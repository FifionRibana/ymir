//! Issue #145 follow-up — **morphological credibility of the closures
//! on the corrected (rigid) transport**.
//!
//! The buoyancy fix proved continents HOLD (global form). This
//! diagnostic asks the orthogonal question the fix does NOT answer:
//! do the closures produce credible *fine* morphology on the good
//! transport, or artefacts (isolated/grid-aligned pics, filaments,
//! uniform shaving)? "Continents hold" ≠ "closures produce good
//! morphology" — measure, do not presume.
//!
//! Protocol (per the user's point):
//!   - Regime: rigid production (`rigid_continental_crust = true`),
//!     closures ON, gallery path (single `run_with_closures`, the
//!     authoritative visual reference — no reclassify).
//!   - **Leave-one-out (LOO) ablation** isolates each closure's
//!     MORPHOLOGICAL contribution: run the full stack, then re-run
//!     with one closure disabled; the delta is what that closure
//!     adds/removes.
//!   - Metric = SPATIAL ([`LandMorphology`]: area, perim/area,
//!     n_components, largest) on the land mask (altitude > 0, the
//!     bipolar render sea level) + VISUAL (PNGs, read by eye).
//!     NOT scalar alone (lesson of the thread).
//!   - Multi-seed {42, 1337, 99} — 42 has cratons; 1337/99 are the
//!     craton-diversity seeds flagged in FOLLOWUPS.md.
//!
//! Closures ablated (the user's five): Davis-Suppe (orogen chains),
//! subduction+accretion (margins), erosion (relief), Stein-Stein
//! (oceans), Track D (subduction+accretion+rifting → evolution
//! trajectory). Equilibrium-height is NOT ablated — it is the
//! load-bearing regulator (disabling it is a known +201% runaway,
//! not a morphology question).
//!
//! Invocation:
//! ```bash
//! cargo test --release -p ymir-core \
//!     --test c1_closure_morphology -- --ignored --nocapture
//! ```
//! Files NOT committed (Phase 1.x + Track A/B/D gallery convention).

use std::path::{Path, PathBuf};

use image::{ImageBuffer, Rgb};

use ymir_core::grid::GridF32;
use ymir_core::morphology::{land_morphology, LandMorphology};
use ymir_core::tectonics_c1::boundary_classification::{classify_boundaries, BoundaryType};
use ymir_core::tectonics::isostasy::{compute_isostasy, IsostasyConfig};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::source_term::apply_stein_stein_bathymetry;
use ymir_core::tectonics_c1::init_r7::{init_c1_state_phase_2_r7, Phase2InitParams};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::state::C1State;
use ymir_core::tectonics_c1::time_loop::{run_with_closures, C1Closures, C1TimeLoopConfig};
use ymir_core::tectonics_v2::boundaries::plate_type::PlateType;
use ymir_core::tectonics_v2::field::Field2D;
use ymir_core::seed::WorldSeed;
use ymir_core::terrain::upscale::FbmUpscaleConfig;
use ymir_core::tectonics_c1::production_upscale::{c1_production_altitude, upscale_from_c1};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::SteinSteinParams;

const GRID_SIZE: usize = 64;
const N_STEPS: usize = 300;
const S_VIZ_MAX: f64 = 3.0;
const ALT_HALF: f32 = 1.13;
const SEEDS: [u64; 4] = [2, 2026, 1988, 4138];

fn output_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/c1_continental_buoyancy/closure_morphology")
}

/// A named LOO ablation: mutate the default (all-ON) closure set.
struct Ablation {
    tag: &'static str,
    mutate: fn(&mut C1Closures),
}

const ABLATIONS: &[Ablation] = &[
    Ablation { tag: "full", mutate: |_| {} },
    Ablation { tag: "no_davis_suppe", mutate: |c| c.davis_suppe.enabled = false },
    Ablation { tag: "no_erosion", mutate: |c| c.erosion.enabled = false },
    Ablation {
        tag: "no_subduction_accretion",
        mutate: |c| {
            c.subduction.enabled = false;
            c.accretion.enabled = false;
        },
    },
    Ablation { tag: "no_stein_stein", mutate: |c| c.oceanic_bathymetry.enabled = false },
    Ablation {
        tag: "no_track_d",
        mutate: |c| {
            c.subduction.enabled = false;
            c.accretion.enabled = false;
            c.rifting.enabled = false;
        },
    },
];

#[test]
#[ignore]
fn closure_morphology_loo_ablation() {
    let dir = output_dir();
    std::fs::create_dir_all(&dir).expect("create output dir");
    let iso_config = IsostasyConfig::c1_default();

    eprintln!("#145 follow-up — closure morphology LOO ablation (rigid production)");
    eprintln!("  grid={GRID_SIZE}² steps={N_STEPS} seeds={SEEDS:?} rigid=true closures=ON");
    eprintln!();
    eprintln!(
        "  {:<26} {:>5} {:>9} {:>7} {:>9} {:>8} | {:>8} {:>8}",
        "config / seed", "land%", "perim/A", "n_comp", "largest", "cont%", "S̃mean", "S̃max"
    );

    for &seed in SEEDS.iter() {
        // Cratonic census per seed (with/without craton diversity).
        let craton_cells = {
            let st = init_c1_state_phase_2_r7(GRID_SIZE, seed, &Phase2InitParams::default());
            st.cratonic_mask.data().iter().filter(|&&b| b).count()
        };
        eprintln!("  --- seed {seed} (cratonic cells at init = {craton_cells}) ---");

        for ab in ABLATIONS {
            let mut closures = C1Closures::default();
            (ab.mutate)(&mut closures);

            let mut state =
                init_c1_state_phase_2_r7(GRID_SIZE, seed, &Phase2InitParams::default());
            let mut kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
            let config = C1TimeLoopConfig {
                rigid_continental_crust: true,
                n_steps: N_STEPS,
                dx: 1.0 / GRID_SIZE as f64,
                dy: 1.0 / GRID_SIZE as f64,
                iso_config: iso_config.clone(),
                drainage_max_distance: 30,
            };

            // Trajectory snapshots only for the full stack at seed 42
            // (Track D evolution credibility — judge the trajectory,
            // not just the final state).
            let traj = ab.tag == "full" && seed == 42;
            let traj_steps: [usize; 4] = [49, 99, 199, 299];

            if traj {
                dump_altitude(&state, 0, &dir, "traj42", &iso_config, &closures);
            }
            run_with_closures(
                &mut state,
                &mut kinematics,
                &config,
                &closures,
                |step, st| {
                    if traj && traj_steps.contains(&step) {
                        dump_altitude(st, step + 1, &dir, "traj42", &iso_config, &closures);
                    }
                },
            );

            // Final-state altitude (post-isostasy + Architecture-C S-S
            // re-apply, matching the gallery render).
            let altitude = render_altitude(&state, &iso_config, &closures);
            let mask: Vec<bool> = altitude.data.iter().map(|&v| v > 0.0).collect();
            let m: LandMorphology = land_morphology(&mask, GRID_SIZE, GRID_SIZE);

            let cont = state
                .plate_type
                .data()
                .iter()
                .filter(|&&t| matches!(t, PlateType::Continental))
                .count();
            let cont_pct = 100.0 * cont as f64 / (GRID_SIZE * GRID_SIZE) as f64;

            let (s_mean, s_max) = {
                let d = state.s.data();
                let sum: f64 = d.iter().sum();
                (sum / d.len() as f64, d.iter().cloned().fold(f64::MIN, f64::max))
            };

            eprintln!(
                "  {:<26} {:>5.1} {:>9.3} {:>7} {:>9.3} {:>8.1} | {:>8.3} {:>8.3}",
                format!("{}/{}", ab.tag, seed),
                100.0 * m.area_fraction,
                m.perimeter_over_area,
                m.n_components,
                m.largest_component_fraction,
                cont_pct,
                s_mean,
                s_max
            );

            // Final PNGs for visual review.
            save_altitude(&altitude, &dir.join(format!("{}_seed{:05}_altitude.png", ab.tag, seed)));
            save_s(&state.s, &dir.join(format!("{}_seed{:05}_s.png", ab.tag, seed)));
        }
        eprintln!();
    }

    eprintln!("  output dir = {}", dir.display());
    eprintln!("  Read the PNGs: full vs each LOO isolates that closure's morphological imprint.");
}

/// H1/H2/H3 discriminator — does the morphology that DS sculpts in S̃
/// become VISIBLE in altitude at higher resolution?
///
/// Runs the full rigid stack at 64² / 128² / 256² (seed 42), holding
/// **physical time** constant (`n_steps ∝ grid`, since CFL `dt ∝ dx ∝
/// 1/grid`). For each grid, measures relief AMPLITUDE restricted to
/// continental cells (where DS thickening lives) — std + p95−p05 — in
/// BOTH S̃ and rendered altitude, plus the conversion ratio.
///
/// - If altitude relief grows richer with resolution → H2 (64² too
///   coarse to carry detail).
/// - If altitude relief stays poor while S̃ relief is rich at every
///   resolution → H1 (S̃→altitude conversion is the bottleneck).
/// PNGs scaled to a constant ~512 px display so the eye compares the
/// SAME physical area at increasing cell counts.
#[test]
#[ignore]
fn resolution_expression_64_128_256() {
    let dir = output_dir().join("resolution");
    std::fs::create_dir_all(&dir).expect("create resolution dir");
    let iso_config = IsostasyConfig::c1_default();
    let seed = 42u64;
    let grids: [usize; 3] = [64, 128, 256];

    eprintln!("#145 follow-up — S̃→altitude expression vs resolution (seed {seed}, full rigid stack)");
    eprintln!("  physical time held constant (n_steps = 300·grid/64)");
    eprintln!(
        "  {:>5} {:>7} {:>6} | {:>9} {:>9} | {:>9} {:>9} | {:>8}",
        "grid", "n_steps", "land%", "S̃ std", "S̃ p95-5", "alt std", "alt p95-5", "alt/S̃"
    );

    for &grid in grids.iter() {
        let n_steps = 300 * grid / 64;
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true,
            n_steps,
            dx: 1.0 / grid as f64,
            dy: 1.0 / grid as f64,
            iso_config: iso_config.clone(),
            drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kinematics, &config, &closures, |_, _| {});

        let altitude = render_altitude(&state, &iso_config, &closures);

        // Restrict to continental cells — where DS thickening lives.
        let mut s_cont: Vec<f64> = Vec::new();
        let mut a_cont: Vec<f64> = Vec::new();
        for j in 0..grid {
            for i in 0..grid {
                if matches!(state.plate_type.get(i, j), PlateType::Continental) {
                    s_cont.push(state.s.get(i, j));
                    a_cont.push(altitude.get(i as i32, j as i32) as f64);
                }
            }
        }
        let (s_std, s_spread) = std_and_spread(&mut s_cont.clone());
        let (a_std, a_spread) = std_and_spread(&mut a_cont.clone());
        let land_pct = 100.0 * s_cont.len() as f64 / (grid * grid) as f64;

        eprintln!(
            "  {:>5} {:>7} {:>6.1} | {:>9.4} {:>9.4} | {:>9.4} {:>9.4} | {:>8.3}",
            format!("{grid}²"),
            n_steps,
            land_pct,
            s_std,
            s_spread,
            a_std,
            a_spread,
            if s_std > 1e-9 { a_std / s_std } else { 0.0 }
        );

        // Constant ~512 px display: scale = 512/grid.
        let scale = (512 / grid).max(1) as u32;
        save_altitude_scaled(&altitude, &dir.join(format!("grid{grid:04}_altitude.png")), scale);
        save_s_scaled(&state.s, &dir.join(format!("grid{grid:04}_s.png")), scale);
    }

    eprintln!();
    eprintln!("  output dir = {}", dir.display());
    eprintln!("  H2 if alt std/spread climb with grid; H1 if alt stays poor while S̃ rich.");
}

/// MESH-CONVERGENCE measure (#145 follow-up — invariance chantier).
///
/// Sweep 64²/128²/256²/512², same seed, full rigid stack, **physical
/// time constant** (`n_steps ∝ grid`). A well-posed model produces the
/// SAME structure at every mesh (finer, not different). We measure
/// whether structure CONVERGES (continents/orogens/bathymetry stable)
/// or DIVERGES (formations appear/disappear with mesh = calibration
/// tuned to 64²).
///
/// Proximity hypothesis (the key mechanism to test): orogens at 64² are
/// EDGE chains born where a convergence boundary falls near a continent
/// coast — with few cells between them, Davis-Suppe piles a chain at
/// that proximity. At 512² the two features are many cells apart, the
/// proximity interaction vanishes, the chain is not born. Test: are
/// orogen cells systematically CLOSER to the coast than the continental
/// average, and does that orogen population COLLAPSE with resolution?
/// All distances reported in PHYSICAL units (cells/grid) for
/// cross-mesh comparison.
#[test]
#[ignore]
fn mesh_convergence_sweep() {
    let dir = output_dir().join("mesh_convergence");
    std::fs::create_dir_all(&dir).expect("create mesh_convergence dir");
    let iso_config = IsostasyConfig::c1_default();
    let seed = 42u64;
    let grids: [usize; 4] = [64, 128, 256, 512];
    // Orogen = continental cell whose S̃ is well above the ~1.0 base
    // continental thickness (clearly DS-thickened, below the ~2.18 cap).
    const OROGEN_S: f64 = 1.5;
    // Physical proximity band: 5 % of the domain.
    const PROX_PHYS: f64 = 0.05;

    eprintln!("#145 follow-up — MESH CONVERGENCE sweep (seed {seed}, full rigid stack)");
    eprintln!("  physical time constant (n_steps = 300·grid/64)");
    eprintln!("  wedge = ANY cell S̃>{OROGEN_S} (DS/accretion thickening, incl. oceanic margin piles)");
    eprintln!(
        "  {:>5} {:>7} {:>6} {:>8} | {:>7} {:>6} {:>11} | {:>9} {:>9} | {:>8}",
        "grid", "n_steps", "land%", "largest", "wedge%", "wedgeN", "wedge d̄coast", "S̃→64 r",
        "alt→64 r", "time"
    );

    // Reference 64² fields (block-averaged to 64² = identity) for the
    // downsample-correlation convergence test.
    let mut ref_s64: Vec<f64> = Vec::new();
    let mut ref_a64: Vec<f64> = Vec::new();

    for &grid in grids.iter() {
        let n_steps = 300 * grid / 64;
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true,
            n_steps,
            dx: 1.0 / grid as f64,
            dy: 1.0 / grid as f64,
            iso_config: iso_config.clone(),
            drainage_max_distance: 30,
        };
        // Snapshot the INIT continental mask (plate geometry before
        // any dynamics) for the boundary-shape A/B test (/btw).
        let mut init_cont = vec![false; grid * grid];
        for j in 0..grid {
            for i in 0..grid {
                if matches!(state.plate_type.get(i, j), PlateType::Continental) {
                    init_cont[j * grid + i] = true;
                }
            }
        }

        let t0 = std::time::Instant::now();
        run_with_closures(&mut state, &mut kinematics, &config, &closures, |_, _| {});
        let elapsed = t0.elapsed();

        let altitude = render_altitude(&state, &iso_config, &closures);

        // Continental (land) mask + wedge mask (ANY plate, S̃>OROGEN_S).
        let mut cont = vec![false; grid * grid];
        let mut wedge = vec![false; grid * grid];
        for j in 0..grid {
            for i in 0..grid {
                let k = j * grid + i;
                if matches!(state.plate_type.get(i, j), PlateType::Continental) {
                    cont[k] = true;
                }
                if state.s.get(i, j) > OROGEN_S {
                    wedge[k] = true;
                }
            }
        }
        let cont_n = cont.iter().filter(|&&b| b).count();
        let wedge_n = wedge.iter().filter(|&&b| b).count();
        let m = land_morphology(&cont, grid, grid);

        // Wedge cells' physical distance to the land/ocean coast.
        let dist = dist_to_coast(&cont, grid);
        let inv_grid = 1.0 / grid as f64;
        let wedge_dcoast = {
            let (mut s, mut c) = (0.0, 0usize);
            for k in 0..grid * grid {
                if wedge[k] && dist[k] != usize::MAX {
                    s += dist[k] as f64 * inv_grid;
                    c += 1;
                }
            }
            if c > 0 { s / c as f64 } else { 0.0 }
        };
        let _ = PROX_PHYS;

        // Downsample S̃ and altitude to 64² (block mean) for the
        // field-to-field convergence correlation vs the 64² reference.
        let s_ds = block_mean_to_64(&|i, j| state.s.get(i, j), grid);
        let a_ds = block_mean_to_64(&|i, j| altitude.get(i as i32, j as i32) as f64, grid);
        let (sr, ar) = if grid == 64 {
            ref_s64 = s_ds.clone();
            ref_a64 = a_ds.clone();
            (1.0, 1.0)
        } else {
            (pearson(&ref_s64, &s_ds), pearson(&ref_a64, &a_ds))
        };

        eprintln!(
            "  {:>5} {:>7} {:>6.1} {:>8.3} | {:>7.1} {:>6} {:>11.4} | {:>9.4} {:>9.4} | {:.1?}",
            format!("{grid}²"),
            n_steps,
            100.0 * cont_n as f64 / (grid * grid) as f64,
            m.largest_component_fraction,
            100.0 * wedge_n as f64 / (grid * grid) as f64,
            m_components(&wedge, grid),
            wedge_dcoast,
            sr,
            ar,
            elapsed,
        );

        // /btw — boundary-shape signature: continental perim/area
        // (final vs init) + how far the coast moved from init (IoU).
        let m_init = land_morphology(&init_cont, grid, grid);
        let (mut inter, mut uni) = (0usize, 0usize);
        for k in 0..grid * grid {
            if cont[k] || init_cont[k] {
                uni += 1;
                if cont[k] && init_cont[k] {
                    inter += 1;
                }
            }
        }
        let iou = if uni > 0 { inter as f64 / uni as f64 } else { 0.0 };
        eprintln!(
            "          boundary: perim/A init={:.3} final={:.3} (Δ {:+.3})   IoU(final,init)={:.3}",
            m_init.perimeter_over_area,
            m.perimeter_over_area,
            m.perimeter_over_area - m_init.perimeter_over_area,
            iou,
        );

        let scale = (512 / grid).max(1) as u32;
        save_altitude_scaled(&altitude, &dir.join(format!("grid{grid:04}_altitude.png")), scale);
        save_s_scaled(&state.s, &dir.join(format!("grid{grid:04}_s.png")), scale);
    }

    eprintln!();
    eprintln!("  output dir = {}", dir.display());
    eprintln!("  CONVERGENCE: large structure (land%, largest, geo) stable across mesh.");
    eprintln!("  PROXIMITY HYPOTHESIS: orogens coast-hugging (oro d̄coast << cont d̄coast,");
    eprintln!("    high oro<5%) at 64² and the orogen population (oro%cont) COLLAPSING with mesh.");
}

/// /btw test #1 — is the 64² boundary irregularity TECTONIC (localised at
/// convergent plate boundaries → Lecture B) or UNIFORM grid crenellation
/// independent of tectonics (→ Lecture A)? Two discriminators at 64²:
///   (1a) coast jaggedness (ocean-neighbour count per coast cell) at
///        convergent-adjacent coast vs elsewhere — equal ⇒ A.
///   (1b) where init≠final coast CHANGES land — enrichment at convergent
///        boundaries vs base rate — ≈1 ⇒ A (dynamics imprint nothing
///        special at convergence).
#[test]
#[ignore]
fn boundary_ab_test_64() {
    let grid = 64usize;
    let seed = 42u64;
    let iso_config = IsostasyConfig::c1_default();
    let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
    let mut kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);

    // Convergent plate-boundary mask from the INIT geometry (where the
    // edge dynamics WOULD imprint over the run).
    let binfo = classify_boundaries(&state.plate_id, &kinematics);
    let mut conv = vec![false; grid * grid];
    for j in 0..grid {
        for i in 0..grid {
            if matches!(binfo.boundary_type.get(i, j), BoundaryType::Convergent) {
                conv[j * grid + i] = true;
            }
        }
    }
    // Dilate convergent mask by 1 cell (a coast cell counts as
    // convergent-adjacent if a convergent boundary is within 1 cell).
    let conv_near = dilate(&conv, grid, 1);

    let mut init_cont = vec![false; grid * grid];
    for j in 0..grid {
        for i in 0..grid {
            if matches!(state.plate_type.get(i, j), PlateType::Continental) {
                init_cont[j * grid + i] = true;
            }
        }
    }

    let config = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / grid as f64,
        dy: 1.0 / grid as f64,
        iso_config: iso_config.clone(),
        drainage_max_distance: 30,
    };
    let closures = C1Closures::default();
    run_with_closures(&mut state, &mut kinematics, &config, &closures, |_, _| {});

    let mut final_cont = vec![false; grid * grid];
    for j in 0..grid {
        for i in 0..grid {
            if matches!(state.plate_type.get(i, j), PlateType::Continental) {
                final_cont[j * grid + i] = true;
            }
        }
    }

    let idx = |i: usize, j: usize| j * grid + i;
    let ocean_neighbours = |mask: &[bool], i: usize, j: usize| -> usize {
        let mut n = 0;
        for (di, dj) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
            let (ni, nj) = (i as i32 + di, j as i32 + dj);
            if ni < 0 || nj < 0 || ni >= grid as i32 || nj >= grid as i32
                || !mask[idx(ni as usize, nj as usize)]
            {
                n += 1;
            }
        }
        n
    };

    // (1a) coast jaggedness, split by convergent-adjacency.
    let (mut rj_conv, mut nc, mut rj_else, mut ne) = (0usize, 0usize, 0usize, 0usize);
    let mut coast_total = 0usize;
    let mut coast_conv = 0usize;
    for j in 0..grid {
        for i in 0..grid {
            if !final_cont[idx(i, j)] {
                continue;
            }
            let on = ocean_neighbours(&final_cont, i, j);
            if on == 0 {
                continue; // interior, not coast
            }
            coast_total += 1;
            if conv_near[idx(i, j)] {
                coast_conv += 1;
                rj_conv += on;
                nc += 1;
            } else {
                rj_else += on;
                ne += 1;
            }
        }
    }
    let rough_conv = if nc > 0 { rj_conv as f64 / nc as f64 } else { 0.0 };
    let rough_else = if ne > 0 { rj_else as f64 / ne as f64 } else { 0.0 };

    // (1b) where the coast CHANGED (init≠final), enrichment at convergent.
    let (mut changed, mut changed_conv) = (0usize, 0usize);
    for k in 0..grid * grid {
        if init_cont[k] != final_cont[k] {
            changed += 1;
            if conv_near[k] {
                changed_conv += 1;
            }
        }
    }
    let base_conv_rate = coast_conv as f64 / coast_total.max(1) as f64;
    let changed_conv_rate = if changed > 0 { changed_conv as f64 / changed as f64 } else { 0.0 };
    let enrichment = if base_conv_rate > 1e-9 { changed_conv_rate / base_conv_rate } else { 0.0 };

    eprintln!("#145 /btw test #1 — boundary origin at 64² (seed {seed}, full rigid stack)");
    eprintln!("  convergent-adjacent coast = {coast_conv}/{coast_total} ({:.1}%)", 100.0 * base_conv_rate);
    eprintln!("  (1a) coast jaggedness (ocean-neighbours/coast-cell):");
    eprintln!("       convergent-adjacent = {rough_conv:.3}   elsewhere = {rough_else:.3}   ratio = {:.3}", if rough_else>1e-9 {rough_conv/rough_else} else {0.0});
    eprintln!("  (1b) coast changes (init≠final) = {changed} cells; {changed_conv} at convergent");
    eprintln!("       changed-conv rate = {:.1}%  vs base {:.1}%  → enrichment = {enrichment:.2}×", 100.0*changed_conv_rate, 100.0*base_conv_rate);
    eprintln!();
    eprintln!("  A (grid noise, dynamics don't imprint): jaggedness ratio ≈1 AND enrichment ≈1.");
    eprintln!("  B (tectonic edge-shaping): jaggedness higher at convergent AND enrichment >1.");
}

/// Dilate a boolean mask by `r` cells (4-neighbour Chebyshev-ish via
/// repeated 4-neighbour expansion).
fn dilate(mask: &[bool], grid: usize, r: usize) -> Vec<bool> {
    let mut cur = mask.to_vec();
    let idx = |i: usize, j: usize| j * grid + i;
    for _ in 0..r {
        let mut next = cur.clone();
        for j in 0..grid {
            for i in 0..grid {
                if cur[idx(i, j)] {
                    continue;
                }
                for (di, dj) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
                    let (ni, nj) = (i as i32 + di, j as i32 + dj);
                    if ni >= 0 && nj >= 0 && ni < grid as i32 && nj < grid as i32
                        && cur[idx(ni as usize, nj as usize)]
                    {
                        next[idx(i, j)] = true;
                        break;
                    }
                }
            }
        }
        cur = next;
    }
    cur
}

/// #147 — contrast counterfactual, variant (γ): FEASIBILITY of the
/// scheme. Question: "can the curtain cede?" Advection-only + rigid
/// no-flux (= the pure curtain, variant C r≈0.045) PLUS a local S̃
/// band-smoothing across the continental/oceanic boundary (softens the
/// sharp 1.0/0.2 contrast the upwind oscillates on). Reuses the REAL
/// `run_advection_only` one step at a time (no churn, no reinvention),
/// smoothing `state.s` between steps. NO closures (no EH) — this
/// isolates "smoothing → curtain", proving scheme FEASIBILITY only,
/// NOT the net system effect (EH already bounds the curtain to ~0.78
/// in the full system).
///
/// Band width = fixed PHYSICAL 2.5 reference-cells (`2.5/64 · grid`),
/// so the band is itself mesh-invariant. λ = smoothing strength per
/// step toward the 4-neighbour mean (0 = off = pure variant C).
///
/// Verdict: λ lifts r→~1 ⇒ NOT a scheme floor (upwind CAN converge on
/// a softened contrast) → proceed to (α) full-system net effect. λ
/// leaves r ~0.045–0.2 ⇒ SCHEME floor → redefine the criterion
/// (r~0.78) without touching the production loop.
#[test]
#[ignore]
fn contrast_counterfactual_gamma() {
    let iso_config = IsostasyConfig::c1_default();
    let seed = 42u64;
    let grids: [usize; 3] = [64, 128, 256];
    let lambdas: [f64; 4] = [0.0, 0.25, 0.5, 1.0];
    const BAND_PHYS_CELLS: f64 = 2.5; // physical band width at ref-64

    eprintln!("#147 — contrast counterfactual (γ): advection-only + S̃ band-smoothing");
    eprintln!("  isolates 'softened contrast → curtain'; proves SCHEME FEASIBILITY only (no EH)");
    eprintln!("  band = {BAND_PHYS_CELLS} physical cells (2.5/64·grid); λ = smooth strength/step");
    eprintln!("  {:>6} | {:>8} {:>8} {:>8}", "λ", "64²", "128²", "256²");

    for &lambda in lambdas.iter() {
        let mut ref64: Vec<f64> = Vec::new();
        let mut row = String::new();
        for &grid in grids.iter() {
            let n_steps = 300 * grid / 64;
            let band_cells = (BAND_PHYS_CELLS / 64.0 * grid as f64).round() as usize;
            let mut state =
                init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
            let kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);

            // Continental/oceanic boundary band (static under advection-only).
            let mut cont = vec![false; grid * grid];
            for j in 0..grid {
                for i in 0..grid {
                    cont[j * grid + i] =
                        matches!(state.plate_type.get(i, j), PlateType::Continental);
                }
            }
            let band = boundary_band(&cont, grid, band_cells);

            let step_cfg = C1TimeLoopConfig {
                rigid_continental_crust: true,
                n_steps: 1,
                dx: 1.0 / grid as f64,
                dy: 1.0 / grid as f64,
                iso_config: iso_config.clone(),
                drainage_max_distance: 30,
            };
            let mut kin = kinematics.clone();
            for _ in 0..n_steps {
                // One real advection step (rigid no-flux), then smooth.
                ymir_core::tectonics_c1::time_loop::run_advection_only(
                    &mut state, &kin, &step_cfg, |_, _| {},
                );
                if lambda > 0.0 {
                    smooth_band(&mut state.s, &band, grid, lambda);
                }
                let _ = &mut kin;
            }

            let s_ds = block_mean_to_64(&|i, j| state.s.get(i, j), grid);
            let r = if grid == 64 {
                ref64 = s_ds.clone();
                1.0
            } else {
                pearson(&ref64, &s_ds)
            };
            row.push_str(&format!(" {r:>8.4}"));
        }
        eprintln!("  {lambda:>6.2} |{row}");
    }
    eprintln!();
    eprintln!("  λ=0 row = pure curtain (variant C, expect ~0.045).");
    eprintln!("  r→~1 with λ ⇒ scheme FEASIBLE (curtain cedes) → go (α). r stuck ⇒ scheme floor.");
}

/// #147 FOLLOWUPS-#6 GATING measurement — is the upscale ROBUST to the
/// S̃ field non-convergence (production r~0.51), or does it diverge?
/// Wires C1→upscale on a throwaway: run production (rigid, full
/// closures) at 64² and 256² (same seed), build the gallery altitude,
/// normalise with the SAME fixed map (sea=0.5) so both are comparable,
/// upscale both to 1024² anisotropic FBM, and compare STRUCTURE
/// (NOT identity — FBM fills more from 256² by design, "different
/// paths" accepted). Criterion = structure convergence: same
/// continents/chains/bathymetry in the same places, detail differing.
///
/// Metrics: (1) coarse-altitude cross-grid r (the field the upscale's
/// slope orientation reads); (2) upscaled-result structure r (both
/// 1024² block-meaned to 64², correlated — FBM averages out, large
/// structure remains); (3) land morphology of each upscaled result
/// (sea_level threshold) — same n_components/largest/area? (4) visual
/// side-by-side PNGs (the judge).
#[test]
#[ignore]
fn upscale_from_c1_structure_converges() {
    let dir = output_dir().join("upscale_robustness");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso_config = IsostasyConfig::c1_default();
    let seed = 42u64;
    const HALF: f32 = ALT_HALF; // fixed normalisation half-range (sea at 0.0 → 0.5)
    let cfg = FbmUpscaleConfig { target_size: 1024, ..Default::default() };

    eprintln!("#147 FOLLOWUPS-#6 contract regression — upscale_from_c1 structure convergence (seed {seed})");
    eprintln!("  structure-convergence (NOT identity); coarse normalised sea=0.5, target 1024²");

    let mut coarse_ds: Vec<Vec<f64>> = Vec::new(); // coarse-altitude →64
    let mut up_ds: Vec<Vec<f64>> = Vec::new();      // upscaled →64
    let mut up_land: Vec<f64> = Vec::new();
    let mut up_largest: Vec<f64> = Vec::new();
    let grids = [64usize, 256];
    for &grid in &grids {
        let n_steps = 300 * grid / 64;
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true,
            n_steps,
            dx: 1.0 / grid as f64,
            dy: 1.0 / grid as f64,
            iso_config: iso_config.clone(),
            drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});

        // Gallery altitude → fixed [0,1] normalisation (sea 0.0 → 0.5).
        let altitude = render_altitude(&state, &iso_config, &closures);
        let mut coarse = altitude.clone();
        for v in coarse.data.iter_mut() {
            *v = ((*v + HALF) / (2.0 * HALF)).clamp(0.0, 1.0);
        }
        coarse_ds.push(block_mean_to_64(&|i, j| coarse.get(i as i32, j as i32) as f64, grid));

        // Upscale through the CONTRACT path (upscale_from_c1 builds the
        // laundered altitude internally — never raw S̃). This is what the
        // regression guards: production uses exactly this entry.
        let up = upscale_from_c1(
            &state,
            &iso_config,
            &closures.oceanic_bathymetry,
            &WorldSeed::new(seed),
            &cfg,
        );
        let up_grid = up.heightmap.width;
        up_ds.push(block_mean_to_64_dim(
            &|i, j| up.heightmap.get(i as i32, j as i32) as f64,
            up_grid,
        ));

        // Land morphology of the upscaled result (land = h > sea 0.5).
        let mask: Vec<bool> = up.heightmap.data.iter().map(|&v| v > 0.5).collect();
        let m = land_morphology(&mask, up.heightmap.width, up.heightmap.height);
        eprintln!(
            "  grid {grid}²→{}²: upscaled land%={:.1} perim/A={:.3} n_comp={} largest={:.3}",
            up_grid, 100.0 * m.area_fraction, m.perimeter_over_area, m.n_components,
            m.largest_component_fraction
        );
        up_land.push(m.area_fraction);
        up_largest.push(m.largest_component_fraction);

        save_heightmap01(&up.heightmap, &dir.join(format!("upscaled_from{grid:04}.png")));
    }

    let coarse_r = pearson(&coarse_ds[0], &coarse_ds[1]);
    let up_r = pearson(&up_ds[0], &up_ds[1]);
    eprintln!();
    eprintln!("  coarse-altitude structure r (64 vs 256, →64) = {coarse_r:.4}");
    eprintln!("  UPSCALED structure r       (64 vs 256, →64) = {up_r:.4}");
    eprintln!("  Visual: {}", dir.display());

    // REGRESSION (Issue #147 #6 contract): the upscale, fed via
    // `upscale_from_c1` (laundered altitude), must stay STRUCTURE-
    // convergent across resolutions. If a future change feeds raw S̃
    // instead, up_r collapses toward the S̃ field r (~0.51) and this
    // fails — the precondition is TESTED, not merely documented.
    assert!(
        up_r >= 0.85,
        "upscaled structure r {up_r:.3} < 0.85 — robustness contract broken \
         (is the upscale reading raw S̃ instead of the laundered altitude? reopens #6)"
    );
    // Same world: one dominant landmass at BOTH resolutions, comparable
    // land fraction (FBM detail differs; large structure must not).
    assert!(
        up_largest[0] > 0.8 && up_largest[1] > 0.8,
        "upscaled largest-component {up_largest:?} — not one dominant landmass at both res"
    );
    assert!(
        (up_land[0] - up_land[1]).abs() < 0.05,
        "upscaled land fractions {up_land:?} diverge > 5 pts across resolution"
    );
}

/// Block-mean an arbitrary square `grid` (multiple of 64) field to 64².
fn block_mean_to_64_dim(get: &dyn Fn(usize, usize) -> f64, grid: usize) -> Vec<f64> {
    block_mean_to_64(get, grid)
}

/// Save a native crop with a CONTINUOUS grayscale palette (no
/// hypsometric color bands), so palette-banding can be told apart from
/// real height (FBM) variation.
fn save_gray01_crop(h: &GridF32, x0: usize, y0: usize, size: usize, path: &Path) {
    let w = (x0 + size).min(h.width) - x0;
    let ht = (y0 + size).min(h.height) - y0;
    let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(w as u32, ht as u32);
    for jj in 0..ht {
        for ii in 0..w {
            let v = h.get((x0 + ii) as i32, (y0 + jj) as i32).clamp(0.0, 1.0);
            let g = (v * 255.0) as u8;
            img.put_pixel(ii as u32, (ht - 1 - jj) as u32, Rgb([g, g, g]));
        }
    }
    img.save(path).expect("save gray crop PNG");
}

/// Hillshade crop (Lambert, light from NW) — shows REAL relief without
/// any palette colour-band artefact. Flat land → uniform mid-gray;
/// FBM ripples → visible shading. dz exaggerated for visibility.
fn save_hillshade_crop(h: &GridF32, x0: usize, y0: usize, size: usize, path: &Path) {
    let w = (x0 + size).min(h.width) - x0;
    let ht = (y0 + size).min(h.height) - y0;
    let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(w as u32, ht as u32);
    let z = 60.0_f64; // vertical exaggeration
    let (lx, ly, lz) = {
        let (a, b, c) = (1.0_f64, 1.0, 2.0);
        let n = (a * a + b * b + c * c).sqrt();
        (a / n, b / n, c / n)
    };
    for jj in 0..ht {
        for ii in 0..w {
            let gx = (x0 + ii) as i32;
            let gy = (y0 + jj) as i32;
            let dzdx = (h.get(gx + 1, gy) - h.get(gx - 1, gy)) as f64 * z;
            let dzdy = (h.get(gx, gy + 1) - h.get(gx, gy - 1)) as f64 * z;
            let nn = (dzdx * dzdx + dzdy * dzdy + 1.0).sqrt();
            let shade = ((-dzdx * lx - dzdy * ly + lz) / nn).clamp(0.0, 1.0);
            let g = (shade * 255.0) as u8;
            img.put_pixel(ii as u32, (ht - 1 - jj) as u32, Rgb([g, g, g]));
        }
    }
    img.save(path).expect("save hillshade crop PNG");
}

/// #151 coastal contour-line PIN — are the faint lines near the coast on
/// land real FBM height, or hypsometric-PALETTE banding? Render the same
/// coastal crop in hypsometric vs continuous grayscale. If the lines
/// vanish in grayscale → palette artefact (cosmetic, not terrain).
#[test]
#[ignore]
fn coast_palette_check() {
    let dir = output_dir().join("coast_palette");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso_config = IsostasyConfig::c1_default();
    let seed = std::env::var("SEED").ok().and_then(|s| s.parse().ok()).unwrap_or(1988u64);
    let mut state = init_c1_state_phase_2_r7(64, seed, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    let closures = C1Closures::default();
    let config = C1TimeLoopConfig {
        rigid_continental_crust: true, n_steps: 300,
        dx: 1.0 / 64.0, dy: 1.0 / 64.0,
        iso_config: iso_config.clone(), drainage_max_distance: 30,
    };
    run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
    let band = std::env::var("FBM_BAND").ok().and_then(|s| s.parse().ok()).unwrap_or(0.20);
    let cfg = FbmUpscaleConfig {
        target_size: 2048, coast_warp_strength: 0.8, coast_warp_frequency: 0.5,
        coastal_amplitude_band: band, amplitude_base: 0.16, submarine_damping: 0.0,
        ..Default::default()
    };
    let up = upscale_from_c1(
        &state, &iso_config, &closures.oceanic_bathymetry, &WorldSeed::new(seed), &cfg,
    );
    // 3-way render of coastal crops: hypso (user's view) vs gray (real
    // height) vs hillshade (relief, no palette bands).
    for (cname, x0, y0) in [("a", 900usize, 1050usize), ("b", 1150, 850), ("c", 600, 1500)] {
        save_heightmap01_crop(&up.heightmap, x0, y0, 460, &dir.join(format!("coast_{cname}_hypso.png")));
        save_gray01_crop(&up.heightmap, x0, y0, 460, &dir.join(format!("coast_{cname}_gray.png")));
        save_hillshade_crop(&up.heightmap, x0, y0, 460, &dir.join(format!("coast_{cname}_hill.png")));
    }
    save_heightmap01(&up.heightmap, &dir.join("full_hypso.png"));
    // Band-width sweep (grayscale crop "a") — does a WIDER coastal taper
    // suppress the near-coast FBM ripples?
    for band in [0.06_f64, 0.15, 0.30] {
        let c2 = FbmUpscaleConfig { coastal_amplitude_band: band, ..cfg.clone() };
        let u2 = upscale_from_c1(
            &state, &iso_config, &closures.oceanic_bathymetry, &WorldSeed::new(seed), &c2,
        );
        save_gray01_crop(&u2.heightmap, 900, 1050, 460,
            &dir.join(format!("band_{:03}_gray.png", (band * 100.0) as i32)));
    }
    // Full continents at a few band widths (hypsometric) for the look call.
    for band in [0.06_f64, 0.20, 0.35] {
        let c2 = FbmUpscaleConfig { coastal_amplitude_band: band, ..cfg.clone() };
        let u2 = upscale_from_c1(
            &state, &iso_config, &closures.oceanic_bathymetry, &WorldSeed::new(seed), &c2,
        );
        save_heightmap01(&u2.heightmap, &dir.join(format!("full_band_{:03}.png", (band * 100.0) as i32)));
    }
    eprintln!("  out = {}", dir.display());
}

/// Bresenham line into an RGB image (clamped).
fn draw_line(img: &mut ImageBuffer<Rgb<u8>, Vec<u8>>, x0: i32, y0: i32, x1: i32, y1: i32, c: [u8; 3]) {
    let (w, h) = (img.width() as i32, img.height() as i32);
    let dx = (x1 - x0).abs();
    let dy = -(y1 - y0).abs();
    let sx = if x0 < x1 { 1 } else { -1 };
    let sy = if y0 < y1 { 1 } else { -1 };
    let (mut x, mut y) = (x0, y0);
    let mut err = dx + dy;
    loop {
        if x >= 0 && y >= 0 && x < w && y < h { img.put_pixel(x as u32, y as u32, Rgb(c)); }
        if x == x1 && y == y1 { break; }
        let e2 = 2 * err;
        if e2 >= dy { err += dy; x += sx; }
        if e2 <= dx { err += dx; y += sy; }
    }
}

/// #151 RELIEF profile — 1D cuts ACROSS the convergent arc (seed 42),
/// sampling raw S̃ and the COARSE production altitude (pre-FBM) vs
/// distance, to read whether a marginal Davis-Suppe orogen pic exists
/// in S̃ (force) and whether it survives the S̃→altitude conversion
/// (expression). Arc located via the DS convergence classifier
/// (`classify_boundaries` → Convergent). Same run/fields as
/// `export_relief_compare`. Outputs CSV (the numbers) + a PNG plot per
/// cut; prints cut centres + inward normals so the cut placement is
/// auditable.
#[test]
#[ignore]
fn profile_convergent_arc_seed42() {
    let dir = output_dir().join("relief_profile");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso_config = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let seed = 42u64;
    let grid = 64usize;

    let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    let closures = C1Closures::default();
    let config = C1TimeLoopConfig {
        rigid_continental_crust: true, n_steps: 300,
        dx: 1.0 / grid as f64, dy: 1.0 / grid as f64,
        iso_config: iso_config.clone(), drainage_max_distance: 30,
    };
    run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});

    // Coarse production altitude (pre-FBM): isostasy + Stein-Stein + despike.
    let alt = c1_production_altitude(&state.s, &state.age, &state.plate_type, &iso_config, &ss);
    let binfo = classify_boundaries(&state.plate_id, &kin);
    let cont = |i: usize, j: usize| matches!(state.plate_type.get(i, j), PlateType::Continental);

    // Convergent cells adjacent to land (the orogen sites along the coast).
    let mut conv: Vec<(usize, usize)> = Vec::new();
    for j in 0..grid { for i in 0..grid {
        if matches!(binfo.boundary_type.get(i, j), BoundaryType::Convergent) {
            conv.push((i, j));
        }
    }}
    // COAST convergences (south half): convergent cells with BOTH a
    // continental and an oceanic 4-neighbour → a real coast convergence
    // with a well-defined inward normal (excludes mid-ocean convergences,
    // which gave degenerate normal=0 cuts).
    let has_both = |ci: usize, cj: usize| {
        let (mut land, mut sea) = (false, false);
        for (di, dj) in [(-1i32,0i32),(1,0),(0,-1),(0,1)] {
            let (a, b) = (ci as i32 + di, cj as i32 + dj);
            if a>=0 && b>=0 && (a as usize)<grid && (b as usize)<grid {
                if cont(a as usize, b as usize) { land = true; } else { sea = true; }
            }
        }
        land && sea
    };
    let mut south: Vec<(usize, usize)> = conv.iter().copied()
        .filter(|&(i, j)| j < grid / 2 && has_both(i, j)).collect();
    south.sort_by_key(|&(i, _)| i);
    eprintln!("seed {seed}: {} convergent cells total, {} south COAST convergences", conv.len(), south.len());
    eprintln!("  south coast-convergent cells (i,j): {south:?}");

    // 3 cut centres spread along the south arc.
    let centres: Vec<(usize, usize)> = if south.len() >= 3 {
        [south.len() / 5, south.len() / 2, 4 * south.len() / 5].iter().map(|&k| south[k]).collect()
    } else { south.clone() };

    // Inward normal at a cell: toward continental neighbours.
    let normal = |ci: usize, cj: usize| -> (f64, f64) {
        let (mut nx, mut ny) = (0.0, 0.0);
        for dj in -1i32..=1 { for di in -1i32..=1 {
            if di == 0 && dj == 0 { continue; }
            let (a, b) = (ci as i32 + di, cj as i32 + dj);
            if a >= 0 && b >= 0 && (a as usize) < grid && (b as usize) < grid {
                let s = if cont(a as usize, b as usize) { 1.0 } else { -1.0 };
                nx += s * di as f64; ny += s * dj as f64;
            }
        }}
        let n = (nx * nx + ny * ny).sqrt().max(1e-9);
        (nx / n, ny / n)
    };

    for (ci, &(cx, cy)) in centres.iter().enumerate() {
        let (nx, ny) = normal(cx, cy);
        eprintln!("  cut {ci}: centre=({cx},{cy}) inward_normal=({nx:.2},{ny:.2})  [t<0 ocean → t=0 convergence → t>0 interior]");
        // Sample t = -10..=12 cells along the inward normal.
        let mut rows: Vec<(i32, f64, f64, u8, char)> = Vec::new();
        for t in -10i32..=12 {
            let i = (cx as f64 + t as f64 * nx).round();
            let j = (cy as f64 + t as f64 * ny).round();
            if i < 0.0 || j < 0.0 || i >= grid as f64 || j >= grid as f64 { continue; }
            let (i, j) = (i as usize, j as usize);
            let s = state.s.get(i, j);
            let a = alt.get(i as i32, j as i32) as f64;
            let cv = matches!(binfo.boundary_type.get(i, j), BoundaryType::Convergent) as u8;
            let pt = if cont(i, j) { 'C' } else { 'O' };
            rows.push((t, s, a, cv, pt));
        }
        // CSV (plate_type: C continental [altitude = isostasy(S̃)] / O
        // oceanic [altitude = Stein-Stein bathymetry, S̃ NOT used]).
        let mut csv = String::from("dist_cells,s_thickness,coarse_altitude,convergent,plate_type\n");
        for &(t, s, a, cv, pt) in &rows {
            csv.push_str(&format!("{t},{s:.4},{a:.4},{cv},{pt}\n"));
        }
        std::fs::write(dir.join(format!("seed{seed:05}_cut{ci}.csv")), &csv).unwrap();
        // Console table.
        eprintln!("    dist  S̃       altitude  pt  conv");
        for &(t, s, a, cv, pt) in &rows {
            eprintln!("    {t:>4}  {s:>6.3}  {a:>+7.3}   {pt}  {}", if cv == 1 { "<-- CONV" } else { "" });
        }

        // PNG plot: S̃ (blue) + altitude (brown), each normalised to its
        // own range; gray vertical marker at t=0; gray baseline.
        let (w, h) = (640i32, 360i32);
        let (ml, mr, mt, mb) = (40i32, 20, 20, 30);
        let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::from_pixel(w as u32, h as u32, Rgb([250, 250, 250]));
        let smin = rows.iter().map(|r| r.1).fold(f64::MAX, f64::min);
        let smax = rows.iter().map(|r| r.1).fold(f64::MIN, f64::max);
        let amin = rows.iter().map(|r| r.2).fold(f64::MAX, f64::min);
        let amax = rows.iter().map(|r| r.2).fold(f64::MIN, f64::max);
        let n = rows.len() as i32;
        let xof = |k: i32| ml + k * (w - ml - mr) / (n - 1).max(1);
        let yof = |v: f64, lo: f64, hi: f64| {
            let t = ((v - lo) / (hi - lo).max(1e-9)).clamp(0.0, 1.0);
            (h - mb) - (t * (h - mt - mb) as f64) as i32
        };
        // t=0 marker.
        if let Some(k0) = rows.iter().position(|r| r.0 == 0) {
            draw_line(&mut img, xof(k0 as i32), mt, xof(k0 as i32), h - mb, [150, 150, 150]);
        }
        for k in 1..rows.len() {
            let (x0, x1) = (xof(k as i32 - 1), xof(k as i32));
            draw_line(&mut img,
                x0, yof(rows[k - 1].1, smin, smax), x1, yof(rows[k].1, smin, smax), [30, 60, 200]);
            draw_line(&mut img,
                x0, yof(rows[k - 1].2, amin, amax), x1, yof(rows[k].2, amin, amax), [160, 90, 30]);
        }
        img.save(dir.join(format!("seed{seed:05}_cut{ci}.png"))).unwrap();
        eprintln!("    S̃ range [{smin:.3},{smax:.3}]  altitude range [{amin:.3},{amax:.3}]");
    }
    eprintln!("  out = {} (BLUE=S̃, BROWN=coarse altitude, GRAY=convergence t=0)", dir.display());
}

/// RELIEF EXPRESSION-vs-GENERATION diagnostic — per seed, three
/// co-registered views from the SAME final C1State: (1) HD altitude
/// (production #151 config, the real product where relief is missing),
/// (2) RAW S̃ thickness (no FBM — does Davis-Suppe orogeny exist in the
/// thickness?), (3) boundary TYPES — convergences (red) reusing the DS
/// convergence classifier (`classify_boundaries`), so the analysis can
/// overlay convergence ↔ altitude ↔ S̃: orogen in S̃ + flat altitude =
/// EXPRESSION; flat in both at a convergence = GENERATION.
#[test]
#[ignore]
fn export_relief_compare() {
    let dir = output_dir().join("relief_compare");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso_config = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let seeds: [u64; 7] = [2, 42, 99, 1337, 1988, 2026, 4138];
    let cfg = FbmUpscaleConfig {
        target_size: 1024, coast_warp_strength: 1.5, coast_warp_frequency: 0.5,
        coastal_amplitude_band: 0.30, amplitude_base: 0.16, submarine_damping: 0.0,
        ..Default::default()
    };
    eprintln!("#151 relief compare — altitude(1024² prod) / raw S̃ / convergence map, per seed");
    for &seed in &seeds {
        let mut state = init_c1_state_phase_2_r7(64, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300,
            dx: 1.0 / 64.0, dy: 1.0 / 64.0,
            iso_config: iso_config.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});

        let up = upscale_from_c1(&state, &iso_config, &ss, &WorldSeed::new(seed), &cfg);
        save_heightmap01(&up.heightmap, &dir.join(format!("seed{seed:05}_altitude.png")));
        save_s_scaled(&state.s, &dir.join(format!("seed{seed:05}_sthickness.png")), 8);

        let binfo = classify_boundaries(&state.plate_id, &kin);
        let grid = 64usize;
        let mut bimg = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(grid as u32 * 8, grid as u32 * 8);
        for j in 0..grid { for i in 0..grid {
            let base = if matches!(state.plate_type.get(i, j), PlateType::Continental) {
                [60, 130, 60]
            } else { [40, 80, 160] };
            let col = match binfo.boundary_type.get(i, j) {
                BoundaryType::Convergent => [220, 30, 30],
                BoundaryType::Divergent => [40, 220, 220],
                BoundaryType::Transform => [230, 220, 40],
                _ => base,
            };
            put_block(&mut bimg, i, grid - 1 - j, 8, col);
        }}
        bimg.save(dir.join(format!("seed{seed:05}_boundaries.png"))).unwrap();
        eprintln!("  seed {seed} done");
    }
    eprintln!("  out = {}", dir.display());
    eprintln!("  legend: convergent=RED divergent=CYAN transform=YELLOW; land=green ocean=blue");
}

/// F3 PIN — the dark dotted lines in the OCEAN. Render, at the COARSE
/// 64² resolution (pre-upscale), the production altitude (Stein-Stein
/// bathymetry), the age field, and the plate boundaries, to test:
/// (1) are the dark lines already in the coarse altitude (→ Stein-Stein,
/// not the upscale)? (2) do they coincide with age discontinuities /
/// plate boundaries (→ age-jump → depth-jump)?
#[test]
#[ignore]
fn pin_f3_ocean_lines() {
    let dir = output_dir().join("f3_pin");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso_config = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    for &seed in &[4138u64, 2] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300,
            dx: 1.0 / grid as f64, dy: 1.0 / grid as f64,
            iso_config: iso_config.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});

        let mut amin = f64::INFINITY; let mut amax = f64::NEG_INFINITY;
        for &v in state.age.data() { amin = amin.min(v); amax = amax.max(v); }
        let arange = (amax - amin).max(1e-9);

        // (1) Coarse production altitude (Stein-Stein), 64² scaled ×8.
        let alt = c1_production_altitude(&state.s, &state.age, &state.plate_type, &iso_config, &ss);
        save_altitude_scaled(&alt, &dir.join(format!("f3_seed{seed:05}_altitude.png")), 8);

        // Counterfactual: 3×3 MEDIAN-filtered age → despikes the
        // pile-up cells, then Stein-Stein. Do the dark dots vanish?
        let mut age_med = state.age.clone();
        for j in 0..grid { for i in 0..grid {
            let mut nb: Vec<f64> = Vec::with_capacity(9);
            for dj in -1i32..=1 { for di in -1i32..=1 {
                let (ni, nj) = (i as i32 + di, j as i32 + dj);
                if ni>=0 && nj>=0 && (ni as usize)<grid && (nj as usize)<grid {
                    nb.push(state.age.get(ni as usize, nj as usize));
                }
            }}
            nb.sort_by(|a,b| a.partial_cmp(b).unwrap());
            age_med.set(i, j, nb[nb.len()/2]);
        }}
        let alt_med = c1_production_altitude(&state.s, &age_med, &state.plate_type, &iso_config, &ss);
        save_altitude_scaled(&alt_med, &dir.join(format!("f3_seed{seed:05}_altitude_median.png")), 8);

        // Despiked age, SAME normalisation as the before-age, to verify the
        // median kills ONLY spikes and preserves any legitimate age gradient.
        let mut aimg2 = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(grid as u32 * 8, grid as u32 * 8);
        for j in 0..grid { for i in 0..grid {
            let g = (((age_med.get(i, j) - amin) / arange).clamp(0.0, 1.0) * 255.0) as u8;
            put_block(&mut aimg2, i, grid - 1 - j, 8, [g, g, g]);
        }}
        aimg2.save(dir.join(format!("f3_seed{seed:05}_age_median.png"))).unwrap();

        // (2a) Age field, grayscale normalised to its range.
        let mut aimg = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(grid as u32 * 8, grid as u32 * 8);
        for j in 0..grid { for i in 0..grid {
            let g = (((state.age.get(i, j) - amin) / arange) * 255.0) as u8;
            put_block(&mut aimg, i, grid - 1 - j, 8, [g, g, g]);
        }}
        aimg.save(dir.join(format!("f3_seed{seed:05}_age.png"))).unwrap();

        // (2b) Plate boundaries (cell 4-adjacent to a different plate_id) in red over ocean.
        let idx = |i: usize, j: usize| state.plate_id.get(i, j);
        let mut bimg = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(grid as u32 * 8, grid as u32 * 8);
        for j in 0..grid { for i in 0..grid {
            let mut bnd = false;
            for (di, dj) in [(-1i32,0i32),(1,0),(0,-1),(0,1)] {
                let (ni, nj) = (i as i32 + di, j as i32 + dj);
                if ni>=0 && nj>=0 && (ni as usize)<grid && (nj as usize)<grid
                    && idx(ni as usize, nj as usize) != idx(i, j) { bnd = true; }
            }
            let oceanic = matches!(state.plate_type.get(i,j), PlateType::Oceanic);
            let col = if bnd { [255,0,0] } else if oceanic { [40,80,160] } else { [60,130,60] };
            put_block(&mut bimg, i, grid - 1 - j, 8, col);
        }}
        bimg.save(dir.join(format!("f3_seed{seed:05}_bounds.png"))).unwrap();
        eprintln!("  seed {seed}: age range [{amin:.3},{amax:.3}]");
    }
    eprintln!("  out = {}", dir.display());
}

/// #151 coast-warp strength sweep — with FBM removed from the coast
/// (band 0.30) the warp carries the coastline irregularity alone, so it
/// reads lighter; sweep stronger warp to restore an irregular coast
/// (watch for fragmentation). Production-ish config, 2048².
#[test]
#[ignore]
fn export_warp_sweep() {
    let dir = output_dir().join("warp_sweep");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso_config = IsostasyConfig::c1_default();
    for &seed in &[1988u64, 2] {
        let mut state = init_c1_state_phase_2_r7(64, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300,
            dx: 1.0 / 64.0, dy: 1.0 / 64.0,
            iso_config: iso_config.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        for warp in [0.8_f64, 1.2, 1.6, 2.0] {
            let cfg = FbmUpscaleConfig {
                target_size: 2048, coast_warp_strength: warp, coast_warp_frequency: 0.5,
                coastal_amplitude_band: 0.30, amplitude_base: 0.16, submarine_damping: 0.0,
                ..Default::default()
            };
            let up = upscale_from_c1(
                &state, &iso_config, &closures.oceanic_bathymetry, &WorldSeed::new(seed), &cfg,
            );
            save_heightmap01(&up.heightmap,
                &dir.join(format!("warp_seed{seed:05}_{:02}.png", (warp * 10.0) as i32)));
        }
        eprintln!("  seed {seed} done");
    }
    eprintln!("  out = {}", dir.display());
}

/// #151 PRODUCTION combined export — coast warp + coastal amplitude taper
/// + mountain amplitude, the full recommended HD config, on the seed set.
#[test]
#[ignore]
fn export_hd_production() {
    let dir = output_dir().join("hd_production");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso_config = IsostasyConfig::c1_default();
    let seeds: [u64; 7] = [2, 42, 99, 1337, 1988, 2026, 4138];
    let cfg = FbmUpscaleConfig {
        target_size: 2048,
        coast_warp_strength: 1.5, // stronger: FBM no longer roughens the coast (band 0.30)
        coast_warp_frequency: 0.5,
        coastal_amplitude_band: 0.30, // FBM only in highland; all lowland smooth
        amplitude_base: 0.16,
        submarine_damping: 0.0, // no FBM in the ocean (smooth bathymetry)
        ..Default::default()
    };
    eprintln!("#151 HD production export — warp 0.8 + band 0.30 + amp 0.16 + subdamp 0.0, 2048²");
    for &seed in &seeds {
        let mut state = init_c1_state_phase_2_r7(64, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300,
            dx: 1.0 / 64.0, dy: 1.0 / 64.0,
            iso_config: iso_config.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        let up = upscale_from_c1(
            &state, &iso_config, &closures.oceanic_bathymetry, &WorldSeed::new(seed), &cfg,
        );
        save_heightmap01(&up.heightmap, &dir.join(format!("hd_seed{seed:05}.png")));
        eprintln!("  seed {seed:>5} done");
    }
    eprintln!("  out = {}", dir.display());
}

/// #151 border-FBM artefact PIN — the FBM streaks the interior near the
/// coast (directional ridges shooting inland). Hypothesis: the coast's
/// extreme coarse-step slope drives FULL anisotropy + slope-amplified
/// amplitude → long ridges along the coast-normal. Counterfactual:
/// disable anisotropy (max_anisotropy=1) and/or the slope amplitude
/// boost (amplitude_slope_factor=0) — does the border streaking vanish?
#[test]
#[ignore]
fn export_border_fbm_pin() {
    let dir = output_dir().join("border_fbm");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso_config = IsostasyConfig::c1_default();
    let seed = 1988u64;
    let mut state = init_c1_state_phase_2_r7(64, seed, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    let closures = C1Closures::default();
    let config = C1TimeLoopConfig {
        rigid_continental_crust: true, n_steps: 300,
        dx: 1.0 / 64.0, dy: 1.0 / 64.0,
        iso_config: iso_config.clone(), drainage_max_distance: 30,
    };
    run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
    // Match the fbm2048 case the artefact was reported in (no coast warp,
    // default amplitude); vary only the slope-driven FBM terms. 1024² for
    // viewability.
    let base = FbmUpscaleConfig {
        target_size: 1024, coast_warp_strength: 0.0, amplitude_base: 0.08,
        ..Default::default()
    };
    let variants: [(&str, f64, f64); 4] = [
        // tag, max_anisotropy, amplitude_slope_factor
        ("default", 3.0, 3.0),
        ("aniso1", 1.0, 3.0),     // isotropic → streaks gone?
        ("slopefac0", 3.0, 0.0),  // no slope amplitude boost
        ("both", 1.0, 0.0),       // both off
    ];
    eprintln!("#151 border-FBM pin (seed {seed}) — aniso/slope-boost counterfactual");
    for (tag, aniso, sfac) in &variants {
        let cfg = FbmUpscaleConfig {
            max_anisotropy: *aniso, amplitude_slope_factor: *sfac, ..base.clone()
        };
        let up = upscale_from_c1(
            &state, &iso_config, &closures.oceanic_bathymetry, &WorldSeed::new(seed), &cfg,
        );
        save_heightmap01(&up.heightmap, &dir.join(format!("border_{tag}.png")));
        eprintln!("  {tag}: aniso={aniso} slope_factor={sfac}");
    }
    eprintln!("  out = {}", dir.display());
}

/// #151 v2 re-validation — render the upscale at v2's DEFAULT workflow
/// target (2048²) with DEFAULT config (no coast warp) to isolate the FBM
/// frequency change. Run before vs after the FBM recalibration (swap
/// upscale.rs) to judge v2's calibration at its real target.
#[test]
#[ignore]
fn export_fbm_2048_isolate() {
    let dir = output_dir().join("fbm_2048");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso_config = IsostasyConfig::c1_default();
    let seed = 1988u64;
    let mut state = init_c1_state_phase_2_r7(64, seed, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    let closures = C1Closures::default();
    let config = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: iso_config.clone(),
        drainage_max_distance: 30,
    };
    run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
    // Default config (coast_warp off) at v2's default target 2048².
    // FBM_ANISO env overrides max_anisotropy (counterfactual: 1 = isotropic).
    let aniso = std::env::var("FBM_ANISO").ok().and_then(|s| s.parse().ok()).unwrap_or(3.0);
    let sfac = std::env::var("FBM_SLOPEFAC").ok().and_then(|s| s.parse().ok()).unwrap_or(3.0);
    let amp = std::env::var("FBM_AMP").ok().and_then(|s| s.parse().ok()).unwrap_or(0.08);
    let band = std::env::var("FBM_BAND").ok().and_then(|s| s.parse().ok()).unwrap_or(0.0);
    let subdamp = std::env::var("FBM_SUBDAMP").ok().and_then(|s| s.parse().ok()).unwrap_or(0.3);
    let cfg = FbmUpscaleConfig {
        target_size: 2048, max_anisotropy: aniso, amplitude_slope_factor: sfac,
        amplitude_base: amp, coastal_amplitude_band: band, submarine_damping: subdamp,
        ..Default::default()
    };
    let up = upscale_from_c1(
        &state, &iso_config, &closures.oceanic_bathymetry, &WorldSeed::new(seed), &cfg,
    );
    // Filename tag from an env var so before/after runs don't clobber.
    let tag = std::env::var("FBM_TAG").unwrap_or_else(|_| "x".into());
    save_heightmap01(&up.heightmap, &dir.join(format!("fbm2048_{tag}.png")));
    // Native-resolution crops (no downscale) so fine anisotropic streaks
    // are visible. 384² windows over the bridge/peninsula + lower interior.
    for (cname, x0, y0) in [("bridge", 820usize, 760usize), ("lower", 520, 1180), ("lowedge", 640, 1560)] {
        save_heightmap01_crop(&up.heightmap, x0, y0, 384, &dir.join(format!("crop_{cname}_{tag}.png")));
    }
    eprintln!("  fbm2048_{tag}: {}² (+3 native crops)", up.heightmap.width);
}

/// #151 coastline warp at the TARGET resolution (4096²) — see the
/// product the eye will actually judge. Seed 1988, off / 0.8 / 1.2 coast
/// warp (isolated), plus 0.8 + raised amplitude (full-product look).
#[test]
#[ignore]
fn export_coast_warp_4096() {
    let dir = output_dir().join("coast_warp_4096");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso_config = IsostasyConfig::c1_default();
    let seed = 1988u64;
    let mut state = init_c1_state_phase_2_r7(64, seed, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    let closures = C1Closures::default();
    let config = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: iso_config.clone(),
        drainage_max_distance: 30,
    };
    run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});

    let variants: [(&str, f64, f64); 4] = [
        // tag, coast_warp_strength, amplitude_base
        ("off", 0.0, 0.08),
        ("c08", 0.8, 0.08),
        ("c12", 1.2, 0.08),
        ("c08_amp", 0.8, 0.16),
    ];
    eprintln!("#151 coast warp @ 4096² (seed {seed}) — 0.8 coarse cell ≈ 51 px at this res");
    for (tag, strength, amp) in &variants {
        let cfg = FbmUpscaleConfig {
            target_size: 4096,
            coast_warp_strength: *strength,
            amplitude_base: *amp,
            ..Default::default()
        };
        let t0 = std::time::Instant::now();
        let up = upscale_from_c1(
            &state, &iso_config, &closures.oceanic_bathymetry, &WorldSeed::new(seed), &cfg,
        );
        save_heightmap01(&up.heightmap, &dir.join(format!("coast4096_{tag}.png")));
        eprintln!("  [{tag:>7}] {}² in {:.2?}", up.heightmap.width, t0.elapsed());
    }
    eprintln!("  out = {}", dir.display());
}

/// #151 coastline-warp variants — does displacing the coarse-altitude
/// sampling (`coast_warp_strength`, in coarse cells) break the blocky
/// 64² coastline (STEP-1 fix)? Eye-judged: credible meander vs
/// procedural-fake regular ripples. All via the contract `upscale_from_c1`.
#[test]
#[ignore]
fn export_coast_warp() {
    let dir = output_dir().join("coast_warp");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso_config = IsostasyConfig::c1_default();
    let seeds: [u64; 2] = [42, 1988];
    // Isolate the coast warp: amplitude/damping at defaults; vary only
    // coast_warp_strength (coarse cells). 0.0 = baseline (blocky).
    let strengths: [(&str, f64); 4] = [
        ("off", 0.0),
        ("c05", 0.5),
        ("c08", 0.8),
        ("c12", 1.2),
    ];
    eprintln!("#151 coastline warp — seeds {seeds:?}, coast_warp_strength sweep (coarse cells)");
    for &seed in &seeds {
        let mut state = init_c1_state_phase_2_r7(64, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true,
            n_steps: 300,
            dx: 1.0 / 64.0,
            dy: 1.0 / 64.0,
            iso_config: iso_config.clone(),
            drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        for (tag, strength) in &strengths {
            let cfg = FbmUpscaleConfig {
                target_size: 1024,
                coast_warp_strength: *strength,
                ..Default::default()
            };
            let up = upscale_from_c1(
                &state, &iso_config, &closures.oceanic_bathymetry, &WorldSeed::new(seed), &cfg,
            );
            save_heightmap01(&up.heightmap, &dir.join(format!("coast_seed{seed:05}_{tag}.png")));
        }
        eprintln!("  seed {seed} done");
    }
    eprintln!("  out = {}", dir.display());
    eprintln!("  EYE: coast meanders credibly (real-coast irregular) vs procedural-fake ripples?");
}

/// Quick HD EXPORT — see the final product (no UI). Runs C1 production
/// (64², rigid, full closures) then `upscale_from_c1` (THE contract
/// function — laundered altitude, NOT raw S̃) → 1024² HD heightmap PNG,
/// for several validated seeds (42, 1988 "best relief", 4138 "oceanic
/// world"). The eye on these decides what (if anything) to enrich next
/// (Phase 3 tectonic morpho / sculpting chantier / nothing).
#[test]
#[ignore]
fn export_hd_upscaled() {
    let dir = output_dir().join("hd_export");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso_config = IsostasyConfig::c1_default();
    // Two upscale configs to SEE the effect of the upscale's own
    // pattern-breaking knobs (default leaves them OFF → blocky coast):
    //  - "default": stock FbmUpscaleConfig (domain_warp_strength = 0.0).
    //  - "warped":  domain warp ON + stronger amplitude — the upscale's
    //    documented "break regular patterns" path. Measure if it
    //    dissolves the coarse 64² coastline blockiness.
    // #147 coastline cause-pinning (contrefactuel A): isolate
    // submarine_damping. If dropping it (0.3→0.0) de-blockifies the
    // coast → STEP 2 (damping suppressed a coast-breaking FBM). If not
    // → STEP 1 (the coast contour follows the bilinear-interpolated 64²
    // altitude; FBM height can't move the level contour). `no_damp_amp`
    // also raises amplitude so "no change" can't be blamed on FBM being
    // too weak to reveal the damping effect.
    let configs: [(&str, FbmUpscaleConfig); 3] = [
        ("default", FbmUpscaleConfig { target_size: 1024, ..Default::default() }),
        (
            "no_damp",
            FbmUpscaleConfig { target_size: 1024, submarine_damping: 0.0, ..Default::default() },
        ),
        (
            "no_damp_amp",
            FbmUpscaleConfig {
                target_size: 1024,
                submarine_damping: 0.0,
                amplitude_base: 0.16,
                ..Default::default()
            },
        ),
    ];
    let seeds: [u64; 3] = [42, 1988, 4138];

    eprintln!("HD export — C1 64² production → upscale_from_c1 (contract) → 1024² PNG");
    for &seed in &seeds {
        let mut state = init_c1_state_phase_2_r7(64, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true,
            n_steps: 300,
            dx: 1.0 / 64.0,
            dy: 1.0 / 64.0,
            iso_config: iso_config.clone(),
            drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});

        for (tag, cfg) in &configs {
            let t0 = std::time::Instant::now();
            let up = upscale_from_c1(
                &state,
                &iso_config,
                &closures.oceanic_bathymetry,
                &WorldSeed::new(seed),
                cfg,
            );
            let dt = t0.elapsed();
            save_heightmap01(
                &up.heightmap,
                &dir.join(format!("hd_seed{seed:05}_{tag}_1024.png")),
            );
            eprintln!("  seed {seed:>5} [{tag:>7}]: {}² HD in {:.2?}", up.heightmap.width, dt);
        }
    }
    eprintln!("  out = {} (cost note: 1024² ~per above; 4096² ≈ 16× the FBM)", dir.display());
}

/// Save a native-resolution `size×size` crop of a `[0,1]` heightmap at
/// `(x0,y0)` (1:1, no scaling) so fine detail (anisotropic streaks) is
/// visible without Read's downscaling.
fn save_heightmap01_crop(h: &GridF32, x0: usize, y0: usize, size: usize, path: &Path) {
    let w = (x0 + size).min(h.width) - x0;
    let ht = (y0 + size).min(h.height) - y0;
    let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(w as u32, ht as u32);
    for jj in 0..ht {
        for ii in 0..w {
            let v = h.get((x0 + ii) as i32, (y0 + jj) as i32).clamp(0.0, 1.0);
            img.put_pixel(ii as u32, (ht - 1 - jj) as u32, Rgb(hypsometric(v, 0.5)));
        }
    }
    img.save(path).expect("save crop PNG");
}

/// Render a normalised `[0,1]` heightmap (sea at 0.5) with the
/// hypsometric palette, 1:1 (no upscale).
fn save_heightmap01(h: &GridF32, path: &Path) {
    let (nx, ny) = (h.width, h.height);
    let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(nx as u32, ny as u32);
    for j in 0..ny {
        for i in 0..nx {
            let v = h.get(i as i32, j as i32).clamp(0.0, 1.0);
            img.put_pixel(i as u32, (ny - 1 - j) as u32, Rgb(hypsometric(v, 0.5)));
        }
    }
    img.save(path).expect("save heightmap01 PNG");
}

/// #147 (C) — visual characterisation of the advection-only
/// decorrelation: is it BULK (interior differs across resolutions too →
/// upwind transport mesh-dependent everywhere) or BOUNDARY-only
/// (interior coherent, only the rim differs → not bulk). Dumps S̃ for
/// advection-only (rigid no-flux, NO closures) at 64² and 256², scaled
/// to the same ~512 px display so the SAME physical area is compared.
#[test]
#[ignore]
fn advection_only_visual_64_256() {
    let dir = output_dir().join("advection_only");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso_config = IsostasyConfig::c1_default();
    let seed = 42u64;
    for &grid in &[64usize, 256] {
        let n_steps = 300 * grid / 64;
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true,
            n_steps,
            dx: 1.0 / grid as f64,
            dy: 1.0 / grid as f64,
            iso_config: iso_config.clone(),
            drainage_max_distance: 30,
        };
        ymir_core::tectonics_c1::time_loop::run_advection_only(
            &mut state, &kin, &config, |_, _| {},
        );
        let scale = (512 / grid).max(1) as u32;
        save_s_scaled(&state.s, &dir.join(format!("advonly_grid{grid:04}_s.png")), scale);
        // Low-range palette (S̃∈[0,0.6] → full ramp) so the OCEAN field
        // (≈0.2 ± advected structure, saturated to flat blue in the
        // [0,3] palette) is visible — that's where advection acts and
        // the decorrelation must live.
        save_s_lowrange(&state.s, &dir.join(format!("advonly_grid{grid:04}_s_ocean.png")), scale);
    }
    eprintln!("advection-only S̃ dumped to {}", dir.display());
}

/// #147 (C-bis) — is the bulk decorrelation in the ADVECTION or already
/// in the INIT? Correlate S̃ at step 0 (no stepping at all) across
/// resolutions, whole-field + continental-only + oceanic-only. If init
/// S̃→64 r is ALREADY ~0.045, the non-convergence is the GRID-DEPENDENT
/// INITIAL CONDITION (R7 per-cell heterogeneity), not the upwind scheme.
#[test]
#[ignore]
fn init_convergence_check() {
    let seed = 42u64;
    let grids: [usize; 3] = [64, 128, 256];
    eprintln!("#147 (C-bis) — INIT S̃ convergence (step 0, no advection, no closures)");
    eprintln!("  {:>5} | {:>9} {:>9} {:>9}", "grid", "all r", "cont r", "ocean r");
    let (mut ref_all, mut ref_c, mut ref_o): (Vec<f64>, Vec<f64>, Vec<f64>) =
        (vec![], vec![], vec![]);
    for &grid in grids.iter() {
        let state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let cont: Vec<bool> = (0..grid * grid)
            .map(|k| matches!(state.plate_type.get(k % grid, k / grid), PlateType::Continental))
            .collect();
        let all = block_mean_to_64(&|i, j| state.s.get(i, j), grid);
        // Continental / oceanic masked fields (non-member cells → 0;
        // correlation still dominated by the masked region's structure).
        let cmask = block_mean_to_64(
            &|i, j| if cont[j * grid + i] { state.s.get(i, j) } else { 0.0 },
            grid,
        );
        let omask = block_mean_to_64(
            &|i, j| if !cont[j * grid + i] { state.s.get(i, j) } else { 0.0 },
            grid,
        );
        let (ra, rc, ro) = if grid == 64 {
            ref_all = all.clone();
            ref_c = cmask.clone();
            ref_o = omask.clone();
            (1.0, 1.0, 1.0)
        } else {
            (pearson(&ref_all, &all), pearson(&ref_c, &cmask), pearson(&ref_o, &omask))
        };
        eprintln!("  {:>4}² | {ra:>9.4} {rc:>9.4} {ro:>9.4}", grid);
    }
    eprintln!();
    eprintln!("  init all r ~0.045 ⇒ NON-CONVERGENCE IS THE INIT (grid-dependent R7), not advection.");
}

/// Render S̃ with a stretched low range (`[0, 0.6]` → full ramp) so
/// sub-continental (ocean) structure is visible. Grayscale-ish ramp.
fn save_s_lowrange(s: &Field2D, path: &Path, scale: u32) {
    let (nx, ny) = (s.nx(), s.ny());
    let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(nx as u32 * scale, ny as u32 * scale);
    for j in 0..ny {
        for i in 0..nx {
            let v = (s.get(i, j) / 0.6).clamp(0.0, 1.0) as f32;
            let g = (v * 255.0) as u8;
            put_block(&mut img, i, ny - 1 - j, scale, [g, g, g]);
        }
    }
    img.save(path).expect("save low-range S̃ PNG");
}

/// Cells within `r` cells of a continental/oceanic boundary (the band
/// the contrast lives on). A boundary cell is a `cont` cell 4-adjacent
/// to non-`cont` or vice versa; dilated by `r`.
fn boundary_band(cont: &[bool], grid: usize, r: usize) -> Vec<bool> {
    let idx = |i: usize, j: usize| j * grid + i;
    let mut seed = vec![false; grid * grid];
    for j in 0..grid {
        for i in 0..grid {
            let c = cont[idx(i, j)];
            for (di, dj) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
                let (ni, nj) = (i as i32 + di, j as i32 + dj);
                if ni >= 0 && nj >= 0 && ni < grid as i32 && nj < grid as i32
                    && cont[idx(ni as usize, nj as usize)] != c
                {
                    seed[idx(i, j)] = true;
                }
            }
        }
    }
    dilate(&seed, grid, r)
}

/// Jacobi smoothing of `s` toward the in-bounds 4-neighbour mean, only
/// on `band` cells, blend factor `lambda`.
fn smooth_band(s: &mut Field2D, band: &[bool], grid: usize, lambda: f64) {
    let old: Vec<f64> = (0..grid * grid).map(|k| s.get(k % grid, k / grid)).collect();
    let idx = |i: usize, j: usize| j * grid + i;
    for j in 0..grid {
        for i in 0..grid {
            if !band[idx(i, j)] {
                continue;
            }
            let (mut sum, mut n) = (0.0, 0.0);
            for (di, dj) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
                let (ni, nj) = (i as i32 + di, j as i32 + dj);
                if ni >= 0 && nj >= 0 && ni < grid as i32 && nj < grid as i32 {
                    sum += old[idx(ni as usize, nj as usize)];
                    n += 1.0;
                }
            }
            let mean = sum / n;
            s.set(i, j, (1.0 - lambda) * old[idx(i, j)] + lambda * mean);
        }
    }
}

/// #147 — counterfactual ATTRIBUTION sweep. Decompose the residual
/// mesh non-convergence (post-Fix-#1: S̃→64 r ~0.47, wedge% ∝1/grid)
/// between the oceanic ACCRETION margin pile and the no-flux CURTAIN,
/// using existing closure toggles only (NO fix code):
///   A. full (all ON, post-Fix-#1) — reference.
///   B. subduction+accretion OFF — isolate the accretion margin pile.
///   C. all closures OFF (advection + rigid no-flux only) — isolate
///      the curtain (the interior speckle is an advection+no-flux
///      artefact on the sharp 1.0/0.2 contrast, present without any
///      closure).
/// Per (variant, grid): S̃→64 r (correlated to that VARIANT's OWN 64²
/// reference) + wedge% (S̃>1.5). Reads the attribution off how each
/// metric moves when a mechanism is removed; coupling shows up as one
/// toggle moving the OTHER mechanism's metric.
#[test]
#[ignore]
fn mesh_convergence_attribution() {
    let iso_config = IsostasyConfig::c1_default();
    let seed = 42u64;
    let grids: [usize; 3] = [64, 128, 256];

    struct Variant {
        tag: &'static str,
        mutate: fn(&mut C1Closures),
    }
    let variants = [
        Variant { tag: "A full", mutate: |_| {} },
        // Split B (#147 step 1) — which of the Track-D pair carries
        // wedge%? Don't assume accretion (DS bet just lost); subduction
        // promotion is per-cell too (#145 finger / ≥2 / Oceanic→Continental).
        Variant { tag: "B1 sub_off", mutate: |c| c.subduction.enabled = false },
        Variant { tag: "B2 acc_off", mutate: |c| c.accretion.enabled = false },
        Variant {
            tag: "B sub+acc_off",
            mutate: |c| {
                c.subduction.enabled = false;
                c.accretion.enabled = false;
            },
        },
        Variant {
            tag: "C advection_only",
            mutate: |c| {
                c.davis_suppe.enabled = false;
                c.equilibrium_height.enabled = false;
                c.erosion.enabled = false;
                c.oceanic_bathymetry.enabled = false;
                c.subduction.enabled = false;
                c.accretion.enabled = false;
                c.rifting.enabled = false;
            },
        },
    ];

    eprintln!("#147 — counterfactual ATTRIBUTION sweep (seed {seed}, rigid, post-Fix-#1)");
    eprintln!("  isolate accretion pile (B) vs curtain (C) in the residual non-convergence");
    eprintln!(
        "  {:<18} {:>5} | {:>8} {:>9}",
        "variant", "grid", "wedge%", "S̃→64 r"
    );

    for v in &variants {
        let mut ref64: Vec<f64> = Vec::new();
        for &grid in grids.iter() {
            let mut closures = C1Closures::default();
            (v.mutate)(&mut closures);
            let n_steps = 300 * grid / 64;
            let mut state =
                init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
            let mut kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
            let config = C1TimeLoopConfig {
                rigid_continental_crust: true,
                n_steps,
                dx: 1.0 / grid as f64,
                dy: 1.0 / grid as f64,
                iso_config: iso_config.clone(),
                drainage_max_distance: 30,
            };
            run_with_closures(&mut state, &mut kinematics, &config, &closures, |_, _| {});

            let wedge_n = (0..grid * grid)
                .filter(|&k| {
                    let (i, j) = (k % grid, k / grid);
                    state.s.get(i, j) > 1.5
                })
                .count();
            let wedge_pct = 100.0 * wedge_n as f64 / (grid * grid) as f64;

            let s_ds = block_mean_to_64(&|i, j| state.s.get(i, j), grid);
            let r = if grid == 64 {
                ref64 = s_ds.clone();
                1.0
            } else {
                pearson(&ref64, &s_ds)
            };

            eprintln!(
                "  {:<18} {:>4}² | {:>8.2} {:>9.4}",
                if grid == 64 { v.tag } else { "" },
                grid,
                wedge_pct,
                r
            );
        }
    }
    eprintln!();
    eprintln!("  Read: wedge% stabilises in B → accretion carries wedge%.");
    eprintln!("        S̃ r climbs in C → curtain carries field decorrelation.");
    eprintln!("        a toggle moving the OTHER metric → coupling (one fix may do both).");
}

/// BFS distance (in cells) to the coast over continental cells.
/// Coast = continental cell 4-adjacent to a non-continental cell or
/// grid edge. Non-continental cells get `usize::MAX`.
fn dist_to_coast(cont: &[bool], grid: usize) -> Vec<usize> {
    let mut dist = vec![usize::MAX; grid * grid];
    let mut q = std::collections::VecDeque::new();
    let idx = |i: usize, j: usize| j * grid + i;
    for j in 0..grid {
        for i in 0..grid {
            if !cont[idx(i, j)] {
                continue;
            }
            let mut coast = i == 0 || j == 0 || i == grid - 1 || j == grid - 1;
            for (di, dj) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
                let (ni, nj) = (i as i32 + di, j as i32 + dj);
                if ni >= 0 && nj >= 0 && ni < grid as i32 && nj < grid as i32
                    && !cont[idx(ni as usize, nj as usize)]
                {
                    coast = true;
                }
            }
            if coast {
                dist[idx(i, j)] = 0;
                q.push_back((i, j));
            }
        }
    }
    while let Some((i, j)) = q.pop_front() {
        let d = dist[idx(i, j)];
        for (di, dj) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
            let (ni, nj) = (i as i32 + di, j as i32 + dj);
            if ni >= 0 && nj >= 0 && ni < grid as i32 && nj < grid as i32 {
                let nk = idx(ni as usize, nj as usize);
                if cont[nk] && dist[nk] == usize::MAX {
                    dist[nk] = d + 1;
                    q.push_back((ni as usize, nj as usize));
                }
            }
        }
    }
    dist
}

fn m_components(mask: &[bool], grid: usize) -> usize {
    land_morphology(mask, grid, grid).n_components
}

/// Block-average a `grid×grid` field down to 64×64 (grid must be a
/// multiple of 64). Returns row-major 64×64. The downsampled field is
/// the high-res structure as seen at the 64² scale — high correlation
/// with the native 64² field ⇒ the LARGE-SCALE structure is
/// mesh-convergent.
fn block_mean_to_64(get: &dyn Fn(usize, usize) -> f64, grid: usize) -> Vec<f64> {
    let b = grid / 64;
    let mut out = vec![0.0; 64 * 64];
    for cj in 0..64 {
        for ci in 0..64 {
            let mut s = 0.0;
            for dj in 0..b {
                for di in 0..b {
                    s += get(ci * b + di, cj * b + dj);
                }
            }
            out[cj * 64 + ci] = s / (b * b) as f64;
        }
    }
    out
}

/// Pearson correlation between two equal-length samples.
fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let (mut cov, mut va, mut vb) = (0.0, 0.0, 0.0);
    for k in 0..a.len() {
        let (da, db) = (a[k] - ma, b[k] - mb);
        cov += da * db;
        va += da * da;
        vb += db * db;
    }
    if va > 1e-12 && vb > 1e-12 { cov / (va.sqrt() * vb.sqrt()) } else { 0.0 }
}

/// Population std and (p95 − p05) of a sample (sorts in place).
fn std_and_spread(v: &mut [f64]) -> (f64, f64) {
    if v.is_empty() {
        return (0.0, 0.0);
    }
    let mean = v.iter().sum::<f64>() / v.len() as f64;
    let var = v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / v.len() as f64;
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let pct = |p: f64| v[((p * (v.len() - 1) as f64).round() as usize).min(v.len() - 1)];
    (var.sqrt(), pct(0.95) - pct(0.05))
}

fn render_altitude(state: &C1State, iso: &IsostasyConfig, closures: &C1Closures) -> GridF32 {
    let isostasy = compute_isostasy(&state.s, iso);
    let mut altitude = isostasy.heightmap;
    apply_stein_stein_bathymetry(
        &mut altitude,
        &state.age,
        &state.plate_type,
        &closures.oceanic_bathymetry,
    );
    altitude
}

fn dump_altitude(
    state: &C1State,
    cycle: usize,
    dir: &Path,
    prefix: &str,
    iso: &IsostasyConfig,
    closures: &C1Closures,
) {
    let altitude = render_altitude(state, iso, closures);
    save_altitude(&altitude, &dir.join(format!("{}_cycle{:03}_altitude.png", prefix, cycle)));
}

/// Default upscale factor for 64² PNGs (nearest-neighbour) so fine
/// morphology is legible by eye. The resolution diagnostic overrides
/// this per-grid to keep the display size ~512 px.
const SCALE: u32 = 8;

fn save_altitude(altitude: &GridF32, path: &Path) {
    save_altitude_scaled(altitude, path, SCALE);
}
fn save_s(s: &Field2D, path: &Path) {
    save_s_scaled(s, path, SCALE);
}

fn put_block(img: &mut ImageBuffer<Rgb<u8>, Vec<u8>>, i: usize, row: usize, scale: u32, rgb: [u8; 3]) {
    for dy in 0..scale {
        for dx in 0..scale {
            img.put_pixel(i as u32 * scale + dx, row as u32 * scale + dy, Rgb(rgb));
        }
    }
}

fn save_altitude_scaled(altitude: &GridF32, path: &Path, scale: u32) {
    let (nx, ny) = (altitude.width, altitude.height);
    let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(nx as u32 * scale, ny as u32 * scale);
    let sea_norm = 0.5_f32;
    for j in 0..ny {
        for i in 0..nx {
            let raw = altitude.get(i as i32, j as i32);
            let t = ((raw + ALT_HALF) / (2.0 * ALT_HALF)).clamp(0.0, 1.0);
            put_block(&mut img, i, ny - 1 - j, scale, hypsometric(t, sea_norm));
        }
    }
    img.save(path).expect("save altitude PNG");
}

fn save_s_scaled(s: &Field2D, path: &Path, scale: u32) {
    let (nx, ny) = (s.nx(), s.ny());
    let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(nx as u32 * scale, ny as u32 * scale);
    let sea_norm = (0.2 / S_VIZ_MAX) as f32;
    for j in 0..ny {
        for i in 0..nx {
            let v = (s.get(i, j) / S_VIZ_MAX).clamp(0.0, 1.0) as f32;
            put_block(&mut img, i, ny - 1 - j, scale, hypsometric(v, sea_norm));
        }
    }
    img.save(path).expect("save S̃ PNG");
}

fn hypsometric(h: f32, sea_norm: f32) -> [u8; 3] {
    let mid = (sea_norm + 1.0) * 0.5;
    let lerp = |t: f32, a: [u8; 3], b: [u8; 3]| -> [u8; 3] {
        let t = t.clamp(0.0, 1.0);
        [
            (a[0] as f32 + t * (b[0] as f32 - a[0] as f32)).round() as u8,
            (a[1] as f32 + t * (b[1] as f32 - a[1] as f32)).round() as u8,
            (a[2] as f32 + t * (b[2] as f32 - a[2] as f32)).round() as u8,
        ]
    };
    if h <= sea_norm * 0.5 {
        lerp(h / (sea_norm * 0.5).max(1e-6), [10, 20, 60], [40, 80, 160])
    } else if h <= sea_norm {
        lerp((h - sea_norm * 0.5) / (sea_norm * 0.5).max(1e-6), [40, 80, 160], [120, 180, 230])
    } else if h <= mid {
        lerp((h - sea_norm) / (mid - sea_norm).max(1e-6), [60, 130, 60], [140, 100, 50])
    } else {
        lerp((h - mid) / (1.0 - mid).max(1e-6), [140, 100, 50], [245, 245, 245])
    }
}

/// #155 amplitude reading (a) — does the 1b-i O-C margin ridge read as a
/// MOUNTAIN CHAIN after the HD upscale (upscale_from_c1 = isostasy +
/// Stein-Stein + bicubic + FBM detail oriented by the altitude gradient),
/// or stay a low soft band? Hillshade (palette-free) of the HD 1024²
/// product on seeds 42 (S arc) + 1988 (N+E margins) — the franc-O-C seeds.
/// Hypothesis: 1b-i turned the dome's MOU gradient into a margin-peaked
/// NET gradient; the FBM amplifies gradients → a low-but-sharp ridge may
/// become a chain in HD where the dome could not. Measure before coding
/// any force/meso fix (reading a before b/c).
#[test]
#[ignore]
fn measure_ridge_hd_amplitude() {
    let dir = output_dir().join("ridge_hd_amplitude");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso_config = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let cfg = FbmUpscaleConfig {
        target_size: 1024, coast_warp_strength: 1.5, coast_warp_frequency: 0.5,
        coastal_amplitude_band: 0.30, amplitude_base: 0.16, submarine_damping: 0.0,
        ..Default::default()
    };
    eprintln!("#155 (a) ridge HD amplitude — hillshade of upscale_from_c1 (1024²)");
    for &seed in &[42u64, 1988u64] {
        let mut state = init_c1_state_phase_2_r7(64, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300,
            dx: 1.0 / 64.0, dy: 1.0 / 64.0,
            iso_config: iso_config.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        let up = upscale_from_c1(&state, &iso_config, &ss, &WorldSeed::new(seed), &cfg);
        let h = &up.heightmap;
        // Full-frame hillshade (palette-free relief).
        save_hillshade_crop(h, 0, 0, h.width, &dir.join(format!("seed{seed:05}_hillshade_full.png")));
        // Also grayscale full (continuous height, no palette bands).
        save_gray01_crop(h, 0, 0, h.width, &dir.join(format!("seed{seed:05}_gray_full.png")));
        eprintln!("  seed {seed}: HD {}² hillshade + gray written", h.width);
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 amplitude reading — quantitative O-C-margin vs passive-margin
/// contrast. Does 1b-i leave a measurable elevation signal at O-C
/// margins above passive (non-convergent) margins, or is it ≈1 (coastal
/// band erases it)? Measured on the COARSE production altitude (where
/// 1b-i imprints directly, before FBM dilution). Bands: cells 1..=6
/// inland of their respective margin. O-C band = typed-wedge is_oc &&
/// 1<=d<=6; passive band = continental, !is_oc, Manhattan dist-to-ocean
/// 1..=6. Ratio >1 (even modest) = sub-visible signal erosion could
/// amplify; ≈1 = 1b-i imprints ~nothing in elevation.
#[test]
#[ignore]
fn measure_oc_vs_passive_margin() {
    use ymir_core::tectonics_c1::boundary_classification::{
        classify_boundaries, oc_override_seed_mask, retarget_upper_plate_continental,
    };
    use ymir_core::tectonics_c1::distance_field::wedge_distance_intra_plate_typed;
    let iso_config = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    eprintln!("#155 O-C vs passive margin elevation contrast (coarse altitude, band 1..=6)");
    for &seed in &[42u64, 1988u64, 1337u64, 2u64] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300,
            dx: 1.0 / 64.0, dy: 1.0 / 64.0,
            iso_config: iso_config.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});

        // Reproduce production geometry on the final state.
        let mut bi = classify_boundaries(&state.plate_id, &kin);
        retarget_upper_plate_continental(&mut bi, &state.plate_id, &state.plate_type);
        let oc_seed = oc_override_seed_mask(&bi, &state.plate_id, &state.plate_type);
        let (wd, is_oc) = wedge_distance_intra_plate_typed(
            &state.plate_id, &bi.upper_plate_mask, &oc_seed, 30.0,
        );
        let alt = c1_production_altitude(&state.s, &state.age, &state.plate_type, &iso_config, &ss);

        let is_cont = |i: usize, j: usize| matches!(state.plate_type.get(i, j), PlateType::Continental);
        // Manhattan dist-to-ocean over continental cells (multi-source BFS).
        let mut d_ocean = vec![usize::MAX; grid * grid];
        let mut q = std::collections::VecDeque::new();
        for j in 0..grid { for i in 0..grid {
            if !is_cont(i, j) { d_ocean[j * grid + i] = 0; q.push_back((i, j)); }
        }}
        while let Some((i, j)) = q.pop_front() {
            let dc = d_ocean[j * grid + i];
            for (di, dj) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
                let (ni, nj) = (i as i32 + di, j as i32 + dj);
                if ni < 0 || nj < 0 || ni as usize >= grid || nj as usize >= grid { continue; }
                let (ni, nj) = (ni as usize, nj as usize);
                if d_ocean[nj * grid + ni] > dc + 1 {
                    d_ocean[nj * grid + ni] = dc + 1;
                    q.push_back((ni, nj));
                }
            }
        }

        // ARTIFACT GUARD (see [[feedback_ratio_across_affine_rescaled_spaces]]):
        // report the O-C−passive GAP in the SAME normalised [0,1] space the
        // HD product uses ((alt+half)/(2·half), sea→0.5), NOT a ratio on raw
        // altitude. A raw-altitude ratio here (1.2-1.5) vs an HD-normalised
        // ratio (1.05-1.12) elsewhere is a false "collapse" (pure +offset
        // arithmetic) — it misled the #155 amplitude diagnostic for a full
        // round. Gap-in-consistent-space is cross-state comparable; ratio is
        // not. `half` mirrors production_upscale::ALTITUDE_NORM_HALF_RANGE.
        let half = 1.13_f64;
        let norm = |a: f64| ((a + half) / (2.0 * half)).clamp(0.0, 1.0);
        let (mut oc_sum, mut oc_n) = (0.0f64, 0usize);
        let (mut pa_sum, mut pa_n) = (0.0f64, 0usize);
        for j in 0..grid { for i in 0..grid {
            if !is_cont(i, j) { continue; }
            let a = norm(alt.get(i as i32, j as i32) as f64);
            let dw = wd.get(i, j);
            if is_oc.get(i, j) && dw >= 1.0 && dw <= 6.0 {
                oc_sum += a; oc_n += 1;
            } else if !is_oc.get(i, j) {
                let do_ = d_ocean[j * grid + i];
                if do_ >= 1 && do_ <= 6 { pa_sum += a; pa_n += 1; }
            }
        }}
        let oc_mean = oc_sum / oc_n.max(1) as f64;
        let pa_mean = pa_sum / pa_n.max(1) as f64;
        eprintln!(
            "  seed {seed:5}: O-C band (norm) = {oc_mean:.4} (n={oc_n:4})  passive band (norm) = {pa_mean:.4} (n={pa_n:4})  GAP = {:+.4}",
            oc_mean - pa_mean
        );
    }
}

/// #155 reading (a-erosion) PROXY — apply v2 droplet hydraulic erosion
/// (erosion::hydraulic::run_erosion — SAME TYPE as the expected phase-4
/// HD erosion, so representative) to the C1 upscaled heightmap, and ask
/// the DIFFERENTIAL question: does erosion amplify the O-C ridge
/// PREFERENTIALLY (O-C/passive elevation ratio RISES → ridge emerges,
/// (a) succeeds) or UNIFORMLY (ratio flat → ridge stays drowned, (a)
/// fails → (b)/(c))? Measures the same O-C-vs-passive contrast on the HD
/// heightmap PRE and POST erosion (coarse band cells mapped to HD
/// blocks), plus a post-erosion hillshade for the eye.
#[test]
#[ignore]
fn measure_aerosion_proxy() {
    use ymir_core::tectonics_c1::boundary_classification::{
        classify_boundaries, oc_override_seed_mask, retarget_upper_plate_continental,
    };
    use ymir_core::tectonics_c1::distance_field::wedge_distance_intra_plate_typed;
    use ymir_core::erosion::hydraulic::{run_erosion, ErosionConfig};
    let dir = output_dir().join("aerosion_proxy");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso_config = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    let cfg = FbmUpscaleConfig {
        target_size: 1024, coast_warp_strength: 1.5, coast_warp_frequency: 0.5,
        coastal_amplitude_band: 0.30, amplitude_base: 0.16, submarine_damping: 0.0,
        ..Default::default()
    };
    let ero_cfg = ErosionConfig { num_droplets: 2_000_000, batch_size: 50_000, ..Default::default() };
    eprintln!("#155 (a-erosion) proxy — droplet hydraulic erosion on C1 HD; O-C/passive ratio pre vs post");
    for &seed in &[42u64, 1988u64, 1337u64] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300,
            dx: 1.0 / 64.0, dy: 1.0 / 64.0,
            iso_config: iso_config.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});

        let mut bi = classify_boundaries(&state.plate_id, &kin);
        retarget_upper_plate_continental(&mut bi, &state.plate_id, &state.plate_type);
        let oc_seed = oc_override_seed_mask(&bi, &state.plate_id, &state.plate_type);
        let (wd, is_oc) = wedge_distance_intra_plate_typed(
            &state.plate_id, &bi.upper_plate_mask, &oc_seed, 30.0,
        );
        let is_cont = |i: usize, j: usize| matches!(state.plate_type.get(i, j), PlateType::Continental);
        // dist-to-ocean (Manhattan) for the passive band.
        let mut d_ocean = vec![usize::MAX; grid * grid];
        let mut q = std::collections::VecDeque::new();
        for j in 0..grid { for i in 0..grid {
            if !is_cont(i, j) { d_ocean[j * grid + i] = 0; q.push_back((i, j)); }
        }}
        while let Some((i, j)) = q.pop_front() {
            let dc = d_ocean[j * grid + i];
            for (di, dj) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
                let (ni, nj) = (i as i32 + di, j as i32 + dj);
                if ni < 0 || nj < 0 || ni as usize >= grid || nj as usize >= grid { continue; }
                let (ni, nj) = (ni as usize, nj as usize);
                if d_ocean[nj * grid + ni] > dc + 1 { d_ocean[nj * grid + ni] = dc + 1; q.push_back((ni, nj)); }
            }
        }
        // Coarse band membership: 0 = none, 1 = O-C, 2 = passive.
        let mut band = vec![0u8; grid * grid];
        for j in 0..grid { for i in 0..grid {
            if !is_cont(i, j) { continue; }
            let dw = wd.get(i, j);
            if is_oc.get(i, j) && dw >= 1.0 && dw <= 6.0 { band[j * grid + i] = 1; }
            else if !is_oc.get(i, j) {
                let do_ = d_ocean[j * grid + i];
                if do_ >= 1 && do_ <= 6 { band[j * grid + i] = 2; }
            }
        }}

        let up = upscale_from_c1(&state, &iso_config, &ss, &WorldSeed::new(seed), &cfg);
        let h0 = up.heightmap.clone();
        let eroded = run_erosion(&h0, &ero_cfg, &WorldSeed::new(seed), |_, _, _| true);
        let h1 = &eroded.heightmap;
        let scale = h0.width / grid; // 16

        // Mean HD height over the coarse band's HD blocks.
        let mean_over_band = |hd: &GridF32, target: u8| -> f64 {
            let (mut s, mut n) = (0.0f64, 0usize);
            for j in 0..grid { for i in 0..grid {
                if band[j * grid + i] != target { continue; }
                for jj in 0..scale { for ii in 0..scale {
                    s += hd.get((i * scale + ii) as i32, (j * scale + jj) as i32) as f64; n += 1;
                }}
            }}
            s / n.max(1) as f64
        };
        // ARTIFACT GUARD ([[feedback_ratio_across_affine_rescaled_spaces]]):
        // h0/h1 are ALREADY normalised [0,1] (upscale_from_c1 output), so
        // report the O-C−passive GAP (cross-state comparable) — NOT a ratio,
        // and NEVER compare this number to a raw-altitude ratio (that mix is
        // the false-collapse trap). The differential verdict is gap_post vs
        // gap_pre: gap RISES → erosion amplifies the ridge preferentially.
        let (oc0, pa0) = (mean_over_band(&h0, 1), mean_over_band(&h0, 2));
        let (oc1, pa1) = (mean_over_band(h1, 1), mean_over_band(h1, 2));
        eprintln!(
            "  seed {seed:5}: PRE  O-C={oc0:.4} passive={pa0:.4} GAP={:+.4}  | POST O-C={oc1:.4} passive={pa1:.4} GAP={:+.4}",
            oc0 - pa0, oc1 - pa1,
        );
        save_hillshade_crop(h1, 0, 0, h1.width, &dir.join(format!("seed{seed:05}_eroded_hillshade.png")));
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 reading (a) CLOSURE — visual confirmation in CONSISTENT
/// normalization. The corrected quantitative verdict: the O-C ridge
/// signal survives upscale + is amplified by erosion (gap +60-76% in
/// consistent space). This triptych confirms by eye whether that
/// preserved-but-modest ridge reads as a distinct CHAIN or stays
/// coastal-dominated (→ méso/structure the next subject).
///
/// CRITICAL: all three states share the SAME normalization AND
/// resolution. The coarse is routed through `upscale_from_c1` with a
/// FLAT config (amplitude_base=0, coast_warp=0 → pure bicubic of the
/// normalized coarse), so it differs from the upscale state ONLY by the
/// FBM detail — not by a normalization offset (the artifact we caught)
/// nor a 64²-vs-1024² hillshade-gradient mismatch (a second cross-space
/// trap). State 2 = prod #151 config; state 3 = state 2 + v2 erosion.
#[test]
#[ignore]
fn triptych_consistent_norm() {
    use ymir_core::erosion::hydraulic::{run_erosion, ErosionConfig};
    let dir = output_dir().join("triptych_norm");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso_config = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    // Flat config: pure bicubic upscale of the normalized coarse (NO FBM,
    // NO warp) → macro ossature at 1024², same normalization as prod.
    let cfg_flat = FbmUpscaleConfig {
        target_size: 1024, amplitude_base: 0.0, coast_warp_strength: 0.0,
        coastal_amplitude_band: 0.0, ..Default::default()
    };
    // Prod #151 config (the real aval).
    let cfg_prod = FbmUpscaleConfig {
        target_size: 1024, coast_warp_strength: 1.5, coast_warp_frequency: 0.5,
        coastal_amplitude_band: 0.30, amplitude_base: 0.16, submarine_damping: 0.0,
        ..Default::default()
    };
    let ero_cfg = ErosionConfig { num_droplets: 2_000_000, batch_size: 50_000, ..Default::default() };
    eprintln!("#155 (a) closure — consistent-norm hillshade triptych (coarse-bicubic / upscale / eroded)");
    for &seed in &[1337u64, 1988u64, 42u64] {
        let mut state = init_c1_state_phase_2_r7(64, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300,
            dx: 1.0 / 64.0, dy: 1.0 / 64.0,
            iso_config: iso_config.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});

        let coarse = upscale_from_c1(&state, &iso_config, &ss, &WorldSeed::new(seed), &cfg_flat);
        let up = upscale_from_c1(&state, &iso_config, &ss, &WorldSeed::new(seed), &cfg_prod);
        let eroded = run_erosion(&up.heightmap, &ero_cfg, &WorldSeed::new(seed), |_, _, _| true);

        let w = coarse.heightmap.width;
        save_hillshade_crop(&coarse.heightmap, 0, 0, w, &dir.join(format!("seed{seed:05}_1_coarse_norm.png")));
        save_hillshade_crop(&up.heightmap, 0, 0, w, &dir.join(format!("seed{seed:05}_2_upscale_norm.png")));
        save_hillshade_crop(&eroded.heightmap, 0, 0, w, &dir.join(format!("seed{seed:05}_3_eroded_norm.png")));
        eprintln!("  seed {seed}: triptych written ({w}²)");
    }
    eprintln!("  out = {}", dir.display());
}
