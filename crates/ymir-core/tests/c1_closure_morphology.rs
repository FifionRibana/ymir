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
use ymir_core::tectonics::isostasy::{compute_isostasy, compute_isostasy_craton, IsostasyConfig};
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

/// #155 MÉSO VIABILITY PROBE (throwaway, before drafting the méso issue
/// or W7-ing the injection — measure, don't deduce). The tectonic-méso
/// approach (deposit parallel ridges in S̃, let erosion orient drainage)
/// is DEDUCED, not measured. Risk: méso structure is FINER than macro;
/// raw S̃ converges only r~0.51 (advection diffuses it) and the upscale
/// reads LAUNDERED altitude precisely because fine S̃ doesn't survive.
/// So synthetic fine ridges in S̃ may be (i) smoothed by upscale, (ii)
/// drowned by FBM, (iii) ignored by erosion (dendritic anyway).
///
/// Inject a SYNTHETIC imbricate stack (cosine in d, wavelength λ) into
/// coarse S̃ at the O-C wedge of a franc-O-C seed, run the FULL aval
/// (upscale + FBM + erosion), and compare to the no-injection baseline:
///   - hillshade (eye): organized parallel ridges + oriented drainage,
///     or smoothed/dendritic?
///   - transverse d-shell profile of HD altitude (normalised): does the
///     injected oscillation (peaks at d≈λ,2λ) survive to HD, or collapse
///     to one smooth bump?
/// Verdict gates #1/#2: structure traverses → tectonic méso VIABLE;
/// smoothed/drowned → méso is an aval-EXPRESSION problem (like S̃ r~0.51),
/// a DIFFERENT chantier — don't draft the tectonic issue.
#[test]
#[ignore]
fn probe_meso_viability() {
    use ymir_core::tectonics_c1::boundary_classification::{
        classify_boundaries, oc_override_seed_mask, retarget_upper_plate_continental,
    };
    use ymir_core::tectonics_c1::distance_field::wedge_distance_intra_plate_typed;
    use ymir_core::erosion::hydraulic::{run_erosion, ErosionConfig};
    let dir = output_dir().join("meso_viability");
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
    let lambda = 3.0_f64; // thrust spacing in coarse cells (~2-3 ridges in d∈0..8)
    let amp = 0.8_f64;     // S̃ modulation amplitude (wedge S̃ ~1-2)
    let band_max = 9.0_f64;
    let half = 1.13_f64;   // production_upscale::ALTITUDE_NORM_HALF_RANGE (consistent norm)
    eprintln!("#155 méso viability — synthetic imbricate stack (λ={lambda} cells, amp={amp}) through the full aval");
    for &seed in &[1988u64, 42u64] {
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

        let _ = half;
        // Aval pass on the current state, tagged. Captures shared refs but
        // NOT `state` (passed as arg), so state.s can be mutated between the
        // two calls (C1State is not Clone → inject in place rather than clone).
        let aval = |st: &C1State, tag: &str| {
            let up = upscale_from_c1(st, &iso_config, &ss, &WorldSeed::new(seed), &cfg);
            let eroded = run_erosion(&up.heightmap, &ero_cfg, &WorldSeed::new(seed), |_, _, _| true);
            save_hillshade_crop(&up.heightmap, 0, 0, up.heightmap.width, &dir.join(format!("seed{seed:05}_{tag}_upscale.png")));
            save_hillshade_crop(&eroded.heightmap, 0, 0, eroded.heightmap.width, &dir.join(format!("seed{seed:05}_{tag}_eroded.png")));
            // Transverse d-shell profile (upscale only, HD already [0,1]):
            // did the injected oscillation (peaks at d≈λ,2λ) survive the
            // upscale, before erosion reorganises?
            let scale = up.heightmap.width / grid;
            let mut shells = vec![(0.0f64, 0usize); 20];
            for j in 0..grid { for i in 0..grid {
                if !is_oc.get(i, j) { continue; }
                let d = wd.get(i, j);
                if d < 0.5 || d > 10.0 { continue; }
                let sh = ((d / 0.5) as usize).min(19);
                let mut s = 0.0; let mut n = 0;
                for jj in 0..scale { for ii in 0..scale {
                    s += up.heightmap.get((i*scale+ii) as i32, (j*scale+jj) as i32) as f64; n += 1;
                }}
                shells[sh].0 += s / n.max(1) as f64; shells[sh].1 += 1;
            }}
            let prof: Vec<String> = shells.iter().enumerate().filter(|(_, c)| c.1 > 0)
                .map(|(k, c)| format!("d{:.1}:{:.3}", k as f64 * 0.5, c.0 / c.1 as f64)).collect();
            eprintln!("  seed {seed} [{tag}] transverse HD profile [0,1]: {}", prof.join(" "));
        };

        aval(&state, "base");
        // Inject the synthetic imbricate stack into S̃ at O-C wedge cells.
        for j in 0..grid { for i in 0..grid {
            if is_oc.get(i, j) {
                let d = wd.get(i, j);
                if d >= 0.5 && d <= band_max {
                    let m = amp * (std::f64::consts::TAU * d / lambda).cos();
                    let v = (state.s.get(i, j) + m).max(0.1);
                    state.s.set(i, j, v);
                }
            }
        }}
        aval(&state, "inj");
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 méso viability LOCALIZE (fast, no upscale/erosion) — where does
/// the injected S̃ structure die? Profiles, per d-shell over O-C cells:
/// (1) raw S̃ (sanity: injection present + oscillating), (2) coarse
/// c1_production_altitude (isostasy(S̃)). If S̃ oscillates but coarse
/// altitude does NOT → isostasy maps the fine structure away (deep). If
/// coarse altitude oscillates but HD (probe_meso_viability) did not →
/// the UPSCALE smooths it (actionable: aval-expression). Decisive for
/// whether méso is tectonic-in-S̃, or an aval-expression chantier.
#[test]
#[ignore]
fn probe_meso_localize() {
    use ymir_core::tectonics_c1::boundary_classification::{
        classify_boundaries, oc_override_seed_mask, retarget_upper_plate_continental,
    };
    use ymir_core::tectonics_c1::distance_field::wedge_distance_intra_plate_typed;
    let iso_config = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    let (lambda, amp, band_max) = (3.0_f64, 0.8_f64, 9.0_f64);
    eprintln!("#155 méso localize — S̃ and COARSE-altitude d-shell profiles, base vs inj");
    for &seed in &[1988u64, 42u64] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0,
            iso_config: iso_config.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        let mut bi = classify_boundaries(&state.plate_id, &kin);
        retarget_upper_plate_continental(&mut bi, &state.plate_id, &state.plate_type);
        let oc_seed = oc_override_seed_mask(&bi, &state.plate_id, &state.plate_type);
        let (wd, is_oc) = wedge_distance_intra_plate_typed(&state.plate_id, &bi.upper_plate_mask, &oc_seed, 30.0);

        let shell_profile = |field: &dyn Fn(usize, usize) -> f64| -> Vec<String> {
            let mut shells = vec![(0.0f64, 0usize); 20];
            for j in 0..grid { for i in 0..grid {
                if !is_oc.get(i, j) { continue; }
                let d = wd.get(i, j);
                if d < 0.5 || d > 10.0 { continue; }
                let sh = ((d / 0.5) as usize).min(19);
                shells[sh].0 += field(i, j); shells[sh].1 += 1;
            }}
            shells.iter().enumerate().filter(|(_, c)| c.1 > 0)
                .map(|(k, c)| format!("d{:.1}:{:.3}", k as f64 * 0.5, c.0 / c.1 as f64)).collect()
        };

        let alt_base = c1_production_altitude(&state.s, &state.age, &state.plate_type, &iso_config, &ss);
        eprintln!("  seed {seed} [base] S̃   : {}", shell_profile(&|i, j| state.s.get(i, j)).join(" "));
        eprintln!("  seed {seed} [base] alt : {}", shell_profile(&|i, j| alt_base.get(i as i32, j as i32) as f64).join(" "));
        // inject
        for j in 0..grid { for i in 0..grid {
            if is_oc.get(i, j) { let d = wd.get(i, j);
                if d >= 0.5 && d <= band_max {
                    let v = (state.s.get(i, j) + amp * (std::f64::consts::TAU * d / lambda).cos()).max(0.1);
                    state.s.set(i, j, v);
                }}
        }}
        let alt_inj = c1_production_altitude(&state.s, &state.age, &state.plate_type, &iso_config, &ss);
        eprintln!("  seed {seed} [inj ] S̃   : {}", shell_profile(&|i, j| state.s.get(i, j)).join(" "));
        eprintln!("  seed {seed} [inj ] alt : {}", shell_profile(&|i, j| alt_inj.get(i as i32, j as i32) as f64).join(" "));
    }
}

/// #155 méso DISTINCTION (fast, no upscale/erosion) — compression (case
/// 1, reversible) vs destruction (case 2/3)? c1_production_altitude's
/// h_raw = S̃·buoyancy is local (structure survives); the [0,1] LAND
/// normalisation divides by raw h_max (line 185), which injected crests
/// self-inflate (the 0.92 cap is applied to sea level ONLY). So the
/// earlier "absorbed" (measured on NORMALISED output) conflates the two.
/// Decisive test: normalise the INJECTED h_raw with the BASE field's
/// fixed range (h_min/h_sea/h_max). If the imbricate oscillation
/// RE-APPEARS under the fixed range → CASE 1 (compression by self-
/// inflated h_max), reversible, structure was never destroyed → the
/// expression fix (percentile-cap the LAND ceiling too) would work and
/// tectonic-méso-in-S̃ is NOT falsified. Measured in CONSISTENT space
/// (same range for both states; gaps, not cross-space ratios).
///
/// ⚠️ GUARD-RAIL — this probe NEARLY MISLED (it over-claimed "+0.587 /
/// cheap land-cap fix"): it hand-rolls the normalisation and OMITS
/// `compute_isostasy`'s gaussian blur (`altitude_smoothing_sigma`),
/// which is the DOMINANT wall — corrected by `probe_sigma_sweep`. So:
/// (1) trust this probe ONLY for the compression-vs-destruction VERDICT
/// (structure survives in h_raw = case 1), NOT for the absolute magnitude
/// (the blur cuts it; the real coarse osc is ~0.05@σ2, not +0.587). (2)
/// Always measure the GAP under a FIXED range, never a ratio across
/// affine-rescaled spaces, and never conclude "absent/smoothed" from a
/// range-normalised output — see [[feedback_measure_structure_not_compressed_output]]
/// and [[feedback_ratio_across_affine_rescaled_spaces]]. (3) Any
/// hand-rolled probe MUST replicate the FULL pipeline (incl. the blur);
/// this one does not — use `probe_sigma_sweep` / `probe_meso_cap` for the
/// real path.
#[test]
#[ignore]
fn probe_meso_distinction() {
    use ymir_core::tectonics_c1::boundary_classification::{
        classify_boundaries, oc_override_seed_mask, retarget_upper_plate_continental,
    };
    use ymir_core::tectonics_c1::distance_field::wedge_distance_intra_plate_typed;
    let iso_config = IsostasyConfig::c1_default();
    let grid = 64usize;
    let (lambda, amp, band_max) = (3.0_f64, 0.8_f64, 9.0_f64);
    let buoy = 1.0 - 2750.0_f64 / 3300.0;
    let sea_norm = 500.0_f64 / (500.0 + 4000.0);
    let sea_frac = 0.4_f64;
    let cap = 0.92_f64;
    eprintln!("#155 méso distinction — inj h_raw under BASE fixed range (buoy={buoy:.4}, sea_norm={sea_norm:.4})");
    for &seed in &[1988u64, 42u64] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0,
            iso_config: iso_config.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        let mut bi = classify_boundaries(&state.plate_id, &kin);
        retarget_upper_plate_continental(&mut bi, &state.plate_id, &state.plate_type);
        let oc_seed = oc_override_seed_mask(&bi, &state.plate_id, &state.plate_type);
        let (wd, is_oc) = wedge_distance_intra_plate_typed(&state.plate_id, &bi.upper_plate_mask, &oc_seed, 30.0);

        // BASE fixed range (replicate compute_isostasy steps 1-2 on base S̃).
        let hraw_base: Vec<f64> = (0..grid*grid).map(|k| state.s.get(k % grid, k / grid) * buoy).collect();
        let hmin = hraw_base.iter().cloned().fold(f64::INFINITY, f64::min);
        let hmax = hraw_base.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mut sorted = hraw_base.clone(); sorted.sort_by(|a,b| a.partial_cmp(b).unwrap());
        let hcap = sorted[((cap * (sorted.len()-1) as f64).round() as usize).min(sorted.len()-1)];
        let hsea = hmin + sea_frac * (hcap - hmin).max(1e-10);
        let normalize = |h: f64| -> f64 {
            if h <= hsea { ((h - hmin) / (hsea - hmin).max(1e-10)) * sea_norm }
            else { sea_norm + ((h - hsea) / (hmax - hsea).max(1e-10)) * (1.0 - sea_norm) }.clamp(0.0, 1.0)
        };

        let shell = |f: &dyn Fn(usize,usize)->f64| -> Vec<f64> {
            let mut sh = vec![(0.0,0usize); 20];
            for j in 0..grid { for i in 0..grid {
                if !is_oc.get(i,j) { continue; }
                let d = wd.get(i,j); if d < 0.5 || d > 10.0 { continue; }
                let k = ((d/0.5) as usize).min(19); sh[k].0 += f(i,j); sh[k].1 += 1;
            }}
            sh.iter().map(|c| if c.1>0 { c.0/c.1 as f64 } else { f64::NAN }).collect()
        };
        // base normalized (base range) vs inj normalized (SAME base range).
        let base_norm = shell(&|i,j| normalize(state.s.get(i,j) * buoy));
        // inject
        for j in 0..grid { for i in 0..grid {
            if is_oc.get(i,j) { let d = wd.get(i,j);
                if d >= 0.5 && d <= band_max {
                    state.s.set(i, j, (state.s.get(i,j) + amp*(std::f64::consts::TAU*d/lambda).cos()).max(0.1));
                }}
        }}
        let inj_fixed = shell(&|i,j| normalize(state.s.get(i,j) * buoy));
        let osc = |v: &[f64]| { // crest(d≈3,6) − trough(d≈4.5,7.5) mean
            let g=|x:f64| v[(x/0.5) as usize];
            ((g(3.0)+g(6.0))/2.0) - ((g(4.5)+g(7.5))/2.0)
        };
        let fmt=|v:&[f64]| v.iter().enumerate().filter(|(_,x)|x.is_finite()).map(|(k,x)| format!("d{:.1}:{:.3}",k as f64*0.5,x)).collect::<Vec<_>>().join(" ");
        eprintln!("  seed {seed} base(baseRange)    : {}", fmt(&base_norm));
        eprintln!("  seed {seed} inj (baseRange FIX): {}", fmt(&inj_fixed));
        eprintln!("  seed {seed} crest-trough osc: base={:+.4}  inj(fixed)={:+.4}  (inj>base ⇒ structure expresses under fixed range = CASE 1)", osc(&base_norm), osc(&inj_fixed));
    }
}

/// #155 méso step 1+2 — blur isolation + upscale survival.
/// Step 1: isolate the coarse wall by toggling the gaussian blur
/// (`altitude_smoothing_sigma` σ=2 vs 0) — does the injected imbricate
/// structure EXPRESS in the normalised coarse altitude once the blur is
/// off? (This probe originally also swept a `land_cap_percentile` lever;
/// that field was MEASURED BAD — over-saturates 14-46% — and REVERTED,
/// so only the σ dimension remains. The cap finding lives in
/// `stage_meso_expression.md`.)
/// Step 2 (the real gate): does the decompressed (σ=0) structure SURVIVE
/// upscale+FBM (case 3, the r~0.51 fine-S̃ wall)?
/// All in consistent normalised space (shell crest−trough gaps, never a
/// cross-space ratio — [[feedback_measure_structure_not_compressed_output]]).
#[test]
#[ignore]
fn probe_meso_cap() {
    use ymir_core::tectonics_c1::boundary_classification::{
        classify_boundaries, oc_override_seed_mask, retarget_upper_plate_continental,
    };
    use ymir_core::tectonics_c1::distance_field::wedge_distance_intra_plate_typed;
    let dir = output_dir().join("meso_cap");
    std::fs::create_dir_all(&dir).expect("create dir");
    let base_iso = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    let (lambda, amp, band_max) = (3.0_f64, 0.8_f64, 9.0_f64);
    let cfg = FbmUpscaleConfig {
        target_size: 1024, coast_warp_strength: 1.5, coast_warp_frequency: 0.5,
        coastal_amplitude_band: 0.30, amplitude_base: 0.16, submarine_damping: 0.0,
        ..Default::default()
    };
    // shell crest(d≈3,6)−trough(d≈4.5,7.5) over O-C cells, on a coarse GridF32.
    let osc_coarse = |alt: &GridF32, isoc: &dyn Fn(usize,usize)->bool, wd: &Field2D| -> f64 {
        let (mut cs, mut cn, mut ts, mut tn) = (0.0, 0usize, 0.0, 0usize);
        for j in 0..grid { for i in 0..grid {
            if !isoc(i, j) { continue; }
            let d = wd.get(i, j); let a = alt.get(i as i32, j as i32) as f64;
            if (2.5..=3.5).contains(&d) || (5.5..=6.5).contains(&d) { cs += a; cn += 1; }
            else if (4.0..=5.0).contains(&d) || (7.0..=8.0).contains(&d) { ts += a; tn += 1; }
        }}
        cs / cn.max(1) as f64 - ts / tn.max(1) as f64
    };
    eprintln!("#155 méso cap — step1 decompress (coarse osc + saturation), step2 upscale survival");
    for &seed in &[2u64, 99u64, 1337u64, 2026u64, 4138u64] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0,
            iso_config: base_iso.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        let mut bi = classify_boundaries(&state.plate_id, &kin);
        retarget_upper_plate_continental(&mut bi, &state.plate_id, &state.plate_type);
        let oc_seed = oc_override_seed_mask(&bi, &state.plate_id, &state.plate_type);
        let (wd, is_oc) = wedge_distance_intra_plate_typed(&state.plate_id, &bi.upper_plate_mask, &oc_seed, 30.0);
        let isoc = |i: usize, j: usize| is_oc.get(i, j);

        // s_base and s_inj (synthetic imbricate stack), no state mutation yet.
        let s_base = state.s.clone();
        let mut s_inj = state.s.clone();
        for j in 0..grid { for i in 0..grid {
            if is_oc.get(i, j) { let d = wd.get(i, j);
                if d >= 0.5 && d <= band_max {
                    s_inj.set(i, j, (s_base.get(i, j) + amp*(std::f64::consts::TAU*d/lambda).cos()).max(0.1));
                }}
        }}

        // Step 1: blur isolation — coarse decompression at σ=2 vs σ=0.
        for sigma in [2.0f32, 0.0] {
            let mut iso = base_iso.clone();
            iso.altitude_smoothing_sigma = sigma;
            let ab = c1_production_altitude(&s_base, &state.age, &state.plate_type, &iso, &ss);
            let ai = c1_production_altitude(&s_inj, &state.age, &state.plate_type, &iso, &ss);
            let (ob, oi) = (osc_coarse(&ab, &isoc, &wd), osc_coarse(&ai, &isoc, &wd));
            eprintln!("  seed {seed} σ={sigma}: coarse osc base={ob:+.4} inj={oi:+.4} (inj≫base@σ0 ⇒ blur was the wall)");
        }

        // Step 2: does the DECOMPRESSED structure survive upscale+FBM? Use
        // the clean coarse expression (σ=0 — blur off, the dominant wall;
        // cap=None — no saturation cost) found in step 1.
        let mut iso = base_iso.clone();
        iso.altitude_smoothing_sigma = 0.0;
        state.s = s_inj; // move injected field into state for upscale
        let coarse_alt = c1_production_altitude(&state.s, &state.age, &state.plate_type, &iso, &ss);
        let oc_coarse = osc_coarse(&coarse_alt, &isoc, &wd);
        let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
        let scale = up.heightmap.width / grid;
        // HD osc: same crest/trough d-bands, mean over each coarse cell's HD block.
        let hd_band = |lo: f64, hi1: f64, lo2: f64, hi2: f64| -> f64 {
            let (mut s, mut n) = (0.0, 0usize);
            for j in 0..grid { for i in 0..grid {
                if !is_oc.get(i, j) { continue; }
                let d = wd.get(i, j);
                if (d>=lo && d<=hi1) || (d>=lo2 && d<=hi2) {
                    for jj in 0..scale { for ii in 0..scale { s += up.heightmap.get((i*scale+ii) as i32,(j*scale+jj) as i32) as f64; n+=1; }}
                }
            }}
            s / n.max(1) as f64
        };
        let oc_hd = hd_band(2.5,3.5,5.5,6.5) - hd_band(4.0,5.0,7.0,8.0);
        eprintln!("  seed {seed} STEP2 σ=0: coarse osc={oc_coarse:+.4}  HD osc={oc_hd:+.4}  (HD≈coarse ⇒ survives upscale; HD≪coarse ⇒ smoothed = CASE 3 upscale wall)");
        save_hillshade_crop(&up.heightmap, 0, 0, up.heightmap.width, &dir.join(format!("seed{seed:05}_sigma0_inj_upscale.png")));
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 méso — altitude_smoothing_sigma sweep on the FULL production path
/// (C1 → isostasy(σ) → upscale+FBM → HD). Measures simultaneously:
/// (A) the méso UNLOCK — injected imbricate coarse osc per σ (consistent
/// normalised space across σ; the threshold where structure expresses);
/// (B) what lowering σ UNMASKS on the REAL (non-injected) production HD —
/// rendered as hillshades to inspect: abrupt tectonic steps at plate
/// boundaries, grid CURTAIN/striping (the invariance-report artefact σ may
/// hide), and coast feathering (the FALSE suspect — should stay absent,
/// it's the upscale coastal_band's job not σ's). σ's stated purpose:
/// "smooths sharp tectonic transitions" — NOT anti-feathering (#151).
#[test]
#[ignore]
fn probe_sigma_sweep() {
    use ymir_core::tectonics_c1::boundary_classification::{
        classify_boundaries, oc_override_seed_mask, retarget_upper_plate_continental,
    };
    use ymir_core::tectonics_c1::distance_field::wedge_distance_intra_plate_typed;
    let dir = output_dir().join("sigma_sweep");
    std::fs::create_dir_all(&dir).expect("create dir");
    let base_iso = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    let (lambda, amp, band_max) = (3.0_f64, 0.8_f64, 9.0_f64);
    let sigmas = [0.0f32, 0.5, 1.0, 1.5, 2.0];
    let cfg = FbmUpscaleConfig {
        target_size: 1024, coast_warp_strength: 1.5, coast_warp_frequency: 0.5,
        coastal_amplitude_band: 0.30, amplitude_base: 0.16, submarine_damping: 0.0,
        ..Default::default()
    };
    eprintln!("#155 σ sweep — (A) méso unlock coarse osc / (B) real-production HD hillshades for artefacts");
    for &seed in &[1988u64, 42u64, 1337u64, 4138u64] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0,
            iso_config: base_iso.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        let mut bi = classify_boundaries(&state.plate_id, &kin);
        retarget_upper_plate_continental(&mut bi, &state.plate_id, &state.plate_type);
        let oc_seed = oc_override_seed_mask(&bi, &state.plate_id, &state.plate_type);
        let (wd, is_oc) = wedge_distance_intra_plate_typed(&state.plate_id, &bi.upper_plate_mask, &oc_seed, 30.0);

        let s_base = state.s.clone();
        let mut s_inj = state.s.clone();
        for j in 0..grid { for i in 0..grid {
            if is_oc.get(i, j) { let d = wd.get(i, j);
                if d >= 0.5 && d <= band_max {
                    s_inj.set(i, j, (s_base.get(i, j) + amp*(std::f64::consts::TAU*d/lambda).cos()).max(0.1));
                }}
        }}
        let osc_coarse = |alt: &GridF32| -> f64 {
            let (mut cs, mut cn, mut ts, mut tn) = (0.0, 0usize, 0.0, 0usize);
            for j in 0..grid { for i in 0..grid {
                if !is_oc.get(i, j) { continue; }
                let d = wd.get(i, j); let a = alt.get(i as i32, j as i32) as f64;
                if (2.5..=3.5).contains(&d) || (5.5..=6.5).contains(&d) { cs += a; cn += 1; }
                else if (4.0..=5.0).contains(&d) || (7.0..=8.0).contains(&d) { ts += a; tn += 1; }
            }}
            cs / cn.max(1) as f64 - ts / tn.max(1) as f64
        };

        for &sigma in &sigmas {
            let mut iso = base_iso.clone(); iso.altitude_smoothing_sigma = sigma;
            // (A) méso unlock — injected coarse osc.
            let ai = c1_production_altitude(&s_inj, &state.age, &state.plate_type, &iso, &ss);
            eprintln!("  seed {seed} σ={sigma}: méso inj coarse osc = {:+.4}", osc_coarse(&ai));
            // (B) REAL production HD hillshade (artefacts: steps / striping / feathering).
            state.s = s_base.clone();
            let up_real = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
            save_hillshade_crop(&up_real.heightmap, 0, 0, up_real.heightmap.width,
                &dir.join(format!("seed{seed:05}_s{sigma:.1}_real.png")));
            // (A-visual) injected HD hillshade for the franc-O-C seeds.
            if seed == 1988 || seed == 1337 {
                state.s = s_inj.clone();
                let up_inj = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
                save_hillshade_crop(&up_inj.heightmap, 0, 0, up_inj.heightmap.width,
                    &dir.join(format!("seed{seed:05}_s{sigma:.1}_inj.png")));
            }
        }
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 méso σ maillon — effect of reduced altitude_smoothing_sigma on
/// the REAL macro O-C ridge (1a+1b-i, NOT synthetic injection — the first
/// test of reduced σ on real structure). O-C-band vs passive-band mean
/// altitude GAP in consistent normalised space (same normalisation across
/// σ; gap not ratio — [[feedback_measure_structure_not_compressed_output]]).
/// Expected: the real macro ridge gets SHARPER (bigger gap) at low σ
/// because σ no longer smooths it. (Real hillshades per σ: sigma_sweep/.)
#[test]
#[ignore]
fn probe_sigma_macro() {
    use ymir_core::tectonics_c1::boundary_classification::{
        classify_boundaries, oc_override_seed_mask, retarget_upper_plate_continental,
    };
    use ymir_core::tectonics_c1::distance_field::wedge_distance_intra_plate_typed;
    let base_iso = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    eprintln!("#155 σ macro — REAL O-C ridge gap (O-C band − passive band, normalised) vs σ");
    for &seed in &[1988u64, 42u64, 1337u64, 4138u64] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0,
            iso_config: base_iso.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        let mut bi = classify_boundaries(&state.plate_id, &kin);
        retarget_upper_plate_continental(&mut bi, &state.plate_id, &state.plate_type);
        let oc_seed = oc_override_seed_mask(&bi, &state.plate_id, &state.plate_type);
        let (wd, is_oc) = wedge_distance_intra_plate_typed(&state.plate_id, &bi.upper_plate_mask, &oc_seed, 30.0);
        let is_cont = |i: usize, j: usize| matches!(state.plate_type.get(i, j), PlateType::Continental);
        // dist-to-ocean for the passive band.
        let mut d_ocean = vec![usize::MAX; grid*grid];
        let mut q = std::collections::VecDeque::new();
        for j in 0..grid { for i in 0..grid { if !is_cont(i,j) { d_ocean[j*grid+i]=0; q.push_back((i,j)); } }}
        while let Some((i,j)) = q.pop_front() { let dc=d_ocean[j*grid+i];
            for (di,dj) in [(-1i32,0i32),(1,0),(0,-1),(0,1)] { let (ni,nj)=(i as i32+di,j as i32+dj);
                if ni<0||nj<0||ni as usize>=grid||nj as usize>=grid {continue;}
                let (ni,nj)=(ni as usize,nj as usize);
                if d_ocean[nj*grid+ni]>dc+1 { d_ocean[nj*grid+ni]=dc+1; q.push_back((ni,nj)); }}}
        for &sigma in &[2.0f32, 0.5, 0.0] {
            let mut iso = base_iso.clone(); iso.altitude_smoothing_sigma = sigma;
            let alt = c1_production_altitude(&state.s, &state.age, &state.plate_type, &iso, &ss);
            let (mut ocs, mut ocn, mut pas, mut pan) = (0.0f64, 0usize, 0.0f64, 0usize);
            for j in 0..grid { for i in 0..grid {
                if !is_cont(i,j) { continue; }
                let a = alt.get(i as i32, j as i32) as f64; let dw = wd.get(i,j);
                if is_oc.get(i,j) && dw>=1.0 && dw<=6.0 { ocs += a; ocn += 1; }
                else if !is_oc.get(i,j) { let d=d_ocean[j*grid+i]; if d>=1 && d<=6 { pas += a; pan += 1; } }
            }}
            let (ocm, pam) = (ocs/ocn.max(1) as f64, pas/pan.max(1) as f64);
            eprintln!("  seed {seed} σ={sigma}: O-C band={ocm:.4} passive={pam:.4} GAP={:+.4}", ocm - pam);
        }
    }
}

/// #155 méso B viability — the central bet: does ISOTROPIC hydraulic
/// erosion on an already-tectonically-ORIENTED FBM texture yield distinct
/// PARALLEL cordilleras, or generic DENDRITIC dissection that erases the
/// orientation? Throwaway (no wiring, no prod change): on the REAL macro
/// ridge (σ=0.5, no synthetic injection), upscale with a BOOSTED
/// anisotropic FBM (amplitude/anisotropy/slope-factor up), then
/// run_erosion. Compare pre- vs post-erosion hillshades:
///   - cordilleras survive/sharpen → B viable, wire HD erosion cleanly;
///   - dendritic dissection (orientation erased) → B needs an ORIENTED
///     erosion (anisotropy input in hydraulic.rs), more work to scope.
/// Visual is authoritative ([[feedback_visual_over_scalar]]).
#[test]
#[ignore]
fn probe_meso_erosion_orientation() {
    use ymir_core::erosion::hydraulic::{run_erosion, ErosionConfig};
    let dir = output_dir().join("meso_erosion_orient");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso = IsostasyConfig::c1_default(); // σ=0.5
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    // BOOSTED anisotropic FBM: stronger oriented ridges on the macro slopes.
    let cfg = FbmUpscaleConfig {
        target_size: 1024,
        amplitude_base: 0.30,        // 0.16 → 0.30 (taller méso ridges)
        max_anisotropy: 8.0,         // 3.0 → 8.0 (long margin-parallel stretch)
        amplitude_slope_factor: 5.0, // 3.0 → 5.0 (amplitude on the ridge flanks)
        coast_warp_strength: 1.5, coast_warp_frequency: 0.5,
        coastal_amplitude_band: 0.30, submarine_damping: 0.0,
        ..Default::default()
    };
    let ero_cfg = ErosionConfig { num_droplets: 2_000_000, batch_size: 50_000, ..Default::default() };
    eprintln!("#155 méso B viability — isotropic erosion on oriented FBM: cordilleras or dendritic?");
    for &seed in &[1988u64, 42u64] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0,
            iso_config: iso.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
        save_hillshade_crop(&up.heightmap, 0, 0, up.heightmap.width,
            &dir.join(format!("seed{seed:05}_1_oriented_fbm.png")));
        let eroded = run_erosion(&up.heightmap, &ero_cfg, &WorldSeed::new(seed), |_, _, _| true);
        save_hillshade_crop(&eroded.heightmap, 0, 0, eroded.heightmap.width,
            &dir.join(format!("seed{seed:05}_2_eroded.png")));
        eprintln!("  seed {seed}: oriented-FBM + eroded hillshades written");
    }
    eprintln!("  out = {}", dir.display());
}

/// Land/sea binary map (Living Landz output): land (h > sea) green, sea blue.
fn save_binarymap(h: &GridF32, sea: f32, path: &Path) {
    let (nx, ny) = (h.width, h.height);
    let mut img = ImageBuffer::<Rgb<u8>, Vec<u8>>::new(nx as u32, ny as u32);
    for j in 0..ny { for i in 0..nx {
        let v = h.get(i as i32, j as i32);
        let c = if v > sea { Rgb([70, 140, 70]) } else { Rgb([40, 80, 160]) };
        img.put_pixel(i as u32, (ny - 1 - j) as u32, c);
    }}
    img.save(path).expect("save binarymap PNG");
}

/// #155 PRODUCT état-des-lieux — the C1 product at the Living Landz target
/// 2048², ALL seeds, WITH sea level visible (hypsometric) + binarymap.
/// Pipeline = C1 (σ=0.5) → upscale_from_c1 (prod #151 FBM) → run_erosion
/// (the decided méso = dendritic dissection; NOTE erosion is NOT YET WIRED
/// into upscale_from_c1, Case 2a — applied manually here pending that
/// maillon). Full world, no crop. Judge: credible mountains on all seeds?
/// coasts natural / relief stays on land? oceans flat (known gap)?
#[test]
#[ignore]
fn product_etat_des_lieux() {
    use ymir_core::erosion::hydraulic::{run_erosion, ErosionConfig};
    let dir = output_dir().join("product_2048");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso = IsostasyConfig::c1_default(); // σ=0.5
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    let cfg = FbmUpscaleConfig {
        target_size: 2048, coast_warp_strength: 1.5, coast_warp_frequency: 0.5,
        coastal_amplitude_band: 0.30, amplitude_base: 0.16, submarine_damping: 0.0,
        ..Default::default()
    };
    let ero_cfg = ErosionConfig { num_droplets: 4_000_000, batch_size: 100_000, ..Default::default() };
    let sea = 0.5f32; // upscale_from_c1 maps sea level to 0.5
    eprintln!("#155 product état-des-lieux — 2048², all seeds, hypsometric + binarymap (C1 σ0.5 + upscale + erosion)");
    for &seed in &[2u64, 42, 99, 1337, 1988, 2026, 4138] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0,
            iso_config: iso.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
        let eroded = run_erosion(&up.heightmap, &ero_cfg, &WorldSeed::new(seed), |_, _, _| true);
        let h = &eroded.heightmap;
        save_heightmap01(h, &dir.join(format!("seed{seed:05}_hypso.png")));
        save_binarymap(h, sea, &dir.join(format!("seed{seed:05}_binary.png")));
        save_hillshade_crop(h, 0, 0, h.width, &dir.join(format!("seed{seed:05}_hillshade.png")));
        let land = h.data.iter().filter(|&&v| v > sea).count();
        eprintln!("  seed {seed}: {}² written, land = {:.1}%", h.width, 100.0 * land as f64 / (h.width*h.height) as f64);
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 HD-production wiring ACCEPTANCE — upscale_from_c1 with the canonical
/// c1_hd_production config now delivers the eroded (dissected) product the
/// 2048² état-des-lieux judged. Confirms: (1) product ≡ état-des-lieux
/// (hypso/hillshade), (2) determinism run×2 byte-identical, (3) cost at
/// 2048², (4) sediment hook present + slope recomputed post-erosion.
#[test]
#[ignore]
fn hd_production_acceptance() {
    let dir = output_dir().join("hd_production");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    let cfg = FbmUpscaleConfig::c1_hd_production(2048);
    eprintln!("#155 HD production acceptance — c1_hd_production(2048), seeds 1988/42/2026");
    eprintln!("  cfg: target={} erosion={}", cfg.target_size,
        cfg.erosion.as_ref().map(|e| e.num_droplets).unwrap_or(0));
    for &seed in &[1988u64, 42, 2026] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0,
            iso_config: iso.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});

        // COST + product: ONE 2048² pass (the heavy, judged config).
        let t0 = std::time::Instant::now();
        let up1 = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
        let dt = t0.elapsed();
        let sed = up1.sediment.is_some();
        let land = up1.heightmap.data.iter().filter(|&&v| v > 0.5).count();
        // DETERMINISM: proven cheaply at 512² (run×2), NOT at the heavy 2048².
        let cfg_det = FbmUpscaleConfig::c1_hd_production(512);
        let d1 = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg_det);
        let d2 = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg_det);
        let det = d1.heightmap.data == d2.heightmap.data;
        eprintln!("  seed {seed}: {}² in {dt:.2?}  det(512²)={}  sediment={}  slope_ok={}  land={:.1}%",
            up1.heightmap.width, if det {"OK"} else {"MISMATCH"}, sed,
            up1.slope.data.iter().any(|&s| s > 0.0),
            100.0 * land as f64 / (up1.heightmap.width*up1.heightmap.height) as f64);
        save_heightmap01(&up1.heightmap, &dir.join(format!("seed{seed:05}_hypso.png")));
        save_hillshade_crop(&up1.heightmap, 0, 0, up1.heightmap.width, &dir.join(format!("seed{seed:05}_hillshade.png")));
        assert!(det, "seed {seed}: HD production not deterministic run×2 (512²)");
        assert!(sed, "seed {seed}: sediment hook missing");
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 INTERIOR relief — economical-hypothesis test. The interior is flat
/// because NOTHING couples craton/age → continental altitude (EH treats
/// cratons uniformly; isostasy = S̃ only; age → oceanic only). Hypothesis:
/// the interior needs no new closure — just LARGE-SCALE altitude variation
/// (from the cratonic_mask we already HAVE) so the aval (σ=0.5 + FBM +
/// erosion) dresses it into relief, the way the macro gradient lets the
/// aval work at margins. Test (throwaway, post-loop = advection-free):
/// bump continental S̃ by Δ·smooth(cratonic_mask) → upscale_from_c1 (erosion
/// ON) → does interior relief APPEAR over cratons vs the flat baseline?
/// Visual authoritative; interior altitude std (continental, away from
/// coast) base vs bumped as a scalar companion.
#[test]
#[ignore]
fn probe_interior_craton() {
    let dir = output_dir().join("interior_craton");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    let cfg = FbmUpscaleConfig::c1_hd_production(1024); // 1024² for speed
    let delta = 1.0f64; // S̃ bump in craton cores (→ ~+0.17 raw altitude)
    eprintln!("#155 interior craton — does the aval dress a gentle craton gradient into relief? (Δ={delta})");
    for &seed in &[42u64, 2026, 1988] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0,
            iso_config: iso.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        let n_craton = (0..grid*grid).filter(|k| state.cratonic_mask.get(k % grid, k / grid)).count();

        // BASELINE (flat interior).
        let up_base = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
        save_hillshade_crop(&up_base.heightmap, 0, 0, up_base.heightmap.width, &dir.join(format!("seed{seed:05}_1_base.png")));

        // Smooth craton field → gentle large-scale gradient.
        let mut cf = GridF32::new(grid, grid, 0.0);
        for j in 0..grid { for i in 0..grid {
            if state.cratonic_mask.get(i, j) { cf.set(i, j, 1.0); }
        }}
        let cf = cf.gaussian_blur(4.0);
        // BUMP continental S̃ post-loop (advection-free, large-scale).
        for j in 0..grid { for i in 0..grid {
            if matches!(state.plate_type.get(i, j), PlateType::Continental) {
                let add = delta * cf.get(i as i32, j as i32) as f64;
                state.s.set(i, j, state.s.get(i, j) + add);
            }
        }}
        let up_bump = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
        save_hillshade_crop(&up_bump.heightmap, 0, 0, up_bump.heightmap.width, &dir.join(format!("seed{seed:05}_2_craton_bumped.png")));
        save_heightmap01(&up_bump.heightmap, &dir.join(format!("seed{seed:05}_2_craton_hypso.png")));

        // Interior relief proxy: std of HD altitude over land cells (>0.5).
        let interior_std = |h: &GridF32| -> f64 {
            let land: Vec<f64> = h.data.iter().filter(|&&v| v > 0.5).map(|&v| v as f64).collect();
            if land.is_empty() { return 0.0; }
            let m = land.iter().sum::<f64>() / land.len() as f64;
            (land.iter().map(|v| (v - m).powi(2)).sum::<f64>() / land.len() as f64).sqrt()
        };
        eprintln!("  seed {seed}: cratons={} cells  land-alt std base={:.4} bumped={:.4}",
            n_craton, interior_std(&up_base.heightmap), interior_std(&up_bump.heightmap));
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 A′ CALIBRATION — judge the now-active (default) craton elevation on
/// RENDERED altitude in METERS (not the raw S̃ ratio). Defaults are A′-on
/// (craton_thickness 1.25 + craton_resist 0.2). Measures craton vs
/// non-craton LAND mean altitude (HD, c1_hd_production) in metres
/// (norm→m: (v−0.5)·2·max_elevation_m, max_elevation_m=4000). Target ~300-
/// 500 m worn shield. Too high → check note (b) (base-erosion melting the
/// non-craton denominator) before any in-band resist nudge.
#[test]
#[ignore]
fn probe_craton_calibration() {
    let dir = output_dir().join("craton_calib");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    let cfg = FbmUpscaleConfig::c1_hd_production(1024);
    let max_elev_m = 4000.0f64;
    eprintln!("#155 A′ calibration — craton vs non-craton land altitude (m), defaults A′-on");
    for &seed in &[42u64, 1988, 2026] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default(); // A′ active
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0,
            iso_config: iso.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
        let scale = up.heightmap.width / grid;
        let (mut cs, mut cn, mut ns, mut nn) = (0.0f64, 0usize, 0.0f64, 0usize);
        for j in 0..grid { for i in 0..grid {
            let craton = state.cratonic_mask.get(i, j);
            for jj in 0..scale { for ii in 0..scale {
                let v = up.heightmap.get((i*scale+ii) as i32, (j*scale+jj) as i32) as f64;
                if v <= 0.5 { continue; } // land only
                let m = (v - 0.5) * 2.0 * max_elev_m;
                if craton { cs += m; cn += 1; } else { ns += m; nn += 1; }
            }}
        }}
        let (cm, nm) = (cs / cn.max(1) as f64, ns / nn.max(1) as f64);
        eprintln!("  seed {seed}: craton land mean = {cm:.0} m  non-craton = {nm:.0} m  craton above = {:.0} m  (target shield ~300-500)", cm - nm);
        save_heightmap01(&up.heightmap, &dir.join(format!("seed{seed:05}_hypso.png")));
        save_hillshade_crop(&up.heightmap, 0, 0, up.heightmap.width, &dir.join(format!("seed{seed:05}_hillshade.png")));
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 POST-A′ état-des-lieux — the product WITH interior cratons, via the
/// canonical c1_hd_production (A′ active by default: thick + erosion-resistant
/// cratons). 7 seeds, 2048², hypsometric (sea level) + hillshade + binarymap.
/// Compares to the pre-A′ snapshot (product_2048): interiors now carry craton
/// shields, margins keep their O-C ridges, oceans still flat (known gap).
#[test]
#[ignore]
fn product_etat_des_lieux_aprime() {
    let dir = output_dir().join("product_2048_aprime");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    let cfg = FbmUpscaleConfig::c1_hd_production(2048); // erosion wired internally
    let sea = 0.5f32;
    eprintln!("#155 post-A′ état-des-lieux — 2048², all seeds, via c1_hd_production (A′ default-on)");
    for &seed in &[2u64, 42, 99, 1337, 1988, 2026, 4138] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default(); // A′ active
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0,
            iso_config: iso.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
        let h = &up.heightmap;
        save_heightmap01(h, &dir.join(format!("seed{seed:05}_hypso.png")));
        save_binarymap(h, sea, &dir.join(format!("seed{seed:05}_binary.png")));
        save_hillshade_crop(h, 0, 0, h.width, &dir.join(format!("seed{seed:05}_hillshade.png")));
        let land = h.data.iter().filter(|&&v| v > sea).count();
        eprintln!("  seed {seed}: {}² written, land = {:.1}%", h.width, 100.0 * land as f64 / (h.width*h.height) as f64);
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 PROMINENCE attribution — why is the orogen/craton differential
/// crushed? Decompose end-of-loop S̃ (raw, consistent space) into three
/// continental zones — margin orogen (O-C wedge), craton, platform — and
/// test the hypothesis (i) orogens cap at EH h_eq=2.0 vs (ii) cratons
/// drift up. Discriminator = a counterfactual run with NO craton
/// differential (init ratio 1.0 + resist 1.0): does the margin/craton
/// ratio re-open? + how many margin cells are above h_eq (EH-collapsed).
#[test]
#[ignore]
fn probe_prominence_attribution() {
    use ymir_core::tectonics_c1::boundary_classification::{
        classify_boundaries, oc_override_seed_mask, retarget_upper_plate_continental,
    };
    use ymir_core::tectonics_c1::distance_field::wedge_distance_intra_plate_typed;
    use ymir_core::tectonics_c1::closures::erosion::params::ErosionParams;
    let iso = IsostasyConfig::c1_default();
    let grid = 64usize;
    let h_eq = 2.0f64;
    let p95 = |mut v: Vec<f64>| -> f64 { if v.is_empty() { return f64::NAN; } v.sort_by(|a,b| a.partial_cmp(b).unwrap()); v[(v.len()*95/100).min(v.len()-1)] };
    let med = |mut v: Vec<f64>| -> f64 { if v.is_empty() { return f64::NAN; } v.sort_by(|a,b| a.partial_cmp(b).unwrap()); v[v.len()/2] };
    eprintln!("#155 prominence attribution — end-of-loop S̃ by zone (margin orogen / craton / platform)");
    for &seed in &[42u64, 2026, 1988] {
        // Geometry (zones) — from a default init (plate_id identical across variants).
        let geom = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let kin0 = PlateKinematics::preset_phase_1_1(geom.num_plates);
        let mut bi = classify_boundaries(&geom.plate_id, &kin0);
        retarget_upper_plate_continental(&mut bi, &geom.plate_id, &geom.plate_type);
        let oc_seed = oc_override_seed_mask(&bi, &geom.plate_id, &geom.plate_type);
        let (wd, is_oc) = wedge_distance_intra_plate_typed(&geom.plate_id, &bi.upper_plate_mask, &oc_seed, 30.0);
        let is_cont = |i: usize, j: usize| matches!(geom.plate_type.get(i, j), PlateType::Continental);
        // EXCLUSIVE zones (de-confound craton∩wedge overlap): 1 pure orogen
        // (wedge ∧ ¬craton), 2 pure craton (craton ∧ ¬wedge), 3 platform,
        // 4 craton∩wedge (reported apart). wedge = O-C && wd<max_d.
        let mut zone = vec![0u8; grid*grid];
        for j in 0..grid { for i in 0..grid {
            if !is_cont(i, j) { continue; }
            let wedge = is_oc.get(i, j) && wd.get(i, j) < 30.0;
            let crat = geom.cratonic_mask.get(i, j);
            zone[j*grid+i] = match (wedge, crat) {
                (true, false) => 1,
                (false, true) => 2,
                (false, false) => 3,
                (true, true) => 4,
            };
        }}
        let zsamp = |st: &C1State, z: u8| -> Vec<f64> {
            let mut v = Vec::new();
            for j in 0..grid { for i in 0..grid { if zone[j*grid+i]==z { v.push(st.s.get(i,j)); } }} v
        };

        let run = |ratio: f64, resist: f64| -> C1State {
            let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams { craton_thickness_ratio: ratio, ..Default::default() });
            let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
            let closures = C1Closures { erosion: ErosionParams { craton_resist: resist, ..ErosionParams::default() }, ..C1Closures::default() };
            let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
            run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
            state
        };

        // craton S̃ at t=0 (A′ init, =1.25×base).
        let c_t0 = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let craton_t0 = p95(zsamp(&c_t0, 2));

        // A′ (default: ratio 1.25, resist 0.2).
        let a = run(1.25, 0.2);
        let (mg, cr, pf) = (p95(zsamp(&a,1)), p95(zsamp(&a,2)), p95(zsamp(&a,3)));
        let mg_med = med(zsamp(&a,1)); let cr_med = med(zsamp(&a,2));
        let above_heq = zsamp(&a,1).iter().filter(|&&v| v > h_eq).count();
        let n_mg = zsamp(&a,1).len();
        eprintln!("  seed {seed} [A′]: margin p95={mg:.3} (med {mg_med:.3}) craton p95={cr:.3} (med {cr_med:.3}, t0={craton_t0:.3}) platform p95={pf:.3}  ratio mg/cr={:.2}  margin>h_eq: {}/{} ({:.0}%)",
            mg/cr.max(1e-6), above_heq, n_mg, 100.0*above_heq as f64/n_mg.max(1) as f64);

        // Counterfactual: NO craton differential (ratio 1.0, resist 1.0).
        let cf = run(1.0, 1.0);
        let (mgc, crc) = (p95(zsamp(&cf,1)), p95(zsamp(&cf,2)));
        eprintln!("  seed {seed} [CF no-craton-diff]: margin p95={mgc:.3} craton p95={crc:.3}  ratio mg/cr={:.2}  (re-opens if > A′)", mgc/crc.max(1e-6));
    }
}

/// #155 OROGEN EQUILIBRIUM — DS targets h_max=2.5 but orogens reach ~1.0.
/// Which agent brakes? Counterfactual chain-of-innocence: margin S̃ p95
/// end-of-loop in {baseline A′, erosion-C1 OFF}. Erosion-OFF rises toward
/// 2.5 → erosion rabote (mech 2). Stays low → DS rate too slow (mech 1)
/// (advection mech 3 a-priori out: continental orogen is rigid/sealed,
/// no-flux bidirectional). Margin = O-C wedge (is_oc && wd<max_d). 42/1988
/// cratons≡margins (orogen on cratonic crust, S̃ starts 1.25×); 2026 =
/// clean non-craton margin. Raw S̃, consistent space.
#[test]
#[ignore]
fn probe_orogen_equilibrium() {
    use ymir_core::tectonics_c1::boundary_classification::{
        classify_boundaries, oc_override_seed_mask, retarget_upper_plate_continental,
    };
    use ymir_core::tectonics_c1::distance_field::wedge_distance_intra_plate_typed;
    use ymir_core::tectonics_c1::closures::erosion::params::ErosionParams;
    let iso = IsostasyConfig::c1_default();
    let grid = 64usize;
    let p95 = |mut v: Vec<f64>| -> f64 { if v.is_empty() { return f64::NAN; } v.sort_by(|a,b| a.partial_cmp(b).unwrap()); v[(v.len()*95/100).min(v.len()-1)] };
    eprintln!("#155 orogen equilibrium — margin S̃ p95: baseline vs erosion-OFF (DS h_max=2.5)");
    for &seed in &[42u64, 1988, 2026] {
        let geom = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let kin0 = PlateKinematics::preset_phase_1_1(geom.num_plates);
        let mut bi = classify_boundaries(&geom.plate_id, &kin0);
        retarget_upper_plate_continental(&mut bi, &geom.plate_id, &geom.plate_type);
        let oc_seed = oc_override_seed_mask(&bi, &geom.plate_id, &geom.plate_type);
        let (wd, is_oc) = wedge_distance_intra_plate_typed(&geom.plate_id, &bi.upper_plate_mask, &oc_seed, 30.0);
        let mut margin = vec![false; grid*grid];
        for j in 0..grid { for i in 0..grid {
            if matches!(geom.plate_type.get(i,j), PlateType::Continental) && is_oc.get(i,j) && wd.get(i,j) < 30.0 { margin[j*grid+i] = true; }
        }}
        let msamp = |st: &C1State| -> Vec<f64> { let mut v=Vec::new(); for j in 0..grid { for i in 0..grid { if margin[j*grid+i] { v.push(st.s.get(i,j)); } }} v };
        let run = |erosion_on: bool| -> C1State {
            let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
            let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
            let mut closures = C1Closures::default();
            if !erosion_on { closures.erosion = ErosionParams { enabled: false, ..ErosionParams::default() }; }
            let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
            run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
            state
        };
        let base = p95(msamp(&run(true)));
        let no_ero = p95(msamp(&run(false)));
        eprintln!("  seed {seed}: margin p95 baseline={base:.3}  erosion-OFF={no_ero:.3}  (rises→2.5 = erosion brakes; stays low = DS rate)  [DS h_max=2.5]");
    }
}

/// #155 WORN-SHIELD sweep — re-calibrate craton_resist DOWN (cratons lower)
/// to reopen prominence + plausible shield height, without re-inverting.
/// Orogens are EH-ceiling-bound ~2.0 → prominence = cratons DOWN. A′
/// resist=0.2 (5×) over-elevates (cratons ~1.3, prominence ~1.6×, height
/// ~600-1100m). Sweep resist {0.2,0.33,0.5,0.7} (5×,3×,2×,1.4×) on RENDERED
/// altitude (the lesson — not raw S̃). ratio_craton=1.25 fixed. 2026 =
/// decisive (separate craton + inversion floor); 42/1988 control
/// (craton≡margin). If no in-band (3-10×, resist 0.1-0.33) value hits
/// worn-height + prominence + no-inversion → STRUCTURAL signal (worn-init
/// needed, not a resist out-of-band).
#[test]
#[ignore]
fn probe_worn_shield_sweep() {
    use ymir_core::tectonics_c1::boundary_classification::{
        classify_boundaries, oc_override_seed_mask, retarget_upper_plate_continental,
    };
    use ymir_core::tectonics_c1::distance_field::wedge_distance_intra_plate_typed;
    use ymir_core::tectonics_c1::closures::erosion::params::ErosionParams;
    let dir = output_dir().join("worn_shield");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    let cfg = FbmUpscaleConfig::c1_hd_production(1024);
    let p95 = |mut v: Vec<f64>| -> f64 { if v.is_empty() { return f64::NAN; } v.sort_by(|a,b| a.partial_cmp(b).unwrap()); v[(v.len()*95/100).min(v.len()-1)] };
    let mean = |v: &[f64]| -> f64 { if v.is_empty() { return f64::NAN; } v.iter().sum::<f64>()/v.len() as f64 };
    eprintln!("#155 worn-shield sweep — RENDERED altitude by zone vs craton_resist (ratio_craton=1.25 fixed)");
    for &seed in &[2026u64, 42, 1988] {
        let geom = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let kin0 = PlateKinematics::preset_phase_1_1(geom.num_plates);
        let mut bi = classify_boundaries(&geom.plate_id, &kin0);
        retarget_upper_plate_continental(&mut bi, &geom.plate_id, &geom.plate_type);
        let oc_seed = oc_override_seed_mask(&bi, &geom.plate_id, &geom.plate_type);
        let (wd, is_oc) = wedge_distance_intra_plate_typed(&geom.plate_id, &bi.upper_plate_mask, &oc_seed, 30.0);
        // zone: 1 orogen(wedge¬craton), 2 craton(craton¬wedge), 3 platform, 4 both
        let mut zone = vec![0u8; grid*grid];
        for j in 0..grid { for i in 0..grid {
            if !matches!(geom.plate_type.get(i,j), PlateType::Continental) { continue; }
            let w = is_oc.get(i,j) && wd.get(i,j) < 30.0; let c = geom.cratonic_mask.get(i,j);
            zone[j*grid+i] = match (w,c) { (true,false)=>1, (false,true)=>2, (false,false)=>3, (true,true)=>4 };
        }}
        eprintln!("  --- seed {seed} ---");
        for &resist in &[0.2f64, 0.33, 0.5, 0.7] {
            let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
            let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
            let closures = C1Closures { erosion: ErosionParams { craton_resist: resist, ..ErosionParams::default() }, ..C1Closures::default() };
            let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
            run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
            let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
            let scale = up.heightmap.width / grid;
            let zalt = |z: u8| -> Vec<f64> { let mut v=Vec::new(); for j in 0..grid { for i in 0..grid { if zone[j*grid+i]==z { for jj in 0..scale { for ii in 0..scale { let a=up.heightmap.get((i*scale+ii) as i32,(j*scale+jj) as i32) as f64; if a>0.5 { v.push(a); } }} } }} v };
            let (oro, crat, plat) = (zalt(1), zalt(2), zalt(3));
            let to_m = |norm: f64| (norm - 0.5) * 8000.0; // conversion-dependent, flagged
            let crat_p95 = p95(crat.clone());
            let prom = p95(oro.clone()) / crat_p95.max(1e-6);
            let inverted = mean(&crat) < mean(&plat); // craton below platform = inverted
            eprintln!("    resist={resist} ({:.0}×): craton p95={crat_p95:.3} (~{:.0} m above sea) prominence oro/cra={prom:.2}  inverted(vs platform)={inverted}",
                1.0/resist, to_m(crat_p95));
        }
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 B-JORDAN sweep — spatial cratonic density lowers thick-craton
/// altitude (compositional isostasy) WITHOUT touching S̃ (orthogonal to
/// resist). Sweep iso.craton_rho_crust {None, 2850, 2900, 2950} on RENDERED
/// altitude. 2026 = decisive (separate craton): craton height plausible?
/// prominence oro/craton → ~3×? altitude-order craton vs platform (Jordan
/// too strong → altitude inversion even though S̃ ordered — resist covers
/// S̃ NOT altitude). 42/1988 control: cratons≡margins → does Jordan LOWER
/// the orogens there (craton∩wedge zone)? Loop S̃ shared (Jordan is
/// altitude-only) → one loop/seed, upscale per density.
#[test]
#[ignore]
fn probe_jordan_sweep() {
    use ymir_core::tectonics_c1::boundary_classification::{
        classify_boundaries, oc_override_seed_mask, retarget_upper_plate_continental,
    };
    use ymir_core::tectonics_c1::distance_field::wedge_distance_intra_plate_typed;
    let iso0 = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    let cfg = FbmUpscaleConfig::c1_hd_production(1024);
    let p95 = |mut v: Vec<f64>| -> f64 { if v.is_empty() { return f64::NAN; } v.sort_by(|a,b| a.partial_cmp(b).unwrap()); v[(v.len()*95/100).min(v.len()-1)] };
    let mean = |v: &[f64]| -> f64 { if v.is_empty() { return f64::NAN; } v.iter().sum::<f64>()/v.len() as f64 };
    let to_m = |n: f64| (n - 0.5) * 8000.0; // conversion-dependent, flagged
    eprintln!("#155 B-Jordan sweep — RENDERED altitude by zone vs craton_rho_crust");
    for &seed in &[2026u64, 42, 1988] {
        let geom = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let kin0 = PlateKinematics::preset_phase_1_1(geom.num_plates);
        let mut bi = classify_boundaries(&geom.plate_id, &kin0);
        retarget_upper_plate_continental(&mut bi, &geom.plate_id, &geom.plate_type);
        let oc_seed = oc_override_seed_mask(&bi, &geom.plate_id, &geom.plate_type);
        let (wd, is_oc) = wedge_distance_intra_plate_typed(&geom.plate_id, &bi.upper_plate_mask, &oc_seed, 30.0);
        let mut zone = vec![0u8; grid*grid];
        for j in 0..grid { for i in 0..grid {
            if !matches!(geom.plate_type.get(i,j), PlateType::Continental) { continue; }
            let w = is_oc.get(i,j) && wd.get(i,j) < 30.0; let c = geom.cratonic_mask.get(i,j);
            zone[j*grid+i] = match (w,c) { (true,false)=>1, (false,true)=>2, (false,false)=>3, (true,true)=>4 };
        }}
        // C1 loop once (defaults A′: resist 0.2, thick 1.25). Jordan is altitude-only.
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso0.clone(), drainage_max_distance: 30 };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        eprintln!("  --- seed {seed} ---");
        for rho in [None, Some(2850.0f32), Some(2900.0), Some(2950.0)] {
            let mut iso = iso0.clone(); iso.craton_rho_crust = rho;
            let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
            let scale = up.heightmap.width / grid;
            let zalt = |z: u8| -> Vec<f64> { let mut v=Vec::new(); for j in 0..grid { for i in 0..grid { if zone[j*grid+i]==z { for jj in 0..scale { for ii in 0..scale { let a=up.heightmap.get((i*scale+ii) as i32,(j*scale+jj) as i32) as f64; if a>0.5 { v.push(a); } }} } }} v };
            let (oro, crat, plat, cw) = (zalt(1), zalt(2), zalt(3), zalt(4));
            // orogen p95 = max of pure-orogen (z1) and craton∩wedge (z4) where present.
            let oro_p95 = { let mut a=oro.clone(); a.extend(cw.clone()); p95(a) };
            let crat_p95 = p95(crat.clone());
            let prom = oro_p95 / crat_p95.max(1e-6);
            let alt_inv = !crat.is_empty() && mean(&crat) < mean(&plat);
            eprintln!("    rho={:?}: craton p95={crat_p95:.3} (~{:.0}m) orogen p95={oro_p95:.3} (~{:.0}m) prom={prom:.2}  craton<platform(altInv)={alt_inv}  [z4 craton∩wedge p95={:.3}]",
                rho, to_m(crat_p95), to_m(oro_p95), p95(cw));
        }
    }
    eprintln!("  out done");
}

/// #155 B-Jordan VISUAL — the density choice + craton/orogen morphology is
/// a read, not a scalar (the rendered ratio is renorm-confounded). Render
/// 2026 (separate craton, decisive) + 1988 (craton≡margin, control) at
/// 2048² via c1_hd_production, sweeping craton_rho_crust {None,2850,2900,
/// 2950}. hypso (height, renorm-confounded across densities — read with
/// care) + hillshade (MORPHOLOGY — plateau-broad-flat vs chain-narrow — the
/// KEY view). Gated (param, not productionized).
#[test]
#[ignore]
fn probe_jordan_render() {
    let dir = output_dir().join("jordan_render");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso0 = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    let cfg = FbmUpscaleConfig::c1_hd_production(2048);
    eprintln!("#155 B-Jordan render — 2026 + 1988, 2048², density {{None,2850,2900,2950}}");
    for &seed in &[2026u64, 1988] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso0.clone(), drainage_max_distance: 30 };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        for (tag, rho) in [("none", None), ("2850", Some(2850.0f32)), ("2900", Some(2900.0)), ("2950", Some(2950.0))] {
            let mut iso = iso0.clone(); iso.craton_rho_crust = rho;
            let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
            save_heightmap01(&up.heightmap, &dir.join(format!("seed{seed:05}_rho{tag}_hypso.png")));
            save_hillshade_crop(&up.heightmap, 0, 0, up.heightmap.width, &dir.join(format!("seed{seed:05}_rho{tag}_hillshade.png")));
            eprintln!("  seed {seed} rho={tag}: written");
        }
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 VERTICAL-SCALE CONTRACT + HIGH-MOUNTAIN CEILING (read-only diagnostic).
///
/// VOLET 1 (contract): the production altitude is DOUBLE-NORMALISED with an
/// internally INCOHERENT metre scale, which is the source of the ~×2 norm→m
/// ambiguity flagged in `project_c1_vertical_scale_contract`:
///  - `compute_isostasy` outputs an ASYMMETRIC [0,1] (sea at
///    `sea_norm = max_depth_m/(max_depth_m+max_elevation_m) = 500/4500 ≈ 0.111`);
///    on LAND the metre anchor is `IsostasyResult.peak_altitude_m` =
///    `(h_max−h_sea)/(h_max−h_min)·max_elevation_m`, i.e. FIELD-RELATIVE
///    (depends on the seed's raw S̃ spread), NOT a fixed constant.
///  - Stein-Stein then OVERWRITES oceanic cells with `−depth_m/depth_scale_m`
///    (depth_scale=5000), a SECOND, FIXED metre scale with sea at 0 (not 0.111).
///  - `upscale_from_c1` then re-normalises `(v+1.13)/2.26` → land tops at 0.943,
///    S-S sea→0.5 but isostasy-continental-sea→0.549 (the two "sea levels" differ).
/// So there is NO single norm→m: land = field-relative `peak_altitude_m`,
/// ocean = fixed 5000 m, spliced at incompatible zeros. The downstream guess
/// `(v−0.5)·2·4000` used in `probe_craton_calibration` is exactly one arbitrary
/// pick — measured here against the model's OWN `peak_altitude_m`.
///
/// VOLET 2 (ceiling attribution): print, per seed, the model's OWN peak metres
/// (`peak_altitude_m`), the S̃ ceiling on the margin orogen vs EH `h_eq=2.0`,
/// and the HD rendered land-max under the downstream guess. Attributes the
/// 6000–8000 m gap: (EH) S̃ capped at 2.0 → (ramp) `peak_altitude_m` is a
/// FRACTION of `max_elevation_m=4000` → (conversion) the codomain itself tops
/// at 4000 m. NOT a fix; scopes the Phase-3 EH/critical-wedge chantier.
#[test]
#[ignore]
fn probe_vertical_scale_ceiling() {
    use ymir_core::tectonics_c1::boundary_classification::{
        classify_boundaries, oc_override_seed_mask, retarget_upper_plate_continental,
    };
    use ymir_core::tectonics_c1::distance_field::wedge_distance_intra_plate_typed;
    let iso = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    let h_eq = 2.0f64;
    let max_elev = iso.max_elevation_m;
    let p95 = |mut v: Vec<f64>| -> f64 { if v.is_empty() { return f64::NAN; } v.sort_by(|a,b| a.partial_cmp(b).unwrap()); v[(v.len()*95/100).min(v.len()-1)] };
    eprintln!("#155 vertical-scale + ceiling — max_elevation_m={max_elev}, sea_norm={:.3}, EH h_eq={h_eq}",
        iso.max_depth_m / (iso.max_depth_m + iso.max_elevation_m));
    for &seed in &[42u64, 1988, 2026] {
        // Zones (margin orogen) — same construction as probe_prominence_attribution.
        let geom = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let kin0 = PlateKinematics::preset_phase_1_1(geom.num_plates);
        let mut bi = classify_boundaries(&geom.plate_id, &kin0);
        retarget_upper_plate_continental(&mut bi, &geom.plate_id, &geom.plate_type);
        let oc_seed = oc_override_seed_mask(&bi, &geom.plate_id, &geom.plate_type);
        let (wd, is_oc) = wedge_distance_intra_plate_typed(&geom.plate_id, &bi.upper_plate_mask, &oc_seed, 30.0);

        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default(); // A′ active
        let config = C1TimeLoopConfig {
            rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0,
            iso_config: iso.clone(), drainage_max_distance: 30,
        };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});

        // S̃ ceiling on the margin orogen (raw, consistent space).
        let mut margin = Vec::new();
        let mut s_global_max = f64::NEG_INFINITY;
        for j in 0..grid { for i in 0..grid {
            let s = state.s.get(i, j);
            s_global_max = s_global_max.max(s);
            if matches!(geom.plate_type.get(i,j), PlateType::Continental)
                && is_oc.get(i,j) && wd.get(i,j) < 30.0 { margin.push(s); }
        }}
        let n_margin = margin.len();
        let n_at_ceiling = margin.iter().filter(|&&s| s >= 1.9).count();
        let margin_p95 = p95(margin.clone());

        // The model's OWN metre statement (IsostasyResult metadata).
        let iso_res = compute_isostasy_craton(&state.s, &iso, state.cratonic_mask.data());
        let peak_m = iso_res.peak_altitude_m;

        // HD rendered land-max under the DOWNSTREAM guess convention (v−0.5)·2·4000.
        let cfg = FbmUpscaleConfig::c1_hd_production(512);
        let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
        let hd_land_max_norm = up.heightmap.data.iter().cloned().fold(0.0f32, f32::max);
        let hd_guess_m = (hd_land_max_norm as f64 - 0.5) * 2.0 * max_elev as f64;

        eprintln!("  seed {seed}:");
        eprintln!("    S̃: global_max={s_global_max:.3}  margin orogen p95={margin_p95:.3}  cells≥1.9: {n_at_ceiling}/{n_margin}  (EH h_eq={h_eq})");
        eprintln!("    model peak_altitude_m = {peak_m:.0} m  (= fraction of max_elevation_m={max_elev:.0}; the ramp's land ceiling)");
        eprintln!("    HD land-max norm={hd_land_max_norm:.3} → guess-convention {hd_guess_m:.0} m  vs target 6000-8000");
    }
    eprintln!("  ATTRIBUTION: EH caps S̃ at 2.0 → peak_altitude_m is a fraction of 4000 → conversion codomain tops at 4000 < 6000-8000.");
}

/// #155 vertical-scale chantier — characterise the LAND-ceiling repair BEFORE
/// coding. The cause: compute_isostasy normalises land against global raw h_max
/// (the oceanic advective spike). Measures, on the production state: raw h by
/// plate_type, so we can choose the terrestrial ceiling (max-continental vs a
/// robust percentile) and quantify over-saturation (how many land cells would
/// clamp to 1.0 under each). Read-only — no isostasy change yet.
#[test]
#[ignore]
fn probe_land_ceiling_repair() {
    let iso = IsostasyConfig::c1_default();
    let grid = 64usize;
    let buoy = 1.0 - iso.rho_crust / iso.rho_mantle;
    let buoy_c = match iso.craton_rho_crust { Some(r) => 1.0 - r / iso.rho_mantle, None => buoy };
    let pct = |v: &[f64], q: f64| -> f64 { if v.is_empty() { return f64::NAN; } let mut s = v.to_vec(); s.sort_by(|a,b| a.partial_cmp(b).unwrap()); s[((q*(s.len()-1) as f64).round() as usize).min(s.len()-1)] };
    eprintln!("#155 land-ceiling repair — raw h by plate_type (buoy cont={buoy:.4} craton={buoy_c:.4})");
    for &seed in &[42u64, 1988, 2026] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        // raw h = S̃ * buoyancy (craton density on cratonic cells), global vs continental.
        let mut h_all = Vec::new();
        let mut h_cont = Vec::new();
        for j in 0..grid { for i in 0..grid {
            let k = j*grid+i;
            let b = if state.cratonic_mask.data()[k] { buoy_c } else { buoy };
            let h = state.s.get(i,j) * b as f64;
            h_all.push(h);
            if matches!(state.plate_type.get(i,j), PlateType::Continental) { h_cont.push(h); }
        }}
        let h_min = h_all.iter().cloned().fold(f64::INFINITY, f64::min);
        let h_max_glob = h_all.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let h_max_cont = h_cont.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        // land = h > h_sea (isostasy's own threshold), under the GLOBAL h_max ramp.
        let h_range = (pct(&h_all, 0.92) - h_min).max(1e-10); // sea uses PercentileCapped{0.92} (the LOW cap)
        let h_sea = h_min + iso.sea_level_fraction as f64 * h_range;
        let land: Vec<f64> = h_all.iter().cloned().filter(|&h| h > h_sea).collect();
        let land_cont: Vec<f64> = h_cont.iter().cloned().filter(|&h| h > h_sea).collect();
        eprintln!("  seed {seed}: h_min={h_min:.3} h_sea={h_sea:.3} | h_max GLOBAL={h_max_glob:.3} CONTINENTAL={h_max_cont:.3} (ratio {:.2})", h_max_glob/h_max_cont.max(1e-6));
        eprintln!("    land cells (h>h_sea): {} total, {} continental | land p50={:.3} p95={:.3} p99={:.3} max={:.3}",
            land.len(), land_cont.len(), pct(&land,0.5), pct(&land,0.95), pct(&land,0.99), pct(&land,1.0));
        // over-saturation: under a candidate ceiling C, land cells with h>=C clamp to 1.0.
        for &(tag, c) in &[("maxCont", h_max_cont), ("landP99", pct(&land,0.99)), ("landP95", pct(&land,0.95))] {
            let clamp = land.iter().filter(|&&h| h >= c).count();
            eprintln!("    ceiling {tag}={c:.3} → {clamp}/{} land cells clamp to 1.0 ({:.1}%)", land.len(), 100.0*clamp as f64/land.len().max(1) as f64);
        }
    }
}

/// #155 land-ceiling repair ACCEPTANCE — contrast the phantom (global h_max)
/// vs repaired (continental ceiling) peak_altitude_m, the over-saturation, and
/// the rendered land change. compute_isostasy_craton (no continental mask) =
/// pre-fix phantom; compute_isostasy_c1 (continental mask) = repaired; the
/// production upscale_from_c1 now uses the repaired path. Renders hypso +
/// hillshade for the visual re-judge (land now spreads to 1.0).
#[test]
#[ignore]
fn probe_land_ceiling_acceptance() {
    use ymir_core::tectonics::isostasy::{compute_isostasy_craton, compute_isostasy_c1};
    let dir = output_dir().join("land_ceiling_accept");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    let cfg = FbmUpscaleConfig::c1_hd_production(1024);
    eprintln!("#155 land-ceiling acceptance — phantom vs repaired peak_altitude_m + render");
    for &seed in &[42u64, 1988, 2026] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        let continental: Vec<bool> = (0..grid*grid).map(|k| matches!(state.plate_type.get(k%grid, k/grid), PlateType::Continental)).collect();
        let phantom = compute_isostasy_craton(&state.s, &iso, state.cratonic_mask.data());
        let repaired = compute_isostasy_c1(&state.s, &iso, Some(state.cratonic_mask.data()), &continental);
        // over-saturation: land cells (coarse) at exactly 1.0 in the repaired map.
        let clamp = repaired.heightmap.data.iter().filter(|&&v| v >= 0.999).count();
        let land = repaired.heightmap.data.iter().filter(|&&v| v > repaired.sea_level_normalized).count();
        // production render (repaired path).
        let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
        let hd_max = up.heightmap.data.iter().cloned().fold(0.0f32, f32::max);
        eprintln!("  seed {seed}: peak_altitude_m phantom={:.0} → repaired={:.0} m  | coarse land-clamp@1.0: {clamp}/{land} ({:.1}%)  | HD land-max norm={hd_max:.3}",
            phantom.peak_altitude_m, repaired.peak_altitude_m, 100.0*clamp as f64/land.max(1) as f64);
        save_heightmap01(&up.heightmap, &dir.join(format!("seed{seed:05}_hypso.png")));
        save_hillshade_crop(&up.heightmap, 0, 0, up.heightmap.width, &dir.join(format!("seed{seed:05}_hillshade.png")));
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 Maillon 2 ACCEPTANCE — sea-level unification + norm→m contract.
/// Checks: (1) continental sea maps to EXACTLY 0.5 (no 0.549) after the
/// production offset; (2) c1_altitude_norm_to_metres round-trips + key anchors;
/// (3) render re-judge (hypso/hillshade). Resolution invariance is preserved by
/// construction (sea-centering is on the 64² coarse, before upscale).
#[test]
#[ignore]
fn probe_sea_unification_acceptance() {
    use ymir_core::tectonics_c1::production_upscale::{
        c1_production_altitude, c1_altitude_norm_to_metres, c1_metres_to_altitude_norm,
    };
    let dir = output_dir().join("sea_unify_accept");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    let half = 1.13f32;
    let cfg = FbmUpscaleConfig::c1_hd_production(1024);
    // norm→m contract anchors (independent of any seed).
    eprintln!("#155 Maillon 2 — norm→m contract: sea(0.5)={:.1}m  norm1.0={:.0}m  norm0.0={:.0}m  roundtrip(2772m)={:.4}",
        c1_altitude_norm_to_metres(0.5, &ss), c1_altitude_norm_to_metres(1.0, &ss),
        c1_altitude_norm_to_metres(0.0, &ss), c1_metres_to_altitude_norm(2772.0, &ss));
    for &seed in &[42u64, 1988, 2026] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        // coarse production altitude (sea-centred metres/scale). Continental
        // sea-level cells should be ~0 → offset (alt+1.13)/2.26 = 0.5 exactly.
        let coarse = c1_production_altitude(&state.s, &state.age, &state.plate_type, &iso, &ss);
        // continental cells: min |alt| ≈ the coast (alt≈0); land = alt>0.
        let mut cont_land = 0usize; let mut cont = 0usize; let mut min_abs = f32::INFINITY;
        for j in 0..grid { for i in 0..grid {
            if matches!(state.plate_type.get(i,j), PlateType::Continental) {
                cont += 1;
                let a = coarse.get(i as i32, j as i32);
                min_abs = min_abs.min(a.abs());
                if a > 0.0 { cont_land += 1; }
            }
        }}
        // production HD render.
        let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
        let off = |a: f32| ((a + half) / (2.0*half)).clamp(0.0,1.0);
        let coast_norm = off(0.0); // continental sea (alt 0) under the offset
        let land_max_norm = up.heightmap.data.iter().cloned().fold(0.0f32, f32::max);
        eprintln!("  seed {seed}: continental coast alt(coarse) min|alt|={min_abs:.4} (→0) | coast→norm {coast_norm:.4} (want 0.5) | HD land-max norm={land_max_norm:.3} = {:.0}m | cont land {cont_land}/{cont}",
            c1_altitude_norm_to_metres(land_max_norm, &ss));
        save_heightmap01(&up.heightmap, &dir.join(format!("seed{seed:05}_hypso.png")));
        save_hillshade_crop(&up.heightmap, 0, 0, up.heightmap.width, &dir.join(format!("seed{seed:05}_hillshade.png")));
    }
    eprintln!("  out = {}", dir.display());
}

/// #155→complétude W7 DRAINAGE — measure what the EXISTING drainage stack
/// (terrain::flow + lakes::detection, already built) produces on the CURRENT
/// C1 HD product (unified scale, sea=0.5, post-erosion). Grounds the W7:
/// rivers by Strahler order, lakes (area/depth/level in METRES via the Maillon-2
/// contract), and the endorheic-vs-exorheic question (trace each lake outlet
/// downstream → reaches sea?). NOT the implementation — a measurement probe.
#[test]
#[ignore]
fn probe_drainage_etat_des_lieux() {
    use ymir_core::terrain::flow::{compute_flow, extract_rivers, FlowConfig, RiverConfig, DIR_NONE, D8_DX, D8_DY};
    use ymir_core::lakes::detection::{detect_lakes, LakeConfig};
    use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
    let dir = output_dir().join("drainage_etat");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let grid = 64usize;
    let target = 512usize;
    let cfg = FbmUpscaleConfig::c1_hd_production(target);
    eprintln!("#155 W7 drainage état-des-lieux — existing stack on C1 HD product ({target}², sea=0.5)");
    for &seed in &[42u64, 1988, 2026] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
        let h = &up.heightmap;
        let (w, ht) = (h.width, h.height);

        // The unified scale puts sea at 0.5.
        let flow = compute_flow(h, &FlowConfig { sea_level: 0.5 });
        let rivers = extract_rivers(&flow, &RiverConfig::default(), w, ht);
        let lakes = detect_lakes(h, &flow.filled, &flow.direction, &flow.basins, &LakeConfig::default());

        // Strahler histogram.
        let mut strahler = [0usize; 12];
        for s in &rivers.segments { strahler[(s.strahler_order as usize).min(11)] += 1; }
        let max_acc = flow.accumulation.data.iter().cloned().fold(0.0f32, f32::max);

        // Endorheic vs exorheic: trace each lake outlet downstream; reaches sea?
        let is_sea = |k: usize| h.data[k] <= 0.5;
        let mut endorheic = 0; let mut exorheic = 0;
        for lk in &lakes.lakes {
            let (mut x, mut y) = (lk.outlet.0 as usize, lk.outlet.1 as usize);
            let mut steps = 0; let mut reached_sea = false;
            loop {
                let k = y*w + x;
                if is_sea(k) { reached_sea = true; break; }
                let d = flow.direction[k];
                if d == DIR_NONE || steps > w*ht { break; }
                x = ((x as i32 + D8_DX[d as usize]).rem_euclid(w as i32)) as usize;
                y = ((y as i32 + D8_DY[d as usize]).rem_euclid(ht as i32)) as usize;
                steps += 1;
            }
            if reached_sea { exorheic += 1; } else { endorheic += 1; }
        }
        // Lake levels/depths in metres (Maillon-2 contract).
        let (mut max_area, mut deepest_m, mut highest_lake_m) = (0usize, 0.0f32, f32::MIN);
        for lk in &lakes.lakes {
            max_area = max_area.max(lk.area);
            deepest_m = deepest_m.max(c1_altitude_norm_to_metres(lk.surface_elevation, &ss) - c1_altitude_norm_to_metres(lk.surface_elevation - lk.max_depth, &ss));
            highest_lake_m = highest_lake_m.max(c1_altitude_norm_to_metres(lk.surface_elevation, &ss));
        }
        let land = h.data.iter().filter(|&&v| v > 0.5).count();
        eprintln!("  seed {seed}: land {:.1}% | basins {} max_acc {max_acc:.0} | rivers {} segs, Strahler {:?}",
            100.0*land as f64/(w*ht) as f64, flow.num_basins, rivers.segments.len(), &strahler[1..7]);
        eprintln!("    lakes {} (exorheic {exorheic} / endorheic {endorheic}) | max area {max_area} cells | deepest {deepest_m:.0} m | highest level {highest_lake_m:.0} m",
            lakes.lakes.len());

        // Render overlay (hypso + rivers blue + lakes).
        let mut img = vec![0u8; w*ht*3];
        for k in 0..w*ht {
            let v = h.data[k].clamp(0.0,1.0);
            let (r,g,b) = if v <= 0.5 { (30,60,(120.0+200.0*v) as u8) } else { let t=(v-0.5)*2.0; (((60.0+150.0*t) as u8),((120.0+80.0*t) as u8),60) };
            img[k*3]=r; img[k*3+1]=g; img[k*3+2]=b;
        }
        let stream = RiverConfig::default().stream_threshold;
        for k in 0..w*ht {
            if lakes.lake_map[k] != 0 { img[k*3]=30; img[k*3+1]=90; img[k*3+2]=180; }
            else if flow.accumulation.data[k] >= stream && h.data[k] > 0.5 { img[k*3]=40; img[k*3+1]=110; img[k*3+2]=230; }
        }
        // Origin-bottom (data row j → image row ht-1-j), matching the
        // product-wide convention (save_heightmap01/binarymap/hillshade_crop +
        // viz hydrology). The extraction reads the same GridF32 indexing as
        // every consumer; only this render's orientation was the Y-flip.
        let mut buf = image::ImageBuffer::new(w as u32, ht as u32);
        for j in 0..ht { for i in 0..w {
            let k = j*w + i;
            buf.put_pixel(i as u32, (ht - 1 - j) as u32, image::Rgb([img[k*3],img[k*3+1],img[k*3+2]]));
        }}
        buf.save(dir.join(format!("seed{seed:05}_drainage.png"))).unwrap();
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 drainage maillon ACCEPTANCE — validate the PRODUCT (c1_drainage) on the
/// real C1 terrain, not the algorithm (that has unit tests). Checks: navigability
/// classes (km² thresholds), lake stats in metres, and RESOLUTION-INDEPENDENCE
/// (512² vs 1024²: km²-anchored network density comparable, where cell-count
/// thresholds would break). Renders the drainage overlay (origin-bottom) for the
/// geometric-coherence re-judge (rivers in valleys / lakes in basins).
#[test]
#[ignore]
fn probe_c1_drainage_acceptance() {
    use ymir_core::tectonics_c1::drainage::{c1_drainage, C1DrainageConfig, Navigability, LakeType};
    use ymir_core::tectonics_c1::production_upscale::c1_cell_area_km2;
    let dir = output_dir().join("drainage_accept");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let dcfg = C1DrainageConfig::default();
    let grid = 64usize;
    eprintln!("#155 drainage acceptance — c1_drainage on C1 product, km² thresholds {:?}", dcfg.thresholds);
    for &seed in &[42u64, 1988, 2026] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});

        // Resolution-independence: 512² and 1024², compare km²-normalised metrics.
        for &target in &[512usize, 1024] {
            let cfg = FbmUpscaleConfig::c1_hd_production(target);
            let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
            let dr = c1_drainage(&up.heightmap, &dcfg, &ss);
            let (w, _h) = (dr.width, dr.height);
            let land_km2 = up.heightmap.data.iter().filter(|&&v| v > 0.5).count() as f32 * c1_cell_area_km2(w);
            let mut nav = [0usize; 4];
            for n in &dr.segment_navigability { nav[match n { Navigability::NonNavigable=>0, Navigability::SmallBoat=>1, Navigability::Barge=>2, Navigability::Ship=>3 }] += 1; }
            let max_drain = dr.segment_drainage_km2.iter().cloned().fold(0.0f32, f32::max);
            let exo = dr.lakes.iter().filter(|l| l.lake_type==LakeType::Exorheic).count();
            let endo = dr.lakes.iter().filter(|l| l.lake_type==LakeType::Endorheic).count();
            let (deepest, highest, biggest) = dr.lakes.iter().fold((0.0f32,f32::MIN,0.0f32), |(d,h,a),l| (d.max(l.depth_m), h.max(l.level_m), a.max(l.area_km2)));
            eprintln!("  seed {seed} @{target}²: land {land_km2:.0} km² | rivers {} (nonNav {} smallBoat {} barge {} ship {}) maxDrain {max_drain:.0} km²",
                dr.rivers.segments.len(), nav[0], nav[1], nav[2], nav[3]);
            eprintln!("    lakes {} (exo {exo}/endo {endo}) deepest {deepest:.0} m, highest level {highest:.0} m, biggest {biggest:.0} km²", dr.lakes.len());

            if target == 1024 {
                // render overlay (origin-bottom, the product convention). Base
                // hypso + lakes + NAVIGABLE rivers by tier (the mappable network:
                // small-boat→barge→ship), NOT every headwater (the full
                // hydrography at stream_km2 is the dense texture, consumer-filtered).
                let mut img = vec![0u8; w*dr.height*3];
                for k in 0..w*dr.height {
                    let v = up.heightmap.data[k].clamp(0.0,1.0);
                    let (r,g,b) = if v <= 0.5 { (30,60,(120.0+200.0*v) as u8) } else { let t=(v-0.5)*2.0; (((60.0+150.0*t) as u8),((120.0+80.0*t) as u8),60) };
                    img[k*3]=r; img[k*3+1]=g; img[k*3+2]=b;
                }
                for (si, seg) in dr.rivers.segments.iter().enumerate() {
                    let col = match dr.segment_navigability[si] {
                        Navigability::Ship => [20u8, 70, 200],
                        Navigability::Barge => [40, 110, 230],
                        Navigability::SmallBoat => [90, 160, 240],
                        Navigability::NonNavigable => continue,
                    };
                    for &(px, py) in &seg.points { let k = py as usize * w + px as usize; img[k*3]=col[0]; img[k*3+1]=col[1]; img[k*3+2]=col[2]; }
                }
                for k in 0..w*dr.height {
                    if dr.lake_map[k]!=0 { img[k*3]=30; img[k*3+1]=90; img[k*3+2]=180; }
                }
                let mut buf = image::ImageBuffer::new(w as u32, dr.height as u32);
                for j in 0..dr.height { for i in 0..w { let k=j*w+i; buf.put_pixel(i as u32,(dr.height-1-j) as u32, image::Rgb([img[k*3],img[k*3+1],img[k*3+2]])); }}
                buf.save(dir.join(format!("seed{seed:05}_drainage.png"))).unwrap();
            }
        }
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 drainage VISUAL validation — overlays of the (numbers-validated)
/// c1_drainage product on the relief, origin-bottom (the Y-flip lesson). Per
/// seed: (1) rivers on HILLSHADE (navigability-tiered: do they sit in valleys?
/// do navigable trunks drain interior→coast?), (2) lakes on HYPSO (are they in
/// real basins, not perched?). Same default c1_drainage config — no re-tuning,
/// this RENDERS the validated product for the eye. 1024² (+ 2048² spot-check).
#[test]
#[ignore]
fn probe_drainage_overlays() {
    use ymir_core::tectonics_c1::drainage::{c1_drainage, C1DrainageConfig, Navigability};
    let dir = output_dir().join("drainage_overlays");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let dcfg = C1DrainageConfig::default();
    let grid = 64usize;

    // hillshade grayscale for one cell (same lighting as save_hillshade_crop).
    let shade_at = |h: &GridF32, i: usize, j: usize| -> f64 {
        let z = 60.0_f64;
        let (lx, ly, lz) = { let (a,b,c)=(1.0_f64,1.0,2.0); let n=(a*a+b*b+c*c).sqrt(); (a/n,b/n,c/n) };
        let (gx, gy) = (i as i32, j as i32);
        let dzdx = (h.get(gx+1, gy) - h.get(gx-1, gy)) as f64 * z;
        let dzdy = (h.get(gx, gy+1) - h.get(gx, gy-1)) as f64 * z;
        let nn = (dzdx*dzdx + dzdy*dzdy + 1.0).sqrt();
        ((-dzdx*lx - dzdy*ly + lz)/nn).clamp(0.0, 1.0)
    };

    for &seed in &[2026u64, 1988, 42, 1337] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});

        let targets: &[usize] = &[1024];
        for &target in targets {
            let cfg = FbmUpscaleConfig::c1_hd_production(target);
            let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
            let h = &up.heightmap;
            let (w, ht) = (h.width, h.height);
            let dr = c1_drainage(h, &dcfg, &ss);

            // ---- (1) rivers on hillshade ----
            let mut riv = vec![0u8; w*ht*3];
            for j in 0..ht { for i in 0..w {
                let k = j*w+i;
                let g = if h.data[k] <= 0.5 { 70u8 } else { (shade_at(h, i, j)*255.0) as u8 };
                let c = if h.data[k] <= 0.5 { [40,70,120] } else { [g,g,g] };
                riv[k*3]=c[0]; riv[k*3+1]=c[1]; riv[k*3+2]=c[2];
            }}
            // Lakes drawn first (the filled flat basins) — rivers are then MASKED
            // under lakes (the viz convention): a flat filled basin is a lake, not
            // a blocky fill-front "river". This separates real valley rivers from
            // the D8-on-flats artifact.
            let put = |img: &mut [u8], k: usize, c: [u8;3]| { img[k*3]=c[0]; img[k*3+1]=c[1]; img[k*3+2]=c[2]; };
            for k in 0..w*ht { if dr.lake_map[k]!=0 { put(&mut riv, k, [35,80,160]); } }
            // navigable tiers + faint non-navigable streams, NOT over lake cells.
            for (si, seg) in dr.rivers.segments.iter().enumerate() {
                let (col, thick) = match dr.segment_navigability[si] {
                    Navigability::Ship => ([10u8,50,170], 2i32),
                    Navigability::Barge => ([30,90,210], 1),
                    Navigability::SmallBoat => ([80,150,235], 0),
                    Navigability::NonNavigable => ([120,170,210], 0), // faint, shows valley-following
                };
                for &(px, py) in &seg.points {
                    for dj in -thick..=thick { for di in -thick..=thick {
                        let ni=(px as i32+di).rem_euclid(w as i32) as usize; let nj=(py as i32+dj).rem_euclid(ht as i32) as usize;
                        let k = nj*w+ni;
                        if dr.lake_map[k]==0 { put(&mut riv, k, col); }
                    }}
                }
            }
            let mut buf = image::ImageBuffer::new(w as u32, ht as u32);
            for j in 0..ht { for i in 0..w { let k=j*w+i; buf.put_pixel(i as u32,(ht-1-j) as u32, image::Rgb([riv[k*3],riv[k*3+1],riv[k*3+2]])); }}
            buf.save(dir.join(format!("seed{seed:05}_{target}_rivers_hillshade.png"))).unwrap();

            // ---- (2) lakes on hypso ----
            let mut buf2 = image::ImageBuffer::new(w as u32, ht as u32);
            for j in 0..ht { for i in 0..w {
                let k=j*w+i;
                let c = if dr.lake_map[k]!=0 { [40,90,190] } else { hypsometric(h.data[k].clamp(0.0,1.0), 0.5) };
                buf2.put_pixel(i as u32,(ht-1-j) as u32, image::Rgb(c));
            }}
            buf2.save(dir.join(format!("seed{seed:05}_{target}_lakes_hypso.png"))).unwrap();

            eprintln!("  seed {seed} @{target}²: {} rivers, {} lakes → overlays written", dr.rivers.segments.len(), dr.lakes.len());
        }
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 flat-routing DIAGNOSTIC — characterise the parallel-bar / 45° fan
/// artifact on flat interiors at 2048². Confirms the root (D8 follows the
/// pit_fill epsilon-tree, NOT real terrain) by a NON-INVASIVE eps=0
/// counterfactual: steepest-descent on the ORIGINAL (unfilled) heightmap =
/// eps=0; flat cells that get a direction ONLY from the eps-fill are the
/// artifact. Splits the artifact: filled depressions (legit drainage to sill)
/// vs native plateaus (cratonic, physically diffuse — channels here are FALSE).
#[test]
#[ignore]
fn probe_flat_routing_diagnostic() {
    use ymir_core::tectonics_c1::drainage::{c1_drainage, C1DrainageConfig};
    use ymir_core::tectonics_c1::production_upscale::c1_cell_area_km2;
    use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};
    let iso = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let dcfg = C1DrainageConfig::default();
    let grid = 64usize;
    // steepest-descent direction on a given heightmap (eps=0 semantics): NONE if
    // no strictly-lower neighbour.
    let d8_eps0 = |h: &GridF32, i: usize, j: usize, w: usize, ht: usize| -> u8 {
        let my = h.data[j*w+i];
        let mut best = 0.0f32; let mut dir = DIR_NONE;
        for d in 0..8 {
            let ni = ((i as i32 + D8_DX[d]).rem_euclid(w as i32)) as usize;
            let nj = ((j as i32 + D8_DY[d]).rem_euclid(ht as i32)) as usize;
            let dist = if d % 2 == 0 { 1.0 } else { 1.414 };
            let s = (my - h.data[nj*w+ni]) / dist;
            if s > best { best = s; dir = d as u8; }
        }
        dir
    };
    let pct = |v: &[f32], q: f32| -> f32 { if v.is_empty() {return f32::NAN;} let mut s=v.to_vec(); s.sort_by(|a,b| a.partial_cmp(b).unwrap()); s[((q*(s.len()-1) as f32) as usize).min(s.len()-1)] };
    eprintln!("#155 flat-routing diagnostic — eps=0 counterfactual + filled/native split");
    for &seed in &[1988u64, 2026] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        for &target in &[1024usize, 2048] {
            let cfg = FbmUpscaleConfig::c1_hd_production(target);
            let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
            let h = &up.heightmap; let (w, ht) = (h.width, h.height);
            let dr = c1_drainage(h, &dcfg, &ss);
            let cell_km2 = c1_cell_area_km2(w);
            let stream_cells = (dcfg.thresholds.stream_km2 / cell_km2).max(1.0);

            // land gradient distribution → flat threshold (15% of median land grad).
            let mut grads = vec![0.0f32; w*ht];
            let mut land_grads = Vec::new();
            for j in 0..ht { for i in 0..w {
                let (gx, gy) = h.gradient_at(i, j); let g = (gx*gx+gy*gy).sqrt();
                grads[j*w+i] = g;
                if h.data[j*w+i] > 0.5 { land_grads.push(g); }
            }}
            let med = pct(&land_grads, 0.5);
            let flat_thr = 0.15 * med;

            let (mut land, mut flat, mut filled_flat, mut native_flat) = (0usize,0usize,0usize,0usize);
            let (mut eps_driven, mut eps_filled, mut eps_native) = (0usize,0usize,0usize);
            let (mut native_above_stream, mut native_flat_total) = (0usize, 0usize);
            for j in 0..ht { for i in 0..w {
                let k = j*w+i;
                if h.data[k] <= 0.5 { continue; }
                land += 1;
                if grads[k] >= flat_thr { continue; }
                flat += 1;
                let raised = dr.flow.filled.data[k] - h.data[k];
                let is_filled = raised > 1e-5;
                if is_filled { filled_flat += 1; } else { native_flat += 1; }
                // eps-driven: compute_flow gave a direction, eps=0 gives NONE.
                let dir_filled = dr.flow.direction[k];
                let dir_orig = d8_eps0(h, i, j, w, ht);
                if dir_filled != DIR_NONE && dir_orig == DIR_NONE {
                    eps_driven += 1;
                    if is_filled { eps_filled += 1; } else { eps_native += 1; }
                }
                // native plateau: does flow accumulation exceed the stream threshold
                // (→ FALSE channels) or stay diffuse (correct)?
                if !is_filled {
                    native_flat_total += 1;
                    if dr.flow.accumulation.data[k] >= stream_cells { native_above_stream += 1; }
                }
            }}
            eprintln!("  seed {seed} @{target}²: land {land} | flat {flat} ({:.1}% land, thr={flat_thr:.2e}) = filled {filled_flat} + native {native_flat}",
                100.0*flat as f64/land.max(1) as f64);
            eprintln!("    eps-driven (dir from fill, NONE at eps=0): {eps_driven} ({:.1}% of flat) = filled {eps_filled} + native {eps_native}",
                100.0*eps_driven as f64/flat.max(1) as f64);
            eprintln!("    native plateau flats above stream-threshold (FALSE channels): {native_above_stream}/{native_flat_total} ({:.1}%)",
                100.0*native_above_stream as f64/native_flat_total.max(1) as f64);
        }
    }
}

/// #155 flat-routing FIX acceptance — re-render the FULL accumulation network
/// (unmasked, directly comparable to the pre-fix blocky image) at 2048² on the
/// artifact seeds 1988/2026, + 1337 control (slope-dominant). After Garbrecht-
/// Martz: the parallel-bars/fans on filled flats should be replaced by
/// convergent drainage to the outlet; slopes unchanged.
#[test]
#[ignore]
fn probe_flat_routing_fix_render() {
    use ymir_core::tectonics_c1::drainage::{c1_drainage, C1DrainageConfig};
    use ymir_core::tectonics_c1::production_upscale::c1_cell_area_km2;
    let dir = output_dir().join("flat_fix");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let dcfg = C1DrainageConfig::default();
    let grid = 64usize;
    let shade = |h: &GridF32, i: usize, j: usize| -> f64 {
        let z=60.0_f64; let (lx,ly,lz)={let(a,b,c)=(1.0_f64,1.0,2.0);let nn=(a*a+b*b+c*c).sqrt();(a/nn,b/nn,c/nn)};
        let dzdx=(h.get(i as i32+1,j as i32)-h.get(i as i32-1,j as i32)) as f64*z;
        let dzdy=(h.get(i as i32,j as i32+1)-h.get(i as i32,j as i32-1)) as f64*z;
        let n=(dzdx*dzdx+dzdy*dzdy+1.0).sqrt(); ((-dzdx*lx-dzdy*ly+lz)/n).clamp(0.0,1.0)
    };
    for &(seed, target) in &[(1988u64,2048usize),(2026,2048),(1337,1024)] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        let cfg = FbmUpscaleConfig::c1_hd_production(target);
        let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
        let h = &up.heightmap; let (w, ht) = (h.width, h.height);
        let dr = c1_drainage(h, &dcfg, &ss);
        let stream = (dcfg.thresholds.stream_km2 / c1_cell_area_km2(w)).max(1.0);
        // FULL accumulation network over hillshade, UNMASKED (so bars would show).
        let mut buf = image::ImageBuffer::new(w as u32, ht as u32);
        for j in 0..ht { for i in 0..w {
            let k=j*w+i;
            let c = if h.data[k] <= 0.5 { [40,70,120] }
                else if dr.flow.accumulation.data[k] >= stream { [40,110,230] }
                else { let g=(shade(h,i,j)*255.0) as u8; [g,g,g] };
            buf.put_pixel(i as u32,(ht-1-j) as u32, image::Rgb(c));
        }}
        buf.save(dir.join(format!("seed{seed:05}_{target}_accum.png"))).unwrap();
        eprintln!("  seed {seed} @{target}²: {} rivers, {} lakes", dr.rivers.segments.len(), dr.lakes.len());
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 flat-routing fix VALIDATION — zoomed overlays on the flat interiors
/// (the decisive zone). Per seed @2048²: (1) crop on the FILLED-depression
/// interior (where bars/fans were) → must be dendritic-convergent now;
/// (2) crop on a NATIVE plateau → must stay diffuse (no invented channels);
/// (3) full lakes/hypso (invariance). 1337 = slope control. Origin-bottom,
/// rivers masked under lakes (product convention). Same c1_drainage config.
#[test]
#[ignore]
fn probe_flat_fix_validation() {
    use ymir_core::tectonics_c1::drainage::{c1_drainage, C1DrainageConfig, Navigability};
    let dir = output_dir().join("flat_fix_validation");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    let dcfg = C1DrainageConfig::default();
    let grid = 64usize;
    let shade = |h: &GridF32, i: usize, j: usize| -> u8 {
        let z=60.0_f64; let (lx,ly,lz)={let(a,b,c)=(1.0_f64,1.0,2.0);let nn=(a*a+b*b+c*c).sqrt();(a/nn,b/nn,c/nn)};
        let dzdx=(h.get(i as i32+1,j as i32)-h.get(i as i32-1,j as i32)) as f64*z;
        let dzdy=(h.get(i as i32,j as i32+1)-h.get(i as i32,j as i32-1)) as f64*z;
        let n=(dzdx*dzdx+dzdy*dzdy+1.0).sqrt(); (((-dzdx*lx-dzdy*ly+lz)/n).clamp(0.0,1.0)*255.0) as u8
    };
    // crop render: hillshade + navigable rivers (masked under lakes), origin-bottom.
    let render_crop = |h:&GridF32, dr:&ymir_core::tectonics_c1::drainage::C1DrainageResult, x0:usize,y0:usize,sz:usize, path:&std::path::Path| {
        let w=h.width;
        let mut buf=image::ImageBuffer::new(sz as u32, sz as u32);
        // base
        for jj in 0..sz { for ii in 0..sz {
            let (i,j)=(x0+ii,y0+jj); let k=j*w+i;
            let c = if h.data[k]<=0.5 {[40,70,120]} else if dr.lake_map[k]!=0 {[35,80,160]} else {let g=shade(h,i,j);[g,g,g]};
            buf.put_pixel(ii as u32,(sz-1-jj) as u32, image::Rgb(c));
        }}
        // navigable rivers (skip lake cells), thicker for higher tiers
        for (si,seg) in dr.rivers.segments.iter().enumerate() {
            let (col,th)=match dr.segment_navigability[si]{Navigability::Ship=>([10u8,50,170],1i32),Navigability::Barge=>([30,90,210],0),Navigability::SmallBoat=>([90,150,235],0),Navigability::NonNavigable=>([140,180,215],0)};
            for &(px,py) in &seg.points {
                let (px,py)=(px as usize,py as usize);
                if px<x0||py<y0||px>=x0+sz||py>=y0+sz {continue}
                for dj in -th..=th { for di in -th..=th {
                    let (ii,jj)=((px as i32-x0 as i32+di),(py as i32-y0 as i32+dj));
                    if ii<0||jj<0||ii>=sz as i32||jj>=sz as i32 {continue}
                    let k=py*w+px; if dr.lake_map[k]!=0 {continue}
                    buf.put_pixel(ii as u32,(sz as i32-1-jj) as u32, image::Rgb(col));
                }}
            }
        }
        buf.save(path).unwrap();
    };
    for &seed in &[1988u64, 2026, 1337] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        let cfg = FbmUpscaleConfig::c1_hd_production(2048);
        let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg);
        let h=&up.heightmap; let (w,ht)=(h.width,h.height);
        let dr=c1_drainage(h,&dcfg,&ss);
        let sz=512usize;
        // centroids: filled-depression flats vs native flats (land, low orig grad, not filled).
        let (mut fx,mut fy,mut fn_,mut nx,mut ny,mut nn)=(0u64,0u64,0u64,0u64,0u64,0u64);
        let mut lg=Vec::new();
        for j in 0..ht { for i in 0..w { if h.data[j*w+i]>0.5 { let (gx,gy)=h.gradient_at(i,j); lg.push((gx*gx+gy*gy).sqrt()); }}}
        lg.sort_by(|a,b|a.partial_cmp(b).unwrap()); let med=lg[lg.len()/2]; let flat_thr=0.15*med;
        for j in 0..ht { for i in 0..w { let k=j*w+i; if h.data[k]<=0.5 {continue}
            let raised=dr.flow.filled.data[k]-h.data[k];
            let (gx,gy)=h.gradient_at(i,j); let g=(gx*gx+gy*gy).sqrt();
            if raised>1e-5 { fx+=i as u64; fy+=j as u64; fn_+=1; }
            else if g<flat_thr { nx+=i as u64; ny+=j as u64; nn+=1; }
        }}
        let clamp=|c:i64,sz:usize,max:usize| (c-(sz as i64)/2).max(0).min(max as i64-sz as i64) as usize;
        if fn_>0 { let (cx,cy)=((fx/fn_) as i64,(fy/fn_) as i64); render_crop(h,&dr,clamp(cx,sz,w),clamp(cy,sz,ht),sz,&dir.join(format!("seed{seed:05}_filled_interior.png"))); }
        if nn>0 { let (cx,cy)=((nx/nn) as i64,(ny/nn) as i64); render_crop(h,&dr,clamp(cx,sz,w),clamp(cy,sz,ht),sz,&dir.join(format!("seed{seed:05}_native_plateau.png"))); }
        // full lakes/hypso (invariance)
        let mut buf=image::ImageBuffer::new(w as u32, ht as u32);
        for j in 0..ht { for i in 0..w { let k=j*w+i; let c= if dr.lake_map[k]!=0 {[40,90,190]} else {hypsometric(h.data[k].clamp(0.0,1.0),0.5)}; buf.put_pixel(i as u32,(ht-1-j) as u32, image::Rgb(c)); }}
        buf.save(dir.join(format!("seed{seed:05}_lakes_hypso_full.png"))).unwrap();
        eprintln!("  seed {seed}: filled-flat {fn_} cells, native-flat {nn} cells, {} lakes", dr.lakes.len());
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 quasi-flat residual DIAGNOSTIC — measure the EXTENT and NATURE of the
/// ladder/parallel residual on near-flat cells NOT resolved by resolve_flats
/// (non-exact-equal, so GM skipped them). The local gradient is the instrument
/// the hillshade lacked: it separates a quasi-flat defect (grad≈0) from natural
/// D8-on-planar-slope parallelism (grad low but franc). Quantify, don't fix.
#[test]
#[ignore]
fn probe_quasi_flat_residual() {
    use ymir_core::terrain::flow::{compute_flow, FlowConfig, RiverConfig, DIR_NONE, D8_DX, D8_DY};
    use ymir_core::tectonics_c1::drainage::{c1_drainage, C1DrainageConfig};
    use ymir_core::tectonics_c1::production_upscale::{c1_cell_area_km2, c1_altitude_norm_to_metres, C1_DOMAIN_KM};
    let iso = IsostasyConfig::c1_default();
    let ss = SteinSteinParams::default();
    // gradient bins in norm/cell; converted to m/km via the contracts.
    let bins = [0.0f32, 1e-5, 1e-4, 1e-3, 1e-2, f32::INFINITY];
    let binname = ["<1e-5","1e-5..1e-4","1e-4..1e-3","1e-3..1e-2",">=1e-2"];
    // local directional coherence: fraction of drainage-neighbours (5x5) sharing dir.
    let coherence = |dir:&[u8], acc:&GridF32, stream:f32, w:usize, ht:usize, k:usize| -> f32 {
        let d=dir[k]; if d==DIR_NONE {return 0.0}
        let (ci,cj)=(k%w,k/w); let (mut same,mut tot)=(0u32,0u32);
        for dj in -2i32..=2 { for di in -2i32..=2 {
            if di==0&&dj==0 {continue}
            let ni=(ci as i32+di).rem_euclid(w as i32) as usize; let nj=(cj as i32+dj).rem_euclid(ht as i32) as usize;
            let m=nj*w+ni; if acc.data[m]>=stream && dir[m]!=DIR_NONE { tot+=1; if dir[m]==d {same+=1} }
        }}
        if tot==0 {0.0} else {same as f32/tot as f32}
    };
    // --- counterfactual: synthetic uniform slopes (parallel-on-planar baseline) ---
    eprintln!("#155 quasi-flat residual — synthetic uniform-slope counterfactual (parallel baseline):");
    for &sg in &[1e-4f32, 1e-3, 5e-3] {
        let n=512usize; let mut hm=GridF32::new(n,n,0.0);
        for j in 0..n { for i in 0..n { hm.set(i,j, (0.75 - i as f32*sg).max(0.2)); }}
        let fl=compute_flow(&hm,&FlowConfig{sea_level:0.5});
        let stream=RiverConfig::default().stream_threshold;
        let (mut cs,mut cn)=(0.0f32,0u32);
        for k in 0..n*n { if hm.data[k]>0.5 && fl.accumulation.data[k]>=stream { cs+=coherence(&fl.direction,&fl.accumulation,stream,n,n,k); cn+=1; }}
        eprintln!("  slope {sg:.0e} norm/cell (~{:.1} m/km): mean dir-coherence over drainage = {:.2} (1.0 = fully parallel)",
            sg*c1_altitude_norm_to_metres_delta(&ss)/(C1_DOMAIN_KM/n as f32), if cn>0 {cs/cn as f32} else {0.0});
    }
    // --- real terrain ---
    let grid=64usize; let dcfg=C1DrainageConfig::default();
    for &seed in &[1988u64, 2026] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let closures = C1Closures::default();
        let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
        run_with_closures(&mut state, &mut kin, &config, &closures, |_, _| {});
        let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &cfg_2048());
        let h=&up.heightmap; let (w,ht)=(h.width,h.height);
        let dr=c1_drainage(h,&dcfg,&ss);
        let f=&dr.flow.filled.data;
        let stream=(dcfg.thresholds.stream_km2/c1_cell_area_km2(w)).max(1.0);
        let m_per_km = c1_altitude_norm_to_metres_delta(&ss)/(C1_DOMAIN_KM/w as f32);
        // recompute exact-flat (GM-resolved) membership on filled.
        let nb=|i:usize,j:usize,d:usize|{let ni=((i as i32+D8_DX[d]).rem_euclid(w as i32))as usize;let nj=((j as i32+D8_DY[d]).rem_euclid(ht as i32))as usize;nj*w+ni};
        let mut bin_tot=[0u64;5]; let mut bin_nonexact=[0u64;5]; let mut bin_coh=[0.0f64;5];
        let (mut drain,mut land)=(0u64,0u64);
        let (mut resid,mut resid_fringe)=(0u64,0u64);
        for j in 0..ht { for i in 0..w { let k=j*w+i;
            if h.data[k]<=0.5 {continue} land+=1;
            if dr.flow.accumulation.data[k] < stream {continue} drain+=1;
            let (gx,gy)=h.gradient_at(i,j); let g=(gx*gx+gy*gy).sqrt();
            let mut b=4; for t in 0..5 { if g>=bins[t]&&g<bins[t+1]{b=t;break} }
            bin_tot[b]+=1; bin_coh[b]+=coherence(&dr.flow.direction,&dr.flow.accumulation,stream,w,ht,k) as f64;
            // exact-flat? no strictly-lower filled neighbour AND an exact-equal one.
            let (mut hl,mut he)=(false,false);
            for d in 0..8 { let mm=nb(i,j,d); if f[mm]<f[k]{hl=true} else if f[mm]==f[k]{he=true} }
            let exact_flat = !hl && he;
            if !exact_flat { bin_nonexact[b]+=1; }
            // residual band = low-gradient (<1e-4) AND not exact-flat (GM skipped it).
            if g<1e-4 && !exact_flat {
                resid+=1;
                // fringe: adjacent to a filled cell (filled>h → partially-filled depression edge)?
                let mut fr=false; for d in 0..8 { let mm=nb(i,j,d); if dr.flow.filled.data[mm]-h.data[mm]>1e-5 {fr=true} }
                if fr {resid_fringe+=1;}
            }
        }}
        eprintln!("  seed {seed} @2048² (m/km factor {m_per_km:.0}): drainage cells {drain} ({:.1}% of land)", 100.0*drain as f64/land.max(1) as f64);
        for b in 0..5 {
            let coh = if bin_tot[b]>0 {bin_coh[b]/bin_tot[b] as f64} else {0.0};
            eprintln!("    grad {:>10} ({:>6.1} m/km lo): {:>7} drainage ({:>4.1}%), non-exact {:>7} ({:>4.1}%), coherence {:.2}",
                binname[b], bins[b]*m_per_km, bin_tot[b], 100.0*bin_tot[b] as f64/drain.max(1) as f64,
                bin_nonexact[b], 100.0*bin_nonexact[b] as f64/bin_tot[b].max(1) as f64, coh);
        }
        eprintln!("    RESIDUAL band (grad<1e-4 & non-exact-flat): {resid} cells = {:.2}% of drainage, {:.3}% of land | fringe-of-depression {:.0}%",
            100.0*resid as f64/drain.max(1) as f64, 100.0*resid as f64/land.max(1) as f64, 100.0*resid_fringe as f64/resid.max(1) as f64);
    }
}
fn cfg_2048() -> FbmUpscaleConfig { FbmUpscaleConfig::c1_hd_production(2048) }
fn c1_altitude_norm_to_metres_delta(ss:&SteinSteinParams)->f32 {
    use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
    c1_altitude_norm_to_metres(1.0,ss)-c1_altitude_norm_to_metres(0.0,ss)
}

/// #155 ladder-residual LOCATOR — confirm the mechanism on the WORST flagged
/// feature before any fix. Finds the largest cluster of (drainage, grad<1e-4,
/// non-exact-flat, high coherence) cells, hard-zooms it (96px) coloring each
/// drainage cell by D8 DIRECTION (8 hues) over faint hillshade: a rectilinear
/// ladder = regular same-hue blocks; correct convergent drainage = mixed hues.
/// Reports those cells: gradient, raised-by-fill?, adjacent-to-resolved-flat?.
#[test]
#[ignore]
fn probe_ladder_locator() {
    use ymir_core::terrain::flow::{DIR_NONE, D8_DX, D8_DY};
    use ymir_core::tectonics_c1::drainage::{c1_drainage, C1DrainageConfig};
    use ymir_core::tectonics_c1::production_upscale::c1_cell_area_km2;
    let dir = output_dir().join("ladder_locator");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso = IsostasyConfig::c1_default(); let ss = SteinSteinParams::default();
    let dcfg = C1DrainageConfig::default(); let grid=64usize;
    // 8 distinct hues for D8 directions.
    let huecol = |d:u8| -> [u8;3] { match d {0=>[230,30,30],1=>[230,140,20],2=>[220,220,30],3=>[60,200,40],4=>[30,200,200],5=>[40,90,230],6=>[150,40,220],7=>[230,40,150],_=>[120,120,120]} };
    let shade=|h:&GridF32,i:usize,j:usize|->u8{let z=40.0_f64;let(lx,ly,lz)={let(a,b,c)=(1.0_f64,1.0,2.0);let nn=(a*a+b*b+c*c).sqrt();(a/nn,b/nn,c/nn)};let dx=(h.get(i as i32+1,j as i32)-h.get(i as i32-1,j as i32))as f64*z;let dy=(h.get(i as i32,j as i32+1)-h.get(i as i32,j as i32-1))as f64*z;let n=(dx*dx+dy*dy+1.0).sqrt();(((-dx*lx-dy*ly+lz)/n).clamp(0.0,1.0)*200.0+30.0)as u8};
    for &seed in &[1988u64, 2026] {
        let mut state=init_c1_state_phase_2_r7(grid,seed,&Phase2InitParams::default());
        let mut kin=PlateKinematics::preset_phase_1_1(state.num_plates);
        let config=C1TimeLoopConfig{rigid_continental_crust:true,n_steps:300,dx:1.0/64.0,dy:1.0/64.0,iso_config:iso.clone(),drainage_max_distance:30};
        run_with_closures(&mut state,&mut kin,&config,&C1Closures::default(),|_,_|{});
        let up=upscale_from_c1(&state,&iso,&ss,&WorldSeed::new(seed),&FbmUpscaleConfig::c1_hd_production(2048));
        let h=&up.heightmap; let(w,ht)=(h.width,h.height); let dr=c1_drainage(h,&dcfg,&ss);
        let f=&dr.flow.filled.data; let dirf=&dr.flow.direction; let acc=&dr.flow.accumulation;
        let stream=(dcfg.thresholds.stream_km2/c1_cell_area_km2(w)).max(1.0);
        let nb=|i:usize,j:usize,d:usize|{let ni=((i as i32+D8_DX[d]).rem_euclid(w as i32))as usize;let nj=((j as i32+D8_DY[d]).rem_euclid(ht as i32))as usize;nj*w+ni};
        let coh=|k:usize|->f32{let d=dirf[k];if d==DIR_NONE{return 0.0}let(ci,cj)=(k%w,k/w);let(mut s,mut t)=(0u32,0u32);for dj in -2i32..=2{for di in -2i32..=2{if di==0&&dj==0{continue}let ni=(ci as i32+di).rem_euclid(w as i32)as usize;let nj=(cj as i32+dj).rem_euclid(ht as i32)as usize;let m=nj*w+ni;if acc.data[m]>=stream&&dirf[m]!=DIR_NONE{t+=1;if dirf[m]==d{s+=1}}}}if t==0{0.0}else{s as f32/t as f32}};
        // suspect mask: drainage, low grad, non-exact-flat, high coherence.
        let mut suspect=vec![false;w*ht];
        for j in 0..ht{for i in 0..w{let k=j*w+i;if h.data[k]<=0.5||acc.data[k]<stream{continue}
            let(gx,gy)=h.gradient_at(i,j);let g=(gx*gx+gy*gy).sqrt();if g>=1e-4{continue}
            let(mut hl,mut he)=(false,false);for d in 0..8{let m=nb(i,j,d);if f[m]<f[k]{hl=true}else if f[m]==f[k]{he=true}}
            let exact=!hl&&he; if exact{continue}
            if coh(k)>=0.7{suspect[k]=true}
        }}
        // largest suspect cluster (BFS), get centroid.
        let mut seen=vec![false;w*ht]; let(mut best_n,mut bx,mut by)=(0usize,0usize,0usize);
        for s0 in 0..w*ht{if !suspect[s0]||seen[s0]{continue}let mut q=std::collections::VecDeque::new();q.push_back(s0);seen[s0]=true;let(mut cnt,mut sx,mut sy)=(0u64,0u64,0u64);
            while let Some(c)=q.pop_front(){cnt+=1;sx+=(c%w)as u64;sy+=(c/w)as u64;let(ci,cj)=(c%w,c/w);for d in 0..8{let m=nb(ci,cj,d);if suspect[m]&&!seen[m]{seen[m]=true;q.push_back(m)}}}
            if cnt as usize>best_n{best_n=cnt as usize;bx=(sx/cnt)as usize;by=(sy/cnt)as usize}}
        let total_suspect:usize=suspect.iter().filter(|&&x|x).count();
        eprintln!("  seed {seed}: suspect cells {total_suspect} ({:.3}% land-ish), largest cluster {best_n} cells @ ({bx},{by})",100.0*total_suspect as f64/(w*ht)as f64);
        if best_n==0 {eprintln!("    no high-coherence low-grad non-exact cluster → ladder NOT a coherent parallel artifact"); continue;}
        // hard zoom 96px on the worst cluster, color drainage by D8 dir.
        let sz=96usize; let x0=(bx as i64-48).max(0).min(w as i64-sz as i64)as usize; let y0=(by as i64-48).max(0).min(ht as i64-sz as i64)as usize;
        let mut buf=image::ImageBuffer::new(sz as u32,sz as u32);
        let(mut raised,mut adj_flat)=(0u64,0u64); let mut cnt=0u64;
        for jj in 0..sz{for ii in 0..sz{let(i,j)=(x0+ii,y0+jj);let k=j*w+i;
            let c= if h.data[k]<=0.5 {[40,70,120]}
                else if acc.data[k]>=stream && dirf[k]!=DIR_NONE { if suspect[k]{cnt+=1; if f[k]-h.data[k]>1e-5{raised+=1} let mut af=false;for d in 0..8{let m=nb(i,j,d);let(mut hl2,mut he2)=(false,false);for e in 0..8{let mm=nb(m%w,m/w,e);if f[mm]<f[m]{hl2=true}else if f[mm]==f[m]{he2=true}}if !hl2&&he2{af=true}}if af{adj_flat+=1}} huecol(dirf[k]) }
                else {let g=shade(h,i,j);[g,g,g]};
            buf.put_pixel(ii as u32,(sz-1-jj)as u32,image::Rgb(c));
        }}
        buf.save(dir.join(format!("seed{seed:05}_ladder_dirhue.png"))).unwrap();
        eprintln!("    zoom @({x0},{y0}) 96px: suspect-in-crop {cnt}, raised-by-fill {raised}, adjacent-to-resolved-flat {adj_flat}");
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 ladder reconcile — render the two candidate identities of the flagged
/// "triangle" so it can be attributed: LAKES (cyan, outlined) vs GM-RESOLVED
/// flats / their convergent fans (yellow tint) vs rivers (blue), over hillshade,
/// 2048², origin-bottom. The triangle's colour identifies it.
#[test]
#[ignore]
fn probe_ladder_reconcile() {
    use ymir_core::terrain::flow::{DIR_NONE, D8_DX, D8_DY};
    use ymir_core::tectonics_c1::drainage::{c1_drainage, C1DrainageConfig};
    use ymir_core::tectonics_c1::production_upscale::c1_cell_area_km2;
    let dir = output_dir().join("ladder_reconcile");
    std::fs::create_dir_all(&dir).expect("create dir");
    let iso = IsostasyConfig::c1_default(); let ss = SteinSteinParams::default();
    let dcfg = C1DrainageConfig::default(); let grid=64usize;
    let shade=|h:&GridF32,i:usize,j:usize|->u8{let z=50.0_f64;let(lx,ly,lz)={let(a,b,c)=(1.0_f64,1.0,2.0);let nn=(a*a+b*b+c*c).sqrt();(a/nn,b/nn,c/nn)};let dx=(h.get(i as i32+1,j as i32)-h.get(i as i32-1,j as i32))as f64*z;let dy=(h.get(i as i32,j as i32+1)-h.get(i as i32,j as i32-1))as f64*z;let n=(dx*dx+dy*dy+1.0).sqrt();(((-dx*lx-dy*ly+lz)/n).clamp(0.0,1.0)*200.0+40.0)as u8};
    for &seed in &[1988u64, 2026] {
        let mut state=init_c1_state_phase_2_r7(grid,seed,&Phase2InitParams::default());
        let mut kin=PlateKinematics::preset_phase_1_1(state.num_plates);
        let config=C1TimeLoopConfig{rigid_continental_crust:true,n_steps:300,dx:1.0/64.0,dy:1.0/64.0,iso_config:iso.clone(),drainage_max_distance:30};
        run_with_closures(&mut state,&mut kin,&config,&C1Closures::default(),|_,_|{});
        let up=upscale_from_c1(&state,&iso,&ss,&WorldSeed::new(seed),&FbmUpscaleConfig::c1_hd_production(2048));
        let h=&up.heightmap;let(w,ht)=(h.width,h.height);let dr=c1_drainage(h,&dcfg,&ss);
        let f=&dr.flow.filled.data; let acc=&dr.flow.accumulation;
        let stream=(dcfg.thresholds.stream_km2/c1_cell_area_km2(w)).max(1.0);
        let nb=|i:usize,j:usize,d:usize|{let ni=((i as i32+D8_DX[d]).rem_euclid(w as i32))as usize;let nj=((j as i32+D8_DY[d]).rem_euclid(ht as i32))as usize;nj*w+ni};
        // exact-flat (GM-resolved) membership.
        let mut exact=vec![false;w*ht];
        for j in 0..ht{for i in 0..w{let k=j*w+i;if h.data[k]<=0.5{continue}let(mut hl,mut he)=(false,false);for d in 0..8{let m=nb(i,j,d);if f[m]<f[k]{hl=true}else if f[m]==f[k]{he=true}}exact[k]=!hl&&he;}}
        let mut buf=image::ImageBuffer::new(w as u32,ht as u32);
        for j in 0..ht{for i in 0..w{let k=j*w+i;
            let lake=dr.lake_map[k]!=0;
            // lake outline: lake cell adjacent to non-lake.
            let lake_edge= lake && (0..8).any(|d|{let m=nb(i,j,d);dr.lake_map[m]==0});
            let c= if h.data[k]<=0.5 {[25,45,90]}
                else if lake_edge {[180,240,255]}
                else if lake {[60,200,220]}
                else if acc.data[k]>=stream {[40,110,235]}      // rivers
                else if exact[k] {let g=shade(h,i,j);[ (g as u16*9/10+60) as u8, (g as u16*9/10+50) as u8, (g/3) as u8 ]} // resolved-flat = yellow tint
                else {let g=shade(h,i,j);[g,g,g]};
            buf.put_pixel(i as u32,(ht-1-j)as u32,image::Rgb(c));
        }}
        buf.save(dir.join(format!("seed{seed:05}_reconcile.png"))).unwrap();
        let nlake:usize=dr.lake_map.iter().filter(|&&x|x!=0).count();
        let nexact:usize=exact.iter().filter(|&&x|x).count();
        eprintln!("  seed {seed}: lake cells {nlake} ({:.2}%), resolved-flat cells {nexact} ({:.2}%) [yellow], rivers blue", 100.0*nlake as f64/(w*ht)as f64, 100.0*nexact as f64/(w*ht)as f64);
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 yellow-fan mechanism — zoom the largest NON-LAKE resolved-flat (the
/// "yellow triangle"), colour by D8 direction, to see WHY it's rectilinear
/// (diamond iso-lines of the tl-dominant distance transform?). Confirm before
/// re-tuning the GM weighting.
#[test]
#[ignore]
fn probe_yellow_fan_mechanism() {
    use ymir_core::terrain::flow::{DIR_NONE, D8_DX, D8_DY};
    use ymir_core::tectonics_c1::drainage::{c1_drainage, C1DrainageConfig};
    let dir = output_dir().join("yellow_fan"); std::fs::create_dir_all(&dir).unwrap();
    let iso=IsostasyConfig::c1_default(); let ss=SteinSteinParams::default(); let dcfg=C1DrainageConfig::default(); let grid=64usize;
    let huecol=|d:u8|->[u8;3]{match d{0=>[230,30,30],1=>[230,140,20],2=>[220,220,30],3=>[60,200,40],4=>[30,200,200],5=>[40,90,230],6=>[150,40,220],7=>[230,40,150],_=>[110,110,110]}};
    for &seed in &[1988u64,2026]{
        let mut state=init_c1_state_phase_2_r7(grid,seed,&Phase2InitParams::default());
        let mut kin=PlateKinematics::preset_phase_1_1(state.num_plates);
        let config=C1TimeLoopConfig{rigid_continental_crust:true,n_steps:300,dx:1.0/64.0,dy:1.0/64.0,iso_config:iso.clone(),drainage_max_distance:30};
        run_with_closures(&mut state,&mut kin,&config,&C1Closures::default(),|_,_|{});
        let up=upscale_from_c1(&state,&iso,&ss,&WorldSeed::new(seed),&FbmUpscaleConfig::c1_hd_production(2048));
        let h=&up.heightmap;let(w,ht)=(h.width,h.height);let dr=c1_drainage(h,&dcfg,&ss);
        let f=&dr.flow.filled.data; let dirf=&dr.flow.direction;
        let nb=|i:usize,j:usize,d:usize|{let ni=((i as i32+D8_DX[d]).rem_euclid(w as i32))as usize;let nj=((j as i32+D8_DY[d]).rem_euclid(ht as i32))as usize;nj*w+ni};
        // exact-flat NON-lake cells.
        let mut ef=vec![false;w*ht];
        for j in 0..ht{for i in 0..w{let k=j*w+i;if h.data[k]<=0.5||dr.lake_map[k]!=0{continue}let(mut hl,mut he)=(false,false);for d in 0..8{let m=nb(i,j,d);if f[m]<f[k]{hl=true}else if f[m]==f[k]{he=true}}ef[k]=!hl&&he;}}
        // largest connected ef cluster.
        let mut seen=vec![false;w*ht];let(mut bn,mut bx,mut by,mut bw,mut bh,mut bi,mut bj)=(0usize,0usize,0usize,0usize,0usize,0usize,0usize);
        for s0 in 0..w*ht{if !ef[s0]||seen[s0]{continue}let mut q=std::collections::VecDeque::new();q.push_back(s0);seen[s0]=true;let(mut c,mut sx,mut sy,mut mni,mut mxi,mut mnj,mut mxj)=(0u64,0u64,0u64,w,0usize,ht,0usize);
            while let Some(z)=q.pop_front(){c+=1;let(zi,zj)=(z%w,z/w);sx+=zi as u64;sy+=zj as u64;mni=mni.min(zi);mxi=mxi.max(zi);mnj=mnj.min(zj);mxj=mxj.max(zj);for d in 0..8{let m=nb(zi,zj,d);if ef[m]&&!seen[m]{seen[m]=true;q.push_back(m)}}}
            if c as usize>bn{bn=c as usize;bx=(sx/c)as usize;by=(sy/c)as usize;bw=mxi-mni+1;bh=mxj-mnj+1;bi=mni;bj=mnj;}}
        eprintln!("  seed {seed}: largest non-lake resolved-flat = {bn} cells, bbox {bw}x{bh} @({bi},{bj}), centroid ({bx},{by})");
        let sz=128usize;let x0=(bx as i64-64).max(0).min(w as i64-sz as i64)as usize;let y0=(by as i64-64).max(0).min(ht as i64-sz as i64)as usize;
        let mut buf=image::ImageBuffer::new(sz as u32,sz as u32);
        for jj in 0..sz{for ii in 0..sz{let(i,j)=(x0+ii,y0+jj);let k=j*w+i;
            let c= if h.data[k]<=0.5{[30,50,90]} else if dr.lake_map[k]!=0{[60,200,220]} else if ef[k]{huecol(dirf[k])} else {[235,235,235]};
            buf.put_pixel(ii as u32,(sz-1-jj)as u32,image::Rgb(c));
        }}
        buf.save(dir.join(format!("seed{seed:05}_yellowfan_dir.png"))).unwrap();
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 non-lake flooded zones DIAGNOSTIC — characterise the residual yellow
/// (flooded to sill, filled>eroded, but NOT in any named lake → depressions
/// whose max depth never reached the 10 m lake-naming threshold). Measure size
/// (km², natural cutoff vs continuum?), depth (m, real water vs films?), shape
/// (compact basin vs filiform artifact), terrain fraction; + marked maps.
#[test]
#[ignore]
fn probe_nonlake_flooded_zones() {
    use ymir_core::tectonics_c1::drainage::{c1_drainage, C1DrainageConfig};
    use ymir_core::tectonics_c1::production_upscale::{c1_cell_area_km2, c1_altitude_norm_to_metres};
    let dir = output_dir().join("nonlake_flooded"); std::fs::create_dir_all(&dir).unwrap();
    let iso=IsostasyConfig::c1_default(); let ss=SteinSteinParams::default(); let dcfg=C1DrainageConfig::default(); let grid=64usize;
    let m_per_norm = c1_altitude_norm_to_metres(1.0,&ss)-c1_altitude_norm_to_metres(0.0,&ss);
    for &seed in &[1988u64,2026,42,1337]{
        let mut state=init_c1_state_phase_2_r7(grid,seed,&Phase2InitParams::default());
        let mut kin=PlateKinematics::preset_phase_1_1(state.num_plates);
        let config=C1TimeLoopConfig{rigid_continental_crust:true,n_steps:300,dx:1.0/64.0,dy:1.0/64.0,iso_config:iso.clone(),drainage_max_distance:30};
        run_with_closures(&mut state,&mut kin,&config,&C1Closures::default(),|_,_|{});
        let up=upscale_from_c1(&state,&iso,&ss,&WorldSeed::new(seed),&FbmUpscaleConfig::c1_hd_production(2048));
        let h=&up.heightmap;let(w,ht)=(h.width,h.height);let dr=c1_drainage(h,&dcfg,&ss);
        let cell_km2=c1_cell_area_km2(w);
        let depth=|k:usize|->f32{(dr.flow.filled.data[k]-h.data[k]).max(0.0)};
        // residual flooded non-lake cells.
        let mut res=vec![false;w*ht]; let mut nres=0u64;
        for k in 0..w*ht{ if h.data[k]>0.5 && dr.lake_map[k]==0 && depth(k)>1e-6 {res[k]=true;nres+=1;} }
        // connected components.
        let mut seen=vec![false;w*ht]; let mut comps:Vec<(u64,f32,f32,usize,usize,usize,usize)>=Vec::new(); // area,maxd_norm,sumd,mni,mxi,mnj,mxj
        for s0 in 0..w*ht{ if !res[s0]||seen[s0]{continue}
            let mut q=std::collections::VecDeque::new();q.push_back(s0);seen[s0]=true;
            let(mut a,mut md,mut sd,mut mni,mut mxi,mut mnj,mut mxj)=(0u64,0.0f32,0.0f32,w,0usize,ht,0usize);
            while let Some(c)=q.pop_front(){a+=1;let d=depth(c);md=md.max(d);sd+=d;let(ci,cj)=(c%w,c/w);mni=mni.min(ci);mxi=mxi.max(ci);mnj=mnj.min(cj);mxj=mxj.max(cj);
                for dy in -1i32..=1{for dx in -1i32..=1{if dx==0&&dy==0{continue}let ni=((ci as i32+dx).rem_euclid(w as i32))as usize;let nj=((cj as i32+dy).rem_euclid(ht as i32))as usize;let m=nj*w+ni;if res[m]&&!seen[m]{seen[m]=true;q.push_back(m)}}}}
            comps.push((a,md,sd,mni,mxi,mnj,mxj));
        }
        comps.sort_by(|x,y|y.0.cmp(&x.0));
        // size histogram (km²) + depth-of-max histogram (m) + shape.
        let szbins=[0.0,1.0,10.0,100.0,1000.0,f32::INFINITY]; let mut szh=[0u64;5]; let mut szarea=[0.0f64;5];
        let dpbins=[0.0,1.0,3.0,6.0,10.0,f32::INFINITY]; let mut dph=[0u64;5];
        for &(a,md,_,mni,mxi,mnj,mxj) in &comps{
            let km2=a as f32*cell_km2; let mdm=md*m_per_norm;
            for t in 0..5{if km2>=szbins[t]&&km2<szbins[t+1]{szh[t]+=1;szarea[t]+=a as f64*cell_km2 as f64;break}}
            for t in 0..5{if mdm>=dpbins[t]&&mdm<dpbins[t+1]{dph[t]+=1;break}}
            let _=(mni,mxi,mnj,mxj);
        }
        let land:u64=(0..w*ht).filter(|&k|h.data[k]>0.5).count() as u64;
        let lake_cells:u64=dr.lake_map.iter().filter(|&&x|x!=0).count() as u64;
        eprintln!("  seed {seed} @2048²: residual non-lake flooded {nres} cells ({:.2}% land, {:.0} km²) in {} components | named-lake cells {} ({:.2}% land)",
            100.0*nres as f64/land.max(1)as f64, nres as f64*cell_km2 as f64, comps.len(), lake_cells, 100.0*lake_cells as f64/land.max(1)as f64);
        eprintln!("    size km² [<1 |1-10|10-100|100-1k|>1k]: count {:?} area {:?}", szh, szarea.iter().map(|x|*x as u64).collect::<Vec<_>>());
        eprintln!("    max-depth m [<1|1-3|3-6|6-10|>=10]: {:?}", dph);
        // top 5 components: area km², max depth m, fill ratio (compact vs filiform).
        for (r,&(a,md,sd,mni,mxi,mnj,mxj)) in comps.iter().take(5).enumerate(){
            let bbox=((mxi-mni+1)*(mxj-mnj+1)) as f32; let fill=a as f32/bbox.max(1.0);
            eprintln!("    #{r}: {:.1} km², max {:.1} m, mean {:.1} m, fill-ratio {:.2} ({}x{})", a as f32*cell_km2, md*m_per_norm, (sd/a as f32)*m_per_norm, fill, mxi-mni+1, mxj-mnj+1);
        }
        // marked map: hypso + residual in magenta.
        let mut buf=image::ImageBuffer::new(w as u32,ht as u32);
        for j in 0..ht{for i in 0..w{let k=j*w+i; let c= if res[k]{[230,40,200]} else if dr.lake_map[k]!=0{[40,90,190]} else {hypsometric(h.data[k].clamp(0.0,1.0),0.5)}; buf.put_pixel(i as u32,(ht-1-j)as u32,image::Rgb(c));}}
        buf.save(dir.join(format!("seed{seed:05}_nonlake_marked.png"))).unwrap();
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 clean PRODUCT render @2048² (no diagnostic coloring) — to check whether
/// the residual non-lake flooded speckle is even VISIBLE in the product (where
/// it renders as terrain, not lake/river), vs only in the diagnostic marking.
/// Side-by-side with the marked map (probe_nonlake_flooded_zones) settles it.
#[test]
#[ignore]
fn probe_clean_product_2048() {
    use ymir_core::tectonics_c1::drainage::{c1_drainage, C1DrainageConfig};
    use ymir_core::tectonics_c1::production_upscale::c1_cell_area_km2;
    let dir = output_dir().join("clean_product"); std::fs::create_dir_all(&dir).unwrap();
    let iso=IsostasyConfig::c1_default(); let ss=SteinSteinParams::default(); let dcfg=C1DrainageConfig::default(); let grid=64usize;
    let shade=|h:&GridF32,i:usize,j:usize|->u8{let z=55.0_f64;let(lx,ly,lz)={let(a,b,c)=(1.0_f64,1.0,2.0);let nn=(a*a+b*b+c*c).sqrt();(a/nn,b/nn,c/nn)};let dx=(h.get(i as i32+1,j as i32)-h.get(i as i32-1,j as i32))as f64*z;let dy=(h.get(i as i32,j as i32+1)-h.get(i as i32,j as i32-1))as f64*z;let n=(dx*dx+dy*dy+1.0).sqrt();(((-dx*lx-dy*ly+lz)/n).clamp(0.0,1.0)*255.0)as u8};
    for &seed in &[1988u64,2026]{
        let mut state=init_c1_state_phase_2_r7(grid,seed,&Phase2InitParams::default());
        let mut kin=PlateKinematics::preset_phase_1_1(state.num_plates);
        let config=C1TimeLoopConfig{rigid_continental_crust:true,n_steps:300,dx:1.0/64.0,dy:1.0/64.0,iso_config:iso.clone(),drainage_max_distance:30};
        run_with_closures(&mut state,&mut kin,&config,&C1Closures::default(),|_,_|{});
        let up=upscale_from_c1(&state,&iso,&ss,&WorldSeed::new(seed),&FbmUpscaleConfig::c1_hd_production(2048));
        let h=&up.heightmap;let(w,ht)=(h.width,h.height);let dr=c1_drainage(h,&dcfg,&ss);
        let stream=(dcfg.thresholds.stream_km2/c1_cell_area_km2(w)).max(1.0);
        // PRODUCT convention: ocean, lakes, rivers (acc>=stream masked under lakes), else hillshade.
        let mut buf=image::ImageBuffer::new(w as u32,ht as u32);
        for j in 0..ht{for i in 0..w{let k=j*w+i;
            let c= if h.data[k]<=0.5 {[30,60,(120.0+200.0*h.data[k]) as u8]}
                else if dr.lake_map[k]!=0 {[40,90,190]}
                else if dr.flow.accumulation.data[k]>=stream {[40,110,230]}
                else {let g=shade(h,i,j);[g,g,g]};
            buf.put_pixel(i as u32,(ht-1-j)as u32,image::Rgb(c));
        }}
        buf.save(dir.join(format!("seed{seed:05}_clean_product.png"))).unwrap();
        eprintln!("  seed {seed}: clean product written");
    }
    eprintln!("  out = {}", dir.display());
}

/// #155 elucidate the 0.78 orogen-render fraction: is it a normalization coupling
/// (peak_altitude_m ties land height to the GLOBAL range incl ocean h_min) or
/// correct? Reconstruct the isostasy raw-h internals at 64² (the sim grid).
#[test]
#[ignore]
fn probe_orogen_fraction() {
    let iso = IsostasyConfig::c1_default(); let grid=64usize;
    let buoy=1.0-iso.rho_crust/iso.rho_mantle;
    let buoy_c=match iso.craton_rho_crust { Some(r)=>1.0-r/iso.rho_mantle, None=>buoy };
    let pct=|v:&[f32],q:f32|->f32{let mut s=v.to_vec();s.sort_by(|a,b|a.partial_cmp(b).unwrap());s[((q*(s.len()-1)as f32)as usize).min(s.len()-1)]};
    eprintln!("#155 orogen-render fraction elucidation (64²)");
    for &seed in &[1988u64,2026,42]{
        let mut state=init_c1_state_phase_2_r7(grid,seed,&Phase2InitParams::default());
        let mut kin=PlateKinematics::preset_phase_1_1(state.num_plates);
        let config=C1TimeLoopConfig{rigid_continental_crust:true,n_steps:300,dx:1.0/64.0,dy:1.0/64.0,iso_config:iso.clone(),drainage_max_distance:30};
        run_with_closures(&mut state,&mut kin,&config,&C1Closures::default(),|_,_|{});
        let mut raw=vec![0.0f32;grid*grid]; let mut raw_cont=Vec::new();
        for j in 0..grid{for i in 0..grid{let k=j*grid+i;let b=if state.cratonic_mask.data()[k]{buoy_c}else{buoy};let h=(state.s.get(i,j) as f32)*b;raw[k]=h;if matches!(state.plate_type.get(i,j),PlateType::Continental){raw_cont.push(h);}}}
        let h_min_g=raw.iter().cloned().fold(f32::INFINITY,f32::min);
        let h_min_c=raw_cont.iter().cloned().fold(f32::INFINITY,f32::min);
        let land_cap=raw_cont.iter().cloned().fold(f32::NEG_INFINITY,f32::max);
        let h_cap=pct(&raw,0.92); let h_range=(h_cap-h_min_g).max(1e-10); let h_sea=h_min_g+0.4*h_range;
        let frac_g=(land_cap-h_sea)/(land_cap-h_min_g);
        let frac_c=(land_cap-h_sea)/(land_cap-h_min_c).max(1e-10);
        // fixed-Airy alternative: orogen elevation = (land_cap - h_sea) directly (raw above sea),
        // vs the range-normalized fraction.
        eprintln!("  seed {seed}: land_cap={land_cap:.3} h_sea={h_sea:.3} h_min global={h_min_g:.3} continental={h_min_c:.3}",);
        eprintln!("    fraction (land_cap-h_sea)/(land_cap-h_min): GLOBAL={frac_g:.3}  CONTINENTAL={frac_c:.3}  (global<continental ⇒ ocean h_min drives the <1)");
    }
}

/// #165 c1_climate ACCEPTANCE — triple validation on the real C1 product @2048²:
/// (1) CONSERVATION (the model's law: evap_in = precip + exit, exact); (2)
/// MAGNITUDE (windward/leeward precip ratio ~3-10×, real ranges); (3) STRUCTURE
/// (visual: west-facing slopes wet, east dry at 45° westerlies). Default 45°.
#[test]
#[ignore]
fn probe_climate_acceptance() {
    use ymir_core::climate::c1_climate;
    use ymir_core::climate::precipitation::{compute_precipitation_with_budget, PrecipParams, SEA_LEVEL_NORM};
    use ymir_core::tectonics_c1::production_upscale::{c1_altitude_norm_to_metres, c1_km_per_cell};
    let dir = output_dir().join("climate"); std::fs::create_dir_all(&dir).unwrap();
    let iso = IsostasyConfig::c1_default(); let ss = SteinSteinParams::default(); let grid = 64usize;
    let lat = 45.0f32; let pp = PrecipParams::default();
    eprintln!("#165 c1_climate acceptance — 2048², latitude {lat}° (westerlies, W→E)");
    for &seed in &[1988u64, 2026, 42] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
        run_with_closures(&mut state, &mut kin, &config, &C1Closures::default(), |_, _| {});
        let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &FbmUpscaleConfig::c1_hd_production(2048));
        let h = &up.heightmap; let (w, ht) = (h.width, h.height);
        let clim = c1_climate(h, &ss, lat, &pp);

        // (1) CONSERVATION
        let temp = &clim.temperature;
        let (_p, evap_in, exit_out, oro_sum) = compute_precipitation_with_budget(h, temp, lat, c1_km_per_cell(w), |n| c1_altitude_norm_to_metres(n, &ss), &pp);
        let psum: f64 = clim.precipitation.data.iter().map(|&v| v as f64).sum();
        // OROGRAPHIC conservation (the convective baseline is a separate source).
        let residual = (evap_in - (oro_sum + exit_out)).abs();
        // (2) MAGNITUDE — windward (W-facing ascending) vs leeward (descending) land precip.
        let dx_m = c1_km_per_cell(w) * 1000.0;
        let (mut wp, mut wn, mut lp, mut ln) = (0.0f64, 0u64, 0.0f64, 0u64);
        for j in 0..ht { for i in 1..w-1 {
            let k = j*w+i; if h.data[k] <= SEA_LEVEL_NORM { continue; }
            // along-wind (eastward) slope: cell to its west (i-1).
            let altw = c1_altitude_norm_to_metres(h.data[k-1], &ss).max(0.0);
            let alt = c1_altitude_norm_to_metres(h.data[k], &ss).max(0.0);
            let asc = (alt - altw)/dx_m; // >0 = ascending eastward = windward (W-facing)
            let pr = clim.precipitation.data[k] as f64;
            if asc > 1e-4 { wp += pr; wn += 1; } else if asc < -1e-4 { lp += pr; ln += 1; }
        }}
        let (wm, lm) = (wp/wn.max(1) as f64, lp/ln.max(1) as f64);
        // T stats
        let tmin = temp.data.iter().cloned().fold(f32::INFINITY, f32::min);
        let tmax = temp.data.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        eprintln!("  seed {seed}: CONSERVE evap={evap_in:.0} = oro_precip {oro_sum:.0} + exit {exit_out:.0} (resid {residual:.2e}) | total precip(incl convective) {psum:.0}");
        eprintln!("    MAGNITUDE windward mean {wm:.3} / leeward mean {lm:.3} = ratio {:.1}× (want ~3-10×) | T [{tmin:.0}, {tmax:.0}]°C", wm/lm.max(1e-9));

        // (3) STRUCTURE — precip (blue) + temperature (cyan→red) maps, origin-bottom.
        // Robust scale: clamp at the land-precip P90 + sqrt, so the windward/leeward
        // GRADIENT is visible (a few steep-coast outliers must not compress it — the
        // measure-the-structure-not-the-compressed-output lesson).
        let mut landp: Vec<f32> = (0..w*ht).filter(|&k| h.data[k] > SEA_LEVEL_NORM).map(|k| clim.precipitation.data[k]).collect();
        landp.sort_by(|a,b| a.partial_cmp(b).unwrap());
        let p90 = if landp.is_empty() { 1e-6 } else { landp[(landp.len()*90/100).min(landp.len()-1)].max(1e-6) };
        let mut pb = image::ImageBuffer::new(w as u32, ht as u32);
        let mut tb = image::ImageBuffer::new(w as u32, ht as u32);
        for j in 0..ht { for i in 0..w { let k=j*w+i;
            let sea = h.data[k] <= SEA_LEVEL_NORM;
            // precip map: ocean dark, land tan(dry)→blue(wet), sqrt of precip/P90.
            let pc = if sea { [20,40,80] } else { let t=(clim.precipitation.data[k]/p90).clamp(0.0,1.0).sqrt(); [ (210.0-180.0*t) as u8, (200.0-120.0*t) as u8, (170.0+85.0*t).min(255.0) as u8 ] };
            pb.put_pixel(i as u32,(ht-1-j) as u32, image::Rgb(pc));
            // temp map: blue(cold)→red(hot), [-20,30].
            let tt = ((temp.data[k]+20.0)/50.0).clamp(0.0,1.0); let tc=[ (40.0+200.0*tt) as u8, 60u8, (240.0-200.0*tt) as u8 ];
            tb.put_pixel(i as u32,(ht-1-j) as u32, image::Rgb(if sea {[30,50,90]} else {tc}));
        }}
        pb.save(dir.join(format!("seed{seed:05}_precip.png"))).unwrap();
        tb.save(dir.join(format!("seed{seed:05}_temp.png"))).unwrap();
    }
    eprintln!("  out = {}", dir.display());
}

/// #165 precip STRUCTURE check (cheap, 512²) — is the W→E rain-shadow gradient
/// real (just visually subtle vs the uniform convective baseline) or flat?
/// Per row, split the land run at its midpoint and compare W-half vs E-half mean
/// precip (at 45° westerlies, W should be wetter). + land precip percentiles.
#[test]
#[ignore]
fn probe_precip_structure() {
    use ymir_core::climate::c1_climate;
    use ymir_core::climate::precipitation::{PrecipParams, SEA_LEVEL_NORM};
    let iso = IsostasyConfig::c1_default(); let ss = SteinSteinParams::default(); let grid = 64usize;
    let pp = PrecipParams::default();
    eprintln!("#165 precip structure (512², 45° westerlies W→E)");
    for &seed in &[1988u64, 2026, 42] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
        run_with_closures(&mut state, &mut kin, &config, &C1Closures::default(), |_, _| {});
        let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &FbmUpscaleConfig::c1_hd_production(512));
        let h = &up.heightmap; let (w, ht) = (h.width, h.height);
        let clim = c1_climate(h, &ss, 45.0, &pp);
        // per-row W-half vs E-half land precip.
        let (mut wsum, mut wn, mut esum, mut en) = (0.0f64, 0u64, 0.0f64, 0u64);
        for j in 0..ht {
            let land: Vec<usize> = (0..w).filter(|&i| h.data[j*w+i] > SEA_LEVEL_NORM).collect();
            if land.len() < 4 { continue; }
            let mid = land[land.len()/2];
            for &i in &land {
                let p = clim.precipitation.data[j*w+i] as f64;
                if i < mid { wsum += p; wn += 1; } else { esum += p; en += 1; }
            }
        }
        let (wm, em) = (wsum/wn.max(1) as f64, esum/en.max(1) as f64);
        // land precip percentiles.
        let mut lp: Vec<f32> = (0..w*ht).filter(|&k| h.data[k] > SEA_LEVEL_NORM).map(|k| clim.precipitation.data[k]).collect();
        lp.sort_by(|a,b| a.partial_cmp(b).unwrap());
        let pct = |q: f32| lp[((q*(lp.len()-1) as f32) as usize).min(lp.len()-1)];
        eprintln!("  seed {seed}: W-half mean {wm:.4} / E-half mean {em:.4} = {:.2}× (>1 = rain shadow eastward) | land precip p10={:.4} p50={:.4} p90={:.4} max={:.4}",
            wm/em.max(1e-9), pct(0.1), pct(0.5), pct(0.9), pct(1.0));
    }
}

/// #165 climate VISUAL validation — precip in LOG + LINEAR + windward-marked, and
/// temperature; 2 seeds @2048², 45° westerlies, origin-bottom. Log reveals the
/// gradient a linear colormap drowns (skewed dist). Reports p50/max + the
/// floor-fraction (land cells at the convective minimum — suspect if a temperate
/// continent is mostly desert-floor).
#[test]
#[ignore]
fn probe_climate_maps() {
    use ymir_core::climate::c1_climate;
    use ymir_core::climate::precipitation::{PrecipParams, SEA_LEVEL_NORM};
    use ymir_core::tectonics_c1::production_upscale::{c1_altitude_norm_to_metres, c1_km_per_cell};
    let dir = output_dir().join("climate_maps"); std::fs::create_dir_all(&dir).unwrap();
    let iso = IsostasyConfig::c1_default(); let ss = SteinSteinParams::default(); let grid = 64usize;
    let pp = PrecipParams::default();
    // hillshade (origin-bottom) for the relief panel.
    let shade=|h:&GridF32,i:usize,j:usize|->u8{let z=55.0_f64;let(lx,ly,lz)={let(a,b,c)=(1.0_f64,1.0,2.0);let nn=(a*a+b*b+c*c).sqrt();(a/nn,b/nn,c/nn)};let dx=(h.get(i as i32+1,j as i32)-h.get(i as i32-1,j as i32))as f64*z;let dy=(h.get(i as i32,j as i32+1)-h.get(i as i32,j as i32-1))as f64*z;let n=(dx*dx+dy*dy+1.0).sqrt();(((-dx*lx-dy*ly+lz)/n).clamp(0.0,1.0)*255.0)as u8};
    for &seed in &[42u64, 99, 1337, 4138, 1988, 2026] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
        run_with_closures(&mut state, &mut kin, &config, &C1Closures::default(), |_, _| {});
        let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &FbmUpscaleConfig::c1_hd_production(2048));
        let h = &up.heightmap; let (w, ht) = (h.width, h.height);
        let clim = c1_climate(h, &ss, 45.0, &pp);
        let pr = &clim.precipitation.data; let temp = &clim.temperature.data;
        let dx_m = c1_km_per_cell(w) * 1000.0;
        // land precip stats.
        let mut lp: Vec<f32> = (0..w*ht).filter(|&k| h.data[k] > SEA_LEVEL_NORM).map(|k| pr[k]).collect();
        lp.sort_by(|a,b| a.partial_cmp(b).unwrap());
        let n = lp.len().max(1);
        let p50 = lp[n/2]; let pmax = lp[n-1].max(1e-6); let pmin = lp[0].max(1e-6);
        let floor = lp.iter().filter(|&&v| v < 0.015).count();
        eprintln!("  seed {seed}: land {n} | precip p50={p50:.4} max={pmax:.4} | floor(<0.015) {:.0}% | log range [{:.2},{:.2}]",
            100.0*floor as f64/n as f64, pmin.log10(), pmax.log10());
        // FIXED-reference colormap (NOT field min/max — those are dominated by the
        // steep-coast windward outliers and squash the uniform temperate base).
        // Absolute anchors: desert P_DRY → tan, temperate base ~0.036 → moderate,
        // heavy rain P_WET → blue. Honest across seeds + vs the floor.
        const P_DRY: f32 = 0.005; const P_WET: f32 = 0.5;
        let lnmin = P_DRY.ln(); let lnmax = P_WET.ln();
        let put = |buf:&mut image::RgbImage, i:usize, j:usize, c:[u8;3]| buf.put_pixel(i as u32, (ht-1-j) as u32, image::Rgb(c));
        let (mut blog, mut blin, mut bslp, mut btmp, mut brel) = (image::RgbImage::new(w as u32,ht as u32), image::RgbImage::new(w as u32,ht as u32), image::RgbImage::new(w as u32,ht as u32), image::RgbImage::new(w as u32,ht as u32), image::RgbImage::new(w as u32,ht as u32));
        for j in 0..ht { for i in 0..w { let k=j*w+i;
            let sea = h.data[k] <= SEA_LEVEL_NORM;
            let p = pr[k];
            // LOG
            let tl = if sea {0.0} else { ((p.max(1e-6).ln()-lnmin)/(lnmax-lnmin)).clamp(0.0,1.0) };
            put(&mut blog,i,j, if sea {[20,40,80]} else {[ (215.0-185.0*tl) as u8,(205.0-125.0*tl) as u8,(165.0+90.0*tl) as u8 ]});
            // LINEAR (fixed P_WET clamp, not /max — so the base isn't squashed)
            let tlin = if sea {0.0} else { (p/P_WET).clamp(0.0,1.0) };
            put(&mut blin,i,j, if sea {[20,40,80]} else {[ (215.0-185.0*tlin) as u8,(205.0-125.0*tlin) as u8,(165.0+90.0*tlin) as u8 ]});
            // SLOPES: windward (ascending eastward) tinted by precip-blue, leeward red.
            let cslp = if sea {[20,40,80]} else {
                let altw = if i>0 {c1_altitude_norm_to_metres(h.data[k-1],&ss).max(0.0)} else {0.0};
                let asc = (c1_altitude_norm_to_metres(h.data[k],&ss).max(0.0)-altw)/dx_m;
                let b = (tl*255.0) as u8;
                if asc > 1e-4 { [ (60.0-40.0*tl) as u8, (120.0+100.0*tl).min(255.0) as u8, b.max(120) ] } // windward, precip-blue
                else if asc < -1e-4 { [ 180, (80.0+60.0*tl) as u8, 60 ] } // leeward orange
                else { [140,140,140] }
            };
            put(&mut bslp,i,j,cslp);
            // TEMP
            let tt=((temp[k]+20.0)/50.0).clamp(0.0,1.0);
            put(&mut btmp,i,j, if sea {[30,50,90]} else {[ (40.0+200.0*tt) as u8,60,(240.0-200.0*tt) as u8 ]});
            // RELIEF (hillshade)
            let g = if sea {70} else { shade(h,i,j) };
            put(&mut brel,i,j, if sea {[40,70,120]} else {[g,g,g]});
        }}
        brel.save(dir.join(format!("seed{seed:05}_relief.png"))).unwrap();
        blog.save(dir.join(format!("seed{seed:05}_precip_LOG.png"))).unwrap();
        blin.save(dir.join(format!("seed{seed:05}_precip_LINEAR.png"))).unwrap();
        bslp.save(dir.join(format!("seed{seed:05}_precip_SLOPES.png"))).unwrap();
        btmp.save(dir.join(format!("seed{seed:05}_temp.png"))).unwrap();
    }
    eprintln!("  out = {}", dir.display());
}

/// #165 climate CONTROL survey — find a NORMAL-relief continent (passive/low
/// WINDWARD (west) margin, interior open to the westerly ocean) to disambiguate
/// dry-interior cause: climatic (no frontal base) vs local (cordillera
/// encirclement). Identify by RELIEF (west-margin height), not climate. 1024².
#[test]
#[ignore]
fn probe_climate_control_survey() {
    use ymir_core::climate::c1_climate;
    use ymir_core::climate::precipitation::{PrecipParams, SEA_LEVEL_NORM};
    use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
    let dir = output_dir().join("climate_control"); std::fs::create_dir_all(&dir).unwrap();
    let iso = IsostasyConfig::c1_default(); let ss = SteinSteinParams::default(); let grid = 64usize;
    let pp = PrecipParams::default();
    let shade=|h:&GridF32,i:usize,j:usize|->u8{let z=55.0_f64;let(lx,ly,lz)={let(a,b,c)=(1.0_f64,1.0,2.0);let nn=(a*a+b*b+c*c).sqrt();(a/nn,b/nn,c/nn)};let dx=(h.get(i as i32+1,j as i32)-h.get(i as i32-1,j as i32))as f64*z;let dy=(h.get(i as i32,j as i32+1)-h.get(i as i32,j as i32-1))as f64*z;let n=(dx*dx+dy*dy+1.0).sqrt();(((-dx*lx-dy*ly+lz)/n).clamp(0.0,1.0)*255.0)as u8};
    eprintln!("#165 climate control survey (1024², 45° westerlies; windward=WEST). margin = mean altitude(m) of the 6 westmost/eastmost land cells per row.");
    for &seed in &[1337u64, 4138, 42, 99] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
        run_with_closures(&mut state, &mut kin, &config, &C1Closures::default(), |_, _| {});
        let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &FbmUpscaleConfig::c1_hd_production(1024));
        let h = &up.heightmap; let (w, ht) = (h.width, h.height);
        let alt = |k:usize| c1_altitude_norm_to_metres(h.data[k], &ss).max(0.0);
        // west/east margin relief: per row, the 6 westmost & 6 eastmost land cells.
        let (mut wsum,mut wn,mut esum,mut en)=(0.0f64,0u64,0.0f64,0u64);
        for j in 0..ht {
            let land:Vec<usize>=(0..w).filter(|&i| h.data[j*w+i]>SEA_LEVEL_NORM).collect();
            if land.len()<12 {continue;}
            for &i in land.iter().take(6) { wsum+=alt(j*w+i) as f64; wn+=1; }
            for &i in land.iter().rev().take(6) { esum+=alt(j*w+i) as f64; en+=1; }
        }
        let (wm,em)=(wsum/wn.max(1) as f64, esum/en.max(1) as f64);
        let clim = c1_climate(h,&ss,45.0,&pp);
        let mut lp:Vec<f32>=(0..w*ht).filter(|&k|h.data[k]>SEA_LEVEL_NORM).map(|k|clim.precipitation.data[k]).collect();
        lp.sort_by(|a,b|a.partial_cmp(b).unwrap()); let nl=lp.len().max(1);
        let floor=100.0*lp.iter().filter(|&&v|v<0.015).count() as f64/nl as f64;
        eprintln!("  seed {seed}: WEST margin {wm:.0} m / EAST margin {em:.0} m | floor {floor:.0}% | (low west = passive windward = control candidate)");
        // renders: hillshade + precip LOG.
        let mut lpmin=f32::MAX; let mut lpmax=0.0f32; for &v in &lp {lpmin=lpmin.min(v.max(1e-6));lpmax=lpmax.max(v);}
        let (lnmin,lnmax)=(lpmin.ln(), lpmax.max(1e-6).ln());
        let mut bh=image::RgbImage::new(w as u32,ht as u32); let mut bp=image::RgbImage::new(w as u32,ht as u32);
        for j in 0..ht{for i in 0..w{let k=j*w+i; let sea=h.data[k]<=SEA_LEVEL_NORM;
            let g=if sea{70}else{shade(h,i,j)};
            bh.put_pixel(i as u32,(ht-1-j)as u32,image::Rgb(if sea{[40,70,120]}else{[g,g,g]}));
            let tl=if sea{0.0}else{((clim.precipitation.data[k].max(1e-6).ln()-lnmin)/(lnmax-lnmin)).clamp(0.0,1.0)};
            bp.put_pixel(i as u32,(ht-1-j)as u32,image::Rgb(if sea{[20,40,80]}else{[(215.0-185.0*tl)as u8,(205.0-125.0*tl)as u8,(165.0+90.0*tl)as u8]}));
        }}
        bh.save(dir.join(format!("seed{seed:05}_relief.png"))).unwrap();
        bp.save(dir.join(format!("seed{seed:05}_precip_log.png"))).unwrap();
    }
    eprintln!("  out = {}", dir.display());
}

/// #165 orographic DEPLETION diagnostic — settle depletion vs low-interior-relief
/// by normalizing the orographic surplus by the along-wind ascent. precip_oro =
/// k_oro·M·ascent → oro/ascent = k_oro·M ∝ the carried moisture flux. Exposes the
/// orographic field ALONE (k_frontal=0, no uniform base masking it) and bins
/// oro/ascent by distance-from-windward-coast (reset per landmass, W→E at 45°).
/// COLLAPSE = depletion (air wrung out); CONSTANT = low relief (not depletion).
#[test]
#[ignore]
fn probe_oro_depletion() {
    use ymir_core::climate::c1_climate;
    use ymir_core::climate::precipitation::{PrecipParams, SEA_LEVEL_NORM};
    use ymir_core::tectonics_c1::production_upscale::{c1_altitude_norm_to_metres, c1_km_per_cell};
    let dir = output_dir().join("climate_maps"); std::fs::create_dir_all(&dir).unwrap();
    let iso = IsostasyConfig::c1_default(); let ss = SteinSteinParams::default(); let grid = 64usize;
    // orographic ONLY: frontal base off (the uniform base masks the depletion profile).
    let pp = PrecipParams { k_frontal: 0.0, ..PrecipParams::default() };
    let bins_km = [0.0f32, 10.0, 30.0, 60.0, 120.0, 250.0, f32::INFINITY];
    eprintln!("#165 orographic depletion — oro/ascent (∝ moisture flux M) by distance-from-windward-coast (km), 1024², 45° W→E");
    for &seed in &[1988u64, 2026, 42, 99] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
        run_with_closures(&mut state, &mut kin, &config, &C1Closures::default(), |_, _| {});
        let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &FbmUpscaleConfig::c1_hd_production(1024));
        let h = &up.heightmap; let (w, ht) = (h.width, h.height);
        let clim = c1_climate(h, &ss, 45.0, &pp); // orographic-only
        let oro = &clim.precipitation.data;
        let km_cell = c1_km_per_cell(w); let dx_m = km_cell * 1000.0;
        let alt = |k:usize| c1_altitude_norm_to_metres(h.data[k], &ss).max(0.0);
        let mut sum = [0.0f64; 6]; let mut cnt = [0u64; 6];
        for j in 0..ht {
            let mut dist_cells = 0u32; let mut in_land = false;
            for i in 0..w { let k=j*w+i;
                let land = h.data[k] > SEA_LEVEL_NORM;
                if land {
                    if !in_land { dist_cells = 0; in_land = true; } else { dist_cells += 1; }
                    if i>0 {
                        let asc = (alt(k) - alt(k-1)) / dx_m; // along +x (wind) ascent, m/m
                        if asc > 1e-4 {
                            let r = oro[k] as f64 / asc as f64; // ∝ M
                            let dkm = dist_cells as f32 * km_cell;
                            for b in 0..6 { if dkm>=bins_km[b] && dkm<bins_km[b+1] { sum[b]+=r; cnt[b]+=1; break; } }
                        }
                    }
                } else { in_land = false; }
            }
        }
        let means: Vec<f64> = (0..6).map(|b| if cnt[b]>0 {sum[b]/cnt[b] as f64} else {f64::NAN}).collect();
        let m0 = means[0];
        // e-folding: first bin whose mean < 0.37*m0.
        let efold = (0..6).find(|&b| means[b].is_finite() && means[b] < 0.37*m0).map(|b| (bins_km[b]+bins_km[b+1].min(300.0))*0.5);
        eprintln!("  seed {seed}: oro/ascent by dist[0-10|10-30|30-60|60-120|120-250|250+]km = {:?}",
            means.iter().map(|v| (v*1000.0).round()/1000.0).collect::<Vec<_>>());
        eprintln!("    -> {} | e-fold(<37% of coastal) ~{:?} km",
            if means.iter().take(4).filter(|v|v.is_finite()).count()>=2 && means[3].is_finite() && means[3] < 0.6*m0 {"COLLAPSE = DEPLETION"} else {"~CONSTANT = low-relief (not depletion)"}, efold);
        if seed == 1988 {
            // render orographic-only LOG.
            let mut lp:Vec<f32>=(0..w*ht).filter(|&k|h.data[k]>SEA_LEVEL_NORM).map(|k|oro[k]).collect();
            lp.sort_by(|a,b|a.partial_cmp(b).unwrap()); let pmax=lp[lp.len()-1].max(1e-6);
            let (lnmin,lnmax)=(1e-4f32.ln(), pmax.ln());
            let mut b=image::RgbImage::new(w as u32,ht as u32);
            for j in 0..ht{for i in 0..w{let k=j*w+i;let sea=h.data[k]<=SEA_LEVEL_NORM;
                let tl=if sea{0.0}else{((oro[k].max(1e-4).ln()-lnmin)/(lnmax-lnmin)).clamp(0.0,1.0)};
                b.put_pixel(i as u32,(ht-1-j)as u32,image::Rgb(if sea{[20,40,80]}else{[(215.0-185.0*tl)as u8,(205.0-125.0*tl)as u8,(165.0+90.0*tl)as u8]}));
            }}
            b.save(dir.join("seed01988_precip_OROGRAPHIC_ONLY.png")).unwrap();
        }
    }
    eprintln!("  out = {}", dir.display());
}

/// #165 BIOMES + mm/yr legend — the judgement instruments. Per seed (2048²,
/// 45°): c1_climate → c1_biomes; renders the categorical BIOME map, precip
/// BANDED by mm/yr (desert/steppe/temperate/oceanic/wet — legend baked into the
/// colormap), temp banded by °C; + prints the biome HISTOGRAM (% per biome —
/// directly judgeable: steppe vs desert interior, etc.). Origin-bottom.
#[test]
#[ignore]
fn probe_climate_biomes() {
    use ymir_core::climate::{c1_climate, c1_biomes};
    use ymir_core::climate::precipitation::{PrecipParams, precip_mm_per_year, SEA_LEVEL_NORM};
    use ymir_core::climate::biomes::Biome;
    let dir = output_dir().join("climate_maps"); std::fs::create_dir_all(&dir).unwrap();
    let iso = IsostasyConfig::c1_default(); let ss = SteinSteinParams::default(); let grid = 64usize;
    let pp = PrecipParams::default();
    // mm/yr bands (legend baked in): desert / steppe / temperate-dry / oceanic / wet.
    let precip_band = |mm: f32| -> [u8;3] {
        if mm < 250.0 {[225,200,140]} else if mm < 500.0 {[200,195,110]} else if mm < 800.0 {[150,180,90]}
        else if mm < 1500.0 {[80,150,200]} else {[30,90,200]}
    };
    let temp_band = |t: f32| -> [u8;3] {
        if t < -10.0 {[230,235,245]} else if t < 0.0 {[150,180,230]} else if t < 10.0 {[120,190,120]}
        else if t < 20.0 {[220,210,120]} else {[220,110,70]}
    };
    eprintln!("#165 biomes + mm/yr (2048², 45° westerlies). precip bands: desert<250 steppe250-500 tempDry500-800 oceanic800-1500 wet>1500 mm/yr");
    for &seed in &[42u64, 99, 1337, 4138, 1988, 2026] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
        run_with_closures(&mut state, &mut kin, &config, &C1Closures::default(), |_, _| {});
        let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &FbmUpscaleConfig::c1_hd_production(2048));
        let h = &up.heightmap; let (w, ht) = (h.width, h.height);
        let clim = c1_climate(h, &ss, 45.0, &pp);
        let biomes = c1_biomes(h, &clim);
        // biome histogram over LAND.
        let order = [Biome::Tundra, Biome::BorealForest, Biome::TemperateGrassland, Biome::TemperateForest, Biome::TemperateRainforest, Biome::Desert, Biome::Savanna, Biome::TropicalSeasonalForest, Biome::TropicalRainforest];
        let mut land = 0u64; let mut cnt = std::collections::HashMap::new();
        for &b in &biomes { if b != Biome::Ocean { land += 1; *cnt.entry(b).or_insert(0u64) += 1; } }
        let hist: Vec<String> = order.iter().filter_map(|b| { let c=*cnt.get(b).unwrap_or(&0); if c>0 {Some(format!("{} {:.0}%", b.name(), 100.0*c as f64/land.max(1) as f64))} else {None} }).collect();
        eprintln!("  seed {seed}: {}", hist.join(" | "));
        // renders.
        let put=|buf:&mut image::RgbImage,i:usize,j:usize,c:[u8;3]|buf.put_pixel(i as u32,(ht-1-j)as u32,image::Rgb(c));
        let (mut bb, mut bpm, mut bt) = (image::RgbImage::new(w as u32,ht as u32), image::RgbImage::new(w as u32,ht as u32), image::RgbImage::new(w as u32,ht as u32));
        for j in 0..ht{for i in 0..w{let k=j*w+i; let sea=h.data[k]<=SEA_LEVEL_NORM;
            put(&mut bb,i,j, biomes[k].color());
            put(&mut bpm,i,j, if sea {[30,50,90]} else {precip_band(precip_mm_per_year(clim.precipitation.data[k]))});
            put(&mut bt,i,j, if sea {[30,50,90]} else {temp_band(clim.temperature.data[k])});
        }}
        bb.save(dir.join(format!("seed{seed:05}_BIOMES.png"))).unwrap();
        bpm.save(dir.join(format!("seed{seed:05}_precip_MMBANDS.png"))).unwrap();
        bt.save(dir.join(format!("seed{seed:05}_temp_CBANDS.png"))).unwrap();
    }
    eprintln!("  out = {}", dir.display());
}

/// #165 A/B VERDICT GRID — the verdict lives in the CROSS windward-coast-height
/// (MEASURED in metres, not eyeballed) × interior-humidity (biome), NOT in
/// biomes or relief alone: a dry interior is CORRECT behind a wall (legitimate
/// Patagonia shadow) and a PROBLEM behind a low coast (rain should penetrate).
/// At 45° westerlies the air enters from the WEST, so "windward coast" = the
/// west margin. Per seed: measure the windward coastal-band altitude (max + p95,
/// metres via the vertical contract), the interior dominant biome + precip
/// (mm/yr), classify into the 4 cases, and render relief (windward height baked
/// into the filename + a red tint on the windward band) / biomes / a W→E transect
/// (altitude + precip + biome strip). MEASURES; does NOT decide A vs B.
#[test]
#[ignore]
fn probe_climate_verdict_grid() {
    use ymir_core::climate::{c1_climate, c1_biomes};
    use ymir_core::climate::precipitation::{PrecipParams, precip_mm_per_year, SEA_LEVEL_NORM};
    use ymir_core::climate::biomes::Biome;
    use ymir_core::tectonics_c1::production_upscale::{c1_altitude_norm_to_metres, c1_km_per_cell};
    let dir = output_dir().join("climate_verdict"); std::fs::create_dir_all(&dir).unwrap();
    let iso = IsostasyConfig::c1_default(); let ss = SteinSteinParams::default(); let grid = 64usize;
    let pp = PrecipParams::default();
    let lat = 45.0f32;
    let alt_m = |n: f32| c1_altitude_norm_to_metres(n, &ss).max(0.0);
    // interior is SEC when the rain-shadow biomes (desert + steppe) dominate it.
    let dry = |b: Biome| matches!(b, Biome::Desert | Biome::TemperateGrassland);
    eprintln!("#165 A/B verdict grid — 2048², 45° westerlies (wind W→E); windward = WEST margin.");
    eprintln!("  coast: HAUTE>3000m / BASSE<2000m (2-3km INTERM). interior SEC if desert+steppe >50%.");
    for &seed in &[42u64, 99, 1337, 4138, 1988, 2026] {
        let mut state = init_c1_state_phase_2_r7(grid, seed, &Phase2InitParams::default());
        let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
        let config = C1TimeLoopConfig { rigid_continental_crust: true, n_steps: 300, dx: 1.0/64.0, dy: 1.0/64.0, iso_config: iso.clone(), drainage_max_distance: 30 };
        run_with_closures(&mut state, &mut kin, &config, &C1Closures::default(), |_, _| {});
        let up = upscale_from_c1(&state, &iso, &ss, &WorldSeed::new(seed), &FbmUpscaleConfig::c1_hd_production(2048));
        let h = &up.heightmap; let (w, ht) = (h.width, h.height);
        let km_cell = c1_km_per_cell(w);
        let clim = c1_climate(h, &ss, lat, &pp);
        let biomes = c1_biomes(h, &clim);
        let w_band = (80.0 / km_cell) as usize;   // windward coastal band: first 80 km of land
        let i_start = (150.0 / km_cell) as usize;  // interior begins 150 km inland of the W coast
        let mut wind_alts: Vec<f32> = Vec::new();
        let mut best_row = (ht/2, -1.0f32);        // row with the tallest windward peak → transect
        let mut interior_cnt: std::collections::HashMap<Biome, u64> = std::collections::HashMap::new();
        let mut interior_land = 0u64; let mut interior_dry = 0u64;
        let mut interior_precip: Vec<f32> = Vec::new();
        for j in 0..ht {
            let mut i0 = None;
            for i in 0..w { if h.data[j*w+i] > SEA_LEVEL_NORM { i0 = Some(i); break; } }
            let Some(i0) = i0 else { continue; };
            let mut row_peak = 0.0f32;
            for i in i0..(i0+w_band).min(w) {
                let k = j*w+i; if h.data[k] > SEA_LEVEL_NORM { let a = alt_m(h.data[k]); wind_alts.push(a); row_peak = row_peak.max(a); }
            }
            if row_peak > best_row.1 { best_row = (j, row_peak); }
            for i in (i0+i_start)..w {
                let k = j*w+i; if h.data[k] > SEA_LEVEL_NORM {
                    interior_land += 1;
                    let b = biomes[k];
                    *interior_cnt.entry(b).or_insert(0) += 1;
                    if dry(b) { interior_dry += 1; }
                    interior_precip.push(precip_mm_per_year(clim.precipitation.data[k]));
                }
            }
        }
        wind_alts.sort_by(|a,b| a.partial_cmp(b).unwrap());
        interior_precip.sort_by(|a,b| a.partial_cmp(b).unwrap());
        let pct = |v:&Vec<f32>, p:f64| if v.is_empty() {0.0} else { v[(((v.len()-1) as f64)*p) as usize] };
        let wind_max = wind_alts.last().copied().unwrap_or(0.0);
        let wind_p95 = pct(&wind_alts, 0.95);
        let int_p50 = pct(&interior_precip, 0.5);
        let int_p90 = pct(&interior_precip, 0.9);
        let dry_frac = if interior_land>0 {100.0*interior_dry as f64/interior_land as f64} else {0.0};
        let dom = interior_cnt.iter().max_by_key(|(_,c)| **c).map(|(b,_)| b.name()).unwrap_or("-");
        let coast = if wind_max > 3000.0 {"HAUTE"} else if wind_max < 2000.0 {"BASSE"} else {"INTERM"};
        let hum = if dry_frac > 50.0 {"SEC"} else {"HUMIDE"};
        let case = match (coast, hum) {
            ("HAUTE","SEC")    => "ombre légitime (Patagonie) [correct]",
            ("BASSE","HUMIDE") => "océanique [correct]",
            ("BASSE","SEC")    => "PROBLÈME (pluie ne pénètre pas)",
            ("HAUTE","HUMIDE") => "improbable (haut + humide)",
            _                  => "intermédiaire (coast 2-3km)",
        };
        eprintln!("  seed {seed}: windward max {wind_max:.0}m p95 {wind_p95:.0}m [{coast}] | interior dom={dom} dry={dry_frac:.0}% precip p50 {int_p50:.0} p90 {int_p90:.0} mm [{hum}] => {case}");
        // renders.
        let put=|buf:&mut image::RgbImage,i:usize,j:usize,c:[u8;3]|buf.put_pixel(i as u32,(ht-1-j)as u32,image::Rgb(c));
        let mut relief = image::RgbImage::new(w as u32, ht as u32);
        let mut bb = image::RgbImage::new(w as u32, ht as u32);
        for j in 0..ht { for i in 0..w { let k=j*w+i;
            if h.data[k] <= SEA_LEVEL_NORM { put(&mut relief,i,j,[30,50,90]); }
            else { let g=((alt_m(h.data[k])/6000.0)*255.0).min(255.0) as u8; put(&mut relief,i,j,[g,g,g]); }
            put(&mut bb,i,j, biomes[k].color());
        }}
        // red tint on the measured windward band (first 80 km of land per row).
        for j in 0..ht {
            let mut i0=None; for i in 0..w { if h.data[j*w+i] > SEA_LEVEL_NORM { i0=Some(i); break; } }
            if let Some(i0)=i0 { for i in i0..(i0+w_band).min(w) { let k=j*w+i; if h.data[k]>SEA_LEVEL_NORM {
                let y=(ht-1-j) as u32; let p=relief.get_pixel(i as u32,y).0;
                relief.put_pixel(i as u32,y, image::Rgb([(p[0] as u16+90).min(255) as u8, p[1]/2, p[2]/2])); } } }
        }
        relief.save(dir.join(format!("seed{seed:05}_relief_wind{:.0}m_{coast}.png", wind_max))).unwrap();
        bb.save(dir.join(format!("seed{seed:05}_biomes.png"))).unwrap();
        // transect on the tallest-windward row: altitude (grey) + precip (cyan) + biome strip.
        let jr = best_row.0; let th = 400usize;
        let mut tr = image::RgbImage::new(w as u32, th as u32);
        for p in tr.pixels_mut() { *p = image::Rgb([20,20,28]); }
        let max_mm = (0..w).map(|i| precip_mm_per_year(clim.precipitation.data[jr*w+i])).fold(1.0f32, f32::max);
        for i in 0..w {
            let k=jr*w+i;
            let a = if h.data[k]>SEA_LEVEL_NORM {alt_m(h.data[k])} else {0.0};
            let ay = (((a/6000.0).min(1.0))*(th as f32 - 40.0)) as usize;
            let mm = precip_mm_per_year(clim.precipitation.data[k]);
            let my = (((mm/max_mm).min(1.0))*(th as f32 - 40.0)) as usize;
            if ay < th { tr.put_pixel(i as u32, (th-1-ay) as u32, image::Rgb([200,200,200])); }
            if my < th { tr.put_pixel(i as u32, (th-1-my) as u32, image::Rgb([60,200,230])); }
            let c = if h.data[k]<=SEA_LEVEL_NORM {[30,50,90]} else {biomes[k].color()};
            for y in 0..30 { tr.put_pixel(i as u32, (th-1-y) as u32, image::Rgb(c)); }
        }
        tr.save(dir.join(format!("seed{seed:05}_transect_row{jr}_maxmm{max_mm:.0}.png"))).unwrap();
    }
    eprintln!("  out = {}", dir.display());
    eprintln!("  relief: grey=alt(0-6000m), red tint=windward 80km band. transect: grey=alt(0-6000m) cyan=precip(0-rowmax mm) bottom=biome strip.");
}
