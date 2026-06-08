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
use ymir_core::tectonics_c1::production_upscale::upscale_from_c1;

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
