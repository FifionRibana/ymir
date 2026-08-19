//! M1 #190 regression diagnostic — seed 42 land collapse + speckled coasts.
//! Read-only measurement (no defaults changed here): runs ONE coarse tectonic
//! pass for seed 42 and evaluates the land topology + hypsometry across
//! target_land_fraction values, to attribute the regression (S1/S2/S4).
//! Run: cargo test -p ymir-core --test m1_regression_seed42 -- --ignored --nocapture

use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::land_topology::land_topology;
use ymir_core::tectonics_c1::production_upscale::{c1_coarse_raw_altitude, c1_normalize_coarse};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};

/// Local slope of the hypsometric curve at altitude `thr`: cells per unit raw
/// altitude within ±`win` (a flat curve = many cells at nearly the same height
/// = the 0-crossing is hypersensitive → speckled coast).
fn slope_at(sorted: &[f32], thr: f32, win: f32) -> f32 {
    let lo = sorted.partition_point(|&v| v < thr - win);
    let hi = sorted.partition_point(|&v| v < thr + win);
    (hi - lo) as f32 / (2.0 * win)
}

fn quantile(sorted: &[f32], q: f32) -> f32 {
    let idx = ((q * (sorted.len() - 1) as f32).round() as usize).min(sorted.len() - 1);
    sorted[idx]
}

fn coarse_raw(seed: u64, grid: usize, init: &Phase2InitParams) -> (ymir_core::grid::GridF32, C1TimeLoopConfig) {
    let run = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / grid as f64,
        dy: 1.0 / grid as f64,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let ss = SteinSteinParams::default();
    let clo = C1Closures::default();
    let mut state = init_c1_state_phase_2_r7(grid, seed, init);
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, &run, &clo, |_, _| {});
    (c1_coarse_raw_altitude(&state, &run.iso_config, &ss), run)
}

fn report(label: &str, raw: &ymir_core::grid::GridF32, tlf: Option<f32>) {
    let norm = c1_normalize_coarse(raw.clone(), tlf);
    let t = land_topology(&norm, 0.5);
    let tlf_s = tlf.map(|f| format!("{f:.2}")).unwrap_or_else(|| "None".into());
    eprintln!(
        "  [{label}] tlf={tlf_s:>4}  emerged {:>5.1}%  largest {:>7.0} km² ({:>4.1}%)  \
         traverse {:>3.0} km  masses {:>3}  wrap x={} y={}",
        t.emerged_fraction * 100.0,
        t.largest_area_km2,
        t.largest_area_frac * 100.0,
        t.bbox_km.0.max(t.bbox_km.1),
        t.num_landmasses,
        t.wraps_x,
        t.wraps_y,
    );
}

#[test]
#[ignore]
fn seed42_land_and_hypsometry() {
    let (seed, grid) = (42u64, 64usize);

    // Pre-M1 tectonic defaults (num_plates 8, seed_cluster_count 1).
    let pre = Phase2InitParams::default();
    eprintln!(
        "pre-M1 tectonic defaults: num_plates={} seed_cluster_count={}",
        pre.num_plates, pre.cluster.seed_cluster_count
    );
    let (raw, _run) = coarse_raw(seed, grid, &pre);

    // STEP 2/3 — land topology across tlf at pre-M1 tectonics (isolates S1/S2).
    // (Erosion sea_level 0.1↔0.5 is an HD-stage coastal-deposition knob; it does
    //  NOT affect the coarse emerged fraction — noted, not varied here.)
    eprintln!("=== seed 42, num_plates 8 / cc 1 (pre-M1 tectonics) — land vs tlf ===");
    for tlf in [None, Some(0.29f32), Some(0.15), Some(0.10), Some(0.08)] {
        report("8pl", &raw, tlf);
    }

    // STEP 4 — hypsometry: where each sea level sits on the coarse curve, and the
    // local slope (cells per raw-altitude unit; low = flat plateau = speckle).
    let mut sorted = raw.data.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let win = 0.02f32;
    eprintln!("=== hypsometry (seed 42) — sea-level placement + local slope ===");
    eprintln!(
        "  isostatic sea (tlf None) = raw 0.00 → slope {:.0} cells/unit",
        slope_at(&sorted, 0.0, win)
    );
    for tlf in [0.29f32, 0.15, 0.10, 0.08] {
        let thr = quantile(&sorted, 1.0 - tlf);
        eprintln!(
            "  tlf {tlf:.2} → threshold raw {thr:+.3} → slope {:.0} cells/unit",
            slope_at(&sorted, thr, win)
        );
    }

    // S4 — current UI default tectonics (island_production: 16 plates / cc 3) at
    // the current default tlf 0.08, for reference (what the user actually sees).
    let mut cur = Phase2InitParams::default();
    cur.num_plates = 16;
    cur.cluster.seed_cluster_count = 3;
    let (raw16, _) = coarse_raw(seed, grid, &cur);
    eprintln!("=== seed 42, num_plates 16 / cc 3 (current island_production) ===");
    report("16pl", &raw16, Some(0.08));
    report("16pl", &raw16, None);
}
