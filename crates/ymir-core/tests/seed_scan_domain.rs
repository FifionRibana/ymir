//! M1 #190 — whole-domain seed scan (the domain IS the map, no crop).
//! For each seed, one cheap coarse pass (pre-M1 defaults: 8 plates, cc 1, tlf
//! None), then the domain-as-map verdict: not a band, no land on the domain
//! border, continent ≤ 60 % of the domain on both axes, deep-enough border ring.
//! Reports the first passing seeds with their bbox %, ocean margins and border
//! depth, plus the bathymetric asymptote distance (the real margin requirement)
//! and the coupled-vs-uncoupled slope telemetry.
//!
//! Run: cargo test -p ymir-core --test seed_scan_domain --release -- --ignored --nocapture

use ymir_core::grid::GridF32;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::land_topology::{DomainMetrics, domain_metrics, slope_shares};
use ymir_core::tectonics_c1::production_upscale::{
    C1_DOMAIN_KM, c1_coarse_raw_altitude, c1_normalize_coarse,
};
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};

const GRID: usize = 64;
const TARGET_SIZE: usize = 8192; // whole-domain export side
const SEA: f32 = 0.5;

fn coarse_norm(seed: u64, init: &Phase2InitParams, ss: &SteinSteinParams) -> GridF32 {
    let run = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / GRID as f64,
        dy: 1.0 / GRID as f64,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let clo = C1Closures::default();
    let mut state = init_c1_state_phase_2_r7(GRID, seed, init);
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, &run, &clo, |_, _| {});
    let raw = c1_coarse_raw_altitude(&state, &run.iso_config, &ss);
    c1_normalize_coarse(raw, None) // tlf None (isostatic sea, crisp coasts)
}

fn row(seed: u64, m: &DomainMetrics) -> String {
    format!(
        "seed {seed:>3}  geom {g}  {v:<24}  bbox {bx:>4.0}×{by:>4.0}% ext {ex:>4.0}×{ey:>4.0}km  \
         margins N{mn:>4.0} S{ms:>4.0} E{me:>4.0} W{mw:>4.0}km  \
         border {bd:>5.0}m ({bf:>3.0}% asymp) deepest {dp:>5.0}m  masses {nm:>3}  \
         emerged {em:>4.1}%(+{dr:.1})",
        g = if m.geometric_pass { "PASS" } else { "fail" },
        v = if m.verdict_pass { "PASS".into() } else { format!("FAIL: {}", m.verdict_reason) },
        bx = m.bbox_frac_x * 100.0,
        by = m.bbox_frac_y * 100.0,
        ex = m.extent_km.0,
        ey = m.extent_km.1,
        mn = m.margin_n_km,
        ms = m.margin_s_km,
        me = m.margin_e_km,
        mw = m.margin_w_km,
        bd = m.border_depth_median_m,
        bf = m.border_depth_frac * 100.0,
        dp = m.deepest_ocean_m,
        nm = m.topo.num_landmasses,
        em = m.emerged_frac * 100.0,
        dr = m.emerged_drift * 100.0,
    )
}

#[test]
#[ignore]
fn scan_domain_map_seeds() {
    let ss = SteinSteinParams::default();
    let init = Phase2InitParams::default(); // pre-M1: 8 plates, cc 1
    eprintln!(
        "pre-M1 defaults: {} plates, cc {}; domain {} km; target {}² → {:.0} m/cell whole-domain",
        init.num_plates,
        init.cluster.seed_cluster_count,
        C1_DOMAIN_KM,
        TARGET_SIZE,
        C1_DOMAIN_KM / TARGET_SIZE as f32 * 1000.0,
    );

    // Seed 42 explicitly (author's reference — expected to FAIL on bbox width).
    let m42 = domain_metrics(&coarse_norm(42, &init, &ss), SEA, &ss, C1_DOMAIN_KM, TARGET_SIZE);
    eprintln!("\n=== seed 42 (author reference) ===\n  {}", row(42, &m42));

    // Scan.
    const N_SCAN: u64 = 200;
    const N_REPORT: usize = 12;
    let mut geom_pass: Vec<(u64, DomainMetrics)> = Vec::new();
    let mut full_pass = 0usize;
    let mut asymp_cells: Vec<usize> = Vec::new();
    let mut deepest: Vec<f32> = Vec::new();
    for seed in 0..N_SCAN {
        let m = domain_metrics(&coarse_norm(seed, &init, &ss), SEA, &ss, C1_DOMAIN_KM, TARGET_SIZE);
        if let Some(c) = m.cells_to_asymptote {
            asymp_cells.push(c);
        }
        deepest.push(m.deepest_ocean_m);
        if m.verdict_pass {
            full_pass += 1;
        }
        if m.geometric_pass {
            geom_pass.push((seed, m));
        }
    }

    eprintln!(
        "\n=== scan of seeds 0..{N_SCAN} — {} geometric-pass ({:.0}%), {} full-pass ({:.0}%) ===",
        geom_pass.len(),
        geom_pass.len() as f32 / N_SCAN as f32 * 100.0,
        full_pass,
        full_pass as f32 / N_SCAN as f32 * 100.0,
    );
    eprintln!("  first {} geometric-pass seeds (not band, no border land, ≤60%/axis):", N_REPORT.min(geom_pass.len()));
    for (seed, m) in geom_pass.iter().take(N_REPORT) {
        eprintln!("  {}", row(*seed, m));
    }
    if geom_pass.is_empty() {
        eprintln!("  (none — every seed fails a geometric clause at these defaults)");
    }

    // The bathymetric ceiling: the theoretical Stein-Stein asymptote vs what the
    // (young) C1 ocean actually reaches. If the asymptote is never reached, the
    // "≥90% of asymptote" bar is a model ceiling, not a seed-selection knob.
    deepest.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let dmed = deepest[deepest.len() / 2];
    let dmax = deepest[deepest.len() - 1];
    eprintln!(
        "\n=== oceanic depth ceiling over {} seeds (Stein-Stein asymptote {:.0} m) ===\n  \
         deepest ocean cell: median {dmed:.0} m ({:.0}% of asymptote), max {dmax:.0} m ({:.0}%)\n  \
         asymptote (≥90%) reached in-domain by {} / {} seeds",
        deepest.len(),
        ss.asymptotic_depth_m,
        dmed / ss.asymptotic_depth_m as f32 * 100.0,
        dmax / ss.asymptotic_depth_m as f32 * 100.0,
        asymp_cells.len(),
        N_SCAN,
    );
    if !asymp_cells.is_empty() {
        asymp_cells.sort_unstable();
        let km = C1_DOMAIN_KM / GRID as f32;
        let med = asymp_cells[asymp_cells.len() / 2];
        eprintln!(
            "  where reached: median {med} coarse cells from the coast ({:.0} km) — the model's real margin need",
            med as f32 * km,
        );
    }

    // Slope telemetry (TASK 5): buildable-land cost of NOT coupling depth_scale to
    // domain_km, measured on seed 42's coarse field at a dramatic 375 km domain.
    let f42 = coarse_norm(42, &init, &ss);
    let base = ss.depth_scale_m as f32;
    let dramatic = 375.0f32;
    let coupled = slope_shares(&f42, SEA, dramatic, base * dramatic / C1_DOMAIN_KM);
    let uncoupled = slope_shares(&f42, SEA, dramatic, base);
    let ref1024 = slope_shares(&f42, SEA, C1_DOMAIN_KM, base);
    eprintln!(
        "\n=== slope telemetry (seed 42 coarse land) — share of LAND above 15°/30°/45° ===\n  \
         domain 1024 km (baseline)        : {:>4.1}% / {:>4.1}% / {:>4.1}%\n  \
         domain  375 km, COUPLED depth    : {:>4.1}% / {:>4.1}% / {:>4.1}%  (slopes preserved)\n  \
         domain  375 km, UNCOUPLED depth  : {:>4.1}% / {:>4.1}% / {:>4.1}%  (×2.73 gradients)",
        ref1024.0 * 100.0, ref1024.1 * 100.0, ref1024.2 * 100.0,
        coupled.0 * 100.0, coupled.1 * 100.0, coupled.2 * 100.0,
        uncoupled.0 * 100.0, uncoupled.1 * 100.0, uncoupled.2 * 100.0,
    );
}
