//! DIAGNOSTIC — is the C1 precipitation transport RESOLUTION-INVARIANT?
//!
//! Motivation (production measurement, arid-hot 25°, seed 10481999410520546993, 400 km):
//! mean precipitation over land is **1219 mm/yr at 2048² and 712 mm/yr at 8192²** (×0.58),
//! and the net runoff `max(0, p − pe)` follows at ×0.33. That collapse propagates into
//! `segment_drainage_km2` (= runoff accumulation / 300 mm), hence into the discharge, the
//! channel width and the navigability class — every hydrological figure the export carries.
//!
//! The production comparison CONFOUNDS two things: the terrain itself differs between
//! resolutions (finer FBM detail, more closed depressions) and the transport is discretised
//! differently. This test removes the confound: the SAME ANALYTIC terrain is built at several
//! resolutions over the SAME physical domain, so any change in the land-mean precipitation
//! is discretisation alone.
//!
//! HYPOTHESIS (mine, from reading the code) — **REFUTED by this test.** The orographic
//! extraction is `oro = k_oro · m · ascent` with `ascent = Δalt / dx_m`, a FRACTION of the
//! carried moisture per cell, so the flux should decay as `(1 − k_oro·S)^N` and compound per
//! CELL instead of per KILOMETRE. The measurement says otherwise: **391 mm/yr at every grid
//! from 512² to 8192², ratio 1.000**, and the interior/coast contrast is identical. The
//! reason the hypothesis fails is that `k_oro` is not what wrings the air out — the CAPACITY
//! CAP (`m > e_sat(T)`, air cannot carry more than saturation) is, and that is a function of
//! the altitude PROFILE, not of how many cells sample it. The ADR already noted this ("the
//! decisive re-measure FALSIFIED lowering the uplift extraction ... the cold coastal peak
//! wrings it out regardless"); this test makes it a resolution statement.
//!
//! CONSEQUENCE: the production gap is NOT in the climate code. It is in the TERRAIN the
//! climate reads — mean land altitude 287 m at 2048² against 693 m at 8192² (×2.42).
//!
//! Run: cargo test -p ymir-core --test precip_resolution_invariance -- --ignored --nocapture

use ymir_core::climate::precipitation::{PrecipParams, compute_precipitation, precip_mm_per_year};
use ymir_core::climate::temperature::compute_temperature_span;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;

const DOMAIN_KM: f32 = 400.0;
const LAT_DEG: f32 = 25.0;
const FULL_RANGE_M: f32 = 11302.0;

/// One analytic continent, sampled at `n²`: ocean on the upwind third, then a coastal plain
/// rising linearly to a 2000 m interior plateau. Identical GEOMETRY at every `n` — the only
/// thing that changes is how many cells discretise it.
fn analytic_terrain(n: usize) -> GridF32 {
    let mut g = GridF32::new(n, n, 0.3); // ocean baseline (below 0.5)
    for j in 0..n {
        for i in 0..n {
            let u = i as f32 / n as f32; // along-wind fraction of the domain
            let v = if u < 0.30 {
                0.30 // ocean
            } else if u < 0.55 {
                // coastal ramp: 0 m → 2000 m over 100 km
                let t = (u - 0.30) / 0.25;
                0.5 + t * (2000.0 / FULL_RANGE_M)
            } else {
                0.5 + 2000.0 / FULL_RANGE_M // interior plateau at 2000 m
            };
            g.set(i, j, v);
        }
    }
    g
}

struct Row {
    n: usize,
    km_per_cell: f32,
    land_mean_mm: f64,
    coastal_mm: f64,
    interior_mm: f64,
}

fn measure(n: usize) -> Row {
    let hm = analytic_terrain(n);
    let alt = |x: f32| (x - 0.5) * FULL_RANGE_M;
    let ss = SteinSteinParams::default();
    let temp = compute_temperature_span(&hm, &ss, LAT_DEG, DOMAIN_KM / 111.0);
    let km_per_cell = DOMAIN_KM / n as f32;
    let p = compute_precipitation(&hm, &temp, LAT_DEG, km_per_cell, alt, &PrecipParams::default());

    let (mut sum, mut cnt) = (0.0f64, 0usize);
    // The ramp toe (just inland of the coast) against the far interior — the rain-shadow
    // signature the mechanism predicts: the finer the grid, the more front-loaded the rain.
    let (mut c_sum, mut c_cnt) = (0.0f64, 0usize);
    let (mut i_sum, mut i_cnt) = (0.0f64, 0usize);
    for j in 0..n {
        for i in 0..n {
            if hm.get(i as i32, j as i32) <= 0.5 {
                continue;
            }
            let mm = precip_mm_per_year(p.get(i as i32, j as i32)) as f64;
            sum += mm;
            cnt += 1;
            let u = i as f32 / n as f32;
            if (0.30..0.40).contains(&u) {
                c_sum += mm;
                c_cnt += 1;
            } else if u >= 0.85 {
                i_sum += mm;
                i_cnt += 1;
            }
        }
    }
    Row {
        n,
        km_per_cell,
        land_mean_mm: sum / cnt.max(1) as f64,
        coastal_mm: c_sum / c_cnt.max(1) as f64,
        interior_mm: i_sum / i_cnt.max(1) as f64,
    }
}

#[test]
#[ignore]
fn precip_transport_resolution_invariance() {
    eprintln!(
        "\n=====  PRECIPITATION TRANSPORT — SAME ANALYTIC CONTINENT, FOUR GRIDS  =====\n\
         domain {DOMAIN_KM} km, lat {LAT_DEG}°, ocean 0–120 km then a 100 km ramp to a \
         2000 m plateau"
    );
    eprintln!(
        "\n{:>6} {:>10} {:>14} {:>14} {:>14} {:>10}",
        "n", "m/cell", "land mean", "coast 120-160", "interior >340", "lee/coast"
    );
    let rows: Vec<Row> = [512, 1024, 2048, 4096, 8192].iter().map(|&n| measure(n)).collect();
    for r in &rows {
        eprintln!(
            "{:>6} {:>10.0} {:>11.0} mm {:>11.0} mm {:>11.0} mm {:>10.3}",
            r.n,
            r.km_per_cell * 1000.0,
            r.land_mean_mm,
            r.coastal_mm,
            r.interior_mm,
            r.interior_mm / r.coastal_mm.max(1e-9)
        );
    }
    let first = &rows[0];
    let last = &rows[rows.len() - 1];
    eprintln!(
        "\nland mean {:.0} → {:.0} mm/yr over a ×{} refinement (ratio {:.3})",
        first.land_mean_mm,
        last.land_mean_mm,
        last.n / first.n,
        last.land_mean_mm / first.land_mean_mm.max(1e-9)
    );
    let (c0, c1) = (
        first.interior_mm / first.coastal_mm.max(1e-9),
        last.interior_mm / last.coastal_mm.max(1e-9),
    );
    let drift = (last.land_mean_mm / first.land_mean_mm.max(1e-9) - 1.0).abs();
    eprintln!(
        "interior/coast contrast {c0:.4} → {c1:.4}\n\nVERDICT: {}",
        if drift > 0.02 {
            "the land mean MOVES with the grid — the transport is not resolution-invariant"
        } else {
            "the transport IS resolution-invariant on fixed geometry (drift < 2 %). The \
             per-cell-extraction hypothesis is REFUTED: `ascent = Δalt/dx_m` is the same \
             physical gradient at every grid, and the CAPACITY CAP `m > e_sat(T)` — not \
             `k_oro` — is what wrings the air out, so the result depends on the altitude \
             profile and not on how many cells discretise it. Any production difference in \
             mean precipitation is therefore caused by the TERRAIN differing between \
             resolutions, not by the climate code."
        }
    );
    // Diagnostic only: report, never gate. A verdict is drawn from the numbers above.
}
