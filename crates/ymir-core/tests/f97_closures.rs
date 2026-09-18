//! ADR 0001 Finding 97 block B — two ISOTROPIC candidates, one per defect, judged with block A's
//! anisotropy column in the table.
//!
//! Finding 96 separated the two defects and named their owners, then had its best candidate
//! refuted by an image: the χ floor is a **D8 path integral**, so clamping to it prints the D8
//! axes into the terrain. Both candidates here drop the path integral.
//!
//!  * **B1 — the area cap**, `h_r_eff = max(h_r, min(h_pre − c·A^q, h_o))`. A pure function of the
//!    pre-incision surface and the drainage area, so it needs **no production change at all**: it
//!    is passed through the `incision_floor` seam Finding 96 already opened. `c` and `q` are read
//!    off Finding 88-D2 and declared below.
//!  * **B2 — the local equilibrium-slope floor**, `h_r_eff = min(h_r + S_eq(A)·dx, h_o)`. The
//!    differential form of χ. ⚠️ This one DOES need a seam (`slope_floor_uk`), because it is
//!    relative to the receiver at run time and cannot be a fixed grid. It is additive and
//!    `None` in production.
//!  * **B3′ — A1 + B2**, the two defects' owners together.
//!
//! Run: cargo test -p ymir-core --release --test f97_closures -- --ignored --nocapture

mod common;

use common::{
    CELL_KM, Crit, Knobs, SEA, aniso, build_field, build_field_with_floor, f95_criteria, land_u16,
    majority, pct, slope_deg, sorted, to_mask, worst_tile,
};
use std::time::Instant;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::coast_metrics::{MIN_SPUR_KM, NECK_KM, coast_spurs};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, breach_monotone};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const ORACLE_P50_M: f32 = 488.1;
const DELIVERED_P50_M: f32 = 424.2;

/// ADR Finding 88-D2, the two published points, and the fit between them. **Population declared**:
/// they are medians over ONE bowl's footprint (234 724 cells, humid, shipped), not over all land —
/// so `c` and `q` are an order-of-magnitude reading of "how the cut scales with area", which is
/// what Finding 88-D2 itself claims (*"the cut scales with drainage area — ×2.0 between the two
/// decades measured"*), and not a calibration.
const F88_CUT_AT_1KM2_M: f32 = 308.9;
const F88_CUT_BELOW_001KM2_M: f32 = 151.6;

fn spurs_km(mask: &[bool], w: usize) -> usize {
    let (m, ww) = majority(mask, w, 1);
    coast_spurs(&marching_squares(&to_mask(&m, ww), 0.5), CELL_KM, MIN_SPUR_KM, NECK_KM).0.len()
}

#[test]
#[ignore]
fn f97_closures() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    eprintln!("\n==========  Finding 97 B . two isotropic candidates  ==========");

    let t0 = Instant::now();
    let pre = build_field(Knobs::no_incision());
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let t1 = Instant::now();
    let delivered = build_field(Knobs::passes(2));
    let t_del = t1.elapsed().as_secs_f64();
    eprintln!(
        "   builds: pre-incision {:.1} s . delivered {t_del:.1} s",
        t1.duration_since(t0).as_secs_f64()
    );

    // ── B1's floor: h_pre − c·A^q, from the PRE-INCISION breached network ────
    let d_pre = c1_drainage_windowed(&pre, None, &on, &ss, DOMAIN_KM);
    let br = breach_monotone(&pre, &d_pre.flow.filled, &d_pre.lake_map, SEA, w, h);
    let d_br = c1_drainage_windowed(&br, None, &on, &ss, DOMAIN_KM);
    let q = (F88_CUT_BELOW_001KM2_M / F88_CUT_AT_1KM2_M).ln() / (0.01f32).ln();
    let c_cap = F88_CUT_AT_1KM2_M;
    eprintln!(
        "\n-- B1 . the area cap, parameters read off Finding 88-D2 --\n   c = {c_cap:.1} m at \
         A = 1 km2, {F88_CUT_BELOW_001KM2_M:.1} m at A = 0.01 km2 => **q = {q:.4}** . at one cell \
         ({cell_km2:.5} km2) the cap allows {:.1} m of cut, at 1000 km2 {:.1} m",
        c_cap * cell_km2.powf(q),
        c_cap * 1000f32.powf(q)
    );
    let floor_b1: Vec<f32> = (0..n)
        .map(|k| {
            if pre.data[k] > SEA {
                let a = (d_br.flow.accumulation.data[k] * cell_km2).max(cell_km2);
                pre.data[k] - (c_cap * a.powf(q)) / n2m
            } else {
                f32::NEG_INFINITY
            }
        })
        .collect();

    // ── B2's one free parameter: the Flint intercept of the DELIVERED channels ──
    //
    // Declared calibration: `uk = median over delivered channel cells (A >= A_c) of S·sqrt(A_km2)`
    // — i.e. the value that reproduces the delivered field's own median channel slope, so the
    // floor cannot be accused of importing a relief target from outside.
    let d_del = c1_drainage_windowed(&delivered, None, &on, &ss, DOMAIN_KM);
    let mut ks: Vec<f32> = Vec::new();
    for k in (0..n).step_by(7) {
        if delivered.data[k] <= SEA {
            continue;
        }
        let a_km2 = d_del.flow.accumulation.data[k] * cell_km2;
        if a_km2 < RELIEF_V1_A_C_KM2 {
            continue;
        }
        let d = d_del.flow.direction[k];
        if d == DIR_NONE {
            continue;
        }
        let nx = ((k % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
        let ny = ((k / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
        let diag = D8_DX[d as usize] != 0 && D8_DY[d as usize] != 0;
        let dx_m = if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M };
        let s = (delivered.data[k] - delivered.data[ny * w + nx]).max(0.0) * n2m / dx_m;
        if s > 0.0 {
            ks.push(s * a_km2.sqrt());
        }
    }
    let uk = pct(&sorted(ks.clone()), 0.50);
    eprintln!(
        "\n-- B2 . the local slope floor, calibrated on the DELIVERED channels --\n   n = {} \
         channel cells (stride 7, A >= A_c) => **(U/K)^(1/n) = k_s = {uk:.4}** . for reference \
         Finding 96's chi calibration was 0.0719 . S_eq at A_c = {:.4} m/m ({:.1} deg), at \
         1000 km2 = {:.5} m/m ({:.2} deg)",
        ks.len(),
        uk / RELIEF_V1_A_C_KM2.sqrt(),
        (uk / RELIEF_V1_A_C_KM2.sqrt()).atan().to_degrees(),
        uk / 1000f32.sqrt(),
        (uk / 1000f32.sqrt()).atan().to_degrees()
    );

    // ── the three candidates ────────────────────────────────────────────────
    let tb1 = Instant::now();
    let b1 = build_field_with_floor(Knobs::passes(2), Some(std::sync::Arc::new(floor_b1)));
    let t_b1 = tb1.elapsed().as_secs_f64();
    let tb2 = Instant::now();
    let b2 = build_field(Knobs { slope_floor_uk: Some(uk), ..Knobs::passes(2) });
    let t_b2 = tb2.elapsed().as_secs_f64();
    let ta1 = Instant::now();
    let a1 = build_field(Knobs { depression_floor: true, ..Knobs::passes(2) });
    let t_a1 = ta1.elapsed().as_secs_f64();
    let tb3 = Instant::now();
    let b3p =
        build_field(Knobs { slope_floor_uk: Some(uk), depression_floor: true, ..Knobs::passes(2) });
    let t_b3p = tb3.elapsed().as_secs_f64();
    eprintln!(
        "\n   built: B1 {t_b1:.1} s (**{:+.1}**) . B2 {t_b2:.1} s (**{:+.1}**) . A1+B2 \
         {t_b3p:.1} s (**{:+.1}**) . A1 {t_a1:.1} s (**{:+.1}**)",
        t_b1 - t_del,
        t_b2 - t_del,
        t_b3p - t_del,
        t_a1 - t_del
    );

    // ── the table, one row per candidate ────────────────────────────────────
    let auth = spurs_km(&land_u16(&pre, &ss), w);
    let cut_of = |f: &GridF32| -> f64 {
        (0..n)
            .filter(|&k| pre.data[k] > SEA && f.data[k] > SEA)
            .map(|k| ((pre.data[k] - f.data[k]) * n2m).max(0.0) as f64)
            .sum()
    };
    let cut_del = cut_of(&delivered);
    let mut rows: Vec<(String, Crit, f64, f64, f32)> = Vec::new();
    for (label, f, cost) in [
        ("delivered", &delivered, t_del),
        ("B1 area cap", &b1, t_b1),
        ("A1 depression skip", &a1, t_a1),
        ("B2 slope floor", &b2, t_b2),
        ("A1+B2", &b3p, t_b3p),
    ] {
        let d = c1_drainage_windowed(f, None, &on, &ss, DOMAIN_KM);
        let bf = breach_monotone(f, &d.flow.filled, &d.lake_map, SEA, w, h);
        let c = f95_criteria(f, &bf, &pre, DELIVERED_P50_M, &ss, &on, cell_km2, n2m, w, h);
        let land: Vec<bool> = (0..n).map(|k| f.data[k] > SEA).collect();
        let a16 = aniso(f, &land, 16);
        let (bx, by, at) = worst_tile(f, &land);
        let sl = sorted(
            (0..n).step_by(37).filter(|&k| land[k]).map(|k| slope_deg(f, k, n2m, w, h)).collect(),
        );
        let cut = cut_of(f);
        let sp_br = spurs_km(&land_u16(&bf, &ss), w);
        let sp_er = spurs_km(&land_u16(f, &ss), w);
        eprintln!(
            "   **{label}**: coast breached+u16 **{sp_br} ({:+})** . eroded u16 {sp_er} ({:+}) . \
             slope p90 {:.1} deg . share>30 {:.2} % . **erosion permitted {:.1} %** . ANISO w16: \
             C p50 {:.3} . C>0.7 {:.2} % . H(theta) {:.4} . R2 {:.4} . **R8 {:.4}** . worst tile \
             ({bx},{by}) R8 **{:.4}**",
            sp_br as i64 - auth as i64,
            sp_er as i64 - auth as i64,
            pct(&sl, 0.90),
            100.0 * sl.iter().filter(|&&x| x > 30.0).count() as f64 / sl.len().max(1) as f64,
            100.0 * cut / cut_del,
            a16.c_p50,
            a16.share_hi,
            a16.h_theta,
            a16.r2,
            a16.r8,
            at.r8
        );
        rows.push((label.to_string(), c, cost, 100.0 * cut / cut_del, a16.r8));
    }

    eprintln!("\n-- C . the table --");
    eprintln!(
        "   | candidate | canyons (rate) | to_nothing | lakes % | paired p50 | vs oracle | \
         erosion permitted | R8 | cost |"
    );
    eprintln!("   |---|---|---|---|---|---|---|---|---|");
    for (label, c, cost, perm, r8) in &rows {
        eprintln!(
            "   | {label} | {}/{} ({:.1} %) | {} | {:.2} | {:.1} m | {:+.1} % | {perm:.1} % | \
             {r8:.4} | {cost:.1} s |",
            c.klass,
            c.scanned,
            100.0 * c.klass as f32 / c.scanned.max(1) as f32,
            c.to_nothing,
            c.lake_pct,
            c.p50,
            100.0 * (c.p50 - ORACLE_P50_M) / ORACLE_P50_M
        );
    }
    eprintln!(
        "   | **ORACLE (F95, 300 passes)** | 0/24 (0.0 %) | 0 | 15.92 | {ORACLE_P50_M:.1} m | - | \
         100 % | n/a | 6 210 s |"
    );
    eprintln!("\n==========  end Finding 97 B  ==========\n");
}
