//! ADR 0001 Finding 98 — is the wall a wall? **No production change.**
//!
//! Finding 97 named a wall: *the defects and the erosion are the same quantity; only integration
//! and deposition redistribute.* It named it BEFORE the two measurements that test it. This is
//! those two, plus the line the table was missing.
//!
//!  * **A** — where the removed erosion went. The round partitions by position relative to χ;
//!    ⚠️ that reference is a property of the FIELD, not of the removal (most land already sits
//!    below χ), so the partition is also drawn **on the closures' own gates**, which is what
//!    distinguishes *preventing* from *suppressing*.
//!  * **B** — the dial: A1+B2 at `k_time` ×1.5 and ×2, with the delivered field at the same
//!    budgets as the control and as the honest denominator for rule 14.
//!  * **C** — ⚠️ the round calls this "anti-comb (F78)". `anti-comb` has **0 hits** in the
//!    dossier; Finding 78's lever is the **shelf clamp**, it runs AFTER the incision
//!    (`incise_lithology` at `production_upscale.rs:451`, `apply_bathymetry_profile` at `:506`),
//!    and Finding 79 deliberately took it from 20 m to **1 m in production** with an
//!    `ALGO_UPSCALE_EROSION` bump. So it CANNOT move R8 on land, the paired relief, or the erosion
//!    permitted — those three are measured anyway, as a negative control on that reasoning.
//!
//! Run: cargo test -p ymir-core --release --test f98_wall -- --ignored --nocapture

mod common;

use common::{
    CELL_KM, Crit, Knobs, SEA, aniso, build_field, build_field_with_floor, f95_criteria, land_u16,
    majority, pct, sorted, to_mask,
};
use std::time::Instant;
use ymir_core::erosion::stream_power::{RELIEF_V1_A_C_KM2, RELIEF_V1_K, RELIEF_V3_K_MULT};
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::coast_metrics::{MIN_SPUR_KM, NECK_KM, coast_spurs};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, breach_monotone, propagation_order};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const ORACLE_P50_M: f32 = 488.1;
const DELIVERED_P50_M: f32 = 424.2;
/// Finding 97 block B's declared calibration: the median of `S·√A` over the DELIVERED channel
/// cells, i.e. the delivered field's own Flint intercept.
const KS: f32 = 0.0451;
/// Finding 88-D2, with the population declared at Finding 97 B.
const CAP_C: f32 = 308.9;

fn spurs(mask: &[bool], w: usize) -> usize {
    let (m, ww) = majority(mask, w, 1);
    coast_spurs(&marching_squares(&to_mask(&m, ww), 0.5), CELL_KM, MIN_SPUR_KM, NECK_KM).0.len()
}

#[test]
#[ignore]
fn f98_wall() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    let k0 = RELIEF_V1_K * RELIEF_V3_K_MULT;
    eprintln!("\n==========  Finding 98 . is the wall a wall?  ==========");
    eprintln!("   shipped k = {k0:.0}, dt = 1, iterations = 2 => k_time = {:.0}", k0 * 2.0);

    let t0 = Instant::now();
    let pre = build_field(Knobs::no_incision());
    let (w, h) = (pre.width, pre.height);
    let n = w * h;

    // the χ reference and B1's floor both need the pre-incision breached network
    let d_pre = c1_drainage_windowed(&pre, None, &on, &ss, DOMAIN_KM);
    let br = breach_monotone(&pre, &d_pre.flow.filled, &d_pre.lake_map, SEA, w, h);
    let d_br = c1_drainage_windowed(&br, None, &on, &ss, DOMAIN_KM);
    let land_pre: Vec<bool> = (0..n).map(|k| br.data[k] > SEA).collect();
    let (order, un) = propagation_order(&d_br.flow.direction, |k| land_pre[k], w, h);
    assert_eq!(un, 0);
    let mut chi = vec![0.0f32; n];
    for &ku in order.iter().rev() {
        let k = ku as usize;
        if !land_pre[k] {
            continue;
        }
        let d = d_br.flow.direction[k];
        if d == DIR_NONE {
            continue;
        }
        let nx = ((k % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
        let ny = ((k / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
        let r = ny * w + nx;
        let diag = D8_DX[d as usize] != 0 && D8_DY[d as usize] != 0;
        let dx_m = if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M };
        let a = (d_br.flow.accumulation.data[k] * cell_km2).max(RELIEF_V1_A_C_KM2);
        chi[k] = if land_pre[r] { chi[r] } else { 0.0 } + (1.0 / a).sqrt() * dx_m;
    }
    let cs = sorted((0..n).filter(|&k| land_pre[k]).map(|k| chi[k]).collect::<Vec<f32>>());
    let c_uk = ORACLE_P50_M / pct(&cs, 0.50).max(1e-6);
    let z_chi: Vec<f32> = (0..n)
        .map(|k| if land_pre[k] { SEA + (c_uk * chi[k]) / n2m } else { br.data[k] })
        .collect();
    // B1's area cap, rebuilt from the same accumulation (Finding 97 B)
    let q = (151.6f32 / CAP_C).ln() / (0.01f32).ln();
    let floor_b1: std::sync::Arc<Vec<f32>> = std::sync::Arc::new(
        (0..n)
            .map(|k| {
                if pre.data[k] > SEA {
                    let a = (d_br.flow.accumulation.data[k] * cell_km2).max(cell_km2);
                    pre.data[k] - (CAP_C * a.powf(q)) / n2m
                } else {
                    f32::NEG_INFINITY
                }
            })
            .collect(),
    );
    // the PRE-INCISION depressions: A1's population, as a declared PROXY (A1 re-evaluates the
    // gate on the CURRENT surface every iteration; this is the surface it sees on iteration 1)
    let in_pit: Vec<bool> =
        (0..n).map(|k| land_pre[k] && d_pre.flow.filled.data[k] > pre.data[k] + 1e-7).collect();
    eprintln!(
        "   prerequisites in {:.1} s . pre-incision pits: {} cells",
        t0.elapsed().as_secs_f64(),
        in_pit.iter().filter(|&&b| b).count()
    );

    // ── the fields ───────────────────────────────────────────────────────────
    let ab = |mult: f32| Knobs {
        slope_floor_uk: Some(KS),
        depression_floor: true,
        k: Some(k0 * mult),
        ..Knobs::passes(2)
    };
    let plain = |mult: f32| Knobs { k: Some(k0 * mult), ..Knobs::passes(2) };
    let mut t = Instant::now();
    let del1 = build_field(plain(1.0));
    let t_del1 = t.elapsed().as_secs_f64();
    t = Instant::now();
    let ab1 = build_field(ab(1.0));
    let t_ab1 = t.elapsed().as_secs_f64();
    t = Instant::now();
    let ab15 = build_field(ab(1.5));
    let t_ab15 = t.elapsed().as_secs_f64();
    t = Instant::now();
    let ab20 = build_field(ab(2.0));
    let t_ab20 = t.elapsed().as_secs_f64();
    t = Instant::now();
    let del15 = build_field(plain(1.5));
    let t_del15 = t.elapsed().as_secs_f64();
    let del20 = build_field(plain(2.0));
    let b1 = build_field_with_floor(Knobs::passes(2), Some(floor_b1));
    let b2 = build_field(Knobs { slope_floor_uk: Some(KS), ..Knobs::passes(2) });
    let shelf = build_field(Knobs { shelf_min_depth_m: Some(20.0), ..Knobs::passes(2) });
    eprintln!(
        "   builds: delivered x1 {t_del1:.1} s . A1+B2 x1 {t_ab1:.1} . x1.5 {t_ab15:.1} . x2 \
         {t_ab20:.1} . delivered x1.5 {t_del15:.1} . (+ delivered x2, B1, B2, shelf-20)"
    );

    // ══ A — the map of the removed erosion ═══════════════════════════════════
    let cut = |f: &GridF32, k: usize| ((pre.data[k] - f.data[k]) * n2m).max(0.0) as f64;
    // B2's gate, as a declared PROXY: the DELIVERED channel cells whose local slope is already at
    // or below the equilibrium slope, i.e. the cells where the floor would have bound.
    let d_del = c1_drainage_windowed(&del1, None, &on, &ss, DOMAIN_KM);
    let mut floor_binds = vec![false; n];
    let mut is_chan = vec![false; n];
    for k in 0..n {
        if del1.data[k] <= SEA {
            continue;
        }
        let a_km2 = d_del.flow.accumulation.data[k] * cell_km2;
        if a_km2 < RELIEF_V1_A_C_KM2 {
            continue;
        }
        is_chan[k] = true;
        let d = d_del.flow.direction[k];
        if d == DIR_NONE {
            continue;
        }
        let nx = ((k % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
        let ny = ((k / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
        let diag = D8_DX[d as usize] != 0 && D8_DY[d as usize] != 0;
        let dx_m = if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M };
        let s = (del1.data[k] - del1.data[ny * w + nx]).max(0.0) * n2m / dx_m;
        floor_binds[k] = s <= KS / a_km2.sqrt();
    }
    eprintln!(
        "\n-- A . the map of the removed erosion (volume in m of cut, land common to both) --\n  \
         populations: channel cells (A >= A_c) {} . pre-incision pit cells {} . delivered cells \
         BELOW chi {} . cells where B2's floor would bind {}",
        is_chan.iter().filter(|&&b| b).count(),
        in_pit.iter().filter(|&&b| b).count(),
        (0..n).filter(|&k| del1.data[k] > SEA && del1.data[k] < z_chi[k]).count(),
        floor_binds.iter().filter(|&&b| b).count()
    );
    for (label, f) in [("A1+B2", &ab1), ("B1 area cap", &b1), ("B2 alone", &b2)] {
        let (mut v_i, mut v_ii, mut v_iii, mut v_gate_pit, mut v_gate_floor, mut v_gate_none) =
            (0f64, 0f64, 0f64, 0f64, 0f64, 0f64);
        let mut v_neg = 0f64;
        for k in 0..n {
            if pre.data[k] <= SEA || f.data[k] <= SEA || del1.data[k] <= SEA {
                continue;
            }
            let d = cut(&del1, k) - cut(f, k);
            if d <= 0.0 {
                v_neg += -d;
                continue;
            }
            // the ROUND's partition, by position relative to chi
            if is_chan[k] && del1.data[k] < z_chi[k] {
                v_i += d;
            } else if in_pit[k] {
                v_ii += d;
            } else {
                v_iii += d;
            }
            // the GATE partition: where the closure could actually have acted
            if in_pit[k] {
                v_gate_pit += d;
            } else if floor_binds[k] {
                v_gate_floor += d;
            } else {
                v_gate_none += d;
            }
        }
        let tot = (v_i + v_ii + v_iii).max(1.0);
        let gt = (v_gate_pit + v_gate_floor + v_gate_none).max(1.0);
        eprintln!(
            "   **{label}**: removed {:.3e} m . ROUND's partition (i) below-chi channel **{:.1} \
             %** . (ii) pit **{:.1} %** . (iii) rest **{:.1} %**  ||  GATE partition: pit \
             **{:.1} %** . floor-binds **{:.1} %** . **NO GATE {:.1} %** . (added back elsewhere: \
             {:.3e} m)",
            tot,
            100.0 * v_i / tot,
            100.0 * v_ii / tot,
            100.0 * v_iii / tot,
            100.0 * v_gate_pit / gt,
            100.0 * v_gate_floor / gt,
            100.0 * v_gate_none / gt,
            v_neg
        );
    }

    // ══ B and C — the dial, and the shelf clamp ══════════════════════════════
    let auth = spurs(&land_u16(&pre, &ss), w);
    let tot_cut = |f: &GridF32| -> f64 {
        (0..n).filter(|&k| pre.data[k] > SEA && f.data[k] > SEA).map(|k| cut(f, k)).sum()
    };
    let (c_del1, c_del15, c_del20) = (tot_cut(&del1), tot_cut(&del15), tot_cut(&del20));
    eprintln!(
        "\n   total cut, delivered: x1 {c_del1:.3e} . x1.5 {c_del15:.3e} (**x{:.2}**) . x2 \
         {c_del20:.3e} (**x{:.2}**) -- the dial does move the budget",
        c_del15 / c_del1,
        c_del20 / c_del1
    );

    let mut rows: Vec<(String, Crit, f64, f64, f32, i64, i64, f64)> = Vec::new();
    for (label, f, cost, denom) in [
        ("delivered x1", &del1, t_del1, c_del1),
        ("A1+B2 x1", &ab1, t_ab1, c_del1),
        ("A1+B2 x1.5", &ab15, t_ab15, c_del15),
        ("A1+B2 x2", &ab20, t_ab20, c_del20),
        ("CONTROL delivered x1.5", &del15, t_del15, c_del15),
        ("C shelf clamp 20 m", &shelf, t_del1, c_del1),
    ] {
        let d = c1_drainage_windowed(f, None, &on, &ss, DOMAIN_KM);
        let bf = breach_monotone(f, &d.flow.filled, &d.lake_map, SEA, w, h);
        let c = f95_criteria(f, &bf, &pre, DELIVERED_P50_M, &ss, &on, cell_km2, n2m, w, h);
        let land: Vec<bool> = (0..n).map(|k| f.data[k] > SEA).collect();
        let r8 = aniso(f, &land, 16).r8;
        let cutf = tot_cut(f);
        let sp_br = spurs(&land_u16(&bf, &ss), w) as i64 - auth as i64;
        let sp_er = spurs(&land_u16(f, &ss), w) as i64 - auth as i64;
        eprintln!(
            "   **{label}**: canyons **{}/{} ({:.1} %)** . coast breached **{sp_br:+}** . eroded \
             u16 {sp_er:+} . relief **{:.1} m ({:+.1} % vs oracle)** . lakes {:.2} % . erosion \
             permitted **{:.1} % of delivered x1** / **{:.1} % of delivered at the SAME k_time** \
             . R8 **{r8:.4}** . {cost:.1} s",
            c.klass,
            c.scanned,
            100.0 * c.klass as f32 / c.scanned.max(1) as f32,
            c.p50,
            100.0 * (c.p50 - ORACLE_P50_M) / ORACLE_P50_M,
            c.lake_pct,
            100.0 * cutf / c_del1,
            100.0 * cutf / denom
        );
        rows.push((label.to_string(), c, cost, cutf, r8, sp_br, sp_er, denom));
    }

    eprintln!("\n-- the decision --");
    for (label, c, _cost, cutf, r8, sp_br, _se, denom) in &rows {
        let dp = 100.0 * (c.p50 - ORACLE_P50_M) / ORACLE_P50_M;
        eprintln!(
            "   {label}: relief {} . canyons {:.1} % . coast {sp_br:+} . R8 {r8:.4} . erosion \
             permitted (same k_time) {:.1} %",
            if dp.abs() <= 10.0 { "INSIDE +-10 %" } else { "**OUTSIDE +-10 %**" },
            100.0 * c.klass as f32 / c.scanned.max(1) as f32,
            100.0 * cutf / denom
        );
    }
    eprintln!(
        "\n==========  end Finding 98 . total {:.1} s  ==========\n",
        t0.elapsed().as_secs_f64()
    );
}
