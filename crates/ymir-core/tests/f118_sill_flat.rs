//! ADR 0001 Finding 118 — **why the law does not bite at the sill: who flattened it, what its
//! receiver is, and what locality costs.**
//!
//! Findings 114–117 left one unknown. The 18 family-1 lakes have a median drainable of 200 m
//! (Finding 117), the retreat wave crosses their outlet paths in 0.2 pass (Finding 115), and
//! nothing moves because `S = 0` at the sill (Finding 114). `E = K·A^m·S^n` is zero at the one
//! cell that would lower the lake. This bench attributes the zero. **No mechanism.**
//!
//! ## Rule 12 / 18, declared
//!
//! The 18 lakes and their sill cells are the closure-ON breached inventory's (Finding 117). Every
//! other stage — pre-incision, eroded OFF, eroded ON, breached OFF — is read **at the same cell**,
//! with that stage's OWN flow direction. On the pre-incision field most of these depressions do not
//! exist yet (Finding 117: the incision creates 91 % of the deep pits), so "S at the sill on S3" is
//! S at a hillslope cell; the bench says, per lake and per stage, whether the depression exists
//! there (`filled > raw` at the lake floor).
//!
//! ## Blocks
//! * **A** — S at the sill on each stage ⇒ who set it to zero: the terrain, the one-pass landing
//!   (Finding 88-D1, 97 % in pass 1; Finding 114, `h ← (h_o + f·h_r)/(1+f)`), or the breach.
//! * **B** — the sill's D8 receiver on the exported stage: a LAKE cell (water — graded onto a
//!   downstream lake) or LAND (the scheme graded the chain onto itself).
//! * **C** — `S_k = (h_sill − h_d@k)/(path k)` for k ∈ {1, 10, 50, 200}; the k at which S first
//!   exceeds 0 is the distance to the knickpoint; and the one-pass lowering the shipped update
//!   WOULD apply if it read S over k cells: `(h − h_k)·f_k/(1+f_k)`, `f_k = k·dt·A^m/(k·dx)`.
//!   **Not a remedy — the price of locality, measured.**
//! * **D** — Δ(sill height) eroded → breached, OFF and ON. Does the breach lower, flatten, or
//!   not touch the sill?
//! * **E** — cited from Finding 117: the eight lakes held by a family-2 basin, drainable to the
//!   basin's level.
//!
//! Run: cargo test -p ymir-core --release --test f118_sill_flat -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, build_field_seed, pct, sorted};
use std::time::Instant;
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::erosion::stream_power::{RELIEF_V1_A_C_KM2, StreamPowerConfig};
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, SegmentKind, c1_drainage_windowed,
};
use ymir_core::tectonics_c1::hd_assembly::assemble_hd_drainage;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowResult, breach_monotone};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const GEO_RATIO: f32 = 7.5;
const S_EQ: f32 = 0.024;
const KS: [usize; 4] = [1, 10, 50, 200];

/// Walk `k` cells down `flow.direction` from `start`; returns (cell reached, path length m).
fn walk(flow: &FlowResult, start: usize, k: usize, w: usize, h: usize) -> (usize, f32) {
    let (mut c, mut len) = (start, 0f32);
    for _ in 0..k {
        let d = flow.direction[c];
        if d == DIR_NONE {
            break;
        }
        let nx = ((c % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
        let ny = ((c / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
        let diag = D8_DX[d as usize] != 0 && D8_DY[d as usize] != 0;
        len += if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M };
        c = ny * w + nx;
    }
    (c, len)
}

struct Stage<'a> {
    name: &'static str,
    g: &'a GridF32,
    flow: &'a FlowResult,
}

#[test]
#[ignore]
fn f118_sill_flat() {
    let ss = SteinSteinParams::default();
    let m = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let sp = StreamPowerConfig::relief_v3(cell_km2, ss.depth_scale_m as f32);
    let kdt = sp.k * sp.dt;
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 118 . why S = 0 at the sill  ==========");

    // ── the five stages ─────────────────────────────────────────────────────
    let pre = build_field_seed(Knobs::no_incision(), PSEED);
    let ero_off = build_field_seed(Knobs::passes(2), PSEED);
    let ero_on = build_field_seed(Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
    let (w, h) = (pre.width, pre.height);
    let d_pre = c1_drainage_windowed(&pre, None, &dcfg, &ss, DOMAIN_KM);
    let d_eoff = c1_drainage_windowed(&ero_off, None, &dcfg, &ss, DOMAIN_KM);
    let d_eon = c1_drainage_windowed(&ero_on, None, &dcfg, &ss, DOMAIN_KM);
    let bre_off = breach_monotone(&ero_off, &d_eoff.flow.filled, &d_eoff.lake_map, SEA, w, h);
    let bre_on = breach_monotone(&ero_on, &d_eon.flow.filled, &d_eon.lake_map, SEA, w, h);
    let d_boff = c1_drainage_windowed(&bre_off, None, &dcfg, &ss, DOMAIN_KM);
    let d_bon = c1_drainage_windowed(&bre_on, None, &dcfg, &ss, DOMAIN_KM);
    let stages = [
        Stage { name: "pre-incision", g: &pre, flow: &d_pre.flow },
        Stage { name: "eroded OFF", g: &ero_off, flow: &d_eoff.flow },
        Stage { name: "eroded ON", g: &ero_on, flow: &d_eon.flow },
        Stage { name: "breached OFF", g: &bre_off, flow: &d_boff.flow },
        Stage { name: "breached ON", g: &bre_on, flow: &d_bon.flow },
    ];

    // ── the inventory: closure-ON breached, as Finding 117 ──────────────────
    let climate = c1_climate_placed(&bre_on, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
    let dclim = DrainageClimate {
        precip_internal: &climate.precipitation,
        temperature: &climate.temperature,
    };
    let pre2 = c1_drainage_windowed(&ero_on, None, &dcfg, &ss, DOMAIN_KM);
    let dr = assemble_hd_drainage(
        &bre_on,
        &dclim,
        Some(pre2),
        &dcfg,
        &ss,
        DOMAIN_KM,
        GEO_RATIO,
        None,
        false,
    )
    .drainage;
    let mut fam1: Vec<_> =
        dr.lakes.iter().filter(|l| l.area_km2 >= 1.0 && l.base.id < 1_000_001).collect();
    fam1.sort_by(|a, b| b.area_km2.total_cmp(&a.area_km2));

    // ── A + D · the sill through the stages ─────────────────────────────────
    eprintln!(
        "\n   ── A/D · the sill cell across the five stages: height (m) · S at d@1 · depression \
         exists at the floor? ──\n   {:>8}  {:>26} {:>26} {:>26} {:>26} {:>26}",
        "area km²", "pre-incision", "eroded OFF", "eroded ON", "breached OFF", "breached ON"
    );
    let (mut s_pre_pos, mut s_eon_zero, mut bre_delta_zero_on, mut dep_pre) = (0, 0, 0, 0);
    let mut water_d1 = 0usize;
    let mut sk_rows: Vec<(f32, [f32; 4], [f32; 4], usize)> = Vec::new();
    for l in &fam1 {
        let (ox, oy) = (l.base.outlet.0 as usize, l.base.outlet.1 as usize);
        let k = oy * w + ox;
        if k >= w * h {
            continue;
        }
        // the lake floor, from the inventory's own cells
        let floor = (0..w * h)
            .filter(|&c| dr.lake_map[c] == l.base.id)
            .min_by(|&a, &b| bre_on.data[a].total_cmp(&bre_on.data[b]))
            .unwrap_or(k);
        let mut line = format!("   {:>8.2}", l.area_km2);
        let mut h_eon = 0f32;
        let mut h_bon = 0f32;
        for (i, st) in stages.iter().enumerate() {
            let (d1, len1) = walk(st.flow, k, 1, w, h);
            let s1 = if len1 > 0.0 { ((m(st.g, k) - m(st.g, d1)).max(0.0)) / len1 } else { 0.0 };
            let dep = st.flow.filled.data[floor] > st.g.data[floor];
            line +=
                &format!("  {:>9.2} S{:>7.4} {}", m(st.g, k), s1, if dep { "dep" } else { " — " });
            match i {
                0 => {
                    if s1 > 1e-4 {
                        s_pre_pos += 1;
                    }
                    if dep {
                        dep_pre += 1;
                    }
                }
                2 => {
                    h_eon = m(st.g, k);
                    if s1 <= 1e-4 {
                        s_eon_zero += 1;
                    }
                }
                4 => {
                    h_bon = m(st.g, k);
                }
                _ => {}
            }
        }
        if (h_bon - h_eon).abs() < 0.005 {
            bre_delta_zero_on += 1;
        }
        eprintln!("{line}");

        // ── B · what is d@1 on the exported stage? ──
        let (d1, _) = walk(&d_bon.flow, k, 1, w, h);
        if dr.lake_map[d1] != 0 {
            water_d1 += 1;
        }

        // ── C · S_k and the one-pass lowering the law WOULD apply at reach k, eroded ON ──
        let a_km2 = (d_eon.flow.accumulation.data[k] * cell_km2).max(RELIEF_V1_A_C_KM2);
        let mut sk = [0f32; 4];
        let mut low = [0f32; 4];
        let mut first_k = 0usize;
        for (i, &kk) in KS.iter().enumerate() {
            let (dk, len) = walk(&d_eon.flow, k, kk, w, h);
            let drop = (m(&ero_on, k) - m(&ero_on, dk)).max(0.0);
            sk[i] = if len > 0.0 { drop / len } else { 0.0 };
            let fk = kdt * a_km2.powf(sp.m) / len.max(CELL_M);
            low[i] = drop * fk / (1.0 + fk);
            if first_k == 0 && sk[i] > 1e-4 {
                first_k = kk;
            }
        }
        sk_rows.push((l.area_km2, sk, low, first_k));
    }
    let n = fam1.len();
    eprintln!(
        "\n   ⇒ A · S > 0 at the (future) sill on PRE-INCISION: **{s_pre_pos} of {n}** · depression \
         already exists there: **{dep_pre} of {n}** · S = 0 at the sill on ERODED ON: **{s_eon_zero} \
         of {n}**\n   ⇒ D · breach leaves the sill height unchanged (|Δ| < 5 mm, ON): **{bre_delta_zero_on} \
         of {n}**\n   ⇒ B · d@1 on the exported stage is a LAKE cell (water): **{water_d1} of {n}** · \
         land: **{} of {n}**",
        n - water_d1
    );

    // ── C · the table ───────────────────────────────────────────────────────
    eprintln!(
        "\n   ── C · S over k cells (eroded ON) and the one-pass lowering the update would apply \
         at that reach ──\n   {:>8} {:>8} {:>8} {:>8} {:>8}   {:>8} {:>8} {:>8} {:>8}  {:>7}",
        "area", "S_1", "S_10", "S_50", "S_200", "low@1", "low@10", "low@50", "low@200", "knick@"
    );
    let (mut s10_pos, mut s1_zero) = (0, 0);
    let mut firsts = Vec::new();
    for (a, sk, low, fk) in &sk_rows {
        if sk[1] > 1e-4 {
            s10_pos += 1;
        }
        if sk[0] <= 1e-4 {
            s1_zero += 1;
        }
        if *fk > 0 {
            firsts.push(*fk as f32);
        }
        eprintln!(
            "   {a:>8.2} {:>8.4} {:>8.4} {:>8.4} {:>8.4}   {:>8.1} {:>8.1} {:>8.1} {:>8.1}  {:>7}",
            sk[0],
            sk[1],
            sk[2],
            sk[3],
            low[0],
            low[1],
            low[2],
            low[3],
            if *fk > 0 { fk.to_string() } else { "none".into() }
        );
    }
    let sf = sorted(firsts.clone());
    eprintln!(
        "   ⇒ S_1 = 0: **{s1_zero} of {n}** · S_10 > 0: **{s10_pos} of {n}** · knickpoint distance \
         (first k with S > 0) p10/p50/p90: **{:.0} / {:.0} / {:.0} cells**",
        if sf.is_empty() { 0.0 } else { pct(&sf, 0.10) },
        if sf.is_empty() { 0.0 } else { pct(&sf, 0.50) },
        if sf.is_empty() { 0.0 } else { pct(&sf, 0.90) }
    );

    // ── controls ────────────────────────────────────────────────────────────
    if let Some(i) = (0..dr.rivers.segments.len())
        .filter(|&i| dr.segment_kind[i] == SegmentKind::Spillway)
        .max_by(|&a, &b| dr.segment_discharge_m3s[a].total_cmp(&dr.segment_discharge_m3s[b]))
    {
        let pts = &dr.rivers.segments[i].points;
        let mut zero = 0usize;
        for p in pts.windows(2) {
            let (a, b) =
                (p[0].1 as usize * w + p[0].0 as usize, p[1].1 as usize * w + p[1].0 as usize);
            if (m(&bre_on, a) - m(&bre_on, b)).abs() < 1e-3 {
                zero += 1;
            }
        }
        eprintln!(
            "\n   CONTROL · Finding 87-A spillway ({:.2} m³/s, {} cells): steps with S = 0 **{zero} of \
             {}** ⇒ {}",
            dr.segment_discharge_m3s[i],
            pts.len(),
            pts.len().saturating_sub(1),
            if zero * 4 < pts.len() {
                "graded, S > 0 almost everywhere"
            } else {
                "**flat runs — check**"
            }
        );
    }
    if let Some(l2) = dr.lakes.iter().find(|l| l.base.id >= 1_000_001 && l.area_km2 >= 1.0) {
        let k = l2.base.outlet.1 as usize * w + l2.base.outlet.0 as usize;
        let mut smax = 0f32;
        for kk in KS {
            let (dk, len) = walk(&d_bon.flow, k, kk, w, h);
            if len > 0.0 {
                smax = smax.max((m(&bre_on, k) - m(&bre_on, dk)).max(0.0) / len);
            }
        }
        eprintln!(
            "   CONTROL · a family-2 sill (basin {}, {:.0} km²): max S over 1..200 cells **{smax:.5}** \
             ⇒ {}",
            l2.base.id,
            l2.area_km2,
            if smax < 1e-4 { "S = 0 by construction, as Finding 114" } else { "**not flat**" }
        );
    }
    eprintln!(
        "\n   E · cited from Finding 117 (no new build): the eight family-1 lakes held by a family-2 \
         basin have drainable-to-basin-level 47.5 / 54.3 / 82.5 / 322.4 / 341.1 / 873.4 / 917.1 / \
         1 206.7 m — all > 0, five above 300 m."
    );
    eprintln!("\n==========  end Finding 118 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
