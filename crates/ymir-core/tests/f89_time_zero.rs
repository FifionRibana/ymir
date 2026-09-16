//! ADR 0001 Finding 89 — the time-zero baseline (C) and the cost of one pass (D).
//!
//! No production change. Block A is a reading of Finding 44, block B writes criteria; this file
//! measures the six criteria on the DELIVERED continent and prices an extra incision pass.
//!
//! ⚠️ The export in `exports/seed10481999410520546993_8192.ymir/` is dated 5 September and predates
//! Findings 83, 86 and 88 (13 lakes, 2 209 km², no `Unresolved`, no base-level bound). It is NOT
//! the delivered product and is not read here.
//!
//! Rule 12: every walk declares its field. Run:
//!   cargo test -p ymir-core --release --test f89_time_zero -- --ignored --nocapture

mod common;

use common::{
    CELL_KM, Knobs, SEA, build_field, declared_flow, land_u16, majority, pct, sorted, to_mask,
};
use std::collections::{HashMap, HashSet};
use std::time::Instant;
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::erosion::stream_power::{RELIEF_V1_A_C_KM2, SHIPPED_K_TIME};
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, LakeType, SegmentKind, SegmentRow, UnresolvedReason,
    apply_lake_water_balance, below_sea_basin_lakes_infil, c1_drainage_windowed,
    clip_rivers_to_lakes, resolve_exorheic_without_outlet, runoff_accumulation, runoff_km2_to_m3s,
};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::coast_metrics::{NECK_KM, coast_spurs};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowConfig, breach_monotone, compute_flow};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
/// Finding 88-D1's control point, by its FLOOR CELL and not by an id (Finding 81).
const DEEP_FLOOR: (usize, usize) = (3487, 5902);

fn median(v: &[f32]) -> f32 {
    pct(&sorted(v.to_vec()), 0.50)
}

fn slope_deg(field: &GridF32, k: usize, n2m: f32, w: usize, h: usize) -> f32 {
    let (x, y) = ((k % w) as i32, (k / w) as i32);
    let mut best = 0.0f32;
    for d in 0..8 {
        let nx = (x + D8_DX[d]).rem_euclid(w as i32) as usize;
        let ny = (y + D8_DY[d]).rem_euclid(h as i32) as usize;
        let diag = D8_DX[d] != 0 && D8_DY[d] != 0;
        let dist = if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M };
        best = best.max((field.data[k] - field.data[ny * w + nx]).abs() * n2m / dist);
    }
    best.atan().to_degrees()
}

#[test]
#[ignore]
fn f89_time_zero() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 89 · C (the baseline) · D (the price of a pass)  ==========");

    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;

    // ── D · the builds, timed ───────────────────────────────────────────────
    let t0 = Instant::now();
    let raw = build_field(Knobs::passes(2));
    let t_p2 = t0.elapsed().as_secs_f64();
    let (w, h) = (raw.width, raw.height);
    let n = w * h;
    let t1 = Instant::now();
    let raw3 = build_field(Knobs::passes(3));
    let t_p3 = t1.elapsed().as_secs_f64();
    let t2 = Instant::now();
    let raw4 = build_field(Knobs::passes(4));
    let t_p4 = t2.elapsed().as_secs_f64();
    let t3 = Instant::now();
    let pre = build_field(Knobs::no_incision());
    let t_pre = t3.elapsed().as_secs_f64();
    // ADR Finding 89-D control — the A/B of the Finding 83 bound at the same pass count, because
    // the raw shipped floor came out at exactly +0.50 m and that is `RELIEF_V3_BASE_LEVEL_M`.
    let raw_nb = build_field(Knobs { base_level_off: true, ..Knobs::passes(2) });

    let it1 = build_field(Knobs::passes(1));
    let pre0 = c1_drainage_windowed(&raw, None, &on, &ss, DOMAIN_KM);
    let field = breach_monotone(&raw, &pre0.flow.filled, &pre0.lake_map, SEA, w, h);
    let wc = water_class(&field, SEA);
    let land: Vec<bool> = (0..n).map(|k| field.data[k] > SEA).collect();
    let land_km2 = land.iter().filter(|&&b| b).count() as f32 * cell_km2;
    let m = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);

    eprintln!(
        "\n── D · wall time, whole pipeline, one core count, release ──\n   \
         iterations 2 (shipped) **{t_p2:.1} s** · 3 **{t_p3:.1} s** · 4 **{t_p4:.1} s** · \
         no incision {t_pre:.1} s"
    );
    let marg23 = t_p3 - t_p2;
    let marg34 = t_p4 - t_p3;
    eprintln!(
        "   marginal pass 2→3 **{marg23:+.1} s** · 3→4 **{marg34:+.1} s** · the incision stage \
         alone (shipped 2 passes vs none) **{:+.1} s** ⇒ ~{:.1} s/pass",
        t_p2 - t_pre,
        (t_p2 - t_pre) / 2.0
    );
    eprintln!(
        "   ⚠️ the constant-cost hypothesis: 2→3 against 3→4 differ by {:+.1} % — {}",
        100.0 * (marg34 - marg23) as f64 / marg23.abs().max(1e-6) as f64,
        if (marg34 - marg23).abs() / marg23.abs().max(1e-6) < 0.15 {
            "**constant within 15 %**, so the extrapolation below is linear"
        } else {
            "**NOT constant**, the extrapolation below is an upper/lower bound only"
        }
    );
    // ══ the control the round demanded — and it FIRED, so it is reported as a finding ══
    //
    // "if the floor is not at −8.67 m before pass 3, it is not the same world". It is not: the
    // ERODED field has that cell at +0.50 m. −8.67 m is the value AFTER `breach_monotone`, and
    // Finding 88-D1 compared passes 0 and 1 on the eroded field against pass 2 on the BREACHED
    // one. Two populations in one row. Both are reported here, separately.
    let kf = DEEP_FLOOR.1 * w + DEEP_FLOOR.0;
    let br = |g: &GridF32| -> GridF32 {
        let d = c1_drainage_windowed(g, None, &on, &ss, DOMAIN_KM);
        breach_monotone(g, &d.flow.filled, &d.lake_map, SEA, w, h)
    };
    let f3 = br(&raw3);
    let f4 = br(&raw4);
    eprintln!(
        "\n   D control at ({},{}) — **the ERODED field, one population**: pass 0 **{:.2} m** → \
         pass 1 **{:.2} m** → pass 2 (shipped) **{:.2} m** → pass 3 **{:.2} m** → pass 4 \
         **{:.2} m**",
        DEEP_FLOOR.0,
        DEEP_FLOOR.1,
        m(&pre, kf),
        m(&it1, kf),
        m(&raw, kf),
        m(&raw3, kf),
        m(&raw4, kf)
    );
    eprintln!(
        "   D control — **the BREACHED field (what is delivered)**: pass 2 **{:.2} m** \
         (Finding 87-B / 88-D1's −8.67) · pass 3 **{:.2} m** · pass 4 **{:.2} m**",
        m(&field, kf),
        m(&f3, kf),
        m(&f4, kf)
    );
    eprintln!(
        "   ⚠️ **the attribution, corrected**: the INCISION takes {:.2} → {:.2} m ({:+.2} m, of \
         which pass 1 carries {:.1} %); `breach_monotone` then takes {:.2} → {:.2} m ({:+.2} m) \
         and it is what puts this cell BELOW SEA LEVEL. With the Finding 83 bound OFF the eroded \
         floor is **{:.2} m** ⇒ the bound {}.",
        m(&pre, kf),
        m(&raw, kf),
        m(&raw, kf) - m(&pre, kf),
        100.0 * (m(&pre, kf) - m(&it1, kf)) / (m(&pre, kf) - m(&raw, kf)).abs().max(1e-6),
        m(&raw, kf),
        m(&field, kf),
        m(&field, kf) - m(&raw, kf),
        m(&raw_nb, kf),
        if (m(&raw, kf) - 0.5).abs() < 0.05 && m(&raw_nb, kf) < 0.45 {
            "**BINDS at this cell, exactly at sea + 0.5 m**"
        } else {
            "does NOT bind at this cell"
        }
    );
    assert!(
        (m(&field, kf) - (-8.67)).abs() < 0.2,
        "D control: the DELIVERED (breached) floor must be −8.67 m, got {:.2}",
        m(&field, kf)
    );

    // per-stage costs for the frequency table
    let ts = Instant::now();
    let _ = compute_flow(&field, &FlowConfig { sea_level: SEA, ..Default::default() });
    let t_flow = ts.elapsed().as_secs_f64();
    let ts = Instant::now();
    let climate = c1_climate_placed(&field, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
    let t_clim = ts.elapsed().as_secs_f64();
    let dclim = DrainageClimate {
        precip_internal: &climate.precipitation,
        temperature: &climate.temperature,
    };
    let ts = Instant::now();
    let mut dr = c1_drainage_windowed(&field, Some(&dclim), &on, &ss, DOMAIN_KM);
    let t_drain = ts.elapsed().as_secs_f64();
    let _ = declared_flow("dr.flow", &dr.flow, &dr.flow);
    eprintln!(
        "   per-stage, once each on the delivered field: `compute_flow` **{t_flow:.1} s** · \
         `c1_climate_placed` **{t_clim:.1} s** · `c1_drainage_windowed` (flow + rivers + lakes) \
         **{t_drain:.1} s**"
    );
    let cfl_steps = 7399.0f64; // Finding 44's own derived count at 8192², A_max 1611 km²
    for (label, per) in
        [("marginal pass (2→3)", marg23 as f64), ("incision only", ((t_p2 - t_pre) / 2.0) as f64)]
    {
        eprintln!(
            "   D extrapolation on {label} = {per:.1} s/pass: 100 passes {:.2} h · 1 000 passes \
             {:.1} h · **Finding 44's derived 7 399 passes {:.1} h = {:.1} days**",
            per * 100.0 / 3600.0,
            per * 1000.0 / 3600.0,
            per * cfl_steps / 3600.0,
            per * cfl_steps / 86400.0
        );
    }

    // ── the chain, to get the delivered lake inventory ──────────────────────
    let pre_d = c1_drainage_windowed(&raw, None, &on, &ss, DOMAIN_KM);
    dr.lakes = pre_d.lakes;
    dr.lake_map = pre_d.lake_map;
    let lakes_in = std::mem::take(&mut dr.lakes);
    dr.lakes = apply_lake_water_balance(
        &field,
        &dr.flow,
        &dclim,
        cell_km2,
        &ss,
        &lakes_in,
        &mut dr.lake_map,
        None,
        w,
        h,
    );
    let racc = runoff_accumulation(&field, &dr.flow, &dclim, cell_km2, None, None, w, h);
    let premerge = dr.lake_map.clone();
    let bs =
        below_sea_basin_lakes_infil(&field, &dclim, &on, &ss, DOMAIN_KM, Some(&premerge), None);
    let t = &bs.termination;
    let mut det_before: HashMap<u32, usize> = HashMap::new();
    for &id in &dr.lake_map {
        if id != 0 && id < 1_000_001 {
            *det_before.entry(id).or_default() += 1;
        }
    }
    let mut det_after: HashMap<u32, usize> = HashMap::new();
    for k in 0..n {
        if bs.lake_map[k] != 0 {
            dr.lake_map[k] = bs.lake_map[k];
        } else if dr.lake_map[k] != 0 && dr.lake_map[k] < 1_000_001 {
            *det_after.entry(dr.lake_map[k]).or_default() += 1;
        }
    }
    let absorbed: HashSet<u32> = det_before
        .keys()
        .copied()
        .filter(|id| det_after.get(id).copied().unwrap_or(0) == 0)
        .collect();
    dr.lakes.retain(|l| l.base.id >= 1_000_001 || !absorbed.contains(&l.base.id));
    dr.lakes.extend(bs.lakes.iter().cloned());
    clip_rivers_to_lakes(&mut dr);
    let inventoried: HashSet<u32> = dr.lakes.iter().map(|l| l.base.id).collect();
    for sw in &bs.spillways {
        let (lx, ly) = *sw.points.last().unwrap_or(&(0, 0));
        dr.push_segment(SegmentRow {
            segment: ymir_core::terrain::flow::RiverSegment {
                points: sw.points.clone(),
                strahler_order: 1,
                avg_flow: 0.0,
                max_flow: 0.0,
                basin_id: 0,
                upstream: vec![],
                downstream: None,
            },
            drainage_km2: sw.drainage_km2,
            navigability: sw.navigability,
            discharge_m3s: sw.discharge_m3s,
            width_m: sw.width_m,
            profile_m: sw.profile_m.clone(),
            catchment_cells: {
                let k = ly as usize * w + lx as usize;
                dr.flow.accumulation.data.get(k).copied().unwrap_or(0.0)
            },
            discharge_profile_m3s: vec![sw.discharge_m3s; sw.points.len()],
            kind: SegmentKind::Spillway,
            source_lake: inventoried.contains(&sw.lake_id).then_some(sw.lake_id),
        });
    }
    let relabelled = resolve_exorheic_without_outlet(&mut dr);

    // ── C1 · the lake fraction, by type and by depth class ──────────────────
    eprintln!(
        "\n── C1 · the lake fraction on the DELIVERED continent (humid, shipped, {} passes of the \
         fixed point, converged {}) ──",
        t.passes_used, !t.not_converged
    );
    let mut cells_of: HashMap<u32, usize> = HashMap::new();
    for k in 0..n {
        if dr.lake_map[k] != 0 {
            *cells_of.entry(dr.lake_map[k]).or_default() += 1;
        }
    }
    let mut by_type: HashMap<String, (usize, f32)> = HashMap::new();
    let mut by_depth: HashMap<&str, (usize, f32, f64)> = HashMap::new();
    let mut water_km2 = 0.0f32;
    let mut water_km3 = 0.0f64;
    for l in &dr.lakes {
        let c = cells_of.get(&l.base.id).copied().unwrap_or(0);
        let a = c as f32 * cell_km2;
        water_km2 += a;
        // volume: sum of (level − h) over the footprint, in km³
        let mut vol = 0.0f64;
        for k in 0..n {
            if dr.lake_map[k] == l.base.id {
                vol += ((l.level_m - m(&field, k)).max(0.0) as f64) * cell_km2 as f64 / 1000.0;
            }
        }
        water_km3 += vol;
        let e = by_type.entry(format!("{:?}", l.lake_type)).or_insert((0, 0.0));
        e.0 += 1;
        e.1 += a;
        let band = if l.depth_m < 5.0 {
            "< 5 m"
        } else if l.depth_m < 30.0 {
            "5–30 m"
        } else if l.depth_m < 100.0 {
            "30–100 m"
        } else {
            "> 100 m"
        };
        let d = by_depth.entry(band).or_insert((0, 0.0, 0.0));
        d.0 += 1;
        d.1 += a;
        d.2 += vol;
    }
    eprintln!(
        "   land {land_km2:.0} km² · water bodies **{}** covering **{water_km2:.0} km²** = \
         **{:.2} % of land** · volume **{water_km3:.1} km³**",
        dr.lakes.len(),
        100.0 * water_km2 / land_km2
    );
    let mut tv: Vec<(&String, &(usize, f32))> = by_type.iter().collect();
    tv.sort_by(|a, b| b.1.1.total_cmp(&a.1.1));
    for (ty, (c, a)) in tv {
        eprintln!(
            "      {c:>3} × {ty:<12} {a:9.1} km² ({:5.1} % of the water)",
            100.0 * a / water_km2
        );
    }
    for band in ["< 5 m", "5–30 m", "30–100 m", "> 100 m"] {
        if let Some((c, a, v)) = by_depth.get(band) {
            eprintln!(
                "      depth {band:<9} {c:>3} bodies · {a:9.1} km² ({:5.1} % of area) · {v:7.1} \
                 km³ ({:5.1} % of volume)",
                100.0 * a / water_km2,
                100.0 * v / water_km3
            );
        }
    }

    // ── C2 · the canyon class, reproduced then extended ──────────────────────
    eprintln!("\n── C2 · the drowned-canyon class (Finding 88-D4's instrument, reproduced) ──");
    let mut rows: Vec<(u32, f32, f32, f32, f64)> = Vec::new();
    let mut scanned = 0usize;
    for l in dr.lakes.iter().filter(|l| l.area_km2 >= 1.0) {
        let cells: Vec<usize> = (0..n).filter(|&k| dr.lake_map[k] == l.base.id).collect();
        if cells.is_empty() {
            continue;
        }
        scanned += 1;
        let cuts: Vec<f32> = cells.iter().map(|&k| m(&pre, k) - m(&field, k)).collect();
        let inside: HashSet<usize> = cells.iter().copied().collect();
        let mut ring: HashSet<usize> = HashSet::new();
        for &k in &cells {
            for d in 0..8 {
                let nx = ((k % w) as i32 + D8_DX[d]).rem_euclid(w as i32) as usize;
                let ny = ((k / w) as i32 + D8_DY[d]).rem_euclid(h as i32) as usize;
                let nk = ny * w + nx;
                if !inside.contains(&nk) {
                    ring.insert(nk);
                }
            }
        }
        let rim: Vec<f32> = ring.iter().map(|&k| slope_deg(&field, k, n2m, w, h)).collect();
        let vol: f64 = cells
            .iter()
            .map(|&k| ((l.level_m - m(&field, k)).max(0.0) as f64) * cell_km2 as f64 / 1000.0)
            .sum();
        rows.push((
            l.base.id,
            cells.len() as f32 * cell_km2,
            median(&cuts),
            if rim.is_empty() { 0.0 } else { median(&rim) },
            vol,
        ));
    }
    let klass: Vec<&(u32, f32, f32, f32, f64)> =
        rows.iter().filter(|r| r.2 > 50.0 && r.3 > 30.0).collect();
    let k_area: f32 = klass.iter().map(|r| r.1).sum();
    let k_vol: f64 = klass.iter().map(|r| r.4).sum();
    let all_area: f32 = rows.iter().map(|r| r.1).sum();
    let all_vol: f64 = rows.iter().map(|r| r.4).sum();
    eprintln!(
        "   of {scanned} bodies ≥ 1 km²: **{}** in the class (median cut > 50 m AND rim p50 > 30°) \
         · they carry **{k_area:.0} km² = {:.1} % of the area** of bodies ≥ 1 km² and \
         **{k_vol:.1} km³ = {:.1} % of their volume**",
        klass.len(),
        100.0 * k_area / all_area.max(1e-6),
        100.0 * k_vol / all_vol.max(1e-9)
    );
    eprintln!("   the cut × rim histogram over all {scanned} bodies:");
    for (clo, chi) in [(0.0, 50.0), (50.0, 200.0), (200.0, 500.0), (500.0, 1e9)] {
        let mut line = String::new();
        for (rlo, rhi) in [(0.0, 15.0), (15.0, 30.0), (30.0, 45.0), (45.0, 90.0)] {
            let c =
                rows.iter().filter(|r| r.2 >= clo && r.2 < chi && r.3 >= rlo && r.3 < rhi).count();
            line.push_str(&format!("{c:>5}"));
        }
        eprintln!(
            "      cut {:>5.0}–{:<5} m | rim 0–15° 15–30° 30–45° 45–90° :{line}",
            clo,
            if chi > 1e8 { "+".to_string() } else { format!("{chi:.0}") }
        );
    }

    // ── C3 · integration, and the `to_nothing` cell named ───────────────────
    eprintln!("\n── C3 · drainage integration of the bodies ≥ 1 km² ──");
    let mut dest: HashMap<&str, (usize, f64)> = HashMap::new();
    for l in dr.lakes.iter().filter(|l| l.area_km2 >= 1.0) {
        let cells: Vec<usize> = (0..n).filter(|&k| dr.lake_map[k] == l.base.id).collect();
        let q = cells.iter().map(|&k| racc[k]).fold(0.0f32, f32::max);
        let key = match l.lake_type {
            LakeType::Exorheic => "Exorheic",
            LakeType::Endorheic => "Endorheic",
            LakeType::Unresolved => "Unresolved",
            LakeType::CraterAcidic | LakeType::CraterNeutral => "Crater",
        };
        let e = dest.entry(key).or_insert((0, 0.0));
        e.0 += 1;
        e.1 += runoff_km2_to_m3s(q) as f64;
    }
    let mut dv: Vec<(&&str, &(usize, f64))> = dest.iter().collect();
    dv.sort_by(|a, b| b.1.0.cmp(&a.1.0));
    for (kk, (c, q)) in dv {
        eprintln!("      {c:>3} bodies · local inflow {q:7.1} m³/s · {kk}");
    }
    let unres: Vec<(u32, Option<UnresolvedReason>)> = dr
        .lakes
        .iter()
        .filter(|l| l.lake_type == LakeType::Unresolved)
        .map(|l| (l.base.id, l.unresolved_reason))
        .collect();
    eprintln!(
        "      end-of-chain relabels {} · Unresolved total {} {:?} · `to_nothing` {} · \
         `to_own_body` {} · q_own_body {:.1} · q_lake_unrouted {:.1}",
        relabelled.len(),
        unres.len(),
        unres,
        t.to_nothing,
        t.to_own_body,
        t.q_own_body_m3s,
        t.q_lake_unrouted_m3s
    );
    // name the to_nothing terminus, cell by cell
    for sw in &bs.spillways {
        let &(lx, ly) = sw.points.last().unwrap();
        let kk = ly as usize * w + lx as usize;
        if wc[kk] != 1 && bs.lake_map[kk] == 0 && premerge[kk] == 0 {
            eprintln!(
                "      **`to_nothing` named**: basin {} spillway ends at ({lx},{ly}) — height \
                 {:.3} m · water_class **{}** · bs.lake_map {} · detected lake_map {} · \
                 dr.lake_map {} · D8 dir {} · discharge {:.1} m³/s",
                sw.lake_id,
                m(&field, kk),
                wc[kk],
                bs.lake_map[kk],
                premerge[kk],
                dr.lake_map[kk],
                if dr.flow.direction[kk] == DIR_NONE {
                    "NONE".to_string()
                } else {
                    dr.flow.direction[kk].to_string()
                },
                sw.discharge_m3s
            );
        }
    }

    // ── C4 · relief, PAIRED ─────────────────────────────────────────────────
    eprintln!("\n── C4 · relief, paired on the cells that are land in BOTH fields ──");
    let common: Vec<usize> = (0..n).filter(|&k| field.data[k] > SEA && pre.data[k] > SEA).collect();
    let del: Vec<f32> = common.iter().map(|&k| m(&field, k)).collect();
    let pr: Vec<f32> = common.iter().map(|&k| m(&pre, k)).collect();
    let cut: Vec<f32> = common.iter().map(|&k| m(&pre, k) - m(&field, k)).collect();
    let (sd, sp) = (sorted(del.clone()), sorted(pr.clone()));
    let work: f64 = cut.iter().map(|&c| c as f64).sum::<f64>() * cell_km2 as f64;
    let chan: f64 = common
        .iter()
        .zip(cut.iter())
        .filter(|(k, _)| dr.flow.accumulation.data[**k] * cell_km2 >= RELIEF_V1_A_C_KM2)
        .map(|(_, &c)| c as f64)
        .sum::<f64>()
        * cell_km2 as f64;
    eprintln!(
        "   paired n = {} · delivered p10/p50/p90 **{:.1} / {:.1} / {:.1} m** · pre-incision \
         {:.1} / {:.1} / {:.1} m · median of the per-cell cuts **{:.1} m** · mean cut {:.1} m",
        common.len(),
        pct(&sd, 0.10),
        pct(&sd, 0.50),
        pct(&sd, 0.90),
        pct(&sp, 0.10),
        pct(&sp, 0.50),
        pct(&sp, 0.90),
        median(&cut),
        cut.iter().map(|&c| c as f64).sum::<f64>() / cut.len() as f64
    );
    eprintln!(
        "   erosion work {:.3e} m·km² · of which on cells with A ≥ {RELIEF_V1_A_C_KM2} km²: \
         **{:.1} %** (the channel share)",
        work,
        100.0 * chan / work.abs().max(1e-9)
    );
    let unpaired_d =
        sorted((0..n).filter(|&k| land[k]).map(|k| m(&field, k)).collect::<Vec<f32>>());
    eprintln!(
        "   ⚠️ unpaired delivered p50 {:.1} m against paired {:.1} m — Findings 63–64, stated so \
         the tolerance below is read on the paired figure",
        pct(&unpaired_d, 0.50),
        pct(&sd, 0.50)
    );

    // ── C5 · the coast, four resolutions ────────────────────────────────────
    eprintln!("\n── C5 · Δ(spurs ≥ 2 cells) delivered vs pre-incision, four resolutions ──");
    let l_del = land_u16(&field, &ss);
    let l_pre = land_u16(&pre, &ss);
    // ⚠️ TWO thresholds, because they are two criteria and the first run conflated them.
    // Findings 80 and 83 measured the acceptance leg at a **1 km** minimum spur length
    // (`coast_spurs(.., CELL_KM, 1.0, NECK_KM)`); "≥ 2 cells" is 98 m at 8192² and counts a
    // different population entirely. Both are reported; the CRITERION is the 1 km column.
    for (label, f) in [("8192", 1usize), ("4096", 2), ("2048", 4), ("1024", 8)] {
        let (md, wd) = majority(&l_del, w, f);
        let (mp, wp) = majority(&l_pre, w, f);
        let ck = CELL_KM * f as f32;
        let pd = marching_squares(&to_mask(&md, wd), 0.5);
        let pp = marching_squares(&to_mask(&mp, wp), 0.5);
        let (k_d, _) = coast_spurs(&pd, ck, 1.0, NECK_KM);
        let (k_p, _) = coast_spurs(&pp, ck, 1.0, NECK_KM);
        let (c_d, _) = coast_spurs(&pd, ck, 2.0 * ck, NECK_KM);
        let (c_p, _) = coast_spurs(&pp, ck, 2.0 * ck, NECK_KM);
        eprintln!(
            "      {label:>5} (cell {:6.1} m): **≥ 1 km (THE CRITERION): delivered {} · \
             pre-incision {} · Δ {:+}** | ≥ 2 cells: {} · {} · Δ {:+}",
            ck * 1000.0,
            k_d.len(),
            k_p.len(),
            k_d.len() as i64 - k_p.len() as i64,
            c_d.len(),
            c_p.len(),
            c_d.len() as i64 - c_p.len() as i64
        );
    }

    // ── C6 · the implicit duration of time zero ─────────────────────────────
    eprintln!(
        "\n── C6 · how many years the delivered continent represents, k_time = {SHIPPED_K_TIME} ──"
    );
    // ⚠️ UNITS. `SHIPPED_K_TIME = 9000` is in Ymir's km² convention; every published `K` is in
    // m². The bridge is `K_ours = 10^(3·2m) · K_lit` — ×1000 at m = 0.5, ×251 at the m = 0.4 the
    // Stock & Montgomery table was fitted with. The first run of this block divided 9000 by a
    // literature `K` directly and reported 4.5e8 yr where the anchoring report says 4.5e5: a
    // factor 10³, and the third unit error in three rounds caught by reproducing a known number.
    eprintln!(
        "   conversion: K_ours = 10^(6m)·K_lit ⇒ **×1000 at m = 0.5**, ×251 at m = 0.4 \
         (the exponent each table was fitted with)"
    );
    for (src, k_lit, mult, note) in [
        (
            "Whipple & Tucker 1999 Table 2, n = 1, m/n = 0.50",
            2.00e-5f64,
            1000.0f64,
            "same (m, n) as Ymir; the authors disown the value as a measurement",
        ),
        (
            "Harel 2016 global mean denudation, back-derived",
            7.5e-7,
            1000.0,
            "from ⟨E⟩ = 242 mm/ky at an ASSUMED A = 1085 km², S = 0.01 — linear in S",
        ),
        (
            "Stock & Montgomery hard rock, LOW",
            1.0e-7,
            251.0,
            "m = 0.4 fit read at Ymir's m = 0.5 — ORDER OF MAGNITUDE ONLY",
        ),
        ("Stock & Montgomery hard rock, HIGH", 1.0e-6, 251.0, "idem"),
    ] {
        let k_ours = k_lit * mult;
        eprintln!(
            "      K_lit {k_lit:.2e} ×{mult:.0} ⇒ K_ours {k_ours:.3e} ⇒ T = k_time/K_ours = \
             **{:.3e} yr = {:.1} Myr** — {src}\n           ({note})",
            SHIPPED_K_TIME as f64 / k_ours,
            SHIPPED_K_TIME as f64 / k_ours / 1e6
        );
    }

    // ══ C1 / C2 / C3 on the ARID bed ═══════════════════════════════════════
    //
    // C4, C5 and C6 are properties of the FIELD and of arithmetic, and the field is
    // climate-independent (Finding 60: the erosion never sees the climate), so they are measured
    // once. C1, C2 and C3 read the LAKE INVENTORY, which is what the bed changes.
    eprintln!("\n── C1 / C2 / C3 · the same three on the ARID bed (lat 25, span 10) ──");
    {
        let clim_a =
            c1_climate_placed(&field, &ss, 25.0, 10.0, &PrecipParams::default(), DOMAIN_KM);
        let dca = DrainageClimate {
            precip_internal: &clim_a.precipitation,
            temperature: &clim_a.temperature,
        };
        let mut da = c1_drainage_windowed(&field, Some(&dca), &on, &ss, DOMAIN_KM);
        let pd_a = c1_drainage_windowed(&raw, None, &on, &ss, DOMAIN_KM);
        da.lakes = pd_a.lakes;
        da.lake_map = pd_a.lake_map;
        let li = std::mem::take(&mut da.lakes);
        da.lakes = apply_lake_water_balance(
            &field,
            &da.flow,
            &dca,
            cell_km2,
            &ss,
            &li,
            &mut da.lake_map,
            None,
            w,
            h,
        );
        let racc_a = runoff_accumulation(&field, &da.flow, &dca, cell_km2, None, None, w, h);
        let pm = da.lake_map.clone();
        let bsa = below_sea_basin_lakes_infil(&field, &dca, &on, &ss, DOMAIN_KM, Some(&pm), None);
        let mut d_before: HashMap<u32, usize> = HashMap::new();
        for &id in &da.lake_map {
            if id != 0 && id < 1_000_001 {
                *d_before.entry(id).or_default() += 1;
            }
        }
        let mut d_after: HashMap<u32, usize> = HashMap::new();
        for k in 0..n {
            if bsa.lake_map[k] != 0 {
                da.lake_map[k] = bsa.lake_map[k];
            } else if da.lake_map[k] != 0 && da.lake_map[k] < 1_000_001 {
                *d_after.entry(da.lake_map[k]).or_default() += 1;
            }
        }
        let abs_a: HashSet<u32> = d_before
            .keys()
            .copied()
            .filter(|id| d_after.get(id).copied().unwrap_or(0) == 0)
            .collect();
        da.lakes.retain(|l| l.base.id >= 1_000_001 || !abs_a.contains(&l.base.id));
        da.lakes.extend(bsa.lakes.iter().cloned());
        clip_rivers_to_lakes(&mut da);
        let rel_a = resolve_exorheic_without_outlet(&mut da);

        let mut cells_a: HashMap<u32, usize> = HashMap::new();
        for k in 0..n {
            if da.lake_map[k] != 0 {
                *cells_a.entry(da.lake_map[k]).or_default() += 1;
            }
        }
        let mut wa = 0.0f32;
        let mut va = 0.0f64;
        let mut ty_a: HashMap<String, (usize, f32)> = HashMap::new();
        let mut deep_a = (0usize, 0.0f32, 0.0f64);
        for l in &da.lakes {
            let a = cells_a.get(&l.base.id).copied().unwrap_or(0) as f32 * cell_km2;
            wa += a;
            let mut vol = 0.0f64;
            for k in 0..n {
                if da.lake_map[k] == l.base.id {
                    vol += ((l.level_m - m(&field, k)).max(0.0) as f64) * cell_km2 as f64 / 1000.0;
                }
            }
            va += vol;
            let e = ty_a.entry(format!("{:?}", l.lake_type)).or_insert((0, 0.0));
            e.0 += 1;
            e.1 += a;
            if l.depth_m >= 100.0 {
                deep_a.0 += 1;
                deep_a.1 += a;
                deep_a.2 += vol;
            }
        }
        eprintln!(
            "      C1 arid: **{}** bodies · **{wa:.0} km² = {:.2} % of land** · volume \
             **{va:.1} km³** · of which depth > 100 m: {} bodies, {:.1} % of area, {:.1} % of \
             volume",
            da.lakes.len(),
            100.0 * wa / land_km2,
            deep_a.0,
            100.0 * deep_a.1 / wa.max(1e-6),
            100.0 * deep_a.2 / va.max(1e-9)
        );
        let mut tva: Vec<(&String, &(usize, f32))> = ty_a.iter().collect();
        tva.sort_by(|a, b| b.1.1.total_cmp(&a.1.1));
        for (t2, (c, a)) in tva {
            eprintln!(
                "         {c:>3} × {t2:<12} {a:9.1} km² ({:5.1} %)",
                100.0 * a / wa.max(1e-6)
            );
        }
        // C2 arid, same instrument
        let mut kl = 0usize;
        let mut sc = 0usize;
        let (mut ka, mut kv, mut aa, mut av) = (0.0f32, 0.0f64, 0.0f32, 0.0f64);
        for l in da.lakes.iter().filter(|l| l.area_km2 >= 1.0) {
            let cells: Vec<usize> = (0..n).filter(|&k| da.lake_map[k] == l.base.id).collect();
            if cells.is_empty() {
                continue;
            }
            sc += 1;
            let cuts: Vec<f32> = cells.iter().map(|&k| m(&pre, k) - m(&field, k)).collect();
            let inside: HashSet<usize> = cells.iter().copied().collect();
            let mut ring: HashSet<usize> = HashSet::new();
            for &k in &cells {
                for d in 0..8 {
                    let nx = ((k % w) as i32 + D8_DX[d]).rem_euclid(w as i32) as usize;
                    let ny = ((k / w) as i32 + D8_DY[d]).rem_euclid(h as i32) as usize;
                    let nk = ny * w + nx;
                    if !inside.contains(&nk) {
                        ring.insert(nk);
                    }
                }
            }
            let rim: Vec<f32> = ring.iter().map(|&k| slope_deg(&field, k, n2m, w, h)).collect();
            let vol: f64 = cells
                .iter()
                .map(|&k| ((l.level_m - m(&field, k)).max(0.0) as f64) * cell_km2 as f64 / 1000.0)
                .sum();
            let a = cells.len() as f32 * cell_km2;
            aa += a;
            av += vol;
            if median(&cuts) > 50.0 && !rim.is_empty() && median(&rim) > 30.0 {
                kl += 1;
                ka += a;
                kv += vol;
            }
        }
        eprintln!(
            "      C2 arid: of {sc} bodies ≥ 1 km², **{kl}** in the canyon class · {ka:.0} km² = \
             {:.1} % of their area · {kv:.1} km³ = {:.1} % of their volume",
            100.0 * ka / aa.max(1e-6),
            100.0 * kv / av.max(1e-9)
        );
        // C3 arid
        let mut da_dest: HashMap<&str, (usize, f64)> = HashMap::new();
        for l in da.lakes.iter().filter(|l| l.area_km2 >= 1.0) {
            let cells: Vec<usize> = (0..n).filter(|&k| da.lake_map[k] == l.base.id).collect();
            let q = cells.iter().map(|&k| racc_a[k]).fold(0.0f32, f32::max);
            let key = match l.lake_type {
                LakeType::Exorheic => "Exorheic",
                LakeType::Endorheic => "Endorheic",
                LakeType::Unresolved => "Unresolved",
                LakeType::CraterAcidic | LakeType::CraterNeutral => "Crater",
            };
            let e = da_dest.entry(key).or_insert((0, 0.0));
            e.0 += 1;
            e.1 += runoff_km2_to_m3s(q) as f64;
        }
        let mut dv2: Vec<(&&str, &(usize, f64))> = da_dest.iter().collect();
        dv2.sort_by(|a, b| b.1.0.cmp(&a.1.0));
        for (kk, (c, q)) in dv2 {
            eprintln!("         C3 arid: {c:>3} bodies · local inflow {q:7.1} m³/s · {kk}");
        }
        eprintln!(
            "         C3 arid: end-of-chain relabels {} · `to_nothing` {} · q_lake_unrouted {:.1} \
             · passes {} (converged {})",
            rel_a.len(),
            bsa.termination.to_nothing,
            bsa.termination.q_lake_unrouted_m3s,
            bsa.termination.passes_used,
            !bsa.termination.not_converged
        );
    }

    eprintln!("\n==========  end Finding 89  ==========\n");
}
