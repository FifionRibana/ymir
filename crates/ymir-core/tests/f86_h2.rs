//! ADR 0001 Finding 86 — the H2 line (A), the non-convergence attributed (B1), the four-term
//! budget over the WHOLE water system (C1/C2), and the two unresolved basins geometrised (D).
//!
//! Every quantity carries its bed, its gate state (relevel OFF / ON) and the PASS COUNT of the
//! fixed point it was read at — Finding 85 showed the shipped loop exits on its 16-pass bound, so
//! a below-sea figure without a pass count is a figure from an unknown state.
//!
//! Run: cargo test -p ymir-core --release --test f86_h2 -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, SEA, build_field, pct, sorted};
use std::collections::{HashMap, HashSet};
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::{PrecipParams, precip_mm_per_year};
use ymir_core::export::height::metric_height_u16;
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, LakeType, MergedUnionRelevel, SegmentKind, SegmentRow,
    Spillway, apply_lake_water_balance, below_sea_basin_lakes_infil, c1_drainage_windowed,
    clip_rivers_to_lakes, exorheic_lakes_missing_outlet, potential_evaporation_mm,
    resolve_exorheic_without_outlet, runoff_accumulation, runoff_km2_to_m3s,
};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{
    D8_DX, D8_DY, DIR_NONE, FlowConfig, RiverSegment, breach_monotone, compute_flow,
};

const DOMAIN_KM: f32 = 400.0;
const RATIO: f32 = 7.5;

#[test]
#[ignore]
fn f86_h2() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 86 · A · B1 · C · D  ==========");

    // ADR Finding 86 B2 -- the relevel SHIPS, so `Default` is ON and the A/B control is `None`.
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    let mut off = on.clone();
    off.merged_union_relevel = None;
    assert!(on.merged_union_relevel.is_some(), "the relevel is not the shipped default");
    let _ = MergedUnionRelevel::default();

    let raw = build_field(Knobs::shipped());
    let (w, h) = (raw.width, raw.height);
    let n = w * h;
    let hl = metric_height_u16(&raw, &ss);
    let u16_step = (hl.max_m - hl.min_m) / 65535.0;
    let pre0 = c1_drainage_windowed(&raw, None, &off, &ss, DOMAIN_KM);
    let field = breach_monotone(&raw, &pre0.flow.filled, &pre0.lake_map, SEA, w, h);
    let flow = compute_flow(&field, &FlowConfig { sea_level: SEA, ..Default::default() });
    let wc = water_class(&field, SEA);
    let land: Vec<bool> = (0..n).map(|k| field.data[k] > SEA).collect();
    eprintln!("   field FNV unchanged by every gate here | u16 step {u16_step:.4} m");

    for (bed, lat, span) in [("humid", 45.0f32, 40.0f32), ("arid-hot", 25.0, 10.0)] {
        let climate =
            c1_climate_placed(&field, &ss, lat, span, &PrecipParams::default(), DOMAIN_KM);
        let dclim = DrainageClimate {
            precip_internal: &climate.precipitation,
            temperature: &climate.temperature,
        };
        let surplus: Vec<f32> = (0..n)
            .map(|k| {
                let p = precip_mm_per_year(climate.precipitation.data[k]);
                let pe = potential_evaporation_mm(climate.temperature.data[k]);
                (p - pe).max(0.0) * cell_km2
            })
            .collect();
        let budget = runoff_km2_to_m3s(
            (0..n).filter(|&k| land[k]).map(|k| surplus[k] as f64).sum::<f64>() as f32,
        );
        eprintln!("\n╔════════ {bed} · breached budget {budget:.1} m³/s ════════╗");

        for (gate, dcfg) in [("GATE OFF (shipped)", &off), ("**RELEVEL ON**", &on)] {
            eprintln!("\n   ── {gate} / {bed} ──");
            let pre = c1_drainage_windowed(&raw, None, dcfg, &ss, DOMAIN_KM);
            let mut dr = c1_drainage_windowed(&field, Some(&dclim), dcfg, &ss, DOMAIN_KM);
            dr.lakes = pre.lakes;
            dr.lake_map = pre.lake_map;
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
            // ── C2 · the DETECTED surface lakes' own balance, before the below-sea stage ──
            let surf: Vec<(u32, f32, f32, LakeType)> = dr
                .lakes
                .iter()
                .filter(|l| l.base.id < 1_000_001)
                .map(|l| {
                    let cells: Vec<usize> =
                        (0..n).filter(|&k| dr.lake_map[k] == l.base.id).collect();
                    let a_km2 = cells.len() as f32 * cell_km2;
                    // net evaporation over the surface, at the floor cell's climate
                    let floor = cells
                        .iter()
                        .copied()
                        .min_by(|&a, &b| field.data[a].total_cmp(&field.data[b]))
                        .unwrap_or(0);
                    let pe = potential_evaporation_mm(climate.temperature.data[floor]);
                    let pr = precip_mm_per_year(climate.precipitation.data[floor]);
                    let ev = runoff_km2_to_m3s((pe - pr).max(0.0) * a_km2);
                    (l.base.id, a_km2, ev, l.lake_type)
                })
                .collect();
            let premerge = dr.lake_map.clone();
            let bs = below_sea_basin_lakes_infil(
                &field,
                &dclim,
                dcfg,
                &ss,
                DOMAIN_KM,
                Some(&premerge),
                None,
            );
            let t = &bs.termination;

            // ── B1 · what the fixed point does while it spins ────────────────
            eprintln!(
                "      B1 passes {} (converged {}) | trace (pass · classes changed · worst m³/s · \
                 on class):",
                t.passes_used, !t.not_converged
            );
            for (p, c, q, r) in &t.pass_trace {
                eprintln!("         pass {p:>2} · changed {c:>3} · worst {q:>10.3} · class {r}");
            }
            let tail: Vec<&(u32, usize, f32, u32)> =
                t.pass_trace.iter().rev().take(6).collect::<Vec<_>>();
            let classes: HashSet<u32> = tail.iter().map(|x| x.3).collect();
            eprintln!(
                "         last 6 passes: changed {:?} · worst {:?} · classes {:?}",
                tail.iter().map(|x| x.1).collect::<Vec<_>>(),
                tail.iter().map(|x| format!("{:.3}", x.2)).collect::<Vec<_>>(),
                classes
            );

            // assemble the tail exactly as production does
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
                    segment: RiverSegment {
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

            // ── A · the H2 population, named, BEFORE the line runs ────────────
            let flagged = exorheic_lakes_missing_outlet(&dr);
            eprintln!(
                "      A: exorheic WITHOUT an outlet reach: **{}** {flagged:?}",
                flagged.len()
            );
            for &id in &flagged {
                let lk = dr.lakes.iter().find(|l| l.base.id == id).unwrap();
                let cells = (0..n).filter(|&k| dr.lake_map[k] == id).count();
                eprintln!(
                    "         id {id} — {} | level {:.2} m | area {:.3} km² | {cells} cells | \
                     type {:?}",
                    if id >= 1_000_001 { "BELOW-SEA basin" } else { "DETECTED surface lake" },
                    lk.level_m,
                    lk.area_km2,
                    lk.lake_type
                );
            }
            let before: Vec<LakeType> = dr.lakes.iter().map(|l| l.lake_type).collect();
            let resolved = resolve_exorheic_without_outlet(&mut dr);
            let unresolved =
                dr.lakes.iter().filter(|l| l.lake_type == LakeType::Unresolved).count();
            let moved = dr.lakes.iter().zip(&before).filter(|(a, b)| a.lake_type != **b).count();
            eprintln!(
                "      A1: `resolve_exorheic_without_outlet` relabelled {} | Unresolved now \
                 **{unresolved}** | lakes whose type moved {moved} | exorheic-without-outlet AFTER \
                 **{}**",
                resolved.len(),
                exorheic_lakes_missing_outlet(&dr).len()
            );

            // ── C1 · the four terms, over the WHOLE system ───────────────────
            let (mut wc_n, mut wc_q) = (0usize, 0.0f64);
            for (i, s) in dr.rivers.segments.iter().enumerate() {
                if dr.segment_kind[i] != SegmentKind::Watercourse || s.downstream.is_some() {
                    continue;
                }
                let &(ex, ey) = s.points.last().unwrap();
                let k = ey as usize * w + ex as usize;
                let mut sea = wc[k] == 1;
                let d = flow.direction[k];
                if !sea && d != DIR_NONE {
                    let nx = ((ex as i32 + D8_DX[d as usize]).rem_euclid(w as i32)) as usize;
                    let ny = ((ey as i32 + D8_DY[d as usize]).rem_euclid(h as i32)) as usize;
                    sea = wc[ny * w + nx] == 1;
                }
                if sea {
                    wc_n += 1;
                    wc_q += dr.segment_discharge_m3s[i] as f64;
                }
            }
            let end_wc = |sw: &Spillway| -> u8 {
                let &(lx, ly) = sw.points.last().unwrap();
                wc[ly as usize * w + lx as usize]
            };
            let sp_ocean: f64 = bs
                .spillways
                .iter()
                .filter(|s| end_wc(s) == 1)
                .map(|s| s.discharge_m3s as f64)
                .sum();
            let to_ocean = wc_q + sp_ocean;
            let evap_below: f64 = bs.basins.iter().map(|b| b.evaporation_m3s as f64).sum();
            let evap_surface: f64 = surf.iter().map(|s| s.2 as f64).sum();
            let mut out_of: HashMap<u32, f32> = HashMap::new();
            for sw in &bs.spillways {
                *out_of.entry(sw.lake_id).or_default() += sw.discharge_m3s;
            }
            // WARNING `inflow_m3s` is Finding 40's CHAINED inflow: summing it counts the same
            // water once per link. `local_inflow_m3s` (Finding 86) is the basin's own catchment
            // and is the only one that sums. Both are printed so the difference is visible.
            let local_in: f64 = bs.basins.iter().map(|b| b.local_inflow_m3s as f64).sum();
            let chained_in: f64 = bs.basins.iter().map(|b| b.inflow_m3s as f64).sum();
            let retained: f64 = bs
                .basins
                .iter()
                .filter(|b| out_of.get(&b.id).copied().unwrap_or(0.0) == 0.0)
                .map(|b| (b.local_inflow_m3s - b.evaporation_m3s).max(0.0) as f64)
                .sum();
            let spill_out: f64 = bs.spillways.iter().map(|s| s.discharge_m3s as f64).sum();
            eprintln!(
                "      C1a the BELOW-SEA system on LOCAL inflow: local {local_in:.1} (chained {chained_in:.1}, x{:.2}) = spillways out {spill_out:.1} + evap {evap_below:.2} + retained {retained:.1} => {:.1} % of local",
                chained_in / local_in.max(1e-9),
                100.0 * (spill_out + evap_below + retained) / local_in.max(1e-9)
            );
            // ── the FIFTH term, measured instead of named: how much runoff physically
            // ARRIVES at the coast, against how much the segment network credits to the sea.
            // The budget is a per-CELL sum; `Watercourse -> sea` is a per-SEGMENT sum over
            // terminal segments only. They are not the same population and this is the gap.
            let racc = runoff_accumulation(&field, &flow, &dclim, cell_km2, None, None, w, h);
            let mut arrives = 0.0f64;
            for k in 0..n {
                if !land[k] || flow.direction[k] == DIR_NONE {
                    continue;
                }
                let d = flow.direction[k] as usize;
                let nx = ((k % w) as i32 + D8_DX[d]).rem_euclid(w as i32) as usize;
                let ny = ((k / w) as i32 + D8_DY[d]).rem_euclid(h as i32) as usize;
                if wc[ny * w + nx] == 1 {
                    arrives += runoff_km2_to_m3s(racc[k]) as f64;
                }
            }
            eprintln!(
                "      C1c the FIFTH term, MEASURED: runoff ARRIVING at the coast (last land cell whose D8 receiver is ocean) = **{arrives:.1} m3/s** = **{:.1} %** of the budget, against `Watercourse -> sea` = {wc_q:.1} ({:.1} %) over {wc_n} terminal segments. The gap **{:.1} m3/s** is runoff that reaches the sea without being credited to a terminal Watercourse segment.",
                100.0 * arrives / budget as f64,
                100.0 * wc_q / budget as f64,
                arrives - wc_q
            );
            let sum5 = arrives + sp_ocean + evap_below + evap_surface + retained;
            eprintln!(
                "      C1d FIVE TERMS: coastal arrival {arrives:.1} + Spillway->ocean {sp_ocean:.1} + evap below {evap_below:.2} + evap surface {evap_surface:.2} + retained {retained:.1} = **{sum5:.1}** = **{:.1} %** of the budget {budget:.1}",
                100.0 * sum5 / budget as f64
            );
            let sum4 = to_ocean + evap_below + evap_surface + retained;
            eprintln!(
                "      C1b FOUR TERMS against the breached budget {budget:.1}: ocean **{to_ocean:.1}** ({:.1} %) | evap BELOW-SEA {evap_below:.2} ({:.1} %) | evap SURFACE LAKES **{evap_surface:.2}** ({:.1} %) | RETAINED in outflowless basins, LOCAL inflow minus evap **{retained:.1}** ({:.1} %) | SUM {sum4:.1} = **{:.1} %** of the budget, missing **{:.1} m3/s**",
                100.0 * to_ocean / budget as f64,
                100.0 * evap_below / budget as f64,
                100.0 * evap_surface / budget as f64,
                100.0 * retained / budget as f64,
                100.0 * sum4 / budget as f64,
                budget as f64 - sum4
            );
            eprintln!(
                "      C2 surface lakes: {} | Σ area {:.1} km² | Σ net evaporation {:.2} m³/s | \
                 types {:?}",
                surf.len(),
                surf.iter().map(|s| s.1 as f64).sum::<f64>(),
                evap_surface,
                {
                    let mut ts: Vec<String> = surf.iter().map(|s| format!("{:?}", s.3)).collect();
                    ts.sort_unstable();
                    ts.dedup();
                    ts
                }
            );

            // ── D · the basins with inflow and no outflow, geometrised ───────
            let mut resid: Vec<&ymir_core::tectonics_c1::drainage::BasinSummary> = bs
                .basins
                .iter()
                .filter(|b| {
                    out_of.get(&b.id).copied().unwrap_or(0.0) == 0.0
                        && b.local_inflow_m3s - b.evaporation_m3s > 0.5
                })
                .collect();
            resid.sort_by(|a, b| {
                (b.local_inflow_m3s - b.evaporation_m3s)
                    .partial_cmp(&(a.local_inflow_m3s - a.evaporation_m3s))
                    .unwrap()
            });
            eprintln!("      D: {} basin(s) with inflow and no outflow", resid.len());
            for b in resid.iter().take(3) {
                let cells: Vec<usize> = (0..n).filter(|&k| bs.lake_map[k] == b.id).collect();
                // the pour point: the lowest cell adjacent to the footprint and not in it
                let inside: HashSet<usize> = cells.iter().copied().collect();
                let mut col: Option<usize> = None;
                for &k in &cells {
                    let (x, y) = ((k % w) as i32, (k / w) as i32);
                    for dy in -1i32..=1 {
                        for dx in -1i32..=1 {
                            let (nx, ny) = (x + dx, y + dy);
                            if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                                continue;
                            }
                            let nk = ny as usize * w + nx as usize;
                            if !inside.contains(&nk)
                                && col.is_none_or(|c| field.data[nk] < field.data[c])
                            {
                                col = Some(nk);
                            }
                        }
                    }
                }
                let Some(c) = col else { continue };
                // the lowest neighbour of the col that is outside the footprint
                let (cx, cy) = ((c % w) as i32, (c / w) as i32);
                let mut esc: Option<usize> = None;
                for dy in -1i32..=1 {
                    for dx in -1i32..=1 {
                        let (nx, ny) = (cx + dx, cy + dy);
                        if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                            continue;
                        }
                        let nk = ny as usize * w + nx as usize;
                        if !inside.contains(&nk)
                            && nk != c
                            && esc.is_none_or(|e| field.data[nk] < field.data[e])
                        {
                            esc = Some(nk);
                        }
                    }
                }
                let drop_steps = esc
                    .map(|e| (field.data[c] - field.data[e]) * n2m / u16_step)
                    .unwrap_or(f32::NAN);
                // where the descent from the escape leads
                let verdict = match esc {
                    None => "no neighbour outside the footprint at all".to_string(),
                    Some(e0) => {
                        let (mut cur, mut steps) = (e0, 0usize);
                        loop {
                            if wc[cur] == 1 {
                                break "the OCEAN".to_string();
                            }
                            if inside.contains(&cur) {
                                break "BACK INTO ITS OWN footprint (the trace is a loop)"
                                    .to_string();
                            }
                            if bs.lake_map[cur] != 0 {
                                break format!("another below-sea body, id {}", bs.lake_map[cur]);
                            }
                            if dr.lake_map[cur] != 0 && dr.lake_map[cur] < 1_000_001 {
                                break format!("a DETECTED surface lake, id {}", dr.lake_map[cur]);
                            }
                            let d = flow.direction[cur];
                            if d == DIR_NONE {
                                break "a FLAT — no D8 direction at all".to_string();
                            }
                            let nx = ((cur % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32);
                            let ny = ((cur / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32);
                            cur = ny as usize * w + nx as usize;
                            steps += 1;
                            if steps > 4 * w {
                                break "nowhere in 4w steps".to_string();
                            }
                        }
                    }
                };
                eprintln!(
                    "         id {} — inflow {:.2} evap {:.2} · a_eq {:.3} a_spill {:.3} · level \
                     {:.2} m sill {:.2} m floor {:.2} m · has_sill {} exorheic {}\n            col \
                     at ({}, {}) h {:.2} m · drop to its escape **{:.2} u16 steps** · the descent \
                     leads to **{verdict}**",
                    b.id,
                    b.local_inflow_m3s,
                    b.evaporation_m3s,
                    b.a_eq_km2,
                    b.a_spill_km2,
                    b.level_m,
                    b.spill_level_m,
                    b.floor_m,
                    b.has_sill,
                    b.exorheic,
                    c % w,
                    c / w,
                    (field.data[c] - SEA) * n2m,
                    drop_steps
                );
                // a_eq == a_spill means a_eq is INFINITE, reported as the sill area
                // (drainage.rs:2170). Inverting it as a rate is meaningless -- Finding 85 read
                // the equality as a knife-edge and it is a reporting convention.
                eprintln!(
                    "            a_eq == a_spill here means a_eq is INFINITE and reported as the sill area (drainage.rs:2170): net_evap = 0, evaporation {:.2} m3/s over {:.3} km2. NOT a knife-edge equality.",
                    b.evaporation_m3s, b.area_km2
                );
            }

            // widths, for the visual block
            let mut max_sp = (0.0f32, 0.0f32);
            for sw in &bs.spillways {
                if sw.discharge_m3s > max_sp.0 {
                    max_sp = (sw.discharge_m3s, sw.width_m * RATIO / 48.83);
                }
            }
            let idx: Vec<usize> = (0..dr.rivers.segments.len())
                .filter(|&i| dr.segment_kind[i] == SegmentKind::Watercourse)
                .collect();
            let mut wcell: Vec<f32> =
                idx.iter().map(|&i| dr.segment_width_m[i] * RATIO / 48.83).collect();
            wcell.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            eprintln!(
                "      widths: max Spillway {:.2} m³/s = {:.2} cells | max Watercourse cells \
                 {:.2} | p50 {:.3}",
                max_sp.0,
                max_sp.1,
                wcell.last().copied().unwrap_or(0.0),
                pct(&wcell, 0.50)
            );
            let _ = sorted(vec![0.0f32]);
        }
    }
}
