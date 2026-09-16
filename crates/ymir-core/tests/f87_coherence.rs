//! ADR 0001 Finding 87 — physical coherence: the spillway drops (A), the 623 m lake's reason (B),
//! the 47 km² marsh (C), lake 55 traced (D), the col contradiction (E), the local outflow and the
//! six-term closure (F), and whether the diffuse coastal arrival should be channelised (G).
//!
//! Every quantity carries its bed, its gate state and the pass count of the fixed point.
//!
//! Run: cargo test -p ymir-core --release --test f87_coherence -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, SEA, build_field, pct, sorted};
use std::collections::{HashMap, HashSet};
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::{PrecipParams, precip_mm_per_year};
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, LakeType, SegmentKind, SegmentRow, Spillway,
    apply_lake_water_balance, below_sea_basin_lakes_infil, c1_drainage_windowed,
    clip_rivers_to_lakes, potential_evaporation_mm, resolve_exorheic_without_outlet,
    runoff_accumulation, runoff_km2_to_m3s, surface_lake_outlet_trace,
};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{
    D8_DX, D8_DY, DIR_NONE, FlowConfig, RiverSegment, breach_monotone, compute_flow,
};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const RATIO: f32 = 7.5;

#[test]
#[ignore]
fn f87_coherence() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 87 · A · B · C · D · E · F · G  ==========");

    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    let mut off = on.clone();
    off.merged_union_relevel = None;
    let head_km2 = on.thresholds.head_km2;

    let raw = build_field(Knobs::shipped());
    let (w, h) = (raw.width, raw.height);
    let n = w * h;
    let pre0 = c1_drainage_windowed(&raw, None, &on, &ss, DOMAIN_KM);
    let field = breach_monotone(&raw, &pre0.flow.filled, &pre0.lake_map, SEA, w, h);
    let flow = compute_flow(&field, &FlowConfig { sea_level: SEA, ..Default::default() });
    let wc = water_class(&field, SEA);
    let land: Vec<bool> = (0..n).map(|k| field.data[k] > SEA).collect();
    // block B needs the PRE-INCISION field at the same cells
    let pre = build_field(Knobs::no_incision());
    eprintln!("   head_threshold (bench) = {head_km2} km² | cell {CELL_M:.2} m");

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

        for (gate, dcfg) in [("**SHIPPED (relevel ON)**", &on), ("A/B control (OFF)", &off)] {
            eprintln!("\n   ── {gate} / {bed} ──");
            let pre_d = c1_drainage_windowed(&raw, None, dcfg, &ss, DOMAIN_KM);
            let mut dr = c1_drainage_windowed(&field, Some(&dclim), dcfg, &ss, DOMAIN_KM);
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
            // ── D · lake 55, and every detected lake the balance calls exorheic ──
            let racc = runoff_accumulation(&field, &flow, &dclim, cell_km2, None, None, w, h);
            let surf_map = dr.lake_map.clone();
            let mut unres_surface: Vec<(u32, f32, f32, &'static str)> = Vec::new();
            for l in dr.lakes.iter().filter(|l| l.base.id < 1_000_001) {
                if l.lake_type != LakeType::Unresolved {
                    continue;
                }
                let cells: Vec<usize> = (0..n).filter(|&k| surf_map[k] == l.base.id).collect();
                let inflow = cells.iter().map(|&k| racc[k]).fold(0.0f32, f32::max);
                // ⚠️ `&dr.flow` and NOT the bench's own `flow`: production routes on the
                // drainage's flow field, which carries `FlatPerturbation`. Tracing on a plain
                // `compute_flow` disagreed with production on 15 of these 19 lakes — the
                // disagreement WAS the perturbation, and reading the reason off the wrong field
                // is reading nothing.
                let reason = surface_lake_outlet_trace(
                    l.base.id,
                    l.base.outlet,
                    &dr.flow,
                    &field,
                    &surf_map,
                    w,
                    h,
                )
                .err()
                // ADR Finding 88-C — the reason is a typed `UnresolvedReason` now.
                .map_or("(the trace SUCCEEDS on the final map — see the note)", |r| r.as_str());
                unres_surface.push((l.base.id, l.area_km2, runoff_km2_to_m3s(inflow), reason));
            }
            eprintln!(
                "      D: detected surface lakes the BALANCE now leaves Unresolved: **{}**",
                unres_surface.len()
            );
            for (id, a, q, why) in &unres_surface {
                eprintln!(
                    "         lake {id} — {a:.3} km² · local inflow **{q:.2} m³/s** · the trace \
                     fails because: **{why}**"
                );
            }

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
            let relabelled = resolve_exorheic_without_outlet(&mut dr);
            eprintln!(
                "      D guard: end-of-chain relabels {} | Unresolved total {} | passes {} \
                 (converged {})",
                relabelled.len(),
                dr.lakes.iter().filter(|l| l.lake_type == LakeType::Unresolved).count(),
                t.passes_used,
                !t.not_converged
            );

            // ── A · every spillway's drop and slope ──────────────────────────
            let level_of: HashMap<u32, f32> =
                dr.lakes.iter().map(|l| (l.base.id, l.level_m)).collect();
            let mut rows: Vec<(f32, u32, u32, f32, f32, f32, usize, f32, f32)> = Vec::new();
            for sw in &bs.spillways {
                let up = level_of.get(&sw.lake_id).copied().unwrap_or(f32::NAN);
                let &(lx, ly) = sw.points.last().unwrap();
                let end = ly as usize * w + lx as usize;
                let (down, dst) = if wc[end] == 1 {
                    (0.0f32, 0u32) // the sea
                } else if dr.lake_map[end] != 0 {
                    (level_of.get(&dr.lake_map[end]).copied().unwrap_or(f32::NAN), dr.lake_map[end])
                } else {
                    (c1_altitude_norm_to_metres(field.data[end], &ss), u32::MAX)
                };
                let mut len_m = 0.0f32;
                for p in sw.points.windows(2) {
                    let (dx, dy) = (
                        (p[1].0 as f32 - p[0].0 as f32).abs(),
                        (p[1].1 as f32 - p[0].1 as f32).abs(),
                    );
                    len_m += (dx * dx + dy * dy).sqrt() * CELL_M;
                }
                let drop = up - down;
                let slope = if len_m > 0.0 { drop / len_m } else { f32::NAN };
                rows.push((
                    sw.discharge_m3s,
                    sw.lake_id,
                    dst,
                    up,
                    down,
                    drop,
                    sw.points.len(),
                    slope,
                    sw.width_m * RATIO / 48.83,
                ));
            }
            rows.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
            let flat = rows.iter().filter(|r| r.0 > 0.0 && r.5 <= 0.0).count();
            eprintln!(
                "      A: {} spillways | with Q > 0 and **drop <= 0: {flat}** (a strait, not a \
                 river)",
                rows.len()
            );
            eprintln!(
                "         {:>9} {:>9} {:>11} {:>9} {:>9} {:>9} {:>6} {:>11} {:>7}",
                "Q m³/s", "from", "-> to", "up m", "down m", "drop m", "cells", "slope", "cells w"
            );
            for r in rows.iter().take(8) {
                eprintln!(
                    "         {:>9.2} {:>9} {:>11} {:>9.2} {:>9.2} {:>9.3} {:>6} {:>11.2e} {:>7.2}",
                    r.0,
                    r.1,
                    if r.2 == 0 {
                        "SEA".to_string()
                    } else if r.2 == u32::MAX {
                        "land".to_string()
                    } else {
                        r.2.to_string()
                    },
                    r.3,
                    r.4,
                    r.5,
                    r.6,
                    r.7,
                    r.8
                );
            }
            let slopes =
                sorted(rows.iter().filter(|r| r.7.is_finite() && r.5 > 0.0).map(|r| r.7).collect());
            if !slopes.is_empty() {
                eprintln!(
                    "         slopes of the {} spillways with a positive drop: p10 {:.2e} p50 \
                     {:.2e} p90 {:.2e} max {:.2e} | anchors: Detroit River 2e-5, rock sill 1e-2, \
                     waterfall > 1e-1",
                    slopes.len(),
                    pct(&slopes, 0.10),
                    pct(&slopes, 0.50),
                    pct(&slopes, 0.90),
                    slopes.last().copied().unwrap_or(0.0)
                );
            }

            // ── F · the closure, on LOCAL inflow AND local outflow ───────────
            let local_in: f64 = bs.basins.iter().map(|b| b.local_inflow_m3s as f64).sum();
            let local_out: f64 = bs.basins.iter().map(|b| b.local_outflow_m3s as f64).sum();
            let evap_below: f64 = bs.basins.iter().map(|b| b.evaporation_m3s as f64).sum();
            let mut out_of: HashMap<u32, f32> = HashMap::new();
            for sw in &bs.spillways {
                *out_of.entry(sw.lake_id).or_default() += sw.discharge_m3s;
            }
            let retained: f64 = bs
                .basins
                .iter()
                .filter(|b| out_of.get(&b.id).copied().unwrap_or(0.0) == 0.0)
                .map(|b| b.local_outflow_m3s as f64)
                .sum();
            eprintln!(
                "      F: below-sea on LOCAL BOTH SIDES: in {local_in:.1} = out {local_out:.1} + \
                 evap {evap_below:.2} (retained, no spillway: {retained:.1}) ⇒ **{:.1} %**",
                100.0 * (local_out + evap_below) / local_in.max(1e-9)
            );
            // the coastal arrival (Finding 86's fifth term) and G in one pass
            let mut arrives = 0.0f64;
            let on_segment: HashSet<usize> = dr
                .rivers
                .segments
                .iter()
                .flat_map(|s| s.points.iter())
                .map(|&(x, y)| y as usize * w + x as usize)
                .collect();
            let (mut ungauged, mut ung_q) = (Vec::<f32>::new(), 0.0f64);
            let mut gauged_q = 0.0f64;
            for k in 0..n {
                if !land[k] || flow.direction[k] == DIR_NONE {
                    continue;
                }
                let d = flow.direction[k] as usize;
                let nx = ((k % w) as i32 + D8_DX[d]).rem_euclid(w as i32) as usize;
                let ny = ((k / w) as i32 + D8_DY[d]).rem_euclid(h as i32) as usize;
                if wc[ny * w + nx] != 1 {
                    continue;
                }
                let q = runoff_km2_to_m3s(racc[k]) as f64;
                arrives += q;
                if on_segment.contains(&k) {
                    gauged_q += q;
                } else {
                    ung_q += q;
                    ungauged.push(flow.accumulation.data[k] * cell_km2);
                }
            }
            let evap_surface: f64 = dr
                .lakes
                .iter()
                .filter(|l| l.base.id < 1_000_001)
                .map(|l| {
                    let cells: Vec<usize> =
                        (0..n).filter(|&k| dr.lake_map[k] == l.base.id).collect();
                    let a = cells.len() as f32 * cell_km2;
                    let fl = cells
                        .iter()
                        .copied()
                        .min_by(|&x, &y| field.data[x].total_cmp(&field.data[y]))
                        .unwrap_or(0);
                    let pe = potential_evaporation_mm(climate.temperature.data[fl]);
                    let pr = precip_mm_per_year(climate.precipitation.data[fl]);
                    runoff_km2_to_m3s((pe - pr).max(0.0) * a) as f64
                })
                .sum();
            let sp_ocean: f64 = bs
                .spillways
                .iter()
                .filter(|s| {
                    let &(lx, ly) = s.points.last().unwrap();
                    wc[ly as usize * w + lx as usize] == 1
                })
                .map(|s| s.discharge_m3s as f64)
                .sum();
            // ADR Finding 87-F — the ocean term in LOCAL units too: the emitted spillway
            // discharge is chained, so summing it re-imports the duplication.
            let ocean_ids: HashSet<u32> = bs
                .spillways
                .iter()
                .filter(|s| {
                    let &(lx, ly) = s.points.last().unwrap();
                    wc[ly as usize * w + lx as usize] == 1
                })
                .map(|s| s.lake_id)
                .collect();
            let sp_ocean_local: f64 = bs
                .basins
                .iter()
                .filter(|b| ocean_ids.contains(&b.id))
                .map(|b| b.local_outflow_m3s as f64)
                .sum();
            let unres_q: f64 = unres_surface.iter().map(|u| u.2 as f64).sum();
            let six = arrives + sp_ocean + evap_below + evap_surface + retained + unres_q;
            let six_local =
                arrives + sp_ocean_local + evap_below + evap_surface + retained + unres_q;
            eprintln!(
                "      F SIX TERMS (chained ocean): coastal {arrives:.1} + Spillway-ocean {sp_ocean:.1} + evap below {evap_below:.2} + evap surface {evap_surface:.2} + retained {retained:.1} + Unresolved lakes {unres_q:.1} = {six:.1} = **{:.1} %**",
                100.0 * six / budget as f64
            );
            eprintln!(
                "      F SIX TERMS (**LOCAL ocean**): {arrives:.1} + **{sp_ocean_local:.1}** + {evap_below:.2} + {evap_surface:.2} + {retained:.1} + {unres_q:.1} = **{six_local:.1}** = **{:.1} %** of {budget:.1}",
                100.0 * six_local / budget as f64
            );

            // ── G · should the ungauged arrival be channelised? ──────────────
            let ug = sorted(ungauged.clone());
            let above = ug.iter().filter(|&&a| a >= head_km2).count();
            eprintln!(
                "      G: coastal arrival {arrives:.1} m³/s, of which **{ung_q:.1} ({:.1} %) on \
                 cells carrying NO segment** ({} cells) · their D8 drainage area km²: p50 {:.4} \
                 p90 {:.3} p99 {:.3} max {:.1} · **{:.1} % are at or above head_threshold \
                 {head_km2}**",
                100.0 * ung_q / arrives.max(1e-9),
                ug.len(),
                pct(&ug, 0.50),
                pct(&ug, 0.90),
                pct(&ug, 0.99),
                ug.last().copied().unwrap_or(0.0),
                100.0 * above as f64 / ug.len().max(1) as f64
            );
            let _ = gauged_q;

            if gate.contains("control") {
                continue; // B, C and E are read on the shipped state only
            }

            // ── B · the deepest below-sea body, and WHY it is deep ───────────
            let deepest = dr
                .lakes
                .iter()
                .filter(|l| l.base.id >= 1_000_001)
                .max_by(|a, b| a.level_m.total_cmp(&b.level_m));
            if let Some(lk) = deepest {
                let cells: Vec<usize> = (0..n).filter(|&k| dr.lake_map[k] == lk.base.id).collect();
                let floor_k = cells
                    .iter()
                    .copied()
                    .min_by(|&a, &b| field.data[a].total_cmp(&field.data[b]))
                    .unwrap();
                let del_floor = (field.data[floor_k] - SEA) * n2m;
                let pre_floor = (pre.data[floor_k] - SEA) * n2m;
                // the same statistic over the whole footprint, not just one cell
                let del = sorted(cells.iter().map(|&k| (field.data[k] - SEA) * n2m).collect());
                let prf = sorted(cells.iter().map(|&k| (pre.data[k] - SEA) * n2m).collect());
                // rim slope: max 8-neighbour gradient at cells just outside the footprint
                let inside: HashSet<usize> = cells.iter().copied().collect();
                let mut rim: Vec<f32> = Vec::new();
                for &k in &cells {
                    let (x, y) = ((k % w) as i32, (k / w) as i32);
                    for dy in -1i32..=1 {
                        for dx in -1i32..=1 {
                            let (nx, ny) = (x + dx, y + dy);
                            if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                                continue;
                            }
                            let nk = ny as usize * w + nx as usize;
                            if inside.contains(&nk) {
                                continue;
                            }
                            let mut g = 0.0f32;
                            for ddy in -1i32..=1 {
                                for ddx in -1i32..=1 {
                                    let (mx, my) = (nx + ddx, ny + ddy);
                                    if mx < 0 || my < 0 || mx as usize >= w || my as usize >= h {
                                        continue;
                                    }
                                    let mk = my as usize * w + mx as usize;
                                    let dist =
                                        ((ddx * ddx + ddy * ddy) as f32).sqrt().max(1.0) * CELL_M;
                                    g = g.max(
                                        ((field.data[mk] - field.data[nk]).abs() * n2m) / dist,
                                    );
                                }
                            }
                            rim.push(g.atan().to_degrees());
                        }
                    }
                }
                let rim = sorted(rim);
                eprintln!(
                    "      B: the deepest below-sea body is id {} — level **{:.1} m**, depth \
                     **{:.1} m**, {} cells ({:.1} km²)",
                    lk.base.id,
                    lk.level_m,
                    lk.depth_m,
                    cells.len(),
                    cells.len() as f32 * cell_km2
                );
                eprintln!(
                    "         the FLOOR cell ({}, {}): DELIVERED {:.2} m · PRE-INCISION **{:.2} \
                     m** ⇒ the incision cut **{:.2} m** there",
                    floor_k % w,
                    floor_k / w,
                    del_floor,
                    pre_floor,
                    pre_floor - del_floor
                );
                eprintln!(
                    "         over the whole footprint: delivered p10 {:.1} p50 {:.1} p90 {:.1} · \
                     pre-incision p10 {:.1} p50 {:.1} p90 {:.1} · median cut **{:.2} m**",
                    pct(&del, 0.10),
                    pct(&del, 0.50),
                    pct(&del, 0.90),
                    pct(&prf, 0.10),
                    pct(&prf, 0.50),
                    pct(&prf, 0.90),
                    pct(&prf, 0.50) - pct(&del, 0.50)
                );
                eprintln!(
                    "         rim slope (degrees, max 8-neighbour gradient just outside): p10 \
                     {:.1} p50 **{:.1}** p90 {:.1} · anchors: crater wall 20–40°, rift flank \
                     5–15°",
                    pct(&rim, 0.10),
                    pct(&rim, 0.50),
                    pct(&rim, 0.90)
                );
            }

            // ── C · the union that is a marsh ────────────────────────────────
            let marsh = dr
                .lakes
                .iter()
                .filter(|l| l.base.id >= 1_000_001 && l.depth_m < 3.0)
                .max_by(|a, b| a.area_km2.total_cmp(&b.area_km2));
            if let Some(lk) = marsh {
                let cells: Vec<usize> = (0..n).filter(|&k| dr.lake_map[k] == lk.base.id).collect();
                let dep = sorted(
                    cells
                        .iter()
                        .map(|&k| lk.level_m - c1_altitude_norm_to_metres(field.data[k], &ss))
                        .collect(),
                );
                let wet = cells.iter().filter(|&&k| bs.wetland[k] == 1).count();
                eprintln!(
                    "      C: the largest shallow union is id {} — {:.3} km², level {:.2} m, \
                     `lake_type` **{:?}** · depth p50 **{:.2} m** p90 {:.2} max {:.2} · **{:.1} % \
                     of its footprint is `wetland`** ({wet} of {} cells)",
                    lk.base.id,
                    lk.area_km2,
                    lk.level_m,
                    lk.lake_type,
                    pct(&dep, 0.50),
                    pct(&dep, 0.90),
                    dep.last().copied().unwrap_or(0.0),
                    100.0 * wet as f64 / cells.len().max(1) as f64,
                    cells.len()
                );
            }

            // ── E · the col, searched the way a flood escapes ────────────────
            let mut resid: Vec<&ymir_core::tectonics_c1::drainage::BasinSummary> = bs
                .basins
                .iter()
                .filter(|b| {
                    out_of.get(&b.id).copied().unwrap_or(0.0) == 0.0 && b.local_outflow_m3s > 0.5
                })
                .collect();
            resid.sort_by(|a, b| b.local_outflow_m3s.total_cmp(&a.local_outflow_m3s));
            for b in resid.iter().take(2) {
                let cells: Vec<usize> = (0..n).filter(|&k| bs.lake_map[k] == b.id).collect();
                let inside: HashSet<usize> = cells.iter().copied().collect();
                // the NAIVE col of Finding 86: the lowest cell adjacent to the footprint
                let mut naive: Option<usize> = None;
                // the REAL saddle: the lowest such cell that HAS a strictly lower neighbour
                // outside the footprint — which is what "the flood escapes here" means
                let mut saddle: Option<(usize, usize)> = None;
                for &k in &cells {
                    let (x, y) = ((k % w) as i32, (k / w) as i32);
                    for dy in -1i32..=1 {
                        for dx in -1i32..=1 {
                            let (nx, ny) = (x + dx, y + dy);
                            if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                                continue;
                            }
                            let c = ny as usize * w + nx as usize;
                            if inside.contains(&c) {
                                continue;
                            }
                            if naive.is_none_or(|b0| field.data[c] < field.data[b0]) {
                                naive = Some(c);
                            }
                            let (cx, cy) = ((c % w) as i32, (c / w) as i32);
                            let mut best: Option<usize> = None;
                            for ey in -1i32..=1 {
                                for ex in -1i32..=1 {
                                    let (ax, ay) = (cx + ex, cy + ey);
                                    if ax < 0 || ay < 0 || ax as usize >= w || ay as usize >= h {
                                        continue;
                                    }
                                    let e = ay as usize * w + ax as usize;
                                    if inside.contains(&e) || e == c {
                                        continue;
                                    }
                                    if field.data[e] < field.data[c]
                                        && best.is_none_or(|b1| field.data[e] < field.data[b1])
                                    {
                                        best = Some(e);
                                    }
                                }
                            }
                            if let Some(e) = best {
                                if saddle.is_none_or(|(s, _)| field.data[c] < field.data[s]) {
                                    saddle = Some((c, e));
                                }
                            }
                        }
                    }
                }
                eprintln!(
                    "      E: basin {} — level {:.2} m · `spill_level_m` **{:.2} m** · a_eq {} \
                     a_spill {:.3} · local in {:.2} out {:.2}",
                    b.id,
                    b.level_m,
                    b.spill_level_m,
                    if b.a_eq_is_infinite {
                        "INFINITE".to_string()
                    } else {
                        format!("{:.3}", b.a_eq_km2)
                    },
                    b.a_spill_km2,
                    b.local_inflow_m3s,
                    b.local_outflow_m3s
                );
                if let Some(c) = naive {
                    eprintln!(
                        "         the NAIVE col (Finding 86's instrument: lowest cell adjacent to \
                         the footprint) is at {:.2} m — **{:.2} m ABOVE the level**, which is what \
                         it must be by construction: the footprint IS every cell at or below the \
                         level",
                        (field.data[c] - SEA) * n2m,
                        (field.data[c] - SEA) * n2m - b.level_m
                    );
                }
                match saddle {
                    None => eprintln!(
                        "         the REAL saddle: **none exists** — no cell adjacent to the \
                         footprint has a strictly lower neighbour outside it. The basin is closed."
                    ),
                    Some((c, e)) => {
                        let mut cur = e;
                        let mut steps = 0usize;
                        let verdict = loop {
                            if wc[cur] == 1 {
                                break "the OCEAN".to_string();
                            }
                            if inside.contains(&cur) {
                                break "BACK into its own footprint".to_string();
                            }
                            if bs.lake_map[cur] != 0 {
                                break format!("below-sea body {}", bs.lake_map[cur]);
                            }
                            let d = flow.direction[cur];
                            if d == DIR_NONE {
                                break "a flat (no D8 direction)".to_string();
                            }
                            let nx = ((cur % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32);
                            let ny = ((cur / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32);
                            cur = ny as usize * w + nx as usize;
                            steps += 1;
                            if steps > 4 * w {
                                break "nowhere in 4w steps".to_string();
                            }
                        };
                        eprintln!(
                            "         the REAL saddle is ({}, {}) at {:.2} m, escaping to ({}, \
                             {}) at {:.2} m — a drop of {:.3} m — and the descent leads to \
                             **{verdict}**",
                            c % w,
                            c / w,
                            (field.data[c] - SEA) * n2m,
                            e % w,
                            e / w,
                            (field.data[e] - SEA) * n2m,
                            (field.data[c] - field.data[e]) * n2m
                        );
                    }
                }
            }
        }
    }
}
