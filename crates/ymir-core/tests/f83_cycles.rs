//! ADR 0001 Finding 83 — block A's hydrological guards **on the production path**, and block C2,
//! the cycle-breaker re-read on the world the bound now makes.
//!
//! Since Finding 83 `Knobs::shipped()` IS the bounded field, so this bench runs the production
//! config and nothing else. The Finding 74 terminal table, the lake invariants and the
//! `BasinSummary` printed here must reproduce Finding 81's BOUNDED column to the digit — if they
//! do not, the production path and the bench have diverged and the promotion does not stand.
//!
//! Run: cargo test -p ymir-core --release --test f83_cycles -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, SEA, build_field, pct, sorted};
use std::collections::{HashMap, HashSet};
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::{PrecipParams, precip_mm_per_year};
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, SegmentKind, SegmentRow, Spillway, apply_lake_water_balance,
    below_sea_basin_lakes_infil, c1_drainage_windowed, clip_rivers_to_lakes,
    exorheic_lakes_missing_outlet, potential_evaporation_mm, runoff_km2_to_m3s,
};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{
    D8_DX, D8_DY, DIR_NONE, FlowConfig, RiverSegment, breach_monotone, compute_flow,
};

const DOMAIN_KM: f32 = 400.0;

#[test]
#[ignore]
fn f83_cycles() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 83 · A-guards + C2 on the PRODUCTION field  ==========");

    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;

    let raw = build_field(Knobs::shipped());
    let (w, h) = (raw.width, raw.height);
    let n = w * h;
    let pre0 = c1_drainage_windowed(&raw, None, &dcfg, &ss, DOMAIN_KM);
    let field = breach_monotone(&raw, &pre0.flow.filled, &pre0.lake_map, SEA, w, h);
    let flow = compute_flow(&field, &FlowConfig { sea_level: SEA, ..Default::default() });
    let wc = water_class(&field, SEA);

    for (bed, lat, span) in [("humid", 45.0f32, 40.0f32), ("arid-hot", 25.0, 10.0)] {
        eprintln!("\n╔═══ PRODUCTION (bounded, eps = 0.5 m) / {bed} ═══╗");
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
            (0..n).filter(|&k| field.data[k] > SEA).map(|k| surplus[k] as f64).sum::<f64>() as f32,
        );

        // ── the production tail, exactly as `base_level_read` ran it ──────────
        let pre = c1_drainage_windowed(&raw, None, &dcfg, &ss, DOMAIN_KM);
        let mut dr = c1_drainage_windowed(&field, Some(&dclim), &dcfg, &ss, DOMAIN_KM);
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
        let premerge = dr.lake_map.clone();
        let bs = below_sea_basin_lakes_infil(
            &field,
            &dclim,
            &dcfg,
            &ss,
            DOMAIN_KM,
            Some(&premerge),
            None,
        );
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

        // ── GUARD · the Finding 74 terminal table ─────────────────────────────
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
        let sp_ocean: f64 =
            bs.spillways.iter().filter(|s| end_wc(s) == 1).map(|s| s.discharge_m3s as f64).sum();
        let n_ocean = bs.spillways.iter().filter(|s| end_wc(s) == 1).count();
        let sp_chain: f64 =
            bs.spillways.iter().filter(|s| end_wc(s) == 2).map(|s| s.discharge_m3s as f64).sum();
        let n_chain = bs.spillways.iter().filter(|s| end_wc(s) == 2).count();
        let terminal = wc_q + sp_ocean;
        eprintln!(
            "   GUARD F74: budget {budget:.1} | W->sea {wc_n}/{wc_q:.1} | S->ocean \
             {n_ocean}/{sp_ocean:.1} | S->chained(wc2) {n_chain}/{sp_chain:.1} | TERMINAL \
             {terminal:.1} (x{:.3}) | HOLE {:.1}",
            terminal / budget as f64,
            budget as f64 - terminal
        );

        // ── GUARD · the lake invariants ───────────────────────────────────────
        let mut fc: HashMap<u32, usize> = HashMap::new();
        let mut fmax: HashMap<u32, f32> = HashMap::new();
        for k in 0..n {
            let id = dr.lake_map[k];
            if id != 0 {
                *fc.entry(id).or_default() += 1;
                let e = fmax.entry(id).or_insert(f32::MIN);
                *e = e.max(field.data[k]);
            }
        }
        let mut ids: HashSet<u32> = HashSet::new();
        let (mut dup, mut empty, mut above, mut amis) = (0, 0, 0, 0);
        for lk in &dr.lakes {
            if !ids.insert(lk.base.id) {
                dup += 1;
            }
            let c = fc.get(&lk.base.id).copied().unwrap_or(0);
            if c == 0 {
                empty += 1;
                continue;
            }
            if fmax.get(&lk.base.id).copied().unwrap_or(f32::MIN) > SEA + lk.level_m / n2m + 1e-6 {
                above += 1;
            }
            if (lk.area_km2 - c as f32 * cell_km2).abs() > 0.5 * cell_km2 {
                amis += 1;
            }
        }
        let dangling: usize = {
            let mut v: Vec<u32> =
                dr.lake_map.iter().copied().filter(|&i| i != 0 && !ids.contains(&i)).collect();
            v.sort_unstable();
            v.dedup();
            v.len()
        };
        eprintln!(
            "   GUARD A1: dup {dup} | empty {empty} | footprint ABOVE level {above} | area \
             mismatch {amis} | dangling ids {dangling} | exorheic WITHOUT outlet **{}**",
            exorheic_lakes_missing_outlet(&dr).len()
        );

        let named = dr
            .segment_source_lake
            .iter()
            .enumerate()
            .filter(|(i, s)| dr.segment_kind[*i] == SegmentKind::Spillway && s.is_some())
            .count();
        let unnamed = bs.spillways.len() - named;
        let inv = dr.lakes.iter().filter(|l| l.base.id >= 1_000_001).count();
        let ar = sorted(bs.basins.iter().map(|b| b.area_km2).collect::<Vec<f32>>());
        let inf = sorted(bs.basins.iter().map(|b| b.inflow_m3s).collect::<Vec<f32>>());
        let evp = sorted(bs.basins.iter().map(|b| b.evaporation_m3s).collect::<Vec<f32>>());
        eprintln!(
            "   GUARD A2 BasinSummary {}: area_km2 p50 {:.3} SUM {:.1} | inflow p50 {:.3} SUM \
             {:.2} | evap SUM {:.2} | exorheic {} | INVENTORIED below-sea lakes {inv} ↔ \
             spillways naming a source {named} (unnamed {unnamed})",
            bs.basins.len(),
            pct(&ar, 0.50),
            ar.iter().map(|&x| x as f64).sum::<f64>(),
            pct(&inf, 0.50),
            inf.iter().map(|&x| x as f64).sum::<f64>(),
            evp.iter().map(|&x| x as f64).sum::<f64>(),
            bs.basins.iter().filter(|b| b.exorheic).count()
        );

        if bed != "humid" {
            continue;
        }

        // ══ C2 · the cycle-breaker on the bounded world ══════════════════════
        // Components of the below-sea water (wc == 2), 8-connected. An id whose footprint
        // spans more than one component is a MERGE: the cycle-breaker glued two basins that
        // do not touch, and the spillway of one can then "leave" into the other's label.
        let mut comp = vec![0u32; n];
        let mut next = 0u32;
        let mut stack: Vec<u32> = Vec::new();
        for s0 in 0..n {
            if wc[s0] != 2 || comp[s0] != 0 {
                continue;
            }
            next += 1;
            comp[s0] = next;
            stack.push(s0 as u32);
            while let Some(kk) = stack.pop() {
                let k = kk as usize;
                let (x, y) = ((k % w) as i32, (k / w) as i32);
                for dy in -1i32..=1 {
                    for dx in -1i32..=1 {
                        let (nx, ny) = (x + dx, y + dy);
                        if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                            let nk = ny as usize * w + nx as usize;
                            if wc[nk] == 2 && comp[nk] == 0 {
                                comp[nk] = next;
                                stack.push(nk as u32);
                            }
                        }
                    }
                }
            }
        }
        let mut per: HashMap<u32, HashMap<u32, usize>> = HashMap::new();
        for k in 0..n {
            if wc[k] == 2 && bs.lake_map[k] != 0 {
                *per.entry(bs.lake_map[k]).or_default().entry(comp[k]).or_default() += 1;
            }
        }
        let mut split: Vec<(u32, Vec<usize>)> = per
            .iter()
            .filter(|(_, c)| c.len() > 1)
            .map(|(&id, c)| {
                let mut v: Vec<usize> = c.values().copied().collect();
                v.sort_unstable_by(|a, b| b.cmp(a));
                (id, v)
            })
            .collect();
        split.sort_by_key(|e| e.0);
        let over100 = split
            .iter()
            .filter(|(_, v)| {
                *v.first().unwrap_or(&1) as f64 / *v.last().unwrap_or(&1) as f64 > 100.0
            })
            .count();
        eprintln!(
            "\n   ── C2 · merged ids under the bound ──\n      below-sea regions {next} | ids \
             spanning > 1 component **{}** | of those, size ratio > 100: **{over100}** ({:.1} %)",
            split.len(),
            100.0 * over100 as f64 / split.len().max(1) as f64
        );
        eprintln!(
            "      {:<12} {:>26} {:>10} {:>14}",
            "id", "component cells", "ratio", "Q out m3/s"
        );
        // discharge leaving each merged id, by its own spillway
        let mut q_by_id: HashMap<u32, f64> = HashMap::new();
        for sw in &bs.spillways {
            *q_by_id.entry(sw.lake_id).or_default() += sw.discharge_m3s as f64;
        }
        let mut q_merged = 0.0f64;
        let mut q_merged_over100 = 0.0f64;
        for (id, v) in &split {
            let ratio = *v.first().unwrap_or(&1) as f64 / *v.last().unwrap_or(&1) as f64;
            let q = q_by_id.get(id).copied().unwrap_or(0.0);
            q_merged += q;
            if ratio > 100.0 {
                q_merged_over100 += q;
            }
            eprintln!("      {id:<12} {:>26} {ratio:>10.1} {q:>14.2}", format!("{v:?}"));
        }
        let q_all: f64 = bs.spillways.iter().map(|s| s.discharge_m3s as f64).sum();
        eprintln!(
            "      merged ids carry {q_merged:.2} m3/s of {q_all:.2} total spillway discharge \
             ({:.1} %); those with ratio > 100 carry {q_merged_over100:.2} ({:.1} % of the merged)",
            100.0 * q_merged / q_all.max(1e-9),
            100.0 * q_merged_over100 / q_merged.max(1e-9)
        );

        // ── the dominant terms of the hole, with the receiver named ───────────
        let footprint: HashMap<u32, usize> = {
            let mut m: HashMap<u32, usize> = HashMap::new();
            for k in 0..n {
                if bs.lake_map[k] != 0 {
                    *m.entry(bs.lake_map[k]).or_default() += 1;
                }
            }
            m
        };
        let mut rank: Vec<&Spillway> = bs.spillways.iter().collect();
        rank.sort_by(|a, b| {
            b.discharge_m3s.partial_cmp(&a.discharge_m3s).unwrap_or(std::cmp::Ordering::Equal)
        });
        eprintln!(
            "\n      the six largest spillways (SegmentKind::Spillway for all of them — they are \
             pushed as such at the tail):"
        );
        eprintln!(
            "      {:>12} {:>12} {:>12} {:>12} {:>12} {:>6}",
            "Q m3/s", "from id", "from cells", "-> id", "-> cells", "wc"
        );
        for sw in rank.iter().take(6) {
            let &(lx, ly) = sw.points.last().unwrap();
            let &(sx, sy) = sw.points.first().unwrap();
            let k = ly as usize * w + lx as usize;
            let k0 = sy as usize * w + sx as usize;
            let rid = bs.lake_map[k];
            // THE decisive column: when a spillway lands in its OWN id, does it land in the
            // SAME below-sea body it left, or in the OTHER body the merge glued to it? The
            // first would be a routing bug; the second is the merge artefact.
            let csize = comp.iter().filter(|&&c| c != 0 && c == comp[k]).count();
            let own: Vec<usize> = per
                .get(&sw.lake_id)
                .map(|m| {
                    let mut v: Vec<usize> = m.values().copied().collect();
                    v.sort_unstable_by(|a, b| b.cmp(a));
                    v
                })
                .unwrap_or_default();
            eprintln!(
                "      {:>12.2} {:>12} {:>12} {:>12} {:>12} {:>6}   start comp {} -> lands in                  comp {} of {} cells; the source id's own bodies are {:?}",
                sw.discharge_m3s,
                sw.lake_id,
                footprint.get(&sw.lake_id).copied().unwrap_or(0),
                if rid == 0 { "ocean/none".to_string() } else { rid.to_string() },
                footprint.get(&rid).copied().unwrap_or(0),
                wc[k],
                comp[k0],
                comp[k],
                csize,
                own
            );
        }
        let self_q: f64 = bs
            .spillways
            .iter()
            .filter(|sw| {
                let &(lx, ly) = sw.points.last().unwrap();
                bs.lake_map[ly as usize * w + lx as usize] == sw.lake_id
            })
            .map(|s| s.discharge_m3s as f64)
            .sum();
        let self_n = bs
            .spillways
            .iter()
            .filter(|sw| {
                let &(lx, ly) = sw.points.last().unwrap();
                bs.lake_map[ly as usize * w + lx as usize] == sw.lake_id
            })
            .count();
        eprintln!(
            "      spillways discharging into their OWN id: **{self_n}**, carrying \
             **{self_q:.2} m3/s** = {:.1} % of the chained total {sp_chain:.2} and {:.1} % of the \
             hole {:.1}",
            100.0 * self_q / sp_chain.max(1e-9),
            100.0 * self_q / (budget as f64 - terminal).max(1e-9),
            budget as f64 - terminal
        );
    }
}
