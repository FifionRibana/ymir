//! ADR 0001 Finding 94 blocks A and B — what the inspection panel reads, and whether the
//! microscope lists reaches or systems. **No production change, no viz change.**
//!
//! The panel's six numeric fields all come from one `Watercourse`, so block A is answered by
//! reading the viz (done in the finding); what the bench must measure is the INPUT that decides
//! `length_km`: whether the trunk can climb. The climb is `segs[cur].upstream` with a `pathlen`
//! tie-break, and it is replicated here in twenty lines — the control on the replication is that
//! the terminal count must match the 347 entries the author sees in the list.
//!
//! Run: cargo test -p ymir-core --release --test f94_panel -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, SEA, build_field, pct, sorted};
use std::collections::{HashMap, HashSet};
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, LakeType, SegmentKind, SegmentRow, apply_lake_water_balance,
    below_sea_basin_lakes_infil, c1_drainage_windowed, clip_rivers_to_lakes,
    resolve_exorheic_without_outlet,
};
use ymir_core::terrain::flow::breach_monotone;

const DOMAIN_KM: f32 = 400.0;
const RATIO: f32 = 7.5;

#[test]
#[ignore]
fn f94_panel() {
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 94 · A/B — the panel's input, and the list  ==========");

    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;

    let raw = build_field(Knobs::passes(2));
    let (w, h) = (raw.width, raw.height);
    let n = w * h;
    let d0 = c1_drainage_windowed(&raw, None, &on, &ss, DOMAIN_KM);
    let field = breach_monotone(&raw, &d0.flow.filled, &d0.lake_map, SEA, w, h);
    let climate = c1_climate_placed(&field, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
    let dclim = DrainageClimate {
        precip_internal: &climate.precipitation,
        temperature: &climate.temperature,
    };
    let pre_d = c1_drainage_windowed(&raw, None, &on, &ss, DOMAIN_KM);
    let mut dr = c1_drainage_windowed(&field, Some(&dclim), &on, &ss, DOMAIN_KM);
    dr.lakes = pre_d.lakes;
    dr.lake_map = pre_d.lake_map;
    let li = std::mem::take(&mut dr.lakes);
    dr.lakes = apply_lake_water_balance(
        &field,
        &dr.flow,
        &dclim,
        cell_km2,
        &ss,
        &li,
        &mut dr.lake_map,
        None,
        w,
        h,
    );
    let detected = dr.lake_map.clone();
    let bs =
        below_sea_basin_lakes_infil(&field, &dclim, &on, &ss, DOMAIN_KM, Some(&detected), None);
    let mut d_before: HashMap<u32, usize> = HashMap::new();
    for &id in &dr.lake_map {
        if id != 0 && id < 1_000_001 {
            *d_before.entry(id).or_default() += 1;
        }
    }
    let mut d_after: HashMap<u32, usize> = HashMap::new();
    for k in 0..n {
        if bs.lake_map[k] != 0 {
            dr.lake_map[k] = bs.lake_map[k];
        } else if dr.lake_map[k] != 0 && dr.lake_map[k] < 1_000_001 {
            *d_after.entry(dr.lake_map[k]).or_default() += 1;
        }
    }
    let absorbed: HashSet<u32> =
        d_before.keys().copied().filter(|id| d_after.get(id).copied().unwrap_or(0) == 0).collect();
    dr.lakes.retain(|l| l.base.id >= 1_000_001 || !absorbed.contains(&l.base.id));
    dr.lakes.extend(bs.lakes.iter().cloned());

    // ── the pre-clip upstream links, to measure what the clip does to them ──
    let pre_segs = dr.rivers.segments.clone();
    let pre_with_up = pre_segs.iter().filter(|s| !s.upstream.is_empty()).count();
    let lake_map_at_clip = dr.lake_map.clone();
    let in_lake = |&(x, y): &(u32, u32)| lake_map_at_clip[y as usize * w + x as usize];
    // parents that emit NO run: every link pointing at them is lost by `filter_map(tail)`
    let mut emits_none = vec![false; pre_segs.len()];
    for (i, s) in pre_segs.iter().enumerate() {
        let (mut a, np) = (0usize, s.points.len());
        let mut runs = 0usize;
        while a < np {
            if in_lake(&s.points[a]) != 0 {
                a += 1;
                continue;
            }
            let mut b = a;
            while b < np && in_lake(&s.points[b]) == 0 {
                b += 1;
            }
            if b - a >= 2 {
                runs += 1;
            }
            a = b + 1;
        }
        emits_none[i] = runs == 0;
    }
    let links_total: usize = pre_segs.iter().map(|s| s.upstream.len()).sum();
    let links_to_swallowed: usize = pre_segs
        .iter()
        .map(|s| s.upstream.iter().filter(|&&u| u < emits_none.len() && emits_none[u]).count())
        .sum();
    eprintln!(
        "   pre-clip: {} segments, {} with a non-empty `upstream`, **{} upstream LINKS** · of \
         those, **{} ({:.1} %) point at a parent that emits NO run** and are dropped by \
         `filter_map(tail)`",
        pre_segs.len(),
        pre_with_up,
        links_total,
        links_to_swallowed,
        100.0 * links_to_swallowed as f64 / links_total.max(1) as f64
    );

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
                let kk = ly as usize * w + lx as usize;
                dr.flow.accumulation.data.get(kk).copied().unwrap_or(0.0)
            },
            discharge_profile_m3s: vec![sw.discharge_m3s; sw.points.len()],
            kind: SegmentKind::Spillway,
            source_lake: inventoried.contains(&sw.lake_id).then_some(sw.lake_id),
        });
    }
    let _ = resolve_exorheic_without_outlet(&mut dr);

    // ── B · the microscope's population: terminals, and can the trunk climb? ──
    let segs = &dr.rivers.segments;
    let ns = segs.len();
    let with_up = segs.iter().filter(|s| !s.upstream.is_empty()).count();
    let terminals: Vec<usize> = (0..ns).filter(|&i| segs[i].downstream.is_none()).collect();
    eprintln!(
        "\n── B · post-clip: {ns} segments · **{with_up} ({:.1} %) have a non-empty `upstream`** · \
         **{} terminals (`downstream == None`)** — the author's list shows 347 entries",
        100.0 * with_up as f64 / ns.max(1) as f64,
        terminals.len()
    );

    // `pathlen`, the viz's own tie-break, replicated
    let mut pathlen = vec![0u32; ns];
    {
        let mut order: Vec<usize> = (0..ns).collect();
        order.sort_by_key(|&i| segs[i].strahler_order);
        for _ in 0..16 {
            let mut changed = false;
            for &i in &order {
                let up = segs[i]
                    .upstream
                    .iter()
                    .copied()
                    .filter(|&u| u < ns)
                    .map(|u| pathlen[u])
                    .max()
                    .unwrap_or(0);
                let v = up + segs[i].points.len() as u32;
                if v != pathlen[i] {
                    pathlen[i] = v;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
    }
    let wc = water_class(&field, SEA);
    // ── Finding 45's chaining, replicated: an inflow reach is linked to the outlet reach of the
    //    EXORHEIC lake it dies on. Without it the entry count is an upper bound, which is why the
    //    first pass of this bench could not use the author's 347 as a control.
    let exo_ids: HashSet<u32> =
        dr.lakes.iter().filter(|l| l.lake_type == LakeType::Exorheic).map(|l| l.base.id).collect();
    // ⚠️ `classify_sink` scans the 3×3 NEIGHBOURHOOD for a lake id, not the mouth cell: the clip
    // keeps runs of NON-lake points, so every reach ending at a lake stops ONE CELL SHORT of it.
    // The first two passes of this bench probed the mouth cell only and reported 2 474 entries
    // against the author's 347 — my replication, not the viz. This is the corrected probe.
    let lake_near = |x: u32, y: u32| -> u32 {
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                    continue;
                }
                let id = dr.lake_map[ny as usize * w + nx as usize];
                if id != 0 {
                    return id;
                }
            }
        }
        0
    };
    // a lake's outlet reach: a segment whose FIRST cell borders it, highest discharge wins
    let mut outlet_of: HashMap<u32, usize> = HashMap::new();
    for (i, sg) in segs.iter().enumerate() {
        let (sx, sy) = sg.points[0];
        let id0 = lake_near(sx, sy);
        let touch = if exo_ids.contains(&id0) { id0 } else { 0 };
        if touch != 0 {
            let q = |j: usize| dr.segment_discharge_m3s.get(j).copied().unwrap_or(0.0);
            if outlet_of.get(&touch).is_none_or(|&j| q(i) > q(j)) {
                outlet_of.insert(touch, i);
            }
        }
    }
    let across_lake = |i: usize| -> Option<usize> {
        let &(mx, my) = segs[i].points.last()?;
        let id = lake_near(mx, my);
        if id == 0 || !exo_ids.contains(&id) {
            return None;
        }
        outlet_of.get(&id).copied().filter(|&o| o != i)
    };
    // the ROOT with chaining: follow downstream, then across an exorheic lake
    let mut root = vec![0usize; ns];
    for i in 0..ns {
        let mut j = i;
        for _ in 0..=ns {
            match segs[j].downstream {
                Some(k) if k < ns => j = k,
                _ => match across_lake(j) {
                    Some(o) => j = o,
                    None => break,
                },
            }
        }
        root[i] = j;
    }
    let mut chained_groups: HashMap<usize, usize> = HashMap::new();
    for i in 0..ns {
        *chained_groups.entry(root[i]).or_default() += 1;
    }
    eprintln!(
        "   **Finding 45's chaining, replicated: {} ENTRIES** (against {} raw terminals) — the \
         author's list shows **347**, so this replication is {} and the counts below are read on it",
        chained_groups.len(),
        terminals.len(),
        if (chained_groups.len() as i64 - 347).abs() <= 40 {
            "**within 40 of it**"
        } else {
            "**NOT the same population — read the gap as the finding**"
        }
    );
    // the trunk climb, with Finding 45's symmetric upstream bridge
    let inflow_to_lake_of = |i: usize, pathlen: &[u32]| -> Option<usize> {
        let (sx, sy) = segs[i].points[0];
        let id = lake_near(sx, sy);
        if id == 0 || !exo_ids.contains(&id) {
            return None;
        }
        segs.iter()
            .enumerate()
            .filter(|&(j, sg)| {
                j != i
                    && sg.downstream.is_none()
                    && sg.points.last().is_some_and(|&(mx, my)| lake_near(mx, my) == id)
            })
            .max_by_key(|&(j, _)| pathlen[j])
            .map(|(j, _)| j)
    };
    let mut trunk_reaches: Vec<f32> = Vec::new();
    let mut trunk_km: Vec<f32> = Vec::new();
    let mut single = 0usize;
    let mut to_sea = 0usize;
    let mut rows: Vec<(usize, f32, f32, usize, f32)> = Vec::new(); // root, km2, q, reaches, km
    // ⚠️ POPULATION. The entries the panel lists are the chained ROOTS (401), not the 2 477 raw
    // terminals — a terminal that gets bridged across an exorheic lake is not an entry. The first
    // three passes of this bench computed the single-reach share over the terminals and reported
    // 81.8 %, which is a statement about fragments, not about what the user sees.
    let roots: Vec<usize> = {
        let mut v: Vec<usize> = chained_groups.keys().copied().collect();
        v.sort_unstable();
        v
    };
    for &r in &roots {
        let mut trunk = vec![r];
        let mut cur = r;
        for _ in 0..=ns {
            match segs[cur].upstream.iter().copied().filter(|&u| u < ns).max_by_key(|&u| pathlen[u])
            {
                Some(next) => {
                    cur = next;
                    trunk.push(cur);
                }
                // Finding 45's symmetric bridge: at the head of an outlet reach, climb the
                // longest reach flowing INTO the same exorheic lake.
                None => match inflow_to_lake_of(cur, &pathlen) {
                    Some(up) if !trunk.contains(&up) => {
                        cur = up;
                        trunk.push(cur);
                    }
                    _ => break,
                },
            }
        }
        let pts: usize = trunk.iter().map(|&s| segs[s].points.len()).sum();
        let km = (pts as f32 - 1.0).max(0.0) * CELL_KM;
        trunk_reaches.push(trunk.len() as f32);
        trunk_km.push(km);
        if trunk.len() == 1 {
            single += 1;
        }
        let &(mx, my) = segs[r].points.last().unwrap();
        if wc[my as usize * w + mx as usize] == 1 {
            to_sea += 1;
        }
        rows.push((
            r,
            dr.segment_drainage_km2.get(r).copied().unwrap_or(0.0),
            dr.segment_discharge_m3s.get(r).copied().unwrap_or(0.0),
            trunk.len(),
            km,
        ));
    }
    let tr = sorted(trunk_reaches.clone());
    let tk = sorted(trunk_km.clone());
    eprintln!(
        "   the TRUNK CLIMB, replicated: reaches per trunk p50 **{:.0}** p90 {:.0} max {:.0} · \
         trunk length km p50 **{:.2}** p90 {:.1} max {:.1}",
        pct(&tr, 0.50),
        pct(&tr, 0.90),
        tr.last().copied().unwrap_or(0.0),
        pct(&tk, 0.50),
        pct(&tk, 0.90),
        tk.last().copied().unwrap_or(0.0)
    );
    eprintln!(
        "   **{single} of {} ENTRIES ({:.1} %) have a trunk of ONE reach** — for those, \
         `length_km` IS the mouth reach's length while `catchment_km2` is its (inherited) \
         catchment, which is the panel's \"1 km, 115 181 km²\"",
        roots.len(),
        100.0 * single as f64 / roots.len().max(1) as f64
    );
    eprintln!(
        "   of the {} terminals, **{to_sea} end on the OCEAN** (`water_class == 1`) and **{} end \
         on a lake or an enclosed cell** — the latter are the inter-water-body fragments Finding 45
         chains only across EXORHEIC lakes",
        terminals.len(),
        terminals.len() - to_sea
    );

    // the top entries by signified catchment, with their trunk
    rows.sort_by(|a, b| b.1.total_cmp(&a.1));
    let sig = RATIO * RATIO;
    eprintln!(
        "\n   the eight largest entries by catchment (signified ×{sig:.2}, as the panel shows them):"
    );
    for (r, km2, q, nr, km) in rows.iter().take(8) {
        let &(mx, my) = segs[*r].points.last().unwrap();
        eprintln!(
            "      root {r:>5}: **{:.0} km² · {:.0} m³/s · {:.1} km** (signified) — trunk of \
             **{nr} reach(es)**, mouth ({mx},{my}) on {} · kind {:?}",
            km2 * sig,
            q * sig,
            km * RATIO,
            match wc[my as usize * w + mx as usize] {
                1 => "the OCEAN",
                2 => "an ENCLOSED below-sea cell",
                _ => "LAND (a lake shore)",
            },
            dr.segment_kind.get(*r)
        );
    }
    // the true source-to-sea systems, if the chain could cross every water body
    // ⚠️ the first pass keyed this by LAKE ID, so it printed 6 where it meant terminals. Count
    // TERMINALS, and report the distinct lakes separately.
    let mut on_exo = 0usize;
    let mut on_other = 0usize;
    let mut distinct: HashSet<u32> = HashSet::new();
    for &r in &terminals {
        let &(mx, my) = segs[r].points.last().unwrap();
        let id = lake_near(mx, my);
        if id == 0 {
            continue;
        }
        distinct.insert(id);
        if exo_ids.contains(&id) {
            on_exo += 1;
        } else {
            on_other += 1;
        }
    }
    eprintln!(
        "\n   TERMINALS dying on a lake footprint: **{}** over **{}** distinct lakes · on an \
         EXORHEIC lake (the only ones Finding 45 bridges): **{on_exo}** · on an \
         Endorheic/Unresolved/below-sea body (never bridged): **{on_other}**",
        on_exo + on_other,
        distinct.len()
    );
    // ══ D — the panel's LAKE claims, as a population ════════════════════════
    //
    // ⚠️ The screenshot's "#341", "#298", "#381" are MICROSCOPE LIST INDICES, not segment ids —
    // the list is sorted by discharge and re-indexed every run, so a bench cannot address them.
    // #381 was identified by its NUMBERS (115 211 km², 1 095 m³/s), not by its index, and the same
    // discipline applies here: D checks the CLAIM the panel makes for every lake instead of three
    // indices that do not survive a re-run.
    eprintln!("\n── D · the panel's `Exutoire` claim, checked for every lake ──");
    let mut claims_outlet = 0usize;
    let mut has_reach = 0usize;
    let mut mismatch: Vec<(u32, f32, LakeType)> = Vec::new();
    for l in &dr.lakes {
        if l.lake_type == LakeType::Endorheic {
            continue; // the panel prints "aucun (fermé)" — no claim to check
        }
        claims_outlet += 1;
        let found = segs.iter().enumerate().any(|(i, sg)| {
            let (sx, sy) = sg.points[0];
            lake_near(sx, sy) == l.base.id && dr.segment_kind.get(i) != Some(&SegmentKind::Spillway)
        });
        let spill = bs.spillways.iter().any(|sw| sw.lake_id == l.base.id);
        if found {
            has_reach += 1;
        } else if !spill {
            mismatch.push((l.base.id, l.area_km2, l.lake_type));
        }
    }
    eprintln!(
        "   lakes the panel would call \"oui → aval\": **{claims_outlet}** · with a real outlet \
         REACH: **{has_reach}** · with neither a reach nor a spillway: **{}** ⇒ that many claims \
         the data does not back",
        mismatch.len()
    );
    for (id, a, t) in mismatch.iter().take(8) {
        eprintln!("      lake **{id}** ({a:.2} km², {t:?}) — claims an outlet, has none");
    }
    eprintln!(
        "   ⚠️ the \"évaporatif\" label: Finding 92-D replaced it with \"cuvette sous-marine SANS \
         exutoire tracé\" in commit e60fa54. A screenshot still showing \"puits sous-marin \
         (évaporatif)\" was taken on a binary that predates it — nothing to re-fix."
    );
    eprintln!("\n==========  end Finding 94 · A/B/D  ==========\n");
}
