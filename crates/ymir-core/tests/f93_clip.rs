//! ADR 0001 Finding 93 — attribution of `clip_rivers_to_lakes`. **No production change.**
//!
//! The log line `rivers clipped to lakes: 14747 -> 7290 segments` is a NET count: Finding 20 built
//! the clip to *split* each parent into its maximal runs of non-lake points, so a parent can emit
//! zero, one or several children. This counts every outcome separately, by cause, with the drop
//! clauses read from the code:
//!
//! ```text
//! if b - a >= 2 && !drop { out.push((a, b)); }
//! let drop = is_outlet && (endorheic.contains(&src_lake) || src_lake >= 1_000_001);
//! ```
//!
//! Run: cargo test -p ymir-core --release --test f93_clip -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, SEA, build_field};
use std::collections::{HashMap, HashSet};
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, LakeType, SegmentKind, SegmentRow, UnresolvedReason,
    apply_lake_water_balance, below_sea_basin_lakes_infil, c1_drainage_windowed,
    clip_rivers_to_lakes, resolve_exorheic_without_outlet, surface_lake_escape,
    surface_lake_escape_trace,
};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::breach_monotone;

const DOMAIN_KM: f32 = 400.0;

#[test]
#[ignore]
fn f93_clip() {
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 93 · the clip, attributed  ==========");

    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;

    let raw = build_field(Knobs::passes(2));
    let (w, h) = (raw.width, raw.height);
    let n = w * h;
    let d0 = c1_drainage_windowed(&raw, None, &on, &ss, DOMAIN_KM);
    let field = breach_monotone(&raw, &d0.flow.filled, &d0.lake_map, SEA, w, h);
    let m = |g: &ymir_core::grid::GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);

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

    // ══ the PRE-CLIP snapshot, the population every count below is against ══
    let pre_segs = dr.rivers.segments.clone();
    let pre_q = dr.segment_discharge_m3s.clone();
    let lake_map = dr.lake_map.clone();
    let endorheic: HashSet<u32> =
        dr.lakes.iter().filter(|l| l.lake_type == LakeType::Endorheic).map(|l| l.base.id).collect();
    let types: HashMap<u32, LakeType> = dr.lakes.iter().map(|l| (l.base.id, l.lake_type)).collect();
    eprintln!(
        "   pre-clip: **{} segments** · {} lakes ({} Endorheic, {} below-sea) · lake_map covers \
         {:.0} km²",
        pre_segs.len(),
        dr.lakes.len(),
        endorheic.len(),
        dr.lakes.iter().filter(|l| l.base.id >= 1_000_001).count(),
        lake_map.iter().filter(|&&v| v != 0).count() as f32 * cell_km2
    );

    // ══ A — the taxonomy, reproducing the clip's own three clauses ══════════
    let in_lake = |&(x, y): &(u32, u32)| lake_map[y as usize * w + x as usize];
    let (mut kept, mut drop_short, mut drop_endo, mut drop_below) =
        (0usize, 0usize, 0usize, 0usize);
    let (mut pts_total, mut pts_in_lake) = (0usize, 0usize);
    let mut parents_zero = 0usize;
    let mut parents_multi = 0usize;
    let mut q_drop_endo = 0.0f64;
    let mut q_drop_below = 0.0f64;
    let mut q_drop_short = 0.0f64;
    // the false positives that matter: an outlet run dropped below a body with NO emitted spillway
    let spillway_of: HashSet<u32> = bs.spillways.iter().map(|s| s.lake_id).collect();
    let mut fp_no_spillway: HashMap<u32, (usize, f64)> = HashMap::new();
    let mut parent_of: Vec<(usize, usize, usize)> = Vec::new(); // child -> (parent, a, b)
    for (i, s) in pre_segs.iter().enumerate() {
        let (mut a, np) = (0usize, s.points.len());
        pts_total += np;
        pts_in_lake += s.points.iter().filter(|p| in_lake(p) != 0).count();
        let mut runs_here = 0usize;
        while a < np {
            if in_lake(&s.points[a]) != 0 {
                a += 1;
                continue;
            }
            let mut b = a;
            while b < np && in_lake(&s.points[b]) == 0 {
                b += 1;
            }
            let is_outlet = a > 0;
            let src_lake = if is_outlet { in_lake(&s.points[a - 1]) } else { 0 };
            let q = pre_q.get(i).copied().unwrap_or(0.0) as f64;
            let drop_by_lake =
                is_outlet && (endorheic.contains(&src_lake) || src_lake >= 1_000_001);
            if drop_by_lake {
                if src_lake >= 1_000_001 {
                    drop_below += 1;
                    q_drop_below += q;
                    if !spillway_of.contains(&src_lake) {
                        let e = fp_no_spillway.entry(src_lake).or_insert((0, 0.0));
                        e.0 += 1;
                        e.1 += q;
                    }
                } else {
                    drop_endo += 1;
                    q_drop_endo += q;
                }
            } else if b - a < 2 {
                drop_short += 1;
                q_drop_short += q;
            } else {
                kept += 1;
                runs_here += 1;
                // the clip reserves a CONTIGUOUS block of new indices per parent, so emission
                // order is parent order: this is the exact child → parent map, and it replaces
                // the catchment-matching heuristic the first pass used (which was ambiguous —
                // several parents share 2 022 km² — and is withdrawn).
                parent_of.push((i, a, b));
            }
            a = b + 1;
        }
        if runs_here == 0 {
            parents_zero += 1;
        } else if runs_here > 1 {
            parents_multi += 1;
        }
    }
    eprintln!(
        "\n── A · the taxonomy, run by run ──\n   points: {pts_total} total, **{pts_in_lake} \
         ({:.1} %) INSIDE a lake footprint**",
        100.0 * pts_in_lake as f64 / pts_total.max(1) as f64
    );
    let all_runs = kept + drop_short + drop_endo + drop_below;
    for (label, c, q, legit) in [
        ("**KEPT** (a run of ≥ 2 non-lake points)", kept, 0.0, "—"),
        (
            "dropped: run of **< 2 points** (`b - a >= 2`)",
            drop_short,
            q_drop_short,
            "the sub-resolution clause, NOT about sinks",
        ),
        (
            "dropped: outlet of an **ENDORHEIC** lake",
            drop_endo,
            q_drop_endo,
            "legitimate — Finding 20: the water dies in the closed basin",
        ),
        (
            "dropped: outlet of a **BELOW-SEA** body (id ≥ 1 000 001)",
            drop_below,
            q_drop_below,
            "legitimate ONLY if a spillway replaces it — see below",
        ),
    ] {
        eprintln!(
            "      {label}: **{c}** ({:.1} % of {all_runs} runs){} — {legit}",
            100.0 * c as f64 / all_runs.max(1) as f64,
            if q > 0.0 {
                format!(", carrying {q:.1} m³/s of parent discharge")
            } else {
                String::new()
            }
        );
    }
    eprintln!(
        "   parents: **{}** emit NO run (their whole polyline is inside lakes or dropped) · \
         **{parents_multi}** emit MORE than one (the clip SPLITS — this is why the net count is \
         not a deletion count)",
        parents_zero
    );
    eprintln!(
        "\n   **THE FALSE POSITIVES**: outlet runs dropped below a below-sea body that emits NO \
         spillway — the water then has no drawn outlet at all: **{} bodies**, {} runs, \
         **{:.1} m³/s**",
        fp_no_spillway.len(),
        fp_no_spillway.values().map(|v| v.0).sum::<usize>(),
        fp_no_spillway.values().map(|v| v.1).sum::<f64>()
    );
    for (id, (c, q)) in &fp_no_spillway {
        eprintln!(
            "      body **{id}** ({:?}, {:.2} km²): {c} dropped outlet run(s), {q:.1} m³/s of \
             parent discharge, and NO spillway emitted",
            types.get(id),
            dr.lakes.iter().find(|l| l.base.id == *id).map(|l| l.area_km2).unwrap_or(0.0)
        );
    }

    // ══ run the real clip and check the arithmetic ══════════════════════════
    clip_rivers_to_lakes(&mut dr);
    eprintln!(
        "\n   the real clip: **{} → {} segments** · my run census predicted {kept} kept ⇒ {}",
        pre_segs.len(),
        dr.rivers.segments.len(),
        if dr.rivers.segments.len() == kept {
            "**EXACT**".to_string()
        } else {
            format!("off by {}", dr.rivers.segments.len() as i64 - kept as i64)
        }
    );
    assert_eq!(
        dr.rivers.segments.len(),
        kept,
        "A: the run census must reproduce the clip exactly, or the taxonomy is not the clip's"
    );

    // ══ C — the ghost river ═════════════════════════════════════════════════
    eprintln!("\n── C · the ghost: a short reach carrying a whole catchment ──");
    let mut ghosts: Vec<(usize, f32, f32, usize)> = dr
        .rivers
        .segments
        .iter()
        .enumerate()
        .map(|(i, s)| {
            (
                i,
                dr.segment_drainage_km2.get(i).copied().unwrap_or(0.0),
                dr.segment_discharge_m3s.get(i).copied().unwrap_or(0.0),
                s.points.len(),
            )
        })
        .filter(|g| g.3 <= 4)
        .collect();
    ghosts.sort_by(|a, b| b.1.total_cmp(&a.1));
    eprintln!(
        "   post-clip reaches of ≤ 4 points: **{}** · the five largest by catchment:",
        ghosts.len()
    );
    for (i, km2, q, np) in ghosts.iter().take(5) {
        let s = &dr.rivers.segments[*i];
        let &(lx, ly) = s.points.last().unwrap();
        let &(fx, fy) = s.points.first().unwrap();
        let src = dr.segment_source_lake.get(*i).copied().flatten();
        eprintln!(
            "      reach {i}: **{km2:.0} km² · {q:.1} m³/s · {np} points ({:.2} km)** · from \
             ({fx},{fy}) to ({lx},{ly}) · `source_lake` {src:?} · kind {:?} · downstream {:?}",
            (*np as f32 - 1.0) * CELL_KM,
            dr.segment_kind.get(*i),
            s.downstream
        );
    }
    let big: Vec<&(usize, f32, f32, usize)> = ghosts.iter().filter(|g| g.1 > 50_000.0).collect();
    eprintln!(
        "   reaches of ≤ 4 points carrying **more than 50 000 km²**: **{}** — Finding 45 named \
         this in 2024: *\"the 1-segment stub inheriting the parent's area (Finding 42), which \
         lives in `clip_rivers_to_lakes` in core and was deliberately not folded in\"*",
        big.len()
    );
    // and the parents they came from, by the EXACT child → parent map
    for (i, km2, _, _) in ghosts.iter().take(5) {
        let (pi, a, b) = parent_of[*i];
        let ps = &pre_segs[pi];
        eprintln!(
            "      reach {i} ({km2:.0} km²) ⇐ pre-clip parent **{pi}** of {} points ({:.2} km), \
             run [{a}..{b}) ⇒ the clip kept **{:.0} %** of it; points inside lakes in that \
             parent: **{}**",
            ps.points.len(),
            (ps.points.len() as f32 - 1.0) * CELL_KM,
            100.0 * (b - a) as f64 / ps.points.len().max(1) as f64,
            ps.points.iter().filter(|p| in_lake(p) != 0).count()
        );
    }
    // ⚠️ the SIGNIFICATION, which is what makes a 4-point reach read as a 115 000 km² river
    let ratio = 7.5f32;
    let sig = ratio * ratio;
    if let Some((i, km2, q, np)) = ghosts.first() {
        eprintln!(
            "\n   ⚠️ **the geographic-scale signification (ratio {ratio}, area ×{sig})**: reach \
             {i} reads {km2:.0} km² and {q:.1} m³/s in GRID units, i.e. **{:.0} km² and \
             {:.0} m³/s signified**, over {:.2} km of grid = **{:.2} km signified**. The \
             screenshot's \"1 095 m³/s, 115 181 km², 1 km\" is this reach, in signified units — \
             and its `downstream` is **{:?}**, so the link is NOT cut in the data.",
            km2 * sig,
            q * sig,
            (*np as f32 - 1.0) * CELL_KM,
            (*np as f32 - 1.0) * CELL_KM * ratio,
            dr.rivers.segments[*i].downstream
        );
    }

    // ══ D — the lake-trace family, three traces read ════════════════════════
    //
    // ⚠️ THE SPILLWAY APPEND FIRST. `resolve_exorheic_without_outlet` asks whether a REACH starts
    // on the lake's footprint, and in `assemble_hd_drainage` the below-sea spillways are pushed
    // as segments BEFORE it runs. The first pass of this bench omitted the append and reported
    // **9** `Unresolved` lakes with 7 end-of-chain relabels — an artefact of my own chain, not a
    // Finding 92 regression. Reading an invariant off an incomplete chain is reading nothing
    // (rule 12's cousin, and the fourth time this campaign has paid for a missing stage).
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
    eprintln!(
        "\n   (the chain completed: {} spillway segments appended before the Finding 86 test)",
        bs.spillways.len()
    );
    let relabelled = resolve_exorheic_without_outlet(&mut dr);
    eprintln!("\n── D · the lake-trace family, after Finding 92-C ──");
    eprintln!("   end-of-chain relabels: {} {:?}", relabelled.len(), relabelled);
    let unres: Vec<(u32, Option<UnresolvedReason>, f32)> = dr
        .lakes
        .iter()
        .filter(|l| l.lake_type == LakeType::Unresolved)
        .map(|l| (l.base.id, l.unresolved_reason, l.area_km2))
        .collect();
    eprintln!("   `Unresolved` lakes now: **{}** {:?}", unres.len(), unres);
    for (id, reason, area) in &unres {
        let cells: Vec<usize> = (0..n).filter(|&k| dr.lake_map[k] == *id).collect();
        if cells.is_empty() {
            eprintln!("      lake {id}: no footprint in the final map");
            continue;
        }
        match surface_lake_escape(*id, &cells, &field, &dr.lake_map, w, h) {
            Ok((sk, ek)) => {
                let t =
                    surface_lake_escape_trace(*id, &cells, &dr.flow, &field, &dr.lake_map, w, h);
                eprintln!(
                    "      lake **{id}** ({area:.2} km², {:?}): saddle ({},{}) {:.2} m → escape \
                     ({},{}) {:.2} m · the trace from the escape: {}",
                    reason,
                    sk % w,
                    sk / w,
                    m(&field, sk),
                    ek % w,
                    ek / w,
                    m(&field, ek),
                    match t {
                        Ok(k) =>
                            format!("**REACHES** ({},{}) at {:.2} m", k % w, k / w, m(&field, k)),
                        Err(e) => format!("**FAILS: {}**", e.as_str()),
                    }
                );
            }
            Err(e) => eprintln!(
                "      lake **{id}** ({area:.2} km², {:?}): no escape — **{}**",
                reason,
                e.as_str()
            ),
        }
    }
    // 1000001's class: a below-sea body with no emitted spillway
    eprintln!("   below-sea bodies with NO emitted spillway:");
    for l in dr.lakes.iter().filter(|l| l.base.id >= 1_000_001) {
        if spillway_of.contains(&l.base.id) {
            continue;
        }
        let b = bs.basins.iter().find(|b| b.id == l.base.id);
        eprintln!(
            "      body **{}** ({:.2} km², level {:.2} m, {:?}): local inflow {:.2} m³/s, local \
             outflow {:.2} m³/s — its carve reach is dropped by the clip and no spillway replaces it",
            l.base.id,
            l.area_km2,
            l.level_m,
            l.lake_type,
            b.map(|b| b.local_inflow_m3s).unwrap_or(f32::NAN),
            b.map(|b| b.local_outflow_m3s).unwrap_or(f32::NAN)
        );
    }
    eprintln!("\n==========  end Finding 93  ==========\n");
}
