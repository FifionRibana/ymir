//! ADR 0001 Finding 88 — evaporation in local terms (A), where the sixth term's water actually
//! goes (B), the Finding 37c port measured (C), the 623.6 m lake's incision history (D1–D4), the
//! `DIR_NONE` of basin 1000001 read on the right field (E), and why 4 850 coastal cells with real
//! catchments carry no segment (F).
//!
//! **Rule 12** (new this round): every block that walks `flow.direction` declares its field and
//! proves by hash that it is the one production routes on, BEFORE reading a reason off it. The
//! declaration is `common::declared_flow`; the failure it guards against cost Finding 87-D three
//! runs and fifteen wrong reasons.
//!
//! Run: cargo test -p ymir-core --release --test f88_routing -- --ignored --nocapture

mod common;

use common::{
    CELL_KM, Knobs, SEA, build_field, declared_flow, flow_field_hash, pct,
    production_k_and_edifices, sorted,
};
use std::collections::{HashMap, HashSet};
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::{PrecipParams, precip_mm_per_year};
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, LakeType, SegmentKind, SegmentRow, UnresolvedReason,
    apply_lake_water_balance, below_sea_basin_lakes_infil, c1_drainage_windowed,
    clip_rivers_to_lakes, potential_evaporation_mm, resolve_exorheic_without_outlet,
    runoff_accumulation, runoff_km2_to_m3s, surface_lake_escape,
};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{
    D8_DX, D8_DY, DIR_NONE, FlowConfig, FlowResult, RiverSegment, breach_monotone, compute_flow,
};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
/// The lake block B and D are about: Finding 87-B's 623.6 m body, by its FLOOR CELL and not by its
/// id (Finding 81's trap — ids are scan-order).
const DEEP_FLOOR: (usize, usize) = (3487, 5902);

/// Follow `flow.direction` to wherever it ends. Returns `(terminus, steps, cycled)`.
///
/// This is NOT a trace that can fail: it answers "where does production put this cell's water",
/// which is the question block B asks and which `walk_to_sink` deliberately refuses to answer
/// (it reports a FAILURE when the descent re-enters the lake, because for a LABEL that is the
/// honest answer). Both are needed and they are different questions.
fn terminus(
    start: usize,
    flow: &FlowResult,
    field: &GridF32,
    w: usize,
    h: usize,
) -> (usize, usize, bool) {
    let mut k = start;
    let mut seen: HashSet<usize> = HashSet::new();
    for step in 0..200_000usize {
        if field.data[k] <= SEA {
            return (k, step, false);
        }
        let d = flow.direction[k];
        if d == DIR_NONE {
            return (k, step, false);
        }
        if !seen.insert(k) {
            return (k, step, true);
        }
        let nx = ((k % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
        let ny = ((k / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
        k = ny * w + nx;
    }
    (k, 200_000, true)
}

fn median(v: &[f32]) -> f32 {
    let s = sorted(v.to_vec());
    pct(&s, 0.50)
}

/// Max-of-8 gradient in degrees at `k`, in the delivered field.
///
/// ⚠️ `n2m` is NOT optional. The first version of this function divided a NORM height difference
/// by a distance in METRES, i.e. it was 11 300× too small, and block D4 consequently found **zero**
/// bodies with a rim above 30° — including the one Finding 87-B had measured at 42.8°. A unit
/// mismatch inside an instrument reads exactly like a result about the world.
fn slope_deg(field: &GridF32, k: usize, n2m: f32, w: usize, h: usize) -> f32 {
    let (x, y) = ((k % w) as i32, (k / w) as i32);
    let mut best = 0.0f32;
    for d in 0..8 {
        let nx = (x + D8_DX[d]).rem_euclid(w as i32) as usize;
        let ny = (y + D8_DY[d]).rem_euclid(h as i32) as usize;
        let diag = D8_DX[d] != 0 && D8_DY[d] != 0;
        let dist = if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M };
        let dz = (field.data[k] - field.data[ny * w + nx]).abs() * n2m;
        best = best.max(dz / dist);
    }
    best.atan().to_degrees()
}

#[test]
#[ignore]
fn f88_routing() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 88 · A · B · C · D · E · F  ==========");

    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    let mut off = on.clone();
    off.merged_union_relevel = None;
    let head_km2 = on.thresholds.head_km2;
    let stream_km2 = on.thresholds.stream_km2;
    eprintln!(
        "   thresholds: head {head_km2} km² · **stream (the EXPORT clip) {stream_km2} km²** · \
         cell {CELL_M:.2} m · 1 norm = {n2m:.1} m"
    );

    let raw = build_field(Knobs::shipped());
    let (w, h) = (raw.width, raw.height);
    let n = w * h;
    let pre0 = c1_drainage_windowed(&raw, None, &on, &ss, DOMAIN_KM);
    let field = breach_monotone(&raw, &pre0.flow.filled, &pre0.lake_map, SEA, w, h);
    let flow_naked = compute_flow(&field, &FlowConfig { sea_level: SEA, ..Default::default() });
    let wc = water_class(&field, SEA);
    let land: Vec<bool> = (0..n).map(|k| field.data[k] > SEA).collect();
    // D needs the field before the incision, and after ONE of its two passes
    let pre = build_field(Knobs::no_incision());
    let it1 = build_field(Knobs::passes(1));
    let (k_field, edifices) = production_k_and_edifices();

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

            // ══ RULE 12 — declare the field, and compare it to the bench's own ══
            let hp = flow_field_hash(&dr.flow);
            let hn = flow_field_hash(&flow_naked);
            let diff_dir =
                (0..n).filter(|&k| dr.flow.direction[k] != flow_naked.direction[k]).count();
            let none_p = (0..n).filter(|&k| dr.flow.direction[k] == DIR_NONE).count();
            let none_n = (0..n).filter(|&k| flow_naked.direction[k] == DIR_NONE).count();
            eprintln!(
                "      RULE 12: production `dr.flow` 0x{hp:016x} | bench `compute_flow` \
                 0x{hn:016x} | **{} D8 directions differ ({:.3} % of cells)** | DIR_NONE: \
                 production {none_p} vs bench {none_n}",
                diff_dir,
                100.0 * diff_dir as f64 / n as f64
            );
            // rule 12: the declaration is the line, not the binding — holding a
            // reference across the chain would borrow `dr` for the rest of the block.
            let _ = declared_flow("dr.flow", &dr.flow, &dr.flow);

            dr.lakes = pre_d.lakes;
            dr.lake_map = pre_d.lake_map;
            let lakes_in = std::mem::take(&mut dr.lakes);
            // How many of these the OLD code labelled Exorheic from the balance alone: that is
            // block C's "before", and it is knowable without a gate because the old line was
            // unconditional (`a_eq >= a_sill` ⇒ Exorheic).
            let n_below_sea_in = lakes_in.iter().filter(|l| l.base.id >= 1_000_001).count();
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
            let surf_map = dr.lake_map.clone();
            // production's own runoff accumulation, on production's own field
            let racc = runoff_accumulation(&field, &dr.flow, &dclim, cell_km2, None, None, w, h);
            // the Finding 87 reading, on the bench field, kept for the comparison
            let racc_bench =
                runoff_accumulation(&field, &flow_naked, &dclim, cell_km2, None, None, w, h);

            // ══ C — the Finding 37c port, measured ═══════════════════════════
            let mut unres: Vec<(u32, f32, f32, UnresolvedReason)> = Vec::new();
            let mut by_reason: HashMap<String, usize> = HashMap::new();
            let n_exo_balance = dr
                .lakes
                .iter()
                .filter(|l| {
                    l.base.id < 1_000_001
                        && matches!(l.lake_type, LakeType::Exorheic | LakeType::Unresolved)
                })
                .count();
            for l in dr.lakes.iter().filter(|l| l.base.id < 1_000_001) {
                if l.lake_type != LakeType::Unresolved {
                    continue;
                }
                let cells: Vec<usize> = (0..n).filter(|&k| surf_map[k] == l.base.id).collect();
                let inflow = cells.iter().map(|&k| racc[k]).fold(0.0f32, f32::max);
                let why = l.unresolved_reason.expect("Finding 88-C: the reason travels with it");
                *by_reason.entry(format!("{why:?}")).or_default() += 1;
                unres.push((l.base.id, l.area_km2, runoff_km2_to_m3s(inflow), why));
            }
            eprintln!(
                "      C: detected lakes the BALANCE calls exorheic: **{n_exo_balance}** → after \
                 the F37c escape trace, **{}** remain Unresolved ({} resolved by the port)",
                unres.len(),
                n_exo_balance - unres.len()
            );
            let mut rs: Vec<(&String, &usize)> = by_reason.iter().collect();
            rs.sort_by(|a, b| b.1.cmp(a.1));
            for (why, c) in rs {
                eprintln!("         {c:>3} × {why}");
            }
            for (id, a, q, why) in unres.iter().take(8) {
                eprintln!("         lake {id} — {a:.3} km² · local inflow {q:.2} m³/s · {:?}", why);
            }
            eprintln!(
                "      C guard: below-sea lakes reaching this function: **{n_below_sea_in}** \
                 (must be 0 — the below-sea labels are set downstream, so none can change here)"
            );
            assert_eq!(
                n_below_sea_in, 0,
                "C stop rule: a below-sea lake reached apply_lake_water_balance"
            );

            // ══ B — where the sixth term's water actually goes ═══════════════
            // On production's field, from the FOOTPRINT, allowing re-entry: this is routing,
            // not labelling.
            let mut cls: HashMap<&str, (usize, f64)> = HashMap::new();
            let mut brows: Vec<(u32, f32, f32, &str, usize)> = Vec::new();
            for (id, a, q, _) in &unres {
                let cells: Vec<usize> = (0..n).filter(|&k| surf_map[k] == *id).collect();
                let deepest = cells
                    .iter()
                    .copied()
                    .min_by(|&x, &y| field.data[x].total_cmp(&field.data[y]))
                    .unwrap_or(0);
                let (term, steps, cycled) = terminus(deepest, &dr.flow, &field, w, h);
                let klass = if field.data[term] <= SEA {
                    "the SEA (already in the coastal term ⇒ DOUBLE COUNTED)"
                } else if surf_map[term] == *id {
                    "inside its OWN footprint (the water stands ⇒ term is REAL)"
                } else if surf_map[term] != 0 {
                    "another SURFACE lake (a chain)"
                } else if cycled {
                    "a cycle outside its footprint"
                } else {
                    "a land pit with no water body (DIR_NONE)"
                };
                let e = cls.entry(klass).or_insert((0, 0.0));
                e.0 += 1;
                e.1 += *q as f64;
                brows.push((*id, *a, *q, klass, steps));
            }
            eprintln!(
                "      B: the {} Unresolved lakes, routed on production's field:",
                unres.len()
            );
            let mut cv: Vec<(&&str, &(usize, f64))> = cls.iter().collect();
            cv.sort_by(|a, b| b.1.0.cmp(&a.1.0));
            for (klass, (c, q)) in cv {
                eprintln!("         {c:>3} lakes · {q:7.1} m³/s → {klass}");
            }
            brows.sort_by(|a, b| b.2.total_cmp(&a.2));
            for (id, a, q, klass, steps) in brows.iter().take(6) {
                eprintln!(
                    "         lake {id} — {a:.1} km² · {q:.1} m³/s · {steps} steps → {klass}"
                );
            }
            let q_real: f64 =
                brows.iter().filter(|r| r.3.contains("REAL")).map(|r| r.2 as f64).sum();
            let q_double: f64 =
                brows.iter().filter(|r| r.3.contains("DOUBLE")).map(|r| r.2 as f64).sum();

            // ══ the chain, to get the below-sea system and the spillways ═════
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
            let bs_types_before: HashMap<u32, LakeType> =
                bs.lakes.iter().map(|l| (l.base.id, l.lake_type)).collect();
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
            // C stop rule: no lake ALREADY labelled by the below-sea path may change
            let changed_below: Vec<u32> = dr
                .lakes
                .iter()
                .filter(|l| l.base.id >= 1_000_001)
                .filter(|l| bs_types_before.get(&l.base.id).is_some_and(|b| *b != l.lake_type))
                .map(|l| l.base.id)
                .collect();
            eprintln!(
                "      C stop rule: below-sea lakes whose label changed after the port: **{}** \
                 {:?} (the round allows at most 2) | end-of-chain relabels {} | passes {} \
                 (converged {})",
                changed_below.len(),
                changed_below,
                relabelled.len(),
                t.passes_used,
                !t.not_converged
            );
            assert!(
                changed_below.len() <= 2,
                "C stop rule fired: {} below-sea lakes changed label",
                changed_below.len()
            );

            // ══ A — evaporation in LOCAL terms ══════════════════════════════
            let local_in: f64 = bs.basins.iter().map(|b| b.local_inflow_m3s as f64).sum();
            let local_out: f64 = bs.basins.iter().map(|b| b.local_outflow_m3s as f64).sum();
            let evap_chained: f64 = bs.basins.iter().map(|b| b.evaporation_m3s as f64).sum();
            let evap_local: f64 = bs.basins.iter().map(|b| b.local_evaporation_m3s as f64).sum();
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
            let close_chained = 100.0 * (local_out + evap_chained) / local_in.max(1e-9);
            let close_local = 100.0 * (local_out + evap_local) / local_in.max(1e-9);
            eprintln!(
                "      A: below-sea local in {local_in:.1} · local out {local_out:.1} · evap \
                 CHAINED {evap_chained:.2} → **{close_chained:.1} %** · evap LOCAL \
                 {evap_local:.2} → **{close_local:.1} %** · chained excess \
                 **{:.2} m³/s** · retained {retained:.1}",
                evap_chained - evap_local
            );
            eprintln!(
                "      ⚠️ A, rule 10: `{close_local:.1} %` is an **IDENTITY**, not a measurement — \
                 local_out + local_evap == local_in per basin by algebra, in both branches. It is \
                 DEDUCED. What is measured here is the chained excess above, and the partition \
                 out/evap = {:.1} / {:.1} %.",
                100.0 * local_out / local_in.max(1e-9),
                100.0 * evap_local / local_in.max(1e-9)
            );
            assert!(
                (close_local - 100.0).abs() < 0.2,
                "A: the identity must hold to the digit, got {close_local}"
            );

            // ══ F/B — the coastal arrival, on production's field ════════════
            let on_segment: HashSet<usize> = dr
                .rivers
                .segments
                .iter()
                .flat_map(|s| s.points.iter())
                .map(|&(x, y)| y as usize * w + x as usize)
                .collect();
            let mut arrives = 0.0f64;
            let mut arrives_bench = 0.0f64;
            let (mut ungauged, mut ung_q) = (Vec::<(f32, usize)>::new(), 0.0f64);
            for k in 0..n {
                if !land[k] || dr.flow.direction[k] == DIR_NONE {
                    continue;
                }
                let d = dr.flow.direction[k] as usize;
                let nx = ((k % w) as i32 + D8_DX[d]).rem_euclid(w as i32) as usize;
                let ny = ((k / w) as i32 + D8_DY[d]).rem_euclid(h as i32) as usize;
                if wc[ny * w + nx] != 1 {
                    continue;
                }
                arrives += runoff_km2_to_m3s(racc[k]) as f64;
                arrives_bench += runoff_km2_to_m3s(racc_bench[k]) as f64;
                if !on_segment.contains(&k) {
                    ung_q += runoff_km2_to_m3s(racc[k]) as f64;
                    ungauged.push((dr.flow.accumulation.data[k] * cell_km2, k));
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
            // ADR Finding 88-B — the ocean term read at the END OF THE CHAIN, not at the first
            // hop. A basin whose spillway lands in ANOTHER below-sea body is credited to neither
            // `sp_ocean_local` nor `retained`, so its local outflow leaves the budget entirely.
            // Follow `chained_into`-style hops on the spillway map until the destination is the
            // ocean, a body with no spillway (retained), or a cycle.
            let mut next_hop: HashMap<u32, Option<u32>> = HashMap::new();
            for sw in &bs.spillways {
                let &(lx, ly) = sw.points.last().unwrap();
                let kk = ly as usize * w + lx as usize;
                // ⚠️ population: `bs.lake_map` holds ONLY the below-sea bodies. A spillway
                // landing on a DETECTED surface lake reads as id 0 there and looked like
                // "nowhere identifiable" on the first run — 117.6 m³/s mislabelled by my own
                // classifier. `premerge` is the detected map at this point in the chain.
                let dest = if wc[kk] == 1 {
                    None // the ocean
                } else if bs.lake_map[kk] != 0 {
                    Some(bs.lake_map[kk])
                } else if premerge[kk] != 0 {
                    Some(u32::MAX - 1) // a DETECTED surface lake
                } else {
                    Some(u32::MAX) // nowhere identifiable
                };
                next_hop.entry(sw.lake_id).or_insert(dest);
            }
            let chain_end = |start: u32| -> &'static str {
                let mut cur = start;
                let mut seen: HashSet<u32> = HashSet::new();
                for _ in 0..64 {
                    if !seen.insert(cur) {
                        return "a CYCLE in the chain";
                    }
                    match next_hop.get(&cur) {
                        None => return "retained (no spillway at the end)",
                        Some(None) => return "the OCEAN",
                        Some(Some(nx)) if *nx == u32::MAX => return "nowhere identifiable",
                        // block B measured it: a detected lake's surface routes to the sea
                        Some(Some(nx)) if *nx == u32::MAX - 1 => {
                            return "a DETECTED LAKE, whose surface routes to the SEA (B)";
                        }
                        Some(Some(nx)) if *nx == cur => return "its OWN body",
                        Some(Some(nx)) => cur = *nx,
                    }
                }
                "over 64 hops"
            };
            let mut by_end: HashMap<&str, f64> = HashMap::new();
            for b in &bs.basins {
                *by_end.entry(chain_end(b.id)).or_default() += b.local_outflow_m3s as f64;
            }
            let ocean_transitive = by_end.get("the OCEAN").copied().unwrap_or(0.0)
                + by_end
                    .get("a DETECTED LAKE, whose surface routes to the SEA (B)")
                    .copied()
                    .unwrap_or(0.0);
            let retained_transitive =
                by_end.get("retained (no spillway at the end)").copied().unwrap_or(0.0);
            eprintln!(
                "      B termination counters (Finding 84's own, on this run): to_ocean {} · \
                 to_other_body {} · to_detected_lake {} · to_own_body {} · **to_nothing {}** · \
                 q_own_body {:.1} · q_detected_lake {:.1} · q_routed_through_lake {:.1} · \
                 q_lake_unrouted {:.1}",
                t.to_ocean,
                t.to_other_body,
                t.to_detected_lake,
                t.to_own_body,
                t.to_nothing,
                t.q_own_body_m3s,
                t.q_detected_lake_m3s,
                t.q_routed_through_lake_m3s,
                t.q_lake_unrouted_m3s
            );
            eprintln!("      B chain, where each basin's LOCAL outflow actually ends:");
            let mut bev: Vec<(&&str, &f64)> = by_end.iter().collect();
            bev.sort_by(|a, b| b.1.total_cmp(a.1));
            for (k, v) in bev {
                eprintln!("         {v:7.1} m³/s → {k}");
            }
            let q_all: f64 = unres.iter().map(|u| u.2 as f64).sum();
            let six_f87 = arrives + sp_ocean_local + evap_chained + evap_surface + retained + q_all;
            let six_f88 = arrives + sp_ocean_local + evap_local + evap_surface + retained + q_real;
            let six_f88c = arrives
                + ocean_transitive
                + evap_local
                + evap_surface
                + retained_transitive
                + q_real;
            eprintln!(
                "      B/F coastal arrival: production field **{arrives:.1}** vs the Finding 87 \
                 bench field {arrives_bench:.1} (Δ {:.1})",
                arrives - arrives_bench
            );
            eprintln!(
                "      ⚠️ B, the PARTITION identity: coastal arrival {arrives:.1} + Σ below-sea \
                 LOCAL inflow {local_in:.1} = **{:.1}** against a budget of {budget:.1} \
                 (Δ {:.2}, {:.2} %). Every drop of land surplus either crosses the coast or enters \
                 a below-sea basin — so a budget built from these two DOES NOT MEASURE \
                 CONSERVATION ALONG THE ROUTING, it re-states the accumulation's own partition.",
                arrives + local_in,
                arrives + local_in - budget as f64,
                100.0 * (arrives + local_in) / budget as f64
            );
            eprintln!(
                "      F87 SIX TERMS re-run: {arrives:.1} + {sp_ocean_local:.1} + \
                 {evap_chained:.2} + {evap_surface:.2} + {retained:.1} + {q_all:.1} = \
                 {six_f87:.1} = **{:.1} %**",
                100.0 * six_f87 / budget as f64
            );
            eprintln!(
                "      F88 six terms, no double count, local evap, ocean term at the FIRST hop: \
                 {arrives:.1} + {sp_ocean_local:.1} + **{evap_local:.2}** + {evap_surface:.2} + \
                 {retained:.1} + **{q_real:.1}** = {six_f88:.1} = **{:.1} %** of {budget:.1} · \
                 (removed from term 6: {q_double:.1} m³/s already counted at the coast)",
                100.0 * six_f88 / budget as f64
            );
            eprintln!(
                "      **F88 SIX TERMS, ocean and retained read at the END OF THE CHAIN**: \
                 {arrives:.1} + **{ocean_transitive:.1}** + {evap_local:.2} + {evap_surface:.2} + \
                 **{retained_transitive:.1}** + {q_real:.1} = **{six_f88c:.1}** = **{:.1} %** of \
                 {budget:.1}",
                100.0 * six_f88c / budget as f64
            );

            // negative control: the biggest EXOREIC detected lake must show up at the coast
            let ctrl = dr
                .lakes
                .iter()
                .filter(|l| l.base.id < 1_000_001 && l.lake_type == LakeType::Exorheic)
                .max_by(|a, b| a.area_km2.total_cmp(&b.area_km2));
            if let Some(c) = ctrl {
                let cells: Vec<usize> = (0..n).filter(|&k| dr.lake_map[k] == c.base.id).collect();
                let q_in = cells.iter().map(|&k| racc[k]).fold(0.0f32, f32::max);
                let deepest = cells
                    .iter()
                    .copied()
                    .min_by(|&x, &y| field.data[x].total_cmp(&field.data[y]))
                    .unwrap_or(0);
                let (term, steps, _) = terminus(deepest, &dr.flow, &field, w, h);
                let q_term = runoff_km2_to_m3s(racc[term]);
                eprintln!(
                    "      B control: biggest Exorheic lake {} ({:.1} km², local inflow {:.1} \
                     m³/s) routes {steps} steps to {} — accumulation there {:.1} m³/s ⇒ \
                     **{}**",
                    c.base.id,
                    c.area_km2,
                    runoff_km2_to_m3s(q_in),
                    if field.data[term] <= SEA { "the SEA" } else { "LAND" },
                    q_term,
                    if field.data[term] <= SEA && q_term >= runoff_km2_to_m3s(q_in) * 0.9 {
                        "PASSES"
                    } else {
                        "FAILS"
                    }
                );
            }

            // ══ F — why 4 850 coastal cells carry no segment ════════════════
            let ug: Vec<f32> = sorted(ungauged.iter().map(|u| u.0).collect());
            let above: Vec<&(f32, usize)> = ungauged.iter().filter(|u| u.0 >= head_km2).collect();
            let mut cause: HashMap<&str, (usize, f32)> = HashMap::new();
            for (a, k) in &above {
                let c = if *a < stream_km2 {
                    "below the EXPORT threshold stream_km2 (Finding 20's clip)"
                } else if dr.lake_map[*k] != 0 {
                    "on a lake footprint (clip_rivers_to_lakes)"
                } else {
                    "at or above stream_km2, off any lake — the SEAM"
                };
                let e = cause.entry(c).or_insert((0, 0.0));
                e.0 += 1;
                e.1 = e.1.max(*a);
            }
            eprintln!(
                "      F: ungauged arrival {ung_q:.1} m³/s on {} cells · area km² p50 {:.4} p90 \
                 {:.3} max {:.2} · **{} cells at or above head_threshold {head_km2}**",
                ug.len(),
                pct(&ug, 0.50),
                pct(&ug, 0.90),
                ug.last().copied().unwrap_or(0.0),
                above.len()
            );
            let mut cvv: Vec<(&&str, &(usize, f32))> = cause.iter().collect();
            cvv.sort_by(|a, b| b.1.0.cmp(&a.1.0));
            for (c, (cnt, mx)) in cvv {
                eprintln!(
                    "         {cnt:>5} cells ({:5.1} %) · max {mx:.2} km² · {c}",
                    100.0 * *cnt as f64 / above.len().max(1) as f64
                );
            }

            if gate.contains("control") {
                continue; // D and E are read on the shipped state only
            }

            // ══ E — basin 1000001's escape, on the RIGHT field ══════════════
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
                if cells.is_empty() {
                    continue;
                }
                match surface_lake_escape(b.id, &cells, &field, &bs.lake_map, w, h) {
                    Ok((sk, ek)) => {
                        let dp = dr.flow.direction[ek];
                        let dn = flow_naked.direction[ek];
                        eprintln!(
                            "      E: basin {} ({:.1} m³/s local) — saddle ({},{}) at {:.4} m, \
                             escape ({},{}) at {:.4} m (Δ {:.4} m = {:.2} u16 steps) · \
                             **DIR_NONE on production {} · on the bench field {}**",
                            b.id,
                            b.local_outflow_m3s,
                            sk % w,
                            sk / w,
                            c1_altitude_norm_to_metres(field.data[sk], &ss),
                            ek % w,
                            ek / w,
                            c1_altitude_norm_to_metres(field.data[ek], &ss),
                            c1_altitude_norm_to_metres(field.data[sk], &ss)
                                - c1_altitude_norm_to_metres(field.data[ek], &ss),
                            (field.data[sk] - field.data[ek]) * 65535.0,
                            dp == DIR_NONE,
                            dn == DIR_NONE
                        );
                        if dp != DIR_NONE {
                            let (term, steps, cyc) = terminus(ek, &dr.flow, &field, w, h);
                            eprintln!(
                                "         and on production's field the descent RUNS: {steps} \
                                 steps → {} (cycled {cyc})",
                                if field.data[term] <= SEA {
                                    "the SEA".to_string()
                                } else {
                                    format!(
                                        "({},{}) at {:.2} m",
                                        term % w,
                                        term / w,
                                        c1_altitude_norm_to_metres(field.data[term], &ss)
                                    )
                                }
                            );
                        }
                    }
                    Err(e) => eprintln!("      E: basin {} — no escape: {}", b.id, e.as_str()),
                }
            }

            // ══ D — the 623.6 m lake ════════════════════════════════════════
            let kf = DEEP_FLOOR.1 * w + DEEP_FLOOR.0;
            let deep_id = dr.lake_map[kf];
            let foot: Vec<usize> = (0..n).filter(|&k| dr.lake_map[k] == deep_id).collect();
            let m = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
            eprintln!(
                "      D1: the floor cell ({},{}) belongs to body {deep_id} ({} cells, {:.1} km²)",
                DEEP_FLOOR.0,
                DEEP_FLOOR.1,
                foot.len(),
                foot.len() as f32 * cell_km2
            );
            eprintln!(
                "      D1 the floor, pass by pass: **pass 0 {:.2} m → pass 1 {:.2} m → pass 2 \
                 (shipped) {:.2} m** ⇒ pass 1 carries **{:.1} %** of the cut",
                m(&pre, kf),
                m(&it1, kf),
                m(&field, kf),
                100.0 * (m(&pre, kf) - m(&it1, kf)) / (m(&pre, kf) - m(&field, kf)).abs().max(1e-6)
            );
            if !foot.is_empty() {
                let cut1: Vec<f32> = foot.iter().map(|&k| m(&pre, k) - m(&it1, k)).collect();
                let cut2: Vec<f32> = foot.iter().map(|&k| m(&pre, k) - m(&field, k)).collect();
                let pre_h: Vec<f32> = foot.iter().map(|&k| m(&pre, k)).collect();
                let del_h: Vec<f32> = foot.iter().map(|&k| m(&field, k)).collect();
                eprintln!(
                    "      D1 the footprint: median cut after pass 1 **{:.2} m**, after pass 2 \
                     **{:.2} m** ⇒ pass 1 carries **{:.1} %**",
                    median(&cut1),
                    median(&cut2),
                    100.0 * median(&cut1) as f64 / median(&cut2).abs().max(1e-6) as f64
                );
                eprintln!(
                    "      ⚠️ D1 population note: median-of-differences **{:.2} m** ≠ \
                     difference-of-medians **{:.2} m** (pre p50 {:.1} − delivered p50 {:.1}) — \
                     Finding 87-B published the SECOND (165.80 m) and called it the median cut. \
                     The paired statistic is the first.",
                    median(&cut2),
                    median(&pre_h) - median(&del_h),
                    median(&pre_h),
                    median(&del_h)
                );
                // D2 — what it relaxes TOWARD, and whether the Finding 83 bound could fire
                let mut bound_would_fire = 0usize;
                let mut landed_on_receiver = 0usize;
                let mut area_of: Vec<(f32, f32)> = Vec::new(); // (log10 area km², cut m)
                for &k in &foot {
                    let d = dr.flow.direction[k];
                    if d == DIR_NONE {
                        continue;
                    }
                    let nx = ((k % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
                    let ny = ((k / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
                    let r = ny * w + nx;
                    // the Finding 83 bound floors the TARGET at sea + 0.5 m
                    if m(&pre, r) < 0.5 {
                        bound_would_fire += 1;
                    }
                    if (m(&field, k) - m(&pre, r)).abs() < 1.0 {
                        landed_on_receiver += 1;
                    }
                    let a = dr.flow.accumulation.data[k] * cell_km2;
                    area_of.push((a.max(1e-6).log10(), m(&pre, k) - m(&field, k)));
                }
                let big: Vec<f32> = area_of.iter().filter(|p| p.0 >= 0.0).map(|p| p.1).collect();
                let small: Vec<f32> = area_of.iter().filter(|p| p.0 < -2.0).map(|p| p.1).collect();
                eprintln!(
                    "      D2: of {} footprint cells — the F83 bound could fire on **{}** \
                     ({:.2} %); **{}** ({:.1} %) landed within 1 m of their receiver's \
                     pre-incision height · cut by drainage area: A ≥ 1 km² median **{:.1} m** \
                     (n={}), A < 0.01 km² median **{:.1} m** (n={})",
                    foot.len(),
                    bound_would_fire,
                    100.0 * bound_would_fire as f64 / foot.len() as f64,
                    landed_on_receiver,
                    100.0 * landed_on_receiver as f64 / foot.len() as f64,
                    if big.is_empty() { f32::NAN } else { median(&big) },
                    big.len(),
                    if small.is_empty() { f32::NAN } else { median(&small) },
                    small.len()
                );
                // D3 — the C-3 and C-2 fields at this place
                if let Some(kf_map) = k_field.as_ref() {
                    // ⚠️ population: the K field covers the whole grid. Comparing a land
                    // footprint against a p50 that includes 55.8 M sea cells is Finding 63's
                    // bias. LAND ONLY.
                    let all = sorted(
                        (0..n).filter(|&k| land[k]).map(|k| kf_map[k]).collect::<Vec<f32>>(),
                    );
                    let here: Vec<f32> = foot.iter().map(|&k| kf_map[k]).collect();
                    eprintln!(
                        "      D3: C-3 erodibility multiplier — map p50 {:.3} · this footprint \
                         p50 **{:.3}** (p10 {:.3}, p90 {:.3}) ⇒ ratio **{:.2}×**",
                        pct(&all, 0.50),
                        median(&here),
                        pct(&sorted(here.clone()), 0.10),
                        pct(&sorted(here.clone()), 0.90),
                        median(&here) / pct(&all, 0.50).max(1e-6)
                    );
                }
                let (fu, fv) = (DEEP_FLOOR.0 as f32 / w as f32, DEEP_FLOOR.1 as f32 / h as f32);
                let mut nearest = f32::INFINITY;
                for e in &edifices {
                    let du = (e.center_uv.0 - fu).abs().min(1.0 - (e.center_uv.0 - fu).abs());
                    let dv = (e.center_uv.1 - fv).abs().min(1.0 - (e.center_uv.1 - fv).abs());
                    nearest = nearest.min((du * du + dv * dv).sqrt() * DOMAIN_KM);
                }
                eprintln!(
                    "      D3: C-2 volcanism — {} edifices placed, nearest to the floor cell \
                     **{:.1} km**",
                    edifices.len(),
                    nearest
                );
                // D4 — how many lakes are in this class
                let mut klass = 0usize;
                let mut scanned = 0usize;
                let mut rows4: Vec<(u32, f32, f32, f32)> = Vec::new();
                for l in dr.lakes.iter().filter(|l| l.area_km2 >= 1.0) {
                    let cells: Vec<usize> =
                        (0..n).filter(|&k| dr.lake_map[k] == l.base.id).collect();
                    if cells.is_empty() {
                        continue;
                    }
                    scanned += 1;
                    let cuts: Vec<f32> = cells.iter().map(|&k| m(&pre, k) - m(&field, k)).collect();
                    let cm = median(&cuts);
                    // ⚠️ the rim is the ring JUST OUTSIDE the footprint — the WALL — which is
                    // how Finding 87-B measured 42.8°. Measuring the inner edge instead (my first
                    // instrument) returned 0 bodies in this class and was the wrong population.
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
                    let rim: Vec<f32> =
                        ring.iter().map(|&k| slope_deg(&field, k, n2m, w, h)).collect();
                    let rm = if rim.is_empty() { 0.0 } else { median(&rim) };
                    if cm > 50.0 && rm > 30.0 {
                        klass += 1;
                        rows4.push((l.base.id, l.area_km2, cm, rm));
                    }
                }
                rows4.sort_by(|a, b| b.2.total_cmp(&a.2));
                // the promoted object's own numbers, on the SAME instrument as the class test
                {
                    let cuts: Vec<f32> = foot.iter().map(|&k| m(&pre, k) - m(&field, k)).collect();
                    let inside: HashSet<usize> = foot.iter().copied().collect();
                    let mut ring: HashSet<usize> = HashSet::new();
                    let mut ring_weighted: Vec<f32> = Vec::new();
                    for &k in &foot {
                        for d in 0..8 {
                            let nx = ((k % w) as i32 + D8_DX[d]).rem_euclid(w as i32) as usize;
                            let ny = ((k / w) as i32 + D8_DY[d]).rem_euclid(h as i32) as usize;
                            let nk = ny * w + nx;
                            if !inside.contains(&nk) {
                                ring.insert(nk);
                                ring_weighted.push(slope_deg(&field, nk, n2m, w, h));
                            }
                        }
                    }
                    let dedup: Vec<f32> =
                        ring.iter().map(|&k| slope_deg(&field, k, n2m, w, h)).collect();
                    eprintln!(
                        "      D4 the promoted body {deep_id}: median cut **{:.1} m** · rim p50 \
                         **{:.1}°** (deduped ring, {} cells) vs **{:.1}°** (Finding 87-B's \
                         multiset weighting, {} entries) · p90 {:.1}° / {:.1}°",
                        median(&cuts),
                        median(&dedup),
                        dedup.len(),
                        median(&ring_weighted),
                        ring_weighted.len(),
                        pct(&sorted(dedup.clone()), 0.90),
                        pct(&sorted(ring_weighted.clone()), 0.90)
                    );
                }
                eprintln!(
                    "      D4: of {scanned} delivered bodies ≥ 1 km², **{klass}** have median cut \
                     > 50 m AND rim p50 (the OUTSIDE ring) > 30° — {}",
                    if klass >= 5 { "**a POPULATION**" } else { "a handful" }
                );
                for (id, a, cm, rm) in rows4.iter().take(8) {
                    eprintln!("         body {id} — {a:7.1} km² · cut {cm:6.1} m · rim {rm:.1}°");
                }
            }
        }
    }
    eprintln!("\n==========  end Finding 88  ==========\n");
}
