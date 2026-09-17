//! The FINAL HD drainage assembly — the climate-dependent tail of the production pipeline.
//!
//! ## Why this lives in `ymir-core` and not in the viz bridge
//!
//! Every step below is core logic; only the caching and the progress events around it were
//! ever viz's business. It nevertheless lived in `ymir-viz::bridge::c1::hd` as a private
//! `build_hd_drainage`, which made it **unreachable from a test**. The consequence is recorded
//! in ADR Finding 56: two separate benches (the first H-1 bench, then the channel-head
//! invariant suite) re-implemented the chain by hand, both dropped the pre-breach lake carry,
//! and both therefore reported `lakes = 0` in every configuration — a metric block that
//! measured nothing, twice.
//!
//! So the assembly is exported here and production calls it. A bench that wants production's
//! lake population now *calls production's code*; it cannot carry the wrong population,
//! because there is only one carry. (ADR method rule 10's sibling: a rule that is not
//! executable does not protect.)
//!
//! ## The order is load-bearing
//!
//! 1. final drainage on the **breached** field, **with** the climate discharge (Finding 22);
//! 2. adopt the **pre-breach** lake geometry — the breach destroys the depressions, so the
//!    post-breach drainage has no lakes to find (this is the step both benches dropped);
//! 3. H-1c: run the surface water balance on that geometry (classification + endorheic
//!    equilibrium level/footprint), mutating `lake_map`;
//! 4. the below-sea basin merge + submerged-lake cleanup;
//! 5. **then** clip the rivers — after the footprint is final, never before, or mouths orphan;
//! 6. append the below-sea spillways as `SegmentKind::Spillway` rows;
//! 7. the geographic-scale post-process.

use crate::grid::GridF32;
use crate::tectonics_c1::cached_product::HdDrainageBundle;
use crate::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use crate::tectonics_c1::drainage::{
    C1DrainageConfig, C1DrainageResult, DrainageClimate, LakeType, SegmentKind, SegmentRow,
    apply_geo_scale_ratio, apply_lake_water_balance, below_sea_basin_lakes_infil,
    c1_drainage_windowed_infil, clip_rivers_to_lakes, resolve_exorheic_without_outlet,
    uncovered_below_sea_components,
};
use crate::terrain::flow::RiverSegment;

/// Build the FINAL HD drainage + wetland mask from the eroded field and the climate.
///
/// `prebreach` is the climate-free drainage of the **un-breached** field: pass `Some` on the
/// relief-v3 path (`eroded` is then the breached field) so the geometrically correct lakes are
/// carried forward; `None` reproduces the legacy path (plain drainage on the raw eroded, no
/// carried lakes, no water balance). `verbose` gates the per-stage `eprintln` accounting.
///
/// A PURE function of exactly its inputs — which is what makes caching it under
/// [`crate::tectonics_c1::cached_product::hd_drainage_key`] stale-proof.
#[allow(clippy::too_many_arguments)]
pub fn assemble_hd_drainage(
    eroded: &GridF32,
    climate: &DrainageClimate<'_>,
    prebreach: Option<C1DrainageResult>,
    dcfg: &C1DrainageConfig,
    ss: &SteinSteinParams,
    window_km: f32,
    geo_scale_ratio: f32,
    // H-1: per-cell infiltrated fraction (None → pre-H-1 balance, byte-identical).
    infiltration: Option<&[f32]>,
    verbose: bool,
) -> HdDrainageBundle {
    // relief-v3: final drainage on the breached field, carrying the pre-breach lakes; legacy path
    // (prebreach None): plain drainage on the raw eroded. Both with the climate discharge (Finding 22).
    let (mut drainage, carried_geometric_lakes) = if let Some(prebreach) = prebreach {
        let mut dr =
            c1_drainage_windowed_infil(eroded, Some(climate), dcfg, ss, window_km, infiltration);
        // The pre-breach lakes are the GEOMETRICALLY correct ones (the breach destroys the
        // depressions) — but they were detected with `climate = None`, so they carry NO
        // water balance. Their geometry is adopted; their classification is fixed below.
        dr.lakes = prebreach.lakes;
        dr.lake_map = prebreach.lake_map;
        (dr, true)
    } else {
        (
            c1_drainage_windowed_infil(eroded, Some(climate), dcfg, ss, window_km, infiltration),
            false,
        )
    };
    // H-1 — THE MISSING LINK, and a CORRECTION, so it runs BY DEFAULT (not gated): lakes
    // classified without ever seeing the climate are simply wrong. The relief-v3 path
    // carries the climate-free pre-breach lakes, so they were exorheic by pure geometry
    // ("the outlet reaches the sea"). Run the balance on them and adopt the CLASSIFICATION
    // ONLY — geometry, levels and footprints stay untouched (adopting an endorheic
    // equilibrium LEVEL would be H-2 by the back door). Crater types (C-2) are never
    // overwritten. The optional H-1 INFILTRATION rides in `infiltration` (measured a
    // secondary term: 0–3 lakes, against 53–100 % from this reclassification alone).
    if carried_geometric_lakes {
        let (w, h) = (eroded.width, eroded.height);
        let cell_km2 = (window_km / w as f32).powi(2);
        let before_n = drainage.lakes.len();
        let before_km2: f32 = drainage.lakes.iter().map(|l| l.area_km2).sum();
        // H-1c — APPLY the balance: classify AND settle endorheic basins at their
        // evaporative equilibrium (level + footprint), draining the exposed floor from
        // `lake_map`. Runs BEFORE `below_sea_basin_lakes_infil` and BEFORE
        // `clip_rivers_to_lakes` so both see the FINAL footprint — river tracks are clipped
        // to the retreated outline instead of ending in the void (the orphaned-mouth defect
        // already fixed twice: enumerate inlets AFTER the footprint is known).
        let lakes_in = std::mem::take(&mut drainage.lakes);
        drainage.lakes = apply_lake_water_balance(
            eroded,
            &drainage.flow,
            climate,
            cell_km2,
            ss,
            &lakes_in,
            &mut drainage.lake_map,
            infiltration,
            w,
            h,
        );
        if verbose {
            let endo = drainage.lakes.iter().filter(|l| l.lake_type == LakeType::Endorheic).count();
            let after_km2: f32 = drainage.lakes.iter().map(|l| l.area_km2).sum();
            eprintln!(
                "[HD] H-1c surface water balance APPLIED: {} → {} lakes | {:.0} → {:.0} km² water ({:.0} km² floor exposed) | {} endorheic",
                before_n,
                drainage.lakes.len(),
                before_km2,
                after_km2,
                (before_km2 - after_km2).max(0.0),
                endo
            );
        }
    }
    let (wetland_mask, below_sea_spillways) = {
        let bs = below_sea_basin_lakes_infil(
            eroded,
            climate,
            dcfg,
            ss,
            window_km,
            Some(&drainage.lake_map),
            infiltration,
        );
        let (mut endo, mut exo) = (0usize, 0usize);
        for lk in &bs.lakes {
            match lk.lake_type {
                LakeType::Endorheic => endo += 1,
                LakeType::Exorheic => exo += 1,
                // Crater types are assigned later (post-drainage), never here.
                // `Unresolved` likewise: ADR Finding 86 posts it at the END of the chain, after
                // the segments are final, so a below-sea lake cannot carry it at this point.
                LakeType::CraterAcidic | LakeType::CraterNeutral | LakeType::Unresolved => {}
            }
        }
        use std::collections::{HashMap, HashSet};
        let mut det_before: HashMap<u32, usize> = HashMap::new();
        for &id in &drainage.lake_map {
            if id != 0 && id < 1_000_001 {
                *det_before.entry(id).or_default() += 1;
            }
        }
        let mut det_after: HashMap<u32, usize> = HashMap::new();
        for k in 0..bs.lake_map.len() {
            if bs.lake_map[k] != 0 {
                drainage.lake_map[k] = bs.lake_map[k];
            } else if drainage.lake_map[k] != 0 && drainage.lake_map[k] < 1_000_001 {
                *det_after.entry(drainage.lake_map[k]).or_default() += 1;
            }
        }
        let absorbed: HashSet<u32> = det_before
            .keys()
            .copied()
            .filter(|id| det_after.get(id).copied().unwrap_or(0) == 0)
            .collect();
        let n_absorbed = absorbed.len();
        drainage.lakes.retain(|lk| lk.base.id >= 1_000_001 || !absorbed.contains(&lk.base.id));
        drainage.lakes.extend(bs.lakes);
        if verbose {
            let (to_sea, chained) = (
                bs.spillways.iter().filter(|s| s.chained_into.is_none()).count(),
                bs.spillways.iter().filter(|s| s.chained_into.is_some()).count(),
            );
            if n_absorbed > 0 {
                eprintln!(
                    "[HD] below-sea cleanup: {n_absorbed} detected lake(s) submerged by a filled below-sea lake -> dropped"
                );
            }
            eprintln!(
                "[HD] below-sea basins: {} lakes ({exo} exorheic, {endo} endorheic); {} spillways ({to_sea} -> sea, {chained} chained)",
                exo + endo,
                bs.spillways.len()
            );
        }
        (bs.wetland, bs.spillways)
    };
    let before_seg = drainage.rivers.segments.len();
    clip_rivers_to_lakes(&mut drainage);
    if verbose {
        eprintln!(
            "[HD] rivers clipped to lakes: {before_seg} -> {} segments (terminate at sinks)",
            drainage.rivers.segments.len()
        );
    }
    // A below-sea basin's outflow is a REAL flow and must stay on the map, but it is not a
    // hierarchised watercourse — giving it a Strahler order produced "order 1 draining
    // 88 468 km²" at the head of the discharge sort. It is tagged `Spillway`: the KIND is
    // authoritative, `strahler_order` must not be read for it, and `width_m` (from
    // discharge) is what to render. `segment_source_lake` is `None` when the source basin
    // sits BELOW the lake-inventory floor, so the consumer is never handed an id absent
    // from `lakes.json`.
    let inventoried: std::collections::HashSet<u32> =
        drainage.lakes.iter().map(|l| l.base.id).collect();
    for sw in &below_sea_spillways {
        // ONE row, every array at once (`push_segment`). Pushing the arrays individually here
        // dropped a field twice — see `SegmentRow`.
        let (lx, ly) = *sw.points.last().unwrap_or(&(0, 0));
        drainage.push_segment(SegmentRow {
            segment: RiverSegment {
                points: sw.points.clone(),
                strahler_order: 1, // MEANINGLESS on a spillway — see `segment_kind`
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
            // A spillway's path is traced over a col, OUTSIDE the accumulation network, so the
            // raster reads ~0 along it (measured). Its contributing area is the basin's, which
            // `drainage_km2` already carries — convert it back to cells rather than reading a
            // raster that does not describe this path.
            catchment_cells: {
                let k = ly as usize * eroded.width + lx as usize;
                drainage.flow.accumulation.data.get(k).copied().unwrap_or(0.0)
            },
            // Per-point discharge (Finding 46). A spillway's discharge is uniform along its
            // traced path — one outflow over a col, not a hierarchy accumulating tributaries.
            discharge_profile_m3s: vec![sw.discharge_m3s; sw.points.len()],
            kind: SegmentKind::Spillway,
            source_lake: inventoried.contains(&sw.lake_id).then_some(sw.lake_id),
        });
    }
    debug_assert!(
        drainage.segment_arrays_aligned(),
        "spillway append desynchronised a parallel array"
    );
    apply_geo_scale_ratio(&mut drainage, geo_scale_ratio, &dcfg.thresholds);
    if verbose && geo_scale_ratio != 1.0 {
        eprintln!(
            "[HD] geographic scale ratio {geo_scale_ratio:.2} -> hydrology signifies x{:.1} area",
            geo_scale_ratio * geo_scale_ratio
        );
    }
    // ADR 0001 Finding 86 — the LAST thing before the bundle leaves: an exorheic lake has an
    // outlet reach, or it is not exorheic. It has to be here and not in the balance, because the
    // question is whether a REACH was emitted, which is only knowable once the segments are final
    // (`clip_rivers_to_lakes` and the spillway append have both run). Finding 29 named the defect
    // as H2, Finding 30 forbade relabelling it Endorheic, Finding 39 closed it on the below-sea
    // path only, and `apply_lake_water_balance` (H-1c, drainage.rs:1078) still posts `Exorheic`
    // from `a_eq >= a_sill` with no trace.
    let unresolved = resolve_exorheic_without_outlet(&mut drainage);
    if verbose && !unresolved.is_empty() {
        eprintln!(
            "[HD] ADR Finding 86: {} lake(s) claimed Exorheic with no outlet reach and are now \
             Unresolved: {unresolved:?}",
            unresolved.len()
        );
    }
    // ADR 0001 Finding 92-B — **THE FINDING 38 INVARIANT, PINNED AT LAST.** Every enclosed
    // below-sea component must be covered by a water body. Finding 38 closed this in 2024 (68
    // river mouths on `water_class == 2`, `lake_map == 0` slivers — "neither lake nor sea") and
    // left no permanent guard; Finding 86's relevel re-opened it for two components carrying
    // 34.8 m³/s, and the thing that found it was a viz screenshot. A production assertion, not a
    // `#[cfg(test)]` one (Finding 72: the ledger that is not read is not a ledger), so it holds on
    // every path that builds a bundle and not only on the ones a bench remembers to call.
    let uncovered = uncovered_below_sea_components(
        eroded,
        &drainage.lake_map,
        crate::tectonics_c1::drainage::C1_SEA_LEVEL_NORM,
        eroded.width,
        eroded.height,
    );
    assert!(
        uncovered.is_empty(),
        "ADR Finding 38/92-B: {} enclosed below-sea component(s) carry no water body. \
         (floor cell, cells): {:?}. This is the invariant Finding 38 closed and Finding 86's \
         `MergedUnionRelevel` re-opened; see `separate_unclaimed_regions`.",
        uncovered.len(),
        uncovered.iter().take(8).collect::<Vec<_>>()
    );
    HdDrainageBundle { drainage, wetland: wetland_mask }
}
