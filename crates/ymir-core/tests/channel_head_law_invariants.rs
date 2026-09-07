//! ADR Finding 56 — the invariant suite for the slope-dependent channel head.
//!
//! `A_c` governs the hillslope/channel split EVERYWHERE, so the law changes the network far
//! beyond the coast. A coastal gain paid for by a degraded network would be a bad trade and the
//! coastline metrics cannot see it. This runs the WHOLE production hydro chain — conditioning,
//! climate, drainage, lakes — with the law OFF and ON, at both resolutions, and reports what
//! moved.
//!
//! ## Two corrections this file carries (Finding 56b)
//!
//! **1. The lake population is production's.** The first version of this bench measured
//! `c1_drainage_windowed` on the BREACHED field, whose lakes production throws away — so it
//! read `lakes 0` in every configuration and the invariants "passed" over an empty set. It is
//! now assembled by [`assemble_hd_drainage`], the same function production calls, so the
//! pre-breach lake carry cannot be got wrong here: there is only one carry.
//!
//! **2. W/D per order is a VALLEY cross-section, not width ÷ relief.** The old version divided
//! a hydraulic channel width (metres, from discharge) by a 1 km-window relief (hundreds of
//! metres) — two different objects, and every order printed `0.0`. The ratio now comes from a
//! transect perpendicular to the flow, walked out to each shoulder; the channel width is
//! reported separately, in metres, so the sub-metre discharge regime stays VISIBLE instead of
//! rounding into a ratio.
//!
//! Every metric block is registered in a [`Sweep`] (method rule 10) and verified at the end:
//! any reported quantity that did not respond, and any invariant checked over an empty
//! population, fails the test loudly.
//!
//! Climate: the arid-hot test bed (25°, span 10) — the standing config for below-sea basins.
//! Note the TERRAIN is climate-independent (`upscale_from_c1_with_progress` takes no climate),
//! so the coastline results elsewhere are the same at any latitude; climate matters HERE,
//! because the lakes and the discharge depend on it.
//!
//! Run: cargo test -p ymir-core --release --test channel_head_law_invariants -- --ignored --nocapture

use std::collections::{HashMap, HashSet};

use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::erosion::stream_power::{CHANNEL_HEAD_S_MIN, CHANNEL_HEAD_S_REF, ChannelHeadLaw};
use ymir_core::grid::GridF32;
use ymir_core::metric_sweep::Sweep;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, C1DrainageResult, DrainageClimate, LakeType, Navigability, SegmentKind,
    c1_drainage_windowed, exorheic_lakes_missing_outlet,
};
use ymir_core::tectonics_c1::hd_assembly::assemble_hd_drainage;
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::{
    c1_altitude_norm_to_metres, upscale_from_c1_with_progress,
};
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::flow::breach_monotone;
use ymir_core::terrain::upscale::{ProductionHdOpts, production_hd_config};

const PSEED: u64 = 10_481_999_410_520_546_993;
const DOMAIN_KM: f32 = 400.0;
const LAT: f32 = 25.0;
const SPAN: f32 = 10.0;
const SEA: f32 = 0.5;
/// Production ships the geographic-scale dial at 1.0 (`Workspace::default`), so that is what
/// the invariants must be measured at. It scales the SIGNIFIED hydrology, never the terrain.
const GEO_SCALE_RATIO: f32 = 1.0;
/// Altitude tolerance, metres. Below the vertical quantisation of the normalised field.
const TOL_M: f32 = 0.5;
/// Valley transect: how far out to look for a shoulder, and the drop that says we crossed it.
const TRANSECT_MAX_KM: f32 = 2.0;
const SHOULDER_DROP_M: f32 = 10.0;

fn terrain(target: usize, law: bool) -> GridF32 {
    let ss = SteinSteinParams::default();
    let run_cfg = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let mut state = init_c1_state_phase_2_r7(64, PSEED, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, &run_cfg, &C1Closures::default(), |_, _| {});
    let seed = WorldSeed::new(PSEED);
    let volc = VolcanismConfig { enabled: true, domain_km: DOMAIN_KM, ..Default::default() };
    let edifices = place_edifices(&state, &kin, &seed, DOMAIN_KM, &volc);
    let mut cfg = production_hd_config(&ProductionHdOpts {
        target_size: target,
        domain_km: DOMAIN_KM,
        depth_scale_m: ss.depth_scale_m as f32,
        sample_origin: [0.0, 0.578_125],
        sample_size: 1.0,
        amplitude_base: 0.04,
        mfd_p: 2.0,
        lithology: LithologyConfig {
            enabled: true,
            soft_multiplier: 10.0,
            volcanic_multiplier: 3.0,
            rift_age_threshold: 1.0,
        },
        fracture: FractureConfig {
            enabled: true,
            amplitude: 6.0,
            decay_km: 25.0,
            domain_km: DOMAIN_KM,
            ..Default::default()
        },
    });
    if law {
        if let Some(sp) = cfg.stream_power.as_mut() {
            sp.a_c_slope_law =
                Some(ChannelHeadLaw { s_ref: CHANNEL_HEAD_S_REF, s_min: CHANNEL_HEAD_S_MIN });
        }
    }
    upscale_from_c1_with_progress(
        &state,
        &run_cfg.iso_config,
        &ss,
        &seed,
        &cfg,
        &edifices,
        &volc,
        Some(&kin),
        &mut |_| {},
        &|| false,
    )
    .0
    .heightmap
}

/// What the hydro chain leaves behind, including the two things the PRE-breach field is the
/// only place to measure.
struct Hydro {
    /// The BREACHED field — what production exports and what the drainage is computed on.
    field: GridF32,
    /// Production's final drainage: pre-breach lakes carried, H-1c applied, below-sea merged,
    /// rivers clipped, spillways appended, geographic scale applied.
    dr: C1DrainageResult,
    /// Closed depressions on the RAW field, i.e. C-1's conditioning TARGET. Measuring these
    /// after the breach is measuring nothing — the breach removes them by construction, which
    /// is the other half of the same mistake as the post-breach lakes.
    pre_pits: usize,
    /// Lakes detected pre-breach, BEFORE the water balance reclassified or settled them.
    pre_lakes: usize,
}

/// The production hydro chain, end to end: pre-breach drainage → breach conditioning → climate
/// → [`assemble_hd_drainage`] (which is the code production runs, not a copy of it).
fn hydro(raw: &GridF32, ss: &SteinSteinParams) -> Hydro {
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let (w, h) = (raw.width, raw.height);

    // The pre-breach drainage: climate-free, and the ONLY place the closed depressions and the
    // geometrically correct lakes exist.
    let pre = c1_drainage_windowed(raw, None, &dcfg, ss, DOMAIN_KM);
    let pre_pits =
        pre.flow.filled.data.iter().zip(raw.data.iter()).filter(|(f, r)| **f > **r + 1e-6).count();
    let pre_lakes = pre.lakes.len();

    let field = breach_monotone(raw, &pre.flow.filled, &pre.lake_map, SEA, w, h);
    let climate = c1_climate_placed(&field, ss, LAT, SPAN, &PrecipParams::default(), DOMAIN_KM);
    let dclim = DrainageClimate {
        precip_internal: &climate.precipitation,
        temperature: &climate.temperature,
    };
    let bundle = assemble_hd_drainage(
        &field,
        &dclim,
        Some(pre),
        &dcfg,
        ss,
        DOMAIN_KM,
        GEO_SCALE_RATIO,
        None,
        false,
    );
    Hydro { field, dr: bundle.drainage, pre_pits, pre_lakes }
}

fn pct(v: &mut Vec<f32>, f: f32) -> f32 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    if v.is_empty() {
        return 0.0;
    }
    v[(((v.len() - 1) as f64 * f as f64).round() as usize).min(v.len() - 1)]
}

/// The VALLEY cross-section at `(x, y)`, measured on a transect PERPENDICULAR to `(dx, dy)`.
///
/// Walks out one cell at a time on each side, tracking the running maximum altitude, and calls
/// the shoulder where the ground has dropped [`SHOULDER_DROP_M`] below that maximum (we are
/// descending into the next valley) or [`TRANSECT_MAX_KM`] is reached. Returns
/// `(width_m, depth_m)` with the width the sum of both half-distances-to-shoulder and the depth
/// the SHALLOWER of the two rises — so one high flank cannot inflate the depth of an
/// asymmetric valley.
///
/// `None` when either side runs off the grid or when a shoulder is under a metre of rise (a
/// transect that never left the flat has no valley to measure, and dividing by it is how the
/// old version produced ratios in the thousands).
fn valley_cross_section(
    field: &GridF32,
    to_m: &dyn Fn(f32) -> f32,
    x: usize,
    y: usize,
    dx: f32,
    dy: f32,
    m_per_cell: f32,
) -> Option<(f32, f32)> {
    let (w, h) = (field.width, field.height);
    let len = (dx * dx + dy * dy).sqrt();
    if len < 1e-6 {
        return None;
    }
    // Perpendicular to the flow direction, unit length.
    let (px, py) = (-dy / len, dx / len);
    let max_steps = (TRANSECT_MAX_KM * 1000.0 / m_per_cell).round().max(2.0) as i32;
    let bed = to_m(field.data[y * w + x]);
    let mut half = [0.0f32; 2];
    let mut rise = [0.0f32; 2];
    for (side, sign) in [(0usize, 1.0f32), (1usize, -1.0f32)] {
        let (mut best, mut best_at) = (bed, 0.0f32);
        for s in 1..=max_steps {
            let fx = x as f32 + 0.5 + sign * px * s as f32;
            let fy = y as f32 + 0.5 + sign * py * s as f32;
            let (ix, iy) = (fx.floor() as i32, fy.floor() as i32);
            if ix < 0 || iy < 0 || ix as usize >= w || iy as usize >= h {
                return None; // the transect left the grid — this segment cannot be measured
            }
            let a = to_m(field.data[iy as usize * w + ix as usize]);
            if a > best {
                best = a;
                best_at = s as f32;
            } else if best - a > SHOULDER_DROP_M {
                break; // over the divide, descending the far side
            }
        }
        half[side] = best_at.max(1.0) * m_per_cell;
        rise[side] = best - bed;
    }
    let depth = rise[0].min(rise[1]);
    if depth < 1.0 {
        return None;
    }
    Some((half[0] + half[1], depth))
}

#[allow(clippy::too_many_arguments)]
fn report(
    label: &str,
    target: usize,
    raw: &GridF32,
    hy: &Hydro,
    ss: &SteinSteinParams,
    sw: &mut Sweep,
) {
    let (field, dr) = (&hy.field, &hy.dr);
    let (w, h) = (field.width, field.height);
    let n = w * h;
    let km = DOMAIN_KM / target as f32;
    let m_per_cell = km * 1000.0;
    let cell_km2 = km * km;
    let to_m = |v: f32| c1_altitude_norm_to_metres(v, ss);
    sw.config(label);
    eprintln!("\n──────── {label} ────────");

    // ── LAKE INVARIANTS, over the WHOLE production population ───────────────────────
    // One pass to index the footprints: the per-lake full-grid scan the old version did is
    // O(lakes × 67 M) at 8192².
    let mut foot: HashMap<u32, Vec<usize>> = HashMap::new();
    for (k, &id) in dr.lake_map.iter().enumerate() {
        if id != 0 {
            foot.entry(id).or_default().push(k);
        }
    }
    let live: HashSet<u32> = dr.lakes.iter().map(|l| l.base.id).collect();
    let (exo, endo) = (
        dr.lakes.iter().filter(|l| l.lake_type == LakeType::Exorheic).count(),
        dr.lakes.iter().filter(|l| l.lake_type == LakeType::Endorheic).count(),
    );
    let crater = dr.lakes.len() - exo - endo;
    let water_km2 = foot.values().map(|v| v.len()).sum::<usize>() as f32 * cell_km2;

    // 1. footprint at or below level, and CONNECTED to the lowest point
    // 2. depth == level − floor, on the exported field AND on the detection field
    // 3. no overlapping / duplicated / dangling footprints
    let mut v_above = 0usize;
    let mut v_split = 0usize; // footprint not 4-connected from its lowest cell
    let mut v_depth_exported = 0usize;
    let mut v_depth_raw = 0usize;
    let mut v_area = 0usize; // declared area_km2 disagrees with the footprint
    let mut v_empty = 0usize; // a live lake with no cells in lake_map
    let mut split_endo = 0usize; // of `v_split`, how many are endorheic
    let mut seen_ids: HashSet<u32> = HashSet::new();
    let mut v_dup = 0usize;
    for lk in &dr.lakes {
        if !seen_ids.insert(lk.base.id) {
            v_dup += 1; // two lakes claiming one id — the only way footprints can overlap
        }
        let cells = match foot.get(&lk.base.id) {
            Some(c) if !c.is_empty() => c,
            _ => {
                v_empty += 1;
                continue;
            }
        };
        let lvl = lk.level_m;
        let mut floor = f32::INFINITY;
        let mut floor_raw = f32::INFINITY;
        let mut argmin = cells[0];
        for &k in cells {
            let a = to_m(field.data[k]);
            if a < floor {
                floor = a;
                argmin = k;
            }
            floor_raw = floor_raw.min(to_m(raw.data[k]));
            if a > lvl + TOL_M {
                v_above += 1;
            }
        }
        if ((lvl - floor) - lk.depth_m).abs() > 1.0 {
            v_depth_exported += 1;
        }
        if ((lvl - floor_raw) - lk.depth_m).abs() > 1.0 {
            v_depth_raw += 1;
        }
        if (lk.area_km2 - cells.len() as f32 * cell_km2).abs() > 1.5 * cell_km2 {
            v_area += 1;
        }
        // 4-connected flood fill from the lowest cell, restricted to the footprint
        let member: HashSet<usize> = cells.iter().copied().collect();
        let mut reached: HashSet<usize> = HashSet::new();
        let mut stack = vec![argmin];
        reached.insert(argmin);
        while let Some(k) = stack.pop() {
            let (x, y) = (k % w, k / w);
            for (nx, ny) in [(x.wrapping_sub(1), y), (x + 1, y), (x, y.wrapping_sub(1)), (x, y + 1)]
            {
                if nx >= w || ny >= h {
                    continue;
                }
                let nk = ny * w + nx;
                if member.contains(&nk) && reached.insert(nk) {
                    stack.push(nk);
                }
            }
        }
        if reached.len() < cells.len() {
            v_split += 1;
            if lk.lake_type == LakeType::Endorheic {
                split_endo += 1;
            }
        }
    }
    // Split by id range. Below-sea basin ids start at 1_000_001; the DETECTED lakes use small
    // ids. `below_sea_basin_lakes_infil` marks EVERY below-sea sink in `lake_map` for sink
    // validity but only inventories those clearing `INVENTORY_MIN_CELLS = 4` — its own unit
    // test asserts exactly that ("marked, not listed"). So a below-sea dangling id is the known
    // residual of Finding 33 Part A; a DETECTED dangling id would be a new defect.
    let v_dangling_belowsea =
        foot.keys().filter(|id| !live.contains(id) && **id >= 1_000_001).count();
    let v_dangling_detected =
        foot.keys().filter(|id| !live.contains(id) && **id < 1_000_001).count();
    let v_dangling = v_dangling_belowsea + v_dangling_detected;

    // 4. level <= inlet arrival altitude — a river must not arrive BELOW the lake surface
    // 5. no mouth without membership — a terminus touching an id absent from `lakes`
    let by_id: HashMap<u32, &ymir_core::tectonics_c1::drainage::C1Lake> =
        dr.lakes.iter().map(|l| (l.base.id, l)).collect();
    let mut v_inlet = 0usize;
    let mut n_inlets = 0usize;
    let mut orphan_mouths = 0usize;
    let mut n_termini = 0usize;
    for s in &dr.rivers.segments {
        let &(ex, ey) = s.points.last().unwrap_or(&(0, 0));
        n_termini += 1;
        let mut touched: HashSet<u32> = HashSet::new();
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let (nx, ny) = (ex as i32 + dx, ey as i32 + dy);
                if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                    continue;
                }
                let id = dr.lake_map[ny as usize * w + nx as usize];
                if id != 0 {
                    touched.insert(id);
                }
            }
        }
        for id in touched {
            if !live.contains(&id) {
                orphan_mouths += 1;
                continue;
            }
            n_inlets += 1;
            let arrival = to_m(field.data[ey as usize * w + ex as usize]);
            if by_id[&id].level_m > arrival + TOL_M {
                v_inlet += 1;
            }
        }
    }

    // 6. exorheic implies a traced outlet reaching a sink
    let missing_outlet = exorheic_lakes_missing_outlet(dr).len();

    // 7. monotone long profiles, FLAT LAKE CROSSINGS TOLERATED — a rise between two points
    //    that both sit inside a lake footprint is the lake surface, not a defect.
    let mut non_monotone = 0usize;
    let mut n_profiled = 0usize;
    for (i, prof) in dr.segment_profile_m.iter().enumerate() {
        let pts = &dr.rivers.segments[i].points;
        if prof.len() < 2 || prof.len() != pts.len() {
            continue;
        }
        n_profiled += 1;
        let in_lake = |p: (u32, u32)| dr.lake_map[p.1 as usize * w + p.0 as usize] != 0;
        if prof
            .windows(2)
            .enumerate()
            .any(|(j, p)| p[1] > p[0] + 1.0 && !(in_lake(pts[j]) && in_lake(pts[j + 1])))
        {
            non_monotone += 1;
        }
    }

    eprintln!(
        "LAKES  {} total ({exo} exo, {endo} endo, {crater} crater) | water {water_km2:.0} km² \
         | pre-breach detected {} | closed depressions on the RAW field {} ({:.0} km²)",
        dr.lakes.len(),
        hy.pre_lakes,
        hy.pre_pits,
        hy.pre_pits as f32 * cell_km2
    );
    eprintln!(
        "  INVARIANTS over {} lakes / {} footprint cells / {n_inlets} inlets / {n_profiled} profiles\n    \
         footprint above level      {v_above}\n    \
         footprint not connected    {v_split} (of which endorheic {split_endo})\n    \
         depth != level-floor (exported field) {v_depth_exported} | (raw pre-breach field) {v_depth_raw}\n    \
         area_km2 != footprint      {v_area}\n    \
         duplicate ids / empty footprints  {v_dup} / {v_empty}\n    \
         dangling lake_map ids      {v_dangling} ({v_dangling_belowsea} below-sea sub-4-cell, \
{v_dangling_detected} detected)\n    \
         level > inlet arrival      {v_inlet}\n    \
         exorheic without outlet    {missing_outlet}\n    \
         orphan mouths              {orphan_mouths} (of {n_termini} termini)\n    \
         non-monotone profiles      {non_monotone}",
        dr.lakes.len(),
        foot.values().map(|v| v.len()).sum::<usize>()
    );
    let n_lakes = dr.lakes.len() as u64;
    let n_foot = foot.values().map(|v| v.len()).sum::<usize>() as u64;
    sw.push("lakes total", dr.lakes.len() as f64);
    sw.push("lakes exorheic", exo as f64);
    sw.push("lakes endorheic", endo as f64);
    sw.push("lake water km2", water_km2);
    sw.push("pre-breach lakes detected", hy.pre_lakes as f64);
    sw.push("pre-breach closed depression cells", hy.pre_pits as f64);
    sw.checked("footprint above level", v_above as u64, n_foot);
    sw.checked("footprint not connected", v_split as u64, n_lakes);
    sw.checked("depth != level-floor (exported)", v_depth_exported as u64, n_lakes);
    sw.checked("depth != level-floor (raw)", v_depth_raw as u64, n_lakes);
    sw.checked("area_km2 != footprint", v_area as u64, n_lakes);
    sw.checked("duplicate lake ids", v_dup as u64, n_lakes);
    sw.checked("empty footprints", v_empty as u64, n_lakes);
    sw.checked("dangling lake_map ids (below-sea)", v_dangling_belowsea as u64, foot.len() as u64);
    sw.checked("dangling lake_map ids (detected)", v_dangling_detected as u64, foot.len() as u64);
    sw.checked("level > inlet arrival", v_inlet as u64, n_inlets as u64);
    sw.checked("exorheic without outlet", missing_outlet as u64, exo as u64);
    sw.checked("orphan mouths", orphan_mouths as u64, n_termini as u64);
    sw.checked("non-monotone profiles", non_monotone as u64, n_profiled as u64);

    // ── SHAPE METRICS ───────────────────────────────────────────────────────────────
    let win = (1000.0 / m_per_cell).round().max(1.0) as usize; // PHYSICAL 1 km window
    let (mut fr, mut slopes) = (Vec::new(), Vec::new());
    let step = if w > 4096 { 4 } else { 2 };
    for y in (win..h - win).step_by(step) {
        for x in (win..w - win).step_by(step) {
            let k = y * w + x;
            if field.data[k] <= SEA {
                continue;
            }
            let gx = 0.5 * (to_m(field.data[k + 1]) - to_m(field.data[k - 1]));
            let gy = 0.5 * (to_m(field.data[k + w]) - to_m(field.data[k - w]));
            slopes.push(((gx * gx + gy * gy).sqrt() / m_per_cell).atan().to_degrees());
            // floor / local ridge over the 1 km window
            let (mut lo, mut hi) = (f32::INFINITY, f32::NEG_INFINITY);
            for j in (y - win)..=(y + win) {
                for i in (x - win)..=(x + win) {
                    let a = to_m(field.data[j * w + i]);
                    lo = lo.min(a);
                    hi = hi.max(a);
                }
            }
            if hi > lo + 1.0 {
                fr.push((to_m(field.data[k]) - lo) / (hi - lo));
            }
        }
    }
    let n_s = slopes.len().max(1);
    let gt30 = 100.0 * slopes.iter().filter(|&&s| s > 30.0).count() as f32 / n_s as f32;
    let gt45 = 100.0 * slopes.iter().filter(|&&s| s > 45.0).count() as f32 / n_s as f32;
    let land_cells = field.data.iter().filter(|&&v| v > SEA).count();
    let chan_cells: usize = dr.rivers.segments.iter().map(|s| s.points.len()).sum();
    let (fr25, fr50, fr75) =
        (pct(&mut fr.clone(), 0.25), pct(&mut fr.clone(), 0.50), pct(&mut fr.clone(), 0.75));
    let smed = pct(&mut slopes.clone(), 0.50);
    let dd = chan_cells as f32 * km / (land_cells as f32 * cell_km2).max(1e-6);
    eprintln!(
        "SHAPE  floor/local-ridge (1 km window) p25 {fr25:.3} | p50 {fr50:.3} | p75 {fr75:.3}\n  \
         slope >30° {gt30:.2} % | >45° {gt45:.2} % | median {smed:.2}°\n  \
         drainage density {dd:.3} km/km² ({chan_cells} channel cells over {:.0} km² of land)",
        land_cells as f32 * cell_km2
    );
    sw.push("floor/ridge p50", fr50);
    sw.push("slope >30 pct", gt30);
    sw.push("slope >45 pct", gt45);
    sw.push("slope median deg", smed);
    sw.push("drainage density", dd);
    sw.push("channel cells", chan_cells as f64);
    let _ = n;

    // ── NETWORK: Strahler, VALLEY W/D per order, confluences, navigability ──────────
    // Spillways carry a MEANINGLESS Strahler order (`segment_kind` is authoritative), so they
    // are excluded from the hierarchy and from the cross-sections — a spillway's path is traced
    // over a col, outside the accumulation network, and has no valley.
    let mut by_order: HashMap<u8, Vec<usize>> = HashMap::new();
    let mut n_spillway = 0usize;
    for (i, s) in dr.rivers.segments.iter().enumerate() {
        if dr.segment_kind[i] == SegmentKind::Spillway {
            n_spillway += 1;
            continue;
        }
        by_order.entry(s.strahler_order).or_default().push(i);
    }
    let confluences = dr
        .rivers
        .segments
        .iter()
        .enumerate()
        .filter(|(i, s)| dr.segment_kind[*i] != SegmentKind::Spillway && s.upstream.len() >= 2)
        .count();
    let mut orders: Vec<u8> = by_order.keys().copied().collect();
    orders.sort_unstable();
    eprint!("NETWORK  confluences {confluences} | spillways {n_spillway} | Strahler");
    for o in &orders {
        eprint!("  S{o}:{}", by_order[o].len());
    }
    eprintln!();
    sw.push("confluences", confluences as f64);
    sw.push("segments", dr.rivers.segments.len() as f64);
    // An INVARIANT claim, not a quantity expected to respond: a channel-head law that changes
    // the DEPTH of the hierarchy has broken the network, and that is precisely how the
    // two-sided form failed (S5 59 → 0 at 8192², i.e. max order 5 → 4).
    sw.invariant("max strahler order", *orders.last().unwrap_or(&0) as f64);
    // S1..=S5 unconditionally, zero when absent — an order that COLLAPSES to nothing (the
    // two-sided law took S5 from 59 to 0 at 8192²) must show as a moving column, not vanish
    // from the report.
    for o in 1u8..=5 {
        sw.push(&format!("strahler S{o}"), by_order.get(&o).map(|v| v.len()).unwrap_or(0) as f64);
    }
    // Rule 10 caught this column TWICE, and the second catch changed the statistic. First
    // pass: `channel width p50` was 0.000 m at every order, both grids, both settings, because
    // MOST of the network carries a discharge of exactly 0.000 m³/s and `w = 5·Q^0.5` is then
    // exactly 0 — the median was measuring the dry majority. Second pass, after switching to
    // p90: still 0.000 for S5 at 8192², because at that grid **not one order-5 reach is wet**.
    //
    // So a percentile of a mostly-zero distribution is not a channel width at all; it is a
    // restatement of the dry fraction. The widths below are therefore CONDITIONAL on Q > 0 and
    // print `n/a` for an order with no wet reach — a printed `n/a` cannot be mistaken for a
    // measured zero, which a `0.000` column can and did.
    eprintln!(
        "  per order:  Q>0 share · channel width m p50/p90/max CONDITIONAL on Q>0  ·  VALLEY \
         width m / depth m / W-D (transect ⊥ flow, p50)"
    );
    let mut dry_orders: Vec<String> = Vec::new();
    for o in &orders {
        let mut widths: Vec<f32> = Vec::new();
        let mut vw: Vec<f32> = Vec::new();
        let mut vd: Vec<f32> = Vec::new();
        let mut wd: Vec<f32> = Vec::new();
        let mut n_wet = 0usize;
        for &i in &by_order[o] {
            widths.push(dr.segment_width_m[i]);
            if dr.segment_discharge_m3s[i] > 0.0 {
                n_wet += 1;
            }
            let pts = &dr.rivers.segments[i].points;
            if pts.len() < 3 {
                continue;
            }
            let mid = pts.len() / 2;
            let (ax, ay) = pts[mid - 1];
            let (bx, by) = pts[mid + 1];
            let (dx, dy) = (bx as f32 - ax as f32, by as f32 - ay as f32);
            let (x, y) = (pts[mid].0 as usize, pts[mid].1 as usize);
            if let Some((wm, dm)) = valley_cross_section(field, &to_m, x, y, dx, dy, m_per_cell) {
                vw.push(wm);
                vd.push(dm);
                wd.push(wm / dm);
            }
        }
        let n_o = by_order[o].len();
        let wet_pct = 100.0 * n_wet as f32 / n_o.max(1) as f32;
        let mut wet: Vec<f32> = widths.iter().copied().filter(|&x| x > 0.0).collect();
        let ch = if wet.is_empty() {
            dry_orders.push(format!("S{o} ({n_o} reaches)"));
            String::from("      n/a /     n/a /      n/a")
        } else {
            format!(
                "{:9.3} / {:7.3} / {:8.3}",
                pct(&mut wet.clone(), 0.5),
                pct(&mut wet.clone(), 0.9),
                pct(&mut wet, 1.0)
            )
        };
        let wmax = widths.iter().copied().fold(0.0f32, f32::max);
        eprintln!(
            "    S{o}  n {n_o:5}  Q>0 {wet_pct:5.1} %  channel w {ch} m   valley w {:7.0} m  \
             d {:6.1} m  W/D {:6.1}  (measured on {} of {n_o})",
            pct(&mut vw.clone(), 0.5),
            pct(&mut vd.clone(), 0.5),
            pct(&mut wd.clone(), 0.5),
            wd.len()
        );
        // Registered: the wet SHARE and the MAXIMUM width. Both are defined for a dry order
        // (0 % and 0 m) yet still respond to the law, so neither can go flat-at-zero silently
        // the way a conditional percentile can go undefined.
        sw.push(&format!("Q>0 share S{o}"), wet_pct);
        sw.push(&format!("channel width max m S{o}"), wmax);
        sw.push(&format!("valley W/D S{o}"), pct(&mut wd.clone(), 0.5));
        sw.push(&format!("valley depth m S{o}"), pct(&mut vd.clone(), 0.5));
    }
    eprintln!(
        "  ORDERS WITH NO WET REACH: {}",
        if dry_orders.is_empty() { String::from("none") } else { dry_orders.join(", ") }
    );
    let mut nav: HashMap<&str, usize> = HashMap::new();
    for kind in &dr.segment_navigability {
        *nav.entry(match kind {
            Navigability::NonNavigable => "non-nav",
            Navigability::SmallBoat => "boat",
            Navigability::Barge => "barge",
            Navigability::Ship => "ship",
        })
        .or_insert(0) += 1;
    }
    let mut q: Vec<f32> = dr.segment_discharge_m3s.clone();
    eprintln!(
        "  NAVIGABILITY  non-nav {} | boat {} | barge {} | ship {}   ·   discharge m³/s p50 \
         {:.3} p90 {:.3} max {:.3}",
        nav.get("non-nav").copied().unwrap_or(0),
        nav.get("boat").copied().unwrap_or(0),
        nav.get("barge").copied().unwrap_or(0),
        nav.get("ship").copied().unwrap_or(0),
        pct(&mut q.clone(), 0.5),
        pct(&mut q.clone(), 0.9),
        pct(&mut q, 1.0)
    );
    sw.pinned(
        "navigable boat reaches",
        nav.get("boat").copied().unwrap_or(0) as f64,
        "the navigability classes are thresholds in m³/s and the discharge sits 2-3 orders \
         below them (Finding 47); the classification cannot respond to a channel-head law, \
         only to the hypsometry",
    );
    sw.push("discharge p90", pct(&mut dr.segment_discharge_m3s.clone(), 0.9));
    sw.push("discharge max", pct(&mut dr.segment_discharge_m3s.clone(), 1.0));

    // ── RULE 7 — the control block: the incision must still be ALIVE ────────────────
    // Not `invariant()`: this quantity is EXPECTED to move a little. It is the block that
    // caught `A_c × 100`, where the incision stopped happening and the hypsometry walked back
    // to the un-eroded field. What must not happen is a LARGE move.
    // BOTH fields. The other benches (`coastal_comb_levers`) report the hypsometry of the RAW
    // eroded field; this one exports the BREACHED field. Reporting both here keeps the two
    // benches comparable instead of leaving a few metres of difference unexplained, and it
    // prices the breach as a measured quantity.
    let alts: Vec<f32> = field.data.iter().filter(|&&v| v > SEA).map(|&v| to_m(v)).collect();
    let raw_alts: Vec<f32> = raw.data.iter().filter(|&&v| v > SEA).map(|&v| to_m(v)).collect();
    // f64 ACCUMULATION, deliberately. An f32 running sum over the 1.8 M land altitudes at
    // 8192² saturates — once the partial sum passes ~10⁹, adding 700 m stops changing it — and
    // it read 672.2 m where the f64 mean is 685.0. `coastal_comb_levers` already sums in f64,
    // so an f32 sum here would also have put the two benches 13 m apart for no physical reason.
    let mean = (alts.iter().map(|&v| v as f64).sum::<f64>() / alts.len().max(1) as f64) as f32;
    let raw_mean =
        (raw_alts.iter().map(|&v| v as f64).sum::<f64>() / raw_alts.len().max(1) as f64) as f32;
    let land_pct = 100.0 * land_cells as f32 / n as f32;
    let raw_land = raw.data.iter().filter(|&&v| v > SEA).count();
    eprintln!(
        "RULE 7  hypsometry BREACHED mean {mean:.1} m p50 {:.0} p90 {:.0} | land {land_pct:.2} % \
         ({:.0} km²)\n        hypsometry RAW      mean {raw_mean:.1} m p50 {:.0} p90 {:.0} | \
         land {:.2} % ({:.0} km²)  ⇒ the breach costs {:.2} m of mean altitude over {} cells",
        pct(&mut alts.clone(), 0.5),
        pct(&mut alts.clone(), 0.9),
        land_cells as f32 * cell_km2,
        pct(&mut raw_alts.clone(), 0.5),
        pct(&mut raw_alts.clone(), 0.9),
        100.0 * raw_land as f32 / n as f32,
        raw_land as f32 * cell_km2,
        raw_mean - mean,
        raw_land as i64 - land_cells as i64
    );
    sw.push("hypsometry mean m (breached)", mean);
    sw.push("hypsometry mean m (raw)", raw_mean);
    sw.push("hypsometry p50 m (breached)", pct(&mut alts.clone(), 0.5));
    sw.push("land pct", land_pct);
}

#[test]
#[ignore]
fn channel_head_law_invariants() {
    let ss = SteinSteinParams::default();
    eprintln!(
        "\n==========  CHANNEL-HEAD LAW — INVARIANT SUITE (arid-hot {LAT}°/span {SPAN})  =========="
    );
    // One sweep PER RESOLUTION: the two grids legitimately differ on almost every quantity
    // (that IS the Findings 42–44 defect), so pooling them would let a resolution difference
    // stand in for the law's effect and rule 10 would pass on nothing.
    let mut sweeps: Vec<(usize, Sweep)> = Vec::new();
    for target in [2048usize, 8192] {
        eprintln!("\n╔══════ {target}² ══════╗");
        let mut sw = Sweep::new(format!("channel-head law @ {target}²"));
        for law in [false, true] {
            let raw = terrain(target, law);
            let hy = hydro(&raw, &ss);
            report(
                if law { "law ON" } else { "shipped (law off)" },
                target,
                &raw,
                &hy,
                &ss,
                &mut sw,
            );
        }
        sweeps.push((target, sw));
    }
    // Rule 10 LAST, so the whole report is on screen when it fires, and over BOTH resolutions
    // so one grid's failure does not hide the other's.
    eprintln!("\n══════ RULE 10 — did every reported quantity respond? ══════");
    let mut bad = 0usize;
    for (target, sw) in &sweeps {
        let f = sw.failures();
        if f.is_empty() {
            eprintln!("  {target}²: OK");
        } else {
            bad += f.len();
            eprintln!("  {target}²: {} column(s) measured NOTHING", f.len());
            for (name, why) in f {
                eprintln!("    · {name}: {why}");
            }
        }
    }
    assert_eq!(bad, 0, "rule 10: {bad} reported quantity/quantities measured nothing (above)");
}
