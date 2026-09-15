//! ADR 0001 Finding 85 blocks A2/A3 — the level/outlet/balance fixed point at 8192², gate OFF
//! and ON, two climate beds, base-level bound ON (production).
//!
//! Beyond the Finding 74 table this carries the two guards Finding 39 earned and the round did not
//! list: **MAX below-sea level** and **footprint above level** (the over-flood net). Finding 39's
//! first cut of the same law produced MAX 613 m and claimed/valid 1.952×, and the fix was to stop
//! absorbing neighbouring regions. This gate absorbs exactly one neighbour — the one the level
//! rule already calls the same water body — and whether that stays bounded is the question.
//!
//! Run: cargo test -p ymir-core --release --test f85_relevel -- --ignored --nocapture

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
    C1DrainageConfig, DrainageClimate, MergedUnionRelevel, SegmentKind, SegmentRow, Spillway,
    apply_lake_water_balance, below_sea_basin_lakes_infil, c1_drainage_windowed,
    clip_rivers_to_lakes, exorheic_lakes_missing_outlet, potential_evaporation_mm,
    runoff_km2_to_m3s,
};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{
    D8_DX, D8_DY, DIR_NONE, FlowConfig, RiverSegment, breach_monotone, compute_flow,
};

const DOMAIN_KM: f32 = 400.0;
const RATIO: f32 = 7.5;

fn fnv(g: &GridF32) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for &v in &g.data {
        h ^= v.to_bits() as u64;
        h = h.wrapping_mul(0x100_0000_01b3);
    }
    h
}

/// Everything the round wants compared between the two gate states, for one (bed, gate).
struct Pass {
    budget: f32,
    wc_n: usize,
    wc_q: f64,
    n_ocean: usize,
    sp_ocean: f64,
    n_chain: usize,
    sp_chain: f64,
    terminal: f64,
    hole: f64,
    evap: f64,
    /// id → (level_m, area_km2, depth_m, footprint cells)
    lakes: HashMap<u32, (f32, f32, f32, usize)>,
    /// cells with a below-sea id, for the drowned-terrain diff
    footprint: Vec<bool>,
    /// the below-sea id per cell — ids are assigned by SCAN ORDER and are NOT comparable between
    /// gate states (Finding 81). Kept so basins can be paired by their FLOOR CELL instead.
    id_map: Vec<u32>,
    /// id → its lowest cell, the geometric identity of the basin
    floor_of: HashMap<u32, usize>,
    max_level_m: f32,
    max_depth_m: f32,
    above_level: usize,
    exorheic_no_outlet: usize,
    max_wc_q: f32,
    max_wc_cells: f32,
    max_sp_q: f32,
    max_sp_cells: f32,
    widths_p50: f32,
    widths_max: f32,
    over_one_cell: f64,
    term: ymir_core::tectonics_c1::drainage::SpillwayTermination,
    /// ADR Finding 85 — the per-basin mass balance: `(id, inflow, evap, out, exorheic, has_sill,
    /// a_eq, a_spill)`. The round's "chaque m³/s résiduel porte sa cuvette et sa raison".
    balance: Vec<(u32, f32, f32, f32, bool, bool, f32, f32)>,
}

#[allow(clippy::too_many_arguments)]
fn one_pass(
    raw: &GridF32,
    field: &GridF32,
    flow: &ymir_core::terrain::flow::FlowResult,
    wc: &[u8],
    dcfg: &C1DrainageConfig,
    ss: &SteinSteinParams,
    dclim: &DrainageClimate<'_>,
    surplus: &[f32],
    land: &[bool],
    n2m: f32,
) -> Pass {
    let (w, h) = (field.width, field.height);
    let n = w * h;
    let cell_km2 = CELL_KM * CELL_KM;
    let budget = runoff_km2_to_m3s(
        (0..n).filter(|&k| land[k]).map(|k| surplus[k] as f64).sum::<f64>() as f32,
    );
    let pre = c1_drainage_windowed(raw, None, dcfg, ss, DOMAIN_KM);
    let mut dr = c1_drainage_windowed(field, Some(dclim), dcfg, ss, DOMAIN_KM);
    dr.lakes = pre.lakes;
    dr.lake_map = pre.lake_map;
    let lakes_in = std::mem::take(&mut dr.lakes);
    dr.lakes = apply_lake_water_balance(
        field,
        &dr.flow,
        dclim,
        cell_km2,
        ss,
        &lakes_in,
        &mut dr.lake_map,
        None,
        w,
        h,
    );
    let premerge = dr.lake_map.clone();
    let bs = below_sea_basin_lakes_infil(field, dclim, dcfg, ss, DOMAIN_KM, Some(&premerge), None);
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

    // ── the Finding 74 table ──────────────────────────────────────────────────
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

    // ── the lakes, their levels, and the over-flood net ───────────────────────
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
    let mut lakes: HashMap<u32, (f32, f32, f32, usize)> = HashMap::new();
    let (mut max_level, mut max_depth, mut above) = (f32::MIN, f32::MIN, 0usize);
    for lk in &dr.lakes {
        if lk.base.id < 1_000_001 {
            continue;
        }
        let cells = fc.get(&lk.base.id).copied().unwrap_or(0);
        lakes.insert(lk.base.id, (lk.level_m, lk.area_km2, lk.depth_m, cells));
        max_level = max_level.max(lk.level_m);
        max_depth = max_depth.max(lk.depth_m);
        if fmax.get(&lk.base.id).copied().unwrap_or(f32::MIN) > SEA + lk.level_m / n2m + 1e-6 {
            above += 1;
        }
    }
    let footprint: Vec<bool> = (0..n).map(|k| bs.lake_map[k] != 0).collect();
    let mut floor_of: HashMap<u32, usize> = HashMap::new();
    for k in 0..n {
        let id = bs.lake_map[k];
        if id != 0 {
            let e = floor_of.entry(id).or_insert(k);
            if field.data[k] < field.data[*e] {
                *e = k;
            }
        }
    }

    // ── widths ───────────────────────────────────────────────────────────────
    let idx: Vec<usize> = (0..dr.rivers.segments.len())
        .filter(|&i| dr.segment_kind[i] == SegmentKind::Watercourse)
        .collect();
    let (mut max_wc_q, mut max_wc_cells) = (0.0f32, 0.0f32);
    for &i in &idx {
        if dr.segment_discharge_m3s[i] > max_wc_q {
            max_wc_q = dr.segment_discharge_m3s[i];
            max_wc_cells = dr.segment_width_m[i] * RATIO / 48.83;
        }
    }
    let mut wcell: Vec<f32> = idx.iter().map(|&i| dr.segment_width_m[i] * RATIO / 48.83).collect();
    wcell.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let (mut max_sp_q, mut max_sp_cells) = (0.0f32, 0.0f32);
    for sw in &bs.spillways {
        if sw.discharge_m3s > max_sp_q {
            max_sp_q = sw.discharge_m3s;
            max_sp_cells = sw.width_m * RATIO / 48.83;
        }
    }

    Pass {
        budget,
        wc_n,
        wc_q,
        n_ocean,
        sp_ocean,
        n_chain,
        sp_chain,
        terminal,
        hole: budget as f64 - terminal,
        evap: bs.basins.iter().map(|b| b.evaporation_m3s as f64).sum(),
        lakes,
        footprint,
        id_map: bs.lake_map.clone(),
        floor_of,
        max_level_m: max_level,
        max_depth_m: max_depth,
        above_level: above,
        exorheic_no_outlet: exorheic_lakes_missing_outlet(&dr).len(),
        max_wc_q,
        max_wc_cells,
        max_sp_q,
        max_sp_cells,
        widths_p50: pct(&wcell, 0.50),
        widths_max: wcell.last().copied().unwrap_or(0.0),
        over_one_cell: 100.0 * wcell.iter().filter(|&&x| x > 1.0).count() as f64
            / wcell.len().max(1) as f64,
        term: bs.termination.clone(),
        balance: {
            let mut out_of: HashMap<u32, f32> = HashMap::new();
            for sw in &bs.spillways {
                *out_of.entry(sw.lake_id).or_default() += sw.discharge_m3s;
            }
            bs.basins
                .iter()
                .map(|b| {
                    (
                        b.id,
                        b.inflow_m3s,
                        b.evaporation_m3s,
                        out_of.get(&b.id).copied().unwrap_or(0.0),
                        b.exorheic,
                        b.has_sill,
                        b.a_eq_km2,
                        b.a_spill_km2,
                    )
                })
                .collect()
        },
    }
}

#[test]
#[ignore]
fn f85_relevel() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 85 · A2/A3 — the level fixed point, gate OFF vs ON  =====");

    let mut off = C1DrainageConfig::default();
    off.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    off.thresholds.full_tree = false;
    let mut on = off.clone();
    on.merged_union_relevel = Some(MergedUnionRelevel::default());

    let raw = build_field(Knobs::shipped());
    let (w, h) = (raw.width, raw.height);
    let n = w * h;
    let hl = metric_height_u16(&raw, &ss);
    let u16_step = (hl.max_m - hl.min_m) / 65535.0;
    eprintln!(
        "\n── the field, which the gate cannot reach ──\n   eroded FNV-1a {:#018x} | u16 step \
         {u16_step:.4} m\n   `merged_union_relevel` is a field of `C1DrainageConfig`; the upscale \
         has no path to it, so this is the SAME field object for both states.",
        fnv(&raw)
    );
    let pre0 = c1_drainage_windowed(&raw, None, &off, &ss, DOMAIN_KM);
    let field = breach_monotone(&raw, &pre0.flow.filled, &pre0.lake_map, SEA, w, h);
    let flow = compute_flow(&field, &FlowConfig { sea_level: SEA, ..Default::default() });
    let wc = water_class(&field, SEA);
    let land: Vec<bool> = (0..n).map(|k| field.data[k] > SEA).collect();
    let ocean_cells = wc.iter().filter(|&&c| c == 1).count();
    eprintln!(
        "   breached FNV-1a {:#018x} | land {} | OCEAN cells {ocean_cells} (the coastal mask the \
         round asks not to move — it is derived from `field`, so it cannot)",
        fnv(&field),
        land.iter().filter(|&&l| l).count()
    );

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
        eprintln!(
            "\n╔════════ {bed} ════════╗  climate FNV precip {:#018x}",
            fnv(&climate.precipitation)
        );

        let a = one_pass(&raw, &field, &flow, &wc, &off, &ss, &dclim, &surplus, &land, n2m);
        let b = one_pass(&raw, &field, &flow, &wc, &on, &ss, &dclim, &surplus, &land, n2m);

        eprintln!("   {:<34} {:>16} {:>16}", "the Finding 74 table", "GATE OFF", "**GATE ON**");
        let row = |nm: &str, x: String, y: String| eprintln!("   {nm:<34} {x:>16} {y:>16}");
        row("budget m³/s", format!("{:.1}", a.budget), format!("{:.1}", b.budget));
        row(
            "Watercourse → sea  n / m³/s",
            format!("{} / {:.1}", a.wc_n, a.wc_q),
            format!("{} / {:.1}", b.wc_n, b.wc_q),
        );
        row(
            "Spillway → ocean   n / m³/s",
            format!("{} / {:.1}", a.n_ocean, a.sp_ocean),
            format!("{} / {:.1}", b.n_ocean, b.sp_ocean),
        );
        row(
            "Spillway → chained n / m³/s",
            format!("{} / {:.1}", a.n_chain, a.sp_chain),
            format!("{} / {:.1}", b.n_chain, b.sp_chain),
        );
        row(
            "TERMINAL (× budget)",
            format!("{:.1} (×{:.3})", a.terminal, a.terminal / a.budget as f64),
            format!("{:.1} (×{:.3})", b.terminal, b.terminal / b.budget as f64),
        );
        row("**HOLE**", format!("{:.1}", a.hole), format!("{:.1}", b.hole));
        row("Σ evaporation m³/s", format!("{:.2}", a.evap), format!("{:.2}", b.evap));
        row(
            "own body   n / m³/s",
            format!("{} / {:.2}", a.term.to_own_body, a.term.q_own_body_m3s),
            format!("{} / {:.2}", b.term.to_own_body, b.term.q_own_body_m3s),
        );
        row(
            "paths ENDING at a lake m³/s",
            format!("{:.2}", a.term.q_detected_lake_m3s),
            format!("{:.2}", b.term.q_detected_lake_m3s),
        );
        row(
            "of which really DROPPED m³/s",
            format!("{:.2}", a.term.q_lake_unrouted_m3s),
            format!("{:.2}", b.term.q_lake_unrouted_m3s),
        );
        row(
            "routed THROUGH a lake m³/s",
            format!("{:.2}", a.term.q_routed_through_lake_m3s),
            format!("{:.2}", b.term.q_routed_through_lake_m3s),
        );
        row(
            "passes / class merges / multi",
            format!(
                "{} / {} / {}",
                a.term.passes_used, a.term.classes_merged, a.term.multi_region_classes
            ),
            format!(
                "{} / {} / {}",
                b.term.passes_used, b.term.classes_merged, b.term.multi_region_classes
            ),
        );
        row(
            "pairs merged AFTER the loop",
            format!("{}", a.term.pairs_merged),
            format!("{}", b.term.pairs_merged),
        );
        row("basins", format!("{}", a.lakes.len()), format!("{}", b.lakes.len()));
        row(
            "MAX below-sea level / depth m",
            format!("{:.1} / {:.1}", a.max_level_m, a.max_depth_m),
            format!("{:.1} / {:.1}", b.max_level_m, b.max_depth_m),
        );
        row("🛑 footprint ABOVE level", format!("{}", a.above_level), format!("{}", b.above_level));
        row(
            "🛑 exorheic WITHOUT outlet",
            format!("{}", a.exorheic_no_outlet),
            format!("{}", b.exorheic_no_outlet),
        );
        row(
            "max Watercourse m³/s / cells",
            format!("{:.3} / {:.2}", a.max_wc_q, a.max_wc_cells),
            format!("{:.3} / {:.2}", b.max_wc_q, b.max_wc_cells),
        );
        row(
            "max Spillway m³/s / cells",
            format!("{:.2} / {:.2}", a.max_sp_q, a.max_sp_cells),
            format!("{:.2} / {:.2}", b.max_sp_q, b.max_sp_cells),
        );
        row(
            "widths p50 / max / >1 cell %",
            format!("{:.3} / {:.2} / {:.1}", a.widths_p50, a.widths_max, a.over_one_cell),
            format!("{:.3} / {:.2} / {:.1}", b.widths_p50, b.widths_max, b.over_one_cell),
        );

        // ── A3 · what moved, basin by basin ───────────────────────────────────
        let drowned: Vec<usize> = (0..n).filter(|&k| b.footprint[k] && !a.footprint[k]).collect();
        let drowned_land = drowned.iter().filter(|&&k| land[k]).count();
        let freed = (0..n).filter(|&k| a.footprint[k] && !b.footprint[k]).count();
        eprintln!(
            "\n   A3: footprint cells {} → {} | NEW cells {} of which LAND (drowned terrain) \
             **{drowned_land}** = **{:.1} km²** | cells freed {freed}",
            a.footprint.iter().filter(|&&x| x).count(),
            b.footprint.iter().filter(|&&x| x).count(),
            drowned.len(),
            drowned_land as f32 * cell_km2
        );
        // ⚠️ PAIRING BY ID IS INVALID HERE. Below-sea ids are assigned by scan order, and a
        // class merge changes the count, so every id after it shifts — Finding 81 withdrew a
        // claim over exactly this. Basins are paired GEOMETRICALLY instead: an OFF basin maps to
        // whichever ON basin covers its FLOOR CELL, which survives a merge (two OFF basins map to
        // one ON basin) and is what "the same place" means.
        let mut rises: Vec<f32> = Vec::new();
        let mut moved = 0usize;
        let mut lost = 0usize;
        let mut rows: Vec<(u32, u32, f32, f32, f32, f32, usize)> = Vec::new();
        for (&id_a, &floor) in &a.floor_of {
            let Some(&(la, aa, _, _ca)) = a.lakes.get(&id_a) else { continue };
            let id_b = b.id_map[floor];
            if id_b == 0 {
                lost += 1;
                continue; // the floor is no longer under water at all
            }
            let Some(&(lb, ab, _, cb)) = b.lakes.get(&id_b) else { continue };
            if (lb - la).abs() > 1e-4 {
                moved += 1;
                rises.push(lb - la);
            }
            rows.push((id_a, id_b, la, lb, aa, ab, cb));
        }
        let rs = sorted(rises.clone());
        eprintln!(
            "   A3: basins paired by FLOOR CELL {} (floors no longer submerged {lost}) | levels              moved on **{moved}** | rise p50 **{:.2} m** p90 **{:.2} m** max **{:.2} m**",
            rows.len(),
            pct(&rs, 0.50),
            pct(&rs, 0.90),
            rs.last().copied().unwrap_or(0.0)
        );
        eprintln!(
            "   A3: distinct ON basins those floors land in: {} of {} (a drop means floors merged)",
            rows.iter().map(|r| r.1).collect::<HashSet<u32>>().len(),
            b.lakes.len()
        );
        // ── where the water goes, basin by basin (the round's last line) ──────
        for (tag, p) in [("GATE OFF", &a), ("**GATE ON**", &b)] {
            let mut bal = p.balance.clone();
            bal.sort_by(|x, y| {
                (y.1 - y.2 - y.3)
                    .partial_cmp(&(x.1 - x.2 - x.3))
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            let tot_in: f64 = bal.iter().map(|r| r.1 as f64).sum();
            let tot_ev: f64 = bal.iter().map(|r| r.2 as f64).sum();
            let tot_out: f64 = bal.iter().map(|r| r.3 as f64).sum();
            eprintln!(
                "
   BALANCE {tag}: {} basins | Σ inflow {tot_in:.1} | Σ evap {tot_ev:.2} | Σ                  out {tot_out:.1} | Σ RESIDUAL **{:.1} m³/s**",
                bal.len(),
                tot_in - tot_ev - tot_out
            );
            eprintln!(
                "      {:>9} {:>10} {:>9} {:>9} {:>10} {:>6} {:>5} {:>11} {:>11}",
                "id",
                "inflow",
                "evap",
                "out",
                "RESIDUAL",
                "exorh",
                "sill",
                "a_eq km²",
                "a_spill km²"
            );
            for (id, inf, ev, o, ex, hs, aeq, asp) in bal.iter().take(6) {
                eprintln!(
                    "      {id:>9} {inf:>10.2} {ev:>9.2} {o:>9.2} {:>10.2} {:>6} {:>5} {:>11.1}                      {asp:>11.3}",
                    inf - ev - o,
                    ex,
                    hs,
                    if aeq.is_finite() { *aeq } else { f32::MAX }
                );
            }
        }
        rows.sort_by(|x, y| (y.3 - y.2).partial_cmp(&(x.3 - x.2)).unwrap());
        eprintln!(
            "   {:>9} {:>9} {:>10} {:>10} {:>9} {:>10} {:>11} {:>11} {:>10}",
            "id OFF",
            "id ON",
            "level OFF",
            "level ON",
            "Δ m",
            "Δ u16",
            "area OFF",
            "area ON",
            "cells ON"
        );
        for (ia, ib, la, lb, aa, ab, cb) in rows.iter().take(8) {
            eprintln!(
                "   {ia:>9} {ib:>9} {la:>10.2} {lb:>10.2} {:>9.2} {:>10.1} {aa:>11.3} {ab:>11.3} {cb:>10}",
                lb - la,
                (lb - la) / u16_step
            );
        }
    }
}
