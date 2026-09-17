//! ADR 0001 Finding 91 — why an interior depression escapes the lake tag. ATTRIBUTION ONLY.
//!
//! Two physically similar interior basins: one is a correct through-flow lake (#1000010), the other
//! reads −1 m with an empty drainage field, a terrestrial biome, and a river falling into it under
//! the name "sub-sea evaporative sink". This measures both at four stages and tests the three
//! candidates the round named plus a fourth of my own.
//!
//! No production change. Run:
//!   cargo test -p ymir-core --release --test f91_orphan -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, SEA, build_field, pct, sorted};
use std::collections::{HashMap, HashSet};
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, apply_lake_water_balance, below_sea_basin_lakes_infil,
    c1_drainage_windowed, runoff_accumulation, runoff_km2_to_m3s,
};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{D8_DX, D8_DY, breach_monotone};

const DOMAIN_KM: f32 = 400.0;

#[test]
#[ignore]
fn f91_orphan() {
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 91 · the untagged interior depression  ==========");

    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    let lake_min_cells = (on.lake_min_area_km2 / cell_km2).ceil().max(1.0) as usize;
    eprintln!(
        "   `lake_min_area_km2` = **{}** ⇒ **{lake_min_cells} cells** at 8192² (the `detect_lakes` \
         filter). ⚠️ the below-sea path discards it: `let _ = cfg.lake_min_area_km2;`",
        on.lake_min_area_km2
    );

    // ── the four stages ─────────────────────────────────────────────────────
    let tect = build_field(Knobs { no_incision: true, bathymetry_off: true, ..Knobs::default() });
    let (w, h) = (tect.width, tect.height);
    let n = w * h;
    let ero_pb = build_field(Knobs { bathymetry_off: true, ..Knobs::passes(2) });
    let raw = build_field(Knobs::passes(2)); // eroded + bathymetry = the pipeline's output
    let d0 = c1_drainage_windowed(&raw, None, &on, &ss, DOMAIN_KM);
    let field = breach_monotone(&raw, &d0.flow.filled, &d0.lake_map, SEA, w, h);
    let m = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);

    // ── the chain, to get the tagging exactly as production does ────────────
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
    let racc = runoff_accumulation(&field, &dr.flow, &dclim, cell_km2, None, None, w, h);

    let wc_t = water_class(&tect, SEA);
    let wc_e = water_class(&ero_pb, SEA);
    let wc_r = water_class(&raw, SEA);
    let wc_f = water_class(&field, SEA);

    // ── the enclosed below-sea components of the DELIVERED field ────────────
    let mut comp_id = vec![0u32; n];
    let mut comps: Vec<(u32, Vec<usize>)> = Vec::new();
    let mut next = 1u32;
    for s in 0..n {
        if wc_f[s] != 2 || comp_id[s] != 0 {
            continue;
        }
        let mut stack = vec![s];
        comp_id[s] = next;
        let mut cells = Vec::new();
        while let Some(k) = stack.pop() {
            cells.push(k);
            // ⚠️ production grows the component with `nb4`, which is BOUNDS-CHECKED and NOT
            // periodic (`region_of`, Finding 40). Matching it here or the populations differ.
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for dd in 0..8 {
                let (nx, ny) = (x + D8_DX[dd], y + D8_DY[dd]);
                if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                    continue;
                }
                let nk = ny as usize * w + nx as usize;
                if wc_f[nk] == 2 && comp_id[nk] == 0 {
                    comp_id[nk] = next;
                    stack.push(nk);
                }
            }
        }
        comps.push((next, cells));
        next += 1;
    }
    eprintln!(
        "\n── the enclosed below-sea (`wc == 2`) components of the delivered field: **{}** ──",
        comps.len()
    );

    // the LOCAL inflow exactly as the below-sea path reads it: max runoff at the LAND 4-neighbours
    let local_inflow = |cells: &[usize]| -> f32 {
        let mut best = 0.0f32;
        for &k in cells {
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for (dx, dy) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
                let nx = (x + dx).rem_euclid(w as i32) as usize;
                let ny = (y + dy).rem_euclid(h as i32) as usize;
                let nk = ny * w + nx;
                if field.data[nk] > SEA {
                    best = best.max(racc[nk]);
                }
            }
        }
        best
    };

    // ── C · the count: components with NO lake tag at all ───────────────────
    let mut orphans: Vec<(u32, usize, usize, f32)> = Vec::new(); // id, size, floor, inflow
    let mut tagged = 0usize;
    for (id, cells) in &comps {
        let covered = cells.iter().any(|&k| bs.lake_map[k] != 0 || detected[k] != 0);
        let floor = cells
            .iter()
            .copied()
            .min_by(|&a, &b| field.data[a].total_cmp(&field.data[b]))
            .unwrap_or(0);
        if covered {
            tagged += 1;
        } else {
            orphans.push((*id, cells.len(), floor, runoff_km2_to_m3s(local_inflow(cells))));
        }
    }
    orphans.sort_by(|a, b| b.1.cmp(&a.1));
    eprintln!(
        "── C · **{} components carry a lake tag · {} carry NONE** (of {} total) · below-sea \
         bodies inventoried: {}",
        tagged,
        orphans.len(),
        comps.len(),
        bs.lakes.len()
    );
    let with_flow: Vec<&(u32, usize, usize, f32)> = orphans.iter().filter(|o| o.3 >= 1.0).collect();
    eprintln!(
        "   of the untagged, **{} have a LOCAL inflow ≥ 1 m³/s** and {} have < 1 · their sizes: \
         p50 {:.0} cells, max {} cells ({:.3} km²)",
        with_flow.len(),
        orphans.len() - with_flow.len(),
        pct(&sorted(orphans.iter().map(|o| o.1 as f32).collect()), 0.50),
        orphans.first().map(|o| o.1).unwrap_or(0),
        orphans.first().map(|o| o.1).unwrap_or(0) as f32 * cell_km2
    );
    // which spillways / rivers land on an untagged component
    let mut q_on_orphan = 0.0f64;
    let mut n_sw = 0usize;
    for sw in &bs.spillways {
        let &(lx, ly) = sw.points.last().unwrap();
        let k = ly as usize * w + lx as usize;
        if wc_f[k] == 2 && bs.lake_map[k] == 0 && detected[k] == 0 {
            n_sw += 1;
            q_on_orphan += sw.discharge_m3s as f64;
            eprintln!(
                "   **a spillway terminates on an untagged component**: basin {} → ({lx},{ly}) \
                 at {:.3} m, {:.1} m³/s, component {} of {} cells",
                sw.lake_id,
                m(&field, k),
                sw.discharge_m3s,
                comp_id[k],
                comps.iter().find(|c| c.0 == comp_id[k]).map(|c| c.1.len()).unwrap_or(0)
            );
        }
    }
    let mut n_riv = 0usize;
    let mut q_riv = 0.0f64;
    for (i, s) in dr.rivers.segments.iter().enumerate() {
        if s.downstream.is_some() {
            continue;
        }
        let &(lx, ly) = s.points.last().unwrap();
        let k = ly as usize * w + lx as usize;
        if wc_f[k] == 2 && bs.lake_map[k] == 0 && detected[k] == 0 {
            let q = dr.segment_discharge_m3s.get(i).copied().unwrap_or(0.0);
            if q >= 1.0 {
                n_riv += 1;
                q_riv += q as f64;
            }
        }
    }
    eprintln!(
        "   terminating there: **{n_sw} spillways carrying {q_on_orphan:.1} m³/s** and \
         **{n_riv} river mouths ≥ 1 m³/s carrying {q_riv:.1} m³/s** — these are the microscope's \
         \"sub-sea evaporative sinks\""
    );

    // ── A · the two basins, at four stages ──────────────────────────────────
    let show = |label: &str, k: usize| {
        eprintln!(
            "      {label:<34} floor ({:>4},{:>4}) · TECTONIC (pre-bath) **{:>9.2} m** · ERODED \
             (pre-bath) **{:>9.2} m** · ERODED+BATH (pipeline out) **{:>8.2} m** · BREACHED \
             (delivered) **{:>8.2} m**",
            k % w,
            k / w,
            m(&tect, k),
            m(&ero_pb, k),
            m(&raw, k),
            m(&field, k)
        );
        eprintln!(
            "      {:<34} water_class {} / {} / {} / {} · lake_map: detected {} · below-sea {}",
            "", wc_t[k], wc_e[k], wc_r[k], wc_f[k], detected[k], bs.lake_map[k]
        );
    };
    eprintln!("\n── A · the two basins, at the four stages, at their FLOOR cell ──");
    // the correct one: the biggest inventoried below-sea body
    if let Some(best) = bs.lakes.iter().max_by(|a, b| a.area_km2.total_cmp(&b.area_km2)) {
        let cells: Vec<usize> = (0..n).filter(|&k| bs.lake_map[k] == best.base.id).collect();
        let floor = cells
            .iter()
            .copied()
            .min_by(|&a, &b| field.data[a].total_cmp(&field.data[b]))
            .unwrap_or(0);
        eprintln!(
            "   TAGGED, the largest inventoried below-sea body: id **{}** · {:.1} km² \
             ({} cells) · level {:.1} m · depth {:.1} m · type {:?}",
            best.base.id,
            best.area_km2,
            cells.len(),
            best.level_m,
            best.depth_m,
            best.lake_type
        );
        show("id", floor);
    }
    for (i, (id, size, floor, q)) in orphans.iter().take(3).enumerate() {
        eprintln!(
            "   UNTAGGED #{i}: component {id} · {size} cells ({:.4} km²) · LOCAL inflow \
             **{q:.2} m³/s** · **{:.1} % of `lake_min_area` ({lake_min_cells} cells)**",
            *size as f32 * cell_km2,
            100.0 * *size as f64 / lake_min_cells as f64
        );
        show("its floor", *floor);
    }

    // ── B4 · the drop rule, tested on the orphans ───────────────────────────
    eprintln!("\n── B4 · Finding 37 TASK 2's drop rule, `let is_dry = inflow <= 0.0;` ──");
    let dry = orphans.iter().filter(|o| o.3 <= 0.0).count();
    eprintln!(
        "   untagged components whose LOCAL inflow is exactly 0 ⇒ dropped as a \"DRY salt flat\": \
         **{dry} of {}** ({:.1} %)",
        orphans.len(),
        100.0 * dry as f64 / orphans.len().max(1) as f64
    );
    let dry_but_fed: Vec<&(u32, usize, usize, f32)> = orphans
        .iter()
        .filter(|o| {
            o.3 <= 0.0
                && bs.spillways.iter().any(|sw| {
                    let &(lx, ly) = sw.points.last().unwrap();
                    comp_id[ly as usize * w + lx as usize] == o.0 && sw.discharge_m3s > 0.0
                })
        })
        .collect();
    eprintln!(
        "   **of those, {} are DRY BY THE RULE AND FED BY A SPILLWAY** — the rule reads the LOCAL \
         inflow, the water arrives CHAINED, so the pocket is dropped while a river visibly \
         terminates on it",
        dry_but_fed.len()
    );
    for o in &dry_but_fed {
        let q: f32 = bs
            .spillways
            .iter()
            .filter(|sw| {
                let &(lx, ly) = sw.points.last().unwrap();
                comp_id[ly as usize * w + lx as usize] == o.0
            })
            .map(|sw| sw.discharge_m3s)
            .sum();
        eprintln!(
            "      component {} · {} cells · floor ({},{}) at {:.3} m · local inflow {:.2} · \
             **spillway inflow {q:.1} m³/s**",
            o.0,
            o.1,
            o.2 % w,
            o.2 / w,
            m(&field, o.2),
            o.3
        );
    }

    // ── B5 · what the two eliminated candidates leave: the region, and the flood ──
    eprintln!("\n── B5 · the orphan against the bodies around it ──");
    for (id, size, floor, q) in orphans.iter().take(3) {
        let cells: Vec<usize> = (0..n).filter(|&k| comp_id[k] == *id).collect();
        let hi = cells.iter().map(|&k| m(&field, k)).fold(f32::MIN, f32::max);
        let lo = cells.iter().map(|&k| m(&field, k)).fold(f32::MAX, f32::min);
        // does the component TOUCH a tagged body (8-conn)? then production saw one region and
        // its flood did not claim these cells.
        let mut touching: HashMap<u32, usize> = HashMap::new();
        for &k in &cells {
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for dd in 0..8 {
                let (nx, ny) = (x + D8_DX[dd], y + D8_DY[dd]);
                if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                    continue;
                }
                let nk = ny as usize * w + nx as usize;
                if bs.lake_map[nk] != 0 {
                    *touching.entry(bs.lake_map[nk]).or_default() += 1;
                }
            }
        }
        let mut tv: Vec<(&u32, &usize)> = touching.iter().collect();
        tv.sort_by(|a, b| b.1.cmp(a.1));
        eprintln!(
            "   component {id}: {size} cells · heights {lo:.2} → {hi:.2} m · local inflow \
             {q:.2} m³/s · **touches {} tagged bodies** {:?}",
            touching.len(),
            tv.iter().take(3).map(|(i, c)| (**i, **c)).collect::<Vec<(u32, usize)>>()
        );
        for (bid, _) in tv.iter().take(2) {
            if let Some(b) = bs.lakes.iter().find(|l| l.base.id == **bid) {
                eprintln!(
                    "      the adjacent body {} sits at level **{:.2} m**, floor {:.2} m — the \
                     orphan's cells are {} that level, so a flood capped at \
                     `surface = level.max(sea)` {} reach them",
                    b.base.id,
                    b.level_m,
                    b.level_m - b.depth_m,
                    if lo > b.level_m { "ABOVE" } else { "below" },
                    if lo > b.level_m { "CANNOT" } else { "should" }
                );
            }
        }
        // how much of the component is a BREACH-dug hole?
        let dug = cells.iter().filter(|&&k| raw.data[k] > SEA && field.data[k] <= SEA).count();
        eprintln!(
            "      **{dug} of {size} cells ({:.1} %) are LAND in the pipeline output and SEA only \
             after `breach_monotone`** — Finding 90-B1's population",
            100.0 * dug as f64 / *size as f64
        );
    }

    // ── B5b · the same count with the relevel OFF, one variable ─────────────
    {
        let mut off = on.clone();
        off.merged_union_relevel = None;
        let bs_off = below_sea_basin_lakes_infil(
            &field,
            &dclim,
            &off,
            &ss,
            DOMAIN_KM,
            Some(&detected),
            None,
        );
        let untagged_off = comps
            .iter()
            .filter(|(_, cells)| {
                !cells.iter().any(|&k| bs_off.lake_map[k] != 0 || detected[k] != 0)
            })
            .count();
        let bs_nodet = below_sea_basin_lakes_infil(&field, &dclim, &on, &ss, DOMAIN_KM, None, None);
        let untagged_nodet = comps
            .iter()
            .filter(|(_, cells)| !cells.iter().any(|&k| bs_nodet.lake_map[k] != 0))
            .count();
        eprintln!(
            "\n── B5b · untagged components under one changed variable ──\n   shipped \
             (relevel ON, detected map Some): **{}** · relevel OFF: **{untagged_off}** · \
             detected map None: **{untagged_nodet}** · bodies: {} / {} / {}",
            orphans.len(),
            bs.lakes.len(),
            bs_off.lakes.len(),
            bs_nodet.lakes.len()
        );
    }

    // ── B1 · the area threshold, on both populations ────────────────────────
    let big_orphans = orphans.iter().filter(|o| o.1 >= lake_min_cells).count();
    eprintln!(
        "\n── B1 · the area threshold ──\n   untagged components at or above \
         `lake_min_area` ({lake_min_cells} cells): **{big_orphans} of {}** · the largest untagged \
         is {} cells = **{:.1} %** of the threshold",
        orphans.len(),
        orphans.first().map(|o| o.1).unwrap_or(0),
        100.0 * orphans.first().map(|o| o.1).unwrap_or(0) as f64 / lake_min_cells as f64
    );
    let mut sizes: HashMap<&str, usize> = HashMap::new();
    sizes.insert("inventoried below-sea bodies", bs.lakes.len());
    sizes.insert("detected surface lakes", {
        let ids: HashSet<u32> = detected.iter().copied().filter(|&i| i != 0).collect();
        ids.len()
    });
    eprintln!("   for scale: {sizes:?}");
    eprintln!("\n==========  end Finding 91  ==========\n");
}
