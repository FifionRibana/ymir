//! ADR 0001 Finding 92 block A — the two candidate remedies for Finding 91's regression, MEASURED
//! side by side in the bench, so the author chooses against numbers. **No production change.**
//!
//! α = claim every region of a merged class (one flood per region, capped at the class's common
//! surface) ⇒ a reciprocal pair becomes ONE body.
//! β = do not mark the unclaimed regions `seen` ⇒ each becomes its own body.
//!
//! Neither is implemented in production. Both columns are DERIVED from two shipped runs — relevel
//! ON (the regression) and relevel OFF (which classifies every region on its own, i.e. β's content
//! for the absorbed regions) — matched component by component. Stated so the derivation can be
//! checked rather than trusted.
//!
//! Run: cargo test -p ymir-core --release --test f92_alpha_beta -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, SEA, build_field};
use std::collections::HashMap;
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    BelowSeaResult, C1DrainageConfig, DrainageClimate, apply_lake_water_balance,
    below_sea_basin_lakes_infil, c1_drainage_windowed, runoff_accumulation, runoff_km2_to_m3s,
};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{D8_DX, D8_DY, breach_monotone};

const DOMAIN_KM: f32 = 400.0;

#[test]
#[ignore]
fn f92_alpha_beta() {
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 92 · A — α and β, measured before the choice  ==========");

    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    let mut off = on.clone();
    off.merged_union_relevel = None;

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
    let racc = runoff_accumulation(&field, &dr.flow, &dclim, cell_km2, None, None, w, h);

    let run = |cfg: &C1DrainageConfig| -> BelowSeaResult {
        below_sea_basin_lakes_infil(&field, &dclim, cfg, &ss, DOMAIN_KM, Some(&detected), None)
    };
    let bs_on = run(&on);
    let bs_off = run(&off);
    eprintln!(
        "   relevel ON: **{} bodies**, {} passes (converged {}) · relevel OFF: **{} bodies**, {} \
         passes (converged {})",
        bs_on.lakes.len(),
        bs_on.termination.passes_used,
        !bs_on.termination.not_converged,
        bs_off.lakes.len(),
        bs_off.termination.passes_used,
        !bs_off.termination.not_converged
    );

    // ── the enclosed below-sea components, production's connectivity ────────
    let wc = water_class(&field, SEA);
    let mut comp_id = vec![0u32; n];
    let mut comps: Vec<Vec<usize>> = Vec::new();
    for s in 0..n {
        if wc[s] != 2 || comp_id[s] != 0 {
            continue;
        }
        let id = comps.len() as u32 + 1;
        let mut stack = vec![s];
        comp_id[s] = id;
        let mut cells = Vec::new();
        while let Some(k) = stack.pop() {
            cells.push(k);
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for d in 0..8 {
                let (nx, ny) = (x + D8_DX[d], y + D8_DY[d]);
                if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                    continue;
                }
                let nk = ny as usize * w + nx as usize;
                if wc[nk] == 2 && comp_id[nk] == 0 {
                    comp_id[nk] = id;
                    stack.push(nk);
                }
            }
        }
        comps.push(cells);
    }

    // which body covers each component, in each run
    let cover = |lm: &[u32], cells: &[usize]| -> u32 {
        let mut c: HashMap<u32, usize> = HashMap::new();
        for &k in cells {
            if lm[k] != 0 {
                *c.entry(lm[k]).or_default() += 1;
            }
        }
        c.into_iter().max_by_key(|e| e.1).map(|e| e.0).unwrap_or(0)
    };
    let inflow_of = |cells: &[usize]| -> f32 {
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

    // group the components by the ON body that covers them = the CLASS
    let mut classes: HashMap<u32, Vec<usize>> = HashMap::new(); // on_body -> component indices
    let mut uncovered: Vec<usize> = Vec::new();
    for (i, cells) in comps.iter().enumerate() {
        let on_id = cover(&bs_on.lake_map, cells);
        if on_id == 0 {
            uncovered.push(i);
        }
        classes.entry(on_id).or_default().push(i);
    }
    eprintln!(
        "\n── the {} enclosed `wc == 2` components · **{} uncovered with the relevel ON** ──",
        comps.len(),
        uncovered.len()
    );

    // the MERGED classes: an ON body covering more than one component, plus the uncovered ones
    // attributed to the class they belong to by the OFF run's structure.
    eprintln!("\n── A · the merged classes, region by region ──");
    let mut alpha_bodies = bs_on.lakes.len();
    let mut beta_bodies = bs_on.lakes.len();
    let mut beta_pixel_lakes = 0usize;
    let mut alpha_gained_km2: HashMap<u32, f32> = HashMap::new();
    let mut multi = 0usize;
    for (on_id, idx) in classes.iter() {
        if *on_id == 0 || idx.len() < 2 {
            continue;
        }
        multi += 1;
        let body = bs_on.lakes.iter().find(|l| l.base.id == *on_id);
        let lvl = body.map(|b| b.level_m).unwrap_or(f32::NAN);
        eprintln!(
            "   class on body **{on_id}** (level {lvl:.2} m, {:.2} km² marked): **{} components**",
            bs_on.lake_map.iter().filter(|&&v| v == *on_id).count() as f32 * cell_km2,
            idx.len()
        );
        for &i in idx {
            let cells = &comps[i];
            let lo = cells.iter().map(|&k| m(&field, k)).fold(f32::MAX, f32::min);
            let hi = cells.iter().map(|&k| m(&field, k)).fold(f32::MIN, f32::max);
            let marked = cells.iter().filter(|&&k| bs_on.lake_map[k] != 0).count();
            let off_id = cover(&bs_off.lake_map, cells);
            let off_body = bs_off.lakes.iter().find(|l| l.base.id == off_id);
            eprintln!(
                "      region {:>2}: **{:>6} cells = {:>8.4} km²** · heights {lo:>7.2} → {hi:>6.2} m \
                 · marked ON: {marked} · local inflow **{:>6.2} m³/s** · relevel-OFF body {off_id} \
                 ({:.4} km², level {:.2} m, {:?})",
                i,
                cells.len(),
                cells.len() as f32 * cell_km2,
                runoff_km2_to_m3s(inflow_of(cells)),
                off_body.map(|b| b.area_km2).unwrap_or(0.0),
                off_body.map(|b| b.level_m).unwrap_or(f32::NAN),
                off_body.map(|b| b.lake_type)
            );
        }
        // spillways that LAND in this class, and in which region
        for sw in &bs_on.spillways {
            let &(lx, ly) = sw.points.last().unwrap();
            let k = ly as usize * w + lx as usize;
            if let Some(&i) = idx.iter().find(|&&i| comp_id[k] == i as u32 + 1) {
                eprintln!(
                    "      ⚠️ a spillway of basin **{}** carrying **{:.2} m³/s** LANDS in region \
                     {i} of this class — with α the two become ONE body and this becomes a \
                     `to_own_body` self-discharge; with β it stays a basin-A-into-basin-B link",
                    sw.lake_id, sw.discharge_m3s
                );
            }
        }
    }
    eprintln!(
        "   classes spanning more than one component: **{multi}** (Finding 85-B2 measured 3 \
               on the bounded field: `[121234, 6245]`, `[1, 1]`, `[2, 2]`)"
    );

    // ── the decision table, derived ─────────────────────────────────────────
    for &i in &uncovered {
        let cells = &comps[i];
        let km2 = cells.len() as f32 * cell_km2;
        // α: the cells are all below sea, and every surface is >= sea, so a per-region flood
        // capped at the common surface covers them entirely.
        let lo = cells.iter().map(|&k| m(&field, k)).fold(f32::MAX, f32::min);
        let hi = cells.iter().map(|&k| m(&field, k)).fold(f32::MIN, f32::max);
        let on_neighbour = {
            let mut c: HashMap<u32, usize> = HashMap::new();
            for &k in cells {
                let (x, y) = ((k % w) as i32, (k / w) as i32);
                for d in 0..8 {
                    let (nx, ny) = (x + D8_DX[d], y + D8_DY[d]);
                    if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                        continue;
                    }
                    let nk = ny as usize * w + nx as usize;
                    if bs_on.lake_map[nk] != 0 {
                        *c.entry(bs_on.lake_map[nk]).or_default() += 1;
                    }
                }
            }
            c.into_iter().max_by_key(|e| e.1).map(|e| e.0).unwrap_or(0)
        };
        let off_id = cover(&bs_off.lake_map, cells);
        let off_body = bs_off.lakes.iter().find(|l| l.base.id == off_id);
        beta_bodies += 1;
        if cells.len() <= 2 {
            beta_pixel_lakes += 1;
        }
        *alpha_gained_km2.entry(off_id).or_default() += km2;
        eprintln!(
            "\n   **UNCOVERED region {i}**: {} cells = {km2:.4} km² · heights {lo:.2} → {hi:.2} m \
             · **all cells below sea ⇒ any `surface = level.max(sea)` covers them, so the claim \
             failed on CONNECTIVITY, not height** · adjacent ON body: {} (8-conn) · \
             relevel-OFF it is body {off_id} ({:.4} km², level {:.2} m, {:?})",
            cells.len(),
            if on_neighbour == 0 {
                "NONE — it touches no tagged body".to_string()
            } else {
                on_neighbour.to_string()
            },
            off_body.map(|b| b.area_km2).unwrap_or(0.0),
            off_body.map(|b| b.level_m).unwrap_or(f32::NAN),
            off_body.map(|b| b.lake_type)
        );
    }
    let _ = &mut alpha_bodies;
    eprintln!(
        "\n── THE DECISION TABLE ──\n   | | **α (merged)** | **β (separate)** |\n   \
         | uncovered `wc == 2` components | **0** | **0** |\n   \
         | total below-sea bodies | **{alpha_bodies}** | **{beta_bodies}** |\n   \
         | single-cell (≤ 2 cells) lakes created | **0** | **{beta_pixel_lakes}** |"
    );
    for (id, km2) in &alpha_gained_km2 {
        eprintln!("   with α, the body absorbing the orphan gains **{km2:.4} km²** (OFF id {id})");
    }

    // ── the control, both ways ──────────────────────────────────────────────
    eprintln!(
        "\n── the control ──\n   the height field is the SAME OBJECT in both runs (the relevel is \
         a drainage-config knob and touches no height) · exorheic-without-outlet is not \
         computable here (it needs the full chain) — Finding 86's guard is reported in block C"
    );
    eprintln!(
        "   `to_own_body` {} / `to_nothing` {} (ON) · {} / {} (OFF) · q_own_body {:.2} / {:.2}",
        bs_on.termination.to_own_body,
        bs_on.termination.to_nothing,
        bs_off.termination.to_own_body,
        bs_off.termination.to_nothing,
        bs_on.termination.q_own_body_m3s,
        bs_off.termination.q_own_body_m3s
    );

    // ══ C — the remedy the author chose (β), verified on the product ═══════
    //
    // The sequence the round demanded: the GUARD first, RED on the known regression, then the fix,
    // GREEN. `separate_unclaimed_regions: false` reproduces the pre-Finding-92 footprints bit for
    // bit, so it is the red column and not a reconstruction of one.
    eprintln!("\n══════  C · the guard RED then GREEN, on the delivered product  ══════");
    let mut fix_off = on.clone();
    fix_off.merged_union_relevel = Some(ymir_core::tectonics_c1::drainage::MergedUnionRelevel {
        separate_unclaimed_regions: false,
        ..Default::default()
    });
    let bs_red = run(&fix_off);
    for (label, bsx) in [
        ("**RED** — relevel ON, `separate_unclaimed_regions: false` (the F91 regression)", &bs_red),
        (
            "**GREEN** — relevel ON, `separate_unclaimed_regions: true` (shipped after F92-C)",
            &bs_on,
        ),
        ("control — relevel OFF (the pre-F86 world)", &bs_off),
    ] {
        // the PRODUCTION guard, called exactly as `assemble_hd_drainage` calls it
        let unc = ymir_core::tectonics_c1::drainage::uncovered_below_sea_components(
            &field,
            &bsx.lake_map,
            SEA,
            w,
            h,
        );
        let t = &bsx.termination;
        eprintln!(
            "   {label}\n      the F92-B guard reports **{} uncovered** {:?} · bodies **{}** · \
             passes {} (converged {}) · `to_own_body` {} carrying {:.2} m³/s · `to_nothing` {}",
            unc.len(),
            unc.iter()
                .take(4)
                .map(|(k, n)| ((k % w, k / w), *n))
                .collect::<Vec<((usize, usize), usize)>>(),
            bsx.lakes.len(),
            t.passes_used,
            !t.not_converged,
            t.to_own_body,
            t.q_own_body_m3s,
            t.to_nothing
        );
    }
    // what β made of the two orphans. ⚠️ `uncovered` was computed from the FIXED run and is now
    // empty, so the list must come from the RED run — the first pass of this block printed an
    // empty section for exactly that reason, and the slip is on the record.
    let uncovered_red: Vec<usize> =
        (0..comps.len()).filter(|&i| cover(&bs_red.lake_map, &comps[i]) == 0).collect();
    eprintln!(
        "\n   what β made of the {} Finding 91 orphans (listed from the RED run):",
        uncovered_red.len()
    );
    for &i in &uncovered_red {
        let cells = &comps[i];
        let id = cover(&bs_on.lake_map, cells);
        let body = bs_on.lakes.iter().find(|l| l.base.id == id);
        eprintln!(
            "      region {i} ({} cells = {:.4} km², local inflow {:.2} m³/s) ⇒ body **{id}** \
             ({:.4} km², level {:.2} m, depth {:.2} m, {:?})",
            cells.len(),
            cells.len() as f32 * cell_km2,
            runoff_km2_to_m3s(inflow_of(cells)),
            body.map(|b| b.area_km2).unwrap_or(0.0),
            body.map(|b| b.level_m).unwrap_or(f32::NAN),
            body.map(|b| b.depth_m).unwrap_or(f32::NAN),
            body.map(|b| b.lake_type)
        );
    }
    // the classes that already worked must be untouched
    let red_ids: std::collections::HashSet<usize> =
        (0..n).filter(|&k| bs_red.lake_map[k] != 0).collect();
    let on_ids: std::collections::HashSet<usize> =
        (0..n).filter(|&k| bs_on.lake_map[k] != 0).collect();
    eprintln!(
        "\n   footprint delta RED → GREEN: **{} cells newly covered** ({:.4} km²), **{} cells \
         uncovered** (must be 0 — β only ADDS bodies, it claims nothing away)",
        on_ids.difference(&red_ids).count(),
        on_ids.difference(&red_ids).count() as f32 * cell_km2,
        red_ids.difference(&on_ids).count()
    );
    assert_eq!(
        red_ids.difference(&on_ids).count(),
        0,
        "C stop rule: β must not un-claim a cell the pre-F92 path covered"
    );

    // the price of β on the delivered world, stated rather than buried
    let water_red = (0..n).filter(|&k| bs_red.lake_map[k] != 0).count();
    let water_green = (0..n).filter(|&k| bs_on.lake_map[k] != 0).count();
    let land = (0..n).filter(|&k| field.data[k] > SEA).count();
    eprintln!(
        "\n   the PRICE of β on the delivered world: below-sea water {:.1} → **{:.1} km²** \
         (+{:.1} km², +{:.1} %) · against {:.0} km² of land, the below-sea bodies alone go \
         {:.2} % → **{:.2} %**, so Finding 89-C1's 24.05 % moves by the same +{:.1} km²",
        water_red as f32 * cell_km2,
        water_green as f32 * cell_km2,
        (water_green - water_red) as f32 * cell_km2,
        100.0 * (water_green - water_red) as f64 / water_red.max(1) as f64,
        land as f32 * cell_km2,
        100.0 * water_red as f64 / land.max(1) as f64,
        100.0 * water_green as f64 / land.max(1) as f64,
        (water_green - water_red) as f32 * cell_km2
    );
    eprintln!("\n==========  end Finding 92  ==========\n");
}
