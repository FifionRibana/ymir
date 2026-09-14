//! ADR Finding 81 — the columns Finding 80 did not read. No production change, no tuning.
//!
//! Two tests, split by what actually depends on the climate:
//!   * `base_level_terrain` — A3's proof, A4's cross-section, B3's population, C's noise bar and
//!     the naming of the two missing kilometre excursions. **All climate-independent**, because
//!     `build_field` has no climate input: the eroded field comes out of
//!     `upscale_from_c1_with_progress` and the climate enters only at `c1_climate_placed`, in the
//!     drainage stage. Printing these twice "per bed" would print the same numbers twice.
//!   * `base_level_hydrology` — A1, A2, B1, B2, which genuinely differ per bed.
//!
//!
//! ⚠️ REPOINTED at ADR Finding 83: the base-level bound SHIPS, so `Knobs::shipped()` is now
//! the BOUNDED field and the pre-Finding-83 "delivered" world is `Knobs::pre83()`. The
//! columns below mean what they did; only the knob that produces them moved.
//!
//! Run: cargo test -p ymir-core --release --test base_level_read -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, SEA, build_field, pct, sorted};
use ymir_core::export::height::metric_height_u16;
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::coast_metrics::{NECK_KM, coast_shape_thresholds, coast_spurs};
use ymir_core::terrain::contour::marching_squares;

#[allow(dead_code)] // ADR Finding 83: shipped, so no longer passed explicitly
const EPS: f32 = 0.5;
const TRANSECTS: usize = 200;
const REACH: i32 = 20;

fn fnv(g: &GridF32) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for &v in &g.data {
        h ^= v.to_bits() as u64;
        h = h.wrapping_mul(0x100_0000_01b3);
    }
    h
}

fn land_u16(f: &GridF32, ss: &SteinSteinParams) -> Vec<bool> {
    let hl = metric_height_u16(f, ss);
    hl.codes
        .iter()
        .map(|&c| hl.min_m + (c as f32 / 65535.0) * (hl.max_m - hl.min_m) > 0.0)
        .collect()
}

fn to_mask(land: &[bool], w: usize) -> GridF32 {
    let mut m = GridF32::new(w, w, 0.0);
    for k in 0..w * w {
        if land[k] {
            m.data[k] = 1.0;
        }
    }
    m
}

/// Median height profile across the shoreline, the Finding 78 instrument.
/// ⚠️ The normal is oriented toward LOWER ground, so **positive offsets are SEAWARD** — the
/// Finding 78 print had that legend inverted and it is corrected here at the source.
fn profile(f: &GridF32, polys: &[Vec<(f32, f32)>], n2m: f32) -> Vec<f32> {
    let (w, h) = (f.width, f.height);
    let mut cols: Vec<Vec<f32>> = vec![Vec::new(); (2 * REACH + 1) as usize];
    let mut idx: Vec<usize> = (0..polys.len()).filter(|&i| polys[i].len() >= 32).collect();
    idx.sort_by_key(|&i| std::cmp::Reverse(polys[i].len()));
    let mut taken = 0usize;
    'outer: for &i in &idx {
        let pl = &polys[i];
        let step = (pl.len() / 16).max(8);
        let mut j = 4usize;
        while j + 4 < pl.len() {
            let (ax, ay) = pl[j - 4];
            let (bx, by) = pl[j + 4];
            let (tx, ty) = (bx - ax, by - ay);
            let len = (tx * tx + ty * ty).sqrt().max(1e-6);
            let (mut nx, mut ny) = (-ty / len, tx / len);
            let (px, py) = pl[j];
            let probe = |dx: f32, dy: f32| -> f32 {
                let (sx, sy) = ((px + dx).round() as i64, (py + dy).round() as i64);
                if sx < 0 || sy < 0 || sx as usize >= w || sy as usize >= h {
                    return f32::NAN;
                }
                f.data[sy as usize * w + sx as usize]
            };
            let (a, b) = (probe(nx * 3.0, ny * 3.0), probe(-nx * 3.0, -ny * 3.0));
            if a.is_nan() || b.is_nan() {
                j += step;
                continue;
            }
            if a > b {
                nx = -nx;
                ny = -ny;
            }
            let mut col = Vec::new();
            let mut ok = true;
            for s in -REACH..=REACH {
                let v = probe(nx * s as f32, ny * s as f32);
                if v.is_nan() {
                    ok = false;
                    break;
                }
                col.push((v - SEA) * n2m);
            }
            if ok {
                for (c, v) in col.into_iter().enumerate() {
                    cols[c].push(v);
                }
                taken += 1;
                if taken >= TRANSECTS {
                    break 'outer;
                }
            }
            j += step;
        }
    }
    assert!(taken >= TRANSECTS / 2, "rule 10: only {taken} transects");
    cols.into_iter()
        .map(|c| {
            let s = sorted(c);
            if s.is_empty() { f32::NAN } else { s[s.len() / 2] }
        })
        .collect()
}

#[test]
#[ignore]
fn base_level_terrain() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    eprintln!("\n==========  Finding 81 · the terrain-side columns  ==========");

    let pre = build_field(Knobs::no_incision());
    let del = build_field(Knobs::pre83());
    let bnd = build_field(Knobs::shipped());
    let (w, n) = (del.width, del.data.len());

    // ── A3 · the climate-independence proof ───────────────────────────────────
    eprintln!(
        "\n── A3 · why there is no second bed here ──\n   `build_field` takes no climate: the \
         eroded field is `upscale_from_c1_with_progress`, and the climate enters at \
         `c1_climate_placed`, in the DRAINAGE stage. So every terrain-side column is identical \
         between beds BY CONSTRUCTION.\n   field FNV-1a: pre {:#018x} | delivered {:#018x} | \
         bounded {:#018x}  — one field per variant, not one per bed.",
        fnv(&pre),
        fnv(&del),
        fnv(&bnd)
    );

    // ── C · the noise bar, measured on the observable it judges ───────────────
    eprintln!("\n── C · the quantisation bar, on the MASK CONTOUR (u16 against f32) ──");
    eprintln!(
        "   {:<16} {:>12} {:>12} {:>10} {:>10}",
        "variant", "f32 km", "u16 km", "delta km", "delta %"
    );
    let mut bars: Vec<f64> = Vec::new();
    for (nm, f) in [("PRE-INCISION", &pre), ("DELIVERED", &del), ("BOUNDED", &bnd)] {
        let l32: Vec<bool> = f.data.iter().map(|&v| v > SEA).collect();
        let l16 = land_u16(f, &ss);
        let c32 = coast_shape_thresholds(
            &marching_squares(&to_mask(&l32, w), 0.5),
            CELL_KM,
            1.0,
            NECK_KM,
        )
        .coast_km;
        let c16 = coast_shape_thresholds(
            &marching_squares(&to_mask(&l16, w), 0.5),
            CELL_KM,
            1.0,
            NECK_KM,
        )
        .coast_km;
        let d = (c16 - c32) as f64;
        bars.push(100.0 * d.abs() / c32 as f64);
        eprintln!("   {nm:<16} {c32:>12.1} {c16:>12.1} {d:>10.2} {:>10.3}", 100.0 * d / c32 as f64);
    }
    let bar = bars.iter().cloned().fold(0.0f64, f64::max);
    eprintln!(
        "   ⇒ the bar, re-declared on this observable: **{bar:.3} %** (the widest of the three)"
    );

    // where the +5 km of Finding 80 come from: the mask XOR against the pre-incision
    let l_pre = land_u16(&pre, &ss);
    let l_bnd = land_u16(&bnd, &ss);
    let l_del = land_u16(&del, &ss);
    let xor_pb = (0..n).filter(|&k| l_pre[k] != l_bnd[k]).count();
    let xor_pd = (0..n).filter(|&k| l_pre[k] != l_del[k]).count();
    let c_pre =
        coast_shape_thresholds(&marching_squares(&to_mask(&l_pre, w), 0.5), CELL_KM, 1.0, NECK_KM);
    let c_bnd =
        coast_shape_thresholds(&marching_squares(&to_mask(&l_bnd, w), 0.5), CELL_KM, 1.0, NECK_KM);
    eprintln!(
        "   the Finding 80 +5 km: contour {:.1} -> {:.1} km ({:+.2} km, {:+.3} %) | mask XOR \
         pre-incision vs BOUNDED **{xor_pb} cells**, vs DELIVERED {xor_pd} cells",
        c_pre.coast_km,
        c_bnd.coast_km,
        c_bnd.coast_km - c_pre.coast_km,
        100.0 * (c_bnd.coast_km - c_pre.coast_km) as f64 / c_pre.coast_km as f64
    );

    // ── C-bis · NAME the two kilometre excursions that go missing ─────────────
    let (sp_pre, _) =
        coast_spurs(&marching_squares(&to_mask(&l_pre, w), 0.5), CELL_KM, 1.0, NECK_KM);
    let (sp_bnd, _) =
        coast_spurs(&marching_squares(&to_mask(&l_bnd, w), 0.5), CELL_KM, 1.0, NECK_KM);
    eprintln!(
        "\n── C-bis · the Δ(>=1 km) = -2, named ──\n   pre-incision {} spurs >= 1 km, bounded {}",
        sp_pre.len(),
        sp_bnd.len()
    );
    let mut missing: Vec<(f32, f32, f32, f32)> = Vec::new(); // x, y, len_pre, nearest len_bnd
    for a in &sp_pre {
        let mut best = (f32::MAX, 0.0f32);
        for b in &sp_bnd {
            let d = (a.mid.0 - b.mid.0).powi(2) + (a.mid.1 - b.mid.1).powi(2);
            if d < best.0 {
                best = (d, b.len_km);
            }
        }
        if best.0.sqrt() > 8.0 {
            missing.push((a.mid.0, a.mid.1, a.len_km, f32::NAN));
        } else if best.1 < 1.0 {
            missing.push((a.mid.0, a.mid.1, a.len_km, best.1));
        }
    }
    for (x, y, lp, lb) in &missing {
        eprintln!(
            "     cell {:.0},{:.0} | length pre-incision {lp:.2} km -> bounded {}",
            x,
            y,
            if lb.is_nan() { "NO spur within 8 cells".to_string() } else { format!("{lb:.2} km") }
        );
    }
    if missing.is_empty() {
        eprintln!("     none matched — the -2 is a re-pairing of the walk, not two named losses");
    }

    // ── B3 · where the 7 738 residual drowned cells are ───────────────────────
    let wc = water_class(&bnd, SEA);
    let resid: Vec<usize> = (0..n).filter(|&k| pre.data[k] > SEA && bnd.data[k] <= SEA).collect();
    let near = |set: &[usize], class: u8, r: i32| -> usize {
        set.iter()
            .filter(|&&k| {
                let (x, y) = ((k % w) as i32, (k / w) as i32);
                (-r..=r).any(|dy| {
                    (-r..=r).any(|dx| {
                        let (nx, ny) = (x + dx, y + dy);
                        nx >= 0
                            && ny >= 0
                            && (nx as usize) < w
                            && (ny as usize) < w
                            && wc[ny as usize * w + nx as usize] == class
                    })
                })
            })
            .count()
    };
    let n_basin = near(&resid, 2, 3);
    let n_ocean = near(&resid, 1, 3);
    eprintln!(
        "\n── B3 · the {} residual drowned cells, by population ──\n   within 3 cells of an \
         ENCLOSED below-sea basin (wc=2): **{n_basin} ({:.1} %)**\n   within 3 cells of the OCEAN \
         (wc=1): **{n_ocean} ({:.1} %)**\n   neither: **{} ({:.1} %)**",
        resid.len(),
        100.0 * n_basin as f64 / resid.len() as f64,
        100.0 * n_ocean as f64 / resid.len() as f64,
        resid.len() - n_basin.max(n_ocean).min(resid.len()),
        100.0 * (resid.len() as f64 - n_basin.max(n_ocean) as f64) / resid.len() as f64
    );

    // ── A4 · the cross-section, and the x6 flattening ─────────────────────────
    eprintln!(
        "\n── A4 · the cross-section, median over {TRANSECTS} transects, metres above sea ──"
    );
    eprintln!(
        "   NEGATIVE offsets are INLAND, positive are SEAWARD (the Finding 78 legend was inverted)"
    );
    eprint!("   {:<16}", "offset");
    for s in (-REACH..=REACH).step_by(2) {
        eprint!("{s:>8}");
    }
    eprintln!();
    let mut inland20: Vec<f32> = Vec::new();
    for (nm, f) in [("PRE-INCISION", &pre), ("DELIVERED", &del), ("BOUNDED", &bnd)] {
        let p = profile(f, &marching_squares(f, SEA), n2m);
        inland20.push(p[0]);
        eprint!("   {nm:<16}");
        for (i, v) in p.iter().enumerate() {
            if (i as i32 - REACH) % 2 == 0 {
                eprint!("{v:>8.1}");
            }
        }
        eprintln!();
    }
    eprintln!(
        "   the Finding 78 x6 coastal flattening, at 20 cells inland: pre {:.1} m -> delivered \
         {:.1} m (x{:.2}) -> bounded {:.1} m (x{:.2}) ⇒ **{:.0} % of the flattening remains**",
        inland20[0],
        inland20[1],
        inland20[0] / inland20[1],
        inland20[2],
        inland20[0] / inland20[2],
        100.0 * (inland20[0] - inland20[2]) as f64 / (inland20[0] - inland20[1]) as f64
    );
}

/// ADR Finding 81 blocks A1, A2, B1, B2 — the hydrology, which is where the bed actually bites.
#[test]
#[ignore]
fn base_level_hydrology() {
    use std::collections::{HashMap, HashSet};
    use ymir_core::climate::c1_climate_placed;
    use ymir_core::climate::precipitation::{PrecipParams, precip_mm_per_year};
    use ymir_core::tectonics_c1::drainage::{
        C1DrainageConfig, DrainageClimate, LakeType, SegmentKind, SegmentRow,
        apply_lake_water_balance, below_sea_basin_lakes_infil, c1_drainage_windowed,
        clip_rivers_to_lakes, exorheic_lakes_missing_outlet, potential_evaporation_mm,
        runoff_km2_to_m3s,
    };
    use ymir_core::terrain::flow::{
        D8_DX, D8_DY, DIR_NONE, FlowConfig, RiverSegment, breach_monotone, compute_flow,
    };

    const DOMAIN_KM: f32 = 400.0;
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 81 · the hydrology, two variants x two beds  ==========");

    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;

    for (vnm, kn) in
        [("DELIVERED (pre-83)", Knobs::pre83()), ("BOUNDED (production)", Knobs::shipped())]
    {
        let raw = build_field(kn);
        let (w, h) = (raw.width, raw.height);
        let n = w * h;
        let pre0 = c1_drainage_windowed(&raw, None, &dcfg, &ss, DOMAIN_KM);
        let field = breach_monotone(&raw, &pre0.flow.filled, &pre0.lake_map, SEA, w, h);
        let flow = compute_flow(&field, &FlowConfig { sea_level: SEA, ..Default::default() });
        let wc = water_class(&field, SEA);

        for (bed, lat, span) in [("humid", 45.0f32, 40.0f32), ("arid-hot", 25.0, 10.0)] {
            eprintln!("\n╔═══ {vnm} / {bed} ═══╗");
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
                (0..n).filter(|&k| field.data[k] > SEA).map(|k| surplus[k] as f64).sum::<f64>()
                    as f32,
            );

            // production tail
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

            // ── B1 · the Finding 74 terminal table ────────────────────────────
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
            let end_wc = |sw: &ymir_core::tectonics_c1::drainage::Spillway| -> u8 {
                let &(lx, ly) = sw.points.last().unwrap();
                wc[ly as usize * w + lx as usize]
            };
            let sp_ocean: f64 = bs
                .spillways
                .iter()
                .filter(|s| end_wc(s) == 1)
                .map(|s| s.discharge_m3s as f64)
                .sum();
            let n_ocean = bs.spillways.iter().filter(|s| end_wc(s) == 1).count();
            let sp_chain: f64 = bs
                .spillways
                .iter()
                .filter(|s| end_wc(s) == 2)
                .map(|s| s.discharge_m3s as f64)
                .sum();
            let n_chain = bs.spillways.iter().filter(|s| end_wc(s) == 2).count();
            let terminal = wc_q + sp_ocean;
            eprintln!(
                "   B1 F74: budget {budget:.1} | W->sea {wc_n}/{wc_q:.1} | S->ocean \
                 {n_ocean}/{sp_ocean:.1} | S->chained(wc2) {n_chain}/{sp_chain:.1} | TERMINAL \
                 {terminal:.1} (x{:.3}) | HOLE {:.1}",
                terminal / budget as f64,
                budget as f64 - terminal
            );

            // ── A1 · the lake invariants ──────────────────────────────────────
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
                if fmax.get(&lk.base.id).copied().unwrap_or(f32::MIN)
                    > SEA + lk.level_m / n2m + 1e-6
                {
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
                "   A1 invariants: dup {dup} | empty {empty} | footprint ABOVE level {above} | \
                 area mismatch {amis} | dangling ids {dangling} | exorheic WITHOUT outlet **{}**",
                exorheic_lakes_missing_outlet(&dr).len()
            );

            // ── A2 · BasinSummary + the F71-A4 reconciliation ─────────────────
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
                "   A2 BasinSummary {}: area_km2 p50 {:.3} SUM {:.1} | inflow p50 {:.3} SUM \
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

            // ── B2 · the eight merged halves ──────────────────────────────────
            if bed == "humid" {
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
                eprintln!(
                    "   B2 merged ids ({} regions total): **{}** span more than one component",
                    next,
                    split.len()
                );
                for (id, v) in &split {
                    eprintln!("      {id} -> components of {v:?} cells");
                }
                for id in [
                    1_000_004u32,
                    1_000_008,
                    1_000_014,
                    1_000_016,
                    1_000_021,
                    1_000_035,
                    1_000_053,
                    1_000_056,
                ] {
                    let e = per.get(&id);
                    eprintln!(
                        "      Finding 75-E id {id}: {}",
                        match e {
                            None => "ABSENT from the below-sea map".to_string(),
                            Some(c) => {
                                let mut v: Vec<usize> = c.values().copied().collect();
                                v.sort_unstable_by(|a, b| b.cmp(a));
                                format!("{} component(s), cells {v:?}", c.len())
                            }
                        }
                    );
                }
            }
        }
    }
}
