//! ADR Finding 82 voie 2 — B (the hole re-read under the bound), C1 (the −5.5 % budget,
//! decomposed), C3 (river widths), C4 (the 6 764 coastal residual cells).
//!
//! Run: cargo test -p ymir-core --release --test base_level_ledger -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, SEA, build_field, pct, sorted};
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::{PrecipParams, precip_mm_per_year};
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, SegmentKind, SegmentRow, apply_lake_water_balance,
    below_sea_basin_lakes_infil, c1_drainage_windowed, clip_rivers_to_lakes,
    potential_evaporation_mm, runoff_km2_to_m3s,
};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{FlowConfig, RiverSegment, breach_monotone, compute_flow};

const EPS: f32 = 0.5;
const DOMAIN_KM: f32 = 400.0;

/// Chebyshev distance to the coast in cells, on a land/sea mask.
fn dist_to_coast(land: &[bool], w: usize) -> Vec<u16> {
    let n = w * w;
    let mut d = vec![u16::MAX; n];
    let mut q: Vec<u32> = Vec::new();
    for k in 0..n {
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        let mut edge = false;
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let (nx, ny) = (x + dx, y + dy);
                if nx >= 0
                    && ny >= 0
                    && (nx as usize) < w
                    && (ny as usize) < w
                    && land[ny as usize * w + nx as usize] != land[k]
                {
                    edge = true;
                }
            }
        }
        if edge {
            d[k] = 0;
            q.push(k as u32);
        }
    }
    let mut i = 0usize;
    while i < q.len() {
        let k = q[i] as usize;
        i += 1;
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let (nx, ny) = (x + dx, y + dy);
                if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= w {
                    continue;
                }
                let nk = ny as usize * w + nx as usize;
                if d[nk] == u16::MAX {
                    d[nk] = d[k] + 1;
                    q.push(nk as u32);
                }
            }
        }
    }
    d
}

#[test]
#[ignore]
fn base_level_ledger() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 82 voie 2 · the ledger under the bound  ==========");

    let del = build_field(Knobs::shipped());
    let bnd = build_field(Knobs { base_level_m: Some(EPS), ..Knobs::shipped() });
    let bnd_pre =
        build_field(Knobs { base_level_m: Some(EPS), bathymetry_off: true, ..Knobs::shipped() });
    let (w, n) = (del.width, del.data.len());

    // ── C1 · the budget, decomposed by distance to coast and by altitude ──────
    // The Finding 81 attribution ("orography") was ASSERTED. This measures it.
    eprintln!("\n── C1 · Δ(P−PE) delivered → bounded, by band ──");
    for (bed, lat, span) in [("humid", 45.0f32, 40.0f32), ("arid-hot", 25.0, 10.0)] {
        let mut per_dist: Vec<[f64; 2]> = vec![[0.0; 2]; 4];
        let mut per_alt: Vec<[f64; 2]> = vec![[0.0; 2]; 4];
        let mut tot = [0.0f64; 2];
        let mut land_n = [0usize; 2];
        for (vi, f) in [&del, &bnd].into_iter().enumerate() {
            let climate = c1_climate_placed(f, &ss, lat, span, &PrecipParams::default(), DOMAIN_KM);
            let land: Vec<bool> = f.data.iter().map(|&v| v > SEA).collect();
            let dist = dist_to_coast(&land, w);
            for k in 0..n {
                if !land[k] {
                    continue;
                }
                land_n[vi] += 1;
                let p = precip_mm_per_year(climate.precipitation.data[k]);
                let pe = potential_evaporation_mm(climate.temperature.data[k]);
                let s = (p - pe).max(0.0) as f64 * cell_km2 as f64;
                tot[vi] += s;
                // distance bands in km: 0-1, 1-5, 5-20, >20
                let dkm = dist[k] as f32 * CELL_KM;
                let di = if dkm < 1.0 {
                    0
                } else if dkm < 5.0 {
                    1
                } else if dkm < 20.0 {
                    2
                } else {
                    3
                };
                per_dist[di][vi] += s;
                let alt = (f.data[k] - SEA) * n2m;
                let ai = if alt < 10.0 {
                    0
                } else if alt < 100.0 {
                    1
                } else if alt < 500.0 {
                    2
                } else {
                    3
                };
                per_alt[ai][vi] += s;
            }
        }
        let m3s = |v: f64| runoff_km2_to_m3s(v as f32) as f64;
        eprintln!(
            "   {bed}: budget {:.1} → {:.1} m³/s ({:+.1} %) | land cells {} → {}",
            m3s(tot[0]),
            m3s(tot[1]),
            100.0 * (tot[1] - tot[0]) / tot[0],
            land_n[0],
            land_n[1]
        );
        eprintln!(
            "      {:<22} {:>12} {:>12} {:>12} {:>10}",
            "band", "delivered", "bounded", "delta m³/s", "share of Δ"
        );
        let dtot = m3s(tot[1]) - m3s(tot[0]);
        for (nm, v) in [
            ("dist 0-1 km", &per_dist[0]),
            ("dist 1-5 km", &per_dist[1]),
            ("dist 5-20 km", &per_dist[2]),
            ("dist > 20 km", &per_dist[3]),
            ("alt < 10 m", &per_alt[0]),
            ("alt 10-100 m", &per_alt[1]),
            ("alt 100-500 m", &per_alt[2]),
            ("alt > 500 m", &per_alt[3]),
        ] {
            let d = m3s(v[1]) - m3s(v[0]);
            eprintln!(
                "      {nm:<22} {:>12.2} {:>12.2} {:>12.2} {:>9.1} %",
                m3s(v[0]),
                m3s(v[1]),
                d,
                100.0 * d / dtot
            );
        }
    }

    // ── C4 · the coastal residual, on the PRE-CLAMP bounded field ─────────────
    let pre = build_field(Knobs::no_incision());
    let wc_b = water_class(&bnd, SEA);
    let resid: Vec<usize> = (0..n).filter(|&k| pre.data[k] > SEA && bnd.data[k] <= SEA).collect();
    let coastal: Vec<usize> = resid
        .iter()
        .copied()
        .filter(|&k| {
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            (-3i32..=3).any(|dy| {
                (-3i32..=3).any(|dx| {
                    let (nx, ny) = (x + dx, y + dy);
                    nx >= 0
                        && ny >= 0
                        && (nx as usize) < w
                        && (ny as usize) < w
                        && wc_b[ny as usize * w + nx as usize] == 1
                })
            })
        })
        .collect();
    let depths =
        sorted(coastal.iter().map(|&k| -(bnd_pre.data[k] - SEA) * n2m).collect::<Vec<f32>>());
    eprintln!(
        "\n── C4 · the {} coastal residual cells, depth on the PRE-CLAMP bounded field ──\n   p10 \
         {:.4} MEDIAN **{:.4}** p90 {:.4} p99 {:.3} max {:.2} m below sea | shallower than 0.1 m: \
         {:.1} % · than 0.5 m: {:.1} %",
        coastal.len(),
        pct(&depths, 0.10),
        pct(&depths, 0.50),
        pct(&depths, 0.90),
        pct(&depths, 0.99),
        depths.last().copied().unwrap_or(f32::NAN),
        100.0 * depths.iter().filter(|&&d| d <= 0.1).count() as f64 / depths.len() as f64,
        100.0 * depths.iter().filter(|&&d| d <= 0.5).count() as f64 / depths.len() as f64
    );
    // concentration: how many 64-cell coastal tiles hold what share of them
    let mut tile: std::collections::HashMap<(usize, usize), usize> = Default::default();
    for &k in &coastal {
        *tile.entry((k % w / 64, k / w / 64)).or_default() += 1;
    }
    let mut counts: Vec<usize> = tile.values().copied().collect();
    counts.sort_unstable_by(|a, b| b.cmp(a));
    let top10: usize = counts.iter().take(counts.len().div_ceil(10)).sum();
    eprintln!(
        "   spread over **{} tiles of 64×64 cells**; the densest 10 % of tiles hold **{:.1} %** \
         of them (uniform would be 10 %)",
        counts.len(),
        100.0 * top10 as f64 / coastal.len() as f64
    );

    // ── B + C3 · the hole and the widths, humid, under the bound ─────────────
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    eprintln!("\n── B · the hole's dominant term, and C3 · the widths (humid, ratio 7.5) ──");
    for (vnm, f) in [("DELIVERED", &del), ("BOUNDED", &bnd)] {
        let pre0 = c1_drainage_windowed(f, None, &dcfg, &ss, DOMAIN_KM);
        let field = breach_monotone(f, &pre0.flow.filled, &pre0.lake_map, SEA, w, w);
        let flow = compute_flow(&field, &FlowConfig { sea_level: SEA, ..Default::default() });
        let wc = water_class(&field, SEA);
        let climate =
            c1_climate_placed(&field, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dclim = DrainageClimate {
            precip_internal: &climate.precipitation,
            temperature: &climate.temperature,
        };
        let pre2 = c1_drainage_windowed(f, None, &dcfg, &ss, DOMAIN_KM);
        let mut dr = c1_drainage_windowed(&field, Some(&dclim), &dcfg, &ss, DOMAIN_KM);
        dr.lakes = pre2.lakes;
        dr.lake_map = pre2.lake_map;
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
            w,
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
        for k in 0..n {
            if bs.lake_map[k] != 0 {
                dr.lake_map[k] = bs.lake_map[k];
            }
        }
        dr.lakes.extend(bs.lakes.iter().cloned());
        clip_rivers_to_lakes(&mut dr);
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
                source_lake: Some(sw.lake_id),
            });
        }
        // B — the dominant terms, with their kind and their receiver's size
        let mut byq: Vec<&ymir_core::tectonics_c1::drainage::Spillway> =
            bs.spillways.iter().collect();
        byq.sort_by(|a, b| b.discharge_m3s.partial_cmp(&a.discharge_m3s).unwrap());
        let size_of = |id: u32| (0..n).filter(|&k| bs.lake_map[k] == id).count();
        eprintln!("   {vnm}: the 4 largest spillways");
        for sw in byq.iter().take(4) {
            let &(lx, ly) = sw.points.last().unwrap();
            let k = ly as usize * w + lx as usize;
            let rcv = bs.lake_map[k];
            eprintln!(
                "      {:>9.2} m³/s [Spillway] from basin {} ({} cells) → id {rcv} ({} cells), \
                 wc {} at the terminus",
                sw.discharge_m3s,
                sw.lake_id,
                size_of(sw.lake_id),
                if rcv == 0 { 0 } else { size_of(rcv) },
                wc[k]
            );
        }
        // C3 — the widths, in render cells at ratio 7.5 (discharge ×56.25 ⇒ width ×7.5)
        let mut wcell: Vec<f32> = (0..dr.rivers.segments.len())
            .filter(|&i| dr.segment_kind[i] == SegmentKind::Watercourse)
            .map(|i| dr.segment_width_m[i] * 7.5 / 48.83)
            .collect();
        let maxq = (0..dr.rivers.segments.len())
            .filter(|&i| dr.segment_kind[i] == SegmentKind::Watercourse)
            .map(|i| dr.segment_discharge_m3s[i])
            .fold(0.0f32, f32::max);
        wcell.sort_by(|a, b| a.partial_cmp(b).unwrap());
        eprintln!(
            "   {vnm}: max Watercourse **{maxq:.3} m³/s** (real) | width in render cells @7.5: \
             p50 {:.3} p90 {:.3} p99 {:.3} max **{:.2}** | share over 1 cell **{:.2} %**",
            pct(&wcell, 0.50),
            pct(&wcell, 0.90),
            pct(&wcell, 0.99),
            wcell.last().copied().unwrap_or(0.0),
            100.0 * wcell.iter().filter(|&&x| x > 1.0).count() as f64 / wcell.len() as f64
        );
    }
}
