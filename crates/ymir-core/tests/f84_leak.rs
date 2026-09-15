//! ADR 0001 Finding 84 block A2/A3 — the Finding 74 leak, measured at 8192² with the
//! `merged_basin_outlet` gate OFF and ON, two climate beds, base-level bound ON (production).
//!
//! The bench calls `below_sea_basin_lakes_infil` directly rather than
//! `hd_assembly::assemble_hd_drainage`, and that is a decision, not laziness: the production carry
//! returns a `HdDrainageBundle`, which is a CACHED artifact with a `RawCodec`. Threading the
//! termination report through it would change the on-disk format for every key, stale or not. The
//! chain below reproduces `assemble_hd_drainage`'s order (it is the same one Findings 81–83 used,
//! and its figures reproduce Finding 81's bounded column to the digit).
//!
//! Run: cargo test -p ymir-core --release --test f84_leak -- --ignored --nocapture

mod common;

use common::{A_C_CELLS, CELL_KM, Knobs, SEA, build_field, pct, sorted};
use std::collections::{HashMap, HashSet};
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::{PrecipParams, precip_mm_per_year};
use ymir_core::export::height::metric_height_u16;
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, MergedBasinOutlet, ReciprocalRule, SegmentKind, SegmentRow,
    Spillway, apply_lake_water_balance, below_sea_basin_lakes_infil, c1_drainage_windowed,
    clip_rivers_to_lakes, exorheic_lakes_missing_outlet, potential_evaporation_mm,
    runoff_km2_to_m3s,
};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{
    D8_DX, D8_DY, DIR_NONE, FlowConfig, RiverSegment, breach_monotone, compute_flow,
};

const DOMAIN_KM: f32 = 400.0;
/// The geographic scale ratio the author reads widths at (Finding 24): discharge ×ratio²,
/// width ×ratio. A render cell is 400 km / 8192 = 48.83 m.
const RATIO: f32 = 7.5;

fn fnv(g: &GridF32) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for &v in &g.data {
        h ^= v.to_bits() as u64;
        h = h.wrapping_mul(0x100_0000_01b3);
    }
    h
}

#[test]
#[ignore]
fn f84_leak() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 84 · A2/A3 — the leak, gate OFF vs ON  ==========");

    let mut dcfg_off = C1DrainageConfig::default();
    dcfg_off.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    dcfg_off.thresholds.full_tree = false;
    let mut dcfg_on = dcfg_off.clone();
    dcfg_on.merged_basin_outlet = Some(MergedBasinOutlet { rule: ReciprocalRule::LevelThenMerge });

    // ── the field, ONCE. The gate lives in `C1DrainageConfig`, which the upscale never sees ──
    let raw = build_field(Knobs::shipped());
    let (w, h) = (raw.width, raw.height);
    let n = w * h;
    let hl = metric_height_u16(&raw, &ss);
    let u16_step = (hl.max_m - hl.min_m) / 65535.0;
    eprintln!(
        "\n── the field and the climate, which the gate cannot reach ──\n   eroded FNV-1a \
         {:#018x} | u16 step {u16_step:.4} m\n   `merged_basin_outlet` is a field of \
         `C1DrainageConfig`; `FbmUpscaleConfig` (which builds the field) has no path to it, so \
         the height field is the SAME OBJECT for both gate states — not two objects that hash \
         alike.",
        fnv(&raw)
    );
    let pre0 = c1_drainage_windowed(&raw, None, &dcfg_off, &ss, DOMAIN_KM);
    let field = breach_monotone(&raw, &pre0.flow.filled, &pre0.lake_map, SEA, w, h);
    let flow = compute_flow(&field, &FlowConfig { sea_level: SEA, ..Default::default() });
    let wc = water_class(&field, SEA);
    eprintln!(
        "   breached FNV-1a {:#018x} (`c1_drainage_windowed` does not call the below-sea stage, so the breach is gate-independent)",
        fnv(&field)
    );

    // ── the control block that is field-derived, so identical for both gates ──
    let land: Vec<bool> = (0..n).map(|k| field.data[k] > SEA).collect();
    let n_land = land.iter().filter(|&&l| l).count();
    let hh = sorted((0..n).filter(|&k| land[k]).map(|k| (field.data[k] - SEA) * n2m).collect());
    let chan = (0..n).filter(|&k| land[k] && flow.accumulation.data[k] >= A_C_CELLS).count();
    let pits = (0..n).filter(|&k| flow.filled.data[k] > field.data[k] + 1e-7).count();
    eprintln!(
        "   CONTROL (field-derived, identical by construction): land {n_land} | hypso mean \
         {:.2} m p50 {:.2} | channel {:.3} % | pits {pits} | ocean cells {}",
        hh.iter().map(|&x| x as f64).sum::<f64>() / hh.len() as f64,
        pct(&hh, 0.50),
        100.0 * chan as f64 / n_land as f64,
        wc.iter().filter(|&&c| c == 1).count()
    );

    for (bed, lat, span) in [("humid", 45.0f32, 40.0f32), ("arid-hot", 25.0, 10.0)] {
        let climate =
            c1_climate_placed(&field, &ss, lat, span, &PrecipParams::default(), DOMAIN_KM);
        eprintln!(
            "\n╔════════ {bed} ════════╗\n   climate FNV-1a: precip {:#018x} | temp {:#018x} \
             (one climate per bed, not one per gate)",
            fnv(&climate.precipitation),
            fnv(&climate.temperature)
        );
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

        for (gate, dcfg) in [("GATE OFF (shipped)", &dcfg_off), ("**GATE ON**", &dcfg_on)] {
            eprintln!("\n   ── {gate} / {bed} ──");
            let pre = c1_drainage_windowed(&raw, None, dcfg, &ss, DOMAIN_KM);
            let mut dr = c1_drainage_windowed(&field, Some(&dclim), dcfg, &ss, DOMAIN_KM);
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
                dcfg,
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

            // ── the Finding 74 terminal table, each column with its population ──
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
                "      F74: budget {budget:.1} | W->sea {wc_n}/{wc_q:.1} | S->ocean \
                 {n_ocean}/{sp_ocean:.1} | S->chained(wc2) {n_chain}/{sp_chain:.1} | TERMINAL \
                 {terminal:.1} (x{:.3}) | HOLE {:.1}",
                terminal / budget as f64,
                budget as f64 - terminal
            );

            // ── the termination report: Finding 74 block 4's states, and the ones it missed ──
            let t = &bs.termination;
            eprintln!(
                "      TERMINATION: ocean {} | other body {} | detected lake {} ({:.2} m³/s, \
                 surplus DROPPED) | **OWN BODY {} ({:.2} m³/s)** | nothing {}",
                t.to_ocean,
                t.to_other_body,
                t.to_detected_lake,
                t.q_detected_lake_m3s,
                t.to_own_body,
                t.q_own_body_m3s,
                t.to_nothing
            );
            eprintln!(
                "      A3 MERGES: pairs by level {} | pairs merged {} | level gaps {} | unions \
                 retraced {} | unions ENDORHEIC {} ({:.2} m³/s unterminated)",
                t.pairs_by_level,
                t.pairs_merged,
                t.merge_level_gaps_m
                    .iter()
                    .map(|g| format!("{g:.4} m = {:.2} u16 steps", g / u16_step))
                    .collect::<Vec<_>>()
                    .join(", "),
                t.unions_retraced,
                t.unions_endorheic,
                t.q_unions_endorheic_m3s
            );

            // ── the basins, and the two structural counts the round asked for ──
            let size_of = |id: u32| (0..n).filter(|&k| bs.lake_map[k] == id).count();
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
            let mut per: HashMap<u32, HashSet<u32>> = HashMap::new();
            for k in 0..n {
                if wc[k] == 2 && bs.lake_map[k] != 0 {
                    per.entry(bs.lake_map[k]).or_default().insert(comp[k]);
                }
            }
            let multi = per.values().filter(|c| c.len() > 1).count();
            let selfspill = bs
                .spillways
                .iter()
                .filter(|sw| {
                    let &(lx, ly) = sw.points.last().unwrap();
                    bs.lake_map[ly as usize * w + lx as usize] == sw.lake_id
                })
                .count();
            let ar = sorted(bs.basins.iter().map(|b| b.area_km2).collect::<Vec<f32>>());
            let named = dr
                .segment_source_lake
                .iter()
                .enumerate()
                .filter(|(i, s)| dr.segment_kind[*i] == SegmentKind::Spillway && s.is_some())
                .count();
            let inv = dr.lakes.iter().filter(|l| l.base.id >= 1_000_001).count();
            eprintln!(
                "      BASINS {}: area p50 {:.3} SUM {:.1} km² | evap SUM {:.2} | ids over >1 \
                 body **{multi}** | spillways into their OWN id **{selfspill}** | inventoried \
                 below-sea {inv} <-> spillways naming a source {named} (unnamed {})",
                bs.basins.len(),
                pct(&ar, 0.50),
                ar.iter().map(|&x| x as f64).sum::<f64>(),
                bs.basins.iter().map(|b| b.evaporation_m3s as f64).sum::<f64>(),
                bs.spillways.len() - named
            );

            // ── the named terms: the six largest spillways ──
            let mut byq: Vec<&Spillway> = bs.spillways.iter().collect();
            byq.sort_by(|a, b| {
                b.discharge_m3s.partial_cmp(&a.discharge_m3s).unwrap_or(std::cmp::Ordering::Equal)
            });
            eprintln!(
                "      {:>10} {:>10} {:>11} {:>12} {:>11} {:>4} {:>10}",
                "Q m³/s", "from id", "from cells", "-> id", "-> cells", "wc", "path len"
            );
            for sw in byq.iter().take(6) {
                let &(lx, ly) = sw.points.last().unwrap();
                let k = ly as usize * w + lx as usize;
                let rcv = bs.lake_map[k];
                eprintln!(
                    "      {:>10.2} {:>10} {:>11} {:>12} {:>11} {:>4} {:>10}",
                    sw.discharge_m3s,
                    sw.lake_id,
                    size_of(sw.lake_id),
                    if rcv == 0 { "ocean/none".to_string() } else { rcv.to_string() },
                    if rcv == 0 { 0 } else { size_of(rcv) },
                    wc[k],
                    sw.points.len()
                );
            }

            // ── the main river, with everything the round demands beside it ──
            let idx: Vec<usize> = (0..dr.rivers.segments.len())
                .filter(|&i| dr.segment_kind[i] == SegmentKind::Watercourse)
                .collect();
            let (mut best, mut bq) = (usize::MAX, 0.0f32);
            for &i in &idx {
                if dr.segment_discharge_m3s[i] > bq {
                    bq = dr.segment_discharge_m3s[i];
                    best = i;
                }
            }
            let mut wcell: Vec<f32> =
                idx.iter().map(|&i| dr.segment_width_m[i] * RATIO / 48.83).collect();
            wcell.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            if best != usize::MAX {
                eprintln!(
                    "      MAIN RIVER: max Watercourse **{bq:.3} m³/s** real | kind {:?} | \
                     catchment **{:.1} km²** real | width {:.1} m -> **{:.2} render cells** @ratio \
                     {RATIO} | carries **{:.2} %** of the budget",
                    dr.segment_kind[best],
                    dr.segment_drainage_km2[best],
                    dr.segment_width_m[best],
                    dr.segment_width_m[best] * RATIO / 48.83,
                    100.0 * bq as f64 / budget as f64
                );
            }
            eprintln!(
                "      F69 WIDTHS (render cells @{RATIO}): p50 {:.3} p90 {:.3} p99 {:.3} max \
                 **{:.2}** | share over 1 cell **{:.2} %** | {} watercourse runs",
                pct(&wcell, 0.50),
                pct(&wcell, 0.90),
                pct(&wcell, 0.99),
                wcell.last().copied().unwrap_or(0.0),
                100.0 * wcell.iter().filter(|&&x| x > 1.0).count() as f64 / wcell.len() as f64,
                wcell.len()
            );

            // ── the lake invariants ──
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
                "      INVARIANTS: dup {dup} | empty {empty} | footprint above level {above} | \
                 area mismatch {amis} | dangling {dangling} | exorheic WITHOUT outlet **{}**",
                exorheic_lakes_missing_outlet(&dr).len()
            );
        }
    }
}
