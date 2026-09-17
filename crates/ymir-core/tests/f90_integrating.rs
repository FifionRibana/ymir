//! ADR 0001 Finding 90 block C — Finding 44's item 2 at the SHIPPED `k_time`.
//!
//! The experiment: hold `k_time = k · dt · iterations = 9000` and raise the pass count, so the
//! integrated erosion budget is identical to the delivered field's and only the Courant number
//! falls. Finding 44 proved `K` and the duration are not separately observable, so holding
//! `k_time` is the only way to move one dial alone.
//!
//! Precondition (Finding 90-A): the +115 is attributed to `breach_monotone`. Recorded here so the
//! two blocks cannot drift apart — the milestones below therefore report the coast criterion on
//! BOTH the eroded and the breached stage, since the stage is where the criterion lives.
//!
//! Run: cargo test -p ymir-core --release --test f90_integrating -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, SEA, build_field, land_u16, majority, pct, sorted, to_mask};
use std::collections::{HashMap, HashSet};
use std::time::Instant;
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::erosion::stream_power::{RELIEF_V1_A_C_KM2, SHIPPED_K_TIME};
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, LakeType, SegmentKind, SegmentRow, apply_lake_water_balance,
    below_sea_basin_lakes_infil, c1_drainage_windowed, clip_rivers_to_lakes,
    resolve_exorheic_without_outlet, runoff_accumulation, runoff_km2_to_m3s,
};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::coast_metrics::{MIN_SPUR_KM, NECK_KM, coast_spurs};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::flow::{D8_DX, D8_DY, breach_monotone};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;

fn median(v: &[f32]) -> f32 {
    pct(&sorted(v.to_vec()), 0.50)
}

fn slope_deg(f: &GridF32, k: usize, n2m: f32, w: usize, h: usize) -> f32 {
    let (x, y) = ((k % w) as i32, (k / w) as i32);
    let mut best = 0.0f32;
    for d in 0..8 {
        let nx = (x + D8_DX[d]).rem_euclid(w as i32) as usize;
        let ny = (y + D8_DY[d]).rem_euclid(h as i32) as usize;
        let diag = D8_DX[d] != 0 && D8_DY[d] != 0;
        let dist = if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M };
        best = best.max((f.data[k] - f.data[ny * w + nx]).abs() * n2m / dist);
    }
    best.atan().to_degrees()
}

fn spurs_km(mask: &[bool], w: usize, factor: usize) -> usize {
    let (m, ww) = majority(mask, w, factor);
    let ck = CELL_KM * factor as f32;
    coast_spurs(&marching_squares(&to_mask(&m, ww), 0.5), ck, MIN_SPUR_KM, NECK_KM).0.len()
}

#[test]
#[ignore]
fn f90_integrating() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 90 · C — Finding 44 item 2 at the shipped k_time  ==========");

    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;

    // ══ C0 — the configuration, read from the code ══════════════════════════
    {
        use ymir_core::erosion::stream_power::StreamPowerConfig;
        // ⚠️ POPULATION. `StreamPowerConfig::default()` is k = 1, dt = 1, iterations = 4 — NOT the
        // shipped relief. The delivered config is `relief_v3(cell_km2, depth_scale_m)`, and the
        // first run of this block read the default and its own assert caught it at k_time = 4.
        let sp = StreamPowerConfig::relief_v3(CELL_KM * CELL_KM, ss.depth_scale_m as f32);
        let a_max = 1611.0f32; // Finding 44's own input: max land accumulation at 8192²
        let plan = sp.timescale_plan(a_max);
        eprintln!(
            "\n── C0 · `timescale_plan()` at the shipped config, A_max = {a_max} km² ──\n   \
             k {} · dt {} · iterations {} · **k_time {}** (shipped {SHIPPED_K_TIME}) · celerity \
             {:.0} m/yr · dt_max **{:.3e} yr** · **cfl_iterations {:.0}** · **Courant {:.0}**",
            sp.k,
            sp.dt,
            sp.iterations,
            sp.k_time(),
            plan.celerity_m_per_yr,
            plan.dt_max_yr,
            plan.cfl_iterations,
            plan.courant
        );
        assert!(
            (sp.k_time() - SHIPPED_K_TIME).abs() < 1e-3,
            "C0 control: k_time must be {SHIPPED_K_TIME}"
        );
    }

    let pre = build_field(Knobs::no_incision());
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let m = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
    let breach_of = |g: &GridF32| -> GridF32 {
        let d = c1_drainage_windowed(g, None, &on, &ss, DOMAIN_KM);
        breach_monotone(g, &d.flow.filled, &d.lake_map, SEA, w, h)
    };
    let br_pre = breach_of(&pre);
    let pre_spurs_unbr = spurs_km(&land_u16(&pre, &ss), w, 1);
    let pre_spurs_br = spurs_km(&land_u16(&br_pre, &ss), w, 1);
    eprintln!(
        "   the authority: pre-incision spurs ≥ 1 km, u16 8192² — **unbreached {pre_spurs_unbr}** \
         (Finding 80's twenty) · breached {pre_spurs_br}"
    );

    // the delivered reference for the paired hypsometry (X = ±10 %, Finding 89-B4)
    let raw2 = build_field(Knobs::passes(2));
    let ref_common: Vec<usize> =
        (0..n).filter(|&k| raw2.data[k] > SEA && pre.data[k] > SEA).collect();
    let ref_p50 = median(&ref_common.iter().map(|&k| m(&raw2, k)).collect::<Vec<f32>>());
    eprintln!("   the delivered reference p50 (paired) = **{ref_p50:.1} m**, tolerance ±10 %");

    for &steps in &[2usize, 10, 100] {
        let t0 = Instant::now();
        let f = if steps == 2 { raw2.clone() } else { build_field(Knobs::integrating(steps)) };
        let secs = t0.elapsed().as_secs_f64();
        let k_used = SHIPPED_K_TIME / steps as f32;
        eprintln!(
            "\n╔══════ MILESTONE {steps} passes · k = {k_used:.2} · k_time = {:.1} · Courant \
             {:.1} · built in {secs:.1} s ══════╗",
            k_used * steps as f32,
            3699.0 * 2.0 / steps as f32
        );
        let bf = breach_of(&f);

        // ── the coast criterion, BOTH stages (Finding 90-A: the stage is the criterion) ──
        let (se, sb) = (
            spurs_km(&(0..n).map(|k| f.data[k] > SEA).collect::<Vec<bool>>(), w, 1),
            spurs_km(&land_u16(&bf, &ss), w, 1),
        );
        let pre_e = spurs_km(&(0..n).map(|k| pre.data[k] > SEA).collect::<Vec<bool>>(), w, 1);
        eprintln!(
            "   COAST ≥ 1 km, 8192²: **eroded {se} (Δ {:+} vs pre-incision {pre_e})** · \
             **breached+u16 {sb} (Δ {:+} vs the {pre_spurs_unbr} authority)**",
            se as i64 - pre_e as i64,
            sb as i64 - pre_spurs_unbr as i64
        );
        // ── the breach's drowning, the mechanism of the added spurs ──
        let drowned = (0..n).filter(|&k| f.data[k] > SEA && bf.data[k] <= SEA).count();
        eprintln!(
            "   the breach drowns **{drowned}** cells ({:.1} km²)",
            drowned as f32 * cell_km2
        );

        // ── paired relief against the DELIVERED reference ──
        let com: Vec<usize> = (0..n).filter(|&k| f.data[k] > SEA && pre.data[k] > SEA).collect();
        let dd = sorted(com.iter().map(|&k| m(&f, k)).collect::<Vec<f32>>());
        let cut: Vec<f32> = com.iter().map(|&k| m(&pre, k) - m(&f, k)).collect();
        let p50 = pct(&dd, 0.50);
        eprintln!(
            "   RELIEF paired (n = {}): p10/p50/p90 **{:.1} / {:.1} / {:.1} m** · median cut \
             {:.1} m · **{:+.1} % against the delivered {ref_p50:.1} m ⇒ {}**",
            com.len(),
            pct(&dd, 0.10),
            p50,
            pct(&dd, 0.90),
            median(&cut),
            100.0 * (p50 - ref_p50) / ref_p50,
            if (100.0 * (p50 - ref_p50) / ref_p50).abs() <= 10.0 {
                "INSIDE ±10 %"
            } else {
                "**OUTSIDE ±10 %**"
            }
        );

        // ── the lake inventory: fraction, canyon class, integration ──
        let climate = c1_climate_placed(&bf, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dclim = DrainageClimate {
            precip_internal: &climate.precipitation,
            temperature: &climate.temperature,
        };
        let pre_d = c1_drainage_windowed(&f, None, &on, &ss, DOMAIN_KM);
        let mut dr = c1_drainage_windowed(&bf, Some(&dclim), &on, &ss, DOMAIN_KM);
        dr.lakes = pre_d.lakes;
        dr.lake_map = pre_d.lake_map;
        let li = std::mem::take(&mut dr.lakes);
        dr.lakes = apply_lake_water_balance(
            &bf,
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
        let racc = runoff_accumulation(&bf, &dr.flow, &dclim, cell_km2, None, None, w, h);
        let pm = dr.lake_map.clone();
        let bs = below_sea_basin_lakes_infil(&bf, &dclim, &on, &ss, DOMAIN_KM, Some(&pm), None);
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
        let absorbed: HashSet<u32> = d_before
            .keys()
            .copied()
            .filter(|id| d_after.get(id).copied().unwrap_or(0) == 0)
            .collect();
        dr.lakes.retain(|l| l.base.id >= 1_000_001 || !absorbed.contains(&l.base.id));
        dr.lakes.extend(bs.lakes.iter().cloned());
        clip_rivers_to_lakes(&mut dr);
        let inv: HashSet<u32> = dr.lakes.iter().map(|l| l.base.id).collect();
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
                source_lake: inv.contains(&sw.lake_id).then_some(sw.lake_id),
            });
        }
        let relab = resolve_exorheic_without_outlet(&mut dr);

        let land_km2 = (0..n).filter(|&k| bf.data[k] > SEA).count() as f32 * cell_km2;
        let water_cells = (0..n).filter(|&k| dr.lake_map[k] != 0).count();
        let water_km2 = water_cells as f32 * cell_km2;
        let mut klass = 0usize;
        let mut scanned = 0usize;
        let mut ty: HashMap<String, usize> = HashMap::new();
        for l in &dr.lakes {
            *ty.entry(format!("{:?}", l.lake_type)).or_default() += 1;
        }
        for l in dr.lakes.iter().filter(|l| l.area_km2 >= 1.0) {
            let cells: Vec<usize> = (0..n).filter(|&k| dr.lake_map[k] == l.base.id).collect();
            if cells.is_empty() {
                continue;
            }
            scanned += 1;
            let cuts: Vec<f32> = cells.iter().map(|&k| m(&pre, k) - m(&bf, k)).collect();
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
            let rim: Vec<f32> = ring.iter().map(|&k| slope_deg(&bf, k, n2m, w, h)).collect();
            if median(&cuts) > 50.0 && !rim.is_empty() && median(&rim) > 30.0 {
                klass += 1;
            }
        }
        let q_unres: f64 = dr
            .lakes
            .iter()
            .filter(|l| l.lake_type == LakeType::Unresolved)
            .map(|l| {
                let c: Vec<usize> = (0..n).filter(|&k| dr.lake_map[k] == l.base.id).collect();
                runoff_km2_to_m3s(c.iter().map(|&k| racc[k]).fold(0.0f32, f32::max)) as f64
            })
            .sum();
        let wc = water_class(&bf, SEA);
        eprintln!(
            "   LAKES: **{} bodies · {water_km2:.0} km² = {:.2} % of {land_km2:.0} km² of \
             land** · {:?} · CANYON CLASS **{klass} of {scanned}** · relabels {} · \
             `to_nothing` {} · Unresolved inflow {q_unres:.1} m³/s · wc==2 cells {}",
            dr.lakes.len(),
            100.0 * water_km2 / land_km2,
            ty,
            relab.len(),
            bs.termination.to_nothing,
            (0..n).filter(|&k| wc[k] == 2).count()
        );
        for (src, k_ours) in
            [("W&T", 2.00e-2f64), ("Harel", 7.50e-4), ("S&M high", 2.51e-4), ("S&M low", 2.51e-5)]
        {
            eprint!("   T({src}) = {:.2e} yr ·", SHIPPED_K_TIME as f64 / k_ours);
        }
        eprintln!(" (unchanged by construction — k_time is held)");
    }
    eprintln!("\n==========  end Finding 90 · C  ==========\n");
}
