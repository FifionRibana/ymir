//! ADR 0001 Finding 96 — a closure for the age of the continent. **No production change.**
//!
//! Direction B measured first, because block B0's grep settles what direction B IS: Finding 6
//! rejected *adding U as a source* and PRESCRIBED *limiting total incision*, so a chi-profile floor
//! is Finding 6's own remedy with a physically-derived limit instead of a smaller K. B1 builds the
//! profile, B2 prices what it loses, B4 prices what it saves; A2 sweeps `talus_passes`, which is a
//! `Knobs` field and therefore free.
//!
//! ⚠️ The ORACLE IS A ROW OF NUMBERS, not a field: rebuilding the 300-pass field costs 1 h 43 of one
//! core and the author has refused it. Finding 95 recorded canyons 0/24, `to_nothing` 0, Δ eroded
//! coast −2, paired relief p50 **488.1 m**, lake fraction 15.92 % — and table C compares against
//! those. Nothing here reconstructs the oracle.
//!
//! Run: cargo test -p ymir-core --release --test f96_closure -- --ignored --nocapture

mod common;

use common::{
    CELL_KM, Knobs, SEA, build_field, build_field_with_floor, land_u16, majority, pct, sorted,
    to_mask,
};
use std::collections::{HashMap, HashSet};
use std::time::Instant;
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
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
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, breach_monotone, propagation_order};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
/// Finding 95's oracle, as numbers.
const ORACLE_P50_M: f32 = 488.1;
const DELIVERED_P50_M: f32 = 424.2;

fn median(v: &[f32]) -> f32 {
    pct(&sorted(v.to_vec()), 0.50)
}

fn spurs_km(mask: &[bool], w: usize, factor: usize) -> usize {
    let (m, ww) = majority(mask, w, factor);
    let ck = CELL_KM * factor as f32;
    coast_spurs(&marching_squares(&to_mask(&m, ww), 0.5), ck, MIN_SPUR_KM, NECK_KM).0.len()
}

/// Max-of-8 slope in degrees, in metres.
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

/// Local standard deviation of altitude over a 3×3 window, in metres — the texture measure B2 needs.
fn local_sigma(f: &GridF32, land: &[bool], n2m: f32, w: usize, h: usize) -> Vec<f32> {
    let mut out = Vec::new();
    // every 37th land cell: a fixed stride, so the population is declared and reproducible
    for k in (0..w * h).step_by(37) {
        if !land[k] {
            continue;
        }
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        let (mut s, mut s2, mut c) = (0.0f64, 0.0f64, 0usize);
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let nx = (x + dx).rem_euclid(w as i32) as usize;
                let ny = (y + dy).rem_euclid(h as i32) as usize;
                let v = (f.data[ny * w + nx] * n2m) as f64;
                s += v;
                s2 += v * v;
                c += 1;
            }
        }
        let m = s / c as f64;
        out.push(((s2 / c as f64 - m * m).max(0.0)).sqrt() as f32);
    }
    out
}

/// Finding 95's criteria chain, **copied VERBATIM** from `f95_long_run.rs:159-310` (commit
/// `751e0bf`) so table C reads the SAME instruments as the oracle row. Deliberately NOT
/// refactored into a shared helper: `f95_long_run` costs 5 h and must not be touched to serve
/// this round. If the two ever drift, this is a COPY and the drift is a defect of this bench.
///
/// `f` = the eroded field before the breach - `bf` = its breach - `pre` = the pre-incision
/// reference the cut is measured against.
#[allow(clippy::too_many_arguments)]
fn f95_criteria(
    f: &GridF32,
    bf: &GridF32,
    pre: &GridF32,
    ref_p50: f32,
    ss: &SteinSteinParams,
    on: &C1DrainageConfig,
    cell_km2: f32,
    n2m: f32,
    w: usize,
    h: usize,
) -> Crit {
    let n = w * h;
    let m = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], ss);
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
    let absorbed: HashSet<u32> =
        d_before.keys().copied().filter(|id| d_after.get(id).copied().unwrap_or(0) == 0).collect();
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
    Crit {
        p50,
        klass,
        scanned,
        to_nothing: bs.termination.to_nothing,
        lake_pct: 100.0 * water_km2 / land_km2,
        unresolved: dr.lakes.iter().filter(|l| l.lake_type == LakeType::Unresolved).count(),
        q_unres,
    }
}

/// Finding 95's criteria, as one row of table C.
struct Crit {
    p50: f32,
    klass: usize,
    scanned: usize,
    to_nothing: usize,
    lake_pct: f32,
    unresolved: usize,
    q_unres: f64,
}

#[test]
#[ignore]
fn f96_closure() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    eprintln!("\n==========  Finding 96 · a closure for the age of the continent  ==========");

    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;

    // ── the three fields the table needs, and nothing more ──────────────────
    let t0 = Instant::now();
    let pre = build_field(Knobs::no_incision());
    let t_pre = t0.elapsed().as_secs_f64();
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let t1 = Instant::now();
    let raw = build_field(Knobs::passes(2));
    let t_del = t1.elapsed().as_secs_f64();
    let m = |g: &GridF32, k: usize| c1_altitude_norm_to_metres(g.data[k], &ss);
    eprintln!(
        "   builds: pre-incision {t_pre:.1} s · delivered (2 passes) {t_del:.1} s · 1 norm unit = \
         {n2m:.1} m · cell {CELL_M:.2} m"
    );

    // the breach of the PRE-INCISION field: the depression-free network chi is built on
    let tb = Instant::now();
    let d_pre = c1_drainage_windowed(&pre, None, &on, &ss, DOMAIN_KM);
    let br_pre = breach_monotone(&pre, &d_pre.flow.filled, &d_pre.lake_map, SEA, w, h);
    let t_breach_pre = tb.elapsed().as_secs_f64();
    let tb2 = Instant::now();
    let d_br = c1_drainage_windowed(&br_pre, None, &on, &ss, DOMAIN_KM);
    let t_flow_br = tb2.elapsed().as_secs_f64();
    eprintln!(
        "   the chi network: breach of the pre-incision field {t_breach_pre:.1} s, its drainage \
         {t_flow_br:.1} s"
    );

    let land_pre: Vec<bool> = (0..n).map(|k| br_pre.data[k] > SEA).collect();
    let del_field = {
        let d = c1_drainage_windowed(&raw, None, &on, &ss, DOMAIN_KM);
        breach_monotone(&raw, &d.flow.filled, &d.lake_map, SEA, w, h)
    };

    // ══ B1 — the chi integral, one ordered descent ══════════════════════════
    //
    // chi(k) = chi(receiver) + (A0/A(k))^(m/n) * dx, with chi = 0 at the sea, so
    // z_chi(k) = z_base + C * chi(k) with C = (U/K)^(1/n). `propagation_order` (Finding 74) gives
    // the traversal: donors before receivers, so ONE pass in reverse order fills chi.
    let tc = Instant::now();
    let (order, unreached) = propagation_order(&d_br.flow.direction, |k| land_pre[k], w, h);
    let m_over_n = 0.5f32; // Harel, anchored (Finding 89-C6's reference)
    let a0_km2 = 1.0f32; // the reference area of the chi integral, declared
    let mut chi = vec![0.0f32; n];
    // `order` is donors-first; chi needs receivers first, so walk it in REVERSE
    for &ku in order.iter().rev() {
        let k = ku as usize;
        if !land_pre[k] {
            continue;
        }
        let d = d_br.flow.direction[k];
        if d == DIR_NONE {
            continue;
        }
        let nx = ((k % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
        let ny = ((k / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
        let r = ny * w + nx;
        let diag = D8_DX[d as usize] != 0 && D8_DY[d as usize] != 0;
        let dx_m = if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M };
        let a_km2 = (d_br.flow.accumulation.data[k] * cell_km2).max(1e-6);
        let up = if land_pre[r] { chi[r] } else { 0.0 };
        chi[k] = up + (a0_km2 / a_km2).powf(m_over_n) * dx_m;
    }
    let t_chi = tc.elapsed().as_secs_f64();
    let chi_land: Vec<f32> = (0..n).filter(|&k| land_pre[k]).map(|k| chi[k]).collect();
    let chi_s = sorted(chi_land.clone());
    eprintln!(
        "\n── B1 · the chi integral ──\n   **built in {t_chi:.2} s** ({} land cells, {} unreached \
         — a cycle if > 0) · chi (m) p10 {:.0} p50 **{:.0}** p90 {:.0} max {:.0}",
        chi_land.len(),
        unreached,
        pct(&chi_s, 0.10),
        pct(&chi_s, 0.50),
        pct(&chi_s, 0.90),
        chi_s.last().copied().unwrap_or(0.0)
    );
    assert_eq!(unreached, 0, "B1: a cycle in the chi network invalidates the integral");

    // calibrate the single free parameter so the PAIRED p50 equals the oracle's 488.1 m
    let common: Vec<usize> =
        (0..n).filter(|&k| br_pre.data[k] > SEA && pre.data[k] > SEA).collect();
    let chi_p50_paired = median(&common.iter().map(|&k| chi[k]).collect::<Vec<f32>>());
    let c_uk = ORACLE_P50_M / chi_p50_paired.max(1e-6);
    eprintln!(
        "   calibration: paired chi p50 = {chi_p50_paired:.1} m of integral ⇒ \
         **(U/K)^(1/n) = {c_uk:.4}** to put the paired relief p50 at the oracle's {ORACLE_P50_M} m"
    );
    // the coast authority, hoisted: B1' and B2 both read it
    let auth = spurs_km(&land_u16(&pre, &ss), w, 1);

    // ══ B1' — THE REPAIR: chi over CHANNEL cells only ══════════════════════
    //
    // The first pass built chi over ALL land, and B2 measured the consequence: 79 % of the field
    // steeper than 30°, local sigma 7× the delivered field. The cause is arithmetic, not a bug:
    // the profile's gradient is (U/K)^(1/n)·(A0/A)^(m/n), so at ONE CELL of drainage area
    // (0.00238 km²) with A0 = 1 km² that is 0.0646 × 20.5 = 1.32 m/m = **53°**, BY CONSTRUCTION.
    // The fluvial law is being applied three decades below its own channel threshold — Finding 65
    // definition 1, `A_c` = 0.1 km², the fluvial/debris-flow break the anchoring report calls
    // ANCHORED (Whipple & Tucker Table 1: 0.059–0.140 km²).
    //
    // The repair is one clamp: integrate chi with A floored at `A_c`, so every cell below the
    // channel threshold inherits the gradient AT the threshold instead of a divergent one.
    let a_c_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    let mut chi_c = vec![0.0f32; n];
    let tc2 = Instant::now();
    for &ku in order.iter().rev() {
        let k = ku as usize;
        if !land_pre[k] {
            continue;
        }
        let d = d_br.flow.direction[k];
        if d == DIR_NONE {
            continue;
        }
        let nx = ((k % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
        let ny = ((k / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
        let r = ny * w + nx;
        let diag = D8_DX[d as usize] != 0 && D8_DY[d as usize] != 0;
        let dx_m = if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M };
        let a_km2 = (d_br.flow.accumulation.data[k] * cell_km2).max(a_c_km2);
        let up = if land_pre[r] { chi_c[r] } else { 0.0 };
        chi_c[k] = up + (a0_km2 / a_km2).powf(m_over_n) * dx_m;
    }
    let t_chi_c = tc2.elapsed().as_secs_f64();
    let chi_c_p50 = median(&common.iter().map(|&k| chi_c[k]).collect::<Vec<f32>>());
    let c_uk_c = ORACLE_P50_M / chi_c_p50.max(1e-6);
    eprintln!(
        "
── B1' · chi CLAMPED at A_c = {a_c_km2} km² ──
   **built in {t_chi_c:.2} s** · paired          chi p50 {chi_c_p50:.1} m ⇒ (U/K)^(1/n) = **{c_uk_c:.4}** · the gradient at one cell is now          {:.3} m/m = **{:.1}°** instead of {:.3} m/m = {:.1}°",
        c_uk_c * (a0_km2 / a_c_km2).sqrt(),
        (c_uk_c * (a0_km2 / a_c_km2).sqrt()).atan().to_degrees(),
        c_uk * (a0_km2 / cell_km2).sqrt(),
        (c_uk * (a0_km2 / cell_km2).sqrt()).atan().to_degrees()
    );
    let z_chi_c = GridF32 {
        width: w,
        height: h,
        data: (0..n)
            .map(|k| if land_pre[k] { SEA + (c_uk_c * chi_c[k]) / n2m } else { br_pre.data[k] })
            .collect(),
    };
    {
        let land_c: Vec<bool> = (0..n).map(|k| z_chi_c.data[k] > SEA).collect();
        let sig = sorted(local_sigma(&z_chi_c, &land_c, n2m, w, h));
        let sl = sorted(
            (0..n)
                .step_by(37)
                .filter(|&k| land_c[k])
                .map(|k| slope_deg(&z_chi_c, k, n2m, w, h))
                .collect(),
        );
        let alt = sorted(common.iter().map(|&k| m(&z_chi_c, k)).collect::<Vec<f32>>());
        let dch = c1_drainage_windowed(&z_chi_c, None, &on, &ss, DOMAIN_KM);
        let pits = (0..n)
            .filter(|&k| land_c[k] && dch.flow.filled.data[k] > z_chi_c.data[k] + 1e-6)
            .count();
        eprintln!(
            "   σ p50 **{:.2} m** (all-land chi 44.73, delivered 6.28) · slope p50 **{:.1}°** p90              **{:.1}°** max {:.1}° · share > 30° **{:.2} %** (all-land chi 79.07, delivered 28.53)              · paired p10/p50/p90 **{:.1} / {:.1} / {:.1} m** · pits **{pits}** · spurs ≥ 1 km              **{} (Δ {:+})**",
            pct(&sig, 0.50),
            pct(&sl, 0.50),
            pct(&sl, 0.90),
            sl.last().copied().unwrap_or(0.0),
            100.0 * sl.iter().filter(|&&s| s > 30.0).count() as f64 / sl.len().max(1) as f64,
            pct(&alt, 0.10),
            pct(&alt, 0.50),
            pct(&alt, 0.90),
            spurs_km(&land_u16(&z_chi_c, &ss), w, 1),
            spurs_km(&land_u16(&z_chi_c, &ss), w, 1) as i64 - auth as i64
        );
    }

    let z_chi = GridF32 {
        width: w,
        height: h,
        data: (0..n)
            .map(|k| if land_pre[k] { SEA + (c_uk * chi[k]) / n2m } else { br_pre.data[k] })
            .collect(),
    };

    // the Flint control: log S vs log A on the CONSTRUCTED field must give theta = m/n
    {
        let (mut sx, mut sy, mut sxx, mut sxy, mut cnt) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0usize);
        for k in (0..n).step_by(17) {
            if !land_pre[k] {
                continue;
            }
            let d = d_br.flow.direction[k];
            if d == DIR_NONE {
                continue;
            }
            let nx = ((k % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
            let ny = ((k / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
            let r = ny * w + nx;
            let diag = D8_DX[d as usize] != 0 && D8_DY[d as usize] != 0;
            let dx_m = if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M };
            let s = (m(&z_chi, k) - m(&z_chi, r)) / dx_m;
            let a = d_br.flow.accumulation.data[k] * cell_km2;
            if s <= 0.0 || a <= 0.0 {
                continue;
            }
            let (lx, ly) = ((a as f64).ln(), (s as f64).ln());
            sx += lx;
            sy += ly;
            sxx += lx * lx;
            sxy += lx * ly;
            cnt += 1;
        }
        let nn = cnt as f64;
        let slope = (nn * sxy - sx * sy) / (nn * sxx - sx * sx);
        let intercept = (sy - slope * sx) / nn;
        eprintln!(
            "   **FLINT CONTROL** on the constructed field ({cnt} samples): log S = {intercept:.4} \
             {slope:+.4}·log A ⇒ **θ = {:.4}** (must be {m_over_n}) · k_s = exp(intercept) = \
             {:.4e} (must be (U/K)^(1/n) = {c_uk:.4}) ⇒ **{}**",
            -slope,
            intercept.exp(),
            if (-slope - m_over_n as f64).abs() < 0.02 { "PASSES" } else { "**FAILS**" }
        );
        assert!(
            (-slope - m_over_n as f64).abs() < 0.05,
            "B1: the Flint control rejects the integral (theta = {:.4})",
            -slope
        );
    }

    // ── the chi field's own numbers ─────────────────────────────────────────
    let chi_alt: Vec<f32> = common.iter().map(|&k| m(&z_chi, k)).collect();
    let chi_alt_s = sorted(chi_alt.clone());
    let pre_alt_s = sorted(common.iter().map(|&k| m(&pre, k)).collect::<Vec<f32>>());
    let del_alt_s = sorted(common.iter().map(|&k| m(&del_field, k)).collect::<Vec<f32>>());
    eprintln!(
        "   paired hypsometry (n = {}): chi p10/p50/p90 **{:.1} / {:.1} / {:.1} m** · delivered \
         {:.1} / {:.1} / {:.1} · pre-incision {:.1} / {:.1} / {:.1}",
        common.len(),
        pct(&chi_alt_s, 0.10),
        pct(&chi_alt_s, 0.50),
        pct(&chi_alt_s, 0.90),
        pct(&del_alt_s, 0.10),
        pct(&del_alt_s, 0.50),
        pct(&del_alt_s, 0.90),
        pct(&pre_alt_s, 0.10),
        pct(&pre_alt_s, 0.50),
        pct(&pre_alt_s, 0.90)
    );
    // depressions: the whole point of the construction
    let d_chi = c1_drainage_windowed(&z_chi, None, &on, &ss, DOMAIN_KM);
    let pits_chi =
        (0..n).filter(|&k| land_pre[k] && d_chi.flow.filled.data[k] > z_chi.data[k] + 1e-6).count();
    let pits_del = {
        let d = c1_drainage_windowed(&raw, None, &on, &ss, DOMAIN_KM);
        (0..n).filter(|&k| raw.data[k] > SEA && d.flow.filled.data[k] > raw.data[k] + 1e-6).count()
    };
    eprintln!(
        "   ⛔ **pits (a cell the priority flood must raise): chi {pits_chi} · delivered \
         {pits_del}** ⇒ chi is {} depression-free",
        if pits_chi * 100 < pits_del { "effectively" } else { "NOT" }
    );

    // ── B2 — what chi loses: texture ────────────────────────────────────────
    let land_del: Vec<bool> = (0..n).map(|k| del_field.data[k] > SEA).collect();
    let sig_chi = sorted(local_sigma(&z_chi, &land_pre, n2m, w, h));
    let sig_del = sorted(local_sigma(&del_field, &land_del, n2m, w, h));
    let sig_pre = sorted(local_sigma(&pre, &land_pre, n2m, w, h));
    eprintln!(
        "\n── B2 · the price: local 3×3 σ of altitude, metres (stride 37, land only) ──\n   \
         chi p50 **{:.2}** p90 {:.2} · delivered p50 **{:.2}** p90 {:.2} · pre-incision p50 \
         {:.2} p90 {:.2} ⇒ chi is **{:.1}× ROUGHER** than the delivered field at the cell          scale — the price is the OPPOSITE of the predicted one (ADR Finding 96 B2)",
        pct(&sig_chi, 0.50),
        pct(&sig_chi, 0.90),
        pct(&sig_del, 0.50),
        pct(&sig_del, 0.90),
        pct(&sig_pre, 0.50),
        pct(&sig_pre, 0.90),
        pct(&sig_chi, 0.50) / pct(&sig_del, 0.50).max(1e-6)
    );
    let sl_chi = sorted(
        (0..n)
            .step_by(37)
            .filter(|&k| land_pre[k])
            .map(|k| slope_deg(&z_chi, k, n2m, w, h))
            .collect(),
    );
    let sl_del = sorted(
        (0..n)
            .step_by(37)
            .filter(|&k| land_del[k])
            .map(|k| slope_deg(&del_field, k, n2m, w, h))
            .collect(),
    );
    eprintln!(
        "   slope° p50/p90/max: chi **{:.1} / {:.1} / {:.1}** · delivered **{:.1} / {:.1} / {:.1}** \
         · share above 30°: chi **{:.2} %** · delivered **{:.2} %**",
        pct(&sl_chi, 0.50),
        pct(&sl_chi, 0.90),
        sl_chi.last().copied().unwrap_or(0.0),
        pct(&sl_del, 0.50),
        pct(&sl_del, 0.90),
        sl_del.last().copied().unwrap_or(0.0),
        100.0 * sl_chi.iter().filter(|&&s| s > 30.0).count() as f64 / sl_chi.len().max(1) as f64,
        100.0 * sl_del.iter().filter(|&&s| s > 30.0).count() as f64 / sl_del.len().max(1) as f64
    );
    // the coast, on the criterion's own instrument
    let sp_chi = spurs_km(&land_u16(&z_chi, &ss), w, 1);
    let sp_del = spurs_km(&land_u16(&del_field, &ss), w, 1);
    let sp_del_ero = spurs_km(&(0..n).map(|k| raw.data[k] > SEA).collect::<Vec<bool>>(), w, 1);
    let sp_pre_ero = spurs_km(&(0..n).map(|k| pre.data[k] > SEA).collect::<Vec<bool>>(), w, 1);
    eprintln!(
        "   coast, spurs ≥ 1 km at 8192²: authority (pre-incision u16) **{auth}** · chi **{sp_chi} \
         (Δ {:+})** · delivered breached+u16 {sp_del} (Δ {:+}) · delivered ERODED {sp_del_ero} \
         (Δ {:+} vs {sp_pre_ero})",
        sp_chi as i64 - auth as i64,
        sp_del as i64 - auth as i64,
        sp_del_ero as i64 - sp_pre_ero as i64
    );

    // ── B4 — what chi saves: the breach has nothing to carve ────────────────
    let tb4 = Instant::now();
    let _ = breach_monotone(&z_chi, &d_chi.flow.filled, &d_chi.lake_map, SEA, w, h);
    let t_breach_chi = tb4.elapsed().as_secs_f64();
    let tb5 = Instant::now();
    let d_del = c1_drainage_windowed(&raw, None, &on, &ss, DOMAIN_KM);
    let _ = breach_monotone(&raw, &d_del.flow.filled, &d_del.lake_map, SEA, w, h);
    let t_breach_del = tb5.elapsed().as_secs_f64();
    eprintln!(
        "\n── B4 · the breach ──\n   on the chi field **{t_breach_chi:.1} s** · on the delivered \
         eroded field **{t_breach_del:.1} s** (the drainage call is inside both) ⇒ **{:+.1} s**",
        t_breach_chi - t_breach_del
    );

    // ══ A2 — MEASURED IN THE FIRST PASS, not repeated (403 s of builds). Results kept in the
    // finding: talus_passes 4 / 16 / 64 → slope p90 44.6 / 37.1 / 33.7°, max 85.2 / 86.4 / 86.3°,
    // share > 30° 28.82 / 31.28 / 32.96 %, paired p50 424.2 / 430.2 / 436.5 m, cost 45.5 / 89.2 /
    // 268.7 s. The tail comes down, the SHARE goes UP, the max never moves, and 64 passes cost
    // 6.5× the delivered incision.

    // == A1 - the closed-depression floor, and B3 - the chi floor ===========
    //
    // Both go through the SAME seam (`stream_power.rs`, the Finding 80/83 site): the relaxation
    // TARGET is bounded at `min(floor, h_o)`, so each can stop the incision and neither can
    // deposit. A1's floor is `flow.filled` - the spill level of the depression the cell sits in,
    // recomputed every iteration and already paid for. B3's floor is the CLAMPED chi profile,
    // which is Finding 6's prescription - *"the fix is to LIMIT total incision, not add U"* -
    // with a limit whose shape comes from a physical profile instead of a smaller global K.
    let ta1 = Instant::now();
    let a1 = build_field(Knobs { depression_floor: true, ..Knobs::passes(2) });
    let t_a1 = ta1.elapsed().as_secs_f64();
    let tb3 = Instant::now();
    let b3 =
        build_field_with_floor(Knobs::passes(2), Some(std::sync::Arc::new(z_chi_c.data.clone())));
    let t_b3 = tb3.elapsed().as_secs_f64();
    // A1 + B3 together. The hypothesis this tests: B3's canyon class survives because chi was
    // built on the BREACHED pre-incision field, where the future lake footprints have been carved
    // through -- so inside them chi is LOW and the floor permits the very cut it should stop.
    // A1's exclusion covers exactly that population, so the union should fix what B3 alone cannot.
    let tab = Instant::now();
    let ab = build_field_with_floor(
        Knobs { depression_floor: true, ..Knobs::passes(2) },
        Some(std::sync::Arc::new(z_chi_c.data.clone())),
    );
    let t_ab = tab.elapsed().as_secs_f64();
    eprintln!(
        "\n-- A1 / B3 . the two floors, built --\n   A1 (depression skip) **{t_a1:.1} s** . B3 \
         (chi floor, (U/K)^(1/n) = {c_uk_c:.4}) **{t_b3:.1} s** . the delivered reference built \
         in {t_del:.1} s => the closures cost **{:+.1} s**, **{:+.1} s** and **{:+.1} s** (union          {t_ab:.1} s)",
        t_a1 - t_del,
        t_b3 - t_del,
        t_ab - t_del
    );

    // the breach of each candidate, then Finding 95's criteria on the pair
    let mut rows: Vec<(&str, Crit, f64, f32, f32, usize)> = Vec::new();
    for (label, fld, cost) in [
        ("delivered", &raw, t_del),
        ("A1 depression skip", &a1, t_a1),
        ("B3 chi floor", &b3, t_b3),
        ("A1+B3 union", &ab, t_ab),
    ] {
        let d = c1_drainage_windowed(fld, None, &on, &ss, DOMAIN_KM);
        let bf = breach_monotone(fld, &d.flow.filled, &d.lake_map, SEA, w, h);
        let c = f95_criteria(fld, &bf, &pre, DELIVERED_P50_M, &ss, &on, cell_km2, n2m, w, h);
        let land: Vec<bool> = (0..n).map(|k| fld.data[k] > SEA).collect();
        let sl = sorted(
            (0..n).step_by(37).filter(|&k| land[k]).map(|k| slope_deg(fld, k, n2m, w, h)).collect(),
        );
        let pits =
            (0..n).filter(|&k| land[k] && d.flow.filled.data[k] > fld.data[k] + 1e-6).count();
        let sp_br = spurs_km(&land_u16(&bf, &ss), w, 1);
        let sp_er = spurs_km(&land_u16(fld, &ss), w, 1);
        eprintln!(
            "   **{label}**: pits **{pits}** . slope p50/p90/max **{:.1} / {:.1} / {:.1} deg** . \
             share > 30 deg **{:.2} %** . spurs >= 1 km on the BREACHED u16 **{sp_br} ({:+})** . \
             ERODED, u16-QUANTISED **{sp_er} ({:+})** -- NOT the same stage as B2's \
             float-threshold reading of the same field, which gives 19 (delta -1). Both are \
             legitimate readings; only the label was ambiguous (ADR Finding 96, stage discipline)",
            pct(&sl, 0.50),
            pct(&sl, 0.90),
            sl.last().copied().unwrap_or(0.0),
            100.0 * sl.iter().filter(|&&x| x > 30.0).count() as f64 / sl.len().max(1) as f64,
            sp_br as i64 - auth as i64,
            sp_er as i64 - auth as i64
        );
        let d30 = 100.0 * sl.iter().filter(|&&x| x > 30.0).count() as f32 / sl.len().max(1) as f32;
        rows.push((label, c, cost, pct(&sl, 0.90), d30, pits));
    }

    // == C - the table that decides, filled =================================
    //
    // The oracle is Finding 95's 300-pass row: canyon class 0, `to_nothing` 0, delta eroded coast
    // -2, paired relief p50 488.1 m, lakes 15.92 %. The gates the round set: within 20 % of the
    // oracle on the numeric criteria, inside +-10 % of relief, under 30 s of cost.
    eprintln!("\n-- C . the decision table --");
    eprintln!(
        "   | candidate | canyons | to_nothing | lakes % | Unresolved | paired p50 | vs oracle | \
         slope p90 | pits | cost |"
    );
    eprintln!("   |---|---|---|---|---|---|---|---|---|---|");
    for (label, c, cost, p90, _d30, pits) in &rows {
        eprintln!(
            "   | {label} | {}/{} | {} | {:.2} | {} | {:.1} m | {:+.1} % | {:.1} deg | {pits} | \
             {cost:.1} s |",
            c.klass,
            c.scanned,
            c.to_nothing,
            c.lake_pct,
            c.unresolved,
            c.p50,
            100.0 * (c.p50 - ORACLE_P50_M) / ORACLE_P50_M,
            p90
        );
    }
    eprintln!(
        "   | **ORACLE (F95, 300 passes)** | **0**/24 | **0** | **15.92** | n/a | \
         **{ORACLE_P50_M:.1} m** | - | n/a | n/a | **6 210 s** |"
    );
    for (label, c, cost, _p90, d30, _pits) in &rows {
        let dp = 100.0 * (c.p50 - ORACLE_P50_M) / ORACLE_P50_M;
        eprintln!(
            "   {label}: relief {} . canyons {} . to_nothing {} . cost {} . texture (share > 30 \
             deg, delivered 28.53 %) {d30:.2} % . Unresolved inflow {:.1} m3/s => **{}**",
            if dp.abs() <= 10.0 { "INSIDE +-10 %" } else { "**OUTSIDE +-10 %**" },
            if c.klass == 0 { "PASS" } else { "**FAIL**" },
            if c.to_nothing == 0 { "PASS" } else { "**FAIL**" },
            if *cost <= 30.0 { "PASS" } else { "**FAIL (> 30 s)**" },
            c.q_unres,
            if dp.abs() <= 10.0 && c.klass == 0 && c.to_nothing == 0 && *cost <= 30.0 {
                "WINS"
            } else {
                "does NOT win"
            }
        );
    }
    eprintln!("\n==========  end Finding 96  ==========\n");
}
