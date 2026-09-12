//! ADR Finding 73/74 — the guards around the accumulation-order fix.
//!
//! Three things live here, and the first one runs BY DEFAULT on purpose:
//!
//! 1. `flat_carry_is_guarded` — a synthetic filled pit fed by a known tributary. Five
//!    accumulations are read across the same flat: the shared topological order, a
//!    deliberately broken order (the pre-fix `filled`-only sort, the NEGATIVE CONTROL — it
//!    must fail), the production geometric accumulation, `runoff_accumulation`, and
//!    `mfd_accumulation`. Cheap, so it protects on every `cargo test`. The 8192² guard below
//!    cannot: "a rule that is not executable does not protect" (ADR, Finding 56's lesson), and
//!    a four-minute test is not executable in a default suite.
//!
//! 2. `forty_two_lake_carry_on_the_production_path` (`#[ignore]`) — the Finding 73 instrument
//!    promoted to a regression guard, on production's own `assemble_hd_drainage`, at the
//!    thresholds declared in Finding 73: median OUT/IN ≥ 1.00 and ≥ 80 % of the population at
//!    ≥ 0.95.
//!
//! 3. `the_height_field_does_not_move` (`#[ignore]`) — the rule-7 block: the fix must change
//!    NOTHING upstream of the drainage. Prints a fingerprint of the eroded and breached
//!    fields plus the depression and below-sea inventories, and pins them.
//!
//! Run the slow ones: cargo test -p ymir-core --release --test accumulation_order -- --ignored --nocapture

use ymir_core::climate::precipitation::precip_mm_per_year;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::drainage::{
    C1_SEA_LEVEL_NORM, DrainageClimate, SegmentKind, potential_evaporation_mm, runoff_accumulation,
};
use ymir_core::terrain::flow::{
    D8_DX, D8_DY, DIR_NONE, FlowConfig, compute_flow, mfd_accumulation, propagation_order,
};

// ── the synthetic fixture ───────────────────────────────────────────────────────────────
//
// A single channel running west→east, with a flat-floored CLOSED BASIN across it. The basin
// floor sits above sea level, so it is land; `pit_fill` raises it to its sill and the result
// is an EXACTLY flat surface (Finding 13: "fills to the EXACT sill — no epsilon increment,
// so flats are truly flat"). That flat is the object under test.

const FW: usize = 96;
const FH: usize = 48;
const SEA: f32 = 0.5;
const PIT_X: std::ops::RangeInclusive<usize> = 40..=55;
const PIT_Y: std::ops::RangeInclusive<usize> = 18..=29;

fn fixture() -> GridF32 {
    let mut g = GridF32::new(FW, FH, 0.0);
    for y in 0..FH {
        for x in 0..FW {
            // a V-shaped valley centred on y = 24, descending eastward
            let v = 0.70 + (y as f32 - 24.0).abs() * 0.0020 - x as f32 * 0.0008;
            g.data[y * FW + x] = v;
        }
    }
    // the closed basin: a flat floor well below its rim, still above sea level
    for y in PIT_Y {
        for x in PIT_X {
            g.data[y * FW + x] = 0.55;
        }
    }
    // the ocean, so there is somewhere for the water to go
    for y in 0..FH {
        for x in 92..FW {
            g.data[y * FW + x] = 0.30;
        }
    }
    g
}

/// The cells `pit_fill` had to RAISE: the filled depression, i.e. the flat.
fn flat_of(raw: &GridF32, filled: &GridF32) -> Vec<bool> {
    (0..raw.data.len()).map(|k| filled.data[k] > raw.data[k] + 1e-7).collect()
}

fn receiver(direction: &[u8], k: usize, w: usize, h: usize) -> Option<usize> {
    let d = direction[k];
    if d == DIR_NONE {
        return None;
    }
    let (i, j) = (k % w, k / w);
    let ni = ((i as i32 + D8_DX[d as usize]).rem_euclid(w as i32)) as usize;
    let nj = ((j as i32 + D8_DY[d as usize]).rem_euclid(h as i32)) as usize;
    Some(nj * w + ni)
}

/// THE METRIC, one definition for every reading: the largest accumulation ENTERING the flat
/// against the largest accumulation LEAVING it. A carrier that works reads ≥ 1 (the exit is
/// downstream, so it also holds the flat's own area); a carrier that drops the water on the
/// flat reads ≈ 0.
fn carry(acc: &[f32], flat: &[bool], direction: &[u8], w: usize, h: usize) -> (f32, f32, f32) {
    let (mut qin, mut qout) = (0.0f32, 0.0f32);
    for k in 0..w * h {
        match receiver(direction, k, w, h) {
            Some(r) if flat[r] && !flat[k] => qin = qin.max(acc[k]),
            Some(r) if flat[k] && !flat[r] => qout = qout.max(acc[r]),
            _ => {}
        }
    }
    (qin, qout, qout / qin.max(1e-9))
}

/// Accumulate unit seeds over `direction` in the given order — the loop every accumulation
/// in this crate runs; ONLY the order differs between the readings below.
fn accumulate(order: &[u32], land: &[bool], direction: &[u8], w: usize, h: usize) -> Vec<f32> {
    let mut acc: Vec<f32> = land.iter().map(|&l| if l { 1.0 } else { 0.0 }).collect();
    for &kk in order {
        let k = kk as usize;
        if let Some(r) = receiver(direction, k, w, h) {
            acc[r] += acc[k];
        }
    }
    acc
}

#[test]
fn flat_carry_is_guarded() {
    let raw = fixture();
    let flow = compute_flow(&raw, &FlowConfig { sea_level: SEA, ..Default::default() });
    let flat = flat_of(&raw, &flow.filled);
    let n_flat = flat.iter().filter(|&&f| f).count();
    let land: Vec<bool> = (0..FW * FH).map(|k| raw.data[k] > SEA).collect();
    let n_land = land.iter().filter(|&&l| l).count();
    eprintln!(
        "\n── synthetic fixture ── {FW}×{FH}, {n_land} land cells, a filled basin of \
         **{n_flat} cells** (the flat under test)"
    );
    assert!(n_flat >= 100, "the fixture must actually contain a flat; got {n_flat}");

    // 1. the shared topological order (the fix)
    let (order, n_cyclic) = propagation_order(&flow.direction, |k| land[k], FW, FH);
    assert_eq!(n_cyclic, 0, "no cycle expected on this fixture");
    assert_eq!(order.len(), n_land, "every land cell must be in the order");
    let acc_ok = accumulate(&order, &land, &flow.direction, FW, FH);
    let (in_ok, out_ok, r_ok) = carry(&acc_ok, &flat, &flow.direction, FW, FH);

    // 2. THE NEGATIVE CONTROL — the pre-fix order: `filled` descending and nothing else.
    //    On a flat every key is equal, so this is an arbitrary order.
    let mut broken: Vec<u32> = (0..FW * FH).filter(|&k| land[k]).map(|k| k as u32).collect();
    broken.sort_unstable_by(|&a, &b| {
        flow.filled.data[b as usize]
            .partial_cmp(&flow.filled.data[a as usize])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let acc_bad = accumulate(&broken, &land, &flow.direction, FW, FH);
    let (in_bad, out_bad, r_bad) = carry(&acc_bad, &flat, &flow.direction, FW, FH);

    // 3. production's GEOMETRIC accumulation — flat-aware since it was written, never broken
    let (in_g, out_g, r_g) = carry(&flow.accumulation.data, &flat, &flow.direction, FW, FH);

    // 4. `runoff_accumulation`, the function the fix changed, on a uniform wet climate
    let t_c = 10.0f32;
    let pe = potential_evaporation_mm(t_c);
    let p_unit = precip_mm_per_year(1.0);
    let precip = GridF32::new(FW, FH, (4.0 * pe.max(1.0)) / p_unit);
    let temperature = GridF32::new(FW, FH, t_c);
    let clim = DrainageClimate { precip_internal: &precip, temperature: &temperature };
    let runoff = runoff_accumulation(&raw, &flow, &clim, 1.0, None, None, FW, FH);
    let (in_r, out_r, r_r) = carry(&runoff, &flat, &flow.direction, FW, FH);

    // 5. MFD — block 4 of the round, read on the same fixture, corrected by nothing
    let mfd = mfd_accumulation(&flow.filled, &flow.direction, SEA, 2.0, FW, FH);
    let (in_m, out_m, r_m) = carry(&mfd.data, &flat, &flow.direction, FW, FH);

    eprintln!("   {:<34} {:>12} {:>12} {:>9}", "reading", "IN", "OUT", "OUT/IN");
    for (name, i, o, r) in [
        ("topological order (the fix)", in_ok, out_ok, r_ok),
        ("`filled`-only order (PRE-FIX)", in_bad, out_bad, r_bad),
        ("flow.accumulation (geometric)", in_g, out_g, r_g),
        ("runoff_accumulation (production)", in_r, out_r, r_r),
        ("mfd_accumulation (incision)", in_m, out_m, r_m),
    ] {
        eprintln!("   {name:<34} {i:>12.2} {o:>12.2} {r:>9.4}");
    }

    // THE GUARD
    assert!(
        r_ok >= 0.95,
        "the topological order must carry the tributary across the flat; read {r_ok:.4}"
    );
    // THE NEGATIVE CONTROL ON THE GUARD: the broken order must FAIL the same threshold, or
    // the guard is not measuring what it claims and would not have caught Finding 73.
    assert!(
        r_bad < 0.95,
        "the deliberately broken order still passed at {r_bad:.4} — this fixture does not \
         exercise the defect, so the guard above is vacuous (rule 10)"
    );
    assert!(
        r_r >= 0.95,
        "`runoff_accumulation` must carry across a flat after the Finding 73 fix; read {r_r:.4}"
    );
    // Production's geometric accumulation was ALREADY flat-aware; pinned so a future
    // "simplification" of its tiebreak fails here rather than twenty findings later.
    assert!(r_g >= 0.95, "flow.accumulation regressed on flats: {r_g:.4}");
    eprintln!(
        "   ⇒ guard green: the fix carries ({r_ok:.4}), the pre-fix order does not \
         ({r_bad:.4}) — a ×{:.0} separation, so the guard is not vacuous",
        r_ok / r_bad.max(1e-9)
    );
    eprintln!("   ⇒ MFD (block 4, DIAGNOSIS ONLY, corrected by nothing): OUT/IN {r_m:.4}");
}

use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::upscale_from_c1_with_progress;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::upscale::{ProductionHdOpts, production_hd_config};

// ── ADR Finding 75 block C · the MFD carrier's NEGATIVE CONTROL ─────────────────────────
//
// The round asked for a bench-only bisection of the MFD carrier. **It cannot be built.**
// `mfd_accumulation` is called at `stream_power.rs:403` INSIDE the incision loop, with no
// parameter and no config field to substitute it, so a bench cannot re-incise with a corrected
// carrier without a production seam — and this round's budget is two production changes, both
// spent (the note, the inventory floor). The seam is specified in the ADR, not opened here.
//
// What CAN be built, and is what the round said must come first: the control that the toggle
// would actually move the reading. Same fixture as Finding 74, same metric — a bench-local copy
// of `mfd_accumulation` whose ONLY difference is the propagation order.

/// `mfd_accumulation` (`flow.rs:866`) with the order corrected and NOTHING else touched.
/// Production sorts `land` by `(filled desc, INDEX desc)` — an index tiebreak, which on a flat
/// is arbitrary. This walks `propagation_order` instead. The weighting, the `cnt == 0` D8
/// fallback and the seeds are copied verbatim.
fn mfd_accumulation_ordered(
    filled: &GridF32,
    direction: &[u8],
    sea_level: f32,
    p: f32,
    w: usize,
    h: usize,
) -> Vec<f32> {
    let n = w * h;
    let is_land: Vec<bool> = (0..n).map(|i| filled.data[i] > sea_level).collect();
    let (order, n_cyclic) = propagation_order(direction, |k| is_land[k], w, h);
    assert_eq!(n_cyclic, 0, "the fixture's D8 graph must be acyclic");
    let mut acc = vec![0.0f32; n];
    for k in 0..n {
        if is_land[k] {
            acc[k] = 1.0;
        }
    }
    let nbr = |i: usize, j: usize, d: usize| -> usize {
        let ni = ((i as i32 + D8_DX[d]).rem_euclid(w as i32)) as usize;
        let nj = ((j as i32 + D8_DY[d]).rem_euclid(h as i32)) as usize;
        nj * w + ni
    };
    for &kk in &order {
        let c = kk as usize;
        let (ci, cj) = (c % w, c / w);
        let zc = filled.data[c];
        let (mut wsum, mut cnt) = (0.0f32, 0usize);
        let (mut wj, mut nj_idx) = ([0.0f32; 8], [0usize; 8]);
        for d in 0..8 {
            let m = nbr(ci, cj, d);
            let drop = zc - filled.data[m];
            if drop > 0.0 {
                let slope = drop / ymir_core::terrain::flow::D8_DIST[d];
                let wgt = slope.powf(p);
                wj[cnt] = wgt;
                nj_idx[cnt] = m;
                wsum += wgt;
                cnt += 1;
            }
        }
        let flow = acc[c];
        if cnt == 0 || wsum <= 0.0 {
            let d = direction[c];
            if d != DIR_NONE {
                acc[nbr(ci, cj, d as usize)] += flow;
            }
            continue;
        }
        for k in 0..cnt {
            acc[nj_idx[k]] += flow * wj[k] / wsum;
        }
    }
    acc
}

#[test]
fn the_mfd_carrier_toggle_moves_the_reading() {
    let raw = fixture();
    let flow = compute_flow(&raw, &FlowConfig { sea_level: SEA, ..Default::default() });
    let flat = flat_of(&raw, &flow.filled);
    let shipped = mfd_accumulation(&flow.filled, &flow.direction, SEA, 2.0, FW, FH);
    let fixed = mfd_accumulation_ordered(&flow.filled, &flow.direction, SEA, 2.0, FW, FH);
    let (i0, o0, r0) = carry(&shipped.data, &flat, &flow.direction, FW, FH);
    let (i1, o1, r1) = carry(&fixed, &flat, &flow.direction, FW, FH);
    eprintln!(
        "
── block C · the MFD carrier, negative control ──
   {:<30} {:>12.2} {:>12.2}          {:>9.4}
   {:<30} {:>12.2} {:>12.2} {:>9.4}",
        "mfd_accumulation (SHIPPED)", i0, o0, r0, "same, order corrected", i1, o1, r1
    );
    assert!(r0 < 0.5, "the shipped carrier must under-read across the flat; got {r0:.4}");
    assert!(
        r1 >= 0.95,
        "the toggle must LIFT the reading to ~1 before it can serve a bisection; got {r1:.4}"
    );
    eprintln!(
        "   ⇒ the toggle moves the fixture {r0:.4} → {r1:.4} (×{:.1}). It is fit to serve a          bisection; the bisection itself needs a production seam at `stream_power.rs:403`.",
        r1 / r0.max(1e-9)
    );
}

// ── the 8192² production-path guard ─────────────────────────────────────────────────────
//
// The Finding 73 instrument, PROMOTED. It runs on `assemble_hd_drainage` — production's one
// carry — at the thresholds declared in Finding 73 BEFORE the population was known: median
// OUT/IN >= 1.00 and >= 80 % of the population at >= 0.95. `#[ignore]` because it costs a
// full 8192² build; the cheap synthetic guard above is the one that runs on every commit,
// and that split is deliberate (a four-minute test does not protect a default suite).

const PSEED: u64 = 10_481_999_410_520_546_993;
const DOMAIN_KM: f32 = 400.0;
const TARGET: usize = 8192;

fn terrain() -> GridF32 {
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
    let cfg = production_hd_config(&ProductionHdOpts {
        target_size: TARGET,
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

#[test]
#[ignore]
fn forty_two_lake_carry_on_the_production_path() {
    use std::collections::{HashMap, HashSet};
    use ymir_core::climate::c1_climate_placed;
    use ymir_core::climate::precipitation::PrecipParams;
    use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, LakeType};
    use ymir_core::tectonics_c1::hd_assembly::assemble_hd_drainage;
    use ymir_core::terrain::flow::breach_monotone;

    let ss = SteinSteinParams::default();
    let raw = terrain();
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let pre0 =
        ymir_core::tectonics_c1::drainage::c1_drainage_windowed(&raw, None, &dcfg, &ss, DOMAIN_KM);
    let field = breach_monotone(
        &raw,
        &pre0.flow.filled,
        &pre0.lake_map,
        C1_SEA_LEVEL_NORM,
        raw.width,
        raw.height,
    );
    let (w, h) = (field.width, field.height);
    // HUMID: the arid world has no exorheic detected surface lake at all (Finding 72), so the
    // population is empty there BY THE WORLD, not by accident. Declared, not discovered.
    let climate = c1_climate_placed(&field, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
    let dclim = DrainageClimate {
        precip_internal: &climate.precipitation,
        temperature: &climate.temperature,
    };
    let pre =
        ymir_core::tectonics_c1::drainage::c1_drainage_windowed(&raw, None, &dcfg, &ss, DOMAIN_KM);
    let dr =
        assemble_hd_drainage(&field, &dclim, Some(pre), &dcfg, &ss, DOMAIN_KM, 1.0, None, false)
            .drainage;

    let mut foot: HashSet<u32> = HashSet::new();
    for &id in &dr.lake_map {
        if id != 0 {
            foot.insert(id);
        }
    }
    let nb_ids = |k: usize| -> Vec<u32> {
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        let mut v = Vec::new();
        for dy in -1i32..=1 {
            for dx in -1i32..=1 {
                let (nx, ny) = (x + dx, y + dy);
                if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                    let id = dr.lake_map[ny as usize * w + nx as usize];
                    if id != 0 && !v.contains(&id) {
                        v.push(id);
                    }
                }
            }
        }
        v
    };
    let (mut out_q, mut in_q): (HashMap<u32, f32>, HashMap<u32, f32>) = Default::default();
    for (i, s) in dr.rivers.segments.iter().enumerate() {
        if dr.segment_kind[i] != SegmentKind::Watercourse {
            continue;
        }
        let q = dr.segment_discharge_m3s[i];
        let (&(sx, sy), &(ex, ey)) = (s.points.first().unwrap(), s.points.last().unwrap());
        let src = nb_ids(sy as usize * w + sx as usize);
        let mouth = nb_ids(ey as usize * w + ex as usize);
        for id in &src {
            if mouth.contains(id) {
                continue;
            }
            let e = out_q.entry(*id).or_insert(0.0);
            *e = e.max(q);
        }
        for id in &mouth {
            if src.contains(id) {
                continue;
            }
            let e = in_q.entry(*id).or_insert(0.0);
            *e = e.max(q);
        }
    }
    let mut spill_fed: HashSet<u32> = HashSet::new();
    for (i, s) in dr.rivers.segments.iter().enumerate() {
        if dr.segment_kind[i] == SegmentKind::Spillway {
            let &(lx, ly) = s.points.last().unwrap();
            let id = dr.lake_map[ly as usize * w + lx as usize];
            if id != 0 {
                spill_fed.insert(id);
            }
        }
    }
    let mut r: Vec<f32> = dr
        .lakes
        .iter()
        .filter(|l| {
            l.base.id < 1_000_001
                && l.lake_type == LakeType::Exorheic
                && foot.contains(&l.base.id)
                && !spill_fed.contains(&l.base.id)
        })
        .filter_map(|l| {
            let qi = in_q.get(&l.base.id).copied().unwrap_or(0.0);
            let qo = out_q.get(&l.base.id).copied()?;
            (qi > 0.0).then_some(qo / qi)
        })
        .collect();
    r.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let np = r.len();
    // Rule 10: an empty or tiny population must FAIL, not pass vacuously.
    assert!(np >= 20, "population collapsed to {np} lakes — the guard would be vacuous");
    let median = r[np / 2];
    let share = r.iter().filter(|&&v| v >= 0.95).count() as f64 / np as f64;
    eprintln!(
        "
── 42-lake carry guard, production path ── n {np} | median          {median:.3} | share >= 0.95 {:.1} % | min {:.3} | max {:.3}",
        100.0 * share,
        r[0],
        r[np - 1]
    );
    assert!(
        median >= 1.00,
        "the clipped network no longer carries flow across an exorheic surface lake: median          OUT/IN {median:.3} over {np} lakes (Finding 73's threshold is 1.00). This is the          accumulation-order regression."
    );
    assert!(share >= 0.80, "only {:.1} % of the population at r >= 0.95", 100.0 * share);
}
