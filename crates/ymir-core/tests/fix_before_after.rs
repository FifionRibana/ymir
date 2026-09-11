//! ADR Finding 74 — the SAME bench run in TWO code states, so every Finding 71 number
//! gets an old -> new pair. Compiles against the PRE-FIX API on purpose (it must run with the
//! fix stashed), so it may not reference `propagation_order`.
//!
//! Every discharge carries its `SegmentKind` **and its population of definition**: local basin
//! only / inherited through the chain / net of evaporation.
//!
//! `assemble_hd_drainage` discards both `BasinSummary` and `Spillway::chained_into`, so this
//! calls `below_sea_basin_lakes_infil` directly after replicating the production sequence up to
//! that point (pre-breach drainage → breach → climate → final drainage → H-1c).
//!
//! Run: cargo test -p ymir-core --release --test fix_before_after -- --ignored --nocapture

use std::collections::{HashMap, HashSet};

use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::{PrecipParams, precip_mm_per_year};
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, LakeType, SegmentKind, SegmentRow, apply_lake_water_balance,
    below_sea_basin_lakes_infil, c1_drainage_windowed, clip_rivers_to_lakes,
    potential_evaporation_mm, runoff_accumulation, runoff_km2_to_m3s,
};
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::upscale_from_c1_with_progress;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::flow::{DIR_NONE, FlowConfig, RiverSegment, breach_monotone, compute_flow};
use ymir_core::terrain::upscale::{ProductionHdOpts, production_hd_config};

const PSEED: u64 = 10_481_999_410_520_546_993;
const DOMAIN_KM: f32 = 400.0;
const SEA: f32 = 0.5;
const TARGET: usize = 8192;
const CELL_KM: f32 = DOMAIN_KM / TARGET as f32;
const BEDS: [(&str, f32, f32); 2] = [("humid", 45.0, 40.0), ("arid-hot", 25.0, 10.0)];
const REFERENCE_RUNOFF_MM: f32 = 300.0;

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

/// The six terminus classes, in priority order. `did` = detected surface lake id.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Terminus {
    Ocean,
    BasinWithSpillway,
    BasinWithoutSpillway,
    SurfaceLake,
    AbsorbedReceiver,
    DryLand,
}

fn cls_name(t: Terminus) -> &'static str {
    match t {
        Terminus::Ocean => "ocean (wc=1)",
        Terminus::BasinWithSpillway => "basin WITH spillway",
        Terminus::BasinWithoutSpillway => "basin WITHOUT spillway",
        Terminus::SurfaceLake => "SURFACE LAKE (did≠0)",
        Terminus::AbsorbedReceiver => "ABSORBED receiver (dead id)",
        Terminus::DryLand => "true DRY LAND",
    }
}

#[test]
#[ignore]
fn fix_before_after() {
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  Finding 72 · the ledger of the below-sea basins  ==========");
    let raw = terrain();
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let pre0 = c1_drainage_windowed(&raw, None, &dcfg, &ss, DOMAIN_KM);
    let field =
        breach_monotone(&raw, &pre0.flow.filled, &pre0.lake_map, SEA, raw.width, raw.height);
    let (w, h) = (field.width, field.height);
    let n = w * h;
    let cell_km2 = CELL_KM * CELL_KM;
    let flow = compute_flow(&field, &FlowConfig { sea_level: SEA, ..Default::default() });
    let wc = water_class(&field, SEA);

    // ── CONTROL BLOCK · the height field must NOT move ────────────────────────────────
    // Rule 7. `runoff_accumulation` is consumed by the segment discharges, the lake water
    // balance and the below-sea basins — NOT by the erosion (which uses `mfd_accumulation`)
    // and NOT by the breach (topological). So the fix must be BIT-IDENTICAL here. Printed as
    // a fingerprint so the two code states can be compared digit by digit; a stop rule in the
    // round, not an assert, because the pre-fix run must also print it.
    let fingerprint = |g: &GridF32, name: &str| {
        let mut hash = 0xcbf2_9ce4_8422_2325u64; // FNV-1a over the raw f32 bits
        let mut sum = 0.0f64;
        let (mut lo, mut hi) = (f32::MAX, f32::MIN);
        let mut land = 0usize;
        for &v in &g.data {
            hash ^= v.to_bits() as u64;
            hash = hash.wrapping_mul(0x100_0000_01b3);
            sum += v as f64;
            lo = lo.min(v);
            hi = hi.max(v);
            if v > SEA {
                land += 1;
            }
        }
        let mut s: Vec<f32> = g.data.iter().copied().filter(|&v| v > SEA).collect();
        s.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        let q = |p: f64| s[((s.len() - 1) as f64 * p) as usize];
        eprintln!(
            "   {name:<10} hash {hash:#018x} | mean {:.9} | min {lo:.6} max {hi:.6} | land \
             {land} | land p10 {:.6} p50 {:.6} p90 {:.6}",
            sum / g.data.len() as f64,
            q(0.10),
            q(0.50),
            q(0.90)
        );
    };
    eprintln!("\n── CONTROL BLOCK · field fingerprints (must be IDENTICAL across the fix) ──");
    fingerprint(&raw, "eroded");
    fingerprint(&field, "breached");
    fingerprint(&pre0.flow.filled, "pre-filled");
    eprintln!(
        "   depressions vs the pre-incision base: {} cells raised by `pit_fill` | below-sea \
         (wc=2) {} cells | ocean (wc=1) {} cells",
        (0..n).filter(|&k| pre0.flow.filled.data[k] > raw.data[k] + 1e-7).count(),
        (0..n).filter(|&k| wc[k] == 2).count(),
        (0..n).filter(|&k| wc[k] == 1).count()
    );

    // A3-four is defined on the HUMID bed (the arid world has no exorheic surface lake at all),
    // so its verdict governs BOTH beds: the arid ledger cannot be more readable than the
    // instrument that would read it. Humid runs first.
    let mut a3_validated = false;
    for (bed, lat, span) in BEDS {
        eprintln!("\n╔═══════════════ {bed} ({lat}°/span {span}) ═══════════════╗");
        let climate =
            c1_climate_placed(&field, &ss, lat, span, &PrecipParams::default(), DOMAIN_KM);
        let dclim = DrainageClimate {
            precip_internal: &climate.precipitation,
            temperature: &climate.temperature,
        };
        // per-cell runoff surplus, mm·km²/yr
        let surplus_mm_km2: Vec<f32> = (0..n)
            .map(|k| {
                let p = precip_mm_per_year(climate.precipitation.data[k]);
                let pe = potential_evaporation_mm(climate.temperature.data[k]);
                (p - pe).max(0.0) * cell_km2
            })
            .collect();

        // ── PERIMETER, both biases, once ───────────────────────────────────────────
        let mut b_land = 0.0f64; // the A0 budget population: field > SEA
        let mut b_lakefoot = 0.0f64; // of which detected-lake footprints
        let mut b_belowsea = 0.0f64; // rain surplus over enclosed below-sea (wc==2), OUT of A0
        let (mut n_land, mut n_lf, mut n_bs) = (0usize, 0usize, 0usize);

        // production sequence up to the below-sea call
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
        for k in 0..n {
            if field.data[k] > SEA {
                n_land += 1;
                b_land += surplus_mm_km2[k] as f64;
                if dr.lake_map[k] != 0 && dr.lake_map[k] < 1_000_001 {
                    n_lf += 1;
                    b_lakefoot += surplus_mm_km2[k] as f64;
                }
            } else if wc[k] == 2 {
                n_bs += 1;
                b_belowsea += surplus_mm_km2[k] as f64;
            }
        }
        let budget = runoff_km2_to_m3s(b_land as f32) as f64;
        eprintln!(
            "\n── PERIMETER ──\n   A0 budget population: {n_land} land cells (field > SEA) ⇒ \
             **{budget:.1} m³/s**\n   of which DETECTED-LAKE footprints: {n_lf} cells = {:.1} \
             m³/s ({:.1} % of budget) — IN the budget (the breach removed the depressions), and \
             open-water PE is NOT deducted ⇒ A0 OVERSTATES by at most this\n   rain surplus over \
             ENCLOSED BELOW-SEA footprints (wc=2): {n_bs} cells = {:.1} m³/s — OUT of the budget",
            runoff_km2_to_m3s(b_lakefoot as f32),
            100.0 * runoff_km2_to_m3s(b_lakefoot as f32) as f64 / budget,
            runoff_km2_to_m3s(b_belowsea as f32)
        );
        eprintln!(
            "   and it is ALSO out of the basins' inflow: `runoff_accumulation` only seeds cells\n   \
             with `heightmap > C1_SEA_LEVEL_NORM`, so a below-sea cell contributes 0. The two\n   \
             perimeter biases are therefore {:.1} m³/s (A0 too high) and {:.1} m³/s (never \
             counted anywhere) — bounded, and both small against the hole.",
            runoff_km2_to_m3s(b_lakefoot as f32),
            runoff_km2_to_m3s(b_belowsea as f32)
        );

        // The map the SPILLWAY TRACE saw: production passes `Some(&drainage.lake_map)` BEFORE the
        // below-sea merge, so `did` inside the trace (`drainage.rs:1631`) is read from this one.
        let premerge_lake_map = dr.lake_map.clone();
        let bs = below_sea_basin_lakes_infil(
            &field,
            &dclim,
            &dcfg,
            &ss,
            DOMAIN_KM,
            Some(&premerge_lake_map),
            None,
        );
        // Same call, same arguments, as inside `below_sea_basin_lakes_infil` (sinks None,
        // infiltration None) ⇒ this array is the one the basin inflows were summed from.
        let runoff = runoff_accumulation(&field, &flow, &dclim, cell_km2, None, None, w, h);

        // ── production's tail, replicated so the network is the CLIPPED one ────────
        // `assemble_hd_drainage` lines 144–235: merge the below-sea water over the detected
        // map, drop the detected lakes it submerged, THEN clip, THEN append the spillways.
        // Replicated here (not called) because `assemble_hd_drainage` discards `BasinSummary`
        // and `chained_into` — the two fields the ledger is made of (the computed-not-exported
        // family, Finding 72). `geo_scale_ratio` stays 1.0, the code's default: identity.
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
        dr.lakes.retain(|lk| lk.base.id >= 1_000_001 || !absorbed.contains(&lk.base.id));
        dr.lakes.extend(bs.lakes.iter().cloned());
        let before_seg = dr.rivers.segments.len();
        clip_rivers_to_lakes(&mut dr);
        let after_clip = dr.rivers.segments.len();
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
        assert!(dr.segment_arrays_aligned(), "the replicated tail desynchronised an array");
        eprintln!(
            "   network: {before_seg} pre-clip → {after_clip} clipped → {} with {} spillways \
             appended | {} detected lake(s) absorbed by the merge",
            dr.rivers.segments.len(),
            bs.spillways.len(),
            absorbed.len()
        );

        // footprints on the FINAL (merged) map — the one the clip ran against
        let mut foot: HashMap<u32, Vec<usize>> = HashMap::new();
        for k in 0..n {
            let id = dr.lake_map[k];
            if id != 0 {
                foot.entry(id).or_default().push(k);
            }
        }
        // 8-neighbour lake ids of a cell, on the final map
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
        // Per lake: max discharge of a WATERCOURSE run whose SOURCE borders it (an outlet, the
        // Finding 37 definition) and of one whose MOUTH borders it (an inlet). A run that is
        // both is counted in NEITHER — declared, and reported.
        let mut out_q: HashMap<u32, f32> = HashMap::new();
        let mut in_q: HashMap<u32, f32> = HashMap::new();
        let (mut n_out, mut n_in, mut n_both) = (0usize, 0usize, 0usize);
        for (i, s) in dr.rivers.segments.iter().enumerate() {
            if dr.segment_kind[i] != SegmentKind::Watercourse {
                continue;
            }
            let q = dr.segment_discharge_m3s[i];
            let (&(sx, sy), &(ex, ey)) = (s.points.first().unwrap(), s.points.last().unwrap());
            let src = nb_ids(sy as usize * w + sx as usize);
            let mouth = nb_ids(ey as usize * w + ex as usize);
            let both: Vec<u32> = src.iter().copied().filter(|id| mouth.contains(id)).collect();
            for id in &src {
                if both.contains(id) {
                    continue;
                }
                n_out += 1;
                let e = out_q.entry(*id).or_insert(0.0);
                *e = e.max(q);
            }
            for id in &mouth {
                if both.contains(id) {
                    continue;
                }
                n_in += 1;
                let e = in_q.entry(*id).or_insert(0.0);
                *e = e.max(q);
            }
            n_both += both.len();
        }
        // lakes a spillway lands in (either map: the trace read the pre-merge one)
        let mut spill_fed: HashSet<u32> = HashSet::new();
        for sw in &bs.spillways {
            let &(lx, ly) = sw.points.last().unwrap();
            let k = ly as usize * w + lx as usize;
            for id in [dr.lake_map[k], premerge_lake_map[k]] {
                if id != 0 {
                    spill_fed.insert(id);
                }
            }
        }

        // ── A3-FOUR · BLOCKING CONTROL, on the SEGMENT network ─────────────────────
        //
        // FOURTH construction, to the specification written in Finding 72 after three of my own
        // failures: discharge of the run LEAVING the lake against the max discharge of the runs
        // ENTERING it — clipped network, `SegmentKind::Watercourse`, population INHERITED.
        //
        // The three that failed, so this is not repeated:
        //   1. one step downstream of `lk.base.outlet` — the PRE-BREACH sill, read on the
        //      BREACHED field where that depression is gone. A saddle carries no accumulation.
        //   2. the footprint's own max-accumulation cell, walked downstream — lands at or below
        //      sea level, where `runoff_accumulation` leaves zero BY CONSTRUCTION.
        //   3. either of those on the arid bed — EMPTY POPULATION, a rule-10 failure of mine.
        //
        // Why the segment network instead: `clip_rivers_to_lakes` (`drainage.rs:756`) makes an
        // outlet run INHERIT `src_q[i]`, the parent's value at the PARENT'S MOUTH, downstream
        // of the lake. That value is only large if the accumulation actually crossed the lake,
        // so the ratio does test the carry — but its MAGNITUDE above 1 is not interpretable
        // (the parent gains catchment below the lake). One-sided test, declared as such.
        //
        // HUMID ONLY, declared in advance: attempt 3 established that the arid bed has no
        // The Finding 73 carry metric, reported BEFORE and AFTER the fix. No stop rule
        // here: this bench exists to produce the old -> new pairs, and Finding 71 already
        // published the old numbers. The ASSERTING guard lives in `accumulation_order.rs`.
        if bed == "humid" {
            eprintln!(
                "\n── A3-four · BLOCKING CONTROL — does the CLIPPED NETWORK carry flow across an \
                 exorheic surface lake? ──"
            );
            eprintln!(
                "   endpoint bookkeeping: {n_out} (run source, lake) pairs · {n_in} (run mouth, \
                 lake) pairs · {n_both} pairs excluded because the SAME run both starts and ends \
                 on the same lake"
            );
            // VALIDATION population: exorheic detected lakes, NOT fed by a spillway, with both
            // an outlet run and an inlet run carrying non-zero discharge.
            let mut rows: Vec<(u32, f32, f32, f32, f32)> = Vec::new(); // id, area, in, out, r
            let mut excluded_spill = 0usize;
            let mut excluded_noinlet = 0usize;
            let mut excluded_nooutlet = 0usize;
            for lk in dr.lakes.iter() {
                if lk.base.id >= 1_000_001 || lk.lake_type != LakeType::Exorheic {
                    continue;
                }
                if !foot.contains_key(&lk.base.id) {
                    continue;
                }
                if spill_fed.contains(&lk.base.id) {
                    excluded_spill += 1;
                    continue;
                }
                let qi = in_q.get(&lk.base.id).copied().unwrap_or(0.0);
                let qo = out_q.get(&lk.base.id).copied().unwrap_or(0.0);
                if qi <= 0.0 {
                    excluded_noinlet += 1;
                    continue;
                }
                if out_q.get(&lk.base.id).is_none() {
                    excluded_nooutlet += 1;
                    continue;
                }
                rows.push((lk.base.id, lk.area_km2, qi, qo, qo / qi));
            }
            eprintln!(
                "   POPULATION {} lakes  (excluded: {excluded_spill} spillway-fed, \
                 {excluded_noinlet} with no inlet run, {excluded_nooutlet} with no outlet run)",
                rows.len()
            );
            if rows.len() < 5 {
                eprintln!(
                    "   ⚠️ POPULATION UNDER 5 — rule 10 on my own control: this cannot be read."
                );
            }
            rows.sort_by(|a, b| a.4.partial_cmp(&b.4).unwrap_or(std::cmp::Ordering::Equal));
            let np = rows.len();
            let median = if np == 0 { 0.0 } else { rows[np / 2].4 };
            let share_ok =
                rows.iter().filter(|r| r.4 >= 0.95).count() as f64 / np.max(1) as f64 * 100.0;
            let pct = |p: f64| -> f32 {
                if np == 0 { 0.0 } else { rows[((np as f64 - 1.0) * p) as usize].4 }
            };
            eprintln!(
                "   ratio OUT/IN  [Watercourse, clipped, discharge population INHERITED at the \
                 outlet run]\n   p10 {:.3} · p25 {:.3} · **median {median:.3}** · p75 {:.3} · p90 \
                 {:.3} · min {:.3} · max {:.3}\n   share with r ≥ 0.95: **{share_ok:.1} %**  \
                 (declared threshold: median ≥ 1.00 AND share ≥ 80 %)",
                pct(0.10),
                pct(0.25),
                pct(0.75),
                pct(0.90),
                rows.first().map(|r| r.4).unwrap_or(0.0),
                rows.last().map(|r| r.4).unwrap_or(0.0)
            );
            eprintln!(
                "   the 6 LOWEST ratios (where a failure would show) and the 3 largest lakes:\n   \
                 {:>8} {:>10} {:>12} {:>12} {:>8}",
                "lake", "area km²", "max IN m³/s", "OUT m³/s", "r"
            );
            for r in rows.iter().take(6) {
                eprintln!("   {:>8} {:>10.2} {:>12.4} {:>12.4} {:>8.3}", r.0, r.1, r.2, r.3, r.4);
            }
            let mut by_area = rows.clone();
            by_area.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
            for r in by_area.iter().take(3) {
                eprintln!(
                    "   {:>8} {:>10.2} {:>12.4} {:>12.4} {:>8.3}  ← by area",
                    r.0, r.1, r.2, r.3, r.4
                );
            }
            a3_validated = np >= 5 && median >= 1.00 && share_ok >= 80.0;

            // SENSITIVITY — the same instrument on an interface Finding 71 knows is unbalanced:
            // a below-sea basin and its spillway. An instrument reading ≈ 1 everywhere measures
            // nothing, so this must read the gap.
            eprintln!(
                "\n   SENSITIVITY — the same ratio at a BELOW-SEA basin / spillway interface \
                 (F71 knows these unbalanced):\n   {:>10} {:>14} {:>16} {:>10} {:>12}",
                "basin", "max IN m³/s", "spillway OUT m³/s", "r", "chained_into"
            );
            let mut sens: Vec<(u32, f32, f32, f32)> = bs
                .spillways
                .iter()
                .map(|sw| {
                    let qi = in_q.get(&sw.lake_id).copied().unwrap_or(0.0);
                    (sw.lake_id, qi, sw.discharge_m3s, sw.discharge_m3s / qi.max(1e-9))
                })
                .collect();
            sens.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
            for s in sens.iter().take(5) {
                let ci = bs
                    .spillways
                    .iter()
                    .find(|sw| sw.lake_id == s.0)
                    .and_then(|sw| sw.chained_into)
                    .map_or("None".to_string(), |v| v.to_string());
                eprintln!(
                    "   {:>10} {:>14.4} {:>16.2} {:>10.1} {:>12}",
                    s.0,
                    s.1,
                    s.2,
                    if s.1 > 0.0 { s.3 } else { f32::INFINITY },
                    ci
                );
            }
            let with_inlet: Vec<f32> =
                sens.iter().filter(|s| s.1 > 0.0).map(|s| s.3).collect::<Vec<_>>();
            let mut wi = with_inlet.clone();
            wi.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            eprintln!(
                "   {} of {} spillway basins have an identifiable inlet run; median r there \
                 **{:.1}** against **{median:.3}** at an ordinary surface lake ⇒ dynamic range \
                 ×{:.0}",
                wi.len(),
                sens.len(),
                wi.get(wi.len() / 2).copied().unwrap_or(0.0),
                wi.get(wi.len() / 2).copied().unwrap_or(0.0) / median.max(1e-9)
            );

            eprintln!(
                "\n   ⇒ {}",
                if a3_validated {
                    "VALIDATES: the clipped network DOES carry flow across an exorheic surface \
                     lake, and the instrument has range. The ledger is readable."
                } else {
                    "DOES NOT VALIDATE"
                }
            );
        }

        // The chain edges, needed by BOTH the cross-check (block 2) and the ledger:
        // Σ discharge of the spillways that NAME a basin in `chained_into`.
        let mut in_chain: HashMap<u32, f64> = HashMap::new();
        for sw in &bs.spillways {
            if let Some(t) = sw.chained_into {
                *in_chain.entry(t).or_default() += sw.discharge_m3s as f64;
            }
        }
        let out_of: HashMap<u32, f64> =
            bs.spillways.iter().map(|s| (s.lake_id, s.discharge_m3s as f64)).collect();
        // ── §3 · BLOCKING · ONE basin reconstructed end to end, independently ─────
        //
        // Target: 1000056, the 240.8 m³/s basin of Finding 71. Admissible error DECLARED BEFORE
        // MEASURING: ±10 % passes, ±10–25 % is marginal and reported as such, **> ±25 % fires
        // stop rule 2** — the ledger cannot then be read from `BasinSummary`.
        //
        // Three routes, kept separate so a gap can be ATTRIBUTED instead of averaged:
        //  (a) the watershed × the measured runoff surplus — geometric, integrated straight from
        //      precipitation and PE, and INDEPENDENT of `runoff_accumulation`;
        //  (b) Σ discharge of the spillways naming it in `chained_into` — the export's chain;
        //  (c) the shoreline sum re-derived exactly as `drainage.rs:1478` does it, from the
        //      bench's own `runoff` array (same call, same arguments ⇒ same array). NOT
        //      independent: it isolates `extra_inflow` = inflow_m3s − (c) so (a) and (b) can be
        //      tested against the right halves.
        eprintln!("\n── §3 · BLOCKING CROSS-CHECK — basin 1000056 reconstructed end to end ──");
        let targets: Vec<u32> = if bs.basins.iter().any(|b| b.id == 1_000_056) {
            vec![1_000_056]
        } else {
            let v: Vec<u32> = bs
                .basins
                .iter()
                .max_by(|a, b| {
                    a.inflow_m3s.partial_cmp(&b.inflow_m3s).unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|b| vec![b.id])
                .unwrap_or_default();
            eprintln!(
                "   ⚠️ basin 1000056 is not present on this bed; falling back to the largest \
                 `inflow_m3s` basin {v:?} — a DIFFERENT object, stated as such"
            );
            v
        };
        let mut xcheck_ok = true;
        for target in targets {
            // (a) the watershed: the D8 basin labels of the SHORELINE ENTRY CELLS — the
            // above-sea cells whose receiver is a footprint cell — then every above-sea cell
            // carrying one of those labels.
            //
            // ⚠️ CORRECTED, and it was my bug, not the machinery's. The first version read
            // `flow.basins` ON the footprint cells. `compute_flow` treats EVERY cell at or
            // below the sea level as `is_ocean`, and `compute_basins` skips those, so a
            // below-sea cell carries label 0 — the collection came back EMPTY and route (a)
            // read 0.00 m³/s. That would have been reported as "the machinery's accounting is
            // not confirmed" when the only thing not working was my label lookup. Fourth
            // instrument failure of this campaign; the pattern is always the same, a raster
            // whose zeros are load-bearing.
            let labels: HashSet<u32> = (0..n)
                .filter(|&k| field.data[k] > SEA)
                .filter(|&k| {
                    let d = flow.direction[k];
                    if d == DIR_NONE {
                        return false;
                    }
                    let (x, y) = ((k % w) as i32, (k / w) as i32);
                    let nx = (x + ymir_core::terrain::flow::D8_DX[d as usize]).rem_euclid(w as i32);
                    let ny = (y + ymir_core::terrain::flow::D8_DY[d as usize]).rem_euclid(h as i32);
                    bs.lake_map[ny as usize * w + nx as usize] == target
                })
                .map(|k| flow.basins[k])
                .filter(|&l| l != 0)
                .collect();
            let mut recon = 0.0f64;
            let mut recon_cells = 0usize;
            for k in 0..n {
                if field.data[k] > SEA && labels.contains(&flow.basins[k]) {
                    recon += surplus_mm_km2[k] as f64;
                    recon_cells += 1;
                }
            }
            let recon_m3s = runoff_km2_to_m3s(recon as f32) as f64;
            // (c) the shoreline sum, the machinery's own definition, re-derived
            let mut shore = 0.0f64;
            for k in 0..n {
                if bs.lake_map[k] != target || wc[k] != 2 {
                    continue;
                }
                let (x, y) = ((k % w) as i32, (k / w) as i32);
                for (dx, dy) in
                    ymir_core::terrain::flow::D8_DX.iter().zip(ymir_core::terrain::flow::D8_DY)
                {
                    let (nx, ny) = (x + dx, y + dy);
                    if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                        continue;
                    }
                    let nk = ny as usize * w + nx as usize;
                    if field.data[nk] <= SEA {
                        continue;
                    }
                    let d = flow.direction[nk];
                    if d == DIR_NONE {
                        continue;
                    }
                    let (tx, ty) = (
                        nx + ymir_core::terrain::flow::D8_DX[d as usize],
                        ny + ymir_core::terrain::flow::D8_DY[d as usize],
                    );
                    if tx >= 0 && ty >= 0 && (tx as usize) < w && (ty as usize) < h {
                        if ty as usize * w + tx as usize == k {
                            shore += runoff[nk] as f64;
                        }
                    }
                }
            }
            let shore_m3s = runoff_km2_to_m3s(shore as f32) as f64;
            let inc = in_chain.get(&target).copied().unwrap_or(0.0);
            let reported = bs
                .basins
                .iter()
                .find(|b| b.id == target)
                .map(|b| b.inflow_m3s as f64)
                .unwrap_or(0.0);
            let extra = reported - shore_m3s; // what the fixpoint added, isolated
            let total = recon_m3s + inc;
            let r = total / reported.max(1e-9);
            eprintln!(
                "   basin {target}   [Spillway source · population: local + chain-inherited]\n   \
                 (a) watershed {recon_cells} above-sea cells over {} D8 label(s) × measured \
                 surplus = **{recon_m3s:.2} m³/s**  ← independent of `runoff_accumulation`\n   \
                 (b) Σ spillways naming it in `chained_into`         = **{inc:.2} m³/s**\n   \
                 (a)+(b)                                            = **{total:.2} m³/s**\n   \
                 `BasinSummary::inflow_m3s`                          = **{reported:.2} m³/s**  ⇒ \
                 ratio (a+b)/reported = **{r:.3}**\n   (c) shoreline sum re-derived from the same \
                 `runoff` array = {shore_m3s:.2} ⇒ `extra_inflow` isolated = **{extra:.2} m³/s** \
                 ({:.1} % of the inflow)",
                labels.len(),
                100.0 * extra / reported.max(1e-9)
            );
            eprintln!(
                "       (a) vs (c) — the same local inflow by two routes: ratio **{:.3}**  \
                 {}\n       (b) vs the isolated `extra_inflow`: ratio **{:.3}**  {}",
                recon_m3s / shore_m3s.max(1e-9),
                if (recon_m3s / shore_m3s.max(1e-9) - 1.0).abs() <= 0.10 {
                    "— within the declared ±10 %"
                } else {
                    "— OUTSIDE ±10 %"
                },
                inc / extra.max(1e-9),
                if (inc / extra.max(1e-9) - 1.0).abs() <= 0.10 {
                    "— within ±10 %: `chained_into` and `chained_region` agree here"
                } else {
                    "— OUTSIDE ±10 %: the DISPLAY id and the ROUTING key disagree (they are \
                     different fields: `chained_into` a lake id, `chained_region` a region label)"
                }
            );
            let dev = (r - 1.0).abs();
            eprintln!(
                "   ⇒ {}",
                if dev <= 0.10 {
                    "PASSES — within the declared ±10 %. The ledger can be read from `BasinSummary`."
                } else if dev <= 0.25 {
                    "MARGINAL — inside ±25 %, outside ±10 %. Read the ledger, carry the reservation."
                } else {
                    "DIVERGES beyond ±25 %"
                }
            );
            if dev > 0.25 {
                xcheck_ok = false;
            }
        }

        // ── terminus classification of every spillway ──────────────────────────────
        let live_bs: HashSet<u32> = bs.lakes.iter().map(|l| l.base.id).collect();
        let has_spill: HashSet<u32> = bs.spillways.iter().map(|s| s.lake_id).collect();
        let classify =
            |sw: &ymir_core::tectonics_c1::drainage::Spillway| -> (Terminus, u32, u32, u8) {
                let &(lx, ly) = sw.points.last().unwrap();
                let k = ly as usize * w + lx as usize;
                // `did_trace` is what the TRACE saw (`drainage.rs:1631` reads the PRE-MERGE
                // map); `did_final` is what survives in the export. The two differ exactly on
                // the absorbed class, which is why it is testable at all.
                let did_trace =
                    if premerge_lake_map[k] < 1_000_001 { premerge_lake_map[k] } else { 0 };
                let did_final = if dr.lake_map[k] < 1_000_001 { dr.lake_map[k] } else { 0 };
                let bid = bs.lake_map[k];
                let cls = if wc[k] == 1 {
                    Terminus::Ocean
                } else if did_final != 0 && inventoried.contains(&did_final) {
                    Terminus::SurfaceLake
                } else if bid != 0 && bid != sw.lake_id {
                    if has_spill.contains(&bid) {
                        Terminus::BasinWithSpillway
                    } else {
                        Terminus::BasinWithoutSpillway
                    }
                } else if did_trace != 0 && !inventoried.contains(&did_trace) {
                    // the trace stopped on a detected lake the below-sea merge then submerged
                    Terminus::AbsorbedReceiver
                } else if let Some(t) = sw.chained_into {
                    if live_bs.contains(&t) || has_spill.contains(&t) {
                        Terminus::BasinWithSpillway
                    } else {
                        Terminus::AbsorbedReceiver
                    }
                } else {
                    Terminus::DryLand
                };
                (cls, did_final.max(did_trace), bid, wc[k])
            };

        // ── 1 · THE BINARY DISCRIMINANT, reported FIRST ───────────────────────────
        let mut sum_surface = 0.0f64;
        let mut per_cls: HashMap<Terminus, (usize, f64)> = HashMap::new();
        for sw in &bs.spillways {
            let (cls, _, _, _) = classify(sw);
            let e = per_cls.entry(cls).or_insert((0, 0.0));
            e.0 += 1;
            e.1 += sw.discharge_m3s as f64;
            if cls == Terminus::SurfaceLake {
                sum_surface += sw.discharge_m3s as f64;
            }
        }
        // THE HOLE, MEASURED ON THIS BED — not borrowed from Finding 71. Terminal = a clipped
        // `Watercourse` run whose mouth is on, or drains directly into, an OCEAN cell (wc == 1)
        // + a spillway ending in the ocean. The F71 figures (82 runs / 46.5 m³/s in humid) are
        // printed beside it as a cross-check of the replicated tail (rule 7).
        let sp_ocean: f64 = per_cls.get(&Terminus::Ocean).map(|e| e.1).unwrap_or(0.0);
        let (mut wc_sea_n, mut wc_sea_q) = (0usize, 0.0f64);
        for (i, s) in dr.rivers.segments.iter().enumerate() {
            if dr.segment_kind[i] != SegmentKind::Watercourse || s.downstream.is_some() {
                continue;
            }
            let &(ex, ey) = s.points.last().unwrap();
            let k = ey as usize * w + ex as usize;
            let mut sea = wc[k] == 1;
            let d = flow.direction[k];
            if !sea && d != DIR_NONE {
                let (nx, ny) = (
                    (ex as i32 + ymir_core::terrain::flow::D8_DX[d as usize]).rem_euclid(w as i32),
                    (ey as i32 + ymir_core::terrain::flow::D8_DY[d as usize]).rem_euclid(h as i32),
                );
                sea = wc[ny as usize * w + nx as usize] == 1;
            }
            if sea {
                wc_sea_n += 1;
                wc_sea_q += dr.segment_discharge_m3s[i] as f64;
            }
        }
        let terminal = wc_sea_q + sp_ocean;
        let hole = budget - terminal;
        eprintln!(
            "\n── the HOLE, measured on THIS bed ──\n   budget {budget:.1} | Watercourse→sea \
             {wc_sea_n} runs = {wc_sea_q:.1} | Spillway→ocean = {sp_ocean:.1} | TERMINAL \
             {terminal:.1} (×{:.3} of budget) ⇒ **HOLE {hole:.1} m³/s**",
            terminal / budget
        );

        // ── the FINDING 71 TABLE, re-derived ──────────────────────────────────────
        let sp_chained: f64 = bs
            .spillways
            .iter()
            .filter(|s| {
                let &(lx, ly) = s.points.last().unwrap();
                wc[ly as usize * w + lx as usize] == 2
            })
            .map(|s| s.discharge_m3s as f64)
            .sum();
        let n_chained = bs
            .spillways
            .iter()
            .filter(|s| {
                let &(lx, ly) = s.points.last().unwrap();
                wc[ly as usize * w + lx as usize] == 2
            })
            .count();
        eprintln!(
            "   F71 classes: Watercourse→sea {wc_sea_n} / {wc_sea_q:.1} | Spillway→OCEAN {} / \
             {sp_ocean:.1} | Spillway→CHAINED (wc=2) {n_chained} / {sp_chained:.1} ({:.1} % of \
             budget)",
            per_cls.get(&Terminus::Ocean).map(|e| e.0).unwrap_or(0),
            100.0 * sp_chained / budget
        );
        // per-kind distribution, the F71 A2 table
        for kind in [SegmentKind::Watercourse, SegmentKind::Spillway] {
            let mut qs: Vec<f32> = (0..dr.rivers.segments.len())
                .filter(|&i| dr.segment_kind[i] == kind)
                .map(|i| dr.segment_discharge_m3s[i])
                .collect();
            if qs.is_empty() {
                continue;
            }
            qs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let q = |p: f64| qs[((qs.len() - 1) as f64 * p) as usize];
            eprintln!(
                "   {:?}: n {} | p50 {:.3} | p90 {:.3} | p99 {:.3} | **max {:.3} m³/s**",
                kind,
                qs.len(),
                q(0.50),
                q(0.90),
                q(0.99),
                qs.last().copied().unwrap_or(0.0)
            );
        }
        // the named top 20, kind beside every discharge (the standing rule)
        let mut top: Vec<usize> = (0..dr.rivers.segments.len()).collect();
        top.sort_by(|&a, &b| {
            dr.segment_discharge_m3s[b]
                .partial_cmp(&dr.segment_discharge_m3s[a])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let n_sp_top20 =
            top.iter().take(20).filter(|&&i| dr.segment_kind[i] == SegmentKind::Spillway).count();
        eprintln!(
            "   the NAMED TOP 20 by discharge — **{n_sp_top20} of 20 are Spillway**:\n   {:>4} \
             {:>12} {:>14} {:>10} {:>8} {:>12}",
            "rank", "Q m³/s", "kind", "width m", "cells", "source_lake"
        );
        for (r, &i) in top.iter().take(20).enumerate() {
            eprintln!(
                "   {:>4} {:>12.3} {:>14} {:>10.2} {:>8} {:>12}",
                r + 1,
                dr.segment_discharge_m3s[i],
                format!("{:?}", dr.segment_kind[i]),
                dr.segment_width_m[i],
                dr.rivers.segments[i].points.len(),
                dr.segment_source_lake[i].map_or("—".to_string(), |v| v.to_string())
            );
        }
        // the lake inventory, so a classification shift is visible
        let (mut n_exo, mut n_endo, mut a_tot) = (0usize, 0usize, 0.0f64);
        for lk in &dr.lakes {
            match lk.lake_type {
                LakeType::Exorheic => n_exo += 1,
                LakeType::Endorheic => n_endo += 1,
                _ => {}
            }
            a_tot += lk.area_km2 as f64;
        }
        let mut lvl: Vec<(u32, f32, f32)> =
            dr.lakes.iter().map(|l| (l.base.id, l.level_m, l.area_km2)).collect();
        lvl.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
        eprintln!(
            "   LAKES: {} total ({n_exo} exorheic, {n_endo} endorheic), Σ area {a_tot:.1} km² | \
             the 5 largest (id, level_m, area_km²): {:?}",
            dr.lakes.len(),
            lvl.iter()
                .take(5)
                .map(|e| (e.0, (e.1 * 10.0).round() / 10.0, (e.2 * 10.0).round() / 10.0))
                .collect::<Vec<_>>()
        );
        eprintln!(
            "\n── 1 · THE BINARY DISCRIMINANT, reported before any table ──\n   Σ discharge of \
             spillways ending on a SURFACE LAKE (did ≠ 0) = **{sum_surface:.1} m³/s**   \
             [Spillway, population: net of evaporation, local + chain-inherited]\n   against the \
             MEASURED hole {:.1} m³/s ⇒ **coverage ×{:.3}** {}",
            hole,
            sum_surface / hole.max(1e-9),
            if (sum_surface / hole.max(1e-9)) >= 1.0 / 1.5 && (sum_surface / hole.max(1e-9)) <= 1.5
            {
                "— WITHIN a factor 1.5: SINGLE MECHANISM CONFIRMED"
            } else {
                "— OUTSIDE a factor 1.5: real but PARTIAL, the complement must be named"
            }
        );

        eprintln!("\n── 2 · TERMINUS TYPOLOGY, six classes ──");
        eprintln!("   {:<30} {:>6} {:>14}", "class", "n", "Σ Q m³/s");
        for t in [
            Terminus::Ocean,
            Terminus::BasinWithSpillway,
            Terminus::BasinWithoutSpillway,
            Terminus::SurfaceLake,
            Terminus::AbsorbedReceiver,
            Terminus::DryLand,
        ] {
            let (c, q) = per_cls.get(&t).copied().unwrap_or((0, 0.0));
            eprintln!("   {:<30} {c:>6} {q:>14.2}", cls_name(t));
        }
        eprintln!(
            "   (spillway-to-ocean Σ {sp_ocean:.1} m³/s is the only class that leaves the continent)"
        );

        eprintln!("\n   the 8 largest spillways, fully classified:");
        let mut byq: Vec<&ymir_core::tectonics_c1::drainage::Spillway> =
            bs.spillways.iter().collect();
        byq.sort_by(|a, b| {
            b.discharge_m3s.partial_cmp(&a.discharge_m3s).unwrap_or(std::cmp::Ordering::Equal)
        });
        eprintln!(
            "   {:>10} {:>10} {:>30} {:>10} {:>12} {:>5} {:>14}",
            "Q m³/s", "source", "terminus class", "did", "below-sea id", "wc", "chained_into"
        );
        for sw in byq.iter().take(8) {
            let (cls, did, bid, cl) = classify(sw);
            eprintln!(
                "   {:>10.2} {:>10} {:>30} {:>10} {:>12} {:>5} {:>14}",
                sw.discharge_m3s,
                sw.lake_id,
                cls_name(cls),
                did,
                bid,
                cl,
                sw.chained_into.map_or("None".to_string(), |v| v.to_string())
            );
        }

        // ── B2 · the `wc = 0` termini, classified ─────────────────────────────────
        // Finding 71 called these "three spillways end on a LAND cell" and Finding 72 corrected
        // the claim in place: `water_class` has no value for a surface lake above sea level, so
        // `wc = 0` means "not ocean and not an enclosed below-sea region" — which INCLUDES every
        // surface lake. This block replaces the guess with the class.
        let mut b2: Vec<(&ymir_core::tectonics_c1::drainage::Spillway, Terminus, u32, u32)> = bs
            .spillways
            .iter()
            .filter_map(|sw| {
                let &(lx, ly) = sw.points.last().unwrap();
                let k = ly as usize * w + lx as usize;
                if wc[k] != 0 {
                    return None;
                }
                let (cls, did, bid, _) = classify(sw);
                Some((sw, cls, did, bid))
            })
            .collect();
        b2.sort_by(|a, b| {
            b.0.discharge_m3s.partial_cmp(&a.0.discharge_m3s).unwrap_or(std::cmp::Ordering::Equal)
        });
        let b2_q: f64 = b2.iter().map(|e| e.0.discharge_m3s as f64).sum();
        eprintln!(
            "\n── B2 · THE `wc = 0` TERMINI, CLASSIFIED ── {} spillway(s), Σ {b2_q:.2} m³/s \
             ({:.1} % of the hole)\n   {:>10} {:>10} {:>30} {:>10} {:>12} {:>14}",
            b2.len(),
            100.0 * b2_q / hole.max(1e-9),
            "Q m³/s",
            "source",
            "terminus class",
            "did",
            "below-sea id",
            "chained_into"
        );
        for (sw, cls, did, bid) in b2.iter() {
            eprintln!(
                "   {:>10.2} {:>10} {:>30} {:>10} {:>12} {:>14}",
                sw.discharge_m3s,
                sw.lake_id,
                cls_name(*cls),
                did,
                bid,
                sw.chained_into.map_or("None".to_string(), |v| v.to_string())
            );
        }
        let mut b2_cls: HashMap<Terminus, usize> = HashMap::new();
        for (_, cls, _, _) in b2.iter() {
            *b2_cls.entry(*cls).or_default() += 1;
        }
        eprintln!(
            "   ⇒ of {} `wc = 0` termini: {} surface lake · {} absorbed receiver · {} basin \
             with spillway · {} basin without · {} TRUE DRY LAND",
            b2.len(),
            b2_cls.get(&Terminus::SurfaceLake).copied().unwrap_or(0),
            b2_cls.get(&Terminus::AbsorbedReceiver).copied().unwrap_or(0),
            b2_cls.get(&Terminus::BasinWithSpillway).copied().unwrap_or(0),
            b2_cls.get(&Terminus::BasinWithoutSpillway).copied().unwrap_or(0),
            b2_cls.get(&Terminus::DryLand).copied().unwrap_or(0)
        );

        // ── 2b · WHAT THE `true DRY LAND` CLASS ACTUALLY IS ──────────────────────
        // Finding 71 said "ends on a LAND cell", Finding 72 corrected it to "a class
        // `water_class` cannot represent", and the class is STILL the largest one here. It is
        // not read as dry land without a test this time. The discriminator: the 8-connected
        // `wc == 2` COMPONENT of the terminal cell against the component(s) of the source
        // basin's own footprint. `break_reciprocal_spill_cycles` MERGES two regions under one
        // id (`lake_map[absorb] = keep`) and resets `chained_into = None`, so a basin whose
        // footprint spans MORE THAN ONE component is a merge survivor — and its spillway can
        // land in the other half while `bs.lake_map` reads its own id.
        let comp_of = {
            let mut lab = vec![0u32; n];
            let mut next = 0u32;
            let mut stack: Vec<u32> = Vec::new();
            for s0 in 0..n {
                if wc[s0] != 2 || lab[s0] != 0 {
                    continue;
                }
                next += 1;
                lab[s0] = next;
                stack.push(s0 as u32);
                while let Some(kk) = stack.pop() {
                    let k = kk as usize;
                    let (x, y) = ((k % w) as i32, (k / w) as i32);
                    for (dx, dy) in
                        ymir_core::terrain::flow::D8_DX.iter().zip(ymir_core::terrain::flow::D8_DY)
                    {
                        let (nx, ny) = (x + dx, y + dy);
                        if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                            continue;
                        }
                        let nk = ny as usize * w + nx as usize;
                        if wc[nk] == 2 && lab[nk] == 0 {
                            lab[nk] = next;
                            stack.push(nk as u32);
                        }
                    }
                }
            }
            eprintln!("\n── 2b · the `true DRY LAND` class, tested ── {next} below-sea regions");
            lab
        };
        let mut comps_of_lake: HashMap<u32, HashSet<u32>> = HashMap::new();
        for k in 0..n {
            if wc[k] == 2 && bs.lake_map[k] != 0 {
                comps_of_lake.entry(bs.lake_map[k]).or_default().insert(comp_of[k]);
            }
        }
        let merged: Vec<u32> =
            comps_of_lake.iter().filter(|(_, c)| c.len() > 1).map(|(&id, _)| id).collect();
        eprintln!(
            "   {} of {} below-sea ids span MORE THAN ONE region (⇒ merge survivors): {:?}",
            merged.len(),
            comps_of_lake.len(),
            {
                let mut m = merged.clone();
                m.sort();
                m
            }
        );
        eprintln!(
            "   {:>10} {:>10} {:>9} {:>9} {:>8} {:>28}",
            "Q m³/s", "source", "own cmp", "end cmp", "#cmps", "verdict"
        );
        let mut dry_loop = 0.0f64;
        let mut dry_merged = 0.0f64;
        let mut dry_other = 0.0f64;
        for sw in byq.iter() {
            let (cls, _, _, _) = classify(sw);
            if cls != Terminus::DryLand {
                continue;
            }
            let &(lx, ly) = sw.points.last().unwrap();
            let k = ly as usize * w + lx as usize;
            let end_c = comp_of[k];
            let own = comps_of_lake.get(&sw.lake_id).cloned().unwrap_or_default();
            // the component holding this basin's own FLOOR: the one with the most cells
            let mut counts: HashMap<u32, usize> = HashMap::new();
            for kk in 0..n {
                if bs.lake_map[kk] == sw.lake_id && wc[kk] == 2 {
                    *counts.entry(comp_of[kk]).or_default() += 1;
                }
            }
            let own_main = counts.iter().max_by_key(|e| *e.1).map(|e| *e.0).unwrap_or(0);
            let verdict = if end_c == 0 {
                dry_other += sw.discharge_m3s as f64;
                "TRUE dry land (wc≠2)"
            } else if end_c == own_main {
                dry_loop += sw.discharge_m3s as f64;
                "LOOP into its own region"
            } else if own.contains(&end_c) {
                dry_merged += sw.discharge_m3s as f64;
                "MERGED half (reset-1932)"
            } else {
                dry_other += sw.discharge_m3s as f64;
                "another region, id mismatch"
            };
            eprintln!(
                "   {:>10.2} {:>10} {:>9} {:>9} {:>8} {:>28}",
                sw.discharge_m3s,
                sw.lake_id,
                own_main,
                end_c,
                own.len(),
                verdict
            );
        }
        eprintln!(
            "   ⇒ Σ by verdict: LOOP {dry_loop:.1} · MERGED-half {dry_merged:.1} · other/true \
             dry {dry_other:.1} m³/s   [Spillway, net of evaporation, local + chain-inherited]"
        );

        // ── A2 · the chain graph, cycles, and the 240.8 chain ─────────────────────
        eprintln!("\n── A2 · CHAIN GRAPH ──");
        let edges: Vec<(u32, u32)> =
            bs.spillways.iter().filter_map(|s| s.chained_into.map(|t| (s.lake_id, t))).collect();
        let mut adj: HashMap<u32, Vec<u32>> = HashMap::new();
        for &(a, b) in &edges {
            adj.entry(a).or_default().push(b);
        }
        eprintln!("   {} chaining edges (source → chained_into):", edges.len());
        let mut sorted_edges = edges.clone();
        sorted_edges.sort();
        for (a, b) in &sorted_edges {
            eprintln!("     {a} → {b}");
        }
        // cycle detection by walking each source with a visited set
        let mut cycles = 0usize;
        for &(start, _) in &edges {
            let mut seen = HashSet::new();
            let mut cur = start;
            loop {
                if !seen.insert(cur) {
                    cycles += 1;
                    eprintln!("   ⛔ CYCLE reached from {start}, revisiting {cur}");
                    break;
                }
                match adj.get(&cur).and_then(|v| v.first()) {
                    Some(&nx) => cur = nx,
                    None => break,
                }
            }
        }
        eprintln!(
            "   cycles found: {cycles} {}",
            if cycles == 0 { "— none, the chain graph is a DAG as the code claims" } else { "" }
        );
        // the biggest spillway's chain, end to end
        if let Some(top) = byq.first() {
            eprintln!(
                "\n   the chain of the largest spillway ({:.2} m³/s), link by link:",
                top.discharge_m3s
            );
            let mut cur = Some(top.lake_id);
            let mut guard = 0;
            while let Some(id) = cur {
                guard += 1;
                if guard > 20 {
                    eprintln!("     … aborted at 20 links");
                    break;
                }
                let sw = bs.spillways.iter().find(|s| s.lake_id == id);
                let sum = bs.basins.iter().find(|b| b.id == id);
                match (sw, sum) {
                    (Some(sw), Some(b)) => {
                        let (cls, did, bid, cl) = classify(sw);
                        eprintln!(
                            "     basin {id}: IN {:.2} (local+chain) − EVAP {:.2} = OUT {:.2} \
                             m³/s | area {:.2} km² | exorheic {} | terminus {} (did {did}, bs \
                             {bid}, wc {cl})",
                            b.inflow_m3s,
                            b.evaporation_m3s,
                            sw.discharge_m3s,
                            b.area_km2,
                            b.exorheic,
                            cls_name(cls)
                        );
                        cur = sw.chained_into;
                    }
                    (None, Some(b)) => {
                        eprintln!(
                            "     basin {id}: IN {:.2} − EVAP {:.2} | area {:.2} km² | exorheic \
                             {} | **NO SPILLWAY — the chain ends here**",
                            b.inflow_m3s, b.evaporation_m3s, b.area_km2, b.exorheic
                        );
                        cur = None;
                    }
                    _ => {
                        eprintln!(
                            "     id {id}: not a below-sea basin — the chain leaves the machinery here"
                        );
                        cur = None;
                    }
                }
            }
        }

        // ── A1 · the ledger, per basin ────────────────────────────────────────────
        eprintln!("\n── A1 · THE LEDGER (source: BasinSummary, post-fixpoint) ──");
        let mut rows: Vec<(u32, f64, f64, f64, f64, f64, f32, bool)> = bs
            .basins
            .iter()
            .map(|b| {
                let inc = in_chain.get(&b.id).copied().unwrap_or(0.0);
                let out = out_of.get(&b.id).copied().unwrap_or(0.0);
                let inn = b.inflow_m3s as f64;
                (
                    b.id,
                    inn,
                    inc,
                    b.evaporation_m3s as f64,
                    out,
                    inn - b.evaporation_m3s as f64 - out,
                    b.area_km2,
                    b.exorheic,
                )
            })
            .collect();
        rows.sort_by(|a, b| b.5.abs().partial_cmp(&a.5.abs()).unwrap_or(std::cmp::Ordering::Equal));
        let (sum_in, sum_ev, sum_out, sum_res): (f64, f64, f64, f64) = rows
            .iter()
            .fold((0.0, 0.0, 0.0, 0.0), |(a, b, c, d), r| (a + r.1, b + r.3, c + r.4, d + r.5));
        eprintln!(
            "   {} basins | Σ IN(local+chain) {sum_in:.1} | Σ EVAP {sum_ev:.1} | Σ OUT \
             {sum_out:.1} | **Σ RESIDUAL {sum_res:.1}** m³/s",
            rows.len()
        );
        eprintln!(
            "   EVAP explains {:.1} % of the hole ({hole:.1} m³/s); residuals explain {:.1} %",
            100.0 * sum_ev / hole.max(1e-9),
            100.0 * sum_res / hole.max(1e-9)
        );
        // THE TWO PERIMETER TERMS AS NAMED LINES, so nobody chases residuals beneath them.
        let pe_open = runoff_km2_to_m3s(b_lakefoot as f32) as f64;
        let rain_below = runoff_km2_to_m3s(b_belowsea as f32) as f64;
        eprintln!(
            "   ── the perimeter, as LINES of the ledger (not residuals) ──\n   \
             {:<52} {:>10.2} m³/s  ({:.1} % of the hole)\n   {:<52} {:>10.2} m³/s  ({:.1} %)\n   \
             {:<52} {:>10.2} m³/s  ({:.1} %)\n   NO residual below {:.2} m³/s is interpretable — \
             it sits under the perimeter floor.",
            "open-water PE never deducted from the A0 budget",
            pe_open,
            100.0 * pe_open / hole.max(1e-9),
            "rain never seeded below sea level (counted NOWHERE)",
            rain_below,
            100.0 * rain_below / hole.max(1e-9),
            "the two together",
            pe_open + rain_below,
            100.0 * (pe_open + rain_below) / hole.max(1e-9),
            pe_open + rain_below
        );
        eprintln!(
            "   the 10 largest |residual| rows:\n   {:>10} {:>12} {:>12} {:>10} {:>10} {:>12} \
             {:>10} {:>9}",
            "id", "IN tot", "of which chain", "EVAP", "OUT", "RESIDUAL", "area km²", "exorheic"
        );
        for r in rows.iter().take(10) {
            eprintln!(
                "   {:>10} {:>12.2} {:>12.2} {:>10.2} {:>10.2} {:>12.2} {:>10.2} {:>9}",
                r.0, r.1, r.2, r.3, r.4, r.5, r.6, r.7
            );
        }

        // ── 4a · the missing handshake, per surface-lake terminus ─────────────────
        eprintln!(
            "\n── 4a · THE MISSING HANDSHAKE — surface lakes that receive a dropped surplus ──"
        );
        let mut recv: HashMap<u32, f64> = HashMap::new();
        for sw in &bs.spillways {
            let (cls, did, _, _) = classify(sw);
            if cls == Terminus::SurfaceLake {
                *recv.entry(did).or_default() += sw.discharge_m3s as f64;
            }
        }
        let mut rl: Vec<(&u32, &f64)> = recv.iter().collect();
        rl.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap_or(std::cmp::Ordering::Equal));
        eprintln!(
            "   {:>10} {:>12} {:>14} {:>14} {:>10}",
            "lake", "dropped Q", "own outlet Q", "should be", "shortfall"
        );
        let mut repaired: Vec<(u32, f64)> = Vec::new();
        for e in rl.iter().take(8) {
            let (did, dropped) = (*e.0, *e.1);
            let own = foot
                .get(&did)
                .map(|c| c.iter().map(|&k| runoff[k]).fold(0.0f32, f32::max))
                .unwrap_or(0.0);
            let own_m3s = runoff_km2_to_m3s(own) as f64;
            eprintln!(
                "   {:>10} {:>12.2} {:>14.3} {:>14.2} {:>10.2}",
                did,
                dropped,
                own_m3s,
                own_m3s + dropped,
                dropped
            );
            repaired.push((did, own_m3s + dropped));
        }

        // ── 4b · arithmetic only: the rivers conservation would produce ───────────
        eprintln!(
            "\n── 4b · ARITHMETIC ONLY — the five biggest 'repaired' outlets, at ratio 7.5 ──"
        );
        repaired.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        eprintln!(
            "   {:>10} {:>14} {:>16} {:>14}",
            "lake", "Q repaired", "width m @7.5", "cells @48.83 m"
        );
        for (did, qq) in repaired.iter().take(5) {
            let wsig = 5.0 * (qq * 56.25).sqrt();
            eprintln!("   {:>10} {:>14.2} {:>16.1} {:>14.2}", did, qq, wsig, wsig / 48.83);
        }
        eprintln!(
            "   (5·√(Q × 56.25) / 48.83. Input to the AUTHOR's conservation-vs-representation \
             decision, not a recommendation.)"
        );

        // ── C · the drainage_km2 denominator, in DIAG ─────────────────────────────
        let measured_mm = b_land / (n_land as f64 * cell_km2 as f64);
        eprintln!(
            "\n── C · `drainage_km2` DENOMINATOR BIAS ── measured runoff {measured_mm:.1} mm/yr \
             against the {REFERENCE_RUNOFF_MM:.0} proxy ⇒ effective area ×{:.3}",
            REFERENCE_RUNOFF_MM as f64 / measured_mm
        );
        let th = &dcfg.thresholds;
        let classify_km2 = |a: f32| -> &'static str {
            if a >= th.ship_km2 {
                "ship"
            } else if a >= th.barge_km2 {
                "barge"
            } else if a >= th.small_boat_km2 {
                "boat"
            } else {
                "non-nav"
            }
        };
        for ratio in [1.0f32, 7.5] {
            let scale = ratio * ratio;
            let mut before: HashMap<&str, usize> = HashMap::new();
            let mut after: HashMap<&str, usize> = HashMap::new();
            for i in 0..dr.rivers.segments.len() {
                let old = dr.segment_drainage_km2[i] * scale;
                let new = old * REFERENCE_RUNOFF_MM / measured_mm as f32;
                *before.entry(classify_km2(old)).or_default() += 1;
                *after.entry(classify_km2(new)).or_default() += 1;
            }
            let g = |m: &HashMap<&str, usize>, k: &str| m.get(k).copied().unwrap_or(0);
            let nav_b = g(&before, "boat") + g(&before, "barge") + g(&before, "ship");
            let nav_a = g(&after, "boat") + g(&after, "barge") + g(&after, "ship");
            eprintln!(
                "   ratio {ratio}: boat {}→{} | barge {}→{} | ship {}→{} | NAVIGABLE {nav_b}→\
                 {nav_a} ({:+.1} %)   [Watercourse population, pre-clip drainage]",
                g(&before, "boat"),
                g(&after, "boat"),
                g(&before, "barge"),
                g(&after, "barge"),
                g(&before, "ship"),
                g(&after, "ship"),
                100.0 * (nav_a as f64 - nav_b as f64) / nav_b.max(1) as f64
            );
        }
    }
}
