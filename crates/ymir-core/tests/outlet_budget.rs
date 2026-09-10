//! ADR Finding 71 — the three worlds for one number: is the 240.8 m³/s spillway a legitimate
//! terminal collector, the terminus of a continuous trunk, or double counting?
//!
//! **Every discharge reported here carries its `SegmentKind`.** No exceptions.
//!
//! Closure target, declared before measuring: the **MEASURED generated runoff**
//! (`Σ max(0, precip − PE) · cell_km2` over land), not the 300 mm reference. The reference is a
//! proxy; the generated runoff is the authority. **Admissible closure error: a factor 2.**
//!
//! Run: cargo test -p ymir-core --release --test outlet_budget -- --ignored --nocapture

use std::collections::HashMap;

use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::{PrecipParams, precip_mm_per_year};
use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, LakeType, SegmentKind, c1_drainage_windowed,
    potential_evaporation_mm, runoff_km2_to_m3s,
};
use ymir_core::tectonics_c1::hd_assembly::assemble_hd_drainage;
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::upscale_from_c1_with_progress;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::flow::{DIR_NONE, FlowConfig, breach_monotone, compute_flow};
use ymir_core::terrain::upscale::{ProductionHdOpts, production_hd_config};

const PSEED: u64 = 10_481_999_410_520_546_993;
const DOMAIN_KM: f32 = 400.0;
const SEA: f32 = 0.5;
const TARGET: usize = 8192;
const CELL_KM: f32 = DOMAIN_KM / TARGET as f32;
const BEDS: [(&str, f32, f32); 2] = [("humid", 45.0, 40.0), ("arid-hot", 25.0, 10.0)];
/// The reference the code uses to turn discharge into an "effective area". A proxy, kept only
/// for comparison against the measured value.
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

fn q(v: &mut Vec<f32>, f: f64) -> f32 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    if v.is_empty() {
        return f32::NAN;
    }
    v[(((v.len() - 1) as f64 * f).round() as usize).min(v.len() - 1)]
}

fn kind_str(k: SegmentKind) -> &'static str {
    match k {
        SegmentKind::Watercourse => "Watercourse",
        SegmentKind::Spillway => "Spillway",
    }
}

#[test]
#[ignore]
fn outlet_budget_and_the_three_worlds() {
    let ss = SteinSteinParams::default();
    eprintln!(
        "\n==========  A · the outlet budget, and the arbitration of the three worlds  =========="
    );
    eprintln!(
        "  closure target = MEASURED generated runoff Σ max(0, precip − PE)·cell_km2 over land.\n  \
         ADMISSIBLE CLOSURE ERROR: factor 2, declared before measuring.\n  \
         every discharge below carries its SegmentKind."
    );
    let raw = terrain();
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let pre0 = c1_drainage_windowed(&raw, None, &dcfg, &ss, DOMAIN_KM);
    let field =
        breach_monotone(&raw, &pre0.flow.filled, &pre0.lake_map, SEA, raw.width, raw.height);
    let (w, h) = (field.width, field.height);
    let cell_km2 = CELL_KM * CELL_KM;
    let flow = compute_flow(&field, &FlowConfig { sea_level: SEA, ..Default::default() });
    // A coastal OUTLET cell: land whose D8 receiver is sea, or which routes off the map.
    let is_coastal_outlet = |k: usize| -> bool {
        if field.data[k] <= SEA {
            return false;
        }
        let d = flow.direction[k];
        if d == DIR_NONE {
            return true;
        }
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        let (nx, ny) = (
            x + ymir_core::terrain::flow::D8_DX[d as usize],
            y + ymir_core::terrain::flow::D8_DY[d as usize],
        );
        if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
            return true;
        }
        field.data[ny as usize * w + nx as usize] <= SEA
    };

    for (bed, lat, span) in BEDS {
        eprintln!("\n╔══════ {bed} ({lat}°/span {span}) ══════╗");
        let climate =
            c1_climate_placed(&field, &ss, lat, span, &PrecipParams::default(), DOMAIN_KM);

        // ── A0 · the budget ────────────────────────────────────────────────────────
        let (mut gen_mm_km2, mut land) = (0.0f64, 0usize);
        for k in 0..w * h {
            if field.data[k] > SEA {
                land += 1;
                let p = precip_mm_per_year(climate.precipitation.data[k]);
                let pe = potential_evaporation_mm(climate.temperature.data[k]);
                gen_mm_km2 += ((p - pe).max(0.0) * cell_km2) as f64;
            }
        }
        let land_km2 = land as f64 * cell_km2 as f64;
        let budget_m3s = runoff_km2_to_m3s(gen_mm_km2 as f32) as f64;
        let depth_mm = gen_mm_km2 / land_km2;
        eprintln!(
            "\n── A0 · BUDGET ── land {:.0} km² | generated runoff {:.1} mm/yr (reference \
             {REFERENCE_RUNOFF_MM:.0}, ×{:.2}) ⇒ **{budget_m3s:.1} m³/s**",
            land_km2,
            depth_mm,
            depth_mm / REFERENCE_RUNOFF_MM as f64
        );

        let dclim = DrainageClimate {
            precip_internal: &climate.precipitation,
            temperature: &climate.temperature,
        };
        let pre = c1_drainage_windowed(&raw, None, &dcfg, &ss, DOMAIN_KM);
        let dr = assemble_hd_drainage(
            &field,
            &dclim,
            Some(pre),
            &dcfg,
            &ss,
            DOMAIN_KM,
            1.0,
            None,
            false,
        )
        .drainage;
        let nseg = dr.rivers.segments.len();
        let endo_ids: Vec<u32> = dr
            .lakes
            .iter()
            .filter(|l| l.lake_type == LakeType::Endorheic)
            .map(|l| l.base.id)
            .collect();

        // ── A1 · the outlets ──────────────────────────────────────────────────────
        // TERMINAL, declared: a segment whose downstream-most cell is a coastal outlet, OR
        // whose downstream-most cell touches an ENDORHEIC lake (the water dies there). A lake
        // outlet upstream of a coastal outlet is NOT terminal and is excluded from the sum, so
        // nothing is counted twice.
        let mut terminal: Vec<(usize, f32, SegmentKind, &'static str)> = Vec::new();
        for i in 0..nseg {
            let &(ex, ey) = dr.rivers.segments[i].points.last().unwrap_or(&(0, 0));
            let k = ey as usize * w + ex as usize;
            let q_i = dr.segment_discharge_m3s[i];
            if is_coastal_outlet(k) {
                terminal.push((i, q_i, dr.segment_kind[i], "sea"));
                continue;
            }
            let mut endo = false;
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let (nx, ny) = (ex as i32 + dx, ey as i32 + dy);
                    if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                        continue;
                    }
                    let id = dr.lake_map[ny as usize * w + nx as usize];
                    if id != 0 && endo_ids.contains(&id) {
                        endo = true;
                    }
                }
            }
            if endo {
                terminal.push((i, q_i, dr.segment_kind[i], "endorheic"));
            }
        }
        terminal.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        let sum_terminal: f64 = terminal.iter().map(|t| t.1 as f64).sum();
        let biggest = terminal.first().map(|t| t.1).unwrap_or(0.0);
        eprintln!(
            "\n── A1 · TERMINAL OUTLETS ── {} of {nseg} segments | Σ {:.1} m³/s | budget \
             {budget_m3s:.1} ⇒ **closure ×{:.3}** {}",
            terminal.len(),
            sum_terminal,
            sum_terminal / budget_m3s.max(1e-9),
            if (sum_terminal / budget_m3s.max(1e-9)) <= 2.0
                && (sum_terminal / budget_m3s.max(1e-9)) >= 0.5
            {
                "— WITHIN the declared factor 2"
            } else {
                "— OUTSIDE the declared factor 2"
            }
        );
        eprintln!(
            "   biggest outlet {biggest:.1} m³/s = {:.1} % of the budget | exceeds budget? {}",
            100.0 * biggest as f64 / budget_m3s.max(1e-9),
            if biggest as f64 > budget_m3s { "⛔ YES — DOUBLE COUNTING" } else { "no" }
        );
        eprintln!("   top 12 terminal outlets:");
        eprintln!(
            "   {:>10} {:>13} {:>12} {:>12} {:>10}",
            "Q m³/s", "KIND", "catch km²", "source lake", "ends in"
        );
        for &(i, qq, kd, endsin) in terminal.iter().take(12) {
            eprintln!(
                "   {qq:>10.2} {:>13} {:>12.1} {:>12} {endsin:>10}",
                kind_str(kd),
                dr.segment_catchment_cells[i] * cell_km2,
                dr.segment_source_lake[i].map_or("-".to_string(), |v| v.to_string())
            );
        }
        // multi-outlet lakes: how many distinct spillways name the same source lake?
        let mut per_lake: HashMap<u32, Vec<usize>> = HashMap::new();
        for i in 0..nseg {
            if dr.segment_kind[i] == SegmentKind::Spillway {
                if let Some(id) = dr.segment_source_lake[i] {
                    per_lake.entry(id).or_default().push(i);
                }
            }
        }
        let multi: Vec<(&u32, &Vec<usize>)> =
            per_lake.iter().filter(|(_, v)| v.len() > 1).collect();
        eprintln!(
            "   lakes with MULTIPLE spillways: {} of {} inventoried spillway sources{}",
            multi.len(),
            per_lake.len(),
            if multi.is_empty() {
                " — the double-counting mechanism is absent".to_string()
            } else {
                format!(
                    " — Σ of their outlets {:.1} m³/s",
                    multi
                        .iter()
                        .flat_map(|(_, v)| v.iter())
                        .map(|&i| dr.segment_discharge_m3s[i] as f64)
                        .sum::<f64>()
                )
            }
        );

        // ── A2 · the dead column, then the discriminating instrument ──────────────
        eprintln!("\n── A2 · the DEAD COLUMN: spillway share by Strahler order ──");
        for o in 1u8..=6 {
            let n: usize = (0..nseg).filter(|&i| dr.rivers.segments[i].strahler_order == o).count();
            let sp: usize = (0..nseg)
                .filter(|&i| {
                    dr.rivers.segments[i].strahler_order == o
                        && dr.segment_kind[i] == SegmentKind::Spillway
                })
                .count();
            if n > 0 {
                eprintln!("   S{o}: {sp} spillways of {n} = {:.2} %", 100.0 * sp as f32 / n as f32);
            }
        }
        eprintln!(
            "   ⇒ as predicted, 100 % of spillways sit at order 1 BY CONSTRUCTION\n     \
             (`strahler_order: 1 // MEANINGLESS on a spillway`). The column measures the\n     \
             constant, not the terrain. Recorded as dead; replaced below."
        );
        let mut order: Vec<usize> = (0..nseg).collect();
        order.sort_by(|&a, &b| {
            dr.segment_discharge_m3s[b]
                .partial_cmp(&dr.segment_discharge_m3s[a])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        eprintln!("\n── A2 · spillway share by DISCHARGE DECILE (count / length) ──");
        for d in 0..10 {
            let lo = d * nseg / 10;
            let hi = ((d + 1) * nseg / 10).min(nseg);
            let slice = &order[lo..hi];
            let sp = slice.iter().filter(|&&i| dr.segment_kind[i] == SegmentKind::Spillway).count();
            let len_all: usize = slice.iter().map(|&i| dr.rivers.segments[i].points.len()).sum();
            let len_sp: usize = slice
                .iter()
                .filter(|&&i| dr.segment_kind[i] == SegmentKind::Spillway)
                .map(|&i| dr.rivers.segments[i].points.len())
                .sum();
            eprintln!(
                "   decile {} (Q {:.4}–{:.2}): {sp:>4} spillways of {:>5} = {:>6.2} % | \
                 {:>6.2} % of length",
                d + 1,
                dr.segment_discharge_m3s[slice[slice.len() - 1]],
                dr.segment_discharge_m3s[slice[0]],
                slice.len(),
                100.0 * sp as f32 / slice.len().max(1) as f32,
                100.0 * len_sp as f32 / len_all.max(1) as f32
            );
        }
        eprintln!("\n── A2 · the TOP 20 SEGMENTS BY DISCHARGE, named ──");
        eprintln!(
            "   {:>4} {:>10} {:>13} {:>12} {:>12} {:>8} {:>9}",
            "#", "Q m³/s", "KIND", "catch km²", "source lake", "order", "cells"
        );
        for (r, &i) in order.iter().take(20).enumerate() {
            eprintln!(
                "   {:>4} {:>10.2} {:>13} {:>12.1} {:>12} {:>8} {:>9}",
                r + 1,
                dr.segment_discharge_m3s[i],
                kind_str(dr.segment_kind[i]),
                dr.segment_catchment_cells[i] * cell_km2,
                dr.segment_source_lake[i].map_or("-".to_string(), |v| v.to_string()),
                dr.rivers.segments[i].strahler_order,
                dr.rivers.segments[i].points.len()
            );
        }
        let sp_top20 =
            order.iter().take(20).filter(|&&i| dr.segment_kind[i] == SegmentKind::Spillway).count();
        eprintln!("   ⇒ {sp_top20} of the top 20 are Spillway");

        // ── A1/A3 quantiles by kind ───────────────────────────────────────────────
        eprintln!("\n── A1b · discharge distribution BY KIND ──");
        for kd in [SegmentKind::Watercourse, SegmentKind::Spillway] {
            let mut v: Vec<f32> = (0..nseg)
                .filter(|&i| dr.segment_kind[i] == kd)
                .map(|i| dr.segment_discharge_m3s[i])
                .collect();
            let len: usize = (0..nseg)
                .filter(|&i| dr.segment_kind[i] == kd)
                .map(|i| dr.rivers.segments[i].points.len())
                .sum();
            if v.is_empty() {
                eprintln!("   {:<13} none", kind_str(kd));
                continue;
            }
            eprintln!(
                "   {:<13} n {:>6} | length {:>8} cells = {:>7.0} km | p50 {:>8.3} | p90 \
                 {:>8.3} | p99 {:>8.3} | max {:>9.3} m³/s",
                kind_str(kd),
                v.len(),
                len,
                len as f32 * CELL_KM,
                q(&mut v.clone(), 0.50),
                q(&mut v.clone(), 0.90),
                q(&mut v.clone(), 0.99),
                q(&mut v.clone(), 1.00)
            );
        }

        // ── A3 · the trunk trace ──────────────────────────────────────────────────
        // From the biggest terminal outlet's cell, walk UPSTREAM choosing at each step the
        // donor with the greatest D8 accumulation. `kind`/`Q` at a cell are read from the
        // segment covering it.
        let mut cover: HashMap<usize, usize> = HashMap::new();
        for i in 0..nseg {
            for &(x, y) in &dr.rivers.segments[i].points {
                cover.insert(y as usize * w + x as usize, i);
            }
        }
        let mut qat: HashMap<usize, f32> = HashMap::new();
        for i in 0..nseg {
            for (j, &(x, y)) in dr.rivers.segments[i].points.iter().enumerate() {
                let k = y as usize * w + x as usize;
                let v = dr.segment_discharge_profile_m3s[i].get(j).copied().unwrap_or(0.0);
                qat.entry(k).and_modify(|e| *e = e.max(v)).or_insert(v);
            }
        }
        if let Some(&(bi, bq, bk, _)) = terminal.first() {
            let &(sx, sy) = dr.rivers.segments[bi].points.last().unwrap_or(&(0, 0));
            let mut k = sy as usize * w + sx as usize;
            let mut path: Vec<usize> = vec![k];
            for _ in 0..2_000_000 {
                // donors of k: neighbours whose D8 receiver is k
                let (x, y) = ((k % w) as i32, (k / w) as i32);
                let mut best: Option<(usize, f32)> = None;
                for d in 0..8usize {
                    let (nx, ny) = (
                        x - ymir_core::terrain::flow::D8_DX[d],
                        y - ymir_core::terrain::flow::D8_DY[d],
                    );
                    if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                        continue;
                    }
                    let nk = ny as usize * w + nx as usize;
                    if field.data[nk] <= SEA || flow.direction[nk] == DIR_NONE {
                        continue;
                    }
                    let rd = flow.direction[nk] as usize;
                    let (rx, ry) = (
                        nx + ymir_core::terrain::flow::D8_DX[rd],
                        ny + ymir_core::terrain::flow::D8_DY[rd],
                    );
                    if rx < 0 || ry < 0 || rx as usize >= w || ry as usize >= h {
                        continue;
                    }
                    if ry as usize * w + rx as usize != k {
                        continue;
                    }
                    let a = flow.accumulation.data[nk];
                    if best.map_or(true, |(_, ba)| a > ba) {
                        best = Some((nk, a));
                    }
                }
                match best {
                    Some((nk, _)) => {
                        path.push(nk);
                        k = nk;
                    }
                    None => break,
                }
            }
            // kind alternations and the longest continuous watercourse run
            let kinds: Vec<Option<SegmentKind>> =
                path.iter().map(|k| cover.get(k).map(|&i| dr.segment_kind[i])).collect();
            let mut alt = 0usize;
            let (mut run, mut best_run) = (0usize, 0usize);
            let mut prev: Option<Option<SegmentKind>> = None;
            for kd in &kinds {
                if let Some(p) = prev {
                    if p != *kd {
                        alt += 1;
                    }
                }
                prev = Some(*kd);
                if *kd == Some(SegmentKind::Watercourse) {
                    run += 1;
                    best_run = best_run.max(run);
                } else {
                    run = 0;
                }
            }
            let covered = kinds.iter().filter(|k| k.is_some()).count();
            eprintln!(
                "\n── A3 · TRUNK TRACE from the biggest terminal outlet ({bq:.2} m³/s, {}) ──\n   \
                 length {} cells = {:.1} km | on a mapped segment: {covered} cells ({:.1} %) | \
                 kind alternations {alt} | longest continuous Watercourse run {best_run} cells \
                 ({:.1} % of the trace)",
                kind_str(bk),
                path.len(),
                path.len() as f32 * CELL_KM,
                100.0 * covered as f32 / path.len() as f32,
                100.0 * best_run as f32 / path.len() as f32
            );
            // discharge profile along the trace, sampled, with the biggest step
            let qs: Vec<f32> = path.iter().map(|k| qat.get(k).copied().unwrap_or(0.0)).collect();
            let (mut jump, mut jump_at) = (1.0f32, 0usize);
            for i in 1..qs.len() {
                let (a, b) = (qs[i], qs[i - 1]); // upstream a, downstream b
                if a > 1e-6 && b / a > jump {
                    jump = b / a;
                    jump_at = i;
                }
            }
            eprintln!(
                "   discharge along the trace (downstream→upstream), every 10 %: {}",
                (0..=10)
                    .map(|t| {
                        let i = (t * (qs.len() - 1) / 10).min(qs.len() - 1);
                        format!("{:.2}", qs[i])
                    })
                    .collect::<Vec<_>>()
                    .join(" → ")
            );
            eprintln!(
                "   BIGGEST STEP downstream/upstream = ×{jump:.1} at {:.1} % along the trace \
                 (Q {:.2} → {:.2}, kinds {:?} → {:?})",
                100.0 * jump_at as f32 / qs.len() as f32,
                qs[jump_at],
                qs[jump_at - 1],
                kinds[jump_at].map(kind_str),
                kinds[jump_at - 1].map(kind_str)
            );
        }

        // ── A4 · reconciling the counts ───────────────────────────────────────────
        let exo = dr.lakes.iter().filter(|l| l.lake_type == LakeType::Exorheic).count();
        let endo = endo_ids.len();
        let below_sea = dr.lakes.iter().filter(|l| l.base.id >= 1_000_001).count();
        let detected = dr.lakes.len() - below_sea;
        let nspill = (0..nseg).filter(|&i| dr.segment_kind[i] == SegmentKind::Spillway).count();
        let spill_named = (0..nseg).filter(|&i| dr.segment_source_lake[i].is_some()).count();
        eprintln!(
            "\n── A4 · COUNTS ──\n   lakes {} = {detected} detected + {below_sea} below-sea | \
             {exo} exorheic, {endo} endorheic\n   spillways {nspill} ({spill_named} naming an \
             inventoried source, {} naming none)\n   ⇒ spillways are produced ONLY by \
             `below_sea_basin_lakes_infil`, so the count to reconcile against is the BELOW-SEA \
             basin count ({below_sea}), not the exorheic lake count ({exo}). An exorheic lake \
             ABOVE sea level overflows through an ordinary Watercourse.",
            dr.lakes.len(),
            nspill - spill_named
        );
    }
}

/// A1 CORRECTED — my first terminal test excluded spillways, which is where the water is.
///
/// A spillway's traced path runs over a col, OUTSIDE the accumulation network, so its last cell
/// has no D8 sea receiver and my `is_coastal_outlet` test rejected it. The code itself gives the
/// right discriminator: `drainage.rs:2024` asserts `wc[end] == 1 || chained_into.is_some()`, so a
/// spillway terminates at the OCEAN iff `water_class` at its last cell is 1. A CHAINED spillway
/// (class 2, an enclosed below-sea region) must NOT be summed — its outflow is re-counted by the
/// spillway of the basin it feeds.
///
/// ⚠️ `chained_into` is computed and then DROPPED: it never reaches `C1DrainageResult` nor
/// `rivers.json`. So neither a consumer nor this bench can build the terminal set from the export
/// alone — it has to be reconstructed from `water_class`. That is an export-contract gap, and it
/// is the reason the water balance could not be closed until now.
#[test]
#[ignore]
fn outlet_budget_corrected() {
    use ymir_core::lakes::connectivity::water_class;
    let ss = SteinSteinParams::default();
    eprintln!(
        "\n==========  A1 CORRECTED · terminal = watercourse-to-sea + spillway-to-OCEAN  =========="
    );
    let raw = terrain();
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let pre0 = c1_drainage_windowed(&raw, None, &dcfg, &ss, DOMAIN_KM);
    let field =
        breach_monotone(&raw, &pre0.flow.filled, &pre0.lake_map, SEA, raw.width, raw.height);
    let (w, h) = (field.width, field.height);
    let cell_km2 = CELL_KM * CELL_KM;
    let flow = compute_flow(&field, &FlowConfig { sea_level: SEA, ..Default::default() });
    let wc = water_class(&field, SEA);

    for (bed, lat, span) in BEDS {
        let climate =
            c1_climate_placed(&field, &ss, lat, span, &PrecipParams::default(), DOMAIN_KM);
        let (mut gen_mm_km2, mut land) = (0.0f64, 0usize);
        for k in 0..w * h {
            if field.data[k] > SEA {
                land += 1;
                let p = precip_mm_per_year(climate.precipitation.data[k]);
                let pe = potential_evaporation_mm(climate.temperature.data[k]);
                gen_mm_km2 += ((p - pe).max(0.0) * cell_km2) as f64;
            }
        }
        let budget = runoff_km2_to_m3s(gen_mm_km2 as f32) as f64;
        let dclim = DrainageClimate {
            precip_internal: &climate.precipitation,
            temperature: &climate.temperature,
        };
        let pre = c1_drainage_windowed(&raw, None, &dcfg, &ss, DOMAIN_KM);
        let dr = assemble_hd_drainage(
            &field,
            &dclim,
            Some(pre),
            &dcfg,
            &ss,
            DOMAIN_KM,
            1.0,
            None,
            false,
        )
        .drainage;
        let nseg = dr.rivers.segments.len();

        let (mut sum_wc_sea, mut n_wc) = (0.0f64, 0usize);
        let (mut sum_sp_ocean, mut n_sp_o) = (0.0f64, 0usize);
        let (mut sum_sp_chained, mut n_sp_c) = (0.0f64, 0usize);
        let mut top: Vec<(f32, SegmentKind, u8, Option<u32>)> = Vec::new();
        for i in 0..nseg {
            let &(ex, ey) = dr.rivers.segments[i].points.last().unwrap_or(&(0, 0));
            let k = ey as usize * w + ex as usize;
            let qq = dr.segment_discharge_m3s[i];
            match dr.segment_kind[i] {
                SegmentKind::Spillway => {
                    top.push((qq, SegmentKind::Spillway, wc[k], dr.segment_source_lake[i]));
                    if wc[k] == 1 {
                        sum_sp_ocean += qq as f64;
                        n_sp_o += 1;
                    } else {
                        sum_sp_chained += qq as f64;
                        n_sp_c += 1;
                    }
                }
                SegmentKind::Watercourse => {
                    // to sea: the D8 receiver is sea, or it routes off the map
                    let d = flow.direction[k];
                    let to_sea = if d == DIR_NONE {
                        true
                    } else {
                        let (x, y) = ((k % w) as i32, (k / w) as i32);
                        let (nx, ny) = (
                            x + ymir_core::terrain::flow::D8_DX[d as usize],
                            y + ymir_core::terrain::flow::D8_DY[d as usize],
                        );
                        nx < 0
                            || ny < 0
                            || nx as usize >= w
                            || ny as usize >= h
                            || field.data[ny as usize * w + nx as usize] <= SEA
                    };
                    if to_sea {
                        sum_wc_sea += qq as f64;
                        n_wc += 1;
                    }
                }
            }
        }
        let terminal = sum_wc_sea + sum_sp_ocean;
        top.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
        eprintln!("\n╔══════ {bed} ══════╗  budget {budget:.1} m³/s");
        eprintln!(
            "   Watercourse → sea      : {n_wc:>4} segments, Σ {sum_wc_sea:>8.1} m³/s ({:>5.1} % \
             of budget)",
            100.0 * sum_wc_sea / budget
        );
        eprintln!(
            "   Spillway → OCEAN (wc=1): {n_sp_o:>4} segments, Σ {sum_sp_ocean:>8.1} m³/s \
             ({:>5.1} % of budget)   ← TERMINAL",
            100.0 * sum_sp_ocean / budget
        );
        eprintln!(
            "   Spillway → CHAINED     : {n_sp_c:>4} segments, Σ {sum_sp_chained:>8.1} m³/s \
             ({:>5.1} % of budget)   ← must NOT be summed",
            100.0 * sum_sp_chained / budget
        );
        eprintln!(
            "   ⇒ TERMINAL TOTAL {terminal:.1} m³/s / budget {budget:.1} = **×{:.3}** {}",
            terminal / budget,
            if (terminal / budget) <= 2.0 && (terminal / budget) >= 0.5 {
                "— WITHIN the declared factor 2"
            } else {
                "— OUTSIDE the declared factor 2"
            }
        );
        eprintln!(
            "   naive sum of ALL spillways {:.1} m³/s = ×{:.3} of budget — the chained ones are \
             the double count",
            sum_sp_ocean + sum_sp_chained,
            (sum_sp_ocean + sum_sp_chained + sum_wc_sea) / budget
        );
        eprintln!("   top 8 spillways, with the water_class of their last cell:");
        for (qq, kd, cls, src) in top.iter().take(8) {
            eprintln!(
                "     {qq:>9.2} m³/s  {}  wc={} ({})  source lake {}",
                kind_str(*kd),
                cls,
                match cls {
                    1 => "OCEAN — terminal",
                    2 => "enclosed below-sea — CHAINED",
                    _ => "land?!",
                },
                src.map_or("-".to_string(), |v| v.to_string())
            );
        }
    }
}
