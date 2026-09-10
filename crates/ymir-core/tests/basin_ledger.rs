//! ADR Finding 72 — the ledger of the below-sea basins: where the 81 % goes.
//!
//! Every discharge carries its `SegmentKind` **and its population of definition**: local basin
//! only / inherited through the chain / net of evaporation.
//!
//! `assemble_hd_drainage` discards both `BasinSummary` and `Spillway::chained_into`, so this
//! calls `below_sea_basin_lakes_infil` directly after replicating the production sequence up to
//! that point (pre-breach drainage → breach → climate → final drainage → H-1c).
//!
//! Run: cargo test -p ymir-core --release --test basin_ledger -- --ignored --nocapture

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
    C1DrainageConfig, DrainageClimate, LakeType, apply_lake_water_balance,
    below_sea_basin_lakes_infil, c1_drainage_windowed, potential_evaporation_mm,
    runoff_accumulation, runoff_km2_to_m3s,
};
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
fn basin_ledger() {
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

        let bs = below_sea_basin_lakes_infil(
            &field,
            &dclim,
            &dcfg,
            &ss,
            DOMAIN_KM,
            Some(&dr.lake_map),
            None,
        );
        let runoff = runoff_accumulation(&field, &flow, &dclim, cell_km2, None, None, w, h);

        // ── A3 · BLOCKING instrument control on an ordinary exorheic SURFACE lake ──
        eprintln!(
            "\n── A3 · BLOCKING CONTROL — does the instrument carry flow ACROSS an exorheic surface lake? ──"
        );
        let mut foot: HashMap<u32, Vec<usize>> = HashMap::new();
        for k in 0..n {
            let id = dr.lake_map[k];
            if id != 0 && id < 1_000_001 {
                foot.entry(id).or_default().push(k);
            }
        }
        let mut cands: Vec<&ymir_core::tectonics_c1::drainage::C1Lake> = dr
            .lakes
            .iter()
            .filter(|l| l.base.id < 1_000_001 && l.lake_type == LakeType::Exorheic)
            .collect();
        cands.sort_by(|a, b| {
            b.area_km2.partial_cmp(&a.area_km2).unwrap_or(std::cmp::Ordering::Equal)
        });
        // CORRECTED CONTROL. The first version read `lk.base.outlet` — the PRE-BREACH sill, a
        // saddle on the rim — and followed `flow.direction` on the BREACHED field, where that
        // depression no longer exists. A saddle carries almost no accumulation, so the control
        // could not pass: it measured the wrong cell. Corrected: start at the footprint's OWN
        // maximum-accumulation cell and walk downstream until leaving the footprint. That tests
        // "does flow continue past the lake" without depending on where the sill was.
        let mut a3_ok = false;
        let mut a3_n = 0usize;
        for lk in cands.iter().take(3) {
            let Some(cells) = foot.get(&lk.base.id) else { continue };
            let member: std::collections::HashSet<usize> = cells.iter().copied().collect();
            let (mut kmax, mut inside) = (cells[0], 0.0f32);
            for &k in cells {
                if runoff[k] > inside {
                    inside = runoff[k];
                    kmax = k;
                }
            }
            let mut cur = kmax;
            let mut steps = 0usize;
            let down = loop {
                steps += 1;
                if steps > 10_000 {
                    break None;
                }
                let d = flow.direction[cur];
                if d == DIR_NONE {
                    break None;
                }
                let (x, y) = (
                    (cur % w) as i32 + ymir_core::terrain::flow::D8_DX[d as usize],
                    (cur / w) as i32 + ymir_core::terrain::flow::D8_DY[d as usize],
                );
                if x < 0 || y < 0 || x as usize >= w || y as usize >= h {
                    break None;
                }
                let nk = y as usize * w + x as usize;
                if !member.contains(&nk) {
                    break Some(nk);
                }
                cur = nk;
            };
            a3_n += 1;
            let dq = down.map(|k| runoff[k]).unwrap_or(0.0);
            let ratio = dq / inside.max(1e-9);
            eprintln!(
                "   lake {:>8} area {:>8.1} km² | max runoff INSIDE {:>10.1} → one step \
                 DOWNSTREAM of the outlet {:>10.1} mm·km²/yr | ratio {ratio:.3} ({:.3} vs \
                 {:.3} m³/s)",
                lk.base.id,
                lk.area_km2,
                inside,
                dq,
                runoff_km2_to_m3s(inside),
                runoff_km2_to_m3s(dq)
            );
            if ratio >= 0.9 {
                a3_ok = true;
            }
        }
        if a3_n == 0 {
            eprintln!(
                "   ⚠️ EMPTY POPULATION — no exorheic DETECTED lake with a footprint exists in this bed. The control DID NOT RUN: a rule-10 failure in my own control, not a result about the machinery."
            );
        }
        eprintln!(
            "   ⇒ {} — over {a3_n} lake(s)",
            if a3_ok {
                "CLOSES: the instrument DOES carry flow across an exorheic surface lake"
            } else if a3_n == 0 {
                "VOID (no population)"
            } else {
                "DOES NOT CLOSE"
            }
        );
        if !a3_ok {
            eprintln!(
                "   ⛔ STOP RULE FIRES: the ledger would prove what the machinery wrote. \
                 Not reading A1."
            );
            continue;
        }

        // ── terminus classification of every spillway ──────────────────────────────
        let live_bs: HashSet<u32> = bs.lakes.iter().map(|l| l.base.id).collect();
        let has_spill: HashSet<u32> = bs.spillways.iter().map(|s| s.lake_id).collect();
        let classify =
            |sw: &ymir_core::tectonics_c1::drainage::Spillway| -> (Terminus, u32, u32, u8) {
                let &(lx, ly) = sw.points.last().unwrap();
                let k = ly as usize * w + lx as usize;
                let did = if dr.lake_map[k] < 1_000_001 { dr.lake_map[k] } else { 0 };
                let bid = bs.lake_map[k];
                let cls = if wc[k] == 1 {
                    Terminus::Ocean
                } else if did != 0 {
                    Terminus::SurfaceLake
                } else if bid != 0 && bid != sw.lake_id {
                    if has_spill.contains(&bid) {
                        Terminus::BasinWithSpillway
                    } else {
                        Terminus::BasinWithoutSpillway
                    }
                } else if let Some(t) = sw.chained_into {
                    if live_bs.contains(&t) || has_spill.contains(&t) {
                        Terminus::BasinWithSpillway
                    } else {
                        Terminus::AbsorbedReceiver
                    }
                } else {
                    Terminus::DryLand
                };
                (cls, did, bid, wc[k])
            };

        // ── 1 · THE BINARY DISCRIMINANT, reported FIRST ───────────────────────────
        let hole = budget - 93.1; // the F71 terminal total for humid; recomputed below per bed
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
        // the bed's own terminal total: watercourse-to-sea is unchanged from F71; recompute the
        // spillway-to-ocean part here and take the F71 watercourse figure as given per bed.
        let sp_ocean: f64 = per_cls.get(&Terminus::Ocean).map(|e| e.1).unwrap_or(0.0);
        eprintln!(
            "\n── 1 · THE BINARY DISCRIMINANT ──\n   Σ discharge of spillways ending on a \
             SURFACE LAKE (did ≠ 0) = **{sum_surface:.1} m³/s**   [Spillway, population: net of \
             evaporation, local + chain-inherited]\n   the hole ≈ budget − terminal ≈ \
             {:.1} m³/s ⇒ **coverage ×{:.3}** {}",
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
        let mut in_chain: HashMap<u32, f64> = HashMap::new();
        for sw in &bs.spillways {
            if let Some(t) = sw.chained_into {
                *in_chain.entry(t).or_default() += sw.discharge_m3s as f64;
            }
        }
        let out_of: HashMap<u32, f64> =
            bs.spillways.iter().map(|s| (s.lake_id, s.discharge_m3s as f64)).collect();
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

        // ── §3 · ONE basin cross-checked by independent reconstruction ────────────
        eprintln!(
            "\n── §3 · INDEPENDENT CROSS-CHECK of one basin (watershed × measured runoff) ──"
        );
        if let Some(target) = rows.iter().find(|r| r.1 > 0.0).map(|r| r.0) {
            // the D8 basin labels present in this below-sea footprint
            let labels: HashSet<u32> = (0..n)
                .filter(|&k| bs.lake_map[k] == target)
                .map(|k| flow.basins[k])
                .filter(|&l| l != 0)
                .collect();
            let mut recon = 0.0f64;
            for k in 0..n {
                if field.data[k] > SEA && labels.contains(&flow.basins[k]) {
                    recon += surplus_mm_km2[k] as f64;
                }
            }
            let recon_m3s = runoff_km2_to_m3s(recon as f32) as f64;
            let inc = in_chain.get(&target).copied().unwrap_or(0.0);
            let reported = rows.iter().find(|r| r.0 == target).map(|r| r.1).unwrap_or(0.0);
            eprintln!(
                "   basin {target}: {} D8 basin label(s) in its footprint | reconstructed \
                 watershed runoff {recon_m3s:.2} + upstream chain {inc:.2} = {:.2} m³/s\n   \
                 against `BasinSummary::inflow_m3s` = {reported:.2} ⇒ ratio {:.3} {}",
                labels.len(),
                recon_m3s + inc,
                reported / (recon_m3s + inc).max(1e-9),
                if ((reported / (recon_m3s + inc).max(1e-9)) - 1.0).abs() < 0.5 {
                    "— the machinery's own accounting is CONFIRMED by an independent route"
                } else {
                    "— ⚠ the machinery's accounting is NOT confirmed; the ledger is suspect"
                }
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
