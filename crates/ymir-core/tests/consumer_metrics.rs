//! Consumer-oriented measurements: river widths in render cells, and buildable connected
//! components. Diagnostic only — nothing is asserted about the terrain, nothing is changed.
//!
//! Target configuration throughout: **8192², 400 km domain, 48.83 m cell**, both climate beds,
//! `geo_scale_ratio` at 1.0 **and** 7.5. The ratio is a post-process on the drainage result
//! (ADR Finding 68), so both values are obtained from ONE drainage by applying
//! `apply_geo_scale_ratio` to a clone — which is exactly what production does, and costs nothing.
//!
//! Run: cargo test -p ymir-core --release --test consumer_metrics -- --ignored --nocapture

use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, DrainageClimate, SegmentKind, apply_geo_scale_ratio, c1_drainage_windowed,
};
use ymir_core::tectonics_c1::hd_assembly::assemble_hd_drainage;
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::{
    c1_altitude_norm_to_metres, upscale_from_c1_with_progress,
};
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::flow::breach_monotone;
use ymir_core::terrain::upscale::{ProductionHdOpts, production_hd_config};

const PSEED: u64 = 10_481_999_410_520_546_993;
const DOMAIN_KM: f32 = 400.0;
const SEA: f32 = 0.5;
const TARGET: usize = 8192;
/// The render cell, in metres — the unit the consumer draws in.
const CELL_M: f32 = DOMAIN_KM / TARGET as f32 * 1000.0;
/// The two climate beds. The terrain is climate-independent; the hydrology is not.
const BEDS: [(&str, f32, f32); 2] = [("arid-hot", 25.0, 10.0), ("humid", 45.0, 40.0)];

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

/// The breached field — what production exports and what a consumer renders.
fn exported_field(raw: &GridF32, ss: &SteinSteinParams) -> GridF32 {
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let pre = c1_drainage_windowed(raw, None, &dcfg, ss, DOMAIN_KM);
    breach_monotone(raw, &pre.flow.filled, &pre.lake_map, SEA, raw.width, raw.height)
}

fn q(v: &mut Vec<f32>, f: f64) -> f32 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    if v.is_empty() {
        return f32::NAN;
    }
    v[(((v.len() - 1) as f64 * f).round() as usize).min(v.len() - 1)]
}

/// MEASURE 1 — river width in RENDER CELLS. Is the hierarchy legible?
#[test]
#[ignore]
fn river_width_in_render_cells() {
    let ss = SteinSteinParams::default();
    eprintln!(
        "\n==========  1 · river width in RENDER CELLS ({TARGET}², {CELL_M:.2} m/cell)  =========="
    );
    eprintln!(
        "  width law: w = CHANNEL_WIDTH_A · Q^CHANNEL_WIDTH_B = 5.0 · Q^0.5 (Leopold & Maddock;\n  \
         the EXPONENT has antecedence, the COEFFICIENT drifted 1.2 → 5.0 with no finding — PROXY)"
    );
    let raw = terrain();
    let field = exported_field(&raw, &ss);
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;

    for (bed, lat, span) in BEDS {
        let climate =
            c1_climate_placed(&field, &ss, lat, span, &PrecipParams::default(), DOMAIN_KM);
        let dclim = DrainageClimate {
            precip_internal: &climate.precipitation,
            temperature: &climate.temperature,
        };
        let pre = c1_drainage_windowed(&raw, None, &dcfg, &ss, DOMAIN_KM);
        let base = assemble_hd_drainage(
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
        for ratio in [1.0f32, 7.5] {
            let mut dr = base.clone();
            apply_geo_scale_ratio(&mut dr, ratio, &dcfg.thresholds);
            eprintln!("\n──────── {bed} · geo_scale_ratio {ratio} ────────");
            // Watercourses only: a spillway is not a rendered river hierarchy member.
            let idx: Vec<usize> = (0..dr.rivers.segments.len())
                .filter(|&i| dr.segment_kind[i] != SegmentKind::Spillway)
                .collect();
            let cells: Vec<f32> = idx.iter().map(|&i| dr.segment_width_m[i] / CELL_M).collect();
            let lens: Vec<f32> =
                idx.iter().map(|&i| dr.rivers.segments[i].points.len() as f32).collect();
            let tot_len: f32 = lens.iter().sum();
            let bucket = |lo: f32, hi: f32| -> (f32, f32) {
                let (mut n, mut l) = (0.0f32, 0.0f32);
                for (c, w) in cells.iter().zip(lens.iter()) {
                    if *c >= lo && *c < hi {
                        n += 1.0;
                        l += *w;
                    }
                }
                (100.0 * n / cells.len().max(1) as f32, 100.0 * l / tot_len.max(1e-6))
            };
            for (lab, lo, hi) in [
                ("< 1 cell", 0.0, 1.0),
                ("1–2 cells", 1.0, 2.0),
                ("2–5 cells", 2.0, 5.0),
                ("> 5 cells", 5.0, f32::INFINITY),
            ] {
                let (bn, bl) = bucket(lo, hi);
                eprintln!("  {lab:<10}  {bn:>6.2} % of segments | {bl:>6.2} % of length");
            }
            eprintln!(
                "  cells: p50 {:.3} | p90 {:.3} | p99 {:.3} | max {:.3}   ({} watercourses)",
                q(&mut cells.clone(), 0.50),
                q(&mut cells.clone(), 0.90),
                q(&mut cells.clone(), 0.99),
                q(&mut cells.clone(), 1.00),
                cells.len()
            );
            // per Strahler order
            eprint!("  by Strahler order (p50 cells):");
            for o in 1u8..=6 {
                let mut v: Vec<f32> = idx
                    .iter()
                    .filter(|&&i| dr.rivers.segments[i].strahler_order == o)
                    .map(|&i| dr.segment_width_m[i] / CELL_M)
                    .collect();
                if !v.is_empty() {
                    eprint!("  S{o}:{:.3} (n {})", q(&mut v, 0.5), v.len());
                }
            }
            eprintln!();
            eprint!("  by Strahler order (max cells):");
            for o in 1u8..=6 {
                let mut v: Vec<f32> = idx
                    .iter()
                    .filter(|&&i| dr.rivers.segments[i].strahler_order == o)
                    .map(|&i| dr.segment_width_m[i] / CELL_M)
                    .collect();
                if !v.is_empty() {
                    eprint!("  S{o}:{:.2}", q(&mut v, 1.0));
                }
            }
            eprintln!();
        }
    }
    eprintln!(
        "\n  CONTROL: the 1.0 → 7.5 shift must be EXACTLY 7.5 (w = a·Q^0.5 with Q × 56.25).\n  \
         If it is not, the shipped law is not `5.0 · Q^0.5`."
    );
}

/// MEASURE 2 — buildable CONNECTED COMPONENTS. Connectivity, not flatness.
///
/// ## The extent measure, declared because the verdict depends on it entirely
///
/// A long thin component of 1 200 cells cannot host a city. Extent is therefore the **diameter of
/// the largest disc that fits inside the component**: an 8-connected chamfer distance transform
/// (weights 1 and √2) from the non-buildable set, then `2 × max` within each component. A ribbon
/// one cell wide scores 2 cells however long it is, which is the intended behaviour.
///
/// TDD §2.2 gives these as VALLEY WIDTHS, which is the same quantity: ravine 50–100 m, hameau
/// 150–250 m, village/bourg 300–600 m, ville 800 m–1.2 km, cité 2 km+.
#[test]
#[ignore]
fn buildable_connected_components() {
    let ss = SteinSteinParams::default();
    eprintln!("\n==========  2 · buildable connected components (8-connectivity)  ==========");
    eprintln!(
        "  criterion: LOCAL SLOPE below a threshold, no constraint on internal relief — a\n  \
         terraced city is acceptable (ADR line 28: the criterion is a LEGIBLE landscape, NOT\n  \
         minimal steep ground). Extent = diameter of the largest inscribed disc."
    );
    let raw = terrain();
    let field = exported_field(&raw, &ss);
    let (w, h) = (field.width, field.height);
    let to_m = |x: f32| c1_altitude_norm_to_metres(x, &ss);

    // dimensionless local slope, central differences, physical units
    let mut slope = vec![0.0f32; w * h];
    for y in 1..h - 1 {
        for x in 1..w - 1 {
            let k = y * w + x;
            if field.data[k] <= SEA {
                continue;
            }
            let gx = 0.5 * (to_m(field.data[k + 1]) - to_m(field.data[k - 1]));
            let gy = 0.5 * (to_m(field.data[k + w]) - to_m(field.data[k - w]));
            slope[k] = (gx * gx + gy * gy).sqrt() / CELL_M;
        }
    }
    let land: usize = field.data.iter().filter(|&&x| x > SEA).count();
    // distance to the coast, in cells, for the altitude/coast cross-tab: cheap proxy via a
    // chamfer transform from the sea.
    let coast_dist = chamfer(&field.data.iter().map(|&x| x > SEA).collect::<Vec<bool>>(), w, h);

    eprintln!("\n  land {land} cells = {:.0} km²", land as f32 * (CELL_M / 1000.0).powi(2));
    for thr in [0.10f32, 0.15, 0.20, 0.25] {
        let mask: Vec<bool> = (0..w * h).map(|k| field.data[k] > SEA && slope[k] < thr).collect();
        let n_build = mask.iter().filter(|&&b| b).count();
        // inscribed radius everywhere inside the mask
        let dist = chamfer(&mask, w, h);
        // 8-connected components, iterative flood fill
        let mut label = vec![0u32; w * h];
        let mut sizes: Vec<u32> = vec![0];
        let mut radii: Vec<f32> = vec![0.0];
        let mut alt_sum: Vec<f64> = vec![0.0];
        let mut coast_min: Vec<f32> = vec![0.0];
        let mut relief: Vec<(f32, f32)> = vec![(0.0, 0.0)];
        let mut next = 1u32;
        let mut stack: Vec<usize> = Vec::new();
        for s in 0..w * h {
            if !mask[s] || label[s] != 0 {
                continue;
            }
            let id = next;
            next += 1;
            sizes.push(0);
            radii.push(0.0);
            alt_sum.push(0.0);
            coast_min.push(f32::INFINITY);
            relief.push((f32::INFINITY, f32::NEG_INFINITY));
            label[s] = id;
            stack.push(s);
            while let Some(k) = stack.pop() {
                let i = id as usize;
                sizes[i] += 1;
                let a = to_m(field.data[k]);
                alt_sum[i] += a as f64;
                if dist[k] > radii[i] {
                    radii[i] = dist[k];
                }
                if coast_dist[k] < coast_min[i] {
                    coast_min[i] = coast_dist[k];
                }
                relief[i].0 = relief[i].0.min(a);
                relief[i].1 = relief[i].1.max(a);
                let (x, y) = (k % w, k / w);
                for dy in -1i32..=1 {
                    for dx in -1i32..=1 {
                        if dx == 0 && dy == 0 {
                            continue;
                        }
                        let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                        if nx < 0 || ny < 0 || nx as usize >= w || ny as usize >= h {
                            continue;
                        }
                        let nk = ny as usize * w + nx as usize;
                        if mask[nk] && label[nk] == 0 {
                            label[nk] = id;
                            stack.push(nk);
                        }
                    }
                }
            }
        }
        let ncomp = (next - 1) as usize;
        let biggest = sizes.iter().skip(1).copied().max().unwrap_or(0);
        let cell_km2 = (CELL_M / 1000.0).powi(2);
        eprintln!(
            "\n  ── slope < {:.0} % ({:.2}°) ── buildable {n_build} cells = {:.1} % of land | \
             {ncomp} components",
            100.0 * thr,
            thr.atan().to_degrees(),
            100.0 * n_build as f32 / land.max(1) as f32
        );
        eprintln!(
            "     LARGEST component {biggest} cells = {:.0} km² = {:.1} % of the buildable land",
            biggest as f32 * cell_km2,
            100.0 * biggest as f32 / n_build.max(1) as f32
        );
        let mut sz: Vec<f32> = sizes.iter().skip(1).map(|&x| x as f32).collect();
        eprintln!(
            "     component size, cells: p50 {:.0} | p90 {:.0} | p99 {:.0} | max {:.0}",
            q(&mut sz.clone(), 0.5),
            q(&mut sz.clone(), 0.9),
            q(&mut sz.clone(), 0.99),
            q(&mut sz.clone(), 1.0)
        );
        // TDD §2.2 classes, on the INSCRIBED DIAMETER
        eprintln!("     components by inscribed diameter (TDD §2.2 valley-width classes):");
        for (lab, lo_m) in [
            ("ravine   ≥  50 m", 50.0f32),
            ("hameau   ≥ 150 m", 150.0),
            ("village  ≥ 300 m", 300.0),
            ("ville    ≥ 800 m", 800.0),
            ("cité     ≥ 2000 m", 2000.0),
        ] {
            let n = (1..=ncomp).filter(|&i| 2.0 * radii[i] * CELL_M >= lo_m).count();
            eprintln!("       {lab}  {n:>6}");
        }
        // the cross-tabs, on the components that clear `ville`
        let mut big: Vec<usize> =
            (1..=ncomp).filter(|&i| 2.0 * radii[i] * CELL_M >= 800.0).collect();
        big.sort_by(|&a, &b| sizes[b].cmp(&sizes[a]));
        eprintln!(
            "     the {} components ≥ ville scale, largest first (top 8):\n       {:>10} {:>10} \
             {:>10} {:>10} {:>12} {:>10}",
            big.len(),
            "cells",
            "km²",
            "⌀ m",
            "alt p50 m",
            "coast km",
            "relief m"
        );
        for &i in big.iter().take(8) {
            eprintln!(
                "       {:>10} {:>10.0} {:>10.0} {:>10.0} {:>12.1} {:>10.0}",
                sizes[i],
                sizes[i] as f32 * cell_km2,
                2.0 * radii[i] * CELL_M,
                alt_sum[i] / sizes[i].max(1) as f64,
                coast_min[i] * CELL_M / 1000.0,
                relief[i].1 - relief[i].0
            );
        }
    }
    eprintln!(
        "\n  CONTROL: the component count and the largest-component share must MOVE across the\n  \
         threshold sweep. A constant column measures the instrument (rule 10)."
    );
}

/// 8-connected chamfer distance (weights 1, √2) from the complement of `mask`, in CELLS.
/// Cells outside the mask get 0; a cell deep inside gets its distance to the nearest boundary,
/// which is the radius of the largest disc centred there that fits inside the mask.
fn chamfer(mask: &[bool], w: usize, h: usize) -> Vec<f32> {
    const D: f32 = std::f32::consts::SQRT_2;
    let big = (w + h) as f32;
    let mut d: Vec<f32> = mask.iter().map(|&b| if b { big } else { 0.0 }).collect();
    for y in 0..h {
        for x in 0..w {
            let k = y * w + x;
            if d[k] == 0.0 {
                continue;
            }
            let mut m = d[k];
            if y > 0 {
                m = m.min(d[k - w] + 1.0);
                if x > 0 {
                    m = m.min(d[k - w - 1] + D);
                }
                if x + 1 < w {
                    m = m.min(d[k - w + 1] + D);
                }
            }
            if x > 0 {
                m = m.min(d[k - 1] + 1.0);
            }
            d[k] = m;
        }
    }
    for y in (0..h).rev() {
        for x in (0..w).rev() {
            let k = y * w + x;
            if d[k] == 0.0 {
                continue;
            }
            let mut m = d[k];
            if y + 1 < h {
                m = m.min(d[k + w] + 1.0);
                if x > 0 {
                    m = m.min(d[k + w - 1] + D);
                }
                if x + 1 < w {
                    m = m.min(d[k + w + 1] + D);
                }
            }
            if x + 1 < w {
                m = m.min(d[k + 1] + 1.0);
            }
            d[k] = m;
        }
    }
    d
}
