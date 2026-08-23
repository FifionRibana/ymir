//! DIAGNOSTIC (read-only) — terraces-following-isolines + missing dendritic
//! valleys. Dumps altitude transects + statistics at each pipeline stage to locate
//! where banding first appears, quantifies the terraces, checks whether erosion
//! carves along the drainage network, runs erosion on a synthetic cone, and
//! compares two resolutions. Reports numbers only; changes NOTHING in the pipeline.
//!
//! Run: cargo test -p ymir-core --test terrace_diagnosis --release -- --ignored --nocapture

use ymir_core::erosion::hydraulic::{ErosionConfig, run_erosion};
use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage};
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::{
    C1_DOMAIN_KM, c1_altitude_norm_to_metres, c1_coarse_normalized_altitude,
};
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::upscale::{FbmUpscaleConfig, upscale_with_fbm};

const GRID: usize = 64;
const SEED: u64 = 42;
const SEA: f32 = 0.5;

fn to_metres(norm: &GridF32, ss: &SteinSteinParams) -> GridF32 {
    GridF32::from_vec(
        norm.width,
        norm.height,
        norm.data.iter().map(|&n| c1_altitude_norm_to_metres(n, ss)).collect(),
    )
}

fn coarse_state(seed: u64) -> (ymir_core::tectonics_c1::state::C1State, C1TimeLoopConfig) {
    let run = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / GRID as f64,
        dy: 1.0 / GRID as f64,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let mut state = init_c1_state_phase_2_r7(GRID, seed, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, &run, &C1Closures::default(), |_, _| {});
    (state, run)
}

/// eroded (from field diff) — the actual carved vs deposited mass (norm units).
fn erosion_balance(before: &GridF32, after: &GridF32) -> (f64, f64) {
    let (mut eroded, mut deposited) = (0.0f64, 0.0f64);
    for (a, b) in before.data.iter().zip(after.data.iter()) {
        let d = (*b - *a) as f64;
        if d < 0.0 {
            eroded += -d;
        } else {
            deposited += d;
        }
    }
    (eroded, deposited)
}

/// Transect stats: min/max (m), step count (|Δ| > 30 m between adjacent cells),
/// plateau widths (runs where |Δ| < 2 m), and whether the LARGE steps cluster at
/// discrete sizes. Land cells only for the altitude histogram.
fn transect_report(label: &str, row_m: &[f32]) {
    let land: Vec<f32> = row_m.iter().copied().filter(|&m| m > 0.0).collect();
    let (mn, mx) = row_m.iter().fold((f32::MAX, f32::MIN), |(a, b), &v| (a.min(v), b.max(v)));
    let diffs: Vec<f32> = row_m.windows(2).map(|w| (w[1] - w[0]).abs()).collect();
    let steps: Vec<f32> = diffs.iter().copied().filter(|&d| d > 30.0).collect();
    // plateau widths (cells): runs where the surface is ~flat (|Δ| < 2 m).
    let mut plateaus = Vec::new();
    let mut run = 1usize;
    for &d in &diffs {
        if d < 2.0 {
            run += 1;
        } else {
            if run >= 3 {
                plateaus.push(run);
            }
            run = 1;
        }
    }
    if run >= 3 {
        plateaus.push(run);
    }
    plateaus.sort_unstable();
    let med_plat = plateaus.get(plateaus.len() / 2).copied().unwrap_or(0);
    let mut steps_sorted = steps.clone();
    steps_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let modal_step = steps_sorted.get(steps_sorted.len() / 2).copied().unwrap_or(0.0);
    eprintln!(
        "  [{label:<22}] range {mn:>7.0}..{mx:>6.0} m (land n={})  \
         steps>30m: {:>4}  modal step {:>5.0} m  plateaus(≥3cell): {:>4}  median plateau {:>4} cells",
        land.len(),
        steps.len(),
        modal_step,
        plateaus.len(),
        med_plat,
    );
}

/// Altitude-VALUE histogram peakiness: are land altitudes clustered at discrete
/// levels (terraced) or spread continuously? Reports the number of populated bins
/// and the ratio of the top bin to the mean bin (a terraced field is peaky).
fn altitude_modality(label: &str, field_m: &GridF32) {
    let land: Vec<f32> = field_m.data.iter().copied().filter(|&m| m > 0.0).collect();
    if land.len() < 100 {
        eprintln!("  [{label}] <100 land cells, skip modality");
        return;
    }
    let (mn, mx) = land.iter().fold((f32::MAX, f32::MIN), |(a, b), &v| (a.min(v), b.max(v)));
    const NB: usize = 120;
    let mut bins = vec![0u32; NB];
    let span = (mx - mn).max(1e-3);
    for &m in &land {
        let b = (((m - mn) / span) * (NB as f32 - 1.0)) as usize;
        bins[b.min(NB - 1)] += 1;
    }
    let populated = bins.iter().filter(|&&c| c > 0).count();
    let maxb = *bins.iter().max().unwrap();
    let mean = land.len() as f32 / populated as f32;
    // count "spikes": bins with > 3× the mean population (discrete levels).
    let spikes = bins.iter().filter(|&&c| c as f32 > 3.0 * mean).count();
    eprintln!(
        "  [{label:<22}] altitude bins populated {populated}/{NB}, peak/mean {:.1}, spikes(>3×mean) {spikes} \
         → {}",
        maxb as f32 / mean,
        if spikes > 3 { "DISCRETE LEVELS (terraced)" } else { "continuous" },
    );
}

fn max_row(field: &GridF32) -> usize {
    let (mut best, mut bi) = (f32::MIN, 0usize);
    for (k, &v) in field.data.iter().enumerate() {
        if v > best {
            best = v;
            bi = k;
        }
    }
    bi / field.width
}

fn row_of(field: &GridF32, r: usize) -> Vec<f32> {
    let w = field.width;
    field.data[r * w..(r + 1) * w].to_vec()
}

fn ero_cfg(target: usize) -> ErosionConfig {
    FbmUpscaleConfig::c1_hd_production(target).erosion.unwrap()
}

fn stages(target: usize, ss: &SteinSteinParams) -> (GridF32, GridF32, GridF32, f64, f64) {
    let (state, _run) = coarse_state(SEED);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), ss, None);
    let seed = WorldSeed::new(SEED);
    let mut cfg = FbmUpscaleConfig::c1_hd_production(target);
    cfg.erosion = None;
    cfg.bathymetry = None;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &cfg).heightmap;
    let ero = ero_cfg(target);
    let eroded = run_erosion(&fbm, &ero, &seed, |_, _, _| true).heightmap;
    let (e, d) = erosion_balance(&fbm, &eroded);
    (coarse, fbm, eroded, e, d)
}

#[test]
#[ignore]
fn diagnose_terraces_and_valleys() {
    let ss = SteinSteinParams::default();
    let target = 2048usize;

    eprintln!("\n=== STEP 1 — where banding first appears (seed {SEED}, {target}²) ===");
    let (coarse, fbm, eroded, e_mass, d_mass) = stages(target, &ss);
    let coarse_m = to_metres(&coarse, &ss);
    let fbm_m = to_metres(&fbm, &ss);
    let eroded_m = to_metres(&eroded, &ss);

    let r = max_row(&eroded);
    let rc = (r * GRID / target).min(GRID - 1);
    eprintln!("(transect: eroded row {r}; coarse row {rc})");
    transect_report("a coarse (post-iso)", &row_of(&coarse_m, rc));
    transect_report("b FBM pre-erosion", &row_of(&fbm_m, r));
    transect_report("c after erosion", &row_of(&eroded_m, r));

    eprintln!("--- altitude-value modality (whole field) ---");
    altitude_modality("a coarse (post-iso)", &coarse_m);
    altitude_modality("b FBM pre-erosion", &fbm_m);
    altitude_modality("c after erosion", &eroded_m);

    // STEP 2 — u16 quantisation refutation.
    let (mn, mx) =
        eroded_m.data.iter().fold((f32::MAX, f32::MIN), |(a, b), &v| (a.min(v), b.max(v)));
    let mpu = (mx - mn) / 65535.0;
    eprintln!(
        "\n=== STEP 2 — quantisation ===\n  eroded range {mn:.1}..{mx:.1} m → u16 step = {mpu:.4} m/unit \
         (author's jumps are 100s of m → u16 quantisation REFUTED as the cause)"
    );

    // STEP 3 — erosion effect + flow-vs-altitude.
    let ero = ero_cfg(target);
    eprintln!(
        "\n=== STEP 3 — erosion ({target}²) ===\n  params: droplets {} ({:.2}/cell), erosion_rate {}, \
         deposition_rate {}, inertia {}, gravity {}, min_slope {}, lifetime {}, radius {}",
        ero.num_droplets,
        ero.num_droplets as f64 / (target * target) as f64,
        ero.erosion_rate,
        ero.deposition_rate,
        ero.inertia,
        ero.gravity,
        ero.min_slope,
        ero.max_lifetime,
        ero.erosion_radius,
    );
    let e_m = e_mass * 11300.0 / (target * target) as f64; // norm→m, per cell avg
    let d_m = d_mass * 11300.0 / (target * target) as f64;
    eprintln!(
        "  mass moved (norm): eroded {e_mass:.1}, deposited {d_mass:.1}, net {:.1}% \
         (avg |Δ|: eroded {e_m:.2} m/cell, deposited {d_m:.2} m/cell)",
        (e_mass - d_mass).abs() / e_mass.max(1.0) * 100.0,
    );

    // Drainage on the eroded field → accumulation vs local altitude.
    let dr = c1_drainage(&eroded, None, &C1DrainageConfig::default(), &ss);
    let acc = &dr.flow.accumulation;
    let (w, h) = (eroded.width, eroded.height);
    // For high-accumulation cells (top 1%), is the cell a local altitude minimum
    // (carved channel) vs its 8 neighbours? Report the fraction.
    let mut accs: Vec<f32> = acc.data.clone();
    accs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let thr = accs[(accs.len() as f64 * 0.99) as usize];
    let (mut hi, mut in_min) = (0usize, 0usize);
    for y in 1..h - 1 {
        for x in 1..w - 1 {
            let k = y * w + x;
            if acc.data[k] < thr || eroded.data[k] <= SEA {
                continue;
            }
            hi += 1;
            let c = eroded.data[k];
            let mut is_min = true;
            for (dx, dy) in
                [(-1i32, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)]
            {
                let nk = ((y as i32 + dy) as usize) * w + (x as i32 + dx) as usize;
                if eroded.data[nk] < c {
                    is_min = false;
                    break;
                }
            }
            if is_min {
                in_min += 1;
            }
        }
    }
    eprintln!(
        "  high-accumulation land cells (top 1%): {hi}; of those {in_min} ({:.0}%) are LOCAL MINIMA \
         (in a carved channel). Low % ⇒ water flows over uncarved terrain.",
        in_min as f32 / hi.max(1) as f32 * 100.0,
    );
    eprintln!(
        "  rivers: {} segments, max Strahler {}",
        dr.rivers.segments.len(),
        dr.rivers.segments.iter().map(|s| s.strahler_order).max().unwrap_or(0),
    );

    // Cross-section perpendicular to the highest-Strahler river channel.
    if let Some(seg) = dr.rivers.segments.iter().max_by_key(|s| s.strahler_order) {
        if let Some(&(px, py)) = seg.points.get(seg.points.len() / 2) {
            let (px, py) = (px as usize, py as usize);
            eprint!(
                "  cross-section ±12 cells across a Strahler-{} channel (m): ",
                seg.strahler_order
            );
            let mut xs = Vec::new();
            for o in -12i32..=12 {
                let x = (px as i32 + o).clamp(0, w as i32 - 1);
                xs.push(eroded_m.get(x, py as i32));
            }
            let cmin = xs.iter().cloned().fold(f32::MAX, f32::min);
            for v in &xs {
                eprint!("{:.0} ", v);
            }
            eprintln!(
                "\n    channel depth below rim: {:.0} m (V/U valley if clearly incised)",
                xs[0].max(xs[24]) - cmin
            );
        }
    }

    // STEP 4 — synthetic smooth cone.
    eprintln!("\n=== STEP 4 — erosion on a synthetic smooth cone (1024²) ===");
    let n = 1024usize;
    let mut cone = GridF32::new(n, n, 0.0);
    let (cx, cy) = (n as f32 / 2.0, n as f32 / 2.0);
    let rmax = n as f32 / 2.0;
    for y in 0..n {
        for x in 0..n {
            let r = ((x as f32 - cx).powi(2) + (y as f32 - cy).powi(2)).sqrt();
            // Smooth cone 0.5..0.8 norm (land), tiny deterministic ripple to seed flow.
            let base = 0.5 + 0.30 * (1.0 - (r / rmax).min(1.0));
            let ripple = 0.002 * ((x as f32 * 0.3).sin() + (y as f32 * 0.3).cos());
            cone.set(x, y, base + ripple);
        }
    }
    let cseed = WorldSeed::new(7);
    let cone_ero = run_erosion(&cone, &ero_cfg(n), &cseed, |_, _, _| true).heightmap;
    let (ce, cd) = erosion_balance(&cone, &cone_ero);
    let dr0 = c1_drainage(&cone, None, &C1DrainageConfig::default(), &ss);
    let dr1 = c1_drainage(&cone_ero, None, &C1DrainageConfig::default(), &ss);
    let maxstr = |d: &ymir_core::tectonics_c1::drainage::C1DrainageResult| {
        d.rivers.segments.iter().map(|s| s.strahler_order).max().unwrap_or(0)
    };
    eprintln!(
        "  cone eroded {ce:.1} / deposited {cd:.1} (norm); rivers before {} (maxS {}) → after {} (maxS {})",
        dr0.rivers.segments.len(),
        maxstr(&dr0),
        dr1.rivers.segments.len(),
        maxstr(&dr1),
    );
    let std_before = std_of(&cone);
    let std_after = std_of(&cone_ero);
    eprintln!(
        "  surface roughness (std of Δ to 4-neighbour mean, norm): before {std_before:.5} → after {std_after:.5} \
         → {}",
        if std_after > std_before * 1.3 {
            "erosion CARVED channels"
        } else {
            "little/no channelisation"
        },
    );

    // STEP 5 — resolution comparison (plateau width in cells vs metres).
    eprintln!("\n=== STEP 5 — resolution comparison (plateau width) ===");
    for t in [512usize, 2048usize] {
        let (_c, _f, er, _e, _d) = stages(t, &ss);
        let erm = to_metres(&er, &ss);
        let r = max_row(&er);
        let km_per_cell = C1_DOMAIN_KM / t as f32; // full-domain scale
        eprint!("  {t}²: ");
        let row = row_of(&erm, r);
        transect_report(&format!("{t}² eroded"), &row);
        eprintln!("        (km/cell {:.3} at full domain)", km_per_cell);
    }
}

/// A proper TERRACE detector on a transect (metres): count flat plateaux (runs of
/// ≥ 5 cells with |Δ| < 5 m) that are bounded by a sharp step (|Δ| > 50 m), and the
/// fraction of the land transect that is "flat" (|Δ| < 5 m). Terraced terrain has a
/// high flat fraction punctuated by big steps; a naturally rough fractal surface
/// has a low flat fraction.
fn terrace_report(label: &str, row_m: &[f32]) {
    let n = row_m.len();
    let diffs: Vec<f32> = row_m.windows(2).map(|w| w[1] - w[0]).collect();
    let land: Vec<bool> = row_m.iter().map(|&m| m > 0.0).collect();
    let mut flat_cells = 0usize;
    let mut land_cells = 0usize;
    for i in 0..diffs.len() {
        if land[i] && land[i + 1] {
            land_cells += 1;
            if diffs[i].abs() < 5.0 {
                flat_cells += 1;
            }
        }
    }
    // Terraces: flat runs ≥5 cells ended by a >50 m step.
    let mut terraces = 0usize;
    let mut run = 1usize;
    let mut widths = Vec::new();
    for i in 0..diffs.len() {
        if diffs[i].abs() < 5.0 {
            run += 1;
        } else {
            if run >= 5 && diffs[i].abs() > 50.0 {
                terraces += 1;
                widths.push(run);
            }
            run = 1;
        }
    }
    widths.sort_unstable();
    let medw = widths.get(widths.len() / 2).copied().unwrap_or(0);
    eprintln!(
        "  [{label:<24}] flat fraction {:>4.0}% of land ({flat_cells}/{land_cells}), \
         terraces(≥5-cell flat + >50m step): {terraces:>3}, median terrace width {medw} cells (of {n})",
        flat_cells as f32 / land_cells.max(1) as f32 * 100.0,
    );
}

#[test]
#[ignore]
fn disentangle_terrace_source() {
    let ss = SteinSteinParams::default();
    let t = 1024usize;
    let (state, _run) = coarse_state(SEED);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(SEED);

    // (i) PURE BILINEAR upscale: FBM amplitude 0 (no noise), no erosion. Isolates
    //     whether the coarse field alone (H3) produces terraces after interpolation.
    let mut bil = FbmUpscaleConfig::c1_hd_production(t);
    bil.erosion = None;
    bil.bathymetry = None;
    bil.amplitude_base = 0.0; // kill the FBM detail entirely
    bil.coast_warp_strength = 0.0;
    let f_bil = upscale_with_fbm(&coarse, SEA, &seed, &bil).heightmap;

    // (ii) FBM full (production amplitude), no erosion. Adds H2.
    let mut fb = FbmUpscaleConfig::c1_hd_production(t);
    fb.erosion = None;
    fb.bathymetry = None;
    let f_fbm = upscale_with_fbm(&coarse, SEA, &seed, &fb).heightmap;

    // (iii) eroded (production). Adds H4.
    let f_ero = run_erosion(&f_fbm, &ero_cfg(t), &seed, |_, _, _| true).heightmap;

    let r = max_row(&f_ero);
    eprintln!("\n=== terrace source disentangle ({t}², transect row {r}) ===");
    terrace_report("i  pure bilinear (H3)", &row_of(&to_metres(&f_bil, &ss), r));
    terrace_report("ii FBM, no erosion (H2)", &row_of(&to_metres(&f_fbm, &ss), r));
    terrace_report("iii eroded (H4)", &row_of(&to_metres(&f_ero, &ss), r));
}

/// Fraction of top-1%-accumulation LAND cells that are local altitude minima
/// (i.e. sit in a carved channel). Higher = erosion incised along the flow.
fn carved_fraction(field: &GridF32, ss: &SteinSteinParams) -> (f32, u8) {
    let dr = c1_drainage(field, None, &C1DrainageConfig::default(), ss);
    let acc = &dr.flow.accumulation;
    let (w, h) = (field.width, field.height);
    let mut accs = acc.data.clone();
    accs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let thr = accs[(accs.len() as f64 * 0.99) as usize];
    let (mut hi, mut minc) = (0usize, 0usize);
    for y in 1..h - 1 {
        for x in 1..w - 1 {
            let k = y * w + x;
            if acc.data[k] < thr || field.data[k] <= SEA {
                continue;
            }
            hi += 1;
            let c = field.data[k];
            let mut is_min = true;
            for (dx, dy) in [(-1i32, 0), (1, 0), (0, -1), (0, 1)] {
                let nk = ((y as i32 + dy) as usize) * w + (x as i32 + dx) as usize;
                if field.data[nk] < c {
                    is_min = false;
                    break;
                }
            }
            if is_min {
                minc += 1;
            }
        }
    }
    let maxstr = dr.rivers.segments.iter().map(|s| s.strahler_order).max().unwrap_or(0);
    (minc as f32 / hi.max(1) as f32, maxstr)
}

/// Read-only VALIDATION of the fix proposal: sweep droplet density × erosion/
/// deposition balance on the synthetic smooth cone and measure channelisation
/// (roughness increase + carved-channel fraction). Confirms which parameters
/// actually incise. Changes NO pipeline default — the config is local to the test.
#[test]
#[ignore]
fn erosion_param_sweep_cone() {
    let ss = SteinSteinParams::default();
    let n = 512usize;
    let mut cone = GridF32::new(n, n, 0.0);
    let (cx, cy, rmax) = (n as f32 / 2.0, n as f32 / 2.0, n as f32 / 2.0);
    for y in 0..n {
        for x in 0..n {
            let r = ((x as f32 - cx).powi(2) + (y as f32 - cy).powi(2)).sqrt();
            let base = 0.5 + 0.30 * (1.0 - (r / rmax).min(1.0));
            let ripple = 0.002 * ((x as f32 * 0.3).sin() + (y as f32 * 0.3).cos());
            cone.set(x, y, base + ripple);
        }
    }
    let rough0 = std_of(&cone);
    let (carved0, _) = carved_fraction(&cone, &ss);
    let cells = (n * n) as f64;
    eprintln!(
        "\n=== erosion param sweep on a smooth cone ({n}²) — baseline roughness {rough0:.5}, carved {:.0}% ===",
        carved0 * 100.0
    );
    eprintln!("  density/cell | ero/dep     | net% | roughness→   | carved% | maxStrahler");
    let base = FbmUpscaleConfig::c1_hd_production(n).erosion.unwrap();
    // (droplets/cell, erosion_rate, deposition_rate)
    let configs = [
        (base.num_droplets as f64 / cells, base.erosion_rate, base.deposition_rate), // production
        (4.0, base.erosion_rate, base.deposition_rate),
        (8.0, base.erosion_rate, base.deposition_rate),
        (4.0, 0.6, 0.15),
        (8.0, 0.6, 0.15),
        (8.0, 0.7, 0.05),
    ];
    for (dens, er, dep) in configs {
        let mut cfg = base.clone();
        cfg.num_droplets = (dens * cells) as usize;
        cfg.erosion_rate = er;
        cfg.deposition_rate = dep;
        let seed = WorldSeed::new(7);
        let out = run_erosion(&cone, &cfg, &seed, |_, _, _| true).heightmap;
        let (e, d) = erosion_balance(&cone, &out);
        let rough = std_of(&out);
        let (carved, maxs) = carved_fraction(&out, &ss);
        eprintln!(
            "  {dens:>10.2} | {er:.2}/{dep:.2}   | {:>4.0} | {rough0:.5}→{rough:.5} (×{:.1}) | {:>5.0}% | {maxs}",
            (e - d) / e.max(1.0) * 100.0,
            rough / rough0.max(1e-9),
            carved * 100.0,
        );
    }
    eprintln!(
        "  (channelisation = roughness ratio ≫1 and carved% rising; production row is the current default)"
    );
}

/// Distance-to-coast field (cells): BFS from all sea cells (≤ sea) over land.
fn coast_distance(field: &GridF32) -> Vec<u32> {
    use std::collections::VecDeque;
    let (w, h) = (field.width, field.height);
    let n = w * h;
    let mut dist = vec![u32::MAX; n];
    let mut q = VecDeque::new();
    for k in 0..n {
        if field.data[k] <= SEA {
            dist[k] = 0;
            q.push_back(k);
        }
    }
    while let Some(k) = q.pop_front() {
        let (x, y) = (k % w, k / w);
        let d = dist[k];
        let push = |nx: usize, ny: usize, dist: &mut Vec<u32>, q: &mut VecDeque<usize>| {
            let nk = ny * w + nx;
            if dist[nk] == u32::MAX {
                dist[nk] = d + 1;
                q.push_back(nk);
            }
        };
        if x + 1 < w {
            push(x + 1, y, &mut dist, &mut q);
        }
        if x > 0 {
            push(x - 1, y, &mut dist, &mut q);
        }
        if y + 1 < h {
            push(x, y + 1, &mut dist, &mut q);
        }
        if y > 0 {
            push(x, y - 1, &mut dist, &mut q);
        }
    }
    dist
}

/// Structure metrics of a drained field: (order histogram, #confluences,
/// carved-fraction of top-1% accumulation cells via 8-neighbours, max order).
fn structure_metrics(field: &GridF32, ss: &SteinSteinParams) -> (Vec<usize>, usize, f32, u8) {
    let dr = c1_drainage(field, None, &C1DrainageConfig::default(), ss);
    let mut hist = vec![0usize; 12];
    let mut confl = 0usize;
    for s in &dr.rivers.segments {
        hist[(s.strahler_order as usize).min(11)] += 1;
        if s.upstream.len() >= 2 {
            confl += 1;
        }
    }
    let maxs = dr.rivers.segments.iter().map(|s| s.strahler_order).max().unwrap_or(0);
    let acc = &dr.flow.accumulation;
    let (w, h) = (field.width, field.height);
    let mut accs = acc.data.clone();
    accs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let thr = accs[(accs.len() as f64 * 0.99) as usize];
    let (mut hi, mut minc) = (0usize, 0usize);
    for y in 1..h - 1 {
        for x in 1..w - 1 {
            let k = y * w + x;
            if acc.data[k] < thr || field.data[k] <= SEA {
                continue;
            }
            hi += 1;
            let c = field.data[k];
            let mut is_min = true;
            for (dx, dy) in
                [(-1i32, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)]
            {
                let nk = ((y as i32 + dy) as usize) * w + (x as i32 + dx) as usize;
                if field.data[nk] < c {
                    is_min = false;
                    break;
                }
            }
            if is_min {
                minc += 1;
            }
        }
    }
    (hist, confl, minc as f32 / hi.max(1) as f32, maxs)
}

/// TASK 3 decider — on the REAL coarse→FBM field (1024²), the sink × density
/// matrix with structure metrics + delta survival. Does sink+density TOGETHER
/// restore incision where sink-alone and density-alone each failed?
#[test]
#[ignore]
fn sink_density_matrix() {
    let ss = SteinSteinParams::default();
    let t = 1024usize;
    let cells = (t * t) as f64;
    let (state, _run) = coarse_state(SEED);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(SEED);
    let mut fcfg = FbmUpscaleConfig::c1_hd_production(t);
    fcfg.erosion = None;
    fcfg.bathymetry = None;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fcfg).heightmap;
    let (_, _, carved_fbm, _) = structure_metrics(&fbm, &ss);
    eprintln!(
        "\n=== TASK 3 — sink × density on the real FBM field ({t}²) ===\n  \
         FBM (no erosion): emerged {:.1}%, carved {:.0}% (the incision floor to beat)",
        emerged_frac(&fbm) * 100.0,
        carved_fbm * 100.0,
    );
    eprintln!("  f    density | net% | carvedΔ (after) | maxS | confl | dep≤5/≤20/≤50 | emerged%");
    for (f, dens) in [(1.0f32, 0.95f64), (0.0, 0.95), (0.25, 4.0), (0.0, 4.0)] {
        let mut ec = ero_cfg(t);
        ec.coastal_deposit_fraction = f;
        ec.num_droplets = (dens * cells) as usize;
        let out = run_erosion(&fbm, &ec, &seed, |_, _, _| true).heightmap;
        let (e, d) = erosion_balance(&fbm, &out);
        let (hist, confl, carved, maxs) = structure_metrics(&out, &ss);
        let (d5, d20, d50) = deposition_locality(&fbm, &out);
        eprintln!(
            "  {f:.2} {dens:>5.2} | {:>4.0} | {:>+4.0}% ({:>3.0}%) | {maxs:>4} | {confl:>5} | {d5:>3.0}/{d20:>3.0}/{d50:>3.0}% | {:.1}",
            (e - d) / e.max(1.0) * 100.0,
            (carved - carved_fbm) * 100.0,
            carved * 100.0,
            emerged_frac(&out) * 100.0,
        );
        let _ = hist;
    }
    eprintln!(
        "  (delta/beach survives while dep≤5c > 0; carvedΔ>0 with a deep histogram = real incision.)"
    );
}

/// PART A — settle the terrace attribution on the AUTHOR'S banded map
/// (seed 10481999410520546993, domain 400 km — domain only relabels km, not shape).
/// The decisive, LL-free test: does the exported height.u16 vary CONTINUOUSLY, or
/// show long runs of identical codes with sharp jumps (data-side plateaux)? Perfect
/// level sets in the DATA ⇒ real terraces; smooth codes ⇒ the banding is a Living
/// Landz DISPLAY issue and deposition/transport rework leaves this chantier's scope.
#[test]
#[ignore]
fn part_a_bigmap_terrace() {
    use ymir_core::export::height::metric_height_u16;
    let ss = SteinSteinParams::default();
    let seed_u = 10481999410520546993u64;
    let t = 2048usize; // proxy for 8192²: u16-code continuity is ~scale-invariant
    let (state, _run) = coarse_state(seed_u);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(seed_u);
    let mut fcfg = FbmUpscaleConfig::c1_hd_production(t);
    fcfg.erosion = None;
    fcfg.bathymetry = None;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fcfg).heightmap;
    let eroded = run_erosion(&fbm, &ero_cfg(t), &seed, |_, _, _| true).heightmap; // f=0.25 default

    eprintln!("\n=== PART A — seed {seed_u} @ {t}² (author's banded map) ===");
    // S̃ closure levels (same test as seed 42 — does it behave the same?).
    let s = &state.s;
    let vals: Vec<f64> = (0..s.ny() * s.nx()).map(|k| s.get(k % s.nx(), k / s.nx())).collect();
    let n = vals.len() as f32;
    let at_heq = vals.iter().filter(|&&v| (v - 2.0).abs() < 1e-4).count();
    eprintln!(
        "  S̃: {:.0}% at h_eq=2.0 (equilibrium clamp), {:.0}% above h_eq — {}",
        at_heq as f32 / n * 100.0,
        vals.iter().filter(|&&v| v > 2.0 + 1e-4).count() as f32 / n * 100.0,
        if at_heq == 0 { "clamp inactive (same as seed 42)" } else { "clamp ACTIVE" },
    );

    // Flat fraction (disentangle) on the eroded field transect.
    let r = max_row(&eroded);
    terrace_report("eroded (bigmap)", &row_of(&to_metres(&eroded, &ss), r));

    // THE decisive test — raw exported u16 codes.
    let mh = metric_height_u16(&eroded, &ss);
    let (w, h) = (eroded.width, eroded.height);
    eprintln!(
        "  height.u16: range {:.0}..{:.0} m over 65535 codes ⇒ {:.3} m/code",
        mh.min_m,
        mh.max_m,
        (mh.max_m - mh.min_m) / 65535.0
    );
    for (label, codes) in [
        ("row", (0..w).map(|x| mh.codes[r * w + x]).collect::<Vec<u16>>()),
        ("col", (0..h).map(|y| mh.codes[y * w + w / 2]).collect::<Vec<u16>>()),
    ] {
        // |Δcode| distribution + longest run of identical codes.
        let (mut z, mut small, mut big) = (0usize, 0usize, 0usize); // Δ=0, 1..=8, >8
        let (mut run, mut maxrun) = (1usize, 1usize);
        for i in 1..codes.len() {
            let d = (codes[i] as i32 - codes[i - 1] as i32).unsigned_abs();
            match d {
                0 => z += 1,
                1..=8 => small += 1,
                _ => big += 1,
            }
            if codes[i] == codes[i - 1] {
                run += 1;
                maxrun = maxrun.max(run);
            } else {
                run = 1;
            }
        }
        let tot = (codes.len() - 1) as f32;
        eprintln!(
            "  [{label}] Δcode: =0 {:.0}%, 1–8 {:.0}%, >8 {:.0}%; longest identical-code run {maxrun} cells",
            z as f32 / tot * 100.0,
            small as f32 / tot * 100.0,
            big as f32 / tot * 100.0,
        );
    }
    eprintln!(
        "  VERDICT: mostly Δ=0 with occasional >8 jumps + long identical runs ⇒ DATA plateaux (real \
         terraces). Mostly 1–8 with short runs ⇒ CONTINUOUS data ⇒ banding is a Living Landz DISPLAY \
         issue (deposition/transport rework leaves scope)."
    );
}

/// PART B — locate the terrace source in the closure chain. Measures the S̃ field
/// (crustal thickness) from the coarse pass for hard-clamp spikes at the global
/// equilibrium level h_eq (2.0) and the Davis-Suppe cap h_max (2.5), and ties the
/// altitude plateaux to the S̃≈h_eq cells. Read-only.
#[test]
#[ignore]
fn terrace_source_closure() {
    let ss = SteinSteinParams::default();
    let (state, _run) = coarse_state(SEED);
    let s = &state.s; // S̃ crustal thickness (the closures' state)
    let (nx, ny) = (s.nx(), s.ny());
    let vals: Vec<f64> =
        (0..ny).flat_map(|j| (0..nx).map(move |i| (i, j))).map(|(i, j)| s.get(i, j)).collect();
    let n = vals.len() as f32;
    let (h_eq, h_max) = (2.0f64, 2.5f64);
    let near = |v: f64, t: f64| (v - t).abs() < 1e-4;
    let at_heq = vals.iter().filter(|&&v| near(v, h_eq)).count();
    let at_hmax = vals.iter().filter(|&&v| near(v, h_max)).count();
    let above_heq = vals.iter().filter(|&&v| v > h_eq + 1e-4).count();
    eprintln!("\n=== PART B — terrace source in the closure chain (seed {SEED}, coarse S̃) ===");
    eprintln!(
        "  S̃ cells: {} total. HARD-CLAMP spike at h_eq=2.0 (equilibrium_height): {} ({:.0}%). \
         at h_max=2.5 (davis_suppe): {} ({:.0}%). above h_eq (still relaxing): {} ({:.0}%).",
        vals.len(),
        at_heq,
        at_heq as f32 / n * 100.0,
        at_hmax,
        at_hmax as f32 / n * 100.0,
        above_heq,
        above_heq as f32 / n * 100.0,
    );

    // Distinct S̃ levels: fine histogram, report bins holding > 2% of cells.
    let (mn, mx) = vals.iter().fold((f64::MAX, f64::MIN), |(a, b), &v| (a.min(v), b.max(v)));
    const NB: usize = 200;
    let mut bins = vec![0u32; NB];
    let span = (mx - mn).max(1e-9);
    for &v in &vals {
        bins[(((v - mn) / span) * (NB as f64 - 1.0)) as usize] += 1;
    }
    eprint!("  S̃ range {mn:.2}..{mx:.2}; dominant levels (>2% of cells): ");
    for (b, &c) in bins.iter().enumerate() {
        if c as f32 / n > 0.02 {
            eprint!(
                "{:.2}({:.0}%) ",
                mn + (b as f64 + 0.5) / NB as f64 * span,
                c as f32 / n * 100.0
            );
        }
    }
    eprintln!();

    // Tie S̃≈h_eq cells to a single altitude plateau (in metres).
    let alt = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let alt_m = to_metres(&alt, &ss);
    let mut plat_m: Vec<f32> = Vec::new();
    for k in 0..vals.len() {
        if near(vals[k], h_eq) {
            plat_m.push(alt_m.data[k]);
        }
    }
    plat_m.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if !plat_m.is_empty() {
        let (lo, hi, md) = (plat_m[0], plat_m[plat_m.len() - 1], plat_m[plat_m.len() / 2]);
        eprintln!(
            "  altitude of the S̃=h_eq plateau: {lo:.0}..{hi:.0} m (median {md:.0} m) — a tight band \
             ⇒ the global h_eq maps to a single altitude LEVEL SET (concentric terrace).",
        );
    }

    // LAND altitude levels (metres) — are there discrete steps, and what makes them?
    let mut land_m: Vec<f32> = alt_m.data.iter().copied().filter(|&m| m > 0.0).collect();
    land_m.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let (lmn, lmx) = (land_m[0], land_m[land_m.len() - 1]);
    const AB: usize = 80;
    let mut abins = vec![0u32; AB];
    let aspan = (lmx - lmn).max(1.0);
    for &m in &land_m {
        abins[(((m - lmn) / aspan) * (AB as f32 - 1.0)) as usize] += 1;
    }
    let amean = land_m.len() as f32 / abins.iter().filter(|&&c| c > 0).count() as f32;
    eprint!(
        "  LAND altitude {lmn:.0}..{lmx:.0} m; spikes (>3× mean bin, width {:.0} m): ",
        aspan / AB as f32
    );
    for (b, &c) in abins.iter().enumerate() {
        if c as f32 > 3.0 * amean {
            eprint!(
                "{:.0}m({:.0}%) ",
                lmn + (b as f32 + 0.5) / AB as f32 * aspan,
                c as f32 / land_m.len() as f32 * 100.0
            );
        }
    }
    eprintln!();

    // Craton correlation: cratons get compute_isostasy_craton (worn-shield altitude)
    // — a candidate discrete level whose boundary is the craton-mask contour.
    let craton = state.cratonic_mask.data();
    let mut cra: Vec<f32> = Vec::new();
    let mut noncra: Vec<f32> = Vec::new();
    for k in 0..alt_m.data.len() {
        if alt_m.data[k] > 0.0 {
            if craton[k] {
                cra.push(alt_m.data[k]);
            } else {
                noncra.push(alt_m.data[k]);
            }
        }
    }
    let med = |v: &mut Vec<f32>| {
        if v.is_empty() {
            return f32::NAN;
        }
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v[v.len() / 2]
    };
    let (cn, ncn) = (cra.len(), noncra.len());
    let _ = med(&mut cra);
    let _ = med(&mut noncra);
    let pct = |v: &[f32], p: f32| v[((v.len() as f32 - 1.0) * p) as usize];
    eprintln!(
        "  cratonic land {cn} cells: altitude p10/median/p90 = {:.0}/{:.0}/{:.0} m; \
         non-cratonic {ncn}: {:.0}/{:.0}/{:.0} m. Cratonic band is WIDE (not a level set) → \
         cratons are NOT the terrace source either.",
        pct(&cra, 0.1),
        pct(&cra, 0.5),
        pct(&cra, 0.9),
        pct(&noncra, 0.1),
        pct(&noncra, 0.5),
        pct(&noncra, 0.9),
    );
    eprintln!(
        "  VERDICT: no coarse discrete level found (h_eq 0%, h_max 0%, craton wide). The flat fraction \
         is dominated by EROSION DEPOSITION (disentangle: bilinear 13% → FBM 6% → eroded 24%), i.e. the \
         same erosion-algorithm limit as the missing valleys — NOT a tectonic/isostasy closure."
    );

    let km_cell_1024 = C1_DOMAIN_KM / 64.0;
    eprintln!(
        "  vertical scale: depth_scale_m {} (uncoupled to domain today); ~120–176 m coarse step > a \
         hamlet valley's relief (150–250 m WIDE at {:.0} km/coarse-cell), not cosmetic. If depth_scale_m \
         were coupled to domain_km the step would scale ∝ domain_km.",
        ss.depth_scale_m, km_cell_1024,
    );
}

/// Median channel incision (m) at top-1% accumulation LAND cells: how much
/// stream-power lowered the channels vs the pre-incision field. For K calibration.
fn median_channel_incision_m(before: &GridF32, after: &GridF32, ss: &SteinSteinParams) -> f32 {
    let dr = c1_drainage(before, None, &C1DrainageConfig::default(), ss);
    let acc = &dr.flow.accumulation;
    let mut accs = acc.data.clone();
    accs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let thr = accs[(accs.len() as f64 * 0.99) as usize];
    let mut inc: Vec<f32> = (0..before.data.len())
        .filter(|&k| acc.data[k] >= thr && before.data[k] > SEA)
        .map(|k| {
            c1_altitude_norm_to_metres(before.data[k], ss)
                - c1_altitude_norm_to_metres(after.data[k], ss)
        })
        .collect();
    if inc.is_empty() {
        return 0.0;
    }
    inc.sort_by(|a, b| a.partial_cmp(b).unwrap());
    inc[inc.len() / 2]
}

/// V/U-ness of the highest-order channel: cross-section incision (m) below the rim.
fn channel_incision_profile(field: &GridF32, ss: &SteinSteinParams) -> f32 {
    let dr = c1_drainage(field, None, &C1DrainageConfig::default(), ss);
    let (w, h) = (field.width, field.height);
    let fm = to_metres(field, ss);
    if let Some(seg) =
        dr.rivers.segments.iter().filter(|s| s.points.len() >= 3).max_by_key(|s| s.strahler_order)
    {
        let mid = seg.points.len() / 2;
        let (ax, ay) = seg.points[mid - 1];
        let (bx, by) = seg.points[mid + 1];
        let (tx, ty) = (bx as f32 - ax as f32, by as f32 - ay as f32);
        let tl = (tx * tx + ty * ty).sqrt().max(1e-6);
        let (px, py) = (-ty / tl, tx / tl);
        let (cx, cy) = seg.points[mid];
        let mut xs = Vec::new();
        for o in -8i32..=8 {
            let sx = (cx as f32 + px * o as f32).round().clamp(0.0, w as f32 - 1.0) as i32;
            let sy = (cy as f32 + py * o as f32).round().clamp(0.0, h as f32 - 1.0) as i32;
            xs.push(fm.get(sx, sy));
        }
        let bottom = xs.iter().cloned().fold(f32::MAX, f32::min);
        return xs[0].max(xs[16]) - bottom;
    }
    0.0
}

/// PART B — striation baseline + FBM knob sweep. Minimise FBM amplitude/anisotropy
/// while keeping drainage ORGANIC (Strahler depth, confluences, no grid/radial
/// alignment). Reports the striation metric + drainage health for each knob value.
/// 1024², seed 42, relief-v1 incision on top of each FBM variant.
#[test]
#[ignore]
fn fbm_striation() {
    use ymir_core::erosion::stream_power::{RELIEF_V1_A_C_KM2, StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let (t, domain) = (1024usize, 400.0f32);
    let cell_km2 = (domain / t as f32).powi(2);
    let base_depth = ss.depth_scale_m as f32;
    let a_c = RELIEF_V1_A_C_KM2 / cell_km2;
    let (state, _run) = coarse_state(SEED);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(SEED);
    let sp = StreamPowerConfig::relief_v1(cell_km2, ss.depth_scale_m as f32);

    let mut report = |label: String, mut fc: FbmUpscaleConfig| {
        fc.erosion = None;
        fc.bathymetry = None;
        let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fc).heightmap;
        // Pre-incision striation (raw FBM, undiluted by stream-power valley walls).
        let slope_fbm = slope_deg_field(&fbm, domain, base_depth);
        let (_, _, asym_pre, wl_pre) = striation_metric(&fbm, &slope_fbm, 20.0);
        let field = incise(&fbm, &sp);
        let slope = slope_deg_field(&field, domain, base_depth);
        let (rg, rc, asym, wl) = striation_metric(&field, &slope, 20.0);
        let _ = (rg, rc, wl);
        eprint!("  [pre-FBM asym {asym_pre:.2} λ {wl_pre:.1}]");
        let dr = c1_drainage(&field, None, &C1DrainageConfig::default(), &ss);
        let maxs = dr.rivers.segments.iter().map(|s| s.strahler_order).max().unwrap_or(0);
        let confl = dr.rivers.segments.iter().filter(|s| s.upstream.len() >= 2).count();
        let nseg = dr.rivers.segments.len();
        let corr = channel_corridor(&field, &ss, a_c);
        let (vf5, _, _, _) = valley_floor(&slope, &corr, &field, cell_km2);
        eprintln!(
            "  {label:<26}: striation asym {asym:.2} (grad {rg:.4} contour {rc:.4}) λ {wl:.1} cells | \
             maxS {maxs} confl {confl} segs {nseg} vfloor<5° {vf5:.0} km²",
        );
    };

    let b = || {
        let mut c = FbmUpscaleConfig::c1_hd_production(t);
        c.erosion = None;
        c.bathymetry = None;
        c
    };
    eprintln!("\n=== PART B — FBM striation baseline + knob sweep (relief-v1 incision, 1024²) ===");
    report("BASELINE (prod FBM)".into(), b());
    eprintln!("-- max_anisotropy (stretch along slope → filaments) --");
    for a in [3.0f64, 2.0, 1.0] {
        report(format!("max_anisotropy={a}"), FbmUpscaleConfig { max_anisotropy: a, ..b() });
    }
    eprintln!("-- amplitude_slope_factor (slope-selective amplitude) --");
    for a in [3.0f64, 1.0, 0.0] {
        report(
            format!("amp_slope_factor={a}"),
            FbmUpscaleConfig { amplitude_slope_factor: a, ..b() },
        );
    }
    eprintln!("-- amplitude_base (overall FBM amplitude) --");
    for a in [0.16f64, 0.08, 0.04, 0.02] {
        report(format!("amplitude_base={a}"), FbmUpscaleConfig { amplitude_base: a, ..b() });
    }
    eprintln!("-- octaves --");
    for o in [7usize, 5, 3] {
        report(format!("octaves={o}"), FbmUpscaleConfig { octaves: o, ..b() });
    }
    eprintln!(
        "  (want: striation asym → 1, λ large, WHILE maxS/confl/segs stay healthy = drainage organic."
    );
    eprintln!("   the floor is the lowest amplitude/anisotropy before the network degenerates.)");
}

/// PART C — valley width/depth as a DISTRIBUTION per Strahler order. Not a target to
/// maximise: the test is whether W/D WIDENS DOWNSTREAM (higher order = wider) and
/// whether wide valleys EXIST (a tail), not a single median. 1024², relief-v1.
#[test]
#[ignore]
fn valley_width_distribution() {
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let (t, domain) = (1024usize, 400.0f32);
    let cell_km2 = (domain / t as f32).powi(2);
    let cell_m = domain / t as f32 * 1000.0;
    let base = ss.depth_scale_m as f32;
    let (state, _run) = coarse_state(SEED);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(SEED);
    let mut fcfg = FbmUpscaleConfig::c1_hd_production(t);
    fcfg.erosion = None;
    fcfg.bathymetry = None;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fcfg).heightmap;
    let sp = incise(&fbm, &StreamPowerConfig::relief_v1(cell_km2, ss.depth_scale_m as f32));
    let slope = slope_deg_field(&sp, domain, base);
    let (w, h) = (t, t);
    let fm: Vec<f32> = sp.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).collect();
    let dr = c1_drainage(&sp, None, &C1DrainageConfig::default(), &ss);
    // Per segment: cross-section W/D at the midpoint, grouped by Strahler order.
    let mut per: std::collections::HashMap<u8, Vec<f32>> = std::collections::HashMap::new();
    for s in dr.rivers.segments.iter().filter(|s| s.points.len() >= 5) {
        let mid = s.points.len() / 2;
        let (ax, ay) = s.points[mid - 1];
        let (bx, by) = s.points[mid + 1];
        let (tx, ty) = (bx as f32 - ax as f32, by as f32 - ay as f32);
        let tl = (tx * tx + ty * ty).sqrt().max(1e-6);
        let (px, py) = (-ty / tl, tx / tl);
        let (cx, cy) = (s.points[mid].0 as usize, s.points[mid].1 as usize);
        let r = 40i32;
        let samp = |o: i32| -> (f32, f32) {
            let sx = (cx as f32 + px * o as f32).round().clamp(0.0, w as f32 - 1.0) as usize;
            let sy = (cy as f32 + py * o as f32).round().clamp(0.0, h as f32 - 1.0) as usize;
            (fm[sy * w + sx], slope[sy * w + sx])
        };
        let (mut bo, mut bmin) = (0i32, f32::MAX);
        for o in -r..=r {
            if samp(o).0 < bmin {
                bmin = samp(o).0;
                bo = o;
            }
        }
        let side = |dir: i32| -> Option<(i32, f32)> {
            let mut o = bo;
            loop {
                o += dir;
                if o.abs() > r {
                    return None;
                }
                let (a, sl) = samp(o);
                if sl > 30.0 {
                    return Some((o, a));
                }
            }
        };
        if let (Some((lo, la)), Some((ro, ra))) = (side(-1), side(1)) {
            let width = (ro - lo) as f32 * cell_m;
            let depth = la.min(ra) - bmin;
            if depth > 1.0 {
                per.entry(s.strahler_order).or_default().push(width / depth);
            }
        }
    }
    eprintln!("\n=== PART C — W/D distribution per Strahler order (relief-v1, 1024²) ===");
    let mut orders: Vec<u8> = per.keys().copied().collect();
    orders.sort();
    let mut medians = Vec::new();
    for o in &orders {
        let v = per.get_mut(o).unwrap();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let (p10, med, p90) = (v[v.len() / 10], v[v.len() / 2], v[v.len() * 9 / 10]);
        medians.push(med);
        eprintln!("  S{o}: n={:>4} W/D p10/median/p90 = {p10:.1}/{med:.1}/{p90:.1}", v.len());
    }
    let widens =
        medians.windows(2).filter(|w| w[1] >= w[0]).count() >= medians.len().saturating_sub(1) / 2;
    eprintln!(
        "  → downstream widening: {} (median W/D {} with order); wide valleys exist if p90 has a high tail",
        if widens { "YES ✓" } else { "NO — gorges everywhere" },
        if widens { "RISES" } else { "does NOT rise" },
    );
}

/// PART D — resolution dependence of incision. E = K·A^m·S^n with A in CELLS and S =
/// Δnorm per cell: a finer cell samples the FBM detail at a steeper local gradient,
/// so S (and incision) rises with resolution. Measure per-order incision at 512 /
/// 1024 / 2048 to size the dependency; the clean fix is a PHYSICAL slope
/// (Δh_m / cell_m) + A in km², not a per-resolution K.
#[test]
#[ignore]
fn incision_resolution() {
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let domain = 400.0f32;
    let (state, _run) = coarse_state(SEED);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(SEED);
    eprintln!(
        "\n=== PART D — incision vs resolution (relief-v1, A_c in km², domain {domain} km) ==="
    );
    for t in [512usize, 1024, 2048] {
        let cell_km2 = (domain / t as f32).powi(2);
        let mut fcfg = FbmUpscaleConfig::c1_hd_production(t);
        fcfg.erosion = None;
        fcfg.bathymetry = None;
        let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fcfg).heightmap;
        let sp = incise(&fbm, &StreamPowerConfig::relief_v1(cell_km2, ss.depth_scale_m as f32));
        let (tab, _) = per_order_incision(&fbm, &sp, &ss);
        eprintln!("  {t}²: {}", fmt_orders(&tab));
    }
    eprintln!(
        "  (rising with resolution ⇒ S measured steeper on finer cells; fix = physical slope Δh_m/cell_m + A km²)"
    );
}

/// Hillshade (NW light, 45° altitude) of an altitude field → grayscale [0,1], the
/// clearest view of fine striations on steep flanks. Uses the effective metre scale.
fn hillshade(field: &GridF32, domain_km: f32, depth_scale: f32) -> GridF32 {
    let (w, h) = (field.width, field.height);
    let cell_m = domain_km / w as f32 * 1000.0;
    let norm_to_m = 2.0 * 1.13 * depth_scale;
    let (lx, ly, lz) = (-0.5f32, 0.5, 0.707); // NW, 45° up
    let mut out = GridF32::new(w, h, 0.0);
    for y in 1..h - 1 {
        for x in 1..w - 1 {
            let dzdx =
                (field.data[y * w + x + 1] - field.data[y * w + x - 1]) * 0.5 * norm_to_m / cell_m;
            let dzdy =
                (field.data[(y + 1) * w + x] - field.data[(y - 1) * w + x]) * 0.5 * norm_to_m
                    / cell_m;
            let inv = 1.0 / (dzdx * dzdx + dzdy * dzdy + 1.0).sqrt();
            let shade = (-dzdx * lx - dzdy * ly + lz) * inv;
            out.data[y * w + x] = shade.clamp(0.0, 1.0);
        }
    }
    out
}

/// TASK 1 — render the FBM amplitude ladder (relief-v1 incision ON) as hillshade
/// PNGs + drainage figures, so the author can SEE whether the striations disappear
/// and the terrain still reads as terrain. Author seed, domain 400 km. Ladder at
/// 2048²; one at 8192² for the recommended value.
#[test]
#[ignore]
fn render_striation_ladder() {
    use std::path::Path;
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let seed_u = 10481999410520546993u64;
    let domain = 400.0f32;
    let base = ss.depth_scale_m as f32;
    let (state, _run) = coarse_state(seed_u);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(seed_u);
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../exports/relief_ladder");
    std::fs::create_dir_all(&dir).unwrap();
    eprintln!("\n=== TASK 1 — striation amplitude ladder (hillshade PNGs) → {} ===", dir.display());

    let render = |t: usize, amp: f64| {
        let cell_km2 = (domain / t as f32).powi(2);
        let mut fc = FbmUpscaleConfig::c1_hd_production(t);
        fc.erosion = None;
        fc.bathymetry = None;
        fc.amplitude_base = amp;
        let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fc).heightmap;
        let sp = incise(&fbm, &StreamPowerConfig::relief_v1(cell_km2, ss.depth_scale_m as f32));
        let hs = hillshade(&sp, domain, base);
        let path = dir.join(format!("amp{amp:.2}_{t}.png"));
        hs.save_png_u8(&path).unwrap();
        // Also a crop (central 1/4) so fine striations are visible without zoom.
        let (w, hh) = (sp.width, sp.height);
        let (x0, y0, cw, ch) = (w / 2 - w / 8, hh / 2 - hh / 8, w / 4, hh / 4);
        let mut crop = GridF32::new(cw, ch, 0.0);
        for j in 0..ch {
            for i in 0..cw {
                crop.data[j * cw + i] = hs.data[(y0 + j) * w + (x0 + i)];
            }
        }
        crop.save_png_u8(&dir.join(format!("amp{amp:.2}_{t}_crop.png"))).unwrap();
        let dr = c1_drainage(&sp, None, &C1DrainageConfig::default(), &ss);
        let maxs = dr.rivers.segments.iter().map(|s| s.strahler_order).max().unwrap_or(0);
        let confl = dr.rivers.segments.iter().filter(|s| s.upstream.len() >= 2).count();
        let slope = slope_deg_field(&sp, domain, base);
        let corr = channel_corridor(
            &sp,
            &ss,
            StreamPowerConfig::relief_v1(cell_km2, ss.depth_scale_m as f32).min_area_cells,
        );
        let (vf5, _, _, _) = valley_floor(&slope, &corr, &sp, cell_km2);
        eprintln!(
            "  amp {amp:.2} @{t}²: maxS {maxs}, confl {confl}, valley floor<5° {vf5:.0} km²  → {}",
            path.display()
        );
    };

    for amp in [0.16f64, 0.08, 0.04, 0.02] {
        render(2048, amp);
    }
    // Recommended value at production scale (provisional — the author's eye decides).
    render(8192, 0.04);
    eprintln!(
        "  RECOMMEND (provisional): amp 0.04 — 4× reduction, drainage healthy; author confirms the visual."
    );
}

/// Recalibrate the PHYSICAL K to reproduce the relief-v1 reference (1024²: drainage
/// relief ~682 m, per-order S2~414 S4~136). The old K=3 was against the normalised
/// slope; the physical law needs a different K.
#[test]
#[ignore]
fn calibrate_k_physical() {
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let (t, domain) = (1024usize, 400.0f32);
    let cell_km2 = (domain / t as f32).powi(2);
    let (state, _run) = coarse_state(SEED);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(SEED);
    let mut fcfg = FbmUpscaleConfig::c1_hd_production(t);
    fcfg.erosion = None;
    fcfg.bathymetry = None;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fcfg).heightmap;
    eprintln!("\n=== physical-K calibration (target: relief ~682 m, S2~414 S4~136) ===");
    for k in [100.0f32, 300.0, 1000.0, 3000.0, 10000.0] {
        let mut cfg = StreamPowerConfig::relief_v1(cell_km2, ss.depth_scale_m as f32);
        cfg.k = k;
        let sp = incise(&fbm, &cfg);
        let rel = drainage_relief_m(&sp, &ss);
        let (tab, _) = per_order_incision(&fbm, &sp, &ss);
        eprintln!("  K={k:>7.0}: relief {rel:>4.0} m, {}", fmt_orders(&tab));
    }
}

/// CLOSURE mosaic — stitch the `closure_*_crop` renders into one comparison image
/// (`exports/sculpt/closure_mosaic.png`) so v1 vs v2 and 2048² vs 8192² sit side by
/// side. Each crop is resampled (box filter) to a common 512² panel; white separators
/// between panels. Layout (2×3):
///   [ v1 2048 | v2 2048 | v2 2048 amp0.01 ]
///   [ (blank) | v2 8192 | v2 8192 amp0.02 ]
/// Run `closure_render` FIRST so the crops exist. Read-only (loads PNGs).
/// Axial (0..π) orientation concentration of a set of angles: `R = |mean(e^{i2θ})|`
/// (1 = all aligned, 0 = uniform) + a 12-bin histogram (fractions). Used to test the
/// striation source: high `R` in the GRID frame ⇒ grid-aligned (D8 routing); high `R`
/// in the LOCAL-GRADIENT frame ⇒ relief-aligned (anisotropic FBM). A grid-frame test
/// ALONE cannot tell them apart, hence both frames.
fn axial_concentration(angles: &[f32]) -> (f32, [f32; 12]) {
    if angles.is_empty() {
        return (0.0, [0.0; 12]);
    }
    let (mut cx, mut cy) = (0.0f64, 0.0f64);
    let mut bins = [0.0f32; 12];
    for &a in angles {
        let t = (a.rem_euclid(std::f32::consts::PI)) as f64;
        cx += (2.0 * t).cos();
        cy += (2.0 * t).sin();
        let b = ((t / std::f64::consts::PI * 12.0) as usize).min(11);
        bins[b] += 1.0;
    }
    let n = angles.len() as f64;
    let r = ((cx / n).powi(2) + (cy / n).powi(2)).sqrt() as f32;
    for v in bins.iter_mut() {
        *v /= angles.len() as f32;
    }
    (r, bins)
}

/// CLOSURE STEP 1 (read-only) — DISCRIMINATE the 8192² striation source: D8 routing vs
/// anisotropic FBM. (1a) two-frame orientation histograms of channel segments at 2048²
/// and 8192² — grid frame (D8 ⇒ peaks at 0/45/90/135°) vs local-gradient frame (FBM ⇒
/// peaks relative to slope). (1b) the FREE ABLATION: relief-v2 at 8192² with isotropic,
/// slope-blind noise (`max_anisotropy=1`, `amplitude_slope_factor=0`) vs the anisotropic
/// default — if the comb dies, it is the FBM and no routing rewrite is needed. Saves
/// iso/aniso 8192² crops for the eye. No flow.rs change.
#[test]
#[ignore]
fn striation_source() {
    use std::path::Path;
    use std::time::Instant;
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let seed_u = 10481999410520546993u64;
    let domain = 400.0f32;
    let base = ss.depth_scale_m as f32;
    let (state, _run) = coarse_state(seed_u);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(seed_u);
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../exports/sculpt");
    std::fs::create_dir_all(&dir).unwrap();

    // Build a relief-v2 field at (res, amp), optionally with ISOTROPIC slope-blind FBM.
    let build = |res: usize, amp: f64, iso: bool| -> GridF32 {
        let mut fc = FbmUpscaleConfig::c1_hd_production(res);
        fc.erosion = None;
        fc.bathymetry = None;
        fc.amplitude_base = amp;
        if iso {
            fc.max_anisotropy = 1.0; // no directional stretching
            fc.amplitude_slope_factor = 0.0; // amplitude independent of slope
        }
        let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fc).heightmap;
        let cell_km2 = (domain / res as f32).powi(2);
        incise(&fbm, &StreamPowerConfig::relief_v2(cell_km2, base))
    };

    // ── (1a) Two-frame channel-segment orientations at 2048² and 8192². ──
    let seg_orient = |f: &GridF32| -> (Vec<f32>, Vec<f32>) {
        let (w, h) = (f.width, f.height);
        let dr = c1_drainage(f, None, &C1DrainageConfig::default(), &ss);
        let (mut grid, mut rel) = (Vec::new(), Vec::new());
        for s in &dr.rivers.segments {
            for win in s.points.windows(2) {
                let (x0, y0) = (win[0].0 as i32, win[0].1 as i32);
                let (dx, dy) = ((win[1].0 as i32 - x0) as f32, (win[1].1 as i32 - y0) as f32);
                if dx == 0.0 && dy == 0.0 {
                    continue;
                }
                let seg = dy.atan2(dx);
                grid.push(seg);
                // local gradient (central diff, clamped to interior).
                let (xi, yi) =
                    (x0.clamp(1, w as i32 - 2) as usize, y0.clamp(1, h as i32 - 2) as usize);
                let k = yi * w + xi;
                let gx = f.data[k + 1] - f.data[k - 1];
                let gy = f.data[k + w] - f.data[k - w];
                if gx == 0.0 && gy == 0.0 {
                    continue;
                }
                rel.push(seg - gy.atan2(gx));
            }
        }
        (grid, rel)
    };

    eprintln!("\n=== CLOSURE STEP 1 — striation source (seed {seed_u}) ===");
    eprintln!("-- (1a) channel-segment orientation concentration R (1=aligned, 0=uniform) --");
    eprintln!("   res  | R_grid | R_gradrel | grid histogram (0°..180°, 15° bins)");
    for res in [2048usize, 8192] {
        let f = build(res, 0.04, false);
        let (grid, rel) = seg_orient(&f);
        let (rg, bins) = axial_concentration(&grid);
        let (rr, _) = axial_concentration(&rel);
        let hb = bins.iter().map(|b| format!("{:.2}", b)).collect::<Vec<_>>().join(" ");
        eprintln!("  {res:>5} | {rg:>6.3} | {rr:>9.3} | {hb}");
    }
    eprintln!(
        "   (D8 ⇒ high R_grid, peaks at bins 0/3/6/9 = 0/45/90/135°; FBM ⇒ R_gradrel > R_grid)"
    );

    // ── (1b) The free ablation at 8192²: isotropic vs anisotropic FBM. ──
    eprintln!("\n-- (1b) ablation @8192²: isotropic slope-blind FBM vs anisotropic default --");
    eprintln!("   variant     | >30° | max° | upper-slope striation contour/gradient power");
    let (fx, fy, fw) = (0.22f32, 0.03f32, 0.34f32);
    for (label, iso) in [("anisotropic", false), ("isotropic", true)] {
        let t0 = Instant::now();
        let f = build(8192, 0.04, iso);
        let ms = t0.elapsed().as_millis();
        let slope = slope_deg_field(&f, domain, base);
        let sm: Vec<f32> = f.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).collect();
        let land: Vec<usize> = (0..f.data.len()).filter(|&k| f.data[k] > SEA).collect();
        let n = land.len().max(1) as f32;
        let a30 = land.iter().filter(|&&k| slope[k] > 30.0).count() as f32 / n * 100.0;
        let mx = land.iter().map(|&k| slope[k]).fold(0.0f32, f32::max);
        let mut le: Vec<f32> = land.iter().map(|&k| sm[k]).collect();
        le.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let e60 = le[le.len() * 6 / 10];
        let upper: Vec<f32> =
            (0..f.data.len()).map(|k| if sm[k] >= e60 { slope[k] } else { 0.0 }).collect();
        let (_, aniso_pow) = striation_spectrum(&f, &upper, 20.0, 48);
        eprintln!("   {label:<11} | {a30:>4.1}% | {mx:.0} | {aniso_pow:.2}  ({ms} ms)");
        // Save the massif crop for the eye.
        let hs = hillshade(&f, domain, base);
        let (w, hh) = (f.width, f.height);
        let (cx0, cy0, cw) =
            ((fx * w as f32) as usize, (fy * hh as f32) as usize, (fw * w as f32) as usize);
        let mut crop = GridF32::new(cw, cw, 0.5);
        for j in 0..cw {
            for i in 0..cw {
                if cx0 + i < w && cy0 + j < hh {
                    crop.data[j * cw + i] = hs.data[(cy0 + j) * w + (cx0 + i)];
                }
            }
        }
        crop.save_png_u8(&dir.join(format!("striation_8192_{label}.png"))).unwrap();
    }
    eprintln!("   → crops: exports/sculpt/striation_8192_anisotropic.png / _isotropic.png");
    eprintln!("   (comb dies with isotropic ⇒ FBM is the source, no routing rewrite needed)");
}

/// CLOSURE STEP 1 (cont.) — isolate WHERE the 8192² contour-terraces are born: the FBM
/// upscale, or the erosion. Renders hillshade crops of the RAW FBM (pre-incision) and
/// the post relief-v2 field at 8192², isotropic slope-blind noise, + the upper-slope
/// striation power of each. If the terraces are already in the raw FBM → the upscale;
/// if they appear only after incision → the erosion coupling. No pipeline change.
#[test]
#[ignore]
fn striation_stage() {
    use std::path::Path;
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let seed_u = 10481999410520546993u64;
    let (domain, base, res) = (400.0f32, ss.depth_scale_m as f32, 8192usize);
    let (state, _run) = coarse_state(seed_u);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(seed_u);
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../exports/sculpt");
    let (fx, fy, fw) = (0.22f32, 0.03f32, 0.34f32);
    let mut fc = FbmUpscaleConfig::c1_hd_production(res);
    fc.erosion = None;
    fc.bathymetry = None;
    fc.amplitude_base = 0.04;
    fc.max_anisotropy = 1.0;
    fc.amplitude_slope_factor = 0.0;
    let raw = upscale_with_fbm(&coarse, SEA, &seed, &fc).heightmap;
    let cell_km2 = (domain / res as f32).powi(2);
    let post = incise(&raw, &StreamPowerConfig::relief_v2(cell_km2, base));

    eprintln!("\n=== CLOSURE STEP 1 — terrace stage @8192² isotropic (upscale vs erosion) ===");
    eprintln!("   stage            | >30° | upper-slope striation contour/gradient power");
    for (label, f) in [("raw FBM (pre)", &raw), ("relief-v2 (post)", &post)] {
        let slope = slope_deg_field(f, domain, base);
        let sm: Vec<f32> = f.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).collect();
        let land: Vec<usize> = (0..f.data.len()).filter(|&k| f.data[k] > SEA).collect();
        let n = land.len().max(1) as f32;
        let a30 = land.iter().filter(|&&k| slope[k] > 30.0).count() as f32 / n * 100.0;
        let mut le: Vec<f32> = land.iter().map(|&k| sm[k]).collect();
        le.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let e60 = le[le.len() * 6 / 10];
        let upper: Vec<f32> =
            (0..f.data.len()).map(|k| if sm[k] >= e60 { slope[k] } else { 0.0 }).collect();
        let (_, pow) = striation_spectrum(f, &upper, 20.0, 48);
        eprintln!("   {label:<16} | {a30:>4.1}% | {pow:.2}");
        let hs = hillshade(f, domain, base);
        let (w, hh) = (f.width, f.height);
        let (cx0, cy0, cw) =
            ((fx * w as f32) as usize, (fy * hh as f32) as usize, (fw * w as f32) as usize);
        let mut crop = GridF32::new(cw, cw, 0.5);
        for j in 0..cw {
            for i in 0..cw {
                if cx0 + i < w && cy0 + j < hh {
                    crop.data[j * cw + i] = hs.data[(cy0 + j) * w + (cx0 + i)];
                }
            }
        }
        let tag = if label.starts_with("raw") { "rawfbm" } else { "post" };
        crop.save_png_u8(&dir.join(format!("striation_stage_{tag}.png"))).unwrap();
    }
    eprintln!("   → crops: exports/sculpt/striation_stage_rawfbm.png / _post.png");
}

/// CLOSURE STEP 1 (verdict) — the terrace source is the FBM's FINEST OCTAVES. With 7
/// octaves + base_frequency 1.0 the finest octave is ~390 m (≈8 px @8192², ≈2 px
/// @2048²): sub-resolution at 2048² (aliased away → clean), fully resolved at 8192²
/// where the drainage carves it into a ~390 m comb. Reducing the octave count coarsens
/// the finest detail (→ 1.5–3 km) so the drainage template gives legible valleys at
/// 8192². Reports raw-FBM + post-relief-v2 steep share & striation for octaves 7/5/4 at
/// 8192², saves post crops. No pipeline change.
#[test]
#[ignore]
fn fbm_octave_ablation() {
    use std::path::Path;
    use std::time::Instant;
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let seed_u = 10481999410520546993u64;
    let (domain, base, res) = (400.0f32, ss.depth_scale_m as f32, 8192usize);
    let (state, _run) = coarse_state(seed_u);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(seed_u);
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../exports/sculpt");
    let (fx, fy, fw) = (0.22f32, 0.03f32, 0.34f32);
    let cell_km2 = (domain / res as f32).powi(2);

    let steep = |f: &GridF32, slope: &[f32]| -> (f32, f32) {
        let land: Vec<usize> = (0..f.data.len()).filter(|&k| f.data[k] > SEA).collect();
        let n = land.len().max(1) as f32;
        (
            land.iter().filter(|&&k| slope[k] > 30.0).count() as f32 / n * 100.0,
            land.iter().map(|&k| slope[k]).fold(0.0f32, f32::max),
        )
    };

    eprintln!("\n=== CLOSURE STEP 1 — FBM octave ablation @8192² (the terrace lever) ===");
    eprintln!("   octaves | finest λ | raw >30° | post >30° | post striation | post max°");
    for oct in [7usize, 5, 4] {
        let mut fc = FbmUpscaleConfig::c1_hd_production(res);
        fc.erosion = None;
        fc.bathymetry = None;
        fc.amplitude_base = 0.04;
        fc.octaves = oct;
        let raw = upscale_with_fbm(&coarse, SEA, &seed, &fc).heightmap;
        let (raw30, _) = steep(&raw, &slope_deg_field(&raw, domain, base));
        let t0 = Instant::now();
        let post = incise(&raw, &StreamPowerConfig::relief_v2(cell_km2, base));
        let ms = t0.elapsed().as_millis();
        let slope = slope_deg_field(&post, domain, base);
        let (p30, mx) = steep(&post, &slope);
        let sm: Vec<f32> = post.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).collect();
        let mut le: Vec<f32> =
            (0..post.data.len()).filter(|&k| post.data[k] > SEA).map(|k| sm[k]).collect();
        le.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let e60 = le[le.len() * 6 / 10];
        let upper: Vec<f32> =
            (0..post.data.len()).map(|k| if sm[k] >= e60 { slope[k] } else { 0.0 }).collect();
        let (_, pow) = striation_spectrum(&post, &upper, 20.0, 48);
        // finest octave wavelength in target px @8192: nscale=base_freq*1024/src_max²,
        // finest freq = nscale*lac^(oct-1) cyc/coarse-cell; coarse cell = res/64 px.
        let nscale = 1.0 * 1024.0 / (64.0f64 * 64.0);
        let finest_cyc = nscale * 2.0f64.powi(oct as i32 - 1);
        let lam_px = (res as f64 / 64.0) / finest_cyc;
        eprintln!(
            "   {oct:>7} | {lam_px:>5.0} px | {raw30:>7.1}% | {p30:>8.1}% | {pow:>13.2} | {mx:.0}  ({ms} ms)"
        );
        let hs = hillshade(&post, domain, base);
        let (w, hh) = (post.width, post.height);
        let (cx0, cy0, cw) =
            ((fx * w as f32) as usize, (fy * hh as f32) as usize, (fw * w as f32) as usize);
        let mut crop = GridF32::new(cw, cw, 0.5);
        for j in 0..cw {
            for i in 0..cw {
                if cx0 + i < w && cy0 + j < hh {
                    crop.data[j * cw + i] = hs.data[(cy0 + j) * w + (cx0 + i)];
                }
            }
        }
        crop.save_png_u8(&dir.join(format!("octave_{oct}_8192.png"))).unwrap();
    }
    eprintln!("   → crops: exports/sculpt/octave_{{7,5,4}}_8192.png");
    eprintln!("   (fewer octaves ⇒ coarser finest detail ⇒ the comb should clear at 8192²)");
}

/// CLOSURE STEP 1 (mechanism) — the 8192² comb is the SMITH–BRETHERTON parallel-rilling
/// instability of detachment-limited stream power (E=K·A^m·S^n, m<1 unstable on smooth
/// slopes), under-damped by hillslope diffusion at fine resolution. Its two physical
/// controls are the diffusion strength (sets the characteristic hillslope length / valley
/// spacing) and A_c (how much surface is channel vs diffusing hillslope). Sweeps both at
/// 8192², reports steep share + upper-slope striation, saves crops. If stronger diffusion
/// and/or larger A_c clear the comb → the fix is a resolution-scaled damping, not a
/// routing rewrite. No pipeline change.
#[test]
#[ignore]
fn rilling_sweep() {
    use std::path::Path;
    use std::time::Instant;
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let seed_u = 10481999410520546993u64;
    let (domain, base, res) = (400.0f32, ss.depth_scale_m as f32, 8192usize);
    let (state, _run) = coarse_state(seed_u);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(seed_u);
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../exports/sculpt");
    let (fx, fy, fw) = (0.22f32, 0.03f32, 0.34f32);
    let cell_km2 = (domain / res as f32).powi(2);
    let mut fc = FbmUpscaleConfig::c1_hd_production(res);
    fc.erosion = None;
    fc.bathymetry = None;
    fc.amplitude_base = 0.04;
    let raw = upscale_with_fbm(&coarse, SEA, &seed, &fc).heightmap;

    eprintln!("\n=== CLOSURE STEP 1 — rilling-instability sweep @8192² (diffusion × A_c) ===");
    eprintln!("   D_mult | A_c km² | >30° | upper-slope striation | max°");
    // (diffusion multiplier vs relief_v2 default 0.15, A_c km²).
    for (dmult, ac) in [(1.0f32, 0.1f32), (4.0, 0.1), (1.0, 0.5), (4.0, 0.5)] {
        let mut cfg = StreamPowerConfig::relief_v2(cell_km2, base);
        cfg.diffusion *= dmult;
        cfg.min_area_cells = ac / cell_km2;
        let t0 = Instant::now();
        let post = incise(&raw, &cfg);
        let ms = t0.elapsed().as_millis();
        let slope = slope_deg_field(&post, domain, base);
        let land: Vec<usize> = (0..post.data.len()).filter(|&k| post.data[k] > SEA).collect();
        let n = land.len().max(1) as f32;
        let a30 = land.iter().filter(|&&k| slope[k] > 30.0).count() as f32 / n * 100.0;
        let mx = land.iter().map(|&k| slope[k]).fold(0.0f32, f32::max);
        let sm: Vec<f32> = post.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).collect();
        let mut le: Vec<f32> = land.iter().map(|&k| sm[k]).collect();
        le.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let e60 = le[le.len() * 6 / 10];
        let upper: Vec<f32> =
            (0..post.data.len()).map(|k| if sm[k] >= e60 { slope[k] } else { 0.0 }).collect();
        let (_, pow) = striation_spectrum(&post, &upper, 20.0, 48);
        eprintln!("   {dmult:>6.1} | {ac:>7.2} | {a30:>4.1}% | {pow:>21.2} | {mx:.0}  ({ms} ms)");
        let hs = hillshade(&post, domain, base);
        let (w, hh) = (post.width, post.height);
        let (cx0, cy0, cw) =
            ((fx * w as f32) as usize, (fy * hh as f32) as usize, (fw * w as f32) as usize);
        let mut crop = GridF32::new(cw, cw, 0.5);
        for j in 0..cw {
            for i in 0..cw {
                if cx0 + i < w && cy0 + j < hh {
                    crop.data[j * cw + i] = hs.data[(cy0 + j) * w + (cx0 + i)];
                }
            }
        }
        crop.save_png_u8(&dir.join(format!("rill_d{dmult:.0}_ac{ac:.2}.png"))).unwrap();
    }
    eprintln!("   → crops: exports/sculpt/rill_d*_ac*.png");
    eprintln!(
        "   (comb clears with higher D and/or larger A_c ⇒ rilling instability, fixable by damping)"
    );
}

/// CLOSURE STEP 1 (clean-room) — discriminate REAL rilling (Smith–Bretherton, follows
/// the slope) from a GRID/solver artifact (follows grid axes), on a SMOOTH plane tilted
/// DIAGONALLY at 30°, no FBM, no coarse pipeline. If the rills run along the 30° downslope
/// direction → real physics (fix = damping/routing). If they snap to 0/45/90° in the grid
/// → D8 routing and/or the Gauss-Seidel sweep artifact (fix = the algorithm). Reports the
/// channel-segment orientation concentration in the GRID frame and RELATIVE to the known
/// 30° slope, for relief-v1 (linear symmetric diffusion) vs relief-v2 (GS diffusion), and
/// saves hillshade crops. 2048², fast.
#[test]
#[ignore]
fn synthetic_slope_rilling() {
    use std::path::Path;
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let (t, domain, base) = (2048usize, 400.0f32, ss.depth_scale_m as f32);
    let theta = 30.0f32.to_radians();
    let (cth, sth) = (theta.cos(), theta.sin());
    // Smooth plane descending along the 30° direction from 0.92 to 0.50, + tiny
    // deterministic hash noise (amp 0.002) to seed rilling without imposing a direction.
    let hash = |x: usize, y: usize| -> f32 {
        let mut h =
            (x as u32).wrapping_mul(374761393).wrapping_add((y as u32).wrapping_mul(668265263));
        h = (h ^ (h >> 13)).wrapping_mul(1274126177);
        ((h ^ (h >> 16)) as f32 / u32::MAX as f32) - 0.5
    };
    let diag_max = (t as f32) * (cth + sth);
    let mut d = vec![0.0f32; t * t];
    for y in 0..t {
        for x in 0..t {
            let proj = (x as f32 * cth + y as f32 * sth) / diag_max; // 0..1 down-slope
            d[y * t + x] = 0.92 - 0.42 * proj + 0.002 * hash(x, y);
        }
    }
    let plane = GridF32::from_vec(t, t, d);
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../exports/sculpt");
    let cell_km2 = (domain / t as f32).powi(2);

    let orient = |f: &GridF32| -> (f32, f32) {
        let (w, h) = (f.width, f.height);
        let dr = c1_drainage(f, None, &C1DrainageConfig::default(), &ss);
        let (mut grid, mut rel) = (Vec::new(), Vec::new());
        for s in &dr.rivers.segments {
            for win in s.points.windows(2) {
                let (dx, dy) = (
                    (win[1].0 as i32 - win[0].0 as i32) as f32,
                    (win[1].1 as i32 - win[0].1 as i32) as f32,
                );
                if dx == 0.0 && dy == 0.0 {
                    continue;
                }
                let seg = dy.atan2(dx);
                grid.push(seg);
                rel.push(seg - theta); // relative to the KNOWN downslope direction
            }
        }
        let _ = (w, h);
        (axial_concentration(&grid).0, axial_concentration(&rel).0)
    };

    eprintln!(
        "\n=== CLOSURE STEP 1 — synthetic 30°-tilted plane (real rilling vs grid artifact) ==="
    );
    eprintln!(
        "   config   | R_grid | R_slope-rel | (high R_slope-rel ⇒ rills follow the 30° slope = REAL;"
    );
    eprintln!("            |        |             |  high R_grid at 0/90° ⇒ grid/solver artifact)");
    for (label, v2) in [("relief-v1", false), ("relief-v2", true)] {
        let cfg = if v2 {
            StreamPowerConfig::relief_v2(cell_km2, base)
        } else {
            StreamPowerConfig::relief_v1(cell_km2, base)
        };
        let post = incise(&plane, &cfg);
        let (rg, rr) = orient(&post);
        eprintln!("   {label:<9} | {rg:>6.3} | {rr:>11.3}");
        hillshade(&post, domain, base)
            .save_png_u8(&dir.join(format!("synthplane_{label}.png")))
            .unwrap();
    }
    eprintln!(
        "   → crops: exports/sculpt/synthplane_relief-v{{1,2}}.png (30° slope; rills should run ↘)"
    );
}

/// STEP 2a — QUANTIFY the Smith–Bretherton stability criterion on the 30° plane, IN
/// METRES (before any fix code). Sweeps the hillslope diffusion `D` (applied EVERYWHERE
/// via `diffuse_channels`, the LEM-correct form; linear-like via a huge `S_c` so the
/// implicit solve stays stable at large D) and reports, at 2048² AND 8192²:
///   - rill wavelength (m) of the unstable comb — must be the SAME in metres at both
///     resolutions (else the resolution-invariance bug is back);
///   - the steep >30° share (rill presence) → `D_crit` = the minimal `D` that damps it;
///   - the emergent channel head `A*` (km²): the smallest drainage area still incised at
///     `D_crit` vs the largest area damped — the gap tells whether legitimate gorges
///     survive the damping (if min-incised-A > max-damped-A, they do).
/// Incision `K` is relief_v1's; lateral OFF (isolate incision vs diffusion). Read-only.
#[test]
#[ignore]
fn rill_stability() {
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let (domain, base) = (400.0f32, ss.depth_scale_m as f32);
    let theta = 30.0f32.to_radians();
    let (cth, sth) = (theta.cos(), theta.sin());
    let hash = |x: usize, y: usize| -> f32 {
        let mut h =
            (x as u32).wrapping_mul(374761393).wrapping_add((y as u32).wrapping_mul(668265263));
        h = (h ^ (h >> 13)).wrapping_mul(1274126177);
        ((h ^ (h >> 16)) as f32 / u32::MAX as f32) - 0.5
    };
    let build_plane = |t: usize| -> GridF32 {
        let diag_max = (t as f32) * (cth + sth);
        let mut d = vec![0.0f32; t * t];
        for y in 0..t {
            for x in 0..t {
                let proj = (x as f32 * cth + y as f32 * sth) / diag_max;
                d[y * t + x] = 0.92 - 0.42 * proj + 0.002 * hash(x, y);
            }
        }
        GridF32::from_vec(t, t, d)
    };

    eprintln!("\n=== STEP 2a — Smith–Bretherton stability on the 30° plane (IN METRES) ===");
    eprintln!("   res  |    D | >30° | rill λ (m) | κΔt (m²) | A* min-incised / max-damped (km²)");
    // 2048² only (fast): the physics + metre-invariance read here; the 8192² confirmation
    // waits until the Gauss-Seidel diffusion is parallelised (single-threaded GS × 240
    // sweeps × 67 M cells is minutes per D). Richer D sweep since it is cheap.
    for t in [2048usize] {
        let plane = build_plane(t);
        let cell_km = domain / t as f32;
        let cell_km2 = cell_km * cell_km;
        let cell_m = cell_km * 1000.0;
        let ds: &[f32] = &[0.0, 0.05, 0.15, 0.4, 1.0, 2.5, 5.0];
        for &dval in ds {
            let mut cfg = StreamPowerConfig::relief_v1(cell_km2, base);
            cfg.critical_slope = 100.0; // huge S_c ⇒ linear-like, but implicit ⇒ stable at any D
            cfg.diffuse_channels = true; // LEM-correct: diffuse everywhere
            cfg.lateral_erosion = 0.0;
            cfg.diffusion = dval;
            let post = incise(&plane, &cfg);
            let slope = slope_deg_field(&post, domain, base);
            let land: Vec<usize> = (0..post.data.len()).filter(|&k| post.data[k] > SEA).collect();
            let n = land.len().max(1) as f32;
            let a30 = land.iter().filter(|&&k| slope[k] > 30.0).count() as f32 / n * 100.0;
            // rill wavelength: contour-direction spectral peak on gentle cells → cells→m.
            let (lam_cells, _) = striation_spectrum(&post, &slope, 5.0, 64);
            let lam_m = lam_cells * cell_m;
            // A*: incised = lowered > 20 m vs the plane. Drainage area (km²) of each cell.
            let dr = c1_drainage(&post, None, &C1DrainageConfig::default(), &ss);
            let norm_to_m = 2.0 * 1.13 * base;
            let (mut min_inc_a, mut max_damp_a) = (f32::MAX, 0.0f32);
            for &k in &land {
                let lowered = (plane.data[k] - post.data[k]) * norm_to_m;
                let a_km2 = dr.flow.accumulation.data[k] * cell_km2;
                if lowered > 20.0 {
                    min_inc_a = min_inc_a.min(a_km2);
                } else {
                    max_damp_a = max_damp_a.max(a_km2);
                }
            }
            let kappa = dval * (400_000.0f32 / 2048.0).powi(2); // physical κΔt (m²)
            let mi = if min_inc_a == f32::MAX { -1.0 } else { min_inc_a };
            eprintln!(
                "  {t:>5} | {dval:>4.2} | {a30:>4.1}% | {lam_m:>10.0} | {kappa:>8.0} | {mi:>9.2} / {max_damp_a:.2}"
            );
        }
    }
    eprintln!(
        "   (D_crit = smallest D where >30% collapses to the planar baseline; rill λ in m must"
    );
    eprintln!(
        "    match across resolutions; gorges survive if min-incised-A stays > max-damped-A)"
    );
}

/// STEP 2b — REMEDY (i) on REAL terrain: relief-v2 + `diffuse_channels=true` (diffusion
/// everywhere, LEM-correct), short D sweep {0.15, 0.25, 0.40, 0.55} + D=1.0 as a
/// non-monotonicity control. The question is not "does the comb die at 0.4?" but "is
/// there a D that kills the comb WITHOUT killing the gorges?". Reports per D: >30° share
/// + upper-slope striation (comb), W/D per Strahler order + floor/ridge + crest curvature
/// (gorge survival — must keep widening downstream: ref v2 W/D 6.2→40.9, floor/ridge 0.32,
/// curv 288 m), Strahler histogram + confluences (health). Plus a GS CONVERGENCE check at
/// D=0.40 (40 vs 80 sweeps). Saves a massif crop per D. 2048², author seed. Read-only.
#[test]
#[ignore]
fn cross_rill_2b() {
    use std::path::Path;
    use std::time::Instant;
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let seed_u = 10481999410520546993u64;
    let (t, domain, base) = (2048usize, 400.0f32, ss.depth_scale_m as f32);
    let cell_km2 = (domain / t as f32).powi(2);
    let (state, _run) = coarse_state(seed_u);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(seed_u);
    let mut fc = FbmUpscaleConfig::c1_hd_production(t);
    fc.erosion = None;
    fc.bathymetry = None;
    fc.amplitude_base = 0.04;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fc).heightmap;
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../exports/sculpt");
    let (fx, fy, fw) = (0.22f32, 0.03f32, 0.34f32);
    let norm_to_m = 2.0 * 1.13 * base;

    // Measure comb + gorge-survival + health for one eroded field.
    let measure = |post: &GridF32, a_c: f32| -> (f32, f32, f32, f32, String, String) {
        let slope = slope_deg_field(post, domain, base);
        let sm: Vec<f32> = post.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).collect();
        let land: Vec<usize> = (0..post.data.len()).filter(|&k| post.data[k] > SEA).collect();
        let n = land.len().max(1) as f32;
        let a30 = land.iter().filter(|&&k| slope[k] > 30.0).count() as f32 / n * 100.0;
        let mut le: Vec<f32> = land.iter().map(|&k| sm[k]).collect();
        le.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let e60 = le[le.len() * 6 / 10];
        let upper: Vec<f32> =
            (0..post.data.len()).map(|k| if sm[k] >= e60 { slope[k] } else { 0.0 }).collect();
        let (_, striation) = striation_spectrum(post, &upper, 20.0, 48);
        // crest curvature on steep upper cells.
        let mut curv = Vec::new();
        for &k in &land {
            let (x, y) = (k % t, k / t);
            if sm[k] >= e60 && slope[k] > 30.0 && x > 0 && x < t - 1 && y > 0 && y < t - 1 {
                let lap = (post.data[k - 1] + post.data[k + 1] - 2.0 * post.data[k]).abs()
                    + (post.data[k - t] + post.data[k + t] - 2.0 * post.data[k]).abs();
                curv.push(lap * norm_to_m);
            }
        }
        curv.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let crest = if curv.is_empty() { 0.0 } else { curv[curv.len() / 2] };
        // floor/local-ridge on channel cells.
        let dr = c1_drainage(post, None, &C1DrainageConfig::default(), &ss);
        let chan: Vec<usize> = (0..post.data.len())
            .filter(|&k| dr.flow.accumulation.data[k] >= a_c && post.data[k] > SEA)
            .collect();
        let mut ratios = Vec::new();
        for &k in chan.iter().step_by(7) {
            let (x, y) = (k % t, k / t);
            let mut ridge = sm[k];
            for dy in -10i32..=10 {
                for dx in -10i32..=10 {
                    let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                    if nx >= 0 && ny >= 0 && nx < t as i32 && ny < t as i32 {
                        ridge = ridge.max(sm[ny as usize * t + nx as usize]);
                    }
                }
            }
            if ridge > 1.0 {
                ratios.push(sm[k] / ridge);
            }
        }
        ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let floor_ridge = if ratios.is_empty() { 0.0 } else { ratios[ratios.len() / 2] };
        // W/D per order + Strahler histogram + confluences.
        let wd = per_order_width_depth(post, domain, base, &ss);
        let wdstr =
            wd.iter().map(|(o, _, _, r)| format!("S{o} {r:.1}")).collect::<Vec<_>>().join(" ");
        let mut hist: std::collections::HashMap<u8, usize> = std::collections::HashMap::new();
        for s in &dr.rivers.segments {
            *hist.entry(s.strahler_order).or_default() += 1;
        }
        let mut ho: Vec<u8> = hist.keys().copied().collect();
        ho.sort();
        let confl = dr.rivers.segments.iter().filter(|s| s.upstream.len() >= 2).count();
        let hstr = format!(
            "{} |confl {confl}",
            ho.iter().map(|o| format!("S{o}:{}", hist[o])).collect::<Vec<_>>().join(" ")
        );
        (a30, striation, floor_ridge, crest, wdstr, hstr)
    };

    eprintln!(
        "\n=== STEP 2b — remedy (i) cross-rill diffusion on real terrain (2048², seed {seed_u}) ==="
    );
    eprintln!(
        "   SUCCESS: comb >30°→~1-2%; gorges: W/D keeps rising S1→S5 (ref 6.2→40.9), floor/ridge~0.32, curv~288m"
    );
    eprintln!(
        "   D    | >30° | striation | floor/ridge | crest curv | W/D per order | Strahler/confl"
    );
    for dval in [0.15f32, 0.25, 0.40, 0.55, 1.0] {
        let mut cfg = StreamPowerConfig::relief_v2(cell_km2, base);
        cfg.diffuse_channels = true;
        cfg.diffusion = dval;
        let t0 = Instant::now();
        let post = incise(&fbm, &cfg);
        let ms = t0.elapsed().as_millis();
        let (a30, stri, fr, crest, wd, hist) = measure(&post, cfg.min_area_cells);
        let tag = if dval > 0.55 { " (control)" } else { "" };
        eprintln!(
            "   {dval:.2} | {a30:>4.1}% | {stri:>9.2} | {fr:>11.2} | {crest:>7.0} m | {wd} | {hist}{tag}  ({ms}ms)"
        );
        let hs = hillshade(&post, domain, base);
        let (cx0, cy0, cw) =
            ((fx * t as f32) as usize, (fy * t as f32) as usize, (fw * t as f32) as usize);
        let mut crop = GridF32::new(cw, cw, 0.5);
        for j in 0..cw {
            for i in 0..cw {
                if cx0 + i < t && cy0 + j < t {
                    crop.data[j * cw + i] = hs.data[(cy0 + j) * t + (cx0 + i)];
                }
            }
        }
        crop.save_png_u8(&dir.join(format!("crossrill_d{dval:.2}.png"))).unwrap();
    }

    // GS convergence check at D=0.40: does doubling the sweeps change the answer?
    eprintln!("\n-- GS convergence @D=0.40: 40 vs 80 implicit sweeps --");
    let build = |iters: usize| -> GridF32 {
        let mut cfg = StreamPowerConfig::relief_v2(cell_km2, base);
        cfg.diffuse_channels = true;
        cfg.diffusion = 0.40;
        cfg.hillslope_implicit_iters = iters;
        incise(&fbm, &cfg)
    };
    let f40 = build(40);
    let f80 = build(80);
    let maxdiff =
        (0..f40.data.len()).map(|k| (f40.data[k] - f80.data[k]).abs()).fold(0.0f32, f32::max);
    let (a40, ..) = measure(&f40, StreamPowerConfig::relief_v2(cell_km2, base).min_area_cells);
    let (a80, ..) = measure(&f80, StreamPowerConfig::relief_v2(cell_km2, base).min_area_cells);
    eprintln!(
        "   max |Δnorm| 40↔80 = {maxdiff:.4} ({:.0} m); >30° 40={a40:.1}% 80={a80:.1}% → {}",
        maxdiff * norm_to_m,
        if maxdiff * norm_to_m < 15.0 {
            "CONVERGED"
        } else {
            "NOT converged — 40 sweeps insufficient"
        }
    );
    eprintln!(
        "   → crops: exports/sculpt/crossrill_d*.png. Recommend the LOWEST D meeting both criteria."
    );
}

/// Count land cells with any downhill neighbour steeper than `sc` (repose violations),
/// and the max slope. Used to test whether talus is a closure (bounded residual) or a
/// solver in disguise (unbounded passes).
fn slope_violations(field: &GridF32, cell_km: f32, depth: f32, sc: f32) -> (usize, f32) {
    use ymir_core::terrain::flow::{D8_DIST, D8_DX, D8_DY};
    let (w, h) = (field.width, field.height);
    let cell_m = cell_km * 1000.0;
    let norm_to_m = 2.0 * 1.13 * depth;
    let (mut viol, mut maxs) = (0usize, 0.0f32);
    for y in 0..h {
        for x in 0..w {
            let k = y * w + x;
            if field.data[k] <= SEA {
                continue;
            }
            let mut bad = false;
            for d in 0..8 {
                let (nx, ny) = (x as i32 + D8_DX[d], y as i32 + D8_DY[d]);
                if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 {
                    continue;
                }
                let j = ny as usize * w + nx as usize;
                let drop = field.data[k] - field.data[j];
                if drop <= 0.0 {
                    continue;
                }
                let s = drop * norm_to_m / (D8_DIST[d] * cell_m);
                if s > maxs {
                    maxs = s;
                }
                if s > sc * 1.02 {
                    bad = true;
                }
            }
            if bad {
                viol += 1;
            }
        }
    }
    (viol, maxs)
}

/// TASK 3 — is talus a CLOSURE or a solver in disguise? Run talus (post relief-v2
/// incision, everywhere) with 1/2/4/8 passes on real 2048² terrain and count the cells
/// still exceeding S_c after each. If residuals fall to ~0 in a few BOUNDED passes it is
/// a closure; if they need unbounded passes it is a solver. Read-only.
#[test]
#[ignore]
fn talus_residual() {
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let seed_u = 10481999410520546993u64;
    let (t, domain, base) = (2048usize, 400.0f32, ss.depth_scale_m as f32);
    let cell_km = domain / t as f32;
    let cell_km2 = cell_km * cell_km;
    let sc = 0.6494f32;
    let (state, _run) = coarse_state(seed_u);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(seed_u);
    let mut fc = FbmUpscaleConfig::c1_hd_production(t);
    fc.erosion = None;
    fc.bathymetry = None;
    fc.amplitude_base = 0.04;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fc).heightmap;
    let land = (0..fbm.data.len()).filter(|&k| fbm.data[k] > SEA).count().max(1);
    let (v0, mx0) = slope_violations(&fbm, cell_km, base, sc);
    eprintln!("\n=== TASK 3 — talus residual vs passes (2048², S_c=tan33°, real terrain) ===");
    eprintln!(
        "   pre-talus (relief-v2 incision only): {v0} violating cells ({:.1}% of land), max slope {:.2}",
        v0 as f32 / land as f32 * 100.0,
        mx0
    );
    eprintln!("   passes | violating cells | % land | max slope");
    for passes in [1usize, 2, 4, 8] {
        let mut cfg = StreamPowerConfig::relief_v2(cell_km2, base);
        cfg.critical_slope = 0.0; // OFF the nonlinear diffusion
        cfg.diffusion = 0.0;
        cfg.talus_slope = sc;
        cfg.talus_passes = passes;
        let post = incise(&fbm, &cfg);
        let (v, mx) = slope_violations(&post, cell_km, base, sc);
        eprintln!("   {passes:>6} | {v:>15} | {:>5.1}% | {mx:.2}", v as f32 / land as f32 * 100.0);
    }
    eprintln!(
        "   (residual → ~0 in a few bounded passes ⇒ CLOSURE; needs unbounded ⇒ solver in disguise)"
    );
}

/// TASK 2 + 2bis — TALUS vs NONLINEAR DIFFUSION head-to-head, at 2048² AND 8192². For
/// each method: comb (>30°, max slope, striation), gorge survival (W/D per order,
/// floor/ridge, crest curvature), drainage health (Strahler/confl), the HEADWATER
/// ramification discriminator (drainage density km/km², S1 segment count, channel-head
/// elevation as % of peak → "rivers start too low" as a number), A_c in km² AND cells,
/// and RUNTIME. The deciding property: does each method give the SAME RESULT IN METRES at
/// both resolutions? Saves a massif crop per method/res. Read-only.
#[test]
#[ignore]
fn talus_vs_diffusion() {
    use std::path::Path;
    use std::time::Instant;
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let seed_u = 10481999410520546993u64;
    let (domain, base) = (400.0f32, ss.depth_scale_m as f32);
    let sc = 0.6494f32;
    let (state, _run) = coarse_state(seed_u);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(seed_u);
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../exports/sculpt");
    let (fx, fy, fw) = (0.22f32, 0.03f32, 0.34f32);

    eprintln!("\n=== TASK 2/2bis — TALUS vs NONLINEAR DIFFUSION (seed {seed_u}) ===");
    for t in [2048usize, 8192] {
        let cell_km = domain / t as f32;
        let cell_km2 = cell_km * cell_km;
        let mut fc = FbmUpscaleConfig::c1_hd_production(t);
        fc.erosion = None;
        fc.bathymetry = None;
        fc.amplitude_base = 0.04;
        let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fc).heightmap;
        let fm: Vec<f32> = fbm.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).collect();
        let peak = fm.iter().cloned().fold(f32::MIN, f32::max);
        let a_c_cells = 0.1 / cell_km2;
        eprintln!(
            "\n-- {t}² (cell {:.0} m, A_c=0.1 km² = {:.1} cells, peak {peak:.0} m) --",
            cell_km * 1000.0,
            a_c_cells
        );
        eprintln!(
            "   method              | ms | >30° | max° | stri | floor/ridge | curv | dens km/km² | S1segs | head%peak | W/D S1→Smax"
        );

        let methods: Vec<(&str, StreamPowerConfig)> = vec![
            ("diffusion D=0.55", {
                let mut c = StreamPowerConfig::relief_v2(cell_km2, base);
                c.diffuse_channels = true;
                c.diffusion = 0.55;
                c
            }),
            ("diffusion D=1.00", {
                let mut c = StreamPowerConfig::relief_v2(cell_km2, base);
                c.diffuse_channels = true;
                c.diffusion = 1.00;
                c
            }),
            ("talus Sc=tan33", {
                let mut c = StreamPowerConfig::relief_v2(cell_km2, base);
                c.critical_slope = 0.0;
                c.diffusion = 0.0;
                c.talus_slope = sc;
                c.talus_passes = 4;
                c
            }),
            ("talus + lin D=0.15", {
                let mut c = StreamPowerConfig::relief_v2(cell_km2, base);
                c.critical_slope = 0.0;
                c.diffusion = 0.15;
                c.diffuse_channels = true;
                c.talus_slope = sc;
                c.talus_passes = 4;
                c
            }),
        ];

        for (label, cfg) in &methods {
            let t0 = Instant::now();
            let post = incise(&fbm, cfg);
            let ms = t0.elapsed().as_millis();
            let sm: Vec<f32> =
                post.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).collect();
            let slope = slope_deg_field(&post, domain, base);
            let land: Vec<usize> = (0..post.data.len()).filter(|&k| post.data[k] > SEA).collect();
            let nland = land.len().max(1);
            let a30 =
                land.iter().filter(|&&k| slope[k] > 30.0).count() as f32 / nland as f32 * 100.0;
            let (_, maxs) = slope_violations(&post, cell_km, base, sc);
            let maxdeg = maxs.atan().to_degrees();
            let mut le: Vec<f32> = land.iter().map(|&k| sm[k]).collect();
            le.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let e60 = le[le.len() * 6 / 10];
            let upper: Vec<f32> =
                (0..post.data.len()).map(|k| if sm[k] >= e60 { slope[k] } else { 0.0 }).collect();
            let (_, stri) = striation_spectrum(&post, &upper, 20.0, 48);
            // floor/ridge + curvature.
            let dr = c1_drainage(&post, None, &C1DrainageConfig::default(), &ss);
            let chan: Vec<usize> = (0..post.data.len())
                .filter(|&k| dr.flow.accumulation.data[k] >= a_c_cells && post.data[k] > SEA)
                .collect();
            let mut ratios = Vec::new();
            for &k in chan.iter().step_by(7) {
                let (x, y) = (k % t, k / t);
                let mut ridge = sm[k];
                for dy in -10i32..=10 {
                    for dx in -10i32..=10 {
                        let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                        if nx >= 0 && ny >= 0 && nx < t as i32 && ny < t as i32 {
                            ridge = ridge.max(sm[ny as usize * t + nx as usize]);
                        }
                    }
                }
                if ridge > 1.0 {
                    ratios.push(sm[k] / ridge);
                }
            }
            ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let fr = if ratios.is_empty() { 0.0 } else { ratios[ratios.len() / 2] };
            let norm_to_m = 2.0 * 1.13 * base;
            let mut curv = Vec::new();
            for &k in &land {
                let (x, y) = (k % t, k / t);
                if sm[k] >= e60 && slope[k] > 30.0 && x > 0 && x < t - 1 && y > 0 && y < t - 1 {
                    let lap = (post.data[k - 1] + post.data[k + 1] - 2.0 * post.data[k]).abs()
                        + (post.data[k - t] + post.data[k + t] - 2.0 * post.data[k]).abs();
                    curv.push(lap * norm_to_m);
                }
            }
            curv.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let crest = if curv.is_empty() { 0.0 } else { curv[curv.len() / 2] };
            // drainage density + S1 segments + channel-head elevation %peak.
            let dens = chan.len() as f32 * cell_km / (nland as f32 * cell_km2);
            let s1segs = dr.rivers.segments.iter().filter(|s| s.strahler_order == 1).count();
            let mut heads: Vec<f32> = dr
                .rivers
                .segments
                .iter()
                .filter(|s| s.upstream.is_empty())
                .filter_map(|s| s.points.first())
                .map(|&(x, y)| sm[y as usize * t + x as usize] / peak * 100.0)
                .collect();
            heads.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let head_p50 = if heads.is_empty() { 0.0 } else { heads[heads.len() / 2] };
            let wd = per_order_width_depth(&post, domain, base, &ss);
            let wdstr =
                wd.iter().map(|(o, _, _, r)| format!("S{o}:{r:.0}")).collect::<Vec<_>>().join(" ");
            eprintln!(
                "   {label:<19} | {ms:>5} | {a30:>4.1}% | {maxdeg:>3.0} | {stri:.2} | {fr:>11.2} | {crest:>4.0} | {dens:>11.2} | {s1segs:>6} | {head_p50:>8.0}% | {wdstr}"
            );
            let hs = hillshade(&post, domain, base);
            let (cx0, cy0, cw) =
                ((fx * t as f32) as usize, (fy * t as f32) as usize, (fw * t as f32) as usize);
            let mut crop = GridF32::new(cw, cw, 0.5);
            for j in 0..cw {
                for i in 0..cw {
                    if cx0 + i < t && cy0 + j < t {
                        crop.data[j * cw + i] = hs.data[(cy0 + j) * t + (cx0 + i)];
                    }
                }
            }
            let tag = label.replace(' ', "").replace('=', "").replace('.', "").replace('+', "");
            crop.save_png_u8(&dir.join(format!("tvd_{tag}_{t}.png"))).unwrap();
        }
    }
    eprintln!(
        "\n   METRE-INVARIANCE: compare each metric 2048² vs 8192² per method — the method whose"
    );
    eprintln!(
        "   numbers match in metres wins (closure property). H-A/H-B: does dens/S1/head reach rise"
    );
    eprintln!(
        "   at 8192² (resolution) or does talus keep more than diffusion at equal res (backfill)?"
    );
    eprintln!("   → crops: exports/sculpt/tvd_*_{{2048,8192}}.png");
}

/// MFD TASK 3 — the CLEAN-ROOM two-sided test on the smooth 30° plane. Sweep the MFD
/// partition exponent p (D8, then 10/4/2/1.1/1) with incision only (no diffusion, no
/// talus). Criterion is BOTH: (1) the comb must NOT appear (>30° low), AND (2) a main
/// channel must STILL emerge (max drainage area stays high — an MFD that disperses so
/// much no channel forms "passes" the comb test while destroying drainage). Saves plane
/// crops. Fast (2048²). Read-only.
#[test]
#[ignore]
fn mfd_plane_sweep() {
    use std::path::Path;
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let (t, domain, base) = (2048usize, 400.0f32, ss.depth_scale_m as f32);
    let cell_km2 = (domain / t as f32).powi(2);
    let theta = 30.0f32.to_radians();
    let (cth, sth) = (theta.cos(), theta.sin());
    let hash = |x: usize, y: usize| -> f32 {
        let mut h =
            (x as u32).wrapping_mul(374761393).wrapping_add((y as u32).wrapping_mul(668265263));
        h = (h ^ (h >> 13)).wrapping_mul(1274126177);
        ((h ^ (h >> 16)) as f32 / u32::MAX as f32) - 0.5
    };
    let diag_max = (t as f32) * (cth + sth);
    let mut d = vec![0.0f32; t * t];
    for y in 0..t {
        for x in 0..t {
            let proj = (x as f32 * cth + y as f32 * sth) / diag_max;
            d[y * t + x] = 0.92 - 0.42 * proj + 0.002 * hash(x, y);
        }
    }
    let plane = GridF32::from_vec(t, t, d);
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../exports/sculpt");
    let a_c = 0.1 / cell_km2;

    eprintln!(
        "\n=== MFD TASK 3 — 30° plane two-sided sweep (comb must die AND a channel must form) ==="
    );
    eprintln!("   p     | >30° | striation | max A (km²) | channel cells | verdict");
    let variants: [(&str, Option<f32>); 6] = [
        ("D8", None),
        ("10", Some(10.0)),
        ("4", Some(4.0)),
        ("2", Some(2.0)),
        ("1.1", Some(1.1)),
        ("1", Some(1.0)),
    ];
    for (label, p) in variants {
        let mut cfg = StreamPowerConfig::relief_v1(cell_km2, base);
        cfg.diffusion = 0.0; // isolate MFD's effect on rilling (no diffusion, no talus)
        cfg.mfd_exponent = p;
        let post = incise(&plane, &cfg);
        let slope = slope_deg_field(&post, domain, base);
        let land: Vec<usize> = (0..post.data.len()).filter(|&k| post.data[k] > SEA).collect();
        let n = land.len().max(1) as f32;
        let a30 = land.iter().filter(|&&k| slope[k] > 30.0).count() as f32 / n * 100.0;
        let (_, stri) = striation_spectrum(&post, &slope, 5.0, 64);
        // channel definition: max drainage area (from the SAME MFD field the incision used)
        // and the number of channel cells (A ≥ A_c).
        let facc = if let Some(pp) = p {
            let fl = ymir_core::terrain::flow::compute_flow(
                &post,
                &ymir_core::terrain::flow::FlowConfig { sea_level: SEA, ..Default::default() },
            );
            ymir_core::terrain::flow::mfd_accumulation(&fl.filled, &fl.direction, SEA, pp, t, t)
        } else {
            ymir_core::terrain::flow::compute_flow(
                &post,
                &ymir_core::terrain::flow::FlowConfig { sea_level: SEA, ..Default::default() },
            )
            .accumulation
        };
        let max_a = facc.data.iter().cloned().fold(0.0f32, f32::max) * cell_km2;
        let chan = (0..t * t).filter(|&k| facc.data[k] >= a_c && post.data[k] > SEA).count();
        let ok_comb = a30 < 6.0;
        let ok_chan = max_a > 50.0; // a real trunk still concentrates
        let verdict = match (ok_comb, ok_chan) {
            (true, true) => "PASS (comb dead, channel lives)",
            (true, false) => "over-dispersed (no channel)",
            (false, true) => "comb survives",
            (false, false) => "worst of both",
        };
        eprintln!(
            "   {label:<5} | {a30:>4.1}% | {stri:>9.2} | {max_a:>11.0} | {chan:>13} | {verdict}"
        );
        hillshade(&post, domain, base)
            .save_png_u8(&dir.join(format!("mfdplane_p{label}.png")))
            .unwrap();
    }
    eprintln!("   → crops: exports/sculpt/mfdplane_p*.png (30° slope; want no comb + a trunk)");
}

/// MFD TASK 4 + 5 — real terrain at 2048² AND 8192². The recommended MFD incision config
/// turns OFF the nonlinear-diffusion SOLVER (critical_slope=0) — MFD prevents the comb, so
/// the expensive Gauss-Seidel is unnecessary — keeping only light linear hillslope + lateral
/// widening. Sweeps p (4/2/1.1) and K (×1/×2, since MFD disperses A → weaker incision) at
/// 2048², runs the recommended pair at 8192². Reports the full metric set + drainage density
/// / S1 / head-reach + RUNTIME, plus TASK 5 D8/MFD ALIGNMENT (are the D8-traced rivers in the
/// MFD-carved valley floors?). Metre-invariance = compare 2048² vs 8192². Read-only.
#[test]
#[ignore]
fn mfd_real() {
    use std::path::Path;
    use std::time::Instant;
    use ymir_core::erosion::stream_power::{RELIEF_V1_K, StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let seed_u = 10481999410520546993u64;
    let (domain, base) = (400.0f32, ss.depth_scale_m as f32);
    let (state, _run) = coarse_state(seed_u);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(seed_u);
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../exports/sculpt");
    let (fx, fy, fw) = (0.22f32, 0.03f32, 0.34f32);
    let norm_to_m = 2.0 * 1.13 * base;

    let mut report = |t: usize, label: &str, p: f32, kmult: f32, fbm: &GridF32, peak: f32| {
        let cell_km = domain / t as f32;
        let cell_km2 = cell_km * cell_km;
        let a_c = 0.1 / cell_km2;
        let mut cfg = StreamPowerConfig::relief_v2(cell_km2, base);
        cfg.critical_slope = 0.0; // MFD prevents the comb → drop the GS solver
        cfg.diffusion = 0.05; // light linear hillslope only
        cfg.diffuse_channels = true;
        cfg.k = RELIEF_V1_K * kmult;
        cfg.mfd_exponent = Some(p);
        let t0 = Instant::now();
        let post = incise(fbm, &cfg);
        let ms = t0.elapsed().as_millis();
        let slope = slope_deg_field(&post, domain, base);
        let sm: Vec<f32> = post.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).collect();
        let land: Vec<usize> = (0..post.data.len()).filter(|&k| post.data[k] > SEA).collect();
        let nland = land.len().max(1);
        let a30 = land.iter().filter(|&&k| slope[k] > 30.0).count() as f32 / nland as f32 * 100.0;
        let (_, maxs) = slope_violations(&post, cell_km, base, 0.6494);
        let maxdeg = maxs.atan().to_degrees();
        let mut le: Vec<f32> = land.iter().map(|&k| sm[k]).collect();
        le.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let e60 = le[le.len() * 6 / 10];
        let upper: Vec<f32> =
            (0..post.data.len()).map(|k| if sm[k] >= e60 { slope[k] } else { 0.0 }).collect();
        let (_, stri) = striation_spectrum(&post, &upper, 20.0, 48);
        // D8 drainage (the exported rivers) on the carved terrain.
        let dr = c1_drainage(&post, None, &C1DrainageConfig::default(), &ss);
        let chan: Vec<usize> = (0..post.data.len())
            .filter(|&k| dr.flow.accumulation.data[k] >= a_c && post.data[k] > SEA)
            .collect();
        let mut ratios = Vec::new();
        for &k in chan.iter().step_by(7) {
            let (x, y) = (k % t, k / t);
            let mut ridge = sm[k];
            for dy in -10i32..=10 {
                for dx in -10i32..=10 {
                    let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                    if nx >= 0 && ny >= 0 && nx < t as i32 && ny < t as i32 {
                        ridge = ridge.max(sm[ny as usize * t + nx as usize]);
                    }
                }
            }
            if ridge > 1.0 {
                ratios.push(sm[k] / ridge);
            }
        }
        ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let fr = if ratios.is_empty() { 0.0 } else { ratios[ratios.len() / 2] };
        let dens = chan.len() as f32 * cell_km / (nland as f32 * cell_km2);
        let s1 = dr.rivers.segments.iter().filter(|s| s.strahler_order == 1).count();
        let confl = dr.rivers.segments.iter().filter(|s| s.upstream.len() >= 2).count();
        let maxo = dr.rivers.segments.iter().map(|s| s.strahler_order).max().unwrap_or(0);
        let mut heads: Vec<f32> = dr
            .rivers
            .segments
            .iter()
            .filter(|s| s.upstream.is_empty())
            .filter_map(|s| s.points.first())
            .map(|&(x, y)| sm[y as usize * t + x as usize] / peak * 100.0)
            .collect();
        heads.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let head_p50 = if heads.is_empty() { 0.0 } else { heads[heads.len() / 2] };
        let wd = per_order_width_depth(&post, domain, base, &ss);
        let wdstr =
            wd.iter().map(|(o, _, _, r)| format!("S{o}:{r:.0}")).collect::<Vec<_>>().join(" ");
        // TASK 5 — D8/MFD alignment: for D8 river cells, vertical offset above the local
        // minimum (±3) of the carved terrain; share NOT sitting in a hollow (offset > 20 m).
        let mut offs = Vec::new();
        let mut off_hi = 0usize;
        for s in &dr.rivers.segments {
            for &(x, y) in &s.points {
                let (x, y) = (x as usize, y as usize);
                let mut lo = sm[y * t + x];
                for dy in -3i32..=3 {
                    for dx in -3i32..=3 {
                        let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                        if nx >= 0 && ny >= 0 && nx < t as i32 && ny < t as i32 {
                            lo = lo.min(sm[ny as usize * t + nx as usize]);
                        }
                    }
                }
                let off = sm[y * t + x] - lo;
                offs.push(off);
                if off > 20.0 {
                    off_hi += 1;
                }
            }
        }
        offs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let off_p50 = if offs.is_empty() { 0.0 } else { offs[offs.len() / 2] };
        let off_share =
            if offs.is_empty() { 0.0 } else { off_hi as f32 / offs.len() as f32 * 100.0 };
        eprintln!(
            "  {t} {label:<12} {ms:>6}ms | >30 {a30:>4.1}% max {maxdeg:>2.0}° stri {stri:.2} | fr {fr:.2} dens {dens:.2} S1 {s1} head {head_p50:.0}% Smax {maxo} confl {confl} | align p50 {off_p50:.0}m off>{off_share:.0}% | W/D {wdstr}"
        );
        let _ = norm_to_m;
        let hs = hillshade(&post, domain, base);
        let (cx0, cy0, cw) =
            ((fx * t as f32) as usize, (fy * t as f32) as usize, (fw * t as f32) as usize);
        let mut crop = GridF32::new(cw, cw, 0.5);
        for j in 0..cw {
            for i in 0..cw {
                if cx0 + i < t && cy0 + j < t {
                    crop.data[j * cw + i] = hs.data[(cy0 + j) * t + (cx0 + i)];
                }
            }
        }
        let tag = label.replace(' ', "").replace('=', "").replace('.', "").replace('×', "x");
        crop.save_png_u8(&dir.join(format!("mfdreal_{tag}_{t}.png"))).unwrap();
    };

    eprintln!("\n=== MFD TASK 4/5 — real terrain (seed {seed_u}, MFD incision, no GS solver) ===");
    for t in [2048usize, 8192] {
        let mut fc = FbmUpscaleConfig::c1_hd_production(t);
        fc.erosion = None;
        fc.bathymetry = None;
        fc.amplitude_base = 0.04;
        let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fc).heightmap;
        let peak =
            fbm.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).fold(f32::MIN, f32::max);
        eprintln!("-- {t}² (peak {peak:.0} m) --");
        if t == 2048 {
            report(t, "p4 Kx1", 4.0, 1.0, &fbm, peak);
            report(t, "p2 Kx1", 2.0, 1.0, &fbm, peak);
            report(t, "p2 Kx2", 2.0, 2.0, &fbm, peak);
            report(t, "p1.1 Kx2", 1.1, 2.0, &fbm, peak);
        } else {
            report(t, "p2 Kx2", 2.0, 2.0, &fbm, peak);
        }
    }
    eprintln!("   (want: comb dead + fr NOT climbing to 0.80 (no under-incision) + W/D widening +");
    eprintln!(
        "    alignment offset small = D8 rivers sit in MFD valleys; metre-invariance 2048↔8192)"
    );
    eprintln!("   → crops: exports/sculpt/mfdreal_*.png");
}

/// Median floor/local-ridge over channel cells, with a PHYSICAL ridge window (`window_km`)
/// instead of a fixed ±10-CELL window — the fixed-cell window measures a 4× smaller ridge
/// in metres at 8192² than 2048², inflating the ratio (the cells-not-metres bug, 3rd time).
fn floor_ridge_phys(
    sm: &[f32],
    acc: &GridF32,
    field: &GridF32,
    a_c: f32,
    t: usize,
    cell_km: f32,
    window_km: f32,
) -> f32 {
    let r = (window_km / cell_km).round() as i32;
    let mut ratios = Vec::new();
    for k in (0..t * t).filter(|&k| acc.data[k] >= a_c && field.data[k] > SEA).step_by(7) {
        let (x, y) = (k % t, k / t);
        let mut ridge = sm[k];
        for dy in -r..=r {
            for dx in -r..=r {
                let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                if nx >= 0 && ny >= 0 && nx < t as i32 && ny < t as i32 {
                    ridge = ridge.max(sm[ny as usize * t + nx as usize]);
                }
            }
        }
        if ridge > 1.0 {
            ratios.push(sm[k] / ridge);
        }
    }
    ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if ratios.is_empty() { 0.0 } else { ratios[ratios.len() / 2] }
}

/// MFD TASK 1 — calibrate K on floor/local-ridge with a PHYSICAL ±1 km ridge window, at
/// 2048² AND 8192². First shows the fixed-cell window is the source of the apparent
/// resolution dependence (0.25 vs 0.72), then sweeps K to hit floor/ridge ≈ 0.4–0.5 at
/// 8192². Reports per K: floor/ridge (phys), per-order incision, W/D per order, >30°,
/// crest curvature, striation. Read-only.
#[test]
#[ignore]
fn mfd_k_sweep() {
    use std::time::Instant;
    use ymir_core::erosion::stream_power::{RELIEF_V1_K, StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let seed_u = 10481999410520546993u64;
    let (domain, base) = (400.0f32, ss.depth_scale_m as f32);
    let (state, _run) = coarse_state(seed_u);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(seed_u);

    eprintln!("\n=== MFD TASK 1 — K calibration on PHYSICAL floor/ridge (±1 km window) ===");
    // First: the metric artifact. Same field, ±10-cell vs ±1 km window, both resolutions.
    eprintln!("-- metric-window check (MFD p2 K×2): fixed ±10 cells vs physical ±1 km --");
    for t in [2048usize, 8192] {
        let cell_km = domain / t as f32;
        let cell_km2 = cell_km * cell_km;
        let mut fc = FbmUpscaleConfig::c1_hd_production(t);
        fc.erosion = None;
        fc.bathymetry = None;
        fc.amplitude_base = 0.04;
        let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fc).heightmap;
        let mut cfg = StreamPowerConfig::relief_v2(cell_km2, base);
        cfg.critical_slope = 0.0;
        cfg.diffusion = 0.05;
        cfg.diffuse_channels = true;
        cfg.k = RELIEF_V1_K * 2.0;
        cfg.mfd_exponent = Some(2.0);
        let post = incise(&fbm, &cfg);
        let sm: Vec<f32> = post.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).collect();
        let dr = c1_drainage(&post, None, &C1DrainageConfig::default(), &ss);
        let a_c = 0.1 / cell_km2;
        // fixed ±10-cell ratio (old metric):
        let r10 =
            floor_ridge_phys(&sm, &dr.flow.accumulation, &post, a_c, t, cell_km, 10.0 * cell_km);
        let rphys = floor_ridge_phys(&sm, &dr.flow.accumulation, &post, a_c, t, cell_km, 1.0);
        eprintln!("   {t}²: ±10 cells = {r10:.2}  |  ±1 km (phys) = {rphys:.2}");
    }

    eprintln!("-- K sweep at 8192² (MFD p2), target physical floor/ridge ≈ 0.4–0.5 --");
    eprintln!("   K×   | ms | floor/ridge(1km) | >30° | curv | per-order incision | W/D S1→S5");
    let t = 8192usize;
    let cell_km = domain / t as f32;
    let cell_km2 = cell_km * cell_km;
    let mut fc = FbmUpscaleConfig::c1_hd_production(t);
    fc.erosion = None;
    fc.bathymetry = None;
    fc.amplitude_base = 0.04;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fc).heightmap;
    let a_c = 0.1 / cell_km2;
    for kmult in [2.0f32, 3.0, 4.0, 6.0] {
        let mut cfg = StreamPowerConfig::relief_v2(cell_km2, base);
        cfg.critical_slope = 0.0;
        cfg.diffusion = 0.05;
        cfg.diffuse_channels = true;
        cfg.k = RELIEF_V1_K * kmult;
        cfg.mfd_exponent = Some(2.0);
        let t0 = Instant::now();
        let post = incise(&fbm, &cfg);
        let ms = t0.elapsed().as_millis();
        let sm: Vec<f32> = post.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).collect();
        let slope = slope_deg_field(&post, domain, base);
        let land: Vec<usize> = (0..post.data.len()).filter(|&k| post.data[k] > SEA).collect();
        let n = land.len().max(1) as f32;
        let a30 = land.iter().filter(|&&k| slope[k] > 30.0).count() as f32 / n * 100.0;
        let dr = c1_drainage(&post, None, &C1DrainageConfig::default(), &ss);
        let fr = floor_ridge_phys(&sm, &dr.flow.accumulation, &post, a_c, t, cell_km, 1.0);
        let norm_to_m = 2.0 * 1.13 * base;
        let mut le: Vec<f32> = land.iter().map(|&k| sm[k]).collect();
        le.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let e60 = le[le.len() * 6 / 10];
        let mut curv = Vec::new();
        for &k in &land {
            let (x, y) = (k % t, k / t);
            if sm[k] >= e60 && slope[k] > 30.0 && x > 0 && x < t - 1 && y > 0 && y < t - 1 {
                let lap = (post.data[k - 1] + post.data[k + 1] - 2.0 * post.data[k]).abs()
                    + (post.data[k - t] + post.data[k + t] - 2.0 * post.data[k]).abs();
                curv.push(lap * norm_to_m);
            }
        }
        curv.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let crest = if curv.is_empty() { 0.0 } else { curv[curv.len() / 2] };
        let (tab, _) = per_order_incision(&fbm, &post, &ss);
        let wd = per_order_width_depth(&post, domain, base, &ss);
        let wdstr =
            wd.iter().map(|(o, _, _, r)| format!("S{o}:{r:.0}")).collect::<Vec<_>>().join(" ");
        eprintln!(
            "   ×{kmult:<3.0} | {ms:>5}ms | {fr:>15.2} | {a30:>4.1}% | {crest:>4.0} | {} | {wdstr}",
            fmt_orders(&tab)
        );
    }
    eprintln!(
        "   (pick the LOWEST K giving floor/ridge ≈ 0.4–0.5 — deep valleys with legible flanks)"
    );
}

/// MFD TASK 2/3 — RIVER MONOTONICITY: the exported long profile of every river segment must
/// decrease toward the sea (water cannot climb). Runs c1_drainage on the MFD-incised field
/// (AFTER incision — see production_upscale.rs:402 then the drainage call), walks each
/// segment's points, and flags any REAL-elevation climb downstream. Attributes each climb to
/// H2/pit-fill (the climbing cell sits on a pit-FILLED flat: filled > real) vs H1/misalignment
/// (the cell is above the local terrain minimum). The permanent acceptance test. Read-only.
#[test]
#[ignore]
fn river_monotonicity() {
    use ymir_core::erosion::stream_power::{RELIEF_V1_K, StreamPowerConfig, incise};
    use ymir_core::terrain::flow::{FlowConfig, compute_flow};
    let ss = SteinSteinParams::default();
    let seed_u = 10481999410520546993u64;
    let (domain, base) = (400.0f32, ss.depth_scale_m as f32);
    let (state, _run) = coarse_state(seed_u);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(seed_u);
    let norm_to_m = 2.0 * 1.13 * base;

    eprintln!("\n=== MFD TASK 2/3 — river long-profile monotonicity (8192², MFD p2 K×3) ===");
    let t = 8192usize;
    let cell_km2 = (domain / t as f32).powi(2);
    let mut fc = FbmUpscaleConfig::c1_hd_production(t);
    fc.erosion = None;
    fc.bathymetry = None;
    fc.amplitude_base = 0.04;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fc).heightmap;
    let mut cfg = StreamPowerConfig::relief_v2(cell_km2, base);
    cfg.critical_slope = 0.0;
    cfg.diffusion = 0.05;
    cfg.diffuse_channels = true;
    cfg.k = RELIEF_V1_K * 3.0;
    cfg.mfd_exponent = Some(2.0);
    let post = incise(&fbm, &cfg);
    let sm: Vec<f32> = post.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).collect();
    let dr = c1_drainage(&post, None, &C1DrainageConfig::default(), &ss);
    // pit-filled surface (to attribute climbs to filled flats).
    let flow = compute_flow(&post, &FlowConfig { sea_level: SEA, ..Default::default() });
    let is_filled = |k: usize| flow.filled.data[k] > post.data[k] + 1e-4;

    let (mut segs, mut viol, mut worst, mut pitfill, mut misalign) =
        (0usize, 0usize, 0.0f32, 0usize, 0usize);
    for s in &dr.rivers.segments {
        segs += 1;
        let mut bad = false;
        for w in s.points.windows(2) {
            let a = sm[w[0].1 as usize * t + w[0].0 as usize];
            let b = sm[w[1].1 as usize * t + w[1].0 as usize];
            let climb = b - a; // downstream elevation change; > 0 = uphill
            if climb > worst {
                worst = climb;
            }
            if climb > 1.0 {
                bad = true;
                let dk = w[1].1 as usize * t + w[1].0 as usize;
                if is_filled(dk) {
                    pitfill += 1;
                } else {
                    misalign += 1;
                }
            }
        }
        if bad {
            viol += 1;
        }
    }
    eprintln!(
        "   segments {segs}, WITH an uphill step: {viol} ({:.1}%), worst climb {worst:.0} m",
        viol as f32 / segs.max(1) as f32 * 100.0
    );
    eprintln!(
        "   climbing steps on pit-FILLED flats (H2/fill): {pitfill}  |  on real terrain (H1/misalign): {misalign}"
    );
    eprintln!(
        "   → verdict: {}",
        if pitfill > misalign {
            "H2-like — rivers cross filled depressions (real floor dips then climbs to the sill)"
        } else {
            "H1-like — D8 segment rides the flank off the MFD thalweg"
        }
    );
    let _ = norm_to_m;
}

#[test]
#[ignore]
fn closure_mosaic() {
    use std::path::Path;
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../exports/sculpt");
    const P: usize = 640; // panel edge
    const S: usize = 12; // separator width
    // Normalised crop window over the MASSIF (upper-left of this seed; the centre is a
    // flat plateau). Same fractional window at every resolution → matched framing.
    let (fx, fy, fw) = (0.22f32, 0.03f32, 0.34f32);
    // Crop the massif window from a full render, box-resample to P×P, then contrast-
    // stretch (hillshade sits near 0.5±small → stretch to [0,1] so structure is visible).
    let to_panel = |g: &GridF32| -> GridF32 {
        let (cx0, cy0) = ((fx * g.width as f32) as usize, (fy * g.height as f32) as usize);
        let cw = (fw * g.width as f32) as usize;
        let mut out = GridF32::new(P, P, 0.5);
        let s = cw as f32 / P as f32;
        for j in 0..P {
            for i in 0..P {
                let (x0, y0) = (cx0 + (i as f32 * s) as usize, cy0 + (j as f32 * s) as usize);
                let (x1, y1) = (
                    (cx0 + ((i + 1) as f32 * s) as usize).min(g.width),
                    (cy0 + ((j + 1) as f32 * s) as usize).min(g.height),
                );
                let (mut acc, mut cnt) = (0.0f32, 0.0f32);
                for y in y0..y1.max(y0 + 1) {
                    for x in x0..x1.max(x0 + 1) {
                        if x < g.width && y < g.height {
                            acc += g.data[y * g.width + x];
                            cnt += 1.0;
                        }
                    }
                }
                out.data[j * P + i] = if cnt > 0.0 { acc / cnt } else { 0.5 };
            }
        }
        // Contrast stretch to the 2–98 percentile.
        let mut v = out.data.clone();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let (lo, hi) = (v[v.len() * 2 / 100], v[v.len() * 98 / 100]);
        let d = (hi - lo).max(1e-4);
        for p in out.data.iter_mut() {
            *p = ((*p - lo) / d).clamp(0.0, 1.0);
        }
        out
    };
    // (row, col, full-render filename) — None cell stays blank/white.
    let cells: [(usize, usize, Option<&str>); 6] = [
        (0, 0, Some("closure_v1_2048.png")),
        (0, 1, Some("closure_v2_2048.png")),
        (0, 2, Some("closure_v2_2048_amp01.png")),
        (1, 0, None),
        (1, 1, Some("closure_v2_8192.png")),
        (1, 2, Some("closure_v2_8192_amp02.png")),
    ];
    let (cols, rows) = (3usize, 2usize);
    let (mw, mh) = (cols * P + (cols - 1) * S, rows * P + (rows - 1) * S);
    let mut mosaic = GridF32::new(mw, mh, 1.0); // white background = separators
    for (r, c, name) in cells {
        let Some(name) = name else { continue };
        let path = dir.join(name);
        let panel = match GridF32::load_png(&path) {
            Ok(g) => to_panel(&g),
            Err(e) => {
                eprintln!("  missing {name}: {e} (run closure_render first)");
                continue;
            }
        };
        let (ox, oy) = (c * (P + S), r * (P + S));
        for j in 0..P {
            for i in 0..P {
                mosaic.data[(oy + j) * mw + (ox + i)] = panel.data[j * P + i];
            }
        }
    }
    let out = dir.join("closure_mosaic.png");
    mosaic.save_png_u8(&out).unwrap();
    eprintln!("\n=== CLOSURE mosaic (massif window, contrast-stretched) → {} ===", out.display());
    eprintln!("  row 0: v1 2048 | v2 2048 | v2 2048 amp0.01");
    eprintln!("  row 1:  (blank) | v2 8192 | v2 8192 amp0.02");
}

/// CLOSURE STEP 1 (read-only) — confirm the missing bounding closure. Reports the
/// land slope distribution (share > 30/35/40/45°, max) for (a) the raw FBM, (b) the
/// relief_v1 sculpt, at the review amp 0.04 AND the default amp 0.16, plus the >30°
/// flank contiguity. If slopes far exceed a plausible angle of repose (~33°), the
/// closure that bounds them (nonlinear hillslope diffusion with critical slope) is
/// missing — the finding that justifies the chantier. 2048², domain 400 km, no
/// pipeline change.
#[test]
#[ignore]
fn step1_slope_distribution() {
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let seed_u = 10481999410520546993u64;
    let (t, domain) = (2048usize, 400.0f32);
    let base = ss.depth_scale_m as f32;
    let cell_km2 = (domain / t as f32).powi(2);
    let (state, _run) = coarse_state(seed_u);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(seed_u);

    let bins = |slope: &[f32], field: &GridF32| -> (f32, f32, f32, f32, f32) {
        let land: Vec<f32> =
            (0..field.data.len()).filter(|&k| field.data[k] > SEA).map(|k| slope[k]).collect();
        let n = land.len().max(1) as f32;
        let frac = |thr: f32| land.iter().filter(|&&s| s > thr).count() as f32 / n * 100.0;
        let mx = land.iter().cloned().fold(0.0f32, f32::max);
        (frac(30.0), frac(35.0), frac(40.0), frac(45.0), mx)
    };

    eprintln!(
        "\n=== CLOSURE STEP 1 — land slope distribution (seed {seed_u}, {t}², {domain} km) ==="
    );
    eprintln!(
        "  stage           amp | >30° | >35° | >40° | >45° |  max° | >30° largest flank (cells)"
    );
    for amp in [0.16f32, 0.04] {
        let mut fc = FbmUpscaleConfig::c1_hd_production(t);
        fc.erosion = None;
        fc.bathymetry = None;
        fc.amplitude_base = amp as f64;
        let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fc).heightmap;
        let cfg = StreamPowerConfig::relief_v1(cell_km2, base);
        let sp = incise(&fbm, &cfg);
        for (label, field) in [("FBM (raw)", &fbm), ("relief_v1 sculpt", &sp)] {
            let slope = slope_deg_field(field, domain, base);
            let (a30, a35, a40, a45, mx) = bins(&slope, field);
            let (_, big, _) = flank_contiguity(&slope, field.width, field.height);
            eprintln!(
                "  {label:<16}{amp:>4.2} | {a30:>4.1} | {a35:>4.1} | {a40:>4.1} | {a45:>4.1} | {mx:>5.1} | {big}",
            );
        }
    }
    eprintln!("  (angle of repose ~33°; slopes far above it with no clustering into arêtes = the");
    eprintln!(
        "   bounding closure — nonlinear hillslope diffusion with critical slope — is absent)"
    );
}

/// CLOSURE STEP 4 — renders for the visual verdict. v1 (no closures) vs v2 (both
/// closures) at amp 0.04, plus v2 at amp 0.02 and 0.01 (does structure now come from
/// the closures, so FBM can drop further?), plus one v2 at 8192². Hillshade + centre
/// crop each, into exports/sculpt/closure_*. Reports max slope / >30% share / relief /
/// trunk W/D so the numbers accompany the picture.
#[test]
#[ignore]
fn closure_render() {
    use std::path::Path;
    use std::time::Instant;
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let seed_u = 10481999410520546993u64;
    let domain = 400.0f32;
    let base = ss.depth_scale_m as f32;
    let (state, _run) = coarse_state(seed_u);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(seed_u);
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../exports/sculpt");
    std::fs::create_dir_all(&dir).unwrap();

    let mut render = |name: &str, t: usize, amp: f64, v2: bool| {
        let cell_km2 = (domain / t as f32).powi(2);
        let mut fc = FbmUpscaleConfig::c1_hd_production(t);
        fc.erosion = None;
        fc.bathymetry = None;
        fc.amplitude_base = amp;
        let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fc).heightmap;
        let cfg = if v2 {
            StreamPowerConfig::relief_v2(cell_km2, base)
        } else {
            StreamPowerConfig::relief_v1(cell_km2, base)
        };
        let t0 = Instant::now();
        let sp = incise(&fbm, &cfg);
        let ms = t0.elapsed().as_millis();
        let hs = hillshade(&sp, domain, base);
        hs.save_png_u8(&dir.join(format!("closure_{name}.png"))).unwrap();
        let (w, hh) = (sp.width, sp.height);
        let (x0, y0, cw, ch) = (w / 2 - w / 8, hh / 2 - hh / 8, w / 4, hh / 4);
        let mut crop = GridF32::new(cw, ch, 0.0);
        for j in 0..ch {
            for i in 0..cw {
                crop.data[j * cw + i] = hs.data[(y0 + j) * w + (x0 + i)];
            }
        }
        crop.save_png_u8(&dir.join(format!("closure_{name}_crop.png"))).unwrap();
        let slope = slope_deg_field(&sp, domain, base);
        let land: Vec<usize> = (0..sp.data.len()).filter(|&k| sp.data[k] > SEA).collect();
        let n = land.len().max(1) as f32;
        let a30 = land.iter().filter(|&&k| slope[k] > 30.0).count() as f32 / n * 100.0;
        let mx = land.iter().map(|&k| slope[k]).fold(0.0f32, f32::max);
        let relief = drainage_relief_m(&sp, &ss);
        let wd = per_order_width_depth(&sp, domain, base, &ss);
        let trunk = wd
            .iter()
            .rev()
            .take(2)
            .map(|(o, w, d, r)| format!("S{o} W/D {r:.1}"))
            .collect::<Vec<_>>()
            .join(" ");
        eprintln!(
            "  {name:<16} {t}² amp{amp:.2} {ms:>6}ms | >30° {a30:>4.1}% max {mx:.0}° relief {relief:>3.0}m | {trunk}"
        );
    };

    eprintln!("\n=== CLOSURE STEP 4 — renders (seed {seed_u}) → exports/sculpt/closure_* ===");
    render("v1_2048", 2048, 0.04, false);
    render("v2_2048", 2048, 0.04, true);
    render("v2_2048_amp02", 2048, 0.02, true);
    render("v2_2048_amp01", 2048, 0.01, true);
    render("v2_8192", 8192, 0.04, true);
    render("v2_8192_amp02", 8192, 0.02, true);
    eprintln!(
        "  (compare closure_v1_2048 vs closure_v2_2048; amp02/amp01 = can FBM drop further?)"
    );
}

/// STEP A+B — render the RECOMMENDED sculpt config (A_c=0.1 km², iters=2, K×0.5) at
/// 2048² + 8192² (where low A_c is resolvable), report upper-slope dissection +
/// floor/ridge + per-order, and the basin-area distribution for navigability (TASK 4).
#[test]
#[ignore]
fn sculpt_render() {
    use std::path::Path;
    use std::time::Instant;
    use ymir_core::erosion::stream_power::{RELIEF_V1_K, StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let seed_u = 10481999410520546993u64;
    let domain = 400.0f32;
    let base = ss.depth_scale_m as f32;
    let (state, _run) = coarse_state(seed_u);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(seed_u);
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../exports/sculpt");
    std::fs::create_dir_all(&dir).unwrap();
    eprintln!("\n=== STEP A+B — recommended sculpt (A_c=0.1 km², iters=2, K×0.5, amp 0.04) ===");
    for t in [2048usize, 8192usize] {
        let cell_km2 = (domain / t as f32).powi(2);
        let mut fc = FbmUpscaleConfig::c1_hd_production(t);
        fc.erosion = None;
        fc.bathymetry = None;
        fc.amplitude_base = 0.04; // reduced FBM (visual review)
        let t0 = Instant::now();
        let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fc).heightmap;
        let mut cfg = StreamPowerConfig::relief_v1(cell_km2, base);
        cfg.min_area_cells = 0.1 / cell_km2;
        cfg.iterations = 2;
        cfg.k = RELIEF_V1_K * 0.5;
        let sp = incise(&fbm, &cfg);
        let sp_ms = t0.elapsed().as_millis();
        let hs = hillshade(&sp, domain, base);
        hs.save_png_u8(&dir.join(format!("recommended_{t}.png"))).unwrap();
        let (w, hh) = (sp.width, sp.height);
        let (x0, y0, cw, ch) = (w / 2 - w / 8, hh / 2 - hh / 8, w / 4, hh / 4);
        let mut crop = GridF32::new(cw, ch, 0.0);
        for j in 0..ch {
            for i in 0..cw {
                crop.data[j * cw + i] = hs.data[(y0 + j) * w + (x0 + i)];
            }
        }
        crop.save_png_u8(&dir.join(format!("recommended_{t}_crop.png"))).unwrap();

        let sm: Vec<f32> = sp.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).collect();
        let peak = sm.iter().cloned().fold(f32::MIN, f32::max);
        let dr = c1_drainage(&sp, None, &C1DrainageConfig::default(), &ss);
        let a_c = cfg.min_area_cells;
        let mut chan_elev: Vec<f32> = (0..w * hh)
            .filter(|&k| dr.flow.accumulation.data[k] >= a_c && sp.data[k] > SEA)
            .map(|k| sm[k])
            .collect();
        chan_elev.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let (p90, mx) = if chan_elev.is_empty() {
            (0.0, 0.0)
        } else {
            (chan_elev[chan_elev.len() * 9 / 10], *chan_elev.last().unwrap())
        };
        let dens = chan_elev.len() as f32 * (domain / t as f32)
            / ((0..w * hh).filter(|&k| sp.data[k] > SEA).count() as f32 * cell_km2).max(1.0);
        let (tab, _) = per_order_incision(&fbm, &sp, &ss);
        eprintln!(
            "  {t}²: SP {sp_ms} ms, peak {peak:.0} m | channels reach p90 {:.0}% / max {:.0}% of peak | \
             drain.dens {dens:.2} km/km² | {}",
            p90 / peak * 100.0,
            mx / peak * 100.0,
            fmt_orders(&tab),
        );
        // TASK 4 — basin-area distribution (mouth segments = downstream None).
        let mut mouths: Vec<f32> = dr
            .rivers
            .segments
            .iter()
            .enumerate()
            .filter(|(_, s)| s.downstream.is_none())
            .map(|(i, _)| dr.segment_drainage_km2[i])
            .collect();
        mouths.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if !mouths.is_empty() {
            let (mmax, p90b, p50b) =
                (mouths[mouths.len() - 1], mouths[mouths.len() * 9 / 10], mouths[mouths.len() / 2]);
            eprintln!(
                "       basin area @mouths: max {mmax:.0} km², p90 {p90b:.0}, p50 {p50b:.0} ({} mouths)",
                mouths.len()
            );
        }
    }
    eprintln!("  renders → {}", dir.display());
}

/// TASK 1+2 — A_c × incision grid on the author's seed: does lowering A_c dissect the
/// UPPER slopes, and does bounding incision stop floors being planed to base level?
/// Metrics: drainage density (km/km²), channel-head reach (p90 channel elevation as %
/// of peak), floor/local-ridge ratio, per-order incision. 1024², domain 400 km.
#[test]
#[ignore]
fn sculpt_grid() {
    use ymir_core::erosion::stream_power::{RELIEF_V1_K, StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let seed_u = 10481999410520546993u64;
    let (t, domain) = (1024usize, 400.0f32);
    let cell_km = domain / t as f32;
    let cell_km2 = cell_km * cell_km;
    let (state, _run) = coarse_state(seed_u);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(seed_u);
    let mut fcfg = FbmUpscaleConfig::c1_hd_production(t);
    fcfg.erosion = None;
    fcfg.bathymetry = None;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fcfg).heightmap;
    let fm: Vec<f32> = fbm.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).collect();
    let peak = fm.iter().cloned().fold(f32::MIN, f32::max);
    let (w, h) = (t, t);
    eprintln!("\n=== TASK 1+2 — A_c × incision grid (seed {seed_u}, 1024², peak {peak:.0} m) ===");
    eprintln!(
        "  A_c km² | iters | K× | drain.dens km/km² | head reach %peak | floor/ridge | relief m | per-order"
    );

    let mut run = |a_c_km2: f32, iters: usize, kmult: f32| {
        let mut cfg = StreamPowerConfig::relief_v1(cell_km2, ss.depth_scale_m as f32);
        cfg.min_area_cells = a_c_km2 / cell_km2;
        cfg.k = RELIEF_V1_K * kmult;
        cfg.iterations = iters;
        let sp = incise(&fbm, &cfg);
        let sm: Vec<f32> = sp.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).collect();
        let dr = c1_drainage(&sp, None, &C1DrainageConfig::default(), &ss);
        let a_c = cfg.min_area_cells;
        // channel cells (A≥A_c land) — density + head reach.
        let chan: Vec<usize> = (0..w * h)
            .filter(|&k| dr.flow.accumulation.data[k] >= a_c && sp.data[k] > SEA)
            .collect();
        let land = (0..w * h).filter(|&k| sp.data[k] > SEA).count();
        let dens = chan.len() as f32 * cell_km / (land as f32 * cell_km2).max(1.0);
        let mut chan_elev: Vec<f32> = chan.iter().map(|&k| sm[k]).collect();
        chan_elev.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let head_p90 = if chan_elev.is_empty() { 0.0 } else { chan_elev[chan_elev.len() * 9 / 10] };
        // floor/local-ridge: for channel cells, ridge = max elev in ±10; median ratio + relief.
        let (mut ratios, mut reliefs) = (Vec::new(), Vec::new());
        for &k in chan.iter().step_by(7) {
            let (x, y) = (k % w, k / w);
            let mut ridge = sm[k];
            for dy in -10i32..=10 {
                for dx in -10i32..=10 {
                    let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                    if nx >= 0 && ny >= 0 && nx < w as i32 && ny < h as i32 {
                        ridge = ridge.max(sm[ny as usize * w + nx as usize]);
                    }
                }
            }
            if ridge > 1.0 {
                ratios.push(sm[k] / ridge);
                reliefs.push(ridge - sm[k]);
            }
        }
        let med = |v: &mut Vec<f32>| {
            if v.is_empty() {
                return 0.0;
            }
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            v[v.len() / 2]
        };
        let (tab, _) = per_order_incision(&fbm, &sp, &ss);
        eprintln!(
            "  {a_c_km2:>6.2} | {iters:>5} | {kmult:>3.1} | {dens:>16.2} | {:>15.0} | {:>11.2} | {:>7.0} | {}",
            head_p90 / peak * 100.0,
            med(&mut ratios),
            med(&mut reliefs),
            fmt_orders(&tab),
        );
    };

    eprintln!(
        "-- TASK 1: A_c sweep (iters=3, K×1.0) — does the network reach the upper slopes? --"
    );
    for a_c in [0.05f32, 0.1, 0.25, 0.5, 1.0] {
        run(a_c, 3, 1.0);
    }
    eprintln!("-- TASK 2: incision bound at A_c=0.1 — do floors stop being planed? --");
    for (it, km) in [(3usize, 1.0f32), (2, 0.5), (1, 1.0), (1, 0.5)] {
        run(0.1, it, km);
    }
    eprintln!(
        "  (want: head reach → high % of peak (upper slopes dissected); floor/ridge NOT ~0 (not planed);"
    );
    eprintln!("   drainage density plausible ~0.5–3 km/km²)");
}

/// TASK 2 — striation metric via a DIRECTIONAL POWER SPECTRUM. The old ±8 roughness
/// ratio failed because the window (~16 cells) was the size of the ~8–11 cell period
/// it tried to resolve. Here: on steep cells, extract a length-`win` (≥32) profile
/// along the contour and along the gradient, endpoint-detrend (removes the slope),
/// DFT, accumulate power per frequency. Striations = a spectral PEAK at a specific
/// wavelength in the contour direction. Returns (dominant wavelength cells,
/// anisotropic power ratio contour/gradient at that wavelength).
fn striation_spectrum(field: &GridF32, slope: &[f32], min_deg: f32, win: usize) -> (f32, f32) {
    let (w, h) = (field.width, field.height);
    let nf = win / 2;
    let mut pc = vec![0.0f64; nf]; // contour power per frequency
    let mut pg = vec![0.0f64; nf]; // gradient power
    let half = win as f32 / 2.0;
    for y in 2..h - 2 {
        for x in 2..w - 2 {
            let k = y * w + x;
            if slope[k] < min_deg || field.data[k] <= SEA {
                continue;
            }
            let gx = field.data[k + 1] - field.data[k - 1];
            let gy = field.data[k + w] - field.data[k - w];
            let gl = (gx * gx + gy * gy).sqrt().max(1e-9);
            let (ux, uy) = (gx / gl, gy / gl);
            let (cxx, cyy) = (-uy, ux);
            let mut accum = |dx: f32, dy: f32, pow: &mut [f64]| {
                let p: Vec<f32> = (0..win)
                    .map(|i| {
                        let t = i as f32 - half;
                        field.sample_bilinear_periodic(x as f32 + dx * t, y as f32 + dy * t)
                    })
                    .collect();
                // endpoint detrend (remove DC + linear slope along the line).
                let (p0, p1) = (p[0], p[win - 1]);
                let d: Vec<f32> = (0..win)
                    .map(|i| p[i] - (p0 + (p1 - p0) * i as f32 / (win as f32 - 1.0)))
                    .collect();
                for kf in 1..nf {
                    let (mut re, mut im) = (0.0f64, 0.0f64);
                    for i in 0..win {
                        let ph = -2.0 * std::f64::consts::PI * kf as f64 * i as f64 / win as f64;
                        re += d[i] as f64 * ph.cos();
                        im += d[i] as f64 * ph.sin();
                    }
                    pow[kf] += re * re + im * im;
                }
            };
            accum(ux, uy, &mut pg);
            accum(cxx, cyy, &mut pc);
        }
    }
    // Dominant frequency = peak of the contour spectrum (where striations show).
    let dom = (1..nf).max_by(|&a, &b| pc[a].partial_cmp(&pc[b]).unwrap()).unwrap_or(1);
    let wavelength = win as f32 / dom as f32;
    let ratio = (pc[dom] / pg[dom].max(1e-12)) as f32;
    (wavelength, ratio)
}

/// TASK 2 validation — does the spectrum metric MOVE when max_anisotropy 3→1 (and
/// amplitude 0.16→0.04)? If yes the metric works + the knob is a real lever; if no,
/// the knobs genuinely are not the striation cause. 1024², relief-v1 incision.
#[test]
#[ignore]
fn striation_spectrum_validate() {
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let (t, domain) = (1024usize, 400.0f32);
    let cell_km2 = (domain / t as f32).powi(2);
    let base = ss.depth_scale_m as f32;
    let (state, _run) = coarse_state(SEED);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(SEED);
    let sp = StreamPowerConfig::relief_v1(cell_km2, ss.depth_scale_m as f32);
    eprintln!("\n=== TASK 2 — striation power-spectrum validation (win=48, steep>20°) ===");
    let run = |label: &str, aniso: f64, amp: f64| {
        let mut fc = FbmUpscaleConfig::c1_hd_production(t);
        fc.erosion = None;
        fc.bathymetry = None;
        fc.max_anisotropy = aniso;
        fc.amplitude_base = amp;
        let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fc).heightmap;
        // pre-incision (raw FBM) + post-incision.
        let sl_f = slope_deg_field(&fbm, domain, base);
        let (wlf, rf) = striation_spectrum(&fbm, &sl_f, 20.0, 48);
        let field = incise(&fbm, &sp);
        let sl = slope_deg_field(&field, domain, base);
        let (wl, r) = striation_spectrum(&field, &sl, 20.0, 48);
        eprintln!(
            "  {label:<28}: pre-FBM λ {wlf:.1} ratio {rf:.2} | post-incision λ {wl:.1} ratio {r:.2}"
        );
    };
    run("baseline aniso=3 amp=0.16", 3.0, 0.16);
    run("aniso=1 (isotropic) amp=0.16", 1.0, 0.16);
    run("aniso=3 amp=0.04", 3.0, 0.04);
    run("aniso=1 amp=0.04", 1.0, 0.04);
    eprintln!(
        "  (ratio ≫1 = striations along contour; if it drops 3→1 the anisotropy knob is the lever)"
    );
}

/// B1 — striation metric. On STEEP cells (slope > `min_deg`), sample the altitude
/// profile ALONG the gradient and ALONG the contour, and measure short-wavelength
/// roughness (RMS of the 1-D Laplacian) in each. Anisotropic FBM elongates noise
/// along the slope → filaments run downslope, wiggling ACROSS (contour) → contour
/// roughness ≫ gradient roughness = striations. Returns (rough_grad, rough_contour,
/// asymmetry contour/grad, dominant wavelength in cells along the contour).
fn striation_metric(field: &GridF32, slope: &[f32], min_deg: f32) -> (f32, f32, f32, f32) {
    let (w, h) = (field.width, field.height);
    let r = 8i32;
    let (mut sg, mut sc, mut n, mut wl_acc, mut wl_n) = (0.0f64, 0.0f64, 0u64, 0.0f64, 0u64);
    for y in 2..h - 2 {
        for x in 2..w - 2 {
            let k = y * w + x;
            if slope[k] < min_deg || field.data[k] <= SEA {
                continue;
            }
            let gx = field.data[k + 1] - field.data[k - 1];
            let gy = field.data[k + w] - field.data[k - w];
            let gl = (gx * gx + gy * gy).sqrt().max(1e-9);
            let (ux, uy) = (gx / gl, gy / gl); // gradient unit
            let (cx, cy) = (-uy, ux); // contour unit
            let prof = |dx: f32, dy: f32| -> Vec<f32> {
                (-r..=r)
                    .map(|i| {
                        field.sample_bilinear_periodic(
                            x as f32 + dx * i as f32,
                            y as f32 + dy * i as f32,
                        )
                    })
                    .collect()
            };
            let rough = |p: &[f32]| -> f32 {
                let mut s = 0.0;
                for i in 1..p.len() - 1 {
                    let l = p[i - 1] - 2.0 * p[i] + p[i + 1];
                    s += (l * l) as f64;
                }
                (s / (p.len() - 2) as f64).sqrt() as f32
            };
            let pg = prof(ux, uy);
            let pc = prof(cx, cy);
            sg += rough(&pg) as f64;
            sc += rough(&pc) as f64;
            n += 1;
            // wavelength along contour: mean-crossings of the detrended profile.
            let mean = pc.iter().sum::<f32>() / pc.len() as f32;
            let cross =
                (1..pc.len()).filter(|&i| (pc[i] - mean) * (pc[i - 1] - mean) < 0.0).count();
            if cross > 0 {
                wl_acc += (2.0 * r as f64) / cross as f64;
                wl_n += 1;
            }
        }
    }
    if n == 0 {
        return (0.0, 0.0, 0.0, 0.0);
    }
    let (rg, rc) = (sg / n as f64, sc / n as f64);
    ((rg) as f32, (rc) as f32, (rc / rg.max(1e-12)) as f32, (wl_acc / wl_n.max(1) as f64) as f32)
}

/// PART A — regression guard for the `ref/relief-streampower-v1` config. NOT a
/// tuning sweep: it asserts the key legibility metrics stay in collapse-catching
/// ranges so a later change cannot silently degrade the relief (droplet-style
/// collapse, lost valleys, over-smoothing). Seed 42, 1024², domain 400 km,
/// uncoupled. #[ignore] (heavy); run before/after any erosion/upscale change.
#[test]
#[ignore]
fn relief_v1_regression() {
    use ymir_core::erosion::stream_power::{RELIEF_V1_A_C_KM2, StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    // 2048² (A_c=0.1 km² is sub-cell at 1024²) + amp 0.04 = the recommended config.
    let (t, domain) = (2048usize, 400.0f32);
    let cell_km2 = (domain / t as f32).powi(2);
    let base = ss.depth_scale_m as f32;
    let a_c = RELIEF_V1_A_C_KM2 / cell_km2;
    let (state, _run) = coarse_state(SEED);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(SEED);
    let mut fcfg = FbmUpscaleConfig::c1_hd_production(t);
    fcfg.erosion = None;
    fcfg.bathymetry = None;
    fcfg.amplitude_base = 0.04;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fcfg).heightmap;
    let sp = incise(&fbm, &StreamPowerConfig::relief_v1(cell_km2, ss.depth_scale_m as f32));

    let relief = drainage_relief_m(&sp, &ss);
    let slope = slope_deg_field(&sp, domain, base);
    let corr = channel_corridor(&sp, &ss, a_c);
    let (vf5_km, _, _, _) = valley_floor(&slope, &corr, &sp, cell_km2);
    let (_, l1, _) = flank_contiguity(&slope, t, t);
    let (tab, _) = per_order_incision(&fbm, &sp, &ss);
    let s2 = tab.iter().find(|(o, _)| *o == 2).map(|(_, m)| *m).unwrap_or(0.0);
    let s4 = tab.iter().find(|(o, _)| *o == 4).map(|(_, m)| *m).unwrap_or(0.0);
    eprintln!(
        "relief-v1 @2048² (A_c=0.1, iters=2, K=1500, amp 0.04): drainage relief {relief:.0} m, \
         valley floor<5° {vf5_km:.0} km², largest flank {l1} cells, per-order {}",
        fmt_orders(&tab),
    );
    // REBASED to the sculpt reference (A_c=0.1 km², iters=2, K=1500, amp 0.04, seed 42,
    // 2048²): relief 397 m, valley floor<5° 5142 km², flank 113909 cells, S2 344, S4
    // 162. Ranges are tight (±~20 %) around these measured values, NOT loosened — they
    // catch silent degradation of the shipped sculpt config.
    assert!(
        (320.0..=480.0).contains(&relief),
        "drainage relief {relief} m out of [320,480] (ref 397)"
    );
    assert!((280.0..=420.0).contains(&s2), "S2 incision {s2} m out of [280,420] (ref 344)");
    assert!((120.0..=210.0).contains(&s4), "S4 incision {s4} m out of [120,210] (ref 162)");
    assert!(
        (4000.0..=6500.0).contains(&vf5_km),
        "valley floor<5° {vf5_km} km² out of [4000,6500] (ref 5142)"
    );
    assert!(l1 > 60000, "largest steep flank {l1} cells (<60000 ⇒ degraded; ref 113909)");
}

/// Per-cell surface slope in DEGREES (central differences; `depth_scale` is the
/// effective vertical scale, i.e. base × coupling factor; `domain_km` sets km/cell).
fn slope_deg_field(field: &GridF32, domain_km: f32, depth_scale: f32) -> Vec<f32> {
    let (w, h) = (field.width, field.height);
    let cell_m = domain_km / w as f32 * 1000.0;
    let norm_to_m = 2.0 * 1.13 * depth_scale;
    let mut out = vec![0.0f32; w * h];
    for y in 0..h {
        for x in 0..w {
            let l = field.data[y * w + (x + w - 1) % w];
            let r = field.data[y * w + (x + 1) % w];
            let d = field.data[((y + h - 1) % h) * w + x];
            let u = field.data[((y + 1) % h) * w + x];
            let sx = (r - l) * 0.5 * norm_to_m / cell_m;
            let sy = (u - d) * 0.5 * norm_to_m / cell_m;
            out[y * w + x] = (sx * sx + sy * sy).sqrt().atan().to_degrees();
        }
    }
    out
}

/// Channel corridor mask: cells with accumulation ≥ A_c, dilated by one ring.
fn channel_corridor(field: &GridF32, ss: &SteinSteinParams, a_c: f32) -> Vec<bool> {
    let dr = c1_drainage(field, None, &C1DrainageConfig::default(), ss);
    let (w, h) = (field.width, field.height);
    let chan: Vec<bool> =
        (0..w * h).map(|k| dr.flow.accumulation.data[k] >= a_c && field.data[k] > SEA).collect();
    let mut corr = chan.clone();
    for y in 0..h {
        for x in 0..w {
            if !chan[y * w + x] {
                continue;
            }
            for (dx, dy) in
                [(-1i32, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)]
            {
                let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                if nx >= 0 && ny >= 0 && nx < w as i32 && ny < h as i32 {
                    corr[ny as usize * w + nx as usize] = true;
                }
            }
        }
    }
    corr
}

/// (a) VALLEY-FLOOR: contiguous flat ground in the channel corridor. Returns
/// (km² <5°, hexes <5°, km² <10°, hexes <10°).
fn valley_floor(
    slope: &[f32],
    corridor: &[bool],
    field: &GridF32,
    cell_km2: f32,
) -> (f32, usize, f32, usize) {
    let (mut c5, mut c10) = (0usize, 0usize);
    for k in 0..slope.len() {
        if corridor[k] && field.data[k] > SEA {
            if slope[k] < 5.0 {
                c5 += 1;
            }
            if slope[k] < 10.0 {
                c10 += 1;
            }
        }
    }
    (c5 as f32 * cell_km2, c5, c10 as f32 * cell_km2, c10)
}

/// (c) FLANK CONTIGUITY: connected components (4-conn, non-periodic) of the >30°
/// mask. Returns (num components ≥ 4 cells, largest, 2nd largest) in cells.
fn flank_contiguity(slope: &[f32], w: usize, h: usize) -> (usize, usize, usize) {
    let steep: Vec<bool> = slope.iter().map(|&s| s > 30.0).collect();
    let mut seen = vec![false; w * h];
    let mut sizes = Vec::new();
    for start in 0..w * h {
        if !steep[start] || seen[start] {
            continue;
        }
        let mut sz = 0usize;
        let mut stack = vec![start];
        seen[start] = true;
        while let Some(k) = stack.pop() {
            sz += 1;
            let (x, y) = (k % w, k / w);
            for (dx, dy) in [(-1i32, 0), (1, 0), (0, -1), (0, 1)] {
                let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                if nx >= 0 && ny >= 0 && nx < w as i32 && ny < h as i32 {
                    let nk = ny as usize * w + nx as usize;
                    if steep[nk] && !seen[nk] {
                        seen[nk] = true;
                        stack.push(nk);
                    }
                }
            }
        }
        sizes.push(sz);
    }
    sizes.sort_unstable_by(|a, b| b.cmp(a));
    let n = sizes.iter().filter(|&&s| s >= 4).count();
    (n, sizes.first().copied().unwrap_or(0), sizes.get(1).copied().unwrap_or(0))
}

/// (b) cliff transition + valley width/depth from perpendicular profiles at trunk
/// probes. Returns (median transition <10°→>30° in m, median valley width m, median
/// depth m). `cell_m` = metres/cell; slope in degrees supplied per cell.
fn profile_metrics(
    field: &GridF32,
    slope: &[f32],
    depth_scale: f32,
    probes: &[((usize, usize), (f32, f32))],
    cell_m: f32,
) -> (f32, f32, f32) {
    let (w, h) = (field.width, field.height);
    // Metres via the effective (coupled) depth scale: (norm − 0.5)·2·1.13·depth.
    let fm: Vec<f32> = field.data.iter().map(|&n| (n - 0.5) * 2.0 * 1.13 * depth_scale).collect();
    let (mut trans, mut widths, mut depths) = (Vec::new(), Vec::new(), Vec::new());
    let r = 40i32; // ±40 cells — enough at fine resolution to reach the flanks
    for &((cx, cy), (px, py)) in probes {
        let sample = |o: i32| -> (f32, f32) {
            let sx = (cx as f32 + px * o as f32).round().clamp(0.0, w as f32 - 1.0) as usize;
            let sy = (cy as f32 + py * o as f32).round().clamp(0.0, h as f32 - 1.0) as usize;
            (fm[sy * w + sx], slope[sy * w + sx])
        };
        // Find the true channel bottom (min altitude) over the profile, then walk
        // OUTWARD to each side to the first >30° flank; the floor edge is the last
        // <10° cell before it. Robust at fine resolution (o=0 may be on a wall).
        let mut bo = 0i32;
        let mut bmin = f32::MAX;
        for o in -r..=r {
            let (a, _) = sample(o);
            if a < bmin {
                bmin = a;
                bo = o;
            }
        }
        let side = |dir: i32| -> Option<(i32, f32, i32)> {
            let (mut floor_edge, mut o) = (bo, bo);
            loop {
                o += dir;
                if o.abs() > r {
                    return None;
                }
                let (a, sl) = sample(o);
                if sl < 10.0 {
                    floor_edge = o;
                }
                if sl > 30.0 {
                    return Some((o, a, floor_edge));
                }
            }
        };
        if let (Some((fl_o, fl_a, fe_l)), Some((fr_o, fr_a, fe_r))) = (side(-1), side(1)) {
            trans.push((((fl_o - fe_l).abs()).min((fr_o - fe_r).abs())) as f32 * cell_m);
            widths.push((fr_o - fl_o).max(1) as f32 * cell_m);
            depths.push(fl_a.min(fr_a) - bmin);
        }
    }
    let med = |v: &mut Vec<f32>| {
        if v.is_empty() {
            return 0.0;
        }
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v[v.len() / 2]
    };
    (med(&mut trans), med(&mut widths), med(&mut depths))
}

/// Median incision (m) per Strahler order + whether the ordering is MONOTONE
/// (S1 < S2 < ... — trunks carve more than headwaters, the physical target).
fn per_order_incision(
    fbm: &GridF32,
    out: &GridF32,
    ss: &SteinSteinParams,
) -> (Vec<(u8, f32)>, bool) {
    let fm = to_metres(fbm, ss);
    let om = to_metres(out, ss);
    let dr = c1_drainage(out, None, &C1DrainageConfig::default(), ss);
    let mut per: std::collections::HashMap<u8, Vec<f32>> = std::collections::HashMap::new();
    for s in &dr.rivers.segments {
        for &(x, y) in &s.points {
            let k = y as usize * out.width + x as usize;
            per.entry(s.strahler_order).or_default().push(fm.data[k] - om.data[k]);
        }
    }
    let mut orders: Vec<u8> = per.keys().copied().collect();
    orders.sort();
    let mut table = Vec::new();
    for o in &orders {
        let v = per.get_mut(o).unwrap();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        table.push((*o, v[v.len() / 2]));
    }
    // Monotone up to the highest 4 orders (ignore sparse top order).
    let mono = table.windows(2).take(3).all(|w| w[0].1 <= w[1].1 + 5.0);
    (table, mono)
}

fn fmt_orders(t: &[(u8, f32)]) -> String {
    t.iter().map(|(o, m)| format!("S{o}={m:.0}m")).collect::<Vec<_>>().join(" ")
}

/// Per-Strahler-order valley WIDTH / DEPTH (m). At every channel-segment point, walk
/// perpendicular to flow counting the contiguous valley-floor cells (slope < 8°) →
/// width; depth = local ridge (max within ±8 perpendicular) − channel. Median per
/// order. The slit pathology = W/D ≪ 1 on trunks; a valley = W/D grows with order.
fn per_order_width_depth(
    field: &GridF32,
    domain_km: f32,
    depth_scale: f32,
    ss: &SteinSteinParams,
) -> Vec<(u8, f32, f32, f32)> {
    let (w, h) = (field.width, field.height);
    let cell_m = domain_km / w as f32 * 1000.0;
    let norm_to_m = 2.0 * 1.13 * depth_scale;
    let slope = slope_deg_field(field, domain_km, depth_scale);
    let dr = c1_drainage(field, None, &C1DrainageConfig::default(), ss);
    let mut wper: std::collections::HashMap<u8, Vec<f32>> = std::collections::HashMap::new();
    let mut dper: std::collections::HashMap<u8, Vec<f32>> = std::collections::HashMap::new();
    for s in &dr.rivers.segments {
        for win in s.points.windows(2) {
            let (x, y) = (win[0].0 as i32, win[0].1 as i32);
            // flow direction from this point to the next; perpendicular = (-dy,dx).
            let (fx, fy) = (win[1].0 as i32 - x, win[1].1 as i32 - y);
            if fx == 0 && fy == 0 {
                continue;
            }
            let (px, py) = (-fy.signum(), fx.signum());
            let zc = field.data[y as usize * w + x as usize];
            // width: contiguous <8° floor cells each way.
            let mut width = 1i32;
            let mut ridge = zc;
            for dir in [-1i32, 1] {
                let mut o = 0i32;
                loop {
                    o += dir;
                    let (nx, ny) = (x + px * o, y + py * o);
                    if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 || o.abs() > 8 {
                        break;
                    }
                    let nk = ny as usize * w + nx as usize;
                    ridge = ridge.max(field.data[nk]);
                    if slope[nk] < 8.0 && field.data[nk] > SEA {
                        width += 1;
                    } else {
                        break;
                    }
                }
            }
            wper.entry(s.strahler_order).or_default().push(width as f32 * cell_m);
            dper.entry(s.strahler_order).or_default().push((ridge - zc) * norm_to_m);
        }
    }
    let med = |v: &mut Vec<f32>| {
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v[v.len() / 2]
    };
    let mut orders: Vec<u8> = wper.keys().copied().collect();
    orders.sort();
    orders
        .iter()
        .filter_map(|o| {
            let (wl, dl) = (wper.get_mut(o)?, dper.get_mut(o)?);
            if wl.is_empty() {
                return None;
            }
            let (wm, dm) = (med(wl), med(dl));
            Some((*o, wm, dm, if dm > 1.0 { wm / dm } else { 0.0 }))
        })
        .collect()
}

/// CLOSURE STEP 3 — shape metrics per configuration: v1 (no closures), +critical
/// slope (a), +lateral (b), +both. Reports the max slope + steep shares (does S_c
/// bound it?), floor/local-ridge, drainage relief + Strahler histogram (health),
/// per-order W/D (slit → valley), the striation ratio ON UPPER SLOPES (has emergent
/// structure replaced the FBM pattern?), and a crest-curvature proxy (arêtes). 2048²,
/// amp 0.04, seed 10481999410520546993. Read-only.
#[test]
#[ignore]
fn closure_grid() {
    use ymir_core::erosion::stream_power::{
        RELIEF_V2_CRITICAL_SLOPE, RELIEF_V2_LATERAL, StreamPowerConfig, incise,
    };
    let ss = SteinSteinParams::default();
    let seed_u = 10481999410520546993u64;
    let (t, domain) = (2048usize, 400.0f32);
    let base = ss.depth_scale_m as f32;
    let cell_km2 = (domain / t as f32).powi(2);
    let (state, _run) = coarse_state(seed_u);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(seed_u);
    let mut fc = FbmUpscaleConfig::c1_hd_production(t);
    fc.erosion = None;
    fc.bathymetry = None;
    fc.amplitude_base = 0.04;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fc).heightmap;
    let fm: Vec<f32> = fbm.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).collect();
    let peak = fm.iter().cloned().fold(f32::MIN, f32::max);

    eprintln!(
        "\n=== CLOSURE STEP 3 — shape metrics (seed {seed_u}, {t}², amp 0.04, peak {peak:.0} m) ==="
    );
    let configs: [(&str, f32, f32); 4] = [
        ("v1 (none)", 0.0, 0.0),
        ("+crit slope (a)", RELIEF_V2_CRITICAL_SLOPE, 0.0),
        ("+lateral (b)", 0.0, RELIEF_V2_LATERAL),
        ("+both (v2)", RELIEF_V2_CRITICAL_SLOPE, RELIEF_V2_LATERAL),
    ];
    for (label, sc, lat) in configs {
        let mut cfg = StreamPowerConfig::relief_v1(cell_km2, base);
        cfg.critical_slope = sc;
        cfg.lateral_erosion = lat;
        if sc > 0.0 {
            cfg.diffusion = ymir_core::erosion::stream_power::RELIEF_V2_DIFFUSION;
        }
        let sp = incise(&fbm, &cfg);
        let sm: Vec<f32> = sp.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).collect();
        let slope = slope_deg_field(&sp, domain, base);
        // slope distribution on land.
        let land: Vec<usize> = (0..sp.data.len()).filter(|&k| sp.data[k] > SEA).collect();
        let nland = land.len().max(1) as f32;
        let frac =
            |thr: f32| land.iter().filter(|&&k| slope[k] > thr).count() as f32 / nland * 100.0;
        let mx = land.iter().map(|&k| slope[k]).fold(0.0f32, f32::max);
        // floor/local-ridge on channel cells.
        let a_c = cfg.min_area_cells;
        let dr = c1_drainage(&sp, None, &C1DrainageConfig::default(), &ss);
        let chan: Vec<usize> = (0..sp.data.len())
            .filter(|&k| dr.flow.accumulation.data[k] >= a_c && sp.data[k] > SEA)
            .collect();
        let mut ratios = Vec::new();
        for &k in chan.iter().step_by(7) {
            let (x, y) = (k % t, k / t);
            let mut ridge = sm[k];
            for dy in -10i32..=10 {
                for dx in -10i32..=10 {
                    let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                    if nx >= 0 && ny >= 0 && nx < t as i32 && ny < t as i32 {
                        ridge = ridge.max(sm[ny as usize * t + nx as usize]);
                    }
                }
            }
            if ridge > 1.0 {
                ratios.push(sm[k] / ridge);
            }
        }
        ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let floor_ridge = if ratios.is_empty() { 0.0 } else { ratios[ratios.len() / 2] };
        // striation ratio on UPPER SLOPES: mask slope to cells above the 60th land
        // elevation pct, then the directional spectrum (contour/gradient power ratio).
        let mut land_elev: Vec<f32> = land.iter().map(|&k| sm[k]).collect();
        land_elev.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let e60 = land_elev[land_elev.len() * 6 / 10];
        let upper_slope: Vec<f32> =
            (0..sp.data.len()).map(|k| if sm[k] >= e60 { slope[k] } else { 0.0 }).collect();
        let (wl_cells, aniso) = striation_spectrum(&sp, &upper_slope, 20.0, 48);
        // crest curvature (m) on steep upper cells: |z_l+z_r−2z| + |z_u+z_d−2z|.
        let norm_to_m = 2.0 * 1.13 * base;
        let mut curv = Vec::new();
        for &k in &land {
            if sm[k] >= e60
                && slope[k] > 30.0
                && k % t > 0
                && k % t < t - 1
                && k >= t
                && k < t * t - t
            {
                let lap = (sp.data[k - 1] + sp.data[k + 1] - 2.0 * sp.data[k]).abs()
                    + (sp.data[k - t] + sp.data[k + t] - 2.0 * sp.data[k]).abs();
                curv.push(lap * norm_to_m);
            }
        }
        curv.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let crest = if curv.is_empty() { 0.0 } else { curv[curv.len() / 2] };
        let relief = drainage_relief_m(&sp, &ss);
        // Strahler histogram (segment count per order) + W/D.
        let mut hist: std::collections::HashMap<u8, usize> = std::collections::HashMap::new();
        for s in &dr.rivers.segments {
            *hist.entry(s.strahler_order).or_default() += 1;
        }
        let mut ho: Vec<u8> = hist.keys().copied().collect();
        ho.sort();
        let hstr = ho.iter().map(|o| format!("S{o}:{}", hist[o])).collect::<Vec<_>>().join(" ");
        let wd = per_order_width_depth(&sp, domain, base, &ss);
        let wdstr = wd
            .iter()
            .map(|(o, w, d, r)| format!("S{o} {w:.0}/{d:.0}={r:.2}"))
            .collect::<Vec<_>>()
            .join("  ");
        eprintln!("\n-- {label} --");
        eprintln!(
            "  slope: >30 {:.1}%  >35 {:.1}%  >40 {:.1}%  MAX {mx:.1}°  | floor/ridge {floor_ridge:.2}  relief {relief:.0} m",
            frac(30.0),
            frac(35.0),
            frac(40.0),
        );
        eprintln!(
            "  upper-slope striation: λ {wl_cells:.1} cells, contour/gradient power {aniso:.2}  | crest curv {crest:.0} m"
        );
        eprintln!("  Strahler {hstr}");
        eprintln!("  W/D per order (m): {wdstr}");
    }
    eprintln!(
        "\n  (want: MAX slope → ~33° with (a); W/D grows with order & >1 on trunks with (b);"
    );
    eprintln!(
        "   striation contour/gradient power → ~1 (isotropic, FBM pattern gone); crest curv up = arêtes)"
    );
}

/// TASK 1+2 — hillslope regime: sweep the incision threshold θ and the critical
/// drainage area A_c, report the per-order table for each, and where the ordering
/// becomes monotone (trunks > headwaters). 1024², seed 42.
#[test]
#[ignore]
fn stream_power_hillslope() {
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    use ymir_core::tectonics_c1::production_upscale::c1_cell_area_km2;
    let ss = SteinSteinParams::default();
    let t = 1024usize;
    let (state, _run) = coarse_state(SEED);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(SEED);
    let mut fcfg = FbmUpscaleConfig::c1_hd_production(t);
    fcfg.erosion = None;
    fcfg.bathymetry = None;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fcfg).heightmap;
    let km2 = c1_cell_area_km2(t);
    let base = |th: f32, ac: f32| StreamPowerConfig {
        k: 1.0,
        m: 0.5,
        n: 1.0,
        dt: 1.0,
        iterations: 3,
        sea_level: SEA,
        diffusion: 0.0,
        diffusion_substeps: 4,
        threshold: th,
        min_area_cells: ac,
        cell_km: 1.0,
        depth_scale_m: 5000.0,
        ..Default::default()
    };

    eprintln!("\n=== TASK 1 — incision threshold θ sweep (K=1,m=.5,n=1,iters=3) ===");
    for th in [0.0f32, 0.02, 0.05, 0.1, 0.15, 0.2] {
        let (tab, mono) = per_order_incision(&fbm, &incise(&fbm, &base(th, 0.0)), &ss);
        eprintln!(
            "  θ={th:>5.1}: {}  {}",
            fmt_orders(&tab),
            if mono { "MONOTONE ✓" } else { "inverted" }
        );
    }
    eprintln!("\n=== TASK 2 — critical drainage area A_c sweep ===");
    for ac in [0.0f32, 5.0, 25.0, 100.0, 400.0] {
        let (tab, mono) = per_order_incision(&fbm, &incise(&fbm, &base(0.0, ac)), &ss);
        eprintln!(
            "  A_c={ac:>5.0} cells ({:>5.1} km²): {}  {}",
            ac * km2,
            fmt_orders(&tab),
            if mono { "MONOTONE ✓" } else { "inverted" }
        );
    }
    eprintln!("  (target: S1 → tens of m, ordering monotone S1<S2<S3<S4)");
}

/// TASK 4 — FROZEN config at production scale (8192², author seed, domain 400 km,
/// UNCOUPLED). One drainage call to keep 67 M cells feasible. Reports the legibility
/// metrics + per-order + runtime + peak RSS + K resolution-stability.
#[test]
#[ignore]
fn stream_power_confirm8192() {
    use std::time::Instant;
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let seed_u = 10481999410520546993u64;
    let t = 8192usize;
    let domain = 400.0f32;
    let cell_m = domain / t as f32 * 1000.0; // 48.8 m/cell
    let km2 = (domain / t as f32).powi(2);
    let base = ss.depth_scale_m as f32;
    // A_c in PHYSICAL km² (resolution-stable): 7.6 km² = the 50-cell value at 1024²
    // that fixed headwaters. In CELLS it must scale with resolution, else the
    // channel head shrinks and headwaters over-carve again at 8192².
    let a_c_km2 = 7.6f32;
    let a_c = a_c_km2 / km2; // cells at THIS resolution
    let (state, _run) = coarse_state(seed_u);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(seed_u);
    let mut fcfg = FbmUpscaleConfig::c1_hd_production(t);
    fcfg.erosion = None;
    fcfg.bathymetry = None;
    let t0 = Instant::now();
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fcfg).heightmap;
    let fbm_ms = t0.elapsed().as_millis();
    let cfg = StreamPowerConfig {
        k: 3.0,
        m: 0.5,
        n: 1.0,
        dt: 1.0,
        iterations: 3,
        sea_level: SEA,
        diffusion: 0.05,
        diffusion_substeps: 4,
        min_area_cells: a_c,
        threshold: 0.0,
        cell_km: 1.0,
        depth_scale_m: 5000.0,
        ..Default::default()
    };
    eprintln!("  A_c = {a_c_km2} km² = {a_c:.0} cells at {t}²");
    let t0 = Instant::now();
    let sp = incise(&fbm, &cfg);
    let sp_ms = t0.elapsed().as_millis();

    // Single drainage on the incised field — reused for corridor, probes, per-order.
    let dr = c1_drainage(&sp, None, &C1DrainageConfig::default(), &ss);
    let (w, h) = (t, t);
    // slope (uncoupled: depth = base).
    let slope = slope_deg_field(&sp, domain, base);
    // corridor (A≥A_c, dilate 1).
    let chan: Vec<bool> =
        (0..w * h).map(|k| dr.flow.accumulation.data[k] >= a_c && sp.data[k] > SEA).collect();
    let mut corr = chan.clone();
    for y in 1..h - 1 {
        for x in 1..w - 1 {
            if chan[y * w + x] {
                for (dx, dy) in [(-1i32, 0), (1, 0), (0, -1), (0, 1)] {
                    corr[((y as i32 + dy) as usize) * w + (x as i32 + dx) as usize] = true;
                }
            }
        }
    }
    let (v5k, v5h, v10k, v10h) = valley_floor(&slope, &corr, &sp, km2);
    let (nc, l1, l2) = flank_contiguity(&slope, w, h);
    // probes: top-Strahler segments.
    let mut segs: Vec<_> = dr.rivers.segments.iter().filter(|s| s.points.len() >= 5).collect();
    segs.sort_by_key(|s| std::cmp::Reverse(s.strahler_order));
    let probes: Vec<((usize, usize), (f32, f32))> = segs
        .iter()
        .take(16)
        .map(|s| {
            let mid = s.points.len() / 2;
            let (ax, ay) = s.points[mid - 1];
            let (bx, by) = s.points[mid + 1];
            let (tx, ty) = (bx as f32 - ax as f32, by as f32 - ay as f32);
            let tl = (tx * tx + ty * ty).sqrt().max(1e-6);
            ((s.points[mid].0 as usize, s.points[mid].1 as usize), (-ty / tl, tx / tl))
        })
        .collect();
    let (tr, wd, dp) = profile_metrics(&sp, &slope, base, &probes, cell_m);
    // per-order incision.
    let fm: Vec<f32> = fbm.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).collect();
    let sm: Vec<f32> = sp.data.iter().map(|&n| c1_altitude_norm_to_metres(n, &ss)).collect();
    let mut per: std::collections::HashMap<u8, Vec<f32>> = std::collections::HashMap::new();
    for s in &dr.rivers.segments {
        for &(x, y) in &s.points {
            let k = y as usize * w + x as usize;
            per.entry(s.strahler_order).or_default().push(fm[k] - sm[k]);
        }
    }
    eprintln!(
        "\n=== TASK 4 — FROZEN config @ {t}² (author seed, domain {domain} km, UNCOUPLED) ==="
    );
    eprintln!("  cfg: K=3 m=0.5 n=1 iters=3 A_c=50 D=0.05 θ=0 droplets=off, cf=1.0 (uncoupled)");
    eprintln!("  runtime: FBM {fbm_ms} ms, stream-power {sp_ms} ms; peak RSS {} MB", peak_ws_mb());
    eprintln!(
        "  (a) valley floor: <5° {v5k:.0} km² ({v5h} hex) | <10° {v10k:.0} km² ({v10h} hex) @ {cell_m:.0} m/hex"
    );
    eprintln!(
        "  (b) cliff transition <10°→>30°: {tr:.0} m; trunk valley width {wd:.0} m, depth {dp:.0} m (W/D {:.1})",
        if dp > 1.0 { wd / dp } else { 0.0 }
    );
    eprintln!("  (c) flanks (>30°): {nc} components, largest {l1}/{l2} cells");
    let mut orders: Vec<u8> = per.keys().copied().collect();
    orders.sort();
    eprint!("  per-order incision: ");
    for o in orders {
        let v = per.get_mut(&o).unwrap();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        eprint!("S{o}={:.0}m ", v[v.len() / 2]);
    }
    eprintln!("\n  (vs 1024²: S1=25 S2=414 S3=307 S4=136 — resolution-stable ⇒ K holds)");
}

/// TASK 1/2/3 (criterion change) — judge LEGIBILITY: (a) valley-floor area,
/// (b) cliff transition sharpness, (c) flank contiguity, + valley width/depth, over
/// a K sweep (at full coupling) and a coupling sweep (at a chosen K). 1024², 400 km.
#[test]
#[ignore]
fn stream_power_legible() {
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let t = 1024usize;
    let domain = 400.0f32;
    let cell_m = domain / t as f32 * 1000.0;
    let km2 = (domain / t as f32).powi(2);
    let base = ss.depth_scale_m as f32;
    let a_c = 50.0f32;
    let (state, _run) = coarse_state(SEED);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(SEED);
    let mut fcfg = FbmUpscaleConfig::c1_hd_production(t);
    fcfg.erosion = None;
    fcfg.bathymetry = None;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fcfg).heightmap;
    let cfg = |k: f32| StreamPowerConfig {
        k,
        m: 0.5,
        n: 1.0,
        dt: 1.0,
        iterations: 3,
        sea_level: SEA,
        diffusion: 0.05,
        diffusion_substeps: 4,
        min_area_cells: a_c,
        threshold: 0.0,
        cell_km: 1.0,
        depth_scale_m: 5000.0,
        ..Default::default()
    };
    let report = |label: String, field: &GridF32, cf: f32| {
        let depth = base * cf;
        let slope = slope_deg_field(field, domain, depth);
        let corr = channel_corridor(field, &ss, a_c);
        let (v5k, v5h, v10k, v10h) = valley_floor(&slope, &corr, field, km2);
        let (nc, l1, l2) = flank_contiguity(&slope, t, t);
        let probes = channel_probes(field, &ss, 12);
        let (tr, wd, dp) = profile_metrics(field, &slope, depth, &probes, cell_m);
        // near-channel vs interfluve mean slope.
        let (mut sc, mut nc_cells, mut si, mut ni) = (0.0f64, 0u64, 0.0f64, 0u64);
        for k in 0..t * t {
            if field.data[k] <= SEA {
                continue;
            }
            if corr[k] {
                sc += slope[k] as f64;
                nc_cells += 1;
            } else {
                si += slope[k] as f64;
                ni += 1;
            }
        }
        let (tab, _) = per_order_incision(&fbm, field, &ss);
        eprintln!(
            "  {label}: vfloor <5° {v5k:.0}km²({v5h}h) <10° {v10k:.0}km²({v10h}h) | flanks {nc} comp (top {l1}/{l2}) | \
             cliff transition {tr:.0}m width {wd:.0}m depth {dp:.0}m (W/D {:.1}) | slope chan {:.0}° inter {:.0}° | {}",
            if dp > 1.0 { wd / dp } else { 0.0 },
            sc / nc_cells.max(1) as f64,
            si / ni.max(1) as f64,
            fmt_orders(&tab),
        );
    };
    eprintln!("\n=== TASK 2 — K sweep (A_c=50, D=0.05), FULL coupling cf=0.39, domain 400 km ===");
    for k in [1.0f32, 2.0, 3.0, 5.0] {
        let f = incise(&fbm, &cfg(k));
        report(format!("K={k:.0} cf=.39"), &f, 0.39);
    }
    eprintln!("\n=== TASK 3 — coupling sweep at K=3 (want steep FLANKS, floors buildable) ===");
    let f3 = incise(&fbm, &cfg(3.0));
    for cf in [0.39f32, 0.55, 0.70, 1.0] {
        report(format!("K=3 cf={cf:.2}"), &f3, cf);
    }
    eprintln!(
        "  (recommend the cf that MAXIMISES flank steepness+sharpness while keeping enough contiguous valley floor;"
    );
    eprintln!(
        "   reject if interfluves go uniformly steep = steep no longer concentrated near channels)"
    );
}

/// TASK 4+5 — recommended regime-split config: coupled-vs-uncoupled per-order +
/// slopes, and the weak-droplet coupling (does a light texture pass preserve the SP
/// valleys, unlike the full pass 323→24 m?). 1024², seed 42.
#[test]
#[ignore]
fn stream_power_recommended() {
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    use ymir_core::tectonics_c1::land_topology::slope_shares;
    let ss = SteinSteinParams::default();
    let t = 1024usize;
    let cells = (t * t) as f64;
    let (state, _run) = coarse_state(SEED);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(SEED);
    let mut fcfg = FbmUpscaleConfig::c1_hd_production(t);
    fcfg.erosion = None;
    fcfg.bathymetry = None;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fcfg).heightmap;

    // Recommended: A_c=50, D=0.05, K=1, m=0.5, n=1, iters=3.
    let cfg = StreamPowerConfig {
        k: 1.0,
        m: 0.5,
        n: 1.0,
        dt: 1.0,
        iterations: 3,
        sea_level: SEA,
        diffusion: 0.05,
        diffusion_substeps: 4,
        min_area_cells: 50.0,
        threshold: 0.0,
        cell_km: 1.0,
        depth_scale_m: 5000.0,
        ..Default::default()
    };
    let sp = incise(&fbm, &cfg);
    let couple = 400.0f32 / 1024.0; // depth_scale ∝ domain (400 km) → relief ×0.39

    eprintln!("\n=== TASK 4 — recommended config, coupled vs uncoupled (domain 400 km) ===");
    let (tab, _) = per_order_incision(&fbm, &sp, &ss);
    eprint!("  per-order UNCOUPLED: ");
    for (o, m) in &tab {
        eprint!("S{o}={m:.0}m ");
    }
    eprintln!(
        "\n  per-order COUPLED  : {}",
        tab.iter().map(|(o, m)| format!("S{o}={:.0}m", m * couple)).collect::<Vec<_>>().join(" ")
    );
    eprintln!(
        "  → coupling scales incision ×{couple:.2}; trunks S3/S4 fall to ~{:.0}/{:.0} m — marginal for city valleys.",
        tab.iter().find(|(o, _)| *o == 3).map(|(_, m)| m * couple).unwrap_or(0.0),
        tab.iter().find(|(o, _)| *o == 4).map(|(_, m)| m * couple).unwrap_or(0.0)
    );

    let base = ss.depth_scale_m as f32;
    let s_un = slope_shares(&sp, SEA, 400.0, base);
    let s_co = slope_shares(&sp, SEA, 400.0, base * couple);
    eprintln!(
        "  slopes >15/30/45°: uncoupled {:.0}/{:.0}/{:.0}% → coupled {:.0}/{:.0}/{:.0}%",
        s_un.0 * 100.0,
        s_un.1 * 100.0,
        s_un.2 * 100.0,
        s_co.0 * 100.0,
        s_co.1 * 100.0,
        s_co.2 * 100.0
    );

    eprintln!("\n=== TASK 5 — weak-droplet coupling (drainage relief must NOT collapse) ===");
    let rel_sp = drainage_relief_m(&sp, &ss);
    let mut weak = ero_cfg(t);
    weak.num_droplets = (0.25 * cells) as usize; // weak hillslope texture (vs 0.95 prod)
    let sp_weak = run_erosion(&sp, &weak, &seed, |_, _, _| true).heightmap;
    let full = run_erosion(&sp, &ero_cfg(t), &seed, |_, _, _| true).heightmap;
    eprintln!(
        "  drainage relief: SP alone {rel_sp:.0} m → SP+weak(0.25/cell) {:.0} m → SP+full(0.95) {:.0} m",
        drainage_relief_m(&sp_weak, &ss),
        drainage_relief_m(&full, &ss),
    );
    eprintln!(
        "  (weak must stay near SP-alone; full collapses it — the reason droplets can't run at full strength)"
    );
}

/// TASK 3 — coupled regime split: hillslope diffusion (A < A_c) + stream power
/// (A ≥ A_c), interleaved. Does it make the per-order ordering monotone AND open the
/// valley cross-sections (V-walls) that plain SP + post-smooth could not? 1024².
#[test]
#[ignore]
fn stream_power_regime() {
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let t = 1024usize;
    let (state, _run) = coarse_state(SEED);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(SEED);
    let mut fcfg = FbmUpscaleConfig::c1_hd_production(t);
    fcfg.erosion = None;
    fcfg.bathymetry = None;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fcfg).heightmap;
    let cfg = |ac: f32, d: f32| StreamPowerConfig {
        k: 1.0,
        m: 0.5,
        n: 1.0,
        dt: 1.0,
        iterations: 3,
        sea_level: SEA,
        diffusion: d,
        diffusion_substeps: 4,
        min_area_cells: ac,
        threshold: 0.0,
        cell_km: 1.0,
        depth_scale_m: 5000.0,
        ..Default::default()
    };
    let probes = channel_probes(&incise(&fbm, &cfg(0.0, 0.0)), &ss, 12);
    let (fd, _) = cross_section_stats(&fbm, &ss, &probes);
    let (sd, _) = cross_section_stats(&incise(&fbm, &cfg(0.0, 0.0)), &ss, &probes);
    eprintln!("\n=== TASK 3 — coupled regime split (hillslope diffusion A<A_c + channel SP) ===");
    eprintln!("  reference cross-sections: FBM {fd:.0} m, SP-only(no diffusion) {sd:.0} m");
    eprintln!("  A_c cells | D    | per-order incision | X-sect depth | relief | ordering");
    for (ac, d) in [(50.0f32, 0.05f32), (100.0, 0.05), (100.0, 0.1), (200.0, 0.05), (100.0, 0.2)] {
        let out = incise(&fbm, &cfg(ac, d));
        let (tab, mono) = per_order_incision(&fbm, &out, &ss);
        let (dep, _) = cross_section_stats(&out, &ss, &probes);
        let rel = drainage_relief_m(&out, &ss);
        eprintln!(
            "  {ac:>5.0}    | {d:.2} | {:<34} | {dep:>4.0} m | {rel:>4.0} m | {}",
            fmt_orders(&tab),
            if mono { "MONOTONE ✓" } else { "inverted" },
        );
    }
    eprintln!(
        "  (want: monotone ordering, S1 tens of m, X-sect depth ABOVE the {sd:.0} m SP-only floor = V-walls open)"
    );
}

#[cfg(windows)]
fn peak_ws_mb() -> u64 {
    #[repr(C)]
    #[derive(Default)]
    struct Pmc {
        cb: u32,
        page_fault_count: u32,
        peak_ws: usize,
        ws: usize,
        qppp: usize,
        qpp: usize,
        qpnpp: usize,
        qnpp: usize,
        pagefile: usize,
        peak_pagefile: usize,
    }
    unsafe extern "system" {
        fn GetCurrentProcess() -> isize;
        fn K32GetProcessMemoryInfo(process: isize, counters: *mut Pmc, cb: u32) -> i32;
    }
    let mut pmc = Pmc { cb: std::mem::size_of::<Pmc>() as u32, ..Default::default() };
    unsafe {
        K32GetProcessMemoryInfo(GetCurrentProcess(), &mut pmc, pmc.cb);
    }
    (pmc.peak_ws / (1024 * 1024)) as u64
}
#[cfg(not(windows))]
fn peak_ws_mb() -> u64 {
    0
}

/// TASK 5 — re-confirm at higher resolution on the AUTHOR'S seed: does the K=1
/// calibration still land major channels in 200–400 m, and what are the real
/// slopes? (8192² deferred: tuning K at production scale is premature while
/// headwaters over-carve — see TASK 3.) 4096², domain 400 km.
#[test]
#[ignore]
fn stream_power_reconfirm() {
    use std::time::Instant;
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    use ymir_core::tectonics_c1::land_topology::slope_shares;
    let ss = SteinSteinParams::default();
    let seed_u = 10481999410520546993u64;
    let t = 4096usize;
    let domain_km = 400.0f32;
    let (state, _run) = coarse_state(seed_u);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(seed_u);
    let mut fcfg = FbmUpscaleConfig::c1_hd_production(t);
    fcfg.erosion = None;
    fcfg.bathymetry = None;
    let ti = Instant::now();
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fcfg).heightmap;
    let fbm_ms = ti.elapsed().as_millis();

    let cfg = StreamPowerConfig {
        k: 1.0,
        m: 0.5,
        n: 1.0,
        dt: 1.0,
        iterations: 3,
        sea_level: SEA,
        diffusion: 0.0,
        diffusion_substeps: 4,
        min_area_cells: 0.0,
        threshold: 0.0,
        cell_km: 1.0,
        depth_scale_m: 5000.0,
        ..Default::default()
    };
    let ti = Instant::now();
    let sp = incise(&fbm, &cfg);
    let sp_ms = ti.elapsed().as_millis();

    eprintln!("\n=== TASK 5 — author seed {seed_u} @ {t}² (domain {domain_km} km) ===");
    eprintln!("  FBM {fbm_ms} ms; stream-power {sp_ms} ms; peak RSS {} MB", peak_ws_mb());
    eprintln!(
        "  drainage relief: FBM {:.0} m → SP {:.0} m",
        drainage_relief_m(&fbm, &ss),
        drainage_relief_m(&sp, &ss)
    );
    let (_, c0, _, m0) = structure_metrics(&fbm, &ss);
    let (_, c1, _, m1) = structure_metrics(&sp, &ss);
    eprintln!("  Strahler/confluences: FBM maxS {m0}/{c0} → SP maxS {m1}/{c1}");

    // Per-order incision — K-scaling check (major channels must stay 200–400 m).
    let fbm_m = to_metres(&fbm, &ss);
    let sp_m = to_metres(&sp, &ss);
    let dr = c1_drainage(&sp, None, &C1DrainageConfig::default(), &ss);
    let mut per: std::collections::HashMap<u8, Vec<f32>> = std::collections::HashMap::new();
    for s in &dr.rivers.segments {
        for &(x, y) in &s.points {
            let k = y as usize * sp.width + x as usize;
            per.entry(s.strahler_order).or_default().push(fbm_m.data[k] - sp_m.data[k]);
        }
    }
    eprint!("  incision by order @ {t}²: ");
    let mut orders: Vec<u8> = per.keys().copied().collect();
    orders.sort();
    for o in orders {
        let v = per.get_mut(&o).unwrap();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        eprint!("S{o}={:.0}m ", v[v.len() / 2]);
    }
    eprintln!(
        "\n  (vs 1024²: S1=379 S2=361 S3=285 S4=130 — resolution-independent ⇒ K holds; scaling ⇒ re-anchor)"
    );

    // TASK 6 at production-ish scale: land slope shares, FBM vs SP, coupled/uncoupled.
    let base = ss.depth_scale_m as f32;
    let sf = slope_shares(&fbm, SEA, domain_km, base);
    let ss_sp = slope_shares(&sp, SEA, domain_km, base);
    let ss_cpl = slope_shares(&sp, SEA, domain_km, base * domain_km / 1024.0);
    eprintln!("  slope >15/30/45° @ domain {domain_km} km:");
    eprintln!(
        "    FBM (uncoupled):      {:>4.1}/{:>4.1}/{:>4.1}%",
        sf.0 * 100.0,
        sf.1 * 100.0,
        sf.2 * 100.0
    );
    eprintln!(
        "    SP  (uncoupled):      {:>4.1}/{:>4.1}/{:>4.1}%",
        ss_sp.0 * 100.0,
        ss_sp.1 * 100.0,
        ss_sp.2 * 100.0
    );
    eprintln!(
        "    SP  (coupled depth):  {:>4.1}/{:>4.1}/{:>4.1}%  (depth_scale ∝ domain — the compounding case)",
        ss_cpl.0 * 100.0,
        ss_cpl.1 * 100.0,
        ss_cpl.2 * 100.0
    );
}

/// Fixed cross-section probes: (centre cell, perpendicular unit) for the N
/// highest-order channels of `field`, reused across configs so profiles compare.
fn channel_probes(
    field: &GridF32,
    ss: &SteinSteinParams,
    n: usize,
) -> Vec<((usize, usize), (f32, f32))> {
    let dr = c1_drainage(field, None, &C1DrainageConfig::default(), ss);
    let mut segs: Vec<_> = dr.rivers.segments.iter().filter(|s| s.points.len() >= 5).collect();
    segs.sort_by_key(|s| std::cmp::Reverse(s.strahler_order));
    segs.iter()
        .take(n)
        .map(|s| {
            let mid = s.points.len() / 2;
            let (ax, ay) = s.points[mid - 1];
            let (bx, by) = s.points[mid + 1];
            let (tx, ty) = (bx as f32 - ax as f32, by as f32 - ay as f32);
            let tl = (tx * tx + ty * ty).sqrt().max(1e-6);
            ((s.points[mid].0 as usize, s.points[mid].1 as usize), (-ty / tl, tx / tl))
        })
        .collect()
}

/// Median cross-section (depth m below rim, floor half-width cells) at fixed probes.
fn cross_section_stats(
    field: &GridF32,
    ss: &SteinSteinParams,
    probes: &[((usize, usize), (f32, f32))],
) -> (f32, f32) {
    let (w, h) = (field.width, field.height);
    let fm = to_metres(field, ss);
    let (mut depths, mut widths) = (Vec::new(), Vec::new());
    for &((cx, cy), (px, py)) in probes {
        let mut xs = Vec::new();
        for o in -8i32..=8 {
            let sx = (cx as f32 + px * o as f32).round().clamp(0.0, w as f32 - 1.0) as i32;
            let sy = (cy as f32 + py * o as f32).round().clamp(0.0, h as f32 - 1.0) as i32;
            xs.push(fm.get(sx, sy));
        }
        let bottom = xs.iter().cloned().fold(f32::MAX, f32::min);
        let rim = xs[0].max(xs[16]);
        depths.push(rim - bottom);
        // floor half-width: cells within 20% of the depth above the bottom.
        let thr = bottom + 0.2 * (rim - bottom);
        widths.push(xs.iter().filter(|&&v| v <= thr).count() as f32 / 2.0);
    }
    let med = |v: &mut Vec<f32>| {
        if v.is_empty() {
            return 0.0;
        }
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v[v.len() / 2]
    };
    (med(&mut depths), med(&mut widths))
}

/// TASK 2/3/6 — stream-power tuning on the real FBM field: diffusion sweep (judged
/// on fixed cross-sections), per-Strahler-order incision + m sweep, and land-slope
/// shares before/after. 1024², read-only.
#[test]
#[ignore]
fn stream_power_tuning() {
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    use ymir_core::tectonics_c1::land_topology::slope_shares;
    let ss = SteinSteinParams::default();
    let t = 1024usize;
    let (state, _run) = coarse_state(SEED);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(SEED);
    let mut fcfg = FbmUpscaleConfig::c1_hd_production(t);
    fcfg.erosion = None;
    fcfg.bathymetry = None;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fcfg).heightmap;
    let sp = |k: f32, m: f32, d: f32| StreamPowerConfig {
        k,
        m,
        n: 1.0,
        dt: 1.0,
        iterations: 3,
        sea_level: SEA,
        diffusion: d,
        diffusion_substeps: 4,
        min_area_cells: 0.0,
        threshold: 0.0,
        cell_km: 1.0,
        depth_scale_m: 5000.0,
        ..Default::default()
    };
    let base = incise(&fbm, &sp(1.0, 0.5, 0.0)); // fix probes on the D=0 carved field
    let probes = channel_probes(&base, &ss, 12);

    eprintln!("\n=== TASK 2 — diffusion sweep (K=1,m=0.5,n=1,iters=3), fixed cross-sections ===");
    let (fd, fw) = cross_section_stats(&fbm, &ss, &probes);
    eprintln!("  FBM (no incision):  depth {fd:>4.0} m, floor half-width {fw:.1} cells");
    for d in [0.0f32, 0.01, 0.03, 0.05, 0.1] {
        let f = incise(&fbm, &sp(1.0, 0.5, d));
        let (dep, wid) = cross_section_stats(&f, &ss, &probes);
        eprintln!(
            "  D={d:>4.2}: cross-section depth {dep:>4.0} m, floor half-width {wid:.1} cells"
        );
    }
    eprintln!(
        "  (want depth preserved AND walls opened: floor width rising a little, not collapsing depth)"
    );

    eprintln!("\n=== TASK 3 — median incision per Strahler order + m sweep ===");
    let fbm_m = to_metres(&fbm, &ss);
    for m in [0.4f32, 0.5, 0.6] {
        let f = incise(&fbm, &sp(1.0, m, 0.0));
        let fm = to_metres(&f, &ss);
        let dr = c1_drainage(&f, None, &C1DrainageConfig::default(), &ss);
        let mut per: std::collections::HashMap<u8, Vec<f32>> = std::collections::HashMap::new();
        for s in &dr.rivers.segments {
            for &(x, y) in &s.points {
                let k = y as usize * f.width + x as usize;
                per.entry(s.strahler_order).or_default().push(fbm_m.data[k] - fm.data[k]);
            }
        }
        eprint!("  m={m:.1}: incision by order ");
        let mut orders: Vec<u8> = per.keys().copied().collect();
        orders.sort();
        for o in orders {
            let v = per.get_mut(&o).unwrap();
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            eprint!("S{o}={:.0}m ", v[v.len() / 2]);
        }
        eprintln!();
    }
    eprintln!("  (small orders should be tens of m, not 300; m governs the large/small split)");

    eprintln!("\n=== TASK 6 — land slope shares (>15/30/45°) before/after (domain 1024 km) ===");
    let base_scale = ss.depth_scale_m as f32;
    let s_fbm = slope_shares(&fbm, SEA, 1024.0, base_scale);
    let s_sp = slope_shares(&incise(&fbm, &sp(1.0, 0.5, 0.05)), SEA, 1024.0, base_scale);
    eprintln!(
        "  FBM:                {:>4.1}% / {:>4.1}% / {:>4.1}%",
        s_fbm.0 * 100.0,
        s_fbm.1 * 100.0,
        s_fbm.2 * 100.0
    );
    eprintln!(
        "  stream-power (D=.05): {:>4.1}% / {:>4.1}% / {:>4.1}%",
        s_sp.0 * 100.0,
        s_sp.1 * 100.0,
        s_sp.2 * 100.0
    );
    eprintln!(
        "  (slopes rise with relief; feeds the cliff/unbuildable + vertical-scale decision — measure only)"
    );
}

/// Drainage relief: median, over top-1% accumulation LAND cells, of (max altitude
/// in an 11×11 window − the cell), in metres. A carved valley network sits well
/// BELOW its interfluves ⇒ this rises. Fixed-location (unlike the single-channel
/// cross-section), so it's comparable across configs. The right incision metric —
/// the local-minimum "carved%" measures PITS (a drained channel is not a local min).
fn drainage_relief_m(field: &GridF32, ss: &SteinSteinParams) -> f32 {
    let dr = c1_drainage(field, None, &C1DrainageConfig::default(), ss);
    let acc = &dr.flow.accumulation;
    let (w, h) = (field.width, field.height);
    let mut accs = acc.data.clone();
    accs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let thr = accs[(accs.len() as f64 * 0.99) as usize];
    let mut rel: Vec<f32> = Vec::new();
    let r = 5i32;
    for y in 0..h {
        for x in 0..w {
            let k = y * w + x;
            if acc.data[k] < thr || field.data[k] <= SEA {
                continue;
            }
            let mut mx = field.data[k];
            for dy in -r..=r {
                for dx in -r..=r {
                    let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                    if nx >= 0 && ny >= 0 && nx < w as i32 && ny < h as i32 {
                        mx = mx.max(field.data[ny as usize * w + nx as usize]);
                    }
                }
            }
            rel.push(
                c1_altitude_norm_to_metres(mx, ss) - c1_altitude_norm_to_metres(field.data[k], ss),
            );
        }
    }
    if rel.is_empty() {
        return 0.0;
    }
    rel.sort_by(|a, b| a.partial_cmp(b).unwrap());
    rel[rel.len() / 2]
}

/// PART B (relief) — the RIGHT incision metric across the coupling, on the real FBM
/// field. Drainage relief = how deep channels sit below their interfluves (median,
/// top-1% flow cells, 11×11 window). Rises = real valleys.
#[test]
#[ignore]
fn stream_power_relief() {
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise};
    let ss = SteinSteinParams::default();
    let t = 1024usize;
    let (state, _run) = coarse_state(SEED);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(SEED);
    let mut fcfg = FbmUpscaleConfig::c1_hd_production(t);
    fcfg.erosion = None;
    fcfg.bathymetry = None;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fcfg).heightmap;
    let sp = StreamPowerConfig { k: 1.0, iterations: 3, sea_level: SEA, ..Default::default() };
    let spd = StreamPowerConfig { diffusion: 0.4, ..sp.clone() };

    eprintln!("\n=== PART B (relief) — drainage relief (m), the incision metric ({t}²) ===");
    let base = drainage_relief_m(&fbm, &ss);
    eprintln!("  FBM baseline:              {base:.0} m");
    let a = drainage_relief_m(&incise(&fbm, &sp), &ss);
    eprintln!("  stream-power alone:        {a:.0} m  (Δ {:+.0})", a - base);
    let b = drainage_relief_m(&incise(&fbm, &spd), &ss);
    eprintln!("  stream-power + diffusion:  {b:.0} m  (Δ {:+.0})", b - base);
    let dro = run_erosion(&fbm, &ero_cfg(t), &seed, |_, _, _| true).heightmap;
    let c = drainage_relief_m(&dro, &ss);
    eprintln!("  droplets alone (prod):     {c:.0} m  (Δ {:+.0})", c - base);
    let both = run_erosion(&incise(&fbm, &spd), &ero_cfg(t), &seed, |_, _, _| true).heightmap;
    let e = drainage_relief_m(&both, &ss);
    eprintln!("  both (SP+diff then droplets): {e:.0} m  (Δ {:+.0})", e - base);
    eprintln!("  (higher = channels sit deeper below interfluves = real dendritic valleys)");
}

/// PART B — routed stream-power incision prototype (Braun & Willett). Calibrates K,
/// then compares stream-power alone / droplets alone / both on the real FBM field
/// with STRUCTURE metrics + runtime, plus the staleness (drainage↔incision) check.
#[test]
#[ignore]
fn stream_power_prototype() {
    use std::time::Instant;
    use ymir_core::erosion::stream_power::{StreamPowerConfig, incise, incise_with_progress};
    let ss = SteinSteinParams::default();
    let t = 1024usize;
    let (state, _run) = coarse_state(SEED);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(SEED);
    let mut fcfg = FbmUpscaleConfig::c1_hd_production(t);
    fcfg.erosion = None;
    fcfg.bathymetry = None;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fcfg).heightmap;
    let (_, c0, carved0, m0) = structure_metrics(&fbm, &ss);
    let prof0 = channel_incision_profile(&fbm, &ss);
    eprintln!("\n=== PART B — stream-power incision ({t}², real FBM field) ===");
    eprintln!(
        "  FBM baseline: carved {:.0}%, maxStrahler {m0}, confluences {c0}, channel profile {prof0:.0} m",
        carved0 * 100.0
    );

    // K calibration — pick K whose median channel incision is plausible (~200–400 m
    // over 4 iterations), NOT chosen for appearance.
    eprintln!(
        "  K calibration (m=0.5, n=1, dt=1, iters=4): K | median channel incision | carvedΔ | maxS"
    );
    let mut chosen_k = 1.0f32;
    let mut best_gap = f32::MAX;
    for k in [0.25f32, 0.5, 1.0, 2.0, 4.0] {
        let cfg = StreamPowerConfig { k, iterations: 4, sea_level: SEA, ..Default::default() };
        let out = incise(&fbm, &cfg);
        let med = median_channel_incision_m(&fbm, &out, &ss);
        let (_, _, carved, maxs) = structure_metrics(&out, &ss);
        eprintln!("    {k:>5.2} | {med:>6.0} m | {:>+3.0}% | {maxs}", (carved - carved0) * 100.0);
        let gap = (med - 300.0).abs();
        if (200.0..=400.0).contains(&med) && gap < best_gap {
            best_gap = gap;
            chosen_k = k;
        }
    }
    eprintln!(
        "  → chosen K = {chosen_k} (median channel incision closest to ~300 m within [200,400])"
    );

    // Coupling comparison: stream-power alone / droplets alone / both.
    let spcfg =
        StreamPowerConfig { k: chosen_k, iterations: 4, sea_level: SEA, ..Default::default() };
    let spcfg_diff = StreamPowerConfig { diffusion: 0.4, ..spcfg.clone() };

    let run = |label: &str, out: &GridF32, dt_ms: u128| {
        let (hist, confl, carved, maxs) = structure_metrics(out, &ss);
        let prof = channel_incision_profile(out, &ss);
        eprintln!(
            "  [{label:<22}] carvedΔ {:>+3.0}% ({:>3.0}%) | maxS {maxs} | confl {confl:>5} | profile {prof:>4.0} m | {dt_ms} ms | hist {:?}",
            (carved - carved0) * 100.0,
            carved * 100.0,
            &hist[1..=maxs.max(1) as usize],
        );
    };

    let ti = Instant::now();
    let sp = incise(&fbm, &spcfg);
    run("stream-power alone", &sp, ti.elapsed().as_millis());

    let ti = Instant::now();
    let sp_d = incise(&fbm, &spcfg_diff);
    run("stream-power + diffusion", &sp_d, ti.elapsed().as_millis());

    let ti = Instant::now();
    let dr = run_erosion(&fbm, &ero_cfg(t), &seed, |_, _, _| true).heightmap;
    run("droplets alone (prod)", &dr, ti.elapsed().as_millis());

    let ti = Instant::now();
    let both = run_erosion(&sp, &ero_cfg(t), &seed, |_, _, _| true).heightmap;
    run("both (SP then droplets)", &both, ti.elapsed().as_millis());

    // Staleness — does the network reorganise across drainage↔incision iterations?
    eprintln!("  staleness (per-iteration maxStrahler / carved% after each incision pass):");
    let stale_cfg =
        StreamPowerConfig { k: chosen_k, iterations: 6, sea_level: SEA, ..Default::default() };
    let mut last_maxs = 0u8;
    let mut converged_at = 0usize;
    incise_with_progress(&fbm, &stale_cfg, &mut |iter, f| {
        let (_, _, carved, maxs) = structure_metrics(f, &ss);
        eprintln!("    iter {}: maxStrahler {maxs}, carved {:.0}%", iter + 1, carved * 100.0);
        if maxs == last_maxs && converged_at == 0 && iter > 0 {
            converged_at = iter;
        }
        last_maxs = maxs;
    });
    eprintln!(
        "  → {}",
        if converged_at > 0 {
            format!("Strahler stable from iteration {}", converged_at + 1)
        } else {
            "still moving at iter 6".into()
        }
    );
    eprintln!(
        "  baselines to beat (ADR 0001): droplets alone carved ~11%, 4/cell ~18% (maxS 6→3 = fragmented)."
    );
}

fn emerged_frac(field: &GridF32) -> f32 {
    field.data.iter().filter(|&&v| v > SEA).count() as f32 / field.data.len() as f32
}

/// Total-deposition-vs-distance-to-coast profile (≤5 / ≤20 / ≤50 cells) as a % of
/// all deposition (before→after), on the exported grid.
fn deposition_locality(before: &GridF32, after: &GridF32) -> (f32, f32, f32) {
    let cd = coast_distance(after);
    let (mut tot, mut d5, mut d20, mut d50) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    for k in 0..before.data.len() {
        let dep = (after.data[k] - before.data[k]).max(0.0) as f64;
        if dep <= 0.0 {
            continue;
        }
        tot += dep;
        let dc = cd[k];
        if dc <= 5 {
            d5 += dep;
        }
        if dc <= 20 {
            d20 += dep;
        }
        if dc <= 50 {
            d50 += dep;
        }
    }
    let t = tot.max(1e-9);
    ((d5 / t * 100.0) as f32, (d20 / t * 100.0) as f32, (d50 / t * 100.0) as f32)
}

/// TASK 1 + 2 — inclined plane (parallel flow, sea at the low edge) at the CURRENT
/// production density (0.95/cell), sweeping the coastal SINK fraction f. Corrected
/// metrics: carved% is reported as the DIFFERENTIAL (after − before on the same
/// field) so it measures INCISION, not the initial-noise floor. Read-only; the sink
/// is a parameter (default f=1.0 = current behaviour). Also reports the emerged-
/// fraction drift on the real coarse→FBM field.
#[test]
#[ignore]
fn inclined_plane_structure() {
    let ss = SteinSteinParams::default();
    let n = 512usize;
    // Tilt along +y: top high land (~0.82), bottom below sea (~0.48) so rivers have
    // an outlet. Small deterministic value-noise seeds channels.
    let mut plane = GridF32::new(n, n, 0.0);
    for y in 0..n {
        for x in 0..n {
            let t = y as f32 / (n as f32 - 1.0);
            let s = (x as f32 * 127.1 + y as f32 * 311.7).sin() * 43758.547;
            let noise = (s - s.floor()) - 0.5;
            plane.set(x, y, 0.48 + 0.34 * (1.0 - t) + 0.008 * noise);
        }
    }
    let cells = (n * n) as f64;
    let mut base = FbmUpscaleConfig::c1_hd_production(n).erosion.unwrap();
    base.num_droplets = (0.95 * cells) as usize; // production density

    let (h0, c0, carved_before, m0) = structure_metrics(&plane, &ss);
    let em0 = emerged_frac(&plane);
    eprintln!(
        "\n=== TASK 1+2 — inclined plane {n}², coastal-sink f sweep @ 0.95 droplets/cell ==="
    );
    eprintln!(
        "  baseline (no erosion): carved(local-min of top-1% flow) {:.0}%, maxStrahler {m0}, \
         confluences {c0}, emerged {:.1}%, order hist {:?}",
        carved_before * 100.0,
        em0 * 100.0,
        &h0[1..=m0.max(1) as usize],
    );
    eprintln!(
        "  f=deposit | net% | carvedΔ (after) | maxS | confl | dep≤5/≤20/≤50 | emerged→ | order hist"
    );
    // f = 1.0 (A: current), 0.5, 0.25, 0.1, 0.0 (total sink).
    for f in [1.0f32, 0.5, 0.25, 0.1, 0.0] {
        let mut cfg = base.clone();
        cfg.coastal_deposit_fraction = f;
        let seed = WorldSeed::new(7);
        let out = run_erosion(&plane, &cfg, &seed, |_, _, _| true).heightmap;
        let (e, d) = erosion_balance(&plane, &out);
        let (hist, confl, carved_after, maxs) = structure_metrics(&out, &ss);
        let (d5, d20, d50) = deposition_locality(&plane, &out);
        let em = emerged_frac(&out);
        eprintln!(
            "  {f:>7.2}  | {:>4.0} | {:>+4.0}% ({:>3.0}%) | {maxs:>4} | {confl:>5} | {d5:>3.0}/{d20:>3.0}/{d50:>3.0}% | {:>4.1}→{:>4.1}% | {:?}",
            (e - d) / e.max(1.0) * 100.0,
            (carved_after - carved_before) * 100.0,
            carved_after * 100.0,
            em0 * 100.0,
            em * 100.0,
            &hist[1..=maxs.max(1) as usize],
        );
    }
    eprintln!(
        "  (carvedΔ = incision above the noise floor; delta/beach survives while dep≤5c > 0. \
         f=1.0 is current production; f=0.0 is the total sink.)"
    );

    // TASK 2 impact on the REAL field: emerged fraction coarse → FBM → eroded, for a
    // few sink fractions (1024², production density), to size the island-shrink risk.
    eprintln!("\n=== emerged-fraction drift on the real coarse→FBM field (1024², 0.95/cell) ===");
    let t = 1024usize;
    let (state, _run) = coarse_state(SEED);
    let coarse = c1_coarse_normalized_altitude(&state, &IsostasyConfig::c1_default(), &ss, None);
    let seed = WorldSeed::new(SEED);
    let mut fcfg = FbmUpscaleConfig::c1_hd_production(t);
    fcfg.erosion = None;
    fcfg.bathymetry = None;
    let fbm = upscale_with_fbm(&coarse, SEA, &seed, &fcfg).heightmap;
    eprintln!(
        "  coarse emerged {:.1}%  →  FBM emerged {:.1}%",
        emerged_frac(&coarse) * 100.0,
        emerged_frac(&fbm) * 100.0,
    );
    for f in [1.0f32, 0.25, 0.0] {
        let mut ec = ero_cfg(t);
        ec.coastal_deposit_fraction = f;
        let out = run_erosion(&fbm, &ec, &seed, |_, _, _| true).heightmap;
        let (e, d) = erosion_balance(&fbm, &out);
        let (_, _, carved, maxs) = structure_metrics(&out, &ss);
        eprintln!(
            "  f={f:.2}: eroded emerged {:.1}%  (net {:>+4.0}%, carved {:.0}%, maxS {maxs})",
            emerged_frac(&out) * 100.0,
            (e - d) / e.max(1.0) * 100.0,
            carved * 100.0,
        );
    }
}

fn std_of(field: &GridF32) -> f32 {
    let (w, h) = (field.width, field.height);
    let mut acc = 0.0f64;
    let mut n = 0u64;
    for y in 1..h as i32 - 1 {
        for x in 1..w as i32 - 1 {
            let c = field.get(x, y);
            let m = (field.get(x - 1, y)
                + field.get(x + 1, y)
                + field.get(x, y - 1)
                + field.get(x, y + 1))
                / 4.0;
            acc += ((c - m) as f64).powi(2);
            n += 1;
        }
    }
    (acc / n as f64).sqrt() as f32
}
