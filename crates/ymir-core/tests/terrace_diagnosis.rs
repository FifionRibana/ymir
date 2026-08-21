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
    let (mn, mx) = eroded_m.data.iter().fold((f32::MAX, f32::MIN), |(a, b), &v| (a.min(v), b.max(v)));
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
            for (dx, dy) in [(-1i32, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)] {
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
            eprint!("  cross-section ±12 cells across a Strahler-{} channel (m): ", seg.strahler_order);
            let mut xs = Vec::new();
            for o in -12i32..=12 {
                let x = (px as i32 + o).clamp(0, w as i32 - 1);
                xs.push(eroded_m.get(x, py as i32));
            }
            let cmin = xs.iter().cloned().fold(f32::MAX, f32::min);
            for v in &xs {
                eprint!("{:.0} ", v);
            }
            eprintln!("\n    channel depth below rim: {:.0} m (V/U valley if clearly incised)", xs[0].max(xs[24]) - cmin);
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
        if std_after > std_before * 1.3 { "erosion CARVED channels" } else { "little/no channelisation" },
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
    eprintln!("  (channelisation = roughness ratio ≫1 and carved% rising; production row is the current default)");
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
            for (dx, dy) in [(-1i32, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)] {
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
    eprintln!("  (delta/beach survives while dep≤5c > 0; carvedΔ>0 with a deep histogram = real incision.)");
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
    let vals: Vec<f64> = (0..ny).flat_map(|j| (0..nx).map(move |i| (i, j))).map(|(i, j)| s.get(i, j)).collect();
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
        at_heq, at_heq as f32 / n * 100.0,
        at_hmax, at_hmax as f32 / n * 100.0,
        above_heq, above_heq as f32 / n * 100.0,
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
            eprint!("{:.2}({:.0}%) ", mn + (b as f64 + 0.5) / NB as f64 * span, c as f32 / n * 100.0);
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
            eprint!("{:.0}m({:.0}%) ", lmn + (b as f32 + 0.5) / AB as f32 * aspan, c as f32 / land_m.len() as f32 * 100.0);
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
        pct(&cra, 0.1), pct(&cra, 0.5), pct(&cra, 0.9),
        pct(&noncra, 0.1), pct(&noncra, 0.5), pct(&noncra, 0.9),
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
    eprintln!("\n=== TASK 1+2 — inclined plane {n}², coastal-sink f sweep @ 0.95 droplets/cell ===");
    eprintln!(
        "  baseline (no erosion): carved(local-min of top-1% flow) {:.0}%, maxStrahler {m0}, \
         confluences {c0}, emerged {:.1}%, order hist {:?}",
        carved_before * 100.0,
        em0 * 100.0,
        &h0[1..=m0.max(1) as usize],
    );
    eprintln!("  f=deposit | net% | carvedΔ (after) | maxS | confl | dep≤5/≤20/≤50 | emerged→ | order hist");
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
            let m = (field.get(x - 1, y) + field.get(x + 1, y) + field.get(x, y - 1) + field.get(x, y + 1)) / 4.0;
            acc += ((c - m) as f64).powi(2);
            n += 1;
        }
    }
    (acc / n as f64).sqrt() as f32
}
