//! ADR Finding 48 — the COASTLINE CONTOUR: barb metric by slope class, before and after
//! smoothing, plus the sea-level offset that must not drift.
//!
//! Marching squares ALREADY interpolates sub-cell along each crossed edge, so the barbs are not
//! blockiness. The diagnosis on record is GRADIENT PINNING: where the field is nearly flat at
//! the iso value, which edges get crossed is decided by tiny fluctuations, so the contour
//! zig-zags. That diagnosis makes a FALSIFIABLE PREDICTION — the fix must improve the LOW-SLOPE
//! shores most. A uniform improvement, or one concentrated on steep shores, refutes it.
//!
//! Both variants are measured on the SAME extracted contour in one run, so the comparison is
//! exact. Run:
//!   cargo test -p ymir-core --release --test coastline_barbs -- --ignored --nocapture

use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::fracture::FractureConfig;
use ymir_core::tectonics_c1::closures::lithology::LithologyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::{
    c1_altitude_norm_to_metres, upscale_from_c1_with_progress,
};
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::contour::{
    COASTLINE_SMOOTH_PASSES, Polyline, SMOOTH_RELAX_BELOW_DEG, marching_squares,
    slope_deg_to_norm_gradient, smooth_polylines_on_isoline,
};
use ymir_core::terrain::upscale::{ProductionHdOpts, production_hd_config};

const PSEED: u64 = 10_481_999_410_520_546_993;
const DOMAIN_KM: f32 = 400.0;
const SEA: f32 = 0.5;

fn production_field(target: usize) -> GridF32 {
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
        target_size: target,
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
    let (up, _) = upscale_from_c1_with_progress(
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
    );
    up.heightmap
}

fn sample(grid: &GridF32, x: f32, y: f32) -> f32 {
    let (w, h) = (grid.width as i32, grid.height as i32);
    let x = x.clamp(0.0, (w - 1) as f32);
    let y = y.clamp(0.0, (h - 1) as f32);
    let (x0, y0) = (x.floor() as i32, y.floor() as i32);
    let (x1, y1) = ((x0 + 1).min(w - 1), (y0 + 1).min(h - 1));
    let (fx, fy) = (x - x0 as f32, y - y0 as f32);
    let g = |xi: i32, yi: i32| grid.data[yi as usize * grid.width + xi as usize];
    let t = g(x0, y0) * (1.0 - fx) + g(x1, y0) * fx;
    let b = g(x0, y1) * (1.0 - fx) + g(x1, y1) * fx;
    t * (1.0 - fy) + b * fy
}

/// Local terrain slope in DEGREES at a sub-cell position — the slope-class axis.
fn slope_deg_at(grid: &GridF32, x: f32, y: f32, m_per_cell: f32, ss: &SteinSteinParams) -> f32 {
    let to_m = |v: f32| c1_altitude_norm_to_metres(v, ss);
    let gx = 0.5 * (to_m(sample(grid, x + 1.0, y)) - to_m(sample(grid, x - 1.0, y)));
    let gy = 0.5 * (to_m(sample(grid, x, y + 1.0)) - to_m(sample(grid, x, y - 1.0)));
    ((gx * gx + gy * gy).sqrt() / m_per_cell).atan().to_degrees()
}

const CLASS_EDGES: [f32; 4] = [0.5, 2.0, 5.0, 15.0];
const CLASS_NAMES: [&str; 5] = ["<0.5", "0.5-2", "2-5", "5-15", ">15"];

fn class_of(deg: f32) -> usize {
    CLASS_EDGES.iter().position(|&e| deg < e).unwrap_or(4)
}

struct Metric {
    /// Per slope class: (vertices, turns > 80 deg).
    per_class: [(usize, usize); 5],
    /// Axial concentration: |mean exp(i*4*theta)| over segment directions. 1 = every segment
    /// axis- or diagonal-aligned (the barb signature), 0 = isotropic.
    axial_r: f64,
    mean_step_cells: f64,
    vertices: usize,
    /// Median |altitude| at the vertices, in metres — 0 for a contour ON the isoline.
    offset_median_m: f32,
    offset_p90_m: f32,
}

/// `class_from` supplies the position each vertex is CLASSIFIED at. It must be the RAW
/// contour for both variants: smoothing moves a vertex, so classifying it where it ends up
/// makes the classes different POPULATIONS between before and after (measured: the <0.5° class
/// went 2903 → 1685 vertices while 0.5–2° went 1333 → 2678 — the worst vertices simply
/// migrated, so "−26 % here, +13 % there" was comparing two different sets). Classifying both
/// variants at the original position makes it a paired comparison.
fn measure(
    grid: &GridF32,
    polys: &[Polyline],
    class_from: &[Polyline],
    m_per_cell: f32,
    ss: &SteinSteinParams,
) -> Metric {
    let mut per_class = [(0usize, 0usize); 5];
    let (mut c4, mut s4, mut n_seg) = (0.0f64, 0.0f64, 0usize);
    let mut step_sum = 0.0f64;
    let mut offsets: Vec<f32> = Vec::new();
    let mut vertices = 0usize;

    for (pi, pl) in polys.iter().enumerate() {
        for w2 in pl.windows(2) {
            let (dx, dy) = (w2[1].0 - w2[0].0, w2[1].1 - w2[0].1);
            let len = (dx * dx + dy * dy).sqrt();
            if len <= 0.0 {
                continue;
            }
            let th = (dy as f64).atan2(dx as f64);
            c4 += (4.0 * th).cos();
            s4 += (4.0 * th).sin();
            n_seg += 1;
            step_sum += len as f64;
        }
        for &(x, y) in pl.iter() {
            vertices += 1;
            offsets.push(c1_altitude_norm_to_metres(sample(grid, x, y), ss));
        }
        for (wi, w3) in pl.windows(3).enumerate() {
            let (ax, ay) = (w3[1].0 - w3[0].0, w3[1].1 - w3[0].1);
            let (bx, by) = (w3[2].0 - w3[1].0, w3[2].1 - w3[1].1);
            let (la, lb) = ((ax * ax + ay * ay).sqrt(), (bx * bx + by * by).sqrt());
            if la <= 0.0 || lb <= 0.0 {
                continue;
            }
            let cos = ((ax * bx + ay * by) / (la * lb)).clamp(-1.0, 1.0);
            let turn = cos.acos().to_degrees();
            // Classify at the ORIGINAL position of this same vertex (index wi+1 of ring pi).
            let (cx, cy) = class_from
                .get(pi)
                .and_then(|r| r.get(wi + 1))
                .copied()
                .unwrap_or((w3[1].0, w3[1].1));
            let k = class_of(slope_deg_at(grid, cx, cy, m_per_cell, ss));
            per_class[k].0 += 1;
            if turn > 80.0 {
                per_class[k].1 += 1;
            }
        }
    }
    let mut abs_off: Vec<f32> = offsets.iter().map(|v| v.abs()).collect();
    abs_off.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let med = |v: &Vec<f32>, f: f32| {
        if v.is_empty() { 0.0 } else { v[((v.len() as f32 - 1.0) * f) as usize] }
    };
    // Signed median, the figure the coherence check reports.
    let mut signed = offsets.clone();
    signed.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    Metric {
        per_class,
        axial_r: if n_seg == 0 {
            0.0
        } else {
            ((c4 / n_seg as f64).powi(2) + (s4 / n_seg as f64).powi(2)).sqrt()
        },
        mean_step_cells: if n_seg == 0 { 0.0 } else { step_sum / n_seg as f64 },
        vertices,
        offset_median_m: med(&signed, 0.5),
        offset_p90_m: med(&abs_off, 0.9),
    }
}

fn report(label: &str, m: &Metric) {
    eprintln!(
        "\n{label}\n  vertices {} | axial R {:.3} | mean step {:.3} cell | sea-level offset \
         median {:+.3} m, |p90| {:.3} m",
        m.vertices, m.axial_r, m.mean_step_cells, m.offset_median_m, m.offset_p90_m
    );
    eprint!("  turns > 80 deg by slope class:");
    for (i, name) in CLASS_NAMES.iter().enumerate() {
        let (n, k) = m.per_class[i];
        eprint!("  {name}: {:.1}% ({n})", 100.0 * k as f32 / n.max(1) as f32);
    }
    eprintln!();
}

fn run(target: usize) {
    let ss = SteinSteinParams::default();
    let grid = production_field(target);
    let m_per_cell = DOMAIN_KM * 1000.0 / grid.width as f32;
    let raw = marching_squares(&grid, SEA);
    let gate = slope_deg_to_norm_gradient(
        SMOOTH_RELAX_BELOW_DEG,
        m_per_cell,
        2.0 * 1.13 * ss.depth_scale_m as f32,
    );
    let smoothed = smooth_polylines_on_isoline(&grid, SEA, &raw, COASTLINE_SMOOTH_PASSES, gate);

    eprintln!(
        "\n==========  COASTLINE BARBS — {target}^2, {:.0} m/cell, {} passes  ==========\n\
         rings/lines: {} (unchanged by smoothing: {})",
        m_per_cell,
        COASTLINE_SMOOTH_PASSES,
        raw.len(),
        smoothed.len() == raw.len()
    );
    let before = measure(&grid, &raw, &raw, m_per_cell, &ss);
    let after = measure(&grid, &smoothed, &raw, m_per_cell, &ss);
    report("BEFORE (raw marching squares)", &before);
    report("AFTER  (smoothed + reprojected)", &after);

    // THE AGGREGATE THE PER-CLASS TABLE CAN HIDE. A per-class table shows where things moved;
    // it does not say whether anything was REMOVED. If barbs merely migrate between classes the
    // total is flat and the "improvement" is a redistribution.
    let tot = |m: &Metric| -> (usize, usize) {
        m.per_class.iter().fold((0, 0), |(n, k), &(a, b)| (n + a, k + b))
    };
    let (nb, kb) = tot(&before);
    let (na, ka) = tot(&after);
    eprintln!(
        "\nTOTAL turns > 80 deg: {kb} of {nb} ({:.2}%) -> {ka} of {na} ({:.2}%) = {:+.1}% relative",
        100.0 * kb as f32 / nb.max(1) as f32,
        100.0 * ka as f32 / na.max(1) as f32,
        100.0 * (ka as f32 - kb as f32) / kb.max(1) as f32
    );

    eprintln!("\nPER-CLASS CHANGE — the falsifiable prediction is that LOW slopes improve MOST");
    eprintln!("  {:<8} {:>10} {:>10} {:>12}", "class", "before %", "after %", "relative");
    for (i, name) in CLASS_NAMES.iter().enumerate() {
        let (nb, kb) = before.per_class[i];
        let (na, ka) = after.per_class[i];
        let b = 100.0 * kb as f32 / nb.max(1) as f32;
        let a = 100.0 * ka as f32 / na.max(1) as f32;
        eprintln!(
            "  {name:<8} {b:>9.1}% {a:>9.1}% {:>11.0}%",
            if b > 0.0 { 100.0 * (a - b) / b } else { 0.0 }
        );
    }
}

#[test]
#[ignore]
fn coastline_barbs_2048() {
    run(2048);
}

#[test]
#[ignore]
fn coastline_barbs_8192() {
    run(8192);
}
