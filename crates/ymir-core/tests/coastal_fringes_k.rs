//! ADR Finding 50 — do the coastal FRINGES coincide with `K` class boundaries?
//!
//! The author's toggle bisection attributed the fringes mainly to C-3 (lithology on, C-3b and
//! H-1 off is markedly more fringed than everything off), with a residue that survives every
//! closure being off. This measures the dominant contributor directly, on the SAME `K` field
//! production incises with (`production_k_field`), so it is not a reconstruction.
//!
//! Two things `K` could be, and they call for different remedies:
//!   - CLASS-WISE / BINARY at HD  → smooth the FIELD over a stated distance in metres;
//!   - already CONTINUOUS at HD   → the cause is elsewhere and the remedy is heavier.
//!
//! Run: cargo test -p ymir-core --release --test coastal_fringes_k -- --ignored --nocapture

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
    production_k_field, upscale_from_c1_with_progress,
};
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::contour::marching_squares;
use ymir_core::terrain::upscale::{FbmUpscaleConfig, ProductionHdOpts, production_hd_config};

const PSEED: u64 = 10_481_999_410_520_546_993;
const DOMAIN_KM: f32 = 400.0;
const SEA: f32 = 0.5;

struct Built {
    height: GridF32,
    k: Vec<f32>,
    w: usize,
}

/// The production chain, returning BOTH the eroded field and the K field it was incised with.
fn build_amp(target: usize, litho: bool, fracture: bool, amplitude: f64) -> Built {
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
    let cfg: FbmUpscaleConfig = production_hd_config(&ProductionHdOpts {
        target_size: target,
        domain_km: DOMAIN_KM,
        depth_scale_m: ss.depth_scale_m as f32,
        sample_origin: [0.0, 0.578_125],
        sample_size: 1.0,
        amplitude_base: amplitude,
        mfd_p: 2.0,
        lithology: LithologyConfig {
            enabled: litho,
            soft_multiplier: 10.0,
            volcanic_multiplier: 3.0,
            rift_age_threshold: 1.0,
        },
        fracture: FractureConfig {
            enabled: fracture,
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
    let w = up.heightmap.width;
    let h = up.heightmap.height;
    let k = production_k_field(&state, Some(&kin), &cfg, &edifices, &volc, w, h)
        .unwrap_or_else(|| vec![1.0; w * h]);
    Built { height: up.heightmap, k, w }
}

/// The shipped amplitude.
fn build(target: usize, litho: bool, fracture: bool) -> Built {
    build_amp(target, litho, fracture, 0.04)
}

fn at(v: &[f32], w: usize, x: f32, y: f32) -> f32 {
    let xi = (x.round() as isize).clamp(0, w as isize - 1) as usize;
    let yi = (y.round() as isize).clamp(0, (v.len() / w) as isize - 1) as usize;
    v[yi * w + xi]
}

/// |∇K| per cell by central difference — the sharpness of the K transition.
fn grad_k(v: &[f32], w: usize, x: f32, y: f32) -> f32 {
    let gx = at(v, w, x + 1.0, y) - at(v, w, x - 1.0, y);
    let gy = at(v, w, x, y + 1.0) - at(v, w, x, y - 1.0);
    0.5 * (gx * gx + gy * gy).sqrt()
}

fn pct(v: &mut Vec<f32>, f: f32) -> f32 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    // ⚠️ INTEGER arithmetic. `(len - 1) as f32 * f` loses precision above f32's 24-bit
    // mantissa: at 8192² the vector holds 67 108 864 entries, 67 108 863 rounds UP to
    // 67 108 864.0, and the index is out of bounds. Found at 8192², invisible at 2048².
    if v.is_empty() {
        return 0.0;
    }
    let idx = ((v.len() - 1) as f64 * f as f64).round() as usize;
    v[idx.min(v.len() - 1)]
}

fn run(target: usize) {
    let km_per_cell = DOMAIN_KM / target as f32;
    eprintln!(
        "\n==========  COASTAL FRINGES vs K — {target}^2, {:.0} m/cell  ==========",
        km_per_cell * 1000.0
    );

    // ── 1. WHAT IS K AT HD? Class-wise or continuous — the question that picks the remedy.
    let b = build(target, true, true);
    let mut all: Vec<f32> = b.k.clone();
    let n = all.len();
    let uniq_1 = b.k.iter().filter(|&&v| (v - 1.0).abs() < 1e-4).count();
    eprintln!(
        "\n1. THE K FIELD OVER THE WHOLE MAP ({n} cells)\n   \
         p1 {:.3} | p25 {:.3} | p50 {:.3} | p75 {:.3} | p99 {:.3} | max {:.3}\n   \
         cells at exactly 1.000: {uniq_1} ({:.1} %)",
        pct(&mut all, 0.01),
        pct(&mut all, 0.25),
        pct(&mut all, 0.50),
        pct(&mut all, 0.75),
        pct(&mut all, 0.99),
        pct(&mut all, 1.0),
        100.0 * uniq_1 as f32 / n as f32
    );

    // ── 2. K ALONG THE SHORE, and the fringe correlation.
    let polys = marching_squares(&b.height, SEA);
    let mut k_shore: Vec<f32> = Vec::new();
    let mut gk_shore: Vec<f32> = Vec::new();
    // Paired per vertex: is this a fringe (turn > 80°) and what is |∇K| there?
    let mut pairs: Vec<(bool, f32)> = Vec::new();
    for pl in &polys {
        for &(x, y) in pl.iter() {
            k_shore.push(at(&b.k, b.w, x, y));
            gk_shore.push(grad_k(&b.k, b.w, x, y));
        }
        for w3 in pl.windows(3) {
            let (ax, ay) = (w3[1].0 - w3[0].0, w3[1].1 - w3[0].1);
            let (bx, by) = (w3[2].0 - w3[1].0, w3[2].1 - w3[1].1);
            let (la, lb) = ((ax * ax + ay * ay).sqrt(), (bx * bx + by * by).sqrt());
            if la <= 0.0 || lb <= 0.0 {
                continue;
            }
            let turn = (((ax * bx + ay * by) / (la * lb)).clamp(-1.0, 1.0)).acos().to_degrees();
            pairs.push((turn > 80.0, grad_k(&b.k, b.w, w3[1].0, w3[1].1)));
        }
    }
    eprintln!(
        "\n2. K ALONG THE SHORE ({} vertices)\n   K      p10 {:.3} | p50 {:.3} | p90 {:.3} | max {:.3}\n   \
         |grad K| p50 {:.5} | p90 {:.5} | p99 {:.5} | max {:.5}  (per cell)",
        k_shore.len(),
        pct(&mut k_shore.clone(), 0.10),
        pct(&mut k_shore.clone(), 0.50),
        pct(&mut k_shore.clone(), 0.90),
        pct(&mut k_shore.clone(), 1.0),
        pct(&mut gk_shore.clone(), 0.50),
        pct(&mut gk_shore.clone(), 0.90),
        pct(&mut gk_shore.clone(), 0.99),
        pct(&mut gk_shore.clone(), 1.0),
    );

    // THE CORRELATION: fringe rate in each |grad K| quintile. If fringes coincide with K
    // transitions, the rate must RISE with the gradient.
    let mut g: Vec<f32> = pairs.iter().map(|p| p.1).collect();
    let edges: Vec<f32> = (1..5).map(|i| pct(&mut g.clone(), i as f32 / 5.0)).collect();
    let bucket = |v: f32| edges.iter().position(|&e| v < e).unwrap_or(4);
    let mut cnt = [(0usize, 0usize); 5];
    for &(fringe, gv) in &pairs {
        let k = bucket(gv);
        cnt[k].0 += 1;
        if fringe {
            cnt[k].1 += 1;
        }
    }
    eprintln!(
        "\n   FRINGE RATE BY |grad K| QUINTILE — rises with the gradient ⇒ the fringes ARE the\n   \
         K transitions; flat ⇒ they are not."
    );
    eprintln!("   {:<10} {:>12} {:>10} {:>12}", "quintile", "|grad K| <", "vertices", "turns>80°");
    for i in 0..5 {
        let hi = if i < 4 { format!("{:.5}", edges[i]) } else { "max".into() };
        eprintln!(
            "   Q{:<9} {:>12} {:>10} {:>11.1}%",
            i + 1,
            hi,
            cnt[i].0,
            100.0 * cnt[i].1 as f32 / cnt[i].0.max(1) as f32
        );
    }

    // ── 3. THE 2×2, IN ABSOLUTE TERMS. A RATE hides this: the fringe rate is ~20 % whatever
    //    the closures do, but the coastline itself gets 26 % longer and gains 39 % more rings,
    //    so the ABSOLUTE number of fringes rises. "More fringed" is a quantity, not a rate —
    //    method rule 2, and I had reached for the rate first.
    struct Shape {
        rings: usize,
        turns: usize,
        fringes: usize,
        len_km: f32,
        tiny_rings: usize,
        main_len_km: f32,
    }
    let shape = |ps: &Vec<Vec<(f32, f32)>>| -> Shape {
        let (mut turns, mut fringes) = (0usize, 0usize);
        let (mut len, mut tiny, mut main) = (0.0f32, 0usize, 0.0f32);
        for pl in ps {
            let mut rl = 0.0f32;
            for w2 in pl.windows(2) {
                rl += ((w2[1].0 - w2[0].0).powi(2) + (w2[1].1 - w2[0].1).powi(2)).sqrt();
            }
            len += rl;
            main = main.max(rl);
            // A "tiny ring" is an islet or a digitation pinched off: under 20 cells of outline.
            if pl.len() < 20 {
                tiny += 1;
            }
            for w3 in pl.windows(3) {
                let (ax, ay) = (w3[1].0 - w3[0].0, w3[1].1 - w3[0].1);
                let (bx, by) = (w3[2].0 - w3[1].0, w3[2].1 - w3[1].1);
                let (la, lb) = ((ax * ax + ay * ay).sqrt(), (bx * bx + by * by).sqrt());
                if la <= 0.0 || lb <= 0.0 {
                    continue;
                }
                let t = (((ax * bx + ay * by) / (la * lb)).clamp(-1.0, 1.0)).acos().to_degrees();
                turns += 1;
                if t > 80.0 {
                    fringes += 1;
                }
            }
        }
        Shape {
            rings: ps.len(),
            turns,
            fringes,
            len_km: len * km_per_cell,
            tiny_rings: tiny,
            main_len_km: main * km_per_cell,
        }
    };

    eprintln!(
        "\n3. THE 2x2, IN ABSOLUTE TERMS (a rate hides it — the coastline itself changes length)"
    );
    eprintln!(
        "   {:<22} {:>7} {:>8} {:>9} {:>10} {:>10} {:>11}",
        "closures", "rings", "tiny", "fringes", "rate %", "coast km", "main ring km"
    );
    let mut base: Option<Shape> = None;
    for (label, litho, frac) in [
        ("all OFF", false, false),
        ("C-3 only", true, false),
        ("C-3b only", false, true),
        ("C-3 + C-3b (shipped)", true, true),
    ] {
        let bb = build(target, litho, frac);
        let sh = shape(&marching_squares(&bb.height, SEA));
        eprintln!(
            "   {:<22} {:>7} {:>8} {:>9} {:>9.2}% {:>10.0} {:>11.0}",
            label,
            sh.rings,
            sh.tiny_rings,
            sh.fringes,
            100.0 * sh.fringes as f32 / sh.turns.max(1) as f32,
            sh.len_km,
            sh.main_len_km
        );
        if label == "all OFF" {
            base = Some(sh);
        } else if let Some(b0) = &base {
            eprintln!(
                "   {:<22} {:>+7} {:>+8} {:>+9} {:>10} {:>+10.0} {:>+11.0}   <- vs all OFF",
                "",
                sh.rings as i64 - b0.rings as i64,
                sh.tiny_rings as i64 - b0.tiny_rings as i64,
                sh.fringes as i64 - b0.fringes as i64,
                "",
                sh.len_km - b0.len_km,
                sh.main_len_km - b0.main_len_km
            );
        }
    }

    // ── 4. SCOPING THE PRE-EXISTING BASELINE — not fixing it. With every K closure off the
    //    coastline is still ~4700 rings of which 99 % are under 20 vertices: a SPECKLE, not a
    //    fringe. Is it the FBM? `amplitude_base = 0` removes the FBM's contribution entirely
    //    (the C-1 relief-budget cap means the base term is inert in production, but zero is
    //    zero), leaving the bilinearly upscaled coarse field plus the incision.
    let no_fbm = build_amp(target, false, false, 0.0);
    let sh0 = shape(&marching_squares(&no_fbm.height, SEA));
    if let Some(b0) = &base {
        eprintln!(
            "
4. SCOPING THE BASELINE (all closures off, and then WITHOUT the FBM)
                {:<24} {:>7} {:>8} {:>9} {:>10}
                {:<24} {:>7} {:>8} {:>9} {:>10.0}
                {:<24} {:>7} {:>8} {:>9} {:>10.0}
                {:<24} {:>+7} {:>+8} {:>+9} {:>+10.0}",
            "",
            "rings",
            "tiny",
            "fringes",
            "coast km",
            "all closures off",
            b0.rings,
            b0.tiny_rings,
            b0.fringes,
            b0.len_km,
            "+ FBM amplitude = 0",
            sh0.rings,
            sh0.tiny_rings,
            sh0.fringes,
            sh0.len_km,
            "  the FBM's share",
            b0.rings as i64 - sh0.rings as i64,
            b0.tiny_rings as i64 - sh0.tiny_rings as i64,
            b0.fringes as i64 - sh0.fringes as i64,
            b0.len_km - sh0.len_km
        );
    }
}

#[test]
#[ignore]
fn coastal_fringes_k_2048() {
    run(2048);
}

#[test]
#[ignore]
fn coastal_fringes_k_8192() {
    run(8192);
}
