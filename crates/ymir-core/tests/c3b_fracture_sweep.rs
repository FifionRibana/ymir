//! C-3b inherited-structure SWEEP (closures roadmap §3b) — density only. The fracture
//! DENSITY field (proximity to plate contacts + sutures, causal) modulates erodibility
//! ISOTROPICALLY: `K = 1 + amplitude · density`. Intact cratonic interior (density → 0)
//! stays at the ×1 reference and RETAINS its relief; fractured belts near contacts erode
//! more. Runs the WHOLE production chain (`upscale_from_c1_with_progress`, export recipe
//! relief-v3), fracture OFF then ON at swept amplitudes, and reports the relief contrast
//! between the two zones + the C-1 invariant.
//!
//! The amplitude is a MEASUREMENT, not a prediction (Molnar's ~1–2 orders is the outer
//! bound). Hard basement (density 0) stays ×1 → global-slowdown nil by construction (the
//! C-3 design survives). C-1 must survive (pits must NOT rise as they did in the failed
//! anisotropic test).
//!
//! Run: cargo test -p ymir-core --test c3b_fracture_sweep c3b_sweep -- --ignored --nocapture

use ymir_core::erosion::stream_power::StreamPowerConfig;
use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::fracture::{FractureConfig, build_hd_density_k};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::VolcanismConfig;
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::upscale_from_c1_with_progress;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::terrain::upscale::FbmUpscaleConfig;

const PSEED: u64 = 10481999410520546993;
const DOMAIN_KM: f32 = 400.0;
const FULL_RANGE_M: f32 = 11302.0;

fn count_closed_depressions(field: &GridF32) -> usize {
    use ymir_core::terrain::flow::{FlowConfig, compute_flow};
    let (w, h) = (field.width, field.height);
    let n = w * h;
    let flow =
        compute_flow(field, &FlowConfig { sea_level: 0.5, flat_perturbation: None, dinf: false });
    let thr = 0.1f32 / FULL_RANGE_M;
    let mut seen = vec![false; n];
    let mut count = 0usize;
    let mut stack = Vec::new();
    for s in 0..n {
        if seen[s] || flow.filled.data[s] - field.data[s] <= thr {
            continue;
        }
        count += 1;
        seen[s] = true;
        stack.push(s);
        while let Some(k) = stack.pop() {
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for dy in -1..=1 {
                for dx in -1..=1 {
                    let (nx, ny) = (x + dx, y + dy);
                    if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                        let nk = ny as usize * w + nx as usize;
                        if !seen[nk] && flow.filled.data[nk] - field.data[nk] > thr {
                            seen[nk] = true;
                            stack.push(nk);
                        }
                    }
                }
            }
        }
    }
    count
}

/// Median local relief (physical ±5-cell window, m) over land cells whose fracture
/// density is in `[lo, hi)` — lets the sweep separate the intact-craton zone (density
/// low, the reference) from the fractured belt (density high).
fn zone_relief_m(field: &GridF32, density: &[f32], lo: f32, hi: f32) -> (f32, u64) {
    let (w, h) = (field.width, field.height);
    let r = 5i32;
    let mut rel = Vec::new();
    for y in (r as usize..h - r as usize).step_by(3) {
        for x in (r as usize..w - r as usize).step_by(3) {
            let k = y * w + x;
            if field.data[k] <= 0.5 || density[k] < lo || density[k] >= hi {
                continue;
            }
            let (mut mn, mut mx) = (f32::MAX, f32::MIN);
            for dy in -r..=r {
                for dx in -r..=r {
                    let v = field.data[(y as i32 + dy) as usize * w + (x as i32 + dx) as usize];
                    mn = mn.min(v);
                    mx = mx.max(v);
                }
            }
            rel.push((mx - mn) * FULL_RANGE_M);
        }
    }
    let n = rel.len() as u64;
    rel.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (rel.get(rel.len() / 2).copied().unwrap_or(0.0), n)
}

fn run_sweep(target: usize) {
    let ss = SteinSteinParams::default();
    let run = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let mut state = init_c1_state_phase_2_r7(64, PSEED, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, &run, &C1Closures::default(), |_, _| {});
    let seed = WorldSeed::new(PSEED);
    let volc = VolcanismConfig::default();

    let km_per_cell = DOMAIN_KM / target as f32;
    let cell_km2 = km_per_cell * km_per_cell;
    let mut base_cfg = FbmUpscaleConfig::c1_hd_production(target);
    base_cfg.sample_origin = [0.09375, 0.578125];
    base_cfg.sample_size = 1.0;
    base_cfg.erosion = None;
    let mut sp = StreamPowerConfig::relief_v3(cell_km2, ss.depth_scale_m as f32);
    sp.mfd_exponent = Some(2.0);
    sp.iterations = 2;
    base_cfg.stream_power = Some(sp);

    // The density field (amplitude 1 → density = K − 1), for zone classification and
    // coverage. Same geometry for every sweep row.
    let dcfg = FractureConfig {
        enabled: true,
        amplitude: 1.0,
        decay_km: 25.0,
        domain_km: DOMAIN_KM,
        ..Default::default()
    };
    let dens_k = build_hd_density_k(
        &state,
        &kin,
        &dcfg,
        None,
        target,
        target,
        base_cfg.sample_origin,
        base_cfg.sample_size,
    );
    let density: Vec<f32> = dens_k.iter().map(|k| k - 1.0).collect();
    let belt = density.iter().filter(|&&d| d >= 0.5).count();
    let n = density.len();
    eprintln!("\n================  C-3b DENSITY SWEEP  {target}²  ================");
    eprintln!(
        "coverage: belt(density≥0.5) {:.1}%  | transition(0.2–0.5) {:.1}%  | craton(<0.2) {:.1}%",
        100.0 * belt as f32 / n as f32,
        100.0 * density.iter().filter(|&&d| (0.2..0.5).contains(&d)).count() as f32 / n as f32,
        100.0 * density.iter().filter(|&&d| d < 0.2).count() as f32 / n as f32,
    );
    eprintln!(
        "{:>10} | {:>18} | {:>18} | {:>10} | {:>8}",
        "amp ×", "CRATON relief m", "BELT relief m", "contrast", "pits"
    );

    for &amp in &[0.0f32, 2.0, 4.0, 8.0, 16.0] {
        let mut cfg = base_cfg.clone();
        cfg.fracture = FractureConfig {
            enabled: amp > 0.0,
            amplitude: amp,
            decay_km: 25.0,
            domain_km: DOMAIN_KM,
            ..Default::default()
        };
        let (res, _craters) = upscale_from_c1_with_progress(
            &state,
            &run.iso_config,
            &ss,
            &seed,
            &cfg,
            &[],
            &volc,
            Some(&kin),
            &mut |_| {},
            &|| false,
        );
        let field = &res.heightmap;
        let (craton, _nc) = zone_relief_m(field, &density, 0.0, 0.2);
        let (beltr, _nb) = zone_relief_m(field, &density, 0.5, 2.0);
        eprintln!(
            "{:>10} | {:>18.0} | {:>18.0} | {:>10.2} | {:>8}",
            if amp > 0.0 { format!("×{amp:.0}") } else { "OFF".into() },
            craton,
            beltr,
            if craton > 0.0 { beltr / craton } else { 0.0 },
            count_closed_depressions(field),
        );
    }
    eprintln!(
        "(CRATON = intact interior, the ×1 reference → relief should stay ~constant across the sweep;\n BELT = fractured, near contacts → relief should DROP with amplitude. pits must not rise = C-1 survives.)"
    );
}

#[test]
#[ignore]
fn c3b_sweep() {
    run_sweep(2048);
}

#[test]
#[ignore]
fn c3b_sweep_8192() {
    run_sweep(8192);
}

/// Δ map: fracture ×4 − OFF at 8192², downsampled ×4 for a viewable PNG. Background =
/// fracture-density zone (craton slate / transition olive / belt warm); RED where the
/// fractured belt eroded DOWN, faint blue where it built up. Ocean dark. Shows WHERE the
/// closure acts (the orogenic belts) leaving the intact craton untouched.
#[test]
#[ignore]
fn c3b_delta_map() {
    let target = 8192usize;
    let ss = SteinSteinParams::default();
    let run = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / 64.0,
        dy: 1.0 / 64.0,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let mut state = init_c1_state_phase_2_r7(64, PSEED, &Phase2InitParams::default());
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, &run, &C1Closures::default(), |_, _| {});
    let seed = WorldSeed::new(PSEED);
    let volc = VolcanismConfig::default();
    let km_per_cell = DOMAIN_KM / target as f32;
    let cell_km2 = km_per_cell * km_per_cell;
    let mut base = FbmUpscaleConfig::c1_hd_production(target);
    base.sample_origin = [0.09375, 0.578125];
    base.sample_size = 1.0;
    base.erosion = None;
    let mut sp = StreamPowerConfig::relief_v3(cell_km2, ss.depth_scale_m as f32);
    sp.mfd_exponent = Some(2.0);
    sp.iterations = 2;
    base.stream_power = Some(sp);

    let dcfg = FractureConfig {
        enabled: true,
        amplitude: 1.0,
        decay_km: 25.0,
        domain_km: DOMAIN_KM,
        ..Default::default()
    };
    let dens_k = build_hd_density_k(
        &state,
        &kin,
        &dcfg,
        None,
        target,
        target,
        base.sample_origin,
        base.sample_size,
    );

    let mk = |amp: f32| {
        let mut cfg = base.clone();
        cfg.fracture = FractureConfig {
            enabled: amp > 0.0,
            amplitude: amp,
            decay_km: 25.0,
            domain_km: DOMAIN_KM,
            ..Default::default()
        };
        upscale_from_c1_with_progress(
            &state,
            &run.iso_config,
            &ss,
            &seed,
            &cfg,
            &[],
            &volc,
            Some(&kin),
            &mut |_| {},
            &|| false,
        )
        .0
        .heightmap
    };
    let off = mk(0.0);
    let on = mk(4.0);

    const STRIDE: usize = 4;
    let (w, h) = (target / STRIDE, target / STRIDE);
    let mut img = image::RgbImage::new(w as u32, h as u32);
    for oy in 0..h {
        for ox in 0..w {
            let (mut dsum, mut nl, mut dens) = (0.0f64, 0u32, 0.0f32);
            let mut any_land = false;
            for dy in 0..STRIDE {
                for dx in 0..STRIDE {
                    let k = (oy * STRIDE + dy) * target + (ox * STRIDE + dx);
                    dens = dens.max(dens_k[k] - 1.0);
                    if on.data[k] > 0.5 || off.data[k] > 0.5 {
                        any_land = true;
                        dsum += (on.data[k] - off.data[k]) as f64;
                        nl += 1;
                    }
                }
            }
            let rgb = if !any_land {
                [12u8, 20, 40]
            } else {
                // Zone tint: craton slate / transition olive / belt warm.
                let base = if dens >= 0.5 {
                    [120i32, 90, 60]
                } else if dens >= 0.2 {
                    [95, 95, 60]
                } else {
                    [80, 82, 88]
                };
                let d = if nl > 0 { (dsum / nl as f64) as f32 * FULL_RANGE_M } else { 0.0 };
                let t = (d.abs() / 150.0).clamp(0.0, 1.0);
                let (tr, tg, tb) = if d < 0.0 { (255, 45, 30) } else { (60, 120, 255) };
                let mix = |c: i32, tg: i32| (c as f32 * (1.0 - t) + tg as f32 * t) as i32;
                [
                    mix(base[0], tr).clamp(0, 255) as u8,
                    mix(base[1], tg).clamp(0, 255) as u8,
                    mix(base[2], tb).clamp(0, 255) as u8,
                ]
            };
            img.put_pixel(ox as u32, oy as u32, image::Rgb(rgb));
        }
    }
    let out = "../../docs/reports/c1_continental_buoyancy/closure_morphology/c3b_delta_map.png";
    img.save(out).unwrap();
    eprintln!("\nΔ map written: {out}  ({w}×{h}; red = belt eroded by fracture, zone-tinted)");
}
