//! C-3 EXPORT comparison — the verdict is the export, not a reconstruction (ADR
//! 0001 / method memory). Decodes two shipped `height.u16` rasters (each with ITS
//! OWN min/max metre scale from its manifest), rebuilds the lithology CLASS map for
//! the export window, and reports the per-class Δ so the "hard basement unchanged,
//! only soft/volcaniclastic move" claim is checked on the actual product.
//!
//! NEW = the lithology-ON export, REF = the C-2 baseline (volcanism on, lithology
//! off). Both: seed 10481999410520546993, 8192², window offset [0.09375, 0.578125],
//! domain 400 km, humid climate.
//!
//! Run: cargo test -p ymir-core --test c3_export_compare -- --ignored --nocapture

use std::path::Path;

use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::lithology;
use ymir_core::tectonics_c1::closures::volcanism::{VolcanismConfig, place_edifices};
use ymir_core::tectonics_c1::init_r7::{Phase2InitParams, init_c1_state_phase_2_r7};
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig, run_with_closures};
use ymir_core::tectonics_v2::boundaries::plate_type::PlateType;

const PSEED: u64 = 10481999410520546993;
const DOMAIN_KM: f32 = 400.0;
const TARGET: usize = 8192;
const SO: [f64; 2] = [0.09375, 0.578125];

const NEW_DIR: &str = "../../exports/seed10481999410520546993_8192.ymir";
const REF_DIR: &str = "../../exports/seed10481999410520546993_8192.volcan.humidymir";

const HARD: u8 = 0;
const SOFT: u8 = 1;
const VOLC: u8 = 2;

/// Decode a `height.u16` (little-endian, linear `min_m..max_m` over the u16 range)
/// to metres, reading `min_m`/`max_m` from the sibling `manifest.json`.
fn load_height_m(dir: &str) -> Vec<f32> {
    let manifest: serde_json::Value =
        serde_json::from_slice(&std::fs::read(Path::new(dir).join("manifest.json")).unwrap())
            .unwrap();
    let layer =
        manifest["layers"].as_array().unwrap().iter().find(|l| l["id"] == "height").unwrap();
    let min_m = layer["min_m"].as_f64().unwrap() as f32;
    let max_m = layer["max_m"].as_f64().unwrap() as f32;
    let bytes = std::fs::read(Path::new(dir).join("height.u16")).unwrap();
    let span = max_m - min_m;
    bytes
        .chunks_exact(2)
        .map(|b| {
            let v = u16::from_le_bytes([b[0], b[1]]) as f32 / 65535.0;
            min_m + v * span
        })
        .collect()
}

/// Rebuild the lithology class map (0 hard / 1 rift-soft / 2 volcaniclastic) for the
/// export window — identical construction to the production K field.
fn class_map() -> Vec<u8> {
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
    let volc = VolcanismConfig { enabled: true, domain_km: DOMAIN_KM, ..Default::default() };
    let edifices = place_edifices(&state, &kin, &seed, DOMAIN_KM, &volc);

    let (w, h) = (TARGET, TARGET);
    let km_per_hd = DOMAIN_KM / TARGET as f32;
    // Coarse rift mask → upscale → threshold.
    let (nx, ny) = (state.plate_id.nx(), state.plate_id.ny());
    let mut rift = GridF32::new(nx, ny, 0.0);
    for j in 0..ny {
        for i in 0..nx {
            if matches!(state.plate_type.get(i, j), PlateType::Continental)
                && (state.age.get(i, j) as f32) < 1.0
            {
                rift.set(i, j, 1.0);
            }
        }
    }
    let up = lithology::upscale_k_to_hd(&rift, w, h, SO, 1.0);
    let mut cls = vec![HARD; w * h];
    for (k, &v) in up.iter().enumerate() {
        if v > 0.5 {
            cls[k] = SOFT;
        }
    }
    let (so32, ss32) = ([SO[0] as f32, SO[1] as f32], 1.0f32);
    for e in &edifices {
        let fx = (e.center_uv.0 - so32[0]).rem_euclid(1.0) / ss32;
        let fy = (e.center_uv.1 - so32[1]).rem_euclid(1.0) / ss32;
        if fx >= 1.0 || fy >= 1.0 {
            continue;
        }
        let (cx, cy) = (fx * w as f32, fy * h as f32);
        let rb = (e.basal_diameter_km * 0.5) / km_per_hd;
        if rb < 1.0 {
            continue;
        }
        let (i0, i1) =
            ((cx - rb).floor().max(0.0) as usize, ((cx + rb).ceil() as usize).min(w - 1));
        let (j0, j1) =
            ((cy - rb).floor().max(0.0) as usize, ((cy + rb).ceil() as usize).min(h - 1));
        for j in j0..=j1 {
            for i in i0..=i1 {
                let d = ((i as f32 + 0.5 - cx).powi(2) + (j as f32 + 0.5 - cy).powi(2)).sqrt();
                if d <= rb {
                    cls[j * w + i] = VOLC;
                }
            }
        }
    }
    cls
}

#[test]
#[ignore]
fn c3_export_compare() {
    if !Path::new(NEW_DIR).exists() || !Path::new(REF_DIR).exists() {
        eprintln!("export dirs missing — skip");
        return;
    }
    let new = load_height_m(NEW_DIR);
    let refh = load_height_m(REF_DIR);
    assert_eq!(new.len(), refh.len(), "raster size mismatch");
    let cls = class_map();
    assert_eq!(cls.len(), new.len(), "class map size mismatch");

    let names = ["HARD", "SOFT(rift)", "VOLC"];
    eprintln!("\n=== C-3 EXPORT COMPARE (NEW − REF, metres) — the verdict is the export ===");
    eprintln!("NEW = {NEW_DIR}");
    eprintln!("REF = {REF_DIR}");
    eprintln!(
        "{:>12} | {:>10} | {:>12} | {:>12} | {:>12} | {:>10}",
        "class", "cells", "mean|Δ| m", "mean Δ m", "max|Δ| m", "changed‰"
    );
    for class in [HARD, SOFT, VOLC] {
        let (mut n, mut sum_abs, mut sum, mut mx, mut changed) =
            (0u64, 0.0f64, 0.0f64, 0.0f32, 0u64);
        for k in 0..new.len() {
            if cls[k] != class {
                continue;
            }
            // land-only (either side above sea) — lithology acts on land incision.
            if new[k] <= 0.0 && refh[k] <= 0.0 {
                continue;
            }
            n += 1;
            let d = new[k] - refh[k];
            sum += d as f64;
            sum_abs += d.abs() as f64;
            mx = mx.max(d.abs());
            if d.abs() > 1.0 {
                changed += 1;
            }
        }
        let nn = n.max(1) as f64;
        eprintln!(
            "{:>12} | {:>10} | {:>12.2} | {:>12.2} | {:>12.1} | {:>10.1}",
            names[class as usize],
            n,
            sum_abs / nn,
            sum / nn,
            mx,
            1000.0 * changed as f32 / n.max(1) as f32,
        );
    }
    eprintln!(
        "(land cells only. If lithology is ON in NEW: HARD mean|Δ| ≪ SOFT/VOLC, and HARD changed‰ small\n = reference intact, contrast localised to the soft classes. If all three are ~0, NEW has lithology OFF.)"
    );

    write_delta_png(
        &new,
        &refh,
        &cls,
        "../../docs/reports/c1_continental_buoyancy/closure_morphology/c3_delta_map.png",
    );
}

/// Δ map, downsampled ×`STRIDE` for a viewable PNG. Background = lithology class
/// (hard slate, rift teal, volcaniclastic violet); RED glow where NEW eroded below
/// REF (Δ<0, intensity ∝ |Δ|/200 m); faint BLUE where it built up. Ocean stays dark.
fn write_delta_png(new: &[f32], refh: &[f32], cls: &[u8], out: &str) {
    const STRIDE: usize = 4;
    let (w, h) = (TARGET / STRIDE, TARGET / STRIDE);
    let mut img = image::RgbImage::new(w as u32, h as u32);
    for oy in 0..h {
        for ox in 0..w {
            // Aggregate the STRIDE² block: mean Δ over land, strongest class present.
            let (mut sd, mut nl, mut best) = (0.0f64, 0u32, HARD);
            let mut any_land = false;
            for dy in 0..STRIDE {
                for dx in 0..STRIDE {
                    let k = (oy * STRIDE + dy) * TARGET + (ox * STRIDE + dx);
                    if cls[k] > best {
                        best = cls[k];
                    }
                    if new[k] > 0.0 || refh[k] > 0.0 {
                        any_land = true;
                        sd += (new[k] - refh[k]) as f64;
                        nl += 1;
                    }
                }
            }
            let rgb = if !any_land {
                [12u8, 20, 40] // ocean
            } else {
                let base = match best {
                    SOFT => [40i32, 110, 110], // rift teal
                    VOLC => [120, 60, 130],    // volcaniclastic violet
                    _ => [80, 82, 88],         // hard slate
                };
                let d = if nl > 0 { (sd / nl as f64) as f32 } else { 0.0 };
                let t = (d.abs() / 200.0).clamp(0.0, 1.0); // 200 m = full glow
                let mix = |c: i32, target: i32| (c as f32 * (1.0 - t) + target as f32 * t) as i32;
                let (tr, tg, tb) = if d < 0.0 { (255, 40, 30) } else { (60, 120, 255) };
                [
                    mix(base[0], tr).clamp(0, 255) as u8,
                    mix(base[1], tg).clamp(0, 255) as u8,
                    mix(base[2], tb).clamp(0, 255) as u8,
                ]
            };
            img.put_pixel(ox as u32, oy as u32, image::Rgb(rgb));
        }
    }
    img.save(out).unwrap();
    eprintln!("\nΔ map written: {out}  ({w}×{h}, red = eroded in NEW, class-tinted background)");
}
