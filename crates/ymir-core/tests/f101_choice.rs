//! ADR 0001 Finding 101 — the panels the AUTHOR needs to choose between ×0.4 and ×0.5.
//! **No production change. Nothing promoted.**
//!
//! ⚠️ Finding 101's own images cannot decide this and said so: on the (2048, 5120) tile the two
//! candidates differ by **0.0003 of R8** and **0.07 m of σ**. The continent water view DID show a
//! difference — and it was rendered for the delivered field and ×0.4 only, not ×0.5.
//!
//! The two candidates' largest relative difference is the **median cut per cell, 13.1 m against
//! 6.9 m (×1.9)**, so the panel that can show it is a **closer** zoom on the tile where they
//! actually differ — not a fixed tile at 1024².
//!
//! Run: cargo test -p ymir-core --release --test f101_choice -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, SEA, build_field, pct, sorted};
use std::path::Path;
use std::time::Instant;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::lakes::connectivity::water_class;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const KS: f32 = 0.0451;
const OUT: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../docs/reports/c1_continental_buoyancy/f101_choice");
/// 512² = 25 km, twice the magnification of Finding 101's panels.
const SIDE: usize = 512;

fn save(g: &GridF32, name: &str) {
    std::fs::create_dir_all(Path::new(OUT)).expect("report dir");
    g.save_png_u8(&Path::new(OUT).join(name)).unwrap_or_else(|e| panic!("save {name}: {e}"));
    eprintln!("   wrote {name}");
}

fn hillshade(t: &GridF32, n2m: f32) -> GridF32 {
    let (w, h) = (t.width, t.height);
    let (az, alt) = (315f32.to_radians(), 45f32.to_radians());
    let (lx, ly, lz) = (alt.cos() * az.sin(), -alt.cos() * az.cos(), alt.sin());
    let mut g = GridF32::new(w, h, 0.0);
    for y in 0..h {
        for x in 0..w {
            let (gx, gy) = t.gradient_at(x, y);
            let (dx, dy) = (gx * n2m / CELL_M, gy * n2m / CELL_M);
            let nrm = (dx * dx + dy * dy + 1.0).sqrt();
            g.data[y * w + x] = ((-dx * lx - dy * ly + lz) / nrm).clamp(0.0, 1.0);
        }
    }
    g
}

fn crop(f: &GridF32, bx: usize, by: usize, side: usize) -> GridF32 {
    let mut g = GridF32::new(side, side, 0.0);
    for y in 0..side {
        for x in 0..side {
            g.data[y * side + x] = f.data[(by + y) * f.width + bx + x];
        }
    }
    g
}

#[test]
#[ignore]
fn f101_choice() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    eprintln!("\n==========  Finding 101 . the panels for the author's choice  ==========");

    let t0 = Instant::now();
    let delivered = build_field(Knobs::passes(2));
    let (w, h) = (delivered.width, delivered.height);
    let n = w * h;
    let kn =
        |m: f32| Knobs { slope_floor_uk: Some(KS * m), depression_floor: true, ..Knobs::passes(2) };
    let f4 = build_field(kn(0.4));
    let f5 = build_field(kn(0.5));
    eprintln!("   three fields in {:.1} s", t0.elapsed().as_secs_f64());

    // ── 1. the missing continent panel, and the two that exist, on ONE ramp ──
    for (label, f) in [("delivered", &delivered), ("x0.4", &f4), ("x0.5", &f5)] {
        let wc = water_class(f, SEA);
        let d = c1_drainage_windowed(f, None, &on, &ss, DOMAIN_KM);
        let mut g = GridF32::new(1024, 1024, 0.0);
        let (mut water, mut land) = (0u64, 0u64);
        for y in 0..1024 {
            for x in 0..1024 {
                let (mut wq, mut lq) = (0u32, 0u32);
                for dy in 0..8 {
                    for dx in 0..8 {
                        let k = (y * 8 + dy) * w + x * 8 + dx;
                        if f.data[k] <= SEA {
                            continue;
                        }
                        lq += 1;
                        if d.lake_map[k] != 0 || wc[k] == 2 {
                            wq += 1;
                        }
                    }
                }
                water += wq as u64;
                land += lq as u64;
                g.data[y * 1024 + x] =
                    if lq == 0 { 0.0 } else { 0.25 + 0.75 * (1.0 - wq as f32 / lq as f32) };
            }
        }
        eprintln!(
            "   continent {label:<10} water {:.2} % of land",
            100.0 * water as f32 / land as f32
        );
        save(&g, &format!("continent_{label}.png"));
    }

    // ── 2. where do x0.4 and x0.5 actually differ? ──────────────────────────
    let mut best = (0usize, 0usize, -1f32);
    for by in (0..h - SIDE).step_by(SIDE) {
        for bx in (0..w - SIDE).step_by(SIDE) {
            let (mut s, mut c) = (0f64, 0usize);
            for y in (by..by + SIDE).step_by(4) {
                for x in (bx..bx + SIDE).step_by(4) {
                    let k = y * w + x;
                    if f4.data[k] <= SEA || f5.data[k] <= SEA {
                        continue;
                    }
                    s += ((f4.data[k] - f5.data[k]) * n2m).abs() as f64;
                    c += 1;
                }
            }
            // require a mostly-land window, else a coastal sliver wins on noise
            if c < (SIDE / 4) * (SIDE / 4) * 8 / 10 {
                continue;
            }
            let m = (s / c as f64) as f32;
            if m > best.2 {
                best = (bx, by, m);
            }
        }
    }
    let (bx, by, dm) = best;
    let dif: Vec<f32> = (0..n)
        .filter(|&k| f4.data[k] > SEA && f5.data[k] > SEA)
        .map(|k| ((f4.data[k] - f5.data[k]) * n2m).abs())
        .collect();
    let ds = sorted(dif);
    eprintln!(
        "\n   |x0.4 - x0.5| over land: p50 **{:.1} m** p90 **{:.1} m** max {:.0} m . the 512² \
         window where they differ MOST: **({bx}, {by})**, mean |delta| **{dm:.1} m**",
        pct(&ds, 0.50),
        pct(&ds, 0.90),
        ds.last().copied().unwrap_or(0.0)
    );

    // ── 3. the zoom, at 512² = 25 km, twice Finding 101's magnification ─────
    for (label, f) in [("delivered", &delivered), ("x0.4", &f4), ("x0.5", &f5)] {
        let c = crop(f, bx, by, SIDE);
        let land: Vec<bool> = c.data.iter().map(|&v| v > SEA).collect();
        let cut: Vec<f32> =
            (0..SIDE * SIDE).filter(|&k| land[k]).map(|k| c.data[k] * n2m).collect();
        eprintln!(
            "   zoom {label:<10} altitude p50 {:.0} m over the window",
            pct(&sorted(cut), 0.50)
        );
        save(&hillshade(&c, n2m), &format!("zoom_{label}.png"));
    }

    // ── 4. the signed difference at that window: WHERE x0.4 cut deeper ──────
    let mut g = GridF32::new(SIDE, SIDE, 0.5);
    for y in 0..SIDE {
        for x in 0..SIDE {
            let k = (by + y) * w + bx + x;
            let d = (f4.data[k] - f5.data[k]) * n2m; // negative = x0.4 is LOWER = cut deeper
            g.data[y * SIDE + x] = (0.5 + d / 200.0).clamp(0.0, 1.0);
        }
    }
    save(&g, "zoom_difference_x0.4_minus_x0.5.png");
    eprintln!(
        "   (mid-grey = equal, BLACK = x0.4 is up to 100 m LOWER than x0.5 = cut deeper, WHITE = \
         higher)"
    );
    eprintln!("\n==========  end . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
