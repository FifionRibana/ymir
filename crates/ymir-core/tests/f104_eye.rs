//! ADR 0001 Finding 104 block B — seed 3's tile, because R8 has been demoted to reported-only and
//! the eye is now the gate. **No production change.**
//!
//! ⚠️ Seed 3 is the seed that decides: Finding 103 measured its R8 at −21.0 % (×0.4), −24.8 %
//! (×0.5) and −31.2 % (×0.6) — the only seed that does not clear a −30 % gate below ×0.6. With R8
//! reported rather than gated, **the question is whether ×0.5's field reads combed or dendritic**,
//! and that cannot be answered by the number that was just demoted.
//!
//! The tile is chosen by **max Δ(R8)** between the delivered field and ×0.6 (Finding 97's
//! correction: selecting by max R8 finds tiles anisotropic in every field and is blind to the
//! difference).
//!
//! Run: cargo test -p ymir-core --release --test f104_eye -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, aniso, build_field_seed, pct, sorted, tile};
use std::path::Path;
use std::time::Instant;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const OUT: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../docs/reports/c1_continental_buoyancy/f104_eye");

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
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

fn sigma_p50(f: &GridF32, land: &[bool], n2m: f32) -> f32 {
    let (w, h) = (f.width, f.height);
    let mut out = Vec::new();
    for y in 1..h - 1 {
        for x in 1..w - 1 {
            let k = y * w + x;
            if !land[k] {
                continue;
            }
            let (mut s, mut s2) = (0f64, 0f64);
            for dy in 0..3 {
                for dx in 0..3 {
                    let v = (f.data[(y + dy - 1) * w + x + dx - 1] * n2m) as f64;
                    s += v;
                    s2 += v * v;
                }
            }
            out.push(((s2 / 9.0 - (s / 9.0) * (s / 9.0)).max(0.0)).sqrt() as f32);
        }
    }
    pct(&sorted(out), 0.50)
}

fn save(g: &GridF32, name: &str) {
    std::fs::create_dir_all(Path::new(OUT)).expect("report dir");
    g.save_png_u8(&Path::new(OUT).join(name)).unwrap_or_else(|e| panic!("save {name}: {e}"));
    eprintln!("   wrote {name}");
}

#[test]
#[ignore]
fn f104_eye() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    let seed = splitmix64(splitmix64(PSEED));
    eprintln!("\n==========  Finding 104 B . seed 3 ({seed}), the eye  ==========");

    let t0 = Instant::now();
    let delivered = build_field_seed(Knobs::passes(2), seed);
    let (w, h) = (delivered.width, delivered.height);
    let n = w * h;
    // seed 3's own Flint intercept, as Finding 102 measured it
    let d_del = c1_drainage_windowed(&delivered, None, &on, &ss, DOMAIN_KM);
    let mut v = Vec::new();
    for k in (0..n).step_by(7) {
        if delivered.data[k] <= SEA {
            continue;
        }
        let a = d_del.flow.accumulation.data[k] * cell_km2;
        if a < RELIEF_V1_A_C_KM2 {
            continue;
        }
        let dir = d_del.flow.direction[k];
        if dir == DIR_NONE {
            continue;
        }
        let nx = ((k % w) as i32 + D8_DX[dir as usize]).rem_euclid(w as i32) as usize;
        let ny = ((k / w) as i32 + D8_DY[dir as usize]).rem_euclid(h as i32) as usize;
        let diag = D8_DX[dir as usize] != 0 && D8_DY[dir as usize] != 0;
        let dx_m = if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M };
        let s = (delivered.data[k] - delivered.data[ny * w + nx]).max(0.0) * n2m / dx_m;
        if s > 0.0 {
            v.push(s * a.sqrt());
        }
    }
    let ks = pct(&sorted(v), 0.50);
    eprintln!("   seed 3 k_s = **{ks:.4}** (Finding 102: 0.0277)");

    let mk = |m: f32| {
        build_field_seed(
            Knobs { slope_floor_uk: Some(ks * m), depression_floor: true, ..Knobs::passes(2) },
            seed,
        )
    };
    let f4 = mk(0.4);
    let f5 = mk(0.5);
    let f6 = mk(0.6);
    let ld: Vec<bool> = (0..n).map(|k| delivered.data[k] > SEA).collect();
    let l6: Vec<bool> = (0..n).map(|k| f6.data[k] > SEA).collect();

    // the tile of maximum delta(R8) between delivered and x0.6
    let (mut bx, mut by, mut bd) = (0usize, 0usize, -1f32);
    for ty in (0..h).step_by(1024) {
        for tx in (0..w).step_by(1024) {
            let (ga, gla) = tile(&delivered, &ld, tx, ty, 1024);
            let (gb, glb) = tile(&f6, &l6, tx, ty, 1024);
            let (ra, rb) = (aniso(&ga, &gla, 16), aniso(&gb, &glb, 16));
            if ra.windows < 64 || rb.windows < 64 {
                continue;
            }
            if (ra.r8 - rb.r8).abs() > bd {
                bd = (ra.r8 - rb.r8).abs();
                bx = tx;
                by = ty;
            }
        }
    }
    eprintln!("   tile of max delta(R8): **({bx}, {by})**, delta {bd:.4}");

    for (label, f) in [("delivered", &delivered), ("x0.4", &f4), ("x0.5", &f5), ("x0.6", &f6)] {
        let (g, gl) = tile(f, &ld, bx, by, 1024);
        eprintln!(
            "   {label:<10} tile R8 **{:.4}** . sigma **{:.2} m**",
            aniso(&g, &gl, 16).r8,
            sigma_p50(&g, &gl, n2m)
        );
        save(&hillshade(&g, n2m), &format!("seed3_{label}_hillshade.png"));
    }
    eprintln!("\n==========  end . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
