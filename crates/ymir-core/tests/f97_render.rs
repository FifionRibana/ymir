//! ADR 0001 Finding 97 — re-render the Finding 96 tile with a CALIBRATED renderer.
//! **No production change.**
//!
//! ⛔ Finding 96 refuted its own best candidate on an image, and Finding 97's block A cannot
//! reproduce the defect that image showed: the χ field's D8-lattice harmonic `R8` is **lower**
//! than the delivered field's, not higher. Before the instrument is called blind, the IMAGE has
//! to be audited — because it is an instrument too, and Finding 96 never calibrated it.
//!
//! The suspicion, stated before the render: Finding 96 exported `save_png_u8` over a ramp of
//! **1 m .. 3 250 m**, i.e. **~13.8 m of altitude per grey level**. Where the χ floor RAISED the
//! terrain and flattened it, adjacent cells fall in the same code and the image shows
//! **quantisation bands along the iso-contours** — long parallel lines on a gentle slope, which
//! is exactly what I described as "striations following the D8 axes".
//!
//! Three renders of the same tile settle it:
//!  * `_u8_global` — Finding 96's own renderer, reproduced, as the control;
//!  * `_u8_local`  — the same, stretched to the TILE's own p1..p99 (~16× finer quantisation);
//!  * `_hillshade` — a shaded relief, which reads the SLOPE and cannot band on altitude at all.
//!
//! If the lines survive the local ramp and the hillshade, they are terrain. If they vanish, they
//! were the renderer, and Finding 96's visual refutation has to be withdrawn.
//!
//! Run: cargo test -p ymir-core --release --test f97_render -- --ignored --nocapture

mod common;

use common::{Knobs, SEA, build_field, build_field_with_floor, pct, sorted};
use std::path::Path;
use std::time::Instant;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, breach_monotone, propagation_order};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const ORACLE_P50_M: f32 = 488.1;
const OUT: &str =
    concat!(env!("CARGO_MANIFEST_DIR"), "/../../docs/reports/c1_continental_buoyancy/f97_render");
/// The tile Finding 96's eye was shown.
const TX: usize = 5120;
const TY: usize = 3072;
const SIDE: usize = 1024;

fn save(g: &GridF32, name: &str) {
    std::fs::create_dir_all(Path::new(OUT)).expect("report dir");
    let p = Path::new(OUT).join(name);
    g.save_png_u8(&p).unwrap_or_else(|e| panic!("save {name}: {e}"));
    eprintln!("   wrote {}", p.display());
}

fn crop(f: &GridF32) -> GridF32 {
    let mut g = GridF32::new(SIDE, SIDE, 0.0);
    for y in 0..SIDE {
        for x in 0..SIDE {
            g.data[y * SIDE + x] = f.data[(TY + y) * f.width + TX + x];
        }
    }
    g
}

fn ramp(t: &GridF32, lo: f32, hi: f32) -> GridF32 {
    let mut g = GridF32::new(t.width, t.height, 0.0);
    for k in 0..t.data.len() {
        g.data[k] = if t.data[k] <= SEA {
            0.0
        } else {
            0.08 + 0.92 * ((t.data[k] - lo) / (hi - lo)).clamp(0.0, 1.0)
        };
    }
    g
}

/// Shaded relief: `cos` of the angle between the surface normal and a light at 315°, 45°.
/// It reads the SLOPE, so a flat step in altitude cannot band it — only real structure shows.
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
            let c = (-dx * lx - dy * ly + lz) / nrm;
            g.data[y * w + x] = c.clamp(0.0, 1.0);
        }
    }
    g
}

/// How many distinct u8 codes a render actually uses, and how many metres one code is worth.
fn codes(g: &GridF32) -> usize {
    let mut seen = [false; 256];
    for &v in &g.data {
        seen[((v.clamp(0.0, 1.0) * 255.0) as usize).min(255)] = true;
    }
    seen.iter().filter(|&&b| b).count()
}

#[test]
#[ignore]
fn f97_render() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = (DOMAIN_KM / 8192.0) * (DOMAIN_KM / 8192.0);
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    eprintln!("\n==========  Finding 97 . auditing the Finding 96 IMAGE  ==========");

    let t0 = Instant::now();
    let pre = build_field(Knobs::no_incision());
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let delivered = build_field(Knobs::passes(2));
    let d_pre = c1_drainage_windowed(&pre, None, &on, &ss, DOMAIN_KM);
    let br = breach_monotone(&pre, &d_pre.flow.filled, &d_pre.lake_map, SEA, w, h);
    let d_br = c1_drainage_windowed(&br, None, &on, &ss, DOMAIN_KM);
    let land_pre: Vec<bool> = (0..n).map(|k| br.data[k] > SEA).collect();
    let (order, un) = propagation_order(&d_br.flow.direction, |k| land_pre[k], w, h);
    assert_eq!(un, 0);
    let mut chi = vec![0.0f32; n];
    for &ku in order.iter().rev() {
        let k = ku as usize;
        if !land_pre[k] {
            continue;
        }
        let d = d_br.flow.direction[k];
        if d == DIR_NONE {
            continue;
        }
        let nx = ((k % w) as i32 + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
        let ny = ((k / w) as i32 + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
        let r = ny * w + nx;
        let diag = D8_DX[d as usize] != 0 && D8_DY[d as usize] != 0;
        let dx_m = if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M };
        let a = (d_br.flow.accumulation.data[k] * cell_km2).max(RELIEF_V1_A_C_KM2);
        chi[k] = if land_pre[r] { chi[r] } else { 0.0 } + (1.0 / a).sqrt() * dx_m;
    }
    let cs = sorted((0..n).filter(|&k| land_pre[k]).map(|k| chi[k]).collect::<Vec<f32>>());
    let c_uk = ORACLE_P50_M / pct(&cs, 0.50).max(1e-6);
    let z_chi: Vec<f32> = (0..n)
        .map(|k| if land_pre[k] { SEA + (c_uk * chi[k]) / n2m } else { br.data[k] })
        .collect();
    let b3 = build_field_with_floor(Knobs::passes(2), Some(std::sync::Arc::new(z_chi)));
    // Finding 97 B's two best candidates. k_s = 0.0451 is block B's declared calibration: the
    // median of S·sqrt(A) over the DELIVERED channel cells, i.e. the delivered field's own Flint
    // intercept, so the floor imports no relief target from outside.
    let b2 = build_field(Knobs { slope_floor_uk: Some(0.0451), ..Knobs::passes(2) });
    let ab = build_field(Knobs {
        slope_floor_uk: Some(0.0451),
        depression_floor: true,
        ..Knobs::passes(2)
    });
    eprintln!("   fields rebuilt in {:.1} s", t0.elapsed().as_secs_f64());

    // Finding 96's own global ramp, reproduced exactly
    let mut lm: Vec<f32> =
        (0..n).filter(|&k| delivered.data[k] > SEA).map(|k| delivered.data[k]).collect();
    lm = sorted(lm);
    let (glo, ghi) = (pct(&lm, 0.01), pct(&lm, 0.99));
    eprintln!(
        "\n   Finding 96's GLOBAL ramp: {:.0} m .. {:.0} m over 235 codes => **{:.1} m per grey \
         level**",
        c1_altitude_norm_to_metres(glo, &ss),
        c1_altitude_norm_to_metres(ghi, &ss),
        (ghi - glo) * n2m / 235.0
    );

    for (label, f) in [("delivered", &delivered), ("b3", &b3), ("b2", &b2), ("a1_b2", &ab)] {
        let t = crop(f);
        let mut tl: Vec<f32> = t.data.iter().copied().filter(|&v| v > SEA).collect();
        tl = sorted(tl);
        let (llo, lhi) = (pct(&tl, 0.01), pct(&tl, 0.99));
        let gr = ramp(&t, glo, ghi);
        let lr = ramp(&t, llo, lhi);
        let hs = hillshade(&t, n2m);
        eprintln!(
            "   {label:<10} tile p1..p99 = {:.0}..{:.0} m . GLOBAL ramp uses **{} codes** \
             ({:.1} m/code) . LOCAL ramp uses {} codes ({:.2} m/code) . hillshade {} codes",
            c1_altitude_norm_to_metres(llo, &ss),
            c1_altitude_norm_to_metres(lhi, &ss),
            codes(&gr),
            (ghi - glo) * n2m / 235.0,
            codes(&lr),
            (lhi - llo) * n2m / 235.0,
            codes(&hs)
        );
        save(&gr, &format!("{label}_u8_global.png"));
        save(&lr, &format!("{label}_u8_local.png"));
        save(&hs, &format!("{label}_hillshade.png"));
    }
    eprintln!("\n==========  end . total {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
