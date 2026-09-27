//! ADR 0001 Finding 124 — **where does the field stop being a translation of itself?**
//!
//! Finding 123 measured that the viz's framing (torus offset x = 0.09375) and the benches' (0.0)
//! give two different worlds: 22.3 % of cells bit-identical after a roll of 768 cells, max |Δ|
//! 383 m. The resolution invariance of Finding 57 has no framing counterpart. This bench walks
//! Finding 117's stage chain (the code's order) under both framings and compares each stage with
//! the other framing's stage ROLLED by 768 cells. The first stage that differs is the one to name.
//!
//! Run: cargo test -p ymir-core --release --test f124_frame -- --ignored --nocapture

mod common;

use common::{Knobs, PSEED, build_field_seed};
use std::time::Instant;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;

const VIZ: [f64; 2] = [0.09375, 0.578_125];
const BENCH: [f64; 2] = [0.0, 0.578_125];

/// (bit-identical cells, cells, max |Δ| m, columns with any difference)
fn compare(
    v: &GridF32,
    b: &GridF32,
    shift: usize,
    ss: &SteinSteinParams,
) -> (usize, usize, f32, usize) {
    let (w, h) = (v.width, v.height);
    let (mut same, mut maxd) = (0usize, 0f32);
    let mut col = vec![false; w];
    for y in 0..h {
        for x in 0..w {
            let a = v.data[y * w + x];
            let c = b.data[y * w + (x + shift) % w];
            if a.to_bits() == c.to_bits() {
                same += 1;
            } else {
                col[x] = true;
                maxd = maxd.max(
                    (c1_altitude_norm_to_metres(a, ss) - c1_altitude_norm_to_metres(c, ss)).abs(),
                );
            }
        }
    }
    (same, w * h, maxd, col.iter().filter(|&&c| c).count())
}

#[test]
#[ignore]
fn f124_frame() {
    let ss = SteinSteinParams::default();
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 124 . the first non-periodic stage  ==========");
    let off = Knobs::no_incision();
    let on = Knobs::passes(2);
    let chain: [(&str, Knobs); 7] = [
        (
            "S1 bicubic only",
            Knobs {
                fbm_amp: Some(0.0),
                volcanism_off: true,
                lithology_off: true,
                fracture_off: true,
                bathymetry_off: true,
                ..off
            },
        ),
        (
            "S2 + FBM",
            Knobs {
                volcanism_off: true,
                lithology_off: true,
                fracture_off: true,
                bathymetry_off: true,
                ..off
            },
        ),
        (
            "S3 + C-2 craters stamped",
            Knobs { lithology_off: true, fracture_off: true, bathymetry_off: true, ..off },
        ),
        (
            "S4 + incision, uniform K",
            Knobs { lithology_off: true, fracture_off: true, bathymetry_off: true, ..on },
        ),
        ("S5 + C-3 lithology in K", Knobs { fracture_off: true, bathymetry_off: true, ..on }),
        ("S6 + C-3b fracture in K", Knobs { bathymetry_off: true, ..on }),
        ("S7 + bathymetry (= delivered)", on),
    ];
    eprintln!(
        "   {:<32} {:>14} {:>9} {:>12} {:>16}",
        "stage", "bit-identical", "share", "max |Δ| m", "columns ≠ / 8192"
    );
    for (label, k) in chain {
        let v = build_field_seed(Knobs { origin: Some(VIZ), ..k }, PSEED);
        let b = build_field_seed(Knobs { origin: Some(BENCH), ..k }, PSEED);
        let shift = (VIZ[0] * v.width as f64).round() as usize;
        let (same, n, maxd, cols) = compare(&v, &b, shift, &ss);
        eprintln!(
            "   {label:<32} {same:>14} {:>8.3}% {maxd:>12.3} {cols:>16}",
            100.0 * same as f64 / n as f64
        );
    }
    eprintln!(
        "\n==========  end Finding 124 frame . {:.1} s  ==========\n",
        t0.elapsed().as_secs_f64()
    );
}
