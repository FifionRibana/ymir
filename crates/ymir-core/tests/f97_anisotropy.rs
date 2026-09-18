//! ADR 0001 Finding 97 block A — the anisotropy instrument, calibrated before any candidate.
//! **No production change.**
//!
//! Finding 96 was refuted by an image no statistic in its table could see: the χ floor prints the
//! D8 axes into the terrain. This bench builds the statistic that should have seen it, and
//! CALIBRATES it against known answers before it is allowed to judge anything.
//!
//! ⚠️ Two antecedents from the rule-11 grep change the design, and both are warnings:
//!
//!  * **L292** — the dossier has already tried an anisotropy metric on striations and it FAILED:
//!    *"the anisotropy knobs do NOT move the striation metric … either the visual striations are
//!    not controlled by these knobs, or the ±8-cell asymmetry metric is too coarse to isolate
//!    them — an honest gap"*. A new instrument has to beat that, so it gets negative controls at
//!    BOTH ends of its scale, not just the striped one the round asked for.
//!  * **Finding 84 / L5015** — the dossier already owns an isotropy statistic, `axis R`, with a
//!    calibrated reference (**0** = isotropic, target **< 0.1**), and Finding 76 measured its
//!    coast version at **0.000** un-incised against **0.879** shipped. But `axis R` is the SECOND
//!    circular harmonic, and the second harmonic is **blind by symmetry** to a 4-fold D8 pattern:
//!    orientations on {0°, 45°, 90°, 135°} give ⟨e^{2iθ}⟩ = (1 + i − 1 − i)/4 = 0. The harmonic
//!    that resonates with a 45° period is the **EIGHTH**. Both are reported, and the gap between
//!    them is itself a result about the dossier's existing instrument.
//!
//! Run: cargo test -p ymir-core --release --test f97_anisotropy -- --ignored --nocapture

mod common;

use common::{
    Knobs, SEA, aniso, build_field, build_field_with_floor, pct, sorted, tile, variogram,
    worst_tile,
};
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

fn row(label: &str, f: &GridF32, land: &[bool]) {
    for win in [8usize, 16, 32] {
        let a = aniso(f, land, win);
        eprintln!(
            "   {label:<22} w={win:<3} windows {:<8} C p50 **{:.3}** p90 {:.3} . C>0.7 \
             **{:.2} %** . H(theta) **{:.4}** . R2 {:.4} . **R8 {:.4}** (C-weighted {:.4})",
            a.windows, a.c_p50, a.c_p90, a.share_hi, a.h_theta, a.r2, a.r8, a.r8_w
        );
    }
}

/// A smooth ramp plus a sinusoid along the 45° D8 diagonal, at Finding 76's own measured fur
/// wavelength of **14 cells**. The instrument must read this near the top of its scale.
fn synth_striped(n: usize, amp: f32) -> GridF32 {
    let mut g = GridF32::new(n, n, 0.0);
    for y in 0..n {
        for x in 0..n {
            let ramp = 0.30 - 0.10 * (y as f32 / n as f32);
            let s = ((x + y) as f32 * std::f32::consts::TAU / 14.0).sin();
            g.data[y * n + x] = ramp + amp * s;
        }
    }
    g
}

/// A provably ISOTROPIC random field: a sum of plane waves with uniformly distributed directions.
/// This fixes the BOTTOM of the scale -- a finite-difference gradient on a square grid is already
/// slightly axis-biased, so "isotropic" is not R8 = 0 exactly, and the instrument must be told
/// what its own floor is before it calls anything anisotropic.
fn synth_isotropic(n: usize, amp: f32) -> GridF32 {
    let mut g = GridF32::new(n, n, 0.0);
    let nw = 64usize;
    let mut st = 0x2545_F491_4F6C_DD1Du64;
    let mut rnd = || {
        st ^= st << 13;
        st ^= st >> 7;
        st ^= st << 17;
        (st >> 11) as f64 / (1u64 << 53) as f64
    };
    let waves: Vec<(f32, f32, f32, f32)> = (0..nw)
        .map(|_| {
            let ang = (rnd() * std::f64::consts::TAU) as f32;
            let k = (2.0 + 26.0 * rnd() as f32) * std::f32::consts::TAU / n as f32;
            (k * ang.cos(), k * ang.sin(), (rnd() * std::f64::consts::TAU) as f32, 1.0)
        })
        .collect();
    for y in 0..n {
        for x in 0..n {
            let mut v = 0.0f32;
            for &(kx, ky, ph, a) in &waves {
                v += a * (kx * x as f32 + ky * y as f32 + ph).sin();
            }
            g.data[y * n + x] = 0.30 + amp * v / (nw as f32).sqrt();
        }
    }
    g
}

#[test]
#[ignore]
fn f97_anisotropy() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = (DOMAIN_KM / 8192.0) * (DOMAIN_KM / 8192.0);
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    eprintln!("\n==========  Finding 97 A . the anisotropy instrument, calibrated  ==========");
    eprintln!("   1 norm unit = {n2m:.1} m . cell {CELL_M:.2} m");

    // ── A2c: the two ends of the scale, on synthetics, BEFORE any real field ──
    let t0 = Instant::now();
    let iso = synth_isotropic(1024, 0.004);
    let strp = synth_striped(1024, 0.002);
    let all = vec![true; 1024 * 1024];
    eprintln!("\n-- A2c . the instrument's own scale, on synthetics (1024^2, all-land) --");
    row("SYNTH isotropic", &iso, &all);
    row("SYNTH striped 45deg/14c", &strp, &all);
    eprintln!("   (synthetics built in {:.1} s)", t0.elapsed().as_secs_f64());

    // ── the three real fields ────────────────────────────────────────────────
    let t1 = Instant::now();
    let pre = build_field(Knobs::no_incision());
    let (w, h) = (pre.width, pre.height);
    let n = w * h;
    let delivered = build_field(Knobs::passes(2));
    eprintln!("\n   builds: pre-incision + delivered in {:.1} s", t1.elapsed().as_secs_f64());

    // the χ floor of Finding 96, rebuilt, so B3 is available as the HIGH anchor
    let t2 = Instant::now();
    let d_pre = c1_drainage_windowed(&pre, None, &on, &ss, DOMAIN_KM);
    let br = breach_monotone(&pre, &d_pre.flow.filled, &d_pre.lake_map, SEA, w, h);
    let d_br = c1_drainage_windowed(&br, None, &on, &ss, DOMAIN_KM);
    let land_pre: Vec<bool> = (0..n).map(|k| br.data[k] > SEA).collect();
    let (order, unreached) = propagation_order(&d_br.flow.direction, |k| land_pre[k], w, h);
    assert_eq!(unreached, 0, "the chi traversal must reach every land cell");
    let a_c = RELIEF_V1_A_C_KM2;
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
        let a_km2 = (d_br.flow.accumulation.data[k] * cell_km2).max(a_c);
        chi[k] = if land_pre[r] { chi[r] } else { 0.0 } + (1.0 / a_km2).sqrt() * dx_m;
    }
    let mut cs: Vec<f32> = (0..n).filter(|&k| land_pre[k]).map(|k| chi[k]).collect();
    cs = sorted(cs);
    let c_uk = ORACLE_P50_M / pct(&cs, 0.50).max(1e-6);
    let z_chi: Vec<f32> = (0..n)
        .map(|k| if land_pre[k] { SEA + (c_uk * chi[k]) / n2m } else { br.data[k] })
        .collect();
    let b3 = build_field_with_floor(Knobs::passes(2), Some(std::sync::Arc::new(z_chi)));
    eprintln!(
        "   the chi floor rebuilt ((U/K)^(1/n) = {c_uk:.4}) and B3 built in {:.1} s",
        t2.elapsed().as_secs_f64()
    );

    // ── A2: the calibration table ────────────────────────────────────────────
    //
    // ⚠️ The round's table asks for the ORACLE as the low anchor. It does not exist as a raster:
    // Finding 95 recorded NUMBERS at 300 and 1 000 passes, and rebuilding the field costs 1 h 43
    // which the author refused. The low anchor is therefore the PRE-INCISION field, which is the
    // right one anyway -- Finding 76 measured its coast parallelism at R = 0.000 against the
    // shipped field's 0.879, so it is the dossier's own isotropic reference.
    eprintln!("\n-- A2 . calibration on the real fields (8192^2, fully-land windows only) --");
    let t3 = Instant::now();
    let l_pre: Vec<bool> = (0..n).map(|k| pre.data[k] > SEA).collect();
    let l_del: Vec<bool> = (0..n).map(|k| delivered.data[k] > SEA).collect();
    let l_b3: Vec<bool> = (0..n).map(|k| b3.data[k] > SEA).collect();
    row("PRE-INCISION (low)", &pre, &l_pre);
    row("DELIVERED", &delivered, &l_del);
    row("B3 chi floor (high)", &b3, &l_b3);
    eprintln!("   (three fields, three windows, in {:.1} s)", t3.elapsed().as_secs_f64());

    // ── the separation, which is the stop rule ───────────────────────────────
    //
    // The spread of R8 across the 8/16/32 windows of the PRE-INCISION field is the instrument's
    // own noise on a field that carries no D8 signal. A candidate is "separated" when it moves by
    // more than 3 of those.
    let pre8 = aniso(&pre, &l_pre, 16);
    let del8 = aniso(&delivered, &l_del, 16);
    let b38 = aniso(&b3, &l_b3, 16);
    eprintln!(
        "\n-- the separation at w=16 --\n   R8: pre {:.4} . delivered {:.4} . B3 {:.4} => B3 \
         moves **{:+.4}** from the delivered field, i.e. **x{:.2}**",
        pre8.r8,
        del8.r8,
        b38.r8,
        b38.r8 - del8.r8,
        b38.r8 / del8.r8.max(1e-6)
    );
    eprintln!(
        "   C p50: pre {:.3} . delivered {:.3} . B3 {:.3} . share C>0.7: {:.2} / {:.2} / {:.2} % \
         . H(theta): {:.4} / {:.4} / {:.4}",
        pre8.c_p50,
        del8.c_p50,
        b38.c_p50,
        pre8.share_hi,
        del8.share_hi,
        b38.share_hi,
        pre8.h_theta,
        del8.h_theta,
        b38.h_theta
    );

    // ── the worst tile per field: Finding 76's rule-9 crop trap, honoured ────
    eprintln!("\n-- the WORST 1024^2 tile per field (chosen per field, not fixed) --");
    for (label, f, l) in [
        ("delivered", &delivered, &l_del),
        ("B3 chi floor", &b3, &l_b3),
        ("pre-incision", &pre, &l_pre),
    ] {
        let (bx, by, a) = worst_tile(f, l);
        eprintln!(
            "   {label:<14} worst tile ({bx}, {by}) . windows {} . C p50 {:.3} . C>0.7 {:.2} % . \
             H(theta) {:.4} . R2 {:.4} . **R8 {:.4}**",
            a.windows, a.c_p50, a.share_hi, a.h_theta, a.r2, a.r8
        );
    }
    // ══ A3 — the three checks the stop rule demands BEFORE it is honoured ═══
    //
    // A2 does not separate B3 from the delivered field. Before concluding "the instrument is
    // false", three things have to be excluded, because each would make A2 a reading of the
    // wrong population rather than a reading of a blind instrument.
    eprintln!("\n-- A3a . the tile the EYE was shown (5120, 3072), the Finding 96 crop --");
    for (label, f, l) in [
        ("delivered", &delivered, &l_del),
        ("B3 chi floor", &b3, &l_b3),
        ("pre-incision", &pre, &l_pre),
    ] {
        let (g, gl) = tile(f, l, 5120, 3072, 1024);
        for win in [4usize, 8, 16] {
            let a = aniso(&g, &gl, win);
            eprintln!(
                "   {label:<14} w={win:<3} windows {:<6} C p50 {:.3} . C>0.7 {:.2} % . H(theta) \
                 {:.4} . R2 {:.4} . **R8 {:.4}**",
                a.windows, a.c_p50, a.share_hi, a.h_theta, a.r2, a.r8
            );
        }
    }

    // A3b -- the IMPRINT itself. The floor did not replace the terrain, it bounded it, so the
    // striping lives in `B3 - delivered` and the terrain's own structure is noise on top of it.
    // Measure the imprint ON the imprint, not on the sum.
    let dif = GridF32 {
        width: w,
        height: h,
        data: (0..n).map(|k| b3.data[k] - delivered.data[k]).collect(),
    };
    let l_dif: Vec<bool> = (0..n).map(|k| l_del[k] && l_b3[k]).collect();
    eprintln!("\n-- A3b . the IMPRINT, B3 minus delivered (the floor's own field) --");
    row("IMPRINT whole grid", &dif, &l_dif);
    let (dg, dgl) = tile(&dif, &l_dif, 5120, 3072, 1024);
    for win in [4usize, 8, 16] {
        let a = aniso(&dg, &dgl, win);
        eprintln!(
            "   IMPRINT tile(5120,3072) w={win:<3} windows {:<6} C p50 {:.3} . C>0.7 {:.2} % . \
             H(theta) {:.4} . R2 {:.4} . **R8 {:.4}**",
            a.windows, a.c_p50, a.share_hi, a.h_theta, a.r2, a.r8
        );
    }

    // A3c -- the directional semivariogram asks a different question: not "which way do the
    // gradients point" but "is the field rougher ACROSS one axis than another, and at what lag".
    // A stripe of wavelength L is a hole effect at L plus an anisotropy ratio above 1.
    eprintln!(
        "\n-- A3c . directional semivariogram on tile (5120,3072), gamma in m^2, dirs \
         [E, NE, N, NW] --"
    );
    for (label, f, l) in
        [("delivered", &delivered, &l_del), ("B3 chi floor", &b3, &l_b3), ("IMPRINT", &dif, &l_dif)]
    {
        let (g, gl) = tile(f, l, 5120, 3072, 1024);
        let v = variogram(&g, &gl, n2m, &[1, 2, 3, 4, 6, 8, 12, 16, 24]);
        let s: Vec<String> = v
            .iter()
            .map(|(lag, g4, r)| {
                format!("L{lag}: {:.0}/{:.0}/{:.0}/{:.0} r={r:.3}", g4[0], g4[1], g4[2], g4[3])
            })
            .collect();
        eprintln!("   {label:<14} {}", s.join(" . "));
    }

    eprintln!("\n==========  end Finding 97 A  ==========\n");
}
