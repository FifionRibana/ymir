//! ADR 0001 Finding 112 — **is the comb in the FLOW FIELD or in the TRACER?**
//!
//! Finding 111 left the drawn network carrying a 4-axis signature of **0.4661 at 1.5 km** on the
//! comb tile, against Finding 97's isotropic calibration of **0.0402** — twelve times the floor —
//! over a terrain with no axis of its own (R2 0.0291, H(θ) 0.9835). The candidate remedy is to
//! stop using D8 for the drawn network. This bench measures whether that *can* work. **It wires
//! nothing.**
//!
//! ## What the grep found, and it answers the round's stop-question with a "yes, but"
//!
//! ⛔ **The remedy was evaluated at the production grid and judged unnecessary** — Finding 11,
//! TASK 5 (ADR L623): *"at 8192² the D8 rivers sit in the MFD valleys (**median offset 6 m, only
//! 18 % above a carved hollow**) — **acceptable**. At 2048² they diverge (30 m, 61 %)… At the
//! 8192² production grid it is a minor issue."* Finding 12 (L668) restates it as *"(ii) snap D8
//! segments to the MFD thalweg (or route rivers on MFD) for the 39 % misalignment"* — the 39 %
//! being the **2048²** figure.
//!
//! ⚠️ **But what was measured there is WHERE the rivers are, not WHAT SHAPE they have.** A river
//! can sit in the right valley to 6 m and still be drawn as a staircase of 45° steps. Finding
//! 111's R8 is a statement about form; Finding 11's 6 m is a statement about position. **The
//! recorded dismissal does not cover this defect**, so the round proceeds — with the dismissal
//! read and quoted, which is what rule 11 asks.
//!
//! ## ⛔ P0 — there are TWO remedies under one name, and only one is reachable without new code
//!
//! `mfd_accumulation` (flow.rs:866) returns a **`GridF32` accumulation field**. There is **no MFD
//! path anywhere in the crate** — `grep "fn .*mfd"` over `ymir-core/src` returns exactly that one
//! function. MFD disperses flow to several downslope neighbours, so it has **no single successor
//! per cell**, and `extract_rivers` — which walks `flow_result.direction` — cannot consume it.
//!
//! 1. **Threshold on MFD, trace on D8** — no new code. Changes *which cells are channels*, not
//!    *what shape a channel has*. The geometry still comes from eight directions.
//! 2. **Threshold AND trace on a single-path field** (steepest descent on the MFD-smoothed
//!    surface, or a real D∞ tracer à la Tarboton) — the only form that can move the geometry, and
//!    it does not exist in this crate.
//!
//! ⇒ **This bench measures the FIELD, not a re-extraction**, because the field is what decides
//! whether remedy 2 is worth writing: `R8` of `log A` under D8 against under MFD, on the same
//! tile with the same instrument. **If the MFD accumulation field is as combed as the D8 one, no
//! tracer will help and the owner is upstream of both.**
//!
//! Run: cargo test -p ymir-core --release --test f112_mfd -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, aniso, build_field_seed};
use std::time::Instant;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::terrain::flow::{breach_monotone, mfd_accumulation};

const DOMAIN_KM: f32 = 400.0;
const S_EQ: f32 = 0.024;
const TILE: usize = 512;
/// Finding 110-B's comb tile.
const COMB: (usize, usize) = (5120, 4096);
/// Finding 11's own MFD exponent, the one that carves the shipped terrain.
const MFD_P: f32 = 2.0;

fn crop(g: &GridF32, x0: usize, y0: usize, w: usize) -> GridF32 {
    let h = g.height;
    let (cw, ch) = (TILE.min(g.width - x0), TILE.min(h - y0));
    let mut o = GridF32::new(cw, ch, 0.0);
    for y in 0..ch {
        for x in 0..cw {
            o.data[y * cw + x] = g.data[(y0 + y) * w + x0 + x];
        }
    }
    o
}

#[test]
#[ignore]
fn f112_mfd() {
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 112 . is the comb in the FIELD or the TRACER?  ==========");

    let f = build_field_seed(Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }, PSEED);
    let (w, h) = (f.width, f.height);
    let n = w * h;
    let dr = c1_drainage_windowed(&f, None, &dcfg, &ss, DOMAIN_KM);
    let bf = breach_monotone(&f, &dr.flow.filled, &dr.lake_map, SEA, w, h);
    // ⚠️ Declared: both accumulations are computed from the SAME `filled` and the SAME
    // `direction` that the drawn network is traced on, so nothing but the partition differs.
    let d8 = &dr.flow.accumulation;
    let mfd = mfd_accumulation(&dr.flow.filled, &dr.flow.direction, SEA, MFD_P, w, h);

    // ── the two channel masks (Finding 65's 1.60 % vs 8.40 %) ───────────────
    let a_c_cells = RELIEF_V1_A_C_KM2 / cell_km2;
    let land: Vec<bool> = (0..n).map(|k| bf.data[k] > SEA).collect();
    let nland = land.iter().filter(|&&b| b).count();
    let cnt = |g: &GridF32| (0..n).filter(|&k| land[k] && g.data[k] >= a_c_cells).count();
    eprintln!(
        "\n   ── the two channel definitions on THIS field (Finding 65 measured 1.60 % / 8.40 %) ──\
         \n      A_c = {RELIEF_V1_A_C_KM2} km² = {a_c_cells:.1} cells · land {nland}\
         \n      **D8  {:>9} cells = {:>5.2} % of land**\
         \n      **MFD {:>9} cells = {:>5.2} % of land**  ⇒ ratio **×{:.2}**",
        cnt(d8),
        100.0 * cnt(d8) as f32 / nland as f32,
        cnt(&mfd),
        100.0 * cnt(&mfd) as f32 / nland as f32,
        cnt(&mfd) as f32 / cnt(d8).max(1) as f32
    );

    // ── THE DISCRIMINANT: R8 of log-accumulation, D8 against MFD ────────────
    // The Finding 97 instrument, unchanged, applied to the accumulation field instead of the
    // height. Calibration carried forward: 0.0402 isotropic, 0.9999 striped.
    let logf = |g: &GridF32| {
        let mut o = GridF32::new(w, h, 0.0);
        for k in 0..n {
            o.data[k] = (1.0 + g.data[k].max(0.0)).ln();
        }
        o
    };
    let (ld8, lmfd) = (logf(d8), logf(&mfd));
    eprintln!("\n   ── the discriminant: R8 of log A, same instrument, same windows ──");
    for (name, x0, y0) in
        [("THE COMB TILE", COMB.0, COMB.1), ("whole continent", usize::MAX, usize::MAX)]
    {
        let (g8, gm, lnd) = if x0 == usize::MAX {
            (ld8.clone(), lmfd.clone(), land.clone())
        } else {
            let c8 = crop(&ld8, x0, y0, w);
            let cm = crop(&lmfd, x0, y0, w);
            let ch = crop(&bf, x0, y0, w);
            let l: Vec<bool> = ch.data.iter().map(|&v| v > SEA).collect();
            (c8, cm, l)
        };
        let a8 = aniso(&g8, &lnd, 16);
        let am = aniso(&gm, &lnd, 16);
        eprintln!(
            "      {name:<16} · **D8  R8 {:.4}** R2 {:.4} θ {:>5.1}° C {:.3} H {:.4}\n      \
             {:<16} · **MFD R8 {:.4}** R2 {:.4} θ {:>5.1}° C {:.3} H {:.4} ⇒ **{}**",
            a8.r8,
            a8.r2,
            a8.theta2_deg,
            a8.c_p50,
            a8.h_theta,
            "",
            am.r8,
            am.r2,
            am.theta2_deg,
            am.c_p50,
            am.h_theta,
            if am.r8 < 0.5 * a8.r8 {
                "MFD is MUCH less combed -- a tracer is worth writing"
            } else if am.r8 < a8.r8 {
                "MFD is somewhat less combed"
            } else {
                "**MFD is NOT less combed -- the owner is upstream of the tracer**"
            }
        );
    }
    // ── the calibration, re-run here so the numbers above are not on trust ──
    let mut iso = GridF32::new(256, 256, 0.0);
    let mut stripe = GridF32::new(256, 256, 0.0);
    for y in 0..256usize {
        for x in 0..256usize {
            let mut v = 0.0f32;
            for i in 0..64 {
                let th = i as f32 * 2.399_963; // golden angle: 64 plane waves, no preferred axis
                v += ((x as f32 * th.cos() + y as f32 * th.sin()) * 0.3).sin();
            }
            iso.data[y * 256 + x] = v;
            stripe.data[y * 256 + x] =
                (((x as f32 + y as f32) / 14.0) * std::f32::consts::TAU).sin();
        }
    }
    let all = vec![true; 256 * 256];
    eprintln!(
        "      CALIBRATION (re-run, not cited): isotropic **{:.4}** · striped 45° **{:.4}** \
         (Finding 97: 0.0402 / 0.9999)",
        aniso(&iso, &all, 16).r8,
        aniso(&stripe, &all, 16).r8
    );

    // ── the control: the biggest watercourse must still be a channel in MFD ─
    let big = (0..n)
        .filter(|&k| land[k])
        .max_by(|&a, &b| d8.data[a].total_cmp(&d8.data[b]))
        .expect("land");
    eprintln!(
        "\n   ── control · the cell of maximum D8 accumulation ──\n      ({}, {}) · D8 **{:.0} \
         km²** · MFD **{:.0} km²** ⇒ still a channel under MFD: **{}**",
        big % w,
        big / w,
        d8.data[big] * cell_km2,
        mfd.data[big] * cell_km2,
        if mfd.data[big] >= a_c_cells { "YES" } else { "**NO -- the remedy would delete it**" }
    );
    // ⚠️ One cell is not a control. The max-D8 cell reading a small MFD value could mean either
    // "MFD has no trunk" or "MFD's trunk is one cell away" — opposite conclusions for whether a
    // tracer is worth writing. Both maxima and the neighbourhood decide it.
    let bigm = (0..n)
        .filter(|&k| land[k])
        .max_by(|&a, &b| mfd.data[a].total_cmp(&mfd.data[b]))
        .expect("land");
    let nbr = |g: &GridF32, k: usize, r: i32| -> f32 {
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        let mut m = 0.0f32;
        for dy in -r..=r {
            for dx in -r..=r {
                let nx = (x + dx).rem_euclid(w as i32) as usize;
                let ny = (y + dy).rem_euclid(h as i32) as usize;
                m = m.max(g.data[ny * w + nx]);
            }
        }
        m
    };
    let far = ((bigm % w) as i64 - (big % w) as i64)
        .abs()
        .max(((bigm / w) as i64 - (big / w) as i64).abs());
    eprintln!(
        "      max MFD cell ({}, {}) · MFD **{:.0} km²** · D8 there {:.0} km² · Chebyshev distance \
         from the max-D8 cell **{far} cells**\n      around the max-D8 cell, max MFD within r = 1 \
         / 4 / 16: **{:.1} / {:.1} / {:.1} km²** ⇒ **{}**",
        bigm % w,
        bigm / w,
        mfd.data[bigm] * cell_km2,
        d8.data[bigm] * cell_km2,
        nbr(&mfd, big, 1) * cell_km2,
        nbr(&mfd, big, 4) * cell_km2,
        nbr(&mfd, big, 16) * cell_km2,
        if nbr(&mfd, big, 16) > 0.25 * d8.data[big] {
            "the trunk EXISTS in MFD, merely displaced -- a tracer is writable"
        } else {
            "**the trunk is DISSOLVED, not displaced -- a thresholded MFD network has no main stem**"
        }
    );
    eprintln!("\n==========  end Finding 112 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
