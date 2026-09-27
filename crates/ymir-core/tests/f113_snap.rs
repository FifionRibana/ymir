//! ADR 0001 Finding 113 — **the snap is refuted before it is written, and by the dossier.**
//!
//! Block B asked to price "D8 for the hierarchy + snap the path to the MFD thalweg" — the half of
//! Finding 12's remedy (ii) that Finding 112 left open. **Rule 11 closed it instead**, and the
//! earliest hit is decisive:
//!
//! ⛔ **Finding 15 (L727) — "The river 'offset' was a metric artifact; rivers already sit in the
//! thalweg":**
//!   1. *"**MFD dominant-flow receiver ≡ D8 steepest** (92.9 % of land cells; the rest are
//!      flat-cell tie-breaks). `argmaxⱼ slopeⱼᵖ = argmaxⱼ slopeⱼ` — the dominant MFD path IS the
//!      D8 line, so **extracting rivers from MFD returns the same polyline. It cannot change the
//!      offset.**"*
//!   2. *"The offset metric was wrong… The correct test is **TRANSVERSE**: is the river cell ≤ its
//!      two banks perpendicular to flow."* ⇒ breached: **94 % in-thalweg, transverse p50 −1.8 m,
//!      p90 0.0 m**.
//!
//! ⇒ **There is nothing to snap to**, and Finding 112's quotation of Finding 11's *"median offset
//! 6 m, acceptable"* cites the metric Finding 15 retired. Finding 15 also retro-explains Finding
//! 112's `×1.04` mask agreement: the two operators pick the same receiver 92.9 % of the time.
//!
//! ⚠️ **Both grep predictions were wrong** — the reviewer's and mine. And my P1 (*"the snap must
//! search PERPENDICULAR to the flow, a disc snap on a monotone field walks downhill"*) is a
//! reinvention of Finding 15's own correction, which I should have grepped before predicting.
//!
//! ## What this bench does instead
//!
//! 1. **Re-measures Finding 15's two numbers on THIS field** — it was 2048², this is 8192² with
//!    Finding 109 ON. A cited number that decides a round gets re-read (rule 12).
//! 2. **Prices the only surviving remedy.** If the line is in the right place and still reads
//!    R8 0.4661 at 1.5 km, the defect is not WHERE the polyline is but WHAT IT IS MADE OF: a
//!    sequence of **cell centres**, so every step is one of eight vectors. The only thing left is
//!    to stop drawing cell centres — resample the polyline in continuous coordinates. That moves
//!    no cell, changes no threshold, touches no concentration, and is a pure export/render change.
//!
//! Run: cargo test -p ymir-core --release --test f113_snap -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, build_field_seed, pct, sorted};
use std::time::Instant;
use ymir_core::climate::c1_climate_placed;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, DrainageClimate, c1_drainage_windowed};
use ymir_core::tectonics_c1::hd_assembly::assemble_hd_drainage;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE, breach_monotone};

const DOMAIN_KM: f32 = 400.0;
const GEO_RATIO: f32 = 7.5;
const S_EQ: f32 = 0.024;
const TILE: usize = 512;
const COMB: (usize, usize) = (5120, 4096);

/// `R2`, `R8` of a set of axial directions (mod π).
fn axial(t: &[f32]) -> (f32, f32) {
    if t.is_empty() {
        return (0.0, 0.0);
    }
    let (mut c2, mut s2, mut c8, mut s8) = (0f64, 0f64, 0f64, 0f64);
    for &a in t {
        let a = a as f64;
        c2 += (2.0 * a).cos();
        s2 += (2.0 * a).sin();
        c8 += (8.0 * a).cos();
        s8 += (8.0 * a).sin();
    }
    let n = t.len() as f64;
    (((c2 * c2 + s2 * s2).sqrt() / n) as f32, ((c8 * c8 + s8 * s8).sqrt() / n) as f32)
}

/// Chords of a polyline given in CONTINUOUS coordinates, at scale `l` samples.
fn chords_f(p: &[(f32, f32)], l: usize) -> Vec<f32> {
    let mut out = Vec::new();
    let mut i = 0usize;
    while i + l < p.len() {
        let (dx, dy) = (p[i + l].0 - p[i].0, p[i + l].1 - p[i].1);
        if dx != 0.0 || dy != 0.0 {
            let mut t = dy.atan2(dx);
            if t < 0.0 {
                t += std::f32::consts::PI;
            }
            out.push(t);
        }
        i += l.max(1);
    }
    out
}

/// A moving average over ±k samples with the endpoints pinned — the cheapest thing that stops a
/// polyline being a sequence of cell centres. It moves no cell and changes no threshold.
fn smooth(p: &[(f32, f32)], k: usize) -> Vec<(f32, f32)> {
    if p.len() < 2 * k + 3 {
        return p.to_vec();
    }
    let mut o = p.to_vec();
    for i in k..p.len() - k {
        let (mut sx, mut sy) = (0f32, 0f32);
        for j in i - k..=i + k {
            sx += p[j].0;
            sy += p[j].1;
        }
        o[i] = (sx / (2 * k + 1) as f32, sy / (2 * k + 1) as f32);
    }
    o
}

#[test]
#[ignore]
fn f113_snap() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let t0 = Instant::now();
    eprintln!("\n==========  Finding 113 . the snap, refuted by Finding 15  ==========");

    // ⚠️ METHOD RULE 18 applied to my own measurement. Control 2 first came back at 50.5 %
    // against Finding 15's 94 %, and TWO things differ from Finding 15 at once: the resolution
    // (8192² against 2048²) and the age closure. Attributing that gap to either without the
    // other's control is exactly what rule 18 forbids — and it matters, because if the closure
    // halves thalweg residence, that belongs on the closure's ledger before it is promoted.
    for (label, knobs) in [
        ("F109 OFF (delivered)", Knobs::passes(2)),
        ("F109 ON  (s_eq 0.024)", Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) }),
    ] {
        eprintln!("\n########  {label}  ########");
        let f = build_field_seed(knobs, PSEED);
        let (w, h) = (f.width, f.height);
        let n = w * h;
        let pre = c1_drainage_windowed(&f, None, &dcfg, &ss, DOMAIN_KM);
        let bf = breach_monotone(&f, &pre.flow.filled, &pre.lake_map, SEA, w, h);
        let climate = c1_climate_placed(&bf, &ss, 45.0, 40.0, &PrecipParams::default(), DOMAIN_KM);
        let dclim = DrainageClimate {
            precip_internal: &climate.precipitation,
            temperature: &climate.temperature,
        };
        let pre2 = c1_drainage_windowed(&f, None, &dcfg, &ss, DOMAIN_KM);
        let dr = assemble_hd_drainage(
            &bf,
            &dclim,
            Some(pre2),
            &dcfg,
            &ss,
            DOMAIN_KM,
            GEO_RATIO,
            None,
            false,
        )
        .drainage;

        // ── CONTROL 1 · Finding 15 claim 1, re-measured at 8192² ────────────────
        // The MFD dominant receiver is `argmax_j slope_j^p`, which for any p > 0 is `argmax_j slope_j`
        // — the steepest descent, i.e. what D8 already picks. Measured against `flow.direction`.
        let fill = &pre.flow.filled;
        let (mut same, mut tot) = (0usize, 0usize);
        for k in 0..n {
            if bf.data[k] <= SEA || pre.flow.direction[k] == DIR_NONE {
                continue;
            }
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            let (mut best, mut bd) = (-1f32, 255u8);
            for d in 0..8u8 {
                let nx = (x + D8_DX[d as usize]).rem_euclid(w as i32) as usize;
                let ny = (y + D8_DY[d as usize]).rem_euclid(h as i32) as usize;
                let diag = D8_DX[d as usize] != 0 && D8_DY[d as usize] != 0;
                let len = if diag { std::f32::consts::SQRT_2 } else { 1.0 };
                let s = (fill.data[k] - fill.data[ny * w + nx]) / len;
                if s > best {
                    best = s;
                    bd = d;
                }
            }
            tot += 1;
            if bd == pre.flow.direction[k] {
                same += 1;
            }
        }
        eprintln!(
            "\n   ── CONTROL 1 · Finding 15 claim 1, re-measured here (it was 2048²) ──\n      MFD \
         dominant receiver == D8 receiver on **{:.2} %** of {tot} routed land cells (Finding 15 at \
         2048²: **92.9 %**) ⇒ **{}**",
            100.0 * same as f32 / tot.max(1) as f32,
            if same as f32 / tot.max(1) as f32 > 0.85 {
                "CONFIRMED -- an MFD tracer returns the same polyline"
            } else {
                "**NOT confirmed at this resolution**"
            }
        );

        // ── CONTROL 2 · Finding 15 claim 2: TRANSVERSE thalweg residence ────────
        let mut off = Vec::new();
        let mut inth = 0usize;
        for s in &dr.rivers.segments {
            for i in 1..s.points.len().saturating_sub(1) {
                let (x, y) = (s.points[i].0 as i32, s.points[i].1 as i32);
                let (px, py) = (
                    s.points[i + 1].0 as i32 - s.points[i - 1].0 as i32,
                    s.points[i + 1].1 as i32 - s.points[i - 1].1 as i32,
                );
                if px == 0 && py == 0 {
                    continue;
                }
                // the perpendicular, quantised to the nearest lattice step
                let (qx, qy) = (-py.signum(), px.signum());
                let at = |dx: i32, dy: i32| -> f32 {
                    let nx = (x + dx).rem_euclid(w as i32) as usize;
                    let ny = (y + dy).rem_euclid(h as i32) as usize;
                    bf.data[ny * w + nx] * n2m
                };
                let (me, a, b) = (at(0, 0), at(qx, qy), at(-qx, -qy));
                off.push(me - a.min(b));
                if me <= a && me <= b {
                    inth += 1;
                }
            }
        }
        let so = sorted(off.clone());
        eprintln!(
            "   ── CONTROL 2 · Finding 15 claim 2, TRANSVERSE (the metric Finding 15 repaired) ──\n   \
         **{:.1} % in-thalweg** of {} river points · transverse offset p50 **{:+.2} m** · p90 \
         **{:+.2} m** (Finding 15, breached, 2048²: **94 % · −1.8 m · 0.0 m**)",
            100.0 * inth as f32 / off.len().max(1) as f32,
            off.len(),
            pct(&so, 0.50),
            pct(&so, 0.90)
        );

        // ── the only surviving remedy: stop drawing CELL CENTRES ────────────────
        eprintln!(
            "\n   ── the only remedy left: the polyline is a sequence of CELL CENTRES ──\n      moving \
         average, endpoints pinned. Moves no cell, changes no threshold, no concentration."
        );
        let in_comb = |p: &(u32, u32)| {
            (p.0 as usize) >= COMB.0
                && (p.0 as usize) < COMB.0 + TILE
                && (p.1 as usize) >= COMB.1
                && (p.1 as usize) < COMB.1 + TILE
        };
        for k in [0usize, 1, 2, 4, 8] {
            let (mut all, mut comb, mut dev) = (Vec::new(), Vec::new(), Vec::new());
            for s in &dr.rivers.segments {
                if s.points.len() < 4 {
                    continue;
                }
                let raw: Vec<(f32, f32)> =
                    s.points.iter().map(|&(x, y)| (x as f32, y as f32)).collect();
                let sm = if k == 0 { raw.clone() } else { smooth(&raw, k) };
                for (a, b) in raw.iter().zip(sm.iter()) {
                    dev.push(((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt());
                }
                all.extend(chords_f(&sm, 32));
                let inside: Vec<(f32, f32)> = s
                    .points
                    .iter()
                    .zip(sm.iter())
                    .filter(|(p, _)| in_comb(p))
                    .map(|(_, q)| *q)
                    .collect();
                comb.extend(chords_f(&inside, 32));
            }
            let (r2a, r8a) = axial(&all);
            let (r2c, r8c) = axial(&comb);
            let sd = sorted(dev);
            eprintln!(
                "      ±{k:>1} samples: continent **R8 {r8a:.4}** (R2 {r2a:.3}) · comb tile **R8 \
             {r8c:.4}** (R2 {r2c:.3}, n = {}) · deviation p50 **{:.2}** p90 **{:.2}** cells",
                comb.len(),
                pct(&sd, 0.50),
                pct(&sd, 0.90)
            );
        }
        eprintln!(
            "      (Finding 111 on the same tile, unsmoothed: **R8 0.4661** at chord 32; isotropic \
         floor **0.031–0.040**)"
        );
    }
    eprintln!("\n==========  end Finding 113 . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
