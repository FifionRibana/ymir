//! ADR 0001 Finding 117-B — **who creates the pits, stage by stage, in the PIPELINE's order.**
//!
//! Finding 116 killed the FBM-amplitude sweep (the knob is inert since C-1, ADR "The DEAD KNOB")
//! but left two anchors: **437** pits > 0.1 m with the FBM off and **6 895** with it on, before
//! incision; **10 867** on the delivered field with the age closure on. This bench walks the
//! production pipeline one module at a time and reports the Δ each one contributes.
//!
//! ## ⛔ Rule 12 — the order is the code's, not the round's
//!
//! `upscale_from_c1` (production_upscale.rs) runs: `upscale_with_fbm` (bicubic **and** FBM in one
//! call, :459) → **C-2 craters stamped** (:461, *"AFTER the FBM"*) → incision with the C-3/C-3b
//! erodibility field (:484–:506; lithology and fracture are **K-modulators inside the incision**,
//! not stages of their own) → **C-2 rim reconstruction** (:579, *"after ALL erosion"*) →
//! bathymetry (:597). The breach is not in `upscale_from_c1` at all; the benches and the HD
//! assembly apply it to the output.
//!
//! The round's list (… C-3 · C-3b · C-2 · incision · breach · bathymetry) differs in three
//! places: craters come BEFORE the incision (and again after it), litho/fracture act only
//! THROUGH the incision, and the breach comes AFTER bathymetry on this path. **The chain below is
//! the code's order.** Two things this chain cannot separate, declared: the incision from the
//! C-2 rim reconstruction (both switch with the stream-power flag), and bicubic overshoot from
//! the coarse field's own 16 pits (Finding 41's number, a different instrument).
//!
//! ## ⛔ The control that can fail
//!
//! The round's control — *"the sum of the Δ recovers the 10 867 to ±1 %"* — is a telescoping sum
//! and holds by construction. The control that MEANS something is that the last cumulative stage
//! is **byte-identical** to the delivered field built in one go; only then are the Δ attributed to
//! the pipeline rather than to a chain that merely ends near the same number.
//!
//! Run: cargo test -p ymir-core --release --test f117_stages -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, build_field_seed, pct, sorted};
use std::time::Instant;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;
use ymir_core::terrain::flow::breach_monotone;

const DOMAIN_KM: f32 = 400.0;
const S_EQ: f32 = 0.024;

fn hash(f: &GridF32) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for v in &f.data {
        for b in v.to_bits().to_le_bytes() {
            h ^= b as u64;
            h = h.wrapping_mul(0x1000_0000_01b3);
        }
    }
    h
}

/// Finding 41's pit definition: 8-connected components the priority flood raises by > `t`.
/// Returns (count, count ≤ 2 cells, depth p50, total area km²).
fn pits(
    fill: &GridF32,
    raw: &GridF32,
    n2m: f32,
    t: f32,
    cell_km2: f32,
) -> (usize, usize, f32, f32) {
    let (w, h) = (raw.width, raw.height);
    let n = w * h;
    let deep: Vec<bool> = (0..n).map(|k| (fill.data[k] - raw.data[k]) * n2m > t).collect();
    let mut seen = vec![false; n];
    let (mut count, mut tiny, mut area) = (0usize, 0usize, 0f32);
    let mut depths = Vec::new();
    for s in 0..n {
        if !deep[s] || seen[s] {
            continue;
        }
        let mut q = vec![s];
        seen[s] = true;
        let (mut cells, mut dmax) = (0usize, 0f32);
        while let Some(k) = q.pop() {
            cells += 1;
            dmax = dmax.max((fill.data[k] - raw.data[k]) * n2m);
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    let nk = (y + dy).rem_euclid(h as i32) as usize * w
                        + (x + dx).rem_euclid(w as i32) as usize;
                    if deep[nk] && !seen[nk] {
                        seen[nk] = true;
                        q.push(nk);
                    }
                }
            }
        }
        count += 1;
        if cells <= 2 {
            tiny += 1;
        }
        area += cells as f32 * cell_km2;
        depths.push(dmax);
    }
    let sd = sorted(depths);
    (count, tiny, if sd.is_empty() { 0.0 } else { pct(&sd, 0.50) }, area)
}

#[test]
#[ignore]
fn f117_stages() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    let cell_km2 = CELL_KM * CELL_KM;
    let mut dcfg = C1DrainageConfig::default();
    dcfg.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    dcfg.thresholds.full_tree = false;
    let t0 = Instant::now();
    eprintln!(
        "\n==========  Finding 117-B . the pits, module by module, in the code's order  ====="
    );

    // The pre-incision base every stage builds on: no stream power, closure irrelevant.
    let off = Knobs::no_incision();
    // The incision stages carry the age closure, because the 10 867 the round names is the
    // closure-ON delivered field of Finding 116.
    let on = Knobs { slope_floor_abs: Some(S_EQ), ..Knobs::passes(2) };
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
            "S4 + incision, uniform K (+ C-2 rims)",
            Knobs { lithology_off: true, fracture_off: true, bathymetry_off: true, ..on },
        ),
        ("S5 + C-3 lithology in K", Knobs { fracture_off: true, bathymetry_off: true, ..on }),
        ("S6 + C-3b fracture in K", Knobs { bathymetry_off: true, ..on }),
        ("S7 + bathymetry  (= delivered?)", on),
    ];

    eprintln!(
        "\n   {:<40} {:>8} {:>7} {:>8} {:>9} {:>7} {:>7} {:>9}",
        "stage", ">0.1 m", "Δ", "≤2cell", "p50 m", ">20 m", "Δ", "area km²"
    );
    let (mut prev1, mut prev20) = (0i64, 0i64);
    let mut last: Option<GridF32> = None;
    for (label, knobs) in chain {
        let f = build_field_seed(knobs, PSEED);
        let dr = c1_drainage_windowed(&f, None, &dcfg, &ss, DOMAIN_KM);
        let (c1, tiny, p50, area) = pits(&dr.flow.filled, &f, n2m, 0.1, cell_km2);
        let (c20, _, _, _) = pits(&dr.flow.filled, &f, n2m, 20.0, cell_km2);
        eprintln!(
            "   {label:<40} {c1:>8} {:>+7} {tiny:>8} {p50:>9.2} {c20:>7} {:>+7} {area:>9.1}",
            c1 as i64 - prev1,
            c20 as i64 - prev20
        );
        prev1 = c1 as i64;
        prev20 = c20 as i64;
        last = Some(f);
    }

    // ── the control that can fail: S7 must BE the delivered field ───────────
    let s7 = last.expect("chain");
    let delivered = build_field_seed(on, PSEED);
    let same = hash(&s7) == hash(&delivered);
    let diff = (0..s7.data.len()).filter(|&k| s7.data[k] != delivered.data[k]).count();
    eprintln!(
        "\n   ⛔ CONTROL · S7 vs the delivered field built in one go: **{}** ({diff} cells differ) \
         ⇒ {}",
        if same { "BIT-IDENTICAL" } else { "**DIFFERENT**" },
        if same {
            "the Δ column is the pipeline's, not the chain's"
        } else {
            "**the toggled order is NOT the pipeline's order; the Δ column is mis-attributed**"
        }
    );

    // ── S8 · the breach, applied to the delivered output as the benches and HD assembly do ──
    let (w, h) = (delivered.width, delivered.height);
    let pre = c1_drainage_windowed(&delivered, None, &dcfg, &ss, DOMAIN_KM);
    let bf = breach_monotone(&delivered, &pre.flow.filled, &pre.lake_map, SEA, w, h);
    let drb = c1_drainage_windowed(&bf, None, &dcfg, &ss, DOMAIN_KM);
    let (c1, tiny, p50, area) = pits(&drb.flow.filled, &bf, n2m, 0.1, cell_km2);
    let (c20, _, _, _) = pits(&drb.flow.filled, &bf, n2m, 20.0, cell_km2);
    eprintln!(
        "   {:<40} {c1:>8} {:>+7} {tiny:>8} {p50:>9.2} {c20:>7} {:>+7} {area:>9.1}",
        "S8 + breach (held lakes stay)",
        c1 as i64 - prev1,
        c20 as i64 - prev20
    );
    eprintln!(
        "\n   anchors from Finding 116: no-FBM pre-incision **437** · FBM pre-incision **6 895** · \
         delivered (closure ON) **10 867** / **921** > 20 m · Finding 41 coarse **16** (other \
         instrument)"
    );
    eprintln!(
        "\n==========  end Finding 117-B . {:.1} s  ==========\n",
        t0.elapsed().as_secs_f64()
    );
}
