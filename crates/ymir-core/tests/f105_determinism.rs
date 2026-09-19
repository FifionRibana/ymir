//! ADR 0001 Finding 105 — is the pipeline bit-reproducible run to run? **The control the
//! promotion gate needed and did not have.**
//!
//! The gate compared production against the bench and found 5.9 M cells differing by a **median
//! of 2 mm**. Before attributing that to the two-pass, one thing has to be excluded: whether
//! building the SAME configuration twice gives the same field at all.
//!
//! ⛔ The project's own design invariant says it must: *"Deterministic: same seed + config =
//! identical output. Rayon batches are sorted before processing."* (CLAUDE.md). If that holds,
//! the residual is the two-pass and is attributable. **If it does not, the Finding 83
//! bit-identity standard is unattainable on this path and every promotion gate written against it
//! has to change.**
//!
//! Run: cargo test -p ymir-core --release --test f105_determinism -- --ignored --nocapture

mod common;

use common::{Knobs, PSEED, SEA, build_field_seed, pct, sorted};
use std::time::Instant;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::production_upscale::c1_altitude_norm_to_metres;

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

#[test]
#[ignore]
fn f105_determinism() {
    let ss = SteinSteinParams::default();
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    eprintln!("\n==========  Finding 105 . is the pipeline bit-reproducible?  ==========");
    let t0 = Instant::now();

    let report = |label: &str, a: &GridF32, b: &GridF32| {
        let n = a.data.len();
        let diff = (0..n).filter(|&k| a.data[k] != b.data[k]).count();
        let d: Vec<f32> = (0..n)
            .filter(|&k| a.data[k] != b.data[k])
            .map(|k| ((a.data[k] - b.data[k]) * n2m).abs())
            .collect();
        let s = sorted(d);
        eprintln!(
            "   {label:<34} hashes **{}** . cells differing **{diff}** ({:.3} %) . |delta| p50 \
             {:.4} m p99 {:.3} m max {:.2} m",
            if hash(a) == hash(b) { "MATCH" } else { "**DIFFER**" },
            100.0 * diff as f64 / n as f64,
            if s.is_empty() { 0.0 } else { pct(&s, 0.50) },
            if s.is_empty() { 0.0 } else { pct(&s, 0.99) },
            s.last().copied().unwrap_or(0.0)
        );
    };

    // 1. the delivered field, twice — the plainest possible control
    let a = build_field_seed(Knobs::passes(2), PSEED);
    let b = build_field_seed(Knobs::passes(2), PSEED);
    report("delivered, built twice", &a, &b);

    // 2. the two-pass production path, twice
    let c = build_field_seed(Knobs { slope_floor_factor: Some(0.6), ..Knobs::passes(2) }, PSEED);
    let d = build_field_seed(Knobs { slope_floor_factor: Some(0.6), ..Knobs::passes(2) }, PSEED);
    report("production two-pass, built twice", &c, &d);

    // 3. and the explicit bench form, twice
    let e = build_field_seed(
        Knobs { slope_floor_uk: Some(0.045095 * 0.6), depression_floor: true, ..Knobs::passes(2) },
        PSEED,
    );
    let f = build_field_seed(
        Knobs { slope_floor_uk: Some(0.045095 * 0.6), depression_floor: true, ..Knobs::passes(2) },
        PSEED,
    );
    report("bench explicit form, built twice", &e, &f);

    eprintln!(
        "\n   ⚠️ a MATCH on all three means the pipeline is deterministic and the promotion \
         gate's residual is the two-pass. A DIFFER means the Finding 83 bit-identity standard \
         cannot be applied on this path, and every gate written against it has to change."
    );
    let land = (0..a.data.len()).filter(|&k| a.data[k] > SEA).count();
    eprintln!("   (land cells: {land})");
    eprintln!("\n==========  end . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
