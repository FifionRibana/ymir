//! ADR 0001 Finding 105 — does PRODUCTION reproduce the bench, bit for bit, and what does it
//! cost? **This is the promotion gate.**
//!
//! Finding 83's pattern: if the production path reproduces the bench field **bit for bit**, every
//! gate measured on the bench field holds by construction. **The hash IS the guard table.**
//!
//! ⛔ And the structural point this round had to raise: the closure cannot ship as a constant,
//! because `k_s` is measured on the DELIVERED field — the output of the incision it modifies. The
//! shipped form is therefore a **two-pass bootstrap** (`FbmUpscaleConfig::slope_floor_factor`),
//! and its price is one extra incision plus one drainage, measured here against the author's
//! "≤ delivered + 30 s".
//!
//! Run: cargo test -p ymir-core --release --test f105_promote -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, build_field_seed, pct, sorted};
use std::time::Instant;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, c1_drainage_windowed};
use ymir_core::tectonics_c1::production_upscale::{c1_altitude_norm_to_metres, flint_intercept};
use ymir_core::terrain::flow::{D8_DX, D8_DY, DIR_NONE};

const DOMAIN_KM: f32 = 400.0;
const CELL_M: f32 = 400_000.0 / 8192.0;
const FACTOR: f32 = 0.6;

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// FNV-1a over the raw f32 bits — the hash Finding 88 rule 12 uses for flow fields, here on the
/// heightmap itself. Two fields with the same hash are the same terrain.
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

/// The bench's own `k_s`, kept here verbatim so the production helper can be checked against the
/// definition it was derived from rather than against itself.
fn bench_k_s(f: &GridF32, on: &C1DrainageConfig, ss: &SteinSteinParams, cell_km2: f32) -> f32 {
    let (w, h) = (f.width, f.height);
    let d = c1_drainage_windowed(f, None, on, ss, DOMAIN_KM);
    let n2m = c1_altitude_norm_to_metres(1.0, ss) - c1_altitude_norm_to_metres(0.0, ss);
    let mut v = Vec::new();
    for k in (0..w * h).step_by(7) {
        if f.data[k] <= SEA {
            continue;
        }
        let a = d.flow.accumulation.data[k] * cell_km2;
        if a < RELIEF_V1_A_C_KM2 {
            continue;
        }
        let dir = d.flow.direction[k];
        if dir == DIR_NONE {
            continue;
        }
        let nx = ((k % w) as i32 + D8_DX[dir as usize]).rem_euclid(w as i32) as usize;
        let ny = ((k / w) as i32 + D8_DY[dir as usize]).rem_euclid(h as i32) as usize;
        let diag = D8_DX[dir as usize] != 0 && D8_DY[dir as usize] != 0;
        let dx_m = if diag { CELL_M * std::f32::consts::SQRT_2 } else { CELL_M };
        let s = (f.data[k] - f.data[ny * w + nx]).max(0.0) * n2m / dx_m;
        if s > 0.0 {
            v.push(s * a.sqrt());
        }
    }
    pct(&sorted(v), 0.50)
}

#[test]
#[ignore]
fn f105_promote() {
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    let mut on = C1DrainageConfig::default();
    on.thresholds.head_km2 = RELIEF_V1_A_C_KM2;
    on.thresholds.full_tree = false;
    let s2 = splitmix64(PSEED);
    let s3 = splitmix64(s2);
    eprintln!("\n==========  Finding 105 . production against the bench, bit for bit  ==========");

    let mut all_ok = true;
    for (si, seed) in [PSEED, s2, s3].into_iter().enumerate() {
        eprintln!("\n========  SEED {} ({seed})  ========", si + 1);

        // ── the DELIVERED reference, and its cost ───────────────────────────
        let t = Instant::now();
        let delivered = build_field_seed(Knobs::passes(2), seed);
        let t_del = t.elapsed().as_secs_f64();

        // ── the BENCH candidate: k_s measured, then the field built from it ─
        let t = Instant::now();
        let ks_bench = bench_k_s(&delivered, &on, &ss, cell_km2);
        let bench = build_field_seed(
            Knobs {
                slope_floor_uk: Some(ks_bench * FACTOR),
                depression_floor: true,
                ..Knobs::passes(2)
            },
            seed,
        );
        let t_bench = t.elapsed().as_secs_f64();

        // ── the production helper, on the same field ────────────────────────
        let ks_prod = flint_intercept(&delivered, &ss, SEA, cell_km2, DOMAIN_KM);
        eprintln!(
            "   k_s: bench **{ks_bench:.6}** . production helper **{ks_prod:.6}** ⇒ **{}**",
            if ks_bench.to_bits() == ks_prod.to_bits() {
                "BIT-IDENTICAL"
            } else {
                "**DIFFERENT -- the specification does not reproduce the measurement**"
            }
        );

        // ── the PRODUCTION path: one call, the factor in the config ─────────
        let t = Instant::now();
        let prod =
            build_field_seed(Knobs { slope_floor_factor: Some(FACTOR), ..Knobs::passes(2) }, seed);
        let t_prod = t.elapsed().as_secs_f64();

        let (hb, hp) = (hash(&bench), hash(&prod));
        let ok = hb == hp;
        all_ok &= ok;
        let diff = (0..bench.data.len()).filter(|&k| bench.data[k] != prod.data[k]).count();
        eprintln!(
            "   hash bench **{hb:016x}** . production **{hp:016x}** ⇒ **{}** ({diff} cells differ)",
            if ok { "BIT-IDENTICAL" } else { "**DIFFERENT -- DO NOT PROMOTE**" }
        );
        eprintln!(
            "   ⏱ delivered **{t_del:.1} s** . bench (k_s + build) {t_bench:.1} s . **production \
             two-pass {t_prod:.1} s** ⇒ **{:+.1} s** against the author's +30 s budget ⇒ **{}**",
            t_prod - t_del,
            if t_prod - t_del <= 30.0 { "WITHIN" } else { "**OVER**" }
        );
    }
    eprintln!(
        "\n==========  end . bit-identity on all three seeds: **{}**  ==========\n",
        if all_ok { "YES" } else { "**NO**" }
    );
}
