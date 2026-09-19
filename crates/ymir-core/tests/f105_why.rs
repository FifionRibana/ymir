//! ADR 0001 Finding 105 — WHY the hashes differ. **Attribution before remedy.**
//!
//! The promotion gate found `k_s` bit-identical between the bench and the production helper, yet
//! the fields differ on more than half the land. So the defect is not in the function, it is in
//! the **stage at which it is called**.
//!
//! ⚠️ The hypothesis, stated before the measurement: the bench measures `k_s` on the **delivered
//! field**, which is the FULL pipeline output — `apply_bathymetry_profile` runs at
//! `production_upscale.rs:506`, AFTER the incision. Production's two-pass measures it on pass 1's
//! RAW output, **before** that re-map. Finding 78 measured what the re-map does: the receivers of
//! coastal channel cells go from a median of **−0.109 m in 3 041 distinct values** to a single
//! clamped **−20.000** (now −1.000 since Finding 79). Those receivers set `S` for the coastal
//! channel cells, so the two stages give different `k_s`.
//!
//! The test: rebuild the bench candidate with `k_s` taken from a **bathymetry-off** delivered
//! field. If that reproduces production, the stage is the cause.
//!
//! Run: cargo test -p ymir-core --release --test f105_why -- --ignored --nocapture

mod common;

use common::{CELL_KM, Knobs, PSEED, SEA, build_field_seed, pct, sorted};
use std::time::Instant;
use ymir_core::erosion::stream_power::RELIEF_V1_A_C_KM2;
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::production_upscale::{c1_altitude_norm_to_metres, flint_intercept};

const DOMAIN_KM: f32 = 400.0;
const FACTOR: f32 = 0.6;

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
fn f105_why() {
    let ss = SteinSteinParams::default();
    let cell_km2 = CELL_KM * CELL_KM;
    let n2m = c1_altitude_norm_to_metres(1.0, &ss) - c1_altitude_norm_to_metres(0.0, &ss);
    eprintln!("\n==========  Finding 105 . why the hashes differ  ==========");

    let t0 = Instant::now();
    let delivered = build_field_seed(Knobs::passes(2), PSEED);
    let del_nb = build_field_seed(Knobs { bathymetry_off: true, ..Knobs::passes(2) }, PSEED);
    let ks_post = flint_intercept(&delivered, &ss, SEA, cell_km2, DOMAIN_KM);
    let ks_pre = flint_intercept(&del_nb, &ss, SEA, cell_km2, DOMAIN_KM);
    // ADR Finding 105 -- the THIRD stage: the incision alone, before droplets AND bathymetry,
    // which is what production pass 1 actually hands to the calibration.
    let inc_only = build_field_seed(
        Knobs { erosion_off: true, bathymetry_off: true, ..Knobs::passes(2) },
        PSEED,
    );
    let ks_inc = flint_intercept(&inc_only, &ss, SEA, cell_km2, DOMAIN_KM);
    eprintln!("   k_s on the INCISION-ONLY field (no droplets, no bathymetry): **{ks_inc:.6}**");
    eprintln!(
        "   k_s on the POST-bathymetry delivered field (what the bench used): **{ks_post:.6}**\n   \
         k_s on the PRE-bathymetry field (what production's pass 1 sees):   **{ks_pre:.6}**\n   \
         ⇒ ratio **{:.4}**, floor difference **{:.6}**",
        ks_pre / ks_post,
        (ks_pre - ks_post) * FACTOR
    );

    // the two candidates, built from the two calibrations
    let mk = |ks: f32| {
        build_field_seed(
            Knobs { slope_floor_uk: Some(ks * FACTOR), depression_floor: true, ..Knobs::passes(2) },
            PSEED,
        )
    };
    let from_post = mk(ks_post);
    let from_pre = mk(ks_pre);
    let prod =
        build_field_seed(Knobs { slope_floor_factor: Some(FACTOR), ..Knobs::passes(2) }, PSEED);

    let cmp = |a: &GridF32, b: &GridF32, label: &str| {
        let n = a.data.len();
        let d: Vec<f32> = (0..n)
            .filter(|&k| a.data[k] > SEA || b.data[k] > SEA)
            .map(|k| ((a.data[k] - b.data[k]) * n2m).abs())
            .collect();
        let s = sorted(d.clone());
        eprintln!(
            "   {label:<34} hashes {} . cells differing {} . |delta| p50 **{:.3} m** p99 {:.2} m \
             max **{:.1} m**",
            if hash(a) == hash(b) { "**MATCH**" } else { "differ" },
            (0..n).filter(|&k| a.data[k] != b.data[k]).count(),
            pct(&s, 0.50),
            pct(&s, 0.99),
            s.last().copied().unwrap_or(0.0)
        );
    };
    cmp(&from_post, &prod, "bench(k_s POST) vs production");
    cmp(&from_pre, &prod, "bench(k_s PRE) vs production");
    cmp(&from_post, &from_pre, "bench POST vs bench PRE");
    let from_inc = mk(ks_inc);
    cmp(&from_inc, &prod, "bench(k_s INCISION-ONLY) vs production");
    eprintln!("\n==========  end . {:.1} s  ==========\n", t0.elapsed().as_secs_f64());
}
