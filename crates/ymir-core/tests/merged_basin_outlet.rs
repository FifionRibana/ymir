//! ADR 0001 Finding 84 — the gate tests for `C1DrainageConfig::merged_basin_outlet`.
//!
//! The behavioural tests live next to the other below-sea fixtures, in
//! `tectonics_c1::drainage`'s own `mod tests` (they need the private flood to compare against).
//! What has to live out here is the CACHE KEY: a seam that changes the drainage but not the key
//! serves a stale `lakes.json` and `rivers.json` from disk, and every measurement after it is a
//! measurement of the old world.
//!
//! Run: cargo test -p ymir-core --release --test merged_basin_outlet

use ymir_core::cache::CacheKey;
use ymir_core::climate::precipitation::PrecipParams;
use ymir_core::tectonics_c1::cached_product::{
    drainage_key_windowed, eroded_key_full, hd_drainage_key,
};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::VolcanismConfig;
use ymir_core::tectonics_c1::drainage::{C1DrainageConfig, MergedBasinOutlet, ReciprocalRule};
use ymir_core::terrain::upscale::FbmUpscaleConfig;

fn cfg(gate: bool) -> C1DrainageConfig {
    C1DrainageConfig {
        merged_basin_outlet: gate
            .then_some(MergedBasinOutlet { rule: ReciprocalRule::LevelThenMerge }),
        ..Default::default()
    }
}

/// `(windowed drainage key, HD drainage key)` for a gate state.
fn keys(gate: bool) -> (String, String) {
    let ss = SteinSteinParams::default();
    let tectonic = CacheKey::root().with("seed", &10_481_999_410_520_546_993u64);
    let ek =
        eroded_key_full(&tectonic, &ss, &FbmUpscaleConfig::default(), &VolcanismConfig::default());
    let dk = drainage_key_windowed(&ek, &cfg(gate), &ss, None, 400.0);
    let hk = hd_drainage_key(
        &ek,
        &cfg(gate),
        &ss,
        45.0,
        &PrecipParams::default(),
        Some(40.0),
        7.5,
        400.0,
    );
    (dk.digest().to_string(), hk.digest().to_string())
}

#[test]
fn the_gate_reaches_both_drainage_keys() {
    let (d_off, h_off) = keys(false);
    let (d_on, h_on) = keys(true);
    eprintln!("   windowed {d_off} -> {d_on}\n   hd       {h_off} -> {h_on}");
    assert_ne!(
        d_off, d_on,
        "the windowed drainage key does not see `merged_basin_outlet`: a stale lake set ships"
    );
    assert_ne!(
        h_off, h_on,
        "the HD drainage key does not see it: lakes.json and rivers.json would be served from cache"
    );
    // NEGATIVE CONTROL — without these two, the test above passes on a key that folds a timestamp
    // or a pointer, and "the key moved" would mean nothing.
    assert_eq!(keys(false).0, d_off, "the drainage key is not stable across identical calls");
    assert_eq!(keys(true).1, h_on, "the HD drainage key is not stable across identical calls");
}

/// The gate is `#[serde(default, skip_serializing_if = "Option::is_none")]`, so `None` must leave
/// the key EXACTLY where it was before the field existed. Without this the promotion of any later
/// finding would invalidate every cached product for nothing.
#[test]
fn none_leaves_the_key_where_it_was() {
    let ss = SteinSteinParams::default();
    let tectonic = CacheKey::root().with("seed", &10_481_999_410_520_546_993u64);
    let ek =
        eroded_key_full(&tectonic, &ss, &FbmUpscaleConfig::default(), &VolcanismConfig::default());
    // `Default` and an explicitly-`None` config must serialise to the same bytes.
    let a = drainage_key_windowed(&ek, &C1DrainageConfig::default(), &ss, None, 400.0);
    let b = drainage_key_windowed(&ek, &cfg(false), &ss, None, 400.0);
    assert_eq!(a.digest(), b.digest(), "`None` is not the pre-Finding-84 serialisation");
    // and the JSON really does not carry the field when it is None
    let json = serde_json::to_string(&C1DrainageConfig::default()).expect("serialise");
    assert!(
        !json.contains("merged_basin_outlet"),
        "`None` is serialised: every cached drainage product would be invalidated for nothing\n{json}"
    );
}
