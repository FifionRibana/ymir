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
use ymir_core::tectonics_c1::drainage::{
    C1DrainageConfig, MergedBasinOutlet, MergedUnionRelevel, ReciprocalRule,
};
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

// ─────────────────────────────────────────────────────────────────────────────
// ADR 0001 Finding 85 — the same two questions for `merged_union_relevel`.
// ─────────────────────────────────────────────────────────────────────────────

/// `true` = the shipped default (Finding 86 B2), `false` = the A/B control.
fn cfg85(gate: bool) -> C1DrainageConfig {
    C1DrainageConfig {
        merged_union_relevel: gate.then(MergedUnionRelevel::default),
        ..Default::default()
    }
}

#[test]
fn the_relevel_gate_reaches_both_drainage_keys() {
    let ss = SteinSteinParams::default();
    let tectonic = CacheKey::root().with("seed", &10_481_999_410_520_546_993u64);
    let ek =
        eroded_key_full(&tectonic, &ss, &FbmUpscaleConfig::default(), &VolcanismConfig::default());
    let key = |gate: bool| {
        let dk = drainage_key_windowed(&ek, &cfg85(gate), &ss, None, 400.0);
        let hk = hd_drainage_key(
            &ek,
            &cfg85(gate),
            &ss,
            45.0,
            &PrecipParams::default(),
            Some(40.0),
            7.5,
            400.0,
        );
        (dk.digest().to_string(), hk.digest().to_string())
    };
    let (d_off, h_off) = key(false);
    let (d_on, h_on) = key(true);
    eprintln!(
        "   windowed {d_off} -> {d_on}
   hd       {h_off} -> {h_on}"
    );
    assert_ne!(d_off, d_on, "the windowed drainage key does not see `merged_union_relevel`");
    assert_ne!(h_off, h_on, "the HD drainage key does not see it: lakes.json would be stale");
    // NEGATIVE CONTROL: without these the test passes on a key that folds a timestamp.
    assert_eq!(key(false).0, d_off, "the drainage key is not stable across identical calls");
    assert_eq!(key(true).1, h_on, "the HD drainage key is not stable across identical calls");
    // and the two gates of Findings 84 and 85 must not collide into one digest
    let both = C1DrainageConfig {
        merged_basin_outlet: Some(MergedBasinOutlet { rule: ReciprocalRule::LevelThenMerge }),
        merged_union_relevel: Some(MergedUnionRelevel::default()),
        ..Default::default()
    };
    let d_both = drainage_key_windowed(&ek, &both, &ss, None, 400.0).digest().to_string();
    assert_ne!(d_both, d_on, "the Finding 84 and 85 gates are indistinguishable in the key");
    assert_ne!(d_both, d_off);
}

/// ⚠️ REWRITTEN at ADR Finding 86 B2: the relevel is PROMOTED, so `Default` now carries it and
/// **must** serialise it — the opposite of what this test asserted one round ago. What the round
/// buys instead is that the A/B control (`None`, the pre-Finding-86 world) is distinguishable in
/// the key, and that `max_passes` still reaches it.
#[test]
fn relevel_is_the_shipped_default_and_the_control_is_distinguishable() {
    let json = serde_json::to_string(&C1DrainageConfig::default()).expect("serialise");
    assert!(
        json.contains("merged_union_relevel"),
        "the promoted relevel is not serialised: cached products would not see it
{json}"
    );
    assert!(
        C1DrainageConfig::default().merged_union_relevel.is_some(),
        "the relevel is not the shipped default"
    );
    // and `max_passes` reaches the key, or a re-bound is served stale
    let ss = SteinSteinParams::default();
    let tectonic = CacheKey::root().with("seed", &10_481_999_410_520_546_993u64);
    let ek =
        eroded_key_full(&tectonic, &ss, &FbmUpscaleConfig::default(), &VolcanismConfig::default());
    let with = |p: u32| {
        let c = C1DrainageConfig {
            merged_union_relevel: Some(MergedUnionRelevel { max_passes: p, ..Default::default() }),
            ..Default::default()
        };
        drainage_key_windowed(&ek, &c, &ss, None, 400.0).digest().to_string()
    };
    assert_ne!(with(16), with(8), "`max_passes` does not reach the key");
    // the A/B control -- explicit `None` -- must be distinguishable from the shipped default
    let ctrl = C1DrainageConfig { merged_union_relevel: None, ..Default::default() };
    assert_ne!(
        drainage_key_windowed(&ek, &ctrl, &ss, None, 400.0).digest().to_string(),
        with(16),
        "the pre-Finding-86 A/B control shares a cache key with the shipped default"
    );
}
