//! ADR Finding 56b — the channel-head law MUST invalidate the eroded cache, and everything
//! derived from it.
//!
//! The author's prerequisite for the visual judgement is "enable the toggle and regenerate at
//! 8192²; the `eroded` and `drainage` caches invalidate since the law touches the field". That
//! is a claim about a content-addressed key, and if it were false the regeneration would return
//! the SHIPPED terrain from cache and the author would judge the wrong picture — silently, with
//! nothing in the log to say so. It is precisely the failure `eroded_key_full`'s own doc comment
//! warns about for volcanism ("the terrain differs, the drainage would not").
//!
//! So it is asserted rather than reasoned about. The law rides inside
//! `FbmUpscaleConfig::stream_power`, which `eroded_key` folds whole via
//! `.with("upscale", upscale_cfg)`; `a_c_slope_law` carries no `skip_serializing_if`, so `None`
//! serialises as `null` and the digest differs. That is the mechanism — this test is the proof,
//! and it will fail if anyone adds a skip attribute or replaces the whole-config fold.
//!
//! Run: `cargo test -p ymir-core --test channel_head_law_cache_key` (fast, always on).

use ymir_core::cache::CacheKey;
use ymir_core::erosion::stream_power::{
    CHANNEL_HEAD_S_MIN, CHANNEL_HEAD_S_REF, ChannelHeadLaw, StreamPowerConfig,
};
use ymir_core::tectonics_c1::cached_product::{drainage_key_windowed, eroded_key_full};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::VolcanismConfig;
use ymir_core::tectonics_c1::drainage::C1DrainageConfig;
use ymir_core::terrain::upscale::FbmUpscaleConfig;

fn cfg(law: bool) -> FbmUpscaleConfig {
    let mut c = FbmUpscaleConfig::default();
    let mut sp = c.stream_power.unwrap_or_else(StreamPowerConfig::default);
    sp.a_c_slope_law =
        law.then_some(ChannelHeadLaw { s_ref: CHANNEL_HEAD_S_REF, s_min: CHANNEL_HEAD_S_MIN });
    c.stream_power = Some(sp);
    c
}

fn keys(law: bool) -> (String, String) {
    let ss = SteinSteinParams::default();
    let volc = VolcanismConfig::default();
    let tectonic = CacheKey::root().with("seed", &10_481_999_410_520_546_993u64);
    let ek = eroded_key_full(&tectonic, &ss, &cfg(law), &volc);
    let dk = drainage_key_windowed(&ek, &C1DrainageConfig::default(), &ss, None, 400.0);
    (ek.digest().to_string(), dk.digest().to_string())
}

/// The law ON and the law OFF must not share an eroded digest, nor a drainage digest.
#[test]
fn the_channel_head_law_invalidates_the_eroded_and_drainage_caches() {
    let (off_e, off_d) = keys(false);
    let (on_e, on_d) = keys(true);
    assert_ne!(
        off_e, on_e,
        "the channel-head law does NOT change the eroded digest — enabling the toggle would \
         return the SHIPPED terrain from cache and the author would judge the wrong picture"
    );
    assert_ne!(
        off_d, on_d,
        "the law changes the terrain but not the drainage digest — the terrain would differ \
         and the drainage would not (the exact defect `eroded_key_full` documents for volcanism)"
    );
    eprintln!("eroded  off {} -> on {}", &off_e[..12], &on_e[..12]);
    eprintln!("drainage off {} -> on {}", &off_d[..12], &on_d[..12]);
}

/// NEGATIVE CONTROL (rule 1): the same call twice must give the SAME digests. Without it, a key
/// that folded a timestamp or a hash-ordered map would pass the test above for the wrong reason
/// — every pair of keys would differ, including two identical configurations.
#[test]
fn the_key_is_stable_across_identical_calls() {
    assert_eq!(keys(false), keys(false), "the eroded/drainage key is not deterministic");
    assert_eq!(keys(true), keys(true), "the eroded/drainage key is not deterministic");
}

/// And the two parameters of the law must BOTH reach the key: a recalibration of `S_ref` or of
/// the `S_min` proxy changes the terrain, so it must not be served from a stale entry.
#[test]
fn recalibrating_either_law_parameter_invalidates_the_cache() {
    let ss = SteinSteinParams::default();
    let volc = VolcanismConfig::default();
    let tectonic = CacheKey::root().with("seed", &1u64);
    let digest = |law: ChannelHeadLaw| {
        let mut c = FbmUpscaleConfig::default();
        let mut sp = c.stream_power.unwrap_or_else(StreamPowerConfig::default);
        sp.a_c_slope_law = Some(law);
        c.stream_power = Some(sp);
        eroded_key_full(&tectonic, &ss, &c, &volc).digest().to_string()
    };
    let base = || ChannelHeadLaw { s_ref: CHANNEL_HEAD_S_REF, s_min: CHANNEL_HEAD_S_MIN };
    assert_ne!(
        digest(base()),
        digest(ChannelHeadLaw { s_ref: 0.1128, ..base() }),
        "`s_ref` does not reach the cache key — the 8192² recalibration would be served stale"
    );
    assert_ne!(
        digest(base()),
        digest(ChannelHeadLaw { s_min: 0.0349, ..base() }),
        "`s_min` does not reach the cache key — replacing the PROXY floor would be served stale"
    );
}
