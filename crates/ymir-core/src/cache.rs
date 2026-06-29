//! Content-addressed cache for generation steps — skip recomputing an upstream
//! step whose inputs are unchanged.
//!
//! **What it accelerates:** downstream calibration on a FIXED relief (climate /
//! biomes / latitude / hypsometry without re-running erosion). It does NOT
//! accelerate the relief loop — changing the relief changes the key, so erosion
//! re-runs (that is invalidation working correctly, not a miss to fix).
//!
//! **Why it cannot serve stale results (for what is in the key):** the cache is
//! *content-addressed*. The on-disk filename embeds `blake3(canonical_json(key))`.
//! A different key → a different digest → a different filename → a guaranteed
//! MISS. There is no "compare then decide" step that could go wrong: the
//! filesystem index IS the validity check.
//!
//! **The one human weak point:** the key hashes the config + an explicit
//! `ALGO_*` version per step, NOT the source code. A code change with no config
//! change is invisible to the digest unless its `ALGO_*` is bumped. See the
//! BUMP rule on each constant below — forgetting it is the only way to serve
//! stale results.
//!
//! Codec: payloads are written via the SHARED `export::raw` codec (the same
//! bytes as the §9 export), distinct policy (disposable, content-addressed).
//! See `docs/design/c1_generation_cache.md`.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde::Serialize;
use serde_json::Value;

use crate::grid::GridF32;

// ── ALGO_* step versions ─────────────────────────────────────────────────
//
// Each cached step has a version integer folded into its cache key. The digest
// hashes (config + this version), NOT the source code, so a CODE change that
// leaves the config untouched is invisible to the cache.
//
//   ⚠️ BUMP this on ANY code change to this step — the cache key hashes config +
//   this version, NOT the source code. Forgetting to bump → cache serves stale
//   results from old code (silent bug). When in doubt, bump (a redundant
//   recompute costs minutes; a stale read validates FALSE results).

/// Version of the 64² tectonic build (`init_c1_state_phase_2_r7` +
/// `run_with_closures`). ⚠️ BUMP on ANY code change to that step (see above).
pub const ALGO_TECTONICS: u32 = 1;

/// Version of the HD altitude build (`upscale_from_c1`: isostasy, Stein-Stein,
/// FBM, hydraulic erosion). ⚠️ BUMP on ANY code change to that step — including
/// the production-altitude / continental-floor code (see above).
pub const ALGO_UPSCALE_EROSION: u32 = 1;

/// Version of the drainage build (`c1_drainage`). ⚠️ BUMP on ANY code change to
/// that step (see above).
pub const ALGO_DRAINAGE: u32 = 2;

// ── Cache key ──────────────────────────────────────────────────────────────

/// An ordered set of named inputs, hashed to a content address.
///
/// Built canonically: entries live in a `BTreeMap` (keys sorted), so the
/// serialized JSON — and hence the digest — is independent of insertion order.
/// Floats serialize to their shortest round-trippable form (deterministic).
/// `std::hash::Hash` is deliberately NOT used: the config structs contain
/// floats, which do not implement `Hash`/`Eq`.
#[derive(Debug, Clone, Default)]
pub struct CacheKey {
    entries: BTreeMap<String, Value>,
}

impl CacheKey {
    /// An empty key.
    pub fn root() -> Self {
        Self { entries: BTreeMap::new() }
    }

    /// A key chained to a parent: it embeds the parent's digest, so any upstream
    /// change propagates to every downstream digest. This materializes the
    /// dependency chain (`k_eroded ⊃ k_state`, `k_drainage ⊃ k_eroded`).
    pub fn derived_from(parent: &CacheKey) -> Self {
        Self::root().with("_parent", &parent.digest())
    }

    /// Add a named input (any `Serialize` config). Panics only if the value is
    /// not representable as JSON (e.g. a NaN/Inf float) — generation configs are
    /// finite, so this is a programming error, not a runtime condition.
    pub fn with<T: Serialize>(mut self, label: &str, value: &T) -> Self {
        let v = serde_json::to_value(value)
            .expect("cache key value must be JSON-serializable (finite floats)");
        self.entries.insert(label.to_string(), v);
        self
    }

    /// Add a named input by its `Debug` representation — for config structs that
    /// do NOT derive `Serialize`. A *derived* `Debug` walks EVERY field, so this
    /// fingerprints the whole struct (completeness by construction). ⚠️ Only use
    /// on types with a DERIVED `Debug`: a hand-written `Debug` that hides fields
    /// would silently drop them from the key (the stale-cache bug).
    pub fn with_debug<T: std::fmt::Debug>(mut self, label: &str, value: &T) -> Self {
        self.entries.insert(label.to_string(), Value::String(format!("{value:?}")));
        self
    }

    /// Fold a step's `ALGO_*` version into the key (the source-code fingerprint).
    pub fn algo(self, version: u32) -> Self {
        self.with("algo_version", &version)
    }

    /// The content address: first 16 hex chars of `blake3(canonical_json)`.
    /// 64 bits — collision-free for any realistic number of cache entries.
    pub fn digest(&self) -> String {
        let json =
            serde_json::to_vec(&self.entries).expect("cache key must serialize (finite floats)");
        let hash = blake3::hash(&json);
        hash.to_hex().as_str()[..16].to_string()
    }
}

// ── Raw payload codec ────────────────────────────────────────────────────

/// A cacheable value: knows its shape (so a reload knows the dimensions) and how
/// to read/write its raw payload via the shared `export::raw` codec.
///
/// `stem` is the path WITHOUT extension (`{dir}/{step}_{digest}`); an impl may
/// write one file (`stem.raw`) or several (`stem_field.raw`).
pub trait RawCodec: Sized {
    /// Shape metadata recorded in the sidecar (e.g. grid dims), passed back to
    /// [`RawCodec::read_raw`] on a HIT.
    fn shape(&self) -> Value;
    /// Write the raw payload(s).
    fn write_raw(&self, stem: &Path) -> Result<(), String>;
    /// Read the raw payload(s), given the recorded shape.
    fn read_raw(stem: &Path, shape: &Value) -> Result<Self, String>;
}

impl RawCodec for GridF32 {
    fn shape(&self) -> Value {
        serde_json::json!({ "width": self.width, "height": self.height })
    }
    fn write_raw(&self, stem: &Path) -> Result<(), String> {
        self.save_raw(&stem.with_extension("raw"))
    }
    fn read_raw(stem: &Path, shape: &Value) -> Result<Self, String> {
        let w = shape["width"].as_u64().ok_or("sidecar: missing width")? as usize;
        let h = shape["height"].as_u64().ok_or("sidecar: missing height")? as usize;
        GridF32::load_raw(&stem.with_extension("raw"), w, h)
    }
}

// ── The wrapper ────────────────────────────────────────────────────────────

#[derive(Serialize, serde::Deserialize)]
struct Sidecar {
    step: String,
    digest: String,
    /// The readable key (debug / provenance) — NOT used for validity (the
    /// filename digest is). Kept so a human can see what produced the file.
    key: BTreeMap<String, Value>,
    shape: Value,
}

/// Default cache directory: `.ymir_cache/` under the current working dir.
/// Internal and disposable — add it to `.gitignore`; deleting it is always safe
/// (content-addressing stays correct regardless of what is or isn't present).
pub fn default_cache_dir() -> PathBuf {
    PathBuf::from(".ymir_cache")
}

/// Read-if-present-valid / compute-and-write-else, content-addressed.
///
/// - HIT: the sidecar `{dir}/{step}_{digest}.json` exists → load the payload and
///   return it (no `compute` call).
/// - MISS: `compute()` runs, the payload + sidecar are written, the value is
///   returned.
///
/// Stale reads are STRUCTURALLY impossible for anything in `key`: a changed key
/// yields a different digest, hence a different filename, hence a MISS. The
/// generation function (`compute`) is untouched and stays pure — the cache only
/// wraps it.
pub fn cached<T: RawCodec>(
    dir: &Path,
    step: &str,
    key: &CacheKey,
    compute: impl FnOnce() -> T,
) -> Result<T, String> {
    let digest = key.digest();
    let stem = dir.join(format!("{step}_{digest}"));
    let sidecar_path = stem.with_extension("json");

    if sidecar_path.exists() {
        let json = std::fs::read_to_string(&sidecar_path)
            .map_err(|e| format!("cache sidecar read error: {e}"))?;
        let meta: Sidecar =
            serde_json::from_str(&json).map_err(|e| format!("cache sidecar parse error: {e}"))?;
        return T::read_raw(&stem, &meta.shape); // HIT
    }

    // MISS — compute, then persist.
    let value = compute();
    std::fs::create_dir_all(dir).map_err(|e| format!("cache dir create error: {e}"))?;
    value.write_raw(&stem)?;
    let meta = Sidecar {
        step: step.to_string(),
        digest: digest.clone(),
        key: key.entries.clone(),
        shape: value.shape(),
    };
    let json =
        serde_json::to_string_pretty(&meta).map_err(|e| format!("sidecar json error: {e}"))?;
    std::fs::write(&sidecar_path, json).map_err(|e| format!("sidecar write error: {e}"))?;
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

    fn tmp(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("ymir_cache_test_{name}"));
        let _ = std::fs::remove_dir_all(&dir); // fresh each run
        dir
    }

    fn sample_grid() -> GridF32 {
        let mut g = GridF32::new(8, 6, 0.0);
        for j in 0..6 {
            for i in 0..8 {
                g.set(i, j, (i as f32) * 0.5 - (j as f32) * 0.25 + 0.1);
            }
        }
        g
    }

    #[test]
    fn digest_is_deterministic_and_order_independent() {
        let a = CacheKey::root().with("seed", &42u64).with("floor", &450.0f32).algo(1);
        let b = CacheKey::root().algo(1).with("floor", &450.0f32).with("seed", &42u64);
        assert_eq!(a.digest(), b.digest(), "digest must be insertion-order independent");
        let c = CacheKey::root().with("seed", &43u64).with("floor", &450.0f32).algo(1);
        assert_ne!(a.digest(), c.digest(), "different input → different digest");
    }

    #[test]
    fn miss_then_hit() {
        let dir = tmp("miss_then_hit");
        let key = CacheKey::root().with("seed", &7u64).algo(ALGO_UPSCALE_EROSION);
        let calls = Cell::new(0u32);

        // First call: MISS — compute runs, files written.
        let g1 = cached(&dir, "eroded", &key, || {
            calls.set(calls.get() + 1);
            sample_grid()
        })
        .unwrap();
        assert_eq!(calls.get(), 1, "first call must compute");
        assert!(dir.join(format!("eroded_{}.raw", key.digest())).exists());
        assert!(dir.join(format!("eroded_{}.json", key.digest())).exists());

        // Second call, same key: HIT — compute must NOT run, value loaded.
        let g2 = cached(&dir, "eroded", &key, || {
            calls.set(calls.get() + 1);
            sample_grid()
        })
        .unwrap();
        assert_eq!(calls.get(), 1, "second call with same key must HIT (no recompute)");

        // Round-trip is byte-exact.
        assert_eq!(g1.width, g2.width);
        assert_eq!(g1.height, g2.height);
        assert_eq!(g1.data, g2.data, "cached field must reload byte-identical");
    }

    #[test]
    fn invalidation_by_config_change() {
        // The relief-fix case: changing `continental_floor_m` changes the digest
        // → MISS → recompute (erosion would re-run on the new relief).
        let dir = tmp("invalidation");
        let calls = Cell::new(0u32);
        let mut compute = |floor: f32| {
            let key = CacheKey::root()
                .with("seed", &7u64)
                .with("continental_floor_m", &floor)
                .algo(ALGO_UPSCALE_EROSION);
            cached(&dir, "eroded", &key, || {
                calls.set(calls.get() + 1);
                sample_grid()
            })
            .unwrap();
            key.digest()
        };

        let d_a = compute(450.0);
        assert_eq!(calls.get(), 1);
        let d_a2 = compute(450.0); // same floor → HIT
        assert_eq!(calls.get(), 1, "unchanged floor must HIT");
        assert_eq!(d_a, d_a2);

        let d_b = compute(200.0); // changed floor → MISS
        assert_eq!(calls.get(), 2, "changed floor must MISS (invalidate)");
        assert_ne!(d_a, d_b, "changed config must change the digest");
    }

    #[test]
    fn algo_version_bump_invalidates() {
        // A code change with no config change: the ONLY protection is bumping
        // ALGO_*. Same config, different algo version → different digest → MISS.
        let dir = tmp("algo_bump");
        let calls = Cell::new(0u32);
        let mut run = |algo: u32| {
            let key = CacheKey::root().with("seed", &7u64).algo(algo);
            cached(&dir, "eroded", &key, || {
                calls.set(calls.get() + 1);
                sample_grid()
            })
            .unwrap();
        };
        run(1);
        run(1);
        assert_eq!(calls.get(), 1, "same algo version → HIT");
        run(2);
        assert_eq!(calls.get(), 2, "bumped algo version → MISS");
    }

    #[test]
    fn chaining_propagates_upstream_change() {
        // k_eroded ⊃ k_state: an upstream change flips the downstream digest.
        let k_state_a = CacheKey::root().with("seed", &7u64).algo(ALGO_TECTONICS);
        let k_state_b = CacheKey::root().with("seed", &8u64).algo(ALGO_TECTONICS);
        let k_eroded_a = CacheKey::derived_from(&k_state_a).algo(ALGO_UPSCALE_EROSION);
        let k_eroded_b = CacheKey::derived_from(&k_state_b).algo(ALGO_UPSCALE_EROSION);
        assert_ne!(
            k_eroded_a.digest(),
            k_eroded_b.digest(),
            "upstream (state) change must change the downstream (eroded) digest"
        );
    }
}
