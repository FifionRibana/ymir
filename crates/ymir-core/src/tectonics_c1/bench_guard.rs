//! ADR 0001 Finding 123 — **the viz / bench identity guard.**
//!
//! Finding 121 declared the guard absent ("the viz run itself was NOT guarded bit-for-bit against
//! the bench") and Finding 123 paid for it: the author's microscope showed a "654 km" trunk that no
//! bench world contains. From here, a world the viz renders is compared to the world a bench
//! measured AT THE SAME SETTINGS before a number of the viz may be read against a Finding.
//!
//! * **The settings** are the eroded product's cache key, [`eroded_key_full`]'s digest — the SAME
//!   function the viz already calls for its cache, so equal digests mean equal configurations
//!   (seed, tectonic run, upscale config incl. the framing origin, volcanism).
//! * **The world** is the FNV-1a hash of the eroded heightmap's `f32` bits ([`field_hash`], the
//!   benches' `hash()` since Finding 105).
//! * **The reference** is `crates/ymir-core/data/bench_field_hashes.json`, written by the
//!   `f123_guard` bench and embedded at compile time: versioned, reviewable, never a cache.
//!
//! [`eroded_key_full`]: crate::tectonics_c1::cached_product::eroded_key_full

use crate::grid::GridF32;
use serde::{Deserialize, Serialize};

/// One bench-measured world.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GuardEntry {
    /// `eroded_key_full(..).digest()` of the bench's configuration.
    pub digest: String,
    /// The state's name in the Findings ("livré", "A1+B2", "C1 nue", "C2 /10", "C2 /10 col").
    pub label: String,
    /// [`field_hash`] of the bench's eroded heightmap, hex.
    pub field_hash: String,
    pub seed: u64,
    pub target: usize,
    pub origin: [f64; 2],
}

const RAW: &str = include_str!("../../data/bench_field_hashes.json");

/// The embedded reference (empty if the file is empty or malformed — which reads NoReference).
pub fn entries() -> Vec<GuardEntry> {
    serde_json::from_str(RAW).unwrap_or_default()
}

/// FNV-1a over the raw `f32` bits — the benches' hash since Finding 105.
pub fn field_hash(g: &GridF32) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for v in &g.data {
        for b in v.to_bits().to_le_bytes() {
            h ^= b as u64;
            h = h.wrapping_mul(0x1000_0000_01b3);
        }
    }
    h
}

/// The guard's verdict for one rendered world.
#[derive(Clone, Debug, PartialEq)]
pub enum GuardStatus {
    /// A bench measured this very configuration and the field is bit-identical.
    Match { label: String },
    /// A bench measured this configuration and the field DIFFERS: no number of this world may be
    /// read against that bench's Finding.
    Mismatch { label: String, bench: String, viz: String },
    /// No bench measured this configuration: the numbers are this world's, unguarded.
    NoReference { viz: String },
}

impl GuardStatus {
    /// `true` when the viz must refuse to display a number.
    pub fn refuses_numbers(&self) -> bool {
        matches!(self, GuardStatus::Mismatch { .. })
    }
}

/// Compare a rendered eroded field with the embedded reference.
pub fn check(digest: &str, eroded: &GridF32) -> GuardStatus {
    check_against(&entries(), digest, eroded)
}

/// [`check`] against an explicit reference.
pub fn check_against(reference: &[GuardEntry], digest: &str, eroded: &GridF32) -> GuardStatus {
    let viz = format!("{:016x}", field_hash(eroded));
    match reference.iter().find(|e| e.digest == digest) {
        Some(e) if e.field_hash == viz => GuardStatus::Match { label: e.label.clone() },
        Some(e) => {
            GuardStatus::Mismatch { label: e.label.clone(), bench: e.field_hash.clone(), viz }
        }
        None => GuardStatus::NoReference { viz },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// ADR Finding 123, rule 13 — the three verdicts, on a reference built in the test.
    #[test]
    fn the_guard_reads_match_mismatch_and_no_reference() {
        let g = GridF32 { width: 2, height: 2, data: vec![0.1, 0.2, 0.3, 0.4] };
        let h = format!("{:016x}", field_hash(&g));
        let mut other = g.clone();
        other.data[3] = 0.41;
        assert_ne!(h, format!("{:016x}", field_hash(&other)), "the hash must see one cell");
        let e = GuardEntry {
            digest: "d".into(),
            label: "x".into(),
            field_hash: h.clone(),
            seed: 1,
            target: 2,
            origin: [0.0, 0.0],
        };
        let r = [e];
        let verdict = |digest: &str, f: &GridF32| check_against(&r, digest, f);
        assert_eq!(verdict("d", &g), GuardStatus::Match { label: "x".into() });
        assert!(verdict("d", &other).refuses_numbers());
        assert!(matches!(verdict("zz", &g), GuardStatus::NoReference { .. }));
        // the embedded reference parses (an empty array is valid)
        let _ = entries();
    }
}
