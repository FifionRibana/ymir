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
//! * **The lakes** (Finding 134) are guarded apart, after the C-2 crater pass, under the final
//!   drainage's digest: [`check_lakes`], reference `data/bench_lake_hashes.json`.
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

// ═══════════════════════════════════════════════════════════════════════════════════════════════
// ADR 0001 Finding 134 — the guard covers the LAKES. The field guard proves the eroded field is the
// bench's; Finding 133 showed the lakes the viz lists may still not be: the viz path listed 26 / 34
// against the bench assembly's 25 / 33 (the C-2 crater lake, added after the assembly). So the lakes
// get their own entry, keyed by the FINAL drainage's cache digest (`hd_drainage_key`: the eroded key
// + the drainage config + the climate) and fingerprinted AFTER the crater pass, on the lake mask and
// on the lake list. One climate is guarded (latitude 45°, span 40°): any other reads NoReference.
// ═══════════════════════════════════════════════════════════════════════════════════════════════

/// One bench-measured lake population.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LakeGuardEntry {
    /// `hd_drainage_key(..).digest()` of the bench's configuration (eroded key + drainage + climate).
    pub digest: String,
    pub label: String,
    /// [`lake_fingerprint`] of the bench's final lakes, `"<lake_map hash>/<lake list hash>"`.
    pub lakes_hash: String,
    /// The number of lakes, for the reader (the hash decides).
    pub lakes: usize,
    pub latitude_deg: f32,
    pub span_deg: f32,
}

const RAW_LAKES: &str = include_str!("../../data/bench_lake_hashes.json");

/// The embedded lake reference (empty if the file is empty or malformed — which reads NoReference).
pub fn lake_entries() -> Vec<LakeGuardEntry> {
    serde_json::from_str(RAW_LAKES).unwrap_or_default()
}

/// FNV-1a over the lake mask's `u32` ids (row-major), then over the list: per lake its id, its type's
/// name, its level and its area (`f32` bits), in list order.
pub fn lake_fingerprint(lake_map: &[u32], lakes: &[crate::tectonics_c1::drainage::C1Lake]) -> String {
    let fnv = |h: &mut u64, bytes: &[u8]| {
        for &b in bytes {
            *h ^= b as u64;
            *h = h.wrapping_mul(0x1000_0000_01b3);
        }
    };
    let mut hm = 0xcbf2_9ce4_8422_2325u64;
    for &v in lake_map {
        fnv(&mut hm, &v.to_le_bytes());
    }
    let mut hl = 0xcbf2_9ce4_8422_2325u64;
    for l in lakes {
        fnv(&mut hl, &l.base.id.to_le_bytes());
        fnv(&mut hl, format!("{:?}", l.lake_type).as_bytes());
        fnv(&mut hl, &l.level_m.to_bits().to_le_bytes());
        fnv(&mut hl, &l.area_km2.to_bits().to_le_bytes());
    }
    format!("{hm:016x}/{hl:016x}")
}

/// One line per lake, in list order: id, type, level, area (their `f32` bits beside), its cell count
/// in the mask and the FNV-1a of its cells' indices (row-major). The bench and the viz print the SAME
/// lines, so a lake-by-lake identity is a text diff (Finding 134, item 0.3).
pub fn lake_listing(lake_map: &[u32], lakes: &[crate::tectonics_c1::drainage::C1Lake]) -> Vec<String> {
    use std::collections::HashMap;
    let mut per: HashMap<u32, (usize, u64)> = HashMap::new();
    for (i, &v) in lake_map.iter().enumerate() {
        if v != 0 {
            let e = per.entry(v).or_insert((0, 0xcbf2_9ce4_8422_2325));
            e.0 += 1;
            for b in (i as u64).to_le_bytes() {
                e.1 ^= b as u64;
                e.1 = e.1.wrapping_mul(0x1000_0000_01b3);
            }
        }
    }
    let mut out: Vec<String> = lakes
        .iter()
        .map(|l| {
            let (c, hc) = per.remove(&l.base.id).unwrap_or((0, 0));
            format!(
                "lake {:>8} {:<13} level {:>9.2} m ({:08x}) area {:>9.3} km² ({:08x}) cells {:>8} mask {:016x}",
                l.base.id,
                format!("{:?}", l.lake_type),
                l.level_m,
                l.level_m.to_bits(),
                l.area_km2,
                l.area_km2.to_bits(),
                c,
                hc
            )
        })
        .collect();
    let mut orphans: Vec<_> = per.into_iter().collect();
    orphans.sort_by_key(|e| e.0);
    for (id, (c, hc)) in orphans {
        out.push(format!("mask id {id:>8} with NO list entry · cells {c:>8} mask {hc:016x}"));
    }
    out
}

/// Compare a rendered lake population with the embedded lake reference.
pub fn check_lakes(
    digest: &str,
    lake_map: &[u32],
    lakes: &[crate::tectonics_c1::drainage::C1Lake],
) -> GuardStatus {
    check_lakes_against(&lake_entries(), digest, lake_map, lakes)
}

/// [`check_lakes`] against an explicit reference.
pub fn check_lakes_against(
    reference: &[LakeGuardEntry],
    digest: &str,
    lake_map: &[u32],
    lakes: &[crate::tectonics_c1::drainage::C1Lake],
) -> GuardStatus {
    let viz = lake_fingerprint(lake_map, lakes);
    match reference.iter().find(|e| e.digest == digest) {
        Some(e) if e.lakes_hash == viz => GuardStatus::Match { label: e.label.clone() },
        Some(e) => {
            GuardStatus::Mismatch { label: e.label.clone(), bench: e.lakes_hash.clone(), viz }
        }
        None => GuardStatus::NoReference { viz },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// ADR Finding 134, rule 13 — the lake guard's verdicts, each with its negative control: ONE mask
    /// cell, ONE level, ONE type, ONE area changed must each read Mismatch; another digest (another
    /// climate) reads NoReference.
    #[test]
    fn the_lake_guard_sees_one_cell_one_level_one_type_one_area() {
        use crate::lakes::detection::Lake;
        use crate::tectonics_c1::drainage::{C1Lake, LakeType};
        let lake = C1Lake {
            base: Lake {
                id: 7,
                surface_elevation: 0.6,
                max_depth: 0.01,
                area: 2,
                basin_id: 0,
                outlet: (0, 0),
                shallow: false,
            },
            level_m: 362.0,
            depth_m: 80.0,
            area_km2: 0.004,
            lake_type: LakeType::Exorheic,
            unresolved_reason: None,
        };
        let map = vec![0u32, 7, 7, 0];
        let lakes = vec![lake.clone()];
        let e = LakeGuardEntry {
            digest: "d".into(),
            label: "x".into(),
            lakes_hash: lake_fingerprint(&map, &lakes),
            lakes: 1,
            latitude_deg: 45.0,
            span_deg: 40.0,
        };
        let r = [e];
        let verdict = |digest: &str, m: &[u32], l: &[C1Lake]| check_lakes_against(&r, digest, m, l);
        assert_eq!(verdict("d", &map, &lakes), GuardStatus::Match { label: "x".into() });
        // negative controls
        let mut m1 = map.clone();
        m1[3] = 7;
        assert!(verdict("d", &m1, &lakes).refuses_numbers(), "one mask cell");
        let mut l1 = lakes.clone();
        l1[0].level_m = 362.01;
        assert!(verdict("d", &map, &l1).refuses_numbers(), "one level");
        let mut l2 = lakes.clone();
        l2[0].lake_type = LakeType::CraterAcidic;
        assert!(verdict("d", &map, &l2).refuses_numbers(), "one type");
        let mut l3 = lakes.clone();
        l3[0].area_km2 = 0.005;
        assert!(verdict("d", &map, &l3).refuses_numbers(), "one area");
        assert!(verdict("d", &map, &[]).refuses_numbers(), "one lake fewer");
        assert!(matches!(verdict("other climate", &map, &lakes), GuardStatus::NoReference { .. }));
        // the embedded reference parses (an empty array is valid)
        let _ = lake_entries();
    }

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
