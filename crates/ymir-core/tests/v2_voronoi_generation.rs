//! Step 6 — Voronoi tessellation: determinism + counts + edge cases.

use ymir_core::tectonics_v2::voronoi::{VoronoiConfig, generate_voronoi};

#[test]
fn determinism_byte_for_byte() {
    let cfg = VoronoiConfig::default();
    let a = generate_voronoi(64, 64, &cfg, 42);
    let b = generate_voronoi(64, 64, &cfg, 42);
    assert_eq!(a.plate_id.data(), b.plate_id.data());
    assert_eq!(a.per_plate_type, b.per_plate_type);
    for (p, q) in a.seed_coords.iter().zip(b.seed_coords.iter()) {
        assert_eq!(p, q);
    }
}

#[test]
fn distinct_plate_count_matches_num_plates() {
    let cfg = VoronoiConfig { num_plates: 8, continental_ratio: 0.3 };
    let v = generate_voronoi(64, 64, &cfg, 42);
    let mut seen = std::collections::HashSet::new();
    for &id in v.plate_id.data() {
        seen.insert(id);
    }
    assert_eq!(seen.len(), 8);
}

#[test]
fn num_plates_equals_one_covers_everything() {
    let cfg = VoronoiConfig { num_plates: 1, continental_ratio: 0.5 };
    let v = generate_voronoi(16, 16, &cfg, 7);
    for &id in v.plate_id.data() {
        assert_eq!(id, 0);
    }
}

#[test]
fn different_seeds_produce_different_tessellations() {
    let cfg = VoronoiConfig::default();
    let a = generate_voronoi(64, 64, &cfg, 42);
    let b = generate_voronoi(64, 64, &cfg, 43);
    assert_ne!(a.plate_id.data(), b.plate_id.data());
}
