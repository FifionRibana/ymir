//! Step 8.6 Phase 4 acceptance test #2.
//!
//! Verifies every embedded preset deserializes into a valid
//! `V2RunSpec` and that the spec round-trips through `build_config`
//! without panicking (i.e. it produces a well-formed `BaselineConfig`
//! the harness would accept). Headless — no Bevy, no run_baseline.

use ymir_viz::bridge_v2::{build_config, presets};

#[test]
fn v2_bridge_preset_load() {
    let names = presets::list();
    assert!(!names.is_empty(), "no embedded presets registered");
    assert!(
        names.iter().any(|n| *n == "active_medley"),
        "active_medley preset must always be present (D4 baseline reference)"
    );

    for name in names {
        let spec = presets::load(name).unwrap_or_else(|e| {
            panic!("preset '{}' failed to deserialize: {}", name, e)
        });
        // Sanity bounds — catches preset typos that would otherwise
        // surface only at run time when build_config is invoked.
        assert!(spec.grid_nx > 0 && spec.grid_nx <= 256, "{} grid_nx", name);
        assert!(spec.grid_ny > 0 && spec.grid_ny <= 256, "{} grid_ny", name);
        assert!(spec.steps > 0 && spec.steps <= 10_000, "{} steps", name);
        assert!(spec.num_plates >= 2, "{} num_plates", name);
        assert!(
            (0.0..=1.0).contains(&spec.continental_ratio),
            "{} continental_ratio out of [0, 1]",
            name
        );
        assert!(spec.bi >= 0.0 && spec.bi <= 5.0, "{} bi out of bounds", name);
        assert!(spec.br >= 0.0 && spec.br <= 5.0, "{} br out of bounds", name);

        // build_config must produce a valid BaselineConfig — this
        // exercises the Voronoi / Boundary / RecyclingConfig::validate
        // path without running the solver.
        let _cfg = build_config::build(&spec);
        // If we got here, the BaselineConfig literal was constructed —
        // the harness side of the integration is contract-clean.
    }
}
