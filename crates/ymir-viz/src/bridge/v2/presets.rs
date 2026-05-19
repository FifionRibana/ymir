//! Step 8.6 Phase 4 — preset registry.
//!
//! Presets are JSON files under `crates/ymir-viz/presets/v2/`. The
//! registry compiles them in via `include_str!` so the binary ships
//! self-contained (no external presets dir to lose at runtime).
//! `load(name)` deserializes the preset's payload into a `V2RunSpec`;
//! `list()` enumerates the available preset names for the UI dropdown.
//!
//! Adding a preset:
//! 1. Drop `name.json` under `crates/ymir-viz/presets/v2/`.
//! 2. Append `("name", include_str!("../../../presets/v2/name.json"))`
//!    to `EMBEDDED_PRESETS` below.
//! 3. The Phase 4 acceptance test (`v2_bridge_preset_load`) will
//!    automatically cover the new entry.

use super::spec::V2RunSpec;

/// Name → JSON payload. Embedded at compile time so presets travel
/// with the binary and tests can run without filesystem I/O.
const EMBEDDED_PRESETS: &[(&str, &str)] = &[
    (
        "quiescent",
        include_str!("../../../presets/v2/quiescent.json"),
    ),
    (
        "single_continent",
        include_str!("../../../presets/v2/single_continent.json"),
    ),
    (
        "convergence",
        include_str!("../../../presets/v2/convergence.json"),
    ),
    (
        "subduction",
        include_str!("../../../presets/v2/subduction.json"),
    ),
    (
        "divergence",
        include_str!("../../../presets/v2/divergence.json"),
    ),
    (
        "active_medley",
        include_str!("../../../presets/v2/active_medley.json"),
    ),
];

/// List the available preset names in registration order.
pub fn list() -> Vec<&'static str> {
    EMBEDDED_PRESETS.iter().map(|(name, _)| *name).collect()
}

/// Deserialize the named preset into a `V2RunSpec`. Returns
/// `Err(String)` with the JSON-parse error path on failure (matches
/// the surface area expected by `tests/v2_bridge_preset_load.rs`).
pub fn load(name: &str) -> Result<V2RunSpec, String> {
    let (_, payload) = EMBEDDED_PRESETS
        .iter()
        .find(|(n, _)| *n == name)
        .ok_or_else(|| format!("unknown preset '{}'", name))?;
    serde_json::from_str::<V2RunSpec>(payload)
        .map_err(|e| format!("preset '{}' deserialize error: {}", name, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Lib-internal smoke: every embedded preset deserializes. Mirrors
    /// the integration test under `tests/` but runs as part of
    /// `cargo test -p ymir-viz --lib`.
    #[test]
    fn every_embedded_preset_loads() {
        for name in list() {
            let spec = load(name).expect(name);
            assert!(spec.grid_nx > 0);
            assert!(spec.grid_ny > 0);
            assert!(spec.steps > 0);
            assert!(spec.num_plates > 0);
        }
    }
}
