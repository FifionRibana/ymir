//! Library facade for `ymir-viz`. Step 8.6 introduced this so the v2
//! bridge can be exercised by integration tests in `tests/` without
//! pulling in the full Bevy plugin graph.
//!
//! Issue #117 — the whole library is gated behind Cargo feature
//! `v2_legacy` because every module (`bridge::v2::*`, the v2 UI
//! panels, the v2 visualization helpers) wires directly into the
//! retired `tectonics_v2/{stokes, mantle, slab, rheology, basal_drag,
//! presets}` subtree. The default build of `ymir-viz` compiles as an
//! empty library; `cargo {build,test} --features v2_legacy` rebuilds
//! the full surface for regression reproducibility. C1's bridge +
//! UI will be reintroduced under a paradigm-agnostic facade in C1
//! Phase 4 (§7 C1.md).

#![cfg(feature = "v2_legacy")]
//!
//! The library scope is narrow: the v2 bridge tree and the v2-only
//! visualization helpers are re-exported. Legacy bridge plugin
//! systems and other crate-local resources stay in the binary
//! entrypoint (`main.rs`).
//!
//! Path layout: the inline `mod bridge { pub mod v2; }` matches the
//! binary's `mod bridge;` tree exactly, so files under `bridge/v2/`
//! that resolve `crate::bridge::v2::*` compile identically in both
//! the bin and the lib. `pub use bridge::v2 as bridge_v2` is a
//! convenience alias for tests / external callers that prefer the
//! flatter name.

pub mod bridge {
    #[path = "../bridge/v2/mod.rs"]
    pub mod v2;
}

pub use bridge::v2 as bridge_v2;

pub mod visualization {
    #[path = "../visualization/colormap.rs"]
    pub mod colormap;

    #[path = "../visualization/overlay.rs"]
    pub mod overlay;

    #[path = "../visualization/v2_viz.rs"]
    pub mod v2_viz;
}

pub mod ui {
    // Mounted on the library facade so the `update_v2_preview`
    // system in `visualization::v2_viz` can reach
    // `parameter_panel_v2::V2EditableSpec` from a `crate::ui::…`
    // path. The bin keeps its own `mod ui;` tree (which also has
    // `metrics_dashboard` and `mod.rs`'s `UiPlugin`); the lib mounts
    // only what `v2_viz` references plus what integration tests
    // under `tests/` need to reach.
    #[path = "../ui/parameter_panel_v2.rs"]
    pub mod parameter_panel_v2;

    // Step 12 Phase 7b — workflow_panel is mounted on the lib facade
    // so its unit tests (`cargo test --lib ui::workflow_panel`) and
    // integration tests under `tests/` can reach
    // `WorkflowCycleHistory`, `CycleMetricsSnapshot`, etc.
    #[path = "../ui/workflow_panel.rs"]
    pub mod workflow_panel;
}

pub mod pipeline;

pub mod phases;
