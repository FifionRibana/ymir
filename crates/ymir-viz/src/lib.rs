//! Library facade for `ymir-viz`. Step 8.6 introduced this so the v2
//! bridge can be exercised by integration tests in `tests/` without
//! pulling in the full Bevy plugin graph.
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

    #[path = "../visualization/v2_viz.rs"]
    pub mod v2_viz;
}
