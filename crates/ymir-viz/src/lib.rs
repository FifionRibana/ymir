//! Library facade for `ymir-viz`.
//!
//! The legacy v2 facade (the `tectonics_v2` Stokes bridge + its UI panels,
//! exposed here for regression tests) was removed with the v2 engine. The viz is
//! now a single-engine (C1) Bevy binary; all logic lives in the bin entrypoint
//! (`main.rs`) and its module tree. This library target is intentionally empty.
