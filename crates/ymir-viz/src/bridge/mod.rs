//! Bridge module — c1 (lightweight dynamic tectonics).
//!
//! The legacy v2 (Stokes) bridge was removed with the v2 engine sunset; the viz
//! now wires a single engine, `bridge::c1::*`, on a background worker thread.

pub mod c1;
