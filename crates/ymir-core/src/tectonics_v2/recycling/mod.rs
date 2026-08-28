//! Conservative mass recycling (Step 6).
//!
//! Subducted mass is recycled through distributive fractions: part
//! goes to immediate arc/collision/rift volcanism, part to delayed
//! mid-ocean spreading (via a mantle-residence ring buffer), and
//! optionally a small fraction is permanently lost to the deep
//! mantle. See [`RecyclingConfig`] for the fractions and
//! [`DelayedRecycler`] for the buffer wrapper.
//!
//! The legacy [`RecyclingBuffer`] is re-exported; it provides the
//! ring-buffer primitive (deposit/advance). Step 6 adds the
//! rollover semantics (mass that emerged from the buffer but cannot
//! be distributed is held at the output without aging further) and
//! the pipeline bookkeeping (arc/coll_v/rift_v immediate
//! accumulators with identical rollover rules).

pub mod buffer;
pub mod config;

pub use buffer::{DelayedRecycler, RecyclingBuffer};
pub use config::{ImmediateAccumulators, RecyclingConfig, RecyclingConfigError};
