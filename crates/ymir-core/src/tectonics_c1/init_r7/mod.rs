//! R7 generalised init for C1 Phase 2 Track B (Issue #131).
//!
//! Sibling of [`super::init`] (Phase 1.1 init, preserved verbatim).
//! Phase 2 Track B introduces three sub-components that compose on
//! top of the v2 Voronoï output:
//!
//! 1. **Boundary displacement** ([`boundary_displacement`]) — apply
//!    Perlin / Simplex noise displacement to the per-cell sampling
//!    position before re-querying the nearest Voronoï seed. Produces
//!    non-rectilinear plate boundaries while preserving the seed-
//!    based plate identity. Resolves the v1 / v2 visual failure mode
//!    of orogenic chains aligned along straight Voronoï edges
//!    (§6.1 design doc). Documented as the "boundary displacement"
//!    option of §6.1 — Lloyd relaxation and multi-scale overlay
//!    remain deferred.
//! 2. **Continental clustering** (Stage E2, separate file) — BFS
//!    cluster-based plate-type assignment producing a cadrable
//!    continental cluster.
//! 3. **Ridge-aligned age = 0** (Stage E3, separate file) — set
//!    `age = 0` on cells adjacent to divergent boundaries at init
//!    time. Resolves the Phase 2 Track A finding that flux-form
//!    advection of `age` produces ~1000× density pile-up at
//!    convergent boundaries (`feedback_age_advection_density_vs_lagrangian`).
//!
//! ## Determinism
//!
//! All R7 init sub-components are deterministic given
//! `(grid_size, params.seed)`. Stochastic elements:
//!
//! - **Boundary displacement**: Perlin / Simplex noise via two
//!   independent `Fbm<Perlin>` instances (one per displacement
//!   component) seeded from `params.seed`. Same seed → same noise
//!   field → same per-cell displacement.
//! - **Continental clustering**: ChaCha8Rng seeded from
//!   `cluster_params.seed` for seed-pick selection. BFS expansion
//!   itself is deterministic given the adjacency graph.
//!
//! No floating-point reproducibility caveats: `f64` arithmetic and
//! `noise::Fbm<Perlin>` are bit-deterministic on a given target
//! triple per the existing Phase 1.1 + Track A precedents.
//!
//! ## Composition with Phase 1.1 init
//!
//! [`super::init::init_c1_state_phase_1_1`] is **preserved
//! verbatim** as the Phase 1.1 regression baseline. Phase 2 Track B
//! ships a parallel entry point (Stage E4) that calls the same v2
//! Voronoï + S̃-init pipeline, then chains the three R7
//! sub-components in order:
//!
//! ```text
//!     generate_voronoi → R7 boundary displacement
//!                      → cluster-BFS plate_type override
//!                      → init_s_field (per overridden plate_type)
//!                      → ridge-aligned age init
//!                      → cratonic mask
//! ```
//!
//! Phase 1.x tests continue to call
//! `init_c1_state_phase_1_1(grid_size, seed)` directly; the Phase
//! 2 entry is a new function call.

pub mod age_init;
pub mod boundary_displacement;
pub mod clustering;
pub mod params;

pub use age_init::{init_age_field_ridge_aligned, AgeInitParams};
pub use boundary_displacement::apply_boundary_displacement;
pub use clustering::{
    assign_continental_clusters, build_plate_adjacency, ContinentalClusterParams,
};
pub use params::R7InitParams;
