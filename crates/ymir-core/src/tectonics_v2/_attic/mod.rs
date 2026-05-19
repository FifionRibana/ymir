//! Issue #117 — v2 Stokes-coupled subtree retired from default build.
//!
//! Modules under `_attic/` were moved here from `tectonics_v2/` during the
//! C1 pivot prep (see `docs/design/c1_lightweight_dynamic_tectonics.md`
//! §4.7 for the rationale and `docs/migrations/v2_to_c1_attic.md` for the
//! audit that drove the move). They remain compilable under the Cargo
//! feature `v2_legacy` so regression baselines (R5b mf sweeps, R6.3 evo
//! sweep, R7.A.2.4 init mode sweep) stay reproducible bit-identically.
//! Default builds exclude this subtree to keep the C1 compile path lean.
//!
//! Retired modules:
//!
//! - [`stokes`] — non-linear Stokes solver subtree (Newton outer, CG/AMG
//!   inner, Picard fallback, continuation, snapshot, parallel reduction)
//! - [`mantle`] — mantle pattern + stream-function builder + evolution
//!   rate phase drift
//! - [`slab`] — slab-pull state ODE + accumulation + convergence
//!   direction (HC1 surfacing during audit — see migration doc)
//! - [`rheology`] — non-linear power-law viscosity + yielding law
//! - [`basal_drag`] — v2-form basal drag (operator-diagonal); any C1
//!   basal drag will be a closure, not an operator term
//! - [`presets`] — `RheologyParams` + `YieldingConfig` + `Preset`
//!   (HC2 surfacing during audit)
//! - [`forcing`] — the two retired body-force implementations
//!   (`SlabPullForce`, `MantleForce`); the preserved `GpeForce`,
//!   `SinusoidalForce`, `ForceSum`, `BodyForce` trait, etc. stay in
//!   `tectonics_v2/forcing/`

#![cfg(feature = "v2_legacy")]

pub mod basal_drag;
pub mod forcing;
pub mod mantle;
pub mod presets;
pub mod rheology;
pub mod slab;
pub mod stokes;
