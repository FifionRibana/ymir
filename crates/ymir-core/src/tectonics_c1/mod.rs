//! # tectonics_c1 — Lightweight dynamic tectonics
//!
//! Successor to `tectonics_v2/` per Step 12 R7 final report
//! (`docs/reports/step12_r7_final_report.md`) and the design document
//! at `docs/design/c1_lightweight_dynamic_tectonics.md`.
//!
//! ## Design in one sentence
//!
//! Advection of `S̃` and `age` under prescribed per-plate velocity,
//! plus closed-form empirical source terms (Davis-Suppe orogeny,
//! Parsons-Sclater bathymetry, stream-power erosion, Airy isostasy)
//! applied per forward-Euler step. No Stokes solver.
//!
//! ## Phase 1.1 status (this module's first issue, #120)
//!
//! Advection only. No closures. The output of a Phase 1.1 run is
//! a transport-correctness check, **not** a plausible continent.
//! Closures land in Phase 1.2 (Davis-Suppe), Phase 1.3 (equilibrium
//! height), Phase 1.4 (erosion + isostasy + downstream).
//!
//! ## Module map
//!
//! - [`state`] — `C1State` (S̃, age, plate_id, plate_type, cratonic
//!   mask). Mirrors v2 field types where possible; introduces a
//!   small `BoolField` for the binary cratonic mask since
//!   `Field2D` is `f64`-only.
//! - [`kinematics`] — per-plate constant translation velocity
//!   (Phase 1.1). Smoothing across plate boundaries and R7-
//!   generalised sampling are Phase 2 work (§6.3 design doc).
//! - [`time_loop`] — forward-Euler advection of `S̃` and `age`,
//!   CFL-trivial step. Closures hook in here in Phase 1.2+.
//! - [`init`] — wraps v2 `generate_voronoi` + `init_s_field` to
//!   produce a `C1State`. Verbatim v2 reuse per §4.8.
//!
//! ## Reused v2 assets (per §4.8 of the design doc)
//!
//! - [`crate::tectonics_v2::advection::step_upwind`] — upwind scheme
//! - [`crate::tectonics_v2::field::Field2D`] and
//!   [`crate::tectonics_v2::field::PeriodicIndex`]
//! - [`crate::tectonics_v2::voronoi::generate_voronoi`] and
//!   [`crate::tectonics_v2::voronoi::VoronoiConfig`]
//! - [`crate::tectonics_v2::init::init_s_field`] and
//!   [`crate::tectonics_v2::init::InitMode`]
//! - [`crate::tectonics_v2::boundaries::plate_type::PlateType`] and
//!   [`crate::tectonics_v2::boundaries::plate_type::PlateTypeField`]
//! - `crate::tectonics::isostasy::compute_isostasy` (legacy
//!   module, paradigm-agnostic, surfaces the heightmap for viz)
//!
//! C1 does **not** reuse v2's `stokes`, `mantle`, `slab`, `rheology`,
//! `basal_drag`, or `presets` — those moved to `tectonics_v2/_attic/`
//! under the `v2_legacy` feature in Issue #117.
//!
//! ## Phase 1.1 outputs
//!
//! See [`docs/reports/c1_phase_1_1_advection/README.md`](../../../docs/reports/c1_phase_1_1_advection/README.md)
//! for the visual + scalar acceptance evidence: 5 cycle snapshots
//! (cycle 0 / 50 / 100 / 200 / 300), per-cycle `S̃` distribution
//! stats, runtime, and the "what this output is not" guide.
//!
//! The integration test that produces those outputs lives in
//! `crates/ymir-core/tests/c1_phase_1_1_advection.rs`. Mass-
//! conservation drift is asserted < `1e-6`; the measured drift
//! is `1.6e-14` (machine precision, well below the threshold).

pub mod closures;
pub mod init;
pub mod kinematics;
pub mod state;
pub mod time_loop;
