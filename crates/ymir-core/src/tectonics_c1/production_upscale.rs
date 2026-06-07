//! Production C1 → anisotropic-FBM upscale wiring (Issue #147
//! FOLLOWUPS-#6 outcome).
//!
//! ## The robustness contract (enforced by the type, not just documented)
//!
//! The S̃ thickness field is mesh-NON-convergent in production
//! (downsample-correlation r~0.51 across 64²/128²/256²/512² — Issue
//! #147). The upscale is nonetheless ROBUST to this — 64²+upscale and
//! 256²+upscale produce the SAME world (structure r~0.90, measured in
//! `stage_upscale_robustness.md`) — **only because it reads the
//! laundered production ALTITUDE** (isostasy + Stein-Stein, which
//! converges at alt r~0.88), NOT the raw S̃ gradient (r~0.51). The
//! non-convergence is washed out by isostasy + Stein-Stein before the
//! FBM orientation ever sees it (same absorption that hides the
//! Davis-Suppe mass swing in production).
//!
//! Feeding raw S̃ to the upscale would silently break this: the FBM
//! would orient on a non-convergent gradient, 64²+FBM and 256²+FBM
//! would DIVERGE (invisible at a single resolution; only a 2-resolution
//! comparison reveals it), and it would reopen FOLLOWUPS #6 (the
//! deferred advection-scheme milestone).
//!
//! [`upscale_from_c1`] makes that illegal state **unrepresentable**: it
//! takes a [`C1State`] and builds the altitude INTERNALLY. There is no
//! S̃/heightmap parameter a caller could substitute. The contract is
//! enforced by the signature; the inline comment is belt-and-braces.
//!
//! The regression test
//! `c1_closure_morphology::upscale_from_c1_structure_converges`
//! asserts the precondition (64²/256² structure convergence) so a
//! future change that breaks it is caught, not just warned against.

use crate::grid::GridF32;
use crate::seed::WorldSeed;
use crate::tectonics::isostasy::{compute_isostasy, IsostasyConfig};
use crate::tectonics_v2::boundaries::plate_type::PlateTypeField;
use crate::tectonics_v2::field::Field2D;
use crate::terrain::upscale::{upscale_with_fbm, FbmUpscaleConfig, UpscaleResult};

use super::closures::oceanic_bathymetry::params::SteinSteinParams;
use super::closures::oceanic_bathymetry::source_term::apply_stein_stein_bathymetry;
use super::state::C1State;

/// The C1 PRODUCTION altitude (Architecture C): `compute_isostasy`
/// followed by the Stein-Stein bathymetry re-apply on oceanic cells.
///
/// **SINGLE SOURCE OF TRUTH** for the production altitude, shared by
/// the viz render (`bridge::c1` / `c1_viz::derive_altitude_field`, which
/// reconstructs these fields from a snapshot and calls this) AND the
/// upscale input ([`upscale_from_c1`]). Factoring it here (Issue #147
/// #6) prevents the upscaled HD terrain from silently diverging from
/// the rendered production altitude: the two agree BY CONSTRUCTION, not
/// by two copies of the same `isostasy + Stein-Stein` sequence drifting
/// apart. The regression test guards robustness; this shared function
/// guards render/upscale equality.
pub fn c1_production_altitude(
    s: &Field2D,
    age: &Field2D,
    plate_type: &PlateTypeField,
    iso: &IsostasyConfig,
    ss: &SteinSteinParams,
) -> GridF32 {
    let isostasy = compute_isostasy(s, iso);
    let mut altitude = isostasy.heightmap;
    apply_stein_stein_bathymetry(&mut altitude, age, plate_type, ss);
    altitude
}

/// Fixed altitude→`[0,1]` normalisation half-range (sea level 0.0 maps
/// to 0.5). This MUST be a constant, NOT a per-call data range:
/// normalising 64² and 256² by their own min/max would make them
/// incomparable and break the measured cross-resolution robustness.
/// Value matches the production gallery palette half-range.
const ALTITUDE_NORM_HALF_RANGE: f32 = 1.13;

/// Upscale a C1 simulation state to a detailed heightmap via the
/// anisotropic-FBM upscale, **through the robustness-preserving
/// altitude path** (Issue #147 FOLLOWUPS-#6 — see module docs).
///
/// Takes the [`C1State`] (not an S̃ field, not a heightmap) BY DESIGN:
/// the altitude is built internally so no caller can feed the
/// non-convergent raw S̃ to the FBM orientation.
///
/// Pipeline (identical to the measured-robust `stage_upscale_robustness`
/// path): `compute_isostasy(iso)` → `apply_stein_stein_bathymetry(ss)`
/// → fixed `[0,1]` normalisation (sea→0.5) → `upscale_with_fbm`.
pub fn upscale_from_c1(
    state: &C1State,
    iso: &IsostasyConfig,
    ss: &SteinSteinParams,
    seed: &WorldSeed,
    cfg: &FbmUpscaleConfig,
) -> UpscaleResult {
    // CONTRACT (Issue #147 #6): the upscale reads the laundered
    // ALTITUDE (isostasy + Stein-Stein, convergent r~0.88), NOT raw S̃
    // (non-convergent r~0.51). This is what makes the upscale robust to
    // S̃ mesh non-convergence. Built via the SHARED `c1_production_altitude`
    // (same source of truth as the viz render) — do NOT inline a raw-S̃
    // path here; it reopens FOLLOWUPS #6 and silently diverges.
    let mut coarse = c1_production_altitude(&state.s, &state.age, &state.plate_type, iso, ss);

    // Fixed normalisation to [0,1] (sea 0.0 → 0.5), resolution-independent.
    let half = ALTITUDE_NORM_HALF_RANGE;
    for v in coarse.data.iter_mut() {
        *v = ((*v + half) / (2.0 * half)).clamp(0.0, 1.0);
    }
    let sea_level_normalized = 0.5_f32;

    upscale_with_fbm(&coarse, sea_level_normalized, seed, cfg)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics_c1::init_r7::{init_c1_state_phase_2_r7, Phase2InitParams};

    /// Smoke: the contract path runs and produces a target-sized
    /// heightmap with finite values (no NaN/Inf), at a small target.
    #[test]
    fn upscale_from_c1_smoke() {
        let state = init_c1_state_phase_2_r7(64, 42, &Phase2InitParams::default());
        let cfg = FbmUpscaleConfig { target_size: 128, ..Default::default() };
        let out = upscale_from_c1(
            &state,
            &IsostasyConfig::c1_default(),
            &SteinSteinParams::default(),
            &WorldSeed::new(42),
            &cfg,
        );
        assert_eq!(out.heightmap.width, 128);
        assert_eq!(out.heightmap.height, 128);
        assert!(out.heightmap.data.iter().all(|v| v.is_finite()));
    }
}
