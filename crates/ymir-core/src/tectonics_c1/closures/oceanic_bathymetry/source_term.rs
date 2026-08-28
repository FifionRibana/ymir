//! Per-time-step application of the Stein-Stein 1992 oceanic
//! bathymetry closure.
//!
//! See the parent module ([`super`]) for the physics derivation,
//! the Architecture C rationale (post-isostasy altitude assignment
//! rather than additive `S̃` source term), and the interaction with
//! Phase 1.4 stream-power erosion through the joint stage-4a
//! altitude preparation.

use crate::grid::GridF32;
use crate::tectonics_v2::boundaries::plate_type::{PlateType, PlateTypeField};
use crate::tectonics_v2::field::Field2D;

use super::params::SteinSteinParams;

/// Continuity offset `C` between the S-S young and old regimes.
///
/// Per Stein & Stein 1992 eq. (4)-(5), the old-regime formula is
/// literally
///
/// ```text
///     d(t) = 5651 - 2473 · exp(-0.0278 · t)         [t ≥ 20 Ma]
/// ```
///
/// where the coefficient `2473 ≈ asymptotic_depth - ridge_depth + 22 m`
/// (the extra `22 m` is a small continuity adjustment so the two
/// regimes match at `t = t_c = 20 Ma` to within ~0.4 m). This
/// constant is hard-coded here rather than exposed as a
/// [`SteinSteinParams`] field because changing it independently of
/// `asymptotic_depth_m` and `ridge_depth_m` would break the
/// formula's published fidelity — the three parameters together
/// define the paper's calibration and must move together if at all.
const SS_OLD_AMPLITUDE: f64 = 2473.0;

/// Compute the Stein-Stein 1992 ocean-floor depth in meters at a
/// given age in Ma.
///
/// Per Stein, C. A. & Stein, S. (1992). "A model for the global
/// variation in oceanic depth and heat flow with lithospheric
/// age." *Nature* 359, 123-129, eq. (4)-(5):
///
/// - Young regime (`age_ma < crossover_age_ma`, default `20 Ma`):
///
///   ```text
///       d(t) = ridge_depth_m + subsidence_rate · √t
///   ```
///
///   Square-root cooling per the half-space conductive model;
///   subsidence rate `b ≈ 365 m / √Ma` in the canonical GDH1
///   parameter set.
/// - Old regime (`age_ma ≥ crossover_age_ma`):
///
///   ```text
///       d(t) = asymptotic_depth_m - SS_OLD_AMPLITUDE · exp(-α · t)
///   ```
///
///   Exponential saturation to the asymptote `d_∞ ≈ 5651 m` per the
///   plate model; thermal time constant `α ≈ 0.0278 Ma⁻¹`. The
///   `2473 m` continuity offset comes from the private
///   `SS_OLD_AMPLITUDE` constant (see module-private const doc for
///   the rationale on hard-coding rather than exposing as a param).
///
/// Negative ages are clamped to `0` (ridge-axis depth). This guards
/// the `√t` in the young regime against `NaN` from a malformed age
/// field; in normal C1 use the age field is non-negative by
/// construction.
///
/// Returns the depth in meters (positive, downward from sea level).
pub fn stein_stein_depth(age_ma: f64, params: &SteinSteinParams) -> f64 {
    let age = age_ma.max(0.0);
    if age < params.crossover_age_ma {
        params.ridge_depth_m + params.subsidence_rate * age.sqrt()
    } else {
        params.asymptotic_depth_m - SS_OLD_AMPLITUDE * (-params.time_constant * age).exp()
    }
}

/// Apply the Stein-Stein 1992 oceanic bathymetry to the altitude
/// field, **in place**.
///
/// For each cell `c` classified as `PlateType::Oceanic`:
///
/// ```text
///     altitude_new(c) = − stein_stein_depth(age_ma(c)) / depth_scale_m
/// ```
///
/// with `age_ma(c) = age(c) · params.age_to_ma`. Continental cells
/// are left unchanged. Plate-type is the canonical oceanic/
/// continental classifier (Stein-Stein semantics — "oceanic
/// lithosphere subsides with age"); no altitude threshold (e.g.,
/// `< sea_level`) is consulted, since lifting a continental cell
/// below sea level does not by itself make it oceanic in the S-S
/// regime, and an oceanic cell uplifted by the upstream
/// tessellation (transient or artefact) should still receive S-S
/// bathymetry to keep the closure's domain coherent.
///
/// **Sign convention.** Oceanic cells receive negative altitude
/// values in `[−1.13, −0.52]` (corresponding to the S-S depth
/// range `[5651, 2600] m` divided by `depth_scale_m = 5000`).
/// "Sea level = 0" in this convention. This breaks the
/// `compute_isostasy` `[0, 1]` normalisation contract on oceanic
/// cells, but the downstream consumers in the C1 stage-4 pipeline
/// (drainage routing operates in `S̃` space; erosion operates on
/// slope magnitudes) are sign-insensitive — they care about
/// gradients, not absolute level.
///
/// **Architecture C**: this mutates the `altitude` heightmap only;
/// `S̃` is not touched. See [`super`] for the rationale and the
/// fallback architectures (A — `S̃` source term; B — hybrid) that
/// would be considered if Stage D visual review surfaces
/// limitations of this design.
///
/// Inputs:
/// - `altitude`: mutable post-isostasy heightmap, updated in place.
///   `f32`-backed [`GridF32`]; the S-S conversion `depth_m /
///   depth_scale_m` is computed in `f64` then cast.
/// - `age`: per-cell non-dim age. `Field2D` (`f64`-backed).
///   Multiplied by `params.age_to_ma` at use to obtain age in Ma.
/// - `plate_type`: per-cell oceanic/continental classification.
///   Only oceanic cells are modified.
/// - `params`: closure tunables (see [`SteinSteinParams`]).
///
/// `params.enabled = false` makes the call a no-op (W4
/// closure-isolation discipline). The early-return is the only
/// guard; no per-cell work runs on the disabled path.
pub fn apply_stein_stein_bathymetry(
    altitude: &mut GridF32,
    age: &Field2D,
    plate_type: &PlateTypeField,
    params: &SteinSteinParams,
) {
    if !params.enabled {
        return;
    }
    let nx = altitude.width;
    let ny = altitude.height;
    debug_assert_eq!(age.nx(), nx, "Stein-Stein: age field nx must match altitude grid");
    debug_assert_eq!(age.ny(), ny, "Stein-Stein: age field ny must match altitude grid");
    debug_assert_eq!(
        plate_type.nx(),
        nx,
        "Stein-Stein: plate_type field nx must match altitude grid"
    );
    debug_assert_eq!(
        plate_type.ny(),
        ny,
        "Stein-Stein: plate_type field ny must match altitude grid"
    );

    let depth_scale = params.depth_scale_m;
    let age_to_ma = params.age_to_ma;

    for j in 0..ny {
        for i in 0..nx {
            if plate_type.get(i, j) != PlateType::Oceanic {
                continue;
            }
            let age_ma = age.get(i, j) * age_to_ma;
            let depth_m = stein_stein_depth(age_ma, params);
            let depth_nondim = (depth_m / depth_scale) as f32;
            altitude.set(i, j, -depth_nondim);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a small `(altitude, age, plate_type)` fixture with the
    /// top half continental, bottom half oceanic. Age uniform at
    /// `30.0` (`age_ma = 30 · 0.667 ≈ 20`, on the boundary between
    /// young and old regimes — close enough to ridge-shallow for
    /// the test to see a depth change).
    fn split_continental_oceanic(nx: usize, ny: usize) -> (GridF32, Field2D, PlateTypeField) {
        let altitude = GridF32::new(nx, ny, 0.5);
        let mut age = Field2D::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                age.set(i, j, 30.0);
            }
        }
        let mut plate_type = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
        for j in 0..ny / 2 {
            for i in 0..nx {
                plate_type.set(i, j, PlateType::Continental);
            }
        }
        (altitude, age, plate_type)
    }

    /// `stein_stein_depth(0, _) = ridge_depth_m`. The young regime
    /// at `t = 0` collapses to the ridge axis: `d_r + b · √0 = d_r`.
    #[test]
    fn stein_stein_depth_at_ridge_axis() {
        let params = SteinSteinParams::default();
        let depth = stein_stein_depth(0.0, &params);
        assert!(
            (depth - params.ridge_depth_m).abs() < 1e-9,
            "ridge-axis depth: got {depth}, expected {} (ridge_depth_m)",
            params.ridge_depth_m
        );
    }

    /// W-T-style scaling check on the young regime: `4× age` →
    /// `2× subsidence`. Verifies the `√t` law is correctly
    /// implemented (not `t`, not `t²`).
    #[test]
    fn stein_stein_depth_young_regime_sqrt_scaling() {
        let params = SteinSteinParams::default();
        // Two ages well inside the young regime
        // (`< crossover_age_ma = 20 Ma`).
        let age1 = 4.0;
        let age2 = 16.0;
        let d1 = stein_stein_depth(age1, &params);
        let d2 = stein_stein_depth(age2, &params);
        let subsidence1 = d1 - params.ridge_depth_m;
        let subsidence2 = d2 - params.ridge_depth_m;
        let ratio = subsidence2 / subsidence1;
        // Expected: √(age2/age1) = √4 = 2.
        let expected = (age2 / age1).sqrt();
        assert!(
            (ratio - expected).abs() < 1e-9,
            "subsidence ratio for 4× age: got {ratio:.6}, expected {expected:.6} (√t law)"
        );
    }

    /// Old regime exponential saturation: as `t → ∞`,
    /// `exp(-α·t) → 0`, so `d → asymptotic_depth_m`.
    #[test]
    fn stein_stein_depth_old_regime_exponential() {
        let params = SteinSteinParams::default();
        // 1000 Ma is well past saturation: `exp(-0.0278 · 1000) ≈
        // 8 × 10⁻¹³`, contribution `< 1e-9 m`.
        let depth_at_large_t = stein_stein_depth(1000.0, &params);
        let residual = (depth_at_large_t - params.asymptotic_depth_m).abs();
        assert!(
            residual < 1e-3,
            "old-regime asymptote: got {depth_at_large_t}, expected {} (asymptotic_depth_m); residual {residual}",
            params.asymptotic_depth_m
        );
    }

    /// Continuity between young and old regimes at the crossover
    /// age. S-S 1992's `2473 m` coefficient is calibrated so the
    /// two formulas agree at `t = 20 Ma` to within ~0.4 m. The
    /// regimes' branch condition is `<` vs `≥`, so we sample
    /// `crossover - ε` (young branch) and `crossover + ε` (old
    /// branch) and check the residual.
    #[test]
    fn stein_stein_depth_continuity_at_crossover() {
        let params = SteinSteinParams::default();
        let t_c = params.crossover_age_ma;
        let eps = 0.01_f64;
        let young_just_below = stein_stein_depth(t_c - eps, &params);
        let old_just_above = stein_stein_depth(t_c + eps, &params);
        let residual = (old_just_above - young_just_below).abs();
        assert!(
            residual < 5.0,
            "young/old continuity at t_c = {t_c} Ma: young({})={young_just_below:.3} vs old({})={old_just_above:.3}, residual {residual:.3} m exceeds 5.0 m tolerance",
            t_c - eps,
            t_c + eps,
        );
    }

    /// Continental cells are untouched by the closure; oceanic
    /// cells are modified. Verifies the `PlateType::Oceanic` filter
    /// in [`apply_stein_stein_bathymetry`].
    #[test]
    fn apply_bathymetry_skips_continental() {
        let (mut altitude, age, plate_type) = split_continental_oceanic(4, 4);
        let initial = altitude.data.clone();
        let params = SteinSteinParams::default();

        apply_stein_stein_bathymetry(&mut altitude, &age, &plate_type, &params);

        // Continental top half (j < 2) must be bit-identical.
        for j in 0..2 {
            for i in 0..4 {
                let idx = j * 4 + i;
                assert_eq!(
                    altitude.data[idx], initial[idx],
                    "continental cell ({i}, {j}) was modified; before={} after={}",
                    initial[idx], altitude.data[idx]
                );
            }
        }
        // Oceanic bottom half (j ≥ 2) must show a change on at
        // least one cell (age = 30 → age_ma ≈ 20 → some depth).
        let mut any_oceanic_changed = false;
        for j in 2..4 {
            for i in 0..4 {
                let idx = j * 4 + i;
                if (altitude.data[idx] - initial[idx]).abs() > 1e-6 {
                    any_oceanic_changed = true;
                }
            }
        }
        assert!(
            any_oceanic_changed,
            "no oceanic cell was modified; S-S apply did not run on the oceanic half"
        );
    }

    /// `enabled = false` → bit-identical altitude regardless of
    /// inputs. Mirrors the W4 closure-isolation discipline used by
    /// Davis-Suppe, equilibrium-height, and erosion.
    #[test]
    fn apply_bathymetry_skips_disabled() {
        let (mut altitude, age, plate_type) = split_continental_oceanic(4, 4);
        let initial = altitude.data.clone();
        let params = SteinSteinParams { enabled: false, ..SteinSteinParams::default() };

        apply_stein_stein_bathymetry(&mut altitude, &age, &plate_type, &params);

        for k in 0..initial.len() {
            assert_eq!(
                altitude.data[k], initial[k],
                "disabled closure modified cell at flat index {k}: before={} after={}",
                initial[k], altitude.data[k]
            );
        }
    }
}
