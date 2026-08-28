//! Per-time-step application of the equilibrium height closure.
//!
//! See the parent module ([`super`]) for the physics derivation
//! (quadratic per Molnar & Lyon-Caen 1988 eq. 2), the global-cell
//! application rationale, and the interaction with Phase 1.2's
//! Davis-Suppe imprint.

use crate::tectonics_v2::field::Field2D;

use super::params::EquilibriumHeightParams;

/// Apply one forward-Euler step of the equilibrium-height closure
/// to the `S̃` field, **in place**.
///
/// For each cell `c`:
///
/// ```text
///     S̃_new(c) = S̃(c) − k_collapse · max(0, S̃(c) − h_eq)² · dt
/// ```
///
/// with a defensive clamp at `h_eq` (see the in-function comment
/// for why this triggers as part of the *intended* threshold
/// behavior for large-excess cells, and is non-physical safety
/// for the small-excess regime).
///
/// Unlike
/// [`super::super::davis_suppe::source_term::apply_davis_suppe_step`],
/// this closure does **not** skip any cell category: boundary
/// cells, cratonic cells, oceanic cells, wedge-body cells — all
/// equally subject to the same gravitational stability criterion.
/// That global application is what makes equilibrium height the
/// first effective sink in C1 (Phase 1.2 Davis-Suppe is source-
/// only, applied selectively on upper-plate wedge bodies).
///
/// Inputs:
/// - `s`: mutable `S̃` field, updated in place
/// - `params`: closure tunables (see
///   [`EquilibriumHeightParams`])
/// - `dt`: time step in the same units as `1 / (k_collapse ·
///   length)`
pub fn apply_equilibrium_height_step(s: &mut Field2D, params: &EquilibriumHeightParams, dt: f64) {
    if !params.enabled {
        return;
    }
    let nx = s.nx();
    let ny = s.ny();
    let h_eq = params.h_eq;
    let k_collapse = params.k_collapse;
    for j in 0..ny {
        for i in 0..nx {
            let s_now = s.get(i, j);
            let excess = s_now - h_eq;
            if excess <= 0.0 {
                continue;
            }
            let mut s_new = s_now - k_collapse * excess.powi(2) * dt;
            // Safety clamp: prevent undershoot of h_eq.
            // For the quadratic formula this is the intended
            // threshold behavior on large-excess cells (k · excess²
            // · dt > excess when k · excess · dt > 1), where one
            // step caps the cell at h_eq.
            // For small-excess wedge cells (k · excess · dt « 1)
            // the clamp never triggers and the formula evolves
            // smoothly.
            if s_new < h_eq {
                s_new = h_eq;
            }
            s.set(i, j, s_new);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Fill a fresh `nx × ny` field with a uniform `value`.
    fn fill_with(nx: usize, ny: usize, value: f64) -> Field2D {
        let mut s = Field2D::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                s.set(i, j, value);
            }
        }
        s
    }

    #[test]
    fn no_op_below_h_eq() {
        // S̃ uniformly under h_eq → asymmetric closure must not act.
        let params = EquilibriumHeightParams::default();
        let initial = 1.0;
        assert!(initial < params.h_eq, "test premise: initial S̃ must be below h_eq");
        let mut s = fill_with(4, 4, initial);
        let before = s.data().to_vec();
        apply_equilibrium_height_step(&mut s, &params, 1.0);
        let after = s.data();
        for k in 0..after.len() {
            assert_eq!(
                before[k], after[k],
                "S̃ below h_eq must not change; mismatch at flat index {k}"
            );
        }
    }

    #[test]
    fn collapse_above_h_eq() {
        // Quadratic formula sanity check with the post-E1.bis
        // defaults (k_collapse = 2.0):
        //   excess     = 1.0
        //   decrement  = 2.0 · 1.0² · 0.1 = 0.2
        //   S̃_new     = 3.0 - 0.2 = 2.8
        // k · excess · dt = 0.2 « 1 → clamp inactive, exact
        // quadratic formula.
        let params = EquilibriumHeightParams::default();
        let initial = 3.0;
        let dt = 0.1;
        let mut s = fill_with(4, 4, initial);
        apply_equilibrium_height_step(&mut s, &params, dt);
        let expected = initial - params.k_collapse * (initial - params.h_eq).powi(2) * dt;
        for j in 0..s.ny() {
            for i in 0..s.nx() {
                let got = s.get(i, j);
                assert!(
                    (got - expected).abs() < 1e-12,
                    "expected S̃ = {expected}, got {got} at ({i},{j})"
                );
                assert!(got < initial, "collapse must reduce S̃");
                assert!(
                    got > params.h_eq,
                    "single subcritical step from S̃=3 must not reach h_eq; got {got}"
                );
            }
        }
    }

    #[test]
    fn smooth_at_h_eq_boundary() {
        // S̃ exactly at h_eq: excess = 0 → no change. Locks the
        // closure against numerical "kicks" through the threshold
        // when the system sits at equilibrium.
        let params = EquilibriumHeightParams::default();
        let mut s = fill_with(4, 4, params.h_eq);
        let before = s.data().to_vec();
        apply_equilibrium_height_step(&mut s, &params, 1.0);
        let after = s.data();
        for k in 0..after.len() {
            assert_eq!(
                before[k], after[k],
                "S̃ at h_eq must not change; mismatch at flat index {k}"
            );
        }
    }

    #[test]
    fn disabled_no_op() {
        // enabled = false → bit-identical pre/post regardless of
        // input values (including cells well above h_eq). Mirrors
        // the W4 closure-isolation discipline used by Davis-Suppe.
        let params =
            EquilibriumHeightParams { enabled: false, ..EquilibriumHeightParams::default() };
        let mut s = fill_with(4, 4, 5.0); // uniformly above h_eq
        // Sprinkle non-uniform values so a forgotten branch can't
        // pass by accident on a uniform field.
        s.set(0, 0, 100.0);
        s.set(1, 2, -3.0);
        s.set(3, 3, params.h_eq);
        let before = s.data().to_vec();
        apply_equilibrium_height_step(&mut s, &params, 1.0);
        let after = s.data();
        for k in 0..after.len() {
            assert_eq!(
                before[k], after[k],
                "disabled closure must not touch any cell; mismatch at flat index {k}"
            );
        }
    }

    #[test]
    fn never_undershoots_h_eq() {
        // Pathological k_collapse · dt = 100 — without the safety
        // clamp the quadratic formula would predict
        //   S̃_new = 3 − 100 · (3 − 2)² · 1 = 3 − 100 = −97.
        // The clamp must hold S̃_new at h_eq = 2.
        //
        // This test locks the *defensive / threshold* clamp; in
        // normal use on small-excess wedge cells `k · excess · dt
        // « 1` and the clamp never triggers. On large-excess
        // boundary outliers the clamp triggering is the intended
        // one-step cap (see module docstring).
        let params =
            EquilibriumHeightParams { k_collapse: 100.0, ..EquilibriumHeightParams::default() };
        let initial = 3.0;
        let dt = 1.0;
        let predicted_unclamped =
            initial - params.k_collapse * (initial - params.h_eq).powi(2) * dt;
        assert!(
            predicted_unclamped < params.h_eq,
            "test premise: unclamped prediction must undershoot h_eq (got {predicted_unclamped})"
        );
        let mut s = fill_with(4, 4, initial);
        apply_equilibrium_height_step(&mut s, &params, dt);
        for j in 0..s.ny() {
            for i in 0..s.nx() {
                let got = s.get(i, j);
                assert_eq!(
                    got, params.h_eq,
                    "safety clamp must hold S̃ at h_eq; got {got}, expected {} at ({i},{j})",
                    params.h_eq
                );
            }
        }
    }
}
