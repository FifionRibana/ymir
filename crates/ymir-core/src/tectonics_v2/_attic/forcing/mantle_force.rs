//! Mantle forcing body-force term (Step 8).
//!
//! Continuous form (D1, §4.9):
//!
//! ```text
//!   f̃_mantle(x) = coupling · S̃(x) · (Mf · v_pattern(x) − v_solved(x))
//! ```
//!
//! Refactor (Step 4 basal-drag pattern) for exact self-consistency
//! at every Newton outer iteration:
//!
//! ```text
//!   RHS part (this BodyForce):         coupling · S̃ · Mf · v_pattern
//!   Operator diagonal (harness sums):  coupling · S̃ · I
//! ```
//!
//! The operator-diagonal part is built by
//! [`super::super::mantle::build_mantle_diagonal_field`] and is
//! added to `drag_diag` at the harness level to form a single
//! `total_diag` passed through `solve_sheet` unchanged. The inner
//! CG does not see any mantle-specific code path — it operates
//! on `A(v;η) + total_diag · I`, an SPD linear augmentation of
//! the Step-7 operator. Newton's own outer convergence drives
//! `v_solved` to the self-consistent balance `f_mantle^* =
//! coupling · S̃ · (Mf · v_pattern − v_solved^*)`.
//!
//! This module only contributes the **constant RHS part**: the
//! `v_solved`-independent contribution `coupling · S̃ · Mf · v_pattern`.
//!
//! Discretisation: `v_pattern` is already at MAC face centres
//! (via [`super::super::mantle::build_mantle_pattern`]); `S̃` is
//! cell-centered and is averaged arithmetically to the face
//! before multiplication, mirroring the GPE pattern (Step 2 D5):
//! `S̃` here plays the role of a *multiplicative weight on a source
//! term*, not a *coefficient of an elliptic operator*. Arithmetic
//! averaging is the correct interpolation for this role; harmonic
//! averaging is reserved for coefficients where it preserves flux
//! continuity across jumps (η at corners, Step 1).
//!
//! Mean-subtraction is **not** performed. A non-zero mean of
//! `f_mantle` — possible when the Fourier basis rounds — is
//! projected out by the null-space-aware preconditioner on `v`
//! (Step 0 `stokes/nullspace.rs`). Subtracting `mean(f_mantle)`
//! here would discard information the solver is equipped to
//! handle. The `v2_mantle_null_space_preservation` test checks
//! that `|mean(v_solved)|` stays below `1e-15` after the solve.

use crate::tectonics_v2::field::Field2D;
use super::super::mantle::MantlePattern;
use crate::tectonics_v2::forcing::body_force::{BodyForce, SimulationState, VectorField};

/// Constant-RHS contribution of the mantle forcing.
///
/// `v_solved` does **not** appear here: its coupling to the
/// mantle force is absorbed into the operator diagonal by the
/// harness (see module doc).
pub struct MantleForce<'a> {
    pub mf: f64,
    pub coupling: f64,
    pub pattern: &'a MantlePattern,
    /// Cell-centered `S̃`. Averaged to faces at assembly.
    pub s_field: &'a Field2D,
}

impl<'a> MantleForce<'a> {
    pub fn new(mf: f64, coupling: f64, pattern: &'a MantlePattern, s_field: &'a Field2D) -> Self {
        debug_assert_eq!(s_field.nx(), pattern.nx());
        debug_assert_eq!(s_field.ny(), pattern.ny());
        Self { mf, coupling, pattern, s_field }
    }
}

impl<'a> BodyForce for MantleForce<'a> {
    fn accumulate(&self, state: &SimulationState, out: &mut VectorField) {
        let nx = state.nx;
        let ny = state.ny;
        debug_assert_eq!(self.pattern.nx(), nx);
        debug_assert_eq!(self.pattern.ny(), ny);
        let scale = self.coupling * self.mf;

        // x-component: vx face (i, j+½) lies between cells (i-1, j) and (i, j).
        // S̃_face = ½(S̃[i-1, j] + S̃[i, j]).
        for j in 0..ny {
            for i in 0..nx {
                let im = state.idx_x.prev(i);
                let s_face = 0.5 * (self.s_field.get(im, j) + self.s_field.get(i, j));
                let vx_pat = self.pattern.v_mantle_x.get(i, j);
                let k = j * nx + i;
                out.fx.data_mut()[k] += scale * s_face * vx_pat;
            }
        }
        // y-component: vy face (i+½, j) lies between cells (i, j-1) and (i, j).
        for j in 0..ny {
            for i in 0..nx {
                let jm = state.idx_y.prev(j);
                let s_face = 0.5 * (self.s_field.get(i, jm) + self.s_field.get(i, j));
                let vy_pat = self.pattern.v_mantle_y.get(i, j);
                let k = j * nx + i;
                out.fy.data_mut()[k] += scale * s_face * vy_pat;
            }
        }
    }

    fn name(&self) -> &'static str { "MantleForce" }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics_v2::field::PeriodicIndex;
    use super::super::super::mantle::{
        build_mantle_pattern, generate_stream_function, StreamFunctionConfig,
    };

    fn env(n: usize) -> (PeriodicIndex, PeriodicIndex, Field2D, f64) {
        let dx = 1.0 / n as f64;
        (PeriodicIndex::new(n), PeriodicIndex::new(n), Field2D::filled(n, n, 1.0), dx)
    }

    fn build_pat(n: usize, seed: u64) -> (MantlePattern, PeriodicIndex, PeriodicIndex, f64) {
        let dx = 1.0 / n as f64;
        let idx_x = PeriodicIndex::new(n);
        let idx_y = PeriodicIndex::new(n);
        let psi = generate_stream_function(n, n, &StreamFunctionConfig { num_modes: 4, seed });
        let pat = build_mantle_pattern(&psi, dx, dx, &idx_x, &idx_y);
        (pat, idx_x, idx_y, dx)
    }

    /// Mf = 0 ⇒ f_mantle = 0.
    #[test]
    fn zero_mf_gives_zero_force() {
        let n = 16;
        let (pat, idx_x, idx_y, dx) = build_pat(n, 42);
        let s = Field2D::filled(n, n, 1.0);
        let mut fx = Field2D::new(n, n);
        let mut fy = Field2D::new(n, n);
        let state = SimulationState { nx: n, ny: n, dx, dy: dx, idx_x: &idx_x, idx_y: &idx_y, s: &s };
        MantleForce::new(0.0, 1.5, &pat, &s).accumulate(
            &state,
            &mut VectorField { fx: &mut fx, fy: &mut fy },
        );
        for v in fx.data().iter().chain(fy.data().iter()) {
            assert_eq!(*v, 0.0);
        }
    }

    /// coupling = 0 ⇒ f_mantle = 0.
    #[test]
    fn zero_coupling_gives_zero_force() {
        let n = 16;
        let (pat, idx_x, idx_y, dx) = build_pat(n, 42);
        let s = Field2D::filled(n, n, 1.0);
        let mut fx = Field2D::new(n, n);
        let mut fy = Field2D::new(n, n);
        let state = SimulationState { nx: n, ny: n, dx, dy: dx, idx_x: &idx_x, idx_y: &idx_y, s: &s };
        MantleForce::new(1.0, 0.0, &pat, &s).accumulate(
            &state,
            &mut VectorField { fx: &mut fx, fy: &mut fy },
        );
        for v in fx.data().iter().chain(fy.data().iter()) {
            assert_eq!(*v, 0.0);
        }
    }

    /// Linearity in Mf.
    #[test]
    fn scales_linearly_with_mf() {
        let n = 16;
        let (pat, idx_x, idx_y, dx) = build_pat(n, 42);
        let s = Field2D::filled(n, n, 1.0);
        let state = SimulationState { nx: n, ny: n, dx, dy: dx, idx_x: &idx_x, idx_y: &idx_y, s: &s };

        let mut fx1 = Field2D::new(n, n);
        let mut fy1 = Field2D::new(n, n);
        MantleForce::new(1.0, 1.0, &pat, &s).accumulate(
            &state,
            &mut VectorField { fx: &mut fx1, fy: &mut fy1 },
        );
        let mut fx3 = Field2D::new(n, n);
        let mut fy3 = Field2D::new(n, n);
        MantleForce::new(3.0, 1.0, &pat, &s).accumulate(
            &state,
            &mut VectorField { fx: &mut fx3, fy: &mut fy3 },
        );
        for k in 0..n * n {
            assert!((fx3.data()[k] - 3.0 * fx1.data()[k]).abs() < 1e-14);
            assert!((fy3.data()[k] - 3.0 * fy1.data()[k]).abs() < 1e-14);
        }
    }

    /// Uniform `S̃` and uniform `v_pattern` aren't both achievable
    /// at once on a periodic grid (v_pattern is div-free with
    /// zero mean), so the most we can cheaply check with analytic
    /// expectations is: uniform `S̃ = s₀`, the face-averaged `S̃`
    /// is also `s₀`, so `fx = coupling · s₀ · Mf · v_mantle_x`
    /// exactly at every face.
    #[test]
    fn uniform_s_factors_cleanly() {
        let n = 16;
        let (pat, idx_x, idx_y, dx) = build_pat(n, 42);
        let s0 = 0.7;
        let s = Field2D::filled(n, n, s0);
        let mut fx = Field2D::new(n, n);
        let mut fy = Field2D::new(n, n);
        let state = SimulationState { nx: n, ny: n, dx, dy: dx, idx_x: &idx_x, idx_y: &idx_y, s: &s };
        let mf = 1.5;
        let c = 2.0;
        MantleForce::new(mf, c, &pat, &s).accumulate(
            &state,
            &mut VectorField { fx: &mut fx, fy: &mut fy },
        );
        let factor = c * s0 * mf;
        for k in 0..n * n {
            let expected_x = factor * pat.v_mantle_x.data()[k];
            let expected_y = factor * pat.v_mantle_y.data()[k];
            assert!((fx.data()[k] - expected_x).abs() < 1e-14);
            assert!((fy.data()[k] - expected_y).abs() < 1e-14);
        }
    }

    /// BodyForce::accumulate contract is additive: calling twice
    /// doubles the output.
    #[test]
    fn accumulation_is_additive() {
        let n = 8;
        let (idx_x, idx_y, s, dx) = env(n);
        let (pat, _, _, _) = build_pat(n, 42);
        let state = SimulationState { nx: n, ny: n, dx, dy: dx, idx_x: &idx_x, idx_y: &idx_y, s: &s };
        let force = MantleForce::new(1.0, 1.0, &pat, &s);
        let mut fx1 = Field2D::new(n, n);
        let mut fy1 = Field2D::new(n, n);
        force.accumulate(&state, &mut VectorField { fx: &mut fx1, fy: &mut fy1 });
        let mut fx2 = Field2D::new(n, n);
        let mut fy2 = Field2D::new(n, n);
        {
            let mut out = VectorField { fx: &mut fx2, fy: &mut fy2 };
            force.accumulate(&state, &mut out);
            force.accumulate(&state, &mut out);
        }
        for k in 0..n * n {
            assert!((fx2.data()[k] - 2.0 * fx1.data()[k]).abs() < 1e-14);
            assert!((fy2.data()[k] - 2.0 * fy1.data()[k]).abs() < 1e-14);
        }
    }
}
