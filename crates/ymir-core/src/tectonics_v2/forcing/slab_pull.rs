//! Slab-pull body force (Step 7).
//!
//! Continuous form (nondim, §4.8):
//! ```text
//!   f̃_slab(x) = Sp · m̃(x) · n̂_convergence(x)
//! ```
//! where `m̃` is the cell-centered slab-mass field maintained by
//! [`super::super::slab::state::SlabState`] and `n̂_convergence` is
//! the cell-centered unit vector produced by
//! [`super::super::slab::convergence_direction::compute_convergence_direction`].
//!
//! Discretisation on the MAC grid — mirror of `GpeForce`.
//! Both operands (`m̃` and `n̂`) live at cell centres; the
//! momentum face `fx[i,j]` at `(i·dx, (j+½)·dy)` sees arithmetic
//! averages of the two cells it bridges:
//!
//! ```text
//!   fx[i,j] += Sp · ½(m[i-1,j] + m[i,j]) · ½(n_x[i-1,j] + n_x[i,j])
//! ```
//! and symmetrically `fy[i,j] +=
//! Sp · ½(m[i,j-1] + m[i,j]) · ½(n_y[i,j-1] + n_y[i,j])`.
//!
//! Why centred averaging. `m̃` is extensive and smooth after
//! advection + decay (no sharp fronts; the drain `-m/τ` smooths
//! with time-scale `τ`). `n̂` is a direction vector bounded on
//! `[-1, 1]`, non-smooth only near epsilon-fallback boundaries,
//! but centred averaging is the only face-interpolation that
//! preserves the rotation-invariance required by the null-space
//! argument below.
//!
//! **No mean subtraction.** A spatially asymmetric slab
//! configuration can leave `mean(f_slab) ≠ 0`. This is a
//! physical signal (net convergent forcing over the tore). The
//! null-space projector installed in the preconditioner `M⁻¹`
//! at Step 0 correctly removes the null-space component of the
//! solved `v`, not of `f`. Subtracting `mean(f_slab)` preemptively
//! would discard information the solver is already equipped to
//! handle. See the `v2_slab_null_space_preservation` test.

use super::super::field::Field2D;
use super::body_force::{BodyForce, SimulationState, VectorField};

/// Slab-pull body-force term.
///
/// Holds references to per-step state (`m_subducted`,
/// `convergence_direction_x/y`) rebuilt by the harness before
/// each Stokes solve. The `Sp` coupling is a scalar copy.
pub struct SlabPullForce<'a> {
    pub sp: f64,
    pub m_subducted: &'a Field2D,
    pub n_x: &'a Field2D,
    pub n_y: &'a Field2D,
}

impl<'a> SlabPullForce<'a> {
    pub fn new(sp: f64, m_subducted: &'a Field2D, n_x: &'a Field2D, n_y: &'a Field2D) -> Self {
        debug_assert_eq!(m_subducted.nx(), n_x.nx());
        debug_assert_eq!(m_subducted.ny(), n_x.ny());
        debug_assert_eq!(m_subducted.nx(), n_y.nx());
        debug_assert_eq!(m_subducted.ny(), n_y.ny());
        Self { sp, m_subducted, n_x, n_y }
    }
}

impl<'a> BodyForce for SlabPullForce<'a> {
    fn accumulate(&self, state: &SimulationState, out: &mut VectorField) {
        let nx = state.nx;
        let ny = state.ny;
        debug_assert_eq!(self.m_subducted.nx(), nx);
        debug_assert_eq!(self.m_subducted.ny(), ny);

        let sp = self.sp;

        // x-component: face (i, j+½) sees cells (i-1, j) and (i, j).
        for j in 0..ny {
            for i in 0..nx {
                let im = state.idx_x.prev(i);
                let m_face = 0.5 * (self.m_subducted.get(im, j) + self.m_subducted.get(i, j));
                let nx_face = 0.5 * (self.n_x.get(im, j) + self.n_x.get(i, j));
                let k = j * nx + i;
                out.fx.data_mut()[k] += sp * m_face * nx_face;
            }
        }
        // y-component: face (i+½, j) sees cells (i, j-1) and (i, j).
        for j in 0..ny {
            for i in 0..nx {
                let jm = state.idx_y.prev(j);
                let m_face = 0.5 * (self.m_subducted.get(i, jm) + self.m_subducted.get(i, j));
                let ny_face = 0.5 * (self.n_y.get(i, jm) + self.n_y.get(i, j));
                let k = j * nx + i;
                out.fy.data_mut()[k] += sp * m_face * ny_face;
            }
        }
    }

    fn name(&self) -> &'static str {
        "SlabPullForce"
    }
}

#[cfg(test)]
mod tests {
    use super::super::super::field::PeriodicIndex;
    use super::*;

    fn env(
        nx: usize,
        ny: usize,
    ) -> (PeriodicIndex, PeriodicIndex, Field2D, Field2D, Field2D, Field2D, Field2D, Field2D) {
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);
        let s = Field2D::filled(nx, ny, 1.0);
        let m = Field2D::new(nx, ny);
        let n_x = Field2D::new(nx, ny);
        let n_y = Field2D::new(nx, ny);
        let fx = Field2D::new(nx, ny);
        let fy = Field2D::new(nx, ny);
        (idx_x, idx_y, s, m, n_x, n_y, fx, fy)
    }

    /// `m = 0` everywhere ⇒ f_slab = 0 everywhere, regardless of n̂.
    #[test]
    fn zero_mass_gives_zero_force() {
        let nx = 8;
        let ny = 8;
        let (idx_x, idx_y, s, m, mut n_x, mut n_y, mut fx, mut fy) = env(nx, ny);
        // Arbitrary non-zero n̂.
        for k in 0..nx * ny {
            n_x.data_mut()[k] = 0.7;
            n_y.data_mut()[k] = -0.3;
        }
        let st =
            SimulationState { nx, ny, dx: 0.125, dy: 0.125, idx_x: &idx_x, idx_y: &idx_y, s: &s };
        let force = SlabPullForce::new(2.0, &m, &n_x, &n_y);
        force.accumulate(&st, &mut VectorField { fx: &mut fx, fy: &mut fy });
        for v in fx.data().iter().chain(fy.data().iter()) {
            assert_eq!(*v, 0.0);
        }
    }

    /// `n̂ = 0` everywhere ⇒ f_slab = 0 regardless of m.
    #[test]
    fn zero_direction_gives_zero_force() {
        let nx = 8;
        let ny = 8;
        let (idx_x, idx_y, s, mut m, n_x, n_y, mut fx, mut fy) = env(nx, ny);
        for v in m.data_mut().iter_mut() {
            *v = 0.5;
        }
        let st =
            SimulationState { nx, ny, dx: 0.125, dy: 0.125, idx_x: &idx_x, idx_y: &idx_y, s: &s };
        let force = SlabPullForce::new(1.0, &m, &n_x, &n_y);
        force.accumulate(&st, &mut VectorField { fx: &mut fx, fy: &mut fy });
        for v in fx.data().iter().chain(fy.data().iter()) {
            assert_eq!(*v, 0.0);
        }
    }

    /// Uniform `m` and `n̂ = (1, 0)` everywhere: the x-component
    /// of the force equals `Sp · m · 1 = Sp · m` on every face.
    /// y-component is 0.
    #[test]
    fn uniform_fields_match_closed_form() {
        let nx = 8;
        let ny = 8;
        let (idx_x, idx_y, s, mut m, mut n_x, n_y, mut fx, mut fy) = env(nx, ny);
        for v in m.data_mut().iter_mut() {
            *v = 0.4;
        }
        for v in n_x.data_mut().iter_mut() {
            *v = 1.0;
        }
        let sp = 2.5;
        let st =
            SimulationState { nx, ny, dx: 0.125, dy: 0.125, idx_x: &idx_x, idx_y: &idx_y, s: &s };
        let force = SlabPullForce::new(sp, &m, &n_x, &n_y);
        force.accumulate(&st, &mut VectorField { fx: &mut fx, fy: &mut fy });
        let expected_fx = sp * 0.4 * 1.0;
        for v in fx.data().iter() {
            assert!((*v - expected_fx).abs() < 1e-14);
        }
        for v in fy.data().iter() {
            assert_eq!(*v, 0.0);
        }
    }

    /// Linearity in `Sp`.
    #[test]
    fn scales_linearly_with_sp() {
        let nx = 8;
        let ny = 8;
        let (idx_x, idx_y, s, mut m, mut n_x, mut n_y, mut fx1, mut fy1) = env(nx, ny);
        for (k, v) in m.data_mut().iter_mut().enumerate() {
            *v = 0.05 * (k as f64 + 1.0);
        }
        for (k, v) in n_x.data_mut().iter_mut().enumerate() {
            *v = (k as f64 * 0.2).sin();
        }
        for (k, v) in n_y.data_mut().iter_mut().enumerate() {
            *v = (k as f64 * 0.2).cos();
        }
        let st =
            SimulationState { nx, ny, dx: 0.125, dy: 0.125, idx_x: &idx_x, idx_y: &idx_y, s: &s };

        SlabPullForce::new(1.0, &m, &n_x, &n_y)
            .accumulate(&st, &mut VectorField { fx: &mut fx1, fy: &mut fy1 });
        let mut fx3 = Field2D::new(nx, ny);
        let mut fy3 = Field2D::new(nx, ny);
        SlabPullForce::new(3.0, &m, &n_x, &n_y)
            .accumulate(&st, &mut VectorField { fx: &mut fx3, fy: &mut fy3 });
        for k in 0..nx * ny {
            assert!((fx3.data()[k] - 3.0 * fx1.data()[k]).abs() < 1e-14);
            assert!((fy3.data()[k] - 3.0 * fy1.data()[k]).abs() < 1e-14);
        }
    }

    /// `BodyForce` contract: accumulate two identical calls →
    /// 2× the single call. Essential so `ForceSum` composes
    /// correctly with `GpeForce`.
    #[test]
    fn accumulation_is_additive() {
        let nx = 8;
        let ny = 8;
        let (idx_x, idx_y, s, mut m, mut n_x, mut n_y, _fx, _fy) = env(nx, ny);
        for (k, v) in m.data_mut().iter_mut().enumerate() {
            *v = 0.1 + 0.01 * k as f64;
        }
        for (k, v) in n_x.data_mut().iter_mut().enumerate() {
            *v = (k as f64 * 0.1).sin();
        }
        for (k, v) in n_y.data_mut().iter_mut().enumerate() {
            *v = (k as f64 * 0.1).cos();
        }
        let sp = 1.5;
        let st =
            SimulationState { nx, ny, dx: 0.125, dy: 0.125, idx_x: &idx_x, idx_y: &idx_y, s: &s };

        let force = SlabPullForce::new(sp, &m, &n_x, &n_y);
        let mut fx_once = Field2D::new(nx, ny);
        let mut fy_once = Field2D::new(nx, ny);
        force.accumulate(&st, &mut VectorField { fx: &mut fx_once, fy: &mut fy_once });

        let mut fx_twice = Field2D::new(nx, ny);
        let mut fy_twice = Field2D::new(nx, ny);
        {
            let mut out = VectorField { fx: &mut fx_twice, fy: &mut fy_twice };
            force.accumulate(&st, &mut out);
            force.accumulate(&st, &mut out);
        }

        for k in 0..nx * ny {
            assert!((fx_twice.data()[k] - 2.0 * fx_once.data()[k]).abs() < 1e-14);
            assert!((fy_twice.data()[k] - 2.0 * fy_once.data()[k]).abs() < 1e-14);
        }
    }
}
