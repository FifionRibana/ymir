//! Slab-mass state and its two update operators.
//!
//! [`SlabState`] wraps the cell-centered scalar field `m̃(x, t)`.
//! It exposes two operators that run at well-defined points in the
//! time-step loop:
//!
//! - [`SlabState::step_ode`] — forward-Euler integration of
//!   `∂m̃/∂t̃ = Q̃_sub_conv − m̃/τ̃_slab`. **Before** the Stokes
//!   solve, so `f_slab` sees the updated `m`.
//! - [`SlabState::advect`] — conservative upwind advection with
//!   the solved velocity. **After** the Stokes solve. Reuses the
//!   Step 0 [`super::super::advection::step_upwind`] scheme
//!   unchanged — `m̃` is an extensive surface density just like
//!   `S̃`.
//!
//! Neither operator clamps or smooths the field. Any saturation
//! in `f_slab` must come from the ODE balance `Q_sub_conv ≈ m/τ`,
//! not from a numerical limiter; the Step 7 spec explicitly rejects
//! the legacy `max_plate_velocity` pattern (§D2) and the
//! smoothing would hide under-/over-shooting from the diagnostic
//! reports.

use super::super::advection::step_upwind;
use super::super::field::{Field2D, PeriodicIndex};

/// Slab-mass field wrapper.
///
/// Owns a single `Field2D` sized `nx × ny`. Advection needs a
/// temporary buffer of the same shape, which we own internally so
/// the caller does not have to plumb it through every step.
#[derive(Clone)]
pub struct SlabState {
    nx: usize,
    ny: usize,
    m: Field2D,
    /// Scratch buffer for conservative advection. Kept here to
    /// avoid a per-step allocation in the time loop.
    scratch: Field2D,
}

impl std::fmt::Debug for SlabState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Field2D does not derive Debug (legacy, keeps the newton
        // transcript readable). Dump the shape and a summary of the
        // slab-mass field so `{:?}` on harness state stays useful.
        let m_data = self.m.data();
        let m_min = m_data.iter().cloned().fold(f64::INFINITY, f64::min);
        let m_max = m_data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let m_sum: f64 = m_data.iter().sum();
        f.debug_struct("SlabState")
            .field("nx", &self.nx)
            .field("ny", &self.ny)
            .field("m_min", &m_min)
            .field("m_max", &m_max)
            .field("m_integrated", &m_sum)
            .finish()
    }
}

impl SlabState {
    /// Construct with `m ≡ 0` everywhere. No mass has been
    /// subducted at `t = 0`; the field fills in dynamically.
    pub fn new_zero(nx: usize, ny: usize) -> Self {
        Self { nx, ny, m: Field2D::new(nx, ny), scratch: Field2D::new(nx, ny) }
    }

    pub fn nx(&self) -> usize {
        self.nx
    }
    pub fn ny(&self) -> usize {
        self.ny
    }

    /// Read-only view of the slab-mass field.
    #[inline]
    pub fn m(&self) -> &Field2D {
        &self.m
    }

    /// Mutable view, used in tests that seed specific `m`
    /// distributions before exercising an operator.
    #[inline]
    pub fn m_mut(&mut self) -> &mut Field2D {
        &mut self.m
    }

    /// Apply the forward-Euler ODE step
    /// `m(t+Δt) = m(t) + Δt · [Q_sub_conv − m(t)/τ_slab]`.
    ///
    /// Stability requires `Δt < τ_slab`. The harness enforces a
    /// defensive `Δt ≤ 0.1 · τ_slab` which sits well above the
    /// CFL-driven bound in practice (`Δt ≈ 0.005`).
    ///
    /// Panics in debug builds if `τ_slab ≤ 0` (division by zero
    /// or nonsensical negative decay).
    pub fn step_ode(&mut self, q_sub_conv: &Field2D, dt: f64, tau_slab: f64) {
        debug_assert_eq!(q_sub_conv.nx(), self.nx);
        debug_assert_eq!(q_sub_conv.ny(), self.ny);
        debug_assert!(tau_slab > 0.0, "τ_slab must be positive, got {}", tau_slab);

        let inv_tau = 1.0 / tau_slab;
        let n = self.nx * self.ny;
        let m_data = self.m.data_mut();
        let q_data = q_sub_conv.data();
        for k in 0..n {
            let m_old = m_data[k];
            let rhs = q_data[k] - m_old * inv_tau;
            m_data[k] = m_old + dt * rhs;
        }
    }

    /// Conservative upwind advection of `m` with the MAC-face
    /// velocity field `(vx, vy)`. Preserves the total slab mass
    /// to machine precision on a periodic domain (see the
    /// `v2_slab_advection_mms` test).
    ///
    /// The caller supplies the same periodic indexers they use
    /// for `S̃` — we reuse the Step 0 scheme unchanged.
    pub fn advect(
        &mut self,
        dx: f64,
        dy: f64,
        dt: f64,
        idx_x: &PeriodicIndex,
        idx_y: &PeriodicIndex,
        vx: &[f64],
        vy: &[f64],
    ) {
        step_upwind(self.nx, self.ny, dx, dy, dt, idx_x, idx_y, &self.m, vx, vy, &mut self.scratch);
        // Swap the buffers: scratch now holds m(t+Δt), the old m
        // becomes the next scratch. One pointer swap, no data
        // copy.
        std::mem::swap(&mut self.m, &mut self.scratch);
    }

    /// Integrated slab mass `Σ m̃` over the domain. Diagnostic
    /// helper for conservation checks and reporting.
    pub fn integrated(&self) -> f64 {
        self.m.data().iter().sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Constant `Q_sub_conv = Q₀`, zero velocity: the ODE
    /// trajectory is the analytic `m(t) = Q₀·τ·(1 − e^{-t/τ})`.
    /// At `t = 5·τ` we should be within 1% of the asymptote
    /// `Q₀·τ`.
    #[test]
    fn ode_converges_to_q_times_tau() {
        let nx = 8;
        let ny = 8;
        let tau = 0.5;
        let dt = 0.01;
        let q0 = 0.3;

        let mut state = SlabState::new_zero(nx, ny);
        let mut q = Field2D::new(nx, ny);
        for v in q.data_mut().iter_mut() {
            *v = q0;
        }

        // Integrate to t = 5·τ = 2.5. 250 steps at dt=0.01.
        for _ in 0..250 {
            state.step_ode(&q, dt, tau);
        }

        let asymptote = q0 * tau;
        for &v in state.m().data().iter() {
            let rel = (v - asymptote).abs() / asymptote;
            assert!(rel < 0.01, "m = {}, expected ≈ {}, rel err = {}", v, asymptote, rel);
        }
    }

    /// With `Q_sub_conv = 0` and non-zero `m`, the ODE is pure
    /// exponential decay. Check half-life ≈ τ·ln(2).
    #[test]
    fn ode_half_life_matches_tau_ln2() {
        let nx = 4;
        let ny = 4;
        let tau = 0.5;
        let dt = 0.001;
        let m0 = 1.0;

        let mut state = SlabState::new_zero(nx, ny);
        for v in state.m_mut().data_mut().iter_mut() {
            *v = m0;
        }
        let q = Field2D::new(nx, ny); // zero

        // Step until half. τ·ln2 = 0.3465... ≈ 347 steps.
        let mut step = 0;
        while state.m().data()[0] > 0.5 * m0 && step < 1000 {
            state.step_ode(&q, dt, tau);
            step += 1;
        }
        let observed_half_life = step as f64 * dt;
        let analytic = tau * std::f64::consts::LN_2;
        let rel = (observed_half_life - analytic).abs() / analytic;
        assert!(
            rel < 0.05,
            "observed half-life = {} τ, expected {} τ·ln2 (rel err {})",
            observed_half_life,
            analytic,
            rel,
        );
    }

    /// Zero velocity advection must leave `m` untouched to
    /// machine precision.
    #[test]
    fn zero_velocity_advection_is_identity() {
        let nx = 8;
        let ny = 8;
        let dx = 1.0 / nx as f64;
        let dt = 0.01;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);

        let mut state = SlabState::new_zero(nx, ny);
        for (k, v) in state.m_mut().data_mut().iter_mut().enumerate() {
            *v = 0.1 + 0.01 * k as f64;
        }
        let snapshot: Vec<f64> = state.m().data().to_vec();

        let vx = vec![0.0; nx * ny];
        let vy = vec![0.0; nx * ny];

        state.advect(dx, dx, dt, &idx_x, &idx_y, &vx, &vy);

        for (a, b) in snapshot.iter().zip(state.m().data().iter()) {
            assert!((a - b).abs() < 1e-15);
        }
    }

    /// Uniform velocity advection preserves the integrated mass
    /// exactly (upwind on a periodic domain is conservative).
    #[test]
    fn uniform_translation_preserves_total_mass() {
        let nx = 16;
        let ny = 16;
        let dx = 1.0 / nx as f64;
        let dt = 0.01;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);

        let mut state = SlabState::new_zero(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                let x = (i as f64 + 0.5) * dx;
                let y = (j as f64 + 0.5) * dx;
                let v =
                    (2.0 * std::f64::consts::PI * x).sin() * (2.0 * std::f64::consts::PI * y).cos();
                state.m_mut().set(i, j, 0.5 + 0.2 * v);
            }
        }
        let initial = state.integrated();

        let vx = vec![0.3; nx * ny];
        let vy = vec![0.2; nx * ny];
        for _ in 0..50 {
            state.advect(dx, dx, dt, &idx_x, &idx_y, &vx, &vy);
        }
        let final_mass = state.integrated();
        let drift = (final_mass - initial).abs() / initial.abs();
        assert!(drift < 1e-12, "mass drift {} exceeds bound", drift);
    }
}
