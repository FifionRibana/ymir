//! Step 10 — non-conservative upwind advection for the age field
//! `A`.
//!
//! # Why a different scheme from `S̃`?
//!
//! The `S̃` advection in [`crate::tectonics_v2::advection::step_upwind`]
//! is **conservative** (flux form):
//!
//! ```text
//!   ∂_t S̃ + ∇·(S̃ ṽ) = 0
//! ```
//!
//! This is correct for `S̃` because it is a thickness (a mass-like
//! density): integrated `S̃` must be conserved. The conservative
//! discretisation guarantees machine-precision conservation even
//! when `∇·ṽ ≠ 0`.
//!
//! The age field `A` follows the **non-conservative** equation
//! (per `solver-scaling.md` §4.11):
//!
//! ```text
//!   ∂A/∂t + ṽ · ∇A = Γ
//! ```
//!
//! `A` is intensive — it is a per-cell scalar (age in `τ*`), not a
//! mass density. A cell moved by plate motion **keeps** its age;
//! there is no notion of "spreading the age over multiple cells"
//! the way mass spreads. Reusing `step_upwind` (flux form) on `A`
//! would expand to
//!
//! ```text
//!   ∂A/∂t + ṽ · ∇A + A · ∇·ṽ = 0
//! ```
//!
//! introducing a spurious `A · ∇·ṽ` source whenever the velocity
//! has non-zero divergence. With mantle forcing or boundary
//! sources making `∇·ṽ` regularly `O(1)` and `A` ranging up to
//! `O(10)`, the spurious term would add `O(10) · ∇·ṽ · dt` per
//! cell per step — dominant over the intended quiescent growth of
//! `1 · dt`.
//!
//! Step 10 therefore implements a **separate** non-conservative
//! upwind scheme for `A`, deviating from the issue's note ("OR
//! call the existing function on A as a separate scalar after S̃
//! advection is complete. The latter is simpler [...]"). The
//! issue's recommendation rests on the implicit assumption
//! `∇·ṽ = 0`, which does not hold under Step 8+ mantle forcing.
//! The §4.11 patch should be updated to clarify the
//! conservative-vs-Lagrangian distinction.
//!
//! # Discretisation
//!
//! For each cell `(i, j)`:
//!
//! 1. Compute cell-centred velocity by averaging the staggered
//!    face velocities:
//!    `vx_c = ½ (vx[i,j] + vx[ip,j])`, `vy_c = ½ (vy[i,j] + vy[i,jp])`.
//! 2. Apply upwind one-sided differences:
//!    `dA/dx = (A[i,j] - A[im,j]) / dx` if `vx_c > 0`, else
//!    `dA/dx = (A[ip,j] - A[i,j]) / dx`.
//!    Likewise for `dA/dy`.
//! 3. Update:
//!    `A_next[i,j] = A[i,j] - dt · (vx_c · dA/dx + vy_c · dA/dy) + dt · Γ_quiescent`.
//!
//! `Γ_quiescent = 1.0` for every cell — this is the "age grows
//! linearly with t" term from §4.11. Boundary-event resets
//! (ridge / arc / collision) overwrite specific cells *after* the
//! advection step in [`super::events::apply_age_events`].
//!
//! # CFL bound
//!
//! Same as the conservative scheme:
//! `dt ≤ cfl_factor · min(dx, dy) / max|v|`. Reuse
//! [`crate::tectonics_v2::advection::cfl_dt`] — caller responsible
//! for honouring the bound.

use crate::tectonics_v2::field::{Field2D, PeriodicIndex};

/// Source term `Γ` per §4.11 quiescent growth: every cell ages at
/// `dA/dt = 1` in the Lagrangian frame, achieved by a uniform
/// source of `+dt` per step before event-driven resets overwrite
/// boundary cells.
pub const QUIESCENT_GROWTH_RATE: f64 = 1.0;

/// Advance `A` by one forward-Euler step using **non-conservative**
/// first-order upwind for the advection term `v · ∇A` plus a
/// uniform source `Γ_quiescent · dt` representing Lagrangian age
/// growth. Output written into `a_next`.
pub fn step_age_advect(
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    dt: f64,
    idx_x: &PeriodicIndex,
    idx_y: &PeriodicIndex,
    a: &Field2D,
    vx: &[f64],
    vy: &[f64],
    a_next: &mut Field2D,
) {
    debug_assert_eq!(vx.len(), nx * ny);
    debug_assert_eq!(vy.len(), nx * ny);
    debug_assert_eq!(a.nx(), nx);
    debug_assert_eq!(a.ny(), ny);
    debug_assert_eq!(a_next.nx(), nx);
    debug_assert_eq!(a_next.ny(), ny);
    let inv_dx = 1.0 / dx;
    let inv_dy = 1.0 / dy;
    let lin = |ii: usize, jj: usize| jj * nx + ii;

    for j in 0..ny {
        let jp = idx_y.next(j);
        let jm = idx_y.prev(j);
        for i in 0..nx {
            let ip = idx_x.next(i);
            let im = idx_x.prev(i);

            // Cell-centred velocity = average of staggered face
            // velocities. This is the "Lagrangian observer's"
            // velocity at the cell centre.
            let vx_c = 0.5 * (vx[lin(i, j)] + vx[lin(ip, j)]);
            let vy_c = 0.5 * (vy[lin(i, j)] + vy[lin(i, jp)]);

            // Upwind one-sided differences for the gradient.
            let da_dx = if vx_c >= 0.0 {
                (a.get(i, j) - a.get(im, j)) * inv_dx
            } else {
                (a.get(ip, j) - a.get(i, j)) * inv_dx
            };
            let da_dy = if vy_c >= 0.0 {
                (a.get(i, j) - a.get(i, jm)) * inv_dy
            } else {
                (a.get(i, jp) - a.get(i, j)) * inv_dy
            };

            let advect = vx_c * da_dx + vy_c * da_dy;
            let next = a.get(i, j) - dt * advect + dt * QUIESCENT_GROWTH_RATE;
            a_next.set(i, j, next);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * a.abs().max(b.abs()).max(1.0)
    }

    #[test]
    fn zero_velocity_grows_uniformly_at_one_per_dt() {
        // The quiescent-growth acceptance #3: with v = 0, A grows
        // linearly at +dt per step, regardless of initial value.
        let nx = 4;
        let ny = 4;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);
        let mut a = Field2D::filled(nx, ny, 3.0);
        let vx = vec![0.0; nx * ny];
        let vy = vec![0.0; nx * ny];
        let dt = 0.1_f64;
        let mut a_next = Field2D::new(nx, ny);
        for k in 1..=10 {
            step_age_advect(nx, ny, 1.0, 1.0, dt, &idx_x, &idx_y, &a, &vx, &vy, &mut a_next);
            std::mem::swap(&mut a, &mut a_next);
            let expected = 3.0 + (k as f64) * dt;
            for v in a.data() {
                assert!(
                    approx(*v, expected, 1e-12),
                    "step {}: got {}, expected {}",
                    k,
                    v,
                    expected
                );
            }
        }
    }

    #[test]
    fn uniform_field_stays_uniform_under_any_velocity() {
        // ∇A = 0 ⇒ no advective contribution. Only quiescent growth.
        let nx = 8;
        let ny = 8;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);
        let mut a = Field2D::filled(nx, ny, 5.0);
        let mut vx = vec![0.0; nx * ny];
        let mut vy = vec![0.0; nx * ny];
        for k in 0..nx * ny {
            vx[k] = ((k as f64 * 0.7).sin()) * 0.3;
            vy[k] = ((k as f64 * 1.3).cos()) * 0.2;
        }
        let dt = 0.05_f64;
        let mut a_next = Field2D::new(nx, ny);
        step_age_advect(nx, ny, 0.125, 0.125, dt, &idx_x, &idx_y, &a, &vx, &vy, &mut a_next);
        std::mem::swap(&mut a, &mut a_next);
        for v in a.data() {
            assert!(approx(*v, 5.0 + dt, 1e-12), "got {}", v);
        }
    }

    #[test]
    fn constant_velocity_advects_a_step_pulse_correctly() {
        // MMS-style: prescribe a step pulse in A, advect with
        // uniform v in +x. After dt the pulse should shift right
        // by v · dt cells. With dt small enough that the pulse
        // moves less than a cell, the upwind scheme produces a
        // small numerical-diffusion broadening but the centre of
        // mass is preserved to first order.
        let nx = 16;
        let ny = 4;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);
        let mut a = Field2D::new(nx, ny);
        // Pulse at i = 8.
        for j in 0..ny {
            a.set(8, j, 1.0);
        }
        let vx = vec![1.0; nx * ny];
        let vy = vec![0.0; nx * ny];
        let dx = 0.1_f64;
        let dt = 0.05_f64; // CFL = 0.5 (vx · dt / dx)
        let mut a_next = Field2D::new(nx, ny);
        // Run a few steps and check the pulse drifts right.
        for _step in 0..3 {
            step_age_advect(nx, ny, dx, dx, dt, &idx_x, &idx_y, &a, &vx, &vy, &mut a_next);
            std::mem::swap(&mut a, &mut a_next);
        }
        // After 3 steps the pulse centre of mass should have
        // moved roughly 3 · vx · dt / dx = 1.5 cells to the right.
        // We just verify the pulse is no longer at i = 8 alone:
        // A[i=9] should now be > A[i=7] (downstream broader than
        // upstream — correct upwind asymmetry).
        let row = 0;
        assert!(
            a.get(9, row) > a.get(7, row),
            "pulse did not drift downstream: a[7]={}, a[9]={}",
            a.get(7, row),
            a.get(9, row)
        );
    }

    #[test]
    fn step_age_advect_is_byte_deterministic_on_same_inputs() {
        // Acceptance #12 / D6 — the structural by-pass requires
        // that any randomness or non-determinism inside the
        // advection step would leak into the regression. Verify
        // that two calls with the same inputs produce byte-equal
        // outputs.
        let nx = 8;
        let ny = 8;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);
        let a = {
            let mut f = Field2D::new(nx, ny);
            for k in 0..nx * ny {
                f.data_mut()[k] = (k as f64 * 0.3).sin();
            }
            f
        };
        let vx: Vec<f64> = (0..nx * ny).map(|k| (k as f64 * 0.7).cos()).collect();
        let vy: Vec<f64> = (0..nx * ny).map(|k| (k as f64 * 1.1).sin()).collect();
        let dt = 0.05;

        let mut out1 = Field2D::new(nx, ny);
        let mut out2 = Field2D::new(nx, ny);
        step_age_advect(nx, ny, 0.125, 0.125, dt, &idx_x, &idx_y, &a, &vx, &vy, &mut out1);
        step_age_advect(nx, ny, 0.125, 0.125, dt, &idx_x, &idx_y, &a, &vx, &vy, &mut out2);
        for (x, y) in out1.data().iter().zip(out2.data().iter()) {
            assert_eq!(x.to_bits(), y.to_bits());
        }
    }
}
