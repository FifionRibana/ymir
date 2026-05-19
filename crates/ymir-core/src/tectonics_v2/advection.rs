//! Conservative first-order upwind advection for the passive
//! thickness field `S̃` on the MAC grid.
//!
//! The equation is
//!
//! ```text
//!   ∂_t̃ S̃ + ∇·(S̃ ṽ) = 0
//! ```
//!
//! written in flux form. The integrated mass `Σ S̃` over the periodic
//! torus is conserved to machine precision by the discretization
//! regardless of whether `ṽ` is exactly divergence-free, because the
//! fluxes leaving a cell through each face are the fluxes entering
//! its neighbour through the same face (opposite signs cancel in the
//! global sum). Higher-order upwind is deferred to a later step.
//!
//! Time integration is explicit forward Euler; the caller is
//! responsible for enforcing the CFL bound (see [`cfl_dt`]).

use super::field::{Field2D, PeriodicIndex};

/// Maximum admissible nondim time step under first-order upwind on
/// the staggered MAC grid: `dt ≤ cfl_factor · dx / max|ṽ|`.
///
/// `cfl_factor ∈ [0.3, 0.5]` is typical; 0.3 is safer with nonlinear
/// source feedback (per `solver-scaling.md` §4.6).
pub fn cfl_dt(dx: f64, dy: f64, vx: &[f64], vy: &[f64], cfl_factor: f64) -> f64 {
    let vmax = vx
        .iter()
        .chain(vy.iter())
        .fold(0.0_f64, |acc, &v| acc.max(v.abs()));
    if vmax == 0.0 {
        f64::INFINITY
    } else {
        cfl_factor * dx.min(dy) / vmax
    }
}

/// Advance `S̃` by one forward-Euler step using first-order upwind
/// conservative flux divergence. Output is written into `s_next`.
///
/// The fluxes are taken at the same staggered faces as `vx` and `vy`.
/// For each face the upwind S value is chosen based on the sign of
/// the face velocity. Formally:
///
/// ```text
///   flux_x[i,j] = vx[i,j] · S̃_up,    S̃_up = S̃[i-1, j] if vx > 0, S̃[i, j] otherwise
///   flux_y[i,j] = vy[i,j] · S̃_up,    S̃_up = S̃[i, j-1] if vy > 0, S̃[i, j] otherwise
///   S̃_next[i, j] = S̃[i, j] - dt · ((flux_x[ip, j] - flux_x[i, j]) / dx
///                                 + (flux_y[i, jp] - flux_y[i, j]) / dy)
/// ```
pub fn step_upwind(
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    dt: f64,
    idx_x: &PeriodicIndex,
    idx_y: &PeriodicIndex,
    s: &Field2D,
    vx: &[f64],
    vy: &[f64],
    s_next: &mut Field2D,
) {
    debug_assert_eq!(vx.len(), nx * ny);
    debug_assert_eq!(vy.len(), nx * ny);
    debug_assert_eq!(s.nx(), nx);
    debug_assert_eq!(s.ny(), ny);
    debug_assert_eq!(s_next.nx(), nx);
    debug_assert_eq!(s_next.ny(), ny);
    let inv_dx = 1.0 / dx;
    let inv_dy = 1.0 / dy;
    let lin = |ii: usize, jj: usize| jj * nx + ii;

    for j in 0..ny {
        let jp = idx_y.next(j);
        let jm = idx_y.prev(j);
        for i in 0..nx {
            let ip = idx_x.next(i);
            let im = idx_x.prev(i);

            let vx_left = vx[lin(i, j)];
            let vx_right = vx[lin(ip, j)];
            let vy_bot = vy[lin(i, j)];
            let vy_top = vy[lin(i, jp)];

            let s_up_left = if vx_left >= 0.0 { s.get(im, j) } else { s.get(i, j) };
            let s_up_right = if vx_right >= 0.0 { s.get(i, j) } else { s.get(ip, j) };
            let s_up_bot = if vy_bot >= 0.0 { s.get(i, jm) } else { s.get(i, j) };
            let s_up_top = if vy_top >= 0.0 { s.get(i, j) } else { s.get(i, jp) };

            let flux_x_left = vx_left * s_up_left;
            let flux_x_right = vx_right * s_up_right;
            let flux_y_bot = vy_bot * s_up_bot;
            let flux_y_top = vy_top * s_up_top;

            let dsdt = -((flux_x_right - flux_x_left) * inv_dx
                + (flux_y_top - flux_y_bot) * inv_dy);
            s_next.set(i, j, s.get(i, j) + dt * dsdt);
        }
    }
}

/// Integrated mass `Σ S̃` over the periodic domain. Used as a
/// conservation diagnostic.
pub fn integrated_mass(s: &Field2D) -> f64 {
    s.data().iter().sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mass_is_conserved_with_nonzero_divergence_velocity() {
        // Even a velocity field that is not exactly divergence-free
        // must preserve total S globally because face fluxes cancel
        // pairwise. This is the property the upwind scheme is chosen
        // for.
        let nx = 8;
        let ny = 6;
        let dx = 0.125;
        let dy = 0.1667;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);
        let mut s = Field2D::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                s.set(i, j, 1.0 + 0.3 * ((i as f64 * 0.7 + j as f64 * 0.9).sin()));
            }
        }
        // A non-divergence-free periodic velocity.
        let mut vx = vec![0.0; nx * ny];
        let mut vy = vec![0.0; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let k = j * nx + i;
                vx[k] = (i as f64 * 0.5).sin();
                vy[k] = (j as f64 * 0.8).cos();
            }
        }
        let mass0 = integrated_mass(&s);
        let dt = cfl_dt(dx, dy, &vx, &vy, 0.3);
        let mut s_next = Field2D::new(nx, ny);
        for _step in 0..50 {
            step_upwind(nx, ny, dx, dy, dt, &idx_x, &idx_y, &s, &vx, &vy, &mut s_next);
            std::mem::swap(&mut s, &mut s_next);
        }
        let mass1 = integrated_mass(&s);
        let rel_drift = (mass1 - mass0).abs() / mass0.abs().max(1.0);
        assert!(rel_drift < 1e-12, "mass drift = {}", rel_drift);
    }

    #[test]
    fn zero_velocity_leaves_s_unchanged() {
        let nx = 4;
        let ny = 4;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);
        let mut s = Field2D::new(nx, ny);
        for k in 0..16 {
            s.data_mut()[k] = (k as f64 * 1.1).sin();
        }
        let s_before = s.data().to_vec();
        let vx = vec![0.0; 16];
        let vy = vec![0.0; 16];
        let mut s_next = Field2D::new(nx, ny);
        step_upwind(nx, ny, 0.25, 0.25, 0.1, &idx_x, &idx_y, &s, &vx, &vy, &mut s_next);
        for (a, b) in s_before.iter().zip(s_next.data().iter()) {
            assert_eq!(*a, *b);
        }
    }

    #[test]
    fn cfl_handles_zero_velocity() {
        let v = vec![0.0; 10];
        let dt = cfl_dt(0.1, 0.1, &v, &v, 0.3);
        assert!(dt.is_infinite());
    }
}
