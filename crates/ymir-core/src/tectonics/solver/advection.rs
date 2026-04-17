//! Upwind advection scheme for crustal thickness transport.

use rayon::prelude::*;

use super::field::Field2D;
use super::grid::StaggeredGrid;

const PAR_THRESHOLD: usize = 64;

/// Compute the divergence of the flux ∇·(S·v) using first-order upwind.
///
/// Each flux uses the value of S from the upwind side, ensuring stability
/// and discrete conservation (sum of div over the periodic grid is exactly 0).
#[allow(clippy::needless_range_loop)]
pub fn compute_divergence_flux(grid: &StaggeredGrid, div: &mut Field2D) {
    let nx = grid.nx();
    let ny = grid.ny();
    let idx_x = grid.idx_x();
    let idx_y = grid.idx_y();
    let inv_dx = 1.0 / grid.dx;

    let process_row = |j: usize, row: &mut [f64]| {
        let ni_fn = |i: usize| idx_x.next(i);
        let pi_fn = |i: usize| idx_x.prev(i);
        let nj = idx_y.next(j);
        let pj = idx_y.prev(j);

        for i in 0..nx {
            let ni = ni_fn(i);
            let pi = pi_fn(i);

            let vx_right = grid.vx.get(ni, j);
            let f_right = if vx_right >= 0.0 {
                vx_right * grid.s.get(i, j)
            } else {
                vx_right * grid.s.get(ni, j)
            };

            let vx_left = grid.vx.get(i, j);
            let f_left = if vx_left >= 0.0 {
                vx_left * grid.s.get(pi, j)
            } else {
                vx_left * grid.s.get(i, j)
            };

            let vy_top = grid.vy.get(i, nj);
            let f_top =
                if vy_top >= 0.0 { vy_top * grid.s.get(i, j) } else { vy_top * grid.s.get(i, nj) };

            let vy_bot = grid.vy.get(i, j);
            let f_bot =
                if vy_bot >= 0.0 { vy_bot * grid.s.get(i, pj) } else { vy_bot * grid.s.get(i, j) };

            row[i] = (f_right - f_left) * inv_dx + (f_top - f_bot) * inv_dx;
        }
    };

    if nx >= PAR_THRESHOLD {
        div.data_mut().par_chunks_mut(nx).enumerate().for_each(|(j, row)| process_row(j, row));
    } else {
        for j in 0..ny {
            let s = j * nx;
            process_row(j, &mut div.data_mut()[s..s + nx]);
        }
    }
}

/// Compute CFL-limited timestep: dt = cfl_factor * dx / max(|v|).
pub fn compute_cfl_dt(grid: &StaggeredGrid, cfl_factor: f64) -> f64 {
    let mut max_v = 0.0_f64;
    let nx = grid.nx();
    let ny = grid.ny();
    for j in 0..ny {
        for i in 0..nx {
            max_v = max_v.max(grid.vx.get(i, j).abs());
            max_v = max_v.max(grid.vy.get(i, j).abs());
        }
    }
    if max_v < 1e-30 {
        return 1e10; // essentially no velocity — return large dt
    }
    cfl_factor * grid.dx / max_v
}

#[cfg(test)]
mod tests {
    use super::*;

    fn deterministic_rand(state: &mut u64) -> f64 {
        *state =
            state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        (*state >> 33) as f64 / (1u64 << 31) as f64 - 0.5
    }

    #[test]
    fn total_divergence_is_zero() {
        let n = 32;
        let mut grid = StaggeredGrid::new(n, n, 1.0);
        let mut state = 12345u64;

        // Fill with arbitrary S and v
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, 1.0 + deterministic_rand(&mut state) * 0.5);
                grid.vx.set(i, j, deterministic_rand(&mut state));
                grid.vy.set(i, j, deterministic_rand(&mut state));
            }
        }

        let mut div = Field2D::new(n, n);
        compute_divergence_flux(&grid, &mut div);

        let total: f64 = div.data().iter().sum();
        assert!(total.abs() < 1e-10, "Total divergence should be ~0, got {total}");
    }

    #[test]
    fn advection_conserves_mass() {
        let n = 32;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, n, dx);
        let mut state = 42u64;

        // Initial S > 0
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, 1.0 + 0.3 * deterministic_rand(&mut state));
                grid.vx.set(i, j, 0.1 * deterministic_rand(&mut state));
                grid.vy.set(i, j, 0.1 * deterministic_rand(&mut state));
            }
        }

        let initial_mass: f64 = grid.s.data().iter().sum();
        let mut div = Field2D::new(n, n);

        for _ in 0..100 {
            let dt = compute_cfl_dt(&grid, 0.5);
            compute_divergence_flux(&grid, &mut div);
            for j in 0..n {
                for i in 0..n {
                    let s = grid.s.get(i, j) - dt * div.get(i, j);
                    grid.s.set(i, j, s);
                }
            }
        }

        let final_mass: f64 = grid.s.data().iter().sum();
        let rel_err = (final_mass - initial_mass).abs() / initial_mass;
        assert!(rel_err < 1e-12, "Mass not conserved: relative error = {rel_err}");
    }

    #[test]
    fn no_negative_thickness() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, n, dx);

        // Small positive S, moderate velocity
        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, 0.01);
                grid.vx.set(i, j, 0.1);
                grid.vy.set(i, j, 0.05);
            }
        }

        let mut div = Field2D::new(n, n);
        for _ in 0..20 {
            let dt = compute_cfl_dt(&grid, 0.4);
            compute_divergence_flux(&grid, &mut div);
            for j in 0..n {
                for i in 0..n {
                    let s = grid.s.get(i, j) - dt * div.get(i, j);
                    grid.s.set(i, j, s);
                }
            }
        }

        for val in grid.s.data() {
            assert!(*val >= -1e-14, "Negative thickness: {val}");
        }
    }
}
