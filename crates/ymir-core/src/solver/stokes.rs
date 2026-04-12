//! Matrix-free Stokes operator for the thin viscous sheet solver.
//!
//! The operator is discretized on a staggered (MAC) grid with periodic BCs.
//! A negative sign convention is used so the operator is symmetric positive definite,
//! enabling solution by conjugate gradient.

use super::field::Field2D;
use super::grid::StaggeredGrid;
use super::plates::PlateField;

/// Apply the Stokes operator A·v in matrix-free fashion.
///
/// `v` and `out` are packed vectors: first N² entries = vx component, next N² = vy.
/// The operator includes a negative sign so that A is SPD when η > 0.
pub fn apply_stokes(v: &[f64], eta: &Field2D, grid: &StaggeredGrid, out: &mut [f64]) {
    let n = grid.n;
    let n2 = n * n;
    let idx = &grid.idx;
    let inv_dx2 = 1.0 / (grid.dx * grid.dx);

    // Helper: index into flat arrays
    let li = |i: usize, j: usize| -> usize { j * n + i };

    // η at cell center (i,j) is eta.get(i,j).
    // We need η interpolated to various staggered positions.
    // η at vertical face midpoint between cells (i,j) and (i, j-1) for the cross term:
    //   averaged from the 4 surrounding cell centers.

    for j in 0..n {
        for i in 0..n {
            let ni = idx.next(i);
            let pi = idx.prev(i);
            let nj = idx.next(j);
            let pj = idx.prev(j);

            // --- vx component at face (i, j) ---
            // ∂/∂x(2η ∂vx/∂x): uses η at cell centers (i,j) and (prev(i),j)
            let eta_right = eta.get(i, j); // cell to the right of face (i,j)
            let eta_left = eta.get(pi, j); // cell to the left

            let dvx_dx_right = (v[li(ni, j)] - v[li(i, j)]) * inv_dx2;
            let dvx_dx_left = (v[li(i, j)] - v[li(pi, j)]) * inv_dx2;
            let term_xx = 2.0 * eta_right * dvx_dx_right - 2.0 * eta_left * dvx_dx_left;

            // ∂/∂y(η(∂vx/∂y + ∂vy/∂x)): η at corners
            // Corner (i, next(j)) — average of 4 cell centers
            let eta_top =
                0.25 * (eta.get(pi, j) + eta.get(i, j) + eta.get(pi, nj) + eta.get(i, nj));
            // Corner (i, j) — average of 4 cell centers
            let eta_bot =
                0.25 * (eta.get(pi, pj) + eta.get(i, pj) + eta.get(pi, j) + eta.get(i, j));

            let dvx_dy_top = (v[li(i, nj)] - v[li(i, j)]) * inv_dx2;
            let dvy_dx_top =
                (v[n2 + li(i, nj)] - v[n2 + li(pi, nj)]) * inv_dx2;
            let dvx_dy_bot = (v[li(i, j)] - v[li(i, pj)]) * inv_dx2;
            let dvy_dx_bot =
                (v[n2 + li(i, j)] - v[n2 + li(pi, j)]) * inv_dx2;

            let term_xy =
                eta_top * (dvx_dy_top + dvy_dx_top) - eta_bot * (dvx_dy_bot + dvy_dx_bot);

            // Negative sign for SPD
            out[li(i, j)] = -(term_xx + term_xy);

            // --- vy component at face (i, j) ---
            // ∂/∂y(2η ∂vy/∂y): uses η at cell centers (i,j) and (i, prev(j))
            let eta_top_vy = eta.get(i, j);
            let eta_bot_vy = eta.get(i, pj);

            let dvy_dy_top = (v[n2 + li(i, nj)] - v[n2 + li(i, j)]) * inv_dx2;
            let dvy_dy_bot = (v[n2 + li(i, j)] - v[n2 + li(i, pj)]) * inv_dx2;
            let term_yy = 2.0 * eta_top_vy * dvy_dy_top - 2.0 * eta_bot_vy * dvy_dy_bot;

            // ∂/∂x(η(∂vx/∂y + ∂vy/∂x)): η at corners
            // Corner (next(i), j) — average of 4 cell centers
            let eta_right_vy =
                0.25 * (eta.get(i, pj) + eta.get(ni, pj) + eta.get(i, j) + eta.get(ni, j));
            // Corner (i, j) — reuse eta_bot from above (same corner position)
            let eta_left_vy =
                0.25 * (eta.get(pi, pj) + eta.get(i, pj) + eta.get(pi, j) + eta.get(i, j));

            let dvy_dx_right_vy = (v[n2 + li(ni, j)] - v[n2 + li(i, j)]) * inv_dx2;
            let dvx_dy_right_vy =
                (v[li(ni, j)] - v[li(ni, pj)]) * inv_dx2;
            let dvy_dx_left_vy = (v[n2 + li(i, j)] - v[n2 + li(pi, j)]) * inv_dx2;
            let dvx_dy_left_vy =
                (v[li(i, j)] - v[li(i, pj)]) * inv_dx2;

            let term_yx = eta_right_vy * (dvy_dx_right_vy + dvx_dy_right_vy)
                - eta_left_vy * (dvy_dx_left_vy + dvx_dy_left_vy);

            out[n2 + li(i, j)] = -(term_yy + term_yx);
        }
    }
}

/// Compute the right-hand side: b = ∇(gravity_factor · S²) + T_plates.
///
/// The gradient of S² is evaluated at the staggered face positions.
pub fn compute_rhs(
    grid: &StaggeredGrid,
    plates: &PlateField,
    gravity_factor: f64,
    rhs: &mut [f64],
) {
    let n = grid.n;
    let n2 = n * n;
    let idx = &grid.idx;
    let inv_dx = 1.0 / grid.dx;
    let li = |i: usize, j: usize| -> usize { j * n + i };

    for j in 0..n {
        for i in 0..n {
            let pi = idx.prev(i);
            let pj = idx.prev(j);

            // vx face (i,j): gradient of S² in x between cells (i,j) and (prev(i),j)
            let s2_right = grid.s.get(i, j) * grid.s.get(i, j);
            let s2_left = grid.s.get(pi, j) * grid.s.get(pi, j);
            // Negative sign: A has -∇·τ, so rhs gets -∇(gS²) to balance
            let dpdx = -gravity_factor * (s2_right - s2_left) * inv_dx;

            // Traction: average of the two adjacent cells
            let tx = 0.5 * (plates.tx.get(pi, j) + plates.tx.get(i, j));

            rhs[li(i, j)] = dpdx + tx;

            // vy face (i,j): gradient of S² in y between cells (i,j) and (i, prev(j))
            let s2_top = grid.s.get(i, j) * grid.s.get(i, j);
            let s2_bot = grid.s.get(i, pj) * grid.s.get(i, pj);
            let dpdy = -gravity_factor * (s2_top - s2_bot) * inv_dx;

            let ty = 0.5 * (plates.ty.get(i, pj) + plates.ty.get(i, j));

            rhs[n2 + li(i, j)] = dpdy + ty;
        }
    }
}

/// Compute the Jacobi preconditioner: 1/diag(A) for each DOF.
pub fn compute_jacobi_precond(eta: &Field2D, grid: &StaggeredGrid, precond: &mut [f64]) {
    let n = grid.n;
    let n2 = n * n;
    let idx = &grid.idx;
    let inv_dx2 = 1.0 / (grid.dx * grid.dx);
    let li = |i: usize, j: usize| -> usize { j * n + i };

    for j in 0..n {
        for i in 0..n {
            let ni = idx.next(i);
            let pi = idx.prev(i);
            let nj = idx.next(j);
            let pj = idx.prev(j);

            // vx diagonal: 2(η_right + η_left) + η_top_corner + η_bot_corner
            let eta_right = eta.get(i, j);
            let eta_left = eta.get(pi, j);
            let eta_top =
                0.25 * (eta.get(pi, j) + eta.get(i, j) + eta.get(pi, nj) + eta.get(i, nj));
            let eta_bot =
                0.25 * (eta.get(pi, pj) + eta.get(i, pj) + eta.get(pi, j) + eta.get(i, j));
            let diag_vx = inv_dx2 * (2.0 * (eta_right + eta_left) + eta_top + eta_bot);
            precond[li(i, j)] = if diag_vx.abs() > 1e-30 {
                1.0 / diag_vx
            } else {
                0.0
            };

            // vy diagonal: η_right_corner + η_left_corner + 2(η_top + η_bot)
            let eta_top_vy = eta.get(i, j);
            let eta_bot_vy = eta.get(i, pj);
            let eta_right_vy =
                0.25 * (eta.get(i, pj) + eta.get(ni, pj) + eta.get(i, j) + eta.get(ni, j));
            let eta_left_vy =
                0.25 * (eta.get(pi, pj) + eta.get(i, pj) + eta.get(pi, j) + eta.get(i, j));
            let diag_vy = inv_dx2 * (eta_right_vy + eta_left_vy + 2.0 * (eta_top_vy + eta_bot_vy));
            precond[n2 + li(i, j)] = if diag_vy.abs() > 1e-30 {
                1.0 / diag_vy
            } else {
                0.0
            };
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn deterministic_rand(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (*state >> 33) as f64 / (1u64 << 31) as f64 - 0.5
    }

    fn make_random_vec(len: usize, state: &mut u64) -> Vec<f64> {
        (0..len).map(|_| deterministic_rand(state)).collect()
    }

    fn dot(a: &[f64], b: &[f64]) -> f64 {
        a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
    }

    #[test]
    fn operator_is_symmetric_uniform_eta() {
        let n = 16;
        let grid = StaggeredGrid::new(n, 1.0 / n as f64);
        let eta = Field2D::filled(n, 1.0);
        let nn2 = 2 * n * n;

        let mut state = 99u64;
        let u = make_random_vec(nn2, &mut state);
        let v = make_random_vec(nn2, &mut state);

        let mut au = vec![0.0; nn2];
        let mut av = vec![0.0; nn2];
        apply_stokes(&u, &eta, &grid, &mut au);
        apply_stokes(&v, &eta, &grid, &mut av);

        let u_av = dot(&u, &av);
        let au_v = dot(&au, &v);
        let rel_err = (u_av - au_v).abs() / u_av.abs().max(1e-14);
        assert!(
            rel_err < 1e-10,
            "Symmetry violated: <u,Av>={u_av}, <Au,v>={au_v}, rel_err={rel_err}"
        );
    }

    #[test]
    fn operator_is_symmetric_variable_eta() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let grid = StaggeredGrid::new(n, dx);
        let mut eta = Field2D::new(n);
        for j in 0..n {
            for i in 0..n {
                let x = (i as f64 + 0.5) * dx;
                let y = (j as f64 + 0.5) * dx;
                eta.set(
                    i,
                    j,
                    1.0 + 0.5 * (2.0 * std::f64::consts::PI * x).sin()
                        * (2.0 * std::f64::consts::PI * y).cos(),
                );
            }
        }

        let nn2 = 2 * n * n;
        let mut state = 777u64;
        let u = make_random_vec(nn2, &mut state);
        let v = make_random_vec(nn2, &mut state);

        let mut au = vec![0.0; nn2];
        let mut av = vec![0.0; nn2];
        apply_stokes(&u, &eta, &grid, &mut au);
        apply_stokes(&v, &eta, &grid, &mut av);

        let u_av = dot(&u, &av);
        let au_v = dot(&au, &v);
        let rel_err = (u_av - au_v).abs() / u_av.abs().max(1e-14);
        assert!(
            rel_err < 1e-10,
            "Symmetry violated with variable η: rel_err={rel_err}"
        );
    }

    #[test]
    fn operator_is_positive_definite() {
        let n = 16;
        let grid = StaggeredGrid::new(n, 1.0 / n as f64);
        let eta = Field2D::filled(n, 1.0);
        let nn2 = 2 * n * n;

        let mut state = 1234u64;
        for _ in 0..10 {
            let v = make_random_vec(nn2, &mut state);
            let mut av = vec![0.0; nn2];
            apply_stokes(&v, &eta, &grid, &mut av);
            let vav = dot(&v, &av);
            assert!(vav > 0.0, "<v, Av> = {vav} should be > 0");
        }
    }

    #[test]
    fn eigenvalue_of_sine_mode() {
        let n = 32;
        let dx = 1.0 / n as f64;
        let grid = StaggeredGrid::new(n, dx);
        let eta = Field2D::filled(n, 1.0);
        let nn2 = 2 * n * n;

        let k = 2.0 * std::f64::consts::PI; // wavenumber for mode fitting in [0,1]

        let mut v = vec![0.0; nn2];
        // vx = sin(k * x_face), vy = 0
        for j in 0..n {
            for i in 0..n {
                let x = i as f64 * dx; // vx face position
                v[j * n + i] = (k * x).sin();
            }
        }

        let mut av = vec![0.0; nn2];
        apply_stokes(&v, &eta, &grid, &mut av);

        // Expected eigenvalue for -∂/∂x(2η ∂vx/∂x) with η=1: 2k²
        // But discrete FD has eigenvalue 2 * (2/dx²) * sin²(k*dx/2)
        // For the ∂/∂y terms, vx doesn't vary in y so that contribution ≈ 0.
        // -∂/∂x(2η·∂vx/∂x) with η=1 has discrete eigenvalue 2·(2/dx²)·2sin²(k·dx/2)
        let expected_eigenvalue = 2.0 * (4.0 / (dx * dx)) * (k * dx / 2.0).sin().powi(2);

        // Check that Av ≈ λ * v for the vx component
        let mut max_err = 0.0_f64;
        let mut max_val = 0.0_f64;
        for j in 0..n {
            for i in 0..n {
                let idx = j * n + i;
                let expected = expected_eigenvalue * v[idx];
                let err = (av[idx] - expected).abs();
                max_err = max_err.max(err);
                max_val = max_val.max(v[idx].abs());
            }
        }

        let rel_err = max_err / (expected_eigenvalue * max_val);
        assert!(
            rel_err < 0.05,
            "Sine eigenvalue test: rel_err = {rel_err}, expected λ = {expected_eigenvalue}"
        );
    }
}
