//! Matrix-free Stokes operator for the thin viscous sheet solver.
//!
//! The operator is discretized on a staggered (MAC) grid with periodic BCs.
//! A negative sign convention is used so the operator is symmetric positive definite,
//! enabling solution by conjugate gradient.

use rayon::prelude::*;

use super::field::Field2D;
use super::grid::StaggeredGrid;
use super::traction::TractionField;

/// Minimum grid dimension for parallel execution. Below this, sequential is faster.
const PAR_THRESHOLD: usize = 64;

/// Lower bound on the absolute value of a Jacobi diagonal before we take
/// its reciprocal. Cells with smaller absolute diagonals are clamped to
/// this floor with their sign preserved, yielding large but finite
/// preconditioner entries. This avoids annihilating information about
/// near-singular cells, which would otherwise force BiCGSTAB to
/// reconstruct their residual contribution indirectly via neighbor
/// coupling.
const JACOBI_DIAG_FLOOR: f64 = 1e-20;

/// Apply the Stokes operator A·v in matrix-free fashion.
///
/// `v` and `out` are packed vectors: first N² entries = vx component, next N² = vy.
/// The operator includes a negative sign so that A is SPD when η > 0.
pub fn apply_stokes(v: &[f64], eta: &Field2D, grid: &StaggeredGrid, out: &mut [f64]) {
    let nx = grid.nx();
    let ny = grid.ny();
    let n2 = nx * ny;
    let idx_x = grid.idx_x();
    let idx_y = grid.idx_y();
    let inv_dx2 = 1.0 / (grid.dx * grid.dx);
    let friction_scaled = grid.basal_friction * inv_dx2;

    let (out_vx, out_vy) = out.split_at_mut(n2);

    // η at cell center (i,j) is eta.get(i,j).
    // We need η interpolated to various staggered positions.
    // η at vertical face midpoint between cells (i,j) and (i, j-1) for the cross term:
    //   averaged from the 4 surrounding cell centers.

    let process_row = |j: usize, row_vx: &mut [f64], row_vy: &mut [f64]| {
        let nj = idx_y.next(j);
        let pj = idx_y.prev(j);

        for i in 0..nx {
            let ni = idx_x.next(i);
            let pi = idx_x.prev(i);
            let li = |ii: usize, jj: usize| -> usize { jj * nx + ii };

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
            let dvy_dx_top = (v[n2 + li(i, nj)] - v[n2 + li(pi, nj)]) * inv_dx2;
            let dvx_dy_bot = (v[li(i, j)] - v[li(i, pj)]) * inv_dx2;
            let dvy_dx_bot = (v[n2 + li(i, j)] - v[n2 + li(pi, j)]) * inv_dx2;

            let term_xy = eta_top * (dvx_dy_top + dvy_dx_top) - eta_bot * (dvx_dy_bot + dvy_dx_bot);

            // Negative sign for SPD
            row_vx[i] = -(term_xx + term_xy);

            // Basal friction: C_b/dx² × S × vx (S interpolated to vx face)
            if friction_scaled > 0.0 {
                let s_face = 0.5 * (grid.s.get(pi, j) + grid.s.get(i, j));
                let s_excess = (s_face - 0.3).max(0.0); // no friction below S=0.3 (thin oceanic)
                row_vx[i] += friction_scaled * s_excess * v[li(i, j)];
            }

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
            let dvx_dy_right_vy = (v[li(ni, j)] - v[li(ni, pj)]) * inv_dx2;
            let dvy_dx_left_vy = (v[n2 + li(i, j)] - v[n2 + li(pi, j)]) * inv_dx2;
            let dvx_dy_left_vy = (v[li(i, j)] - v[li(i, pj)]) * inv_dx2;

            let term_yx = eta_right_vy * (dvy_dx_right_vy + dvx_dy_right_vy)
                - eta_left_vy * (dvy_dx_left_vy + dvx_dy_left_vy);

            row_vy[i] = -(term_yy + term_yx);

            // Basal friction: C_b/dx² × S × vy (S interpolated to vy face)
            if friction_scaled > 0.0 {
                let s_face = 0.5 * (grid.s.get(i, pj) + grid.s.get(i, j));
                let s_excess = (s_face - 0.3).max(0.0); // no friction below S=0.3 (thin oceanic)
                row_vy[i] += friction_scaled * s_excess * v[n2 + li(i, j)];
            }
        }
    };

    if nx >= PAR_THRESHOLD {
        out_vx.par_chunks_mut(nx).zip(out_vy.par_chunks_mut(nx)).enumerate().for_each(
            |(j, (row_vx, row_vy))| {
                process_row(j, row_vx, row_vy);
            },
        );
    } else {
        for j in 0..ny {
            let vx_start = j * nx;
            let (row_vx, row_vy) =
                (&mut out_vx[vx_start..vx_start + nx], &mut out_vy[vx_start..vx_start + nx]);
            // SAFETY: we need non-overlapping mutable borrows of the two halves
            // which split_at_mut already guarantees. The sequential inner slices
            // are also non-overlapping because we index by j*nx.
            process_row(j, row_vx, row_vy);
        }
    }
}

/// Compute the right-hand side: b = -∇(GPE) + T_plates.
///
/// The GPE is density-corrected: `Φ = ρ × (1 - ρ/ρ_mantle) × S²`.
/// When `rho_mantle` is 0 (or rho field is all zeros), falls back to the
/// simple `S²` formulation for backward compatibility.
pub fn compute_rhs(
    grid: &StaggeredGrid,
    plates: &TractionField,
    gravity_factor: f64,
    rho_continental: f64,
    rho_mantle: f64,
    rhs: &mut [f64],
) {
    let nx = grid.nx();
    let ny = grid.ny();
    let n2 = nx * ny;
    let idx_x = grid.idx_x();
    let idx_y = grid.idx_y();
    let inv_dx = 1.0 / grid.dx;
    let use_density = rho_mantle > 0.0;

    let (rhs_vx, rhs_vy) = rhs.split_at_mut(n2);

    // GPE potential at a cell: ρ × (1 - ρ/ρ_m) × S², or just S² if no density
    let gpe = |i: usize, j: usize| -> f64 {
        let s = grid.s.get(i, j);
        if use_density {
            let rho = grid.rho.get(i, j);
            let buoyancy = rho * (1.0 - rho / rho_mantle);
            // Normalize by continental reference so that continental cells
            // have GPE ≈ S² (same as the non-density-corrected formula).
            // Oceanic cells get GPE ≈ 0.60 × S² (40% less spreading pressure).
            let ref_buoyancy = rho_continental * (1.0 - rho_continental / rho_mantle);
            (buoyancy / ref_buoyancy) * s * s
        } else {
            s * s
        }
    };

    let process_row = |j: usize, row_vx: &mut [f64], row_vy: &mut [f64]| {
        let pj = idx_y.prev(j);
        for i in 0..nx {
            let pi = idx_x.prev(i);

            // vx face (i,j): gradient of GPE in x
            let dpdx = -gravity_factor * (gpe(i, j) - gpe(pi, j)) * inv_dx;
            let tx = 0.5 * (plates.tx.get(pi, j) + plates.tx.get(i, j));
            row_vx[i] = dpdx + tx;

            // vy face (i,j): gradient of GPE in y
            let dpdy = -gravity_factor * (gpe(i, j) - gpe(i, pj)) * inv_dx;
            let ty = 0.5 * (plates.ty.get(i, pj) + plates.ty.get(i, j));
            row_vy[i] = dpdy + ty;
        }
    };

    if nx >= PAR_THRESHOLD {
        rhs_vx.par_chunks_mut(nx).zip(rhs_vy.par_chunks_mut(nx)).enumerate().for_each(
            |(j, (row_vx, row_vy))| {
                process_row(j, row_vx, row_vy);
            },
        );
    } else {
        for j in 0..ny {
            let s = j * nx;
            process_row(j, &mut rhs_vx[s..s + nx], &mut rhs_vy[s..s + nx]);
        }
    }
}

/// Compute the Jacobi preconditioner: 1/diag(A) for each DOF.
pub fn compute_jacobi_precond(eta: &Field2D, grid: &StaggeredGrid, precond: &mut [f64]) {
    let nx = grid.nx();
    let ny = grid.ny();
    let n2 = nx * ny;
    let idx_x = grid.idx_x();
    let idx_y = grid.idx_y();
    let inv_dx2 = 1.0 / (grid.dx * grid.dx);
    let friction_scaled = grid.basal_friction * inv_dx2;

    let (pre_vx, pre_vy) = precond.split_at_mut(n2);

    let process_row = |j: usize, row_vx: &mut [f64], row_vy: &mut [f64]| {
        let nj = idx_y.next(j);
        let pj = idx_y.prev(j);
        for i in 0..nx {
            let ni = idx_x.next(i);
            let pi = idx_x.prev(i);

            // vx diagonal: 2(η_right + η_left) + η_top_corner + η_bot_corner
            let eta_right = eta.get(i, j);
            let eta_left = eta.get(pi, j);
            let eta_top =
                0.25 * (eta.get(pi, j) + eta.get(i, j) + eta.get(pi, nj) + eta.get(i, nj));
            let eta_bot =
                0.25 * (eta.get(pi, pj) + eta.get(i, pj) + eta.get(pi, j) + eta.get(i, j));
            let mut diag_vx = inv_dx2 * (2.0 * (eta_right + eta_left) + eta_top + eta_bot);
            if friction_scaled > 0.0 {
                let s_face = 0.5 * (grid.s.get(pi, j) + grid.s.get(i, j));
                let s_excess = (s_face - 0.3).max(0.0); // no friction below S=0.3 (thin oceanic)
                diag_vx += friction_scaled * s_excess;
            }
            // Sign-preserving reciprocal with a floor. Near-singular cells
            // produce a large but finite preconditioner entry instead of
            // the annihilating 0.0 the old ternary yielded — see #50.
            row_vx[i] = diag_vx.signum() / diag_vx.abs().max(JACOBI_DIAG_FLOOR);

            // vy diagonal: η_right_corner + η_left_corner + 2(η_top + η_bot)
            let eta_top_vy = eta.get(i, j);
            let eta_bot_vy = eta.get(i, pj);
            let eta_right_vy =
                0.25 * (eta.get(i, pj) + eta.get(ni, pj) + eta.get(i, j) + eta.get(ni, j));
            let eta_left_vy =
                0.25 * (eta.get(pi, pj) + eta.get(i, pj) + eta.get(pi, j) + eta.get(i, j));
            let mut diag_vy =
                inv_dx2 * (eta_right_vy + eta_left_vy + 2.0 * (eta_top_vy + eta_bot_vy));
            if friction_scaled > 0.0 {
                let s_face = 0.5 * (grid.s.get(i, pj) + grid.s.get(i, j));
                let s_excess = (s_face - 0.3).max(0.0); // no friction below S=0.3 (thin oceanic)
                diag_vy += friction_scaled * s_excess;
            }
            row_vy[i] = diag_vy.signum() / diag_vy.abs().max(JACOBI_DIAG_FLOOR);
        }
    };

    if nx >= PAR_THRESHOLD {
        pre_vx.par_chunks_mut(nx).zip(pre_vy.par_chunks_mut(nx)).enumerate().for_each(
            |(j, (row_vx, row_vy))| {
                process_row(j, row_vx, row_vy);
            },
        );
    } else {
        for j in 0..ny {
            let s = j * nx;
            process_row(j, &mut pre_vx[s..s + nx], &mut pre_vy[s..s + nx]);
        }
    }
}

/// Stencil coefficients for the 5-point Stokes operator, stored per DOF.
///
/// For each DOF, stores [center, left, right, bottom, top] — the absolute
/// values of the matrix entries. The center (diagonal) is positive; the
/// neighbors appear with a negative sign in the operator.
pub struct StencilCoeffs {
    pub vx: Vec<[f64; 5]>,
    pub vy: Vec<[f64; 5]>,
}

impl StencilCoeffs {
    /// Extract stencil coefficients from the current viscosity field.
    pub fn compute(eta: &Field2D, grid: &StaggeredGrid) -> Self {
        let nx = grid.nx();
        let ny = grid.ny();
        let n2 = nx * ny;
        let idx_x = grid.idx_x();
        let idx_y = grid.idx_y();
        let inv_dx2 = 1.0 / (grid.dx * grid.dx);
        let friction_scaled = grid.basal_friction * inv_dx2;

        let mut vx_coeffs = vec![[0.0; 5]; n2];
        let mut vy_coeffs = vec![[0.0; 5]; n2];

        for j in 0..ny {
            let nj = idx_y.next(j);
            let pj = idx_y.prev(j);

            for i in 0..nx {
                let ni = idx_x.next(i);
                let pi = idx_x.prev(i);
                let k = j * nx + i;

                // vx stencil (same η interpolation as apply_stokes)
                let eta_right = eta.get(i, j);
                let eta_left = eta.get(pi, j);
                let eta_top =
                    0.25 * (eta.get(pi, j) + eta.get(i, j) + eta.get(pi, nj) + eta.get(i, nj));
                let eta_bot =
                    0.25 * (eta.get(pi, pj) + eta.get(i, pj) + eta.get(pi, j) + eta.get(i, j));

                let mut diag = inv_dx2 * (2.0 * (eta_right + eta_left) + eta_top + eta_bot);
                if friction_scaled > 0.0 {
                    let s_face = 0.5 * (grid.s.get(pi, j) + grid.s.get(i, j));
                    let s_excess = (s_face - 0.3).max(0.0); // no friction below S=0.3 (thin oceanic)
                    diag += friction_scaled * s_excess;
                }
                let c_left = inv_dx2 * 2.0 * eta_left;
                let c_right = inv_dx2 * 2.0 * eta_right;
                let c_bot = inv_dx2 * eta_bot;
                let c_top = inv_dx2 * eta_top;

                vx_coeffs[k] = [diag, c_left, c_right, c_bot, c_top];

                // vy stencil
                let eta_top_vy = eta.get(i, j);
                let eta_bot_vy = eta.get(i, pj);
                let eta_right_vy =
                    0.25 * (eta.get(i, pj) + eta.get(ni, pj) + eta.get(i, j) + eta.get(ni, j));
                let eta_left_vy =
                    0.25 * (eta.get(pi, pj) + eta.get(i, pj) + eta.get(pi, j) + eta.get(i, j));

                let mut diag_vy =
                    inv_dx2 * (eta_right_vy + eta_left_vy + 2.0 * (eta_top_vy + eta_bot_vy));
                if friction_scaled > 0.0 {
                    let s_face = 0.5 * (grid.s.get(i, pj) + grid.s.get(i, j));
                    let s_excess = (s_face - 0.3).max(0.0); // no friction below S=0.3 (thin oceanic)
                    diag_vy += friction_scaled * s_excess;
                }
                let c_left_vy = inv_dx2 * eta_left_vy;
                let c_right_vy = inv_dx2 * eta_right_vy;
                let c_bot_vy = inv_dx2 * 2.0 * eta_bot_vy;
                let c_top_vy = inv_dx2 * 2.0 * eta_top_vy;

                vy_coeffs[k] = [diag_vy, c_left_vy, c_right_vy, c_bot_vy, c_top_vy];
            }
        }

        Self { vx: vx_coeffs, vy: vy_coeffs }
    }
}

/// Apply SSOR (Symmetric Gauss-Seidel) preconditioner: z = M_SSOR⁻¹ · r.
///
/// Block-diagonal: sweeps vx and vy independently on the 5-point stencil.
/// Forward sweep (row-major order), then backward sweep (reverse order).
/// With periodic BCs, all neighbors use the latest available z values
/// (symmetric Gauss-Seidel, not strict triangular SSOR).
pub fn apply_ssor(
    r: &[f64],
    coeffs: &StencilCoeffs,
    nx: usize,
    ny: usize,
    omega: f64,
    z: &mut [f64],
) {
    let n2 = nx * ny;
    let scale = omega * (2.0 - omega);

    // Process vx block, then vy block independently
    ssor_sweep(&r[..n2], &coeffs.vx, nx, ny, omega, &mut z[..n2]);
    ssor_sweep(&r[n2..], &coeffs.vy, nx, ny, omega, &mut z[n2..]);

    // Scale by ω(2-ω)
    for val in z.iter_mut() {
        *val *= scale;
    }
}

fn ssor_sweep(r: &[f64], coeffs: &[[f64; 5]], nx: usize, ny: usize, omega: f64, z: &mut [f64]) {
    let wrap_x = |i: i32| -> usize { ((i % nx as i32) + nx as i32) as usize % nx };
    let wrap_y = |j: i32| -> usize { ((j % ny as i32) + ny as i32) as usize % ny };

    // Initialize z = 0
    z.iter_mut().for_each(|v| *v = 0.0);

    // Forward sweep (Gauss-Seidel, row by row, left to right)
    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            let [diag, c_left, c_right, c_bot, c_top] = coeffs[k];

            let pi = wrap_x(i as i32 - 1);
            let ni = wrap_x(i as i32 + 1);
            let pj = wrap_y(j as i32 - 1);
            let nj = wrap_y(j as i32 + 1);

            // Use latest z values for all neighbors (symmetric GS on periodic grid)
            let neighbor_sum = c_left * z[j * nx + pi]
                + c_right * z[j * nx + ni]
                + c_bot * z[pj * nx + i]
                + c_top * z[nj * nx + i];

            z[k] = (omega / diag) * (r[k] + neighbor_sum);
        }
    }

    // Backward sweep (reverse order)
    for j in (0..ny).rev() {
        for i in (0..nx).rev() {
            let k = j * nx + i;
            let [diag, c_left, c_right, c_bot, c_top] = coeffs[k];

            let pi = wrap_x(i as i32 - 1);
            let ni = wrap_x(i as i32 + 1);
            let pj = wrap_y(j as i32 - 1);
            let nj = wrap_y(j as i32 + 1);

            let neighbor_sum = c_left * z[j * nx + pi]
                + c_right * z[j * nx + ni]
                + c_bot * z[pj * nx + i]
                + c_top * z[nj * nx + i];

            z[k] = (omega / diag) * (r[k] + neighbor_sum);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn deterministic_rand(state: &mut u64) -> f64 {
        *state =
            state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
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
        let grid = StaggeredGrid::new(n, n, 1.0 / n as f64);
        let eta = Field2D::filled(n, n, 1.0);
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
        let grid = StaggeredGrid::new(n, n, dx);
        let mut eta = Field2D::new(n, n);
        for j in 0..n {
            for i in 0..n {
                let x = (i as f64 + 0.5) * dx;
                let y = (j as f64 + 0.5) * dx;
                eta.set(
                    i,
                    j,
                    1.0 + 0.5
                        * (2.0 * std::f64::consts::PI * x).sin()
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
        assert!(rel_err < 1e-10, "Symmetry violated with variable η: rel_err={rel_err}");
    }

    /// Symmetry must hold on rectangular grids too. With nx != ny, any
    /// bug that confuses idx_x with idx_y in the stencil, or uses the
    /// wrong stride for linear indexing, breaks <Lu, v> = <u, Lv>.
    /// A variable viscosity field is used so coefficient handling is
    /// exercised on both axes.
    #[test]
    fn operator_is_symmetric_rectangular_variable_eta() {
        let nx = 16;
        let ny = 12;
        let dx = 1.0 / nx as f64;
        let grid = StaggeredGrid::new(nx, ny, dx);
        let mut eta = Field2D::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                let x = (i as f64 + 0.5) * dx;
                let y = (j as f64 + 0.5) * dx;
                eta.set(
                    i,
                    j,
                    1.0 + 0.5
                        * (2.0 * std::f64::consts::PI * x).sin()
                        * (2.0 * std::f64::consts::PI * y).cos(),
                );
            }
        }

        let n2 = nx * ny;
        let nn2 = 2 * n2;
        let mut state = 4242u64;
        let mut u = make_random_vec(nn2, &mut state);
        let mut v = make_random_vec(nn2, &mut state);

        // Null-space projection: Stokes with periodic BCs has a rank-2
        // null space (constant mode per velocity component). Project it
        // out of both vectors so the symmetry check isn't dominated by
        // numerical noise in the null direction.
        let project = |w: &mut [f64]| {
            let (wx, wy) = w.split_at_mut(n2);
            let mean_x: f64 = wx.iter().sum::<f64>() / n2 as f64;
            let mean_y: f64 = wy.iter().sum::<f64>() / n2 as f64;
            for a in wx.iter_mut() {
                *a -= mean_x;
            }
            for a in wy.iter_mut() {
                *a -= mean_y;
            }
        };
        project(&mut u);
        project(&mut v);

        let mut au = vec![0.0; nn2];
        let mut av = vec![0.0; nn2];
        apply_stokes(&u, &eta, &grid, &mut au);
        apply_stokes(&v, &eta, &grid, &mut av);

        let u_av = dot(&u, &av);
        let au_v = dot(&au, &v);
        let rel_err = (u_av - au_v).abs() / u_av.abs().max(1e-14);
        assert!(
            rel_err < 1e-10,
            "Symmetry violated on {nx}×{ny} grid: <u,Lv>={u_av}, <Lu,v>={au_v}, rel_err={rel_err}"
        );
    }

    #[test]
    fn operator_is_positive_definite() {
        let n = 16;
        let grid = StaggeredGrid::new(n, n, 1.0 / n as f64);
        let eta = Field2D::filled(n, n, 1.0);
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
        let grid = StaggeredGrid::new(n, n, dx);
        let eta = Field2D::filled(n, n, 1.0);
        let nn2 = 2 * n * n;

        let k = 2.0 * std::f64::consts::PI;

        let mut v = vec![0.0; nn2];
        for j in 0..n {
            for i in 0..n {
                let x = i as f64 * dx;
                v[j * n + i] = (k * x).sin();
            }
        }

        let mut av = vec![0.0; nn2];
        apply_stokes(&v, &eta, &grid, &mut av);

        // -∂/∂x(2η·∂vx/∂x) with η=1 has discrete eigenvalue 2·(2/dx²)·2sin²(k·dx/2)
        let expected_eigenvalue = 2.0 * (4.0 / (dx * dx)) * (k * dx / 2.0).sin().powi(2);

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

    #[test]
    fn parallel_apply_stokes_matches_sequential() {
        let n = 64; // Above PAR_THRESHOLD to test parallel path
        let dx = 1.0 / n as f64;
        let grid = StaggeredGrid::new(n, n, dx);
        let eta = Field2D::filled(n, n, 1.0);
        let nn2 = 2 * n * n;

        let mut state = 42u64;
        let v: Vec<f64> = (0..nn2).map(|_| deterministic_rand(&mut state)).collect();

        // Sequential reference
        let mut out_seq = vec![0.0; nn2];
        apply_stokes_sequential(&v, &eta, &grid, &mut out_seq);

        // Parallel (default apply_stokes with n >= PAR_THRESHOLD)
        let mut out_par = vec![0.0; nn2];
        apply_stokes(&v, &eta, &grid, &mut out_par);

        for i in 0..nn2 {
            let err = (out_seq[i] - out_par[i]).abs();
            assert!(
                err < 1e-14,
                "Mismatch at {i}: seq={}, par={}, diff={err}",
                out_seq[i],
                out_par[i]
            );
        }
    }

    /// Sequential reference implementation for regression testing.
    fn apply_stokes_sequential(v: &[f64], eta: &Field2D, grid: &StaggeredGrid, out: &mut [f64]) {
        let nx = grid.nx();
        let ny = grid.ny();
        let n2 = nx * ny;
        let idx_x = grid.idx_x();
        let idx_y = grid.idx_y();
        let inv_dx2 = 1.0 / (grid.dx * grid.dx);
        let li = |i: usize, j: usize| -> usize { j * nx + i };

        for j in 0..ny {
            for i in 0..nx {
                let ni = idx_x.next(i);
                let pi = idx_x.prev(i);
                let nj = idx_y.next(j);
                let pj = idx_y.prev(j);

                let eta_right = eta.get(i, j);
                let eta_left = eta.get(pi, j);
                let dvx_dx_right = (v[li(ni, j)] - v[li(i, j)]) * inv_dx2;
                let dvx_dx_left = (v[li(i, j)] - v[li(pi, j)]) * inv_dx2;
                let term_xx = 2.0 * eta_right * dvx_dx_right - 2.0 * eta_left * dvx_dx_left;

                let eta_top =
                    0.25 * (eta.get(pi, j) + eta.get(i, j) + eta.get(pi, nj) + eta.get(i, nj));
                let eta_bot =
                    0.25 * (eta.get(pi, pj) + eta.get(i, pj) + eta.get(pi, j) + eta.get(i, j));
                let dvx_dy_top = (v[li(i, nj)] - v[li(i, j)]) * inv_dx2;
                let dvy_dx_top = (v[n2 + li(i, nj)] - v[n2 + li(pi, nj)]) * inv_dx2;
                let dvx_dy_bot = (v[li(i, j)] - v[li(i, pj)]) * inv_dx2;
                let dvy_dx_bot = (v[n2 + li(i, j)] - v[n2 + li(pi, j)]) * inv_dx2;
                let term_xy =
                    eta_top * (dvx_dy_top + dvy_dx_top) - eta_bot * (dvx_dy_bot + dvy_dx_bot);
                out[li(i, j)] = -(term_xx + term_xy);

                let eta_top_vy = eta.get(i, j);
                let eta_bot_vy = eta.get(i, pj);
                let dvy_dy_top = (v[n2 + li(i, nj)] - v[n2 + li(i, j)]) * inv_dx2;
                let dvy_dy_bot = (v[n2 + li(i, j)] - v[n2 + li(i, pj)]) * inv_dx2;
                let term_yy = 2.0 * eta_top_vy * dvy_dy_top - 2.0 * eta_bot_vy * dvy_dy_bot;

                let eta_right_vy =
                    0.25 * (eta.get(i, pj) + eta.get(ni, pj) + eta.get(i, j) + eta.get(ni, j));
                let eta_left_vy =
                    0.25 * (eta.get(pi, pj) + eta.get(i, pj) + eta.get(pi, j) + eta.get(i, j));
                let dvy_dx_right_vy = (v[n2 + li(ni, j)] - v[n2 + li(i, j)]) * inv_dx2;
                let dvx_dy_right_vy = (v[li(ni, j)] - v[li(ni, pj)]) * inv_dx2;
                let dvy_dx_left_vy = (v[n2 + li(i, j)] - v[n2 + li(pi, j)]) * inv_dx2;
                let dvx_dy_left_vy = (v[li(i, j)] - v[li(i, pj)]) * inv_dx2;
                let term_yx = eta_right_vy * (dvy_dx_right_vy + dvx_dy_right_vy)
                    - eta_left_vy * (dvy_dx_left_vy + dvx_dy_left_vy);
                out[n2 + li(i, j)] = -(term_yy + term_yx);
            }
        }
    }

    #[test]
    fn operator_with_friction_is_symmetric() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, n, dx);
        grid.basal_friction = 2.0;

        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, 0.5 + 0.3 * ((i + j) as f64 / n as f64));
            }
        }

        let eta = Field2D::filled(n, n, 1.0);
        let n_dof = 2 * n * n;

        let mut u = vec![0.0; n_dof];
        let mut w = vec![0.0; n_dof];
        let mut au = vec![0.0; n_dof];
        let mut aw = vec![0.0; n_dof];

        use rand::Rng;
        let mut rng = rand::thread_rng();
        for i in 0..n_dof {
            u[i] = rng.r#gen::<f64>() - 0.5;
            w[i] = rng.r#gen::<f64>() - 0.5;
        }

        apply_stokes(&u, &eta, &grid, &mut au);
        apply_stokes(&w, &eta, &grid, &mut aw);

        let u_dot_aw: f64 = u.iter().zip(aw.iter()).map(|(a, b)| a * b).sum();
        let au_dot_w: f64 = au.iter().zip(w.iter()).map(|(a, b)| a * b).sum();

        let rel_diff = (u_dot_aw - au_dot_w).abs() / (u_dot_aw.abs() + au_dot_w.abs()).max(1e-30);
        assert!(rel_diff < 1e-12, "Operator with friction should be symmetric: {rel_diff}");
    }

    #[test]
    fn zero_friction_is_backward_compatible() {
        let n = 8;
        let dx = 1.0 / n as f64;
        let mut grid = StaggeredGrid::new(n, n, dx);
        grid.basal_friction = 0.0;

        for j in 0..n {
            for i in 0..n {
                grid.s.set(i, j, 1.0);
            }
        }

        let eta = Field2D::filled(n, n, 1.0);
        let n_dof = 2 * n * n;

        let mut v = vec![0.0; n_dof];
        for i in 0..n_dof {
            v[i] = (i as f64) * 0.01;
        }

        let mut out = vec![0.0; n_dof];
        apply_stokes(&v, &eta, &grid, &mut out);

        for &val in &out {
            assert!(val.is_finite(), "Output should be finite");
        }
    }

    #[test]
    fn jacobi_precond_handles_near_singular_cells() {
        // Uniform η at 1e-25 simulates an early-continuation or
        // pathological cratonic-core state where the whole vx stencil
        // collapses. With dx = 0.1 and the 2·(η_R + η_L) + η_T + η_B
        // stencil that gives a raw diagonal of ~6e-23, well below the
        // 1e-20 JACOBI_DIAG_FLOOR. Under the old ternary (abs > 1e-30)
        // the reciprocal would survive unclamped and blow up to ~1.6e+22,
        // but a slightly smaller η (1e-40) would have fallen below 1e-30
        // and been annihilated to 0.0. The new floored reciprocal always
        // yields a finite, sign-preserving, bounded-magnitude entry.
        let n = 8;
        let dx = 0.1;
        let grid = StaggeredGrid::new(n, n, dx);
        let mut eta = Field2D::new(n, n);
        for k in 0..n * n {
            eta.data_mut()[k] = 1e-40;
        }

        let mut precond = vec![0.0; 2 * n * n];
        compute_jacobi_precond(&eta, &grid, &mut precond);

        for (i, &val) in precond.iter().enumerate() {
            assert!(val.is_finite(), "non-finite preconditioner entry at index {i}: {val}");
        }

        // Every entry should saturate to 1/JACOBI_DIAG_FLOOR = 1e20
        // instead of the old code's silent 0.0.
        for (i, &val) in precond.iter().enumerate() {
            assert!(
                val.abs() > 1e10,
                "expected near-singular cell to produce large preconditioner entry at {i}, got {val}"
            );
        }
    }

    #[test]
    fn jacobi_precond_unchanged_on_normal_cells() {
        // On a uniform well-conditioned η field every Jacobi entry should
        // be a finite positive number of reasonable magnitude — the floor
        // change from #50 must not perturb the normal regime.
        let n = 4;
        let dx = 0.25;
        let grid = StaggeredGrid::new(n, n, dx);
        let mut eta = Field2D::new(n, n);
        for k in 0..n * n {
            eta.data_mut()[k] = 1.0;
        }

        let mut precond = vec![0.0; 2 * n * n];
        compute_jacobi_precond(&eta, &grid, &mut precond);

        for &val in precond.iter() {
            assert!(val.is_finite());
            assert!(val > 0.0 && val < 10.0, "expected normal preconditioner magnitude, got {val}");
        }
    }
}
