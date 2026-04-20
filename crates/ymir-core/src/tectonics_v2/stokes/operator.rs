//! Discrete Stokes operator on a MAC (staggered) grid with periodic BCs.
//!
//! Layout (same convention as the legacy `tectonics/solver/grid.rs`):
//! - p, η, S at cell centers (i, j) — physical position `((i+0.5)dx, (j+0.5)dy)`.
//! - vx at left vertical face of cell (i, j) — position `(i dx, (j+0.5)dy)`.
//! - vy at bottom horizontal face of cell (i, j) — position `((i+0.5)dx, j dy)`.
//! - ε̇_xy and σ_xy live at corners (i·dx, j·dy); η there is computed by
//!   **harmonic averaging** of the four surrounding cell centers, so
//!   Step 1's variable-η rheology is a data change, not an assembly change.
//!
//! Indexing uses [`PeriodicIndex`] for wraparound. All flat arrays are
//! row-major with stride `nx`: `arr[j * nx + i]`.

use super::super::field::{Field2D, PeriodicIndex};

/// Geometry and viscosity field needed by the momentum and divergence
/// operators. Borrows η non-mutably — callers update η out of band.
pub struct StokesGrid {
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub dy: f64,
    pub idx_x: PeriodicIndex,
    pub idx_y: PeriodicIndex,
}

impl StokesGrid {
    pub fn new(nx: usize, ny: usize, dx: f64, dy: f64) -> Self {
        Self {
            nx,
            ny,
            dx,
            dy,
            idx_x: PeriodicIndex::new(nx),
            idx_y: PeriodicIndex::new(ny),
        }
    }

    #[inline]
    pub fn n_cells(&self) -> usize {
        self.nx * self.ny
    }
}

/// Harmonic average of four positive quantities.
///
/// Reduces to the arithmetic value when all four inputs are equal
/// (constant-η case at Step 0). Guards against division by zero if any
/// argument is non-positive by falling back to the arithmetic mean.
#[inline]
fn harmonic4(a: f64, b: f64, c: f64, d: f64) -> f64 {
    if a > 0.0 && b > 0.0 && c > 0.0 && d > 0.0 {
        4.0 / (1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d)
    } else {
        0.25 * (a + b + c + d)
    }
}

/// η at the four corners surrounding cell face/center, harmonic-averaged.
#[inline]
fn eta_corner(eta: &Field2D, im: usize, i: usize, jm: usize, j: usize) -> f64 {
    harmonic4(eta.get(im, jm), eta.get(i, jm), eta.get(im, j), eta.get(i, j))
}

/// Apply the momentum operator `A v = -∇·(2η ε̇(v))` on the MAC grid.
///
/// Inputs are flat arrays of length `nx*ny` each; outputs are the same
/// shape. Uses `(grid.dx, grid.dy)` spacing and the supplied η field.
/// The sign convention makes A symmetric positive semi-definite (SPD
/// modulo the 2-D null space of constant velocities per component).
pub fn apply_momentum(
    grid: &StokesGrid,
    eta: &Field2D,
    vx: &[f64],
    vy: &[f64],
    out_vx: &mut [f64],
    out_vy: &mut [f64],
) {
    let nx = grid.nx;
    let ny = grid.ny;
    let dx = grid.dx;
    let dy = grid.dy;
    let inv_dx = 1.0 / dx;
    let inv_dy = 1.0 / dy;
    debug_assert_eq!(vx.len(), nx * ny);
    debug_assert_eq!(vy.len(), nx * ny);

    for j in 0..ny {
        let jp = grid.idx_y.next(j);
        let jm = grid.idx_y.prev(j);
        for i in 0..nx {
            let ip = grid.idx_x.next(i);
            let im = grid.idx_x.prev(i);
            let lin = |ii: usize, jj: usize| jj * nx + ii;

            // ---------- x-momentum at vx(i, j) ----------
            //   A_x v = -∂x(2 η ∂x vx) - ∂y(η (∂y vx + ∂x vy))
            // Normal-stress flux at cell centers (i, j) and (im, j):
            let eta_cc_right = eta.get(i, j);
            let eta_cc_left = eta.get(im, j);
            let dvx_dx_right = (vx[lin(ip, j)] - vx[lin(i, j)]) * inv_dx;
            let dvx_dx_left = (vx[lin(i, j)] - vx[lin(im, j)]) * inv_dx;
            let sigma_xx_right = 2.0 * eta_cc_right * dvx_dx_right;
            let sigma_xx_left = 2.0 * eta_cc_left * dvx_dx_left;
            let d_sigma_xx_dx = (sigma_xx_right - sigma_xx_left) * inv_dx;

            // Shear-stress flux at corners (i, j) and (i, jp):
            let eta_corner_top = eta_corner(eta, im, i, j, jp);
            let eta_corner_bot = eta_corner(eta, im, i, jm, j);
            let dvx_dy_top = (vx[lin(i, jp)] - vx[lin(i, j)]) * inv_dy;
            let dvx_dy_bot = (vx[lin(i, j)] - vx[lin(i, jm)]) * inv_dy;
            let dvy_dx_top = (vy[lin(i, jp)] - vy[lin(im, jp)]) * inv_dx;
            let dvy_dx_bot = (vy[lin(i, j)] - vy[lin(im, j)]) * inv_dx;
            let sigma_xy_top = eta_corner_top * (dvx_dy_top + dvy_dx_top);
            let sigma_xy_bot = eta_corner_bot * (dvx_dy_bot + dvy_dx_bot);
            let d_sigma_xy_dy = (sigma_xy_top - sigma_xy_bot) * inv_dy;

            out_vx[lin(i, j)] = -(d_sigma_xx_dx + d_sigma_xy_dy);

            // ---------- y-momentum at vy(i, j) ----------
            //   A_y v = -∂x(η (∂y vx + ∂x vy)) - ∂y(2 η ∂y vy)
            // Normal-stress flux at cell centers (i, j) and (i, jm):
            let eta_cc_top = eta.get(i, j);
            let eta_cc_bot = eta.get(i, jm);
            let dvy_dy_top = (vy[lin(i, jp)] - vy[lin(i, j)]) * inv_dy;
            let dvy_dy_bot = (vy[lin(i, j)] - vy[lin(i, jm)]) * inv_dy;
            let sigma_yy_top = 2.0 * eta_cc_top * dvy_dy_top;
            let sigma_yy_bot = 2.0 * eta_cc_bot * dvy_dy_bot;
            let d_sigma_yy_dy = (sigma_yy_top - sigma_yy_bot) * inv_dy;

            // Shear-stress flux at corners (ip, j) and (i, j):
            let eta_corner_right = eta_corner(eta, i, ip, jm, j);
            let eta_corner_left = eta_corner(eta, im, i, jm, j);
            let dvx_dy_right = (vx[lin(ip, j)] - vx[lin(ip, jm)]) * inv_dy;
            let dvx_dy_left = (vx[lin(i, j)] - vx[lin(i, jm)]) * inv_dy;
            let dvy_dx_right = (vy[lin(ip, j)] - vy[lin(i, j)]) * inv_dx;
            let dvy_dx_left = (vy[lin(i, j)] - vy[lin(im, j)]) * inv_dx;
            let sigma_xy_right = eta_corner_right * (dvx_dy_right + dvy_dx_right);
            let sigma_xy_left = eta_corner_left * (dvx_dy_left + dvy_dx_left);
            let d_sigma_xy_dx = (sigma_xy_right - sigma_xy_left) * inv_dx;

            out_vy[lin(i, j)] = -(d_sigma_xy_dx + d_sigma_yy_dy);
        }
    }
}

/// Diagonal of the momentum operator for Jacobi preconditioning.
///
/// For the vx(i, j) DOF the diagonal coefficient reads
/// `2(η[i-1,j]+η[i,j])/dx² + (η_c[i,j]+η_c[i,jp])/dy²` (and symmetrically
/// for vy), where `η_c` are harmonic-corner viscosities. Reduces to
/// `4η/dx² + 2η/dy²` (= 6 for η=dx=dy=1) when η is constant — the
/// expected value for the discrete variable-coefficient Stokes momentum
/// operator on this MAC layout.
pub fn momentum_diagonal(
    grid: &StokesGrid,
    eta: &Field2D,
    diag_vx: &mut [f64],
    diag_vy: &mut [f64],
) {
    let nx = grid.nx;
    let ny = grid.ny;
    let dx = grid.dx;
    let dy = grid.dy;
    let inv_dx2 = 1.0 / (dx * dx);
    let inv_dy2 = 1.0 / (dy * dy);

    for j in 0..ny {
        let jp = grid.idx_y.next(j);
        let jm = grid.idx_y.prev(j);
        for i in 0..nx {
            let ip = grid.idx_x.next(i);
            let im = grid.idx_x.prev(i);
            let lin = |ii: usize, jj: usize| jj * nx + ii;

            let eta_right_cc = eta.get(i, j);
            let eta_left_cc = eta.get(im, j);
            let eta_c_top = eta_corner(eta, im, i, j, jp);
            let eta_c_bot = eta_corner(eta, im, i, jm, j);
            diag_vx[lin(i, j)] =
                2.0 * (eta_right_cc + eta_left_cc) * inv_dx2 + (eta_c_top + eta_c_bot) * inv_dy2;

            let eta_top_cc = eta.get(i, j);
            let eta_bot_cc = eta.get(i, jm);
            let eta_c_right = eta_corner(eta, i, ip, jm, j);
            let eta_c_left = eta_corner(eta, im, i, jm, j);
            diag_vy[lin(i, j)] =
                (eta_c_right + eta_c_left) * inv_dx2 + 2.0 * (eta_top_cc + eta_bot_cc) * inv_dy2;
        }
    }
}

/// Discrete divergence `(B v)[i,j] = (vx[ip,j] - vx[i,j])/dx + (vy[i,jp] - vy[i,j])/dy`.
pub fn apply_divergence(
    grid: &StokesGrid,
    vx: &[f64],
    vy: &[f64],
    out_p: &mut [f64],
) {
    let nx = grid.nx;
    let ny = grid.ny;
    let inv_dx = 1.0 / grid.dx;
    let inv_dy = 1.0 / grid.dy;

    for j in 0..ny {
        let jp = grid.idx_y.next(j);
        for i in 0..nx {
            let ip = grid.idx_x.next(i);
            let lin = |ii: usize, jj: usize| jj * nx + ii;
            out_p[lin(i, j)] = (vx[lin(ip, j)] - vx[lin(i, j)]) * inv_dx
                + (vy[lin(i, jp)] - vy[lin(i, j)]) * inv_dy;
        }
    }
}

/// Discrete adjoint of `apply_divergence`: `B^T p` lives at velocity faces.
///
/// Derived from `⟨B v, p⟩ = ⟨v, B^T p⟩` under the dot product that
/// weights every DOF uniformly (no volume factors, per standard MAC
/// Stokes convention). The adjoint acts as a discrete negative gradient.
pub fn apply_divergence_transpose(
    grid: &StokesGrid,
    p: &[f64],
    out_vx: &mut [f64],
    out_vy: &mut [f64],
) {
    let nx = grid.nx;
    let ny = grid.ny;
    let inv_dx = 1.0 / grid.dx;
    let inv_dy = 1.0 / grid.dy;

    for j in 0..ny {
        let jm = grid.idx_y.prev(j);
        for i in 0..nx {
            let im = grid.idx_x.prev(i);
            let lin = |ii: usize, jj: usize| jj * nx + ii;
            // (B^T p) at vx(i,j) = (p[im,j] - p[i,j]) / dx
            out_vx[lin(i, j)] = (p[lin(im, j)] - p[lin(i, j)]) * inv_dx;
            // (B^T p) at vy(i,j) = (p[i,jm] - p[i,j]) / dy
            out_vy[lin(i, j)] = (p[lin(i, jm)] - p[lin(i, j)]) * inv_dy;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dot(a: &[f64], b: &[f64]) -> f64 {
        a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
    }

    #[test]
    fn harmonic_of_equal_values_equals_the_value() {
        assert!((harmonic4(2.5, 2.5, 2.5, 2.5) - 2.5).abs() < 1e-14);
        assert!((harmonic4(1.0, 1.0, 1.0, 1.0) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn momentum_on_zero_is_zero() {
        let grid = StokesGrid::new(8, 8, 0.125, 0.125);
        let eta = Field2D::filled(8, 8, 1.0);
        let vx = vec![0.0; 64];
        let vy = vec![0.0; 64];
        let mut out_vx = vec![9.9; 64];
        let mut out_vy = vec![9.9; 64];
        apply_momentum(&grid, &eta, &vx, &vy, &mut out_vx, &mut out_vy);
        for v in out_vx.iter().chain(out_vy.iter()) {
            assert_eq!(*v, 0.0);
        }
    }

    #[test]
    fn momentum_diagonal_for_constant_eta_is_6() {
        // With η = dx = dy = 1, the diagonal of A at every velocity DOF is 6.
        let grid = StokesGrid::new(8, 8, 1.0, 1.0);
        let eta = Field2D::filled(8, 8, 1.0);
        let mut dvx = vec![0.0; 64];
        let mut dvy = vec![0.0; 64];
        momentum_diagonal(&grid, &eta, &mut dvx, &mut dvy);
        for (k, (&a, &b)) in dvx.iter().zip(dvy.iter()).enumerate() {
            assert!((a - 6.0).abs() < 1e-12, "diag_vx[{}] = {}", k, a);
            assert!((b - 6.0).abs() < 1e-12, "diag_vy[{}] = {}", k, b);
        }
    }

    #[test]
    fn momentum_is_symmetric() {
        // A symmetric ⇒ ⟨A u, w⟩ = ⟨u, A w⟩ for all u, w.
        let nx = 8;
        let ny = 8;
        let grid = StokesGrid::new(nx, ny, 0.13, 0.17);
        // Non-uniform but positive η to exercise harmonic averaging.
        let mut eta = Field2D::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                let e = 1.0 + 0.5 * ((i * 3 + j * 7) % 5) as f64 / 5.0;
                eta.set(i, j, e);
            }
        }
        let n2 = nx * ny;
        let mut ux = vec![0.0; n2];
        let mut uy = vec![0.0; n2];
        let mut wx = vec![0.0; n2];
        let mut wy = vec![0.0; n2];
        for k in 0..n2 {
            ux[k] = ((k as f64 * 1.7).sin()) * 1.1;
            uy[k] = ((k as f64 * 2.3).cos()) * 0.7;
            wx[k] = ((k as f64 * 0.9).sin()) * 0.5;
            wy[k] = ((k as f64 * 1.3).cos()) * 1.3;
        }
        let mut aux_x = vec![0.0; n2];
        let mut aux_y = vec![0.0; n2];
        let mut awx = vec![0.0; n2];
        let mut awy = vec![0.0; n2];
        apply_momentum(&grid, &eta, &ux, &uy, &mut aux_x, &mut aux_y);
        apply_momentum(&grid, &eta, &wx, &wy, &mut awx, &mut awy);
        let lhs = dot(&aux_x, &wx) + dot(&aux_y, &wy);
        let rhs = dot(&ux, &awx) + dot(&uy, &awy);
        let rel = (lhs - rhs).abs() / lhs.abs().max(rhs.abs()).max(1.0);
        assert!(rel < 1e-12, "symmetry broken: |lhs-rhs|/max = {}", rel);
    }

    #[test]
    fn divergence_transpose_is_adjoint_of_divergence() {
        // ⟨B v, p⟩ = ⟨v, B^T p⟩ for random v and p.
        let nx = 7;
        let ny = 5;
        let grid = StokesGrid::new(nx, ny, 0.19, 0.23);
        let n2 = nx * ny;
        let mut vx = vec![0.0; n2];
        let mut vy = vec![0.0; n2];
        let mut p = vec![0.0; n2];
        for k in 0..n2 {
            vx[k] = ((k as f64 * 1.1).sin()) * 1.3;
            vy[k] = ((k as f64 * 1.7).cos()) * 0.9;
            p[k] = ((k as f64 * 2.1).sin()) * 1.7;
        }
        let mut div_v = vec![0.0; n2];
        apply_divergence(&grid, &vx, &vy, &mut div_v);
        let mut bt_x = vec![0.0; n2];
        let mut bt_y = vec![0.0; n2];
        apply_divergence_transpose(&grid, &p, &mut bt_x, &mut bt_y);
        let lhs = dot(&div_v, &p);
        let rhs = dot(&vx, &bt_x) + dot(&vy, &bt_y);
        let rel = (lhs - rhs).abs() / lhs.abs().max(rhs.abs()).max(1.0);
        assert!(rel < 1e-12, "B^T not adjoint of B: rel = {}", rel);
    }

    #[test]
    fn divergence_of_constant_is_zero() {
        let grid = StokesGrid::new(6, 6, 0.25, 0.25);
        let vx = vec![3.14; 36];
        let vy = vec![-2.71; 36];
        let mut d = vec![9.9; 36];
        apply_divergence(&grid, &vx, &vy, &mut d);
        for v in d {
            assert!(v.abs() < 1e-14);
        }
    }
}
