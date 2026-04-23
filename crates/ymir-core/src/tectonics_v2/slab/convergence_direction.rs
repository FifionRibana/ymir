//! Per-cell unit vector pointing "towards convergence".
//!
//! Physical meaning (§4.8): the slab-pull traction at a surface
//! cell should point along the direction in which the subducted
//! material is hanging. The surface proxy chosen here is
//!
//! ```text
//!   n̂_convergence(x) = -∇(div v)(x) / |∇(div v)(x)|
//! ```
//!
//! `-∇(div v)` points **from divergent regions toward convergent
//! regions** — a cell adjacent to a subduction zone (where
//! `div v < 0`) feels a vector pointing into the zone, which is
//! the horizontal component of the slab's pulling direction.
//!
//! In quiescent cells (|∇(div v)| below a user-supplied epsilon)
//! we fall back to the zero vector rather than divide a noisy
//! gradient and get a unit vector that flips randomly every
//! step. This is the "undefined direction" case; the matching
//! `f_slab = Sp · m · n̂` contribution is zero, which is the
//! physically right answer.
//!
//! The output is **cell-centered**; the face interpolation happens
//! inside [`super::super::forcing::slab_pull::SlabPullForce`]
//! (same pattern as `GpeForce`).

use super::super::boundaries::source_sink::div_v_cell;
use super::super::field::{Field2D, PeriodicIndex};

/// Knobs controlling `n̂_convergence` computation.
#[derive(Clone, Copy, Debug)]
pub struct ConvergenceDirectionConfig {
    /// Threshold below which `|∇(div v)|` is treated as "no
    /// definable direction"; `n̂` is set to `(0, 0)` on such
    /// cells.
    pub epsilon: f64,
}

impl Default for ConvergenceDirectionConfig {
    fn default() -> Self {
        Self { epsilon: super::EPSILON_DEFAULT }
    }
}

/// Compute `n̂_convergence` in-place, writing the x-component into
/// `n_x` and y-component into `n_y` (both cell-centered, shape
/// `nx × ny`).
///
/// Steps:
/// 1. `div_v` at each cell from the MAC faces (delegates to
///    [`div_v_cell`]).
/// 2. Cell-centered centred differences of `div_v` →
///    `g = ∇(div v)`.
/// 3. Normalise `-g`, falling back to `(0, 0)` when `|g| < ε`.
///
/// The caller owns the output buffers and a scratch `div_v` buffer
/// to avoid per-step allocations in the time loop.
pub fn compute_convergence_direction(
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    idx_x: &PeriodicIndex,
    idx_y: &PeriodicIndex,
    vx: &[f64],
    vy: &[f64],
    div_scratch: &mut Field2D,
    n_x: &mut Field2D,
    n_y: &mut Field2D,
    config: &ConvergenceDirectionConfig,
) {
    debug_assert_eq!(div_scratch.nx(), nx);
    debug_assert_eq!(div_scratch.ny(), ny);
    debug_assert_eq!(n_x.nx(), nx);
    debug_assert_eq!(n_x.ny(), ny);
    debug_assert_eq!(n_y.nx(), nx);
    debug_assert_eq!(n_y.ny(), ny);

    // Pass 1: div(v) at each cell centre.
    div_v_cell(nx, ny, dx, dy, idx_x, idx_y, vx, vy, div_scratch);

    // Pass 2: centred difference ∇(div v), normalise −g.
    let inv_2dx = 0.5 / dx;
    let inv_2dy = 0.5 / dy;
    let eps = config.epsilon;
    let eps2 = eps * eps;

    for j in 0..ny {
        let jp = idx_y.next(j);
        let jm = idx_y.prev(j);
        for i in 0..nx {
            let ip = idx_x.next(i);
            let im = idx_x.prev(i);

            let gx = (div_scratch.get(ip, j) - div_scratch.get(im, j)) * inv_2dx;
            let gy = (div_scratch.get(i, jp) - div_scratch.get(i, jm)) * inv_2dy;

            let mag2 = gx * gx + gy * gy;
            if mag2 > eps2 {
                let inv_mag = 1.0 / mag2.sqrt();
                // n̂ = -g / |g|
                n_x.set(i, j, -gx * inv_mag);
                n_y.set(i, j, -gy * inv_mag);
            } else {
                n_x.set(i, j, 0.0);
                n_y.set(i, j, 0.0);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Uniform velocity ⇒ div(v) = 0 ⇒ ∇(div v) = 0 ⇒ n̂ = 0.
    #[test]
    fn uniform_flow_gives_zero_direction() {
        let nx = 8;
        let ny = 8;
        let dx = 1.0 / nx as f64;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);

        let vx = vec![0.3; nx * ny];
        let vy = vec![-0.1; nx * ny];

        let mut div = Field2D::new(nx, ny);
        let mut n_x = Field2D::new(nx, ny);
        let mut n_y = Field2D::new(nx, ny);
        compute_convergence_direction(
            nx,
            ny,
            dx,
            dx,
            &idx_x,
            &idx_y,
            &vx,
            &vy,
            &mut div,
            &mut n_x,
            &mut n_y,
            &ConvergenceDirectionConfig::default(),
        );
        for k in 0..nx * ny {
            assert_eq!(n_x.data()[k], 0.0);
            assert_eq!(n_y.data()[k], 0.0);
        }
    }

    /// A 1D sinusoidal convergence field has `div(v) = cos(2πx)`
    /// (after suitable choice of `vx`). Its gradient is
    /// `-2π sin(2πx)` in x, zero in y. The unit vector `n̂`
    /// should have x-component `+sign(sin(2πx))` (since
    /// `n̂ = -∇(div)/|∇(div)|`) and y-component 0 — modulo the
    /// epsilon fallback near `sin = 0`.
    #[test]
    fn one_d_convergence_points_along_x() {
        let nx = 32;
        let ny = 8;
        let dx = 1.0 / nx as f64;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);

        // Build vx on vertical faces s.t. div(v) = sin(2πx) at cell centre.
        // On the MAC grid with vx at face (i, j+½), div = (vx[i+1] - vx[i])/dx.
        // We want this to equal sin(2π (i+0.5) dx). Choose
        // vx(x_face) = -cos(2π x_face) / (2π); then (vx[i+1]-vx[i])/dx
        // approximates sin(2π (i+0.5) dx) at O(dx²).
        let mut vx = vec![0.0; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64 * dx;
                vx[j * nx + i] =
                    -(2.0 * std::f64::consts::PI * x).cos() / (2.0 * std::f64::consts::PI);
            }
        }
        let vy = vec![0.0; nx * ny];

        let mut div = Field2D::new(nx, ny);
        let mut n_x = Field2D::new(nx, ny);
        let mut n_y = Field2D::new(nx, ny);
        compute_convergence_direction(
            nx,
            ny,
            dx,
            dx,
            &idx_x,
            &idx_y,
            &vx,
            &vy,
            &mut div,
            &mut n_x,
            &mut n_y,
            &ConvergenceDirectionConfig::default(),
        );

        // At a cell where sin(2π (i+0.5) dx) is far from zero,
        // ∇(div) is dominated by its x-component (a cosine),
        // so n̂ should have |n_x| ≈ 1 and |n_y| ≈ 0.
        let i_probe = 8; // x ≈ 0.265625 → cos(2π x) ≈ 0.195 ≠ 0
        let j_probe = 3;
        let nx_val = n_x.get(i_probe, j_probe);
        let ny_val = n_y.get(i_probe, j_probe);
        assert!(nx_val.abs() > 0.99, "|n_x| = {}, expected ≈ 1", nx_val.abs());
        assert!(ny_val.abs() < 1e-10, "|n_y| = {}, expected ≈ 0", ny_val.abs());
    }

    /// With a huge `epsilon`, every cell falls back to zero.
    #[test]
    fn large_epsilon_zeros_direction() {
        let nx = 16;
        let ny = 16;
        let dx = 1.0 / nx as f64;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);

        let mut vx = vec![0.0; nx * ny];
        for (k, v) in vx.iter_mut().enumerate() {
            *v = (k as f64 * 0.1).sin();
        }
        let vy = vec![0.0; nx * ny];

        let mut div = Field2D::new(nx, ny);
        let mut n_x = Field2D::new(nx, ny);
        let mut n_y = Field2D::new(nx, ny);
        compute_convergence_direction(
            nx,
            ny,
            dx,
            dx,
            &idx_x,
            &idx_y,
            &vx,
            &vy,
            &mut div,
            &mut n_x,
            &mut n_y,
            &ConvergenceDirectionConfig { epsilon: 1.0e10 },
        );
        for k in 0..nx * ny {
            assert_eq!(n_x.data()[k], 0.0);
            assert_eq!(n_y.data()[k], 0.0);
        }
    }

    /// Unit-vector property: wherever fallback did not trigger,
    /// |n̂|² should equal 1 exactly.
    #[test]
    fn nonzero_direction_is_unit_magnitude() {
        let nx = 16;
        let ny = 16;
        let dx = 1.0 / nx as f64;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);

        // Arbitrary smooth but non-uniform 2D field: vx = sin(2πx)·cos(2πy)
        // on faces, giving a non-trivial div(v).
        let mut vx = vec![0.0; nx * ny];
        let mut vy = vec![0.0; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64 * dx;
                let y = (j as f64 + 0.5) * dx;
                vx[j * nx + i] =
                    (2.0 * std::f64::consts::PI * x).sin() * (2.0 * std::f64::consts::PI * y).cos();
                let xp = (i as f64 + 0.5) * dx;
                let yp = j as f64 * dx;
                vy[j * nx + i] = -(2.0 * std::f64::consts::PI * xp).cos()
                    * (2.0 * std::f64::consts::PI * yp).sin();
            }
        }

        let mut div = Field2D::new(nx, ny);
        let mut n_x = Field2D::new(nx, ny);
        let mut n_y = Field2D::new(nx, ny);
        compute_convergence_direction(
            nx,
            ny,
            dx,
            dx,
            &idx_x,
            &idx_y,
            &vx,
            &vy,
            &mut div,
            &mut n_x,
            &mut n_y,
            &ConvergenceDirectionConfig { epsilon: 1.0e-8 },
        );

        for k in 0..nx * ny {
            let mag2 = n_x.data()[k] * n_x.data()[k] + n_y.data()[k] * n_y.data()[k];
            // Each cell: either fallback (mag² = 0) or unit (mag² ≈ 1).
            let ok = mag2 < 1.0e-20 || (mag2 - 1.0).abs() < 1.0e-12;
            assert!(ok, "cell {}: mag² = {}", k, mag2);
        }
    }
}
