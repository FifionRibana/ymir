//! Gravitational potential energy (GPE) spreading force.
//!
//! Continuous form (nondim, thin viscous sheet):
//! ```text
//!   f̃_GPE = -Ar · ∇(½ S̃²) = -Ar · S̃ · ∇S̃
//! ```
//! Discrete form on the MAC grid (staggered standard, consistent
//! with Step 0's `apply_momentum`):
//! ```text
//!   f_x[vx face (i dx, (j+0.5) dy)]
//!     = -Ar · ½ · (S(i, j)² - S(i-1, j)²) / dx
//!     = -Ar · ½ · (S(i, j) + S(i-1, j)) · (S(i, j) - S(i-1, j)) / dx
//! ```
//! and symmetrically for `f_y` at the horizontal faces. The product-
//! and-sum expansion of the central term reveals that the operator
//! is exactly a divided-difference of the potential `½ S²` between
//! the two cell centres on either side of the face.
//!
//! # What this formulation does NOT solve
//!
//! Issue #78 (GPE gradient spike at sharp material interfaces) is
//! **not** addressed here. At Step 2 the thickness field stays in a
//! narrow band around `S ≈ 1` (advected by a gentle placeholder
//! forcing with no boundary sources), so `|∇S|` remains `O(1)` and
//! the staggered differencing is well-behaved. When Steps 5/6
//! introduce sharp oceanic/continental interfaces, `|∇S|` locally
//! blows up to `O(1/dx)` and this formulation — like any pointwise
//! discretisation of `S·∇S` — will spike. Real cures (wider stencil,
//! separate smoothed-`S` for GPE, semi-implicit spreading) are
//! explicit Step-5/6 work and live out of scope here.

use super::super::scales::Scales;
use super::body_force::{BodyForce, SimulationState, VectorField};

/// Nondim GPE forcing term. The Argand number is **derived** from
/// the 4 primary scales via [`Scales::argand_number`]; it is not an
/// independent knob. Ar is stored on construction and re-used.
#[derive(Clone, Copy, Debug)]
pub struct GpeForce {
    pub ar: f64,
}

impl GpeForce {
    /// Build a `GpeForce` whose `Ar` comes from the solver's scales.
    pub fn from_scales(scales: &Scales) -> Self {
        Self { ar: scales.argand_number() }
    }

    /// Explicit constructor. Prefer `from_scales`; this is provided
    /// for tests that want to probe Ar-sensitivity.
    pub fn with_ar(ar: f64) -> Self {
        Self { ar }
    }
}

impl BodyForce for GpeForce {
    fn accumulate(&self, state: &SimulationState, out: &mut VectorField) {
        let nx = state.nx;
        let ny = state.ny;
        let inv_dx = 1.0 / state.dx;
        let inv_dy = 1.0 / state.dy;
        let s = state.s;
        let ar = self.ar;

        // x-component: vx face (i, j+½) sees potentials at cell
        // centres (i-1, j) and (i, j).
        for j in 0..ny {
            for i in 0..nx {
                let im = state.idx_x.prev(i);
                let s_right = s.get(i, j);
                let s_left = s.get(im, j);
                let f = -ar * 0.5 * (s_right + s_left) * (s_right - s_left) * inv_dx;
                let k = j * nx + i;
                out.fx.data_mut()[k] += f;
            }
        }
        // y-component: vy face (i+½, j) sees potentials at cell
        // centres (i, j-1) and (i, j).
        for j in 0..ny {
            for i in 0..nx {
                let jm = state.idx_y.prev(j);
                let s_top = s.get(i, j);
                let s_bot = s.get(i, jm);
                let f = -ar * 0.5 * (s_top + s_bot) * (s_top - s_bot) * inv_dy;
                let k = j * nx + i;
                out.fy.data_mut()[k] += f;
            }
        }
    }

    fn name(&self) -> &'static str {
        "GpeForce"
    }
}

#[cfg(test)]
mod tests {
    use super::super::super::field::{Field2D, PeriodicIndex};
    use super::*;

    fn linear_idx(nx: usize, i: usize, j: usize) -> usize {
        j * nx + i
    }

    fn state_env(
        nx: usize,
        ny: usize,
        _dx: f64,
        _dy: f64,
        s_fill: impl Fn(usize, usize) -> f64,
    ) -> (PeriodicIndex, PeriodicIndex, Field2D) {
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);
        let mut s = Field2D::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                s.set(i, j, s_fill(i, j));
            }
        }
        (idx_x, idx_y, s)
    }

    #[test]
    fn gpe_zero_when_thickness_uniform() {
        let (idx_x, idx_y, s) = state_env(8, 8, 0.125, 0.125, |_, _| 1.0);
        let mut fx = Field2D::new(8, 8);
        let mut fy = Field2D::new(8, 8);
        let st = SimulationState {
            nx: 8,
            ny: 8,
            dx: 0.125,
            dy: 0.125,
            idx_x: &idx_x,
            idx_y: &idx_y,
            s: &s,
        };
        let mut out = VectorField { fx: &mut fx, fy: &mut fy };
        GpeForce::with_ar(2.0).accumulate(&st, &mut out);
        for v in fx.data().iter().chain(fy.data().iter()) {
            assert_eq!(*v, 0.0);
        }
    }

    #[test]
    fn gpe_is_linear_in_ar() {
        let nx = 16;
        let ny = 16;
        let dx = 1.0 / nx as f64;
        let (idx_x, idx_y, s) = state_env(nx, ny, dx, dx, |i, j| {
            1.0 + 0.1
                * ((i as f64 / nx as f64) * std::f64::consts::TAU).sin()
                * ((j as f64 / ny as f64) * std::f64::consts::TAU).cos()
        });
        let st = SimulationState { nx, ny, dx, dy: dx, idx_x: &idx_x, idx_y: &idx_y, s: &s };
        let mut fx1 = Field2D::new(nx, ny);
        let mut fy1 = Field2D::new(nx, ny);
        let mut v1 = VectorField { fx: &mut fx1, fy: &mut fy1 };
        GpeForce::with_ar(1.0).accumulate(&st, &mut v1);

        let mut fx3 = Field2D::new(nx, ny);
        let mut fy3 = Field2D::new(nx, ny);
        let mut v3 = VectorField { fx: &mut fx3, fy: &mut fy3 };
        GpeForce::with_ar(3.0).accumulate(&st, &mut v3);

        for k in 0..(nx * ny) {
            let a = fx1.data()[k];
            let b = fx3.data()[k];
            assert!((b - 3.0 * a).abs() < 1e-14, "fx scale broken at k={}", k);
            let a2 = fy1.data()[k];
            let b2 = fy3.data()[k];
            assert!((b2 - 3.0 * a2).abs() < 1e-14, "fy scale broken at k={}", k);
        }
    }

    #[test]
    fn gpe_accumulates_additively() {
        let nx = 8;
        let ny = 8;
        let dx = 0.125;
        let (idx_x, idx_y, s) =
            state_env(nx, ny, dx, dx, |i, j| 1.0 + 0.2 * (i as f64 + j as f64).sin());
        let st = SimulationState { nx, ny, dx, dy: dx, idx_x: &idx_x, idx_y: &idx_y, s: &s };

        let gpe = GpeForce::with_ar(2.0);

        let mut fx_once = Field2D::new(nx, ny);
        let mut fy_once = Field2D::new(nx, ny);
        gpe.accumulate(&st, &mut VectorField { fx: &mut fx_once, fy: &mut fy_once });

        let mut fx_twice = Field2D::new(nx, ny);
        let mut fy_twice = Field2D::new(nx, ny);
        {
            let mut out = VectorField { fx: &mut fx_twice, fy: &mut fy_twice };
            gpe.accumulate(&st, &mut out);
            gpe.accumulate(&st, &mut out);
        }

        for k in 0..(nx * ny) {
            assert!((fx_twice.data()[k] - 2.0 * fx_once.data()[k]).abs() < 1e-14);
            assert!((fy_twice.data()[k] - 2.0 * fy_once.data()[k]).abs() < 1e-14);
        }
    }

    /// The integral of the GPE force over the periodic domain is
    /// zero — the potential `½ S²` is single-valued, so the divided
    /// difference sums to zero. Equivalently, `GpeForce` adds no net
    /// momentum to the domain, which is what "conservative" means at
    /// this granularity.
    #[test]
    fn gpe_integral_is_zero_on_periodic_domain() {
        let nx = 16;
        let ny = 16;
        let dx = 1.0 / nx as f64;
        let (idx_x, idx_y, s) = state_env(nx, ny, dx, dx, |i, j| {
            1.0 + 0.15 * ((i as f64 * 2.0).sin() + (j as f64 * 1.3).cos())
        });
        let st = SimulationState { nx, ny, dx, dy: dx, idx_x: &idx_x, idx_y: &idx_y, s: &s };
        let mut fx = Field2D::new(nx, ny);
        let mut fy = Field2D::new(nx, ny);
        GpeForce::with_ar(2.0).accumulate(&st, &mut VectorField { fx: &mut fx, fy: &mut fy });
        let sum_fx: f64 = fx.data().iter().sum();
        let sum_fy: f64 = fy.data().iter().sum();
        assert!(sum_fx.abs() < 1e-12, "Σfx = {}", sum_fx);
        assert!(sum_fy.abs() < 1e-12, "Σfy = {}", sum_fy);
    }

    /// For a smooth S field the staggered formulation agrees with the
    /// analytic `-Ar·S·∇S` at each face to 2nd order in dx — check
    /// the pointwise values at a single face.
    #[test]
    fn gpe_matches_analytic_on_smooth_field() {
        let nx = 64;
        let ny = 64;
        let dx = 1.0 / nx as f64;
        let k = 2.0 * std::f64::consts::PI;
        let s_of = |x: f64, y: f64| 1.0 + 0.1 * (k * x).sin() * (k * y).cos();
        let (idx_x, idx_y, s) = state_env(nx, ny, dx, dx, |i, j| {
            let x = (i as f64 + 0.5) * dx;
            let y = (j as f64 + 0.5) * dx;
            s_of(x, y)
        });
        let st = SimulationState { nx, ny, dx, dy: dx, idx_x: &idx_x, idx_y: &idx_y, s: &s };
        let mut fx = Field2D::new(nx, ny);
        let mut fy = Field2D::new(nx, ny);
        let ar = 2.0;
        GpeForce::with_ar(ar).accumulate(&st, &mut VectorField { fx: &mut fx, fy: &mut fy });

        // Probe at a face away from zero-gradient points.
        let i0 = 10;
        let j0 = 15;
        let xf = i0 as f64 * dx;
        let yf = (j0 as f64 + 0.5) * dx;
        // -Ar * S * ∂_x S at (xf, yf).
        let s_val = s_of(xf, yf);
        let dsdx = 0.1 * k * (k * xf).cos() * (k * yf).cos();
        let analytic = -ar * s_val * dsdx;
        let numeric = fx.data()[linear_idx(nx, i0, j0)];
        let rel = (analytic - numeric).abs() / analytic.abs().max(numeric.abs()).max(1e-12);
        assert!(rel < 2e-2, "analytic={}, numeric={}", analytic, numeric);
    }
}
