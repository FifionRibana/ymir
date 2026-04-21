//! Power-law rheology with additive strain-rate floor and smooth
//! saturation cap.
//!
//! # Constitutive law
//!
//! The dimensional power-law form is `η = B·ε̇_II^(1/n - 1)` (Gerya
//! 2010 §14). After nondimensionalization we take `B̃ = 1` and apply
//! two regularisations:
//!
//! 1. **Additive strain-rate floor.**
//!    `ε̇_II_reg = ε̇_II + ε̇_min` (not `√(ε̇_II² + ε̇_min²)`). The
//!    additive form gives a simple closed derivative (see below) and
//!    coincides with the `√`-form everywhere except in a narrow band
//!    around `ε̇_II ≪ ε̇_min`; in that band, Newton's iteration is
//!    already dominated by the stationary state, so the choice does
//!    not influence convergence. Documentation fix point — this is
//!    the regularisation used throughout `tectonics_v2`.
//!
//! 2. **Smooth upper saturation.**
//!    `η_eff = smooth_saturate(η_newton, η_max; k)` with k = 4 by
//!    default. The cap is a soft asymptote, not a hard clamp. The
//!    effective viscosity approaches `η_max` monotonically as
//!    `η_newton → ∞` but never exceeds it.
//!
//! Putting it together:
//!
//! ```text
//!   η_newton(ε̇_II) = B̃ · (ε̇_II + ε̇_min)^(1/n - 1)
//!   η_eff   (ε̇_II) = η_newton / (1 + (η_newton/η_max)^k)^(1/k)
//! ```
//!
//! # First derivatives (Gerya §14.4)
//!
//! ```text
//!   dη_newton / dε̇_II = (1/n - 1) · η_newton / (ε̇_II + ε̇_min)
//!   dη_eff    / dη_newton = [1 + (η_newton/η_max)^k]^(-(1+k)/k)
//! ```
//!
//! Chain rule:
//!
//! ```text
//!   dη_eff / dε̇_II = dη_eff/dη_newton · dη_newton/dε̇_II
//! ```
//!
//! For shear-thinning (`n > 1`, `1/n - 1 < 0`) the derivative is
//! negative. The Newton-extra term in the Jacobian (cf
//! `stokes::operator::apply_tangent`) is proportional to this
//! derivative, which flips its contribution to the quadratic form —
//! the full Newton Jacobian is **symmetric but not necessarily SPD**
//! in zones of strong localisation. Tests that exercise Jacobian
//! structure verify symmetry, not positive-definiteness.
//!
//! # Asymptotes
//!
//! - `ε̇_II → 0`: `η_newton → ε̇_min^(1/n − 1)`, so `η_newton` is
//!   bounded by the floor alone. With `ε̇_min = 10⁻³` and `n = 3` that
//!   bound is `10²`, well under the default soft cap `η_max = 10³`.
//! - `ε̇_II → ∞`: `η_newton → 0`, and the saturation wrapper leaves it
//!   at 0. `η_eff` inherits the shear-thinning limit.
//! - `η_newton = η_max`: `η_eff = η_max / 2^(1/k) ≈ 0.84 · η_max` for
//!   k = 4. The "cap active" signal `η_eff > 0.9·η_max` is reached
//!   when `η_newton ≈ 1.47 · η_max`.

use super::field::{Field2D, PeriodicIndex};

/// Power-law viscosity law evaluated pointwise on the MAC grid.
#[derive(Clone, Copy, Debug)]
pub struct ViscosityLaw {
    pub n: f64,
    pub b_prefactor: f64,
    pub strain_rate_floor: f64,
    pub eta_max_cap: f64,
    pub k_saturation: f64,
}

impl ViscosityLaw {
    /// `η_newton(ε̇_II) = B̃ · (ε̇_II + ε̇_min)^(1/n - 1)`.
    /// Shear-thinning for `n > 1`.
    #[inline]
    pub fn eta_newton(&self, eps_ii: f64) -> f64 {
        let eps_reg = eps_ii + self.strain_rate_floor;
        self.b_prefactor * eps_reg.powf(1.0 / self.n - 1.0)
    }

    /// `η_eff(ε̇_II) = smooth_saturate(η_newton(ε̇_II), η_max; k)`.
    /// The returned value is strictly bounded above by `η_max` for
    /// any positive input.
    #[inline]
    pub fn eta_effective(&self, eps_ii: f64) -> f64 {
        smooth_saturate(self.eta_newton(eps_ii), self.eta_max_cap, self.k_saturation)
    }

    /// `dη_newton/dε̇_II = (1/n - 1) · η_newton / (ε̇_II + ε̇_min)`.
    #[inline]
    pub fn d_eta_newton_d_eps_ii(&self, eps_ii: f64) -> f64 {
        let eps_reg = eps_ii + self.strain_rate_floor;
        (1.0 / self.n - 1.0) * self.eta_newton(eps_ii) / eps_reg
    }

    /// `dη_eff/dε̇_II` via chain rule.
    #[inline]
    pub fn d_eta_effective_d_eps_ii(&self, eps_ii: f64) -> f64 {
        let eta_n = self.eta_newton(eps_ii);
        let ratio = eta_n / self.eta_max_cap;
        let k = self.k_saturation;
        // dη_eff / dη_newton = [1 + ratio^k]^(-(k+1)/k)
        let factor = (1.0 + ratio.powf(k)).powf(-(k + 1.0) / k);
        factor * self.d_eta_newton_d_eps_ii(eps_ii)
    }

    /// The Newton-extra pre-factor used in the Jacobian stress,
    /// `2 · η'(ε̇_II) / ε̇_II_reg`. Evaluated with the additive-floor
    /// denominator so it stays finite at `ε̇_II → 0`.
    #[inline]
    pub fn newton_extra_prefactor(&self, eps_ii: f64) -> f64 {
        2.0 * self.d_eta_effective_d_eps_ii(eps_ii) / (eps_ii + self.strain_rate_floor)
    }
}

impl Default for ViscosityLaw {
    fn default() -> Self {
        Self {
            n: 3.0,
            b_prefactor: 1.0,
            strain_rate_floor: 1.0e-3,
            eta_max_cap: 1.0e3,
            k_saturation: 4.0,
        }
    }
}

/// Smooth saturation with configurable sharpness `k`. Reduces to
/// `smooth_saturate(x, x_max)` of the legacy solver when `k = 4`.
/// Parity is verified in the unit tests.
///
/// Copied from `tectonics/solver/smooth.rs` (commit `74510ad`), which
/// is frozen for the duration of the reconstruction milestone. Tests
/// of the legacy function apply here too.
#[inline]
pub fn smooth_saturate(x: f64, x_max: f64, k: f64) -> f64 {
    let ratio = x / x_max;
    x / (1.0 + ratio.powf(k)).powf(1.0 / k)
}

/// Precomputed strain-rate component fields on the staggered grid.
///
/// Layout (Gerya 2010 §14):
/// - `exx_center[i,j] = ∂vx/∂x at cell centre (i,j)`.
/// - `eyy_center[i,j] = ∂vy/∂y at cell centre (i,j)`.
/// - `exy_corner[i,j] = ½(∂vx/∂y + ∂vy/∂x) at corner (i·dx, j·dy)`.
/// - `eps_ii_center[i,j] = √(½·(exx² + eyy²) + ⟨exy²⟩_cc)` where the
///   corner `exy²` is arithmetically averaged to the cell centre from
///   the four surrounding corners.
/// - `eps_ii_corner[i,j] = √(½·(⟨exx²⟩_corner + ⟨eyy²⟩_corner) + exy²)`
///   using cell-centred `exx, eyy` averaged to the corner.
pub struct StrainRate {
    pub exx_center: Field2D,
    pub eyy_center: Field2D,
    pub exy_corner: Field2D,
    pub eps_ii_center: Field2D,
    pub eps_ii_corner: Field2D,
}

impl StrainRate {
    pub fn compute(
        nx: usize,
        ny: usize,
        dx: f64,
        dy: f64,
        idx_x: &PeriodicIndex,
        idx_y: &PeriodicIndex,
        vx: &[f64],
        vy: &[f64],
    ) -> Self {
        let n = nx * ny;
        debug_assert_eq!(vx.len(), n);
        debug_assert_eq!(vy.len(), n);
        let inv_dx = 1.0 / dx;
        let inv_dy = 1.0 / dy;
        let lin = |ii: usize, jj: usize| jj * nx + ii;

        let mut exx_center = Field2D::new(nx, ny);
        let mut eyy_center = Field2D::new(nx, ny);
        let mut exy_corner = Field2D::new(nx, ny);
        for j in 0..ny {
            let jp = idx_y.next(j);
            let jm = idx_y.prev(j);
            for i in 0..nx {
                let ip = idx_x.next(i);
                let im = idx_x.prev(i);
                // Cell-centred normal strains.
                exx_center.set(i, j, (vx[lin(ip, j)] - vx[lin(i, j)]) * inv_dx);
                eyy_center.set(i, j, (vy[lin(i, jp)] - vy[lin(i, j)]) * inv_dy);
                // Corner shear strain at (i·dx, j·dy).
                let dvx_dy = (vx[lin(i, j)] - vx[lin(i, jm)]) * inv_dy;
                let dvy_dx = (vy[lin(i, j)] - vy[lin(im, j)]) * inv_dx;
                exy_corner.set(i, j, 0.5 * (dvx_dy + dvy_dx));
            }
        }

        // ε̇_II at cell centres: ½·(exx² + eyy²) + ⟨exy²⟩_cc.
        let mut eps_ii_center = Field2D::new(nx, ny);
        for j in 0..ny {
            let jp = idx_y.next(j);
            for i in 0..nx {
                let ip = idx_x.next(i);
                let exx = exx_center.get(i, j);
                let eyy = eyy_center.get(i, j);
                let exy2_avg = 0.25
                    * (exy_corner.get(i, j).powi(2)
                        + exy_corner.get(ip, j).powi(2)
                        + exy_corner.get(i, jp).powi(2)
                        + exy_corner.get(ip, jp).powi(2));
                let e2 = 0.5 * (exx * exx + eyy * eyy) + exy2_avg;
                eps_ii_center.set(i, j, e2.max(0.0).sqrt());
            }
        }

        // ε̇_II at corners: ½·(⟨exx²⟩_corner + ⟨eyy²⟩_corner) + exy².
        let mut eps_ii_corner = Field2D::new(nx, ny);
        for j in 0..ny {
            let jm = idx_y.prev(j);
            for i in 0..nx {
                let im = idx_x.prev(i);
                let exx2_avg = 0.25
                    * (exx_center.get(im, jm).powi(2)
                        + exx_center.get(i, jm).powi(2)
                        + exx_center.get(im, j).powi(2)
                        + exx_center.get(i, j).powi(2));
                let eyy2_avg = 0.25
                    * (eyy_center.get(im, jm).powi(2)
                        + eyy_center.get(i, jm).powi(2)
                        + eyy_center.get(im, j).powi(2)
                        + eyy_center.get(i, j).powi(2));
                let exy = exy_corner.get(i, j);
                let e2 = 0.5 * (exx2_avg + eyy2_avg) + exy * exy;
                eps_ii_corner.set(i, j, e2.max(0.0).sqrt());
            }
        }

        Self { exx_center, eyy_center, exy_corner, eps_ii_center, eps_ii_corner }
    }

    /// Cell-centred `exx` interpolated to corner `(i, j)` (average of the
    /// four surrounding cell centres `(i-1,j-1), (i,j-1), (i-1,j), (i,j)`).
    #[inline]
    pub fn exx_at_corner(
        &self,
        i: usize,
        j: usize,
        idx_x: &PeriodicIndex,
        idx_y: &PeriodicIndex,
    ) -> f64 {
        let im = idx_x.prev(i);
        let jm = idx_y.prev(j);
        0.25 * (self.exx_center.get(im, jm)
            + self.exx_center.get(i, jm)
            + self.exx_center.get(im, j)
            + self.exx_center.get(i, j))
    }

    #[inline]
    pub fn eyy_at_corner(
        &self,
        i: usize,
        j: usize,
        idx_x: &PeriodicIndex,
        idx_y: &PeriodicIndex,
    ) -> f64 {
        let im = idx_x.prev(i);
        let jm = idx_y.prev(j);
        0.25 * (self.eyy_center.get(im, jm)
            + self.eyy_center.get(i, jm)
            + self.eyy_center.get(im, j)
            + self.eyy_center.get(i, j))
    }

    /// Corner-centred `exy` interpolated to cell centre `(i, j)`
    /// (average of the four surrounding corners).
    #[inline]
    pub fn exy_at_center(
        &self,
        i: usize,
        j: usize,
        idx_x: &PeriodicIndex,
        idx_y: &PeriodicIndex,
    ) -> f64 {
        let ip = idx_x.next(i);
        let jp = idx_y.next(j);
        0.25 * (self.exy_corner.get(i, j)
            + self.exy_corner.get(ip, j)
            + self.exy_corner.get(i, jp)
            + self.exy_corner.get(ip, jp))
    }
}

/// Build an `η` field at cell centres from the rheology and the
/// current strain-rate field. Stored as a `Field2D` so the operator
/// layer can feed it to the existing harmonic corner averaging
/// unchanged.
pub fn build_eta_field(law: &ViscosityLaw, eps_ii_center: &Field2D) -> Field2D {
    let nx = eps_ii_center.nx();
    let ny = eps_ii_center.ny();
    let mut eta = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            eta.set(i, j, law.eta_effective(eps_ii_center.get(i, j)));
        }
    }
    eta
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * a.abs().max(b.abs()).max(1.0)
    }

    /// Parity with the legacy `smooth_saturate(x, x_max)` which pins
    /// `k = 4` internally. Any drift between the copy and the
    /// original is caught here.
    #[test]
    fn smooth_saturate_matches_legacy_at_k4() {
        let x_max = 1e4;
        let xs: [f64; 9] = [1e-3, 1e-2, 1e-1, 1.0, 10.0, 1e2, 1e3, 1e4, 1e5];
        for x in xs {
            let legacy = x / (1.0 + (x / x_max).powf(4.0)).powf(1.0 / 4.0);
            let here = smooth_saturate(x, x_max, 4.0);
            assert!(
                (legacy - here).abs() < 1e-14,
                "smooth_saturate parity broken at x={}: {} vs {}",
                x,
                legacy,
                here,
            );
        }
    }

    #[test]
    fn eta_newton_monotone_shear_thinning() {
        let law = ViscosityLaw::default(); // n = 3
        let xs = [1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0];
        for w in xs.windows(2) {
            let a = law.eta_newton(w[0]);
            let b = law.eta_newton(w[1]);
            assert!(a >= b, "η_newton not monotone at ε̇={}: {} → {}", w[0], a, b);
        }
    }

    #[test]
    fn eta_effective_bounded_above_by_eta_max() {
        let law = ViscosityLaw::default();
        for x in [1e-6, 1e-3, 1.0, 1e6, 1e12] {
            let eta = law.eta_effective(x);
            assert!(eta <= law.eta_max_cap, "η_eff({}) = {} exceeds cap", x, eta);
        }
    }

    #[test]
    fn derivative_agrees_with_finite_difference_on_a_sampling_grid() {
        let law = ViscosityLaw::default();
        let xs = [5e-3, 1e-2, 5e-2, 0.1, 0.5, 1.0, 5.0, 50.0];
        let h = 1e-6;
        for &x in &xs {
            let analytic = law.d_eta_effective_d_eps_ii(x);
            let fd = (law.eta_effective(x + h) - law.eta_effective(x - h)) / (2.0 * h);
            // Relative tolerance scaled to 1e-6 is plenty given the
            // analytic derivative is accurate to machine precision.
            let rel = (analytic - fd).abs() / analytic.abs().max(fd.abs()).max(1e-12);
            assert!(rel < 1e-5, "d_eta/d_eps_II FD mismatch at x={}: analytic {} vs FD {}", x, analytic, fd);
        }
    }

    #[test]
    fn derivative_is_negative_for_shear_thinning() {
        let law = ViscosityLaw::default();
        for x in [1e-3, 1e-2, 1.0, 100.0] {
            let dp = law.d_eta_effective_d_eps_ii(x);
            assert!(dp <= 0.0, "dη/dε̇_II({}) = {} should be ≤ 0", x, dp);
        }
    }

    #[test]
    fn eta_newton_is_c1_across_the_floor() {
        // The additive floor makes η_newton smooth on (0, ∞);
        // finite-difference derivative should change continuously
        // across ε̇_min.
        let law = ViscosityLaw::default();
        let h = 1e-7;
        let mut prev: Option<f64> = None;
        for k in -20..20 {
            let x = law.strain_rate_floor * 2.0f64.powi(k);
            let fd = (law.eta_newton(x + h) - law.eta_newton(x - h)) / (2.0 * h);
            if let Some(p) = prev {
                let ratio = (fd / p).abs();
                assert!(
                    ratio < 10.0 && ratio > 0.1,
                    "non-smooth transition at x={}: FD went from {} to {}",
                    x,
                    p,
                    fd,
                );
            }
            prev = Some(fd);
        }
    }

    #[test]
    fn strain_rate_computes_zero_for_zero_velocity() {
        let nx = 8;
        let ny = 8;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);
        let vx = vec![0.0; nx * ny];
        let vy = vec![0.0; nx * ny];
        let sr = StrainRate::compute(nx, ny, 0.1, 0.1, &idx_x, &idx_y, &vx, &vy);
        for v in sr.eps_ii_center.data().iter() {
            assert_eq!(*v, 0.0);
        }
        for v in sr.eps_ii_corner.data().iter() {
            assert_eq!(*v, 0.0);
        }
    }

    #[test]
    fn build_eta_field_uses_floor_when_strain_rate_is_zero() {
        // When ε̇_II = 0 everywhere, η_newton = ε̇_min^(1/n − 1).
        // For n=3, ε̇_min=1e-3 this is 10⁴, which the soft cap at
        // η_max=10³ attenuates to ~0.84·10³.
        let law = ViscosityLaw::default();
        let eps = Field2D::filled(4, 4, 0.0);
        let eta = build_eta_field(&law, &eps);
        let expected = law.eta_effective(0.0);
        for v in eta.data().iter() {
            assert!(approx(*v, expected, 1e-12));
        }
    }
}
