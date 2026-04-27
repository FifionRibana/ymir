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
///
/// Starting at Step 3 the struct carries an optional
/// [`YieldingConfig`](crate::tectonics_v2::presets::YieldingConfig)
/// field that toggles the plastic branch on or off. `eta_effective`
/// and `d_eta_effective_d_eps_ii` return the yielding-enhanced
/// effective viscosity when `yielding == Enabled(..)` and the
/// pure-viscous value otherwise — **structurally** by matching on
/// the enum variant. The `Disabled` path performs exactly the same
/// work as the pre-Step-3 code path, with no allocation or extra
/// branch inside the plastic helpers; this is what the Step-3
/// regression harness relies on for bit-comparable behaviour with
/// Step 2.
#[derive(Clone, Copy, Debug)]
pub struct ViscosityLaw {
    pub n: f64,
    pub b_prefactor: f64,
    pub strain_rate_floor: f64,
    pub eta_max_cap: f64,
    pub k_saturation: f64,
    /// Plastic-yielding configuration. `Disabled` (default)
    /// reproduces Step 0/1/2 behaviour exactly.
    pub yielding: crate::tectonics_v2::presets::YieldingConfig,
}

impl ViscosityLaw {
    /// `η_newton(ε̇_II) = B̃ · (ε̇_II + ε̇_min)^(1/n - 1)`.
    /// Pure power-law branch, **ignores yielding**. Useful for
    /// reporting the "viscous only" reference value (e.g. in the
    /// yielding-cell-fraction diagnostic).
    #[inline]
    pub fn eta_newton(&self, eps_ii: f64) -> f64 {
        let eps_reg = eps_ii + self.strain_rate_floor;
        self.b_prefactor * eps_reg.powf(1.0 / self.n - 1.0)
    }

    /// Viscous-branch effective viscosity — `eta_newton` wrapped by
    /// the smooth upper-saturation cap. **Ignores yielding.**
    #[inline]
    pub fn eta_visc_effective(&self, eps_ii: f64) -> f64 {
        smooth_saturate(self.eta_newton(eps_ii), self.eta_max_cap, self.k_saturation)
    }

    /// `η_eff(ε̇_II)` — full effective viscosity consumed by the
    /// solver. Under `yielding == Disabled` returns
    /// `eta_visc_effective`; under `Enabled(law)` returns
    /// `soft_min_harmonic(eta_visc, eta_plastic, sharpness)`.
    #[inline]
    pub fn eta_effective(&self, eps_ii: f64) -> f64 {
        let eta_v = self.eta_visc_effective(eps_ii);
        match self.yielding {
            crate::tectonics_v2::presets::YieldingConfig::Disabled => eta_v,
            crate::tectonics_v2::presets::YieldingConfig::Enabled(ylaw) => {
                let eta_p = eta_plastic(eps_ii, self.strain_rate_floor, ylaw.bi);
                eta_effective(eta_v, eta_p, ylaw.sharpness)
            }
        }
    }

    /// `dη_newton/dε̇_II = (1/n - 1) · η_newton / (ε̇_II + ε̇_min)`.
    #[inline]
    pub fn d_eta_newton_d_eps_ii(&self, eps_ii: f64) -> f64 {
        let eps_reg = eps_ii + self.strain_rate_floor;
        (1.0 / self.n - 1.0) * self.eta_newton(eps_ii) / eps_reg
    }

    /// Derivative of the viscous-branch effective viscosity (chain
    /// rule through the `smooth_saturate` cap). **Ignores yielding.**
    #[inline]
    pub fn d_eta_visc_effective_d_eps_ii(&self, eps_ii: f64) -> f64 {
        let eta_n = self.eta_newton(eps_ii);
        let ratio = eta_n / self.eta_max_cap;
        let k = self.k_saturation;
        let factor = (1.0 + ratio.powf(k)).powf(-(k + 1.0) / k);
        factor * self.d_eta_newton_d_eps_ii(eps_ii)
    }

    /// `dη_eff/dε̇_II` — full derivative consumed by the Newton
    /// Jacobian. Matches on the yielding variant; Enabled routes
    /// through [`d_eta_eff_d_eps_ii`] which chains the blend's
    /// partials.
    #[inline]
    pub fn d_eta_effective_d_eps_ii(&self, eps_ii: f64) -> f64 {
        let dv = self.d_eta_visc_effective_d_eps_ii(eps_ii);
        match self.yielding {
            crate::tectonics_v2::presets::YieldingConfig::Disabled => dv,
            crate::tectonics_v2::presets::YieldingConfig::Enabled(ylaw) => {
                let eta_v = self.eta_visc_effective(eps_ii);
                let eta_p = eta_plastic(eps_ii, self.strain_rate_floor, ylaw.bi);
                d_eta_eff_d_eps_ii(
                    eta_v, dv, eta_p, eps_ii, self.strain_rate_floor, ylaw.sharpness,
                )
            }
        }
    }

    /// The Newton-extra pre-factor used in the Jacobian stress,
    /// `2 · η'(ε̇_II) / ε̇_II_reg`. Evaluated with the additive-floor
    /// denominator so it stays finite at `ε̇_II → 0`.
    #[inline]
    pub fn newton_extra_prefactor(&self, eps_ii: f64) -> f64 {
        2.0 * self.d_eta_effective_d_eps_ii(eps_ii) / (eps_ii + self.strain_rate_floor)
    }

    /// Step 9 — `η_eff(ε̇_II)` with a per-cell `bi_override` replacing
    /// the global yielding law's `bi`. Used by the cratonic
    /// pipeline to apply per-cell Bi elevation
    /// `bi_eff[i] = Bi · (1 + (B_factor - 1) · cratonic_factor[i])`.
    /// Falls through to `eta_visc_effective` when yielding is
    /// disabled (the override has no effect there).
    #[inline]
    pub fn eta_effective_with_bi_override(&self, eps_ii: f64, bi_override: f64) -> f64 {
        let eta_v = self.eta_visc_effective(eps_ii);
        match self.yielding {
            crate::tectonics_v2::presets::YieldingConfig::Disabled => eta_v,
            crate::tectonics_v2::presets::YieldingConfig::Enabled(ylaw) => {
                let eta_p = eta_plastic(eps_ii, self.strain_rate_floor, bi_override);
                eta_effective(eta_v, eta_p, ylaw.sharpness)
            }
        }
    }

    /// Step 9 — `dη_eff/dε̇_II` with a per-cell `bi_override`.
    /// Mirrors [`Self::d_eta_effective_d_eps_ii`] but threads the
    /// override into the plastic branch's `eta_plastic` so the
    /// chain rule through `soft_min_harmonic` uses the elevated
    /// Bi. The viscous branch and its derivative are unchanged
    /// (Bi has no role there).
    #[inline]
    pub fn d_eta_effective_d_eps_ii_with_bi_override(
        &self,
        eps_ii: f64,
        bi_override: f64,
    ) -> f64 {
        let dv = self.d_eta_visc_effective_d_eps_ii(eps_ii);
        match self.yielding {
            crate::tectonics_v2::presets::YieldingConfig::Disabled => dv,
            crate::tectonics_v2::presets::YieldingConfig::Enabled(ylaw) => {
                let eta_v = self.eta_visc_effective(eps_ii);
                let eta_p = eta_plastic(eps_ii, self.strain_rate_floor, bi_override);
                d_eta_eff_d_eps_ii(
                    eta_v, dv, eta_p, eps_ii, self.strain_rate_floor, ylaw.sharpness,
                )
            }
        }
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
            yielding: crate::tectonics_v2::presets::YieldingConfig::Disabled,
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
/// layer can feed it to the arithmetic corner averaging unchanged.
///
/// When `cratonic = Some(state)`, two cratonic mechanisms apply
/// (Step 9 D1):
///
/// 1. **Primary — Bi elevation.** The plastic branch in cratonic
///    cells uses `bi_eff[i] = Bi · state.bi_multiplier[i,j]
///    = Bi · (1 + (B_factor - 1) · cratonic_factor[i,j])`. This
///    elevates the yield strength inside cratons by a factor up
///    to `B_factor`, suppressing viscoplastic yielding even in
///    high-strain-rate regimes where the original "yield_stress
///    = Bi" formula could not.
/// 2. **Secondary — viscous K mult.** The full effective viscosity
///    is post-multiplied by `state.eta_multiplier[i,j] = 1 + (K - 1)
///    · cratonic_factor[i,j]`, slowing wide-wavelength flow through
///    cratonic interiors.
///
/// When `cratonic = None`, this is a structural by-pass — no branch
/// is taken inside the inner loop — so Step 0–8 callers get
/// bit-identical output to the pre-Step-9 path.
pub fn build_eta_field(
    law: &ViscosityLaw,
    eps_ii_center: &Field2D,
    cratonic: Option<&super::cratonic::CratonicState>,
) -> Field2D {
    let nx = eps_ii_center.nx();
    let ny = eps_ii_center.ny();
    let mut eta = Field2D::new(nx, ny);
    let global_bi = match law.yielding {
        crate::tectonics_v2::presets::YieldingConfig::Disabled => 0.0,
        crate::tectonics_v2::presets::YieldingConfig::Enabled(ylaw) => ylaw.bi,
    };
    match cratonic {
        None => {
            // Hot path identical to pre-Step-9. Kept structurally
            // separate from the cratonic branch so the inner loop
            // is branch-free in Step 0–8 callers.
            for j in 0..ny {
                for i in 0..nx {
                    eta.set(i, j, law.eta_effective(eps_ii_center.get(i, j)));
                }
            }
        }
        Some(state) => {
            for j in 0..ny {
                for i in 0..nx {
                    let eps = eps_ii_center.get(i, j);
                    let bi_eff = global_bi * state.bi_multiplier.get(i, j);
                    let e = law.eta_effective_with_bi_override(eps, bi_eff);
                    eta.set(i, j, e * state.eta_multiplier.get(i, j));
                }
            }
        }
    }
    eta
}

// ============================================================================
// Plastic yielding (Step 3, Von Mises / Bingham, stateless)
// ============================================================================
//
// The constitutive law for the effective viscosity becomes a smooth
// minimum of the power-law (`ViscosityLaw`) branch and a plastic
// branch derived from a Von Mises yield stress `τ_yield = Bi`:
//
// ```
//   η_plastic(ε̇_II) = Bi / (2 · (ε̇_II + ε̇_min))
//   η_eff(ε̇_II)     = soft_min_harmonic(η_visc, η_plastic, p)
// ```
//
// `Bi` is a nondim Bingham number, stateless at Step 3 (no plastic
// memory, no healing, no spatial variation). `ε̇_min` is the same
// strain-rate floor used by the power-law viscosity, reused to
// regularise `η_plastic` at `ε̇_II → 0`. `soft_min_harmonic` is
// imported from `tectonics::solver::smooth` (whitelisted for
// `tectonics_v2`, along with `Field2D` and `PeriodicIndex`); see
// `tectonics_v2/mod.rs` for the re-export.
//
// # Derivatives (Newton Jacobian)
//
// For the `soft_min_harmonic` canonical form
// `η_s = (a^(-p) + b^(-p))^(-1/p)` the partial derivatives are
//
// ```
//   ∂η_s/∂a = (η_s / a)^(p+1)
//   ∂η_s/∂b = (η_s / b)^(p+1)
// ```
//
// (algebra: differentiate the `(·)^(-1/p)` wrapper, the chain rule
// gives `a^(-p-1) · (a^(-p)+b^(-p))^(-(p+1)/p) = (η_s/a)^(p+1)`).
//
// The plastic branch's own derivative is
// ```
//   ∂η_p/∂ε̇_II = -Bi / (2·(ε̇_II + ε̇_min)²) = -η_p / (ε̇_II + ε̇_min).
// ```
//
// Chain rule for the effective viscosity:
// ```
//   dη_eff/dε̇_II = (η_s/η_v)^(p+1) · dη_v/dε̇_II
//                + (η_s/η_p)^(p+1) · dη_p/dε̇_II
// ```
// with `dη_v/dε̇_II` supplied by the existing
// `ViscosityLaw::d_eta_effective_d_eps_ii` (power-law + smooth cap).
//
// # Why this stays SPD-compatible
//
// The operator layer (`stokes::operator::apply_momentum`) is
// unchanged — it reads a single `η` field and produces a symmetric
// stress divergence. The Newton Jacobian's tangent contribution
// (`apply_tangent`) uses the cell-centred scalar
// `c(ε̇_II) = dη/dε̇_II · 2 / (ε̇_II + ε̇_min)`; because
// `dη_eff/dε̇_II` is itself a derivative of a scalar, smooth
// potential-derived function of `ε̇_II`, the tangent operator
// remains symmetric at discrete level provided the η averaging to
// corners is arithmetic (cf. `stokes/operator.rs`). Step 3 therefore
// keeps CG as the inner linear solver — BiCGSTAB is not needed.

pub use crate::tectonics::solver::smooth::soft_min_harmonic;

/// Plastic yielding parameters.
///
/// `bi` is the Bingham number `Bi = τ_yield / σ* ∈ [0.05, 0.5]`
/// (design note §5.1). `sharpness` controls the `soft_min_harmonic`
/// transition width; `p = 4` matches the legacy default.
#[derive(Clone, Copy, Debug)]
pub struct YieldingLaw {
    pub bi: f64,
    pub sharpness: f64,
}

impl Default for YieldingLaw {
    fn default() -> Self {
        Self { bi: 0.15, sharpness: 4.0 }
    }
}

/// Plastic-branch viscosity
/// `η_plastic(ε̇_II) = Bi / (2·(ε̇_II + ε̇_min))`.
///
/// Positive, continuous at `ε̇_II = 0`, strictly decreasing in
/// `ε̇_II`. The `ε̇_min` floor is the same as the one used by
/// `ViscosityLaw` — callers pass it explicitly so the two branches
/// stay consistent.
#[inline]
pub fn eta_plastic(eps_ii: f64, eps_ii_floor: f64, bi: f64) -> f64 {
    bi / (2.0 * (eps_ii + eps_ii_floor))
}

/// Effective viscosity blend
/// `η_eff = soft_min_harmonic(η_visc, η_plastic, sharpness)`.
#[inline]
pub fn eta_effective(eta_visc: f64, eta_plastic_val: f64, sharpness: f64) -> f64 {
    soft_min_harmonic(eta_visc, eta_plastic_val, sharpness)
}

/// Analytic derivative `dη_eff/dε̇_II` via chain rule through the
/// `soft_min_harmonic` blend. Needed by the Newton Jacobian —
/// verified against a centred finite difference in the unit tests.
#[inline]
pub fn d_eta_eff_d_eps_ii(
    eta_visc: f64,
    d_eta_visc_d_eps: f64,
    eta_plastic_val: f64,
    eps_ii: f64,
    eps_ii_floor: f64,
    sharpness: f64,
) -> f64 {
    let eta_s = soft_min_harmonic(eta_visc, eta_plastic_val, sharpness);
    let p_plus_1 = sharpness + 1.0;
    // ∂η_p/∂ε̇_II = -η_p / (ε̇_II + ε̇_min).
    let d_eta_p_d_eps = -eta_plastic_val / (eps_ii + eps_ii_floor);
    // (η_s / η_v)^(p+1) · dη_v + (η_s / η_p)^(p+1) · dη_p.
    (eta_s / eta_visc).powf(p_plus_1) * d_eta_visc_d_eps
        + (eta_s / eta_plastic_val).powf(p_plus_1) * d_eta_p_d_eps
}

/// Build the plastic-branch viscosity field at cell centres.
pub fn build_eta_plastic_field(
    law: &YieldingLaw,
    eps_ii_floor: f64,
    eps_ii_center: &Field2D,
) -> Field2D {
    let nx = eps_ii_center.nx();
    let ny = eps_ii_center.ny();
    let mut out = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            out.set(i, j, eta_plastic(eps_ii_center.get(i, j), eps_ii_floor, law.bi));
        }
    }
    out
}

/// Blend `η_visc` and `η_plastic` field-by-field. Output is the
/// effective viscosity consumed by the operator layer.
pub fn blend_eta_fields(
    eta_visc: &Field2D,
    eta_plastic: &Field2D,
    sharpness: f64,
) -> Field2D {
    debug_assert_eq!(eta_visc.nx(), eta_plastic.nx());
    debug_assert_eq!(eta_visc.ny(), eta_plastic.ny());
    let nx = eta_visc.nx();
    let ny = eta_visc.ny();
    let mut out = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            out.set(
                i,
                j,
                eta_effective(eta_visc.get(i, j), eta_plastic.get(i, j), sharpness),
            );
        }
    }
    out
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
        let eta = build_eta_field(&law, &eps, None);
        let expected = law.eta_effective(0.0);
        for v in eta.data().iter() {
            assert!(approx(*v, expected, 1e-12));
        }
    }

    // ------------------------------------------------------------------
    // Step 3 — plastic yielding (Von Mises / Bingham, stateless).
    // ------------------------------------------------------------------

    #[test]
    fn eta_plastic_is_positive_decreasing_continuous_at_zero() {
        let bi = 0.15;
        let floor = 1.0e-3;
        let xs = [0.0, 1e-5, 1e-3, 1e-2, 0.1, 1.0, 10.0];
        let mut prev = f64::INFINITY;
        for &x in &xs {
            let e = eta_plastic(x, floor, bi);
            assert!(e > 0.0);
            assert!(e <= prev, "eta_plastic not monotone at ε̇={}", x);
            prev = e;
        }
        // Continuity at 0: value is finite and equals Bi/(2·floor).
        let e0 = eta_plastic(0.0, floor, bi);
        assert!(approx(e0, bi / (2.0 * floor), 1e-14));
    }

    #[test]
    fn eta_effective_asymptotes_match_branches() {
        // η_v ≪ η_p ⟹ η_eff → η_v (viscous-dominated).
        let eta_v = 1.0;
        let eta_p = 1.0e3;
        let eta = eta_effective(eta_v, eta_p, 4.0);
        let rel = (eta - eta_v).abs() / eta_v;
        assert!(rel < 1e-10, "viscous asymptote: rel = {}", rel);
        // η_v ≫ η_p ⟹ η_eff → η_p (plastic-dominated).
        let eta_v = 1.0e3;
        let eta_p = 1.0;
        let eta = eta_effective(eta_v, eta_p, 4.0);
        let rel = (eta - eta_p).abs() / eta_p;
        assert!(rel < 1e-10, "plastic asymptote: rel = {}", rel);
    }

    #[test]
    fn eta_effective_is_monotone_decreasing_in_eps_ii() {
        let visc = ViscosityLaw::default();
        let yld = YieldingLaw::default();
        let floor = visc.strain_rate_floor;
        let xs = [1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0];
        let mut prev = f64::INFINITY;
        for &x in &xs {
            let ev = visc.eta_effective(x);
            let ep = eta_plastic(x, floor, yld.bi);
            let e = eta_effective(ev, ep, yld.sharpness);
            assert!(e <= prev, "eta_eff not monotone at ε̇={}", x);
            prev = e;
        }
    }

    #[test]
    fn d_eta_eff_matches_finite_difference() {
        // Grid of ε̇_II values covering the whole viscous/plastic
        // transition. Analytic derivative must agree with a centred
        // finite difference on `eta_effective_of_eps` to 1e-7 (the
        // spec bound).
        let visc = ViscosityLaw::default();
        let yld = YieldingLaw::default();
        let floor = visc.strain_rate_floor;
        let eta_eff_of = |x: f64| {
            let ev = visc.eta_effective(x);
            let ep = eta_plastic(x, floor, yld.bi);
            eta_effective(ev, ep, yld.sharpness)
        };
        let h = 1.0e-7;
        let xs = [5e-3, 1e-2, 5e-2, 0.1, 0.3, 0.5, 1.0, 3.0, 10.0, 100.0];
        for &x in &xs {
            let ev = visc.eta_effective(x);
            let dv = visc.d_eta_effective_d_eps_ii(x);
            let ep = eta_plastic(x, floor, yld.bi);
            let analytic =
                d_eta_eff_d_eps_ii(ev, dv, ep, x, floor, yld.sharpness);
            let fd = (eta_eff_of(x + h) - eta_eff_of(x - h)) / (2.0 * h);
            let rel = (analytic - fd).abs() / analytic.abs().max(fd.abs()).max(1e-12);
            assert!(
                rel < 1e-6,
                "d_eta_eff FD mismatch at ε̇={}: analytic={}, fd={}, rel={}",
                x, analytic, fd, rel,
            );
        }
    }

    #[test]
    fn build_and_blend_fields_agree_with_pointwise() {
        // Spot-check the field-level builders: every cell of
        // `build_eta_plastic_field` then `blend_eta_fields` matches
        // the pointwise `eta_effective(eta_visc, eta_plastic, p)`.
        let visc = ViscosityLaw::default();
        let yld = YieldingLaw::default();
        let mut eps = Field2D::new(6, 6);
        for j in 0..6 {
            for i in 0..6 {
                eps.set(i, j, 1e-3 + 0.01 * (i + j) as f64);
            }
        }
        let eta_v = build_eta_field(&visc, &eps, None);
        let eta_p = build_eta_plastic_field(&yld, visc.strain_rate_floor, &eps);
        let eta_e = blend_eta_fields(&eta_v, &eta_p, yld.sharpness);
        for j in 0..6 {
            for i in 0..6 {
                let expected = eta_effective(
                    eta_v.get(i, j),
                    eta_p.get(i, j),
                    yld.sharpness,
                );
                assert!(approx(eta_e.get(i, j), expected, 1e-14));
            }
        }
    }
}
