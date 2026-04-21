//! Body-force abstraction used as the Stokes momentum RHS.
//!
//! At Step 0 only a sinusoidal placeholder is driven; the trait is the
//! future integration point for GPE spreading (Step 2), slab pull
//! (Step 7), and mantle forcing (Step 8).
//!
//! The trait samples the force at face locations on the MAC grid (x-face
//! for the x-component, y-face for the y-component). Coordinates passed
//! to the trait are nondimensional, in [0, L̃x) × [0, L̃y).

use std::f64::consts::PI;

use super::field::Field2D;

/// A vector body force sampled at grid face locations.
pub trait BodyForce: Send + Sync {
    /// x-component sampled at the left vertical face of cell (i, j),
    /// coordinates ((i·dx), (j+0.5)·dx).
    fn fx(&self, x: f64, y: f64, s: f64) -> f64;
    /// y-component sampled at the bottom horizontal face of cell (i, j),
    /// coordinates ((i+0.5)·dx, (j·dx)).
    fn fy(&self, x: f64, y: f64, s: f64) -> f64;
}

/// Null body force. Used in tests where flow should decay to v = 0.
#[derive(Clone, Copy, Debug, Default)]
pub struct ZeroForce;

impl BodyForce for ZeroForce {
    fn fx(&self, _x: f64, _y: f64, _s: f64) -> f64 { 0.0 }
    fn fy(&self, _x: f64, _y: f64, _s: f64) -> f64 { 0.0 }
}

/// Sinusoidal placeholder forcing: `f̃ = ε · sin(2π x̃ / L̃x) · ê_x`.
///
/// In the thin-sheet formulation this force **produces flow** — there
/// is no pressure unknown available to balance it as a gradient, so
/// the momentum operator must deform the velocity field to cancel it.
///
/// # Amplitude calibration
///
/// Step 0 used `ε = 0.1`, giving `peak|ṽ| ≈ 1.3·10⁻³` and
/// `peak ε̇_II ≈ 8·10⁻³` — well above the default `ε̇_min = 10⁻³`
/// in Step 1's rheology, which matters because at Step 1 `ε̇_II` has
/// to clear that floor to activate the power-law regime. For
/// Step 1 the rheology-saturated effective viscosity pushes the
/// linear response down by a factor `η_eff(ε̇≈0) ≈ 100`, so an
/// equivalent Step-0 amplitude of 0.1 yields `peak|ṽ| ≈ 1.3·10⁻⁵`
/// and `peak ε̇ ≈ 10⁻⁴` — **below** `ε̇_min`. The floor dominates,
/// `η_max/η_min → 1`, and the power-law is never exercised.
///
/// To keep the placeholder diagnostic at Step 1 (`η_max/η_min` in
/// the single-digit-to-~30 range, Newton genuinely iterating), we
/// raise the default amplitude to `ε = 10.0`. That restores the
/// Step-0-equivalent `ε̇` range once the ~100× viscosity drag from
/// the floor is applied. The placeholder disappears at Step 2 when
/// GPE provides the physical driving term; amplitude is a
/// step-local knob, not a milestone-wide parameter.
///
/// (Note: in an incompressible-Stokes reading, the same force is a
/// pure gradient and produces `ṽ = 0` — one of the reasons the
/// thin-sheet vs incompressible-Stokes architecture distinction
/// matters.)
#[derive(Clone, Copy, Debug)]
pub struct SinusoidalForce {
    pub amplitude: f64,
    pub lx: f64,
}

impl SinusoidalForce {
    pub fn new(amplitude: f64, lx: f64) -> Self {
        Self { amplitude, lx }
    }
}

impl Default for SinusoidalForce {
    fn default() -> Self {
        Self { amplitude: 10.0, lx: 1.0 }
    }
}

impl BodyForce for SinusoidalForce {
    fn fx(&self, x: f64, _y: f64, _s: f64) -> f64 {
        self.amplitude * (2.0 * PI * x / self.lx).sin()
    }
    fn fy(&self, _x: f64, _y: f64, _s: f64) -> f64 {
        0.0
    }
}

/// Sample a `BodyForce` onto MAC-staggered face-located fields.
///
/// `s` is the crustal thickness at cell centers and may be consulted by
/// state-dependent forces (no-op for the Step 0 placeholders).
pub fn sample_to_faces<F: BodyForce>(
    force: &F,
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    s: &Field2D,
    fx: &mut Field2D,
    fy: &mut Field2D,
) {
    for j in 0..ny {
        for i in 0..nx {
            // vx lives at (i·dx, (j+0.5)·dy); S for this face interpolates
            // between cells (i-1, j) and (i, j).
            let xfx = (i as f64) * dx;
            let yfx = (j as f64 + 0.5) * dy;
            let i_left = if i == 0 { nx - 1 } else { i - 1 };
            let s_fx = 0.5 * (s.get(i_left, j) + s.get(i, j));
            fx.set(i, j, force.fx(xfx, yfx, s_fx));

            // vy lives at ((i+0.5)·dx, j·dy); S interpolates (i, j-1) and (i, j).
            let xfy = (i as f64 + 0.5) * dx;
            let yfy = (j as f64) * dy;
            let j_below = if j == 0 { ny - 1 } else { j - 1 };
            let s_fy = 0.5 * (s.get(i, j_below) + s.get(i, j));
            fy.set(i, j, force.fy(xfy, yfy, s_fy));
        }
    }
}
