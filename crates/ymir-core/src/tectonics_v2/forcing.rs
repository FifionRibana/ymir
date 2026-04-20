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

/// Sinusoidal placeholder forcing driving a Kolmogorov-like shear
/// flow: `f̃ = ε · sin(2π ỹ / L̃y) · ê_x`.
///
/// **Deviation from prompt literal.** The Step 0 spec names
/// `ε·sin(2πx̃/L̃x)·ê_x`, but that force is a pure gradient
/// (`∇(-ε·cos(2πx̃/L̃x)/(2π/L̃x))`), so under incompressible Stokes
/// the momentum balance is satisfied with `ṽ = 0` and `p̃` absorbing
/// the force. The smoke test would then be trivial and the advection
/// coupling would not be exercised. The y-periodic variant is a
/// minimal, rotational force (`∇×f ≠ 0`) that still produces a simple
/// analytic flow and the prompt-specified "~3% S displacement over
/// 300 steps" behaviour.
///
/// Analytic reference: with η = 1 and domain `[0, L̃y]`, the steady
/// solution is `ṽx = ε · sin(2π ỹ / L̃y) / (2π/L̃y)²`, `ṽy = 0`,
/// `p̃ = 0`. Peak `|ṽ|` with `ε = 0.1`, `L̃y = 1` is `≈ 2.5e-3`.
#[derive(Clone, Copy, Debug)]
pub struct SinusoidalForce {
    pub amplitude: f64,
    pub ly: f64,
}

impl SinusoidalForce {
    pub fn new(amplitude: f64, ly: f64) -> Self {
        Self { amplitude, ly }
    }
}

impl Default for SinusoidalForce {
    fn default() -> Self {
        Self { amplitude: 0.1, ly: 1.0 }
    }
}

impl BodyForce for SinusoidalForce {
    fn fx(&self, _x: f64, y: f64, _s: f64) -> f64 {
        self.amplitude * (2.0 * PI * y / self.ly).sin()
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
