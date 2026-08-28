//! Sinusoidal placeholder forcing.
//!
//! `f̃ = ε · sin(2π x̃ / L̃x) · ê_x`. See
//! [`super::body_force::BodyForce`] for the accumulation contract.
//!
//! The amplitude default is `ε = 10.0` (Step-1 regression setting).
//! `ε = 0.1` produces `peak ε̇` below `ε̇_min` when combined with
//! Step 1's power-law rheology, silently reducing the solver to a
//! linear regime — see the commit history of `tectonics_v2`.

use std::f64::consts::PI;

use super::body_force::{BodyForce, SimulationState, VectorField};

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
    fn accumulate(&self, state: &SimulationState, out: &mut VectorField) {
        let nx = state.nx;
        let ny = state.ny;
        let dx = state.dx;
        let lx = self.lx;
        let amp = self.amplitude;
        for j in 0..ny {
            for i in 0..nx {
                // vx face at (i·dx, (j+0.5)·dy).
                let x = i as f64 * dx;
                let k = j * nx + i;
                out.fx.data_mut()[k] += amp * (2.0 * PI * x / lx).sin();
            }
        }
        // fy = 0: nothing to add.
    }

    fn name(&self) -> &'static str {
        "SinusoidalForce"
    }
}

/// Null body force. Useful as a test baseline and as a placeholder
/// in `ForceSum` slots that are currently disabled.
#[derive(Clone, Copy, Debug, Default)]
pub struct ZeroForce;

impl BodyForce for ZeroForce {
    fn accumulate(&self, _state: &SimulationState, _out: &mut VectorField) {}
    fn name(&self) -> &'static str {
        "ZeroForce"
    }
}
