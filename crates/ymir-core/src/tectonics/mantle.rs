//! Mantle convection proxy — large-scale flow field driving plate motion.
//!
//! Generates a smooth, divergence-free velocity field from a stream function
//! decomposed into low-frequency Fourier modes on the torus. This provides
//! the continuous energy source that keeps plates moving indefinitely.

use super::solver::field::Field2D;

/// Configuration for the mantle flow field.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct MantleConfig {
    /// Whether mantle flow is enabled. Default: true.
    pub enabled: bool,
    /// Number of Fourier modes in the stream function. More modes = more
    /// complex flow patterns. Default: 6.
    pub num_modes: usize,
    /// Amplitude of the mantle flow velocity (grid units / timestep). Default: 0.5.
    pub amplitude: f64,
    /// Coupling strength: how much the mantle drags the lithosphere.
    /// Thick crust couples more strongly. Default: 1.0.
    pub coupling: f64,
    /// Rate at which convection pattern evolves over time.
    /// 0.0 = static pattern, 0.01 = slow evolution. Default: 0.0.
    pub evolution_rate: f64,
}

impl Default for MantleConfig {
    fn default() -> Self {
        Self { enabled: true, num_modes: 6, amplitude: 0.5, coupling: 1.0, evolution_rate: 0.0 }
    }
}

/// A single Fourier mode of the stream function.
#[derive(Debug, Clone)]
struct MantleMode {
    kx: i32,
    ky: i32,
    amplitude: f64,
    phase: f64,
}

/// The mantle flow field.
#[derive(Clone)]
pub struct MantleFlow {
    modes: Vec<MantleMode>,
    pub vx: Field2D,
    pub vy: Field2D,
    n: usize,
}

impl MantleFlow {
    /// Generate a mantle flow field from a seed.
    ///
    /// Uses low-frequency Fourier modes (wave numbers 1-3) to create
    /// large-scale convection cells. The stream function ensures the
    /// resulting velocity field is divergence-free.
    pub fn generate(n: usize, seed: u64, config: &MantleConfig) -> Self {
        use std::f64::consts::TAU;

        // xorshift64-like RNG mapped to [-1, 1]
        let mut state = seed.wrapping_add(0x9E3779B97F4A7C15);
        let mut next_random = move || -> f64 {
            state ^= state >> 12;
            state ^= state << 25;
            state ^= state >> 27;
            let val = state.wrapping_mul(0x2545F4914F6CDD1D);
            (val as f64 / u64::MAX as f64) * 2.0 - 1.0
        };

        let mut modes = Vec::with_capacity(config.num_modes);
        let max_k = 3i32;

        for _ in 0..config.num_modes {
            let kx = ((next_random() * max_k as f64).round() as i32).clamp(-max_k, max_k);
            let ky = ((next_random() * max_k as f64).round() as i32).clamp(-max_k, max_k);

            if kx == 0 && ky == 0 {
                continue;
            }

            let amplitude = config.amplitude * next_random();
            let phase = next_random() * TAU;

            modes.push(MantleMode { kx, ky, amplitude, phase });
        }

        let mut flow = MantleFlow { modes, vx: Field2D::new(n), vy: Field2D::new(n), n };
        flow.recompute_field();
        flow
    }

    /// Recompute the cached velocity field from modes.
    ///
    /// Stream function: ψ(x,y) = Σ A_i × sin(kx_i × x + ky_i × y + φ_i)
    /// Velocity:        vx = ∂ψ/∂y,  vy = -∂ψ/∂x  (div-free)
    fn recompute_field(&mut self) {
        use std::f64::consts::TAU;
        let n = self.n;
        let nf = n as f64;

        for j in 0..n {
            for i in 0..n {
                let x = TAU * i as f64 / nf;
                let y = TAU * j as f64 / nf;

                let mut vx_sum = 0.0;
                let mut vy_sum = 0.0;

                for mode in &self.modes {
                    let kx = mode.kx as f64;
                    let ky = mode.ky as f64;
                    let arg = kx * x + ky * y + mode.phase;
                    let cos_arg = arg.cos();
                    vx_sum += mode.amplitude * ky * cos_arg;
                    vy_sum -= mode.amplitude * kx * cos_arg;
                }

                self.vx.set(i, j, vx_sum);
                self.vy.set(i, j, vy_sum);
            }
        }
    }

    /// Optionally evolve the convection pattern by perturbing mode phases.
    pub fn evolve(&mut self, evolution_rate: f64, _step: usize) {
        if evolution_rate <= 0.0 {
            return;
        }
        for (idx, mode) in self.modes.iter_mut().enumerate() {
            let rate = evolution_rate * (1.0 + 0.3 * (idx as f64).sin());
            mode.phase += rate;
        }
        self.recompute_field();
    }
}
