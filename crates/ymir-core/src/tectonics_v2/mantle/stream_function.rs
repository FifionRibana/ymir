//! Fourier stream-function generation for the mantle pattern.
//!
//! Generates `ψ(x, y) = Σ_k a_k · sin(kx_k · 2π x + φx_k) · sin(ky_k · 2π y + φy_k)`
//! sampled at **grid nodes** `(i · dx, j · dy)` for
//! `i = 0..nx, j = 0..ny` (periodic). The separable product form
//! matches D2 of the Step 8 spec; it preserves div-freeness when
//! combined with the staggered curl in `pattern.rs`.
//!
//! Nodal semantics
//! ----------------
//! The returned `Field2D` has shape `nx × ny` and is the standard
//! row-major container used throughout `tectonics_v2`, **but**
//! entry `[i, j]` represents `ψ(i · dx, j · dy)` — a grid corner,
//! **not** a cell centre. This is a deliberate semantic choice: the
//! discrete curl on the MAC-staggered velocity grid is exactly
//! divergence-free (`div v ≡ 0` by algebraic cancellation) only
//! when `ψ` is nodal. A cell-centered `ψ` would leave O(dx²)
//! residual divergence and violate the Step 8 acceptance
//! `div_v_mantle_max < 10⁻¹⁰` (at 256², dx² ≈ 1.5e-5 ≫ 1e-10).
//!
//! The two interpretations share the same storage layout; only
//! the spatial mapping differs. The `pattern.rs` staggered curl
//! consumes this nodal view directly.
//!
//! Normalisation
//! -------------
//! After sampling, `ψ` is rescaled so `max|ψ| = 1` (or left at
//! zero if all modes cancel). The final velocity magnitude is
//! then controlled by `Mf`: `v_mantle = Mf · curl(ψ)`.

use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use super::super::field::Field2D;

/// Knobs for [`generate_stream_function`].
#[derive(Clone, Copy, Debug)]
pub struct StreamFunctionConfig {
    /// Number of Fourier modes summed into `ψ`. Default `6`.
    pub num_modes: usize,
    /// Seed for mode generation (wave numbers, amplitudes,
    /// phases).
    pub seed: u64,
}

impl Default for StreamFunctionConfig {
    fn default() -> Self {
        Self { num_modes: super::NUM_MODES_DEFAULT, seed: 42 }
    }
}

/// One Fourier mode of the separable stream function.
#[derive(Clone, Copy, Debug)]
struct Mode {
    kx: u32,
    ky: u32,
    amplitude: f64,
    phase_x: f64,
    phase_y: f64,
}

/// Generate a nodal stream function field.
///
/// Samples `ψ(i · dx, j · dy)` at grid nodes; returns a
/// `Field2D` sized `nx × ny` with the nodal interpretation
/// documented in the module comment. The field is normalised so
/// `max|ψ| = 1` (or left zero if all modes cancelled).
///
/// Wave-number convention: `kx, ky ∈ {1, 2, 3}` drawn uniformly
/// per D2. Phases uniform on `[0, 2π)`. Amplitudes drawn from a
/// standard normal (Box–Muller via two uniforms), giving varied
/// Fourier weights while keeping the normalised output in `[-1,
/// 1]`.
pub fn generate_stream_function(
    nx: usize,
    ny: usize,
    config: &StreamFunctionConfig,
) -> Field2D {
    let mut rng = ChaCha8Rng::seed_from_u64(config.seed);
    let modes = draw_modes(&mut rng, config.num_modes);
    sample_nodal(nx, ny, &modes)
}

fn draw_modes(rng: &mut ChaCha8Rng, num_modes: usize) -> Vec<Mode> {
    use std::f64::consts::TAU;
    let mut modes = Vec::with_capacity(num_modes);
    for _ in 0..num_modes {
        // Wave numbers ∈ {1, 2, 3}. Convert from a uniform [0,1)
        // via `random::<f64>()` rather than `gen_range`, matching
        // the rand API used elsewhere in tectonics_v2 (voronoi.rs).
        let kx = 1 + ((rng.random::<f64>() * 3.0).floor() as u32).min(2);
        let ky = 1 + ((rng.random::<f64>() * 3.0).floor() as u32).min(2);
        // Box–Muller on two uniforms for a Gaussian amplitude so
        // mode strength varies naturally — some dominant, some
        // weak, emergent from the seed.
        let u1: f64 = rng.random::<f64>().max(1.0e-12);
        let u2: f64 = rng.random::<f64>();
        let amplitude = (-2.0 * u1.ln()).sqrt() * (TAU * u2).cos();
        let phase_x = rng.random::<f64>() * TAU;
        let phase_y = rng.random::<f64>() * TAU;
        modes.push(Mode { kx, ky, amplitude, phase_x, phase_y });
    }
    modes
}

fn sample_nodal(nx: usize, ny: usize, modes: &[Mode]) -> Field2D {
    use std::f64::consts::TAU;
    let mut psi = Field2D::new(nx, ny);
    let nxf = nx as f64;
    let nyf = ny as f64;
    let mut max_abs = 0.0_f64;
    for j in 0..ny {
        let y = j as f64 / nyf;
        for i in 0..nx {
            let x = i as f64 / nxf;
            let mut sum = 0.0_f64;
            for m in modes {
                let sx = (m.kx as f64 * TAU * x + m.phase_x).sin();
                let sy = (m.ky as f64 * TAU * y + m.phase_y).sin();
                sum += m.amplitude * sx * sy;
            }
            psi.set(i, j, sum);
            let abs = sum.abs();
            if abs > max_abs {
                max_abs = abs;
            }
        }
    }
    if max_abs > 0.0 {
        let inv = 1.0 / max_abs;
        for v in psi.data_mut().iter_mut() {
            *v *= inv;
        }
    }
    psi
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_seed_gives_same_field() {
        let cfg = StreamFunctionConfig { num_modes: 6, seed: 42 };
        let a = generate_stream_function(32, 32, &cfg);
        let b = generate_stream_function(32, 32, &cfg);
        for (va, vb) in a.data().iter().zip(b.data().iter()) {
            assert_eq!(va, vb);
        }
    }

    #[test]
    fn different_seeds_give_different_fields() {
        let a = generate_stream_function(32, 32, &StreamFunctionConfig { num_modes: 6, seed: 42 });
        let b = generate_stream_function(32, 32, &StreamFunctionConfig { num_modes: 6, seed: 43 });
        let mut differ = false;
        for (va, vb) in a.data().iter().zip(b.data().iter()) {
            if (va - vb).abs() > 1e-12 { differ = true; break; }
        }
        assert!(differ, "seeds 42 vs 43 produced identical fields");
    }

    /// After normalisation, `max|ψ| = 1` (or exactly 0 if every
    /// mode cancelled, which is degenerate).
    #[test]
    fn normalisation_hits_unit_max() {
        let psi = generate_stream_function(
            64, 64,
            &StreamFunctionConfig { num_modes: 6, seed: 42 },
        );
        let max = psi.data().iter().cloned().map(f64::abs).fold(0.0_f64, f64::max);
        assert!((max - 1.0).abs() < 1e-12, "max |ψ| = {}, expected 1", max);
    }

    /// On a periodic grid sampled at exactly nx points per
    /// period, the mean of a sum of pure-sine modes is **not**
    /// generally zero after normalisation (the normalisation
    /// rescales the field based on its L∞ norm, and a rescaled
    /// sum of sines has mean equal to (rescale factor) × 0 = 0
    /// analytically, but discretely the sine samples on a
    /// periodic grid also sum to 0 exactly via periodicity).
    /// Verify mean ≈ 0.
    #[test]
    fn mean_is_zero_on_periodic_grid() {
        let psi = generate_stream_function(
            32, 32,
            &StreamFunctionConfig { num_modes: 6, seed: 42 },
        );
        let mean: f64 = psi.data().iter().sum::<f64>() / (32 * 32) as f64;
        assert!(mean.abs() < 1e-12, "mean ψ = {}, expected ≈ 0", mean);
    }
}
