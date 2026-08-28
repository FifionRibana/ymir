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
//!
//! Time evolution (Step 12 R6 — Phys.A phase drift)
//! ------------------------------------------------
//! At Step 8 the pattern was static (D6, "Out of scope at Step 8").
//! Step 12 R6 wires `MantleConfig::Enabled.evolution_rate` by drifting
//! all modes' phases linearly with non-dimensional time:
//! ```
//!   φx_k(t) = φx_k(0) + ω · t
//!   φy_k(t) = φy_k(0) + ω · t       with  ω = evolution_rate · TAU
//! ```
//! Wave numbers and amplitudes are frozen at init. The normalisation
//! `init_norm = max|ψ(t=0)|` is also frozen, so `|ψ(t)|` may drift
//! around 1 instead of being clamped each step — this avoids jitter
//! on the argmax position. See [`StreamFunctionBuilder`].
//!
//! Div-freeness of `v_mantle = curl(ψ)` survives the drift by
//! construction: the staggered-curl cancellation in `pattern.rs` is
//! algebraic in the four corner values of ψ, independent of what
//! phases were used to compute them.

use std::f64::consts::TAU;

use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use crate::tectonics_v2::field::Field2D;

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

/// Builder for nodal stream functions, with optional Step 12 R6
/// phase-drift evolution.
///
/// The builder draws and freezes the Fourier modes (wave numbers,
/// amplitudes, base phases) and the t=0 normalisation factor
/// `init_norm = max|ψ(t=0)|`. Subsequent samples at non-dimensional
/// time `t > 0` apply a phase offset `ω · t` (with
/// `ω = evolution_rate · TAU`) to every mode's `(φx, φy)` and
/// divide by the *frozen* `init_norm`. The amplitude `max|ψ(t)|`
/// is therefore allowed to drift around 1 (no runtime renorm) —
/// preferred over per-step renormalisation because the argmax
/// position would otherwise jitter step-to-step.
///
/// For Step 8 callers that only need the static t=0 pattern,
/// [`generate_stream_function`] is a thin wrapper that builds and
/// samples in one call — and produces output bit-identical to the
/// pre-R6 implementation.
#[derive(Clone, Debug)]
pub struct StreamFunctionBuilder {
    base_modes: Vec<Mode>,
    init_norm: f64,
}

impl StreamFunctionBuilder {
    /// Draw the base modes from `config.seed` and compute the t=0
    /// normalisation factor. The grid size enters via the discrete
    /// `max|ψ|` sample (different `(nx, ny)` may select different
    /// argmax cells).
    pub fn new(nx: usize, ny: usize, config: &StreamFunctionConfig) -> Self {
        let mut rng = ChaCha8Rng::seed_from_u64(config.seed);
        let base_modes = draw_modes(&mut rng, config.num_modes);
        let psi_t0 = sample_nodal_unscaled(nx, ny, &base_modes);
        let raw_max = psi_t0.data().iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        let init_norm = if raw_max > 0.0 { raw_max } else { 1.0 };
        Self { base_modes, init_norm }
    }

    /// Sample ψ at non-dimensional time `t_nondim` with phase-drift
    /// rate `evolution_rate`. `t_nondim = 0` or `evolution_rate = 0`
    /// returns the static t=0 pattern bit-identical to
    /// [`generate_stream_function`].
    pub fn sample_at_time(
        &self,
        nx: usize,
        ny: usize,
        t_nondim: f64,
        evolution_rate: f64,
    ) -> Field2D {
        let phase_offset = evolution_rate * TAU * t_nondim;
        let mut psi = if phase_offset == 0.0 {
            sample_nodal_unscaled(nx, ny, &self.base_modes)
        } else {
            let drifted = drift_modes(&self.base_modes, phase_offset);
            sample_nodal_unscaled(nx, ny, &drifted)
        };
        let inv = 1.0 / self.init_norm;
        for v in psi.data_mut().iter_mut() {
            *v *= inv;
        }
        psi
    }

    /// Frozen normalisation factor `max|ψ(t=0)|` (or `1.0` if all
    /// modes degenerated to zero). Exposed for diagnostics.
    pub fn init_norm(&self) -> f64 {
        self.init_norm
    }
}

/// Generate a nodal stream function field at t=0.
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
///
/// Output is bit-identical to the pre-R6 implementation. For
/// multi-step rebuilds with phase drift, use
/// [`StreamFunctionBuilder`] directly.
pub fn generate_stream_function(nx: usize, ny: usize, config: &StreamFunctionConfig) -> Field2D {
    StreamFunctionBuilder::new(nx, ny, config).sample_at_time(nx, ny, 0.0, 0.0)
}

/// Generate a nodal stream function field at non-dimensional time
/// `t_nondim` with phase-drift rate `evolution_rate` (Step 12 R6).
///
/// One-shot convenience equivalent to building a
/// [`StreamFunctionBuilder`] and calling
/// [`StreamFunctionBuilder::sample_at_time`] once. For per-step
/// rebuilds in the harness loop, prefer the builder so the modes
/// and `init_norm` are drawn only once.
pub fn generate_stream_function_at_time(
    nx: usize,
    ny: usize,
    config: &StreamFunctionConfig,
    t_nondim: f64,
    evolution_rate: f64,
) -> Field2D {
    StreamFunctionBuilder::new(nx, ny, config).sample_at_time(nx, ny, t_nondim, evolution_rate)
}

fn draw_modes(rng: &mut ChaCha8Rng, num_modes: usize) -> Vec<Mode> {
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

/// Apply a uniform phase offset to every mode's `(phase_x, phase_y)`.
/// Wave numbers and amplitudes are passed through unchanged.
fn drift_modes(base: &[Mode], phase_offset: f64) -> Vec<Mode> {
    base.iter()
        .map(|m| Mode {
            kx: m.kx,
            ky: m.ky,
            amplitude: m.amplitude,
            phase_x: m.phase_x + phase_offset,
            phase_y: m.phase_y + phase_offset,
        })
        .collect()
}

/// Sample ψ at grid nodes without normalisation. Caller divides by
/// the desired scaling factor afterwards (`init_norm` for the
/// builder, or `max|ψ|` for the legacy renormalise-each-call path).
fn sample_nodal_unscaled(nx: usize, ny: usize, modes: &[Mode]) -> Field2D {
    let mut psi = Field2D::new(nx, ny);
    let nxf = nx as f64;
    let nyf = ny as f64;
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
            if (va - vb).abs() > 1e-12 {
                differ = true;
                break;
            }
        }
        assert!(differ, "seeds 42 vs 43 produced identical fields");
    }

    /// After normalisation, `max|ψ| = 1` (or exactly 0 if every
    /// mode cancelled, which is degenerate).
    #[test]
    fn normalisation_hits_unit_max() {
        let psi =
            generate_stream_function(64, 64, &StreamFunctionConfig { num_modes: 6, seed: 42 });
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
        let psi =
            generate_stream_function(32, 32, &StreamFunctionConfig { num_modes: 6, seed: 42 });
        let mean: f64 = psi.data().iter().sum::<f64>() / (32 * 32) as f64;
        assert!(mean.abs() < 1e-12, "mean ψ = {}, expected ≈ 0", mean);
    }
}
