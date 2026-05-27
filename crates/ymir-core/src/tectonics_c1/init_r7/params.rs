//! Tunables for R7 boundary displacement (Phase 2 Track B
//! sub-component 1, Issue #131).
//!
//! See [`super`] for the module-level rationale and the
//! [`boundary_displacement`](super::boundary_displacement) module
//! for the algorithm.

/// Phase 2 Track B R7 boundary displacement parameters.
///
/// Defaults calibrated for the `grid_size = 64`, 8-plate Voronoï
/// configuration produced by Phase 1.1 init. For other grid sizes
/// the caller should typically rescale `amplitude` proportionally
/// (e.g. `amplitude = grid_size as f64 / 8.0`) — the noise
/// frequency stays grid-relative because the apply function
/// normalises the sample position into `[0, 1)` before passing to
/// the FBM stack.
///
/// ## Default rationale
///
/// - `amplitude = 8.0` — magnitude of the per-cell displacement
///   vector in cell units. At `grid_size = 64` this is 1/8 of the
///   domain, producing a typical boundary deviation of `~1–3`
///   cells (since FBM output is approximately `[−1, +1]` modulo
///   octave normalisation). Larger amplitudes produce more curved
///   but less coherent boundaries; Stage E1 unit test 5 enforces
///   that the per-cell reassignment count stays in
///   `[1, 20 %]` of total cells.
/// - `frequency = 4.0` — noise frequency multiplier. The apply
///   function samples noise at `(i / nx, j / ny) × frequency`, so
///   `frequency = 4.0` means roughly 4 wavelengths across the
///   domain.
/// - `octaves = 3` — fractal Brownian motion octave count. More
///   octaves produce finer-scale boundary roughness.
/// - `persistence = 0.5` — amplitude ratio between successive
///   octaves. Standard FBM value.
/// - `seed = 0` — caller should typically pass the same seed
///   used for `generate_voronoi` so the same input
///   `(grid_size, seed)` produces a bit-identical Phase 2 R7
///   init.
///
/// `enabled` follows the same W4 closure-isolation discipline as
/// Phase 1.2 / 1.3 / 1.4 / Track A closures: when `false`,
/// [`super::boundary_displacement::apply_boundary_displacement`]
/// is a no-op and the run reproduces the pre-displacement Voronoï
/// state bit-identically.
#[derive(Clone, Copy, Debug)]
pub struct R7InitParams {
    /// Master enable/disable. When `false`,
    /// [`super::boundary_displacement::apply_boundary_displacement`]
    /// is a no-op (W4 closure-isolation discipline).
    pub enabled: bool,
    /// Displacement amplitude in cell units. For other grid sizes,
    /// rescale via `amplitude = grid_size as f64 / 8.0`.
    pub amplitude: f64,
    /// Noise frequency multiplier (wavelengths per domain).
    pub frequency: f64,
    /// FBM octave count.
    pub octaves: u32,
    /// FBM amplitude ratio between successive octaves.
    pub persistence: f64,
    /// Deterministic seed. Same `(grid_size, seed)` → same
    /// displaced `plate_id` byte-for-byte. The apply function
    /// derives two independent `u32` channels from this `u64`
    /// (one per displacement component).
    pub seed: u64,
}

impl Default for R7InitParams {
    fn default() -> Self {
        Self {
            enabled: true,
            amplitude: 8.0,
            frequency: 4.0,
            octaves: 3,
            persistence: 0.5,
            seed: 0,
        }
    }
}
