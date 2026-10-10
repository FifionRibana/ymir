//! Seeded noise generator with anisotropic FBM support.

use noise::{NoiseFn, OpenSimplex};

/// A seeded noise generator that supports anisotropic sampling.
pub struct SeededNoise {
    /// One noise source per octave (different seeds for independence).
    sources: Vec<OpenSimplex>,
}

impl SeededNoise {
    /// Create a noise generator with `n_octaves` independent sources.
    pub fn new(base_seed: u32, n_octaves: usize) -> Self {
        let sources = (0..n_octaves)
            .map(|i| OpenSimplex::new(base_seed.wrapping_add(i as u32 * 7919)))
            .collect();
        Self { sources }
    }

    /// Sample a single octave at (x, y).
    pub fn sample(&self, octave: usize, x: f64, y: f64) -> f64 {
        self.sources[octave].get([x, y])
    }

    /// Sample FBM (sum of octaves) at (x, y) with the given parameters.
    /// Returns a value in approximately [-1, 1].
    pub fn fbm(&self, x: f64, y: f64, octaves: usize, lacunarity: f64, persistence: f64) -> f64 {
        let mut value = 0.0;
        let mut amplitude = 1.0;
        let mut frequency = 1.0;
        let mut max_amplitude = 0.0;

        for i in 0..octaves.min(self.sources.len()) {
            value += self.sample(i, x * frequency, y * frequency) * amplitude;
            max_amplitude += amplitude;
            frequency *= lacunarity;
            amplitude *= persistence;
        }

        value / max_amplitude
    }

    /// ADR Finding 169-A1 -- [`Self::fbm`] made PERIODIC on a torus of `period` × `period` (in the units of `x`, `y`,
    /// i.e. wavelengths): each octave samples the same OpenSimplex source on a 4D Clifford torus, (x, y) ↦
    /// (R cos 2πx/P, R sin 2πx/P, R cos 2πy/P, R sin 2πy/P) with R = P·f/2π, so its local wavelength stays 1/f and the
    /// field wraps exactly at P.
    pub fn fbm_torus(&self, x: f64, y: f64, period: f64, octaves: usize, lacunarity: f64, persistence: f64) -> f64 {
        let tau = std::f64::consts::TAU;
        let (ax, ay) = (tau * x / period, tau * y / period);
        let (cx, sx, cy, sy) = (ax.cos(), ax.sin(), ay.cos(), ay.sin());
        let mut value = 0.0;
        let mut amplitude = 1.0;
        let mut frequency = 1.0;
        let mut max_amplitude = 0.0;
        for i in 0..octaves.min(self.sources.len()) {
            let r = period * frequency / tau;
            value += self.sources[i].get([r * cx, r * sx, r * cy, r * sy]) * amplitude;
            max_amplitude += amplitude;
            frequency *= lacunarity;
            amplitude *= persistence;
        }
        value / max_amplitude
    }

    /// Sample anisotropic FBM: the noise is stretched along a direction.
    ///
    /// `angle` is the direction of compression (radians, 0 = x-axis).
    /// `ratio` is the anisotropy ratio: 1.0 = isotropic, 3.0 = stretched
    /// 3x perpendicular to `angle` (creating elongated features
    /// perpendicular to the compression direction).
    #[allow(clippy::too_many_arguments)]
    pub fn fbm_anisotropic(
        &self,
        x: f64,
        y: f64,
        angle: f64,
        ratio: f64,
        octaves: usize,
        lacunarity: f64,
        persistence: f64,
    ) -> f64 {
        // Rotate coordinates so that the compression direction aligns with x
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        let rx = x * cos_a + y * sin_a;
        let ry = -x * sin_a + y * cos_a;

        // Compress along the rotated x (slope direction) by the ratio
        let sx = rx * ratio;
        let sy = ry;

        self.fbm(sx, sy, octaves, lacunarity, persistence)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn noise_fbm_is_zero_mean() {
        let noise = SeededNoise::new(42, 7);
        let n = 1000;
        let mut sum = 0.0;
        for i in 0..n {
            for j in 0..n {
                sum += noise.fbm(i as f64 * 0.01, j as f64 * 0.01, 7, 2.0, 0.5);
            }
        }
        let mean = sum / (n * n) as f64;
        assert!(mean.abs() < 0.05, "FBM should be approximately zero-mean: {mean}");
    }
}
