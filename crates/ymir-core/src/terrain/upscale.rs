//! Anisotropic FBM upscaling: transforms a coarse isostatic heightmap into
//! detailed terrain by adding fractal noise modulated by local slope.

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use super::noise::SeededNoise;
use crate::grid::GridF32;
use crate::seed::WorldSeed;

/// Configuration for anisotropic FBM upscaling.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FbmUpscaleConfig {
    /// Target resolution (side length in pixels). Default: 1024.
    pub target_size: usize,
    /// Number of FBM octaves. More = finer detail. Default: 7.
    pub octaves: usize,
    /// Frequency multiplier per octave. Default: 2.0.
    pub lacunarity: f64,
    /// Amplitude multiplier per octave. Default: 0.5.
    pub persistence: f64,
    /// Base noise amplitude as fraction of altitude range. Default: 0.08.
    pub amplitude_base: f64,
    /// Extra amplitude scaling on steep slopes. Default: 3.0.
    pub amplitude_slope_factor: f64,
    /// Maximum anisotropy ratio on steep slopes. Default: 3.0.
    pub max_anisotropy: f64,
    /// Amplitude reduction below sea level. Default: 0.3.
    pub submarine_damping: f64,
    /// Base frequency of the first octave, in cycles per source pixel. Default: 1.0.
    pub base_frequency: f64,
}

impl Default for FbmUpscaleConfig {
    fn default() -> Self {
        Self {
            target_size: 1024,
            octaves: 7,
            lacunarity: 2.0,
            persistence: 0.5,
            amplitude_base: 0.08,
            amplitude_slope_factor: 3.0,
            max_anisotropy: 3.0,
            submarine_damping: 0.3,
            base_frequency: 1.0,
        }
    }
}

/// Result of FBM upscaling.
pub struct UpscaleResult {
    /// The upscaled heightmap at target resolution.
    pub heightmap: GridF32,
    /// The slope magnitude field (useful for erosion and viz).
    pub slope: GridF32,
}

/// Upscale a coarse heightmap using anisotropic FBM.
///
/// The coarse heightmap (from isostasy, typically 64-128) is bilinearly
/// interpolated to the target resolution, then fractal noise is added with
/// amplitude and direction controlled by the local terrain properties.
///
/// `sea_level` is the normalized value [0,1] that separates ocean from land.
pub fn upscale_with_fbm(
    coarse: &GridF32,
    sea_level: f32,
    seed: &WorldSeed,
    config: &FbmUpscaleConfig,
) -> UpscaleResult {
    let src_w = coarse.width;
    let src_h = coarse.height;
    let dst = config.target_size;
    let scale_x = (src_w - 1) as f64 / (dst - 1) as f64;
    let scale_y = (src_h - 1) as f64 / (dst - 1) as f64;

    // Create noise generator
    let noise_seed = seed.derive_seed("fbm_upscale") as u32;
    let noise = SeededNoise::new(noise_seed, config.octaves);

    // Precompute slope and direction on the coarse grid
    let (slope_map, direction_map) = compute_terrain_analysis(coarse);

    // Base frequency: one cycle per source cell mapped to target coords
    let freq = config.base_frequency / (src_w as f64);

    // Process each output row in parallel
    let row_data: Vec<(Vec<f32>, Vec<f32>)> = (0..dst)
        .into_par_iter()
        .map(|j| {
            let mut h_row = vec![0.0f32; dst];
            let mut s_row = vec![0.0f32; dst];

            for i in 0..dst {
                // Source coordinates in coarse pixel space
                let sx = i as f64 * scale_x;
                let sy = j as f64 * scale_y;

                // 1. Bilinear interpolation of the coarse heightmap
                let base_height = coarse.sample_bilinear(sx as f32, sy as f32);

                // 2. Sample terrain properties from coarse analysis
                let slope_mag = slope_map.sample_bilinear(sx as f32, sy as f32);
                let slope_dir = direction_map.sample_bilinear(sx as f32, sy as f32);

                // 3. Compute amplitude modulation
                let altitude_factor =
                    if base_height > sea_level { 1.0 } else { config.submarine_damping };

                let amplitude = config.amplitude_base
                    * (1.0 + slope_mag as f64 * config.amplitude_slope_factor)
                    * altitude_factor;

                // 4. Compute anisotropy ratio (proportional to slope)
                let anisotropy = 1.0 + (config.max_anisotropy - 1.0) * (slope_mag as f64).min(1.0);

                // 5. Sample anisotropic FBM
                let nx = i as f64 * freq;
                let ny = j as f64 * freq;

                let noise_val = if anisotropy > 1.01 {
                    noise.fbm_anisotropic(
                        nx,
                        ny,
                        slope_dir as f64,
                        anisotropy,
                        config.octaves,
                        config.lacunarity,
                        config.persistence,
                    )
                } else {
                    noise.fbm(nx, ny, config.octaves, config.lacunarity, config.persistence)
                };

                // 6. Combine: base + noise
                let final_height = (base_height as f64 + amplitude * noise_val).clamp(0.0, 1.0);

                h_row[i] = final_height as f32;
                s_row[i] = slope_mag;
            }

            (h_row, s_row)
        })
        .collect();

    // Copy into GridF32
    let mut heightmap = GridF32::new(dst, dst, 0.0);
    let mut slope_out = GridF32::new(dst, dst, 0.0);

    for (j, (h_row, s_row)) in row_data.into_iter().enumerate() {
        for (i, (h, s)) in h_row.into_iter().zip(s_row).enumerate() {
            heightmap.set(i, j, h);
            slope_out.set(i, j, s);
        }
    }

    UpscaleResult { heightmap, slope: slope_out }
}

/// Compute slope magnitude and direction on the coarse grid.
fn compute_terrain_analysis(heightmap: &GridF32) -> (GridF32, GridF32) {
    let w = heightmap.width;
    let h = heightmap.height;
    let mut slope = GridF32::new(w, h, 0.0);
    let mut direction = GridF32::new(w, h, 0.0);

    for j in 0..h {
        for i in 0..w {
            let (gx, gy) = heightmap.gradient_at(i, j);
            let mag = (gx * gx + gy * gy).sqrt();
            let dir = gy.atan2(gx);
            slope.set(i, j, mag);
            direction.set(i, j, dir);
        }
    }

    (slope, direction)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn upscale_preserves_mean_altitude() {
        let coarse = GridF32::new(16, 16, 0.5);
        let seed = WorldSeed::new(42);
        let config = FbmUpscaleConfig { target_size: 128, ..Default::default() };

        let result = upscale_with_fbm(&coarse, 0.1, &seed, &config);
        let mean = result.heightmap.mean();

        assert!((mean - 0.5).abs() < 0.05, "Mean should be near 0.5, got {mean}");
    }

    #[test]
    fn upscale_is_deterministic() {
        let coarse = GridF32::new(16, 16, 0.5);
        let seed = WorldSeed::new(42);
        let config = FbmUpscaleConfig { target_size: 64, ..Default::default() };

        let r1 = upscale_with_fbm(&coarse, 0.1, &seed, &config);
        let r2 = upscale_with_fbm(&coarse, 0.1, &seed, &config);

        assert_eq!(
            r1.heightmap.data, r2.heightmap.data,
            "Same seed should produce identical output"
        );
    }

    #[test]
    fn upscale_adds_detail_on_slopes() {
        let n = 32;
        let mut coarse = GridF32::new(n, n, 0.0);
        for j in 0..n {
            for i in 0..n {
                coarse.set(i, j, i as f32 / n as f32 * 0.6 + 0.2);
            }
        }

        let seed = WorldSeed::new(42);
        let config = FbmUpscaleConfig { target_size: 256, ..Default::default() };

        let result = upscale_with_fbm(&coarse, 0.1, &seed, &config);
        let mean = result.heightmap.mean();

        let variance = result
            .heightmap
            .data
            .iter()
            .map(|&v| {
                let dev = v - mean;
                dev * dev
            })
            .sum::<f32>()
            / result.heightmap.data.len() as f32;

        assert!(variance > 1e-5, "Should have measurable noise variance: {variance}");
    }

    #[test]
    fn upscale_output_in_range() {
        let mut coarse = GridF32::new(16, 16, 0.0);
        for j in 0..16 {
            for i in 0..16 {
                coarse.set(i, j, (i + j) as f32 / 30.0);
            }
        }

        let seed = WorldSeed::new(123);
        let config = FbmUpscaleConfig { target_size: 128, ..Default::default() };

        let result = upscale_with_fbm(&coarse, 0.1, &seed, &config);

        for &v in &result.heightmap.data {
            assert!((0.0..=1.0).contains(&v), "Height out of range: {v}");
        }
    }
}
