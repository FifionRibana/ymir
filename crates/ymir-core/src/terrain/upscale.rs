//! Anisotropic FBM upscaling: transforms a coarse isostatic heightmap into
//! detailed terrain by adding fractal noise modulated by local slope.

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use super::noise::SeededNoise;
use crate::grid::GridF32;
use crate::seed::WorldSeed;

/// Configuration for anisotropic FBM upscaling.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
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
    /// Domain warp strength as fraction of base noise frequency. Default: 0.4.
    /// 0.0 = no warping, 0.5 = moderate, 1.0 = heavy distortion.
    pub domain_warp_strength: f64,
    /// Frequency of the warp noise relative to base frequency. Default: 0.5.
    pub domain_warp_frequency: f64,
    /// Number of FBM octaves for the warp noise itself. Default: 3.
    pub domain_warp_octaves: usize,
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
            domain_warp_strength: 0.0,
            domain_warp_frequency: 0.5,
            domain_warp_octaves: 3,
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

/// Hermite smoothstep: 0 at lo, 1 at hi, smooth transition.
fn smoothstep(x: f64, lo: f64, hi: f64) -> f64 {
    let t = ((x - lo) / (hi - lo)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

// Thresholds for isotropic/anisotropic blending (slope magnitude).
const ISOTROPY_LOW: f64 = 0.01;
const ISOTROPY_HIGH: f64 = 0.05;

// Maximum angular perturbation (radians, ~17 degrees).
const ANGLE_PERTURBATION_MAX: f64 = 0.3;
// Frequency of the angle perturbation noise (low = spatially coherent).
const ANGLE_PERTURBATION_FREQ: f64 = 0.02;

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

    // Create noise generators
    let noise_seed = seed.derive_seed("fbm_upscale") as u32;
    let noise = SeededNoise::new(noise_seed, config.octaves);
    // Separate single-octave source for angle perturbation (different seed)
    let angle_noise = SeededNoise::new(noise_seed.wrapping_add(99991), 1);
    // Domain warp: two independent FBM fields for X and Y displacement
    let warp_noise_x = SeededNoise::new(noise_seed.wrapping_add(55555), config.domain_warp_octaves);
    let warp_noise_y = SeededNoise::new(noise_seed.wrapping_add(77777), config.domain_warp_octaves);

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

                // 4. Compute anisotropy ratio with sigmoid rolloff
                let slope_f64 = slope_mag as f64;
                let aniso_t = smoothstep(slope_f64, 0.0, 1.0);
                let anisotropy = 1.0 + (config.max_anisotropy - 1.0) * aniso_t;

                // 5. Blend factor: isotropic at low slopes, anisotropic at high slopes
                let aniso_blend = smoothstep(slope_f64, ISOTROPY_LOW, ISOTROPY_HIGH);

                // 6. Angular perturbation to break long-range parallelism
                let angle_offset = angle_noise.sample(
                    0,
                    i as f64 * ANGLE_PERTURBATION_FREQ,
                    j as f64 * ANGLE_PERTURBATION_FREQ,
                ) * ANGLE_PERTURBATION_MAX;
                let perturbed_dir = slope_dir as f64 + angle_offset;

                // 7. Domain warping: distort noise coordinates to break regular patterns
                let raw_nx = i as f64 * freq;
                let raw_ny = j as f64 * freq;

                let (nx, ny) = if config.domain_warp_strength > 0.0 {
                    let warp_freq = freq * config.domain_warp_frequency;
                    let wx = i as f64 * warp_freq;
                    let wy = j as f64 * warp_freq;
                    let inv_freq = config.domain_warp_strength / freq;
                    let warp_dx = warp_noise_x.fbm(
                        wx,
                        wy,
                        config.domain_warp_octaves,
                        config.lacunarity,
                        config.persistence,
                    ) * inv_freq;
                    let warp_dy = warp_noise_y.fbm(
                        wx,
                        wy,
                        config.domain_warp_octaves,
                        config.lacunarity,
                        config.persistence,
                    ) * inv_freq;
                    (raw_nx + warp_dx, raw_ny + warp_dy)
                } else {
                    (raw_nx, raw_ny)
                };

                // 8. Sample FBM
                let noise_val = if aniso_blend < 0.001 {
                    // Pure isotropic — skip aniso sample
                    noise.fbm(nx, ny, config.octaves, config.lacunarity, config.persistence)
                } else if aniso_blend > 0.999 {
                    // Pure anisotropic — skip iso sample
                    noise.fbm_anisotropic(
                        nx,
                        ny,
                        perturbed_dir,
                        anisotropy,
                        config.octaves,
                        config.lacunarity,
                        config.persistence,
                    )
                } else {
                    // Blend both
                    let fbm_iso =
                        noise.fbm(nx, ny, config.octaves, config.lacunarity, config.persistence);
                    let fbm_aniso = noise.fbm_anisotropic(
                        nx,
                        ny,
                        perturbed_dir,
                        anisotropy,
                        config.octaves,
                        config.lacunarity,
                        config.persistence,
                    );
                    fbm_iso + (fbm_aniso - fbm_iso) * aniso_blend
                };

                // 9. Combine: base + noise
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

    #[test]
    fn flat_interior_uses_isotropic_noise() {
        // On a flat heightmap, the interior pixels (away from edges) should produce
        // identical noise regardless of max_anisotropy, because slope_mag ≈ 0
        // and the blend factor forces pure isotropic sampling.
        // Edge pixels may differ because GridF32::gradient_at returns non-zero
        // gradients at boundaries (clamp behavior).
        let coarse = GridF32::new(16, 16, 0.5);
        let seed = WorldSeed::new(42);

        let config_iso =
            FbmUpscaleConfig { target_size: 64, max_anisotropy: 1.0, ..Default::default() };
        let config_aniso =
            FbmUpscaleConfig { target_size: 64, max_anisotropy: 3.0, ..Default::default() };

        let r_iso = upscale_with_fbm(&coarse, 0.1, &seed, &config_iso);
        let r_aniso = upscale_with_fbm(&coarse, 0.1, &seed, &config_aniso);

        // Check interior pixels only (skip outer 25% to avoid edge effects)
        let margin = 16; // 25% of 64
        let mut max_diff = 0.0f32;
        for j in margin..(64 - margin) {
            for i in margin..(64 - margin) {
                let idx = j * 64 + i;
                let diff = (r_iso.heightmap.data[idx] - r_aniso.heightmap.data[idx]).abs();
                max_diff = max_diff.max(diff);
            }
        }

        assert!(
            max_diff < 1e-6,
            "Flat interior should be isotropic regardless of max_anisotropy, max_diff={max_diff}"
        );
    }
}
