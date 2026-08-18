//! Anisotropic FBM upscaling: transforms a coarse isostatic heightmap into
//! detailed terrain by adding fractal noise modulated by local slope.

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use super::noise::SeededNoise;
use crate::erosion::hydraulic::ErosionConfig;
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
    /// **Sub-domain window ORIGIN** in normalized coarse coordinates `[0,1]`
    /// `(x, y)`. The upscale samples only the window `[origin, origin+size]²`
    /// of the coarse field instead of the whole (periodic) torus, so a fraction
    /// of the domain is rendered at the full `target_size` (→ finer km/cell).
    /// Default `[0.0, 0.0]`. The coarse grid is sampled periodically, so a window
    /// whose `origin+size` exceeds 1 simply wraps.
    ///
    /// `skip_serializing_if` on the full-domain default keeps the eroded cache
    /// key byte-identical for non-windowed configs (existing `.raw` stay valid).
    #[serde(default, skip_serializing_if = "is_full_domain_origin")]
    pub sample_origin: [f64; 2],
    /// **Sub-domain window SIZE** as a fraction of the coarse field `(0,1]`.
    /// `1.0` (default) = the whole torus → byte-identical to pre-window. `< 1`
    /// zooms into a window rendered at `target_size`.
    #[serde(default = "default_sample_size", skip_serializing_if = "is_full_sample_size")]
    pub sample_size: f64,
    /// Domain warp strength as fraction of base noise frequency. Default: 0.4.
    /// 0.0 = no warping, 0.5 = moderate, 1.0 = heavy distortion.
    pub domain_warp_strength: f64,
    /// Frequency of the warp noise relative to base frequency. Default: 0.5.
    pub domain_warp_frequency: f64,
    /// Number of FBM octaves for the warp noise itself. Default: 3.
    pub domain_warp_octaves: usize,
    /// **Coastline warp** (Issue #151 follow-up). Displacement applied
    /// to the COARSE-ALTITUDE sampling coordinates `(sx, sy)` BEFORE
    /// bilinear interpolation, in units of COARSE pixels. This makes
    /// the `altitude = sea_level` contour (the coastline) meander
    /// instead of following the blocky interpolated coarse polygon —
    /// the domain warp above only warps the NOISE, never the coarse
    /// sampling, so it cannot move the coast (cause pinned: blocky
    /// coast = STEP 1, the interpolated contour). Amplitude ~0.5–1.0
    /// coarse cells breaks the 1-cell stairstep without making the
    /// coast chaotic. Default: 0.0 (OFF → byte-identical to pre-#151;
    /// the existing v2 pipeline is unaffected).
    pub coast_warp_strength: f64,
    /// Frequency of the coastline-warp noise, in cycles per COARSE
    /// pixel. Low = coherent meander over several coarse cells.
    /// Default: 0.5.
    pub coast_warp_frequency: f64,
    /// **Coastal amplitude taper** (Issue #151 follow-up). FBM amplitude
    /// is multiplied by `smoothstep(|height − sea_level|, 0,
    /// coastal_amplitude_band)`, so the noise tapers to ~0 AT the
    /// sea-level contour and rises to full beyond `band` (altitude
    /// units). Without it, `±amplitude·noise` flips cells land↔ocean at
    /// the coast → a feathered/combed coastline, worst on flat
    /// near-sea-level regions (the interior of strong relief is
    /// unaffected — so a global amplitude cut is the wrong fix, it would
    /// flatten mountains). This is a LOCAL fix: kills coastal feathering,
    /// keeps inland mountain detail. The coast warp supplies the macro
    /// coast irregularity; this removes the micro-feathering. Default:
    /// 0.0 (OFF → byte-identical to pre-#151; v2 unaffected).
    pub coastal_amplitude_band: f64,
    /// **HD hydraulic erosion** (#155 méso). When `Some`, [`upscale_from_c1`]
    /// (tectonics_c1::production_upscale) applies droplet hydraulic erosion
    /// to the HD heightmap AFTER the FBM step — the dendritic dissection that
    /// turns the C1 macro ridge into credible eroded mountains (the méso bar
    /// for Living Landz; validated by the 2048² état-des-lieux). Default
    /// `None` → no erosion (byte-identical; smoke/regression/probe tests and
    /// v2's `upscale_with_fbm` path are unaffected — only `upscale_from_c1`
    /// reads this field). The canonical C1 HD product config that turns it
    /// ON is [`FbmUpscaleConfig::c1_hd_production`].
    #[serde(default)]
    pub erosion: Option<ErosionConfig>,
    /// **Submarine bathymetry re-map** (#submarine). When `Some`,
    /// [`upscale_from_c1`](crate::tectonics_c1::production_upscale::upscale_from_c1)
    /// re-maps the ocean floor toward the plateau→slope→abyss envelope AFTER the
    /// FBM + erosion (see [`crate::terrain::bathymetry`]) — gives the missing
    /// continental shelf + abyssal plain the Stein-Stein slab lacked. Touches
    /// ONLY sub-sea cells (coastline / land invariant). Default `None` →
    /// byte-identical (only `upscale_from_c1` reads it; v2's `upscale_with_fbm`
    /// path is unaffected). The canonical C1 HD product config that turns it ON is
    /// [`FbmUpscaleConfig::c1_hd_production`].
    #[serde(default)]
    pub bathymetry: Option<crate::terrain::bathymetry::BathymetryProfile>,
}

/// serde default for [`FbmUpscaleConfig::sample_size`] (a missing field must
/// deserialize to the full domain `1.0`, not `f64`'s `0.0`).
fn default_sample_size() -> f64 {
    1.0
}

/// A full-domain origin serializes to nothing (keeps the pre-window cache key).
fn is_full_domain_origin(o: &[f64; 2]) -> bool {
    o[0] == 0.0 && o[1] == 0.0
}

/// A full-domain size (`1.0`) serializes to nothing (keeps the pre-window key).
#[allow(clippy::trivially_copy_pass_by_ref)]
fn is_full_sample_size(s: &f64) -> bool {
    *s == 1.0
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
            sample_origin: [0.0, 0.0],
            sample_size: 1.0,
            domain_warp_strength: 0.0,
            domain_warp_frequency: 0.5,
            domain_warp_octaves: 3,
            coast_warp_strength: 0.0,
            coast_warp_frequency: 0.5,
            coastal_amplitude_band: 0.0,
            erosion: None,
            bathymetry: None,
        }
    }
}

impl FbmUpscaleConfig {
    /// Canonical C1 HD product config (#155) — the #151 coastline params
    /// + HD hydraulic erosion ON (the validated 2048² état-des-lieux
    /// product). `target_size` parameterised; `num_droplets` is scaled at
    /// the judged density (4M at 2048²) ∝ `target_size²` so the dissection
    /// is resolution-coherent (#151 spirit), 2048² == 4M exactly (the judged
    /// point). The path that DELIVERS the HD product must use THIS (else the
    /// erosion-off product — judged non-deliverable — ships).
    ///
    /// COST (product characteristic): the 4M-droplet erosion at 2048² is
    /// ~2 min/seed (measured 105–125 s). This is an OFFLINE / background
    /// export path, NOT interactive — the Living Landz pipeline must generate
    /// the HD product in the background, not on-the-fly.
    ///
    /// NOTE (#155 follow-up): `ErosionConfig::sea_level` defaults to 0.1 and
    /// is PRESERVED here as-judged — the état-des-lieux rendered with that
    /// (mismatched) value. The normalised sea level is actually 0.5; setting
    /// `sea_level = 0.5` is the CORRECT value but CHANGES coastal deposition
    /// → a separate follow-up maillon to re-judge the coastal render, NOT a
    /// silent fix here. Do not "correct" it in passing.
    #[must_use]
    pub fn c1_hd_production(target_size: usize) -> Self {
        let num_droplets = (4_000_000u64 * (target_size as u64).pow(2) / (2048u64).pow(2)) as usize;
        Self {
            target_size,
            coast_warp_strength: 1.5,
            coast_warp_frequency: 0.5,
            coastal_amplitude_band: 0.30,
            amplitude_base: 0.16,
            submarine_damping: 0.0,
            erosion: Some(ErosionConfig {
                num_droplets,
                batch_size: 100_000,
                // sea_level 0.1 preserved as-judged (see note above).
                ..Default::default()
            }),
            // #submarine — re-map the ocean floor to the plateau→slope→abyss
            // envelope (the diagnostic found a uniform ~−2600 m slab with almost
            // no shelf and no abyss). Anchored default profile; touches only
            // sub-sea cells. None on the plain Default keeps v2/byte-identical.
            bathymetry: Some(crate::terrain::bathymetry::BathymetryProfile::default()),
            ..Default::default()
        }
    }
}

/// Result of FBM upscaling.
pub struct UpscaleResult {
    /// The upscaled heightmap at target resolution.
    pub heightmap: GridF32,
    /// The slope magnitude field (useful for erosion and viz).
    pub slope: GridF32,
    /// Sediment / water-passage map from HD erosion (#155) — `Some` only
    /// when `cfg.erosion` ran. Forwarded from `ErosionResult.sediment`; the
    /// hook for a future rivers/lakes chantier (NOT consumed yet).
    pub sediment: Option<GridF32>,
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
    let target = config.target_size;

    // Output dimensions: `target_size` is applied to the longer source
    // axis; the other axis is derived to preserve the source aspect ratio.
    // Must stay in sync with `export::upscale_dims`, which documents the
    // same contract for the exported artefacts.
    let (dst_w, dst_h) = if src_w == src_h {
        (target, target)
    } else if src_w >= src_h {
        let h = ((target as f64) * src_h as f64 / src_w as f64).round() as usize;
        (target, h.max(1))
    } else {
        let w = ((target as f64) * src_w as f64 / src_h as f64).round() as usize;
        (w.max(1), target)
    };

    // Sub-domain window: sample only `[origin, origin+size]` of the coarse
    // field. `size == 1` and `origin == 0` (defaults) → the full domain, and
    // `scale`/`origin` collapse to the pre-window `src/dst` mapping (byte-
    // identical). A smaller `size` maps `dst` pixels onto fewer coarse cells →
    // finer effective km/cell. Coarse sampling is periodic (`_periodic`), so a
    // window crossing the torus edge wraps.
    let scale_x = config.sample_size * src_w as f64 / dst_w as f64;
    let scale_y = config.sample_size * src_h as f64 / dst_h as f64;
    let origin_x = config.sample_origin[0] * src_w as f64;
    let origin_y = config.sample_origin[1] * src_h as f64;

    // Create noise generators
    let noise_seed = seed.derive_seed("fbm_upscale") as u32;
    let noise = SeededNoise::new(noise_seed, config.octaves);
    // Separate single-octave source for angle perturbation (different seed)
    let angle_noise = SeededNoise::new(noise_seed.wrapping_add(99991), 1);
    // Domain warp: two independent FBM fields for X and Y displacement
    let warp_noise_x = SeededNoise::new(noise_seed.wrapping_add(55555), config.domain_warp_octaves);
    let warp_noise_y = SeededNoise::new(noise_seed.wrapping_add(77777), config.domain_warp_octaves);
    // Dedicated coastline-warp noise (Issue #151 follow-up) — distinct
    // seeds from the noise/domain-warp so the coast displacement is
    // independent of the FBM detail.
    let coast_octaves = config.domain_warp_octaves.max(1);
    let coast_warp_noise_x = SeededNoise::new(noise_seed.wrapping_add(13331), coast_octaves);
    let coast_warp_noise_y = SeededNoise::new(noise_seed.wrapping_add(24421), coast_octaves);

    // Precompute slope and direction on the coarse grid
    let (slope_map, direction_map) = compute_terrain_analysis(coarse);

    // #151: RESOLUTION-INDEPENDENT noise frequency. The FBM is sampled in
    // COARSE-CELL coordinate space (`sx, sy`), NOT target-pixel space, so
    // the noise feature size is fixed relative to the terrain regardless of
    // `target_size`. (Previously `freq = base/src_w` applied to the TARGET
    // pixel index → the feature count scaled with the upscale ratio, e.g.
    // 4× finer at 4096² than 1024² for the same config — see #151.)
    //
    // The coefficients reference the prior 1024² calibration
    // (`NOISE_REF_TARGET`) so that a 1024² render is BYTE-IDENTICAL to
    // pre-#151 (old `i·base/src` ≡ `sx·nscale` when `target == 1024`); all
    // other target sizes now match the 1024² feature size instead of
    // diverging.
    const NOISE_REF_TARGET: f64 = 1024.0;
    let src_max = src_w.max(src_h) as f64;
    let nscale = config.base_frequency * NOISE_REF_TARGET / (src_max * src_max);
    let ascale = ANGLE_PERTURBATION_FREQ * NOISE_REF_TARGET / src_max;

    // Process each output row in parallel
    let row_data: Vec<(Vec<f32>, Vec<f32>)> = (0..dst_h)
        .into_par_iter()
        .map(|j| {
            let mut h_row = vec![0.0f32; dst_w];
            let mut s_row = vec![0.0f32; dst_w];

            for i in 0..dst_w {
                // Source coordinates in coarse pixel space (offset into the
                // sub-domain window; see `origin_x`/`scale_x` above).
                let sx = origin_x + i as f64 * scale_x;
                let sy = origin_y + j as f64 * scale_y;

                // 0. Coastline warp (Issue #151 follow-up): displace the
                // COARSE-ALTITUDE sampling position so the sea-level
                // contour meanders instead of following the blocky 64²
                // polygon. In coarse-pixel units; the domain warp (step
                // 7) never touches this. OFF when strength == 0 →
                // (csx, csy) == (sx, sy), byte-identical to pre-#151.
                let (csx, csy) = if config.coast_warp_strength > 0.0 {
                    let cwf = config.coast_warp_frequency;
                    let (wcx, wcy) = (sx * cwf, sy * cwf);
                    let cdx = coast_warp_noise_x.fbm(
                        wcx,
                        wcy,
                        coast_octaves,
                        config.lacunarity,
                        config.persistence,
                    ) * config.coast_warp_strength;
                    let cdy = coast_warp_noise_y.fbm(
                        wcx,
                        wcy,
                        coast_octaves,
                        config.lacunarity,
                        config.persistence,
                    ) * config.coast_warp_strength;
                    (sx + cdx, sy + cdy)
                } else {
                    (sx, sy)
                };

                // 1. Bilinear interpolation of the coarse heightmap (at the
                // coast-warped position — this is what bends the coastline).
                let base_height = coarse.sample_bilinear_periodic(csx as f32, csy as f32);

                // 2. Sample terrain properties from coarse analysis (same
                // warped position so amplitude/orientation follow the
                // displaced coarse terrain coherently).
                let slope_mag = slope_map.sample_bilinear_periodic(csx as f32, csy as f32);
                let slope_dir = direction_map.sample_bilinear_periodic(csx as f32, csy as f32);

                // 3. Compute amplitude modulation
                let altitude_factor =
                    if base_height > sea_level { 1.0 } else { config.submarine_damping };

                let mut amplitude = config.amplitude_base
                    * (1.0 + slope_mag as f64 * config.amplitude_slope_factor)
                    * altitude_factor;

                // #151: coastal amplitude taper — damp the FBM to ~0 AT the
                // sea-level contour so it doesn't feather the coastline by
                // flipping near-sea cells land↔ocean; full amplitude inland
                // (mountains preserved). OFF (band 0) → no-op, byte-identical.
                if config.coastal_amplitude_band > 0.0 {
                    let dist_from_sea = (base_height as f64 - sea_level as f64).abs();
                    amplitude *= smoothstep(dist_from_sea, 0.0, config.coastal_amplitude_band);
                }

                // 4. Compute anisotropy ratio with sigmoid rolloff
                let slope_f64 = slope_mag as f64;
                let aniso_t = smoothstep(slope_f64, 0.0, 1.0);
                let anisotropy = 1.0 + (config.max_anisotropy - 1.0) * aniso_t;

                // 5. Blend factor: isotropic at low slopes, anisotropic at high slopes
                let aniso_blend = smoothstep(slope_f64, ISOTROPY_LOW, ISOTROPY_HIGH);

                // 6. Angular perturbation to break long-range parallelism.
                // #151: sampled in coarse-cell space (sx·ascale) →
                // resolution-independent, byte-identical at 1024².
                let angle_offset =
                    angle_noise.sample(0, sx * ascale, sy * ascale) * ANGLE_PERTURBATION_MAX;
                let perturbed_dir = slope_dir as f64 + angle_offset;

                // 7. Domain warping: distort noise coordinates to break regular
                // patterns. #151: all in coarse-cell space (sx·nscale).
                let raw_nx = sx * nscale;
                let raw_ny = sy * nscale;

                let (nx, ny) = if config.domain_warp_strength > 0.0 {
                    let warp_freq = nscale * config.domain_warp_frequency;
                    let wx = sx * warp_freq;
                    let wy = sy * warp_freq;
                    let inv_freq = config.domain_warp_strength / nscale;
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
    let mut heightmap = GridF32::new(dst_w, dst_h, 0.0);
    let mut slope_out = GridF32::new(dst_w, dst_h, 0.0);

    for (j, (h_row, s_row)) in row_data.into_iter().enumerate() {
        for (i, (h, s)) in h_row.into_iter().zip(s_row).enumerate() {
            heightmap.set(i, j, h);
            slope_out.set(i, j, s);
        }
    }

    UpscaleResult { heightmap, slope: slope_out, sediment: None }
}

/// Compute slope magnitude and direction on the coarse grid.
fn compute_terrain_analysis(heightmap: &GridF32) -> (GridF32, GridF32) {
    let w = heightmap.width;
    let h = heightmap.height;
    let mut slope = GridF32::new(w, h, 0.0);
    let mut direction = GridF32::new(w, h, 0.0);

    for j in 0..h {
        for i in 0..w {
            let (gx, gy) = heightmap.gradient_at_periodic(i, j);
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
    fn full_domain_window_is_byte_identical_to_default() {
        // An explicit full-domain window ([0,0], 1.0) must reproduce the default
        // (no-window) output bit-for-bit — the byte-identical guard.
        let n = 24;
        let mut coarse = GridF32::new(n, n, 0.0);
        for j in 0..n {
            for i in 0..n {
                coarse.set(i, j, i as f32 / n as f32 * 0.6 + 0.2);
            }
        }
        let seed = WorldSeed::new(7);
        let base = FbmUpscaleConfig { target_size: 96, ..Default::default() };
        let explicit =
            FbmUpscaleConfig { sample_origin: [0.0, 0.0], sample_size: 1.0, ..base.clone() };
        let r_def = upscale_with_fbm(&coarse, 0.1, &seed, &base);
        let r_exp = upscale_with_fbm(&coarse, 0.1, &seed, &explicit);
        assert_eq!(r_def.heightmap.data, r_exp.heightmap.data);
    }

    #[test]
    fn sub_domain_window_zooms_into_region() {
        // Noise off (amplitude 0) → heightmap == bilinear coarse samples, so the
        // window maps deterministically onto its coarse sub-region. A 0.5-wide
        // window at origin 0.25 renders coarse x∈[0.25,0.75] across the target.
        let n = 32;
        let mut coarse = GridF32::new(n, n, 0.0);
        for j in 0..n {
            for i in 0..n {
                coarse.set(i, j, i as f32 / n as f32); // 0..1 gradient in x
            }
        }
        let seed = WorldSeed::new(1);
        let cfg = FbmUpscaleConfig {
            target_size: n,
            amplitude_base: 0.0, // kill the FBM → pure coarse resampling
            sample_origin: [0.25, 0.0],
            sample_size: 0.5,
            ..Default::default()
        };
        let r = upscale_with_fbm(&coarse, 0.1, &seed, &cfg);
        let left = r.heightmap.get(0, 0); // sx = 0.25·n = 8 → 0.25
        let right = r.heightmap.get((n - 1) as i32, 0); // sx = 8 + 0.5·31 = 23.5 → ~0.734
        assert!((left - 0.25).abs() < 0.02, "window left edge ~0.25, got {left}");
        assert!((right - 0.734).abs() < 0.03, "window right edge ~0.734, got {right}");
        // Half-domain window → about half the x-span a full render would show.
        assert!((0.35..0.6).contains(&(right - left)), "half-domain span, got {}", right - left);
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

    /// The upscale must preserve the source aspect ratio: target_size
    /// becomes the longer output axis, the shorter axis is derived.
    /// Source 16×10 (landscape, 8:5) → 128×80 at target_size=128.
    #[test]
    fn upscale_preserves_aspect_ratio_landscape() {
        let mut coarse = GridF32::new(16, 10, 0.5);
        for j in 0..10 {
            for i in 0..16 {
                coarse.set(i, j, (i + j) as f32 / 26.0);
            }
        }
        let seed = WorldSeed::new(42);
        let config = FbmUpscaleConfig { target_size: 128, ..Default::default() };
        let result = upscale_with_fbm(&coarse, 0.1, &seed, &config);

        assert_eq!(result.heightmap.width, 128);
        assert_eq!(result.heightmap.height, 80);
        assert_eq!(result.slope.width, 128);
        assert_eq!(result.slope.height, 80);
        assert_eq!(result.heightmap.data.len(), 128 * 80);
    }

    /// Same as above with a portrait source: 10×16 → 80×128.
    #[test]
    fn upscale_preserves_aspect_ratio_portrait() {
        let coarse = GridF32::new(10, 16, 0.5);
        let seed = WorldSeed::new(42);
        let config = FbmUpscaleConfig { target_size: 128, ..Default::default() };
        let result = upscale_with_fbm(&coarse, 0.1, &seed, &config);

        assert_eq!(result.heightmap.width, 80);
        assert_eq!(result.heightmap.height, 128);
        assert_eq!(result.slope.width, 80);
        assert_eq!(result.slope.height, 128);
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
