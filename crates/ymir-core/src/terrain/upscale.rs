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
    /// **Flow-conditioning relief budget** β (C-1, closures roadmap §1). The
    /// additive isotropic FBM is the sole creator of closed depressions (16
    /// tectonic → 90 682 after the FBM): a pit forms wherever the noise out-slopes
    /// the bed and reverses the descent. When `> 0`, two coupled mechanisms
    /// suppress that:
    ///
    /// 1. **Amplitude cap** (this knob). The noise amplitude is bounded by a relief
    ///    budget `amplitude ≤ β · slope_mag / (nscale·S)`, where the divisor
    ///    `nscale·S` (S = Σ(persistence·lacunarity)ᵒ) converts a bed slope into the
    ///    FBM amplitude whose summed per-octave downslope rise equals it. Smaller β
    ///    = stricter (fewer pits, less detail); the limit β → 0 recovers the smooth
    ///    coarse bed (its 16 depressions). On a flat (slope → 0) the cap → 0 and
    ///    removes the fabricated pit; on a steep flank it is generous (texture
    ///    kept).
    /// 2. **Downslope stretch** (fixed [`FLOW_STRETCH`], not tuned here). The noise
    ///    is elongated along the bed gradient so its along-flow derivative shrinks —
    ///    a bed displaced laterally without a counter-slope, while the contour axis
    ///    still carries relief.
    ///
    /// Both depend only on the coarse slope field and config, so they are identical
    /// at every `target_size` (band policy unchanged; low bands still bit-identical
    /// across resolutions). Default `0.0` (OFF → byte-identical to the pre-C-1
    /// additive noise; all determinism/byte guards stay green). The canonical HD
    /// product ([`FbmUpscaleConfig::c1_hd_production`]) turns it ON.
    #[serde(default, skip_serializing_if = "is_zero_f64")]
    pub flow_conditioning: f64,
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
    /// **Routed stream-power incision** (ADR 0001, prototype). When `Some`,
    /// [`upscale_from_c1`](crate::tectonics_c1::production_upscale::upscale_from_c1)
    /// carves valleys along the drainage network (Braun & Willett) AFTER the FBM and
    /// BEFORE droplet erosion — deterministic, hierarchy by construction, ~13× faster
    /// than the droplet pass, and (unlike droplets) it RAISES drainage relief instead
    /// of collapsing it. Default `None` → skipped, byte-identical, OFF in production
    /// until confirmed at 8192². Pair it with a WEAK droplet pass (reduced
    /// `ErosionConfig.num_droplets`) for hillslope texture — the full droplet pass
    /// erases the carved valleys (measured relief 323 → 24 m).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stream_power: Option<crate::erosion::stream_power::StreamPowerConfig>,
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
    /// **Sea-level calibration on a target LAND-AREA fraction** (M1). When
    /// `Some(f)`, [`upscale_from_c1`](crate::tectonics_c1::production_upscale::upscale_from_c1)
    /// shifts the COARSE altitude so exactly `f` of cells stay above 0 m (the
    /// `1−f` quantile becomes 0 m → "0 m = coastline" by construction), instead
    /// of leaving sea at the isostatic level (which emerges ~55–60 % land). Read
    /// ONLY by `upscale_from_c1`; `None` (default) → no shift, byte-identical.
    /// The canonical HD config sets `Some(0.29)` (Earth-like ocean fraction).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_land_fraction: Option<f32>,
    /// **C-3 lithological heterogeneity** (closures roadmap §3). When
    /// `lithology.enabled`, `upscale_from_c1` builds a per-cell erodibility
    /// MULTIPLIER (hard basement 1.0, softer on rift + volcaniclastic footprints,
    /// derived causally from the tectonic state — never noise/geometry) and threads
    /// it into the stream-power incision. Default (`enabled = false`) → uniform K,
    /// byte-identical. Only read when `stream_power` is also `Some`. See
    /// [`crate::tectonics_c1::closures::lithology`]. Skipped from serialization when
    /// disabled, so a lithology-off eroded cache key is byte-identical to pre-C-3.
    #[serde(
        default,
        skip_serializing_if = "crate::tectonics_c1::closures::lithology::LithologyConfig::is_disabled"
    )]
    pub lithology: crate::tectonics_c1::closures::lithology::LithologyConfig,
    /// **C-3b inherited structure** (closures roadmap §3b). When `fracture.enabled`,
    /// `upscale_from_c1` builds a DIRECTIONAL fracture field (density × strike, causal
    /// from plate kinematics) and threads it into the stream-power incision so valleys
    /// align on the tectonic fabric (anisotropic K, no relief added → C-1 survives).
    /// Default (`enabled = false`) → byte-identical. Only read when `stream_power` is
    /// `Some`. See [`crate::tectonics_c1::closures::fracture`].
    #[serde(
        default,
        skip_serializing_if = "crate::tectonics_c1::closures::fracture::FractureConfig::is_disabled"
    )]
    pub fracture: crate::tectonics_c1::closures::fracture::FractureConfig,
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

/// A zero (OFF) flow-conditioning strength serializes to nothing (keeps the
/// pre-C-1 eroded cache key byte-identical for un-conditioned configs).
#[allow(clippy::trivially_copy_pass_by_ref)]
fn is_zero_f64(v: &f64) -> bool {
    *v == 0.0
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
            flow_conditioning: 0.0,
            sample_origin: [0.0, 0.0],
            sample_size: 1.0,
            domain_warp_strength: 0.0,
            domain_warp_frequency: 0.5,
            domain_warp_octaves: 3,
            coast_warp_strength: 0.0,
            coast_warp_frequency: 0.5,
            coastal_amplitude_band: 0.0,
            erosion: None,
            stream_power: None,
            bathymetry: None,
            target_land_fraction: None,
            lithology: crate::tectonics_c1::closures::lithology::LithologyConfig::default(),
            fracture: crate::tectonics_c1::closures::fracture::FractureConfig::default(),
        }
    }
}

/// The two competing terms of the per-cell FBM amplitude. The effective amplitude is the
/// SMALLER: a configured LEVEL versus the C-1 relief BUDGET. Named because which one binds
/// is a load-bearing fact — in production the cap binds at EVERY cell, which makes
/// `amplitude_base` inert (ADR "The DEAD KNOB"); a silent `min` hid that for three closures.
#[derive(Debug, Clone, Copy)]
pub struct AmplitudeTerms {
    /// The configured level: `amplitude_base · (1 + slope·factor) · submarine_damping`.
    pub base: f64,
    /// The C-1 relief budget `β · slope / divisor`. `None` when `flow_conditioning == 0`
    /// (the pre-C-1 regime, where `base` genuinely drives the terrain).
    pub cap: Option<f64>,
}

impl AmplitudeTerms {
    /// The amplitude actually applied — the smaller of the two terms.
    #[inline]
    #[must_use]
    pub fn effective(&self) -> f64 {
        match self.cap {
            Some(c) => self.base.min(c),
            None => self.base,
        }
    }
    /// `true` when the relief budget is what limits this cell — i.e. `amplitude_base` has
    /// NO effect here. In production this is true everywhere.
    #[inline]
    #[must_use]
    pub fn cap_binds(&self) -> bool {
        matches!(self.cap, Some(c) if c < self.base)
    }
}

/// Inputs to [`production_hd_config`] — everything the SHIPPED HD path varies per run.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductionHdOpts {
    pub target_size: usize,
    /// Geometric domain span (km) — the domain IS the map. NEVER `geo_scale_ratio`.
    pub domain_km: f32,
    /// Vertical scale for the stream-power incision (`SteinSteinParams::depth_scale_m`).
    pub depth_scale_m: f32,
    /// Framing roll of the coarse sampling origin, and the sampled fraction (always 1.0 —
    /// the whole torus renders).
    pub sample_origin: [f64; 2],
    pub sample_size: f64,
    /// FBM base amplitude. ⚠️ INERT on this path — the C-1 relief-budget cap binds at every
    /// cell (proven byte-identical at 4×; see ADR "The DEAD KNOB"). Kept so the config
    /// states what it *would* be, and so the guard test can assert the cap still binds.
    pub amplitude_base: f64,
    /// MFD partition exponent for relief-v3.
    pub mfd_p: f32,
    pub lithology: crate::tectonics_c1::closures::lithology::LithologyConfig,
    pub fracture: crate::tectonics_c1::closures::fracture::FractureConfig,
}

/// **THE single source of truth for the SHIPPED HD config.** Returns exactly what
/// production runs — relief-v3 stream-power incision, droplets off, the C-1 conditioning,
/// the closures' configs — with NO further mutation expected from the caller.
///
/// This exists because `c1_hd_production` is NOT production: the viz used to build it and
/// then mutate amplitude, sampling, erosion, stream-power, lithology and fracture, so any
/// bench calling `c1_hd_production` got something else than what ships. That divergence
/// occurred SEVEN times, most recently hiding an inert `amplitude_base` for three closures.
/// A rule that must be remembered seven times is a design flaw: the viz and the benches now
/// call THIS, and `production_upscale_config_is_the_shipped_one` guards it.
#[must_use]
pub fn production_hd_config(o: &ProductionHdOpts) -> FbmUpscaleConfig {
    let mut cfg = FbmUpscaleConfig::c1_hd_production(o.target_size);
    cfg.amplitude_base = o.amplitude_base;
    cfg.sample_origin = o.sample_origin;
    cfg.sample_size = o.sample_size;
    // relief-v3 replaces the droplet pass: droplets collapse the stream-power valleys.
    cfg.erosion = None;
    let km_per_cell = o.domain_km / o.target_size as f32;
    let cell_km2 = km_per_cell * km_per_cell;
    let mut sp =
        crate::erosion::stream_power::StreamPowerConfig::relief_v3(cell_km2, o.depth_scale_m);
    sp.mfd_exponent = Some(o.mfd_p);
    sp.iterations = 2;
    cfg.stream_power = Some(sp);
    cfg.lithology = o.lithology.clone();
    cfg.fracture = o.fracture.clone();
    cfg
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
    /// M1: `ErosionConfig::sea_level` is set to **0.5** — the real normalised sea
    /// level (`C1_SEA_LEVEL_NORM`). It was 0.1 (a mismatch preserved "as-judged"
    /// for the earlier état-des-lieux); M1 IS the judged change. Coastal impact:
    /// hydraulic erosion's coastal deposition triggers at the TRUE waterline now,
    /// so beach/delta deposition lands at the actual coast (0 m) instead of ~0.1
    /// norm (deep ocean) — deposition that previously vanished offshore now
    /// shapes the real shoreline. Combined with `target_land_fraction` (below).
    ///
    /// M1 #190 REGRESSION FIX: `target_land_fraction` defaults to `None` (no
    /// sea-level calibration). Quantile calibration moves sea level onto the FLAT
    /// continental-shelf plateau (measured slope ~2 000 cell/unit at tlf 0.08 vs
    /// ~12 350 at the isostatic level), so the 0-crossing becomes hypersensitive
    /// → speckled coasts + marginal land (the seed-42 regression). The isostatic
    /// sea level sits on the STEEP part of the hypsometric curve → crisp coasts.
    /// Calibration stays available as an OPT-IN knob (set `target_land_fraction`
    /// on the returned config); it is just not defaulted on. Bounding the
    /// landmass without drowning it is a separate lever (continental_fraction).
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
            // C-1 (closures roadmap §1): flow-condition the FBM so it stops being
            // the sole creator of closed depressions. β = 0.1 cuts the post-FBM pit
            // population ~13× at 8192² / ~28× at 2048² while PRESERVING mountain
            // morphology (steep-slope shares and local relief are held or sharpened
            // — the downslope stretch concentrates relief into coherent valleys).
            // See the C-1 section of docs/adr/0001 for the trajectory table.
            flow_conditioning: 0.1,
            // M1 #190 regression fix: NO sea-level calibration by default (it lands
            // on the flat shelf → speckled coasts). Opt-in only. See docstring.
            target_land_fraction: None,
            erosion: Some(ErosionConfig {
                num_droplets,
                batch_size: 100_000,
                // M1: real normalised sea level (was 0.1 — see docstring).
                sea_level: 0.5,
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

// C-1: fixed downslope elongation applied by the flow conditioning. The noise is
// stretched ×`FLOW_STRETCH` along the bed gradient (its along-flow frequency is
// divided by this), so the along-flow derivative shrinks and the noise cannot
// reverse the descent, while the across-flow (contour) axis keeps full-frequency
// relief. This is a structural part of the method, not a tuning knob (stretching
// only lowers a frequency, so it never crosses Nyquist — the band policy holds).
// The tuning knob is the relief budget `flow_conditioning` (β) on the amplitude.
const FLOW_STRETCH: f64 = 8.0;

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

    // C-1 relief budget: the divisor `nscale · S` converts a slope (Δh per coarse
    // cell) into the FBM amplitude whose summed per-octave downslope rise equals
    // that slope. `S = Σ (persistence·lacunarity)ᵒ` is the octave gradient sum
    // (≈ `octaves` at the p·l = 1 default). Precomputed once — config-only, so it
    // is identical at every resolution. See `flow_conditioning` doc.
    let flow_budget_divisor = if config.flow_conditioning > 0.0 {
        let ratio = config.persistence * config.lacunarity;
        let s: f64 = (0..config.octaves).map(|o| ratio.powi(o as i32)).sum();
        nscale * s
    } else {
        0.0 // unused when conditioning is OFF
    };

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

                // AMPLITUDE COMPOSITION, made EXPLICIT (see ADR "The DEAD KNOB"). Two terms
                // compete and the SMALLER wins; the silent `min` that used to live here is
                // exactly what hid, for three closures, that the cap wins EVERYWHERE in
                // production — making `amplitude_base` (and the viz amplitude selector, and
                // the striation ladder) inert. Named terms so the fact is legible and
                // `amplitude_cap_binds_everywhere_in_production` can assert it.
                let terms = AmplitudeTerms {
                    // (1) the configured LEVEL, slope-boosted and submarine-damped.
                    base: config.amplitude_base
                        * (1.0 + slope_mag as f64 * config.amplitude_slope_factor)
                        * altitude_factor,
                    // (2) the C-1 relief BUDGET cap (β = flow_conditioning): on a flat
                    // (slope → 0) no flow direction exists, so any additive bump is a
                    // fabricated pit; the cap → 0 there and removes it. On a slope it scales
                    // with the bed slope. `None` when conditioning is off (pre-C-1 regime,
                    // where the base term genuinely drives the terrain).
                    cap: (config.flow_conditioning > 0.0)
                        .then(|| config.flow_conditioning * slope_mag as f64 / flow_budget_divisor),
                };
                let mut amplitude = terms.effective();

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
                let noise_val = if config.flow_conditioning > 0.0 {
                    // C-1 flow-aligned sample: STRETCH features DOWNSLOPE by
                    // `flow_conditioning = E` — pass the slope direction with
                    // `ratio = 1/E < 1`, which LOWERS the along-slope frequency
                    // (wavelength ×E) while leaving the across-slope (contour)
                    // frequency untouched. The along-flow derivative shrinks by ≈E
                    // so the noise cannot reverse the descent; the contour axis
                    // still carries full-frequency relief (downslope flutes, not
                    // transverse ridges). Stretching only ever LOWERS a frequency,
                    // so it never crosses Nyquist — the band policy holds (compressing
                    // the contour axis instead would alias into salt-and-pepper pits).
                    // Applied at ALL slopes (unlike the legacy `aniso_blend`, which
                    // stays isotropic on gentle slopes — exactly where pits form); on
                    // a true flat the amplitude cap has already zeroed the noise, so
                    // the ill-defined slope_dir there is harmless.
                    noise.fbm_anisotropic(
                        nx,
                        ny,
                        perturbed_dir,
                        1.0 / FLOW_STRETCH,
                        config.octaves,
                        config.lacunarity,
                        config.persistence,
                    )
                } else if aniso_blend < 0.001 {
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

    /// Count interior cells strictly lower than all 8 neighbours (closed pits).
    fn local_minima(g: &GridF32) -> usize {
        let (w, h) = (g.width, g.height);
        let mut n = 0;
        for y in 1..h - 1 {
            for x in 1..w - 1 {
                let v = g.data[y * w + x];
                let mut is_min = true;
                for dy in -1i32..=1 {
                    for dx in -1i32..=1 {
                        if dx == 0 && dy == 0 {
                            continue;
                        }
                        if g.data[(y as i32 + dy) as usize * w + (x as i32 + dx) as usize] <= v {
                            is_min = false;
                        }
                    }
                }
                if is_min {
                    n += 1;
                }
            }
        }
        n
    }

    /// C-1 permanent guard: on a uniform tilted ramp the un-conditioned additive
    /// FBM fabricates many closed pits (noise out-slopes the gentle bed); the
    /// flow-conditioned FBM must fabricate far fewer, because the amplitude cap and
    /// downslope stretch keep the noise from reversing the descent. Falsifiable and
    /// resolution-cheap — the C-1 acceptance in miniature.
    #[test]
    fn flow_conditioning_suppresses_fabricated_pits() {
        let n = 48;
        let mut coarse = GridF32::new(n, n, 0.0);
        for j in 0..n {
            for i in 0..n {
                // a gentle land ramp well above sea (0.5): 0.6 → 0.9 across x.
                coarse.set(i, j, 0.6 + i as f32 / n as f32 * 0.3);
            }
        }
        let seed = WorldSeed::new(42);
        let base =
            FbmUpscaleConfig { target_size: 384, amplitude_base: 0.04, ..Default::default() };
        let off = upscale_with_fbm(&coarse, 0.5, &seed, &base).heightmap;
        let on = upscale_with_fbm(
            &coarse,
            0.5,
            &seed,
            &FbmUpscaleConfig { flow_conditioning: 0.1, ..base.clone() },
        )
        .heightmap;
        let (m_off, m_on) = (local_minima(&off), local_minima(&on));
        assert!(
            m_on * 4 < m_off,
            "conditioning must cut fabricated pits >4×: off {m_off}, on {m_on}"
        );
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

#[cfg(test)]
mod production_config_guards {
    use super::*;

    fn opts() -> ProductionHdOpts {
        ProductionHdOpts {
            target_size: 2048,
            domain_km: 400.0,
            depth_scale_m: 5000.0,
            sample_origin: [0.09375, 0.578125],
            sample_size: 1.0,
            amplitude_base: 0.04,
            mfd_p: 2.0,
            lithology: Default::default(),
            fracture: Default::default(),
        }
    }

    /// NON-REGRESSION: `production_hd_config` must return the SHIPPED recipe. It exists so
    /// the viz and the benches consume ONE config; if either side drifts, this breaks.
    /// Guards the seventh bench/production divergence (ADR "The DEAD KNOB").
    #[test]
    fn production_upscale_config_is_the_shipped_one() {
        let c = production_hd_config(&opts());
        // relief-v3: stream-power ON, droplets OFF (they collapse the SP valleys).
        let sp = c.stream_power.as_ref().expect("relief-v3 ships stream-power incision");
        assert!(c.erosion.is_none(), "droplets must be OFF on the shipped path");
        assert_eq!(sp.mfd_exponent, Some(2.0), "relief-v3 is MFD");
        assert_eq!(sp.iterations, 2, "production runs 2 incision iterations");
        // C-1 conditioning ON — this is what makes the amplitude cap bind (below).
        assert!(c.flow_conditioning > 0.0, "C-1 flow conditioning ships ON");
        // The framing is the whole torus, rolled — never a crop.
        assert_eq!(c.sample_size, 1.0, "the domain IS the map: no crop");
        assert_eq!(c.amplitude_base, 0.04, "the shipped amplitude value (inert, see below)");
    }

    /// THE INERTNESS, ASSERTED. In production the C-1 relief budget is smaller than the
    /// configured level at every plausible slope, so `amplitude_base` changes nothing —
    /// proven byte-identical at 4× in `tests/amplitude_anomaly.rs`. This pins the fact: if
    /// β, the divisor or the base ever make the LEVEL bind again, the knob comes back to
    /// life and this test tells whoever changed it.
    #[test]
    fn amplitude_cap_binds_everywhere_in_production() {
        let c = production_hd_config(&opts());
        // `flow_budget_divisor` in the loop is nscale·S; use the same order of magnitude by
        // sweeping the ratio over a wide range of slopes and divisors.
        for divisor in [1.0f64, 4.0, 16.0, 64.0] {
            for slope in [0.0f64, 0.001, 0.01, 0.05, 0.2, 0.5, 1.0] {
                let terms = AmplitudeTerms {
                    base: c.amplitude_base * (1.0 + slope * c.amplitude_slope_factor),
                    cap: Some(c.flow_conditioning * slope / divisor),
                };
                assert!(
                    terms.cap_binds() || terms.base == 0.0 || slope == 0.0,
                    "the relief budget must bind (slope {slope}, divisor {divisor}): \
                     base {} vs cap {:?}",
                    terms.base,
                    terms.cap
                );
                assert_eq!(terms.effective(), terms.base.min(terms.cap.unwrap()));
            }
        }
    }

    /// With conditioning OFF (the pre-C-1 regime) the configured level DOES drive the
    /// terrain — which is why Finding 5's amplitude sweep was valid when it was made.
    #[test]
    fn without_conditioning_the_base_term_drives() {
        let t = AmplitudeTerms { base: 0.16, cap: None };
        assert_eq!(t.effective(), 0.16);
        assert!(!t.cap_binds());
    }
}
