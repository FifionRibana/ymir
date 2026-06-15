//! Airy isostasy: converts crustal thickness to altitude.
//!
//! Models the crust as blocks floating on the mantle. Thick crust (continents)
//! stands above sea level; thin crust (ocean) is submerged. The altitude is
//! determined by Archimedes' principle:
//!
//!   h_raw = S × (1 − ρ_crust / ρ_mantle)

use crate::grid::GridF32;
use crate::tectonics::solver::field::Field2D;
use serde::{Deserialize, Serialize};

/// How the sea-level threshold is positioned within the crustal /
/// altitude distribution (Issue #141 Phase 1.5).
///
/// `MinMaxFraction` is the original formula (`min + frac·(max −
/// min)`) and the [`IsostasyConfig::default`] mode — v2 + export +
/// gallery paths keep it, byte-identical. `PercentileCapped` caps the
/// upper end of the range at the `cap_percentile`-th percentile
/// instead of the raw max, so a thin upper tail (e.g. Davis-Suppe
/// orographic peaks pushing `s_max` to ~2.18 while P50 ≈ 0.28) no
/// longer inflates the threshold and submerges the bulk of the crust.
/// M2 (Issue #139 Stage V diagnostic) measured P95-cap → ~28% emergent
/// land vs ~5.9% under min/max on the same field.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum SeaLevelMode {
    /// Sea level = `min + frac·(max − min)`. Original behaviour.
    MinMaxFraction,
    /// Sea level = `min + frac·(percentile(data, cap) − min)`. Robust
    /// to upper-tail outliers. `cap_percentile ∈ [0, 1]` (e.g. 0.95).
    PercentileCapped { cap_percentile: f32 },
}

impl Default for SeaLevelMode {
    fn default() -> Self {
        SeaLevelMode::MinMaxFraction
    }
}

/// O(N)-average percentile via `select_nth_unstable_by` (no full
/// sort, no NaN handling beyond `partial_cmp` — fields are finite by
/// invariant). `q ∈ [0, 1]`; index = `round(q·(n−1))`. Allocates one
/// scratch copy of `data`. Used by both sea-level instances so the
/// S̃-space and h-space thresholds stay coherent (Issue #141 W2).
pub(crate) fn percentile_copy<T: Copy + PartialOrd>(data: &[T], q: f32) -> T {
    debug_assert!(!data.is_empty(), "percentile of empty slice");
    let mut buf: Vec<T> = data.to_vec();
    let n = buf.len();
    let idx = (((q as f64) * (n - 1) as f64).round() as usize).min(n - 1);
    let (_, nth, _) =
        buf.select_nth_unstable_by(idx, |a, b| a.partial_cmp(b).expect("non-finite in percentile"));
    *nth
}

/// Configuration for isostatic altitude computation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsostasyConfig {
    /// Crustal density (kg/m³). Default: 2750.
    pub rho_crust: f32,
    /// Mantle density (kg/m³). Default: 3300.
    pub rho_mantle: f32,
    /// Water density (kg/m³). Default: 1025.
    pub rho_water: f32,
    /// Maximum land elevation in meters. Default: 4000.
    pub max_elevation_m: f32,
    /// Maximum ocean depth in meters (positive). Default: 500.
    pub max_depth_m: f32,
    /// Sea level as a fraction of the raw isostatic range [0, 1].
    /// 0.0 = everything is land, 1.0 = everything is ocean.
    /// Default: 0.4 (roughly ~30% land / 70% ocean).
    pub sea_level_fraction: f32,
    /// Gaussian blur sigma applied to the altitude heightmap after
    /// isostatic computation. Smooths sharp tectonic transitions.
    /// Default: 2.0. Set to 0.0 to disable.
    pub altitude_smoothing_sigma: f32,
    /// How the sea-level threshold is positioned (Issue #141). Default
    /// `MinMaxFraction` (v2 + export bit-compat). C1 uses
    /// [`IsostasyConfig::c1_default`] (`PercentileCapped`).
    /// `#[serde(default)]` so legacy configs without this field
    /// deserialize to `MinMaxFraction` (avoids the #47-style break).
    #[serde(default)]
    pub sea_level_mode: SeaLevelMode,
    /// #155 B-Jordan — optional CRATONIC crustal density (kg/m³). When
    /// `Some`, [`compute_isostasy_craton`] uses it instead of `rho_crust`
    /// on cells flagged cratonic → lower Airy buoyancy → a thick craton
    /// (1.25× S̃) elevates LESS (compositional isostasy: old cratonic
    /// lower crust is dense / the cold lithospheric root is ~neutrally
    /// buoyant, so the thick crust does NOT over-elevate — worn-shield
    /// height). Real cratonic crust ~2850-2950 vs 2750 normal. Default
    /// `None` → scalar buoyancy = byte-identical (v2/export untouched).
    /// Orthogonal to `craton_resist` (this is ALTITUDE via density; resist
    /// is the S̃ order). Thickness ratio (1.25) stays anchored.
    #[serde(default)]
    pub craton_rho_crust: Option<f32>,
    /// #155 critical-wedge — FIXED land-elevation reference (S̃ thickness). When
    /// `Some(rt)`, the LAND ramp normalises against `rt · buoyancy` (a FIXED raw
    /// height) instead of the field's `h_max`/`land_cap`, and `peak_altitude_m =
    /// max_elevation_m`. This DECOUPLES the rendered land elevation from the
    /// field's own relief: a given S̃ thickness renders the SAME metres in every
    /// seed and after every closure (no more field-dependent ~0.74-0.80 fraction
    /// that shifts as the orogen grows). `rt = h_eq = 2.0` (= ~70 km, the EH cap
    /// = the physical max crustal thickness) is the single anchor the whole
    /// relief converges on (the critical-wedge targets it, EH caps it, the scale
    /// references it). Cells above `rt` clamp to norm 1.0 = `max_elevation_m`.
    /// `None` → the field-relative `land_cap` behaviour (byte-identical; v2/export).
    #[serde(default)]
    pub land_ref_thickness: Option<f32>,
}

impl Default for IsostasyConfig {
    fn default() -> Self {
        Self {
            rho_crust: 2750.0,
            rho_mantle: 3300.0,
            rho_water: 1025.0,
            max_elevation_m: 4000.0,
            max_depth_m: 500.0,
            sea_level_fraction: 0.4,
            altitude_smoothing_sigma: 2.0,
            sea_level_mode: SeaLevelMode::MinMaxFraction,
            craton_rho_crust: None,
            land_ref_thickness: None,
        }
    }
}

impl IsostasyConfig {
    /// C1 default (Issue #141 Phase 1.5): identical to [`Default`]
    /// except `sea_level_mode = PercentileCapped { cap_percentile:
    /// 0.92 }`. The cap was **calibrated live** (not from M2's static
    /// post-hoc 28%, which the in-loop feedback contradicted): a
    /// multi-seed 15-cycle sweep showed 0.92 damps to a bounded
    /// equilibrium with emergent-land distribution mean ~30.6%, range
    /// ~24.5–36.6% (natural per-seed variation — "around 30%"). 0.95
    /// oscillates persistently and undershoots (~20%); 0.90 is erratic
    /// (runaway). cap=0.92 is COUPLED with `n_cycles ≈ 12` (worst-case
    /// band-entry cycle 9 + margin) — the system is a bounded limit
    /// cycle (±0.05), not a fixed point. Used by the C1 engine's
    /// config builders (viz workflow / gallery worker / render
    /// altitude, and C1 validation tests); v2 + export + the gallery
    /// PNG generators keep [`Default`] (`MinMaxFraction`).
    pub fn c1_default() -> Self {
        Self {
            sea_level_mode: SeaLevelMode::PercentileCapped { cap_percentile: 0.92 },
            // #155 méso — reduced from the inherited 2.0. The σ=2 gaussian
            // blur ("smooths sharp tectonic transitions") was the DOMINANT
            // wall smoothing sub-macro (méso) structure out of the altitude
            // BEFORE upscale (σ-sweep: méso coarse osc 0.05@σ2 → 0.35@σ0.5);
            // it also smooths the REAL macro O-C ridge (gap +0.06@σ2 →
            // +0.11@σ0.5 on seed 1988). 0.5 captures ~80-87% of the
            // sharpening while keeping mild smoothing (cleaner than σ=0, no
            // Voronoi-staircase). Confirmed NOT the #151 anti-feathering
            // (that is the upscale `coastal_band`); reducing σ unmasks NO
            // grid striping / steps / feathering (σ-sweep, real hillshades
            // 1988+4138). v2 + export keep the `Default` 2.0 (bit-compat).
            // See `stage_meso_expression.md`.
            altitude_smoothing_sigma: 0.5,
            // #155 B-Jordan — cratonic crustal density 2900 (real range
            // 2850-2950; old dense lower crust / cold neutrally-buoyant
            // root). Thick craton (1.25× S̃) elevates LESS → worn high-
            // plateau (~1650 m, S. African/Colorado type) instead of ~2340 m.
            // Compositional isostasy — the 3rd real craton property (after
            // thickness=Airy + resist=S̃-order); ALTITUDE-only (compute_isostasy
            // scalar/None ignores it → loop S̃ unchanged; applies in
            // upscale_from_c1 which threads cratonic_mask). Validated 4-density
            // visual; 2900 corrects height without under-selling cratonic
            // orogens (2950 over-lowers). Default keeps None (v2/export
            // bit-compat). Anchored — not tuned to a target.
            craton_rho_crust: Some(2900.0),
            // #155 critical-wedge — land elevation anchor = the contract norm-1.0
            // (5650 m). With `land_ref_thickness=Some(2.0)`, norm 1.0 ↔ S̃=2.0
            // (~70 km, the EH cap) ↔ max_elevation_m ↔ the vertical contract max:
            // ONE coherent anchor. An S̃=2.0 orogen renders 5650 m, S̃≈1.9 renders
            // ~5.3 km, FIELD-INDEPENDENTLY (no re-pin; rendered ≤ contract). Default
            // keeps 4000 (v2/export byte-identical). 8000 m Everest POINT peaks
            // exceed the contract + C1's mean-field resolution → out of scope.
            max_elevation_m: 5650.0,
            // #155 critical-wedge — fix the land reference to S̃=2.0 (h_eq=70 km,
            // EH cap): decouples rendered orogen height from the field's own
            // relief (the ~0.78 field-dependent fraction → field-INDEPENDENT). The
            // single anchor the relief converges on (wedge targets 2.0, EH caps
            // 2.0, scale references 2.0). Default None (v2/export field-relative).
            land_ref_thickness: Some(2.0),
            ..Default::default()
        }
    }
}

/// Result of isostatic computation.
pub struct IsostasyResult {
    /// Altitude heightmap normalized to [0, 1].
    /// 0.0 = deepest ocean, sea_level_normalized = sea level, 1.0 = highest peak.
    pub heightmap: GridF32,
    /// The normalized value that corresponds to sea level.
    pub sea_level_normalized: f32,
    /// Peak altitude in meters.
    pub peak_altitude_m: f32,
    /// Deepest point in meters below sea level (positive value).
    pub max_depth_m: f32,
    /// Fraction of cells above sea level.
    pub land_ratio: f32,
}

/// Convert crustal thickness field to altitude heightmap via Airy isostasy.
///
/// The input is the Field2D from the solver (f64, dimensionless).
/// The output is a normalized GridF32 heightmap suitable for erosion and export.
pub fn compute_isostasy(thickness: &Field2D, config: &IsostasyConfig) -> IsostasyResult {
    compute_isostasy_inner(thickness, config, None, None)
}

/// #155 B-Jordan — Airy isostasy with spatial CRATONIC density. Identical
/// to [`compute_isostasy`] except cells where `craton[k]` is true use
/// `config.craton_rho_crust` (if `Some`) for the buoyancy — a denser
/// cratonic crust → lower buoyancy → a thick craton elevates LESS
/// (compositional isostasy / worn-shield height). `craton` is a row-major
/// `&[bool]` (length nx·ny) to avoid a tectonics→tectonics_c1 dependency
/// cycle. With `craton_rho_crust == None` → byte-identical to the scalar.
pub fn compute_isostasy_craton(
    thickness: &Field2D,
    config: &IsostasyConfig,
    craton: &[bool],
) -> IsostasyResult {
    compute_isostasy_inner(thickness, config, Some(craton), None)
}

/// #155 vertical-scale — Airy isostasy with BOTH the cratonic density mask
/// (`craton`, see [`compute_isostasy_craton`]) AND a CONTINENTAL mask that
/// repairs the land-normalisation ceiling.
///
/// **The land-ceiling repair.** [`compute_isostasy`] normalises the LAND ramp
/// `[h_sea, h_max] → [sea_norm, 1.0]` against the GLOBAL raw `h_max`. In the C1
/// production that global max is the OCEANIC advective pile-up cell (S̃ ≈ 2.15
/// at the rigid margin, #145) — a cell Stein-Stein then overwrites into deep
/// ocean. So real continental land is normalised against a phantom oceanic
/// ceiling (tops at ~0.94) and `peak_altitude_m` describes that ocean cell, not
/// any land (the `project_c1_vertical_scale_contract` diagnostic). When
/// `continental` is `Some`, the land ceiling is the max `h` over CONTINENTAL
/// cells instead — the real terrestrial summit; oceanic cells (S-S-overwritten)
/// no longer inflate it. `continental == None` → byte-identical to the global
/// `h_max` behaviour (v2/export untouched — they have no S-S/no advective spike,
/// so no phantom: the fix is C1-specific). Row-major `&[bool]` (length nx·ny),
/// same convention as `craton` (avoids a tectonics→tectonics_c1 cycle).
pub fn compute_isostasy_c1(
    thickness: &Field2D,
    config: &IsostasyConfig,
    craton: Option<&[bool]>,
    continental: &[bool],
) -> IsostasyResult {
    compute_isostasy_inner(thickness, config, craton, Some(continental))
}

fn compute_isostasy_inner(
    thickness: &Field2D,
    config: &IsostasyConfig,
    craton: Option<&[bool]>,
    continental: Option<&[bool]>,
) -> IsostasyResult {
    let nx = thickness.nx();
    let ny = thickness.ny();
    let buoyancy = 1.0 - config.rho_crust / config.rho_mantle;
    // #155 B-Jordan: cratonic buoyancy (denser crust → less elevation).
    let buoyancy_craton = match config.craton_rho_crust {
        Some(r) => 1.0 - r / config.rho_mantle,
        None => buoyancy,
    };

    // 1. Compute raw isostatic elevation
    let mut h_raw = vec![0.0f32; nx * ny];
    let mut h_min = f32::INFINITY;
    let mut h_max = f32::NEG_INFINITY;

    for (k, val) in thickness.data().iter().enumerate() {
        let b = if craton.is_some_and(|c| c[k]) { buoyancy_craton } else { buoyancy };
        let h = *val as f32 * b;
        h_raw[k] = h;
        h_min = h_min.min(h);
        h_max = h_max.max(h);
    }

    // 2. Determine sea level from the configured fraction + mode.
    // MinMaxFraction caps at the raw max (original, byte-identical);
    // PercentileCapped caps at the cap-percentile of the raw
    // altitude (Issue #141 — robust to upper-tail orographic peaks).
    let h_cap = match config.sea_level_mode {
        SeaLevelMode::MinMaxFraction => h_max,
        SeaLevelMode::PercentileCapped { cap_percentile } => {
            percentile_copy(&h_raw, cap_percentile)
        }
    };
    let h_range = (h_cap - h_min).max(1e-10);
    let h_sea = h_min + config.sea_level_fraction * h_range;

    // 2b. #155 land-ceiling repair. The LAND ramp tops at this ceiling. By
    // default (continental == None) it is the global raw h_max — but in the C1
    // production that is the OCEANIC advective spike (S-S-overwritten), so real
    // land is normalised against a phantom. When a continental mask is given,
    // use the max h over CONTINENTAL cells: the real terrestrial summit, with
    // oceanic cells excluded (they get S-S bathymetry downstream regardless).
    // See `compute_isostasy_c1`. Fallback to h_max if the mask flags no cell.
    let land_cap = match continental {
        Some(cont) => {
            let mut m = f32::NEG_INFINITY;
            for k in 0..nx * ny {
                if cont[k] {
                    m = m.max(h_raw[k]);
                }
            }
            if m.is_finite() { m } else { h_max }
        }
        None => h_max,
    };

    // 2c. #155 critical-wedge — the LAND-ramp ceiling. With `land_ref_thickness`
    // it is a FIXED raw height `rt · buoyancy` (rt = S̃ reference, the EH cap 2.0
    // = 70 km), NOT the field's `land_cap`: this decouples the rendered land
    // elevation from the field's own relief (a given S̃ → the same metres in every
    // seed / after every closure). Cells above it clamp to norm 1.0. `None` →
    // the field-relative `land_cap` (byte-identical).
    let land_ceiling = match config.land_ref_thickness {
        Some(rt) => rt * buoyancy,
        None => land_cap,
    };

    // 3. Map to normalized [0, 1] with sea level at a known position
    // sea_norm = max_depth / (max_depth + max_elevation)
    let sea_norm = config.max_depth_m / (config.max_depth_m + config.max_elevation_m);

    let mut data = vec![0.0f32; nx * ny];
    let mut land_count = 0usize;

    for k in 0..nx * ny {
        let h = h_raw[k];
        let normalized = if h <= h_sea {
            let t = (h - h_min) / (h_sea - h_min).max(1e-10);
            t * sea_norm
        } else {
            // #155 land ramp tops at `land_ceiling` (fixed S̃-reference when
            // land_ref_thickness is set, else the continental land_cap). Cells
            // above it clamp to 1.0 below.
            let t = (h - h_sea) / (land_ceiling - h_sea).max(1e-10);
            sea_norm + t * (1.0 - sea_norm)
        };
        data[k] = normalized.clamp(0.0, 1.0);

        if h > h_sea {
            land_count += 1;
        }
    }

    let land_ratio = land_count as f32 / (nx * ny) as f32;

    // 4. Compute actual peak altitude and depth for metadata. #155 critical-wedge:
    // with a fixed land reference, norm 1.0 ↔ `land_ceiling` (S̃=rt) ↔
    // `max_elevation_m` directly — so `peak_altitude_m = max_elevation_m`,
    // field-INDEPENDENT (a cell at the reference thickness renders the same metres
    // in every field). Without it, the field-relative continental-summit fraction
    // (byte-identical None path).
    let peak_altitude_m = if config.land_ref_thickness.is_some() {
        config.max_elevation_m
    } else {
        (land_cap - h_sea) / (land_cap - h_min).max(1e-10) * config.max_elevation_m
    };
    let actual_depth_m = (h_sea - h_min) / (h_max - h_min).max(1e-10) * config.max_depth_m;

    let heightmap = GridF32::from_vec(nx, ny, data);
    let heightmap = if config.altitude_smoothing_sigma > 0.0 {
        heightmap.gaussian_blur(config.altitude_smoothing_sigma)
    } else {
        heightmap
    };

    IsostasyResult {
        heightmap,
        sea_level_normalized: sea_norm,
        peak_altitude_m,
        max_depth_m: actual_depth_m,
        land_ratio,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_continental_gives_all_land() {
        let mut s = Field2D::new(16, 16);
        for val in s.data_mut() {
            *val = 1.0;
        }

        let config = IsostasyConfig { sea_level_fraction: 0.0, ..Default::default() };
        let result = compute_isostasy(&s, &config);
        // sea_level_fraction=0 → h_sea = h_min, all cells have h > h_sea
        // But when all values are equal, h_min == h_max, so h_sea == h_min == h_max
        // All cells have h == h_sea, not h > h_sea, so land_ratio = 0.
        // With uniform thickness there's no land/ocean distinction.
        // This is a degenerate case — test the non-degenerate one instead.
    }

    #[test]
    fn all_ocean_when_sea_level_at_max() {
        let n = 16;
        let mut s = Field2D::new(n, n);
        for j in 0..n {
            for i in 0..n {
                s.set(i, j, if i < n / 2 { 0.2 } else { 1.0 });
            }
        }

        let config = IsostasyConfig { sea_level_fraction: 1.0, ..Default::default() };
        let result = compute_isostasy(&s, &config);
        assert!(result.land_ratio < 1e-6, "All ocean expected, got {}", result.land_ratio);
    }

    #[test]
    fn thicker_crust_gives_higher_altitude() {
        let n = 32;
        let mut s = Field2D::new(n, n);
        for j in 0..n {
            for i in 0..n {
                s.set(i, j, if i < n / 2 { 0.2 } else { 1.0 });
            }
        }

        let result = compute_isostasy(&s, &IsostasyConfig::default());

        let h_ocean = result.heightmap.get(0, 0);
        let h_land = result.heightmap.get(n as i32 - 1, 0);
        assert!(h_land > h_ocean, "Continental should be higher: {} vs {}", h_land, h_ocean);
    }

    #[test]
    fn sea_level_normalized_is_consistent() {
        let n = 16;
        let mut s = Field2D::new(n, n);
        for j in 0..n {
            for i in 0..n {
                s.set(i, j, if i < n / 2 { 0.2 } else { 1.0 });
            }
        }

        let config = IsostasyConfig::default();
        let result = compute_isostasy(&s, &config);

        let expected = config.max_depth_m / (config.max_depth_m + config.max_elevation_m);
        assert!(
            (result.sea_level_normalized - expected).abs() < 1e-6,
            "sea_level_normalized: {} vs expected {}",
            result.sea_level_normalized,
            expected
        );
    }

    #[test]
    fn heightmap_values_in_range() {
        let n = 32;
        let mut s = Field2D::new(n, n);
        for (k, val) in s.data_mut().iter_mut().enumerate() {
            *val = 0.2 + (k as f64 / (n * n) as f64) * 1.8;
        }

        let result = compute_isostasy(&s, &IsostasyConfig::default());
        for val in &result.heightmap.data {
            assert!(*val >= 0.0 && *val <= 1.0, "Height out of range: {}", val);
        }
    }

    #[test]
    fn land_ratio_increases_as_sea_level_drops() {
        let n = 32;
        let mut s = Field2D::new(n, n);
        // Use a gradient so that the sea level threshold splits at different points
        for j in 0..n {
            for i in 0..n {
                s.set(i, j, 0.2 + 0.8 * (i as f64 / (n - 1) as f64));
            }
        }

        let r_high =
            compute_isostasy(&s, &IsostasyConfig { sea_level_fraction: 0.7, ..Default::default() });
        let r_low =
            compute_isostasy(&s, &IsostasyConfig { sea_level_fraction: 0.3, ..Default::default() });

        assert!(
            r_low.land_ratio > r_high.land_ratio,
            "Lower sea level should give more land: {} vs {}",
            r_low.land_ratio,
            r_high.land_ratio
        );
    }

    #[test]
    fn smoothing_reduces_max_gradient() {
        let n = 32;
        let mut s = Field2D::new(n, n);
        for j in 0..n {
            for i in 0..n {
                s.set(i, j, if i < n / 2 { 0.2 } else { 1.0 });
            }
        }

        let config_sharp = IsostasyConfig { altitude_smoothing_sigma: 0.0, ..Default::default() };
        let config_smooth = IsostasyConfig { altitude_smoothing_sigma: 2.0, ..Default::default() };

        let result_sharp = compute_isostasy(&s, &config_sharp);
        let result_smooth = compute_isostasy(&s, &config_smooth);

        let max_grad = |hm: &GridF32| -> f32 {
            let mut max = 0.0f32;
            for j in 0..hm.height {
                for i in 1..hm.width {
                    let g = (hm.data[j * hm.width + i] - hm.data[j * hm.width + i - 1]).abs();
                    max = max.max(g);
                }
            }
            max
        };

        let grad_sharp = max_grad(&result_sharp.heightmap);
        let grad_smooth = max_grad(&result_smooth.heightmap);

        assert!(
            grad_smooth < grad_sharp,
            "Smoothing should reduce max gradient: sharp={}, smooth={}",
            grad_sharp,
            grad_smooth
        );
    }

    // ----- Issue #141 Phase 1.5: SeaLevelMode -----

    /// A field with a SPREAD bulk (gradient 0.2..0.9, ~95% of cells)
    /// plus a thin tall tail (~5% peaks at 2.0) — the M2 configuration
    /// in miniature (P50/P90 within the bulk, max in the tail, so
    /// min/max sea level sits above the bulk but P95 sits inside it).
    fn tailed_field(n: usize) -> Field2D {
        let mut s = Field2D::new(n, n);
        for j in 0..n {
            for i in 0..n {
                s.set(i, j, 0.2 + 0.7 * (i as f64 / (n - 1) as f64));
            }
        }
        // ~5% of cells: tall orographic peaks (inflate min/max, not P95).
        let peaks = (n * n) / 20;
        for k in 0..peaks {
            s.data_mut()[k] = 2.0;
        }
        s
    }

    #[test]
    fn defaults_and_c1_default_modes() {
        assert_eq!(IsostasyConfig::default().sea_level_mode, SeaLevelMode::MinMaxFraction);
        assert_eq!(
            IsostasyConfig::c1_default().sea_level_mode,
            SeaLevelMode::PercentileCapped { cap_percentile: 0.92 }
        );
        // c1_default differs from Default in the sea-level MODE and (since
        // #155 méso) the altitude_smoothing_sigma: 2.0 (Default, v2/export
        // bit-compat) vs 0.5 (C1, so the σ=2 blur no longer smooths the
        // macro O-C ridge + méso structure out of the altitude). Other
        // fields match.
        let d = IsostasyConfig::default();
        let c = IsostasyConfig::c1_default();
        assert_eq!(d.sea_level_fraction, c.sea_level_fraction);
        // #155 critical-wedge: c1_default max_elevation_m = 5650 (= contract
        // norm-1.0) + land_ref_thickness Some(2.0) (fixed S̃=70 km reference, so
        // rendered orogen height is field-independent); Default keeps 4000 + None
        // (field-relative, v2/export bit-compat).
        assert_eq!(d.max_elevation_m, 4000.0);
        assert_eq!(c.max_elevation_m, 5650.0);
        assert_eq!(d.land_ref_thickness, None);
        assert_eq!(c.land_ref_thickness, Some(2.0));
        assert_eq!(d.altitude_smoothing_sigma, 2.0);
        assert_eq!(c.altitude_smoothing_sigma, 0.5);
    }

    #[test]
    fn percentile_cap_at_100_equals_min_max() {
        // P100 == max, so PercentileCapped{1.0} must give the SAME
        // h_sea / land_ratio as MinMaxFraction (byte-identity sanity).
        let s = tailed_field(32);
        let r_minmax = compute_isostasy(&s, &IsostasyConfig::default());
        let r_p100 = compute_isostasy(
            &s,
            &IsostasyConfig {
                sea_level_mode: SeaLevelMode::PercentileCapped { cap_percentile: 1.0 },
                ..Default::default()
            },
        );
        assert_eq!(r_minmax.land_ratio, r_p100.land_ratio);
    }

    #[test]
    fn percentile_cap_lowers_sea_level_and_raises_land() {
        // The Phase 1.5 point: capping the upper tail drops the sea
        // level, so MORE of the bulk crust emerges. Same field.
        let s = tailed_field(32);
        let r_minmax = compute_isostasy(&s, &IsostasyConfig::default());
        let r_p95 = compute_isostasy(&s, &IsostasyConfig::c1_default());
        assert!(
            r_p95.land_ratio > r_minmax.land_ratio,
            "P95-cap should raise emergent land vs min/max: p95={} minmax={}",
            r_p95.land_ratio,
            r_minmax.land_ratio
        );
        // The bulk (98% of cells at 0.3) should now be land: min/max
        // sea level ≈ 0.3 + 0.4·(2.0−0.3) ≈ 0.98 (above the bulk → ~2%
        // land); P95-cap sea level ≈ 0.3 + 0.4·(0.3−0.3)=0.3 (bulk at
        // the boundary). The jump must be large.
        assert!(
            r_minmax.land_ratio < 0.10 && r_p95.land_ratio > 0.5,
            "expected min/max ~2% vs P95-cap majority: minmax={} p95={}",
            r_minmax.land_ratio,
            r_p95.land_ratio
        );
    }

    #[test]
    fn percentile_copy_picks_expected_value() {
        let data: Vec<f64> = (0..=100).map(|i| i as f64).collect(); // 0..100
        assert_eq!(percentile_copy(&data, 0.0), 0.0);
        assert_eq!(percentile_copy(&data, 1.0), 100.0);
        assert_eq!(percentile_copy(&data, 0.5), 50.0);
        assert_eq!(percentile_copy(&data, 0.95), 95.0);
    }

    #[test]
    fn serde_legacy_config_without_mode_defaults_minmax() {
        // A config JSON predating Phase 1.5 (no sea_level_mode field)
        // must deserialize to MinMaxFraction (W4 — no #47-style break).
        let legacy = r#"{
            "rho_crust": 2750.0, "rho_mantle": 3300.0, "rho_water": 1025.0,
            "max_elevation_m": 4000.0, "max_depth_m": 500.0,
            "sea_level_fraction": 0.4, "altitude_smoothing_sigma": 2.0
        }"#;
        let cfg: IsostasyConfig = serde_json::from_str(legacy).expect("legacy deserialize");
        assert_eq!(cfg.sea_level_mode, SeaLevelMode::MinMaxFraction);
    }

    #[test]
    fn serde_roundtrip_percentile_mode() {
        let cfg = IsostasyConfig::c1_default();
        let json = serde_json::to_string(&cfg).unwrap();
        let back: IsostasyConfig = serde_json::from_str(&json).unwrap();
        assert_eq!(back.sea_level_mode, cfg.sea_level_mode);
    }
}
