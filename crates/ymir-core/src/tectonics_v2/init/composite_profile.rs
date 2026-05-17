//! Step 12 R7.A.2 — composite continental S̃ initialisation.
//!
//! Combines [`super::radial_profile`] (radial dome) with the orogenic
//! ridge formula from [`super::orogenic_profile`] additively, with a
//! cap. Each continental plate gets a smooth radial-symmetric dome
//! peaked at its interior plus a linear ridge along its PCA principal
//! axis. The two contributions superpose to produce the
//! "continent-wide elevated terrain + linear chain on top" morphology
//! we associate with e.g. South America (continental dome + Andes
//! along the western coast).
//!
//! ## Algorithm (R7.A.2.1 §3)
//!
//! 1. Run [`super::radial_profile::build`] with the user's
//!    `continental_value`, `oceanic_value`, `profile_shape` → get the
//!    **dome field** with the Step 13 RadialProfile semantics (oceanic
//!    cells flat at `oceanic_value`, continental cells smooth from
//!    boundary → interior).
//! 2. Run the same per-plate geometry pipeline as Orogenic-seul
//!    (periodic centroid + PCA principal axis + `L_plate` from BFS).
//!    Reuses [`super::orogenic_profile::PlateAccum`] verbatim so the
//!    geometry is byte-identical to the Orogenic-seul mode.
//! 3. For each continental cell, compute `ridge_amount ∈ [0, 1]`
//!    from the orogenic formula with an optional
//!    `offset_along_axis_ratio` (R7.A.2 addition; the standalone
//!    Orogenic mode is hardcoded to offset 0):
//!    ```text
//!    d_along_adjusted = d_along − offset · half_length
//!    t_along          = (1 − |d_along_adjusted|/half_length).clamp(0, 1)
//!    long_mod         = 3·t_along² − 2·t_along³  (smoothstep)
//!    trans            = exp(−(d_perp / width_sigma)²)
//!    ridge_amount     = long_mod · trans
//!    ```
//! 4. Combine additively with cap:
//!    ```text
//!    S̃ = clamp(dome[i,j] + (peak − base) · ridge_amount, 0, cap)
//!    ```
//!    for continental cells. Oceanic cells get `oceanic_value`
//!    uniform (bypass the formula).
//!
//! ## Cap conventions
//!
//! - [`CompositeCap::UsePeakOrogenic`] — cap at
//!   `orogenic_ridge.peak_value`. Composite max equals Orogenic-seul
//!   max; ~7 % of the absolute crest amplitude is clipped at the
//!   exact axis × centroid intersection (~1-2 cells per plate).
//! - [`CompositeCap::Fixed { value }`] — explicit user cap.
//!
//! ## Degenerate cases
//!
//! Same as Orogenic-seul:
//! - Plates with < 5 cells / rank-1 PCA fall back to
//!   `OrogenicOrientation::Fixed { angle_rad: 0.0 }`.
//! - Plates with `L_plate ≤ 0` (single-plate-on-torus / fully-on-
//!   boundary): ridge degenerates to amplitude 0 ⇒ `S̃ = dome` for
//!   that plate. Same convention as Orogenic-seul and RadialProfile.
//!
//! ## R7.A.1 byte-equal preservation
//!
//! This module only **reads** the `pub(crate)` `PlateAccum` API of
//! `orogenic_profile.rs`. It never invokes `orogenic_profile::build`
//! nor mutates its internals. Existing R7.A.1 `InitMode::Orogenic`
//! callers are unaffected; the orogenic test suite still passes
//! byte-for-byte.

use serde::{Deserialize, Serialize};

use super::super::boundaries::PlateType;
use super::super::field::Field2D;
use super::super::voronoi::compute_dist_to_inter_plate_boundary;
use super::orogenic_profile::{min_image_delta, OrogenicOrientation, PlateAccum};
use super::radial_profile::{self, ProfileShape};
use super::PlateInitData;

/// Radial dome parameters for the composite mode. Mirrors the subset
/// of [`super::radial_profile`] inputs the composite actually uses
/// (continental peak + profile shape). `oceanic_value` is shared at
/// the `Composite` variant top level since both the dome and the
/// ridge agree on the same flat oceanic baseline.
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq)]
pub struct CompositeRadialParams {
    pub continental_value: f64,
    pub profile_shape: ProfileShape,
}

/// Orogenic ridge parameters for the composite mode. Same fields as
/// [`super::orogenic_profile`] except that `oceanic_value` lives at
/// the `Composite` top level (shared with the dome), and a new
/// `offset_along_axis_ratio` displaces the ridge along the PCA axis
/// (R7.A.2 Q.f.3).
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq)]
pub struct CompositeOrogenicRidgeParams {
    pub peak_value: f64,
    pub base_continental_value: f64,
    pub half_length_ratio: f64,
    pub width_sigma_ratio: f64,
    pub orientation: OrogenicOrientation,
    /// R7.A.2 Q.f.3 — shift the ridge along its PCA principal axis
    /// by `offset · half_length` from the centroid. Default `0.0`
    /// keeps the ridge centred on the centroid (R7.A.1 behaviour).
    /// Range `[-1, 1]`: at `±1` the ridge sits at the longitudinal
    /// boundary of the plateau; beyond that the ridge effectively
    /// disappears (long_mod clamped to 0).
    pub offset_along_axis_ratio: f64,
}

/// Composite cap mode (R7.A.2 Q.f.2). `UsePeakOrogenic` (default)
/// makes the composite peak match the Orogenic-seul peak; `Fixed`
/// is a user override for experimentation.
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum CompositeCap {
    UsePeakOrogenic,
    Fixed { value: f64 },
}

impl Default for CompositeCap {
    fn default() -> Self {
        CompositeCap::UsePeakOrogenic
    }
}

/// Default `continental_value` for the composite dome. Inherits the
/// Step 13 RadialProfile default (0.95).
pub const COMPOSITE_RADIAL_CONTINENTAL_DEFAULT: f64 = 0.95;

/// Build the S̃ field for `InitMode::Composite`. See the module
/// docstring for the algorithm and R7.A.2.1 spec doc under
/// `docs/reports/step12_r7_a_composite_profile/` for the full
/// rationale.
///
/// # Panics
///
/// - Any of `radial.continental_value`, `orogenic.peak_value`,
///   `orogenic.base_continental_value`, `oceanic_value` is non-finite.
/// - `orogenic.half_length_ratio <= 0` or
///   `orogenic.width_sigma_ratio <= 0` (degenerate ridge).
pub fn build(
    nx: usize,
    ny: usize,
    p: &PlateInitData<'_>,
    radial: CompositeRadialParams,
    orogenic_ridge: CompositeOrogenicRidgeParams,
    oceanic_value: f64,
    cap: CompositeCap,
) -> Field2D {
    assert!(
        radial.continental_value.is_finite()
            && orogenic_ridge.peak_value.is_finite()
            && orogenic_ridge.base_continental_value.is_finite()
            && oceanic_value.is_finite(),
        "composite build: continental/peak/base/oceanic must be finite",
    );
    assert!(
        orogenic_ridge.half_length_ratio > 0.0
            && orogenic_ridge.half_length_ratio.is_finite(),
        "composite build: half_length_ratio must be > 0",
    );
    assert!(
        orogenic_ridge.width_sigma_ratio > 0.0
            && orogenic_ridge.width_sigma_ratio.is_finite(),
        "composite build: width_sigma_ratio must be > 0",
    );

    // ---- 1. Radial dome via the existing Step 13 path. ----
    let dome = radial_profile::build(
        nx,
        ny,
        p,
        radial.continental_value,
        oceanic_value,
        radial.profile_shape,
    );

    // ---- 2. Per-plate geometry (same pipeline as Orogenic-seul). ----
    let mut plate_meta: Vec<PlateAccum> = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let pid = p.plate_id.get(i, j) as usize;
            if pid >= plate_meta.len() {
                plate_meta.resize_with(pid + 1, PlateAccum::new);
            }
            plate_meta[pid].accumulate_pos(i, j, nx, ny);
        }
    }
    for acc in plate_meta.iter_mut() {
        acc.finalise_centroid(nx, ny);
    }
    for j in 0..ny {
        for i in 0..nx {
            let pid = p.plate_id.get(i, j) as usize;
            plate_meta[pid].accumulate_covariance(i, j, nx, ny);
        }
    }
    let bfs = compute_dist_to_inter_plate_boundary(nx, ny, p.plate_id);
    let mut plate_l: Vec<f64> = vec![0.0; plate_meta.len()];
    for j in 0..ny {
        for i in 0..nx {
            let pid = p.plate_id.get(i, j) as usize;
            let d = bfs.distance.get(i, j);
            if d.is_finite() && d > plate_l[pid] {
                plate_l[pid] = d;
            }
        }
    }
    // Resolve PCA axes (or fallback per orientation).
    let fallback_axis = match orogenic_ridge.orientation {
        OrogenicOrientation::PlateMainAxisPca => (1.0, 0.0),
        OrogenicOrientation::Fixed { angle_rad } => (angle_rad.cos(), angle_rad.sin()),
    };
    let mut plate_axis: Vec<(f64, f64)> = vec![(1.0, 0.0); plate_meta.len()];
    for (pid, acc) in plate_meta.iter().enumerate() {
        plate_axis[pid] = match orogenic_ridge.orientation {
            OrogenicOrientation::Fixed { angle_rad } => (angle_rad.cos(), angle_rad.sin()),
            OrogenicOrientation::PlateMainAxisPca => {
                acc.principal_axis().unwrap_or(fallback_axis)
            }
        };
    }

    // ---- 3. Compose: dome + (peak − base) · ridge_amount, capped. ----
    let cap_value = match cap {
        CompositeCap::UsePeakOrogenic => orogenic_ridge.peak_value,
        CompositeCap::Fixed { value } => value,
    };
    let ridge_amplitude =
        orogenic_ridge.peak_value - orogenic_ridge.base_continental_value;
    let nxf = nx as f64;
    let nyf = ny as f64;
    let mut s = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            // Oceanic cells: flat, no dome adjustment, no ridge.
            if matches!(p.plate_type.get(i, j), PlateType::Oceanic) {
                s.set(i, j, oceanic_value);
                continue;
            }
            let dome_value = dome.get(i, j);
            let pid = p.plate_id.get(i, j) as usize;
            let l_plate = plate_l[pid];
            if l_plate <= 0.0 {
                // Same degenerate convention as Orogenic-seul: no
                // ridge scale ⇒ dome-only value.
                s.set(i, j, dome_value);
                continue;
            }
            let half_length = (orogenic_ridge.half_length_ratio * l_plate).max(1e-12);
            let width_sigma = (orogenic_ridge.width_sigma_ratio * l_plate).max(1e-12);
            let acc = &plate_meta[pid];
            let (cx, cy) = match acc.centroid() {
                Some(c) => c,
                None => {
                    s.set(i, j, dome_value);
                    continue;
                }
            };
            let cell_x = i as f64 + 0.5;
            let cell_y = j as f64 + 0.5;
            let dx = min_image_delta(cell_x - cx, nxf);
            let dy = min_image_delta(cell_y - cy, nyf);
            let (ux, uy) = plate_axis[pid];
            let d_along_raw = dx * ux + dy * uy;
            let d_perp = -dx * uy + dy * ux;
            // R7.A.2 Q.f.3 — shift the ridge along the axis.
            let d_along = d_along_raw - orogenic_ridge.offset_along_axis_ratio * half_length;
            let t_along = (1.0 - d_along.abs() / half_length).clamp(0.0, 1.0);
            let long_mod = t_along * t_along * (3.0 - 2.0 * t_along);
            let dp_over_sigma = d_perp / width_sigma;
            let trans = (-(dp_over_sigma * dp_over_sigma)).exp();
            let ridge_amount = long_mod * trans;
            let composite = dome_value + ridge_amplitude * ridge_amount;
            s.set(i, j, composite.min(cap_value).max(0.0));
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics_v2::voronoi::{generate_voronoi, VoronoiConfig};

    fn build_plates(
        nx: usize,
        ny: usize,
        seed: u64,
        num_plates: usize,
        continental_ratio: f64,
    ) -> crate::tectonics_v2::voronoi::VoronoiPlates {
        let cfg = VoronoiConfig { num_plates, continental_ratio };
        generate_voronoi(nx, ny, &cfg, seed)
    }

    fn make_init_data<'a>(
        plates: &'a crate::tectonics_v2::voronoi::VoronoiPlates,
    ) -> PlateInitData<'a> {
        PlateInitData {
            plate_id: &plates.plate_id,
            plate_type: &plates.plate_type,
            seed_coords: Some(&plates.seed_coords),
        }
    }

    fn default_radial() -> CompositeRadialParams {
        CompositeRadialParams {
            continental_value: 0.95,
            profile_shape: ProfileShape::Smoothstep,
        }
    }

    fn default_orogenic_ridge() -> CompositeOrogenicRidgeParams {
        CompositeOrogenicRidgeParams {
            peak_value: 1.20,
            base_continental_value: 0.85,
            half_length_ratio: 0.40,
            width_sigma_ratio: 0.10,
            orientation: OrogenicOrientation::PlateMainAxisPca,
            offset_along_axis_ratio: 0.0,
        }
    }

    /// Acceptance #1 — every oceanic cell holds exactly
    /// `oceanic_value` regardless of dome / ridge configuration.
    #[test]
    fn composite_oceanic_uniform() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);
        let s = build(
            nx, ny, &p,
            default_radial(), default_orogenic_ridge(),
            0.20,
            CompositeCap::UsePeakOrogenic,
        );
        for j in 0..ny {
            for i in 0..nx {
                if matches!(plates.plate_type.get(i, j), PlateType::Oceanic) {
                    let v = s.get(i, j);
                    assert!(
                        (v - 0.20).abs() < 1e-12,
                        "oceanic cell ({i},{j}) S̃ = {v}, expected 0.20 exact",
                    );
                }
            }
        }
    }

    /// Acceptance #2 — with the ridge amplitude nullified
    /// (peak == base ⇒ ridge_amplitude = 0), composite collapses to
    /// pure RadialProfile output.
    #[test]
    fn composite_dome_visible_without_ridge() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);
        let mut ridge = default_orogenic_ridge();
        ridge.peak_value = 0.85;
        ridge.base_continental_value = 0.85;
        let s_composite = build(
            nx, ny, &p,
            default_radial(), ridge,
            0.20,
            // Use Fixed cap = 1.0 so a hypothetical residual ridge
            // (none expected here) wouldn't be silently clamped to a
            // value tied to the now-zeroed peak.
            CompositeCap::Fixed { value: 1.0 },
        );
        let s_radial = radial_profile::build(
            nx, ny, &p,
            0.95, 0.20,
            ProfileShape::Smoothstep,
        );
        for (a, b) in s_composite.data().iter().zip(s_radial.data().iter()) {
            assert!(
                (a - b).abs() < 1e-12,
                "composite with ridge_amplitude=0 diverged from \
                 RadialProfile: {a} vs {b}",
            );
        }
    }

    /// Acceptance #3 — with the dome collapsed (radial continental_value
    /// = oceanic_value ⇒ dome flat at 0.20), composite is essentially
    /// `0.20 + (peak-base) · ridge_amount` on continental cells. The
    /// max value across the field must reach `0.20 + (peak-base) ≈ 0.55`
    /// at the ridge axis (without ever needing the cap).
    #[test]
    fn composite_ridge_visible_without_dome() {
        let nx = 128;
        let ny = 128;
        let plates = build_plates(nx, ny, 42, 4, 0.6);
        let p = make_init_data(&plates);
        let flat_radial = CompositeRadialParams {
            continental_value: 0.20,
            profile_shape: ProfileShape::Smoothstep,
        };
        let s = build(
            nx, ny, &p,
            flat_radial, default_orogenic_ridge(),
            0.20,
            // Cap high so the test reads the raw additive value.
            CompositeCap::Fixed { value: 5.0 },
        );
        // peak - base = 1.20 - 0.85 = 0.35 ⇒ uncapped max ≈ 0.20 + 0.35
        // = 0.55. Allow a small tolerance for cell discrete sampling.
        let max_s = s.data().iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(
            max_s >= 0.45 && max_s <= 0.56,
            "max(S̃) = {max_s}, expected ≈ 0.55 (= 0.20 + 0.35 ridge amplitude)",
        );
        // All cells must lie in [oceanic, dome + ridge_amplitude],
        // i.e. [0.20, 0.55] under this degenerate dome.
        for &v in s.data() {
            assert!(
                v >= 0.20 - 1e-12 && v <= 0.55 + 1e-9,
                "cell S̃ = {v} out of expected range [0.20, 0.55]",
            );
        }
    }

    /// Acceptance #4 — the cap is strictly enforced. Use a low cap
    /// (1.10) below the uncapped composite max (1.30) and assert no
    /// cell exceeds it.
    #[test]
    fn composite_cap_respected() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);
        let s = build(
            nx, ny, &p,
            default_radial(), default_orogenic_ridge(),
            0.20,
            CompositeCap::Fixed { value: 1.10 },
        );
        let max_s = s.data().iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(
            max_s <= 1.10 + 1e-12,
            "max(S̃) = {max_s} > cap 1.10",
        );
        // At least one cell must hit (or come very close to) the cap
        // — otherwise the cap is moot for these defaults and the test
        // is not exercising the clip path.
        assert!(
            max_s >= 1.10 - 1e-9,
            "max(S̃) = {max_s} far below the cap — clip path not exercised",
        );
    }

    /// Acceptance #5 — determinism. Same seed → byte-equal output.
    #[test]
    fn composite_determinism() {
        let nx = 32;
        let ny = 32;
        let plates_a = build_plates(nx, ny, 42, 8, 0.4);
        let plates_b = build_plates(nx, ny, 42, 8, 0.4);
        let pa = make_init_data(&plates_a);
        let pb = make_init_data(&plates_b);
        let sa = build(
            nx, ny, &pa, default_radial(), default_orogenic_ridge(),
            0.20, CompositeCap::UsePeakOrogenic,
        );
        let sb = build(
            nx, ny, &pb, default_radial(), default_orogenic_ridge(),
            0.20, CompositeCap::UsePeakOrogenic,
        );
        assert_eq!(sa.data(), sb.data());
    }

    /// Acceptance #6 — somewhere on the field, S̃ approaches the cap
    /// at a plate centroid × ridge axis intersection. Same
    /// resolution-sensitivity caveat as orogenic_peak_at_centroid:
    /// use a 128² × 4-plates × 0.6-cont_ratio grid so the largest
    /// continental plate has L_plate ≥ 30 cells and σ ≥ 3 cells, well
    /// above the half-cell rounding floor.
    #[test]
    fn composite_peak_at_centroid_on_axis() {
        let nx = 128;
        let ny = 128;
        let plates = build_plates(nx, ny, 42, 4, 0.6);
        let p = make_init_data(&plates);
        let s = build(
            nx, ny, &p,
            default_radial(), default_orogenic_ridge(),
            0.20,
            CompositeCap::UsePeakOrogenic, // cap = 1.20
        );
        let max_s = s.data().iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        // Cap = 1.20. With L_plate ≥ 30 the axis-centroid cell sits
        // within 0.5 cell of the continuous centroid, giving
        // trans ≥ exp(-(0.5/(0.1·30))²) ≈ 0.97 and long_mod ≈ 1, so
        // the uncapped value at that cell is ≈ 0.95 + 0.35·0.97 = 1.29
        // → capped to 1.20. Within 1 % of the cap.
        assert!(
            max_s >= 1.19,
            "composite max(S̃) = {max_s}, expected ≥ 1.19 (cap = 1.20, large-plate setup)",
        );
        assert!(
            max_s <= 1.20 + 1e-12,
            "composite max(S̃) = {max_s} > cap 1.20",
        );
    }
}
