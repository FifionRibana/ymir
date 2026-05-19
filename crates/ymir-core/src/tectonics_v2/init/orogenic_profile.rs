//! Step 12 R7.A.1.2 — orogenic continental S̃ initialisation.
//!
//! Builds a continuous S̃ field where each continental plate carries a
//! single **linear ridge** (orogenic chain) along its principal axis,
//! instead of the radial peak pattern of [`super::radial_profile`].
//!
//! ## Algorithm (R7.A.1.1 §3)
//!
//! 1. For each continental plate `p` compute:
//!    - centroid `(cx, cy)` via **periodic-aware circular mean** of the
//!      cells' coordinates (the cells of a plate may wrap around the
//!      torus; a naive Cartesian average would place the centroid on
//!      the wrong side of the domain)
//!    - principal axis `(ux, uy)` via PCA on the **unwrapped**
//!      coordinates relative to the centroid (minimum-image
//!      convention). Returns `None` (→ fallback) for plates with
//!      fewer than 5 cells or with rank-1 covariance (cells colinear
//!      / point-like).
//!    - `L_plate` = the max BFS distance to the inter-plate boundary
//!      over cells of `p` (reused from [`super::radial_profile`]).
//! 2. For each continental cell `(x, y)`:
//!    - unwrap `(x − cx, y − cy)` to minimum-image vector `(dx, dy)`
//!    - project: `d_along = dx·ux + dy·uy`,
//!      `d_perp  = dx·(−uy) + dy·ux`
//!    - longitudinal modulation:
//!      `t_along = (1 − |d_along| / half_length).clamp(0, 1)`
//!      `long_mod = 3·t_along² − 2·t_along³` (smoothstep)
//!    - transversal Gaussian: `trans = exp(−(d_perp / width_sigma)²)`
//!    - `ridge_amount = long_mod · trans` in `[0, 1]`
//!    - `S̃ = base + (peak − base) · ridge_amount`
//! 3. Oceanic cells: `S̃ = oceanic_value` uniform (independent of
//!    plate orientation; the orogenic mechanism only acts on
//!    continental crust).
//!
//! ## Scale conventions
//!
//! `half_length = half_length_ratio · L_plate` and
//! `width_sigma = width_sigma_ratio · L_plate`, both in cell units
//! (the same unit as the BFS distance). With the default ratios
//! `(0.40, 0.08)` the ridge spans roughly 80% of the plate diameter
//! longitudinally and has a transversal sigma of ~8% of L_plate.
//!
//! ## Degenerate cases
//!
//! - **Plate with < 5 cells or rank-1 covariance** → fall back to the
//!   user-supplied [`OrogenicOrientation::Fixed`] direction (default
//!   `angle_rad = 0.0`, i.e. axis +x). PCA failure is silent at
//!   runtime; tests on synthetic plates cover the fallback branch.
//! - **`L_plate = 0`** (plate is entirely on the inter-plate boundary)
//!   → `S̃ = base_continental_value` uniform on that plate (ridge
//!   degenerates to zero amplitude). Same convention as
//!   `radial_profile` for this corner case.
//! - **Single-plate-on-torus** (BFS never reaches an inter-plate
//!   boundary, all distances are `INFINITY`): same as `L_plate = 0`
//!   above — the plate has no characteristic boundary distance, so
//!   the ridge has no scale; we degenerate to the base value.
//!
//! ## Side notes
//!
//! - The output is a **scalar field**, so no div-free constraint
//!   applies (unlike `MantlePattern`, which is built as `curl(ψ)` and
//!   needs the staggered grid trick to preserve `div = 0`).
//! - Determinism: the only stochastic element is the Voronoï seed
//!   (upstream of this module). Same plates ⇒ byte-equal output.

use serde::{Deserialize, Serialize};

use super::super::boundaries::PlateType;
use super::super::field::Field2D;
use super::super::voronoi::compute_dist_to_inter_plate_boundary;
use super::PlateInitData;

/// Default `peak_value` — S̃ at the ridge axis. Continental highland
/// (Himalaya / Andes proxy) sits roughly 5–8 km above the reference
/// continental thickness `~ 1.0`.
pub const OROGENIC_PEAK_VALUE_DEFAULT: f64 = 1.20;

/// Default `base_continental_value` — S̃ far from the ridge but still
/// inside the continent. Continental lowland (`~ 28` km crust)
/// corresponds to ~0.85 in dimensionless thickness.
pub const OROGENIC_BASE_VALUE_DEFAULT: f64 = 0.85;

/// Default `oceanic_value` — same as [`super::OCEANIC_S_DEFAULT`] /
/// [`super::radial_profile::OCEANIC_VALUE_DEFAULT`].
pub const OROGENIC_OCEANIC_VALUE_DEFAULT: f64 = 0.20;

/// Default `half_length_ratio` — the ridge spans `2 · 0.40 = 80 %` of
/// the plate's characteristic distance `L_plate`.
pub const OROGENIC_HALF_LENGTH_RATIO_DEFAULT: f64 = 0.40;

/// Default `width_sigma_ratio` — `σ = 0.08 · L_plate`. At 32² this
/// resolves to ~1–2 cells (thin ridge, may need a bigger σ — see
/// R7.A.1.3 visual check). At 64² it resolves to ~3–5 cells.
pub const OROGENIC_WIDTH_SIGMA_RATIO_DEFAULT: f64 = 0.08;

/// Minimum number of cells a plate must hold for PCA-based orientation
/// to be considered reliable. Below this, the orogenic build falls
/// back to the configured fixed orientation.
const ORO_MIN_CELLS_FOR_PCA: usize = 5;

/// Minimum determinant of the covariance matrix for PCA to be
/// considered well-conditioned (avoids degenerate rank-1 inputs like a
/// 1-cell-wide line). Cell coordinates are in `[0, nx) × [0, ny)`, so
/// this scale-independent threshold matches the "two non-degenerate
/// axes of extent > sqrt(threshold) cells" intuition.
const ORO_PCA_MIN_DET: f64 = 0.5;

/// Orientation strategy for the orogenic ridge.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum OrogenicOrientation {
    /// Compute the principal axis of the plate's cells via PCA, with
    /// periodic-torus handling. Falls back to
    /// `Fixed { angle_rad: 0.0 }` for plates with fewer than 5 cells
    /// or with degenerate (rank-1) covariance.
    PlateMainAxisPca,
    /// Use a constant orientation for every plate. Useful for unit
    /// tests, presets that want a global "trend" orientation, or as
    /// a fallback path the user can prove is hit on a given seed.
    Fixed { angle_rad: f64 },
}

impl Default for OrogenicOrientation {
    fn default() -> Self {
        OrogenicOrientation::PlateMainAxisPca
    }
}

/// Build the S̃ field for `InitMode::Orogenic`. See module docstring
/// for the full algorithm; the parameter rationale and design choices
/// are documented in `docs/reports/step12_r7_a_orogenic_profile/`.
///
/// # Panics
///
/// - Any of `peak_value`, `base_continental_value`, `oceanic_value`,
///   `half_length_ratio`, `width_sigma_ratio` is non-finite.
/// - `half_length_ratio <= 0.0` or `width_sigma_ratio <= 0.0` (the
///   ridge would degenerate to a point/line of zero amplitude — this
///   is a user-config error, not a runtime case).
pub fn build(
    nx: usize,
    ny: usize,
    p: &PlateInitData<'_>,
    peak_value: f64,
    base_continental_value: f64,
    oceanic_value: f64,
    half_length_ratio: f64,
    width_sigma_ratio: f64,
    orientation: OrogenicOrientation,
) -> Field2D {
    assert!(
        peak_value.is_finite()
            && base_continental_value.is_finite()
            && oceanic_value.is_finite(),
        "orogenic build: peak/base/oceanic must be finite (got {peak_value}, \
         {base_continental_value}, {oceanic_value})",
    );
    assert!(
        half_length_ratio > 0.0 && half_length_ratio.is_finite(),
        "orogenic build: half_length_ratio must be > 0 (got {half_length_ratio})",
    );
    assert!(
        width_sigma_ratio > 0.0 && width_sigma_ratio.is_finite(),
        "orogenic build: width_sigma_ratio must be > 0 (got {width_sigma_ratio})",
    );

    // ---- Pass 1 — per-plate accumulators (circular sums for the
    // periodic-aware centroid, cell count, plate type tag). ----
    let mut plate_meta: Vec<PlateAccum> = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let pid = p.plate_id.get(i, j) as usize;
            if pid >= plate_meta.len() {
                plate_meta.resize_with(pid + 1, PlateAccum::new);
            }
            plate_meta[pid].accumulate_pos(i, j, nx, ny);
            plate_meta[pid].plate_type = Some(p.plate_type.get(i, j));
        }
    }

    // ---- Per-plate geometry — centroid, PCA axis, L_plate. ----
    let bfs = compute_dist_to_inter_plate_boundary(nx, ny, p.plate_id);
    for acc in plate_meta.iter_mut() {
        acc.finalise_centroid(nx, ny);
    }
    // Second pass to compute covariance now that centroids are known.
    for j in 0..ny {
        for i in 0..nx {
            let pid = p.plate_id.get(i, j) as usize;
            plate_meta[pid].accumulate_covariance(i, j, nx, ny);
        }
    }
    // Per-plate L_plate from BFS (only meaningful for continental
    // plates but computed for all uniformly; oceanic entries unused).
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
    // Resolve PCA axes (or fallback).
    let mut plate_axis: Vec<(f64, f64)> = vec![(1.0, 0.0); plate_meta.len()];
    let fallback_axis = match orientation {
        OrogenicOrientation::PlateMainAxisPca => (1.0, 0.0),
        OrogenicOrientation::Fixed { angle_rad } => {
            (angle_rad.cos(), angle_rad.sin())
        }
    };
    for (pid, acc) in plate_meta.iter().enumerate() {
        plate_axis[pid] = match orientation {
            OrogenicOrientation::Fixed { angle_rad } => {
                (angle_rad.cos(), angle_rad.sin())
            }
            OrogenicOrientation::PlateMainAxisPca => {
                acc.principal_axis().unwrap_or(fallback_axis)
            }
        };
    }

    // ---- Pass 3 — build the field. ----
    let mut s = Field2D::new(nx, ny);
    let nxf = nx as f64;
    let nyf = ny as f64;
    for j in 0..ny {
        for i in 0..nx {
            // Oceanic cells: flat value, no ridge.
            if matches!(p.plate_type.get(i, j), PlateType::Oceanic) {
                s.set(i, j, oceanic_value);
                continue;
            }
            let pid = p.plate_id.get(i, j) as usize;
            let acc = &plate_meta[pid];
            let l_plate = plate_l[pid];
            if l_plate <= 0.0 {
                // Single-plate-on-torus or fully-on-boundary plate
                // (see module docstring). Ridge has no scale; fall
                // back to the base continental value.
                s.set(i, j, base_continental_value);
                continue;
            }
            let half_length = (half_length_ratio * l_plate).max(1e-12);
            let width_sigma = (width_sigma_ratio * l_plate).max(1e-12);
            let cx = acc.centroid_x.expect("continental cell ⇒ centroid set");
            let cy = acc.centroid_y.expect("continental cell ⇒ centroid set");
            let cell_x = i as f64 + 0.5;
            let cell_y = j as f64 + 0.5;
            let dx = min_image_delta(cell_x - cx, nxf);
            let dy = min_image_delta(cell_y - cy, nyf);
            let (ux, uy) = plate_axis[pid];
            let d_along = dx * ux + dy * uy;
            let d_perp = -dx * uy + dy * ux;
            let t_along = (1.0 - d_along.abs() / half_length).clamp(0.0, 1.0);
            let long_mod = t_along * t_along * (3.0 - 2.0 * t_along);
            let dp_over_sigma = d_perp / width_sigma;
            let trans = (-(dp_over_sigma * dp_over_sigma)).exp();
            let ridge_amount = long_mod * trans;
            let value = base_continental_value
                + (peak_value - base_continental_value) * ridge_amount;
            s.set(i, j, value);
        }
    }
    s
}

/// Minimum-image delta on a torus of length `n`. Returns the signed
/// difference in `[-n/2, n/2)` that is closest to zero (modulo `n`).
///
/// Exposed at `pub(crate)` from R7.A.2 so [`super::composite_profile`]
/// can share the same wrap convention as the orogenic ridge.
#[inline]
pub(crate) fn min_image_delta(raw: f64, n: f64) -> f64 {
    let half = 0.5 * n;
    let mut d = raw;
    while d > half {
        d -= n;
    }
    while d < -half {
        d += n;
    }
    d
}

/// Per-plate scratch accumulator used by the two preparatory passes
/// over the grid. After [`Self::accumulate_pos`] runs over every
/// cell, [`Self::finalise_centroid`] derives the periodic-aware
/// centroid via the circular-mean trick. After
/// [`Self::accumulate_covariance`] runs (which needs the centroid
/// already), [`Self::principal_axis`] returns the unit eigenvector
/// of the larger eigenvalue of the 2×2 covariance, or `None` for
/// degenerate inputs.
///
/// Exposed at `pub(crate)` from R7.A.2 so
/// [`super::composite_profile`] can reuse the exact same plate
/// geometry pipeline as the standalone Orogenic mode — same
/// centroid, same PCA axis, same fallback semantics. Visibility is
/// the only change; the algorithm is unchanged so
/// `InitMode::Orogenic` output stays byte-equal to R7.A.1.
pub(crate) struct PlateAccum {
    n: usize,
    // Circular-mean accumulators for x and y (unit-circle projections).
    cos_x_sum: f64,
    sin_x_sum: f64,
    cos_y_sum: f64,
    sin_y_sum: f64,
    // Centroid resolved by finalise_centroid (cell-coord units, in
    // [0, nx) / [0, ny)). None if the plate had no cells (skipped
    // pid) — should not happen in practice since accumulate_pos is
    // called exactly once per cell.
    centroid_x: Option<f64>,
    centroid_y: Option<f64>,
    // Covariance accumulators on unwrapped coords relative to the
    // centroid. Filled by accumulate_covariance.
    sum_dx2: f64,
    sum_dy2: f64,
    sum_dxdy: f64,
    plate_type: Option<PlateType>,
}

impl PlateAccum {
    pub(crate) fn new() -> Self {
        Self {
            n: 0,
            cos_x_sum: 0.0,
            sin_x_sum: 0.0,
            cos_y_sum: 0.0,
            sin_y_sum: 0.0,
            centroid_x: None,
            centroid_y: None,
            sum_dx2: 0.0,
            sum_dy2: 0.0,
            sum_dxdy: 0.0,
            plate_type: None,
        }
    }

    /// Centroid `(x, y)` in cell coordinates after
    /// [`Self::finalise_centroid`] has been called. `None` while the
    /// accumulator is still in the "positions-only" phase.
    pub(crate) fn centroid(&self) -> Option<(f64, f64)> {
        match (self.centroid_x, self.centroid_y) {
            (Some(x), Some(y)) => Some((x, y)),
            _ => None,
        }
    }

    pub(crate) fn accumulate_pos(&mut self, i: usize, j: usize, nx: usize, ny: usize) {
        use std::f64::consts::TAU;
        self.n += 1;
        // Project each cell-center coord onto the unit circle of its
        // axis. Circular mean of these projections gives a centroid
        // angle ∈ [0, 2π) which we map back to the [0, n) cell range.
        let theta_x = TAU * (i as f64 + 0.5) / nx as f64;
        let theta_y = TAU * (j as f64 + 0.5) / ny as f64;
        self.cos_x_sum += theta_x.cos();
        self.sin_x_sum += theta_x.sin();
        self.cos_y_sum += theta_y.cos();
        self.sin_y_sum += theta_y.sin();
    }

    pub(crate) fn finalise_centroid(&mut self, nx: usize, ny: usize) {
        use std::f64::consts::TAU;
        if self.n == 0 {
            return;
        }
        let inv_tau = 1.0 / TAU;
        // atan2 returns [-π, π]; +TAU if negative to keep [0, TAU).
        let angle_x = self.sin_x_sum.atan2(self.cos_x_sum);
        let angle_y = self.sin_y_sum.atan2(self.cos_y_sum);
        let ax = if angle_x < 0.0 { angle_x + TAU } else { angle_x };
        let ay = if angle_y < 0.0 { angle_y + TAU } else { angle_y };
        // Map back from [0, TAU) → [0, n).
        let cx = ax * inv_tau * nx as f64;
        let cy = ay * inv_tau * ny as f64;
        // Snap into [0, n) defensively (floating arithmetic may
        // produce exactly `n` for cells right at the wrap).
        let cx = if cx >= nx as f64 { 0.0 } else { cx };
        let cy = if cy >= ny as f64 { 0.0 } else { cy };
        self.centroid_x = Some(cx);
        self.centroid_y = Some(cy);
    }

    pub(crate) fn accumulate_covariance(&mut self, i: usize, j: usize, nx: usize, ny: usize) {
        let (cx, cy) = match (self.centroid_x, self.centroid_y) {
            (Some(a), Some(b)) => (a, b),
            _ => return,
        };
        let dx = min_image_delta(i as f64 + 0.5 - cx, nx as f64);
        let dy = min_image_delta(j as f64 + 0.5 - cy, ny as f64);
        self.sum_dx2 += dx * dx;
        self.sum_dy2 += dy * dy;
        self.sum_dxdy += dx * dy;
    }

    /// Unit eigenvector of the larger eigenvalue of the 2×2 sample
    /// covariance matrix. Returns `None` for plates with too few
    /// cells (PCA unreliable) or with rank-1 / near-isotropic
    /// covariance (caller falls back).
    pub(crate) fn principal_axis(&self) -> Option<(f64, f64)> {
        if self.n < ORO_MIN_CELLS_FOR_PCA {
            return None;
        }
        let inv_n = 1.0 / self.n as f64;
        let cxx = self.sum_dx2 * inv_n;
        let cyy = self.sum_dy2 * inv_n;
        let cxy = self.sum_dxdy * inv_n;
        // Determinant test: if both extents are small or one is much
        // larger than the other we still want to return that one
        // axis. Use the discriminant of the characteristic
        // polynomial instead — λ₊ − λ₋ = sqrt((cxx − cyy)² + 4cxy²).
        // The "well-conditioned" check is on this spread; if
        // eigenvalues are within ORO_PCA_MIN_DET of each other in
        // absolute terms, treat as near-isotropic and fall back.
        let diff = cxx - cyy;
        let disc_sq = diff * diff + 4.0 * cxy * cxy;
        if disc_sq < ORO_PCA_MIN_DET {
            // Eigenvalues too close → axis indeterminate. Caller
            // falls back to the configured fixed orientation.
            return None;
        }
        let disc = disc_sq.sqrt();
        let lambda_plus = 0.5 * (cxx + cyy + disc);
        // Eigenvector of `λ₊` satisfies `(C - λ₊·I) v = 0`. Pick the
        // row with the larger off-diagonal magnitude to avoid the
        // (cxy = 0) degenerate path.
        let (vx, vy) = if cxy.abs() > 1e-12 {
            (lambda_plus - cyy, cxy)
        } else if cxx >= cyy {
            (1.0, 0.0)
        } else {
            (0.0, 1.0)
        };
        let norm = (vx * vx + vy * vy).sqrt();
        if norm < 1e-12 {
            return None;
        }
        Some((vx / norm, vy / norm))
    }
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

    /// Acceptance #3 — every oceanic cell holds exactly `oceanic_value`
    /// (independent of plate orientation; ridge formula is gated
    /// behind the continental classification).
    #[test]
    fn orogenic_oceanic_uniform() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);
        let s = build(
            nx, ny, &p,
            1.20, 0.85, 0.20,
            0.40, 0.08,
            OrogenicOrientation::PlateMainAxisPca,
        );
        for j in 0..ny {
            for i in 0..nx {
                if matches!(plates.plate_type.get(i, j), PlateType::Oceanic) {
                    let v = s.get(i, j);
                    assert!(
                        (v - 0.20).abs() < 1e-12,
                        "oceanic cell ({i},{j}) = {v}, expected 0.20 exact",
                    );
                }
            }
        }
    }

    /// Acceptance #2 — the continental classification is independent
    /// of the init mode. Counting "cells whose S̃ is at or above the
    /// base value" (= continental cells) on the orogenic build must
    /// match the same count on the RadialProfile build for the same
    /// Voronoï layout.
    #[test]
    fn orogenic_continental_fraction_correct() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);
        let s_oro = build(
            nx, ny, &p,
            1.20, 0.85, 0.20,
            0.40, 0.08,
            OrogenicOrientation::PlateMainAxisPca,
        );
        let count_cont_from_plate_type = (0..ny)
            .flat_map(|j| (0..nx).map(move |i| (i, j)))
            .filter(|&(i, j)| matches!(plates.plate_type.get(i, j), PlateType::Continental))
            .count();
        let count_cont_from_oro = s_oro
            .data()
            .iter()
            .filter(|&&v| v > 0.20 + 1e-9)
            .count();
        // Orogenic continental cells all have S̃ ≥ base = 0.85 > 0.20.
        // Oceanic cells have S̃ = 0.20. So `S̃ > 0.20 + eps` strictly
        // identifies the continental subset.
        assert_eq!(
            count_cont_from_plate_type, count_cont_from_oro,
            "orogenic continental cell count ({count_cont_from_oro}) diverged \
             from the Voronoï classification ({count_cont_from_plate_type})",
        );
    }

    /// Acceptance #4 — somewhere on a continental plate, S̃ reaches
    /// (or comes within numerical tolerance of) the configured
    /// `peak_value`.
    ///
    /// **Resolution-sensitivity caveat** (acknowledged in R7.A.1.1
    /// caveats): with the default `width_sigma_ratio = 0.08`, small
    /// plates (`L_plate < ~12` cells) have `σ < 1` cell so even the
    /// cell closest to the ridge axis sits at ≳ 1 σ from the axis,
    /// and the half-cell rounding offset relative to the continuous
    /// centroid cuts the peak attainment to ~75 %. The contract
    /// "max(S̃) approaches `peak_value`" therefore requires a plate
    /// big enough that σ ≥ ~2 cells.
    ///
    /// This test uses a 128² grid with 4 plates @ 0.6 continental
    /// ratio so the largest continental plate hosts `L_plate ≥ 30`
    /// cells (σ ≥ 2.4 cells); under these conditions, at least one
    /// cell along the ridge axis must reach within 5 % of `peak`.
    #[test]
    fn orogenic_peak_at_centroid() {
        let nx = 128;
        let ny = 128;
        let plates = build_plates(nx, ny, 42, 4, 0.6);
        let p = make_init_data(&plates);
        let s = build(
            nx, ny, &p,
            1.20, 0.85, 0.20,
            0.40, 0.08,
            OrogenicOrientation::PlateMainAxisPca,
        );
        // Global max over all cells — independent of which exact
        // (i, j) cell falls on the ridge axis. Robust to half-cell
        // rounding and PCA orientation.
        let max_s = s.data().iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(
            max_s >= 1.14,
            "max(S̃) = {max_s}, expected ≥ 1.14 (peak = 1.20, within 5 % \
             tolerance for half-cell offset + cell-discrete sampling)",
        );
        // Sanity: max must not exceed peak (formula is bounded by
        // construction: ridge_amount ∈ [0, 1] ⇒ S̃ ∈ [base, peak]).
        assert!(
            max_s <= 1.20 + 1e-12,
            "max(S̃) = {max_s} exceeded peak = 1.20",
        );
    }

    /// Acceptance #5 — far from the ridge (longitudinally outside
    /// `half_length`) the field decays to `base_continental_value`.
    /// Iterate over continental cells and find at least one that is
    /// at `t_along = 0` (i.e. `|d_along| ≥ half_length`); assert S̃
    /// equals base there.
    #[test]
    fn orogenic_decays_far_longitudinal() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);
        let s = build(
            nx, ny, &p,
            1.20, 0.85, 0.20,
            0.40, 0.08,
            OrogenicOrientation::PlateMainAxisPca,
        );
        // For each continental cell, value can be no lower than the
        // base when ridge_amount = 0. Conversely, ridge_amount > 0
        // requires the cell to be inside `half_length` longitudinally
        // AND within a few σ perpendicularly. So we just check: every
        // continental cell holds S̃ ≥ base − ε, S̃ ≤ peak + ε.
        let mut found_base = false;
        for j in 0..ny {
            for i in 0..nx {
                if !matches!(plates.plate_type.get(i, j), PlateType::Continental) {
                    continue;
                }
                let v = s.get(i, j);
                assert!(
                    v >= 0.85 - 1e-12 && v <= 1.20 + 1e-12,
                    "continental cell ({i},{j}) S̃ = {v} out of [base, peak] = [0.85, 1.20]",
                );
                if (v - 0.85).abs() < 1e-9 {
                    found_base = true;
                }
            }
        }
        // At least one continental cell must be far enough from any
        // ridge that ridge_amount underflows to 0 (exp(-something
        // huge) = 0 in f64 around 745). With sigma = 0.08 · L_plate
        // and 64² × 8 plates there will always be some such cell.
        assert!(
            found_base,
            "no continental cell decayed exactly to the base value 0.85; \
             ridge may be wider than the plate (check defaults)",
        );
    }

    /// Acceptance #6 — same seed → byte-equal output. Determinism is
    /// not a contract by accident: the centroid + PCA pipeline reads
    /// floats and uses `atan2`, both of which are deterministic on
    /// IEEE-754 hardware. The test pins it.
    #[test]
    fn orogenic_determinism() {
        let nx = 32;
        let ny = 32;
        let plates_a = build_plates(nx, ny, 42, 8, 0.4);
        let plates_b = build_plates(nx, ny, 42, 8, 0.4);
        let pa = make_init_data(&plates_a);
        let pb = make_init_data(&plates_b);
        let sa = build(
            nx, ny, &pa, 1.20, 0.85, 0.20, 0.40, 0.08,
            OrogenicOrientation::PlateMainAxisPca,
        );
        let sb = build(
            nx, ny, &pb, 1.20, 0.85, 0.20, 0.40, 0.08,
            OrogenicOrientation::PlateMainAxisPca,
        );
        assert_eq!(sa.data(), sb.data());
    }

    /// Periodic-aware centroid: a synthetic plate split across the
    /// torus wrap (cells at i ∈ {0, 1, nx−1}) should produce a
    /// centroid at i ≈ 0 (the cluster's geometric centre), NOT at
    /// i ≈ (nx − 1) / 3 (the naive Cartesian mean).
    #[test]
    fn periodic_centroid_handles_wrap() {
        let nx = 32;
        let ny = 32;
        let mut acc = PlateAccum::new();
        // Cluster around i = 0 with one cell at the right edge.
        acc.accumulate_pos(nx - 1, 16, nx, ny);
        acc.accumulate_pos(0, 16, nx, ny);
        acc.accumulate_pos(1, 16, nx, ny);
        acc.finalise_centroid(nx, ny);
        let cx = acc.centroid_x.unwrap();
        // The naive Cartesian mean of {nx-1=31, 0, 1} is 32/3 ≈
        // 10.67 (wrong). The periodic mean should be ≈ 0.5 (the
        // centroid of the cluster crossing the wrap). Pin within ±1.
        assert!(
            cx < 1.5 || cx > nx as f64 - 1.5,
            "periodic centroid x = {cx}, expected ≈ 0 (wrap handling broken)",
        );
    }

    /// PCA fallback path: a plate with fewer than 5 cells should
    /// produce `principal_axis() == None`, prompting the build to
    /// use the configured fixed orientation.
    #[test]
    fn pca_falls_back_below_min_cells() {
        let nx = 32;
        let ny = 32;
        let mut acc = PlateAccum::new();
        acc.accumulate_pos(10, 10, nx, ny);
        acc.accumulate_pos(11, 10, nx, ny);
        acc.accumulate_pos(12, 10, nx, ny);
        acc.finalise_centroid(nx, ny);
        for i in 10..13 {
            acc.accumulate_covariance(i, 10, nx, ny);
        }
        assert!(acc.principal_axis().is_none());
    }

    /// PCA fallback path: even with enough cells, a nearly-isotropic
    /// distribution (square cluster) should fall back since no axis
    /// dominates.
    #[test]
    fn pca_falls_back_on_isotropic() {
        let nx = 32;
        let ny = 32;
        let mut acc = PlateAccum::new();
        for j in 14..18 {
            for i in 14..18 {
                acc.accumulate_pos(i, j, nx, ny);
            }
        }
        acc.finalise_centroid(nx, ny);
        for j in 14..18 {
            for i in 14..18 {
                acc.accumulate_covariance(i, j, nx, ny);
            }
        }
        // Eigenvalues for a 4×4 square cluster about its centre:
        // var_x = var_y = mean((i − 15.5)²) = (0.5² + 1.5²) · 2 / 4
        //                = (0.25 + 2.25) / 2 = 1.25 each, cxy = 0.
        // Discriminant = (0)² + 4·0 = 0 → falls under ORO_PCA_MIN_DET.
        assert!(acc.principal_axis().is_none());
    }

    /// PCA detects an elongated cluster: synthetic plate elongated
    /// along the x axis must yield an axis close to ±(1, 0).
    #[test]
    fn pca_detects_x_elongation() {
        let nx = 32;
        let ny = 32;
        let mut acc = PlateAccum::new();
        // 10-cell-long line along x at j = 16.
        for i in 5..15 {
            acc.accumulate_pos(i, 16, nx, ny);
        }
        acc.finalise_centroid(nx, ny);
        for i in 5..15 {
            acc.accumulate_covariance(i, 16, nx, ny);
        }
        let (ux, uy) = acc.principal_axis().expect("non-degenerate axis");
        // Direction is unsigned (axis), so test |ux| ≈ 1, |uy| ≈ 0.
        assert!(
            ux.abs() > 0.95 && uy.abs() < 0.05,
            "expected x-axis, got ({ux}, {uy})",
        );
    }

    /// PCA detects an elongated cluster along y.
    #[test]
    fn pca_detects_y_elongation() {
        let nx = 32;
        let ny = 32;
        let mut acc = PlateAccum::new();
        for j in 5..15 {
            acc.accumulate_pos(16, j, nx, ny);
        }
        acc.finalise_centroid(nx, ny);
        for j in 5..15 {
            acc.accumulate_covariance(16, j, nx, ny);
        }
        let (ux, uy) = acc.principal_axis().expect("non-degenerate axis");
        assert!(
            uy.abs() > 0.95 && ux.abs() < 0.05,
            "expected y-axis, got ({ux}, {uy})",
        );
    }

    /// Fixed orientation path: regardless of plate shape, the
    /// returned axis is exactly `(cos θ, sin θ)`. Verify by passing a
    /// plate that PCA would happily resolve (10-cell x-line) but
    /// requesting `Fixed { angle_rad: π/2 }` (y axis).
    #[test]
    fn fixed_orientation_overrides_pca() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);
        let s_pca = build(
            nx, ny, &p, 1.20, 0.85, 0.20, 0.40, 0.08,
            OrogenicOrientation::PlateMainAxisPca,
        );
        let s_fix = build(
            nx, ny, &p, 1.20, 0.85, 0.20, 0.40, 0.08,
            OrogenicOrientation::Fixed { angle_rad: std::f64::consts::FRAC_PI_2 },
        );
        // The two fields must differ on at least one continental
        // cell — same plates, but different orientation strategy.
        let mut diff = false;
        for (a, b) in s_pca.data().iter().zip(s_fix.data().iter()) {
            if (a - b).abs() > 1e-9 {
                diff = true;
                break;
            }
        }
        assert!(diff, "PCA and Fixed produced identical S̃ field");
    }
}
