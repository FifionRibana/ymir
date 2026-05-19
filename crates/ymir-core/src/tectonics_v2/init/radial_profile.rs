//! Step 13 Phase 2 — radial-profile S̃ initialisation.
//!
//! Produces a continuous S̃ field that decreases smoothly from each
//! continental plate's interior (peak `continental_value`) to its
//! inter-plate boundary (where it meets `oceanic_value`). Oceanic
//! plates receive a flat `oceanic_value`.
//!
//! Motivation: Step 11 visual exploration revealed that all pre-
//! Step-13 init modes produce continents with rigid polygonal
//! Voronoi outlines and quasi-uniform interior thickness. The
//! `RadialProfile` mode addresses the first limitation by reusing
//! the Phase 1 BFS distance-to-inter-plate-boundary utility to ramp
//! the value smoothly from boundary to interior. The second
//! limitation (interior heterogeneity) is addressed in Phase 3
//! (`RadialProfileWithFBM`).
//!
//! ## Algorithm (issue D1)
//!
//! 1. Compute `dist_to_boundary[i, j]` via
//!    [`super::super::voronoi::compute_dist_to_inter_plate_boundary`]
//!    (Chebyshev 8-conn, periodic torus).
//! 2. For each plate `p`, compute
//!    `L_plate = max(dist_to_boundary[i, j])` over cells of `p`.
//!    This is the per-plate characteristic distance from boundary
//!    to plate interior, the same normalisation used in
//!    [`super::super::cratonic::factor`] (see Step 9 module
//!    docstring for the geometric reading).
//! 3. For each continental cell, normalise
//!    `t = dist / L_plate` ∈ `[0, 1]`. Apply `profile(t)` with the
//!    chosen [`ProfileShape`]. Set
//!    `S̃ = oceanic_value + (continental_value - oceanic_value) ·
//!    profile(t)`.
//! 4. Oceanic cells: `S̃ = oceanic_value` (uniform). The radial
//!    profile only applies to continental cells, since oceanic
//!    plates are physically thinner overall (issue D1).
//!
//! ## Profile shapes
//!
//! - [`ProfileShape::Smoothstep`] — cubic `3t² − 2t³`. C¹-continuous,
//!   slope at midpoint = 1.5. Default; gives a natural-looking
//!   margin.
//! - [`ProfileShape::Linear`] — slope 1 everywhere. Sharper at
//!   boundary and interior ends than Smoothstep.
//! - [`ProfileShape::Pow`] — `t^exponent`. Default exponent 1.0
//!   (= Linear). Exponent > 1 keeps most of the plate close to
//!   `oceanic_value` with a steeper rise near the interior; exponent
//!   < 1 keeps most of the plate close to `continental_value` with a
//!   steeper drop near the boundary. The UI clamps to `[0.3, 3.0]`
//!   in Phase 5; the algorithm here accepts any positive exponent.
//!
//! ## Degenerate cases
//!
//! - **Single continental plate spans the entire torus** (no
//!   inter-plate boundary): `bfs.distance` is `INFINITY` everywhere
//!   on that plate; the per-plate max-distance pass leaves
//!   `L_plate = 0`. Cells with non-finite distance saturate to
//!   `t = 1.0` (interior), so `S̃ = continental_value` uniformly —
//!   the same fall-through convention as `cratonic::factor`.
//! - **Plate consists entirely of boundary cells** (every cell at
//!   `d = 0`): `L_plate = 0` → guarded as `t = 0` →
//!   `S̃ = oceanic_value`. No division-by-zero.

use serde::{Deserialize, Serialize};

use super::super::boundaries::PlateType;
use super::super::field::Field2D;
use super::super::voronoi::compute_dist_to_inter_plate_boundary;
use super::PlateInitData;

/// Default `continental_value` for `InitMode::RadialProfile`
/// (= 0.95 ≈ 33 km, near the reference continental thickness).
pub const CONTINENTAL_VALUE_DEFAULT: f64 = 0.95;

/// Default `oceanic_value` for `InitMode::RadialProfile`
/// (= 0.20 ≈ 7 km, the reference oceanic thickness).
pub const OCEANIC_VALUE_DEFAULT: f64 = 0.20;

/// Default `Pow` exponent — equivalent to Linear, neutral starting
/// point so users see the smoothstep/linear/pow distinction by
/// adjusting only this slider in the UI.
pub const POW_EXPONENT_DEFAULT: f64 = 1.0;

/// Radial-profile shape applied to the normalised distance
/// `t ∈ [0, 1]`. See module docstring for visual interpretation.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum ProfileShape {
    /// Cubic smoothstep `3t² − 2t³`. Default.
    Smoothstep,
    /// Linear ramp `t`.
    Linear,
    /// Power profile `t^exponent`. UI clamps to `[0.3, 3.0]`;
    /// the algorithm itself accepts any positive exponent.
    Pow { exponent: f64 },
}

impl Default for ProfileShape {
    fn default() -> Self {
        ProfileShape::Smoothstep
    }
}

/// Build the S̃ field for `InitMode::RadialProfile`. See module
/// docstring for the algorithm.
///
/// # Panics
///
/// - `continental_value` or `oceanic_value` not finite.
///
/// `Pow` with non-positive exponent does not panic but produces
/// values that may exceed `[oceanic_value, continental_value]`
/// (negative exponent: blow-up at `t = 0`; zero exponent:
/// `profile = 1` everywhere, S̃ = continental_value uniformly).
/// Callers wanting a strict range should clamp the exponent to a
/// positive interval (UI does this in Phase 5).
pub fn build(
    nx: usize,
    ny: usize,
    p: &PlateInitData<'_>,
    continental_value: f64,
    oceanic_value: f64,
    profile_shape: ProfileShape,
) -> Field2D {
    assert!(
        continental_value.is_finite() && oceanic_value.is_finite(),
        "RadialProfile requires finite continental_value ({continental_value}) \
         and oceanic_value ({oceanic_value})"
    );

    let bfs = compute_dist_to_inter_plate_boundary(nx, ny, p.plate_id);

    // Per-plate L_plate = max BFS distance over cells of that plate.
    // Only used for continental cells, but we compute over all
    // plates for simplicity (oceanic entries are unused).
    // INFINITY entries (single-plate-on-torus degenerate case) are
    // skipped so `plate_max_dist[pid]` stays at 0 and the per-cell
    // branch picks up the saturation case explicitly.
    let mut plate_max_dist: Vec<f64> = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let pid = p.plate_id.get(i, j) as usize;
            if pid >= plate_max_dist.len() {
                plate_max_dist.resize(pid + 1, 0.0);
            }
            let d = bfs.distance.get(i, j);
            if d.is_finite() && d > plate_max_dist[pid] {
                plate_max_dist[pid] = d;
            }
        }
    }

    let mut s = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            if matches!(p.plate_type.get(i, j), PlateType::Oceanic) {
                s.set(i, j, oceanic_value);
                continue;
            }

            // Continental cell.
            let pid = p.plate_id.get(i, j) as usize;
            let d = bfs.distance.get(i, j);
            let l_plate = plate_max_dist[pid];
            let t = if !d.is_finite() {
                // Single-plate-on-torus: BFS never seeded inside this
                // plate. Saturate to interior (the same
                // fall-through used by cratonic::factor).
                1.0
            } else if l_plate <= 0.0 {
                // Plate consists entirely of boundary cells (all at
                // d=0). t=0 → S̃ = oceanic_value. Guards against the
                // 0/0 = NaN case.
                0.0
            } else {
                (d / l_plate).clamp(0.0, 1.0)
            };
            let p_t = apply_profile(t, profile_shape);
            let value = oceanic_value + (continental_value - oceanic_value) * p_t;
            s.set(i, j, value);
        }
    }
    s
}

#[inline]
fn apply_profile(t: f64, shape: ProfileShape) -> f64 {
    let t = t.clamp(0.0, 1.0);
    match shape {
        ProfileShape::Smoothstep => t * t * (3.0 - 2.0 * t),
        ProfileShape::Linear => t,
        ProfileShape::Pow { exponent } => t.powf(exponent),
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

    /// Acceptance #2: cells at the centroid (max BFS distance within
    /// the plate) of a continental plate hold `S̃ = continental_value`.
    #[test]
    fn continental_at_center() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);
        let s = build(nx, ny, &p, 0.95, 0.20, ProfileShape::Smoothstep);

        let bfs = compute_dist_to_inter_plate_boundary(nx, ny, p.plate_id);
        let mut found_one = false;
        for pid in 0..plates.num_plates {
            if !matches!(plates.per_plate_type[pid], PlateType::Continental) {
                continue;
            }
            let mut max_d = 0.0_f64;
            let mut max_at: Option<(usize, usize)> = None;
            for j in 0..ny {
                for i in 0..nx {
                    if plates.plate_id.get(i, j) as usize == pid {
                        let d = bfs.distance.get(i, j);
                        if d.is_finite() && d > max_d {
                            max_d = d;
                            max_at = Some((i, j));
                        }
                    }
                }
            }
            let Some((i, j)) = max_at else { continue };
            if max_d <= 0.0 {
                continue;
            }
            // At the deepest cell, t = 1, profile(1) = 1 →
            // S̃ = continental_value exactly.
            assert!(
                (s.get(i, j) - 0.95).abs() < 1e-12,
                "continental plate {} centre at ({},{}) max_d={} expected 0.95, got {}",
                pid,
                i,
                j,
                max_d,
                s.get(i, j)
            );
            found_one = true;
        }
        assert!(
            found_one,
            "no continental plate with non-degenerate L_plate found in 64² × 8 plates @ 40%"
        );
    }

    /// Acceptance #3: continental cells exactly on an inter-plate
    /// boundary (`d = 0`) hold `S̃ = oceanic_value` (since
    /// `profile(0) = 0`).
    #[test]
    fn oceanic_at_boundary() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);
        let s = build(nx, ny, &p, 0.95, 0.20, ProfileShape::Smoothstep);

        let bfs = compute_dist_to_inter_plate_boundary(nx, ny, p.plate_id);
        let mut found_continental_boundary = false;
        for j in 0..ny {
            for i in 0..nx {
                if !matches!(plates.plate_type.get(i, j), PlateType::Continental) {
                    continue;
                }
                if bfs.distance.get(i, j) == 0.0 {
                    assert!(
                        (s.get(i, j) - 0.20).abs() < 1e-12,
                        "continental boundary cell ({},{}) expected 0.20, got {}",
                        i,
                        j,
                        s.get(i, j)
                    );
                    found_continental_boundary = true;
                }
            }
        }
        assert!(
            found_continental_boundary,
            "no continental boundary cell found in 64² × 8 plates @ 40%"
        );
    }

    /// Acceptance #1 / smoothness: cell-to-cell deltas are bounded
    /// by the full amplitude `(continental - oceanic)` (a generous
    /// trivial bound; the smoothstep slope is much tighter at
    /// `1.5 · amp / L_plate` per cell). All values lie in
    /// `[oceanic_value, continental_value]` (no overshoot or NaN).
    #[test]
    fn smoothness() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);
        let s = build(nx, ny, &p, 0.95, 0.20, ProfileShape::Smoothstep);
        let amp = 0.95 - 0.20;

        for j in 0..ny {
            for i in 0..nx {
                let here = s.get(i, j);
                assert!(
                    here.is_finite()
                        && (0.20 - 1e-12..=0.95 + 1e-12).contains(&here),
                    "cell ({},{}) value {} out of [oceanic=0.20, continental=0.95]",
                    i,
                    j,
                    here
                );
                let ip = (i + 1) % nx;
                let jp = (j + 1) % ny;
                let dr = (here - s.get(ip, j)).abs();
                let dd = (here - s.get(i, jp)).abs();
                assert!(
                    dr <= amp + 1e-12,
                    "horizontal Δ at ({},{}) = {} > amplitude {}",
                    i,
                    j,
                    dr,
                    amp
                );
                assert!(
                    dd <= amp + 1e-12,
                    "vertical Δ at ({},{}) = {} > amplitude {}",
                    i,
                    j,
                    dd,
                    amp
                );
            }
        }
    }

    /// Acceptance #8: `Pow { exponent: 2.0 }` keeps more continental
    /// cells close to `oceanic_value` (the interior peak rises
    /// later, transition steeper near interior end) than
    /// `Pow { exponent: 0.5 }` (interior peak rises early,
    /// transition steep near boundary). Measured by cell count
    /// below the midpoint value.
    #[test]
    fn pow_steeper() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.4);
        let p = make_init_data(&plates);

        let s_steep = build(nx, ny, &p, 0.95, 0.20, ProfileShape::Pow { exponent: 2.0 });
        let s_gentle = build(nx, ny, &p, 0.95, 0.20, ProfileShape::Pow { exponent: 0.5 });

        let mid = 0.5 * (0.95 + 0.20);
        let mut count_steep_below = 0;
        let mut count_gentle_below = 0;
        for j in 0..ny {
            for i in 0..nx {
                if !matches!(plates.plate_type.get(i, j), PlateType::Continental) {
                    continue;
                }
                if s_steep.get(i, j) < mid {
                    count_steep_below += 1;
                }
                if s_gentle.get(i, j) < mid {
                    count_gentle_below += 1;
                }
            }
        }
        assert!(
            count_steep_below > count_gentle_below,
            "Pow{{2}} ({} cells below midpoint) is not steeper than Pow{{0.5}} ({} cells)",
            count_steep_below,
            count_gentle_below
        );
    }

    /// D1: every cell in an oceanic plate has `S̃ = oceanic_value`
    /// exactly (no radial profile applied).
    #[test]
    fn oceanic_plates_uniform() {
        let nx = 64;
        let ny = 64;
        let plates = build_plates(nx, ny, 42, 8, 0.3);
        let p = make_init_data(&plates);
        let s = build(nx, ny, &p, 0.95, 0.20, ProfileShape::Smoothstep);

        let mut count_oceanic = 0;
        for j in 0..ny {
            for i in 0..nx {
                if matches!(plates.plate_type.get(i, j), PlateType::Oceanic) {
                    assert_eq!(
                        s.get(i, j),
                        0.20,
                        "oceanic cell ({},{}) expected 0.20, got {}",
                        i,
                        j,
                        s.get(i, j)
                    );
                    count_oceanic += 1;
                }
            }
        }
        assert!(
            count_oceanic > 0,
            "no oceanic cells found in 64² × 8 plates @ 30% continental"
        );
    }

    /// Linear and Pow{1.0} produce identical fields (Pow with
    /// exponent 1.0 is exactly the identity).
    #[test]
    fn pow_one_equals_linear() {
        let nx = 32;
        let ny = 32;
        let plates = build_plates(nx, ny, 42, 6, 0.4);
        let p = make_init_data(&plates);
        let s_linear = build(nx, ny, &p, 0.95, 0.20, ProfileShape::Linear);
        let s_pow1 = build(nx, ny, &p, 0.95, 0.20, ProfileShape::Pow { exponent: 1.0 });
        for j in 0..ny {
            for i in 0..nx {
                let l = s_linear.get(i, j);
                let q = s_pow1.get(i, j);
                assert!(
                    (l - q).abs() < 1e-15,
                    "Linear vs Pow{{1.0}} differ at ({},{}): {} vs {}",
                    i,
                    j,
                    l,
                    q
                );
            }
        }
    }

    /// Determinism: same `(plate_id, plate_type, params)` →
    /// byte-identical output.
    #[test]
    fn deterministic_same_inputs() {
        let nx = 32;
        let ny = 32;
        let plates_a = build_plates(nx, ny, 42, 6, 0.4);
        let plates_b = build_plates(nx, ny, 42, 6, 0.4);
        let s_a = build(nx, ny, &make_init_data(&plates_a), 0.95, 0.20, ProfileShape::Smoothstep);
        let s_b = build(nx, ny, &make_init_data(&plates_b), 0.95, 0.20, ProfileShape::Smoothstep);
        assert_eq!(s_a.data(), s_b.data());
    }
}
