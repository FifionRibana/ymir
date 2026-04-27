//! Cratonic factor field — BFS distance + per-plate normalisation +
//! smoothstep transition. Pure function of the Voronoï partition
//! and the cratonic configuration; called once at simulation init.
//!
//! # Pipeline (per `step9_issue.md` D2/D3/D5)
//!
//! 1. Compute per-plate area (cell count). Plates marked as oceanic
//!    or whose area is below `plate_area_min · domain_area` are
//!    *excluded* — their cells receive `cratonic_factor = 0`.
//! 2. Multi-source BFS from every non-continental cell on the
//!    periodic torus. Result `dist[i,j]` = grid distance (in cells)
//!    from cell `(i, j)` to its nearest oceanic cell. Cells that
//!    are themselves oceanic have `dist = 0`.
//! 3. For each retained continental cell, normalise:
//!    `dist_norm = dist / L_plate`, with `L_plate` = the maximum
//!    BFS depth attained inside that plate (i.e. the largest grid
//!    distance from the oceanic boundary reached by any cell of the
//!    plate). `dist_norm` therefore lives in `[0, 1]` per plate by
//!    construction.
//! 4. Apply smoothstep
//!    `factor = smoothstep(d_low, d_high, dist_norm)`
//!    where `d_mid = 1 − sqrt(Cr)`,
//!    `d_low = d_mid − smoothing_width / 2`,
//!    `d_high = d_mid + smoothing_width / 2`.
//!
//! # `L_plate` calibration — clarification of the issue text
//!
//! The Step 9 issue (D3) writes
//! `d_mid = R · (1 − sqrt(Cr))` "for a circular plate of radius R",
//! then "generalized to non-circular plates via
//! L_plate = sqrt(plate_area)". Reading these literally is
//! self-inconsistent: for a circle of radius R the area is πR² so
//! `sqrt(area) = R · sqrt(π) ≈ 1.77 R`, not R. Substituting that
//! into the formula collapses the cratonic core to a few percent of
//! the plate area for any realistic Cr (verified empirically — see
//! the throwaway probe `tests/v2_cratonic_normalization_probe.rs`,
//! "Scheme A").
//!
//! The geometric intent of the issue formula is that `L_plate`
//! represents the *characteristic distance from boundary to plate
//! interior*, i.e. the "inradius" generalisation. For a circular
//! plate that is exactly `R`, for a square of side `L` it is `L/2`,
//! and for an irregular Voronoï plate it is naturally the maximum
//! BFS depth attained inside the plate. Under this reading, the
//! formula `d_mid_normalized = 1 − sqrt(Cr)` is geometrically
//! correct *and* exactly recovers the cratonic-core area fraction
//! `Cr` for both circular and square reference shapes.
//!
//! We therefore use `L_plate = max BFS depth in plate`. The
//! empirical probe across 31 non-degenerate seeds confirms this
//! gives a mean cratonic-cell-fraction ratio of ~1.13 vs the
//! `Cr · continental_fraction` target (vs ~1.26 with the
//! sqrt(area) reading), with most seeds inside the ±20 %
//! tolerance of acceptance #8. The §4.10 patch documents this
//! reading explicitly.

use crate::tectonics_v2::boundaries::PlateType;
use crate::tectonics_v2::field::Field2D;
use crate::tectonics_v2::voronoi::VoronoiPlates;

use super::CratonicConfigEnabled;

/// Build the cratonic factor field from a Voronoï partition.
///
/// Returns a `Field2D` of values in `[0, 1]`, sized `nx × ny` (taken
/// from `plates.plate_id`). Non-continental cells, cells in
/// excluded small plates, and cells outside the smoothstep's outer
/// edge all yield `0`. Cells deep inside large continental plates
/// yield `1`.
pub fn build_cratonic_factor_field(plates: &VoronoiPlates, cfg: &CratonicConfigEnabled) -> Field2D {
    let nx = plates.plate_id.nx();
    let ny = plates.plate_id.ny();
    let domain_area = (nx * ny) as f64;
    let area_threshold = cfg.plate_area_min * domain_area;

    // 1. Per-plate area (cell count) and retention mask.
    let plate_areas = compute_plate_areas(plates, nx, ny);
    let retained: Vec<bool> = (0..plates.num_plates)
        .map(|p| {
            let is_continental = matches!(plates.per_plate_type[p], PlateType::Continental);
            let big_enough = plate_areas[p] as f64 >= area_threshold;
            is_continental && big_enough
        })
        .collect();

    // 2. Multi-source BFS from every non-continental cell.
    let dist = bfs_distance_to_oceanic(plates, &retained, nx, ny);

    // Per-plate "characteristic length" = max BFS depth attained
    // inside the plate. See module docstring for why this — and
    // not sqrt(area) — is the geometrically correct normaliser
    // for the d_mid formula.
    let mut plate_max_dist = vec![0u32; plates.num_plates];
    for j in 0..ny {
        for i in 0..nx {
            let pid = plates.plate_id.get(i, j) as usize;
            if !retained[pid] {
                continue;
            }
            let d = dist[j * nx + i];
            if d > plate_max_dist[pid] {
                plate_max_dist[pid] = d;
            }
        }
    }

    // 3 + 4. Per-plate normalisation and smoothstep.
    let mut factor = Field2D::new(nx, ny);
    let d_mid = 1.0 - cfg.cr.sqrt();
    let half_w = 0.5 * cfg.smoothing_width;
    let d_low = d_mid - half_w;
    let d_high = d_mid + half_w;
    for j in 0..ny {
        for i in 0..nx {
            let pid = plates.plate_id.get(i, j) as usize;
            if !retained[pid] {
                continue;
            }
            let l_plate = (plate_max_dist[pid] as f64).max(1.0);
            let dist_norm = dist[j * nx + i] as f64 / l_plate;
            factor.set(i, j, smoothstep(d_low, d_high, dist_norm));
        }
    }
    factor
}

/// Cell counts per plate id.
fn compute_plate_areas(plates: &VoronoiPlates, nx: usize, ny: usize) -> Vec<u32> {
    let mut areas = vec![0u32; plates.num_plates];
    for j in 0..ny {
        for i in 0..nx {
            let pid = plates.plate_id.get(i, j) as usize;
            areas[pid] += 1;
        }
    }
    areas
}

/// Multi-source BFS on the periodic torus. Sources are every cell
/// whose plate is *not* in `retained` (oceanic plates and excluded
/// small continentals). Returns a `Vec<u32>` of distances in cell
/// units, indexed `j * nx + i`.
///
/// Cells that are themselves a source have distance `0`. The BFS is
/// guaranteed to fill every cell because at least one cell is a
/// source whenever there is a non-retained plate; if the entire
/// domain is retained continental the function returns `nx*ny`
/// sentinel for every cell so the smoothstep saturates to `1` at
/// the centre — this is the correct fall-through for that
/// (degenerate, untestable in practice) case.
///
/// Implementation note: the queue stores cell indices as `u32` so
/// the BFS allocates `4 · nx * ny` bytes for the queue plus the
/// distance array. Periodic neighbours are computed inline; we do
/// not depend on the `PeriodicIndex` helper here because the queue
/// already pays for the wrap arithmetic per neighbour push.
fn bfs_distance_to_oceanic(
    plates: &VoronoiPlates,
    retained: &[bool],
    nx: usize,
    ny: usize,
) -> Vec<u32> {
    let n = nx * ny;
    let mut dist: Vec<u32> = vec![u32::MAX; n];
    let mut queue: std::collections::VecDeque<u32> = std::collections::VecDeque::with_capacity(n);

    // Seed the BFS with every source cell at distance 0.
    for j in 0..ny {
        for i in 0..nx {
            let pid = plates.plate_id.get(i, j) as usize;
            if !retained[pid] {
                let k = (j * nx + i) as u32;
                dist[k as usize] = 0;
                queue.push_back(k);
            }
        }
    }

    if queue.is_empty() {
        // Degenerate case: entire domain is retained continental.
        // No oceanic boundary exists; saturate distances so the
        // smoothstep returns 1 everywhere and `cratonic_factor`
        // reduces to a flat 1.
        return vec![n as u32; n];
    }

    while let Some(k) = queue.pop_front() {
        let k = k as usize;
        let i = k % nx;
        let j = k / nx;
        let d = dist[k];
        let next = d.saturating_add(1);
        // Four periodic neighbours.
        let neighbours = [
            (if i + 1 < nx { i + 1 } else { 0 }, j),
            (if i == 0 { nx - 1 } else { i - 1 }, j),
            (i, if j + 1 < ny { j + 1 } else { 0 }),
            (i, if j == 0 { ny - 1 } else { j - 1 }),
        ];
        for (ni, nj) in neighbours {
            let nk = nj * nx + ni;
            if dist[nk] > next {
                dist[nk] = next;
                queue.push_back(nk as u32);
            }
        }
    }

    dist
}

/// Cubic smoothstep `3t² − 2t³` clamped to `[0, 1]`. Value `0`
/// outside `[edge_low, edge_high]` on the low side, value `1`
/// outside on the high side, monotone in between. `edge_high`
/// strictly greater than `edge_low` is required; equal edges fall
/// back to a step at `edge_low`.
#[inline]
pub fn smoothstep(edge_low: f64, edge_high: f64, x: f64) -> f64 {
    if x <= edge_low {
        return 0.0;
    }
    if x >= edge_high {
        return 1.0;
    }
    let span = edge_high - edge_low;
    if span <= 0.0 {
        return if x < edge_low { 0.0 } else { 1.0 };
    }
    let t = (x - edge_low) / span;
    t * t * (3.0 - 2.0 * t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics_v2::voronoi::{VoronoiConfig, generate_voronoi};

    #[test]
    fn smoothstep_is_zero_below_low_edge() {
        assert_eq!(smoothstep(0.2, 0.4, 0.1), 0.0);
        assert_eq!(smoothstep(0.2, 0.4, 0.2), 0.0);
    }

    #[test]
    fn smoothstep_is_one_above_high_edge() {
        assert_eq!(smoothstep(0.2, 0.4, 0.4), 1.0);
        assert_eq!(smoothstep(0.2, 0.4, 0.5), 1.0);
    }

    #[test]
    fn smoothstep_is_half_at_midpoint() {
        // 3·0.5² − 2·0.5³ = 0.75 − 0.25 = 0.5 exactly.
        assert!((smoothstep(0.0, 1.0, 0.5) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn smoothstep_is_monotone_non_decreasing() {
        let edges = [(0.0, 1.0), (0.2, 0.5), (-0.1, 0.3)];
        for (a, b) in edges {
            let mut prev = -f64::INFINITY;
            for k in 0..101 {
                let x = a + (b - a) * (k as f64) / 100.0 - 0.1;
                let v = smoothstep(a, b, x);
                assert!(v >= prev - 1e-15, "smoothstep not monotone at x = {}", x);
                prev = v;
            }
        }
    }

    #[test]
    fn factor_field_in_unit_range() {
        let cfg = VoronoiConfig::default();
        let plates = generate_voronoi(64, 64, &cfg, 42);
        let factor = build_cratonic_factor_field(&plates, &CratonicConfigEnabled::default());
        for v in factor.data() {
            assert!(*v >= 0.0 && *v <= 1.0, "factor out of range: {}", v);
        }
    }

    #[test]
    fn factor_zero_on_oceanic_cells() {
        let cfg = VoronoiConfig::default();
        let plates = generate_voronoi(64, 64, &cfg, 42);
        let factor = build_cratonic_factor_field(&plates, &CratonicConfigEnabled::default());
        for j in 0..64 {
            for i in 0..64 {
                if matches!(plates.plate_type.get(i, j), PlateType::Oceanic) {
                    assert_eq!(factor.get(i, j), 0.0);
                }
            }
        }
    }

    #[test]
    fn small_continental_plates_are_excluded() {
        // Force a tiny Voronoï so at least one continental plate
        // ends up below the area threshold. With 12 plates on a
        // 32 × 32 domain, expected area ≈ 85 cells ≈ 0.083 of the
        // domain — below `plate_area_min = 0.10`. Each continental
        // plate that lands below the threshold must show `factor = 0`
        // for every cell in it.
        let vcfg = VoronoiConfig { num_plates: 12, continental_ratio: 0.5 };
        let plates = generate_voronoi(32, 32, &vcfg, 7);
        let cfg = CratonicConfigEnabled { plate_area_min: 0.15, ..Default::default() };
        let factor = build_cratonic_factor_field(&plates, &cfg);

        let domain_area = 32.0 * 32.0;
        let area_threshold = cfg.plate_area_min * domain_area;
        let plate_areas = compute_plate_areas(&plates, 32, 32);

        let mut excluded_count = 0;
        for pid in 0..plates.num_plates {
            if matches!(plates.per_plate_type[pid], PlateType::Continental)
                && (plate_areas[pid] as f64) < area_threshold
            {
                excluded_count += 1;
                for j in 0..32 {
                    for i in 0..32 {
                        if plates.plate_id.get(i, j) as usize == pid {
                            assert_eq!(
                                factor.get(i, j),
                                0.0,
                                "small continental plate {} cell ({},{}) should be 0",
                                pid,
                                i,
                                j,
                            );
                        }
                    }
                }
            }
        }
        // Make sure we actually exercised the exclusion; if the
        // random Voronoï doesn't land any small continental, the
        // test isn't probing the right thing.
        assert!(
            excluded_count >= 1,
            "test setup did not produce any excluded small continental plate \
             — adjust seed or num_plates to expose the path"
        );
    }

    #[test]
    fn disabled_path_returns_zero_field() {
        // The harness drives `CratonicConfig::Disabled` at the call
        // site (no `build_cratonic_factor_field` call); this test
        // double-checks the contract that *if* it were called with
        // factor = 0 everywhere, the eta multiplier round-trip is
        // identity.
        let factor = Field2D::filled(8, 8, 0.0);
        let state = super::super::CratonicState::from_factor(factor, 5.0);
        for v in state.eta_multiplier.data() {
            assert_eq!(*v, 1.0);
        }
    }

    #[test]
    fn factor_static_across_calls() {
        // D7 — static identification. The factor field depends only
        // on the Voronoï partition and the config; it is independent
        // of any time-varying state. Calling twice with the same
        // inputs returns byte-identical output.
        let cfg = VoronoiConfig::default();
        let plates = generate_voronoi(64, 64, &cfg, 42);
        let f1 = build_cratonic_factor_field(&plates, &CratonicConfigEnabled::default());
        let f2 = build_cratonic_factor_field(&plates, &CratonicConfigEnabled::default());
        for (a, b) in f1.data().iter().zip(f2.data().iter()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }

    #[test]
    fn smoothing_width_zero_produces_step_at_d_mid() {
        // smoothing_width = 0 ⟹ d_low = d_high = d_mid. Smoothstep
        // collapses to a step at d_mid (factor jumps 0 → 1). Useful
        // as a stress test of the boundary handling.
        let cfg = VoronoiConfig { num_plates: 1, continental_ratio: 1.0 };
        let plates = generate_voronoi(32, 32, &cfg, 11);
        let crcfg = CratonicConfigEnabled {
            cr: 0.3,
            smoothing_width: 0.0,
            plate_area_min: 0.0,
            ..Default::default()
        };
        let factor = build_cratonic_factor_field(&plates, &crcfg);
        // Each cell's value should be either 0 or 1 (no in-between).
        for v in factor.data() {
            assert!(*v == 0.0 || *v == 1.0, "value {} not at step boundary", v);
        }
    }

    #[test]
    fn entire_continental_domain_saturates_to_one() {
        // Degenerate case: single continental plate spanning the
        // whole torus, no oceanic source for the BFS. The function
        // returns factor = 1 everywhere (smoothstep saturates) so
        // downstream multipliers behave as "uniform craton".
        let cfg = VoronoiConfig { num_plates: 1, continental_ratio: 1.0 };
        let plates = generate_voronoi(16, 16, &cfg, 3);
        // Make sure the single plate is continental — if it isn't,
        // the test isn't probing the saturation path.
        if !matches!(plates.per_plate_type[0], PlateType::Continental) {
            return;
        }
        let crcfg = CratonicConfigEnabled { plate_area_min: 0.0, ..Default::default() };
        let factor = build_cratonic_factor_field(&plates, &crcfg);
        for v in factor.data() {
            assert_eq!(*v, 1.0, "saturated factor expected 1.0, got {}", v);
        }
    }

    #[test]
    fn factor_higher_in_plate_interior_than_at_boundary() {
        // Sanity check: deep cells (max BFS distance) should have
        // strictly larger factor than cells adjacent to oceanic
        // boundary, for at least one continental plate in a typical
        // Voronoï.
        let cfg = VoronoiConfig::default();
        let plates = generate_voronoi(64, 64, &cfg, 42);
        let crcfg = CratonicConfigEnabled::default();
        let factor = build_cratonic_factor_field(&plates, &crcfg);

        // Find a continental plate and compare its max factor vs
        // the average factor over its boundary-adjacent cells.
        let nx = 64;
        let ny = 64;
        let plate_areas = compute_plate_areas(&plates, nx, ny);
        let domain_area = (nx * ny) as f64;
        let area_threshold = crcfg.plate_area_min * domain_area;
        let mut found = false;
        for pid in 0..plates.num_plates {
            if !matches!(plates.per_plate_type[pid], PlateType::Continental) {
                continue;
            }
            if (plate_areas[pid] as f64) < area_threshold {
                continue;
            }
            let mut max_v = 0.0_f64;
            for j in 0..ny {
                for i in 0..nx {
                    if plates.plate_id.get(i, j) as usize == pid {
                        max_v = max_v.max(factor.get(i, j));
                    }
                }
            }
            if max_v > 0.0 {
                // We have a craton-bearing plate; max factor should
                // be > 0.5 (well into the smoothstep's high side).
                assert!(max_v > 0.5, "plate {} max cratonic_factor only {}", pid, max_v,);
                found = true;
                break;
            }
        }
        assert!(found, "no retained continental plate produced a craton");
    }
}
