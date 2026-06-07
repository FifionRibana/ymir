//! Per-time-step application of the stream-power erosion
//! closure.
//!
//! See the parent module ([`super`]) for the physics derivation
//! (Whipple-Tucker 1999 SPIM with Lague 2014 calibration
//! discipline), the interaction with Phase 1.3's equilibrium-
//! height sink, and the rationale for eroding `S̃` directly
//! (implicit Airy isostatic compensation).

use crate::grid::GridF32;
use crate::tectonics_v2::field::Field2D;
use crate::tectonics_v2::workflow::drainage::DrainageMap;

use super::params::ErosionParams;

/// Apply one forward-Euler step of the stream-power erosion
/// closure to the `S̃` field, **in place**.
///
/// For each cell `c`:
///
/// ```text
///     S̃_new(c) = max(floor,
///                    S̃(c) − k · A(c)^m · |∇h(c)|^n · dt )
/// ```
///
/// where `A(c)` is the transitive drainage area at cell `c`
/// (counts upstream cells routing through `c`, including `c`
/// itself), and `|∇h(c)|` is the centered-difference slope
/// magnitude of the altitude heightmap with periodic wraparound.
///
/// Cells with `slope ≤ 0` (flat / sink) or `A ≤ 0` are skipped —
/// no erosion where no transport gradient or no flow accumulates.
///
/// Inputs:
/// - `s`: mutable `S̃` field, updated in place. `f64`-backed.
/// - `altitude`: post-isostasy heightmap. `f32`-backed
///   [`GridF32`], converted to `f64` per cell for the gradient
///   computation (precision loss is bounded by the heightmap's
///   `f32` resolution — fine for slope magnitudes well above
///   `1e-7`).
/// - `drainage_areas`: per-cell transitive drainage areas
///   (flat row-major slice, length `nx * ny`). Computed via the
///   private `compute_drainage_areas` helper from a
///   [`DrainageMap`] that the per-step time loop already needs
///   for other downstream passes.
/// - `params`: closure tunables (see [`ErosionParams`]).
/// - `dt`: time step.
/// - `dx`: cell size (non-dim; `1.0 / nx` at unit-domain).
pub fn apply_erosion_step(
    s: &mut Field2D,
    altitude: &GridF32,
    drainage_areas: &[u32],
    params: &ErosionParams,
    dt: f64,
    dx: f64,
) {
    if !params.enabled {
        return;
    }
    let nx = s.nx();
    let ny = s.ny();
    debug_assert_eq!(
        drainage_areas.len(),
        nx * ny,
        "drainage_areas length must match grid (nx · ny)"
    );
    debug_assert_eq!(
        altitude.width, nx,
        "altitude width must match S̃ grid"
    );
    debug_assert_eq!(
        altitude.height, ny,
        "altitude height must match S̃ grid"
    );

    let k = params.k;
    let m = params.m;
    let n_exp = params.n;
    let floor = params.floor;

    for j in 0..ny {
        for i in 0..nx {
            let slope = compute_local_slope(altitude, i, j, dx);
            if slope <= 0.0 {
                continue;
            }
            let a = drainage_areas[j * nx + i] as f64;
            if a <= 0.0 {
                continue;
            }
            let erosion_rate = k * a.powf(m) * slope.powf(n_exp);
            let delta = erosion_rate * dt;
            let s_old = s.get(i, j);
            // Issue #145 — clean non-injecting removal. The legacy
            // `max(floor, s_old - delta)` RAISED any cell already below
            // `floor` up to `floor`, injecting mass (+247% standalone on
            // rigid transport; the known Phase-1.4 non-conservation). Erode
            // toward `floor` but NEVER above `s_old`. The ONLY behavioural
            // change vs legacy is sub-floor cells (s_old < floor): legacy
            // injected them up to floor; we leave them untouched. For
            // s_old >= floor the result is identical (floor returned
            // exactly, so no FP drift vs the old `max`). Conservative
            // deposition (sediment routing) is a separate follow-up, NOT #145.
            let s_new = if s_old <= floor {
                s_old // already at/below floor — never inject
            } else if delta >= s_old - floor {
                floor // fully eroded to the floor (exact, no FP drift)
            } else {
                s_old - delta
            };
            s.set(i, j, s_new);
        }
    }
}

/// Centered-difference slope magnitude on the altitude heightmap,
/// with periodic wraparound. Returns `|∇h|` in non-dim altitude /
/// non-dim length units.
///
/// `GridF32::get(x, y)` clamps to `0.0` out-of-bounds; the
/// periodic indices computed here always sit in-bounds so the
/// clamping never triggers. The `f32 → f64` casts at the four
/// neighbour samples introduce no measurable error for the slope
/// magnitudes the erosion closure responds to (≥ `1e-6`).
fn compute_local_slope(altitude: &GridF32, i: usize, j: usize, dx: f64) -> f64 {
    let nx = altitude.width;
    let ny = altitude.height;
    let im1 = if i == 0 { nx - 1 } else { i - 1 };
    let ip1 = if i == nx - 1 { 0 } else { i + 1 };
    let jm1 = if j == 0 { ny - 1 } else { j - 1 };
    let jp1 = if j == ny - 1 { 0 } else { j + 1 };

    let h_e = altitude.get(ip1 as i32, j as i32) as f64;
    let h_w = altitude.get(im1 as i32, j as i32) as f64;
    let h_n = altitude.get(i as i32, jp1 as i32) as f64;
    let h_s = altitude.get(i as i32, jm1 as i32) as f64;

    let dz_dx = (h_e - h_w) / (2.0 * dx);
    let dz_dy = (h_n - h_s) / (2.0 * dx);

    (dz_dx * dz_dx + dz_dy * dz_dy).sqrt()
}

/// Compute per-cell **transitive** drainage areas from a
/// [`DrainageMap`].
///
/// W-T's `A` is the number of upstream cells whose flow paths
/// pass through this cell (including itself). The v2
/// [`compute_drainage_targets`](crate::tectonics_v2::workflow::drainage::compute_drainage_targets)
/// produces per-cell drainage *targets* (where each cell drains
/// to), not areas. This function inverts targets into transitive
/// areas in `O(N log N)` via path-length-descending accumulation:
///
/// 1. Sort cell indices by `path_length` descending. Cells
///    further from their final destination are visited first.
/// 2. For each cell `i` in sorted order, add its accumulated
///    area to its target: `areas[target_idx[i]] += areas[i]`.
///    The guard `target != i` skips sinks (oceanic cells with
///    `target_idx[i] = i`).
///
/// Correctness: a cell with `path_length = k` is only ever
/// accumulated INTO by cells with `path_length ≥ k + 1` (their
/// drainage paths route through it). Descending-order processing
/// ensures all upstream accumulations land before downstream
/// accumulation reads them. The
/// `drainage_area_transitive_accumulation_linear_chain` test
/// locks this against the naive in-degree-only variant.
///
/// An `O(N + max_path_length)` counting-sort variant is possible
/// (bucket by `path_length: u8`), but at 64²×4096 cells the sort
/// cost is sub-microsecond and the simpler `sort_by_key` is
/// retained. Promote if Phase 2+ profiling identifies this as a
/// hot path.
///
/// Crate-visible (`pub(crate)`) since Phase 1.4 Stage E2: the
/// C1 time loop consumes it to prepare the `drainage_areas`
/// argument for [`apply_erosion_step`]. Promote to `pub`
/// (potentially in `workflow::drainage`) if Phase 2+ needs it
/// from outside `ymir-core`.
pub(crate) fn compute_drainage_areas(map: &DrainageMap) -> Vec<u32> {
    let n = map.target_idx.len();
    let mut areas = vec![1u32; n];

    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by_key(|&i| std::cmp::Reverse(map.path_length[i]));

    for i in indices {
        let target = map.target_idx[i];
        if target != i {
            areas[target] += areas[i];
        }
    }
    areas
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flat_altitude(nx: usize, ny: usize, value: f32) -> GridF32 {
        GridF32 {
            width: nx,
            height: ny,
            data: vec![value; nx * ny],
        }
    }

    fn ramp_altitude_x(nx: usize, ny: usize) -> GridF32 {
        // Altitude increases linearly with i: h[i, j] = i / nx.
        // Slope magnitude = 1 (after dx = 1/nx normalisation).
        let mut data = vec![0.0f32; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                data[j * nx + i] = i as f32 / nx as f32;
            }
        }
        GridF32 {
            width: nx,
            height: ny,
            data,
        }
    }

    fn filled_field(nx: usize, ny: usize, value: f64) -> Field2D {
        let mut f = Field2D::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                f.set(i, j, value);
            }
        }
        f
    }

    /// Sanity: flat altitude → slope = 0 → no erosion.
    #[test]
    fn no_erosion_when_slope_zero() {
        let nx = 4;
        let ny = 4;
        let mut s = filled_field(nx, ny, 1.0);
        let altitude = flat_altitude(nx, ny, 0.5);
        let drainage_areas = vec![10u32; nx * ny];
        let params = ErosionParams::default();
        let before = s.data().to_vec();
        apply_erosion_step(
            &mut s,
            &altitude,
            &drainage_areas,
            &params,
            1.0,
            1.0 / nx as f64,
        );
        for k in 0..before.len() {
            assert_eq!(
                before[k],
                s.data()[k],
                "flat altitude → no erosion; mismatch at flat index {k}"
            );
        }
    }

    /// Sanity: drainage area zero → no erosion, even on a slope.
    #[test]
    fn no_erosion_when_drainage_zero() {
        let nx = 4;
        let ny = 4;
        let mut s = filled_field(nx, ny, 1.0);
        let altitude = ramp_altitude_x(nx, ny);
        let drainage_areas = vec![0u32; nx * ny]; // every cell A = 0
        let params = ErosionParams::default();
        let before = s.data().to_vec();
        apply_erosion_step(
            &mut s,
            &altitude,
            &drainage_areas,
            &params,
            1.0,
            1.0 / nx as f64,
        );
        for k in 0..before.len() {
            assert_eq!(
                before[k],
                s.data()[k],
                "A = 0 → no erosion; mismatch at flat index {k}"
            );
        }
    }

    /// W-T eq. (1) validation: doubling `A` increases erosion by
    /// `2^m = 2^0.5 ≈ 1.414` for the default `m = 0.5`.
    #[test]
    fn erosion_scales_with_drainage_sqrt() {
        let nx = 4;
        let ny = 4;
        let altitude = ramp_altitude_x(nx, ny);
        let dx = 1.0 / nx as f64;
        let dt = 1.0;
        let params = ErosionParams { k: 0.1, ..ErosionParams::default() };

        // Run 1 — A = 10
        let mut s1 = filled_field(nx, ny, 1.0);
        let areas1 = vec![10u32; nx * ny];
        apply_erosion_step(&mut s1, &altitude, &areas1, &params, dt, dx);

        // Run 2 — A = 20
        let mut s2 = filled_field(nx, ny, 1.0);
        let areas2 = vec![20u32; nx * ny];
        apply_erosion_step(&mut s2, &altitude, &areas2, &params, dt, dx);

        // Sample an interior cell (avoid the periodic wrap edge
        // where the ramp creates a discontinuity).
        let (i, j) = (2, 2);
        let delta1 = 1.0 - s1.get(i, j);
        let delta2 = 1.0 - s2.get(i, j);
        let ratio = delta2 / delta1;
        let expected = 2.0_f64.powf(params.m); // 2^0.5 ≈ 1.4142
        assert!(
            (ratio - expected).abs() / expected < 1e-9,
            "erosion ratio for A doubling: got {ratio:.4}, expected 2^m = {expected:.4} (m = {})",
            params.m
        );
    }

    /// W-T eq. (1) validation: doubling slope increases erosion
    /// linearly for the default `n = 1`.
    #[test]
    fn erosion_scales_linearly_with_slope() {
        let nx = 4;
        let ny = 4;
        let dx = 1.0 / nx as f64;
        let dt = 1.0;
        let params = ErosionParams { k: 0.1, ..ErosionParams::default() };

        // Slope 1× via the standard ramp.
        let altitude1 = ramp_altitude_x(nx, ny);
        let mut s1 = filled_field(nx, ny, 1.0);
        let areas = vec![10u32; nx * ny];
        apply_erosion_step(&mut s1, &altitude1, &areas, &params, dt, dx);

        // Slope 2× via a doubled-amplitude ramp.
        let mut data2 = vec![0.0f32; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                data2[j * nx + i] = 2.0 * (i as f32 / nx as f32);
            }
        }
        let altitude2 = GridF32 { width: nx, height: ny, data: data2 };
        let mut s2 = filled_field(nx, ny, 1.0);
        apply_erosion_step(&mut s2, &altitude2, &areas, &params, dt, dx);

        let (i, j) = (2, 2);
        let delta1 = 1.0 - s1.get(i, j);
        let delta2 = 1.0 - s2.get(i, j);
        let ratio = delta2 / delta1;
        let expected = 2.0_f64.powf(params.n); // 2^1 = 2
        assert!(
            (ratio - expected).abs() / expected < 1e-9,
            "erosion ratio for slope doubling: got {ratio:.4}, expected 2^n = {expected:.4} (n = {})",
            params.n
        );
    }

    /// `enabled = false` → bit-identical pre/post regardless of
    /// inputs. Mirrors the W4 closure-isolation discipline used
    /// by Davis-Suppe and equilibrium-height.
    #[test]
    fn disabled_no_op() {
        let nx = 4;
        let ny = 4;
        let mut s = filled_field(nx, ny, 5.0);
        // Mix in some non-uniform values so a forgotten branch
        // can't pass by accident on a uniform field.
        s.set(0, 0, 100.0);
        s.set(1, 2, -3.0);
        s.set(3, 3, 0.5);
        let altitude = ramp_altitude_x(nx, ny);
        let drainage_areas = vec![50u32; nx * ny];
        let params = ErosionParams {
            enabled: false,
            k: 100.0, // pathologically large but disabled
            ..ErosionParams::default()
        };
        let before = s.data().to_vec();
        apply_erosion_step(
            &mut s,
            &altitude,
            &drainage_areas,
            &params,
            1.0,
            1.0 / nx as f64,
        );
        for k in 0..before.len() {
            assert_eq!(
                before[k],
                s.data()[k],
                "disabled closure must not touch any cell; mismatch at flat index {k}"
            );
        }
    }

    /// Pathological `k · A^m · S^n · dt` would drive `S̃`
    /// arbitrarily negative without the floor clamp. The floor
    /// must hold `S̃` at the oceanic baseline.
    ///
    /// This test locks the defensive clamp; in normal use the
    /// product is `≪ S̃` and the floor never triggers.
    #[test]
    fn floor_at_oceanic_baseline() {
        let nx = 4;
        let ny = 4;
        let mut s = filled_field(nx, ny, 1.0); // continental
        let altitude = ramp_altitude_x(nx, ny);
        let drainage_areas = vec![1_000u32; nx * ny];
        let params = ErosionParams {
            k: 100.0, // pathological, would predict S̃ < 0 without floor
            ..ErosionParams::default()
        };
        let dx = 1.0 / nx as f64;
        // Verify the test premise: unclamped formula does predict undershoot.
        let predicted = 1.0
            - params.k
                * (1000.0_f64).powf(params.m)
                * (1.0_f64).powf(params.n)
                * 1.0;
        assert!(
            predicted < params.floor,
            "test premise: unclamped prediction {predicted} must undershoot floor {}",
            params.floor
        );
        apply_erosion_step(
            &mut s,
            &altitude,
            &drainage_areas,
            &params,
            1.0,
            dx,
        );
        // Every interior cell with slope > 0 should be clamped
        // at the floor (the ramp gives slope ≈ 1 everywhere).
        let (i, j) = (2, 2);
        assert_eq!(
            s.get(i, j),
            params.floor,
            "floor clamp must hold S̃ at oceanic baseline; got {} expected {}",
            s.get(i, j),
            params.floor
        );
    }

    /// Transitive drainage-area accumulation validation. Builds a
    /// 4-cell linear chain `A → B → C → D` and verifies the
    /// `compute_drainage_areas` algorithm produces
    /// `[1, 2, 3, 4]`, not the in-degree-only `[1, 2, 2, 2]`.
    ///
    /// This is the algorithm-correctness lock: any future
    /// "simplification" of `compute_drainage_areas` to a single
    /// `areas[target] += 1` pass will fail this test.
    #[test]
    fn drainage_area_transitive_accumulation_linear_chain() {
        // Cells laid out at flat indices 0, 1, 2, 3:
        //   0 (A) drains to 1 (B), path_length = 3
        //   1 (B) drains to 2 (C), path_length = 2
        //   2 (C) drains to 3 (D), path_length = 1
        //   3 (D) drains to itself (sink), path_length = 0
        let map = DrainageMap {
            target_idx: vec![1, 2, 3, 3],
            path_length: vec![3, 2, 1, 0],
        };
        let areas = compute_drainage_areas(&map);
        assert_eq!(
            areas,
            vec![1, 2, 3, 4],
            "transitive accumulation must produce [1, 2, 3, 4], not in-degree [1, 1, 1, 2]"
        );
    }
}
