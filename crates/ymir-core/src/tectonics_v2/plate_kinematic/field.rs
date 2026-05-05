//! Step 11 — Plate kinematic drift field construction.
//!
//! Builds the `(vx_drift, vy_drift)` row-major buffers added by the
//! harness to every solver output post-solve (`v_total = v_solver +
//! v_drift`) from a per-plate velocity assignment + a Voronoï plate
//! map.
//! The shape follows the [`super::super::init::InitMode::Uniform`]
//! pattern introduced at Step 8.6 Phase 8a:
//!
//! 1. Set per-cell own velocity from `velocities[plate_id[i, j]]`.
//! 2. Seed a BFS at every cell whose 4-neighbour (NESW, periodic)
//!    belongs to a different plate. Distance starts at 0; target
//!    velocity = the across-boundary neighbour's plate velocity
//!    (deterministic — first match in NESW order).
//! 3. Chebyshev BFS (8-neighbour, periodic). Each hop adds `+1` to
//!    the distance.
//! 4. For each cell, blend the own velocity with the *midpoint*
//!    between own and across-boundary neighbour with cubic
//!    smoothstep weight on `d / boundary_smoothing_width`. Outside
//!    the band (`d >= width`), the value is the own velocity
//!    exactly.
//!
//! ## Why midpoint, not direct neighbour, at the boundary
//!
//! The issue's literal smoothstep formula
//! `v = weight_a · vel[p_a] + weight_b · vel[p_b]` (with
//! `weight_a = smoothstep(t)`) produces `v = vel[p_b]` exactly at the
//! boundary cell on plate A's side, and `v = vel[p_a]` on plate B's
//! side — i.e. an *inverted* discontinuity at the boundary itself,
//! which contradicts the smoothness contract from D8 ("if smoothstep
//! interpolation produces visible discontinuities … diagnose").
//!
//! Switching to the midpoint convention (matches `init::Uniform`):
//!
//! ```text
//! result = own · st + ((own + other) / 2) · (1 - st)
//! ```
//!
//! gives `result = midpoint` on both sides at `d = 0`, so the field
//! is *continuous* across the boundary, and `result = own` at
//! `d = width`. Same overall family of behaviour, but smooth by
//! construction.

use super::super::voronoi::{compute_dist_to_inter_plate_boundary, PlateIdField};

/// Build the `(vx, vy)` initial velocity buffers from a per-plate
/// assignment.
///
/// Returns row-major `Vec<f64>` buffers of length `nx · ny`.
///
/// # Panics
///
/// - `boundary_smoothing_width` not finite or `≤ 0`.
/// - Some `plate_id[i, j]` falls outside `0..velocities.len()`.
///
/// # Determinism
///
/// Output is deterministic for fixed `(plate_id, velocities,
/// boundary_smoothing_width)` — the BFS seeds in row-major order and
/// neighbour ties in BFS pick the first NESW match deterministically.
pub fn build(
    nx: usize,
    ny: usize,
    plate_id: &PlateIdField,
    velocities: &[(f64, f64)],
    boundary_smoothing_width: f64,
) -> (Vec<f64>, Vec<f64>) {
    assert!(
        boundary_smoothing_width.is_finite() && boundary_smoothing_width > 0.0,
        "boundary_smoothing_width must be a positive finite scalar; got {}",
        boundary_smoothing_width
    );
    assert!(
        !velocities.is_empty(),
        "velocities must not be empty (one entry per plate)",
    );

    // Step 13 Phase 1: BFS distance computation delegated to the
    // shared `compute_dist_to_inter_plate_boundary` utility (also
    // used by `init::Uniform` and `init::radial_profile`).
    // Bit-identical with the pre-refactor implementation by
    // construction: per-plate properties (here the (vx, vy) entry of
    // `velocities`) are constant within a plate, so propagating
    // `target_plate_id` through the BFS and indexing `velocities`
    // at the end is equivalent to propagating `(target_vx,
    // target_vy)` directly. See
    // `tectonics_v2::voronoi::distance` module docstring for why
    // `cratonic::factor` is *not* unified into the same utility
    // (different connectivity, different sources, different
    // physical meaning).
    let bfs = compute_dist_to_inter_plate_boundary(nx, ny, plate_id);

    let n = nx * ny;
    let mut own_vx = vec![0.0_f64; n];
    let mut own_vy = vec![0.0_f64; n];

    // Per-cell own velocity from plate_id lookup. Bounds-checked so
    // an out-of-range plate_id panics early with a clear message
    // (would otherwise silently corrupt the field).
    for j in 0..ny {
        for i in 0..nx {
            let pid = plate_id.get(i, j) as usize;
            let (vx, vy) = *velocities.get(pid).unwrap_or_else(|| {
                panic!(
                    "plate_id[{}, {}] = {} out of range for velocities of length {}",
                    i, j, pid, velocities.len()
                )
            });
            own_vx[j * nx + i] = vx;
            own_vy[j * nx + i] = vy;
        }
    }

    let w = boundary_smoothing_width;
    let mut out_vx = vec![0.0_f64; n];
    let mut out_vy = vec![0.0_f64; n];
    for j in 0..ny {
        for i in 0..nx {
            let idx = j * nx + i;
            let d = bfs.distance.get(i, j);
            let own_x = own_vx[idx];
            let own_y = own_vy[idx];
            let tpid = bfs.target_plate_id[idx];
            let (other_x, other_y) = if tpid == u16::MAX {
                // BFS never reached this cell (degenerate single-
                // plate-on-torus). dist=INFINITY → the `d >= w`
                // branch short-circuits to `own`; the value picked
                // here is unused. Defensive default = own.
                (own_x, own_y)
            } else {
                velocities[tpid as usize]
            };
            if d >= w {
                out_vx[idx] = own_x;
                out_vy[idx] = own_y;
            } else {
                let t = (d / w).clamp(0.0, 1.0);
                let st = t * t * (3.0 - 2.0 * t);
                let mid_x = 0.5 * (own_x + other_x);
                let mid_y = 0.5 * (own_y + other_y);
                out_vx[idx] = own_x * st + mid_x * (1.0 - st);
                out_vy[idx] = own_y * st + mid_y * (1.0 - st);
            }
        }
    }

    (out_vx, out_vy)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics_v2::voronoi::{generate_voronoi, VoronoiConfig};
    use std::collections::VecDeque;

    fn build_test_plates(nx: usize, ny: usize, seed: u64, num_plates: usize) {
        let _ = (nx, ny, seed, num_plates);
    }

    fn make_plates(nx: usize, ny: usize, seed: u64, num_plates: usize) -> crate::tectonics_v2::voronoi::VoronoiPlates {
        build_test_plates(nx, ny, seed, num_plates);
        let cfg = VoronoiConfig { num_plates, continental_ratio: 0.3 };
        generate_voronoi(nx, ny, &cfg, seed)
    }

    /// All-zero per-plate input → all-zero output. Sanity check that
    /// the algorithm produces identically zero velocities when the
    /// inputs are zero, regardless of plate geometry. (Bit-identity
    /// to `Zero` is enforced at the harness level by the structural
    /// short-circuit; here we only verify the algorithm doesn't
    /// introduce numerical noise from its own computation.)
    #[test]
    fn all_zero_input_produces_zero_field() {
        let nx = 32;
        let ny = 32;
        let plates = make_plates(nx, ny, 42, 6);
        let velocities = vec![(0.0, 0.0); plates.num_plates];
        let (vx, vy) = build(nx, ny, &plates.plate_id, &velocities, 1.5);
        assert_eq!(vx, vec![0.0_f64; nx * ny]);
        assert_eq!(vy, vec![0.0_f64; nx * ny]);
    }

    /// Cells far from any inter-plate boundary
    /// (`dist > boundary_smoothing_width`) hold their plate's
    /// velocity exactly. Acceptance criterion #2 from the issue.
    #[test]
    fn interior_uniform() {
        let nx = 64;
        let ny = 64;
        let plates = make_plates(nx, ny, 42, 4);
        // Non-trivial per-plate velocities so a wrong assignment is
        // distinguishable.
        let velocities: Vec<(f64, f64)> = (0..plates.num_plates)
            .map(|p| (0.1 * (p as f64 + 1.0), -0.05 * (p as f64 + 1.0)))
            .collect();
        let width = 1.5;
        let (vx, vy) = build(nx, ny, &plates.plate_id, &velocities, width);

        // Replicate the BFS distance to know which cells qualify
        // as "interior" (`d >= width`).
        let mut dist = vec![f64::INFINITY; nx * ny];
        let mut q: VecDeque<(usize, usize)> = VecDeque::new();
        for j in 0..ny {
            for i in 0..nx {
                let id = plates.plate_id.get(i, j);
                let ip = (i + 1) % nx;
                let im = (i + nx - 1) % nx;
                let jp = (j + 1) % ny;
                let jm = (j + ny - 1) % ny;
                if [(ip, j), (im, j), (i, jp), (i, jm)]
                    .iter()
                    .any(|&(ni, nj)| plates.plate_id.get(ni, nj) != id)
                {
                    dist[j * nx + i] = 0.0;
                    q.push_back((i, j));
                }
            }
        }
        while let Some((i, j)) = q.pop_front() {
            let d = dist[j * nx + i];
            for dj in [-1_i32, 0, 1] {
                for di in [-1_i32, 0, 1] {
                    if di == 0 && dj == 0 {
                        continue;
                    }
                    let ni = ((i as i32 + di).rem_euclid(nx as i32)) as usize;
                    let nj = ((j as i32 + dj).rem_euclid(ny as i32)) as usize;
                    let nidx = nj * nx + ni;
                    let nd = d + 1.0;
                    if nd < dist[nidx] {
                        dist[nidx] = nd;
                        q.push_back((ni, nj));
                    }
                }
            }
        }

        let mut interior_count = 0usize;
        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;
                if dist[idx] >= width {
                    let pid = plates.plate_id.get(i, j) as usize;
                    let (ex_vx, ex_vy) = velocities[pid];
                    assert!(
                        (vx[idx] - ex_vx).abs() < 1e-15,
                        "interior cell ({},{}) plate {}: vx {} != expected {}",
                        i, j, pid, vx[idx], ex_vx
                    );
                    assert!(
                        (vy[idx] - ex_vy).abs() < 1e-15,
                        "interior cell ({},{}) plate {}: vy {} != expected {}",
                        i, j, pid, vy[idx], ex_vy
                    );
                    interior_count += 1;
                }
            }
        }
        assert!(
            interior_count > 0,
            "no interior cells found in 64² with 4 plates and width=1.5 — \
             grid likely too small or width too wide"
        );
    }

    /// Across any boundary the cell-to-cell velocity step is bounded
    /// by the smoothstep slope. With `width=1.5`, the discrete
    /// midpoint convention gives a per-cell delta ≤ 0.5·|v_a - v_b|
    /// at worst (the midpoint vs. own jump at `d=0`). We test a
    /// looser bound (= |v_a - v_b|) to leave headroom for the BFS
    /// "first NESW neighbour" determinism choice.
    #[test]
    fn boundary_smoothing_bounded_step() {
        let nx = 64;
        let ny = 64;
        let plates = make_plates(nx, ny, 42, 4);
        let velocities: Vec<(f64, f64)> = (0..plates.num_plates)
            .map(|p| (0.5 * ((p % 2) as f64) - 0.25, 0.0))
            .collect();
        let max_delta_per_component = velocities
            .iter()
            .flat_map(|&(a, b)| [a, b])
            .fold(0.0_f64, f64::max)
            - velocities
                .iter()
                .flat_map(|&(a, b)| [a, b])
                .fold(0.0_f64, f64::min);

        let (vx, vy) = build(nx, ny, &plates.plate_id, &velocities, 1.5);

        let mut max_step = 0.0_f64;
        for j in 0..ny {
            for i in 0..nx {
                let here = (vx[j * nx + i], vy[j * nx + i]);
                let ip = (i + 1) % nx;
                let jp = (j + 1) % ny;
                let right = (vx[j * nx + ip], vy[j * nx + ip]);
                let down = (vx[jp * nx + i], vy[jp * nx + i]);
                let dx_r = (here.0 - right.0).abs().max((here.1 - right.1).abs());
                let dx_d = (here.0 - down.0).abs().max((here.1 - down.1).abs());
                max_step = max_step.max(dx_r).max(dx_d);
            }
        }
        // Without smoothing the step at a boundary would be the
        // raw `|v_a - v_b|`. Smoothstep midpoint convention halves
        // that (the boundary itself sits at the midpoint), then the
        // 1-cell hop into the interior advances by another half.
        // We assert the step stays at or below the unsmoothed
        // jump — i.e. smoothing didn't make things worse.
        assert!(
            max_step <= max_delta_per_component + 1e-12,
            "max cell-to-cell step {} exceeds raw inter-plate jump {}",
            max_step,
            max_delta_per_component
        );
    }

    /// Max magnitude of the constructed field is bounded by the max
    /// magnitude over the input plate velocities (no overshoot from
    /// the smoothstep blend). Acceptance criterion #4.
    #[test]
    fn no_overshoot() {
        let nx = 64;
        let ny = 64;
        let plates = make_plates(nx, ny, 42, 5);
        let velocities = vec![
            (0.5, 0.0),
            (-0.5, 0.0),
            (0.0, 0.5),
            (0.0, -0.5),
            (0.3, 0.4), // |v| = 0.5
        ];
        let max_input_mag = velocities
            .iter()
            .map(|&(a, b): &(f64, f64)| (a * a + b * b).sqrt())
            .fold(0.0_f64, f64::max);

        let (vx, vy) = build(nx, ny, &plates.plate_id, &velocities, 1.5);

        let max_out_mag = vx
            .iter()
            .zip(vy.iter())
            .map(|(&a, &b)| (a * a + b * b).sqrt())
            .fold(0.0_f64, f64::max);

        assert!(
            max_out_mag <= max_input_mag + 1e-12,
            "max output magnitude {} > max input magnitude {} (overshoot)",
            max_out_mag,
            max_input_mag
        );
    }

    /// Determinism: same `(plate_id, velocities, width)` → identical
    /// `(vx, vy)` byte-for-byte.
    #[test]
    fn deterministic_same_inputs() {
        let nx = 32;
        let ny = 32;
        let plates_a = make_plates(nx, ny, 42, 4);
        let plates_b = make_plates(nx, ny, 42, 4);
        let velocities = vec![(0.5, 0.0), (-0.5, 0.0), (0.0, 0.3), (0.0, -0.3)];
        let (ax, ay) = build(nx, ny, &plates_a.plate_id, &velocities, 1.5);
        let (bx, by) = build(nx, ny, &plates_b.plate_id, &velocities, 1.5);
        assert_eq!(ax, bx);
        assert_eq!(ay, by);
    }
}
