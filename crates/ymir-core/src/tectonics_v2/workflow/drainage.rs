//! Phase R1 — low-res steepest-descent drainage utility.
//!
//! For each continental cell `(i, j)`, walks 4-NESW steepest descent
//! up to `max_distance` steps and records its drainage destination:
//!
//! - **Oceanic reach** — descent crossed `S̃ ≤ sea_level_reference`:
//!   `target_idx = oceanic neighbour`, `path_length = N` where `N` is
//!   the number of moves taken.
//! - **Pit** — start cell or some descent step has no strictly-lower
//!   neighbour: `target_idx = current cell`, `path_length = steps so
//!   far`. Pits anchored on the start cell record `path_length = 0`.
//! - **Cycle** — the next steepest neighbour is already in the visited
//!   path (rare on smooth fields, possible on flat plateaus): aborts
//!   at the pre-cycle cell. `path_length = steps taken so far`.
//! - **Max-distance** — `path_length = max_distance`, `target_idx =
//!   last reached cell`.
//!
//! ## Conventions
//!
//! - **Periodic boundaries (toroidal)** — same convention as the Stokes
//!   solver and `low_res_erosion`. Neighbours of `(i, j)` are
//!   `(i, j-1), (i+1, j), (i, j+1), (i-1, j)` with wrap.
//! - **Tie-break NESW** — when two neighbours have equal `S̃`, the
//!   one earlier in NESW order wins. Deterministic; no randomness.
//! - **Flat index `j * nx + i`** — row-major, identical to
//!   [`Field2D::set`] and `state.s_field`. Same convention as
//!   [`super::low_res_erosion`]. (The Step 12 brief proposed
//!   `i * ny + j` but every consumer in `tectonics_v2` is row-major
//!   `j * nx + i`; following codebase convention avoids a layout
//!   mismatch in Phase R2's macro-redistribution lookups.)
//! - **Oceanic cells stay quiet** — cells with `S̃ ≤ sea_level_reference`
//!   record `target_idx = self`, `path_length = 0`. Phase R2 only
//!   reads drainage destinations for cells that get eroded, which by
//!   construction are continental, so the oceanic-self anchor never
//!   matters semantically — but having it spelled out lets the test
//!   suite assert `target_idx[oceanic] == oceanic`.
//!
//! ## Complexity
//!
//! `O(nx · ny · max_distance²)` worst case — `Vec::contains` for cycle
//! detection. At Phase R5 defaults (`max_distance ≈ 10`, `nx · ny =
//! 64²`) this is ~410 k comparisons, well under one millisecond per
//! call. A `HashSet`-based path-visited check would lower the constant
//! but is over-engineering at this stage.

use crate::tectonics_v2::field::Field2D;

/// Result of a drainage pass — per-cell destination index and the
/// number of NESW moves taken to reach it.
///
/// Both buffers are flat row-major (`j * nx + i`) of length `nx · ny`.
#[derive(Clone, Debug)]
pub struct DrainageMap {
    /// `target_idx[k]` is the flat index of the drainage destination
    /// for cell `k`. For oceanic cells (`S̃[k] ≤ sea_level_reference`)
    /// and pits anchored at the start cell, `target_idx[k] == k`.
    pub target_idx: Vec<usize>,
    /// `path_length[k]` is the number of NESW moves taken from cell
    /// `k` to its destination. `0` for oceanic cells and start-anchored
    /// pits; positive for everything else; bounded above by
    /// `max_distance` (so `u8` is large enough for any practical
    /// `max_distance ≤ 255`).
    pub path_length: Vec<u8>,
}

/// Compute the drainage destination of every cell in `s_field`.
///
/// See module docstring for the algorithm and conventions. `max_distance
/// == 0` is a valid no-op (every cell drains to itself with length 0).
pub fn compute_drainage_targets(
    s_field: &Field2D,
    sea_level_reference: f64,
    max_distance: usize,
) -> DrainageMap {
    let nx = s_field.nx();
    let ny = s_field.ny();
    let n = nx * ny;
    let data = s_field.data();

    // Periodic neighbour LUTs (NESW).
    let prev_x: Vec<usize> = (0..nx).map(|i| (i + nx - 1) % nx).collect();
    let next_x: Vec<usize> = (0..nx).map(|i| (i + 1) % nx).collect();
    let prev_y: Vec<usize> = (0..ny).map(|j| (j + ny - 1) % ny).collect();
    let next_y: Vec<usize> = (0..ny).map(|j| (j + 1) % ny).collect();

    // Default: every cell drains to itself (length 0). Continental
    // cells overwrite below; oceanic cells stay quiet.
    let mut target_idx: Vec<usize> = (0..n).collect();
    let mut path_length = vec![0_u8; n];

    if max_distance == 0 {
        return DrainageMap { target_idx, path_length };
    }

    // Single reused path buffer — reset per start cell instead of
    // reallocated. Capacity bounded by `max_distance + 1`.
    let mut path: Vec<usize> = Vec::with_capacity(max_distance + 1);

    for j in 0..ny {
        for i in 0..nx {
            let start = j * nx + i;
            if data[start] <= sea_level_reference {
                continue; // oceanic — keep self / 0 anchor
            }

            let (target, length) = walk_steepest_descent(
                start,
                data,
                nx,
                sea_level_reference,
                max_distance,
                &prev_x,
                &next_x,
                &prev_y,
                &next_y,
                &mut path,
            );
            target_idx[start] = target;
            path_length[start] = length;
        }
    }

    DrainageMap { target_idx, path_length }
}

#[allow(clippy::too_many_arguments)]
fn walk_steepest_descent(
    start: usize,
    data: &[f64],
    nx: usize,
    sea_level_reference: f64,
    max_distance: usize,
    prev_x: &[usize],
    next_x: &[usize],
    prev_y: &[usize],
    next_y: &[usize],
    path: &mut Vec<usize>,
) -> (usize, u8) {
    path.clear();
    path.push(start);

    let mut current = start;
    for step in 0..max_distance {
        let ci = current % nx;
        let cj = current / nx;

        // 4 NESW neighbours in N→E→S→W order. The order matters:
        // ties (equal `S̃`) fall through to the first-found neighbour
        // since `if n_h < best_h` is strict.
        let neighbours = [
            prev_y[cj] * nx + ci, // N
            cj * nx + next_x[ci], // E
            next_y[cj] * nx + ci, // S
            cj * nx + prev_x[ci], // W
        ];
        let current_h = data[current];

        let mut best_lin: Option<usize> = None;
        let mut best_h = current_h;
        for &n_lin in &neighbours {
            let n_h = data[n_lin];
            if n_h < best_h {
                best_h = n_h;
                best_lin = Some(n_lin);
            }
        }

        let Some(next) = best_lin else {
            // Pit — no strictly-lower neighbour. `step` counts the
            // NESW moves already taken (0 if pit at start).
            return (current, step as u8);
        };

        // Crossed sea level → oceanic neighbour is the basin sink.
        if data[next] <= sea_level_reference {
            return (next, (step + 1) as u8);
        }

        // Cycle — `next` is somewhere in our visited path. Abort at
        // the current cell so the caller's redistribution never feeds
        // a loop. Cycle detection is O(path.len()) per step; at
        // realistic `max_distance` this is irrelevant. Rare in
        // practice (smooth Phase A fields) — happens only on flat
        // plateaus where multiple NESW siblings tie.
        if path.contains(&next) {
            return (current, step as u8);
        }

        current = next;
        path.push(current);
    }

    // Walked the full budget without reaching ocean or a pit.
    (current, max_distance as u8)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a symmetric ridge: oceanic rows at both top and bottom of
    /// the grid, continental peak at `j ≈ ny/2`. Each row is uniform
    /// so the NESW tie-break is exercised at the peak (N and S
    /// neighbours have equal S̃; N wins by NESW priority).
    ///
    /// Why a ridge rather than a simple top→bottom slope: the Stokes
    /// solver uses periodic boundaries (toroidal grid). On a linear
    /// `j=0` top → `j=ny-1` bottom slope, cell `(_, 0)` sees
    /// `(_, ny-1)` as its N neighbour — which is oceanic — and drains
    /// in *one* step via wrap. The ridge layout deliberately makes
    /// the wrap N-S boundary both-oceanic so it cannot be used as a
    /// drainage shortcut for any continental cell.
    fn ridge_field(nx: usize, ny: usize, sea_level: f64) -> Field2D {
        let mut s = Field2D::new(nx, ny);
        let peak_j = ny / 2;
        let h_peak = sea_level + 0.5;
        let h_edge = sea_level - 0.3;
        let max_dist = peak_j.max(ny - peak_j - 1) as f64;
        for j in 0..ny {
            let dist = ((j as i64) - (peak_j as i64)).unsigned_abs() as f64;
            let t = dist / max_dist; // 0 at peak, 1 at the farthest edge
            let val = h_peak + t * (h_edge - h_peak);
            for i in 0..nx {
                s.set(i, j, val);
            }
        }
        s
    }

    #[test]
    fn drainage_simple_slope_reaches_ocean() {
        // 8×8 ridge. With sea_level=0.5, h_peak=1.0, h_edge=0.2,
        // peak_j=4, max_dist=4: val(j) = 1.0 - 0.8 · |j-4|/4.
        //
        //   j |  0    1    2    3    4    5    6    7
        // val | 0.2  0.4  0.6  0.8  1.0  0.8  0.6  0.4
        //
        // Oceanic (val ≤ 0.5): rows {0, 1, 7}.
        // Continental (val > 0.5): rows {2, 3, 4, 5, 6}.
        let nx = 8;
        let ny = 8;
        let sea_level = 0.5;
        let s = ridge_field(nx, ny, sea_level);
        let map = compute_drainage_targets(&s, sea_level, 20);

        // Peak (j=4): N and S neighbours both at 0.8, NESW tie-break
        // picks N → descent through rows 3, 2, 1 (oceanic). length=3.
        for i in 0..nx {
            let lin = 4 * nx + i;
            assert_eq!(
                map.target_idx[lin], nx + i,
                "peak ({i}, 4) should drain N to ({i}, 1)"
            );
            assert_eq!(map.path_length[lin], 3);
        }

        // North slope (j=3): N=row 2 (continental, 0.6), then row 1
        // (oceanic). length=2, target=(_, 1).
        for i in 0..nx {
            let lin = 3 * nx + i;
            assert_eq!(map.target_idx[lin], nx + i);
            assert_eq!(map.path_length[lin], 2);
        }
        // North coast (j=2): N=row 1 oceanic, length=1.
        for i in 0..nx {
            let lin = 2 * nx + i;
            assert_eq!(map.target_idx[lin], nx + i);
            assert_eq!(map.path_length[lin], 1);
        }

        // South slope (j=5): N=row 4 peak (higher), so descent goes
        // S → row 6 (0.6), then row 7 (oceanic). length=2.
        for i in 0..nx {
            let lin = 5 * nx + i;
            assert_eq!(map.target_idx[lin], 7 * nx + i);
            assert_eq!(map.path_length[lin], 2);
        }
        // South coast (j=6): S=row 7 oceanic, length=1.
        for i in 0..nx {
            let lin = 6 * nx + i;
            assert_eq!(map.target_idx[lin], 7 * nx + i);
            assert_eq!(map.path_length[lin], 1);
        }

        // Oceanic rows {0, 1, 7} target self with length 0.
        for &j in &[0_usize, 1, 7] {
            for i in 0..nx {
                let lin = j * nx + i;
                assert_eq!(map.target_idx[lin], lin);
                assert_eq!(map.path_length[lin], 0);
            }
        }
    }

    #[test]
    fn drainage_pit_stays_local() {
        // Plateau at S̃=1.0 with a single pit at (3, 3) sunk to S̃=0.6.
        // All cells stay continental (sea_level=0.4 < 0.6 < 1.0); the
        // pit is surrounded by cells higher than itself.
        let nx = 8;
        let ny = 8;
        let sea_level = 0.4;
        let mut s = Field2D::filled(nx, ny, 1.0);
        let pit = (3, 3);
        s.set(pit.0, pit.1, 0.6);
        let pit_lin = pit.1 * nx + pit.0;

        let map = compute_drainage_targets(&s, sea_level, 10);

        // The pit drains to itself with length 0.
        assert_eq!(map.target_idx[pit_lin], pit_lin);
        assert_eq!(map.path_length[pit_lin], 0);

        // The 4 NESW neighbours of the pit (all at S̃=1.0) drain to
        // the pit in one step.
        let pit_neighbours = [
            (3, 2), // N
            (4, 3), // E
            (3, 4), // S
            (2, 3), // W
        ];
        for (i, j) in pit_neighbours {
            let lin = j * nx + i;
            assert_eq!(map.target_idx[lin], pit_lin, "neighbour ({i}, {j})");
            assert_eq!(map.path_length[lin], 1);
        }

        // Cells far from the pit are on a flat plateau (all neighbours
        // equal at S̃=1.0) → pits at start cell. `target = self`,
        // `length = 0`. Pick (0, 0) — far enough that the pit's
        // descent gradient never reaches it.
        let far = 0_usize;
        assert_eq!(map.target_idx[far], far);
        assert_eq!(map.path_length[far], 0);
    }

    #[test]
    fn drainage_long_distance_capped() {
        // 16×16 ridge, sea_level=0.5, h_peak=1.0, h_edge=0.2, peak=8,
        // max_dist=8: val(j) = 1.0 - 0.8 · |j-8|/8.
        //
        //   j  |  0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
        // val  | 0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1.0  0.9  0.8  0.7  0.6  0.5  0.4  0.3
        //
        // Oceanic (val ≤ 0.5): rows {0, 1, 2, 3, 13, 14, 15}.
        // Continental: rows {4..=12}.
        //
        // From the peak (j=8), descent N (tie-break) takes 5 hops to
        // hit oceanic row 3. With max_distance=3 the path saturates
        // mid-slope at row 5.
        let nx = 16;
        let ny = 16;
        let sea_level = 0.5;
        let s = ridge_field(nx, ny, sea_level);
        let cap = 3_usize;
        let map = compute_drainage_targets(&s, sea_level, cap);

        for i in 0..nx {
            let lin = 8 * nx + i; // peak row
            assert_eq!(
                map.path_length[lin], cap as u8,
                "peak row should saturate at cap={cap}"
            );
            assert_eq!(map.target_idx[lin], 5 * nx + i);
        }

        // No `path_length` ever exceeds the cap.
        for &p in &map.path_length {
            assert!(p as usize <= cap, "path_length {} exceeds cap {}", p, cap);
        }
    }

    #[test]
    fn drainage_deterministic() {
        // Re-running on the same field must produce bit-identical
        // outputs (same Vecs). Includes a non-trivial field with
        // ties so the NESW tie-break path is exercised.
        let nx = 12;
        let ny = 10;
        let sea_level = 0.5;
        let s = ridge_field(nx, ny, sea_level);
        let m1 = compute_drainage_targets(&s, sea_level, 10);
        let m2 = compute_drainage_targets(&s, sea_level, 10);
        assert_eq!(m1.target_idx, m2.target_idx);
        assert_eq!(m1.path_length, m2.path_length);
    }

    #[test]
    fn drainage_oceanic_cells_anchor_to_self() {
        // Mixed field. Verify every oceanic cell has
        // target_idx == self, length == 0 (the documented contract).
        let nx = 8;
        let ny = 8;
        let sea_level = 0.5;
        let s = ridge_field(nx, ny, sea_level);
        let map = compute_drainage_targets(&s, sea_level, 10);
        let data = s.data();
        for k in 0..nx * ny {
            if data[k] <= sea_level {
                assert_eq!(map.target_idx[k], k, "oceanic cell {k}");
                assert_eq!(map.path_length[k], 0, "oceanic cell {k}");
            }
        }
    }

    #[test]
    fn drainage_max_distance_zero_is_self_everywhere() {
        let nx = 4;
        let ny = 4;
        let s = ridge_field(nx, ny, 0.5);
        let map = compute_drainage_targets(&s, 0.5, 0);
        for k in 0..nx * ny {
            assert_eq!(map.target_idx[k], k);
            assert_eq!(map.path_length[k], 0);
        }
    }

    /// Phase R1 diagnostic — distribution of `path_length` over a
    /// realistic 32² active-medley-like Voronoï tessellation using
    /// `InitMode::RadialProfile` (Step 13 shape: peak at plate centre,
    /// descent to plate edge). Prints a histogram to stdout.
    ///
    /// Run with:
    /// ```bash
    /// cargo test --release -p ymir-core --lib \
    ///   tectonics_v2::workflow::drainage::tests::drainage_path_length_distribution_diagnostic \
    ///   -- --ignored --nocapture
    /// ```
    ///
    /// Expected (per the Step 12 refactor brief): coastal continental
    /// cells path_length 1–2 (immediate descent across coast), interior
    /// cells 5–10 (full radial slope to plate edge), pits (rare on
    /// smooth fields) path_length 0. If all cells land at length 1,
    /// the algorithm is degenerate (likely a periodic wrap shortcut in
    /// the test field — but a real Voronoï tessellation has
    /// continental patches embedded in oceanic surroundings so this
    /// should not occur).
    #[test]
    #[ignore]
    fn drainage_path_length_distribution_diagnostic() {
        use crate::tectonics_v2::init::{
            init_s_field, InitContext, InitMode, PlateInitData, ProfileShape,
        };
        use crate::tectonics_v2::voronoi::{generate_voronoi, VoronoiConfig};

        let nx = 32;
        let ny = 32;
        let seed = 42;
        let cfg = VoronoiConfig { num_plates: 8, continental_ratio: 0.3 };
        let plates = generate_voronoi(nx, ny, &cfg, seed);

        let plate_data = PlateInitData {
            plate_id: &plates.plate_id,
            plate_type: &plates.plate_type,
            seed_coords: Some(&plates.seed_coords),
        };
        let init_ctx = InitContext { nx, ny, seed, amplitude: 0.2, plate_data: Some(plate_data) };
        let init_mode = InitMode::RadialProfile {
            continental_value: 1.0,
            oceanic_value: 0.2,
            profile_shape: ProfileShape::Smoothstep,
        };
        let s = init_s_field(init_mode, &init_ctx);

        let sea_level = 0.5;
        let max_distance = 10;
        let map = compute_drainage_targets(&s, sea_level, max_distance);

        let data = s.data();
        let mut hist = vec![0_usize; max_distance + 1];
        let mut continental_count = 0_usize;
        for k in 0..nx * ny {
            if data[k] > sea_level {
                continental_count += 1;
                hist[map.path_length[k] as usize] += 1;
            }
        }

        println!(
            "\n=== R1 drainage stats — 32² active-medley-like, max_distance={}, sea_level={} ===",
            max_distance, sea_level
        );
        println!(
            "Total cells: {}  continental: {} ({:.1} %)",
            nx * ny,
            continental_count,
            100.0 * continental_count as f64 / (nx * ny) as f64
        );
        println!("Path length histogram (continental cells):");
        for (len, count) in hist.iter().enumerate() {
            let pct = 100.0 * (*count) as f64 / continental_count.max(1) as f64;
            let bar = "#".repeat((pct / 2.0).round() as usize);
            println!("  len {:>2}: {:>5} ({:>5.1} %) {}", len, count, pct, bar);
        }
        let mean = hist
            .iter()
            .enumerate()
            .map(|(len, count)| (len * count) as f64)
            .sum::<f64>()
            / continental_count.max(1) as f64;
        let max_len = hist.iter().rposition(|&c| c > 0).unwrap_or(0);
        println!("Mean path_length over continental cells: {:.2}", mean);
        println!("Max path_length: {}", max_len);

        assert!(continental_count > 0, "expected some continental cells");
        assert!(
            mean > 0.5,
            "drainage should produce non-trivial paths (mean > 0.5), got {}",
            mean
        );
    }

    /// Phase R1 diagnostic — same scenario as
    /// [`drainage_path_length_distribution_diagnostic`] but on a 64²
    /// grid with 4 large continental plates (`continental_ratio=0.5`,
    /// `Pow` profile with exponent 1.5 — sharper than Smoothstep so
    /// the slope is non-zero everywhere on the continent, not just at
    /// the rim). This is the layout where long-distance drainage
    /// (path_length 5–10) should appear.
    ///
    /// ```bash
    /// cargo test --release -p ymir-core --lib \
    ///   tectonics_v2::workflow::drainage::tests::drainage_long_distance_diagnostic \
    ///   -- --ignored --nocapture
    /// ```
    #[test]
    #[ignore]
    fn drainage_long_distance_diagnostic() {
        use crate::tectonics_v2::init::{
            init_s_field, InitContext, InitMode, PlateInitData, ProfileShape,
        };
        use crate::tectonics_v2::voronoi::{generate_voronoi, VoronoiConfig};

        let nx = 64;
        let ny = 64;
        let seed = 42;
        let cfg = VoronoiConfig { num_plates: 4, continental_ratio: 0.5 };
        let plates = generate_voronoi(nx, ny, &cfg, seed);

        let plate_data = PlateInitData {
            plate_id: &plates.plate_id,
            plate_type: &plates.plate_type,
            seed_coords: Some(&plates.seed_coords),
        };
        let init_ctx = InitContext { nx, ny, seed, amplitude: 0.2, plate_data: Some(plate_data) };
        // Pow with exponent 1.5 produces a non-flat continental
        // interior — every continental cell has a strictly-positive
        // outward slope, so steepest-descent runs farther than on a
        // Smoothstep plateau.
        let init_mode = InitMode::RadialProfile {
            continental_value: 1.0,
            oceanic_value: 0.2,
            profile_shape: ProfileShape::Pow { exponent: 1.5 },
        };
        let s = init_s_field(init_mode, &init_ctx);

        let sea_level = 0.5;
        let max_distance = 15;
        let map = compute_drainage_targets(&s, sea_level, max_distance);

        let data = s.data();
        let mut hist = vec![0_usize; max_distance + 1];
        let mut continental_count = 0_usize;
        for k in 0..nx * ny {
            if data[k] > sea_level {
                continental_count += 1;
                hist[map.path_length[k] as usize] += 1;
            }
        }

        println!(
            "\n=== R1 drainage stats — 64² large continents (Pow exp=1.5), max_distance={} ===",
            max_distance
        );
        println!(
            "Total cells: {}  continental: {} ({:.1} %)",
            nx * ny,
            continental_count,
            100.0 * continental_count as f64 / (nx * ny) as f64
        );
        println!("Path length histogram (continental cells):");
        for (len, count) in hist.iter().enumerate() {
            let pct = 100.0 * (*count) as f64 / continental_count.max(1) as f64;
            let bar = "#".repeat((pct / 2.0).round() as usize);
            println!("  len {:>2}: {:>5} ({:>5.1} %) {}", len, count, pct, bar);
        }
        let mean = hist
            .iter()
            .enumerate()
            .map(|(len, count)| (len * count) as f64)
            .sum::<f64>()
            / continental_count.max(1) as f64;
        let max_len = hist.iter().rposition(|&c| c > 0).unwrap_or(0);
        println!("Mean path_length over continental cells: {:.2}", mean);
        println!("Max path_length: {}", max_len);
    }
}
