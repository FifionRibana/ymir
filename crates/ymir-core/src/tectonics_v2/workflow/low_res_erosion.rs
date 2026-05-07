//! Phase A low-res parametric erosion (D2 of `step12_issue.md`).
//!
//! Per-cycle diffusive erosion algorithm applied to all continental
//! cells:
//!
//! ```text
//! For each cell i with S̃[i] > sea_level_reference:
//!     slope_i = max |S̃[i] - S̃[neighbor]| over the 4 NESW neighbors
//!     Δh = α · slope_i · (S̃[i] - sea_level_reference)
//!     S̃[i] -= Δh
//!     if β > 0:
//!         downslope_neighbor.S̃ += β · Δh
//! ```
//!
//! Conventions:
//!
//! - **Boundaries**: fully periodic (toroidal), matching the Stokes
//!   solver. The 4 NESW neighbors of `(i, j)` are
//!   `(i, j-1), (i+1, j), (i, j+1), (i-1, j)` with wrap.
//! - **Continental mask**: strict `S̃[i] > sea_level_reference`. The
//!   filter is applied at cell-selection time, so `Δh` is positive by
//!   construction (never need to clamp at zero).
//! - **`sea_level_reference`**: passed in by the caller. Phase 3
//!   orchestrator extracts it from
//!   [`crate::tectonics::isostasy::compute_isostasy`] —
//!   `IsostasyResult::sea_level_normalized` — for the adaptive
//!   threshold contract (Option 2 of Phase 0 finding E). Phase 2
//!   tests pass `0.5` for simplicity.
//! - **Redistribution `β > 0`**: deterministic NESW priority on ties.
//!   The downslope neighbor is the strictly-lowest of the 4; on
//!   equality the first wins in the order N, E, S, W. No randomness,
//!   no hash-dependent ordering.
//! - **Two-pass implementation**: Δh is computed read-only against
//!   the original `S̃` field (pass 1), then applied in-place with
//!   optional redistribution (pass 2). This avoids cross-cell
//!   contamination — every cell's Δh is computed from the same
//!   coherent snapshot.
//! - **Mass algebra**: total mass change is `-(1-β) · Σ Δh`. So
//!   `β = 1.0` conserves exactly (modulo IEEE-754 rounding —
//!   `O(ε · N · Δh̄)` ≈ 1e-16 at 32², well below the 1e-6 acceptance
//!   threshold) and `β = 0.0` decreases monotonically each cycle.

use super::PhaseAParams;
use crate::tectonics_v2::field::Field2D;

/// Per-call summary of the erosion pass. The orchestrator (Phase 3)
/// accumulates these into per-cycle metrics
/// (`erosion_volume_removed_per_cycle`).
///
/// `volume_removed` is the **gross** integrated `Σ Δh` over all
/// eroded cells, independent of `β`. Net mass change of the grid is
/// `-(1 - β) · volume_removed`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ErosionStats {
    pub volume_removed: f64,
    pub peak_delta_h: f64,
}

/// Apply one pass of low-res parametric erosion in-place on `s`.
///
/// See module docstring for the algorithm and conventions.
pub fn apply(
    s: &mut Field2D,
    params: &PhaseAParams,
    sea_level_reference: f64,
) -> ErosionStats {
    let nx = s.nx();
    let ny = s.ny();
    let alpha = params.alpha;
    let beta = params.beta;

    // Periodic neighbour LUTs (NESW). Allocations are O(nx + ny),
    // dominated by the O(nx · ny) computation pass — negligible.
    let prev_x: Vec<usize> = (0..nx).map(|i| (i + nx - 1) % nx).collect();
    let next_x: Vec<usize> = (0..nx).map(|i| (i + 1) % nx).collect();
    let prev_y: Vec<usize> = (0..ny).map(|j| (j + ny - 1) % ny).collect();
    let next_y: Vec<usize> = (0..ny).map(|j| (j + 1) % ny).collect();

    // Pass 1: compute Δh and downslope linear-index for every
    // continental cell, reading `s` only (no in-place mutation yet).
    let n_cells = nx * ny;
    let mut delta_h = vec![0.0_f64; n_cells];
    let mut downslope_lin = vec![0_usize; n_cells];
    {
        let s_data = s.data();
        for j in 0..ny {
            for i in 0..nx {
                let lin = j * nx + i;
                let s_i = s_data[lin];
                if s_i <= sea_level_reference {
                    continue;
                }
                // NESW order — N first, then E, S, W. Used both for
                // the slope magnitude and the downslope tie-breaker.
                let neighbors = [
                    (i, prev_y[j]),  // N
                    (next_x[i], j),  // E
                    (i, next_y[j]),  // S
                    (prev_x[i], j),  // W
                ];

                let mut max_slope = 0.0_f64;
                let mut best_idx = 0_usize;
                let mut best_h = s_data[neighbors[0].1 * nx + neighbors[0].0];

                for (k, &(ni, nj)) in neighbors.iter().enumerate() {
                    let s_n = s_data[nj * nx + ni];
                    let mag = (s_i - s_n).abs();
                    if mag > max_slope {
                        max_slope = mag;
                    }
                    // Strict-less: ties keep the NESW-priority winner
                    // (k=0 = N is initialised; k=1..4 only displaces
                    // when strictly lower).
                    if k > 0 && s_n < best_h {
                        best_h = s_n;
                        best_idx = k;
                    }
                }

                let dh = alpha * max_slope * (s_i - sea_level_reference);
                delta_h[lin] = dh;
                let down = neighbors[best_idx];
                downslope_lin[lin] = down.1 * nx + down.0;
            }
        }
    }

    // Pass 2: erode source + redistribute downslope (in-place).
    // Order-of-iteration does not affect total mass conservation
    // (the algebra holds independent of order); within-cycle
    // intra-cycle deposit-then-erode interactions are intentional —
    // they emerge from the diffusive reading.
    let mut volume_removed = 0.0_f64;
    let mut peak_dh = 0.0_f64;
    {
        let s_mut = s.data_mut();
        for lin in 0..n_cells {
            let dh = delta_h[lin];
            if dh <= 0.0 {
                continue;
            }
            s_mut[lin] -= dh;
            volume_removed += dh;
            if dh > peak_dh {
                peak_dh = dh;
            }
            if beta > 0.0 {
                let down_lin = downslope_lin[lin];
                s_mut[down_lin] += beta * dh;
            }
        }
    }

    ErosionStats { volume_removed, peak_delta_h: peak_dh }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flat_continental_field(nx: usize, ny: usize) -> Field2D {
        // Smooth periodic continental relief, all cells > 0.5
        // (continental). Range [0.6, 0.8], slope non-zero almost
        // everywhere — exercises the per-cell algorithm without any
        // coast for unit-test simplicity.
        let mut s = Field2D::new(nx, ny);
        let two_pi = 2.0 * std::f64::consts::PI;
        for j in 0..ny {
            for i in 0..nx {
                let sx = (i as f64 / nx as f64) * two_pi;
                let sy = (j as f64 / ny as f64) * two_pi;
                s.set(i, j, 0.7 + 0.05 * (sx.sin() + sy.cos()));
            }
        }
        s
    }

    #[test]
    fn empty_grid_returns_zero_stats() {
        let mut s = Field2D::filled(8, 8, 0.3);  // all oceanic
        let params = PhaseAParams { alpha: 0.05, beta: 0.0, ..Default::default() };
        let stats = apply(&mut s, &params, 0.5);
        assert_eq!(stats.volume_removed, 0.0);
        assert_eq!(stats.peak_delta_h, 0.0);
    }

    #[test]
    fn flat_continental_does_not_erode() {
        // Constant S̃ = 0.7 → slope = 0 everywhere → Δh = 0.
        let mut s = Field2D::filled(8, 8, 0.7);
        let mass_before: f64 = s.data().iter().sum();
        let params = PhaseAParams { alpha: 0.05, beta: 0.0, ..Default::default() };
        let stats = apply(&mut s, &params, 0.5);
        let mass_after: f64 = s.data().iter().sum();
        assert_eq!(stats.volume_removed, 0.0);
        assert_eq!(mass_before, mass_after);
    }

    #[test]
    fn relief_continental_erodes_with_zero_beta() {
        let mut s = flat_continental_field(16, 16);
        let mass_before: f64 = s.data().iter().sum();
        let params = PhaseAParams { alpha: 0.05, beta: 0.0, ..Default::default() };
        let stats = apply(&mut s, &params, 0.5);
        assert!(stats.volume_removed > 0.0);
        let mass_after: f64 = s.data().iter().sum();
        // β = 0.0 → mass strictly decreases by exactly volume_removed
        // (modulo IEEE-754 summation order).
        assert!(mass_before - mass_after > 0.0);
        assert!((mass_before - mass_after - stats.volume_removed).abs() < 1e-12);
    }

    #[test]
    fn full_redistribution_conserves_mass() {
        let mut s = flat_continental_field(16, 16);
        let mass_before: f64 = s.data().iter().sum();
        let params = PhaseAParams { alpha: 0.05, beta: 1.0, ..Default::default() };
        let _stats = apply(&mut s, &params, 0.5);
        let mass_after: f64 = s.data().iter().sum();
        assert!(
            (mass_after - mass_before).abs() < 1e-12,
            "β=1.0 must conserve mass exactly (modulo float rounding): \
             before={mass_before} after={mass_after} delta={}",
            mass_after - mass_before
        );
    }
}
