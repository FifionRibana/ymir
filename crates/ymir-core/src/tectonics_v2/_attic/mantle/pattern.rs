//! MAC-staggered mantle velocity pattern `v_pattern = curl(ψ ẑ)`.
//!
//! Given a **nodal** stream function `ψ` (stored as a `Field2D`
//! `nx × ny` with entry `[i, j]` representing `ψ(i·dx, j·dy)`,
//! see `stream_function.rs`), this module builds the two
//! face-centered velocity components on the standard MAC
//! staggered grid used throughout `tectonics_v2`:
//!
//! - `vx[i, j]` lives at `(i·dx, (j+½)·dy)` (left vertical face
//!   of cell `(i, j)`).
//! - `vy[i, j]` lives at `((i+½)·dx, j·dy)` (bottom horizontal
//!   face of cell `(i, j)`).
//!
//! Staggered curl formulas
//! ------------------------
//!
//! ```text
//!   vx[i, j] = (ψ[i, j+1] − ψ[i, j]) / dy          // ∂ψ/∂y at x = i·dx
//!   vy[i, j] = −(ψ[i+1, j] − ψ[i, j]) / dx         // −∂ψ/∂x at y = j·dy
//! ```
//!
//! Both indices wrap periodically. The formulas are the first
//! differences of nodal ψ across the corners bounding each face:
//! a vertical face at `x = i·dx` is bounded by nodes `(i, j)` and
//! `(i, j+1)`; a horizontal face at `y = j·dy` by `(i, j)` and
//! `(i+1, j)`.
//!
//! Exact discrete div-freeness
//! ---------------------------
//!
//! The staggered divergence at cell `(i, j)` is
//!
//! ```text
//!   div(v)[i,j] = (vx[i+1,j] − vx[i,j]) / dx + (vy[i,j+1] − vy[i,j]) / dy
//! ```
//!
//! Substituting the formulas above and expanding:
//!
//! ```text
//!   (ψ[i+1,j+1] − ψ[i+1,j] − ψ[i,j+1] + ψ[i,j]) / (dx · dy)
//! + (−ψ[i+1,j+1] + ψ[i,j+1] + ψ[i+1,j] − ψ[i,j]) / (dx · dy)
//! = 0
//! ```
//!
//! — the four corner terms cancel algebraically, independent of
//! grid resolution. `div(v) ≡ 0` to f64 precision at 64², 128²,
//! 256², …, which satisfies the Step 8 strict acceptance
//! `div_v_mantle_max < 10⁻¹⁰`.
//!
//! This is the reason the stream function is stored **nodally**
//! and not cell-centered: a cell-centered discretisation would
//! leave a residual `div ≈ O(dx²)` which at 256² is already
//! `~1.5e-5`, far above the acceptance.

use crate::tectonics_v2::field::{Field2D, PeriodicIndex};

/// MAC-staggered mantle velocity pattern.
///
/// Both components are `nx × ny` `Field2D`s storing face values.
#[derive(Clone)]
pub struct MantlePattern {
    pub v_mantle_x: Field2D,
    pub v_mantle_y: Field2D,
}

impl std::fmt::Debug for MantlePattern {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Field2D does not derive Debug (legacy convention in this
        // crate — cf. SlabState). Dump shape + peak magnitude.
        f.debug_struct("MantlePattern")
            .field("nx", &self.nx())
            .field("ny", &self.ny())
            .field("peak_magnitude", &self.peak_magnitude())
            .finish()
    }
}

impl MantlePattern {
    pub fn nx(&self) -> usize {
        self.v_mantle_x.nx()
    }
    pub fn ny(&self) -> usize {
        self.v_mantle_x.ny()
    }

    /// Peak `|v_mantle|` over the domain (sampled by taking the
    /// max of `|v_x|` and `|v_y|` cell-by-cell). Used by the
    /// report as a telemetry anchor `peak|v_mantle| = Mf · peak|v_pattern|`.
    pub fn peak_magnitude(&self) -> f64 {
        let mut peak = 0.0_f64;
        let nx = self.nx();
        let ny = self.ny();
        for j in 0..ny {
            for i in 0..nx {
                let vx = self.v_mantle_x.get(i, j);
                let vy = self.v_mantle_y.get(i, j);
                let mag = (vx * vx + vy * vy).sqrt();
                if mag > peak {
                    peak = mag;
                }
            }
        }
        peak
    }

    /// Recompute `v_mantle = curl(ψ ẑ)` in place from a new nodal
    /// `ψ`. Same algebra as [`build_mantle_pattern`] but writes
    /// into the existing buffers (no allocation). Used by Step 12
    /// R6 to rebuild the pattern each step under phase drift —
    /// the buffer footprint stays constant across the run.
    pub fn rebuild_from_psi(
        &mut self,
        psi_nodal: &Field2D,
        dx: f64,
        dy: f64,
        idx_x: &PeriodicIndex,
        idx_y: &PeriodicIndex,
    ) {
        debug_assert_eq!(psi_nodal.nx(), self.nx());
        debug_assert_eq!(psi_nodal.ny(), self.ny());
        let nx = self.nx();
        let ny = self.ny();
        let inv_dx = 1.0 / dx;
        let inv_dy = 1.0 / dy;
        for j in 0..ny {
            let jp = idx_y.next(j);
            for i in 0..nx {
                let ip = idx_x.next(i);
                let vx = (psi_nodal.get(i, jp) - psi_nodal.get(i, j)) * inv_dy;
                self.v_mantle_x.set(i, j, vx);
                let vy = -(psi_nodal.get(ip, j) - psi_nodal.get(i, j)) * inv_dx;
                self.v_mantle_y.set(i, j, vy);
            }
        }
    }
}

/// Build `v_mantle = curl(ψ ẑ)` from a nodal stream function on
/// a MAC staggered velocity grid.
///
/// Both `dx` and `dy` are passed in to match the harness's grid
/// spacing; on square domains they are equal, but keeping them
/// separate is harmless and matches the rest of the crate.
pub fn build_mantle_pattern(
    psi_nodal: &Field2D,
    dx: f64,
    dy: f64,
    idx_x: &PeriodicIndex,
    idx_y: &PeriodicIndex,
) -> MantlePattern {
    let nx = psi_nodal.nx();
    let ny = psi_nodal.ny();
    let mut v_mantle_x = Field2D::new(nx, ny);
    let mut v_mantle_y = Field2D::new(nx, ny);
    let inv_dx = 1.0 / dx;
    let inv_dy = 1.0 / dy;
    for j in 0..ny {
        let jp = idx_y.next(j);
        for i in 0..nx {
            let ip = idx_x.next(i);
            // vx[i, j] = (ψ[i, j+1] − ψ[i, j]) / dy
            let vx = (psi_nodal.get(i, jp) - psi_nodal.get(i, j)) * inv_dy;
            v_mantle_x.set(i, j, vx);
            // vy[i, j] = −(ψ[i+1, j] − ψ[i, j]) / dx
            let vy = -(psi_nodal.get(ip, j) - psi_nodal.get(i, j)) * inv_dx;
            v_mantle_y.set(i, j, vy);
        }
    }
    MantlePattern { v_mantle_x, v_mantle_y }
}

/// Compute the max absolute value of the staggered divergence of
/// the pattern. Returns `0` (to f64 precision) by construction;
/// the helper exists so tests and diagnostics can sample it
/// empirically.
pub fn pattern_div_max(
    pattern: &MantlePattern,
    dx: f64,
    dy: f64,
    idx_x: &PeriodicIndex,
    idx_y: &PeriodicIndex,
) -> f64 {
    let nx = pattern.nx();
    let ny = pattern.ny();
    let inv_dx = 1.0 / dx;
    let inv_dy = 1.0 / dy;
    let mut max = 0.0_f64;
    for j in 0..ny {
        let jp = idx_y.next(j);
        for i in 0..nx {
            let ip = idx_x.next(i);
            let dvx = (pattern.v_mantle_x.get(ip, j) - pattern.v_mantle_x.get(i, j)) * inv_dx;
            let dvy = (pattern.v_mantle_y.get(i, jp) - pattern.v_mantle_y.get(i, j)) * inv_dy;
            let div = dvx + dvy;
            if div.abs() > max {
                max = div.abs();
            }
        }
    }
    max
}

#[cfg(test)]
mod tests {
    use super::super::stream_function::{StreamFunctionConfig, generate_stream_function};
    use super::*;

    /// Core contract: curl of nodal ψ gives discrete-div-free
    /// velocity on the staggered grid at any resolution.
    #[test]
    fn curl_is_exactly_div_free() {
        for &n in &[16, 32, 64, 128] {
            let dx = 1.0 / n as f64;
            let idx_x = PeriodicIndex::new(n);
            let idx_y = PeriodicIndex::new(n);
            let psi =
                generate_stream_function(n, n, &StreamFunctionConfig { num_modes: 6, seed: 42 });
            let pattern = build_mantle_pattern(&psi, dx, dx, &idx_x, &idx_y);
            let div = pattern_div_max(&pattern, dx, dx, &idx_x, &idx_y);
            assert!(div < 1.0e-10, "div max at N={} is {} (expected < 1e-10)", n, div,);
        }
    }

    /// ψ ≡ 0 ⇒ v_mantle ≡ 0.
    #[test]
    fn zero_psi_gives_zero_pattern() {
        let n = 8;
        let dx = 1.0 / n as f64;
        let idx_x = PeriodicIndex::new(n);
        let idx_y = PeriodicIndex::new(n);
        let psi = Field2D::new(n, n);
        let pattern = build_mantle_pattern(&psi, dx, dx, &idx_x, &idx_y);
        for &v in pattern.v_mantle_x.data().iter().chain(pattern.v_mantle_y.data().iter()) {
            assert_eq!(v, 0.0);
        }
    }

    /// Linear ψ = x ⇒ vx = 0, vy = -1 (uniform). On a periodic
    /// domain "ψ = x" can't exist literally — but ψ = sin(2πx)
    /// at nodes gives predictable behaviour we can probe.
    #[test]
    fn pattern_is_nontrivial_when_psi_nontrivial() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let idx_x = PeriodicIndex::new(n);
        let idx_y = PeriodicIndex::new(n);
        let psi = generate_stream_function(n, n, &StreamFunctionConfig { num_modes: 4, seed: 7 });
        let pattern = build_mantle_pattern(&psi, dx, dx, &idx_x, &idx_y);
        let peak = pattern.peak_magnitude();
        assert!(peak > 0.01, "peak |v_mantle| = {}, expected nontrivial", peak);
    }

    /// Scaling ψ by α scales v_mantle by α (curl is linear).
    #[test]
    fn pattern_is_linear_in_psi() {
        let n = 16;
        let dx = 1.0 / n as f64;
        let idx_x = PeriodicIndex::new(n);
        let idx_y = PeriodicIndex::new(n);
        let mut psi_a =
            generate_stream_function(n, n, &StreamFunctionConfig { num_modes: 4, seed: 3 });
        let mut psi_b = psi_a.clone();
        for v in psi_b.data_mut().iter_mut() {
            *v *= 3.0;
        }
        let pat_a = build_mantle_pattern(&psi_a, dx, dx, &idx_x, &idx_y);
        let pat_b = build_mantle_pattern(&psi_b, dx, dx, &idx_x, &idx_y);
        for k in 0..n * n {
            let ax = pat_a.v_mantle_x.data()[k];
            let bx = pat_b.v_mantle_x.data()[k];
            assert!((bx - 3.0 * ax).abs() < 1e-14);
            let ay = pat_a.v_mantle_y.data()[k];
            let by = pat_b.v_mantle_y.data()[k];
            assert!((by - 3.0 * ay).abs() < 1e-14);
        }
        let _ = &mut psi_a;
    }
}
