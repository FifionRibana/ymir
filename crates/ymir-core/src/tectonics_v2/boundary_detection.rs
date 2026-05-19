//! Dynamic boundary-flag detection (Step 6).
//!
//! Every time step the `boundary_flag` field is recomputed from the
//! current velocity field. Flags change naturally as plates press
//! together, pull apart, or slide past one another.
//!
//! # Classification rule
//!
//! For each cell `(i, j)`, compute `div(v)_cell` using the Step 5
//! staggered finite-difference (no face reconstruction). Then:
//!
//! - `div > +threshold` → `Rift`. Source lookup (`Q_spread` vs
//!   `Q_rift_v`) is resolved downstream via `plate_type`.
//! - `div < -threshold` → convergent; classify by comparing `cell`'s
//!   plate_type to the argmax-neighbour (the neighbour contributing
//!   the most convergent motion into this cell):
//!   - `(Oceanic, Continental)` → `Subduction` on this cell.
//!   - `(Continental, Oceanic)` → `None` (the continental side
//!     does not drain; the arc it receives comes from the
//!     neighbour-lookup pattern installed in Step 5's
//!     `compute_source_sink_terms`).
//!   - `(Oceanic, Oceanic)` → `OceanicSubduction`.
//!   - `(Continental, Continental)` → `ContinentalCollision`.
//! - `|div| ≤ threshold` → `None`.
//!
//! # Threshold value
//!
//! `DetectionConfig::default().threshold = 1e-4`. This is a
//! **numerical** threshold, not derived from physical scales: above
//! the machine-noise floor of `div(v)_cell` in a floor-dominated
//! Ar = 0.1 solve, below the divergence signatures that should
//! activate boundary mechanisms in a developed flow.

use super::boundaries::boundary_flag::BoundaryFlag;
use super::boundaries::plate_type::{PlateType, PlateTypeField};
use super::field::PeriodicIndex;
use super::voronoi::PlateIdField;

#[derive(Clone, Copy, Debug)]
pub struct DetectionConfig {
    pub threshold: f64,
}

impl Default for DetectionConfig {
    fn default() -> Self {
        Self { threshold: 1.0e-4 }
    }
}

/// Classify each cell's `boundary_flag` from the current velocity
/// field. The divergence is computed on-the-fly from `(vx, vy)` on
/// the MAC staggered grid. The `plate_id` is accepted but not
/// consulted directly — plate-type lookup is sufficient for
/// classification. It is kept in the signature so callers can
/// extend the rule later without changing the API.
pub fn detect_boundaries(
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    idx_x: &PeriodicIndex,
    idx_y: &PeriodicIndex,
    vx: &[f64],
    vy: &[f64],
    plate_type: &PlateTypeField,
    _plate_id: &PlateIdField,
    config: &DetectionConfig,
    out: &mut super::boundaries::boundary_flag::BoundaryFlagField,
) {
    debug_assert_eq!(plate_type.nx(), nx);
    debug_assert_eq!(plate_type.ny(), ny);
    debug_assert_eq!(out.nx(), nx);
    debug_assert_eq!(out.ny(), ny);
    debug_assert_eq!(vx.len(), nx * ny);
    debug_assert_eq!(vy.len(), nx * ny);

    let inv_dx = 1.0 / dx;
    let inv_dy = 1.0 / dy;
    let lin = |i: usize, j: usize| j * nx + i;
    let thr = config.threshold;

    for j in 0..ny {
        let jp = idx_y.next(j);
        let jm = idx_y.prev(j);
        for i in 0..nx {
            let ip = idx_x.next(i);
            let im = idx_x.prev(i);

            let div_v = (vx[lin(ip, j)] - vx[lin(i, j)]) * inv_dx
                + (vy[lin(i, jp)] - vy[lin(i, j)]) * inv_dy;

            let flag = if div_v > thr {
                BoundaryFlag::Rift
            } else if div_v < -thr {
                // Convergent: find the neighbour contributing the
                // most convergent motion into this cell.
                //   From E (ip,j): convergent contribution is
                //     `max(0, -vx[ip, j])` (vx positive = outflow at
                //     the right face, so -vx = inflow from east).
                //   From W (im,j): `max(0, +vx[i, j])`.
                //   From N (i,jp): `max(0, -vy[i, jp])`.
                //   From S (i,jm): `max(0, +vy[i, j])`.
                let c_e = (-vx[lin(ip, j)]).max(0.0);
                let c_w = vx[lin(i, j)].max(0.0);
                let c_n = (-vy[lin(i, jp)]).max(0.0);
                let c_s = vy[lin(i, j)].max(0.0);
                let mut best_c = c_e;
                let mut best_idx: u8 = 0;
                if c_w > best_c { best_c = c_w; best_idx = 1; }
                if c_n > best_c { best_c = c_n; best_idx = 2; }
                if c_s > best_c { best_idx = 3; }
                let _ = best_c;

                let neighbor = match best_idx {
                    0 => plate_type.get(ip, j),
                    1 => plate_type.get(im, j),
                    2 => plate_type.get(i, jp),
                    _ => plate_type.get(i, jm),
                };
                let cell_t = plate_type.get(i, j);
                classify_convergent(cell_t, neighbor)
            } else {
                BoundaryFlag::None
            };

            out.set(i, j, flag);
        }
    }
}

/// Decide the flag for a cell given its type and the type of the
/// neighbour contributing the most convergent motion. See module
/// doc for the rule.
#[inline]
fn classify_convergent(cell: PlateType, neighbor: PlateType) -> BoundaryFlag {
    match (cell, neighbor) {
        (PlateType::Oceanic, PlateType::Continental) => BoundaryFlag::Subduction,
        (PlateType::Continental, PlateType::Oceanic) => BoundaryFlag::None,
        (PlateType::Oceanic, PlateType::Oceanic) => BoundaryFlag::OceanicSubduction,
        (PlateType::Continental, PlateType::Continental) => BoundaryFlag::ContinentalCollision,
    }
}

#[cfg(test)]
mod tests {
    use super::super::boundaries::boundary_flag::BoundaryFlagField;
    use super::super::boundaries::plate_type::{PlateType, PlateTypeField};
    use super::super::voronoi::PlateIdField;
    use super::*;
    use std::f64::consts::PI;

    fn run_detect(
        nx: usize,
        ny: usize,
        vx: &[f64],
        vy: &[f64],
        plate_type: &PlateTypeField,
        threshold: f64,
    ) -> BoundaryFlagField {
        let dx = 1.0 / nx as f64;
        let dy = 1.0 / ny as f64;
        let idx_x = PeriodicIndex::new(nx);
        let idx_y = PeriodicIndex::new(ny);
        let pid = PlateIdField::new(nx, ny);
        let mut out = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
        let cfg = DetectionConfig { threshold };
        detect_boundaries(
            nx, ny, dx, dy, &idx_x, &idx_y, vx, vy, plate_type, &pid, &cfg, &mut out,
        );
        out
    }

    #[test]
    fn rift_fires_where_divergence_positive() {
        // v = (sin(2πx), 0) on 32x8. div = 2π cos(2πx).
        // Rift where cos(2πx) > 0 and large enough.
        let nx = 32;
        let ny = 8;
        let dx = 1.0 / nx as f64;
        let mut vx = vec![0.0; nx * ny];
        let vy = vec![0.0; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let x_face = i as f64 * dx;
                vx[j * nx + i] = (2.0 * PI * x_face).sin();
            }
        }
        let plate_type = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
        let out = run_detect(nx, ny, &vx, &vy, &plate_type, 1e-4);
        // Row j=0: look at cells where x_c = (i+0.5)/nx. The
        // analytical div at cell centre is 2π cos(2π · x_c). So
        // cells with cos(2π · x_c) > 0 AND large enough magnitude →
        // Rift; cells with cos < 0 → convergent → OceanicSubduction
        // (since all oceanic).
        let mut rift_count = 0;
        let mut sub_count = 0;
        let mut none_count = 0;
        for i in 0..nx {
            match out.get(i, 0) {
                BoundaryFlag::Rift => rift_count += 1,
                BoundaryFlag::OceanicSubduction => sub_count += 1,
                BoundaryFlag::None => none_count += 1,
                _ => panic!("unexpected flag on oceanic-only domain"),
            }
        }
        assert!(rift_count > 0);
        assert!(sub_count > 0);
        assert!(rift_count + sub_count + none_count == nx);
    }

    #[test]
    fn classify_convergent_respects_plate_type_pair() {
        assert!(matches!(
            classify_convergent(PlateType::Oceanic, PlateType::Continental),
            BoundaryFlag::Subduction
        ));
        assert!(matches!(
            classify_convergent(PlateType::Continental, PlateType::Oceanic),
            BoundaryFlag::None
        ));
        assert!(matches!(
            classify_convergent(PlateType::Oceanic, PlateType::Oceanic),
            BoundaryFlag::OceanicSubduction
        ));
        assert!(matches!(
            classify_convergent(PlateType::Continental, PlateType::Continental),
            BoundaryFlag::ContinentalCollision
        ));
    }

    #[test]
    fn below_threshold_maps_to_none() {
        // Zero velocity → div = 0 → all cells None.
        let nx = 8;
        let ny = 8;
        let vx = vec![0.0; nx * ny];
        let vy = vec![0.0; nx * ny];
        let plate_type = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
        let out = run_detect(nx, ny, &vx, &vy, &plate_type, 1e-4);
        for j in 0..ny {
            for i in 0..nx {
                assert!(matches!(out.get(i, j), BoundaryFlag::None));
            }
        }
    }
}
