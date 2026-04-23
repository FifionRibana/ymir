//! Unified crust-geometry structure (Step 6).
//!
//! Step 5 carried `plate_types` and `flags` as separate `Arc`-wrapped
//! fields on `BoundaryConfig::Enabled`. Step 6 introduces dynamic
//! per-step detection of `boundary_flag` and a `plate_id` field for
//! Voronoi tessellations; the three fields become tightly coupled
//! (flag detection reads plate_type; stats group by plate_id; etc.).
//! Packing them into [`CrustGeometry`] is cheap and makes the
//! coupling explicit in the type system.
//!
//! The struct is the **initial** state. The run-local harness
//! clones the mutable part (`boundary_flag`) and updates it each
//! step when `geometry_kind == Voronoi`. The shared `Arc` wrapper
//! [`CrustGeometryShared`] lets multiple configs / tests reference
//! the same underlying tessellation without copying.

use std::sync::Arc;

use super::super::boundary_detection::DetectionConfig;
use super::super::voronoi::{generate_voronoi, PlateIdField, VoronoiConfig};
use super::boundary_flag::{BoundaryFlag, BoundaryFlagField};
use super::plate_type::{PlateType, PlateTypeField};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GeometryKind {
    /// Flags are fixed at init (Step 5 static layouts).
    Static,
    /// Flags are recomputed each step from `div(v)` (Step 6).
    Voronoi,
}

/// Initial crust geometry. The `boundary_flag` field is the initial
/// value; the harness owns a mutable run-local copy for dynamic
/// geometries.
#[derive(Clone, Debug)]
pub struct CrustGeometry {
    pub plate_id: PlateIdField,
    pub plate_type: PlateTypeField,
    pub boundary_flag: BoundaryFlagField,
    pub detection_config: DetectionConfig,
    pub geometry_kind: GeometryKind,
    pub layout_name: String,
}

pub type CrustGeometryShared = Arc<CrustGeometry>;

impl CrustGeometry {
    /// Build from a Step 5 static layout: plate_id filled with zero,
    /// boundary_flag copied from the layout, `geometry_kind = Static`.
    pub fn from_static(
        plate_type: PlateTypeField,
        boundary_flag: BoundaryFlagField,
        layout_name: impl Into<String>,
    ) -> Self {
        let nx = plate_type.nx();
        let ny = plate_type.ny();
        Self {
            plate_id: PlateIdField::filled(nx, ny, 0),
            plate_type,
            boundary_flag,
            detection_config: DetectionConfig::default(),
            geometry_kind: GeometryKind::Static,
            layout_name: layout_name.into(),
        }
    }

    /// Build from a Voronoi tessellation. The boundary_flag is
    /// initialized to `None` everywhere — it will be populated at
    /// run time by the first `detect_boundaries` call, after the
    /// Stokes solve at step 1 produces a non-trivial velocity
    /// field. `geometry_kind = Voronoi`.
    pub fn from_voronoi(nx: usize, ny: usize, config: &VoronoiConfig, seed: u64) -> Self {
        let plates = generate_voronoi(nx, ny, config, seed);
        Self {
            plate_id: plates.plate_id,
            plate_type: plates.plate_type,
            boundary_flag: BoundaryFlagField::filled(nx, ny, BoundaryFlag::None),
            detection_config: DetectionConfig::default(),
            geometry_kind: GeometryKind::Voronoi,
            layout_name: format!(
                "voronoi_seed{}_n{}",
                seed, config.num_plates,
            ),
        }
    }

    pub fn nx(&self) -> usize { self.plate_type.nx() }
    pub fn ny(&self) -> usize { self.plate_type.ny() }

    pub fn is_dynamic(&self) -> bool {
        matches!(self.geometry_kind, GeometryKind::Voronoi)
    }

    /// Count cells whose plate_type is Continental. Used by the
    /// physics report's type-distribution metric.
    pub fn continental_cell_count(&self) -> usize {
        self.plate_type
            .data()
            .iter()
            .filter(|&&t| matches!(t, PlateType::Continental))
            .count()
    }

    /// Count distinct plate ids observed in the grid (typically
    /// equal to `num_plates` when all plates host ≥ 1 cell).
    pub fn distinct_plate_count(&self) -> usize {
        let mut seen = std::collections::HashSet::new();
        for &id in self.plate_id.data() {
            seen.insert(id);
        }
        seen.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn static_geometry_has_zero_plate_id() {
        let pt = PlateTypeField::filled(4, 4, PlateType::Continental);
        let bf = BoundaryFlagField::filled(4, 4, BoundaryFlag::None);
        let g = CrustGeometry::from_static(pt, bf, "test");
        for j in 0..4 {
            for i in 0..4 {
                assert_eq!(g.plate_id.get(i, j), 0);
            }
        }
        assert!(!g.is_dynamic());
    }

    #[test]
    fn voronoi_geometry_is_dynamic() {
        let cfg = VoronoiConfig::default();
        let g = CrustGeometry::from_voronoi(32, 32, &cfg, 42);
        assert!(g.is_dynamic());
        assert_eq!(g.distinct_plate_count(), 8);
        assert_eq!(g.layout_name, "voronoi_seed42_n8");
    }
}
