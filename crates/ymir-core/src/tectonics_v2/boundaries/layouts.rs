//! Synthetic boundary-layout generators for Step 5.
//!
//! These layouts are **prescribed statically**: they produce a
//! plate-type field and a boundary-flag field that do not change
//! over the run (Step 6 handles dynamic reclassification). Each
//! layout is designed to isolate one physical effect of the
//! source/sink machinery.

use super::boundary_flag::{
    BoundaryConfig, BoundaryFlag, BoundaryFlagField, BoundaryRates,
};
use super::plate_type::{PlateType, PlateTypeField};

/// A named pair (plate_types, flags) that the harness consumes.
///
/// `name` is reproducibility metadata: the physics report prints it
/// in the config dump so the layout is unambiguous.
#[derive(Clone, Debug)]
pub struct BoundaryLayout {
    pub plate_types: PlateTypeField,
    pub flags: BoundaryFlagField,
    pub name: &'static str,
}

impl BoundaryLayout {
    pub fn nx(&self) -> usize { self.plate_types.nx() }
    pub fn ny(&self) -> usize { self.plate_types.ny() }

    /// Wrap into a [`BoundaryConfig::Enabled`] with the given rates,
    /// using the Step 5 back-compat builder
    /// [`BoundaryConfig::enabled_static`]. The default recycling
    /// mode is `Open`, preserving the Step 5 per-cell rate-based
    /// source/sink pipeline. Callers wanting Step 6's Closed mode
    /// should build via [`BoundaryConfig::enabled_voronoi_closed`]
    /// or manually swap `recycling_mode` after construction.
    pub fn into_config(self, rates: BoundaryRates) -> BoundaryConfig {
        BoundaryConfig::enabled_static(self.plate_types, self.flags, self.name, rates)
    }
}

/// Horizontal oceanic strip — the canonical Step 5 baseline layout.
///
/// ```text
///    continental
///   -------- Rift row ---------     (north edge of the strip)
///    oceanic (S̃ ≈ 0.2 at steady state)
///   -- OceanicSubduction row --     (south edge of the strip)
///    continental (arc side)
/// ```
///
/// The strip occupies the central third of the grid (rows
/// `[ny/3, 2ny/3)`). The **north row** of the strip is flagged
/// `Rift` (spreading centre), the **south row** is flagged
/// `OceanicSubduction`. Cells one row north of the rift and one row
/// south of the subduction are continental interior. This layout
/// produces a mass source at the rift and a mass sink at the
/// subduction, separated by oceanic interior with no source.
pub fn horizontal_oceanic_strip(nx: usize, ny: usize) -> BoundaryLayout {
    let mut plate_types = PlateTypeField::filled(nx, ny, PlateType::Continental);
    let mut flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);

    let strip_south = ny / 3;
    let strip_north = 2 * ny / 3;
    for j in strip_south..strip_north {
        for i in 0..nx {
            plate_types.set(i, j, PlateType::Oceanic);
        }
    }
    // Rift at the north edge of the strip; subduction at the south
    // edge. Spanning the full x-range of each row.
    if strip_north > 0 {
        let rift_j = strip_north - 1;
        for i in 0..nx {
            flags.set(i, rift_j, BoundaryFlag::Rift);
        }
    }
    if strip_south < ny {
        let sub_j = strip_south;
        for i in 0..nx {
            flags.set(i, sub_j, BoundaryFlag::OceanicSubduction);
        }
    }

    BoundaryLayout {
        plate_types,
        flags,
        name: "horizontal_oceanic_strip",
    }
}

/// Vertical rift line through an otherwise continental domain.
/// Tests intracontinental rifting (`Q_rift-v`, not `Q_spread`).
pub fn vertical_rift_line(nx: usize, ny: usize) -> BoundaryLayout {
    let plate_types = PlateTypeField::filled(nx, ny, PlateType::Continental);
    let mut flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
    let rift_i = nx / 2;
    for j in 0..ny {
        flags.set(rift_i, j, BoundaryFlag::Rift);
    }
    BoundaryLayout {
        plate_types,
        flags,
        name: "vertical_rift_line",
    }
}

/// Balanced sub/spread layout: one rift row and one oceanic-
/// subduction row with **matched count** of cells, so if the
/// continental arc wrapping is absent and the rates are equal then
/// `∫ Q_sub + ∫ Q_spread = 0` exactly.
///
/// For exact mass balance at the run level, the caller is expected
/// to disable `Q_arc`, `Q_coll-v`, `Q_rift-v` via the rates (set
/// them to zero) and to tune `k_sub · |Δv_conv|_mean ≈ k_spread`.
/// The v2_balanced_mass_balance test relies on this structure plus
/// a zero-velocity initial condition that the advection preserves.
pub fn balanced_sub_spread(nx: usize, ny: usize) -> BoundaryLayout {
    let plate_types = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
    let mut flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
    // Rift on row ny/3, subduction on row 2ny/3 — same cell count
    // per row on the square periodic domain.
    let rift_j = ny / 3;
    let sub_j = (2 * ny) / 3;
    for i in 0..nx {
        flags.set(i, rift_j, BoundaryFlag::Rift);
        flags.set(i, sub_j, BoundaryFlag::OceanicSubduction);
    }
    BoundaryLayout {
        plate_types,
        flags,
        name: "balanced_sub_spread",
    }
}

/// Continental collision band at the middle of the grid.
/// All-continental domain with a single `ContinentalCollision` row.
pub fn continental_collision_band(nx: usize, ny: usize) -> BoundaryLayout {
    let plate_types = PlateTypeField::filled(nx, ny, PlateType::Continental);
    let mut flags = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
    let coll_j = ny / 2;
    for i in 0..nx {
        flags.set(i, coll_j, BoundaryFlag::ContinentalCollision);
    }
    BoundaryLayout {
        plate_types,
        flags,
        name: "continental_collision_band",
    }
}

/// CLI token → layout builder. Unknown tokens error out with a
/// list of available names.
pub fn build_layout(name: &str, nx: usize, ny: usize) -> Result<BoundaryLayout, String> {
    match name {
        "horizontal_oceanic_strip" => Ok(horizontal_oceanic_strip(nx, ny)),
        "vertical_rift_line" => Ok(vertical_rift_line(nx, ny)),
        "balanced_sub_spread" => Ok(balanced_sub_spread(nx, ny)),
        "continental_collision_band" => Ok(continental_collision_band(nx, ny)),
        other => Err(format!(
            "unknown layout '{}'; expected one of horizontal_oceanic_strip, \
             vertical_rift_line, balanced_sub_spread, continental_collision_band",
            other,
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn horizontal_strip_has_oceanic_band() {
        let l = horizontal_oceanic_strip(12, 9);
        // Rows [3, 6) are oceanic.
        for i in 0..12 {
            assert_eq!(l.plate_types.get(i, 3), PlateType::Oceanic);
            assert_eq!(l.plate_types.get(i, 4), PlateType::Oceanic);
            assert_eq!(l.plate_types.get(i, 5), PlateType::Oceanic);
            assert_eq!(l.plate_types.get(i, 0), PlateType::Continental);
            assert_eq!(l.plate_types.get(i, 8), PlateType::Continental);
        }
    }

    #[test]
    fn horizontal_strip_has_rift_on_north_edge() {
        let l = horizontal_oceanic_strip(12, 9);
        // strip_north = 6; rift at j = 5.
        for i in 0..12 {
            assert_eq!(l.flags.get(i, 5), BoundaryFlag::Rift);
            assert_eq!(l.flags.get(i, 3), BoundaryFlag::OceanicSubduction);
        }
    }

    #[test]
    fn vertical_rift_is_single_column() {
        let l = vertical_rift_line(8, 8);
        let c = 4;
        for j in 0..8 {
            assert_eq!(l.flags.get(c, j), BoundaryFlag::Rift);
            // Column c-1 is None.
            if c > 0 {
                assert_eq!(l.flags.get(c - 1, j), BoundaryFlag::None);
            }
        }
        for j in 0..8 {
            for i in 0..8 {
                assert_eq!(l.plate_types.get(i, j), PlateType::Continental);
            }
        }
    }

    #[test]
    fn balanced_layout_has_matched_rift_and_subduction_counts() {
        let l = balanced_sub_spread(16, 12);
        let mut rift = 0usize;
        let mut sub = 0usize;
        for j in 0..12 {
            for i in 0..16 {
                match l.flags.get(i, j) {
                    BoundaryFlag::Rift => rift += 1,
                    BoundaryFlag::OceanicSubduction => sub += 1,
                    _ => {}
                }
            }
        }
        assert_eq!(rift, sub);
        assert_eq!(rift, 16);
    }

    #[test]
    fn continental_collision_has_single_row() {
        let l = continental_collision_band(8, 8);
        for i in 0..8 {
            assert_eq!(l.flags.get(i, 4), BoundaryFlag::ContinentalCollision);
        }
        for j in 0..8 {
            if j != 4 {
                for i in 0..8 {
                    assert_eq!(l.flags.get(i, j), BoundaryFlag::None);
                }
            }
        }
    }

    #[test]
    fn build_layout_dispatches() {
        assert!(build_layout("horizontal_oceanic_strip", 8, 8).is_ok());
        assert!(build_layout("balanced_sub_spread", 8, 8).is_ok());
        assert!(build_layout("garbage", 8, 8).is_err());
    }
}
