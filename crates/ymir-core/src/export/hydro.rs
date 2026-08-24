//! JSON vector layers for the `.ymir` container — rivers and lakes.
//!
//! Serializes the EXISTING [`crate::tectonics_c1::drainage`] outputs (no
//! re-computation) with `serde_json`. Coordinates are in erosion-grid cells,
//! the same row-major space as the container rasters (see
//! [`crate::export::container`] for the orientation invariant).

use serde::Serialize;

use crate::tectonics_c1::drainage::{C1DrainageResult, Navigability};
use crate::terrain::flow::RiverSegment;

/// Coordinate space tag written into `rivers.json` so a consumer never guesses.
pub const COORDINATE_SPACE: &str = "erosion_grid_cells";

/// A river segment enriched with its parallel drainage-area / navigability
/// (which live in `C1DrainageResult` as arrays parallel to `rivers.segments`).
#[derive(Serialize)]
struct RiverSegmentView<'a> {
    #[serde(flatten)]
    segment: &'a RiverSegment,
    /// Upstream drainage area in km² (`segment_drainage_km2[i]`).
    drainage_km2: f32,
    /// Navigability class (`segment_navigability[i]`).
    navigability: Navigability,
    /// Mean discharge in m³/s (`segment_discharge_m3s[i]`) — runoff × catchment.
    discharge_m3s: f32,
    /// Bankfull channel width in metres (`segment_width_m[i]`), hydraulic geometry
    /// `w = 5·Q^0.5` on the DISCHARGE. A rendering hint (stroke width) — topology
    /// stays continuous regardless of any display cutoff, and at production scale
    /// most reaches are sub-cell, so the consumer must render this as a stroke.
    width_m: f32,
    /// Long profile — bed elevation (m) at each `segment.points`, upstream→downstream
    /// (`segment_profile_m[i]`).
    profile_m: &'a [f32],
}

#[derive(Serialize)]
struct RiversView<'a> {
    /// `"erosion_grid_cells"` — segment points are `(x, y)` cell indices.
    coordinate_space: &'static str,
    segments: Vec<RiverSegmentView<'a>>,
}

/// Serialize `drainage.rivers` as `rivers.json` bytes, attaching each segment's
/// parallel `drainage_km2` and `navigability`. Deterministic (segment order).
pub fn rivers_json(drainage: &C1DrainageResult) -> Vec<u8> {
    let segments = drainage
        .rivers
        .segments
        .iter()
        .enumerate()
        .map(|(i, segment)| RiverSegmentView {
            segment,
            // The arrays are parallel by contract; `get` stays defensive.
            drainage_km2: drainage.segment_drainage_km2.get(i).copied().unwrap_or(0.0),
            navigability: drainage
                .segment_navigability
                .get(i)
                .copied()
                .unwrap_or(Navigability::NonNavigable),
            discharge_m3s: drainage.segment_discharge_m3s.get(i).copied().unwrap_or(0.0),
            width_m: drainage.segment_width_m.get(i).copied().unwrap_or(0.0),
            profile_m: drainage.segment_profile_m.get(i).map(|v| v.as_slice()).unwrap_or(&[]),
        })
        .collect();
    let view = RiversView { coordinate_space: COORDINATE_SPACE, segments };
    serde_json::to_vec_pretty(&view).expect("rivers view serialization is infallible")
}

/// Serialize `drainage.lakes` (the enriched `C1Lake` list) as `lakes.json`
/// bytes. Deterministic (lake order).
pub fn lakes_json(drainage: &C1DrainageResult) -> Vec<u8> {
    serde_json::to_vec_pretty(&drainage.lakes).expect("lakes serialization is infallible")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::GridF32;
    use crate::lakes::detection::Lake;
    use crate::tectonics_c1::drainage::{C1Lake, LakeType};
    use crate::terrain::flow::{FlowResult, RiverNetwork, RiverSegment};
    use serde_json::Value;

    fn synthetic_drainage() -> C1DrainageResult {
        let seg = RiverSegment {
            points: vec![(1, 1), (2, 2), (3, 3)],
            strahler_order: 4,
            avg_flow: 12.0,
            max_flow: 30.0,
            basin_id: 7,
            upstream: vec![],
            downstream: None,
        };
        let lake = C1Lake {
            base: Lake {
                id: 3,
                surface_elevation: 0.6,
                max_depth: 0.05,
                area: 42,
                basin_id: 7,
                outlet: (5, 6),
                shallow: false,
            },
            level_m: 1130.0,
            depth_m: 90.0,
            area_km2: 12.5,
            lake_type: LakeType::Endorheic,
        };
        C1DrainageResult {
            flow: FlowResult {
                filled: GridF32::new(4, 4, 0.0),
                direction: vec![0; 16],
                accumulation: GridF32::new(4, 4, 1.0),
                basins: vec![0; 16],
                num_basins: 1,
            },
            rivers: RiverNetwork { segments: vec![seg] },
            segment_drainage_km2: vec![55_000.0],
            segment_navigability: vec![Navigability::Ship],
            segment_discharge_m3s: vec![520.0],
            segment_width_m: vec![114.0],
            segment_profile_m: vec![vec![1200.0, 1150.0, 1130.0]],
            lakes: vec![lake],
            lake_map: vec![0; 16],
            width: 4,
            height: 4,
        }
    }

    /// rivers.json carries the coordinate space and, per segment, the strahler
    /// order zipped with its navigability class.
    #[test]
    fn rivers_view_round_trips_strahler_and_navigability() {
        let d = synthetic_drainage();
        let v: Value = serde_json::from_slice(&rivers_json(&d)).unwrap();
        assert_eq!(v["coordinate_space"], "erosion_grid_cells");
        let seg = &v["segments"][0];
        assert_eq!(seg["strahler_order"], 4, "flattened segment keeps strahler_order");
        assert_eq!(seg["navigability"], "Ship", "navigability zipped onto the segment");
        assert_eq!(seg["drainage_km2"], 55000.0);
        assert_eq!(seg["basin_id"], 7);
    }

    /// lakes.json round-trips back into `Vec<C1Lake>` with its typed fields.
    #[test]
    fn lakes_json_round_trips() {
        let d = synthetic_drainage();
        let lakes: Vec<C1Lake> = serde_json::from_slice(&lakes_json(&d)).unwrap();
        assert_eq!(lakes.len(), 1);
        assert_eq!(lakes[0].lake_type, LakeType::Endorheic);
        assert_eq!(lakes[0].base.id, 3);
        assert!((lakes[0].area_km2 - 12.5).abs() < 1e-6);
    }
}
