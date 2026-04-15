//! Lake detection from pit-filled vs original heightmap comparison.

use std::collections::VecDeque;

use serde::{Deserialize, Serialize};

use crate::grid::GridF32;
use crate::terrain::flow::DIR_NONE;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct LakeConfig {
    /// Minimum depth (filled - original) to classify as lake. Default: 0.001.
    pub min_depth: f32,
    /// Minimum area in pixels to keep a lake. Default: 20.
    pub min_area: usize,
}

impl Default for LakeConfig {
    fn default() -> Self {
        Self { min_depth: 0.001, min_area: 20 }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Lake {
    pub id: u32,
    pub surface_elevation: f32,
    pub max_depth: f32,
    pub area: usize,
    pub basin_id: u32,
    pub outlet: (u32, u32),
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct LakeResult {
    pub lake_map: Vec<u32>,
    pub lakes: Vec<Lake>,
    pub width: usize,
    pub height: usize,
}

/// Detect lakes by comparing eroded heightmap with pit-filled heightmap.
pub fn detect_lakes(
    eroded: &GridF32,
    filled: &GridF32,
    direction: &[u8],
    basins: &[u32],
    config: &LakeConfig,
) -> LakeResult {
    let w = eroded.width;
    let h = eroded.height;
    let n = w * h;

    // Step 1: identify potential lake cells
    let is_lake_cell: Vec<bool> =
        (0..n).map(|i| (filled.data[i] - eroded.data[i]) > config.min_depth).collect();

    // Step 2: connected component labeling (BFS, 8-connectivity, periodic)
    let mut lake_map = vec![0u32; n];
    let mut next_id = 1u32;

    for start in 0..n {
        if !is_lake_cell[start] || lake_map[start] != 0 {
            continue;
        }

        let id = next_id;
        next_id += 1;
        lake_map[start] = id;

        let mut queue = VecDeque::new();
        queue.push_back(start);

        while let Some(cur) = queue.pop_front() {
            let cx = cur % w;
            let cy = cur / w;

            for dy in -1i32..=1 {
                for dx in -1i32..=1 {
                    if dx == 0 && dy == 0 {
                        continue;
                    }
                    let nx = ((cx as i32 + dx) % w as i32 + w as i32) as usize % w;
                    let ny = ((cy as i32 + dy) % h as i32 + h as i32) as usize % h;
                    let nidx = ny * w + nx;

                    if is_lake_cell[nidx] && lake_map[nidx] == 0 {
                        lake_map[nidx] = id;
                        queue.push_back(nidx);
                    }
                }
            }
        }
    }

    // Step 3: compute per-lake stats and filter by min_area
    let num_labels = (next_id - 1) as usize;
    let mut areas = vec![0usize; num_labels];
    let mut surface_sums = vec![0.0f64; num_labels];
    let mut min_eroded = vec![f32::INFINITY; num_labels];
    let mut basin_ids = vec![0u32; num_labels];

    for i in 0..n {
        let lid = lake_map[i];
        if lid == 0 {
            continue;
        }
        let idx = (lid - 1) as usize;
        areas[idx] += 1;
        surface_sums[idx] += filled.data[i] as f64;
        if eroded.data[i] < min_eroded[idx] {
            min_eroded[idx] = eroded.data[i];
        }
        if basin_ids[idx] == 0 {
            basin_ids[idx] = basins[i];
        }
    }

    // Find outlets: lake cell whose D8 direction points to a non-lake cell (or different lake)
    let mut outlets: Vec<Option<(u32, u32, f32)>> = vec![None; num_labels]; // (x, y, filled_h)

    for i in 0..n {
        let lid = lake_map[i];
        if lid == 0 {
            continue;
        }
        let d = direction[i];
        if d == DIR_NONE {
            continue;
        }
        let cx = i % w;
        let cy = i / w;
        let dx = crate::terrain::flow::D8_DX[d as usize];
        let dy = crate::terrain::flow::D8_DY[d as usize];
        let nx = ((cx as i32 + dx) % w as i32 + w as i32) as usize % w;
        let ny = ((cy as i32 + dy) % h as i32 + h as i32) as usize % h;
        let nidx = ny * w + nx;

        if lake_map[nidx] != lid {
            // This cell is an outlet (drains out of the lake)
            let idx = (lid - 1) as usize;
            let h_here = filled.data[i];
            if outlets[idx].is_none() || h_here < outlets[idx].unwrap().2 {
                outlets[idx] = Some((cx as u32, cy as u32, h_here));
            }
        }
    }

    // Build lake list, filtering by min_area
    let mut lakes = Vec::new();
    let mut id_remap = vec![0u32; num_labels]; // old_idx -> new_id (0 = removed)
    let mut new_id = 1u32;

    for idx in 0..num_labels {
        if areas[idx] < config.min_area {
            continue;
        }

        let surface_elevation = (surface_sums[idx] / areas[idx] as f64) as f32;
        let max_depth = surface_elevation - min_eroded[idx];
        let outlet = outlets[idx].map(|(x, y, _)| (x, y)).unwrap_or_else(|| {
            // Fallback: pick the cell with highest accumulation
            // (not passed here, so just pick center of mass)
            (0, 0)
        });

        id_remap[idx] = new_id;
        lakes.push(Lake {
            id: new_id,
            surface_elevation,
            max_depth,
            area: areas[idx],
            basin_id: basin_ids[idx],
            outlet,
        });
        new_id += 1;
    }

    // Remap lake_map: remove filtered lakes, renumber survivors
    for v in &mut lake_map {
        if *v > 0 {
            let idx = (*v - 1) as usize;
            *v = id_remap[idx];
        }
    }

    LakeResult { lake_map, lakes, width: w, height: h }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pit_detected_as_lake() {
        let n = 32;
        let mut eroded = GridF32::new(n, n, 0.3);
        // Create a pit at center
        for j in 14..18 {
            for i in 14..18 {
                eroded.set(i, j, 0.1);
            }
        }
        // Filled version: pit raised to surrounding level
        let mut filled = eroded.clone();
        for j in 14..18 {
            for i in 14..18 {
                filled.set(i, j, 0.3);
            }
        }

        let direction = vec![0u8; n * n]; // dummy directions
        let basins = vec![1u32; n * n];
        let config = LakeConfig { min_depth: 0.001, min_area: 1 };

        let result = detect_lakes(&eroded, &filled, &direction, &basins, &config);
        assert!(!result.lakes.is_empty(), "Should detect a lake");
        assert_eq!(result.lakes[0].area, 16); // 4x4 pit
        assert!((result.lakes[0].max_depth - 0.2).abs() < 0.01);
    }

    #[test]
    fn flat_terrain_no_lakes() {
        let n = 32;
        let hmap = GridF32::new(n, n, 0.5);
        let direction = vec![0u8; n * n];
        let basins = vec![1u32; n * n];
        let config = LakeConfig::default();

        let result = detect_lakes(&hmap, &hmap, &direction, &basins, &config);
        assert!(result.lakes.is_empty(), "Flat terrain should have no lakes");
    }

    #[test]
    fn min_area_filter() {
        let n = 32;
        let mut eroded = GridF32::new(n, n, 0.3);
        let mut filled = eroded.clone();

        // Small pit: 4 cells (should be filtered)
        for j in 5..7 {
            for i in 5..7 {
                eroded.set(i, j, 0.1);
                filled.set(i, j, 0.3);
            }
        }
        // Large pit: 36 cells (should survive)
        for j in 20..26 {
            for i in 20..26 {
                eroded.set(i, j, 0.1);
                filled.set(i, j, 0.3);
            }
        }

        let direction = vec![0u8; n * n];
        let basins = vec![1u32; n * n];
        let config = LakeConfig { min_depth: 0.001, min_area: 10 };

        let result = detect_lakes(&eroded, &filled, &direction, &basins, &config);
        assert_eq!(result.lakes.len(), 1, "Only the large pit should survive filtering");
        assert_eq!(result.lakes[0].area, 36);
    }
}
