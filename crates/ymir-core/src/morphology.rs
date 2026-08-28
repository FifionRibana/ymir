//! Land-mask morphology metrics — distinguish continental **masses**
//! from drainage **filaments** (piste 4).
//!
//! A scalar land *fraction* does NOT validate spatial structure: 30 %
//! land as a few compact landmasses and 30 % land as a filamentous
//! drainage network are identical by area. The Issue #141 cap=0.92
//! calibration passed every fraction/convergence gate yet produced
//! filaments (perim/area ~1.1, ~48 components); only a visual caught
//! it. These metrics are the gate that would have caught it, and are
//! the permanent morphological acceptance for C1 land output.
//!
//! All metrics use **4-neighbour** connectivity on a row-major
//! `nx × ny` boolean mask (`true` = land).

/// Morphology of a boolean land mask. Reference values (C1 gallery
/// production, seed 42, 64²): masses ≈ `{ area 0.33, perim/area 0.515,
/// n_components 11, largest 0.625 }`; filaments (reclassify workflow)
/// ≈ `{ 0.30, 1.10, 48, 0.378 }`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LandMorphology {
    /// Land cells / total cells.
    pub area_fraction: f64,
    /// Boundary edges / land cells (an edge is a land-cell face whose
    /// 4-neighbour is non-land or off-grid). **The key filament
    /// discriminator**: compact masses → low (~0.4–0.5); thin
    /// filaments → high (~1.0+). Small/sparse masks read high too
    /// (boundary-dominated), so compare at a fixed scale.
    pub perimeter_over_area: f64,
    /// Number of 4-connected land components.
    pub n_components: usize,
    /// Largest component's cells / total land cells. One dominant
    /// landmass → high (~0.6+); fragmented → low (~0.3). `0.0` when
    /// there is no land.
    pub largest_component_fraction: f64,
}

/// Compute [`LandMorphology`] for a row-major `nx × ny` land mask
/// (`mask[j*nx + i]`). O(N) time, O(N) scratch. Deterministic.
pub fn land_morphology(mask: &[bool], nx: usize, ny: usize) -> LandMorphology {
    debug_assert_eq!(mask.len(), nx * ny, "mask length must be nx*ny");
    let total = nx * ny;
    let area = mask.iter().filter(|&&b| b).count();
    if area == 0 {
        return LandMorphology {
            area_fraction: 0.0,
            perimeter_over_area: 0.0,
            n_components: 0,
            largest_component_fraction: 0.0,
        };
    }

    let idx = |i: usize, j: usize| j * nx + i;
    // Perimeter: count land-cell faces adjacent to non-land or off-grid.
    let mut perimeter = 0usize;
    for j in 0..ny {
        for i in 0..nx {
            if !mask[idx(i, j)] {
                continue;
            }
            for (di, dj) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
                let (ni, nj) = (i as i32 + di, j as i32 + dj);
                if ni < 0
                    || nj < 0
                    || ni >= nx as i32
                    || nj >= ny as i32
                    || !mask[idx(ni as usize, nj as usize)]
                {
                    perimeter += 1;
                }
            }
        }
    }

    // Connected components (4-neighbour, iterative flood fill).
    let mut seen = vec![false; total];
    let mut n_components = 0usize;
    let mut largest = 0usize;
    let mut stack: Vec<usize> = Vec::new();
    for start in 0..total {
        if !mask[start] || seen[start] {
            continue;
        }
        n_components += 1;
        let mut size = 0usize;
        stack.push(start);
        seen[start] = true;
        while let Some(c) = stack.pop() {
            size += 1;
            let (i, j) = (c % nx, c / nx);
            for (di, dj) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
                let (ni, nj) = (i as i32 + di, j as i32 + dj);
                if ni >= 0 && nj >= 0 && ni < nx as i32 && nj < ny as i32 {
                    let k = idx(ni as usize, nj as usize);
                    if mask[k] && !seen[k] {
                        seen[k] = true;
                        stack.push(k);
                    }
                }
            }
        }
        largest = largest.max(size);
    }

    LandMorphology {
        area_fraction: area as f64 / total as f64,
        perimeter_over_area: perimeter as f64 / area as f64,
        n_components,
        largest_component_fraction: largest as f64 / area as f64,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build an `nx × ny` mask from a closure `(i,j) -> bool`.
    fn mask_from(nx: usize, ny: usize, f: impl Fn(usize, usize) -> bool) -> Vec<bool> {
        (0..nx * ny).map(|k| f(k % nx, k / nx)).collect()
    }

    #[test]
    fn empty_mask() {
        let m = land_morphology(&vec![false; 64], 8, 8);
        assert_eq!(m.area_fraction, 0.0);
        assert_eq!(m.n_components, 0);
        assert_eq!(m.perimeter_over_area, 0.0);
        assert_eq!(m.largest_component_fraction, 0.0);
    }

    #[test]
    fn single_cell() {
        // One land cell in an 8×8 grid: area 1, all 4 faces are
        // boundary, one component.
        let m = land_morphology(&mask_from(8, 8, |i, j| i == 3 && j == 3), 8, 8);
        assert_eq!(m.area_fraction, 1.0 / 64.0);
        assert_eq!(m.perimeter_over_area, 4.0);
        assert_eq!(m.n_components, 1);
        assert_eq!(m.largest_component_fraction, 1.0);
    }

    #[test]
    fn compact_block_beats_scattered_on_every_axis() {
        // THE discrimination test: 16 land cells as one 4×4 block vs
        // 16 isolated cells, same grid + same area.
        let block = land_morphology(&mask_from(16, 16, |i, j| i < 4 && j < 4), 16, 16);
        // 16 isolated cells on a coarse lattice (every 4th cell) so
        // none touch.
        let scattered =
            land_morphology(&mask_from(16, 16, |i, j| i % 4 == 0 && j % 4 == 0), 16, 16);

        assert_eq!(block.area_fraction, scattered.area_fraction); // same area
        // Block: 1 component, scattered: 16.
        assert_eq!(block.n_components, 1);
        assert_eq!(scattered.n_components, 16);
        // Block perim/area = 16/16 = 1.0; scattered = 64/16 = 4.0.
        assert!(block.perimeter_over_area < scattered.perimeter_over_area);
        assert_eq!(block.perimeter_over_area, 1.0);
        assert_eq!(scattered.perimeter_over_area, 4.0);
        // Block: one dominant mass; scattered: 1/16 each.
        assert_eq!(block.largest_component_fraction, 1.0);
        assert!((scattered.largest_component_fraction - 1.0 / 16.0).abs() < 1e-12);
    }

    #[test]
    fn two_separate_blobs() {
        // Two 2×2 blocks far apart → 2 components, largest = half.
        let m = land_morphology(
            &mask_from(16, 16, |i, j| (i < 2 && j < 2) || (i >= 14 && j >= 14)),
            16,
            16,
        );
        assert_eq!(m.n_components, 2);
        assert!((m.largest_component_fraction - 0.5).abs() < 1e-12);
    }

    #[test]
    fn larger_compact_blob_has_lower_perim_over_area() {
        // Scale sanity: a bigger compact blob has lower perim/area
        // (boundary grows slower than area) — why the gallery's 64²
        // masses read ~0.5, not the 1.0 of a 4×4.
        let small = land_morphology(&mask_from(32, 32, |i, j| i < 4 && j < 4), 32, 32);
        let big = land_morphology(&mask_from(32, 32, |i, j| i < 16 && j < 16), 32, 32);
        assert!(big.perimeter_over_area < small.perimeter_over_area);
        // 16×16 block: perimeter 64, area 256 → 0.25.
        assert_eq!(big.perimeter_over_area, 0.25);
    }
}
