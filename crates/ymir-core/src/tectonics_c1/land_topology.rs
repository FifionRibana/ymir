//! Land-topology diagnostics (M1) — judge a seed as an "island continent".
//!
//! Connected-component labelling of the land mask (`altitude_norm > sea_level`)
//! on the coarse C1 field, with **periodic wrap in both axes** (the domain is a
//! torus). Reports the number of landmasses and, for the largest, its area, its
//! bounding box, and whether it connects to itself across the x/y seam. A
//! landmass that WRAPS the torus is not an island surrounded by ocean and the
//! seed should be rejected for that use case.
//!
//! Resolution-independent: run on the coarse altitude field (before FBM/erosion),
//! area via [`crate::tectonics_c1::production_upscale::c1_cell_area_km2`].

use crate::grid::GridF32;

use super::production_upscale::{C1_DOMAIN_KM, c1_cell_area_km2};

/// Metrics for the land mask of a coarse altitude field.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LandTopology {
    /// Number of distinct landmasses (periodic 4-connected components of land).
    pub num_landmasses: usize,
    /// Largest landmass — cell count.
    pub largest_cells: usize,
    /// Largest landmass — area in km².
    pub largest_area_km2: f32,
    /// Largest landmass — fraction of the whole domain area `[0,1]`.
    pub largest_area_frac: f32,
    /// Largest landmass connects to itself across the x seam (wraps the torus in
    /// x). A `true` here means the mass spans the domain → NOT an island.
    pub wraps_x: bool,
    /// Largest landmass connects to itself across the y seam.
    pub wraps_y: bool,
    /// Bounding box of the largest landmass in cells (inclusive): `(min, max)`.
    /// Meaningless as an extent when the corresponding axis wraps.
    pub bbox_min: (usize, usize),
    pub bbox_max: (usize, usize),
    /// Bounding box `(width_km, height_km)` of the largest landmass.
    pub bbox_km: (f32, f32),
    /// Total emerged fraction of the field (all land cells / all cells).
    pub emerged_fraction: f32,
}

impl LandTopology {
    /// An empty (all-ocean) field.
    fn empty() -> Self {
        Self {
            num_landmasses: 0,
            largest_cells: 0,
            largest_area_km2: 0.0,
            largest_area_frac: 0.0,
            wraps_x: false,
            wraps_y: false,
            bbox_min: (0, 0),
            bbox_max: (0, 0),
            bbox_km: (0.0, 0.0),
            emerged_fraction: 0.0,
        }
    }
}

// ── Union-find (disjoint set) ───────────────────────────────────────────────

struct DisjointSet {
    parent: Vec<usize>,
}

impl DisjointSet {
    fn new(n: usize) -> Self {
        Self { parent: (0..n).collect() }
    }
    fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]]; // path halving
            x = self.parent[x];
        }
        x
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.parent[ra] = rb;
        }
    }
}

/// Compute the land-topology metrics of an altitude field. Land = cell value
/// strictly `> sea_level_norm`. 4-connectivity with periodic wrap in x and y.
/// The km² scale uses the field's own width (a square coarse torus of
/// [`C1_DOMAIN_KM`] on a side).
pub fn land_topology(altitude_norm: &GridF32, sea_level_norm: f32) -> LandTopology {
    let (w, h) = (altitude_norm.width, altitude_norm.height);
    let n = w * h;
    if n == 0 {
        return LandTopology::empty();
    }

    let is_land = |k: usize| altitude_norm.data[k] > sea_level_norm;
    let land_count = (0..n).filter(|&k| is_land(k)).count();
    let emerged_fraction = land_count as f32 / n as f32;
    if land_count == 0 {
        return LandTopology { emerged_fraction, ..LandTopology::empty() };
    }

    // Union land cells with their right and down neighbours (periodic). Unioning
    // both directions over every cell covers all 4-adjacencies, seam included.
    let mut ds = DisjointSet::new(n);
    for y in 0..h {
        for x in 0..w {
            let k = y * w + x;
            if !is_land(k) {
                continue;
            }
            let rx = (x + 1) % w;
            let kr = y * w + rx;
            if is_land(kr) {
                ds.union(k, kr);
            }
            let dy = (y + 1) % h;
            let kd = dy * w + x;
            if is_land(kd) {
                ds.union(k, kd);
            }
        }
    }

    // Component sizes.
    let mut sizes: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for k in 0..n {
        if is_land(k) {
            *sizes.entry(ds.find(k)).or_insert(0) += 1;
        }
    }
    let num_landmasses = sizes.len();
    // Largest by size; tie-break on the smallest root index for determinism.
    let largest_root =
        sizes.iter().max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0))).map(|(&root, _)| root).unwrap();
    let largest_cells = sizes[&largest_root];

    // Bounding box + wrap detection for the largest component.
    let (mut xmin, mut xmax, mut ymin, mut ymax) = (usize::MAX, 0usize, usize::MAX, 0usize);
    let (mut wraps_x, mut wraps_y) = (false, false);
    for y in 0..h {
        for x in 0..w {
            let k = y * w + x;
            if is_land(k) && ds.find(k) == largest_root {
                xmin = xmin.min(x);
                xmax = xmax.max(x);
                ymin = ymin.min(y);
                ymax = ymax.max(y);
            }
        }
    }
    // Wrap = the seam edge is INTERNAL to the largest mass: both the last and the
    // first cell of a row/column are land and belong to it (they are unioned).
    for y in 0..h {
        let (k_last, k_first) = (y * w + (w - 1), y * w);
        if is_land(k_last) && is_land(k_first) && ds.find(k_last) == largest_root {
            wraps_x = true;
            break;
        }
    }
    for x in 0..w {
        let (k_last, k_first) = ((h - 1) * w + x, x);
        if is_land(k_last) && is_land(k_first) && ds.find(k_last) == largest_root {
            wraps_y = true;
            break;
        }
    }

    let cell_km2 = c1_cell_area_km2(w);
    let km_per_cell = C1_DOMAIN_KM / w as f32;
    let domain_km2 = C1_DOMAIN_KM * C1_DOMAIN_KM;
    LandTopology {
        num_landmasses,
        largest_cells,
        largest_area_km2: largest_cells as f32 * cell_km2,
        largest_area_frac: (largest_cells as f32 * cell_km2) / domain_km2,
        wraps_x,
        wraps_y,
        bbox_min: (xmin, ymin),
        bbox_max: (xmax, ymax),
        bbox_km: ((xmax - xmin + 1) as f32 * km_per_cell, (ymax - ymin + 1) as f32 * km_per_cell),
        emerged_fraction,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a coarse field with a rectangular land block `[x0,x1)×[y0,y1)` (norm
    /// 0.8) on an ocean background (norm 0.2). Ranges may exceed the grid to make
    /// a block that touches/crosses an edge (caller keeps them in-bounds).
    fn field_with_block(w: usize, h: usize, x0: usize, x1: usize, y0: usize, y1: usize) -> GridF32 {
        let mut d = vec![0.2f32; w * h];
        for y in y0..y1 {
            for x in x0..x1 {
                d[y * w + x] = 0.8;
            }
        }
        GridF32::from_vec(w, h, d)
    }

    #[test]
    fn compact_blob_is_one_nonwrapping_landmass() {
        let f = field_with_block(20, 20, 6, 12, 6, 12); // interior 6×6 block
        let t = land_topology(&f, 0.5);
        assert_eq!(t.num_landmasses, 1);
        assert!(!t.wraps_x && !t.wraps_y, "compact interior blob must not wrap");
        assert_eq!(t.largest_cells, 36);
        assert_eq!((t.bbox_min, t.bbox_max), ((6, 6), (11, 11)));
        assert!((t.emerged_fraction - 36.0 / 400.0).abs() < 1e-6);
    }

    #[test]
    fn block_touching_both_x_edges_wraps_x() {
        // Land in columns 0..3 and 17..20 but only rows 3..7 (interior in y), so
        // it crosses the x seam but leaves y-edge rows ocean → wraps_x, not y.
        let mut f = field_with_block(20, 10, 0, 3, 3, 7);
        for y in 3..7 {
            for x in 17..20 {
                f.data[y * 20 + x] = 0.8;
            }
        }
        let t = land_topology(&f, 0.5);
        assert!(t.wraps_x, "land spanning the x seam must be flagged wraps_x");
        assert!(!t.wraps_y, "y-edge rows are ocean → must not wrap y");
        assert_eq!(t.num_landmasses, 1, "the two edge strips are one mass via the seam");
    }

    #[test]
    fn block_touching_both_y_edges_wraps_y() {
        // Rows 0..2 and 8..10 land, with an ocean gap in x so it can't wrap x.
        let mut f = field_with_block(12, 10, 3, 9, 0, 2);
        for y in 8..10 {
            for x in 3..9 {
                f.data[y * 12 + x] = 0.8;
            }
        }
        let t = land_topology(&f, 0.5);
        assert!(t.wraps_y, "land spanning the y seam must be flagged wraps_y");
        assert!(!t.wraps_x, "interior-x block must not wrap x");
    }

    #[test]
    fn all_ocean_reports_no_land() {
        let f = GridF32::new(8, 8, 0.2);
        let t = land_topology(&f, 0.5);
        assert_eq!(t.num_landmasses, 0);
        assert_eq!(t.largest_cells, 0);
        assert_eq!(t.emerged_fraction, 0.0);
    }
}
