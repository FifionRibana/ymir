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

// ── Island-continent acceptance (M1 #190 geometric budget) ──────────────────

/// Acceptance criteria for "a continent surrounded by ocean, framed by the
/// export window". A seed passes if its LARGEST landmass does not wrap the torus
/// and fits inside `window_km` with an ocean margin on every side.
#[derive(Debug, Clone, Copy)]
pub struct IslandCriteria {
    /// Export window side (km) the mass must fit inside.
    pub window_km: f32,
    /// Ocean margin (km) required on EACH side (so the coast never touches the
    /// window edge → every basin terminates inside the window, no halo needed).
    pub ocean_margin_km: f32,
    /// Minimum acceptable traverse (km) — reject a speck of land.
    pub min_traverse_km: f32,
}

impl IslandCriteria {
    /// Max landmass traverse (km) that still leaves `ocean_margin_km` per side.
    pub fn max_traverse_km(&self) -> f32 {
        self.window_km - 2.0 * self.ocean_margin_km
    }
}

/// Does this land topology describe an island continent that fits the window?
/// Exactly one criterion per clause: at least one landmass, the largest does not
/// wrap either seam, and its bounding-box traverse fits the window with margin.
pub fn is_island_fit(t: &LandTopology, c: &IslandCriteria) -> bool {
    if t.num_landmasses == 0 || t.wraps_x || t.wraps_y {
        return false;
    }
    let traverse = t.bbox_km.0.max(t.bbox_km.1);
    traverse >= c.min_traverse_km && traverse <= c.max_traverse_km()
}

/// Border-clean island evaluation (M1 #190 reframe): the window is SIZED to the
/// largest landmass (`window_km = traverse + 2·margin`) and centred on its bbox
/// centre. The binding property is that the continent is surrounded by ocean —
/// so tentacular shape is fine; what matters is that NO land (of any mass) lies
/// on the window border ring. Compactness is reported, never selected on.
#[derive(Debug, Clone, Copy)]
pub struct IslandEval {
    /// Largest-mass metrics (from [`land_topology`]).
    pub topo: LandTopology,
    /// Window side (km) = `traverse + 2·margin`.
    pub window_km: f32,
    /// Resulting cell size at HD 8192² (m).
    pub m_per_cell: f32,
    /// `m_per_cell` within the 30–50 m band.
    pub resolution_ok: bool,
    /// No land cell (of ANY landmass) lies on the window border ring.
    pub border_clean: bool,
    /// Number of OTHER landmasses fully inside the window (satellite islets — a
    /// feature: independently placeable continents from one tectonic run).
    pub satellites_inside: usize,
    /// Border ring thickness in km (`ring_cells · km_per_cell`).
    pub ring_km: f32,
    /// Window centre (largest-mass bbox centre) in coarse cells.
    pub center_cell: (usize, usize),
    /// Reported only, never a selector: equiv-disc-diameter / traverse.
    pub compactness: f32,
}

impl IslandEval {
    /// Accept: largest mass does not wrap, resolution in band, border clean.
    pub fn accepted(&self) -> bool {
        !self.topo.wraps_x && !self.topo.wraps_y && self.resolution_ok && self.border_clean
    }
}

/// Evaluate the border-clean island predicate on a coarse normalized field, with
/// the window sized to the largest landmass and centred on its bbox centre. The
/// border ring is `ring_cells` thick on the (periodic) coarse grid. Land = value
/// `> sea`.
pub fn evaluate_island(coarse: &GridF32, sea: f32, margin_km: f32, ring_cells: usize) -> IslandEval {
    let topo = land_topology(coarse, sea);
    let (w, h) = (coarse.width, coarse.height);
    let km_per_cell = C1_DOMAIN_KM / w as f32;
    let traverse_km = topo.bbox_km.0.max(topo.bbox_km.1);
    let window_km = traverse_km + 2.0 * margin_km;
    let m_per_cell = window_km / 8192.0 * 1000.0;
    let resolution_ok = (30.0..=50.0).contains(&m_per_cell);
    let ring_km = ring_cells as f32 * km_per_cell;
    let disc = 2.0 * (topo.largest_area_km2 / std::f32::consts::PI).sqrt();
    let compactness = if traverse_km > 0.0 { disc / traverse_km } else { 0.0 };
    let cx = (topo.bbox_min.0 + topo.bbox_max.0) as f32 / 2.0;
    let cy = (topo.bbox_min.1 + topo.bbox_max.1) as f32 / 2.0;
    let center_cell = ((cx.round() as usize) % w, (cy.round() as usize) % h);

    let mut eval = IslandEval {
        topo,
        window_km,
        m_per_cell,
        resolution_ok,
        border_clean: false,
        satellites_inside: 0,
        ring_km,
        center_cell,
        compactness,
    };
    if topo.num_landmasses == 0 {
        return eval; // no land → not a clean island
    }

    // Periodic signed offset of integer coord `p` from float centre `c` (axis `len`).
    let offset = |p: usize, c: f32, len: usize| -> f32 {
        let l = len as f32;
        let mut d = p as f32 - c;
        d -= (d / l).round() * l;
        d
    };

    // CCL to label components (satellites need per-component membership).
    let n = w * h;
    let is_land = |k: usize| coarse.data[k] > sea;
    let mut ds = DisjointSet::new(n);
    for y in 0..h {
        for x in 0..w {
            let k = y * w + x;
            if !is_land(k) {
                continue;
            }
            let kr = y * w + (x + 1) % w;
            if is_land(kr) {
                ds.union(k, kr);
            }
            let kd = ((y + 1) % h) * w + x;
            if is_land(kd) {
                ds.union(k, kd);
            }
        }
    }
    let mut sizes: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for k in 0..n {
        if is_land(k) {
            *sizes.entry(ds.find(k)).or_insert(0) += 1;
        }
    }
    let largest_root =
        sizes.iter().max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0))).map(|(&r, _)| r).unwrap();

    let half = (window_km / km_per_cell).round().max(1.0) / 2.0;
    let ring = ring_cells as f32;
    let mut border_clean = true;
    // Per-component: stays strictly inside the window (no cell on/beyond the ring)?
    let mut inside: std::collections::HashMap<usize, bool> = std::collections::HashMap::new();
    for y in 0..h {
        for x in 0..w {
            let k = y * w + x;
            if !is_land(k) {
                continue;
            }
            let dx = offset(x, cx, w).abs();
            let dy = offset(y, cy, h).abs();
            let in_window = dx <= half && dy <= half;
            let on_ring = in_window && (dx > half - ring || dy > half - ring);
            if on_ring {
                border_clean = false;
            }
            let strictly_inside = in_window && !on_ring;
            let e = inside.entry(ds.find(k)).or_insert(true);
            if !strictly_inside {
                *e = false;
            }
        }
    }
    eval.border_clean = border_clean;
    eval.satellites_inside =
        inside.iter().filter(|&(&r, &ins)| r != largest_root && ins).count();
    eval
}

/// One landmass evaluated as a placeable continent with its OWN window.
#[derive(Debug, Clone, Copy)]
pub struct IslandCandidate {
    pub area_km2: f32,
    pub traverse_km: f32,
    pub window_km: f32,
    pub m_per_cell: f32,
    pub resolution_ok: bool,
    pub border_clean: bool,
    pub wraps: bool,
}

impl IslandCandidate {
    /// Placeable on a globe: non-wrapping, resolution in band, border clean.
    pub fn placeable(&self) -> bool {
        !self.wraps && self.resolution_ok && self.border_clean
    }
}

/// Evaluate EVERY landmass as its own border-clean continent (M1 #190 multi-
/// continent harvest, report-only): each mass gets a window sized to itself and
/// is placeable when non-wrapping, resolution in band, and no OTHER land lies on
/// its ring. One tectonic run can thus yield several placeable continents that
/// share a consistent tectonic history. Sorted by area (desc), deterministic.
pub fn harvest_islands(
    coarse: &GridF32,
    sea: f32,
    margin_km: f32,
    ring_cells: usize,
) -> Vec<IslandCandidate> {
    use std::collections::{HashMap, HashSet};
    let (w, h) = (coarse.width, coarse.height);
    let n = w * h;
    let km_per_cell = C1_DOMAIN_KM / w as f32;
    let cell_km2 = km_per_cell * km_per_cell;
    let is_land = |k: usize| coarse.data[k] > sea;

    let mut ds = DisjointSet::new(n);
    for y in 0..h {
        for x in 0..w {
            let k = y * w + x;
            if !is_land(k) {
                continue;
            }
            let kr = y * w + (x + 1) % w;
            if is_land(kr) {
                ds.union(k, kr);
            }
            let kd = ((y + 1) % h) * w + x;
            if is_land(kd) {
                ds.union(k, kd);
            }
        }
    }
    // Per-root area + bbox.
    let mut agg: HashMap<usize, (usize, usize, usize, usize, usize)> = HashMap::new();
    for y in 0..h {
        for x in 0..w {
            let k = y * w + x;
            if !is_land(k) {
                continue;
            }
            let r = ds.find(k);
            let e = agg.entry(r).or_insert((0, x, x, y, y));
            e.0 += 1;
            e.1 = e.1.min(x);
            e.2 = e.2.max(x);
            e.3 = e.3.min(y);
            e.4 = e.4.max(y);
        }
    }
    // Roots that wrap a seam (not placeable).
    let mut wraps: HashSet<usize> = HashSet::new();
    for y in 0..h {
        let (a, b) = (y * w + (w - 1), y * w);
        if is_land(a) && is_land(b) {
            wraps.insert(ds.find(a));
        }
    }
    for x in 0..w {
        let (a, b) = ((h - 1) * w + x, x);
        if is_land(a) && is_land(b) {
            wraps.insert(ds.find(a));
        }
    }

    let offset = |p: usize, c: f32, len: usize| -> f32 {
        let l = len as f32;
        let mut d = p as f32 - c;
        d -= (d / l).round() * l;
        d
    };
    let ring = ring_cells as f32;
    let mut roots: Vec<usize> = agg.keys().copied().collect();
    roots.sort_unstable(); // deterministic scan order
    let mut out = Vec::with_capacity(roots.len());
    for r in roots {
        let (area, xmin, xmax, ymin, ymax) = agg[&r];
        let traverse = ((xmax - xmin + 1).max(ymax - ymin + 1)) as f32 * km_per_cell;
        let window_km = traverse + 2.0 * margin_km;
        let m_per_cell = window_km / 8192.0 * 1000.0;
        let resolution_ok = (30.0..=50.0).contains(&m_per_cell);
        let wr = wraps.contains(&r);
        let (cx, cy) = ((xmin + xmax) as f32 / 2.0, (ymin + ymax) as f32 / 2.0);
        let half = (window_km / km_per_cell).round().max(1.0) / 2.0;
        // Border-clean: no land (any mass) on this window's ring.
        let mut border_clean = !wr;
        if border_clean {
            'scan: for y in 0..h {
                for x in 0..w {
                    if !is_land(y * w + x) {
                        continue;
                    }
                    let (dx, dy) = (offset(x, cx, w).abs(), offset(y, cy, h).abs());
                    if dx <= half && dy <= half && (dx > half - ring || dy > half - ring) {
                        border_clean = false;
                        break 'scan;
                    }
                }
            }
        }
        out.push(IslandCandidate {
            area_km2: area as f32 * cell_km2,
            traverse_km: traverse,
            window_km,
            m_per_cell,
            resolution_ok,
            border_clean,
            wraps: wr,
        });
    }
    out.sort_by(|a, b| b.area_km2.partial_cmp(&a.area_km2).unwrap());
    out
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

    /// Border-clean: a centred island with ocean margin passes; a second mass
    /// intruding on the window ring breaks it.
    #[test]
    fn border_clean_predicate() {
        let (w, h) = (64usize, 64usize);
        // Centred 18-cell block on a 64² torus (16 km/cell) → window ~338 km
        // (~41 m/cell), block inside with margin → border-clean, 0 satellites.
        let mut d = vec![0.2f32; w * h];
        for y in 23..41 {
            for x in 23..41 {
                d[y * w + x] = 0.8;
            }
        }
        let e = evaluate_island(&GridF32::from_vec(w, h, d.clone()), 0.5, 25.0, 1);
        assert!((30.0..=50.0).contains(&e.m_per_cell), "cell {} m in band", e.m_per_cell);
        assert!(e.resolution_ok);
        assert!(e.border_clean, "isolated centred block is border-clean");
        assert_eq!(e.satellites_inside, 0);
        assert!(e.accepted());

        // A distinct landmass on the window ring breaks border-clean.
        let mut d2 = d;
        d2[32 * w + 42] = 0.8; // inside the window but on its ring
        let e2 = evaluate_island(&GridF32::from_vec(w, h, d2), 0.5, 25.0, 1);
        assert!(!e2.border_clean, "a mass touching the ring is not border-clean");
        assert!(!e2.accepted());
    }

    /// Island acceptance: only a non-wrapping mass that fits the window with an
    /// ocean margin passes.
    #[test]
    fn island_fit_predicate() {
        let crit =
            IslandCriteria { window_km: 328.0, ocean_margin_km: 25.0, min_traverse_km: 80.0 };
        // max traverse = 328 − 50 = 278 km.
        let base = LandTopology {
            num_landmasses: 3,
            wraps_x: false,
            wraps_y: false,
            bbox_km: (250.0, 240.0),
            ..LandTopology::empty()
        };
        assert!(is_island_fit(&base, &crit), "250 km mass fits a 278 km budget");

        // Too big (traverse 300 > 278).
        assert!(!is_island_fit(&LandTopology { bbox_km: (300.0, 100.0), ..base }, &crit));
        // Wrapping → rejected regardless of size.
        assert!(!is_island_fit(&LandTopology { wraps_x: true, ..base }, &crit));
        assert!(!is_island_fit(&LandTopology { wraps_y: true, ..base }, &crit));
        // Speck below the floor.
        assert!(!is_island_fit(&LandTopology { bbox_km: (40.0, 30.0), ..base }, &crit));
        // No land.
        assert!(!is_island_fit(&LandTopology::empty(), &crit));
    }
}
