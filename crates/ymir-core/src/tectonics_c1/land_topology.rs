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
    /// TRUE BAND in x: the largest mass circumnavigates the torus (there is NO
    /// empty column) — a real east-west band with no coast. This is the ONLY case
    /// that should reject a mass as "not an island". A finite mass that merely
    /// STRADDLES the seam is NOT a band (see `straddles_x`).
    pub wraps_x: bool,
    /// TRUE BAND in y (no empty row).
    pub wraps_y: bool,
    /// The largest mass touches both x edges but is FINITE (seam-straddle) — a
    /// debug flag, NOT a rejection criterion. Such a mass is surrounded by ocean;
    /// it just renders split by an arbitrary origin.
    pub straddles_x: bool,
    pub straddles_y: bool,
    /// Bounding box of the largest mass in the ROLLED frame (unrolled so the mass
    /// is contiguous): `(min, max)` cells. `min` = the unroll origin.
    pub bbox_min: (usize, usize),
    pub bbox_max: (usize, usize),
    /// CIRCULAR extent `(width_km, height_km)` of the largest mass — correct
    /// across the seam (`domain − largest empty run`), NOT the wrapped bbox.
    pub bbox_km: (f32, f32),
    /// Torus centre (cells) of the largest mass — the unrolled mid-point, valid
    /// for a seam-straddling mass (periodic).
    pub center_cell: (usize, usize),
    /// Cyclic shift (cells) that makes the largest mass contiguous; compose with
    /// `window_offset_in_torus` to locate the exported window in the torus.
    pub roll_origin: (usize, usize),
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
            straddles_x: false,
            straddles_y: false,
            bbox_min: (0, 0),
            bbox_max: (0, 0),
            bbox_km: (0.0, 0.0),
            center_cell: (0, 0),
            roll_origin: (0, 0),
            emerged_fraction: 0.0,
        }
    }
}

/// Circular extent of a 1-D occupancy (a set of occupied slots on a ring of size
/// `n`). Returns `(band, extent, unroll)`:
/// - `band`: every slot occupied → the set circumnavigates (a true band);
/// - `extent`: slots in the TIGHTEST circular arc covering the set
///   (`n − largest empty run`);
/// - `unroll`: index of the slot just AFTER the largest empty run — the start of
///   the arc, i.e. the cyclic shift that makes the set contiguous.
fn circular_extent(occ: &[bool]) -> (bool, usize, usize) {
    let n = occ.len();
    let count = occ.iter().filter(|&&b| b).count();
    if count == 0 {
        return (false, 0, 0);
    }
    if count == n {
        return (true, n, 0); // circumnavigates → band
    }
    // Largest run of EMPTY slots, allowing the run to wrap the seam (scan 2n).
    let (mut best_gap, mut best_end) = (0usize, 0usize);
    let mut cur = 0usize;
    for i in 0..2 * n {
        if !occ[i % n] {
            cur += 1;
            if cur > best_gap && cur <= n {
                best_gap = cur;
                best_end = (i + 1) % n; // slot after this empty run
            }
        } else {
            cur = 0;
        }
    }
    (false, n - best_gap, best_end)
}

/// Circular extent of a mass from its per-axis occupancy (which columns / rows it
/// occupies). Correct across the seam: a finite mass straddling the seam reports
/// its true extent + a torus centre, and only a circumnavigating mass is a band.
struct MassExtent {
    band_x: bool,
    band_y: bool,
    straddles_x: bool,
    straddles_y: bool,
    extent_x: usize,
    extent_y: usize,
    roll_x: usize,
    roll_y: usize,
    center_x: usize,
    center_y: usize,
}

fn mass_extent(occ_x: &[bool], occ_y: &[bool]) -> MassExtent {
    let (w, h) = (occ_x.len(), occ_y.len());
    let (band_x, extent_x, roll_x) = circular_extent(occ_x);
    let (band_y, extent_y, roll_y) = circular_extent(occ_y);
    MassExtent {
        band_x,
        band_y,
        straddles_x: !band_x && occ_x[0] && occ_x[w - 1],
        straddles_y: !band_y && occ_y[0] && occ_y[h - 1],
        extent_x,
        extent_y,
        roll_x,
        roll_y,
        // Torus centre = middle of the circular arc, unrolled back onto the torus.
        center_x: (roll_x + extent_x / 2) % w,
        center_y: (roll_y + extent_y / 2) % h,
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

    // CIRCULAR extent of the largest mass (correct across the seam). Occupancy
    // per axis → largest empty run → finite extent + unroll origin + true-band.
    let mut occ_x = vec![false; w];
    let mut occ_y = vec![false; h];
    for y in 0..h {
        for x in 0..w {
            let k = y * w + x;
            if is_land(k) && ds.find(k) == largest_root {
                occ_x[x] = true;
                occ_y[y] = true;
            }
        }
    }
    let m = mass_extent(&occ_x, &occ_y);

    let cell_km2 = c1_cell_area_km2(w);
    let km_per_cell = C1_DOMAIN_KM / w as f32;
    let domain_km2 = C1_DOMAIN_KM * C1_DOMAIN_KM;
    LandTopology {
        num_landmasses,
        largest_cells,
        largest_area_km2: largest_cells as f32 * cell_km2,
        largest_area_frac: (largest_cells as f32 * cell_km2) / domain_km2,
        wraps_x: m.band_x,
        wraps_y: m.band_y,
        straddles_x: m.straddles_x,
        straddles_y: m.straddles_y,
        bbox_min: (m.roll_x, m.roll_y),
        bbox_max: (m.roll_x + m.extent_x.saturating_sub(1), m.roll_y + m.extent_y.saturating_sub(1)),
        bbox_km: (m.extent_x as f32 * km_per_cell, m.extent_y as f32 * km_per_cell),
        center_cell: (m.center_x, m.center_y),
        roll_origin: (m.roll_x, m.roll_y),
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
pub fn evaluate_island(
    coarse: &GridF32,
    sea: f32,
    margin_km: f32,
    ring_cells: usize,
) -> IslandEval {
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
    // Torus centre of the largest mass (correct across the seam; the periodic
    // `offset` below measures distances the short way round).
    let center_cell = topo.center_cell;
    let (cx, cy) = (center_cell.0 as f32, center_cell.1 as f32);

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
    eval.satellites_inside = inside.iter().filter(|&(&r, &ins)| r != largest_root && ins).count();
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
    use std::collections::HashMap;
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
    // Per-root: area + occupancy per axis (for the circular extent, seam-correct).
    let mut agg: HashMap<usize, (usize, Vec<bool>, Vec<bool>)> = HashMap::new();
    for y in 0..h {
        for x in 0..w {
            let k = y * w + x;
            if !is_land(k) {
                continue;
            }
            let r = ds.find(k);
            let e = agg.entry(r).or_insert_with(|| (0, vec![false; w], vec![false; h]));
            e.0 += 1;
            e.1[x] = true;
            e.2[y] = true;
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
        let (area, occ_x, occ_y) = &agg[&r];
        let m = mass_extent(occ_x, occ_y);
        let traverse = m.extent_x.max(m.extent_y) as f32 * km_per_cell;
        let window_km = traverse + 2.0 * margin_km;
        let m_per_cell = window_km / 8192.0 * 1000.0;
        let resolution_ok = (30.0..=50.0).contains(&m_per_cell);
        let wr = m.band_x || m.band_y; // only a TRUE band is unplaceable
        let (cx, cy) = (m.center_x as f32, m.center_y as f32);
        let area = *area;
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
    fn seam_straddle_x_is_finite_not_a_band() {
        // Land in columns 0..3 and 17..20 (rows 3..7): ONE finite mass straddling
        // the x seam (6 columns wide), NOT a circumnavigating band. The bug flagged
        // this as wraps_x; the fix reports straddles_x + a correct 6-column extent.
        let mut f = field_with_block(20, 10, 0, 3, 3, 7);
        for y in 3..7 {
            for x in 17..20 {
                f.data[y * 20 + x] = 0.8;
            }
        }
        let t = land_topology(&f, 0.5);
        assert!(!t.wraps_x, "a finite seam-straddle is NOT a band");
        assert!(t.straddles_x, "it does straddle the x seam");
        assert!(!t.wraps_y && !t.straddles_y);
        assert_eq!(t.num_landmasses, 1, "the two edge strips are one mass via the seam");
        // extent_x = 6 columns (17,18,19,0,1,2), not the whole width.
        assert!((t.bbox_km.0 - 6.0 * (1024.0 / 20.0)).abs() < 1e-3, "extent {} km", t.bbox_km.0);
    }

    #[test]
    fn full_band_x_is_a_band() {
        // Every column occupied (rows 4..6) → a real circumnavigating band.
        let f = field_with_block(20, 10, 0, 20, 4, 6);
        let t = land_topology(&f, 0.5);
        assert!(t.wraps_x, "a mass occupying every column IS a band");
    }

    #[test]
    fn seam_straddle_y_is_finite_not_a_band() {
        // Rows 0..2 and 8..10 land (cols 3..9): finite mass straddling the y seam.
        let mut f = field_with_block(12, 10, 3, 9, 0, 2);
        for y in 8..10 {
            for x in 3..9 {
                f.data[y * 12 + x] = 0.8;
            }
        }
        let t = land_topology(&f, 0.5);
        assert!(!t.wraps_y && t.straddles_y, "finite y-seam-straddle, not a band");
        assert!(!t.wraps_x);
    }

    /// STEP 2 regression guard: the SAME blob placed across the seam must report
    /// the same extent / traverse / area as placed in the middle of the domain.
    #[test]
    fn seam_straddle_matches_centred_blob() {
        let (w, h) = (32usize, 32usize);
        let block = |cx: i32, cy: i32| {
            let mut d = vec![0.2f32; w * h];
            for dy in -3..=3 {
                for dx in -3..=3 {
                    let x = (cx + dx).rem_euclid(w as i32) as usize;
                    let y = (cy + dy).rem_euclid(h as i32) as usize;
                    d[y * w + x] = 0.8;
                }
            }
            land_topology(&GridF32::from_vec(w, h, d), 0.5)
        };
        let mid = block(16, 16); // centred
        let seam = block(0, 0); // straddles both seams
        assert_eq!(mid.largest_cells, seam.largest_cells, "same area");
        assert_eq!(mid.bbox_km, seam.bbox_km, "same circular extent");
        assert!(!seam.wraps_x && !seam.wraps_y, "a 7×7 blob is never a band");
        assert!(seam.straddles_x && seam.straddles_y, "the seam blob straddles both seams");
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
