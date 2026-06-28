//! D8 flow routing, flow accumulation, basin labeling, and river extraction.

use std::cmp::Ordering as CmpOrd;
use std::collections::BinaryHeap;

use serde::{Deserialize, Serialize};

use crate::grid::GridF32;
use crate::terrain::noise::SeededNoise;

// ── D8 direction encoding ───────────────────────────────────────────────

/// D8 direction: 0-7 for the 8 neighbors, 255 = no outlet.
///   7  0  1
///   6  x  2
///   5  4  3
pub const DIR_NONE: u8 = 255;
pub const D8_DX: [i32; 8] = [0, 1, 1, 1, 0, -1, -1, -1];
pub const D8_DY: [i32; 8] = [-1, -1, 0, 1, 1, 1, 0, -1];
pub const D8_DIST: [f32; 8] = [1.0, 1.414, 1.0, 1.414, 1.0, 1.414, 1.0, 1.414];

// ── Configuration ───────────────────────────────────────────────────────

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct FlowConfig {
    pub sea_level: f32,
    /// #drainage fix A — optional micro-relief restored on the pit-filled FLATS
    /// before routing. The priority-flood fills depressions to a PERFECTLY flat
    /// sill; Garbrecht-Martz then imposes a uniform distance-to-outlet gradient →
    /// D8 routes cardinal-straight (straight channels + 90° junctions). Real
    /// plains are never perfectly flat (old meanders, deposits). This restores a
    /// little coherent micro-relief on those flats so the routing WANDERS down it
    /// instead. It touches ONLY the routing surface — `filled` (and thus every
    /// lake level) keeps the exact original sill, so the hydrology is unchanged.
    /// `None` → no perturbation (byte-identical legacy routing).
    pub flat_perturbation: Option<FlatPerturbation>,
}

impl Default for FlowConfig {
    fn default() -> Self {
        Self { sea_level: 0.1, flat_perturbation: None }
    }
}

/// #drainage fix A — coherent value-noise added to the Garbrecht-Martz flat
/// gradient so the routing wanders on the flats instead of running
/// cardinal-straight down the uniform distance-to-outlet gradient.
///
/// `amplitude` is a FRACTION of the G-M descent step (the per-cell `tl` weight,
/// `fhmax+1`). It MUST stay `< 0.5`: at `< 0.5` the toward-outlet neighbour is
/// still strictly lower after the noise, so drainage is mathematically guaranteed
/// (no spurious pit, network stays connected). Within that bound the noise
/// dominates the lateral tie-break → the river meanders. `frequency` (cells⁻¹)
/// sets the wavelength.
///
/// FREQUENCY vs the LARGE-FLAT combs (#fix A suite). The residual diagonal combs
/// sit on the LARGE flats: there the noise must flip the D8 lateral choice OFTEN
/// to break a long straight run into a short staircase. Counter-intuitively a LOW
/// frequency makes it WORSE (a smooth large-scale tilt = a long straight descent,
/// measured: large-flat straightness 59 %→72 %); a HIGHER frequency (shorter
/// wavelength) flips the choice more often → long combs collapse into short
/// staircases (λ≈14→6: large-flat 59 %→48 %, 17+ runs 6 %→1 %). So `frequency` is
/// the lever, raised, not a low-frequency octave.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct FlatPerturbation {
    pub seed: u32,
    pub amplitude: f32,
    pub frequency: f64,
    pub octaves: usize,
}

impl Default for FlatPerturbation {
    fn default() -> Self {
        // amplitude 0.45 of the G-M step (just under the 0.5 no-pit bound);
        // frequency 0.17 (λ≈6) — high enough to collapse the large-flat combs.
        Self { seed: 0xF1A7_5EED, amplitude: 0.45, frequency: 0.17, octaves: 4 }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct RiverConfig {
    pub stream_threshold: f32,
    pub river_threshold: f32,
    pub major_river_threshold: f32,
}

impl Default for RiverConfig {
    fn default() -> Self {
        Self { stream_threshold: 500.0, river_threshold: 2000.0, major_river_threshold: 10000.0 }
    }
}

// ── Results ─────────────────────────────────────────────────────────────

/// Result of the heavy flow computation (Phase A).
pub struct FlowResult {
    pub filled: GridF32,
    pub direction: Vec<u8>,
    pub accumulation: GridF32,
    pub basins: Vec<u32>,
    pub num_basins: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RiverSegment {
    pub points: Vec<(u32, u32)>,
    pub strahler_order: u8,
    pub avg_flow: f32,
    pub max_flow: f32,
    pub basin_id: u32,
    pub upstream: Vec<usize>,
    pub downstream: Option<usize>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct RiverNetwork {
    pub segments: Vec<RiverSegment>,
}

// ── Priority queue cell ─────────────────────────────────────────────────

#[derive(Clone, Copy)]
struct PqCell {
    height: f32,
    idx: usize,
}

impl PartialEq for PqCell {
    fn eq(&self, other: &Self) -> bool {
        self.height == other.height && self.idx == other.idx
    }
}
impl Eq for PqCell {}

// Min-heap: reverse ordering so BinaryHeap pops the lowest.
impl PartialOrd for PqCell {
    fn partial_cmp(&self, other: &Self) -> Option<CmpOrd> {
        Some(self.cmp(other))
    }
}
impl Ord for PqCell {
    fn cmp(&self, other: &Self) -> CmpOrd {
        other
            .height
            .partial_cmp(&self.height)
            .unwrap_or(CmpOrd::Equal)
            .then_with(|| other.idx.cmp(&self.idx))
    }
}

// ── Phase A: Heavy computation ──────────────────────────────────────────

/// Compute flow directions, accumulation, and basins from a heightmap.
pub fn compute_flow(heightmap: &GridF32, config: &FlowConfig) -> FlowResult {
    let w = heightmap.width;
    let h = heightmap.height;
    let n = w * h;

    // Step 1: identify ocean
    let mut is_ocean = vec![false; n];
    for i in 0..n {
        is_ocean[i] = heightmap.data[i] <= config.sea_level;
    }

    // Step 2: priority flood pit filling (fills depressions to the EXACT sill —
    // no epsilon increment, so flats are truly flat for the Garbrecht-Martz step).
    let filled = pit_fill(heightmap, &is_ocean, w, h);

    // Step 2b: flat resolution (Garbrecht-Martz 1997). Replaces the old eps-fill
    // micro-gradient (which followed the flood TREE → parallel-bar / 45°-fan
    // artifacts on flat interiors, #155 diagnostic 2ec0348). Imposes a convergent
    // drainage gradient on the EXACT-equal flats (= the pit-filled depressions;
    // native FBM-textured plateaus are never exactly flat, so they keep their
    // real micro-gradient and are untouched). Returns a per-cell `flat_grad`
    // (f64, 0 on non-flat cells) used ONLY for routing — `filled` keeps the exact
    // sill so lake levels are unchanged.
    //
    // #drainage fix A — with `flat_perturbation`, a coherent value-noise term is
    // ADDED to `flat_grad` on the flats (bounded below the G-M descent step) so D8
    // wanders laterally instead of running cardinal-straight down the uniform
    // distance gradient. It cannot create a pit (the bound keeps the toward-outlet
    // neighbour strictly lower) and never touches `filled` → the hydrology (lakes,
    // water balance, endorheic basins) is unchanged; only the TRACÉ wanders.
    let flat_grad =
        resolve_flats(&filled, &is_ocean, config.flat_perturbation.as_ref(), w, h);

    // Step 3: D8 flow direction (steepest descent on `filled`; flat cells routed
    // down the `flat_grad` toward the outlet).
    let direction = compute_d8(&filled, &flat_grad, &is_ocean, w, h);

    // Step 4: flow accumulation (topological order by filled, then flat_grad).
    let accumulation = compute_accumulation(&filled, &flat_grad, &direction, &is_ocean, w, h);

    // Step 5: basin labeling
    let (basins, num_basins) = compute_basins(&direction, &is_ocean, w, h);

    FlowResult { filled, direction, accumulation, basins, num_basins }
}

fn pit_fill(heightmap: &GridF32, is_ocean: &[bool], w: usize, h: usize) -> GridF32 {
    let n = w * h;
    let mut filled = heightmap.clone();
    let mut processed = vec![false; n];
    let mut pq = BinaryHeap::new();

    // Seed: ocean cells and land cells adjacent to ocean
    for j in 0..h {
        for i in 0..w {
            let idx = j * w + i;
            if is_ocean[idx] {
                processed[idx] = true;
                // Check neighbors for land cells to seed
                for d in 0..8 {
                    let ni = ((i as i32 + D8_DX[d]) % w as i32 + w as i32) as usize % w;
                    let nj = ((j as i32 + D8_DY[d]) % h as i32 + h as i32) as usize % h;
                    let nidx = nj * w + ni;
                    if !is_ocean[nidx] && !processed[nidx] {
                        processed[nidx] = true;
                        pq.push(PqCell { height: heightmap.data[nidx], idx: nidx });
                    }
                }
            }
        }
    }

    // Flood fill
    while let Some(cell) = pq.pop() {
        let ci = cell.idx % w;
        let cj = cell.idx / w;

        for d in 0..8 {
            let ni = ((ci as i32 + D8_DX[d]) % w as i32 + w as i32) as usize % w;
            let nj = ((cj as i32 + D8_DY[d]) % h as i32 + h as i32) as usize % h;
            let nidx = nj * w + ni;

            if processed[nidx] || is_ocean[nidx] {
                continue;
            }

            // Fill to the exact sill (no epsilon): a depression becomes a TRUE
            // flat, which the Garbrecht-Martz step then drains. (Was `+ eps`, the
            // flood-tree micro-gradient that caused the flat-routing artifact.)
            let fill_h = heightmap.data[nidx].max(filled.data[cell.idx]);
            filled.data[nidx] = fill_h;
            processed[nidx] = true;
            pq.push(PqCell { height: fill_h, idx: nidx });
        }
    }

    filled
}

fn compute_d8(filled: &GridF32, flat_grad: &[f64], is_ocean: &[bool], w: usize, h: usize) -> Vec<u8> {
    let n = w * h;
    let mut direction = vec![DIR_NONE; n];

    for j in 0..h {
        for i in 0..w {
            let idx = j * w + i;
            if is_ocean[idx] {
                continue;
            }

            let my_h = filled.data[idx];
            let mut best_dir = DIR_NONE;
            let mut best_slope = 0.0f32;

            // Pass 1 — steepest descent on `filled` (real terrain + the sill exits
            // of flats). Unchanged from the original; slopes/real drainage untouched.
            for d in 0..8 {
                let ni = ((i as i32 + D8_DX[d]) % w as i32 + w as i32) as usize % w;
                let nj = ((j as i32 + D8_DY[d]) % h as i32 + h as i32) as usize % h;
                let nidx = nj * w + ni;

                let slope = (my_h - filled.data[nidx]) / D8_DIST[d];
                if slope > best_slope {
                    best_slope = slope;
                    best_dir = d as u8;
                }
            }

            // Pass 2 — flat resolution. Only reached when no strictly-lower
            // neighbour exists (a pit-filled flat interior). Route to the
            // EQUAL-elevation neighbour with the smallest `flat_grad` (down the
            // Garbrecht-Martz gradient toward the outlet). `flat_grad` is 0 on
            // non-flat cells, so the flat's spill cells (which DO drain via pass 1)
            // attract the flow. Guaranteed to descend (the toward-lower term
            // dominates `flat_grad`), so no spurious sink is introduced.
            if best_dir == DIR_NONE {
                let my_g = flat_grad[idx];
                let mut best_g = my_g;
                for d in 0..8 {
                    let ni = ((i as i32 + D8_DX[d]) % w as i32 + w as i32) as usize % w;
                    let nj = ((j as i32 + D8_DY[d]) % h as i32 + h as i32) as usize % h;
                    let nidx = nj * w + ni;
                    if !is_ocean[nidx] && filled.data[nidx] == my_h && flat_grad[nidx] < best_g {
                        best_g = flat_grad[nidx];
                        best_dir = d as u8;
                    }
                }
            }

            direction[idx] = best_dir;
        }
    }

    direction
}

/// Garbrecht-Martz 1997 flat resolution. Imposes a convergent drainage gradient
/// over the TRULY-flat regions (cells whose `filled` value exactly equals a
/// neighbour's and which have no strictly-lower neighbour = the pit-filled
/// depression interiors). Returns a per-cell `flat_grad` (0 on non-flat cells);
/// a flat cell drains to the equal-elevation neighbour with the smallest value.
///
/// The gradient combines two distance transforms over the flat:
/// - **toward-lower** `tl`: graph distance from the flat's spill cells (the
///   equal-elevation cells adjacent to lower terrain / ocean). 0 at the spill,
///   growing inward.
/// - **away-from-higher** `fh`: graph distance from the flat's high edges (cells
///   adjacent to HIGHER terrain — where inflow enters). 0 at the high edge.
///
/// `flat_grad = tl·(fhmax+1) + (fhmax − fh)`. The `tl` term DOMINATES (weight
/// `fhmax+1` exceeds the `(fhmax − fh)` range), so every flat cell has a
/// strictly-smaller-`flat_grad` neighbour toward a spill → guaranteed drainage,
/// no interior minimum. The `(fhmax − fh)` term breaks ties between equal-`tl`
/// neighbours toward the convergent (away-from-inflow) direction — this is what
/// dissolves the parallel-bar / fan pattern the old eps-fill produced.
///
/// Native FBM-textured plateaus are never EXACTLY flat (micro-relief gives a
/// strictly-lower neighbour), so they are not classified as flat here and keep
/// their real diffuse drainage — the fix touches only the pit-filled flats.
///
/// #drainage fix A — with `perturb`, a coherent value-noise term bounded to
/// `amplitude·(fhmax+1)` with `amplitude < 0.5` is ADDED to `flat_grad` on the
/// flats. Since a single `tl` step is `(fhmax+1)` and the noise difference
/// between two cells is `< 2·0.5·(fhmax+1) = (fhmax+1)`, the toward-spill
/// neighbour (one `tl` lower) stays strictly lower → drainage still guaranteed,
/// no interior minimum. But the noise (up to nearly half a `tl` step) overrides
/// the `(fhmax − fh)` tie-break and small `tl` differences → the route WANDERS
/// instead of converging cardinal-straight.
fn resolve_flats(
    filled: &GridF32,
    is_ocean: &[bool],
    perturb: Option<&FlatPerturbation>,
    w: usize,
    h: usize,
) -> Vec<f64> {
    use std::collections::VecDeque;
    let n = w * h;
    let f = &filled.data;
    let nb = |i: usize, j: usize, d: usize| -> usize {
        let ni = ((i as i32 + D8_DX[d]) % w as i32 + w as i32) as usize % w;
        let nj = ((j as i32 + D8_DY[d]) % h as i32 + h as i32) as usize % h;
        nj * w + ni
    };

    // 1. Classify flat cells: land, exact-equal to a neighbour, no lower neighbour.
    let mut needs = vec![false; n];
    for j in 0..h {
        for i in 0..w {
            let c = j * w + i;
            if is_ocean[c] {
                continue;
            }
            let (mut has_lower, mut has_equal) = (false, false);
            for d in 0..8 {
                let m = nb(i, j, d);
                if is_ocean[m] || f[m] < f[c] {
                    has_lower = true;
                } else if f[m] == f[c] {
                    has_equal = true;
                }
            }
            needs[c] = !has_lower && has_equal;
        }
    }

    // 2. toward-lower BFS — seeds: flat cells adjacent to a spill (ocean, lower,
    //    or an equal NON-flat cell that itself drains).
    let mut tl = vec![-1i32; n];
    let mut q = VecDeque::new();
    for j in 0..h {
        for i in 0..w {
            let c = j * w + i;
            if !needs[c] {
                continue;
            }
            let mut is_spill = false;
            for d in 0..8 {
                let m = nb(i, j, d);
                if is_ocean[m] || f[m] < f[c] || (f[m] == f[c] && !needs[m]) {
                    is_spill = true;
                }
            }
            if is_spill {
                tl[c] = 0;
                q.push_back(c);
            }
        }
    }
    while let Some(c) = q.pop_front() {
        let (ci, cj) = (c % w, c / w);
        for d in 0..8 {
            let m = nb(ci, cj, d);
            if needs[m] && f[m] == f[c] && tl[m] < 0 {
                tl[m] = tl[c] + 1;
                q.push_back(m);
            }
        }
    }

    // 3. away-from-higher BFS — seeds: flat cells adjacent to higher terrain.
    let mut fh = vec![-1i32; n];
    let mut q2 = VecDeque::new();
    for j in 0..h {
        for i in 0..w {
            let c = j * w + i;
            if !needs[c] {
                continue;
            }
            let mut is_high = false;
            for d in 0..8 {
                let m = nb(i, j, d);
                if !is_ocean[m] && f[m] > f[c] {
                    is_high = true;
                }
            }
            if is_high {
                fh[c] = 0;
                q2.push_back(c);
            }
        }
    }
    while let Some(c) = q2.pop_front() {
        let (ci, cj) = (c % w, c / w);
        for d in 0..8 {
            let m = nb(ci, cj, d);
            if needs[m] && f[m] == f[c] && fh[m] < 0 {
                fh[m] = fh[c] + 1;
                q2.push_back(m);
            }
        }
    }
    let fhmax = fh.iter().copied().max().unwrap_or(0).max(0) as f64;

    // 4. Combined gradient (tl dominant → guaranteed descent; fh breaks bars).
    //    With a perturbation, add coherent noise bounded to amplitude·(fhmax+1)
    //    (amplitude < 0.5 keeps the toward-spill neighbour strictly lower → still
    //    no interior minimum) so the lateral choice wanders → meandering rivers.
    let noise = perturb.map(|p| (SeededNoise::new(p.seed, p.octaves.max(1)), p));
    let mut flat_grad = vec![0.0f64; n];
    for c in 0..n {
        if needs[c] {
            let t = if tl[c] < 0 { 0 } else { tl[c] } as f64;
            let hh = if fh[c] < 0 { fhmax } else { fh[c] as f64 };
            let mut g = t * (fhmax + 1.0) + (fhmax - hh);
            if let Some((ng, p)) = &noise {
                let (i, j) = (c % w, c / w);
                let bound = (p.amplitude as f64).min(0.49) * (fhmax + 1.0);
                let v = ng.fbm(i as f64 * p.frequency, j as f64 * p.frequency, p.octaves, 2.0, 0.5);
                g += v * bound;
            }
            flat_grad[c] = g;
        }
    }
    flat_grad
}

fn compute_accumulation(
    filled: &GridF32,
    flat_grad: &[f64],
    direction: &[u8],
    is_ocean: &[bool],
    w: usize,
    h: usize,
) -> GridF32 {
    let n = w * h;

    // Sort land cells by decreasing height, then decreasing flat_grad so flat
    // cells are processed from the inflow side down to the outlet (correct
    // topological order on flats where `filled` ties).
    let mut land_cells: Vec<usize> = (0..n).filter(|&i| !is_ocean[i]).collect();
    land_cells.sort_unstable_by(|&a, &b| {
        filled.data[b]
            .partial_cmp(&filled.data[a])
            .unwrap_or(CmpOrd::Equal)
            .then(flat_grad[b].partial_cmp(&flat_grad[a]).unwrap_or(CmpOrd::Equal))
    });

    let mut acc = GridF32::new(w, h, 0.0);
    // Initialize land cells to 1.0
    for &idx in &land_cells {
        acc.data[idx] = 1.0;
    }

    // Propagate downstream (highest to lowest)
    for &idx in &land_cells {
        let dir = direction[idx];
        if dir == DIR_NONE {
            continue;
        }
        let i = idx % w;
        let j = idx / w;
        let ni = ((i as i32 + D8_DX[dir as usize]) % w as i32 + w as i32) as usize % w;
        let nj = ((j as i32 + D8_DY[dir as usize]) % h as i32 + h as i32) as usize % h;
        let nidx = nj * w + ni;
        acc.data[nidx] += acc.data[idx];
    }

    acc
}

fn compute_basins(direction: &[u8], is_ocean: &[bool], w: usize, h: usize) -> (Vec<u32>, u32) {
    let n = w * h;
    let mut basins = vec![0u32; n];
    let mut next_basin = 1u32;

    // Trace each land cell downstream to find its basin
    for start in 0..n {
        if is_ocean[start] || basins[start] != 0 {
            continue;
        }

        // Trace path to ocean or to an already-labeled cell
        let mut path = Vec::new();
        let mut cur = start;
        loop {
            if is_ocean[cur] {
                break;
            }
            if basins[cur] != 0 {
                break;
            }
            path.push(cur);
            let dir = direction[cur];
            if dir == DIR_NONE {
                break;
            }
            let ci = cur % w;
            let cj = cur / w;
            let ni = ((ci as i32 + D8_DX[dir as usize]) % w as i32 + w as i32) as usize % w;
            let nj = ((cj as i32 + D8_DY[dir as usize]) % h as i32 + h as i32) as usize % h;
            cur = nj * w + ni;
        }

        // Determine basin ID
        let basin_id = if !is_ocean[cur] && basins[cur] != 0 {
            basins[cur]
        } else {
            let id = next_basin;
            next_basin += 1;
            id
        };

        // Label entire path
        for &idx in &path {
            basins[idx] = basin_id;
        }
    }

    (basins, next_basin - 1)
}

// ── Phase B: Interactive river extraction ────────────────────────────────

/// Extract the river network from flow data. Fast enough for main thread.
pub fn extract_rivers(
    flow_result: &FlowResult,
    config: &RiverConfig,
    width: usize,
    height: usize,
) -> RiverNetwork {
    let n = width * height;
    let acc = &flow_result.accumulation;
    let dir = &flow_result.direction;
    let basins = &flow_result.basins;

    // Identify river cells
    let is_river: Vec<bool> = (0..n).map(|i| acc.data[i] >= config.stream_threshold).collect();

    // Count upstream river neighbors for each cell
    let mut upstream_count = vec![0u8; n];
    for idx in 0..n {
        if !is_river[idx] || dir[idx] == DIR_NONE {
            continue;
        }
        let i = idx % width;
        let j = idx / width;
        let d = dir[idx] as usize;
        let ni = ((i as i32 + D8_DX[d]) % width as i32 + width as i32) as usize % width;
        let nj = ((j as i32 + D8_DY[d]) % height as i32 + height as i32) as usize % height;
        let nidx = nj * width + ni;
        if is_river[nidx] {
            upstream_count[nidx] += 1;
        }
    }

    // Find headwater cells (river cells with 0 upstream river neighbors)
    let headwaters: Vec<usize> =
        (0..n).filter(|&i| is_river[i] && upstream_count[i] == 0).collect();

    // Trace segments from each headwater
    let mut segments = Vec::new();
    let mut cell_to_segment: Vec<Option<usize>> = vec![None; n];

    for &start in &headwaters {
        trace_segment(
            start,
            &is_river,
            &upstream_count,
            dir,
            acc,
            basins,
            width,
            height,
            &mut segments,
            &mut cell_to_segment,
        );
    }

    // Link upstream/downstream
    for seg_idx in 0..segments.len() {
        let last = segments[seg_idx].points.last().copied();
        if let Some((lx, ly)) = last {
            let lidx = ly as usize * width + lx as usize;
            let d = dir[lidx];
            if d != DIR_NONE {
                let ni = ((lx as i32 + D8_DX[d as usize]) % width as i32 + width as i32) as usize
                    % width;
                let nj = ((ly as i32 + D8_DY[d as usize]) % height as i32 + height as i32) as usize
                    % height;
                let nidx = nj * width + ni;
                if let Some(ds_seg) = cell_to_segment[nidx] {
                    if ds_seg != seg_idx {
                        segments[seg_idx].downstream = Some(ds_seg);
                        segments[ds_seg].upstream.push(seg_idx);
                    }
                }
            }
        }
    }

    // Strahler ordering (process from headwaters)
    compute_strahler(&mut segments);

    RiverNetwork { segments }
}

fn trace_segment(
    start: usize,
    is_river: &[bool],
    upstream_count: &[u8],
    dir: &[u8],
    acc: &GridF32,
    basins: &[u32],
    w: usize,
    h: usize,
    segments: &mut Vec<RiverSegment>,
    cell_to_segment: &mut [Option<usize>],
) {
    let seg_idx = segments.len();
    let mut points = Vec::new();
    let mut flow_sum = 0.0f32;
    let mut max_flow = 0.0f32;

    let mut cur = start;
    loop {
        let ci = cur % w;
        let cj = cur / w;
        points.push((ci as u32, cj as u32));
        cell_to_segment[cur] = Some(seg_idx);
        let f = acc.data[cur];
        flow_sum += f;
        max_flow = max_flow.max(f);

        let d = dir[cur];
        if d == DIR_NONE {
            break;
        }

        let ni = ((ci as i32 + D8_DX[d as usize]) % w as i32 + w as i32) as usize % w;
        let nj = ((cj as i32 + D8_DY[d as usize]) % h as i32 + h as i32) as usize % h;
        let nidx = nj * w + ni;

        if !is_river[nidx] {
            break;
        }

        // Stop at junctions (cells with >=2 upstream river neighbors)
        // unless this is the first cell of the segment
        if upstream_count[nidx] >= 2 && points.len() > 1 {
            // Include the junction cell as last point of this segment
            points.push((ni as u32, nj as u32));
            cell_to_segment[nidx] = Some(seg_idx);
            let f = acc.data[nidx];
            flow_sum += f;
            max_flow = max_flow.max(f);

            // Start a new segment from the junction if not already traced
            if cell_to_segment[nidx].is_none() || cell_to_segment[nidx] == Some(seg_idx) {
                let avg = flow_sum / points.len() as f32;
                let basin_id = basins[start];
                segments.push(RiverSegment {
                    points,
                    strahler_order: 1,
                    avg_flow: avg,
                    max_flow,
                    basin_id,
                    upstream: Vec::new(),
                    downstream: None,
                });

                // Continue tracing from junction as a new segment
                // only if this junction hasn't already started a segment
                let d2 = dir[nidx];
                if d2 != DIR_NONE {
                    let ni2 = ((ni as i32 + D8_DX[d2 as usize]) % w as i32 + w as i32) as usize % w;
                    let nj2 = ((nj as i32 + D8_DY[d2 as usize]) % h as i32 + h as i32) as usize % h;
                    let nidx2 = nj2 * w + ni2;
                    if is_river[nidx2] && cell_to_segment[nidx2].is_none() {
                        trace_segment(
                            nidx,
                            is_river,
                            upstream_count,
                            dir,
                            acc,
                            basins,
                            w,
                            h,
                            segments,
                            cell_to_segment,
                        );
                    }
                }
                return;
            }
            break;
        }

        // Cell already part of another segment — stop
        if cell_to_segment[nidx].is_some() {
            break;
        }

        cur = nidx;
    }

    if points.is_empty() {
        return;
    }

    let avg = flow_sum / points.len() as f32;
    let basin_id = basins[start];
    segments.push(RiverSegment {
        points,
        strahler_order: 1,
        avg_flow: avg,
        max_flow,
        basin_id,
        upstream: Vec::new(),
        downstream: None,
    });
}

fn compute_strahler(segments: &mut [RiverSegment]) {
    // Iterative: keep processing until no changes
    let mut changed = true;
    while changed {
        changed = false;
        for i in 0..segments.len() {
            if segments[i].upstream.is_empty() {
                if segments[i].strahler_order != 1 {
                    segments[i].strahler_order = 1;
                    changed = true;
                }
                continue;
            }

            let upstream_orders: Vec<u8> =
                segments[i].upstream.iter().map(|&u| segments[u].strahler_order).collect();

            let max_order = *upstream_orders.iter().max().unwrap_or(&1);
            let count_max = upstream_orders.iter().filter(|&&o| o == max_order).count();

            let new_order = if count_max >= 2 { max_order + 1 } else { max_order };

            if segments[i].strahler_order != new_order {
                segments[i].strahler_order = new_order;
                changed = true;
            }
        }
    }
}

// ── Tests ───────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cone_flow_accumulation() {
        let n = 64;
        let mut hmap = GridF32::new(n, n, 0.0);
        let center = n as f32 / 2.0;
        for j in 0..n {
            for i in 0..n {
                let dx = i as f32 - center;
                let dy = j as f32 - center;
                let dist = (dx * dx + dy * dy).sqrt();
                hmap.set(i, j, (1.0 - dist / (center * 0.9)).max(0.0));
            }
        }

        let config = FlowConfig { sea_level: 0.05, ..Default::default() };
        let result = compute_flow(&hmap, &config);

        // Max accumulation should be large (many cells drain to edges)
        let max_acc = result.accumulation.data.iter().cloned().fold(0.0f32, f32::max);
        assert!(max_acc > 10.0, "Should have significant flow accumulation: {max_acc}");
    }

    #[test]
    fn pit_filling_resolves_pit() {
        let mut hmap = GridF32::new(32, 32, 0.3);
        // Create a pit at (16,16) surrounded by higher terrain
        hmap.set(16, 16, 0.1);
        // Ocean border
        for i in 0..32 {
            hmap.set(i, 0, 0.0);
            hmap.set(i, 31, 0.0);
            hmap.set(0, i, 0.0);
            hmap.set(31, i, 0.0);
        }

        let config = FlowConfig { sea_level: 0.05, ..Default::default() };
        let result = compute_flow(&hmap, &config);

        // The pit cell should have a valid flow direction
        let pit_idx = 16 * 32 + 16;
        assert_ne!(result.direction[pit_idx], DIR_NONE, "Pit should have a flow direction");

        // Filled height should be >= surrounding minimum
        assert!(
            result.filled.data[pit_idx] >= 0.1,
            "Pit should be filled: {}",
            result.filled.data[pit_idx]
        );
    }

    #[test]
    fn tilted_plane_flow_increases_downhill() {
        let n = 32;
        let mut hmap = GridF32::new(n, n, 0.0);
        for j in 0..n {
            for i in 0..n {
                // High on left, ocean on right
                let t = i as f32 / n as f32;
                hmap.set(i, j, 0.5 - t * 0.6);
            }
        }

        let config = FlowConfig { sea_level: 0.05, ..Default::default() };
        let result = compute_flow(&hmap, &config);

        // Flow accumulation should generally increase toward the right (downhill)
        // Compare average accumulation of left quarter vs right quarter (land only)
        let quarter = n / 4;
        let left_avg: f32 = (0..n)
            .flat_map(|j| (0..quarter).map(move |i| (i, j)))
            .filter(|&(i, j)| hmap.data[j * n + i] > config.sea_level)
            .map(|(i, j)| result.accumulation.data[j * n + i])
            .sum::<f32>()
            / (n * quarter) as f32;

        let right_avg: f32 = (0..n)
            .flat_map(|j| (quarter * 2..quarter * 3).map(move |i| (i, j)))
            .filter(|&(i, j)| hmap.data[j * n + i] > config.sea_level)
            .map(|(i, j)| result.accumulation.data[j * n + i])
            .sum::<f32>()
            / (n * quarter) as f32;

        assert!(
            right_avg > left_avg,
            "Downhill should have more flow: left={left_avg}, right={right_avg}"
        );
    }

    #[test]
    fn flow_all_land_cells_have_direction() {
        let n = 64;
        let mut hmap = GridF32::new(n, n, 0.0);
        let center = n as f32 / 2.0;
        for j in 0..n {
            for i in 0..n {
                let dx = i as f32 - center;
                let dy = j as f32 - center;
                let dist = (dx * dx + dy * dy).sqrt();
                hmap.set(i, j, (1.0 - dist / (center * 0.9)).max(0.0));
            }
        }

        let config = FlowConfig { sea_level: 0.05, ..Default::default() };
        let result = compute_flow(&hmap, &config);

        // Every land cell must have a valid flow direction (pit filling guarantees this)
        let land_without_dir = (0..n * n)
            .filter(|&i| hmap.data[i] > config.sea_level && result.direction[i] == DIR_NONE)
            .count();
        assert_eq!(land_without_dir, 0, "All land cells should have a valid flow direction");

        // All land cells should have accumulation >= 1.0
        let min_land_acc = (0..n * n)
            .filter(|&i| hmap.data[i] > config.sea_level)
            .map(|i| result.accumulation.data[i])
            .fold(f32::INFINITY, f32::min);
        assert!(
            min_land_acc >= 1.0,
            "All land cells should have accumulation >= 1.0, got {min_land_acc}"
        );
    }
}
