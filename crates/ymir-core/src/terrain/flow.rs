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
    /// #drainage fix A — D∞ (Tarboton) on the FLATS: route along the CONTINUOUS
    /// angle of the `flat_grad` gradient (split fractionally between the two
    /// bracketing D8 neighbours) instead of snapping to one of 8 directions. The
    /// slope pass stays mono-D8. `direction` carries the PRIMARY (larger-fraction)
    /// neighbour, so lakes/basins/B use a valid mono direction; only
    /// `accumulation` becomes fractional. `false` → mono D8 (byte-identical).
    ///
    /// ⚠️ KEPT OFF — measured INFERIOR (`probe_dinf_compare`). The rendered river
    /// trace follows the primary (one D8 neighbour per cell), so it re-quantises;
    /// and the `flat_grad` distance field (8-conn BFS ≈ Chebyshev) has a CARDINAL
    /// gradient, so the continuous angle points cardinally → the trace collapses
    /// to CARDINAL combs (diagonal 67 %→1 %, local orientation R2 0.40→0.91 —
    /// WORSE than D8's diagonal). The grid-rendered parallelism is a floor of
    /// drawing rivers as grid cells; no flat-routing tweak removes it. Default D8.
    pub dinf: bool,
}

impl Default for FlowConfig {
    fn default() -> Self {
        Self { sea_level: 0.1, flat_perturbation: None, dinf: false }
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
    /// ADR 0001 Finding 37 — extend each RETAINED watercourse (one that reaches `stream_threshold`
    /// somewhere) UPSTREAM to this accumulation, in the SAME cell units as `stream_threshold`. The
    /// natural value is the erosion regime split's critical area A_c (where the fluvial regime, and
    /// thus the channel, begins). `<= 0` OR `>= stream_threshold` ⟹ NO extension (the first exported
    /// point stays at `stream_threshold`) — the byte-identical default. The watercourse COUNT is
    /// unchanged (retention still decides WHICH exist); only their upstream extent grows.
    pub head_threshold: f32,
    /// With extension on: `true` ramifies the WHOLE upstream tree down to `head_threshold` (dense
    /// dendritic network — filter by Strahler order at render time); `false` extends only the MAIN
    /// STEM (the max-accumulation branch at each confluence — one headwater tail per watercourse).
    pub full_tree: bool,
}

impl Default for RiverConfig {
    fn default() -> Self {
        Self {
            stream_threshold: 500.0,
            river_threshold: 2000.0,
            major_river_threshold: 10000.0,
            head_threshold: 0.0, // no upstream extension → byte-identical
            full_tree: true,
        }
    }
}

// ── Results ─────────────────────────────────────────────────────────────

/// Result of the heavy flow computation (Phase A).
#[derive(Debug, Clone)]
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
    let flat_grad = resolve_flats(&filled, &is_ocean, config.flat_perturbation.as_ref(), w, h);

    // Step 3+4: flow direction + accumulation. D8 (mono) by default; with `dinf`,
    // the flat pass routes on the CONTINUOUS flat_grad gradient (fractional split)
    // → `accumulation` is fractional, `direction` is the primary neighbour.
    let (direction, accumulation) = if config.dinf {
        let (dir, dir2, frac1) = compute_dinf(&filled, &flat_grad, &is_ocean, w, h);
        let acc =
            compute_accumulation_dinf(&filled, &flat_grad, &dir, &dir2, &frac1, &is_ocean, w, h);
        (dir, acc)
    } else {
        let dir = compute_d8(&filled, &flat_grad, &is_ocean, w, h);
        let acc = compute_accumulation(&filled, &flat_grad, &dir, &is_ocean, w, h);
        (dir, acc)
    };

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

fn compute_d8(
    filled: &GridF32,
    flat_grad: &[f64],
    is_ocean: &[bool],
    w: usize,
    h: usize,
) -> Vec<u8> {
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

/// #drainage fix A — D∞ (Tarboton 1997) flow on the FLATS. The slope pass is the
/// SAME mono-D8 steepest descent as `compute_d8` (slopes already route well). On a
/// flat (no strictly-lower `filled` neighbour) the flow follows the CONTINUOUS
/// angle of the `flat_grad` gradient and is split between the two BRACKETING D8
/// neighbours by angular proximity. Returns `(primary_dir, secondary_dir, frac1)`:
/// `frac1` of the flow goes to `primary_dir`, `1−frac1` to `secondary_dir`
/// (`DIR_NONE` if none). Both candidates are valid downstream (equal `filled`,
/// strictly-lower `flat_grad`) → drainage still guaranteed, no pit. The primary is
/// the larger fraction, so `direction` stays a valid mono field for lakes/basins/B.
fn compute_dinf(
    filled: &GridF32,
    flat_grad: &[f64],
    is_ocean: &[bool],
    w: usize,
    h: usize,
) -> (Vec<u8>, Vec<u8>, Vec<f32>) {
    use std::f64::consts::FRAC_PI_4;
    let n = w * h;
    let mut dir = vec![DIR_NONE; n];
    let mut dir2 = vec![DIR_NONE; n];
    let mut frac1 = vec![1.0f32; n];
    // D8 direction angles (atan2(dy, dx)).
    let mut ang = [0.0f64; 8];
    for d in 0..8 {
        ang[d] = (D8_DY[d] as f64).atan2(D8_DX[d] as f64);
    }
    let angdiff = |a: f64, b: f64| -> f64 {
        let mut x = (a - b).abs() % (2.0 * std::f64::consts::PI);
        if x > std::f64::consts::PI {
            x = 2.0 * std::f64::consts::PI - x;
        }
        x
    };

    for j in 0..h {
        for i in 0..w {
            let idx = j * w + i;
            if is_ocean[idx] {
                continue;
            }
            let my_h = filled.data[idx];
            // Pass 1 — mono steepest descent on `filled` (slopes). Unchanged.
            let mut best_dir = DIR_NONE;
            let mut best_slope = 0.0f32;
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
            if best_dir != DIR_NONE {
                dir[idx] = best_dir;
                continue;
            }

            // Pass 2 — D∞ on the flat. Valid downstream = equal `filled`, strictly
            // lower `flat_grad`. Gradient of `flat_grad` over equal-filled
            // neighbours; flow = −∇flat_grad (toward the outlet).
            let my_g = flat_grad[idx];
            let nbr = |d: usize| -> usize {
                let ni = ((i as i32 + D8_DX[d]) % w as i32 + w as i32) as usize % w;
                let nj = ((j as i32 + D8_DY[d]) % h as i32 + h as i32) as usize % h;
                nj * w + ni
            };
            let g_at = |d: usize| -> f64 {
                let m = nbr(d);
                if !is_ocean[m] && filled.data[m] == my_h { flat_grad[m] } else { my_g }
            };
            // central differences (E=2, W=6, S=4, N=0); y grows downward.
            let gx = (g_at(2) - g_at(6)) * 0.5;
            let gy = (g_at(4) - g_at(0)) * 0.5;
            let (fx, fy) = (-gx, -gy);

            // valid downstream candidates.
            let valid = |d: usize| -> bool {
                let m = nbr(d);
                !is_ocean[m] && filled.data[m] == my_h && flat_grad[m] < my_g
            };

            if fx == 0.0 && fy == 0.0 {
                // degenerate gradient → mono to the min-flat_grad valid neighbour.
                let mut bg = my_g;
                for d in 0..8 {
                    if valid(d) && flat_grad[nbr(d)] < bg {
                        bg = flat_grad[nbr(d)];
                        dir[idx] = d as u8;
                    }
                }
                continue;
            }
            let theta = fy.atan2(fx);
            // bracketing pair = valid dirs within 45° of θ, weighted linearly.
            let mut wd = [0.0f64; 8];
            let mut wsum = 0.0;
            for d in 0..8 {
                if valid(d) {
                    let ad = angdiff(theta, ang[d]);
                    if ad < FRAC_PI_4 {
                        let wt = FRAC_PI_4 - ad;
                        wd[d] = wt;
                        wsum += wt;
                    }
                }
            }
            if wsum <= 0.0 {
                // θ points away from every valid neighbour → mono closest valid.
                let mut best = f64::MAX;
                for d in 0..8 {
                    if valid(d) {
                        let ad = angdiff(theta, ang[d]);
                        if ad < best {
                            best = ad;
                            dir[idx] = d as u8;
                        }
                    }
                }
                continue;
            }
            // top two weights → primary + secondary.
            let (mut p1, mut p2) = (usize::MAX, usize::MAX);
            for d in 0..8 {
                if wd[d] > 0.0 {
                    if p1 == usize::MAX || wd[d] > wd[p1] {
                        p2 = p1;
                        p1 = d;
                    } else if p2 == usize::MAX || wd[d] > wd[p2] {
                        p2 = d;
                    }
                }
            }
            dir[idx] = p1 as u8;
            if p2 != usize::MAX {
                dir2[idx] = p2 as u8;
                frac1[idx] = (wd[p1] / (wd[p1] + wd[p2])) as f32;
            }
        }
    }
    (dir, dir2, frac1)
}

/// #drainage fix A — fractional flow accumulation for D∞. Each cell sends `frac1`
/// of its accumulation to `direction` and `1−frac1` to `dir2` (if any), in
/// topological order (filled desc, flat_grad desc). The split conserves water
/// (`frac1 + (1−frac1) = 1`), so the total runoff is preserved exactly.
fn compute_accumulation_dinf(
    filled: &GridF32,
    flat_grad: &[f64],
    direction: &[u8],
    dir2: &[u8],
    frac1: &[f32],
    is_ocean: &[bool],
    w: usize,
    h: usize,
) -> GridF32 {
    let n = w * h;
    let mut land_cells: Vec<usize> = (0..n).filter(|&i| !is_ocean[i]).collect();
    land_cells.sort_unstable_by(|&a, &b| {
        filled.data[b]
            .partial_cmp(&filled.data[a])
            .unwrap_or(CmpOrd::Equal)
            .then(flat_grad[b].partial_cmp(&flat_grad[a]).unwrap_or(CmpOrd::Equal))
    });
    let mut acc = GridF32::new(w, h, 0.0);
    for &idx in &land_cells {
        acc.data[idx] = 1.0;
    }
    let push = |acc: &mut GridF32, from: usize, d: u8, frac: f32| {
        if d == DIR_NONE || frac <= 0.0 {
            return;
        }
        let (i, j) = (from % w, from / w);
        let ni = ((i as i32 + D8_DX[d as usize]) % w as i32 + w as i32) as usize % w;
        let nj = ((j as i32 + D8_DY[d as usize]) % h as i32 + h as i32) as usize % h;
        acc.data[nj * w + ni] += acc.data[from] * frac;
    };
    for &idx in &land_cells {
        let d1 = direction[idx];
        if d1 == DIR_NONE {
            continue;
        }
        let f1 = frac1[idx];
        push(&mut acc, idx, d1, f1);
        push(&mut acc, idx, dir2[idx], 1.0 - f1);
    }
    acc
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

/// MULTIPLE-FLOW-DIRECTION accumulation (Freeman 1991 / Quinn 1991) on the pit-filled
/// surface, for the STREAM-POWER INCISION ONLY (rivers/lakes keep D8 — this reads
/// `filled`/`direction` from a D8 [`FlowResult`] and returns a separate accumulation
/// grid; the caller decides which consumer uses it). Each cell spreads its accumulation
/// to ALL lower neighbours weighted by `slopeᵖ / Σ slopeᵖ`. `p → ∞` recovers D8 (single
/// steepest); `p → small` disperses. Dispersing the drainage area breaks the positive
/// feedback (a rill captures area → incises → captures more) that drives the
/// Smith–Bretherton parallel-rilling comb (ADR 0001 Finding 10), so the comb never forms
/// — attacking the CAUSE rather than smoothing it away. Flat cells (no lower `filled`
/// neighbour) fall back to the single D8 `direction` (which the Garbrecht–Martz flat pass
/// already resolved). Deterministic (descending-`filled` order, index tiebreak).
pub fn mfd_accumulation(
    filled: &GridF32,
    direction: &[u8],
    sea_level: f32,
    p: f32,
    w: usize,
    h: usize,
) -> GridF32 {
    let n = w * h;
    let is_ocean: Vec<bool> = (0..n).map(|i| filled.data[i] <= sea_level).collect();
    let mut land: Vec<usize> = (0..n).filter(|&i| !is_ocean[i]).collect();
    land.sort_unstable_by(|&a, &b| {
        filled.data[b].partial_cmp(&filled.data[a]).unwrap_or(CmpOrd::Equal).then(b.cmp(&a))
    });
    let mut acc = GridF32::new(w, h, 0.0);
    for &i in &land {
        acc.data[i] = 1.0;
    }
    let nbr = |i: usize, j: usize, d: usize| -> usize {
        let ni = ((i as i32 + D8_DX[d]) % w as i32 + w as i32) as usize % w;
        let nj = ((j as i32 + D8_DY[d]) % h as i32 + h as i32) as usize % h;
        nj * w + ni
    };
    let (mut wj, mut nj_idx) = ([0.0f32; 8], [0usize; 8]);
    for &c in &land {
        let (ci, cj) = (c % w, c / w);
        let zc = filled.data[c];
        let (mut wsum, mut cnt) = (0.0f32, 0usize);
        for d in 0..8 {
            let m = nbr(ci, cj, d);
            let drop = zc - filled.data[m];
            if drop > 0.0 {
                let slope = drop / D8_DIST[d];
                let wgt = slope.powf(p);
                wj[cnt] = wgt;
                nj_idx[cnt] = m;
                wsum += wgt;
                cnt += 1;
            }
        }
        let flow = acc.data[c];
        if cnt == 0 || wsum <= 0.0 {
            // Flat / no lower neighbour → single D8 fallback (already resolved).
            let d = direction[c];
            if d != DIR_NONE {
                acc.data[nbr(ci, cj, d as usize)] += flow;
            }
            continue;
        }
        for k in 0..cnt {
            acc.data[nj_idx[k]] += flow * wj[k] / wsum;
        }
    }
    acc
}

/// DRAINAGE CARVE (breach) — guarantee the exported river long profile is MONOTONICALLY
/// NON-INCREASING toward the sea (water cannot climb; ADR 0001 Finding 12, DEFECT 2). The
/// uphill artefact came from rivers traced on the pit-FILLED surface crossing filled
/// depressions whose REAL floor dips then climbs to the sill. This lowers each downstream
/// receiver to at most its donor's carved height, in donor→receiver order (decreasing
/// `filled`, so a receiver is finalised only after all its donors) — O(n) after the sort,
/// deterministic. **Lake cells are SKIPPED**: a genuine closed depression is a LAKE (flat
/// surface), not a trench, so the carve leaves it at its level and the monotonicity test
/// tolerates the flat crossing. Returns the carved heightmap; ocean is untouched.
pub fn carve_monotone(
    height: &GridF32,
    filled: &GridF32,
    direction: &[u8],
    lake_map: &[u32],
    sea_level: f32,
    w: usize,
    h: usize,
) -> GridF32 {
    let n = w * h;
    let mut carved = height.clone();
    // 1. FILL lakes to their flat surface (the pit-fill sill = `filled`), so a river
    //    crossing a lake runs on flat water (real depression floor → surface): the
    //    long profile is level there, not a climb out of the hollow. A lake is water,
    //    not a trench — this is why we fill (raise) rather than carve (lower) them.
    for k in 0..n {
        if lake_map[k] != 0 {
            carved.data[k] = filled.data[k];
        }
    }
    // 2. CARVE the non-lake reversals: lower each downstream receiver to at most its
    //    donor, in donor→receiver order (decreasing `filled`), so the real profile is
    //    monotone non-increasing along the network between lakes.
    let mut land: Vec<usize> = (0..n).filter(|&i| height.data[i] > sea_level).collect();
    land.sort_unstable_by(|&a, &b| {
        filled.data[b].partial_cmp(&filled.data[a]).unwrap_or(CmpOrd::Equal).then(b.cmp(&a))
    });
    for &k in &land {
        let d = direction[k];
        if d == DIR_NONE {
            continue;
        }
        let (i, j) = (k % w, k / w);
        let ni = ((i as i32 + D8_DX[d as usize]) % w as i32 + w as i32) as usize % w;
        let nj = ((j as i32 + D8_DY[d as usize]) % h as i32 + h as i32) as usize % h;
        let r = nj * w + ni;
        // Lakes are already flat (step 1); don't carve into/out of them.
        if lake_map[k] != 0 || lake_map[r] != 0 {
            continue;
        }
        if carved.data[r] > carved.data[k] {
            carved.data[r] = carved.data[k]; // receiver never above its donor → monotone
        }
    }
    carved
}

/// PRIORITY-FLOOD COMPLETE BREACHING (Lindsay 2016), lakes excepted — the GUARANTEED
/// one-pass monotone conditioning for DEFECT 2 (ADR 0001 Finding 13). Detected lakes are
/// pre-filled to their flat sill (a lake is water, not a trench); every OTHER cell is given
/// a monotone-descending path to the ocean or a lake by CARVING a trench along its outlet
/// path (never filling). Result: no non-lake pit remains, so the exported river long profile
/// cannot climb (flat across lakes only). Deterministic (min-heap by height, index tiebreak).
///
/// Per pit, the trench is carved backwards along the priority-flood tree with a small
/// `EPS` descent; carving stops as soon as the path is already low enough, so it is
/// near-linear in practice (Lindsay's complete breaching).
pub fn breach_monotone(
    height: &GridF32,
    filled: &GridF32,
    lake_map: &[u32],
    sea_level: f32,
    w: usize,
    h: usize,
) -> GridF32 {
    breach_monotone_protected(height, filled, lake_map, sea_level, w, h, None)
}

/// [`breach_monotone`] with a `protect` mask: cells set `true` are kept at their
/// ORIGINAL height (neither carved nor filled) and treated as inert barriers, so a
/// legitimate closed depression survives the breach — used for ACTIVE volcanic
/// crater bowls (C-2), which must stay closed for the climate-dependent crater-lake
/// stage to fill them (the generic breach otherwise flattens or breaches every pit,
/// climate-independently, in a cache the crater lakes cannot live in). `None` →
/// exactly [`breach_monotone`].
pub fn breach_monotone_protected(
    height: &GridF32,
    filled: &GridF32,
    lake_map: &[u32],
    sea_level: f32,
    w: usize,
    h: usize,
    protect: Option<&[bool]>,
) -> GridF32 {
    let n = w * h;
    let prot = |k: usize| protect.is_some_and(|p| p[k]);
    let mut z = height.data.clone();
    // Lakes → flat sill surface (base level, never breached). Protected cells keep
    // their original bowl (not flattened).
    for k in 0..n {
        if lake_map[k] != 0 && !prot(k) {
            z[k] = filled.data[k];
        }
    }
    let is_base = |k: usize, z: &[f32]| z[k] <= sea_level || lake_map[k] != 0;
    // Descent per carved step: ~0.1 m in norm (norm_to_m ≈ 11300 → 1e-5 ≈ 0.11 m).
    const EPS: f32 = 1e-5;
    let mut visited = vec![false; n];
    let mut backlink = vec![usize::MAX; n];
    let mut heap: BinaryHeap<PqCell> = BinaryHeap::new();
    // Protected cells are inert: marked visited (never carved/filled), not sources.
    for k in 0..n {
        if prot(k) {
            visited[k] = true;
        }
    }
    for k in 0..n {
        if !visited[k] && is_base(k, &z) {
            visited[k] = true;
            heap.push(PqCell { height: z[k], idx: k });
        }
    }
    while let Some(c) = heap.pop() {
        let ci = c.idx;
        let (cx, cy) = (ci % w, ci / w);
        for d in 0..8 {
            let nx = ((cx as i32 + D8_DX[d]) % w as i32 + w as i32) as usize % w;
            let ny = ((cy as i32 + D8_DY[d]) % h as i32 + h as i32) as usize % h;
            let nb = ny * w + nx;
            if visited[nb] {
                continue;
            }
            visited[nb] = true;
            backlink[nb] = ci;
            if height.data[nb] < z[ci] {
                // `nb` is lower than its only outlet `ci` → a pit. Breach: carve the outlet
                // path (ci → backlink → … → base) down to a descending ramp under nb.
                let mut target = height.data[nb] - EPS;
                let mut cur = ci;
                while cur != usize::MAX && !is_base(cur, &z) && z[cur] > target {
                    z[cur] = target;
                    target -= EPS;
                    cur = backlink[cur];
                }
            }
            // `nb` keeps its real elevation (it now drains via the carved/real path).
            z[nb] = height.data[nb];
            heap.push(PqCell { height: z[nb], idx: nb });
        }
    }
    // FILL mop-up (priority-flood, lakes excepted): the breach carves the deep/long
    // depressions into drained channels but leaves a small residual of micro-pits it cannot
    // fully connect in one pass; a standard priority-flood FILL raises each remaining non-lake
    // cell to at least its outlet, GUARANTEEING a monotone-descending path by construction.
    // It cannot undo the breach (carved channels already sit at/above their outlet), so it
    // only touches the residual — mostly-carved + a negligible flat fill = guaranteed monotone.
    let mut visited2 = vec![false; n];
    let mut heap2: BinaryHeap<PqCell> = BinaryHeap::new();
    for k in 0..n {
        if prot(k) {
            visited2[k] = true; // inert: never filled
        }
    }
    for k in 0..n {
        if !visited2[k] && is_base(k, &z) {
            visited2[k] = true;
            heap2.push(PqCell { height: z[k], idx: k });
        }
    }
    while let Some(c) = heap2.pop() {
        let ci = c.idx;
        let (cx, cy) = (ci % w, ci / w);
        for d in 0..8 {
            let nx = ((cx as i32 + D8_DX[d]) % w as i32 + w as i32) as usize % w;
            let ny = ((cy as i32 + D8_DY[d]) % h as i32 + h as i32) as usize % h;
            let nb = ny * w + nx;
            if visited2[nb] {
                continue;
            }
            visited2[nb] = true;
            if z[nb] < z[ci] {
                z[nb] = z[ci]; // raise a residual pit to its outlet → monotone by construction
            }
            heap2.push(PqCell { height: z[nb], idx: nb });
        }
    }
    GridF32::from_vec(w, h, z)
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

    let stream = config.stream_threshold;
    // Downstream cell of `k` (toroidal D8), if any.
    let ds_idx = |k: usize| -> Option<usize> {
        let d = dir[k];
        if d == DIR_NONE {
            return None;
        }
        let (i, j) = (k % width, k / width);
        let ni = ((i as i32 + D8_DX[d as usize]) % width as i32 + width as i32) as usize % width;
        let nj = ((j as i32 + D8_DY[d as usize]) % height as i32 + height as i32) as usize % height;
        Some(nj * width + ni)
    };

    // Identify river cells. Finding 37: with `head_threshold` in (0, stream) the network extends
    // UPSTREAM from `stream` down to `head_threshold` for every watercourse that reaches `stream`.
    let head = if config.head_threshold > 0.0 { config.head_threshold.min(stream) } else { stream };
    let is_river: Vec<bool> = if head >= stream {
        // No extension → the original mask (byte-identical).
        (0..n).map(|i| acc.data[i] >= stream).collect()
    } else {
        let dense: Vec<bool> = (0..n).map(|i| acc.data[i] >= head).collect();
        // Dense cells in DECREASING accumulation → downstream processed before upstream.
        let mut order: Vec<usize> = (0..n).filter(|&i| dense[i]).collect();
        order.sort_by(|&a, &b| {
            acc.data[b].partial_cmp(&acc.data[a]).unwrap_or(std::cmp::Ordering::Equal)
        });
        // RETENTION: keep a dense cell only if its downstream reaches `stream` (so a sub-`stream`
        // gully draining to a small pit is NOT exported — the watercourse count stays stable).
        let mut retained = vec![false; n];
        for &k in &order {
            if acc.data[k] >= stream {
                retained[k] = true;
            } else if let Some(ds) = ds_idx(k) {
                if dense[ds] && retained[ds] {
                    retained[k] = true;
                }
            }
        }
        if config.full_tree {
            (0..n).map(|i| dense[i] && retained[i]).collect()
        } else {
            // MAIN STEM: keep every ≥ stream cell; in the extension, keep only the max-accumulation
            // upstream branch feeding each cell (one headwater tail per watercourse).
            let mut best_up = vec![usize::MAX; n];
            for &k in &order {
                if acc.data[k] >= stream || !retained[k] {
                    continue;
                }
                if let Some(ds) = ds_idx(k) {
                    if best_up[ds] == usize::MAX || acc.data[k] > acc.data[best_up[ds]] {
                        best_up[ds] = k;
                    }
                }
            }
            let mut ms: Vec<bool> =
                (0..n).map(|i| dense[i] && retained[i] && acc.data[i] >= stream).collect();
            for &k in &order {
                if acc.data[k] >= stream || !retained[k] {
                    continue;
                }
                if let Some(ds) = ds_idx(k) {
                    if ms[ds] && best_up[ds] == k {
                        ms[k] = true;
                    }
                }
            }
            ms
        }
    };

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

    /// PERMANENT acceptance guard (ADR 0001 Finding 13, DEFECT 2): after `breach_monotone`,
    /// NO non-lake land cell may be a strict local minimum — every one has a not-higher
    /// neighbour, i.e. a monotone-descending path to the sea exists, so no exported river
    /// can climb. A river that climbs is a bug, not a trade-off.
    #[test]
    fn breach_leaves_no_interior_pit() {
        let (w, h) = (48usize, 48usize);
        // A bowl: high rim, deep interior pit, with an ocean strip on the left edge.
        let mut d = vec![0.6f32; w * h];
        for j in 0..h {
            for i in 0..w {
                let (dx, dy) = (i as f32 - w as f32 / 2.0, j as f32 - h as f32 / 2.0);
                let r = (dx * dx + dy * dy).sqrt();
                d[j * w + i] = 0.35 + 0.02 * r; // rises outward → central depression
                if i == 0 {
                    d[j * w + i] = 0.0; // ocean outlet on the left edge
                }
            }
        }
        let height = GridF32::from_vec(w, h, d);
        let flow = compute_flow(&height, &FlowConfig { sea_level: 0.05, ..Default::default() });
        // No lakes for this test (empty mask) → everything must breach to the ocean.
        let lake_map = vec![0u32; w * h];
        let carved = breach_monotone(&height, &flow.filled, &lake_map, 0.05, w, h);
        // No interior (non-edge, non-ocean) cell may be a strict local minimum.
        let mut pits = 0;
        for j in 1..h - 1 {
            for i in 1..w - 1 {
                let k = j * w + i;
                if carved.data[k] <= 0.05 {
                    continue; // ocean
                }
                let mut has_not_higher = false;
                for dd in 0..8 {
                    let ni = ((i as i32 + D8_DX[dd]) % w as i32 + w as i32) as usize % w;
                    let nj = ((j as i32 + D8_DY[dd]) % h as i32 + h as i32) as usize % h;
                    if carved.data[nj * w + ni] <= carved.data[k] {
                        has_not_higher = true;
                        break;
                    }
                }
                if !has_not_higher {
                    pits += 1;
                }
            }
        }
        assert_eq!(pits, 0, "breach must leave no interior pit (found {pits})");
    }
}
