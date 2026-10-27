//! ADR Finding 182 -- **the sediment deposition at the physics level** (gated: behind the continent profile's flag, `None`
//! everywhere else, so production and every earlier bench are unchanged).
//!
//! The declaration is `docs/reports/relief_method/f182_deposition/f182_declared.md`. In short:
//!
//! * **Fluvial** (Yuan et al. 2019, JGR, eq. 13–15, p. 7): the deposit at cell i is (G / p̃) / Ñ_i × Σ_{j ∈ ups(i)} (the net
//!   erosion of j), Ñ_i the cells of i's catchment, i included (eq. 14), p̃ = 1 (a uniform precipitation). G = 1.55, the
//!   middle of the natural range « ~0 to 3.1 » (Guerit et al. 2019, *Geology* 47(9), p. 855). **The coupling is split
//!   explicitly** (the author's decision): the implicit incision unchanged, then ONE pass upstream → downstream with the
//!   step's eroded matter; a cell's deposit never exceeds the flux entering it. Yuan's iterative coupling (eq. 19–21) is
//!   the bench's reference only, through [`FluvialHook`]; it never runs in the cascade.
//! * **Lakes** (A2): a closed depression of the routing surface (priority flood, `filled > z`) takes no fluvial deposit;
//!   what reaches its outlet fills it by volume, from the bottom up, to its spill; the remainder continues downstream.
//! * **Marine** (Yuan et al. 2019, EPSL, eq. 5–6 and Table 1): one class, the silt, K_M = 400 m²/yr (the Ogooué inversion's
//!   300–500), surface porosity 0.63; the sand is omitted (no sourced value: a declared simplification). Explicit sub-steps
//!   Δt_sub ≤ Δx² / (4 K_M) at the real cell (F162-0: no conversion); no iteration, no linear system. The transport moves
//!   only the deposited sediment: the C1 bathymetry is the tectonics' and the marine diffusion never erodes it.
//! * **The mass** (A4): eroded (the incision + the hillslope's loss to the sea) = land + lakes + sea (solid) + exported,
//!   per step, in f64, on the torus (nothing leaves the map).

use crate::grid::GridF32;
use crate::terrain::flow::{D8_DIST, D8_DX, D8_DY, DIR_NONE, FlowConfig, compute_flow};

/// The marine part (Yuan et al. 2019, EPSL): the silt's transport coefficient and its surface porosity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Marine {
    /// K_M (m²/yr), applied unconverted on the real cell (F162-0).
    pub k_m2_yr: f32,
    /// φ₁,₀, the silt's surface porosity: a solid volume V deposits a thickness V / (1 − φ).
    pub porosity: f32,
}

impl Marine {
    /// The declaration: K_M = 400 m²/yr, the middle of the Ogooué inversion's best fit « in the range of 300−500 m²/yr »
    /// (abstract, p. 1; the sensitivity 300 and 500), φ = 0.63 (Table 1, p. 6).
    pub fn declared() -> Self {
        Self { k_m2_yr: 400.0, porosity: 0.63 }
    }
}

/// ADR Finding 182 -- the deposition at the physics level. `lakes: false` routes the lakes like any cell (no fill);
/// `marine: None` exports what reaches the sea (deleted, as before F182).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Deposition {
    /// G (dimensionless; p̃ = 1).
    pub g: f32,
    pub lakes: bool,
    pub marine: Option<Marine>,
    /// ADR Finding 183-A1 -- the sediments' load through the rebound's flexure (`None` = F182: the rebound reads the incision).
    pub load: Option<Load>,
    /// ADR Finding 183-A2 -- the compaction (`None` = F182: the land deposits without porosity, the marine at φ₀).
    pub compaction: Option<Compaction>,
    /// ADR Finding 183-A3 -- the recycling by subduction, read from C1's consumption (`false` = F182).
    pub recycle: bool,
    /// ADR Finding 184-A1 -- the marine domain re-read at each sub-step (`false` = F182 / F183: fixed for the step).
    pub marine_live: bool,
    /// ADR Finding 184-A2 -- the fluvial stage in this many equal sub-steps (1 = F183).
    pub fluvial_substeps: usize,
    /// ADR Finding 185 -- the fluvial sub-step bounded by a duration (years): a step of duration d runs ⌈d / Δt_f⌉ sub-steps;
    /// `f32::INFINITY` = F184 (`fluvial_substeps`).
    pub fluvial_dt_max: f32,
    /// ADR Finding 186-A1 -- the fluvial deposition as a closure: `Some(G)` makes the incision use K / (1 + G (Ñ − 1)/Ñ) and the
    /// pass deposit nothing fluvial (it still routes, fills the lakes and delivers to the sea); `None` = F185.
    /// **Not adopted (F187, the author: « oui on supprime k_eff »): double counting with an effective K.**
    pub fluvial_keff: Option<f32>,
    /// ADR Finding 188-B -- **the sediment's erodibility**: a cell with a cover (H > 0) erodes at K2's soft cover's K (K_couv,
    /// ρ_K times the basement's, Stock & Montgomery 1999) while the lowering stays in the cover. `false` = F187.
    pub soft_cover: bool,
    /// ADR Finding 186-B -- the deposition handed to the amplifications (the chain on the basement, the burial closure
    /// [`crate::cascade::amplify::Transmit`]); read by the chain's caller. `false` = the amplifications as today.
    pub transmit: bool,
    /// ADR Finding 186 (the author's fallback) -- the live marine diffusion on 2 × 2 blocks (128² at 256²); `false` = per cell.
    pub marine_blocks: bool,
}

/// ADR Finding 183-A1 -- the densities of the load (kg/m³).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Load {
    /// ρ_g, the silt's grains: Sclater & Christie 1980, JGR 85(B7), Table A1a p. 3732, « Shale » (f₀ 0.63, c 0.51 km⁻¹,
    /// ρ_sg 2.72 g/cm³), the row matching Yuan et al. 2019 EPSL's silt.
    pub rho_g: f32,
    /// ρ_w, sea water: Whitehead & Clift 2009, JGR, p. 7 §[24] (1 030 kg/m³).
    pub rho_w: f32,
}

impl Load {
    pub fn declared() -> Self {
        Self { rho_g: 2720.0, rho_w: 1030.0 }
    }
}

/// ADR Finding 183-A2 -- Athy's law (Yuan et al. 2019 EPSL eq. 7, p. 3; Table 1, p. 6): φ(z) = φ₀ e^(−z/z₀).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Compaction {
    pub phi0: f32,
    /// z₀ (m); `f32::INFINITY` = no compaction (the A5 control).
    pub z0_m: f32,
}

impl Compaction {
    /// φ₀ = 0.63, z₀ = 1 960 m (Yuan EPSL Table 1, after Sclater & Christie 1980 Table A1a: 1 / 0.51 km⁻¹).
    pub fn declared() -> Self {
        Self { phi0: 0.63, z0_m: 1960.0 }
    }
    /// T(H) (m): the solid thickness of a column of height H in Athy's equilibrium, H − φ₀ z₀ (1 − e^(−H/z₀)).
    pub fn solid(&self, h_m: f64) -> f64 {
        let (p, z0) = (self.phi0 as f64, self.z0_m as f64);
        if z0.is_infinite() { (1.0 - p) * h_m } else { h_m - p * z0 * (1.0 - (-h_m / z0).exp()) }
    }
}

/// ADR Finding 184-A1 -- **the marine step with its domain re-read at each sub-step** (h < `sea`, norm): the sources enter in
/// equal parts per sub-step; a part on a cell at or above the sea walks down the steepest descent (D8, periodic) to the sea,
/// filling each sea cell it meets up to its room (`sea` − h) and walking on with the rest; at a closed depression the rest stays
/// (stranded, returned in thickness). The diffusion runs between the sub-step's sea cells, each donor limited to its supply
/// (`sup`, thickness) and each receiver's inflows capped at its room. Returns the sub-steps and the stranded thickness.
#[allow(clippy::too_many_arguments)]
pub fn marine_step_live(z: &mut GridF32, src: &[f64], dt_yr: f64, m: &Marine, cell_m: f64, sup: &mut [f64], torus: bool, sea: f64) -> (usize, f64) {
    marine_step_live_sea(z, src, dt_yr, m, cell_m, sup, torus, sea, None, None, None)
}

/// ADR Finding 187-A1/A2 -- [`marine_step_live`] with the ocean by connectivity (`Some(torus origin)`): at each sub-step the
/// domain is the largest 4-connected component of h < h_sea (relabelled only when a cell crossed the sea, so exactly the
/// labelling of every sub-step), the sources fill only the ocean's cells, and an unconnected below-sea cell is a closed
/// depression (a walk ending there strands its rest). `None` is [`marine_step_live`] exactly.
///
/// ADR Finding 189 2.3 -- with `exported_out`, **the walk circulates only in the sea, on a filled surface**: it moves on the
/// traversable sea (the 4-connected component of h ≤ `sea` holding the marine domain, the cells filled to the sea included,
/// fixed for the call: the step only raises cells and never above the sea), deposits on a marine-domain cell at most its room
/// and goes on down the steepest descent; with no lower traversable neighbour it goes to the nearest marine-domain cell with
/// room through the traversable cells (4-connected; the progradation over the filled cells). With no such cell the rest leaves
/// the map (added to `*exported_out`, thickness: the open window). A source on a cell out of the traversable sea is stranded
/// in place (none in the adopted profile: the cascade routes those by the fluvial pass first). `None` keeps
/// F184's walk.
#[allow(clippy::too_many_arguments)]
pub fn marine_step_live_sea(
    z: &mut GridF32,
    src: &[f64],
    dt_yr: f64,
    m: &Marine,
    cell_m: f64,
    sup: &mut [f64],
    torus: bool,
    sea: f64,
    connected: Option<Option<[usize; 2]>>,
    mut strand_out: Option<&mut [f64]>,
    mut exported_out: Option<&mut f64>,
) -> (usize, f64) {
    let (w, h) = (z.width, z.height);
    let n = w * h;
    let k_m = m.k_m2_yr as f64;
    let thick = 1.0 / (1.0 - m.porosity as f64);
    let srcs: Vec<usize> = (0..n).filter(|&k| src[k] > 0.0).collect();
    if srcs.is_empty() && !sup.iter().any(|&v| v > 0.0) {
        return (0, 0.0);
    }
    let nsub = if k_m > 0.0 && dt_yr > 0.0 { (dt_yr * 4.0 * k_m / (cell_m * cell_m)).ceil().max(1.0) as usize } else { 1 };
    let lambda = if k_m > 0.0 { k_m * (dt_yr / nsub as f64) / (cell_m * cell_m) } else { 0.0 };
    let at = |k: usize, dx: i32, dy: i32| -> usize {
        let (x, y) = ((k % w) as i32 + dx, (k / w) as i32 + dy);
        if torus {
            y.rem_euclid(h as i32) as usize * w + x.rem_euclid(w as i32) as usize
        } else if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 {
            usize::MAX
        } else {
            y as usize * w + x as usize
        }
    };
    let nb4: Vec<[usize; 4]> = (0..n).map(|k| [at(k, 1, 0), at(k, -1, 0), at(k, 0, 1), at(k, 0, -1)]).collect();
    let mut hz: Vec<f64> = z.data.iter().map(|&v| v as f64).collect();
    let mut stranded = 0f64;
    let (mut sc, mut cap, mut delta) = (vec![0f64; n], vec![1f64; n], vec![0f64; n]);
    // F187-A1: the ocean and the below-sea mask it was labelled from
    let (mut ocean, mut below): (Vec<bool>, Vec<bool>) = (Vec::new(), Vec::new());
    let refresh = |hz: &[f64], ocean: &mut Vec<bool>, below: &mut Vec<bool>| {
        if let Some(t) = connected {
            let b: Vec<bool> = hz.iter().map(|&v| v < sea).collect();
            if b != *below {
                *ocean = super::amplify::largest_component(&b, w, t);
                *below = b;
            }
        }
    };
    // F189 2.3: the traversable sea, and a stamp for the progradation's breadth-first search
    let filled_walk = exported_out.is_some();
    let mut trav: Vec<bool> = Vec::new();
    let (mut stamp, mut epoch, mut bfs) = (vec![0u32; n], 0u32, std::collections::VecDeque::new());
    if filled_walk {
        refresh(&hz, &mut ocean, &mut below);
        trav = sea_traversable(&hz, sea, if connected.is_some() { Some(&ocean) } else { None }, w, torus);
    }
    for _ in 0..nsub {
        refresh(&hz, &mut ocean, &mut below);
        if filled_walk {
            let in_sea = |c: usize, hz: &[f64], ocean: &[bool]| hz[c] < sea && (connected.is_none() || ocean[c]);
            for &k0 in &srcs {
                let mut rem = src[k0] * thick / nsub as f64;
                if !trav[k0] {
                    hz[k0] += rem;
                    sup[k0] += rem;
                    stranded += rem;
                    if let Some(so) = strand_out.as_deref_mut() {
                        so[k0] += rem;
                    }
                    continue;
                }
                let mut c = k0;
                for _ in 0..8 * n {
                    if in_sea(c, &hz, &ocean) {
                        let room = sea - hz[c];
                        if rem < room {
                            hz[c] += rem;
                            sup[c] += rem;
                            rem = 0.0;
                            break;
                        }
                        hz[c] = sea;
                        sup[c] += room;
                        rem -= room;
                        if rem <= 0.0 {
                            break;
                        }
                    }
                    // the steepest descent among the traversable neighbours
                    let mut best: Option<(usize, f64)> = None;
                    for d in 0..8 {
                        let j = at(c, D8_DX[d], D8_DY[d]);
                        if j == usize::MAX || !trav[j] {
                            continue;
                        }
                        let s = (hz[c] - hz[j]) / D8_DIST[d] as f64;
                        if s > 0.0 && best.is_none_or(|b| s > b.1) {
                            best = Some((j, s));
                        }
                    }
                    if let Some((j, _)) = best {
                        c = j;
                        continue;
                    }
                    // the progradation: the nearest marine cell with room through the traversable cells
                    epoch += 1;
                    bfs.clear();
                    bfs.push_back(c);
                    stamp[c] = epoch;
                    let mut found = None;
                    while let Some(q) = bfs.pop_front() {
                        if in_sea(q, &hz, &ocean) {
                            found = Some(q);
                            break;
                        }
                        for &j in &nb4[q] {
                            if j != usize::MAX && trav[j] && stamp[j] != epoch {
                                stamp[j] = epoch;
                                bfs.push_back(j);
                            }
                        }
                    }
                    match found {
                        Some(j) => c = j,
                        None => break,
                    }
                }
                if rem > 0.0 {
                    // the ocean filled with no deeper cell left: out of the map
                    if let Some(e) = exported_out.as_deref_mut() {
                        *e += rem;
                    }
                }
            }
        }
        // the sources: placed under the sea, walking seaward from a cell at or above it (F184's walk)
        for &k0 in srcs.iter().filter(|_| !filled_walk) {
            let mut rem = src[k0] * thick / nsub as f64;
            let mut c = k0;
            for _ in 0..n {
                if hz[c] < sea && (connected.is_none() || ocean[c]) {
                    let put = rem.min(sea - hz[c]);
                    hz[c] += put;
                    sup[c] += put;
                    rem -= put;
                    if rem <= 0.0 {
                        break;
                    }
                }
                let mut best: Option<(usize, f64)> = None;
                for d in 0..8 {
                    let j = at(c, D8_DX[d], D8_DY[d]);
                    if j == usize::MAX {
                        continue;
                    }
                    let s = (hz[c] - hz[j]) / D8_DIST[d] as f64;
                    if s > 0.0 && best.is_none_or(|b| s > b.1) {
                        best = Some((j, s));
                    }
                }
                match best {
                    Some((j, _)) => c = j,
                    None => {
                        hz[c] += rem;
                        sup[c] += rem;
                        stranded += rem;
                        if let Some(so) = strand_out.as_deref_mut() {
                            so[c] += rem;
                        }
                        rem = 0.0;
                        break;
                    }
                }
            }
            if rem > 0.0 {
                hz[c] += rem;
                sup[c] += rem;
                stranded += rem;
                if let Some(so) = strand_out.as_deref_mut() {
                    so[c] += rem;
                }
            }
        }
        if lambda <= 0.0 {
            continue;
        }
        // the sub-step's domain, the donors' supply limit, the receivers' room
        let dom: Vec<bool> = if connected.is_some() {
            refresh(&hz, &mut ocean, &mut below);
            ocean.clone()
        } else {
            hz.iter().map(|&v| v < sea).collect()
        };
        for k in 0..n {
            sc[k] = 0.0;
            if !dom[k] {
                continue;
            }
            let mut out = 0.0;
            for &j in &nb4[k] {
                if j != usize::MAX && dom[j] && hz[k] > hz[j] {
                    out += lambda * (hz[k] - hz[j]);
                }
            }
            if out > 0.0 {
                sc[k] = (sup[k].max(0.0) / out).min(1.0);
            }
        }
        for k in 0..n {
            cap[k] = 1.0;
            if !dom[k] {
                continue;
            }
            let mut inflow = 0.0;
            for &j in &nb4[k] {
                if j != usize::MAX && dom[j] && hz[j] > hz[k] {
                    inflow += lambda * (hz[j] - hz[k]) * sc[j];
                }
            }
            if inflow > 0.0 {
                cap[k] = ((sea - hz[k]).max(0.0) / inflow).min(1.0);
            }
        }
        for k in 0..n {
            delta[k] = 0.0;
            if !dom[k] {
                continue;
            }
            let mut d = 0.0;
            for &j in &nb4[k] {
                if j == usize::MAX || !dom[j] {
                    continue;
                }
                if hz[j] > hz[k] {
                    d += lambda * (hz[j] - hz[k]) * sc[j] * cap[k];
                } else if hz[k] > hz[j] {
                    d -= lambda * (hz[k] - hz[j]) * sc[k] * cap[j];
                }
            }
            delta[k] = d;
        }
        for k in 0..n {
            if delta[k] != 0.0 {
                hz[k] += delta[k];
                sup[k] += delta[k];
            }
        }
    }
    for k in 0..n {
        z.data[k] = hz[k] as f32;
    }
    (nsub, stranded)
}

/// ADR Finding 189 2.3 -- the traversable sea: the cells at or below `sea` 4-connected to the marine domain (`ocean`, the
/// connected ocean; `None`: every cell at or below the sea), the cells filled up to the sea included.
pub fn sea_traversable(hz: &[f64], sea: f64, ocean: Option<&[bool]>, w: usize, torus: bool) -> Vec<bool> {
    let n = hz.len();
    let h = n / w;
    let at = |k: usize, dx: i32, dy: i32| -> usize {
        let (x, y) = ((k % w) as i32 + dx, (k / w) as i32 + dy);
        if torus {
            y.rem_euclid(h as i32) as usize * w + x.rem_euclid(w as i32) as usize
        } else if x < 0 || y < 0 || x >= w as i32 || y >= h as i32 {
            usize::MAX
        } else {
            y as usize * w + x as usize
        }
    };
    let nb4: Vec<[usize; 4]> = (0..n).map(|k| [at(k, 1, 0), at(k, -1, 0), at(k, 0, 1), at(k, 0, -1)]).collect();
    match ocean {
        None => hz.iter().map(|&v| v <= sea).collect(),
        Some(o) => {
            let mut t = o.to_vec();
            let mut st: Vec<usize> = (0..n).filter(|&k| o[k]).collect();
            while let Some(k) = st.pop() {
                for &j in &nb4[k] {
                    if j != usize::MAX && !t[j] && hz[j] <= sea {
                        t[j] = true;
                        st.push(j);
                    }
                }
            }
            t
        }
    }
}

/// ADR Finding 186 (the author's fallback) -- **the live marine step on 2 × 2 blocks**:
/// - the restriction: a block's source is the sum of its 4 cells' (as a thickness over the block, ÷ 4); its room under the sea
///   is the mean of its cells' rooms, so the block's height is h_sea − that room when it has some, else its cells' mean;
/// - the block's supply is the mean of its cells';
/// - [`marine_step_live`] runs on the blocks (Δx doubled, so the stability limit and the sub-steps at 128²);
/// - the prolongation (F186 step 2): the blocks' thickness change interpolated bilinearly on the cells, a gain capped by the
///   cell's room under the sea (the author's rule), a loss by the cell's sediment (its symmetric extension, not in the author's
///   rule; the loss per block in proportion to the sediment was measured worse, F186), each cell's excess shared among the
///   block's cells and its ring that still have room (or sediment). (The first prolongation shared a gain « in proportion to the
///   room »: it sent the deposit to the deepest cells and kept the shelf from forming, the reviewer's error, F186.)
///
/// The mass is exact: the last cell takes the remainder. Returns the sub-steps and the stranded thickness (per cell).
#[allow(clippy::too_many_arguments)]
pub fn marine_step_live_blocks(z: &mut GridF32, src: &[f64], dt_yr: f64, m: &Marine, cell_m: f64, sup: &mut [f64], torus: bool, sea: f64) -> (usize, f64) {
    let n = z.width;
    let nb = n / 2;
    let kids = |b: usize| -> [usize; 4] {
        let (bx, by) = (b % nb, b / nb);
        let (x, y) = (2 * bx, 2 * by);
        [y * n + x, y * n + x + 1, (y + 1) * n + x, (y + 1) * n + x + 1]
    };
    let mut zb = GridF32::new(nb, nb, 0.0);
    let mut srcb = vec![0f64; nb * nb];
    let mut supb = vec![0f64; nb * nb];
    for b in 0..nb * nb {
        let c = kids(b);
        let room: f64 = c.iter().map(|&k| (sea - z.data[k] as f64).max(0.0)).sum::<f64>() / 4.0;
        let mean: f64 = c.iter().map(|&k| z.data[k] as f64).sum::<f64>() / 4.0;
        zb.data[b] = if room > 0.0 { (sea - room) as f32 } else { mean as f32 };
        srcb[b] = c.iter().map(|&k| src[k]).sum::<f64>() / 4.0;
        supb[b] = c.iter().map(|&k| sup[k]).sum::<f64>() / 4.0;
    }
    let z0b = zb.clone();
    let sup0b = supb.clone();
    let (ns, strd) = marine_step_live(&mut zb, &srcb, dt_yr, m, 2.0 * cell_m, &mut supb, torus, sea);
    let _ = &z0b;
    // ADR Finding 186 step 2 -- the prolongation (the author's): the blocks' thickness change interpolated bilinearly on the
    // cells (the block centres at 2b + 1, so the weights 3/4 and 1/4; on the torus the interpolation keeps the total exactly),
    // a gain capped by the cell's room under the sea, a loss by the cell's sediment; each cell's excess goes to the cells of its
    // block and the ring around it that still have room (or sediment), in proportion; what no cell can take stays on the cell.
    let db: Vec<f64> = (0..nb * nb).map(|bk| supb[bk] - sup0b[bk]).collect();
    let wrap = |i: i64, m: usize| -> usize { i.rem_euclid(m as i64) as usize };
    let nc = n * n;
    let mut t = vec![0f64; nc];
    for y in 0..n {
        for x in 0..n {
            let (bx0, fx) = if x % 2 == 0 { (x as i64 / 2 - 1, 0.75) } else { (x as i64 / 2, 0.25) };
            let (by0, fy) = if y % 2 == 0 { (y as i64 / 2 - 1, 0.75) } else { (y as i64 / 2, 0.25) };
            let at = |bx: i64, by: i64| -> f64 {
                if torus {
                    db[wrap(by, nb) * nb + wrap(bx, nb)]
                } else {
                    let (cx, cy) = (bx.clamp(0, nb as i64 - 1) as usize, by.clamp(0, nb as i64 - 1) as usize);
                    db[cy * nb + cx]
                }
            };
            let v0 = at(bx0, by0) * (1.0 - fx) + at(bx0 + 1, by0) * fx;
            let v1 = at(bx0, by0 + 1) * (1.0 - fx) + at(bx0 + 1, by0 + 1) * fx;
            t[y * n + x] = v0 * (1.0 - fy) + v1 * fy;
        }
    }
    // the interpolation's total against the blocks' (exact up to rounding on the torus; the rounding residue goes to the cells
    // in proportion to |t| so the mass is kept to the last bit the f64 allows)
    let tot_b: f64 = db.iter().sum::<f64>() * 4.0;
    let tot_c: f64 = t.iter().sum();
    let abs_c: f64 = t.iter().map(|v| v.abs()).sum();
    if abs_c > 0.0 {
        let r = tot_b - tot_c;
        for v in t.iter_mut() {
            *v += r * v.abs() / abs_c;
        }
    }
    // the caps and the excess
    let room = |k: usize, z: &GridF32| (sea - z.data[k] as f64).max(0.0);
    let mut excess = vec![0f64; nc];
    let mut put = vec![0f64; nc];
    for k in 0..nc {
        let v = t[k];
        if v > 0.0 {
            let r = room(k, z);
            let p = v.min(r);
            put[k] = p;
            excess[k] = v - p;
        } else if v < 0.0 {
            let a = sup[k].max(0.0);
            let p = v.max(-a);
            put[k] = p;
            excess[k] = v - p;
        }
    }
    let ring = |k: usize| -> Vec<usize> {
        let (x, y) = ((k % n) as i64, (k / n) as i64);
        let (bx, by) = (x - x.rem_euclid(2), y - y.rem_euclid(2));
        let mut v = Vec::with_capacity(16);
        for dy in -1..3 {
            for dx in -1..3 {
                let (cx, cy) = (bx + dx, by + dy);
                if torus {
                    v.push(wrap(cy, n) * n + wrap(cx, n));
                } else if cx >= 0 && cy >= 0 && (cx as usize) < n && (cy as usize) < n {
                    v.push(cy as usize * n + cx as usize);
                }
            }
        }
        v
    };
    for k in 0..nc {
        let e = excess[k];
        if e == 0.0 {
            continue;
        }
        let rg = ring(k);
        // what each ring cell can still take: room left (a gain) or sediment left (a loss)
        let cap_of = |j: usize, put: &[f64]| -> f64 {
            if e > 0.0 { (room(j, z) - put[j].max(0.0)).max(0.0) } else { (sup[j].max(0.0) + put[j].min(0.0)).max(0.0) }
        };
        let caps: Vec<f64> = rg.iter().map(|&j| cap_of(j, &put)).collect();
        let cs: f64 = caps.iter().sum();
        let mag = e.abs();
        if cs <= 0.0 {
            put[k] += e;
            continue;
        }
        let take = mag.min(cs);
        let mut left = take;
        let last = caps.iter().rposition(|&c| c > 0.0).unwrap();
        for (i, &j) in rg.iter().enumerate() {
            if caps[i] <= 0.0 {
                continue;
            }
            let d = if i == last { left } else { take * caps[i] / cs };
            left -= d;
            put[j] += d * e.signum();
        }
        // beyond the ring's capacity: on the cell itself
        put[k] += (mag - take) * e.signum();
    }
    for k in 0..nc {
        if put[k] != 0.0 {
            sup[k] += put[k];
            z.data[k] = (z.data[k] as f64 + put[k]) as f32;
        }
    }
    (ns, strd * 4.0)
}

/// ADR Finding 183-A1/A2 -- the solid a lowering `ec` (≤ h) removes from a cover of solid `s` and thickness `h` (norm): with
/// the compaction, S − T(H − e_c) (the column kept in Athy's equilibrium); without, S e_c / H (proportional).
pub fn cover_removed(tab: Option<&CompactionTable>, s: f64, h: f64, ec: f64, n2m: f64) -> f64 {
    if ec <= 0.0 || h <= 0.0 {
        return 0.0;
    }
    match tab {
        Some(t) => (s - t.c.solid((h - ec) * n2m) / n2m).clamp(0.0, s.max(0.0)),
        None => s * (ec / h).min(1.0),
    }
}

/// ADR Finding 183-A2 -- T⁻¹ as a table computed once: H at uniform solid steps `ds_m`, filled by one monotone sweep of a
/// fine H grid; read by linear interpolation; beyond the table, the end slope 1 / T'(H_max). No iteration at run time.
#[derive(Clone, Debug)]
pub struct CompactionTable {
    pub c: Compaction,
    pub ds_m: f64,
    pub h: Vec<f64>,
    /// The worst round trip |T(T⁻¹(S)) − S| (m) over the table's mid-points, measured at construction.
    pub max_err_m: f64,
}

impl CompactionTable {
    pub fn new(c: Compaction) -> Self {
        let (ds, h_max) = (0.5f64, 30_000.0f64);
        let s_max = c.solid(h_max);
        let ns = (s_max / ds).floor() as usize + 1;
        // the fine H grid (0.05 m), swept once
        let dh = 0.05f64;
        let mut h = vec![0f64; ns];
        let (mut hk, mut tk) = (0f64, 0f64);
        for (j, hj) in h.iter_mut().enumerate().skip(1) {
            let target = j as f64 * ds;
            while tk < target {
                hk += dh;
                tk = c.solid(hk);
            }
            // linear inside the fine cell
            let (h0, t0) = (hk - dh, c.solid(hk - dh));
            *hj = h0 + (target - t0) / (tk - t0) * dh;
        }
        let mut t = Self { c, ds_m: ds, h, max_err_m: 0.0 };
        let mut err = 0f64;
        for j in 0..ns - 1 {
            let s = (j as f64 + 0.5) * ds;
            err = err.max((c.solid(t.inv(s)) - s).abs());
        }
        t.max_err_m = err;
        t
    }
    /// H = T⁻¹(S) (m).
    pub fn inv(&self, s_m: f64) -> f64 {
        if s_m <= 0.0 {
            return 0.0;
        }
        let x = s_m / self.ds_m;
        let j = x.floor() as usize;
        if j + 1 < self.h.len() {
            let f = x - j as f64;
            self.h[j] + f * (self.h[j + 1] - self.h[j])
        } else {
            let hl = *self.h.last().unwrap();
            let sl = (self.h.len() - 1) as f64 * self.ds_m;
            let slope = 1.0 - self.c.phi0 as f64 * (-hl / self.c.z0_m as f64).exp();
            hl + (s_m - sl) / slope
        }
    }
}

impl Deposition {
    /// G of the declaration: Guerit et al. 2019 p. 855, natural systems « ~0 to 3.1 » → the middle, 1.55 (the bounds as
    /// sensitivity).
    pub const G_DECLARED: f32 = 1.55;
    /// The declaration: G = 1.55, the lakes filled, the marine silt.
    pub fn declared() -> Self {
        Self { g: Self::G_DECLARED, lakes: true, marine: Some(Marine::declared()), load: None, compaction: None, recycle: false, marine_live: false, fluvial_substeps: 1, fluvial_dt_max: f32::INFINITY, fluvial_keff: None, transmit: false, marine_blocks: false, soft_cover: false }
    }
    /// ADR Finding 188-H -- **the adopted deposition profile** (the author: « ok B adopté »): F187's configuration without
    /// K_eff (the load, the compaction, the recycling, the live marine step, G = 0, no K_eff) with F188's soft cover (B), the
    /// deposition handed to the amplifications. It goes with [`crate::cascade::Sea`] `{ connected: true }` at the physics level,
    /// C1's consumption for the recycling, and in the chain `AmpConfig::sea_strict` with `Transmit::bury_smooth`.
    /// Known defects (F188): the fluvial and marine stranding on dead ends (G1 bis: 3 487 km³ on seed 8), `bury_smooth`'s
    /// read through the parent cell, the step dependence (E3, V2 bis), the hillslopes' unloading (D2a).
    pub fn continent_default() -> Self {
        Self {
            g: 0.0,
            load: Some(Load::declared()),
            compaction: Some(Compaction::declared()),
            recycle: true,
            marine_live: true,
            soft_cover: true,
            transmit: true,
            ..Self::declared()
        }
    }
    /// Nothing deposited: the F181 level bit for bit (the A5 control).
    pub fn nothing() -> Self {
        Self { g: 0.0, lakes: false, marine: None, load: None, compaction: None, recycle: false, marine_live: false, fluvial_substeps: 1, fluvial_dt_max: f32::INFINITY, fluvial_keff: None, transmit: false, marine_blocks: false, soft_cover: false }
    }
}

/// One fluvial pass: what each cell received (norm; a height on one cell).
#[derive(Clone, Debug, Default)]
pub struct Pass {
    /// The fluvial deposit (eq. 13–15, capped by the incoming flux).
    pub fluvial: Vec<f64>,
    /// The lake fill (A2).
    pub lake: Vec<f64>,
    /// What reached a land base node (no receiver, no lake capacity): deposited there.
    pub stranded: Vec<f64>,
    /// What reached a sea cell (the marine source, norm × cells).
    pub to_sea: Vec<f64>,
    /// Σ of the eroded matter routed (norm × cells).
    pub eroded: f64,
    /// The lakes' outlets (one per sub-basin) and how many of them filled to the spill.
    pub outlets: usize,
    pub outlets_full: usize,
    /// ADR Finding 184 -- the fluvial deposit made under the cap (G / Ñ ≥ 1: the whole incoming flux), solid.
    pub capped: f64,
}

impl Pass {
    /// The deposit on cell `k` (norm): fluvial + lake + stranded.
    pub fn deposit(&self, k: usize) -> f64 {
        self.fluvial[k] + self.lake[k] + self.stranded[k]
    }
    /// ADR Finding 184-A2: adds another sub-step's pass (the arrays, the eroded matter, the cap; the last outlets).
    pub fn absorb(&mut self, o: &Pass) {
        for (a, b) in [(&mut self.fluvial, &o.fluvial), (&mut self.lake, &o.lake), (&mut self.stranded, &o.stranded), (&mut self.to_sea, &o.to_sea)] {
            for (x, y) in a.iter_mut().zip(b) {
                *x += y;
            }
        }
        self.eroded += o.eroded;
        self.capped += o.capped;
        self.outlets = o.outlets;
        self.outlets_full = o.outlets_full;
    }
    pub fn sums(&self) -> [f64; 4] {
        let s = |v: &[f64]| v.iter().sum::<f64>();
        [s(&self.fluvial), s(&self.lake), s(&self.stranded), s(&self.to_sea)]
    }
}

/// ADR Finding 182 A1–A2 -- **the fluvial pass**: route `e` (the step's eroded matter per cell, norm) over `surface`'s D8
/// steepest descent on the priority-flood surface (the incision's routing, periodic on a torus), upstream → downstream,
/// depositing min((G / Ñ_i) · Q_in, Q_in) on every non-lake cell and filling each lake by volume at its outlet. A base node
/// is a sea cell (`sea_mask`, else `surface ≤ sea`): its flux goes to the sea; a land base node keeps it.
///
/// ADR Finding 183-A2: `e` and every output are solid volumes; a new deposit's thickness is its solid / (1 − `p_dep`), so a
/// lake's capacity in solid is its volume × (1 − `p_dep`). `p_dep` = 0 is F182 exactly.
#[allow(clippy::too_many_arguments)]
pub fn pass(surface: &GridF32, e: &[f64], sea: f32, sea_mask: Option<&[bool]>, g: f64, lakes: bool, torus: bool, p_dep: f64) -> Pass {
    pass_routed(surface, e, sea, sea_mask, false, g, lakes, torus, p_dep)
}

/// ADR Finding 187-A2/A3 -- [`pass`] with `routed`: the `sea_mask` (the ocean by connectivity) is also the routing's sea --
/// the flood's seeds, the fixed receivers -- and a below-sea cell outside it is a closed depression, a lake filled by volume
/// up to its sill. `false` is [`pass`] exactly.
#[allow(clippy::too_many_arguments)]
pub fn pass_routed(surface: &GridF32, e: &[f64], sea: f32, sea_mask: Option<&[bool]>, routed: bool, g: f64, lakes: bool, torus: bool, p_dep: f64) -> Pass {
    let wsol = 1.0 - p_dep;
    let (w, h) = (surface.width, surface.height);
    let n = w * h;
    let ocean: Option<&[bool]> = if routed { Some(sea_mask.expect("the routed pass needs the ocean")) } else { None };
    let flow = match ocean {
        Some(o) => crate::terrain::flow::compute_flow_ocean(surface, &FlowConfig { sea_level: sea, ..Default::default() }, o),
        None => compute_flow(surface, &FlowConfig { sea_level: sea, ..Default::default() }),
    };
    let base = |k: usize| ocean.map_or(surface.data[k] <= sea, |o| o[k]);
    let mut receiver = vec![0usize; n];
    for k in 0..n {
        let d = flow.direction[k];
        if d == DIR_NONE || base(k) {
            receiver[k] = k;
            continue;
        }
        let nx = (k % w) as i32 + D8_DX[d as usize];
        let ny = (k / w) as i32 + D8_DY[d as usize];
        receiver[k] = if torus {
            ny.rem_euclid(h as i32) as usize * w + nx.rem_euclid(w as i32) as usize
        } else if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 {
            k
        } else {
            ny as usize * w + nx as usize
        };
    }
    // the stack, receivers before donors (stream_power's construction)
    let mut ndonor = vec![0u32; n];
    for k in 0..n {
        if receiver[k] != k {
            ndonor[receiver[k]] += 1;
        }
    }
    let mut offset = vec![0usize; n + 1];
    for k in 0..n {
        offset[k + 1] = offset[k] + ndonor[k] as usize;
    }
    let mut donors = vec![0usize; offset[n]];
    let mut fill = offset.clone();
    for k in 0..n {
        let r = receiver[k];
        if r != k {
            donors[fill[r]] = k;
            fill[r] += 1;
        }
    }
    let mut stack = Vec::with_capacity(n);
    let mut work = Vec::new();
    for k in 0..n {
        if receiver[k] == k {
            work.push(k);
            while let Some(node) = work.pop() {
                stack.push(node);
                for d in offset[node]..offset[node + 1] {
                    work.push(donors[d]);
                }
            }
        }
    }
    // Ñ_i (eq. 14): the catchment's cells, i included
    let mut cnt = vec![1f64; n];
    for &k in stack.iter().rev() {
        let r = receiver[k];
        if r != k {
            cnt[r] += cnt[k];
        }
    }
    // A2: the lakes (filled > z on land), each cell's outlet (the lake cell whose receiver leaves the lake), the
    // sub-basins sorted from the bottom up, their capacity
    let is_lake: Vec<bool> =
        (0..n).map(|k| lakes && ocean.map_or(surface.data[k] > sea, |o| !o[k]) && flow.filled.data[k] > surface.data[k]).collect();
    let mut outlet = vec![usize::MAX; n];
    for &k in &stack {
        if is_lake[k] {
            let r = receiver[k];
            outlet[k] = if r != k && is_lake[r] { outlet[r] } else { k };
        }
    }
    let mut groups: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
    for k in 0..n {
        if is_lake[k] {
            groups.entry(outlet[k]).or_default().push(k);
        }
    }
    for v in groups.values_mut() {
        v.sort_by(|&a, &b| surface.data[a].total_cmp(&surface.data[b]).then(a.cmp(&b)));
    }
    let mut p = Pass {
        fluvial: vec![0.0; n],
        lake: vec![0.0; n],
        stranded: vec![0.0; n],
        to_sea: vec![0.0; n],
        eroded: e.iter().sum(),
        outlets: groups.len(),
        outlets_full: 0,
        capped: 0.0,
    };
    let mut q = vec![0f64; n];
    for &k in stack.iter().rev() {
        let r = receiver[k];
        let qin = q[k];
        if r == k {
            let tot = qin + e[k];
            let is_sea = sea_mask.map_or(surface.data[k] <= sea, |m| m[k]);
            if is_sea {
                p.to_sea[k] += tot;
            } else {
                p.stranded[k] += tot;
            }
            continue;
        }
        let out = if is_lake[k] {
            let f = qin + e[k];
            if outlet[k] == k {
                // the sub-basin fills from the bottom up to the level L holding f, at most to its spill
                let cells = &groups[&k];
                let lvl = |c: usize| flow.filled.data[c] as f64;
                let cap: f64 = cells.iter().map(|&c| lvl(c) - surface.data[c] as f64).sum::<f64>() * wsol;
                let put = if f >= cap {
                    p.outlets_full += 1;
                    let mut s = 0.0;
                    for &c in cells {
                        let d = (lvl(c) - surface.data[c] as f64) * wsol;
                        p.lake[c] = d;
                        s += d;
                    }
                    s
                } else if f > 0.0 {
                    // the volume to hold, in thickness
                    let f = f / wsol;
                    let m = cells.len();
                    let (mut sum, mut level) = (0f64, 0f64);
                    let mut top = m;
                    for j in 0..m {
                        sum += surface.data[cells[j]] as f64;
                        let next = if j + 1 < m { (surface.data[cells[j + 1]] as f64).min(lvl(cells[j])) } else { lvl(cells[j]) };
                        if (j + 1) as f64 * next - sum >= f {
                            level = (f + sum) / (j + 1) as f64;
                            top = j + 1;
                            break;
                        }
                    }
                    if top == m && level == 0.0 {
                        level = (f + sum) / m as f64;
                    }
                    let mut s = 0.0;
                    for &c in &cells[..top] {
                        let d = (level.min(lvl(c)) - surface.data[c] as f64).max(0.0) * wsol;
                        p.lake[c] = d;
                        s += d;
                    }
                    s
                } else {
                    0.0
                };
                f - put
            } else {
                f
            }
        } else {
            let d = (g / cnt[k] * qin).min(qin).max(0.0);
            // F184: the cap's share (G / Ñ ≥ 1 deposits the whole incoming flux)
            if g >= cnt[k] {
                p.capped += d;
            }
            p.fluvial[k] = d;
            qin + e[k] - d
        };
        q[r] += out;
    }
    p
}

/// ADR Finding 182 A3 -- **the marine step**: the step's source `src` (solid volume, norm × cells, on the `marine` cells)
/// deposited at a thickness src / (1 − φ) and diffused by K_M in explicit sub-steps over `dt_yr` (Δt_sub ≤ Δx² / (4 K_M),
/// λ = K_M Δt_sub / Δx² ≤ 1/4), between marine cells only (no flux across the coast), each donor's outgoing flux limited
/// to its sediment `sed` (norm). Periodic on a torus. Returns the sub-steps run. Conservative: each link's flux is one
/// f64 expression read on both sides.
pub fn marine_step(z: &mut GridF32, marine: &[bool], src: &[f64], dt_yr: f64, m: &Marine, cell_m: f64, sed: &mut [f64], torus: bool) -> usize {
    let (w, h) = (z.width, z.height);
    let n = w * h;
    let k_m = m.k_m2_yr as f64;
    let thick = 1.0 / (1.0 - m.porosity as f64);
    let idx: Vec<usize> = (0..n).filter(|&k| marine[k]).collect();
    let any = idx.iter().any(|&k| src[k] > 0.0 || sed[k] > 0.0);
    if !any {
        return 0;
    }
    if k_m <= 0.0 || dt_yr <= 0.0 {
        for &k in &idx {
            let t = src[k] * thick;
            sed[k] += t;
            z.data[k] = (z.data[k] as f64 + t) as f32;
        }
        return 0;
    }
    let nsub = (dt_yr * 4.0 * k_m / (cell_m * cell_m)).ceil().max(1.0) as usize;
    let lambda = k_m * (dt_yr / nsub as f64) / (cell_m * cell_m);
    // the four neighbours of each marine cell that are marine (usize::MAX = none)
    let nb = |k: usize, d: usize| -> usize {
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        let (dx, dy) = [(1, 0), (-1, 0), (0, 1), (0, -1)][d];
        let (nx, ny) = (x + dx, y + dy);
        let j = if torus {
            ny.rem_euclid(h as i32) as usize * w + nx.rem_euclid(w as i32) as usize
        } else if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 {
            return usize::MAX;
        } else {
            ny as usize * w + nx as usize
        };
        if marine[j] { j } else { usize::MAX }
    };
    let links: Vec<[usize; 4]> = idx.iter().map(|&k| [nb(k, 0), nb(k, 1), nb(k, 2), nb(k, 3)]).collect();
    let mut hz: Vec<f64> = z.data.iter().map(|&v| v as f64).collect();
    let add: Vec<f64> = idx.iter().map(|&k| src[k] * thick / nsub as f64).collect();
    let mut sc = vec![0f64; n];
    let mut delta = vec![0f64; idx.len()];
    for _ in 0..nsub {
        for (i, &k) in idx.iter().enumerate() {
            hz[k] += add[i];
            sed[k] += add[i];
        }
        // the donors' supply limit (serial: a sub-step is ~10⁵ flops, below rayon's overhead)
        for (i, &k) in idx.iter().enumerate() {
            let mut out = 0.0;
            for &j in &links[i] {
                if j != usize::MAX && hz[k] > hz[j] {
                    out += lambda * (hz[k] - hz[j]);
                }
            }
            sc[k] = if out > 0.0 { (sed[k].max(0.0) / out).min(1.0) } else { 0.0 };
        }
        for (i, &k) in idx.iter().enumerate() {
            let mut d = 0.0;
            for &j in &links[i] {
                if j == usize::MAX {
                    continue;
                }
                if hz[j] > hz[k] {
                    d += lambda * (hz[j] - hz[k]) * sc[j];
                } else if hz[k] > hz[j] {
                    d -= lambda * (hz[k] - hz[j]) * sc[k];
                }
            }
            delta[i] = d;
        }
        for (i, &k) in idx.iter().enumerate() {
            hz[k] += delta[i];
            sed[k] += delta[i];
        }
    }
    for &k in &idx {
        z.data[k] = hz[k] as f32;
    }
    nsub
}

/// ADR Finding 182 A4 -- the weights spreading the hillslope's loss to the sea over the coastal sea cells: Σ over a sea
/// cell's land 4-neighbours of (z_land − z_sea)⁺, the gradient that drives the linear diffusion and the talus across the
/// coast.
pub fn coast_weights(z: &GridF32, land: &[bool], torus: bool) -> Vec<f64> {
    let (w, h) = (z.width, z.height);
    (0..w * h)
        .map(|k| {
            if land[k] {
                return 0.0;
            }
            let (x, y) = ((k % w) as i32, (k / w) as i32);
            let mut s = 0.0;
            for (dx, dy) in [(1, 0), (-1, 0), (0, 1), (0, -1)] {
                let (nx, ny) = (x + dx, y + dy);
                let j = if torus {
                    ny.rem_euclid(h as i32) as usize * w + nx.rem_euclid(w as i32) as usize
                } else if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 {
                    continue;
                } else {
                    ny as usize * w + nx as usize
                };
                if land[j] {
                    s += (z.data[j] as f64 - z.data[k] as f64).max(0.0);
                }
            }
            s
        })
        .collect()
}

/// The step's fluvial stage as a [`FluvialHook`] sees it: the field after the uplift, the incision (the cascade's own,
/// the rock law included), the pass (the cascade's own), and Yuan's tolerance base U·Δt (the step's mean |U Δt|, norm).
pub struct FluvialCtx<'a> {
    pub before: &'a GridF32,
    pub incise: &'a dyn Fn(&GridF32) -> GridF32,
    pub pass: &'a dyn Fn(&GridF32, &[f64]) -> Pass,
    pub u_dt: f64,
    pub step: usize,
    /// ADR Finding 183-A2: the porosity of a new deposit (its thickness = solid / (1 − p)); 0 = F182.
    pub p_dep: f64,
}

/// A fluvial stage's result: the field, the gross incision per cell (norm, for the rebound and Σerosion), the deposits.
pub struct FluvialOut {
    pub z: GridF32,
    pub e: Vec<f64>,
    pub pass: Pass,
    pub iterations: usize,
}

/// ADR Finding 182 -- a replacement of the fluvial stage, for the bench's reference only (Yuan's iterative coupling).
pub type FluvialHook<'h> = dyn Fn(&FluvialCtx<'_>) -> FluvialOut + 'h;

/// ADR Finding 182 A1 -- **the explicit splitting** (the cascade's stage): the incision unchanged, then one pass routing
/// its eroded matter on the incised surface; the field is the incised one plus the deposits.
pub fn splitting(c: &FluvialCtx) -> FluvialOut {
    let hs = (c.incise)(c.before);
    let e: Vec<f64> = c.before.data.iter().zip(&hs.data).map(|(&b, &a)| (b as f64 - a as f64).max(0.0)).collect();
    let p = (c.pass)(&hs, &e);
    let mut z = hs;
    let w = 1.0 - c.p_dep;
    for k in 0..z.data.len() {
        let d = p.deposit(k);
        if d != 0.0 {
            z.data[k] = (z.data[k] as f64 + d / w) as f32;
        }
    }
    FluvialOut { z, e, pass: p, iterations: 1 }
}

/// ADR Finding 182 -- the level's deposition record (m unless stated).
#[derive(Clone, Debug, Default)]
pub struct SedRecord {
    /// The sediment thickness at the end (m): the deposits, less what the incision and the hillslope lowered since.
    pub thickness_m: Vec<f32>,
    /// ADR Finding 188-B2 -- Σ of the marine walk's stranded thickness per cell (m; empty unless the live marine step ran).
    pub marine_stranded_m: Vec<f32>,
    /// ADR Finding 189 2.3 -- Σ of what left the map because the ocean was filled with no deeper cell (solid, norm × cells).
    pub marine_exported: f64,
    /// ADR Finding 189 2.3 -- Σ of the marine sources out of the traversable sea, routed by the fluvial pass (solid, norm ×
    /// cells; none in the adopted profile).
    pub marine_rerouted: f64,
    /// Σ over the level of the fluvial deposits, the lake fills and the marine thickness change (m per cell).
    pub fluvial_m: Vec<f32>,
    pub lake_m: Vec<f32>,
    pub marine_m: Vec<f32>,
    /// The budget over the level (norm × cells, solid): [eroded by the incision, the hillslope's loss to the sea, fluvial,
    /// lakes, stranded, marine solid, exported; F183: subducted, the cover eroded by the incision, the cover lowered by the
    /// hillslope, without A2 the solid created by counting a porous cover's height as solid].
    pub budget: [f64; 11],
    /// ADR Finding 183-A4 -- the worst per-step relative residual of the stock control (in = ΔStock + the cover eroded + the
    /// cover lowered + subducted + exported, ΔStock measured on the cover's solid field).
    pub max_stock_residual: f64,
    /// ADR Finding 183 -- the cover's solid at the end (m of solid per cell); Σ |ΔL_cover| over the level (A1, kg/m² / g in
    /// norm units); Σ |ΔH| of the compaction (norm); the compaction table's worst round trip (m).
    pub solid_m: Vec<f32>,
    pub load_cover_abs: f64,
    /// ADR Finding 188 D2 (a record only) -- the uniform mode (k = 0) of F183 A1's flexure, cumulated over the steps (norm): its
    /// erosional part, mean(ρ_c e_basement / ρ_m), and its load part, mean(−ΔL / ρ_m); and Σ e_basement (norm × cells), the
    /// basement the incision removed (the hillslope's lowering is not in it).
    pub w0_erosion_n: f64,
    pub w0_load_n: f64,
    pub basement_incised_n: f64,
    pub compaction_change: f64,
    pub table_err_m: f64,
    /// ADR Finding 184 -- the fluvial deposit made under the cap (solid, norm × cells); the marine source stranded in a closed
    /// depression by the live domain (solid).
    pub capped: f64,
    pub marine_stranded: f64,
    /// The worst per-step relative residual of the budget, and the worst f32-field residual (the deposits as written in z
    /// against the pass's f64).
    pub max_rel_residual: f64,
    pub max_field_residual: f64,
    pub marine_substeps: usize,
    /// Per step: its duration (yr), the marine sub-steps, the marine source (norm × cells), the sea cells, the eroded matter,
    /// the sediment on the sea cells before the step (norm × cells).
    pub steps: Vec<(f64, usize, f64, usize, f64, f64)>,
    /// The fluvial stage's iterations per step (1 for the splitting).
    pub iterations: Vec<usize>,
    /// Wall seconds: the fluvial stage beyond the incision alone is not separable; the marine step.
    pub secs_marine: f64,
    /// The lakes' outlets and those full at the last step.
    pub last_outlets: (usize, usize),
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A5: G = 0 deposits nothing on land; everything reaches the sea or a lake, the mass exact.
    #[test]
    fn the_pass_conserves_and_g0_deposits_nothing_on_land() {
        let n = 32;
        let mut z = GridF32::new(n, n, 0.0);
        for k in 0..n * n {
            let (x, y) = ((k % n) as f32, (k / n) as f32);
            z.data[k] = 0.5 + 0.02 * (1.0 + (x * 0.4).sin() * (y * 0.3).cos()) * ((x - 16.0).abs().min(16.0) / 16.0);
        }
        for x in 0..n {
            z.data[x] = 0.4;
        }
        let e: Vec<f64> = z.data.iter().map(|&v| if v > 0.5 { 1e-4 } else { 0.0 }).collect();
        for g in [0.0, 1.55, 3.1] {
            for lakes in [false, true] {
                let p = pass(&z, &e, 0.5, None, g, lakes, true, 0.0);
                let s = p.sums();
                let tot: f64 = s.iter().sum();
                assert!((tot - p.eroded).abs() <= 1e-12 * p.eroded, "g {g} lakes {lakes}: {tot} vs {}", p.eroded);
                if g == 0.0 {
                    assert_eq!(s[0], 0.0);
                }
                for k in 0..n * n {
                    assert!(p.fluvial[k] >= 0.0 && p.lake[k] >= 0.0);
                }
            }
        }
    }

    /// ADR Finding 186 (the fallback): the block step conserves the source exactly, and no cell gains above the sea but the
    /// stranded part.
    #[test]
    fn the_block_marine_step_conserves() {
        let n = 32;
        let sea = 0.5f64;
        let mut z = GridF32::new(n, n, 0.0);
        for k in 0..n * n {
            let x = (k % n) as f32;
            z.data[k] = if x < 10.0 { 0.53 + 0.01 * (10.0 - x) } else { 0.495 - 0.002 * (x - 10.0) };
        }
        let mut src = vec![0f64; n * n];
        src[6 * n + 4] = 0.04;
        src[12 * n + 12] = 0.01;
        let tot: f64 = src.iter().sum();
        let m = Marine::declared();
        let mut sup = vec![0f64; n * n];
        let z0 = z.clone();
        let (ns, strd) = marine_step_live_blocks(&mut z, &src, 2.0e5, &m, 1562.5, &mut sup, true, sea);
        assert!(ns >= 1);
        let solid: f64 = sup.iter().sum::<f64>() * (1.0 - m.porosity as f64);
        assert!((solid - tot).abs() < 1e-12 * tot, "{solid} against {tot}");
        let above: f64 = (0..n * n).filter(|&k| z0.data[k] < 0.5).map(|k| (z.data[k] as f64 - sea).max(0.0)).sum();
        assert!(above <= strd + 1e-6, "above the sea {above} against stranded {strd}");
        // F186 step 2: the bilinear prolongation keeps the mass off the torus too (the edge clamp, the residue redistributed)
        let (mut z2, mut sup2) = (z0.clone(), vec![0f64; n * n]);
        marine_step_live_blocks(&mut z2, &src, 2.0e5, &m, 1562.5, &mut sup2, false, sea);
        let solid2: f64 = sup2.iter().sum::<f64>() * (1.0 - m.porosity as f64);
        assert!((solid2 - tot).abs() < 1e-9 * tot, "{solid2} against {tot} (bounded)");
    }

    /// ADR Finding 184-A1: the live marine step conserves the source, fills the sea cells only up to the sea, walks a source
    /// landing on land down to the sea, and strands it in a closed depression.
    #[test]
    fn the_live_marine_step_fills_to_the_sea_and_walks() {
        let n = 24;
        let sea = 0.5f64;
        // a ramp: land on x < 8 (rising inland), sea beyond, and a closed land pit at (3, 3)
        let mut z = GridF32::new(n, n, 0.0);
        for k in 0..n * n {
            let x = (k % n) as f32;
            z.data[k] = if x < 8.0 { 0.52 + 0.01 * (8.0 - x) } else { 0.49 - 0.002 * (x - 8.0) };
        }
        z.data[3 * n + 3] = 0.51;
        let mut src = vec![0f64; n * n];
        src[5 * n + 2] = 0.05; // on land: walks to the sea
        src[3 * n + 3] = 1e-3; // in the pit: stranded
        src[10 * n + 9] = 0.02; // on the sea
        let tot_src: f64 = src.iter().sum();
        let m = Marine::declared();
        let mut sup = vec![0f64; n * n];
        let z0 = z.clone();
        let (ns, strd) = marine_step_live(&mut z, &src, 2.0e5, &m, 1562.5, &mut sup, true, sea);
        assert!(ns > 1);
        let solid: f64 = sup.iter().sum::<f64>() * (1.0 - m.porosity as f64);
        assert!((solid - tot_src).abs() < 1e-12 * tot_src, "{solid} against {tot_src}");
        assert!(strd > 0.0, "the pit strands its source");
        for k in 0..n * n {
            let gained = z.data[k] as f64 - z0.data[k] as f64;
            if z0.data[k] < 0.5 && gained > 0.0 {
                assert!(z.data[k] as f64 <= sea + 1e-6, "a sea cell {k} filled above the sea: {}", z.data[k]);
            }
        }
        assert!(sup[5 * n + 2] == 0.0, "the land source walked away");
    }

    /// ADR Finding 183-A2: the table inverts T to < 1 mm; z₀ = ∞ is linear; a porous pass conserves its solid.
    #[test]
    fn the_compaction_table_and_the_porous_pass() {
        let t = CompactionTable::new(Compaction::declared());
        assert!(t.max_err_m < 1e-3, "{}", t.max_err_m);
        for h in [0.0, 1.0, 100.0, 1960.0, 8000.0, 29_000.0, 40_000.0] {
            let s = t.c.solid(h);
            assert!((t.inv(s) - h).abs() < 2e-3, "h {h}: {}", t.inv(s));
        }
        let lin = CompactionTable::new(Compaction { phi0: 0.63, z0_m: f32::INFINITY });
        for s in [0.5, 10.0, 1234.5] {
            assert!((lin.inv(s) - s / 0.37).abs() < 1e-6 * s);
        }
        let n = 32;
        let mut z = GridF32::new(n, n, 0.6);
        for x in 0..n {
            z.data[x] = 0.4;
        }
        z.data[10 * n + 10] = 0.55;
        z.data[10 * n + 11] = 0.56;
        let e: Vec<f64> = z.data.iter().map(|&v| if v > 0.5 { 1e-4 } else { 0.0 }).collect();
        let p = pass(&z, &e, 0.5, None, 1.55, true, true, 0.63);
        assert!((p.sums().iter().sum::<f64>() - p.eroded).abs() <= 1e-12 * p.eroded);
    }

    /// A3: the marine step conserves the solid volume and never erodes the bathymetry below its sediment.
    #[test]
    fn the_marine_step_conserves_and_moves_only_sediment() {
        let n = 24;
        let mut z = GridF32::new(n, n, 0.0);
        for k in 0..n * n {
            z.data[k] = 0.2 + 0.01 * ((k % n) as f32);
        }
        let marine = vec![true; n * n];
        let mut src = vec![0f64; n * n];
        src[5 * n + 20] = 1e-3;
        let mut sed = vec![0f64; n * n];
        let z0 = z.clone();
        let m = Marine::declared();
        let ns = marine_step(&mut z, &marine, &src, 2.0e5, &m, 1562.5, &mut sed, true);
        assert!(ns > 1);
        let solid: f64 = sed.iter().sum::<f64>() * (1.0 - m.porosity as f64);
        assert!((solid - 1e-3).abs() < 1e-12, "{solid}");
        for k in 0..n * n {
            assert!(sed[k] > -1e-15);
            assert!(z.data[k] >= z0.data[k] - 1e-7);
        }
    }
}
