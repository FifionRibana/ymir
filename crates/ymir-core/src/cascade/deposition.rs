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
use crate::terrain::flow::{D8_DX, D8_DY, DIR_NONE, FlowConfig, compute_flow};

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
        Self { g: Self::G_DECLARED, lakes: true, marine: Some(Marine::declared()), load: None, compaction: None, recycle: false }
    }
    /// Nothing deposited: the F181 level bit for bit (the A5 control).
    pub fn nothing() -> Self {
        Self { g: 0.0, lakes: false, marine: None, load: None, compaction: None, recycle: false }
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
}

impl Pass {
    /// The deposit on cell `k` (norm): fluvial + lake + stranded.
    pub fn deposit(&self, k: usize) -> f64 {
        self.fluvial[k] + self.lake[k] + self.stranded[k]
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
    let wsol = 1.0 - p_dep;
    let (w, h) = (surface.width, surface.height);
    let n = w * h;
    let flow = compute_flow(surface, &FlowConfig { sea_level: sea, ..Default::default() });
    let mut receiver = vec![0usize; n];
    for k in 0..n {
        let d = flow.direction[k];
        if d == DIR_NONE || surface.data[k] <= sea {
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
        (0..n).map(|k| lakes && surface.data[k] > sea && flow.filled.data[k] > surface.data[k]).collect();
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
    pub compaction_change: f64,
    pub table_err_m: f64,
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
