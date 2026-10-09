//! ADR Finding 159 -- **the multi-scale cascade prototype** (gated: nothing in production calls it; the viz's Cascade
//! window and the `f159_*` benches do).
//!
//! The author's LOD track (2026-10-08): « raffinant une carte d'abord à 64², puis 128², puis 256² … Chaque niveau
//! produirait la forme ou le relief répondant à la résolution ». Each level is upscaled from the previous one (bicubic
//! ×2), then run to a steady state under tectonic uplift, the existing implicit stream power, the talus and the linear
//! diffusion. The design and its every number are declared in
//! `docs/reports/relief_method/f159_cascade/f159_declared.md`; the short version:
//!
//! * **U(x) = U₀ · max(0, h_iso) / 1 000 m** (m/yr), h_iso the coarse C1 isostatic altitude; constant in time, so a
//!   level HAS a steady state and the steady state is the stop (D7). U₀ is calibrated once at the first level so the
//!   steady state's mean land altitude is h_iso's (n = 1 makes the steady state linear in U).
//! * **K is physical and the same at every level** (Kwang & Parker 2017: with m/n = 0.5 and no diffusion the law is
//!   scale-free, so K cannot carry a length). **The length is D's** (Perron 2009, `L_c = √(D/K)`):
//!   `D_L = max(D_phys, α·K_SI·cell²)`, so each level carries valleys down to 2.6–5.1 of its cells.
//! * The sea is the base level: during a level its cells are held at sea level and its land mask is fixed; the
//!   bathymetry is restored at the end of the level.
//!
//! Read-only on its input: [`Cascade::new`] takes the coarse field by value and nothing here touches a tectonic state.

use crate::erosion::stream_power::{
    RELIEF_V1_A_C_KM2, StreamPowerConfig, incise_with_floor, linear_diffusion, talus_sweep,
};
use crate::grid::GridF32;
use std::time::Instant;

/// The four sub-steps of a level, in the order the author reads them (the viz frieze).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SubStep {
    /// « Agrandissement »: the bicubic ×2 of the previous level.
    Upscale,
    /// « Soulèvement »: the upscaled field plus the level's whole uplift (Σ U·dt).
    Uplift,
    /// « Érosion »: plus the stream power's whole change (incision + lateral erosion).
    Erosion,
    /// « Diffusion »: plus the hillslope closures' (talus + linear diffusion) = the level's result.
    Diffusion,
}

impl SubStep {
    pub const ALL: [SubStep; 4] = [SubStep::Upscale, SubStep::Uplift, SubStep::Erosion, SubStep::Diffusion];
    pub fn label(self) -> &'static str {
        match self {
            SubStep::Upscale => "Agrandissement",
            SubStep::Uplift => "Soulèvement",
            SubStep::Erosion => "Érosion",
            SubStep::Diffusion => "Diffusion",
        }
    }
}

/// The declared parameters (`f159_declared.md` D1–D7). [`CascadeConfig::declared`] is the round's.
#[derive(Clone, Debug, PartialEq)]
pub struct CascadeConfig {
    /// The cascade's levels after the 64² tectonics, each the double of the previous one.
    pub levels: Vec<usize>,
    /// The domain's side (km): the témoin's 400 km.
    pub domain_km: f32,
    /// `SteinSteinParams::depth_scale_m` (the norm ↔ metres contract).
    pub depth_scale_m: f32,
    /// Sea level (norm).
    pub sea_level: f32,
    /// K in the code's units (E [m/yr] = K · A_km²^0.5 · S); K_SI = K / 1000 for m = 0.5.
    pub k: f32,
    /// Years per step.
    pub dt_yr: f32,
    /// The sub-grid diffusivity's factor: `D = α · K_SI · cell²` gives `L_c = √α` cells.
    pub alpha: f32,
    /// The physical diffusivity (m²/yr), PROXY: `K_SI · A_c` (Perron's `L_c² ≈ A_c`).
    pub d_phys_m2_yr: f32,
    /// U₀ (m/yr per 1 000 m of h_iso). The trial value when `calibrate` is on.
    pub u0_m_yr: f32,
    /// Calibrate U₀ at the first level (D2).
    pub calibrate: bool,
    /// Equilibrium: mean |Δz| per step < `eq_ratio` × mean U·dt, over land …
    pub eq_ratio: f32,
    /// … for this many consecutive steps.
    pub eq_window: usize,
    /// The step cap per level (same length as `levels`).
    pub max_steps: Vec<usize>,
    /// The linear diffusion's explicit sub-steps.
    pub diffusion_substeps: usize,
}

impl CascadeConfig {
    /// The F159 declaration (D1–D7).
    pub fn declared(depth_scale_m: f32) -> Self {
        let k = 2.0e-3;
        Self {
            levels: vec![128, 256, 512],
            domain_km: 400.0,
            depth_scale_m,
            sea_level: 0.5,
            k,
            dt_yr: 1.0e5,
            alpha: 0.16,
            d_phys_m2_yr: k / 1000.0 * RELIEF_V1_A_C_KM2 * 1.0e6,
            u0_m_yr: 1.0e-4,
            calibrate: true,
            eq_ratio: 0.01,
            eq_window: 5,
            max_steps: vec![400, 300, 200],
            diffusion_substeps: 4,
        }
    }
    /// Metres per norm unit (`c1_altitude_norm_to_metres`' slope).
    pub fn norm_to_m(&self) -> f32 {
        2.0 * 1.13 * self.depth_scale_m
    }
    pub fn cell_km(&self, n: usize) -> f32 {
        self.domain_km / n as f32
    }
    /// D at a level of `n` cells (m²/yr): `max(D_phys, α·K_SI·cell²)`.
    pub fn diffusivity(&self, n: usize) -> f32 {
        let cell_m = self.cell_km(n) * 1000.0;
        (self.alpha * self.k / 1000.0 * cell_m * cell_m).max(self.d_phys_m2_yr)
    }
    /// The linear diffusion's dimensionless weight per step, `D·dt/cell²`.
    pub fn diffusion_weight(&self, n: usize) -> f32 {
        let cell_m = self.cell_km(n) * 1000.0;
        self.diffusivity(n) * self.dt_yr / (cell_m * cell_m)
    }
    /// The solver's config at a level: `relief_v3` at the level's cell, with this round's K and dt, one iteration, and
    /// the talus and diffusion OFF inside the call (they are the « Diffusion » sub-step).
    pub fn solver(&self, n: usize) -> StreamPowerConfig {
        let cell_km = self.cell_km(n);
        let mut sp = StreamPowerConfig::relief_v3(cell_km * cell_km, self.depth_scale_m);
        sp.k = self.k;
        sp.dt = self.dt_yr;
        sp.iterations = 1;
        sp.sea_level = self.sea_level;
        sp.diffusion = 0.0;
        sp.talus_slope = 0.0;
        sp
    }
    /// The hillslope config at a level: the talus and the linear diffusion at the level's weight.
    pub fn hillslope(&self, n: usize) -> StreamPowerConfig {
        let cell_km = self.cell_km(n);
        let mut sp = StreamPowerConfig::relief_v3(cell_km * cell_km, self.depth_scale_m);
        sp.sea_level = self.sea_level;
        sp.diffusion = self.diffusion_weight(n);
        sp.diffusion_substeps = self.diffusion_substeps;
        sp
    }
}

/// Progress of a running level (the viz's bar, the benches' log).
#[derive(Clone, Copy, Debug)]
pub struct CascadeProgress {
    /// `None` = the calibration run.
    pub level: Option<usize>,
    pub n: usize,
    pub step: usize,
    pub max_steps: usize,
    /// mean |Δz| / mean U·dt of the last step.
    pub ratio: f32,
}

/// One computed level, with its additive decomposition: `upscaled + Σuplift + Σerosion + Σdiffusion = z` on land
/// (exact up to f32 rounding; the sea is restored to `upscaled`).
#[derive(Clone, Debug)]
pub struct LevelRecord {
    pub n: usize,
    pub cell_km: f32,
    /// D (m²/yr) and the weight per step.
    pub d_m2_yr: f32,
    pub weight: f32,
    /// The bicubic of the previous level (norm), sea included.
    pub upscaled: GridF32,
    /// U at this level (m/yr).
    pub uplift: GridF32,
    /// The land mask the level ran on (fixed during the level).
    pub land: Vec<bool>,
    /// Σ of each operator's change over the level (norm), 0 on the sea.
    pub sum_uplift: Vec<f32>,
    pub sum_erosion: Vec<f32>,
    pub sum_diffusion: Vec<f32>,
    /// The talus's part of `sum_diffusion` (norm), for the talus activity report.
    pub sum_talus: Vec<f32>,
    /// The level's result (norm), the bathymetry restored.
    pub z: GridF32,
    pub steps: usize,
    pub equilibrium: bool,
    /// The ratio per step.
    pub ratios: Vec<f32>,
    /// Wall seconds: upscale, uplift, erosion, diffusion.
    pub secs: [f64; 4],
    /// Land cells the talus moved at the last step (share of land).
    pub talus_last_share: f32,
}

impl LevelRecord {
    /// The field after a sub-step (norm): the additive decomposition, cumulated.
    pub fn substep_field(&self, s: SubStep) -> GridF32 {
        let mut g = self.upscaled.clone();
        let parts: &[&Vec<f32>] = match s {
            SubStep::Upscale => &[],
            SubStep::Uplift => &[&self.sum_uplift],
            SubStep::Erosion => &[&self.sum_uplift, &self.sum_erosion],
            SubStep::Diffusion => return self.z.clone(),
        };
        for p in parts {
            for (v, d) in g.data.iter_mut().zip(p.iter()) {
                *v += d;
            }
        }
        g
    }
    pub fn total_secs(&self) -> f64 {
        self.secs.iter().sum()
    }
    /// Retained bytes (the memory estimate of D8.9): 4 grids + 4 sums + the mask.
    pub fn retained_bytes(&self) -> usize {
        let n = self.n * self.n;
        n * 4 * 8 + n
    }
}

/// The calibration run's record (D2).
#[derive(Clone, Debug)]
pub struct Calibration {
    pub u0_trial: f32,
    pub target_mean_m: f32,
    pub trial_mean_m: f32,
    pub u0: f32,
    pub steps: usize,
    pub equilibrium: bool,
    pub secs: f64,
    pub talus_last_share: f32,
}

/// The cascade: the coarse field (64², rolled to the framing), the uplift, and the levels computed so far.
pub struct Cascade {
    pub cfg: CascadeConfig,
    pub coarse: GridF32,
    /// U at 64² (m/yr), at the CALIBRATED U₀ once the calibration ran.
    pub u_coarse: GridF32,
    pub calibration: Option<Calibration>,
    pub levels: Vec<LevelRecord>,
}

impl Cascade {
    /// `coarse` is the 64² normalised C1 altitude (sea 0.5), ALREADY rolled to the framing ([`roll`]).
    pub fn new(coarse: GridF32, cfg: CascadeConfig) -> Self {
        let u_coarse = uplift_field(&coarse, cfg.u0_m_yr, &cfg);
        Self { cfg, coarse, u_coarse, calibration: None, levels: Vec::new() }
    }

    /// Is every declared level computed?
    pub fn done(&self) -> bool {
        self.levels.len() >= self.cfg.levels.len()
    }

    /// Compute the next level (running the calibration first when it is due). `None` when all are done or when
    /// `cancel` fired (nothing is kept from a cancelled level).
    pub fn next_level(
        &mut self,
        progress: &mut dyn FnMut(CascadeProgress),
        cancel: &dyn Fn() -> bool,
    ) -> Option<&LevelRecord> {
        let li = self.levels.len();
        if li >= self.cfg.levels.len() {
            return None;
        }
        let n = self.cfg.levels[li];
        if li == 0 && self.cfg.calibrate && self.calibration.is_none() {
            let t = Instant::now();
            let z0 = upsample_to(&self.coarse, n);
            let u = upsample_to(&self.u_coarse, n).map_nonneg();
            let target = mean_land_m(&z0, &land_mask(&z0, self.cfg.sea_level), &self.cfg);
            let trial = run_level(z0, u, &self.cfg, n, self.cfg.max_steps[0], None, progress, cancel)?;
            let trial_mean = mean_land_m(&trial.z, &trial.land, &self.cfg);
            let u0 = self.cfg.u0_m_yr * target / trial_mean.max(1e-3);
            self.calibration = Some(Calibration {
                u0_trial: self.cfg.u0_m_yr,
                target_mean_m: target,
                trial_mean_m: trial_mean,
                u0,
                steps: trial.steps,
                equilibrium: trial.equilibrium,
                secs: t.elapsed().as_secs_f64(),
                talus_last_share: trial.talus_last_share,
            });
            self.u_coarse = uplift_field(&self.coarse, u0, &self.cfg);
        }
        let t = Instant::now();
        let (z0, u) = match self.levels.last() {
            None => (upsample_to(&self.coarse, n), upsample_to(&self.u_coarse, n).map_nonneg()),
            Some(prev) => (upsample2(&prev.z), upsample2(&prev.uplift).map_nonneg()),
        };
        let up_secs = t.elapsed().as_secs_f64();
        let mut rec = run_level(z0, u, &self.cfg, n, self.cfg.max_steps[li], Some(li), progress, cancel)?;
        rec.secs[0] = up_secs;
        self.levels.push(rec);
        self.levels.last()
    }

    /// Every remaining level.
    pub fn run_all(&mut self, progress: &mut dyn FnMut(CascadeProgress), cancel: &dyn Fn() -> bool) {
        while !self.done() {
            if self.next_level(progress, cancel).is_none() {
                return;
            }
        }
    }
}

trait NonNeg {
    fn map_nonneg(self) -> Self;
}
impl NonNeg for GridF32 {
    fn map_nonneg(mut self) -> Self {
        for v in &mut self.data {
            *v = v.max(0.0);
        }
        self
    }
}

/// U = U₀ · max(0, h_iso) / 1 000 m (m/yr) on a normalised field.
pub fn uplift_field(z: &GridF32, u0: f32, cfg: &CascadeConfig) -> GridF32 {
    let n2m = cfg.norm_to_m();
    let mut u = z.clone();
    for v in &mut u.data {
        *v = u0 * ((*v - cfg.sea_level) * n2m).max(0.0) / 1000.0;
    }
    u
}

/// The framing roll: pixel (i, j) of the result is pixel (i + ox, j + oy) of `g`, wrapped (D1: the coarse coordinate
/// of level pixel i is `ox + i·64/N`).
pub fn roll(g: &GridF32, offset: [usize; 2]) -> GridF32 {
    let (w, h) = (g.width, g.height);
    let mut out = GridF32::new(w, h, 0.0);
    for y in 0..h {
        for x in 0..w {
            out.data[y * w + x] = g.data[((y + offset[1]) % h) * w + (x + offset[0]) % w];
        }
    }
    out
}

/// Periodic Catmull-Rom bicubic ×2 on the cell-corner mapping: even pixels are `g`'s values exactly, odd pixels the
/// cubic at the midpoint, `(−a + 9b + 9c − d) / 16`. Separable (rows, then columns).
pub fn upsample2(g: &GridF32) -> GridF32 {
    let (w, h) = (g.width, g.height);
    let mid = |a: f32, b: f32, c: f32, d: f32| (-a + 9.0 * b + 9.0 * c - d) / 16.0;
    let (w2, h2) = (2 * w, 2 * h);
    let mut rows = vec![0f32; w2 * h];
    for y in 0..h {
        let r = &g.data[y * w..(y + 1) * w];
        for x in 0..w {
            rows[y * w2 + 2 * x] = r[x];
            rows[y * w2 + 2 * x + 1] = mid(r[(x + w - 1) % w], r[x], r[(x + 1) % w], r[(x + 2) % w]);
        }
    }
    let mut out = GridF32::new(w2, h2, 0.0);
    for y in 0..h {
        for x in 0..w2 {
            let at = |yy: usize| rows[(yy % h) * w2 + x];
            out.data[2 * y * w2 + x] = at(y);
            out.data[(2 * y + 1) * w2 + x] = mid(at(y + h - 1), at(y), at(y + 1), at(y + 2));
        }
    }
    out
}

/// Repeated [`upsample2`] up to `n` (a power-of-two multiple of `g`'s size).
pub fn upsample_to(g: &GridF32, n: usize) -> GridF32 {
    let mut out = g.clone();
    while out.width < n {
        out = upsample2(&out);
    }
    assert_eq!(out.width, n, "the cascade's levels double: {} does not reach {n}", g.width);
    out
}

pub fn land_mask(z: &GridF32, sea: f32) -> Vec<bool> {
    z.data.iter().map(|&v| v > sea).collect()
}

/// Mean land altitude (m) over a mask.
pub fn mean_land_m(z: &GridF32, land: &[bool], cfg: &CascadeConfig) -> f32 {
    let (mut s, mut c) = (0f64, 0usize);
    for (v, &l) in z.data.iter().zip(land) {
        if l {
            s += ((v - cfg.sea_level) * cfg.norm_to_m()) as f64;
            c += 1;
        }
    }
    (s / c.max(1) as f64) as f32
}

/// One level to its steady state (D3, D7). `None` if cancelled.
#[allow(clippy::too_many_arguments)]
fn run_level(
    upscaled: GridF32,
    uplift: GridF32,
    cfg: &CascadeConfig,
    n: usize,
    max_steps: usize,
    level: Option<usize>,
    progress: &mut dyn FnMut(CascadeProgress),
    cancel: &dyn Fn() -> bool,
) -> Option<LevelRecord> {
    let cells = n * n;
    let land = land_mask(&upscaled, cfg.sea_level);
    let n2m = cfg.norm_to_m();
    let solver = cfg.solver(n);
    let hill = cfg.hillslope(n);
    let cell_m = cfg.cell_km(n) * 1000.0;
    // the sea is the base level: held at sea level during the level (D3)
    let mut z = upscaled.clone();
    for k in 0..cells {
        if !land[k] {
            z.data[k] = cfg.sea_level;
        }
    }
    let du: Vec<f32> = (0..cells).map(|k| if land[k] { uplift.data[k] * cfg.dt_yr / n2m } else { 0.0 }).collect();
    let n_land = land.iter().filter(|&&l| l).count().max(1);
    let mean_du = du.iter().map(|&v| v as f64).sum::<f64>() / n_land as f64;
    let (mut s_up, mut s_ero, mut s_dif, mut s_tal) =
        (vec![0f32; cells], vec![0f32; cells], vec![0f32; cells], vec![0f32; cells]);
    let mut secs = [0f64; 4];
    let (mut ratios, mut calm, mut steps, mut equilibrium, mut talus_last) = (Vec::new(), 0usize, 0usize, false, 0f32);
    while steps < max_steps {
        if cancel() {
            return None;
        }
        let start = z.data.clone();
        // 1. uplift
        let t = Instant::now();
        for k in 0..cells {
            z.data[k] += du[k];
            s_up[k] += du[k];
        }
        secs[1] += t.elapsed().as_secs_f64();
        // 2. erosion (the implicit stream power, one iteration, talus and diffusion off inside)
        let t = Instant::now();
        let before = z.data.clone();
        z = incise_with_floor(&z, &solver, None, None, &mut |_, _| {});
        for k in 0..cells {
            s_ero[k] += z.data[k] - before[k];
        }
        secs[2] += t.elapsed().as_secs_f64();
        // 3. diffusion: the talus, then the linear diffusion (production's order)
        let t = Instant::now();
        let before = z.data.clone();
        talus_sweep(&mut z, &hill, cell_m, n2m);
        let mut moved = 0usize;
        for k in 0..cells {
            let d = z.data[k] - before[k];
            if d != 0.0 {
                s_tal[k] += d;
                if land[k] {
                    moved += 1;
                }
            }
        }
        talus_last = moved as f32 / n_land as f32;
        linear_diffusion(&mut z, &hill, None);
        // the sea stays the base level (the talus may have shed onto it)
        for k in 0..cells {
            if !land[k] {
                z.data[k] = cfg.sea_level;
            }
            s_dif[k] += z.data[k] - before[k];
        }
        secs[3] += t.elapsed().as_secs_f64();
        steps += 1;
        let change: f64 = (0..cells).filter(|&k| land[k]).map(|k| (z.data[k] - start[k]).abs() as f64).sum();
        // no uplift: no steady state to reach, the level runs to its cap
        let ratio = if mean_du > 0.0 { (change / n_land as f64 / mean_du) as f32 } else { f32::INFINITY };
        ratios.push(ratio);
        progress(CascadeProgress { level, n, step: steps, max_steps, ratio });
        calm = if ratio < cfg.eq_ratio { calm + 1 } else { 0 };
        if calm >= cfg.eq_window {
            equilibrium = true;
            break;
        }
    }
    for k in 0..cells {
        if !land[k] {
            z.data[k] = upscaled.data[k];
            s_dif[k] = 0.0;
            s_tal[k] = 0.0;
        }
    }
    Some(LevelRecord {
        n,
        cell_km: cfg.cell_km(n),
        d_m2_yr: cfg.diffusivity(n),
        weight: cfg.diffusion_weight(n),
        upscaled,
        uplift,
        land,
        sum_uplift: s_up,
        sum_erosion: s_ero,
        sum_diffusion: s_dif,
        sum_talus: s_tal,
        z,
        steps,
        equilibrium,
        ratios,
        secs,
        talus_last_share: talus_last,
    })
}

/// The viz's « Ombrage » for a normalised field, as RGBA rows **north up** (the field is stored south-first): the
/// hypsometric palette (`hypsometric_bipolar` at sea 0.5) times `0.25 + 0.75 · hillshade` (sun az 315°, alt 45°),
/// exactly the workspace's `ReliefView::Shade` formula, copied so the benches' images and the Cascade window share it.
pub fn shade_rgba(z: &GridF32, cell_m: f32, depth_scale_m: f32) -> Vec<u8> {
    let (w, h) = (z.width, z.height);
    let n2m = 2.0 * 1.13 * depth_scale_m;
    let (az, alt) = (315f32.to_radians(), 45f32.to_radians());
    let l = [alt.cos() * az.sin(), alt.cos() * az.cos(), alt.sin()];
    let mut out = vec![0u8; w * h * 4];
    for y in 0..h {
        let (ym, yp) = (y.saturating_sub(1), (y + 1).min(h - 1));
        for x in 0..w {
            let (xm, xp) = (x.saturating_sub(1), (x + 1).min(w - 1));
            let dzdx = (z.data[y * w + xp] - z.data[y * w + xm]) * n2m / (cell_m * (xp - xm).max(1) as f32);
            let dzdy = (z.data[yp * w + x] - z.data[ym * w + x]) * n2m / (cell_m * (yp - ym).max(1) as f32);
            let nn = (dzdx * dzdx + dzdy * dzdy + 1.0).sqrt();
            let s = ((-dzdx * l[0] - dzdy * l[1] + l[2]) / nn).max(0.0) / l[2];
            let [r, g, b] = hypsometric(z.data[y * w + x]);
            let f = (0.25 + 0.75 * s).min(1.35);
            let o = ((h - 1 - y) * w + x) * 4;
            out[o] = (r as f32 * f).min(255.0) as u8;
            out[o + 1] = (g as f32 * f).min(255.0) as u8;
            out[o + 2] = (b as f32 * f).min(255.0) as u8;
            out[o + 3] = 255;
        }
    }
    out
}

/// A signed difference (m) as RGBA rows north up: blue = lowered, red = raised, saturating at `sat_m`.
pub fn diff_rgba(d_m: &[f32], w: usize, h: usize, sat_m: f32) -> Vec<u8> {
    let mut out = vec![0u8; w * h * 4];
    for y in 0..h {
        for x in 0..w {
            let t = (d_m[y * w + x] / sat_m).clamp(-1.0, 1.0);
            let (r, g, b) = if t >= 0.0 {
                (255.0, 255.0 * (1.0 - t), 255.0 * (1.0 - t))
            } else {
                (255.0 * (1.0 + t), 255.0 * (1.0 + t), 255.0)
            };
            let o = ((h - 1 - y) * w + x) * 4;
            out[o] = r as u8;
            out[o + 1] = g as u8;
            out[o + 2] = b as u8;
            out[o + 3] = 255;
        }
    }
    out
}

/// `ymir-viz`'s `hypsometric_bipolar(t, 0.5)`, copied (the core cannot depend on the viz).
fn hypsometric(t: f32) -> [u8; 3] {
    let sea = 0.5f32;
    let mid = (sea + 1.0) * 0.5;
    let lerp = |t: f32, a: [u8; 3], b: [u8; 3]| -> [u8; 3] {
        let t = t.clamp(0.0, 1.0);
        [
            (a[0] as f32 + t * (b[0] as f32 - a[0] as f32)).round() as u8,
            (a[1] as f32 + t * (b[1] as f32 - a[1] as f32)).round() as u8,
            (a[2] as f32 + t * (b[2] as f32 - a[2] as f32)).round() as u8,
        ]
    };
    if t <= sea * 0.5 {
        lerp(t / (sea * 0.5), [10, 20, 60], [40, 80, 160])
    } else if t <= sea {
        lerp((t - sea * 0.5) / (sea * 0.5), [40, 80, 160], [120, 180, 230])
    } else if t <= mid {
        lerp((t - sea) / (mid - sea), [60, 130, 60], [140, 100, 50])
    } else {
        lerp((t - mid) / (1.0 - mid), [140, 100, 50], [245, 245, 245])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bump(n: usize) -> GridF32 {
        // an island: a Gaussian dome on a sea at 0.45 (norm)
        let mut g = GridF32::new(n, n, 0.45);
        let c = n as f32 / 2.0;
        for y in 0..n {
            for x in 0..n {
                let r2 = ((x as f32 - c).powi(2) + (y as f32 - c).powi(2)) / (n as f32 * 0.22).powi(2);
                g.data[y * n + x] = 0.45 + 0.15 * (-r2).exp();
            }
        }
        g
    }

    /// The bicubic keeps the coarse values on the even pixels and reproduces a cubic exactly (Catmull-Rom is exact on
    /// polynomials up to degree 2 at the midpoint; checked on a periodic sine within its error bound).
    #[test]
    fn the_upsample_keeps_the_coarse_pixels() {
        let g = bump(32);
        let u = upsample2(&g);
        for y in 0..32 {
            for x in 0..32 {
                assert_eq!(u.data[2 * y * 64 + 2 * x], g.data[y * 32 + x]);
            }
        }
        // a periodic sine of wavelength 16 cells: the midpoints within 1 % of the amplitude
        let mut s = GridF32::new(32, 32, 0.0);
        for y in 0..32 {
            for x in 0..32 {
                s.data[y * 32 + x] = (x as f32 * std::f32::consts::TAU / 16.0).sin();
            }
        }
        let su = upsample2(&s);
        for x in 0..64 {
            let want = (x as f32 * 0.5 * std::f32::consts::TAU / 16.0).sin();
            assert!((su.data[5 * 64 + x] - want).abs() < 0.01, "x {x}: {} vs {want}", su.data[5 * 64 + x]);
        }
    }

    /// ADR Finding 159 -- the cascade on a synthetic island: each level reaches its equilibrium, the land mask and the
    /// sea hold, and the decomposition adds up to the result. Negative control: with the uplift off the land only
    /// lowers (no operator lifts the mean).
    #[test]
    fn the_cascade_reaches_equilibrium_and_its_substeps_add_up() {
        let mut cfg = CascadeConfig::declared(5000.0);
        cfg.levels = vec![32, 64];
        cfg.max_steps = vec![600, 400];
        cfg.domain_km = 100.0;
        let coarse = bump(16);
        let mut c = Cascade::new(coarse.clone(), cfg.clone());
        c.run_all(&mut |_| {}, &|| false);
        assert_eq!(c.levels.len(), 2);
        let cal = c.calibration.as_ref().unwrap();
        assert!(cal.equilibrium, "the calibration run did not reach equilibrium in {} steps", cal.steps);
        for l in &c.levels {
            assert!(l.equilibrium, "level {} stopped at its cap ({} steps, last ratio {:?})", l.n, l.steps, l.ratios.last());
            let full = l.substep_field(SubStep::Diffusion);
            let mut sum = l.substep_field(SubStep::Erosion);
            for k in 0..l.n * l.n {
                if l.land[k] {
                    sum.data[k] += l.sum_diffusion[k];
                    assert!((sum.data[k] - full.data[k]).abs() < 1e-4, "the decomposition does not add up at {k}");
                    assert!(full.data[k] >= cfg.sea_level, "a land cell went under the sea at {k}");
                } else {
                    assert_eq!(full.data[k], l.upscaled.data[k], "the bathymetry was not restored at {k}");
                }
            }
        }
        // the calibration hit its target within 10 % (talus inert on this gentle dome)
        let l0 = &c.levels[0];
        let m = mean_land_m(&l0.z, &l0.land, &cfg);
        assert!((m - cal.target_mean_m).abs() < 0.1 * cal.target_mean_m, "mean {m} vs target {}", cal.target_mean_m);
        // the input is untouched
        assert_eq!(c.coarse.data, coarse.data);
        // negative control: no uplift → the mean land altitude can only fall
        let mut off = cfg.clone();
        off.calibrate = false;
        off.u0_m_yr = 0.0;
        off.levels = vec![32];
        off.max_steps = vec![50];
        let mut c0 = Cascade::new(coarse, off.clone());
        c0.run_all(&mut |_| {}, &|| false);
        let l = &c0.levels[0];
        assert!(mean_land_m(&l.z, &l.land, &off) < mean_land_m(&l.upscaled, &l.land, &off));
    }
}
