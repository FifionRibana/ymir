//! ADR Findings 159–161 -- the « Cascade » window: the multi-scale cascade (`ymir_core::cascade`), level by level.
//!
//! F161: N1 by default (the implicit solver at n = 1, the area capped, the smooth retargeting), physics at 256²,
//! « k du niveau » with « Caler sur la Corse », and the « Référence Corse » view (Copernicus DEM GLO-30, produced using
//! Copernicus WorldDEM-30 © DLR e.V. 2010-2014 and © Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS
//! by the European Union and ESA; all rights reserved).
//!
//! The author (2026-10-09): « on aurait clairement les différentes résolutions qui seront balayées, soit d'un seul
//! trait, soit résolution par résolution, par action utilisateur … visualiser les résultats intermédiaires » ; and at F160:
//! « je n'ai rien pour les visualiser sur le viz, ni leurs formes et affluents. Pareil pour les lacs. »
//!
//! A floating window, opened from the top bar, on the world the workspace frames (its seed and its framing roll). It
//! never reads or writes the HD world: its own worker thread runs the 64² tectonics, the PHYSICS level (F160-D1:
//! uplift, the implicit solver, diffusion, to equilibrium, calibrated on the peak), then the AMPLIFICATION levels
//! (Schott 2024 transposed: agrandissement, érosion, talus, dépôt), and the final retargeting.
//!
//! Each level shows its sub-steps and their times; the rivers (Strahler order, the trunks ≥ 100 km² darker) and the
//! lakes (filled to their spill) are drawn on it; the hover reads altitude, drainage area, order and lake. The budgets
//! of the next level are editable, and each process has a button that adds iterations to the selected level on top of
//! its current state (« sur la base de l'existant »); the levels after it are then discarded, and the window says so.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

use bevy::prelude::*;
use bevy_egui::{EguiContexts, EguiPrimaryContextPass, egui};
use crossbeam_channel::{Receiver, Sender, unbounded};
use egui::Color32 as C;
use ymir_core::cascade::amplify::{AmpConfig, AmpLevel, Budget, Chain, Variant, Warp, WarpPlace};
use ymir_core::cascade::hydro::{Hydro, draw, hydrology};
use ymir_core::cascade::{
    CascadeConfig, CascadeProgress, LevelRecord, PeakCalibration, SubStep, diff_rgba, physics_level, roll, shade_m_rgba,
};
use ymir_core::grid::GridF32;
use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::init_r7::init_c1_state_phase_2_r7;
use ymir_core::tectonics_c1::kinematics::PlateKinematics;
use ymir_core::tectonics_c1::production_upscale::c1_coarse_normalized_altitude;
use ymir_core::tectonics_c1::time_loop::{C1TimeLoopConfig, run_with_closures};

use crate::bridge::c1::C1RunSpec;

const COPPER: C = C::from_rgb(0xB8, 0x73, 0x33);
const INK: C = C::from_rgb(0x1A, 0x12, 0x0A);
const DIM: C = C::from_rgb(0x8A, 0x8A, 0x8A);
const GREEN: C = C::from_rgb(0x5F, 0xA5, 0x5F);
const AMBER: C = C::from_rgb(0xC0, 0x9A, 0x4A);
const DOMAIN_KM: f32 = 400.0;
/// F160-D1: the author's « type Corse » peak (the middle of 2 700–3 000 m).
const PEAK_TARGET_M: f32 = 2850.0;
const LAST_LEVEL: usize = 8192;
/// F161: N1's k calibrated on Corsica by `f161_cascade` for the témoin's world (the closest trial at 512² and 1 024²,
/// where the target was out of reach; within ±20 % at 2 048²).
const F161_K: [(usize, f32); 3] = [(512, 1.0e3), (1024, 1.0e3), (2048, 177.8)];

pub struct CascadePlugin;

impl Plugin for CascadePlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<CascadeUi>();
        app.add_systems(EguiPrimaryContextPass, draw_cascade);
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Process {
    Erosion,
    Talus,
    Deposit,
}

impl Process {
    fn label(self) -> &'static str {
        match self {
            Process::Erosion => "érosion",
            Process::Talus => "talus",
            Process::Deposit => "dépôt",
        }
    }
}

enum Cmd {
    /// The physics level at n cells with F163's variants (resets the amplification).
    Physics(usize, PhysOpts),
    /// The amplification's settings: variant, budget scale, k multiplier (resets the amplification, keeps the physics).
    Config(Variant, f32, f32, Remedies),
    /// The next level at this budget (and, for N1, this k).
    Next(Budget, Option<f32>),
    /// F161-V: calibrate N1's k for the next level on Corsica's 2–4 cell octave (the target in m).
    Calibrate(f32),
    /// Levels up to n cells at the declared budgets.
    All(usize),
    /// Add iterations of one process to amplification level i (the levels after it are discarded).
    Add(usize, Process, usize),
    Retarget,
}

/// A computed field, ready to display: its result (m) and its hydrology.
#[derive(Clone)]
struct Shown {
    n: usize,
    result: Arc<GridF32>,
    hydro: Arc<Hydro>,
}

enum Evt {
    Coarse(Shown, f64),
    Physics(Shown, Arc<LevelRecord>, PeakCalibration),
    Level(usize, Shown, Arc<AmpLevel>),
    Truncate(usize),
    Retargeted(Shown, usize),
    /// F161-V: the calibrated k for a level, and its trials (k, octave 0).
    Calibrated(usize, f32, Vec<(f32, f32)>),
    Progress(String, f32),
    Idle,
}

struct Worker {
    tx: Sender<Cmd>,
    rx: Receiver<Evt>,
    cancel: Arc<AtomicBool>,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum View {
    Shade,
    Hypso,
    Diff,
}

/// What the frieze selects: the 64² tectonics, the physics level, an amplification level, the retargeted field.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Sel {
    Coarse,
    Physics,
    Amp(usize),
    Retargeted,
}

/// The window's state. `seed` / `offset_cells` are the workspace's world, set by its top-bar button.
#[derive(Resource)]
pub struct CascadeUi {
    pub open: bool,
    pub seed: u64,
    pub offset_cells: [i64; 2],
    worker: Option<Worker>,
    world: Option<(u64, [i64; 2])>,
    coarse: Option<(Shown, f64)>,
    physics: Option<(Shown, Arc<LevelRecord>, PeakCalibration)>,
    levels: Vec<(Shown, Arc<AmpLevel>)>,
    retargeted: Option<(Shown, usize)>,
    busy: bool,
    progress: Option<(String, f32)>,
    note: Option<String>,
    phys_n: usize,
    variant: Variant,
    budget_scale: f32,
    /// F160 diagnostic A1: the k multiplier (1 = the declaration).
    k_boost: f32,
    /// F162-V: the remedies against the predictability (ρ, π) and the deficit (R4), with their amplitudes.
    rho_on: bool,
    rho_amp: f32,
    pi_on: bool,
    pi_factor: f32,
    r4_on: bool,
    /// F164-V: the warp (none / W-phys / W-tous), A and L in previous-level cells.
    warp_place: Option<WarpPlace>,
    warp_amp: f32,
    warp_corr: f32,
    /// F163-V: the physics level's variants, read when « Calculer la physique » is pressed.
    rr_on: bool,
    rr_tau: f32,
    rough_on: bool,
    rough_m: f32,
    mfd_p: Option<f32>,
    /// F161-V: « k du niveau » (N1) for the next level, and the level it was set for.
    n1_k: f32,
    n1_k_for: usize,
    /// F161-V: the « Référence Corse » view, and its loaded grid (cell count, shown, octave 0).
    corse: bool,
    corse_cache: Option<(usize, egui::TextureHandle)>,
    calib_note: Option<String>,
    next_budget: Budget,
    /// The resolution `next_budget` was last filled for (it follows the declared budget until edited).
    budget_for: usize,
    add_iters: [usize; 3],
    sel: Sel,
    sub: usize,
    view: View,
    rivers: bool,
    lakes: bool,
    min_order: u8,
    texture: Option<egui::TextureHandle>,
    tex_key: Option<String>,
    diff_sat: Option<f32>,
}

impl Default for CascadeUi {
    fn default() -> Self {
        Self {
            open: false,
            seed: 0,
            offset_cells: [0, 0],
            worker: None,
            world: None,
            coarse: None,
            physics: None,
            levels: Vec::new(),
            retargeted: None,
            busy: false,
            progress: None,
            note: None,
            phys_n: 256,
            variant: Variant::N1,
            budget_scale: 1.0,
            k_boost: 1.0,
            rho_on: false,
            rho_amp: 0.4,
            pi_on: false,
            pi_factor: 0.15,
            r4_on: false,
            warp_place: None,
            warp_amp: 0.5,
            warp_corr: 4.0,
            rr_on: false,
            rr_tau: 0.5,
            rough_on: false,
            rough_m: 42.4,
            mfd_p: Some(2.0),
            n1_k: 100.0,
            n1_k_for: 0,
            corse: true,
            corse_cache: None,
            calib_note: None,
            next_budget: Budget::moyen(256),
            budget_for: 0,
            add_iters: [100, 200, 50],
            sel: Sel::Coarse,
            sub: 3,
            view: View::Shade,
            rivers: true,
            lakes: true,
            min_order: 2,
            texture: None,
            tex_key: None,
            diff_sat: None,
        }
    }
}

/// The 64² normalised C1 altitude of `seed`, rolled to the framing: the tectonic run of `run_hd` (`hd_run_config`, the
/// default spec) and production's `target_land_fraction` (`None`).
pub fn cascade_coarse(seed: u64, offset_cells: [i64; 2]) -> GridF32 {
    let spec = C1RunSpec { seed, ..C1RunSpec::default() };
    let run = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: spec.n_steps,
        dx: 1.0 / spec.grid_size as f64,
        dy: 1.0 / spec.grid_size as f64,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: spec.drainage_max_distance,
    };
    let mut state = init_c1_state_phase_2_r7(spec.grid_size, seed, &spec.init_params);
    let mut kin = PlateKinematics::preset_phase_1_1(state.num_plates);
    run_with_closures(&mut state, &mut kin, &run, &spec.closures, |_, _| {});
    let ss = SteinSteinParams::default();
    let coarse = c1_coarse_normalized_altitude(&state, &run.iso_config, &ss, None);
    let g = coarse.width as i64;
    roll(&coarse, [offset_cells[0].rem_euclid(g) as usize, offset_cells[1].rem_euclid(g) as usize])
}

fn n2m() -> f32 {
    2.0 * 1.13 * SteinSteinParams::default().depth_scale_m as f32
}
fn to_m(g: &GridF32) -> GridF32 {
    let k = n2m();
    GridF32 { width: g.width, height: g.height, data: g.data.iter().map(|v| (v - 0.5) * k).collect() }
}
fn shown(g: GridF32) -> Shown {
    let n = g.width;
    let hy = hydrology(&g, DOMAIN_KM / n as f32);
    Shown { n, result: Arc::new(g), hydro: Arc::new(hy) }
}

fn spawn(seed: u64, offset_cells: [i64; 2], ccfg: CascadeConfig) -> Worker {
    let (tx, cmd_rx) = unbounded::<Cmd>();
    let (etx, rx) = unbounded::<Evt>();
    let cancel = Arc::new(AtomicBool::new(false));
    let stop = cancel.clone();
    std::thread::Builder::new()
        .name("ymir-cascade".into())
        .spawn(move || {
            let t = Instant::now();
            let coarse = cascade_coarse(seed, offset_cells);
            let _ = etx.send(Evt::Coarse(shown(to_m(&coarse)), t.elapsed().as_secs_f64()));
            // F161: N1 by default, as the window
            let (mut variant, mut scale, mut boost) = (Variant::N1, 1.0f32, 1.0f32);
            // F162: the talus at a quarter from 1 024² by default (decided at F161's commit); ρ, π, R4 as ticked
            let amp_cfg = |v: Variant, f: f32, b: f32, r: &Remedies| AmpConfig {
                k_boost: b,
                talus_fine_scale: 0.25,
                n1_rho: r.rho,
                n1_pi_m: r.pi_m.clone(),
                n1_recal_depth: if r.r4 { 2 } else { 1 },
                warp: r.warp,
                ..AmpConfig::declared(seed, v, f)
            };
            let mut remedies = Remedies::default();
            let mut chain: Option<Chain> = None;
            let send_level = |etx: &Sender<Evt>, ch: &Chain, i: usize| {
                let l = &ch.levels[i];
                let _ = etx.send(Evt::Level(i, shown(l.result()), Arc::new(l.clone())));
            };
            for cmd in cmd_rx.iter() {
                stop.store(false, Ordering::Relaxed);
                let cancel = || stop.load(Ordering::Relaxed);
                match cmd {
                    Cmd::Physics(n, o) => {
                        let cap = if n <= 128 { 400 } else { 300 };
                        // F163-V: the physics level's variants (the random receiver, the roughness, the MFD exponent)
                        let ccfg = CascadeConfig {
                            receiver_tau: o.tau,
                            roughness_m: o.rough_m,
                            mfd_exponent: o.mfd,
                            seed,
                            ..ccfg.clone()
                        };
                        let ptx = etx.clone();
                        let mut prog = |p: CascadeProgress| {
                            let what = if p.level.is_none() { "physique, essai" } else { "physique, calée" };
                            let text = format!("{what} {}² · pas {} / {} · |Δz| / U·dt = {:.4}", p.n, p.step, p.max_steps, p.ratio);
                            drop(ptx.send(Evt::Progress(text, p.step as f32 / p.max_steps as f32)));
                        };
                        if let Some((rec, cal)) = physics_level(&coarse, n, &ccfg, PEAK_TARGET_M, cap, &mut prog, &cancel) {
                            let pm = to_m(&rec.z);
                            chain = Some(Chain::new(pm.clone(), amp_cfg(variant, scale, boost, &remedies)));
                            let _ = etx.send(Evt::Truncate(0));
                            let _ = etx.send(Evt::Physics(shown(pm), Arc::new(rec), cal));
                        }
                    }
                    Cmd::Config(v, f, b, r) => {
                        (variant, scale, boost) = (v, f, b);
                        remedies = r;
                        if let Some(ch) = chain.as_mut() {
                            ch.cfg = amp_cfg(v, f, b, &remedies);
                            ch.levels.clear();
                            let _ = etx.send(Evt::Truncate(0));
                        }
                    }
                    Cmd::Calibrate(target) => {
                        if let Some(ch) = chain.as_mut() {
                            let n = ch.last_result().width * 2;
                            let _ = etx.send(Evt::Progress(format!("{n}² : calage de k sur la Corse (≤ 8 essais)"), 0.0));
                            // F162-D1: the octave is not monotonic in k; calibrate on its peak
                            let (k, trials) = ch.calibrate_peak(target, &cancel);
                            let _ = etx.send(Evt::Calibrated(n, k, trials));
                        }
                    }
                    Cmd::Next(b, k) => {
                        if let Some(ch) = chain.as_mut() {
                            if let Some(k) = k {
                                let n = ch.last_result().width * 2;
                                ch.cfg.set_n1_k(n, k);
                            }
                            if run_level(ch, b, &etx, &cancel) {
                                send_level(&etx, ch, ch.levels.len() - 1);
                            }
                        }
                    }
                    Cmd::All(top) => {
                        if let Some(ch) = chain.as_mut() {
                            while ch.last_result().width < top && !cancel() {
                                let b = ch.cfg.budget(ch.last_result().width * 2);
                                if !run_level(ch, b, &etx, &cancel) {
                                    break;
                                }
                                send_level(&etx, ch, ch.levels.len() - 1);
                            }
                        }
                    }
                    Cmd::Add(i, p, iters) => {
                        if let Some(ch) = chain.as_mut().filter(|c| i < c.levels.len()) {
                            ch.levels.truncate(i + 1);
                            let _ = etx.send(Evt::Truncate(i + 1));
                            let cfg = ch.cfg.clone();
                            let l = &mut ch.levels[i];
                            let _ = etx.send(Evt::Progress(format!("{}² : + {iters} itérations de {}", l.n, p.label()), 0.0));
                            match p {
                                Process::Erosion => l.erode(iters, &cfg, &cancel),
                                Process::Talus => l.talus(iters, &cfg, &cancel),
                                Process::Deposit => l.deposit(iters, &cfg, &cancel),
                            }
                            send_level(&etx, ch, i);
                        }
                    }
                    Cmd::Retarget => {
                        if let Some((g, nc)) = chain.as_ref().and_then(|c| c.retargeted()) {
                            let _ = etx.send(Evt::Retargeted(shown(g), nc));
                        }
                    }
                }
                let _ = etx.send(Evt::Idle);
            }
        })
        .expect("spawn the cascade worker");
    Worker { tx, rx, cancel }
}

/// One amplification level at `b`: U, then E, T, D (and S+U's 2×2 retargeting), with progress. `false` if cancelled.
fn run_level(ch: &mut Chain, b: Budget, etx: &Sender<Evt>, cancel: &dyn Fn() -> bool) -> bool {
    let prev = ch.last_result();
    let cfg = ch.cfg.clone();
    let n = prev.width * 2;
    let _ = etx.send(Evt::Progress(format!("{n}² : agrandissement"), 0.0));
    let total = (b.erosion + b.talus + b.deposit).max(1) as f32;
    let mut done = 0usize;
    let l = ch.start_level();
    for (p, count) in [(Process::Erosion, b.erosion), (Process::Talus, b.talus), (Process::Deposit, b.deposit)] {
        let mut left = count;
        while left > 0 {
            if cancel() {
                ch.levels.pop();
                return false;
            }
            let k = left.min(50);
            match p {
                Process::Erosion => l.erode(k, &cfg, cancel),
                Process::Talus => l.talus(k, &cfg, cancel),
                Process::Deposit => l.deposit(k, &cfg, cancel),
            }
            left -= k;
            done += k;
            let text = format!("{n}² : {} {}/{count}", p.label(), count - left);
            let _ = etx.send(Evt::Progress(text, done as f32 / total));
        }
    }
    if cfg.variant == Variant::SU {
        l.recalage_2x2(&prev, &cfg);
    }
    true
}

impl CascadeUi {
    fn start(&mut self) {
        let ss = SteinSteinParams::default();
        if let Some(w) = &self.worker {
            w.cancel.store(true, Ordering::Relaxed);
        }
        self.worker = Some(spawn(self.seed, self.offset_cells, CascadeConfig::declared(ss.depth_scale_m as f32)));
        self.world = Some((self.seed, self.offset_cells));
        self.coarse = None;
        self.physics = None;
        self.levels.clear();
        self.retargeted = None;
        self.busy = true;
        self.progress = None;
        self.note = None;
        self.sel = Sel::Coarse;
        self.tex_key = None;
    }

    fn drain(&mut self) {
        let Some(w) = &self.worker else { return };
        let evts: Vec<Evt> = w.rx.try_iter().collect();
        for e in evts {
            match e {
                Evt::Coarse(s, t) => {
                    self.coarse = Some((s, t));
                    self.busy = false;
                }
                Evt::Physics(s, r, c) => {
                    self.physics = Some((s, r, c));
                    self.sel = Sel::Physics;
                    self.note = None;
                }
                Evt::Level(i, s, l) => {
                    if i < self.levels.len() {
                        self.levels[i] = (s, l);
                    } else {
                        self.levels.push((s, l));
                    }
                    self.sel = Sel::Amp(i);
                    self.sub = 3;
                    self.retargeted = None;
                }
                Evt::Truncate(len) => {
                    if self.levels.len() > len {
                        self.note = Some(format!(
                            "les niveaux au-delà de {}² sont écartés : « Niveau suivant » les recalcule sur cette base",
                            self.level_n(len) / 2
                        ));
                    }
                    self.levels.truncate(len);
                    self.retargeted = None;
                    if let Sel::Amp(i) = self.sel {
                        if i >= len {
                            self.sel = if len == 0 { Sel::Physics } else { Sel::Amp(len - 1) };
                        }
                    }
                }
                Evt::Retargeted(s, nc) => {
                    self.retargeted = Some((s, nc));
                    self.sel = Sel::Retargeted;
                }
                Evt::Calibrated(n, k, trials) => {
                    self.n1_k = k;
                    self.n1_k_for = n;
                    self.calib_note = Some(format!(
                        "{n}² : k calé = {k:.3e} en {} essais ({})",
                        trials.len(),
                        trials.iter().map(|(k, v)| format!("{k:.2e} → {v:.1} m")).collect::<Vec<_>>().join(", ")
                    ));
                }
                Evt::Progress(t, f) => self.progress = Some((t, f)),
                Evt::Idle => {
                    self.busy = false;
                    self.progress = None;
                }
            }
            self.tex_key = None;
        }
    }

    fn send(&mut self, c: Cmd) {
        if let Some(w) = &self.worker {
            if w.tx.send(c).is_ok() {
                self.busy = true;
            }
        }
    }

    /// The resolution of amplification level i (0-based), given the physics level.
    fn level_n(&self, i: usize) -> usize {
        self.physics.as_ref().map_or(self.phys_n, |p| p.0.n) << (i + 1)
    }

    fn shown(&self, sel: Sel) -> Option<&Shown> {
        match sel {
            Sel::Coarse => self.coarse.as_ref().map(|c| &c.0),
            Sel::Physics => self.physics.as_ref().map(|p| &p.0),
            Sel::Amp(i) => self.levels.get(i).map(|l| &l.0),
            Sel::Retargeted => self.retargeted.as_ref().map(|r| &r.0),
        }
    }
}

fn draw_cascade(mut contexts: EguiContexts, mut cu: ResMut<CascadeUi>) {
    if !cu.open {
        return;
    }
    let Ok(ctx) = contexts.ctx_mut() else { return };
    cu.drain();
    if cu.busy {
        ctx.request_repaint();
    }
    let mut open = true;
    egui::Window::new("Cascade multi-échelle (F159–F161)")
        .open(&mut open)
        .default_size([860.0, 820.0])
        .resizable(true)
        .show(ctx, |ui| body(ui, &mut cu));
    cu.open = open;
}

fn button(ui: &mut egui::Ui, text: &str, enabled: bool) -> bool {
    let b = egui::Button::new(egui::RichText::new(text).color(INK).strong().size(12.0)).fill(COPPER).corner_radius(5.0);
    ui.add_enabled(enabled, b).clicked()
}

fn small(text: String) -> egui::RichText {
    egui::RichText::new(text).color(DIM).size(11.0)
}

fn body(ui: &mut egui::Ui, cu: &mut CascadeUi) {
    ui.label(small(format!(
        "Monde : graine {} · cadrage ({}, {}) · 400 km · prototype, la production n'est pas touchée",
        cu.seed, cu.offset_cells[0], cu.offset_cells[1]
    )));
    let stale = cu.world.is_some_and(|w| w != (cu.seed, cu.offset_cells));
    let ready = cu.worker.is_some() && cu.coarse.is_some() && !cu.busy && !stale;
    // ── 1. the world and the physics level ──
    ui.horizontal(|ui| {
        if cu.worker.is_none() || stale {
            let label = if stale { "Recommencer sur le monde courant" } else { "Démarrer (tectonique 64²)" };
            if button(ui, label, !cu.busy || stale) {
                cu.start();
            }
        }
        ui.label(small("Niveau de physique :".into()));
        for n in [128usize, 256] {
            ui.selectable_value(&mut cu.phys_n, n, format!("{n}²"));
        }
        if button(ui, "Calculer la physique", ready) {
            let n = cu.phys_n;
            let o = phys_opts(cu);
            cu.send(Cmd::Physics(n, o));
        }
        if cu.busy && cu.coarse.is_some() && ui.button("Annuler").clicked() {
            if let Some(w) = &cu.worker {
                w.cancel.store(true, Ordering::Relaxed);
            }
        }
    });
    // F163-V: the physics level's variants (taken when « Calculer la physique » is pressed)
    ui.horizontal(|ui| {
        ui.label(small("Variantes de la physique :".into()));
        ui.checkbox(&mut cu.rr_on, "récepteur aléatoire").on_hover_text(
            "P-mfd, F163-D : à chaque pas, le récepteur est tiré parmi les voisins plus bas dont la pente ≥ τ × la plus forte, avec une probabilité ∝ la pente ; les plats gardent le D8",
        );
        ui.add_enabled(cu.rr_on, egui::DragValue::new(&mut cu.rr_tau).range(0.05..=1.0).speed(0.01).prefix("τ "));
        ui.checkbox(&mut cu.rough_on, "rugosité initiale").on_hover_text(
            "P-bruit, F163-D : un bruit de bande (8 et 4 cellules) ajouté aux terres avant la physique ; 42,4 m = 0,5 × l'octave 3–6 km de la Corse à 1 563 m",
        );
        ui.add_enabled(cu.rough_on, egui::DragValue::new(&mut cu.rough_m).range(0.0..=300.0).speed(0.5).suffix(" m"));
        ui.separator();
        ui.label(small("MFD p :".into())).on_hover_text("P-dissection, F163-D : l'exposant de l'aire MFD de l'incision (2 = F159–F162 ; D8 = sans MFD)");
        for (name, p) in [("1", Some(1.0f32)), ("1,5", Some(1.5)), ("2", Some(2.0)), ("3", Some(3.0)), ("4", Some(4.0)), ("6", Some(6.0)), ("D8", None)] {
            ui.selectable_value(&mut cu.mfd_p, p, name);
        }
    });
    // ── 2. the amplification's settings and the next level's budget ──
    ui.horizontal(|ui| {
        ui.label(small("Amplification :".into()));
        let (mut v, mut f) = (cu.variant, cu.budget_scale);
        for x in [Variant::N1, Variant::N1NoCap, Variant::N1NoRecal, Variant::S, Variant::SRho, Variant::SU] {
            ui.selectable_value(&mut v, x, x.label());
        }
        ui.separator();
        for (name, x) in [("faible", 0.5f32), ("moyen", 1.0), ("fort", 2.0)] {
            ui.selectable_value(&mut f, x, name);
        }
        ui.separator();
        let mut b = cu.k_boost;
        ui.add(egui::DragValue::new(&mut b).range(0.1..=1000.0).speed(0.5).prefix("intensité k × "))
            .on_hover_text("Diagnostic F160-A1 : 1 = la déclaration ; ×10 creuse peu, ×100 rabote le relief (mesuré)");
        if (v, f, b) != (cu.variant, cu.budget_scale, cu.k_boost) {
            (cu.variant, cu.budget_scale, cu.k_boost) = (v, f, b);
            cu.budget_for = 0;
            let r = remedies(cu);
            cu.send(Cmd::Config(v, f, b, r));
        }
    });
    // F162-V: ρ, π and R4, ticked with their amplitudes (N1 only); a change re-sends the settings (the levels reset)
    ui.horizontal(|ui| {
        ui.label(small("Contre la prévisibilité :".into()));
        let before = (cu.rho_on, cu.rho_amp, cu.pi_on, cu.pi_factor, cu.r4_on);
        ui.checkbox(&mut cu.rho_on, "ρ dureté bruitée").on_hover_text("k(p) = k·(1 − ρ)/(1 − 0,5), ρ = 0,5 + a·bruit fractal (25 km → 2 cellules), F162-D2");
        ui.add_enabled(cu.rho_on, egui::DragValue::new(&mut cu.rho_amp).range(0.0..=0.49).speed(0.01).prefix("a "));
        ui.checkbox(&mut cu.pi_on, "π perturbation").on_hover_text("Un bruit à l'octave nouvelle (4 cellules) ajouté avant l'érosion, RMS = facteur × l'octave 0 de la Corse, F162-D2");
        ui.add_enabled(cu.pi_on, egui::DragValue::new(&mut cu.pi_factor).range(0.0..=1.0).speed(0.01).prefix("× "));
        ui.checkbox(&mut cu.r4_on, "R4 recalage 4×4").on_hover_text("Le recalage lisse contre le niveau n − 2 : chaque niveau amplifie aussi l'octave du dessus, F162-D3");
        if before != (cu.rho_on, cu.rho_amp, cu.pi_on, cu.pi_factor, cu.r4_on) {
            let r = remedies(cu);
            let (v, f, b) = (cu.variant, cu.budget_scale, cu.k_boost);
            cu.budget_for = 0;
            cu.send(Cmd::Config(v, f, b, r));
        }
    });
    // F164-V: the warp of the upscaled field; a change re-sends the settings (the levels reset)
    ui.horizontal(|ui| {
        ui.label(small("Déformation :".into())).on_hover_text(
            "F164-D1 : le champ agrandi est rééchantillonné en x + d(x), d lisse (|d| ≤ A, longueur L, en cellules du niveau précédent) ; le recalage se fait sur le champ déformé. PROXY de l'hétérogénéité géologique absente",
        );
        let before = (cu.warp_place, cu.warp_amp, cu.warp_corr);
        for (name, w) in [("aucune", None), ("W-phys", Some(WarpPlace::Phys)), ("W-tous", Some(WarpPlace::All))] {
            ui.selectable_value(&mut cu.warp_place, w, name);
        }
        let on = cu.warp_place.is_some();
        ui.add_enabled(on, egui::DragValue::new(&mut cu.warp_amp).range(0.0..=2.0).speed(0.01).prefix("A "));
        ui.add_enabled(on, egui::DragValue::new(&mut cu.warp_corr).range(1.0..=16.0).speed(0.1).prefix("L "));
        if before != (cu.warp_place, cu.warp_amp, cu.warp_corr) {
            let r = remedies(cu);
            let (v, f, b) = (cu.variant, cu.budget_scale, cu.k_boost);
            cu.budget_for = 0;
            cu.send(Cmd::Config(v, f, b, r));
        }
    });
    let have_phys = cu.physics.is_some();
    let next_n = cu.level_n(cu.levels.len());
    if cu.budget_for != next_n {
        cu.next_budget = Budget::moyen(next_n).scaled(cu.budget_scale);
        cu.budget_for = next_n;
    }
    if cu.n1_k_for != next_n {
        // F161: the k the bench calibrated on Corsica for the témoin's world (PSEED, 256² physics); another world
        // recalibrates with « Caler sur la Corse »
        cu.n1_k = F161_K.iter().find(|(n, _)| *n == next_n).map_or(100.0, |p| p.1);
        cu.n1_k_for = next_n;
    }
    ui.horizontal(|ui| {
        ui.label(small(format!("Budget du prochain niveau ({next_n}²) :")));
        ui.add(egui::DragValue::new(&mut cu.next_budget.erosion).prefix("érosion ").range(0..=20000));
        ui.add(egui::DragValue::new(&mut cu.next_budget.talus).prefix("talus ").range(0..=20000));
        ui.add(egui::DragValue::new(&mut cu.next_budget.deposit).prefix("dépôt ").range(0..=20000));
        if ui.small_button("déclaré").on_hover_text("Le budget déclaré (F160-D3) × faible / moyen / fort").clicked() {
            cu.next_budget = Budget::moyen(next_n).scaled(cu.budget_scale);
        }
    });
    if cu.variant.is_n1() {
        ui.horizontal(|ui| {
            ui.label(small(format!("k du niveau ({next_n}², N1) :")));
            let speed = cu.n1_k * 0.02;
            ui.add(egui::DragValue::new(&mut cu.n1_k).range(1e-3..=1e8).speed(speed).custom_formatter(|v, _| format!("{v:.3e}")));
            let target = corse_target(next_n);
            let can = ready && have_phys && target.is_some();
            let hint = match target {
                Some(t) => format!("Cale k pour que l'octave 2–4 cellules atteigne celle de la Corse ({t:.1} m) à ±20 % (≤ 8 essais, F161-D3)"),
                None => "Données Corse absentes (data/corsica, prep_corse.py)".to_string(),
            };
            if ui.add_enabled(can, egui::Button::new("Caler sur la Corse")).on_hover_text(hint).clicked() {
                if let Some(t) = target {
                    cu.send(Cmd::Calibrate(t));
                }
            }
        });
        if let Some(n) = &cu.calib_note {
            ui.label(small(n.clone()));
        }
    }
    ui.horizontal(|ui| {
        let can = ready && have_phys && next_n <= LAST_LEVEL;
        if button(ui, "Niveau suivant", can) {
            let b = cu.next_budget;
            let k = cu.variant.is_n1().then_some(cu.n1_k);
            cu.send(Cmd::Next(b, k));
        }
        if button(ui, "Tout (→ 1024²)", can && next_n <= 1024) {
            cu.send(Cmd::All(1024));
        }
        if button(ui, "Jusqu'à 2048²", can) {
            cu.send(Cmd::All(LAST_LEVEL));
        }
        if button(ui, "Recalage final", ready && !cu.levels.is_empty()) {
            cu.send(Cmd::Retarget);
        }
    });
    if let Some((t, f)) = &cu.progress {
        ui.add(egui::ProgressBar::new(*f).text(t.clone()));
    } else if cu.busy && cu.coarse.is_none() {
        ui.label("tectonique 64² en cours…");
    }
    if let Some(n) = &cu.note {
        ui.label(egui::RichText::new(n).color(AMBER).size(11.0));
    }
    ui.separator();
    // ── 3. the frieze ──
    let mut nodes: Vec<(Sel, String, String, C, bool)> = Vec::new();
    nodes.push((
        Sel::Coarse,
        "64² tectonique".into(),
        cu.coarse.as_ref().map_or("—".into(), |c| format!("✓ {:.1} s", c.1)),
        if cu.coarse.is_some() { GREEN } else { DIM },
        cu.coarse.is_some(),
    ));
    nodes.push((
        Sel::Physics,
        format!("{}² physique", cu.physics.as_ref().map_or(cu.phys_n, |p| p.0.n)),
        cu.physics.as_ref().map_or("—".into(), |p| format!("{} {} pas · {:.1} s", if p.1.equilibrium { "✓" } else { "⚠" }, p.1.steps, p.2.secs)),
        if have_phys { GREEN } else { DIM },
        have_phys,
    ));
    for (i, (s, l)) in cu.levels.iter().enumerate() {
        nodes.push((Sel::Amp(i), format!("{}² {}", s.n, cu.variant.label()), format!("✓ {:.1} s", l.secs.iter().sum::<f64>()), GREEN, true));
    }
    if let Some((s, _)) = &cu.retargeted {
        nodes.push((Sel::Retargeted, format!("{}² recalé", s.n), "✓".into(), GREEN, true));
    }
    egui::ScrollArea::horizontal().id_salt("cascade_frieze").show(ui, |ui| {
        ui.horizontal(|ui| {
            for (i, (sel, title, status, col, avail)) in nodes.into_iter().enumerate() {
                if i > 0 {
                    ui.label(egui::RichText::new("→").color(DIM));
                }
                let on = cu.sel == sel;
                let text = egui::RichText::new(format!("{title}\n{status}")).size(11.0).color(if on { INK } else { col });
                let b = egui::Button::new(text).fill(if on { COPPER } else { C::TRANSPARENT }).min_size(egui::vec2(118.0, 36.0));
                if ui.add_enabled(avail, b).clicked() {
                    cu.sel = sel;
                    cu.sub = 3;
                }
            }
        });
    });
    // ── 4. the selected level ──
    match cu.sel {
        Sel::Physics => {
            if let Some((_, r, c)) = cu.physics.clone() {
                ui.horizontal(|ui| {
                    for (j, s) in SubStep::ALL.iter().enumerate() {
                        let on = cu.sub == j;
                        let text = egui::RichText::new(format!("{} · {:.2} s", s.label(), r.secs[j])).size(11.0).color(if on { INK } else { DIM });
                        if ui.add(egui::Button::new(text).fill(if on { COPPER } else { C::TRANSPARENT })).clicked() {
                            cu.sub = j;
                        }
                    }
                });
                ui.label(small(format!(
                    "calage sur le point culminant {:.0} m : essai U₀ {:.1e} → {:.0} m, U₀ = {:.3e} m/an par km de h_iso · {} pas ({}) · D {:.3} m²/an · \
                     {} cellules intérieures sous 0 m rendues à la terre (D5, profondeur moyenne {:.0} m)",
                    c.target_peak_m,
                    c.u0_trial,
                    c.trial_peak_m,
                    c.u0,
                    r.steps,
                    if r.equilibrium { "équilibre" } else { "plafond" },
                    r.d_m2_yr,
                    c.lifted,
                    c.lifted_mean_depth_m
                )));
            }
        }
        Sel::Amp(i) => {
            if let Some((_, l)) = cu.levels.get(i).cloned() {
                ui.horizontal(|ui| {
                    let names = ["Agrandissement", "Érosion", "Talus", "Dépôt"];
                    let its = [0, l.done.erosion, l.done.talus, l.done.deposit];
                    for j in 0..4 {
                        let on = cu.sub == j;
                        let label = if j == 0 {
                            format!("{} · {:.2} s", names[j], l.secs[j])
                        } else {
                            format!("{} · {} it · {:.1} s", names[j], its[j], l.secs[j])
                        };
                        let text = egui::RichText::new(label).size(11.0).color(if on { INK } else { DIM });
                        if ui.add(egui::Button::new(text).fill(if on { COPPER } else { C::TRANSPARENT })).clicked() {
                            cu.sub = j;
                        }
                    }
                });
                ui.horizontal(|ui| {
                    ui.label(small("Ajouter à ce niveau :".into()));
                    for (j, p) in [Process::Erosion, Process::Talus, Process::Deposit].into_iter().enumerate() {
                        ui.add(egui::DragValue::new(&mut cu.add_iters[j]).range(1..=20000));
                        let b = ui
                            .add_enabled(ready, egui::Button::new(format!("+ {}", p.label())))
                            .on_hover_text("Sur la base de l'existant ; les niveaux suivants sont écartés");
                        if b.clicked() {
                            let k = cu.add_iters[j];
                            cu.send(Cmd::Add(i, p, k));
                        }
                    }
                });
                let c2 = (l.cell_km * l.cell_km) as f64;
                let carved: f64 =
                    l.sum_erosion.iter().zip(&l.ocean).filter(|(_, o)| !**o).map(|(v, _)| -*v as f64).sum::<f64>() * c2 * 1e-3;
                ui.label(small(format!(
                    "cellule {:.3} km · k {:.2e} · pente de référence (p90) {:.4} · d_L {:.0} m · creusé {:.2} km³ · s ≥ s_max sur {:.2} % des \
                     mises à jour · D5 : {} cellules sous 0 m rendues à la terre",
                    l.cell_km,
                    l.k,
                    l.s_ref,
                    l.depth_m,
                    carved,
                    100.0 * l.smax_bind.0 as f64 / l.smax_bind.1.max(1) as f64,
                    l.lifted
                )));
            }
        }
        Sel::Retargeted => {
            if let Some((_, nc)) = &cu.retargeted {
                ui.label(small(format!(
                    "recalage final (Schott 2024 §5.1) : {nc} cellules de crête (aire < 2 cellules) remises sur le niveau de physique, erreur diffusée \
                     sur 500 itérations"
                )));
            }
        }
        Sel::Coarse => {}
    }
    // ── 5. the view ──
    ui.horizontal(|ui| {
        for (v, name) in [(View::Shade, "Ombrage"), (View::Hypso, "Hypsométrie"), (View::Diff, "Différence")] {
            if ui.selectable_label(cu.view == v, name).clicked() {
                cu.view = v;
            }
        }
        ui.separator();
        ui.checkbox(&mut cu.rivers, "Rivières");
        ui.add(egui::DragValue::new(&mut cu.min_order).range(1..=8).prefix("ordre ≥ "));
        ui.checkbox(&mut cu.lakes, "Lacs");
        ui.separator();
        ui.checkbox(&mut cu.corse, "Référence Corse")
            .on_hover_text("Le relief de la Corse (Copernicus DEM GLO-30) à la taille de cellule du niveau affiché, même ombrage");
    });
    let key = format!(
        "{:?}/{}/{:?}/{}/{}/{}/{}/{}",
        cu.sel,
        cu.sub,
        cu.view,
        cu.rivers,
        cu.lakes,
        cu.min_order,
        cu.levels.len(),
        cu.retargeted.is_some()
    );
    if cu.tex_key.as_deref() != Some(key.as_str()) {
        if let Some(img) = build_image(cu) {
            cu.texture = Some(ui.ctx().load_texture("cascade", img, egui::TextureOptions::NEAREST));
            cu.tex_key = Some(key);
        }
    }
    let shown_n = cu.shown(cu.sel).map(|s| s.n);
    let corse_tex = match (cu.corse, shown_n) {
        (true, Some(n)) => corse_texture(ui.ctx(), &mut cu.corse_cache, n),
        _ => None,
    };
    if let (Some(tex), Some(s)) = (&cu.texture, cu.shown(cu.sel)) {
        let room = if corse_tex.is_some() { ui.available_width() / 1.5 } else { ui.available_width() };
        let side = room.min(ui.available_height() - 36.0).max(160.0);
        let resp = ui
            .horizontal(|ui| {
                let resp = ui.add(egui::Image::new((tex.id(), egui::vec2(side, side))).sense(egui::Sense::hover()));
                if let Some(ct) = &corse_tex {
                    // 200 km against our 400 km: half the side, the same km per pixel
                    ui.add(egui::Image::new((ct.id(), egui::vec2(side / 2.0, side / 2.0))));
                }
                resp
            })
            .inner;
        let n = s.n;
        if let Some(pos) = resp.hover_pos() {
            let r = resp.rect;
            let (u, v) = ((pos.x - r.left()) / r.width(), (pos.y - r.top()) / r.height());
            let x = ((u * n as f32) as usize).min(n - 1);
            let y = n - 1 - ((v * n as f32) as usize).min(n - 1); // north up
            let k = y * n + x;
            let hy = &s.hydro;
            let lake = match hy.lake_id[k] {
                0 => "pas de lac".to_string(),
                id => {
                    let l = &hy.lakes[id as usize - 1];
                    format!("lac n° {id} : {:.1} km², déversoir {:.0} m", l.area_km2, l.spill_m)
                }
            };
            let cell = DOMAIN_KM / n as f32;
            let order = if hy.strahler[k] == 0 { "—".to_string() } else { hy.strahler[k].to_string() };
            ui.label(
                small(format!(
                    "({x}, {y}) · {:.1} km, {:.1} km · altitude {:.0} m · aire drainée {:.1} km² · Strahler {order} · {lake}",
                    x as f32 * cell,
                    y as f32 * cell,
                    s.result.data[k],
                    hy.acc_km2[k]
                ))
                .monospace(),
            );
        }
        let mut legend =
            format!("{} lacs · rivières : aire ≥ {:.1} km², troncs ≥ 100 km² en bleu foncé", s.hydro.lakes.len(), s.hydro.river_min_km2);
        if let (View::Diff, Some(sat)) = (cu.view, cu.diff_sat) {
            let what = match (cu.sel, cu.sub) {
                (Sel::Retargeted, _) => "recalé − dernier niveau",
                (_, 0) => "changement total du niveau",
                (Sel::Physics, 1) => "Σ soulèvement",
                (_, 1) => "Σ érosion",
                (Sel::Physics, 2) => "Σ érosion",
                (_, 2) => "Σ talus",
                (Sel::Physics, _) => "Σ diffusion (talus + diffusion)",
                _ => "Σ dépôt",
            };
            legend.push_str(&format!(" · {what} : bleu = abaissé, rouge = relevé, saturé à {sat:.0} m (p98)"));
        }
        if cu.sub < 3 && matches!(cu.sel, Sel::Amp(_) | Sel::Physics) {
            legend.push_str(" · rivières et lacs : ceux du résultat du niveau");
        }
        ui.label(small(legend));
    }
}

/// The texture of the selection: Ombrage / hypsometry / difference, north up, with the rivers and the lakes.
fn build_image(cu: &mut CascadeUi) -> Option<egui::ColorImage> {
    let depth = SteinSteinParams::default().depth_scale_m as f32;
    let s = cu.shown(cu.sel)?.clone();
    let n = s.n;
    let cell_m = DOMAIN_KM / n as f32 * 1000.0;
    // the sub-step's field, and what the « Différence » view shows: the level's total change on the first sub-step,
    // then each process's own Σ
    let (field, diff): (GridF32, Option<Vec<f32>>) = match cu.sel {
        Sel::Amp(i) => {
            let l = &cu.levels.get(i)?.1;
            let d = match cu.sub {
                0 => s.result.data.iter().zip(&l.upscaled.data).map(|(a, b)| a - b).collect(),
                1 => l.sum_erosion.clone(),
                2 => l.sum_talus.clone(),
                _ => l.sum_deposit.clone(),
            };
            (l.substep_field(cu.sub), Some(d))
        }
        Sel::Physics => {
            let (_, r, _) = cu.physics.as_ref()?;
            let k = n2m();
            let d = match cu.sub {
                0 => r.z.data.iter().zip(&r.upscaled.data).map(|(a, b)| (a - b) * k).collect(),
                1 => r.sum_uplift.iter().map(|v| v * k).collect(),
                2 => r.sum_erosion.iter().map(|v| v * k).collect(),
                _ => r.sum_diffusion.iter().map(|v| v * k).collect(),
            };
            (to_m(&r.substep_field(SubStep::ALL[cu.sub.min(3)])), Some(d))
        }
        Sel::Retargeted => {
            let last = cu.levels.last()?.0.result.clone();
            let d = s.result.data.iter().zip(&last.data).map(|(a, b)| a - b).collect();
            ((*s.result).clone(), Some(d))
        }
        Sel::Coarse => ((*s.result).clone(), None),
    };
    let mut rgba = match (cu.view, &diff) {
        (View::Diff, Some(d)) => {
            let mut a: Vec<f32> = d.iter().map(|v| v.abs()).collect();
            a.sort_by(|x, y| x.partial_cmp(y).unwrap());
            let sat = a[((a.len() - 1) as f32 * 0.98) as usize].max(1.0);
            cu.diff_sat = Some(sat);
            diff_rgba(d, n, n, sat)
        }
        (View::Hypso, _) => shade_m_rgba(&field, f32::MAX, depth),
        _ => shade_m_rgba(&field, cell_m, depth),
    };
    draw(&mut rgba, &s.hydro, cu.min_order, cu.rivers, cu.lakes);
    Some(egui::ColorImage::from_rgba_unmultiplied([n, n], &rgba))
}

/// F163-D: the physics level's variants.
#[derive(Clone, Copy, Debug, PartialEq)]
struct PhysOpts {
    /// The weighted random receiver's τ (`None` = D8).
    tau: Option<f32>,
    /// The initial roughness (m RMS, 0 = none).
    rough_m: f32,
    /// The incision area's MFD exponent (`None` = D8 area).
    mfd: Option<f32>,
}

impl Default for PhysOpts {
    /// P0, F159–F162's level.
    fn default() -> Self {
        Self { tau: None, rough_m: 0.0, mfd: Some(2.0) }
    }
}

fn phys_opts(cu: &CascadeUi) -> PhysOpts {
    PhysOpts { tau: cu.rr_on.then_some(cu.rr_tau), rough_m: if cu.rough_on { cu.rough_m } else { 0.0 }, mfd: cu.mfd_p }
}

/// F162-V: the remedies the worker builds its `AmpConfig` from.
#[derive(Clone, Debug, Default)]
struct Remedies {
    rho: Option<f32>,
    pi_m: Vec<(usize, f32)>,
    r4: bool,
    /// F164-V: the warp of the upscaled field.
    warp: Option<Warp>,
}

/// The remedies as ticked; π's amplitudes are the factor × Corsica's octave 0 at each level (512² … 8 192²).
fn remedies(cu: &CascadeUi) -> Remedies {
    Remedies {
        rho: cu.rho_on.then_some(cu.rho_amp),
        pi_m: if cu.pi_on {
            [512usize, 1024, 2048, 4096, 8192].iter().filter_map(|&n| corse_target(n).map(|t| (n, cu.pi_factor * t))).collect()
        } else {
            Vec::new()
        },
        r4: cu.r4_on,
        warp: cu.warp_place.map(|place| Warp { amp: cu.warp_amp, corr: cu.warp_corr, place }),
    }
}

/// F161-V: the Corsica grid at a cascade level of `n` cells (its 200 km grid of n / 2 cells, the same cell), from
/// `YMIR_CORSICA_DIR` or `data/corsica` (`prep_corse.py`; Copernicus DEM GLO-30, licence in `docs/refs`).
fn corse_grid(n: usize) -> Option<GridF32> {
    let dir = std::env::var("YMIR_CORSICA_DIR").unwrap_or_else(|_| "data/corsica".into());
    let b = std::fs::read(format!("{dir}/corse_{}.bin", n / 2)).ok()?;
    let w = u32::from_le_bytes(b.get(0..4)?.try_into().ok()?) as usize;
    let h = u32::from_le_bytes(b.get(4..8)?.try_into().ok()?) as usize;
    let data: Vec<f32> = b[8..].chunks_exact(4).map(|c| f32::from_le_bytes(c.try_into().unwrap())).collect();
    (data.len() == w * h).then_some(GridF32 { width: w, height: h, data })
}

/// F161-D3: Corsica's 2–4 cell octave at a cascade level of `n` cells (the calibration's target), `None` without data.
fn corse_target(n: usize) -> Option<f32> {
    let g = corse_grid(n)?;
    let ocean = ymir_core::cascade::amplify::ocean_mask(&g);
    let zm: Vec<f32> = g.data.iter().zip(&ocean).map(|(v, o)| if *o { 0.0 } else { *v }).collect();
    let land = ocean.iter().filter(|o| !**o).count();
    ymir_core::cascade::measure::octave_rms(&zm, g.width, land).first().copied()
}

/// F161-V: the « Référence Corse » texture for a cascade level of `n` cells (built once per level), the same shading,
/// rivers (order ≥ 2) and lakes.
fn corse_texture(ctx: &egui::Context, cache: &mut Option<(usize, egui::TextureHandle)>, n: usize) -> Option<egui::TextureHandle> {
    if cache.as_ref().map(|c| c.0) != Some(n) {
        let g = corse_grid(n)?;
        let m = g.width;
        let hy = hydrology(&g, 200.0 / m as f32);
        let depth = SteinSteinParams::default().depth_scale_m as f32;
        let mut rgba = shade_m_rgba(&g, 200_000.0 / m as f32, depth);
        draw(&mut rgba, &hy, 2, true, true);
        let tex = ctx.load_texture(format!("corse{n}"), egui::ColorImage::from_rgba_unmultiplied([m, m], &rgba), egui::TextureOptions::NEAREST);
        *cache = Some((n, tex));
    }
    cache.as_ref().map(|c| c.1.clone())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// ADR Findings 159–160 -- the window's worker: the coarse is the framed C1 altitude (the roll is a pure
    /// relabelling); « Calculer la physique » gives a physics level calibrated on the peak; « Niveau suivant » one
    /// amplification level at the given budget; « + érosion » adds iterations on top of the level's state; « Recalage
    /// final » retargets. Tiny budgets keep it a unit test.
    #[test]
    fn the_cascade_worker_runs_physics_amplification_and_retargeting() {
        let w = spawn(1, [6, 37], CascadeConfig::declared(5000.0));
        let wait = |want: fn(&Evt) -> bool| -> Vec<Evt> {
            let mut got = Vec::new();
            loop {
                let e = w.rx.recv_timeout(std::time::Duration::from_secs(300)).expect("the worker answers");
                let stop = want(&e);
                got.push(e);
                if stop {
                    return got;
                }
            }
        };
        let first = wait(|e| matches!(e, Evt::Coarse(..)));
        let Some(Evt::Coarse(s, _)) = first.last() else { unreachable!() };
        assert_eq!(s.n, 64);
        let unrolled = cascade_coarse(1, [0, 0]);
        let framed = cascade_coarse(1, [6, 37]);
        assert_eq!(framed.data[0], unrolled.data[37 * 64 + 6], "the roll is the framing's");
        w.tx.send(Cmd::Physics(128, PhysOpts::default())).unwrap();
        let e = wait(|e| matches!(e, Evt::Idle));
        let phys = e.iter().find_map(|e| if let Evt::Physics(s, _, _) = e { Some(s.clone()) } else { None }).expect("a physics level");
        assert_eq!(phys.n, 128);
        let peak = phys.result.data.iter().fold(f32::MIN, |a, &v| a.max(v));
        assert!((peak - PEAK_TARGET_M).abs() < 0.2 * PEAK_TARGET_M, "the calibrated peak {peak}");
        w.tx.send(Cmd::Next(Budget { erosion: 3, talus: 3, deposit: 3 }, None)).unwrap();
        let e = wait(|e| matches!(e, Evt::Idle));
        let l = e.iter().find_map(|e| if let Evt::Level(0, _, l) = e { Some(l.clone()) } else { None }).expect("one level");
        assert_eq!((l.n, l.done.erosion, l.done.talus, l.done.deposit), (256, 3, 3, 3));
        w.tx.send(Cmd::Add(0, Process::Erosion, 2)).unwrap();
        let e = wait(|e| matches!(e, Evt::Idle));
        let l2 = e.iter().find_map(|e| if let Evt::Level(0, _, l) = e { Some(l.clone()) } else { None }).expect("the level updated");
        assert_eq!(l2.done.erosion, 5, "« + érosion » adds on top of the existing iterations");
        w.tx.send(Cmd::Retarget).unwrap();
        let e = wait(|e| matches!(e, Evt::Idle));
        assert!(e.iter().any(|e| matches!(e, Evt::Retargeted(..))));
    }
}
