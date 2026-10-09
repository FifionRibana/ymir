//! ADR Finding 159 -- the « Cascade » window: the multi-scale cascade prototype (`ymir_core::cascade`), level by level.
//!
//! The author (2026-10-09): « on aurait clairement les différentes résolutions qui seront balayées, soit d'un seul
//! trait, soit résolution par résolution, par action utilisateur … visualiser les résultats intermédiaires. Ainsi, si
//! quelque chose ne va pas, on ne va pas plus loin. »
//!
//! A floating window, opened from the top bar, on the world the workspace frames (its seed and its framing roll). It
//! never reads or writes the HD world: its own worker thread runs the 64² tectonics and the cascade, and the window only
//! displays. « Niveau suivant » computes one level, « Tout » the remaining ones; each level and each of its sub-steps
//! (agrandissement, soulèvement, érosion, diffusion) can be viewed once computed, as Ombrage, hypsometry, or the
//! difference view (what the level, or the selected sub-step's process, changed), with the closed depressions' outline.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

use bevy::prelude::*;
use bevy_egui::{EguiContexts, EguiPrimaryContextPass, egui};
use crossbeam_channel::{Receiver, Sender, unbounded};
use egui::Color32 as C;
use ymir_core::cascade::{Calibration, Cascade, CascadeConfig, CascadeProgress, LevelRecord, SubStep, diff_rgba, roll, shade_rgba};
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

pub struct CascadePlugin;

impl Plugin for CascadePlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<CascadeUi>();
        app.add_systems(EguiPrimaryContextPass, draw_cascade);
    }
}

enum Cmd {
    Next,
    All,
}

enum Evt {
    Coarse(Arc<GridF32>, f64),
    Progress(CascadeProgress),
    Calibrated(Calibration),
    Level(Arc<LevelRecord>),
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

/// The window's state. `seed` / `offset_cells` are the workspace's world, set by its top-bar button.
#[derive(Resource)]
pub struct CascadeUi {
    pub open: bool,
    pub seed: u64,
    pub offset_cells: [i64; 2],
    worker: Option<Worker>,
    /// The world the worker was started on.
    world: Option<(u64, [i64; 2])>,
    coarse: Option<(Arc<GridF32>, f64)>,
    calibration: Option<Calibration>,
    levels: Vec<Arc<LevelRecord>>,
    busy: bool,
    progress: Option<CascadeProgress>,
    /// 0 = the 64² tectonics, i = cascade level i − 1.
    sel_level: usize,
    sel_sub: SubStep,
    view: View,
    outline: bool,
    texture: Option<egui::TextureHandle>,
    tex_key: Option<(usize, SubStep, View, bool, usize)>,
    /// The field and cell size of the displayed texture, for the hover readout.
    shown: Option<(GridF32, f32)>,
    /// The difference view's saturation (m), the p98 of |Δ|.
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
            calibration: None,
            levels: Vec::new(),
            busy: false,
            progress: None,
            sel_level: 0,
            sel_sub: SubStep::Diffusion,
            view: View::Shade,
            outline: true,
            texture: None,
            tex_key: None,
            shown: None,
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

fn spawn(seed: u64, offset_cells: [i64; 2], cfg: CascadeConfig) -> Worker {
    let (tx, cmd_rx) = unbounded::<Cmd>();
    let (evt_tx, rx) = unbounded::<Evt>();
    let cancel = Arc::new(AtomicBool::new(false));
    let stop = cancel.clone();
    std::thread::Builder::new()
        .name("ymir-cascade".into())
        .spawn(move || {
            let t = Instant::now();
            let coarse = cascade_coarse(seed, offset_cells);
            let _ = evt_tx.send(Evt::Coarse(Arc::new(coarse.clone()), t.elapsed().as_secs_f64()));
            let mut cas = Cascade::new(coarse, cfg);
            for cmd in cmd_rx.iter() {
                stop.store(false, Ordering::Relaxed);
                let all = matches!(cmd, Cmd::All);
                while !cas.done() {
                    let had_cal = cas.calibration.is_some();
                    let ptx = evt_tx.clone();
                    let level = cas
                        .next_level(&mut |p| drop(ptx.send(Evt::Progress(p))), &|| stop.load(Ordering::Relaxed))
                        .map(|l| Arc::new(l.clone()));
                    if !had_cal && let Some(c) = &cas.calibration {
                        let _ = evt_tx.send(Evt::Calibrated(c.clone()));
                    }
                    match level {
                        Some(l) => drop(evt_tx.send(Evt::Level(l))),
                        None => break,
                    }
                    if !all {
                        break;
                    }
                }
                let _ = evt_tx.send(Evt::Idle);
            }
        })
        .expect("spawn the cascade worker");
    Worker { tx, rx, cancel }
}

impl CascadeUi {
    fn start(&mut self) {
        let ss = SteinSteinParams::default();
        self.worker = Some(spawn(self.seed, self.offset_cells, CascadeConfig::declared(ss.depth_scale_m as f32)));
        self.world = Some((self.seed, self.offset_cells));
        self.coarse = None;
        self.calibration = None;
        self.levels.clear();
        self.busy = true;
        self.progress = None;
        self.sel_level = 0;
        self.texture = None;
        self.tex_key = None;
    }

    fn drain(&mut self) {
        let Some(w) = &self.worker else { return };
        let evts: Vec<Evt> = w.rx.try_iter().collect();
        for e in evts {
            match e {
                Evt::Coarse(g, s) => {
                    self.coarse = Some((g, s));
                    self.busy = false;
                }
                Evt::Progress(p) => self.progress = Some(p),
                Evt::Calibrated(c) => self.calibration = Some(c),
                Evt::Level(l) => {
                    self.levels.push(l);
                    self.sel_level = self.levels.len();
                    self.sel_sub = SubStep::Diffusion;
                }
                Evt::Idle => {
                    self.busy = false;
                    self.progress = None;
                }
            }
        }
    }

    fn send(&mut self, c: Cmd) {
        if let Some(w) = &self.worker {
            if w.tx.send(c).is_ok() {
                self.busy = true;
            }
        }
    }

    /// The displayed field and its cell (m), for (level, sub-step); `None` when not computed.
    fn field(&self, level: usize, sub: SubStep) -> Option<(GridF32, f32)> {
        if level == 0 {
            return self.coarse.as_ref().map(|(g, _)| ((**g).clone(), 400_000.0 / g.width as f32));
        }
        let l = self.levels.get(level - 1)?;
        Some((l.substep_field(sub), l.cell_km * 1000.0))
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
    egui::Window::new("Cascade multi-échelle (F159)")
        .open(&mut open)
        .default_size([780.0, 720.0])
        .resizable(true)
        .show(ctx, |ui| body(ui, &mut cu));
    cu.open = open;
}

fn button(ui: &mut egui::Ui, text: &str, enabled: bool) -> bool {
    let b = egui::Button::new(egui::RichText::new(text).color(INK).strong().size(12.0)).fill(COPPER).corner_radius(5.0);
    ui.add_enabled(enabled, b).clicked()
}

fn body(ui: &mut egui::Ui, cu: &mut CascadeUi) {
    ui.horizontal(|ui| {
        ui.label(
            egui::RichText::new(format!(
                "Monde : graine {} · cadrage ({}, {}) · 400 km · prototype, la production n'est pas touchée",
                cu.seed, cu.offset_cells[0], cu.offset_cells[1]
            ))
            .color(DIM)
            .size(11.0),
        );
    });
    let stale = cu.world.is_some_and(|w| w != (cu.seed, cu.offset_cells));
    ui.horizontal(|ui| {
        if cu.worker.is_none() || stale {
            let label = if stale { "Recommencer sur le monde courant" } else { "Démarrer (tectonique 64²)" };
            if button(ui, label, !cu.busy || stale) {
                if let Some(w) = &cu.worker {
                    w.cancel.store(true, Ordering::Relaxed);
                }
                cu.start();
            }
        }
        let ready = cu.worker.is_some() && cu.coarse.is_some() && !cu.busy && !stale;
        let remaining = cu.levels.len() < CascadeConfig::declared(5000.0).levels.len();
        if button(ui, "Niveau suivant", ready && remaining) {
            cu.send(Cmd::Next);
        }
        if button(ui, "Tout", ready && remaining) {
            cu.send(Cmd::All);
        }
        if cu.busy && cu.coarse.is_some() && ui.button("Annuler").clicked() {
            if let Some(w) = &cu.worker {
                w.cancel.store(true, Ordering::Relaxed);
            }
        }
    });
    if let Some(p) = cu.progress {
        let what = if p.level.is_none() { "calibration de U₀" } else { "niveau" };
        ui.add(
            egui::ProgressBar::new(p.step as f32 / p.max_steps as f32)
                .text(format!("{what} {}² · pas {} / {} · |Δz| / U·dt = {:.4}", p.n, p.step, p.max_steps, p.ratio)),
        );
    } else if cu.busy && cu.coarse.is_none() {
        ui.label("tectonique 64² en cours…");
    }
    ui.separator();
    // ── the frieze of levels ──
    let cfg = CascadeConfig::declared(SteinSteinParams::default().depth_scale_m as f32);
    ui.horizontal(|ui| {
        let n_nodes = 1 + cfg.levels.len();
        for i in 0..n_nodes {
            let (title, status, col) = if i == 0 {
                match &cu.coarse {
                    Some((_, s)) => ("64² tectonique".to_string(), format!("✓ {s:.1} s"), GREEN),
                    None => ("64² tectonique".to_string(), "—".to_string(), DIM),
                }
            } else {
                let n = cfg.levels[i - 1];
                match cu.levels.get(i - 1) {
                    Some(l) => (
                        format!("{n}²"),
                        format!(
                            "{} {} pas · {:.1} s",
                            if l.equilibrium { "✓" } else { "⚠ plafond" },
                            l.steps,
                            l.total_secs()
                        ),
                        if l.equilibrium { GREEN } else { AMBER },
                    ),
                    None => {
                        let running = cu.progress.is_some_and(|p| p.n == n && p.level.is_some());
                        (format!("{n}²"), if running { "en cours…".into() } else { "—".into() }, DIM)
                    }
                }
            };
            let avail = i == 0 && cu.coarse.is_some() || i > 0 && cu.levels.len() >= i;
            let sel = cu.sel_level == i;
            let text = egui::RichText::new(format!("{title}\n{status}")).size(11.0).color(if sel { INK } else { col });
            let b = egui::Button::new(text).fill(if sel { COPPER } else { C::TRANSPARENT }).min_size(egui::vec2(140.0, 36.0));
            if ui.add_enabled(avail, b).clicked() {
                cu.sel_level = i;
            }
            if i + 1 < n_nodes {
                ui.label(egui::RichText::new("→").color(DIM));
            }
        }
    });
    // ── the sub-steps of the selected level ──
    if cu.sel_level > 0 {
        if let Some(l) = cu.levels.get(cu.sel_level - 1).cloned() {
            ui.horizontal(|ui| {
                for (j, s) in SubStep::ALL.iter().enumerate() {
                    let sel = cu.sel_sub == *s;
                    let text = egui::RichText::new(format!("{} · {:.2} s", s.label(), l.secs[j]))
                        .size(11.0)
                        .color(if sel { INK } else { DIM });
                    if ui.add(egui::Button::new(text).fill(if sel { COPPER } else { C::TRANSPARENT })).clicked() {
                        cu.sel_sub = *s;
                    }
                }
            });
            let mut info = format!(
                "cellule {:.3} km · D {:.3} m²/an · L_c {:.2} cellule · poids {:.4} · {} pas, dernier |Δz| / U·dt {:.4} · talus {:.2} % \
                 des terres au dernier pas",
                l.cell_km,
                l.d_m2_yr,
                (l.d_m2_yr / (cfg.k / 1000.0)).sqrt() / (l.cell_km * 1000.0),
                l.weight,
                l.steps,
                l.ratios.last().copied().unwrap_or(f32::NAN),
                100.0 * l.talus_last_share
            );
            if cu.sel_level == 1 {
                if let Some(c) = &cu.calibration {
                    info.push_str(&format!(
                        "\ncalibration : U₀ {:.2e} m/an par km de h_iso (essai {:.0} m contre {:.0} m, {} pas, {:.1} s)",
                        c.u0, c.trial_mean_m, c.target_mean_m, c.steps, c.secs
                    ));
                }
            }
            ui.label(egui::RichText::new(info).color(DIM).size(11.0));
        }
    }
    ui.horizontal(|ui| {
        for (v, name) in [(View::Shade, "Ombrage"), (View::Hypso, "Hypsométrie"), (View::Diff, "Différence (niveau − agrandissement)")] {
            if ui.selectable_label(cu.view == v, name).clicked() {
                cu.view = v;
            }
        }
        ui.checkbox(&mut cu.outline, "Contour des dépressions");
    });
    // ── the image ──
    let sub = if cu.sel_level == 0 { SubStep::Upscale } else { cu.sel_sub };
    let key = (cu.sel_level, sub, cu.view, cu.outline, cu.levels.len());
    if cu.tex_key != Some(key) {
        if let Some((g, cell_m)) = cu.field(cu.sel_level, sub) {
            // ADR Finding 159 (after run 1) -- the difference view shows what the selected sub-step DID: the level's
            // total change for « Agrandissement », then Σ uplift / Σ erosion / Σ diffusion (m), saturating at the p98
            let diff = if cu.view == View::Diff && cu.sel_level > 0 {
                cu.levels.get(cu.sel_level - 1).map(|l| {
                    let n2m = 2.0 * 1.13 * SteinSteinParams::default().depth_scale_m as f32;
                    let d: Vec<f32> = match sub {
                        SubStep::Upscale => l.z.data.iter().zip(&l.upscaled.data).map(|(a, b)| (a - b) * n2m).collect(),
                        SubStep::Uplift => l.sum_uplift.iter().map(|v| v * n2m).collect(),
                        SubStep::Erosion => l.sum_erosion.iter().map(|v| v * n2m).collect(),
                        SubStep::Diffusion => l.sum_diffusion.iter().map(|v| v * n2m).collect(),
                    };
                    let mut a: Vec<f32> = d.iter().map(|v| v.abs()).collect();
                    a.sort_by(|x, y| x.partial_cmp(y).unwrap());
                    let sat = a[((a.len() - 1) as f32 * 0.98) as usize].max(1.0);
                    (d, sat)
                })
            } else {
                None
            };
            cu.diff_sat = diff.as_ref().map(|d| d.1);
            let img = render(&g, cell_m, cu.view, diff.as_ref(), cu.outline);
            cu.texture = Some(ui.ctx().load_texture("cascade", img, egui::TextureOptions::NEAREST));
            cu.tex_key = Some(key);
            cu.shown = Some((g, cell_m));
        }
    }
    if let (Some(tex), Some((g, cell_m))) = (&cu.texture, &cu.shown) {
        let side = ui.available_width().min(ui.available_height() - 20.0).max(160.0);
        let resp = ui.add(egui::Image::new((tex.id(), egui::vec2(side, side))).sense(egui::Sense::hover()));
        if let Some(pos) = resp.hover_pos() {
            let r = resp.rect;
            let (u, v) = ((pos.x - r.left()) / r.width(), (pos.y - r.top()) / r.height());
            let (w, h) = (g.width, g.height);
            let x = ((u * w as f32) as usize).min(w - 1);
            // north up: the screen's top row is the field's last
            let y = h - 1 - ((v * h as f32) as usize).min(h - 1);
            let m = (g.data[y * w + x] - 0.5) * 2.0 * 1.13 * SteinSteinParams::default().depth_scale_m as f32;
            ui.label(
                egui::RichText::new(format!("cellule ({x}, {y}) · {:.1} km, {:.1} km · {m:.0} m", x as f32 * cell_m / 1000.0, y as f32 * cell_m / 1000.0))
                    .color(DIM)
                    .size(11.0)
                    .monospace(),
            );
        }
        if let (View::Diff, Some(sat)) = (cu.view, cu.diff_sat) {
            let what = match cu.sel_sub {
                SubStep::Upscale => "le changement total du niveau (niveau − agrandissement)",
                SubStep::Uplift => "Σ soulèvement",
                SubStep::Erosion => "Σ érosion",
                SubStep::Diffusion => "Σ diffusion (talus + diffusion)",
            };
            ui.label(egui::RichText::new(format!("{what} · bleu = abaissé, rouge = relevé, saturé à {sat:.0} m (p98)")).color(DIM).size(11.0));
        }
    }
}

/// The texture: Ombrage / hypsometry / difference against `base` (m), north up, with the closed depressions outlined
/// (cells the priority flood fills by more than 0.5 m; the cascade has no lake inventory, so this is its stand-in for
/// the workspace's « Contour des lacs »).
fn render(g: &GridF32, cell_m: f32, view: View, diff: Option<&(Vec<f32>, f32)>, outline: bool) -> egui::ColorImage {
    let (w, h) = (g.width, g.height);
    let depth = SteinSteinParams::default().depth_scale_m as f32;
    let n2m = 2.0 * 1.13 * depth;
    let mut rgba = match (view, diff) {
        (View::Diff, Some((d, sat))) => diff_rgba(d, w, h, *sat),
        (View::Hypso, _) => {
            // the shade at a flat sun = the bare palette: a constant field has hillshade 1 → factor 1
            let mut flat = shade_rgba(g, f32::MAX, depth);
            for p in flat.chunks_exact_mut(4) {
                p[3] = 255;
            }
            flat
        }
        _ => shade_rgba(g, cell_m, depth),
    };
    if outline {
        let flow = ymir_core::terrain::flow::compute_flow(g, &ymir_core::terrain::flow::FlowConfig { sea_level: 0.5, ..Default::default() });
        let pit: Vec<bool> = flow.filled.data.iter().zip(&g.data).map(|(f, z)| (f - z) * n2m > 0.5).collect();
        for y in 0..h {
            for x in 0..w {
                let k = y * w + x;
                if !pit[k] {
                    continue;
                }
                let edge = (x == 0 || !pit[k - 1]) || (x + 1 == w || !pit[k + 1]) || (y == 0 || !pit[k - w]) || (y + 1 == h || !pit[k + w]);
                if edge {
                    let o = ((h - 1 - y) * w + x) * 4;
                    rgba[o..o + 4].copy_from_slice(&[0x4F, 0xC3, 0xF7, 255]);
                }
            }
        }
    }
    egui::ColorImage::from_rgba_unmultiplied([w, h], &rgba)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// ADR Finding 159 -- the window's worker answers « Niveau suivant » with ONE level and « Tout » with the rest, then
    /// idles; the coarse field it starts from is the framed C1 altitude (the roll is a pure relabelling). A short config
    /// (two levels, a few steps, no calibration) keeps it a unit test.
    #[test]
    fn the_cascade_worker_runs_level_by_level_then_all() {
        let mut cfg = CascadeConfig::declared(5000.0);
        cfg.levels = vec![128, 256];
        cfg.max_steps = vec![3, 2];
        cfg.calibrate = false;
        let w = spawn(1, [6, 37], cfg);
        let wait = |want: fn(&Evt) -> bool| -> Vec<Evt> {
            let mut got = Vec::new();
            loop {
                let e = w.rx.recv_timeout(std::time::Duration::from_secs(120)).expect("the worker answers");
                let stop = want(&e);
                got.push(e);
                if stop {
                    return got;
                }
            }
        };
        let first = wait(|e| matches!(e, Evt::Coarse(..)));
        let Some(Evt::Coarse(g, _)) = first.last() else { unreachable!() };
        assert_eq!(g.data, cascade_coarse(1, [6, 37]).data);
        let unrolled = cascade_coarse(1, [0, 0]);
        assert_eq!(g.data[0], unrolled.data[37 * 64 + 6], "the roll is the framing's");
        w.tx.send(Cmd::Next).unwrap();
        let e = wait(|e| matches!(e, Evt::Idle));
        assert_eq!(e.iter().filter(|e| matches!(e, Evt::Level(_))).count(), 1, "« Niveau suivant » computes one level");
        w.tx.send(Cmd::All).unwrap();
        let e = wait(|e| matches!(e, Evt::Idle));
        let lv: Vec<usize> = e.iter().filter_map(|e| if let Evt::Level(l) = e { Some(l.n) } else { None }).collect();
        assert_eq!(lv, vec![256], "« Tout » computes the rest");
    }
}
