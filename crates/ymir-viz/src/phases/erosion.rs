//! Erosion phase — particle hydraulic erosion (Beyer 2015) on the
//! UpscaleFbm output.
//!
//! Reads [`super::upscale_fbm::FbmCache`], dispatches
//! `ymir_core::erosion::hydraulic::run_erosion` on a background
//! thread, and stores the eroded heightmap + sediment field when the
//! worker finishes. Threading is required because synchronous
//! compute at the default droplet counts (200k–500k) blocks the
//! Bevy main thread for 1–10 s — long enough that Windows marks
//! the app "Not Responding" and the user thinks the click did
//! nothing. The legacy bridge used the same worker-thread pattern.

use bevy::prelude::*;
use bevy_egui::egui;
use crossbeam_channel::{bounded, Receiver};
use std::time::Instant;
use ymir_core::erosion::hydraulic::{run_erosion, ErosionConfig, ErosionResult};
use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;

use super::isostasy::IsostasyCache;
use super::upscale_fbm::FbmCache;
use crate::bridge::v2::V2SolverBridge;
use crate::pipeline::{ActivePhase, PipelinePhase};
use crate::visualization::v2_viz::{V2VizSprite, V2VizState};

#[derive(Resource, Clone, Debug)]
pub struct ErosionParams {
    pub num_droplets: usize,
    pub deposition_rate: f32,
    pub erosion_rate: f32,
    pub inertia: f32,
    pub evaporation_rate: f32,
    pub max_lifetime: usize,
    pub erosion_radius: usize,
    pub coastal_deposition_range: usize,
    pub seed: u64,
    pub recompute_requested: bool,
    /// Step 8.6 follow-up — set by the panel's "↻ Continue" button.
    /// Consumed by `handle_erosion_compute` like `recompute_requested`,
    /// but starts the worker from the existing eroded heightmap rather
    /// than the FBM input — letting the user pile additional droplets
    /// on top of a prior run.
    pub continue_requested: bool,
}

impl Default for ErosionParams {
    fn default() -> Self {
        let cfg = ErosionConfig::default();
        Self {
            // 1 M droplets ≈ 1 per cell at 1024², enough to carve
            // valleys + deposit visible sediment. The core default
            // is 5 M (production), the previous viz default was
            // 200 k (under-eroded — barely visible). Threaded
            // dispatch makes 1 M painless on the UI thread (≈ 5–8 s
            // wall clock, no freeze).
            num_droplets: 1_000_000,
            deposition_rate: cfg.deposition_rate,
            erosion_rate: cfg.erosion_rate,
            inertia: cfg.inertia,
            evaporation_rate: cfg.evaporation_rate,
            max_lifetime: cfg.max_lifetime,
            erosion_radius: cfg.erosion_radius,
            coastal_deposition_range: cfg.coastal_deposition_range,
            seed: 42,
            recompute_requested: false,
            continue_requested: false,
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ErosionState {
    Idle,
    Running { dispatched_at: Instant },
    Completed,
}

/// Messages the erosion worker thread streams back to the main
/// thread. `Snapshot` fires once per `run_erosion` batch (default
/// every 50 k droplets, ~5–10× per total run at the 1 M default)
/// so the render system can paint the in-progress heightmap and
/// the panel can show a live progress bar — same UX shape as the
/// per-step v2 progress events. `Done` fires once at the end with
/// the full `ErosionResult` (heightmap + sediment + stats).
pub enum ErosionMessage {
    Snapshot {
        heightmap: GridF32,
        completed: usize,
        total: usize,
    },
    Done(ErosionResult),
}

#[derive(Resource)]
pub struct ErosionCache {
    pub state: ErosionState,
    pub result: Option<ErosionResult>,
    /// In-progress heightmap snapshot (latest worker batch). The
    /// render system prefers this over `result.heightmap` when
    /// `state == Running` so the user sees the erosion accumulate
    /// in real time. Cleared on `Done`.
    pub preview_heightmap: Option<GridF32>,
    /// `(droplets_completed, total_droplets)` from the latest
    /// snapshot. Drives the panel's progress bar.
    pub progress: Option<(usize, usize)>,
    pub last_status: Option<Result<String, String>>,
    last_signature: Option<u64>,
    /// Receiver from the worker thread. Drained by
    /// `poll_erosion_result` each frame.
    receiver: Option<Receiver<ErosionMessage>>,
}

impl Default for ErosionCache {
    fn default() -> Self {
        Self {
            state: ErosionState::Idle,
            result: None,
            preview_heightmap: None,
            progress: None,
            last_status: None,
            last_signature: None,
            receiver: None,
        }
    }
}

impl ErosionCache {
    pub fn mark_dirty(&mut self) {
        self.last_signature = None;
    }
}

pub fn draw_section(
    ui: &mut egui::Ui,
    params: &mut ErosionParams,
    cache: &ErosionCache,
    fbm_ready: bool,
) {
    let is_running = matches!(cache.state, ErosionState::Running { .. });
    let can_run = fbm_ready && !is_running;
    ui.add_space(4.0);
    ui.heading("Erosion (hydraulic)");
    ui.label(
        egui::RichText::new(
            "Particle hydraulic erosion (Beyer 2015) on the upscaled \
             heightmap. Synchronous — UI freezes briefly at high \
             droplet counts.",
        )
        .small()
        .weak(),
    );
    ui.add_space(4.0);

    ui.add(
        egui::Slider::new(&mut params.num_droplets, 50_000..=10_000_000)
            .text("droplets")
            .logarithmic(true),
    );
    ui.add(
        egui::Slider::new(&mut params.erosion_rate, 0.0..=1.0)
            .text("erosion rate")
            .step_by(0.01),
    );
    ui.add(
        egui::Slider::new(&mut params.deposition_rate, 0.0..=1.0)
            .text("deposition rate")
            .step_by(0.01),
    );
    ui.add(
        egui::Slider::new(&mut params.inertia, 0.0..=0.5)
            .text("inertia")
            .step_by(0.01),
    );
    ui.add(
        egui::Slider::new(&mut params.evaporation_rate, 0.0..=0.05)
            .text("evaporation")
            .step_by(0.001),
    );
    ui.add(
        egui::Slider::new(&mut params.max_lifetime, 30..=400)
            .text("max lifetime"),
    );
    ui.add(
        egui::Slider::new(&mut params.erosion_radius, 1..=8)
            .text("brush radius (cells)"),
    );
    ui.add(
        egui::Slider::new(&mut params.coastal_deposition_range, 0..=30)
            .text("coastal deposit range"),
    );
    ui.add(egui::DragValue::new(&mut params.seed).prefix("seed = "));

    ui.add_space(4.0);
    let can_continue = cache.result.is_some() && !is_running;
    ui.horizontal(|ui| {
        if ui
            .add_enabled(can_run, egui::Button::new("⚙ Run Erosion"))
            .on_hover_text(
                "Run hydraulic erosion on the upscaled heightmap. \
                 Threaded — UI stays responsive while the worker \
                 runs.",
            )
            .clicked()
        {
            params.recompute_requested = true;
        }
        if ui
            .add_enabled(can_continue, egui::Button::new("\u{21bb} Continue"))
            .on_hover_text(
                "Run another batch of droplets on top of the previous \
                 eroded heightmap. Each click stacks more erosion onto \
                 the existing landscape.",
            )
            .clicked()
        {
            params.continue_requested = true;
        }
    });

    if let ErosionState::Running { dispatched_at } = cache.state {
        let elapsed = dispatched_at.elapsed().as_secs_f64();
        ui.label(
            egui::RichText::new(format!("⏳ Computing… ({:.1}s elapsed)", elapsed))
                .small()
                .color(egui::Color32::YELLOW),
        );
        if let Some((completed, total)) = cache.progress {
            let frac = if total > 0 {
                (completed as f32 / total as f32).clamp(0.0, 1.0)
            } else {
                0.0
            };
            ui.add(
                egui::ProgressBar::new(frac)
                    .text(format!("{} / {} droplets", completed, total))
                    .desired_width(f32::INFINITY),
            );
        } else {
            ui.label(
                egui::RichText::new("(awaiting first batch snapshot)")
                    .small()
                    .weak(),
            );
        }
    } else if let Some(result) = &cache.result {
        ui.label(
            egui::RichText::new(format!(
                "✓ {}×{} eroded · {:.2e} deposited · avg lifetime {:.1}",
                result.heightmap.width,
                result.heightmap.height,
                result.stats.total_deposited,
                result.stats.avg_lifetime,
            ))
            .small()
            .color(egui::Color32::LIGHT_GREEN),
        );
    } else if let Some(Err(err)) = &cache.last_status {
        ui.label(
            egui::RichText::new(format!("Error: {}", err))
                .small()
                .color(egui::Color32::LIGHT_RED),
        );
    } else if !fbm_ready {
        ui.label(
            egui::RichText::new(
                "No erosion data — Run UpscaleFbm first, then click \
                 ⚙ Run Erosion.",
            )
            .small()
            .weak(),
        );
    } else {
        ui.label(
            egui::RichText::new(
                "Ready — click ⚙ Run Erosion to dispatch the worker.",
            )
            .small()
            .weak(),
        );
    }
}

pub fn handle_erosion_compute(
    fbm_cache: Res<FbmCache>,
    iso_cache: Res<IsostasyCache>,
    mut params: ResMut<ErosionParams>,
    mut cache: ResMut<ErosionCache>,
) {
    let do_run = params.recompute_requested;
    let do_continue = params.continue_requested;
    if !do_run && !do_continue {
        return;
    }
    params.recompute_requested = false;
    params.continue_requested = false;

    if matches!(cache.state, ErosionState::Running { .. }) {
        // A worker is already in flight — ignore the new request
        // (the user double-clicked or held down the button).
        return;
    }

    // Continue starts from the prior eroded heightmap; Run resets to
    // the FBM input. If Continue is requested but no prior result
    // exists, fall back to a Run.
    let start_heightmap = if do_continue {
        if let Some(prev) = cache.result.as_ref() {
            prev.heightmap.clone()
        } else {
            cache.last_status = Some(Err(
                "Continue: no prior eroded heightmap. Click Run Erosion first.".to_string(),
            ));
            return;
        }
    } else {
        let Some(fbm) = fbm_cache.result.as_ref() else {
            cache.last_status = Some(Err(
                "No upscaled heightmap — run UpscaleFbm first.".to_string(),
            ));
            return;
        };
        fbm.heightmap.clone()
    };

    let cfg = ErosionConfig {
        num_droplets: params.num_droplets,
        deposition_rate: params.deposition_rate,
        erosion_rate: params.erosion_rate,
        inertia: params.inertia,
        evaporation_rate: params.evaporation_rate,
        max_lifetime: params.max_lifetime,
        erosion_radius: params.erosion_radius,
        coastal_deposition_range: params.coastal_deposition_range,
        sea_level: iso_cache
            .result
            .as_ref()
            .map(|r| r.sea_level_normalized)
            .unwrap_or(0.4),
        ..ErosionConfig::default()
    };
    let world_seed = WorldSeed::new(params.seed);
    let heightmap = start_heightmap;

    // Channel buffered at 8 — the worker fires one snapshot per
    // batch (~50 k droplets) and one final Done. At 1 M droplets
    // that's ≤ 20 messages over ~5–8 s; the main thread polls at
    // 60 fps so the channel never approaches saturation.
    let (tx, rx) = bounded::<ErosionMessage>(8);
    cache.state = ErosionState::Running {
        dispatched_at: Instant::now(),
    };
    cache.preview_heightmap = None;
    cache.progress = None;
    cache.receiver = Some(rx);
    cache.last_status = Some(Ok(format!(
        "dispatched {} droplets ({})",
        params.num_droplets,
        if do_continue { "continue" } else { "fresh" },
    )));

    info!(
        "[erosion] dispatching worker — {} droplets @ {}×{} ({})",
        cfg.num_droplets,
        heightmap.width,
        heightmap.height,
        if do_continue { "continue" } else { "fresh" },
    );
    std::thread::Builder::new()
        .name("ymir-erosion-worker".into())
        .spawn(move || {
            let t0 = Instant::now();
            let result = run_erosion(&heightmap, &cfg, &world_seed, |completed, total, hmap| {
                // Stream the in-progress heightmap so the main
                // thread can paint partial erosion. Clone is
                // ≈ width·height·4 B (16 MB at 2048², 4 MB at
                // 1024²) — well under one frame at memcpy speed.
                let _ = tx.send(ErosionMessage::Snapshot {
                    heightmap: hmap.clone(),
                    completed,
                    total,
                });
                true
            });
            info!(
                "[erosion] worker finished in {:.2}s — {:.2e} deposited",
                t0.elapsed().as_secs_f64(),
                result.stats.total_deposited
            );
            let _ = tx.send(ErosionMessage::Done(result));
        })
        .expect("failed to spawn erosion worker thread");
}

/// Drain the worker channel — runs each frame. Each `Snapshot`
/// updates `preview_heightmap` + `progress` so the render system
/// paints the in-progress heightmap and the panel shows a live
/// progress bar; `Done` finalises by storing the full result and
/// transitioning to `Completed`.
pub fn poll_erosion_result(mut cache: ResMut<ErosionCache>) {
    let still_running = matches!(cache.state, ErosionState::Running { .. });
    if !still_running {
        return;
    }
    let Some(rx) = cache.receiver.as_ref() else {
        cache.state = ErosionState::Idle;
        return;
    };
    // Drain everything queued this frame so we land on the latest
    // snapshot rather than stale ones if the worker outpaces the
    // poll.
    let mut last_snapshot: Option<(GridF32, usize, usize)> = None;
    let mut done: Option<ErosionResult> = None;
    while let Ok(msg) = rx.try_recv() {
        match msg {
            ErosionMessage::Snapshot { heightmap, completed, total } => {
                last_snapshot = Some((heightmap, completed, total));
            }
            ErosionMessage::Done(result) => {
                done = Some(result);
            }
        }
    }
    if let Some((heightmap, completed, total)) = last_snapshot {
        cache.preview_heightmap = Some(heightmap);
        cache.progress = Some((completed, total));
        cache.last_signature = None;
    }
    if let Some(result) = done {
        cache.last_status = Some(Ok(format!(
            "Eroded {}×{}, {:.2e} deposited",
            result.heightmap.width, result.heightmap.height, result.stats.total_deposited
        )));
        cache.result = Some(result);
        cache.preview_heightmap = None;
        cache.progress = None;
        cache.state = ErosionState::Completed;
        cache.receiver = None;
        cache.last_signature = None;
    }
}

#[allow(clippy::too_many_arguments)]
pub fn render_erosion_phase(
    active: Res<ActivePhase>,
    mut cache: ResMut<ErosionCache>,
    iso_cache: Res<IsostasyCache>,
    fbm_cache: Res<FbmCache>,
    bridge: Res<V2SolverBridge>,
    viz: Res<V2VizState>,
    mut images: ResMut<Assets<Image>>,
    mut sprite_q: Query<&mut Sprite, With<V2VizSprite>>,
) {
    if active.0 != PipelinePhase::Erosion {
        return;
    }
    let _ = iso_cache; // sea_level is baked into hypsometric_colormap's [0,1] mapping
    // Prefer the in-progress snapshot, then the cached result, then
    // walk upstream so the user sees the FBM / isostasy / V2 state
    // before erosion has run at all.
    let Some((grid, source_tag)) = super::select_grid_for_phase(
        PipelinePhase::Erosion,
        &bridge,
        &iso_cache,
        &fbm_cache,
        &cache,
    ) else {
        return;
    };
    let nx = grid.width;
    let ny = grid.height;
    // Signature includes the source tag and a `progress_token` so
    // each in-progress snapshot still triggers a re-render and a
    // fallback → primary upgrade does too.
    let progress_token = cache.progress.map(|(c, _)| c as u64).unwrap_or(0);
    let signature =
        ((nx as u64) << 48) | ((ny as u64) << 32) | (progress_token & 0x00FF_FFFF) | ((source_tag as u64) << 24);
    if cache.last_signature == Some(signature) {
        return;
    }
    super::paint_grid_to_v2_sprite(&grid, &viz, &mut images, &mut sprite_q);
    cache.last_signature = Some(signature);
}
