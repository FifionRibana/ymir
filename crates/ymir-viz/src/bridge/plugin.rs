//! Bevy plugin that manages the solver bridge and polls events each frame.

use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::time::Duration;

use bevy::prelude::*;
use crossbeam_channel::{Receiver, Sender, bounded};

use ymir_core::tectonics::solver::config::{
    ContinuationConfig, NewtonConfig, PicardConfig, TectonicsConfig,
};
use ymir_core::tectonics::solver::tectonics::DynamicPlateContext;
use ymir_core::tectonics::solver::workspace::StepStats;

use super::commands::SolverCommand;
use super::events::SolverEvent;
use super::thread::spawn_solver_thread;
use crate::state::{
    DynamicPlateIds, ErosionCache, ErosionState, FbmState, FlowCache, FlowState, SolverConfig,
    TectonicState, UiActions, UpscaleCache,
};
use crate::visualization::render::TerrainDisplay;

/// Current state of the solver.
#[derive(Default, Debug, Clone)]
pub enum SolverState {
    #[default]
    Idle,
    Running {
        step: usize,
        total_steps: usize,
        stats: Option<StepStats>,
    },
    Completed {
        elapsed: Duration,
    },
    Failed {
        error: String,
    },
}

/// Resource holding the channels to communicate with the solver thread.
#[derive(Resource)]
pub struct SolverBridge {
    pub commands_tx: Sender<SolverCommand>,
    pub events_rx: Receiver<SolverEvent>,
    pub state: SolverState,
    pub cancel_flag: Arc<AtomicBool>,
}

pub struct TectonicsBridgePlugin;

impl Plugin for TectonicsBridgePlugin {
    fn build(&self, app: &mut App) {
        let (cmd_tx, cmd_rx) = bounded::<SolverCommand>(4);
        let (evt_tx, evt_rx) = bounded::<SolverEvent>(256);
        let cancel = Arc::new(AtomicBool::new(false));

        spawn_solver_thread(cmd_rx, evt_tx, cancel.clone());

        app.insert_resource(SolverBridge {
            commands_tx: cmd_tx,
            events_rx: evt_rx,
            state: SolverState::Idle,
            cancel_flag: cancel,
        });

        app.init_resource::<crate::state::UiActions>();
        app.init_resource::<DynamicPlateIds>();
        app.add_systems(
            Update,
            (
                poll_solver_events,
                handle_step,
                dispatch_fbm_upscale,
                dispatch_erosion,
                dispatch_flow,
                update_river_extraction,
                super::export_system::handle_export,
                super::export_system::handle_load,
            ),
        );
    }
}

fn poll_solver_events(
    mut bridge: ResMut<SolverBridge>,
    mut terrain_display: ResMut<TerrainDisplay>,
    mut isostasy_cache: ResMut<crate::state::IsostasyCache>,
    mut dynamic_plates: ResMut<DynamicPlateIds>,
    mut upscale_cache: ResMut<UpscaleCache>,
    mut erosion_cache: ResMut<ErosionCache>,
    mut flow_cache: ResMut<FlowCache>,
) {
    while let Ok(event) = bridge.events_rx.try_recv() {
        match event {
            SolverEvent::Progress { step, total_steps, stats } => {
                bridge.state = SolverState::Running { step, total_steps, stats: Some(stats) };
            }
            SolverEvent::Snapshot { s_field, plate_ids, plates, .. } => {
                let grid_size = s_field.n();
                terrain_display.update_field(s_field);
                isostasy_cache.valid = false;

                if let Some(ids) = plate_ids {
                    dynamic_plates.grid_size = grid_size;
                    dynamic_plates.ids = Some(ids);
                }
                if let Some(pl) = plates {
                    dynamic_plates.active_count = pl.iter().filter(|p| p.active).count();
                    dynamic_plates.plates = Some(pl);
                }
            }
            SolverEvent::Completed { s_field, plate_ids, plates, elapsed, .. } => {
                let grid_size = s_field.n();
                terrain_display.update_field(s_field);
                isostasy_cache.valid = false;
                bridge.state = SolverState::Completed { elapsed };

                if let Some(ids) = plate_ids {
                    dynamic_plates.grid_size = grid_size;
                    dynamic_plates.ids = Some(ids);
                }
                if let Some(pl) = plates {
                    dynamic_plates.active_count = pl.iter().filter(|p| p.active).count();
                    dynamic_plates.plates = Some(pl);
                }
            }
            SolverEvent::FbmCompleted { heightmap, slope, elapsed } => {
                upscale_cache.heightmap = Some(heightmap);
                upscale_cache.slope = Some(slope);
                upscale_cache.state = FbmState::Completed { elapsed };
            }
            SolverEvent::ErosionProgress { completed, total } => {
                erosion_cache.state = ErosionState::Running { completed, total };
            }
            SolverEvent::ErosionSnapshot { heightmap, completed, total } => {
                erosion_cache.heightmap = Some(heightmap);
                erosion_cache.state = ErosionState::Running { completed, total };
            }
            SolverEvent::ErosionCompleted { heightmap, sediment, stats, elapsed } => {
                erosion_cache.heightmap = Some(heightmap);
                erosion_cache.sediment = Some(sediment);
                erosion_cache.stats = Some(stats);
                erosion_cache.state = ErosionState::Completed { elapsed };
            }
            SolverEvent::FlowCompleted { result, elapsed } => {
                flow_cache.result = Some(result);
                flow_cache.state = FlowState::Completed { elapsed };
                flow_cache.rivers_dirty = true;
            }
            SolverEvent::Failed { error } => {
                bridge.state = SolverState::Failed { error };
            }
        }
    }
}

fn handle_step(
    mut ui_actions: ResMut<UiActions>,
    mut bridge: ResMut<SolverBridge>,
    terrain_display: Res<TerrainDisplay>,
    dynamic_plates: Res<DynamicPlateIds>,
    tectonic_state: Option<Res<TectonicState>>,
    solver_config: Res<SolverConfig>,
) {
    if !ui_actions.step_requested {
        return;
    }
    ui_actions.step_requested = false;

    // Don't send if solver is already running
    if matches!(bridge.state, SolverState::Running { .. }) {
        return;
    }

    let Some(ref tecto) = tectonic_state else {
        return;
    };

    let grid_size = tecto.init.grid_size;
    let dx = 1.0 / grid_size as f64;

    // Get current s_field from terrain display (last snapshot), or from initial thickness
    let s_field = if let Some(ref field) = terrain_display.s_field {
        field.clone()
    } else {
        let mut field = ymir_core::tectonics::solver::field::Field2D::new(grid_size);
        for j in 0..grid_size {
            for i in 0..grid_size {
                field.set(i, j, tecto.init.thickness.data[j * grid_size + i] as f64);
            }
        }
        field
    };

    // Get current plate context from dynamic state, or from initial
    let plate_ctx = if let (Some(ids), Some(plates)) = (&dynamic_plates.ids, &dynamic_plates.plates)
    {
        let traction = ymir_core::tectonics::plates::rebuild_traction(ids, plates, grid_size);
        DynamicPlateContext { ids: ids.clone(), plates: plates.clone(), traction }
    } else {
        let traction = tecto.init.to_traction_field();
        DynamicPlateContext {
            ids: tecto.init.plate_ids.clone(),
            plates: tecto.init.plates.clone(),
            traction,
        }
    };

    let config = build_tectonics_config(&solver_config);

    let _ = bridge.commands_tx.send(SolverCommand::SingleStep {
        config,
        plate_ctx,
        s_field,
        grid_size,
        dx,
    });

    bridge.state = SolverState::Running { step: 0, total_steps: 1, stats: None };
}

fn build_tectonics_config(sc: &SolverConfig) -> TectonicsConfig {
    TectonicsConfig {
        num_timesteps: sc.num_timesteps,
        gravity_factor: sc.gravity_factor,
        cfl_factor: sc.cfl_factor,
        s_min: 0.1,
        s_max: 2.5,
        nonlinear_solver: sc.nonlinear_solver,
        picard: PicardConfig {
            power_law_n: sc.power_law_n,
            relaxation: sc.picard_relaxation,
            strain_rate_min: sc.strain_rate_min,
            eta_max: sc.eta_max,
            ..PicardConfig::default()
        },
        newton: NewtonConfig {
            preconditioner: sc.preconditioner,
            inexact: sc.inexact_newton,
            ..NewtonConfig::default()
        },
        continuation: ContinuationConfig {
            enabled: sc.continuation_enabled,
            ..ContinuationConfig::default()
        },
        boundaries: sc.boundaries.clone(),
        dynamic_boundaries: sc.dynamic_boundaries,
        cratonic: sc.cratonic.clone(),
        yielding: sc.yielding.clone(),
    }
}

/// Dispatch a pending FBM upscale command to the solver thread.
///
/// The UI sets `upscale_cache.pending_config` when the user clicks "Run FBM".
/// This system picks it up, grabs the isostasy heightmap from the cache,
/// and sends the command to the background thread.
fn dispatch_fbm_upscale(
    mut upscale_cache: ResMut<UpscaleCache>,
    bridge: ResMut<SolverBridge>,
    isostasy_cache: Res<crate::state::IsostasyCache>,
) {
    let (Some(config), Some(seed), Some(sea_level)) = (
        upscale_cache.pending_config.take(),
        upscale_cache.pending_seed.take(),
        upscale_cache.pending_sea_level.take(),
    ) else {
        return;
    };

    let Some(ref heightmap) = isostasy_cache.heightmap else {
        upscale_cache.state = FbmState::Idle;
        return;
    };

    let _ = bridge.commands_tx.send(SolverCommand::RunFbmUpscale {
        coarse: heightmap.clone(),
        sea_level,
        seed,
        config,
    });
}

/// Dispatch a pending erosion command to the solver thread.
fn dispatch_erosion(
    mut erosion_cache: ResMut<ErosionCache>,
    bridge: ResMut<SolverBridge>,
    upscale_cache: Res<UpscaleCache>,
) {
    let (Some(config), Some(seed)) =
        (erosion_cache.pending_config.take(), erosion_cache.pending_seed.take())
    else {
        return;
    };

    let Some(ref heightmap) = upscale_cache.heightmap else {
        erosion_cache.state = ErosionState::Idle;
        return;
    };

    let _ = bridge.commands_tx.send(SolverCommand::RunErosion {
        heightmap: heightmap.clone(),
        config,
        seed,
    });
}

/// Dispatch a pending flow computation to the solver thread.
fn dispatch_flow(
    mut flow_cache: ResMut<FlowCache>,
    bridge: ResMut<SolverBridge>,
    erosion_cache: Res<ErosionCache>,
    upscale_cache: Res<UpscaleCache>,
) {
    let Some(config) = flow_cache.pending_config.take() else {
        return;
    };

    // Prefer eroded heightmap, fall back to upscale
    let heightmap = erosion_cache.heightmap.as_ref().or(upscale_cache.heightmap.as_ref());

    let Some(heightmap) = heightmap else {
        flow_cache.state = FlowState::Idle;
        return;
    };

    let _ = bridge
        .commands_tx
        .send(SolverCommand::RunFlowComputation { heightmap: heightmap.clone(), config });
}

/// Re-extract rivers when thresholds change (runs on main thread, fast).
fn update_river_extraction(mut flow_cache: ResMut<FlowCache>) {
    if !flow_cache.rivers_dirty {
        return;
    }

    let Some(ref result) = flow_cache.result else {
        return;
    };

    let w = result.accumulation.width;
    let h = result.accumulation.height;
    let config = flow_cache.river_config.clone();

    let rivers = ymir_core::terrain::flow::extract_rivers(result, &config, w, h);
    flow_cache.rivers = Some(rivers);
    flow_cache.rivers_dirty = false;
}
