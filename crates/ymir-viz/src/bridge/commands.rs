//! Commands sent from the Bevy main thread to the solver thread.

use ymir_core::erosion::hydraulic::ErosionConfig;
use ymir_core::grid::GridF32;
use ymir_core::seed::WorldSeed;
use ymir_core::tectonics::solver::config::TectonicsConfig;
use ymir_core::tectonics::solver::field::Field2D;
use ymir_core::tectonics::solver::tectonics::DynamicPlateContext;
use ymir_core::terrain::flow::FlowConfig;
use ymir_core::terrain::upscale::FbmUpscaleConfig;

#[allow(dead_code, clippy::large_enum_variant)]
pub enum SolverCommand {
    RunTectonics {
        config: TectonicsConfig,
        plate_ctx: DynamicPlateContext,
        initial_s: Field2D,
        grid_size: usize,
        dx: f64,
    },
    SingleStep {
        config: TectonicsConfig,
        plate_ctx: DynamicPlateContext,
        s_field: Field2D,
        grid_size: usize,
        dx: f64,
    },
    RunFbmUpscale {
        coarse: GridF32,
        sea_level: f32,
        seed: WorldSeed,
        config: FbmUpscaleConfig,
    },
    RunErosion {
        heightmap: GridF32,
        config: ErosionConfig,
        seed: WorldSeed,
    },
    RunFlowComputation {
        heightmap: GridF32,
        config: FlowConfig,
    },
    Cancel,
    Shutdown,
}
