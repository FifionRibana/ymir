//! Commands sent from the Bevy main thread to the solver thread.

use ymir_core::tectonics::solver::config::TectonicsConfig;
use ymir_core::tectonics::solver::field::Field2D;
use ymir_core::tectonics::solver::tectonics::DynamicPlateContext;

#[allow(dead_code, clippy::large_enum_variant)]
pub enum SolverCommand {
    RunTectonics {
        config: TectonicsConfig,
        plate_ctx: DynamicPlateContext,
        initial_s: Field2D,
        grid_size: usize,
        dx: f64,
    },
    Cancel,
    Shutdown,
}
