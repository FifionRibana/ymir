//! Commands sent from the Bevy main thread to the solver thread.

use ymir_core::tectonics::plates::Plate;
use ymir_core::tectonics::solver::config::TectonicsConfig;
use ymir_core::tectonics::solver::field::Field2D;
use ymir_core::tectonics::solver::traction::TractionField;

#[allow(dead_code, clippy::large_enum_variant)]
pub enum SolverCommand {
    RunTectonics {
        config: TectonicsConfig,
        plates: TractionField,
        plate_ids: Vec<usize>,
        plate_info: Vec<Plate>,
        initial_s: Field2D,
        grid_size: usize,
        dx: f64,
    },
    Cancel,
    Shutdown,
}
