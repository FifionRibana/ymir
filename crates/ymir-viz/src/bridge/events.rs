//! Events sent from the solver thread back to Bevy.

use std::time::Duration;

use ymir_core::erosion::hydraulic::ErosionStats;
use ymir_core::grid::GridF32;
use ymir_core::tectonics::plates::Plate;
use ymir_core::tectonics::solver::field::Field2D;
use ymir_core::tectonics::solver::workspace::StepStats;
use ymir_core::terrain::flow::FlowResult;

#[allow(dead_code)]
pub enum SolverEvent {
    Progress {
        step: usize,
        total_steps: usize,
        stats: StepStats,
    },
    Snapshot {
        step: usize,
        s_field: Field2D,
        plate_ids: Option<Vec<usize>>,
        plates: Option<Vec<Plate>>,
    },
    Completed {
        s_field: Field2D,
        plate_ids: Option<Vec<usize>>,
        plates: Option<Vec<Plate>>,
        elapsed: Duration,
        total_steps: usize,
    },
    FbmCompleted {
        heightmap: GridF32,
        slope: GridF32,
        elapsed: Duration,
    },
    ErosionProgress {
        completed: usize,
        total: usize,
    },
    ErosionSnapshot {
        heightmap: GridF32,
        completed: usize,
        total: usize,
    },
    ErosionCompleted {
        heightmap: GridF32,
        sediment: GridF32,
        stats: ErosionStats,
        elapsed: Duration,
    },
    FlowCompleted {
        result: FlowResult,
        elapsed: Duration,
    },
    Failed {
        error: String,
    },
}
