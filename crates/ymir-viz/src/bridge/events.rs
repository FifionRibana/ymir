//! Events sent from the solver thread back to Bevy.

use std::time::Duration;

use ymir_core::grid::GridF32;
use ymir_core::tectonics::plates::Plate;
use ymir_core::tectonics::solver::field::Field2D;
use ymir_core::tectonics::solver::workspace::StepStats;

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
    Failed {
        error: String,
    },
}
