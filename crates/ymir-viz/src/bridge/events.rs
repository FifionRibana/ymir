//! Events sent from the solver thread back to Bevy.

use std::time::Duration;
use ymir_core::tectonics::solver::field::Field2D;
use ymir_core::tectonics::solver::workspace::StepStats;

#[allow(dead_code)]
pub enum SolverEvent {
    Progress { step: usize, total_steps: usize, stats: StepStats },
    Snapshot { step: usize, s_field: Field2D },
    Completed { s_field: Field2D, elapsed: Duration, total_steps: usize },
    Failed { error: String },
}
