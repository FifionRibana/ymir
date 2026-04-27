//! Startup-only continuation on the power-law exponent `n`.
//!
//! A Newton solver initialised from `v₀ = 0` on a strongly
//! shear-thinning system (`n = 3`) can diverge if the initial
//! residual sits far enough from the attractor. The classical
//! remedy is to warm-start by solving a sequence of gentler
//! problems: run the solver at `n = 1, 1.5, 2, 2.5, 3`, using each
//! converged iterate as the initial guess for the next.
//!
//! This is **startup-only**: applied exactly once, at the first
//! macro time step. Subsequent time steps inherit the previous
//! step's velocity and skip the ramp. Per-timestep continuation or
//! adaptive retry are out of scope for Step 1 — the fallback for
//! mid-run Newton trouble is the dual-track Picard path.

use super::super::cratonic::CratonicState;
use super::super::field::Field2D;
use super::super::presets::ContinuationConfig;
use super::super::rheology::ViscosityLaw;
use super::nonlinear_solver::{NewtonSolver, NonlinearOutcome, NonlinearSolver};
use super::operator::StokesGrid;
use super::solver::LinearSolver;

/// Result record for a full continuation ramp.
#[derive(Clone, Debug)]
pub struct ContinuationOutcome {
    /// One outcome per `n` value in the ramp, in order.
    pub sub_outcomes: Vec<(f64, NonlinearOutcome)>,
    /// True if every sub-solve converged.
    pub all_converged: bool,
    /// Total linear iterations across the ramp.
    pub linear_iters_total: u32,
}

/// Run the continuation ramp once and leave the final `n` equal to
/// the last value of `schedule.n_steps`. The input `law` is
/// temporarily mutated per sub-solve — callers get back the final `n`
/// through `law_out.n` on return.
///
/// If any sub-solve does not converge, the ramp stops early and
/// `all_converged` is `false`; the partial trace is still returned
/// for diagnostics.
pub fn run_continuation(
    grid: &StokesGrid,
    law_final: &ViscosityLaw,
    drag_diag: Option<&Field2D>,
    cratonic: Option<&CratonicState>,
    schedule: &ContinuationConfig,
    rhs_x: &[f64],
    rhs_y: &[f64],
    vx: &mut [f64],
    vy: &mut [f64],
    newton: &NewtonSolver,
    linear_solver: &dyn LinearSolver,
) -> ContinuationOutcome {
    let mut sub_outcomes = Vec::with_capacity(schedule.n_steps.len());
    let mut all_converged = true;
    let mut linear_iters_total = 0u32;

    for &n_current in &schedule.n_steps {
        let mut law_k = *law_final;
        law_k.n = n_current;
        let outcome = newton.solve(grid, &law_k, drag_diag, cratonic, rhs_x, rhs_y, vx, vy, linear_solver);
        if let NonlinearOutcome::Converged { linear_iters_total: lit, .. } = &outcome {
            linear_iters_total = linear_iters_total.saturating_add(*lit);
        }
        let converged = outcome.converged();
        sub_outcomes.push((n_current, outcome));
        if !converged {
            all_converged = false;
            break;
        }
    }

    ContinuationOutcome { sub_outcomes, all_converged, linear_iters_total }
}
