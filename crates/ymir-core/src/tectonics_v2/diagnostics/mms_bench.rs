//! Runtime MMS convergence benchmarks included in the Step-1 report.
//!
//! The baseline scenario at Step 1 does not fully exercise the
//! power-law rheology (the placeholder forcing is too gentle to
//! spread `ε̇_II` far above `ε̇_min`). To provide visible evidence
//! that the Newton path and the variable-η assembly are correct,
//! the binary runs two lightweight MMS convergence checks at report
//! time and reports the observed order-of-accuracy.
//!
//! The two manufactured solutions are identical to those used in
//! `tests/v2_stokes_mms_variable_eta.rs` and
//! `tests/v2_stokes_mms.rs` — factored out here so the binary doesn't
//! duplicate their assembly.

use std::f64::consts::PI;

use crate::tectonics_v2::field::Field2D;
use crate::tectonics_v2::rheology::{self, StrainRate, ViscosityLaw};
use crate::tectonics_v2::stokes::nonlinear_solver::{
    NewtonConfig, NewtonSolver, NonlinearSolver,
};
use crate::tectonics_v2::stokes::operator::apply_momentum;
use crate::tectonics_v2::stokes::solver::ConjugateGradient;
use crate::tectonics_v2::stokes::{solve_sheet, Grid, SheetConfig};

#[derive(Clone, Debug)]
pub struct MmsResults {
    /// Linear thin-sheet MMS with `η = 1` everywhere.
    pub const_eta: MmsSeries,
    /// Linear thin-sheet MMS with prescribed variable η.
    pub variable_eta: MmsSeries,
    /// Nonlinear Newton convergence test on a target-generated RHS
    /// at `n = 3` (single grid; reports the iteration trace tail).
    pub newton_tail: NewtonTail,
    /// Step-2 addition: GPE staggered-flux convergence on smooth S.
    pub gpe: MmsSeries,
}

#[derive(Clone, Debug, Default)]
pub struct MmsSeries {
    pub sizes: Vec<usize>,
    pub errors: Vec<f64>,
    pub slopes: Vec<f64>,
}

impl MmsSeries {
    pub fn final_slope(&self) -> Option<f64> {
        self.slopes.last().copied()
    }
}

#[derive(Clone, Debug, Default)]
pub struct NewtonTail {
    pub size: usize,
    pub residuals: Vec<f64>,
    pub outer_iters: u32,
    pub converged: bool,
}

/// Run all MMS series and the Newton-tail check. All computations
/// are small (128² max) and deterministic.
pub fn run_all() -> MmsResults {
    MmsResults {
        const_eta: mms_const_eta(&[16, 32, 64, 128]),
        variable_eta: mms_variable_eta(&[32, 64, 128]),
        newton_tail: newton_tail_at_n3(32),
        gpe: mms_gpe(&[32, 64, 128]),
    }
}

fn mms_gpe(sizes: &[usize]) -> MmsSeries {
    let mut series = MmsSeries::default();
    series.sizes = sizes.to_vec();
    for &n in sizes {
        series.errors.push(gpe_smooth_error(n));
    }
    for w in series.errors.windows(2) {
        series.slopes.push((w[0] / w[1]).log2());
    }
    series
}

fn gpe_smooth_error(n: usize) -> f64 {
    use std::f64::consts::TAU;
    use crate::tectonics_v2::field::PeriodicIndex;
    use crate::tectonics_v2::forcing::{BodyForce, GpeForce, SimulationState, VectorField};

    let nx = n;
    let ny = n;
    let dx = 1.0 / nx as f64;
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);
    let ar = 2.0;
    let alpha = 0.1;
    let k = TAU;

    let mut s = Field2D::new(nx, ny);
    for j in 0..ny {
        for i in 0..nx {
            let x = (i as f64 + 0.5) * dx;
            let y = (j as f64 + 0.5) * dx;
            s.set(i, j, 1.0 + alpha * (k * x).sin() * (k * y).cos());
        }
    }
    let mut fx = Field2D::new(nx, ny);
    let mut fy = Field2D::new(nx, ny);
    let st = SimulationState { nx, ny, dx, dy: dx, idx_x: &idx_x, idx_y: &idx_y, s: &s };
    GpeForce::with_ar(ar).accumulate(
        &st,
        &mut VectorField { fx: &mut fx, fy: &mut fy },
    );

    let mut sq = 0.0_f64;
    let mut count = 0_usize;
    for j in 0..ny {
        for i in 0..nx {
            let xf = i as f64 * dx;
            let yf = (j as f64 + 0.5) * dx;
            let s_at = 1.0 + alpha * (k * xf).sin() * (k * yf).cos();
            let dsdx = alpha * k * (k * xf).cos() * (k * yf).cos();
            let analytic = -ar * s_at * dsdx;
            sq += (fx.data()[j * nx + i] - analytic).powi(2);
            count += 1;
        }
    }
    (sq / count as f64).sqrt()
}

fn mms_const_eta(sizes: &[usize]) -> MmsSeries {
    let mut series = MmsSeries::default();
    series.sizes = sizes.to_vec();
    for &n in sizes {
        series.errors.push(const_eta_error(n));
    }
    for w in series.errors.windows(2) {
        series.slopes.push((w[0] / w[1]).log2());
    }
    series
}

fn const_eta_error(n: usize) -> f64 {
    let nx = n;
    let ny = n;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let grid = Grid::new(nx, ny, dx, dy);
    let eta = Field2D::filled(nx, ny, 1.0);
    let mut fx = vec![0.0; nx * ny];
    let mut fy = vec![0.0; nx * ny];
    let mut vx_ex = vec![0.0; nx * ny];
    let mut vy_ex = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let xfx = i as f64 * dx;
            fx[j * nx + i] = 8.0 * PI * PI * (2.0 * PI * xfx).sin();
            vx_ex[j * nx + i] = (2.0 * PI * xfx).sin();
            let yfy = j as f64 * dy;
            fy[j * nx + i] = 8.0 * PI * PI * (2.0 * PI * yfy).sin();
            vy_ex[j * nx + i] = (2.0 * PI * yfy).sin();
        }
    }
    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    let mut cfg = SheetConfig::default();
    cfg.tol = 1.0e-12;
    cfg.max_iter = 10_000;
    let _ = solve_sheet(&grid, &eta, &fx, &fy, &mut vx, &mut vy, &cfg);
    rms_err(&vx, &vy, &vx_ex, &vy_ex)
}

fn mms_variable_eta(sizes: &[usize]) -> MmsSeries {
    let mut series = MmsSeries::default();
    series.sizes = sizes.to_vec();
    for &n in sizes {
        series.errors.push(variable_eta_error(n));
    }
    for w in series.errors.windows(2) {
        series.slopes.push((w[0] / w[1]).log2());
    }
    series
}

fn variable_eta_error(n: usize) -> f64 {
    let nx = n;
    let ny = n;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let grid = Grid::new(nx, ny, dx, dy);
    let mut eta = Field2D::new(nx, ny);
    let mut fx = vec![0.0; nx * ny];
    let mut fy = vec![0.0; nx * ny];
    let mut vx_ex = vec![0.0; nx * ny];
    let mut vy_ex = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let xc = (i as f64 + 0.5) * dx;
            let yc = (j as f64 + 0.5) * dy;
            eta.set(
                i, j,
                1.0 + 0.5 * (2.0 * PI * xc).sin() * (2.0 * PI * yc).cos(),
            );
            let xf = i as f64 * dx;
            let yf = (j as f64 + 0.5) * dy;
            vx_ex[j * nx + i] = (2.0 * PI * xf).sin() * (2.0 * PI * yf).cos();
            fx[j * nx + i] = 8.0 * PI * PI * (2.0 * PI * xf).sin() * (2.0 * PI * yf).cos()
                - 4.0 * PI * PI * (4.0 * PI * xf).cos() * (2.0 * PI * yf).cos().powi(2);
            let xf2 = (i as f64 + 0.5) * dx;
            let yf2 = j as f64 * dy;
            vy_ex[j * nx + i] = -(2.0 * PI * xf2).cos() * (2.0 * PI * yf2).sin();
            fy[j * nx + i] = -8.0 * PI * PI * (2.0 * PI * xf2).cos() * (2.0 * PI * yf2).sin()
                - 2.0 * PI * PI * (4.0 * PI * xf2).sin() * (4.0 * PI * yf2).sin();
        }
    }
    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    let mut cfg = SheetConfig::default();
    cfg.tol = 1.0e-12;
    cfg.max_iter = 10_000;
    let _ = solve_sheet(&grid, &eta, &fx, &fy, &mut vx, &mut vy, &cfg);
    rms_err(&vx, &vy, &vx_ex, &vy_ex)
}

fn newton_tail_at_n3(n: usize) -> NewtonTail {
    let nx = n;
    let ny = n;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let grid = Grid::new(nx, ny, dx, dy);
    let mut law = ViscosityLaw::default();
    law.n = 3.0;
    // Target velocity and RHS (same layout as v2_newton_convergence).
    let mut vx_t = vec![0.0; nx * ny];
    let mut vy_t = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let xf = i as f64 * dx;
            let yf = (j as f64 + 0.5) * dy;
            vx_t[j * nx + i] = 0.3 * (2.0 * PI * xf).sin() * (2.0 * PI * yf).cos();
            let xf2 = (i as f64 + 0.5) * dx;
            let yf2 = j as f64 * dy;
            vy_t[j * nx + i] = 0.2 * (2.0 * PI * xf2).cos() * (2.0 * PI * yf2).sin();
        }
    }
    crate::tectonics_v2::stokes::nullspace::project_velocity(&mut vx_t, &mut vy_t);
    let sr = StrainRate::compute(
        nx, ny, dx, dy, &grid.idx_x, &grid.idx_y, &vx_t, &vy_t,
    );
    let eta = rheology::build_eta_field(&law, &sr.eps_ii_center);
    let mut rhs_x = vec![0.0; nx * ny];
    let mut rhs_y = vec![0.0; nx * ny];
    apply_momentum(&grid, &eta, &vx_t, &vy_t, &mut rhs_x, &mut rhs_y);
    crate::tectonics_v2::stokes::nullspace::project_velocity(&mut rhs_x, &mut rhs_y);

    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    let mut cfg = NewtonConfig::default();
    cfg.rel_tol = 1.0e-8;
    cfg.linear_tol = 1.0e-10;
    cfg.max_outer_iters = 15;
    let newton = NewtonSolver::new(cfg);
    let cg = ConjugateGradient::new(cfg.linear_tol, cfg.linear_max_iter);
    let outcome = newton.solve(&grid, &law, &rhs_x, &rhs_y, &mut vx, &mut vy, &cg);
    NewtonTail {
        size: n,
        residuals: outcome.trace().residuals.clone(),
        outer_iters: outcome.outer_iters(),
        converged: outcome.converged(),
    }
}

fn rms_err(vx: &[f64], vy: &[f64], vx_ex: &[f64], vy_ex: &[f64]) -> f64 {
    let n = (vx.len() + vy.len()) as f64;
    let s = vx
        .iter()
        .zip(vx_ex.iter())
        .chain(vy.iter().zip(vy_ex.iter()))
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>();
    (s / n).sqrt()
}

/// Render the MMS results as a markdown section.
pub fn render_markdown(results: &MmsResults) -> String {
    let mut s = String::new();
    s.push_str("## Discretisation validation (MMS convergence at report time)\n\n");
    s.push_str("The baseline run at Step 1 does not fully exercise power-law behaviour (the placeholder forcing saturates against `ε̇_min`). The following MMS convergence checks are run **every time the report is generated** so the discretisation remains visibly verified.\n\n");

    s.push_str("### Constant η (Picard path, Step 0 operator)\n\n");
    s.push_str("Manufactured solution `v = (sin(2πx), sin(2πy))`, `η = 1`.\n\n");
    s.push_str("| N | v_err RMS | slope to next |\n|---|---|---|\n");
    for (k, n) in results.const_eta.sizes.iter().enumerate() {
        let slope_str = results
            .const_eta
            .slopes
            .get(k)
            .map(|s| format!("{:.3}", s))
            .unwrap_or_else(|| "—".into());
        s.push_str(&format!(
            "| {} | {:.3e} | {} |\n",
            n, results.const_eta.errors[k], slope_str,
        ));
    }
    s.push_str(&format!(
        "\nFinal slope: `{:.3}` (expected ≥ 1.7; quadratic target = 2.0).\n\n",
        results.const_eta.final_slope().unwrap_or(f64::NAN),
    ));

    s.push_str("### Variable η (linear, prescribed η field)\n\n");
    s.push_str("Manufactured solution `v = (sin(2πx)cos(2πy), -cos(2πx)sin(2πy))`, `η(x,y) = 1 + 0.5·sin(2πx)·cos(2πy)`. Validates the η-variable Picard path used under Newton.\n\n");
    s.push_str("| N | v_err RMS | slope to next |\n|---|---|---|\n");
    for (k, n) in results.variable_eta.sizes.iter().enumerate() {
        let slope_str = results
            .variable_eta
            .slopes
            .get(k)
            .map(|s| format!("{:.3}", s))
            .unwrap_or_else(|| "—".into());
        s.push_str(&format!(
            "| {} | {:.3e} | {} |\n",
            n, results.variable_eta.errors[k], slope_str,
        ));
    }
    s.push_str(&format!(
        "\nFinal slope: `{:.3}` (expected ≥ 1.7).\n\n",
        results.variable_eta.final_slope().unwrap_or(f64::NAN),
    ));

    s.push_str("### GPE force (staggered `-Ar·∇(½S²)`, smooth S)\n\n");
    s.push_str("Smooth manufactured `S = 1 + 0.1·sin(2πx)·cos(2πy)` at `Ar = 2`. Validates the GPE discretisation introduced at Step 2 against the analytic `-Ar·S·∇S`.\n\n");
    s.push_str("| N | v_err RMS | slope to next |\n|---|---|---|\n");
    for (k, n) in results.gpe.sizes.iter().enumerate() {
        let slope_str = results
            .gpe
            .slopes
            .get(k)
            .map(|s| format!("{:.3}", s))
            .unwrap_or_else(|| "—".into());
        s.push_str(&format!(
            "| {} | {:.3e} | {} |\n",
            n, results.gpe.errors[k], slope_str,
        ));
    }
    s.push_str(&format!(
        "\nFinal slope: `{:.3}` (expected ≥ 1.7).\n\n",
        results.gpe.final_slope().unwrap_or(f64::NAN),
    ));

    s.push_str("### Nonlinear Newton tail (n = 3)\n\n");
    s.push_str(&format!(
        "Target-generated RHS on {}² grid. Newton outer iterations: `{}`{}.\n\n",
        results.newton_tail.size,
        results.newton_tail.outer_iters,
        if results.newton_tail.converged { "" } else { " — did NOT converge" },
    ));
    s.push_str("Residual trail:\n\n");
    s.push_str("```\n");
    for r in &results.newton_tail.residuals {
        s.push_str(&format!("  {:.3e}\n", r));
    }
    s.push_str("```\n\n");
    if results.newton_tail.residuals.len() >= 3 {
        let last = results.newton_tail.residuals.len();
        let r_km2 = results.newton_tail.residuals[last - 3];
        let r_km1 = results.newton_tail.residuals[last - 2];
        let r_k = results.newton_tail.residuals[last - 1];
        let r1 = if r_km1 > 0.0 { r_km2 / r_km1 } else { f64::NAN };
        let r2 = if r_k > 0.0 { r_km1 / r_k } else { f64::NAN };
        s.push_str(&format!(
            "Tail reductions: `{:.1}×` then `{:.1}×` (super-linear target: both ≥ 100×; strict quadratic requires an exact inner solve).\n\n",
            r1, r2,
        ));
    }
    s
}
