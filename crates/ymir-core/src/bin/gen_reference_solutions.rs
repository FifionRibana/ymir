//! `gen_reference_solutions` — Step 8.5a Phase 4.2.
//!
//! Generates high-precision reference solutions `x_ref` for
//! step0/3/6/7 benchmark snapshots. The scalar-parity tests in
//! [`v2_amg_scalar_parity`] load these and compare both
//! `JacobiCG(tol=1e-8)` and `AmgCG(tol=1e-8)` outputs against them
//! with a threshold derived from the triangle-inequality formula
//! `C · κ · (tol_test + tol_ref)`.
//!
//! Reference strategy (per reviewer-approved option d.3):
//!
//! 1. Run AMG at `tol=1e-12, max_iter=10000` to get the finest
//!    solution achievable within budget.
//! 2. Record the achieved residual norm `tol_ref = ‖r‖/‖b‖` —
//!    for well-conditioned cases this reaches ~1e-12; for poorly
//!    conditioned cases it plateaus earlier.
//! 3. Empirically estimate κ(A) by running AMG at the test
//!    tolerance `tol_coarse = 1e-8` and computing
//!    `κ ≈ rel_diff(x_coarse, x_ref) / (tol_coarse + tol_ref)`.
//!
//! `AMG` is chosen over dense LU because (i) dense LU on the full
//! coupled `2N × 2N` Stokes operator would be 512 MB + 15-60 min
//! per case, infeasible; (ii) decoupled u-u/v-v LU would solve a
//! DIFFERENT system than production CG (which includes the
//! `apply_tangent` Newton tangent contribution).
//!
//! Usage:
//! ```bash
//! cargo run --release --bin gen_reference_solutions -- --all
//! ```

use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;

use ymir_core::tectonics_v2::field::Field2D;
use ymir_core::tectonics_v2::stokes::amg::{AmgConfig, AmgPreconditioner};
use ymir_core::tectonics_v2::stokes::nullspace;
use ymir_core::tectonics_v2::stokes::operator::{
    apply_momentum, apply_tangent, StokesGrid, TangentContext,
};
use ymir_core::tectonics_v2::stokes::precond::VelocityJacobi;
use ymir_core::tectonics_v2::stokes::snapshot::{
    field_from_vec, LinearStokesSnapshot, ReferenceSolution, SNAPSHOT_FORMAT_VERSION,
};
use ymir_core::tectonics_v2::stokes::solver::{ConjugateGradient, LinearSolver, SolverStats};
use ymir_core::tectonics_v2::stokes::sparse_assembly::assemble_picard_csr;

const ALL_CASES: &[&str] = &[
    "step0_quiescent",
    "step3_floor_yielding",
    "step6_voronoi",
    "step7_slab_off",
];
// step8_activated and step8_activated_128 intentionally omitted
// — α merge scope excludes them; scalar-parity testing of step8
// waits for Step 8.5a.2.

const REF_TOL: f64 = 1.0e-12;
const REF_MAX_ITER: usize = 10_000;
const KAPPA_TOL_COARSE: f64 = 1.0e-8;
const KAPPA_MAX_ITER: usize = 2_000;

struct Args {
    cases: Vec<String>,
    output_dir: PathBuf,
}

fn parse_args() -> Result<Args, String> {
    let mut a = Args {
        cases: Vec::new(),
        output_dir: PathBuf::from("bench_data").join("reference_solutions"),
    };
    let raw: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < raw.len() {
        match raw[i].as_str() {
            "--all" => a.cases = ALL_CASES.iter().map(|s| s.to_string()).collect(),
            "--case" => {
                i += 1;
                a.cases.push(raw[i].clone());
            }
            "--output-dir" => {
                i += 1;
                a.output_dir = PathBuf::from(&raw[i]);
            }
            "--help" | "-h" => {
                eprintln!("Usage: gen_reference_solutions [--all | --case NAME] [--output-dir PATH]");
                for c in ALL_CASES {
                    eprintln!("  {c}");
                }
                std::process::exit(0);
            }
            other => return Err(format!("unknown argument: {other}")),
        }
        i += 1;
    }
    if a.cases.is_empty() {
        return Err("pass --all or --case <name>; see --help".into());
    }
    Ok(a)
}

fn snapshot_path(case: &str) -> PathBuf {
    PathBuf::from("bench_data").join(format!("{case}.bin"))
}

struct Replay {
    grid: StokesGrid,
    eta_center: Field2D,
    drag_diag: Option<Field2D>,
    ctx: TangentContext,
    rhs_pack: Vec<f64>,
    b_norm: f64,
    n_cells: usize,
}

fn build_replay(snap: &LinearStokesSnapshot) -> Replay {
    let grid = StokesGrid::new(snap.nx, snap.ny, snap.dx, snap.dy);
    let eta_center = field_from_vec(snap.eta_center.clone(), snap.nx, snap.ny);
    let drag_diag = snap
        .drag_diag
        .as_ref()
        .map(|v| field_from_vec(v.clone(), snap.nx, snap.ny));
    let ctx = TangentContext {
        eta_center: eta_center.clone(),
        c_center: field_from_vec(
            snap.tangent_c_center.clone().unwrap_or_else(|| vec![0.0; snap.n_cells()]),
            snap.nx,
            snap.ny,
        ),
        exx_center: field_from_vec(
            snap.tangent_exx_center.clone().unwrap_or_else(|| vec![0.0; snap.n_cells()]),
            snap.nx,
            snap.ny,
        ),
        eyy_center: field_from_vec(
            snap.tangent_eyy_center.clone().unwrap_or_else(|| vec![0.0; snap.n_cells()]),
            snap.nx,
            snap.ny,
        ),
        exy_corner: field_from_vec(
            snap.tangent_exy_corner.clone().unwrap_or_else(|| vec![0.0; snap.n_cells()]),
            snap.nx,
            snap.ny,
        ),
    };
    let mut rhs_pack = Vec::with_capacity(2 * snap.n_cells());
    rhs_pack.extend_from_slice(&snap.rhs_vx);
    rhs_pack.extend_from_slice(&snap.rhs_vy);
    let b_norm: f64 = rhs_pack.iter().map(|v| v * v).sum::<f64>().sqrt();
    Replay {
        grid,
        eta_center,
        drag_diag,
        ctx,
        rhs_pack,
        b_norm,
        n_cells: snap.n_cells(),
    }
}

fn solve_jacobi(replay: &Replay, diag_vx: &[f64], diag_vy: &[f64], diag_floor: f64, tol: f64, max_iter: usize) -> (Vec<f64>, SolverStats) {
    let vjac = VelocityJacobi::from_diagonal(diag_vx, diag_vy, diag_floor);
    let n = replay.n_cells;
    let cg = ConjugateGradient::new(tol, max_iter);
    let mut x_pack = vec![0.0f64; 2 * n];
    let b_scratch = replay.rhs_pack.clone();
    let mut tmp_ax = vec![0.0; n];
    let mut tmp_ay = vec![0.0; n];
    let mut matvec = |v: &[f64], out: &mut [f64]| {
        let (vx_in, vy_in) = v.split_at(n);
        let (out_x, out_y) = out.split_at_mut(n);
        apply_momentum(
            &replay.grid,
            &replay.eta_center,
            replay.drag_diag.as_ref(),
            vx_in,
            vy_in,
            &mut tmp_ax,
            &mut tmp_ay,
        );
        out_x.copy_from_slice(&tmp_ax);
        out_y.copy_from_slice(&tmp_ay);
        apply_tangent(&replay.grid, &replay.ctx, vx_in, vy_in, out_x, out_y);
    };
    let mut precond = |r: &[f64], z: &mut [f64]| vjac.apply(r, z);
    let stats = cg.solve(&mut matvec, &mut precond, &b_scratch, &mut x_pack);
    {
        let (vx, vy) = x_pack.split_at_mut(n);
        nullspace::project_velocity(vx, vy);
    }
    (x_pack, stats)
}

fn solve_amg(replay: &Replay, snap: &LinearStokesSnapshot, tol: f64, max_iter: usize) -> (Vec<f64>, SolverStats) {
    let a_picard =
        assemble_picard_csr(&replay.grid, &replay.eta_center, replay.drag_diag.as_ref());
    let precond = AmgPreconditioner::build(&a_picard, snap.n_cells(), AmgConfig::default());
    let n = replay.n_cells;
    let cg = ConjugateGradient::new(tol, max_iter);
    let mut x_pack = vec![0.0f64; 2 * n];
    let b_scratch = replay.rhs_pack.clone();
    let mut tmp_ax = vec![0.0; n];
    let mut tmp_ay = vec![0.0; n];
    let mut matvec = |v: &[f64], out: &mut [f64]| {
        let (vx_in, vy_in) = v.split_at(n);
        let (out_x, out_y) = out.split_at_mut(n);
        apply_momentum(
            &replay.grid,
            &replay.eta_center,
            replay.drag_diag.as_ref(),
            vx_in,
            vy_in,
            &mut tmp_ax,
            &mut tmp_ay,
        );
        out_x.copy_from_slice(&tmp_ax);
        out_y.copy_from_slice(&tmp_ay);
        apply_tangent(&replay.grid, &replay.ctx, vx_in, vy_in, out_x, out_y);
    };
    let mut precond_fn = |r: &[f64], z: &mut [f64]| precond.apply(r, z);
    let stats = cg.solve(&mut matvec, &mut precond_fn, &b_scratch, &mut x_pack);
    // Project to canonical zero-mean gauge.
    {
        let (vx, vy) = x_pack.split_at_mut(n);
        nullspace::project_velocity(vx, vy);
    }
    (x_pack, stats)
}

fn rel_diff_max(a: &[f64], b: &[f64]) -> f64 {
    let norm_b = b.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if norm_b == 0.0 {
        return 0.0;
    }
    let diff = a
        .iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0_f64, f64::max);
    diff / norm_b
}

fn main() -> ExitCode {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("argument error: {e}");
            return ExitCode::from(2);
        }
    };

    if let Err(e) = std::fs::create_dir_all(&args.output_dir) {
        eprintln!("mkdir {}: {}", args.output_dir.display(), e);
        return ExitCode::FAILURE;
    }

    for case in &args.cases {
        let snap_path = snapshot_path(case);
        let snap = match LinearStokesSnapshot::load(&snap_path) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("[{case}] load snapshot failed: {e}");
                return ExitCode::FAILURE;
            }
        };
        let replay = build_replay(&snap);

        // --- Reference solve: AMG at tol=1e-12, max_iter=10000.
        let t0 = Instant::now();
        let (x_ref, stats_ref) = solve_amg(&replay, &snap, REF_TOL, REF_MAX_ITER);
        let t_ref = t0.elapsed().as_secs_f64();
        let tol_ref_achieved = stats_ref.final_residual / replay.b_norm.max(1.0);
        eprintln!(
            "[{case}] reference AMG: iters={}, converged={}, r_final/‖b‖={:.3e} ({:.1}s)",
            stats_ref.iterations,
            stats_ref.converged(),
            tol_ref_achieved,
            t_ref,
        );

        // --- κ estimate: run both AMG and Jacobi at tol=1e-8, take the max
        // rel_diff vs reference. Using the worst-case solver yields an
        // empirical κ that bounds both solvers' scalar-parity drift
        // uniformly, eliminating the need for different thresholds per
        // solver. (Theory: ‖x_cg − x*‖ ≲ C · κ(A) · tol, and the hidden
        // constant C depends on the Krylov path — different for Jacobi
        // and AMG. Measuring rel_diff/tol empirically absorbs C into κ.)
        let t0 = Instant::now();
        let (x_coarse_a, stats_coarse_a) = solve_amg(&replay, &snap, KAPPA_TOL_COARSE, KAPPA_MAX_ITER);
        let (x_coarse_j, stats_coarse_j) = solve_jacobi(
            &replay,
            &snap.diag_vx,
            &snap.diag_vy,
            snap.diag_floor,
            KAPPA_TOL_COARSE,
            KAPPA_MAX_ITER,
        );
        let t_coarse = t0.elapsed().as_secs_f64();
        let rel_a = rel_diff_max(&x_coarse_a, &x_ref);
        let rel_j = rel_diff_max(&x_coarse_j, &x_ref);
        let rel = rel_a.max(rel_j);
        let denom = (KAPPA_TOL_COARSE + tol_ref_achieved).max(1e-20);
        let kappa = rel / denom;
        eprintln!(
            "[{case}] κ estimate: AMG(it={}) rel={:.3e}  Jacobi(it={}) rel={:.3e}  max={:.3e}  κ≈{:.3e} ({:.1}s)",
            stats_coarse_a.iterations, rel_a,
            stats_coarse_j.iterations, rel_j,
            rel, kappa, t_coarse,
        );

        // --- Save.
        let n = snap.n_cells();
        let ref_sol = ReferenceSolution {
            format_version: SNAPSHOT_FORMAT_VERSION,
            case_label: case.clone(),
            nx: snap.nx,
            ny: snap.ny,
            x_vx: x_ref[..n].to_vec(),
            x_vy: x_ref[n..].to_vec(),
            tol_ref_achieved,
            kappa_estimated: kappa,
            ref_solver: "amg".into(),
            ref_iters: stats_ref.iterations,
            ref_converged: stats_ref.converged(),
        };
        let out_path = args.output_dir.join(format!("{case}.bin"));
        if let Err(e) = ref_sol.save(&out_path) {
            eprintln!("[{case}] save failed: {e}");
            return ExitCode::FAILURE;
        }
        eprintln!("[{case}] → {}  ✓", out_path.display());
    }

    ExitCode::SUCCESS
}
