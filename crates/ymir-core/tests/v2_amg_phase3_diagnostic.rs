//! Step 8.5a Phase 3 diagnostic — mandatory measurements per
//! reviewer contract.
//!
//! Phase 3.0: formal multi-run AMG gates for step0/3/6/7
//! (re-measures Phase 2.7 bonus numbers with 5-run wallclock
//! mean+std; iter count is D9-deterministic so the headline
//! figure is a single value).
//!
//! Phase 3.1: step8 diagnostic "carte du territoire":
//!
//!  - η profile comparison step6 vs step8 (min, max, contrast,
//!    10-bin histogram). Establishes whether step8's operator is
//!    structurally harder than step6.
//!
//!  - AMG hierarchy structure comparison (levels built, per-level
//!    size ratio, per-level C-point fraction). Identifies whether
//!    coarsening is working.
//!
//!  - Instrumented V-cycle on step8: per-level residual norm
//!    before/after pre-smooth, at the coarse solve, after post-
//!    smooth. V-cycle reduction ratio ‖r_after‖/‖r_before‖ —
//!    > 0.5 signals V-cycle broken, < 0.1 signals V-cycle works
//!    but iters too few.
//!
//! No tuning, no fixes — just the map. Output goes to stderr via
//! `eprintln!` — invoke with `cargo test -- --nocapture`. Phase
//! 3.2 decisions follow the reviewer's interpretation of these
//! measurements.

use std::path::PathBuf;
use std::time::Instant;

use ymir_core::tectonics_v2::field::Field2D;
use ymir_core::tectonics_v2::stokes::amg::setup::{build_hierarchy, extract_diagonal_block};
use ymir_core::tectonics_v2::stokes::amg::smoother::sgs_sweep;
use ymir_core::tectonics_v2::stokes::amg::splitting::{classical_rs_splitting, CfType};
use ymir_core::tectonics_v2::stokes::amg::strong_connections::compute_strong_connections;
use ymir_core::tectonics_v2::stokes::amg::{AmgConfig, AmgPreconditioner};
use ymir_core::tectonics_v2::stokes::nullspace;
use ymir_core::tectonics_v2::stokes::operator::{apply_momentum, apply_tangent, StokesGrid, TangentContext};
use ymir_core::tectonics_v2::stokes::snapshot::{field_from_vec, LinearStokesSnapshot};
use ymir_core::tectonics_v2::stokes::solver::{ConjugateGradient, LinearSolver, SolverStats};
use ymir_core::tectonics_v2::stokes::sparse_assembly::{assemble_picard_csr, CsrMatrix};

// ---------------------------------------------------------------------------
// Shared fixtures
// ---------------------------------------------------------------------------

fn bench_data_path(case: &str) -> PathBuf {
    let manifest = std::env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR");
    PathBuf::from(manifest).join("..").join("..").join("bench_data").join(format!("{}.bin", case))
}

fn load_or_skip(case: &str) -> Option<LinearStokesSnapshot> {
    let path = bench_data_path(case);
    match LinearStokesSnapshot::load(&path) {
        Ok(s) => Some(s),
        Err(e) => {
            eprintln!("[skip] {}: {}", case, e);
            None
        }
    }
}

struct StokesReplay {
    grid: StokesGrid,
    eta_center: Field2D,
    drag_diag: Option<Field2D>,
    ctx: TangentContext,
    rhs_pack: Vec<f64>,
    cfg_tol: f64,
    cfg_max_iter: usize,
}

fn build_replay(snap: &LinearStokesSnapshot) -> StokesReplay {
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
    StokesReplay {
        grid,
        eta_center,
        drag_diag,
        ctx,
        rhs_pack,
        cfg_tol: snap.tol,
        cfg_max_iter: snap.max_iter,
    }
}

fn solve_amg(replay: &StokesReplay, snap: &LinearStokesSnapshot, cfg: AmgConfig) -> SolverStats {
    let a_picard =
        assemble_picard_csr(&replay.grid, &replay.eta_center, replay.drag_diag.as_ref());
    let precond = AmgPreconditioner::build(&a_picard, snap.n_cells(), cfg);
    let n = replay.grid.n_cells();
    let cg = ConjugateGradient::new(replay.cfg_tol, replay.cfg_max_iter);
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
    cg.solve(&mut matvec, &mut precond_fn, &b_scratch, &mut x_pack)
}

// ---------------------------------------------------------------------------
// 3.0 — Formal multi-run gates for step0/3/6/7
// ---------------------------------------------------------------------------

/// Target iter-count caps per reviewer contract (Phase 3):
///   step0_quiescent ≤ 10, step3_floor_yielding ≤ 15,
///   step6_voronoi ≤ 40, step7_slab_off ≤ 40.
/// Margin ×1.5-2 over Phase 2.7 bonus (4 / 9 / 9 / 8) to absorb
/// any measurement variance.
const STEP0_GATE: usize = 10;
const STEP3_GATE: usize = 15;
const STEP6_GATE: usize = 40;
const STEP7_GATE: usize = 40;

fn measure_multi_run(
    case: &str,
    gate_iters: usize,
    n_runs: usize,
) -> Option<(usize, f64, f64)> {
    let snap = load_or_skip(case)?;
    let replay = build_replay(&snap);
    let cfg = AmgConfig::default();

    // Iter count (D9-deterministic so one measurement is enough,
    // but we verify determinism by comparing the first vs last run).
    let mut iter_counts = Vec::with_capacity(n_runs);
    let mut wallclocks = Vec::with_capacity(n_runs);
    for _ in 0..n_runs {
        let t0 = Instant::now();
        let stats = solve_amg(&replay, &snap, cfg);
        let dt = t0.elapsed().as_secs_f64();
        iter_counts.push(stats.iterations);
        wallclocks.push(dt);
    }
    // D9 sanity: all iter counts identical.
    for &k in &iter_counts[1..] {
        assert_eq!(
            k, iter_counts[0],
            "{}: AMG iter count drifted across runs ({:?})",
            case, iter_counts
        );
    }
    let iters = iter_counts[0];
    let mean: f64 = wallclocks.iter().sum::<f64>() / n_runs as f64;
    let var: f64 = wallclocks.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n_runs as f64;
    let std = var.sqrt();
    eprintln!(
        "[gate] {:<24} AMG iters = {:>3}  wallclock {:.3} ± {:.3} ms  (gate ≤ {})",
        case,
        iters,
        mean * 1000.0,
        std * 1000.0,
        gate_iters,
    );
    assert!(
        iters <= gate_iters,
        "{} AMG iter count {} exceeds gate {}",
        case,
        iters,
        gate_iters,
    );
    Some((iters, mean, std))
}

#[test]
fn phase3_0_formal_gates_step0_step3_step6_step7() {
    eprintln!("=== Phase 3.0 — formal multi-run gates (5 runs each) ===");
    measure_multi_run("step0_quiescent", STEP0_GATE, 5);
    measure_multi_run("step3_floor_yielding", STEP3_GATE, 5);
    measure_multi_run("step6_voronoi", STEP6_GATE, 5);
    measure_multi_run("step7_slab_off", STEP7_GATE, 5);
}

// ---------------------------------------------------------------------------
// 3.1 — Diagnostic "carte du territoire" for step8
// ---------------------------------------------------------------------------

/// η field profile: min, max, contrast (max/min), and a 10-bin
/// histogram on log scale to reveal spatial distribution.
fn report_eta_profile(case: &str, snap: &LinearStokesSnapshot) {
    let eta = &snap.eta_center;
    let min = eta.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = eta.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mean: f64 = eta.iter().sum::<f64>() / eta.len() as f64;
    let contrast = if min > 0.0 { max / min } else { f64::INFINITY };
    eprintln!(
        "[eta ] {:<22} min={:.3e}  max={:.3e}  mean={:.3e}  contrast={:.2e}",
        case, min, max, mean, contrast
    );
    // Log-scale 10-bin histogram.
    if min > 0.0 && max > min {
        let log_min = min.ln();
        let log_max = max.ln();
        let step = (log_max - log_min) / 10.0;
        let mut bins = [0usize; 10];
        for &v in eta.iter() {
            let k = ((v.ln() - log_min) / step) as usize;
            let k = k.min(9);
            bins[k] += 1;
        }
        eprint!("       hist (log10 bins): ");
        for b in bins.iter() {
            eprint!("{:>5} ", b);
        }
        eprintln!();
    }
}

/// Report the hierarchy structure: level sizes, coarsening
/// ratios, C-point fractions per level.
fn report_hierarchy_structure(case: &str, label: &str, a: &CsrMatrix, cfg: AmgConfig) {
    let h = build_hierarchy(a.clone(), cfg);
    eprintln!("[hier] {:<22} {} hierarchy — {} levels", case, label, h.levels.len());
    for (idx, lvl) in h.levels.iter().enumerate() {
        let n = lvl.a.n_rows;
        let ratio = if idx == 0 {
            1.0
        } else {
            n as f64 / h.levels[idx - 1].a.n_rows as f64
        };
        // If we can still compute a splitting on this level (it's
        // not the coarsest), report C/F fractions too.
        let c_frac = if lvl.p.is_some() {
            // Re-compute splitting to get the stats (cheap at
            // setup scale); inherits determinism from the module.
            let strong = compute_strong_connections(&lvl.a, cfg.strong_connection_threshold);
            let cf = classical_rs_splitting(&strong);
            let c = cf.iter().filter(|&&c| c == CfType::Coarse).count();
            Some(c as f64 / n as f64)
        } else {
            None
        };
        let coarsest = if lvl.coarse_lu.is_some() { "  LU" } else { "" };
        eprintln!(
            "       level {:>2}: n={:>6}  ratio-to-prev={:.3}  C-frac={}{}",
            idx,
            n,
            ratio,
            c_frac.map(|f| format!("{:.3}", f)).unwrap_or_else(|| "  -- ".into()),
            coarsest,
        );
    }
}

fn infinity_norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x.abs()).fold(0.0_f64, f64::max)
}

/// Manually replay a V-cycle on a scalar hierarchy, instrumenting
/// residual norms at every transition. Uses only the public API
/// of the AMG submodules (no modification to production code).
fn instrumented_vcycle(
    case: &str,
    label: &str,
    snap: &LinearStokesSnapshot,
    cfg: AmgConfig,
) {
    // Build the Picard CSR, extract the u-u block, build hierarchy.
    let grid = StokesGrid::new(snap.nx, snap.ny, snap.dx, snap.dy);
    let eta = field_from_vec(snap.eta_center.clone(), snap.nx, snap.ny);
    let drag = snap
        .drag_diag
        .as_ref()
        .map(|v| field_from_vec(v.clone(), snap.nx, snap.ny));
    let a_picard = assemble_picard_csr(&grid, &eta, drag.as_ref());
    let n = snap.n_cells();
    let a_uu = extract_diagonal_block(&a_picard, 0, n);
    let hierarchy = build_hierarchy(a_uu.clone(), cfg);

    // Build a RHS that matches the vx half of the snapshot's RHS
    // (after gauge fixing). This is what the real CG sees on the
    // u half of the first AMG apply.
    let mut b = snap.rhs_vx.clone();
    nullspace::subtract_mean(&mut b);

    let n_levels = hierarchy.levels.len();
    eprintln!("[vcyc] {:<22} {} — V-cycle on u-u block, {} levels", case, label, n_levels);

    // x starts at zero.
    let x = vec![0.0f64; n];

    // Initial residual.
    let mut ax = vec![0.0f64; n];
    a_uu.apply(&x, &mut ax);
    let mut r = vec![0.0f64; n];
    for k in 0..n {
        r[k] = b[k] - ax[k];
    }
    let r0_norm = infinity_norm(&r);
    eprintln!("       level  0 | residual before ||r||_∞ = {:.3e}", r0_norm);

    // Manual V-cycle down: pre-smooth each level, restrict residual.
    let mut residuals: Vec<Vec<f64>> = vec![r.clone()];
    let mut xs: Vec<Vec<f64>> = vec![x.clone()];
    for k in 0..(n_levels - 1) {
        // Pre-smooth.
        for _ in 0..cfg.pre_smooth_sweeps {
            sgs_sweep(&hierarchy.levels[k].a, &residuals[k], &mut xs[k]);
        }
        // Fresh residual at this level.
        let mut ax_k = vec![0.0f64; residuals[k].len()];
        hierarchy.levels[k].a.apply(&xs[k], &mut ax_k);
        let fine_residual: Vec<f64> = residuals[k].iter().zip(ax_k.iter()).map(|(b, a)| b - a).collect();
        let r_pre = infinity_norm(&fine_residual);
        // Restrict to next level.
        let r_op = hierarchy.levels[k].r.as_ref().expect("non-coarsest has R");
        let mut b_next = vec![0.0f64; r_op.n_rows];
        r_op.apply(&fine_residual, &mut b_next);
        eprintln!(
            "       level {:>2} | after pre-smooth  ||r||_∞ = {:.3e}   restricted ||b_next||_∞ = {:.3e}",
            k,
            r_pre,
            infinity_norm(&b_next),
        );
        residuals.push(b_next);
        xs.push(vec![0.0f64; residuals[k + 1].len()]);
    }

    // Coarsest: LU solve.
    let coarsest = n_levels - 1;
    let lu = hierarchy.levels[coarsest]
        .coarse_lu
        .as_ref()
        .expect("coarsest has LU");
    let mut x_coarse = vec![0.0f64; residuals[coarsest].len()];
    lu.solve(&residuals[coarsest], &mut x_coarse);
    xs[coarsest] = x_coarse;
    eprintln!(
        "       level {:>2} | coarse LU solve — ||x_coarse||_∞ = {:.3e}",
        coarsest,
        infinity_norm(&xs[coarsest]),
    );

    // Up: prolongate + correct + post-smooth.
    for k in (0..coarsest).rev() {
        let p = hierarchy.levels[k].p.as_ref().expect("non-coarsest has P");
        let mut correction = vec![0.0f64; xs[k].len()];
        p.apply(&xs[k + 1], &mut correction);
        for i in 0..xs[k].len() {
            xs[k][i] += correction[i];
        }
        for _ in 0..cfg.post_smooth_sweeps {
            sgs_sweep(&hierarchy.levels[k].a, &residuals[k], &mut xs[k]);
        }
        // Residual after full V-cycle up through level k.
        let mut ax_post = vec![0.0f64; xs[k].len()];
        hierarchy.levels[k].a.apply(&xs[k], &mut ax_post);
        let post_r: Vec<f64> = residuals[k].iter().zip(ax_post.iter()).map(|(b, a)| b - a).collect();
        let r_post = infinity_norm(&post_r);
        eprintln!(
            "       level {:>2} | after prolongate+post-smooth ||r||_∞ = {:.3e}",
            k,
            r_post,
        );
    }

    // Final residual at level 0.
    let mut ax_final = vec![0.0f64; n];
    a_uu.apply(&xs[0], &mut ax_final);
    let r_final: Vec<f64> = b.iter().zip(ax_final.iter()).map(|(b, a)| b - a).collect();
    let r_final_norm = infinity_norm(&r_final);
    let reduction = r_final_norm / r0_norm;
    eprintln!(
        "       >>>> V-cycle reduction ratio ‖r_after‖/‖r_before‖ = {:.3e}  ({})",
        reduction,
        if reduction > 0.5 {
            "V-cycle INEFFICIENT (> 0.5)"
        } else if reduction < 0.1 {
            "V-cycle WORKS (< 0.1) — may just need more iters"
        } else {
            "V-cycle marginal"
        },
    );
}

#[test]
fn phase3_1_diagnostic_step6_vs_step8() {
    eprintln!("\n=== Phase 3.1 — diagnostic carte du territoire step8 ===\n");
    let snap6 = match load_or_skip("step6_voronoi") {
        Some(s) => s,
        None => return,
    };
    let snap8 = match load_or_skip("step8_activated") {
        Some(s) => s,
        None => return,
    };

    // η profile side by side.
    eprintln!("--- η profile ---");
    report_eta_profile("step6_voronoi", &snap6);
    report_eta_profile("step8_activated", &snap8);

    // Hierarchy structure side by side (u-u block — same pattern
    // as v-v by symmetry).
    eprintln!("\n--- AMG hierarchy structure (u-u block) ---");
    let cfg = AmgConfig::default();
    {
        let grid = StokesGrid::new(snap6.nx, snap6.ny, snap6.dx, snap6.dy);
        let eta = field_from_vec(snap6.eta_center.clone(), snap6.nx, snap6.ny);
        let drag = snap6
            .drag_diag
            .as_ref()
            .map(|v| field_from_vec(v.clone(), snap6.nx, snap6.ny));
        let a = assemble_picard_csr(&grid, &eta, drag.as_ref());
        let a_uu = extract_diagonal_block(&a, 0, snap6.n_cells());
        report_hierarchy_structure("step6_voronoi", "u-u", &a_uu, cfg);
    }
    {
        let grid = StokesGrid::new(snap8.nx, snap8.ny, snap8.dx, snap8.dy);
        let eta = field_from_vec(snap8.eta_center.clone(), snap8.nx, snap8.ny);
        let drag = snap8
            .drag_diag
            .as_ref()
            .map(|v| field_from_vec(v.clone(), snap8.nx, snap8.ny));
        let a = assemble_picard_csr(&grid, &eta, drag.as_ref());
        let a_uu = extract_diagonal_block(&a, 0, snap8.n_cells());
        report_hierarchy_structure("step8_activated", "u-u", &a_uu, cfg);
    }

    // Instrumented V-cycle on step6 (works — expect ratio < 0.1)
    // and step8 (plateau — expect ratio > 0.5 if V-cycle broken,
    // or similar to step6 if the problem is CG outer-loop depth
    // rather than V-cycle efficacy).
    eprintln!("\n--- V-cycle per-level residual trace ---");
    instrumented_vcycle("step6_voronoi", "reference (converges)", &snap6, cfg);
    eprintln!();
    instrumented_vcycle("step8_activated", "plateau case", &snap8, cfg);
}
