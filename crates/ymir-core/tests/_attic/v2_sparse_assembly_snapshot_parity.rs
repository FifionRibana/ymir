//! Step 8.5a Phase 1 — sparse vs matrix-free parity on real
//! Stokes snapshots.
//!
//! The unit tests inside `sparse_assembly.rs` verify the stencil
//! on synthetic η fields. This integration test pins the invariant
//! end-to-end: load a real captured snapshot, rebuild the same
//! `eta_center` + `drag_diag` the physics-run Newton solve saw, and
//! verify that applying the assembled CSR matches `apply_momentum`
//! on multiple seeded test vectors to the same relative parity
//! threshold (~f64 epsilon after accumulating 9 products per row).
//!
//! Covers four of the six Stokes snapshots (the fast-loading ones);
//! the two large step8 cases are skipped to keep this test light
//! (full benchmark validates them).

use std::path::PathBuf;

use ymir_core::tectonics_v2::stokes::operator::{apply_momentum, StokesGrid};
use ymir_core::tectonics_v2::stokes::snapshot::{field_from_vec, LinearStokesSnapshot};
use ymir_core::tectonics_v2::stokes::sparse_assembly::assemble_picard_csr;

fn bench_data_path(case: &str) -> PathBuf {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")
        .expect("CARGO_MANIFEST_DIR not set");
    PathBuf::from(manifest_dir).join("..").join("..").join("bench_data").join(format!("{}.bin", case))
}

fn seeded_zero_mean(seed: u64, n_cells: usize) -> Vec<f64> {
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let mut v = vec![0.0f64; 2 * n_cells];
    for k in 0..2 * n_cells {
        v[k] = rng.random::<f64>() * 2.0 - 1.0;
    }
    let mx: f64 = v[..n_cells].iter().sum::<f64>() / n_cells as f64;
    let my: f64 = v[n_cells..].iter().sum::<f64>() / n_cells as f64;
    for k in 0..n_cells {
        v[k] -= mx;
        v[n_cells + k] -= my;
    }
    v
}

fn max_abs(x: &[f64]) -> f64 {
    x.iter().map(|v| v.abs()).fold(0.0_f64, f64::max)
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0_f64, f64::max)
}

fn check_snapshot_parity(case: &str) {
    let path = bench_data_path(case);
    let snap = match LinearStokesSnapshot::load(&path) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("[skip] {}: {}", case, e);
            return; // silently skip when snapshot absent (CI without bench_data)
        }
    };
    let grid = StokesGrid::new(snap.nx, snap.ny, snap.dx, snap.dy);
    let eta = field_from_vec(snap.eta_center.clone(), snap.nx, snap.ny);
    let drag = snap
        .drag_diag
        .as_ref()
        .map(|v| field_from_vec(v.clone(), snap.nx, snap.ny));

    let csr = assemble_picard_csr(&grid, &eta, drag.as_ref());
    let n = snap.n_cells();
    assert_eq!(csr.n_rows, 2 * n);
    assert_eq!(csr.n_cols, 2 * n);

    for seed in 0..10u64 {
        let x = seeded_zero_mean(seed, n);
        let mut y_csr = vec![0.0; 2 * n];
        csr.apply(&x, &mut y_csr);

        let mut y_mf_vx = vec![0.0; n];
        let mut y_mf_vy = vec![0.0; n];
        apply_momentum(
            &grid,
            &eta,
            drag.as_ref(),
            &x[..n],
            &x[n..],
            &mut y_mf_vx,
            &mut y_mf_vy,
        );
        let norm = max_abs(&y_mf_vx).max(max_abs(&y_mf_vy)).max(1e-300);
        let diff = max_abs_diff(&y_csr[..n], &y_mf_vx)
            .max(max_abs_diff(&y_csr[n..], &y_mf_vy));
        let rel = diff / (norm * 9.0);
        assert!(
            rel < 1e-14,
            "{} seed {}: rel parity {:.3e} (diff {:.3e}, norm {:.3e})",
            case,
            seed,
            rel,
            diff,
            norm
        );
    }
}

#[test]
fn parity_on_step0_quiescent() {
    check_snapshot_parity("step0_quiescent");
}

#[test]
fn parity_on_step3_floor_yielding() {
    check_snapshot_parity("step3_floor_yielding");
}

#[test]
fn parity_on_step6_voronoi() {
    check_snapshot_parity("step6_voronoi");
}

#[test]
fn parity_on_step7_slab_off() {
    check_snapshot_parity("step7_slab_off");
}
