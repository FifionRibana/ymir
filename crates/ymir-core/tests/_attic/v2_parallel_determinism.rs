//! Step 8.5b Phase 1 — cross-thread-count determinism of the
//! `parallel_reduce` primitives.
//!
//! Each test builds a **local** rayon thread pool via
//! `ThreadPoolBuilder::new().num_threads(n).build()` and executes
//! the reduction inside its `install()` scope. This avoids
//! touching the global pool (Q3 answer: "tests auto-contained,
//! pas besoin de --test-threads=1 global").
//!
//! The contract (stronger than D1 requires, achievable thanks to
//! the chunk-in-index-order pattern): `par_dot`, `par_norm2`,
//! `par_axpy`, `par_max_abs` produce `f64` outputs with identical
//! bit patterns across any of `num_threads ∈ {1, 2, 4, 8}`.

use rayon::ThreadPoolBuilder;
use ymir_core::tectonics_v2::stokes::parallel_reduce::{par_axpy, par_dot, par_max_abs, par_norm2};

const THREAD_COUNTS: &[usize] = &[1, 2, 4, 8];

/// Deterministic pseudo-input — `sin/cos`-based, not RNG, so this
/// test has no seed plumbing and its fixtures are byte-reproducible
/// on any machine.
fn fixture(n: usize) -> (Vec<f64>, Vec<f64>) {
    let a: Vec<f64> = (0..n).map(|k| (k as f64 * 0.011).sin()).collect();
    let b: Vec<f64> = (0..n).map(|k| (k as f64 * 0.017).cos() + 1.5).collect();
    (a, b)
}

fn run_on_pool<F, R>(num_threads: usize, f: F) -> R
where
    F: FnOnce() -> R + Send,
    R: Send,
{
    let pool =
        ThreadPoolBuilder::new().num_threads(num_threads).build().expect("build local rayon pool");
    pool.install(f)
}

#[test]
fn par_dot_is_bit_identical_across_thread_counts() {
    let n = 4_096;
    let (a, b) = fixture(n);
    let reference = run_on_pool(1, || par_dot(&a, &b));
    for &nt in THREAD_COUNTS {
        let got = run_on_pool(nt, || par_dot(&a, &b));
        assert_eq!(
            got.to_bits(),
            reference.to_bits(),
            "par_dot diverged at {nt} threads: got {got:.20e} vs reference {reference:.20e}"
        );
    }
}

#[test]
fn par_norm2_is_bit_identical_across_thread_counts() {
    let n = 4_096;
    let (a, _) = fixture(n);
    let reference = run_on_pool(1, || par_norm2(&a));
    for &nt in THREAD_COUNTS {
        let got = run_on_pool(nt, || par_norm2(&a));
        assert_eq!(got.to_bits(), reference.to_bits(), "par_norm2 diverged at {nt} threads");
    }
}

#[test]
fn par_axpy_is_bit_identical_across_thread_counts() {
    let n = 4_096;
    let (x, base_y) = fixture(n);
    let alpha = 0.3;
    let mut ref_y = base_y.clone();
    run_on_pool(1, || par_axpy(alpha, &x, &mut ref_y));
    for &nt in THREAD_COUNTS {
        let mut y = base_y.clone();
        run_on_pool(nt, || par_axpy(alpha, &x, &mut y));
        for (i, (got, expected)) in y.iter().zip(ref_y.iter()).enumerate() {
            assert_eq!(
                got.to_bits(),
                expected.to_bits(),
                "par_axpy diverged at {nt} threads cell {i}: got {got} vs reference {expected}"
            );
        }
    }
}

#[test]
fn par_max_abs_is_bit_identical_across_thread_counts() {
    let n = 4_096;
    let (a, _) = fixture(n);
    let reference = run_on_pool(1, || par_max_abs(&a));
    for &nt in THREAD_COUNTS {
        let got = run_on_pool(nt, || par_max_abs(&a));
        assert_eq!(got.to_bits(), reference.to_bits(), "par_max_abs diverged at {nt} threads");
    }
}

#[test]
fn par_dot_is_bit_identical_across_runs_on_same_pool() {
    // Sanity: repeatedly firing the same pool on the same input
    // must produce the same bits. Cross-pool is the stronger
    // claim above; this is a minimum floor.
    let n = 10_000;
    let (a, b) = fixture(n);
    let pool = ThreadPoolBuilder::new().num_threads(8).build().unwrap();
    let r0 = pool.install(|| par_dot(&a, &b));
    for _ in 0..50 {
        assert_eq!(pool.install(|| par_dot(&a, &b)).to_bits(), r0.to_bits());
    }
}

#[test]
fn par_dot_covers_small_and_boundary_sizes() {
    // Exercise the size classes that hit the CHUNK_COUNT=16 edge
    // cases: chunks of size 1, partial final chunk, exactly fills
    // a chunk. Each size must be bit-identical across thread
    // counts.
    for n in [15_usize, 16, 17, 31, 32, 33, 256, 257] {
        let (a, b) = fixture(n);
        let reference = run_on_pool(1, || par_dot(&a, &b));
        for &nt in THREAD_COUNTS {
            let got = run_on_pool(nt, || par_dot(&a, &b));
            assert_eq!(
                got.to_bits(),
                reference.to_bits(),
                "par_dot diverged at n={n}, threads={nt}"
            );
        }
    }
}
