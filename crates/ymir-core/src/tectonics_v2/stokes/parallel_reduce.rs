//! Step 8.5b D1 — deterministic parallel reductions.
//!
//! The design goal is a single invariant: **same input produces the
//! exact same `f64` bit pattern regardless of the rayon thread
//! count** (1, 2, 4, 8, 16, or any other). This is stronger than
//! what D1 of the Step 8.5b spec requires (which tolerates
//! variation across thread counts down to 1e-10 relative) and it
//! is achieved by the pattern below:
//!
//! 1. Work is split into a fixed number of **index-defined chunks**
//!    ([`CHUNK_COUNT`] = 16). A chunk is `[chunk_idx * chunk_size,
//!    (chunk_idx + 1) * chunk_size)`; the last chunk clamps to `n`.
//! 2. Each chunk's partial is computed **sequentially** in left-to-
//!    right index order (deterministic floating-point accumulation
//!    within the chunk, regardless of which rayon worker picks it).
//! 3. Rayon's `.collect::<Vec<_>>()` over an `IndexedParallelIterator`
//!    preserves index order, so the partials vector is
//!    deterministic.
//! 4. The final reduction over the partials is a **sequential**
//!    left-to-right sum (or max, for `par_max_abs`).
//!
//! Because the chunk boundary set and the final accumulation order
//! are fixed by index (not by which thread processed them), the
//! result does not depend on work-stealing order. The only choice
//! left to rayon is "which worker takes which chunk", and that
//! choice does not influence the numeric output.
//!
//! `CHUNK_COUNT` is fixed to 16 rather than derived from
//! `std::thread::available_parallelism()` so that reductions stay
//! cross-machine bit-identical. Scalar-parity across machines at
//! 1e-10 relative is guaranteed by this constant; bit-parity is
//! guaranteed on a given machine by the indexed scan.
//!
//! # Bit-parity versus Step 8.5a
//!
//! Switching from the 8.5a scalar `dot(a, b) = a.iter().zip(b).map(..).sum()`
//! to `par_dot` changes the *floating-point accumulation order* (16
//! chunk sub-sums rather than one running left-to-right sum). So
//! the output is **not** bit-identical to 8.5a output at the last
//! ULP. This is intentional — per D5 of the Step 8.5b spec,
//! bit-parity vs 8.5a is abandoned; scalar-parity at 1e-10 relative
//! is the new contract. The determinism the helpers preserve is
//! Step-8.5b-internal (runs of the same build on the same machine
//! agree byte-for-byte, regardless of thread count).

use rayon::prelude::*;

/// Fixed chunk count for deterministic parallel reductions.
///
/// Intentionally not a function of core count. See module docs for
/// the rationale.
pub const CHUNK_COUNT: usize = 16;

/// Deterministic parallel dot product.
///
/// Returns `sum_i a[i] * b[i]` computed as a chunked reduction that
/// is bit-identical regardless of thread count.
pub fn par_dot(a: &[f64], b: &[f64]) -> f64 {
    debug_assert_eq!(a.len(), b.len(), "par_dot: length mismatch");
    let n = a.len();
    if n == 0 {
        return 0.0;
    }
    let chunk_size = n.div_ceil(CHUNK_COUNT);
    let partials: Vec<f64> = (0..CHUNK_COUNT)
        .into_par_iter()
        .map(|chunk_idx| {
            let start = chunk_idx * chunk_size;
            let end = (start + chunk_size).min(n);
            let mut sum = 0.0_f64;
            for i in start..end {
                sum += a[i] * b[i];
            }
            sum
        })
        .collect();
    let mut total = 0.0_f64;
    for partial in &partials {
        total += *partial;
    }
    total
}

/// Deterministic parallel Euclidean norm (`sqrt(dot(a, a))`).
pub fn par_norm2(a: &[f64]) -> f64 {
    par_dot(a, a).sqrt()
}

/// Deterministic parallel AXPY: `y[i] += alpha * x[i]` for all `i`.
///
/// Each cell is computed independently (no reduction), so the
/// outcome is order-independent and therefore bit-identical across
/// thread counts by construction. `par_iter_mut` is preferred over
/// the chunk pattern here because there is nothing to accumulate.
pub fn par_axpy(alpha: f64, x: &[f64], y: &mut [f64]) {
    debug_assert_eq!(x.len(), y.len(), "par_axpy: length mismatch");
    y.par_iter_mut().zip(x.par_iter()).for_each(|(yi, xi)| {
        *yi += alpha * *xi;
    });
}

/// Deterministic parallel `max |a[i]|`.
///
/// Uses the same chunk pattern as [`par_dot`]. Returns `0.0` on an
/// empty slice.
pub fn par_max_abs(a: &[f64]) -> f64 {
    let n = a.len();
    if n == 0 {
        return 0.0;
    }
    let chunk_size = n.div_ceil(CHUNK_COUNT);
    let partials: Vec<f64> = (0..CHUNK_COUNT)
        .into_par_iter()
        .map(|chunk_idx| {
            let start = chunk_idx * chunk_size;
            let end = (start + chunk_size).min(n);
            let mut m = 0.0_f64;
            for i in start..end {
                let v = a[i].abs();
                if v > m {
                    m = v;
                }
            }
            m
        })
        .collect();
    let mut m = 0.0_f64;
    for partial in &partials {
        if *partial > m {
            m = *partial;
        }
    }
    m
}

#[cfg(test)]
mod tests {
    //! Unit tests for the reduction primitives. Cross-thread-count
    //! determinism is covered by the dedicated integration test in
    //! [`crates/ymir-core/tests/v2_parallel_determinism.rs`]; the
    //! tests here exercise correctness and edge cases only.

    use super::*;

    fn ramp(n: usize) -> Vec<f64> {
        (0..n).map(|k| 1.0 + 0.001 * (k as f64)).collect()
    }

    #[test]
    fn par_dot_matches_scalar_within_eps() {
        let n = 10_000;
        let a = ramp(n);
        let b: Vec<f64> = (0..n).map(|k| 2.0 - 0.0003 * (k as f64)).collect();
        let par = par_dot(&a, &b);
        let scalar: f64 = a.iter().zip(b.iter()).map(|(x, y)| x * y).sum();
        // Different accumulation orders — expect agreement at f64
        // precision but not ulp-for-ulp.
        let rel = ((par - scalar) / scalar).abs();
        assert!(rel < 1e-12, "par={par} scalar={scalar} rel={rel:.3e}");
    }

    #[test]
    fn par_dot_edge_cases() {
        assert_eq!(par_dot(&[], &[]), 0.0);
        assert_eq!(par_dot(&[3.0], &[4.0]), 12.0);
        let a = vec![1.0; 15];
        let b = vec![1.0; 15];
        assert!((par_dot(&a, &b) - 15.0).abs() < 1e-14);
        let c = vec![1.0; 17];
        let d = vec![1.0; 17];
        assert!((par_dot(&c, &d) - 17.0).abs() < 1e-14);
    }

    #[test]
    fn par_norm2_matches_scalar_within_eps() {
        let n = 1_024;
        let a = ramp(n);
        let par = par_norm2(&a);
        let scalar: f64 = a.iter().map(|x| x * x).sum::<f64>().sqrt();
        let rel = ((par - scalar) / scalar).abs();
        assert!(rel < 1e-12, "par={par} scalar={scalar}");
    }

    #[test]
    fn par_axpy_computes_y_plus_alpha_x() {
        let x = vec![1.0, 2.0, 3.0, 4.0];
        let mut y = vec![10.0, 20.0, 30.0, 40.0];
        par_axpy(0.5, &x, &mut y);
        assert_eq!(y, vec![10.5, 21.0, 31.5, 42.0]);
    }

    #[test]
    fn par_axpy_empty_noop() {
        let x: Vec<f64> = Vec::new();
        let mut y: Vec<f64> = Vec::new();
        par_axpy(2.0, &x, &mut y);
        assert!(y.is_empty());
    }

    #[test]
    fn par_max_abs_finds_extremum() {
        let a = vec![0.5, -3.2, 1.7, -0.4, 2.9];
        assert!((par_max_abs(&a) - 3.2).abs() < 1e-14);
    }

    #[test]
    fn par_max_abs_edge_cases() {
        assert_eq!(par_max_abs(&[]), 0.0);
        assert!((par_max_abs(&[-0.7]) - 0.7).abs() < 1e-14);
        let zeros = vec![0.0; 100];
        assert_eq!(par_max_abs(&zeros), 0.0);
    }

    #[test]
    fn par_dot_deterministic_across_repeated_invocations() {
        // Same input called 100 times must produce identical bits.
        let n = 4_096;
        let a = ramp(n);
        let b: Vec<f64> = (0..n).map(|k| (k as f64).sin()).collect();
        let r0 = par_dot(&a, &b);
        for _ in 0..99 {
            assert_eq!(par_dot(&a, &b).to_bits(), r0.to_bits());
        }
    }

    #[test]
    fn par_max_abs_deterministic_across_repeated_invocations() {
        let n = 4_096;
        let a: Vec<f64> = (0..n).map(|k| (k as f64).sin()).collect();
        let r0 = par_max_abs(&a);
        for _ in 0..99 {
            assert_eq!(par_max_abs(&a).to_bits(), r0.to_bits());
        }
    }
}
