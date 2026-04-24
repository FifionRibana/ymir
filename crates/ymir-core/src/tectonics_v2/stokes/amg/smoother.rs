//! Red-Black Gauss-Seidel smoother — Step 8.5b Phase 3 (replaces
//! Step 8.5a's sequential SGS).
//!
//! # Why replace SGS?
//!
//! Sequential SGS is inherently serial — each row update reads
//! neighbours already updated in the current sweep. RBGS breaks
//! the dependency by **partitioning the unknowns into colours**
//! such that no two unknowns of the same colour share a non-zero
//! off-diagonal entry (see [`super::coloring`]). Within each
//! colour all rows can then be updated in parallel because their
//! updates don't depend on each other; between colours the sweep
//! proceeds sequentially.
//!
//! On structured grids with a 5-point / 9-point stencil this is a
//! textbook acceleration of Gauss-Seidel (Trottenberg-Oosterlee-
//! Schüller §5.3). On coarser AMG levels the greedy coloring adapts
//! automatically — any number of colours is valid, only parallel
//! efficiency degrades.
//!
//! # Convergence equivalence
//!
//! The iteration matrix of **colour-ordered** Gauss-Seidel differs
//! from row-ordered Gauss-Seidel only by a row permutation. On SPD
//! M-matrices (our regime) the spectral radius is unchanged, so
//! convergence per sweep is preserved to within the noise of a
//! different update order. The "within 5 %" agreement on the
//! Poisson gate (D2) is the measured invariant.
//!
//! # Bit-determinism
//!
//! Within a colour, parallel writes target **disjoint cells** (by
//! the coloring invariant), so the numeric result does not depend
//! on which worker processed which cell. The sequential colour
//! loop gives an implicit barrier. Same machine + same code → same
//! bits, independent of thread count.
//!
//! # Implementation note — unsafe raw-pointer write
//!
//! `&mut [f64]` cannot be shared across rayon workers under
//! Rust's borrow model, even when the writes target distinct
//! indices. The helper [`SyncSlicePtr`] wraps a raw `*mut f64`
//! and promises `Sync`; the safety argument sits in the coloring
//! invariant (same colour ⇒ no shared non-zero entry) and the
//! barrier-like semantics of `par_iter.for_each` (all workers for
//! colour `c` complete before colour `c+1` starts). The caller
//! of `rbgs_sweep` discharges this invariant by passing a colour
//! partition produced by [`super::coloring::greedy_coloring`].

use rayon::prelude::*;

use super::super::sparse_assembly::CsrMatrix;

/// Shared raw pointer into a mutable slice. Safe to send/share
/// across threads **if and only if** the caller guarantees that
/// simultaneous writers target distinct indices (and that no
/// thread writes while another is reading the same index).
///
/// Used by [`rbgs_sweep`]: within a colour the greedy partition
/// gives disjoint writer indices, and between colours the rayon
/// `for_each` barrier separates read and write phases.
struct SyncSlicePtr {
    ptr: *mut f64,
    len: usize,
}

unsafe impl Send for SyncSlicePtr {}
unsafe impl Sync for SyncSlicePtr {}

impl SyncSlicePtr {
    fn new(s: &mut [f64]) -> Self {
        Self { ptr: s.as_mut_ptr(), len: s.len() }
    }

    /// # Safety
    /// Caller guarantees no concurrent writer targets `i` and no
    /// concurrent writer targets any cell another reader is
    /// accessing within the same parallel block.
    #[inline(always)]
    unsafe fn read(&self, i: usize) -> f64 {
        debug_assert!(i < self.len);
        unsafe { *self.ptr.add(i) }
    }

    /// # Safety
    /// Caller guarantees no other thread is writing to or reading
    /// from `i` simultaneously.
    #[inline(always)]
    unsafe fn write(&self, i: usize, v: f64) {
        debug_assert!(i < self.len);
        unsafe { *self.ptr.add(i) = v; }
    }
}

/// Apply one symmetric RBGS sweep on `A · x = b`: forward over
/// colours in ascending order, then backward over colours in
/// descending order. Symmetry preserves SPD structure which CG
/// relies on.
///
/// `colors[c]` holds the sorted row indices of colour `c`, as
/// produced by [`super::coloring::greedy_coloring`]. The caller
/// is expected to have stored the partition alongside the level's
/// matrix so it is recomputed at most once per hierarchy build.
pub fn rbgs_sweep(a: &CsrMatrix, colors: &[Vec<usize>], b: &[f64], x: &mut [f64]) {
    let n = a.n_rows;
    debug_assert_eq!(a.n_cols, n, "RBGS requires a square matrix");
    debug_assert_eq!(b.len(), n);
    debug_assert_eq!(x.len(), n);

    if colors.is_empty() {
        // Fallback: a single-colour degenerate case (e.g. 1×1
        // matrix). Forward + backward passes collapse to a single
        // sequential update.
        sequential_gs_sweep(a, b, x);
        return;
    }

    let sync = SyncSlicePtr::new(x);

    // --- Forward pass: colours 0, 1, ..., C−1 ---
    for group in colors {
        group.par_iter().for_each(|&i| {
            update_row(a, b, &sync, i);
        });
    }

    // --- Backward pass: colours C−1, C−2, ..., 0 ---
    for group in colors.iter().rev() {
        group.par_iter().for_each(|&i| {
            update_row(a, b, &sync, i);
        });
    }
}

/// Single Gauss-Seidel update of row `i`: `x[i] ← (b[i] − Σ_{j≠i} a[i, j] x[j]) / a[i, i]`.
/// Accesses `x` through the shared pointer; see [`SyncSlicePtr`]
/// safety contract.
#[inline(always)]
fn update_row(a: &CsrMatrix, b: &[f64], sync: &SyncSlicePtr, i: usize) {
    let start = a.row_ptr[i];
    let end = a.row_ptr[i + 1];
    let mut acc = b[i];
    let mut diag = 0.0_f64;
    for k in start..end {
        let j = a.col_idx[k];
        if j == i {
            diag = a.values[k];
        } else {
            // Safe: j ≠ i (handled above) and by the coloring
            // invariant j has a different colour from i — no
            // concurrent writer targets j during this for_each.
            let xj = unsafe { sync.read(j) };
            acc -= a.values[k] * xj;
        }
    }
    if diag.abs() > 1e-300 {
        // Safe: all workers in this for_each write to distinct
        // indices (disjoint colour partition).
        unsafe { sync.write(i, acc / diag); }
    }
}

/// Sequential Gauss-Seidel fallback: one forward + one backward
/// sweep in index order, no parallelism. Used when no colouring
/// is provided (`colors.is_empty()`) or by tests that need a
/// stable reference behaviour.
pub fn sequential_gs_sweep(a: &CsrMatrix, b: &[f64], x: &mut [f64]) {
    let n = a.n_rows;
    debug_assert_eq!(a.n_cols, n);
    debug_assert_eq!(b.len(), n);
    debug_assert_eq!(x.len(), n);

    for i in 0..n {
        sequential_update_row(a, b, x, i);
    }
    for i in (0..n).rev() {
        sequential_update_row(a, b, x, i);
    }
}

#[inline]
fn sequential_update_row(a: &CsrMatrix, b: &[f64], x: &mut [f64], i: usize) {
    let start = a.row_ptr[i];
    let end = a.row_ptr[i + 1];
    let mut acc = b[i];
    let mut diag = 0.0_f64;
    for k in start..end {
        let j = a.col_idx[k];
        if j == i {
            diag = a.values[k];
        } else {
            acc -= a.values[k] * x[j];
        }
    }
    if diag.abs() > 1e-300 {
        x[i] = acc / diag;
    }
}

#[cfg(test)]
mod tests {
    use super::super::coloring::greedy_coloring;
    use super::super::strong_connections::build_poisson_laplacian_csr;
    use super::*;

    fn norm2(v: &[f64]) -> f64 {
        v.iter().map(|x| x * x).sum::<f64>().sqrt()
    }

    fn matvec(a: &CsrMatrix, x: &[f64], y: &mut [f64]) {
        a.apply(x, y);
    }

    #[test]
    fn rbgs_reduces_residual_monotonically_on_poisson() {
        use rand::{Rng, SeedableRng};
        use rand_chacha::ChaCha8Rng;
        let n = 8;
        let a = build_poisson_laplacian_csr(n);
        let colors = greedy_coloring(&a);
        let total = n * n;
        let mut rng = ChaCha8Rng::seed_from_u64(42);
        let mut b: Vec<f64> = (0..total).map(|_| rng.random::<f64>() * 2.0 - 1.0).collect();
        let mean: f64 = b.iter().sum::<f64>() / total as f64;
        for v in b.iter_mut() {
            *v -= mean;
        }

        let mut x = vec![0.0_f64; total];
        let mut ax = vec![0.0_f64; total];
        matvec(&a, &x, &mut ax);
        let r0: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
        let r0_norm = norm2(&r0);

        let mut prev = r0_norm;
        for sweep in 0..30 {
            rbgs_sweep(&a, &colors, &b, &mut x);
            matvec(&a, &x, &mut ax);
            let r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
            let r_norm = norm2(&r);
            assert!(
                r_norm <= prev + 1e-12,
                "sweep {}: residual increased from {:.3e} to {:.3e}",
                sweep,
                prev,
                r_norm,
            );
            prev = r_norm;
        }
        assert!(
            prev < r0_norm / 10.0,
            "30 RBGS sweeps gave only {:.3}× reduction",
            r0_norm / prev,
        );
    }

    #[test]
    fn rbgs_is_exact_on_diagonal_system() {
        let d_vals: Vec<f64> = (1..=5).map(|k| k as f64).collect();
        let b: Vec<f64> = (0..5).map(|k| (k as f64) * 0.5 + 1.0).collect();
        let a = CsrMatrix {
            n_rows: 5,
            n_cols: 5,
            row_ptr: vec![0, 1, 2, 3, 4, 5],
            col_idx: vec![0, 1, 2, 3, 4],
            values: d_vals.clone(),
        };
        let colors = greedy_coloring(&a);
        // Diagonal matrix ⇒ no off-diagonals ⇒ 1 colour, all rows.
        assert_eq!(colors.len(), 1);
        let mut x = vec![0.0; 5];
        rbgs_sweep(&a, &colors, &b, &mut x);
        for k in 0..5 {
            assert!((x[k] - b[k] / d_vals[k]).abs() < 1e-15);
        }
    }

    #[test]
    fn rbgs_damps_random_high_frequency_error() {
        // Smooths the error equation `A·x = 0` seeded with a RANDOM
        // high-freq-weighted field. A 2-colour RBGS on the 5-pt
        // Laplacian is an exact eigen-invariant for the pure
        // checkerboard mode (so we do NOT use it here — see report
        // §Phase 3 smoother note). A random seed averaged over many
        // modes decays as expected: ≥ 40 % reduction in 3 sweeps.
        use rand::{Rng, SeedableRng};
        use rand_chacha::ChaCha8Rng;
        let n = 8;
        let a = build_poisson_laplacian_csr(n);
        let colors = greedy_coloring(&a);
        let total = n * n;
        let mut rng = ChaCha8Rng::seed_from_u64(11);
        let mut x: Vec<f64> = (0..total).map(|_| rng.random::<f64>() * 2.0 - 1.0).collect();
        let m: f64 = x.iter().sum::<f64>() / total as f64;
        for v in x.iter_mut() {
            *v -= m;
        }
        let x0_norm = norm2(&x);
        let b = vec![0.0; total];
        for _ in 0..3 {
            rbgs_sweep(&a, &colors, &b, &mut x);
        }
        let xn = norm2(&x);
        assert!(
            xn < x0_norm * 0.6,
            "3 RBGS sweeps reduced random high-freq error by only {:.3}× ({:.3e} → {:.3e})",
            x0_norm / xn.max(1e-300),
            x0_norm,
            xn,
        );
    }

    #[test]
    fn rbgs_is_deterministic_across_runs() {
        let n = 6;
        let a = build_poisson_laplacian_csr(n);
        let colors = greedy_coloring(&a);
        let total = n * n;
        let b: Vec<f64> = (0..total).map(|k| (k as f64).sin()).collect();
        let mut x_a = vec![0.1_f64; total];
        let mut x_b = x_a.clone();
        for _ in 0..5 {
            rbgs_sweep(&a, &colors, &b, &mut x_a);
            rbgs_sweep(&a, &colors, &b, &mut x_b);
        }
        for k in 0..total {
            assert_eq!(x_a[k].to_bits(), x_b[k].to_bits());
        }
    }
}
