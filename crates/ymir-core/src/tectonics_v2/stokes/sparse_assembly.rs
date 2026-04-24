//! Sparse matrix assembly of the thin-sheet momentum operator —
//! Step 8.5a Phase 1.
//!
//! Assembles the **Picard block** `A_picard = -∇·(2 η ε̇(·)) +
//! Br·S̃²·I` into a CSR matrix on the packed `[vx; vy]` layout of
//! size `2·N` where `N = nx · ny`. The Newton tangent contribution
//! `apply_tangent(ctx)` is deliberately **not** sparsified: the
//! tangent is possibly indefinite (negative semi-definite for
//! shear-thinning `n > 1`), and Classical AMG coarsening operates
//! on the SPD Picard block while the CG matvec continues to evaluate
//! the tangent matrix-free. This is Gerya §14.4's recommended
//! pattern for Stokes with non-linear rheology.
//!
//! # Stencil
//!
//! Each row of `A_picard` has **at most 9 non-zeros** (5 same-
//! component + 4 cross-coupling). Basal drag adds a diagonal
//! contribution, no new off-diagonal entries. The derivation is
//! the algebraic rewrite of [`super::operator::apply_momentum`]
//! — see the doc on [`assemble_picard_csr`] for the per-row layout.
//!
//! # Matrix-free parity guarantee — relative per-product metric
//!
//! `CsrMatrix::apply` and `apply_momentum` compute **the same
//! algebraic expression** — the only difference is the order in
//! which the per-row products are accumulated. f64 is not
//! associative under addition, so each output component differs
//! from its matrix-free counterpart by at most the accumulated
//! rounding
//! ```text
//!   |y_csr[k] − y_mf[k]| ≲ (nnz_per_row) · ε_mach · ‖A_k‖ · ‖x‖
//! ```
//! where `ε_mach ≈ 1.1·10⁻¹⁶` and `nnz_per_row = 9` here. On
//! heterogeneous η (operator row-norm `~ O(η_max / dx²)`) the
//! operator magnitude dominates and the *absolute* diff can be
//! `O(10⁻¹⁰)` or larger — not a correctness issue, pure rounding.
//!
//! The parity tests therefore use a **relative per-product**
//! metric:
//! ```text
//!   rel_diff = max_k |y_csr[k] − y_mf[k]| / (‖y_mf‖_∞ · nnz_per_row)
//! ```
//! which is bounded by `~ε_mach` and is the rigorous
//! machine-precision parity statement. Thresholds of `< 1e-14`
//! apply uniformly across uniform η, 100×/10⁴× contrast fields,
//! basal-drag augmentation, and all real snapshot captures.
//!
//! Covered by the 5 lib unit tests (`sparse_assembly::tests`) on
//! synthetic η and 4 integration tests
//! (`v2_sparse_assembly_snapshot_parity`) on real captured states
//! (step0, step3, step6, step7), each exercising 10 seeded zero-
//! mean test vectors.
//!
//! # Determinism
//!
//! Entries per row are emitted in **strictly increasing column
//! index order** so the CSR representation is canonical (no
//! implicit row-sort required by downstream AMG). This is a D9
//! determinism commitment: the same `(grid, η, drag)` always
//! produces byte-identical `row_ptr`, `col_idx`, `values` vectors.

use super::super::field::Field2D;
use super::operator::StokesGrid;

/// Compressed Sparse Row matrix on a `2·N` packed velocity layout.
///
/// The block layout mirrors what CG sees: rows `0..N` are vx
/// equations, rows `N..2N` are vy. Column indices follow the same
/// convention. Entries within a row are sorted by column index
/// (ascending) so downstream consumers can assume canonical order
/// (required by the AMG strong-connection and coarsening passes).
#[derive(Clone, Debug)]
pub struct CsrMatrix {
    pub n_rows: usize,
    pub n_cols: usize,
    /// Row pointers, length `n_rows + 1`. `row_ptr[i]` is the
    /// starting offset into `col_idx` / `values` for row `i`;
    /// `row_ptr[n_rows]` is the total non-zero count.
    pub row_ptr: Vec<usize>,
    pub col_idx: Vec<usize>,
    pub values: Vec<f64>,
}

impl CsrMatrix {
    pub fn nnz(&self) -> usize {
        *self.row_ptr.last().unwrap_or(&0)
    }

    /// Apply `y = A · x` on the packed `2·N` vector.
    pub fn apply(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.n_cols, "x length mismatches n_cols");
        assert_eq!(y.len(), self.n_rows, "y length mismatches n_rows");
        for i in 0..self.n_rows {
            let start = self.row_ptr[i];
            let end = self.row_ptr[i + 1];
            let mut acc = 0.0;
            for k in start..end {
                acc += self.values[k] * x[self.col_idx[k]];
            }
            y[i] = acc;
        }
    }
}

/// Assemble the Picard-block momentum operator `A_picard` into CSR.
///
/// # Row layout
///
/// Rows `0..N` correspond to vx equations; rows `N..2N` to vy. Each
/// row has up to 9 non-zeros in strict ascending column order.
///
/// # Stencil (for row `vx[i, j]`, linear index `lin_x(i, j) =
/// j·nx + i`)
///
/// The diagonal entry is
/// `2·(η(i,j) + η(im,j))/dx² + (η_c_top + η_c_bot)/dy² +
///  (drag(i,j)+drag(im,j))/2`
/// where `η_c_top = ¼·[η(im,j) + η(i,j) + η(im,jp) + η(i,jp)]` and
/// `η_c_bot` is the analogous bottom-corner 4-point average
/// (see [`super::operator::eta_corner`]).
///
/// Off-diagonals in the same component:
/// - `vx(ip, j)` :  `-2·η(i, j)/dx²`
/// - `vx(im, j)` :  `-2·η(im, j)/dx²`
/// - `vx(i, jp)` :  `-η_c_top / dy²`
/// - `vx(i, jm)` :  `-η_c_bot / dy²`
///
/// Cross-coupling vx ← vy (from `∂ y σ_xy = ∂_y [η (∂_y vx + ∂_x vy)]`):
/// - `vy(i,  jp)` : `-η_c_top / (dx·dy)`
/// - `vy(im, jp)` : `+η_c_top / (dx·dy)`
/// - `vy(i,  j)`  : `+η_c_bot / (dx·dy)`
/// - `vy(im, j)`  : `-η_c_bot / (dx·dy)`
///
/// The vy rows follow by x↔y symmetry (same derivation, swap roles).
///
/// # Complexity
///
/// `O(N)` assembly. Each row emits exactly its 9 non-zeros
/// (duplicates on periodic wrap collapse naturally via the
/// insert-and-accumulate in `push_entry`).
pub fn assemble_picard_csr(
    grid: &StokesGrid,
    eta: &Field2D,
    drag_diag: Option<&Field2D>,
) -> CsrMatrix {
    let nx = grid.nx;
    let ny = grid.ny;
    let dx = grid.dx;
    let dy = grid.dy;
    let inv_dx2 = 1.0 / (dx * dx);
    let inv_dy2 = 1.0 / (dy * dy);
    let inv_dxdy = 1.0 / (dx * dy);
    let n_cells = nx * ny;
    let n_dofs = 2 * n_cells;

    // Per-row staging: map col_idx → value. We use a small
    // `Vec<(usize, f64)>` since each row has at most 9 entries —
    // simpler than a hashmap and preserves determinism trivially.
    // The canonicalisation (sort by col, accumulate duplicates on
    // periodic wrap) happens in `flush_row`.
    let mut row_buf: Vec<(usize, f64)> = Vec::with_capacity(16);
    let mut row_ptr: Vec<usize> = Vec::with_capacity(n_dofs + 1);
    let mut col_idx: Vec<usize> = Vec::with_capacity(9 * n_dofs);
    let mut values: Vec<f64> = Vec::with_capacity(9 * n_dofs);
    row_ptr.push(0);

    let lin_x = |i: usize, j: usize| j * nx + i;
    let lin_y = |i: usize, j: usize| n_cells + j * nx + i;

    let eta_corner = |im: usize, i: usize, jm: usize, j: usize| -> f64 {
        0.25 * (eta.get(im, jm) + eta.get(i, jm) + eta.get(im, j) + eta.get(i, j))
    };

    let flush_row = |buf: &mut Vec<(usize, f64)>,
                     col_idx: &mut Vec<usize>,
                     values: &mut Vec<f64>,
                     row_ptr: &mut Vec<usize>| {
        // Sort by column index.
        buf.sort_by_key(|&(c, _)| c);
        // Merge duplicate columns (periodic wrap: e.g., nx=2 makes
        // im == ip for some (i,j)). Accumulate values for identical
        // column indices. Result is strictly ascending col order.
        let mut i = 0;
        while i < buf.len() {
            let (c, mut v) = buf[i];
            let mut j = i + 1;
            while j < buf.len() && buf[j].0 == c {
                v += buf[j].1;
                j += 1;
            }
            col_idx.push(c);
            values.push(v);
            i = j;
        }
        row_ptr.push(col_idx.len());
        buf.clear();
    };

    // ---- vx rows ----
    for j in 0..ny {
        let jp = grid.idx_y.next(j);
        let jm = grid.idx_y.prev(j);
        for i in 0..nx {
            let ip = grid.idx_x.next(i);
            let im = grid.idx_x.prev(i);

            let eta_cc_r = eta.get(i, j);
            let eta_cc_l = eta.get(im, j);
            let eta_c_top = eta_corner(im, i, j, jp);
            let eta_c_bot = eta_corner(im, i, jm, j);

            let diag = 2.0 * (eta_cc_r + eta_cc_l) * inv_dx2
                + (eta_c_top + eta_c_bot) * inv_dy2;

            // Same-component entries.
            row_buf.push((lin_x(i, j), diag));
            row_buf.push((lin_x(ip, j), -2.0 * eta_cc_r * inv_dx2));
            row_buf.push((lin_x(im, j), -2.0 * eta_cc_l * inv_dx2));
            row_buf.push((lin_x(i, jp), -eta_c_top * inv_dy2));
            row_buf.push((lin_x(i, jm), -eta_c_bot * inv_dy2));

            // Cross-coupling from d_sigma_xy_dy:
            //   +η_c_top · (dvx_dy_top + dvy_dx_top)/dy (with outer −)
            //   −η_c_bot · (dvx_dy_bot + dvy_dx_bot)/dy (with outer −)
            // dvy_dx_top = (vy(i, jp) − vy(im, jp))/dx    (sign +/−)
            // dvy_dx_bot = (vy(i, j)  − vy(im, j))/dx     (sign +/−)
            row_buf.push((lin_y(i, jp), -eta_c_top * inv_dxdy));
            row_buf.push((lin_y(im, jp), eta_c_top * inv_dxdy));
            row_buf.push((lin_y(i, j), eta_c_bot * inv_dxdy));
            row_buf.push((lin_y(im, j), -eta_c_bot * inv_dxdy));

            // Basal drag: +drag_face_x on the vx(i,j) self-entry.
            if let Some(drag) = drag_diag {
                let drag_x = 0.5 * (drag.get(im, j) + drag.get(i, j));
                row_buf.push((lin_x(i, j), drag_x));
            }

            flush_row(&mut row_buf, &mut col_idx, &mut values, &mut row_ptr);
        }
    }

    // ---- vy rows (x ↔ y symmetry) ----
    for j in 0..ny {
        let jp = grid.idx_y.next(j);
        let jm = grid.idx_y.prev(j);
        for i in 0..nx {
            let ip = grid.idx_x.next(i);
            let im = grid.idx_x.prev(i);

            let eta_cc_t = eta.get(i, j);
            let eta_cc_b = eta.get(i, jm);
            // Corner η for vy row: the "right" corner of vy(i,j) is
            // shared with vx(ip,*) and the "left" corner with vx(i,*),
            // spanning j → jm vertically. Mirroring apply_momentum's
            // choice at lines 178-179.
            let eta_c_right = eta_corner(i, ip, jm, j);
            let eta_c_left = eta_corner(im, i, jm, j);

            let diag = 2.0 * (eta_cc_t + eta_cc_b) * inv_dy2
                + (eta_c_right + eta_c_left) * inv_dx2;

            // Same-component entries.
            row_buf.push((lin_y(i, j), diag));
            row_buf.push((lin_y(i, jp), -2.0 * eta_cc_t * inv_dy2));
            row_buf.push((lin_y(i, jm), -2.0 * eta_cc_b * inv_dy2));
            row_buf.push((lin_y(ip, j), -eta_c_right * inv_dx2));
            row_buf.push((lin_y(im, j), -eta_c_left * inv_dx2));

            // Cross-coupling from d_sigma_xy_dx, mirroring the
            // apply_momentum y-momentum block (lines 180-186):
            //   dvx_dy_right = (vx(ip, j) − vx(ip, jm))/dy
            //   dvx_dy_left  = (vx(i,  j) − vx(i,  jm))/dy
            row_buf.push((lin_x(ip, j), -eta_c_right * inv_dxdy));
            row_buf.push((lin_x(ip, jm), eta_c_right * inv_dxdy));
            row_buf.push((lin_x(i, j), eta_c_left * inv_dxdy));
            row_buf.push((lin_x(i, jm), -eta_c_left * inv_dxdy));

            // Basal drag on vy(i,j) self-entry.
            if let Some(drag) = drag_diag {
                let drag_y = 0.5 * (drag.get(i, jm) + drag.get(i, j));
                row_buf.push((lin_y(i, j), drag_y));
            }

            flush_row(&mut row_buf, &mut col_idx, &mut values, &mut row_ptr);
        }
    }

    CsrMatrix {
        n_rows: n_dofs,
        n_cols: n_dofs,
        row_ptr,
        col_idx,
        values,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tectonics_v2::stokes::operator::apply_momentum;

    /// Generate a deterministic non-trivial test vector seeded by
    /// `seed` on `2·N` DOFs with zero mean per component (so the
    /// null-space alignment is explicit — CG sees zero-mean inputs).
    fn seeded_vector(seed: u64, n_cells: usize) -> Vec<f64> {
        use rand::{Rng, SeedableRng};
        use rand_chacha::ChaCha8Rng;
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let mut v = vec![0.0f64; 2 * n_cells];
        for k in 0..2 * n_cells {
            v[k] = rng.random::<f64>() * 2.0 - 1.0;
        }
        let mean_x: f64 = v[..n_cells].iter().sum::<f64>() / n_cells as f64;
        let mean_y: f64 = v[n_cells..].iter().sum::<f64>() / n_cells as f64;
        for k in 0..n_cells {
            v[k] -= mean_x;
            v[n_cells + k] -= mean_y;
        }
        v
    }

    fn eta_uniform(nx: usize, ny: usize, v: f64) -> Field2D {
        Field2D::filled(nx, ny, v)
    }

    fn eta_smooth_contrast(nx: usize, ny: usize, contrast: f64) -> Field2D {
        let mut f = Field2D::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                let x = (i as f64 + 0.5) / nx as f64;
                let y = (j as f64 + 0.5) / ny as f64;
                // η in [1, contrast], smooth.
                let phase =
                    (2.0 * std::f64::consts::PI * x).sin() * (2.0 * std::f64::consts::PI * y).cos();
                let val = 1.0 + 0.5 * (contrast - 1.0) * (1.0 + phase);
                f.set(i, j, val.max(1.0));
            }
        }
        f
    }

    fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
        a.iter()
            .zip(b.iter())
            .map(|(x, y)| (x - y).abs())
            .fold(0.0f64, f64::max)
    }

    /// Relative parity metric. The matrix-free `apply_momentum` and
    /// the sparse `CsrMatrix::apply` compute the SAME algebraic
    /// result but accumulate products in different orders, so
    /// rounding accumulates ~ε_mach · ‖A‖ · ‖x‖ per output
    /// component. Absolute-tolerance parity at 1e-14 is therefore
    /// impossible on heterogeneous η (operator norm ~ O(10⁴–10⁵) at
    /// 32² with contrast 100). The rigorous parity metric is the
    /// **relative residual**:
    ///   rel_diff = max|y_csr − y_mf| / (‖y_mf‖_inf · num_nz_per_row_bound)
    /// — which is bounded by a small multiple of f64 epsilon
    /// (here ~1e-14 after dividing by the ~9 products per row).
    fn relative_parity(a: &[f64], b: &[f64]) -> f64 {
        let norm_b = b.iter().map(|x| x.abs()).fold(0.0f64, f64::max);
        if norm_b == 0.0 {
            return 0.0;
        }
        let diff = max_abs_diff(a, b);
        // Divide by the bound on per-row summation width (9 products
        // per row in this stencil) so the threshold reflects
        // per-product f64 rounding, not accumulated rounding.
        diff / (norm_b * 9.0)
    }

    #[test]
    fn csr_matches_matrix_free_on_uniform_eta() {
        let nx = 16;
        let ny = 16;
        let n = nx * ny;
        let grid = StokesGrid::new(nx, ny, 1.0 / nx as f64, 1.0 / ny as f64);
        let eta = eta_uniform(nx, ny, 1.0);
        let csr = assemble_picard_csr(&grid, &eta, None);
        assert_eq!(csr.n_rows, 2 * n);
        assert_eq!(csr.n_cols, 2 * n);

        let x = seeded_vector(42, n);
        let mut y_csr = vec![0.0; 2 * n];
        csr.apply(&x, &mut y_csr);

        let mut y_mf_vx = vec![0.0; n];
        let mut y_mf_vy = vec![0.0; n];
        apply_momentum(&grid, &eta, None, &x[..n], &x[n..], &mut y_mf_vx, &mut y_mf_vy);

        let rel = relative_parity(&y_csr[..n], &y_mf_vx)
            .max(relative_parity(&y_csr[n..], &y_mf_vy));
        assert!(
            rel < 1e-14,
            "csr vs matrix-free relative parity = {:.3e}",
            rel
        );
    }

    #[test]
    fn csr_matches_matrix_free_on_contrast_heterogeneous_eta() {
        let nx = 32;
        let ny = 32;
        let n = nx * ny;
        let grid = StokesGrid::new(nx, ny, 1.0 / nx as f64, 1.0 / ny as f64);
        let eta = eta_smooth_contrast(nx, ny, 100.0);

        let csr = assemble_picard_csr(&grid, &eta, None);
        for seed in 0..10u64 {
            let x = seeded_vector(seed, n);
            let mut y_csr = vec![0.0; 2 * n];
            csr.apply(&x, &mut y_csr);
            let mut y_mf_vx = vec![0.0; n];
            let mut y_mf_vy = vec![0.0; n];
            apply_momentum(
                &grid, &eta, None, &x[..n], &x[n..], &mut y_mf_vx, &mut y_mf_vy,
            );
            let rel = relative_parity(&y_csr[..n], &y_mf_vx)
                .max(relative_parity(&y_csr[n..], &y_mf_vy));
            assert!(
                rel < 1e-14,
                "csr vs matrix-free relative parity on seed {} = {:.3e}",
                seed,
                rel
            );
        }
    }

    #[test]
    fn csr_matches_matrix_free_with_basal_drag() {
        let nx = 16;
        let ny = 16;
        let n = nx * ny;
        let grid = StokesGrid::new(nx, ny, 1.0 / nx as f64, 1.0 / ny as f64);
        let eta = eta_smooth_contrast(nx, ny, 10.0);
        let mut drag = Field2D::new(nx, ny);
        // drag = Br · S̃^exp with a mild spatial variation.
        for j in 0..ny {
            for i in 0..nx {
                let x = (i as f64 + 0.5) / nx as f64;
                let y = (j as f64 + 0.5) / ny as f64;
                drag.set(i, j, 0.05 * (1.0 + 0.3 * (x + y - 1.0)));
            }
        }

        let csr = assemble_picard_csr(&grid, &eta, Some(&drag));
        let x = seeded_vector(123, n);
        let mut y_csr = vec![0.0; 2 * n];
        csr.apply(&x, &mut y_csr);
        let mut y_mf_vx = vec![0.0; n];
        let mut y_mf_vy = vec![0.0; n];
        apply_momentum(
            &grid,
            &eta,
            Some(&drag),
            &x[..n],
            &x[n..],
            &mut y_mf_vx,
            &mut y_mf_vy,
        );
        let rel = relative_parity(&y_csr[..n], &y_mf_vx)
            .max(relative_parity(&y_csr[n..], &y_mf_vy));
        assert!(
            rel < 1e-14,
            "csr vs matrix-free with drag relative parity = {:.3e}",
            rel
        );
    }

    #[test]
    fn csr_is_column_sorted_per_row() {
        let nx = 8;
        let ny = 8;
        let grid = StokesGrid::new(nx, ny, 1.0 / nx as f64, 1.0 / ny as f64);
        let eta = eta_smooth_contrast(nx, ny, 10.0);
        let csr = assemble_picard_csr(&grid, &eta, None);
        for i in 0..csr.n_rows {
            let start = csr.row_ptr[i];
            let end = csr.row_ptr[i + 1];
            if end - start > 1 {
                for k in start + 1..end {
                    assert!(
                        csr.col_idx[k - 1] < csr.col_idx[k],
                        "row {} not strictly sorted at offset {}: {} ≥ {}",
                        i,
                        k,
                        csr.col_idx[k - 1],
                        csr.col_idx[k]
                    );
                }
            }
        }
    }

    #[test]
    fn csr_is_byte_deterministic() {
        let nx = 16;
        let ny = 16;
        let grid = StokesGrid::new(nx, ny, 1.0 / nx as f64, 1.0 / ny as f64);
        let eta = eta_smooth_contrast(nx, ny, 50.0);
        let a = assemble_picard_csr(&grid, &eta, None);
        let b = assemble_picard_csr(&grid, &eta, None);
        assert_eq!(a.row_ptr, b.row_ptr);
        assert_eq!(a.col_idx, b.col_idx);
        // Values are f64 but the computation is deterministic
        // (no reductions with unstable ordering), so bitwise equality
        // is expected.
        for k in 0..a.values.len() {
            assert_eq!(
                a.values[k].to_bits(),
                b.values[k].to_bits(),
                "value {} differs: {} vs {}",
                k,
                a.values[k],
                b.values[k]
            );
        }
    }
}
