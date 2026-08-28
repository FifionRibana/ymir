//! Classical Algebraic Multigrid (Ruge-Stüben 1987) preconditioner
//! — Step 8.5a Phase 2+.
//!
//! # Architecture — Option B' (unknown-based, two scalar hierarchies)
//!
//! The thin-sheet momentum operator on the packed `[vx; vy]` 2·N
//! layout is scalar-per-unknown but has cross-coupling (`vx ↔ vy`
//! through the shear term of `∇·(2η ε̇)`). Per the pre-
//! implementation review (Q3 decision): **start with two
//! independent scalar AMG hierarchies**, one on the `u-u` block
//! (rows/cols `0..N`) and one on the `v-v` block (`N..2N`). The
//! cross-coupling is handled inside CG's matvec but not in the
//! preconditioner — the smoother's Gauss-Seidel sweeps see it
//! implicitly through the residual.
//!
//! If Phase 2.7 benchmarks show `step8_activated` cannot meet
//! its gate under Option B', escalate to a dedicated follow-up
//! step implementing point-based 2×2-block AMG (explicit scope
//! creep, out-of-scope for 8.5a).
//!
//! # Hierarchy setup — Classical Ruge-Stüben two-pass
//!
//! Per design review Q2:
//!
//! - Strong connections at threshold θ = 0.25 (configurable).
//! - C/F splitting via two-pass original algorithm (Briggs-
//!   Henson-McCormick Ch. 8.8). Not CLJP / PMIS (those target
//!   parallel execution; AMG 8.5a is single-threaded).
//! - Prolongation weights via the Ruge-Stüben classical formula.
//! - Deterministic tie-breaking by **lowest index wins**.
//! - Coarsen until `n_coarse ≤ min_coarse_unknowns` (default 50)
//!   or `max_levels` reached (default 7) — whichever first.
//!
//! # Smoother and coarse-grid solve
//!
//! - Pre- and post-smoother: Symmetric Gauss-Seidel, 2 sweeps
//!   each (D8 default). Symmetric (forward + backward) preserves
//!   SPD structure required by CG.
//! - Coarsest grid solve: direct LU (Doolittle with partial
//!   pivoting, reimplemented here per QA3 — no new `nalgebra`
//!   dependency).
//!
//! # Null-space handling (Option α)
//!
//! Every `AmgPreconditioner::apply` wraps a `project_velocity`
//! before entering the V-cycle and another after exit, mirroring
//! the Jacobi convention in [`super::precond::VelocityJacobi`].
//! This preserves the SPD structure CG relies on without
//! projecting the operator at setup (which would densify the
//! CSR). See the Phase 0 report §QA4 for the motivation.
//!
//! # Phase 2 implementation status — WIP
//!
//! - **2.0 (this commit)**: skeleton, types, module layout.
//! - **2.1**: `strong_connections.rs`.
//! - **2.2**: `splitting.rs` (Classical RS two-pass) with
//!   100-run byte-determinism test.
//! - **2.3**: `prolongation.rs` + `restriction.rs`.
//! - **2.4**: `smoother.rs` (SGS).
//! - **2.5**: `coarse_solve.rs` (Doolittle LU).
//! - **2.6**: `vcycle.rs`.
//! - **2.7**: CG dispatch (`LinearSolverConfig::AmgCG`) +
//!   benchmark gate (poisson_contrast_10000 ≤ 100 iters).

pub mod coarse_solve;
pub mod coloring;
pub mod fmg;
pub mod prolongation;
pub mod restriction;
pub mod setup;
pub mod smoother;
pub mod splitting;
pub mod strong_connections;
pub mod vcycle;

use super::sparse_assembly::CsrMatrix;

/// AMG configuration — reviewer-validated defaults per design Q2.
#[derive(Clone, Copy, Debug)]
pub struct AmgConfig {
    /// Strong-connection threshold `θ` in
    /// `|a_ij| ≥ θ · max_{k≠i} |a_ik|`. Standard Classical RS
    /// default is 0.25.
    pub strong_connection_threshold: f64,
    /// Hard cap on hierarchy depth; coarsening stops when reached
    /// even if `min_coarse_unknowns` not yet met. Guard against
    /// pathological inputs that refuse to coarsen.
    pub max_levels: usize,
    /// Coarsening stops when the coarse-grid size drops to this
    /// threshold. Default 50 — small enough that the LU direct
    /// solve is ~O(50³) = 125 k ops, negligible per V-cycle.
    pub min_coarse_unknowns: usize,
    /// Pre-smoother symmetric Gauss-Seidel sweep count per level.
    /// One sweep = one forward + one backward pass. Default 1
    /// (= 2 "passes" in the D8 wording, SGS-symmetric).
    pub pre_smooth_sweeps: u32,
    /// Post-smoother symmetric GS sweep count per level. Default 1.
    pub post_smooth_sweeps: u32,
}

impl Default for AmgConfig {
    fn default() -> Self {
        Self {
            strong_connection_threshold: 0.25,
            max_levels: 7,
            min_coarse_unknowns: 50,
            pre_smooth_sweeps: 1,
            post_smooth_sweeps: 1,
        }
    }
}

/// Single scalar AMG hierarchy, applied to an `N × N` SPD block.
///
/// Phase 2.0: struct skeleton only; the level vector is populated
/// by the `setup` module in Phase 2.6 once the C/F splitting,
/// prolongation, and coarse-solve machinery lands.
#[derive(Debug)]
pub struct AmgHierarchy {
    /// Level 0 is the fine-grid operator (a scalar block of the
    /// original 2N·2N `A_picard`). Level `k+1` is constructed by
    /// `R_k · A_k · P_k`. The last entry is the coarse-grid
    /// factorisation's source matrix.
    pub levels: Vec<AmgLevel>,
}

/// One level of the AMG hierarchy.
///
/// - `a` is the operator at this level.
/// - `p` / `r` are the prolongation (to finer) and restriction
///   (from finer) operators; both `None` on the coarsest level.
/// - `coarse_lu` is the direct-solve factorisation, `Some` only
///   on the coarsest level.
#[derive(Debug)]
pub struct AmgLevel {
    pub a: CsrMatrix,
    pub p: Option<CsrMatrix>,
    pub r: Option<CsrMatrix>,
    pub coarse_lu: Option<coarse_solve::LuFactorisation>,
    /// Step 8.5b Phase 3: algebraic greedy coloring of `a`, stored
    /// as `colors[c] = sorted row indices of colour c`. Enables
    /// Red-Black Gauss-Seidel smoothing (parallel within each
    /// colour). Empty on the coarsest level (where the direct LU
    /// solve bypasses the smoother).
    pub colors: Vec<Vec<usize>>,
}

/// AMG preconditioner on the `[vx; vy]` 2·N layout — Option B'.
///
/// Holds two independent scalar hierarchies. Each `apply` extracts
/// the vx and vy residual halves, runs V-cycle on each in isolation,
/// and re-packs. Null-space projection wraps the entry and exit.
#[derive(Debug)]
pub struct AmgPreconditioner {
    pub n_cells: usize,
    pub u_hierarchy: AmgHierarchy,
    pub v_hierarchy: AmgHierarchy,
    pub cfg: AmgConfig,
}

impl AmgPreconditioner {
    /// Construct the two-hierarchy preconditioner from the full
    /// `2N × 2N` `A_picard` by extracting the `u-u` and `v-v`
    /// scalar blocks and running `setup::build_hierarchy` on each.
    ///
    /// The cross-coupling blocks (`u-v`, `v-u`) are **discarded**
    /// for the preconditioner — Option B' of the design review
    /// Q3. They remain in the full CG matvec (via the sparse
    /// matrix directly); only the preconditioner is block-
    /// diagonal.
    pub fn build(a_picard: &CsrMatrix, n_cells: usize, cfg: AmgConfig) -> Self {
        debug_assert_eq!(a_picard.n_rows, 2 * n_cells);
        debug_assert_eq!(a_picard.n_cols, 2 * n_cells);

        // Extract the diagonal scalar blocks.
        let a_uu = setup::extract_diagonal_block(a_picard, 0, n_cells);
        let a_vv = setup::extract_diagonal_block(a_picard, n_cells, n_cells);

        let u_hierarchy = setup::build_hierarchy(a_uu, cfg);
        let v_hierarchy = setup::build_hierarchy(a_vv, cfg);

        Self { n_cells, u_hierarchy, v_hierarchy, cfg }
    }

    /// Apply the AMG preconditioner: `z = M⁻¹ r` for CG's inner
    /// iteration. Two independent V-cycles (Option B'), wrapped
    /// with `project_velocity` to preserve the 2-D velocity null
    /// space (Option α — mirrors the Jacobi precond convention).
    pub fn apply(&self, r: &[f64], z: &mut [f64]) {
        debug_assert_eq!(r.len(), 2 * self.n_cells);
        debug_assert_eq!(z.len(), 2 * self.n_cells);

        let n = self.n_cells;
        // Pre-project the residual to remove null-space content
        // (constant-per-component).
        let mut r_proj = r.to_vec();
        {
            let (r_x, r_y) = r_proj.split_at_mut(n);
            super::nullspace::project_velocity(r_x, r_y);
        }

        // V-cycle on each half independently.
        let (z_x, z_y) = z.split_at_mut(n);
        for v in z_x.iter_mut() {
            *v = 0.0;
        }
        for v in z_y.iter_mut() {
            *v = 0.0;
        }
        vcycle::v_cycle(&self.u_hierarchy, &self.cfg, &r_proj[..n], z_x);
        vcycle::v_cycle(&self.v_hierarchy, &self.cfg, &r_proj[n..], z_y);

        // Post-project to clean residual null-space drift.
        super::nullspace::project_velocity(z_x, z_y);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn config_defaults_are_documented() {
        let c = AmgConfig::default();
        assert_eq!(c.strong_connection_threshold, 0.25);
        assert_eq!(c.max_levels, 7);
        assert_eq!(c.min_coarse_unknowns, 50);
        assert_eq!(c.pre_smooth_sweeps, 1);
        assert_eq!(c.post_smooth_sweeps, 1);
    }

    #[test]
    fn preconditioner_build_produces_non_empty_hierarchies() {
        // Post-Phase-2.6 invariant: the setup pipeline builds at
        // least one level per hierarchy, with an LU-factored
        // coarsest level at the end.
        use crate::tectonics_v2::field::Field2D;
        use crate::tectonics_v2::stokes::operator::StokesGrid;
        use crate::tectonics_v2::stokes::sparse_assembly::assemble_picard_csr;
        let nx = 16;
        let ny = 16;
        let grid = StokesGrid::new(nx, ny, 1.0 / nx as f64, 1.0 / ny as f64);
        let eta = Field2D::filled(nx, ny, 1.0);
        let a = assemble_picard_csr(&grid, &eta, None);
        let p = AmgPreconditioner::build(&a, nx * ny, AmgConfig::default());
        assert_eq!(p.n_cells, nx * ny);
        assert!(!p.u_hierarchy.levels.is_empty(), "u hierarchy is empty");
        assert!(!p.v_hierarchy.levels.is_empty(), "v hierarchy is empty");
        // Coarsest level carries LU; no finer does.
        assert!(p.u_hierarchy.levels.last().unwrap().coarse_lu.is_some());
        assert!(p.v_hierarchy.levels.last().unwrap().coarse_lu.is_some());
    }
}
