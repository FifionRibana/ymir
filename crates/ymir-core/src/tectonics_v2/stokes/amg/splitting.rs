//! C/F splitting via the Classical Ruge-Stüben two-pass algorithm
//! — Step 8.5a Phase 2.2.
//!
//! Partitions the fine-grid DOFs into **C-points** (kept on the
//! coarser level) and **F-points** (interpolated from C-points).
//! Reference: Briggs-Henson-McCormick, "A Multigrid Tutorial",
//! 2nd ed., §8.8 — the original Ruge-Stüben 1987 algorithm.
//!
//! # Determinism — reviewer's top-priority invariant
//!
//! Phase 2.2 MUST pass a 100-run byte-determinism test. If the
//! splitting is not deterministic, the hierarchy is not
//! deterministic, and scalar-parity between AMG runs (D9) is
//! impossible. Tie-breaking: lowest column index wins at every
//! decision point (Pass 1 "which point has most strong
//! dependencies", Pass 2 "which F-F strong pair needs
//! promotion").
//!
//! # Phase 2.2 status — stub
//!
//! Real implementation + 100-run determinism test land in the
//! next commit.

/// Label for a fine-grid DOF after splitting.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CfType {
    Coarse,
    Fine,
    /// Used only transiently during the splitting passes; no
    /// output cell should carry `Undecided`.
    Undecided,
}

/// Compute the C/F labelling for the given strong-connection
/// structure. Returned vector has one entry per fine-grid DOF.
///
/// Phase 2.2 stub.
pub fn classical_rs_splitting(_strong: &[Vec<usize>]) -> Vec<CfType> {
    panic!("classical_rs_splitting — lands in Phase 2.2 (Classical RS two-pass)");
}
