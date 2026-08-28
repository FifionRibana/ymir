//! Per-cell plate-boundary classification.
//!
//! ## What this module does
//!
//! For each cell of the grid, look at its 4-connected neighbours and
//! classify the relative motion across any boundary into one of five
//! categories ([`BoundaryType`]):
//!
//! - [`BoundaryType::Internal`] — all 4 neighbours share the cell's
//!   `plate_id` (not on a plate boundary at all).
//! - [`BoundaryType::Convergent`] — at least one neighbour belongs to
//!   a different plate **and** the relative velocity has a clear
//!   net component pointing into the boundary.
//! - [`BoundaryType::Divergent`] — same but the relative velocity
//!   points away from the boundary.
//! - [`BoundaryType::Transform`] — relative velocity is essentially
//!   tangent to the boundary (cosine threshold below).
//! - [`BoundaryType::Ambiguous`] — multiple differing neighbours give
//!   conflicting verdicts (typical of triple junctions). Per W2 of
//!   Issue #123, these are excluded from the orogenic source-term
//!   application.
//!
//! ## Algorithm
//!
//! For each cell `c` and each 4-neighbour `n`:
//!
//! 1. If `plate_id(n) == plate_id(c)`, skip — no boundary on this edge.
//! 2. Compute relative velocity `v_rel = v_c − v_n` from the per-plate
//!    kinematics.
//! 3. Compute the unit normal `n̂` from `c` toward `n` in the integer
//!    cell grid (one of `(±1, 0)` or `(0, ±1)`).
//! 4. Project: `d = v_rel · n̂`. The cosine of the angle between
//!    relative motion and the boundary normal is `d / |v_rel|`.
//! 5. Verdict for this neighbour:
//!    - `cos > +COSINE_THRESHOLD` → Convergent
//!    - `cos < −COSINE_THRESHOLD` → Divergent
//!    - else → Transform
//!    - `|v_rel| < V_REL_FLOOR` (plates essentially co-moving) →
//!      Transform (the projection is meaningless).
//!
//! Aggregating across all differing neighbours: if all non-Transform
//! verdicts agree, that's the cell verdict; if Convergent and
//! Divergent are both present on different neighbours, the cell is
//! [`BoundaryType::Ambiguous`]. Transform never blocks an agreeing
//! Convergent/Divergent group.
//!
//! ## Upper-plate heuristic (W2 watchpoint, Issue #123)
//!
//! For cells classified as Convergent, the cell is marked "upper
//! plate" if **|v_c| > |v_n|** for at least one of its converging
//! neighbours. Heuristic — in real subduction (Andes, Java) the
//! *slower* continental plate overrides the *faster* oceanic plate
//! because of density contrast. Phase 1.2 ignores plate type and
//! goes with the simpler velocity-magnitude heuristic; Phase 3
//! Lallemand subduction-arc closure will refine using
//! continental/oceanic information.
//!
//! Edge case: when `|v_c| == |v_n|` (symmetric collision), neither
//! cell is marked upper. This produces no asymmetric orogenic wedge
//! on those boundaries. The Phase 1.1 hand-tuned kinematics preset
//! has several symmetric pairs (cardinal vs cardinal, |v| = 0.01);
//! the diagonal-vs-cardinal pairs are the asymmetric ones that
//! will carry the orogenic signal in Stage 5.
//!
//! ## Connectivity choice
//!
//! 4-neighbour only (N, S, E, W). 8-neighbour (with diagonals) is
//! more precise on oblique boundaries but doubles per-cell work
//! and complicates the normal-vector arithmetic. Phase 4
//! performance + UI work can revisit if the orogenic visual at
//! 64²-512² shows staircase artefacts on oblique boundaries.

use crate::tectonics_v2::boundaries::plate_type::{PlateType, PlateTypeField};
use crate::tectonics_v2::field::PeriodicIndex;
use crate::tectonics_v2::voronoi::PlateIdField;

use super::kinematics::PlateKinematics;
use super::state::BoolField;

/// Per-cell classification of plate-boundary geometry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryType {
    /// All 4 neighbours share this cell's plate id.
    Internal,
    /// Relative velocity points into the boundary across at least
    /// one differing neighbour, with no conflicting Divergent
    /// verdict on a different neighbour.
    Convergent,
    /// Relative velocity points away from the boundary; no
    /// conflicting Convergent verdict.
    Divergent,
    /// Relative velocity is essentially tangent to the boundary at
    /// every differing neighbour.
    Transform,
    /// Conflicting Convergent / Divergent verdicts on different
    /// neighbours (triple junctions, complex geometries).
    Ambiguous,
}

/// `BoundaryType` field stored row-major matching the v2
/// `Field2D` / `PlateTypeField` conventions. Custom type because the
/// v2 `Field2D` is `f64`-specialised.
#[derive(Clone, Debug)]
pub struct BoundaryTypeField {
    nx: usize,
    ny: usize,
    data: Vec<BoundaryType>,
}

impl BoundaryTypeField {
    pub fn filled(nx: usize, ny: usize, value: BoundaryType) -> Self {
        Self { nx, ny, data: vec![value; nx * ny] }
    }

    #[inline]
    pub fn nx(&self) -> usize {
        self.nx
    }

    #[inline]
    pub fn ny(&self) -> usize {
        self.ny
    }

    #[inline]
    pub fn get(&self, i: usize, j: usize) -> BoundaryType {
        self.data[j * self.nx + i]
    }

    #[inline]
    pub fn set(&mut self, i: usize, j: usize, value: BoundaryType) {
        self.data[j * self.nx + i] = value;
    }

    pub fn data(&self) -> &[BoundaryType] {
        &self.data
    }
}

/// Result of [`classify_boundaries`]: per-cell type plus the upper-
/// plate mask used by the orogenic source term.
pub struct BoundaryInfo {
    pub boundary_type: BoundaryTypeField,
    pub upper_plate_mask: BoolField,
}

impl BoundaryInfo {
    /// Count cells per boundary type. Returns `[Internal, Convergent,
    /// Divergent, Transform, Ambiguous]` counts. Useful for the
    /// cycle-0 stats remontée required at Stage 2 end.
    pub fn counts(&self) -> [usize; 5] {
        let mut counts = [0_usize; 5];
        for &b in self.boundary_type.data() {
            let idx = match b {
                BoundaryType::Internal => 0,
                BoundaryType::Convergent => 1,
                BoundaryType::Divergent => 2,
                BoundaryType::Transform => 3,
                BoundaryType::Ambiguous => 4,
            };
            counts[idx] += 1;
        }
        counts
    }

    /// Count of cells in the upper-plate mask. Sanity check that the
    /// mask is non-empty when at least one asymmetric convergent pair
    /// exists in the kinematics.
    pub fn upper_plate_count(&self) -> usize {
        self.upper_plate_mask.data().iter().filter(|&&b| b).count()
    }
}

/// Cosine of the angle between `v_rel` and the boundary normal at
/// which a verdict flips between Convergent/Divergent and Transform.
///
/// `0.1` is approximately `cos(84°)`. Below this, the relative
/// motion is treated as essentially perpendicular to the boundary
/// normal — i.e. tangent to the boundary itself — and the cell is
/// classified as Transform.
const COSINE_THRESHOLD: f64 = 0.1;

/// Magnitude floor under which `v_rel` is considered negligible (no
/// boundary verdict possible). Below this, the neighbour edge is
/// classified Transform by default.
///
/// Choice rationale: typical Phase 1.1 preset has plate speeds of
/// `~ 0.01` non-dim. Symmetric pairs produce `|v_rel| = 0.02`. A
/// floor at `1e-12` only triggers for numerically-identical
/// velocities, not for the natural `(0.01, 0)` vs `(−0.01, 0)`
/// case.
const V_REL_FLOOR: f64 = 1e-12;

/// Classify every cell in the grid by 4-neighbour plate-boundary
/// geometry under the given per-plate kinematics.
///
/// See module docstring for the algorithm and the upper-plate
/// heuristic.
pub fn classify_boundaries(plate_id: &PlateIdField, kinematics: &PlateKinematics) -> BoundaryInfo {
    let nx = plate_id.nx();
    let ny = plate_id.ny();
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);

    let mut boundary_type = BoundaryTypeField::filled(nx, ny, BoundaryType::Internal);
    let mut upper_plate_mask = BoolField::filled(nx, ny, false);

    // 4-neighbour offsets and their outward normals from the central
    // cell. Each entry: (neighbour_offset_i, neighbour_offset_j,
    // normal_x, normal_y). The normal points from the central cell
    // toward the neighbour, which is what the dot-product convention
    // expects (positive dot = velocity component pointing toward
    // neighbour = convergent if also into the neighbour's volume).
    let neighbours: [(i32, i32, f64, f64); 4] = [
        (1, 0, 1.0, 0.0),   // East
        (-1, 0, -1.0, 0.0), // West
        (0, 1, 0.0, 1.0),   // North
        (0, -1, 0.0, -1.0), // South
    ];

    for j in 0..ny {
        for i in 0..nx {
            let pid_c = plate_id.get(i, j);
            let (vx_c, vy_c) = kinematics.velocities[pid_c as usize];
            let v_mag_c = (vx_c * vx_c + vy_c * vy_c).sqrt();

            // Per-cell aggregation state.
            let mut has_differing_neighbour = false;
            let mut has_convergent = false;
            let mut has_divergent = false;
            let mut is_upper = false;

            for (di, dj, nx_norm, ny_norm) in neighbours.iter() {
                let ni = if *di > 0 {
                    idx_x.next(i)
                } else if *di < 0 {
                    idx_x.prev(i)
                } else {
                    i
                };
                let nj = if *dj > 0 {
                    idx_y.next(j)
                } else if *dj < 0 {
                    idx_y.prev(j)
                } else {
                    j
                };

                let pid_n = plate_id.get(ni, nj);
                if pid_n == pid_c {
                    continue; // same plate — no boundary on this edge
                }

                has_differing_neighbour = true;

                let (vx_n, vy_n) = kinematics.velocities[pid_n as usize];
                let v_mag_n = (vx_n * vx_n + vy_n * vy_n).sqrt();
                let vrel_x = vx_c - vx_n;
                let vrel_y = vy_c - vy_n;
                let vrel_mag = (vrel_x * vrel_x + vrel_y * vrel_y).sqrt();

                if vrel_mag < V_REL_FLOOR {
                    // Plates co-move — no normal motion. Treat as
                    // Transform contribution.
                    continue;
                }

                let dot = vrel_x * nx_norm + vrel_y * ny_norm;
                let cos = dot / vrel_mag;

                if cos > COSINE_THRESHOLD {
                    has_convergent = true;
                    // Upper-plate heuristic — the strictly-faster
                    // plate overrides on this convergent edge.
                    if v_mag_c > v_mag_n {
                        is_upper = true;
                    }
                } else if cos < -COSINE_THRESHOLD {
                    has_divergent = true;
                }
                // else: Transform contribution — no flag flipped.
            }

            let verdict = if !has_differing_neighbour {
                BoundaryType::Internal
            } else if has_convergent && has_divergent {
                BoundaryType::Ambiguous
            } else if has_convergent {
                BoundaryType::Convergent
            } else if has_divergent {
                BoundaryType::Divergent
            } else {
                // Has differing neighbours but no normal-motion
                // verdict on any of them: purely tangential.
                BoundaryType::Transform
            };

            boundary_type.set(i, j, verdict);
            if matches!(verdict, BoundaryType::Convergent) && is_upper {
                upper_plate_mask.set(i, j, true);
            }
        }
    }

    BoundaryInfo { boundary_type, upper_plate_mask }
}

/// #155 maillon 1a — retarget the Davis-Suppe "upper plate" to the
/// CONTINENTAL plate at **O-C (subduction)** convergences.
///
/// [`classify_boundaries`] marks `upper_plate_mask` by the velocity
/// heuristic `is_upper = v_mag_c > v_mag_n` (the FASTER plate, plate_type-
/// blind). Real orogeny thickens the overriding **continental** plate
/// (Andes on South America). This post-process overrides ONLY O-C
/// convergent cells:
/// - a Convergent **Continental** cell with an Oceanic differing-plate
///   neighbour → `upper = true` (overriding plate, gets the DS wedge);
/// - a Convergent **Oceanic** cell with a Continental differing-plate
///   neighbour → `upper = false` (subducting / lower plate).
///
/// **Strict fallback (critical):** C-C (collision) and O-O cells are
/// **NOT touched** — the velocity-based `upper_plate_mask` from
/// `classify_boundaries` stays, so any seed WITHOUT O-C subduction is
/// byte-identical. Reuses the same 4-neighbour periodic convention as
/// `classify_boundaries`.
///
/// Scope is QUI is thickened, NOT the wedge geometry (dome vs chain) nor
/// C-C collision orogeny — those are maillon 1b. Called only on the
/// production DS path (`time_loop`), where `plate_type` is in hand; the
/// 28 test/age-init callers of `classify_boundaries` are unaffected.
pub fn retarget_upper_plate_continental(
    boundary: &mut BoundaryInfo,
    plate_id: &PlateIdField,
    plate_type: &PlateTypeField,
) {
    let nx = plate_id.nx();
    let ny = plate_id.ny();
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);
    for j in 0..ny {
        for i in 0..nx {
            if !matches!(boundary.boundary_type.get(i, j), BoundaryType::Convergent) {
                continue;
            }
            let pid_c = plate_id.get(i, j);
            let cont_c = matches!(plate_type.get(i, j), PlateType::Continental);
            // Differing-plate 4-neighbours = the boundary partners.
            let mut has_oceanic_diff = false;
            let mut has_continental_diff = false;
            for (ni, nj) in
                [(idx_x.next(i), j), (idx_x.prev(i), j), (i, idx_y.next(j)), (i, idx_y.prev(j))]
            {
                if plate_id.get(ni, nj) == pid_c {
                    continue;
                }
                match plate_type.get(ni, nj) {
                    PlateType::Oceanic => has_oceanic_diff = true,
                    PlateType::Continental => has_continental_diff = true,
                }
            }
            if cont_c && has_oceanic_diff {
                // O-C: continental side overrides → gets the wedge.
                boundary.upper_plate_mask.set(i, j, true);
            } else if !cont_c && has_continental_diff {
                // O-C: oceanic side subducts → lower plate, no wedge.
                boundary.upper_plate_mask.set(i, j, false);
            }
            // else C-C / O-O → strictly untouched (velocity fallback).
        }
    }
}

/// #155 maillon 1b-i — companion to [`retarget_upper_plate_continental`]:
/// the mask of cells that ARE O-C continental-override seeds (the cells
/// `retarget` sets `upper = true` via the O-C branch). Used to route the
/// Davis-Suppe wedge GEOMETRY by convergence type: cells whose nearest
/// upper-plate seed is O-C get the margin-peaked ridge profile (Andes),
/// C-C / velocity-fallback seeds keep the rising-to-plateau dome (Tibet).
///
/// Same O-C criterion and 4-neighbour periodic convention as `retarget`:
/// a Convergent **Continental** cell with an Oceanic differing-plate
/// neighbour. C-C / O-O / non-Convergent cells → `false`.
pub fn oc_override_seed_mask(
    boundary: &BoundaryInfo,
    plate_id: &PlateIdField,
    plate_type: &PlateTypeField,
) -> BoolField {
    let nx = plate_id.nx();
    let ny = plate_id.ny();
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);
    let mut mask = BoolField::filled(nx, ny, false);
    for j in 0..ny {
        for i in 0..nx {
            if !matches!(boundary.boundary_type.get(i, j), BoundaryType::Convergent) {
                continue;
            }
            if !matches!(plate_type.get(i, j), PlateType::Continental) {
                continue;
            }
            let pid_c = plate_id.get(i, j);
            let has_oceanic_diff =
                [(idx_x.next(i), j), (idx_x.prev(i), j), (i, idx_y.next(j)), (i, idx_y.prev(j))]
                    .into_iter()
                    .any(|(ni, nj)| {
                        plate_id.get(ni, nj) != pid_c
                            && matches!(plate_type.get(ni, nj), PlateType::Oceanic)
                    });
            if has_oceanic_diff {
                mask.set(i, j, true);
            }
        }
    }
    mask
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a `PlateIdField` from a closure that maps `(i, j)` →
    /// plate id. Helper for tests.
    fn build_plate_id<F: Fn(usize, usize) -> u16>(nx: usize, ny: usize, f: F) -> PlateIdField {
        let mut p = PlateIdField::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                p.set(i, j, f(i, j));
            }
        }
        p
    }

    /// #155 maillon 1a — the BOUNDED proof. The retarget helper must
    /// flip ONLY O-C (subduction) convergent cells and leave C-C
    /// (collision) and O-O cells byte-identical to the velocity-based
    /// `upper_plate_mask`. Geography-independent (no kinematics, no
    /// seed): three synthetic 4×1 boundaries, one per plate-type pair.
    ///
    /// This is the test that proves the C-C/O-O fallback leaks nothing —
    /// the real-seed runs all carry O-C boundaries (seed 2026 retargets
    /// 13044 cells, NOT zero, falsifying the earlier "pure-intraplate"
    /// premise), so the boundedness must be proven structurally here.
    #[test]
    fn retarget_bounded_to_oc_subduction() {
        use crate::tectonics_c1::state::BoolField;
        // 4×1 (ny=1 ⇒ periodic y-neighbour is self ⇒ x-axis only).
        // plate_id [0,0,1,1]; the 0|1 boundary is between cell 1 & 2.
        let nx = 4;
        let ny = 1;
        let plate_id = build_plate_id(nx, ny, |i, _| if i < 2 { 0 } else { 1 });

        // Common: only cells 1 & 2 are Convergent (the boundary pair).
        let make_boundary = |m1: bool, m2: bool| {
            let mut bt = BoundaryTypeField::filled(nx, ny, BoundaryType::Internal);
            bt.set(1, 0, BoundaryType::Convergent);
            bt.set(2, 0, BoundaryType::Convergent);
            let mut mask = BoolField::filled(nx, ny, false);
            mask.set(1, 0, m1);
            mask.set(2, 0, m2);
            BoundaryInfo { boundary_type: bt, upper_plate_mask: mask }
        };
        let pt = |a: PlateType, b: PlateType| {
            // cells 0,1 = `a`; cells 2,3 = `b`.
            let mut p = PlateTypeField::filled(nx, ny, a);
            p.set(2, 0, b);
            p.set(3, 0, b);
            p
        };

        // --- O-C: continental (cells 0,1) overriding oceanic (2,3). ---
        // Seed the mask "wrong" (cont=false, ocean=true) to see the flip.
        let mut oc = make_boundary(false, true);
        retarget_upper_plate_continental(
            &mut oc,
            &plate_id,
            &pt(PlateType::Continental, PlateType::Oceanic),
        );
        assert!(oc.upper_plate_mask.get(1, 0), "O-C: continental cell must become upper=true");
        assert!(!oc.upper_plate_mask.get(2, 0), "O-C: oceanic cell must become lower=false");

        // --- C-C: both continental → strictly untouched. ---
        let mut cc = make_boundary(false, true);
        retarget_upper_plate_continental(
            &mut cc,
            &plate_id,
            &pt(PlateType::Continental, PlateType::Continental),
        );
        assert!(
            !cc.upper_plate_mask.get(1, 0),
            "C-C: cell 1 mask must stay false (velocity fallback)"
        );
        assert!(
            cc.upper_plate_mask.get(2, 0),
            "C-C: cell 2 mask must stay true (velocity fallback)"
        );

        // --- O-O: both oceanic → strictly untouched. ---
        let mut oo = make_boundary(true, false);
        retarget_upper_plate_continental(
            &mut oo,
            &plate_id,
            &pt(PlateType::Oceanic, PlateType::Oceanic),
        );
        assert!(
            oo.upper_plate_mask.get(1, 0),
            "O-O: cell 1 mask must stay true (velocity fallback)"
        );
        assert!(
            !oo.upper_plate_mask.get(2, 0),
            "O-O: cell 2 mask must stay false (velocity fallback)"
        );
    }

    #[test]
    fn convergent_boundary_classified() {
        // Two plates split at column nx/2. Plate 0 on the left
        // moves east; plate 1 on the right moves west. They
        // converge along the boundary at column nx/2.
        let nx = 16;
        let ny = 8;
        let plate_id = build_plate_id(nx, ny, |i, _j| if i < nx / 2 { 0 } else { 1 });
        // Plate 0 faster so it's the upper plate on the boundary.
        let kinematics = PlateKinematics { velocities: vec![(0.02, 0.0), (-0.01, 0.0)] };

        let info = classify_boundaries(&plate_id, &kinematics);

        // Cells at columns nx/2 - 1 (last of plate 0) and nx/2
        // (first of plate 1) sit on the boundary; everything else
        // should be Internal.
        let left_edge_col = nx / 2 - 1;
        let right_edge_col = nx / 2;
        for j in 0..ny {
            assert_eq!(
                info.boundary_type.get(left_edge_col, j),
                BoundaryType::Convergent,
                "left-edge cell ({left_edge_col}, {j}) should be Convergent"
            );
            assert_eq!(
                info.boundary_type.get(right_edge_col, j),
                BoundaryType::Convergent,
                "right-edge cell ({right_edge_col}, {j}) should be Convergent"
            );
            // Faster plate (0) cells on the boundary are upper.
            assert!(
                info.upper_plate_mask.get(left_edge_col, j),
                "plate 0 cell ({left_edge_col}, {j}) should be upper-plate"
            );
            // Slower plate (1) cells on the boundary are NOT upper.
            assert!(
                !info.upper_plate_mask.get(right_edge_col, j),
                "plate 1 cell ({right_edge_col}, {j}) should NOT be upper-plate"
            );
        }

        // Periodic wrap-around (col nx-1 ↔ col 0) inverts the
        // boundary normal, so that line classifies as Divergent
        // in this scenario. We only assert on the "interior"
        // boundary (cols 7/8) which is the test's intent. Global
        // counts should show non-zero Internal + Convergent +
        // Divergent (interior + wrap).
        let counts = info.counts();
        assert!(counts[0] > 0, "Internal count should be > 0 (interior of plates)");
        assert!(
            counts[1] >= 2 * ny,
            "Convergent count should cover both sides of interior boundary"
        );
        assert!(counts[2] >= 2 * ny, "Divergent count should cover the wrap-around boundary");
    }

    #[test]
    fn divergent_boundary_classified() {
        // Mirror of the convergent test: plate 0 moves west, plate 1
        // moves east. They diverge along the column nx/2 boundary.
        let nx = 16;
        let ny = 8;
        let plate_id = build_plate_id(nx, ny, |i, _j| if i < nx / 2 { 0 } else { 1 });
        let kinematics = PlateKinematics { velocities: vec![(-0.02, 0.0), (0.01, 0.0)] };

        let info = classify_boundaries(&plate_id, &kinematics);

        let left_edge_col = nx / 2 - 1;
        let right_edge_col = nx / 2;
        for j in 0..ny {
            assert_eq!(info.boundary_type.get(left_edge_col, j), BoundaryType::Divergent,);
            assert_eq!(info.boundary_type.get(right_edge_col, j), BoundaryType::Divergent,);
            // No cell on the designed divergent boundary should be
            // upper-plate (upper-plate only set on Convergent). The
            // wrap-around boundary at cols 15/0 may carry an
            // upper-plate flag on its Convergent face, so we don't
            // assert globally — only on these two columns.
            assert!(
                !info.upper_plate_mask.get(left_edge_col, j),
                "designed Divergent boundary cell ({left_edge_col}, {j}) should not be upper-plate"
            );
            assert!(
                !info.upper_plate_mask.get(right_edge_col, j),
                "designed Divergent boundary cell ({right_edge_col}, {j}) should not be upper-plate"
            );
        }
    }

    #[test]
    fn transform_boundary_classified() {
        // Two plates split at column nx/2. Plate 0 moves north,
        // plate 1 moves south. Relative motion is along the
        // boundary, not across it → Transform.
        let nx = 16;
        let ny = 8;
        let plate_id = build_plate_id(nx, ny, |i, _j| if i < nx / 2 { 0 } else { 1 });
        let kinematics = PlateKinematics { velocities: vec![(0.0, 0.01), (0.0, -0.01)] };

        let info = classify_boundaries(&plate_id, &kinematics);

        let left_edge_col = nx / 2 - 1;
        let right_edge_col = nx / 2;
        for j in 0..ny {
            assert_eq!(info.boundary_type.get(left_edge_col, j), BoundaryType::Transform);
            assert_eq!(info.boundary_type.get(right_edge_col, j), BoundaryType::Transform);
        }
        // Transform should not produce upper-plate cells.
        assert_eq!(info.upper_plate_count(), 0);
    }

    /// Cycle-0 stats on the actual Phase 1.1 init state. Calibration
    /// sanity check before the Stage 4 source term is wired up:
    /// roughly what fraction of cells will the orogenic source term
    /// affect? Runs `cargo test -- --nocapture` to see the numbers.
    #[test]
    fn phase_1_1_init_state_cycle_0_counts() {
        use crate::tectonics_c1::init::init_c1_state_phase_1_1;
        use crate::tectonics_c1::kinematics::PlateKinematics;

        let state = init_c1_state_phase_1_1(64, 42);
        let kinematics = PlateKinematics::preset_phase_1_1(state.num_plates);
        let info = classify_boundaries(&state.plate_id, &kinematics);

        let counts = info.counts();
        let total = (state.nx() * state.ny()) as f64;
        let upper = info.upper_plate_count();
        let pct = |n: usize| 100.0 * n as f64 / total;

        eprintln!("Phase 1.1 cycle-0 boundary classification on 64² init (seed 42):");
        eprintln!("  Internal    = {:>4}  ({:>5.1} %)", counts[0], pct(counts[0]));
        eprintln!("  Convergent  = {:>4}  ({:>5.1} %)", counts[1], pct(counts[1]));
        eprintln!("  Divergent   = {:>4}  ({:>5.1} %)", counts[2], pct(counts[2]));
        eprintln!("  Transform   = {:>4}  ({:>5.1} %)", counts[3], pct(counts[3]));
        eprintln!("  Ambiguous   = {:>4}  ({:>5.1} %)", counts[4], pct(counts[4]));
        eprintln!(
            "  Upper-plate = {:>4}  ({:>5.1} %)  — orogenic source-term coverage",
            upper,
            pct(upper)
        );
        eprintln!("  Total cells = {} (= 64²)", state.nx() * state.ny());

        // Sanity: counts sum to total grid size, Internal dominates,
        // at least some Convergent / Divergent / Ambiguous detected.
        assert_eq!(counts.iter().sum::<usize>(), state.nx() * state.ny());
        assert!(
            counts[0] > counts[1] + counts[2] + counts[3] + counts[4],
            "Internal should dominate over boundary cells on the v2 default 8-plate layout"
        );
        assert!(counts[1] > 0, "expected some Convergent cells");
        assert!(counts[2] > 0, "expected some Divergent cells");
    }

    #[test]
    fn triple_junction_ambiguous() {
        // Three-plate layout designed so that ONE cell has a
        // Convergent verdict on one neighbour and a Divergent
        // verdict on another.
        //
        // Layout:
        //   Plate 0 — left half (i < 4)
        //   Plate 1 — right half lower (i >= 4 and j < 4)
        //   Plate 2 — right half upper (i >= 4 and j >= 4)
        //
        // Velocities chosen so that:
        //   Plate 0 vs Plate 1 — Convergent (E/W across i=4)
        //   Plate 1 vs Plate 2 — Divergent  (N/S across j=4)
        //
        // The cell at (4, 3) (plate 1, just below the j=4 line)
        // has west neighbour plate 0 (Convergent) and north
        // neighbour plate 2 (Divergent) → Ambiguous.
        let nx = 9;
        let ny = 9;
        let plate_id = build_plate_id(nx, ny, |i, j| {
            if i < 4 {
                0
            } else if j < 4 {
                1
            } else {
                2
            }
        });
        // Plate 0 east (converges with plate 1 which moves west).
        // Plate 1 south (away from plate 2 which moves north) →
        // divergent across the j=4 line. The faster-than-zero y
        // components for plates 1 and 2 give a clear cosine = 1
        // on the divergent edge.
        let kinematics =
            PlateKinematics { velocities: vec![(0.02, 0.0), (-0.01, -0.02), (0.0, 0.02)] };

        let info = classify_boundaries(&plate_id, &kinematics);

        let counts = info.counts();
        assert!(
            counts[4] > 0,
            "expected at least one Ambiguous cell at the triple junction; counts: {:?}",
            counts
        );

        // The designed Ambiguous cell is plate 1's (4, 3) and its
        // mirrors. Spot-check: cell (4, 3) should be Ambiguous.
        assert_eq!(
            info.boundary_type.get(4, 3),
            BoundaryType::Ambiguous,
            "cell (4, 3) should see convergent west + divergent north"
        );
    }
}
