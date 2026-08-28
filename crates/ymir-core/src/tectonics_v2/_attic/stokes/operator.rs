//! Discrete thin-viscous-sheet momentum operator on a MAC (staggered)
//! grid with periodic BCs.
//!
//! Following England & McKenzie (1982): the depth-integrated
//! horizontal momentum balance reads
//! ```text
//!   -∇·(2 η ε̇(v)) = f_ext
//! ```
//! where `v` is the 2-D horizontal velocity and `f_ext` gathers GPE,
//! plate traction, slab pull, and basal drag. **The 2-D velocity is
//! NOT divergence-free**: `∇·v ≠ 0` is physically meaningful — it is
//! the rate at which the column thickens (`∂_t S + ∇·(Sv) = 0`).
//! There is no incompressibility constraint and no pressure unknown.
//!
//! The discrete operator is
//! ```text
//!   A v ≡ -∇·(2 η ε̇(v))
//! ```
//! expanded per component. For constant η this reduces to
//! `-η (∇² v + ∇(∇·v))`, i.e. Laplacian + grad-div. The grad-div term
//! is real and must be in the discretization — a "pure Laplacian"
//! approximation would drop the physics that couples `v_x` and `v_y`
//! through normal strain.
//!
//! Layout (same as legacy `tectonics/solver/grid.rs`):
//! - `η`, `S` at cell centres `((i+0.5)dx, (j+0.5)dy)`.
//! - `vx` at left vertical face of cell (i, j) — `(i dx, (j+0.5)dy)`.
//! - `vy` at bottom horizontal face of cell (i, j) — `((i+0.5)dx, j dy)`.
//! - `ε̇_xy` and `σ_xy` at nodal corners `(i dx, j dy)`, with η there
//!   computed by **arithmetic 4-point averaging** of the four
//!   surrounding cell centres. (Step 0 documented harmonic averaging;
//!   Step 1 switched to arithmetic so the Newton Jacobian is exactly
//!   symmetric — see the `eta_corner` doc-comment.)

use rayon::prelude::*;

use super::super::rheology::{StrainRate, ViscosityLaw};
use crate::tectonics_v2::field::{Field2D, PeriodicIndex};

/// Geometry needed by the momentum operator. Borrows nothing; η is
/// passed as an argument to the apply/diagonal routines.
pub struct StokesGrid {
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub dy: f64,
    pub idx_x: PeriodicIndex,
    pub idx_y: PeriodicIndex,
}

impl StokesGrid {
    pub fn new(nx: usize, ny: usize, dx: f64, dy: f64) -> Self {
        Self { nx, ny, dx, dy, idx_x: PeriodicIndex::new(nx), idx_y: PeriodicIndex::new(ny) }
    }

    #[inline]
    pub fn n_cells(&self) -> usize {
        self.nx * self.ny
    }
}

/// Arithmetic 4-point average of `η` at a node (corner), built from
/// the four surrounding cell centres.
///
/// # Why not harmonic averaging?
///
/// Step 0 documented "harmonic 4-point averaging" — appropriate for
/// the thin viscous sheet against sharp material contrasts (Gerya
/// 2010 §14.3). Step 1 introduces power-law rheology, which makes
/// `η` a smooth function of `ε̇_II(v)`. The Newton Jacobian of
/// `apply_momentum` w.r.t. `v` then carries a
/// `dη_corner / dη_cell` chain-rule factor; for harmonic averaging
/// this factor is `η_corner² / (4 η_cell²)`, which breaks the
/// natural symmetry of the Newton-extra stress when assembled at
/// mismatched staggered locations.
///
/// With arithmetic averaging, `dη_corner / dη_cell = ¼`, a constant.
/// The Newton-extra stress at the corner becomes
/// `⟨c·contract⟩_cc→co · ε̇_xy(v_k)`, which is the adjoint of the
/// cell-centre contribution and makes the full Jacobian exactly
/// symmetric at discrete level.
///
/// For Step 1's smooth η fields the two averagings are
/// indistinguishable at O(dx²); for sharp viscosity contrasts
/// (cratonic rigidity, Step 9) the choice may be revisited.
#[inline]
fn eta_corner(eta: &Field2D, im: usize, i: usize, jm: usize, j: usize) -> f64 {
    0.25 * (eta.get(im, jm) + eta.get(i, jm) + eta.get(im, j) + eta.get(i, j))
}

/// Apply `A v = -∇·(2 η ε̇(v))  (+ Br · S̃² · v)` on the MAC grid.
///
/// The discretization assembles normal stresses `σ_αα = 2η ∂_α v_α`
/// at cell centres and shear stresses `σ_xy = η (∂_y v_x + ∂_x v_y)`
/// at corners (η harmonic-averaged over the four surrounding cells).
/// The stress divergence is then differenced into face-centred
/// outputs. The sign convention makes `A` SPD on the zero-mean
/// velocity subspace.
///
/// # Basal drag (Step 4)
///
/// If `drag_diag` is `Some(&Br·S̃²)` (a cell-centered field, built
/// by [`crate::tectonics_v2::basal_drag::build_drag_diagonal_field`]),
/// a positive diagonal contribution is added **after** the viscous
/// stencil:
/// ```text
///   out_vx[i, j] += drag_face_x(i, j) · vx[i, j]
///   out_vy[i, j] += drag_face_y(i, j) · vy[i, j]
/// ```
/// where `drag_face_*` is the arithmetic 2-point cell-to-face average
/// of the cell-centered `drag_diag` — the same convention the viscous
/// stencil uses for η at vx/vy faces (see the `d_sigma_xx_dx` block
/// above: `eta_cc_right/left` are cells `(i,j)` and `(im,j)` averaged
/// across the vx face). When `drag_diag` is `None` (basal drag
/// disabled) the augmentation loop is skipped entirely — the fast
/// path matches the Step 0/1/2/3 behaviour bit-for-bit.
pub fn apply_momentum(
    grid: &StokesGrid,
    eta: &Field2D,
    drag_diag: Option<&Field2D>,
    vx: &[f64],
    vy: &[f64],
    out_vx: &mut [f64],
    out_vy: &mut [f64],
) {
    let nx = grid.nx;
    let ny = grid.ny;
    let dx = grid.dx;
    let dy = grid.dy;
    let inv_dx = 1.0 / dx;
    let inv_dy = 1.0 / dy;
    debug_assert_eq!(vx.len(), nx * ny);
    debug_assert_eq!(vy.len(), nx * ny);

    // Step 8.5b: parallelise over rows j via `par_chunks_mut(nx)`. Each
    // cell's output depends only on read-only inputs and the cell
    // index, so execution order is irrelevant to the numeric result
    // (bit-identical across thread counts).
    out_vx.par_chunks_mut(nx).zip(out_vy.par_chunks_mut(nx)).enumerate().for_each(
        |(j, (row_vx, row_vy))| {
            let jp = grid.idx_y.next(j);
            let jm = grid.idx_y.prev(j);
            for i in 0..nx {
                let ip = grid.idx_x.next(i);
                let im = grid.idx_x.prev(i);
                let lin = |ii: usize, jj: usize| jj * nx + ii;

                // ---------- x-momentum at vx(i, j) ----------
                let eta_cc_right = eta.get(i, j);
                let eta_cc_left = eta.get(im, j);
                let dvx_dx_right = (vx[lin(ip, j)] - vx[lin(i, j)]) * inv_dx;
                let dvx_dx_left = (vx[lin(i, j)] - vx[lin(im, j)]) * inv_dx;
                let sigma_xx_right = 2.0 * eta_cc_right * dvx_dx_right;
                let sigma_xx_left = 2.0 * eta_cc_left * dvx_dx_left;
                let d_sigma_xx_dx = (sigma_xx_right - sigma_xx_left) * inv_dx;

                let eta_corner_top = eta_corner(eta, im, i, j, jp);
                let eta_corner_bot = eta_corner(eta, im, i, jm, j);
                let dvx_dy_top = (vx[lin(i, jp)] - vx[lin(i, j)]) * inv_dy;
                let dvx_dy_bot = (vx[lin(i, j)] - vx[lin(i, jm)]) * inv_dy;
                let dvy_dx_top = (vy[lin(i, jp)] - vy[lin(im, jp)]) * inv_dx;
                let dvy_dx_bot = (vy[lin(i, j)] - vy[lin(im, j)]) * inv_dx;
                let sigma_xy_top = eta_corner_top * (dvx_dy_top + dvy_dx_top);
                let sigma_xy_bot = eta_corner_bot * (dvx_dy_bot + dvy_dx_bot);
                let d_sigma_xy_dy = (sigma_xy_top - sigma_xy_bot) * inv_dy;

                row_vx[i] = -(d_sigma_xx_dx + d_sigma_xy_dy);

                // ---------- y-momentum at vy(i, j) ----------
                let eta_cc_top = eta.get(i, j);
                let eta_cc_bot = eta.get(i, jm);
                let dvy_dy_top = (vy[lin(i, jp)] - vy[lin(i, j)]) * inv_dy;
                let dvy_dy_bot = (vy[lin(i, j)] - vy[lin(i, jm)]) * inv_dy;
                let sigma_yy_top = 2.0 * eta_cc_top * dvy_dy_top;
                let sigma_yy_bot = 2.0 * eta_cc_bot * dvy_dy_bot;
                let d_sigma_yy_dy = (sigma_yy_top - sigma_yy_bot) * inv_dy;

                let eta_corner_right = eta_corner(eta, i, ip, jm, j);
                let eta_corner_left = eta_corner(eta, im, i, jm, j);
                let dvx_dy_right = (vx[lin(ip, j)] - vx[lin(ip, jm)]) * inv_dy;
                let dvx_dy_left = (vx[lin(i, j)] - vx[lin(i, jm)]) * inv_dy;
                let dvy_dx_right = (vy[lin(ip, j)] - vy[lin(i, j)]) * inv_dx;
                let dvy_dx_left = (vy[lin(i, j)] - vy[lin(im, j)]) * inv_dx;
                let sigma_xy_right = eta_corner_right * (dvx_dy_right + dvy_dx_right);
                let sigma_xy_left = eta_corner_left * (dvx_dy_left + dvy_dx_left);
                let d_sigma_xy_dx = (sigma_xy_right - sigma_xy_left) * inv_dx;

                row_vy[i] = -(d_sigma_xy_dx + d_sigma_yy_dy);
            }
        },
    );

    // ---------- Basal drag augmentation (Step 4) ----------
    //
    // `drag_diag` is cell-centered `Br · S̃^exp`. Augment by `drag ·
    // v` on both components, with the drag coefficient averaged from
    // cells to faces by the same arithmetic 2-point convention that
    // the viscous stencil uses for η across the vx/vy faces.
    //
    // Done as a separate pass so the viscous code above stays
    // verbatim from Step 0 — this preserves readability and lets
    // `momentum_diagonal` extract the drag-augmented diagonal with
    // the matching arithmetic.
    if let Some(drag) = drag_diag {
        debug_assert_eq!(drag.nx(), nx);
        debug_assert_eq!(drag.ny(), ny);
        out_vx.par_chunks_mut(nx).zip(out_vy.par_chunks_mut(nx)).enumerate().for_each(
            |(j, (row_vx, row_vy))| {
                let jm = grid.idx_y.prev(j);
                for i in 0..nx {
                    let im = grid.idx_x.prev(i);
                    let lin = |ii: usize, jj: usize| jj * nx + ii;
                    let drag_x = 0.5 * (drag.get(im, j) + drag.get(i, j));
                    let drag_y = 0.5 * (drag.get(i, jm) + drag.get(i, j));
                    row_vx[i] += drag_x * vx[lin(i, j)];
                    row_vy[i] += drag_y * vy[lin(i, j)];
                }
            },
        );
    }
}

/// Diagonal of `A` at each velocity DOF, for Jacobi preconditioning.
///
/// For η constant and `drag_diag = None` and `dx = dy = 1` this
/// returns 6 at every DOF — the expected diagonal of the discrete
/// thin-sheet momentum operator (`4 η/dx² + 2 η/dy²` for `vx`,
/// symmetric for `vy`).
///
/// # Case (B) — analytical reconstruction
///
/// This function is an **analytical rewrite** of the viscous stencil
/// in [`apply_momentum`] — `stokes/precond.rs` consumes the diagonal
/// as an external slice and does NOT rebuild it. Therefore any new
/// diagonal contribution (Step 4 basal drag, future Step 7/8 spike
/// operators, …) must be added **explicitly** here with the same
/// cell-to-face averaging that `apply_momentum` uses, or CG's
/// preconditioner drifts silently from the assembled operator.
/// Consistency with `apply_momentum` is enforced by the unit test
/// `v2_precond_drag_diagonal` which probes `A · e_k[k]` against the
/// analytical diagonal at 1e-14.
pub fn momentum_diagonal(
    grid: &StokesGrid,
    eta: &Field2D,
    drag_diag: Option<&Field2D>,
    diag_vx: &mut [f64],
    diag_vy: &mut [f64],
) {
    let nx = grid.nx;
    let ny = grid.ny;
    let dx = grid.dx;
    let dy = grid.dy;
    let inv_dx2 = 1.0 / (dx * dx);
    let inv_dy2 = 1.0 / (dy * dy);

    // Step 8.5b: parallelise over rows — same rationale as
    // `apply_momentum` (cell-local writes, read-only inputs).
    diag_vx.par_chunks_mut(nx).zip(diag_vy.par_chunks_mut(nx)).enumerate().for_each(
        |(j, (row_vx, row_vy))| {
            let jp = grid.idx_y.next(j);
            let jm = grid.idx_y.prev(j);
            for i in 0..nx {
                let ip = grid.idx_x.next(i);
                let im = grid.idx_x.prev(i);

                let eta_right_cc = eta.get(i, j);
                let eta_left_cc = eta.get(im, j);
                let eta_c_top = eta_corner(eta, im, i, j, jp);
                let eta_c_bot = eta_corner(eta, im, i, jm, j);
                row_vx[i] = 2.0 * (eta_right_cc + eta_left_cc) * inv_dx2
                    + (eta_c_top + eta_c_bot) * inv_dy2;

                let eta_top_cc = eta.get(i, j);
                let eta_bot_cc = eta.get(i, jm);
                let eta_c_right = eta_corner(eta, i, ip, jm, j);
                let eta_c_left = eta_corner(eta, im, i, jm, j);
                row_vy[i] = (eta_c_right + eta_c_left) * inv_dx2
                    + 2.0 * (eta_top_cc + eta_bot_cc) * inv_dy2;
            }
        },
    );

    // Basal drag: augment the diagonal with `drag_face_*` (arithmetic
    // 2-point cell-to-face average of the cell-centered `drag_diag`),
    // matching `apply_momentum`'s convention 1-to-1. Disabled path:
    // no loop, no branch beyond the `if let` at top level.
    if let Some(drag) = drag_diag {
        debug_assert_eq!(drag.nx(), nx);
        debug_assert_eq!(drag.ny(), ny);
        diag_vx.par_chunks_mut(nx).zip(diag_vy.par_chunks_mut(nx)).enumerate().for_each(
            |(j, (row_vx, row_vy))| {
                let jm = grid.idx_y.prev(j);
                for i in 0..nx {
                    let im = grid.idx_x.prev(i);
                    let drag_x = 0.5 * (drag.get(im, j) + drag.get(i, j));
                    let drag_y = 0.5 * (drag.get(i, jm) + drag.get(i, j));
                    row_vx[i] += drag_x;
                    row_vy[i] += drag_y;
                }
            },
        );
    }
}

/// Cached pieces of the current Newton iterate used to form the
/// tangent Jacobian. Built once per Newton iteration so the linear
/// solver (CG) reuses it across matrix-vector products.
///
/// # Symmetric-preserving discretisation
///
/// Arithmetic averaging of `η` to corners (see `eta_corner`) gives
/// `dη_corner / dη_cell = ¼`. With this choice the Newton-extra
/// stress can be written using a **single cell-centred scalar**
/// `S(δv) = c · contract` (contract evaluated at cell centres) and
/// an averaging to corners for the shear component:
/// ```text
///   σ^N_xx[cc] = S(δv) · ε̇_xx(v_k)
///   σ^N_yy[cc] = S(δv) · ε̇_yy(v_k)
///   σ^N_xy[co] = ⟨S(δv)⟩_co · ε̇_xy(v_k)
/// ```
/// where `⟨·⟩_co` averages the four cells around a corner.
/// The pairing of cell-centre stress with `ε̇_{xx,yy}(w)|_cc` and
/// corner stress with `ε̇_xy(w)|_co` then collapses (via the
/// adjoint of the averaging operator) to
/// `∑_cc c·contract(u)·contract(w)` — symmetric in `(u, w)`.
///
/// `contract(δv)` at cell centre is
/// ```text
///   (ε̇(v_k):ε̇(δv))|_cc = ε̇_xx(v_k)·ε̇_xx(δv) + ε̇_yy(v_k)·ε̇_yy(δv)
///                      + 2·⟨ε̇_xy(v_k)·ε̇_xy(δv)⟩_cc
/// ```
/// where the shear product is averaged **after** multiplication (not
/// before) — this consistency with the discrete definition of `ε̇_II`
/// is what lets the tangent be the exact Jacobian of the discrete
/// residual.
///
/// The Newton extra is negative semi-definite for shear-thinning
/// (`n > 1`, η' < 0), so the full Jacobian `J = A_picard + A_tangent`
/// is symmetric but **not necessarily SPD** in zones of strong
/// localisation (Gerya §14.4). Symmetry is verified by unit test;
/// positive-definiteness is not required and not tested.
pub struct TangentContext {
    /// η(ε̇_II(v_k)) at cell centres. Also feeds the Picard part via
    /// `apply_momentum` (with the same arithmetic corner averaging).
    pub eta_center: Field2D,
    /// Scalar prefactor `c_cc = η'(ε̇_II_cc) / (ε̇_II_cc + ε̇_min)` at
    /// cell centres.
    pub c_center: Field2D,
    /// Native strain-rate components of `v_k`.
    pub exx_center: Field2D,
    pub eyy_center: Field2D,
    pub exy_corner: Field2D,
}

impl TangentContext {
    /// Build the Newton tangent context from the current strain rate.
    ///
    /// When `cratonic = Some(state)`, two scalings apply consistently
    /// to both `eta_center` and the per-cell tangent
    /// `c_center = dη/dε̇ / (ε̇ + floor)`:
    /// 1. The plastic branch's Bi is replaced by `bi_eff[i,j] =
    ///    Bi · state.bi_multiplier[i,j]` (Step 9 D1 primary mechanism)
    ///    inside both `eta_effective_with_bi_override` and
    ///    `d_eta_effective_d_eps_ii_with_bi_override`. The chain
    ///    rule through the soft-min blend uses the elevated
    ///    `η_p = bi_eff/(2(ε̇+floor))`.
    /// 2. The full result is post-multiplied by
    ///    `state.eta_multiplier[i,j]` (Step 9 D1 secondary mechanism
    ///    — K viscous mult). Both `eta_center` and `c_center` scale
    ///    by the same factor because `m(x)` has no `ε̇` dependence.
    ///
    /// When `cratonic = None`, structural by-pass — bit-identical
    /// to the pre-Step-9 path.
    pub fn from_strain_rate(
        grid: &StokesGrid,
        law: &ViscosityLaw,
        sr: &StrainRate,
        cratonic: Option<&crate::tectonics_v2::cratonic::CratonicState>,
    ) -> Self {
        let nx = grid.nx;
        let ny = grid.ny;
        let mut eta_center = Field2D::new(nx, ny);
        let mut c_center = Field2D::new(nx, ny);
        let global_bi = match law.yielding {
            crate::tectonics_v2::presets::YieldingConfig::Disabled => 0.0,
            crate::tectonics_v2::presets::YieldingConfig::Enabled(ylaw) => ylaw.bi,
        };
        match cratonic {
            None => {
                for j in 0..ny {
                    for i in 0..nx {
                        let eps = sr.eps_ii_center.get(i, j);
                        eta_center.set(i, j, law.eta_effective(eps));
                        c_center.set(
                            i,
                            j,
                            law.d_eta_effective_d_eps_ii(eps) / (eps + law.strain_rate_floor),
                        );
                    }
                }
            }
            Some(state) => {
                for j in 0..ny {
                    for i in 0..nx {
                        let eps = sr.eps_ii_center.get(i, j);
                        let m = state.eta_multiplier.get(i, j);
                        let bi_eff = global_bi * state.bi_multiplier.get(i, j);
                        eta_center.set(i, j, law.eta_effective_with_bi_override(eps, bi_eff) * m);
                        c_center.set(
                            i,
                            j,
                            m * law.d_eta_effective_d_eps_ii_with_bi_override(eps, bi_eff)
                                / (eps + law.strain_rate_floor),
                        );
                    }
                }
            }
        }
        Self {
            eta_center,
            c_center,
            exx_center: sr.exx_center.clone(),
            eyy_center: sr.eyy_center.clone(),
            exy_corner: sr.exy_corner.clone(),
        }
    }
}

/// Apply **only the Newton-extra** part of the Jacobian:
/// `δv ↦ -∇·[c · (ε̇(v_k):ε̇(δv)) · ε̇(v_k)]`.
///
/// For shear-thinning (`n > 1`) this operator is negative
/// semi-definite — consequently the full Jacobian
/// `J = A_picard + apply_tangent` is symmetric but not necessarily
/// positive-definite in zones of strong localisation (Gerya §14.4).
/// Tests that probe structure check **symmetry only**.
pub fn apply_tangent(
    grid: &StokesGrid,
    ctx: &TangentContext,
    dvx: &[f64],
    dvy: &[f64],
    out_vx: &mut [f64],
    out_vy: &mut [f64],
) {
    let nx = grid.nx;
    let ny = grid.ny;
    let dx = grid.dx;
    let dy = grid.dy;
    let inv_dx = 1.0 / dx;
    let inv_dy = 1.0 / dy;
    let idx_x = &grid.idx_x;
    let idx_y = &grid.idx_y;
    let lin = |ii: usize, jj: usize| jj * nx + ii;

    // Step 8.5b: each of the four loops below writes cell-local
    // values and is safely parallelised by row (`par_chunks_mut(nx)`
    // on the Field2D's flat buffer). The loops run sequentially
    // relative to each other because they pipeline intermediate
    // buffers (dexx/deyy/dexy → s_cc → sigma → divergence).

    // --- 1. δv's native strain-rate components ---
    let mut dexx_cc = Field2D::new(nx, ny);
    let mut deyy_cc = Field2D::new(nx, ny);
    let mut dexy_co = Field2D::new(nx, ny);
    dexx_cc
        .data_mut()
        .par_chunks_mut(nx)
        .zip(deyy_cc.data_mut().par_chunks_mut(nx))
        .zip(dexy_co.data_mut().par_chunks_mut(nx))
        .enumerate()
        .for_each(|(j, ((dexx_row, deyy_row), dexy_row))| {
            let jp = idx_y.next(j);
            let jm = idx_y.prev(j);
            for i in 0..nx {
                let ip = idx_x.next(i);
                let im = idx_x.prev(i);
                dexx_row[i] = (dvx[lin(ip, j)] - dvx[lin(i, j)]) * inv_dx;
                deyy_row[i] = (dvy[lin(i, jp)] - dvy[lin(i, j)]) * inv_dy;
                let dvx_dy = (dvx[lin(i, j)] - dvx[lin(i, jm)]) * inv_dy;
                let dvy_dx = (dvy[lin(i, j)] - dvy[lin(im, j)]) * inv_dx;
                dexy_row[i] = 0.5 * (dvx_dy + dvy_dx);
            }
        });

    // --- 2. Cell-centre scalar S(δv) = c_cc · contract_cc(δv) ---
    //     contract_cc(δv) = ε̇_xx(v_k)·ε̇_xx(δv) + ε̇_yy(v_k)·ε̇_yy(δv)
    //                     + 2·⟨ε̇_xy(v_k)·ε̇_xy(δv)⟩_cc
    // Average-of-products on the shear term keeps it consistent with
    // the definition of `ε̇_II_cc`.
    let mut s_cc = Field2D::new(nx, ny);
    s_cc.data_mut().par_chunks_mut(nx).enumerate().for_each(|(j, s_row)| {
        let jp = idx_y.next(j);
        for i in 0..nx {
            let ip = idx_x.next(i);
            let pxy_avg = 0.25
                * (ctx.exy_corner.get(i, j) * dexy_co.get(i, j)
                    + ctx.exy_corner.get(ip, j) * dexy_co.get(ip, j)
                    + ctx.exy_corner.get(i, jp) * dexy_co.get(i, jp)
                    + ctx.exy_corner.get(ip, jp) * dexy_co.get(ip, jp));
            let contract = ctx.exx_center.get(i, j) * dexx_cc.get(i, j)
                + ctx.eyy_center.get(i, j) * deyy_cc.get(i, j)
                + 2.0 * pxy_avg;
            s_row[i] = ctx.c_center.get(i, j) * contract;
        }
    });

    // --- 3. Newton-extra stress components ---
    //   σ^N_xx[cc] = S(δv) · ε̇_xx(v_k)
    //   σ^N_yy[cc] = S(δv) · ε̇_yy(v_k)
    //   σ^N_xy[co] = ⟨S(δv)⟩_co · ε̇_xy(v_k)  (adjoint-consistent
    //                averaging from 4 surrounding cells)
    let mut sigma_xx_cc = Field2D::new(nx, ny);
    let mut sigma_yy_cc = Field2D::new(nx, ny);
    let mut sigma_xy_co = Field2D::new(nx, ny);
    sigma_xx_cc
        .data_mut()
        .par_chunks_mut(nx)
        .zip(sigma_yy_cc.data_mut().par_chunks_mut(nx))
        .zip(sigma_xy_co.data_mut().par_chunks_mut(nx))
        .enumerate()
        .for_each(|(j, ((sxx_row, syy_row), sxy_row))| {
            let jm = idx_y.prev(j);
            for i in 0..nx {
                let im = idx_x.prev(i);
                sxx_row[i] = s_cc.get(i, j) * ctx.exx_center.get(i, j);
                syy_row[i] = s_cc.get(i, j) * ctx.eyy_center.get(i, j);
                let s_avg =
                    0.25 * (s_cc.get(im, jm) + s_cc.get(i, jm) + s_cc.get(im, j) + s_cc.get(i, j));
                sxy_row[i] = s_avg * ctx.exy_corner.get(i, j);
            }
        });

    // --- 4. Divergence: adds to existing out (caller placed Picard there). ---
    out_vx.par_chunks_mut(nx).zip(out_vy.par_chunks_mut(nx)).enumerate().for_each(
        |(j, (row_vx, row_vy))| {
            let jp = idx_y.next(j);
            let jm = idx_y.prev(j);
            for i in 0..nx {
                let ip = idx_x.next(i);
                let im = idx_x.prev(i);
                let d_sigma_xx_dx = (sigma_xx_cc.get(i, j) - sigma_xx_cc.get(im, j)) * inv_dx;
                let d_sigma_xy_dy = (sigma_xy_co.get(i, jp) - sigma_xy_co.get(i, j)) * inv_dy;
                row_vx[i] += -(d_sigma_xx_dx + d_sigma_xy_dy);
                let d_sigma_xy_dx = (sigma_xy_co.get(ip, j) - sigma_xy_co.get(i, j)) * inv_dx;
                let d_sigma_yy_dy = (sigma_yy_cc.get(i, j) - sigma_yy_cc.get(i, jm)) * inv_dy;
                row_vy[i] += -(d_sigma_xy_dx + d_sigma_yy_dy);
            }
        },
    );
}

/// Apply the full Newton Jacobian `J δv = A_picard δv + A_tangent δv`.
///
/// `ctx.eta_center` is used for the Picard part via harmonic averaging
/// at corners, so both pieces share the same η field and the operator
/// collapses to a single symmetric (possibly indefinite) linear map.
///
/// # Basal drag
///
/// Basal drag's Jacobian is `+Br · S̃² · I` (identity-scaled,
/// diagonal), which lives entirely in the Picard block. Forwarding
/// `drag_diag` to [`apply_momentum`] is sufficient; [`apply_tangent`]
/// stays unaware of drag since `S̃` is frozen during a Newton solve.
pub fn apply_jacobian(
    grid: &StokesGrid,
    ctx: &TangentContext,
    drag_diag: Option<&Field2D>,
    dvx: &[f64],
    dvy: &[f64],
    out_vx: &mut [f64],
    out_vy: &mut [f64],
) {
    apply_momentum(grid, &ctx.eta_center, drag_diag, dvx, dvy, out_vx, out_vy);
    apply_tangent(grid, ctx, dvx, dvy, out_vx, out_vy);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dot(a: &[f64], b: &[f64]) -> f64 {
        a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
    }

    #[test]
    fn eta_corner_is_arithmetic_average() {
        let mut eta = Field2D::new(4, 4);
        eta.set(0, 0, 1.0);
        eta.set(1, 0, 2.0);
        eta.set(0, 1, 3.0);
        eta.set(1, 1, 4.0);
        // Corner (1, 1) averages cells (0, 0), (1, 0), (0, 1), (1, 1).
        let v = eta_corner(&eta, 0, 1, 0, 1);
        assert!((v - 2.5).abs() < 1e-14, "arithmetic average at corner = {}", v);
    }

    #[test]
    fn momentum_on_zero_is_zero() {
        let grid = StokesGrid::new(8, 8, 0.125, 0.125);
        let eta = Field2D::filled(8, 8, 1.0);
        let vx = vec![0.0; 64];
        let vy = vec![0.0; 64];
        let mut out_vx = vec![9.9; 64];
        let mut out_vy = vec![9.9; 64];
        apply_momentum(&grid, &eta, None, &vx, &vy, &mut out_vx, &mut out_vy);
        for v in out_vx.iter().chain(out_vy.iter()) {
            assert_eq!(*v, 0.0);
        }
    }

    #[test]
    fn momentum_diagonal_for_constant_eta_is_6() {
        let grid = StokesGrid::new(8, 8, 1.0, 1.0);
        let eta = Field2D::filled(8, 8, 1.0);
        let mut dvx = vec![0.0; 64];
        let mut dvy = vec![0.0; 64];
        momentum_diagonal(&grid, &eta, None, &mut dvx, &mut dvy);
        for (k, (&a, &b)) in dvx.iter().zip(dvy.iter()).enumerate() {
            assert!((a - 6.0).abs() < 1e-12, "diag_vx[{}] = {}", k, a);
            assert!((b - 6.0).abs() < 1e-12, "diag_vy[{}] = {}", k, b);
        }
    }

    #[test]
    fn momentum_is_symmetric() {
        // A symmetric ⇒ ⟨A u, w⟩ = ⟨u, A w⟩ for all u, w. This is the
        // property that justifies CG on the thin-sheet operator.
        let nx = 8;
        let ny = 8;
        let grid = StokesGrid::new(nx, ny, 0.13, 0.17);
        let mut eta = Field2D::new(nx, ny);
        for j in 0..ny {
            for i in 0..nx {
                let e = 1.0 + 0.5 * ((i * 3 + j * 7) % 5) as f64 / 5.0;
                eta.set(i, j, e);
            }
        }
        let n2 = nx * ny;
        let mut ux = vec![0.0; n2];
        let mut uy = vec![0.0; n2];
        let mut wx = vec![0.0; n2];
        let mut wy = vec![0.0; n2];
        for k in 0..n2 {
            ux[k] = ((k as f64 * 1.7).sin()) * 1.1;
            uy[k] = ((k as f64 * 2.3).cos()) * 0.7;
            wx[k] = ((k as f64 * 0.9).sin()) * 0.5;
            wy[k] = ((k as f64 * 1.3).cos()) * 1.3;
        }
        let mut aux_x = vec![0.0; n2];
        let mut aux_y = vec![0.0; n2];
        let mut awx = vec![0.0; n2];
        let mut awy = vec![0.0; n2];
        apply_momentum(&grid, &eta, None, &ux, &uy, &mut aux_x, &mut aux_y);
        apply_momentum(&grid, &eta, None, &wx, &wy, &mut awx, &mut awy);
        let lhs = dot(&aux_x, &wx) + dot(&aux_y, &wy);
        let rhs = dot(&ux, &awx) + dot(&uy, &awy);
        let rel = (lhs - rhs).abs() / lhs.abs().max(rhs.abs()).max(1.0);
        assert!(rel < 1e-12, "symmetry broken: |lhs-rhs|/max = {}", rel);
    }

    /// A non-divergence-free test input must produce a non-zero output
    /// through the grad-div term. This catches the common bug of
    /// discretizing only the Laplacian and dropping the coupling
    /// introduced by `∇(∇·v)`.
    #[test]
    fn momentum_includes_grad_div_coupling() {
        let nx = 8;
        let ny = 8;
        let grid = StokesGrid::new(nx, ny, 1.0, 1.0);
        let eta = Field2D::filled(nx, ny, 1.0);
        // vx = +1 on left half, -1 on right half — cell-wise
        // deliberately divergent flow. The normal-strain contribution
        // drives the nonzero diagonal response.
        let mut vx = vec![0.0; nx * ny];
        let vy = vec![0.0; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                vx[j * nx + i] = if i < nx / 2 { 1.0 } else { -1.0 };
            }
        }
        let mut out_vx = vec![0.0; nx * ny];
        let mut out_vy = vec![0.0; nx * ny];
        apply_momentum(&grid, &eta, None, &vx, &vy, &mut out_vx, &mut out_vy);
        let peak = out_vx.iter().fold(0.0_f64, |a, &v| a.max(v.abs()));
        assert!(peak > 0.1, "grad-div coupling missing: peak={}", peak);
    }
}
