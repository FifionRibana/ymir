//! End-to-end relaxation contract for mantle forcing.
//!
//! Setup: mantle forcing in **isolation** — GPE, yielding, basal
//! drag, boundary sources, slab-pull all Disabled. Uniform `S̃ = 1`,
//! small constant η (chosen so `η · k² ≪ coupling` at the pattern's
//! dominant wave numbers). Solve Stokes once; since η is constant
//! and there is no nonlinearity, the Picard/Newton loop converges
//! in one step — equivalent to "fully relaxed steady state".
//!
//! At steady state with constant η and coupling-augmented diagonal
//! `(−η∇² + coupling) v = coupling · Mf · v_mantle`, the solution
//! for each Fourier mode `k` of `v_mantle` is
//!
//! ```text
//!   v_k = coupling · Mf · v_mantle_k / (η · k² + coupling)
//! ```
//!
//! We measure **magnitude alignment** `α = <v, v_m> / |v_m|²`
//! (inner product with the driver, normalised). This is the
//! scalar fraction of `Mf · v_mantle` that `v` reproduces. Large
//! coupling → α ≈ 1 (strong tracking). Small coupling → α small
//! (loose tracking). The spec thresholds (> 0.95 at c=10,
//! > 0.3 at c=0.3) are met at η = 0.01, with the dominant Fourier
//! wave numbers k ∈ {2π, 4π, 6π} giving `η · k² ∈ [0.4, 3.6]`.

use ymir_core::tectonics_v2::field::{Field2D, PeriodicIndex};
use ymir_core::tectonics_v2::forcing::{BodyForce, MantleForce, SimulationState, VectorField};
use ymir_core::tectonics_v2::mantle::pattern::build_mantle_pattern;
use ymir_core::tectonics_v2::mantle::{MantleConfig, build_mantle_diagonal_field};
use ymir_core::tectonics_v2::stokes::{Grid, SheetConfig, solve_sheet};

struct RelaxationRun {
    alignment: f64,
    peak_v: f64,
    peak_v_mantle_scaled: f64,
}

fn run_relaxation(coupling: f64, mf: f64) -> RelaxationRun {
    let nx = 32;
    let ny = 32;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let grid = Grid::new(nx, ny, dx, dy);

    // Small η so coupling dominates the viscous drag at the
    // pattern's single-mode wave number. The test probes the
    // mantle coupling contract, not the general Stokes
    // response; keeping the pattern to a pure k=1 mode gives
    // a predictable alignment coefficient
    // `coupling / (2·η·k² + coupling)` — uniform across all
    // samples. The `2η·k²` factor (not `η·k²`) comes from the
    // thin-viscous-sheet operator: for a div-free velocity,
    // `-∇·(2η·ε̇) = -η·∇²v` in x AND a coupled term in y,
    // doubling the effective diagonal vs a plain Laplacian.
    // η = 0.005 with k² = (2π)² ≈ 39.5 gives 2η·k² ≈ 0.4,
    // well below the coupling strengths probed below.
    let eta = Field2D::filled(nx, ny, 0.005);

    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);
    let s = Field2D::filled(nx, ny, 1.0);

    // Hand-crafted nodal ψ = sin(2π·x) · sin(2π·y). Single
    // Fourier mode, k=1 in each direction → k² = (2π)² ≈ 39.5
    // is the ONLY spatial frequency in v_mantle.
    use std::f64::consts::TAU;
    let mut psi = Field2D::new(nx, ny);
    for j in 0..ny {
        let y = j as f64 / ny as f64;
        for i in 0..nx {
            let x = i as f64 / nx as f64;
            psi.set(i, j, (TAU * x).sin() * (TAU * y).sin());
        }
    }
    let pattern = build_mantle_pattern(&psi, dx, dy, &idx_x, &idx_y);

    // Mantle diagonal: coupling · S̃. No drag to sum — total_diag
    // is just the mantle contribution.
    let mantle_cfg =
        MantleConfig::Enabled { mf, coupling, num_modes: 1, seed: 42, evolution_rate: 0.0 };
    let mantle_diag = build_mantle_diagonal_field(&mantle_cfg, &s).expect("enabled → Some");

    // Constant-RHS contribution: coupling · S̃ · Mf · v_pattern.
    let mut fx = Field2D::new(nx, ny);
    let mut fy = Field2D::new(nx, ny);
    let state = SimulationState { nx, ny, dx, dy, idx_x: &idx_x, idx_y: &idx_y, s: &s };
    MantleForce::new(mf, coupling, &pattern, &s)
        .accumulate(&state, &mut VectorField { fx: &mut fx, fy: &mut fy });

    let mut vx = vec![0.0; nx * ny];
    let mut vy = vec![0.0; nx * ny];
    let cfg = SheetConfig::default();
    let stats =
        solve_sheet(&grid, &eta, Some(&mantle_diag), fx.data(), fy.data(), &mut vx, &mut vy, &cfg);
    assert!(stats.converged, "CG failed to converge: {:?}", stats);

    // Alignment metric: α = <v, v_m_scaled> / |v_m_scaled|²
    // where v_m_scaled = Mf · v_pattern. This is the scalar
    // fraction of the driver reproduced by v. Direction-only
    // alignment (cosine similarity) is always ~1 in this
    // single-force setup because nothing else pushes v; we need
    // a magnitude-aware metric to distinguish coupling regimes.
    let mut dot = 0.0_f64;
    let mut norm_m_sq = 0.0_f64;
    let mut peak_v = 0.0_f64;
    let mut peak_vm = 0.0_f64;
    for k in 0..nx * ny {
        let vmx = mf * pattern.v_mantle_x.data()[k];
        let vmy = mf * pattern.v_mantle_y.data()[k];
        let vx_k = vx[k];
        let vy_k = vy[k];
        dot += vx_k * vmx + vy_k * vmy;
        norm_m_sq += vmx * vmx + vmy * vmy;
        let v_mag = (vx_k * vx_k + vy_k * vy_k).sqrt();
        let vm_mag = (vmx * vmx + vmy * vmy).sqrt();
        if v_mag > peak_v {
            peak_v = v_mag;
        }
        if vm_mag > peak_vm {
            peak_vm = vm_mag;
        }
    }
    let alignment = if norm_m_sq > 0.0 { dot / norm_m_sq } else { 0.0 };
    RelaxationRun { alignment, peak_v, peak_v_mantle_scaled: peak_vm }
}

#[test]
fn strong_coupling_tracks_mantle_closely() {
    let r = run_relaxation(10.0, 1.0);
    eprintln!(
        "coupling=10: alignment = {:.3}, peak|v| = {:.3e}, peak|Mf·v_m| = {:.3e}",
        r.alignment, r.peak_v, r.peak_v_mantle_scaled,
    );
    assert!(r.alignment > 0.95, "alignment at coupling=10 = {:.3} (spec: > 0.95)", r.alignment,);
}

#[test]
fn weak_coupling_tracks_mantle_loosely() {
    let r = run_relaxation(0.3, 1.0);
    eprintln!(
        "coupling=0.3: alignment = {:.3}, peak|v| = {:.3e}, peak|Mf·v_m| = {:.3e}",
        r.alignment, r.peak_v, r.peak_v_mantle_scaled,
    );
    assert!(r.alignment > 0.3, "alignment at coupling=0.3 = {:.3} (spec: > 0.3)", r.alignment,);
    // Sanity: weak coupling should NOT clamp v to v_mantle —
    // alignment strictly below strong-coupling result.
    assert!(
        r.alignment < 0.95,
        "alignment at coupling=0.3 = {:.3} matches strong-coupling expectation (setup suspicious)",
        r.alignment,
    );
}

/// Monotonicity: alignment increases with coupling.
#[test]
fn alignment_is_monotonic_in_coupling() {
    let mut prev = -1.0_f64;
    for &c in &[0.3_f64, 1.0, 3.0, 10.0] {
        let r = run_relaxation(c, 1.0);
        eprintln!("coupling={}: alignment = {:.3}", c, r.alignment);
        assert!(
            r.alignment > prev - 1e-6,
            "alignment non-monotonic at coupling={} (prev={}, current={})",
            c,
            prev,
            r.alignment,
        );
        prev = r.alignment;
    }
}
