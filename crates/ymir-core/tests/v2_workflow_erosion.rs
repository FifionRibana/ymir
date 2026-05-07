//! Step 12 Phase 2 acceptance — low-res parametric erosion (D2).
//!
//! Three tests pinned to the issue's acceptance criteria:
//!
//! - **#2 / `v2_workflow_erosion_mass_balanced`** — `β = 1.0` conserves
//!   total mass within `1e-6`. Algebraically the residual is
//!   `O(ε · N · Δh̄)` ≈ machine-precision; the `1e-6` threshold is the
//!   issue spec, intentionally generous.
//! - **#3 / `v2_workflow_erosion_diffusive`** — `β = 0.0` makes total
//!   continental mass decrease monotonically over `N` consecutive
//!   cycles.
//! - **`v2_workflow_erosion_applied_everywhere`** — interior cells of
//!   a continental disk (far from the coast) lose mass over 10
//!   cycles. This is the **counter-isostasy contract**: erosion must
//!   reach interior cratons, otherwise isostasy bulges them
//!   indefinitely (the whole motivation for the per-cycle Phase A
//!   erosion mechanism).

use ymir_core::tectonics_v2::field::Field2D;
use ymir_core::tectonics_v2::workflow::{low_res_erosion, PhaseAParams};

/// 32² periodic continental field with a small two-frequency relief.
/// Range `[0.6, 0.8]` — strictly above the test sea-level reference
/// `0.5`, so every cell is continental. Slope non-zero almost
/// everywhere → Phase A erosion engages on every cell.
fn periodic_continental_relief(nx: usize, ny: usize) -> Field2D {
    let mut s = Field2D::new(nx, ny);
    let two_pi = 2.0 * std::f64::consts::PI;
    for j in 0..ny {
        for i in 0..nx {
            let sx = (i as f64 / nx as f64) * two_pi;
            let sy = (j as f64 / ny as f64) * two_pi;
            s.set(i, j, 0.7 + 0.05 * (sx.sin() + sy.cos()));
        }
    }
    s
}

/// 32² continental disk in oceanic surround. Continental cells follow
/// `S̃ = 0.55 + 0.35 · (1 - r/R)²` for `r < R = 12`, i.e. peak `0.9`
/// at the centre, edge ≈ `0.55` (just above the `0.5` continental
/// threshold). Outside the disk: `S̃ = 0.3` (oceanic).
///
/// Used by `v2_workflow_erosion_applied_everywhere` to verify that
/// interior craton cells (far from the coast) do erode — *not* just
/// coastal cells.
fn continental_disk_with_oceanic_surround(nx: usize, ny: usize) -> Field2D {
    let mut s = Field2D::filled(nx, ny, 0.3);
    let cx = (nx as f64 - 1.0) * 0.5;
    let cy = (ny as f64 - 1.0) * 0.5;
    let radius = (nx.min(ny) as f64) * 0.375;
    for j in 0..ny {
        for i in 0..nx {
            let dx = i as f64 - cx;
            let dy = j as f64 - cy;
            let r = (dx * dx + dy * dy).sqrt();
            if r < radius {
                let t = 1.0 - r / radius;
                s.set(i, j, 0.55 + 0.35 * t * t);
            }
        }
    }
    s
}

#[test]
fn v2_workflow_erosion_mass_balanced() {
    let mut s = periodic_continental_relief(32, 32);
    let mass_before: f64 = s.data().iter().sum();

    let params = PhaseAParams { alpha: 0.05, beta: 1.0, ..Default::default() };
    let stats = low_res_erosion::apply(&mut s, &params, 0.5);

    let mass_after: f64 = s.data().iter().sum();
    assert!(
        stats.volume_removed > 0.0,
        "erosion must engage on the relief field"
    );
    let drift = (mass_after - mass_before).abs();
    assert!(
        drift < 1.0e-6,
        "β=1.0 mass conservation residual {drift:.3e} exceeds 1e-6 \
         (before={mass_before:.6}, after={mass_after:.6}); investigate \
         floating-point summation, do NOT relax the threshold"
    );
}

#[test]
fn v2_workflow_erosion_diffusive() {
    let mut s = periodic_continental_relief(32, 32);
    let params = PhaseAParams { alpha: 0.05, beta: 0.0, ..Default::default() };

    let mut prev_mass: f64 = s.data().iter().sum();
    let n_cycles = 5;
    for cycle in 1..=n_cycles {
        let stats = low_res_erosion::apply(&mut s, &params, 0.5);
        assert!(
            stats.volume_removed > 0.0,
            "cycle {cycle}: must keep eroding while continental mass \
             remains (volume_removed = 0 means the relief flattened, \
             unexpected this early)"
        );
        let mass: f64 = s.data().iter().sum();
        assert!(
            mass < prev_mass,
            "cycle {cycle}: continental mass must decrease monotonically \
             with β=0.0 (got {mass} >= prev {prev_mass})"
        );
        prev_mass = mass;
    }
}

#[test]
fn v2_workflow_erosion_applied_everywhere() {
    // Continental disk centred on (16, 16), radius 12 → centre cell
    // is at distance 12 from the coast. We measure `Δmass` on the
    // centre cell + a near-centre cell, both > 5 cells from the coast.
    let nx = 32;
    let ny = 32;
    let mut s = continental_disk_with_oceanic_surround(nx, ny);

    // Sanity: probe cells must be continental at start.
    let probes = [(16, 16), (15, 16), (16, 15), (14, 14)];
    for &(i, j) in &probes {
        assert!(
            s.get(i, j) > 0.5,
            "probe ({i}, {j}) must start continental, got {}",
            s.get(i, j)
        );
    }
    let initial: Vec<f64> = probes.iter().map(|&(i, j)| s.get(i, j)).collect();

    let params = PhaseAParams { alpha: 0.05, beta: 0.0, ..Default::default() };
    for _ in 0..10 {
        low_res_erosion::apply(&mut s, &params, 0.5);
    }

    for (k, &(i, j)) in probes.iter().enumerate() {
        let s_after = s.get(i, j);
        let delta = initial[k] - s_after;
        assert!(
            delta > 0.0,
            "interior cell ({i}, {j}) must erode over 10 cycles \
             (counter-isostasy contract): before={}, after={}, Δ={delta}",
            initial[k],
            s_after
        );
    }
}
