//! Phase R2 — macro mass redistribution.
//!
//! Combines three Phase A mechanisms into one atomic call:
//!
//! - **M1 — erosion** per continental cell: `Δh = α · slope · (S̃ -
//!   sea_level)`. `slope = max |S̃[i] − S̃[neighbour]|` over the 4 NESW
//!   neighbours (same formula as the legacy
//!   [`super::low_res_erosion::apply`] — preserved for continuity).
//! - **M2 — drainage + deposition** via
//!   [`super::drainage::compute_drainage_targets`]. Each cell's eroded
//!   mass is transported to its drainage destination (oceanic cell,
//!   pit, or max-distance terminus).
//! - **M3 — isostatic rebound** as an *implicit* compensation: a
//!   fraction `rebound_ratio` of the eroded mass stays "behind" in
//!   the source cell — not in S̃, but in altitude (the crust rises by
//!   isostatic equilibrium when sediment leaves). Only the
//!   non-rebound fraction `1 - rebound_ratio` is actually moved as
//!   crustal thickness. With the Earth-calibrated default
//!   `rebound_ratio = 0.80` (ρ_crust / ρ_mantle ≈ 2700/3300 ≈ 0.82),
//!   20 % of the eroded volume migrates downstream per cycle and
//!   80 % is rebound-compensated.
//!
//! ## Atomic arithmetic
//!
//! For every cell `i` with `eroded[i] > 0`:
//!
//! ```text
//! net_loss[i]         = eroded[i] · (1 - rebound_ratio)
//! deposited[target_i] = eroded[i] · (1 - rebound_ratio)
//! S̃[i]        -= net_loss[i]
//! S̃[target_i] += deposited[i]
//! ```
//!
//! The implicit rebound `eroded[i] · rebound_ratio` is *not* written to
//! S̃ — it represents the crust that didn't migrate. Total mass is
//! conserved exactly modulo IEEE-754 rounding (one subtraction + one
//! addition per eroded cell ≈ ε per operation; 64² grid → drift
//! ≤ 1e-14 in absolute terms, far below the 1e-10 contract).
//!
//! When `target_i == i` (pit anchor, oceanic cell, or max-distance
//! terminus where the wrap pointed back to start — last is
//! degenerate but well-defined), the subtract / add at the same
//! index cancels arithmetically: net effect on `S̃` is zero. The
//! cell-physical interpretation is "closed continental basin —
//! sediment accumulates locally, isostasy restores altitude through
//! the implicit rebound term".
//!
//! ## Two-pass implementation
//!
//! Pass 1 reads a *consistent snapshot* of `S̃` to compute `eroded[]`
//! and `drainage targets`. Pass 2 mutates `S̃` in-place. This avoids
//! intra-cycle order-dependence — every cell's eroded volume is
//! computed from the *same* `S̃` state, even when its target cell is
//! also a source that gets eroded.

use super::drainage::compute_drainage_targets;
use super::PhaseAParams;
use crate::tectonics_v2::field::Field2D;

/// Per-call diagnostics emitted alongside the in-place mutation.
///
/// `total_eroded` is the gross sum of `α · slope · (S̃ − sea_level)`
/// over all continental cells. `total_deposited` is the fraction that
/// actually migrated as S̃; `total_rebound` is the implicit fraction
/// that compensated locally. `mass_balance_check` is the absolute
/// difference `|sum_after − sum_before|` — should sit at ≤ 1e-10 for
/// any non-pathological run (the R2 conservation acceptance test
/// pins this).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RedistributionStats {
    pub total_eroded: f64,
    pub total_deposited: f64,
    pub total_rebound: f64,
    pub mass_balance_check: f64,
    pub peak_delta_h: f64,
    pub max_path_length: u8,
}

/// Apply one cycle of macro mass redistribution in-place on `s`.
///
/// See module docstring for the algorithm. Returns
/// [`RedistributionStats`] for the orchestrator's per-cycle metrics
/// (R3 wires this into [`super::CycleOutput`]).
///
/// Only `alpha`, `isostatic_rebound_ratio`, and `max_drainage_distance`
/// from `params` are consumed — `n_cycles` and `k_cycle` are
/// orchestrator-level knobs that the caller layers on top.
pub fn apply(
    s: &mut Field2D,
    params: &PhaseAParams,
    sea_level_reference: f64,
) -> RedistributionStats {
    let nx = s.nx();
    let ny = s.ny();
    let n_cells = nx * ny;
    let alpha = params.alpha;
    let factor = 1.0 - params.isostatic_rebound_ratio;

    // Periodic NESW LUTs (same convention as drainage + low_res_erosion).
    let prev_x: Vec<usize> = (0..nx).map(|i| (i + nx - 1) % nx).collect();
    let next_x: Vec<usize> = (0..nx).map(|i| (i + 1) % nx).collect();
    let prev_y: Vec<usize> = (0..ny).map(|j| (j + ny - 1) % ny).collect();
    let next_y: Vec<usize> = (0..ny).map(|j| (j + 1) % ny).collect();

    // Pass 1 — compute eroded[] from a snapshot of `s` (read-only).
    let mut eroded = vec![0.0_f64; n_cells];
    {
        let data = s.data();
        for j in 0..ny {
            for i in 0..nx {
                let lin = j * nx + i;
                let s_i = data[lin];
                if s_i <= sea_level_reference {
                    continue;
                }
                let n_lins = [
                    prev_y[j] * nx + i,
                    j * nx + next_x[i],
                    next_y[j] * nx + i,
                    j * nx + prev_x[i],
                ];
                let mut max_slope = 0.0_f64;
                for &nl in &n_lins {
                    let mag = (s_i - data[nl]).abs();
                    if mag > max_slope {
                        max_slope = mag;
                    }
                }
                eroded[lin] = alpha * max_slope * (s_i - sea_level_reference);
            }
        }
    }

    // Pass 1bis — drainage targets on the same snapshot.
    let drainage =
        compute_drainage_targets(s, sea_level_reference, params.max_drainage_distance);

    // Snapshot mass for the conservation diagnostic.
    let mass_before: f64 = s.data().iter().sum();

    // Pass 2 — apply eroded[] · factor as paired (subtract, add) on the
    // same index pair. Multiple sources sharing a target accumulate.
    let mut total_eroded = 0.0_f64;
    let mut peak_dh = 0.0_f64;
    let mut max_pl = 0_u8;
    {
        let data_mut = s.data_mut();
        for lin in 0..n_cells {
            let e = eroded[lin];
            if e <= 0.0 {
                continue;
            }
            let net_loss = e * factor;
            data_mut[lin] -= net_loss;
            data_mut[drainage.target_idx[lin]] += net_loss;

            total_eroded += e;
            if e > peak_dh {
                peak_dh = e;
            }
            if drainage.path_length[lin] > max_pl {
                max_pl = drainage.path_length[lin];
            }
        }
    }

    let total_deposited = total_eroded * factor;
    let total_rebound = total_eroded * params.isostatic_rebound_ratio;
    let mass_after: f64 = s.data().iter().sum();
    let mass_balance_check = (mass_after - mass_before).abs();

    RedistributionStats {
        total_eroded,
        total_deposited,
        total_rebound,
        mass_balance_check,
        peak_delta_h: peak_dh,
        max_path_length: max_pl,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a field with smooth periodic continental relief, all
    /// continental cells > sea_level; same helper as the legacy
    /// `low_res_erosion` test suite so the two modules are tested on
    /// equivalent inputs (Phase 2 contract).
    fn flat_continental_field(nx: usize, ny: usize) -> Field2D {
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

    /// Mixed field: continental patch in the middle, oceanic surround.
    /// Used for the transport-to-oceanic and craton-preservation tests.
    fn continent_in_ocean(nx: usize, ny: usize, sea_level: f64) -> Field2D {
        let mut s = Field2D::filled(nx, ny, sea_level - 0.2); // oceanic floor
        // Continental peak at the centre, descending linearly to coast
        let cx = nx as i64 / 2;
        let cy = ny as i64 / 2;
        let radius = (nx.min(ny) as i64 / 3).max(2);
        for j in 0..ny {
            for i in 0..nx {
                let dx = (i as i64 - cx).abs();
                let dy = (j as i64 - cy).abs();
                let d = dx.max(dy); // Chebyshev — gives a square continent
                if d <= radius {
                    let t = d as f64 / radius as f64; // 0 at peak, 1 at edge
                    let val = (sea_level + 0.5) - t * 0.4; // peak=1.0 → edge=0.6
                    s.set(i, j, val);
                }
            }
        }
        s
    }

    #[test]
    fn macro_conservation_mass_balanced() {
        // Default rebound (0.80) on smooth continental relief — total
        // mass change should be near IEEE-754 floor.
        let nx = 16;
        let ny = 16;
        let mut s = flat_continental_field(nx, ny);
        let mass_before: f64 = s.data().iter().sum();
        let params = PhaseAParams::default();
        let stats = apply(&mut s, &params, 0.5);
        let mass_after: f64 = s.data().iter().sum();
        let drift = (mass_after - mass_before).abs();
        let tol = 1e-10_f64.max(mass_before * 1e-12);
        assert!(
            drift < tol,
            "mass drift {drift} exceeds tolerance {tol} (before={mass_before}, after={mass_after})"
        );
        // Stats consistency
        assert!(stats.total_eroded > 0.0);
        assert!((stats.total_deposited - stats.total_eroded * (1.0 - params.isostatic_rebound_ratio)).abs() < 1e-12);
        assert!((stats.total_rebound - stats.total_eroded * params.isostatic_rebound_ratio).abs() < 1e-12);
        assert!(stats.mass_balance_check < tol);
    }

    #[test]
    fn macro_conservation_at_full_rebound_is_noop() {
        // rebound = 1.0 → factor = 0 → s_field unchanged.
        let nx = 16;
        let ny = 16;
        let s_init = flat_continental_field(nx, ny);
        let mut s = s_init.clone();
        let params = PhaseAParams { isostatic_rebound_ratio: 1.0, ..Default::default() };
        let stats = apply(&mut s, &params, 0.5);

        // Bit-identical comparison: every cell unchanged.
        for k in 0..nx * ny {
            assert_eq!(
                s.data()[k],
                s_init.data()[k],
                "cell {k} should be unchanged when rebound=1.0"
            );
        }
        assert_eq!(stats.total_deposited, 0.0);
        assert!(stats.mass_balance_check < 1e-14);
    }

    #[test]
    fn macro_conservation_at_zero_rebound_redistributes_fully() {
        // rebound = 0.0 → factor = 1.0 → full eroded mass redistributed,
        // still mass-conserving (deposit balances loss exactly).
        let nx = 16;
        let ny = 16;
        let mut s = flat_continental_field(nx, ny);
        let mass_before: f64 = s.data().iter().sum();
        let params = PhaseAParams { isostatic_rebound_ratio: 0.0, ..Default::default() };
        let stats = apply(&mut s, &params, 0.5);
        let mass_after: f64 = s.data().iter().sum();
        let drift = (mass_after - mass_before).abs();
        let tol = 1e-10_f64.max(mass_before * 1e-12);
        assert!(drift < tol, "drift {drift} > tol {tol} at rebound=0");
        assert!((stats.total_deposited - stats.total_eroded).abs() < 1e-12);
        assert_eq!(stats.total_rebound, 0.0);
    }

    #[test]
    fn macro_continental_to_oceanic_transport() {
        // Continental patch + oceanic surround. After one apply pass,
        // the oceanic cells immediately adjacent to the coast should
        // have GAINED mass (sediment deposition), while the peak
        // should have LOST mass.
        let nx = 16;
        let ny = 16;
        let sea_level = 0.4;
        let mut s = continent_in_ocean(nx, ny, sea_level);
        let s_init = s.clone();

        // Use a reasonable α; default rebound 0.80 keeps the effect
        // small but observable (factor=0.20 of eroded goes to sediment).
        let params = PhaseAParams { alpha: 0.05, ..Default::default() };
        let _stats = apply(&mut s, &params, sea_level);

        // Sum mass over oceanic cells (S̃_init ≤ sea_level) — should
        // strictly increase (sediment received). Sum mass over
        // continental cells — should strictly decrease (mass lost).
        let mut oceanic_before = 0.0_f64;
        let mut oceanic_after = 0.0_f64;
        let mut continental_before = 0.0_f64;
        let mut continental_after = 0.0_f64;
        for k in 0..nx * ny {
            if s_init.data()[k] <= sea_level {
                oceanic_before += s_init.data()[k];
                oceanic_after += s.data()[k];
            } else {
                continental_before += s_init.data()[k];
                continental_after += s.data()[k];
            }
        }
        assert!(
            oceanic_after > oceanic_before,
            "oceanic mass should increase by sediment: {oceanic_before} → {oceanic_after}"
        );
        assert!(
            continental_after < continental_before,
            "continental mass should decrease: {continental_before} → {continental_after}"
        );
        // Total mass conservation
        let drift = ((oceanic_after + continental_after)
            - (oceanic_before + continental_before))
            .abs();
        assert!(drift < 1e-10, "drift {drift}");
    }

    #[test]
    fn macro_rebound_preserves_thickness_over_cycles() {
        // 5 cycles with rebound=0.80 on a continental patch — central
        // peak cells should retain ≥ 50 % of their initial thickness
        // (target for the R4 craton-preservation acceptance is "≥ 50 %
        // continental cells with S̃ > 0.8 retained after 5 cycles").
        let nx = 24;
        let ny = 24;
        let sea_level = 0.4;
        let mut s = continent_in_ocean(nx, ny, sea_level);
        let s_init = s.clone();
        let params = PhaseAParams { alpha: 0.05, ..Default::default() };

        for _ in 0..5 {
            let stats = apply(&mut s, &params, sea_level);
            assert!(stats.total_eroded.is_finite());
            // No runaway: total_eroded per cycle must stay bounded.
            assert!(
                stats.peak_delta_h < 0.5,
                "peak_dh {} exploded in cycle",
                stats.peak_delta_h
            );
        }

        // Count cells originally above S̃ = 0.8 that are still above
        // 0.5 × 0.8 = 0.4. (Light contract — the strong R4 contract
        // gates on the *post-tectonic* state, not the standalone
        // R2 acceptance.)
        let mut kept = 0_usize;
        let mut total = 0_usize;
        for k in 0..nx * ny {
            if s_init.data()[k] >= 0.8 {
                total += 1;
                if s.data()[k] >= 0.4 {
                    kept += 1;
                }
            }
        }
        assert!(total > 0, "expected some high-thickness cells");
        let ratio = kept as f64 / total as f64;
        assert!(
            ratio >= 0.5,
            "rebound should keep ≥ 50 % of high-thickness cells, got {:.1} % ({}/{})",
            100.0 * ratio,
            kept,
            total
        );
    }

    #[test]
    fn macro_no_runaway_over_20_cycles() {
        // 20 cycles on smooth continental relief — no NaN / Inf, mass
        // stays conserved across the long run, peak Δh stays bounded.
        let nx = 16;
        let ny = 16;
        let mut s = flat_continental_field(nx, ny);
        let mass_init: f64 = s.data().iter().sum();
        let params = PhaseAParams { alpha: 0.05, ..Default::default() };

        for cycle in 0..20 {
            let stats = apply(&mut s, &params, 0.5);
            for &v in s.data() {
                assert!(v.is_finite(), "non-finite S̃ at cycle {cycle}: {v}");
            }
            assert!(
                stats.peak_delta_h < 1.0,
                "peak_dh runaway at cycle {cycle}: {}",
                stats.peak_delta_h
            );
        }

        let mass_final: f64 = s.data().iter().sum();
        let drift = (mass_final - mass_init).abs();
        // Cumulative drift over 20 cycles — still tight.
        let tol = 1e-9_f64;
        assert!(drift < tol, "cumulative drift {drift} > {tol} over 20 cycles");
    }

    #[test]
    fn macro_pit_target_is_self_no_net_change() {
        // Sanity check on the pit branch — a single uniform-plateau
        // field (one cell pulled down to make it a continent above a
        // virtual sea_level) where every cell's target is itself
        // (perfectly flat → no slope → eroded=0 → no mutation). With
        // strictly-flat S̃, max_slope = 0 → eroded[i] = 0 → no
        // mutation regardless of pit logic.
        let nx = 8;
        let ny = 8;
        let mut s = Field2D::filled(nx, ny, 0.7);
        let mass_before: f64 = s.data().iter().sum();
        let stats = apply(&mut s, &PhaseAParams::default(), 0.5);
        let mass_after: f64 = s.data().iter().sum();
        assert_eq!(stats.total_eroded, 0.0);
        assert_eq!(mass_before, mass_after);
    }

    /// R2 diagnostic — run 5 cycles on a 32² active-medley-like INIT,
    /// print per-cycle drift / total_eroded / total_rebound / max_path
    /// stats so the reviewer can sanity-check the magnitudes.
    ///
    /// ```bash
    /// cargo test --release -p ymir-core --lib \
    ///   tectonics_v2::workflow::macro_redistribution::tests::macro_active_medley_diagnostic \
    ///   -- --ignored --nocapture
    /// ```
    #[test]
    #[ignore]
    fn macro_active_medley_diagnostic() {
        use crate::tectonics_v2::init::{
            init_s_field, InitContext, InitMode, PlateInitData, ProfileShape,
        };
        use crate::tectonics_v2::voronoi::{generate_voronoi, VoronoiConfig};

        let nx = 32;
        let ny = 32;
        let seed = 42;
        let cfg = VoronoiConfig { num_plates: 8, continental_ratio: 0.3 };
        let plates = generate_voronoi(nx, ny, &cfg, seed);
        let plate_data = PlateInitData {
            plate_id: &plates.plate_id,
            plate_type: &plates.plate_type,
            seed_coords: Some(&plates.seed_coords),
        };
        let init_ctx = InitContext { nx, ny, seed, amplitude: 0.2, plate_data: Some(plate_data) };
        let init_mode = InitMode::RadialProfile {
            continental_value: 1.0,
            oceanic_value: 0.2,
            profile_shape: ProfileShape::Smoothstep,
        };
        let mut s = init_s_field(init_mode, &init_ctx);
        let sea_level = 0.5;
        let params = PhaseAParams::default();

        println!(
            "\n=== R2 macro_redistribution diagnostic — 32² active-medley-like INIT ==="
        );
        println!(
            "params: α={}, rebound={}, max_distance={}",
            params.alpha, params.isostatic_rebound_ratio, params.max_drainage_distance
        );
        let m0: f64 = s.data().iter().sum();
        println!("mass init: {:.6}", m0);

        let mut cumulative_drift = 0.0_f64;
        for cycle in 1..=5 {
            let mb: f64 = s.data().iter().sum();
            let stats = apply(&mut s, &params, sea_level);
            let ma: f64 = s.data().iter().sum();
            let drift = (ma - mb).abs();
            cumulative_drift += drift;
            println!(
                "cycle {}: mass {:.6} → {:.6}  drift={:.3e}  eroded={:.4e}  deposited={:.4e}  rebound={:.4e}  peak_dh={:.4e}  max_path={}",
                cycle,
                mb,
                ma,
                drift,
                stats.total_eroded,
                stats.total_deposited,
                stats.total_rebound,
                stats.peak_delta_h,
                stats.max_path_length
            );
        }
        let mf: f64 = s.data().iter().sum();
        println!("mass final: {:.6} (cumulative drift {:.3e})", mf, cumulative_drift);
    }
}
