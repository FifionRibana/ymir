//! Bisection calibration of `k_spread` (issue #89 D4).
//!
//! `k_spread` is not a user-facing knob — it is a **closure property
//! of the chosen layout** such that `s_oceanic_mean` at steady state
//! lands in a target band (canonical target `[0.18, 0.22]`, matching
//! `S̃_oceanic ≈ 0.2 ± 10%` from `solver-scaling.md` §4.7).
//!
//! The calibration is a bisection on `k_spread` over a user-supplied
//! bracket (default `[0.1, 1.0]`). At each probe, it runs a single
//! short simulation via a user-supplied closure and reads back
//! `s_oceanic_mean`. The search terminates as soon as
//! `s_oceanic_mean` lands inside `target_range`, or `max_iters` is
//! reached (→ `MaxItersReached`), or the target is
//! unreachable within the bracket (→ `OutOfRange`).
//!
//! The closure signature keeps this module independent of the full
//! harness: anyone with a way to turn `k_spread` into an observed
//! `s_oceanic_mean` can call it. The Step-5 CLI wires in the
//! `BaselineConfig` runner; the `v2_k_spread_calibration` test
//! substitutes a deterministic analytic response so the harness
//! robustness is testable without running real simulations.

/// Canonical `k_spread` bracket to probe.
///
/// Narrowed to `[0.05, 1.0]` after empirical observation on the
/// Step 5 baseline: with type-aware initial `S̃` (oceanic cells at
/// 0.2, continental at 1.0) and 6·τ* simulation time, `s_oceanic_mean`
/// crosses the target band `[0.18, 0.22]` below `k_spread ≈ 0.1`,
/// and sits at ≈ 0.23 at `k_spread = 0.1` itself. The bracket must
/// start below 0.1 to bracket the target from below, but a full-zero
/// lower bound would let the calibration pick `k_spread = 0` and
/// inactivate the rift layer entirely (dropping
/// `boundary_type_diversity` from 2 to 1). A 0.05 floor keeps the
/// spreading mechanism physically active and the resulting
/// `s_oceanic_mean` within ~5% of the target band's upper bound.
///
/// This is the same family of observations as Step 3's
/// `yielding_cell_fraction = 0` and Step 4's `drag/visc ≈ 10⁻⁷`:
/// a quantitative consequence of the honest `Ar = 0.1` thin-sheet
/// scaling. The GPE-only regime at Step 5 has `|Δṽ_conv|` so small
/// (`peak|v| ≈ 5e-5`, `Q_sub ≈ k_sub · 5e-5` per step) that any
/// sizeable `k_spread` sur-fills the oceanic strip. The calibrated
/// value is **expected to evolve** once Steps 7 (slab pull) and 8
/// (mantle forcing) amplify the convergent motion — the `k_spread`
/// of today is not the `k_spread` of tomorrow; it is an evolving
/// closure property of the active-mechanism set, not a fixed
/// constant. Recalibration post-Step 7 / post-Step 8 is
/// anticipated.
pub const K_SPREAD_BRACKET: (f64, f64) = (0.05, 1.0);

/// Input parameters for [`calibrate_k_spread`].
#[derive(Clone, Copy, Debug)]
pub struct KSpreadCalibration {
    pub bracket: (f64, f64),
    pub target_range: (f64, f64),
    pub convergence_tol: f64,
    pub simulation_time: f64,
    pub max_iters: usize,
}

impl KSpreadCalibration {
    /// Canonical Step-5 configuration: bracket `[0.1, 1.0]`, target
    /// `[0.18, 0.22]`, tol `0.01`, simulated time `3·τ*`, up to 20
    /// bisection iterations.
    pub fn step5_default() -> Self {
        Self {
            bracket: K_SPREAD_BRACKET,
            target_range: (0.18, 0.22),
            convergence_tol: 0.01,
            simulation_time: 3.0,
            max_iters: 20,
        }
    }
}

/// Single-iteration record — what the bisection tried and what came
/// back. Reported verbatim in the physics report so the calibration
/// is reproducible.
#[derive(Clone, Copy, Debug)]
pub struct CalibrationIter {
    pub k_spread: f64,
    pub s_oceanic_mean: f64,
}

#[derive(Clone, Debug)]
pub struct CalibrationResult {
    pub k_spread: f64,
    pub iterations: Vec<CalibrationIter>,
    pub final_s_oceanic_mean: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum CalibrationError {
    /// Hit `max_iters` without landing inside `target_range`.
    MaxItersReached { last_k: f64, last_s_oceanic: f64 },
    /// Evaluated the bracket endpoints and found that the target is
    /// not straddled — no monotone bisection can reach it. `bound`
    /// is the endpoint (lo or hi) that was evaluated last; its
    /// `s_oceanic_at_bound` is the best approximation available.
    OutOfRange { bound: f64, s_oceanic_at_bound: f64 },
}

/// Bisect `k_spread` over `cfg.bracket` until `probe(k)` lands
/// inside `cfg.target_range`.
///
/// Assumes `probe(k)` is **monotonically increasing in k** on the
/// bracket (larger spreading → thicker oceanic steady state). That's
/// the physical expectation at Step 5; if the closure violates it,
/// the bisection may converge to a wrong `k`. The assumption is
/// cheap to validate at the endpoints (and the `OutOfRange` detection
/// does that implicitly: if `s_lo` and `s_hi` do not straddle the
/// target, we abort).
pub fn calibrate_k_spread<F>(
    cfg: &KSpreadCalibration,
    mut probe: F,
) -> Result<CalibrationResult, CalibrationError>
where
    F: FnMut(f64) -> f64,
{
    let (mut lo, mut hi) = cfg.bracket;
    let (t_lo, t_hi) = cfg.target_range;
    let mut iterations: Vec<CalibrationIter> = Vec::new();

    // Pre-bracket check: evaluate both endpoints. Abort
    // `OutOfRange` if the target band is not bracketed.
    let s_lo = probe(lo);
    iterations.push(CalibrationIter { k_spread: lo, s_oceanic_mean: s_lo });
    if s_lo >= t_lo && s_lo <= t_hi {
        return Ok(CalibrationResult { k_spread: lo, iterations, final_s_oceanic_mean: s_lo });
    }
    let s_hi = probe(hi);
    iterations.push(CalibrationIter { k_spread: hi, s_oceanic_mean: s_hi });
    if s_hi >= t_lo && s_hi <= t_hi {
        return Ok(CalibrationResult { k_spread: hi, iterations, final_s_oceanic_mean: s_hi });
    }
    // Target-straddling requires the two endpoints to sit on
    // opposite sides of the target band. The spec's canonical case
    // is `s_lo < t_lo` (too little spreading) and `s_hi > t_hi`
    // (too much spreading). Deviation from that → OutOfRange.
    let below_lo = s_lo < t_lo && s_hi < t_lo;
    let above_hi = s_lo > t_hi && s_hi > t_hi;
    if below_lo {
        return Err(CalibrationError::OutOfRange { bound: hi, s_oceanic_at_bound: s_hi });
    }
    if above_hi {
        return Err(CalibrationError::OutOfRange { bound: lo, s_oceanic_at_bound: s_lo });
    }
    // We might also have the pathological ordering `s_lo > t_hi`
    // and `s_hi < t_lo` (anti-monotone). Treat it the same as
    // OutOfRange — the monotonicity assumption is violated, so the
    // bisection can't run.
    if s_lo > t_hi && s_hi < t_lo {
        return Err(CalibrationError::OutOfRange { bound: hi, s_oceanic_at_bound: s_hi });
    }

    // Normalise so that mid-value response is compared monotonically:
    // after the endpoint check we know lo gives `s_lo` and hi gives
    // `s_hi`, with the target band strictly between them. Whichever
    // endpoint is below the band is `lo`.
    let target_mid = 0.5 * (t_lo + t_hi);
    if s_lo > s_hi {
        std::mem::swap(&mut lo, &mut hi);
        // Caller expects bisection on `k_spread`, not on flipped
        // endpoints; the swap is purely local to track which
        // endpoint is "below target".
    }

    let mut last_mid = 0.5 * (lo + hi);
    let mut last_s = f64::NAN;
    for _iter in 0..cfg.max_iters {
        let mid = 0.5 * (lo + hi);
        let s_mid = probe(mid);
        iterations.push(CalibrationIter { k_spread: mid, s_oceanic_mean: s_mid });
        last_mid = mid;
        last_s = s_mid;

        if s_mid >= t_lo && s_mid <= t_hi {
            return Ok(CalibrationResult {
                k_spread: mid,
                iterations,
                final_s_oceanic_mean: s_mid,
            });
        }

        // Monotone-increasing case: s_mid < target → raise `lo`;
        // s_mid > target → lower `hi`.
        if s_mid < target_mid {
            lo = mid;
        } else {
            hi = mid;
        }

        // Additional early termination if the k-interval collapses
        // below convergence_tol without hitting the target (pathology
        // like a step-function response). Surface this as
        // MaxItersReached, not OutOfRange — the bracket was valid, the
        // search just did not find a k inside the target band.
        if (hi - lo).abs() < cfg.convergence_tol {
            break;
        }
    }
    Err(CalibrationError::MaxItersReached { last_k: last_mid, last_s_oceanic: last_s })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn converges_on_linear_response() {
        // s_oceanic(k) = 0.4 * k — monotone, so target
        // `[0.18, 0.22]` → k ∈ [0.45, 0.55].
        let cfg = KSpreadCalibration::step5_default();
        let result = calibrate_k_spread(&cfg, |k| 0.4 * k).unwrap();
        assert!(result.final_s_oceanic_mean >= 0.18);
        assert!(result.final_s_oceanic_mean <= 0.22);
        assert!(result.k_spread >= 0.45 && result.k_spread <= 0.55);
    }

    #[test]
    fn max_iters_with_tight_tolerance_hits_max_iters_reached() {
        // Force `MaxItersReached` deterministically: straddle a
        // target band that is narrow enough for 2 bisections not to
        // reach. `s(k) = 0.2·k + 0.03` has `s(0.1)=0.05` (below
        // band) and `s(1.0)=0.23` (above band), so the bracket
        // straddles the target. With `target = [0.199, 0.201]` and
        // `max_iters = 2` the bisection visits `k ∈ {0.55, 0.775}`
        // with `s ∈ {0.14, 0.185}` — neither lands in the 0.002-wide
        // band, so termination is through `MaxItersReached`.
        let cfg = KSpreadCalibration {
            bracket: K_SPREAD_BRACKET,
            target_range: (0.199, 0.201),
            convergence_tol: 0.001,
            simulation_time: 3.0,
            max_iters: 2,
        };
        let response = |k: f64| 0.2 * k + 0.03;
        let result = calibrate_k_spread(&cfg, response);
        assert!(
            matches!(result, Err(CalibrationError::MaxItersReached { .. })),
            "expected MaxItersReached, got {:?}",
            result,
        );
    }

    #[test]
    fn target_above_bracket_high_is_out_of_range() {
        // s(k) = 0.1 * k always in [0.01, 0.10] over bracket
        // [0.1, 1.0]. Target (0.99, 1.01) is unreachable: s_hi is
        // well below t_lo, so OutOfRange fires. The error's `bound`
        // is the endpoint that produced the closest-to-target
        // response (`hi = 1.0`).
        let cfg = KSpreadCalibration {
            bracket: K_SPREAD_BRACKET,
            target_range: (0.99, 1.01),
            convergence_tol: 0.01,
            simulation_time: 3.0,
            max_iters: 20,
        };
        let response = |k: f64| 0.1 * k;
        match calibrate_k_spread(&cfg, response) {
            Err(CalibrationError::OutOfRange { bound, s_oceanic_at_bound }) => {
                assert!((bound - 1.0).abs() < 1e-12);
                assert!((s_oceanic_at_bound - 0.1).abs() < 1e-12);
            }
            other => panic!("expected OutOfRange, got {:?}", other),
        }
    }

    #[test]
    fn lower_endpoint_already_in_band_returns_immediately() {
        // A responder that returns 0.20 for any probe → the
        // bracket's low endpoint already lands in target, so the
        // bisection terminates on the first probe. The expected
        // `k_spread` is whatever `K_SPREAD_BRACKET.0` happens to be
        // — read it from the constant instead of hardcoding a stale
        // value.
        let cfg = KSpreadCalibration {
            bracket: K_SPREAD_BRACKET,
            target_range: (0.18, 0.22),
            convergence_tol: 0.01,
            simulation_time: 3.0,
            max_iters: 20,
        };
        let response = |_: f64| 0.20;
        let r = calibrate_k_spread(&cfg, response).unwrap();
        assert_eq!(r.iterations.len(), 1);
        assert_eq!(r.k_spread, K_SPREAD_BRACKET.0);
    }
}
