//! Step 12 R6.2 — `MantleConfig::Enabled.evolution_rate` wiring tests.
//!
//! Four regression contracts:
//!
//! 1. **Static when `evolution_rate = 0`**: `sample_at_time(t=N·dt, evo=0)`
//!    is bit-identical to `sample_at_time(t=0, evo=0)` for any N. The
//!    `evolution_rate == 0` early-out in `StreamFunctionBuilder` avoids
//!    the `drift_modes` allocation, so output bits are exact.
//!
//! 2. **Bit-identical public wrapper**: `generate_stream_function(seed)`
//!    keeps the pre-R6 output bit-for-bit. Internally it now routes
//!    through `StreamFunctionBuilder::new(...).sample_at_time(0.0, 0.0)`;
//!    the byte equality with a freshly-built builder confirms the
//!    refactor preserved the wire-format.
//!
//! 3. **Measurable evolution under non-zero `evolution_rate`**: at
//!    `evolution_rate = 0.5` and `t = 1.0` the phase offset is `π` and
//!    every mode's sines flip sign — the L∞ deviation must exceed a
//!    comfortable floor (0.1 here, well above the discretisation noise).
//!
//! 4. **Div-freeness preserved across the drift**: `MantlePattern`
//!    rebuilt via `rebuild_from_psi` from a phase-drifted ψ keeps
//!    `pattern_div_max < 1e-10` for every step in a 20-step sweep —
//!    the Step 8 strict acceptance survives at every t.

use ymir_core::tectonics_v2::field::PeriodicIndex;
use ymir_core::tectonics_v2::mantle::pattern::pattern_div_max;
use ymir_core::tectonics_v2::mantle::{
    StreamFunctionBuilder, StreamFunctionConfig, build_mantle_pattern, generate_stream_function,
    generate_stream_function_at_time,
};

const SEED: u64 = 42;
const NUM_MODES: usize = 6;

/// Contract 1 — `evolution_rate = 0` keeps the pattern strictly
/// constant across any number of step rebuilds.
#[test]
fn evolution_rate_zero_pattern_constant_multistep() {
    let n: usize = 32;
    let cfg = StreamFunctionConfig { num_modes: NUM_MODES, seed: SEED };
    let builder = StreamFunctionBuilder::new(n, n, &cfg);
    let psi_t0 = builder.sample_at_time(n, n, 0.0, 0.0);
    // Sweep ten "step-equivalent" times. With evolution_rate = 0,
    // every sample must be byte-equal to t=0.
    for k in 1..=10 {
        let t = k as f64 * 0.02; // dt_nondim ≈ 0.02 in the harness
        let psi_t = builder.sample_at_time(n, n, t, 0.0);
        assert_eq!(
            psi_t.data(),
            psi_t0.data(),
            "evolution_rate = 0 at t = {} should be bit-identical to t = 0",
            t,
        );
    }
}

/// Contract 2 — the public `generate_stream_function` wrapper still
/// produces the same field bytes as a builder built and sampled in
/// the new code path. This protects the pre-R6 bit-identical contract
/// for every existing caller of the legacy entry point.
#[test]
fn generate_stream_function_matches_builder_sample_at_t_zero() {
    for &(nx, ny) in &[(16_usize, 16_usize), (32, 32), (64, 64), (24, 40)] {
        let cfg = StreamFunctionConfig { num_modes: NUM_MODES, seed: SEED };
        let legacy = generate_stream_function(nx, ny, &cfg);
        let via_builder = StreamFunctionBuilder::new(nx, ny, &cfg).sample_at_time(nx, ny, 0.0, 0.0);
        assert_eq!(
            legacy.data(),
            via_builder.data(),
            "generate_stream_function diverged from builder.sample_at_time(0, 0) \
             at ({nx}, {ny})",
        );
        // Also confirm the t=0 convenience function matches.
        let one_shot = generate_stream_function_at_time(nx, ny, &cfg, 0.0, 0.0);
        assert_eq!(
            legacy.data(),
            one_shot.data(),
            "generate_stream_function_at_time(t=0, evo=0) diverged at ({nx}, {ny})",
        );
        // Sanity: pre-R6 contract `max|ψ| = 1`.
        let max_abs = legacy.data().iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        assert!(
            (max_abs - 1.0).abs() < 1e-12,
            "max|ψ| = {} at ({nx}, {ny}), expected 1 ± 1e-12",
            max_abs,
        );
    }
}

/// Contract 3 — non-zero `evolution_rate` produces a *measurably*
/// different field after enough simulated time.
///
/// Phys.A structural note: because the same offset shifts both
/// `(φx, φy)`, the separable product `sin(arg_x + φ_x) · sin(arg_y + φ_y)`
/// is invariant at `phase_offset ∈ {π, 2π, …}` (both sines flip; their
/// product doesn't). The drifted field's effective period in phase is
/// therefore `π`, not `TAU`. To stay clear of the identity points we
/// probe at `phase_offset = π / 2` (= quadrature shift, `sin → cos`):
/// `evolution_rate = 0.25`, `t = 1.0` ⇒ `0.25 · TAU · 1.0 = π/2`.
///
/// At quadrature, no clean algebraic identity holds; we assert the
/// raw L∞ distance is above a comfortable noise floor.
#[test]
fn evolution_rate_nonzero_evolves_measurably() {
    let n: usize = 32;
    let cfg = StreamFunctionConfig { num_modes: NUM_MODES, seed: SEED };
    let builder = StreamFunctionBuilder::new(n, n, &cfg);
    let psi_t0 = builder.sample_at_time(n, n, 0.0, 0.25);
    let psi_t1 = builder.sample_at_time(n, n, 1.0, 0.25);
    let l_inf = psi_t0
        .data()
        .iter()
        .zip(psi_t1.data().iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        l_inf > 0.1,
        "L∞ distance between ψ(t=0) and ψ(t=1, evo=0.25) is {l_inf}, \
         expected > 0.1 — quadrature shift should give visibly different field",
    );

    // Equally important: confirm `phase_offset = π` IS a degenerate
    // identity point (documents the Phys.A structural caveat). With
    // evolution_rate = 0.5, t = 1.0 ⇒ phase_offset = π, and both
    // sines flip sign so the product is preserved.
    let psi_pi = builder.sample_at_time(n, n, 1.0, 0.5);
    let psi_zero = builder.sample_at_time(n, n, 0.0, 0.5);
    let identity_l_inf = psi_zero
        .data()
        .iter()
        .zip(psi_pi.data().iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        identity_l_inf < 1e-10,
        "expected phase_offset = π to be a degenerate identity of the \
         separable product (Phys.A structural caveat), but L∞ = {identity_l_inf}",
    );
}

/// Contract 4 — phase drift preserves the strict div-freeness
/// acceptance from Step 8 (`div_v_mantle_max < 1e-10`) at every step
/// in a 20-step sweep. Exercises `MantlePattern::rebuild_from_psi`
/// (the R6 in-place rebuild path used by the harness loop).
#[test]
fn evolution_rate_preserves_div_freeness_multistep() {
    let n: usize = 32;
    let dx = 1.0 / n as f64;
    let idx_x = PeriodicIndex::new(n);
    let idx_y = PeriodicIndex::new(n);
    let cfg = StreamFunctionConfig { num_modes: NUM_MODES, seed: SEED };
    let evolution_rate = 0.1_f64;
    let dt = 0.02_f64;
    let builder = StreamFunctionBuilder::new(n, n, &cfg);
    let psi0 = builder.sample_at_time(n, n, 0.0, evolution_rate);
    let mut pattern = build_mantle_pattern(&psi0, dx, dx, &idx_x, &idx_y);
    let div0 = pattern_div_max(&pattern, dx, dx, &idx_x, &idx_y);
    assert!(div0 < 1e-10, "div_v_mantle_max at t=0 = {div0:.3e}, expected < 1e-10",);
    for step in 1..=20 {
        let t = step as f64 * dt;
        let psi_t = builder.sample_at_time(n, n, t, evolution_rate);
        pattern.rebuild_from_psi(&psi_t, dx, dx, &idx_x, &idx_y);
        let div = pattern_div_max(&pattern, dx, dx, &idx_x, &idx_y);
        assert!(
            div < 1e-10,
            "step {step} (t = {t:.3}): div_v_mantle_max = {div:.3e}, \
             expected < 1e-10 — phase drift broke div-freeness",
        );
    }
}
