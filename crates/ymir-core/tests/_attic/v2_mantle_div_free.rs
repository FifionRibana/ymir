//! Discrete div-free verification of the mantle pattern
//! construction.
//!
//! `v_mantle = curl(ψ ẑ)` built from a nodal `ψ` on the MAC
//! staggered grid must produce `div(v_mantle) ≡ 0` to f64
//! precision, at any resolution. This is the guarantee that
//! mantle forcing does not pollute Step 6's boundary detection
//! from `div(v_solved)`.
//!
//! Acceptance (Step 8 spec): `div_v_mantle_max < 10⁻¹⁰` at
//! N ∈ {64, 128, 256}.

use ymir_core::tectonics_v2::field::PeriodicIndex;
use ymir_core::tectonics_v2::mantle::pattern::pattern_div_max;
use ymir_core::tectonics_v2::mantle::{
    StreamFunctionConfig, build_mantle_pattern, generate_stream_function,
};

#[test]
fn div_is_below_strict_threshold_at_64_128_256() {
    for &n in &[64_usize, 128, 256] {
        let dx = 1.0 / n as f64;
        let idx_x = PeriodicIndex::new(n);
        let idx_y = PeriodicIndex::new(n);
        let psi = generate_stream_function(n, n, &StreamFunctionConfig { num_modes: 6, seed: 42 });
        let pattern = build_mantle_pattern(&psi, dx, dx, &idx_x, &idx_y);
        let div = pattern_div_max(&pattern, dx, dx, &idx_x, &idx_y);
        eprintln!("N = {}, div_v_mantle_max = {:.3e}", n, div);
        assert!(
            div < 1.0e-10,
            "N={}: div_v_mantle_max = {:.3e} exceeds 1e-10 — \
             mantle pattern would pollute div(v) diagnostic",
            n,
            div,
        );
    }
}

/// The construction is seed-independent: changing the seed
/// generates a different pattern but div must still be zero.
#[test]
fn div_free_under_several_seeds() {
    let n = 128;
    let dx = 1.0 / n as f64;
    let idx_x = PeriodicIndex::new(n);
    let idx_y = PeriodicIndex::new(n);
    for &seed in &[1_u64, 7, 42, 100, 9999] {
        let psi = generate_stream_function(n, n, &StreamFunctionConfig { num_modes: 6, seed });
        let pattern = build_mantle_pattern(&psi, dx, dx, &idx_x, &idx_y);
        let div = pattern_div_max(&pattern, dx, dx, &idx_x, &idx_y);
        assert!(div < 1.0e-10, "seed={}: div = {:.3e}", seed, div,);
    }
}

/// The number of modes does not affect div-freeness (curl of a
/// sum of stream functions is still div-free).
#[test]
fn div_free_under_various_mode_counts() {
    let n = 64;
    let dx = 1.0 / n as f64;
    let idx_x = PeriodicIndex::new(n);
    let idx_y = PeriodicIndex::new(n);
    for &num_modes in &[1_usize, 3, 6, 12, 20] {
        let psi = generate_stream_function(n, n, &StreamFunctionConfig { num_modes, seed: 42 });
        let pattern = build_mantle_pattern(&psi, dx, dx, &idx_x, &idx_y);
        let div = pattern_div_max(&pattern, dx, dx, &idx_x, &idx_y);
        assert!(div < 1.0e-10, "num_modes={}: div = {:.3e}", num_modes, div,);
    }
}
