//! Mantle convection forcing — Step 8.
//!
//! Installs the **initiator** of the mechanism hierarchy. Unlike
//! slab-pull (Step 7, an amplifier that cannot bootstrap out of
//! floor-domination), mantle forcing imposes a velocity bias
//! `v_mantle(x) = Mf · v_pattern(x)` independently of the
//! dynamic state. This is what breaks the quiescent fixed point
//! and lets the yielding checkpoint (transported since Step 3)
//! resolve.
//!
//! Formulation (§4.9 + D1):
//!
//! ```text
//!   f_mantle = coupling · S̃ · (Mf · v_pattern − v_solved)
//! ```
//!
//! Refactored into operator + RHS contributions for exact
//! self-consistency at every Newton outer iteration, following
//! the Step 4 basal-drag pattern:
//!
//! ```text
//!   RHS part (handled by MantleForce):      coupling · S̃ · Mf · v_pattern
//!   Operator diagonal (handled by harness): coupling · S̃ · I
//! ```
//!
//! The operator-diagonal part sums with the basal-drag diagonal
//! (`Br · S̃²`) into a single `total_diag` field passed through
//! `solve_sheet` unchanged. The inner CG does not see any
//! mantle-specific dispatch — it solves `(A(v;η) + total_diag · I) δv = -r`
//! where `total_diag = drag_diag + mantle_diag`. Self-consistency
//! of `v_solved` with `f_mantle` follows from Newton's own
//! convergence on the augmented operator — no callback-based RHS
//! recomputation is required.
//!
//! Pattern construction (D2):
//!
//! The mantle velocity field is built from a Fourier stream
//! function `ψ(x, y)` via the staggered discrete curl
//! `v_mantle = (∂ψ/∂y, −∂ψ/∂x)`. `ψ` is sampled at **grid nodes**
//! `(i·dx, j·dy)` — i.e. cell corners — not at cell centres.
//! This is essential because the discrete curl-of-nodal-ψ on a
//! MAC-staggered velocity grid is **exactly divergence-free** by
//! algebraic cancellation, independent of grid resolution. A
//! cell-centered `ψ` would leave O(dx²) residual divergence and
//! violate the strict acceptance `div_v_mantle_max < 10⁻¹⁰`.
//!
//! See [`stream_function::generate_stream_function`] and
//! [`pattern::build_mantle_pattern`] for the construction.

pub mod pattern;
pub mod stream_function;

pub use pattern::{build_mantle_pattern, MantlePattern};
pub use stream_function::{
    generate_stream_function, generate_stream_function_at_time, StreamFunctionBuilder,
    StreamFunctionConfig,
};

use crate::tectonics_v2::field::Field2D;

/// Default forcing amplitude. Sits in the middle of the §4.9
/// target band `Mf ∈ [0.3, 2.0]`.
pub const MF_DEFAULT: f64 = 1.0;

/// Default coupling `c ∈ [0.1, 10.0]`. `c = 1.0` means a plate
/// takes roughly one τ* to match the mantle pattern under pure
/// relaxation (no other forces).
pub const COUPLING_DEFAULT: f64 = 1.0;

/// Default number of Fourier modes in the stream function.
pub const NUM_MODES_DEFAULT: usize = 6;

/// Default evolution rate. `0.0` = static pattern (Step 8 baseline,
/// bit-identical with pre-Step-12-R6 behaviour). Non-zero activates
/// the Step 12 R6 Phys.A phase-drift evolution: every Fourier mode's
/// `(φx, φy)` shifts by `ω · t_nondim` with `ω = evolution_rate · TAU`,
/// while wave numbers, amplitudes and the t=0 normalisation stay
/// frozen.
pub const EVOLUTION_RATE_DEFAULT: f64 = 0.0;

/// Mantle forcing configuration.
///
/// `Disabled` bypasses the whole pipeline (pattern not generated,
/// `MantleForce` not pushed into the `ForceSum`, diagonal
/// contribution not summed). The Step 8 regression verifies
/// zero-cost-when-disabled in scalar parity with Step 7 physics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MantleConfig {
    Disabled,
    Enabled {
        /// §4.9 amplitude, target band `[0.3, 2.0]`.
        mf: f64,
        /// Coupling `c`, target band `[0.1, 10.0]`.
        coupling: f64,
        /// Number of Fourier modes in the stream function.
        num_modes: usize,
        /// Seed for mode tuple generation (wave numbers,
        /// amplitudes, phases). Independent of the baseline's
        /// main `WorldSeed` so the mantle pattern can be pinned
        /// independently for sweeps.
        seed: u64,
        /// Evolution rate. `0.0` = static pattern (Step 8
        /// baseline). Non-zero = Step 12 R6 Phys.A phase drift —
        /// each step rebuilds ψ with all mode phases shifted by
        /// `evolution_rate · TAU · t_nondim`. Wave numbers,
        /// amplitudes and t=0 normalisation are frozen.
        evolution_rate: f64,
    },
}

impl MantleConfig {
    pub fn label(&self) -> &'static str {
        match self {
            MantleConfig::Disabled => "disabled",
            MantleConfig::Enabled { .. } => "enabled",
        }
    }

    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "disabled" | "off" => Ok(MantleConfig::Disabled),
            "enabled" | "on" => Ok(MantleConfig::Enabled {
                mf: MF_DEFAULT,
                coupling: COUPLING_DEFAULT,
                num_modes: NUM_MODES_DEFAULT,
                seed: 42,
                evolution_rate: EVOLUTION_RATE_DEFAULT,
            }),
            other => Err(format!(
                "unknown --mantle-config value '{}'; expected disabled|enabled",
                other,
            )),
        }
    }

    pub fn describe(&self) -> String {
        match self {
            MantleConfig::Disabled => "Disabled".to_string(),
            MantleConfig::Enabled { mf, coupling, num_modes, seed, evolution_rate } => format!(
                "Enabled (Mf = {:.3}, coupling = {:.3}, num_modes = {}, seed = {}, evolution_rate = {:.3})",
                mf, coupling, num_modes, seed, evolution_rate,
            ),
        }
    }

    pub fn mf_or_zero(&self) -> f64 {
        match self {
            MantleConfig::Disabled => 0.0,
            MantleConfig::Enabled { mf, .. } => *mf,
        }
    }

    pub fn coupling(&self) -> f64 {
        match self {
            MantleConfig::Disabled => 0.0,
            MantleConfig::Enabled { coupling, .. } => *coupling,
        }
    }
}

impl Default for MantleConfig {
    fn default() -> Self {
        MantleConfig::Disabled
    }
}

/// Build the cell-centered `coupling · S̃` diagonal field.
///
/// Step 4 basal-drag pattern: the field is summed into the
/// momentum operator diagonal (and preconditioner). `None` when
/// disabled → no cost; `Some(&field)` when enabled → `coupling`
/// times the cell-centered `S̃` (linear in `S̃`, NOT `S̃²` — see
/// D5 rationale: `S̃` is a weight on a body-force term, not a
/// viscous coefficient).
///
/// The caller sums `drag_diag + mantle_diag` into a single
/// `total_diag` slice before passing to `solve_sheet`. When one
/// is `None` the other passes through unchanged — this preserves
/// the Step 7 bit-identity in the regression run.
pub fn build_mantle_diagonal_field(cfg: &MantleConfig, s: &Field2D) -> Option<Field2D> {
    match cfg {
        MantleConfig::Disabled => None,
        MantleConfig::Enabled { coupling, .. } => {
            let nx = s.nx();
            let ny = s.ny();
            let mut out = Field2D::new(nx, ny);
            let c = *coupling;
            for j in 0..ny {
                for i in 0..nx {
                    out.set(i, j, c * s.get(i, j));
                }
            }
            Some(out)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_is_disabled() {
        assert_eq!(MantleConfig::default(), MantleConfig::Disabled);
    }

    #[test]
    fn parse_disabled_and_enabled() {
        assert_eq!(MantleConfig::parse("disabled").unwrap(), MantleConfig::Disabled);
        match MantleConfig::parse("enabled").unwrap() {
            MantleConfig::Enabled { mf, coupling, num_modes, evolution_rate, .. } => {
                assert_eq!(mf, MF_DEFAULT);
                assert_eq!(coupling, COUPLING_DEFAULT);
                assert_eq!(num_modes, NUM_MODES_DEFAULT);
                assert_eq!(evolution_rate, EVOLUTION_RATE_DEFAULT);
            }
            MantleConfig::Disabled => panic!("expected Enabled"),
        }
    }

    #[test]
    fn baseline_mf_and_coupling_in_spec_bands() {
        assert!(MF_DEFAULT >= 0.3 && MF_DEFAULT <= 2.0);
        assert!(COUPLING_DEFAULT >= 0.1 && COUPLING_DEFAULT <= 10.0);
    }

    #[test]
    fn diagonal_is_none_when_disabled() {
        let s = Field2D::filled(4, 4, 1.0);
        assert!(build_mantle_diagonal_field(&MantleConfig::Disabled, &s).is_none());
    }

    #[test]
    fn diagonal_scales_linearly_with_s_and_coupling() {
        let mut s = Field2D::new(3, 3);
        for j in 0..3 {
            for i in 0..3 {
                s.set(i, j, 0.1 * (i + j + 1) as f64);
            }
        }
        let cfg = MantleConfig::Enabled {
            mf: 1.0, coupling: 2.5, num_modes: 1, seed: 1, evolution_rate: 0.0,
        };
        let d = build_mantle_diagonal_field(&cfg, &s).expect("enabled → Some");
        for j in 0..3 {
            for i in 0..3 {
                let expected = 2.5 * s.get(i, j);
                assert!((d.get(i, j) - expected).abs() < 1e-14);
            }
        }
    }
}
