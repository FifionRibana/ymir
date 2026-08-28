//! Basal drag — velocity damping via the momentum-operator diagonal.
//!
//! Step 4 introduces the basal drag force
//!
//! ```text
//!   f_drag(i, j) = -Br · S̃(i, j)^exp · ṽ(i, j),        exp = 2 by default.
//! ```
//!
//! Unlike GPE (Step 2) or future slab-pull (Step 7), this term couples
//! linearly to `ṽ` and therefore joins the **operator**, not the RHS
//! (design decision D2 in the Step 4 spec / issue #87). Per-cell
//! contribution is a positive scalar `Br · S̃²` added to the diagonal
//! of the momentum operator `A`, interpolated from cell centres to the
//! MAC faces by the natural arithmetic 2-point cell-to-face average
//! already used by the viscous stencil.
//!
//! # Why `S̃²` and not `S̃`
//!
//! The design note §4.4 presents the canonical dimensionless form
//! `-Br · S̃ · ṽ`, but also notes that a `S̃²` weighting is physically
//! cleaner: it is `C∞`, vanishes naturally as `S̃ → 0`, and removes the
//! ad-hoc `S ≥ 0.3` hard threshold that the legacy solver used to skip
//! drag on thin oceanic crust. Step 4 adopts `S̃²` uniformly
//! (decision D1). The `s_exponent` field is kept configurable so the
//! test suite can probe against the canonical `S̃` variant if needed;
//! the physical default is `2.0`.
//!
//! # Jacobian
//!
//! The residual contribution `f_drag_residual(v) = Br · S̃² · v`
//! (sign flipped from the body-force form, since `A v = … - f_ext` and
//! drag sits on the operator side) has derivative with respect to `v`
//!
//! ```text
//!   ∂(Br · S̃² · v) / ∂v = Br · S̃² · I        (identity-scaled, diagonal).
//! ```
//!
//! `S̃` is **not** differentiated through: the thickness field is
//! advected separately and treated as frozen during a Newton solve.
//! So drag contributes nothing to `apply_tangent`; its entire effect
//! lives in the Picard block, which `apply_momentum` and
//! `momentum_diagonal` already handle once the drag contribution is
//! wired into their diagonals.
//!
//! # Preconditioner interaction (Case B)
//!
//! `stokes/precond.rs::VelocityJacobi::from_diagonal` takes the
//! diagonal as an **external slice**: it does not reconstruct. The
//! analytical reconstruction lives in
//! `stokes/operator.rs::momentum_diagonal`, which is a symbolic rewrite
//! of the viscous stencil parallel to `apply_momentum`. This is
//! **Case B** per the Step 4 spec (diagonal rebuilt analytically, not
//! extracted by probing the assembled operator). The drag augmentation
//! is therefore added in **both** `apply_momentum` and
//! `momentum_diagonal` with identical face-interpolation conventions;
//! consistency is enforced by `tests/v2_precond_drag_diagonal.rs`.
//!
//! # Zero cost when disabled
//!
//! [`BasalDragConfig::Disabled`] causes [`build_drag_diagonal_field`]
//! to return `None`. Callers propagate `None` as
//! `Option<&Field2D>` into the operator routines, which short-circuit
//! before any face loop. No allocation, no fma on `Br = 0.0` fields,
//! no dead branch inside the hot path. This is what the Step 4
//! regression ratio `[0.95, 1.05]` vs Step 3 requires (decision D4).

use crate::tectonics_v2::field::Field2D;

/// Basal-drag law parameters.
///
/// - `br`: the nondimensional basal drag number. Design target
///   `Br ∈ [0.01, 0.3]` (solver-scaling §4.4); default `0.05`
///   (decision D3).
/// - `s_exponent`: the exponent on `S̃`. Default `2.0` (decision D1);
///   configurable for MMS and diagnostic experiments.
#[derive(Clone, Copy, Debug)]
pub struct BasalDragLaw {
    pub br: f64,
    pub s_exponent: f64,
}

impl Default for BasalDragLaw {
    fn default() -> Self {
        Self { br: 0.05, s_exponent: 2.0 }
    }
}

/// Basal-drag configuration.
///
/// `Disabled` by-passes the contribution structurally (see the module
/// doc on zero-cost-when-disabled). `Enabled(law)` activates the
/// `Br · S̃^exp` diagonal augmentation.
#[derive(Clone, Copy, Debug)]
pub enum BasalDragConfig {
    Disabled,
    Enabled(BasalDragLaw),
}

impl BasalDragConfig {
    /// Short stable label for reports and logs.
    pub fn label(&self) -> &'static str {
        match self {
            BasalDragConfig::Disabled => "disabled",
            BasalDragConfig::Enabled(_) => "enabled",
        }
    }

    /// Parse a CLI token. Mirrors
    /// [`crate::tectonics_v2::presets::YieldingConfig::parse`] — the
    /// physics binary exposes `--basal-drag-config
    /// {enabled|disabled}`, and the `Enabled` variant receives the Br
    /// and exponent values through a separate `--br` flag so the two
    /// CLI flags compose cleanly.
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "disabled" | "off" => Ok(BasalDragConfig::Disabled),
            "enabled" | "on" => Ok(BasalDragConfig::Enabled(BasalDragLaw::default())),
            other => Err(format!(
                "unknown --basal-drag-config value '{}'; expected disabled|enabled",
                other,
            )),
        }
    }

    /// Human-readable descriptor threaded into `SolverConfigDump`.
    pub fn describe(&self) -> String {
        match self {
            BasalDragConfig::Disabled => "Disabled".to_string(),
            BasalDragConfig::Enabled(law) => {
                format!("Enabled (Br = {:.3}, S exponent = {:.1})", law.br, law.s_exponent,)
            }
        }
    }
}

/// Per-cell basal-drag diagonal coefficient `Br · S̃^exp`.
///
/// Exposed for the unit test that pins the algebraic form; normal
/// callers use [`build_drag_diagonal_field`] and let the operator
/// layer consume the resulting `Option<&Field2D>`.
#[inline]
pub fn drag_diagonal_at_cell(br: f64, s: f64, exponent: f64) -> f64 {
    br * s.powf(exponent)
}

/// Assemble the cell-centered basal-drag diagonal field `Br · S̃^exp`.
///
/// Returns `None` when `cfg == Disabled` — the caller then propagates
/// `None: Option<&Field2D>` into the operator routines, which
/// structurally short-circuit before any face loop. When enabled, the
/// returned field has the same shape as `s`.
///
/// The field is cell-centered; `apply_momentum` and `momentum_diagonal`
/// take care of the cell-to-face arithmetic averaging at matvec time.
pub fn build_drag_diagonal_field(cfg: &BasalDragConfig, s: &Field2D) -> Option<Field2D> {
    match cfg {
        BasalDragConfig::Disabled => None,
        BasalDragConfig::Enabled(law) => {
            let nx = s.nx();
            let ny = s.ny();
            let mut out = Field2D::new(nx, ny);
            for j in 0..ny {
                for i in 0..nx {
                    out.set(i, j, drag_diagonal_at_cell(law.br, s.get(i, j), law.s_exponent));
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
    fn default_law_is_step4_baseline() {
        let law = BasalDragLaw::default();
        assert_eq!(law.br, 0.05);
        assert_eq!(law.s_exponent, 2.0);
    }

    #[test]
    fn drag_diagonal_at_cell_matches_algebraic_form() {
        let br = 0.05;
        let s = 1.2_f64;
        let exp = 2.0;
        let expected = br * s * s;
        assert!((drag_diagonal_at_cell(br, s, exp) - expected).abs() < 1e-14);
    }

    #[test]
    fn drag_diagonal_annihilates_at_zero_s() {
        // Oceanic annihilation: S̃ = 0 ⇒ drag = 0 regardless of Br.
        assert_eq!(drag_diagonal_at_cell(0.3, 0.0, 2.0), 0.0);
        assert_eq!(drag_diagonal_at_cell(0.3, 0.0, 1.0), 0.0);
    }

    #[test]
    fn build_drag_diagonal_field_returns_none_when_disabled() {
        let s = Field2D::filled(6, 6, 1.0);
        assert!(build_drag_diagonal_field(&BasalDragConfig::Disabled, &s).is_none());
    }

    #[test]
    fn build_drag_diagonal_field_populates_per_cell_values() {
        let mut s = Field2D::new(4, 4);
        for j in 0..4 {
            for i in 0..4 {
                s.set(i, j, 1.0 + 0.1 * (i + j) as f64);
            }
        }
        let law = BasalDragLaw { br: 0.1, s_exponent: 2.0 };
        let cfg = BasalDragConfig::Enabled(law);
        let out = build_drag_diagonal_field(&cfg, &s).expect("enabled → Some");
        for j in 0..4 {
            for i in 0..4 {
                let expected = 0.1 * s.get(i, j).powi(2);
                assert!(
                    (out.get(i, j) - expected).abs() < 1e-14,
                    "cell ({},{}) = {}, expected {}",
                    i,
                    j,
                    out.get(i, j),
                    expected,
                );
            }
        }
    }

    #[test]
    fn label_matches_variant() {
        assert_eq!(BasalDragConfig::Disabled.label(), "disabled");
        assert_eq!(BasalDragConfig::Enabled(BasalDragLaw::default()).label(), "enabled");
    }

    #[test]
    fn parse_accepts_enabled_and_disabled_tokens() {
        assert!(matches!(BasalDragConfig::parse("disabled").unwrap(), BasalDragConfig::Disabled));
        assert!(matches!(BasalDragConfig::parse("off").unwrap(), BasalDragConfig::Disabled));
        assert!(matches!(BasalDragConfig::parse("enabled").unwrap(), BasalDragConfig::Enabled(_)));
        assert!(matches!(BasalDragConfig::parse("on").unwrap(), BasalDragConfig::Enabled(_)));
    }

    #[test]
    fn parse_rejects_garbage() {
        assert!(BasalDragConfig::parse("truthy").is_err());
        assert!(BasalDragConfig::parse("").is_err());
    }

    #[test]
    fn describe_emits_br_and_exponent_when_enabled() {
        let law = BasalDragLaw { br: 0.05, s_exponent: 2.0 };
        let desc = BasalDragConfig::Enabled(law).describe();
        assert!(desc.contains("Enabled"));
        assert!(desc.contains("0.050"));
        assert!(desc.contains("2.0"));
        assert_eq!(BasalDragConfig::Disabled.describe(), "Disabled");
    }
}
