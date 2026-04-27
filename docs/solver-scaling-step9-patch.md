# Step 9 — §4.10 amendments and clarifications

This document collects three formal patches to `solver-scaling.md`
§4.10 produced during Step 9 implementation. They will be folded
into `solver-scaling.md` itself when the Step 9 PR merges.

## Patch 1 — `Cr` parameter formalisation

`solver-scaling.md` §4.10 references the cratonic regions
qualitatively but does not define a parameter that controls the
cratonic *extent*. Step 9 introduces such a parameter,
documented here.

> **`Cr` — cratonic fraction parameter.** A nondimensional handle
> introduced for the Step 9 implementation, **not a number from
> the geodynamics literature**. It is constructed in the spirit
> of §5.1 dimensionless numbers but is local to the Step 9
> design.
>
> Operationally, `Cr ∈ [0.1, 0.5]` (default `0.3`) sets the target
> fraction of plate area occupied by the cratonic core for plates
> large enough to host one. The cratonic factor field
> `cratonic_factor[i] ∈ [0, 1]` is built per Step 9 D2: a
> distance-from-boundary BFS, normalized per-plate, and a
> smoothstep transition centred at `d_mid = 1 - sqrt(Cr)`.
>
> **Note on `L_plate`.** §4.10 (per Step 9 issue D2/D3) writes
> `d_mid = R · (1 - sqrt(Cr))` for a circular plate of radius R,
> then "generalised to non-circular plates via L_plate =
> sqrt(plate_area)". Reading these together literally is
> self-inconsistent: for a circle of radius R, area is πR² so
> `sqrt(plate_area) = R · sqrt(π) ≈ 1.77 R`, not R.
> Substituting that into the formula collapses the cratonic core
> to a few percent of plate area for any realistic `Cr` (verified
> by the calibration probe `tests/v2_cratonic_normalization_probe`).
>
> The geometric intent is that `L_plate` represents the
> *characteristic distance from boundary to plate interior*,
> i.e. the inradius generalisation. For a circular plate that is
> exactly `R`; for a square of side `L` it is `L/2`; for an
> irregular Voronoï plate it is naturally the maximum BFS depth
> attained inside the plate. The Step 9 implementation uses
> `L_plate = max BFS depth in plate`, under which the formula
> `d_mid_normalized = 1 - sqrt(Cr)` is geometrically correct
> and exactly recovers the cratonic-core area fraction `Cr` for
> both circular and square reference shapes.
>
> **Empirical realisation.** The calibration probe across 31
> random Voronoï seeds at 64² with default settings produces a
> realised `cratonic_cell_fraction / (Cr · continental_fraction)`
> ratio of ~1.13 mean (vs ~1.26 with the literal `sqrt(area)`
> reading). Per-seed variation can reach up to 1.5× because
> Voronoï plates are neither circular nor square — irregular
> blob shapes give larger cratonic fractions than the regular-
> shape derivation predicts. Step 9 acceptance criterion #8
> ("`cratonic_cell_fraction ≈ Cr · continental_fraction` within
> ±20 %") is therefore evaluated **in aggregate across multiple
> seeds**, not strictly per-instance — the user-acknowledged
> reformulation is "mean across seeds within ±20 %, per-seed
> dispersion up to 1.5× documented".

## Patch 2 — `B_factor` primary-mechanism amendment

The original Step 9 issue D1 articulated the primary plastic-
immunity mechanism via the formula

```text
yield_stress[i] = Bi · (cratonic_factor[i]
                      + (1 - cratonic_factor[i]) · weakening(plastic_strain[i]))
```

which assumes a `weakening` function depending on accumulated
plastic strain. In the current milestone, plastic memory is **not**
implemented (deferred) and `weakening` is implicitly the constant
`1.0` everywhere. Substituting this collapses the formula to
`yield_stress = Bi · 1 = Bi` for any `cratonic_factor`,
making the primary mechanism a no-op in the current scope.

Diagnostic on the Step 8-shape immunity test (32², mantle on,
slab off) with K = 5 (secondary mechanism) and the original D1
formula showed `peak_yielding_in_craton = 0.99` — cratons yielded
essentially everywhere. The K viscous mult cannot suppress
viscoplastic yielding `η_p = Bi/(2(ε̇+ε̇_min))` in saturated regimes
because `η_p` becomes the soft-min branch when ε̇ is large
regardless of how much the secondary K stiffens the viscous
branch. Acceptance #6 (`peak_yielding_in_craton ≤ 0.01`) cannot
be reached with the original formula in this milestone.

> **Amendment — operational form of the primary mechanism.** The
> primary mechanism "cratons immune to plastic weakening"
> generalises to "cratons have an *elevated* yield strength"
> via a new parameter `B_factor ∈ [3, 10]`, default `8`,
> multiplying `Bi`. The default is derived from the analytical
> threshold `B > η_v / (2·K·η_p_default) ≈ 6.1` in activated
> regimes (`peak|v| ~ O(1)`, ε̇ large) and validated empirically
> by the `B_factor` sweep on the Step 8 shape immunity test.
> Lower values (B_factor = 3–5) are sufficient in non-activated
> regimes but produce yielding leakage in cratons under saturated
> conditions. The implementation provides full configurability;
> defaults are chosen to satisfy immunity acceptance across all
> tested regimes.
>
> ```text
> yield_stress[i] = Bi · (1 + (B_factor - 1) · cratonic_factor[i])
>                      · weakening(plastic_strain[i])
> ```
>
> Limits:
>
> - `cratonic_factor = 0` (mobile cells): `yield_stress = Bi · weakening`
>   (= Bi today; plastic-memory-modulated when implemented).
> - `cratonic_factor = 1` (full cratonic core): `yield_stress =
>   B_factor · Bi · weakening` (= B_factor · Bi today).
>
> In the absence of plastic memory (current milestone), `weakening = 1`
> everywhere and this elevation **is** the operational form of the
> primary mechanism. When plastic memory is later implemented, the
> formula retains `weakening` modulating mobile belts; cratons'
> `plastic_strain` stays zero by design D1 so `weakening(0) = 1`
> and `B_factor · Bi` survives unmodified.
>
> `B_factor = 1` reduces the formula to the identity (the
> pre-amendment behaviour) and is preserved as a degenerate
> diagnostic configuration that exercises the secondary K
> mechanism alone.
>
> **Range justification.** Lower bound `B_factor ≥ 3` matches
> the K range for symmetric primary/secondary configurability.
> Upper bound `B_factor ≤ 10` keeps the κ(A) impact
> manageable: `B_factor` enters multiplicatively in `η_p`, not
> on the full η_eff, so its effective contribution to η_eff is
> bounded by `min(η_v, B_factor · η_p_default)`. Beyond 10 the
> cratonic branch saturates (η_v starts capping) and additional
> Bi elevation gives diminishing returns while pushing the
> conditioning further. If a sweep at the top of range fails
> acceptance #6, a remontée for architectural review is
> required — extending B_factor beyond 10 silently is an anti-
> pattern.

## Patch 3 — Combined-mechanism characterisation

> **Combined effect.** With both primary (`B_factor`) and
> secondary (`K`) mechanisms active in a fully cratonic cell
> (`cratonic_factor = 1`), the effective viscosity is
> `η_eff = K · soft_min(η_v, B_factor · η_p_default)`. The
> conditioning impact stays manageable because:
>
> 1. **Primary**: Bi elevation acts only on `η_p` *pre-blend*.
>    Once the soft-min selects the viscous branch (`η_v` wins),
>    further Bi increase has zero effect.
> 2. **Secondary**: K acts post-blend on the result. Its effect
>    on κ(A) scales linearly with K and is the dominant
>    contributor to viscosity contrast across the cratonic
>    boundary.
>
> The smoothstep transition `cratonic_factor[i]` is the same
> field for both mechanisms, so the cratonic-mobile transition
> is geometrically consistent. Smoothness acceptance is checked
> via the `eta_multiplier` ratio at adjacent cells crossing the
> boundary (acceptance #3, target `≤ K · 1.05`); this
> isolates the cratonic-induced contrast from the underlying
> `η_law(ε̇_II)` gradient that exists at any boundary between
> dynamically-different regions and would confound an `η_eff`
> ratio.
