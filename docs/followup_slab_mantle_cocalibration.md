# Follow-up: Slab+Mantle co-calibration

**Status:** draft for GitHub issue, post-Step 8 merge.

**Depends on:** Step 8 (issue #95) shipped.
**Blocks:** N/A (Step 9 cratonic immunity can proceed on Step 8's
mantle-only base; this issue is orthogonal to that step's scope).

## Problem

At Step 7 we established that slab-pull alone cannot bootstrap the
system out of floor-domination (§4.8 revision, `slab-pull is an
amplifier, not an initiator`). At Step 8 we established that mantle
forcing provides the initiator and resolves the yielding checkpoint
transported since Step 3 — **but only with slab-pull held Disabled**.

Running the nominal spec configuration — Step 7 physics setup (slab
Enabled at §4.8 baseline) plus mantle Enabled at §4.9 baseline —
produces catastrophic numerical divergence within ~15 timesteps at
64², Mf=1.0, coupling=1.0:

| steps | peak\|v_solved\| | peak\|f_slab\| | alignment |
|---|---|---|---|
| 5 | `9.6e0` | `9.8e0` | `+0.22` |
| 10 | `3.3e1` | `5.5e1` | `+0.23` |
| 15 | `1.5e7` | `1.0e6` | `−48` |
| 20 | `7.9e14` | `4.0e13` | `−1.9e9` |

Full diagnostic in
[`docs/reports/step8_physics_report.md` §Slab+Mantle interaction
instability finding](reports/step8_physics_report.md). Reproducible
via the `#[ignore]`-d test
[`crates/ymir-core/tests/v2_mantle_runaway_diagnostic.rs`](../crates/ymir-core/tests/v2_mantle_runaway_diagnostic.rs).

## Root cause (physical, not implementation)

Closed-loop gain analysis in the mantle-activated regime. Once mantle
forcing pulls `v ~ O(Mf) = O(1)`, the power-law rheology exits the
floor band, `η_newton → O(1)`, and the viscous diagonal at wave
number `k=1` on a 64² grid is `2·η·k² ≈ 80`. In the same regime, the
discrete divergence operator in `Q_sub_conv = k_slab · max(0, −div v)`
amplifies by `2/dx = 128`. Then

```text
m_subducted ≈ Q · τ_slab = 64 · v
f_slab      = Sp · m ≈ 1.5 · 64 · v = 96 · v

G_activated = 96 / 80 ≈ 1.2   >   1
```

— linear instability in the activated regime. The §4.8 target band
`Sp ∈ [0.5, 3]` was calibrated against quiescent-regime balance
assumptions and is **not co-calibrated** with §4.9's `Mf ∈ [0.3, 2]`
band in the activated regime.

Per the spec language in
[`docs/solver-scaling.md` §4.8 activation-regime stability
constraint](solver-scaling.md), stability requires

```text
Sp · k_slab_accum · τ_slab · (2/dx) / (2·η_op · k²)  <  1
```

which at 64² baseline reduces to `Sp < 0.6` — below the §4.8 band.

This is the **second refutation of a §4.x design-note prediction** in
this milestone. The first was Step 7 (slab-pull as amplifier, not
initiator). Both findings are revisions of implicit assumptions in
`solver-scaling.md`, not implementation bugs.

## Resolution paths (none selected; issue asks which)

### (a) Recalibrate Sp in the activated regime

Reset the §4.8 target band based on the activated-regime operator
balance. Straightforward algebra gives a new upper bound
`Sp < 2·η_op·k² · dx / (2·k_slab·τ_slab) = 0.6` at 64² baseline.

- **Pros:** minimal code change — just the §4.8 text and the default
  `Sp` constant.
- **Cons:** breaks the §4.8 band semantics. Sp now has a different
  meaning (no longer "natural slab-pull strength" but "maximum Sp
  stable on the 64² Newton operator"). Scale-dependent: refining to
  128² or 256² changes the bound by `dx`.

### (b) Modify the discrete divergence operator used in `Q_sub_conv`

The `1/dx` amplification is a discretisation choice, not a physical
requirement. A smoothed or gradient-bounded variant would reduce the
loop gain without altering the §4.8 `Sp` band:

- Average `|div v|` over a 3×3 stencil before applying `k_slab`.
- Use the L²-projection of `div v` onto a coarse grid.
- Cap `max(0, −div v)` at a physically-motivated ceiling tied to the
  thermal advection time scale.

- **Pros:** preserves §4.8 semantics. Grid-independent.
- **Cons:** modifies `slab/accumulation.rs` behaviour. Needs its own
  MMS + regression cycle. Picks a specific variant out of several
  plausible ones — additional design decision.

### (c) Physical saturation of `m_subducted`

Introduce an upper bound or nonlinear growth law on the `m` ODE so
`m_steady` does not scale linearly with `|div v|`. Options:

- `dm/dt = k · Q − m/τ − β · m²` (quadratic self-limit).
- `dm/dt = k · Q · (1 − m/m_max) − m/τ` (logistic saturation).
- Slab detachment events when `m > m_max`: remove the mass
  instantly, matching the geological picture of slab break-off.

- **Pros:** physically motivated. Matches the original §4.8
  bounded-growth intent (the τ exponential decay was already a step
  in this direction).
- **Cons:** most invasive. Changes Step 7's slab-pull contract.
  Requires its own MMS validation and re-running Step 7 physics to
  check the mantle-off regime still behaves.

## Decision requested

Which path? The reviewer's call. Once selected, a dedicated step
implementing the chosen resolution can start; Step 9 (cratonic
immunity) is not blocked by this issue and can run in parallel on
Step 8's mantle-only base.

## Deliverables of the resolution step

- Implementation of (a), (b), or (c).
- MMS tests for any new code path.
- Promote `v2_mantle_runaway_diagnostic` from `#[ignore]`-d to a
  non-ignored regression guard — the test currently reproduces the
  instability on baseline parameters; after resolution, the same test
  must pass (bounded `peak|v_solved|`) or be updated to a new
  stability bound consistent with the chosen path.
- Re-open the §Regression run convention to full scope (Step 8's
  exception is removed; regressions once again track "Step N physics
  − new-mechanism").
- Update `solver-scaling.md` §4.8 to state the resolved constraint.
- Update `tectonics_v2/README.md` to remove the Step 8 exception.

## Non-goals

- Time evolution of the mantle pattern (`evolution_rate > 0`) —
  orthogonal, deferred per D6.
- Coupling sweep — the `coupling` parameter was held at `1.0`
  throughout Step 8 and the instability scales with `Sp · k_slab`
  rather than with `coupling` at leading order.
- Thermal Arrhenius viscosity (§4.10 Q3) — would change `η_op` in
  the activated regime and in principle relaxes the stability bound,
  but introducing it to resolve this issue would conflate two
  unrelated pieces of physics.
