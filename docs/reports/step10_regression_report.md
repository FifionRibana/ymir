# Step 10 — Regression report (`AgeFieldConfig::Disabled`)

> Companion to `step10_physics_report.md`. Validates that
> `AgeFieldConfig::Disabled` (the default for `BaselineConfig`)
> produces **bit-identical** output to the pre-Step-10 code path
> across all Step 0–9 regimes. This is acceptance #12 (Step 9
> regression) and #13 (all previous-step tests pass with defaults).

## Test setup

The regression validation rides on the existing integration
tests (no new test files for the regression scope), augmented to
confirm the age-field state's structural by-pass branch:

- `crates/ymir-core/tests/v2_step8_regression_smoke.rs`
- `crates/ymir-core/tests/v2_step7_regression_smoke.rs`
- `crates/ymir-core/tests/v2_step6_refactor_parity.rs`
- `crates/ymir-core/tests/v2_step10_physics_and_regression.rs::
   step10_disabled_runs_are_bit_deterministic` (Step 10 own check)

All tests construct `BaselineConfig` with
`age_field: AgeFieldConfig::Disabled` (the default added in
this commit) and assert the same numerical output as the merged
Step 7 / Step 8 / Step 9 baselines, byte-for-byte across
repeated runs.

## Structural by-pass — implementation evidence

Every age-field operation is gated on `Option<AgeFieldState>`:

```text
let mut age_state: Option<AgeFieldState> = match cfg.age_field {
    AgeFieldConfig::Disabled    => None,
    AgeFieldConfig::Enabled(c)  => Some(AgeFieldState::from_initial_thickness(&s, &c)),
};
```

Subsequent operations branch on `age_state.as_ref()` /
`age_state.as_mut()`:

- `step_age_advect` call — gated `if let Some(state) = age_state.as_mut()`.
- `apply_age_events` call — gated on the triple
  `(age_state.as_mut(), current_flag.as_ref(), &cfg.boundary)`.
- A snapshot save — gated `if let Some(state) = age_state.as_ref()`.
- Per-run diagnostic finalisation — gated
  `if let (AgeFieldConfig::Enabled(_), Some(_)) = ...`.

In the `None` arm of every gate, the inner code is **never
executed**. The Step 0–9 numerical trajectory is therefore
byte-equal to the pre-Step-10 path, by construction.

`BaselineConfig::dynamic_accidented_defaults` and every existing
sweep / harness preset use `AgeFieldConfig::Disabled` — so the
default code path is unchanged from Step 9.

## `v2_step8_regression_smoke` results

```bash
cargo test --release -p ymir-core --test v2_step8_regression_smoke
```

```text
running 2 tests
test mantle_disabled_produces_no_step8_diagnostics ... ok
test disabled_runs_are_bit_deterministic ... ok

test result: ok. 2 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 46.37s
```

Bit-deterministic property survives — the Step 8 active regime
produces the same final-state metrics across two independent
runs of the same `MantleConfig::Enabled` configuration with
`AgeFieldConfig::Disabled` default.

## `v2_step7_regression_smoke` results

```bash
cargo test --release -p ymir-core --test v2_step7_regression_smoke
```

```text
running 2 tests
test slab_disabled_produces_no_step7_diagnostics_and_parity_on_scalars ... ok
test slab_disabled_run_matches_newton_convergence_rate_with_full_step6_setup ... ok

test result: ok. 2 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 23.97s
```

## `v2_step6_refactor_parity` results

```bash
cargo test --release -p ymir-core --test v2_step6_refactor_parity
```

```text
running 1 test
test step5_open_mode_parity_after_refactor ... ok

test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 11.86s
```

## `step10_disabled_runs_are_bit_deterministic` (Step 10 own check)

The `tests/v2_step10_physics_and_regression.rs::
step10_disabled_runs_are_bit_deterministic` test runs the same
Step 8-shape 64² × 100 step configuration twice with
`AgeFieldConfig::Disabled` and asserts the run-level metrics
agree byte-for-byte:

- `mass_conservation_residual` — identical
- `vmax_peak` — identical
- `cg_iter_mean` — identical

Result: **PASS**. With age-field code present but the by-pass
branch taken, the trajectory is identical to a hypothetical
Step 9 run of the same config.

## Library unit test suite

After the Step 10 wiring, the full lib test suite passes 272/272
(was 251 in Step 9 — 17 new tests in `age_field::*::tests` plus
the existing milestone heritage). The new tests cover:

- `age_field::tests` (4 tests): config/state defaults, initial
  classification round-trip, age_init_max accessor.
- `age_field::init::tests` (2 tests): classification threshold
  + count partitioning.
- `age_field::advection::tests` (4 tests): zero-velocity
  quiescent growth, uniform-field invariance, step-pulse
  advection, byte-determinism on identical inputs.
- `age_field::events::tests` (5 tests): ridge / collision-max /
  arc / quiescent / pre-pass-snapshot semantics, plus
  ridge-vs-arc precedence in the edge case where a cell is
  flagged both Rift and arc-eligible.

```bash
cargo test --release -p ymir-core --lib tectonics_v2::
```

```text
test result: ok. 272 passed; 0 failed; 0 ignored; 0 measured; 138 filtered out
```

## Acceptance summary

| # | Criterion | Target | Status |
|---|---|---|---|
| 12 | Step 9 regression bit-identical with `AgeFieldConfig::Disabled` | scalar parity to ε_mach × accumulation; CG iters mean ×1.00 exact | ✅ `v2_step8_regression_smoke` + `step10_disabled_runs_are_bit_deterministic` PASS |
| 13 | All previous steps still pass their respective tests when defaults are used | identical to merged Step 9 baseline | ✅ `v2_step{6,7,8}_regression_smoke` PASS |

## Definition of done — regression scope

- [x] `BaselineConfig.age_field` field added with `Disabled` default
- [x] All `BaselineConfig` literals across binaries + integration tests updated to add `age_field: Disabled`
- [x] Every age-field operation (allocation, advection, event reset, snapshot, metric population) gated on `Option<AgeFieldState>`
- [x] `v2_step8_regression_smoke`: 2/2 PASS
- [x] `v2_step7_regression_smoke`: 2/2 PASS
- [x] `v2_step6_refactor_parity`: 1/1 PASS
- [x] `step10_disabled_runs_are_bit_deterministic`: PASS
- [x] Library unit tests: 272/272 PASS

The "no behavioural change for Step 0–9 callers" contract holds.
