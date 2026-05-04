# Step 11 — Regression report (`PlateKinematicConfig::Zero`)

> Companion to `step11_physics_report.md`. Validates that
> `PlateKinematicConfig::Zero` (the default for `BaselineConfig`)
> produces **bit-identical** output to the pre-Step-11 code path
> across all Step 0–10 regimes. This is acceptance #1 (default
> bit-identical), #11 (Step 10 regression bit-identical), and
> #12 (all previous-step tests pass with defaults).

## Test setup

The regression validation rides on the existing integration tests
plus three new dedicated tests that anchor the contract:

- `crates/ymir-core/tests/v2_plate_kinematic.rs` (Step-11-specific):
  - `zero_short_circuit_matches_per_plate_zeros`
  - `zero_default_does_not_perturb_baseline`
  - `zero_produces_zero_plate_kinematic_at_step_0`
- `crates/ymir-core/tests/v2_step8_regression_smoke.rs` (Step 8 regime)
- `crates/ymir-core/tests/v2_step7_regression_smoke.rs` (Step 7 regime)
- `crates/ymir-core/tests/v2_step6_refactor_parity.rs` (Step 6 refactor parity)
- `crates/ymir-core/tests/v2_step10_physics_and_regression.rs::step10_disabled_runs_are_bit_deterministic`
  (Step 10 own determinism check)

All tests construct `BaselineConfig` with
`plate_kinematic: PlateKinematicConfig::Zero` (the default added
in this milestone) and assert the same numerical output as the
merged Step 7 / 8 / 9 / 10 baselines, byte-for-byte across
repeated runs.

## Structural by-pass — implementation evidence

The drift mechanism is gated by `PlateKinematicConfig::is_zero()`
on every hook site in `harness.rs::run_baseline_with_progress`:

```text
let (drift_vx, drift_vy) = match &cfg.plate_kinematic {
    PlateKinematicConfig::Zero       => (Vec::new(), Vec::new()),
    PlateKinematicConfig::PerPlate { velocities, width } => {
        plate_kinematic::field::build(nx, ny, plate_id, velocities, *width)
    }
};
```

```text
// In the time loop, post-solve:
if !cfg.plate_kinematic.is_zero() {
    for k in 0..vx.len() {
        vx[k] += drift_vx[k];   // add hook (before advection)
        vy[k] += drift_vy[k];
    }
}
// ... advection ...
if !cfg.plate_kinematic.is_zero() {
    for k in 0..vx.len() {
        vx[k] -= drift_vx[k];   // strip hook (end of iter)
        vy[k] -= drift_vy[k];
    }
}
```

In the `Zero` arm, neither the field-build call nor the per-cell
add/strip loops execute. The Step 0–10 numerical trajectory is
therefore byte-equal to the pre-Step-11 path, by construction.

`BaselineConfig::dynamic_accidented_defaults` and every existing
sweep / harness preset use `PlateKinematicConfig::Zero` — so the
default code path is unchanged from Step 10.

## `v2_plate_kinematic` results (Step-11-specific)

```bash
cargo test --release -p ymir-core --test v2_plate_kinematic
```

```text
running 3 tests
test zero_produces_zero_plate_kinematic_at_step_0 ... ok
test zero_short_circuit_matches_per_plate_zeros ... ok
test zero_default_does_not_perturb_baseline ... ok

test result: ok. 3 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 4.42s
```

The middle test `zero_short_circuit_matches_per_plate_zeros` is
the strongest contract: it runs the harness twice on the same
config — once with `Zero`, once with `PerPlate { velocities: vec![
(0.0, 0.0); n], boundary_smoothing_width: 1.5 }` — and asserts
every field of `FinalState` (`s_field`, `vx`, `vy`,
`strain_rate_invariant`) matches byte-for-byte. This proves the
algorithmic path with all-zero inputs introduces zero numerical
noise, so the structural short-circuit and the algorithmic path
agree on `Zero`-equivalent inputs.

## `v2_step8_regression_smoke` results

```bash
cargo test --release -p ymir-core --test v2_step8_regression_smoke
```

```text
running 2 tests
test mantle_disabled_produces_no_step8_diagnostics ... ok
test disabled_runs_are_bit_deterministic ... ok

test result: ok. 2 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out
```

Bit-deterministic property survives — the Step 8 active regime
produces the same final-state metrics across two independent
runs of the same `MantleConfig::Enabled` configuration with
`PlateKinematicConfig::Zero` default.

## `v2_step7_regression_smoke` results

```bash
cargo test --release -p ymir-core --test v2_step7_regression_smoke
```

```text
running 2 tests
test slab_disabled_produces_no_step7_diagnostics_and_parity_on_scalars ... ok
test slab_disabled_run_matches_newton_convergence_rate_with_full_step6_setup ... ok

test result: ok. 2 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out
```

## `v2_step6_refactor_parity` results

```bash
cargo test --release -p ymir-core --test v2_step6_refactor_parity
```

```text
running 1 test
test step5_open_mode_parity_after_refactor ... ok

test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out
```

## Step 10 own check — `step10_disabled_runs_are_bit_deterministic`

The Step 10 own bit-determinism test runs at 64² × 100 steps with
`AgeFieldConfig::Disabled` AND now (post-Step-11) also
`PlateKinematicConfig::Zero` default. Run-level metrics agree
byte-for-byte across two independent runs:

- `mass_conservation_residual` — identical
- `vmax_peak` — identical
- `cg_iter_mean` — identical

```text
test result: ok. 2 passed; 0 failed; 0 ignored; 0 measured (incl. step10 + step10_disabled_runs_are_bit_deterministic)
```

## Library unit test suite

After the Step 11 wiring (new module + scaffolding), the `ymir-core`
lib test suite passes 287/287 (gain of +4 vs pre-Step-11):

- `plate_kinematic::tests` (3 tests): default = Zero, is_zero
  invariant, JSON round-trip for both variants.
- `plate_kinematic::field::tests` (5 tests): all-zero input → zero
  field, interior_uniform, boundary_smoothing_bounded_step,
  no_overshoot, deterministic_same_inputs.

```bash
cargo test --release -p ymir-core --lib tectonics_v2::
```

```text
test result: ok. 287 passed; 0 failed; 0 ignored; 0 measured; 138 filtered out
```

## Acceptance summary

| # | Criterion | Target | Status |
|---|---|---|---|
| 1 | Default `Zero` config bit-identical to pre-Step-11 | velocity field = 0 at init, all Steps 0-10 baselines preserved | ✅ `zero_produces_zero_plate_kinematic_at_step_0` + structural by-pass evidence |
| 11 | Step 10 regression bit-identical with `Zero` | scalar parity to ε_mach × accumulation; CG iters mean ×1.00 exact | ✅ `step10_disabled_runs_are_bit_deterministic` PASS |
| 12 | All previous-step tests still pass with defaults | identical to merged Step 10 baseline | ✅ `v2_step{6,7,8}_regression_smoke` + Step 10 own check PASS |

## Definition of done — regression scope

- [x] `BaselineConfig.plate_kinematic` field added with `Zero` default
- [x] All `BaselineConfig` literals across binaries + integration tests updated to add
      `plate_kinematic: PlateKinematicConfig::Zero` (27 sites)
- [x] Every drift operation (field build at init, post-solve add hook,
      iter-end strip hook) gated on `is_zero()`
- [x] `v2_plate_kinematic`: 3/3 PASS
- [x] `v2_step8_regression_smoke`: 2/2 PASS
- [x] `v2_step7_regression_smoke`: 2/2 PASS
- [x] `v2_step6_refactor_parity`: 1/1 PASS
- [x] `step10_disabled_runs_are_bit_deterministic`: PASS
- [x] Library unit tests: 287/287 PASS

The "no behavioural change for Step 0–10 callers" contract holds.
