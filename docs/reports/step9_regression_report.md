# Step 9 — Regression report (`CratonicConfig::Disabled`)

> Companion to `step9_physics_report.md`. Validates that
> `CratonicConfig::Disabled` (the default for `BaselineConfig`)
> produces **bit-identical** output to the pre-Step-9 code path
> across all Step 0–8 regimes. This is acceptance #13 (Step 8
> regression) and #14 (Step 7 regression).

## Test setup

The regression validation rides on two existing integration
tests, augmented to confirm the cratonic state's structural
by-pass branch:

- `crates/ymir-core/tests/v2_step8_regression_smoke.rs`
- `crates/ymir-core/tests/v2_step7_regression_smoke.rs`

Both tests construct `BaselineConfig` with
`cratonic: CratonicConfig::Disabled` (the default added in
Phase 4+5+6 — see commit `cb454e5`) and assert the same numerical
output as the merged Step 7 / Step 8 baselines, byte-for-byte
across repeated runs.

## Structural by-pass — implementation evidence

Every eta-build site in the solver pipeline matches on
`Option<&CratonicState>`:

```text
match cratonic {
    None    => /* Step 0–8 hot path, no per-cell branch */,
    Some(_) => /* Step 9 cratonic-aware path */,
}
```

Sites:
- `rheology::build_eta_field` (cell-centre η field)
- `stokes::operator::TangentContext::from_strain_rate`
  (Newton tangent + Picard η)
- `stokes::nonlinear_solver::compute_residual` and
  `evaluate_residual_norm`
- `stokes::picard::compute_residual`
- `stokes::continuation::run_continuation`

In the `None` arm, the inner loop is byte-equal to the
pre-Step-9 implementation (one preceding `match cratonic` ahead
of the loop, no per-cell branch). In the `Some` arm, the per-cell
multipliers (`bi_multiplier`, `eta_multiplier`) are applied; this
arm is **only** entered when `BaselineConfig.cratonic =
CratonicConfig::Enabled(_)`.

`BaselineConfig::dynamic_accidented_defaults` and every existing
sweep / harness preset use `CratonicConfig::Disabled` — so the
default code path is unchanged from Step 8.

## `v2_step8_regression_smoke` results

```bash
cargo test --release -p ymir-core --test v2_step8_regression_smoke
```

```text
running 2 tests
test mantle_disabled_produces_no_step8_diagnostics ... ok
test disabled_runs_are_bit_deterministic ... ok

test result: ok. 2 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 46.06s
```

Both tests pass after the Phase 4+5+6 wiring (commit `cb454e5`)
and after the Phase 7b `B_factor` plumbing (commit `db58ec8`,
`e48c9a3`). The bit-deterministic property survives because:

1. `Disabled` config takes the structural by-pass at every eta-
   build site (no extra arithmetic in the inner loop).
2. The `cratonic_state: Option<CratonicState>` field is
   `None` whenever `cfg.cratonic = Disabled` OR
   `cfg.boundary = Disabled` (no Voronoï partition to compute
   the factor field from).
3. No `f64` operation in the new code path is ever evaluated
   on the regression run — the `match` is the only addition and
   it has no numerical effect.

## `v2_step7_regression_smoke` results

```bash
cargo test --release -p ymir-core --test v2_step7_regression_smoke
```

```text
running 2 tests
test slab_disabled_run_matches_newton_convergence_rate_with_full_step6_setup ... ok
test slab_disabled_produces_no_step7_diagnostics_and_parity_on_scalars ... ok

test result: ok. 2 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 22.29s
```

Same conclusion as Step 8: `Disabled` cratonic config preserves
the Step 7 behaviour byte-for-byte.

## `v2_step6_refactor_parity` (Step 5/6) results

```bash
cargo test --release -p ymir-core --test v2_step6_refactor_parity
```

```text
running 1 test
test step5_open_mode_parity_after_refactor ... ok

test result: ok. 1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out; finished in 11.86s
```

Step 5/6 regression (Open mode boundary, no cratonic) preserved.

## Library unit test suite

After all Step 9 wiring (Phases 1–7b), the full lib test suite
passes 251 / 251 — the cratonic-aware additions are 24 of those
tests; the rest is the Step 0–8 inheritance unchanged.

```bash
cargo test --release -p ymir-core --lib tectonics_v2::
```

```text
test result: ok. 251 passed; 0 failed; 0 ignored; 0 measured; 138 filtered out
```

## Acceptance summary

| # | Criterion | Target | Status |
|---|---|---|---|
| 13 | Step 8 regression bit-identical with `CratonicConfig::Disabled` | scalar parity to ε_mach × accumulation; CG iters mean ×1.00 | ✅ `v2_step8_regression_smoke` PASS |
| 14 | Step 7 regression preserved | identical to Step 7 baseline | ✅ `v2_step7_regression_smoke` PASS |

## Definition of done — regression scope

- [x] `BaselineConfig.cratonic` field added with `Disabled` default
- [x] All `BaselineConfig` literals across binaries + integration tests updated to add `cratonic: Disabled`
- [x] Every eta-build call site takes `Option<&CratonicState>` and uses the structural by-pass when `None`
- [x] `v2_step8_regression_smoke`: 2/2 PASS
- [x] `v2_step7_regression_smoke`: 2/2 PASS
- [x] `v2_step6_refactor_parity`: 1/1 PASS
- [x] Library unit tests: 251/251 PASS

The "no behavioural change for Step 0–8 callers" contract holds.
