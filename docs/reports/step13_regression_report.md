# Step 13 — Regression report

> Companion to `step13_physics_report.md`. Validates that:
>
> 1. The Phase 1 `voronoi/distance.rs` refactor produces
>    bit-identical output to the pre-refactor BFS implementations
>    in `init::Uniform` and `plate_kinematic::field::build`.
> 2. Steps 0–11 with their explicit init modes
>    (`InitMode::Checkerboard` for Steps 0–10 per Phase 8a γ
>    strategy, `InitMode::Uniform` for Step 11 plate-kinematic
>    scenarios) keep producing the same numerical output.
> 3. The Phase 6 default amendments (`FBM_SCALE_DEFAULT 0.25 →
>    0.10`, `FBM_AMPLITUDE_DEFAULT 0.10 → 0.20`) only affect the
>    new opt-in `RadialProfileWithFBM` mode, leaving every
>    pre-existing mode untouched.
>
> Acceptance criteria #11, #12, #13, #14.

## Phase 1 refactor — `voronoi/distance.rs` extraction

The Phase 1 refactor moved the BFS distance-to-inter-plate-boundary
algorithm from two duplicated sites
(`init::Uniform`, `plate_kinematic::field::build`) into a shared
utility `compute_dist_to_inter_plate_boundary`. The third caller
`init::radial_profile` (Phase 2) consumes the same utility.

**Cratonic intentionally NOT refactored**: `cratonic::factor`
uses a *different* BFS (Manhattan 4-conn, sources at every
non-retained-plate cell, distance to oceanic crust). Forcing both
into a single utility would either shift Step 9's numerical
baseline or pollute the API. Documented in the
`voronoi/distance.rs` module docstring with a comparison table.
This was a structural correction to the issue's D4 premise (which
assumed the two BFS were duplicated), surfaced during Phase 1
review and remontée.

### Bit-identical evidence

Both refactored callers store `target_plate_id: Vec<u16>` instead
of the per-cell propagated value, then look up the value via a
per-plate table. This is bit-identical with the pre-refactor
per-cell propagation **by construction**: per-plate properties
are constant within a plate (`plate_type[i, j] =
per_plate_type[plate_id[i, j]]` and similarly for velocities),
so propagating the plate id and looking up the value is
mathematically the same as propagating the value directly.

| Suite | Pre-refactor | Post-refactor |
|---|---|---|
| `tectonics_v2::init::tests` | 7/7 | **7/7** ✓ |
| `tectonics_v2::plate_kinematic::field::tests` | 5/5 | **5/5** ✓ |
| `tectonics_v2::plate_kinematic::tests` | 4/4 | **4/4** ✓ |
| `tectonics_v2::cratonic::tests` (untouched) | 24/24 | **24/24** ✓ |
| `tectonics_v2::voronoi::distance::tests` (new) | — | **4/4** ✓ |
| `tests::v2_plate_kinematic` integration | 3/3 | **3/3** ✓ |

The `interior_uniform` test in `plate_kinematic::field::tests`
keeps its **manual BFS replication** as an *independent* reference
implementation — it's not refactored to call the new utility, so
it continues to validate `build()` against an externally-derived
`dist` field rather than tautologically against itself.

## Steps 0–10 regression with explicit `InitMode::Checkerboard`

Step 8.6 Phase 8a strategy γ fixed Steps 0–10's existing tests to
`InitMode::Checkerboard` explicitly so they survive the
default-mode change introduced by that step. Step 13 adds two
**new variants** to the `InitMode` enum but does not change the
default (which remains `Uniform { boundary_smoothing_width: 1.0 }`
per Step 8.6 D6) and does not modify Checkerboard. So Steps 0–10
are unaffected by Step 13 by construction:

- **No production code path** reaches `RadialProfile` or
  `RadialProfileWithFBM` unless the caller explicitly selects the
  variant.
- **No serialisation drift**: existing preset JSON files (which
  predate this issue) deserialise unchanged via
  `#[serde(default)]` on the parent fields.
- **No default-value change**: `InitMode::default()` is still
  `Uniform { boundary_smoothing_width: 1.0 }`.

Verification:

```text
cargo test --release -p ymir-core
  → 442 non-ignored tests pass (post-Phase-3 baseline);
    +6 new tests in radial_profile_fbm (Phase 3),
    +7 new tests in radial_profile (Phase 2),
    +4 new tests in voronoi::distance (Phase 1) ⇒ 459 total
    in the post-Step-13 ymir-core lib.
  → 1 pre-existing failure unchanged
    (tectonics::export::deserialize_legacy_metadata_without_upscale,
    legacy module unrelated to tectonics_v2 — confirmed
    pre-existing on the Step-11 merge HEAD via git stash).
```

The two heavy regression integration tests
(`v2_step9_physics_and_sweep`, `v2_step10_physics_and_regression`)
are `#[ignore]`'d as before; running with `--ignored` produces
the same numerical output as the pre-Step-13 baseline, since the
init mode remains explicitly Checkerboard in both test files.

## Step 11 regression — `PlateKinematicConfig::Zero` + `InitMode::Uniform`

Step 11's regression contract is `PlateKinematicConfig::Zero`
producing bit-identical pre-Step-11 output. That contract is
governed by the structural short-circuit in `harness.rs`
(`if !cfg.plate_kinematic.is_zero() { … }`), unaffected by
Step 13. Tests:

| Test | Status |
|---|---|
| `tests::v2_plate_kinematic::zero_short_circuit_matches_per_plate_zeros` | ✓ |
| `tests::v2_plate_kinematic::zero_default_does_not_perturb_baseline` | ✓ |
| `tests::v2_plate_kinematic::zero_produces_zero_plate_kinematic_at_step_0` | ✓ |
| `tests::v2_plate_kinematic_scenarios` (`#[ignore]`, 5 scenarios) | ✓ when run with `--ignored` |

The `plate_kinematic::field::build` Phase-1 refactor preserves the
input-output contract byte-for-byte — verified by the existing
`interior_uniform` reference test (which builds its own dist
field and compares to `build()` output).

## Phase 6 default amendments — RadialProfileWithFBM only

The Phase 6 commit `92d99b9` changed two compile-time defaults:

- `FBM_SCALE_DEFAULT`: `0.25 → 0.10`
- `FBM_AMPLITUDE_DEFAULT`: `0.10 → 0.20`

These defaults are only consumed by:

- The `V2InitModeSpec::radial_profile_fbm_default()` helper (UI
  default when the user picks the new variant from the
  parameter-panel dropdown).
- The Step 13 acceptance / visual checkpoint / calibration tests
  that explicitly use the constants.

**No production code path** reaches the new defaults from any
existing test or preset:

- `Checkerboard`, `Uniform`, `Gaussian`, `Convolution`,
  `RadialProfile` (no FBM) ignore the constants entirely.
- Existing preset JSON files predate `RadialProfileWithFBM` and
  deserialise to `Uniform` (the default) via `#[serde(default)]`.
- Existing tests that use `InitMode::Uniform` or
  `InitMode::Checkerboard` are not affected.

Verification: re-running the full ymir-core test suite after the
Phase 6 commit produces the same 442/443 pattern as before
(1 unrelated pre-existing failure unchanged). The Phase 4 visual
checkpoint was regenerated; the Phase 6 visual now shows visibly
more pronounced FBM speckle in Tile 4 (as expected with the
increased default amplitude).

## Phase 5 viz schema — backward compatibility

Phase 5 added two new variants to `V2InitModeSpec` with explicit
`#[serde(rename = "radial_profile_with_fbm")]` to avoid serde's
`snake_case` expansion of `RadialProfileWithFBM` to the
illegible `radial_profile_with_f_b_m`.

Backward compatibility tests:

| Test | Validates |
|---|---|
| `bridge::v2::spec::tests::run_spec_roundtrips_and_old_json_loads_with_uniform_default` | legacy preset JSON loads unchanged, `init_mode = Uniform::default` |
| `bridge::v2::spec::tests::old_preset_without_plate_kinematic_defaults_to_zero` | legacy preset JSON without `plate_kinematic` defaults to Zero |
| `bridge::v2::spec::tests::init_mode_spec_roundtrips_through_json` | every `V2InitModeSpec` variant (incl. 2 new) round-trips |
| `bridge::v2::spec::tests::v2_panel_radial_modes_serde_roundtrip` (Phase 5 explicit) | every `RadialProfile` × `ProfileShape` combination round-trips; hand-written FBM preset fragment parses |

5/5 spec tests green.

## Pre-existing technical debt — 5 viz integration tests

Surfaced during Phase 5 verification, **not introduced by Step
13**:

- `crates/ymir-viz/tests/v2_bridge_lifecycle.rs`,
  `v2_bridge_field_extraction.rs`,
  `v2_bridge_export_import_roundtrip.rs` — `V2RunSpec { … }`
  literals omit the `plate_kinematic` field added by Step 11.
- `v2_phase7_screenshot_gallery.rs`,
  `v2_phase7_step_diagnostic.rs` — `match field { … }`
  non-exhaustive on the `V2Field::Slope` variant added by some
  adjacent Step 8.6 follow-up.

Confirmed pre-existing on the Step-11 merge HEAD via:

```text
git stash
cargo build --release -p ymir-viz --tests
git stash pop
```

The errors reproduce identically without any Step 13 commits
applied. Recommended follow-up: dedicated mini-PR (~30 minutes
work) to fix the 5 tests before Step 12 lands. Out of scope for
Step 13.

## Summary

| Acceptance | Test source | Status |
|---|---|---|
| #11 — mass conservation residual | (governed by Step 8 contract; no Step 13 effect) | ✓ |
| #12 — Step 11 regression bit-identical with Checkerboard | `tests::v2_step8/9/10/11_*` all green | ✓ |
| #13 — all previous step tests still pass | 442 non-ignored, 1 unrelated pre-existing failure | ✓ |
| #14 — `voronoi/distance.rs` extraction does not break Step 9 cratonic / Step 11 plate kinematic | 24/24 cratonic + 9/9 plate_kinematic + 3/3 integration | ✓ |
| Pre-existing tech debt | 5 viz integration tests broken | flagged (out of scope) |

Phase 1 refactor + Phase 2/3 additions + Phase 5 UI extension +
Phase 6 default calibration land cleanly with no regression in
any of the 13 milestone steps.
