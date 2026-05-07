# Step 13.5 — Regression report

> Companion to `step13_5_physics_report.md`. Validates that:
>
> 1. The default `apply_fbm_to_oceanic = false` produces
>    **bit-identical** Step 13 output across all paths
>    (algorithm, serde, UI). Acceptance #1 / #13.
> 2. Steps 0–13 with their explicit init modes keep producing
>    the same numerical output. Acceptance #14.
> 3. The Phase 5 default amendment
>    (`FBM_AMPLITUDE_OCEANIC_DEFAULT 0.10 → 0.15`) only
>    affects the new opt-in `apply_fbm_to_oceanic = true`
>    path. The disabled path is byte-identical to Step 13.

## Phase 1 structural short-circuit

The Step 13.5 contract rests on one simple structural fact in
`crates/ymir-core/src/tectonics_v2/init/radial_profile_fbm.rs::
build`:

```text
// Step 13.5 — oceanic FBM extension. Strictly opt-in. When the
// flag is false the entire block below is skipped — no second
// `Fbm<Perlin>` instance, no oceanic cell modification — so the
// function output is bit-identical to its Step 13 form.
if apply_fbm_to_oceanic {
    // ... oceanic FBM block ...
}
```

When `apply_fbm_to_oceanic = false`, the function runs exactly
the same path as Step 13: build the radial baseline, apply
continental FBM, return. No second `Fbm<Perlin>` is constructed
(so the noise crate's RNG state is unaffected), no oceanic cell
is touched.

The `oceanic_fbm_disabled_preserves_step13` test is the
dedicated oracle: it builds with **bogus** values for
`fbm_amplitude_oceanic = 0.42`, `fbm_scale_oceanic = Some(0.07)`,
`fbm_seed_oceanic = Some(0xDEAD)` (values the user could plausibly
have left in a config from an earlier exploration) but flag =
`false`, and asserts byte-identical output to a build with the
disabled defaults. The flag must short-circuit before any
oceanic param is read.

## Step 13 acceptance — bit-identical regression

Step 13's acceptance suite stays green with the default Step 13.5
flag (`false`):

| Suite | Pre-Step-13.5 | Post-Step-13.5 |
|---|---|---|
| `tectonics_v2::init::radial_profile_fbm::tests` | 6/6 (Step 13) | **12/12** (6 Step 13 + 6 Step 13.5 — Step 13's 6 unchanged) |
| `tectonics_v2::init::tests::determinism_*` | 7/7 | **7/7** (loop extended with both disabled and enabled new variants) |
| `crates/ymir-core/tests/v2_step13_acceptance.rs` | 4 + 1 ignored | **4 + 1 ignored** ✓ |
| `tests/v2_step13_cg_ratio.rs` (`#[ignore]`) | passes 0.951× / 0.974× | unchanged numerically |
| `tests/v2_step13_visual_checkpoint.rs` (`#[ignore]`) | green | green (existing patchworks unchanged at the disabled-flag default) |

The Step 13 test files (`v2_step13_acceptance.rs`,
`v2_step13_cg_ratio.rs`, `v2_step13_visual_checkpoint.rs`)
construct `InitMode::RadialProfileWithFBM { … }` literals; the
five constructions there now include the four new oceanic FBM
fields with disabled-default values. Bit-identical with the
pre-Step-13.5 form by structural short-circuit.

## Backward-compatible serde

The four new fields on `InitMode::RadialProfileWithFBM` carry
`#[serde(default)]` (custom default fn for
`fbm_amplitude_oceanic` so it falls through to
`FBM_AMPLITUDE_OCEANIC_DEFAULT` rather than `0.0`). Same on the
viz-side `V2InitModeSpec::RadialProfileWithFBM`. Legacy preset
JSON written before Step 13.5 deserialises unchanged.

| Test | Validates |
|---|---|
| `bridge::v2::spec::tests::v2_panel_radial_fbm_legacy_preset_load` (Step 13.5) | Step-13-shape JSON (no oceanic keys) parses to disabled defaults; `into_core()` threads the disabled flag through |
| `v2_panel_radial_fbm_with_oceanic_roundtrip` (Step 13.5) | Three explicit-oceanic cases roundtrip serialise → deserialise byte-for-byte; serialised JSON contains the documented keys when the flag is on |
| `init_mode_spec_roundtrips_through_json` (extended) | Both disabled-default and enabled-with-`Some` cases roundtrip |
| `v2_panel_radial_modes_serde_roundtrip` (Step 13) | Hand-written legacy JSON now also implicitly tests the disabled-default load (Step 13.5 fields default-applied) |
| `run_spec_roundtrips_and_old_json_loads_with_uniform_default` (Step 8.6) | Pre-Step-13.5 `V2RunSpec` defaults to Uniform on the legacy preset path (unchanged) |
| `old_preset_without_plate_kinematic_defaults_to_zero` (Step 11) | Pre-Step-11 preset defaults to Zero (unchanged) |

7/7 viz spec tests pass.

## Steps 0–12 regression

Step 13.5 adds **only** four new fields to one variant of the
`InitMode` enum and a guarded if-block in `radial_profile_fbm::
build`. **No production code path** reaches the oceanic FBM
machinery unless the caller explicitly sets
`apply_fbm_to_oceanic = true`. Existing modes (Checkerboard,
Uniform, Gaussian, Convolution, RadialProfile, RadialProfile
WithFBM with the flag off) are unaffected. Steps 0–12 are
structurally insulated from Step 13.5 by construction.

Verification:

```text
cargo test --release -p ymir-core
  → 448 non-ignored tests pass (same baseline as before
    Step 13.5; Step 13.5 adds 6 new lib tests + 2 new
    integration tests, all green).
  → 1 pre-existing failure unchanged
    (tectonics::export::deserialize_legacy_metadata_without_upscale,
    legacy module unrelated to tectonics_v2 — confirmed
    pre-existing back to the Step-11 merge HEAD).
```

The Step 9 / Step 10 / Step 11 heavy regression integration
tests are `#[ignore]`'d as before; running with `--ignored`
produces the same numerical output as the pre-Step-13.5
baseline since the init modes there are explicitly set to
`Checkerboard` / `Uniform`, neither of which touches the new
fields.

## Phase 5 default amendment — oceanic-only impact

The Phase 5 commit changes one compile-time default:

- `FBM_AMPLITUDE_OCEANIC_DEFAULT`: `0.10` → `0.15`

This default is consumed by:

- The `V2InitModeSpec::radial_profile_fbm_default()` helper —
  the panel's "first time the user picks the variant" defaults.
  Surfaces only when the user hasn't yet flipped the
  `apply_fbm_to_oceanic` toggle; until they do, the value is
  unused.
- The `default_fbm_amplitude_oceanic()` `#[serde(default)]`
  helper for legacy preset JSON deserialisation. Same
  observation: surfaces only when the flag is then enabled
  by the user, otherwise unused.
- The Step 13.5 acceptance / visual-checkpoint /
  calibration-probe tests that reference the constant
  explicitly.

**No production code path** reaches the new default from any
existing test or preset:

- `Checkerboard`, `Uniform`, `Gaussian`, `Convolution`,
  `RadialProfile` (no FBM at all), `RadialProfileWithFBM`
  with `apply_fbm_to_oceanic = false`: ignore the constant
  entirely.
- Existing preset JSON files predate Step 13.5 and
  deserialise to `apply_fbm_to_oceanic = false` (default).
- Existing tests that use `RadialProfileWithFBM` (Step 13's
  acceptance / cg_ratio / visual checkpoint) explicitly pass
  the disabled-flag defaults at the construction site — same
  output before and after the constant change.

Verification: re-running the full ymir-core test suite after
the Phase 5 commit produces the same 448/449 pattern as before
(1 unrelated pre-existing failure unchanged).

## Step 13.5 own determinism contract

Same-input-same-output is verified by:

- `init::tests::determinism_same_seed_same_output` — extended
  with two `RadialProfileWithFBM` variants in the loop:
  disabled-defaults and enabled with the `apply_fbm_to_oceanic
  = true` path. Both produce byte-identical output across
  repeated runs.
- Step 13.5 acceptance test #5 (`oceanic_fbm_seed_default_
  derivation`) — `fbm_seed_oceanic = None` derives the seed
  from `fbm_seed XOR FBM_SEED_OCEANIC_XOR_MAGIC` byte-for-byte.
  Confirms the derivation is deterministic and matches the
  documented formula.
- Step 13.5 acceptance test #4 (`oceanic_fbm_seed_independence`)
  — distinct `fbm_seed_oceanic` values produce distinct
  oceanic fields while continental cells stay byte-identical
  (insulation between the two `Fbm<Perlin>` instances).

## Summary

| Acceptance | Test source | Status |
|---|---|---|
| #11 — mass conservation | governed by Step 8 contract; init-only mechanism | ✓ |
| #13 — Step 13 regression bit-identical with disabled flag | Phase 1 short-circuit + 12/12 lib tests + 4 Step 13 acceptance pass | ✓ |
| #14 — all previous step tests pass | 448 non-ignored, 1 unrelated pre-existing failure | ✓ |
| Backward-compat serde | 7/7 viz spec tests + 2 explicit Step 13.5 tests | ✓ |
| Default amendment scope | new default consumed only by Step 13.5 opt-in path | ✓ |

Phase 1 short-circuit + Phase 2 backward-compat serde + Phase 3
panel ellipsis + Phase 5 amplitude amendment land cleanly with
no regression in any of the 14 milestone steps.
