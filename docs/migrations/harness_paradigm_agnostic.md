# Issue #125 Phase 1.3 Stage H1 — Harness paradigm-agnostic refactor audit

**Status:** Phase 1.3 Stage H1 audit only. Read-only. No source modifications,
no `cfg` changes, no `git mv`. H2 / H3 land the implementation if the
recommendation below is accepted.

**Scope:** Cartograph the v2 harness and its workflow callers, identify the
load-bearing paradigm-coupling points, evaluate two design options (trait
abstraction vs two separate functions), and recommend one based on the
explicit decision criteria from the Phase 1.3 spec.

**Audit cut-off commit:** `ecf2921` (Phase 1.3 Stage E1.bis), branch
`125-c1-phase-13-equilibrium-height-closure-harness-paradigm-agnostic-refactor`,
off `3746da7` (Phase 1.2 merge into `milestone/c1-lightweight-dynamic-tectonics`).

**HC4 reference:** `docs/migrations/v2_to_c1_attic.md` § HC4 (lines 389-393)
explicitly defers the paradigm-agnostic workflow refactor to Phase 1.3+. This
is that refactor.

---

## Summary of findings

- **Harness is monolithically v2-coupled.** Single file `harness.rs`, 3 public
  functions, `BaselineConfig` embeds 6 v2-specific config types
  (`MantleConfig`, `SlabPullConfig`, `BasalDragConfig`, `LinearSolverConfig`,
  `YieldingConfig`, `RheologyConfig` via `Preset`). Whole file gated under
  `v2_legacy`.
- **The leak point is `CycleOutput.baseline: BaselineResult`** in
  `workflow/mod.rs:230`. This is the *single* type-system bridge between
  the otherwise paradigm-agnostic workflow scaffolding and the v2 harness.
  Everything else in `workflow/mod.rs` (`WorkflowConfig`, `WorkflowParams`,
  `PhaseAParams`, `PhaseBParams`, `PhaseBOutput`) is data-only and reusable.
- **Phase B is already paradigm-agnostic at runtime.** It consumes `Field2D
  + sea_level_normalized + iso_config`, not `BaselineResult` — no refactor
  needed for Phase B.
- **C1 has zero coupling to harness.** `tectonics_c1/time_loop.rs` is a
  self-contained per-cycle runner that already satisfies the HC4 statement
  "C1 provides its own per-cycle tectonic runner" — it just needs a
  workflow wrapper.
- **Structural divergence is material.** v2's `BaselineConfig` carries 23
  fields including 6 nested config types; C1's runner takes 4 args
  (`&mut C1State`, `&PlateKinematics`, `&C1TimeLoopConfig`, `&C1Closures`).
  Per-step payloads (`StepProgress<'_>` vs `&C1State`) are structurally
  different, not just type-renamed.
- **Surface area to preserve.** 63 bin/test entries in `crates/ymir-core/
  Cargo.toml` carry `required-features = ["v2_legacy"]`. Any signature
  change to the v2 path triples the regression risk.

**Recommendation:** **Option B — Two separate functions.** All three
decision criteria from the Phase 1.3 spec point the same way, and the
effort estimate (1-2 days) sits under the W3 trait fallback threshold
(3 days).

---

## §1 — Harness inventory

### §1.1 — Public surface of `tectonics_v2/diagnostics/harness.rs`

Single file. Module-level gated by `#[cfg(feature = "v2_legacy")]` in
`diagnostics/mod.rs:32-39`.

| Item | Signature | v2-coupling severity |
|---|---|---|
| `pub fn build_force` | `(kind: ForceKind, scales: &Scales, sin_amplitude: f64, domain_lx: f64) -> Box<dyn BodyForce>` | high — returns Stokes body-force trait |
| `pub fn run_baseline` | `(cfg: &BaselineConfig) -> BaselineResult` | high — only entry into the v2 pipeline |
| `pub fn run_baseline_with_progress<F>` | `(cfg: &BaselineConfig, mut on_progress: F) -> BaselineResult where F: FnMut(&StepProgress<'_>) -> bool` | high — adds early-stop callback |
| `enum ForceKind { Gpe, Sinusoidal }` | + `pub fn label()` | data |
| `enum NonlinearChoice { Newton, Picard }` | + `pub fn label()`, `parse()` | data, but Stokes-specific |
| `struct BaselineConfig` (23 fields) | embeds 6 v2 configs | **bridgehead** |
| `struct ContinuationState` | velocity warm-start | v2-Stokes-specific |
| `struct HarnessCaptureSpec` | benchmark snapshot spec | v2-bench-specific |
| `struct BaselineResult` | `{ metrics, config_dump, final_state }` | composite v2 result |
| `struct FinalState` | 5 v2-optional fields | composite v2 result |
| `struct StepProgress<'a>` | per-step callback payload | composite v2 progress |

### §1.2 — v2-specific types embedded in `BaselineConfig`

| Field | Type | Step |
|---|---|---|
| `mantle` | `MantleConfig` | Step 8 |
| `slab_pull` | `SlabPullConfig` | Step 7 |
| `basal_drag` | `BasalDragConfig` | Step 4 |
| `linear_solver` | `LinearSolverConfig` | Step 8.5a |
| `yielding` | `YieldingConfig` | Step 3 (via preset.rheology) |
| `preset` | `Preset` | carries rheology / scales |

These configs have no semantic analog in C1 (no Stokes solver, no mantle
coupling, no slab-pull body force). A paradigm-agnostic abstraction over
`BaselineConfig` would have to be opaque — workflow code could never
reach inside.

### §1.3 — v2-specific optional fields in `FinalState`

| Field | Type | Step |
|---|---|---|
| `age_field` | `Option<Field2D>` | Step 10 |
| `cratonic_factor` | `Option<Field2D>` | Step 9 |
| `plate_id` | `Option<PlateIdField>` | from boundary config |
| `plate_type` | `Option<PlateTypeField>` | from boundary config |
| `boundary_flag` | `Option<BoundaryFlagField>` | from boundary config |

C1's analog is `C1State { s, age, plate_id, plate_type, cratonic_mask, num_plates }`. Field shapes
overlap (S̃ field, age field, plate_id field, plate_type field) — but C1
has no `boundary_flag` field, and its `cratonic_mask` is `BoolField` while
v2's `cratonic_factor` is `Field2D`. The runtime per-cycle product is
*similar enough that a workflow consumer (Phase B erosion) doesn't care*,
but the typed envelope is different.

---

## §2 — Workflow callers

### §2.1 — `tectonics_v2/workflow/phase_a.rs` — the load-bearing call site

Currently the only piece of the workflow that runs the harness. Module
gated by `#[cfg(feature = "v2_legacy")]` in `workflow/mod.rs:42`. Public
entry points re-exported on lines 54-58:

- `run_phase_a_cycle(cfg: &BaselineConfig, wf: &WorkflowConfig) -> CycleOutput`
- `run_phase_a_cycle_with_progress<F>(cfg, wf, on_progress: F) -> CycleOutput`
- `run_phase_a_loop(cfg, wf) -> PhaseAOutput`
- `final_state_to_continuation(FinalState) -> ContinuationState`

All four take or return v2-typed values.

### §2.2 — `tectonics_v2/workflow/phase_b.rs` — already paradigm-agnostic at runtime

Currently gated for symmetry (`workflow/mod.rs:44`), but the **runtime
input it actually reads** is `Field2D + sea_level_normalized + iso_config
+ FbmUpscaleConfig + ErosionConfig`. Of these, only `Field2D` is shared
between v2 and C1 (both use `tectonics_v2::field::Field2D`, which the C1
state exposes via `state.s`). The other inputs are paradigm-agnostic
domain types.

`PhaseBOutput { heightmap, sediment, slope, deviation, deviation_p95 }`
is `GridF32`-based — no v2 types.

**Implication:** under either Option A or B, Phase B's wrapper can stay
verbatim; only its callers (Phase A loop) need re-pointing.

### §2.3 — `workflow/mod.rs` — the bridgehead struct

```text
// workflow/mod.rs:230
#[cfg(feature = "v2_legacy")]
pub struct CycleOutput {
    pub baseline: crate::tectonics_v2::diagnostics::harness::BaselineResult,
    pub erosion_volume_removed: f64,
    pub erosion_peak_delta_h: f64,
    pub sea_level_normalized: f64,
    pub mass_drift: f64,
    pub craton_recomputation_change: Option<f64>,
}
```

Only the `baseline: BaselineResult` field forces the gate. The other six
fields (erosion volume + peak Δh, sea level, mass drift, craton-change)
are paradigm-agnostic Phase A diagnostics that C1 can populate equally
well.

---

## §3 — C1 cross-paradigm coupling

**Status: none.** Confirmed by Grep against `crates/ymir-core/src/
tectonics_c1/`:

- No `use crate::tectonics_v2::diagnostics::harness` imports.
- No `use crate::tectonics_v2::workflow` imports.
- `tectonics_c1::time_loop` is the self-contained C1 per-cycle runner.
  Its public surface is:
  - `run_advection_only(state, kinematics, config, on_step)` — Phase 1.1
  - `run_with_closures(state, kinematics, config, closures, on_step)` —
    Phase 1.2 + 1.3 (now wires Davis-Suppe + equilibrium-height)
- `C1Closures { davis_suppe, equilibrium_height }` parallels the
  per-closure flags of `BaselineConfig.mantle/slab_pull/basal_drag`.

**Implication:** the C1 side of the abstraction already exists. H2 only
needs to wrap it in a workflow entry that produces the same paradigm-
agnostic scalars (`erosion_volume_removed`, `mass_drift`, …) as v2 does.

---

## §4 — `v2_legacy` feature surface area

Declared at `crates/ymir-core/Cargo.toml:32` (`v2_legacy = []`).

| Surface | Count | Notes |
|---|---|---|
| Module-level `#[cfg(feature = "v2_legacy")]` gates | 28 | `diagnostics/mod.rs` (9), `tectonics_v2/mod.rs` (2), `workflow/mod.rs` (6), `forcing/mod.rs` (2), `_attic/mod.rs` (7), inline (2) |
| Cargo.toml entries `required-features = ["v2_legacy"]` | 63 | bin + bench + test entries that don't compile without the flag |
| Existing `workflow` items gated transitively | 4 fns + 2 structs | `run_phase_a_*`, `run_phase_b`, `CycleOutput`, `PhaseAOutput` |

Any change that touches the v2 path signature multiplies risk by 63
(every bin/test entry).

---

## §5 — Structural comparison v2 vs C1

The decision criteria from the Phase 1.3 spec call for an honest
divergence table:

| Aspect | v2 harness | C1 time_loop | divergence |
|---|---|---|---|
| Config bundle | 23-field `BaselineConfig` with 6 nested v2-only configs | 4 args (`&mut C1State`, `&PlateKinematics`, `&C1TimeLoopConfig`, `&C1Closures`); no nested v2 configs | **structural** |
| Mutable state | Built inside `run_baseline` from `BaselineConfig` | Pre-built `C1State` passed in by caller | structural |
| Per-step progress | `StepProgress<'_>` — Stokes residuals, nonlinear iter count, mantle/slab residuals, age/cratonic snapshots | `usize, &C1State` — step index + simple state ref | **structural** |
| Result envelope | `BaselineResult { metrics: Metrics, config_dump: SolverConfigDump, final_state: FinalState }` | None — mutates `&mut C1State` in place | **structural** |
| Time stepping | Variable Δt from nonlinear continuation; multi-iteration per visible step | Fixed CFL Δt; single Forward-Euler step per visible step | structural |
| Mass balance | Source / sink balance through Stokes velocity + boundary fluxes | Pure advection of S̃ + per-cell additive closures | semantic |
| Early-stop | `on_progress -> bool` allows aborting mid-run | `on_step` has no return value | semantic |

The two paths share **no common signature modulo types**. They share
**no common return envelope**. They share **no common per-step payload
shape**. The only thing they share is the *concept* "advance a tectonic
state by N visible steps".

---

## §6 — Design options

### Option A — Trait abstraction

#### A.1 Sketch

```rust
pub trait TectonicRunner {
    type Config;
    type Progress<'a>;
    type Result;
    fn run<F>(cfg: &Self::Config, on_progress: F) -> Self::Result
    where
        F: FnMut(&Self::Progress<'_>) -> ControlFlow;
}

pub enum ControlFlow { Continue, Stop }

pub struct V2Runner;
impl TectonicRunner for V2Runner {
    type Config = BaselineConfig;
    type Progress<'a> = StepProgress<'a>;
    type Result = BaselineResult;
    fn run<F>(cfg, on_progress) -> Self::Result { /* delegate to harness::run_baseline_with_progress */ }
}

pub struct C1Runner;
impl TectonicRunner for C1Runner {
    type Config = (/* C1State + kinematics + config + closures bundle */);
    type Progress<'a> = (usize, &'a C1State);
    type Result = ();  // mutates in place
    fn run<F>(cfg, on_progress) -> Self::Result { /* delegate to time_loop::run_with_closures */ }
}

// workflow becomes generic
pub fn run_phase_a_cycle<R: TectonicRunner>(
    cfg: &R::Config,
    wf: &WorkflowConfig,
) -> CycleOutput<R::Result> { … }
```

#### A.2 Pros

- Single conceptual entry point: workflow code calls `R::run(…)`.
- Future paradigms (Phase 2.x evolution) plug in via new `impl
  TectonicRunner for FooRunner` — extensibility argument.
- Reduces lexical duplication in `workflow/` (one `run_phase_a_cycle`
  generic).

#### A.3 Cons

- **3 associated types** (`Config`, `Progress`, `Result`) → workflow
  becomes generic over `R`, which propagates through:
  - `CycleOutput<R::Result>` — generic
  - `PhaseAOutput<R::Result>` — generic
  - All 63 bin/test entries — every callsite must monomorphise the
    parameter
- The two runners' `Result`s are fundamentally different shapes
  (`BaselineResult` vs `()`). Workflow code that wants e.g. "the
  final S̃ field" must access it via *different* paths depending on
  the runner — abstraction *leaks* through the workflow regardless.
- C1's `Progress` is `(usize, &'a C1State)` — that's a tuple with a
  lifetime. v2's `Progress<'a>` is `StepProgress<'a>` — a struct
  with a lifetime. GAT (generic associated types with lifetimes)
  works since edition 2024, but error messages and lifetime
  inference become harder for downstream callers.
- `on_progress` semantics differ — v2's callback returns `bool`
  (early stop), C1's callback returns `()`. Reconciling forces a
  `ControlFlow` enum and a wrapper for the C1 callback that always
  returns `Continue`. Subtle correctness risk on misuse.
- The bin/test surface (63 entries) is touched by signature changes
  to `run_phase_a_*`. Even if behaviour-preserving, all 63 entries
  need re-checking under both `--features v2_legacy` and default.

#### A.4 Effort estimate H2

- Trait definition + GAT setup: 0.5 d
- `V2Runner` impl: 0.5 d
- `C1Runner` impl + the wrapping `(&mut C1State, …)` Config bundle: 0.5 d
- Refactor `workflow/phase_a.rs` to be generic: 1.5-2 d (the most
  fragile step — generic over `Result` makes `CycleOutput` generic,
  which cascades)
- Update 63 bin/test entries: 1-2 d (mostly mechanical, but several
  need explicit `::<V2Runner>` annotations)
- Coverage tests for both runners: 0.5 d
- Doc + migration notes: 0.5 d

**Total: 4-5 days.** Past the 3-day W3 fallback threshold.

### Option B — Two separate functions

#### B.1 Sketch

Module structure after H2:

```text
crates/ymir-core/src/tectonics_v2/workflow/
├── mod.rs              ← WorkflowConfig, WorkflowParams, PhaseAParams,
│                          PhaseBParams, CycleOutputCommon (paradigm-agnostic)
├── phase_a_v2.rs       ← (rename of current phase_a.rs)
│                          gated #[cfg(feature = "v2_legacy")]
│                          re-exports run_phase_a_cycle_v2, _with_progress, _loop
├── phase_a_c1.rs       ← NEW, default-features-on
│                          run_phase_a_cycle_c1(...) calls time_loop::run_with_closures
├── phase_b.rs          ← unchanged (already paradigm-agnostic at runtime;
│                          may need to drop the cfg gate)
├── drainage.rs         ← unchanged
└── macro_redistribution.rs  ← unchanged
```

Shared per-cycle output:

```rust
// workflow/mod.rs — promoted to paradigm-agnostic
pub struct CycleOutputCommon {
    pub erosion_volume_removed: f64,
    pub erosion_peak_delta_h: f64,
    pub sea_level_normalized: f64,
    pub mass_drift: f64,
    pub craton_recomputation_change: Option<f64>,
}

// workflow/phase_a_v2.rs — keeps the baseline field
#[cfg(feature = "v2_legacy")]
pub struct CycleOutputV2 {
    pub baseline: harness::BaselineResult,
    pub common: CycleOutputCommon,
}

// workflow/phase_a_c1.rs — paradigm-specific surface
pub struct CycleOutputC1 {
    pub final_state: C1State,        // moved out of run, owned
    pub common: CycleOutputCommon,
}
```

Phase A entry points after H2:

- v2 path: `run_phase_a_cycle_v2(cfg: &BaselineConfig, wf: &WorkflowConfig) -> CycleOutputV2` (gated).
- C1 path: `run_phase_a_cycle_c1(state, kinematics, c1_cfg, closures, wf) -> CycleOutputC1` (default).

Phase B path is paradigm-agnostic — same `phase_b::run_phase_b(s_field,
wf, seed)` works for both. The only change is dropping its `v2_legacy`
gate.

#### B.2 Pros

- **Honest about structural divergence.** Two paths because the
  paradigms actually are two different things. No forced
  abstraction that would let one leak into the other.
- **Zero risk to the 63 existing v2 bin/test entries.** A pure rename
  + cfg-gate keeps every existing call site working bit-identically
  under `--features v2_legacy`.
- **No GATs, no `ControlFlow` enum wrapping, no `R::Result`
  monomorphisation.** Each path stays plain Rust.
- **Each path owns its full type richness.** v2 keeps
  `BaselineResult` with its 23 fields and 5 optional final-state
  fields. C1 keeps its mutate-in-place `&mut C1State` semantics.
- **Phase B becomes default-features-on for free.** Just drop its
  cfg gate. PhaseBOutput was already paradigm-agnostic.
- **Lower bin/test surface area.** No bin/test currently exercises
  C1 through the workflow (C1 is too new). The new C1 path only
  needs ~2-3 smoke tests added in H3.

#### B.3 Cons

- Lexical duplication of orchestration logic between the two
  `phase_a_*.rs`. Mitigatable by extracting the common bits (erosion
  pass, isostasy call, craton recompute) into helpers in
  `workflow/phase_a_common.rs` — the per-cycle erosion + isostasy +
  craton logic is paradigm-agnostic by construction.
- Future Phase 2.x paradigm (if it materialises) requires a third
  `phase_a_p2.rs`. Acceptable — Phase 2 is far away and its
  paradigm shape is unknown.
- Two test surfaces to keep in step (v2 regression tests + new C1
  smoke tests). Mitigated by shared `CycleOutputCommon` invariants
  (e.g. `mass_drift` semantics).

#### B.4 Effort estimate H2

- `git mv phase_a.rs phase_a_v2.rs` + re-export gating in
  `workflow/mod.rs`: 0.25 d
- Extract `CycleOutputCommon` + `CycleOutputV2` split: 0.5 d
- Create `phase_a_c1.rs` parallel structure calling
  `time_loop::run_with_closures` + paradigm-agnostic post-pass
  (erosion + isostasy + craton update via existing helpers): 0.5-0.75 d
- Drop `v2_legacy` gate from `phase_b.rs` (it's already paradigm-
  agnostic at runtime): 0.1 d
- H3 smoke tests on the C1 path (default features): 0.5 d
- Doc + migration update (this file → "implemented" + brief change
  note in `v2_to_c1_attic.md` HC4 § "resolved"): 0.25 d

**Total: 1.5-2 days.** Under the W3 trait fallback threshold (3 days).

---

## §7 — Decision criteria applied

The Phase 1.3 spec lists three criteria. Applied to the audit findings:

| # | Criterion | Verdict | Reasoning |
|---|---|---|---|
| 1 | Same signature modulo types? | **NO** | Configs (23 nested-v2 fields vs 4 paradigm-agnostic args), per-step payloads (`StepProgress<'_>` vs `(usize, &C1State)`), result envelopes (`BaselineResult` vs `()`), and time-stepping semantics (variable continuation Δt vs fixed CFL Δt) all diverge structurally. |
| 2 | Configs / returns structurally different? | **YES** | See §5 divergence table. The signature differences are not type-renames; they reflect that v2 owns its state internally while C1 mutates an externally-owned `C1State`. |
| 3 | Trait would need type erasure or complex associated types? | **YES** | Option A requires 3 associated types incl. a GAT (`Progress<'a>`) and a `ControlFlow` enum to reconcile the callback-return mismatch. Not impossible, but explicitly the case the criterion warns against. |

All three criteria point the same way → **Option B is the indicated
choice** per the criterion-driven decision rule.

---

## §8 — Recommendation + open questions

### §8.1 Recommendation

**Option B — Two separate functions.**

Decisive factors, in order:

1. Effort estimate (1.5-2 d) sits clearly below the W3 trait-fallback
   threshold (3 d). Option A is 4-5 d.
2. All three criteria from the spec point to B.
3. Risk surface is materially lower: 63 existing bin/test entries
   stay untouched under `v2_legacy`, vs Option A which touches every
   single one of them via the generic parameter.
4. Honesty argument: the two paradigms ARE structurally different;
   forcing a trait abstraction misrepresents the difference and
   makes the abstraction itself a maintenance liability.

The future-paradigm extensibility argument (Option A's main pro) is
not load-bearing — Phase 2.x is far enough that we don't know what
shape it will take.

### §8.2 Risks identified

| ID | Risk | Mitigation |
|---|---|---|
| R1 | C1 path has no `with_progress` variant initially → workflow API asymmetry | Ship only `run_phase_a_cycle_c1(..., on_step: F)` from day one; C1's `time_loop::run_with_closures` already takes a per-step callback. Symmetry trivially restored. |
| R2 | `phase_b.rs` un-gating reveals a hidden `v2_legacy`-only call we missed | Verify in H2 first commit: `cargo build` without `v2_legacy` after the gate drop. Grep already confirmed only `Field2D + isostasy + upscale + erosion` consumed. |
| R3 | The shared `CycleOutputCommon` scalars (mass_drift, erosion_volume_removed) need to be **computed identically** in both paths or the diagnostics will diverge silently | Extract the post-tectonic-step erosion / isostasy / craton block into `workflow/phase_a_common.rs` and have both `phase_a_v2.rs` and `phase_a_c1.rs` call it. One implementation, two callers — no silent divergence. |
| R4 | Workflow tests live under `v2_legacy` today; H3 needs to move some of them to default-features so the C1 path has coverage | Smoke tests under `cfg(not(feature = "v2_legacy"))` referencing only the C1 entry points. The v2 regression test stays gated. |
| R5 | The `Disabled` regression contract (`workflow_disabled = run_baseline` bit-identical, currently asserted by `tests/v2_workflow_disabled_regression.rs`) is v2-specific. The C1 path needs its own equivalent: `phase_a_c1` with `WorkflowConfig::Disabled` must equal a direct `time_loop::run_with_closures` call. | Bit-identical regression test in H3, modelled on the v2 one. |
| R6 | If Option B's lexical duplication grows (e.g. Phase 2 adds a third path), pressure to retrofit a trait will appear | Document the criterion-driven decision in this file so the *next* refactor opens with explicit history rather than re-deriving. |

### §8.3 Open questions for user before GO H2

1. **Gating default for `phase_a_c1.rs`:** propose default-features-on (no
   feature flag). C1 is mainline as of Phase 1.3; v2 is the legacy path.
   Confirm or specify a `c1` feature if you want symmetry with
   `v2_legacy`.
2. **Naming:** `run_phase_a_cycle_v2` vs keeping `run_phase_a_cycle` (the
   current name) under `v2_legacy`. Renaming is more honest but touches
   the 63 bin/test entries. Recommended: **keep `run_phase_a_cycle` as
   the v2 name** (no rename), add `run_phase_a_cycle_c1` as the new
   default path. Zero churn in existing bin/test entries.
3. **`phase_b.rs` un-gating:** recommended to drop the cfg gate in H2.
   Already paradigm-agnostic at runtime per §2.2. Confirm.
4. **Migration doc update:** propose updating `v2_to_c1_attic.md` HC4
   section to "resolved by Phase 1.3 H2" with a back-link to this file.
   Confirm.

---

## §9 — Out of scope (not H1, not H2, not H3)

- Migrating individual workflow callers in `crates/ymir-viz/src/bridge/v2/`
  to the C1 path. Bridge layer stays v2-only for now (per Issue #117).
- Retiring the v2 harness itself. The plan is to keep it gated; Phase
  2.x decides whether to retire or keep as a regression anchor.
- Re-implementing v2's continuation warm-start in C1. C1 has no
  Stokes solver and no notion of velocity warm-start; the C1 cycle
  loop is stateless across cycles modulo the `C1State` it carries.
- Adding new C1 closures (Phase 1.4 erosion is a separate phase).

---

## §10 — Implementation roadmap (sketch only, H2 scope)

This section is descriptive only — the *audit* recommends Option B,
the *implementation* lands in H2. Listed here so the user can sanity-
check the proposal end-to-end.

1. **H2 Commit 1** (~0.5 d) — Carve out `CycleOutputCommon`. Move the
   paradigm-agnostic fields out of `CycleOutput`, repackage as
   `CycleOutputV2 { baseline, common }`. Run v2 regression tests under
   `v2_legacy`.
2. **H2 Commit 2** (~0.5-0.75 d) — `git mv phase_a.rs phase_a_v2.rs`.
   Update `workflow/mod.rs` re-exports. Run v2 regression tests.
3. **H2 Commit 3** (~0.5 d) — Drop `v2_legacy` gate from `phase_b.rs`.
   Verify `cargo build` (default features) passes. Smoke test on
   the GridF32 surface.
4. **H2 Commit 4** (~0.5 d) — Create `workflow/phase_a_common.rs`
   with the shared post-tectonic-step pass (isostasy + craton +
   erosion volume). Call from `phase_a_v2.rs`.
5. **H2 Commit 5** (~0.5 d) — Create `workflow/phase_a_c1.rs`,
   default-features-on, calls `time_loop::run_with_closures`,
   delegates the post-pass to `phase_a_common.rs`.
6. **H3** (~0.5 d) — Smoke tests for the C1 path under default
   features, including the Disabled-equals-direct-call regression
   parallel.
7. **Doc tail** (~0.25 d) — Update this file's § Status to
   "implemented", update `v2_to_c1_attic.md` HC4 to "resolved".

Total H2 + H3: **1.5-2 days**, matching the §6.B.4 estimate.

---

## §11 — Hidden-coupling check (W7-style)

Per the Phase 1.3 Stage H1 W7 watchpoint, surfaced findings beyond the
documented HC4 debt:

- **None.** The HC4 statement in `v2_to_c1_attic.md` lines 389-393 is
  the only paradigm-coupling debt and it is precisely the
  `CycleOutput.baseline: BaselineResult` leak point identified in
  §2.3. C1 has no existing coupling to harness (§3). Workflow Phase B
  is already paradigm-agnostic at runtime (§2.2). The 63 v2_legacy
  bin/test entries are scope-bounded by Cargo.toml and require no
  refactor beyond the H2 plan above.

No STOP-and-surface required. GO/STOP on H2 is the explicit user gate
that follows this audit.
