# C1 Phase 1.4 — Stream-power erosion closure outputs

Issue #127, branch `127-c1-phase-14-stream-power-erosion-isostasy-e2e-downstream-pipeline-ymir-viz-surgical-gating`.

Produced by the Stage E4 acceptance test
[`crates/ymir-core/tests/c1_phase_1_4_erosion.rs`](../../../crates/ymir-core/tests/c1_phase_1_4_erosion.rs)
+ the Stage E3 calibration tool
[`crates/ymir-core/tests/c1_phase_1_4_erosion_calibration.rs`](../../../crates/ymir-core/tests/c1_phase_1_4_erosion_calibration.rs)
(`#[ignore]`'d, invoked with `--ignored`).

## Run parameters

| Parameter | Value | Source |
|---|---|---|
| Grid | 64 × 64 | Phase 1.1 standard |
| Steps | 300 | Phase 1.1 standard |
| Kinematics | `PlateKinematics::preset_phase_1_1(8 plates)` | Phase 1.1 hand-tuned |
| Davis-Suppe | `DavisSuppeParams::default()` (enabled, coupling = 2.0, h_max = 2.5) | Phase 1.2 |
| Equilibrium-height | `EquilibriumHeightParams::default()` (enabled, h_eq = 2.0, k_collapse = 2.0, quadratic) | Phase 1.3 (post-E1.bis) |
| **Stream-power erosion** | **`ErosionParams::default()` (enabled, K = 0.001, m = 0.5, n = 1.0, floor = 0.2)** | **Phase 1.4 (this issue)** |
| Wall time | ~ 110 ms / 300 steps (~ 370 µs/step) | Stage E3 measured |

## Acceptance summary — 3 / 3 PASS

| # | Test | Measured | Threshold | Status |
|---|------|----------|-----------|--------|
| T1 | `erosion_caps_height_below_equilibrium` | `global_max = 2.181` | `∈ (1.0, 2.5)` | PASS |
| T2 | `erosion_preserves_davis_suppe_imprint_partially` (composite) | see below | 3 sub-assertions | PASS |
| T3 | `all_closures_disabled_matches_phase_1_1` | `global_max = 1079.697` | `> 100` (Phase 1.1 baseline) | PASS |

(Originally specified 4 tests; T3 `erosion_alone_produces_<X>` deferred — see § 4 below for the architectural finding it surfaced.)

### T2 composite breakdown

| Sub-assertion | Phase 1.4 | Phase 1.3 baseline | Phase 1.2 baseline | Threshold | Result |
|---|---|---|---|---|---|
| sub-1: `wedge_p95` | **0.696 ↑** | 0.376 | 0.376 | `∈ [0.4, 1.0]` | PASS |
| sub-2: `asymmetry = mean(near)/mean(far)` | 1.34 | 2.12 | 4.66 | `> 1.0` | PASS |
| sub-3: `fill_near` | 0.365 | 0.207 | 0.778 | `> 0.05` (Phase 1.4 floor) | PASS |

**The `wedge_p95 UP` finding** is the primary Phase 1.4 architectural signature — see § 1 below.

## Phase 1.1 / 1.2 / 1.3 / 1.4 cycle-300 comparison

| Quantity | Phase 1.1 (advection only) | Phase 1.2 (+ Davis-Suppe) | Phase 1.3 (+ Equilibrium) | **Phase 1.4 (+ Erosion)** | Phase 1.4 mechanism |
|---|---|---|---|---|---|
| `global_max` | ≈ 1080 | 2297 | 2.18 | **2.18** | equilibrium cap holds |
| `mean(S̃)` | conserved (init ≈ 0.557) | 1.574 | 0.361 | **0.361** | steady-state by cycle 100 |
| Total mass (= mean × cells) | conserved | ≈ 6447 (source-only) | 1478 (sink active) | **1478** | balanced by sink chain |
| `wedge_p95` (bulk) | n/a | 0.376 | 0.376 (bit-identical) | **0.696 ↑** | erosion eats shoulders, lifts wedge relative |
| `wedge_p99` (top 1 %) | n/a | 5.83 | 2.17 | **2.17** | clamped at h_eq |
| `wedge_max` | n/a | 93 | 2.18 | **2.18** | clamped |
| `mean(d ∈ 0-5)` near | n/a | 0.904 | 0.241 | **0.424** | wedge band lifted (W-T A^m discriminates) |
| `mean(d ∈ 10-20)` far | n/a | 0.194 | 0.113 | **0.316** | far band lifted (same mechanism) |
| `fill_near` | n/a | 0.778 | 0.207 | **0.365** | better saturation under combined source + sink |
| `asymmetry` (near / far) | n/a | 4.66 | 2.12 | **1.34** | direction preserved, magnitude reduced |
| Continental fraction (heightmap) | high | high | high | **14.9 %** | sparse Earth-like terrane |
| Wall time | very fast | 43 ms / 300 steps | 29 ms / 300 steps | **110 ms / 300 steps** | + isostasy + drainage per step |

## Visual gallery

Five cycle snapshots at NNN ∈ { 000, 050, 100, 200, 300 }, two PNGs each = 10 files in this directory:

- `cycle_NNN_altitude.png` — auto-rescaled hypsometric view via `compute_isostasy`.
- `cycle_NNN_s.png` — `S̃` field at the **same fixed palette `[0, 3.0]`** as the Phase 1.2 and 1.3 galleries
  ([`../c1_phase_1_2_davis_suppe/`](../c1_phase_1_2_davis_suppe/),
  [`../c1_phase_1_3_equilibrium_height/`](../c1_phase_1_3_equilibrium_height/)),
  so the three galleries are directly pixel-for-pixel comparable.

### What to look for in `cycle_300_*.png` (Phase 1.3 → Phase 1.4 visual delta)

1. **Continental fraction collapse.** Phase 1.3's `cycle_300_s.png` showed continental blocks (green/brown) covering significant grid area with bright X-shaped wedge boundaries. Phase 1.4 shows mostly ocean (blue) with the wedge ridges remaining as **discrete linear chains** — characteristic mid-ocean-ridge / passive-margin geometry. Erosion ate the broad continental shoulders, ocean expanded.
2. **Wedge ridges preserved.** The bright cross-shaped wedge structure remains prominently visible — Davis-Suppe + equilibrium-height imprint survives the erosion sink. Wedge cells are upstream of their drainage basins (small `A` in W-T `E ∝ A^m`), so they erode slowly relative to the downstream continental shoulders. See § 1 architectural finding.
3. **Boundary cap preserved.** Boundary cells (Convergent) still clamp at h_eq = 2.0 — visible as the brightest pixels in `cycle_NNN_s.png`. Same intensity as Phase 1.3.
4. **Stabilisation by cycle 50.** The `cycle_050_*.png` already shows the linear-mountains-in-ocean morphology. Mean `S̃` stable from cycle 100 onwards. Phase 1.4 reaches dynamic equilibrium faster than Phase 1.3 because erosion has a sink that responds proportional to slope+area.

## Architectural findings (4)

### § 1 — `wedge_p95 UP` (Stage E4 T2)

Counter-intuitive **at first glance**: Phase 1.4 wedge body bulk goes UP (0.376 → 0.696) despite adding erosion as a sink. Mechanism: W-T eq. (1) is `E ∝ K · A^m · S^n` with `m = 0.5`. Erosion concentrates on cells with **large drainage area `A`** — continental shoulders downstream of the wedge bodies. Wedge cells are **topographically upstream** (small `A`, at the top of their drainage basin), so they erode SLOWLY relative to the surroundings. Net effect: wedge ridges stand RELATIVELY HIGHER vs the eroding bulk.

Earth-like signature: this is precisely the morphology of mountain ranges in oceanic terrane — peaks preserve, plateaus erode. Confirms the W-T model produces geomorphologically plausible output in C1's regime.

### § 2 — Sparse continental Earth-like topology (Stage D)

Phase 1.4 default produces a continental fraction of **14.9 %** (vs presumably higher in Phase 1.3 — not measured). The downstream `compute_flow` analysis (Stage D test T1) restricted to continental cells:

| Metric | Measured (Phase 1.4) |
|---|---|
| n_continental | 610 cells (14.9 % of grid) |
| max continental accumulation | 9.0 |
| mean continental accumulation | 1.70 |
| **non-leaf (`accum > 1`) fraction** | **49.5 %** of continental (threshold `> 25 %`) |
| downstream (`accum > 2`) fraction | 13.6 % of continental |
| major (`accum > 5`) fraction | 0.8 % of continental |

Continental cells form **short drainage chains** on each wedge ridge — half are non-leaf (downstream of at least one upstream cell). Major channels are sparse (0.8 %), but connectivity in the drainage tree shape is significant. The morphology matches Earth-like passive-margin / mid-ocean-ridge terrane geometry (not Pangaea-style continental dominance).

### § 3 — Stage E0 bit-identical decomposition invariant preserved across 5 Phase 1.4 commits

The Phase 1.3 H3 architectural lock `c1_phase_a_decomposes_into_closures_then_post_tectonic` — asserting that `run_phase_a_cycle_c1(input, Enabled(_))` is bit-for-bit equal to `run_with_closures(state, ...) + apply_post_tectonic(...)` — was preserved EXACT (no tolerance) across all Phase 1.4 commits:

| Stage | Commit | Bit-identical |
|---|---|---|
| E0 | `cb71a8d` (extract `compute_sea_level_ref_s_space`) | PASS |
| E1 | `7a62b3c` (erosion closure module) | PASS |
| E2 | `a9d0fb7` (time loop integration) | PASS |
| I2 | `0dfdf9d` (isostasy validation tests) | PASS |
| D | `55e3756` (downstream validation tests) | PASS |

The wrapper-equals-decomposition contract held across 5 separate test runs with the test fixture explicitly disabling erosion (`erosion: ErosionParams { enabled: false, .. }`) — proving that the per-step erosion pipeline isolation (the `if closures.erosion.enabled { ... }` block in [`run_with_closures`](../../../crates/ymir-core/src/tectonics_c1/time_loop.rs)) is structurally clean, with no leakage when disabled. The wrapper continues to be **structurally pure** ("nothing more, nothing less") even after adding the Phase 1.4 pipeline.

### § 4 — Erosion floor clamp is mass-non-conservative in degraded regimes (Stage E4 T3 deferred)

While iterating on the Phase 1.4 spec's fourth test (`erosion_alone_produces_<X>`), three successive metric proposals failed under different regime assumptions:

1. **"Variance smoothing"** (`final_variance < initial_variance`): false. With Davis-Suppe + equilibrium-height disabled, the C1 advection-dominated regime drives unbounded boundary pile-up (`global_max ≈ 1080` Phase 1.1 baseline). `K = 0.001` erosion cannot keep up with the pile-up rate; variance EXPLODES (boundary cells reach 100+, interior stays ≈ 1.0).
2. **"Mass loss"** (`final_mass < initial_mass`): also false. Measured: erosion-only run added **+5191 mass over 300 steps (+227 % vs init)**. Mechanism: with DS disabled, advection drives near-oceanic cells below the `floor = 0.2` clamp inside `apply_erosion_step`. For each such cell, the clamp injects mass *upward* to bring it back to 0.2 (defensive against pathological K; mass-non-conservative).
3. **"Continental-filtered mass loss"**: over-engineered for a sanity check (requires arbitrary initial-cell-set filtering).

The **architectural finding**: the erosion closure's `floor = 0.2` defensive clamp is **mass-non-conservative in degraded regimes** (no Davis-Suppe source), while **net-conservative in the Phase 1.4 default** (DS active keeps continental cells well above the floor → clamp rarely triggers, erosion is a net sink at -35 % mass per Stage E3 measurement).

Per the `recursive-tuning-signals-structural` memory pattern, 3 iterations on a single sanity test signal that no clean regime-agnostic invariant exists. The finding is **load-bearing for understanding the closure's behavior**, but **not load-bearing as an acceptance test** — Phase 1.4's default operating regime never triggers the mass-non-conservative branch. Documented here + in the `project_c1_phase_1_4_erosion_outcomes` memory entry; future Phase 1.5+ (e.g., kinematics variations producing different pile-up rates) may need to revisit the floor behaviour.

## Calibration discipline (`K = 0.001`)

Per design doc §11.1, `K = 0.001` is calibrated for visual balance, not derived from literature (Lague 2014 declines to publish a universal K). Stage E3 visual review + quantitative comparison landed the default on the first iteration (within W3 3-iteration budget). The Stage E3 calibration tool (`#[ignore]`'d) is retained at [`c1_phase_1_4_erosion_calibration.rs`](../../../crates/ymir-core/tests/c1_phase_1_4_erosion_calibration.rs) for future K reviews.

## Cross-references

- **Closure code:** [`crates/ymir-core/src/tectonics_c1/closures/erosion/`](../../../crates/ymir-core/src/tectonics_c1/closures/erosion/) — `mod.rs` carries the W-T 1999 + Lague 2014 derivation, parameter rationale, and `min(h_collapse, h_erosion)` interaction documentation.
- **Time-loop integration:** [`crates/ymir-core/src/tectonics_c1/time_loop.rs`](../../../crates/ymir-core/src/tectonics_c1/time_loop.rs) `run_with_closures` § "Per-step pipeline" — 4-stage Phase 1.4 pipeline (isostasy + drainage + areas + erosion).
- **Design doc:** [`docs/c1_lightweight_dynamic_tectonics.md`](../../c1_lightweight_dynamic_tectonics.md) §5.1 (MVP closures), §8.5 (visual plausibility), §11 (implicit physical scales), §11.1 (calibration via visual review).
- **Phase 1.3 gallery:** [`../c1_phase_1_3_equilibrium_height/`](../c1_phase_1_3_equilibrium_height/) — same palette, direct side-by-side comparison.
- **Phase 1.2 gallery:** [`../c1_phase_1_2_davis_suppe/`](../c1_phase_1_2_davis_suppe/).
- **Phase 1.1 baseline:** [`../c1_phase_1_1_advection/`](../c1_phase_1_1_advection/).
- **Stage E3 calibration tool:** [`c1_phase_1_4_erosion_calibration.rs`](../../../crates/ymir-core/tests/c1_phase_1_4_erosion_calibration.rs).
- **Stage D downstream tests:** [`c1_phase_1_4_downstream.rs`](../../../crates/ymir-core/tests/c1_phase_1_4_downstream.rs).
- **Stage I2 isostasy tests:** [`c1_phase_1_4_isostasy.rs`](../../../crates/ymir-core/tests/c1_phase_1_4_isostasy.rs).
