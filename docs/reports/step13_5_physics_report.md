# Step 13.5 — Physics report (oceanic FBM extension)

> Companion to `step13_5_regression_report.md`. Validates that the
> opt-in oceanic FBM extension produces measurable bathymetric
> heterogeneity at the calibrated default amplitude, that the
> strict `OCEANIC_CLAMP_MAX = 0.49` threshold prevents accidental
> threshold-crossing, and that solver health is preserved.
> Acceptance criteria #1–#9, #11, #15, #16.

## Methodological note

Step 13.5 is a focused extension of Step 13 (continental FBM
heterogeneity) that adds opt-in FBM perturbation on oceanic
cells. Init-only mechanism, no dynamics changes. Default
behaviour preserves Step 13 bit-identical via a structural
short-circuit on the `apply_fbm_to_oceanic` flag — no second
`Fbm<Perlin>` instance, no oceanic cell touched when the flag
is `false`.

The Phase 4 calibration sweep produced one notable finding
that justifies a focused §4.13 amendment (the milestone's
**7th** design-driven amendment):

- **`noise::Fbm<Perlin>::σ ≈ 0.27 × amplitude` reproduced on
  the oceanic side**: the same auto-normalisation pattern
  Step 13 Phase 6 surfaced for the continental FBM holds for
  the oceanic FBM. Sweep confirmed σ scales linearly with
  amplitude across the calibration grid (`amp ∈ {0.05, 0.10,
  0.15, 0.20, 0.25}`), and is essentially insensitive to
  `fbm_scale_oceanic` over `[0.05, 0.20]` (variation < 5 %)
  — the wavelength-vs-`L_plate` ratio is comfortably under 1
  on every column. The `None` default for
  `fbm_scale_oceanic` (= reuse continental `fbm_scale = 0.10`)
  is therefore empirically justified, not just by parsimony.

## Test setup

All measurements at 64² unless noted; 32² spot-checks confirm
the mechanism scales without surprise. Voronoï layouts:
`single_continent` (seed=12, 4 plates, 50 % continental) for
the per-acceptance probes — same as Step 13's calibration so
the two sides compose cleanly. Multi-preset galerie additionally
covers `convergence` (seed=23, 6 plates, 40 % continental).

| Source | Tests |
|---|---|
| `crates/ymir-core/tests/v2_step13_5_acceptance.rs` | acceptance #7 (32² + 64²), `fbm_oceanic_calibration_probe` |
| `crates/ymir-core/tests/v2_step13_5_cg_ratio.rs` | acceptance #9 |
| `crates/ymir-core/tests/v2_step13_visual_checkpoint.rs` | acceptance #16 (single-preset amplitude sweep + multi-preset disabled/enabled gallery) |
| `crates/ymir-core/src/tectonics_v2/init/radial_profile_fbm.rs` (test mod) | acceptance #1, #2, #3, #4, #5, #6 |

## Phase 4 calibration sweep

`single_continent` 64², per (`fbm_amplitude_oceanic`,
`fbm_scale_oceanic`): `σ_fbm_oceanic_isolated / max(S̃_oceanic) /
clip%`.

```text
amp \ scale  | 0.05         0.10         0.15         0.20
amp=0.05     | σ=0.014/0.24/0%  σ=0.013/0.24/0%  σ=0.013/0.23/0%  σ=0.014/0.24/0%
amp=0.10     | σ=0.028/0.28/0%  σ=0.027/0.27/0%  σ=0.026/0.27/0%  σ=0.028/0.28/0%
amp=0.15     | σ=0.042/0.32/0%  σ=0.040/0.31/0%  σ=0.039/0.30/0%  σ=0.042/0.31/0%
amp=0.20     | σ=0.055/0.36/0%  σ=0.053/0.34/0%  σ=0.052/0.34/0%  σ=0.056/0.35/0%
amp=0.25     | σ=0.069/0.40/0%  σ=0.067/0.38/0%  σ=0.065/0.37/0%  σ=0.070/0.39/0%
```

Findings:

1. **Clip fraction = 0 % across the entire grid** — even the
   most aggressive `(amp=0.25, scale=0.05)` pair gives
   `max(S̃_oceanic) ≈ 0.40 < OCEANIC_CLAMP_MAX = 0.49`. The
   threshold protection has comfortable headroom. The strict
   `0.49` bound never fires at sane amplitudes.
2. **σ scales linearly with amplitude** (`σ ≈ 0.27 × amp`),
   reproducing the Step 13 continental finding for
   auto-normalised `noise::Fbm<Perlin>`.
3. **σ insensitive to scale** in `[0.05, 0.20]` (< 5 %
   variation between columns).
4. **Target `σ_fbm_oceanic_isolated ∈ [0.02, 0.08]`** is hit
   by `amp ∈ [0.10, 0.25]`. `amp = 0.15` (Phase 5 default
   choice) lands at `σ ≈ 0.040`, **mid-band with margin on
   both sides**.

## Acceptance #1–#6 — algorithmic correctness

Established by the Phase 1 unit tests in
`init::radial_profile_fbm::tests` (six new tests covering
acceptance #1–#6):

| Test | Acceptance | Result |
|---|---|---|
| `oceanic_fbm_disabled_preserves_step13` | #1 — disabled flag short-circuits, output byte-identical to Step 13 | ✓ |
| `oceanic_fbm_enabled_varies` | #2 — enabled flag → oceanic variance > 0 | ✓ |
| `oceanic_fbm_no_threshold_crossing` | #3 — sweep `amp ∈ {0.05, 0.10, 0.20, 0.30, 0.40}`, every oceanic cell stays in `[0, 0.49]` | ✓ |
| `oceanic_fbm_seed_independence` | #4 — distinct `fbm_seed_oceanic` → distinct oceanic; continental insulated | ✓ |
| `oceanic_fbm_seed_default_derivation` | #5 — `None` derives from `fbm_seed XOR 0xC0FFEE` byte-for-byte | ✓ |
| `oceanic_fbm_scale_independence` | #6 — distinct `fbm_scale_oceanic` → distinct oceanic spectral content | ✓ |

The disabled-flag test is the structural-short-circuit oracle:
it builds with bogus `fbm_amplitude_oceanic = 0.42`,
`fbm_scale_oceanic = Some(0.07)`, `fbm_seed_oceanic = Some(0xDEAD)`
but flag = `false`, and asserts byte-identical output to a
build with the disabled defaults. The flag must short-circuit
before any oceanic param is read.

## Acceptance #7 — oceanic FBM contribution measurable

`σ_fbm_oceanic_isolated ∈ [0.02, 0.08]` over oceanic cells with
the Phase 5 calibrated defaults
(`FBM_AMPLITUDE_OCEANIC_DEFAULT = 0.15`, `fbm_scale_oceanic =
None` ⇒ reuse `0.10`, `fbm_seed_oceanic = None` ⇒ XOR derive):

| Grid | oceanic cells | σ_fbm_oceanic_isolated | max(S̃_oceanic) | clip% |
|------|---------------|------------------------|-----------------|-------|
| 64²  | 928           | **0.0400** ✓           | 0.307           | 0 %   |
| 32²  | 232           | **0.0398** ✓           | 0.295           | 0 %   |

PASS at both grids with margin (mid-band of `[0.02, 0.08]`,
max well under the 0.49 clamp). σ is essentially
grid-independent (Δσ ≈ 0.0002 between 32² and 64²) — the 232
oceanic cells at 32² are well above the 150-cell sample-noise
threshold introduced in Step 13's caveats, so per-aggregate
`σ` is statistically reliable here.

## Acceptance #8 — bathymetric range bounded

By construction: the algorithm's `clamp(perturbed, 0,
OCEANIC_CLAMP_MAX)` guarantees `min(S̃_oceanic) ≥ 0` and
`max(S̃_oceanic) ≤ 0.49`. Validated by acceptance #3's
amplitude sweep test (every cell, every amplitude) and
empirically by the calibration sweep (clip % = 0 at the
default amplitude).

## Acceptance #9 — CG iters ratio

`tests/v2_step13_5_cg_ratio.rs::oceanic_fbm_cg_ratio_acceptance`
runs the same 64² × 30-step mantle-on shape twice, varying
only `apply_fbm_to_oceanic`:

| Mode                              | cg_iter_mean |
|-----------------------------------|--------------|
| oceanic_disabled (Step 13 path)   | 1384.42      |
| oceanic_enabled (Step 13.5 path)  | 1384.87      |

`ratio = 1.000×` — solver health is **transparent** to the
oceanic FBM. PASSES the [0.90, 1.10] band with maximum
margin.

The disabled-flag baseline (1384.42) reproduces Step 13's
Phase 7 RadialProfileWithFBM CG mean (1384.4) — independent
confirmation that Phase 1's structural short-circuit is
bit-identical.

### Aside — why oceanic FBM is solver-transparent

Step 13's continental FBM extension reduced CG iters by
≈ 5 % vs the `Uniform` baseline (smoother init = better
Newton warm-start). The oceanic FBM extension is neutral
(1.000×). Plausible explanation: oceanic cells sit at low
`S̃ ≈ 0.20`, where the Stokes operator's stiffness coefficient
(`η ∝ S̃²`) is small — the preconditioner traverses that
low-stiffness band regardless of bathymetric perturbation.
Conditioning is dominated by continental cells (`S̃ ≈ 0.95`),
which are unchanged when the oceanic flag flips.

Operationally: enabling oceanic FBM on an existing run carries
no measurable solver-health cost. Cheap to ship behind the
opt-in flag.

## Acceptance #11 — mass conservation

Governed by Step 8's contract; oceanic FBM is an **init-only**
modification of the initial S̃ field, the time-loop dynamics
are unchanged, and the mass conservation invariant rides on
the existing `step_upwind` / `boundary_q` machinery. No new
breakage path; covered by the Step 11 / Step 12 regression
tests (which both pass — see regression report).

## Acceptance #15 + #16 — UI controls and visual output

UI: the parameter panel's `RadialProfileWithFBM` block now
includes a separate "Oceanic FBM noise (Step 13.5)" section
behind the `apply_fbm_to_oceanic` toggle. Conditional sliders
for amplitude, scale (with "Use continental scale" checkbox
mapping to `None`), seed (with "Derive from continental seed
XOR 0xC0FFEE" checkbox mapping to `None`), and a randomize
button reusing the SplitMix64 mixer from the continental
side. **Acceptance #16 tooltip warning** — when
`fbm_amplitude_oceanic > OCEANIC_CLAMP_MAX − oceanic_value`,
an italic informational message tells the user that oceanic
cells may saturate at the 0.49 clamp and that volcanic
islands are a separate Step 13.6 if pursued. Not a hard
block.

Visual: three patchworks under
`docs/reports/step13_5_visual_checkpoint/`.

### Phase 4 amplitude sweep (single-preset)

`single_continent` 64² with four oceanic-FBM configurations.

![patchwork standard](step13_5_visual_checkpoint/patchwork_oceanic_amp_sweep_64sq.png)

Layout: `Step 13 default (uniform)` | `amp=0.10` | `amp=0.20` |
`amp=0.40 (clipping demo)`.

Per-tile oceanic stats (mean / std / range / clip%):

| Tile | oceanic mean | oceanic std | range | clip% |
|------|--------------|-------------|-------|-------|
| Step 13 default            | 0.200 | 0.000 | [0.20, 0.20] | 0 % |
| amp = 0.10                 | 0.200 | 0.027 | [0.13, 0.27] | 0 % |
| amp = 0.20                 | 0.199 | 0.053 | [0.07, 0.34] | 0 % |
| amp = 0.40 (clipping demo) | 0.198 | 0.106 | [0.00, 0.49] | 0 % |

The `amp=0.40` tile reaches `max ≈ 0.485` — right at the edge
of the upper clamp without crossing, with the **lower** clamp
at `0.0` activating on FBM negative tails (visible as
dark patches in the standard view). At Phase 5's calibrated
`amp = 0.15` the perturbation stays comfortably inside the
linear-σ regime.

### Phase 4 amplitude sweep — oceanic-zoomed view

Continental cells blanked to mid-grey, oceanic cells remapped
from `[0, 0.49] → [0, 1]` so the FBM signature is directly
readable.

![patchwork oceanic zoom](step13_5_visual_checkpoint/patchwork_oceanic_zoom_64sq.png)

The progression Step 13 default → amp=0.10 → amp=0.20 →
amp=0.40 reads as: uniform mid-grey, light speckle, clear
bathymetric texture, pronounced texture with floor saturation.

### Phase 7 multi-preset gallery (disabled vs enabled)

Two presets × two modes (oceanic disabled / oceanic enabled
with the Phase 5 default `amp = 0.15`).

![galerie disabled vs enabled](step13_5_visual_checkpoint/galerie_oceanic_disabled_vs_enabled_64sq.png)

Oceanic-zoomed:

![galerie zoom](step13_5_visual_checkpoint/galerie_oceanic_zoom_disabled_vs_enabled_64sq.png)

Layout: rows = `single_continent` | `convergence`; cols =
`oceanic_disabled (Step 13)` | `oceanic_enabled (Step 13.5,
amp=0.15)`.

Cross-preset oceanic stats with the Phase 5 default amplitude:

| Preset             | oceanic mean | oceanic std | max(S̃_oceanic) |
|--------------------|--------------|-------------|------------------|
| single_continent   | 0.199        | 0.040       | 0.307            |
| convergence        | 0.201        | 0.040       | 0.311            |

σ is essentially identical across presets — the mechanism is
preset-shape-insensitive. Same finding pattern as Step 13's
continental FBM at this resolution.

## Caveats

### 1. Volcanic islands deferred (Step 13.6 if pursued)

Threshold protection at `OCEANIC_CLAMP_MAX = 0.49` is
**strict** — oceanic cells cannot cross the `0.5` continental
classification threshold via FBM perturbation regardless of
amplitude. Volcanic islands (cells crossing the threshold via
deliberate emergence) involve game-design decisions that
benefit from a dedicated step (Step 13.6 if pursued) with its
own calibration. Documented in §4.13 (D7) and surfaced via
the UI tooltip when the user pushes amplitude beyond the
expected linear-σ regime.

### 2. `fbm_scale_oceanic` grid scaling

`fbm_scale_oceanic` is in domain fractions (same convention
as `fbm_scale`), so wavelength scales linearly with grid
resolution. The default `None` (= reuse continental scale
0.10) gives `wavelength ≈ 6 cells` at 64² and `≈ 3 cells` at
32², a comfortable wavelength-vs-`L_plate` ratio at both
grids. For grids ≥ 256² the user should raise the scale to
maintain the same physical wavelength relative to plate
sizes — same recommendation as Step 13.

### 3. Reformulation of acceptance #7 inherited from Step 13

The original Step 13 issue had `σ_total ∈ [0.04, 0.10]` for
oceanic acceptance #7 — a vacuous-truth shape that Step 13
Phase 6 corrected to `σ_fbm_isolated ≥ 0.040` (lower bound
only, FBM contribution measured directly via subtraction of
the FBM-disabled baseline). Step 13.5 inherits this
reformulation: the test asserts `σ_fbm_oceanic_isolated ∈
[0.02, 0.08]` (issue D5 target) measured the same way. With
the Step 13 oceanic baseline being a uniform constant, the
FBM-isolated and total-variance forms agree numerically here,
but the FBM-isolated form is robust to any future change of
the oceanic baseline (e.g., per-plate variation).

### 4. 7th design-driven amendment of the milestone

Step 13.5 produced one finding worthy of an amendment patch:
the `noise::Fbm<Perlin>::σ ≈ 0.27 × amplitude` auto-
normalisation. This reproduces Step 13's continental-side
finding on the oceanic side; same coefficient, same
implication for default-amplitude calibration. Folded into
the §4.13 amendment alongside the Step 13.5 design notes.
The pattern of "issue D-decision lands → empirical probing
during implementation either confirms or revises the premise
→ amendment patch" has now held for 14 steps.

### 5. Solver-transparency finding

The CG ratio = 1.000× neutrality (vs Step 13's 0.95×
continental gain) is a *finding*, not a target. Documented
above. Does not influence the default-amplitude choice — that
was anchored on the σ target band — but worth tracing for
the next consumer of these results.
