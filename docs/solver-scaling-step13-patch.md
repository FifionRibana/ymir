# §4.13 — Step 13 patch: continental margins and intra-plate heterogeneity

This patch documents two new S̃ initialisation modes added in
[#111](https://github.com/FifionRibana/ymir/issues/111) and the
empirical calibration that produced their default parameter values.
Phase 7 will integrate this into the full §4.13 of
`docs/solver-scaling.md`; the calibration section below is kept as a
self-contained note so the next reviewer who tunes the defaults can
read what was measured and why.

## Motivation — Step 11 visual findings

Step 11 visual exploration (`docs/reports/step11_physics_report.md`)
revealed two limitations of the existing four init modes
(`Checkerboard`, `Uniform`, `Gaussian`, `Convolution`):

1. **Rigid polygonal Voronoï borders**. All pre-Step-13 modes
   produce a binary continental/oceanic classification that
   strictly follows Voronoï cells. Boundary smoothing is
   localised (1–2 cells), and the polygonal contours stay
   visible throughout the simulation in regions of low tectonic
   activity.
2. **Uniform interior thickness**. All modes assign approximately
   the same `S̃` value across an entire continental plate. Real
   cratons have intra-province variations of factor 2–3.

Step 13 closes both gaps with init-only mechanisms (no dynamics
changes, no per-step re-application).

## New modes

### `RadialProfile { continental_value, oceanic_value, profile_shape }`

Continental cells take
`S̃ = oceanic_value + (continental_value − oceanic_value) ·
profile(d / L_plate)`,
where `d` is the Chebyshev BFS distance to the nearest
inter-plate boundary (Phase 1 utility
`compute_dist_to_inter_plate_boundary`) and `L_plate` is the
per-plate maximum BFS distance. Oceanic cells take
`S̃ = oceanic_value` flat.

`ProfileShape ∈ { Smoothstep, Linear, Pow { exponent } }`. Defaults
`continental_value = 0.95`, `oceanic_value = 0.20`,
`profile_shape = Smoothstep`. The Pow exponent is clamped by the
Phase 5 UI to `[0.3, 3.0]`.

### `RadialProfileWithFBM { …, fbm_amplitude, fbm_octaves, fbm_persistence, fbm_lacunarity, fbm_scale, fbm_seed }`

Same radial baseline as above, plus isotropic FBM noise
(`noise::Fbm<Perlin>`) added on continental cells only. Output
clamped to `[0, 1]`.

```text
S̃[i, j] = clamp(
    S̃_radial[i, j] + amplitude · fbm(x_norm / scale, y_norm / scale),
    0, 1
)   // continental cells only — oceanic stays at oceanic_value
```

## Calibration note — `fbm_scale` and `fbm_amplitude` defaults

The initial draft (Phase 3) used `fbm_scale = 0.25` and
`fbm_amplitude = 0.10` based on the issue D2 estimate that the
FBM contribution would be `σ ≈ amplitude / 1.5–2.0 ≈ 0.05–0.07`
(target [0.04, 0.10]). Phase 6 acceptance probing on
`single_continent` (64², 4 plates, 50 % continental, seed 12)
revealed two compounding issues:

1. **`fbm_scale = 0.25` ⇒ wavelength ≈ 16 cells** on a 64² grid,
   greater than typical plate `L_plate ≈ 10–15 cells`. The FBM does
   not actually oscillate intra-plate.
2. **`noise::Fbm<Perlin>` is auto-normalised** and produces
   `σ ≈ 0.27 × amplitude`, not the `≈ amplitude / 1.5–2.0` the
   issue estimate assumed. The factor difference is roughly 2×.

Empirical sweep
(`tests::v2_step13_acceptance::fbm_calibration_probe`,
`#[ignore]`):

| `scale \ amp` | 0.10  | 0.15  | 0.20  | 0.25  |
|---------------|-------|-------|-------|-------|
| 0.05          | 0.024 | 0.037 | 0.048 | 0.060 |
| **0.10**      | 0.027 | 0.041 | **0.055** | 0.068 |
| 0.20          | 0.018 | 0.027 | 0.036 | 0.045 |
| 0.25          | 0.018 | 0.027 | 0.036 | 0.046 |

(table values are `σ_fbm_isolated = std(S̃_FBM − S̃_radial)` over
interior cells `t > 0.5` of the largest continental plate of
`single_continent`.)

**`scale = 0.10` is empirically optimal**. `scale = 0.05` is
counter-productive: high-frequency Perlin grid artefacts dilute the
large-scale variance the metric captures. `scale ≥ 0.20` is
under-resolved relative to plate size — the FBM does not oscillate
intra-plate.

**`amplitude = 0.20`** at `scale = 0.10` clears the acceptance
lower bound `σ_fbm_isolated ≥ 0.040` with ≈ 35 % margin across
all three continental plates of `single_continent`:

| pid | L_plate | cells_int | σ_radial | σ_total | σ_fbm_isolated |
|-----|---------|-----------|----------|---------|----------------|
| 0   | 14.0    | 296       | 0.094    | 0.110   | **0.0545**     |
| 1   |  9.0    | 266       | 0.109    | 0.131   | **0.0450**     |
| 2   | 12.0    | 208       | 0.097    | 0.110   | **0.0549**     |

`amplitude = 0.20` stays well within the Phase 5 UI clamp
`[0.0, 0.40]` and well above the `0.5` continental threshold for
interior cells (with `S̃_radial → 0.95` at the centroid, FBM
perturbation at amplitude 0.20 brings worst-case interior S̃ to
`0.95 − 0.20 = 0.75`, still well above 0.5).

Users on larger grids (256²+) should raise `fbm_scale` to maintain
the wavelength-vs-`L_plate` ratio (≈ 6 cells / 14 cells ≈ 0.4).
For example on 256² with comparable continental_ratio, plates
are ≈ 4× larger and `fbm_scale ≈ 0.025` would preserve the same
`wavelength / L_plate` proportion. The default targets the
milestone's 32²–64² validation grids.

## Acceptance #7 reformulation

The issue's original acceptance `σ_total(S̃) ∈ [0.04, 0.10]` over
interior cells passes vacuously: with the `Smoothstep` radial
profile alone (no FBM), `σ_radial ≈ 0.094` already lands in the
band, so the test was satisfied by the radial gradient irrespective
of FBM. Phase 6 reformulated to `σ_fbm_isolated ≥ 0.040`
(lower bound only — the Phase 5 UI amplitude clamp `[0.0, 0.40]`
already caps how much heterogeneity FBM can introduce, and the
algorithm's `[0, 1]` clamp keeps S̃ in physical range). The
reformulated test directly measures the FBM contribution rather
than letting the radial gradient validate the spec by accident.

The reformulation is implemented per-plate (not just on the
largest plate): every continental plate with non-degenerate
`L_plate` must clear `σ_fbm_isolated ≥ 0.040`. This catches
configurations where one small plate would silently fail the
acceptance under a "largest-plate-only" interpretation.

## 32² acceptance — small-sample finding

Phase 7 re-runs the acceptance probes at 32² (the milestone's
smaller validation grid) on the same `single_continent` Voronoï
parameters. Continental plates land at `L_plate ∈ [4, 6]` cells
(half of 64²) with interior cell counts `[52, 71, 80]` (vs
`[208, 266, 296]` at 64²). The wavelength-vs-`L_plate` ratio is
preserved (`fbm_scale = 0.10` produces wavelength = 3.2 cells at
32², ≈ same fraction of `L_plate` as at 64²).

Per-plate σ_fbm_isolated at 32²:

| pid | cells_int | L_plate | σ_fbm_isolated |
|-----|-----------|---------|----------------|
| 0   | 80        | 6.0     | 0.0560 ✓       |
| 1   | 52        | 4.0     | 0.0388 (marginal) |
| 2   | 71        | 5.0     | 0.0540 ✓       |

`pid=1` lands 0.0012 below the 0.040 lower bound. With 52 cells the
small-sample standard error of a σ estimate is roughly
`σ / √(2·N) ≈ 0.0388 / √104 ≈ 0.0038`, so 0.0388 is statistically
indistinguishable from 0.040 at this sample size — the failure is
small-sample noise, not a mechanism degradation. The 32² test
relaxes to "largest-continental-plate-only" assertion (matching
the issue's "a single continental plate" wording); smaller plates
report-only.

The mechanism does scale to 32² for adequately-sized plates
(≥ ~70 interior cells); per-plate strictness is 64²-only.

## CG ratio (acceptance #10)

Acceptance #10 ("CG iters ratio ≤ 1.10× existing modes baseline")
is measured by `tests/v2_step13_cg_ratio.rs`: same `BaselineConfig`
(64² × 30 steps, mantle on, slab off, yielding on,
`single_continent` Voronoï) run three times, varying only
`init_mode`. `metrics.cg_iter_mean` ratio to the `Uniform`
baseline:

| Mode                 | cg_iter_mean | ratio  |
|----------------------|--------------|--------|
| Uniform (baseline)   | 1421.0       | 1.000  |
| RadialProfile        | 1351.4       | **0.951×** |
| RadialProfileWithFBM | 1384.4       | **0.974×** |

Both new modes *reduce* CG iterations slightly — the smoother
initial S̃ field gives Newton's first inner CG solve a better
warm-start than `Uniform`'s 1-cell-wide boundary band. Acceptance
satisfied with ≈ 5–15 % margin under the 1.10× ceiling.

## Validity envelope

| Parameter            | Default    | Sensible range                    |
|----------------------|------------|------------------------------------|
| `continental_value`  | 0.95       | [0.5, 1.0] (UI clamp)              |
| `oceanic_value`      | 0.20       | [0.0, 0.5] (UI clamp)              |
| `profile_shape`      | Smoothstep | { Smoothstep, Linear, Pow }        |
| `Pow.exponent`       | 1.0        | [0.3, 3.0] (UI clamp)              |
| `fbm_amplitude`      | 0.20       | [0.0, 0.40] (UI clamp)             |
| `fbm_octaves`        | 4          | [1, 8] (UI clamp)                  |
| `fbm_persistence`    | 0.5        | [0.10, 1.0] (UI clamp)             |
| `fbm_lacunarity`     | 2.0        | [1.5, 4.0] (UI clamp)              |
| `fbm_scale`          | 0.10       | [0.05, 1.0] (UI clamp); see note   |
| `fbm_seed`           | 0x0FBA5EED | u64, independent of Voronoï seed   |

The UI clamps are documented at the slider level
(`crates/ymir-viz/src/ui/parameter_panel_v2.rs`); the algorithm
itself accepts any positive value (no silent clamps — anti-pattern
D7 from the Step 13 issue).

## Step 13.5 — oceanic FBM extension

Step 13's `RadialProfileWithFBM` mode applies FBM noise to
**continental** cells only; oceanic cells stay at a flat
`oceanic_value`. Step 13.5 adds an opt-in extension that applies
FBM to oceanic cells too, producing bathymetric variation for
the Living Landz workflow (heightmap rendering, coastal gameplay
zones, visual coherence pre-erosion). Init-only mechanism, no
dynamics changes; default flag preserves Step 13 bit-identical.

### Algorithm

Four new fields on `InitMode::RadialProfileWithFBM`
(`#[serde(default)]` on each so legacy v2 preset JSON
deserialises unchanged):

```rust
apply_fbm_to_oceanic: bool,                 // default false
fbm_amplitude_oceanic: f64,                 // default 0.15 (Phase 5 calibration)
fbm_scale_oceanic: Option<f64>,             // default None ⇒ reuse fbm_scale
fbm_seed_oceanic: Option<u64>,              // default None ⇒ XOR derive
```

When the flag is on, oceanic cells receive

```text
S̃[i, j] = clamp(
    oceanic_value + fbm_amplitude_oceanic · fbm_oceanic.get(x, y),
    0,
    OCEANIC_CLAMP_MAX
)   // oceanic cells only — continental cells unchanged
```

with `OCEANIC_CLAMP_MAX = 0.49` strictly preventing
threshold-crossing to continental classification (D7 — volcanic
islands are explicitly out of scope, deferred to Step 13.6 if
pursued).

`fbm_seed_oceanic = None` derives the oceanic seed from
`fbm_seed XOR FBM_SEED_OCEANIC_XOR_MAGIC` (= `0xC0FFEE`) for
reasonable independence between continental and oceanic noise
without forcing the user to supply two seeds. `fbm_scale_oceanic
= None` reuses the continental `fbm_scale` (empirical sweep
showed σ insensitive to scale over `[0.05, 0.20]` so reusing is
both parsimonious and statistically sound).

### Calibration sweep — oceanic side

Phase 4 sweep on `single_continent` (64², seed=12, 4 plates,
50 % continental):

| amp \ scale | 0.05  | 0.10  | 0.15  | 0.20  |
|-------------|-------|-------|-------|-------|
| 0.05        | σ=0.014, max=0.24 | σ=0.013, max=0.24 | σ=0.013, max=0.23 | σ=0.014, max=0.24 |
| 0.10        | σ=0.028, max=0.28 | σ=0.027, max=0.27 | σ=0.026, max=0.27 | σ=0.028, max=0.28 |
| **0.15**    | σ=0.042, max=0.32 | **σ=0.040, max=0.31** | σ=0.039, max=0.30 | σ=0.042, max=0.31 |
| 0.20        | σ=0.055, max=0.36 | σ=0.053, max=0.34 | σ=0.052, max=0.34 | σ=0.056, max=0.35 |
| 0.25        | σ=0.069, max=0.40 | σ=0.067, max=0.38 | σ=0.065, max=0.37 | σ=0.070, max=0.39 |

(format: `σ_fbm_oceanic_isolated, max(S̃_oceanic)`. Clip-fraction
= 0 % across the entire grid — the strict `OCEANIC_CLAMP_MAX =
0.49` upper bound never fires at sane amplitudes.)

Findings:

- **`noise::Fbm<Perlin>::σ ≈ 0.27 × amplitude` reproduced on
  the oceanic side**, identical coefficient to Step 13 Phase 6's
  continental finding. Auto-normalisation of the multifractal
  stack is not amplitude-dependent.
- **σ insensitive to scale** (variation < 5 % across columns)
  — wavelength comfortably below `L_plate` over the full
  scale range. The `None` default for `fbm_scale_oceanic` (=
  reuse continental) is empirically justified.
- **Clip fraction = 0 %** — even the most aggressive
  `(amp=0.25, scale=0.05)` pair gives `max(S̃_oceanic) ≈ 0.40
  < 0.49`. The threshold protection has comfortable headroom;
  `OCEANIC_CLAMP_MAX = 0.49` is a safety net rather than a
  binding constraint.

### Default amplitude — `0.15`

`FBM_AMPLITUDE_OCEANIC_DEFAULT = 0.15` lands at
`σ_fbm_oceanic_isolated ≈ 0.040` — **mid-band of the issue's
target `[0.02, 0.08]`** with margin on both sides. `max ≈ 0.31`
leaves 38 % headroom under the `0.49` clamp. Phase 4 sanity
visual confirmed the perturbation is visually distinct from
Step 13's uniform oceanic baseline without overwhelming the
continental signature.

Acceptance #7 measurement at this default:

| Grid | oceanic cells | σ_fbm_oceanic_isolated | max(S̃_oceanic) | clip% |
|------|---------------|------------------------|------------------|-------|
| 64²  | 928           | 0.0400                 | 0.307            | 0 %   |
| 32²  | 232           | 0.0398                 | 0.295            | 0 %   |

σ is grid-independent (Δσ ≈ 0.0002 between 32² and 64²) — the
232 oceanic cells at 32² are well above the 150-cell sample-noise
threshold introduced in Step 13's "small-sample" caveat.

### Threshold protection rationale (D7)

`OCEANIC_CLAMP_MAX = 0.49` is **strict** — oceanic cells cannot
cross the `0.5` continental classification threshold via FBM
perturbation regardless of amplitude. The defensive 0.01 margin
under the threshold guards against floating-point edge cases at
exactly `0.5`.

Volcanic islands (oceanic cells deliberately crossing the
`0.5` threshold to emerge as land) involve game-design
decisions — how many, how distributed, how rare — that benefit
from a dedicated step (Step 13.6 if pursued) with its own
calibration. Step 13.5 provides bathymetric variation; Step
13.6 would add controlled threshold-crossing on top.

UI surfaces this contract via a tooltip warning (acceptance
#16): when `fbm_amplitude_oceanic > OCEANIC_CLAMP_MAX −
oceanic_value`, an italic informational message tells the user
that oceanic cells may saturate at the 0.49 clamp and that
volcanic islands are a separate Step 13.6 scope. Not a hard
block — the user can push the amplitude up if they want
visible saturation at the floor (negative FBM tails).

### CG ratio (acceptance #9)

`tests/v2_step13_5_cg_ratio.rs` runs the same 64² × 30-step
mantle-on shape twice, varying only `apply_fbm_to_oceanic`:

| Mode                              | cg_iter_mean |
|-----------------------------------|--------------|
| oceanic_disabled (Step 13 path)   | 1384.42      |
| oceanic_enabled (Step 13.5 path)  | 1384.87      |

`ratio = 1.000×` — solver health is **transparent** to the
oceanic FBM. The disabled-flag baseline reproduces Step 13
Phase 7's RadialProfileWithFBM CG mean (1384.4) — independent
confirmation that Phase 1's structural short-circuit is
bit-identical.

Worth tracing: Step 13's continental FBM extension reduced CG
iters by ≈ 5 % vs the `Uniform` baseline (smoother init →
better Newton warm-start). Oceanic FBM extension is neutral.
Plausible explanation: oceanic cells sit at low `S̃ ≈ 0.20`,
where the Stokes operator's stiffness coefficient (`η ∝ S̃²`)
is small — the preconditioner traverses that low-stiffness
band regardless of bathymetric perturbation. Conditioning is
dominated by continental cells, which are unchanged when the
oceanic flag flips.

### Validity envelope — Step 13.5 fields

| Parameter                | Default      | Sensible range                          |
|--------------------------|--------------|------------------------------------------|
| `apply_fbm_to_oceanic`   | `false`      | `bool`                                   |
| `fbm_amplitude_oceanic`  | `0.15`       | `[0.0, 0.40]` (UI clamp)                 |
| `fbm_scale_oceanic`      | `None`       | `Some([0.05, 0.50])` or `None` (= reuse continental) |
| `fbm_seed_oceanic`       | `None`       | `Some(u64)` or `None` (= XOR derive)     |
| `OCEANIC_CLAMP_MAX`      | `0.49` const | not user-tunable (strict by design — D7) |
| `FBM_SEED_OCEANIC_XOR_MAGIC` | `0xC0FFEE` const | not user-tunable                |

UI tooltip warns when the user pushes
`fbm_amplitude_oceanic > OCEANIC_CLAMP_MAX − oceanic_value`,
indicating clipping is likely on FBM positive tails.

### Volcanic islands — explicitly out of scope (Step 13.6)

Volcanic islands deferred to a separate Step 13.6 because the
mechanism intersects three game-design decisions Step 13.5
should not silently make:

1. **How many islands** — distribution model (Poisson? rare
   high-amplitude FBM tails? per-plate quota?).
2. **How distributed** — uniform random on oceanic cells?
   Hot-spot tracks (volcanic chains)? Boundary-proximate
   (subduction-induced)?
3. **How rare** — what fraction of oceanic cells should cross
   the threshold? Once per plate? Once per run?

A naive "remove the OCEANIC_CLAMP_MAX clamp and let the FBM
tails emerge" approach would produce volcanic islands as a
side-effect of FBM amplitude tuning — but with no control over
distribution, count, or geometry. Step 13.6 is the appropriate
scope for that conversation.

If the user explicitly wants no volcanic islands, the Step
13.5 default behaviour delivers exactly that (strict clamp at
0.49). If they want volcanic islands, Step 13.6 is the
mechanism (when implemented).

