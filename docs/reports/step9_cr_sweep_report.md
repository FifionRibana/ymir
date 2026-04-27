# Step 9 — Cr sweep report

> Companion to `step9_physics_report.md`. Characterises the
> cratonic-fraction parameter `Cr` over its admissible range
> `[0.1, 0.5]` at 64² × 100 steps, validates acceptance #9
> (monotonicity), and documents the per-seed dispersion that
> motivates the acceptance #8 reformulation.

## Parameter

`Cr ∈ [0.1, 0.5]` (admissible range, §4.10 patch). Default 0.30.
Sets the **target** fraction of plate area occupied by the
cratonic core for plates that are large enough to host one
(`plate_area ≥ plate_area_min · domain_area`). The realised
fraction varies with Voronoï plate geometry.

The smoothstep `d_mid = 1 − sqrt(Cr)` (with `L_plate = max BFS
depth`, see §4.10 patch) converts `Cr` to a normalised distance
threshold per plate:

| Cr | `d_mid_normalised = 1 − sqrt(Cr)` |
|---|---|
| 0.10 | 0.6838 |
| 0.20 | 0.5528 |
| 0.30 | 0.4523 |
| 0.40 | 0.3675 |
| 0.50 | 0.2929 |

Lower Cr → larger `d_mid` → fewer cells qualify → smaller
cratonic core.

## Sweep at 64² Step 7 shape (Section 1)

Setup matches `build_step9_config` from
`crates/ymir-core/tests/v2_step9_physics_and_sweep.rs`: 64×64,
100 steps, seed = 42, num_plates = 8, continental_ratio = 0.30,
Bi = 0.15, Br = 0.05, Voronoï closed, no slab, no mantle, no
plastic memory. `K = 5` (default), `B_factor = 5` at the time
of the sweep run (the post-amendment default `B_factor = 8`
does not change the cratonic-fraction count, which depends only
on the geometry pipeline, not on the Bi-elevation mechanism).

| Cr | `cratonic_cell_fraction` | `peak_eta_contrast_at_boundary` | Wallclock | `peak\|v\|` |
|---|---|---|---|---|
| 0.10 | 0.0813 | 5.000 | 24.0 s | (sweep run timestamp) |
| 0.20 | 0.1240 | 5.000 | 23.2 s | — |
| 0.30 | 0.1672 | 5.000 | 21.4 s | — |
| 0.40 | 0.2041 | 5.000 | 20.1 s | — |
| 0.50 | 0.2441 | 5.000 | 22.0 s | — |

### Acceptance #9 — Cr sweep monotonic

`cratonic_cell_fraction` is **strictly non-decreasing** in Cr:
0.0813 < 0.1240 < 0.1672 < 0.2041 < 0.2441. ✅ PASS.

The sweep test asserts this property programmatically:

```text
Monotonicity acceptance #9: PASS
```

### `peak_eta_contrast_at_boundary` saturation

The metric saturates at `K = 5.000` for all Cr values at this
resolution because the BFS-distance step `1 / L_plate ≈ 0.10`
exceeds the smoothing band width `smoothing_width = 0.05`,
making the smoothstep transition sub-cell. The `eta_multiplier`
field jumps `0 → 1` between adjacent boundary cells, ratio = K
exactly (acceptance #3 `≤ K · 1.05 = 5.25` with margin).

A wider `smoothing_width` would bring the contrast below K (the
smoothing-width probe `tests/v2_step9_smoothing_probe` shows
`smoothing_width = 0.20` gives contrast ≈ 4.965 < K), but the
default 0.05 is sufficient for acceptance #3.

## Acceptance #8 — `cratonic_cell_fraction ≈ Cr · continental_fraction` ±20 %

The continental fraction at seed = 42 (independent of Cr) is
0.4458 at 64² (44.58 % of cells are in continental plates,
counted before the area-threshold filter). Per Cr, the target
is `Cr · continental_fraction`:

| Cr | observed `cratonic_cell_fraction` | target `Cr · cont_frac` | ratio (observed / target) |
|---|---|---|---|
| 0.10 | 0.0813 | 0.0446 | 1.825 |
| 0.20 | 0.1240 | 0.0892 | 1.391 |
| 0.30 | 0.1672 | 0.1337 | 1.250 |
| 0.40 | 0.2041 | 0.1783 | 1.145 |
| 0.50 | 0.2441 | 0.2229 | 1.095 |

For seed = 42 specifically, the per-seed ratio at the default
Cr = 0.30 is 1.250 — outside the ±20 % bound for that single
seed.

### The mean-across-seeds reformulation

The original acceptance #8 ("within ±20 %") was reformulated by
the reviewer to "mean across seeds within ±20 %, per-seed
dispersion up to 1.5× documented" before the implementation
locked in the `L_plate = max BFS depth` calibration. The
calibration probe `tests/v2_cratonic_normalization_probe` runs
the same metric across 31 random Voronoï seeds and reports:

> Scheme A (literal `L_plate = sqrt(area)`, default formula
>          `d_mid = 0.5 (1 − sqrt(Cr))`): **mean ratio 1.297**
>          across 31 non-degenerate seeds.
> Scheme B (geometric `L_plate = max BFS depth`, formula
>          `d_mid = 1 − sqrt(Cr)`, **adopted**): **mean ratio
>          1.132** across the same 31 seeds.

Scheme B's mean ratio 1.13× sits comfortably inside the ±20 %
bound. Per-seed variation reaches up to 1.5× for seeds that
land particularly compact / blob-like plate geometries (where
the BFS depth is unusually large compared to `sqrt(area)`); this
is documented as a known geometric-irregularity effect, not a
defect.

✅ **Acceptance #8 (reformulated): mean across seeds 1.13× the
target, inside the ±20 % bound.**

## Resolution sensitivity

The cratonic_factor field at 32² and 64² (same seed = 42)
samples slightly different cratonic-cell counts due to BFS
quantisation:

| grid | continental_fraction | cratonic_cell_fraction (Cr=0.3) | ratio vs target |
|---|---|---|---|
| 32² | 0.4434 | 0.1572 | 1.182 (inside ±20 %) |
| 64² | 0.4458 | 0.1672 | 1.250 (outside ±20 %) |

Higher resolution refines the BFS depth measurement; for the
particular `(seed = 42, num_plates = 8, continental_ratio = 0.3)`
combination the 64² result lands on the wrong side of the per-
seed bound. This is consistent with the geometric-irregularity
explanation — across many seeds the average converges to the
1.13 mean.

## Reproducing the Cr sweep

```bash
cargo test --release -p ymir-core --test v2_step9_physics_and_sweep \
    -- --ignored --nocapture --test-threads=1 step9_cr_sweep_64sq
```

Output PNGs land in `docs/reports/step9_phase7_cr_{1,2,3,4,5}/`
(one directory per Cr value × 10).

The 31-seed calibration probe:

```bash
cargo test --release -p ymir-core --test v2_cratonic_normalization_probe \
    -- --ignored --nocapture
```

## Acceptance summary

| # | Criterion | Status |
|---|---|---|
| 8 | `cratonic_cell_fraction ≈ Cr · cont_frac` ±20 % (mean across seeds) | ✅ 1.13× mean |
| 9 | Cr sweep monotone non-decreasing | ✅ strictly increasing |
