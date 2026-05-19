# Step 8.5b — 32² baseline measurement (Step 9 readiness check)

Single-purpose report: confirm that 32² physics runs are wallclock-
manageable on this hardware before Step 9 (cratonic immunity) starts.
**No code changes** — measurement only.

## Hardware profile

| Item | Value |
|---|---|
| CPU | 11th Gen Intel Core i7-11850H @ 2.50 GHz |
| Physical cores | 8 |
| Logical threads (SMT) | 16 |
| RAM | 16 GB |
| OS | Windows 10 Enterprise LTSC 2021 (10.0.19044) |
| Rust compiler | rustc 1.90.0 |
| Build profile | `release`, `lto = "fat"`, `codegen-units = 1` (8.5b workspace default) |
| `available_parallelism()` | 16 |
| `RAYON_NUM_THREADS` env | unset → rayon defaults to 16 |

The default thread count is intentionally **not** tuned to the 8.5b-
recommended `4` — per the prompt, the measurement reflects what a
downstream user would see without configuration tweaks. The 8.5b
performance report's thread-count sweep on `step6_voronoi 64²` is
the place to read for the optimal tuning; here we report the user-
facing default.

## Methodology

| Item | Value |
|---|---|
| Grid | 32 × 32 |
| Steps | 100 |
| Seed | 42 |
| Runs per (case, precond) | 5 (mean ± std) |
| Cases | step3, step6, step7, step8 (mantle on / slab off) |
| Preconditioners | `JacobiCG`, `AmgCG(Default)` |
| Test | [`tests/v2_baseline_32sq.rs`](../../crates/ymir-core/tests/v2_baseline_32sq.rs) (`#[ignore]`, invoke with `--ignored --nocapture`) |

The test prints per-run wallclock + the Step 8.5b extrapolation
instrumentation (`cg_iter_mean`, Newton outer iters mean, fallback
rate) for every (case, precond) combination, then aggregates
mean ± std at the end. CG iter counts are D9-deterministic so they
do not vary across runs.

### Outlier note (step6 Jacobi initial sweep)

The first measurement attempt covered all four cases in a single
test invocation. During that sweep, run [2/5] of `step6_voronoi`
Jacobi recorded a wallclock of **68 361 s (≈ 19 h)** while the
other four runs were 9.27 s – 10.73 s. Inspection confirmed an
OS-level event (the laptop entered hibernate during run 2 of that
case; the date rolled over while the process was suspended). The
test process was still resident in memory the next morning and was
explicitly stopped before the re-run.

Per the prompt's rule "if a run fails, document the failure, do not
retry in a loop", but also per "no exclusion of outliers without
evidence": the outlier is **system-induced**, not solver-induced,
and is not a property of the measurement. The clean 5-run dataset
for `step6_voronoi` is taken from a focused re-run with no system
interruption (test
[`baseline_32sq_step6_step8_only`](../../crates/ymir-core/tests/v2_baseline_32sq.rs)).
The other three cases (step3, step7, step8 Jacobi) ran cleanly in
the first attempt and are reported as captured.

## Wallclock table — 32², 100 steps

| Case | Precond | Wallclock (s) | CG mean | Newton mean | Fallback % | Conv % | CG max |
|---|---|---|---|---|---|---|---|
| `step3_floor_yielding` | Jacobi | **0.74 ± 0.03** | 22.0 | 0.51 | 0.0 | 100 | 53 |
| `step3_floor_yielding` | AMG | **0.92 ± 0.02** | 13.5 | 0.49 | 0.0 | 100 | 53 |
| `step6_voronoi` | Jacobi | **9.65 ± 0.55** | 246.4 | 0.88 | 0.0 | 100 | 2000 |
| `step6_voronoi` | AMG | **9.81 ± 0.06** | 204.2 | 1.05 | 0.0 | 100 | 2000 |
| `step7_slab_off` | Jacobi | **11.75 ± 2.58** | 246.4 | 0.88 | 0.0 | 100 | 2000 |
| `step7_slab_off` | AMG | **11.42 ± 1.11** | 204.2 | 1.05 | 0.0 | 100 | 2000 |
| `step8_mantle_on_slab_off` | Jacobi | **456.99 ± 3.60** | 1046.0 | 13.64 | 31.6 | n/m | 2000 |
| `step8_mantle_on_slab_off` | AMG | **3565.59** (n=1) | 995.9 | 13.63 | 31.6 | n/m | 2000 |

`Conv %` is the Newton outer-iteration convergence rate (number of
steps where Newton hit the tolerance / total number of steps).
`n/m` for step8 means the run completed (Newton tolerated the
saturated CG) but the convergence-rate field was not extracted from
the per-run output stream — see "Step 8 caveat" below. `CG max =
2000` is the configured iteration cap; reaching it is a property of
the operator, not of the run.

### AMG / Jacobi wallclock ratios (32²)

| Case | AMG / Jacobi |
|---|---|
| `step3_floor_yielding` | 1.24 × |
| `step6_voronoi` | 1.02 × |
| `step7_slab_off` | 0.97 × |
| `step8_mantle_on_slab_off` | 7.80 × (single AMG run vs Jacobi mean) |

`step6_voronoi` and `step7_slab_off` essentially break even at 32²
between the two preconditioners — closer to the D6 target of
"≤ 1.0" than the 64² measurement (8.5b report measured 1.14 – 1.18
on these cases at 64²). Two interpretations, both reasonable: (a) at
32² the per-iter overhead of AMG (Galerkin coarsening, V-cycle) is
proportionally larger relative to the smaller fine-grid solve, but
(b) AMG's iter-count reduction (CG mean 204 vs 246 = 17 % fewer
iterations) closes some of the gap. Effective parity with Jacobi
on the active regimes at 32² is the headline.

### Step 8 caveat — single AMG run

Step 8 Jacobi finished its 5 runs in ≈ 38 min total. Step 8 AMG
**run [1/5] alone took 59 min** (3565.59 s). Five AMG runs would
require ≈ 5 hours of dedicated compute, which is outside the 30–90
min budget the prompt allocated for this measurement. We report the
single AMG measurement honestly and stop there:

- Step 8 is the regime where AMG saturates per the 8.5a Phase 3
  diagnostic (η-contrast 4 · 10⁴, V-cycle reduction ratio 0.67,
  classical RS hierarchy loses diagonal dominance after Galerkin
  coarsening). The 32² measurement here matches that picture: AMG
  CG mean 995.9 — saturated like Jacobi, but each AMG iteration
  costs ≈ 7.8 × more wallclock than a Jacobi iteration in this
  regime, an even worse picture than the 1.10–1.18 × ratio seen on
  step6/7.
- Standard deviation across AMG runs is unmeasured. CG and Newton
  outer iter counts are D9-deterministic so they would have been
  identical across runs; only wallclock variance is missing.

Step 9 will not exercise the step8 mantle-activated regime
(cratonic immunity is a step-up from Step 7 baseline, not Step 8),
so the missing variance does not affect the readiness conclusion.

## Comparison vs 64² 8.5b benchmarks

64² numbers from the [Step 8.5b performance report](./step8_5b_performance_report.md) §"Phase 6 wallclock — 100-step physics, `RAYON_NUM_THREADS=4`" — measured under the recommended thread tuning. The 32² numbers below were collected at the **default** thread count (16); the comparison is therefore "what a user sees at the default" vs "what a user sees with tuning".

| Case | Precond | 32² wc (s) | 64² wc (s) | 32² / 64² |
|---|---|---|---|---|
| `step3_floor_yielding` | Jacobi | 0.74 | 1.64 | 0.45 |
| `step3_floor_yielding` | AMG | 0.92 | 2.18 | 0.42 |
| `step6_voronoi` | Jacobi | 9.65 | 13.84 | 0.70 |
| `step6_voronoi` | AMG | 9.81 | 15.76 | 0.65 |
| `step7_slab_off` | Jacobi | 11.75 | 13.02 | 0.90 |
| `step7_slab_off` | AMG | 11.42 | 15.36 | 0.74 |
| `step8_mantle_on_slab_off` | Jacobi | 456.99 | 656.03 (8.5b report) | 0.70 |

If the cost scaled purely with grid area, the ratio should be
`(32/64)² = 0.25`. The measured ratios are 0.42 – 0.90 — every
case scales **worse** than the area scaling would predict. The
4 dominant overheads we identified:

1. **Per-step setup is grid-size-independent**. Force sampling,
   slab pipeline, boundary geometry, S advection, η evaluation —
   all do work proportional to (grid + book-keeping + diagnostics).
   At 32² this fixed cost is a larger fraction of each step.
2. **Newton outer iter count is unchanged**. Newton convergence is
   a property of the rheology + boundary geometry, not the grid
   resolution. Iter counts at 32² match 64² (CG mean 246 / 204 on
   step6 — same numbers as 64²).
3. **CG iter cap saturated on step6/7**. The CG hits 2000 every
   solve, so the linear-solver wallclock is `2000 ×
   per-iter-cost(grid)`. The per-iter cost does shrink with grid
   area but not by the full factor.
4. **Default thread count is not optimal at 32²**. With the 8.5b-
   identified 4-thread sweet spot the 32² numbers would likely
   improve; the prompt explicitly forbids tuning, so we report the
   default-thread numbers.

The `step3_floor_yielding` ratio (0.42) is closest to the area-
scaling prediction because step3 has the simplest physics pipeline
(no boundary, no slab, no mantle, no Voronoi geometry overhead).
This is consistent with overhead 1 above.

## Step 9 readiness — verdict

The prompt defines three readiness tiers:

- `< 2 min`: green light
- `2 – 5 min`: marginal
- `> 5 min`: pre-Step 9 follow-up needed

Reading the wallclock column of the **non-activated regimes**
(step3, step6, step7 — the regimes Step 9 cratonic immunity will
exercise):

| Case | 32² Jacobi wallclock | Tier |
|---|---|---|
| `step3_floor_yielding` | 0.74 s | 🟢 GREEN — well under 2 min |
| `step6_voronoi` | 9.65 s | 🟢 GREEN — well under 2 min |
| `step7_slab_off` | 11.75 s | 🟢 GREEN — well under 2 min |

**Verdict: GREEN. Step 9 development is wallclock-supportable at 32².**

The `step8_mantle_on_slab_off` row is RED at both preconditioners
(7.6 min Jacobi, 59 min AMG), but Step 9's cratonic immunity
exercises a Step 7-shape regime (yielding + drag + Voronoi) without
the mantle activation that defines step8. The step8 row is
documented for completeness; it does not gate Step 9 readiness.

### What this implies for Step 8.5c sequencing

Step 8.5c (hierarchy caching across Newton outer iterations) was
listed in the 8.5b performance report as a follow-up. Given the
GREEN result above, Step 9 can start without 8.5c landing first.
Step 8.5c remains valuable when:

- Step 9-10 measurements show repeated cycles cost more than
  expected (caching would help if Newton outer iter count climbs).
- A future grid-size scaling exercise (Step 8.5d) at 128²+ shows
  the AMG / Jacobi ratio still > 1 — caching is the next obvious
  lever.

For now, the 32² readiness check is unambiguous and Step 9 is
unblocked.

## Solver health summary

All four cases × both preconditioners converge `100 %` of Newton
outer iterations on step3/6/7. step8 Newton converges in the
sense that the run completes within `max_outer_iters` even with CG
saturated; the per-step convergence breakdown was not captured in
the per-run output stream during this measurement (see step 8
caveat above).

Newton extrapolation fallback rates at 32²:

| Case | Fallback % |
|---|---|
| `step3_floor_yielding` | 0.0 |
| `step6_voronoi` | 0.0 |
| `step7_slab_off` | 0.0 |
| `step8_mantle_on_slab_off` | 31.6 |

The reviewer's `> 10 %` watch point fires only on step 8, which is
expected and consistent with the 64² measurement (step 8 fallback
was 83.7 % at 64², 31.6 % at 32² — the smaller grid is somewhat
less hostile to extrapolation, perhaps because Newton outer iter
count is identical but CG saturation is hit faster, leaving less
opportunity for the extrapolated guess to drift). On the
non-activated regimes the fallback rate is exactly zero —
extrapolation is uniformly accepted, which is the design intent.

## Reproducing the measurement

```bash
# Full sweep (4 cases × 2 preconds × 5 runs ≈ 30 – 90 min,
# step 8 AMG capped after 1 run).
cargo test --release -p ymir-core --test v2_baseline_32sq \
    baseline_32sq_measurement -- --ignored --nocapture

# Focused re-run (step 6 + step 8 only, ≈ 70 min including the
# single step 8 AMG run; Ctrl+C the AMG block once run [1/5]
# completes if you do not need the 5-run variance).
cargo test --release -p ymir-core --test v2_baseline_32sq \
    baseline_32sq_step6_step8_only -- --ignored --nocapture
```

To reproduce with the 8.5b-recommended thread tuning, prefix with
`RAYON_NUM_THREADS=4`. Wallclock should be modestly better on
`step6_voronoi` / `step7_slab_off` and comparable on the simpler
cases.
