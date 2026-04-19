# Slab-Pull Reformulation — Phase 1-bis Diagnostic

Issue #75 — *Reformulate slab-pull as an auto-regulated operator term instead of RHS forcing.*
Branch: `75-reformulate-slab-pull-...` (same branch as Phase 1). Date: 2026-04-19.

This report answers the empirical questions left open by the Phase 1
diagnostic: **where** the slab-pull cost actually goes, **whether** the
RHS is really spiked, and **what** intervention is most promising. Four
reference scenarios on 64²/seed 42/300 steps, three reps each.
Instrumentation lives in `tectonics/solver/diagnostics.rs` and routes
tracing debug events through four dedicated targets.

**Top-line result.** B vs C wallclock ratio is **4.75×**, not 1.5× as
the parent task estimated. Slab-pull is the dominant cost driver, but
*via an η-contrast cascade*, not via RHS spiking. T_plates is not
spiked at all. The Phase 2 refactor is justified, but for a different
reason than the issue originally gave.

---

## 1. Scenario wallclock comparison

300 macro steps, adaptive dt, Newton solver, `dt_target = 2.0`.

| Scenario | Description               | mean (s) | min (s) | max (s) | Δ vs A     |
|----------|---------------------------|---------:|--------:|--------:|-----------:|
| A        | bare thin-sheet           |     10.6 |    10.5 |    10.9 | —          |
| B        | all extensions on         |    841.0 |   836.0 |   846.5 | +830 s     |
| C        | all on except slab-pull   |    177.1 |   173.2 |   179.3 | +167 s     |
| D        | all on except mantle flow |    851.6 |   792.7 |   891.1 | +841 s     |

Derived:

- **Slab-pull premium (B − C):** 663.9 s per 300-step run = **79% of B's wallclock.**
- **Mantle premium (B − D):** −10.6 s. Mantle is free (within noise; D is even slightly slower).
- **Cost of everything *except* slab-pull/mantle (C − A):** 166.5 s.

Variance across 3 reps is ≤ 1% for A/B/C and 6% for D; the ordering is robust.

---

## 2. Per-phase wallclock breakdown

Mean per-sub-step elapsed (microseconds) from `phase_timings`, steps
100–200 (steady regime), averaged across 3 reps.

| Phase        | A      | B          | C         | D          |
|--------------|-------:|-----------:|----------:|-----------:|
| boundaries   |      0 |        285 |       237 |        245 |
| **solve**    | 20,524 |  1,372,010 |   630,639 |  1,317,110 |
| advection    |     90 |        119 |       113 |        113 |
| recycling    |      0 |         35 |        30 |         30 |
| plates       |      0 |        377 |       329 |        309 |

**Solve dominates in every scenario (≥ 99.95%).** Per-substep:

- B → C saves **741 ms/step** in the solve phase alone — a 54% reduction.
- B → D saves **55 ms/step** — 4% reduction. Mantle flow adds almost nothing to the linear solve.
- Boundaries, advection, recycling, plates together account for < 0.1% of wallclock in B. Optimising them is pointless.

Per-substep t_solve times (B = 1.37 ms vs C = 0.63 ms) differ by 2.2×,
but the wallclock ratio (841/177 = 4.75×) is larger. The missing
factor is **Newton iteration count and adaptive sub-step multiplicity**:
B's η contrast forces both more Newton iters and more sub-steps per
macro step. Both are symptoms of the same ill-conditioning.

---

## 3. RHS spike analysis

From `rhs_breakdown`, steps 100–200, 3 reps merged. `spike_p95 =
max_abs / p95` was used instead of `max_abs / p50` because periodic
Stokes has a rank-2 null space projected out of the RHS, pushing p50
to near zero and making the p50-based ratio degenerate (hits
`DENOM_FLOOR = 1e-20`). `max_abs / p95` stays bounded and expresses
how far the peak exceeds "typical large values."

| Scenario | `gpe_rhs_spike_p95` | `tplates_rhs_spike_p95` | `gpe_max_abs` | `tp_max_abs` | `gpe_norm` | `tp_norm` |
|----------|--------------------:|------------------------:|--------------:|-------------:|-----------:|----------:|
| A        |                8.93 |                    1.39 |         12.15 |         2.14 |       86.4 |      88.3 |
| B        |                3.64 |                    1.18 |         15.05 |         5.92 |      176.2 |     312.7 |
| C        |                2.18 |                    1.46 |          5.12 |         2.84 |       89.0 |      92.7 |
| D        |                3.77 |                    1.01 |         14.88 |         4.99 |      171.2 |     311.9 |

**Which term dominates the spike?** GPE, unambiguously. `gpe_spike_p95`
is 2.2–8.9 across scenarios; `tplates_spike_p95` never exceeds 1.5 —
T_plates is essentially flat. In B, `gpe_max_abs = 15` is 2.5× larger
than `tp_max_abs = 5.9`.

**Which term dominates the L2 norm?** T_plates in B/D (312 vs 176),
GPE in A (86 vs 88, tied) and C (89 vs 92). The slab-pull velocity
boost doubles `tp_norm` from A to B without creating any spike — it
rescales a smooth-piecewise-constant field.

**Key finding.** The Phase 1 report's hypothesis stands: the RHS spike
lives in the GPE gradient across thin-oceanic / thick-continental
thickness jumps. T_plates is never the spike. Moving slab-pull from
RHS to operator does not fix a problem that exists in the RHS.

---

## 4. Viscosity distribution analysis

From `eta_breakdown`, steps 100–200, 3 reps merged. One event per
Newton outer iteration.

| Scenario | `eta_ratio` (mean) | `eta_ratio` (max) | `yielding_cells_fraction` | `saturated_cells_count` |
|----------|-------------------:|------------------:|--------------------------:|------------------------:|
| A        |                2.5 |               2.6 |                      0.00 |                       0 |
| B        |               61.8 |              73.2 |                      1.00 |                       0 |
| C        |               11.3 |              11.6 |                      1.00 |                       0 |
| D        |               57.1 |              66.7 |                      1.00 |                       0 |

**The cascade is real and measurable.** Disabling slab-pull (B→C)
drops `eta_ratio` from 62 to 11 — a **5.5× reduction in viscosity
contrast.** This is the mechanism behind the solve-time reduction:

```
slab-pull boost → plate velocities grow → strain rates grow →
  → η drops toward η_min on yielded cells (and possibly toward
    η_max on cratonic cells that were barely saturating before)
    → eta_ratio blows up → BiCGSTAB condition number blows up
    → more linear iters per Newton, more Newton iters per solve.
```

**Caveat on `yielding_cells_fraction = 1.00`.** The metric as defined
in the task fires whenever `eta_final < 1.01 × eta_plastic`. Because
`apply_yielding` uses `soft_min_harmonic(eta_visc, eta_plastic)`, the
post-yielding value is **always** ≤ `eta_plastic` wherever strain rate
is above `1e-20` (which is most of the domain). So the criterion
captures "yielding was active in the pipeline" rather than "the
plastic branch dominated." The informative signal is `eta_ratio`, not
`yield_frac`. Saturation never fires because `eta_max = 1e4` is never
reached in these runs.

Mantle flow (D) perturbs `eta_ratio` only marginally (62 → 57), which
lines up with its near-zero wallclock effect.

---

## 5. Residual localization analysis

From `residual_spatial`, steps 100–200, 3 reps merged. `F(v_converged)`
comes from the last Newton iteration's `ws.jfnk_f_v`. Boundary cells
are everything classified as `!= BoundaryType::None`.

| Scenario | `residual_localization` (mean) | max  | boundary cell fraction |
|----------|-------------------------------:|-----:|-----------------------:|
| A        |                          0.000 | 0.00 |                    n/a |
| B        |                          0.277 | 0.74 |                  ~0.15 |
| C        |                          0.478 | 0.99 |                  ~0.15 |
| D        |                          0.264 | 0.69 |                  ~0.15 |

(A reports 0 because `boundaries.enabled = false` disables the whole
classification — no boundary mask, so no localization computed.)

**Interpretation.** The boundary cells hold ~15% of the grid but
**27–48% of the residual L2 energy.** Concentration factor is 1.8×
(B), 3.2× (C), 1.8× (D). Residual is mildly concentrated at
boundaries, not overwhelmingly so — the remaining 52–72% lives in the
plate interior.

**Counterintuitive observation.** C (no slab-pull) has *higher*
localization (48%) than B (28%). Interpretation: slab-pull's η
cascade spreads the residual over the whole domain (ill-conditioned
interior), whereas without slab-pull the η field is flatter and the
remaining residual concentrates at the only places where the RHS is
spiked — the boundaries where GPE gradient dominates. This is
consistent with the RHS breakdown in §3.

---

## 6. Hypothesis validation

1. **Is the RHS actually spiked?** *Yes, but only the GPE component.*
   `gpe_spike_p95` = 3.6 (B) and `gpe_max_abs` = 15. T_plates is
   never spiked (`tp_spike_p95` ≤ 1.5 in every scenario). The issue's
   original claim that slab-pull spikes T_plates is empirically
   falsified.

2. **Does slab-pull trigger widespread yielding?** *Yielding is "on"
   everywhere in B, C, D indiscriminately* — the criterion as
   specified does not discriminate. But **η contrast jumps 5.5×**
   (62 vs 11) between B and C. Slab-pull drives ill-conditioning
   through elevated strain rates, not through enabling yielding per
   se. The cascade hypothesis is confirmed by `eta_ratio`, which is
   the meaningful signal.

3. **Where does the slab-pull cost go?** Not into RHS assembly. Into
   the linear solve, via η contrast. BiCGSTAB effective iteration
   count grows with condition number, and the condition number of
   the Stokes operator at fixed grid resolution grows roughly with
   `eta_max / eta_min`. Going from 11 to 62 roughly predicts a
   factor ~5 more BiCGSTAB iterations per solve, which matches the
   observed wallclock ratio.

4. **Where does the mantle flow cost go?** Nowhere. It adds < 1% to
   wallclock and < 10% to `eta_ratio`. The task's estimate of
   "46-second mantle flow cost" is not reproducible at this config
   — mantle is effectively free.

5. **Most promising intervention?**

   | Option | Expected impact | Effort | Rationale |
   |--------|-----------------|--------|-----------|
   | (a) slab-pull → operator | **High** — breaks the η cascade at its source | Medium | γ_slab·(v·n̂)·n̂ is auto-regulating: it bounds v locally at the margin, preventing the velocity-boost-driven strain-rate growth that drives η contrast |
   | (b) smooth GPE gradient | Marginal — affects only max_abs, not the η cascade; distorts ridge-push | Low | GPE spike is physical, not a bug |
   | (c) regularize yielding | Medium — flattens the plastic branch, possibly caps `eta_ratio`; but yielding is physical | Medium | Addresses a symptom |
   | (d) multigrid preconditioner | **High** — addresses conditioning regardless of cause; future-proofs for 256²/512² | High | Already on the block 2 roadmap |
   | (e) do nothing | Not viable — 14 min per 300-step run at 64² is already over budget; 512² will be prohibitive | — | — |

   **(a) and (d) are complementary.** (a) removes the feedback loop
   from slab-pull; (d) absorbs whatever conditioning the other
   extensions still produce.

---

## 7. Recommendation for Phase 2

**Proceed with Phase 2 as scoped, with a revised justification.** The
reformulation `γ_slab · (v·n̂) · n̂` is the correct next step not
because it removes a T_plates spike (there is none), but because it
removes the *physical mechanism* that drives the η-contrast cascade:
the unbounded per-plate velocity boost in `apply_slab_pull`. Replacing
it with a bounded, velocity-coupled operator term caps subduction
velocity locally, which flattens strain rates, which flattens η,
which improves BiCGSTAB conditioning — the observed 4.75× wallclock
slowdown should largely disappear. The GPE spike will remain and
should be addressed separately (or left alone, since it represents
real physics).

After Phase 2 lands, re-run this diagnostic. The expected signature
of success is `eta_ratio` in the B regime dropping from 62 toward
C's 11, and B wallclock dropping from 841 s toward C's 177 s. If
`eta_ratio` stays high, the operator term is not actually replacing
the velocity boost and the refactor needs to be revisited.

---

## Appendix A — Instrumentation

New module: [`crates/ymir-core/src/tectonics/solver/diagnostics.rs`](crates/ymir-core/src/tectonics/solver/diagnostics.rs)
— four public helpers (`emit_rhs_breakdown`, `emit_eta_breakdown`,
`emit_residual_spatial`) plus a private `scalar_dist` that computes L2
norm, max, and percentiles.

Emit sites:

- `stokes.rs::compute_rhs` end → `rhs_breakdown` (once per call; Newton
  calls it multiple times per macro step).
- `newton.rs::solve_velocity_newton` after `compute_nonlinear_residual`
  → `eta_breakdown` (once per Newton outer iter).
- `tectonics.rs::execute_tectonic_pass` after successful Newton →
  `residual_spatial`; at phase boundaries → `phase_timings` (once per
  sub-step).

Expensive per-field percentile computations are guarded with
`tracing::enabled!(target, Level::DEBUG)`, so the production solver
pays nothing when the targets are filtered out. The
`phase_timings` info log always fires — it is two `Instant::elapsed()`
subtractions per phase, negligible.

Scenario runner: [`crates/ymir-core/examples/phase1bis_scenarios.rs`](crates/ymir-core/examples/phase1bis_scenarios.rs).
Aggregator: [`scripts/phase1bis_aggregate.sh`](scripts/phase1bis_aggregate.sh). Logs in
`logs/phase1bis_<S>_<NN>.log`; per-scenario wallclock in
`logs/summary.txt`.

## Appendix B — Reproducing

```bash
cargo build --release --example phase1bis_scenarios
for s in A B C D; do
  for r in 1 2 3; do
    ./target/release/examples/phase1bis_scenarios.exe "$s" "$r" logs \
      >> logs/summary.txt
  done
done
./scripts/phase1bis_aggregate.sh
```

Total wallclock on the reference machine: 4 × 3 runs × (A=10 / B=840 /
C=177 / D=850 s) ≈ **59 minutes**. A single rep of each is enough for
the ordering; 3 reps are needed only for the ± variance bounds in §1.

## Appendix C — Caveats and known instrumentation gaps

1. `yielding_cells_fraction` criterion (`eta_final < 1.01 × eta_plastic`) fires
   too broadly because of the soft-min blend — saturates at 1.0 for all
   scenarios with yielding enabled. `eta_ratio` is the informative signal.
2. `spike_ratio` using p50 as denominator is degenerate due to null-space
   projection pushing the median to zero. The report uses `max_abs / p95`
   instead, which is well-defined.
3. `emit_residual_spatial` is only meaningful for Newton (Picard does not
   populate `ws.jfnk_f_v`). The scenario runner forces `NonlinearSolver::Newton`.
4. `phase_timings` clocks each call to `execute_tectonic_pass`, which is
   once per adaptive *sub-step*, not once per macro step. With B averaging
   ~2 sub-steps per macro step (inferable from `841 s / 300 steps /
   1.37 ms per solve ≈ 2.04`), per-sub-step times under-report the
   macro-step wallclock. Aggregating over 100 macro steps smooths this
   out.
5. Wallclocks are ~13× larger than the task description's reference (65 s
   for B at 64²/300 steps). Likely explanations: continuation enabled by
   default (n_steps = [1.0, 1.5, 2.0, 2.5, 3.0] cold-starts the first
   step with 5 successive solves), Newton as solver (vs Picard default),
   different defaults from whatever config produced the 65 s reference.
   The *relative* ordering of A/B/C/D is the load-bearing finding and is
   unaffected by the absolute scale.
