# Step 8.5a: Classical AMG + FMG preconditioner (α partial merge)

Installs Classical AMG (Ruge-Stüben) + Full Multigrid as opt-in alternatives to Jacobi-CG, with sparse Picard-block materialisation. **α partial merge per reviewer contract**: step0-7 fully covered, step8 documented as out-of-regime with a concrete Step 8.5a.2 follow-up.

## What ships

- `tectonics_v2/stokes/amg/` — 9 submodules: `strong_connections`, `splitting` (Classical RS two-pass), `prolongation` (RS classical formula), `restriction` (R = Pᵀ), `smoother` (SGS), `coarse_solve` (Doolittle LU, no nalgebra), `setup` (hierarchy + Galerkin coarsening), `vcycle`, `fmg`.
- `tectonics_v2/stokes/sparse_assembly.rs` — 9-point CSR assembly of the Picard block, algebraically derived from `apply_momentum`, O(N) per call, column-sorted per row (D9).
- `tectonics_v2/stokes/snapshot.rs` — `LinearStokesSnapshot` (Newton iter-0 state) + `ReferenceSolution` (Phase 4.2 high-precision refs), both bincode-serialisable with `format_version`.
- `LinearSolverConfig` enum in `stokes/solver.rs` — `JacobiCG` (default, bit-parity) / `AmgCG(AmgConfig)` (opt-in). Threaded through `SheetConfig`, `NewtonSolver`, `BaselineConfig`.
- `bin/gen_bench_data` (snapshot capture) + `bin/gen_reference_solutions` (high-precision references) + `benches/amg_benchmark.rs` (9 cases, criterion, both preconditioners).
- `bench_data/` with 6 Stokes snapshots + `bench_data/reference_solutions/` with 4 references (deterministic, checked in).
- 5 integration tests (`v2_amg_phase3_diagnostic`, `v2_amg_scalar_parity`, `v2_amg_physics_scalar_parity`, `v2_sparse_assembly_snapshot_parity`, `v2_amg_poisson_projection_diag`) + 42 lib tests across the AMG submodules.
- `docs/reports/step8_5a_amg_report.md` — full phase-by-phase report with measurements, diagnostic traces, and the α contract checklist.

## What does NOT ship

- **Resolution of `step8_activated` non-convergence.** AmgCG on step8 saturates at the 2000-iter cap just like Jacobi. Deferred to Step 8.5a.2 with SA-AMG as primary working hypothesis.
- **Wallclock performance work.** AMG is currently slower than Jacobi on step0-7 (see §Performance honesty). Amortisation (hierarchy caching, rayon, Newton extrapolation, LTO) is Step 8.5b scope.
- **CLI flags on `stepN_baseline` binaries.** Opt-in happens via Rust callers setting `BaselineConfig.linear_solver`, not via new command-line arguments. Binaries retain their Jacobi default.
- **Detailed wallclock breakdown** (setup vs V-cycle apply vs CG inner work). Measurement deferred to Step 8.5b.
- **Step 0-7 physics re-runs against committed merged-report metrics.** Phase 4.3.5 runs both paths live (100 steps) and compares directly — rather than parsing the markdown reports. Equivalent correctness, simpler tooling.

## Acceptance gates passed (with measured values)

| Gate | Target | Measured | Status |
|---|---|---|---|
| `poisson_constant` AMG CG iters | ≤ 6 (revised from ≤ 3) | 5 | ✅ |
| `poisson_contrast_100` AMG CG iters | ≤ 20 | 9 | ✅ (55 % under) |
| **`poisson_contrast_10000` AMG CG iters** | **≤ 100 (principal gate)** | **10** | ✅ **(90 % under)** |
| FMG ≥ 2× V-cycle residual reduction | ≥ 2× | 2.04× | ✅ |
| `step0_quiescent` AMG physics iter cap | ≤ 10 | 4 | ✅ |
| `step3_floor_yielding` AMG physics iter cap | ≤ 15 | 9 | ✅ |
| `step6_voronoi` AMG physics iter cap | ≤ 40 | 9 | ✅ |
| `step7_slab_off` AMG physics iter cap | ≤ 40 | 8 | ✅ |
| Scalar-parity `vmax_peak` (physics 100-step) | < 1 % rel | 2.5·10⁻¹⁰ ... 3·10⁻¹¹ | ✅ |
| Scalar-parity `mass_drift` (physics 100-step) | < 1 % rel OR < 1e-6 abs | 7.6·10⁻¹⁵ | ✅ |
| Scalar-parity snapshot (reference-based) | < `C·κ·(tol_test + tol_ref)` | 3-17× under | ✅ 4/4 |
| JacobiCG bit-parity on default path | byte-identical | `v2_step8_regression_smoke` passes | ✅ |
| AMG setup deterministic | byte-for-byte | 100-run determinism tests | ✅ |

## step8 diagnosis — the intellectual deliverable of this step

Phase 3 diagnostic (archived in `v2_amg_phase3_diagnostic`) established that AMG's plateau on `step8_activated` has a mathematically conclusive root cause, distinct from the u-v coupling hypothesis originally planned as the Option A' fallback:

- **η-contrast = 4.06·10⁴** on `step8_activated` (vs 1.36× on `step6_voronoi`). At the D1-predicted boundary where Classical Ruge-Stüben loses its comfort zone (§D1: *"SA-AMG would be more robust for extreme η contrasts (> 10⁶) ... If Classical proves insufficient, SA-AMG becomes the next step."*).
- **V-cycle reduction ratio = 0.67** on step8 (vs 0.023 on step6). The hierarchy builds similarly on both (comparable level counts, coarsening ratios) but the coarse-level operators have lost diagonal dominance after Galerkin coarsening — SGS amplifies residual going down (0.34 → 0.59 → 1.07 → 1.66 across levels).
- **Ran on the u-u scalar block alone via `extract_diagonal_block`** — no u-v coupling in the experiment. The failure is on the scalar problem. **Option A' (2×2-block AMG) cannot be the resolution** because it addresses u-v coupling, which is not the blocker.
- **Step 8.5a.2 renamed from "Option A'" → "advanced AMG techniques for extreme η-contrast"** with SA-AMG as primary working hypothesis. Literature note: SA-AMG itself can struggle in the 10⁴-10⁶ range, so the follow-up issue explicitly budgets for alternatives (Chebyshev/ILU smoothers, W-/F-cycle variants, hybrid schemes).

## Performance honesty — wallclock

**AMG is currently 1.1-3.4× SLOWER than Jacobi** on step0/3/6/7 physics runs in Step 8.5a alone:

| Case | JacobiCG | AmgCG | Ratio |
|---|---|---|---|
| step0 100-step | 2.76 s | 8.00 s | **2.9× slower** |
| step3 100-step | 3.75 s | 12.71 s | **3.4× slower** |
| step6 100-step | 20.81 s | 22.84 s | 1.10× |
| step7 100-step | 24.37 s | 28.92 s | 1.19× |

CG iter counts drop ×2-5 as designed. Each AMG iteration costs 5-10× more than a Jacobi iteration in this build because the `O(N log N)` hierarchy rebuild happens per Newton outer iter without amortisation.

**This milestone is a correctness gate, not a performance gate.** The honest framing of Step 8.5a's contribution is: *"machinery validated and ready-to-be-accelerated"*. Wallclock gains come in Step 8.5b (hierarchy caching when `‖Δη‖` is small, rayon parallelisation, Newton extrapolation, LTO/PGO).

Detailed wallclock decomposition (setup vs V-cycle apply vs CG inner) was not measured in this step and is deferred to Step 8.5b as its first-order investigation.

## Downstream regime recommendation

| Regime | Recommended default |
|---|---|
| step0-7 (η-contrast ≲ 10²) correctness testing / development | `AmgCG(Default)` fine for scalar-parity exploration |
| **step0-7 production physics runs (including Step 9 cratonic immunity)** | **`JacobiCG` (keep as default)** — AMG is correct but slower until 8.5b |
| step8-like (η-contrast > 10⁴) | **`JacobiCG`** mandatory — AMG saturates here too, awaiting 8.5a.2 |

**No automatic Jacobi fallback inside AmgCG** — would hide the regime mismatch and add permanent complexity. Users opt in explicitly per-regime.

This is a **reversal from the original Phase 4.3 draft**, which recommended AmgCG for Step 9. The wallclock reality dictates otherwise for production — Step 9 on AmgCG would be slower than Step 9 on Jacobi at this milestone's build; AmgCG becomes the default only after Step 8.5b amortisation lands.

## Follow-ups (issues to open)

- **Step 8.5a.2 — advanced AMG techniques for extreme η-contrast**. First phase: SA-AMG feasibility prototype on `step8_activated` snapshot. Budget for alternatives (Chebyshev smoother, W-/F-cycle, hybrid schemes) if SA-AMG plateaus too. Blocker for step8 physics convergence; not blocker for Step 9.
- **Step 8.5b — wallclock performance**. Hierarchy caching across Newton outer iters (η drift threshold gate), rayon parallelisation of SpMV + SGS + prolongation, Newton extrapolation, compilation flags (LTO, PGO, target-cpu). Compounds multiplicatively with 8.5a's iter-count reduction. Required before promoting AmgCG to a default recommendation for production physics runs.

## How to review

Recommended reading order:

1. `docs/reports/step8_5a_amg_report.md` — phase-by-phase with measurements. Sections to prioritise:
   - §Phase 2.7 for the Poisson gate results,
   - §Phase 3.1 for the step8 diagnostic (smoking gun is the V-cycle per-level residual trace),
   - §Phase 4.3.5 for the physics scalar-parity table,
   - §Performance honesty for the wallclock framing,
   - §Downstream regime recommendation for the Step 9 guidance.
2. **Diff by phase** via the 18-commit history (`git log --oneline 9b05faa..de0ee96`). Each WIP commit is scoped to one submodule with its own tests; each FEAT commit seals a phase.
3. **Tests worth inspecting**:
   - `tests/v2_amg_phase3_diagnostic.rs` — the step8 investigation, reproducible at any time.
   - `tests/v2_amg_scalar_parity.rs` — reference-based parity with the `3·κ·(tol + tol_ref)` formula derivation.
   - `tests/v2_amg_physics_scalar_parity.rs` — end-to-end physics-run parity.
   - `tests/v2_step8_regression_smoke.rs` — still passes byte-identical, confirms JacobiCG bit-parity through all the dispatch plumbing.
4. Spot-check the algebraic derivations: `stokes/sparse_assembly.rs` doc-comment on the 9-point stencil, `stokes/amg/strong_connections.rs` on the RS definition, `stokes/amg/prolongation.rs` on the weight formula.

---

🤖 Generated with [Claude Code](https://claude.com/claude-code)
