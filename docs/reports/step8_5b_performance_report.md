# Step 8.5b — Performance (rayon + Newton extrapolation + LTO)

**Status**: Phase 0 landed (hardware profile captured, LTO enabled).
Subsequent phases will extend this report inline.

**Branch**: `99-step-85b-performance-rayon-newton-extrapolation-lto`
**Target**: `milestone/solver-reconstruction`
**Prior step**: [Step 8.5a — Classical AMG + FMG](./step8_5a_amg_report.md)
**Motivation**: §Performance honesty of Step 8.5a report documented
that AMG is 1.1–3.4× slower than Jacobi on step0–7 despite iter-count
reductions of ×10–24. Step 8.5b installs the three compound
accelerators (rayon, Newton extrapolation, LTO) that amortise the
AMG machinery into an actual wallclock win.

## Hardware profile (D7)

| Item | Value |
|---|---|
| CPU | 11th Gen Intel Core i7-11850H @ 2.50 GHz |
| Physical cores | 8 |
| Logical threads (SMT) | 16 |
| RAM | 16 GB (16 325 677 056 bytes) |
| OS | Windows 10 Enterprise LTSC 2021 (10.0.19044) |
| Rust compiler | rustc 1.90.0 (1159e78c4 2025-09-14) |
| Cargo | 1.90.0 (840b83a10 2025-07-30) |
| Build profile | `release`, `lto = "fat"`, `codegen-units = 1` |

All wallclock numbers in later phases are meaningful **only relative
to this hardware**. Cross-machine portability is explicitly not a
property of the reported gains.

## Phase 0 — LTO enabled, report skeleton

`Cargo.toml` now carries:

```toml
[profile.release]
lto = "fat"
codegen-units = 1
```

Cold release-build time penalty: ~20–60 s on this hardware (measured
in Phase 0 smoke test, exact number recorded below). Debug/test
profiles unchanged. `target-cpu=native` and PGO deliberately excluded
per D4 (portability + ceremony/gain ratio).

### Phase 0 smoke — LTO compiles cleanly

| Measurement | Value |
|---|---|
| Cold release build (`cargo build --release -p ymir-core`) | 33.5 s |
| Incremental test-profile rebuild | 8.8 s |
| `v2_step8_regression_smoke` (2 tests) | 14.1 s runtime, both pass |
| `v2_amg_scalar_parity` (4 tests) | 0.08 s runtime, all pass |

The 8.5a bit-parity smoke test (`disabled_runs_are_bit_deterministic`)
passes with LTO enabled — confirms LTO does not perturb the
Jacobi-CG default path's byte-identical guarantee on this hardware /
thread-count combination.

Cold-build penalty is within the 20–60 s budget of D4.

## Phase 1 — parallel_reduce.rs helpers (TBD)

Deterministic chunk-sequential reduction helpers. D1 reference pattern
(chunk → sequential per chunk → sequential sum of chunks in index
order) applied to dot product, axpy, norm, and max-abs.

## Phase 2 — Jacobi path parallelised (TBD)

Matrix-free `apply_momentum`, Jacobi preconditioner apply, CG inner
products.

## Phase 3 — RBGS replaces SGS (TBD)

Red-Black coloring at level 0 (structured grid), algebraic greedy
coloring at coarser levels (C/F splitting is algebraic).

## Phase 4 — AMG setup parallelised (TBD)

Galerkin R·A·P, prolongation/restriction row construction, sparse
assembly.

## Phase 5 — Newton extrapolation order 2 (TBD)

`v_guess(t) = 2·v(t-Δt) − v(t-2Δt)`, with step-0/step-1 initialisation
fallbacks and residual-regression safeguard.

## Phase 6 — Benchmarks (TBD)

Multi-run measurement (5 runs, mean ± std per D6). Four configurations
per case: Jacobi-8.5a-baseline / Jacobi-serial-LTO / Jacobi-parallel /
AMG-parallel.

## Phase 7 — Physics re-runs + downstream default (TBD)

Step 0–8 physics runs (scalar-parity + wallclock targets), README
default recommendation updated.

## Gate summary (to be filled)

| Gate | Target | Measured | Status |
|---|---|---|---|
| Jacobi step0-7 wallclock gain | ≥ ×4 | — | — |
| Jacobi step8 wallclock gain | ≥ ×3 | — | — |
| AMG step0-7 wallclock gain | ≥ ×2 | — | — |
| AMG/Jacobi ratio on step0-7 | ≤ 1.0 | — | — |
| Scalar-parity step0-7 physics (Jacobi) | < 1e-10 rel | — | — |
| Scalar-parity step8 physics (Jacobi) | < 1e-10 rel | — | — |
| Determinism across thread counts (fixed count) | byte-identical | — | — |
| RBGS vs SGS convergence | within 5 % | — | — |
| Newton convergence preserved | ≥ 95 % | — | — |
