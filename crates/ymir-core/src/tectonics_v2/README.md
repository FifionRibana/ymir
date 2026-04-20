# tectonics_v2 — Step 0

Incremental rebuild of the tectonic solver. Step 0 ships a linear,
constant-viscosity Stokes solver coupled to a passive advected
thickness field `S̃` on a fully periodic toroidal domain, plus the
diagnostics framework the rest of the milestone extends.

Physical reference: [`docs/solver-scaling.md`](../../../../docs/solver-scaling.md).
Milestone tracker: [`docs/solver-reconstruction-roadmap.md`](../../../../docs/solver-reconstruction-roadmap.md).

## Entry-condition decisions (archived for the milestone)

1. **Viscosity scale `η*`.** The design note's formula
   `η* = ρ*·g·τ*·S*` evaluated at default primary scales
   (`ρ* = 3300 kg/m³`, `g = 9.81 m/s²`, `τ* = 30 Myr`, `S* = 35 km`)
   gives
   `η* = 3300 × 9.81 × (30·3.156·10¹³) × 3.5·10⁴ ≈ 1.07·10²⁴ Pa·s`.
   The previous handoff's `~10²³` figure was an arithmetic slip.
   Adopted value: **`η* = 1.073·10²⁴ Pa·s`** (rounded display `~10²⁴`).
   Used everywhere in `scales.rs`; reported at solver startup.

2. **`Field2D` / `PeriodicIndex` audit.** The legacy types in
   `crates/ymir-core/src/tectonics/solver/field.rs` are self-contained
   (zero imports into `tectonics/`), covered by stride- and
   wrap-aware unit tests (square, rectangular, coprime), and carry no
   tectonic state. **Decision: re-export directly** via
   [`field.rs`](./field.rs) (`pub use` only). No duplication until
   the legacy module is retired at milestone end. No other symbol is
   imported from `tectonics/`.

3. **Gauge-fixing strategy.** The periodic torus admits a 1-D
   pressure null space (constant mode) and a 2-D velocity null space
   (constant `vx` and `vy`). The strategy is: **mean-subtract
   `P`, `vx`, `vy` both before and after every preconditioner
   application (`M⁻¹`), plus once more on the final iterate.**
   Implemented in [`stokes/nullspace.rs`](./stokes/nullspace.rs) and
   wrapped by [`stokes/precond.rs`](./stokes/precond.rs). Verified by
   `tests/v2_nullspace.rs` — a solve with RHS deliberately carrying
   nonzero mean in every component returns `|mean(P)|`,
   `|mean(vx)|`, `|mean(vy)| < 1·10⁻¹⁰`.

## Discretization choice

MAC (staggered) grid with periodic BCs, following the legacy layout:
- `p`, `η`, `S` at cell centres `((i+0.5)dx, (j+0.5)dy)`.
- `vx` at left vertical faces `(i dx, (j+0.5)dy)`.
- `vy` at bottom horizontal faces `((i+0.5)dx, j dy)`.
- `ε̇_xy` and `σ_xy` at nodal corners `(i dx, j dy)`, with η there
  computed by **harmonic 4-point averaging** of the surrounding cell
  centres. Harmonic averaging reduces trivially to `η` for constant
  fields, so the assembly is identical at Step 0 and Step 1 — the
  variable-η rheology in Step 1 is a data change, not an assembly
  change.

The MMS convergence test verifies this discretization is 2nd-order in
`v`:
```
N=16:  v_err=6.475e-3
N=32:  v_err=1.609e-3  (slope 2.008)
N=64:  v_err=4.018e-4  (slope 2.002)
N=128: v_err=1.004e-4  (slope 2.001)
```

## Linear solver

Pressure Schur-complement with nested **conjugate gradient**
(`ConjugateGradient` implementing the [`LinearSolver`][trait] trait).
The outer CG applies `S = B·A⁻¹·B^T` to pressure iterates; each
application invokes an inner CG on the momentum block `A`. Both blocks
are SPD on their respective zero-mean subspaces; CG is mathematically
appropriate throughout.

The `LinearSolver` trait is the integration point for BiCGSTAB at
Step 3 when plastic yielding makes the global system non-symmetric.
Direct calls to `ConjugateGradient::solve` outside this trait are
forbidden by convention.

[trait]: ./stokes/solver.rs

### Preconditioner

Block-diagonal:
- **Velocity block** `M_v = diag(A)⁻¹` (Jacobi) with a configurable
  floor on the diagonal (`cfg.diag_floor`, default `1e-20`).
- **Pressure block** `M_p = diag(1/η)` (viscosity-scaled mass matrix).

Both wrap the mean-projection in the `apply` entry points.

## Transverse T1 — absorbed

The roadmap previously listed "T1 — Null-space-aware preconditioner"
as an independent transverse task. Step 0 now ships the null-space-
aware preconditioner as part of the core solver, so T1 is absorbed.
[`docs/solver-reconstruction-roadmap.md`](../../../../docs/solver-reconstruction-roadmap.md)
is updated to remove T1 from the open transverse list.

## Dormant metrics

The [`Metrics`](./diagnostics/metrics.rs) struct declares the full
list of metrics the milestone will eventually track. Those not
applicable at Step 0 are `Option<_>` and stay `None`:

| metric | activated at |
|---|---|
| `s_eq` (active-orogen mean thickness) | Step 5+ |
| `boundary_type_diversity` | Step 5 |
| `yielding_cell_fraction` | Step 3 |
| `cratonic_stability` | Step 9 |
| `newton_outcome_distribution` | Step 1 |
| `age_field_stats` | Step 10 |

## Files

```
tectonics_v2/
├── README.md                 ← this file
├── mod.rs                    ← public re-exports
├── scales.rs                 ← Scales + dim ↔ nondim conversions
├── field.rs                  ← re-exports Field2D, PeriodicIndex
├── forcing.rs                ← BodyForce trait + ZeroForce + SinusoidalForce
├── advection.rs              ← conservative first-order upwind S advection
├── stokes/
│   ├── mod.rs                ← pressure Schur-complement coordinator
│   ├── nullspace.rs          ← mean projectors for P, vx, vy
│   ├── operator.rs           ← MAC momentum + divergence + adjoint
│   ├── precond.rs            ← block-diag preconditioners with null-space wrapping
│   └── solver.rs             ← LinearSolver trait + CG
└── diagnostics/
    ├── mod.rs                ← re-exports
    ├── metrics.rs            ← Metrics + histogram + SolverConfigDump
    ├── report.rs             ← markdown writer
    └── harness.rs            ← baseline runner (seed-parametrisable)
```

Plus the binary entry point at
[`crates/ymir-core/src/bin/step_baseline.rs`](../../bin/step_baseline.rs).

## Running the baseline

```bash
cargo run --release --bin step_baseline -- \
    --seed 42 --grids 64,128 --steps 300 \
    --output docs/reports/step0_report.md
```

Writes the markdown report and a heightmap set under
`docs/reports/step0_heightmaps/`.

## Note on the placeholder body force

The Step 0 spec names `f̃ = ε · sin(2π x̃ / L̃x) ê_x`. That force is a
pure gradient, so under incompressible Stokes it is balanced entirely
by pressure and produces `ṽ ≡ 0`. With `ṽ = 0`, the advection step is
a no-op and the coupled smoke test is trivial. This crate ships the
**rotational variant** `f̃ = ε · sin(2π ỹ / L̃y) · ê_x` instead,
which drives a simple Kolmogorov-like shear flow (analytic steady
solution `ṽx = ε · sin(2π ỹ) / (2π/L̃y)²`), still periodic, with
non-trivial motion and non-zero `peak |ṽ|`. The deviation is
documented in [`forcing.rs`](./forcing.rs); the spec's intent
("~3% S displacement over 300 steps, measurable but small") is
preserved.
