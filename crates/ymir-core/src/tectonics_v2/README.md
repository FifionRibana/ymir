# tectonics_v2 — Step 0

Incremental rebuild of the tectonic solver. Step 0 ships a linear,
constant-viscosity **thin viscous sheet** solver (England & McKenzie
1982) coupled to a conservative upwind advection of the crustal
thickness `S̃` on a fully periodic toroidal domain, plus the
diagnostics framework the rest of the milestone extends.

Physical reference: [`docs/solver-scaling.md`](../../../../docs/solver-scaling.md).
Milestone tracker: [`docs/solver-reconstruction-roadmap.md`](../../../../docs/solver-reconstruction-roadmap.md).

## Formulation: thin viscous sheet, NOT incompressible Stokes

The equations are:

```
Momentum:   -∇·(2 η̃ ε̇̃(ṽ)) = Ar·∇Φ̃ + f̃_ext       (1) [elliptic, SPD]
Thickness:  ∂_t̃ S̃ + ∇·(S̃ ṽ) = Q̃                  (2) [mass balance]
```

Two properties that distinguish this from incompressible Stokes and
are foundational for every Step 1–10 that follows:

- **No incompressibility constraint.** `∇·ṽ ≠ 0` in 2-D is physically
  meaningful — it is the rate at which the crustal column thickens
  (plug it into (2) and it becomes the source of `∂_t̃ S̃`). Enforcing
  `∇·ṽ = 0` would make it impossible for orogens to thicken under
  convergence or for rifts to thin under divergence.
- **No pressure unknown.** `Φ̃` is the gravitational potential
  energy and enters (1) as a physical driving term weighted by the
  Argand number `Ar` — not as a Lagrange multiplier dual to a
  constraint. Nothing in the solver carries a pressure field.

Consequence for the linear algebra: (1) assembles an SPD operator on
the velocity. A **single preconditioned conjugate-gradient solve per
time step** suffices — no saddle point, no Schur complement, no
nested iteration.

### Early faux-départ note

The first Step 0 commit (`a8c8f3a`) shipped an incompressible-Stokes
saddle-point solver (pressure Schur complement with nested CG). That
was an implementation error, not a design choice: the spec never
called for incompressibility. It passed the original unit tests
because the placeholder forcing happened to produce a div-free
velocity (so the constraint did not bite) and the original MMS
manufactured solution was div-free by construction. The current tree
replaces that solver with the correct thin-sheet elliptic operator;
the faux-départ commit remains in the history as documentation of
the correction trajectory.

## Entry-condition decisions (archived for the milestone)

1. **Viscosity scale `η*`.**
   `η* = ρ*·g·τ*·S* = 3300 × 9.81 × (30·3.156·10¹³) × 3.5·10⁴
   ≈ 1.073·10²⁴ Pa·s`. The previous handoff's `~10²³` was an
   arithmetic slip. Adopted value: **`η* = 1.073·10²⁴ Pa·s`** (rounded
   display `~10²⁴`). Reported at solver startup.

2. **`Field2D` / `PeriodicIndex` audit.** The legacy types in
   `crates/ymir-core/src/tectonics/solver/field.rs` are self-contained
   (zero external imports), tested for stride and wrap (square,
   rectangular, coprime), and carry no tectonic state. **Decision:
   re-export directly** via [`field.rs`](./field.rs) (`pub use` only).
   No other symbol is imported from `tectonics/`.

3. **Gauge-fixing strategy.** The fully periodic torus admits a **2-D
   velocity null space** (constant `vx` and `vy`, rigid-body
   translation). The strategy is: **mean-subtract `vx` and `vy`
   before and after every preconditioner application `M⁻¹`, plus once
   more on the final iterate.** Verified by `tests/v2_nullspace.rs`
   with RHS deliberately carrying a nonzero mean in each component.

## Discretization choice

MAC (staggered) grid with periodic BCs:
- `η`, `S` at cell centres `((i+0.5)dx, (j+0.5)dy)`.
- `vx` at left vertical faces `(i dx, (j+0.5)dy)`.
- `vy` at bottom horizontal faces `((i+0.5)dx, j dy)`.
- `ε̇_xy` at nodal corners `(i dx, j dy)`, with η there by
  **harmonic 4-point averaging** of surrounding cells. Harmonic
  averaging reduces to `η` for constant fields, so the assembly is
  identical at Step 0 and Step 1 — the power-law rheology in Step 1
  is a data change on `η`, not an assembly change.

For constant η the discrete operator reduces to
`A v = -η (∇² v + ∇(∇·v))`. The grad-div part is essential: it
couples `vx` and `vy` through normal strain, and dropping it would
reduce the 2-D thin-sheet to two decoupled scalar Laplacians — wrong
physics. The MMS test uses a deliberately non-div-free manufactured
solution (`v = (sin(2πx), sin(2πy))`) so that both terms are
exercised.

MMS convergence (rel tol 1e-12, manufactured non-div-free solution):

```
N=16:  v_err=9.16e-3
N=32:  v_err=2.28e-3  (slope 2.008)
N=64:  v_err=5.68e-4  (slope 2.002)
N=128: v_err=1.42e-4  (slope 2.001)
```

## Linear solver

Preconditioned **conjugate gradient** (`ConjugateGradient`
implementing the [`LinearSolver`][trait] trait). The trait is the
integration point for BiCGSTAB at Step 3 when plastic yielding makes
the system non-symmetric. Direct calls to `ConjugateGradient::solve`
outside the trait are forbidden by convention.

Preconditioner: `M = diag(A)⁻¹` (Jacobi), wrapped with the velocity
mean-projection (see [`precond.rs`](./stokes/precond.rs)). A
configurable floor on `|diag|` protects against degenerate-η cells.

[trait]: ./stokes/solver.rs

## Transverse T1 — absorbed

The roadmap previously listed "T1 — Null-space-aware preconditioner"
as an independent transverse task. Step 0 ships the null-space-aware
preconditioner as part of the core solver, so T1 is absorbed.
[`docs/solver-reconstruction-roadmap.md`](../../../../docs/solver-reconstruction-roadmap.md)
reflects this.

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
│   ├── mod.rs                ← thin-sheet solver entry point
│   ├── nullspace.rs          ← mean projectors for vx, vy
│   ├── operator.rs           ← MAC momentum operator
│   ├── precond.rs            ← velocity Jacobi with null-space wrapping
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

`f̃ = ε · sin(2π x̃ / L̃x) · ê_x` per the Step 0 spec. In the
thin-sheet formulation this force **produces flow** (no pressure
available to absorb it as a gradient); the analytic steady solution
with constant η is `ṽx = ε · sin(2π x̃ / L̃x) / (8 π² η / L̃x²)`,
`ṽy = 0`, giving `peak|ṽ| ≈ 1.27·10⁻³` at `ε = 0.1`, `L̃x = 1`,
`η = 1`. The single-Fourier-mode character of this solution is why
CG converges in very few iterations at Step 0; the iteration count
becomes meaningful once Step 1 (power-law η) and Step 2 (GPE)
introduce real structure in the solution.
