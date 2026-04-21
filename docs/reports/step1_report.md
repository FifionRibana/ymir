# Step 1 — Power-law rheology + Newton solver (baseline)

> **Step 1 reference run for milestone "Solver reconstruction".**
> Compared against Step 0 (`docs/reports/step0_report.md`).
> Subsequent steps' reports will diff against this one.

- Seed: `42`
- Formulation: thin viscous sheet + power-law rheology (n > 1). Linear internal system per Newton iteration is symmetric (variational structure); CG suffices.

## Physical scales

```
Scales:
- L* = 3.500e5 m (350 km)
- S* = 3.500e4 m (35 km)
- τ* = 9.467e14 s (30.00 Myr)
- ρ* = 3300.0 kg/m³
- v* = 3.697e-10 m/s (1.167 cm/yr)
- ε̇* = 1.056e-15 1/s
- η* = 1.073e24 Pa·s
- σ* = 1.133e9 Pa
- p* = 1.133e9 Pa
- f* = 3.237e3 N/m³
- Ar = 0.100
- 2π check (not used, informational): 6.283185
```

## Grid 64×64

### Solver configuration

| field | value |
|---|---|
| formulation | thin viscous sheet (elliptic, no pressure) with power-law rheology |
| discretization | MAC staggered (v face / η S cell-centre / ε̇_xy corner) |
| η averaging to corners | arithmetic 4-point at corners (see operator.rs) |
| preconditioner | velocity Jacobi (Picard-block diagonal), null-space wrapped |
| gauge fixing | mean(vx), mean(vy) projected before & after every M⁻¹ + post-solve |
| preset | `dynamic-accidented` |
| nonlinear solver | `newton` |
| rheology `n` (after continuation) | 3.00 |
| rheology `ε̇_min` | 1.0e-3 |
| rheology `η_max` (soft cap) | 1.0e3 |
| continuation schedule | `[1.0, 1.5, 2.0, 2.5, 3.0]` |
| Newton rel tol | 1.0e-6 |
| Newton max outer iters | 20 |
| CG tolerance | 1.0e-8 |
| CG max iter | 2000 |
| CFL factor | 0.30 |
| grid spacing (nondim) | 0.015625 |
| body force | SinusoidalForce(ε=0.1) |
| seed | 42 |

### Timing

- wallclock total: `0.408 s`
- wallclock per step (mean): `1.361 ms`
- steps: `300`

### Linear-solver health (CG inside Newton)

- κ(A) estimate from CG iterations (per Newton step): `2.46e0`
- CG iterations per Newton step — mean: `15.2`, max: `17`
- CG iteration histogram (5 bins):

  | bin ≤ | count |
  |---|---|
  | 4 | 1 |
  | 7 | 0 |
  | 10 | 0 |
  | 13 | 0 |
  | 17 | 22 |

### Newton (nonlinear) health

- outcome distribution — Converged: `100.0%`, Stalled: `0.0%`, Diverged: `0.0%`, CappedIters: `0.0%`
- Newton outer iters per timestep — mean: `0.1`, max: `6`
- effective η_max/η_min over run — mean: `1.04`, max: `1.04`
- cap-activation fraction (η_eff > 0.9·η_max) — during ramp: `0.000%`; steady state: `0.000%` (spec target < 1%)
- continuation ramp: ✅ all 5 sub-solves converged

### Mass conservation of S

- initial mass: `4.096000000e3`
- final mass: `4.096000000e3`
- relative drift: `-8.882e-16`

### Null-space health

- max |mean(vx)| across solves: `3.102e-25`
- max |mean(vy)|: `0.000e0`

### Velocity magnitude

- peak |v|: `1.306e-5`

### Heightmaps of S

- `docs/reports/step1_heightmaps/s_64x64_t0000.png`
- `docs/reports/step1_heightmaps/s_64x64_t0150.png`
- `docs/reports/step1_heightmaps/s_64x64_t0300.png`

### Comparison with Step 0

#### Grid 64×64 — comparison vs Step 0

| metric | previous | current | ratio / note |
|---|---|---|---|
| wallclock (s) | 0.152 | 0.408 | ×2.69 |
| CG iters / linear solve (mean) | 0.0 (solver-trivial) | 15.2 | N/A — no denominator; report absolute |
| S mass drift (relative) | 2.331e-15 | -8.882e-16 | gate 1e-10 |
| max \|mean(vx)\| | 2.647e-23 | 3.102e-25 | bruit machine |
| max \|mean(vy)\| | 0.000e0 | 0.000e0 | bruit machine |

### Dormant metrics (inactive at Step 1)

| metric | activated at |
|---|---|
| S̃_eq (active-orogen mean thickness) | Step 5+ |
| boundary type diversity | Step 5 |
| yielding cell fraction | Step 3 |
| cratonic stability | Step 9 |
| age field stats | Step 10 |

## Grid 128×128

### Solver configuration

| field | value |
|---|---|
| formulation | thin viscous sheet (elliptic, no pressure) with power-law rheology |
| discretization | MAC staggered (v face / η S cell-centre / ε̇_xy corner) |
| η averaging to corners | arithmetic 4-point at corners (see operator.rs) |
| preconditioner | velocity Jacobi (Picard-block diagonal), null-space wrapped |
| gauge fixing | mean(vx), mean(vy) projected before & after every M⁻¹ + post-solve |
| preset | `dynamic-accidented` |
| nonlinear solver | `newton` |
| rheology `n` (after continuation) | 3.00 |
| rheology `ε̇_min` | 1.0e-3 |
| rheology `η_max` (soft cap) | 1.0e3 |
| continuation schedule | `[1.0, 1.5, 2.0, 2.5, 3.0]` |
| Newton rel tol | 1.0e-6 |
| Newton max outer iters | 20 |
| CG tolerance | 1.0e-8 |
| CG max iter | 2000 |
| CFL factor | 0.30 |
| grid spacing (nondim) | 0.007812 |
| body force | SinusoidalForce(ε=0.1) |
| seed | 42 |

### Timing

- wallclock total: `1.655 s`
- wallclock per step (mean): `5.518 ms`
- steps: `300`

### Linear-solver health (CG inside Newton)

- κ(A) estimate from CG iterations (per Newton step): `1.05e1`
- CG iterations per Newton step — mean: `30.6`, max: `34`
- CG iteration histogram (5 bins):

  | bin ≤ | count |
  |---|---|
  | 7 | 1 |
  | 14 | 0 |
  | 20 | 0 |
  | 27 | 0 |
  | 34 | 22 |

### Newton (nonlinear) health

- outcome distribution — Converged: `100.0%`, Stalled: `0.0%`, Diverged: `0.0%`, CappedIters: `0.0%`
- Newton outer iters per timestep — mean: `0.1`, max: `6`
- effective η_max/η_min over run — mean: `1.04`, max: `1.04`
- cap-activation fraction (η_eff > 0.9·η_max) — during ramp: `0.000%`; steady state: `0.000%` (spec target < 1%)
- continuation ramp: ✅ all 5 sub-solves converged

### Mass conservation of S

- initial mass: `1.638400000e4`
- final mass: `1.638400000e4`
- relative drift: `-1.221e-15`

### Null-space health

- max |mean(vx)| across solves: `8.484e-22`
- max |mean(vy)|: `0.000e0`

### Velocity magnitude

- peak |v|: `1.305e-5`

### Heightmaps of S

- `docs/reports/step1_heightmaps/s_128x128_t0000.png`
- `docs/reports/step1_heightmaps/s_128x128_t0150.png`
- `docs/reports/step1_heightmaps/s_128x128_t0300.png`

### Comparison with Step 0

#### Grid 128×128 — comparison vs Step 0

| metric | previous | current | ratio / note |
|---|---|---|---|
| wallclock (s) | 0.527 | 1.655 | ×3.14 |
| CG iters / linear solve (mean) | 0.0 (solver-trivial) | 30.6 | N/A — no denominator; report absolute |
| S mass drift (relative) | -3.553e-15 | -1.221e-15 | gate 1e-10 |
| max \|mean(vx)\| | 1.423e-22 | 8.484e-22 | bruit machine |
| max \|mean(vy)\| | 0.000e0 | 0.000e0 | bruit machine |

### Dormant metrics (inactive at Step 1)

| metric | activated at |
|---|---|
| S̃_eq (active-orogen mean thickness) | Step 5+ |
| boundary type diversity | Step 5 |
| yielding cell fraction | Step 3 |
| cratonic stability | Step 9 |
| age field stats | Step 10 |

---
*Generated by `cargo run --release --bin step_baseline`.*
