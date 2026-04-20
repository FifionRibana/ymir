# Step 0 — Nondim Stokes core + S advection (baseline)

> **Baseline reference for milestone "Solver reconstruction".**
> Do NOT compare against earlier reports — this is the first report.
> Subsequent steps' reports will diff against this one.

- Seed: `42`
- Entry-condition decisions archived in `tectonics_v2/README.md`.

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
| discretization | MAC staggered (v face / P η S cell-centre) |
| harmonic averaging | harmonic 4-point for η at corners |
| preconditioner | block-diag Jacobi (v) + diag(1/η) mass (P), null-space wrapped |
| gauge fixing | mean(P), mean(vx), mean(vy) projected before & after every M^-1 and once post-solve |
| outer CG tolerance | 1.0e-8 |
| inner CG tolerance | 1.0e-10 |
| outer CG max iter | 200 |
| inner CG max iter | 500 |
| CFL factor | 0.30 |
| grid spacing (nondim) | 0.015625 |
| body force | SinusoidalForce(ε=0.1) |
| seed | 42 |

### Timing

- wallclock total: `0.382 s`
- wallclock per step (mean): `1.273 ms`
- steps: `300`

### Solver health

- κ(A) estimate: N/A — outer CG converged in 0 iterations (the Kolmogorov-like placeholder forcing produces an exactly divergence-free velocity from A⁻¹f, so the Schur complement problem is trivially satisfied by p=0). The framework slot is exercised; real κ estimates come online at Step 2 when GPE spreading makes the Schur-complement nontrivial.
- effective η_max/η_min over run: `1.000` (placeholder; trivially 1.0 at Step 0)
- outer CG iterations — mean: `0.0`, max: `0`
- outer CG iteration histogram (5 bins):

  | bin ≤ | count |
  |---|---|
  | 0 | 300 |
  | 0 | 0 |
  | 0 | 0 |
  | 0 | 0 |
  | 0 | 0 |

- inner CG iterations (per inner solve) — mean: `1.0`, max: `1`

### Mass conservation of S

- initial mass: `4.096000000e3`
- final mass: `4.096000000e3`
- relative drift: `3.331e-16`

### Null-space health (post-solve means)

- max |mean(P)| across solves: `0.000e0`
- max |mean(vx)|: `3.895e-20`
- max |mean(vy)|: `0.000e0`

### Velocity magnitude

- peak |v|: `2.532e-3`

### Heightmaps of S

- `docs/reports/step0_heightmaps/s_64x64_t0000.png`
- `docs/reports/step0_heightmaps/s_64x64_t0150.png`
- `docs/reports/step0_heightmaps/s_64x64_t0300.png`

### Dormant metrics (inactive at Step 0)

| metric | activated at |
|---|---|
| S̃_eq (active-orogen mean thickness) | Step 5+ |
| boundary type diversity | Step 5 |
| yielding cell fraction | Step 3 |
| cratonic stability | Step 9 |
| Newton outcome distribution | Step 1 |
| age field stats | Step 10 |

## Grid 128×128

### Solver configuration

| field | value |
|---|---|
| discretization | MAC staggered (v face / P η S cell-centre) |
| harmonic averaging | harmonic 4-point for η at corners |
| preconditioner | block-diag Jacobi (v) + diag(1/η) mass (P), null-space wrapped |
| gauge fixing | mean(P), mean(vx), mean(vy) projected before & after every M^-1 and once post-solve |
| outer CG tolerance | 1.0e-8 |
| inner CG tolerance | 1.0e-10 |
| outer CG max iter | 200 |
| inner CG max iter | 500 |
| CFL factor | 0.30 |
| grid spacing (nondim) | 0.007812 |
| body force | SinusoidalForce(ε=0.1) |
| seed | 42 |

### Timing

- wallclock total: `1.376 s`
- wallclock per step (mean): `4.586 ms`
- steps: `300`

### Solver health

- κ(A) estimate: N/A — outer CG converged in 0 iterations (the Kolmogorov-like placeholder forcing produces an exactly divergence-free velocity from A⁻¹f, so the Schur complement problem is trivially satisfied by p=0). The framework slot is exercised; real κ estimates come online at Step 2 when GPE spreading makes the Schur-complement nontrivial.
- effective η_max/η_min over run: `1.000` (placeholder; trivially 1.0 at Step 0)
- outer CG iterations — mean: `0.0`, max: `0`
- outer CG iteration histogram (5 bins):

  | bin ≤ | count |
  |---|---|
  | 0 | 300 |
  | 0 | 0 |
  | 0 | 0 |
  | 0 | 0 |
  | 0 | 0 |

- inner CG iterations (per inner solve) — mean: `1.0`, max: `1`

### Mass conservation of S

- initial mass: `1.638400000e4`
- final mass: `1.638400000e4`
- relative drift: `0.000e0`

### Null-space health (post-solve means)

- max |mean(P)| across solves: `0.000e0`
- max |mean(vx)|: `2.626e-20`
- max |mean(vy)|: `0.000e0`

### Velocity magnitude

- peak |v|: `2.533e-3`

### Heightmaps of S

- `docs/reports/step0_heightmaps/s_128x128_t0000.png`
- `docs/reports/step0_heightmaps/s_128x128_t0150.png`
- `docs/reports/step0_heightmaps/s_128x128_t0300.png`

### Dormant metrics (inactive at Step 0)

| metric | activated at |
|---|---|
| S̃_eq (active-orogen mean thickness) | Step 5+ |
| boundary type diversity | Step 5 |
| yielding cell fraction | Step 3 |
| cratonic stability | Step 9 |
| Newton outcome distribution | Step 1 |
| age field stats | Step 10 |

---
*Generated by `cargo run --release --bin step_baseline`.*
