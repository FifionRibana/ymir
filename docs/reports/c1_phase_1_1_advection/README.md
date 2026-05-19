# C1 Phase 1.1 — Advection-only prototype outputs

Issue #120 visual + scalar outputs. Regenerate with:

```
cargo test -p ymir-core --test c1_phase_1_1_advection -- --nocapture
```

## Run configuration

- Grid: 64²
- Seed: 42
- Steps: 300 (forward Euler)
- Kinematics: hand-tuned 8-plate preset
  ([`crates/ymir-core/src/tectonics_c1/kinematics.rs::PlateKinematics::preset_phase_1_1`](../../../crates/ymir-core/src/tectonics_c1/kinematics.rs))
  — two convergent cardinal pairs (E/W, N/S), one divergent
  diagonal pair (NE/SW), two diagonal rounders (SE, NW).
- CFL: `Δt = 0.5 · dx / max|v|` with `max|v| = √2 · 0.008 ≈
  0.0113`, giving `Δt ≈ 0.69` non-dim per step.
- Init: v2 `generate_voronoi` (default 8 plates, 30 % continental
  ratio) + v2 `init_s_field(InitMode::default())`
  (`Uniform { boundary_smoothing_width = 1.0 }`). Continental
  cells `S̃ = 1.0`, oceanic `S̃ = 0.2`, smoothstep blend on the
  1-cell collar.

## Snapshots

Each cycle dumps two PNGs:

| File | What it shows | Palette |
|---|---|---|
| `cycle_NNN_altitude.png` | Airy-isostasy heightmap of `S̃` (Phase B viewport) | **per-frame** auto-rescale through `tectonics::isostasy::compute_isostasy` |
| `cycle_NNN_s.png` | Direct `S̃` heightmap | **absolute** palette `[0, 2.0]`, saturating above |

The fixed-palette `cycle_NNN_s.png` series is the transport-
correctness check; the auto-rescale `altitude` PNGs are useful
to anticipate how the downstream isostasy phase will frame the
output as the C1 pipeline grows in Phase 1.4.

Snapshots saved: cycle 0, 50, 100, 200, 300.

## Visual reading guide

- **cycle 000 — initial state**: a single coherent Voronoï
  landmass (brown on dark blue ocean) with the 1-cell
  smoothstep transitions at inter-plate boundaries.
- **cycle 050 — visible transport signal**: continental
  material has been dragged by the per-plate velocities in
  different directions, breaking the coherent landmass; bright
  saturated points appear at convergence corners where mass
  piles up because no closure absorbs it.
- **cycles 100 / 200 / 300 — pile-up dominance**: most cells
  drained to ~0 (divergence zones have grown to dominate grid
  coverage); thin lines of accumulated mass remain at persistent
  convergence boundaries. `S̃` max climbs to ~1080 × initial,
  intentionally — Phase 1.1 has no closure to bound this.

## Acceptance signal

Per Issue #120:

| Criterion | Result |
|---|---|
| Mass conservation `< 1e-6` drift over 300 steps | ✓ drift `1.6e-14` (machine precision) |
| No NaN anywhere | ✓ asserted per-step |
| Wall time `<< 1 minute` at 64² | ✓ `245 ms` total, `820 µs/step` |
| Convergent boundary thickening visible | ✓ bright points at convergence corners in `cycle_050_s.png` |
| Divergent boundary thinning visible | ✓ ocean-blue spreading across most cells by `cycle_100_s.png` |
| Cratonic cells transport rigidly | Partial — covered by the mass-conservation + no-NaN invariants; full spatial check is deferred to a later phase that adds a cratonic-mask overlay |

## Per-cycle `S̃` distribution

From `cargo test --nocapture` output:

```
cycle 000: S̃ in [0.20,    1.00], mean = 0.557, std = 3.78e-1
cycle 050: S̃ in [0.00,  234.0],  mean = 0.557, std = 5.06e0
cycle 100: S̃ in [0.00,  627.2],  mean = 0.557, std = 1.23e1
cycle 200: S̃ in [0.00, 1077.7],  mean = 0.557, std = 1.89e1
cycle 300: S̃ in [0.00, 1079.7],  mean = 0.557, std = 1.89e1
```

The constant mean confirms exact mass conservation. The growing
max captures the "no closure to absorb pile-up" Phase 1.1
characteristic. The std stabilising between cycle 200 and 300
hints that first-order upwind's implicit numerical dissipation
eventually balances the convergent build-up.

## What this output is **not**

A plausible continent. The PNGs at cycles 100-300 are not meant
to look like a planet — they're meant to confirm that mass
moves with plate velocity, accumulates at convergence
boundaries, and drains at divergence boundaries. Mountain
morphology, isostatic balance, erosion, climate, biomes — none
of those exist in this code yet.

Phase 1.2 adds the Davis-Suppe orogenic-profile closure that
absorbs the convergence pile-up into mountain elevation,
bounding the time-asymptotic `S̃` max well below the Phase 1.1
unbounded value. Phase 1.3 adds equilibrium height. Phase 1.4
adds erosion + isostasy + the end-to-end heightmap pipeline.
See [`docs/design/c1_lightweight_dynamic_tectonics.md`](../../design/c1_lightweight_dynamic_tectonics.md)
§7 for the full plan.
