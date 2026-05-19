# Step 8.6 Phase 8g — visual revalidation post-corrections

Re-run of the Phase 7 4-preset gallery
(`active_medley`, `convergence`, `divergence`, `subduction`) with every
Step 8.6 augmenté correction live:

- Phase 8a `InitMode::Uniform` (default) — flat per-plate-type
  initialisation, no more sinusoidal perturbation in S̃.
- Phase 8b Voronoï plate-boundary overlay (black) and per-plate
  velocity-vector overlay (yellow arrows) burned into every
  captured frame.
- Phase 8c real-time metrics dashboard (UI surface, not in PNGs).
- Phase 8d Init mode + missing physics knobs exposed (no preset
  changes — all four still load with their pre-Phase-8d JSON via
  the `#[serde(default)]` shims).
- Phase 8e export/import system (UI surface, not in PNGs).

## Run parameters

| Setting           | Value |
|-------------------|-------|
| Grid              | 32²   |
| Steps             | 100   |
| Capture interval  | every 10 steps (+ step 0) |
| Init mode         | `Uniform { boundary_smoothing_width: 1.0 }` |
| Overlays          | Voronoï boundaries + velocity vectors |
| Wallclock total   | 2515 s ≈ 42 min (4 presets, sequential) |
| Per-preset wallclock | active_medley 518 s · convergence 726 s · divergence 753 s · subduction 518 s |

A 64² verification run is optional follow-up; the Phase 8f
equilibrium analysis indicates the physics-relevant metrics
(`peak|v|`, yielding pattern, cratonic structure) reach a steady
state by step 100 at this scale, so the rendered patterns won't
change meaningfully.

## Reviewer judgment (per acceptance criterion)

| Criterion | active_medley | convergence | divergence | subduction |
|-----------|---------------|-------------|------------|------------|
| No residual sinusoidal / checkerboard pattern in S̃ | yes | yes | yes | yes |
| Continents emerge from dynamics, not from init | yes | yes | yes | yes |
| Voronoï overlay clarifies plate structure | yes | yes | yes | yes |
| Velocity arrows show plausible directions for the regime | yes | yes | yes | yes |
| ε̇_II patterns coherent (cratons dark, mobile belts active) | yes | yes | yes | yes |
| age / cratonic / \|v\| comparable to Phase 7 baseline | yes | yes | yes | yes |

Per-preset notes:

- **active_medley**: S̃ rows show flat per-plate fills with the
  Voronoï mesh in black; continental cells (green) and oceanic cells
  (teal) are clearly delineated by the overlay rather than by the
  sinusoidal artefact that Phase 7 flagged. Strain rate shows the
  expected mobile-belt + craton pattern. Velocity arrows are visible
  in the cratonic/strain rows where the field is dark enough for
  yellow to read.
- **convergence**: continental masses are smaller than active_medley
  (the spec carries a different `continental_ratio`); yielding bands
  in ε̇_II concentrate at plate-collision interfaces, exactly the
  shape this preset is meant to demonstrate.
- **divergence**: per-plate flat S̃, Voronoï clean, yielding patterns
  trace the rifting boundaries. The "L-shape craton" the Phase 7
  follow-up flagged is preserved with the new init.
- **subduction**: ε̇_II shows pronounced yielding near subduction
  edges; arrows point along the expected subduction direction
  through the dark cratonic cells. No `vmax` runaway — the regime
  stays within the Step 8 envelope.

Each per-field patchwork (`_<field>_patchwork.png`) tiles the 11
captured steps in a roughly square grid sorted by step ascending,
with 2-px black gutters between cells. Each `_all.png` stacks the
five field patchworks vertically (S̃, age, cratonic, strain, vmag).

## Verdict

**All four presets pass the Phase 8g acceptance criteria.** The
Phase 8a `Uniform` init removed the sinusoidal artefact that
Phase 7 reviewer flagged; the Phase 8b overlays make the plate
structure and motion legible without crowding the field rendering;
ε̇_II behaves physically (cratons dark, mobile belts active); no
preset shows runaway or NaN/artefact symptoms.

**Phase 8h sunset of the legacy bridge is authorised.**

## Artefacts

Patchwork images (committed):
- `active_medley/_all.png` + 5 per-field patchworks
- `convergence/_all.png` + 5 per-field patchworks
- `divergence/_all.png` + 5 per-field patchworks
- `subduction/_all.png` + 5 per-field patchworks

Per-step source PNGs (`step_NNNN_<field>.png`) and `final/<field>.png`
are gitignored — regenerable via:

```text
cargo test --release -p ymir-viz --test v2_phase8g_visuals \
    --jobs 1 -- --ignored --nocapture
```

Override grid / step count / capture interval via env vars:
`YMIR_PHASE8G_GRID`, `YMIR_PHASE8G_STEPS`, `YMIR_PHASE8G_INTERVAL`.
