# Step 8.6 Phase 7 — Reviewer checkpoint

> **Hard checkpoint** per `step8_6_issue.md` D7 / §"Reviewer Phase 7 est
> un checkpoint dur". Phase 8 sunset is **not authorised** until you
> validate the visuals here.

## What's in this directory

For each of the six v2 presets registered in
`crates/ymir-viz/src/bridge/v2/presets.rs`, the harness ran a
`32² × 50` step baseline and dumped all five raster fields as PNG:

```
step8_6_phase7_gallery/
├── quiescent/
│   ├── s.png           ← S̃ thickness (hypsometric)
│   ├── age.png         ← geological age A (linear teal→orange)
│   ├── cratonic.png    ← cratonic factor (linear grayscale, 0-1)
│   ├── strain.png      ← ε̇_II (log scale, purple→yellow)
│   └── vmag.png        ← |v| velocity magnitude (log scale)
├── single_continent/   ← same 5 fields
├── convergence/
├── subduction/
├── divergence/
└── active_medley/
```

The grid was deliberately reduced to `32²` for the gallery so it
builds in ~5 minutes total. The full-resolution UI (run via
`YMIR_BRIDGE=v2 cargo run --release -p ymir-viz`) lets you re-run
each preset at `64²` × 100 steps for higher-detail spot checks.

## What to validate

Per D2 (visual coherence acceptance) on each preset:

1. **No NaN / Inf artifacts.** No black holes, no banding, no
   uniform-coloured patches that look like overflow / underflow.
2. **Continents identifiable.** `s.png` should show distinct
   bright (continental) and dark (oceanic) regions matching the
   Voronoï layout — not a uniform grey wash.
3. **Age field gradient (preset-dependent).**
   - Quiescent / single_continent: age field mostly grows uniformly
     (no boundary events) — large continental block stays bright,
     oceanic stays dark.
   - active_medley: heavy ridge / arc / collision activity →
     scattered teal "freshly reset" cells, ochre interior bands.
4. **Cratonic factor smooth.** White at the cratonic-core cells,
   grading through grey at the outer boundary, black on oceanic
   cells. Hard-edge stripes mean the BFS / smoothstep is broken.
5. **ε̇_II contrast.** When the preset enables yielding (every
   preset except `quiescent`), `strain.png` should show *localised*
   hot bands at boundary cells and quiet (purple) interiors. A
   uniformly hot map suggests every cell is yielding (regime
   pathology).
6. **|v| dynamic range.** `vmag.png` log scale should show a
   bright streamline pattern when mantle is enabled (every preset
   except `quiescent` / `single_continent`). Quiescent presets
   should be uniformly dark — that's the §7 expected behaviour
   (peak\|v\| ~ 1e-5).

## Per-preset expectation cheatsheet

| Preset | Expected character |
|---|---|
| `quiescent` | Tutorial. Static-looking S, dark vmag, low strain. Sanity. |
| `single_continent` | One large bright Voronoï cell on `s`, otherwise dark. Mantle off → vmag near zero. Cratonic-on shows white interior on `cratonic.png`. |
| `convergence` | Mantle pattern visible on vmag. Strain peaks at plate convergent margins (where boundary_flag = ContinentalCollision). |
| `subduction` | Mantle Mf=1.5 (stronger), more flow visible. Strain band at oceanic-continental interface; arcs / ridges on `age` should show as resets. |
| `divergence` | Same as convergence but pulling apart — ridges in age field, basins on S. |
| `active_medley` | Step 8 / Step 9 baseline. The "reference" — should match `docs/reports/step9_phase7_baseline/` qualitatively. |

## How to also run the live UI

```bash
# Default — bridge mode is v2.
cargo run --release -p ymir-viz

# Force legacy (acceptance #11 regression check):
YMIR_BRIDGE=legacy cargo run --release -p ymir-viz
```

In the v2 panel:

- **Preset dropdown** at the top — switch between the six presets.
- **Display** section — switch between the five field rasters.
- **▶ Run** kicks off the simulation. Wait ~10-30 s at 32², longer
  at 64².
- **📷 Capture** dumps the current field as PNG under
  `<output_dir>/screenshots/`. Filename includes preset + field +
  unix timestamp.

## Checkpoint outcome

When you've reviewed:

- ✅ **All six presets visually coherent** → proceed to Phase 8
  (sunset). Tell the agent "phase 7 ok, lance la phase 8" or similar.
- ❌ **One or more presets look broken** → tell the agent which
  preset + which field + what's wrong. The agent diagnoses (v2
  solver vs visualisation) before tuning anything (D8 anti-pattern).
