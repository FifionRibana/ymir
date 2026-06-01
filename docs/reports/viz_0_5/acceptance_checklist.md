# Viz-0.5 — Issue #139 Acceptance Checklist

Hover-to-inspect + workflow-mode pipeline + hypsometric lens + continuation-ready structure. Manual UI/visual checklist (the headless part is automated; see bottom). Viz-only — `tectonics_c1` untouched (9th bit-identical guard at Final).

## How to run

```bash
cargo run --release -p ymir-viz --features v2_legacy
```

C1 engine selected by default; "C1 Engine Controls" panel top-left; "C1 Cell Inspector" panel bottom-right.

## Checklist

### 1. Hover-to-inspect (all views)

- [ ] Hover the map in **each** field view (S̃, Age, PlateId, PlateType, Altitude, Cratonic). The bottom-right "C1 Cell Inspector" updates with `cell (i,j)`, `S̃`, `age`, `altitude (non-dim)`, `altitude (m)`.
- [ ] **Non-dim altitude is shown FIRST** (the verification value); meters SECOND, labelled "after hypsometric curve".
- [ ] Cross-check a known cell: hover the **centre** of the map → `(i,j)` ≈ `(nx/2, ny/2)`. Hover an **off-centre** cell (e.g. upper-left) → small `i`, small `j` (Y-axis is top-down).
- [ ] Hover works in views OTHER than Altitude (the altitude cache derives independent of the active field).
- [ ] Moving off the map / over a panel clears the readout ("Hover over the map…").

### 2. Workflow mode — continent + migrating coast

- [ ] Switch **Pipeline → Workflow**. Press ▶ Run (default 64² seed 42, 5×20).
- [ ] The continent does NOT vanish — a stable landmass remains (Altitude / PlateType views). See the finding below: emergent land at seed 42 is ~5–6%, the honest above-sea fraction (NOT collapsed).
- [ ] The coast **migrates in discrete jumps** at each cycle boundary (every k_cycle steps) as `apply_post_tectonic` reclassifies — visible in PlateType / Altitude.

### 3. Gallery ↔ Workflow toggle (same seed)

- [ ] Run **Gallery** then **Workflow** at the same seed (64² seed 42). In Gallery the PlateType coast is static (Issue #137 contract); in Workflow it migrates per cycle. The difference should be visible.

### 4. Cadence sliders + calibration warning

- [ ] In Workflow mode, open "Workflow cadence". `n_cycles` (default 5) and `k_cycle` (default 20) sliders show, with the live "total = N×K" readout.
- [ ] Move a slider off the default → an **amber calibration warning** appears ("alpha=0.01 calibrated for k_cycle=20, n_cycles=5; raising without lowering alpha over-erodes"). Running off-default visibly over-erodes (expected — the A1-c lesson).

### 5. Hypsometric lens (meters only)

- [ ] Open "Hypsometric (hover meters)". Adjust **land scale** / **land gamma** → the hover **meters** value changes; the **non-dim** value does NOT (W3: cosmetic lens, never the verification value or the model).
- [ ] Ocean scale/gamma are shown read-only ("anchored Stein-Stein").

### 6. Gallery + v2 unchanged

- [ ] Gallery mode behaves exactly as Issue #137 (static coast, per-step animation, 6 fields, overlays).
- [ ] Switch Engine → v2 (legacy): v2 works unchanged; C1 controls gray out.

## Headless acceptance (automated)

```bash
cargo test --release -p ymir-viz --features v2_legacy --bin ymir-viz \
  bridge::c1::thread -- --nocapture
```

- `workflow_mode_produces_calibrated_cycle_count` — cycle-0 + n_cycles×(k_cycle per-step + 1 post-cycle) snapshots, monotone indices, total = n_cycles×k_cycle.
- `workflow_mode_coast_reclassifies` — `apply_post_tectonic` reclassify changes `plate_type` at a cycle boundary.
- `worker_retains_state_for_continuation` — `run_workflow_cycles` resumes a carried state (Stage E4 continuation core).
- `workflow_mode_continent_preserved` — **the A1-c regression guard** (qualitative signature): mass-conserving (`|Δmass|/mass < 1e-3`/cycle), converged (`|Δfrac| < 0.01` last cycle), and workflow `iso_land ≥ 0.9 × gallery iso_land` at the same 100 steps.
- `workflow_acceptance_continent_preserved_seed_42` — **product acceptance**: coast migrated substantially vs init (> 500 cells reclassified; measured 1205) AND emergent land nonzero (> 0.02; measured 0.0588).
- `workflow_continent_diagnostic` (`#[ignore]`) — per-cycle evidence trajectory.

## Stage V diagnostic findings (the ~0.058 continental fraction)

The calibrated workflow converges to a continental/emergent fraction of **~0.058** at seed 42 (NOT the ~0.20–0.45 originally estimated). The Stage V diagnostic established this is **correct**, not a regression:

1. **The ~0.27 init "continent" is a geometric LABEL, not emergent land.** Init `plate_type` is assigned by Voronoï continental seeds (0.27 of cells). The *above-sea-level* land (`compute_isostasy(s).land_ratio`, S̃-only) is ~0.05 even in the pure-tectonic gallery (0.0449 at 100 steps). The [0.20,0.45] target compared the geometric label to a final emergent fraction — apples/oranges, confirmed on evidence.
2. **Not macro erosion.** `macro_redistribution` is mass-conserving (`Δmass ≈ 0` every cycle, `isostatic_rebound_ratio = 0.80`). The continental loss is **declassification**, not mass transport.
3. **Adaptive `sea_level_ref` drift is the mechanism.** `sea_level_ref = s_min + 0.4·(s_max − s_min)` jumps 0.52 → 0.98 in cycle 1 as `s_max` doubles (1.0 → 2.18) from Davis-Suppe orographic peaks — declassifying 753 continental cells by threshold rise. The workflow then converges (last-cycle Δ 0.0007) and sits ABOVE the gallery isostatic floor (0.0588 > 0.0449).
4. **Pre-existing C1 debt, revealed not created.** The `min/max`-based `sea_level_ref` is fragile to Davis-Suppe peak outliers. A robust **percentile** sea level is a Phase 1.5 follow-up (core change, out of viz-only scope). This also explains "little emergent land at 256²" from the Viz-0 thread: C1 intrinsically produces little emergent land at seed 42 under the current sea-level formula — model behavior, not a viz bug. **Relevant for Phase 3 validation**: "little emergent land" is the normal seed-42 state; an arc on a thin band is consistent.
