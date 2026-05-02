# ymir-viz

Interactive Bevy + bevy_egui visualisation for Ymir's `tectonics_v2`
solver. Run a thin viscous sheet simulation, watch the rasters
update live, capture screenshots, and round-trip completed runs to
JSON for offline review.

After Step 8.6 Phase 8h sunset this binary is **v2-only**. The
legacy `tectonics::` bridge and every plugin / panel that drove its
pipeline phases have been removed.

## Run

```bash
cargo run -p ymir-viz --release
```

Always use `--release` — debug builds are 10–20× slower on the
solver thread.

The window opens with the v2 top bar (status badge + step counter
+ progress bar + wallclock), a left-side metrics dashboard, a
right-side parameter editor, and the centre sprite that paints the
currently-selected raster field.

A single env knob is read at startup:

| Var | Effect |
|-----|--------|
| `RUST_LOG` | Bevy / `bevy_egui` log filter override (defaults to `warn,ymir_core::tectonics_v2=info,ymir_viz=info`). |

## Workflow

### 1. Pick a preset

Six presets ship under `presets/v2/`:

| Preset            | What it demonstrates |
|-------------------|----------------------|
| `quiescent`       | Mantle off + boundary off + perturbation amplitude reduced. Steady, near-zero velocity — sanity / regression baseline. |
| `single_continent`| One large continental plate surrounded by oceanic. Useful for craton-immunity inspection. |
| `convergence`     | Multiple plates with strong convergent forcing — yielding bands at collision interfaces. |
| `divergence`     | Strong divergent forcing — yielding traces rifting boundaries. |
| `subduction`     | Subduction-flavoured boundary fluxes — pronounced ε̇_II near subduction edges. |
| `active_medley`   | Step 8 mantle on + cratonic on + age on; the §4.11 validated regime. **Default.** |

Pick from the dropdown in the right panel. The dropdown clobbers
the editable spec but preserves the user's `output_dir` and
`capture_endpoints` choices.

### 2. Choose an init mode (Phase 8a)

The "Initialisation (S̃)" section of the right panel exposes four
init modes for the crustal-thickness field:

| Mode             | Pattern |
|------------------|---------|
| `Checkerboard`   | Legacy sinusoidal perturbation around per-plate-type means (oceanic 0.2, continental 1.0). Required for Steps 0–10 numerical regression baselines; not recommended for visual review (introduces global bumps not present in the underlying physics). |
| `Uniform`        | **Default.** Flat per-plate-type with smoothstep transition over `boundary_smoothing_width` cells across inter-plate edges. Aligns with TDD §4.2. Recommended for visual review and for any visualisation downstream of the milestone. |
| `Gaussian`       | Per-plate Gaussian peaked at the Voronoï seed coordinate, decaying with periodic minimum-image distance. `sigma_continental` and `sigma_oceanic` measured in cells. Useful for emphasising plate centres. |
| `Convolution`    | Periodic Gaussian blur of the binary classification mask. `sigma` in cells. Smooth everywhere; loses sharp inter-plate edges. |

`Uniform` is the recommended default for both research and
illustration. `Checkerboard` is the regression-preserving fallback;
the Steps 0–10 lib tests under `crates/ymir-core/tests/v2_*.rs`
opt into it explicitly.

### 3. Tune physics knobs

Collapsible sections in the right panel expose:

- **Grid & seed** — resolution (32² / 64² / 128²), seed, plate count
  (3–15), continental ratio, total step count.
- **Yielding & drag** — Bi (yield), Br (basal drag).
- **Initialisation (S̃)** — init mode + per-mode parameters
  (Phase 8a).
- **Mantle** — Mf, modes, mantle seed, evolution rate.
- **Cratonic immunity** — Cr, K (viscous), B_factor, smoothing
  width, plate area minimum.
- **Age field** — continental / oceanic init ages.
- **Solver options** — slab toggle (forward-compat), linear solver
  pick (JacobiCG / AmgCG), force kind (GPE / Sinusoidal).

### 4. Run

Hit ▶ Run. The bridge thread builds the full `BaselineConfig` from
the spec, spawns the harness, and streams `Started` / `Progress` /
`Completed` events back to the main thread. The metrics dashboard
on the left fills in live (peak |v|, ⟨S̃⟩, max ε̇_II, cratonic
fraction) and turns into the full final-state summary on
completion.

⏹ Cancel signals the worker thread to bail at the next step
boundary; the run state stays at the last `peek_state` so the
sprite doesn't blank out.

### 5. Display options

The "Display" section of the right panel:

- **Field selector** — S̃ / Age / Cratonic factor / ε̇_II / |v|.
  Hovering each option surfaces a one-line legend caption; the
  legend bar below renders 32 stops of the matching colormap.
  S̃ uses a fixed `[0, 2.5]` band; ε̇_II uses log `[1e-3, 1e2]`;
  |v| uses log `[1e-5, 1e1]`; Age uses dynamic min/max; Cratonic
  is fixed `[0, 1]`.
- **Voronoï boundaries** (Phase 8b) — paints inter-plate edges in
  black on top of the field colour. Default off.
- **Velocity vectors** (Phase 8b) — one yellow arrow per plate at
  the periodic-aware centroid; length proportional to mean
  per-plate velocity. Quiescent plates (`|v̄| < 1` cell)
  skipped. Default off.

### 6. Capture (Phase 6)

📷 Capture saves the currently-selected field to PNG under
`<output_dir>/screenshots/`. Filename format:
`{preset_label}_{field}_{unix_ts}.png`. The status line below the
button confirms the path.

If you want overlays in the saved PNG, the dedicated test
`v2_phase8g_visuals` already burns Voronoï + velocity into every
captured frame; the bare Capture button writes raw colormap.

### 7. Export / import a run (Phase 8e)

💾 Export run saves the current `Completed` run to JSON under
`<output_dir>/snapshots/{preset_label}_{unix_ts}.json`. The
schema:

```json
{
  "format_version": 1,
  "exported_at": "2026-04-30T12:34:56Z",
  "elapsed_seconds": 42.3,
  "spec": { /* full V2RunSpec */ },
  "scalar_metrics": { /* dashboard-relevant scalars */ },
  "final_state": { /* every raster field */ }
}
```

The export populates the import path field below the button so a
quick Export → Import round-trip is one click each side.

📂 Import run reads a JSON path from the text field and replaces
the bridge state with `V2RunState::Imported`. Renders identically
to a fresh `Completed` run from the sprite's perspective. Useful
for offline review without re-running the solver, or for
regression baselines (a snapshot from a known-good revision can
be reloaded and visually compared against a current run).

The loader rejects unknown `format_version` — bumps reserved for
schema changes. Forward-compat probe lives in the integration test
`v2_bridge_export_import_roundtrip`.

## Programmatic harness reproduction

The integration tests under `crates/ymir-viz/tests/v2_*.rs` can be
run with `--ignored` for the longer ones:

```bash
# Phase 8f equilibrium analysis (~26 min at 32²; bump grid via
# YMIR_PHASE8F_GRID=64 for ~60–90 min canonical reference).
cargo test --release -p ymir-viz --test v2_phase8f_equilibrium \
    --jobs 1 -- --ignored --nocapture

# Phase 8g visual revalidation (4 presets × 32² × 100 steps,
# ≈ 42 min total).
cargo test --release -p ymir-viz --test v2_phase8g_visuals \
    --jobs 1 -- --ignored --nocapture

# Phase 7 follow-up frame-by-frame strip + patchwork composition
# (per-preset diagnostic, env-driven).
cargo test --release -p ymir-viz --test v2_phase7_step_diagnostic \
    --jobs 1 -- --ignored --nocapture
cargo test --release -p ymir-viz --test v2_phase7_patchwork \
    --jobs 1 -- --ignored --nocapture
```

Override knobs via env vars:

| Var | Default | Effect |
|-----|---------|--------|
| `YMIR_DIAG_PRESET`   | `active_medley` | Preset for `v2_phase7_step_diagnostic`. |
| `YMIR_DIAG_STEPS`    | `30`            | Step count. |
| `YMIR_DIAG_INTERVAL` | `1`             | Capture every Nth step. |
| `YMIR_DIAG_GRID`     | `32`            | Grid edge. |
| `YMIR_PATCHWORK_PRESETS` | (all) | Comma-list of presets to re-tile. |
| `YMIR_PHASE8F_GRID`  | `32`            | Equilibrium-test grid edge. |
| `YMIR_PHASE8G_GRID`  | `32`            | Visual-revalidation grid edge. |
| `YMIR_PHASE8G_STEPS` | `100`           | Visual-revalidation step count. |
| `YMIR_PHASE8G_INTERVAL` | `10`         | Capture interval in steps. |

## Layout

```
crates/ymir-viz/
├── presets/v2/                  — six preset JSONs (loaded via include_str!)
├── src/
│   ├── main.rs                  — Bevy app entry (v2-only post-sunset)
│   ├── lib.rs                   — library facade (re-exports for tests)
│   ├── camera.rs                — pan / zoom / cursor world position
│   ├── bridge/
│   │   └── v2/                  — solver bridge (commands, events, plugin,
│   │                              spec, snapshot, presets, build_config)
│   ├── visualization/
│   │   ├── colormap.rs          — hypsometric / age / log_hot / cratonic
│   │   ├── overlay.rs           — Voronoï boundaries + velocity vectors
│   │   └── v2_viz.rs            — sprite update + screenshot path
│   └── ui/
│       ├── parameter_panel_v2.rs— right-side editor (Phase 8d sliders)
│       └── metrics_dashboard.rs — left-side dashboard (Phase 8c)
└── tests/                       — integration tests, all `v2_*.rs`
```

## Reports

The Step 8.6 augmenté milestone produced these artefacts under
`docs/reports/`:

- `step8_6_viz_portage_report.md` — closing report (this milestone).
- `step8_6_phase7_gallery/` — Phase 7 reviewer-validated checkpoint
  (the gallery that surfaced the ε̇_II + init-mode bugs).
- `step8_6_phase8f_equilibrium/active_medley_32sq.md` — equilibrium
  verdict.
- `step8_6_phase8g_visuals/` — post-corrections visual revalidation
  (24 patchwork PNGs + REPORT.md).

## License

Proprietary — same terms as the workspace root.
