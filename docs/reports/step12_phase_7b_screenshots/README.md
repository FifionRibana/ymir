# Step 12 Phase 7b.7 — UI screenshots capture procedure

Status: **scaffold + placeholders** (commit *Phase 7b.7*). The four
PNGs in this directory are 1280×720 placeholders rendered by
`scripts/generate_phase_7b_placeholders.ps1` so Phase 8 reports can
reference them with non-broken Markdown links. Real captures will be
swapped in by a follow-up commit (`DOC : capture Phase 7b.7
screenshots`) once the workflow panel is exercised in an interactive
viz session — the Bevy app is graphical and cannot be driven from a
headless agent loop.

## Expected captures

| File | Bridge state | What the screenshot proves |
|---|---|---|
| `01_panel_idle_pre_run.png` | `V2RunState::Idle` (post-init, no run yet) | Workflow panel is mounted under the parameter editor. `Run Phase A` enabled, `Stop` / `Continue` / `Run Phase B` / `Export HD heightmap` disabled. Step 12.X scope-limit caveat visible above the section headers. |
| `02_phase_a_running_cycle_3_of_5.png` | `V2RunState::Running` with `spec.workflow == V2WorkflowSpec::On` and `step == 3, total == 5` | Top-bar progress bar at 3/5; left dashboard streams `peek_state` metrics; right panel status badge "Phase A running — cycle 3/5 · Xs"; cycle-history table populated with 2 prior cycles' metrics (cycle 1, cycle 2). `Stop` enabled. |
| `03_phase_a_completed.png` | `V2RunState::WorkflowPhaseACompleted` | Status badge "Phase A done — 5 cycles in Xs" in light-green. Cycle-history table shows all 5 rows (erosion vol, sea level, mass drift, Δ craton). `Run Phase A`, `Continue`, `Run Phase B (HD)` all enabled. |
| `04_phase_b_hd_output_ready.png` | `V2RunState::WorkflowPhaseBCompleted` | Status badge "Phase B done — N×N HD in Xs · p95 = …" in light-blue. Left dashboard shows D5 metric line. `Export HD heightmap` enabled; if the user has clicked Export prior to capture, the success "Saved → …" line is visible. |

Resolution target: 1280×720 minimum (text on the panel must be
legible — this is the 320 px slider column rendered at default DPI;
2× HiDPI captures are also acceptable). PNG format.

## Reproduction procedure

The capture session below is the canonical sequence. Times in
parentheses are approximate wallclock at 64² × default cycle
parameters on a development laptop; actual numbers vary.

1. **Build + launch viz**
   ```powershell
   cargo run --release -p ymir-viz
   ```
   The window opens with the Tectonics phase active by default.

2. **Configure the workflow run**
   - Right panel → **Preset** dropdown → select `single_continent`.
     (`single_continent` produces a clean polygonal continent —
     visually unambiguous when documenting that Phase A preserves
     plate shape and Phase B HD adds the FBM-driven complexity.)
   - Right panel → expand **Workflow (Step 12)** section.
   - Toggle **Enable interleaved tectonic-erosion workflow** ON.
   - Confirm the D8 defaults are set (the toggle resets to
     `V2PhaseAParams::default()` + `V2PhaseBParams::default()`):
     `N_cycles = 5`, `k_cycle = 20`, `α = 0.01`, `β = 0`,
     `hd_grid_size = 2048²`, `droplets = 5×10⁶`.

3. **Capture #01 — `01_panel_idle_pre_run.png`**
   - Bridge state at this point is `Idle`; the workflow panel
     buttons reflect that (`Run Phase A` enabled, others disabled).
   - Screenshot the Bevy window with both the right panel and the
     left dashboard visible. The viewport sprite shows the
     `single_continent` initial state.

4. **Start Phase A**
   - Click **▶ Run Phase A**.
   - The bridge transitions to `Running` carrying
     `spec.workflow == On`. The top bar progress bar starts; the
     dashboard begins streaming `peek_state` metrics.

5. **Capture #02 — `02_phase_a_running_cycle_3_of_5.png`**
   - Watch the right panel status badge for "Phase A running —
     cycle 3/5". This typically lands ~10–20 s into the run on
     64² × 20 steps × 5 cycles. The cycle-history table at this
     moment shows 2 rows (cycle 1 and 2 already completed).
   - Screenshot.

6. **Wait for Phase A completion**
   - The bridge transitions to `WorkflowPhaseACompleted` (~30–60 s
     total). Status badge turns light-green; cycle history shows
     all 5 rows.

7. **Capture #03 — `03_phase_a_completed.png`**
   - Screenshot with the cycle-history table fully populated.

8. **Start Phase B**
   - Click **▶ Run Phase B (HD)**.
   - 30–90 s wallclock at 2048² × 5×10⁶ droplets. The bridge
     state stays `WorkflowPhaseACompleted` during the HD run
     (a follow-up could add a `WorkflowPhaseBRunning` variant; not
     in 7b scope), then transitions to `WorkflowPhaseBCompleted`.

9. **Click Export HD heightmap** (optional but documented)
   - The success line "Saved → C:\Users\…\AppData\Local\Temp\ymir_v2_phase_b_<ts>.png"
     appears below the buttons. Verify the PNG opens in an external
     viewer and shows the HD heightmap.

10. **Capture #04 — `04_phase_b_hd_output_ready.png`**
    - Screenshot with the Phase B status badge + Export success
      line visible. The viewport sprite still shows the Phase A S̃
      field (the v2 viewport is wired to LR rasters; HD display
      would require a separate sprite, future scope).

## Caveat for Phase 8 reports

Reports referencing this directory should include a footnote /
caveat:

> UI screenshots (manual capture): the Bevy viz cannot be driven from
> the CLI agent loop, so the captures in
> `docs/reports/step12_phase_7b_screenshots/` were taken by an
> interactive operator following the procedure documented in that
> directory's `README.md`. The procedure is reproducible and
> deterministic given the recorded preset + workflow defaults; only
> the wallclock-dependent status fields (elapsed time, exact step
> within a cycle) vary between captures.

## Placeholder lifecycle

The four `*.png` files in this directory currently render the text
"PENDING MANUAL CAPTURE" over a 1280×720 grey background. They were
produced by:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/generate_phase_7b_placeholders.ps1
```

The script is checked in under `scripts/` so the placeholders can be
regenerated if accidentally deleted. Real captures should overwrite
the `*.png` files in place; the README and procedure stay valid.
