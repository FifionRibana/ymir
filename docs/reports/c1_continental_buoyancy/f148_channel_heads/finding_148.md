# Finding 148 — rivers, extension 1 of 2: channel heads by convergence (measured; the stop rule fires, nothing built), Hack at the T3b junctions, the lake 1000011

**Status: measured. B is NOT built (the stop rule fired); J is built (the export); nothing committed.**

**Files:**
- predictions: `f148_predictions.md`, written before any measurement;
- instruments: `f148_declared.md`;
- raw output: `f148_rivers_raw.txt` (J, K, M) and `f148_diag_raw.txt` (the diagnoses; `f148_diag_raw_run1.txt` is the
  first run, without the catchment ratios);
- checks: `checks_before.txt` and `checks_after.txt`.

## Part 0

- F147 committed (`254c66f`, not pushed).
- **Recorded, proposed by the reviewer and not contested by the author:**
  - cut the receiving reach at the junctions that violate Hack's convention (the 1 : 1 match with `rivers.json` is
    dropped, as Living Landz reads `rivers_ll.json`);
  - amend the ridge instrument for zero-length corner touches, declared;
  - a = 5 by default;
  - the grid's radius and origin are verified on Living Landz's side, against the reference file.
- **Checks before** = F147's after-checks at `254c66f` (the same tree).
- **Checks after** (`checks_after.txt`):
  - `cargo check --workspace` clean;
  - lib 608 passed (605 + 3), 0 failed;
  - viz 30 passed;
  - **the guard 6 / 6, field and lakes = banc through `run_hd`** (C2 /10 col's field is the bench hash
    `a8d2d538d692c2f0`);
  - `rivers_ll` takes 0.24–0.44 s per state.

## G — how a river is born today (rules 11 / 11c)

**By area only:**
- `hd.rs:725–731` sets `head_km2` = A_c = 0.1 km² and `full_tree = false`;
- `flow.rs:1200–1255` (`extract_rivers`) makes a cell a river cell when:
  - its accumulation is ≥ 20 km² (`stream_km2`);
  - or it lies on the main stem of the largest sub-20 km² branch feeding such a cell, up to 0.1 km².
- No slope, no convergence.
- Erosion's own slope-dependent channel-head law (F56, `stream_power.rs`, Montgomery & Dietrich's A·S^α) does not
  enter the designation.
- The grep found no planform curvature or convergence index in the code (only test-side curvature benches, F84 /
  terrace diagnosis).

**The designation reaches more than the three consumers:**
- `resolve_exorheic_without_outlet` (`hd_assembly.rs:255`, Finding 86) relabels an Exorheic lake as Unresolved when
  no segment's source borders its footprint, so removing reaches can change a lake's type;
- `clip_rivers_to_lakes` and the spillway append work on the segments;
- the viz lake panel and the river microscope read them.
- The valley construction's skeleton (it works on the accumulation) and the climate (biomes) do not read the
  segments.

**Hence the declaration**: the head rule, as measured, lives in the `rivers_ll` export only and never touches
`dr.rivers`.

## J — Hack and the ridges (the export, built)

**J1**: `rivers_ll` now chains PIECES. Each segment is cut wherever a reach joins it below its head, and at each
junction the largest catchment continues.

| measure | F147 | F148 |
|---|---|---|
| rivers | 10 246 | 10 246 |
| T3b joins | 2 261 | 2 261 |
| of them, the joining stream larger (the receiver cut) | "831" (wrong, below) | **465** |
| Hack violations (a tributary larger than its receiver above the junction by > 0.1 %, per-vertex catchments) | — | **0** of 7 984 confluences |
| crossings · through a shared run · self-intersections | 0 · 0 · 0 | **0 · 0 · 0** |
| lateral p99 | 0.133 cell | 0.133 cell |

- **F147's 831 was a units error.**
  - The témoin runs at a geographic scale ratio of 7.5. `segment_catchment_cells` is signified (× 56.25, measured
    on all 16 226 segments: p1 = p99 = 56.25), while `flow.accumulation` is real.
  - F147 compared one to the other. F148 scales the accumulation to each segment's catchment, so both sides are
    signified.
  - The ADR's F147 line is corrected in place.
- The format goes to **0.4.0**: `segments` lists the segments a river covers wholly or in part.
- A permanent test with a negative control: the larger joiner continues; the smaller does not.

**J2**: the ridge instrument amended as declared. It cuts each final edge exactly at the cell boundaries, and a
crossing is a stretch of length > 0.001 cell inside an off-outlet cell.
- F147's sample instrument, same build: **23**.
- Amended: **134**.
- **The amendment did not do what it was meant to.** It removes the zero-length corner touches, but it also sees the
  short stretches the 0.25-cell sampling stepped over.
- The 134 stretches (`f148_diag_raw.txt`):
  - length p10 / p50 / p90 / max = 0.003 / 0.012 / 0.042 / **0.098 cell** (max 4.8 m);
  - all ≤ 0.1 cell; total 2.5 cells (122 m) over 21 541 km.
- The smoothed line clips the corners of neighbouring cells that drain elsewhere.
- **A length tolerance would be a second amendment, the reviewer's to decide.** Not made.

## K — the lake 1000011 (attribution)

| what | measured |
|---|---|
| the lake | Unresolved, reason `NoOutletReach` (F86), 494 km², level 459.3 m; id ≥ 1 000 001: a below-sea-path basin |
| its outlet | **traced**: spillway segment 16217, `source_lake = 1000011`, 33 points, **1 625 m³/s**. Its first point (4424, 3975), at 459.32 m, is the col; D8 runs away from the lake (459.3 → 454.4 m in 5 steps) |
| the gap | between the col and the footprint lies a rim cell at **459.3 m, the lake's own level**, outside the footprint. So the spillway starts **2 cells** from the shore |
| the core's check | `exorheic_lakes_missing_outlet` accepts a source within **1** cell (the 8-neighbourhood) → no outlet → `resolve_exorheic_without_outlet` relabels it Unresolved |
| the panel's search | `workspace.rs:5317–5387` accepts a segment source within **2** cells → finds the spillway's watercourse ("#3"). It calls the lake "exoréique" because it tests only `Endorheic` |

**At which stage is the outlet lost?** None. It is traced and exported. The **core's check misses it**, by a rim cell
at the lake's level.

**Why does the panel say the opposite of the log?**
- It searches farther (2 cells against 1);
- and it treats Unresolved as not-closed.

**Here the panel is right and the core's label is wrong.**

**Not a false display, so nothing is corrected.** A display fix would show the core's false verdict faithfully.

**Queued, a named defect:** `exorheic_lakes_missing_outlet` misses a traced spillway whose col sits behind a rim cell
at the lake's level.
- A candidate fix: count a segment tagged `source_lake = id`, or a source within 2 cells.
- In this world it is the only Unresolved lake.

## M — convergence against A·S² (the reviewer's hypothesis, measured)

**The index**: C = div(∇z/|∇z|) on the conditioned field, Gaussian σ = 2 cells. Unit-tested on a pit, a peak, a
valley, and a plane as the negative control.

**Populations (the J build):**
- (a) **873** fusion candidates (69 in the control, removed);
- (b) **705** micro-rivers into a lake < 1 km (21 in the control);
- (a) ∩ (b) = 0;
- **(c) 2 828** control rivers: 255 ≥ 10 km; 2 695 ≥ 50 % on a construction floor.

**AUC** = P(score(c) > score(a ∪ b)):

| field | max | head (first 10 cells) | median | head median (a) / (b) / (c) |
|---|---|---|---|---|
| convergence σ = 1 | 0.796 | **0.857** | 0.837 | −0.11 / −0.05 / +3.20 km⁻¹ |
| **convergence σ = 2 (declared)** | **0.736** | **0.784** | 0.753 | −0.12 / −0.07 / +1.26 km⁻¹ |
| convergence σ = 4 | 0.713 | 0.748 | 0.722 | −0.13 / −0.09 / +0.58 km⁻¹ |
| A·S² | 0.462 | 0.418 | 0.306 | 639 / 966 / 292 m² |

- **The heads of (a) and (b) are planar to slightly divergent; the control's are convergent.** The reviewer's
  hypothesis holds at the head.
- **A·S² separates the wrong way**: the bundles and the micro-rivers sit on STEEPER ground than the control.

**The rule simulated** (σ = 2; a reach is channel from its first cell with A ≥ 0.1 km² and C ≥ θ, lake outlets and
spillways exempt; downstream closure; J build: 10 246 rivers, 21 541 km, 1 660 pairs, 726 micro-lake rivers):

| θ (km⁻¹) | rivers | km | pairs | micro-lake | (a) removed | (b) removed | (a)∪(b) | **(c) lost** |
|---|---|---|---|---|---|---|---|---|
| (c) p0.5 = −1.24 | 10 223 | 21 249 | 1 636 | 732 | 0.8 % | 4.4 % | 2.4 % | **0.99 %** |
| (c) p1 = −0.95 | 10 199 | 21 106 | 1 615 | 729 | 0.9 % | 6.5 % | 3.4 % | 1.87 % |
| (c) p2 = −0.68 | 10 131 | 20 869 | 1 593 | 721 | 2.7 % | 11.2 % | 6.5 % | 3.36 % |
| (c) p5 = −0.29 | 9 913 | 19 710 | 1 375 | 740 | 14.3 % | 24.3 % | 18.8 % | 7.85 % |
| (c) p10 = −0.03 | 9 454 | 17 679 | 950 | 660 | 34.7 % | 49.2 % | 41.2 % | 14.85 % |
| 0 | 9 400 | 17 448 | 915 | 655 | 36.8 % | 51.2 % | 43.2 % | 16.27 % |
| (c) p20 = +0.30 | 8 676 | 14 916 | 563 | 497 | 61.1 % | 72.5 % | 66.2 % | 29.74 % |

A·S² does worse at every θ (at its p0.5, 0.0 % of (a)∪(b) for 0.92 % of (c); at its p20, 12.5 % for 26 %).

**STOP RULE: no θ removes ≥ 50 % of (a) ∪ (b) while losing ≤ 1 % of (c), for either criterion. Reported; B not
built** (no gated rule in the viz, no flag, no display toggle), so part P does not apply.

**Why** (a hypothesis to check, not measured): **16 % of the control has no convergent cell at all along its trace**
(θ = 0 loses 16.3 %).
- 2 695 of the 2 828 control rivers lie ≥ 50 % on a construction floor. The construction's floors are planar strips
  (W(A) from 178 m up to 2.9 km), where C ≈ 0 along the centreline.
- So the rule's statistic (one eligible cell along the trace) cannot tell a river in a wide built valley from a rill
  on a planar slope, even though their HEADS differ (AUC 0.78–0.86).
- In the simulation, the micro-lake count even rises at some θ (726 → 740): trimming shortens rivers below 1 km.

**The cost** if built: the convergence field takes **0.17 s**; the build with the rule 0.36 s (0.39 s without). That
is +0.17 s against 249.8 s.

## Predictions

**The reviewer's (hypotheses to check):**
- "AUC > 0.85, better than A·S²":
  - at the declared σ = 2, REFUTED (max 0.736, head 0.784);
  - at σ = 1, the head score reaches 0.857 (a sensitivity, not the decision);
  - "better than A·S²" HELD (A·S² < 0.5).
- "A threshold removes ≥ 70 % of the pairs without losing a control river": REFUTED (1.5 % of the pairs at 0.99 %
  control loss).
- G-micro > 80 % and G-lakes: not reached (B not built).
- **K, "the panel reads an outlet declared by the balance that the trace did not produce"**: REFUTED. The trace
  produced it (1 625 m³/s); the core's 1-cell check misses it.
- Cost < 2 s: HELD (+0.17 s).
- Meta: HELD.

**Mine:**
- **P-J**:
  - Hack 0, crossings 0, river count < 1 %: HELD;
  - "amended ridges ≤ 3": REFUTED (134).
- **P-M**:
  - "max AUCs > 0.85 by a length confound": REFUTED (0.74 / 0.46; the confound did not dominate);
  - "head AUC 0.60–0.75": REFUTED (0.78);
  - **"the bundles are incised channels, not sheet flow": REFUTED at the head** (median −0.12 against +1.26);
  - "< 30 % of the pairs at ≤ 1 % control loss": HELD (1.5 %);
  - "(b) removed > 60 %": REFUTED (4.4 %);
  - "the stop rule does not fire": REFUTED.
- **P-K**:
  - "a display defect": REFUTED (a core false negative);
  - "the 2-cell against 1-cell search": HELD;
  - "the stage that loses it is the trace or the clip": REFUTED (nothing loses it).
- **P-cost, 0.5–2 s**: REFUTED (+0.17 s).
- **Meta**: HELD.

## For the author

**Round 2 of the extension is yours. Proposals, none built:**
1. **The display toggle on its own** (« Masquer faisceaux et micro-rivières »).
   - It needs the F147 flags plus a `micro_lake_inflow` flag; it does not depend on any head rule.
   - The stop rule kept it out of this round.
2. A head rule on the **head score** (the mean over the first cells), where the separation lives (AUC 0.78–0.86),
   instead of one eligible cell anywhere. To measure with the same stop rule.
3. The K fix in the core (queued above).
4. The ridge instrument's length tolerance (≤ 0.1 cell), if the reviewer wants it.

**Then the geology.**
