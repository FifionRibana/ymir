# Finding 148 — predictions, written 2026-10-06 BEFORE any measurement of this round

**Reading declaration, unfavourable.** Non-blind on Findings 73–147, the reviewer's predictions, and the code read for
parts G and K before writing this.

**G, read in the code:**
- Production designates a river cell by AREA ONLY (`hd.rs:725–731`, `flow.rs:1200–1255`):
  - every cell with an accumulation ≥ 20 km² (`stream_km2`);
  - plus, for each such cell, the main stem of its largest sub-20 km² inflow, extended upstream to A_c = 0.1 km²
    (`head_km2`, `full_tree = false`).
- No slope, no convergence. Erosion has its own slope-dependent channel-head law (F56, `stream_power.rs`); the river
  designation does not use it.
- **The designation reaches the lakes**:
  - `resolve_exorheic_without_outlet` (hd_assembly.rs:255) relabels an Exorheic lake with no reach starting on its
    shore as Unresolved;
  - `clip_rivers_to_lakes` and the spillway append work on the segments.
- The valley construction's skeleton and the climate do not read the segments.

**K, read in the code**: the viz lake panel (`workspace.rs:5277–5419`) differs from the core in two ways:
- it calls a lake "closed" only when it is `Endorheic`, so an `Unresolved` lake reads "exoréique · Exutoire : oui →
  aval";
- it searches the outlet reach within **2** cells of the shore, where the core's check uses **1**.

## Predictions

- **P-J**:
  - by construction 0 violations of Hack's convention; the crossings and self-intersections stay at 0;
  - the amended ridge instrument (the length of the polyline in an off-trace cell of another outlet, a zero-length
    corner touch excluded) gives **≤ 3** where F147's sample instrument gave 23 (the 2 samples off the corners by
    0.004–0.037 cell may be real millimetre passes);
  - the river count moves by < 1 %.
- **P-M**:
  - **a confound first**: a per-river MAXIMUM of any cell index along the trace grows with the trace's length, and the
    control (c) is long by selection. So the max-score AUCs of (a)+(b) against (c) will be high for a trivial reason:
    **> 0.85 for both convergence and A·S²** (A grows with length too).
  - The head-reach score, which does not see length: AUC **0.60–0.75** for convergence.
  - **The bundles (a) are NOT sheet flow**: F146 put 92 % of the pairs between 10–100 km² rivers, 88 % on uncarved
    terrain; they are incised channels on a regional tilt, with convergence. **No threshold that loses ≤ 1 % of (c)
    removes ≥ 50 % of the parallel pairs (I expect < 30 %)**, against the reviewer's 70 %.
  - **The micro-rivers into lakes (b) are removed easily** (flat shores, divergent): > 60 % at ≤ 1 % control loss.
  - Whether the stop rule fires depends on the sizes of (a) and (b). I predict it **does not fire** (the (b) removal
    carries (a)+(b) past 50 %), but G-parallel fails if B is built.
- **P-G-lakes**: holds bit for bit, because the rule lives in the `rivers_ll` export and never touches `dr.rivers`
  (see G).
- **P-K**:
  - a display defect: the panel shows the Unresolved lake 1000011 as exorheic with an outlet, because of its own
    closed-or-not test and its 2-cell search;
  - the "#3" it names is a reach rising 2 cells off the shore, or one whose source touches the lake on the other
    side;
  - the core's stage that loses the outlet is the below-sea spillway trace or `clip_rivers_to_lakes`.
  - Against the reviewer's "two sources" in part: the "exorheic" label does come from the balance, but the outlet the
    panel shows is its own geometric search.
- **P-cost**: convergence field (a Gaussian at σ = 2 cells on 8192² + a divergence) plus the rule: **0.5–2 s per
  world**. With the reviewer (< 2 s).

## Meta

My disagreements with the reviewer: the bundles' nature (M, G-parallel), K's mechanism. **At least one of my
predictions is wrong.**
