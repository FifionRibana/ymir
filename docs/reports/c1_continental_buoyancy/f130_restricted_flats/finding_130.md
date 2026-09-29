# Finding 130 — the restricted clause returns the raised cells to the floor but not θ to its CI; the parabolic valley's own D8 tree does not give its D8 errors; the skeleton's flats are the lakes', and resolving them makes LTD cardinal

**Nothing is promoted; no clause is activated.**
- The benches are in `crates/ymir-core/tests/f126_coast.rs`, Finding 130's section.
- The raw outputs are beside this report.
- mur ↔ mer is ON in every world. bruit ↔ pied is ON where there is a wall profile (B2 + B1 + foot).
- The hydrology is D8 everywhere.

**Units.**
- 1 cell = 48.8 m (domain), 8192² over 400 km.
- Heights are in metres. Areas are domain km². Distances are in m or in W.

**Reading declaration, unfavourable.** My predictions (`f130_predictions.md`) were written before any grep,
reading or measurement. They are non-blind on Findings 73–129 and on the reviewer's predictions. Every
definition and gate was declared before its run:
- B0's geometry, q and tolerance, after the reading and before any B run;
- A, P and E, before any of their runs.

After B failed, one diagnostic was declared before it ran, and it is reported, not gated.

## 0 — Finding 129 committed; the guard; the greps

- **Finding 129 was committed** as `a35f67b`. The guard through `run_hd` reads **6 / 6 "= banc"** (`guard_six_states.txt`).
- **This round's library code** is gated and bit-neutral when off (see "State"):
  - `ValleyConstruction::trunk_band_downstream` (A2);
  - `Skeleton::line_parent` (each line's parent line and its junction sample);
  - `carve_diag` (what `carve` decided per cell);
  - `ltd_directions_masked` (the LTD core on `f32` or `f64` elevations, with a closed-boundary mask; `f32`
    is bit-identical);
  - `flow::flat_resolution` (`resolve_flats`, exposed with its flat set).
- **A permanent test with a negative control**:
  `the_restricted_clause_takes_the_band_only_downstream_of_the_junction`.
  - Control: the full clause re-lays cells upstream of the junction.
  - Restricted: those cells are bit-identical to no clause, and the cells downstream equal the full clause's.

**The greps (rule 11/11b/11c), earliest hit read.**
- **`FlatPerturbation`** first appears in Finding 87: "a flat-routing and quantisation item — `FlatPerturbation`
  and D∞ exist in the config for this class — and it is named, not fixed".
- **`resolve_flats`** first appears in Finding 74, where it was recorded as NOTHING FOUND. **`Garbrecht`**
  appears only in Finding 128.
- **Barnes** first appears in Finding 13, as the priority-flood conditioning.
- **Finding 79 C: NOTHING FOUND** under that name. Finding 79's section "C / C-bis" (L9229) is the
  deposition class, "refuted as a short path" for the coastal fringe.
- `alluvial` appears once, in Finding 111's rule-11 list. `min(field` appears only in Finding 129.
- **The answer to the round's question.** The repository HAS a flat resolution: `resolve_flats`
  (Garbrecht–Martz 1997, `flow.rs:581`), with `FlatPerturbation` on top.
  - It is applied to the field the skeleton reads, because `skeleton` calls `compute_flow` with the
    drainage chain's `flat_perturbation`.
  - But it resolves D8's **pointer**, through a routing surface `flat_grad`, and never the **elevation**.
  - LTD reads the elevation, so it sees the lakes' surfaces as flats and falls back (P0).

## B — [O03]'s parabolic valley, as its figures define it: B0 fails on two pointers, B1 fails; B stops

**B0 — what the article defines.**
- **The text does not define the valley.** §3.1 [11] (p. TNN 1-3) gives only "a parabolic valley", its 3-D
  view in Fig. 2d, and "theoretical areas … obtained using the analytical integrals of the drainage lines
  through the corners of the draining cells". There is no equation, no dimensions and no spacing.
- [9] (p. TNN 1-3): "the D8-LAD method with λ = 0 reproduces the classical D8 method".
- **Its figures define more.** The PDF is vector, and I read its coordinates:
  - Fig. 2h / 2l draw a **25 × 25** grid (625 cell outlines).
  - D8's pattern is identical in every column, so the down-valley slope is uniform, and it is symmetric
    about row 12.
  - D8's diagonals are exactly rows 7–10 and 14–17, columns 0–23 (192 cells). Rows 11–13 flow east, and
    every other cell flows toward the axis.
  - There is one outlet, east of (24, 12), and a closed boundary elsewhere.
- **The parameter q is bounded by D8's figure.** For z = s·x + a (y − 12.5)², D8's discrete criterion on
  drops makes that pattern hold iff q = a/s ∈ (0.2195, 0.268).
- **Declared, before any run:**
  - q = 0.244, the midpoint, chosen from D8's figure only; 0.2205 and 0.2670 are reported as sensitivity;
  - the analytic flow-tube areas (C = (y − yc)·e^{2qx} is constant on a drainage line);
  - F129's traced areas beside them.
- **Fig. 3j–l, read from the vector markers:**
  - D8 (λ = 0): ME −0.038, MAE **0.566**, RMSE **2.043**;
  - LTD (λ = 1): ME −0.261, MAE **0.282**, RMSE **0.330**;
  - the plane's markers read 0.504 / 1.907 as a check, as in Finding 128.
- **Tolerance**: ± half a marker (ME ± 0.017, MAE ± 0.034, RMSE ± 0.064).

| run | D8 vs Fig. 2h | D8: ME / MAE / RMSE | LTD: ME / MAE / RMSE | LTD's diagonals vs Fig. 2l | verdict |
|---|---|---|---|---|---|
| **as declared, q = 0.244** | **623 / 625** | −0.030 / **0.453** / **0.952** | −0.234 / 0.267 / 0.314 | 204 of 252 | **B0 FAILS · B1 FAILS · LTD NOT JUDGED** |
| q = 0.2205 (sensitivity) | 623 / 625 | −0.054 / 0.450 / 0.890 | −0.231 / 0.272 / 0.318 | **252 of 252** (ours 256) | — |
| q = 0.2670 (sensitivity) | 623 / 625 | −0.007 / 0.457 / 1.012 | −0.199 / 0.247 / 0.298 | 148 of 252 | — |
| diagnostic, one-inflow outlet, q = 0.244 | **625 / 625** | −0.028 / **0.452** / **0.951** | −0.234 / 0.267 / 0.314 | 204 | not a gate |

- **B0 fails on MY outlet.** It admitted diagonal inflow, so (24, 11) and (24, 13) went into it; the figure
  sends them to (24, 12) (`b_where_d8_leaves_fig2h.txt`).
- **B1 fails.** D8's MAE is 0.453 against 0.566 ± 0.034, and its RMSE 0.952 against 2.043 ± 0.064.
- **By the round's rule, B stops.** "La géométrie n'est toujours pas la sienne": LTD is not judged, and C
  does not run.
- **The diagnostic** was declared before it ran, and is not a gate. With the outlet admitting (24, 12) only,
  our D8 reproduces Fig. 2h's **625 pointers**, and its errors do not move (0.452 / 0.951).
  - **Our D8 tree is the paper's D8 tree, and it does not give the paper's D8 errors.** The difference cannot
    be in the tree, and our analytic and traced areas agree within 0.001.
  - So it lies in the paper's theoretical areas, which it describes in one sentence and does not publish.
- **As the article publishes itself, its D8 is not reproducible.** Finding 128's plane showed the same
  thing: our D8 errors below the paper's while LTD's absolute values matched.
- **Seen, not used.**
  - LTD's absolute errors at the declared q fall inside the tolerance of the paper's (0.267 / 0.314 against
    0.282 ± 0.034 / 0.330 ± 0.064).
  - At q = 0.2205 LTD contains all 252 of Fig. 2l's diagonals.
  - Neither is a verdict: LTD is not judged when D8 is not reproduced, and q was not re-chosen.

## P — the flats of the field the skeleton reads

**P0 — where the 6.05 % are (`p_flats.txt`).** The population is S1 breached: land cells with no strictly
lower neighbour, 684 893 = Finding 129's count.

| exclusive class | cells | share |
|---|---|---|
| a rounding tie of the metres conversion | 0 | 0.0 % |
| **on a pre-incision lake footprint** | **684 790** | **100.0 %** |
| a cell the breach changed | 85 | 0.0 % |
| coastal (≤ 976 m of the sea) | 0 | 0.0 % |
| other | 18 | 0.0 % |

- **They are the surfaces of the pre-incision's lakes**: 36 connected flats. The breached field keeps the held
  lakes' filled surfaces flat.
- **All 684 893 are in `resolve_flats`' flat set.** Garbrecht–Martz saw them and resolved D8's pointer
  there.
- **That pointer is CARDINAL on 20.1 % only.** `FlatPerturbation` makes it mostly diagonal. Finding 129's
  "the cardinal Garbrecht–Martz pointer" was wrong, and is corrected here.
- **Why the flats remain for LTD**: the resolution never touches the elevation, and LTD reads the elevation.

**P1 — the Garbrecht–Martz gradient laid into the elevation, on the skeleton's field only.**
- The method is `resolve_flats`' own `flat_grad`, with the chain's `FlatPerturbation`, added per flat to
  `filled` in f64. It is bounded below a tenth of the flat's smallest rise.
- The hydrology's D8 is untouched.

| | reference (LTD on `filled`, no increment) | **P1** |
|---|---|---|
| flat fallbacks | 684 893 (6.046 % of land) | **0 (0.0000 %)** |
| CARDINAL pointers on the former flats | 137 539 (20.1 %) | **666 441 (97.3 %)** |
| non-flat cells, steepest facet unchanged | — | 99.971 % |
| non-flat cells, LTD pointer unchanged | — | 99.959 % |

- **The flats fall to zero, and C's condition (< 0.1 %) is met.**
- **But on the 684 893 lake-surface cells LTD then goes cardinal**: 97.3 %, against 20.1 % for the perturbed
  D8.
  - The Garbrecht–Martz field is an 8-connected BFS distance, whose gradient is cardinal. Its D∞ use was
    rejected for the same reason (`FlowConfig::dinf`, `flow.rs:45-51`).
  - Read on it, LTD reproduces the lattice it was meant to leave.
  - **P1 solves the fallback count and worsens the thing the fallbacks stood for.** Named, not treated.

## C — NOT RUN

C requires B2 to pass and P1 < 0.1 %. P1 < 0.1 % holds. **B2 was not judged (B1 failed).** No gate was
declared, and nothing of C was measured. `f129_c_refs` was not run.

## A — the restricted clause, on the trunks

### A0 — where the full clause's raised cells are (`a0_a1_a2_law.txt`)

The population: trunk ≥ 10 km² cells more than 1 m above their law with the full clause and not without it,
**1 252 cells** (Finding 129 counted the difference of totals, 1 209).

| | cells | share |
|---|---|---|
| laid by a larger line's sample | 1 252 | 100 % |
| the owner's chain never meets the larger line (NO junction) | 4 | 0.3 % |
| **UPSTREAM of the junction** | **1 247** | **99.9 %** of those with one |
| at or downstream | 1 | 0.1 % |
| **within 1 W of the junction** | 1 246 | 99.8 % |

- **The distance along the larger line is p10 / p50 / p90 19 / 24 / 72 m, i.e. 0.05 / 0.10 / 0.24 W.** The
  laying sample is the larger line's FIRST sample above the junction.
- **Its floor is higher than the junction's in 99.9 %**: p50 +4.3 m (p10 +2.3, p90 +8.5). It is +3.2 m (p50)
  above the cell's own law.
- **⇒ The deformation is the larger line laying its floor from a point just UPSTREAM of the junction**, i.e.
  higher, over the tributary's last cells.

### A1 — the references, with a link bootstrap (1 000 resamples, 95 %)

| world | raised > 1 m above the law | θ laid [95 %] |
|---|---|---|
| **témoin C2/10 col** | 82 | **0.493** [0.492, 0.495] |
| B2 → A_c, no clause (the floor) | **86** | **0.493** [0.491, 0.496] |
| B2 → A_c, the full clause | 1 295 | 0.469 [0.464, 0.473] |
| B2 (A_c) + B1 + foot, no clause | 66 | 0.492 [0.489, 0.495] |
| B2 (A_c) + B1 + foot, the full clause | 1 278 | 0.468 [0.464, 0.472] |

The i.i.d. bootstrap ignores the links' spatial correlation, so these intervals are a LOWER bound on the
uncertainty.

### A2 — the restricted clause

| world | raised > 1 m | θ laid [95 %] | over-dug (canyons) | col 7 | col 14 | 13's dam cell | canyon 13 | lakes |
|---|---|---|---|---|---|---|---|---|
| **témoin C2/10 col** | 82 | **0.493** [0.492, 0.495] | **0** | **623.6 m** | **119.6 m** | **299.7 m** | — | 25 |
| B2 → A_c, no clause | 86 | 0.493 [0.491, 0.496] | 4 (3, 7, 13, 14) | 787.6 m | 493.7 m | 546.2 m | dammed 52 cells down | 33 |
| B2 → A_c, full clause | 1 295 | 0.469 [0.464, 0.473] | 0 | 714.2 m | 231.6 m | 362.8 m | gone | 25 |
| **B2 → A_c, RESTRICTED** | **86** | **0.485** [0.482, 0.488] | **4 (the same bodies)** | **787.5 m** | **493.7 m** | **546.2 m** | **dammed, 52 cells down** | 33 |
| B2 (A_c) + B1 + foot, no clause | 66 | 0.492 [0.489, 0.495] | 3 (Finding 127) | — | — | — | — | — |
| B2 (A_c) + B1 + foot, full clause | 1 278 | 0.468 [0.464, 0.472] | 0 | 704.7 m | 250.7 m | 351.6 m | gone | 26 |
| **B2 (A_c) + B1 + foot, RESTRICTED** | **67** | **0.484** [0.480, 0.487] | **3** (3, 13, 14) | 772.8 m | 500.4 m | 531.8 m | **dammed, 52 cells down** | 35 |

- **The raised cells return to the floor**: 86 against 86, and 67 against 66.
- **θ laid does NOT return to the no-clause CI.** 0.485 [0.482, 0.488] is outside [0.491, 0.496], and 0.484
  is outside [0.489, 0.495].
- **Every canyon comes back.**
  - B2 → A_c has the same four bodies, and the cols of 7 and 14 and 13's dam cell stand at the no-clause
    heights to 0.1 m (787.5 / 493.7 / 546.2 m).
  - B2 + B1 + foot has Finding 127's three.
  - The lakes return to the no-clause count (33; 35).
  - The coast, σ and both R8 values stay where the other B2 worlds are (`a2_worlds.txt`).
- **The pre-written table has no row for this case** (raised at the floor, θ outside the CI, canyons > 0).
  Rows 1 and 2 need θ inside the CI; rows 3 and 4 need the raised cells above the floor. It is stated as such,
  and no row is forced.
- **What the facts say, beside the table.** The restriction removes the clause's takeovers upstream of the
  junction. A0 finds 99.9 % of the deformation there, and the canyons return with them. So **the dams and the
  deformation are the SAME cells**: the tributary's last cells, laid from the larger line's first sample above
  the junction.
  - This is the opposite of row 2's reading ("barrages et déformation ne sont pas les mêmes cellules"). It
    comes closest to row 4's verdict: the restriction does nothing useful, and A0 already locates the cells.
  - The θ that remains off (0.485) comes from the takeovers the restriction keeps, at or downstream of the
    junction. Which cells carry it is NOT measured.
- **⇒ Neither form of the clause passes.** The full one removes the canyons by raising the tributaries' last
  cells above their law. The restricted one keeps the law's cells and keeps the canyons. **The clause is
  kept gated and off; activation remains the author's.**
- Named, not built: a remedy must lay those tributary cells lower than the larger line's upstream sample
  without damming the trunk. For example, from the junction's own floor rather than the nearest covering
  sample.

## E — the one-sided construction, measured only (`e_below_the_law.txt`)

The population is the témoin's trunk ≥ 10 km² cells the construction leaves uncarved, with the field
below its law: **22 775** (Finding 129's 32.1 %). The drainage is that of the as-built témoin (stage ii).

| | ALL trunk ≥ 10 km² (base rate) | **below their law** |
|---|---|---|
| cells | 70 904 | 22 775 |
| **under a lake after drainage** | 6 258 (8.8 %) | **6 258 (27.5 %)** |
| ≤ 100 m above the sea | 32.9 % | 46.8 % |
| low plain (≤ 100 m and slope ≤ 0.01) | 11.3 % | 13.5 % |
| coastal (≤ 2 km of the sea) | 10.9 % | 17.7 % |
| D8 distance to the base, p10 / p50 / p90 | 2.3 / 16.9 / 42.2 km | 1.3 / **10.1** / 27.4 km |
| altitude, p10 / p50 / p90 | 31 / 170 / 484 m | 20 / 108 / 374 m |
| drained area, p10 / p50 / p90 | 12 / 29 / 247 km² | 13 / 38 / 343 km² |

- **Every trunk cell under a lake is one of them.** They make up 27.5 % of the population.
- **The other 72.5 % are not lakes, and mostly not plains** (13.5 %). They are the lower, larger reaches:
  10 km from the base at p50 against 17 km, and 108 m high against 170 m.
- **The law passes 128 m above the terrain there** (p50; p10 21 m, p90 343 m).
- Read, not treated: the construction's floor law (k·χ from the sea, k a PROXY) climbs faster than these
  lower trunks do. A remblai (F79-C's deposition class) would fill mostly non-lake reaches by a hundred
  metres. That is the datum for the author's later decision, and nothing was built.

## The clause table, updated

| unwanted effect | Finding 129's status | what Finding 130 measured | **Finding 130's status** |
|---|---|---|---|
| **B1 spurs, Finding 38** | mur ↔ mer PROVED, ON | ON everywhere; base rate held in every world; F38 holds | **PROVED, kept ON** |
| **B1 teeth** | bruit ↔ pied ON; the profile's +38 % named | ON in the B1 worlds | unchanged |
| **B2 → A_c canyons: the full clause** | deforms the law at construction | A0: 1 252 raised cells, 99.9 % upstream of the junction, at 0.10 W, laid +4.3 m above the junction's floor | **the deformation is located: the tributary's last cells** |
| **the restricted clause (A2)** | not run | raised at the floor (86); θ 0.485 outside [0.491, 0.496]; canyons 4 / 3 (all back) | **does nothing useful: dams and deformation are the same cells.** Gated and off |
| **B2 terrain R8** × 2.3 (LTD) | NOT MEASURED; B1 failed | B0 fails on MY outlet; B1 fails; the paper's own D8 tree (625 / 625) gives D8 errors 0.452 / 0.951 against 0.566 / 2.043 | **the reference is not reproducible from what [O03] publishes.** C not run |
| **the skeleton's flats** (new) | 6.05 %, "cardinal fallback" | 100 % on the pre-incision lakes' surfaces; the fallback pointer cardinal on 20.1 % only; P1 → 0 fallbacks but 97.3 % cardinal on them | **fallbacks solved, lattice worsened.** Named, not treated |
| **the construction's uncarved third** | named | 27.5 % under a lake (base rate 8.8 %); the rest lower, larger reaches, 10 km from the base; the law 128 m above the terrain (p50) | **measured, nothing built** |

## Predictions scored

**Mine** (`f130_predictions.md`):
- **P0 HELD.**
- **A**
  - **P-A0 HALF**: within 1 W HELD (99.8 %); the laying sample upstream and higher HELD (99.9 %);
    "35–65 % upstream" REFUTED (99.9 %); "≥ 15 % without a junction" REFUTED (0.3 %).
  - **P-A1 HALF**: 86 and 1 295 HELD; the CIs not overlapping HELD; their half-width ±0.004–0.010 REFUTED
    (±0.002).
  - **P-A2 HALF**: θ 0.475–0.488, outside the CI, HELD (0.485); 300–800 raised REFUTED (86); canyons 0 / 0
    REFUTED (4 / 3); cols within ±30 m of the full clause REFUTED (at the no-clause heights); "row 3" REFUTED
    (no row applies).
- **B**
  - **P-B0 HELD** (the text does not define it).
  - **P-B1 HALF**: D8 not reproduced HELD; "MAE in, RMSE out" REFUTED (both out).
  - **P-Bdiag REFUTED**: the errors do not move with the outlet.
- **P**
  - **P-P0 HALF**: on the lakes HELD (100 %); coastal < 15 % HELD (0); FlatPerturbation applied HELD; "the
    cardinal fallback" REFUTED (20.1 %).
  - **P-P1 HALF**: < 0.1 % HELD (0); the facet unchanged on 100 % REFUTED narrowly (99.971 %); LTD
    changing 1–5 % REFUTED (0.041 %).
- **C: NOT RUN**, as predicted.
- **P-E HALF**: < 30 % under a lake HELD (27.5 %); "≥ 50 % low (< 100 m)" REFUTED narrowly (46.8 %); the
  downstream side HELD in sense (p50 10 km against the base rate's 17 km).
- **Meta ("at least two wrong") HELD.**

**The reviewer's:**
- **A0** ("≥ 85 % upstream, within a W"): HELD (99.9 %, 99.8 %).
- **A2**:
  - "raised at the floor": HELD;
  - "θ laid within the CI of 0.493": REFUTED;
  - "canyons 0 on B2 + B1 + foot": REFUTED (3);
  - "one return on B2 → A_c": REFUTED (four);
  - "row 2": REFUTED.
- **B** ("the exact geometry finds the paper's D8; LTD passes in absolute"): REFUTED on D8; LTD not judged.
- **P0**: "mostly lake surfaces" HELD; "FlatPerturbation not applied to the skeleton's field" REFUTED.
  **P1** ("under 0.1 %"): HELD.
- **C**: not run.
- **E** ("> 60 % under a lake"): REFUTED (27.5 %).
- **Meta HELD.**

## Limitations, stated

1. **B's geometry is read from figures**:
   - the grid and D8's pointers exactly (vector coordinates);
   - q only within the interval D8's pattern allows;
   - the paper's theoretical-area computation not at all.

   Since the paper's D8 tree reproduces exactly while its errors do not, B cannot pass on what [O03]
   publishes. The next reference would need the authors' data, or a benchmark that publishes its areas.
2. **B0's failure was MY outlet.** The one-inflow diagnostic shows it does not explain B1's failure.
3. **The θ intervals are i.i.d. link bootstraps**, a lower bound on the uncertainty (spatially correlated
   links).
4. **P1's increment** (a tenth of each flat's rise, the flat_grad normalised per flat) is my adaptation of
   `resolve_flats` into the elevation. Garbrecht & Martz 1997 was not read (no PDF in `docs/refs`); the
   repository's method was used, as the round allows.
5. **A0's "junction"** is the chain of `line_parent` (each line ends on its parent's cell). The 4 cells without
   a junction belong to lines whose chain never meets the laying line.
6. **E's classes** (lake, low plain ≤ 100 m and ≤ 0.01, coastal ≤ 2 km) are mine, declared before the run.

## State

**Uncommitted**, awaiting the next round's go-ahead.
- `valley_construction.rs`: `trunk_band_downstream`, `Skeleton::line_parent`, `densify_index`, `carve_diag` /
  `CarveDiag`, `ltd_directions_masked`, and the new permanent test.
- `flow.rs`: `flat_resolution`.
- `cached_product.rs`: the key test's list.
- `f126_coast.rs`: Finding 130's benches.
- This report, and ADR Finding 130.

**Checks after the last change:**
- The guard reads 6 / 6 "= banc" (`guard_after_this_round.txt`).
- The definition's fresh build hashes to `a8d2d538d692c2f0`, the reference (`identity.txt`).
- `f129_b1` and `f129_b2` reproduce Finding 129 to the digit (0.368 / 0.440 / 0.403; 6 124 / 684 893 /
  1 407 / 77): `ltd_directions_masked`'s `f32` path is bit-identical.
- `cargo test -p ymir-core --release --lib`: 587 passed, 0 failed.
- `cargo check --workspace --tests --release`: clean.
