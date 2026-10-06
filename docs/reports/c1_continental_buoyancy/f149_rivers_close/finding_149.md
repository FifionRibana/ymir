# Finding 149 — rivers, extension 2 of 2 (the last round): the masking toggle, the K fix, one attempt at the head score (the stop rule fires again), and the rivers closed

**Status: T and K built; H measured and NOT built (the stop rule fired). The bundles and micro-rivers become an
accepted defect, masked by T. Nothing committed.**

**Files:**
- predictions: `f149_predictions.md`, written before any measurement;
- instruments: `f149_declared.md`;
- raw output: `f149_rivers_raw.txt` (T and H) and `f149_lake_guard_raw.txt` (K, the guard regenerated);
- the lake reference before the regeneration: `bench_lake_hashes_BEFORE.json`;
- checks: `checks_before.txt` and `checks_after.txt`.

## Part 0

- F148 committed (`d937aec`, not pushed).
- `docs/adr/accepted_defects.md` § 2: the smoothed polylines clip off-valley cell corners (134 stretches < 0.1 cell,
  122 m over 21 541 km).
- **Rule 15, recorded**: F148's stop rule is not relaxed after the fact. Part H used it unchanged.
- Checks before = F148's after-checks.
- **Checks after** (`checks_after.txt`, run 2):
  - `cargo check --workspace` clean;
  - lib 610 passed (608 + 2), 0 failed;
  - viz 30 passed;
  - **the guard 6 / 6, field and lakes = banc through `run_hd`** (C2 /10 col's field `a8d2d538d692c2f0`);
  - `rivers_ll` takes 0.22–0.40 s per state.
- **Run 1 of the after-checks failed** on lakes for the four regenerated states (`checks_after_run1_STALE_CACHE.txt`).
  - The cause is my error: K changes the HD drainage bundle and `ALGO_HD_DRAINAGE` was not bumped, so `run_hd`
    served the pre-fix bundle from its cache.
  - Amended as declared (`f149_declared.md`): `ALGO_HD_DRAINAGE` 11 → 12.
  - The bump changes the HD drainage key, the lake guard's `digest`, so the reference was regenerated a second time.
  - Against the first regeneration, only the six `digest`s changed and every `lakes_hash` is identical
    (`bench_lake_hashes_REGEN1.json`, `f149_lake_guard_raw_regen2.txt`).

## T — the masking toggle (built; display only)

- **`micro_lake_inflow`**: shorter than 1 km and ending in a lake. It is F148-M (b)'s definition; the bench checks
  the flag equals it on every river.
- `rivers_ll.json` goes to **0.5.0** (11.32 MB); the doc is updated.
- The viz box « Masquer faisceaux et micro-rivières » hides `fusion_candidate || micro_lake_inflow`. The box
  « masquer les déversoirs retracés » stays separate. Each shows its count and km.

On the témoin (10 246 rivers, 21 541 km):

| box | rivers | km |
|---|---|---|
| fusion candidates | 942 | 3 430 |
| micro_lake_inflow | 726 | 351 |
| both (overlap) | 0 | — |
| **« Masquer faisceaux et micro-rivières »** | **1 668 (16.3 %)** | **3 781 (17.6 %)** |
| + the retracing spillways | 5 | 53 |
| all boxes | 1 673 | 3 834 (8 573 rivers, 17 707 km left) |

T touches only the export's flags and the viz. `rivers.json` is byte-identical (the existing test); the guard is
below.

## K — the lake 1000011 (the core fix, built)

`exorheic_lakes_missing_outlet`: a lake has an outlet if a segment carries `segment_source_lake = Some(id)`, in
addition to the 1-cell border test. No wider distance.
- A permanent test: a tagged spillway 2 cells off the shore keeps its lake Exorheic.
- Its negative controls: the same segment untagged, and a lake with no segment, both go Unresolved.
- The F86 tests stay green.

**The type changes, per guard state** (`f149_lake_guard_raw.txt`). They are exactly the Exorheic lakes whose only
outlet is a tagged segment, since the fix changes nothing else.

| state | lakes kept Exorheic by the fix (formerly Unresolved) | Unresolved left | `lakes_hash` |
|---|---|---|---|
| livré | none | 3 (not moved by the fix) | unchanged |
| A1+B2 | 1000016 (404 km², spillway 11803) | 0 | changed |
| C1 nue | 1000017 (325 km², spillway 15051) | 2 (not moved by the fix) | changed |
| **C2 /10 col (the témoin)** | **1000011 only** (494 km², spillway 16217) | **0** | changed |
| C2 /3 | none | 0 | unchanged |
| C2 /10 niveau mer | 1000017 (374 km², spillway 12793) | 0 | changed |
| (C2 /10 col, lake_base ON; not guarded) | 1000011 (493 km², spillway 13468) | 0 | not written |

**The regeneration, declared** (`bench_lake_hashes_BEFORE.json` → `crates/ymir-core/data/bench_lake_hashes.json`):
- the diff touches the four states with a non-empty list, and only their `lakes_hash`;
- only the hash's second half (the lakes' attributes) changes; the first half (the footprints) and every lake count
  are unchanged;
- the six fields reproduce (C2 /10 col `a8d2d538d692c2f0`).
- **Gate: held.** On the témoin only lake 1000011 changes type; one lake changes in each of three other states, by
  the same mechanism (a below-sea-path basin whose tagged spillway starts off its shore).

## H — one attempt: the head score (measured; the stop rule fires; NOT built)

**The score**: the mean of C (F148's convergence) over the first 10 trace cells. N = 10 is F148's head score, NOT
blind.

**The rule** judges a head once: a segment with no channel input is a channel from its start if its head score is ≥
θ, else not at all. Lake outlets and spillways are exempt; a channel input carries the channel downstream.

**The populations** (F148's definitions, same build): (a) 873 · (b) 705 · (a)∪(b) 1 578 · (c) 2 828.

**σ = 2 (declared)**: head AUC 0.784. (c)'s head quantiles p0.5 / p1 / p2 / p5 = −2.48 / −1.60 / −1.16 / −0.71 km⁻¹.

| θ | rivers | km | pairs | micro-lake | (a) | (b) | **(a)∪(b) removed** | **(c) lost** |
|---|---|---|---|---|---|---|---|---|
| (c) p0.5 = −2.48 | 10 070 | 21 058 | 1 625 | 715 | 1.6 % | 1.6 % | 1.6 % | 0.39 % |
| (c) p1 = −1.60 | 9 893 | 20 673 | 1 602 | 700 | 2.5 % | 3.7 % | **3.0 %** | **0.60 %** |
| (c) p2 = −1.16 | 9 703 | 20 268 | 1 574 | 686 | 3.9 % | 5.7 % | 4.7 % | 1.45 % |
| (c) p5 = −0.71 | 9 265 | 19 323 | 1 481 | 631 | 9.7 % | 13.2 % | 11.3 % | 3.29 % |
| (c) p10 = −0.33 | 8 326 | 16 971 | 1 111 | 541 | 28.9 % | 26.0 % | 27.6 % | 6.40 % |
| (c) p20 = −0.02 | 6 796 | 13 239 | 604 | 348 | 58.9 % | 53.2 % | 56.3 % | 12.94 % |
| 0 | 6 665 | 12 988 | 555 | 329 | 61.4 % | 55.9 % | 58.9 % | 13.44 % |

**σ = 1 (NOT blind, reported beside)**: head AUC 0.857.
- At ≤ 1 % of (c): at best 7.2 % of (a)∪(b) (θ = −0.94, 0.64 %).
- Past 50 %: θ ≈ 0, (c) loses 3.6–3.7 %.

**THE STOP RULE FIRES for σ = 2, and for σ = 1 too.** Nothing is built: no rule in `rivers_ll`, no head-rule toggle.
**The bundles and the micro-rivers become an accepted defect, masked by T** (`accepted_defects.md` § 3).

**The control rivers lost at the θ nearest the 1 % bound** (σ = 2, θ = p1, 17 rivers):
- length p10 / p50 / p90 / max 0.54 / **8.17** / 16.4 / 17.9 km; **8 of 17 ≥ 10 km**;
- 11 ≥ 50 % on a construction floor;
- **3 pair members**;
- heads: **13 on uncarved ground**, 3 on a construction floor, 1 on a carved wall;
- ends: 12 confluences, 4 lakes, 1 sea.

**They are real rivers, not bundle members the control let through.** Their HEAD lies on planar uncarved ground,
and their valley is built downstream. Losing them is what costs the rule its threshold.

At σ = 1, θ = p2: 18 lost, p50 5.8 km, 7 ≥ 10 km, 4 pair members, 10 heads uncarved.

**The reason for the defect**: the separation is real at the head (AUC 0.78; 0.86 at σ = 1, non-blind), but the
control's tail overlaps it.
- About 1–2 % of the control's rivers rise on ground as planar as the bundles'.
- A head rule cannot drop the bundles without dropping them.

## Predictions

**The reviewer's (hypotheses to check):**
- H, "the stop rule fires at σ = 2": HELD.
- H, "at best 15–30 % of (a)∪(b) for 1 % of (c)": REFUTED (3.0 % at 0.60 %; 4.7 % at 1.45 %).
- H, "at least a third of the lost control rivers are parallel traces on the built floors": REFUTED (3 of 17 pair
  members; 3 of 17 heads on a floor).
- K, "only lake 1000011 changes on the témoin": HELD.
- T, "~1 600 rivers": HELD (1 668).
- Meta: HELD.

**Mine:**
- P-H:
  - "the stop rule fires": HELD;
  - "5–20 % at ≤ 1 %": REFUTED (3.0 %);
  - "σ = 1 also fires, ≤ 35 %": HELD (7.2 %).
- P-H, the lost control:
  - "mostly short (< 3 km) on a floor": REFUTED (p50 8.2 km; 11 of 17 on a floor ≥ 50 % though);
  - "fewer than a third pair members": HELD (3 of 17);
  - "heads mostly on a floor or a wall": REFUTED (13 of 17 uncarved).
- P-K:
  - "only lake 1000011 on the témoin": HELD;
  - "0–2 lakes in each other state": HELD (0, 1, 1, 0, 1);
  - "the field guard untouched": HELD.
- P-T, "~1 650 rivers, ~3 800 km": HELD (1 668, 3 781).
- P-cost, ~0: HELD (a flag; a set lookup).
- Meta: HELD.

## F — the rivers closed: a summary for the author

**What `rivers_ll.json` 0.5.0 delivers** (a new container layer; `rivers.json` unchanged):
- **every river of the drainage as a main stem** (Hack's convention at every confluence, the T3b junctions included);
- **smoothed polylines in their valley floor**:
  - within 0.13 cell of the D8 trace (p99);
  - no crossings between rivers, no self-intersections;
  - a tributary ends exactly on a vertex of its receiver;
  - spillways that retrace a watercourse share its vertices.
- **per river**: length, catchment, Strahler, discharge, the end (sea, lake, endorheic lake, confluence, terminal),
  the lakes touched, the length rank in its basin, and the flags `parallel_of`, `fusion_candidate`,
  `retraces_spillway`, `micro_lake_inflow`;
- **per vertex**: catchment, accumulated discharge, slope, valley-floor width W(A), bed width a·Q^0.5 (a = 5 in the
  header, b = 0.5 anchored);
- ~0.4 s per world.

**What is accepted** (`docs/adr/accepted_defects.md`):
- § 2, the corner clips: 134 stretches under 0.1 cell;
- § 3, the bundles and the micro-rivers into lakes: not removable by a head rule without losing real rivers; flagged
  and masked by the toggle.

**What Living Landz does on its side**:
- **the hex grid**: compute the edges from the polylines and verify the radius reading (centre to corner) and the
  origin (hex (0,0) centred on the bottom-left corner) against `docs/rivers_ll_hex_reference.json`;
- **the ford / obstacle / bonus rules on hex edges**, from `vertex.bed_width_m` (≥ 69 m is wider than one hex on
  ~1 % of the river km at a = 5), the discharge and the slope;
- **the selection** (length, catchment, Strahler, the flags);
- **the masking**: hide `fusion_candidate || micro_lake_inflow` (and `retraces_spillway`) as the viz toggle does;
- **the SDF rendering**.

**What stays queued**:
- the breach's ramp anchored on each pit's floor (F141's γ);
- the Unresolved lakes the fix leaves (livré 3, C1 nue 2). The fix does not move them, so they are not this round's
  mechanism. Their reasons were not read this round.
- K is closed: lake 1000011's outlet is recognised.

**The next work is the geology v1.**
