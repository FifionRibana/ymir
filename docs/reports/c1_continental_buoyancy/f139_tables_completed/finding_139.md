# Finding 139 — the tables completed: T3 is calibrated on the D8 lakes; "steepen" removes every shortage fall (m_fit ≤ 51, never at the cap); a per-lake retreat spreads the emptying over half the selector; the merged spill recovers the below-sea lakes' level

**Raw output, in this folder:** `f139_t.txt` (bench `f139_t`), `checks_before.txt`, `checks_after.txt`.

The predictions were written before the measurement (`f139_predictions.md`, C1 declared not blind). The instruments
were declared before it (`f139_declared.md`). **The specification stays v2; its decisions stay open.** No production
code.

**World**: the témoin, ON extended, run_hd tail. F138's instruments, with their amendments.

## Partie 0

- **The "FNV" label** in F138's raw output (`f138_t.txt`) is corrected: the code draws φ with splitmix64. The
  correction is marked in the file, in `finding_138.md` and in the bench's source.
- **F138 committed** (`e4802c1`). Added to ADR F138, the review of v2:
  1. one r per world empties every lake at once;
  2. "fixed slope, the rest falls" makes falls at large outlets;
  3. the rim clamp on `carve_diag`'s output risks a one-cell ridge;
  4. at ×1, r = 1 reduces the lakes by only 14 %;
  5. three columns were missing.
- **Checks**: before = F138's after-checks (no production code changed since, declared). After = `checks_after.txt`.

## C1 — T3's calibration (T3's area at r = 1 against the final area, run_hd mask)

- **The 14 D8 lakes: ratio 0.982–1.000.** T3's hypsometry on the input reproduces the final lakes within 2 %.
- **The spillway lakes**:
  - 4 of 8 within 11 %: 1000012 0.999, 1000013 0.971, 1000014 0.949, 1000018 0.898;
  - **4 off by > 20 %**: 1000011 0.409, 1000016 0.265, 1000017 0.120, 1000019 0.568. Their final level is the
    merged one, above their body's ring (C5).
  - 1000001 has no body.
- **4 of 22 off by > 20 %**: under a third, so no stop.

## C2 — T1b, the minimal steepness m_fit (23 descents)

`S_req` = the least slope from L that meets the law below before the next base; `m_fit = S_req / S_loi(A)`.
- **m_fit: min 0.5, p50 6.6, max 50.7.** **No descent needs more than the 28° cap**: the largest S_req is 6.57°
  (1000018: 203 m over 1.76 km).
- The highest m_fit: 1000014 50.7, 1000018 42.5, lake 10 28.3, 1000011 27.8, 1000016 26.8, 1000019 20.2, lake 1
  18.8.
- The D8 lakes' S_req is 0.29–4.08°; the steepest is lake 10, 898 m over 12.6 km.

## C3 — the rule "steepen" (m = 10)

| rule | shortage falls ≥ 1 m | median | max |
|---|---|---|---|
| fixed slope (F138) | 10 | 129 m | 581 m |
| **steepen** (S = max(10·S_loi, S_req), capped) | **0** | — | — |

The steepen rule printed 3 "falls" of 0 m. They are float ties: S_req sits exactly on the fit test's limit. All
three are < 0.5 m.

**Under "steepen", no shortage fall exists in this world.** The steepened gorges are 1.1–6.6° (lake 10 4.1°,
1000018 6.6°).

## C4 — T4's variant (m = 10; the head fall φ·D_g₀ first, the gorge from L − φ·D_g₀)

| rule | fit / no room | shortage falls |
|---|---|---|
| fixed slope | **16 / 7** (F138: 13 / 10) | 7, median 138 m, max 581 m |
| steepen | 21 / 2 (both float ties, 0 m) | 0 |

- With the head fall taken first, **lakes 2, 3 and 9 now fit** under the fixed slope (φ 0.25 / 0.25 / 0.34, head
  falls 199 / 181 / 101 m).
- Lakes 1 and 10, and 1000016 and 1000019, have φ = 0 and keep their shortage.

## R — retreat per lake: `r_lake = min(2, r_world·(A/A_ref)^p)`, 14 D8 lakes, A_ref = 418.7 km²

| p | r_world 0 | 0.5 | 1 | 1.5 | 2 |
|---|---|---|---|---|---|
| 0 (v2) | 14 lakes, 1 308.5 km² | 14, 1 223.6 | 14, 1 145.7 | 14, 367.3 | **0**, 0 |
| 0.25 | 14, 1 308.5 | 14, 1 214.3 | 14, 845.5 | 12, 135.4 | 3, 5.2 |
| 0.5 | 14, 1 308.5 | 14, 1 202.7 | 14, **521.3** | 9, 103.0 | **4**, 19.0 |

**The emptying age** (map (a): r 0 → ×0.7, 1 → ×1, 2 → ×1.4, piecewise linear):
- **p = 0**: all 14 at ×1.40 (spread 0 %). This is v2's defect, reproduced.
- **p = 0.25**: 8 lakes from ×1.19 to ×1.40 (spread **30 %**); 6 never within the range.
- **p = 0.5**: 8 lakes from **×1.03 to ×1.40** (spread **53 %**); 6 never within the range.
  - The largest outlets empty first: 19 at ×1.03, 1 at ×1.05, 13 ×1.12, 11 ×1.14, 10 ×1.17, 3 ×1.22.
  - Never within the range: 2, 4, 7, 8, 9, 20.

At ×1 (r_world = 1) no lake is empty at any p. At p = 0.5 the total area there is already down to 521 km² (−55 %
against r = 0), against −12 % at p = 0.

**p's provenance**: a DECISION. The link to the knickpoint celerity is ANCHORED in form (F115: c = K·A^m), not in
value.

## C5 — the below-sea lakes: L_in = the merged spill

`L_in_merged` = the median, over the lake's final footprint, of `ocean_flood`'s spill on the breached input (F38's
"one water body to its col").
- **8 of 9 within 1 m of today's level**:
  - 1000011 459.3 m (against the body ring's 325.8);
  - 1000012, 1000013 and 1000014 exact;
  - 1000016 115.3 (ring −1.1);
  - 1000017 268.8 (ring 126.1);
  - 1000018 −0.6 m;
  - 1000019 614.9 (ring 570.2).
- **The definition holds** for these 8. L(1) = min(L_in_merged, today) = today.
- **1000001 fails**: merged spill 3 046.6 m against a level of 176.8 m. It is F135-K4's anomaly, already named:
  a 1.2 km² below-sea-merge lake in a hole the construction dug into a ~3 000 m high (ring cut 2 968 m in OFF
  and ON). It stays excluded, named, not explained.

## Predictions

**Mine:**
- **P-C1** held (18 of 22 within 20 %; C1 was declared not blind).
- **P-C2** held: p50 6.6 in (3, 10]; max 50.7 in 30–60; none beyond the cap. The "min ~1" came out 0.5.
- **P-C3** held: 0 falls.
- **P-C4** held: 3 more fit, lakes 2, 3, 9.
- **P-R** held: 53 % spread; 0 lakes empty at r_world = 1.
- **P-C5** held: 8 of 9.
- **Meta** ("at least one of all these is wrong"): **refuted**, none is.

**The reviewer's:**
- **C1** held (18 of 22).
- **C2** held (p50 6.6; none beyond the cap).
- **C3** held (≤ 1: 0).
- **C4** held (≥ 2: 3).
- **R**: "the emptying spreads over more than half of the range" held (53 %). "At r_world = 1, ≥ 3 lakes empty"
  **refuted** (0).
- **C5** held (8 of 9).
- **Meta** held.

## What the tables say for the open decisions (no decision taken)

- **"Steepen" makes the shortage fall disappear in this world.** With it, falls come only from φ (the head),
  which is the author's « les deux » reduced to one cause here: no room is ever lacking below 6.6°.
  - Under "fixed slope" the shortage falls are the large outlets' (lake 10 581 m), which review point 2 called
    unphysical.
- **A per-lake retreat (p > 0)** spreads the emptying (53 % of the selector at p = 0.5). It keeps 4 lakes at the
  oldest age, and shrinks the lakes at ×1 by −55 % instead of −12 %.
- **The below-sea lakes** get L_in from the merged spill (8 of 9). 1000001 stays out.

## Limitations, stated

1. **m_fit and steepen** use today's levels (r = 1) and the law along the ON path with the world's accumulation.
2. **The "0 m falls"** are float ties at S_req, counted as none.
3. **R's areas** are T3's hypsometry. It is calibrated on the D8 lakes (C1), not under retreat.
4. **The emptying age** uses map (a); map (b) shifts it, and is not tabulated.
5. **1000001** stays an unexplained construction artefact.

## State

**Uncommitted**, awaiting the go-ahead:
- `f126_coast.rs` (`f139_t`; F138's label correction is already committed with F138);
- this folder;
- ADR Finding 139.

**No production code. The specification stays v2, its decisions open.**
