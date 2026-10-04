# Finding 145 — the base's first defect measured: `carve_diag`'s cones lay 18–23 % of the carved cells from a line of another basin, almost all deepened (the ridges replaced by the walls' intersection, divides lowered p50 ~420 m); the artificial captures are small and not separated from the measure's noise; C1 (the own-basin carve) costs nothing, keeps θ, R8 and canyons, but raises the walls across the divides ×1.5–3.9

**Status: measured, a candidate in a bench. No production code. Nothing promoted.**

**Raw outputs, in this folder:**
- `f145_grep.txt` (Part G);
- `f145_m_raw.txt` (M1–M3, C);
- `f145_cap_raw.txt` (M3's refined captures, an amendment);
- `checks_before.txt`, `checks_after.txt`.

**Also:** `etat_au_F144.md`, the one-page state. The predictions were written before any measurement
(`f145_predictions.md`, with the non-blind items declared). The instruments and the basin label were declared before
(`f145_declared.md`, with one non-blind amendment).

## Partie 0

1. **F144 committed** (`9c4deb5`).
2. **ADR corrected in place, flagged "[CORRECTED in F145, after F144.]"**:
   - F143's "H_f is probably computed with 10·S_loi" is false. H_f already used the real slope; G-tag's defect was
     its two-cell span.
   - **Method rule, recorded**: a reviewer's hypothesis enters the ADR labelled **"hypothesis to check"**, never as a
     fact (two false lines in four rounds: F140, F143).
3. **Recorded (ADR F145)**:
   - the author's decisions (2026-10-04): **the gorge is paused** (committed, gated, toggle hidden); the base is
     treated, starting with the cones;
   - not contested: the scope gives 18 bodies; the cap rises to ~240 m on small outlets;
   - the cost principle (249.8 s reference);
   - the paused gorge's state.
4. **`etat_au_F144.md`** is written: the real state, the rules earned since F131, the author's decisions in force,
   dated.
   - **The "v2 note" is not in the repository** (searched in `docs/`, `docs/handoff`, `docs/design`, the reports and
     the handoff archive). The page stands alone, meant to replace that note's Part 0; the rest of the note is not
     rewritten.
5. **Checks**: before = F144's after-checks (`9c4deb5`; no file under `src` changed since); after = `checks_after.txt`.

## G — grep (rules 11, 11c; `f145_grep.txt`)

- **NOTHING FOUND**: «cône», «échantillon le plus proche», «minimum entre deux lignes», «ligne de partage», «priorité
  par aire», and the tags «F135-C», «F141-Z», «F144-Z». The ADR writes them as sections ("Finding 135 … the nearest
  sample lays 99 %…", F141's and F144's "Z —").
- Earliest hits read:
  - **«nearest sample»**: Finding 127-C, «`carve` lays a cell from its NEAREST sample…»;
  - **«cross-line minimum»**: Finding 134;
  - **«largest drained area»**: Finding 128-A, the confluence clause (off in the témoin);
  - **«divide»**: Finding 37c (the spillways' escape);
  - **«watershed»**: Finding 73;
  - **«basin id»**: Finding 63;
  - **`carve_diag`**: Finding 130's gated diagnostics.
- In the code: `valley_construction.rs:139` (the nearest-sample doc), `flow.rs:1155` («Determine basin ID»).

## Q — the code, before any measurement

1. **How `carve_diag` assigns a cell to a sample (`who`)** (`valley_construction.rs`):
   - every sample seeds the cell it stands in if its cone lies below the terrain (2191–2201; condition at 2196);
   - then a Dijkstra-like front propagates each sample's Euclidean distance to the 8 neighbours. A cell keeps the
     NEAREST sample whose cone lies below its terrain (2202–2223; condition at 2216);
   - the value is that sample's cone, `zf + tan(28°)·max(0, d − W/2)` (`geo`, 2167–2187), lowered to the minimum over
     the 8 neighbours' samples **of other lines** (2306–2322);
   - and the terrain is only ever lowered (2348–2355).
   - The confluence band (2227) is off in the témoin.
2. **What bounds a cone's reach**: **no radius and no basin.** A front stops only where the cone rises above the
   terrain (2196, 2216). W only sets the flat floor's half-width. Two valleys' 28° walls therefore meet wherever both
   lie below the terrain, whatever the divide.
3. **A basin label at this stage**: yes, and free. `compute_flow` on the breached input (`valley_construction.rs:620`)
   returns `basins` (`flow.rs:136`, from `compute_basins` at 237 / 1122). The skeleton keeps `direction` and the
   areas but drops `basins`.
   - The bench labels a cell by the last land cell on its D8 path over the skeleton's directions (declared): the same
     outlet basin, cost O(n).

## M — the defect measured (the construction; the témoin = C2 /10 col; ON extended; B2 → A_c for M2)

**The control**: the test-only copy of `carve_diag` (this configuration), with its filter off, equals production bit for
bit in all three worlds (0 cells differ).

### M1 — the foreign-laid cells, against the own-basin carve

| world | carved cells | laid by a foreign-basin line | **real defect (Δz < −1 m)** | legitimate (≤ 2 cells from a divide, \|Δz\| ≤ 1 m) | Δz p5 / p50 / p95 (m) |
|---|---|---|---|---|---|
| témoin | 2 697 596 | 494 061 (18.3 %) | **490 220 = 18.17 %** | 1 162 (0.2 % of the foreign-laid) | −2 316 / −568 / −25 |
| ON extended | 2 164 855 | 502 253 (23.2 %) | **498 123 = 23.01 %** | 1 190 (0.2 %) | −2 309 / −559 / −25 |
| B2 → A_c | 6 485 290 | 661 789 (10.2 %) | **608 317 = 9.38 %** | 22 532 (3.4 %) | −850 / −159 / 0 |

By the distance to a divide (témoin):

| band | foreign | deepened | \|Δz\| ≤ 1 m |
|---|---|---|---|
| 0–2 cells | 315 628 | 314 088 | 1 162 |
| 3–10 | 122 150 | 120 823 | 732 |
| > 10 | 56 283 | 55 309 | 667 |

- **Almost every foreign-laid cell is deepened, by hundreds of metres.** 64 % lie within 2 cells of a divide.
- **What it is, by the code (Q2)**: the nearest valley's 28° wall lays a cell even across the D8 divide. Where the
  own valley's wall would lie higher, or not reach the cell at all (then the terrain stays), the cell is laid lower.
  **The ridge between two valleys becomes the intersection of their walls.**
- **The "legitimate" category, defined as |Δz| ≤ 1 m, is nearly empty (0.2 %).** The two cones never agree to a
  metre at a divide.

### M2 — F127's dams

In B2 → A_c, **none of F127's four true cols is a defect cell.** All four are laid by a line of their own basin, with
Δz = 0:
- (3852, 2337) at 492.4 m;
- (4551, 3280) at 623.5 m;
- (2253, 4495) at 538.3 m;
- (1805, 4580) at 119.5 m.

**F127's dams are not cross-basin cones.** At the construction these cols are not raised: F127 measured the dams on the
full world. **F144-Z's share in G-drained, recalled**: 99.3 % / 94.0 % (×1.4, p .5 / p 0), by a line-catchment test,
not this basin test.

### M3 — visible consequences

- **Divides lowered** (divide cells laid > 1 m below what the own-basin carve gives):
  - témoin: **225 961, p50 425.4 m, p90 1 701.6 m**;
  - ON extended: 229 241, p50 416.1 m, p90 1 695.4 m;
  - B2 → A_c: 584 618, p50 147.7 m, p90 655.3 m.
- **Artificial captures** (the route changes under the témoin's carve, and the own-basin carve keeps it):
  - **the raw count** (2.4–2.7 M cells, 5 700–6 600 km², thousands of pairs) **is an upper bound, not captures.** The
    route label separates 37 099 input basins, and only 184 of them reach 10 km² (holding 75 % of the land). A carved
    mouth shifted by a few cells moves a cell into a neighbouring coastal mini-basin;
  - **refined** (`f145_cap`, the amendment: both basins ≥ 10 km²):
    - témoin: 171.5 km² in 30 pairs, against a noise floor of 141.8 km² (the own-basin carve's own moves);
    - ON extended: 197.8 km² in 37 pairs, noise floor 161.1 km²;
  - **the largest refined pairs are still adjacent mouths**: 78.7 km² from (3282, 2955) to (3281, 2956), 20.1 km² and
    12.4 km², all one cell apart;
  - **real captures between distant outlets exist but are small**: 9.2 km² from the 255 km² basin of (2651, 5908) into
    the 1 805 km² basin of (2758, 5210); 7.0 km² from the 1 373 km² basin of (3288, 2697); 6.3 km² from the 449 km²
    basin of (4058, 4319); in all about 70 km², ~0.3 % of the land;
  - **the measure does not separate them from its noise floor**, which is not decomposed.
- **The stop rule does not trigger**: a real defect of 9–23 % ≫ 1 %, and captures exist.

## C — C1, the own-basin carve, against the témoin's carve (rule 18; at the construction)

**C1**: a sample lays a cell only if its line's basin is the cell's basin. A cell no own-basin sample reaches **keeps
its terrain**: the restriction must not invent a floor.

| world | carve | **walls across the divides** (> 28°) | θ | R8 terrain | R8 network | canyons | relief p50 (m) | coast (spurs /km) | carve time |
|---|---|---|---|---|---|---|---|---|---|
| témoin | témoin's (negative control) | 182 091 | 0.493 [0.492, 0.495] | 0.0579 | 0.6904 | 0 | 565.9 | 0.0096 | 2.2 s |
| | **C1** | **271 169 (×1.49)** | 0.495 [0.493, 0.496] | 0.0571 | 0.6930 | 0 | 577.0 | 0.0096 | 2.2 s |
| ON extended | témoin's | 182 568 | 0.439 [0.431, 0.446] | 0.0508 | 0.6637 | 0 | 592.7 | 0.0096 | 1.6 s |
| | **C1** | **271 534 (×1.49)** | 0.440 [0.432, 0.448] | 0.0498 | 0.6650 | 0 | 605.4 | 0.0096 | 1.6 s |
| B2 → A_c | témoin's | 269 364 | 0.493 [0.491, 0.496] | 0.1447 | 0.6664 | 0 | 430.2 | 0.0067 | 11.4 s |
| | **C1** | **1 037 566 (×3.85)** | 0.495 [0.492, 0.497] | **0.1698** | 0.6687 | 0 | 443.6 | 0.0067 | 9.2 s |

- **Canyons stay 0.** θ stays within the declared ±0.01 (+0.001 / +0.002). R8 network moves by < 0.003. The coast is
  unchanged.
- **The relief rises +2–3 %**: the divides are kept higher.
- **F127's four cols are unchanged** (all own-basin).
- **Walls across the divides rise ×1.49 (témoin, ON) and ×3.85 (B2 → A_c)**: the restriction leaves each side of a
  divide at its own wall, or at the terrain, so the divide becomes a step. R8 terrain rises at B2 → A_c (0.145 →
  0.170).
- **The cost: none.** The filtered carve takes the same time (2.2 s; 1.6 s; 9.2 against 11.4 s). In production the
  basin labels already come out of `compute_flow`. Against run_hd 249.8 s: ±0.
- **Captures under C1** (the route changed against the input, raw): 717 905 / 743 627 / 1 308 435 cells. That is the
  raw measure's noise, the same mouth-shift artefact.

## Predictions

**The reviewer's (hypotheses to check):**
- **M1**:
  - "the real defect 3–10 % of the carved cells": **refuted** (18.2 % témoin, 23.0 % ON; 9.4 % only in B2 → A_c);
  - "> 60 % of the foreign-laid within 2 cells of a divide with |Δz| ≤ 1 m": **refuted** (0.2 %). 64 % lie within
    2 cells, but almost all deepened.
- **M2** (≥ 3 of F127's 4 dams are real-defect cells): **refuted** (0).
- **M3** (≥ 1 artificial capture on the témoin): **held in kind** (distant-outlet captures of 6–9 km²), but not
  separated from the measure's noise.
- **C1**:
  - canyons stay 0: **held**;
  - F127's dams decrease: **refuted** (unchanged);
  - walls across the divides > ×2: **refuted** on the témoin and ON (×1.49), **held** on B2 → A_c (×3.85);
  - cost < +5 s: **held** (±0).
- **Meta**: **held**.

**Mine** (`f145_predictions.md`):
- **P-M1**:
  - "foreign-laid 8–20 %": **held** (18.3; 10.2), except ON (23.2);
  - "most deepened": **held**;
  - "real defect 5–15 %": **refuted** on the témoin and ON (18.2 / 23.0), **held** on B2 (9.4);
  - "legit < 40 %": **held** (0.2 %).
- **P-M2** (≤ 1 of 4): **held** (0).
- **P-M3**:
  - "divides lowered, p50 > 5 m": **held** (p50 ~420 m);
  - "≥ 1 capture": **held** in kind;
  - "the stop rule does not trigger": **held**.
- **P-C1**:
  - canyons 0: **held**;
  - dams unchanged: **held**;
  - walls > ×2: **refuted** on the témoin and ON, held on B2;
  - captures fall to ~0: **refuted** (the raw measure's noise stays);
  - θ ±0.01: **held**;
  - cost < 5 s: **held**.
- **Meta**: **held**.

## What the round says (no decision taken; nothing promoted)

1. **The cones' cross-basin laying is large and systematic**: 18–23 % of the carved cells in the témoin and ON, almost
   all deepened by hundreds of metres. The divides fall to the intersection of two 28° walls. That is the construction
   as written (Q2: nothing bounds a cone but the terrain), not a rare edge case.
2. **The captures it causes are small**: ~70 km² between distant outlets in the témoin, below the measure's own noise.
   So the defect changes the ridges, not (much) the drainage.
3. **F127's dams are not cross-basin.**
4. **C1 is free and neutral on θ, R8 network, canyons and the coast.** It trades lowered divides for steps at the
   divides (walls ×1.5, ×3.9 at A_c), because the restriction keeps the terrain where no own wall reaches.
   **A hard basin boundary is a seam**: F143 / F144's pattern again.
5. **Candidates for the author** (none built):
   - (i) C1 with a divide transition: own-basin first, the foreign wall allowed within a few cells of the divide and
     blended;
   - (ii) a cone limited to its basin AND its own wall's reach, so that a divide sits where the two own walls meet;
   - (iii) keep the cones as they are and accept the walls' intersection as the ridge (the témoin's look), now that
     the captures are measured small.

## Limitations, stated

1. **The construction stage only.** No production code, so C1 was not run through the light pass, the droplets or the
   breach. Canyons, the relief and R8 are measured on the construction field.
2. **One world, one seed.**
3. **The capture measure** needs an outlet merge (adjacent mouths). Its noise floor is not decomposed.
4. **The basin label** is the input's D8 outlet basin. A line's basin is the majority over its samples.
5. **M2's cols** are at the construction; F127's dams are a full-world effect.

## State

**Uncommitted**, awaiting the go-ahead:
- `f126_coast.rs`: `f145_m`, `f145_cap`;
- ADR: F143's line corrected in place, Finding 145;
- this folder (with `etat_au_F144.md`).

**No production code.**
