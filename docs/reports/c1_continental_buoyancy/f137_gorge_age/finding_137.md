# Finding 137 — the age sets how far the gorge has retreated: the measurements before the design, and its specification (nothing built)

## The author's criterion (recorded as given)

- « Je veux réduire la taille des lacs avec l'âge, car au départ ils étaient tous remplis à la hauteur
  d'équilibre, ce qui donnait des lacs énormes. » (2026-09-30)
- « Un lac présent est un niveau de base pour les rivières qui s'y jettent ; un lac vidé ne l'est plus. »
  (2026-09-29) « Les lacs des bassins sous la mer sont des lacs présents. » (2026-09-30)
- « Chute ou gorge, les deux sont valides. » « Le cas des chutes lié au type de roche n'est pas encore
  activable […]. Ce sera donc intéressant de faire apparaître des chutes d'eau (et donc de les tagger comme tel)
  là où c'est pertinent. Ça donnera de la diversité au monde. » (2026-10-02)
- « C'est bien d'avoir l'âge qui fixe jusqu'où la gorge a reculé. » (2026-10-02)
  - The reviewer's reading, accepted by the author: a young continent has a full lake and a gorge downstream; an
    old continent has a notched rim and a lowered lake.
- **Cost**: a world in minutes, not hours.

**Raw outputs, in this folder:**
- `grep_rule11_age.txt`;
- the bench `f137_q.txt`;
- the map `q2_k_field_map.png`;
- `checks_before.txt` and `checks_after.txt`.

The predictions were written before any grep or measurement (`f137_predictions.md`). The instruments were
declared before the measurements (`f137_declared.md`). **The specification is `spec_gorge_age.md`.**

**Units**: 1 cell = 48.8 m (domain).

## Partie 0

- F136 committed (`c209c9c`). Added to ADR F136:
  - K5 is not blind (both readings kept);
  - G = F114's family 1;
  - M's split (59.9 / 33.7 %);
  - K7's one-sided invariant.
- **ADR:259 corrected in place**, with the correction marked: the C-3 erodibility field exists.
- Checks: before = F136's after-checks (no code changed since, declared). After = `checks_after.txt`. F137
  changes no production code (one bench).

## Q1 — F114's family 1, and the delivered worlds today

**Where it was measured.** F114's bench (`tests/f114_sill.rs:85`) builds `Knobs { slope_floor_abs: Some(0.024),
..passes(2) }`:
- **the delivered stream-power world with the age closure (A1+B2)**, before the construction existed (F121);
- at the framing of its time (Findings 107–123's origin);
- read on the eroded and breached fields.

**Today, with F136-K6b's amended instrument**, at the canonical framing, through run_hd's tail:

| world | final lakes | ≥ 1 km² with a D8 outlet (without) | steps > 50 / > 200 / > 500 m | drop over 5 km > 200 m |
|---|---|---|---|---|
| **livré** (delivered stream power) | 60 | 43 (11) | **9 / 6 / 3** | 19 lakes |
| **A1+B2** (F114's world) | 38 | 18 (11) | **5 / 2 / 1** | 6 lakes |
| **ON extended** (construction) | 34 | 14 (12) | **8 / 6 / 2** | 7 lakes |

- The step is **not an artefact of the construction alone**: the delivered world has as many steps > 200 m (6),
  and many outlet descents the instrument cannot close. Its "drop over 5 km" reaches 1 259, 1 373 and 1 148 m
  with no graded reach within 10 km.
- The steepest: livré lake 1 at 28.6°.
- **The OFF construction HID it**: OFF cuts through the cols (F135-K1: 11 cols > 100 m), so it has neither the
  lakes nor their steps. The delivered world and ON extended both have them.

## Q2 — C-3's `production_k_field`

**Its source** (`production_upscale.rs:662`, `lithology/mod.rs`, `fracture`):
- `build_coarse_k`: the coarse 64² grid (6.25 km cells). Hard continental basement ×1; **rift-soft ×10** where a
  continental cell's age < 1. Upscaled **bilinear** at the altitude's mapping.
- `stamp_volcanic_k`: `max(K, 3)` on the edifices' basal discs.
- Times **the C-3b fracture density factor** (`build_hd_density_k`).

**Its range on land** (11.3 M cells):
- K: p1 1.04, p10 1.06, **p50 1.27**, p90 5.46, p99 23.8, max 56.7.
- The fracture factor alone: 1.03–7.0, **≠ 1 everywhere**.
- Histogram: [1, 1.5) 61.8 %, [1.5, 2.5) 16.8 %, [2.5, 3.5) 5.0 %, [3.5, 6) 7.5 %, [6, 10) 3.3 %, > 10 5.7 %.
- Volcanic discs 4.29 % of land; rift-soft 4.87 %.

**Its spatial structure** (autocorrelation, first lag under 1/e):

| field | x | y |
|---|---|---|
| K | 12.5–25 km | > 100 km |
| log K | 25–50 km | > 100 km |
| lithology classes | 12.5–25 km | 25–50 km |
| fracture factor | 50–100 km | > 100 km |
| altitude (reference) | 12.5–25 km | 6.25–12.5 km |

- K's ACF at 64 cells (3.1 km) is still 0.94.
- The map (`q2_k_field_map.png`, north up) shows a smooth ramp (the fracture density, rising to the south-east),
  one coarse rift strip in the east (blocky 6 km cells), and the volcanic discs.

**Its correlations** (Pearson r of log K on land) with:

| variable | r |
|---|---|
| the fracture factor | **0.746** |
| the rift class | **0.619** |
| the volcanic mask | **0.451** |
| the slope | 0.330 |
| the altitude | 0.142 |
| the cratonic mask | −0.178 |
| the precipitation | **−0.006** |

It is tectonic and volcanic, weakly related to the relief, and unrelated to the climate.

**The stages that read it**:
- the stream-power incision (`incise_with_floor`, `stream_power.rs:751`): **the delivered world AND the
  construction's light pass** (`production_upscale.rs:509–553`);
- **`carve` does not read it** (the construction proper);
- the infiltration uses its own conductivity matrix, not this field.

**⇒ It is not a hardness map and cannot place local waterfalls alone.** Its only local structure is its
contrasts: the volcanic discs' edges and the rift strip's edges. The specification uses it only there, as an
optional and gated criterion (spec § 5).

## Q3 — the length and class of F136's steps (non-blind: F136's numbers)

| lake | declared H / L → mean slope | class | amended drop / distance → mean slope | class |
|---|---|---|---|---|
| 1 | 136 m / 0.33 km → 22.4° | steep | 290 m / 1.94 km → 8.5° | gorge |
| 2 | 413 / 0.61 → 34.1° | steep | 862 / 3.31 → 14.6° | steep |
| 4 | 119 / 0.39 → 17.0° | steep | 119 / 0.39 → 16.8° | steep |
| 7 | 158 / 0.30 → 27.7° | steep | 317 / 1.85 → 9.7° | gorge |
| 8 | 58 / 0.22 → 14.8° | steep | 51 / 0.12 → 23.5° | steep |
| 9 | 281 / 1.55 → 10.3° | steep | 266 / 1.50 → 10.1° | steep |
| 10 | 819 / 3.29 → 14.0° | steep | 811 / 3.02 → 15.0° | steep |
| 11 | 276 / 1.60 → 9.8° | gorge | 260 / 1.55 → 9.5° | gorge |
| 12 | 88 / 2.95 → 1.7° | gorge | 0 | — |

- **No wall (> 45°) under either instrument.**
- The six > 200 m (amended): **3 gorges (8.5–9.7°), 3 steep (10.1–15.0°)**.
- The delivered world's steps (Q1): 3.6–28.6°, one gorge, the rest steep, no wall.
- The mean slope hides the local cell-to-cell slopes inside a zone; those are not read here.

## Q4 — the 12 below-sea-basin lakes' spillways (ON extended)

Their outflow is a `Spillway` segment, not a D8 outlet.
- **9 of 12 have a spillway path.** 1000002 and 1000008 (both at 76.3 m) and 1000010 (0.4 m) have no spillway
  segment of their own.
- The instrument is F136-K6b from the lake's level, along the spillway, its downstream chain, then the D8.
- Results:
  - **1000018**: level 471.8 m; drop to the first graded reach **307 m** at 4.46 km (**3.9°, gorge**); total to
    the sea 473 m.
  - **1000019**: level 614.9 m (F114's resistant "614.9 m lake"); **348 m** at 3.22 km (**6.2°, gorge**).
  - **1000001**: level 176.8 m; 87 m at 0.90 km (5.5°, gorge).
  - **1000011** (Unresolved): 0.8 m, but 438 m in total over 30 km (graded).
  - 1000012, 1000013, 1000017: 0 (graded from the col); 1000014, 1000016: the path reaches the sea before 2 km.
- **Drops > 50 / > 200 / > 500 m: 3 / 2 / 0, all gorges.** The spillways' descents are long and gentle, not
  steps.

## A — the age in the code (`grep_rule11_age.txt`)

**1. How the age enters today.**
- **In the construction**:
  - through **k**: `floor_m = base + age_k·χ` (`valley_construction.rs:371`), F121's k = 0.07186 → 0.07183 at the
    canonical framing (PROXY, calibrated on F95's oracle p50 488.1 m);
  - through **W(k)'s γ**: `width_m` (`:325`, `(k/k₀)^γ`, gated, `None` by default);
  - through **the light pass's slope floor S_eq** when the closure is set. The C2 modes force 0.024
    (`workspace.rs:945–949`).
- **In the delivered world**: through **S_eq** ("Âge du continent — constante (F109)", absolute slope floor + A1)
  and the incision's **k_time** (k·dt·iterations).
  - K and the duration are not separately observable (F44 / F90).
  - **F121-C3**: "the factor is not the age": ×1.71 of k buys ×1.079 of relief.

**2. A break position or a retreat?**
- **None.** No knickpoint position and no retreat state exist anywhere.
- F115 measures a celerity (`stream_power.rs:1028`, `K·A_max^m` = 180 618 m/yr at the shipped K, Courant 3 699)
  as a diagnostic.
- H-2 (sill incision, ADR:2667) was computed as zero by F114 and never built.

**3. How the age is shown to the author.**
- "âge k": a segmented selector **×0.7 / ×1 / ×1.4** × `F121_AGE_K` (`workspace.rs:1210–1220`), "the three ages
  MEASURED in F121-C3".
- "Âge du continent — constante (F109)": a checkbox with `s_eq` ∈ **{0.021, 0.024, 0.027}**, the measured window
  only, no free slider (`:1116–1145`).
- W(k)'s **γ ∈ {0, 0.5, 1}** (`:1231`).

## S — the specification

`spec_gorge_age.md`, for the author's and the reviewer's review **before** any construction:
- a retreat parameter r ∈ [0, 2]: r = 0 young (full lake, gorge from the col); r = 1 today's ON; r = 2 emptied;
- its link to age, a DECISION: the selector ×0.7 / ×1 / ×1.4 → 0 / 1 / 2, or an independent selector;
- the lake level L(r) between L_in, L_bed and L_floor;
- **the rim invariant**, closing F136-K7's one-sided invariant: `floor_outlet = max(law_below, z_gorge)`;
- the descent's slope S_g = 14.5° (MEASURED: the delivered world's median, a PROXY);
- waterfalls on a drop ≥ H_f within ~100 m (DECISION), C-3 only at its contrasts, and the `features` tag in
  `rivers.json`;
- the spillways;
- gates from negative controls, the first being **r = 1 = today's ON bit for bit**;
- the viz toggle and binary questions;
- a cost estimate < 2 s per world.

## Predictions

**Mine:**
- **P-Q1**: F114's world and stage held. "≥ 3 steps > 200 m in today's delivered worlds" half: livré 6 held,
  A1+B2 2 refuted. "Fewer than ON extended's 6" refuted (livré = 6). "OFF hid it" held.
- **P-Q2**: an erodibility multiplier, tectonic and volcanic, weak with the relief, not with the climate, read by
  the incision and the light pass, not by `carve`: held. "Correlation length > 20 km" held for log K and in y; for
  K in x it is 12.5–25 km, so undecided.
- **P-Q3**: held. No wall; the amended ones 8.5–15.0°.
- **P-Q4**: held (2 of 12 > 200 m; 6 of 9 measured < 100 m).
- **P-A**: held (k, γ, S_eq; no break position).
- **Meta**: of my four disagreements (Q2's bound, Q3, Q4, A), none is clearly wrong. **"At least two wrong" not
  held.**

**The reviewer's:**
- **Q1** held: livré has 6 steps > 200 m. F114's own world today has 2.
- **Q2**: "not a lithology, unusable alone" held. "Correlation length > 50 km" half: y > 100 km, x 12.5–25 km for
  K.
- **Q3** ("≥ 4 of 6 walls") **refuted**: 0 walls.
- **Q4** ("≥ 3 of 12 > 200 m") **refuted**: 2.
- **A**: "only through k" **refuted** (γ and S_eq); "no break position" held.
- **Meta** held.

## Limitations, stated

1. **Q1's instrument cannot close long steep descents.** "No graded reach within 10 km" reads NaN. The "drop over
   5 km" column is given beside it.
2. **Q3 reads mean slopes**; the steepest single cells inside a zone are not read, so a wall a few cells high
   inside a "steep" zone is possible.
3. **Q2's correlation lengths** are brackets (lags double).
4. **3 below-sea lakes have no spillway segment of their own** (chained or merged, not traced).
5. **The specification's L_floor** (r > 1) is not measured yet.

## State

**Uncommitted**, awaiting the go-ahead: `f126_coast.rs` (`f137_q`), this folder (with `spec_gorge_age.md`), and
ADR Finding 137. **No production code. Nothing built.**
