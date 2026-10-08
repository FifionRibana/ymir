# Finding 157 — the gorge paused (closing summary); the lakes' outline on the hillshade; the lake base alone ready for the author's look: F133's table reproduced to the digit, the same 11 new lakes and crops, 6 steps > 200 m below the based lakes (two above 40° over 100 m), both worlds exported for Living Landz, no measurable cost; and the production world's defects inventoried (planar walls, anisotropy and comb teeth are all born at the construction)

**Status: the viz outline built; the rest measured; nothing promoted; the gorge paused. Nothing committed.**

**Raw outputs, in this folder:**
- `f157_viz_raw.txt` (B2, B3, B5);
- `f157_b_raw.txt` (the cost, B4);
- `f157_b7_raw.txt` (B7);
- `f157_b1_f135_w3m_raw.txt` (B1);
- the 24 `f157_{crop}_{off,on}.png`;
- `etat_gorge_pause.md`, `breach_spill_anchor.patch`;
- `checks_before.txt`, `checks_after.txt`.

**Process:**
- The predictions were written before any measurement (`f157_predictions.md`, with the non-blind items declared).
- The instruments were declared before the benches (`f157_declared.md`).
- **One amendment came after a run** (B7 split out of `f157_b` after a memory abort), declared there as NOT blind.

## The author's criterion, recorded as given

- **On the gorge, after the T / T-all / P0 / NP images**: « Je ne suis pas convaincu. Je ne vois pas où s'arrête le
  lac, mais je vois un rebord qui apparaît en T-all par rapport à P0 (et en T d'ailleurs). »
- **Decisions of 2026-10-08**:
  - the gorge: « Pause, avec bilan de clôture »;
  - the lake base alone: « La valider pour la promouvoir (crops restants + Living Landz) »;
  - seed 20261008002's crash: « Le laisser dans la file ».
- **The reviewer's reading** (a hypothesis to check): the visible rim is the ring's clamp (G-ring). P0 erases it
  because the light pass erodes the ring freely. Once the design is protected (T, T-all) the edge stands. No setting
  of the light pass gives a held rim without an edge: the defect is in the construction.
- **Limit**: 2 rounds for the lake base. This one prepares; F158 promotes after the author's look.

## Partie 0

- **F156 committed** (`e796eb1`).
- **The breach's production patch was not applied.** It is saved as
  [`breach_spill_anchor.patch`](breach_spill_anchor.patch) (5 files: `RampAnchor`, the skeleton's explicit floor
  anchor, the three ALGO bumps, `viz_hd_lakes_on`, `hd.rs`).
  - It was generated on the clean tree and the files restored.
  - It was not compiled, and its doc comments still say "confirmed blind", which F156 refuted.
- **Checks before** = F156's after-checks: no `src` change between them and this round's start.
- **Checks after**: `cargo check` clean; lib 621; viz **32** (+1, the outline's test); **the guard 6 / 6, field and
  lakes = banc** (C2 /10 col `a8d2d538d692c2f0`).

## C — the gorge paused

- **[`etat_gorge_pause.md`](etat_gorge_pause.md)**: what is built (v3–v6, the falls' rule, the 18 bodies, the plain,
  the retreat per lake, T and T-all), what holds, what fails with its cause, the eight candidates set aside, the
  track for a resumption (a hypothesis), and the queue.
- **`accepted_defects.md` is unchanged.** The gorge is paused, not accepted: no defect enters.
- An ADR closing entry heads ADR Finding 157.

## V — the lakes' outline on the hillshade (viz)

- **Relief → Ombrage** draws the outline of `run_hd`'s final lake mask:
  - a **cyan** line on the lake cells touching another cell (4-neighbours);
  - a **dark navy** line on the shore cells touching a lake.
- Two-tone, so it reads on a light slope and on a dark one alike. It is drawn before the river overlay, which stays
  on top.
- **A « Contour des lacs » box** in the « 🗻 Ombrage ▾ » menu, on by default, hides it.
- A view: the world does not change.
- Permanent test `the_outline_rings_each_lake_inside_and_out`, with its negative control (no lake, no outline).
- **The images below are drawn with it**: the hillshade paints a lake like land, and the outline shows where it stops.

## B — the lake base alone

### B1 — F133's table, replayed on the current code (`f135_w3m`, unchanged)

| | OFF | ON F132 | ON extended |
|---|---|---|---|
| canyons | 0 | 0 | 0 |
| spurs near a coastal wall | 0.0142 /km | 0.0142 | 0.0142 |
| R8 terrain | 0.0452 | 0.0395 | 0.0411 |
| relief p50 | 515.4 m | 551.9 m | 546.6 m |
| Δz between flat resolutions (cells) | 1 150 588 | 235 602 | 12 509 |
| lakes (area) | 26 (3 492.8 km²) | 34 (4 454.0 km²) | 34 (4 446.2 km²) |
| θ on the carved links | 0.493 [0.492, 0.495] | 0.436 [0.427, 0.445] | 0.439 [0.431, 0.446] |

- **Identical to F133 and F135-W to the digit**, lakes 26 / 34 included.
- OFF reads the guard's field and lakes (Match). **Nothing to attribute.**

**The cost** (`f157_b`, single timings):
- `build_world` OFF 151 s, ON 147 s;
- the run_hd tail OFF 100 s, ON 98 s.
- **ON costs nothing measurable** (−3 s, within the noise).

### B2 — the crops (`f157_viz`, on `run_hd` itself)

- **The viz path lists 26 / 34 lakes.**
- **The same 11 new lakes as F133v, at the same centres**: lakes 1, 2, 3, 4, 7, 8, 9, 10, 11, 12, 19.
- 44 268 342 conditioned cells differ between OFF and ON.
- **No 615² window has zero changed cells.** With the declared |Δz| ≤ 1 m threshold, the control crop is at origin
  (2688, 5888), 78 % land, 3 051 changed cells.

| crop | centre (data) | origin (615²) | lakes OFF | lakes ON | changed cells |
|---|---|---|---|---|---|
| lake 1 (132.8 km²) | (3491, 2500) | (3184, 2193) | 0 | 1 | 198 226 |
| lake 2 (35.5 km²) | (4019, 2839) | (3712, 2532) | 0 | 2 | 305 223 |
| lake 3 (159.7 km²) | (4124, 3096) | (3817, 2789) | 0 | 2 | 271 278 |
| lake 4 (45.8 km²) | (2220, 3584) | (1913, 3277) | 7 | 9 | 82 772 |
| lake 7 (5.3 km²) | (1906, 3697) | (1599, 3390) | 1 | 3 | 130 034 |
| lake 8 (12.3 km²) | (4060, 4008) | (3753, 3701) | 3 | 4 | 189 167 |
| lake 9 (79.7 km²) | (4331, 4107) | (4024, 3800) | 4 | 5 | 129 053 |
| lake 10 (249.0 km²) | (2951, 4213) | (2644, 3906) | 0 | 1 | 275 562 |
| lake 11 (163.0 km²) | (1429, 4287) | (1122, 3980) | 1 | 1 | 208 717 |
| lake 12 (10.2 km²) | (3468, 4693) | (3161, 4386) | 3 | 3 | 109 537 |
| lake 19 (72.0 km²) | (3342, 5798) | (3035, 5491) | 2 | 2 | 171 536 |
| control | (2995, 6195) | (2688, 5888) | 1 | 1 | 3 051 |

### B3 — the images for the author

- `f157_{lake1,…,lake19,control}_{off,on}.png`:
  - 615 × 615 cells (30 km);
  - north up;
  - Relief → Ombrage with the lakes' outline, the viz's own buffer.

### B4 — the steps below the based lakes (ON)

The amended instrument (F136-K6b completed). For each lake, the drop to the first graded reach ≥ 2 km, its length,
the mean slope, and the max slope over 100 m.

| lake | km² | level | drop | length | mean | max over 100 m | crop |
|---|---|---|---|---|---|---|---|
| 1 | 132.8 | 327.2 m | **290.1 m** | 1.94 km | 8.5° | 27.3° | lake 1 |
| 2 | 35.5 | 1 292.0 m | **862.1 m** | 3.31 km | 14.6° | **43.6°** | lake 2 / 3 |
| 4 | 45.8 | 168.2 m | 119.0 m | 0.39 km | 16.8° | 22.4° | lake 4 / 7 |
| 7 | 5.3 | 383.4 m | **316.6 m** | 1.85 km | 9.7° | 31.3° | lake 7 |
| 8 | 12.3 | 520.5 m | 51.3 m | 0.12 km | 23.5° | 23.5° | lake 8 / 9 |
| 9 | 79.7 | 418.0 m | **265.9 m** | 1.50 km | 10.1° | 27.0° | lake 9 |
| 10 | 249.0 | 974.5 m | **810.7 m** | 3.02 km | 15.0° | **47.1°** | lake 10 |
| 11 | 163.0 | 362.0 m | **260.4 m** | 1.55 km | 9.5° | 27.8° | lake 11 |
| 20 | 33.2 | 127.6 m | 23.9 m | 0.12 km | 11.5° | 11.5° | — |
| 3, 12, 13, 19 | | | 0 (graded at once) | | | | |
| 15 | | | no graded reach within 10 km | | | | |

- **6 steps > 200 m**: lakes 1, 2, 7, 9, 10, 11. They match F136-K6b (290 / 862 / 260 m for lakes 1, 2, 11).
- Two exceed 40° over 100 m: lake 10 at 47.1° and lake 2 at 43.6°.
- The mean slopes over the steps are 8.5–15.0°.
- 12 lakes have no D8 receiver: below-sea basins draining by spillway.

### B5 — the Living Landz exports

- **OFF**: `exports/f157/off/seed10481999410520546993_8192.ymir/`
- **ON**: `exports/f157/on/seed10481999410520546993_8192.ymir/`
- Both are written by `run_hd` itself, as the production world: the guard's C2 /10 col literal, container 1.1.0 with
  rivers_ll and the geology, ≈ 1.1 GB each.
- `exports/` is git-ignored.

### B6 — what the promotion would change (not made)

- **The production field**: ON's eroded field differs (guard key `834fa66b41d6e05d`, field `a620a9aa0d5882a0`, against
  OFF's `a8d2d538d692c2f0`). The construction's χ stops at the based lakes.
- **The guard**:
  - promoted as the default state's lake base: the « C2 /10 col (défaut) » entry, field and lakes, regenerated
    (`f123_guard`, `f134_lake_guard`);
  - promoted inside `ValleyConstruction::new`: also C1 nue, C2 /3 and C2 /10 niveau mer.
  - livré and A1+B2 do not move.
- **The lake references**: 26 → 34 lakes, 3 492.8 → 4 446.2 km².
  - The F134 lake guard's hashes change.
  - The lake ids the dossier quotes are ON's (F139's numbering), as now.
- **The exports**: the lakes, the rivers (rivers_ll), the biomes, and the geology (the lake margins' loose deposits,
  the evaporites) move with the lakes.
- **The cost**: none measurable (B1).

### B7 — the production world's defects (OFF, the 12 crops, stage by stage)

**The stages**: S1 (tectonics + FBM) → the construction → C1 bare (no light pass) → OFF (production).
- Production has **no droplet pass** (relief-v3 sets `cfg.erosion = None`; `upscale.rs:410`).
- So C1 bare equals the construction on land, bit for bit on every crop. The light pass is the only stage after the
  construction that acts on land.

| (mean over the 12 crops) | S1 | construction | C1 bare | OFF |
|---|---|---|---|---|
| planar walls (land at 28° ± 0.5°) | 0.10 % | **17.15 %** | 17.15 % | 9.07 % |
| terrain R8 | 0.0428 | **0.0933** | 0.0933 | 0.0876 |
| comb teeth (F124: a reach ≥ 80 % on the construction's walls) / reaches | — | — | 428 / 6 684 | 371 / 7 267 |

- **Planar walls and flat facets**: born at the **construction** (28° planar walls, F120).
  - From 0.1 % to 17 % of the land per crop on average, up to 48 % (lake 2's crop).
  - **The light pass halves them** (9.1 %).
- **Axis alignment** (the terrain R8): born at the **construction** (×2.2 from S1). The light pass barely moves it
  (−6 %). F126 traced the sub-valleys' stripes to the D8 lattice of the construction's trunks.
- **Comb teeth**: by definition on the construction's walls, so born there.
  - In OFF, in every crop: 1–127 per crop, 5.1 % of the reaches.
  - The light pass removes some (428 → 371).
- **The network R8 on the 1–3 km² reaches** is not read: too few such chords per crop (NaN on 5 of 12 in OFF).
- No remedy. **The inventory of a possible work**, all at the construction stage: planar walls, the trunks' D8
  lattice, the teeth on the walls.

## The author's questions for the look (copied as given)

Par crop, dans le viz puis dans Living Landz :
1. **Amont** : la vallée arrive-t-elle au lac à hauteur de sa surface ?
2. **Aval** : en sortie de lac, la descente te semble-t-elle acceptable (gorge, rapides), ou vois-tu un mur ou une
   marche artificielle ?
3. **Non-régression** : vois-tu en ON un défaut absent en OFF ?
4. **Les lacs ajoutés** te conviennent-ils ?

## Predictions

**The reviewer's (hypotheses to check; judged):**
- B1, "the table is reproduced to the digit (26 / 34); nothing to attribute": **held**.
- B4:
  - "≥ 4 steps > 200 m": **held** (6);
  - "≥ 1 beyond 40° over 100 m": **held** (47.1°, 43.6°);
  - "all classed 'steep' or 'gorge' by mean slope": **not judged**. The classes are not defined; the means are
    8.5–15.0°.
- V, "the guard stays 6 / 6": **held**.
- Meta (« au moins une est fausse »): **refuted**. All four judged predictions held; B4's classes were not judged.

**Mine (`f157_predictions.md`):**
- **V**: P-V1 (the guard 6 / 6): **held**.
- **B**:
  - P-B1 (to the digit, 26 / 34): **held**;
  - P-B2 (the same 11; no zero-change control; one at ≤ 1 m): **held**;
  - P-B4 (8 / 6 / 2; ≥ 4 > 200 m; ≥ 1 > 40°; means 2–15°):
    - **held for > 200 m (6)**, > 40° and the means (8.5–15.0°);
    - "> 50 m: 8": 8 measured over 50 m (51.3, 119.0 and the six), **held**;
    - "> 500 m: 2": **held** (lakes 2, 10);
  - P-B5 (the containers with rivers_ll and the geology): **held**;
  - P-B6:
    - "the field hash changes": **held**;
    - "4 of 6 states move": **held only if** the promotion goes inside `ValleyConstruction::new`; 1 if only the
      default state's;
    - "0–10 s": **held** (−3 s);
  - P-B7:
    - "teeth in ≥ 6 crops": **held** (all 12);
    - "planar 5–20 % per crop": **refuted in its range** (1.5–22.6 % in OFF; 2.7–48 % at the construction);
    - "R8 0.03–0.10 per crop": **refuted in its range** (0.034–0.120);
    - "born at the construction, the light pass cutting the alignment": **held for the walls and the teeth**; the R8
      is cut only −6 %.

## Limitations, stated

1. **One world, one seed** (PSEED). The crops are the témoin's.
2. **The cost** is single timings (±30 s).
3. **B7's network R8** is not readable on these crops.
4. **The patch** is unapplied and uncompiled.
5. **The memory**: `f157_b`'s B7 aborted at 512 MB; B7 was rerun alone (amendment).

## State

**Uncommitted**, awaiting the feu vert.
- `crates/ymir-viz/src/ui/workspace.rs`:
  - the lakes' outline (`lake_outline_mask`, `draw_lake_outline`, `ReliefView::Shade { outline }`);
  - the « Contour des lacs » box and its test;
  - the bench `f157_viz` (with `run_with`).
- `crates/ymir-core/tests/f126_coast.rs` (`f157_b`), `tests/f157_b7.rs`.
- ADR Finding 157 with the gorge's closing entry; this folder.
- The exports in `exports/f157/` (git-ignored).
