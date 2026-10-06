# Finding 151 — geology v1, round 2 of 3: the specification, the fossil-belt prototype, and the viz aligned on production

**Status: specified and prototyped; built only the viz's two small changes and a read-only C1 observer; nothing
committed.**

**Files:**
- predictions: `f151_predictions.md`, written before any measurement;
- instruments: `f151_declared.md`;
- raw output: `f151_geology_raw.txt` (bench `crates/ymir-core/tests/f151_geology.rs`);
- the maps: `map_f_*.png`, `map_s_class_sources.png`;
- the specification: `spec_geologie_v1.md`;
- the rules proposal: `geology_rules_v0.toml`;
- checks: `checks_before.txt` and `checks_after.txt`.

## The author's criterion and decisions, recorded as given

- Usage : « Chutes + ressources. […] Là, on va faire du zoning : ici, il y a de l'or possible ; ici, du fer, etc. »
- **The decisions of 2026-10-06**:
  - **production**: « comme la garde, car je l'active manuellement » (volcanism, lithology and fracture on);
  - **fossil belts**: « Oui, enregistrer les anciennes collisions »;
  - **the v1 resources**: fer, or, argent, cuivre, étain, pierre, argile, sable et gravier, sel, tourbe, soufre et
    obsidienne, gemmes;
  - **zoning**: « Je pense qu'une règle dans un fichier est plus versatile, et donc une grille des roches. La roche
    n'est pas une variable du monde. Les ressources, elles, le sont. J'ai d'ailleurs une mécanique de
    "réinitialisation" de certaines parties lorsque remasquées par le mist. Donc la répartition n'est qu'une image à
    l'instant t des chances. »;
  - **viz**: « Dans Ymir-viz, je veux avoir une visualisation des roches (comme les biomes), et pour les chances, juste
    générer une vue des chances (avec opacité = % chance). Comme ça, on verra une image avec des patchs de couleurs
    pour les ressources. »
- **The consequence (the reviewer's)**: Ymir exports CHANCES, never deposits. Living Landz draws the deposits in the
  chances at each reset.
- The limit: round 2 of 3. The cost: seconds per world against 249.8 s.

## Part 0

1. F150 committed (`267f868`, not pushed).
2. **The viz starts as production.**
   - C-2, C-3 and C-3b are checked by default (`workspace.rs`); H-1 is not.
   - A permanent test, `startup_defaults::the_workspace_starts_with_the_guards_closures`, pins the defaults to the
     guard's literal.
3. **The panel's three fields are wired.**
   - « Type de plaque » and « Craton » come from `HdResult.tectonic`, sampled with the overlay's own mapping.
   - « Épaisseur crustale / S̃ » comes from a new `crust_s` in the debug labels. It is non-dimensional: no km scale
     exists in the code.
4. **Checks after** (`checks_after.txt`):
   - `cargo check --workspace` clean;
   - lib 610 passed;
   - viz 31 passed (30 + the defaults test);
   - **the guard 6 / 6, field and lakes = banc** (C2 /10 col `a8d2d538d692c2f0`).
   - With the defaults test, the viz's startup world is the guard's.

## F — the fossil belts (prototype)

- **The observer**: `time_loop::run_with_closures_observed`, whose callback also reads the kinematics.
  `run_with_closures` delegates to it.
  - **A production-code change of the signature only, declared**: the brief asked for no production code in F, but the
    plate velocities change at every merge and split, so a past boundary cannot be classified without them.
  - **The final state is bit-identical** with and without the recording (`s`, `age`, `plate_id`, `plate_type`,
    `cratonic_mask`). The recording costs **+0.02 s** (0.20 → 0.22 s).
  - The guard: below.

**On the témoin (300 steps):**

| measure | value |
|---|---|
| cells ever convergent | 520 (162 continental at the end) |
| cells ever in collision (C-C) | **45**, all continental, every one for exactly 48 steps, ending at step 48 |
| cells ever on a subduction upper plate | 123 |
| cells ever divergent | 345 |
| convergent at the end · collision at the end | 154 · 0 |
| plate merges (sutures) | **6**: five at step 49 (3 → 2, 5 → 1, 7 → 0, 6 → 0, **4 → 0**: 51 cells, 100 % continental), one at step 99 (1 → 0, 226 cells, 17 % continental) |
| rift splits | 0 |
| suture cells | 360, 106 continental at the end |

**What the map shows** (`map_f_ever_collision_and_sutures.png`): the continent's north-western block is an **accreted
terrane**.
- The collision that ended at step 48 merged plate 4 into plate 0, and its suture runs as a diagonal band across the
  continent.
- The rest of the sutures lie at sea.

**The proposed rule (PROXY)**: a **fossil belt** cell was in collision for ≥ 10 steps (≈ 7 Ma at ≈ 0.67 Ma per step,
PROXY from `state.rs`), OR it is a continental suture cell, AND it is not convergent at the end.
- On the témoin: **106 coarse cells**, all continental, 2 on the craton.
- The 45 collision cells are all among the continental suture cells: the collision ended in the merge.
- The duration threshold does not bite on this world: 1 to 20 steps keep the same 45 cells, 50 steps none.

## S — the measurements the spec needs

**The land shares of the rock classes' sources**, with this measurement's thresholds (`map_s_class_sources.png`):

| class | share of the land |
|---|---|
| evaporites (endorheic lake) | 0.95 % |
| loose deposits (valley floor, lake, flat with A ≥ 1 km²) | 11.59 % |
| volcanic (edifice discs) | 4.04 % |
| rift fill | 4.78 % |
| belt (fossil, or within 2 coarse cells of a convergent cell) | 14.25 % |
| craton | 7.70 % |
| **default basement** | **56.69 %** |

**The export sizes**: 12 PROXY chance grids (a field mask × a chance falling 10 points per km), u8.

| resolution | raw | PNG (deflate) |
|---|---|---|
| full (49 m) | 805 MB | 15.0 MB |
| ½ (98 m) | 201 MB | 3.5 MB |
| **¼ (195 m)** | **50 MB** | **1.4 MB** |

- The rock-class grid at full resolution: 67 MB raw, 0.52 MB in PNG.
- The 12 grids with their distance fields take 21.6 s here: bench code, unoptimised. Round 3 measures the real cost.
- **The spec's choice**: chances at ¼, raw in the container. One chance cell covers ~9 hexes, finer than any zone. A
  PNG encoding (1.4 MB) is offered as the fallback.

## The specification

`spec_geologie_v1.md` gives:
- the 9 rock classes, with sources, hardness, colours and priority;
- the contacts: bilinear, Gaussian 3 km, iso-line; volcanic rock bound to the relief, not a disc;
- the rules file: TOML, AND within a rule, MAX across rules, placers downstream;
- the export: rocks at full resolution, chances at ¼, the rules' SHA-256 in the manifest;
- the viz's two views;
- the round-3 gates: G-field, G-blocks, G-circles, G-nonempty, G-placers, G-rules, and the cost.

`geology_rules_v0.toml` gives **34 rules for the 12 resources**:
- 3 placers: gold, tin, gems;
- every chance a PROXY;
- every source a real deposit context cited from memory, to check;
- parsed (TOML 1.0).

## Predictions

**The reviewer's (hypotheses to check):**
- **F**, "≥ 3 sutures, ≥ 100 fossil cells, bit-identical, < 1 s": **HELD** (6 sutures, 106 cells, identical, +0.02 s).
- **S**, "the default class > 50 % of the land": **HELD** (56.7 %).
- **Export**, "¼ compressed < 50 MB": **HELD** (1.4 MB).
- **Viz**, "the guard reads 6/6 at startup": **HELD** (6/6; the defaults equal the guard's literal, tested).
- **Meta**: HELD.

**Mine:**
- **P-F**:
  - "ever convergent 600–1 500": REFUTED (520);
  - "ever collision 50–300": REFUTED (45);
  - "≥ 4 sutures": HELD (6);
  - "bit-identical, < 1 s": HELD.
- **P-S**, "the default 35–55 %": REFUTED (56.7 %).
- **P-export**:
  - "¼ ≤ 10 MB compressed": HELD (1.4);
  - "full compressed 50–150 MB": REFUTED (15).
- **P-viz**: HELD.
- **P-cost**, "< 0.5 s": HELD (+0.02 s).
- **Meta**: HELD.

## For the author (before round 3)

1. **The rules file**: edit `geology_rules_v0.toml` (the chances, the conditions, the colours). Round 3 builds the
   reader for it.
2. **The export**: chances at ¼ raw (50 MB), or PNG (1.4 MB)?
3. **The fossil belt**: keep "collision ≥ 10 steps or a continental suture"?
4. **TOML** as the rules format (a new dependency)?
