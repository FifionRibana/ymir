# Finding 152 — geology v1, round 3 of 3: the rocks, the favourabilities and the views built; the gates; and what the author will see

**Status: built.** Two gates of note:
- **G-blocks failed in run 1** and holds after a declared, non-blind amendment (the spec's terrain snapping);
- **G-veins makes no districts**: reported, not fixed.

Nothing committed.

**Files:**
- predictions: `f152_predictions.md`, written before any measurement;
- instruments and amendment: `f152_declared.md`;
- raw output: `f152_geology_raw.txt` (run 2, final) and `f152_geology_raw_run1.txt` (run 1: G-blocks failed);
- the views: `view_*.png`;
- checks: `checks_before.txt` and `checks_after.txt`.
- Code: `crates/ymir-core/src/geology/`. The rules: `crates/ymir-core/data/geology_rules_v0.toml`. The format:
  `docs/geology_format.md`.

## The author's criterion and decisions, recorded as given

- **F151's decisions** (2026-10-06):
  - production = the guard's world;
  - fossil belts recorded;
  - the 12 resources;
  - a rock grid plus a rules file;
  - the viz's « Roches » and « Chances » views.
- **The decisions of 2026-10-06**:
  - « Favorabilité dans Ymir, rareté dans Living Landz »;
  - « Oui, moduler les filons par la fracturation (districts) ».
- **Proposed by the reviewer, not contested by the author**:
  - a PNG export (lossless);
  - TOML rules;
  - the fossil-belt rule kept as a PROXY (its duration threshold is untested on the témoin, where every collision
    lasted 48 steps);
  - a build with the v0 rules, which the author then adjusts by eye.
- **The consequence**: the exported value is a **relative favourability** (0–100), not a probability. Living Landz
  applies a global rarity per resource.
- **The limit**: the last round of geology v1. After the author's look the work is closed; the author then adjusts the
  rules file alone.

## Part 0

- F151 committed (`5a124bc`), the observer hook included.
- **Checks after** (`checks_after.txt`):
  - `cargo check --workspace` clean;
  - lib 616 passed (610 + 6 geology tests);
  - viz 31 passed;
  - **the guard 6 / 6, field and lakes = banc, through `run_hd` with the geology stage running**: G-field holds.
  - In `run_hd` the stage takes **2.97–3.32 s** per state (the bench 2.58 s; `run_hd` adds the inputs'
    conversions).

## B — what was built

| part | where | what |
|---|---|---|
| fossil belts in production | `geology::history` | the read-only observer records collisions and sutures in a geology pass of its own (0.35 s). The terrain's cached run is untouched. F151's PROXY rule gives 106 fossil coarse cells on the témoin |
| the structural density | `TectonicHistory::structural_density` | `max(C-3b, exp(−d / 25 km))` from the sutures and the past collisions; the zoning only. **A test pins that C-3b's own density is unchanged by building it.** |
| the rock grid | `geology::rocks` | 9 classes, u8 per HD cell, with hardness and colours. Contacts bilinear → Gaussian 3 km → **snapped to the terrain** (the amendment) → 0.5. Volcanic rock bound to the cone's relief |
| the rules reader | `geology::rules` | TOML; unknown keys, rocks, belts and settings refused; `modulate = "structural"` with `[modulation.structural]`. **d_ref = 0.994, MEASURED** (the p90 of d over the belts and the arc zone; p50 0.664). **g_min = 0.2, DECISION** |
| the zoning | `geology::zoning` | the context on the ¼ grid, built once per world. AND within a rule, MAX across rules, placers downstream along the D8 rivers |
| the export | `geology::export`, the container 1.1.0 | `geologie_roches.png` (full resolution), `favorabilite_<id>.png` (¼), `geologie.json` (the rules' SHA-256, the classes, the resources). The existing layers are unchanged |
| the zoning cache | `GeologyView` in `HdResult` | the context stays with the world; « Recharger les règles » reruns the zoning only |
| the viz | `workspace.rs` | the « ⛏ Géologie » menu: « Roches » (with legend and hardness), « Chances » (one resource or all, the highest wins, opacity = favourability, legend), the rules file's path, « Recharger les règles » with the re-zoning time, a refused file's error |

The vein rules carrying `modulate`: gold (belt, arc), silver (arc, belt), tin (belt), copper (arc). Not the
stratiform deposits, the placers or the surface resources.

**The rock classes on the témoin** (share of the non-water cells):
- basement 62.6 %;
- belt 14.2 %;
- craton 8.4 %;
- basaltic volcanic 5.5 %;
- rift fill 3.8 %;
- loose deposits 2.9 %;
- evaporites 2.4 %;
- arc volcanic 0.3 %.

The loose deposits are smaller than F151's 11.6 %: the construction's floor is not available in `run_hd`, and the
floodplain proxy is narrower.

## G — the gates (run 2, final)

| gate | measure | verdict |
|---|---|---|
| **G-field** | the 6/6 guard (field and lakes) through `run_hd`, which now runs the geology stage; the C-3b test; the geology reads the world and writes its own buffers | **PASS** (6/6 = banc) |
| **G-blocks** | the share within 1 cell of a coarse grid line, against the land control (0.0274): all rock contacts **1.10 ×**; the coarse-sourced ones (craton / belt / rift / basement) **1.05 ×**; the ¼ favourability edges **0.98 ×** | **PASS** (≤ 1.5). **Run 1 failed**: coarse-sourced 4.01 × |
| **G-circles** | 14 volcanic patches ≥ 20 cells; Q / Q_disc max **0.704** | **PASS** (≤ 0.8) |
| **G-nonempty** | land ¼ cells with favourability ≥ 10 (table below) | **PASS**: all 12 non-empty |
| **G-placers** | the ¼ cells a placer raised that hold a river of ≥ 5 km²: gold 1 745 / 1 745, tin 471 / 471, gems 1 325 / 1 325 | **PASS** (100 %) |
| **G-rules** | the permanent test `changing_one_rule_changes_only_its_resource`, with a negative control | **PASS** |
| **G-veins** | below | **districts: NOT formed** |
| **the cost** | the stage 2.58 s in the bench (history 0.35, rocks 0.94, context 1.23, zoning 0.05), **2.97–3.32 s in `run_hd`**; re-zoning alone **0.05 s** | +1.2–1.3 % of 249.8 s |

**G-nonempty** (land ¼ cells with favourability ≥ 10):

| resource | area | max | placer HD cells |
|---|---|---|---|
| stone | 19 409 km² | 60 | — |
| gold | 6 404 km² | 35 | 4 876 |
| silver | 5 692 km² | 25 | — |
| copper | 5 525 km² | 40 | — |
| iron | 4 227 km² | 40 | — |
| gems | 3 568 km² | 20 | 3 482 |
| sulphur and obsidian | 1 721 km² | 60 | — |
| tin | 1 452 km² | 30 | 1 342 |
| sand and gravel | 735 km² | 50 | — |
| salt | 600 km² | 80 | — |
| clay | 483 km² | 50 | — |
| **peat** | **165 km²** | 60 | — |

**G-veins**: the modulated favourability inside the belts.

| resource, belts | p10 / p50 / p90 | share ≥ half the max | the control (C-3b alone): mean change |
|---|---|---|---|
| gold, fossil belts (37 847 ¼ cells) | 33 / 34 / 35 | **100 %** | **−76.4 %** (mean 34.1 → 8.0) |
| tin, fossil belts | 28 / 29 / 30 | 100 % | −76.5 % |
| gold, all belts (94 444) | 23 / 32 / 35 | 100 % | −45.3 % |
| tin, all belts | 0 / 0 / 30 | 40.1 % | −76.5 % |

- **What the sutures bring**: without them, the fossil-belt gold and tin fall to g_min. The sutures are what make the
  fossil belts favourable at all; C-3b's density there is near 0, since today's contacts are far.
- **No districts**: inside a belt the favourability is a flat sheet (100 % above half its maximum). The belts are
  drawn from proximity to the same contacts and sutures that make the density, and the density varies on a 25 km
  scale, wider than a belt.
- **Ymir has no structure finer than a belt** (no fault lineaments), so the modulation separates belts from one
  another but cannot make districts within one.
- **Not fixed** (a change of method, not of a parameter). Candidates for a later round: a fine-scale structural
  proxy, for example the HD drainage's straight segments as lineaments; or accept flat belts.

## The run history

- **Run 1**: G-blocks failed on the coarse-sourced contacts (4.01 ×).
  - The cause is geometric: a long straight edge of a coarse mask stays straight, and on the grid line, under any
    symmetric smoothing.
- **Amendment** (declared before run 2, NOT blind): the spec's optional terrain snapping. Hard classes take `+β·a`,
  the soft rift `−β·a`, with `a` the terrain anomaly at ≈ 1.6 km and β = 0.25, so a contact moves within ≈ 2 km and
  runs along the valleys (PROXY: differential erosion).
- **Run 2**: G-blocks 1.05 ×; every other gate unchanged; the cost +0.26 s.

## The export (the témoin)

- `geologie_roches.png` 0.65 MB;
- the 12 `favorabilite_*.png` **0.57 MB**;
- `geologie.json` 3.7 kB.

## R — what the author will see (the views, bench renders with the viz's colours)

- The captures are `view_rocks.png`, `view_chances_all.png`, `view_chances_{gold,tin,copper,salt}.png`, and
  `view_chances_gold_CONTROL_c3b_only.png`.
- They use the same `geology::render` functions as the viz, over a grey shaded relief. The viz draws them over its
  hypsometric relief.
- **The viz itself was not screenshotted** (no GUI automation here): the author's look is in the viz.

The binary questions, answered as I see the renders (to confirm by the author):
- **Blocks or perfect circles?**
  - Blocks: **no** since the amendment (in run 1 the craton had straight east–west edges).
  - Circles: **a few coastal basaltic patches still read round**, although their Q ratio passes (≤ 0.42).
- **Vein resources in patches rather than sheets?** **No**: gold and tin are uniform bands along the belts (G-veins).
  The placers are thin lines along the rivers.
- **A resource absent, or everywhere?** **Stone is nearly everywhere** (19 409 km², favourability 40–60), so in the
  « all » view it hides the others wherever it wins. **Peat is the smallest** (165 km²). None is absent.

## Predictions

**The reviewer's (hypotheses to check):**
- **G-veins' control, "−80 % or more"**: REFUTED (−76.4 %). Declared before measuring: with g_min = 0.2 the control
  cannot pass −80 %.
- **G-nonempty**: "all 12" HELD; "copper the smallest" REFUTED (peat, 165 km²; copper 5 525 km²).
- **"G-blocks and G-circles hold"**: G-circles HELD; G-blocks FAILED in run 1 and held after the declared amendment.
- **Export < 5 MB**: HELD (0.57 MB).
- **Cost**:
  - "< 3 s": HELD in the bench (2.58 s), REFUTED in `run_hd` (2.97–3.32 s);
  - "re-zoning < 2 s": HELD (0.05 s).
- **Meta**: HELD.

**Mine:**
- **P-veins**:
  - "the control −65 to −80 %": HELD (−76.4 %);
  - "20–40 % above half: districts": REFUTED (100 %, no districts).
- **P-nonempty**:
  - "all 12": HELD;
  - "salt or sulphur the smallest": REFUTED (peat).
- **P-blocks, "holds"**: REFUTED in run 1.
- **P-circles, "at least one young cone fails"**: REFUTED (none).
- **P-export, "1–3 MB"**: REFUTED (0.57).
- **P-cost**:
  - "2–6 s": HELD (2.58 in the bench, 2.97–3.32 in `run_hd`);
  - "re-zoning 0.3–1.5 s": REFUTED (0.05).
- **Meta**: HELD.

## For the author

1. **Look in the viz**: « ⛏ Géologie » → « Roches », « Chances ».
2. **Tune the rules file**: copy `crates/ymir-core/data/geology_rules_v0.toml` to `geology_rules.toml` in the viz's
   working directory, edit it, and press « Recharger les règles ». Without that file the viz uses the shipped v0. Two
   first candidates:
   - lower stone (it hides everything in the « all » view);
   - raise peat's reach.
3. **Districts**: accept flat belts, or ask for a fine-scale structural proxy in a later round.
4. **The geology v1 work closes after your look.**
