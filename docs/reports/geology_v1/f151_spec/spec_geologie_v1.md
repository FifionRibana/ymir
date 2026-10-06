# Geology v1 — specification (Finding 151; built in round 3)

**What this is**: what round 3 builds, and the gates it must pass.

**The provenance of every number**:
- **MEASURED** on the témoin (`f151_geology_raw.txt`, `finding_150.md`);
- **ANCHORED** on a cited source;
- **PROXY**, a stand-in to replace;
- **DECISION**, the author's or this spec's choice, changeable.

**The author's decisions (2026-10-06), recorded as given:**
- **production**: « comme la garde, car je l'active manuellement » (C-2, C-3, C-3b ON). The viz starts that way since
  F151.
- **fossil belts**: « Oui, enregistrer les anciennes collisions ».
- **the v1 resources**: fer, or, argent, cuivre, étain, pierre, argile, sable et gravier, sel, tourbe, soufre et
  obsidienne, gemmes.
- **zoning**: « Je pense qu'une règle dans un fichier est plus versatile, et donc une grille des roches. La roche n'est
  pas une variable du monde. Les ressources, elles, le sont. J'ai d'ailleurs une mécanique de "réinitialisation" de
  certaines parties lorsque remasquées par le mist. Donc la répartition n'est qu'une image à l'instant t des chances. »
- **viz**: « Dans Ymir-viz, je veux avoir une visualisation des roches (comme les biomes), et pour les chances, juste
  générer une vue des chances (avec opacité = % chance). »

**The consequence: Ymir exports CHANCES, never deposits.** Living Landz draws the deposits inside the chances at each
reset.

## 0. The pipeline (round 3)

```
eroded world (unchanged) ──► rock grid (u8, fixed per world)   [cached with the world]
                               │
rules file (TOML) ────────────►└─► chance grids (u8 0–100 per resource)   [cached on (world, rules hash)]
```

- **Read-only**: neither stage changes the terrain, the lakes or the rivers (G-field).
- **The rock grid** depends on the world only. It is recomputed only when the world is.
- **The chances** depend on the rock grid, the context fields and the rules. Changing the rules reruns this stage
  alone (G-rules; its time is measured in round 3).

## 1. The rock classes (`geologie_roches`, u8 per HD cell)

| id | class | sources (Ymir fields) | hardness, for the falls at lake sills | viz colour |
|---|---|---|---|---|
| 0 | none (sea, lake water) | `water_class`, `lake_map` | — | transparent |
| 1 | **craton** (granite, gneiss) | `cratonic_mask` (64²) | **hard** | `#C9A55A` gold |
| 2 | **belt** (metamorphic; fossil or active) | the fossil-belt record (§1.1) + the convergent cells at the end | hard to medium | `#B5487F` plum |
| 3 | **basaltic volcanic** (shields: hotspots, rifts) | edifices with `setting` hotspot or rift (`kind` shield) | **hard** | `#4A4A52` slate |
| 4 | **arc volcanic** (andesite, tuff) | edifices with `setting` arc (`kind` stratocone) | medium to soft | `#9B59B6` violet |
| 5 | **rift fill** | continental `age < 1` (C-3's rift-soft class) | soft | `#2EB8A8` teal |
| 6 | **loose recent deposits** (alluvium, lake beds, slope feet): the author's « meuble » | valley floors, `lake_map` margins and drained floors, low slope with A ≥ 1 km², foot of steep slopes | **soft** (a fall does not survive) | `#E2C896` sand |
| 7 | **evaporites** | endorheic lake footprints and their arid margins | soft | `#F2F0E6` white |
| 8 | **undifferentiated continental basement** | everything else on land (DECISION: the default, for want of sedimentary basins) | **hard** | `#8C8C8C` grey |

**The hardness**: ANCHORED in its class contrast, PROXY in its values.
- Stock & Montgomery 1999 (C-3's source) measure crystalline rock 1–5 orders harder than mudstone and volcaniclastic
  rock. Hard / medium / soft is that contrast.
- The values the gorge would use are for the gorge's work, which is paused.

**The priority at an overlap** (DECISION): 7 > 6 > 4 > 3 > 5 > 2 > 1 > 8. A surface deposit covers the rock under it,
and an edifice covers its basement.

**The land shares** (MEASURED, F151-S, with that measurement's thresholds, not round 3's):
- default basement 56.7 %;
- belt 14.3 %;
- loose deposits 11.6 %;
- craton 7.7 %;
- rift fill 4.8 %;
- volcanic 4.0 %;
- evaporites 0.95 %.

### 1.1 The fossil belts (prototype MEASURED in F151-F)

- **The recording**: during the C1 loop, through `run_with_closures_observed` (read-only; the final state is
  bit-identical, MEASURED; +0.02 s).
  - Per coarse cell, the steps convergent / collision / subduction upper / divergent, and the last step of each.
  - The sutures: at each merge, the previous step's touching cells of the two plates.
- **The rule (PROXY)**: a cell is **fossil belt** if it was in collision for ≥ 10 steps (≈ 7 Ma at ≈ 0.67 Ma per step,
  PROXY), OR it is a suture cell on continental crust, AND it is not convergent at the end (else it is an **active**
  belt).
  - On the témoin: **106 coarse cells**, all continental.
  - The 45 collision cells are all among the suture cells: the collisions ended in the merge at step 49.
  - The 10-step threshold does not bite here: every collision lasted 48 steps.
- **The active belts**: the convergent cells at the end, plus 2 coarse cells (12.5 km, PROXY) on the overriding /
  continental side.

### 1.2 The contacts: never on the 64² grid

- A coarse source (craton, belt, rift) is up-sampled **bilinearly** to HD (C-3's `upscale_k_to_hd` sampler), smoothed
  by a Gaussian of **σ = 0.5 coarse cell ≈ 3 km** (PROXY), and cut at the **0.5 iso-line**. This gives a curved,
  grid-free contact.
- **Snapping to the terrain** (DECISION, optional in v1): within 2 km of the iso-line, a contact may follow the nearest
  divide or valley floor. It is measured in round 3 before it is kept.
- **Volcanic rock follows the relief, not a disc.** The class covers the cells within 1.5 basal radii of an edifice
  whose altitude stands above the edifice's base level by > 10 % of its height (PROXY). The base level is the median
  altitude of the ring at 1.3–1.6 basal radii. The eroded cone's flanks and valleys then shape the contact (G-circles).
- **Loose deposits and evaporites** come from HD fields (49 m) and need no smoothing; thresholds only (PROXY, round 3):
  - floodplain: slope < 0.02 and A ≥ 1 km²;
  - slope foot: slope < 0.05 within 150 m below a slope > 0.3;
  - lake margin: within 100 m of a lake footprint.

## 2. The rules file (`geology_rules.toml`)

**Declarative, read by Ymir** (TOML: comments allowed, editable by hand; a new dependency, `toml`, DECISION). The v0
proposal is `geology_rules_v0.toml`, beside this spec.

**The vocabulary of conditions**, every field computed per HD cell by Ymir:

| condition | meaning |
|---|---|
| `rock = ["craton", "belt", …]` | the rock class is one of these |
| `belt = "fossil" \| "active"` | the belt record (§1.1) |
| `arc_within_km = d` | within `d` km of an active-margin (subduction upper plate) cell |
| `volcano_within_km = d`, `volcano_active = true`, `volcano_setting = "arc" \| "hotspot" \| "rift"` | edifices |
| `wetland = true` | the HD wetland mask |
| `slope_min`, `slope_max` | m/m |
| `valley_floor = true` | the construction's floor |
| `lake = true`, `endorheic_lake_within_km = d` | lakes |
| `river_min_area_km2 = a` | flow accumulation |
| `precip_min_mm`, `precip_max_mm`, `temp_min_c`, `temp_max_c` | climate |
| `coast_within_km = d` | distance to the sea |

**The semantics** (DECISION):
- a rule's conditions all hold (AND), and the rule gives `chance`;
- a resource's chance at a cell is the **max** of its rules (no sum: two settings do not make a richer deposit);
- `chance` is 0–100, an integer.

**The placers** (gold, tin, gems): a `[resource.placer]` block. From the cells where the resource's own non-placer
chance is ≥ `source_min_chance`:
- the chance travels **downstream along the D8 rivers** (cells with A ≥ `river_min_area_km2`);
- it starts at `start_chance` and loses `decay_per_km` per km of flow path;
- it is set only where the slope is < `max_slope` (gravels settle where the gradient drops).
- A placer is a line along the drainage, never a disc (G-placers).
- **The decay is PROXY**: real placers concentrate within a few km to tens of km of their source (Evans 1993, to
  check).

**The header**: `format`, `version` (semver), and a free `notes`. Ymir records the file's SHA-256 in the export
(§3).

## 3. The export (a new part of the container; nothing added to the existing layers)

| layer | dtype | resolution | raw size | PNG size (MEASURED, F151-S) | reason |
|---|---|---|---|---|---|
| `geologie_roches` | u8 class id | **full** (49 m, as `biome`) | 67 MB | 0.52 MB | the falls at sills and the rock view need the cell; it matches the biome layer |
| `chances_<resource>` × 12 | u8 0–100 | **¼ (195 m)**, DECISION | 4.2 MB each, 50 MB for 12 | 0.04–0.27 MB each, **1.4 MB for 12** | a hex covers ~1.75 cells: one chance cell is ~9 hexes, finer than any zone. The chances are a zoning, not a placement. A placer along a river becomes a 195 m strip, still on its river |

- **The PROXY grids measured**: the full resolution would cost 805 MB raw or 15 MB in PNG for the 12; ½, 201 MB raw or
  3.5 MB in PNG. The real grids' entropy may differ; round 3 re-measures.
- **The encoding** (DECISION, for the author): the container's rasters are raw today.
  - The chances can stay raw (50 MB at ¼).
  - Or the container gains a per-layer `encoding = "png"` (1.4 MB), a small writer and reader change.
  - This spec proposes **raw at ¼**: one convention, no new codec. PNG is the fallback if the size matters.
- **The header**, in `manifest.json` (the per-layer metadata):
  - `rules_format`, `rules_version`, `rules_sha256`;
  - the resource list with each resource's colour;
  - the rock classes with their colours.
- **Downsampling**: the max of each 4 × 4 block. A chance never vanishes by averaging.
- **The documentation**: `docs/geology_format.md` (round 3).

## 4. The viz views

- **« Roches »**: a categorical layer like « Biomes », with the colours of §1 and a legend (class, hardness).
- **« Chances »**: one colour per resource (from the rules file), **opacity = chance %** over the shaded relief.
  - **One resource**: a picker.
  - **All**: at a cell, the **highest chance wins** (ties: the rules file's order). The opacity is that chance. The
    legend lists the resources with their colours.
- The cell panel adds « Roche » and the chances at the cell.

## 5. The gates of round 3 (written now)

| gate | instrument | pass |
|---|---|---|
| **G-field** | the guard 6/6 (field and lakes) and the C2 /10 col hash `a8d2d538d692c2f0` | bit-identical (the geology is read-only) |
| **G-blocks** | the land contact cells (class change between 4-neighbours); the share of them within 1 HD cell of a coarse grid line (every 128 HD cells, both axes) | ≤ **1.5 ×** the share a random placement gives (2/128 per axis ≈ 3.1 %) (DECISION) |
| **G-circles** | each volcanic patch's isoperimetric ratio 4πA/P² against that of a rasterised disc of the same area | the ratio of the two ≤ **0.8** for each patch (DECISION) |
| **G-nonempty** | for each resource, the land cells with chance ≥ 10 on the témoin | > 0, or its absence explained. Tin and gems hang on the belts: the témoin has a fossil belt (F151-F) |
| **G-placers** | the cells whose chance comes from a placer rule | 100 % on the river network (A ≥ the rule's area), and each downstream of a source |
| **G-rules** | change one rule's chance and rerun the zoning | only that resource's grid changes (the others bit-identical), and the field hash is unchanged |
| **cost** | seconds per world for the rock grid + the zoning; seconds for a re-zoning alone | reported against 249.8 s (no threshold) |
