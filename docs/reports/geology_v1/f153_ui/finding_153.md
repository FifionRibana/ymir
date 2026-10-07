# Finding 153 — geology v1: the interface finishes the author asked for, and the work closed

**Status: built (viz only, plus one core function for the inspector); nothing committed.** The world, the exports and
the guard are unchanged.

## The author's criterion, recorded as given

- « À mon avis, ça semble très bien déjà. Avant de clore, j'aimerais faire des updates UI. »
- « La roche devrait avoir une étape en haut, vu que c'est une vue au même titre que les biomes. Ce n'est pas un
  masque qui s'affiche sur les autres vues. Ou sinon il faut que ce soit un masque. »
- « Les menus de la toolbar ne bloquent pas le clic : lorsqu'on clique sur un élément du menu pour l'activer, le menu
  se ferme. Ce n'est pas simplement activable. »
- « Il faut que l'inspecteur affiche plus d'informations : une section Roches au même titre que les biomes, ainsi que
  les ressources et leur favorabilité sur le pixel survolé. »
- This round touches the viz only: the world, the exports and the guard do not change.

## Part 0

- F152 committed (`f85f431`), with its ADR entry completed ("Recorded at commit (F153)"):
  - G-blocks held only after the declared amendment;
  - no districts is a known v1 limit;
  - the volcanic favourability discs come from the distance rules, not the code.
- **Checks after** (`checks_after.txt`):
  - `cargo check --workspace` clean;
  - lib 617 passed (+1, the inspector = export test);
  - **the guard 6 / 6, field and lakes = banc** (C2 /10 col `a8d2d538d692c2f0`).
  - The viz tests: 31 passed, **after one fix**. The first run failed `c1_worker_runs_hd_chain_and_ships_product`,
    which counted 6 worker phases; the new `Geology` phase makes 7. I had updated its expected ORDER but not its two
    counts.

## P — « Roches » is a step, « Chances » a mask

- **The frieze gains an eighth step, « Roches »**, after « Biomes ».
  - A worker phase `HdPhase::Geology` reports the geology stage's start and end.
  - The step shows the categorical rock view and its legend (class, hardness).
- **« Chances » is no longer a view but a MASK**: a « Chances » box in the toolbar, like « Réseau » or « Rivières LL ».
  - It draws over any step.
  - Its menu « ⛏ ressources ▾ » holds: the box again; the resource (one, or all with the strongest winning, as radio
    buttons, each in its colour); the rules file's path; « Recharger les règles » with the re-zoning time; a refused
    file's error.
  - Opacity = favourability.
  - The mask has its own legend, top left of the map (the resources shown, their maxima, "grille ¼ (195 m)").

## M — the toolbar menus stay open

- **The four toolbar menus** (rivers' filters, « ⛏ ressources », the relief view, « Tectonique ») open through
  `tool_menu`: egui's `MenuButton` with `PopupCloseBehavior::CloseOnClickOutside`.
  - A click on an item (a box, a slider, a radio button) uses it **without closing the menu**.
  - The menu closes on a click outside it, or on Escape (egui's popups close on Escape).
- The resource picker was a combo box, a second popup whose clicks would count as "outside" and close the menu. It is
  now radio buttons inside the menu.
- **No click on a menu reaches the map**: the menus are egui popups on a higher layer than the map canvas
  (`ui.interact`), and egui gives a click to the topmost layer only.
- **Not verified interactively** (no GUI driving here): the author's look confirms it.

## I — the inspector

- **« ROCHES »** (beside « BIOME »): the hovered cell's class, with its colour and its hardness (dure / moyenne /
  tendre).
- **« RESSOURCES »**: every resource of non-zero favourability at the cell, strongest first:
  - its favourability (0–100);
  - **the rule that gave it**, its `source` from the file, numbered in the file's order;
  - the structural modulation factor « × g » when the rule is modulated, or « placer : source » when a placer raised
    it.
  - The section says the favourability is read on the ¼ grid (195 m), as exported.
- **The values are exactly the export's**:
  - `geology::zoning::explain` READS the favourability from the zoning grid that the export writes;
  - only the attribution is recomputed, by the same evaluation as the zoning (shared `best_rule`). A tie goes to the
    rule, as the MAX keeps it.
  - **A permanent test**, `the_inspector_shows_the_exports_values`, decodes the exported PNGs and compares them with the
    inspector's values at 7 cells: a modulated rule, a plain rule, a placer source, a placer cell, a tie, an empty
    cell. It also checks each attribution.
- The inspector reads the ACTIVE zoning: after « Recharger les règles », the reloaded rules.

## Predictions

- **The reviewer's (hypotheses to check)**:
  - "the guard stays 6/6, the hash identical": HELD;
  - "the inspector = export test passes on every tested cell": HELD (7 cells × 2 resources).
- **Mine**:
  - P-guard: HELD;
  - P-inspector = export: HELD.
  - **P-attribution**: "ties go to the first rule". **REFUTED in its detail**: my first version of the test expected a
    placer at a cell where the placer (30) and the modulated rule (35 × 0.85 = 30) tie, and the tie went to the rule.
    The test was corrected to name the tie, and to test the placer at a cell where it wins (40 > 28). The témoin's
    count of ties was not measured.
  - P-menus: built as predicted; not verified interactively.
- **Meta**: HELD.

---

## C — the geology v1, closed: a summary for the author

**What Ymir delivers** (`.ymir` container 1.1.0, `docs/geology_format.md`):
- **`geologie_roches.png`**: one rock class per cell (49 m), fixed for a world. The classes are craton, belt
  (metamorphic, fossil or active), basaltic volcanic, arc volcanic, rift fill, recent loose deposits, evaporites, and
  undifferentiated basement. Each class carries its hardness, for the falls at the lake sills.
- **`favorabilite_<resource>.png`**: a **relative favourability** (0–100) per resource on the ¼ grid (195 m, ≈ 9
  hexes). It is not a probability and not a deposit. Covered: iron, gold, silver, copper, tin, stone, clay, sand and
  gravel, salt, peat, sulphur and obsidian, gems.
- **`geologie.json`**: the legend and the SHA-256 of the rules file the favourabilities came from.
- The geology stage costs ≈ 3 s per world; a re-zoning costs 0.05 s.

**How to tune the rules:**
1. Copy `crates/ymir-core/data/geology_rules_v0.toml` to `geology_rules.toml` in the viz's working directory.
2. Edit it: a rule's conditions all hold, and it gives its `chance`; the strongest rule wins; placers run downstream.
   The format is in `docs/geology_format.md`.
3. In the viz: « ⛏ ressources ▾ » → « Recharger les règles ». The map and the inspector update in 0.05 s. A file that
   fails to parse is refused with its error.

**The first tunings suggested:**
- **lower stone**: at 40–60 nearly everywhere, it hides the others in the « toutes » view;
- **widen peat**: 165 km² on the témoin;
- **replace the volcanic discs** (sulphur and obsidian's `volcano_within_km`) with the volcanic ROCK class
  (`rock = ["basaltic_volcanic", "arc_volcanic"]`), which follows the cone's relief.

**The known limits:**
- **no districts inside the belts**: Ymir has no structure finer than a belt. A finer structural proxy could lift it
  later;
- **no sedimentary basins**: erosion deposits nothing, so there is no coal, no limestone and no sandstone;
- the fossil-belt rule is a PROXY: its duration threshold is untested on the témoin, where every collision lasted 48
  steps.

**What Living Landz does on its side:**
- a **rarity** per resource (Ymir gives the favourability, never the rarity);
- **the distribution** of the deposits according to the favourability;
- **the draw at each reset of the mist** (the resources are not world state; the rock is).

**The next work, in the agreed order: lakes, gorges and falls resumed.** The rock hardness at the lake sills is now
available (`geologie_roches`, the `hardness` of each class).
