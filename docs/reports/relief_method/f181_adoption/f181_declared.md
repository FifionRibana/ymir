# Finding 181 — declared before any measurement

**The order of work (the rules):**
1. Part B (the intermediate and old ages) and Part C (the « équerres ») are measured first, with F180's flags as committed
   (`Margin::seaward_oo_collision_f180()`).
2. Part A (the adoption) becomes effective only if Part B holds.
3. **If Part B changes the cordilleras' aspect (width, disappearance, plateaus), the plates go to the author before the
   adoption is made effective** (the author's complement). No closure is modified, and no constant is calibrated.

## The author's words (condition v; to be copied verbatim into the ADR)

> « Non rien ne semble me gêner pour l'instant. Les graines montrent de la variabilité. »

They concern F180's diagnostic maps (`formations_15_F180.png`, `origin_15_F180.png`, `relief_active_15_F180.png`,
`transect_F180_seed7.png`). « Pour l'instant »: the validation holds for the young profile measured at F180.

## Part B — the definitions (F180's diagnostic, plus)

- The ages: `AGES[2]` (intermediate, 100 Myr) and `AGES[3]` (old, 300 Myr), with the physics level's settings unchanged
  (K2 rock, the torus, N_PHYS steps). F176 is the reference and F180 the candidate; the 15 seeds of F176's selection.
- **The formations and condition (v)**: F177's complement definitions (`formations178`). An Andean cordillera is an
  8-connected group of Andean-origin mountain cells touching an A1-active coast cell (« any » rule), as at F180.
- **A cordillera's width perpendicular to the coast** = its area / its length along the coast (the touched active coast cells
  × 6.25 km), in world km, with the terrestrial equivalent (× 7.5). This is parameter-free: the mean extent across the coast.
- **W₁₀**: global, active and passive (F178's `wall_index`, `w10_by_type` on the A1 map).
- **The macro**: the Kapos classes, the median, the land.
- **The coast's dimension**: `coast_torus` at 256² (box sizes 2–64 cells, its own definition), with the block lattice at 4 cells
  (one C1 cell).
- **The circumnavigation at the age**: `circumnavigation` on the 256² land at that age (r* in 256² cells). C1's end land does not
  depend on the age; F180's tier 1 gave 1 lost (seed 23) there.
- **Tier 2 (2 048², old)**: the unclipped transect of seed 7 (cells chosen as at F180), and the whole map of the seed with the
  longest cordillera at the old age (seed 7 excluded if another ties; if seed 7 is the longest, its map).

## Part C — the « équerres » (diagnostic, no correction)

- **The offshore active coasts**: A1's majority-active coast cells whose C1 cell and all 8 neighbours are sea at the end of C1
  (the Airy land, F177's B4 definition). The « équerres » are their 8-connected groups.
- **Reported**: their number on the 15 seeds; their share of the majority-active coast cells; and the origin of their
  continental type:
  - continental at t₀ (the initial continent or the band's continental part);
  - re-typed during the run (oceanic at t₀, continental at the end);
  - in the seaward band at t₀;
  - their s at the end.
- **The candidate, from the code before measuring**: `subduction/source_term.rs:266–268`. The continent-ocean floor-triggered
  reassignment gives a consumed oceanic cell below s = 0.05, with ≥ 2 convergent continental neighbours, to the continental
  plate **and re-types it continental**: a continental-typed cell with almost no crust, at sea. Its own comment mentions
  « grid-aligned false-land lines ». The measure says whether it is the source.

## The budget (announced; < 45 min wall per tier; resumable per world)

| tier | content | estimate |
|---|---|---|
| 1 | 256², 15 seeds × F176 / F180 × intermediate / old (60 physics runs, 30 C1 runs) | ~30 min |
| C | C1, 15 seeds, F180 (the équerres) | ~3 min |
| 2 | 2 048², old: seed 7's transect and one whole map | ~5 min |

## The adoption rule (written before)

The adoption becomes effective only if, at the intermediate and the old ages:
- (v) holds: ≥ 5 of 15 seeds with a cordillera, the separated highland never the only formation;
- there is no rampart: W₁₀ < 0.43 on every panel seed;
- the circumnavigation holds: lost on fewer than 2 seeds.

Otherwise: stop, report, the flags stay off. Stops: any closure change or solver; a guard state that changes.
