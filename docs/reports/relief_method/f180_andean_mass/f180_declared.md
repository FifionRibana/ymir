# Finding 180 — declared before any measurement

**Two changes to one chain, Davis-Suppe's source, behind the continent profile's flag: (A1) the Andean seed only behind the
coasts A1 calls active by the majority rule; (A2) a per-seed mass bound on Davis-Suppe's increments. Each is measured alone at
C1. No solver, no iteration; production unchanged.**

## Part 0 — the defect measured first (C1, the 15 seeds of F176's selection, F179 as committed)

- **The Andean coverage** of a world = the share of its coast cells covered by the Andean seed. The coast cells are A1's coast
  cells (continental type with an oceanic 4-neighbour) of the continents (4-connected, ≥ 20 cells). A coast cell is covered if
  it was seeded at any step of the run, or if a cell whose dominant Davis-Suppe family over the run is Andean lies at a periodic
  Euclidean distance < 3 cells.
- It is compared seed by seed with the A1 active share, both rules (« any », « majority »), at the end of C1.
- R0 (the reviewer's) and my prediction are in `f180_predictions.md`. No code of the closures changes for Part 0 (bench only).

## A1 — the Andean seed behind the active coasts only

- **Which coasts**: every coast cell that A1's walk classifies active **by the majority rule**: at the walk's first depth,
  strictly more active than passive boundary cells (an active boundary cell = an overriding trench cell, `overriding_trench_cells`,
  the rule of F178).
  - A coast cell on a plate boundary (A1's 3 × 3 case) has no walk; it keeps the existing boundary seeds, unchanged.
  - Passive coasts, and mixed coasts without a majority, get no Andean seed.
- **The walk** is A1's own: from the coast cell through P's oceanic cells (4-connected, periodic, breadth-first, the neighbours
  in the order E, W, N, S) to the first depth with boundary cells. **The walk's trench** is the first active boundary cell of
  that depth in the walk's order. **Its path** is the breadth-first predecessor chain from that trench cell back to the coast
  cell; d_front = the path's length (the coast's depth from the trench).
- **The seed** is P's continental cell on this path at depth d = max(d_arc, d_front) from the trench, if d ≤ D_gap:
  - d_front ≥ d_arc → the seed is the coast cell itself;
  - d_front < d_arc → the path is continued straight inland past the coast, in the direction of its last step (from the path's
    oceanic cell to the coast cell), for d_arc − d_front cells; the seed is the cell reached, if it is a continental cell of P
    (otherwise no seed for that coast).
- d_arc = 4 (F179: 180 km, Syracuse & Abers 2006 p. 10, ÷ 7.5 ÷ 6.25 = 3.84 → 4); D_gap = 10 (F178: 470 km).
- The seed is skipped by Davis-Suppe like a boundary seed (F179's skip mask). No reading of the sea level.
- It replaces F179's arc seed when the flag is set (`C1Ext::andean_active`).

## A2 — the per-seed mass bound (Dahlen 1990, p. 59, eq. 1)

- **The source**: Dahlen 1990, p. 59: « The mass flux per unit length along strike into the toe of the wedge is ρhV », and eq. 1,
  ∂/∂t[½ ρ W² tan(α+β)] = ρhV — the wedge's cross-section grows by exactly hV per unit length per unit time, whatever the
  erosion (Willett et al. 2001 p. 456: the tectonics is kinematic, independent of the erosion; p. 459: U = V_c h / L). No L is
  used: the bound is on the whole wedge of a seed.
- **The attachment**: every Davis-Suppe cell is attached to its nearest seed (Davis-Suppe's metric: 8-connected, √2 diagonals,
  within the plate, up to its reach); ties go to the seed reached first (the lower index).
- **Per seed σ and per step**:
  - the demand R_σ = Σ of the Davis-Suppe increments of σ's cells (the step computed as now);
  - the supply S_σ = (v_rel · n̂) · dt · 64 · s_in, the section hV of one cell length of trench, in s × cells.
- **If R_σ > S_σ**, every increment of σ is multiplied by S_σ / R_σ, so the wedge's shape is kept. There is no iteration; a seed
  with R_σ ≤ S_σ keeps its increments bit for bit.
- **The sharing**: an incoming (boundary) cell or a trench cell that feeds k seeds gives each S / k.
- **s_in, v_rel · n̂**:
  - a boundary seed: across its convergent edge of largest v_rel · n̂, the other plate's cell (s_in = its s);
  - an Andean seed: the sinking partner of its walk's trench cell (option 5's partner, the largest v_rel · n̂ of the subducting
    pair), s_in its s, v_rel · n̂ that contact's. For F179's arc seeds (the profile « F179 + A2 »), the walk is replaced by the
    nearest overriding trench cell of the same plate through P's cells.
  - **s_in is the whole thickness, an upper bound** (C1 has no sediment layer; the overriding plate's shortening is not sourced
    here).
- The bound applies to **every** Davis-Suppe family (collision included); its effect on the collisions is measured.
- **The mass record**: C1's step stats carry, per step, the seeds, the bounded seeds, the mass removed per cell, and the
  largest R_σ / S_σ after the bound.

## A3 — the controls (permanent)

- Both flags off: F179 bit for bit.
- Production's hash (`f175_c1_hash`).
- An infinite supply (`ds_mass_bound = ∞`): the result without the bound, bit for bit.
- The mass control: on a test world, after the bound, Σ increments ≤ S_σ (1 + 10⁻⁹) for every seed at every step.

## A4 — the viz

Two boxes: « Graine andine derrière les côtes actives (C1) » and « Borne de masse de Davis-Suppe (C1) ».

## Part B — the measures, by tiers

- **Tier 1 (C1, the 15 seeds)**: F176, F179, F179 + A1, F179 + A2, **F180 = F179 + A1 + A2**:
  - the Andean coverage against A1;
  - the land, the circumnavigation, r*;
  - the land added over F176 (in the band at t₀, the type at t₀);
  - the uplift by origin by distance band (F178's definitions);
  - s_src by family (the dominant family's cells), against F176;
  - the share of bounded seed-steps and the mass removed, by family.
- **Tier 2 (256², young, the 15 seeds)**, F176 and F180: F177's B3 maps and the formations (the complement's definitions),
  condition (v), W₁₀ (global, active, passive), the macro.
- **Tier 3 (2 048², young)**: the unclipped transects of seed 18 and of the seed whose longest cordillera is new; a whole map
  of a seed with a cordillera.

## The budget (announced; < 45 min wall per tier; resumable per world)

| tier | content | estimate |
|---|---|---|
| 0 | C1, F179, the 15 seeds (coverage) | ~5 min |
| 1 | C1, 15 seeds × 5 profiles | ~15 min |
| 2 | 256², 15 seeds × F176, F180 | ~20 min |
| 3 | 2 048², two transects and a map | ~5 min |

## The rules (written before)

- **Stops**:
  - the rampart (W₁₀ young ≥ 0.43 on a panel seed, or an Andean coverage above 70 % on ≥ 3 seeds);
  - the selection fails, or the circumnavigation is lost, on ≥ 2 seeds;
  - the mass control fails;
  - any other closure change, a solver or an iteration;
  - any value changed on the way (d_arc, D_gap, s_in).
- **Adoption** only if all hold: (i) the A1 median in 10–50 % with both rules; (ii) ≥ 70 % of the continents with an active
  margin; (iii) condition (v); (iv) no rampart; (v) the author validates B3's maps, the formations' plate and the transects.
- Forbidden: as in F165–F179.
