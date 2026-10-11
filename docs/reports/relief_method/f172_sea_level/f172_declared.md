# F172 — declared before any measurement (2026-10-11)

**The sea-level rise, and a screening of the seeds.**

**Production is unchanged**; the guard 6 / 6 is not touched.

**Frozen**: T, Airy, the sources' drive, the torus, K2 (ρ_K 10, d_c 1.5 km, σ_lnK 0.5), 16 plates with the 20 / 24 fallback,
F164's cascade.

**Not blind**: F171's readings and the code. No F172 variant has been run.

**Part 0**: F171 committed (`3c0b653`) with (a) K2 as the continent profile's rock, (b) the change of criterion (mechanisms by
their effect, the macro layout by a distribution of seeds, A the author's slider with a provisional default, ρ_K 10), (c) the
coast against Corsica in absolute value, (d) the walls, facets and rampart measured here. The author's view and position are
recorded verbatim.

## The two sea constants (read in the sources)

| constant | value | source |
|---|---|---|
| **Δ** (the lowstand) | **134 m** | Lambeck, K., Rouby, H., Purcell, A., Sun, Y., Sambridge, M. (2014). Sea level and global ice volumes from the Last Glacial Maximum to the Holocene. PNAS 111(43), 15296–15303. doi:10.1073/pnas.1411762111 |
| **T_q** (the lowstand's duration) | **2.58 Myr**, the base of the Quaternary | ICS, the GSSP list: the Gelasian GSSP (Monte San Nicola, Sicily), 2.58 Ma, ratified in 2009 as the base of the Pleistocene and the Quaternary |

- **Lambeck et al. 2014**, read on PubMed Central (PMC4217469); the PNAS site refused the fetch, HTTP 403.
  - The abstract: « a slow fall to −134 m from 29 to 21 ka BP ».
  - The Discussion: the eustatic sea level « reaches its lowest value of ∼134 m », « greater than the frequently cited −125 m ».
  - The brief's « ~120 à 130 m » is the frequently cited figure. **The article prevails: Δ = 134 m.**
  - The printed page of each quote is not available in the scraped text (abstract and Discussion sections).
- **T_q**: stratigraphy.org/gssps, read 2026-10-11.
- **The sensitivity**: Δ = 60 m, a more realistic mean Quaternary level than a permanent lowstand. The témoin, one age (old),
  measure only.

## A — the diagnostics: the rampart, the walls, the facets (F171's worlds, K2, no sea-level rise)

### A1 — the perimeter (64²)

- **The perimeter**: the end-of-C1 land cells (h_src > 0) 4-adjacent to the sea.
- **Each perimeter cell's boundary**: the nearest class within Chebyshev 1, in priority convergent > divergent > transform >
  none (internal). It comes from C1's final `classify_boundaries`, retargeted.
- Reported: the four shares, and Σ Δs_DS against the distance to the coast (64² land, bins of 1 cell, 0–7 and ≥ 8).

### A2 — where the mountains are, and the rampart index (declared)

- On the 256² physics level, the periodic EDT of the land to the sea.
  - f₁₀ and f₂₀: the mountain class's share of the land within 10 and 20 km of the coast.
  - f_in: its share of the land farther than 20 km.
- **The rampart index R = f₂₀ / max(f_in, 0.05).**
- **A rampart** = f₂₀ ≥ 0.40 and R ≥ 2: the coastal band is at least 40 % mountain, and twice the interior.
- Reported at B0, young, intermediate and old.

### A3 — where the walls and the facets are

- At 2 048², in the mountain class (the 256² classes enlarged):
  - **the walls** = the cells whose slope (`physio::slope_deg`) exceeds 28°;
  - **the facets** = the battery's `read_cell` facet measure on the sub-mask.
- Each is reported by the distance to the coast (0–10, 10–20, > 20 km) and by the rock at the surface (cover, basement, shield).
- At B0 and old, the 4 seeds. Measured only.

## B — the sea-level rise

### B1 — the mechanism

- **The physics level**: during the last T_q of the history (the steps whose remaining time, the sum of the later steps'
  durations, is ≤ T_q; under F169's age mapping), the base level is −Δ.
  - Every cell ≤ −Δ is the sea for the solver, the talus and the diffusion. This is the history level's existing rule, so an
    inland cell below −Δ is a base level too; the count is reported.
- **The amplification levels**: the whole chain runs on the field raised by Δ (so its sea at 0 is the real −Δ: the ocean, D5,
  the land floor), and the result is lowered by Δ at the end.
- **At the very end**, the sea is back at 0 m.
  - The final ocean is the largest periodic sea component of the cells ≤ 0 (`ocean_mask_t`): every connected cell under 0 is
    the sea, which drowns the valleys.
  - The closed depressions under 0 not connected to it stay land and are handled by the existing hydrology (lakes).
- **No deposition**: the drowned valleys are not filled. A known limit.

### B3 — the negative control (permanent)

Δ = 0 is F171 (K2) bit for bit.

### B4 — the shelf

- **The width of the 0 to −Δ zone**: for each coastal sea cell (a final sea cell 4-adjacent to land), the distance (EDT, 2 048²)
  to the nearest cell ≤ −Δ. Reported: the median and p90 (km).
- **The mean depth against the distance from the coast**, 0–20 km by 1 km, plotted.

### B5 — viz

In the Cascade window: a checkbox « Remontée de la mer » (the history level and the chain), and a view « Paysage de bas niveau /
après remontée » (the same field shown with the sea at −Δ or at 0).

## C — the effect (témoin, 42, 1, 9; K2; B0, young, intermediate, old; physics 256², the chain to 2 048², torus)

### The coast against Corsica (declared bands)

| measure | Corsica | « close to Corsica » |
|---|---|---|
| dimension | 1.12 | 1.06–1.18 |
| length ratio | 1.37 | 1.15–1.60 |
| aligned share | 0.27 | 0.22–0.32 |
| block index | 0.78 | 0.58–0.98 |

- **A coast « gets closer to Corsica »** if, after the rise, the dimension and the length ratio both move towards Corsica's
  values, and the number of measures inside the bands does not fall.
- Reported: before / after on every world, with « before » = F171 K2, the same world without the rise.

### The drowned valleys (2 048², declared detector)

- **D** = the final ocean cells with z > −Δ, i.e. the ocean the rise added. Its periodic 8-connected components are the
  candidates.
- **The length inland** = the BFS geodesic distance inside the component from its cells adjacent to the lowstand ocean
  (final-ocean cells ≤ −Δ).
- **A drowned valley** = a component of length ≥ 2 km.
- **Its river's size** = the largest D8 drainage area (the final field's `hydrology`) among the land cells adjacent to it.
- Reported:
  - the count, the lengths' p50 / p90 / max, Spearman(length, area);
  - the islands created: the land components ≥ 1 km² after the rise minus those of the lowstand land, the cells > −Δ;
  - the land lost against F171 K2: the share of F171's land cells that are the sea after the rise.

### The rest

- Macro (the classes must move < 3 points against F171 K2), A's diagnostic, the lakes, the texture and the regularity as in
  F171.
- **Images**: the whole maps before and after; the 25 km coast crops Corsica | F171 | F172; the shelf profile; the A2 map.
- The cost.

## D — the screening (physics 256² only, no chain)

- **The seeds**: 1 to 32, and the témoin: 33 seeds.
- **The profile**: 16 plates with the fallback, K2, the sea-level rise. Two ages: young (20 Myr) and old (300 Myr).
- **The descriptors**, per seed and age:
  - the N retained, or none;
  - the land % and r*;
  - the classes, the median and the peak;
  - R and the rampart verdict;
  - the perimeter's convergent share;
  - the mountain regions > 100 km² (8-connected, 256²);
  - the cover / basement / shield shares;
  - the lakes; the cost.
- **The outputs**:
  - the catalogue (`catalogue.tsv` and a plate of 256² thumbnails, shaded and classes; one row per seed, columns young | old);
  - each descriptor's median and deciles;
  - **the share of seeds within Europe's bands** (mountains 16.5–37.2 %, plain + plateau 37.5–84.5 %, median 98–390 m, all three)
    for at least one age;
  - the share with a rampart (the A2 threshold) at young;
  - Spearman(the perimeter's convergent share, R) at young and at old.
- **The cost**: the total, per seed, and the selection's share.

## R — the stop rules (the brief's)

- **The rise's success** (a mechanism): on ≥ 3 of C's 4 seeds and for ≥ 2 ages, the coast gets closer to Corsica (the
  declared sense), and the classes move < 3 points.
- **No effect** (the dimension rises < 0.01 everywhere): report the valleys' width at the mouth and the coastal profile, and stop
  before the screening.
- **The screening has no success rule.** A rampart on > 80 % of the seeds at young is flagged as a probable C1 bias.
- **Forbidden**: as F165–F171.
