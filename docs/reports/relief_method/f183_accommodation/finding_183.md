# Finding 183 — the sediments' accommodation (the load, the compaction, the recycling by subduction): built behind their flags; the combination keeps the sea on the panel, but the declared stop fires at tier 2 (the circumnavigation lost on 6 of 14 seeds at the intermediate age, the old age collapses), with a scheme defect of F182 at its root

**Status: Part A is built and its controls pass. Tier 0 and tier 1 ran. Tier 2 and tier 3 ran, then the brief's stop fired:
the combination A1 + A2 + A3 loses the circumnavigation on 6 of the 14 seeds that had it at the intermediate age, and on all
4 at the old age.** Nothing is committed. The reference profile is F181's default (`Margin::continent_default()`).

Declared before any measurement: `f183_declared.md` (Part 0, the resumption R1–R6). Predictions: `f183_predictions.md`.

**The author's decision (verbatim)**: « Oui ok pour ajouter les 3 points. »

## Part 0 — the sources

- **ρ_g = 2 720 kg/m³**: Sclater & Christie 1980, JGR 85(B7), **Table A1a, p. 3732** (the PDF's page 22, read on the scanned
  image: the text layer garbles the digits). The row « Shale » gives f₀ 0.63, c 0.51 × 10⁻⁵ cm⁻¹ (1 / c = 1 961 m) and
  ρ_sg 2.72 g/cm³. It is Yuan EPSL's silt exactly (φ₁,₀ 0.63, z₁,₀ 1 960 m). The row « Sand » (0.49, 3 704 m, 2.65) is Yuan's
  sand. The reviewer's memory was confirmed by the PDF.
- **ρ_w = 1 030 kg/m³**: Whitehead & Clift 2009 p. 7 §[24].
- **ρ_c = 2 700, ρ_m = 3 300 kg/m³**: the rebound's (F165, no page). Not harmonised with Airy's 2 800 here (queued).
- **The compaction**: Yuan EPSL eq. 7 (p. 3), Table 1 (p. 6); the one-class column T(H) = H − φ₀ z₀ (1 − e^(−H/z₀)) (Appendix A).
- **The flexure**: the rebound's (α = 64.4 km).
- **An order of magnitude for the deltas** (not a criterion): the Ogooué's post-Early-Pliocene sediments reach « ∼2 km at the
  location of the present-day shelf break » in 5.3 Myr (Yuan EPSL p. 4).

## Part A — what is built (behind `Deposition::{load, compaction, recycle}`, all off by default)

- **A1 the load** (`Load`, ρ_g 2 720, ρ_w 1 030):
  - The cover's solid S and thickness H are tracked per cell.
  - Its load is S (ρ_g − ρ_w) under water and S ρ_g + (H − S) ρ_w on land.
  - The step's net ΔL = Δ(the cover's load) − ρ_c e_b passes once through the rebound's filter, Δz = −G ∗ ΔL / ρ_m, at the
    step's end. It replaces the rebound when on.
  - The incision takes the cover first (at ρ_b, inside the cover's load change), then the basement (ρ_c), as the author
    specified.
- **A2 the compaction** (`Compaction`, φ₀ 0.63, z₀ 1 960 m):
  - every deposit enters at the surface porosity;
  - at the step's end every cover re-compacts, H = T⁻¹(S);
  - T⁻¹ is a table computed once (0.5 m steps to 30 km; its worst round trip is 7.3 × 10⁻⁵ m);
  - the eroded cover leaves its solid S − T(H − e_c);
  - the pass carries solids (a new deposit's thickness is solid / (1 − φ₀)).
- **A3 the recycling** (`recycle`, C1's `Consumption`):
  - C1's run now records its consumption per snapshot interval. This is an observation: C1's results are unchanged.
  - A cell is consumed when it was oceanic on plate P, joins a plate that existed before the step, and P survives.
  - On a plate with no continental cell, every step removes r = 1 − (1 − F_j)^share of the cover. All of it subducts (the
    declared upper bound).
- **A4 the amended mass control**: in = ΔStock (measured on the cover's solid) + the cover eroded + the cover lowered by the
  hillslope + subducted + exported, at 10⁻⁹. F182's flux control is kept.
- **A5 the controls** (permanent: `the_accommodation_controls`, `the_compaction_table_and_the_porous_pass`):
  - the three flags off through the new entry are F182 (the same field);
  - zero densities give no cover load;
  - z₀ = ∞ compacts nothing;
  - no consumption recycles nothing, bit for bit;
  - with the three on, both controls hold < 10⁻⁹, the recycling removes solid, and no solid is created.
- **F182 bit for bit, checked across the trees**: the same synthetic history at HEAD and at the commit c00817b gives the same
  field and rebound hash (84aad857019d8a19).

## Tier 0 — C1's consumption (15 seeds; C1 does not depend on the physics level's age)

- **The recyclable plates** (no continental cell) cover **8.4–28.3 % of the map (median 14.9 %)**. The plates carrying
  continental cells, with their bands and margins, cover the rest.
- **Their area-weighted consumed share over the history**: median 23.7 %, range 2.5–188 %. Above 100 %, a plate receives cells
  (C1's reassignment to the overriding plate) that are consumed again.
- So the recycling's capacity is bounded by the small recyclable area before C1's rate.

## Tier 1 — the panel, the five variants (256², 50 worlds)

| variant | young: the sea kept | circumnavigable | the sea at the end (median) | intermediate: the sea kept | circumnavigable | the sea at the end |
|---|---|---|---|---|---|---|
| F182 | 1 / 5 | 0 / 5 | 0 % | 0 / 5 | 0 / 5 | 0 % |
| + A1 | 5 / 5 | 0 / 5 | 37.4 % | 0 / 5 | 0 / 5 | 0 % |
| + A2 | 5 / 5 | 0 / 5 | 22.1 % | 0 / 5 | 0 / 5 | 0 % |
| + A3 | 2 / 5 | 0 / 5 | 0 % | 0 / 5 | 0 / 5 | 0 % |
| **A1 + A2 + A3** | **5 / 5** | **4 / 5** | **47.0 %** | **5 / 5** | **4 / 5** | **47.1 %** |

- **The baseline without the deposition**: the sea is 38–51 %, circumnavigable on 4 / 5 at both ages. Seed 23 is not
  circumnavigable without the deposition either: F181's « lost on 1 seed ». So on the panel the combination loses nothing.
- **Subducted / input** (the combination): median 5.4 % young, **9.1 % intermediate** (range 0.9–19.7 %).
- **The mass controls**: the flux control at worst 3.1 × 10⁻¹¹, the stock control at worst 6.7 × 10⁻¹¹, over the 50 worlds.
- **The capacity at t₀** (local Airy columns, an upper bound on the load's effect): eroded / capacity = 0.71–1.67 (F182's),
  0.34–0.81 (A2), 0.52–1.21 (A1), **0.09–0.22 (A1 + A2)**. The load and the compaction together multiply the ocean's solid
  capacity by ~8.
- **The thickest cover** (the combination): median 4.9 km young, 5.6 km intermediate (max 7.0 km), against the Ogooué's ~2 km
  in 5.3 Myr.
- **With A1 alone**, a one-cell spike reaches 47 km at z = 41 km (seed 29 young). It is the defect below, amplified.
- Images: `sea_panel_young_top_intermediate_bottom.png` (grey F182, blue A1, green A2, orange A3, red the combination) and
  `thickness_A123_young_top_intermediate_bottom.png`.

## Tier 2 — F182's tiers on the combination, and the stop

### The 15 seeds at the intermediate age (F182's tier 1)

| measure | without | with (A1 + A2 + A3) |
|---|---|---|
| circumnavigable | **14 / 15** | **8 / 15** (lost on seeds 1, 3, 6, 7, 10, 11) |
| large lakes (≥ 64 cells), median | 3 | 4 (per-seed change median +67 %) |
| dust (1–2 cells), median | 4 | 16 |
| lake area, median | 1 477 km² | 2 112 km² |
| shelf, passive / active (median) | 2.2 / 3.1 km | 5.3 / 10.0 km |
| shelf in Harris ÷ 7.5 (4.65–16.4 km) | — | 5 / 15 |
| shelf in Harris unconverted (34.9–123.0 km) | — | 1 / 15 |
| plains / mountain | 44.0 / 23.1 % | 15.7 / 51.3 % |
| median altitude | 115 m | 518 m (per-seed change median +396 m) |
| the coast's dimension / length ratio | 1.180 / 1.257 | 1.064 / 1.062 (\|Δ\| up to 0.163) |
| deltas (land gained at the coast) | — | 28 per seed, 8.1 % of the land |
| the sediment: land / lakes / sea | — | 37.7 / 1.0 / 61.3 % |
| (v) | — | 15 / 15 |
| W₁₀ \|Δ\| max | — | 0.066 |

### The panel, young and old (F182's tier 2)

- **Young**:
  - circumnavigable 4 / 5 → 4 / 5;
  - median altitude 101 → 277 m; plains 47.8 → 28.3 %;
  - shelf passive 2.2 → 6.2 km; (v) 5 / 5.
- **Old (300 Myr): the collapse.**
  - The land covers 94–100 % of the map on 4 of 5 seeds (seeds 19, 23, 29 entirely land).
  - The circumnavigation is lost on 4 / 4.
  - The median altitude is 53 → 5 142 m; the mountain share 15.6 → 99.9 %.

### What drives it: F182's marine scheme fixes the sea for a whole step

- The marine step's domain is the step's sea cells (`!land` after the uplift). A sea cell keeps receiving the river's source
  and the diffusion's inflow for the whole step, even once it stands above sea level. Yuan's eq. 5 holds only for h < h_sea.
- **The defect grows with the step's duration**:
  - at the young age, in the combination, 1–15 % of the marine deposit stands above 100 m;
  - at the intermediate age, measured on the two losing seeds and one keeper (`f183_diag_losers`):

    | seed | above 100 m | above 500 m | the sea over the history |
    |---|---|---|---|
    | 7 (lost) | 54 % | 9.5 % | 55 → 19 % |
    | 1 (lost) | 7.8 % | 3.0 % | 52 → 38 % |
    | 18 (kept) | 0.9 % | 0.5 % | 60 → 52 % |

  - at the old age the longest steps last 5.7 Myr.
- **Not all of it is the defect.** On seed 1 the sea settles at 38 % from 74 Myr on. Losing the circumnavigation there means a
  narrow seaway filling (F177: the straits are 1-cell rift strips). A shallow strait filling with sediment is physically
  plausible.
- **The defect was already in F182's tier-1 numbers.** The « thickest cover » of 8–23 km and the sea's fast loss there partly
  rest on it.

### The reference (Yuan's iterative coupling, the bench only, 2 panel seeds)

- **It never converged**: 50 iterates on 300 / 300 steps for both seeds. Yuan proves convergence for G / p̃ ≤ 1 only (JGR p. 8),
  and G = 1.55.
- Its numbers are not a converged implicit reference:
  - the deposit map's r = 0.22–0.34;
  - the total deposit +159–254 %;
  - the median altitude 369 → 2 267 m and 1 123 → 2 184 m.
- **P7** holds formally: the deposit gap is 1.09–1.10 at the longest step (1.47 Myr) against 0.26–0.29 at the shortest
  (18.7 kyr), with a Spearman (duration, gap) of 0.97–0.99. It is uninterpretable for the same reason.

### The sensitivity (the panel, intermediate, on the combination)

| variant | large lakes | shelf passive | plains | median altitude | the sediment at sea |
|---|---|---|---|---|---|
| G = 0 | 3 | 6.6 km | 50.8 % | 199 m | 98.5 % |
| G = 3.1 | 4 | 23.0 km | 7.9 % | 902 m | 48.8 % |
| K_M = 300 | 1 | 21.1 km | 23.3 % | 412 m | 61.1 % |
| K_M = 500 | 3 | 16.7 km | 25.3 % | 382 m | 62.0 % |

G moves the relief most: the fluvial deposit on land raises it. K_M changes little.

### Tier 3 (2 048², seed 23, intermediate)

- The land grows 102 224 → 112 087 km²; the median altitude 73 → 393 m.
- The shelf transect at a passive coast (the 256² width is 138 km, the break 3 973 m deep) ranges −5 022 … 307 m.
- Images: `map2048_off/on_seed23.png`, `crop2048_off/on_seed23.png`, `transect_shelf_seed23.png`.

## The predictions

- **The reviewer's**:
  - **Q1** (A1 alone keeps the sea on ≥ 3 / 10): held, 5 / 10 (young 5, intermediate 0).
  - **Q2** (the combination's circumnavigation on ≥ 9 / 10): refuted, 8 / 10. The two are seed 23, not circumnavigable without
    the deposition.
  - **Q3** (the recycling removes ≥ 20 % at the intermediate age): refuted, median 9.1 %, max 19.7 %.
  - **Q4** (> 5 km under the large deltas): the thickest cover is 4.4–7.0 km; its location against the deltas was not checked.
    Undecided.
  - **Q5** (the amended mass control holds): held, ≤ 6.7 × 10⁻¹¹.
- **F182's, on the combination (15 seeds, intermediate)**:
  - P1 (the large lakes −50 %): refuted (+67 %).
  - P2 (the shelf in Harris ÷ 7.5 on ≥ 10 / 15): refuted, 5 / 15; Harris unconverted 1 / 15.
  - P3 (the plains +5 points, the median +20 m): refuted (−29 points, +396 m).
  - P5 ((v) 15 / 15, W₁₀ < 0.03): (v) held, W₁₀ refuted (0.066).
  - P6 (the dimension < 0.03): refuted (up to 0.163).
  - P7: held formally on a non-converged reference.
- **Mine**:
  - N0: held narrowly (8 / 15 under 30 %).
  - N1 and N2: refuted (A1 alone and A2 alone keep the sea young on 5 / 5).
  - N3: held (< 20 %).
  - N4: refuted (the combination keeps the sea on 5 / 5 at both ages on the panel).
  - N5: held (4.4–7.0 km).
  - N6: refuted on its number (6.7 × 10⁻¹¹, not < 10⁻¹²), held on the brief's 10⁻⁹.
  - N7: held (25–43 s per world).

## The stop and the author's questions

**The rule (written before)**: « Si la combinaison perd encore la circumnavigation sur au moins 2 mondes : arrête-toi après le
palier 1 et rapporte ».
- **Read at tier 1** as a loss caused by the deposition against the profile without it, the rule did not fire: 8 / 8 kept on
  the panel, seed 23 lacking it from the start. This reading is mine; the author may read the rule literally, against all 10
  worlds, in which case it fired at tier 1.
- **At tier 2 it fires on any reading**: 6 of 14 seeds at the intermediate age, and 4 / 4 at the old age.

**For the author**:
1. **F182's marine domain is fixed for a step** (the root of the old-age collapse and of most of the losses). A closure-level
   correction: the domain re-read at each marine sub-step (h < h_sea). A cell at or above sea level takes no further source;
   its share goes to its lowest marine 4-neighbour, else stays (a stranded deposit). It is explicit, with no iteration, but it
   changes F182 (behind a flag, so F182 stays bit for bit).
2. **The design question the brief reserved**: the map's land / sea ratio (0.96–1.63 against the Earth's 0.41), or an export to
   an ocean off the map. It cannot be judged before 1, since the defect inflates the filling.
3. **The relief**: on the combination the land's median altitude rises by +185 (young) to +396 m (intermediate), driven by G
   (the sensitivity: G = 0 → 199 m, G = 3.1 → 902 m). This is a large change of the profile's hypsometry. Is it acceptable in
   principle, once 1 is corrected?
4. **The reference**: Yuan's iteration does not converge at G = 1.55 > 1. Keep it as a diagnostic, or drop it?

**The queue (the author's, F183)**: harmonise ρ_c between the rebound and Airy, on a source with a page; the erodibility of the
sediment against the rock.

**Not changed, waiting for the decision**: the code stays behind its flags (off by default, production unchanged). The tiers
are resumable in `tests/f183_cascade.rs`.

## The checks (after; the three flags off by default)

- `cargo check --workspace`: clean.
- 667 lib tests pass (F182's 665 plus the two new controls); 33 viz tests pass.
- The guard reads Match on all six states (field and lakes).
- Production's C1 hash is unchanged (`checks_f183.txt`).

**At commit** -- the author's words (verbatim): « Feu vert pour commiter le F183 tel quel (les trois drapeaux coupés par défaut,
les contrôles, le rapport avec l'arrêt). Ajoute à l'ADR trois erreurs du relecteur : 1. F182 : le schéma marin n'exigeait pas
que le domaine soit relu à chaque sous-pas (h < h_mer), alors que l'éq. 5 de Yuan ne vaut que sous la mer ; 2. F183, Q3 : 20 %
de recyclage prédits sans voir que les plaques recyclables ne couvrent que 8 à 28 % de la carte ; 3. F183 : la règle « perd la
circumnavigation » était ambiguë. Ta lecture (une perte due au dépôt par rapport au même monde sans dépôt) était la bonne. Elle
sera désormais écrite ainsi. La référence de Yuan est abandonnée : sa convergence n'est prouvée que pour G ≤ 1 (JGR p. 8), et
G = 1,55. Garde son code dans le banc, avec cette note. »
- **The reviewer's errors, recorded (the author's)**:
  1. **F182**: the marine scheme did not require the domain to be re-read at each sub-step (h < h_sea), although Yuan's eq. 5
     holds only under the sea.
  2. **F183, Q3**: 20 % of recycling was predicted without seeing that the recyclable plates cover only 8–28 % of the map.
  3. **F183**: the rule « perd la circumnavigation » was ambiguous. The reading applied at tier 1 was the right one: a loss
     caused by the deposition, against the same world without it. **From now on the rule is written that way.**
- **Yuan's reference is abandoned**: its convergence is proven only for G ≤ 1 (JGR p. 8), and G = 1.55. Its code stays in the
  benches (`yuan_ref` in `tests/f182_cascade.rs`, `yuan_ref183` in `tests/f183_cascade.rs`) with that note.
- **Committed as is**:
  - the three flags off by default;
  - the controls;
  - the report with the stop.
- Production is unchanged: the C1 hash, the guard 6/6 Match, 667 lib and 33 viz tests.
