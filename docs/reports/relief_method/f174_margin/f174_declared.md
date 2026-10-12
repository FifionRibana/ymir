# Finding 174 — declared before any measurement

Written before the bench `f174_cascade.rs` was run. The windows were built and looked at (`prep_coasts_output.txt`, a
montage checked for orientation) before this file was written, but **none of the four coast measures was computed on them**.
The known facts at writing time are listed where they matter (the Brittany window's maximum is 297 m).

## The configuration (all of Part A)

- **The panel**: the seeds 8, 18, 19, 23, 29 (F173). Each seed's continent profile is `select_continent_default` (16 plates,
  then 20, then 24), exactly as F173.
- **The worlds**: the C1 64² history (Airy sources, T_F166), the physics level at 256² (300 steps, rebound, the torus), the
  rock K2 (F171: K1's 1.5 km cover on a basement 10× harder, plus the heterogeneous K, σ_lnK 0.5, L 25 km), the frozen F164
  chain to 2 048² in torus mode, **no sea** (F173's « none » branch, bit for bit the same calls).
- **The ages**: B0 (linear time) and old (A = 300 Myr), as the brief asks. A1, A2 (our worlds) and A3 use the same 10 worlds.
- The old témoin (10481999410520546993) is not measured; it serves only the bit-for-bit controls.

## Part 0 — the Spratt re-zero

- `SeaCurve::from_noaa_spratt2016` subtracts the series' value at 0 ka from every level, so the curve reads 0 m today
  (Spratt & Lisiecki's stack reads +8.96 m at 0 ka). An alignment of reference, not a calibration: the shape and the
  amplitudes are untouched.
- The value before the series (`before_m`, 0 m) is unchanged.
- The test: the parsed curve reads exactly 0 at 0 ka, and every other level is shifted by the same constant.
- F173's measures used the un-zeroed series; they are not re-run (the sea is out of the profile, decision (a)).

## A1 — is the coastline inherited from C1?

- **The C1 contour**: s = 0.915 is where Airy's altitude is 0 (b·35·s = b·35 − 0.45 km, b = 500 / 3 300). The contour is taken
  on the field the cascade receives: C1's last Airy snapshot (C1's σ = 0.5 blur, oceanic cells read at s = 0.2), normalised
  (`C1History::normalized`), upsampled by `upsample_to` to 256² — the physics level's drive, exactly — and on to 2 048² by
  `upsample_to` again. Land = the value above the cascade's sea level (0.5 normalised).
  - The chain's warp (W-tous) is **not** applied to the C1 contour: the warp is part of what may move the coast, and is
    therefore part of what is measured.
- **The final coastline**: at 2 048², the land on the torus (`land_t`: the ocean is the largest periodic sea component); a
  coastline cell is a land cell with a 4-neighbour ocean. The C1 coastline cell is defined the same way on the C1 field.
- **The distance**: for each final coastline cell, the Euclidean distance (periodic EDT) to the nearest C1 coastline cell, in
  cells of 256² (2 048² cells ÷ 8). The share of the coast's length is the share of the coastline cells.
- **Reported**: the shares within 1, 2 and 5 cells of 256², and the median and p90 distance. Also, as a split of the cause,
  the same distances at 256² between the physics level's final coast and the C1 contour (no chain, no warp involved).
- **P1 is judged** on the 2 048² share within 1 cell: P1 holds at an age if that share is ≥ 90 % on ≥ 4 of the 5 seeds.
  **P1 is refuted** (the stop rule « la côte n'est pas héritée de C1 ») if it fails at both ages. If it holds at one age only,
  that is reported as mixed.

## A2 — the reference coasts by type

### The rule (computable on the DEM windows and on our worlds)

- A coastline cell is a land cell with a 4-neighbour sea cell, where the sea is the **largest 4-connected sea component**
  (of the window, or of the torus on our worlds). Water not connected to it (a closed lagoon, a lake) counts as land, as on our
  worlds.
- For each coastline cell, the p90 of the altitude of the land cells within 5 km (a disk of radius 25.6 cells at 195.3 m).
  **The segment is rocky if that p90 > 200 m, low otherwise** (the brief's example, taken as is).
- A window's type is the type of ≥ 80 % of its coastline cells; otherwise it is mixed.
- **Validation on the real windows**: the rule's type is compared with the declared (geographic) type of the 8 windows. The
  bands below are built from the **declared** type. If the rule disagrees on a window, that is reported, the rule is flagged
  for the author, and it is not tuned here. Known before writing: the Brittany window's maximum is 297 m, so its p90 near the
  coast may fall under 200 m and the rule may call it low.

### The windows (50 km, 256 × 256 cells of 195.3125 m, the same DEM and pipeline as F161's Corsica)

Copernicus DEM GLO-30 (1″), projected to a local metric frame and averaged into 48.83 m cells (the sea at −50 m, land the mean
of its land pixels), then block-averaged ×4 — `prep_corse.py`'s method exactly (`prep_coasts.py`). A window is a mosaic of
the tiles its box touches; a tile Copernicus does not publish (HTTP 404) has no land and is sea (Jutland's E007 tiles).

| window | centre | declared type | why |
|---|---|---|---|
| Galicia, Rías Baixas | 42.40 N, 8.70 W | rocky | drowned river valleys (rías) in granite, a classic rocky coast |
| Dalmatia, Split–Trogir | 43.55 N, 16.35 E | rocky | karst coast with the islands Brač, Šolta, Čiovo; a drowned folded margin |
| North Brittany (Trégor, Morlaix) | 48.72 N, 3.65 W | rocky | granite and gneiss coast with skerries; low relief (max 297 m) |
| Corsica west (Sagone, Ajaccio) | 42.20 N, 8.65 E | rocky | F161's reference; a crop of `corse_1024` (other islands removed by F161) |
| Landes (Mimizan) | 44.20 N, 1.35 W | low | a straight sandy coast with dunes, the textbook smooth coast |
| Vendée (Les Sables-d'Olonne) | 46.55 N, 1.65 W | low | sandy beaches, marshes and small rocky points (the brief's Vendée-Charente) |
| West Jutland (Ringkøbing) | 56.05 N, 8.30 E | low | a North Sea barrier coast with a lagoon (the brief's Frisia / North Sea) |
| Golfe du Lion (Agde, Thau) | 43.25 N, 3.45 E | low | lagoons behind a sandy lido, a deltaic plain |

- Frisia / the Netherlands was not taken: much of its land is below 0 m, and the « land = DEM > 0 » rule would drown the
  polders. West Jutland stands for the North Sea barrier coast instead.
- The windows were chosen by their geology and by eye from maps, not by any coast measure.

### The measures (physio::coast's four, through a window variant)

- The box-counting dimension over boxes of 2–64 cells (0.39–12.5 km, the closest to the brief's 200 m–10 km on the same
  lattice as Corsica's F161 reference), the length ratio L(2)/L(32), the aligned share of the 8-step chords, and the block
  index (lattice 32), exactly `coast_t`'s code.
- **The window variant**: the window's border is never a coast. Only land/sea edges between two cells inside the window
  count; a coastline cell needs an in-window sea neighbour; the crack contours stop at the border. On our worlds the same
  variant is applied to 256² windows cut from the 2 048² field, so real and generated coasts are measured identically.

### The bands and our worlds

- **A band per type** = [min, max] of the type's 4 real windows, per measure, not widened.
- **Our worlds**: the 2 048² field cut into 64 non-overlapping windows of 256² (50 km). A window is eligible if it has
  ≥ 128 coastline cells (≈ 25 km) and a type by the rule (≥ 80 %). Per world and type, the median of each measure over its
  eligible windows, counted against that type's band (0–4). Also reported: the worlds' coast with `coast_torus` (F173's
  measure) for continuity.

## A3 — convergence in the time step

- The physics level at N_phys = 300 (the normal step) against N_phys = 600 (the step halved), the same T, the same age, the
  same rock, everything else identical. Under the age mapping, step_durations solves r for 600 steps, so each step is
  ≈ halved (exactly halved at A = T).
- **Reported**, at 256² on the physics level's output (m): |Δz| p50, p99 and max over all cells and over land; the classes'
  largest move (points) and the share of cells whose class changed; the network: the cells of Strahler order ≥ 3
  (`hydro::hydrology`, the declared threshold) in one run with no order-≥ 3 cell within 1 cell (3 × 3) in the other, as a
  share of the order-≥ 3 cells, the mean of both directions.
- **P3 is judged** on max |Δz| > 200 m on ≥ 3 seeds (at either age) and the classes' largest move < 1 point everywhere.

## A4 — why every coast was active

- **The per-step shares**: along the continent's perimeter at 64² (continental cells with a 4-neighbour oceanic cell, by
  `plate_type` at that step), the share whose 3 × 3 neighbourhood holds a convergent, else divergent, else transform boundary,
  else no boundary (F173's priority), after `retarget_upper_plate_continental`; for every C1 step. Reported: the mean over the
  first 50 steps, over all steps, and the series every 25 steps.
- **The same-plate share**: the share of perimeter cells whose oceanic neighbour belongs to the **same** plate. A passive
  margin (continent and ocean on one plate) can exist only there.
- **P4 is judged** on the mean convergent share over the first 50 steps > 60 % (on the panel's mean and per seed).
- **The code reading** (not a measure): where the continental velocity is zeroed, where plate_id and the boundary types come
  from, and what follows. **No correction**; proposals in order of preference.

## Part B — status at declaration time

- **B2's sources could not be read from here**: Péron-Pinvidic et al. 2013 (Marine and Petroleum Geology 43, 21–47) is behind
  HAL's bot protection; Harris et al. 2014 (Marine Geology 352, 4–24) is an Elsevier page that does not render. Neither
  refusal was circumvented. Christensen & Mooney 1995 (JGR 100, 9761–9788) was not fetched. None is in `docs/refs`.
- So the margin's profile, w_n and B-hétéro's σ cannot be declared from a read source. **Part B stops before it starts**,
  whatever P1 says, and the author is asked to add the three articles to `docs/refs`. Part C is not run.
- What would be declared once the sources are read: the profile shape (proximal / necking / distal domains), w_n and the shelf
  width expected on the map (the layout is compressed ×7.5, the slopes are real), the ×2 and ÷2 variants, B-hétéro's σ and
  correlation length.

## The checks

- `checks154.sh` before (`checks_before.txt`) and after (`checks_after.txt`): cargo check, the lib tests, the viz tests, the
  guard 6/6. Production unchanged; no cache key touched (no ALGO_* bump expected).
- Cost reported in CPU seconds per world.

# Parts B and C — declared before any measurement (2026-10-12, after Part A's commit `31384f5`)

Written after reading the three sources and C1's code, before the bench `f174bc` was written or run. Journal page numbers
are the articles' own.

## The constraint (the author's)

« attention à ne pas retomber dans le biais d'intégrer un solveur dans les mécanismes sauf si réellement nécessaire. Mais c'est
justement ce qu'on a essayé de fuir avec C1! » So: **initial conditions only** (continents inside their plates, the rim, a
prescribed margin profile, the B-hétéro field), then C1 runs **unchanged**. No advection of `plate_id`, no simulated rifting
or stretching, no solver. C1's code is not touched: the margin is applied to the initial state the continent profile builds
(`cascade::margin`), so with the margin off C1 is bit for bit by construction.

## A finding of the code reading, declared before measuring (it shapes B1's verification)

- **The subduction closure cannot consume a rim of the continent's own plate.** `apply_subduction_step`
  (`closures/subduction/source_term.rs`, l. 169–228) acts on an oceanic convergent cell only if a 4-neighbour of a **different
  plate is continental**; otherwise the boundary is « oceanic-oceanic convergent, not handled by this closure ». With an oceanic
  rim of the same plate between the coast and the plate boundary, the boundary is ocean against ocean: nothing is consumed and
  nothing is promoted. The arc mass is placed on continental cells within `arc_distance` of the consumed cell, so it cannot
  reach the continent either.
- `retarget_upper_plate_continental` and `oc_override_seed_mask` (Davis-Suppe's O-C seeds) likewise need a continental cell on
  the boundary; on an O-O boundary the velocity-based upper plate stays and Davis-Suppe takes its C-C / fallback « dome »
  geometry from those seeds.
- So the expectation is that the rim is **not** consumed and the coasts stay passive by construction. This is measured (B1's
  verification), not assumed. If it holds, the stop rule « no rim consumption » reports it; changing the closure's partner rule
  would be a change to C1's closures, not an initial condition, and is left to the author.

## B1 — continents inside their plates

- The continent profile builds C1's initial state as before (`C1Profile::init_state`). Then, for every cell of a continental
  plate, d = the Euclidean distance (cell centres, periodic, in 64² cells) to the nearest cell of an **oceanic plate**
  (boundaries between two continental plates are interior to the continent and get no rim).
- **The rim**: the cells with d ≤ w_rim become `Oceanic` **on the same plate** (`plate_id` unchanged), with C1's oceanic init
  values: s = 0.2 (`init_s_field`'s oceanic reference, 7 km), the age by C1's own rule (`init_age_field_ridge_aligned` re-run on
  the new types: ridge value on divergent cells, the oceanic baseline elsewhere), and out of the cratonic mask.
- **w_rim = 2.5 cells of 64² (15.6 km on the map), a design value.** Why: the rim must hold more than C1's 1-cell boundary
  smoothstep (`InitMode::Uniform`, width 1) and the Airy record's σ = 0.5 cell blur, so the coast never sits on the boundary
  cell; 2.5 is the middle of the brief's 2–3 cells. At the ×7.5 layout compression it stands for ~117 km of real outer domain
  and ocean (Péron-Pinvidic's distal + outer domains are ~100–150 km wide, p. 34–35). Variants ×0.5 (1.25) and ×2 (5),
  measured only.
- **The selection is redone with this initial condition** (`select_continent_default_margin`): the same rules (16 plates, then
  20, then 24; k ∈ [⌈0.35 N⌉, ⌊0.75 N⌋]; land 50–62 %; circumnavigable with r* ≥ 1), each candidate's C1 run starting from the
  margined state. Run per variant.
- **The verification** (every C1 step): the continental perimeter (continental cells with an oceanic 4-neighbour, any plate);
  **active** = a convergent boundary cell in its 3 × 3 (A4's rule); **rim consumed** = a 4-neighbour of another plate (no rim
  left between the coast and the boundary). Reported: both shares at the end, the series every 25 steps, the first step a
  share exceeds 5 %, and the corresponding age of the history (the age mapping of each profile).

## B2 — the thinned margin (on the continental cells, inside the rim)

**The source**: Péron-Pinvidic, Manatschal & Osmundsen 2013, Marine and Petroleum Geology 43, 21–47.
- The domains, seawards: proximal, necking, distal, outer, oceanic (p. 24, Fig. 3).
- The basement thins oceanwards from > 30 km in proximal settings to ~0–5 km at the most distal (p. 22, Fig. 2).
- The necking domain: drastic thinning from ~30 km to < 10 km (p. 24); its outer end, the taper break, is the nearest point to
  the coast where the crust is ≤ 10 km (p. 25). The distal domain: basement < 10 km (p. 25); 3–7 km under the Angola sag
  basin (p. 34).
- Widths: the proximal domain is a > 200 km wide corridor at Iberia (p. 27; > 500 km at the Grand Banks, p. 27); the necking
  domain ~50 km (Møre) to > 150 km (Vøring), ~60 (Deep Galicia) to ~150 km (Southern Iberia), ~100 km (southern
  Angola–Gabon) (p. 35); the distal sag basin ~150 km wide (Angola, p. 34).

**The profile** (C1's convention s = 1 ↔ 35 km), piecewise linear in x = the distance inward from the continent's edge (the
inner side of the rim, x = d − w_rim):

| domain | thickness (km) | s | real width (km) | map width at ×7.5 (km) | 64² cells |
|---|---|---|---|---|---|
| distal (from the edge) | 7 → 10 | 0.200 → 0.286 | 150 (p. 34) | 20.0 | 3.2 |
| necking | 10 → 30 | 0.286 → 0.857 | 100 (p. 35) | 13.3 | 2.1 |
| proximal | 30 → 35 (or s_craton) | 0.857 → s_int | 200 (p. 27) | 26.7 | 4.3 |
| interior | unchanged | unchanged | | | |

- The outer and oceanic domains are the rim (oceanic cells, s = 0.2).
- s_int = 1, or C1's craton ratio on the cratonic shield; the profile is s(x) = 0.2 + (s_int − 0.2)·g(x), g through 0, 0.107,
  0.821, 1 at the knots. Cells beyond the proximal domain keep their initial s.
- **The scale**: the layout is compressed ×7.5 (the map is 400 km for a continent of Earth size), so the main widths are the
  real ones ÷ 7.5. The altitudes are Airy's, so the margin's slopes on the map are ×7.5 the real ones. Variants: all three
  widths ×2 and ÷2, measured only.
- **The expected shelf on the map** (pure Airy on the initial profile, before C1, erosion and the chain): the coast
  (s = 0.915) lies in the proximal domain, 41 % of its width from the necking; Airy puts the proximal / necking knot at
  −448 m, the first marked slope increase. So the shelf from the coast to that knot is **10.9 km on the map** (82 km real),
  its break at **~450 m**; ×2 → 21.8 km, ÷2 → 5.5 km.

## B3 — the shelf reference

**The source**: Harris, Macmillan-Lawler, Rupp & Baker 2014, Marine Geology 352, 4–24.
- The shelf width (the shelf break to the nearest land, every 0.1°): passive margins **88.2 ± 0.7 km**, active 31 ± 0.4 km, all
  57 ± 0.41 km (Table 2, p. 11; the abstract rounds to 88 and 31, the caption gives 84 ± 0.66 for passive — recorded as an
  inconsistency of the article). Per ocean, the passive means span **34.9 km (North Pacific) to 123.0 km (South Atlantic)**
  (Table 2).
- The shelf break is drawn where the slope increases markedly (p. 6); **the article gives no statistic of its depth.** The
  ~140 m value (Pinet 2003, Gross 1972, via Wikipedia) is therefore not verified by a read source; the break depth is
  reported, without a band.
- **The band (declared)**: the per-ocean passive means ÷ 7.5 = **4.65–16.4 km on the map**, centre 88.2 / 7.5 = 11.8 km. The
  margin's success reads the median shelf width against it.

**The shelf measure** (2 048², the final field, sea at 0 m): from each coastline cell (land with an ocean 4-neighbour), the
steepest-descent path through the ocean (D8, strictly descending), up to 60 km or −3 000 m. The break is the first point,
at least 1 km out, where the mean slope over the next 1 km is ≥ 3 × the mean slope from the coast to that point (Harris's
shelves are ~10 × gentler than the slopes; 3 is a conservative detection). No such point: no shelf (width 0). Reported: the
median and p90 of the width (km), the median depth of the break, and the share of coastline cells with a shelf.

## B4 — B-hétéro (measured only)

**The source**: Christensen & Mooney 1995, JGR 100(B7), 9761–9788.
- Continental crust: 39.2 km ± 8.5 km over all observations (p. 9764, Fig. 2); weighted 41.0 ± 6.2 km; **shields and platforms
  41.5 ± 5.8 km** (Table 3).
- **σ**: the initial continent is a platform before C1 builds its orogens, so the shields-and-platforms spread: σ_rel =
  5.8 / 41.5 = **0.14**.
- **The field**: n = a periodic fBm (`noise01_t`, 1 octave, wavelength **40 km on the map** ≈ 300 km real, a design value of
  the margin's own scale), standardised to zero mean and unit variance over the continental cells; s_int ← s_int · (1 + 0.14 n),
  under the main margin (B-hétéro = the main margin + the field). The ramp then reads the perturbed s_int, so the coast's level
  curve moves with it.
- Its effect on the coast and on the macro (classes, median, land) is reported.

## B5 — the controls (permanent tests)

- The margin off (`Margin::off()`, w_rim = 0, no ramp, no field): the initial state is the profile's bit for bit, and the
  default selection is `select_continent_default` bit for bit.
- C1 untouched (no C1 source changed).
- w_rim = 0 and a ramp of zero width: the initial state is unchanged except the ramp's knots collapse — tested as identical to
  the off state.

## Part C — the measures

**The worlds**: the panel 8, 18, 19, 23, 29; K2; the torus; no sea; the four ages (B0, young 20 Myr, intermediate 100 Myr,
old 300 Myr); the physics level at 256², the chain to 2 048². Per seed and age: **Part A's world** (no margin, F173's « none »
exactly) and **the main margin** (its own selection). **The variants** (w_n ×2, ÷2, w_rim ×0.5, ×2, B-hétéro), each with its own
selection, at **young and old** only (cost).

**Reported**:
- The selection: N, k, land, r*, found on 5 (per variant).
- The margins: the active and rim-consumed shares at the end; the series; the first step past 5 %.
- The rampart: W₁₀ at young and rampart_W (F173's `wall_index`).
- The shelf (B3's measure): median, p90, break depth, share with a shelf, against 4.65–16.4 km.
- The coast:
  - A1 redone (the share within 1 cell of 256² of the C1 contour of that world);
  - **by type, by the rock at the coastline** (the decision at Part A's commit): the 256² rock at the surface (F171's
    cover / basement / shield, nearest-upsampled to 2 048²); a 256² window of the 2 048² field is rocky if ≥ 80 % of its
    coastline cells are basement or shield, low if ≥ 80 % are cover, otherwise mixed; eligible with ≥ 128 coastline cells.
    The rocky medians against the rocky reference (dimension 1.15–1.22, length ratio 1.43–1.56), the low medians against the
    sandy one (dimension 0.95–1.00) and the lagoon one (1.25–1.35);
  - the seaside slope: the median slope (°) of the 2 048² land coastline cells' 3 × 3 (Horn);
  - the dendritic sea arms: the 8-connected components of the ocean at 2 048² lying inside C1's land by ≥ 2 cells of 256²,
    with ≥ 1 km²; their number and their penetration (the largest distance from C1's coast reached, km).
- The macro: the classes, the median, the land, against Part A's world (no rule; the rampart was a bias).
- **Statistical convergence** (seeds 8 and 23, the main margin, B0 and old; the physics level at 300 against 600 steps): the
  drainage density and λ (`read_cell`), the sinuosity S₅ (`planform`), the classes. Q8 is judged on < 5 % relative.
- Lakes (F173's count and area); texture (relief p50, slope p50 / p90 at 2 048²).
- Images (north up): the whole maps (Part A | margin | B-hétéro, at young); the active / passive perimeter at the end of C1
  (64²); 25 km coast crops (Corsica west | the Landes | Part A | margin); a transect across a passive and, if any, an active
  margin (altitude and s).
- Cost: CPU seconds per world, the selection included.

**The success** (the brief's, on ≥ 4 of 5 seeds): W₁₀ at young down ≥ 40 % against Part A; the active share at the end in
10–60 %; the median shelf width in 4.65–16.4 km. The coast's dimension is reported, not part of the success.
**The stops**: the selection fails on ≥ 2 seeds → stop and report; < 5 % active on every seed → reported (« the subduction
closure does not do what is expected of it »).
