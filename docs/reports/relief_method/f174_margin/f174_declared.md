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
