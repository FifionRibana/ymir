# F168 — declared before any measurement (2026-10-10)

**C1 (b), part 1: a continent that fills the map and can be sailed around.**

**Production is unchanged**: C1's production default and the guard 6 / 6 are not touched. The cascade uses a **C1
parameter profile** (« continent profile ») behind a flag.

**Frozen**:
- T = 5.62 Myr;
- the Airy mapping (the oceanic plates at S̃ = 0.2);
- the sources' drive;
- F164's cascade;
- the viz's offset rule;
- **Davis-Suppe's parameters, not touched**.

**The one constant set on our worlds is a design target, not a physical calibration**: the initial continental
fraction, chosen once on the témoin for 55–60 % land at the end, then frozen for the four seeds.

**Not blind**:
- F167's readings (the témoin's land 31 091 km² = 19.4 % of the 400 × 400 km map);
- the code read below.

No profile other than production has been run.

## A1 — how C1 makes its land

- **`ContinentalClusterParams`** (`tectonics_c1/init_r7/clustering.rs:58–78, 150–200`):
  - **`continental_fraction`** (default 0.29) is **a fraction of the PLATES**, not of the area. The target is
    round(`num_plates` · f) continental plates, clamped to [1, `num_plates`];
  - `seed_cluster_count` (default 1) BFS seeds pick them across the plate adjacency, from `seed` (ChaCha8);
  - **so the land is quantised in whole plates**: with `num_plates` = 8 (default, `init_r7/mod.rs:130–135`),
    f ∈ [0.1875, 0.3125) gives 2 plates, [0.3125, 0.4375) 3, [0.4375, 0.5625) 4, [0.5625, 0.6875) 5, [0.6875, 0.8125) 6;
  - the area follows each seed's Voronoi cells.
- **`craton_shield_fraction`** (Some(0.15), `init_r7/mod.rs:136–150`): the share of the cratonic area kept as thick
  high shield (S̃ × 1.25); the rest is normal continental crust.
- **`num_plates`** (8): the plate count, a fragmentation lever (and the quantum above).
- **From the initial to the final land under Airy** (land = S̃_src ≥ 0.915, S̃_src = s − C1's erosion):
  - the continental plates start at S̃ = 1.0 (land) and stay in place;
  - Davis-Suppe thickens upper-plate wedges;
  - the subduction consumes oceanic crust and adds arc mass;
  - the rift thinning can sink thinned continental crust below 0.915;
  - the type flips add land at trenches.
  - **Measured per seed and per size**: the land at t₀ and at the end, and the cells lost (rift, others) and gained.
- **The continental cells' velocity is zero** (`apply_continental_rigidity`, F165-A1): **the continents never move**.
  Their position and outline are fixed at init by the Voronoi plates; only the boundaries' classification and the
  sources depend on the plate velocities.

## A2 — Davis-Suppe's wedge (measured only, nothing changed)

- **The parameters** (`davis_suppe/source_term.rs:53–84, 133–141`): coupling 2.0, h_max 2.5, **l_taper 4, l_decay 6,
  max_distance 30 cells** (O-C rate ×6 in `C1Closures::default`, #155).
  - **Their origin**: « Defaults selected for Phase 1.2 visual demonstration at 64² » (`source_term.rs:53`), a visual
    calibration (doc §11.1), not a literature value.
  - **Their metric**: cells of the reference 64² grid (`DS_REFERENCE_GRID`, rescaled ∝ nx, #147), i.e. a fraction of
    the domain. On our 400 km map that is 25 / 37.5 / 187.5 km; ×7.5 « signified », 187 / 281 / 1 406 km.
- **The profile**: per seed and per size, the mean Σ Δs_DS of the final land cells (64²) by Chebyshev distance
  (periodic) to the nearest current convergent cell or suture: bins 0, 1, 2 … 10 cells and > 10 cells.
- **The share of the land with Σ Δs_DS > 0.05.**

## A3 — the offset and the border

- **The viz's offset** (`ymir-viz/src/bridge/c1/hd.rs:451–455`): `land_topology(coarse, sea)` (`tectonics_c1/land_topology.rs`,
  periodic 4-connected components, the largest mass's circular extent per axis) gives its torus centre `center_cell`;
  the offset is centre − 32 per axis, so the largest mass is centred.
  - It is applied to the HD run, the export window (`sample_origin`) and the Cascade window (`workspace.rs:936`).
  - **It guarantees nothing by itself**: `land_topology` flags a « band » (`wraps_x/y`) when the largest mass occupies
    every column or row. The offset only centres.
  - **The benches until F167 used the témoin's canonical roll (6, 37) for every seed**, not the viz's offset. That is
    why seed 9's land touched the map border (the F166 rim).
- **The cascade works on the periodic field**:
  - `upsample2`, `restrict` and the warp sample periodically;
  - `compute_flow` routes D8 across the wrap (`terrain/flow.rs:256`);
  - but the solver makes a cell whose D8 receiver leaves the map its own fixed base (`erosion/stream_power.rs:651`).
  - So land on the map border after the roll is never incised and rises under U: **the rim of F166**.
  - It disappears only where the border is sea; that is A3's link between the offset and the rim.
- **The circumnavigation test** (declared):
  1. on the torus, the largest land mass (periodic 8-connected) **does not wrap**: a periodic BFS tracking the
     unwrapped coordinates never reaches a cell at two different unwrapped offsets (both winding numbers zero);
  2. **the sea width along the best loop**: r* = the largest r such that the component of the r-dilated land (all the
     land, Chebyshev, periodic) that contains the main mass still does not wrap.
     - A loop of sea at least 2r* + 1 cells wide then surrounds it. Reported in cells and km.
  - It is read at the end of C1 (64², h_src > 0 under Airy) and at 2 048² (the final field > 0).
  - Also reported: `land_topology`'s band flags, and **the land cells touching the window border after the viz's
    offset** (islets included), without a rule.
- **The offset in this bench**: the viz's rule applied to the cascade's own final 64² field (the Airy h_src, the last
  snapshot), per seed and per size.

## B — the variants

- **The trials (the témoin only)**: f = 0.29 (L0, 2 plates), 0.375 (3), 0.5 (4), 0.625 (5), 0.75 (6).
  - **The land share** = the cells of the history's 256² physics level above 0 m, over the whole torus
    (offset-independent).
- **L1** = the trial nearest 40 %; **L2** = the trial inside 55–60 %, else the nearest.
  - **If no trial lies in 50–65 %, or the response is not monotone in f, the « unreachable » stop fires.**
- **Frozen for the four seeds** (the témoin, 42, 1, 9).
- **The negative control** (permanent): the continent profile at L0 (f = 0.29) is `Phase2InitParams::default()` and
  gives C1 bit for bit. In the bench, L0 with the canonical roll reproduces F167's témoin peak (2 670 m).
- **The coast reference** (« degraded against F167 »): L0 under the viz's offset (F167's configuration, the viz's
  framing).
- **The viz**: « Profil C1 : production / continent L1 / continent L2 ».

## C — the measures (4 seeds × L0, L1, L2)

- **The land**:
  - its share of the torus; the masses > 1 000 km² (periodic 8-connected, at 256²);
  - the circumnavigation test (64² end of C1, 2 048²);
  - the border cells.
- **The interior**: the share of the land farther than 20 km and 40 km from the coast (Euclidean distance transform
  on the 256² land, × 1.5625 km).
- **Macro**: the classes against Europe.
  - **The geometric hypothesis**: the mountain share against the land area, L0 → L2 (the slope per 10 000 km²).
- **C1's tectonics**:
  - the plates at the end;
  - the continent–continent collision cells (Convergent with a continental cell on both sides, F165's rule);
  - the active margins (upper-plate cells);
  - A2's share;
  - the cratons: cells, altitude p50 at 256², and their blockiness (F165's block index of the craton mask's boundary,
    lattice 4 at 256²).
- **The peaks and the slope guard.**
- **The lakes.**
- **The border rim**: the border land cells after the offset and their mean altitude against the land's.
- **The texture** (whole and mountain class), F per class, **the coast** (F165's B4), **the regularity** (F167's
  baseline) for L2.
- **The images** (north up, the viz's offset):
  - L0 | L1 | L2 per seed, shaded and classes, at 256² and 2 048²;
  - the sea loop drawn on the 256² views (the boundary of the r*-dilated main blob);
  - the Σ Δs_DS profile;
  - the 25 km coast crops against Corsica.
- **The cost**: CPU seconds per world, C1 included.

## R — the stop rules (the brief's)

- **Macro success in L2**, on ≥ 3 of 4 seeds:
  - the mountains in 16.5–37.2 %;
  - plain + plateau in 37.5–84.5 %;
  - circumnavigable;
  - no coast degraded by more than 10 % against F167.
  - Then F169 (the sea-level rise) on L2.
- **The mountains still above 37.2 % in L2**: decompose them (Davis-Suppe by distance to the margin, the cratons, the
  arcs), and stop. F169 is then the wedge's reach.
- **The circumnavigation lost on ≥ 1 seed in L2**: report which size loses it, and propose the most neutral guarantee.
  Implement nothing.
- **The target unreachable** (no trial in 50–65 %, or not monotone): stop and report.
- **Forbidden**: as F165–F167. No constant set on our worlds but the land target.
