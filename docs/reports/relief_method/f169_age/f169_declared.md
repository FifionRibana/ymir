# F169 — declared before any measurement (2026-10-10)

**A world that does not depend on its framing, a continent that can be sailed around, and the continent's age.**

**Production is unchanged**: everything goes through cascade-only flags (`None` = F168 bit for bit) and the C1 « continent »
profile. The guard 6 / 6 is not touched.

**Frozen**: T = 5.62 Myr, K, the Airy mapping, the sources' drive, F164's cascade, Davis-Suppe. The ages A are the author's
and the reviewer's **design values**, not calibrations.

**Not blind**: F168's readings and the code read below. No F169 variant has been run.

**Part 0**: F168 committed (`3271a44`). Its ADR carries the five requested points; the « At commit » block makes the invariance
defect and the reviewer's error explicit, with the author's three answers verbatim.

## A1 — the fully periodic cascade (`torus`, a cascade flag)

### The inventory of what depends on the window (read in the code)

| where | what | today | under `torus` |
|---|---|---|---|
| `stream_power.rs:650` | a D8 receiver leaving the map | the cell becomes a fixed base (the rim) | wraps; the base level is the sea cells only |
| `stream_power.rs:1126` `linear_diffusion` | the border rows / columns | never diffused | periodic 5-point Laplacian on every cell |
| `stream_power.rs:1046` `talus_sweep` | off-map neighbours skipped; the sweep order ties broken by index | | periodic neighbours; ties by the cell's **terrain index** |
| `stream_power.rs:941` `receiver_hash` (random receiver, off in the frozen cascade) | hashed on the window index | | hashed on the terrain index; periodic |
| `stream_power.rs:850` lateral erosion | stops at the window edge | | periodic |
| `flow.rs:878` `mfd_accumulation` | float sums, order ties by index | | ties by the terrain index |
| `amplify.rs:346` `ocean_mask` | the ocean = sea **4-connected to the window border** (decides the D5 basins) | | the ocean = sea periodically 4-connected to the **largest** sea component (ties: size, then the smallest terrain index) |
| `amplify.rs` `norms`, `gather`, `steepest`, `talus`, `deposit` | off-map neighbours skipped | | periodic |
| `amplify.rs:575` `exact_area` | sort ties by index | | ties by the terrain index |
| `amplify.rs:403` `noise01` (π, talus ν, hardness, the warp's d, the physics roughness) | OpenSimplex FBM at the **window's** km coordinates, not periodic | | **periodic on the torus** (a 4D torus embedding of the same OpenSimplex sources, radius per octave so the local wavelength is kept), sampled at the **terrain's** coordinates |
| `cascade.rs:688`, `history.rs:367` the roughness; `amplify.rs:749` π | the land mean / RMS summed in window order (f64) | | summed in terrain order |
| `cascade.rs:559` the rebound (`filter_isotropic`, FFT) | FFT of the window-ordered field | | FFT in the terrain frame (rolled to it, filtered, rolled back) |
| `compute_flow` (pit fill, Garbrecht-Martz flats, D8, D8 accumulation) | already periodic; its outputs do not depend on the scan order (the fill value is unique, the BFS distances are unique, the D8 accumulation adds integers) | | unchanged |
| `upsample2`, `restrict`, `sample_cr`, `warp_grid`, `jacobian` | already periodic | | unchanged |

- **The terrain frame**: the input carries its origin, the shift (in 64² cells) between the given grid and C1's own frame. The
  terrain index of a cell is its index in C1's frame at the level's resolution.
- So every tie-break, hash, noise, ordered sum and FFT is anchored on the terrain. Every neighbour operation is periodic.

### The invariance test (permanent; criteria declared here)

- Roll the 64² input by (dx, dy) and give the cascade that origin. The 2 048² output must equal the reference output rolled by
  (32·dx, 32·dy).
- **Tolerance: bit for bit.** Every source of a difference listed above is anchored, so none is expected: no FFT order, no
  tie-break, no sum order.
- **Shifts**: (0, 0) reference, (5, 0), (0, 23), (37, 50), the last crossing both edges.
- **Where**:
  - the bench runs the full chain (history physics 256², the frozen chain to 2 048²) on the 4 seeds;
  - the permanent unit test runs the same code at 64 → 256 (a physics history and one amplification level on a synthetic
    island).
- **If it is not bit for bit**: report the max |Δz| (m) and the share of differing D8 receivers. Then find and fix the source,
  or stop (the brief's rule): no physical measure on a non-invariant cascade.

## A2 — the offset at the end only

- The cascade computes in C1's frame (origin 0). The viz's offset rule (`hd.rs:451–455`, on the world's final 64² field)
  rolls the outputs for display and images only.
- **The rim test**: the mean altitude of the land touching the display window's border over the mean land altitude should be
  ≈ 1. The reading rule: 0.8–1.25 means no rim. Also reported on C1's frame border.

## A3 — the coast on the torus

`physio::coast_torus`: the same four measures with periodic land / sea edges. The box counts run on the torus (n divisible by
every box), and the crack contours close across the seam. The window border is never a coast. Corsica stays the reference
(non-periodic; its window border is sea).

## A4 — the continental plates' selection (an initial condition)

- **The candidates**, for f ∈ {0.5, 0.625}: the BFS cluster of round(8 f) plates grown from each of the 8 start plates (C1's
  BFS, unchanged), i.e. 16 candidates.
- **The criterion**: C1 run to the end under the profile, with the Airy record.
  - The end land (64², h_src > 0, whole torus) is within 50–62 %.
  - Circumnavigable (F168's test at 64²).
  - The largest r*.
- **Ties**: f = 0.625 first (the larger land, the author's direction), then the order in which C1's RNG would draw the start
  plates (its draw sequence; the first draw is today's start plate).
- **None qualifies**: reported per seed. The author then goes to more plates.
- **The negative control** (permanent): the selection off is F168 bit for bit (`init_c1_state_phase_2_r7` unchanged; the start
  plate goes through a separate entry point that, at today's RNG start plate, is the same state bit for bit).
- **Production untouched**: `ContinentalClusterParams` is not changed (it enters the tectonic cache key through `Debug`).

## A5 — B0

B0 = the selected profile + `torus` + the offset at the end, with linear time (A = T). It is measured with D's quantities and
compared with F168 L2. It is C's base.

## B — when C1 thickens (diagnostic)

- **B1**: per selected seed:
  - Σ Δs per term (Davis-Suppe, subduction, rift thinning, transport) per C1 step, summed over the end land;
  - per land cell, **the thickening age** = the first step at which its cumulative Σ Δs_DS reaches 50 % of its total (cells
    with Σ Δs_DS ≤ 0.05: no age), mapped;
  - the share of the land's Σ Δs_DS received in the last 10 %, 25 % and 50 % of the steps;
  - the C-C collision cells (F165's rule) at each snapshot.
- **B2 — « the rest »** (mountain class outside the cratons, the Davis-Suppe zones and the arcs, at 256², B0): altitude, R7,
  distance to the coast, its share of Σ rebound, and its altitude at t₀ (h_src snapshot 0).

## C — the continent's age A

### C1 — C1's steps in geological time

- N = 300 C1 steps; snapshots every 10 (30 intervals); one physics step per C1 step (N_phys = 300).
- d_N = T/N = 18 743 yr; d_k = d_N · r^(N−k), k = 1 … N; r ≥ 1 is solved (f64 bisection) so that Σ d_k = A.
- Physics step i ↔ C1 step i + 1, with dt_i = d_(i+1). On interval j, D_j = Σ d_k over its 10 steps, and the step's uplift
  du_i = Δh_j · d_(i+1) / D_j.
- **A = T gives r = 1 exactly** (special-cased, so no rounding enters): dt = T/N and du = Δh_j / 10, which is B0 bit for bit.
  A permanent test.

| A | r | d_1 (the oldest step) | the first half of the steps | the last 10 % of the steps |
|---|---|---|---|---|
| 5.62 Myr (control) | 1 | 18.7 kyr | 2.8 Myr | 0.56 Myr |
| 20 Myr | 1.007252 | 0.163 Myr | 14.9 Myr | 0.63 Myr |
| 100 Myr | 1.014687 | 1.47 Myr | 89.9 Myr | 0.70 Myr |
| 300 Myr | 1.019306 | 5.70 Myr | 283.9 Myr | 0.75 Myr |

### C2 — the time step

- **The linear diffusion** (explicit, 5-point): sub-steps s_i = max(4, ⌈w_i / 0.2⌉), w_i = D·dt_i / cell². The stability
  limit is 0.25; 0.2 is the margin. At A = T, w ≪ 0.8, so s = 4 as today (bit for bit).
- **The talus** is a repose closure (a threshold, no rate): once per step whatever dt, no sub-steps.
- **The rebound** is linear in the step's erosion: no limit.
- **The implicit stream power** is unconditionally stable. Reported: the largest dt (d_1) and the largest K·A^m·dt / dx over
  the land at that step (the « Courant » number, indicative).
- The cost is reported.

### C3 — the ages (design values)

The control A = T (B0); young 20 Myr; intermediate 100 Myr; **old 300 Myr (the proposed default)**.

### C4 — viz

In the Cascade window: a slider « Âge du continent (Ma) » with presets jeune / intermédiaire / vieux, and a view « Âge
d'épaississement ».

## D — the measures (4 seeds: témoin, 42, 1, 9; B0 and the 3 ages; physics 256², the frozen chain to 2 048²)

- **Infrastructure**: the invariance test, the rim, the circumnavigation and r*, the land %, the selected candidate.
- **Macro**:
  - the classes against Europe;
  - the land hypsometry's deciles against Europe's (40 / 81 / 117 / 153 / 195 / 262 / 383 / 580 / 889 m);
  - the peaks;
  - the mountain class decomposed by thickening age (the first, second and last third of the steps), the cratons, the arcs
    and the rest.
- **Plateaus and plains**: where (distance to the coast, to a convergent cell or suture), their altitude, their R7.
- **Lakes.**
- **Texture**: Corsica's battery in the mountain class and on the whole world; F by class; the slope guard.
- **The coast**: on the torus against Corsica. « Degraded » is measured against F167's coast, recomputed on the torus.
- **Regularity**: F167's baseline on the old profile.
- **Images** (north up, checked; the offset at the end):
  - the 4 seeds B0 | young | intermediate | old, shaded and classes, at 256² and 2 048²;
  - the thickening-age map;
  - the hypsometric curves against Europe;
  - a cross-section of an old belt and a recent belt (old profile).
- **Cost**: CPU seconds per world, per age.

## R — the stop rules (the brief's)

- **Infrastructure**:
  - the invariance test fails → stop, fix or report;
  - the selection fails on ≥ 2 seeds → stop and report (next: more plates).
- **Macro success (old profile)**, on ≥ 3 of 4 seeds:
  - mountains 16.5–37.2 %;
  - plain + plateau 37.5–84.5 %;
  - circumnavigable;
  - no coast degraded > 10 % against F167 (on the torus).
  - Then F170.
- **The old profile erases the high mountains** (peak < 1.5 km on ≥ 3): report B1's chronology. Next: the collisions' timing
  in C1. A is not tuned.
- **Mountains still > 37.2 % on the old profile**: decompose by thickening age and stop.
- **Forbidden**: as F165–F168. No constant set on our worlds.
