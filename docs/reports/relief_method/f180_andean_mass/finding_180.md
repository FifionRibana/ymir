# Finding 180 — the Andean seed behind the active coasts only, and Davis-Suppe's per-seed mass bound: the bound works (active coasts 12 → 4.9 km, collisions −7 %, mass control holds), but the Andean wedge still covers > 70 % of the coast on 11 seeds — the declared coverage stop fires after tier 1

**Status: built behind the continent profile's flag (two changes to Davis-Suppe's source chain: where the Andean seed sits, A1;
a per-seed mass bound, A2; no solver, no iteration), measured on C1; the coverage stop fired, so tiers 2–3 (256², 2 048²) were
not run; nothing adopted, nothing committed without the feu vert.** Production unchanged (`f175_c1_hash`).

- Declared before any measurement: `f180_declared.md` (A1's walk and seed, A2 with Dahlen p. 59, the sharing, s_in).
  Predictions: `f180_predictions.md`.
- The bench: `crates/ymir-core/tests/f180_cascade.rs`. The observer replays the loop's seeds and records the bound's statistics
  by family. Results in `results/`, assembled in `f180_bench_output.txt`.

## Part 0 — the defect measured first (F179 as committed, C1, 15 seeds)

- **The Andean coverage** (the coast cells seeded at some step, or within < 3 cells of an Andean-dominant Davis-Suppe cell):
  **68–100 %, median 94 %; above 70 % on 14 of 15 seeds.** That is 2.5–13 × the majority A1 share (median ratio 5.15).
- **R0 held on 15 of 15** (`images/coverage_F179_panel.png`).

## Part A — built

- **A1** (`andean_seed.rs`: `overriding_trench_partners`, `andean_seeds_active`, `andean_feeds_nearest`):
  - the coasts A1's walk calls active by the majority rule, the walk's trench and path;
  - the seed at d = max(4, d_front) on the path (the coast cell, or the path continued straight inland), if d ≤ 10;
  - the seed is skipped like a boundary seed;
  - behind `C1Ext::andean_active`.
- **A2** (`ds_mass_bound.rs`, `apply_ds_mass_bound`):
  - each Davis-Suppe cell attached to its nearest seed;
  - per seed and step, R_σ = Σ increments against S_σ = (v_rel · n̂) dt · 64 · s_in (Dahlen 1990, p. 59, eq. 1), shared S/k
    among the seeds a feeding cell feeds;
  - R_σ > S_σ → every increment × S_σ / R_σ (one pass);
  - every family; behind `C1Ext::ds_mass_bound`;
  - its record in `C1StepStats::ds_bound`.
- **A3** (permanent, `the_f180_controls`, all pass):
  - both flags off are F179 bit for bit;
  - an infinite supply is the unbounded run bit for bit;
  - **the mass control: after the bound, Σ increments / S_σ ≤ 1 + 10⁻⁹ for every seed at every step**;
  - production's hash unchanged.
- **A4**: the viz boxes « Graine andine derrière les côtes actives (C1) » and « Borne de masse de Davis-Suppe (C1) ».

## Tier 1 — C1, the 15 seeds of F176's selection

| | F176 | F179 | F179 + A1 | F179 + A2 | **F180** |
|---|---|---|---|---|---|
| Andean coverage, median (> 70 %) | 0 % | 94.4 % (14) | 91.8 % (14) | 83.8 % (14) | **75.9 % (11)** |
| coverage / majority A1, median (in 0.8–1.5) | — | 5.15 (0) | 5.12 (0) | 4.17 (0) | 3.67 (0) |
| A1 any / majority | 27.4 / 17.9 % | 24.9 / 18.3 % | 23.9 / 17.9 % | 25.5 / 19.0 % | 26.1 / 19.6 % |
| land, median (range) | 51.1 % (50.0–57.5) | 61.5 % (53.3–65.7) | 59.2 % (52.9–65.3) | 51.7 % (44.8–58.4) | 49.7 % (44.8–57.7) |
| circumnavigation lost | 0 | 1 | 1 | 0 | 1 (seed 23) |
| active 0–6 cells, median | 2.4 km | 12.0 km | 13.3 km | 5.0 km | **4.9 km** |
| passive 0–6 cells, median | 1.6 km | 4.9 km | 4.4 km | 1.6 km | 1.4 km |
| collision s_src / F176 (median) | 1 | 1.054 | 1.055 | 0.924 | **0.927** |
| bounded seed-steps: collision / O-C / Andean | — | — | — | 22 / 77 / 70 % | 23 / 86 / 63 % |

- **A2 (the bound) does what Dahlen's balance asks**:
  - it brings the active coasts' tectonic uplift from 12–13 km back to ~5 km, and the passive coasts to F176's level;
  - the land comes back to F176's (49.7 % median; the land added over F176 is 68–359 cells, 92–100 % in the band);
  - the collisions' s_src falls by 7 %;
  - the Andean seeds are bounded on 63–70 % of their steps.
- **A1 (the seed behind the active coasts only) barely moves the coverage**: F179 → F179 + A1 gives 94.4 → 91.8 %. With the
  bound added, it reaches 75.9 %.

## The stop, and its mechanism (measured, `f180_stop_mechanism`)

**The rule « an Andean coverage above 70 % of the coast on ≥ 3 seeds » fires on 11 seeds.** Tiers 2–3 were not run.

- **The covered passive coasts are close to a seed.** Their distance to the nearest cell ever seeded is a median of 1.4–5
  cells per seed (p90 6–13).
- **They receive a real Andean wedge.** The cumulative Andean Davis-Suppe in their 3-cell disc is 8–65 km of crust (the active
  coasts: 10–235 km). So these are not traces.
- **The cause.** The seed is now placed only behind majority-active coasts. But Davis-Suppe's wedge spreads from each seed
  isotropically within the plate, out to its reach of 30 cells (`wedge_distance_intra_plate_typed`: 8-connected, any
  direction). The wedge therefore runs along the coast into the neighbouring passive coasts. A real Andean wedge spreads inland
  across the trench's strike.
- **Where the coverage comes from**: the seed's position (A1) sets where the wedge starts, but not where it spreads. That is
  why A1 alone leaves the coverage at 92 %, and why the bound lowers it only through the smaller increments.

## Predictions

**The reviewer's**:

| | prediction | outcome |
|---|---|---|
| R0 | coverage ≥ 2 × majority on ≥ 10 / 15 under F179 | held (15) |
| R1 | coverage 0.8–1.5 × majority on ≥ 12 / 15 under F180 | refuted (0; median 3.67) |
| R2 | the bound cuts the collisions' s_src by < 10 % | held (−7 %) |
| R3 | active 0–6 uplift median 2–6 km | held (4.9 km) |
| R4 | the circumnavigation holds on 15 / 15 | refuted (lost on seed 23) |
| R5 | condition (v) | not measured (stop) |
| R6 | no rampart (W₁₀) | not measured (stop) |

**Mine**:
- held:
  - R0 (I predicted 3–5 ×; it came out 2.5–13 ×);
  - R1 refuted (I predicted 1.3–2.5 ×; it came out 2.4–11.6 ×, median 3.67);
  - R2 (2–10 %: −7 %);
  - R3 (3–7 km: 4.9 km);
  - « 0–1 seed loses the circumnavigation » (1);
  - the mass control.
- refuted: the land, 50–58 % predicted against 44.8–57.7 %.

## Open, for the author (nothing chosen)

- (a) **The wedge's direction**: Davis-Suppe's attachment is isotropic within the plate. An Andean wedge spreading across the
  trench's strike only (inland from its seed, not along the coast) would be a third change to this chain, for the author's
  next finding. The source would be needed (e.g. a wedge orthogonal to the trench, Dahlen's 2-D cross-section).
- (b) The coverage measure counts any Andean-dominant cell within 3 cells. The bound already caps the mass. Whether the rule
  should read the coverage by mass is the author's call; the rule as written fired.
- (c) Tiers 2–3 (~25 min) were not run and are resumable. The bound alone (« F179 + A2 ») keeps the circumnavigation on
  15 / 15.

## F180 — diagnostic (the author's complement; diagnostic only, no adoption)

The stop played and stays valid: **nothing is adopted in this finding.** The author's complement adds the measures below. The
budget, announced: ~30 min wall (C1 diagnostic ~3 min, 256² tier ~20 min, 2 048² tier ~5 min); actual: C1 diagnostic and 256²
tier 743 s, 2 048² tier 153 s.

**The reviewer's error (recorded)**:
- The rule « Andean coverage > 70 % » measures the presence of the Andean wedge as the dominant family, not the danger it
  aimed at, the rampart.
- Tier 1's numbers already said so: the passive coasts (0–6 cells) rise 1.4 km under F180 against 1.6 km under F176; the
  active coasts 4.9 km.
- The direct measure of the rampart is W₁₀ at 256², which had not run.

### 256², young, the 15 seeds (F176 the reference)

| | F176 | F180 |
|---|---|---|
| seeds with ≥ 1 Andean cordillera | 9 / 15 | **15 / 15** |
| collision the only mountain formation | 0 | **0** |
| W₁₀ young, panel (global): 8 / 18 / 19 / 23 / 29 | 0.248 / 0.349 / 0.045 / 0.160 / 0.218 | **0.193 / 0.184 / 0.064 / 0.283 / 0.267** |
| W₁₀ active / passive, panel | 0.20/0.25, 0.38/0.32, 0.13/0.02, 0.24/0.12, 0.25/0.21 | 0.40/0.17, 0.22/0.16, 0.13/0.03, 0.49/0.16, 0.58/0.12 |
| active ≥ 2 × passive (W₁₀), panel seeds | 1 | **4 (8, 19, 23, 29)** |
| mountain % at young, range | 30–83 % | 25–55 % |

- **The cordilleras under F180** (count, total length along the coast, km world): 1–6 per seed, 38–606 km.
  - The six seeds without one under F176 gain one: 23 (269 km), 29 (244), 6 (262), 7 (606), 10 (131), 11 (531).
- **No rampart**: the panel's global W₁₀ stays 0.06–0.28 (under 0.30), and the 0.43 stop does not fire.
  - The active coasts alone reach W₁₀ 0.40–0.58 on seeds 8, 23 and 29: steep Andean fronts, against passive coasts at
    0.12–0.17.
- The mountain share falls under F180 (25–55 % against 30–83 %): the bound also cuts the collisions (s_src −7 % at C1).
- Maps (to validate, condition (v) of the rules):
  - `images/origin_15_F180.png` (B3);
  - `images/formations_15_F180.png` (red: cordillera; pink: Andean not touching an active coast; blue: collision; navy:
    separated collision; orange: arc; purple: other);
  - `images/relief_active_15_F180.png` (the relief with A1's majority-active coasts in red).
    A few straight red lines offshore (seeds 18, 2) are type boundaries on the band's submerged shelf, not shorelines.

### 2 048², young (the transects without a clamp)

- **Seed 18**: the transect spans −4 939 … 3 716 m. The active coast peaks at 3 716 m, the passive one at 2 163 m; the world's
  maximum is 5 632 m (`images/transect_F180_seed18.png`).
- **Seed 7** (the longest new cordillera, 606 km): −5 168 … 4 693 m. The active coast carries a ~4.7 km range; the passive one
  stays at 197 m; the world's maximum is 5 964 m (`images/transect_F180_seed7.png`, the whole map
  `images/map2048_F180_seed7.png`).

### C1, the 15 seeds, F180

- **The fragmentation of the majority-active coasts**: 5–19 segments per seed (8-connected). Their median length is 1–2 cells
  (median over seeds 2.0); the p90 is 2–10 and the longest 3–35 cells. The gaps between neighbouring segments along the coast
  are a median of 1–14 coast cells.
- **Option 5's polarity**: over the subducting pairs' overriding trench cells, **0.55 changes per 10 cells** (median over seeds;
  0–1.27). Per pair the median is 0 (most pairs keep one polarity; the changes are concentrated on a few).
- **The decisions on |Δage| < 0.01**: 11–50 % of option 5's age-decided consumptions over the run (median 35 %).
- **The net Andean uplift along the coasts** (C1's net Δh after the bound and C1's erosion, × the Andean share of the
  Davis-Suppe in the 3-cell disc; no threshold):
  - **active coasts**: median −0.8 km (median over seeds), but a long positive tail: p90 +2.0 to +8.4 km on 9 of 15 seeds;
  - **passive coasts**: median −0.8 km, p90 ≈ 0.
  - In C1's own s (with its erosion) the coast cells mostly lose height; the Andean wedge's net gain sits on a minority of
    active coast cells, the cordilleras of the 256² tier.

### Dahlen 1990, the 2-D section (read, nothing built)

- The wedge is posed in a vertical section perpendicular to the strike.
- p. 60: « Let (x, z) be a system of Cartesian coordinates with x aligned along the top of the wedge and z pointing down
  (Figure 6) ».
- p. 61: the force balance « on an infinitesimal segment of the wedge lying between x and x+dx … per unit length along strike »
  (eq. 5).
- p. 59: the influx ρhV is « per unit length along strike ».

### The reviewer's predictions (diagnostic)

| | prediction | outcome |
|---|---|---|
| D1 | W₁₀ young < 0.30 on every panel seed | held (0.06–0.28) |
| D2 | active ≥ 2 × passive (W₁₀) on ≥ 3 panel seeds | held (4) |
| D3 | condition (v): ≥ 5 of 15 seeds with a cordillera, the separated highland never alone | held (15 / 15; never alone) |
| D4 | median active segment ≤ 3 cells **and** ≥ 1 polarity change per 10 cells | refuted on the second half (segments 2 cells, held; changes 0.55 per 10 cells) |

**Nothing adopted.** The flags stay off. Waiting for the author's feu vert to commit F180 with this diagnostic.

**At commit (2026-10-18), decided**: the author's words: « Feu vert pour commiter le F180 avec son diagnostic, drapeaux coupés. »
- **F180 is committed with its flags off**: `Margin::andean_active = None` and `Margin::ds_mass_bound = None` by default, so
  production and F176–F179 are unchanged bit for bit (`the_f180_controls`; production's C1 hash).
  `Margin::seaward_oo_collision_f180()` and the two viz boxes are only used on request.
- Committed: the bench (`f180_cascade.rs`, with the diagnostic); the active-coast seed (`andean_seeds_active`,
  `overriding_trench_partners`, `andean_feeds_nearest`); the mass bound (`ds_mass_bound.rs`, `C1StepStats::ds_bound`); the
  controls and the viz boxes. `config.rs` and `workspace.rs` stay out.
- **A note of method (the author's; not an error): two defects compensate each other.**
  - The trenches are fragmented: the active segments are 2 cells long (median), and 35 % of the polarities are decided on
    |Δage| < 0.01.
  - Davis-Suppe's wedge spreads in every direction and stitches these pieces back into continuous cordilleras.
  - Correcting the wedge alone (Dahlen's 2-D section, p. 60–61) would give comb-like cordilleras. Both must be corrected
    together, starting with the trenches' coherence.
  - **Queued after the adoption and the deposition.**
