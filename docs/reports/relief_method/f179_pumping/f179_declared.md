# Finding 179 — declared before any measurement

**Part A measures the Davis-Suppe pumping; it bounds nothing (A3 is a separate recomputation that feeds nothing back). Part B
changes one closure — where the Andean wedge starts (at the arc) — behind the continent profile's flag, with the seed skipped
like a boundary seed; no solver; production unchanged. No adoption in this finding.**

## Part 0

- **The `f165` test**: nothing was removed or renamed. `compose` and `load_grid` exist in `tests/common/cascade_bench.rs`.
  - `tests/f165/main.rs` is a submodule of `f165_cascade.rs` (`#[path = "f165/main.rs"] mod f165_main;`, its imports
    through `super::`). But cargo also auto-discovers every `tests/<dir>/main.rs` as a test target of its own, « f165 ».
  - Compiled alone, its `super::` paths do not resolve: hence the 51 errors.
  - It has been there since it was added (`ee841a8`, F165, 2026-10-10).
  - **The repair, test-only**: the file is renamed `tests/f165/f165_main.rs` and the `#[path]` updated. Every `ymir-core`
    test now compiles.
- **Production does not use the sources mode.**
  - `run_hd`'s eroded product runs C1 with `run_with_closures` (`tectonics_c1/cached_product.rs:281`), with no term
    excluded;
  - it then upscales the final state (`:305`, `upscale_from_c1_with_progress(&state, …)`), whose coarse field is
    `c1_coarse_normalized_altitude(state, …)` (`tectonics_c1/production_upscale.rs:455`) → `c1_coarse_raw_altitude(state,
    …)` (`:293`): C1's own `s`, erosion included.
  - The sources (`run_c1_sources*`, `cascade/history.rs`) serve only the cascade, its benches and the viz's cascade tab. The
    pumping does not reach production.

## A2 — the kinematic supply, from the sources

- **Dahlen 1990** (Annu. Rev. Earth Planet. Sci. 18, 55–99):
  - p. 59: « The mass flux per unit length along strike into the toe of the wedge is ρhV », and the wedge's growth is the mass
    conservation law ∂/∂t[½ ρ W² tan(α+β)] = ρhV (eq. 1): the wedge's cross-section grows by exactly hV per unit time;
  - p. 60, eq. 4: « An eroding wedge will attain a dynamic steady state when the accretionary influx rate of fresh material
    into the toe is balanced by the erosive efflux … ē W sec(α+β) ≈ ē W = hV »;
  - Taiwan, p. 84: V = 70 km/Myr, h = 7 km, W = 90 km, ē = 5.5 km/Myr.
- **Willett, Slingerland & Hovius 2001** (Am. J. Sci. 301, 455–485):
  - p. 459: « At steady state, for accretion of a constant thickness of new material, h, conservation of mass requires that
    U = V_c h / L » (with the velocity field u = U, v = V_c x / L, eq. 1, p. 458);
  - p. 456: the tectonic model « is purely kinematic without feedback; … independent of erosion rate or total erosion »;
  - p. 455: the balance is « not as simple as an equating of surface uplift rates and erosion rates on a pointwise basis ».
- **The formula retained**: the mean thickening rate a wedge can receive is the influx over its width, **U = V_c h / L**
  (Willett p. 459; the same as Dahlen's ē = hV / W, p. 60). Per cell and per C1 step:

  **Δs_kin = (v_rel · n̂) · dt · 64 · s_in / L**

  - **V_c** ↔ `v_rel · n̂` (> 0) on the convergent edge of the wedge's seed, C1's domain units per time; × dt × 64 = cells per
    step;
  - **h** ↔ `s_in`, the s of the incoming cell across that edge: the other plate's cell for a boundary seed, and the sinking
    partner of the trench for an Andean seed;
  - **L (W)** ↔ the wedge width = Davis-Suppe's reach, `max_distance` = **30 cells** at 64² (the brief's choice: the reach in
    cells Davis-Suppe uses);
  - each Davis-Suppe cell takes the supply of its nearest seed (the same metric as `wedge_distance_intra_plate_typed`: 8-
    connected, √2 diagonals, within the plate).
- **Terms of the sources with no C1 equivalent, and what is done** (no value outside the sources):
  - ρ (constant, it cancels);
  - the taper α + β and sec(α+β) ≈ 1 (Dahlen's narrow-taper approximation, used as stated);
  - **the accreted fraction**: Dahlen's and Willett's h is the thickness of material accreted, i.e. scraped off. C1 has no
    sediment layer, so h is taken as **the incoming cell's whole s** — an upper bound of the accreted thickness (the bound is
    therefore generous);
  - the split between frontal accretion and underplating (Willett) is not needed: the bound counts the total influx;
  - the erosion rate ē is not used: the bound is kinematic, which is the point.
- **The reviewer's reading checked against the sources**: both articles make the influx kinematic (hV, V_c h / L), and Dahlen
  makes the steady-state width the adjusting variable (W = hV / ē). They **support** « the flux is set by the kinematics and
  the relief adjusts »; Willett adds that the balance is not pointwise. No contradiction to report.

## A1 and A3 — the definitions

- **The families of Davis-Suppe cells** (per step, by the nearest seed, as C1's typed wedge does):
  - **collision**: a non-O-C seed on a continental cell with a continental neighbour of another plate;
  - **O-C boundary**: the existing O-C override seeds;
  - **Andean**: the F178 / F179 seeds;
  - **other**: any other seed (the velocity fallback).
  A cell's family over the run is the one that gave it the most Davis-Suppe.
- **A1, the pumping factor** of a cell = Σ Davis-Suppe increments over the run / the net change of C1's s (s_end − s₀) on that
  cell. It is reported for the cells with a positive net change (median, p90, max per family and profile); the cells with a
  net ≤ 0 are counted apart, since their factor is unbounded. The family's aggregate ratio ΣΣ / ΣΣ is also given.
- **A2**: per Davis-Suppe cell and step with an increment: the share of steps where the increment exceeds Δs_kin, and the
  cumulative excess Σ max(0, inc − Δs_kin), in s and in km of crust (× 35), by family and profile.
- **A3**: per cell, s_src,bounded = s_src − Σ max(0, inc − Δs_kin), with s_src = s₀ + Σ every term but C1's erosion (F178's
  definition). Its distribution against s_src on F176's collision belts (the collision-family cells) and on the active coasts
  (the A1-active 0–3 cell band, F178's definition). An estimate of what F180 would change, nothing more.
- Profiles: F175 (`seaward_oo`), F176 (`seaward_oo_collision`), F178 (`seaward_oo_collision_andean`) and F179 (the seed at the
  arc); the 15 seeds of F176's selection (F177's cached profiles).

## Part B — the Andean seed at the arc (declared)

- **B1, d_arc = 4 cells**: the lower bound of the middle 50 % of arc-trench distances, 180 km (Syracuse & Abers 2006 p. 10), ÷
  7.5 ÷ 6.25 = 3.84 → the first whole depth beyond.
- **B2, the seed**:
  - for each active trench on which P overrides (F178's rule, `overriding_trench_cells`), grouped into 8-connected segments of
    one plate, a walk (4-connected) from the segment through **P's cells of any type**;
  - d_front = the depth of the first continental cell of P reached; d = max(d_arc, d_front);
  - the seeds are P's continental cells at exactly depth d, if d ≤ D_gap = 10;
  - no reading of the sea level.
- **B3, the seed is skipped like a boundary seed**: the Davis-Suppe step receives a skip mask (the seeds) and gives them no
  source, so the wedge starts at d = 1 of the seed. Davis-Suppe's formula is unchanged (the mask only extends the existing
  boundary skip to these seeds).
- **B4, the controls (permanent)**: the seed off is F176 bit for bit; production's hash (`f175_c1_hash`); a reach of 0 has no
  effect.
- **B5**: nothing else changes (C1's step is out of scope).
- **The measures (C1 only, the 15 seeds)**, for F176, F178 and F179:
  - the land at the end of C1 and r*;
  - the origin of the land added over F176 (in the band at t₀ or not, the type at t₀);
  - the uplift by origin by distance band (F178's definitions);
  - the A1 share with both rules;
  - the Andean seeds' pumping factor;
  - **the seeds' position**: the distribution of d, and the share of seeds on cells emerged at t₀ (Airy, s₀ ≥ 0.915) —
    measured, not a criterion.

## The budget (announced; resumable per world)

C1 only, the 15 seeds × 4 profiles (F175, F176, F178, F179) with the pumping observer (three typed wedges and the supply per
step), parallel by seed: **~10–15 min wall**.

## The rules (written before)

- No adoption. Stops: any other closure change, or a solver; any bound applied (A3 feeds nothing back).
- If A1 shows that the cascade's relief accepted since F166 rests mostly on pumping (P3 true on the collisions), it goes at the
  head of the report.
- Forbidden: as in F165–F178. No constant calibrated (d_arc and D_gap from Syracuse & Abers; L is Davis-Suppe's reach; h the
  incoming s).
