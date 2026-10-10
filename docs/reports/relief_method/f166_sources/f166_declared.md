# F166 — declared before any measurement (2026-10-10)

**The physics level driven by C1's tectonic sources only.** The tectonics supplies the forcing; the erosion is done
once, by the physics level (the author's decision at F165's commit).

**Production is unchanged**: the guard 6 / 6 (field and lakes) runs before and after. C1's outputs stay bit-identical
(the term observer is read-only; a permanent test). The cascade is frozen (p6 · W-tous W1, F164's budgets and k).

**Not blind**:
- F165's readings and images (the subsidence, the drowned margin, T = 1.155 Myr);
- the code read for the inventory below.

Nothing about the per-term decomposition, the sources' drive or its worlds has been measured.

## A1 — the terms of C1's step (`tectonics_c1/time_loop.rs::run_with_closures_observed`)

| # | term | file:line (call) | formula (read from the code) | physical meaning (the code's doc) | class |
|---|---|---|---|---|---|
| 1 | **transport** | time_loop.rs:677 / 691 (`step_upwind_masked`) | first-order upwind advection of `s` and `age` with one velocity per plate; continental cells have v = 0 and a no-flux boundary | plate motion | **tectonic** |
| 2 | **Davis-Suppe** | time_loop.rs:739; davis_suppe/source_term.rs:286–292 | on upper-plate wedge cells: ds/dt = coupling · \|v\| · max(0, h_crit − s) · e^(−d/l_decay); h_crit = h_max·(1 − e^(−d/l_taper)) (C-C) or h_max·e^(−d/l_taper) (O-C, rate ×6, clamped at h_crit) | convergent orogenic thickening to the critical taper | **tectonic** |
| 3 | **equilibrium height** | time_loop.rs:757; equilibrium_height/source_term.rs:51–68 | where s > h_eq = 2.0 (~70 km): s ← max(h_eq, s − k_collapse·(s − h_eq)²·dt), k_collapse = 2.0 | « gravitational sink », the collapse of crust thicker than the Tibet value | **AMBIGUOUS**: a collapse of thick crust (tectonic by the brief's rule), but it removes mass instead of spreading it, ends in a numerical clamp at h_eq, and k_collapse is « calibrated, not literature-derived » (doc §11.1). **Measured both ways (variant).** |
| 4 | **erosion** | time_loop.rs:838; erosion/source_term.rs:116–143 | s ← s − K·A^m·S^n·dt (K 0.001, m 0.5, n 1; ×0.2 on cratons), toward a **floor of 0.2** (oceanic thickness), never below it | stream-power erosion of the crust | **surface** |
| 5 | **subduction** | time_loop.rs:884; subduction/source_term.rs:236–266 | on convergent oceanic cells: s ← s − 0.5·conv·dt; 0.5 of it redistributed as arc mass on continental cells within `arc_distance`; the type flips to continental below a threshold | slab consumption and arc magmatism | **tectonic** |
| 6 | **rift thinning** | time_loop.rs:898; rifting/source_term.rs:135–138 | on divergent continental cells: s ← s − 1.0·div·dt | rift extension | **tectonic** |
| 7 | **accretion merge** | time_loop.rs:919 | `plate_id` only | collision merge | **tectonic** (does not touch s) |
| 8 | **rift split** | time_loop.rs:933; rifting/split.rs:277 | `plate_id`, and `age` = 0 on the rift strip | new plate | **tectonic** (age only: oceanic h_iso through Stein-Stein) |
| — | the craton mask | init only | never recomputed | | — |
| — | **the isostatic datum** | `isostasy.rs::compute_isostasy_inner` (lines 281–292) | the sea level h_sea = h_min + 0.4 · (p92(h_raw) − h_min) is a percentile of the WHOLE field; the land ceiling is fixed (S̃ = 2.0) | the M1 sea-level convention | **other: a mapping, not a term on s**. A land cell of constant s changes altitude when the field's distribution moves. It cannot be removed through s_src; **its drift is measured in A2** and reported, not corrected |

Stein-Stein (time_loop.rs:819) rewrites only the transient altitude inside the step, not `s`.

**The main drive (« sources »)**: terms 1, 2, 3, 5, 6 (7, 8 do not touch s). **The variant « sources − EH »** excludes
term 3. The surface term (4) is excluded in both.

## A2 — the decomposition of F165's history (the témoin)

- The observer records each term's per-cell Δs at every step (read-only).
- **Σ Δs per term** over the run, on:
  - **the interior**: the final land cells (h_iso > 0) farther than 3 cells (Chebyshev) from the sea, a current
    convergent cell or a suture;
  - **the margin**: the other final land cells.
- **Σ Δh per term**:
  - by Airy's linear land slope at the final datum: Δh = Δs · b · peak_m / (land_ceiling − h_sea), where b is the
    cell's buoyancy (craton or not);
  - plus **the datum drift**: the change of a land cell's altitude at constant s between the initial and the final
    h_sea.
- **The question**: the share of the interior's negative Σ Δh that comes from the surface term (4) and the sinks (3),
  against the tectonic terms and the datum.

## A3 — the initial state

- `init_c1_state_phase_2_r7`:
  - continental s = 1.0 (~35 km), oceanic s = 0.2;
  - the cratonic cells × 1.25 (`craton_thickness_ratio`), 15 % of them shields (`craton_shield_fraction`);
  - craton density 2 900 (altitude only).
- **h_iso(t₀) on each seed's land** (deciles, mean): measured and reported.

## B — the mechanism

- **B1, the record**: per cell, S(t) = Σ of the excluded terms' Δs (term 4, plus term 3 in the variant).
  - **s_src(t) = s(t) − S(t)**. On the land this is s(t₀) + Σ tectonic terms in the matter's frame (continental v = 0,
    F165).
  - With nothing excluded, s_src = s − 0 = s exactly.
- **B2, the uplift**:
  - h_src(t) = C1's isostasy (`c1_production_altitude_craton`) applied to s_src(t), with C1's `age`, `plate_type` and
    craton mask at t, normalised like F165;
  - U = Δh_src / Δt over the 30 intervals, upsampled as F165;
  - **everything else is F165's**: the erosion, the talus, the diffusion, MFD p = 6, the rebound (ρc/ρm, α = 64.4 km),
    300 steps, no equilibrium, no deposition;
  - the initial state h_src(t₀) = h_iso(t₀), plus the roughness (F163).
- **B3, T**:
  - calibrated **once, on the témoin, with the main drive** (sources), so its peak at 256² falls in [2 700, 3 000] m;
  - F165's method: trials at log₁₀ T = 5, 6, 7, 8, 9, then bisection on the bracketing pair, ≤ 12 trials; the cost
    reported;
  - **T is then frozen** for the four seeds and the variant. The peaks are not recalibrated or bounded;
  - T is reported against C1's nominal 10–30 Ma.
- **The slope guard** (followed, not blocking): each world's peak, and the slope p50 / p90 of its mountain class at
  2 048² against Corsica's mountain class.
- **B4, the variants and the controls**:
  - **« sources − EH »** (term 3 excluded), the same T;
  - **F164 and F165's history** as references.
    - F165 saved no grids, so they are recomputed: deterministic, and pinned by the steady and history tests.
    - The check is the témoin's F165 history peak (2 852 m) and land (20 376 km²).
  - **The negative control** (a permanent test): with nothing excluded and F165's T, the drive is F165's history bit
    for bit.

## C — the measures (seeds: the témoin, 42, 1, 9; physics 256²; the frozen cascade to 2 048²)

- **Macro**:
  - F165's classifier, unchanged, against Europe;
  - plus the count of mountain regions (8-connected) of more than 100 km², reported without a rule.
- **The peaks and the slope guard.**
- **The land**:
  - the shares lost and gained against F164 at 256²;
  - the interior at ≤ 0 m (land not connected to the ocean).
- **The lakes at 256²** (`cascade::hydro`): count, area, and the share of the lake area on plateau-class cells.
- **The texture**: the Corsica battery, whole and mountain class, at 2 048², with the facets and the 28° walls
  followed; F per class.
- **The coast**: F165's B4 against Corsica.
- **The blocks**: F165's block index (lattice 4 at 256², the 64² blocks) applied to the boundary of {Σ U > 0} of the
  sources' drive at 256².
- **The warp's depressions**:
  - the new depressions of each warp step (F165's footprint rule);
  - the share whose footprint still holds a lake (`hydro` lake cells) in the final 2 048² field.
- **Seed 9's rim**: the cause (file:line), a diagnosis only.
- **The cost**: in CPU seconds per world.
  - The bench reads its own process's total CPU time (all threads, user + kernel) before and after each run, through
    `powershell (Get-Process -Id <pid>).TotalProcessorTime`.
  - A run whose wall time exceeds its CPU time by more than ×5 is flagged as suspended.

## R — the stop rules (the brief's)

- **Macro success**, on ≥ 3 of 4 seeds:
  - mountains within ×1.5 of Europe (16.5–37.2 %);
  - plain + plateau within ×1.5 of Europe (37.5–84.5 %);
  - the land lost against F164 under 15 %;
  - no B4 measure degraded by more than 10 % against F164.
  - The facets, the walls, the peaks and the slopes are followed, not blocking.
- **The « attribution » stop**: Q1 refuted (the subsidence does not come from the surface terms). Then I stop and
  report, with no other drive.
- **Macro success with no plateaus** (< 1 %): reported, and the plateau question goes to C1's work (b).
- **Forbidden**: as F165. No calibrated constant but T, and T is the same for every world.
