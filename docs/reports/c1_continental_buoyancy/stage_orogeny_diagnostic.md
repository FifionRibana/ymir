# Orogeny diagnostic — why C1 continents are plateaus without mountain chains

The understanding that must precede any relief-chantier fix (the analogue of the
#145 buoyancy diagnostic). Symptom established by the quantitative cut profiles
(`profile_convergent_arc_seed42`); cause established by code-reading + ONE
counterfactual. Both MEASURED, not presumed.

## Symptom (profiles across the convergent arc, seed 42, 3 cuts, consistent)
At a coast convergence: S̃ has a **sharp spike at t=0** (1.4–2.1) but on the
**OCEANIC** cell (plate_type=O) → altitude there = Stein-Stein bathymetry
(deep, ~−0.55) → invisible. The **CONTINENTAL** side (t≥1) is a flat plateau
(S̃≈1.0, altitude≈+0.4). No orogen.

## Code-read verdict — (A): an orogeny mechanism EXISTS but is mis-targeted
Three mechanisms act at a convergence; none builds a continental orogen:

1. **Davis-Suppe** (`apply_davis_suppe_step`) — thickens the "upper plate"
   wedge body toward `h_critical(d) = h_max·(1−exp(−d/l_taper))`. But "upper"
   is defined in `classify_boundaries` by **`is_upper = v_mag_c > v_mag_n`** —
   the strictly **FASTER** plate, **plate_type-BLIND**. In real orogeny the
   overriding (thickening) plate is the CONTINENTAL one (Andes on South
   America). The model targets the faster plate → if the ocean converges
   faster, DS thickens OCEANIC cells (lost to Stein-Stein); the continent gets
   nothing. DS also SKIPS the convergent boundary cells themselves (acts on the
   interior wedge body).
2. **Subduction arc-mass** (`distribute_arc_mass` via `apply_subduction_step`) —
   DOES add S̃ to continental cells (arc volcanism): `arc_efficiency=0.5` of the
   consumed oceanic mass, spread by BFS over `arc_distance=3` onto ALL
   continental cells reached (`per_cell = arc_mass/n`). Real continental
   thickening, but DIFFUSE + weak (small consumed mass, spread wide, capped by
   equilibrium-height, eroded) → a soft broad bump, not a chain.
3. **Advective pile-up at the rigid margin** — the t=0 oceanic S̃ spike is
   NEITHER DS NOR a closure deposit: oceanic crust advects into the rigid
   no-flux continental margin (#145 `step_upwind_masked`) and accumulates on the
   oceanic boundary cell. Subduction consumes it; advection refills faster.

**Causal chain:** convergence → ocean advects into the rigid margin → no-flux
piles S̃ on the OCEANIC boundary cell → `apply_stein_stein_bathymetry`
overwrites that cell's altitude with age-based depth (S̃ **ignored**) → no relief.
Meanwhile DS thickens the FASTER plate (not the continent), subduction arc-mass
is diffuse, accretion only merges plates → continent stays a plateau.

**Altitude branch (confirmed):** `apply_stein_stein_bathymetry` overwrites
**Oceanic** cells (`altitude = −stein_stein_depth(age)`, S̃ unused);
**Continental** cells keep `isostasy(S̃)`. So a high-S̃ OCEANIC cell never
becomes relief; a high-S̃ CONTINENTAL cell does (the dome counterfactual + the
existing cratonic domes prove this path works).

## Counterfactual — DOUBLE verdict (both measured)
Throwaway: force `upper = continental` (the convergent CONTINENTAL cells seed
the DS wedge), re-run seed 42 (env-gated patch in `time_loop`, reverted after).

1. **Lever (A) is REAL** (not an innocenced suspect this time — a true lever):
   the continent thickens and gains relief — cut-2 continental side S̃ 1.13→1.39,
   altitude +0.43→+0.55; S̃ render interior visibly browner; altitude render a
   brown highland massif where the baseline was green/flat. **DS mis-targeting
   (faster vs continental) was a dominant cause of the flat continent.**
2. **Geometry is wrong**: the relief is a **broad interior DOME (Tibet)**, NOT a
   **linear marginal chain (Andes)**. Cause: `h_critical(d)` RISES with distance
   from the boundary + `max_distance≈30` ≈ the whole 64² continent → the wedge
   FILLS the interior instead of BORDERING the margin.

## Chantier macro-border = TWO distinct components
- **Retarget** (`upper → continental`): necessary, proven. NOT 1 line — requires
  threading `plate_type` into `classify_boundaries` (19 callers; via optional
  param or a continental-aware path). The upper-plate criterion
  (`is_upper = v_mag_c > v_mag_n`) is the mis-target.
- **Geometry** (concentrate the wedge near the margin): necessary for a chain,
  distinct from retarget.

## Physical nuance — dome vs ridge are BOTH real (do not "fix" dome uniformly)
Tibet (broad plateau) and the Andes (narrow chain) are two real orogen types by
convergence type:
- **Subduction (ocean under continent, O-C)** → narrow chain along the margin
  (Andes).
- **Collision (continent-continent, C-C)** → wide thickened plateau (Tibet).

So the geometry fix is NOT "always concentrate near the margin" (that would make
Andes everywhere, wrong at collisions) but **"concentrate at SUBDUCTIONS, fill at
COLLISIONS"** — routed by convergence type. The diversity the analysis wants is
HAVING BOTH per tectonics, not replacing a uniform defect (dome) with another
(uniform ridge).

**Does the model distinguish subduction from collision? (code-read, verified):**
- `classify_boundaries` is plate_type-BLIND (Convergent/Divergent/Transform by
  velocity normal only) — does NOT distinguish O-C from C-C.
- BUT the **subduction** closure already keys on the O-C geometry (Filter 2 =
  Oceanic cell + Continental neighbour) → the INFORMATION (plate_type of the two
  plates at the boundary) IS available and already used there.
- There is NO collision (C-C) orogeny: `accretion` only merges plates.
- ⇒ The model HAS the ingredients to route by convergence type; DS simply does
  not consume them. A geometry fix can branch on the boundary's plate-type pair:
  narrow marginal wedge at O-C, broad fill at C-C.

## Injection points (surface only)
- Retarget: the `is_upper` criterion in `classify_boundaries` (thread `plate_type`).
- Geometry by type: branch the DS wedge profile/extent (`l_taper`/`max_distance`,
  or which side it deposits) on the O-C vs C-C plate-type pair at the boundary.
- Existing continental-deposit hook: `distribute_arc_mass` (subduction) already
  puts S̃ on the continental side — concentrating it would make an O-C chain.

## NOT done here
The fix is not coded. This is the diagnostic. The chantier (retarget + type-routed
geometry) scopes from this. The throwaway counterfactual patch was reverted.
