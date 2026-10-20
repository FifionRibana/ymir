# Finding 184 — declared before any measurement

The reference is F183's combination (`Deposition::declared()` with `load`, `compaction` and `recycle`, the fluvial pass in one
sub-step) on F181's default profile (`Margin::continent_default()`), at 256². Every change is behind a new flag of
`Deposition`. With the flags off the deposition is F183 bit for bit.

## A1 — the marine domain re-read at each sub-step (`Deposition::marine_live`)

The marine step of F182/F183 takes the step's sea cells as its domain for the whole step. F184 re-reads it at every marine
sub-step:
- **The domain** at a sub-step is the set of cells with h < h_sea at that sub-step (h the current surface, h_sea the step's base
  level). The links (4-neighbours, periodic) join two cells of the domain. Yuan et al. 2019 EPSL eq. 5 holds for h < h_sea only.
- **The sources** (the rivers' outflow, the hillslope's loss to the sea) enter in equal parts at each sub-step, as in F182.
  - A part landing on a cell with h ≥ h_sea walks down the steepest descent (D8, periodic, the current surface) to the first
    cell with h < h_sea. On its way it is placed in that cell up to its room (h_sea − h), and the remainder walks on from there.
    This is the delta's progradation, by a deterministic walk (bounded by the map's cells), not an iteration to a tolerance.
  - A part that reaches a cell with no lower neighbour before finding the sea (a closed depression) stays there: **stranded**,
    reported apart.
- **The diffusive flux** into a domain cell is capped at its room under the sea:
  - with I_r the sum of the sub-step's inflows into r (each computed as in F182: λ Δh × the donor's supply factor), every inflow
    into r is multiplied by c_r = min(1, (h_sea − h_r) / I_r);
  - the capped part stays in the donor;
  - each link's flux is one f64 expression read on both sides, so the mass is conserved;
  - the explicit bound is unchanged (λ ≤ 1/4); the cap only lowers the fluxes.
- **The supply limit** is unchanged (the transport moves only the deposited sediment).
- **The control**: `marine_live = false` is F183 bit for bit.

## A2 — the fluvial stage in k equal sub-steps (`Deposition::fluvial_substeps` = k)

- The step's uplift is split in k equal parts. Each sub-step adds its part, re-reads the land mask, then runs the fluvial stage at
  dt / k: the implicit incision with the rock law (K from E_cum at the step's start), the fluvial pass and the lakes (one pass
  each), and the cover's update.
- **Once per step**, unchanged:
  - the rebound (on the step's summed incision; with A1 of F183 the load's flexure at the step's end);
  - the hillslope and the sea reset;
  - the marine step;
  - the recycling;
  - the compaction;
  - the rock law's E_cum;
  - the random receivers (`receiver_tau`, `None` in the profile).
- **k = 1 is F183 bit for bit.**
- **Measured**: k = 1, 4, 16.

## The convergence criterion (the project's, adopted at F174)

- **The source**:
  - F174, at the commit of Parts B′, C′ and D (ADR): « The time-step question is closed: the physics level's step is not halved;
    the criterion is statistical convergence (density, λ, sinuosity, classes), which holds (F174 A3; B′'s R8 within 2.6 % on
    4 / 4) »;
  - declared in `docs/reports/relief_method/f174_margin/f174_declared.md` (A3; the Q8 / R8 rule, « judged on < 5 % relative »);
  - computed in `tests/f174bc_cascade.rs:2394`.
- **The rule** applied here, between k and 4k, on the same world:
  - the drainage density and λ (`read_cell`) and the sinuosity S₁₀ (`planform`, the 10 km window) each move by **< 5 %
    relative**;
  - the Kapos classes' largest move is **< 1 point** (F174 A3's bar for the classes).
- **k is retained** as the smallest of {1, 4} whose comparison with 4k holds on every panel world (5 seeds × intermediate and
  old).
  - If 4 against 16 fails, k = 16 is not shown converged: **the brief's stop**, with the cost reported.
- **Also reported, not judged**: the fluvial volume deposited and the median altitude (the reviewer's S3).

## The cones' diagnostic

- **The share of the fluvial deposit made under the cap**: in the pass, the deposit of a cell where G / Ñ ≥ 1 is the whole
  incoming flux. It is summed apart (`Pass::capped`).
- **The fluvial thickness on land**: the map of the cumulated fluvial deposit.
- **The lakes born of a deposit**, on the final 256² field:
  - a lake (an 8-connected component of land cells with `filled > z`) whose spill cell has a cumulated fluvial deposit > 1 m;
  - the spill cell is the receiver, outside the lake, of the lake's outlet cell (the routing's D8 on the filled surface).

## The measures (Part B)

- **Tier 1** (the panel, intermediate and old):
  - the variants: without the deposition; F183's combination (k = 1); + A1; + A1 + A2 at k = 4 and k = 16;
  - the quantities:
    - the sea over the history;
    - the circumnavigation **against the same world without the deposition**;
    - the marine deposit above 0 m and above 100 m (a share of the marine Σ), and the stranded share;
    - the median altitude; the Kapos plains and mountain;
    - the lakes (large ≥ 64 cells, dust 1–2 cells, area, born of a deposit);
    - the share under the cap; the largest thickness;
    - the mass controls (flux and stock, F183's); the cost per world;
  - the convergence criterion between k = 1 / 4 and 4 / 16 (+ A1 in both).
- **Tier 2** (the converged configuration):
  - the 15 seeds at the intermediate age, and the panel at the young and old ages;
  - F182's and F183's quantities: the lakes, the shelf (Harris unconverted and ÷ 7.5), Kapos, the coasts, the deltas, (v), W₁₀,
    the circumnavigation, the sediment's shares (land, lakes, sea).
- **Tier 3** (2 048², seed 23, intermediate): the whole map, a crop on a delta, the shelf's transect.

## What a physical deposit should produce (written before; checked, nothing calibrated on it)

On the converged configuration, against the same world without the deposition:
- **E1**: the plains' share does not fall;
- **E2**: the mountain gains ≤ 2 points;
- **E3**: the area of the large closed lakes falls.

A contradicted expectation is diagnosed before any adoption proposal.

## The rules (the brief's)

- **The circumnavigation is lost** only against the same world without the deposition.
- **The stops**:
  - the mass control fails;
  - an iteration to convergence or a solver;
  - k = 16 does not converge;
  - the deposition loses the circumnavigation on ≥ 2 seeds, on the panel or on the 15.
- **The adoption** needs:
  - the mass control;
  - the circumnavigation;
  - (v);
  - E1–E3 held or their departure explained;
  - the author's look at the plates.
