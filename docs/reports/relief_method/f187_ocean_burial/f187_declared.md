# Finding 187 — declared before any measure

The ocean defined by connectivity (A, flag `Sea::connected`) and the burial bounded by the 256² surface (B), before the
deposition's adoption. Written before the measures; amended only above an « Amendment » line with its reason.

## Part 0 — the inventory (read in the code at 5e40634)

### The sea in the physics level (a history level, `hold_sea = false`, the sea at `sea`, 0 m)

Every test is **per cell, with no connectivity**: a cell below the sea in the middle of the land is the sea.

| what | where | the test |
|---|---|---|
| the incision's base nodes | `erosion/stream_power.rs:670` | `d == DIR_NONE \|\| field ≤ sea_level` → fixed node |
| the incision's routing (priority flood seeds, D8) | `terrain/flow.rs:198` (`compute_flow`), called at `stream_power.rs:651` | `is_ocean = h ≤ sea_level`; the flood starts from every such cell, so a below-sea basin is never filled |
| the incision's MFD area | `terrain/flow.rs:902` | `is_ocean = filled ≤ sea_level` |
| the land mask of the step | `cascade.rs:642`, `:698` | `land = z > sea` |
| the hillslope's base level | `cascade.rs:841-846` | every `!land` cell keeps its pre-hillslope value |
| the fluvial pass (deposition) | `cascade.rs:716` (`sea_mask = !land`), `cascade/deposition.rs:533` (`compute_flow` at `sea`), `:537` (receiver fixed if `surface ≤ sea`), `:595` (a lake needs `surface > sea`), `:628` (an outlet is the sea if `sea_mask`) | per cell |
| the marine domain (F184 A1) | `cascade/deposition.rs:186` | `dom = hz < sea`, re-read at each sub-step |
| the marine sources' walk | `cascade/deposition.rs:145` | a source fills the first cell with `hz < sea` it reaches |

Two comparisons coexist: ≤ (the base nodes, the routing, the pass) and < (the marine domain, the walk).

### The sea in the amplification chain

- `amplify.rs:491-526`, `ocean_mask_t`: **the largest 4-connected component (periodic on the torus) of the cells with
  z ≤ 0**; ties go to the component holding the smallest terrain index.
- `amplify.rs:875-884`: the level's ocean; the D5 cells are the non-ocean cells ≤ 0 (and the previous level's D5).
- `amplify.rs:885-890`: the previous level's working field puts **the ocean at 0 m** (the cause of F186's family ii).
- `amplify.rs:903` → `amplify.rs:442`, `lift_interior`: **every non-ocean cell ≤ 0.5 m is lifted to the land floor, 0.5 m**
  (the cause of family i).

### The difference with A1, and A4's alignment (declared)

- A1 is the chain's definition (the largest 4-connected periodic component) with a **strict** comparison, h < h_sea, where the
  chain uses z ≤ 0.
- They differ only on the cells exactly at the sea. Those are common in the physics level's output: the marine step fills a
  cell to the sea exactly (`deposition.rs:146`).
- **A4**: behind the flag, the chain takes A1's strict comparison (`ocean_mask_t` with z < 0). Flag off, the chain is unchanged.
- The tie rule (the smallest terrain index) is the chain's, kept for A1.

### The measure of Part 0 (on the panel, intermediate, no deposition, flag off)

- At each step of the physics level, the components of `h < h_sea` other than the largest one (4-connected, periodic): their
  number, their area, their greatest depth below the sea.
- Reported at the history's end and over time: the steps with at least one, the maximum over the history, the first step.
- **How**: an observer on the level's steps, `Drive::observe` (`None` everywhere else, so production and every earlier bench
  are unchanged), reached through a new entry point of the history level. It reads the field; it changes nothing.

## Part A — the ocean by connectivity (`Sea::connected`)

- **A1**: at each step (and at each marine sub-step), the ocean is the largest 4-connected periodic component of `h < h_sea`.
  It has no parameter.
- **A2**: only the ocean's cells are:
  - the incision's base nodes;
  - the marine domain (the sources, the walk, the diffusion);
  - the cells held at their pre-hillslope value.
- **A3**: the other below-sea cells are closed depressions under the lakes' rules:
  - the routing fills them to their sill (the priority flood starts from the ocean only);
  - the deposition fills them by volume, as a lake;
  - they receive their catchment's flux.
- **Where it goes**:
  - the incision and its routing get an ocean mask as an argument (a new entry point; the old one passes `None`);
  - `compute_flow` and the MFD area get the same mask through new variants;
  - the pass gets the mask for its receivers, its lakes and its outlets;
  - the marine step gets the connected domain at each sub-step;
  - `StreamPowerConfig` gets no new field, so the cache digests stay as they are (F186's lesson).
- **A5's controls**:
  - the flag off = F186 bit for bit;
  - on a world where every below-sea cell is connected, the same result with and without the flag;
  - the mass within 10⁻⁹;
  - production's hash bit for bit.

## Part B — the burial bounded by the 256² surface

- **B1**: on each covered 256² cell the volume is max(0, z_256 − b̄_chain) × n, where b̄_chain is the mean of the fine basement
  the chain left (after `lift_interior` and the ocean at 0 m). L follows as at F186: the exact sort, the bilinear
  interpolation, z = max(fine basement, L).
- **B2**: so the fill never exceeds the 256² surface.
  - Shown: f(L) = Σ max(0, L − b) ≥ n (L − b̄), so f(z_256) ≥ n (z_256 − b̄), the target, hence L ≤ z_256 on each cell.
  - The bilinear L then stays below the bilinear z_256.
- **B3's controls**:
  - H = 0: no change;
  - a chain that leaves the coarse basement (b̄ = z_256 − H) gives F186's volume;
  - the highest filled cell against the interpolated 256² surface, measured (bilinear, the interpolation the burial uses).
- **Where it goes**: `Transmit::bury_bounded`; F186's `bury` stays as it is.

## Part C — the measures (on mains; the budget announced before each tier; < 45 min of wall time each)

- **Tier 1**, the panel at 256², intermediate:
  - the runs: without the deposition (F186); without the deposition with A; with the deposition with A (the per-cell marine
    step, F186's configuration A);
  - the quantities:
    - the unconnected below-sea depressions (number, area, fate);
    - the sea, and the circumnavigation against the same world without the deposition, with A;
    - (v) and W₁₀; the plains, the mountain, the median altitude;
    - the large lakes (E3, reported); the mass control; the cost;
  - family (iii): on seed 23, the rock at the cells F186 found raised, with and without A.
- **Tier 2**, the chain to 2 048², seed 23, with A and B:
  - the filled basins' greatest height against the surrounding land;
  - the fine volume against the coarse;
  - the plates (the whole map; the lake, delta and plateau crops; without | with).

## The rules (the author's, written before)

- **The stops**:
  - the mass control fails;
  - a solver or an iteration to convergence;
  - the circumnavigation lost to the deposition on ≥ 2 seeds;
  - A moves the world without the deposition by more than 5 points of plains on one seed (stop and report before going on).
- **The adoption of the deposition** (with A and B) only if:
  - the mass, the circumnavigation and (v) hold;
  - no burial bump is left;
  - the author validates the plates.
- **Then** F188: the chain's budget, with the per-cell marine step's acceleration validated against the reference.
- **Forbidden**: as at F165–F186.
