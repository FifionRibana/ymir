# Finding 182 — predictions, written before any measurement

The reference is F181's default (the seaward margin, option 5, the merge at collision, the Andean seed at the arc and
Davis-Suppe's mass bound) at 256², without deposition. « With » is `Deposition::declared()` (G = 1.55, the lakes, the marine
silt K_M = 400 m²/yr, φ = 0.63).

## The reviewer's (the brief's, verbatim in substance)

- **P1**: the large lakes fall by 50 % (with the declared class: ≥ 64 cells).
- **P2**: the shelf width falls in the Harris ÷ 7.5 band on ≥ 10 of the 15 seeds. (Kept as written. Since the author's
  complement it is reported against both Harris unconverted and ÷ 7.5, with no criterion.)
- **P3**: the plains gain 5 points and the median altitude 20 m.
- **P4**: the mass control holds.
- **P5**: (v) holds on 15 / 15, and W₁₀ changes by < 0.03.
- **P6**: the coast's dimension changes by < 0.03.
- **P7** (the complement): the gap between the explicit splitting and Yuan's method is larger at the longest step (1.47 Myr)
  than at the shortest (18.7 kyr).

## Mine

- **M1 (the mass)**: the control holds at < 10⁻¹² relative on every world (f64 bookkeeping; the marine solid measured on the
  field). The f32 field residual is ~10⁻⁶–10⁻⁵ relative.
- **M2 (the lakes)**: the lakes on the final field fall in count AND in area. The large lakes fall by more than 50 %, because
  a step's fill reaches the spill wherever the inflow exceeds the capacity. The dust falls less (1–2-cell pits re-form each
  step from the incision and the roughness).
- **M3 (where the sediment goes)**: at the intermediate age the sea takes most of it (> 60 %), the lakes 10–30 %, the fluvial
  deposit on land < 15 % (G / Ñ falls fast downstream).
- **M4 (the relief)**: the median altitude rises by < 20 m and the plains by < 3 points. The fluvial deposit is small on land;
  the lakes are a small share of the land.
- **M5 (the coast)**: the deltas add < 2 % of land; the coast's dimension changes by < 0.03 (P6 holds).
- **M6 (the shelf)**: the deposition widens the passive shelf by a few cells at most. Its median stays well under Harris
  unconverted (88 km), near or under the ÷ 7.5 band's top (16.4 km): with K_M = 400 m²/yr over 100 Myr the diffusion length is
  ~200 km, but the supply is small against the ocean's volume.
- **M7 (the protection)**: (v) holds on 15 / 15; W₁₀ changes by < 0.03; the circumnavigation is lost on at most one more seed.
- **M8 (the reference)**: Yuan's iterative coupling converges on most steps in < 12 iterates (Yuan's Fig. 2 range). The
  deposit maps correlate at r > 0.8 and the land deposit mass differs by < 30 %. The median altitude differs by < 5 m.
  P7 holds: the gap grows with the step's duration.
- **M9 (the cost)**: the marine step costs ~10–40 s per world at 100 Myr and ~3× that at 300 Myr; the fluvial pass adds
  < 20 % to the level's CPU.
