# Finding 183 — predictions, written before any measurement

The reference is F182's deposition (`Deposition::declared()`) on F181's default profile, at 256², the panel (8, 18, 19, 23, 29)
at the young (20 Myr) and intermediate (100 Myr) ages. The variants are F182, + A1, + A2, + A3, and + A1 + A2 + A3.

## The reviewer's (hypotheses; at least one is false)

- **Q1**: A1 alone keeps the sea on ≥ 3 of the 10 worlds.
- **Q2**: A1 + A2 + A3 keep the circumnavigation on ≥ 9 of the 10.
- **Q3**: at the intermediate age, the recycling (A3) removes ≥ 20 % of the eroded volume.
- **Q4**: the largest sediment thickness exceeds 5 km under the large deltas.
- **Q5**: the amended mass control holds at every step.

## Mine

- **N0 (tier 0, C1)**: the oceanic plates C1 consumes are few. The area-weighted cumulative consumed fraction of the recyclable
  plates over the history is < 30 % on most seeds, so the recycling's capacity is bounded by C1's slow consumption.
- **N1 (A1 alone)**: the load takes back ~(ρ_b − ρ_w) / ρ_m ≈ 0.2–0.4 of a marine deposit's thickness through the flexure, but
  α = 64 km spreads it. A1 alone delays the sea's disappearance and keeps it on ≤ 2 of the 10 worlds (Q1 refuted).
- **N2 (A2 alone)**: the compaction reduces a km-thick column's thickness by 20–40 %. It delays the filling but keeps the sea on
  ≤ 2 of the 10 worlds.
- **N3 (A3 alone)**: it removes < 20 % of the eroded volume at the intermediate age (Q3 refuted), so it does not save the sea
  alone.
- **N4 (the combination)**: A1 + A2 + A3 keep the sea, at a reduced area, on the young worlds whose ratio is ≤ 1 (8, 18, 19),
  and lose it on ≥ 2 intermediate worlds. Q2 is refuted, and the brief's stop after tier 1 fires.
- **N5 (the thickness)**: the largest thickness is 2–6 km, at the mouths of the largest rivers. Q4 is undecided at ± 1 km.
- **N6 (the mass)**: the amended control holds at < 10⁻¹² relative at every step (Q5 held). F182's porosity term is non-zero
  without A2, < 1 % of the input.
- **N7 (the cost)**: A2's table and A1's one FFT per step add < 10 % to a world. The marine step's cost now scales with the sea
  that survives (up to ~65 000 sub-steps at the intermediate age, about a minute).
