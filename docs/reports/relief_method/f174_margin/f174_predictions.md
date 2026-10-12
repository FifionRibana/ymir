# Finding 174 — predictions, written before any measurement

## The reviewer's (hypotheses to check; at least one is false)

- **P1**: with no margin, ≥ 90 % of the coastline lies within 1 cell of 256² of C1's contour, on ≥ 4 of 5 seeds.
- **P2**: the real low coasts have a dimension in 1.02–1.08, the rocky ones in 1.10–1.20.
- **P3**: with the step halved, max |Δz| exceeds 200 m on ≥ 3 seeds, but the classes move by less than a point.
- **P4**: early in C1, more than 60 % of the perimeter is convergent (the mean over the first 50 steps). The immobile
  continents are the probable cause.
- **P5**: with the margin, the shelf's median width falls in the expected band (B2), and the coastline's share within 1 cell
  of C1's contour drops by ≥ 20 points. Yet the dimension stays under the rocky band: a smooth ramp gives a smooth coast.
- **P6**: B-hétéro raises the dimension by ≥ 0.03 but moves the classes by more than 3 points.
- **P7**: W₁₀ at young moves by less than 20 % with the margin.

P5–P7 belong to Parts B and C, which are not run this round (the sources are unreadable, `f174_declared.md`); they stay open.

## Mine

- **Spratt**: the re-zeroed curve reads 0 at 0 ka; nothing else moves (F173 is not re-run).
- **A1 (P1)**:
  - At B0, P1 **holds**: ≥ 90 % within 1 cell on ≥ 4 / 5 seeds. The physics level moves the coast by less than a 256² cell at
    B0 (it erodes ~2 km of relief over 5.6 Myr, the coast at 0 m stays where the uplift put it), and the warp moves it by a
    fraction of a cell.
  - At old, P1 **fails**: 60–90 % within 1 cell; the coastal lowlands, graded to sea level over 300 Myr, are partly cut
    below 0 m and the coast retreats by 1–5 cells in places. So P1 is mixed and the stop rule does not fire.
  - ≥ 98 % within 5 cells at both ages: the coast's shape at > 8 km is C1's.
  - The physics level's own share within 1 cell (256²) is higher than the 2 048² share by 2–10 points (the chain adds the
    rest).
- **A2 (P2)**:
  - The rule calls 7 of the 8 windows right; North Brittany is called low (its p90 near the coast < 200 m).
  - Landes: dimension ≤ 1.02 and length ratio ≤ 1.05 (a near-straight line). Jutland and the Golfe du Lion, with their lagoons
    and lidos, have a dimension > 1.08, so the low band reaches past P2's 1.08: **P2's low range refuted** on its upper end.
  - Rocky: Galicia and Dalmatia (rías, islands) above 1.20, so **P2's rocky range refuted** on its upper end; the rocky band
    ~1.10–1.30.
  - The two bands overlap in dimension; the length ratio separates them better (low ≤ 1.3 except the lagoon coasts, rocky
    ≥ 1.3).
  - Our worlds: their coast windows are mostly low by the rule at old, mostly rocky at B0; their dimension (~1.0) sits under
    both bands' lower ends except the Landes, i.e. our coasts are **as smooth as the smoothest low coast**, and smoother
    than every rocky one.
- **A3 (P3)**: max |Δz| > 200 m on ≥ 3 seeds (local channel-head and knickpoint shifts, as F173's 593 m); p99 < 50 m; the
  classes move < 1 point; order-≥ 3 channels displaced 5–20 %. **P3 holds.**
- **A4 (P4)**:
  - The convergent share over the first 50 steps is 35–60 % (under P4's 60 % on the panel's mean); divergent 20–35 %,
    transform the rest; **no boundary ≈ 0 %**.
  - The same-plate share ≈ 0: the continent profile makes continents out of whole plates, so every coast is a plate
    boundary from t₀, and since `plate_id` is never advected the boundaries never move; a passive margin cannot form.
  - So the zero continental velocity is **not** the main cause: the classification reads the plates' velocities, not the
    cells'. The cause is the plate geometry (whole-plate continents, static plates). P4's number **refuted** (or barely held),
    its cause **refuted**.
- **The stop**: Part B stops anyway (the sources); the round ends with A0–A4 and the request for the articles.
- **Cost**: A1–A2 add nothing to the world (measures only); A3 doubles the physics level (+100 % of its CPU, ~+15–25 % of a
  world).
