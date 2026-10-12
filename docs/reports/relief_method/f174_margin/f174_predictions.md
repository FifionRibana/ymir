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

# Parts B and C — predictions, written before any measurement

## The reviewer's (hypotheses to check; at least one is false)

- **Q1**: W₁₀ at young drops by ≥ 40 % on ≥ 4 of 5 seeds, but rampart_W still fires on ≥ 1 seed, where margins became active.
- **Q2**: at the end, 10–50 % of the perimeter is active, on ≥ 4 seeds.
- **Q3**: the median shelf width is in the B3 band, on ≥ 4 seeds.
- **Q4**: at B0, the coastline's share within 1 cell of C1's contour drops by ≥ 20 points.
- **Q5**: the coast's dimension rises by ≥ 0.03 but stays under the rocky band on ≥ 4 seeds. A smooth ramp gives a smooth coast.
- **Q6**: the selection finds a continent on 5 of 5 seeds, but with more fallbacks to 20 or 24 plates than in Part A.
- **Q7**: B-hétéro raises the dimension by ≥ 0.05 but moves the classes by > 3 points.
- **Q8**: with the step halved, the drainage density, λ and the sinuosity move by < 5 %.

## Mine

- **The rim is not consumed** (the code reading in `f174_declared.md`): the active share at the end < 5 % on every seed, the
  rim-consumed share ≈ 0. **Q2 refuted**; the stop « no consumption » is reported. Q1's second half (rampart_W where margins
  became active) refuted with it.
- **Q1's first half held**: W₁₀ at young drops ≥ 40 % on ≥ 4 seeds — the O-C Andes wedge no longer seeds at the coast; what
  remains of Davis-Suppe is the O-O fallback dome, set back behind the rim and the ramp. rampart_W fires on 0–1 seed.
- **Q6 held**: the main margin moves the coast ~9–10 cells of 64² inward of the plate boundary (rim 2.5 + ~7 cells of ramp),
  so the continents shrink and the selection needs larger k: found on 4–5 seeds, with more fallbacks than Part A (which needed
  20 plates on 2 seeds). The ×2 variants fail the selection on ≥ 2 seeds.
- **Q3 held**: the median shelf width 6–14 km on ≥ 4 seeds (the pure-Airy expectation is 10.9 km; Davis-Suppe and the erosion
  narrow it a little). The break depth 250–600 m: deeper than the textbook 140 m, since there is no sediment.
- **Q4 refuted**: at B0 the share within 1 cell drops by 5–15 points only — the physics level lowers the coastal land by tens
  of metres at B0, which moves a coast on a ~2° ramp by < 1 cell of 256².
- **Q5**: the dimension stays 0.98–1.02 (under the rocky band, held), but rises by < 0.03 (refuted on the rise).
- **Q7**: B-hétéro raises the dimension by 0.01–0.04 (refuted on 0.05) and moves the classes by > 3 points (held).
- **Q8 held**: density, λ and S₅ within 5 %.
- **The macro**: the mountain class at young falls by 20–50 % relative to Part A (the rampart was much of it); land in 50–62 %
  by the selection.
- **The coast by rock type**: the margin's coastline is mostly cover (low type) at every age; its windows' dimension 0.98–1.02,
  inside the sandy reference (0.95–1.00) or just above, far under the lagoon one.
- **The arms**: fewer at old than in Part A? No — the arms come from valleys worn to 0 m, unchanged by the margin; similar
  counts.
- **Success**: 0 / 5 (the active share fails); the margin does not enter the profile.
- **Cost**: the selection adds 30–120 s CPU per world (~100–300 C1 runs at < 0.5 s); the world itself unchanged.
