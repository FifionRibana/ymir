# Finding 140 — predictions, written 2026-10-02 BEFORE any code or measurement of this round

**Reading declaration, unfavourable.** Non-blind on Findings 73–139, the reviewer's predictions, and in particular:
- F139-R's table: p = 0.5, 14 D8 lakes: 1 308.5 → 1 202.7 → 521.3 → 103.0 → 19.0 km²;
- F138-T4: 11 tagged head falls at today's levels with m = 10;
- F139-C2 / C3: S_req ≤ 6.6°, "steepen" leaves no shortage fall.

**Not blind on how I will build it.** The design below is mine and is written before coding. The construction does
not know the final world's lakes, so it must derive L_bed (the bed sill) and the outlet area itself. That is why
I expect it to differ from the benches' tables.

## The design I expect to build (stated so the predictions can be read against it)

- **The bodies**: the construction's present lakes: input lakes and closed-depression components ≥ 1 km².
- **L_in**: an input lake's fill level, or a depression's spill (the merged spill, F139-C5).
- **L_bed**: the minimax height from the bed's lowest cell to the outlet.
- **L(r)** sets the lake's base in the χ walk (instead of today's L_in), so rivers in are based on the lowered
  lake and a drained lake (r = 2) bases on L_floor.
- **The gorge** raises the outlet path's sample floors: `zf = max(law, z_gorge)`.
- **The ring** is clamped ≥ L(r).

## Predictions

- **G-inert**: holds to the bit (with the reviewer).
- **×1, p = 0.5**: the 14 D8 lakes' total area is within **20 %** of 521.3 km², **not within 5 %** (against the
  reviewer). The construction's own L_bed and outlet area differ from the tables', which read the final world.
- **G-levels**: at r ≤ 1, ≥ 70 % of the present lakes sit at L(r) ± 1 m. It fails more beyond r = 1, because the bed
  sill is cut only along the trunk through the lake.
- **G-rim**: holds on ≥ 90 % of the lakes with r ≤ 1. The light pass may still lower a few ring cells.
- **G-ring**: fails on ≥ 1 lake (with the reviewer). The clamp leaves ring cells above an outside cell another valley
  carved: > 28°.
- **G-slope100**: holds on ≥ 90 % of the untagged descents; at most 2 fail (at the gorge / law junction or under
  the light pass).
- **G-drained**: holds where measurable (×1.4: ~10 lakes drained per F139-R).
- **G-area**: holds (non-increasing).
- **G-tag** at ×1, p = 0.5: **6–12** tagged falls (overlapping the reviewer's 8–14); no gorge cell over H_f.
- **Non-regression**: canyons 0, coast unchanged, θ on the carved links within OFF's CI (with the reviewer).
- **Cost**: < 5 s added per world (with the reviewer).

## Meta
- Disagreements: the 5 % bound at ×1 and the low end of G-tag. **At least one of my predictions is wrong.**
