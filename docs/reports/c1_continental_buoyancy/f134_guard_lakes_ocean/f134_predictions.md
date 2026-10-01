# Finding 134 — predictions, written 2026-10-01 BEFORE any grep, reading or measurement of this round

**Reading declaration, unfavourable.** Non-blind on Findings 73–133, F133v and the reviewer's predictions.
**Non-blind on Δ**: F133-F already printed the extended residual's |Δz| p50 6.50 / p90 68.80 / max 631.2 m (témoin).
Facts held: the conditioned z differs OFF/ON on 44.3 M cells (land ≈ 11 M); the construction only LOWERS land
cells; `run_hd` breaches with `breach_monotone_protected` (active crater bowls kept closed) and then adds the C-2
crater lakes; the bench assembly uses plain `breach_monotone` and no crater pass.

## 0.3 — the bench with the crater pass
- **P-0.3**: the crater pass ALONE does not give the bench the crater lake: the bench's plain breach opens the active
  crater's bowl, so the pass finds it dry. The bench needs the PROTECTED breach (with the crater mask) too; with both,
  the lake masks are identical lake by lake, OFF and ON (agreeing with the reviewer's "identity", with that caveat).

## O — the ocean
- **P-O**: the first OFF/ON difference on a sea cell appears at the BATHYMETRY stage (the last of the upscale), through
  a statistic computed on the whole grid; the light pass and the construction leave sea cells bit-identical. If the
  statistic is a global mean depth, z_ON/z_OFF is constant within 1e-3 on > 90 % of the sea cells (agreeing with the
  reviewer in kind, at 90 % not 95 %).

## D — the p90 of 1 177 m upstream of lake 2 (disagreeing with the reviewer)
- **P-D**: it was read on the CONDITIONED field (no water), on cells OUTSIDE lake 2's ON footprint (the F133v instrument
  excluded lm_on cells): it is NOT the lake. Sign: z_ON > z_OFF (OFF's χ ran through the lake to the sea and carved
  the trunk floors deep; ON lays them from the lake's col). The large values sit on cells CARVED by OFF's
  construction near the lake's outlet trunk; their order is (lake level − OFF's floor), comparable to the col's
  altitude above the OFF trunk floor, not to the lake's depth.

## T
- **P-T**: θ on the carved links upstream of a based lake ≤ 0.43; downstream 0.47–0.50; untouched basins 0.49 ± 0.01. The
  carved cells > 1 m off their law sit within ~2 W of a based lake's shore (the shore's base jump).

## U (disagreeing with the reviewer)
- **P-U**: the added under-law cells are long reaches, not delta places: median distance to a based shore 5–15 km,
  < 30 % within 2 km, law − terrain p50 > 30 m.

## Δ
- **P-Δ (non-blind on the numbers)**: p50 6.5 / p90 68.8 / max 631 m (F133-F); the function is `carve` (its nearest-sample
  propagation, `who[]`), the SAME code as Finding 121/127's "laid from the NEAREST skeleton sample". The reviewer's
  "p90 < 2 m" is refuted by the existing output.

## Meta
Disagreements: P-0.3 (the protected breach), P-D, P-U. At least two are wrong.
