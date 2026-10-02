# F135 — instruments, declared BEFORE the measurements (2026-10-02)

**The realigned témoin** is the viz state "C2 /10 col (défaut)" itself:
- the field: `ValleyConstruction::new(F121_AGE_K, Some(0.1))`, slope floor 0.024, `Knobs::passes(2)`, canonical
  framing (guard field `a8d2d538d692c2f0`);
- the assembly: `run_hd`'s tail (`common::viz_hd_lakes_on`: protected breach, climate 45° / 40°, H-1
  infiltration, assembly, crater pass; guard lakes `b7e502b2c1e126c3/9960d2b71e538692`).
- "F132" adds `lake_base: InputLakes`; "extended" adds `InputLakesAndBasins`.
- Every world prints its field hash and its lake fingerprint against the guard's reference.

## W1 — the configuration differences
- **Field**: `serde_json` of `bench_cfg(témoin knobs)` against `bench_cfg(viz knobs)`, every key that
  differs, recursively; and the two eroded digests (the viz's must read `cf1d4539c7afe588`).
- **Drainage config**: `dcfg()` against `viz_dcfg()`, the same way.
- **Assembly steps** (breach, crater pass, infiltration field, climate): read from the code, cited by line.

## W3 — F133's table on the realigned témoin (OFF / F132 / extended)
F133's instruments are kept, on the realigned world and its assembly:
- **canyons**: `f95_criteria` + `over_dug_depression` (the drainage config is `viz_dcfg()`; the breached
  field is the CONDITIONED one);
- **coast**: spurs near a coastal wall, per km, with F126's wall mask on PRE's carve;
- **R8**: terrain (window 16) and network (chord 8, on the run_hd-tail watercourses);
- **lakes** and their km²;
- **relief p50**: `Crit::p50`;
- **Δz between flat resolutions**: F132 / F133-F's A/B patches on S1, cells outside the present lakes
  (input lakes ∪ closed depressions) with z_A ≠ z_B, |Δz| p50 / p90 / max. OFF is also given outside the
  input lakes alone.

## M — trunk metrics without the lakes
- **Lake mask**: the UNION of the three worlds' final lake maps (run_hd tail, after the crater pass), so
  that every state is read on the same cells. The world's own mask is given beside it for (iii).
- **θ (i) / (ii) / (iii)**: F133's `theta_sa_ci` instrument (the OFF skeleton's trunk links ≥ 10 km²) on:
  - (i) the built field;
  - (ii) its breach: the run_hd breach, protected, on the built world's craters;
  - (iii) the world.

  The fit and the bootstrap are F133's, on the link list; a link is dropped if either end lies in the
  mask. Without a mask it must reproduce F133's values on the realigned world.
- **θ on the carved links**: `theta_laid_ci`'s population, the same mask rule.
- **Under-law share**: trunk ≥ 10 km² cells not in the mask, uncarved, with z_S1 < floor_m.

## K — the cuts through the cols
- **Items**:
  - every input lake (S1's `c1_drainage_windowed`, default config, as F133);
  - every `basin_base` closed-depression component (on S1's breached field);
  - (K4) every unbased final lake of the extended world.
- **Level and col**: the lowest cell of the item's outer 8-ring, on S1 for lakes and on the breached S1 for
  depressions; the level is that cell's z.
- **The built floor**: `carve(S1)`, OFF and extended.
  - **cut = level − (the lowest built z on the outer ring)**; the cut at the col cell itself is beside it.
  - Counts of cuts > 10, > 100 and > 500 m per state.
  - Whether a W3 canyon floor lies inside the item.
- **K2**: the path runs on the OFF skeleton (S1):
  - downstream of the col along `sk.direction` (60 km at most);
  - upstream through the lake by the max-area donor (40 km at most).

  It shows z_S1, the OFF and extended CONDITIONED fields, the extended lake level over its span, and the col
  marked. SVG, one per lake (1, 2, 11 = F133v / F134 ids, found by centre).
- **K3**: lake 2's attributed upstream cells (F133v's attribution).
  - Distance = the OFF drainage's D8 path length to the col (within 3 cells).
  - 2 km bins to 30 km: count and p50 Δz (conditioned), on all cells and on cells carved in both worlds.
- **K4**: each extended cut > 10 m is labelled:
  - a base on an absent lake (an item whose cells are < 50 % under the extended final lake map);
  - an unbased final lake (final lake < 50 % covered by the present set);
  - other.

## C — carve_diag
- **Laid by the nearest propagation**: a carved cell whose value is its `who` sample's cone.
- **Laid by the cross-line minimum**: another line's cone, from a neighbour's `who`, is lower.
- The cone is rebuilt in the bench, `v = zf + tan(28°)·max(d − hw, 0)` (`wall_profile` None, foot 0), with
  `geo`'s f32 arithmetic. **The rebuild must reproduce `out` bit for bit** on every carved cell, or the split
  is not reported.
- **Populations**:
  - the realigned témoin's constructed cells (OFF);
  - the residual of path dependence (extended, A/B patches; class in A and in B);
  - Finding 127's dams: the 4 true cols (3852, 2337), (4551, 3280), (2253, 4495), (1805, 4580), each 3 × 3,
    in B2 → A_c realigned (`a_min_km2 = 0.1`). The laying sample's line and its area are given.

## V
Two viz layers that read the world and never write it:
- "Ombrage": hillshade, azimuth 315°, altitude 45°, on the CONDITIONED field in metres, multiplied into
  the hypsometric colour;
- "Différence": current − a stored reference world, in metres, a signed blue–white–red scale with its
  saturation chosen and shown.

The field hash is checked unchanged (the guard). The "capture" is the layer's own buffer, from the same
`layer_color_image` the viz draws, cropped on lake 2, written by a test. It is NOT a screen grab, and is said
so.
