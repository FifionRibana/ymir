# Finding 123-A — the teeth are on the constructed walls, but not as the hypothesis states them

Draft for the ADR. Seed 1, 8192², the C2/10 world (Finding 121), the export layer
(`full_tree = false`). Bench `f123_teeth.rs`, output `a_teeth.txt`. Blocks B, C2 and E wait on the
author's block D (the round orders D before B).

## Grep (rule 11/11b)

`Genevaux hierarchical` **0** · `hillslope profile` **0** · `mur` **0** · `river id` / `river_id`
**0** · `longueur de tronc` **0** · `wall` 34 (earliest L187, the V-walls of the regime split) ·
`convex` 3 (L654: *"talus for straight repose flanks … or nonlinear critical-slope for convex flanks —
but that is the GS solver again"*) · `stream_threshold` 6 (L877: segments are selected by an
accumulation threshold, not by traversal) · `trunk length` 1 (L13935, Finding 93: *"trunk length p50
1.37 km · p90 14.3 · max 37.1"*) · `Finding 120` **0** (not yet in the ADR).

## Instrument

- **Walls**: `carved && !floor` of the PRE carve at `F121_AGE_K` (Findings 121/122's approximation).
- **A tooth**: a drawn `Watercourse` segment with ≥ 80 % of its cells on a wall; "on walls" is
  ≥ 50 %.
- **A**: the D8 drained area in-domain at the segment's OWN cells. ⚠️ The first run read the last
  point, which is the junction ON the parent, and so read the trunk's area (teeth p50 29.4 km²,
  and teeth counted among the "trunks" of the control). Corrected before any reading below.
- **Angle to the parent**: the acute angle between the tooth's last 8 cells and the parent's cells
  around the junction.

## The teeth

| | value |
| --- | --- |
| drawn segments | 12 609 (16 981 km in-domain) |
| on walls (≥ 50 %) | 2 902 |
| **teeth (≥ 80 %)** | **1 787 = 61.6 % of the on-wall segments** |
| teeth length | 2 094 km = **12.3 %** of the drawn network |
| teeth A p10 / p50 / p90 | **0.122 / 0.269 / 2.277 km²** (A_c = 0.1) |
| teeth length p50 | 0.93 km |
| angle to the parent p10 / p50 / p90 | **0 / 56 / 86°** (off-wall segments: p50 23°) |
| teeth by Strahler | 1: 1 075 · 2: 673 · 3: 35 · 4: 4 |
| CONTROL — trunk segments A ≥ 10 km² (4 886) | wall share p50 **0.0 %**, p90 **0.0 %** |

**Verdict against the round's three clauses** (≥ 70 % of on-wall segments, ~90°, A median < A_c):
**none holds.** 61.6 %, p50 56°, A p50 0.27 km². The teeth ARE on the walls — 1 787 segments at
≥ 80 % wall — so the hypothesis does not fall on location, but they are not the sub-A_c
perpendicular D8 lines it describes. They cannot be sub-A_c: a drawn segment exists only above the
extraction threshold, and Finding 110 measured 0 of 11 848 below A_c.

## The "S2, 580 km², 654 km" entry — NOT reproduced

The microscope's `aggregate_watercourses` is ported line for line (lake chaining included; the
shared-terminus catchment branch omitted). No S2 watercourse of this world comes near: the closest
are 459 km² / 300 km and 798 km² / 299 km. The longest S2 trunk here is about 300 km signified. The
author's viz world differs from this bench's in a way not identified (resolution, framing,
parameters); the entry's mouth coordinates would settle it.

## C1 and C2 read, not applied

| drawn layer at | segments kept | teeth kept | S2 watercourses to the sea that vanish |
| --- | --- | --- | --- |
| A ≥ 0.1 km² | 12 609 | 1 787 | 0 of 55 |
| A ≥ 0.5 km² | 7 308 | 616 | 0 of 55 |
| A ≥ 1 km² | 6 299 | 370 | 0 of 55 |
| Strahler ≥ 2 | 6 646 | 712 | 0 of 55 |

Identity (the trunk as the max-A path, lake chaining kept): over 647 watercourses, 557 trunks
unchanged, **90 shorter**, none longer; length ratio p10 **0.21**. E.g. the 459 km² S2 entry:
shipped trunk 300 km, max-A trunk 22 km. The shipped rule climbs the LONGEST path in points, which
through lake chaining can follow a long minor branch; whether the 90 shortened trunks are Finding
93's lake-crossing trunks is not checked.
