# Finding 175 — option 5: an ocean-ocean subduction closure in C1, with the seaward margin

**Status: built behind a flag (no solver, two closures changed as declared, A1–A7), measured by tiers; the success fails on
0 / 5 seeds (the active share); one point is stopped because it needs a closure change beyond A1–A7; nothing committed without
the feu vert.** Production unchanged: C1's production output is bit for bit (a hash of three seeds, HEAD against the F175 tree);
`ALGO_TECTONICS` is **not** bumped — see « The guard and a stale cache » below.

- Declared before any measurement: `f175_declared.md` (with the amendments after the smoke run). Predictions:
  `f175_predictions.md`.
- The bench: `crates/ymir-core/tests/f175_cascade.rs`, four resumable tiers (`f175_tier1_c1` … `f175_tier4_cascade`); each
  world's lines in `results/`, assembled in `f175_bench_output.txt`. Images in `images/` (north up, checked; the Copernicus notice
  in `images/NOTICE.md`).

## The bench by tiers (the author's budget rule)

The first bench (100 worlds, all at 2 048², 4 ages, 3–5 h) was disproportionate — the reviewer's error, recorded — and was
interrupted on seed 8; its output was discarded. The revised bench, announced at ~25 min wall / ~1 h CPU:

| tier | content | wall | CPU (est.) |
|---|---|---|---|
| 1 | C1 only: 5 seeds × Part A / B′ / B′ + option 5, the fixed continent (F174 B′'s plates); the controls | 21 s | ~1 min |
| 2 | the selection redone with option 5, once per seed (parallel) | 127 s | ~6 min |
| 3 | the physics level at 256²: 15 worlds × young and old | 472 s | ~11 min |
| 4 | the cascade to 2 048²: seeds 8 and 18 (the most Andean cells after 8), young, B′ + option 5 | 140 s | ~16 min |

Total **~12.6 min wall** (tier 4 was run twice: the shelf was added after its first run, +2.3 min).

## What was built (A1–A7)

- **The flag**: `OoSubductionParams` (D_init = 100 km ÷ 7.5 = 2.13 cells of 64²; Stern & Gerya 2018 p. 44 after Hall et al.
  2003), carried by `Margin::oo` and passed to C1 through `run_with_closures_terms_oo` / `run_c1_sources_mapped_oo`;
  `C1Closures` untouched (its `Debug` is in the tectonic cache key).
- **A1** the older boundary cell sinks (ties: the velocity flag, then the lower plate id), Δs = the existing rate;
  **A2** the arc on the overriding plate's continental cells, else its oceanic cells; **A3** the floor-triggered reassignment
  to the overriding plate (the existing mechanism; the type stays oceanic); **A4** the per-pair D_init tracker, no rate
  threshold; **A5** the existing C-O rule takes over at the continent; **A6** no spontaneous collapse; **A7** Davis-Suppe's
  seeds cleared on O-O convergent cells (`clear_oo_seeds`).
- **The controls** (permanent, `the_oo_subduction_path_controls`; and per seed at tier 1): the flag off is C1 bit for bit
  (5 / 5); D_init infinite consumes nothing O-O and leaves plate_id unchanged (5 / 5), and differs from B′ by A7 alone (|Δs| max
  1.4–2.1, mean 0.02–0.14); the arc mass is conserved (error ≤ 2 × 10⁻¹³); a second A7 cut finds nothing.
- **The viz**: « Subduction océan-océan (C1) », and « Marges à la fin (actives 3, passives 2, arcs ½) » in the terms' frieze.

## Tier 1 — C1 (the fixed continent)

| seed | active at the end: A / B′ / B′ + 5 | Andean margins (B′ + 5) | O-O cells consumed | reassigned | pairs subducting | arc mass on continent / ocean | DS on the band: B′ → B′ + 5 | r* (B′ + 5) |
|---|---|---|---|---|---|---|---|---|
| 8 | 30.9 / 0.0 / **4.4 %** | 1 (4 cells) | 20 711 | 491 | 30 | 1.57 / 15.07 | 1 837 → 564 | 4 |
| 18 | 3.7 / 1.9 / **2.3 %** | 1 (3 cells) | 6 649 | 121 | 22 | 0.78 / 6.72 | 520 → 234 | 2 |
| 19 | 29.8 / 0.0 / **0.0 %** | 0 | 10 293 | 12 | 27 | 0.26 / 13.51 | 1 055 → 355 | 4 |
| 23 | 0.0 / 0.0 / **0.0 %** | 0 | 10 844 | 104 | 31 | 2.30 / 14.92 | 1 187 → 386 | 1 |
| 29 | 17.9 / 0.0 / **0.0 %** | 0 | 18 925 | 30 | 29 | 0.31 / 25.60 | 1 009 → 446 | 3 |

- **The O-O path works and consumes a lot, but almost never reaches a coast**: 0–4.4 % of the perimeter active at the end; 0–1
  Andean margin per seed.
- **87–99 % of the arc mass lands on oceanic cells** (island arcs, in orange on `f175_margins_end.png`: lines along the
  O-O boundaries in mid-ocean) — which have **no altitude** under Airy (declared in A2). No island at 64² or 256².
- **Davis-Suppe on the band falls by 55–69 %** (A7): the O-O fallback dome is gone.
- **The validation criteria**:
  - **V1**: 100 % of the initiations (22–31 per world) lie within 2 cells of a continent — the O-O boundaries the continents
    meet are the bands' outer edges; far above the note's ≤ 15 %.
  - **V2**: 96–100 % within 200 km of the previous step's consumption.
  - **V3**: 22–31 initiations per world, 59–78 % of them in the first 25 steps, the rest at the merges (steps 53–64) and a few
    late (74–157).
  - **V4**: by construction the older side sinks on every age-decided consumption (0 violations); 0–54 exact ties (velocity
    rule); but **10–36 % of the decisions rest on |Δage| < 0.01** — the advected age's smearing, not a clock.
  - **V5**: the arc mass conserved; one-cell fingers at the end against B′: 10 / 9, 0 / 0, 1 / 1, 2 / 2, 10 / 8 (above B′ on
    seeds 8 and 29).
  - **V6**: see the active share and the Andean margins above.

**Why the trench does not reach the continent** (a reading of the code and of the numbers, a hypothesis to check): the
consumption is at most 0.5 × v_rel · n̂ × dt ≤ 0.5 × dx per step (the CFL), so a distal band cell (s 0.20–0.29) needs ≥ 19 steps
— ~40 at the typical convergence — to fall under the floor (0.05); the distal band is ~3 cells deep. But the **accretion
closure merges a plate pair after 50 consecutive convergent steps** (`AccretionParams::merge_time_threshold`): the pair merges
before the band is consumed, the boundary disappears, and the activity decays with the plate count (16 → 2–4 by step ~100–150).
Making the margins turn Andean would need the O-O consumption to outpace the merge (or the merge to wait for it) — **a change
of a closure beyond A1–A7: stopped, not implemented.**

## Tier 2 — the selection redone with option 5

Found on **5 / 5** seeds (N 16 on 4, 20 on seed 23; land 52.5–58.6 %, r* 2–3); a different continent from B′'s on 4 seeds.
Nothing else was measured on the reselected continents.

## Tier 3 — the physics level at 256² (the fixed continent)

| seed | mountain at young: A / B′ / **B′ + 5** | W₁₀ young: A / B′ / B′ + 5 | W₁₀ active / passive (B′ + 5, young) | islands | lakes (B′ + 5, young) |
|---|---|---|---|---|---|
| 8 | 71.0 / 89.7 / **71.5 %** | 0.666 / 0.695 / 0.139 | 0.621 / 0.120 (190 / 4 614 cells) | 0 | 11 (6 636 km²) |
| 18 | 47.1 / 29.6 / **22.6 %** | 0.397 / 0.051 / 0.044 | — / 0.044 (86 cells active, none steep) | 0 | 3 (156 km²) |
| 19 | 82.3 / 58.4 / **40.3 %** | 0.649 / 0.202 / 0.021 | — / 0.021 | 0 | 8 (3 135 km²) |
| 23 | 59.4 / 62.6 / **48.8 %** | 0.420 / 0.111 / 0.039 | — / 0.039 | 0 | 1 (435 km²) |
| 29 | 69.8 / 48.6 / **43.2 %** | 0.572 / 0.097 / 0.105 | — / 0.105 | 0 | 5 (166 km²) |

- **Seed 8's plateau is undone in part**: 89.7 → 71.5 % mountain at young (−18.2 points; −26.8 at old, 78.4 → 51.6 %); the median
  1 778 → 971 m. No world above 80 % mountain with option 5.
- **W₁₀ falls further** (option 5 removes the band's Davis-Suppe wall): 0.021–0.139 at young.
- **Active against passive**: only seed 8 has active coastal land at 256² (190 cells): W₁₀ 0.621 against 0.120 (×5.2).
- **No island** at 256² on any world (the arcs are oceanic, invisible).

## Tier 4 — the cascade to 2 048² (seeds 8 and 18, young, B′ + option 5)

- The coast: dimension 0.979 / 0.985, length ratio 0.897 / 0.892 — as smooth as before; low (cover) windows (20 / 28), their
  dimension 1.002 / 0.988 (the sandy reference 0.95–1.00); the seaside slope p50 2.1° / 0.7°.
- **The shelf** (F174's break rule): median 0 km on both (p90 1.1 / 1.4 km; break 117 / 90 m). **S7 held.**
- Islands ≥ 1 km²: 0 / 2 (5 km², not on arc deposits).
- The transects (`f175_transect_seed8.png`, `…seed18.png`): the Andean margin (red) carries a higher coastal relief and a steeper
  drop than the passive one (green); both fall to −5 km in ~15 km.

## The success and the stops

- **The stops**: the selection fails on 0 seeds; the circumnavigation is lost on 0 → neither fires.
- **The success** (≥ 4 / 5 at the fixed continent): DS O-O seeds 0 ✓ (5 / 5); no world above 80 % mountain ✓ (5 / 5); **10–60 % of
  the perimeter active ✗ (0 / 5)**; V4 ✓; V5 ✓ on 3 (fingers above B′ on 8 and 29); circumnavigable ✓ (5 / 5). **0 / 5: the
  seaward margin and option 5 do not enter the profile.**
- **The point stopped** (a closure change beyond A1–A7): the margins turn Andean only if the O-O consumption reaches the continent
  before the accretion merge — the merge threshold (50 steps) against the consumption rate (≥ 19 steps per cell, the band ~3 cells
  deep). And island arcs would need the arc-loaded oceanic cells to emerge (a promotion rule or a mapping change, declared in A2).

## Predictions

- **The reviewer's**: S1 refuted (−18.2 points, not ≥ 30); S2 refuted (0 / 5 in 10–60 %); S3 refuted on the islands, held on the
  circumnavigation (5 / 5); S4 refuted (only seed 8 has active segments, ×5.2 there); S5 refuted (2 seeds with one Andean
  margin); S6 held (5 / 5); S7 held (no shelf).
- **Mine**: held: S3's islands absent (the Airy reason), S6, S7, V3's count and timing, V4, the arc conservation; refuted: S1's
  −25 to −40 (−18), the ties ≥ 70 % (0–54 exact ties; the age decides, often on < 0.01), S2's 5–25 % (0–4.4 %), S4 and S5 on 2–3
  seeds (1 and 2), V1's 20–60 % (100 %), V2's « none » early (most initiations follow a consumption within 200 km), V5's fingers
  ≤ B′ (2 seeds above), the success 1–3 (0).

## The guard and a stale cache (found by this round, not caused by it)

- A first after-check ran with `ALGO_TECTONICS` bumped 1 → 2, to make the guard exercise the new time loop. Every cached state
  was recomputed (C1, then the HD upscale and erosion, chained to the tectonic key; 3 189 s) and the guard found **all six fields
  and lakes different** (`checks_after_with_algo_bump.txt`).
- **F175 is not the cause**: C1's production run hashes identically at HEAD (`b080a3b`, a separate worktree) and in the F175 tree
  (`tests/f175_c1_hash.rs`; seeds 42, 8 and the old témoin: c69acd85a0f7f484, 3565285d65e4c341, fd23c309498a9761 on both).
- **So the guard's references match cached products that the current code no longer reproduces**: a stale cache, older than
  F175, somewhere in the tectonic or HD steps. The bump was reverted (the hash check proves the bit-for-bit directly) and the
  after-checks rerun on the existing cache (`checks_after.txt`).
- **For the author**: which earlier change left the cache stale is not investigated here; re-baselining the guard means bumping
  the ALGO_* of the stale step and regenerating `bench_field_hashes.json` (`f123_guard`), which changes the guard's reference.

## Open, for the author (nothing chosen)

- (a) The accretion merge against the O-O consumption: a merge that waits for an active subduction, or a faster O-O
  consumption — a closure change beyond A1–A7.
- (b) Island arcs: promoting an arc-loaded oceanic cell to continental, or reading its s under Airy.
- (c) The age: C1's age is advected, not aged; 10–36 % of the O-O polarity decisions rest on < 0.01 of difference.
- (d) Whether to keep option 5 behind its flag: it removes the band's Davis-Suppe dome (seed 8's plateau −18 points, W₁₀ ×0.2–1) at
  no cost (C1 +3 s per world), but it does not make the coasts active.
- (e) The stale cache revealed by the guard (above): find the step, bump it, re-baseline the guard.

## At commit (2026-10-17)
- **Option 5 is committed behind its flag, off by default** (`Margin::oo = None`; `Margin::seaward_oo()` only on request).
- It **removes Davis-Suppe's dome at ocean-ocean boundaries** (seed 8: 89.7 → 71.5 % mountain at young; no world above 80 %),
  for +3 s of C1 per world.
- It **does not make the coasts active** (0–4.4 % of the perimeter): the accretion closure merges the plates before the trench
  reaches the continent, so active margins would need a change to a second closure.
- **The reviewer's error, recorded**: the success criterion « 10–60 % of the perimeter active » assumed option 5 alone would
  suffice. The reviewer proposes to close the margins workstream and accept passive margins — **the author's decision is
  awaited**.
- **The budget rule applied**: the bench ran by tiers (~13 min wall, resumable).
- **Outside F175**: the guard's references match a cache the current code no longer reproduces. C1's hash is identical (HEAD
  against the F175 tree); the gap is downstream. To be diagnosed in F176's Part 0.
- **Queued**: (a) the accretion merge against the ocean-ocean consumption; (b) emerged arcs (reading the s of arc-loaded cells
  in Airy); (c) the polarity decided by less than 0.01 of age difference.
