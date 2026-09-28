# Finding 127 — the foot clause returns B1's teeth to the profile's own; the off-grid tracé fails its calibration because the D8 tree it must keep is itself the lattice; the canyons are dammed by a floor, not dug

**One clause (A) and one method (B) were built, both gated and off. Nothing is promoted.** The benches
are in `crates/ymir-core/tests/f126_coast.rs`, Finding 127's section; the raw outputs are beside this
report. **mur ↔ mer is ON in every world** (proved in Finding 126).

**Units.** 1 cell = 48.8 m (domain). "Signified" = ×7.5 for lengths, ×56.25 for areas. Heights are
metres. Areas are domain km².

**Reading declaration, unfavourable.** My predictions (`f127_predictions.md`) were written before any
grep or measurement. They are non-blind on Findings 73–126 and on the reviewer's predictions. The B0
and B1 gates were declared before any run.

## 0 — Finding 126 committed; rule 11c; the guard

- Finding 126 and **method rule 11c** were committed as `8a5d1fc`. The rule, in the ADR right after
  rule 11b: every grep on the dossier covers its typographic variants (`−` `–` `-`, `×` `x`, `≥` `>=`,
  `≤` `<=`, `’` `'`, `∞` `inf`, thin spaces). The table counts each variant in the ADR.
- The guard through `run_hd` then read **6 / 6 "= banc"**.
- This round's code (`foot_quiet`, `skeleton_trace`, and the polyline builder split into `line_samples`
  and `densify`) is **bit-neutral when off**: a fresh uncached build of the definition hashes to
  `a8d2d538d692c2f0`, the reference.

**The rule-11c greps.**
- "Tarboton": NOTHING FOUND (the expanded grep).
- "D∞": 5 hits. The earliest (Finding 10, L583) is about ROUTING.
- **A continuous-angle method was tried once, and for routing only.** `FlowConfig::dinf` was kept off,
  measured inferior (`probe_dinf_compare`, `flow.rs:45`): its trace re-quantised to the primary D8
  neighbour, and on flats its BFS distance gradient is cardinal (R2 0.40 → 0.91).
- Finding 112 (L17552) records "no MFD path and no D∞ tracer in the crate".
- The true col is Finding 119's (L17867): the receiver of `Lake::outlet`, which is itself a water
  cell.

## A — bruit ↔ pied: the foot clause removes the interaction and returns B1's teeth to the profile's own level

**The clause, gated: `WallProfile::foot_quiet = Some(τ)`, τ = 0.5 (PROXY).** The detail's amplitude is
0 on the concave foot (u < 300 m beyond the floor edge), full beyond (1 + τ)·300 = 450 m, with a
smoothstep between; the transition is declared as a fraction of the foot's own length. The profile
stays. `carve` now carries the `u` of the sample it keeps through the cross-line minimum. Permanent
test `the_foot_clause_silences_the_detail_on_the_foot_only`:
- the weight is 0 on the foot, 1 beyond, and monotone;
- negative control: the clause changes wall cells;
- it never raises a cell;
- without a foot it is a bit-exact no-op.

**The two identity controls, read first.**
- The témoin + mur ↔ mer hashes to `a8d2d538d692c2f0`, the definition's reference: the coast clause
  never binds on a planar wall.
- B1b (no foot) + mur ↔ mer hashes to `bbfad1625373130b` both without and with the foot clause:
  **bit-identical**, as the round required.

| world (mur ↔ mer ON) | teeth | mouths / heads / points at the FOOT | spurs near a coastal wall | σ p50 | R8 terrain | R8 network | over-dug | F38 |
|---|---|---|---|---|---|---|---|---|
| **témoin C2/10 col** | **1 440** | 28.3 / 8.5 / 21.4 % | 0.0142 /km | 5.573 m | 0.0452 | 0.5613 | 0 | holds |
| B1 | 4 805 | 89.7 / 59.6 / 49.3 % | 0.0095 | 5.489 m | 0.0492 | 0.5407 | 0 | holds |
| **B1 + FOOT CLAUSE** | **1 986 (−58.7 %)** | **48.4 / 23.0 / 27.8 %** | 0.0095 | 5.516 m (**+0.5 %**) | 0.0494 | 0.5476 | 0 | holds |
| B1a, profile only (Finding 126, no mur ↔ mer) | 2 130 | 46.4 / 22.2 / 27.8 % | — | — | 0.0497 | 0.5513 | 0 | holds |
| B2 (A_c) + B1 | 1 453 | 23.2 / 41.9 / 32.2 % | 0.0066 | 9.925 m | 0.1034 | 0.5054 | 2 | holds |
| **B2 (A_c) + B1 + FOOT CLAUSE** | **1 084 (−25.4 %)** | 18.9 / 45.8 / 27.8 % | 0.0066 | 10.252 m (**+3.3 %**) | 0.1036 | 0.5230 | **3** | holds |

- **On B1 the clause removes the interaction exactly.** The teeth fall to 1 986, at B1a's level (the
  profile alone, 2 130 in Finding 126 without mur ↔ mer). Their positions become B1a's to the percent:
  27.8 % of points at the foot against 27.8 %, mouths 48 % against 46 %. What remains over the témoin
  (+38 %) is the profile's own share. This clause does not address it, and it was not asked to.
- **It touches nothing else on B1**: spurs at the base rate, R8 terrain +0.0002, no canyon, F38 holds.
  **The texture's price is small: σ +0.5 %.** Removing the detail from the foot does not smooth the
  land; the foot is a small share of it.
- **On B2 (A_c) + B1 the teeth fall 25 %, and the clause REVEALS a canyon** (rule 18). A third
  over-dug body appears: body 3, 1.18 km² (domain), floor (3895, 2306). That is where B2 → A_c's canyon
  3 sits (floor (3892, 2310), Finding 125). My reading, NOT measured: C shows the dams are floors;
  near a dam the foot's noise on the neighbouring walls breached it, and without the noise the dam
  stands. **σ rises 3.3 % there.** The ~ 55 % of wall area that is foot
  at A_c loses its texture.
- **⇒ bruit ↔ pied is PROVED for what it was built for** (the interaction, seed 1, 8192²). Its cost is
  measured: +1 canyon on B2 → A_c's skeleton, and σ +3.3 % there.

## B — the off-grid tracé FAILS its own calibration (B0). B1–B3 were not run

**What was built, gated and off.** `ValleyConstruction::skeleton_trace: Option<SkeletonTrace>`, with
`SkeletonTrace::f127(interp)`: the 1–10 km² segments, a flat threshold of 1 m/km (**PROXY**). The
accumulation, the junctions and the choice of trunks stay D8. The MEDIAN LINE of each segment with
`a_min ≤ A < 10 km²` is retraced by continuous steepest descent on the construction's input field:
- the field is interpolated bicubically (Catmull-Rom with its analytic gradient), or bilinearly as the
  control;
- a step is half a cell, from the segment's source until the path enters its D8 receiver's corridor
  (± W/2, checked within 6 cells);
- on a flat (|∇| < 1 m/km), one step follows the D8 path;
- a trace that first enters ANOTHER line's corridor, or the sea, or runs out of steps (16 × its D8
  length) keeps its D8 geometry, and is counted;
- the floor and the width follow the D8 prefix by arclength fraction.

Permanent test `the_trace_leaves_the_lattice_and_keeps_its_receiver`. Its negative control: on a
rounded V valley at 22.5°, every one-cell chord of the D8 trunk sits 22.5° off the axis. The trace
stays within 8° and reaches its receiver. When off, the definition's world is bit-identical (fresh
build `a8d2d538d692c2f0`).

**Not Finding 113's snap** (written so a grep does not close the question): that moved an EXPORTED
river toward a thalweg it already occupied. This replaces the skeleton's geometry before the
construction digs, at a continuous angle.

**Not the `dinf` of `FlowConfig`** either. The rule-11c grep found that D∞ ROUTING was tried and kept
off (`probe_dinf_compare`). Its trace re-quantised to the primary D8 neighbour, and on flats its BFS
distance gradient is cardinal (R2 0.40 → 0.91). This tracé keeps continuous positions and never takes
a flat gradient: it falls back to D8 instead.

**B0 — the calibration, BLOCKING.** The gate was declared before the run: bicubic skeleton R8 ≤ 0.08 in
the 1–3 and 3–10 km² bands at all four roughnesses. The instruments are Finding 126's isotropic
synthetics and 8-point chords, with the band read at the nearest trunk cell. D8 is re-read with the
same rule.

| roughness ρ (ratio) | D8 raw, 3–10 / 1–3 | D8 + 250 m | **tracé bicubic** | tracé bilinear | bicubic: SAME receiver / other / lost |
|---|---|---|---|---|---|
| 0.1 (0.071) | 0.879 / 0.986 | 0.786 / 0.977 | **0.768 / 0.974** | 0.762 / 0.970 | **12.9 %** / 83.8 % / 3.4 % |
| 0.3 (0.209) | 0.749 / 0.797 | 0.553 / 0.663 | **0.511 / 0.574** | 0.500 / 0.543 | **36.9 %** / 55.1 % / 8.0 % |
| 1 (0.613) | 0.320 / 0.363 | 0.117 / 0.178 | **0.104 / 0.110** | 0.107 / 0.098 | **42.6 %** / 13.9 % / 43.5 % |
| 3 (1.071) | 0.186 / 0.240 | 0.061 / 0.126 | **0.062 / 0.111** | 0.058 / 0.114 | **23.6 %** / 5.3 % / 71.1 % |

**⇒ B0 FAILS** (worst 0.974 against 0.08). By the round's rule, **B1–B3 were not run**, and the real
skeleton was never retraced.

**Where it fails — the attribution, not a remedy.**
- **The method's constraint breaks, not its interpolant.** Only 13–43 % of the traces reach the
  receiver the D8 tree assigns them. The rest keep their D8 geometry, so the skeleton's R8 stays near
  D8's.
- **On smooth terrain the D8 TREE is itself the lattice.** At ρ 0.1, 84 % of the traces enter ANOTHER
  line's corridor first. A continuous descent follows the true (radial) slope, while D8's parallel runs
  within each 45° sector define receivers that no continuous path reaches. Keeping the D8 tree while
  freeing the geometry is inconsistent exactly where the lattice dominates, at a roughness ratio below
  ~0.4, which is where the real terrain sits (0.379, Finding 126-D1).
- **On rough terrain the descent is trapped.** At ρ 3, 71 % are lost: they run out of steps in the
  closed depressions of the UNBREACHED input field, which the D8 tree crosses by its breach.
- **Does the descent print its own lattice?** Mostly not, but noisily. The traced PATHS themselves, all
  outcomes, read R8 0.039–0.179 at 1–3 km² and 0.029–0.716 at 3–10 km². The counts are small where the
  values are high (733 chords at 3–10 km², ρ 0.1).
- **The bilinear control does not separate from the bicubic** at any roughness (skeleton R8 within 0.03;
  paths within 0.07). The interpolant is not what fails.
- **Named, not built:** a tracé whose TREE is continuous too, i.e. accumulation along continuous paths.
  That is a routing change, and Finding 112 warns it trades away the concentration. Or tracing on the
  BREACHED field (the D8 tree's own), which removes the pits but puts back the breach's D8-staircase
  channels.

## C — the canyons at their true col: they are DAMMED, not dug, and the dam is a FLOOR laid higher

**Instrument.** The true col is Finding 119's: the receiver, on the HD assembly's own flow, of
`Lake::outlet`, the lake's last WATER cell. It is read in the canyons' world (B2 → A_c) and at the same
cells in the témoin and in B2 ≥ 1 km², where no canyon forms. "Owner" is the construction's masks on
PRE (a FLOOR, a WALL, or UNCARVED = the terrain, a rim) with the band of the nearest trunk. These are
planar worlds, so mur ↔ mer never binds and they are Findings 125–126's worlds bit for bit. Heights
are metres.

| canyon (km², domain) | PRE floor / col | **floor**: témoin · B2 ≥ 1 · A_c | **true col**: témoin · B2 ≥ 1 · A_c | col's owner, témoin → A_c | S at the col (A_c) |
|---|---|---|---|---|---|
| 3 (1.49) ⚠ | 2 179 / 1 731 | 272.6 · 296.3 · 297.7 | 294.5 · 295.0 · **405.2** | trunk WALL → **FLOOR of a < 1 km² valley** | 0.000 |
| 7 (10.55) | 2 322 / 1 819 | 646.7 · 646.7 · 647.0 | 623.6 · 623.6 · **787.6** | trunk FLOOR → **FLOOR of a < 1 km² valley** | 0.021 |
| 13 (15.37) | 1 479 / 1 620 | 328.4 · 328.4 · 329.0 | 970.2 · 771.8 · 537.3 | trunk WALL → WALL of a < 1 km² valley (a SIDE spill, see below) | 0.000 |
| 14 (41.97) | 1 223 / 1 018 | 134.2 · 134.2 · 135.0 | 119.6 · 119.6 · **493.7** | trunk FLOOR → **FLOOR** (trunk band) | 0.037 |

⚠ Body 3 fails the instrument check: its recorded outlet is not a cell of its lake on the assembly's
map. Its row is read with that caveat. Bodies 7, 13 and 14 pass: floor and outlet in the lake, col
outside.

- **No canyon is a pre-existing depression.** PRE's own fill is 0.0 m at all four floors.
- **No canyon's floor is dug.** The floors sit at the same height in all three worlds: within 1 m for
  7, 13 and 14, +25 m for 3.
- **The COL rises.** In bodies 7 and 14 the col cell lies BELOW the floor in the témoin (623.6 < 646.7;
  119.6 < 134.2): the valley drains. Under B2 → A_c the same cell sits +164 m and +374 m higher, and
  the valley is closed. Body 3's col rises +110 m.
- **What holds the dam is a FLOOR, not a wall or a rim.** In 7 and 3 the col becomes the floor of a
  < 1 km² valley; in 14 it stays a floor in the trunk's band, laid 374 m higher.
- **The construction does not RAISE anything.** Every col is still far below PRE (493.7 against
  1 017.9 m). The dam is RELATIVE: at A_c the col cell is lowered LESS than at 1 km².
- **The mechanism, by the code** (read, not measured cell by cell). `carve` lays each cell from its
  NEAREST skeleton sample ("nearest, not lowest", by design, so a downstream sample cannot flatten an
  upstream floor). At A_c a 0.1 km² tributary's own floor rises fast upstream: χ grows by
  (A₀/A)^0.5 · dx ≈ 3.16 × 48.8 m per cell, so zf gains ≈ 11 m per cell. A trunk-floor cell whose
  nearest sample now belongs to such a tributary takes that tributary's higher floor. The sample
  positions are smoothed (±250 m) and densified, so a tributary's samples can lie nearer to a trunk cell
  than the trunk's own. The continuity minimum across lines only looks at the 8 neighbours' samples.
  The instrument does not expose which sample carved the col (`who`), so "whose floor" is the masks'
  band, not the sample's line.
- **Body 13 is a side spill.** Its recorded col (537.3 m) is a < 1 km² valley's WALL, LOWER than in the
  témoin (970.2 m). The lake spills sideways because its trunk route (the témoin's drainage) is dammed
  elsewhere, above 537.3 m. That dam is not located here. It is the next read: walk the témoin's D8
  path downstream of the floor and compare the heights.
- **The clause this points to, named and not built — the confluence clause:** "a tributary's floor may
  not stand above its trunk's floor across the trunk's corridor". Build the trunk first, then let a
  tributary sample take a trunk cell only where it is lower.

## The clause table, updated

| unwanted effect | Finding 126's status | what Finding 127 measured | **Finding 127's status** |
|---|---|---|---|
| **B1 spurs** | mur ↔ mer: PROVED | ON in every world; the base rate holds in every one (0.0066–0.0142 /km) | **PROVED, kept ON** |
| **Finding 38** | the same clause: PROVED | holds in every world with mur ↔ mer ON | **PROVED** |
| **B1 teeth** | located in the foot; *bruit ↔ pied* named | the foot clause: B1 4 805 → 1 986 (−58.7 %), back to the profile's own level and its positions to the percent; B2 (A_c) + B1 −25 %; B1b bit-identical | **PROVED for the interaction.** Cost measured: σ +0.5 % (B1) / +3.3 % (B2 + B1), and +1 canyon revealed on B2 → A_c. The profile's own +38 % over the témoin is left, and not addressed |
| **B2 terrain R8** × 2.3 | the moving average REFUTED; the alignment is the D8 lattice's | the off-grid tracé on the D8 tree fails its own calibration (B0): 13–43 % of traces reach their D8 receiver; divergence on smooth terrain, pits on rough | **"a continuous geometry on the D8 tree": REFUTED at calibration** (B1–B3 not run). What remains is a CONTINUOUS TREE (routing, which Finding 112 warns trades the concentration): named, not built |
| **B2 → A_c canyons** | "intersection" REFUTED; look at the sill | the dam is a FLOOR: the col rises +110 / +164 / +374 m at unchanged floors; body 13 is a side spill | **LOCATED at the sill.** The candidate is named, not built: the *confluence clause* (a tributary's floor may not stand above its trunk's floor across the trunk's corridor) |

## Predictions scored

**Mine** (written before the round):
- **P0 HELD.**
- **A.** **P-A1 REFUTED, narrowly**: −58.7 % against my −45 to −56 %. **P-A2 HELD.** **P-A3 HALF**: σ
  +0.5 % on B1, +3.3 % on B2 + B1. **P-A4 HALF**: canyons +1 on B2 + B1.
- **B.** **P-B0 REFUTED** (FAIL; bilinear does not separate). **P-B1 / B2 / B3 NOT MEASURED** (stopped
  by B0).
- **C.** **P-C REFUTED**: the floor is not dug, the col rises.
- **Meta HELD.**

**Reviewer's** (in the prompt):
- **A**: −60 to −75 % **REFUTED, narrowly** (−58.7 %); σ −3 % **REFUTED** in sign (+0.5 / +3.3 %);
  spurs at the base rate and R8 unchanged **HELD**.
- **B0**: bicubic 0.05–0.10 and bilinear 0.15–0.30 **REFUTED** (bicubic 0.11–0.97; the two agree).
- **B1–B3 NOT MEASURED.**
- **C**: a barrage in 3 of 4 **HELD in count**, but it is held by a FLOOR, not a wall or rim; "the 4th a
  pre-existing depression deepened" **REFUTED** (PRE's fill is 0 at all four; the 4th, body 13, is a
  side spill).
- **Meta HELD.**

## Limitations, stated

1. **B0's gate reads the skeleton, divergent traces included, which keep their D8 geometry.** The paths'
   own R8 is a diagnostic added after the first partial run. It changes no verdict. It is noisy where
   chords are few (94–733 chords in some rows).
2. **B0's "other receiver" counts any non-upstream corridor first entered**, a sibling tributary's near
   the junction included. At ρ = 3 "other" is only 5 %, so this cannot carry the failure. "Lost" is a
   step budget (16 × the D8 length); a larger budget was not tried.
3. **B1–B3 were not run**, by the round's rule, so the real terrain's topology share is NOT measured.
   The synthetics bracket it (the real roughness ratio is 0.379, between ρ 0.3 and 1: 37–43 % reach
   their receiver there).
4. **C's owner is the construction's masks on PRE (floor / wall / uncarved) with the nearest-trunk
   band**, not the sample that carved the cell (`who` is not exposed). Body 3 fails the instrument check
   (its recorded outlet is not in its lake), so it is read with that caveat. Body 13's dam is not
   located.
5. **A's canyon revealed on B2 + B1 is located, not attributed.** "The noise was breaching the dam" is a
   reading.
6. **One seed, one resolution (8192²).**

## Uncommitted state (awaiting the author)

- `valley_construction.rs`:
  - `WallProfile::foot_quiet` (the foot clause) and `foot_weight`;
  - `ValleyConstruction::skeleton_trace`, with `SkeletonTrace`, `TraceInterp`, `TraceStats` (and its
    diagnostic `paths`);
  - `interp_grad` (bicubic / bilinear with gradient) and `traced_polylines`;
  - the polyline builder split into `line_samples` + `densify`, byte-identically;
  - two permanent tests with negative controls.
- `cached_product.rs`: the key test lists `skeleton_trace` and the nested `foot_quiet`.
- `tests/f126_coast.rs`: Finding 127's benches (`f127_a`, `f127_b0`, `f127_b12` — not run, `f127_c`).
- This report.

**Nothing promoted.**
