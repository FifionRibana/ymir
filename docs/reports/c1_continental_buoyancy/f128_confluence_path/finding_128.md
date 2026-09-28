# Finding 128 — the confluence clause removes every canyon but fails its concavity control; the path-based directions are read, and our D8-LTD fails its first reproduction on a tie, then reproduces the source once the tie is honoured

**One clause (A) and one method (C / D) were built, both gated and off. Nothing is promoted.**
- The benches are in `crates/ymir-core/tests/f126_coast.rs`, Finding 128's section.
- The reading is `reading_path_based.md`.
- The raw outputs are beside this report.
- **mur ↔ mer and bruit ↔ pied are ON in every world.**

**Units.** 1 cell = 48.8 m (domain). Heights are metres. Areas are domain km².

**Reading declaration, unfavourable.** My predictions (`f128_predictions.md`) were written before any
grep, reading or measurement. They are non-blind on Findings 73–127, on the reviewer's predictions and
on the addendum. I declared my memory of the papers as a belief to be tested. The C0 targets were
declared in the reading (§7) before the run. **No gate threshold for C was declared, because C did not
run** (rule 15).

## 0 — Finding 127 committed; the guard

- Finding 127 was committed as `0ae6243`. The guard through `run_hd` reads **6 / 6 "= banc"**.
- This round's code (`trunk_band`, `ltd_directions`, `accumulate_cells`, the exposed
  `Skeleton::direction`) is **bit-neutral when off**: a fresh build of the definition hashes to
  `a8d2d538d692c2f0` (`identity.txt`).
- `docs/refs/{orlandini2003,orlandini2014,paik2008}.pdf` stay uncommitted (rights).

**The rule-11c greps.**
- Orlandini, LTD / LAD, transversal deviation, path-based, nondispersive, Garbrecht, "nearest sample",
  "min des cônes": **NOTHING FOUND** in the ADR (every variant).
- D∞ appears 7 times, all about ROUTING.
- `FlowConfig` appears once (Finding 127).
- **The `dinf` history is not in the ADR but in commit `4f47495`** (2026-06-29), read whole in the
  reading §6. Its "Tarboton" is in the commit message, not in the dossier.
- The flat resolution the dossier already has is Garbrecht–Martz (`resolve_flats`, `flow.rs:220`).

## B — the reading (`reading_path_based.md`, with pages)

- **The parameter exists.** It is λ, a dampening factor, in [O03] eqs. (1)–(2), p. 1-3, not in the
  abstract. λ = 1 is advocated ([O03] §4). [O14] is parameter-free (p. 528). So λ = 1 is the method's
  definition, not a PROXY.
- **LTD over LAD, for the reason the sources give**: the transverse deviation is the one with an
  analytical basis ([O14] Theorem 3.1, Lemma 3.3 |δ⁺| ≤ h; Remark 3.8 excludes LAD). The sources do not
  say "LAD too sensitive on flats".
- **Paik 2008's two critiques.**
  - (a) λ has no analytical basis: answered by [O14].
  - (b) Asymmetry and false confluences on divergent terrain (a cone): NOT answered. Finding 126's
    synthetic is a dome, which is exactly that case.
  - Paik's numbers: the plane at 4 % of D8's lateral deviation, the cone at 60 %.
- **`FlowConfig::dinf` failed for two reasons.** It re-quantised to the primary D8 neighbour, which
  LTD cannot do. And the Garbrecht–Martz flat gradient is cardinal, which LTD would inherit. So LTD
  corrects the first, not the second. On flats the fallback is `compute_flow`'s pointer, counted.
- **B does not stop the round**: the methods do what the round needs, with two named caveats (flats,
  divergent terrain).

## A — the confluence clause removes every canyon, and FAILS Finding 121's concavity control

**The clause, gated: `ValleyConstruction::trunk_band`.**
- On the FLOOR BAND of any primitive (within W/2 of one of its samples), the LINE with the largest
  drained area covering the cell lays it, from that line's nearest covering sample. W is monotone in A,
  so the largest covering half-width ranks the lines.
- It is **not a minimum**, and it never acts within a line: a line's own downstream sample never takes
  over its upstream band.
- Permanent test `the_confluence_clause_lets_the_trunk_lay_its_band`, on a hand-made trunk and a
  higher-floored tributary. Negative control first: without the clause the tributary dams some trunk-band
  cells. With it, every band cell is at the trunk's floor; cells away from the junction are
  bit-identical; nothing is raised above the field.
- When off, the definition's world is bit-identical (`identity.txt`).

mur ↔ mer and bruit ↔ pied are ON in every world. The témoin reproduces Finding 127 to the digit.

| world | over-dug (canyons) | col 7 (4551,3280) | col 14 (1805,4580) | teeth | spurs near a coastal wall | σ p50 | R8 terrain | R8 network | F38 |
|---|---|---|---|---|---|---|---|---|---|
| **témoin C2/10 col** | **0** | **623.6 m** | **119.6 m** | 1 440 | 0.0142 /km | 5.573 m | 0.0452 | 0.5613 | holds |
| B2 → A_c, no clause | 4 (bodies 3, 7, 13, 14) | 787.6 m | 493.7 m | 73 | 0.0067 | 12.245 m | 0.1253 | 0.5356 | holds |
| **B2 → A_c + CLAUSE** | **0** | 714.2 m | 231.6 m | 22 | 0.0067 | 12.296 m | 0.1266 | 0.5263 | holds |
| B2 ≥ 1 km² + CLAUSE | 0 (0 without, Finding 126) | 623.6 m | 119.6 m | 275 (262 without) | 0.0056 | 9.409 m | 0.1063 (0.1049) | 0.5091 | holds |
| B2 (A_c) + B1 + foot + CLAUSE | **0** (3 without, Finding 127) | 704.7 m | 250.7 m | 538 (1 084 without) | 0.0066 | 10.240 m (10.252) | 0.1098 (0.1036) | 0.5131 | holds |

- **Every canyon goes.** 4 → 0 on B2 → A_c, 3 → 0 on B2 (A_c) + B1 + foot, and the lake count falls
  back to the témoin's (25 bodies, from 33). The coast, σ and both R8s barely move.
- **The cols of 7 and 14 do NOT come back to the témoin's height.** They fall (787.6 → 714.2 m,
  493.7 → 231.6 m) but stay +90.6 m and +112 m above it. Yet no depression closes: those cells are no
  longer the sills, and the water leaves another way. The col is the wrong read once the dam is
  removed; the over-dug census is the right one.
- **Canyon 13 was a confluence dam too.** Without the clause, the témoin's drainage path from its
  floor is dammed 52 cells down, at (2315, 4614), at 546.2 m against the témoin's 299.7 m (+246 m).
  With the clause, the floor carries no lake and the path is open.
- **Canyon 3's instrument, repaired and declared.** Its recorded `Lake::outlet` (3851, 2336) is not a
  cell of its lake on the assembly's map, which is why Finding 127 flagged it. The col is now taken as
  the receiver of the lake's LOWEST exit on the assembly's own flow (13 exits): (3907, 2310) at
  405.2 m, against the témoin's 394.6 m (+10.6 m). Under the clause its floor carries no lake.

**Finding 121's trap, the control the round required.** It is read on the témoin's 70 897 trunk
≥ 10 km² cells (skeleton on PRE), slopes to each cell's D8 receiver, with the concavity θ fitted from
S ∝ A^−θ over the cells whose slope is > 10⁻⁴ in both worlds.

| | \|Δz\| vs the témoin, p50 / p90 | slope ratio world / témoin, p50 | concavity θ, témoin → world |
|---|---|---|---|
| B2 → A_c, no clause | 1.47 / 43.42 m | 1.008 | 0.346 → 0.366 (+5.8 %) |
| **B2 → A_c + CLAUSE** | 1.06 / 32.40 m | **1.000** | 0.350 → **0.402 (+14.9 %)** |

- **The slope holds; the concavity does not.** The trunk floors come CLOSER to the témoin (p90
  43 → 32 m) and their median slope is unchanged. But θ rises 15 %, well outside the ±2 % the round set.
  My reading, NOT measured: a ≥ 10 km² tributary of a larger trunk now has its lower end laid at the
  larger trunk's floor across that trunk's band, which is a step at the band's edge and flat inside.
  That is Finding 121's flattening at the scale of W/2 of the RECEIVING line.
- **The source cones are touched.** The clause moves 388 470 cells at the construction (on PRE, the
  A_c skeleton):
  - 244 930 are raised: where a smaller line's cone had cut below the larger line's floor, or had carved
    a cell the larger line's floor leaves above the terrain (so "all on a floor band" is false);
  - 143 540 are lowered: the dams removed;
  - **90 097 lie within 3 cells of a line's head.**
- **⇒ The confluence clause is PROVED for the canyons, and FAILS the round's own trap control**
  (concavity +15 %, source cones moved). It is kept gated and off, and not advanced as a remedy until a
  form passes both.
- **Named, not built.** Restrict the takeover to the larger line's band DOWNSTREAM of the junction point
  (where the dams sit), leaving the tributary's own profile upstream of it. That keeps "no minimum"
  while not re-laying the tributary's approach.

## C0 — the implementation against the source: the FIRST run FAILS (the code); C and D stop this round

**What was built, gated and off.**
- `ValleyConstruction::ltd_directions` makes the skeleton's pointers D8-LTD (§7 of the reading), and
  the skeleton only.
- `ltd_directions` is the tree itself. Cells are processed in descending elevation. Each takes
  Tarboton's steepest facet, and between the facet's cardinal and diagonal pointers keeps the one with
  the smaller |cumulative transverse deviation|. The deviation is `q × t`, inherited from the donor
  with the largest area. A candidate that is not strictly lower is refused. On a flat, the pointer is
  `compute_flow`'s, with the deviation reset, and the cell is counted.
- `accumulate_cells` is a Kahn accumulation.
- Permanent test `the_ltd_tree_follows_a_planar_slope_off_the_lattice`: on a plane 14.04° off north,
  D8 paths are 14° off (the negative control) and LTD paths within 1.5°; every land cell drains once.
- When off, the definition's world is bit-identical (`identity.txt`).

**C0, as declared in `reading_path_based.md` §7 before the run.** The geometry is a 1:4 plane (TDD
14.04° off the cardinal axis), 30² and 60², sink cells around it. Theoretical areas are the analytic
strip of each cell's drainage lines, on 8 × 8 sub-points. The pass bands, around the sources' ratios:
- MAE ratio in [0.2, 0.9] ([O03] Fig. 3b ≈ 0.44);
- RMSE ratio in [0.05, 0.4] (Fig. 3c ≈ 0.13);
- Paik's cumulative lateral deviation ratio in [0.01, 0.15] ([P08] Table 1: 0.04).

| run | plane | D8: ME / MAE / RMSE | LTD: ME / MAE / RMSE | LTD confluences | ratios MAE / RMSE / lateral | verdict |
|---|---|---|---|---|---|---|
| **1 (as declared)** | 30² | +0.029 / 0.361 / 0.983 | −0.193 / **0.514** / 0.594 | **137** | **1.424** / 0.604 / 0.182 | **FAIL** |
| **1 (as declared)** | 60² | +0.088 / 0.431 / 1.515 | −0.197 / **0.539** / 0.625 | **490** | **1.251** / 0.413 / 0.096 | **FAIL** |
| diagnostic: no 1 000 m offset | 30² / 60² | same | 0.537 / 0.514 (MAE) | 149 / 495 | 1.487 / 1.193 (MAE) | FAIL |
| **2: the tie rule honoured** | 30² | +0.029 / 0.361 / 0.983 | −0.191 / **0.215** / **0.258** | 18 | 0.595 / 0.262 / 0.116 | PASS |
| **2: the tie rule honoured** | 60² | +0.088 / 0.431 / 1.515 | −0.194 / 0.299 / 0.403 | 240 | 0.695 / 0.266 / 0.073 | PASS |
| 2, no offset | 60² | same | −0.194 / 0.224 / 0.268 | 74 | 0.520 / 0.177 / 0.059 | PASS |

**⇒ By the addendum's rule, the first run fails and C and D do not run this round: "c'est le code, pas la
méthode".**

**Where the code fails.**
- The first run's LTD ME (−0.19) is already [O03]'s (≈ −0.21 at λ = 1). That is the analytic error of a
  confluence-free LTD on this plane: each path collects ≈ 0.97 cells per unit of flow length, against a
  strip 1.21 cells wide.
- An MAE of 0.5 therefore means the paths MERGE: 137 cells with two donors on a 30² plane, where D8
  has 0.
- **The declared 1:4 slope is exactly degenerate.** With tan r = ¼, the cumulative deviation reaches an
  EXACT tie at the second step of every four-step cycle (|2 sin r| = |cos r − 2 sin r| ⟺ tan r = ¼).
  [O14] eq. (5)'s "≤" gives that tie to the cardinal.
- In f32, the facet slopes carry a rounding that breaks the tie differently from cell to cell.
  Neighbouring paths fall out of phase, and out-of-phase paths merge. Removing the 1 000 m offset does
  not help (149 / 495), so it is the tie, not the magnitude.

**The fix, applied, and run 2 as its evidence (not as a licence for C).** Ties within 1e-4 cell go to
the cardinal. Confluences fall to 18 (30²), and LTD then reproduces [O03]'s ABSOLUTE values at 30²: MAE
0.215 against ≈ 0.22, RMSE 0.258 against ≈ 0.25, ME −0.191 against ≈ −0.21. All ratios pass.
- Our D8 errors are smaller than [O03]'s (0.36 / 0.98 against 0.50 / 1.9). Their plane's orientation and
  size are not in the text, so the ratios differ while LTD's absolute values agree.
- The residual confluences at 60² with the offset (240, against 74 without) say f32 rounding still
  reaches some ties there.

**C and D are left for the next round**, on the fixed implementation: the two references, then the
gate, then the method, then D if C passes. Nothing of them was measured, so every C and D prediction is
unscored.

## The clause table, updated

| unwanted effect | Finding 127's status | what Finding 128 measured | **Finding 128's status** |
|---|---|---|---|
| **B1 spurs, Finding 38** | mur ↔ mer: PROVED | ON everywhere; the base rate holds in every world (0.0056–0.0142 /km); F38 holds | **PROVED, kept ON** |
| **B1 teeth** | bruit ↔ pied: PROVED for the interaction | ON everywhere | **kept ON**. The profile's own +38 % is named, not treated (the round's instruction) |
| **B2 → A_c canyons** | dammed by a FLOOR, located at the sill | the confluence clause: 4 → 0 and 3 → 0; canyon 13's dam (+246 m) removed; canyon 3's instrument repaired (the lowest exit: +10.6 m) and gone | **PROVED for the canyons. Its Finding 121 control FAILS**: trunk concavity θ +14.9 %, 90 097 moved cells near line heads. Kept gated and off. A narrower form (downstream of the junction only) is named |
| **B2 terrain R8** × 2.3 | "continuous geometry on the D8 tree": REFUTED at calibration | the path-based tree read (B); **our D8-LTD fails its first reproduction (C0, the code: an exact tie broken by rounding)**; with the tie honoured it reproduces [O03]'s absolute values | **NOT MEASURED this round.** C and D stop by the addendum's rule. The implementation is fixed and ready for C |

## Predictions scored

**Mine:**
- **P0 HELD.**
- **B**: P-B1, P-B2 and P-B4 **HELD**. P-B3 **HALF**: Paik's second critique is divergent terrain, not
  the path's start.
- **A**: **P-A1 REFUTED** (4 → 0, not 4 → 1). **P-A2 REFUTED** (the cols stay 91 / 112 m high, but are
  no longer sills). **P-A3 HALF** (slope held, concavity refuted, source cones touched).
- **C0**: **P-C0 REFUTED** (run 1).
- **C and D: NOT MEASURED.**
- **Meta HELD.**

**The reviewer's:**
- **A**: the canyons → 0 **HELD**; the cols back within ±5 m **REFUTED**; slope and concavity ±2 %
  **HALF**; canyon 13 a confluence **HELD**.
- **B**: LTD the right choice **HELD**; "LAD too sensitive on flats" **REFUTED** as the reason; "its
  damping parameter stays PROXY" **REFUTED**; "D∞ failed for lack of memory, which LTD corrects" **HALF**.
- **C and D NOT MEASURED.**
- **Meta HELD.**

**The addendum's amended C** (chord 1: LTD ≈ D8; chords 16–32: 0.10–0.20, independent of the
roughness): **NOT MEASURED.** Its chord-1 half is an identity: every D8-like step is one of eight
directions, so R8 = 1.000 for both.

## Limitations, stated

1. **C0's first run failed on MY declared geometry.** A 1:4 plane is exactly degenerate for LTD (a tie
   every fourth step), and [O03]'s plane is not specified in its text. The fix is the paper's own tie
   rule made robust to rounding (1e-4 cell). Run 2's PASS is the evidence for the fix, not a licence for
   C this round. At 60² with the 1 000 m offset, 240 confluences remain: f32 rounding still reaches some
   ties there.
2. **Our D8 errors are below [O03]'s** (MAE 0.36 against 0.50), so its plane differs from ours; only
   LTD's absolute values agree (0.215 against 0.22).
3. **The concavity θ is a fit over single-cell slopes** (> 10⁻⁴ in both worlds). The témoin's own θ
   reads 0.346 and 0.350 on the two subsets. The +15 % is well beyond that spread, but the metric is one
   number.
4. **"Raised" cells include cells the clause leaves UNCARVED.** The larger line's floor there is above
   the terrain, so they return to the field. They are counted, not attributed further.
5. **The cols of 7 and 14 after the clause are not the sills**, and the new outlets are not located.
6. **Canyon 3's repaired instrument is declared** (the lowest exit on the assembly's flow); it is not
   Finding 119's definition.
7. **One seed, one resolution (8192²).**

## Uncommitted state (awaiting the author)

- `valley_construction.rs`:
  - `ValleyConstruction::trunk_band` (the confluence clause) and `ltd_directions`;
  - `ltd_directions()` (D8-LTD, with the tie honoured) and `accumulate_cells()`;
  - `Skeleton::direction` and `ltd_flat_cells`;
  - two permanent tests with negative controls.
- `cached_product.rs`: the key test lists the two new fields.
- `tests/f126_coast.rs`: Finding 128's benches (`f128_c0`, `f128_a`).
- This report and `reading_path_based.md`.

**Nothing promoted.**
