# Finding 129 — the construction lays its law exactly where it cuts, and the confluence clause bends it there; the source cones are moved both ways (A2 not triggered); our D8-LTD misses the reserved case on one ratio, so C and D stop

**Nothing was built for production and nothing is promoted.**
- The confluence clause stays gated and off. Its activation is the author's decision, with A in hand.
- The benches are in `crates/ymir-core/tests/f126_coast.rs`, Finding 129's section. The raw outputs are
  beside this report.
- **mur ↔ mer is ON in every world.** bruit ↔ pied is a clause of the wall PROFILE. These are Finding
  128's B2 → A_c worlds, with planar walls, so bruit ↔ pied has nothing to act on (declared).

**Units.**
- 1 cell = 48.8 m (domain), 8192² over 400 km.
- Heights and drops are in metres. Areas are domain km².
- θ is the exponent of S ∝ A^−θ.

**Reading declaration, unfavourable.** My predictions (`f129_predictions.md`) were written before any
grep or measurement. They are non-blind on Findings 73–128 and on the reviewer's predictions.
- The instrument and B1's case were declared before any measurement.
- **Two A benches were added after reading the first stage values.** They were declared non-blind, and
  their predictions were written before their outputs:
  - `f129_a_law`: the law's reference measured, and where the as-built leaves it;
  - `f129_a_cones`: the cones, before and after separated.
- They add references and a partition, not variants.

## 0 — Finding 128 committed; the guard; the greps

- **Finding 128 was committed** as `5effc92`. The guard through `run_hd` reads **6 / 6 "= banc"** (`guard_six_states.txt`).
- `docs/refs/*.pdf` stay uncommitted (rights).
- **This round's library change is `LtdStats`.** `ltd_directions` now wraps `ltd_directions_stats`
  byte-identically, and it is used by the skeleton only under `ltd_directions` (off).
- **After every change of this round**, the guard reads 6 / 6 "= banc" again (`guard_after_this_round.txt`).
  The fresh build of the definition hashes to `a8d2d538d692c2f0`, the reference (`identity.txt`).

**The rule-11/11b/11c greps (`grep_rule11.txt`).** The earliest hit is the one read.
- **"source cone" / "cône de source": NOTHING FOUND** in the dossier.
- **concavity** first appears in Finding 128 (its control). **slope-area** appears once, in Finding
  128's "θ (S ∝ A^−θ) goes 0.350 → 0.402". So the dossier never measured a concavity before Finding 128.
- **Harel** first appears in Finding 89, as the B6 duration reference. The construction's law cites
  Harel 2016 eqs. (4)–(7), PDF p. 15 (`valley_construction.rs`, module table).
- **Paik** and **false confluence** appear only in Finding 128. **dome** appears first in Finding 126
  D1's isotropic synthetics.
- **Finding 121's trap is at ADR L17966**: "a minimum over every cone flattened the long profile". That
  is why the construction lays each cell from the NEAREST sample.

## A — θ by stage, and what the clause does to the cones

**The instrument, declared before any measurement.**
- θ is the slope-area OLS slope, −d ln S / d ln A.
- It is fitted over the TRUNK ≥ 10 km² cells of the construction-input skeleton (S1). The A range runs
  from 10 km² to the largest trunk.
- S is the stage's drop to the cell's D8 receiver on that skeleton, over the link length. Only cells with
  S > 10⁻⁴ in that stage are kept.
- Finding 128's own form (the PRE skeleton, cells with S > 10⁻⁴ in both worlds) is reported beside it,
  for continuity.
- **Harel 2016's 0.51 ± 0.12 is a χ-integral estimate**, not a slope-area one. It is a reference of the
  same quantity by another estimator.

### References, measured first (rule 15)

| reference | θ | cells |
|---|---|---|
| **the instrument on the law itself** (z = base + k·χ on the S1 skeleton, k = 0.07183) | **0.500** (S > 0: 0.500) | 70 720 links, none across two bases |
| Harel et al. 2016, ⟨m/n⟩ (PDF p. 23; 0.51 ± 0.14 at p. 36) | 0.51 ± 0.12 | — |
| S1, the construction's input | 0.274 | 58 791 |
| PRE (no incision, bathymetry) | 0.274 | 58 791 |
| the delivered (livré) | 0.259 | 46 276 |
| A1+B2 | 0.300 | 50 107 |

- **The instrument reads the law exactly.** The 0.5 reference is now measured, not written.
- **S1 and PRE are the same field on land.** Production ships no droplet pass (`upscale.rs:410`,
  `cfg.erosion = None`), and the bathymetry only touches sub-sea cells. Their identity is expected, not
  a coincidence.
- The shipped incision (0.259) and its closure (0.300) sit near the input's own 0.274, far below the law
  and below Harel's band.

### The stages (`a_theta_by_stage.txt`)

| world | (i) as built | (ii) + drainage / breach | (iii-a) + the light pass | (iii-b) the world |
|---|---|---|---|---|
| **témoin C2/10 col** | **0.338** (65 802) | 0.369 (63 054) | 0.351 (62 950) | 0.351 |
| B2 → A_c, no clause | 0.335 (62 621) | 0.367 (59 870) | 0.364 (59 257) | 0.364 |
| **B2 → A_c + CLAUSE** | **0.317** (62 430) | 0.346 (59 666) | **0.401** (61 031) | **0.401** |

- **(iii-a) = (iii-b) on land by construction.** There is no droplet pass, and the bathymetry is
  sub-sea only.
- **Finding 128's paired form reproduces to the digit.** It reads 0.346 → 0.366 and 0.350 → 0.402.
  These are Finding 128's worlds.

**Where the as-built construction leaves its law (`a_law_reference.txt`).**
- These are the 70 904 trunk ≥ 10 km² cells, from `carve` on S1, the pipeline's own call.
- (i) is reproduced: 0.338 / 0.335 / 0.317.

| world | as built − law, p10 / p50 / p90 | > 1 m above the law | > 1 m below the law | UNCARVED (min kept the field) | θ, links lowered at both ends | θ, the other links | θ, links within ±1 m of the law at both ends |
|---|---|---|---|---|---|---|---|
| **témoin C2/10 col** | −192.2 / −0.0 / +0.0 m | 82 | 23 937 (33.8 %) | 22 776 (32.1 %) | **0.493** (46 885) | 0.161 (18 917) | 0.497 (45 593) |
| B2 → A_c, no clause | −192.2 / −0.2 / +0.0 m | 86 | 26 608 (37.5 %) | 22 778 (32.1 %) | **0.493** (43 699) | 0.161 (18 921) | 0.485 (39 443) |
| **B2 → A_c + CLAUSE** | −192.2 / −0.0 / +0.0 m | **1 295** | 26 528 (37.4 %) | 22 794 (32.1 %) | **0.469** (43 491) | 0.161 (18 939) | 0.483 (38 860) |

- **Where the construction cuts, it lays its law** (0.493 against the instrument's 0.500 on the law).
- **(i)'s 0.338 is a mixture.** A third of the trunk ≥ 10 km² cells are UNCARVED. There the input field
  lies BELOW its law (22 775 of 22 776), so `min(field, V)` keeps the field; those links read 0.161.
- The construction does not act on that third, by its own rule ("it can only lower a cell"). This is
  named, not treated.
- **B2's sub-valleys cut 2 671 more trunk cells below their law** (the cross-line minimum). The laid
  links still read 0.493.
- **The clause bends the laid floor, 0.493 → 0.469.** On the links within ±1 m of the law at both ends,
  θ is unchanged (0.485 → 0.483).
  - What moves is the cells the clause pulls OFF the law: 1 209 more trunk cells sit more than 1 m above
    their own law.
  - They can only be cells of a larger line's floor band, by the clause's definition, laid from that
    line's nearest sample.
  - **Where they sit relative to the junction: NOT measured.**
  - This is Finding 128's "step at the band's edge", measured for the first time as a count and a θ.

**The pre-written reading, applied.**
- The round's rule: if the no-clause floor is already at 0.5 at (i) and the clause moves it away, "**la
  clause déforme**, quel que soit le rapprochement de Harel aux étages suivants".
- The no-clause constructed floor is at the law: 0.493, and 0.485 within ±1 m.
- The clause moves it away: 0.469, and the whole trunk 0.335 → 0.317.
- **⇒ The clause deforms the law at construction.** The restoring branch is refuted: at no stage does
  the clause bring (i) toward 0.5.
- **The world's +0.037 comes from the light pass** acting on the clause's field (0.346 → 0.401). The
  other worlds move by −0.018 and −0.003 over that stage. Why is NOT measured.
- The +0.037 places the clause's world inside Harel's band (0.39–0.63) where the témoin (0.351) is not.
  **By the pre-written reading, that does not count for the clause.**

### The source cones (`a_theta_by_stage.txt`, `a_cones_before_after.txt`)

**The census on the pipeline's S2, B2 → A_c.**
- The clause moves 388 282 cells; **90 069 lie within 3 cells of a line's head** (Finding 128: 90 097,
  on PRE).
- **Δz: RAISED 55 185 (61 %), lowered 34 884.** p10 −51.2 / p50 +0.5 / p90 **+13.9 m**.
- **Position:**
  - on the larger line's floor band: 82 895 (**92.0 %**);
  - left uncarved, back at the field: 7 163;
  - on a wall: 11.

**The touched cones' profiles** (head → 10th cell, 10 cells = 488 m). There are 3 956 lines ≥ 11 cells
with a moved cell among their first four, out of 103 075 lines (3.8 %). The before / after split is on
`carve(S1)`, with the census's 3 953 / 200 reproduced.

| | before (no clause) | after (clause) | change |
|---|---|---|---|
| **MEAN** slope | 0.1017 m/m | 0.1127 m/m | **+10.8 %** |
| p50 | 0.1174 | 0.1122 | −4.4 % |
| p25 | 0.0500 | 0.0407 | **−19 %** |

| what the clause does to a touched line | lines |
|---|---|
| reverses its head (drop < 0 after, ≥ 0 before) | **67** (16 with a lowered head) |
| un-reverses it (< 0 before, ≥ 0 after) | **158** (291 reversed before, 133 still) |
| flattens it below half | **341** |
| steepens it above × 1.5 | 291 |

- The head cell itself does not move: Δz p10 −3.6, p50 0.0, p90 +61.2 m. 61 % of the heads are on a
  floor band.
- **The clause's reversals come from a RAISED downstream reach**, not a lowered head (51 of 67 heads
  are not lowered). It is the same raising as on the trunks.
- 5 388 of the 21 514 lines shorter than 11 cells are touched. They are not profiled; their head Δz is
  p10 −10.2, p50 0.0, p90 +1.7 m.

**Finding 121's trap (a flattening near sources), looked for explicitly: NOT found as a systematic
effect.**
- The mean rises 10.8 %, and the clause un-reverses more heads than it reverses (158 against 67).
- The lower quartile flattens 19 %, and 341 lines are halved. That is a flattening of the flattest
  quarter, offset elsewhere.

**⇒ A2 NOT triggered.** The round did not quantify "the cone cells are a defect", and neither did I
before the data. So this is my reading after the data, stated as such:
- the cones are moved both ways, with a net fewer reversed heads;
- the mean steepens;
- 92 % of the moved cells lie on the larger line's band, where the clause is meant to act.

**The measured deformation is on the TRUNKS** (the laid floor 0.493 → 0.469, the 1 209 cells raised off
the law), not on the cones. A2's restricted form was designed for the cones. Whether it would remove
the trunk deformation depends on where the 1 209 cells sit relative to the junctions, which is not
measured. **The author decides.**

## B1 — D8-LTD on the reserved case: FAIL, narrowly, on one criterion of four; C and D stop

**The case, declared before the run** (`f129_predictions.md`).
- It is [O03]'s PARABOLIC VALLEY (Fig. 2d, Fig. 3j–l, TNN 1-4/1-5), not the plane that found Finding
  128's tie bug.
- The geometry is not in [O03]'s text. Declared:
  - 30² interior cells of 10 m, with sink cells around them;
  - the flow along +x to the east edge;
  - z = (N − x) + (y − yc)² / 30 in cell units, so the cross gradient at the valley's side edge equals
    the longitudinal one.
- The theoretical areas come from the exact descent lines of 8 × 8 sub-points per cell.
- The paper's values, read on the render: D8 (λ 0) MAE ≈ 0.57, RMSE ≈ 2.1; LTD (λ 1) MAE ≈ 0.29, RMSE
  ≈ 0.35.
- **PASS iff all four hold**, else STOP before C:
  - LTD MAE in [0.19, 0.39];
  - LTD RMSE in [0.20, 0.60];
  - MAE ratio in [0.3, 0.8];
  - RMSE ratio in [0.08, 0.40].

| | ME | MAE | RMSE | [O03] Fig. 3j–l |
|---|---|---|---|---|
| D8 | −0.227 | 0.647 | 1.092 | ME ≈ 0.0 · MAE ≈ 0.57 · RMSE ≈ 2.1 |
| **D8-LTD** | −0.283 | **0.368** ✓ | **0.440** ✓ | ME ≈ −0.26 · MAE ≈ 0.29 · RMSE ≈ 0.35 |
| **ratio LTD / D8** | | **0.569** ✓ (paper ≈ 0.51) | **0.403** ✗ (band ≤ 0.40; paper ≈ 0.17) | |

- **⇒ FAIL.** One criterion of four misses, by 0.003. LTD takes no flat fallback here.
- **I do not re-run it on another geometry.** A second valley chosen after seeing this one would be a
  forking path.
- **By the round's rule, C and D do not run, not even C's references.**
- **Reading, not a verdict.**
  - LTD's own errors are inside their bands but +27 % (MAE) and +26 % (RMSE) above the paper's.
  - Our D8 RMSE is half the paper's (1.09 against ≈ 2.1). The ratio fails because our D8 is better than
    [O03]'s on our geometry, not because our LTD leaves its own band.
  - The gate was declared on the ratio, and the ratio fails.

## B2 — the ties on the real pre-incision (`b2_ties_pre_incision.txt`)

The surface is the one the skeleton reads: S1 breached by the pipeline's own drainage.

| quantity | value |
|---|---|
| land cells | 11 327 660 |
| D8 exact ties (the two steepest slopes equal, f32) | 6 124 (**0.054 %**) |
| D8 ties within 1e-6 relative | 6 152 (0.054 %) |
| D8 exact ties after the u16 export (0.115 m / code) | 15 316 (**0.135 %**) |
| LTD choices the tie rule decides (\|c1\| vs \|c2\| within 1e-4 cell) | 1 407 (**0.012 %**), 77 exact in f32 |
| where they sit, by the land's slope quartile (Q1 ≤ 0.040 · Q2 ≤ 0.088 · Q3 ≤ 0.161 m/m) | 23.1 / 26.0 / 26.2 / 24.7 % |
| **LTD flat fallbacks** (no strictly lower neighbour: `compute_flow`'s pointer) | **684 893 (6.05 %)** |

- **Ties do not matter on the real terrain.**
  - The tie rule decides one LTD choice in 8 000, uniformly across the slope quartiles, not on smooth
    terrain.
  - The u16 export multiplies D8's ties by 2.5 and keeps them below 0.2 %.
  - Finding 128's bug belonged to an exactly degenerate synthetic, not to this surface.
- **The number that matters is elsewhere.** 6.05 % of land cells have no strictly lower neighbour on the
  breached field. LTD takes `compute_flow`'s pointer there, and its flat gradient is Garbrecht–Martz and
  cardinal (Finding 128 B, cause 2).
- Any LTD skeleton would carry that share of D8 lattice. Named, not treated.

## C and D — NOT RUN

"Pas de C si B1 échoue ; pas de D si C échoue."
- **No gate for C was declared**, because its references were not measured (rule 15).
- **Nothing of C or D was measured**: not the free continuous trace, LTD's chords, Paik's dome, or the
  skeleton on LTD directions.
- The bench `f129_c_refs` (both references, chords {1, 8, 16, 32}, the dome) is in the file,
  `#[ignore]`, and was NOT run.

## The clause table, updated

| unwanted effect | Finding 128's status | what Finding 129 measured | **Finding 129's status** |
|---|---|---|---|
| **B1 spurs, Finding 38** | mur ↔ mer PROVED, kept ON | ON in every world | **PROVED, kept ON** |
| **B1 teeth** | bruit ↔ pied kept ON; the profile's +38 % named | not re-measured (planar walls here) | **kept ON**, unchanged |
| **B2 → A_c canyons** (the confluence clause) | PROVED for the canyons; its concavity control FAILS | see the next four rows | **DEFORMS the law at construction** (the pre-written reading). Kept gated and off; activation is the author's |
| — the law at construction | — | the instrument reads 0.500 on the law; the laid floor is at it (0.493); the clause bends it to 0.469 and raises 1 209 trunk cells > 1 m above their law | position relative to the junction NOT measured |
| — the world's +15 % | concavity 0.350 → 0.402 | made by the light pass on the clause's field (0.346 → 0.401) | why: NOT measured |
| — the source cones | 90 097 moved cells near heads | 92 % on the larger band, 61 % raised, p90 +13.9 m; mean slope +10.8 %, p25 −19 %; 67 heads reversed against 158 un-reversed | **not a defect in Finding 121's sense (my reading): A2 NOT triggered** |
| **the construction's uncarved third** (new, named) | — | 32.1 % of the trunk ≥ 10 km² cells lie below their law; `min` keeps the field (θ 0.161); hence (i) = 0.338 | **named, not treated** |
| **B2 terrain R8** × 2.3 (LTD) | NOT MEASURED; the implementation fixed | B1 FAILS (RMSE ratio 0.403 against ≤ 0.40); ties negligible (0.012 %); flat fallbacks 6.05 % | **NOT MEASURED.** C and D stop by the round's rule |

## Predictions scored

**Mine** (`f129_predictions.md`):
- **P0 HELD.**
- **A**
  - **P-A0a HELD**: the law reads 0.500.
  - **P-A0b HALF**: the gap is the uncarved third (held); it reads 0.161, not ≈ 0.27 (refuted); the
    laid links read ≥ 0.45 (held).
  - **P-A1 HALF**:
    - "without the clause ≈ 0.45–0.50 at (i)": REFUTED on the declared instrument (0.335);
    - "the clause moves it away": HELD; "to ≥ 0.52": REFUTED (it moves down);
    - "the clause distorts": HELD;
    - "the drop comes at (ii)–(iii)": REFUTED (it is at (i); (ii) raises every world).
  - **P-A2 HALF**: pre-incision and Harel held; the delivered and A1+B2 "0.35–0.45" refuted (0.259 /
    0.300).
  - **P-A3 HALF**: raised, p90 ≥ 10 m and in the band all held; "the mean slope falls, A2 triggered"
    REFUTED (the mean rises 10.8 %).
  - **P-A3b REFUTED**: 67 of 200 reversals are the clause's; 16 of those 67 have a lowered head.
- **B**
  - **P-B1 REFUTED** (FAIL on the RMSE ratio).
  - **P-B2 HALF**: D8 ties < 1 % held; LTD 1–5 % refuted (0.012 %); "on smooth terrain" refuted.
- **C and D NOT MEASURED.**
- **Meta ("at least two wrong") HELD.**

**The reviewer's:**
- **A**
  - "(i) no clause 0.40–0.43": REFUTED (0.335).
  - "with the clause 0.48–0.50, it restores the law": REFUTED (0.317, away).
  - "(ii) and (iii) both fall as much": REFUTED. (ii) raises both; (iii) raises the clause's by 0.055.
  - Cones "Δz mostly negative and small (p90 ≤ 5 m)": REFUTED (61 % raised, p90 +13.9 m).
  - "in the trunk's band": HELD (92 %).
  - "not Finding 121's flattening": HELD on the mean and the median, with the lower quartile's −19 %
    named.
  - "A2 not triggered": HELD, on my reading.
- **B1** ("within ±10 % of the paper"): REFUTED (LTD +27 % / +26 %).
- **B2** ("2–8 % of land cells tie, on smooth terrain"): REFUTED on both halves.
- **C and D NOT MEASURED.**
- **Meta ("at least one of the four wrong") HELD.**

## Limitations, stated

1. **B1's geometry is mine.** [O03] does not give its valley, and its values were read on a render.
   - Our D8 RMSE is half the paper's, so the failing ratio has our D8 in its denominator.
   - The code builds z in metres at 10 m per cell, where the declaration says 1 m per cell. D8 and
     LTD's pointers are invariant under a uniform positive scaling of z: the slope ORDER, Tarboton's
     facet angle (a ratio), the strictly-lower test, and a deviation counted in cells. So it is the
     declared case.
2. **θ is an OLS over single-link slopes** of spatially correlated cells, with no confidence interval.
   - Differences ≤ 0.005 are not read.
   - The clause's shifts at (i) (−0.018 on the whole trunk, −0.024 on the laid links) stand against a
     0.000 difference between the témoin and B2 on the laid links.
3. **The partitions of `f129_a_law` were chosen after seeing (i).**
   - The partitions are: lowered at both ends / the others / within ±1 m of the law.
   - They were declared as such, with their predictions written first.
   - The ±1 m tolerance is mine; the lowered-at-both-ends partition needs none.
4. **The 1 209 cells raised off the law are located only by the clause's definition** (on a larger line's
   floor band). Their position relative to the junction is not measured, and that is what would tell
   whether A2's form reaches them.
5. **The light pass's +0.055 on the clause's world** is measured, not explained.
6. **A2's trigger was not quantified**, by the round or by me, before the data. The call is mine, and
   the author's to overrule.
7. **The cones.**
   - Lines shorter than 11 cells (5 388 touched) are not profiled.
   - The census is on S2, the before / after on `carve(S1)`. They differ by the active rims only, and
     reproduce each other: 3 953 / 3 956 lines, and 200 / 200 reversals.
8. **The machine was suspended during `f129_a`.** Its 9 276 s for the témoin is not a cost.

## State

**Uncommitted**, awaiting the next round's go-ahead.
- `valley_construction.rs`: `LtdStats` and `ltd_directions_stats`.
- `f126_coast.rs`: Finding 129's benches.
- This report, and ADR Finding 129.

**Checks after the last change:**
- `cargo test -p ymir-core --release --lib`: 586 passed, 0 failed, 6 ignored.
- `cargo check --workspace --tests --release`: clean.
- The guard reads 6 / 6, and the identity hash is `a8d2d538d692c2f0`.
