# Channel-head law `A_c(S)` — what to judge, and what not to look for

ADR 0001, Findings 56 / 56b / 56c. **The law is gated OFF** (`StreamPowerConfig::a_c_slope_law =
None`, byte-identical to the shipped relief). It stays off until this visual judgement is made.

Seed 10 481 999 410 520 546 993 · 400 km domain · production config (`production_hd_config`) ·
`geo_scale_ratio = 1.0`. **Two climate beds**, because the hydrology is a climate result:
`humid` (45°, span 40) is the production default; `arid-hot` (25°, span 10) is the only bed that
produces endorheic below-sea basins. The terrain is climate-independent, so the coastline
figures hold at any latitude.

**Prerequisite, verified mechanically** (`cargo test -p ymir-core --test
channel_head_law_cache_key`): enabling the toggle DOES invalidate the caches — eroded digest
`5a735ac0 → a84ac294`, drainage `9770ded0 → f1d20f47`, and both `S_ref` and `S_min` reach the key
so a recalibration cannot be served stale. Regenerating with the law on cannot silently return
the shipped terrain.

---

## 1. The coastline — the panel AND the absolute count, together

Method rule 9 as amended: a fixed crop chosen on one configuration can render flatteringly clean
for a remedy that *displaces* a defect instead of removing it, so **the panel and an absolute
count must agree** before a visual verdict is accepted.

**Panels** — `exports/coastal_fringes/`, 640×640 at 8192², crop fixed on the densest spur window
**of the SHIPPED panel**:

| file | what it is |
|---|---|
| `law_SHIPPED_8192.png` | the shipped relief — the defect |
| `law_A_c_S__LAW_8192.png` | the law ON |
| `law_REFERENCE_coarse_8192.png` | **the anchor**: what "no comb" looks like (coarse field, no FBM, no incision) |

**Absolute counts on the same fields** (absolute, not per 100 km — the rate has a moving
denominator, which is the trap that cost a round):

| | 2048² | 8192² |
|---|---|---|
| spurs, shipped → law *(reference target)* | 1 254 → **560** *(15)* | 1 849 → **544** *(19)* |
| share of the gap to the reference closed | **56 %** | **71 %** |
| coastline km, shipped → law *(target)* | 4 860 → 3 333 *(1 581)* | 6 206 → **2 692** *(1 602)* |
| local parallelism R *(target)* | 0.515 → 0.513 *(0.330)* | 0.525 → **0.430** *(0.363)* |

**The question.** Is the residual comb in `law_A_c_S__LAW_8192.png` closer to the reference panel
than to the shipped one? The counts say 56 % / 71 % of the way. If the picture disagrees, the
picture wins and the crop should be re-chosen — the defect may have moved.

**Known residual, and its cause.** `S_ref` is 0.3319 at 2048² against 0.1128 at 8192², so **no
single channel-head calibration can hold at both grids** until the hypsometry converges. The
remaining gap is inherited from that, not from the law.

---

## 2. The lakes — look at them, the basins moved more than the coast did

`A_c` moves the hillslope/channel split everywhere, so the basins moved:

| | 2048² off → on | 8192² off → on |
|---|---|---|
| lakes, **humid** (production default) | 80 → **68** (−15 %) | 86 → **76** (−12 %) |
| lakes, arid-hot | 41 → 39 (−5 %) | 56 → **48** (−14 %) |
| water km², humid | 6 847 → 6 670 | 5 971 → **4 981** (−17 %) |
| water km², arid-hot | 3 795 → 3 788 | 1 638 → **1 323** (−19 %) |
| endorheic, humid | 5 → 4 | 2 → 2 |
| endorheic, arid-hot | 27 → **23** | 28 → **24** |
| drainage density km/km², humid | 0.431 → 0.447 (+4 %) | 0.337 → **0.501** (+49 %) |
| drainage density km/km², arid-hot | 0.611 → 0.619 (+1 %) | 0.540 → **0.721** (+33 %) |
| confluences, humid | 1 769 → **1 680** (−5 %) | 2 281 → 3 165 (+39 %) |
| confluences, arid-hot | 3 102 → **3 002** (−3 %) | 5 086 → 6 106 (+20 %) |
| Strahler S5, humid | 36 → 104 | 39 → 99 |
| Strahler S5, arid-hot | 54 → **128** | 51 → **131** |

Read at 2048²: the law **reorganises** the network (order-5 reaches ×2.4 to ×2.7) **without
densifying** it — the confluences actually FALL, by 3 % in arid and 5 % in humid. At 8192² it
does both. **The lake population falls in every bed and at every resolution**, which is the
expected direction (a denser channel network gives more basins an outlet).

**So the survivors are the question**, and it is a visual one: are the 68 / 76 remaining lakes
credible basins, or residue? The metrics cannot answer that.

### Two things to look for specifically

**(a) Rivers ending in nothing.** Orphan mouths — a terminus touching a lake id absent from
`lakes.json` — go **1 → 46 at 2048² and 2 → 26 at 8192²**, in *both* beds. Cause:
`below_sea_basin_lakes_infil` marks every below-sea sink in `lake_map` for sink validity but only
inventories those clearing 4 cells, and the denser network sends far more termini into the
sub-4-cell ones. In the microscope this reads as a watercourse stopping in open ground. A
**reporting** defect in a known place, fixable independently of the channel head — but it must be
closed before the law ships, whatever the visual verdict.

**(b) Lake depths against the height raster, in the HUMID bed.** `depth == level − floor` fails
for **59 of 80 lakes at 2048² and 44 of 86 at 8192²** in humid, against 0–4 in arid. The depths
match the RAW pre-breach field there instead (2–5 failures). Mechanism: H-1c recomputes
`level_m`/`depth_m` on the breached field only for the lakes it re-settles, i.e. the endorheic
ones — and humid has almost none. So **in the production climate most lake depths in
`lakes.json` are pre-breach values read against a post-breach raster.** Pre-existing, unchanged
by the law, and the largest export-coherence problem this suite found.

---

## 3. Navigability — ⚠️ CORRECTED: it DOES move, and downwards, in the production bed

An earlier version of this note said the navigability spread would not move. **That was measured
on the arid bed only and it is wrong for the production default.**

| | 2048² off → on | 8192² off → on |
|---|---|---|
| small boat, **humid** | 214 → **175** (−18 %) | 35 → **11** (−69 %) |
| barge, humid | 2 → 3 | 3 → **2** |
| small boat, arid-hot | 108 → 103 | 2 → 2 |
| barge, arid-hot | 0 → 1 | 0 → 0 |
| discharge p90 m³/s, humid | 1.759 → 1.983 | 1.017 → 1.386 |

**In the humid bed the law cuts navigable reaches by 18 % at 2048² and 69 % at 8192².** The
discharge p90 *rises*, so this is not less water — it is the same water spread over a denser
network, so fewer individual reaches clear a threshold. Whether that is a defect or a truer
picture is a judgement call, but it is not "no change", and barges do exist in humid (2–3).

⚠️ **And ADR Finding 60 gives a reason to distrust the classification itself, independently of
the law.** Navigability classifies on the **climatic discharge**, while the channels it classifies
were carved by the **geometric drainage area** — the incision provably takes no climatic input.
Two quantities that ought to be one, and here they move in opposite directions. So this table may
be measuring the decoupling rather than the channel-head law. **Do not use it to accept or reject
the law**; §1, §2 and §4 are the criteria.

In the arid bed navigability genuinely does not move (108 → 103, 2 → 2), and there the
classification is pinned: the thresholds are in m³/s and the discharge sits 2–3 orders below
them (Finding 47). At 8192² arid **not one of the 51 order-5 reaches carries any discharge at
all**, and the widest channel on the grid is 5.1 m.

---

## 4. The mountain arêtes — what the two-sided form destroyed

This is the check the rejected form failed:

| | 2048² | 8192² |
|---|---|---|
| slope > 45°, shipped → **two-sided** | 1.21 → **0.00 %** | 8.70 → **1.70 %** |
| Strahler S5, shipped → two-sided | 117 → 237 | 59 → **0 (gone)** |
| slope > 45°, shipped → **raise-only (shipped form)** | 1.21 → **1.11 %** | 8.70 → **8.27 %** |
| slope > 30°, raise-only | 14.50 → 12.84 % | 23.68 → 21.56 % |

With the clamp the steep shares move 5–9 % relative and the hierarchy deepens instead of
collapsing. **If the crests look intact, the clamp did its job.** (Note the S5 figures you have
— 59 → 139 — were measured on the unclipped network; on production's clipped network it is
51 → 131. Same story, different absolutes.)

---

## 5. What the metrics already certify, so you need not look for it

- The incision is **alive, and doing 3.2× less work at 8192²** — the right measure is the WORK,
  not the 23 m rise in the mean. Erosion work (pre-incision field minus eroded field, paired
  over the cells that are land in both): **633 m at 2048² against 199 m at 8192²**, lowering
  98.80 % of common land at 2048² but only 66.52 % at 8192². A mean that goes UP is precisely
  the signature already on file for a mechanism that clears nothing, so quoting the rise was
  the wrong evidence even though the conclusion held. The pre-incision reference is 865 m at
  BOTH grids (resolution-invariant to 0.01 %, block A).
- Lake invariants that hold over the real population in both beds: footprint at or below its
  level, monotone long profiles, no duplicate id, no empty footprint, `area_km2` == footprint.
- Three that do not, all pre-existing and none made worse by the law except (a) above:
  **footprint not 4-connected** from its lowest cell for 51–65 of 68–86 lakes (and in humid these
  are almost all EXORHEIC, so H-1c's altitude-sorted truncation does not explain them —
  unattributed); 6–18 `lake_map` ids absent from the inventory, all below-sea sub-4-cell; and
  **2–5 exorheic lakes with no traced outlet in the HUMID bed** (0 in arid), which is a
  mass-balance violation since an exorheic label requires an outlet.

---

## 6. The decision

- **Accept** → the law is un-gated *after* the orphan-mouth regression in §2(a) is closed.
- **Reject** → the coastal target stands at density ≈ 1 per 100 km, p90 ≈ 24 km, local R ≈ 0.33
  (the coarse field's own values), and the next lever is C-4 coastal erosion, which must lower
  the indentation density while raising the p90 length.
- **Undecided on the picture** → say which crop, and it will be re-chosen and re-rendered.

Reproduce: `cargo test -p ymir-core --release --test coastal_comb_levers -- --ignored --nocapture`
(coast + panels) and `--test channel_head_law_invariants` (lakes, network, rule-7, rule-10, both
beds, both resolutions).
