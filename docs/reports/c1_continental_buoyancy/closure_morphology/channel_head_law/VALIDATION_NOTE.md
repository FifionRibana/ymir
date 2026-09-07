# Channel-head law `A_c(S)` — what to judge, and what not to look for

ADR 0001, Findings 56 and 56b. **The law is gated OFF** (`StreamPowerConfig::a_c_slope_law =
None`, byte-identical to the shipped relief). It stays off until this visual judgement is made.

Seed 10 481 999 410 520 546 993 · 400 km domain · production config (`production_hd_config`) ·
`geo_scale_ratio = 1.0` · arid-hot test bed 25°/span 10 for the hydrology.

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
denominator):

| | 2048² | 8192² |
|---|---|---|
| spurs, shipped → law *(reference target)* | 1 254 → **560** *(15)* | 1 849 → **544** *(19)* |
| share of the gap to the reference closed | **56 %** | **71 %** |
| coastline km, shipped → law *(target)* | 4 860 → 3 333 *(1 581)* | 6 206 → **2 692** *(1 602)* |
| local parallelism R *(target)* | 0.515 → 0.513 *(0.330)* | 0.525 → **0.430** *(0.363)* |

**The question.** Is the residual comb in `law_A_c_S__LAW_8192.png` closer to the reference panel
than to the shipped one? The counts say 56 % / 71 % of the way. If the picture disagrees with
that, the picture wins and the crop should be re-chosen — the defect may have moved.

**Known residual, and its cause.** `S_ref` is 0.3319 at 2048² against 0.1128 at 8192², so **no
single channel-head calibration can hold at both grids** until the hypsometry converges. The
remaining gap is inherited from that, not from the law.

---

## 2. The lakes — look at them, the network moved more than the coast did

`A_c` moves the hillslope/channel split everywhere, so the basins moved, and **at 8192² the lake
population moved by more than the coastline did**:

| | 2048² off → on | 8192² off → on |
|---|---|---|
| lakes in the inventory | 41 → 39 | 56 → **48** (−14 %) |
| water surface km² | 3 795 → 3 788 | 1 638 → **1 323** (−19 %) |
| exorheic / endorheic | 14/27 → 16/23 | 28/28 → 24/24 |
| drainage density km/km² | 0.611 → 0.619 (+1 %) | 0.540 → **0.721** (+33 %) |
| confluences | 3 102 → **3 002** (−3 %) | 5 086 → 6 106 (+20 %) |
| Strahler S5 reaches | 54 → **128** | 51 → **131** |

Read at 2048²: the law **reorganises** the network (order-5 reaches ×2.4) **without densifying
it** — the confluences actually fall 3 %. At 8192² it does both.

### One regression to look for specifically

**Rivers ending in nothing.** Orphan mouths — a river terminus touching a lake id that is not in
`lakes.json` — go **1 → 46 at 2048² and 2 → 26 at 8192²**. Cause: `below_sea_basin_lakes_infil`
marks every below-sea sink in `lake_map` for sink validity but only inventories those clearing
4 cells, and the denser network sends far more termini into the sub-4-cell ones. In the
microscope this reads as a watercourse stopping in open ground.

This is a **reporting** defect in a known place, not a terrain defect, and it is fixable
independently of the channel head. But it must be closed before the law ships even with a
favourable visual verdict.

---

## 3. Do NOT look for barges — the navigability spread will not move

| | 2048² off → on | 8192² off → on |
|---|---|---|
| small boat | 108 → 103 | **2 → 2** |
| barge | 0 → **1** | 0 → 0 |
| ship | 0 → 0 | 0 → 0 |
| discharge m³/s p90 | 0.283 → 0.205 | 0.005 → 0.021 |

The navigability classes are thresholds in m³/s and the discharge sits two to three orders of
magnitude below them (Finding 47). **The classification cannot respond to a channel-head law** —
only to the hypsometry.

How dry it is: **at 8192², not one of the 51 order-5 reaches carries any discharge at all** in
the shipped configuration (the law raises that to 6.9 % of 131 reaches). Across the network only
14–25 % of reaches are wet at all, and the widest channel anywhere on the fine grid is 5.1 m.
This is the arid test bed doing what it should — a desert continent has no exotic trunk rivers —
but it is also why no width or navigability calibration can be anchored here yet.

So the barge/ship question is untouched by this change. It is blocked on the discharge, hence on
the hypsometry, which the roadmap now ranks as priority 1 ahead of everything that depends on
slope (`docs/roadmap_closures.md`).

---

## 4. What the metrics already certify, so you need not look for it

- The incision is **alive** — hypsometry mean 282 → 305 m at 2048² and 685 → 697 m at 8192²
  (raw eroded field), against 865 m un-eroded. That is the control block that caught
  `A_c × 100`, where the coast looked cured because erosion had stopped happening.
- The **arêtes survive**: slope > 45° 1.21 → 1.11 % at 2048², 8.70 → 8.27 % at 8192². The
  two-sided form of the law erased them (1.21 → 0.00 %, and order 5 vanished at 8192²) and was
  rejected for it; the shipped form is clamped **raise-only**.
- Lake invariants that hold over the real population: footprint at or below level, exorheic
  implies a traced outlet, monotone long profiles, no duplicate id, no empty footprint, and lake
  depths coherent with the **exported** height raster.
- Two that do not, both pre-existing and unchanged by the law: **more than half of every lake
  population has a footprint that is not 4-connected** from its lowest cell (partly H-1c's
  altitude-sorted truncation, partly unattributed), and 6–18 `lake_map` ids are absent from the
  inventory (all below-sea sub-4-cell).

---

## 5. The decision

- **Accept** → the law is un-gated *after* the orphan-mouth regression in §2 is closed.
- **Reject** → the coastal target stands at density ≈ 1 per 100 km, p90 ≈ 24 km, local R ≈ 0.33
  (the coarse field's own values), and the next lever is C-4 coastal erosion, which must lower
  the indentation density while raising the p90 length.
- **Undecided on the picture** → say which crop, and it will be re-chosen and re-rendered.

Reproduce: `cargo test -p ymir-core --release --test coastal_comb_levers -- --ignored --nocapture`
(coast + panels) and `--test channel_head_law_invariants` (lakes, network, rule-7, rule-10).
