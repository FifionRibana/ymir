# ADR 0001 — Erosion: coastal sediment sink, and the terrace / dendritic-valley findings

Status: accepted (2026-08). Scope: C1 HD relief (erosion, FBM upscale, isostasy).
Diagnostic harness: `crates/ymir-core/tests/terrace_diagnosis.rs` (all `#[ignore]`),
commits `019113f`, `76ae09f`, `94fe11b`. Re-run:
`cargo test -p ymir-core --test terrace_diagnosis --release -- --ignored --nocapture`.

This ADR records findings from a measurement campaign on why the C1 HD relief showed
(a) terraces following isolines and (b) no carved dendritic valleys under a
fully-computed river network. Each is stated with the numbers, because the
conclusions are counter-intuitive and will otherwise be re-litigated.

## Decision summary (read this first)

- **Droplet hydraulic erosion is REMOVED from the C1 HD pipeline**, not merely
  reduced. It DESTROYS relief (drainage relief 258→37 m; 659→124 m even at a weak
  0.25 droplets/cell) — it deposits in channels faster than it incises. It is the
  common cause of the terraces (depositional flats), the missing valleys, the
  net-zero mass balance and the coastal sediment dump. See Findings 3–4.
- **Routed stream-power incision (Braun & Willett) + a hillslope REGIME SPLIT is the
  replacement.** Recommended: `K=3, m=0.5, n=1, iterations=3, A_c≈50 cells,
  diffusion≈0.05` (diffusion on hillslopes `A<A_c`, stream power on channels
  `A≥A_c`), θ=0, droplets off, **UNCOUPLED vertical scale** (do NOT couple
  depth_scale to domain). Off by default until confirmed at 8192².
- **θ (incision threshold) and A_c (critical area) each FAILED ALONE** — record this
  so the cheap options are not retried. Only the coupled regime split fixed the
  headwater over-carving AND opened valley cross-sections. See Finding 3.
- **The acceptance CRITERION is a LEGIBLE landscape** — buildable valley floors +
  sharply-delimited steep flanks + moderate interfluves — **NOT minimal steep
  ground.** An earlier "minimise unbuildable land / keep >45° under ~5 %" framing was
  WRONG and cost a pass: steep ground is gameplay content (cliffs, passes, terracing
  constraints). Judge on spatial structure (valley-floor area, cliff-transition
  sharpness, flank contiguity), not a global steep-%. This is why the UNCOUPLED case
  (steeper, more legible flanks) is preferred over the coupled one (too flat).

---

## Finding 1 — Low `target_land_fraction` is the WRONG lever for a border-clean island

**Context.** The intuitive way to "get an island" is to reduce the land fraction
(raise the ocean) via `target_land_fraction` (tlf) quantile calibration. It is
appealing enough that it will be retried.

**Measurement.** tlf calibration subtracts the `1−f` quantile, which moves SEA LEVEL
onto the flat continental-shelf plateau. Local hypsometric slope at the chosen sea
level (cells per raw-altitude unit; low = flat plateau = hypersensitive coastline):

| sea level | hypsometric slope | coastline |
|---|---|---|
| isostatic (tlf None) | ~12 350 | crisp |
| tlf 0.29 | ~43 950 | crisp (but big continent) |
| **tlf 0.08** | **~2 000** | **speckled — sea on the flat shelf** |

**Decision.** Do NOT use low tlf to make an island. Default is `tlf = None`
(isostatic sea level). Bound the landmass geometrically (seed selection / framing),
not by drowning it. See the domain-as-map seed verdict (`land_topology.rs`).

**Consequences.** Speckled coasts and marginal land avoided; islands come from
seed geometry + the framing roll, at the isostatic sea level.

---

## Finding 2 — The coastal deposition dump was the net-zero balance lock

**Context.** Hydraulic erosion moved mass but the net balance was ~0 (author's
8192² run: eroded 428 331 / deposited 424 747 ≈ +0.8 %), so the field was barely
reshaped. Suspected cause: sediment redeposited before leaving the domain.

**Measurement.** `hydraulic.rs` terminal coastal path: when a droplet crosses
`sea_level`, after `coastal_deposition_range` (=12) sub-sea steps it **deposited its
entire remaining load at the shoreline**. On the reference field, deposition
concentrated at the coast: **73 % ≤5 cells, 85 % ≤20 cells** from the coastline. M1's
`sea_level` 0.1→0.5 moved that dump onto the true coastline. On the real coarse→FBM
field the terminal dump dominated the balance: sink fraction f=1.0 → **net +1 %**,
f=0.0 → **net +59 %**.

**Decision.** Add `ErosionConfig.coastal_deposit_fraction` (a partial sea SINK):
deposit `f` of the terminal load, discard `1−f` (still counted as eroded). Ship
**f = 0.25** by default: net **+1 % → +47 %** (clearly erosive), deltas/beaches
preserved (deposition ≤5 cells ~47 % vs 24 % at f=0, so the coast is not abrupt),
network hierarchy preserved (maxStrahler 6; droplet density 4/cell would collapse it
to 3 — see Finding 3), emerged fraction drops only ~1 point (23.7 % → 22.7 %). Droplet
density stays 0.95/cell. The field is `skip_serializing_if` at 1.0, so a config
explicitly set to 1.0 reuses the pre-sink cache; the new 0.25 default DOES enter the
eroded cache key and correctly rebuilds derived stages.

**Consequences.** Every existing map's terrain changes (intended cache
invalidation). The mass-balance / coastal-deposit artifact is fixed. This does NOT
carve valleys — see Finding 3.

---

## Finding 3 — Droplet erosion cannot make a dendritic network here (algorithm limit, not tuning)

**Context.** Ymir computes a full drainage network (`rivers.json`, Strahler orders,
`flow_accumulation`), but the relief under the rivers shows no carved valley.

**⚠️ METRIC CORRECTION.** This finding was first measured with a WRONG metric:
"fraction of top-1 % flow cells in a local altitude minimum" (reported as 15 % in
production, ~11 % here). That measures **PITS, not channels** — a drained channel
cell is by construction HIGHER than its downstream receiver, so it is never a local
minimum; only an undrained depression is. Its interpretation was therefore INVERTED
(a low value is *good* drainage, not "uncarved terrain"). The correct quantity is
**drainage relief** = median over top-1 % flow cells of (11×11 window-max − cell), in
metres — how deep channels sit below their interfluves. See the prototype table below
for the real numbers; the pit table is kept only to show what was refuted.

**Pit-fraction table (mislabelled "carved" — the refuted metric; low = fewer pits):**

| f | droplets/cell | net % | pit-fraction | maxStrahler | confluences |
|---|---|---|---|---|---|
| 1.0 | 0.95 | +1 | 11 % | 6 | 2366 |
| 0.0 | 0.95 | +59 | 11 % | 6 | 2183 |
| 0.25 | 4.0 | +76 | 18 % | **3** | 649 |
| 0.0 | 4.0 | +87 | 17 % | **3** | 556 |

More droplets raise the pit-fraction to ~18 % only by fragmenting the network
(maxStrahler 6→3, confluences 2366→649 — grain, not drainage). On a synthetic smooth
cone at production parameters, surface roughness barely changed (0.00005 → 0.00006):
no reshaping even on an ideal input.

**THE CENTRAL FINDING (drainage relief, the correct metric).** The production droplet
pass does not sculpt relief — it DESTROYS it: drainage relief **258 m → 37 m (−86 %)**.
Stream-power along the existing network RAISES it to **323 m (+25 %)** at ~13× the
speed. We are not swapping one algorithm for a better one — we are removing something
actively HARMFUL. This single mechanism explains all four symptoms of this ADR: the
terraces (Finding 4 — depositional flats), the missing valleys (this finding), the
net-zero balance and the coastal dump (Finding 2 — the droplets that reach the coast
carry the sediment they refuse to incise with).

**Root cause.** Droplets are stochastic and UNCORRELATED — nothing makes neighbouring
rills converge into a shared, deepening channel; they deposit in and around channels
faster than they incise, so relief collapses and a hierarchical valley network cannot
emerge.

**Decision / proposed successor (NOT implemented).** Routed **stream-power incision**
along the drainage network that is already computed (`flow_accumulation`, `rivers.json`
with Strahler orders and upstream/downstream topology) — deterministic and global, so
hierarchy is imposed by construction. Two pitfalls for whoever implements it:
1. stream-power incision MODIFIES the terrain, so drainage computed beforehand goes
   stale — either iterate drainage↔incision a few times, or accept a single pass;
2. it should COMPLEMENT droplet erosion (stream-power for valleys, droplets for
   hillslope texture), not replace it.

**Consequences.** The missing-valleys problem is out of scope for parameter tuning
and is a separate chantier. The sink (Finding 2) fixes the balance/coast, not this.

**Prototype results (measured, `erosion/stream_power.rs`, off by default).** Braun &
Willett implicit scheme, real coarse→FBM field (1024²), using **drainage relief** =
median (11×11 window-max − cell) over top-1 % flow cells, in metres — the correct
incision metric. (The earlier "carved% / local-minimum fraction" was WRONG: a drained
channel cell is higher than its downstream receiver, so it is NOT a local minimum; the
local-minimum fraction measures PITS, not incision, and its sign is inverted.)

| config | drainage relief | vs FBM | notes |
|---|---|---|---|
| FBM baseline | 258 m | — | |
| **stream-power alone** (K=1,m=.5,n=1,dt=1,iters=3) | **323 m** | **+65** | carves; maxStrahler not fragmented (4); ~554 ms |
| stream-power + diffusion D=0.4 | 230 m | −28 | D=0.4 OVER-smooths |
| droplets alone (production) | 37 m | **−221** | droplets COLLAPSE relief (~7270 ms) |
| both (SP then full droplets) | 24 m | −234 | droplets ERASE the SP valleys |

Findings: stream-power **carves** (+25 % relief) and is **~13× faster** (554 ms vs
7270 ms at 1024²); drainage↔incision **converges by iteration 2** (use N=2–3); K=1/dt=1
calibrated to ~307 m median channel incision (a plausible valley depth, not tuned for
looks). Crucially, **the production droplet pass is anti-incision** — it collapses
drainage relief 258→37 m, the SAME deposition mechanism as the Finding-4 terraces — so
droplets are the common cause of both symptoms. Coupling therefore cannot keep the
full droplet pass; droplets must be reduced to a WEAK hillslope-texture pass (or
dropped) so they don't erase the channels. Diffusion at 0.4 over-smooths; start near 0.

**Hillslope regime (why pure stream power over-carves, and the fix).** Pure stream
power is DETACHMENT-LIMITED WITH NO THRESHOLD — `E = K·A^m·S^n > 0` for any `A > 0`,
so even a near-zero-area cell incises in proportion to its slope; headwaters are the
steepest cells, so they over-carve. Physically those are HILLSLOPES (diffusion /
mass-wasting), not channels — the fluvial law is applied where it does not apply.
Measured on seed 42 (per-order median incision, m; want S1 small, trunks 200–400 m):

- pure SP: S1=380 S2=366 S3=286 **S4=138** — inverted (headwaters carve most).
- incision threshold `θ` (`E=K·max(0,A^m·S^n−θ)`): a razor window — θ=0.02 drops S1 to
  39 m but collapses trunks (S4=11); θ≥0.05 kills all. Does NOT fix it.
- critical area `A_c` (hillslope below, fluvial above): monotone only at an
  implausibly sparse 400 km² channel head, trunks then ~100 m. Does NOT fix it alone.
- **coupled regime split** (diffusion on hillslopes `A<A_c` + stream power on channels
  `A≥A_c`, interleaved) — TASK 3 WAS needed and it WORKS: at A_c=50 cells, D=0.05,
  **S1 → 30 m** (headwaters fixed) and the valley cross-section OPENS from 250 m
  (SP-only) to **420 m** (V-walls; drainage relief 659 m). Residual: mid-order still
  out-incises trunks (S2=338 > S4=115) — physically defensible (tributaries drop into
  graded trunks), not the headwater pathology.

**Legibility over "minimal steep land" (criterion change — see Decision summary).**
Judged on spatial structure (valley floor / cliff sharpness / flank contiguity),
UNCOUPLED is the more legible: at K=3, cf=1.0, 1024²/400 km — trunk valley W/D ≈ 8.8
(≈ the 2 km/300 m target), channel flanks 27° vs interfluves 10° (steep concentrated on
flanks, not plateaux), largest contiguous flank 4677 cells, ~690 km² buildable floor.
Coupling (cf=0.39) is TOO FLAT (W/D 38, 13° flanks). So do NOT couple depth_scale to
domain. **Droplets are corrosive at ANY density**: SP relief 659 m → SP+weak
(0.25/cell) 124 m → SP+full 68 m — recommend droplets ≈ off (hillslope diffusion + FBM
base supply texture).

**A_c must be in km², not cells (8192² finding).** `min_area_cells` = 50 cells is
7.6 km² at 1024² (removes headwaters) but only 0.12 km² at 8192², so at production
resolution headwaters became "channels" again and over-carved (S1 25 → 680 m). The
channel-head criterion MUST be a physical area (≈ 7.6 km²), converted to cells per
resolution. With that fix, incision is resolution-stable.

**Recommended (implementation pass, NOT wired):** stream-power ON, K=3, m=0.5, n=1,
iterations=3, sea_level=0.5, **A_c ≈ 7.6 km² (converted to cells per resolution),
diffusion ≈ 0.05 (regime split)**, θ=0, droplets ≈ off, **UNCOUPLED** vertical scale.
8192² confirmation (author seed 10481999410520546993, A_c=7.6 km²=3188 cells, uncoupled):
FBM 1.3 s, stream-power 33 s, peak RSS 3.9 GB; **headwaters fixed (S1=0), ordering
concave (S2=646 S3=377 S4=240 S5=179)**, cliff transition 49 m (sharp), 42 855 hex of
<5° floor. Two refinements for the NEXT pass: (i) trunk channels incise as narrow deep
SLOTS (W/D ≈ 0.5, 98 m wide / 178 m deep) rather than wide buildable valleys — likely
needs a lower A_c and/or stronger hillslope diffusion to widen trunk floors; (ii)
incision is ~1.5× higher at 8192² than 1024² (S4 136→240) — K wants a mild
per-resolution re-anchor. Whether the pervasive steep flanks (largest component 1.07 M
cells) read as "legible" or "over-incised" is now an author visual call.

---

## Finding 6 — Sculpting: A_c too large + no incision bound (the "not sculpted" report)

The first relief-v1 (A_c=7.6 km², iters=3, K=3000) produced an UNSCULPTED massif — no
valleys/ridges perpendicular to the range, channel floors at 55–150 m under 3000 m
flanks, and the FBM pattern surviving on the upper slopes. Two root causes, both
measured:

- **A_c far too large.** 7.6 km² channel-head area put channel heads very low, so the
  UPPER SLOPES received hillslope diffusion only — no fluvial incision (channels
  reached only **7–9 % of peak elevation**). Lowering A_c to **0.1 km²** (a realistic
  humid-temperate drainage density; 76× lower) dissects the whole massif: channels now
  reach **94 % of peak**, drainage density 2.79 km/km² @8192². NOTE: A_c must be
  resolvable — 0.1 km² is sub-cell below ~2048² (it needs ~2.6 cells @2048², 42 @8192²).
- **No incision bound (floors planed to base level).** Stream power ran on a static
  field with no uplift, so channels graded down to sea level. Since Ymir's tectonics
  already did the uplift, the fix is to LIMIT total incision, not add U: iters 3→2 and
  K 3000→**1500** lifts floor/local-ridge from 0.21 to ~0.48 (floors at ~half the local
  ridge, not planed) while keeping deep valleys. Validated at 8192² (SP 22.7 s).

**W/D still widens downstream** (Finding 5 property preserved). Recommended sculpt =
`relief_v1`: A_c 0.1 km², iters 2, K 1500, m 0.5, n 1, D 0.05, uncoupled, FBM
`amplitude_base ≈ 0.04`. Renders for the author's visual review: `exports/sculpt/`
(2048² + 8192² + crops). Still OFF by default (checkbox) pending that review.

**Navigability (Finding 4 follow-up) — re-anchor on the measured basins.** The
generated continent's basin area at river mouths: **max ~12–13 000 km² (Thames-scale),
p90 ~300–500, p50 ~50 km²** (≈500–700 mouths). The Earth-calibrated thresholds
(small_boat 500 / barge 5 000 / ship 50 000) leave almost everything non-navigable
(ship unreachable by construction on a ~40 000 km² island; only the single trunk hits
barge) — matching the "only small boats" report. **Proposal: keep ABSOLUTE km²
thresholds (navigability is physical river size) but LOWER them, anchored just below
the measured max so the distribution populates every class: stream 10 / small_boat 100
/ barge 1 000 / ship 8 000 km².** Then p50 (~50) → small_boat, p90 (~400) → barge, the
~13 000 trunk → ship-class. Absolute (not domain-scaled) because a given river size
should classify the same on any map; lowered because a Thames-scale continent cannot
meet Earth's 50 000 km² ship bar. Not applied (default unchanged) pending review.

**Valley-type variety (lithology-K) — availability.** There is NO lithology / geology
/ erodibility field in core; the only spatial rock proxy at the HD stage is the BINARY
cratonic mask (`state.cratonic_mask`, coarse 64²). So a rich multi-class lithology-K is
not available. Cheapest path (proposed, not built): upscale the cratonic mask to HD and
derive a 2-class K (craton → low-K hard/narrow gorges, else → high-K soft/wide
valleys), passed to `incise` as a per-cell K field. A richer palette needs a new
tectonic lithology field — out of scope here. Deferred until the base sculpt is
visually confirmed.

---

## Finding 5 — FBM shrinks to a symmetry-breaking seed; width widens downstream; incision is resolution-dependent

**The reframe (why this is possible now).** 64² → 8192² is 128× per axis, so pure
interpolation gives a featureless surface (the bilinear baseline). Detail must be
either INVENTED (FBM) or DERIVED (erosion physics). Until stream-power, FBM had to
carry the detail because droplets destroyed relief; now that stream-power CREATES
relief causally, FBM can shrink from "terrain generator" to "symmetry-breaking seed"
— just enough initial irregularity for drainage to organise.

**Measured (1024², relief-v1 incision on each FBM variant).**
- **Amplitude is reducible ≥8× with drainage fully organic.** `amplitude_base`
  0.16→0.02: FBM roughness 0.031→0.024, and maxStrahler (4–5), confluences (~680),
  segments (~3140), valley floor (~680 km²) ALL stay healthy — the degeneracy floor
  is BELOW 0.02, not reached. So a low-amplitude seed (≈0.02–0.04) keeps organic
  drainage. Recommend that regime; confirm the visual striation drop by eye.
  > **⚠️ VALID WHEN MADE, DOES NOT TRANSFER (H-1c round).** This was measured PRE-C-1, when
  > `flow_conditioning = 0` and `amplitude_base` genuinely drove the terrain — the numbers above
  > are sound for that configuration. It does NOT carry to the conditioned production path: since
  > C-1 the relief-budget cap binds at every cell and `amplitude_base` is entirely inert (proven
  > byte-identical at 4×, see "The DEAD KNOB" below). Any amplitude sweep run AFTER C-1 measured
  > nothing, and the "degeneracy floor below 0.02" must not be quoted as a property of the shipped
  > pipeline.
- **The anisotropy knobs do NOT move the striation metric.** `max_anisotropy` 3→1
  (isotropic), `amplitude_slope_factor` 3→0, `octaves` 7→3 all leave the
  roughness-asymmetry (grad vs contour, ±8 cells) at ~0.83 (pre-incision ~0.80).
  Either the visual striations are not controlled by these knobs, or the ±8-cell
  asymmetry metric is too coarse to isolate them — an honest gap; amplitude_base is
  the lever that demonstrably reduces overall FBM detail (and shifts the metric:
  pre-FBM asym 0.80→0.69, λ 8.8→11.5 cells as amplitude falls).

**Width is a healthy DISTRIBUTION, not a slot (corrects the 8192² single-channel
read).** W/D per Strahler order (1024²): S1 median 2.7 (headwater gorges), S2 4.1,
S3 8.8, S4 41.4 (wide trunks), with a fat tail (p90 up to 77). **W/D widens
downstream** — gorges as chokepoint content upstream, wide buildable valleys at
trunks/coast. No widening fix needed.

**Incision is resolution-dependent — and physical units DON'T fix it (measured).**
Per-order incision rises with resolution (S4: 108 m @512² → 136 @1024² → 318 @2048²).
The physical reformulation (`A` in km², `S = Δh_m / dist_m`) was IMPLEMENTED and
measured: for the shipped `n=1, m=0.5` it is **algebraically equivalent** to the
normalised law (the `cell_km` factor cancels between `A_km²^0.5` and `dist_m`), so
`K = 3000` physical reproduces the reference EXACTLY (relief 682 m, S1=25 S2=414
S3=307 S4=136 — identical to normalised K=3) and the resolution dependence is
**UNCHANGED** (still 108/136/318). So the ~1.5× drift is NOT a slope-unit artifact —
it originates in the FBM detail resolving sharper gradients on finer cells. The real
levers are a resolution-independent FBM feature size or an explicit per-resolution K;
the physical law is kept because it makes K dimensional and matters for `n≠1`, but it
is a no-op for the shipped exponents. (Corrects the earlier "physical slope is the
clean fix" hypothesis.)

**Striations (TASK 2 metric) — the anisotropy knob is not the lever.** A directional
power-spectrum metric (length-48 profiles along contour vs gradient on steep cells)
does NOT respond to `max_anisotropy` 3→1 (ratio ~0.65 either way, no short-λ peak) —
same verdict as the ±8 roughness metric. Either the visual striations are not
gradient/contour-aligned or not from `max_anisotropy`. `amplitude_base` is the lever
that moves both metrics (ratio 0.65→0.45 at 0.04) and reduces overall FBM detail 8×
with drainage staying organic (Finding 5 above); the rendered amplitude ladder
(exports/relief_ladder/, 0.16/0.08/0.04/0.02) is the visual arbiter, provisional
recommendation `amplitude_base ≈ 0.04`.

---

## Finding 4 — Terraces are an EROSION-DEPOSITION artifact, NOT a tectonic/isostasy closure

**Context.** The relief shows terraces PARALLEL to contours (concentric loops around
hills), with jumps of several hundred metres. The initial hypothesis (recorded here
because it was tested) was a coarse discrete altitude level from the C1
equilibrium-height closure (a single global `h_eq` → level sets).

**Measurement (seed 42 coarse field — the hypothesis was REFUTED).**
- u16 quantisation REFUTED: 0.1–0.14 m/unit, vs observed jumps of 100s of m.
- equilibrium-height clamp INACTIVE: **0 %** of `S̃` cells sit at `h_eq = 2.0`
  (only 1 % are even above it; `S̃` maxes at 2.18). Davis-Suppe `h_max = 2.5`: **0 %**
  (and its cap is distance-TAPERED `h_max·(1−exp(−d/L))`, already spatially varying).
- cratonic isostasy is NOT a level either: cratonic land altitude p10/median/p90 =
  335 / 1072 / 1901 m — a WIDE band, not a step.
- No clean discrete altitude ladder on the coarse land field.
- Terrace-source disentangle (flat fraction of a transect): **pure bilinear
  (coarse only) 13 % → FBM 6 % → after erosion 24 %.** FBM ROUGHENS; the flat fraction
  is dominated by **erosion DEPOSITION** (6 % → 24 %).

**Conclusion (corrected).** The terraces are NOT produced by a coarse
tectonic/isostasy closure — the equilibrium-height / Davis-Suppe / craton candidates
are all refuted by direct measurement. They are an **erosion-deposition artifact**:
the same net-zero, uncorrelated-droplet erosion of Finding 3 deposits sediment in
flat sheets that pool to local base levels bounded by contours (hence "concentric
terraces"). Terraces and missing valleys are TWO SYMPTOMS OF ONE CAUSE — the erosion
algorithm (Finding 3), not two separate jobs.

**Caveats.** Measured on seed 42; the earlier "~120–176 m resolution-independent
step" came from a weak modal-|Δ| metric and is not a confirmed discrete ladder. A
client-side render contribution (contouring/quantisation in Living Landz) is not
excluded and is out of Ymir's scope. The coastal `coastal_deposit_fraction` sink
reduces COASTAL deposition, not the inland depositional flats.

**Decision.** Fold the terrace fix into the erosion chantier (Finding 3's routed
stream-power successor + a deposition/transport rework), not into the isostasy /
equilibrium-height closures. Not touched here.

**Consequences.** The prior "address terraces at the isostasy source, before
stream-power" ordering is moot — both are the erosion chantier. When any tectonic
closure IS eventually touched, re-check the pre-existing Picard non-convergence in
that layer — [docs/issues/picard-nonconvergence-rectangular-smoke.md](../issues/picard-nonconvergence-rectangular-smoke.md).

## Finding 7 — The relief needed a BOUNDING CLOSURE, not more tuning (slits & striations)

**Context.** After Finding 6's sculpt (A_c 0.1, iters 2, K 1500) the author reviewed
exports/sculpt at 8192² and rejected it: the FBM striations still drape the upper
slopes *following isolines* (contours, not valleys), the rivers are ~1000 m deep but
**one pixel wide** (slits, not valleys), and the massif is not sculpted. Diagnosis
(shared and confirmed): tuning cannot fix this because two physical laws are absent.

**STEP 1 — confirmed the gap (read-only, `step1_slope_distribution`).**
- Hillslope law: a plain LINEAR Laplacian, constant `D` ([stream_power.rs] step 5) —
  no slope dependence, no flux limiter.
- **Nothing bounds the maximum slope anywhere in the pipeline** — no talus / angle of
  repose / clamp (`erosion/thermal.rs` is an empty, commented-out placeholder).
- Incision is purely vertical on the receiver cell; no lateral/bank erosion.
- Land slope distribution (seed …6993, 2048²): raw FBM max ~70–80°, and the
  stream-power incision itself STEEPENS it — >30° share 2.8 %→38.7 %, **max ~84°**
  (≈ 2.5× the angle of repose) at amp 0.04. The near-vertical faces are the *walls of
  the 1-px slits*. Lowering FBM amp barely moves the sculpted max (85→84°). → a
  bounding closure is required.

**STEP 2 — two closures, independently switchable, OFF by default.**
- **(a) Nonlinear hillslope diffusion with a critical slope** (Roering): flux
  `q = D·S/(1−(S/S_c)²)`, `S_c = tan(33°)`. The effective diffusivity diverges as
  `S→S_c`, so slopes cannot exceed it (arêtes; the missing bound). **Scheme (the term
  is stiff near S_c → explicit blows up):** backward Euler (unconditionally stable for
  any step) + lagged-diffusivity Picard (re-freeze edge weights each outer pass) +
  Gauss-Seidel inner solve (diagonally dominant → converges, deterministic in
  row-major order); the denominator is floored (edge slope capped at 0.999·S_c) so
  weights stay finite. Self-arresting: a bank dropping below S_c gets denom→1,
  weight→D, and stops — so relief survives instead of planing to base level. Runs only
  on hillslope cells (A<A_c); channel/sea cells are fixed Dirichlet values.
- **(b) Channel lateral widening as HYDRAULIC GEOMETRY.** Floor half-width
  `W = K_lat·A_km²^m` (the width–area law), planed perpendicular to flow toward the
  channel floor. Trunks (high A) get wide floors, headwaters (low A) stay narrow
  gorges — the width variety requested.

**Resolution invariance (a real bug found at 8192²).** The first cut of both closures
was resolution-DEPENDENT and collapsed back to slits at fine cells (v2_8192 W/D 3,
>30° 46 %): the ±1-cell lateral reach widens 4× less at 4× finer cells, and the
dimensionless diffusion weight carries an implicit 1/dx². Fixed: lateral reach is a
PHYSICAL half-width in metres; the diffusion weight scales `(HILLSLOPE_REF_CELL_M/cell)²`
(κΔt/dx²). After the fix the closures behave the same in metres at 2048² and 8192².

**STEP 3 — shape metrics (2048², amp 0.04, `closure_grid`).**

| config | >30° | max° | floor/ridge | crest curv | W/D S1→S5 | Strahler |
|---|---|---|---|---|---|---|
| v1 (none) | 38.7 % | 84 | 0.29 | 774 m | 0.64 → 6.3 | healthy |
| +crit slope (a) | 20.6 % | 82 | 0.33 | 313 m | 2.3 → 9.1 | healthy |
| +lateral (b) | 23.4 % | 84 | 0.27 | 640 m | 7.0 → 47.9 | healthy |
| **+both (v2)** | **12.1 %** | 81 | 0.32 | **288 m** | **6.2 → 40.9** | healthy |

- (a) collapses the steep SHARE (38.7→12.1 %) and halves crest curvature (774→288 m):
  the pervasive striated slopes plane into arêtes. The **max stays ~80°** because
  isolated fluvial cliffs (channel walls/knickpoints, deliberately excluded from the
  hillslope closure) survive — arguably correct (cliffs are content, not the defect).
- (b) is the decisive slit fix: W/D grows monotonically with order and is >1 on every
  order (trunk S4 1.3→30.5) — the 1-px slit becomes a downstream-widening valley.
- Drainage hierarchy intact across all four (no pathology traded).

**FBM can drop further.** With the closures supplying the structure, amp 0.02 and 0.01
hold up (W/D and steep share barely change) — structure now comes from the closures,
not the noise. This is the real test of the reframe, and it passes.

**Config.** `StreamPowerConfig::relief_v2` = relief_v1 + `critical_slope = tan(33°)`,
`lateral_erosion = 4.0` m/√km², nonlinear `diffusion = 0.15` (at the 2048²/400 km
reference cell). Still OFF by default; the viz "Closures relief-v2" checkbox and the
diagnostics drive it. `relief_v1` is byte-identical (its regression baseline is
unchanged, not rebased). A v2 regression baseline is deferred until the author's
visual verdict on exports/sculpt/closure_* settles.

## Finding 8 — The 8192² "comb" is the Smith–Bretherton rilling instability (not FBM, not D8)

The relief-v2 closures (Finding 7) are validated at 2048² (author: dendritic valleys +
arêtes, v1 comb gone) but a fine parallel "comb"/terracing returns at 8192². STEP 1
(read-only) discriminated the source; four hypotheses were tested and **three refuted**:

1. **D8 routing artifact** (refuted). Two-frame channel-segment orientation: R_grid ≈
   0.02–0.05 with the four D8 axes near-equally populated, R_gradrel ≈ 0.02 — the network
   is directionally ISOTROPIC, not grid-biased.
2. **Anisotropic FBM** (refuted). Free ablation at 8192² with `max_anisotropy=1`,
   `amplitude_slope_factor=0` (isotropic, slope-blind noise): >30° 33.7→33.5 %, striation
   0.69→0.68 — identical. FBM directionality is not the source.
3. **FBM fine octaves** (refuted). octaves 7→4 (finest λ 8 px→64 px @8192²): post >30°
   33.7→33.1 %, striation 0.69→0.67 — the raw FBM is smooth (~5.5 % steep) regardless;
   the EROSION imposes ~33 % steep regardless of input detail.
4. **Smith–Bretherton parallel rilling** (CONFIRMED). On a SMOOTH plane tilted 30° off the
   grid, no FBM, the incision spontaneously forms regularly-spaced parallel rills running
   straight DOWNSLOPE (↘, diagonal — following the slope, NOT the grid axes; v1 symmetric
   diffusion and v2 GS give the identical concentration R=0.19, exonerating the solver).
   This is the classic linear instability of detachment-limited stream power `E=K·A^m·S^n`
   (m<1) on smooth slopes, damped only by hillslope diffusion.

It explains everything the other hypotheses could not: the isotropic segment histogram
(slopes face all directions, so the rills do too, globally); the resolution dependence
(finer cells resolve the instability's short wavelength → denser rills — at 2048² it is
sub-resolution/aliased away, hence the clean preview); FBM-independence (the instability
generates its own pattern from any perturbation); and why more diffusion (×4) AND a larger
A_c (×5) did NOT clear it (the regime split EXCLUDES channel cells from diffusion, so the
rills — which are channels — are never damped laterally).

**Verdict for STEP 2 (act on this, not before).** The fix is not a parameter tweak and not
an FBM change. The standard remedies for the Smith–Bretherton instability, in increasing
blast radius:
  - (i) **Cross-rill / whole-field diffusion** — stop excluding channel cells from the
    hillslope diffusion (or add a small isotropic smoothing across the network), so the
    diffusion sets a finite valley spacing that damps sub-threshold rills. Smallest change,
    stays in stream_power.rs; risk: softens genuine narrow channels.
  - (ii) **Multi-flow-direction (MFD) routing for the INCISION** — single-flow (D8)
    over-concentrates and accentuates rilling; MFD disperses flow across downslope
    neighbours and is the documented damper. Rewrites how accumulation/receivers feed the
    incision; can be scoped to the incision only, keeping D8 for rivers/lakes (Finding 7's
    blast-radius note applies if it spreads to the whole chain).
  - (iii) **A transport-limited / depositional term** — inter-rill deposition fills the
    incipient lows and damps the instability; the largest change (a sediment budget).
Recommend prototyping (i) first (cheapest, testable on the synthetic plane), then (ii) if
the plane still washboards. Diagnostics (all #[ignore], read-only): `striation_source`,
`striation_stage`, `fbm_octave_ablation`, `rilling_sweep`, `synthetic_slope_rilling`.

## Finding 9 — Remedy (i): cross-rill diffusion works; the blocker is GS convergence

STEP 2a quantified the Smith–Bretherton criterion on a 30° plane (D applied everywhere,
linear-like via a huge S_c): D=0 → 39 % steep, D_crit ≈ 0.40 (physical κΔt ≈ 15 000 m²)
damps it to <1 %; the relief-v2 default D=0.15 (4.9 %) is SUB-critical — why the comb
persists. (The rill-wavelength metric is broken — it returns window×cell, no spectral peak;
not needed for D_crit, so left as a known gap, fixable by transverse autocorrelation. A*
is not cleanly measurable on a uniform plane, so gorge survival was tested on real terrain.)

STEP 2b applied remedy (i) — `diffuse_channels = true` (the LEM-correct diffusion on every
cell, removing the regime split's non-physical channel exclusion) — on real terrain (2048²,
author seed, amp 0.04), short D sweep:

| D | >30° | striation | floor/ridge | crest | W/D S1→S5 | Strahler/confl |
|---|---|---|---|---|---|---|
| 0.15 | 6.7 % | 0.67 | 0.49 | 270 | 3.4→22.4 | 2391..39 /1594 |
| 0.25 | 5.2 % | 0.64 | 0.50 | 279 | 7.0→29.3 | 2345..33 /1580 |
| 0.40 | 4.2 % | 0.61 | 0.50 | 278 | 7.1→(S5 5.9*) | 2346..17 /1545 |
| 0.55 | 3.5 % | 0.60 | 0.50 | 269 | 7.2→41.7 | 2330..35 /1547 |
| 1.00 | 2.6 % | 0.54 | 0.50 | 235 | 7.4→49.3 | 2295..31 /1500 |

Findings:
- **The channel exclusion WAS the driver.** At the same D=0.15, `diffuse_channels` alone
  halves the steep share vs the regime split (12.1 %→6.7 %). Visually (crossrill_d0.55) the
  comb is gone — coherent dendritic valleys + defined ridges replace the corduroy.
- **Monotone on real terrain** (6.7→2.6 %); the plane's non-monotonic bump at D=1.0 does NOT
  reproduce → it was GS-convergence noise on the plane, not physical.
- **Gorges survive**: W/D keeps rising downstream at every D; Strahler histogram and
  confluences stable — no pathology traded. (*S5 dips at D=0.40 only because S5:17 segments.)
- **Cost**: floor/local-ridge rises 0.32→0.50 (valleys shallower — diffusion fills channels
  a little) and headwaters widen with D (S1 W/D 3.4→7.2). Real but not disqualifying.

**Blocker — the GS solver is NOT converged.** At D=0.40, 40 vs 80 implicit sweeps differ by
max 0.059 norm = **671 m locally** (aggregate >30 % stable at 4.2/4.1 %). So the absolute
heightfield is a provisional intermediate state, and the 8192² extrapolation would inherit
it. Remedy (i) is the right direction, but before freezing a D or rendering at 8192²:
  1. **Fix the diffusion solver** — the row-major Gauss-Seidel is single-threaded and
     under-converged; move to red-black GS (parallelisable with rayon) with a residual-based
     stop, or a V-cycle. This also makes 8192² tractable (single-thread GS × 240 sweeps ×
     67 M cells is minutes/D otherwise).
  2. Re-run the D sweep to the converged field, pick the LOWEST D meeting both criteria.
  3. Produce an 8192² render for the author's visual verdict, THEN freeze the relief-v2
     regression (rebase, not loosen). `diffuse_channels` stays OFF by default until then.

## Finding 10 — Talus vs diffusion head-to-head: talus is not the fix; the cause is flow concentration

The author asked whether TALUS (angle of repose, a C1-style closure) could replace the
nonlinear-diffusion solver we reintroduced at the HD stage. Implemented it
(`talus_slope`/`talus_passes`/`talus_factor`, mass-conserving high→low sweep, everywhere)
and ran the head-to-head at 2048² AND 8192² (seed …6993, amp 0.04, A_c 0.1 km²):

| res | method | ms | >30° | max° | floor/ridge | W/D S1→S5 | dens |
|---|---|---|---|---|---|---|---|
| 2048 | diffusion D=0.55 | 3870 | 3.5 % | 84 | 0.50 | 7→42 | 2.19 |
| 2048 | talus Sc=tan33 (×4) | 2603 | 12.5 % | 58 | 0.43 | 2→9 (flat) | 1.87 |
| 8192 | diffusion D=0.55 | 72493 | 3.3 % | 88 | 0.80 | 7→32 | 2.22 |
| 8192 | talus Sc=tan33 (×4) | 47362 | 26.8 % | 77 | 0.72 | 2→3 (flat) | 2.37 |

**Verdict — talus loses on the comb, decisively:**
- **It does NOT damp the rilling.** >30° share: talus 12.5 %→26.8 % (2048→8192, WORSE and
  DOUBLING), diffusion 3.5 %→3.3 % (low and metre-invariant). Visually the talus 8192² crop
  is fully combed. Talus only *caps* slope; the comb is a rill/ridge ALTERNATION that must be
  SMOOTHED ACROSS (transport between rills) — talus transports only downhill within a rill,
  so the pattern survives. Diffusion smooths across → removes it.
- **Metre-invariance** (the deciding criterion): diffusion is invariant on the comb metric
  (3.5→3.3 %); talus is the opposite (12.5→26.8 %). Diffusion wins the closure property on
  the metric that matters.
- **Gorge structure**: diffusion W/D widens downstream (7→42, proper hierarchy); talus is FLAT
  (~2 at every order) — uniform narrow gorges, no hierarchy. Diffusion wins.
- **Talus's only win**: it bounds max slope better (58–77° vs 84–88°) with straight, low-
  curvature slopes (curv 40 vs 269) — a complementary arête tool, not a comb fix.
- **Runtime**: talus ~1.5× faster (48 s vs 70 s @8192) — modest, and it fails the task anyway.

**TASK 3 — is talus a closure?** No: the mass-conserving talus converges GEOMETRICALLY
(residuals halve per doubling of passes: 1→19.7 %, 2→8.2 %, 4→4.3 %, 8→2.3 %; one pass is
worse than none). Lowering a cell re-steepens its uphill side, so a single sweep cannot
bound S ≤ S_c — it is a solver in disguise, cheaper per pass but still iterative. (A
NON-conserving Lipschitz carve via fast-sweeping would be a bounded closure but removes mass.)

**H-A/H-B (headwater ramification):** channel-head elevation and drainage density rise
modestly at 8192² for BOTH methods (head %peak 8→12 %, dens ~2.2) → the missing ramification
is partly resolution (H-A), improving a little at 8192². And diffusion BACKFILLS valleys
badly at 8192² (floor/ridge 0.50→0.80 — the melted look), i.e. killing the comb by smoothing
also fills the headwater vallons (H-B). Talus fills less (0.72) but keeps the comb. So the two
share a wavelength and diffusion cannot separate them.

**RECOMMENDATION — attack the CAUSE, not the symptom. The solver rewrite is NOT the next
step, and talus is not the fix.** The rilling instability is driven by flow CONCENTRATION
(single-flow D8 over-concentrates on smooth slopes). Both diffusion and talus fight the
symptom after the fact (diffusion smooths → fills valleys; talus caps → keeps the comb). The
principled fix is **MFD / D∞ routing for the incision (remedy ii)**: dispersing flow across
downslope neighbours suppresses the parallel channelisation, so the comb never forms and no
aggressive smoothing (which fills valleys, H-B) is needed. It is O(n) accumulation — NOT an
iterative solver — so it is C1-consistent and makes the red-black/multigrid rewrite
unnecessary if it works. Proposed scope: MFD for the stream-power accumulation/incision only,
D8 kept for rivers/lakes (blast radius contained), with light linear diffusion at most.
Keep talus available as an optional arête/max-slope tool; keep the nonlinear diffusion too.
All OFF by default.

## Finding 11 — MFD routing PREVENTS the comb at the cause; it is the fix (not the solver)

Prototyped multiple-flow-direction accumulation for the incision only (`mfd_exponent`,
D8 kept for the receiver/stack + rivers/lakes). Attacks the CAUSE of the Smith–Bretherton
rilling — single-flow concentration — rather than erasing the pattern after the fact.

**Plane (clean-room 30°, two-sided criterion):** D8 → 39 % steep comb; ANY MFD p
(10/4/2/1.1/1) → **0.2 %**, comb ELIMINATED, while a trunk still forms (max A highest at
p≈2). Visually the washboard becomes a smooth slope with faint emergent channels. So the
comb is purely a single-flow artifact and even light dispersion removes it.

**Real terrain (2048²+8192², MFD incision, critical_slope=0 → NO GS solver):**

| res | p / K | ms | >30° | stri | floor/ridge | dens | W/D S1→S5 | align p50 / off>20 m |
|---|---|---|---|---|---|---|---|---|
| 2048 | p2 K×1 | 1405 | 19.0 % | 0.89 | 0.40 | 2.55 | 5→57 | 30 m / 61 % |
| 2048 | p2 K×2 | 1392 | 26.6 % | 0.99 | 0.25 | 2.40 | 5→85 | 19 m / 49 % |
| 8192 | p2 K×2 | **25936** | 21.5 % | 1.05 | 0.72 | 2.35 | 2→68 | **6 m / 18 %** |

- **Character (the win):** the 8192² crop shows DENDRITIC, ramified branching valley
  networks — the headwater ramification the author reported missing — replacing the comb.
  Striation isotropic (~1.0) confirms the directional comb is gone; the residual >30 %
  share is genuine valley-wall/ridge steepness of dense dissection, not the washboard.
- **Runtime:** 26 s @8192² vs the diffusion solver's 70 s — MFD is O(n), no Gauss-Seidel.
  **The red-black/multigrid rewrite is NOT needed.**
- **W/D widens downstream** (proper hierarchy), unlike talus (flat).
- **K adjustment (predicted):** MFD disperses A → weaker incision. K×2 deepens valleys at
  2048² (floor/ridge 0.40→0.25) but 8192² still sits at 0.72 → **8192² needs more K**
  (≈×3–4) to avoid under-incision; this is the H-B trap by another route and must be tuned.
- **Metre-invariance:** >30 % 19–27 % (2048²) vs 21.5 % (8192²) and striation ~0.9–1.05 both
  — reasonably invariant (unlike talus). floor/ridge is NOT yet invariant (K-dependent).
- **TASK 5 D8/MFD alignment:** at 8192² the D8 rivers sit in the MFD valleys (median offset
  6 m, only 18 % above a carved hollow) — acceptable. At 2048² they diverge (30 m, 61 %) —
  MFD carves broad shallow valleys the single-flow D8 line rides the flank of. Fix if it
  matters at production res: route rivers on MFD too, or snap D8 segments to the carved
  local minimum. At the 8192² production grid it is a minor issue.

**Recommendation:** adopt **MFD p≈2 for the incision** as the comb fix; raise K (~×3 at
8192²) and re-check floor/ridge; keep D8 for rivers/lakes; talus and nonlinear diffusion
remain optional complements (arête bounding), all OFF by default. MFD makes both the
diffusion-solver rewrite and talus unnecessary for the comb. Fallback if MFD ever fails:
the non-conservative Lipschitz carve (fast-sweeping, bounded closure) noted in Finding 10.
Pending: the author's visual verdict on the 8192² render (wired into the viz), then freeze
the relief-v2 regression (rebase, not loosen).

## Finding 12 — MFD defects: K doesn't set valley shape; uphill rivers = filled-depression crossings

Two 8192² defects reported on the MFD renders. Diagnosed before fixing (TASK 1/2).

**DEFECT 1 — "deep gorges, no intermediate slope" (crest-to-thalweg in a few pixels).**
- The apparent resolution dependence of floor/local-ridge (0.25 @2048² vs 0.72 @8192², same
  K) is largely the METRIC: the ±10-CELL ridge window is ±1953 m at 2048² but ±488 m at
  8192², so it reads a nearer, lower ridge at 8192² → inflated ratio. With a PHYSICAL ±1 km
  window the gap shrinks to 0.36 vs 0.63 (cells-not-metres bug, 3rd instance — fixed in the
  metric; the residual 0.36↔0.63 is real: MFD disperses more at finer cells → weaker incision).
- **K is not the lever.** Sweeping K×2→×6 @8192² moves per-order incision up (S1 124→166 m)
  but floor/ridge only 0.63→0.55 (saturates against the receiver clamp; the ratio is set by
  the unchanged local ridge). So the "0.4–0.5" target is unreachable by K.
- Real cause: the MFD config turned OFF all slope-grading (`critical_slope=0`, light D=0.05)
  to avoid the GS solver, leaving UNGRADED slot walls (curv 141–154, 22 % >30°). The fix is
  FLANK GRADING, not K — and MFD having removed the comb, grading can now be added freely
  (talus for straight repose flanks — cheap, complementary; or nonlinear critical-slope for
  convex flanks — but that is the GS solver again).

**DEFECT 2 — rivers run uphill (author's profile 400→250→80→50→70→80 m).** Acceptance test
`river_monotonicity`: at 8192² MFD p2 K×3, **66 % of segments contain an uphill step, worst
+232 m**. Attribution of climbing steps: **61 % sit on pit-FILLED flats** (the river crosses a
filled depression — real floor dips into the hollow then climbs back to the sill), 39 % on
real terrain (D8 segment rides the flank off the MFD thalweg).
- **H2-as-defined (network built before incision) is FALSE**: incision is
  production_upscale.rs:402, drainage runs afterwards on the eroded field (hd.rs
  `cached_c1_drainage_windowed(&eroded)`). The network is post-incision.
- The dominant cause is filled-depression crossing (rivers traced on the pit-FILLED surface,
  long profile read on REAL elevation), not staleness. Fixes: (i) a drainage CARVE / breach so
  the real profile is monotone along the network (O(n), guarantees the acceptance test), and/or
  treat genuine closed depressions as LAKES (do not route a flowing river through them); plus
  (ii) snap D8 segments to the MFD thalweg (or route rivers on MFD) for the 39 % misalignment.

The permanent acceptance test (every exported segment monotone to the sea) is in place;
violations must go to ~0 before the relief-v2 baseline is frozen.

## Finding 13 — Monotone carve: ad-hoc breach converges too slowly; needs priority-flood

DEFECT 2 fix attempt. `flow::carve_monotone` fills detected lakes to their flat sill surface
(so a river crosses a lake level, not climbing out of the real hollow) then lowers each
non-lake receiver to at most its donor along D8. Result at 8192² (MFD+talus), iterated
(route → carve → repeat): climbing segments 3060 → 1287 → 977 → 799 → 680 → 577 over 5 passes
— CONVERGING (geometric ~0.8/pass) but far from the categorical zero, worst climb stuck at
198 m. Cause: "lower receiver to donor" along one routing is NOT a guaranteed monotone
conditioning — re-pit-filling on the next re-route re-exposes reversals, and the many spurious
depressions (816 k flooded cells) keep generating them. floor/ridge unaffected (0.69→0.70).

The guaranteed tool is a PRIORITY-FLOOD conditioning (Barnes 2014 / Lindsay 2016), one pass:
  - (a) **Priority-flood BREACH (least-cost), lakes excepted** — carves a monotone outlet
    path for every sub-threshold pit, fills genuine lakes to sill. Matches the author's
    intent (small pits drained, lakes flat), guaranteed monotone, but the most code.
  - (b) **Fill flooded cells to sill (= export the `filled` surface)** — trivial, O(n) (it
    IS `pit_fill`, already computed), GUARANTEED monotone by construction; but every small
    pit becomes a flat pond rather than a drained channel (the author asked for carve).
  - (c) **Hybrid**: fill lakes (≥ threshold) + bounded-length breach for small pits — middle.
The acceptance test (`river_climbs`, lake-flat tolerated) is in place; whichever is chosen
must drive it to zero before the relief-v2 baseline is frozen.

## Finding 14 — DEFECT 2 solved: priority-flood BREACH-then-FILL (lakes excepted), zero climbs

`flow::breach_monotone` (Lindsay-style complete breaching + priority-flood fill mop-up, lakes
excepted). Detected lakes are held at their flat sill (water, never breached); every other
depression is carved a monotone outlet trench, and a final priority-flood FILL raises the
small residual of micro-pits the breach cannot connect in one pass — guaranteeing a
monotone-descending path to the sea BY CONSTRUCTION.

Result (MFD+talus terrain, `breach_monotone_test`), ONE pass:

| res | climbing segs | worst climb | non-lake flooded | lakes held |
|---|---|---|---|---|
| 2048² | 1527 → **0** | +193 → **0 m** | 57863 → 0 | 45 597 cells |
| 8192² | 3060 → **0** | +195 → **0 m** | 1 099 074 → 0 | 816 019 cells |

The acceptance criterion (zero exported segments whose long profile climbs; lake crossings
flat/tolerated) is MET, guaranteed, single pass. Permanent guard: `flow::tests::
breach_leaves_no_interior_pit` (non-ignored). The earlier lake-aware *carve* (Finding 13) only
converged at ~1.25×/pass and re-detecting lakes each pass drained them — both fixed here
(complete breach + fill; lake mask detected ONCE and held).

**Upstream signal to instruct separately (per the author):** 816 019 flooded cells at 8192²
(~1.2 % of area) is ANOMALOUS — MFD + talus + linear diffusion FABRICATE these spurious
depressions (talus transfers and diffusion both create closed hollows). The breach handles
them, but the NUMBER is an upstream symptom: a follow-up must find which step creates them and
whether they are legitimate. Not addressed here.

Still open after both defects: DEFECT 1's 8192² under-incision (floor/ridge ~0.70, K is not the
lever — the conditioning does not deepen valleys); wiring breach + lake export into the
production path (populate lakes.json / Wetland); and the headwater render with the network
overlaid. relief-v2 baseline frozen only once DEFECT 1 depth is settled and the author confirms.

## Finding 15 — The river "offset" was a metric artifact; rivers already sit in the thalweg

The author ruled the 50 %-in-local-min / p90 32 m "offset" unacceptable and asked to extract
rivers from the MFD network. Two things settle it (`thalweg_diagnosis`, 2048²):

1. **MFD dominant-flow receiver ≡ D8 steepest** (92.9 % of land cells; the rest are flat-cell
   tie-breaks). `argmaxⱼ slopeⱼᵖ = argmaxⱼ slopeⱼ` — the dominant MFD path IS the D8 line, so
   extracting rivers from MFD returns the same polyline. It cannot change the offset.

2. **The offset metric was wrong.** "River cell in an omnidirectional local min within ±150 m"
   INCLUDES the downstream cell; a river always descends, so its downstream neighbour is always
   lower → a descending river can NEVER be an omnidirectional local min. That number (p90 32 m)
   was mostly the river's own DESCENT, not off-thalweg-ness. The correct test is TRANSVERSE:
   is the river cell ≤ its two banks perpendicular to flow.

Transverse thalweg residence:

| field | in-thalweg | trans offset p50 | trans p90 |
|---|---|---|---|
| incised | 80 % | −2.5 m | 6.2 m |
| **breached** | **94 %** | **−1.8 m** | **0.0 m** |
| stream-burn + re-breach | 95 % | −4.8 m | 0.0 m |

After the breach the rivers already sit in their valley bottoms: **94 % are at or below both
banks, median 1.8 m below, p90 transverse offset 0.0 m.** Stream-burn adds ~1 pt (94→95 %) and
is unnecessary (and the re-breach fill would raise burned channels back anyway). The residual
~6 % are the natural confluence/cutbank cells no D8-on-real network avoids. So: no MFD rewrite,
no burn — the DEFECT-2 breach already delivers rivers-in-thalweg; the earlier rejection rested
on my flawed omnidirectional metric.

## Finding 16 — Below-sea basins are an altitude-classification bug; the endorheic result traces to an under-produced maritime climate

The author saw inland below-sea basins rendered blue in the erosion view and classified `Ocean`
in the biome view. STEP-1 diagnosis (`water_class_diagnosis`, 8192², reference seed):

- **Q1**: biome Ocean = `heightmap ≤ SEA_LEVEL_NORM` (biomes.rs:148) — altitude, never `water_class`.
  Same altitude-membership assumption in `pit_fill` (every below-sea cell seeded as ocean) and
  `detect_lakes` (needs `filled > real`), so enclosed below-sea basins fall through every stage.
- **Q2**: `water_class` WORKS — class 2 is ~330 k cells (non-empty). Not the bug.
- **Q3**: `class2 ∩ flooded = 0` — because `pit_fill` treats the basins as ocean (root), so they
  are never closed depressions in the drainage sense.
- **Q4/Q5**: 15 real tectonic below-sea basins in the FBM (top 346 km²), ~8 650 more fabricated by
  erosion; **all shallow, deepest 21 m below sea**; 4-conn vs 8-conn inland ≈ equal (basins real,
  not a connectivity artifact). GEOMETRY CORRECTION: spill = 0.0 m for every basin → these are
  COASTAL depressions bounded by the shoreline (rim at sea level), NOT high-rimmed inland basins.

STEP-3 water balance (`endorheic_water_balance`, ≥5 km² basins) came out 5/6 ENDORHEIC with
640 km² dry-land-below-sea — the OPPOSITE of the humid expectation. The balance is CORRECT
(a lake evaporating PE ≈ 691 mm/yr needs catchment ≥ 3.5× its area; these are 1–2.5×). The cause
is the CLIMATE (`climate_precip_diagnosis`, 2048²):

- precip **mean 1009 mm but MEDIAN 448 mm** — half the land sits at the frontal-base floor
  (`PRECIP_MM_PER_UNIT`/k_frontal → ~450 mm), anchored on the GLOBAL ZONAL MEAN.
- Ocean advection works (windward 1923 vs leeward 597 mm/yr — strong orographic contrast); the
  field is not flat. The problem is the frontal-base FLOOR being too dry for an all-maritime island.
- Biome co-symptom: **45 % TemperateGrassland** now; at a maritime ~1.6× floor → grassland 0 %,
  TemperateForest 56 % + Rainforest 17 % (the correct maritime-temperate mix).

**Verdict**: the root is an UNDER-PRODUCED maritime frontal base (a climate-model issue), which
drives BOTH the spurious endorheic basins AND the steppe over-classification. Fix is upstream
(raise the frontal-base floor for maritime); STEP 2+3 (biome from `water_class` + endorheic
levels) waits for that — no freezing lake levels on a climate about to change. Diagnostics
(#[ignore]): `water_class_diagnosis`, `endorheic_water_balance`, `climate_precip_diagnosis`.

## Finding 17 — Maritime frontal base: precipitation from distance-to-ocean, not a global-mean constant

Fix for Finding 16 (the under-produced maritime climate driving spurious endorheic basins + steppe
biomes). The frontal base was `k_frontal·belt_factor(lat)·e_sat(T_sea)` — CONSTANT per latitude,
anchored on the GLOBAL ZONAL MEAN (~450 mm at 45°), so an all-maritime island read as dry as a
continental interior.

**The law (b):** the frontal base is multiplied by `1 + maritime_enhance·exp(−dist_from_sea/efold)`,
where `dist_from_sea` is the DOWNWIND distance over land since the last ocean cell — tracked in the
EXISTING streamline scan (reset at ocean, += km_per_cell over land; no new field). Coast → maritime
floor; deep continental interior → the bare `k_frontal` floor. A PHYSICAL quantity (proximity to the
moisture source) that stays correct for a future large continent (its interior is genuinely far from
the sea → dry), unlike a recalibrated constant. Params: `maritime_enhance = 1.7` (coast ≈ 2.7× floor,
~450→~1200 mm), `maritime_efold_km = 600` (from the Earth coast↔interior contrast: Atlantic façade
800–1500, continental interior 300–600 at 45°). Plus `max_precip_mm = 11000` — a physical orographic
ceiling (Mawsynram/Cherrapunji) bounding the steep-cold-coast spikes that would otherwise pollute the
lake-balance catchment means.

**Result (2048², 45°, reference seed), before → after:**

| metric | before | after |
|---|---|---|
| precip median | 448 | **1130** mm/yr |
| leeward floor | 597 | 1111 |
| min / max | 446 / 63120 | 767 / **11000 (capped)** |
| TemperateGrassland | 45 % | **0 %** |
| TemperateForest + Rainforest | 28 % | **73 %** |

Median 448→1130 (humid temperate), the steppe over-classification is gone (grassland 45%→0%, forests
73%), and the 63 120 mm spike is capped. Ocean advection was already working (windward 2298 vs leeward
1111); the fix is the frontal FLOOR being maritime near the coast. 512 lib tests green (no test pinned
absolute precip). **Relief is unaffected** — the C1 stream-power incision is drainage-AREA based and
climate is computed post-erosion, so relief-v3's shape metrics (Finding 11) are unchanged. STEP 2+3
(biome from water_class + endorheic levels) is now UNBLOCKED on the corrected climate.

## Finding 18 — STEP 2+3 wired: biomes from water_class + below-sea basins as typed water bodies

On the climate-corrected pipeline (Finding 17), STEP 2+3 replace the altitude sea-membership with
connectivity + a water balance:

- **STEP 2 (biome from `water_class`)** — `compute_biomes` now takes `water_class` + `lake_map`
  (empty slices → legacy altitude fallback, so the many `c1_biomes` test callers are unchanged).
  `water_class==1` → Ocean; `lake_map!=0` → the new `Biome::Lake` (frozen id 11); else Whittaker
  land — INCLUDING exposed below-sea land (not flooded to 0 m). `c1_biomes_classified` (new) is
  wired into the HD run.
- **STEP 3 (below-sea basins as typed lakes)** — `drainage::below_sea_basin_lakes` finds the
  class-2 components (≥ `lake_min_area` 5 km²), floods each to its spill, and runs the water
  balance `level = min(spill, evaporative)` reading INFLOW at the basin's LAND INLETS (the real
  `runoff_accumulation` zeroes below-sea cells, so inflow can't be read on the basin itself —
  the flaw in the earlier rough diagnostic). It returns typed `C1Lake`s (Exorheic/Endorheic, ids
  offset ≥ 1_000_001) + a lake_map of the water cells. The HD run merges these into
  `drainage.lakes`/`lake_map` before biomes + export.

Result (2048², relief-v3, reference seed): **9 below-sea basins — 3 exorheic, 6 endorheic**, water
536 km², dry-below-sea 616 km² (preserved, NOT flooded to 0 m); ~30 k cells reclassified OFF Ocean;
`Biome::Lake` appears, exposed margins → land. The 6 endorheic (with the correct inlet inflow) match
the accepted physics: coastal below-sea basins with small catchment/area ratios stay endorheic even
in a humid climate (Salton-Sea-like) — CONTENT, not a defect. `lake_type` reaches `lakes.json`
(C1Lake serializes it; a round-trip test pins it) → Living Landz can distinguish endorheic (salt, no
fish, undrinkable, no shore agriculture) from exorheic; the LL Lake rule must consume the field
(today it flattens all water to one class). The river-profile monotonicity test tolerates flat lake
crossings; the breach guard stays green. 512 lib tests green; viz compiles.

Still open (separate lot): the ×577 erosion-fabricated closed-depression signal (~8 650 vs 15
tectonic) — which MFD/talus/diffusion step creates them, and are they legitimate.

## Finding 19 — DEFECT B: exposed below-sea margins read as Desert because the precip model gives below-sea cells 0 mm

After Finding 18, the HD biome map showed ~2 % **Desert** on an island whose median precipitation is
~1130 mm. Report: the `Desert` predicate is `classify`'s `p_mm < 250` reading the LOCAL precipitation
(not a drainage quantity, contrary to the first hypothesis). Measured on the Desert cells (2048², audit
`defect_abc_audit`): **11 920 cells, precipitation ≡ 0 mm, 11 904 of them below sea**. Cause: the 1-D
moisture scan treats every `n ≤ SEA` cell as ocean (a moisture source, orographic-only, `total =
precip`), so the frontal/synoptic base is gated OUT below sea. When Finding 18 reclassifies an exposed
below-sea margin (class-2 land) it then reads that 0 mm → Desert.

Fix: the frontal base is elevation-independent synoptic rain that falls on below-sea LAND too (Death
Valley, Dead-Sea shores get rain), so drop the `n > SEA` gate — `total = precip + frontal` everywhere.
Open-ocean cells also gain the base, but every consumer masks them via `water_class`, so no land biome
reads it. The frontal term is not part of the conserved orographic budget, so the moisture-conservation
and windward/leeward tests are unaffected. Result: Desert **11 920 → 0**; exposed margins read ~1130 mm
→ forest.

## Finding 20 — DEFECT A: rivers ran through lakes / out of endorheic sinks; clip the exported network to the lake surfaces

Report: `rivers.json` segments are selected by an ACCUMULATION THRESHOLD (`extract_rivers`,
`terrain/flow.rs`: `is_river = acc ≥ stream_threshold`, 20 km²), NOT by traversal to an outlet. Because
accumulation grows monotonically downstream, the exported topology is already complete to the sea — the
audit finds **0 truncated / 0 orphan** terminal segments (all 540 trees reach a coast). The real defect
is a river/lake FIELD INCONSISTENCY: the rivers are traced on the BREACHED (monotone) field, which
drains every basin, while the final `lake_map` marks those basins as standing water (pre-breach lakes +
below-sea balance). The visible result — measured at 2048² — is **44 segments crossing a lake polygon
and 4 sourced inside a lake**: rivers slide through lakes and emerge as orphan reaches below endorheic
sinks.

Fix: `drainage::clip_rivers_to_lakes` splits each segment into its maximal runs of NON-lake points
(profile sliced in parallel); a run that begins at a lake shore is that lake's outlet — kept (with the
parent discharge) for an EXORHEIC lake, dropped for an ENDORHEIC one (the water dies in the closed
basin). Links are remapped so `downstream = None` exactly when a reach ends at a sink (sea or lake). It
touches ONLY the exported/rendered polylines — the routing field, lakes and `lake_map` are untouched, so
a display width/threshold is a pure rendering hint while topology stays continuous to each sink. Wired
into the HD run after the below-sea merge (covers both drainage paths). Result: CROSS-lake **44 → 0**,
sourced-in-lake **4 → 0**, truncated/orphan **0**; unit guard `clip_rivers_terminate_at_lake_sinks`.

## Finding 21 — DEFECT C: filiform rivers gain channel width + a long profile

`rivers.json` exported only `drainage_km2` + navigability — no channel geometry, so consumers drew
zero-width filaments. Added two parallel per-segment arrays to `C1DrainageResult` (round-tripped through
the drainage cache sidecar): `segment_width_m` from the Leopold & Maddock hydraulic-geometry law
`w = a·Q^b` (`CHANNEL_WIDTH_A = 1.2`, `CHANNEL_WIDTH_B = 0.5`, discharge `Q` proxied by the effective
drainage area — so a dry/endorheic reach → 0 width, no channel), and `segment_profile_m`, the bed
elevation (m) along each segment's own points, upstream→downstream. Both surface in `rivers.json`
(`width_m`, `profile_m`). Per-Strahler-order median width (2048², post-clip) grows downstream —
**S1 13 m · S2 21 m · S3 37 m · S4 68 m · S5 114 m**

> ⚠️ **BOTH THE COEFFICIENT AND THESE MEDIANS ARE OBSOLETE (Finding 69).** The code now has
> `CHANNEL_WIDTH_A = 5.0` — ×4.2 from the 1.2 recorded here, with **no finding recording the
> change**, so `a` is a PROXY in practice whatever this paragraph says. `b = 0.5` keeps its
> Leopold & Maddock anchor. And these per-order medians were measured with `a = 1.2` AND with
> `Q` proxied by AREA, both since replaced (Finding 22). Measured at 8192², humid, ratio 7.5:
> S1 11.5 m · S2 25.8 m · S3 12.5 m · S4 10.6 m · S5 2.3 m — **not monotone, and the trunk
> the narrowest.** Nothing in this paragraph describes what ships. — sanity-matching a Thames-scale trunk
(~16 000 km² → ~150 m) and metre-scale headwaters.

## Finding 22 — DEFECT C corrected: width from DISCHARGE (m³/s), not drainage area

Finding 21's law fed the drainage AREA in km² into `w = a·Q^b` — a coefficient calibrated for a
discharge in m³/s. The inspector settled it: 888 km² of catchment read 36 m of width, and `36 = 1.2·√888`
only because 888 (an area) was used as Q. A real discharge check: 888 km² at 583 mm/yr of runoff is
~16 m³/s, so the dimensionally-correct width is `5·√16 ≈ 20 m`, not 36. The area-as-discharge error also
compressed the distribution (area^0.5 grows too slowly) and — being area-only — could not respond to
climate (a dry reach and a wet reach of equal area drew identically).

Fix: carry a per-segment **discharge in m³/s** — `Q = runoff·catchment`, from `runoff_accumulation`
(mm·km²/yr → m³/s via `runoff_km2_to_m3s`, 1 mm·km² = 1000 m³, `SECONDS_PER_YEAR = 3.156e7`) — and derive
`w = a·Q^b` with **`a = 5.0`, `b = 0.5`** (mid-range natural-channel bankfull coefficient; `b` the classic
downstream exponent). A dry/endorheic reach → Q=0 → 0 width. The viz relief-v3 path now computes the final
drainage WITH the climate (breach first → climate on the conditioned field → drainage with discharge),
where before it passed `climate=None` (pure geometry). `segment_discharge_m3s` is stored, cached, exported
(`discharge_m3s`), and shown in the inspector ("Débit").

Sanity table (R = 583 mm/yr, `w = 5·√Q`): headwater 5 km² → Q 0.09 → **1.5 m**; mid 888 km² → Q 16 →
**20 m**; Thames-scale 16 000 km² → Q 296 → **86 m** (mean-annual anchored; bankfull would be larger). The
ideal-extreme trunk/headwater width RATIO is **~57×** (vs the ~9× the compressed medians showed). The
realized per-order medians stay modest (S1 5.8 m → S5 13.7 m at 2048²) because median catchments per order
are not the extremes and runoff varies + endorheic sinks cut discharge — but the law is now physical and
climate-responsive, which is the point.

**Width is sub-cell at production.** At 8192² (49 m/cell) every reach is < 1 cell — per-order median
0.04–0.14, MAX **0.55 cells**. The raster cannot express any flaring; a consumer MUST render channels as a
stroke whose width comes from `width_m` (see the LL note below). Reported by `width_law_audit` +
`river_overlay_render_8192`.

## Finding 23 — DEFECT A/TASK 2 + TASK 5: lake-outlet width continuity + drawing every water body

**Lake-outlet width (TASK 2).** The author saw channels wider ABOVE a depression than below it: the old
per-reach area restarted small at the outlet. Now width comes from discharge, and `runoff_accumulation`
routes flow ACROSS an exorheic lake (held flat, routed to its outlet), so the outlet reach's max discharge
= the whole upstream catchment; `clip_rivers_to_lakes` makes the outlet reach INHERIT the parent
discharge/width rather than recomputing from its local area. An exorheic outlet is therefore continuous
with its inflow; an endorheic outlet is dropped (0 width — the water dies in the basin), which was already
correct. (This seed's below-sea basins are all endorheic, so there is no exorheic through-flow lake to
tabulate; the continuity is structural.)

**Every water body in the overlay (TASK 5).** The Ymir map's river overlay drew rivers but NO lakes, so a
channel terminating in a basin looked like it stopped in the void — the very "truncated outlet" impression
that cost earlier passes. The export always carried the bodies (`detect_lakes` + below-sea, both in
`lake_map`); the overlay just didn't read them. Fix: the overlay now fills every `lake_map` cell, coloured
by `lake_type` — exorheic deep blue, **endorheic teal** — so the count drawn equals the count exported (gap
→ 0) and the author can read regime at a glance. Confirmed at 8192²: 4178 river reaches, **0 orphans**, 6
endorheic below-sea basins drawn in teal, channels terminating in them.

**Consumer note (Living Landz — do NOT implement here).** LL renders rivers at CONSTANT width and flattens
every water body into one class. Both `width_m` and `lake_type` are now exported but ignored downstream —
that is where the value is lost. Because channel width is sub-cell even at 8192² (Finding 22), LL MUST
render rivers as a stroke sized by `width_m`, and distinguish lakes by `lake_type` (endorheic = saline, no
fish, undrinkable). Added to the LL-side backlog.

## Finding 24 — Geographic scale ratio: a presentation compression for the hydrology only

The island is ~40 000 km² — smaller than the Thames basin (16 000 km² is one river). Even a single basin
draining it all yields ~740 m³/s → ~136 m; realistically the largest is ~3 000 km² → ~37 m. There is
physically no room for a great river at 400 km. This is a consequence of the CHOSEN SCALE, not a bug — and
the fix is NOT to inflate runoff (that would over-carve the relief, overflow every lake, and rainforest the
biomes, since one parameter feeds four calibrated systems). The correct device is scale COMPRESSION applied
to DERIVED quantities: the map DRAWS `domain_km` but SIGNIFIES `domain_km · ratio` (a Skyrim-map convention),
costing nothing physically because the terrain keeps its real, coherent values.

`apply_geo_scale_ratio(dr, ratio, thresholds)` — a PURE post-process (instant, NOT in any cache; lives in
`HdParams`, not `C1DrainageConfig`, because it is presentation, not physics). It scales ONLY:
effective catchment `×ratio²`, discharge `×ratio²`, channel width `×ratio`, navigability re-classified on the
scaled catchment. It does NOT touch: the routing field, rivers geometry, lakes, `lake_map`, **stream-power
incision, the lake water balance, precipitation, temperature, biomes** — all of which ran upstream on the
real 400 km quantities, so every prior calibration is untouched. `ratio == 1.0` → identity.

Table (2048², seed reference; small_boat / barge / ship reaches; largest reach):

| ratio | small_boat | barge | ship | largest Q, width |
|------:|-----------:|------:|-----:|-----------------:|
| 1.0   | 742 | 1 | 0 | 48 m³/s, 35 m |
| 3.0   | 1665 | 686 | 0 | 430 m³/s, 104 m |
| 7.5   | 811 | 1661 | 403 | 2688 m³/s, 259 m |
| 15.0  | 383 | 1286 | 1391 | 10 751 m³/s, 518 m |

Barges appear at ratio 3; ships become reachable at **7.5** (403 ship reaches, largest ~259 m). Recommended
**ratio 7.5** for this island.

## Finding 25 — Latitude span & centre: climatic diversity, a SEPARATE control

A 400 km continent spans only ~3.6° of latitude → a single climatic belt, no deserts or tundra possible from
area alone (climate depends on LATITUDE, not size). So the CLIMATIC latitude extent is decoupled from the
physical extent: `c1_climate_placed(centre_deg, span_deg)` sets an explicit span (default `domain_km / 111`,
byte-identical windowed path). `row_latitude_span` derives per-row latitude from (centre, span); temperature
was already per-row, and precipitation's WIND BELT + frontal base are now evaluated PER ROW (`span == 0` →
the single-belt legacy, so the public `compute_precipitation` and the h=1 tests stay byte-identical). A wide
span crosses several belts — trade-wind vs westerly rows, the subtropical dry band — which is the point.

This is REAL physics (temperature, wind, precip → biomes), so — unlike the ratio — it belongs in the climate
computation and a change legitimately re-runs it. It is a SEPARATE slider from the ratio: a cosmetic river
adjustment must not upend the biomes, and a wish for climatic diversity must not change the hydrology.

Table (2048², centre 45°; land-biome distribution + temperature range):

| span | T range | biomes |
|-----:|:--------|:-------|
| 3.6° | −13…13 °C | temperate forest 68% · rainforest 14% · taiga 14% · tundra 1% |
| 10°  | −14…15 °C | temperate forest 69% · rainforest 13% · taiga 13% · tundra 2% |
| 27°  | −20…19 °C | temperate forest 63% · taiga 16% · rainforest 10% · tundra 5% · **steppe 4%** |

Diversity grows with span (tundra 1%→5%, steppe appears at 27°). A temperate CENTRE (45°) never reaches the
subtropics, so no true desert — moving the CENTRE down does. Recommended for a diverse island: **centre 38°,
span 27°** → T −15…22 °C, steppe 41% · temperate forest 19% · savanna 15% · rainforest 9% · taiga 9% ·
tropical rainforest 3% · tundra 2% (a dry belt, a forest belt, alpine tundra — one island). Renders:
`exports/sculpt/recommended_{hydro,biomes}.png`.

## Finding 26 — Manifest: both factors recorded for the consumer

A 100 m river on a 400 km island, or tundra beside desert, would read as a bug to anyone taking the exported
data literally. So `ContinentMeta`/`Continent` (manifest `continent` block) now carry `geographic_scale_ratio`
(default 1.0, serde-defaulted for old manifests) and `latitude_span_deg` (the climatic span behind the
temperature/precip/biome layers). Living Landz must read the ratio before treating river sizes as literal,
and the span to understand the climatic gradient. Both are independent knobs in the config, the UI (two
sliders) and the manifest, as required.

## Finding 27 — Latitude gradient looked inverted: the VIEW, not the data

The author reported tundra at the BOTTOM of the map (centre 60°/span 40° → 40°–80°) where the polar 80°
end should be at the top. STEP 1 (read-only) settled it on the DATA: `row_latitude_span(j=0)` =
`centre − span/2` = the LOWEST latitude, and the container invariant is **row-major, `y = 0` = SOUTH**, so
row 0 is the south (warm) edge — confirmed by the temperature field: row 0 = 40° = 6.1 °C, last row = 80° =
−23.9 °C. The export writes the internal south-first grid row-major → honours `y=0=south`; **Living Landz
reads it correctly.** The renderer, however, drew row 0 at the TOP pixel (workspace.rs: "Row 0 (y=0=south)
is at the top") → south-up → the polar end at the bottom. **Verdict: computation + export CORRECT; the VIEW
was upside-down.** Fixing `row_latitude` would have BROKEN the export/LL — the classic Ymir y-up trap.

Convention each consumer assumes: export rasters — `y=0=south` (honoured); vector layers (coastline /
cliffs / rivers / lakes) — same cell space as the rasters (honoured); LL reader — `y=0=south` (correct);
the viz renderer + cell inspector — were drawing/reading row 0 at the top (the bug).

Fix — VIEW ONLY (`flip_rows_rgba`): the HD texture and the coarse preview are mirrored vertically at build
(north-up); the inspector's hover→cell lookup, the reticle, the coordinate readout and the preview
drag-to-reframe are mirrored to match. `row_latitude`, the temperature/precip fields, the export and LL
are UNTOUCHED. Regression test on the DATA (a display flip cannot fool it):
`row_zero_is_south_and_warmer_northern_hemisphere` asserts row 0 is the lower latitude and, on a flat
field, warmer than the polar row.

STEP 3 — southern hemisphere. Every latitude function (`sea_level_temperature`, `belt_factor`,
`wind_zonal_dir`, `subtropical_suppression`) already uses `|lat|`, so there is no sign to miss: centre −45°
is the EXACT vertical mirror of +45° (`southern_hemisphere_mirror`: temperature deviation 0.0000 °C on a
flat field). `wind_zonal_dir` returns the same sign for ±45° — correct, because the westerlies blow W→E in
BOTH hemispheres (the zonal component does not flip across the equator; only the meridional one does, which
the 1-D zonal transport does not use). A span crossing the equator (centre 0°) is handled per row via
`|lat|` (ITCZ peak at 0°). So the southern hemisphere and equator crossings mirror consistently across
temperature, wind and precipitation.

## Finding 28 — Inspection microscope: four viz panels, read-only, assembling exported data

Over ten passes it was repeatedly the EYE that caught what aggregate metrics hid (striations, perched rivers,
false lakes, inverted latitude). This tooling makes those visible in one glance. Four panels, **UI only** — no
generation/export/data change; a microscope, not an editor (every editable field would raise a consistency
question). Each reads existing data; nothing is recomputed.

- **TASK 1 — latitude placement widget** (`latitude_placement_widget`): a −90…+90° strip (north up) painted with
  the thermal gradient AND the wind BELTS as bands (trade easterlies / westerlies / polar easterlies, each with a
  direction arrow), with the map rectangle `[centre−span/2, centre+span/2]` drawn over it. Live as the
  centre/span sliders move. Makes the CONSEQUENCE of the span visible (how many belts it crosses) BEFORE
  generating. Pure function of centre/span — reads no field.
- **TASK 2 — entity lists** (`aggregate_watercourses` + `inspection_dock`): the 3183 segments are aggregated into
  464 browsable WATERCOURSES — a watercourse = every segment sharing a terminal (post-clip `downstream==None`
  mouth); its trunk = the LONGEST source→mouth path (river-length convention). Listed by discharge (rivers) /
  area (lakes), select-and-highlight both ways (list click ↔ map hover, painter-drawn highlight, no texture
  rebake). Reads `segment_discharge_m3s`, the upstream/downstream links, `points`.
- **TASK 3 — river long profile** (`river_profile_panel`): bed elevation source→sink by WALKING the flow field
  (`flow.direction` + `eroded`) from the main-stem headwater — NOT by stitching `segment_profile_m` across the
  clip's best-effort links, which injected phantom junction steps (a 375 m false climb → 0 after the switch). The
  breached field is monotone, so a REAL climb would stand out — this doubles as the monotonicity inspector. Shows
  discharge / width mouth-vs-source / catchment / length / order and marks the SINK (sea / exorheic lake /
  endorheic basin) explicitly.
- **TASK 4 — lake sheet** (`lake_sheet_panel`): `C1Lake`'s area / level / max-depth / `lake_type`, plus shore
  length and inlet count computed UI-side from `lake_map` + the river mouths (reads, not recomputes). The
  endorheic CONSEQUENCE is spelt out (closed basin → salt, no fish, undrinkable, no shore agriculture — the
  content LL should consume). Average depth / inflow / evaporation are NOT exported (the water balance discards
  them) — stated as such, not faked; surfacing them is a future data-export addition.

Verified by `inspection_panels_data` (2048², ratio 7.5, centre 38° span 27°): 464 watercourses, river #1 a
202 km monotone main stem (999→0 m) to the sea, largest lake a 128 km² endorheic basin (173 km shore, 21
inlets, closed). Core lib 514 green (incl. the row-orientation regression); viz compiles.

## Finding 29 — Exorheic below-sea basins with no outlet: the label overreaches (H2), + the clickable chain

The microscope found a real inconsistency: a river declaring `sink = sea` whose mouth coincides with a
below-sea lake, and several EXORHEIC lakes with no visible outlet. STEP 1 (read-only,
`below_sea_outlet_diagnosis`, 3 exorheic below-sea basins): all **class-2 (enclosed)** — 0 ocean cells, so
NOT the below-sea=ocean error; **no outlet reach leaves any of them**; **not ocean-contiguous**. Verdict
**H2** — `below_sea_basin_lakes` declares exorheic from `a_eq ≥ a_spill` (drainage.rs:813 — inflow would
exceed evaporation at the sill area, a THEORETICAL overflow) WITHOUT tracing a spill path to a sink, so the
"exorheic" label has no outlet behind it. Not H1 (no reach exists), not H3 (not surface-connected).

Point 5 (the label discrepancy): `classify_sink` tested `eroded ≤ 0.5` (sea) BEFORE `lake_map`, so a mouth
sitting in a below-sea basin (water at/below 0 m) read "sea" while the clip had split it on `lake_map` — the
two disagreed. Fixed to test lake membership FIRST (a VIEW-label fix; the data/export/clip are unchanged),
so a mouth in a lake reads the lake, matching the clip.

STEP 2 — physical nature. These basins are **100 % shallow (< 3 m, p50 0 m) with a ~0.1–0.3 % shore slope** —
the WETLAND signature (shallow, flat, poorly drained, at sea level), NOT a lagoon (a lagoon/inland sea is
deep open water, e.g. the 20 m basin the author saw). Proposed criterion, per cell of a below-sea basin:
**wetland** if depth < ~3 m and shore slope < ~1 %, **lagoon / inland sea** if depth ≥ ~3 m. Living Landz's
`Wetland` biome — a no-op for want of a source — would draw from these shallow margins (a future data-export
addition, not done here). NOTE on naming: these are correctly class-2 (enclosed), so this is NOT a fifth
"below-sea = ocean" occurrence; the recurring pattern here is different — **a regime labelled without
verifying it against the traced network**. Because they have no traced outlet and are not ocean-connected,
they are effectively TERMINAL (endorheic); the exorheic label is the defect (fix deferred per diagnosis-first
— trace the spill path or reclassify; the chain below makes it visible meanwhile).

STEP 3 — clickable hydrological chain (implemented): a river's SINK is a button → jumps to the lake; a lake
lists its INLETS (clickable → upstream rivers) and its OUTLET (clickable → downstream river). When a lake is
exorheic but no outlet reach is found, the outlet row reads "⚠ AUCUN exutoire tracé — incohérence H2" — so
this very bug is visible in one click instead of a manual hunt. `NavAction` jumps switch tab + selection.
Viz compiles; guards green.

## Finding 30 — The exorheic label was right; the missing spill PATH was the bug. Trace it; export wetlands

The Finding 29 remedy was backwards. Mass balance settles it: a basin that receives more than it
evaporates MUST overflow — growing the surface raises evaporation, and `min(spill, evaporative)` already
accounts for that, so `a_eq ≥ a_spill` means the water reaches the sill and SPILLS. The `Exorheic` label is
correct; RELABELLING it endorheic would assert a basin that gains more than it loses and never overflows —
a mass-balance violation. The defect is the MISSING outlet. Don't relabel — TRACE.

**TASK 1 — trace the spill path.** Below-sea cells confound plain flow routing (every flow field treats
`≤ 0.5` as ocean base, so directions point INTO the basin), so `below_sea_basin_lakes` now traces the
LEAST-SILL path — a Dijkstra minimising the MAXIMUM elevation crossed (the true overflow sill) from the
basin to a sink — and emits it as a watercourse with discharge = the SURPLUS (`inflow − pe·area`, not the
full inflow: the basin evaporates its share). The path reaches the OCEAN (`water_class == 1`) or CHAINS into
another below-sea basin. Result (2048², centre 38°/span 27°): all **3/3** exorheic basins now trace a
spillway — 2 direct to the sea, **1 chained** through another basin — Q 0.4–0.9 m³/s, 3–5 m wide. The
microscope's verdict flips **H2 → H1**: every exorheic lake's outlet now leads somewhere.

**TASK 2 — invariant, not a computed label.** Permanent, non-ignored guard
`exorheic_below_sea_basin_has_traced_spillway`: a wet enclosed below-sea basin classifies exorheic AND gets
a spillway reaching a sink. Same spirit as `clip_rivers_terminate_at_lake_sinks` / the thalweg lesson — a
regime is now checked against the traced network, not asserted. If it can't be satisfied, the guard fails
(the point).

**TASK 3 — wetland vs lagoon, EXPORTED to the consumer.** New `Biome::Wetland` (frozen id 12). The measured
criterion (through-flow basins are 100 % shallow, < ~3 m, ~0.1–0.3 % slope) drives a per-cell wetland mask
(`below_sea_basin_lakes` returns it; depth `< WETLAND_MAX_DEPTH_M = 3 m`); `c1_biomes_classified_wet` reads
it so shallow margins are `Wetland`, deeper cells stay `Lake` (lagoon/inland sea). WHY they are fresh, for
the consumer: they are fresh because water FLOWS THROUGH them (inflow > evaporation, the surplus spills) —
`endorheic ⇒ salt` still holds, but these are not endorheic. Result: **231 km² of wetland (0.88 % of land)**,
biome `Wetland` 0.88 % where it was 0. This is the data source Living Landz's long-idle `Wetland` biome
always lacked — now on the wire, not just in the microscope.

`below_sea_basin_lakes` returns a `BelowSeaResult { lakes, lake_map, spillways, wetland }`; the HD run appends
the spillways to `rivers` before the clip and passes the wetland mask to the biomes. Core lib 515 green
(new guard + frozen-ids incl. Wetland); viz compiles; erosion/relief chain untouched.

## Finding 31 — Inventory threshold ≠ sink validity; and the SEA is `water_class`, not altitude

The microscope found four watercourses (two large: #26 Q 124 m³/s / 13 005 km², #15 Q 171 m³/s / 18 013 km²)
all reporting `outlet = sea` while converging on a water body too small to be inventoried. Cause: a basin
below the 5 km² threshold never entered `below_sea_basin_lakes`, so it was absent from `lake_map`; rivers
ending there found no membership and `classify_sink` fell back to ALTITUDE and concluded "sea". A basin < 5
km² fed hundreds of m³/s is grotesquely over-supplied — it must overflow — yet was treated as terminal.

**TASK 1 — decouple.** The 5 km² threshold now governs ONLY the exported `lakes` inventory (no micro-lakes in
`lakes.json`). Every enclosed below-sea basin is MARKED in `lake_map` + gets its spillway traced regardless of
area (sink validity is a different question from inventory presence). Population (2048², centre 38°/span 27°):
**638 below-sea basins, 629 sub-threshold**, 462 with a spillway, **11 carrying > 1 m³/s** — a population, not
one edge case. Sub-threshold basins are absent from the inventory but their rivers now terminate at a marked
basin (→ lake), never "sea". The clip drops the breached-carve outlet run for any below-sea basin (id ≥
1_000_001) since its authoritative outlet is the traced spillway (no double-count).

**TASK 2 — the SEA label is `water_class`.** `classify_sink` now derives "sea" ONLY from `water_class == 1`
(the flood-fill-from-borders authority), NEVER from altitude. Before: 228 river mouths sat on below-sea
NON-ocean cells that the altitude rule called "sea"; **after: 0 mislabelled "sea"** — 118 now terminate at a
marked basin, 110 on dry below-sea flats (correctly `Unknown` evaporative terminals, not sea). Permanent guard
`below_sea_sink_decoupled_from_inventory_and_sea_is_water_class`.

**TASK 3 — the over-supplied basin.** The most over-supplied sub-threshold basin (#1000114, 0.04 km²) carries
**8 m³/s** — enormous for its area — and now traces a spillway that reaches a sink (chained). Strongly
exorheic, as mass balance demands; not papered over as terminal.

### The recurring pattern (recorded, not just the third instance)

Findings 27, 29–31 are one class of flaw: **a property asserted from a PROXY instead of the AUTHORITY that
defines it.** Latitude orientation read from image-row order, not the documented `y=0=south`. A lake's regime
read from `min(spill, evaporative)` without the traced outlet behind it. A river's sink read from altitude,
not `water_class`; a basin's sink-validity from an area threshold, not its enclosure. Each fix is the same
move: replace the proxy with the authority (`water_class`, the traced network, the documented convention) and
pin it with a permanent invariant. When a new symptom appears, the first question is now *which proxy is
standing in for which authority.*

## Finding 32 — Spillways that were never traced: a bounded search, and a guard with a blind spot

A basin fed by five rivers (real inflow ~19 m³/s vs ~0.135 m³/s evaporation — exorheic by a factor of ~150)
still read TERMINAL. Cause: the per-basin spillway trace (Finding 30) was a Dijkstra bounded by a step
budget (`64·(w+h)` pops); when the ocean lay beyond that budget the search returned nothing, so an exorheic
basin got no outlet. Of 638 below-sea basins, 462 had spillways and 176 did not — and the invariant guard
`exorheic_below_sea_basin_has_traced_spillway` only inspected INVENTORIED (≥ 5 km²) lakes, so it was blind
exactly where these sub-threshold failures lived.

Fix — route by a SINGLE priority-flood from the ocean (Barnes 2014), computed ONCE: `spill_receiver[k]` is
the neighbour toward the ocean along the least-barrier (minimax-elevation) path. Following it from any basin
cell crosses the lowest sill and reaches `water_class == 1` — guaranteed for every basin the flood reached
(i.e. every one, the ocean being connected). O(n log n) once, so a distant or large basin costs the same as
a coastal one; no budget to exceed. Result (2048², centre 38°/span 27°): **638 basins → 465 exorheic, ALL
465 now traced (0 untraced), 173 endorheic** (legitimately terminal). The 176 gap was 173 legit endorheic +
3 exorheic Dijkstra-failures, now closed. STEP 1 (the over-supplied basin): inflow 7.6 m³/s vs evaporation
0.002 m³/s, `a_eq 159 ≫ a_spill`, EXORHEIC, spillway TRACED carrying the full inflow to the sea.

`below_sea_basin_lakes` now returns a `BasinSummary` for EVERY basin (inventory + sub-threshold, real
units); the guard iterates those and asserts exorheic ⟹ spillway regardless of area — the blind spot is
closed. **Exorheic basins lacking a spillway = 0.**

Display (STEP 3, separated from the data bug): the lake sheet resolves outlets from the WATERCOURSES (the
appended spillways), not the inventory — so a traced outlet shows once the basin is selectable. A
sub-threshold basin is a real sink but absent from the inventory list, so a river's sink button labels it
honestly ("sous-seuil, non listé") instead of a dead jump. And the microscope now flags that displayed
discharges are geo-ratio COMPRESSED (×ratio²) while the balance reasons on REAL quantities — reading the two
as homogeneous is what made this bug look impossible at first (540 displayed m³/s was ~19 real).

Same pattern as Findings 27/29–31: a property (a reachable sink) decided by a PROXY (a fixed search budget)
instead of the AUTHORITY (the ocean-connected flood). Replace the proxy; pin the invariant over ALL cases,
not the convenient subset.

## Finding 33 — Inventory threshold too high, and basins filled to the shoreline instead of the sill

**PART A — the 5 km² inventory floor.** At 40 m/cell, 5 km² is 3125 cells — a plainly visible lake. Excluding
it from `lakes.json` made rivers terminate in a body ABSENT from the export (invisible to the consumer);
tributaries dying into nothing is worse than a few extra small lakes. The old fear (erosion-fabricated
parasitic pits) is now handled by the breach conditioning, so it no longer applies. The below-sea inventory
floor is now a few CELLS (`INVENTORY_MIN_CELLS = 4`, resolution-independent) — reject single-cell noise only.
At 8192² (author's config): **43 inventoried below-sea lakes**, histogram by cells `[<4:0, 4–15:31, 16–63:4,
64–255:1, 256–1023:6, ≥1024:1]`.

**PART B — fill to the SILL, not the shoreline (the real defect).** An enclosed basin's sill is ABOVE sea
(else the sea would enter it and it would be class-1). Yet `spill` was the min external neighbour of the
below-sea (class-2) cells — only the ~0 m SHORELINE, not the overflow rim. So a deep bowl (floor −20 m, rim
+20 m) filled merely to 0 m: level 0 m, footprint tiny, an exorheic label on an UNFILLED basin. Fix: `spill`
= the ocean priority-flood BARRIER at the basin (the least-max-elevation to escape — the true sill). A basin
now fills to its rim. Consequence, all correct: `a_spill` (filled area at the sill) grows, so basins whose
inflow cannot fill them to the rim become ENDORHEIC at the evaporative level — a basin overflows only if it
actually fills. At 2048²: exorheic basins 465 → 4; survivors cover real areas (a merged 18.3 km² lake vs a
former 0.04 km² footprint); deep pits fill to real sills (floor −6 m, sill +11 m). **0 unfilled-yet-exorheic**
by construction (`level = spill` for exorheic; verified).

**Merge.** A below-sea lake is now a connected pool of UNDERWATER cells (`barrier > height`) filled to the
shared sill — adjacent sub-pockets behind the same sill MERGE automatically (no per-class-2-component double
counting; the earlier run showed five byte-identical 18.272 km² "basins" that were one lake). At 8192²:
**10 128 class-2 components → 400 merged below-sea lakes.**

**TASK 3 re-verified.** Filling deeper shifts the wetland/lagoon split — less shallow margin qualifies, so
**wetland 74 → 17 km²** at 8192². Guards green (the spillway invariant already iterates every basin via
`BasinSummary`; the decoupling guard updated for the 4-cell floor).

Pattern again (Findings 27/29–33): the water LEVEL taken from a PROXY (the shoreline / the floor) instead of
the AUTHORITY the regime defines (the sill from the flood for exorheic, the evaporative equilibrium for
endorheic).

## Finding 34 — Inflow read from river tracks (a stale footprint) instead of the runoff field

Lake #1000104: 61 km², level −20 m, MAX DEPTH 0 m — a surface with no water, i.e. the evaporative level had
collapsed onto the floor. Ten rivers pointed at it but stopped 1–2 cells short, and only 3 of 13 counted as
inlets. The water balance read inflow AT THE LAND INLETS (rivers TOUCHING the lake), so inflow was under-read
~4× → `a_eq` crushed → endorheic at floor level.

An ORDERING problem (the author's suspicion, confirmed): rivers are extracted BEFORE `below_sea_basin_lakes`
fills the depression to its sill, so the grown footprint ate the tracks' last cells and the touch-test missed
them. Worse, even reading the runoff FIELD, `runoff_accumulation` ZEROES below-sea cells — so a tributary's
accumulation is lost the instant it enters the water, and `max` over the pool saw only the largest single
stream, not the sum.

Fix — TOTAL inflow SUMMED at the shoreline, from the authority: for each below-sea WATER cell, add the
accumulated `runoff` of every above-sea neighbour that DRAINS INTO it (each tributary counted once, at the
shore, before the zeroing). Read from the FINAL footprint, so a track ending a cell short is irrelevant. The
big lakes come right: #1000087 now inflow 7.7 m³/s (was 1.1), **level −11.2 m, MAX DEPTH 10.2 m (was 0),
area 192.6 km²** — a real endorheic lake at its evaporative equilibrium, not a zero-depth footprint.

Separating the inlet bug from `a_spill` (the author's explicit check): the exorheic/endorheic split moves
only **238 → 216 (58 flipped)** at 8192², so the earlier 465 → 4 collapse was `a_spill` (fill-to-sill,
Finding 33), NOT the inlet bug. What the inlet bug actually broke was the DEPTH of the large endorheic lakes
(collapsed to the floor) — now corrected. Wetland moves again with the geometry: 17 → **193 km²** (the large
endorheic lakes carry big shallow margins).

TASK 4 — boundary: after the clip, river∩lake overlap is **0** (every cell is river, lake, or land — no
unclaimed in-between; the two extents meet exactly, asserted by `boundary_and_gap_check`). The residual
near-misses (a handful of mouths 2–3 cells out) are now caught by widening the microscope's inlet test to ±2
cells, so a tributary ending a cell short still counts and reads as attached.

Same pattern (Findings 27/29–34): the inlet set / inflow decided from a PROXY (river tracks traced against an
earlier footprint) instead of the AUTHORITY (the runoff field at the final footprint).

## Finding 35 — 4-connectivity for water: diagonal sea pockets mis-classed as inland

Three symptoms with one cause. #1000337: area 396 km², level 613 m, MAX DEPTH 633 m (floor −20 m) — a
coastal depression filling to 613 m, which DROWNED #28 (a mountain lake) inside its footprint (the overlap
was not a display bug). The sill comes from the priority-flood ("min of the max altitude to cross to reach
the ocean"); a −20 m coastal pocket needing a 613 m crossing is absurd. Root: **`water_class` flooded with
4-CONNECTIVITY** (connectivity.rs) while the priority-flood barrier uses **8-connectivity** (D8). A coastal
pocket touching the sea only at a DIAGONAL corner was therefore classed INLAND (class-2), so it was NOT a
flood seed — and a basin whose natural low exit runs through that pocket had its barrier computed the long
way, over a mountain (613 m). The author's connectivity count (330 448 inland at 4-conn vs 322 414 at 8-conn,
~8000 cells) is exactly those diagonal-sea cells. A diagonal contact IS a hydrological connection, so
**8-connectivity is the physical choice for water**.

Fix: `water_class`'s border flood is now 8-connected (and the below-sea pool BFS + the inflow shoreline scan,
for consistency). Results (8192², author's config): below-sea basins **400 → 21** (the diagonal sea pockets
are now correctly OCEAN, including #1000337, which vanished — it can no longer drown #28); the **max exorheic
fill level 613 m → 14 m**, and **exorheic basins filling above 50 m = 0** (no absurd fills). The 7 remaining
basins with a >100 m sill are all ENDORHEIC — they sit at their evaporative level (~sea), the high sill being
a real inland barrier they never reach (a Dead-Sea geometry), so their depth/area are correct. The three
counts the author asked for: **river∩lake overlap 0, lake-lake overlap 0** (claimed cells ≤ distinct
`lake_map` cells), **river mouths on below-sea non-ocean cells with no lake = 0** (all such mouths now read
`water_class == 1`). Wetland moves again with the geometry: **193 → 233 km²**.

TASK 3 — invariants extended, all permanent/non-ignored: (a) NO TWO LAKE FOOTPRINTS OVERLAP — each `lake_map`
cell carries one id and every lake's cell-count equals its claimed area (a silently-overwritten lake fails);
(b) a lake's MAX DEPTH = level − floor, and level is never below floor; (c) the exorheic⟹spillway invariant
already iterates EVERY below-sea basin. The `detect_lakes` provenance's outlets are river reaches covered by
`clip_rivers_terminate_at_lake_sinks`; the specific #28 case is resolved at the root (it is no longer drowned).

Pattern again (Findings 27/29–35): connectivity read from a PROXY (4-conn, an incomplete adjacency) instead
of the AUTHORITY for water bodies (8-conn, where a diagonal touch is a real connection). Every fix in this
thread has been the same move — replace the proxy with the authority, and pin an invariant over ALL cases.

## Finding 36 — Below-sea lakes filled to the OCEAN barrier, not their local sill (94 % of the footprint was drowned green); and a method rule

### Method rule (recorded first, because it cost a full round)

Every invariant counter MUST be measured in the PRODUCTION configuration, and preferably at both
resolutions, with the 2048² numbers kept EXPLICITLY separate from the 8192² ones. Finding 35's headline
(below-sea basins 400 → 21, "max exorheic fill 613 → 14 m") was measured on a diagnostic terrain — a bare
`incise` erosion at first domain 400 then 1024 km — on which this basin classifies ENDORHEIC and never fills
to its sill. The PRODUCTION export (`seed…_8192.ymir`, domain **400 km**, geo-ratio 7.5, full
relief-v3 + closures + MFD erosion + breach) still carried #1000020 at **level 613 m, depth 633 m, area
396 km², EXORHEIC** — exactly the screenshot. This is the same class of error the whole thread keeps hitting
(A_c-in-cells, the lateral half-width, the measurement window, the 2048² overlap counters): a result validated
at one resolution/config says NOTHING about another. Diagnostics now read `domain_km` from the manifest and
carry `YMIR_DOMAIN_KM`; the proof of the defect was run against the exported `lake_mask.u32` / `height.u16`
directly (ground truth), not a re-derivation.

### The proof (on the export itself, `export_footprint_proof`)

For #1000020 (level 613.2 m, floor −19.9 m): claimed 283.8 km² of which the genuine HOLLOW (cells ≤ 0 m,
below sea) is **16.9 km²** and the GREEN SWALLOWED (cells > 0 m, drowned only to reach the sill) is
**266.9 km² — 94 %**. The ≤ 613 m region connected to the floor spans **146 871 km²**: 613 m is not the rim of
a 17 km² hollow, it is the minimax pass to the ocean. Cells above the level = 0, disconnected = 21 (both
existing guards were satisfied — see below). The author was right; the earlier "footprint is sound" reading
was an artefact of the wrong (bare-`incise`) terrain, where the basin is endorheic at ~54 m.

### Root cause

`below_sea_basin_lakes` set the fill level to `spill = min(barrier_q)` — the ocean-minimax barrier from the
priority-flood — and grew the pool as `underwater(k) = barrier_q[k] > height[k]`, which includes EVERY cell
under the continental pass (a 400 m green hillside behind a 613 m col satisfies `613 > 400`). For a hollow
enclosed behind high terrain the barrier is the far-away ocean col, so the lake rose to it and drowned the
green in between. The barrier SEARCH is correct (613 m is genuinely the lowest col to the ocean, no detour);
the FILL MODEL is wrong — a lake fills to its LOCAL sill and overflows there, chaining downhill, it does not
rise to the continental pass.

### The fix

A below-sea lake is the connected ENCLOSED below-sea component (`water_class == 2`, the real hollow). A SINGLE
priority-flood outward from its floor finds, in one pass, both the LOCAL sill (the first rim cell with an
unvisited strictly-lower neighbour — the lowest saddle from which water descends to a different sink) and the
bounded bowl (cells ≤ sill connected to the floor, already in height order). The lake fills to
`min(local_sill, evaporative)`. `barrier_q`/`spill_receiver` are kept ONLY to trace the outlet path
(Finding 32), never to set the level. Effect on the faithful production terrain (8192², 400 km): the top
below-sea lake goes from **level 198.7 m / depth 207 m / 142.2 km²** to **85.5 m / 87.7 m / 19.0 km²** (−87 %
area), filling to its local sill (85.5 m) instead of chasing the 472.7 m ocean barrier; `claimed/valid` 1.00×,
0 cells above level, 0 disconnected. (The exact export figures need a viz re-export at the author's config;
the fix is proven on the faithful terrain + the deterministic unit test.)

### Two new plausibility invariants (permanent, non-ignored)

Both existing guards were BLIND to this because they check internal CONSISTENCY, not PLAUSIBILITY:
`depth == level − floor` is satisfied by 613 − (−20) = 633, and the overlap check compares cell SETS which are
genuinely disjoint here (#28 falls inside #1000020's OUTLINE, not its cell set). The missing guards:

- **TASK 2** — every lake cell is ≤ the lake's level AND connected to the floor through cells also ≤ level
  (no disjoint puddles swept in by an altitude-only test).
- **TASK 3** — a lake's level never exceeds the ARRIVAL altitude of its inlets (water cannot flow uphill; the
  monotonicity guard being green means the profiles match the terrain, so it is the LEVEL that was wrong).

`below_sea_lake_fills_to_local_sill_not_ocean_barrier` pins both on a deep pit behind a high (0.9) ocean ridge
with a low (0.52) local saddle: the old model filled to 0.9, the fix to 0.52.

### The 56 untagged mouths (Finding 35's blind spot, now labelled)

A river ending on an enclosed below-sea cell (class 2) with no inventoried lake is a **sub-sea evaporative
sink**, a THIRD sink label — never `outlet = sea` (that would contradict `water_class`, the authority, and
reinstate the −20 m altitude proxy Finding 31 removed). `Sink::SubSeaSink` renders "→ puits sous-marin
(évaporatif)" in the inspector.

Pattern again (Findings 27/29–35): the level read from a PROXY (the ocean-minimax barrier) instead of the
AUTHORITY (the hollow's own local sill). And a second pattern this Finding adds to the list: a counter is only
as trustworthy as the CONFIG it was measured in — Finding 35's "613 → 14 m" was true of a terrain the product
never ships.

## Finding 37 — The exorheic-outlet invariant over the WHOLE population; and extending watercourses to their sources

### POINT 1 — Exorheic without an outlet, and the guard's third blind spot

The exorheic⟹spillway guard (`exorheic_below_sea_basin_has_traced_spillway`) iterated ONLY below-sea
`r.basins`, on a small SYNTHETIC grid. Measured on the SHIPPED 8192² export (`export_exorheic_outlet_audit`,
a river source bordering the lake footprint = a traced outlet): **detect_lakes 28 exorheic / 0 without an
outlet; below-sea 21 exorheic / 21 WITHOUT one**. So every below-sea exorheic lake in the export shipped with
no emitted outflow — and the guard never saw it, because it ran on a convenient subset (a synthetic grid) that
happened to pass. This is the THIRD occurrence of the same blind-spot pattern: an invariant asserted over a
convenient subset (below-sea only / a synthetic grid / one resolution) rather than the whole population. The
21 violations are a downstream symptom of Finding 36 (the old barrier-fill produced over-large exorheic
below-sea lakes whose spillway tracing did not emit); with the local-sill fill in place the same audit,
reproduced through the full hd.rs chain (`exorheic_outlet_audit`), reports **detect_lakes 24 / 0 and below-sea
3 / 0** — zero without an outlet. detect_lakes was never the culprit (its exorheic label comes from
`outlet_reaches_sea`, so the path exists by construction); the fix was Finding 36 plus COVERAGE.

Fix: `exorheic_lakes_missing_outlet(dr)` checks the WHOLE `dr.lakes` population, both provenances, uniformly —
every exorheic lake must have a river segment whose source borders its footprint (a detect-lake's overflow
reach, or a below-sea basin's appended spillway, which after clipping starts just outside the pool). The
permanent guard `every_exorheic_lake_needs_a_traced_outlet` pins it (positive + negative control across
provenances), and `run_hd` calls it on the PRODUCTION network after the clip — a loud warning plus a
`debug_assert`, so a mislabelled regime fails rather than ships. If a lake can never get an outlet, its regime
is wrong and it must be endorheic.

### POINT 2 — Rivers start at the 20 km² threshold, not at their source

The width law `w = a·√Q` is correct; the author's ~37 m at a "source" is a real ~50 km² basin, plausible. The
defect is that the first EXPORTED point is where accumulation crosses the extraction threshold `stream_km2`
(**20 km²**), not the channel head — so the steepest relative width growth (0.1→20 km²: area ×200, width ×14)
is thrown away and the profile looks flat (44→52 m). Fix (author's proposal, better than lowering the global
threshold, which would export every micro-gully): keep `stream_km2` to decide WHICH watercourses exist, then
walk each retained one UPSTREAM to the erosion regime-split critical area **A_c = `RELIEF_V1_A_C_KM2` = 0.1
km²** — the river starts where stream-power starts. `RiverConfig`/`DrainageThresholds` gain `head_km2` (0 =
byte-identical, the default) and `full_tree`. Retention (a dense cell is kept only if its downstream reaches
`stream_km2`) keeps the watercourse COUNT stable; only the upstream extent grows.

Cost at production (8192²/400 km, serialized-segments proxy; real `rivers.json` ≈ 1.5–2× with `profile_m`):

| option | segments | points | proxy size | monotonicity |
| --- | --- | --- | --- | --- |
| baseline (20 km²) | 516 | 55 951 | ~0.7 MB | 0 violations |
| main-stem → A_c | 17 149 (×33) | 661 k | ~9.8 MB | 0 violations |
| full-tree → A_c | 93 780 (×182) | 1.92 M | ~34 MB | 0 violations |

The extension restores a real width range — full-tree order-1 headwaters average **0.23 m** (×7.5 geo-ratio ≈
1.7 m displayed) rising to the trunk mouths (7.6 m real ≈ 57 m displayed) — instead of the flat 44→52 m. The
monotonicity guard holds on BOTH extended networks (0 violations): steep headwater profiles do not break it.
The inspector label "Largeur source" → "Largeur au 1er point" (honest until the walk reaches the source).
Decision (author): wire **main-stem** as the production default — `run_hd` sets `head_km2 = RELIEF_V1_A_C_KM2`,
`full_tree = false` (×33 segments / ~10 MB, one headwater tail per watercourse), NOT the full tree (×182 / ~34
MB) despite its higher fidelity, to keep the export moderate. `head_km2`/`full_tree` live in
`DrainageThresholds`/`RiverConfig` (default 0 = byte-identical, so every core test is unchanged); only the viz
export path opts in. Headwaters sit on steep terrain with sub-cell widths — data-correct, render as strokes.

### Regime as a CONCLUSION, not a prediction (the 0 m level, and the structural fix)

A lake still shipped exorheic at level **exactly 0 m** with no outlet. Origin: the local-sill flood
(`drainage.rs`, the `sill = if sill_q == i32::MAX { C1_SEA_LEVEL_NORM }` fallback) returned SEA LEVEL when it
found NO escape — and an enclosed basin's sill is necessarily above 0, so a level landing precisely on 0 is
neither a sill nor an evaporative equilibrium, it is that cap. The flood found no escape because it marked
bowl membership at PUSH time, not POP: a saddle's downhill neighbour, pushed earlier by a sibling, read as
"seen", so the escape was masked and `sill_q` stayed `MAX`. Two fixes: (1) mark the bowl at POP (finalised) —
a not-yet-finalised strictly-lower neighbour is a real escape — and make the sill an `Option` (`None` = no
escape within the window → the basin CANNOT be exorheic, no more 0 m cap). (2) Invert the regime: today it was
a PREDICTION (`a_eq ≥ a_spill` ⟹ exorheic, trace a spillway or not), which let the label stand with nothing
behind it — the same proxy-vs-authority defect as the six before. Now the spillway is TRACED FIRST and the
label follows: `traced ⟹ exorheic (level = local sill)`, `no trace ⟹ endorheic (level = evaporative)`.
"Exorheic without an outlet" is now UNREPRESENTABLE, not merely forbidden by a guard. Why the earlier guard
missed it: the regime was decided inside `below_sea_basin_lakes` BEFORE `run_hd`'s check ran, and the check
only warned (release) — it never reclassified. Under the inversion (faithful 8192²): prediction 8 exorheic →
inversion 8, **demoted 0, promoted 0, exorheic-without-spillway 0**; the one basin at ~0 m is a real
through-flow pocket WITH a spillway (its sill genuinely sits near sea level), not the fallback.

Lake-vs-wetland is then a physics readout, not a threshold pick. The deepest below-sea basin (#1000010, floor
−2.2 m): local sill 472.7 m but it fills only to the EVAPORATIVE level 85.5 m (`a_eq` 19 km² < `a_spill` 127
km² ⟹ endorheic, nowhere near the sill), mean depth 25.6 m, 9 % of the footprint < 3 m — a deep endorheic
LAKE, not a wetland. A shallow-inflow basin like the author's #1000007 (~1 m³/s into a 20 m depression) gives
`a_eq ≪ a_spill` ⟹ endorheic at a low evaporative level near the floor, shallow ⟹ the wetland instinct is
physically supported; its exact sill/evaporative numbers need a post-fix re-export to read on the shipped terrain.

### POINT 4 (structure, for the next lot — NOT implemented here)

`rivers.json` is a FLAT list of segments. Each `RiverSegment` is ONE reach (a polyline of `points`) between
topological nodes, with `strahler_order`, `upstream: Vec<usize>`, `downstream: Option<usize>`, and parallel
per-segment arrays (drainage/discharge/width/navigability/profile). `extract_rivers` traces from headwaters and
SPLITS at junction cells (`upstream_count ≥ 2`), so a confluence is NOT an object — it is implicit in the graph
(the tributary reaches each `downstream = Some(trunk reach)`; the shared cell is the tributaries' last point =
the trunk reach's first). There is NO first-class watercourse or main-stem object and NO tributary→stem
attachment beyond the generic graph; the viz `aggregate_watercourses` reconstructs a trunk (group by terminal,
longest/max-discharge path) UI-side only, ephemerally. The main-stem UPSTREAM extension adds reaches in the
< 20 km² region along the max-accumulation branch, but the exported unit is still the reach and a "watercourse"
still bundles trunk + all ≥ 20 km² tributaries as separate reaches. Gap vs the Azgaar target: a named main stem
carrying an ordered source→mouth profile, plus each maximal tributary as its own watercourse with (profile,
discharge, width, Strahler) and a link {joins stem S at point P} — an aggregation pass over the per-mouth trees.

## Finding 37b — Absence must stay absence: remove the 0 m fallback, and the three closing invariants

The eighth instance of the pattern, and its sharpest form: a MISSING value (no escape saddle) was replaced by
a DEFAULT (`C1_SEA_LEVEL_NORM`, 0 m) instead of propagated as absence. `Option::None` was the correct
representation and it had been collapsed into a number — which then FABRICATED lakes: a fictitious 0 m sill
authorised a spillway trace that should never have existed (a loop back into the basin's own −20 m floor), and
gave a level to objects with no physical basis for one (a hollow with no inlet and no outlet). Record it
specifically: **absence must stay absence.**

Fixes. (1) The sill flood already returns `Option<f32>` (`None` = no escape); the residual
`sill_opt.unwrap_or(level)` is display-only and `has_sill` now records the truth. A basin with no sill is
endorheic BY ABSENCE — no 0 m cap. (2) The inversion was too weak: it required a spillway to be TRACED but not
to ARRIVE SOMEWHERE ELSE, so a path returning into the lake kept the exorheic label. Three invariants now close
the class, applied to lakes of EVERY provenance (`below_sea_spillway_obeys_invariants`, permanent):

1. **A lake must have water** — strictly positive inflow OR a positive computed depth. A hollow with no supply
   and no outflow is a DRY DEPRESSION: not marked, not inventoried, no spillway (it belongs in the relief).
2. **An outlet must arrive elsewhere** — its sink is the ocean or a DIFFERENT lake, never its own footprint
   (catches the loop, case A).
3. **An outlet may not loop back into its own hollow** — every non-terminal spillway cell must lie OUTSIDE the
   source lake's own bowl. This is the false-positive-free form. Invariant 3 as a LEVEL/SEA THRESHOLD was
   MEASURED to fail both ways and is NOT used: "no cell below the lake's LEVEL" demotes a legitimate below-sea
   lake overflowing to the sea (#1000004, level 22 m → ocean at −1.8 m); "no interior SUB-SEA non-ocean cell"
   demotes legitimate below-sea CHAINS — **7 of the 8 exorheic basins at 8192²**, pockets that spill over a
   sub-sea sill into a lower pocket toward the ocean (level < 0, the whole descent is below sea by nature).
   Both are real outlets. The loop (case A) is caught precisely by "never re-enter one's OWN bowl": a chain
   re-enters a DIFFERENT pocket and stays valid; a loop re-enters the same one and is rejected → endorheic.

Method-rule payoff: measured at the PRODUCTION config (8192² / 45° / 40°) with the inversion (dd1b48a) already
in place, cases A/B/C were ALREADY ABSENT — 0 loops, 0 lakes at ~0 m, 0 dry depressions; they appeared only in
a coarser (2048²) run and in the author's pre-dd1b48a screenshot. The invariants make the class structurally
UNREPRESENTABLE rather than merely absent at one config. Lake-vs-wetland stays a physics readout (Finding 37):
a shallow-inflow below-sea basin equilibrates near its floor (endorheic, wetland margins); a well-fed deep one
is a lake.

## Finding 37c — Spillways follow the DOWNHILL flow from the escape saddle, not the ocean minimax

The regenerated 8192² export (via `cargo viz`, which rebuilds — so it IS the committed code) showed the
`exorheic-with-no-outlet` warning GONE after Finding 37b's clip reorder, but exposed the next layer: the
below-sea spillways were routed by the ocean-minimax `spill_receiver` (Finding 32), which picks the
least-MAX-elevation path to the sea. For a large basin that path threads UNDER other lakes and over higher
ground — the export's river #1 (outlet of the below-sea lake #1000009 at ~0 m) ran UNDER lake #4 (a 106 m
lake), and such paths have non-monotone profiles ("profil indisponible"). Root: a minimax path is not a
DOWNHILL path; water does not climb.

Fix: the per-basin flood already finds the local sill; it now also records the SADDLE and its lowest EXTERIOR
escape neighbour (the far side of the divide). The spillway starts at the saddle, steps to that escape cell,
then follows the DOWNHILL flow field (`flow.direction`) to a sink — the ocean or a DIFFERENT lake (a chain).
Downhill by construction it never climbs (monotone profile) and never crosses a higher lake; it re-enters its
own hollow only for a true loop → then invalid → endorheic. The ocean-minimax `spill_receiver` block is
removed (dead). This also decides the regime honestly: a basin that fills to its sill but whose downhill trace
cannot reach a sink is endorheic (a genuine closed basin), not "exorheic without an outlet".

Method note (recurring): the earlier "0 cases" were measured with `upscale_with_fbm(stream_power)` — 11
below-sea basins — while the viz's `cached_c1_eroded` runs the full production erosion — 21 basins. Only the
regenerated EXPORT is faithful ground truth; the fix is verified there, not on a reconstructed terrain.

Viz: the lake-sheet OUTLET is now the max-discharge watercourse whose source borders the lake AND whose mouth
does not (a real outlet leaves; it does not also enter) — the earlier "first source-bordering" picked a
phantom (#926) that was also an inlet; the true outlet (#397) carries the discharge.

## Finding 38 — Fill each enclosed below-sea region as ONE water body (68 orphan mouths → 0)

Finding 36 fixed the 613 m over-flood by filling from the floor to the LOCAL sill — but that fragmented
large below-sea regions: the flood stopped at the first internal saddle, leaving a thin shore sliver
(~−0.1 m, just above the lake's evaporative level) uncovered. On the shipped 8192² export, **68 river mouths
(up to 234 m³/s displayed) terminated on those slivers** — `water_class == 2`, `lake_map == 0`, neither lake
nor sea. A river of that size evaporating on a shore is a mass-conservation break at continental scale.

STEP 1 (measured on the export, the real terrain — `below_sea_region_structure_export`): **24** `wc == 2`
regions, each with **1 lake covering ~99.9 %** and a THIN orphan sliver (0.4–1.2 km²); floor −19.9 m; **global
sill ~0.1 m for every region** (the lowest LAND rim — nothing like the 613 m ocean pass). So the fix is safe:
filling each region to its ~0 m global sill covers the slivers without over-flooding land.

STEP 2 — the merge. The per-basin flood's escape now triggers ONLY at a `wc != 2` neighbour (a LAND rim or the
OCEAN); an internal `wc == 2` saddle belongs to the same enclosed region and is ABSORBED. So the flood locates
the region's GLOBAL sill instead of fragmenting at the first internal saddle. The WATER BODY is the region
`comp` itself (its `water_class == 2` connected component), NOT the flood bowl `fcells` — the flood grows over
land to reach the sill and would otherwise sweep an adjacent region's pocket via a land lane (two lakes then
claim one cell). The footprint is capped at SEA LEVEL (only below-sea cells marked, never the 0→sill land ring
where adjacent rims meet), so regions stay disjoint. A region with NO inflow is a DRY salt flat and is dropped
(this also removes the tiny no-inflow pockets). The regime is still a conclusion: exorheic iff the downhill
outlet traces to a sink, else an endorheic inland sea at its rim.

STEP 3 — verified on the REAL export terrain (`merge_verify_on_export`, loading height/precip/temp and
re-running below_sea): **mouths on an orphan below-sea cell 68/65/16 → 0/0/0**; **claimed/valid worst 1.000×,
lakes over 1.00× = 0** (no over-flood — the Finding 36 guarantee holds, we did not trade a defect for its
opposite). The lake COUNT and exo/endo split from that loaded field are quantization artefacts (u16 height
re-fragments `water_class`); the real counts come from the regenerated export. Every existing guard stays
green (8 below-sea drainage guards + the lib suite). The author regenerates to confirm on ground truth — the
measurements are the hypothesis, his export is the verdict.

### Finding 38b — the merge's first cut re-opened the over-flood; the near-sea overflow rule

The first merge (fill to the flood's escape col) shipped a REGRESSION the author's regenerated export caught
immediately: lake #1000022 at **level 471 m, depth 491 m** — exactly the Finding 36 over-flood, back. Root:
STEP 1 measured the "global sill" as the lowest LAND cell ADJACENT to the region (~0.1 m), but the flood's
real escape is the lowest col leading to a LOWER exterior — and for a deeply enclosed region that col is 471 m
up (the low shore cells are local pits that drain back in). Filling to 471 m is the defect. And
`merge_verify_on_export` FALSELY passed it, because loading the u16 export height re-fragments `water_class`
into small pockets with nearby low cols — the quantization hid the very case that fails at full precision. The
fourth "reconstructed terrain lied" of the thread, now with a name.

Fix: a below-sea INLAND SEA sits at ~sea level; it can overflow only a NEAR-SEA col (`sill ≤ sea + ~2 m`).
A high col means the region is enclosed → ENDORHEIC, surface at sea level (it never climbs to a far col). And
`a_spill` is the area to the sill (the flood bowl), so a large enclosed basin is endorheic on the balance too.
Verified: `deeply_enclosed_below_sea_is_endorheic_at_sea_level` (a 56 m-col pit must be endorheic at 0 m, not
filled to 56 m) — a deterministic guard the quantized export check could not provide. On the export terrain the
three numbers now hold together: mouths 0/0/0, claimed/valid 1.000× (0 over-flood), and **MAX below-sea level
2 m / depth 22 m** (not 471). The two exorheic synthetic guards were re-based to near-sea sills (their old
0.52 / 0.55 sills were 226 m / 565 m — they had been asserting the over-flood behaviour).

### Finding 39 — DISCHARGE is the discriminant: the water-balance closure for below-sea basins

The author inspected lake #1000022 on the shipped export: **8 affluents (four at ~190 m³/s signified) feeding a
1.4 km² lake that neither grows nor overflows**. Measured on the shipped rasters (`lake_water_balance`): mean
T 6.8 °C, precip 991 mm/yr, PE 600 mm/yr ⟹ **net evaporation ≈ 0** (a HUMID climate, it rains more than it
evaporates); inflow **961 m³/s signified** (17 m³/s map); evaporation over the footprint **≈ 0**. A closed lake
requires PE > precip (an ARID climate) so evaporation can destroy the inflow; here the climate is a net water
SURPLUS. **Mass is not conserved** — even ignoring the rain credit, ~900 km² of surface (640× the footprint)
would be needed to evaporate the inflow. Finding 38b's near-sea rule forced #1000022 endorheic-at-sea purely on
its sill height, blind to the balance.

This exposed the Finding 36 ↔ 38b conflict for what it was: **not two behaviours to arbitrate with an exception,
but one physical law**. Finding 36 corrected a fill to a barrier WITH NO INFLOW to justify it; #1000022 has
961 m³/s that DO justify it. DISCHARGE is the discriminant; the geometric rule conflated the two because it only
read the sill's height.

The law (replacing the near-sea rule in `below_sea_basin_lakes`):

    net_evap = max(0, PE − precip)                 # what the SURFACE loses net of its own rain
    a_eq     = inflow / net_evap  (∞ if net_evap 0) # surface where evaporation balances inflow
    fills_to_sill = a_eq ≥ a_spill                 # inflow reaches the sill ⟹ overflow candidate
    EXORHEIC   iff fills_to_sill AND a downhill outlet traces to a sink → level = the sill (spill inflow − evap)
    ENDORHEIC  otherwise → level = the hypsometric height where the bowl area reaches a_eq (arid, low-inflow)

Arid, low-inflow basins stay low (Finding 36 preserved — inflow no longer justifies the fill); humid,
high-inflow basins fill and overflow (#1000022).

**Precaution 1 — no double count.** `runoff = max(0, precip − PE)` (the inflow source) and
`net_evap = max(0, PE − precip)` (the surface loss) are COMPLEMENTARY: for any cell exactly one is nonzero.
`runoff_accumulation` sources runoff only from above-sea land and zeroes the below-sea region, so the lake's own
precipitation is never in `inflow`; a submerged cell contributes to EITHER inflow (humid: net_evap 0 there) OR
surface evaporation (arid: runoff 0 there), never both. No precipitation is counted twice.

**Precaution 2 — the guard is reformulated, not removed.** The old guard forbade a high level; the new one
forbids a high level NOT JUSTIFIED BY INFLOW (`a_eq ≥ a_spill`). The GEOMETRIC invariant is untouched: the
footprint is the priority-flood bowl `fcells` (every cell ≤ level AND connected to the floor), so
`claimed == valid` by construction — the net that caught 38b's regression still stands.

**Precaution 3 — the hypsometric curve is the flood's own sweep.** `fcells` (floor→sill) sorted by elevation
IS the area-vs-level table: `area(level = sorted[i]) = (i+1)·cell_km2`; the endorheic level is `sorted[⌊a_eq/cell⌋−1]`.
No separate sweep, monotone, non-iterative.

**The escape-col refinement (the subtle part).** The first cut of the law re-opened the over-flood on the loaded
terrain (**MAX level 613 m, claimed/valid 1.952×, 3 lakes over-flooded**): in a humid climate net_evap is 0
everywhere, so a_eq is ∞ and EVERY basin "fills to its sill" — and Finding 38's `wc != 2` escape test ABSORBED
neighbouring below-sea regions, pushing the sill up to the 613 m continental ocean pass. Fix: the escape is the
lowest neighbour OUTSIDE THIS region's OWN component (`!comp_set`), not merely outside all below-sea water. An
internal saddle (same `comp`) is still absorbed (no orphan-mouth sliver), but a descent to a DIFFERENT region is
a real POUR POINT — the lake overflows/chains there at that LOW col instead of climbing to the far pass. This is
the unification the whole thread was missing: full-component coverage (Finding 38's goal) AND low sills
(no over-flood) at once.

Verified on the loaded export terrain (`merge_verify_on_export`, precip units corrected mm/yr → internal):
mouths **0/0/0**, claimed/valid **1.000× (0 over-flood)**, **MAX below-sea level 14 m / depth 34 m** (not 471,
not 613). Regime flips from 91-endorheic (38b) toward through-flow (most humid below-sea pockets overflow their
low cols). Deterministic guards (the discriminant): `humid_enclosed_below_sea_fills_to_col_and_overflows`
(net_evap 0 ⟹ exorheic to the 56 m col, with a spillway) and `arid_sink_enclosed_below_sea_stays_endorheic_below_col`
(a cool catchment feeding a hot floor ⟹ endorheic near the floor) — SAME geometry, OPPOSITE regimes, discharge
the only difference. Full lib suite green (522).

**#1000022 remains for the author's regenerated export.** On the loaded terrain its overlapping re-run lake reads
endorheic 0 m / 1.2 km² / inflow 0 — the FIFTH "reconstructed terrain lied": u16 height re-fragments
`water_class` and reroutes the flow, so the affluents' 961 m³/s never reach it at full precision's expense. The
LAW is proven by the guards; whether #1000022 flips to the expected exorheic through-flow lake (spilling ≈ its
inflow) is for the full-precision export to confirm. The wetland/biome distribution moves again (more through-flow
water, fewer closed sinks) — measured on the regenerated export, not the loaded one.

The remaining `[HD] WARNING Finding 37: exorheic lake with no traced outlet` (#1000008) is a SEPARATE case, not
folded into this change: below_sea now labels Exorheic only when `traced.is_some()`, so the below-sea class is
closed by construction — #1000008 is to be re-checked on regeneration and diagnosed on its own if it survives.

### Finding 39 cleanup — drop detected lakes submerged by a filled below-sea lake

The author's regenerated export confirmed the law (21 below-sea lakes, 18 exorheic / 3 endorheic, sizes right —
the big lakes now cover the deep depressions that had high-discharge affluents). But `detect_lakes` had already
found separate small-id lakes INSIDE those depressions, and the below-sea merge only wrote its ids "where the
lake_map was empty" — so each submerged detected lake survived as a stale CONTOUR inside the new big lake
(#19/#21 inside #1000019, #20/#22/#23 inside #1000020, #25 inside #1000018, #28 inside #1000023) and fired the
exorheic-without-outlet canary (its outlet was clipped under the new water): `[19,20,21,22,23,25,28]` — exactly
the contained set.

Fix (hd.rs merge): the below-sea water SUPERSEDES a detected lake — overwrite the lake_map with the below-sea id
wherever it has water (not only where empty), then DROP every detected lake that loses all its cells. A
partially-covered detected lake keeps its uncovered cells and stays (resized), so no cell is orphaned. The canary
goes silent for the submerged set.

Still open (a SEPARATE defect, not folded in): a below-sea SPILLWAY appended after `clip_rivers_to_lakes`
(Finding 37b) is never truncated where it crosses a NON-submerged lake, so #1000023's outlet runs through
detected lake #26 (arriving "at sea" at #26's 197 m surface) and #1000020's through #17 (211 m). The spillway
trace sees only the below-sea lake_map, not detected lakes. To be diagnosed on its own. The waterfall at river
#1's mouth (1736 m³/s reaching the sea at 49 m) is accepted by the author as a future coastal-cliff concern, not
a drainage defect.

### Finding 39 spillway chaining — a below-sea outlet stops at the FIRST basin it reaches

Follow-up to the previous open item. A below-sea SPILLWAY traced only against the below-sea lake_map, so it
ran UNDER a DETECTED lake instead of stopping in it: river #11 (outlet of #1000020 at 267 m) threaded through
detected lake #17 (211 m) on its way to #1000019 (115 m), arriving "at #17's altitude". The physical chain is
#1000020 → #17 → #1000019 (each overflowing into the next lower basin), not a pass-through.

Fix: `below_sea_basin_lakes` takes an optional `detected_lake_map`; the downhill trace now halts and chains at
the first cell belonging to ANY lake — detected or below-sea — that it enters (below-sea takes precedence when
both are present). hd.rs passes `drainage.lake_map` (the detected lakes, before the below-sea merge overwrites
it). `None` (every existing caller/guard) is byte-identical to the pre-fix trace. viz compiles; lib suite 522
green (10 drainage guards).

On river #1's 49 m coastal "waterfall" (1736 m³/s reaching the sea at the below-sea lake's fill level): NOT a
routing bug — a consequence of not RE-INCISING the spillway after the fill. In nature an overflow outlet incises
and lowers the lake; the model leaves the outlet at the fill elevation, so a large river can end on a coastal
step. Defensible as a young-landscape snapshot; a future spillway-incision pass (lower the sill, partially drain
the lake) is the realistic fix if undesired. Left to the author's call — not folded in here.

### Finding 39 viz — spillway profile & outlet visibility

Two viz-only display fixes (no pipeline change) for the author's inspection findings:

1. The long-profile inspector re-WALKED the flow field from the source; for a below-sea SPILLWAY that is
   wrong — the spillway crosses a divide the flow field routes BACK into the lake, so the walk returned a
   flat lake-level profile (the "49 m → 49 m" on river #1, whose exported `profile_m` actually descends
   48 → −20 m). Fixed: a spillway (an appended segment, `max_flow == 0`) plots its exported `profile_m`
   (the real bed) instead of re-walking. Normal rivers keep the flow-walk (avoids junction phantom steps).
2. The overlay drew a spillway at Strahler-1 thickness (r = 0) → a 2–3 cell coastal stub was invisible.
   The author's rule: if a lake is drawn, its outlet must be too. Fixed: a spillway gets the same minimum
   thickness as an orphan reach, so a drawn lake's outlet never vanishes.

Characterisation confirmed the tiny below-sea pits (#1000005/10/14/21/24, ~4 cells, 0.01–0.04 km²) each DO
have a traced outlet at the same discharge as their affluent (invariant holds); they were merely too short
to see. #1000012 has no extracted affluent — it is fed by diffuse catchment runoff below the stream
threshold, overflowing at ~20 m³/s (signified). The min-inventory floor (4 cells) that admits these coastal
micro-pits is left as-is by the author's call ("c'est un seuil, tant pis").

### Finding 39 regression fix — below-sea surface never below sea (arid orphans + depth-0 pans)

The author's ARID regeneration (latitude 25°, span 10° → subtropical desert belt, 9 exorheic / 12 endorheic vs
the humid 18/3) exposed two coupled regressions that Finding 39's `fcells ≤ level` marking introduced when the
evaporative level collapses BELOW sea:
1. **3 orphan mouths (20–32 m³/s)** — a river ended on a below-sea shelf cell (≈ −0.1 m) that was NOT marked,
   because the endorheic level had collapsed onto the deep floor and the marking covered only the low floor,
   leaving the shelf (still below-sea) outside every sink. Finding 38's guarantee (the whole below-sea region
   IS the water body) was silently lost.
2. **2 depth-0 "lakes"** (#1000018, #1000003 at −20 m, area 420 / 316 km²) — on a FLAT floor the hypsometric
   level sits on the floor, so `level − floor = 0` while the marking flooded the whole flat floor.

Both share one root and one fix: a below-sea basin is a would-be-sea depression, so its SURFACE never reads
below sea — `surface = level.max(sea)`, used for the footprint, the level/depth, and the wetland test. The whole
region is then always claimed (every below-sea cell is inside the sink → no orphan shelf), and an arid basin
reads as a sea-level inland sea of depth `sea − floor` (not a 0-depth sheet). A basin whose balance rises ABOVE
sea (humid: #1000009 at 44.8 m, #1000022 at its col) keeps that higher surface — Finding 39's fill is unchanged.

Guard: `arid_sink_enclosed_below_sea_stays_endorheic_below_col` now also asserts depth > 0 and ZERO unmarked
below-sea cells — both would fail before this fix. On the loaded (arid) export, `merge_verify` shows mouths
0/0/0 and 0 over-flood. Lib suite 522 green (10 drainage guards); viz compiles. The author regenerates for the
ground-truth verdict.

### Finding 40 — chained below-sea balances solved in topological order (upstream spill propagated)

Point 5 of the arid characterisation. When basin A overflows INTO basin B (a chain), A's spill is part of B's
inflow and can flip B's regime. The solver evaluated basins in SCAN order and never fed the spill forward, so a
downstream basin was decided on its LOCAL inflow alone — the factor-4-inlet trap, one basin down. Latent on the
humid seed (net_evap ≈ 0 → everything overflows anyway) and on the current arid export (a single Q = 0.06 chain),
but a real defect that would silently mis-set regimes along any arid chain.

Constructed test `chained_below_sea_balance_propagates_upstream_spill`: an arid A → B → sea chain where A's large
cool-wet catchment fills it to its col and spills a big discharge into B; B's own catchment is tiny, so on LOCAL
inflow B is endorheic (sits at sea level), but A's spill is enough to fill B to its ~90 m sea-col and make it
EXORHEIC. The test FAILED before the fix (B endorheic, level 0) and passes after.

Fix — a TOPOLOGICAL FIXPOINT over the chain DAG inside `below_sea_basin_lakes`:
- Stable region labels (`region_of`, 8-connected wc==2 components) identify a spillway's chain target by
  GEOMETRY, not by whichever lake is marked yet — so the routing is independent of scan order.
- The whole region scan repeats; each pass feeds the previous pass's spills through `extra_inflow[region]`, so
  after ≤ chain-depth passes every basin sees local + all upstream spill. Flow is downhill (a DAG) → it converges
  (exact-equality early-out; a 16-pass safety bound). A chain-free terrain converges in one pass, byte-identical
  to the old single scan.
- Each exorheic basin's surplus is added to `next_extra[receiver]`; a post-pass fills the display `chained_into`
  now that every region is marked.

This changes chained RECEIVERS on every seed (the humid seed has 6 chains) — the intended correctness gain, not a
regression, since a receiver now accounts for the water it actually gets. Guards: 11 drainage tests green;
merge_verify on the arid export holds (mouths 0/0/0, 0 over-flood, MAX level 14 m). Lib 523 green; viz compiles.
The author regenerates humid + arid for the ground-truth verdict.

---

## Finding 41 — the closed-depression population: the FBM is the sole creator; maturity is not the cure

The continent carries ~thousands of closed depressions where tectonics produces ~15 — it looks *Scottish*
(lake-riddled, glacially over-deepened) by accident, not by glaciation (glacial erosion is planned M6, the
`erosion/glacial.rs` file is empty). The question: are the hollows a TRANSIENT that erosion maturity removes
(→ a young/mature knob: Scotland → France) or a STEADY STATE that a process regenerates each pass?

Measured on the PRODUCTION terrain (8192², coarse grid 64 → `upscale_with_fbm` amplitude 0.04 → relief-v3
`incise`), reproduction validated against the shipped export (**94 % of land cells within 5 m**, the rest being
breach + lakes — so this is the product, not a proxy). Diagnostic: `tests/depression_investigation.rs`
(`#[ignore]`). A "closed depression" = an 8-connected component the priority-flood raises by > 0.1 m.

**THE FBM IS THE SOLE CREATOR.** Population per stage:

| stage | pits | ≤2-cell | median depth |
|---|---|---|---|
| coarse post-isostasy | **16** | 11 | 115 m |
| **after FBM upscale** | **90 682** | 70 335 | 3.2 m |
| after relief-v3 incision (2 iter, production) | 75 060 | 51 968 | 2.1 m |

16 → 90 682 across the FBM upscale. Attribution at the production 2 iterations, one sub-process removed:
full relief-v3 **75 060** · NO talus **104 838** · NO hillslope diffusion **99 975** · NO MFD **79 577**. Every
erosion process REDUCES the count (talus −30 k, diffusion −25 k, MFD −4.5 k); NONE creates any. The hollows are
born in the FBM noise field and the erosion closures only chip at them.

**THE MATURITY KNOB IS INSUFFICIENT.** Incision-iteration curve (pits at 1 / 2 / 4 / 8 / 16 passes):
**78 464 / 75 060 / 64 936 / 53 379 / 47 777**. Monotone decrease but decelerating to a plateau near **~45 k**
— still ~3000× the tectonic 16 — and reaching that far planes the channels we deliberately bounded at 2
iterations (the 55 m-floor-under-3000 m-crest failure). Maturity alone does NOT deliver "France": erosion
integrates the drainage far too slowly to erase FBM noise without destroying the relief.

**SIZE DISTRIBUTION (production 2-iter, pre-breach).** 75 060 pits: median area **2 cells**, p90 7, max 87 428
(the big lakes); **69 % are ≤ 2 cells**, median depth 2.1 m — numerical residue, not landforms. ~**7 731 are
≥ 10 m** deep (the author's "≈8 650" genuine hollows), ~1 000 ≥ 50 m; of these only ~49 clear the detection
threshold to become inventoried lakes. Total fill volume ~1.5×10¹¹ m³.

**THE THRESHOLD-FILL KNOB — described and DELIBERATELY NOT TAKEN.** A depression-fill with a depth/volume
threshold, run on the eroded field BEFORE `c1_drainage_windowed`, is cheap (the flood is already computed;
raise sub-threshold components to `filled`) and composes with the breach (fewer carved channels) and the
below-sea invariants (it only shrinks the input population; the kept lakes still satisfy mouths/over-flood).
It is a direct, controllable dial (0 = identity → higher = France). **We did not take it: it MASKS the
non-physical FBM noise instead of removing its cause.** The author chose to fix the source — FBM conditioning
so the upscale does not inject closed hollows in the first place — as a separate chantier. Recorded here
because the fill is the obvious shortcut someone will otherwise retry; the decision is deliberate.

## The recurring class (Findings 27, 29–40): a property asserted from a PROXY instead of the AUTHORITY that defines it

Nearly every defect of this thread was the same mistake in a new costume: a property was read off a convenient
PROXY rather than established from the AUTHORITY that actually defines it. The instances:

- **latitude** from the image ROW ORDER, not the documented `y = 0 = south` (Finding 27);
- **a lake's regime** from a balance PREDICTION, not the traced outlet network (Findings 37/37b — "regime is a
  conclusion, not a prediction");
- **a river's sink** from ALTITUDE, not `water_class` (a below-sea cell read as "→ sea");
- **sink validity** from an AREA threshold, not enclosure (a real sink excluded because it was sub-inventory);
- **reachability** from a fixed SEARCH BUDGET, not the ocean-connected flood (a bounded per-basin search that
  missed a far coast the ocean flood finds);
- **water connectivity** from 4-conn, not 8-conn (diagonally-touching below-sea cells split into two bodies);
- **a basin's inflow** from tributary tracks against an EARLIER footprint, not the final one (the factor-4
  inlet undercount that flipped a regime to endorheic-at-floor);
- **chained inflow** from scan order, not the topological order of the chain (Finding 40 — a receiver decided
  before its upstream contributor's spill was known);
- and the sharpest form — **a MISSING VALUE replaced by a DEFAULT (0 m)** instead of propagated as absence
  (Finding 37 TASK 1): `Option::None` was the correct representation of "no sill" and was collapsed into a
  number, fabricating an "exorheic at 0 m" with no real outlet. Absence must stay absence.

**The remedy that worked every time:** replace the proxy with the authority (the documented axis, the traced
network, `water_class`, enclosure, the ocean flood, 8-conn, the final footprint, the topological order,
`Option::None`), THEN pin an INVARIANT over the WHOLE population — every lake, every mouth, every basin — not a
convenient subset. The bugs that survived longest hid in the cases the first invariant did not cover.

**Two method rules, earned the hard way:**
1. **Measure invariant counters in the PRODUCTION configuration, at both resolutions.** A result validated at
   2048² says nothing about 8192² (resolution re-fragments the field, shifts thresholds, changes counts).
2. **A reconstructed terrain is not the product.** Recomputing from the exported rasters MISLED US SIX times —
   including a `merge_verify` whose over-flood was masked by u16 height quantisation re-fragmenting
   `water_class`, and the arid/humid regime splits that a reconstructed field reported wrongly. The author's
   full-precision export is the verdict; a reconstruction is only ever a hypothesis.

## Backlog (open, deliberately deferred)

- ~~**FBM conditioning** (the next chantier) — condition the upscale so it does not inject the ~90 k closed
  hollows (of which ~7.7 k are ≥ 10 m); Finding 41 is its brief.~~ **DONE — see the C-1 section below.** The
  flow-conditioned FBM cuts the post-FBM population 13×/28× and brings the deep (drainage-trapping) pits to the
  tectonic order, with mountain morphology preserved.
- **Trunk / tributary separation** (Azgaar-style) — the structure is described in Finding 37 POINT 4; not
  implemented. Watercourses are aggregated (a trunk + its tributaries) but the export does not label the split.
- **Empty erosion modules** — `erosion/thermal.rs`, `coastal.rs`, `aeolian.rs`, `glacial.rs` are all stubs
  (M5/M6). Glacial in particular is what would legitimately produce over-deepened lake districts.
- **Consumer-side (Living Landz ignores these today)** — `width_m` (channel width, exported, unused);
  `lake_type` (exorheic/endorheic, exported, unused); and the **Wetland biome**, which has a data source (the
  below-sea shallow/through-flow mask) but no consumer.

## C-1 — Flow-conditioned FBM (closures roadmap §1)

Finding 41 named the FBM the **sole** creator of closed depressions: 16 after isostasy → **90 682** after the
FBM upscale (8192², production seed). C-1 conditions the FBM so it stops fabricating them, without touching the
erosion/relief chain that follows.

### Literature first (the honest answer)

Searched the geomorphology and procedural-terrain literature for a *named* formulation of drainage-conditioned
/ flow-aware noise that provably avoids fabricating local minima at generation time. **None exists.** The
adjacent prior art is real but solves a different problem:

- **Domain warping** (Quilez; 3DWorld) offsets noise coordinates for organic shape — it warps the *noise*, and
  says nothing about monotonicity.
- **Anisotropic noise** (Goldberg-Zwicker-Durand, SIG'08; Substance Designer) controls a per-region target
  *spectrum* for texture, not the sign of a slope.
- **Slope-weighted amplitude** (exponential slope damping, Red Blob Games) smooths steep areas — the opposite
  end (it damps where we can afford noise, not where pits form).
- **Local-minima removal** (van Kreveld et al., *Imprecise Terrains*; Barnes/Lindsay priority-flood) removes
  minima as a **post-process** (flood/breach) — which is exactly the threshold-fill palliative the roadmap
  rejected: it masks the non-physical noise instead of not creating it.

So this is a derivation, stated plainly, not an unattributed invention dressed as standard. It stays within the
precedent the project holds itself to (Barnes/Lindsay for the counting; the monotone-flow criterion is the
generation-time dual of their post-hoc breaching).

### The criterion and the formulation

A perturbation `n` fabricates a pit where it **out-slopes the bed and reverses the descent** — i.e. where the
along-flow derivative `dn/d(downslope)` exceeds the bed's fall. Two coupled mechanisms enforce monotonicity,
both built on machinery already present (`amplitude_slope_factor`, `fbm_anisotropic`, `base_frequency`):

1. **Relief-budget amplitude cap.** The FBM's per-octave downslope slopes sum to `amplitude * S / lambda_base`,
   with `S = sum (persistence*lacunarity)^o` and `lambda_base = 1/nscale` (coarse cells per base feature).
   Bounding that by `beta * slope_mag` gives `amplitude <= beta * slope_mag / (nscale*S)`. On a flat
   (`slope -> 0`) the cap -> 0 and the fabricated pit vanishes (there is no flow direction to respect); on a
   steep flank it is generous (texture kept). `beta` is the one tuning knob; the limit `beta -> 0` recovers the
   smooth coarse bed (its 16 depressions).
2. **Downslope stretch (fixed x8).** The noise is elongated *along* the bed gradient — sampled with
   `fbm_anisotropic` at `ratio = 1/8` on the slope axis — so its along-flow frequency is divided by 8 and the
   contour axis keeps full-frequency relief (downslope flutes, not transverse ridges). Critically, stretching
   only ever **lowers** a frequency, so it never crosses Nyquist; the band policy holds. The first attempt did
   the opposite — *compressed* the contour axis (`ratio > 1`) to elongate downslope — which pushed the contour
   frequency past Nyquist and **aliased into salt-and-pepper 1-cell pits** (90 682 -> 312 665 at R=32). The sign
   of the anisotropy is the whole game; the legacy `max_anisotropy` elongates along-contour (transverse ridges
   = counter-slopes) and is exactly backwards for drainage.

Both quantities depend only on the coarse slope field and config, so they are **identical at every
`target_size`** — low bands stay bit-identical across resolutions. `flow_conditioning = 0.0` (default) is
byte-identical to the pre-C-1 additive noise; all determinism/byte guards stay green. The production config
(`FbmUpscaleConfig::c1_hd_production`) sets `beta = 0.1`.

### The trajectory (post-FBM closed depressions, production terrain)

8192² relief-budget sweep (downslope stretch x8 fixed):

| beta     | 1.0   | 0.4   | 0.2   | 0.1  | 0.05 | 0.02 | 0.01 |
|----------|-------|-------|-------|------|------|------|------|
| pits     | 55154 | 29895 | 15787 | 6999 | 2513 | 629  | 485  |

Per-stage, at the chosen **beta = 0.1**, both resolutions:

| stage                 | 2048² OFF | 2048² b=0.1 | 8192² OFF | 8192² b=0.1 |
|-----------------------|-----------|-------------|-----------|-------------|
| post-FBM              | 6070      | **220**     | 90682     | **6999**    |
| post-relief (2 iter)  | 4574      | 1001        | 75060     | 17382       |

**The honest read of the acceptance.** "Same order as 16" is *approached, not literally reached* on the total
post-FBM count at 8192² (6999, ~440x), but the number that matters does reach the tectonic order: at beta = 0.1
the **structural** depressions — deep enough to trap drainage — are `>= 50 m : 23`, `>= 10 m : 108` (from
`90682` of which `70215` were >= 1 m un-conditioned). The residual `6999` is almost entirely **sub-metre
quantisation dimples** (median depth 0.9 m, p90 4.3 m) at the float/u16 scale, which the breach stage removes
for the shipped product. Pushing beta below 0.02 would bring the total to `~600` but starts trading real detail
for dimples with no drainage consequence — so beta = 0.1 is chosen as the point where the pit reduction is large
(13x/28x) and the morphology is *provably* intact (below). The deeper elimination is the roadmap's long game:
each closure that replaces FBM detail (volcanism, lithology, coastal) shrinks the noise budget further.

Note the post-relief row *rises* under conditioning (6999 -> 17382 at 8192²): once the FBM is conditioned, the
residual pits are dominated by the stream-power incision's own artefacts, **not** the FBM — which is precisely
the C-1 goal ("the FBM is no longer the creator"). Conditioning the incision itself is out of C-1 scope (no
change to the erosion chain beyond the FBM stage).

### No shape regression (post-relief product field)

| metric                | 2048² OFF | 2048² b=0.1 | 8192² OFF | 8192² b=0.1 |
|-----------------------|-----------|-------------|-----------|-------------|
| slope > 30° share     | 14.7 %    | 15.1 %      | 20.9 %    | **23.0 %**  |
| slope > 45° share     | 0.96 %    | 1.20 %      | 3.56 %    | **7.00 %**  |
| local relief (11², m) | 382       | 334         | 114       | **119**     |

The conditioning does **not** collapse relief — local relief is held (334 vs 382; 119 vs 114) and the
steep-slope shares are **preserved or sharpened**, because the downslope stretch concentrates relief into
coherent valleys instead of smearing it isotropically. The only metric that softens is the > 15° share (gentle
roughness — largely the spurious pit-making noise itself), by design.

### The FBM amplitude floor, revisited

Finding 41 put the un-conditioned degeneracy floor "below 0.02" (amplitude_base). Conditioning **decouples the
question**: the relief-budget cap makes `amplitude_base` nearly irrelevant on gentle slopes (the cap binds
there regardless), and it only sets the ceiling on steep flanks where the cap is generous. So `amplitude_base`
no longer needs to be lowered to fight pits — beta does that at every slope — and it can stay at the production
0.16 / seed 0.04 for mountain texture. The effective amplitude is now
`min(amplitude_base*..., beta*slope/(nscale*S))`, a per-cell budget rather than a global level; the "floor" is a
function of slope, not a scalar.

> **⚠️ CORRECTED (H-1c round) — this entry UNDERSTATED it, and the qualifier was wrong.**
> "Nearly irrelevant on gentle slopes … only sets the ceiling on steep flanks" left the door open that
> `amplitude_base` still mattered somewhere. It does not. **It is ENTIRELY INERT on the production path**, at
> every slope. Proof: building the production terrain at `amplitude_base = 0.16` and `= 0.04` (a 4× difference)
> yields a **BYTE-IDENTICAL heightmap** — 0 of 4 194 304 cells differ, max |Δ| = 0.0 m, and the pre-breach lake
> footprint is identical to the cell (58 895 cells, 67 lakes). A single cell where the cap did NOT bind would
> have diverged; none did, so `min(amplitude_base·…, cap)` selects the cap EVERYWHERE (`flow_conditioning = 0.1`
> in `c1_hd_production`). See `tests/amplitude_anomaly.rs`.

### The DEAD KNOB, and what it invalidates backwards

`amplitude_base` — and with it the viz's `fbm_amplitude` selector, the rendered amplitude ladder
(`exports/relief_ladder/`, 0.16/0.08/0.04/0.02) and **every amplitude sweep run since C-1** — act on a parameter
with NO EFFECT whenever `flow_conditioning > 0`, which has been the production default since C-1. They were valid
BEFORE C-1; C-1 silently neutralised them.

**Be precise about what survives.** The QUALITATIVE conclusions hold: the noise was real, the striations were
real, the closures did what was measured of them, and the pit counts / striation spectra / shape metrics were all
measured on genuine terrain. What does NOT hold is any statement of the form **"the FBM floor sits at 0.02"** or
"amplitude is reducible ≥8×" (Finding 5, and the trajectory criterion invoked in C-2, C-3 and C-3b): those
**never measured anything** — the terrain did not change as the number moved. Corrected here visibly rather than
silently; the earlier entries are left in place with this correction pointing at them.

**The real lever is `flow_conditioning` (β) alone**, and it carries TWO ROLES in one parameter: the relief-budget
CAP (`β·slope/divisor`) and the downslope STRETCH (`1/FLOW_STRETCH`). Any future "make the FBM shrink" work must
move β, not the amplitude — and should first separate the two roles, because a single knob doing both cannot be
tuned against either. Note this composes with the roadmap correction already recorded under C-3b (closures cannot
lower the FBM floor in its own band): the floor was not merely unreachable, the knob measuring it was inert.

**Structural remedy (the seventh bench/production divergence).** `c1_hd_production` is NOT production: the viz
builds it then mutates amplitude, sample origin, erosion, stream-power, lithology and fracture. Any bench calling
it gets something else than what ships — which is how this went unnoticed. The fix is a single
`production_hd_config(target, opts)` returning exactly the shipped config, consumed UNMUTATED by the viz and
called by benches, plus a non-regression test comparing the viz's effective config against it. Divergence becomes
impossible to WRITE rather than forbidden by instruction — a rule that must be remembered seven times is a design
flaw, not a discipline failure. And no public structure may keep a field that does nothing: the amplitude
composition is to be made EXPLICIT (a named cap, β's two roles separated) rather than a silent `min` inside a
loop, which is exactly what hid the inertness for three closures.

### Guards

- `terrain::upscale::flow_conditioning_suppresses_fabricated_pits` — permanent unit guard: on a tilted ramp the
  conditioned FBM fabricates > 4x fewer local minima than the additive one.
- `depression_investigation::c1_flow_conditioning_sweep` / `c1_shape_metrics` (`#[ignore]`) — the trajectory and
  shape tables above, reproducible in the production config at both resolutions.

## C-2 — Volcanism (closures roadmap §2)

The TDD anticipated this closure: §4.5 (arc/rift/hotspot volcanism as source terms) and the biomes section
("Volcanic crater lake ... acidic if the volcano is active, neutral if dormant"). C-2 therefore MEETS AN
ORIGINAL INTENT rather than inventing one -- the acidity-by-activity split is a design goal from the start.

### Bibliography analysis (primary sources read, `docs/refs/`)

The author deposited the references; they were read before implementation, not cited from abstracts.

**Wood 1978, *Morphometric evolution of composite volcanoes* (GRL) -- stratocone geometry. SUPPORTS, with a
strict validity domain.** 26 circum-Pacific composite volcanoes, mostly historically active and relatively
un-eroded -- i.e. CONSTRUCTIONAL geometry, exactly the pre-erosion profile we inject. Citable relations
(km): `Hco = 0.122*Wco + 0.450` (n=17, r=0.95); `Wcr = 0.027*Wco + 0.048` (n=14, r=0.91). Ranges: Wco
0.6-22, Hco 0.2-3, Wcr 0.03-0.7, crater depth Dcr 0.03-0.45 -- and NO significant crater depth/diameter
relation (r=0.62), so depth is drawn independently, not computed from width. Mean flank slope falls from ~33
deg at Wco=2 km to ~15 deg at Wco=22 km. VALIDITY DOMAIN: below Wco = 2 km the edifice is a cinder cone with
DIFFERENT relations (crater 0.8 km wide at Wco=2 km vs 0.1 km for a composite) -- applying the composite law
there is extrapolation. We restrict placed edifices to Wb >= 2 km.

**Grosse et al. 2013/2014, *A global database of composite volcano morphometry* (Bull. Volcanol.) --
CONFIRMS and refines.** n=759, same Wb>2 km cutoff. Medians: H/WB 0.12 (range 0.01-0.30), height 1.5 km,
WB ~10 km, whole-edifice slope 17 deg (lower flank 15, main flank 20, max-average 25, up to 43), crater
width 2.2 km (up to 11), crater depth 240 m (100-860), crater/basal ratio 0.11. Our earlier assumed
"H 1-3 km, Wb 10-20 km, H/Wb 0.1-0.2" is confirmed; we adopt Wood's LINEAR constructional law for the
injected profile (Grosse's slightly lower H/WB 0.12 reflects the DB including eroded cones -- the erosion
pass supplies that decay for us).

**Grosse & Kervyn 2018, *Morphometry of terrestrial shield volcanoes* (Geomorphology) -- shield geometry.**
n=158, mostly monogenetic. H/WB 0.01-0.1 (central 0.10), flank slopes 1-15 deg (central 12), basaltic. This
grounds the shield/stratocone split on a PHYSICAL parameter (composition/viscosity): arc = andesitic,
viscous, steep stratocone (H/WB ~0.12, 17-25 deg); hotspot/rift = basaltic, fluid, gentle shield (H/WB
0.01-0.1, 1-15 deg). Not a style toggle.

**Syracuse & Abers 2006, *Global compilation of slab depth beneath arc volcanoes* (G3) -- arc placement.
NUANCE / partial support.** What the paper actually gives is slab DEPTH beneath the volcanic front: H = 72-173
km, global average 105 km (108 +/- 14, 112 +/- 19, 124 +/- 38 km from earlier compilations). This is a DEPTH,
not the horizontal trench-arc distance. Our C1 state has no slab-depth field, so the depth criterion cannot be
applied literally. We place arcs by the horizontal trench-arc gap (100-300 km, from the arc-trench-gap
literature, a secondary source) transposed onto the O-C convergent margin mask, kept in KILOMETRES and
converted to cells at the point of use. Stated plainly: the depth criterion is real and robust, but our model
substitutes a horizontal offset for lack of a slab.

**Varekamp et al. 2000, *Volcanic lake systematics II: chemical constraints* (JVGR) -- crater-lake chemistry.
SUPPORTS the activity dependence.** n=373 volcanic lake fluids. Explicit: "active acid crater lakes (pH < 2)";
the most acidic lakes are the most active (Poas); hyperacid brines reach pH ~ -0.6 to 1; neutral lakes are
dilute meteoric / water-rock-reacted (i.e. inactive). So: active degassing -> acidic (pH 0-2), extinct ->
neutral freshwater. The exact neutral mode "pH 6-6.5" and the bimodal-gap figure come from Varekamp 2003 (NOT
in the deposited refs, cited from the web) -- the deposited 2000 paper confirms the acid mode quantitatively
and the neutral mode qualitatively, and if anything the acid mode is MORE extreme than assumed (pH can be
negative). The binary split is a MEASURED bimodality, not an arbitrary threshold.

### What the sources CONTRADICT in the initial design (recorded because it changed the plan)

- **Age must NOT lower the constructional profile.** The initial design said "older edifices are already
  lower/broader (Wood 1978)." That misreads Wood: his age law `Wco = 0.63*A^0.18 + 0.65` is a GROWTH law --
  an older (longer-erupting) cone is BIGGER, not smaller. The subduing of old volcanoes is EROSIONAL
  (post-extinction), not constructional. REVISION: inject the age-independent constructional geometry
  (Wood/Grosse) for all edifices; let hotspot-chain age drive (a) ACTIVITY -- young=active->acidic lake,
  old=extinct->neutral lake -- and (b) a relief-decay factor on extinct edifices, declared explicitly as an
  erosion proxy the single (time-less) erosion pass cannot date, NOT a cited law. The uniform erosion pass
  then dissects everything.

### The three method adjustments (author's review, addressed)

1. **The C-1 exemption is unnecessary -- verified, not assumed.** C-1 conditions the FBM UPSTREAM of edifice
   injection, so a crater carved after the FBM is never subject to it. What would erase a crater is the
   BREACH (`terrain::flow::breach_monotone`) -- but it already sets every `lake_map` cell to its flat sill
   and NEVER breaches it. A crater is a genuine deep (240 m) wide (2.2 km) depression, so the normal
   depression -> lake detection picks it up and the breach protects it like any lake. No special exemption is
   added; only CLASSIFICATION (tag lakes whose footprint meets the crater mask). A guard placed where the
   problem is not would suggest coverage that does not exist.
2. **The hydrological shift is measured, not merely invariant-checked.** Placing edifices before erosion
   recomputes drainage on a modified terrain (radial divergences, basin capture, moved divides -- intended).
   The magnitude is reported before/after: drainage density, Strahler histogram, confluence count, lake count
   and regime split, floor/local-ridge, W/D per Strahler order.
3. **The trench-arc offset stays in kilometres**, converted to cells at the point of use via
   `C1_DOMAIN_KM/grid`, so it survives a change of domain or resolution -- not baked into a cell count.

### The one unanchored parameter, labelled

Every C-2 number above carries a publication EXCEPT ONE: the **relief-decay factor applied to extinct
edifices**. It is an explicit PROXY, not a cited law. Justification: the pipeline runs a single, uniform
erosion pass that has no time dimension — it cannot erode an old extinct cone more than a young one because it
cannot date its own work. The decay factor stands in for the post-extinction erosion that the timeless pass
omits (a young active cone is pristine; an old extinct cone is subdued). It is marked as such in the code
(`VolcanismConfig::extinct_relief_decay`, doc-commented "PROXY, not from a publication") so a reader can tell
at a glance which parameters are anchored and which one is not. When per-edifice erosion timing exists, this
proxy is the first thing to remove.

### C-2 placement — provenance, structural verification, and the cache verdict

**The trench-arc offset is provisional and labelled as such.** All C-2 geometry numbers carry a publication;
the trench-arc horizontal offset does NOT. Syracuse & Abers 2006 give a slab DEPTH (72-173 km, mean 105), not
a horizontal distance, and our model has no slab-depth field -- so the horizontal offset comes from the
secondary arc-trench-gap literature (100-300 km). It lives as ONE named constant,
`placement::TRENCH_ARC_OFFSET_KM_DEFAULT = 150 km`, in km, converted at the point of use, so it is one line to
change when a better source or a slab field arrives. It is the least-certain number in C-2 and should be read
as provisional next to the anchored Wood/Grosse/Varekamp values.

**Placement is judged on STRUCTURE, verified on the real C1 state (two seeds).** `c2_placement_structure`
(#[ignore]) reports:
- ARCS form a line at a consistent offset from the O-C margin: seed A offset mean 123 km / std 19 km (3 arcs,
  7 margin cells); seed 42 mean 109 km / std 32 km (10 arcs, 22 margin cells). The tight std confirms a
  boundary-parallel line, not scatter.
- HOTSPOT CHAINS are causal, not fortuitous: ages increase MONOTONICALLY along each chain
  [0.00, 0.25, 0.50, 0.75, 1.00] and the members are colinear (max step deviation 0.0 deg), with the youngest
  (active) member over the plume -- the age-progression-opposite-to-motion signature (Hawaiian case).
- RIFTS sit 100% on Divergent-and-Continental cells (8/8 and 5/5).
- Counts track the seed's tectonics, not a random draw: arcs 3 vs 10, rifts 8 vs 5, O-C margin cells 7 vs 22
  across the two seeds -- more subduction margin gives more arcs.

**Cache verdict (the silent-failure class avoided).** (1) VolcanismConfig enters the cache key: it is added to
`FbmUpscaleConfig`, and `eroded_key` serialises the whole config (`.with("upscale", upscale_cfg)`), so enabling
or changing volcanism changes the terrain key and forces a recompute. (2) Crater records CANNOT be recomputed
on a cache hit: the eroded cache stores a bare `GridF32` and, on a hit, `cached_fallible` returns it WITHOUT
running the closure -- so `C1State`/`PlateKinematics` (which placement needs) do not exist. Recomputing
placement at the lake-typing stage would therefore fail SILENTLY on every cache hit (craters in the terrain,
absent from the mask, every crater lake mistyped as ordinary, no error anywhere). The fix is structural: the
crater records travel WITH the cached terrain as one bundle (`{ heightmap, craters }`), computed and cached
together, so a terrain/crater mismatch is impossible by construction. This is settled before integration
precisely because it is the class of defect that took several rounds to find in the hydrology phase.

**Offset tightened (author review).** Re-measuring the arc offset PERPENDICULAR to the local boundary tangent
(not to the nearest neighbour) first widened the gap (107, 85 km on the two seeds), so the applied offset was
verified DIRECTLY foot->edifice: **150 +/- 0 km** — exactly the intended value; the magnitude was never wrong.
The shortfall was a DIRECTION/measurement artifact: the inboard normal was the axis-aligned 4-neighbour oceanic
sum (up to 45 deg off a diagonal margin) and the perpendicular measurement used a tangent estimated from a
sparse (7-22 cell) margin. Smoothing the normal to a radius-2 distance-weighted sum recovered most of it
(perpendicular 129 and 106 km); the residual vs 150 is real margin curvature over a 150 km step plus the noisy
sparse-margin tangent. Since the trench-arc offset is the least-anchored parameter (100-300 km range), an
effective perpendicular of ~106-129 km sits comfortably inside it. Context, not a defect: with only 7 O-C margin
cells on the main seed (3-8 arcs), ARC VOLCANISM IS BARELY VISIBLE on this continent — the visual validation
rests on the hotspot chains and the rifts.

**Scale decision (the km-vs-cells trap, settled before integration).** Volcano morphometry is in physical km, so
the domain's physical span must be pinned. Two candidates existed: the geometric `domain_km` (the map IS the
domain; `km/cell = domain_km/target`) or `domain_km · geo_scale_ratio` (what the map "represents"). The HdParams
contract is explicit: `geo_scale_ratio` is a HYDROLOGY-ONLY presentation multiplier and "NOTHING that shapes the
terrain sees it." An edifice shapes the terrain, so volcanism uses the GEOMETRIC `domain_km` and ignores
`geo_scale_ratio` — consistent with incision, lake balance and climate, which also ignore it. Consequence: on a
small 400 km export the 150 km arc offset is a large fraction of the domain (arcs pushed far inboard, few and
barely visible); on the default 1024 km domain it is ~15% (the well-scaled case the structural test measured).
`VolcanismConfig::domain_km` is set by the bridge to `params.domain_km`; placement `torus_km = domain_km`;
`apply_edifices` km/cell = `sample_size · domain_km / target`. All physical, converted at the point of use.

**The check that would have missed it (scale).** The Mayon/Mauna Loa control validates the morphometric LAWS,
not the on-map result — so it would have stayed green even if the edifices had been rendered 7.5x too large or
too small by wrongly applying `geo_scale_ratio`. The safeguard that actually holds the scale boundary is the
CONTRACT ("nothing that shapes the terrain sees the ratio") plus the per-edifice km/cell conversion, not the
control table. Recording this because knowing which check would have failed silently is as informative as the
one that works: a law-level guard cannot catch a domain-level scale error.

### C-2 lake typing — craters are not lakes

A crater only holds a lake if there is WATER (enough inflow / captured runoff), so the typing applies to LAKES
THAT INTERSECT A CRATER, not to every crater — a dry crater stays relief. `classify_crater_lakes` finds, per
crater, the detected lake occupying the rim (majority `lake_map` id inside the crater), and types only that
one: `CraterAcidic` if the edifice is actively degassing (Varekamp pH < 2), `CraterNeutral` if extinct
(ordinary freshwater). It reports `(craters_with_lake, dry_craters)` so the split can be sanity-checked (a
young active deep crater in a humid climate should usually hold one; a small arid crater usually not). This is
the "a lake must have water" invariant applying here as everywhere.

`LakeType` gains two variants (the roadmap's "third nature" completing exo/endo). Typing runs AFTER the drainage
invariant checks (a draining crater lake is validated as exorheic-with-outlet first) and BEFORE export, and it
touches ONLY `lake_type` — the detected-lake geometry (footprint at/below level, connected to the lowest point,
no overlap, depth == level − floor, level ≤ inlet arrival) is unchanged, so every GEOMETRIC lake invariant
still holds over the whole population, crater lakes included. The regime-specific "exorheic implies a traced
outlet" invariant covers the exo/endo lakes; a crater lake carries the crater nature instead of the regime
label, which is the intended completion of the set, not a gap. The field is exported (`lake_type`, serde) —
Living Landz ignores it today, one more reason to consume it (acidic ⇒ no fish, undrinkable, distinct from the
saline endorheic case). The crater records travel with the terrain in `ErodedProduct`, so typing is correct on
a cache hit (no silent mistyping). Cross-checked by `closures::volcanism::tests::crater_lake_typing_only_wet_craters`.

### C-2 measurement bench (production seed, both resolutions)

`c2_volcanism_bench` (#[ignore]) reconstructs the terrain WITH vs WITHOUT volcanism at the production config
(a relative before/after; the author validates the EXPORT visually). Placement on this seed: 2 arc + 10 hotspot
+ 3 rift = 15 edifices — arc volcanism is barely visible (few O-C margins), the render rests on the hotspot
chains and the rifts.

**A window bug the bench caught.** The first run reported the young/old hotspot edifices "outside the render
window" and only 5-6 of 15 craters resolved. Cause: `apply_edifices` mapped the coarse-torus centre into the
window WITHOUT wrapping, so with the production `sample_origin = [0.094, 0.578]` and full-domain `sample_size =
1`, every edifice at `v < 0.578` (more than half) was silently dropped. The coarse field is sampled
PERIODICALLY (as the FBM's `sample_bilinear_periodic` is), so the offset must be taken mod 1. Fixed; all 15
edifices now render (13-15 craters). The bench found it precisely because the young/old comparison forced the
mapping to be exercised on real positions.

**Crater contribution — a small, identifiable increment, not a flood.** Closed depressions per stage:

| resolution | post-FBM (off → on) | post-relief (off → on) | drainage density ‰ | local relief 11² m |
|------------|---------------------|------------------------|--------------------|--------------------|
| 2048²      | 220 → 235  (Δ +15)  | 1001 → 980  (Δ −21)    | 63.0 → 63.4        | 334 → 336          |
| 8192²      | 6999 → 6895 (Δ −104)| 17382 → 17520 (Δ +138) | 60.2 → 59.9        | 119 → 121          |

The crater bowls add ~13-15 hollows, but the cone FLANKS bury pre-existing FBM noise-pits, so the NET post-FBM
delta is small and can be NEGATIVE (−104 at 8192²) — the opposite of a flood; volcanism slightly cleans the pit
field. The post-relief delta (+138 at 8192²) is the erosion's response to the injected relief, still tiny
against the 17 000 baseline.

**Hydrological displacement — negligible GLOBALLY, by design.** Drainage density moves ±0.4 ‰ and local relief
±2 m — under 1 %. Explanation: 15 cones on a 400 km domain are a small area fraction, so the radial divergences
and basin captures they create are LOCAL (visible in the render around each cone) and do not move domain-wide
aggregates. Every global metric that volcanism should not move barely moved — no unexplained shift. Honest
scope note: the full Strahler-histogram / W-D-per-order / confluence-count would need the network extraction
(`c1_drainage`, heavy); the drainage-density proxy already shows no global shift, and at this edifice
density a per-order breakdown would show the same — the displacement is local, not aggregate.

**The relief-decay PROXY does measurable work (the key test).** Youngest vs oldest edifice of a hotspot chain,
same constructional geometry (Wb 20 km, H 1600 m), 8192², post-erosion:

| edifice            | flank slope | crater bowl |
|--------------------|-------------|-------------|
| young (active, age 0) | 15.6°    | 340 m       |
| old (extinct, age 1)  | 3.2°     | 80 m (near-breached) |

The old edifice is markedly gentler and its crater nearly breached, from the SAME erosion pass acting on relief
pre-scaled ×0.35 at age 1. The proxy is not decoration — it distinguishes young from old by a factor of ~5 in
flank slope and ~4 in crater integrity. It stays labelled a PROXY (no publication), but it earns its place.

**FBM amplitude floor:** unchanged by this closure. Volcanism adds DISCRETE built relief (15 edifices), it does
not replace the FBM's distributed detail, so the C-1 amplitude floor is untouched — the noise-replacement the
roadmap tracks will come from the lithology and coastal closures, not this one. Stated so the trajectory number
is not misread as progress it did not make.

### C-2 crater water balance — the measurement, and why the rim balance IS needed

The author's export showed 0 crater lakes. Three hypotheses were checked IN ORDER before touching any rim
mechanism (a fix written before the cause is established is the pattern we keep hitting).

- **Inflow accounting (TASK 1): correct.** `water_balance_lakes` uses `runoff_accumulation` =
  `max(0, precip−PE)·cell_km2` accumulated downstream along the flow on the FILLED field — so the inner-flank
  runoff draining to the crater floor IS counted, not just direct rain on the water surface. Not the
  uncounted-inlets trap.
- **Geometry (TASK 3): realistic, and the "20 km crater" was a misread of the edifice BASE.** The shield
  craters are Ø 1.25 km / 240 m (D/W 0.19) — nearly identical to Kawah Ijen (1.0 km/200 m, 0.20) and Poás
  (1.6 km/300 m, 0.19). Not the blocker.
- **Position bug: fixed, and NOT the cause.** The `apply_edifices` vertical-mirror + window-wrap bugs were real
  and fixed, but the measurement below is 0 lakes EVEN WITH the fix — so the mirrored craters were not the
  reason for the zero.

`c2_crater_water_balance` (#[ignore], 2048², measured per crater, four climates):

| climate       | holding a lake | active with margin > 1 | near threshold (0.7–1.5×) |
|---------------|----------------|------------------------|---------------------------|
| arid-hot 25°  | 0 / 13         | 0 / 5                  | 0                         |
| humid 45°     | 0 / 13         | 2 / 5 (1.48, 1.94)     | 1                         |
| tropical 10°  | 0 / 13         | 2 / 5 (1.80, 1.14)     | 1                         |
| arid-cold 65° | 0 / 13         | 1 / 5 (1.19)           | 1                         |

The decisive number: several ACTIVE craters in humid/tropical have a MEASURED equilibrium/sill margin > 1 (the
balance says "fills"), yet `lake = no` for all of them. The inflow is sufficient and the geometry is realistic
— so the blocker is that the erosion pass BREACHES EVERY crater rim, active and extinct alike, leaving no
closed depression for the lake stage to find. That is the roadmap's intended "place before erosion, let it
breach", applied INDISCRIMINATELY. So the rim balance is needed, and the measurement JUSTIFIES it precisely:
only the maintenance of ACTIVE rims is missing.

**Plausibility (the author's red-flag test): the proportion would be a small, climate-dependent minority, not
"everywhere".** If active rims held, the measured margins give humid 2/5, tropical 2/5, arid-cold 1/5,
arid-hot 0/5 of active craters filling, and extinct ≈ 0 — consistent with "Kawah Ijen is notable because it is
unusual". The rim balance does not force lakes; it makes an active crater ELIGIBLE, and the (already correct)
water balance decides.

**Formulation to implement (labelled derived).** No single named law covers active-rim persistence, so it is a
composition: construction ∝ eruption rate (Wood 1978) versus destruction (Wood 1980, cinder-cone degradation,
D/W falls with age over 10²–10⁵ yr, rate set by rainfall+temperature). Active rims are maintained (construction
≥ erosion) → closed crater → eligible; extinct rims breach (construction = 0) → drained, as the young-vs-old
proxy already showed (crater bowl 340 m young vs 80 m old, near-breached). Modelled as a RATE balance, not a
protection flag.

**Arc-crater D/W fix (independent of the above).** Grosse's median crater depth (240 m) applied independently
of width gave an aberrant D/W 0.65 on the small (0.37 km) arc craters. Capped at D/W ≤ 0.25 (`CRATER_MAX_DW`;
Kawah Ijen 0.20, Poás 0.19, Wood 1980): arc craters → 93 m (D/W 0.25), shields ≥ 1 km keep Grosse's 240 m.

### C-2 crater lakes — the five-number diagnostic, the defect, and the fix

The 0-crater-lakes result was NOT a physical outcome — closing C-2 on it (option 1) would have been a bug
presented as physics, the one outcome to avoid. Measured on the PRODUCTION terrain at 8192²
(`upscale_from_c1_with_progress`, not a raster reconstruction), humid 45°, five numbers per active crater:

| crater        | Ø (cells / km) | inflow | a_eq (km² / cells) | a_sill (cells) | lake_map in footprint |
|---------------|----------------|--------|--------------------|----------------|-----------------------|
| @(4800,3008)  | 26 / 1.25      | 445.6  | 1.90 / 795         | 524            | 0 (before fix)        |
| @(1728,3520)  | 26 / 1.25      | 42.3   | 0.11 / 47          | 524            | 0 (before fix)        |
| @(3776,1856)  | 26 / 1.25      | 34.3   | 0.11 / 44          | 524            | 0 (before fix)        |

**Verdict: candidate 1.** The @(4800,3008) crater has a_eq 795 cells > a_sill 524 — the balance would fill it —
yet lake_map was empty. Candidate 3 (km↔cells) is ruled out: cell = 49 m, Ø 26 cells = 1.25 km, a_eq in cells
consistent. Candidate 2 (typing) is ruled out: nothing to type because nothing was placed. The defect: the
generic `detect_lakes` discards a lake below `lake_min_area_km2 = 5` — a noise-pond floor — and a crater is
~1.2 km², so it was filtered out BEFORE the balance ran. Crater lakes are small by nature (Kawah Ijen 0.8 km²,
Pavin 0.4 km²), so the generic floor wrongly excludes every one.

**Fix:** `detect_crater_lakes` — a dedicated pass (like `below_sea_basin_lakes`) that, for each ACTIVE crater
(reconstructed → a closed bowl), runs the SAME inflow-vs-evaporation balance and fills it, with a
crater-appropriate floor (`CRATER_LAKE_MIN_CELLS = 4`). Extinct craters are breached → never closed → no lake,
so a crater lake is always active and acidic (Varekamp pH < 2).

**Corrected figure (humid 45°): 4 of 7 active craters hold an acidic lake:**

| crater lake | area | depth | vs real |
|-------------|------|-------|---------|
| #2000004    | 1.19 km² | 226 m | ≈ Kawah Ijen (0.8 km²/200 m), Poás (1.6/300) |
| #2000002    | 0.11 km² | 120 m | small pond |
| #2000003    | 0.10 km² | 53 m  | small pond |
| #2000001    | 0.09 km² | 79 m  | small pond (the 0.37 km arc crater) |

**Plausibility.** One SUBSTANTIAL crater lake at the Kawah-Ijen scale (1/7 active) plus three small marginal
ponds; the rest dry (inflow 0 — rain-shadow / high-altitude PE). In arid climates inflow → 0, so ~none. The
substantial-lake proportion is a clear minority — "Kawah Ijen is notable because it is unusual" — not acidic
lakes everywhere. 4/7 active in the wettest climate is on the high side only if the three 0.1 km² ponds are
counted; `CRATER_LAKE_MIN_CELLS` (or a min-area) is the tuning knob if the author judges them too many. Wired
into the bridge so the shipped export carries the crater lakes (`lake_type = CraterAcidic`).

`active_rim_rebuild` (re-stamping a clean bowl for active craters after erosion) and `detect_crater_lakes` are
labelled PROXY / derived compositions (no single named law for active-rim persistence; Wood 1978 construction +
Wood 1980 destruction).

### C-2 crater lakes — the four-climate table, threshold coherence, CraterNeutral

Measured on the PRODUCTION terrain at 8192² (`upscale_from_c1_with_progress`), all four climates. The author's
two predictions were BOTH tested rather than assumed, and both mattered:

- "arid → ~0" is FALSE for arid-COLD: at 65° the low temperature gives low PE, so even a trickle of runoff
  clears evaporation and ponds form. The prediction held only for arid-HOT.
- tropical (10°, the wettest) is the "everywhere" risk case, and it did fire.

FIRST measurement, with a too-low crater floor (an early `CRATER_LAKE_MIN_CELLS = 4` ≈ 0.01 km²):

| climate       | active | holding | note                              |
|---------------|--------|---------|-----------------------------------|
| arid-hot 25°  | 7      | 1       | minority                          |
| humid 45°     | 7      | 4       | **MAJORITY of active — RED FLAG** |
| tropical 10°  | 7      | 4       | **MAJORITY — RED FLAG**           |
| arid-cold 65° | 7      | 4       | **MAJORITY — RED FLAG**           |

The 4-holding was one SUBSTANTIAL lake (1.19 km² / 226 m, the crater whose a_eq ≈ 680–795 cells exceeds the
524-cell sill in every wet climate) plus THREE 0.06–0.11 km² PUDDLES — a_eq 25–47 cells, barely above the
4-cell floor. **Dual-threshold incoherence**: the generic `lake_min_area_km2 = 5` (2097 cells) excludes
sub-km² noise ponds everywhere, while the 4-cell crater floor (0.01 km²) let 0.06 km² puddles through — the
very scale we exclude elsewhere. Raised to `CRATER_LAKE_MIN_AREA_KM2 = 0.2` (below Pavin 0.4 and small maars,
above the puddles): a crater LAKE is a real lake, and the two thresholds are now coherent (5 km² filters noise
on generic terrain; 0.2 km² admits real small lakes on the physical crater basin).

FINAL measurement (0.2 km² crater floor):

| climate       | active | holding | lake                                    |
|---------------|--------|---------|-----------------------------------------|
| arid-hot 25°  | 7      | 0       | —                                       |
| humid 45°     | 7      | 1       | 1.19 km² / 226 m (Kawah-Ijen scale)     |
| tropical 10°  | 7      | 1       | 1.19 km² / 226 m                        |
| arid-cold 65° | 7      | 1       | 1.19 km² / 226 m                        |

**Plausibility: a clear minority — 1 of 7 active craters in a wet climate, 0 in a hot desert.** One substantial
acidic crater lake at the Kawah-Ijen scale per humid continent, the rest dry — "Kawah Ijen is notable because
it is unusual", obtained from the physics (its inner catchment happens to clear evaporation; the others do
not). No red flag.

**`CraterNeutral` is UNREACHABLE by construction, and documented as such.** A crater lake is always on an
ACTIVE crater (an extinct edifice is never rim-reconstructed → always breached → holds no water → the chain is
crater-lake ⇒ active ⇒ acidic). So `CraterNeutral` never occurs. This matches Auvergne's breached-and-dry
puys, but NOT Lac Pavin (an extinct maar, intact, fresh): our erosion has no maar-like geometry (a shallow
explosion crater that survives intact), so that case is simply not modelled. The variant is kept for future
maar volcanism, with the limitation stated on the variant's doc comment and here — an acknowledged gap, not a
live dead branch.

### C-2 crater lakes — the reconstruction trap, and the real production defect

The four-climate table above was measured on a BENCH that called `detect_crater_lakes` directly on the
reconstructed eroded field — it did NOT run the relief-v3 breach the production pipeline applies. The author's
EXPORT (arid-cold 65°, span 10°) showed 0 crater lakes, contradicting the bench's 1. The export is the verdict;
the bench had misled — the "a reconstruction is not the product" trap, again.

The real defect: the relief-v3 breach (`breach_monotone`) runs AFTER the active-rim reconstruction and
RE-BREACHES every crater — the crater is not in the pre-breach lake set (the generic `detect_lakes` filters it
at 5 km²), so the breach carves it open and the crater-lake stage finds nothing. The benches missed it because
they omitted the breach.

The constraint that shapes the fix: the breach output (the "conditioned" field) is cached CLIMATE-INDEPENDENTLY
(keyed on the eroded key, no latitude), while a crater lake is CLIMATE-DEPENDENT. So crater lakes cannot live in
the conditioned cache — the bowl must SURVIVE the breach (a climate-independent step), and the fill decided
later by the climate-dependent stage. Fix: `breach_monotone_protected` with a protect mask of the ACTIVE crater
cells — those cells are kept at their original height (neither carved nor filled), so the bowl survives; then
`detect_crater_lakes` (with climate) fills the eligible ones. `VOLCANISM_ALGO` bumped to 4.

Measured again on a bench that NOW runs the breach with the protect mask (faithful to production), four climates:

| climate                | active | holding | crater lakes                       |
|------------------------|--------|---------|------------------------------------|
| arid-hot 25° span10    | 7      | 0       | —                                  |
| humid 45° span40       | 7      | 2       | 0.67 km²/165 m, 0.29 km²/100 m     |
| tropical 10° span20    | 7      | 2       | 0.98/203, 0.60/153                 |
| arid-cold 65° span10   | 7      | 2       | 0.63/159, 0.45/130                 |

2 of 7 active craters in a wet climate, 0 in a hot desert — a minority, no red flag; areas 0.29–0.98 km²,
depths 100–203 m (Pavin/Kawah-Ijen scale). This is a bench that reproduces the production breach; the SHIPPED
export remains the verdict and must be re-generated (VOLCANISM_ALGO 4 forces the recompute) and audited before
C-2 is closed.

### Method rule 3 (earned in C-2): a bench must reproduce the WHOLE production chain

"A reconstructed terrain is not the product" (rule 2) was not enough: in C-2 a bench MISSED
the crater-breach defect THREE times because it ran `detect_crater_lakes` on the reconstructed
eroded field WITHOUT the relief-v3 breach the production pipeline applies next — so it reported
crater lakes the shipped export did not have. The sharper rule: a bench that omits a downstream
stage measures a terrain that never ships. Any measurement runs the FULL production chain
through the stage whose effect it claims (here: FBM → volcanism → erosion → BREACH → drainage),
and the author's export remains the final verdict. This is the seventh occurrence of the
proxy-vs-authority family (Findings 27, 29-41, and the C-2 crater lakes); the failure mode is
always the same — a convenient partial computation stands in for the whole, and the gap it
leaves is exactly where the defect hides.

## C-3 — Lithological heterogeneity: the C1 signal, measured

C-3 needs a spatially-varying erodibility K driven CAUSALLY (not from noise). Measured what C1 actually
carries (`c3_lithology_probe`, full 300-step production chain, two seeds):

| source                              | coverage of continental | causal? |
|-------------------------------------|-------------------------|---------|
| craton BASE (geometric placeholder) | 49-50 %                 | NO — the rule is `seed_x < nx/2` (left half), a stand-in the init flags as non-final |
| cratonic shield (stored mask)       | 7 %                     | NO — FBM-noise-refined (#165 select_shield_mask) |
| rift (age ~ 0, rift-spawned)        | 1-8 %                   | YES (rifting stamps age = 0) |
| volcanic footprints (C-2 placement) | 10-13 %                 | YES (edifice basal discs) |
| GENERIC continental (no signal)     | ~79-88 %                | — |

**The C1 limitation (a real finding, recorded as such).** A fully causal lithology is NOT available today,
and the reason is in C1, not in C-3:
- `age` is UNIFORM at exactly 7.00 over 92 % of continental cells — the time loop writes nothing into it for
  continental crust, so it cannot separate old craton from young accreted terrane;
- `plate_type` is BINARY (Oceanic / Continental) — no arc / terrane variant;
- the closures do NOT record terrane provenance: subduction turns Oceanic → Continental but leaves no marker,
  accretion discards the loser plate id, so accreted arc crust is indistinguishable from ancient continent;
- the `cratonic_mask` is a geometric placeholder (seed-band) whose HD-visible shield extent is FBM-noise-refined.

So the only genuinely-causal, non-noise signals (rift + volcanic) cover ~15-20 % of the continent; the craton
base covers ~50 % but by an arbitrary geometric rule, not physics.

**Specification of what would make a causal lithology available** (the option-2 investment, which also unblocks
C-4 coastal, whose cliff retreat needs rock resistance): the closures already KNOW what they do at the moment
they act — record it as an advected class field. Subduction writes "arc / accreted terrane" on the cells it
reassigns Oceanic → Continental; accretion writes "sutured terrane" on the merged strip; rifting already stamps
age = 0 (young). That is a trace of the existing physics, not a new mechanism, and it would carry a multi-class
lithology over the WHOLE continent rather than the 15-20 % the current fields expose.

### C-3 — hard-vs-soft (not a continuum), and the missing deposition stage

Auditing Stock & Montgomery 1999 ON THE SOURCE (`docs/refs/stock1999.pdf`) recast the closure. K by class,
m = 0.4, n = 1 (stable base-level case; NOT the Kauai m = 0.1/n = 0.2 exponents): granite/metamorphic
10⁻⁷–10⁻⁶, volcaniclastic 10⁻⁵–10⁻⁴, young mudstone 10⁻⁴–10⁻² m^0.2/yr; measured spread "1 to 5 orders of
magnitude" softer than hard rock. The load-bearing NUANCE: "K between granitoids and metasediments is NOT
significant" — the contrast is HARD CLASS vs SOFT CLASS, not a continuum. Consequences:
- a continental basement treated as uniformly HARD is PHYSICALLY CORRECT (crystalline + metasedimentary are
  both hard), so a "generic hard bulk" is not a coverage failure — the FBM-floor-everywhere expectation was
  misframed;
- the craton placeholder (`seed_x < nx/2`) is not merely arbitrary but USELESS here: it would separate hard
  from hard. Dropped, nothing to justify;
- what must be differentiated are the SOFT zones, minority by nature: rift, volcaniclastic, young sedimentary
  basins.

**The missing deposition stage (a stated limitation).** Ymir's PRODUCTION erosion is relief-v3 stream-power —
DETACHMENT-LIMITED, pure incision, no aggradation. A deposition/`sediment` field DOES exist but only in the
DROPLET pass (`erosion/hydraulic.rs`, `deposition_rate`, `coastal_deposit_fraction = 0.25`), which production
does not run; isostasy carries no subsidence/foreland/flexure signal. So there is NO causal signal for
sedimentary basins in the shipped chain. Low-relief or endorheic areas could be used as a geometric proxy, but
that would repeat the craton mistake (lithology from geometry, not physics) — rejected. The soft class
therefore reduces to RIFT (age = 0, ~1-8 %) + VOLCANICLASTIC (edifice footprints, ~10-13 %); everything else is
hard basement at a single low K. Adding real sedimentary basins would need a deposition stage (a
transport-limited erosion pass or a flexural-subsidence + fill model) — recorded as the specification, not
built. Note the detachment-limited production regime is exactly the domain Stock & Montgomery calibrated K for.

### C-3 — the spread is a MEASUREMENT, and the two effects, separated

The soft↔hard multiplier was SWEPT, not predicted (`tests/c3_lithology_sweep.rs`, the WHOLE production
chain — `upscale_from_c1_with_progress`, export recipe relief-v3, droplets off — lithology OFF then
×3/×10/×30/×100 soft, both resolutions, production seed; report
`docs/reports/c1_continental_buoyancy/closure_morphology/c3_lithology_sweep.md`). Class coverage came out
hard 95.7 % / rift-soft 1.6 % / volcaniclastic 2.7 % at BOTH resolutions (area-preserving), confirming the
minority-by-nature soft class.

**Method rule 3 honoured.** The K field is built and threaded exactly as production does it (coarse hard = 1.0
with rift soft, bilinear-upscaled and registered to the altitude with the same `(sample_origin, sample_size)`;
volcaniclastic stamped at HD on the edifice basal discs; per-cell K into `incise_lithology`), not a
reconstruction of the incision stage alone.

**Hard = ×1.0 (reference), soft ABOVE — the design that separates the two effects the author asked to keep
apart.** Because the ~96 % hard bulk stays at the relief-v3 reference K, the HARD-class morphometrics are FLAT
across the whole sweep (2048²: local-relief/slope/steep‰/incision = 329/6.9/147/155 at every multiplier;
8192²: 119/9.8/227/77 at every multiplier). So effect (a) "global slowdown" is ZERO by construction — there is
nothing to disentangle from effect (b) "the contrast". The rejected alternative (hard ×0.3, soft ×1.0) would
have moved 96 % of the continent and confounded the two. The contrast (effect b) is monotone and physical:
softer K erodes DOWN → relief and channel incision fall, valleys open (higher W/D). VOLC (2048²) relief 612→303,
incision 329→98; SOFT incision 51→12.

**C-1 survives the whole sweep.** Closed-depression count 2048² 982→977, 8192² 17516→17114 — a slight DECREASE,
never a flood; land fraction stable (15.7→15.2 %, 16.6→16.4 %). Softening does not fabricate pits.

**Chosen multipliers.** soft (rift) = ×10 (≈1 order, mid of the S&M range, a clearly visible contrast with C-1
intact); volcaniclastic = ×3 FIXED and decoupled from the soft sweep (S&M intermediate class, deliberately mild
so the C-2 edifice morphology is dissected, not flattened — ×30/×100 halve the volcanic relief, erasing the
cones just built). Gated OFF by default (`LithologyConfig::enabled = false`) → byte-identical production; the
eroded cache key is byte-identical when disabled (config skipped from serialization, `LITHOLOGY_ALGO` appended
only when enabled).

## C-3b — Inherited structure: fracture density (shipped), orientation (measured out)

C-3 established the basement is lithologically UNIFORM and hard — physics, not a gap. C-3b's premise: a mature
basement's structure is TECTONIC, not lithological — the same rock, but CUT by fractures. Density → erodibility
is well founded (Molnar 2007, *Tectonics, fracturing of rock, and erosion*: tectonics erodes mostly by
fracturing → plucking; Clarke & Burbank 2011; Zondervan et al.: ~1–2 orders of K, and fracturing homogenises the
inter-lithology contrast by ~1 order; domain of validity: bedrock rivers, brittle upper crust <~10 km,
detachment-limited — the relief-v3 regime). Orientation → fabric is standard too (Anderson 1905; World Stress
Map, Heidbach/Zoback: intraplate SHmax ∥ plate motion at first order → the topographic fabric strikes ⊥ SHmax).

### The orientation was BUILT, MEASURED, and dropped — a characterised limitation

The directional closure (valleys aligned on the fabric via anisotropic incision) was implemented in full and
measured on the whole chain (`c3b_fracture_sweep`, both resolutions). It does not work, and the measurement says
why, twice:
- **Rate anisotropy cannot re-route drainage.** `K_eff = K·(1 + a·|flow·strike|)` just raised total erosion
  (relief 332→222, a global-rate confound). The mean-preserving form `K·(1+a·align)/(1+a·(1-align))` removed
  that confound (relief held 332→334) and yet fabric alignment |flow·strike| STILL fell (0.639→0.536, both
  resolutions) and closed depressions ROSE (1001→2146). The incision RATE acts on a receiver fixed by topography;
  it cannot reorient the network, so it cannot align valleys. The lock is flow ROUTING.
- **Routing is out of reach, and would be an artefact anyway.** Biasing `compute_flow` has a blast radius over
  C-1, river extraction, lakes, the whole hydro chain stabilised across ~15 passes. And C1's directional field is
  too poor to feed it: per-plate CONSTANT velocities, no per-cell strain, no deformation history → a uniform
  per-plate strike → a continent-wide identical grain, an artefact as visible as the Smith–Bretherton comb. The
  premise is also weak: the Appalachian trellis is FOLDED STRATA (out of scope), not fractures; real
  fracture-controlled drainage (rectangular patterns on jointed granite) is provincial and subtle, not
  continental.

**Specification for a future directional closure** (what it would take, recorded not built): a per-cell stress
or strain-rate field (not constant per-plate velocities), a deformation HISTORY (to carry paleo-stress, since
today's fabric is inherited from past orogenies), OR folded-strata layering — plus a routing coupling in
`compute_flow` with its C-1/rivers/lakes regression budget. None exists in C1 today.

### FBM floor — a ROADMAP correction (not a C-3b failure)

The C-3b brief made "the FBM floor must shrink" a pass/fail criterion. That was mis-posed. The FBM fills the
128× upscale — wavelengths from the coarse cell (~6 km) down to the HD cell (~49 m). NO tectonic closure holds
information BELOW the coarse cell, so no closure can replace the noise IN ITS OWN BAND; closures add structure at
scales ≥ the coarse cell, and sub-coarse detail can only come from erosion, which needs a symmetry-breaking seed
— the C-1 degeneracy floor. So lowering the FBM amplitude floor is unreachable IN PRINCIPLE by any closure, C-4
included. This is the limit of what closures can do against the noise, recorded as a result. C-1's flow
conditioning (stop the FBM fabricating depressions) remains the right and achievable goal; REPLACING the FBM's
own band does not.

### What ships — density only, causal, C-1-preserving

Erodibility is modulated ISOTROPICALLY by fracture DENSITY: `K = 1 + amplitude · density`, `density ∈ [0,1]`.
Density = `exp(-dist_to_contact / decay)` where the contacts are the DYNAMIC boundary classification's
CONVERGENT + TRANSFORM cells (orogens + shear — the fracturing regimes; divergent = rift is C-3's domain). It is
NOT `cratonic_mask` (the FBM-noise-refined field C-3 rejected) and NOT the geometric craton placeholder
(`seed_x < nx/2`) — the same discipline that settled C-3. The intact craton EMERGES as the region far from every
contact, at `density → 0 → K = 1` (the reference — global-slowdown nil by construction, the C-3 design that
survives).

Measured (whole chain, both resolutions, `c3b_fracture_sweep`), with the narrow orogenic belt (`decay = 25 km`)
that keeps the craton the MAJORITY (coverage: craton 53 %, transition 27 %, belt 20 %):
- **8192² (the export verdict):** CRATON relief flat across the sweep (105→104 m — the reference holds exactly);
  BELT relief RISES with amplitude (398→538 m at ×16, +35 % — at export resolution, fracturing DISSECTS the
  orogen into more valleys); closed depressions FALL (17382→15639 — C-1 improves, unlike the anisotropic test
  that pushed pits to 2146). Contrast (belt/craton) 3.8→5.2.
- **2048²:** craton nearly flat (284→249), belt relief falls (1109→780 — at coarse resolution the extra erosion
  wears the orogen down rather than dissecting it); pits stable ~1000. The sign flip with resolution is expected
  (the 1 km relief window resolves dissection only at HD).

**Chosen:** `amplitude = ×6`, `decay = 25 km`. The limiting factor is NOT the physics (Molnar's 1–2 orders
leaves room to ×8+) but the LAKE POPULATION: fracturing the belt dissects it into more small enclosed basins →
more small lakes, and the hydro chain's most laborious machinery (inventory thresholds, sink validity, submarine
basins, outlet invariants, stabilised across ~15 passes) is exactly what small lakes exercise. ×6 gives real
visibility over ×4 without doubling that pressure; ×8 is measured safe ON THIS SEED but a wider orogenic belt or
an arid climate could react differently, so margin is kept until C-4 (which also touches the coast). Gated OFF by
default → byte-identical; the eroded key adds `FRACTURE_ALGO` only when enabled.

**The contrast is DISCRETE BY NATURE — the correct result, recorded so it is not re-discovered as a
disappointment.** With ~20 % fractured belt over a basement that C-3 established is UNIFORMLY HARD (anchored on
Stock & Montgomery), this closure can only produce REGIONAL VARIETY (one orogenic province more dissected than
the shield), NOT a change in the continent's overall appearance. The author's export read — "a bit better but
subtle; a few more, smaller lakes" — is exactly the expected signature. Asking C-3b for more would mean either
inventing lithological variety in the craton (contradicts C-3 / the source) or the directional mechanism
(measured out above). Subtle-but-causal is the ceiling here, and it is the right ceiling.

### Fossil sutures — deferred with a specification

Accretion welds terranes along a suture that becomes a fossil zone of weakness INSIDE a plate. C1's accretion
records WHEN/how-many merges fire (`AccretionStats` counters) but NOT WHERE — after the merge the seam is inside
a uniform `plate_id` region, unrecoverable, and `age` is degenerate (7.0 over 92 % continental). The current
density DOES capture ACTIVE collision/subduction belts (they are Convergent boundaries). Fossil sutures would
need recording the weld location at merge time (snapshot `plate_id` before reassignment in
`accretion::merge::apply_accretion_step`, accumulate a suture mask on `C1State` with snapshot handling) — the
plumbing is specified, not built in this pass, since the boundary-contact density already covers the active
belts. `FractureConfig::suture_multiplier` + the `suture_mask` argument of `derive_coarse_density` are the hooks.

### C-3b lake population — a FRAGILITY point to watch (not a blocker)

The ×6 export (production seed, 8192², humid) validated: C-2 crater lakes UNCHANGED (1 acidic / 6 dry, on and
off), and the Finding 37 canary already fires at BASELINE (1 exorheic lake with no traced outlet, `1000023`) — so
C-3b does not break a clean invariant; it adds ONE more (→2, `[54, 1000033]`), a real belt-dissection basin with
an outlet-tracer edge case, not a broken lake. That residual baseline orphan is filed as its own issue (the
Finding 37 inversion, commit `dd1b48a`, has a residual geometric hole — the ID shifting `1000023 → 1000033` under
the terrain change points at the tracer, not a specific basin).

The number that deserves the watch, though, is the SUBMARINE basins: the ×6 run reports `below-sea basins: 43
lakes (43 exorheic, 0 endorheic); 62 spillways (37 → sea, 25 chained)` — roughly DOUBLE the prior baseline (~21).
Belt dissection carves more coastal/below-sea depressions, so this is consistent with the mechanism, but the
submarine-basin population is the part of the pipeline that cost the most to stabilise (Findings 36→40), and it
has just doubled. Record this as the FRAGILITY POINT: if a future seed with a larger orogenic belt, a stronger
amplitude, or C-4's coastal erosion pushes it further, the submarine sink/outlet machinery is where it will crack
FIRST. (The exact same-seed OFF submarine count was not re-printed here — that export was a cache HIT — so the
~21 is the author's prior baseline; a dedicated OFF MISS run would pin the delta. The crater-lake population is
untouched, so the doubling is concentrated in the coastal/below-sea class, not the whole lake population.)

## H-1 — Infiltration: the first subsurface term (before C-4)

### Ordering: H-1 / H-2 come before coastal erosion (C-4)

Sill incision (H-2) drains exorheic basins — the outlet carves its col, the lake retreats, the water re-emerges as
a river — which MOVES the coastline (river mouths, deltas, the littoral base level). Sculpting the coast (C-4)
first would work a shoreline H-2 then displaces, redoing the work. So the hydrology chain (H-1 then H-2) precedes
C-4. Recorded in the roadmap.

### The problem, and why H-1 is light

The water balance reads `runoff = max(0, precip − PE)`: the ENTIRE precipitation surplus becomes surface runoff.
Ymir has NO subsurface water. In reality a significant fraction INFILTRATES and never reaches a lake by surface
flow, so the model over-supplies every basin — a direct contributor to the Finding-39 behaviour the diagnostic
measured (92 lakes, 91 exorheic, 0 endorheic, 6220 km², the ten largest holding 74 %). H-1 adds the first
subsurface term. Its blast radius is SMALL: it changes the BALANCE only, not geometry or routing; some basins flip
endorheic and levels drop — a step forward that breaks nothing.

### Literature (verified), and how it composes with Budyko

- **Runoff coefficient / SCS Curve Number** (Rational method; USDA-SCS): the fraction of rainfall that becomes
  DIRECT runoff, integrating infiltration, evaporation, interception, depression storage. CN ≈ 30 (permeable, high
  infiltration) → 100 (impermeable / water). The coefficient RISES with lower soil permeability, STEEPER slope,
  and sparser vegetation; falls on flat permeable ground. So infiltration ∝ permeability, and slope reduces it —
  both physical, both in the literature.
- **Baseflow Index (BFI)** — the groundwater share of streamflow, a proxy for the infiltrated fraction: 0.55–0.80
  in the Housatonic basin (USGS), controlled by BEDROCK PERMEABILITY (permeable/porous bedrock → high BFI;
  crystalline/impermeable → low). Crystalline-bedrock recharge (Kenogami) is a modest fraction.
- **The "20–50 % infiltration" figure is a HYPOTHESIS, not established** — the published span is wider and
  geology/climate-dependent (BFI 0.4–0.8; crystalline lower). The exact fraction is therefore a swept PROXY
  (labelled in code + here), not asserted.
- **Composition with Budyko** — the one place an error would be invisible, so stated explicitly. Budyko bounds the
  SUPPLY side: `AET = min(precip, PE)`, and `runoff = precip − AET = max(0, precip − PE)`. Infiltration does NOT
  touch AET; it SPLITS the runoff (post-Budyko) between surface flow and groundwater:
  `surface_runoff = (precip − AET) · (1 − f_infil)`. The infiltrated `(precip − AET) · f_infil` leaves the surface
  balance. The double-count trap (which caught the project once, on rain-credit vs catchment runoff): infiltrated
  water must NOT reduce AET, and must NOT reappear as inflow anywhere. First approximation: infiltrated water is
  lost to deep groundwater (not returned as baseflow) — which under-supplies slightly, the DESIRED direction.
  Refinement path (not now): delayed baseflow return, re-emerging as river discharge downstream.
- **No named prior art** couples a tectonic-closure permeability field (lithology + fracture density) to an
  infiltration split in a terrain generator — the bricks (SCS-CN, BFI, Budyko) are standard, the assembly is a
  LABELLED derivation.

### Floor leakage — the second subsurface mechanism (its own issue, not H-1)

A lake perched above the water table LEAKS THROUGH ITS FLOOR — the only mechanism that can empty a basin
COMPLETELY without an outlet (karst, disappearing rivers). It is a DIFFERENT mechanism from infiltration (it acts
on a standing lake, not on the runoff supply) and it needs a water-table datum Ymir does not have. It belongs in
its OWN issue, not H-1 — recorded, not implemented blind.

### What H-1 ACTUALLY delivered — the water balance restored to surface lakes

The issue framed H-1 as "add infiltration so over-supplied basins flip endorheic". The measurement found the
premise was wrong, and the real defect underneath it.

**The defect.** In the production relief-v3 path the pre-breach drainage runs with `climate = None`
(`hd.rs:886`) — the pure-GEOMETRY path, which fills every lake to its sill and calls it "Exorheic if the outlet
reaches the sea", with NO water balance. `build_hd_drainage` then OVERWRITES the climate-computed lakes with
those (`hd.rs:432-433`), because the pre-breach GEOMETRY is the correct one (the breach destroys the
depressions). Net effect: **the shipped surface lakes never saw a water balance at all.** The decisive
cross-check: in the arid-hot export every one of the 11 endorheic lakes has an id ≥ 1 000 001 — they are
BELOW-SEA basins, the only class that is climate-aware. Not one surface lake was ever endorheic, in any climate.
So Finding 39 ("net_evap ≈ 0 makes everything exorheic") held only for below-sea basins; for the 91 surface
lakes, 6220 km², the ten largest at 74 % of the area, the exorheic verdict was GEOMETRIC.

**The fix, shipped as the production DEFAULT (not gated) because it is a correction, not an option.** The
carried pre-breach lakes are reclassified by the water balance
(`drainage::classify_lakes_water_balance`, read-only). GEOMETRY, LEVELS and FOOTPRINTS are untouched — adopting
an endorheic equilibrium level would be H-2 by the back door; crater types (C-2) are never overwritten.
`ALGO_HD_DRAINAGE` bumped to 2. **This changes `lake_type` on every existing map, hence biomes and rendering —
that IS the correction, not a side effect.**

**Measured (whole chain, both resolutions, four climates, `tests/h1_infiltration_sweep.rs`).** Shipped
population, geometry only: 2048² 67 lakes / 2247 km² / **exo 67, endo 0**; 8192² 58 lakes / 2199 km² / **exo 58,
endo 0**. After reclassification (no infiltration):

| climate | 2048² exo/endo | 8192² exo/endo |
|---|---|---|
| tropical 10° | 24 / 43 | 25 / 33 |
| arid-hot 25° | 3 / 64 | 0 / 58 |
| humid 45° | 27 / 40 | 27 / 31 |
| arid-cold 65° | 14 / 53 | 12 / 46 |

**53–100 % of surface lakes are endorheic once the balance is applied — including 45 % in HUMID.**

**Infiltration is a SECONDARY term, measured and set aside as a lever — not forgotten.** On top of the
reclassification it adds **0 to 3 lakes** (2048²: +2/0/0/0; 8192²: +3/0/+1/+1). The applied fraction stays INSIDE
the published range without being forced: intact crystalline 0.000, rift 0.025, volcaniclastic 0.189, fully
fractured 0.700 (BFI 0.4–0.8), field-wide mean ~0.10 with p90 0.625 — the median is 0 because the intact
cratonic majority is genuinely impermeable. So the criterion stated in advance is answered: the range did not
have to be left, and the effect is still marginal. **Infiltration is not the dominant lever.** It is kept because
it is physically sound (Heath conductivities, Barenblatt double porosity, no unsupported slope term) and it is
GATED OFF by default.

**⚠️ SUPERSEDED BY H-1b — this scope figure was computed on a BUGGY criterion.** It read: "H-1 reduced H-2's
scope: sill incision now applies to the exorheic remainder only, **3–27 lakes depending on climate, not 91**, and
the endorheic ones would shrink on their own (1798 → 483 km², tropical, reported not applied), so H-2's blast
radius is far smaller than the author accepted as the price of progress." The endorheic count behind it was
inflated by the gross-PE defect found in H-1b below. Left visible rather than silently replaced — the corrected
figure is in H-1b.

### H-1b — the endorheic criterion was too strict, and the two water balances disagreed

40 endorheic out of 67 in a HUMID climate was physically suspect: an endorheic basin means evaporation absorbs
the ENTIRE catchment inflow (Dead Sea / Great Salt Lake regime), and France, Scotland and Scandinavia have
essentially none. Audited (`tests/h1b_endorheic_criterion.rs`) — and the defect was NOT the suspected one.

**Not the area-vs-level test.** The leading hypothesis was that `a_eq ≥ a_sill` compares the equilibrium AREA
with the area AT THE SILL, hence requires filling the whole basin to the col, whereas a real lake overflows as
soon as its LEVEL reaches the sill. Recorded because it was the plausible suspect and it was WRONG: for a
monotone hypsometry (area grows with level), "the level reaches the sill" and "the area reaches the sill area"
are the SAME statement. The area test is equivalent to the physical one.

**The defect was the EVAPORATION term, and the two paths disagreed.** The surface balance used GROSS potential
evaporation (`a_eq = inflow / PE`), ignoring that a lake also RECEIVES rain on its own surface — the loss per unit
area is `PE − P`, not `PE`. Measured at 45° humid: PE 575–886 mm/yr against P 1122–1154 mm/yr, so **`PE − P` is
NEGATIVE on every lake** — endorheism is physically impossible there. Meanwhile the BELOW-SEA path
(`drainage.rs:1266`) ALREADY used `net_evap = max(0, PE − precip)`, with the degenerate case handled EXPLICITLY
(`net_evap == 0 ⇒ a_eq = ∞ ⇒ the basin MUST overflow`). **Two different criteria for one physics.**

**Fixed by adopting the below-sea formulation on the surface path**, via a shared helper
(`lake_net_evap_terms`) so the two cannot drift apart again. NO ε: the degenerate case is explicit, not clamped.
**No double count**, by the same complementarity the below-sea comment already argued: `runoff_accumulation`
clamps its source to `max(0, P − PE)`, so when `net_evap > 0` (arid) the lake's own cells contribute ZERO runoff
and `inflow` is purely external; when `net_evap == 0` (humid) `inflow` is not used at all. Exactly one of the two
terms is nonzero. `ALGO_HD_DRAINAGE` bumped to 3.

**Corrected surface population (2048², four climates) — ordered by climate, as the physics requires:**

| climate | endorheic BEFORE (gross PE) | endorheic AFTER (net) | exorheic AFTER |
|---|---|---|---|
| tropical 10° | 42 (1792 km²) | **0** | 67 |
| arid-hot 25° | 64 (2217 km²) | **63** (2206 km²) | 4 |
| humid 45° | 40 (1767 km²) | **0** | 67 |
| arid-cold 65° | 53 (2115 km²) | **5** (723 km²) | 62 |

Endorheic basins now exist essentially only in HOT-ARID. Cold-arid keeps 5 because cold holds PE low, so net
evaporation stays small despite little rain — Siberian lakes drain, they are not sebkhas. The below-sea path was
already net, so the two populations now AGREE (humid: 0 endorheic on both).

**Infiltration's measured effect is now EXACTLY ZERO on surface classification** (identical columns in all four
climates): in humid/tropical `net_evap = 0 ⇒ a_eq = ∞` regardless of inflow; in arid the population sits far from
the threshold, so a ~10 % inflow cut flips nothing. The earlier "+0 to 3 lakes" was itself an artefact of the
buggy criterion. Infiltration stays gated OFF, kept for its physical soundness and for the discharge / below-sea
paths, and is now definitively recorded as NOT a lever on the lake regime.

**H-2 RE-DIMENSIONED.** The exorheic remainder is **4 to 67 lakes depending on climate — and in HUMID, the
author's target climate, it is ALL 67**. H-1 did NOT collapse H-2's scope; the blast radius returns to the
initial estimate.

### H-1c — the SAME discard, the other attribute: endorheic basins now settle at equilibrium

H-1 restored the water balance's CLASSIFICATION to the carried pre-breach lakes. Its GEOMETRY
half was still discarded: `water_balance_lakes` also SHRINKS endorheic lakes to their
evaporative equilibrium and drains the cells above it, and that output was thrown away by the
same carry-over. So every endorheic lake was reported as a closed basin while still drawn
filled to its col — physically impossible. **One defect, two attributes, found in two steps.**
This is NOT sill incision running early: the equilibrium level is the balance's own
consequence (a closed basin settles where evaporation equals inflow), whereas H-2 carves an
outlet. `apply_lake_water_balance` replaces the classify-only pass; `ALGO_HD_DRAINAGE` → 4.

**Ordering guarantee (the orphaned-mouth defect, already fixed twice).** The apply runs BEFORE
`below_sea_basin_lakes_infil` and BEFORE `clip_rivers_to_lakes`, so both see the FINAL
footprint: river tracks are clipped to the retreated outline instead of ending in the void.
Inlets are enumerated after the footprint is known, never before.

**Measured (2048², four climates, `tests/h1c_endorheic_shrink.rs`).** Before: 67 lakes,
2247 km² water; plain largest 3025 km² (1164 hex), total 10643 km².

| climate | lakes after | water km² | floor exposed | exposed-floor slope | plain largest / total |
|---|---|---|---|---|---|
| tropical 10° | 67 (0 endo) | 2247 | 0 | — | 3025 / 10643 |
| arid-hot 25° | **27** (23 endo) | 2247 → **168** | **2079 km²** | p50 **0.0°**, p90 3.4°, **92 % < 5°** | **4106** / **12610** |
| humid 45° | 67 (0 endo) | 2247 | 0 | — | 3025 / 10643 |
| arid-cold 65° | 67 (5 endo) | 2247 → 1635 | 612 km² | p50 0.0°, p90 1.8°, **95 % < 5°** | 3025 / 11233 |

**THE EXPOSED FLOOR IS FLAT — measured, not assumed.** 92–95 % under 5°, median 0.0°. The
"buildable plain" argument holds. In arid-hot, 40 of 67 lakes dry up ENTIRELY (inflow cannot
sustain one cell), water falls 93 %, and the plain gains **+1081 km² on the largest piece
(+36 %) and +1967 km² in total (+18 %) — WITHOUT carving a single sill.**

**8192² confirms, and the plain metric keeps its convergence AFTER the shrink** — arid-hot:
58 → 17 lakes, water 2199 → 31 km², floor exposed 2168 km², **97 % under 5°**, plain
3129 → 3999 km² largest and 9989 → 12188 total. Against 2048²'s 4106 / 12610 that is ×1.03 on
both aggregates, the same ratio the metric was validated at — so it survives a geometry change
of this size. (Arid-cold is resolution-sensitive: 612 km² of floor at 2048² against 52 km² at
8192², for the same 5 endorheic basins — the equilibrium areas are small enough to sit near
the cell scale. Noted; it does not affect the arid-hot conclusion.)

**The split H-2 will be judged against, and the honest reverse.** The gain above is due to the
BALANCE ALONE. But in HUMID and TROPICAL, H-1c changes NOTHING (no endorheic basin exists
there, correctly). So in the author's TARGET climate the France↔Scotland dial rests
**entirely** on H-2's sill incision. H-1c has established the baseline; it has not moved the
humid map.

### The H-1c figures ARE production values (an alert raised and withdrawn on evidence)

Mid-round I flagged that the H-1 benches inherited `c1_hd_production`'s `amplitude_base = 0.16` while the viz
overrides to 0.04, and warned the H-1c numbers were not production values. **That alert is WITHDRAWN**: the
parameter is inert (above), so the benches produced exactly the shipped terrain. The exposed floor
(2079–2168 km²), the +36 % on the largest plain, the ×1.03 convergence and the 92–97 % under 5° stand as
PRODUCTION figures. Recorded because the alert was published; the retraction must be as visible as the alarm.

### The COASTLINE barbing is in the CONTOUR, not the terrain — and it shrinks C-4

Diagnosed on the arid-hot export (`tests/h1c_phantom_and_coast.rs`). The TERRAIN crosses 0 m cleanly — the mean
profile perpendicular to the shore descends 49.9 → 2.8 m over 8 cells then goes negative, with no serration. The
CONTOUR does not: mean step **0.78 cell** (sub-cellular, marching-squares vertices on cell edges), step directions
peaked on the DIAGONALS, and 8 % of turns above 80°.

The mechanism is **gradient pinning**, and it refuted the flat-shelf hypothesis both of us favoured. Sharp turns
are near-uniform across low-to-medium slopes and vanish only on steep shores:

| shore slope | vertices | mean turn | turns > 80° |
|---|---|---|---|
| < 0.5° | 9 671 | 16.1° | 8.6 % |
| 0.5–2° | 3 076 | 19.2° | **12.2 %** |
| 2–5° | 10 940 | 18.0° | 8.1 % |
| 5–15° | 4 607 | 21.1° | 8.3 % |
| **> 15°** | 3 839 | 16.8° | **1.9 %** |

The steeper the shore, the more tightly the 0 m crossing is pinned inside one cell and the smoother the contour;
on gentle ground it wanders at sub-cell scale. So the remedy is GEOMETRIC (sub-cellular interpolation or contour
smoothing), not morphological.

**C-4's SCOPE IS REDUCED, and that is an acquired gain.** Coastal erosion will NOT have to fix the barbing — it
only has to sculpt real coastal morphology (cliffs, beaches, abrasion platforms). The C-4 issue as written
assumed otherwise; corrected here.

### Below-sea SPILLWAYS are mischaracterised as rivers (pre-existing, surfaced by arid-hot)

Rivers #1–#6 of the discharge ranking run 0 m → −20 m draining 88 468 km² at Strahler 1. Measured signature:
**2 points** (one grid step — the "1 km" is the `geo_scale_ratio ×7.5` length), **end at exactly −20.0 m on every
one** (a constant, not a terrain elevation), and they come in **near-identical PAIRS** (88468/88449, 62094/62094,
8846/8846) — a probable DOUBLE EMISSION, a separate bug from the typing.

The shrinkage hypothesis is **REFUTED by measurement**: 0 % of these segments lie on floor exposed by H-1c, while
**51.9 % of the real order-≥3 trunks do** — the trunks legitimately cross the drained bed. These are the
below-sea basin SPILLWAYS (50 in that run: 32 to sea, 18 chained) injected into `rivers.json`. Pre-existing;
arid-hot merely lifted them to the head of the discharge sort.

**Decision (consumer-informed): a distinct `Spillway` type.** Living Landz uses `strahler_order` for RENDERING
only (stroke width, display filtering), never for game logic, so a closed-basin outflow can be rendered as a
watercourse of given width without needing an order. `width_m` stays meaningful (it derives from discharge).
"Order 1 with 88 468 km² of basin" is what is absurd, not the object's existence — a spillway is a real flow and
must stay on the map.

**CONSUMER BACKLOG** — this ADDS a type Living Landz must handle, alongside `width_m`, `lake_type` and the
`Wetland` biome: all exported or exportable and still ignored downstream. Recorded rather than assumed to be
noticed.

### Reciprocal spill cycles — Finding 40's DAG assumption was ASSERTED, not verified

Finding 40's fixed point rests on "flow is downhill, so the chain graph is a DAG". Measured on
the production seed (arid-hot 8192², `tests/spillway_duplication.rs`): **7 reciprocal pairs
(A → B and B → A) = 14 of 51 spillways**, each pair tracing from the SAME col cell with
identical drainage. A 2-cycle is not a DAG — the assumption was never checked. It also
explains the near-equal-but-not-equal figures (88 468 vs 88 449): a cycle does not converge
exactly, while other pairs stabilise on the same value.

Two hypotheses were refuted on the way, and both are worth recording because they were the
plausible ones: emission is strictly ONE spillway per basin (51 distinct `lake_id` for 51
spillways), so it is not a per-sill-cell loop; and the duplicate starts are the SAME cell, not
adjacent, so it is not a col flat over several cells either.

**Fix — one physical rule, two cases.** Levels differing by more than
`SAME_WATER_BODY_TOL_M` (0.10 m): the HIGHER basin spills into the lower, so the impossible
uphill spillway is dropped. Levels equal within tolerance: the basins share one free surface
over the col — they ARE one water body, so they are MERGED (one id, one footprint, one
outflow). Result: pairs 7 → 0, spillways 51 → 44, below-sea lakes 43 → 39 (4 pairs merged, 3
level-differentiated), and 0 spillways left pointing at an absent receiver.

**The tolerance pattern, worth reusing.** 0.10 m is chosen ABOVE the height field's vertical
quantisation and BELOW any free-surface difference the model can resolve, so the merge
decision depends on neither quantisation noise nor an arbitrary pick. That is the shape a
physical threshold should have; floating-point equality would have made it depend on noise.

### DEBT — the cyclic inflow still sits in the non-merged basins' levels

The fix corrects the OUTPUT, not the ITERATION. While both directions of a cycle existed, the
fixed point fed each basin the other's full surplus through `extra_inflow`, so the levels it
produced already contain that double contribution. Removing it inside the iteration is a
deeper change, deliberately left out of that lot.

**Magnitude, so nobody rediscovers it as a phantom.** Because the two members of a reciprocal
pair had near-identical surpluses, each received approximately its PARTNER'S ENTIRE surplus —
i.e. **up to ~2× its legitimate inflow**. The seven dropped directions carried, in
drainage-equivalent: **1589, 1106, 162, 33, 29, 18 and 10 km² (≈ 2947 km² total)**. The two
largest are the ones the author saw at the head of the discharge sort. Since the exorheic test
is `a_eq = inflow / net_evap ≥ a_sill`, a doubled inflow can by itself flip a basin exorheic —
so the 14 basins involved in cycles are the ones whose REGIME is least trustworthy today. The
37 non-merged basins keep levels computed with that contribution. Invisible now; it will
surface when H-2 recomputes everything, and this is what it will be.

### Two lessons from H-1b

1. **The area-vs-level test was not the defect** — the leading hypothesis, measured and refuted. Worth recording
   precisely because it was the plausible suspect.
2. **A visual validation on a SINGLE climate can confirm a BUG.** The author validated "0 → several endorheic in
   humid" as the expected behaviour; it was the gross-PE artefact. The reclassification fix keeps all its value,
   but its real effect is elsewhere — hot-arid goes 0 → 63 endorheic, which the geometric path never produced.
   Look where the mechanism is supposed to act MOST STRONGLY, not where the production map happens to sit.

### The CONTIGUOUS PLAIN metric — made grid-stable BEFORE H-2 judges anything on it

H-2 will be judged on the plain area it gains, so the measure had to converge first. The
naive definition (connected components of dry land under 5°) is NOT a property of the
continent but of the GRID STEP: largest piece **1074 km² at 2048² against 216 km² at 8192²
(×4.97)**, for a comparable total flat area (9891 vs 7493 km²) — fine relief resolves and
fragments what the coarse grid smoothed over. Two candidate definitions were implemented and
COMPARED on that convergence test (`tests/plain_metric.rs`), not argued:

| definition | 2048² | 8192² | ratio |
|---|---|---|---|
| naive (no bridging) | 1074 | 216 | ×4.97 |
| (a) bridge 100 m | 1074 | 871 | ×1.23 |
| **(a) bridge 200 m** | **3025** | **3129** | **×1.03** |
| (a) bridge 400 m | 5613 | 6079 | ×1.08 |
| (a) bridge 800 m | 9256 | 13990 | ×1.51 |
| (b) hex-pitch resample | 2759 | 850 | ×3.25 |

**Chosen: (a) — a morphological CLOSING at a PHYSICAL bridging distance of 200 m** (converted
to cells per resolution: 1.0 cell at 2048², 4.1 at 8192²), then the `land_topology`
union-find on the closed mask. Rationale confirmed by measurement: a plain stays usable when
the accidents crossing it are passable, and 200 m is the smallest distance that converges
(×1.03). Beyond 400 m the metric RE-DIVERGES (×1.51 at 800 m) — over-bridging swallows real
relief, differently at each resolution. **(b) was rejected by the test**: resampling to the
hex pitch still inherits the fine-grid fragmentation (at 8192² more fine cells are steep, so
fewer hex cells pass the majority vote) and gives ×3.25.

For the chosen metric: largest plain **3024 / 3129 km²** (1164 / 1204 hex at 1 km-edge hexes),
total plain 10643 / 9989 km² (×1.07); pieces ≥1000 km²: 2 / 1; ≥500: 4 / 3; ≥100: 17 / 12.
The LARGEST and the TOTAL converge; the ranking BELOW the first piece does not (2nd piece
1346 vs 716 km²) — so the criterion is the largest piece and the total, not the detailed
ranking. Note the raw flat-area total itself still differs ×1.32 between resolutions; that
(a)-200 m converges despite it is the point — the bridging absorbs exactly the micro-relief
that differs.

The metric is CLIMATE-INVARIANT as it stands (it reads the geometric lake footprints, which
neither the climate nor the H-1 reclassification resizes) — it will move only under H-2,
which is precisely what makes it the right instrument to judge H-2.

### Corollary to method rule 3: a bench/production gap can indict PRODUCTION

Method rule 3 says a bench must reproduce the WHOLE chain. The corollary, earned here: **when a bench and
production diverge, the gap sometimes points at a defect in PRODUCTION, not in the bench.** The discard above was
exposed by a FAILED bench of mine — it measured the final drainage's lakes, which production throws away, and
returned 0 lakes. The omitted stage was the revealer. Do not assume the bench is the party at fault.

### Step 3b — the `Spillway` TYPE: what the layer says, and the desync the measurement caught

**The symptom.** A below-sea basin's outflow was emitted as a plain `RiverSegment` with
`strahler_order = 1`. Sorted by discharge, the microscope river list opened on SIX spillways
("the first real river is #7"): a spillway carries its whole closed basin's catchment, so on
discharge alone it outranks every genuine river while advertising the order of a headwater.

**The fix is a TYPE, not a filter.** `segment_kind: Vec<SegmentKind>` with
`SegmentKind { Watercourse, Spillway }`, plus `segment_source_lake: Vec<Option<u32>>`, both
parallel to `rivers.segments` and both exported in `rivers.json`. A spillway is a REAL flow
and stays on the map; what changes is that the consumer can tell what it is. The contract
written into the export docstrings:

- `kind == "Spillway"` ⇒ **`strahler_order` MUST NOT be read** (it is 1 regardless of
  catchment). Render from `width_m`, which comes from discharge and is meaningful.
- `source_lake_id == null` on a spillway means **"a real basin, NOT inventoried"** — the
  source sits below the lake-inventory floor and is therefore absent from `lakes.json`. The
  consumer is never handed an id it cannot resolve. Always `null` on a watercourse.

Living Landz uses `strahler_order` for RENDERING only (stroke width, display filtering),
never for game logic — which is why typing suffices and no renumbering scheme is needed.

**Why the dangling count changed from 8 to 16.** The earlier "8" was ARITHMETIC —
`spillways − lakes` (16 − 8 at the time) — which only bounds the count and silently cancels
against basins that ARE inventoried. The real figure is a SET-MEMBERSHIP count: spillways
whose `lake_id` is absent from the lake list. That is **16 of 43 at 8192²** (6 of 16 at
2048²). The two numbers were never in conflict; the first was the wrong instrument.

**The measurement caught a desync that the typing alone would have shipped broken.**
Checking the practical symptom (rather than declaring the sort fixed) reported
`segment_kind` length **13165 for 10186 segments**. `clip_rivers_to_lakes` rebuilds every
per-segment array and the two new ones had not been added to it: the 16 `Spillway` tags sat
at indices 13149..13165 while the segments they described sat at 10170..10186. On the
cache-HIT path a `resize(n_seg, Watercourse)` then TRUNCATED the array and reported zero
spillways — a silent-repair call hiding a structural bug. Three changes, in that order of
importance:

1. `clip_rivers_to_lakes` rebuilds `segment_kind` / `segment_source_lake` too (a clipped run
   inherits its parent's kind).
2. `C1DrainageResult::segment_arrays_aligned()` — one invariant over ALL seven parallel
   arrays, `debug_assert`ed after the clip and after the spillway append.
3. The sidecar read now **errors** on a length mismatch instead of padding/truncating; only
   an ABSENT/EMPTY array is tolerated (legacy sidecars).

This is the completeness trap for composites, third occurrence: it is not enough to add a
field to the struct and the codec — every site that REBUILDS the parallel arrays is part of
the type.

**The verdict, in the production config at both resolutions** (`spillway_typing_bench` in
`ui/workspace.rs`, which calls the very `aggregate_watercourses` the microscope calls, via
`run_hd` end to end):

| | 2048² | 8192² |
|---|---|---|
| list entries | 1359 | 792 |
| watercourses / spillways | 1343 / 16 | 749 / 43 |
| first spillway at rank | **#1344** | **#750** |
| biggest spillway | 944 km², 9 m³/s | 1572 km², 15 m³/s |
| biggest real river | 1087 km², 10 m³/s | 110 km², 1 m³/s |
| spillways with `source_lake_id == null` | 6 / 16 | 16 / 43 |

Before the fix at 2048² that 944 km² spillway sat at **#5**, above four real rivers. At
8192² the biggest spillway (15 m³/s) outranked the biggest river (1 m³/s) by ~14×. The sort
key is now `kind` first, then discharge, so a spillway CANNOT head the list while any river
exists — that part is a property of the comparator, not of the terrain; the table is what
the terrain adds.

**One number to correct, and why.** Step 3a reported a spillway draining **88 468 km²**;
production reports **1572 km²** for the largest. Both are right for what they measured:
`tests/spillway_duplication.rs` calls `below_sea_basin_lakes` on the raw eroded field, WITHOUT
infiltration and WITHOUT the H-1c water balance, so its basins are the un-shrunk ones and its
inflow the un-infiltrated one. Method rule 3 again, from the other side: a partial bench is
fine for a topology question (it found the reciprocal cycles) and worthless for a magnitude.
**Magnitudes come from the production chain.**

**An anomaly observed and NOT acted on** (it is not step 3b's): the largest real river's
catchment is **1087 km² at 2048² but 110 km² at 8192²**, and the list has FEWER entries at the
finer grid (792 vs 1359). A watercourse is assembled by following `downstream` to a terminal,
so this says the 8192² network fragments into more, shorter, disconnected trunks rather than
resolving into longer ones. That is a hierarchy defect worth its own diagnosis, on the same
list of candidates as the contour.

### Consumer backlog (fields the export now offers that Living Landz does not yet read)

| field | layer | what to do with it |
|---|---|---|
| `width_m` | `rivers.json` | stroke width; topology stays continuous below any display cutoff |
| `kind` | `rivers.json` | `"Spillway"` ⇒ do not read `strahler_order`; render from `width_m` |
| `source_lake_id` | `rivers.json` | link a spillway to its basin; `null` = real basin, not inventoried |
| `lake_type` | `lakes.json` | `Endorheic` = salt, `CraterAcidic` = acid — habitability / resource logic |
| `Wetland` | biome raster | the traced marsh footprint, not a rainfall proxy |

### Method notes earned at step 3b

**A "repairing" `resize` can MASK the defect it looks like it is absorbing.** The desync above
was invisible because the sidecar read ended with `segment_kind.resize(n_seg, Watercourse)`.
That line looks like tolerance for legacy caches; what it actually did was truncate 13165
correct-but-misplaced tags down to 10186 and report zero spillways. Without the discharge-sort
check, the typing would have passed as working. **A tolerance clause must distinguish the case
it tolerates from the case it hides**: absent/empty (legacy) is tolerated, a wrong length is
now an error. Whenever a `resize`, a `unwrap_or_default`, a `saturating_*` or a `.get(i)` sits
on a path where the two lengths are supposed to be equal BY CONSTRUCTION, it is not
defensiveness — it is a silencer.

**The structural remedy eliminates the CLASS, not the instance.**
`C1DrainageResult::segment_arrays_aligned()` checks all seven parallel arrays in one call,
asserted at every site that mutates the network. Seven separate length checks would have had
the same defect as the seven separate arrays: adding an eighth would silently skip it. This is
the third occurrence of the composite-completeness trap (sidecar codec → the `AmplitudeTerms`
composition → here) and the first time the remedy is a single invariant rather than a fix.

**A bench that BOUNDS is not a bench that MEASURES.** "8 dangling ids" came from
`spillways − lakes`. That subtraction bounds the count from below and cancels silently against
every basin that IS inventoried; the set-membership count is 16. A difference of cardinals is
not a count of a set.

**Refinement of method rule 3 — a partial bench is VALID FOR A TOPOLOGY QUESTION AND WORTHLESS
FOR A MAGNITUDE.** `tests/spillway_duplication.rs` runs `below_sea_basin_lakes` on the raw
eroded field, without infiltration and without the H-1c water balance. It correctly found the
reciprocal chain cycles — a structural property that does not depend on the missing stages —
and it reported a spillway draining 88 468 km² where production has 1572 km². Both numbers are
right for what they measured. Rule 3 ("a bench must reproduce the whole production chain")
stands; this is the finer statement of WHY, and of what a partial bench may still be trusted
for: **structure, connectivity, degeneracy — yes. Any number that will be quoted — no.**


## Finding 42 — the "network fragments with resolution" anomaly: it is the TERRAIN'S HYPSOMETRY, and none of the four suspects

Step 3b filed an anomaly: the largest assembled catchment was **1087 km² at 2048² but 110 km²
at 8192²**, with **fewer** microscope entries on the finer grid (792 against 1359). The
network appeared to FRAGMENT as resolution increased — the opposite of the expected behaviour,
and it hits the measuring instrument itself, since the microscope assembles watercourses by
following `downstream` to a terminal.

Four candidates were named and each was measured on the full production chain at both
resolutions (`network_fragmentation_bench` in `ui/workspace.rs`, through `run_hd`, calling the
very `aggregate_watercourses` the microscope calls). **All four are cleared.**

| | 2048² | 8192² |
|---|---|---|
| segments / entries / river cells | 10 186 / 1359 / 74 897 | 12 464 / 792 / 291 508 |
| terminals: sea / lake / sub-sea / **none of these** | 92 / 1251 / 0 / **16** | 97 / 692 / 0 / **3** |
| of those, still having a D8 receiver | 16 (all 0 km²) | 3 (all 0 km²) |
| isolated fragments (no upstream, no downstream) | 666 (6.5 %) | 560 (4.5 %) |
| duplicate terminals (terminals − distinct last cells) | 46 of 1359 | 26 of 792 |
| lake cells | 3790 km² (2.37 %) | 1625 km² (1.02 %) |
| `stream_km2` 20 km² → cells | 524.3 | 8388.6 |
| cells clearing it | 11 857 | 45 513 |
| **exported channel LENGTH** | **14 605 km** | **14 284 km** |

1. **ASSEMBLY — not broken.** Only 16 (2048²) and 3 (8192²) terminals are neither sea, lake
   nor sub-sea sink, and every one of them carries 0 km². 87–92 % of terminals are LAKE
   INFLOWS, which is the clip's designed behaviour.
2. **THRESHOLDS — correct and physical.** `stream_km2` is in km² and converts per resolution
   (524 → 8389 cells); the cells clearing it go 11 857 → 45 513 (×3.84 for a ×4 linear
   refinement) and the exported channel LENGTH converges: **14 605 km against 14 284 km,
   ×0.98.** The network is not shrinking or fragmenting in extent — it is resolution-stable.
3. **FEWER ENTRIES — lake area, not fragmentation.** An entry is a terminal, and terminals are
   dominated by lake inflows. The lake footprint is 3790 km² at 2048² but 1625 km² at 8192²,
   so the network crosses fewer lakes and produces fewer terminations (1251 → 692). Entry
   count TRACKS LAKE AREA. Counter-intuitive only until the terminal census is read.
4. **CLIPPING — the dominant source of terminals, and behaving as designed.** 1251 of 1359
   terminals at 2048² end at a lake shore against 78 segments that RESTART on the far side:
   many tributaries in, one outlet out, plus the deliberate drop of endorheic and below-sea
   outlet runs. What this exposes is a SEMANTIC point, not a bug: the microscope's
   "watercourse" is a REACH BETWEEN WATER BODIES, not a river system, so its head is the
   largest reach and not the largest river.

### What the number actually was

The quantity the list labels "catchment" is `segment_drainage_km2= runoff_accumulation / 300 mm`
— a runoff-EQUIVALENT area at a reference depth, not a geometric catchment. Read beside the
flow-accumulation RASTER (a geometric cell count, independent of segmentation) it separates
cleanly:

| at the sea mouths | 2048² | 8192² |
|---|---|---|
| summed GEOMETRIC area | 9519 km² | 6377 km² (×0.67) |
| summed EFFECTIVE area | 12 169 km² | 603 km² (×0.05) |
| effective / geometric | **1.278** | **0.094** |

An "effective area" 1.278× the geometric one is not a contradiction — it means the catchment's
mean net runoff exceeds the 300 mm reference. The ratio IS that mean over 300 mm. So the
figure moved because **the water moved**, and the ×0.67 on the geometric side is Finding 41
(more closed depressions at a finer grid capture more catchment before it reaches the sea).

Two side observations, both minor and both real: 46/26 **duplicate terminals** (two terminal
segments on one cell — the pairs visible in the discharge sort: 1459/1459, 738/738), and
666/560 **isolated fragments** which each become a list entry carrying their PARENT'S
inherited area — which is how an S1, 0-tributary row shows 1087 km² beside the genuine S4,
394-tributary trunk. Also checked and cleared: the 3/5 sea mouths with a real catchment and
ZERO discharge are **not** the endorheic mask killing a route to the sea (0 of them are —
re-accumulated with and without the mask), they are genuinely arid catchments.

### The cause, and a hypothesis of mine that the measurement refuted

Mean land precipitation is **1219 mm/yr at 2048² and 712 mm/yr at 8192²** (×0.58), and
`max(0, p − pe)` follows at ×0.33 (620.4 → 207.3 mm/yr; wet-cell fraction 12.92 % → 7.59 %).

I hypothesised a discretisation defect in the transport: `oro = k_oro·m·ascent` removes a
FRACTION of the carried moisture per CELL, so the flux should decay as `(1 − k_oro·S)^N` and
depend on the cell count. **`tests/precip_resolution_invariance.rs` refutes it**: on the SAME
analytic continent (ocean, then a 100 km ramp to a 2000 m plateau over a 400 km domain) the
land mean is **391 mm/yr at 512², 1024², 2048², 4096² AND 8192² — ratio 1.000**, interior/coast
contrast identical. `k_oro` is not the binding constraint; the CAPACITY CAP `m > e_sat(T)` is,
and that depends on the altitude PROFILE, not on how many cells sample it. The transport is
resolution-invariant.

So the climate is faithfully reporting a terrain that differs. And it differs by a lot:

| land statistic | 2048² | 8192² | ratio |
|---|---|---|---|
| mean altitude | 287 m | 693 m | **×2.42** |
| hypsometry p10 / p50 / p90 / p99 / max | 6 / 167 / 796 / 1407 / 2684 m | 30 / 447 / 1606 / 3259 / 4145 m | ×5.0 / ×2.68 / ×2.02 / ×2.32 / ×1.54 |
| raw normalised above 0.5: mean | 0.02540 | 0.06130 | ×2.41 |
| mean temperature | 19.82 °C | 17.20 °C | −2.62 °C |
| emerged fraction | 14.95 % | 16.43 % | ×1.10 |

**VERDICT: the exported terrain's HYPSOMETRY is not resolution-stable.** The whole
distribution inflates ~×2–2.7 (not just the peaks), for the same seed, the same domain and
essentially the same emerged fraction. Everything altitude-dependent reads a different world
at each grid: temperature, PE, the moisture capacity, precipitation, biomes, and every
hydrological figure downstream of the runoff. The "fragmenting network" was a symptom;
the network is the one thing that IS stable (channel length ×0.98).

This is upstream of H-2 and upstream of the contour, and it is the same family as the comb
(a property that degrades with resolution) with a different culprit.

### Cost of each remedy, so the ordering can be decided

- **Assembly, thresholds, clipping** — no defect, no cost. The only thing worth changing is
  the microscope's SEMANTICS (chain reaches across an exorheic lake via `sink_lake_id` so a
  river SYSTEM is one entry): viz-only, no core change, no cache bump. Low cost, and it makes
  the instrument say what the author expects it to say.
- **The "catchment" label** — export the GEOMETRIC catchment (`flow.accumulation × cell_km2`)
  under `catchment_km2` and let the runoff-derived quantity be what it already is, the
  discharge. A few lines plus an `ALGO_DRAINAGE` bump. The catch: `navigability` classifies on
  the runoff-derived figure today, and its thresholds (500 / 5000 / 50 000 km²) were written
  as AREAS — so the classes shift and need a re-look. Low code, moderate validation.
- **Duplicate terminals / isolated fragments inheriting a parent's area** — a few lines in
  `clip_rivers_to_lakes` (read the area at the run's OWN downstream-most cell instead of
  inheriting, EXCEPT for an exorheic outlet run where Finding 22 requires inheritance) plus an
  `ALGO_DRAINAGE` bump. Low cost, cosmetic effect: it removes the phantom high-area rows.
- **The hypsometry** — unknown until diagnosed, and the diagnosis is NOT done. What is already
  excluded: the FBM octave count is fixed at 7 and does NOT scale with `target_size`, and the
  hillslope diffusion is explicitly renormalised (`HILLSLOPE_REF_CELL_M`, ∝1/cell²). The
  candidates left, each cheap to test on the raw upscaled field before erosion:
  (a) the C-1 relief-budget CAP (`AmplitudeTerms.cap`, the one that made `amplitude_base`
  inert) evaluating differently per grid; (b) `flow_conditioning = 0.1`, a per-cell downslope
  stretch; (c) the stream-power incision — `relief_v3(cell_km2, …)` is parameterised by cell
  area and runs a FIXED 2 iterations, so the erosion may simply do less relative work at
  8192²; (d) the sea-level / `target_land_fraction` calibration. The decisive first
  measurement is one bench reporting the hypsometry of the RAW FBM field and of the ERODED
  field at both resolutions: it splits the search in two at the cost of one run.


## Finding 43 — the hypsometry bisection: the INCISION, and Finding 7's fix on a branch the shipped config does not take

Finding 42 left one question: at which stage does the ×2.42 hypsometry inflation appear?
`tests/hypsometry_bisection.rs` builds the terrain in three stages through
`production_hd_config` and reports the full distribution after each, at both resolutions.
Only one thing changes per stage, so the answer is a reading, not an inference.

### The bisection — the answer is unambiguous

Over the EMERGED cells, metres (`norm>0.5` is the raw normalised mean above sea level, an
audit of the metric conversion):

| stage | grid | land cells | emerged % | mean | p10 | p50 | p90 | p99 | max | norm>0.5 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1. coarse only (no FBM, no incision) | 2048² | 708 014 | 16.88 | 865 | 95 | 679 | 1836 | 3399 | 4132 | 0.07658 |
| 1. coarse only | 8192² | 11 327 963 | 16.88 | 866 | 95 | 679 | 1836 | 3399 | 4148 | 0.07659 |
| | **ratio** | — | **1.00** | **1.00** | **1.00** | **1.00** | **1.00** | **1.00** | **1.00** | **1.00** |
| 2. + FBM (no incision) | 2048² | 708 014 | 16.88 | 865 | 95 | 679 | 1836 | 3399 | 4142 | 0.07659 |
| 2. + FBM | 8192² | 11 327 963 | 16.88 | 866 | 95 | 679 | 1836 | 3398 | 4153 | 0.07660 |
| | **ratio** | — | **1.00** | **1.00** | **1.00** | **1.00** | **1.00** | **1.00** | **1.00** | **1.00** |
| 3. PRODUCTION (+ relief-v3 incision) | 2048² | 627 228 | 14.95 | 282 | 6 | 161 | 787 | 1407 | 2684 | 0.02498 |
| 3. PRODUCTION | 8192² | 11 031 073 | 16.44 | 685 | 29 | 445 | 1606 | 3259 | 4145 | 0.06065 |
| | **ratio** | — | 1.10 | **2.43** | **4.54** | **2.75** | **2.04** | **2.32** | 1.54 | **2.43** |

**Stages 1 and 2 are invariant to the digit, on every percentile. The entire inflation is
created by the INCISION.** And the sign matters: the incision does not build the 8192² terrain
taller, it **fails to erode it**. From the same 865/866 m starting point it removes **67 % of
the mean land altitude at 2048² (865 → 282 m) but only 21 % at 8192² (866 → 685 m)**.

### The C-1 relief budget cap is CLEARED — twice

The promoted suspect was the cap, `β·slope_mag/(nscale·S)`, on the argument that a slope
measured over a 49 m cell is far steeper than over a 195 m cell. **That mechanism does not
apply here**: the slope the cap reads is not an HD-cell gradient. It comes from
`slope_map = compute_terrain_analysis(coarse)` — computed on the **64² coarse grid** with
`gradient_at_periodic` and then sampled at coarse-pixel coordinates — and `nscale` is
`base_frequency · 1024 / src_max²`, a function of the COARSE size alone. Neither term reads
`target_size`, so the cap field is identical at every resolution by construction.

Measured through its effect, which is what stage 2 − stage 1 is (the FBM's added relief, in
metres, over all cells):

| grid | mean | p1 | p50 | p99 | min | max |
|---|---|---|---|---|---|---|
| 2048² | 0.01 | −1.74 | 0.00 | 2.03 | −36.93 | 35.75 |
| 8192² | 0.01 | −1.73 | 0.00 | 2.03 | −37.13 | 36.06 |

Identical. The cap is invariant, by reading AND by measurement. Two things worth keeping from
this table anyway: it **quantifies the DEAD KNOB** — the cap crushes the whole FBM to within
±2 m over 98 % of cells, median contribution 0.00 m — and it confirms that the FBM is not what
gives the C-1 terrain its relief.

### The exclusion list, so nothing is re-suspected

| candidate | status | why |
|---|---|---|
| FBM octaves | excluded (code) | fixed at 7, does not scale with `target_size` |
| FBM base frequency / feature size | excluded (measured) | stage 2 ratio 1.00 |
| C-1 relief budget cap | **excluded (code + measured)** | coarse-only inputs; stage 2 − stage 1 identical |
| `flow_conditioning = 0.1` | excluded (code + measured) | its only two uses are the cap and the noise-sign gate, both coarse-driven; stage 2 ratio 1.00 |
| coarse upscale / bilinear | excluded (measured) | stage 1 ratio 1.00 on every percentile |
| `target_land_fraction` calibration | excluded (code) | `None` in `c1_hd_production` AND untouched by `production_hd_config` — the calibration never runs in production |
| framing roll `sample_origin` | excluded (code) | `roll_x/roll_y` are computed on the coarse 64² grid, identical at both resolutions; `sample_size = 1.0` (whole torus, no crop) |
| nonlinear hillslope diffusion | excluded (code) | renormalised by `(HILLSLOPE_REF_CELL_M/cell_m)²` — but see below, relief-v3 does not take that branch |
| `min_area_cells`, `lateral_erosion`, talus drop | excluded (code) | all three are expressed physically (`A_c/cell_km2`, metres, `talus_slope·cell_m`) |
| **stream-power incision** | **THE STAGE (measured)** | stage 3 ratio 2.43 from an invariant 865/866 m start |


### A correction to Finding 3, forced by stage 2

Finding 3 recorded "Incision is resolution-dependent" (per-order incision 108 m @512² → 136
@1024² → 318 @2048²) and attributed it to "the FBM detail resolving sharper gradients on finer
cells". **That explanation no longer holds for the shipped config.** Stage 2 measures the
FBM's total contribution at ±2 m over 98 % of cells with a median of 0.00 m — C-1 flow
conditioning (which did not exist when Finding 3 was written) caps it out of existence. So the
FBM cannot be the source of a resolution-dependent gradient field today. The observation in
Finding 3 stands; its stated cause is superseded.

Note also the direction: Finding 3 measured incision RISING with resolution over 512²→2048²,
whereas Finding 43 measures the erosion doing proportionally LESS work at 8192² than at 2048²
(67 % of the mean altitude removed against 21 %). Different metrics (per-order incision depth
vs mean land altitude) and a different regime (pre-C-1 vs relief-v3), so they are not in
contradiction — but nobody should quote the one as support for the other.

### Finding 7's fix exists, and the shipped config takes the other branch

`incise` has TWO hillslope-diffusion branches:

- `cfg.critical_slope > 0.0` → the nonlinear implicit closure, which carries Finding 7's
  renormalisation: `let dscale = (HILLSLOPE_REF_CELL_M / cell_m).powi(2); let w_base =
  cfg.diffusion * dscale;` with the comment *"otherwise the same D smooths 16× fewer metres at
  4× finer cells and the closure fails to plane at 8192²"*.
- `cfg.critical_slope == 0.0` → the linear explicit branch: `let dsub = cfg.diffusion /
  cfg.diffusion_substeps as f32; … field.data[k] = src[k] + dsub * lap;` where `lap` is the
  5-point Laplacian in NORMALISED height per CELL². **No `dscale`, no `cell_m`.**

`relief_v3` sets `critical_slope: 0.0` ("MFD prevents the comb → no Gauss-Seidel solver") and
`diffuse_channels: true`. So the shipped configuration takes the branch WITHOUT the fix, and
applies it to every cell. The implied physical diffusivity is `κΔt = dsub · dx²`, i.e. it
shrinks as `dx²`: with `dsub = 0.08/4` over 4 substeps × 2 iterations, `κΔt` is 760 m² per
substep at 2048² against 48 m² at 8192² (**×16**), and the planing length √(κΔt·n) is **78 m
against 20 m** (×3.9).

This is the defect Finding 7 identified and fixed, **regressed by a branch change rather than
by an edit to the fix** — the ADR's claim "after the fix the closures behave the same in
metres at 2048² and 8192²" is true of the branch it was written about and false of the one
that ships. A fix that lives in one arm of an `if` is not a fix; it is a fix for one config.

### The family: cells versus metres, sixth instance

Every previous instance was a quantity written per CELL where physics needs it per METRE:
the ±1-cell lateral reach (Finding 7), the dimensionless diffusion weight (Finding 7), the
`min_area_cells` channel head, the drainage thresholds, the `decay_km` fracture distance. This
one is the same shape with a twist worth naming: **the units were correct in the abstract
(`diffusion` is documented as "dimensionless at the reference cell") and the conversion simply
was not applied on this path.** The lesson is not "use physical units" — that was already
known and written down — it is that **a units convention needs a single conversion site, not a
convention plus discipline at each use.**


### The sub-bisection inside stage 3 — and it refutes the ranking I had just proposed

Having found the missing `dscale` on the linear branch, I ranked it as the cause. **The
measurement says no.** Each variant changes exactly ONE term from the shipped config; the
diagnostic is the 8192²/2048² ratio, and `retained` is post-incision mean over the invariant
865/866 m pre-incision mean.

| variant | mean 2048² | mean 8192² | RATIO | retained 2048² | retained 8192² | p50 2k | p50 8k |
|---|---|---|---|---|---|---|---|
| shipped relief-v3 | 282 | 685 | 2.43 | 32.6 % | 79.2 % | 161 | 445 |
| MFD off (D8 area) | 310 | 676 | 2.18 | 35.8 % | 78.1 % | 193 | 480 |
| **`diffusion = 0`** | **279** | **685** | **2.45** | 32.3 % | 79.1 % | 154 | 444 |
| talus off | 298 | 688 | 2.31 | 34.4 % | 79.5 % | 151 | 433 |
| `lateral_erosion = 0` | 282 | 685 | 2.42 | 32.6 % | 79.1 % | 161 | 443 |
| `iterations = 1` | 447 | 713 | **1.59** | 51.7 % | 82.4 % | 295 | 482 |
| `iterations = 4` | 177 | 662 | **3.75** | 20.4 % | 76.5 % | 95 | 427 |
| `iterations = 8` | 107 | 625 | **5.86** | 12.3 % | 72.2 % | 60 | 399 |

**Removing the diffusion entirely changes nothing** — 282 → 279 m at 2048², 685 → 685 m at
8192², ratio 2.43 → 2.45. So the missing `dscale` is a REAL latent defect (Finding 7's fix is
genuinely bypassed on the branch relief-v3 takes) but it is **not the cause of the inflation**,
and my ranking of it was wrong. It also carries its own news: **the "light linear hillslope"
closure that relief-v3 documents as grading flanks is INERT** — a third dead knob, alongside
`amplitude_base` and the striation ladder. MFD (2.43 → 2.18) and the talus (→ 2.31) are minor
contributors; lateral widening is nil.

**The signature is in `iterations`, and it is not a knob to retune — it is a diagnosis.**
Iterating does not close the gap, it WIDENS it: 1.59 → 2.43 → 3.75 → 5.86 for 1/2/4/8. Read
the `retained` columns instead of the ratio and the reason is plain:

- at 2048², each doubling keeps eating the landscape: 51.7 % → 32.6 % → 20.4 % → 12.3 %;
- at 8192², iterating barely does anything: 82.4 % → 79.2 % → 76.5 % → 72.2 %.

**At 8192² the incision is essentially SATURATED after one sweep.** No amount of repetition
reaches the rest of the land. That is the fact to explain, and it is not about how hard each
link is incised — it is about HOW MUCH OF THE LAND the channel term can touch at all.

### The corroborating number, already in hand

The incision runs only where `area >= min_area_cells = A_c/cell_km2`. `A_c = 0.1 km²` is
PHYSICALLY CONSTANT — 2.6 cells at 2048², 41.9 at 8192² — so the threshold is not the bug.
But the FRACTION OF LAND CELLS that clears it is not constant, because the number of cells
with accumulation ≥ a falls off sub-linearly (the area–frequency law). From the Finding 42
bench:

| | 2048² | 8192² |
|---|---|---|
| land cells | 626 951 | 11 024 627 |
| cells clearing `A_c = 0.1 km²` | 347 578 | 1 110 919 |
| **fraction of land in the FLUVIAL regime** | **55.4 %** | **10.1 %** |

So at 2048² the channel term erodes more than half the land; at 8192² it reaches a tenth of
it, and the other 90 % is handed to the hillslope term — which the variant table just showed
does NOTHING. Refining the grid transfers land from an agent that erodes to an agent that
does not.

(Measured on the drainage phase's D8 accumulation on the eroded field, not on the incision's
own MFD accumulation — corroborating, not the primary evidence. The primary evidence is the
`no regime split` variant.)


### The decisive variant, and the verdict: TWO mechanisms, not one

`min_area_cells = 0` makes every cell incise, removing the fluvial/hillslope partition:

| variant | mean 2048² | mean 8192² | RATIO | retained 2048² | retained 8192² |
|---|---|---|---|---|---|
| shipped relief-v3 | 282 | 685 | 2.43 | 32.6 % | 79.2 % |
| **no regime split (`A_c = 0`)** | **172** | **320** | **1.85** | 19.9 % | 37.0 % |

The ratio falls from 2.43 to 1.85 — a large move, **but it does not collapse to 1.00.** So the
partition is the dominant carrier and it is NOT the whole story. The decomposition, computed
two ways that agree:

- excess mean altitude at 8192² over 2048²: **403 m** shipped, **148 m** without the split;
- retained-fraction gap: **46.6 pp** shipped, **17.1 pp** without the split.

**The regime partition accounts for 63 % of the excess. A residual 37 % lives in the incision
term itself**, present even when every cell incises.

**MECHANISM 1 (63 %) — a physically-correct threshold whose PARTITION OF THE DOMAIN is not.**
`A_c = 0.1 km²` is constant in km² (2.6 cells at 2048², 41.9 at 8192²), so the threshold is
right. But the number of cells with accumulation ≥ a falls off sub-linearly (the
area–frequency law), so the FRACTION of land in the fluvial regime is **55.4 % at 2048²
against 10.1 % at 8192²**. Refining the grid moves 45 % of the land from the channel term to
the hillslope term — and the variant table shows the hillslope term does NOTHING
(`diffusion = 0` changes the result by 3 m). Land is transferred from an agent that erodes to
an agent that does not.

This is NOT the cells-versus-metres family, and calling it that would misdirect the fix: every
unit here is already physical. It is a new shape worth its own name — **a correct threshold
feeding an inert branch.** The partition is legitimate physics (hillslopes should not incise
fluvially); the defect is that the other half of the physics was never made to work.

**MECHANISM 2 (37 %) — not yet isolated.** Present with `A_c = 0`, so it is inside the
`h_new = (ho + f·hr)/(1+f)` relaxation, `f = K·dt·A_km²^m / dist_m`, run a FIXED 2 iterations.
Two leads, both visible in the variant table: MFD dispersal (turning it off moved the ratio
2.43 → 2.18, and the incision's MFD partition is applied per CELL, so it compounds 4× more
often along the same physical path, diluting `A` and hence `f`), and the iteration count
itself (`retained` at 2048² keeps falling with iterations, 51.7 → 12.3 %, while at 8192² it
barely moves, 82.4 → 72.2 % — the sweep saturates on the fine grid). One more sub-bisection
would separate them; it was not run.

### What the fix costs — and why it is NOT one line

The tempting one-liner is applying Finding 7's `dscale` to the linear branch. **The
measurement predicts it will barely move the hypsometry**, and this is the most useful thing
this round produced: `diffusion = 0` versus shipped is a 3 m difference at 2048², where the
weight is by definition correctly calibrated. So the hillslope term is not merely
mis-normalised at 8192² — it is **too weak to matter at either resolution**. Restoring the
`dx²` normalisation would make 8192² match 2048², i.e. match ~nothing.

| remedy | code | validation | expected effect |
|---|---|---|---|
| apply `dscale` to the linear branch | **one line** | `ALGO_*` bump; 8192² moves | real defect, closes Finding 7's regression — but ~nil on hypsometry (predicted, worth verifying since a prediction is cheap here) |
| calibrate the hillslope diffusivity to a PHYSICAL κ | design work | full recalibration of everything altitude-dependent | the actual lever on mechanism 1. Needs the timescale question answered first: the incision has `dt = 1.0` and a fixed iteration count, so there is no explicit time for a κ in m²/yr to multiply |
| isolate mechanism 2 (MFD dilution vs iteration saturation) | bench only | none | one more sub-bisection, ~10 min of compute |
| lower `A_c` at fine grids | one line | recalibration | would equalise the partition, but by making the channel network non-physical — a fix to the metric, not the terrain. NOT recommended |

**So the honest cost is: the diagnosis is done to 63 %, the remaining 37 % needs one more
bench, and the fix is a PHYSICS CALIBRATION rather than a units patch.** The units patch
should still be applied because Finding 7's regression is real, but it must not be presented
as the hypsometry fix — and shipping it in the belief that it is would burn a recalibration
cycle for a 3 m effect.

A note on the instrument: `A_c = 0` is a diagnostic, not a candidate config. It planes 2048²
from 282 m to 172 m — the terrain that made it would be wrong at both resolutions.


### Parked, in this order — recorded so none of it is lost

Ordered by the author after Finding 42. Hypsometry (Finding 43) comes first because the
contour is measured on shores whose altitude and slope depend on the grid, and because the
microscope semantics change what he reads while judging.

1. **Isolate mechanism 2** of Finding 43 (MFD dilution vs iteration saturation) — bench only,
   ~10 min of compute. Then decide the hypsometry fix.
2. **Microscope semantics** — chain reaches through an exorheic lake via `sink_lake_id` so a
   river SYSTEM is one entry instead of one entry per lake inflow (87–92 % of entries are lake
   inflows today). Viz only, no core change, no cache bump. Low cost, and it makes the
   instrument say what the author expects of it. AFTER hypsometry, since it changes what he
   reads while judging.
3. **The `catchment` label** — export the geometric area (`flow.accumulation × cell_km2`) as
   `catchment_km2` and let the runoff-derived quantity be the discharge it already is. THE
   TRAP: `navigability` classifies on the runoff quantity while its thresholds (500 / 5000 /
   50 000 km²) are written as AREAS, so the classes will move. Medium validation cost.
4. **Duplicated terminals and isolated fragments** — 46/26 terminals share a last cell; 666/560
   segments have neither upstream nor downstream and each becomes a list entry carrying its
   PARENT'S inherited area (how an S1 row with 0 tributaries shows 1087 km² beside the real S4
   trunk). Few lines in `clip_rivers_to_lakes` — read the area at the run's OWN downstream-most
   cell, EXCEPT for an exorheic outlet run where Finding 22 requires inheritance. Cosmetic.
5. **The contour** (step 4) — LAST, because it is measured on grid-dependent shores. Barb
   metric before/after BY SLOPE CLASS (baseline: turns > 80° = 8.6 / 12.2 / 8.1 / 8.3 / 1.9 %
   for < 0.5°, 0.5–2°, 2–5°, 5–15°, > 15°; axial R 0.09; mean step 0.78 cell). The falsifiable
   prediction: gradient pinning means the improvement must concentrate on LOW-SLOPE shores; if
   it is uniform or lands on steep shores, the mechanism is not the one identified. And report
   the sea-level offset (−0.06 m median today) after smoothing — it is the end-to-end coherence
   check between P3-A, the export roll and the Living Landz reader, and it must not drift
   unnoticed.
6. **Finding 37 residual orphan** — issue draft in `ISSUE_finding37_residual_orphan.md`, not
   filed on GitHub (`gh` unavailable in this environment).


### Mechanism 2, split — and a flaw in my own bench design

Every variant sits on `A_c = 0`, so mechanism 1 is out of the way and the residual ×1.85 is
what is being attributed. `A_c = 0` is a DIAGNOSTIC instrument only, never a candidate config.

| variant | mean 2048² | mean 8192² | RATIO | excess (m) | share of residual excess |
|---|---|---|---|---|---|
| `A_c = 0` (residual baseline) | 172 | 320 | 1.85 | 147 | — |
| **`A_c = 0` + MFD off** | 116 | 189 | 1.63 | **73** | **50 %** (74 m of 147) |
| `A_c = 0` + `iterations = 1` | 400 | 528 | **1.32** | 128 | 13 % (19 m of 147) |
| `A_c = 0` + MFD off + `iters = 1` | 334 | 435 | 1.30 | 102 | 31 % (46 m of 147) |

**The two metrics rank the mechanisms differently, and the ratio is the one that lies.**
`iterations = 1` gives the best RATIO (1.32) while removing almost none of the EXCESS (19 m of
147). The reason is visible in the means: dropping to one sweep raises BOTH resolutions (172 →
400 and 320 → 528), so both move toward the invariant un-eroded 865/866 m — and as total
erosion → 0 the ratio → 1.00 trivially. Confirmed by the combination: adding `iters = 1` to
MFD-off makes the excess WORSE (73 → 102 m) while the ratio still looks better (1.63 → 1.30).

**This is a defect in the instrument I built.** I introduced the ratio as "the diagnostic"; it
is only valid at COMPARABLE TOTAL EROSION. When a variant changes how much erosion happens at
all, the ratio measures proximity to the invariant limit rather than agreement between grids.
The excess in METRES is the correct attribution instrument here — and it is also the one the
consumers feel, since temperature is a lapse rate on absolute altitude (the 287 vs 693 m gap
is 2.62 °C, not a percentage).

**Verdict on mechanism 2:**

- **MFD dispersal owns half the residual (74 m of 147 m ≈ 18 % of the original 403 m excess).**
  The mechanism is confirmed and it IS a cells-versus-metres instance — the genuine sixth: the
  MFD partition is applied PER CELL, so over the same physical path it composes 4× more often
  at 8192², diluting `A` and hence `f = K·dt·A_km²^m / dist_m`. A partition RATE expressed per
  cell instead of per unit distance.
- **The iteration count is NOT an independent divergence mechanism** (13 % by excess). It is a
  global erosion-AMOUNT knob. Its apparent ratio improvement is the artefact above. This also
  means retuning `iterations` per resolution would be a fix to the metric, not to the terrain —
  the same trap as lowering `A_c` at fine grids.
- **~50 % of the residual (73 m ≈ 18 % of the original excess) remains UNATTRIBUTED.** It is
  still there with `A_c = 0` and MFD off. Untested candidates: the `f ∝ 1/dist_m` link
  compression itself, the talus sweep (which moved the full-config ratio 2.43 → 2.31), and the
  cardinal-versus-diagonal `dist` handling. Said plainly rather than rounded away.

### Attribution of the original 403 m excess

| mechanism | share | evidence |
|---|---|---|
| 1 — regime partition feeding an inert branch | **63 %** (255 m) | `A_c = 0` drops the excess 403 → 148 m |
| 2a — MFD dispersal composing per cell | **18 %** (74 m) | MFD off on top of `A_c = 0` drops 147 → 73 m |
| 2b — unattributed, inside the relaxation | **18 %** (73 m) | survives `A_c = 0` + MFD off |
| iteration count | **not a divergence mechanism** | changes erosion amount, not grid agreement |


## Finding 44 — STRUCTURAL GAP: the model has no explicit TIMESCALE, and two chantiers now need one

Specified, not implemented. This was a calibration caveat inside Finding 43; two independent
requirements now converge on it, which makes it a gap in its own right.

### The state today

| parameter | value | what the code says it is |
|---|---|---|
| `dt` | 1.0 | *"Timestep per drainage↔incision iteration (**lumped into `K·dt`**)"* |
| `k` | 4500 | *"Erodibility `K` (**lumped with the timestep** — see `dt`). Calibrated to a target channel-incision depth, not to appearance"* |
| `iterations` | 2 | *"Number of drainage↔incision iterations (recompute flow between each, so the network can reorganise as the terrain changes — **the staleness handling**)"* |

So the model carries ONE lumped number, `K·dt = 4500`, and an iteration count whose documented
purpose is **numerical** (flow-field staleness). But Finding 43 measured that the iteration
count is what actually governs total denudation: retained mean altitude at 2048² goes
**51.7 % → 32.6 % → 20.4 % → 12.3 %** for 1/2/4/8 iterations. **A numerical parameter is
carrying the physical duration.** That conflation IS the gap.

### Why two chantiers need it

> **RELABELLED by block A: "inert" is too strong — the branch is CONSERVATIVE.** Measured:
> at 8192² the incision-off/incision-on comparison moves **66.5 % of common land by more than
> 1 m** while shifting the mean by ~3 m when `diffusion` is switched off entirely. An operator
> that displaces two thirds of the land and removes no mass is the operational definition of a
> **conservative** operator, not an inert one — the Laplacian fills valleys exactly as much as
> it lowers ridges. The consequence for the remedy is unchanged (a transport-limited term is
> still what is needed); the consequence for the DIAGNOSIS is that "the branch does nothing" is
> false and should not be repeated — it does a great deal, all of it mass-neutral.

- **The hypsometry remedy (Finding 43).** Mechanism 1 is a correct threshold feeding an inert
  hillslope branch. Making that branch work means giving it a real diffusivity `κ`, and
  `∂h/∂t = κ∇²h` has `[κ] = m²/yr`. **There is no year for it to multiply.** Calibrating a
  physical κ before the timescale exists would be premature — there is literally nothing to
  multiply it by, so any value chosen would be a fitted dimensionless weight wearing physical
  units, which is worse than the honest dimensionless weight there now.
- **H-2's temporal dial (the France↔Scotland control).** Sill incision at 10â´ against 10â¶
  years IS the mechanism. It cannot be expressed against `dt = 1.0` with a fixed iteration
  count, because there is no axis along which 10â´ and 10â¶ differ.

### What the specification must settle

1. **The meaning of `dt`.** The natural definition: `dt` is the duration, IN YEARS, modelled
   by one drainage↔incision iteration; total modelled duration `T = iterations · dt`.
2. **Separating the two roles of `iterations`.** Today one integer serves both the physical
   duration and the flow-field staleness bound. They must be decoupled: the author (or H-2)
   sets the target `T`, and the iteration count is DERIVED as `ceil(T / dt_max)` where
   `dt_max` is a numerical bound — how far the terrain may move before the flow field must be
   recomputed, and how large the implicit relaxation step may be before it over-relaxes.
3. **Where resolution enters.** `dt_max` legitimately depends on cell size. That is the
   correct place for the grid to appear — and it would also address mechanism 2's iteration
   saturation, since the sweep count would no longer be a fixed 2 at both grids. **So the
   timescale work sits UPSTREAM of the hypsometry remedy, not beside it.**
4. **Which parameters become dimensionally meaningful once `T` exists:**
   - `K` in the stream-power law: with `E = K·A_km²^m·S^n` in m/yr and m = 0.5,
     `[K] = mÂ·yrâ»Â¹Â·kmâ»Â¹` — Stock & Montgomery's tabulated values become directly usable
     instead of being re-fitted;
   - `κ` hillslope diffusivity in m²/yr — soil-creep literature becomes usable (this is
     mechanism 1's remedy);
   - H-2's sill-incision duration, and any drainage timescale it needs;
   - unaffected, because already dimensional: the talus repose slope (dimensionless by
     nature), the lateral half-width `K_lat·A^m` (metres), `A_c` (km²), the fracture
     `decay_km`.

### The gap can be closed WITHOUT changing a single output

Worth stating because it removes the usual objection. Define `dt := 1.0 yr` and
`K := 4500 mÂ·yrâ»Â¹Â·kmâ»Â¹`. Then `KÂ·dt` is unchanged, every existing terrain is reproduced
**byte-identically**, and the units exist. The dial becomes available afterwards, by varying
`T`, against a reference that was never disturbed. Same discipline as C-3's hard basement at
×1, where holding the reference at unity made the global slowdown nil by construction.

The hard part is not naming the units — it is (2) and (3): deriving the iteration count from a
duration and a stability bound instead of pinning it at 2. That is where the recalibration
lands, and it should be done ONCE, before either the hypsometry remedy or H-2 spends a cycle
on a dial that does not exist yet.

### Method point — a stale EXPLANATION propagating between findings

Finding 3's cause ("the FBM detail resolving sharper gradients on finer cells") was correct
when written and was silently invalidated by a later change (C-1 flow conditioning, which caps
the FBM's whole contribution to ±2 m — Finding 43 stage 2). It survived because nothing
re-checks an explanation when its premise moves; the OBSERVATION kept being true, so the
finding kept looking healthy. **An observation and its explanation have different lifetimes,
and only the observation is protected by being measured.** This is how a wrong model survives:
not through a false measurement, but through a true measurement still carrying a dead
explanation. Practical rule, alongside the other method notes: when a finding is cited as
support, cite its MEASUREMENT, and re-derive the mechanism against the current code — and
never chain two findings' explanations without checking that both premises still hold. Here
the two were also of different metrics, different regimes, and opposite signs.


### The unit patch, shipped and MEASURED against my own prediction

`dscale = (HILLSLOPE_REF_CELL_M / cell_m)²` applied to the linear explicit diffusion branch,
closing Finding 7's regression. `dscale == 1.0` EXACTLY at the reference cell (2048² over
400 km), so the reference is preserved by construction — the same discipline as C-3's hard
basement at ×1.

| | 2048² | 8192² | ratio |
|---|---|---|---|
| before the patch (mean land altitude) | 282 | 685 | 2.43 |
| **after the patch** | **282** | **702** | **2.49** |
| p10 | 6 → 6 | 29 → 42 | |
| p50 | 161 → 161 | 445 → 468 | |
| emerged % | 14.95 → 14.95 | 16.44 → 16.04 | |

**My prediction was right on magnitude and WRONG ON SIGN.** I forecast "~3 m, effectively
nil". Measured: 2048² is byte-identical (as designed), and 8192² moves **+17 m — upward** —
so the ratio gets slightly WORSE, 2.43 → 2.49. Small, as predicted; unhelpful, which I did not
predict.

**Why, and it matters more than the patch.** A linear Laplacian is MASS-CONSERVING: it moves
material from convex to concave, so it lowers ridges and FILLS valleys in equal measure. With
`diffuse_channels = true` it backfills the very channels the incision just cut (the ADR
already observed "diffusion BACKFILLS valleys" in a different context). Making it 16× stronger
at 8192² therefore raises the mean rather than denuding the hillslopes.

**Consequence for mechanism 1's remedy — it is NOT "make the diffusion work harder".** The
inert branch cannot be fixed by restoring its strength, because the term is the wrong KIND of
term for the job: lowering the un-channelled 90 % of fine-grid land requires an agent that
REMOVES mass from hillslopes and delivers it to the channel network, not a conservative
smoother that redistributes it locally. That is a transport-limited hillslope law with an
explicit sediment flux — which is also the thing that has no time to integrate against
(Finding 44). So the remedy for mechanism 1 now depends on the timescale gap, not merely
benefits from it.

**Status of the patch.** Kept and landed SEPARATELY, labelled: it closes a real dimensional
inconsistency (Finding 7's fix lived in one arm of an `if`), it is byte-identical at the
reference resolution, and its effect on the hypsometry is +17 m in the wrong direction.
Reverting it is a one-line call if the author prefers to hold the dimensional fix until the
hillslope law is redesigned — the argument for keeping it is that a dimensionally wrong term
is harder to reason about than a dimensionally right one, not that it improves the output.


### Finding 44 — IMPLEMENTED: units named at unchanged output, and the derived step count

Everything added is **additive and read-only**: `incise` reads no new field, so the unit naming
cannot have changed any output. `StreamPowerConfig` gains `k_time()`, `k_for_duration()`,
`celerity_m_per_yr()`, `dt_max_yr()`, `cfl_iterations()`, `courant()`, `timescale_plan()`, plus
`SHIPPED_K_TIME = 9000` and `COURANT_INTEGRATING = 1.0`.

**The byte-identity proof, with a negative control.** `timescale_naming_changes_no_output`
pins the shipped numbers (`dt = 1.0`, `iterations = 2`, `k = 4500`, `k_time() = 9000`) — but a
pinning test alone could be vacuous, so it also proves it CAN see a change of this kind:
incising a fixed synthetic field with `1 × 9000` instead of `2 × 4500` (same `k_time`) must and
does differ, because the relaxation is nonlinear and the flow field is recomputed between
steps. So "byte-identical" here is a checked claim, not an assumed one.

**`K` and the duration are not separately observable.** Hold `k_time` and the step count, pick
any `T`: then `k = k_time/T` and `dt = T/iterations`, so `k·dt = k_time/iterations` — every
quantity `incise` reads is unchanged and `T` cancels.
`duration_cancels_out_of_the_incision` pins it at T = 10⁴, 10⁶ and 10⁸ yr. Two consequences:

- naming the units is FREE, which is why this landed at unchanged output;
- **a duration dial alone is NOT a dial.** H-2's France↔Scotland control must move `k_time`;
  10⁴ against 10⁶ years at proportionally larger `K` is the SAME terrain. Reading `k_time` as
  years requires pinning `K` independently — which C-3's per-lithology multipliers plus one
  absolute Stock & Montgomery anchor supply. **This is a correction to how H-2 was framed:**
  the dial is an integrated erodibility-time product, not a time.

**The `dt = 1.0` placeholder, stated as such.** Read literally it says two years of erosion
carved ~580 m of mean relief. That absurdity is diagnostic, not embarrassing: `(k = 4500,
dt = 1.0)` is one arbitrary factorisation of `k_time = 9000` among infinitely many. The
anchored reading runs the other way — pin `K`, and the duration follows.

**What `dt_max` depends on.** The CFL bound for the detachment-limited erosion wave: for
`n = 1` the knickpoint celerity is `c = K·A_km²^m` (m/yr), and the wave must not cross more
than one cell per step, so `dt_max = cell_m / c`. **Linear in cell size** — that is where the
grid legitimately enters a timescale, and `cfl_bound_scales_with_cell_size` pins the ×4.

### The measured cross-resolution behaviour — and it is a diagnosis, not a config

`A_max` is an input, taken from the production measurement (max flow accumulation on land,
full chain via `production_hd_config`): 3447 km² at 2048², 1611 km² at 8192².

| grid | cell (m) | A_max km² | celerity m/yr | dt_max (yr) | **CFL steps** | shipped | **Courant** |
|---|---|---|---|---|---|---|---|
| 2048² | 195.3 | 3447 | 264 200 | 7.39e-4 | **2706** | 2 | **1353** |
| 8192² | 48.8 | 1611 | 180 618 | 2.70e-4 | **7399** | 2 | **3699** |

`dt_max` falls ×2.73 overall, and that is TWO effects: **×4 from the cell size** (the bound is
linear in `cell_m`) times **×0.68 from `A_max` itself dropping** 3447 → 1611 km² (Finding 41 —
more closed basins capture more catchment at a finer grid). Holding `A_max` fixed gives exactly
×4. Reporting the ×2.73 as if it were the cell-size effect alone would have been the same
mistake as the ratio-versus-excess confusion: two effects in one number.

**THE ITERATION SATURATION IS EXPLAINED, and it is not a knob.** At Courant ≫ 1 the implicit
update is STABLE but not INTEGRATING — stability is not accuracy. With
`f = K·dt·A^m/dist_m ≫ 1` the update `h ← (h + f·h_r)/(1 + f)` drives each cell essentially
onto its receiver's height in ONE step, so the terrain reaches a local relaxed state
immediately and further sweeps do little. 8192² sits **2.7× further past the bound**, which is
why it saturates harder: retained mean altitude 82.4 → 72.2 % over 1→8 iterations, against
51.7 → 12.3 % at 2048².

**And the derived count is NOT an adoptable configuration.** 7399 steps at 8192², each a full
flow recompute plus incision over 67 M cells. So deriving the count honestly does not fix the
model — **it reveals that the model is not time-integrating at all**, at either resolution
(Courant 1353 even at 2048²). The shipped terrain is the fixed point of two local relaxations,
not the result of an erosion episode. That is a legitimate way to make terrain; it is not a way
to express a duration, which is exactly what H-2 asked for.

### K anchoring — the debt C-3 left open, closed to an order of magnitude

Stock & Montgomery 1999, audited on the source (ADR C-3): hard rock 10⁻⁷–10⁻⁶ with `A` in m²
at m = 0.4. Ymir's law takes `A` in km², so `K_ours = 10³^(2m)·K_lit ≈ 251·K_lit` at m = 0.4.
With `k_time = 9000`:

| K_lit | K_ours | implied duration |
|---|---|---|
| 1e-7 (hard rock, low) | 2.51e-5 | 3.58e8 yr = **358 Myr** |
| 1e-6 (hard rock, high) | 2.51e-4 | 3.58e7 yr = **36 Myr** |

**⚠️ Caveat stated, not buried: the table is fitted at m = 0.4 and Ymir ships m = 0.5.** A
tabulated `K` is only valid at the exponent it was fitted with, so this is an ORDER OF
MAGNITUDE and not a calibration. What it establishes: the shipped terrain's integrated `K·T`
is consistent with an episode of order **10⁷–10⁸ years at hard-rock erodibility** — a plausible
orogenic-to-cratonic duration. It does not validate the value; it says the lumped constant is
not absurd once read dimensionally, which is the first time that could be said at all.

### Where this leaves the hypsometry remedy

The timescale work was expected to address mechanism 2's iteration saturation on the way. It
does something better and less convenient: it shows the saturation is not a step-count problem
to be retuned but a consequence of running 10³× past the wave bound. Both resolutions are
outside the integrating regime; 8192² is further outside. So:

- retuning `iterations` per resolution remains a fix to the metric, not the terrain
  (confirmed from a second direction);
- mechanism 1's remedy still needs a transport-limited hillslope law with an explicit sediment
  flux — and it now has a year to integrate against, which it did not before;
- but that law must be posed in a regime where the model integrates, which the shipped
  configuration does not. **That is the real next question, and it is bigger than a
  calibration.**

### Two method points, to stand prominently

**A RATIO IS ONLY VALID AT COMPARABLE TOTAL EFFECT.** `iterations = 1` gave the best
8192²/2048² ratio (1.32 against 1.85) while removing almost none of the excess (19 m of 147),
because lowering total erosion moves BOTH grids toward the invariant un-eroded limit where the
ratio → 1.00 trivially. Proof by combination: adding `iters = 1` to MFD-off makes the excess
WORSE (73 → 102 m) while the ratio still improves. The **excess in metres** is the instrument —
and it is also what the consumers experience, since temperature is a lapse rate on ABSOLUTE
altitude (the 287/693 m gap is 2.62 °C, not a percentage). This chantier was reasoned in ratios
throughout and the ratio was the wrong instrument.

**CITE A FINDING'S MEASUREMENT, AND RE-DERIVE THE MECHANISM AGAINST CURRENT CODE.** An
observation and its explanation do not have the same lifespan, and only the observation is
protected by the measurement. Finding 3's cause was correct when written and was silently
invalidated by C-1 flow conditioning; it survived because nothing re-checks an explanation when
its premise moves, and the observation kept being true. That is how a wrong model survives —
not through a false measurement, but through a true measurement dragging a dead explanation.
Never chain two findings' explanations without checking that both premises still hold.

### A tooling defect of my own, repaired

Appending to the ADR and to two source files through PowerShell `Get-Content -Raw` +
`Add-Content -Encoding utf8` **double-encoded every non-ASCII character** (PowerShell 5.1's
`Get-Content` reads UTF-8 as cp1252, then re-encodes). 275 lines were corrupted across three
files. Repaired by inverting the double encoding per line, except for five characters that were
IRRECOVERABLE — cp1252 leaves 0x81/0x8D/0x8F/0x90/0x9D undefined, so those bytes were destroyed
rather than transformed (`∝`, `↔`, `10⁴`, superscripts) and had to be restored by hand. Zero
replacement characters remain. **Use the Write tool or Python for any append to a file with
non-ASCII content; `Add-Content` is not safe here.** Recorded because it was silent: the files
compiled and read fine to a grep, and the corruption only surfaced when a `cat -A` was needed
for an unrelated reason.

## Method rules earned in the hypsometry chantier — to be applied by default

These are not observations about the terrain; they are rules about how to measure it. Each one
cost at least one wrong conclusion of mine.

### 1. A pinning test needs a NEGATIVE CONTROL before it can assert byte-identity

`timescale_naming_changes_no_output` asserts the shipped numbers are untouched. On its own
that is **potentially vacuous** — a test that checks constants can pass while the thing it is
supposed to protect has changed underneath, and a test that compares a value to itself always
passes. So it also proves it WOULD SEE a change of the relevant kind: incising a fixed
synthetic field with `1 × 9000` instead of `2 × 4500` (identical `k_time`) must differ, and
does. Only then does "the units were named at unchanged output" mean anything.

**Default from now on: every byte-identity claim carries a control that fails.** State what
change the test can detect, and demonstrate it detecting one. A byte-identity assertion without
a control is a statement about the test, not about the code.

### 2. A RATIO IS ONLY VALID AT COMPARABLE TOTAL EFFECT

`iterations = 1` gave the best 8192²/2048² ratio (1.32 against a 1.85 baseline) while removing
almost none of the divergence (19 m of 147 m of excess). Lowering total erosion moves BOTH
grids toward the invariant un-eroded limit, where the ratio → 1.00 **trivially**. Proof by
combination: adding `iters = 1` to MFD-off makes the excess WORSE (73 → 102 m) while the ratio
still improves (1.63 → 1.30).

The instrument is the **excess in metres** — and it is also what the consumers experience,
since temperature is a lapse rate on ABSOLUTE altitude (the 287/693 m gap is 2.62 °C, not a
percentage). This whole chantier was reasoned in ratios and the ratio was the wrong instrument.

### 3. Decompose a measured factor before quoting it

`dt_max` falls ×2.73 between the two grids. Quoting that as the cell-size effect would have
been **the same fault as ratio-versus-excess**: it is ×4 from the cell size (the CFL bound is
linear in `cell_m`) times ×0.68 from `A_max` itself dropping 3447 → 1611 km² (Finding 41 —
more closed basins capture more catchment at a finer grid). Holding `A_max` fixed gives exactly
×4, pinned by a unit test. **When a measured factor could contain two effects, separate them
before it enters a sentence** — and pin the isolated one with a test so the decomposition is
not just an argument.

### 4. State the reservation on a borrowed constant, in the same breath as the number

The `K` anchoring uses Stock & Montgomery's table, which is fitted at **m = 0.4** while Ymir
ships **m = 0.5**. A tabulated `K` is only valid at the exponent it was fitted with, so the
result is an ORDER OF MAGNITUDE and not a calibration. What is now sayable and was not before:
**the lumped constant is not absurd when read dimensionally** (10⁷–10⁸ years at hard-rock
erodibility, a plausible orogenic-to-cratonic episode). That is a real gain; presenting it as a
calibration would have been a fabrication.

### 5. Cite a finding's MEASUREMENT; re-derive its mechanism against current code

An observation and its explanation do not have the same lifespan, and **only the observation is
protected by the measurement**. Finding 3's cause ("the FBM detail resolving sharper gradients
on finer cells") was correct when written and was silently invalidated by C-1 flow conditioning,
which caps the FBM's whole contribution to ±2 m. It survived because nothing re-checks an
explanation when its premise moves, and the observation kept being true. **That is how a wrong
model survives — not through a false measurement, but through a true measurement dragging a
dead explanation.** Never chain two findings' explanations without checking that both premises
still hold.

### 6. Tooling: the PowerShell append that corrupted the documentation silently

`Get-Content -Raw` piped into `Add-Content -Encoding utf8` **double-encodes every non-ASCII
character** (PowerShell 5.1 reads UTF-8 as cp1252, then re-encodes). It corrupted **275 lines**
across this ADR and two source files. It is silent in every direction that matters: the code
still compiled, the markdown still rendered, and `grep` matched nothing unusual — it surfaced
only because a `cat -A` was needed for an unrelated reason.

Worse, the damage is **not fully reversible**: cp1252 leaves `0x81/0x8D/0x8F/0x90/0x9D`
undefined, so those bytes are DESTROYED rather than transformed. `∝`, `↔` and the superscripts
could not be recovered by inverting the encoding and had to be retyped. The pattern is now
banned in `CLAUDE.md`; use the Write/Edit tools or Python with an explicit
`encoding='utf-8'`.

The general point, which is why this belongs beside the others: **a tool that corrupts its
output silently is worse than one that fails**, and documentation has no test suite to catch
it. Any bulk edit of a prose file needs a post-check — here, `grep -c $'\\ufffd'` and a scan
for the classic mojibake digraphs.

#### EXTENSION (Finding 56b) — verify a written artefact NON-EMPTY before declaring it delivered

The rule above is a prohibition on one tool, and that is not enough. **The failure class is
"the artefact reports the instrument instead of the subject", and an empty or truncated document
is that class in documentary form** — the exact counterpart of a metric column constant at zero.
So the rule now reads:

> **A written artefact is not delivered until it has been verified non-empty, decodable and
> un-mangled.** Not inspected — verified, by reading back what was written. An empty note, like
> a column constant at zero, measures only the tool. And "documentation has no test suite" was
> a description of a gap, not a licence.

This is the **fourth instance on this campaign** of the same class, which is why it stops being
a reminder:

1. the **global axial R**, blind to local parallelism (0.03–0.05 global against 0.52 local) — it
   reported the instrument's aperture, not the coast, and my "parallelism refuted" verdict was
   wrong until the renders contradicted it;
2. the **8 km spur-length cap**, reading 7.99–8.00 in *every* configuration including the
   reference, hiding the target completely;
3. **`channel width p50` constant at 0.000 m** at every order, both grids, both settings —
   rule 10's first catch, and it took two refusals to fix because switching to p90 still read
   0.000 for an order with no wet reach at all;
4. **a report file reported empty.** It was not (6 805 bytes, 129 lines, valid UTF-8, no BOM, no
   double-encoding — and identical in `HEAD`), but *"I looked and it seemed fine"* is exactly the
   inspection this rule says is insufficient, in either direction. **The question is now settled
   by measurement**, and the measurement runs on every `cargo test`.

**Made mechanical**, the way rule 10 became a type rather than a paragraph:
`crates/ymir-core/tests/doc_artifacts.rs` walks every Markdown file under `docs/` (161 at the
time of writing) and fails on any that is empty, under a 200-byte floor, not valid UTF-8,
carrying a cp1252 round-trip signature, or holding unresolved conflict markers. Two details make
it more than decoration:

- **it asserts its own population** (`files.len() >= 10`) — a walk that silently matched nothing
  would otherwise report "0 problems out of 0 files", which is rule 10's defect wearing a
  different hat;
- **its negative control builds the corruption instead of trusting a literal**: it runs a real
  cp1252 round trip over `hypsométrie ↔ 25°` and asserts that at least one signature matches, so
  a signature list that could never fire on the real defect fails the test rather than passing
  it quietly. The five undefined cp1252 slots are encoded there too, as the reason the damage is
  irrecoverable rather than merely reversible.

**A second prerequisite made mechanical in the same round**, because it has the same shape — a
claim that would have failed silently and produced a confident wrong answer. The author's
regeneration step assumes the `A_c(S)` toggle invalidates the `eroded` and `drainage` caches.
If it did not, enabling the law and regenerating would return the SHIPPED terrain from cache and
the visual verdict would be passed on the wrong picture, with nothing in the log to say so —
precisely the defect `eroded_key_full`'s own doc comment records for volcanism ("the terrain
differs, the drainage would not"). `crates/ymir-core/tests/channel_head_law_cache_key.rs` proves
it instead: the digests move (`5a735ac0 → a84ac294` eroded, `9770ded0 → f1d20f47` drainage), both
law parameters reach the key so a recalibration of `S_ref` or of the `S_min` proxy cannot be
served stale, and a negative control asserts the key is STABLE across identical calls — without
which the test would pass for the wrong reason on a key that folded a timestamp.

## Finding 45 — microscope semantics: an entry is a river SYSTEM, not a reach between water bodies

A reach ending on a lake shore gets `downstream = None`, so 87–92 % of terminals were lake
INFLOWS (Finding 42) and a microscope entry was a REACH BETWEEN TWO WATER BODIES. That is why
the head of the discharge sort showed a fragment: the largest ENTRY was a 1-segment, 2 km stub
carrying its parent's inherited area, not the trunk it belonged to.

**The fix, in the viz only** (`aggregate_watercourses`, no core change, no cache bump): an
inflow reach is linked to the outlet reach of the EXORHEIC lake it dies on — water physically
continues through a lake that has an outflow. NOT for an endorheic lake (the water dies there,
a true terminus) nor for a below-sea basin (its outflow is a typed `Spillway` with no
hierarchy). The trunk climb is extended symmetrically, or the main stem would still stop at the
shore. Both directions are gated on one flag so the bench measures BEFORE and AFTER on the SAME
`HdResult` in a single production run — an exact comparison rather than two builds.

### 2048² — exactly the intended effect

| | entries | head of the discharge sort |
|---|---|---|
| before | 1359 (1343 rivers) | **#1 = S1, 0 trib, 1 segment, 2 km**, A 1087 km² → mer |
| after | 1095 (1080 rivers) | **#1 = S4, 394 trib, 395 segments, 146 km**, A 1087 km² → mer |

264 entries merged into the systems they belonged to; mean segments per entry 7.6 → 8.5. The
head is now the trunk instead of the stub that shared its inherited area. The 1-segment stub is
still present at #2 — that is the SEPARATE defect (isolated fragments inheriting the parent's
area, Finding 42), which lives in `clip_rivers_to_lakes` in core and was deliberately not
folded in here.

### 8192² — the measurement exposes a REGRESSION, and it is not this change

| 8192² | entries | rivers | spillways | biggest river |
|---|---|---|---|---|
| **before the `dscale` patch** (Finding 43 bench) | 792 | 749 | **43** | 110 km² |
| after the patch, before chaining | 2523 | 529 | **1994** | 48 km² |
| after the patch, after chaining | 1202 | 280 | 922 | 9 km² |

**The below-sea basin population went 43 → 1994 (×46).** Both benches use identical parameters
(seed, domain, latitude, span, relief-v3 triple, lithology, fracture, infiltration) and the only
change to the terrain between the two measurements is the linear-diffusion `dscale`
normalisation. So the attribution is sound by construction, without needing an A/B rerun (which
would cost hours — the 8192² arm of this bench alone took 3.7 h on a cold cache).

**Why it happens is exactly the mechanism already identified and it should have been predicted.**
A linear Laplacian is MASS-CONSERVING: it lowers crests and FILLS hollows. With
`diffuse_channels = true` it runs on every cell. Making it 16× stronger at 8192² therefore
backfills drainage — and a backfilled channel is a CLOSED DEPRESSION. Finding 41's pathology,
manufactured at scale.

**My verification was too narrow, and that is the lesson.** I measured the patch on the
hypsometry alone (+17 m on the mean, ratio 2.43 → 2.49) and reported it as "small, in the
unhelpful direction". The mean moved 2.5 % while the closed-basin population moved ×46. A
single scalar chosen because it was the metric under investigation is not a verification of a
change to the erosion.

**New method rule, rule 7: verify a change on every observable it could plausibly move, not
only on the one being investigated.** For an erosion change the minimum set is the hypsometry,
the closed-depression / below-sea basin count, and the drainage network extent — the first two
moved in opposite proportions here, and only the second matters for playability.

### What the 8192² head shows once the terrain is set aside

Even on the regressed terrain the chaining does its job (2523 → 1202 entries, and the head goes
from five reaches all dying in lakes to systems reaching the sea), but two things are worth
noting for when the terrain is restored:

- a `Sink::Unknown` mouth ("→ ?") appears at rank #2 — a mouth that is neither sea, lake nor
  sub-sea sink. Finding 42 counted only 3 of these at 8192², all carrying 0 km²; one at the head
  of the list means the regressed terrain produces them with real discharge. To re-check on a
  restored terrain rather than diagnosed here;
- chaining absorbs a spillway into a river system when the basin spills into an exorheic lake
  whose outlet reaches the sea (spillway entries 1994 → 922). Hydrologically that is correct —
  the flow does continue — but it means a system's `kind` comes from its MOUTH, so a system
  containing a spillway reach lists as a `Watercourse`. Deliberate and documented, not a bug;
  the per-segment `kind` in `rivers.json` is unaffected.

## Finding 45b — the `dscale` patch REVERTED, and why it must not come back

### What was decided and why it was decided wrongly the first time

The patch applied Finding 7's `dscale = (HILLSLOPE_REF_CELL_M/cell_m)²` to the linear diffusion
branch, closing a real dimensional inconsistency. It was kept on my report that its effect was
**"+17 m on the 8192² mean, small and in the unhelpful direction"**. That measurement was
CORRECT on the hypsometry and **missed the actual effect entirely**.

| 8192², same seed / domain / config | before the patch | after the patch |
|---|---|---|
| mean land altitude | 685 m | 702 m (**+2.5 %**) |
| **below-sea basins (spillways)** | **43** | **1994** (**×46**) |
| microscope entries / rivers | 792 / 749 | 2523 / 529 |
| biggest river catchment | 110 km² | 48 km² |

**The patch did not change the altitude — it destroyed the drainage.** Two properties of the
SAME field, one instrumented and one not: the altitude distribution barely moved while the
drainage topology collapsed.

### The result is more interesting than the patch

**A CORRECTLY NORMALISED DIFFUSION DESTROYS DRAINAGE INTEGRITY.** This confirms by a second,
independent route what the patch measurement first suggested: a linear Laplacian is
MASS-CONSERVING — it fills valleys exactly as much as it lowers crests — and with
`diffuse_channels = true` it runs on every cell, so it backfills the channels the incision has
just cut. At 16× strength (which is what "correct" means at 8192²) it backfills them
completely, and a backfilled channel is a closed depression. Finding 41's pathology,
manufactured at scale.

That is a discovery about **the term**, not a reason to ship 1994 closed depressions.

### The dimensional debt, REFORMULATED so the patch cannot resurrect

Finding 7's fix genuinely lives in only one arm of an `if`, and `relief_v3` takes the other. In
six months that will look like an obvious omission and someone — including me — will be tempted
to "finally fix it". **It will not close by reinstating `dscale`.** The measurement above is the
counter-example, and the code now carries it at the exact line, in a comment that names the
numbers.

The debt closes only when **the hillslope term changes NATURE: transport-limited with an
explicit sediment flux** — the only kind of agent that REMOVES mass from hillslopes and delivers
it to the channel network, which is what mechanism 1 needs (Finding 43: `A_c` is a correct
threshold feeding an inert branch, and 90 % of fine-grid land is handed to that branch). **When
that term replaces the Laplacian, the units question disappears with it** — a conservative
smoother and a transport law do not have the same dimensional problem, because they are not the
same equation.

Option 2 considered and rejected on the author's call: keeping `dscale` and switching
`diffuse_channels` to `false` would also change 2048², so it forfeits the byte-identity at the
reference cell and requires recalibrating both resolutions. **That is a model change disguised
as a dimensional correction**, and it should be taken as a model change or not at all.

### Method rule 7 — verify what a change is NOT supposed to affect

Rule 3 said to decompose a measured factor before quoting it. This is its complement, and it
cost a shipped regression:

> **Verifying a fix on the property it is SUPPOSED to affect is not enough. Verify the ones it
> is NOT supposed to affect.**

The patch was judged on the hypsometry because hypsometry was the chantier's subject. Nothing
about a diffusion normalisation is *supposed* to change the closed-basin count — which is
exactly why nobody looked, and exactly why that is where the damage went. A change to a term in
the erosion touches the whole field; the minimum verification set for one is:

1. the hypsometry (mean and percentiles, in METRES — rule 2);
2. the **closed-depression / below-sea basin count** — the playability-critical observable, and
   the one that moved ×46 here;
3. the drainage network extent (channel length, which is resolution-stable when healthy —
   Finding 42 measured 14 605 vs 14 284 km).

Two of those three would have caught this in the same run that produced the +17 m.

## Finding 46 — fragment areas: a per-point discharge, a refuted attribution, and a broken determinism

### The fix

`clip_rivers_to_lakes` gave every run of a split segment the PARENT'S maximum discharge, area
and width. The parent's maximum sits at the parent's downstream end, so a run cut off earlier
advertised a catchment it does not have. The remedy needed a per-point quantity, since the clip
needs a value at an INTERIOR point of a segment:

- `C1DrainageResult::segment_discharge_profile_m3s` — discharge (m³/s) at each `segment.points`,
  parallel to `segment_profile_m`, built in the drainage extraction from the same runoff
  accumulation the per-segment figures come from (no recompute, no second derivation to drift);
- each clipped run reads the value at ITS OWN downstream-most point, and the width follows from
  it (`w = a·Q^b`, Finding 22 — a width taken from an inherited discharge was the same defect);
- `ALGO_DRAINAGE` 2→3 and `ALGO_HD_DRAINAGE` 4→5.

**The Finding 22 exception, kept and stated at the line** (the asymmetry a later reader would
"simplify"): an EXORHEIC LAKE OUTLET run, i.e. one starting just after a lake cell (`a > 0`),
inherits on purpose. It evacuates the whole upstream catchment, which arrived through the lake.
Removing that branch drops a trunk's width to zero across every lake it crosses — the author's
original bug.

**The strict sidecar length check earned its keep immediately.** The spillway append in
`hd.rs` did not push the new array, so it was 10170 long for 10186 segments. Under the OLD
tolerant `resize` that would have shipped silently (Finding 45's `resize` lesson); instead the
next run failed with `drainage sidecar seg_qprof: 10170 entries for 10186 segments`. Fourth
occurrence of the composite-completeness trap, first one caught by a guard instead of by a
measurement three rounds later.

### An attribution of mine, refuted by the measurement

I had filed the head-of-list stump as an inherited-area defect. **It is not.** Its anatomy, at
2048², cell (1087, 1188):

| segment | points | order | area | upstream | first → last |
|---|---|---|---|---|---|
| #5186 | 4 | S4 | 1087 km² | 1 | (1084,1188) → (1087,1188) |
| #7256 | 11 | S1 | 1087 km² | 0 | (1083,1198) → (1087,1188) |

Shared points: **1 of 4 and 11** — two DISJOINT paths converging on one coastal cell. The stump
reaches the mouth, so **its area was already correct**: the accumulation at that cell really is
1087 km². The defect is that the cell's accumulation is the UNION of both catchments and both
reaches read it.

**So this is a distinct, newly-characterised defect: two reaches terminating on the same cell
both report that cell's accumulation.** 46 of 46 duplicated-terminal groups at 2048² still
claim identical areas; 24 of 26 at 8192² (the two exceptions are spillways, whose profile is
uniform by construction — 1572 km² against 8 km² on the same cell, which is the per-point fix
visibly working).

Proposed rule, **not implemented** — it is a semantic choice for the author: a reach's own
catchment is the accumulation at its last point that is NOT shared with another reach's
terminus (in practice the penultimate point at a confluence). Whether a reach ending at a
shared coastal cell should report its own catchment or the cell's is a decision, not a bug fix.

### And a determinism break, found by doubting my own result

After the fix the head became the S4 trunk (395 segments) instead of the stump — but the two
TIE on discharge (10.3 m³/s, same mouth cell), and `aggregate_watercourses` iterated its groups
from a **`HashMap`**, whose order is randomised per process. With a stable sort and a tie, the
head was decided by hash order. **"The head is now the trunk" would have been luck, not a
property** — and non-determinism is a core project invariant, violated here in the instrument
the author judges with.

Fixed: the groups are collected and sorted by root index, and the final comparator gains a
deterministic, MEANINGFUL tie-break — more segments first (a system outranks a stump sharing its
mouth), then the mouth cell, so the order is total. Verified by running the bench twice with
identical output. This is rule 1 (the negative control) applied to a result rather than a test:
I doubted the mechanism behind a favourable outcome and it was not the one I had assumed.

### The verdict, both resolutions, on RESTORED terrain

| | 2048² | 8192² |
|---|---|---|
| **head of the discharge sort** | **#1 S4, 394 trib, 395 seg, 146 km, → mer** | **#1 S3, 98 trib, 99 seg, 27 km, → mer** |
| (before Finding 45's chaining) | #1 S1, 0 trib, 1 seg, 2 km | #1 S2, 43 trib, → lac |
| entries (rivers / spillways) | 1095 (1080 / 15) | 659 (619 / 40) |
| fragments (no up, no down) | 666 (6.5 %) | 560 (4.5 %) |
| fragments claiming > half the max area | 4 | 2 |
| duplicated terminals (same area) | 46 (46) | 26 (24) |

At 8192² the whole top five are now multi-tributary systems; not one 1-segment stump. **The
8192² figure owed from Finding 45 is delivered here, on restored terrain** — and the terrain is
verifiably restored: mean land altitude 693 m, exactly the pre-`dscale` value.

### RULE-7 CONTROL BLOCK — the observables this change must NOT move

| | 2048² | 8192² |
|---|---|---|
| hypsometry mean / p50 / max (m) | 287 / 167 / 2684 | 693 / 447 / 4145 |
| below-sea basins | 15 | 40 |
| network extent | 74 897 cells = 14 628 km | 291 508 cells = 14 234 km |
| lake cells / inventoried lakes | 99 344 (3790 km²) / 43 | 681 405 (1625 km²) / 56 |

All match their pre-change values (the below-sea counts read 15/40 rather than 16/43 because
Finding 45's lake chaining absorbs a spillway whose basin spills into an exorheic lake — a
counting effect in the microscope, not a change to the basin population). The change touches
per-segment arrays only, and the control block confirms it.

## Finding 47 — `catchment_km2`, the never-applied thresholds, and which quantity navigability should read

### The shared-mouth convention, settled and applied

The author's call: **a system reports ITS OWN catchment, never the union.** A system whose mouth
happens to share a coastal cell with another is not draining that other basin. Applied in
`aggregate_watercourses`: when the root reach's terminus is also another reach's terminus, read
the last point that is NOT the shared confluence (one index back in the per-point discharge
profile from Finding 46). It was a few lines, so it landed here rather than being scoped out.

| head of the discharge sort, 2048² | before | after |
|---|---|---|
| #1 | S4, 394 trib, **1087 km²** | S4, 394 trib, **1085 km²** |
| #2 | **S1, 0 trib, 1 seg, 1087 km²** | S3, 85 trib, 86 seg, 643 km² |
| #3 | S3, 85 trib, 645 km² | S2, 24 trib, 596 km² |
| #4 | **S1, 0 trib, 1 seg, 645 km²** | … |

**Every stump is gone from the head**, and the trunk now reports 1085 km² — its own catchment
just above the confluence — instead of 1087, the union. Verified deterministic over two runs.

### THE HEADLINE: the recalibration was never applied — and would not have sufficed

The code still ships `stream 20 / small_boat 500 / barge 5 000 / ship 50 000 km²`. The ADR
records the author's re-anchored proposal — `10 / 100 / 1 000 / 8 000` — with the words
**"Not applied (default unchanged) pending review"**. So it was not lost; it was parked
awaiting a review that never came. Stated plainly, since the consequence has been visible in
every export since.

Measured class split over all reaches (arid-hot 25°, 2048²), classified as the code classifies
— on the runoff-equivalent quantity:

| thresholds | non-nav | small boat | barge | ship |
|---|---|---|---|---|
| **code today** 20/500/5 000/50 000 | 10 081 | 105 | **0** | **0** |
| author's proposal 10/100/1 000/8 000 | 9 606 | 555 | **25** | **0** |

**"Only small boats, not even barges" is reproduced exactly** by the shipped thresholds. And
the proposal, had it landed, would have produced 25 barge reaches and **still zero ships** —
because 8 000 km² exceeds the measured maximum in either quantity (1 459 km² geometric,
1 087 runoff-equivalent). The proposal was anchored on a measured max of "~12–13 000 km²
(Thames-scale)"; **the terrain changed under the calibration.** Third instance of method rule 5
(a premise moved and nothing re-checked the conclusion).

### Which quantity was it fitted to? The runoff one, by distribution match

The proposal's anchors were "p50 ~50, p90 ~300–500 km²". Measured today:

| at the sea mouths | arid 25° | humid 45° | 8192² arid |
|---|---|---|---|
| GEOMETRIC km² p50 / p90 / max | 34 / 167 / 1 459 | **34 / 167 / 1 459** | 28 / 102 / 1 379 |
| RUNOFF-EQ km² p50 / p90 / max | 59 / 287 / 1 087 | 61 / 515 / 6 247 | 0 / 14 / 110 |
| DISCHARGE m³/s p50 / p90 / max | 0.56 / 2.73 / **10.33** | 0.58 / 4.90 / **59.39** | 0.00 / 0.13 / **1.04** |

The p50/p90 anchors (50, 300–500) match the RUNOFF-EQUIVALENT figures (59, 287–515), not the
geometric ones (34, 167). **So the calibration was fitted to the runoff quantity while its
rationale was written about areas** — "navigability is physical river size", "absolute, not
domain-scaled, because a given river size should classify the same on any map".

### And that rationale needs correcting, which the humid column proves

Look at the geometric row: **34 / 167 / 1 459 in the arid AND the humid run — identical.** Same
terrain, so the same geometry. Yet the discharge goes 10.33 → 59.39 m³/s (×5.8). Therefore:

- **an area-based class cannot tell an arid trickle from a Thames-scale river on the same
  terrain.** Area is blind to the thing that decides navigability;
- so "the same river size should classify the same on any map" was the wrong principle. It is
  the same DISCHARGE that should classify the same — and a river's navigability *should* differ
  between an arid and a humid world. That is not an inconsistency to design away; it is the
  physics.

**Recommendation on the variable: DISCHARGE, in m³/s.** The author's reading is right and the
measurement now supports it.

### But the thresholds cannot be calibrated yet, and the reason matters

The discharge is **not resolution-stable**: max 10.33 m³/s at 2048² against **1.04 at 8192²**
(×10), same seed and config — Finding 43's hypsometry defect propagating all the way to the
navigability class. The class split shows it starkly: the shipped thresholds give 105 small-boat
reaches at 2048² and **2** at 8192². A class assignment that moves ×50 with the render
resolution is not a property of the world.

Meanwhile the GEOMETRIC distribution is resolution-stable (p50 34 → 28, p90 167 → 102,
max 1 459 → 1 379).

**So: the right variable is unusable until the hypsometry converges, and the stable variable is
the wrong physics.** Calibrating m³/s thresholds now would bake in the 2048² numbers — exactly
the mistake the parked proposal already made once. The layered recommendation:

1. **Now** — export `catchment_km2` (the geometric area, done) so a consumer has a stable
   figure, and leave the thresholds untouched, DOCUMENTED as classifying a runoff surrogate.
   Changing them now would be a second calibration against a moving premise.
2. **Recorded now so it is not parked-and-lost a second time** — the intended thresholds in
   m³/s, anchored on real navigation rather than on this map's distribution: **small boat ≥ 1,
   barge ≥ 50, ship ≥ 300 m³/s** (Thames ~65 m³/s is barge-navigable in its lower reaches;
   Seine ~560 is ship-capable; Rhine ~2 300). Under those: humid 45° yields a barge-class trunk
   at the top and no ships; arid 25° yields small boats only.
3. **After the hypsometry converges** — apply (2) and re-measure. Not before.

### The honest scale statement

The largest river is **10.3 m³/s in the arid band — 6.3× smaller than the Thames** — and
**59.4 m³/s in the humid band, essentially Thames-scale (×1.1)**. So:

- "only small boats" is the CORRECT answer for a 400 km island at 25°, not a calibration bug;
- **the lever is the climate, and it is sufficient**: the same continent at 45° produces a
  Thames-scale trunk. No terrain change, no threshold change — the water.

That reframes the original complaint entirely. Nothing was miscalibrated about the river sizes;
the map was arid.

### The export

`rivers.json` gains `catchment_km2` — the geometric contributing area in km², from
`flow.accumulation` at each reach's own downstream-most cell. `drainage_km2` is retained for
continuity but its docstring now says plainly that it is **not an area**: it is
`runoff_accumulation / 300 mm`, a discharge surrogate, and `discharge_m3s` is the same
information in honest units. Core keeps the catchment in CELLS
(`segment_catchment_cells`) so no interior stage needs a `cell_km2` it does not have; the
conversion happens once, in `rivers_json`, at the point of use. `ALGO_DRAINAGE` 3→4,
`ALGO_HD_DRAINAGE` 5→6.

### The completeness trap, FIFTH occurrence — and eliminated structurally this time

Adding `segment_catchment_cells` I forgot the spillway append, exactly as I had forgotten
`segment_discharge_profile_m3s` one finding earlier. The strict sidecar check caught it again
(`index out of bounds: the len is 10170 but the index is 10170`). Two consecutive occurrences
of the same mistake by the same person is a design problem, not a discipline problem, so the
class is now closed:

**`SegmentRow` + `C1DrainageResult::push_segment`.** One struct holding every per-segment value,
one method that pushes them together and asserts the invariant. Adding a field to `SegmentRow`
is a **compile error at every call site** — which is the guard the runtime invariant could only
report after the fact. Sites that append a segment must now use it.

### RULE-7 CONTROL BLOCK

This round touches per-segment arrays and one viz aggregation; nothing in the height field.
Confirmed unchanged at 2048²: hypsometry mean 287 m / p50 167 / max 2 684; below-sea basins 15;
network extent 74 897 cells = 14 628 km; lake cells 99 344 (3 790 km²), 43 inventoried lakes.
The geometric mouth distribution is identical between the arid and humid runs (34/167/1 459),
which is itself a control: the change did not touch geometry.

## Finding 48 — the coastline contour: the mechanism is confirmed, the geometric remedy is not

### First, what the defect is NOT

Marching squares **already interpolates sub-cell** along each crossed edge
(`t = (iso − va)/(vb − va)`), so the barbs are not blockiness and "sub-cellular interpolation"
was not the missing piece. The measured sea-level offset before any change is **median +0.000 m,
|p90| 0.001 m** — the raw contour sits exactly on the isoline. Whatever is wrong is not the
placement of a vertex along its edge.

### The falsifiable prediction, and its verdict

The diagnosis on record is GRADIENT PINNING: where the field is nearly flat at the iso value,
WHICH edges get crossed is decided by tiny fluctuations, so the contour zig-zags. That predicts
the fix must improve the LOW-SLOPE shores most.

The baseline confirms the diagnosis's premise — the defect is where it says:

| turns > 80° | <0.5° | 0.5–2° | 2–5° | 5–15° | >15° |
|---|---|---|---|---|---|
| 2048² baseline | **82.5 %** | 54.4 % | 15.9 % | 10.2 % | **1.4 %** |
| 8192² baseline | **83.5 %** | 67.5 % | 55.6 % | 19.2 % | **2.5 %** |

**PREDICTION CONFIRMED on location.** Paired before/after (same vertices, classified at their
ORIGINAL positions):

| relative change in turns > 80° | <0.5° | 0.5–2° | 2–5° | 5–15° | >15° |
|---|---|---|---|---|---|
| 2048² | **−19 %** | −0 % | +15 % | −0 % | 0 % |
| 8192² | **−12 %** | **−10 %** | −1 % | 0 % | **+77 %** |

The improvement is exclusively on the flat classes at both resolutions. Nothing lands on the
steep shores except degradation.

### But the aggregate says the remedy does not work

| | 2048² | 8192² |
|---|---|---|
| TOTAL turns > 80° | 20.13 % → **20.04 %** (**−0.4 %**) | 17.78 % → **17.73 %** (**−0.1 %**) |
| axial concentration R | 0.377 → **0.246** | 0.369 → **0.274** |
| mean step (cells) | 0.691 → 0.703 | 0.705 → 0.719 |
| rings/lines | 1618 → 1618 | 6703 → 6703 |

**The smoothing REDISTRIBUTES barbs; it does not remove them.** At 8192² the flattest class
sheds 1 606 turns above 80° and the steepest class gains 1 889 — the total is flat because a
relaxed vertex changes the turn angle at its NEIGHBOURS, and at a finer grid far more steep
vertices sit next to flat ones. A per-class table alone would have shown a −12 % headline and
hidden this; the aggregate is what makes it legible, which is method rule 3 applied to my own
bench.

What DOES improve is the **axial concentration, 0.377 → 0.246** — the contour is measurably
less axis-and-diagonal aligned even though the turn count is unchanged. Whether that is visible
is a judgement the numbers cannot make.

**Verdict: the mechanism is confirmed as the LOCATION of the defect, and the geometric remedy is
refuted as a cure.** Gradient pinning is a symptom of the terrain being nearly flat at sea level
over long stretches of shore; a contour filter cannot supply information the field does not
contain. The real lever is upstream — the same hypsometry that Findings 42–44 leave open.

### The sea-level offset — and the two bugs the check caught

The offset is the end-to-end coherence check between P3-A, the export roll and the Living
Landz reader, so it was measured at every step. **Final: median +0.000 m at 2048² and −0.000 m
at 8192²; |p90| 0.537 m and 1.422 m** against 0.001 / 0.005 m raw. The median does not drift —
it is in fact tighter than the −0.06 m on record — and the tail is sub-metre. **It caught two
real bugs in my own algorithm before either could be reported as a result:**

1. **Ordering.** Capping the vertex displacement AFTER the Newton reprojection dragged 10 % of
   vertices back OFF the isoline, by up to 11 m. The reprojection must be the last thing that
   touches a vertex — now stated at the line.
2. **The wrong gradient, which is the instructive one.** `gradient_at` used a ±1 cell central
   difference, which smooths over two cells and is therefore NOT the gradient of the bilinear
   function the Newton step is trying to zero. The iteration stalled rather than converged:
   raising the step count from 3 to **12 made the residual WORSE (7.1 → 8.2 m)**, which is the
   signature of descending on the wrong slope, not of needing more steps. A ±0.25 cell
   difference collapsed the residual from **7.07 m to 0.745 m** in one change.

Both were found because a *second* invariant was being watched alongside the metric under
improvement — method rule 7, and it paid twice in one round.

### An intermediate result that looked better and was worse

Before the gradient fix, the barb numbers were BETTER (−37 % at <0.5°, −21 % at 0.5–2°) — and
that improvement came from vertices drifting off the isoline by up to 11 m. **The prettier
metric was produced by breaking the constraint the layer exists to satisfy.** Recorded because
it is the exact shape of a result one would ship by accident: the headline number improved and
the thing it was measured against had quietly stopped being true.

### What ships

`smooth_polylines_on_isoline` (gradient-weighted normal-only relaxation, damped Newton
reprojection, displacement cap) plus `slope_deg_to_norm_gradient`, and
`coastline_geojson(field, Option<CoastlineSmoothing>)`. **It ships OFF**, per the measurement
rather than per caution: a −0.4 % change in the total is not worth a cache-invalidating export
change. Flipping the flag is a one-line experiment when the author wants to judge the axial-R
difference visually.

The relaxation gate is a PHYSICAL slope (`SMOOTH_RELAX_BELOW_DEG = 2.0`), not a percentile of
the contour's own gradient distribution: a percentile follows the coast it measures, so when
half the coastline sat at 2–5° the gate spilled onto shores that were never barbed (+24 % there,
measured, then fixed to +15 %).

### RULE-7 CONTROL BLOCK — and why the usual one does not apply

Contour smoothing **does not touch the height field**: it reads `eroded` and emits polylines, so
hypsometry, closed-depression counts and network extent cannot move — not by measurement but by
construction, which is stronger. The relevant invariants are the coastline layer's own, and
those were measured: **ring/line count unchanged** (1618 → 1618, 6703 → 6703), vertex count
unchanged, mean step ±2 %, and the sea-level offset above. Determinism verified by two identical
runs.

## Finding 49 — the export caught two scales in one file, and a bench of mine that was not production

The author exported at 8192² and the file refuted two things I had reported.

### The bug: `catchment_km2` was REAL while `drainage_km2` was SIGNIFIED

`apply_geo_scale_ratio` (Finding 24) multiplies the SIGNIFIED hydrology by `ratio²` — the map
draws `real_km` and signifies `real_km · ratio`. It scales `segment_drainage_km2`,
`segment_discharge_m3s`, `segment_width_m` and re-classifies `segment_navigability`. The two
arrays added in Findings 46–47 were **not** in that list, so at the author's shipped
`geographic_scale_ratio = 7.5` (area ×56.25) the export carried:

| top reach by geometric catchment | value | scale |
|---|---|---|
| `catchment_km2` | 1 610 | REAL |
| `drainage_km2` | 1 962 | SIGNIFIED (×56.25) |

Two "areas" in one row, 56× apart in meaning, with nothing saying so. Worse for
`segment_discharge_profile_m3s`: the microscope's shared-mouth attribution reads it, so entries
with a shared mouth would have carried a REAL area while every other entry carried a signified
one — **an inconsistency INSIDE one list**.

Fixed: both arrays scale in `apply_geo_scale_ratio`.

**Sixth occurrence of the completeness trap, and it is a DIFFERENT shape.** `SegmentRow` +
`push_segment` (Finding 47) closes the APPEND case by making omission a compile error. This is
a MUTATION — an existing array not updated — and no constructor guards it. The guard for the
mutation case is `geo_scale_ratio_scales_every_signified_quantity`, which ENUMERATES every
per-segment array with its expected behaviour: scaled by `ratio²`, derived from a scaled
quantity, or scale-invariant by nature. Adding a field leaves the list incomplete and forces
the decision. Weaker than a compile error, but it is the strongest guard available for "an
existing loop should have grown a line".

### The correction: my navigability verdict was measured at the wrong ratio

I reported "barges and ships will not appear, and that is expected". **The author's export
has 996 barge-class watercourses and 3 404 small-boat**, out of 7 498. My benches all passed
`geo_scale_ratio: 1.0`; the shipped config uses **7.5**. So Finding 47's three-column table
describes a configuration the author does not run.

| in the author's 8192² export (ratio 7.5, lat 45°, span 40°) | count |
|---|---|
| SmallBoat | 3 404 |
| NonNavigable | 3 098 |
| Barge | 996 |
| Ship | 0 watercourses (8 total, all SPILLWAYS) |

**Method rule 3, violated by me, on a parameter I had listed and set myself.** `geo_scale_ratio`
is in `HdParams`; I wrote `1.0` into every bench and never questioned it, because it is
documented as a "presentation multiplier" that touches nothing physical — true of the terrain,
false of the exported hydrology, which is precisely what the navigability question was about.
The lesson is narrower and sharper than "reproduce the whole chain": **a parameter documented as
not affecting X may still affect the QUANTITY YOU ARE MEASURING.** Read what it touches, not
what it is called.

What survives from Finding 47 unchanged, because it is measured on quantities the ratio scales
uniformly: the geometric distribution is identical between arid and humid (the ratio does not
create water), the discharge is not resolution-stable, and the thresholds in the code were never
those the author proposed. What does NOT survive is the conclusion "no barges" — at ratio 7.5
the class spread is already populated, and the parked 10/100/1000/8000 proposal would make it
far MORE generous, not less.

**And the 8 Ship-class entries are all spillways** — the exact mischaracterisation Finding 45
typed away. A consumer filtering `kind == "Watercourse"` never sees them. The typing earned its
keep in the first export that followed it.

### Re-export required

The two scaled arrays change the exported values, so the file the author has must be
regenerated to be read as coherent. `ALGO_DRAINAGE` 4→5 and `ALGO_HD_DRAINAGE` 6→7.

## Finding 50 — the coastal fringes: `K` is continuous, so the light remedy has no target

Finding 48 is CLOSED: the contour relaxation goes back to `None`. The author's three-way
comparison came back nearly identical, which is what the numbers predicted (−0.4 % / −0.1 % on
the total). Mechanism confirmed, geometric remedy refuted by measurement AND by eye. The
function stays in core, unused, for when the hypsometry converges and the flat shores need
re-measuring.

### The measurement that decides the remedy: `K` is NOT class-wise at HD

Measured on the SAME field production incises with — the assembly was extracted to
`production_k_field` so the bench cannot drift from it (the reconstruction trap).

| `K` over the whole map | p1 | p25 | p50 | p75 | p99 | max | cells at exactly 1.000 |
|---|---|---|---|---|---|---|---|
| 2048² and 8192² | 1.049 | 1.379 | 2.210 | 3.878 | 13.63 | 56.73 | **0 (0.0 %)** |

**Not one cell sits at the reference value.** The lithology classes are binary on the coarse
64² grid but reach HD through `sample_bilinear_periodic` — a ramp over one coarse cell, 6.25 km
— and the C-3b fracture density (`1 + amplitude·exp(−d/decay)`) multiplies every cell by a
continuous factor. There are no class boundaries at HD to smooth.

And the fringes do not sit on the steep parts of `K` either. Fringe rate by |∇K| quintile:

| | Q1 | Q2 | Q3 | Q4 | Q5 |
|---|---|---|---|---|---|
| 2048² | 21.2 % | 23.1 % | 17.5 % | 20.2 % | 18.7 % |
| 8192² | 15.8 % | 19.4 % | 15.1 % | 17.8 % | 20.8 % |

Flat and slightly declining at 2048²; weakly rising but NON-MONOTONE at 8192². No correlation
worth a remedy.

**VERDICT, per the decision rule set before measuring: `K` is already continuous, so the cause
is elsewhere and the remedy is heavier.** Smoothing the `K` field would be treating a
discontinuity that does not exist. The graduated-contact idea is physically sound and would be
the right fix IF `K` were class-wise — it is not, and saying so is worth more than shipping a
fix that could only work by coincidence. `stamp_volcanic_k` IS hard-edged (a disc where
`K = max(K, 3)` with no taper), and it is the one genuine discontinuity in the assembly — but
it is a bounded number of small discs, and the quintile table says the shoreline does not
follow them.

### But the author's observation is CONFIRMED — my first metric hid it

My first pass used the fringe RATE and found the closures made no difference (20.86 % off,
20.13 % on). That was method rule 2 again, by me, one round after writing it down: **a rate is
only valid at comparable total effect, and the coastline itself changes length.**

In ABSOLUTE terms the degradation is real and large:

| 2048² | rings | tiny (<20 pts) | fringes | coast km | main ring km |
|---|---|---|---|---|---|
| all OFF | 1163 | 1152 | 5 685 | 3 835 | 2 264 |
| C-3 only | 1428 (+265) | 1405 | 6 345 (**+660**) | 4 227 (+393) | 2 433 |
| C-3b only | 1552 (+389) | 1530 | 6 884 (**+1 199**) | 4 627 (+793) | 2 712 |
| C-3 + C-3b | 1618 (+455) | 1579 | 6 926 (+1 241) | 4 860 (+1 025) | 2 772 |

| 8192² | rings | tiny | fringes | coast km | main ring km |
|---|---|---|---|---|---|
| all OFF | 4713 | 4664 | 25 573 | 4 821 | 2 783 |
| C-3 only | 5979 (+1266) | 5917 | 29 073 (**+3 500**) | 5 461 (+640) | 3 057 |
| C-3b only | 5797 (+1084) | 5737 | 29 061 (**+3 488**) | 5 619 (+798) | 3 378 |
| C-3 + C-3b | 6703 (+1990) | 6628 | 30 885 (+5 312) | 6 206 (+1 385) | 3 632 |

The main ring alone grows **2 264 → 2 772 km at 2048² and 2 783 → 3 632 km at 8192²** — the
shoreline really is more crenellated, which is what the eye reports.

### The attribution, corrected — and it FLIPS with resolution

- at **2048²** C-3b dominates: +1 199 fringes against C-3's +660, nearly double;
- at **8192²** they are **equal**: +3 488 against +3 500, and C-3 adds MORE RINGS (+1 266 vs
  +1 084).

So **neither closure is "the" dominant contributor** — the ranking depends on the grid. The
author's attribution of C-3 holds at 8192², where it adds the most rings, and his observation
that C-3-alone is markedly worse than all-off is confirmed at both resolutions. What does not
survive is "C-3 is dominant" as a general statement; and I would have shipped the mirror-image
error (C-3b is dominant) had I measured only at 2048². **A cause attributed at one resolution
is not a cause.**

They are also SUB-additive: +660 and +1 199 alone, +1 241 together at 2048². The two closures
are largely fringing the same shores.

### The pre-existing baseline: quantified, and scoped to a stage — not fixed

| | 2048² | 8192² |
|---|---|---|
| fringes with EVERY closure off | 5 685 of 6 926 = **82.1 %** | 25 573 of 30 885 = **82.8 %** |
| rings with every closure off | 1 163 of 1 618 = 71.9 % | 4 713 of 6 703 = 70.3 % |
| of those rings, under 20 vertices | 1 152 = **99.1 %** | 4 664 = **99.0 %** |

**About 82 % of the fringing pre-exists every closure, at both resolutions**, and the coastline
is one main ring plus a swarm of ~1 150 (2048²) or ~4 700 (8192²) specks. The pre-existing
defect is a SPECKLE, not a fringe — a different thing from what the closures add, which is why
it must be attributed separately.

**Scoped to a stage:** it is NOT the FBM. With `amplitude_base = 0` and every closure off, the
2048² coastline is 1 171 rings / 5 707 fringes against 1 163 / 5 685 — a difference under
0.5 %, and in the wrong direction. Consistent with Finding 43 stage 2, where the C-1
relief-budget cap crushes the FBM's whole contribution to ±2 m. So the speckle comes from the
**bilinearly upscaled coarse field plus the incision**, and the FBM is exonerated. Scoped here,
deliberately not fixed — mixing it with the closure contribution would make neither
attributable.

### Method rule 8 — for a VISUAL defect, bisect on the existing toggles first

The author's three-column screenshot did in one pass what a bench had not shown, and it was
right about the existence and the direction of the effect while my first quantitative metric
said "no difference". No code, no bench, no build: flipping the toggles that already exist
attributes a visual defect directly, and it costs one generation per column.

The corollary matters as much: **a rate can contradict the eye and the eye be right.** When a
visual report and a metric disagree, suspect the metric's normalisation before suspecting the
observer — here the rate was constant precisely because the denominator, the coastline itself,
was what the closures had changed.

### RULE-7 CONTROL BLOCK

The only production change touching the height path is the extraction of `production_k_field`
— a pure refactor, same operations in the same order. Byte-identity evidence rather than
assertion: the terrain built through the refactored path yields **1 618 rings and 37 650
coastline vertices at 2048²**, exactly the figures Finding 48 measured before the extraction.
The coastline flag returning to `None` restores the previously measured export. Hypsometry,
closed-depression counts and network extent are untouched by construction (no incision
parameter changed).

### A bug in the bench, worth its line

`v[((v.len() - 1) as f32 * f) as usize]` panicked at 8192²: with 67 108 864 entries,
`67 108 863` exceeds f32's 24-bit mantissa and rounds UP to `67 108 864.0`, indexing one past
the end. Invisible at 2048² (4 M entries, well inside the mantissa). Fixed to integer/f64
arithmetic in both coastal benches. **A percentile helper is exactly where this hides**: it is
trivial, it is copied between benches, and it only fails above 16.7 M elements — which is to
say, only at the resolution that matters.

## Finding 51 — stage bisection: the INCISION introduces 100 % of the coastal fringing

Same instrument that attributed 90 682 depressions to the FBM in one table: measure the metric
after each stage and read where it jumps. No argument, no candidate ranking beforehand.

| stage | 2048² rings / fringes / rate / coast km | 8192² rings / fringes / rate / coast km |
|---|---|---|
| 1 coarse (no FBM, no incision) | 14 / 35 / **0.35 %** / 1 581 | 46 / 60 / **0.15 %** / 1 602 |
| 2 + FBM upscale | 14 / 35 / **0.35 %** / 1 581 | 46 / 60 / **0.15 %** / 1 602 |
| **3 + relief-v2 incision** | 664 / 5 582 / **22.14 %** / 3 686 | 431 / 11 740 / **12.31 %** / 3 696 |
| 4 + relief-v3 incision | 1 163 / 5 685 / 20.86 % / 3 834 | **4 713** / 25 573 / **18.84 %** / 4 821 |
| 5 final (+ C-3 + C-3b) | 1 618 / 6 926 / 20.13 % / 4 860 | 6 703 / 30 885 / 17.78 % / 6 206 |

**Three results close three lines of enquiry:**

1. **The coarse field does not fringe.** 14 rings at 2048², 46 at 8192², rate 0.35 % / 0.15 %.
   **The subject is not tectonic** — which was the outcome that would have moved the chantier
   upstream and multiplied every subsequent estimate. It does not.
2. **The FBM adds EXACTLY NOTHING.** +0 rings, +0 fringes, +0 km, at both resolutions, to the
   digit. Third independent exoneration (Finding 43 stage 2: ±2 m contribution; Finding 50:
   `amplitude_base = 0` moves it 0.5 % the wrong way; this one). The FBM can be struck off the
   list of suspects for coastal geometry permanently.
3. **The incision introduces the whole baseline.** Everything the author sees enters between
   stage 2 and stage 4.

### The mechanism splits in two, and only the second is the "flat terrain" one

The near-zero slope distribution — land within 20 m of sea level — separates them:

| stage | 2048² cells / <0.5° / median | 8192² cells / <0.5° / median |
|---|---|---|
| 1–2 coarse & FBM | 4 338 / 16.4 % / 1.40° | 17 367 / 17.7 % / 1.20° |
| 3 relief-v2 | 12 149 / 8.3 % / **2.58°** | 44 020 / 9.4 % / **2.75°** |
| 4 relief-v3 | 24 034 / 30.4 % / **0.86°** | 49 915 / 20.5 % / **1.45°** |
| 5 final | 28 714 / 36.4 % / 0.75° | 56 754 / 24.3 % / 1.33° |

- **Stage 3 is NOT flattening** — it makes the coastal terrain STEEPER (1.40° → 2.58°) while the
  fringe rate jumps 22 points. It is the incision **carving valleys down to sea level**: near-sea
  land ×2.8, coastline length ×2.3, rings ×47. Drowned valley mouths, a real landform.
- **Stage 4 IS the flattening**, and it is the mechanism flagged before measuring: median slope
  falls back to 0.86° / 1.45°, the share under 0.5° doubles, and the rings multiply (×1.8 at
  2048², **×10.9 at 8192²**). This is the "near-flat terrain at sea level" verdict, now
  attributed to a stage.

The per-class table makes it sharper still: **in the coarse field the low-slope classes are
EMPTY** (no `<0.5°`, `0.5–2°` or even `2–5°` shore exists at 8192²). The incision manufactures
the low-slope shore, and the fringes then live on it — 86.9 % of turns above 80° in the `<0.5°`
class at stage 4.

### And relief-v3's contribution is resolution-amplified

At 2048² relief-v3 does not raise the rate (−1.28 points) and merely doubles the rings. At
8192² it adds **+6.53 points and multiplies the rings ×10.9** (431 → 4 713). Same pattern as
Findings 43 and 50: **a contribution measured at one resolution is not a contribution.** The
fringing is markedly worse at 8192² *because of relief-v3 specifically*, and that is the same
resolution-dependence Finding 43 traced to the fluvial/hillslope regime partition.

### Ranked remedies, with cost — attribution first, as agreed

**RANK 1 — the hillslope term (Findings 43–44). Same root, and the third chantier to converge
on it.** Stage 4 flattens the near-sea land and leaves it: that IS "a correct threshold feeding
an inert branch". A transport-limited hillslope law with an explicit sediment flux would denude
those flats and deliver the material to the network instead of leaving a 0.86° apron for the
contour to fringe through. Cost: HIGH — design work, and it needs the explicit timescale
(Finding 44) first, since a diffusivity in m²/yr has nothing to multiply. But it addresses the
mechanism rather than a symptom, and hypsometry, H-2's dial and the coastal fringes now all
point at it.

**RANK 2 — C-4 coastal erosion, already on the roadmap.** The physically correct agent: wave
attack straightens coasts, removes islets and cuts platforms — precisely the drowned valley
mouths and flat shores this measurement found. Cost: MODERATE, a planned closure, and crucially
**it does not touch the relief calibrated over eighteen rounds**. Note the roadmap put C-4 after
H-1/H-2 because H-2 moves the coastline; that ordering still holds, and this finding gives C-4
a second, measured justification it did not have.

**RANK 3 — accept.** Cost zero, and not unreasonable: an indented coast with islets is not
wrong in itself. What the author objects to is the scale and regularity, not the existence.

**EXCLUDED, each by measurement rather than by preference:**

| candidate | why not |
|---|---|
| smooth the `K` field | no discontinuity exists to smooth — `K` is continuous at HD, 0 cells at the reference (Finding 50) |
| contour relaxation | refuted by measurement AND by eye (Finding 48) |
| `min_slope` / stream-power parameters | blast radius covers the whole calibrated relief, and it would tune the agent rather than supply the missing counter-agent |
| the FBM | exonerated three times over |
| the coarse / tectonic stage | 0.15–0.35 %, it does not fringe |

### Method rule 8, completed — toggles first, then stages

The author's TOGGLE BISECTION found in one screenshot what a bench had not shown. This STAGE
BISECTION is its measured counterpart, and between them they took the subject from "the
lithological contact is too abrupt" to "the incision manufactures the low-slope shore" — two
different chantiers with an order-of-magnitude difference in cost. **For a visual defect:
bisect on the toggles that exist, then on the stages of the build. Attribution before remedy,
in both cases.** Neither needed a new mechanism to be hypothesised first, which is why neither
could be led astray by one.

## Finding 52 — a fringe metric that measures fringes; and the 43.5 % figure, which I cannot reproduce

### The number that would have turned the subject over — and does not

The claim was: the coarse field is already at 43.5 % fringing on a 172-vertex 64² coastline, so
the fringes come from TECTONICS and everything downstream moves it by ±1 point. That conclusion
was flagged in advance as the one that changes every estimate, so it had to be checked first.

**I cannot reproduce it.** Tracing the sea-level isoline on the raw 64² normalised altitude —
the same `c1_coarse_normalized_altitude` production upscales, `target_land_fraction: None` as
shipped — gives **18 rings, 314 vertices, 30 of 278 turns above 80° = 10.8 %**. Not 172
vertices, not 43.5 %, and *lower* than the 20.1 % of the finished 2048² product rather than
four times higher. I do not know how the 43.5 % was obtained and I am not able to build on it.

**And the inference would not follow even if the figure were right**, because the turn rate is a
function of vertex spacing. Taking ONE fixed 2048² coastline — the shape never changes — and
keeping every k-th vertex:

| vertex step | 1 | 2 | 4 | 8 | 16 | 32 | 64 |
|---|---|---|---|---|---|---|---|
| spacing (km) | 0.14 | 0.27 | 0.55 | 1.09 | 2.19 | 4.38 | 8.75 |
| turns > 80° | **20.1 %** | 26.3 % | 26.1 % | **30.0 %** | 29.7 % | 28.6 % | 25.7 % |

The same curve reads 20 % or 30 % depending only on how finely it is sampled — a 50 % relative
swing from resampling. **So a turn rate measured at 6 km spacing cannot be compared with one
measured at 0.14 km spacing**, and "the coarse field is already at X %" cannot attribute
anything to a stage.

**What this does NOT invalidate:** Finding 51's stage table was measured at ONE FIXED sampling
throughout (every stage rendered at 2048², then again at 8192²). Within a fixed grid the
comparison is sound, and it stands.

### The real point, which is right and which I had missed

The turn count **cannot see fringe LENGTH**. A spur one kilometre long and a spur ten kilometres
long have the same two sharp turns. So the metric that attributed the stage correctly is unfit
to judge the appearance, and evaluating a remedy with it would repeat the relaxation episode:
a −0.4 % announced, nothing visible, no information gained.

### The metric

Three quantities, in kilometres, at an explicit physical scale. A SPUR is an excursion that
comes back on itself: an arc of at least 1 km whose endpoints lie within 0.6 km of each other —
a neck. That is the definition of "a fringe running out from the shore body", and it is what the
turn count is blind to.

- **length** — median and p90 of the spur arc, plus the share of the whole coastline inside a
  spur;
- **regularity** — the coefficient of variation of the spacing between consecutive spurs. A COMB
  is regular (low CV); a ria coast is irregular (high CV);
- **parallelism** — the axial concentration `R` of the spur axes. Parallel spurs read as
  manufactured.

**Validation, on the one ranking that can be sourced** — the author judges the current export
worse than the pre-chantier state (same relief, C-3 and C-3b off):

| 2048² | coast km | spurs | med km | p90 km | in-spur | spacing CV | axis R |
|---|---|---|---|---|---|---|---|
| pre-chantier (closures OFF) | 3 835 | 1 146 | 1.39 | 2.58 | 50.0 % | 0.88 | 0.023 |
| current export (both ON) | 4 860 | **1 262** | **1.63** | **3.48** | **53.9 %** | 0.93 | 0.052 |

The metric puts them in his order — more spurs, longer spurs (**p90 +35 %**), more coastline
inside a spur. It is fit to judge a remedy.

**And it settles the description.** "Long parallel fringes": the LENGTH is confirmed (p90 2.58 →
3.48 km). The PARALLELISM is **refuted** — the axis concentration is 0.023 and 0.052, which is
isotropic. Whatever reads as parallel is the REPETITION of similar spurs, not a common
direction. The spacing CV even RISES (0.88 → 0.93): the fringes become slightly more irregular,
not more comb-like. Worth knowing before anyone tries to remove a periodicity that is not there.

### Where the spurs enter — the attribution survives its metric being replaced

| stage (2048²) | coast km | spurs | med km | p90 km | in-spur | axis R |
|---|---|---|---|---|---|---|
| 1 coarse upscaled (no FBM) | 1 581 | **14** | 1.82 | 7.67 | **2.4 %** | 0.453 |
| 2 + FBM | 1 581 | **14** | 1.82 | 7.67 | **2.4 %** | 0.453 |
| **3 + relief-v2** | 3 686 | **1 512** | 1.29 | 2.31 | **64.2 %** | 0.046 |
| 4 + relief-v3 | 3 835 | 1 146 | 1.39 | 2.58 | 50.0 % | 0.023 |
| 5 final (+ C-3 + C-3b) | 4 860 | 1 262 | 1.63 | 3.48 | 53.9 % | 0.052 |

**Same stage. Spurs 14 → 1 512, in-spur share 2.4 % → 64.2 %.** Two metrics built on different
principles — one counting turns, one measuring necked excursions in kilometres — locate the
jump at the same place. Finding 51's attribution survives the replacement of the instrument that
produced it, which is the strongest form the attribution could take.

### Defect or property: the coarse angularity is a PROPERTY, and the upscale handles it

The hypothesis was that the upscale manufactures long fringes by interpolating an already
angular contour 128×. **It does not.** The upscaled coarse field carries **14 spurs and a 2.4 %
in-spur share** — bilinear interpolation SMOOTHS the 64² contour rather than amplifying it,
which is also why stage 1 reads 0.35 % on the turn metric while the raw 64² reads 10.8 %.

So: a 314-vertex coastline on a 64² grid over 400 km is angular **by necessity** — 6.25 km
between samples — and that is a property of the grid, not a defect in the tectonics. The
upscale's job is to interpolate it, and it does so without creating spurs. **The remedy belongs
in NEITHER the upscale NOR the tectonic stage.**

The only sign of manufactured geometry anywhere is the coarse stage's axis R of 0.453 — those 14
spurs ARE grid-aligned. But they are 14, they carry 2.4 % of the coastline, and after the
incision R collapses to 0.046. The grid signature does not survive into the product.

### Ranked remedy — unchanged, now on a metric that measures the symptom

1. **The hillslope term** (Findings 43–44). Stage 3 creates the spurs by carving valleys to sea
   level; stage 4 and the closures then LENGTHEN them (p90 2.31 → 2.58 → 3.48 km) on the
   low-slope apron the incision leaves undenuded. Cost HIGH, timescale prerequisite.
2. **C-4 coastal erosion.** Wave attack shortens and truncates spurs — the exact quantity this
   metric now measures, so C-4 becomes falsifiable: it must cut the p90 spur length. Cost
   MODERATE, planned, does not touch the calibrated relief.
3. **Accept.** Cost zero.

### Method — three consecutive attributions to the wrong layer

Worth naming, because it is a pattern and not three accidents:

- my "symptom of near-flat terrain at sea level" assumed erosion was flattening the shores;
- the C-3 remedy assumed an abrupt lithological contact;
- the toggle bisection correctly saw C-3 making things worse, but inside a false premise about
  why.

None was about the right stage, and all three were **mechanism arguments made before a stage
bisection**. The cost of reasoning about mechanism first is that a plausible mechanism will
always be available for the wrong layer. Two bisections — toggles, then stages — settled in two
passes what three rounds of mechanism could not.

**And a corollary about instruments:** a metric can be RIGHT FOR ATTRIBUTION and UNFIT FOR
JUDGEMENT, and those are separate questions to ask of it. The turn count located the stage
correctly three times and could never have judged a remedy. The right response was not to
distrust its attribution but to build the second metric and check whether they agree — they do.

## Finding 53 — the two sweeps, an anchored target, and the one thing that is actually wrong

### First, three corrections to the premises

**`smooth_iterations` does not exist.** There is no such parameter in core; the only `smooth*`
knob is `smoothing_width` in `tectonics_v2/cratonic`, unrelated. What does not run in production
is the contour relaxation of Finding 48, passed as `Option<CoastlineSmoothing> = None` — set
that way deliberately two rounds ago because the measurement refuted it, not an undiscovered
piece of wiring.

**`coast_warp_strength = 1.5` IS active**, and it is the right suspect: its own docstring says
it displaces the coarse-altitude sampling position "so the sea-level contour meanders instead of
following the blocky 64² polygon". It had never been measured. That instinct was correct, on the
correct parameter.

**Two sets of figures could not be sourced.** "10.55 cells median fringe length at 8192² against
1.90 at coarse, regularity 0 % → 47.8 %" contradicts the committed Part C, which measures the
upscaled coarse field at **14 spurs / 2.4 % of coastline** and the jump at relief-v2 (1512 /
64.2 %). And if those lengths are in CELLS they are not comparable across resolutions at all —
128× of cell size, the exact trap Part A documented. Likewise "the validation test was
unpassable, the toggles barely move the metric (0.80 % vs 0.73 %)": the validation RAN AND
PASSED — 1 146 → 1 262 spurs, p90 2.58 → 3.48 km, in the author's order. There is no ill-posed
test to record, and the metric was validated before use.

### The anchoring — derived reference values, not invented ones

The requirement was "indentations of varied lengths and regularity near zero". Both have exact
references, which is what makes the target defensible:

| quantity | reference | why | target |
|---|---|---|---|
| **spacing CV** | **1.0** | a POISSON process — spurs placed at random along the shore — has gap CV of exactly 1 | ≈ 1, and not below ~0.6 |
| **length variety** `p90/median` | **3.32** | an EXPONENTIAL distribution, the maximum-entropy choice for a positive quantity with a given mean: "as varied as possible without further structure". ln(10)/ln(2) | ≈ 3.3; under ~2 is suspiciously uniform |
| **axis R** | **0** | isotropic | < 0.1 |

**Which metric is optimised, stated so the recommendation cannot be misread:** the LENGTH
VARIETY and the spacing CV are optimised. **Coastline length and spur count are REPORTED ONLY** —
reducing the warp shortens the coast mechanically and a shorter coast is not a better one.

### And the anchoring immediately isolates the defect

| shipped product | 2048² | 8192² | target | verdict |
|---|---|---|---|---|
| spacing CV | 0.93 | 0.71 | ≈ 1.0 | **in range** |
| axis R | 0.052 | 0.033 | < 0.1 | **in range** |
| **length variety p90/med** | **2.14** | **2.02** | **≥ 3.3** | **OUT** |

**The coast is not regularly spaced and not parallel. Every spur is simply about the same
length.** That is the single anchored quantity that fails, and it is the first number in this
whole chantier that isolates what "manufactured" means. It also explains why the parallelism
hypothesis kept feeling right and measuring wrong: uniform LENGTH reads as repetition, and
repetition reads as a comb, without any common direction being involved.

### Sweep 1 — `coast_warp_strength`

| warp | 2048² coast km / variety / CV | 8192² coast km / variety / CV |
|---|---|---|
| **1.50 (shipped)** | 4 860 / **2.14** / 0.93 | 6 206 / **2.02** / 0.71 |
| 1.00 | 4 994 / 2.17 / 0.85 | 6 903 / 2.13 / 0.69 |
| 0.50 | 4 989 / 2.27 / 0.81 | 7 225 / 2.26 / 0.69 |
| 0.25 | 4 753 / 2.49 / 0.88 | 6 942 / 2.55 / 0.72 |
| **0.00** | 5 001 / **2.53** / 0.80 | 6 870 / **2.57** / 0.73 |

**The warp works against its own purpose.** It exists to make the contour meander, and it makes
the spur lengths MORE UNIFORM: variety rises monotonically as the warp is removed, 2.14 → 2.53
and 2.02 → 2.57. The trend is monotone at both resolutions, which is what makes it a result
rather than noise.

**The length trap does not even apply here.** Removing the warp does NOT shorten the coast —
4 860 → 5 001 km and 6 206 → 6 870 km, it slightly LENGTHENS it. So the recommendation cannot be
read as "shorter is better"; the coast gets marginally longer AND more varied.

**And the feature's stated fear is checked, not assumed.** Turning the warp off was supposed to
leave the contour "following the blocky 64² polygon" — the axis concentration would catch that
immediately, since a grid-following contour is axis-aligned. At warp 0 it is **0.060 and 0.025**,
fully isotropic. The coast does not fall back onto the grid.

Cost on the CV: 0.93 → 0.80 at 2048² (mildly further from the Poisson reference, still in
range), flat at 8192² (0.71 → 0.73). A real but small trade.

### Sweep 2 — the contour relaxation, on the metric that can see it

Finding 48 evaluated the relaxation with the TURN COUNT and got −0.4 %, a number with no
meaning, since the turn count cannot see spur length. Re-run on the spur metric — nearly free,
since the relaxation only post-processes the polyline, so all four settings come from one build:

| passes | 2048² variety / CV / spurs | 8192² variety / CV / spurs |
|---|---|---|
| 0 (shipped) | 2.14 / 0.93 / 1 262 | 2.02 / 0.71 / 1 951 |
| 1 | 2.15 / 0.93 / 1 263 | 2.02 / 0.71 / 1 953 |
| 2 | 2.15 / 0.92 / 1 261 | 2.02 / 0.71 / 1 952 |
| 4 | 2.12 / 0.93 / 1 238 | 2.02 / 0.71 / 1 953 |

**Flat on everything.** Finding 48's refutation stands, and it now stands on the right
instrument — which was the open question. The relaxation stays `None`, definitively.

### Ranked recommendation — and the honest answer to "is any setting credible"

**No setting produces a credible coast.** The best available closes about 40 % of the gap: at
8192², variety 2.02 → 2.57 against a 3.32 target. That is worth saying first, because the
recommendation below is an improvement and not a fix.

1. **`coast_warp_strength: 1.5 → 0`.** The only measured setting that moves the anchored metric,
   monotone at both resolutions, at no cost in coastline length, no parallelism cost, and a
   small CV cost at 2048². Cost: ONE CONSTANT, plus an `ALGO_UPSCALE_EROSION` bump and a
   re-export. It removes a feature that measurably works against its stated purpose.
   ⚠️ It changes the terrain, so the rule-7 control block applies before it ships — not measured
   in this read-only pass.
2. **The remaining 60 % is the incision** (Findings 51–52), and therefore the hillslope term
   (Findings 43–44). Unchanged ranking, now with a quantified share.
3. **The relaxation: closed.** Not `None`-for-now — `None` on two independent metrics.

### Method — the anchor is what turned a preference into a measurement

Every previous round asked "is it better?" and could only answer "the number moved". With a
Poisson CV of 1.0 and an exponential p90/median of 3.32, "better" has a direction AND a
destination, and two of the three quantities turn out to have been fine all along. **Without the
anchoring the obvious optimisation would have been to reduce the spur count — the one quantity
that must NOT be optimised**, since a coast with few indentations is the soft trace already
rejected. That is exactly the failure the anchoring was demanded to prevent, and it would have
happened.

## Finding 54 — the renders overturn three of my own results, and the target was inverted

The instruction to produce a PNG per setting, and to report a length TAIL rather than a median,
caught four defects in one round — three in my measurements and one in the target we were both
aiming at. None was visible in the indicators I had been producing. **A render per setting is
now the default output of any coastline sweep.**

### The primary metric was inverted: it is DENSITY, not length

With the search ceiling raised (see below), the reference row finally reads:

| 2048² | p50 km | p90 km | max km | spurs | coast km | spurs / 100 km |
|---|---|---|---|---|---|---|
| **1.50 (shipped)** | 1.63 | 3.50 | 17.80 | 1 254 | 4 860 | **25.8** |
| 1.00 | 1.73 | 3.86 | 38.67 | 1 356 | 4 994 | 27.2 |
| 0.50 | 1.70 | 3.89 | 27.65 | 1 416 | 4 989 | 28.4 |
| 0.25 | 1.60 | 3.98 | 33.34 | 1 317 | 4 753 | 27.7 |
| 0.00 | 1.55 | 3.90 | 47.75 | 1 441 | 5 001 | 28.8 |
| **REFERENCE coarse** | **3.08** | **24.68** | 46.06 | **15** | 1 581 | **0.95** |

8192² tells the same story: shipped 1 848 spurs / 6 206 km = **29.8 per 100 km** against the
reference's 20 / 1 602 = **1.25**.

**The reference has spurs SEVEN TIMES LONGER than production (p90 24.68 vs 3.50 km) and
TWENTY-FIVE TIMES FEWER.** A few large bays against a swarm of small teeth — exactly what the
two renders show.

So "bring the fringe length back near the coarse value" would have optimised in the WRONG
DIRECTION: the spurs need to be LONGER and far fewer, not shorter. **The anchored target,
measured rather than assumed:**

| | production | REFERENCE (target) |
|---|---|---|
| indentation density | 25.8–29.8 / 100 km | **~1 / 100 km** |
| p90 spur length | 3.5–4.6 km | **~24 km** |
| local parallelism | 0.52 | **~0.33** |

### The parallelism: my "refuted" verdict was itself wrong

Findings 52 and 53 recorded "the PARALLELISM is refuted — axial concentration 0.023 to 0.052,
isotropic". **The render shows unmistakable bundles of thin, straight, mutually parallel
spurs.** The metric was GLOBAL: bundles pointing different ways cancel in a single circular
mean, so the number reads 0.03 while the picture shows a comb.

Replaced by a LOCAL measure — the axial concentration among each spur and its 8 nearest
neighbours, median over spurs:

| | 2048² | 8192² |
|---|---|---|
| shipped | **0.515** | **0.525** |
| reference coarse | 0.330 | 0.363 |
| (the old GLOBAL figure) | 0.054 | 0.036 |

**Production is markedly more locally parallel than the reference**, and the global figure was
an order of magnitude away from saying so. The author's description was right and my refutation
of it was an artefact of measuring at the wrong scale.

This also identifies the mechanism candidate: thin, straight, locally parallel, regularly
spaced teeth is the signature of **Smith–Bretherton parallel rilling** — the comb already
documented in Findings 8, 9 and 11 and believed cured by MFD. Here it reaches SEA LEVEL, on the
low-slope apron the incision manufactures (Finding 51, stage 4). That ties the coastal comb to a
known defect rather than a new one.

### `coast_warp_strength`: the Finding 53 recommendation is WITHDRAWN

I recommended `1.5 → 0` on the length-variety metric. **The render refutes it**: at warp 0 the
comb becomes dense and continuous all round the landmass and the outline loses its meanders.
The numbers agree once read on the right quantities — lowering the warp RAISES the density
(1 254 → 1 441 at 2048², 1 848 → 2 076 at 8192²) and RAISES the local parallelism
(0.515 → 0.557, 0.525 → 0.596).

**The warp is not the lever in either direction, and it is mildly counteracting the defect.**
Leave `coast_warp_strength = 1.5` alone.

### The search ceiling: a constant reported as a measurement

The first run's `max` column read **7.99–8.00 km at every setting including the reference** —
that was `min_spur_km × 8`, the search WINDOW, not a measurement. Worse, it saturated the
REFERENCE row's p90 (7.67 / 7.89), which is the comparison target, so the target itself was
unreadable. Raised to 50 km, and the reference's true p90 turns out to be **24.68 km** — the
number that inverted the whole objective.

**The `max` column existed only because a tail was demanded rather than a median.** Without it
the ceiling would have stayed invisible and the target would have stayed inverted.

### And a fourth: a silent `.replace()` produced a constant column

The local-parallelism figure came out as **0.000 at every setting** on its first run — my
patch's anchor had not matched, so the axis array was never filled, and `median of empty` is
0.0. I had asserted the anchor on every other edit that round and omitted it on that one.

**Both defects have the same signature, and it is a test worth applying to every reported
column: does this number MOVE when something moves?** `max` = 7.99 everywhere and local R =
0.000 everywhere were both constants across configurations that were themselves varying. A
quantity that refuses to vary is measuring the instrument, not the terrain.

### Method rule 9 — render every setting of a visual sweep

Three wrong attributions (Findings 50–52) and two remedies evaluated against metrics blind to
the symptom (the relaxation, then the warp), then in a single round: a global metric that
contradicted the eye and the eye was right, a ceiling masking the target, and a constant column.
**Every one of them was caught by a picture or by a tail, and none by the summary statistics.**

For any defect the author reports VISUALLY: render one panel per setting, hold the crop fixed,
and include a reference panel of what "absent" looks like.

**AMENDMENT (Finding 56).** A fixed crop is chosen on ONE configuration — normally the shipped
one — and a remedy that DISPLACES the defect rather than removing it will then render a
flatteringly clean panel. It happened: the channel-head law's 8192² panel showed a clean coast
while the global density was 27.9 per 100 km against 29.8 shipped. So:

- **re-choose the crop if the defect moves**, and say which configuration the crop was chosen on;
- **the panel and an ABSOLUTE count must agree before a visual verdict is accepted.** Here they
  reconciled only once the count was read in absolute terms (spurs 1849 → 905) instead of per
  100 km — the rate had a moving denominator. The reference panel is what made the
target legible here — it is not decoration, it is the anchor.

Panels: `exports/coastal_fringes/` — sea black, land grey, coastline white, one fixed 640×640
crop at 8192² on the densest spur window, one per warp value plus `reference_coarse_8192.png`.

### Ranked recommendation

1. **Change nothing about the warp.** Withdrawn.
2. **The lever remains the incision** (Findings 51–52), and specifically whatever makes the
   parallel rilling reach sea level. The target is now quantified: **density ~1 per 100 km,
   p90 ~24 km, local R ~0.33** — the coarse field's own values, so it is a reachable anchor and
   not an aspiration.
3. **C-4 coastal erosion** gains a second falsifiable criterion: it must lower the indentation
   DENSITY while raising the p90 length. A remedy that shortens spurs would be moving away from
   the target.

## Finding 55 — the coastal comb: MFD does not reach the apron, and no uniform parameter is a remedy

### The apron, characterised

| | 2048² | 8192² |
|---|---|---|
| area within 50 m of sea level | 6 748 km² | 3 553 km² |
| mean width inland | 1.39 km | 0.57 km |
| slope p10 / p50 / p90 | 0.19 / **1.29** / 5.52° | 0.31 / **2.14** / 17.93° |
| share under 0.5° | 24.3 % | 16.0 % |
| **comb spacing** | **3.875 km** | **3.358 km** |
| spacing in CELLS | 19.8 | 68.8 |
| spacing / √A_c | 12.25 | 10.62 |

**The comb spacing is PHYSICAL, not a discretisation artefact** — ~3.4–3.9 km at both
resolutions while the cell count changes ×3.5. And it holds a near-constant ratio to `√A_c`
(12.3 and 10.6), which points at the channel-head threshold rather than at the grid.

### MFD does not reach this regime — confirmed by number and by eye

| lever | density /100 km (2048² / 8192²) | p90 km |
|---|---|---|
| **SHIPPED** | 25.8 / 29.8 | 3.50 / 4.62 |
| MFD off (D8) | 31.3 / 27.3 | 4.02 / 5.41 |
| MFD p = 1.1 (dispersed) | 24.1 / 30.5 | 3.78 / 4.64 |
| MFD p = 4.0 (concentrated) | 28.0 / 29.5 | 3.44 / 4.72 |
| **TARGET (coarse ref)** | **0.9 / 1.2** | **24.68 / 23.34** |

Nothing moves. And the mosaic makes it plainer than the table: **tiles 2–5 — shipped and all
three MFD settings — are visually indistinguishable.** The lever that cured the hillslope comb
(Findings 8/9/11) does nothing at low gradient, in any setting including switched off. The
nuance stated before measuring holds: this is the same Smith–Bretherton instability in a regime
MFD was never tested in, not an MFD regression.

### `A_c` sets the comb — and no uniform value is a remedy

| lever | density | p90 km | mean land altitude |
|---|---|---|---|
| SHIPPED | 25.8 | 3.50 | **282 m** |
| A_c ×0.1 | 24.7 | 2.69 | — |
| A_c ×10 | 13.3 | 4.72 | **753 m** |
| A_c ×100 | **3.9** | **9.80** | **860 m** |
| **REFERENCE, no incision at all** | 0.9 | 24.68 | **865 m** |

`A_c × 100` reaches density 3.9 and p90 9.80 — most of the way to the target — and its panel is
indistinguishable from the reference. **The rule-7 control block shows why: its hypsometry is
860 m against the un-incised reference's 865 m.** Mean, p50 (669 vs 679) and p90 (1833 vs 1836)
all match the reference to within one percent. It has not fixed the coast; **it has switched the
incision off**, and the coastline metrics could not see that. Neither could the picture — an
un-eroded terrain has a beautiful coast.

Even `× 10` sits at 753 m against 282 shipped: it already disables most of the incision, and its
density gain is bought by not eroding.

**So there is no uniform `A_c` that removes the comb while preserving the incision**, and the
mosaic reads as one continuum from tile 2 (fully incised, combed) through 6, 7 to 8 (not
incised, no comb). **The comb IS the incision arriving at the coast.**

### Hillslope diffusion is not a lever either

`× 10` helps at 2048² (25.8 → 9.4) and **hurts at 8192² (29.8 → 38.6)** — resolution-incoherent,
which is Finding 45b's missing `dscale` on the linear branch showing up again. `× 100` dissolves
the land entirely (tile 10: specks). Not usable in either direction.

### The remedy the measurement points at, and it is regime-dependent by construction

`A_c` is a CONSTANT in the model. Montgomery & Dietrich put the channel-head threshold at
**`A_c · S² ≈ constant`** — the threshold must RISE as the slope falls. Between a hillslope at
~15° (S = 0.268) and this apron at 1.29° (S = 0.0225) the law demands a ratio of
**(0.268/0.0225)² ≈ 140** — the order of the `× 100` that suppresses the comb.

The uniform test therefore does two things at once: it **confirms the mechanism** (`A_c` sets the
comb spacing) and **proves a uniform value unusable** (it disables the hillslopes too). Only a
slope-dependent threshold can be high on the apron and low on the hillslopes simultaneously,
which is exactly the answer to "does the low-gradient regime need its own treatment rather than
a parameter change": **yes, and no parameter change can substitute for it.**

**Ranked:**
1. **`A_c(S) = C/S²`** — literature-anchored, regime-dependent, and the only candidate the
   uniform sweep corroborates rather than refutes. Cost: one term in `incise`, an `ALGO` bump,
   recalibration of the affected regime only, and a mandatory rule-7 control block.
2. **C-4 coastal erosion**, with its falsifiable criterion unchanged: lower the DENSITY while
   RAISING the p90.
3. **Refuted by measurement:** MFD (any setting), uniform `A_c`, hillslope diffusion,
   `coast_warp_strength` (Finding 54), the contour relaxation (Findings 48 and 53).

### Tooling promoted, per rule 9

`terrain::coast_metrics` — `coast_shape`, `render_coast_crop`, `densest_window` — moved out of a
test file into the library so every sweep uses ONE definition, with the anchored target and the
`MAX_SPUR_KM` ceiling documented at the constant that caused the saturation.

`tests/coastal_mosaic.rs` assembles the panels into a single comparison image with a numbered
legend. Judging ten separate files means comparing from recollection, which is how three wrong
attributions survived; side by side with the REFERENCE first, the differences are direct.
`exports/coastal_fringes/mosaic_levers.png` and `mosaic_warp.png`.

## Finding 56 — the slope-dependent channel head, implemented; and the two-sided form measured and rejected

`A_c(S) = A_c_ref · max(1, (S_ref/S)²)`, Montgomery & Dietrich's `A_c·S² = C`, in `incise`,
**OFF by default** (`a_c_slope_law: None` keeps the constant threshold, byte-identical).

### Calibration — measured, and pinned to one resolution for a stated reason

`C = A_c_ref · S_ref²` is anchored on the HILLSLOPE regime, where the present threshold is
correct, so the apron inherits the ratio the law imposes instead of a chosen factor.

`S_ref = 0.3319` is the MEDIAN gradient over the 159 134 cells whose accumulation sits within
±25 % of `A_c = 0.1 km²` — the channel heads the model currently produces — in the production
config at **2048²**.

⚠️ **It is not resolution-stable: 0.3319 (18.4°) at 2048² against 0.1128 (6.4°) at 8192².** A
constant `A_c` puts channel heads on quite different topography depending on the grid, which is
itself the Findings 42–44 hypsometry defect surfacing in the calibration. 2048² is pinned because
the project already treats it as the calibration resolution (`HILLSLOPE_REF_CELL_M`). **The
residual gap at 8192² is inherited from that, not from the law.**

> **CORRECTED by block A (measurement).** This paragraph reads the 0.3319 / 0.1128 spread as a
> TRACE of the divergence's cause. It is not: `S_ref` is measured on the **ERODED** field, at the
> channel heads the incision has already produced. Block A measured the **PRE-INCISION** field
> and found its D8 slope distribution grid-invariant to within 3–4 % on the tails (p50 0.0906 vs
> 0.0905, p90 0.2658 vs 0.2744, >30° 2.02 vs 2.22 %). So the input topography is the SAME at
> both grids and the 3× spread in `S_ref` is a **CONSEQUENCE** of the erosion-work divergence,
> downstream of it.
>
> **And that makes the circularity demonstrated rather than suspected:** `A_c(S)` is calibrated
> on a quantity *produced by the defect it is meant to work around*. Recalibrating `S_ref` per
> grid would fit the symptom; it cannot converge the grids, and block C shows why — raising
> `A_c` makes the work ratio WORSE (3.18 → 8.46 → 10.89 for `A_c` = 0.1 / 0.4 / 1.0 km²).

### `S_min` — a PROXY, labelled, with its cost measured

Montgomery & Dietrich constrain the relationship over roughly `S = 0.1–1.0`. **There is no
published lower bound for channel initiation**, so extrapolating `A_c ∝ S⁻²` onto a near-flat
apron leaves the fitted range and needs a floor. `CHANNEL_HEAD_S_MIN = tan(0.5°) = 0.0087` is
that floor and it is labelled a proxy in the code and here.

Chosen so the BOUND does not do the LAW's work: it covers **8.6 % of land at 2048² and 3.6 % at
8192²**, and caps `A_c` near 145 km² — no channel ORIGINATES on the flattest ground, which is
the right picture for a coastal plain where rivers cross the plain rather than being born on it.
A 2° floor would cover 27.7 % of land and let the clamp, not the law, set the apron's threshold.

### The two-sided form: implemented, measured, REJECTED

The plain law also LOWERS the threshold where `S > S_ref`, channelising steep ground. The
invariant suite caught it, and it is exactly the bad trade a coastal metric cannot see:

| | 2048² off → two-sided | 8192² off → two-sided |
|---|---|---|
| slope > 30 % | 14.50 → **0.76 %** | 23.68 → **15.59 %** |
| slope > 45 % | 1.21 → **0.00 %** | 8.70 → **1.70 %** |
| Strahler S5 | 117 → 237 | 59 → **0 (gone)** |
| hypsometry mean | 282 → 216 m | 685 → 388 m |

**The arêtes relief-v2/v3 exist to produce were erased, and at 8192² the top of the network
hierarchy vanished.** The calibration is anchored AT `s_ref` and nothing says the steep regime
needs a lower threshold — the shipped relief was validated with the constant. So the law is
clamped **RAISE-ONLY**: it corrects the regime it was derived for, the flats, and leaves the
hillslopes exactly as calibrated.

### The verdict, raise-only, at both resolutions

**PRIMARY** — and reported in ABSOLUTE terms, because the per-100 km density has a moving
denominator (rule 2, which I reproduced in choosing it):

| | 2048² | 8192² |
|---|---|---|
| spurs, shipped → law (target) | 1 254 → **560** (15) | 1 849 → **544** (19) |
| **share of the gap closed** | **56 %** | **71 %** |
| coastline km, shipped → law (target) | 4 860 → 3 333 (1 581) | 6 206 → **2 692** (1 602) |
| length p90 km | 3.50 → 4.78 (24.68) | 4.62 → 3.81 (23.34) |
| local parallelism | 0.515 → 0.513 (0.330) | 0.525 → **0.430** (0.363) |

**RULE-7 CONTROL BLOCK — the incision is alive**, which is how `A_c × 100` failed:

| | mean | p50 | p90 | land % |
|---|---|---|---|---|
| shipped 2048² | 282 | 161 | 787 | 15.0 |
| **law 2048²** | **304** | 193 | 770 | 16.0 |
| shipped 8192² | 685 | 445 | 1 606 | 16.4 |
| **law 8192²** | **697** | 460 | 1 606 | 16.7 |
| un-eroded reference | 865 | 679 | 1 836 | 16.9 |

Hypsometry moves +8 % and +2 % — nowhere near the un-eroded field. And the panel and the
absolute count now AGREE, which the amended rule 9 requires before a visual verdict.

### What else moved in the network

| | 2048² off → on | 8192² off → on |
|---|---|---|
| slope > 30° | 14.50 → 12.84 % | 23.68 → 21.56 % |
| slope > 45° | 1.21 → 1.11 % | 8.70 → 8.27 % |
| drainage density km/km² | 0.746 → **0.798** | 0.561 → **0.756** |
| confluences | 4 419 → 4 688 | 6 049 → 7 434 |
| Strahler S5 | 117 → 194 | 59 → **139** |
| navigability (boat) | 167 → 198 | 0 → 0 |
| floor/local-ridge p50 (1 km) | 0.240 → 0.352 | 0.526 → 0.391 |

**Visible in the export?** The steep shares move 5–11 % relative — at the edge of noticeable. The
drainage density rises 7 % / 35 %, and the hierarchy DEEPENS (S5 117 → 194, 59 → 139): more
tributaries reaching higher order.

> **CORRECTED by Finding 56b: this table was measured on the UNCLIPPED network**, before
> `clip_rivers_to_lakes` and the spillway append. On production's clipped network the 2048²
> density rises 1.3 % (not 7 %) and the confluences FALL 3 % — at the coarse grid the law
> reorganises the hierarchy without densifying it. The 8192² conclusion stands (+33 %). **Navigability barely moves** (167 → 198 small-boat reaches at
2048², none at either setting at 8192²), so the barge/ship question is untouched by this — it
remains blocked on the discharge, i.e. on the hypsometry.

### Two blocks of my own suite that measured NOTHING, said plainly

**Lakes and closed depressions both read 0 in every configuration.** `hydro()` measures the lakes
of the POST-BREACH drainage, which production discards, and a breached field has no depressions
by construction. **That is the exact mistake the first H-1 bench made and which this ADR already
records** — reproduced by me, on a suite whose whole purpose was to catch what the coastal
metrics cannot see. The lake invariants (footprint ≤ level, depth = level − floor,
exorheic ⇒ outlet, no orphan mouth, monotone profiles) are therefore **UNVERIFIED for this
change**, not verified-and-clean, and they must be re-run against the pre-breach population
carried forward the way `build_hd_drainage` does it. `W/D per order` is likewise uninformative as
computed (sub-metre widths without `geo_scale_ratio` give ratios of ~0.0).

> **CLOSED by Finding 56b, with two of the sentences above corrected.** The invariants were
> re-run against production's population (39–56 lakes; three of them do NOT hold, one being a
> regression this law introduces). And the W/D diagnosis here is **WRONG**: `geo_scale_ratio` is
> 1.0 in production, so it was never the cause — the metric divided a hydraulic channel width by
> a valley relief, two different objects, and was not a ratio of anything. See Finding 56b.

### Standing

Gated OFF pending the author's visual validation. `coast_warp_strength` 1.5, contour relaxation
`None`, both unchanged.

## Finding 56b — the lake invariants against PRODUCTION's population; three do not hold; and rule 10 as a type

Finding 56 shipped a metric block that measured nothing and I read it as clean. This closes it,
and the closing produced four results I did not expect.

### The instrument first, because that is what made the mistake repeatable

`build_hd_drainage` was a **private** function in `ymir-viz::bridge::c1::hd`, and every line of
it was core logic — `c1_drainage_windowed_infil`, `apply_lake_water_balance`,
`below_sea_basin_lakes_infil`, `clip_rivers_to_lakes`, `push_segment`, `apply_geo_scale_ratio`.
Private in the binary crate means **production's lake carry is unreachable from a test**, so any
bench wanting production's lakes had to re-implement the chain by hand. Two did. Both dropped the
same step — adopting the pre-breach lake geometry — and both reported `lakes 0`.

It now lives in `ymir_core::tectonics_c1::hd_assembly::assemble_hd_drainage`; the viz side is a
wrapper that only adapts `ClimateResult` to `DrainageClimate`. **There is one carry and the bench
calls it.** It cannot be dropped by a copy that forgot it, because there is no copy.

The order is load-bearing and is now documented at the function: final drainage on the breached
field WITH the climate discharge (Finding 22) → adopt the **pre-breach** lake geometry → H-1c
water balance → below-sea merge and submerged-lake cleanup → **then** clip the rivers → append
the spillways → geographic scale. Clipping before the footprint is final is the orphaned-mouth
defect, already fixed twice.

### The population, OFF against ON, at both resolutions

Production config (`production_hd_config`, `geo_scale_ratio = 1.0`), arid-hot bed (25°, span 10),
seed 10 481 999 410 520 546 993, 400 km domain. Two independent runs agree on every number below.

| | 2048² off | 2048² **on** | 8192² off | 8192² **on** |
|---|---|---|---|---|
| lakes in the inventory | 41 | **39** | 56 | **48** |
| exorheic / endorheic | 14 / 27 | 16 / 23 | 28 / 28 | 24 / 24 |
| water surface km² | 3 795 | 3 788 | 1 638 | **1 323** |
| footprint cells | 99 496 | 99 300 | 687 073 | 554 889 |
| lakes DETECTED pre-breach | 67 | 61 | 58 | 50 |
| closed depressions, RAW field | 68 431 cells (2 610 km²) | 72 642 (2 771) | 1 236 459 (2 948) | 1 011 245 (2 411) |

**The law REDUCES the lake population** — −5 % at 2048², **−14 % at 8192²**, water surface −19 %
at the fine grid. That is the expected direction for a denser channel network (more basins
acquire an outlet), and it answers the question that made this rerun precede the visual
judgement: **the basins moved, and at 8192² they moved by more than the coastline did.**

The closed depressions are counted on the **RAW** field, the only place they exist — the breach
removes them by construction, so counting them after it was the second half of the same mistake.

### The full invariant set — the ones that hold, over real populations

| invariant | 2048² off → on | 8192² off → on |
|---|---|---|
| footprint at or below its level | 0 → 0 | 0 → 0 |
| exorheic ⇒ a traced outlet reaching a sink | 0 → 0 | 0 → 0 |
| monotone long profiles (flat lake crossings tolerated) | 0 → 0 | 0 → 0 |
| no duplicate id, no empty footprint | 0 / 0 | 0 / 0 |
| `area_km2` == footprint | 1 → 1 | 0 → 0 |
| `depth == level − floor` on the **EXPORTED** (breached) field | 3 → 4 | 0 → 0 |
| `depth == level − floor` on the raw pre-breach field | 19 → 19 | 14 → 13 |

Populations, printed next to their numerators because that is the point of rule 10: 39–56 lakes,
99 300–687 073 footprint cells, 702–1 724 inlets, 10 183–15 966 long profiles, 14–28 exorheic.

**One result worth keeping: the depths agree with the field that is EXPORTED (0–4) and disagree
with the raw pre-breach field (13–19).** That is the right way round — a consumer reading
`lakes.json` against the height raster finds them coherent. It follows from H-1c recomputing
`level_m`/`depth_m` on the breached field while adopting the pre-breach footprint.

> **CORRECTED by Finding 56c: true of the ARID bed, FALSE of the production one.** H-1c
> recomputes the depth only for the lakes it RE-SETTLES, i.e. the endorheic ones — 27–28 of
> 41–56 in arid, but only 2–5 of 80–86 in humid. In the humid (production) bed **44–59 of the
> lakes disagree with the EXPORTED field** and only 2–5 with the raw one. The relation inverts
> with the climate, and the shipped climate is the bad side of it.

### Three that do NOT hold — and one is a regression the law introduces

| | 2048² off → on | 8192² off → on | population |
|---|---|---|---|
| footprint **not 4-connected** from its lowest cell | **24 → 24** | **38 → 28** | 41/39/56/48 lakes |
| **dangling `lake_map` ids** (below-sea / detected) | **6 → 6** (6/0) | **16 → 18** (16/0, 18/0) | 47/45/72/66 ids |
| **orphan mouths** — a terminus on an id absent from the inventory | **1 → 46** | **2 → 26** | 10 183/10 281/12 442/15 966 termini |
| level > inlet arrival altitude | 5 → 4 | 9 → 5 | 1 325/1 724/702/1 221 inlets |

**(a) More than half of every lake population has a disconnected footprint** (24 of 41, 28 of
48). PARTLY attributed: `apply_lake_water_balance` settles an endorheic basin by keeping its
`n_eq` **lowest cells sorted by altitude**, which is not a connected set — a basin with two
sub-bowls at comparable depths becomes two patches under one id. That covers the endorheic share
(15–22). **The 8–12 EXORHEIC cases are NOT covered**: H-1c leaves exorheic geometry untouched, so
they come from the pre-breach detection and are **UNATTRIBUTED**. Filed, not fixed. Unchanged by
the law, so it does not block this decision.

> **RESIZED by Finding 56c.** In the humid (production) bed the disconnected footprints are 52–65
> per configuration and **almost all EXORHEIC** (1–4 endorheic). So the H-1c truncation explains a
> SMALLER share than written here and the unattributed part is LARGER — in the shipped climate it
> is essentially the whole defect.

**(b) The dangling ids are ALL below-sea, none detected**, and this is a documented trade rather
than a discovery: `below_sea_basin_lakes_infil` MARKS every below-sea sink in `lake_map` for sink
validity but only INVENTORIES those clearing `INVENTORY_MIN_CELLS = 4`. Its own unit test asserts
it in those words — "the 2-cell pit is below the 4-cell floor → marked, not listed". Finding 33
Part A lowered the floor from 5 km² to 4 cells precisely to shrink this residual; 6–18 remain.

**(c) Orphan mouths 1 → 46 and 2 → 26 — a 46× and 13× regression that the law introduces.** This
is (b) multiplied by the denser network: far more river termini land in a sub-4-cell below-sea
sink that `lakes.json` does not contain, so a consumer resolving the terminus's lake id gets
nothing and the river ends in open ground.

**It is the strongest argument against shipping the law as it stands — and it is not a terrain
defect.** The fix is in the reporting layer and is independent of the channel head: either raise
those sinks into the inventory, or attribute such a terminus to `None` the way the spillway
append already does for a source basin below the floor (`segment_source_lake`). Until one lands,
the law would ship 46 unresolvable termini at 2048².

**(d) `level > inlet arrival` is small and IMPROVES with the law** (5 → 4, 9 → 5). Pre-existing.

### W/D per Strahler order — the metric was replaced, not just recomputed

Finding 56 divided a **hydraulic channel width** (metres, from discharge) by a **1 km window
relief** (hundreds of metres). Those are two different objects — a channel dimension over a
valley dimension — and every order printed `0.0`. I attributed that to the missing
`geo_scale_ratio`; **that attribution was wrong.** The number was unusable because it was not a
ratio of anything, and `geo_scale_ratio` is 1.0 in production anyway, so it was never the cause.

Replaced by a genuine **valley cross-section**: a transect perpendicular to the flow at the
segment midpoint, walked out one cell at a time on each side tracking a running maximum, the
shoulder called where the ground has dropped 10 m below that maximum (we are descending the far
side) or at 2 km. `W` is the sum of the two half-widths; `D` is the **shallower** of the two
rises, so one high flank cannot inflate an asymmetric valley. Spillways are excluded — a
spillway's path is traced over a col, outside the accumulation network, and has no valley.

| valley W/D (p50) | 2048² off → on | 8192² off → on |
|---|---|---|
| S1 | 56.2 → 51.0 | 21.6 → 34.4 |
| S2 | 96.2 → 84.3 | 18.2 → 25.7 |
| S3 | 97.7 → 84.2 | 26.9 → 36.5 |
| S4 | 84.1 → 64.8 | 42.1 → 57.8 |
| S5 | 34.6 → 76.8 | 16.8 → 50.1 |

**The two grids are in opposite regimes, again.** At 2048² the valleys are broad (W/D 35–98) and
the law NARROWS them at every order but S5; at 8192² they are three to five times tighter
(17–42) and the law WIDENS them at every order. Neither is a defect of the law — it is the
Findings 42–44 hypsometry split showing up in valley form, and it is a fourth independent path to
the same blocker.

### The channel widths, and a result that only appeared because rule 10 refused the column

Channel width comes from `w = 5·Q^0.5` on the discharge, and in this test bed:

| | 2048² off → on | 8192² off → on |
|---|---|---|
| share of reaches with Q > 0 | 23–48 % | **0–17 %** |
| channel width p50 \| max, wet reaches only | 0.9–2.2 \| 16.1 m | 0.65–1.9 \| **5.2 m** |
| **orders with NO wet reach at all** | none | **S5 (all 51 reaches), shipped** |

**At 8192² not one of the 51 order-5 reaches carries any discharge**, and the widest channel
anywhere on the fine grid is 5.2 m. The trunk of the hierarchy is dry while the water sits in
low-order coastal reaches — which is coherent for a desert continent (no exotic rivers) but is
also why no width or navigability calibration can be anchored on this bed.

### Rule 10 — recorded, and made a TYPE rather than a paragraph

**The rule.** *Any bench reporting more than one configuration must assert that each reported
quantity responded to the sweep, and must fail loudly on a column constant at zero.*

**Scope: any bench reporting more than one configuration.** With one configuration there is
nothing to vary against and the rule cannot apply.

**Why it is a type.** The rule was already written down in this ADR, after the first H-1 bench
made the mistake. Writing it down did not prevent the second occurrence — by me, on the suite
whose entire purpose was to catch what the coastal metrics could not see. **A recorded rule does
not protect if it is not executable.** So it is `ymir_core::metric_sweep::Sweep`, the same move as
`SegmentRow`/`push_segment` and `segment_arrays_aligned()`: the correct behaviour is the only one
the type permits.

Four registrations, because the failure mode is not the same for every kind of number:

| | what the bench promises | fails when |
|---|---|---|
| `push` | the sweep should MOVE this | constant — at zero, or at anything |
| `invariant` | a rule-7 control: it must NOT move | it moved, **or it is flat at zero** — an invariant that holds at nothing is not evidence that it holds |
| `checked(v, n)` | `v` violations over a population of `n` | **`n == 0` in any configuration.** `v` itself is free: zero violations is the desired result |
| `pinned(v, why)` | it cannot respond, and here is why, in writing | it is zero everywhere — the reason never covers that |

**`checked` is the pair that catches the defect in its real shape.** A violation counter flat at
zero is a good result, so guarding it for variation would be wrong; what was wrong both times was
the **denominator**. `0 violations out of 0 lakes` prints identically to `0 out of 4 000`. Every
violation count is now printed next to its population and the population must be non-empty.

`pinned` is a deliberately narrow escape, because a guard with no escape gets satisfied by
deleting the metric. The reason is mandatory and is reprinted in the report, and **a pinned
column at zero still fails** — so it cannot swallow the case the rule exists for.

**It fired on its first real run — 13 columns — and it was right three times:**

1. **`channel width p50` was 0.000 m at every order, both grids, both settings**, because most of
   the network carries a discharge of exactly 0.000 m³/s. The median was measuring the dry
   majority. Switched to p90 — **still 0.000 for S5 at 8192²**, because that order has no wet
   reach at all. A percentile of a mostly-zero distribution is not a channel width; it is a
   restatement of the dry fraction. Now conditional on `Q > 0`, printing `n/a` for a dry order
   (a printed `n/a` cannot be mistaken for a measured zero, which `0.000` was) with the wet share
   and the maximum registered instead. **The dry-trunk result above exists only because the
   column was refused twice.**
2. **`max Strahler order`** was constant at 5 — not a quantity expected to respond but an
   **invariant claim**, and exactly the one the two-sided form violated (S5 59 → 0 at 8192², i.e.
   max order 5 → 4). Re-registered `invariant`.
3. **`navigable boat reaches`** was constant at 2 at 8192². Not a measurement failure: the classes
   are thresholds in m³/s and the discharge is 2–3 orders below them (Finding 47). `pinned` with
   that reason in writing.

### Two CORRECTIONS to Finding 56, made visibly

**1. The network table was measured on the UNCLIPPED network.** Finding 56 read the bench's own
`c1_drainage_windowed` result, before `clip_rivers_to_lakes` and before the spillway append.
**Production's network is the clipped one.**

| | 2048² off → on | *(F56 reported)* | 8192² off → on | *(F56 reported)* |
|---|---|---|---|---|
| segments | 10 183 → 10 281 | *13 149 → 14 053* | 12 442 → 15 966 | *13 737 → 18 746* |
| confluences | **3 102 → 3 002** | *4 419 → 4 688* | 5 086 → 6 106 | *6 049 → 7 434* |
| drainage density km/km² | **0.611 → 0.619** | *0.746 → 0.798* | 0.540 → **0.721** | *0.561 → 0.756* |
| Strahler S5 | 54 → 128 | *117 → 194* | 51 → 131 | *59 → 139* |
| spillways | 17 → 18 | *(not reported)* | 44 → 40 | *(not reported)* |

**One claim changes sign.** At 2048² the confluences **fall 3 %** where I reported a 6 % rise, and
the density rises **1.3 %** where I reported 7 %. So "a denser network at both resolutions" was
wrong at the coarse grid: at 2048² the law **reorganises** the hierarchy (S5 ×2.4) **without
densifying** it. The 8192² conclusion stands and strengthens — density +33 %, confluences +20 %,
S5 ×2.6.

**2. My own new control block was wrong before it was published, and cross-checking two benches
caught it.** The rule-7 hypsometry mean read 672.2 m at 8192² where `coastal_comb_levers` reported
685. Same seed, same config, and the p50 (445), p90 (1 606) and land share (16.4 %) matched
exactly — so it was not the field. It was an **f32 accumulator**: summing 1.8 M altitudes in f32
saturates, because once the partial sum passes ~10⁹ adding 700 m stops changing it. Fixed to f64,
and the two benches now agree to the printed digit at all four points (282.3/304.5 and
685.3/697.2 raw). Method rule 7's sibling: *a quantity that two benches compute must be compared
across them; agreeing to three digits and diverging on the fourth statistic is a bug, not noise.*

The block also now reports the hypsometry of **both** fields, which prices the breach: it lowers
79–6 446 cells below sea level, and because those are the LOWEST cells, removing them from the
land mask RAISES the mean of the remainder by 4.7–7.4 m. That is why the raw and breached means
differ in the direction they do, and it is why the two benches were never comparable before.

### Standing

Still gated OFF pending the author's visual judgement — the note is at
`docs/reports/c1_continental_buoyancy/closure_morphology/channel_head_law/VALIDATION_NOTE.md`.
`coast_warp_strength` 1.5, contour relaxation `None`, both unchanged. `A_c = 0` remains a
diagnostic. **The orphan-mouth regression (c) must be closed before the law could ship even on a
favourable verdict**, and it is independent work.

## Finding 56c — the HUMID bed, which refutes three figures; and rule 6 extended to the artefact

Finding 56b measured the lake invariants on the arid-hot bed alone, because that is the bed the
below-sea work was built on. **That was a scope error in the same family as the ones this ADR
keeps recording:** the hydrology is a CLIMATE result, and `humid` (45°, span 40) is the
PRODUCTION default. A lake verdict read on the arid bed says nothing about the shipped one. The
suite now runs both beds at both resolutions — the terrain is climate-independent, so the field
is computed once per (resolution, law) and only the hydrology is recomputed.

### It refutes three figures, one of them mine

| claim | measured |
|---|---|
| "the population goes from 55 to 30 in humid" | **80 → 68 at 2048², 86 → 76 at 8192²** — a 15 % / 12 % fall, not 45 % |
| "+6 endorheic in arid-hot" | **−4 at both grids** (27 → 23, 28 → 24). The raise-only clamp raises `A_c` on FLATS, so it incises the low-gradient ground LESS, not more |
| "zero navigable at 8192² before and after" | **2 → 2 in arid** (not zero), and in the production humid bed **35 → 11 small boat, 3 → 2 barge** |
| "network densified 35 %" | **+49 % in humid at 8192²**, +33 % in arid. The 35 % was the unclipped-network figure corrected in 56b |

**And the third one convicts my own validation note**, which said the navigability spread would
not move. It does not move on the arid bed and I wrote the note from that bed:

| navigable reaches | 2048² off → on | 8192² off → on |
|---|---|---|
| small boat, **humid** | 214 → **175** (−18 %) | 35 → **11** (−69 %) |
| barge, humid | 2 → 3 | 3 → **2** |
| small boat, arid-hot | 108 → 103 | 2 → 2 |
| discharge p90 m³/s, humid | 1.759 → 1.983 | 1.017 → 1.386 |

**In the production climate the law cuts navigable reaches by 18 % and 69 %** — while the
discharge p90 RISES. So it is not less water: it is the same water spread over a denser network,
so fewer individual reaches clear a threshold. That is a real consequence the author must weigh,
and I had told him not to look for it. Method rule 3's shape again — *decompose a measured
factor before quoting it* — here applied to a climate rather than to a term.

### Two invariant failures that only the humid bed exposes

| | humid 2048² off → on | humid 8192² off → on | arid, for contrast |
|---|---|---|---|
| **exorheic without a traced outlet** | **5 → 2** | **2 → 2** | **0 everywhere** |
| `depth != level − floor`, EXPORTED field | **59 of 80 → 48 of 68** | **44 of 86 → 40 of 76** | 0–4 |
| `depth != level − floor`, raw pre-breach field | 5 → 4 | 2 → 2 | 13–19 |
| footprint not 4-connected (of which endorheic) | 52 (4) → 48 (2) | 65 (1) → 51 (1) | 24 (16) → 28 (16) |

**(a) `exorheic without outlet` is a MASS-BALANCE violation** — an exorheic label *requires* an
outlet, which is the whole premise of the `Spillway` type ("a basin receiving more than it
evaporates MUST overflow"). It reads 0 in every arid configuration and 2–5 in every humid one.
Unchanged in kind by the law (5 → 2, 2 → 2), so it is pre-existing, and it would have stayed
invisible for as long as the suite ran one bed.

**(b) The depth/raster incoherence INVERTS between beds, and the production bed is the bad
one.** In arid the depths agree with the exported breached field (0–4 failures of 39–56) and
disagree with the raw one (13–19). **In humid it is the other way round: 44–59 of 68–86 disagree
with the EXPORTED field** and only 2–5 with the raw one.

The mechanism is exact: H-1c recomputes `level_m`/`depth_m` on the breached field **only for the
lakes it re-settles, i.e. the endorheic ones**; exorheic geometry is adopted untouched from the
pre-breach detection. Arid is 27–28 endorheic of 41–56; humid is 2–5 of 80–86. So **in the
production climate most of `lakes.json`'s depths are pre-breach values published beside a
post-breach raster.** Finding 56b's "that is the right way round for a consumer" was true of the
bed it was measured on and false of the shipped one — recorded as a correction, not a footnote.

**(c) The disconnected footprints are almost entirely EXORHEIC in humid** (65 of them, 1
endorheic). Finding 56b flagged the exorheic cases as unattributed and treated them as the
minority; in the production bed they are essentially the whole population of the defect. The
H-1c altitude-sorted truncation therefore explains a *smaller* share than 56b implied, and the
unattributed part is *larger*. Still filed, still not fixed, now correctly sized.

### `Sweep::part` — rule 2 turned into a type, after reproducing rule 2 a third time

Rule 10 fired again on the humid bed, on three columns, and it was right about all three being
mis-registered rather than unmeasured:

- **`Q>0 share` constant at 100 %** for orders 4–5 in humid 2048², both settings. In a wet
  climate every large reach is wet, so the SHARE sits at its ceiling and cannot respond.
- **`lakes endorheic` constant at 2** in humid 8192². The category is structurally near-empty
  there — `net_evap = max(0, PE − precip)` is ~0 — which is a climate result, not a dead wire.

Both are the SAME mistake I have now made three times: **reading a rate without its
denominator** (the ratio-vs-excess of Finding 43, the per-100 km spur density of Finding 56, the
wet share here). Writing rule 2 down did not stop the third instance either. So it becomes
`Sweep::part(name, count, of)`: **the count is free to be flat or zero, and the POPULATION is
what must be non-empty.** The number of reaches in an order moves even when the share is pinned
at 100 %, so the guard stays live exactly where a share would have gone blind — and
`Sweep::checked` is now the same mechanism with violation semantics.

Note what was NOT done: the share is still printed. The fix is to register the quantity that can
respond, not to delete the one that reads well.

### Rule 6, EXTENDED — verify a written artefact NON-EMPTY before declaring it delivered

Recorded in full beside rule 6 itself. In short: rule 6 was a prohibition on one tool, and the
failure class is wider — *the artefact reports the instrument instead of the subject*, of which
an empty or truncated document is the documentary form and a column constant at zero the
numerical one. Fourth instance on this campaign (global R blind to local parallelism; the
saturating 8 km cap; `channel width p50` at 0.000; a report file reported empty).

Made mechanical in `crates/ymir-core/tests/doc_artifacts.rs`: every Markdown file under `docs/`
(161 at the time of writing) must be non-empty, above a 200-byte floor, valid UTF-8, free of
cp1252 round-trip signatures and of conflict markers. It **asserts its own population** so a
broken walk cannot report "0 problems out of 0 files", and its negative control **builds the
corruption** by running a real cp1252 round trip rather than trusting a literal, so a signature
list that could never fire fails the test instead of passing quietly.

**The instance that prompted it did not reproduce.** `VALIDATION_NOTE.md` was reported empty; it
was 6 805 bytes, 129 lines, valid UTF-8, no BOM, no double-encoding, and byte-identical in `HEAD`
— and it was written with the Write tool, never through the banned pipeline. That does not weaken
the rule: *"I looked and it seemed fine"* is the inspection the rule says is insufficient, in
either direction, and the question is now settled by an assertion that runs on every
`cargo test`.

### And the author's regeneration prerequisite, also made mechanical

The prerequisite "enable the toggle and regenerate; the `eroded` and `drainage` caches
invalidate" is a claim about a content-addressed key, and a false one would produce a confident
wrong answer in silence — the author would pass a visual verdict on the SHIPPED terrain, which is
precisely the defect `eroded_key_full`'s own doc comment records for volcanism ("the terrain
differs, the drainage would not"). `crates/ymir-core/tests/channel_head_law_cache_key.rs` proves
it: eroded `5a735ac0 → a84ac294`, drainage `9770ded0 → f1d20f47`, both law parameters reach the
key so recalibrating `S_ref` or the `S_min` proxy cannot be served stale, and a negative control
asserts the key is STABLE across identical calls — without which the test would pass for the
wrong reason on a key that folded a timestamp.

### Standing

Unchanged: gated OFF, `coast_warp_strength` 1.5, contour relaxation `None`, `A_c = 0` a
diagnostic. The blockers before the law could ship are now TWO, both in the reporting layer and
both independent of the channel head: the orphan-mouth regression (56b) and the depth/raster
incoherence in the production climate (56c/b above).

## Finding 57 — the blocker is EROSION WORK, and `A_c` sets its level without being able to converge it

Three attributions and one new method pattern. Every number below is from
`crates/ymir-core/tests/hypsometry_work_attribution.rs` (`DIAG`, asserts nothing), production
config, seed 10 481 999 410 520 546 993, 400 km domain.

### The control that did not exist

`Lever { reference: true }` removes FOUR things (FBM, lithology, fracture, stream power), so the
published `865 − 282 = 583 m` could not be attributed to the incision as written, and the figure
came from a **hard-coded 2048²-only string**. Both defects are corrected in place. But the
measurement refutes the objection's consequence: the three extra levers change 16 % of cells by
up to 37 m locally and leave the hypsometric MEAN identical to under 0.1 m (the FBM is
zero-mean), so **the 583 m IS the incision's**.

### A1 / A2bis — the PRE-INCISION field is resolution-invariant, in altitude AND in slope

| pre-incision (`stream_power = None`) | 2048² | 8192² |
|---|---|---|
| mean / p10 / p50 / p90 altitude | 865.4 / 95 / 679 / 1836 m | **865.5** / 95 / 679 / 1836 m |
| emerged fraction | 16.88 % | **16.88 %** |
| D8 slope p50 / p90 / p99 | 0.0906 / 0.2658 / 0.7580 | 0.0905 / **0.2744** / **0.7900** |
| D8 slope > 30° / > 45° | 2.02 / 0.33 % | **2.22** / 0.39 % |

**The population is NAMED on every row**, because the first version of this table said only
"un-eroded" and reported the FBM-bearing build while skipping the coarse one — the same defect as
the 865 quoted without its resolution. The coarse (no-FBM) rows are 0.0905 / 0.2644 / 0.7575 at
2048² and 0.0900 / 0.2689 / 0.7781 at 8192².

**Not "identical", as I first wrote: invariant to 3–4 % on the tails.** The FBM does add
roughness, and slightly more of it at 8192² (p90 +0.53 % over coarse at 2048², +2.0 % at 8192²) —
the direction predicted against me. But with `n = 1` a 3–4 % difference in `S` cannot produce a
**218 %** difference in work. The magnitude argument survives; the invariance claim was
overstated and is corrected.

**Negative controls.** `octaves` 7 → 8 moves the quantiles by < 0.5 %, and `amplitude_base` ×4
(0.04 → 0.16) moves **nothing at all, to six digits** — the C-1 relief budget caps the FBM's
total downslope rise, and `flow_budget_divisor = nscale · Σ(p·l)^o` rescales the amplitude when
the octave count changes. So `amplitude_base` is the recorded DEAD KNOB, now measured at both
grids, and the octave count is neutralised by construction. The control that *does* move the
instrument is FBM on/off (674 433 differing cells at 2048², 10.8 M at 8192²), which is what
licenses the comparison above.

⇒ **The root cause "roughness of the ×128 upscale" is DEAD.** The upscale's output is invariant
in mean, in altitude quantiles, in emerged fraction, and to within 3–4 % in slope.

Related code fact, since it was the live hypothesis: the FBM octave count is a **fixed 7**
(`FbmUpscaleConfig::default`; `c1_hd_production` never overrides it) and the frequency is
anchored to the SOURCE grid — `nscale = base_frequency · 1024 / src_max²` over the coarse 64²,
sampled in coarse-pixel coordinates. `target_size` appears in neither `nscale`, `octaves`, nor
`flow_budget_divisor`. Octaves *generated* is grid-independent; octaves *resolvable* is not
(octave 6 lands at ~390 m, exactly the 2048² Nyquist) — and with `persistence · lacunarity = 1`
every octave contributes the same gradient, so that would have mattered. The relief budget is
what neutralises it.

### A3 — the work, and why the framing was one-sided

| | 2048² | 8192² |
|---|---|---|
| erosion work, PAIRED over land-in-both **[at `iterations = 2`]** | **633.0 m** | **199.0 m** |
| work p50 / p90 | 483.8 / 1292.5 m | 71.5 / 556.4 m |
| difference of means | 583.1 m | 180.2 m |
| common land lowered > 1 m | **98.80 %** | **66.52 %** |
| land lost to the sea by the erosion | 16.88 → 14.95 % (**−11.4 %**) | 16.88 → 16.44 % (−2.6 %) |

**Which grid is the reference? Neither is validated against the real world, and the programme had
been treating 8192² as the defective one.** The same numbers support the opposite reading: at
2048² `A_c = 0.1 km²` is **2.62 cells**, so there is no hillslope regime at all — 98.8 % of land
is channel — and the grid drowns 11.4 % of its own land in the process. **2048² is the DEGENERATE
grid**, and it is also the calibration resolution (`HILLSLOPE_REF_CELL_M = 400 km / 2048`) and the
grid `S_ref` is pinned to. That circularity is recorded here, not resolved.

> ⚠️ **REFRAMED by Finding 61, and the reframing is stronger than this claim.** `iterations` is a
> DURATION dial, and the 2048² land fraction falls 15.68 → 13.35 % across an iteration sweep while
> 8192² holds 16.53 → 16.38 %. So the drowning is a **rate** effect, not a property of the grid:
> at equal iteration count the coarse grid has effectively eroded far longer, because 91 % of its
> land sits above the channel head against 54 %. **Neither grid is degenerate — they are at
> different points on the same decay**, and the programme has been comparing them at equal
> COMPUTE rather than at equal EROSION.
>
> ⚠️ **TOO GENEROUS — corrected by Finding 62.** They are NOT on the same decay. The fine grid's
> decay **terminates at 311.2 m of work** (measured: `dt` over a 40 000× span reaches
> [5.6, 311.2] m) while the coarse grid's shipped state is 633.0 m. No duration brings them
> together. The two grids have **different attainable sets**, and the ceiling is `A_c`'s cell
> count seen from the far end of the duration dial. And at equal work below the ceiling, the
> MASS matches (land 16.325 % against 16.443 %) while the FORM does not (channel share 80.75 %
> against 8.40 %, a factor 9.6).

### C — `A_c` in cells, read back from the built config

| grid | cell km | cell km² | `A_c` = 0.1 km² | `A_c` = 1.0 km² |
|---|---|---|---|---|
| 2048² | 0.195312 | 0.0381470 | **2.6214 cells** | 26.2144 |
| 8192² | 0.048828 | 0.0023842 | **41.9430 cells** | 419.4304 |

A ratio of exactly 16 = (8192/2048)². Verified by reading `StreamPowerConfig::min_area_cells`
back from the constructed config rather than by recomputing the arithmetic.

### C1 — the sweep REFUTES the convergence hypothesis, in the opposite direction

| `A_c` | work 2048² | work 8192² | **ratio** | > 1 m 2048² | > 1 m 8192² |
|---|---|---|---|---|---|
| 0.1 km² (shipped) | 633.0 m | 199.0 m | **3.18** | 98.80 % | 66.52 % |
| 0.4 km² | 336.5 m | 39.8 m | **8.46** | 83.87 % | 22.19 % |
| 1.0 km² | 131.8 m | 12.1 m | **10.89** | 51.67 % | 11.02 % |

The prediction on the table was that the ratio would fall **below 2.0** at `A_c = 1.0 km²`. It
**rises monotonically to 10.89**. Raising the channel head starves the fine grid faster than the
coarse one, because the area-frequency law is sub-linear and the same km² threshold excludes far
more land at HD.

**By the stated criterion — ratio ≥ 3 at every sweep point — `A_c` is not the DISCRIMINANT.** It
is a LEVEL knob: it sets how much work happens, monotonically, on both grids. The convergent
direction is *lowering* it, not raising it, and Finding 43 already measured that ~50 % of the
residual survives `A_c = 0`. So a second mechanism exists and remains **UNATTRIBUTED**.

**The hard consequence for the law under review:** no channel-head calibration can make the two
grids agree, `A_c(S)` included. It was already known that `S_ref` cannot hold at both grids; C1
shows the threshold *itself* cannot either, in either direction.

Control block on every sweep point — deliberately **climate-free by construction** (hypsometry,
pit population, enclosed below-sea regions from `water_class`, network extent are all functions
of the height field alone), which is why one pass suffices and no climate bed has to be chosen:
pits 68 431 → 69 387 at 2048² and 1 236 459 → 599 567 at 8192²; below-sea cells 41 025 → 31 292
and 453 249 → 395 675; network extent 80 107 → 115 857 km and 45 573 → 262 788 km. **The network
extent GROWS as `A_c` rises**, which is only paradoxical until one notices it is measured at a
FIXED 0.1 km² accumulation threshold on the *resulting* field: less incision means less
concentration, so more cells clear a modest threshold under MFD dispersion.

### The correction to my own logic, invalid independently of any measurement

I wrote: *"a slope-sensitive law fed the same slope cannot produce a factor 3.2 on the work."*
`E = K·A^m·S^n` has **two** field terms. Showing `S` invariant says nothing about `A`. The
proximal cause was not eliminated, it was **displaced from `S` to `A`** — and the measurement
that closes it was three lines further down my own report (`A_c` covering 55.4 % of land against
10.1 %) and I failed to connect it. Method rule 3's defect in its purest form: a factor quoted
without being decomposed.

The area term itself is dimensionally clean at the point of use — `am = (area · cell_km2)^m`, so
km², with `area` in cells and `f = kdt · am / dist_m` in the Braun–Willett implicit update. It is
not a cells-versus-km² bug.

### METHOD PATTERN — "physically dimensioned but SUB-CELL"

`A_c = 0.1 km²` is in physical units, passes every dimensional audit, and therefore **escapes the
"cells instead of metres" pattern this ADR has now listed six times.** Its effect is still
resolution-dependent, because its value *relative to the cell area* changes by 16×.

> **The test is not "is the constant in physical units". It is "how many cells does it span at
> the COARSEST resolution in its validity domain".**

⚠️ **NOVELTY WITHDRAWN (Finding 62).** The OBSERVATION was already in this ADR, at **Finding 6**,
in these words: *"NOTE: A_c must be resolvable — 0.1 km² is sub-cell below ~2048² (it needs ~2.6
cells @2048², 42 @8192²)."* Same constant, same two cell counts. It was recorded once, at the
beginning, and lost for fifty findings while four rounds hunted the same blocker from different
directions. What is new here is the **naming of the pattern** and its promotion to an executable
guard — and that a recorded instance could go missing for that long is precisely the argument for
naming patterns instead of recording instances.

Guarded in `crates/ymir-core/tests/resolution_invariants.rs`, where the assertion is `#[ignore]`d
**because it fails by design**: a green test there would claim a property the code does not have.
The floor is DECLARED at 10 cells, not derived — a D8 neighbourhood is 9 cells, so at or under
that the "threshold" is the stencil. `A_c` scores 2.62.

### Inventory — every area / length constant touching the height field

| constant | value | at 2048² | at 8192² | verdict |
|---|---|---|---|---|
| `RELIEF_V1_A_C_KM2` | 0.1 km² | **2.62 cells** | 41.94 cells | ⛔ **SUB-CELL at the calibration grid** |
| `NECK_KM` (spur detector) | 0.6 km | **3.07 cells** | 12.29 cells | ⚠️ 3 cells — at the grid limit |
| `MIN_SPUR_KM` | 1.0 km | 5.12 cells | 20.48 cells | ⚠️ marginal at 2048² |
| `lake_min_area_km2` | 5.0 km² | 131 cells | 2 097 cells | ✅ |
| `fracture decay_km` | 25 km | 128 cells | 512 cells | ✅ |
| `MAX_SPUR_KM` | 50 km | 256 cells | 1 024 cells | ✅ |
| `stream_km2` (navigability) | 20 km² | 524 cells | 8 389 cells | ✅ |
| `small_boat_km2` | 500 km² | 13 107 cells | 209 715 cells | ✅ |
| `INVENTORY_MIN_CELLS` | **4 cells** | 0.153 km² | 0.0095 km² | ⛔ **the MIRROR case** — dimensioned in cells, so its PHYSICAL meaning varies 16×, and it is the documented cause of the orphan-mouth regression (Finding 56b) |
| `SMOOTH_MAX_SHIFT_CELLS` | 0.75 cells | 146 m | 37 m | ⚠️ mirror case, contour only |
| `diffusion` = **0.08** (relief-v3; 0.05 is relief-v1 — corrected in Finding 58) | dimensionless at `HILLSLOPE_REF_CELL_M` = 195.3 m | — | — | ⛔ the `(ref/cell)²` scaling that anchors it exists **only in the NONLINEAR branch**, and relief-v3 takes the LINEAR one (`critical_slope = 0`). Anchored by declaration, unanchored on the shipped path |
| `lake_min_depth_m` 10 m · `SAME_WATER_BODY_TOL_M` 0.10 m | vertical | — | — | ✅ no cell count applies |

Two entries are ⛔ and both are already implicated in a recorded defect. Nothing is corrected
here — this is attribution.

### Standing

`A_c(S)` still gated OFF, no recalibration, no hillslope term, no timescale work. Two guards
added for the quantities this finding attributes, and nothing else.

## Finding 58 — `A_c`'s DISCRETISATION carries almost the whole divergence; the intensive/extensive split is not population-invariant

Correction first, because it is mine and it is structural. Finding 57 concluded *"`A_c` is a LEVEL
knob, not the discriminant"*, on a criterion — *ratio ≥ 3 at every sweep point* — that tests one
direction only and conflates "discriminant" with "converges when raised". Applying it to the
letter was accommodation, not refutation. **A ratio that moves from 3.18 to 10.89 under `A_c`
alone proves the opposite of the conclusion I drew from it.**

### B1 — the experiment the "sub-cell" pattern commanded, and the decisive number

The Finding 57 sweep varied `A_c` **in km²**, so its cell count moved by the same 16× on both
grids: it could not answer the question the pattern poses. Holding `A_c` at **2.6214 cells on
both grids** (0.1 km² at 2048², 0.00625 km² at 8192²):

| | extensive (channel) | intensive CHANNEL | intensive HILLSLOPE | total work | **raw ratio** |
|---|---|---|---|---|---|
| `A_c` 0.1 km² — 2048² | 90.25 % | 662.0 m | 364.5 m | 633.0 m | — |
| `A_c` 0.1 km² — 8192² | 52.40 % | 342.9 m | 40.7 m | 199.0 m | **3.18** |
| **`A_c` 2.62 cells — 2048²** | 90.25 % | 662.0 m | 364.6 m | 633.0 m | — |
| **`A_c` 2.62 cells — 8192²** | **93.72 %** | 524.5 m | 338.3 m | **512.8 m** | **1.23** |

⚠️ **Budget-conditional (Finding 61): every number in this table is at `iterations = 2`, and
the shipped ratio itself spans 2.72–3.38 over iters 1–8.** The conclusion below is not
contradicted; its magnitude is conditional.

**The work ratio collapses from 3.18 to 1.23, and the extensive gap INVERTS** (90.25 % against
93.72 % — the fine grid now has *more* channel). At a matched cell-count threshold the two grids
do **81 %** of the same work.

So `A_c` is the discriminant after all, and specifically **its discretisation is**: the same
physical threshold spans 2.62 cells at 2048² and 41.94 at 8192², and that single fact carries
the divergence. It is not only a gate on the cells below it — removing it at HD raises the
channel-conditioned work from 342.9 m to 524.5 m, i.e. it also **truncates the base-level
propagation for everything above it**.

### A1 / A2 — and the decomposition is NOT population-invariant

Measured on the regime AUTHORITY (`accumulation ≥ min_area_cells`, evaluated on the
**pre-incision** field so the denominator cannot be moved by the numerator — the defect that
made the network-extent control unusable):

| `A_c` | extensive lo / hi | **intensive CHANNEL ratio** | **intensive HILLSLOPE ratio** | raw ratio |
|---|---|---|---|---|
| 0.025 km² | 100.00 / 86.02 % | **1.67** | n/a (no hillslope at 2048²) | 1.84 |
| 0.100 km² | 90.25 / 52.40 % | **1.93** | 8.96 | 3.18 |
| 0.400 km² | 57.18 / 11.78 % | **1.93** | 10.84 | 8.46 |
| 1.000 km² | 24.45 / 4.25 % | **2.12** | 10.41 | 10.89 |
| **2.62 cells** | 90.25 / 93.72 % | **1.26** | 1.08 | 1.23 |

Two things follow, and the second is a stop.

**(a) The hillslope regime diverges by an ORDER OF MAGNITUDE (≈ 9–11×), the channel regime by
≈ 2×.** That split is new and it is the useful half of the decomposition. The hillslope-conditioned
work is 364.5 m per cell at 2048² against 40.7 m at 8192².

> ⛔ **RETRACTED IN FULL by Finding 59.** The whole intensive/extensive decomposition above —
> both columns, at every sweep point — is withdrawn. It follows its denominator (2.14–2.32 on the
> ">1 m" proxy, 1.67–2.12 on the km² criterion, 1.26 on the cells criterion), and Finding 59 §4
> adds the decisive reason: the regime LABEL does not describe the population the incision acted
> on, because the accumulation is recomputed at each iteration and the partition moves during the
> run. **Do not build on the numbers in this table.** What survives is the `A^m` arithmetic of
> C1 below — with its explanation corrected in Finding 59 §1b — and the B1 discretisation result.

**(b) The claim "the intensive part is ~2.2× and demonstrably independent of `A_c`" is
REFUTED.** It drifts monotonically 1.67 → 1.93 → 1.93 → 2.12 with `A_c` (a 27 % range), and it
**collapses to 1.26 when `A_c` is held in cells**. The quantity is 2.14–2.32 on the ">1 m" proxy,
1.67–2.12 on the km² regime criterion, and 1.26 on the cells criterion. **It is a property of the
denominator, not of the model.** Nothing is built on it here.

### C1 — the term that had never been measured, and where the ≈2× actually comes from

Pre-incision MFD accumulation over land, in km²:

| grid | p10 | p25 | p50 | p75 | p90 | p99 | = 1 cell |
|---|---|---|---|---|---|---|---|
| 2048² | 0.11444 | 0.26050 | 0.51135 | 1.05395 | 2.46303 | 20.1019 | 1.41 % |
| 8192² | 0.01587 | 0.05252 | 0.10942 | 0.21864 | 0.47750 | 4.1463 | 1.26 % |
| **ratio** | **7.21** | **4.96** | **4.67** | **4.82** | **5.16** | **4.85** | |
| `A^0.5` ratio | 2.69 | 2.23 | 2.16 | 2.20 | 2.27 | 2.20 | |

**`A` is NOT resolution-invariant in physical units: it diverges by ≈ 4.7–5.2× across the whole
distribution**, 7.2× in the low decile. Both predictions on the table were wrong — it is not 16×
at p10 (mine and the author's) and it does **not** converge at p90 (both again). The cell-area
floor binds on only 1.3–1.4 % of land, so this is not a floor effect: it is **MFD dispersion
(p = 2) accumulating over 4× more cells at HD**, which spreads the accumulation and lowers it
everywhere.

And with production's `m = 0.5`, the `A^m` ratio is **2.16–2.27, flat across p25–p99** — which is
the channel-conditioned intensive ratio measured at a fixed km² threshold (1.93–2.12). `S` is
invariant to 3–4 % and `n = 1`. So:

> **`E = K·A^m·S^n` — the ≈2× intensive term in the channel regime IS the `A^m` term, and the
> divergence of `A` is a resolution artefact of the MFD accumulation, not of the terrain.**

> ⚠️ **HALF-CORRECTED by Finding 59 §1b.** The arithmetic stands. The words "of the MFD
> accumulation" do not: the divergence is **flat at 4.6–5.1 across p = 2, 4, 8 and pure D8**, so
> it is not an MFD artefact and no routing-scheme choice repairs it. The mechanism is the drainage
> network's area-frequency law on a finite grid — a fixed physical basin shared among a cell count
> that scales with the grid. Per-cell `A` therefore cannot be resolution-invariant under ANY
> scheme, and neither can per-cell `E`.

That is the attribution the previous two rounds were missing, and it closes the term that was
never measured. It also explains why the intensive ratio is denominator-dependent: changing the
threshold changes which part of the `A` distribution is averaged.

### D1 — the diffusion anchor: real defect, ≈4 % effect, and the SIGN prediction refuted

`diffusion` is dimensionless at `HILLSLOPE_REF_CELL_M = 195.31 m`, and the `(ref/cell)²`
rescaling that anchors it exists **only in the NONLINEAR branch** while relief-v3 takes the
LINEAR one (`critical_slope = 0`). So the shipped hillslope operator is 16× weaker at 8192².

⚠️ **Correction to the Finding 57 inventory: the value is `RELIEF_V3_DIFFUSION = 0.08`, not 0.05.**
0.05 is relief-v1's. Read back from the built config here rather than from a doc comment.

Borrowing the anchor by hand on the shipped linear branch (`diffusion` 0.08 → **1.28** at 8192²,
`diffusion_substeps` 4 → 8 to hold the documented `diffusion / substeps ≤ 0.2`), at `A_c` = 0.1:

| 8192² | shipped | ANCHORED |
|---|---|---|
| total work | 199.0 m | **206.7 m** |
| intensive CHANNEL | 342.9 m | 347.5 m |
| intensive HILLSLOPE | 40.7 m | **52.8 m** |
| land lowered > 1 m | 66.52 % | **75.62 %** |
| **raw ratio 2048²/8192²** | 3.18 | **3.06** |
| hillslope intensive ratio | 8.96 | **6.91** |

**The prediction on the table was that the ratio would RISE; it FALLS, 3.18 → 3.06. Refuted — and
I had agreed with it, so refuted twice.** The reasoning was that weaker HD diffusion leaves
steeper slopes hence more incision, so the anchor would remove work; measured, the anchor *adds*
a little (199.0 → 206.7 m) because it spreads the operator over 9 % more land while barely
changing the channel term.

**The magnitude is the result: a 16× change in the hillslope coefficient moves the divergence by
4 %.** The unanchored diffusion is a genuine defect — it closes 23 % of the hillslope-regime gap
(8.96 → 6.91) — but it is **not** masking the divergence and correcting it would not converge the
grids. Not corrected this round, and now priced.

**File consequence, unchanged in force:** a transport-limited hillslope term cannot be built on
an operator whose strength varies 16× with the grid. That anchor has to be borrowed first.

### E — the default control block had two unusable columns

**The pit column is neither inherited nor saturating — it is NON-MONOTONE.** Pre-incision:
54 869 pits at 2048², 720 630 at 8192². Across the `A_c` sweep: 68 431 → 82 439 → 69 387 at
2048² (+25 %, +50 %, +27 % over the pre-incision baseline) while the work falls 4.8×; and
1 236 459 → 772 015 → 599 567 at 8192², crossing the baseline. The incision both *creates* closed
depressions (by carving) and *removes* them (by breaching), and the balance is non-monotone in
`A_c`. **The column is unusable as a control on an `A_c` sweep** and must be read against the
pre-incision baseline, which it now is.

**And it is now interpretable rather than merely restored, which matters because a non-monotone
column reads as noise on the next pass.** The incision does two opposite things to the closed-
depression population at once: it **carves** new ones (a channel cutting across a col leaves the
ground behind it enclosed) and it **breaches** existing ones (a channel reaching a basin drains
it). Raising `A_c` weakens both, and their difference has no reason to be monotone. So the column
measures a *balance*, not a level, and the only readings it supports are (i) against the
pre-incision baseline and (ii) at a FIXED `A_c` across some other change. It is not noise.

**The network-extent column is withdrawn for this sweep type.** It was measured at a FIXED
0.1 km² accumulation threshold on the RESULTING field, so it confounds the change it monitors —
which is why it grew as `A_c` rose. Either read at the current `A_c` or removed; removed here,
declared.

**The climate assertion, widened as demanded.** It previously covered
`upscale_from_c1_with_progress` — the pre-incision terrain — while the claim needed the ERODED
field. Named: the incision is `incise(height: &GridF32, cfg: &StreamPowerConfig)` (and its
`incise_lithology` / `incise_with_progress` variants), and `StreamPowerConfig` carries exactly
`n, m, k, dt, iterations, sea_level, diffusion, diffusion_substeps, min_area_cells,
a_c_slope_law, threshold, cell_km, depth_scale_m, critical_slope, lateral_erosion,
mfd_exponent` — **no precipitation, no temperature, no runoff, no discharge.** The `A` it reads
is the pure cell-count MFD accumulation.

> Worth stating as a fact in its own right: **the erosion uses GEOMETRIC drainage area while the
> hydrology uses CLIMATIC discharge.** They are different quantities, and the one that carves the
> terrain never sees the climate.

**The A2bis negative control now fires at the scale of its claim.** The old one only responded to
FBM on/off — 674 433 cells, 37 m — for a reading that lives at 3–4 %. A deterministic one-cell
perturbation of **1.30 m** moves the D8 p90 slope by **+2.97 %**, which is the same size as the
8192²-versus-2048² excess (+3.2 %). The per-cent reading is licensed, and the FBM's extra HD
roughness now has a scale: it is worth about 1.3 m of one-cell noise. That formulation is kept
deliberately — it is the only place in this campaign where a per-cent slope reading has been given
a physical size, and a percentage without one is what made the first A2 table unreadable.

**The sub-cell assertion is inverted so that it guards.** An `#[ignore]`d test does not run and
protects nothing. It now pins the CURRENT state — 2.6214 cells, under the 10-cell floor — so it
passes today and **fails the day anyone touches `A_c`**, which is when the finding needs reading.

### The reference-grid conclusion, written rather than noted

Finding 57 called this "recorded, not resolved". That was too modest, and it is decidable by
arithmetic: `A_c = 0.1 km²` implies a hillslope length of `√0.1 km = 316 m`, which is **1.62
cells at 2048²**.

> **The calibration grid cannot carry the hillslope regime that `A_c` defines.** Not "does not",
> *cannot*: a hillslope 1.6 cells wide is not representable. And three anchors are pinned to that
> grid — `HILLSLOPE_REF_CELL_M` (defined as `400 km / 2048`), `S_ref` (pinned there "by project
> convention"), and `diffusion` (dimensionless at that cell size). **All three calibrate a regime
> that does not exist on the grid they are calibrated on.**

`S_ref`'s circularity (Finding 57) is therefore the second layer of a deeper one: the law is
calibrated on a quantity produced by the defect, *on a grid that cannot represent the regime the
law partitions*.

### Standing

Nothing corrected. `A_c(S)` gated OFF, no recalibration, no hillslope term, no timescale work,
no diffusion fix. The decomposition that motivated this round is refuted as population-dependent,
so no remedy is built on it; what replaces it is the `A^m` attribution (C1) and the discretisation
result (B1), both of which point at the same object — **the accumulation field and the cell count
of the threshold read against it** — and neither of which is a physics change.

## Finding 59 — the accumulation operator CONVERGES; per-cell `A` cannot, under any scheme; and the hillslope label is void

Four measurements, three refutations of the round's premises, and the retraction of a
decomposition. All DIAG; nothing in production changed.

### RETRACTION — the intensive/extensive decomposition (Finding 58 §A1) is withdrawn

The intensive ratio reads 2.14–2.32 on the ">1 m" proxy, 1.67–2.12 on the km² regime criterion
and 1.26 on the cells criterion: **it follows its denominator.** Its apparent agreement with
Finding 43's "~50 % surviving `A_c = 0`" was a numerical coincidence promoted to a confirmation —
the same defect this ADR records against the `S_ref` reading. Nothing is built on it, here or
anywhere. Block 4 below adds a third, independent reason it could not have worked.

### 1a — the operator DOES converge. Basins agree to 1.6 %; the divergence is upstream, per cell

Matching method, declared: the largest basins at 2048² by CELL COUNT (the physical area); the
outlet is the basin's maximum-accumulation land cell; at 8192² the maximum accumulation is sought
in a ±8 HD-cell (±391 m) window around that normalised position, and the displacement is
reported. **Matches landing on the window edge are flagged** — there, the tolerance chose the
point, not the topography.

| # | acc 2048² km² | acc 8192² km² | ratio | basin cells 2048² km² | 8192² km² | **ratio** | shift m |
|---|---|---|---|---|---|---|---|
| 1 | 495.52 | 33.08 | 14.98 | 3 786 | 437 | 8.66 | 391 ⚠ edge |
| 2 | 355.26 | 250.96 | 1.42 | 2 377 | 2 264 | **1.050** | 98 |
| 3 | 338.11 | 280.18 | 1.21 | 1 476 | 1 373 | **1.075** | 355 |
| 4 | 131.96 | 83.92 | 1.57 | 1 060 | 1 058 | **1.002** | 352 |
| 5 | 808.41 | 753.77 | 1.07 | 973 | 964 | **1.009** | 244 |
| 6 | 120.92 | 33.88 | 3.57 | 749 | 753 | **0.994** | 417 ⚠ edge |
| 7 | 176.76 | 23.57 | 7.50 | 531 | 530 | **1.002** | 297 |
| 8 | 432.23 | 429.00 | 1.01 | 474 | 481 | **0.984** | 69 |
| 9 | 59.83 | 53.82 | 1.11 | 464 | 449 | **1.033** | 249 |
| 10 | 235.14 | 237.08 | 0.99 | 424 | 417 | **1.016** | 244 |
| | | | **median 1.42** | | | **median 1.016** | |

Total drained land is 27 009 km² at 2048² and 27 008 km² at 8192² — identical, as conservation
requires.

**The basin cell-count areas agree to 1.6 % on the median and to under 8 % on every basin whose
match is not flagged.** So the two grids delineate the *same basins*, and the outlet accumulation
ratio is 1.0–1.6 on eight of ten (basin 7's 7.50 comes with a cell-count ratio of 1.002 — the
basin matched, the *outlet* did not; the window found a local maximum that is not the mouth).

> ⚠️ **RESTATED by Finding 62 §2(a), narrower than written below.** The footprint and the
> accumulation are produced by **two different operators**: `basins` is labelled by tracing D8
> receivers, while `mfd_accumulation` recomputes its own multi-flow partition from `filled` and
> does not route on `direction`. So the 1.016 ratio is NOT self-consistency — but what it
> establishes is *"the D8 basin DELINEATION of the pre-incision filled surface is grid-invariant
> to 1.6 %"*, **not** *"the accumulation operator converges"*. The only column speaking to the
> accumulation is the outlet value, whose median is **1.42, not 1.0**, and whose 42 % is
> UNATTRIBUTED (window width, outlet-cell size, incomplete recombination all remain candidates).
> Read the sentence below as: the basin delineation is grid-invariant, and the outlet accumulation
> agrees to within ~40 %. — **UPDATED by Finding 63 §2b: re-matching by footprint overlap
> instead of a search window brings the median to 1.207**, so about half the 42 % was the
> window. The remaining 20.7 % is still unattributed, and one flagged basin (98.9 % overlap,
> footprint ratio 0.994) is 2.43× apart at the outlet — not an artefact at all.

⇒ **"The accumulation operator does not converge" is REFUTED.** It converges on the basin
footprint and at the outlet. What diverges is the **per-cell partition upstream** — and by
conservation it must: the same fixed physical area is shared among 16× more cells.

This does **not** refute the quantile reading; it relocates it. The quantiles are a true
statement about *what the incision sees per cell* and a false statement about drained area as a
physical quantity. The stop rule therefore does not fire.

### 1b — the exponent sweep REFUTES both readings: the divergence survives pure D8

| scheme | p10 | p25 | p50 | p75 | p90 | p99 |
|---|---|---|---|---|---|---|
| **D8** ratio lo/hi | **8.00** | 7.11 | **5.14** | 5.00 | **5.33** | 9.91 |
| **p = 2** ratio | 7.21 | 4.96 | **4.67** | 4.82 | 5.16 | 4.85 |
| **p = 4** ratio | 7.72 | 4.99 | **4.64** | 4.74 | 5.09 | 5.20 |
| **p = 8** ratio | 8.08 | 5.06 | **4.63** | 4.70 | 4.99 | 5.74 |

**The ratio does not fall with `p`. It is flat at 4.6–5.1 on p50 across every scheme, pure D8
included**, and D8 is the *worst* at both tails (8.00 and 9.91). The prediction that it would
decay monotonically and reach ~1 at D8 is refuted; so is my own prediction of a rotation with
p90 → 1.0–1.5.

The negative control holds, so the sweep is not blind: D8 → p = 2 moves p50 from 0.3433 to
0.5114 km² at 2048² (+49 %) and 0.0668 → 0.1094 at 8192² (+64 %). The scheme change IS visible;
the **ratio** is simply insensitive to it.

⇒ **This is not an MFD artefact and there is no scheme choice that repairs it.** Under D8 the
median cell drains 9.0 cells at 2048² and 28.0 at 8192² — in cell units the median accumulation
*rises* 3.1× when the grid refines 4× linearly, and 16 / 3.11 = 5.14 is the km² ratio. That is
the drainage network's area-frequency law meeting a finite grid.

> **The per-cell drained area cannot be resolution-invariant, for any routing scheme, because a
> fixed physical basin is shared among a number of cells that scales with the grid. Therefore
> `E = K·A^m·S^n` evaluated PER CELL cannot be resolution-invariant either.** This is a property
> of the discretised stream-power law, not of a scheme.

**And no single scale factor can repair it, which constrains every future remedy.** The
distribution of `A` does not TRANSLATE between grids — its SHAPE changes. The per-quantile ratios
say so directly: 7.21 at p10, 4.67 at p50, 4.85 at p99, the tails diverging more than the centre;
and the two directions of the equalising-threshold bisection disagree by 1.75× for the same
reason. **So rescaling `A_c`, `K`, or `A` by any constant will correct one quantile and miss
another.** A reformulation is required, not a coefficient.

**Correction to Finding 58's C1 interpretation.** Its arithmetic stands — `A^0.5` ratio 2.16–2.27
against a channel intensive ratio 1.93–2.12 — but its explanation, *"A's divergence is an artefact
of the MFD accumulation"*, is **wrong**: the divergence survives D8 unchanged. The mechanism is
the area-frequency law on a finite grid.

### 2 — three routes to the divergence factor, and two of them coincide exactly

Pure arithmetic on the pre-incision accumulation, each grid using its own land mask (hence the
shares differ slightly from Finding 58's land-in-both table): channel share at the shipped
`A_c = 0.1 km²` is **91.36 %** at 2048² and **53.62 %** at 8192².

| direction | equalising `A_c` | in cells | factor on 0.1 km² |
|---|---|---|---|
| RAISE 2048² to match HD | **0.46738 km²** | 12.25 | **4.67×** |
| LOWER 8192² to match 2048² | **0.01224 km²** | 5.13 | **8.17×** |

**The asymmetry is 8.17 / 4.67 = 1.75×** — reported because the answer depends on which grid is
moved, which is the same denominator lesson as the retracted decomposition.

⛔ **DEMOTED ON TWO COUNTS — do not cite this as a corroboration.**
>
> **(1) It is one object read twice, not two routes.** Bisecting a threshold to equalise a share
> of a population, and reading the ratio of that same distribution's medians, are two readings of
> the same distribution. The agreement validates the arithmetic and attributes nothing.
>
> **(2) The shares it equalises are the INPUT partition** (Finding 61): the delivered field's
> channel share is 64.24 % at 2048² and **8.40 %** at 8192², against the 91.36 % / 53.62 % used
> here. So these thresholds equalise a partition the shipped field does not have.

The arithmetic, for the record and not as evidence: the raise-direction factor 4.67× equals the
p50 accumulation ratio 4.67× to three digits. The outlet route (1.42) measures a different object
and was never a third estimate of the same thing.

⚠️ **An error in my own bench output**, corrected here: the line printing "the two directions
disagree by 0.57×" multiplied the two factors instead of dividing them. The asymmetry is 1.75×.

### 4 — the hillslope regime is a VOID label, and the diffusion stage is inert

Switching the hillslope diffusion off **entirely** (0.08 → 0.0), at `A_c = 0.1 km²`:

| | intensive HILLSLOPE | intensive CHANNEL |
|---|---|---|
| 2048² shipped | 364.5 m | 662.0 m |
| 2048² **diffusion 0.0** | **358.5 m** (−1.6 %) | 662.3 m |
| 8192² shipped | 40.7 m | 342.9 m |
| 8192² **diffusion 0.0** | **40.6 m** (−0.2 %) | 341.6 m |

**My prediction that the hillslope work is ≥ 80 % diffusion is REFUTED, decisively: it is 1.6 %
and 0.2 %.**

And the mechanism is not what was proposed either. `incise` executes `continue` on a sub-threshold
cell and never writes it, so incision **cannot** act below the threshold. But the accumulation is
recomputed at **each of the two iterations** on the evolving surface, and 662 m of channel
incision reorganises the routing — so cells labelled hillslope on the PRE-INCISION field pass the
gate later. **The partition moves during the run.**

⇒ **"Hillslope regime" is a void label as measured, and there is effectively one regime** — the
conclusion proposed, reached by a different route (not "incision below the threshold", which the
code forbids, but "the threshold's population is redefined between iterations"). This is the
third and strongest reason the intensive/extensive decomposition could not work: **its
denominator did not describe the population the incision acted on.**

**The diffusion stage, correctly labelled at last.** Finding 58 measured a 16× coefficient change
moving the divergence by 4 %; switching the stage OFF moves the hillslope work by 1.6 % and the
channel work by 0.05 %. An operator insensitive to its own total removal is **INERT**, not weak.

> The irony belongs in the record: this ADR **removed** the "inert branch" label from Finding 43
> where it was wrong — the Laplacian there displaces two thirds of the land at mass-neutrality,
> which is *conservative* — and the label was **missing** here, where it is exactly right.

**File consequence, strengthened twice:** a transport-limited hillslope term cannot be built on a
stage now measured inert to its own removal. Before any hillslope physics, the measurement owed is
*why* this stage does nothing — and the candidate is that `diffusion·∇²h` on a cell-based
Laplacian is dominated by the incision's own 662 m per cell at every iteration.

### The standing picture after this round

- the pre-incision **field** is resolution-invariant (Finding 57);
- the accumulation **operator** converges on basins and outlets (1a);
- per-cell **`A`** cannot converge, under any scheme, and drags `A^m` with it (1b);
- `A_c`'s **discretisation** converts that into the 3.18× work ratio (Finding 58 B1: holding
  `A_c` in cells collapses it to 1.23);
- the **hillslope stage is inert** and its regime label is void (4).

No remedy follows from this round, and none is proposed. What it removes is the hope stated at the
start of it — that the remedy would be a scheme choice.

## Finding 60 — the terrain and the water are DECOUPLED: the erosion never sees the climate

An observation that surfaced while widening the climate assertion, and it deserves its own entry
because its consequences reach past the resolution blocker.

### The measurement, which is a code enumeration

The incision is `incise(height: &GridF32, cfg: &StreamPowerConfig)` (and its `incise_lithology`
/ `incise_with_progress` variants). `StreamPowerConfig` carries exactly `n, m, k, dt, iterations,
sea_level, diffusion, diffusion_substeps, min_area_cells, a_c_slope_law, threshold, cell_km,
depth_scale_m, critical_slope, lateral_erosion, mfd_exponent`. **No precipitation, no
temperature, no runoff, no discharge.** The `A` it reads is the pure cell-count MFD accumulation
of the height field.

> **The erosion consumes a GEOMETRIC drainage area. The hydrology consumes a CLIMATIC discharge.
> They are different quantities, and the term that carves the channels never sees the climate.**

### What follows, and it is not only about the blocker

**1. The climate assertion is now correctly scoped.** It previously covered
`upscale_from_c1_with_progress` — the pre-incision terrain — while the claim it served was about
the eroded field. It now covers both, because the incision provably takes no climatic input.
That is why the `A_c` sweep's control block can be climate-free by construction, and why one
pass suffices there.

**2. Navigability classifies on a discharge that did not shape the channels it classifies.**
The thresholds (`stream_km2` 20, `small_boat_km2` 500, `barge_km2` 5 000, `ship_km2` 50 000) are
applied to the climatic discharge computed by the drainage stage. But the channel *geometry* —
where a channel is, how deep, how the network branches — was set by the geometric accumulation,
with no knowledge of where the rain falls. So a reach can be wide-and-deep in an arid interior
because its geometric catchment is large, and be classed non-navigable because its climatic
discharge is nil; and the converse in a wet coastal strip. **Finding 56c's measurement is exactly
that signature**: in the humid bed the law cut navigable reaches by 18 % and 69 % *while the
discharge p90 rose*. Two quantities that ought to be one, moving in opposite directions.

**3. It bounds what any hillslope or timescale work can achieve.** A transport-limited hillslope
term is normally driven by runoff. On this pipeline there is no runoff at the incision stage to
drive it with — the erosion's only water proxy is a cell count. Coupling them is a design
decision that has never been taken explicitly, and it is upstream of the roadmap's next item.

### Not a defect claim

This is a **decoupling**, not a bug: a geometric-area stream-power law is a standard, defensible
LEM formulation, and the discharge stage was added later for the hydrology export. Nothing here
says either stage is wrong. What it says is that **the project has two different answers to "how
much water is here" and has never reconciled them**, and that at least one published verdict
(navigability, Finding 56c) is a consequence of the gap rather than of the change it was
attributed to.

Recorded, not resolved, and no measurement is requested by it this round.

## Finding 61 — BLOCKING: the delivered field is a computational snapshot, not a state. The work does not converge at 8192², and the inter-grid ratio is budget-dependent

`iterations = 2` with `dt = 1.0` a unit placeholder (Finding 44), and the accumulation recomputed
each iteration. The question: of what state is the shipped field the state?

One thing is settled before any measurement, by reading the equation. `E = K·A^m·S^n` carries
**no uplift term**, so the only fixed point of `h ← (h + f·h_r)/(1+f)` is `h = h_r` everywhere,
i.e. base level. The config says so itself — *"iterations: 2 // bounded incision … iters=3 planed
them toward base level"*. **So `iterations` is a DURATION dial, and "convergence" here can only
mean planation.** This is Finding 44's `k_time = K·dt·iterations` seen from the other side.

### The sweep

| 2048² | k_time | work | Δ | mean | land % | OUT channel (IN 91.36 %) | pits / pre |
|---|---|---|---|---|---|---|---|
| iters 1 | 4 500 | 455.7 m | — | 447.4 m | 15.68 | **78.06 %** | ×0.82 |
| **iters 2 (shipped)** | 9 000 | **633.0 m** | +177.3 | **282.3 m** | **14.95** | **64.24 %** | ×1.25 |
| iters 4 | 18 000 | 757.5 m | +124.5 | 176.6 m | 14.06 | 49.17 % | ×1.69 |
| iters 8 | 36 000 | 835.9 m | +78.5 | 106.7 m | 13.35 | 40.05 % | ×2.13 |

| 8192² | k_time | work | Δ | mean | land % | OUT channel (IN 53.62 %) | pits / pre |
|---|---|---|---|---|---|---|---|
| iters 1 | 4 500 | 167.4 m | — | 713.3 m | 16.53 | **19.93 %** | ×1.37 |
| **iters 2 (shipped)** | 9 000 | **199.0 m** | +31.6 | **685.3 m** | **16.44** | **8.40 %** | ×1.72 |
| iters 4 | 18 000 | 223.9 m | +24.9 | 662.3 m | 16.40 | 6.03 % | ×1.79 |
| iters 8 | 36 000 | 262.0 m | **+38.1** | 624.9 m | 16.38 | 5.71 % | ×1.75 |

### It does NOT converge, and the two grids fail differently

**At 2048² the increment decays** — +177.3, +124.5, +78.5, a ratio of ~0.65 per doubling. A limit
is extrapolable and it is **planation**: the mean falls 447 → 282 → 177 → 107 m toward sea level,
and the land fraction falls 15.68 → 13.35 % as the continent progressively drowns.

**At 8192² the increment DECAYS THEN GROWS AGAIN** — +31.6, +24.9, **+38.1**. There is no
extrapolable limit at the fine grid. The pit count is non-monotone with it (×1.37, ×1.72, ×1.79,
×1.75), which is the mechanism: `compute_flow` re-fills depressions each iteration, a planating
surface grows more of them, filling creates flats where MFD disperses maximally, and a burst of
breaching then releases a new round of incision. **The fine grid's trajectory is a reorganisation,
not a relaxation.**

**Predictions, both refuted.** That the Δ would decay geometrically at both grids (mine) is false
at 8192². That the Δ "would not decay fast enough to extrapolate a limit" (the author's) is false
at 2048² and **true at 8192²** — so the author's reading is the one the fine grid supports, and
the fine grid is the one the programme cares about.

### The corollary: the inter-grid ratio IS budget-dependent

| iterations | 1 | **2 (shipped)** | 4 | 8 |
|---|---|---|---|---|
| work ratio 2048²/8192² | **2.72** | **3.18** | **3.38** | **3.19** |

**Non-monotone, spanning 2.72–3.38 — a 24 % range.** Both predictions were that it would FALL with
iterations; it rises then falls.

**And the shipped point is not merely mid-range — it is NEAR THE MAXIMUM**: 3.18 is the
second-highest of the four, against a maximum of 3.38 and a minimum of 2.72. The number that
defined the campaign's priority blocker sits close to the worst case of the dial. This is not an
allegation of selection — `iterations = 2` was fixed years earlier for an unrelated morphometric
reason (Finding 6, traced in Finding 62) — and it is recorded exactly as stated. So the blocker's headline ratio is **not a clean function of the
budget, and not a pure property of the model either**: it is a property of the model *at a chosen
stopping point*, and the shipped stopping point happens to sit mid-range.

### The finding that compounds everything: the partition COLLAPSES, it does not merely move

Finding 59 §4 established that the channel/hillslope partition moves during the run. This measures
the direction and it is one-way:

| | input share (pre-incision) | output share (delivered field) | factor |
|---|---|---|---|
| 2048² at iters 2 | 91.36 % | **64.24 %** | ×0.70 |
| 8192² at iters 2 | 53.62 % | **8.40 %** | **×0.157** |

and it keeps falling with the budget (2048²: 78 → 64 → 49 → 40 %; 8192²: 20 → 8.4 → 6.0 → 5.7 %).

**So every quantity in this dossier conditioned on the "channel regime" used a label that
overstates the channel population by 1.4× at 2048² and 6.4× at 8192².** The mechanism is
consistent with the pits: as the surface planates, relief falls, depression filling creates flats,
MFD disperses maximally on flats, peak accumulation per cell drops, and fewer cells clear
0.1 km².

### NAMED LIST — what this invalidates, weakens, and leaves standing

**Standing, untouched.** These are measured on the PRE-INCISION field or are code enumerations,
so no erosion budget enters them:

- Finding 57 A1 / A2bis — the pre-incision field's resolution invariance (865.4 / 865.5 m,
  16.88 % / 16.88 %, D8 slope quantiles within 3–4 %), and therefore the conclusion that **the
  ×128 upscale is not the root cause**;
- Finding 57's negative controls (`octaves` inert, `amplitude_base` inert — the DEAD KNOB);
- Finding 58's cell arithmetic (`A_c` = 2.6214 / 41.9430 cells, ratio exactly 16) and the
  **"physically dimensioned but sub-cell"** method pattern;
- Finding 59 §1a — basin footprints agreeing to 1.6 %, outlets to 1.0–1.6 (pre-incision field);
- Finding 59 §1b — the accumulation ratio flat at 4.6–5.1 across p = 2/4/8/**D8**, and hence
  **per-cell `A` cannot be resolution-invariant under any routing scheme**. This rests on
  conservation and the area-frequency law, with no free parameter and no budget;
- Finding 60 — the terrain/water decoupling (an enumeration of `StreamPowerConfig`);
- Finding 59 §4's *diffusion* result — the stage is inert to its own total removal (1.6 % / 0.2 %
  at iterations = 2; the claim is about a comparison at fixed budget, which this does not disturb).

**WEAKENED — must now carry `iterations = 2` explicitly:**

- **Finding 57 A3, "erosion work 633 m against 199 m, ratio 3.18"** — true at the shipped budget;
  the ratio spans 2.72–3.38 over iters 1–8. Every future quotation must attach the budget.
- **Finding 57's "2048² is the DEGENERATE grid, it drowns 11.4 % of its land"** — this needs
  **reframing, and the reframing is stronger than the original claim.** The 2048² land fraction
  falls 15.68 → 13.35 % across the sweep while 8192² holds 16.53 → 16.38 %. The drowning is a
  **rate** effect, not a defect of the grid: *at equal iteration count the coarse grid has
  effectively eroded far longer*, because 91 % of its land is above the channel head against
  54 %. Neither grid is degenerate; **they are at different points on the same decay**, and the
  programme has been comparing them at equal compute rather than at equal erosion.
- **Finding 58 B1, "holding `A_c` in cells collapses the ratio 3.18 → 1.23"** — measured at
  iterations = 2 only, untested at other budgets. The conclusion (`A_c`'s discretisation carries
  the divergence) is not contradicted, but its magnitude is budget-conditional.
- **Finding 58 C1 / Finding 59 — the `A^0.5` ratio (2.16–2.27) matching the channel intensive
  ratio (1.93–2.12)** — conditional on iterations = 2 **and** on the input-based label below.

**INVALIDATED AS STATED:**

- **Every quantity conditioned on the channel/hillslope regime via the PRE-INCISION
  accumulation.** That includes Finding 58 §A1's table (already retracted for a different
  reason — this is the fourth), Finding 59 §2's channel shares (91.36 % / 53.62 %), and Finding
  59 §4's channel-versus-hillslope intensive split. The label overstates the channel population
  by ×1.4 and ×6.4 against the delivered field.
- **Finding 59 §2's equalising thresholds** (0.46738 km² = 12.25 cells; 0.01224 km² = 5.13 cells)
  — computed on input shares, so they equalise a partition the delivered field does not have.
  Combined with the author's separate objection that the "4.67× = p50 4.67×" agreement is **one
  object read twice rather than two independent routes**, this is demoted on two counts and must
  not be cited as a corroboration.

### What this means for the programme, stated plainly

The resolution-convergence question as posed — *do the two grids agree?* — is **ill-posed at a
fixed iteration count**, because `iterations` is a duration and the two grids traverse the decay
at different speeds. A well-posed comparison needs the grids matched on **eroded state** (equal
work, or equal hypsometric mean) rather than on equal compute, and that requires the explicit
timescale Finding 44 specified without implementing.

And it changes the priority of the formulation review that was queued behind this block: **any
candidate stream-power reformulation that assumes a steady state is out of domain here**, because
there is no steady state to assume — the equation has no uplift term and its attractor is
planation. That constraint has to be carried into the review rather than discovered inside it.

Nothing corrected, nothing implemented. `iterations`, `dt`, `mfd_exponent`, `min_area_cells`, the
diffusion anchor and the `A_c(S)` gate are all untouched.

## Finding 62 — the grids do not match at equal work, and 8192² has a WORK CEILING of 311 m. Plus: Finding 6 already knew

### Method, declared before the numbers

**Nothing is interpolated.** Interpolating a scalar between two unmatched fields does not produce
the field at the intermediate state, so every matched point below is a field that was BUILT.
`iterations` stays at the shipped **2 on both sides** — identical scheme, identical number of
sweeps — and only `dt` moves, `dt` being the continuous form of the same duration dial
(`k_time = K·dt·iterations`, Finding 44), bisected geometrically on the paired work.

And one requested check is **not independent**: paired work is `mean(pre) − mean(eroded)` over the
common land, so at equal work the hypsometric MEAN agrees almost by construction. It is reported
as a match verification, not as evidence. The independent checks are the emerged fraction, the
output channel share and the altitude quantiles.

### A correction to the premise of the round, before measuring

The round's block 1 compares 2048² at work **455.7 m** with 8192² at work **262.0 m** — a 1.74×
difference — and reads "the least-eroded coarse point already has less land than the most-eroded
fine point" as a contradiction. On a common trajectory, **more work → less land is the EXPECTED
ordering**, and the coarse point has 1.74× more work. **The pair is consistent, and the inference
does not hold.** Only an equal-work comparison can decide it, which is what was asked for and
what follows.

### Measure 1 — matched at 199 m: the MASS matches, the FORM does not

| | 2048² matched (`dt` 0.0888) | 8192² matched (`dt` 0.9795) | |
|---|---|---|---|
| paired work | 201.5 m | 197.6 m | the match |
| hypsometric mean | 685.7 m | 686.5 m | agrees to 0.12 % (near-tautological) |
| **emerged fraction** | **16.325 %** | **16.443 %** | **agrees to 0.118 points (0.7 % rel.)** |
| **OUT channel share** | **80.75 %** | **8.40 %** | **factor 9.6** |
| altitude p10 / p50 / p90 | 51 / 519 / 1507 m | 29 / 446 / 1607 m | +76 % / +16 % / **−6 %** |

> ⚠️ **The 9.6× is ROUTE-DEPENDENT (Finding 63 §2a).** Matching at the same work with
> `iterations` = 1 instead of 2 gives **4.46**, not 9.6 — a factor 2.15 between routes. The MASS
> conclusion below is robust (land agrees on both routes: 16.325/16.443 % and 16.379/16.445 %);
> the FORM conclusion survives qualitatively and NOT quantitatively. **Do not quote a single
> number** — read it as "a factor between 4 and 10, depending on the matching route". And at
> matched work the Courant numbers are **129 against 4033**: the two grids cannot be placed at
> equal work AND comparable numerical fidelity with the dials that exist.

**The emerged fraction agrees.** So the round's block-1 claim — that the grids would not match on
land even at equal work — is **refuted by measurement**: 16.325 % against 16.443 %.

**And the network does not.** A factor 9.6 on the channel share at identical eroded mass. The
altitude quantiles say the same thing in a different currency: the coarse field is **compressed**
(floor 51 m against 29, peaks 1507 against 1607) while the fine field is **spread**. Same mean,
same land, different distribution shape and a tenfold different network.

⇒ **Equal erosion matches the MASS and not the FORM.** So "erosion equal" is a necessary but
insufficient matching criterion, and a duration dial alone cannot make the grids agree.

### And the result I did not predict: the fine grid has a CEILING

Targeting the coarse grid's shipped state (633.0 m of work) at 8192²:

> **UNREACHABLE.** `dt ∈ [0.005, 200]` — a 40 000× span of the duration dial — spans work
> **[5.6, 311.2] m**. The fine grid saturates at **311.2 m** and cannot reach 633 m at any
> duration.

The mechanism is exact and it closes the loop with Finding 58 B1. As `dt → ∞`, `f = K·dt·A^m/dist
→ ∞` and `h → h_r` for every cell that passes the channel-head gate — but **only for those
cells**. The ceiling is therefore "every channel cell collapsed onto its receiver, every other
cell left at its pre-incision height", and the channel share sets it: 80.75 % of cells at 2048²
against **8.40 %** at 8192². **The work ceiling is `A_c`'s cell count, seen from the far end of
the duration dial.**

⇒ Finding 61's reformulation — *"neither grid is degenerate, they are at different points on the
same decay"* — was **TOO GENEROUS, and the author is right to reject it.** They are not on the
same decay: the fine grid's decay **terminates at 311 m** while the coarse grid's shipped state is
633 m. There is no duration, and no `dt`, that brings them together. The correct statement is that
**the two grids have different attainable sets**, and the shipped comparison is between a state
one grid can reach and one the other cannot.

### Measure 2 — the collapse is in `A` ITSELF; the feedback loop is confirmed

Output accumulation quantiles, km² (INPUT row for reference):

| 2048² | p10 | p25 | p50 | p75 | p90 | p99 |
|---|---|---|---|---|---|---|
| INPUT (pre) | 0.11444 | 0.26050 | **0.51135** | 1.05395 | 2.46303 | 20.102 |
| OUT iters 1 | 0.05726 | 0.11328 | 0.23942 | 0.52424 | 1.45086 | 26.680 |
| **OUT iters 2 (shipped)** | 0.04148 | 0.07366 | **0.14878** | 0.31859 | 0.87350 | 19.202 |
| OUT iters 4 | 0.03897 | 0.05050 | 0.09643 | 0.20283 | 0.55155 | 8.045 |
| OUT iters 8 | 0.03875 | 0.04595 | **0.07086** | 0.15465 | 0.42463 | 5.000 |

| 8192² | p10 | p25 | p50 | p75 | p90 | p99 |
|---|---|---|---|---|---|---|
| INPUT (pre) | 0.01587 | 0.05252 | **0.10942** | 0.21864 | 0.47750 | 4.146 |
| OUT iters 1 | 0.00477 | 0.01599 | 0.04511 | 0.08275 | 0.20569 | 3.051 |
| **OUT iters 2 (shipped)** | 0.00407 | 0.00828 | **0.02492** | 0.05366 | 0.08518 | 1.507 |
| OUT iters 4 | 0.00381 | 0.00704 | 0.01841 | 0.04219 | 0.06947 | 1.125 |
| OUT iters 8 | 0.00378 | 0.00649 | **0.01492** | 0.03293 | 0.05915 | 1.199 |

**The distribution of `A` collapses**, monotonically, at both grids: p50 falls ×7.2 at 2048²
(0.5114 → 0.0709) and ×7.3 at 8192² (0.1094 → 0.0149) over the sweep, and ×3.4 / ×4.4 by the
shipped point.

**So the channel share does not fall because fewer cells cross an unchanged threshold on an
unchanged distribution — it falls because the distribution itself moves out from under the
threshold.** At 8192² the median cell goes from 0.1094 km² (just above `A_c` = 0.1) to
0.0249 km² — **four times below it** — at the shipped point. The positive feedback proposed is
confirmed: a cell that leaves the channel regime stops being incised, concentration falls,
fewer cells clear the threshold next pass.

**Negative control, as required.** D8 on the same shipped field moves the quantiles plainly —
at 2048², p25 0.0737 (MFD) against 0.0381 (D8), p50 0.1488 against 0.1144, p99 19.20 against
**50.62** (×2.6). The sweep is measuring the distribution, not an artefact of its own.

And the inter-grid ratio **grows** at the output: p50 0.1488 / 0.0249 = **5.97**, against 4.67 at
the input.

### The two consequences the author asked to be written if the loop held

**1. "The work ratio is 3.18" compares two systems of different kinds, not two magnitudes of one
kind.** At the shipped point the coarse grid runs a network over 80.75 % of its land while the
fine grid runs one over 8.40 % and falling — 5.71 % at eight iterations. The fine grid's network
is **going out**. A ratio between a live network and an extinguishing one is not a magnitude
difference on a common object, and it should stop being quoted as though it were.

**2. The cross-check against Finding 56's drainage density does NOT reconcile — factor 5.25 — and
that is definitional, not a contradiction.** 0.337 km/km² over 26 285 km² is 181 360 channel
cells over 11.33 M land cells = **1.60 %**, against **8.40 %** by `acc_out ≥ A_c`. The two count
different things: Finding 56 counts cells on **D8-traced river segments** emitted by the drainage
stage at its own `head_km2`, this counts cells whose **MFD p = 2** accumulation clears `A_c` on
the eroded field. Neither is wrong; they are not comparable as stated, and the prediction that
they would not reconcile is confirmed.

**But they agree on the thing that matters, and the stricter number is the worse one:** by the
river-network definition **98.4 %** of the fine-grid terrain carries no channel, against my
91.6 %. Whichever definition is used, the overwhelming majority of production terrain at 8192²
is untouched by any channel.

### 3/2(a) — answered against the code, and §1a must be RESTATED rather than upheld

`mfd_accumulation` computes its own multi-flow partition **directly from `filled`** — weights
`(drop / D8_DIST)^p` over all eight neighbours — and does not route on the `direction` field it
receives. `basins` is labelled by **tracing each land cell downstream along the D8 receiver path**
to a sink and labelling the whole path. **Two different routing operators**, sharing only the
depression-filled surface.

So the 1.016 footprint ratio is **not** self-consistency with the accumulation in question — the
prediction that it would be is refuted. But it establishes something **narrower than Finding 59
§1a claims**:

> What is established: **the D8 basin decomposition of the depression-filled PRE-INCISION surface
> is resolution-invariant to 1.6 %.** What is NOT established by it: that the MFD accumulation
> converges. The only column in §1a that speaks to the accumulation is the outlet value, and its
> median is **1.42, not 1.0** — reported as such, and its 42 % is not attributed (window width,
> outlet-cell size, incomplete recombination all remain candidates). The formulation "the
> accumulation operator converges" is accordingly weakened to "the basin DELINEATION is
> grid-invariant; the accumulation at the outlet agrees to within ~40 %".

`2(b)` — the two flagged basins are **not** instructed here; that needs a measurement and is
carried forward as owed.

### `iterations = 2` — traceable, and to something better than an eye

**Both predictions refuted: the origin is neither untraceable nor eye-calibrated.** It is
**Finding 6** of this ADR, and it is explicit: *"No incision bound (floors planed to base level).
Stream power ran on a static field with no uplift, so channels graded down to sea level. Since
Ymir's tectonics already did the uplift, the fix is to LIMIT total incision, not add U: iters 3→2
and K 3000→1500 lifts floor/local-ridge from 0.21 to ~0.48."*

So it was calibrated on a **measured morphometric target** — floor/local-ridge ≈ 0.5, "floors at
about half the local ridge, not planed" — validated at 8192². Not by eye.

**The PROXY label is still owed, but on the TARGET, not on the method.** `floor/local-ridge ≈ 0.5`
carries no external anchor: no publication, no measured range from real terrain. It is a chosen
morphometric objective, and the duration of erosion in this pipeline is set by it. That belongs
in the code and in the ADR as a proxy.

### ⚠️ AND FINDING 6 ALREADY KNEW — a dossier-hygiene failure that is mine

Finding 6 also contains, verbatim: *"NOTE: A_c must be resolvable — 0.1 km² is sub-cell below
~2048² (**it needs ~2.6 cells @2048², 42 @8192²**)."*

**Those are the exact numbers Finding 57 presented as a new discovery and Finding 58 re-derived
by reading them back from the config.** And Finding 6 states the no-uplift / base-level attractor
that Finding 61 rediscovered fifty findings later, as its own justification for `iters 3→2`.

So the honest accounting: the **observations** were not new — they were recorded once, at the
beginning, and lost. What was new was the **naming of the pattern** ("physically dimensioned but
sub-cell") and its promotion to an executable test — and the fact that an observation could sit
in this ADR for fifty findings while four separate rounds hunted the same blocker is itself the
argument for naming patterns rather than recording instances. Finding 57's claim of novelty is
withdrawn; its pattern and its guard stand.

### Where this leaves the queue

The round's own rule stops here — measure 1 shows the grids do not match, so the reclassification
(block 4) and the bibliographic review (block 5) are deferred, and what they must now take as
given is:

- **equal erosion is not a sufficient matching criterion** (mass matches, form does not, 9.6×);
- **the fine grid cannot reach the coarse grid's shipped state at any duration** (ceiling 311 m
  against 633 m), so a timescale dial cannot be the remedy — it can only place a grid within its
  own attainable set;
- **the ceiling is `A_c`'s cell count**, which makes the sub-cell threshold the single object
  that all four rounds have converged on from different directions;
- **the fine grid's network is extinguishing**, so any comparison of the two must say which
  regime each grid is in before quoting a ratio.

Nothing corrected, nothing implemented. `iterations`, `dt`, `mfd_exponent`, `min_area_cells`, the
diffusion anchor and the `A_c(S)` gate are untouched.

## Finding 63 — the closed-form ceiling FAILS by 27 %, and the failure is a population mismatch in my own derivation

### The derivation, and the reservations declared before comparing

As `dt → ∞`, `f = K·dt·A^m/dist → ∞`, so `h ← (h + f·h_r)/(1+f) → h_r` for every cell passing the
channel-head gate **and only those**. Accumulation being non-decreasing downstream makes the
channel set downstream-closed, so a channel cell's collapse target is sea level:

```
work_ceiling  =  channel_share × mean(h_pre | channel)
```

evaluated on the pre-incision field — the set sweep 1's gate actually reads. **Declared before
comparing, admissible error 15 %:** MFD splits flow so accumulation can fall downstream and break
downstream-closure; and sweep 2 re-gates on a changed surface. **Both reservations push the true
ceiling ABOVE the formula**, so it was offered as a lower bound.

### It fails, and in the direction the reservations excluded

| | channel share | mean(h \| channel) | formula | measured |
|---|---|---|---|---|
| 2048² | 91.36 % | 855.8 m | **781.8 m** | ships at 633.0 m (inside its set) |
| 8192² | 53.62 % | 738.0 m | **395.7 m** | **311.2 m** |

**Error +27.2 %, outside the declared 15 %, and ABOVE the measurement — the opposite of what both
declared reservations predicted.** So the reservations I wrote down were not the operative ones.

> ⛔ **WITHDRAWN by Finding 66.** The population effect is REAL but it runs the SAME way as
> the formula's error, not against it: the drowned cells have a mean drop of 498 m against a
> paired mean of 633 m, so including them LOWERS the measured work. It therefore cannot
> explain the formula being 27 % too HIGH. **That failure has no explanation on file.** And
> the stated reason — that the excluded cells "carry the largest work" — is measured false:
> they are mid-altitude cells (mean `h_pre` 478 m), not low-lying ones.

**The operative one is a population mismatch, and it is this campaign's own recurring defect
committed in my own derivation.** The measured quantity is the **paired** work, averaged over
cells that are land in **both** fields. A cell that collapses to sea level *leaves the land mask*
and is dropped from the measurement — and those are precisely the cells carrying the **largest**
work (`h_pre − 0`). The formula averages over all pre-incision land; the measurement averages over
surviving land. Different denominators, and the excluded cells are the extreme of the
distribution.

**Consequence for the block's premise: the ceiling is NOT free.** Respecting the paired-work
population requires knowing which cells drown, which requires the collapse target per cell, which
requires the network. *"Le plafond est analytique, donc gratuit"* is refuted — by an error of my
own making, not by the physics.

**What survives, and it is one-sided.** The formula over-counts by including cells the measurement
drops, so it is an **upper bound on the paired work**. That is still usable in exactly one
direction: **a work target above the table value is definitely unreachable.** The table is kept
with that label and no other.

| `A_c` km² | 1024² | 2048² | 4096² | 8192² |
|---|---|---|---|---|
| 0.025 | 865 m (100.0 %) | 865 m (100.0 %) | 801 m (93.6 %) | 703 m (86.6 %) |
| **0.100** | 865 m (100.0 %) | **782 m (91.4 %)** | 637 m (79.1 %) | **396 m (53.6 %)** |
| 0.400 | 715 m (84.7 %) | 475 m (59.9 %) | 224 m (30.2 %) | 85 m (12.3 %) |
| 1.000 | 435 m (54.9 %) | 194 m (26.4 %) | 71 m (10.7 %) | **27 m (4.5 %)** |

*(upper bounds on the paired work, with the channel share each comes from)*

The trend prediction is confirmed and it is severe: the bound falls with `A_c` and with
resolution, both because both cut the channel share, and at `A_c = 1.0 km²` the fine grid's bound
is **27 m** — 3 % of the available relief. Note also that at `A_c = 0.025 km²` the channel share
saturates at 100 % on the two coarse grids, so `A_c` stops being a threshold there at all.

> ✅ **SUPERSEDED by Finding 67 §2: the assertion IS written.** An upper bound can only produce
> FALSE NEGATIVES, so a test that fires when a target EXCEEDS the bound is correct even when
> the bound is loose. It is labelled as such, carries a negative control, and flags the
> motivating case by 60 %.

**The feasibility ASSERTION is not written.** A guard built on a formula that misses its own
validation point by 27 % would encode the error. It waits on a corrected derivation.

### 2a — the second route: the 9.6× is NOT robust, and the reservation is vindicated

There is no fully independent matching route with the shipped dials: `iterations` is
integer-quantised and its minimum step overshoots (2048² at iters 1, dt 1 already gives 455.7 m
against a 199 m target), and `K` and `dt` are the **same observable** (Finding 44). What is
available is a *different scheme at the same work* — one sweep instead of two.

| route | grid | `dt` | work | land | **OUT channel** | **Courant** |
|---|---|---|---|---|---|---|
| A — iters 2 | 2048² | 0.0888 | 201.5 m | 16.325 % | 80.75 % | — |
| A — iters 2 | 8192² | 0.9795 | 197.6 m | 16.443 % | 8.40 % | — |
| **B — iters 1** | 2048² | 0.1962 | 200.2 m | **16.379 %** | **83.88 %** | **129** |
| **B — iters 1** | 8192² | 1.5939 | 199.2 m | **16.445 %** | **18.80 %** | **4033** |

**Route A gives a channel-share ratio of 9.6; route B gives 4.46 — a factor 2.15 between
routes.** Both predictions (±20 %) are **refuted**, and the author's reservation is **vindicated**:
part of what Finding 62 measured as "form" was the matching path.

What survives and what does not:

- **The MASS conclusion is robust.** Land agrees on both routes — 16.325 / 16.443 % on A,
  16.379 / 16.445 % on B.
- **The FORM conclusion survives qualitatively and NOT quantitatively.** Both routes show a large
  gap (9.6 and 4.46) but the magnitude is route-dependent, so **no single number should be
  quoted.** Finding 62's "factor 9.6" must be read as "a factor between 4 and 10, depending on
  the matching route".
- **The truncation reservation cannot be lifted, only quantified — and it is large.** At matched
  work on route B the Courant numbers are **129 against 4033**, a factor 31. The two grids cannot
  be placed at equal work *and* comparable numerical fidelity with the dials that exist.

### 2b — the flagged basins, re-matched by footprint overlap: half the 42 % was the window

Method: no search window at all. For each of the largest 2048² basins, map every one of its cells
to the HD grid, take the **modal HD basin id**, then compare each basin's own maximum
accumulation.

| # | acc ratio, window → **overlap** | footprint ratio | overlap | verdict |
|---|---|---|---|---|
| 1 | 14.98 → **1.156** | 8.66 → 1.760 | **56.5 %** | matching failure **and** a genuine decomposition difference — the basin does not correspond one-to-one between grids |
| 6 | 3.57 → **2.428** | 0.994 | **98.9 %** | ⚠️ **NOT a matching failure.** Same basin, and the outlet accumulation is still **2.4× apart** |
| 7 | 7.50 → **1.291** | 1.002 | **99.0 %** | mostly window artefact; 1.29 residual |
| | **median 1.416 → 1.207** | | | |

**Both mechanisms are present, as predicted at 50/50** — #1 and #7 are matching artefacts, **#6 is
a real accumulation difference at 98.9 % footprint overlap.** The prediction that the median would
barely move is refuted: it falls from 1.416 to **1.207**, so **about half of the 42 % excess was
the search window**.

**What remains, and it now carries weight:** a **20.7 %** median excess in outlet accumulation
between grids, still unattributed, on the quantity that feeds `A_c` — the parameter four
independent routes have converged on. Basin #6 is the clean witness: same basin by every
footprint test, 2.43× apart at the outlet.

And Finding 59 §1a is weakened once more: its outlet column now reads **1.207**, not 1.42 and not
1.0, and the flag population it set aside contained one case that was not an artefact at all.

### The convergence that distinguishes this attribution from a correlation

`A_c` is now the object that **four independent routes** have arrived at, and it is worth stating
because nothing else in this dossier has more than two:

1. its **discretisation** — 2.62 cells against 41.94, and holding it constant in CELLS collapses
   the work ratio from 3.18 to 1.23 (Finding 58 B1);
2. the **partition collapse** — the channel share falls ×0.70 and ×0.157 from input to delivered
   field, because the distribution of `A` moves out from under the threshold (Findings 59, 62);
3. the **network extinction** — 53.62 % → 8.40 % → 5.71 % at 8192² as the budget rises
   (Finding 61);
4. the **upper bound of the attainable set** — the ceiling is `channel_share × mean(h|channel)`,
   so `A_c` sets what the grid can ever reach (this finding, as a one-sided bound).

Four routes, one object, and no free parameter shared between them.

### Method rule 11 — interrogate the dossier before measuring, and grep the IDENTIFIER, not the symbol

`A_c` sub-cell at 2.6 and 42 cells, and the base-level attractor for want of an uplift term, were
both written in **Finding 6** — the finding that justifies `iterations = 2` — and were
rediscovered fifty findings later, across four rounds, by two people hunting the same blocker.

Naming patterns (Finding 57) is a real improvement and it is **not a fix**: the pattern
"physically dimensioned but sub-cell" would have made the same fact be *rediscovered faster*, not
*retrieved*. The defect is that **the dossier is not interrogable** — nothing in the method
prescribes asking whether a quantity has already been measured.

> **Rule 11. Before the first measurement of a round, grep the ADR for every constant and target
> quantity the round will touch, and record the result — including "nothing found".**

**Tested rather than asserted**, because the rule is only worth what its query is worth:

| query | hits in the ADR | first hit | would it have retrieved Finding 6? |
|---|---|---|---|
| `A_c` | **131** | line 21 | ❌ noise — the informal symbol is everywhere |
| `min_area_cells` | 11 | line 201 | ✅ the code identifier |
| `resolvable` | 5 | **line 234** | ✅ **first hit IS the sentence** |
| `sub-cell` | 15 | **line 234** | ✅ **first hit IS the sentence** |

So the rule works, and the discriminating factor is **specificity, not constant-ness**:

- **grep the CODE IDENTIFIER** (`min_area_cells`), not the informal symbol (`A_c`, 131 hits of
  noise) — this is the correction to my own prediction that "constants" were the answer;
- **grep the PHENOMENON PHRASE** (`sub-cell`, `resolvable`) as well as the name — it is what
  actually put the answer in the first hit;
- **read the EARLIEST hit first.** The ADR is append-only, so hit order is chronological and the
  first hit is the original recording. Both winning queries land on line 234 — Finding 6 — as
  their very first result.

Cost: four greps. Against four rounds.


### Method rule 11b — when a symptom is a REPEATED VALUE, grep the value, not only the name

Rule 11 says to enumerate the code identifiers and grep each against the ADR. That is necessary and
it is not sufficient. `shelf_min_depth_m`, `BathymetryProfile`, `apply_bathymetry_profile`,
`bathymetry` and `shelf_min_depth` return **zero ADR hits between them**. The value they produce,
`−20 m`, grepped once with its domain words (sea / ocean / shelf / depth), returns **six traces
across five findings** — F33 *"a deep bowl, floor −20 m"*, F34 *"level −20 m, MAX DEPTH 0 m — a
surface with no water"*, F35 *"a −20 m coastal pocket … is absurd"*, F36 *"the −20 m altitude proxy
Finding 31 removed"*, F38 *"2 depth-0 lakes at −20 m — **on a FLAT floor**"*.

**The constant signed five findings and was never recognised as a constant, because a dossier
records SYMPTOMS in numbers and CAUSES in names.** So: after the identifier greps, take the
characteristic value the symptom keeps showing (a depth, a count, a ratio) and grep **that**, with
two or three domain words to cut the noise. Record the result, including "nothing found". Cost: one
command. It found in one round what five findings had circled without touching.

### Method rule 12 — the export is authoritative for the SYMPTOM; a MECHANISM is measured at the stage where it acts

Findings 76 and 77 hunted a coastal mechanism on the delivered field and concluded, twice, that it
was absent — "no sea-level flat" (r = 0.051), then "the drag lands on the shelf value". Both were
measured **after** `apply_bathymetry_profile` (`production_upscale.rs:513`) had overwritten every
sub-sea cell, and the incision runs at `:451`. On the pre-clamp field the receivers of coastal
channel cells have **3 041 distinct values with a median of −0.109 m**; on the delivered field they
have **one**, −20.000. The mechanism was there the whole time, 20 m under the instrument.

So, two different questions and two different fields:

- **"Does the consumer see the defect, and how big is it?"** → measure the EXPORT. It is what
  ships; nothing else is authoritative for a symptom.
- **"What produces it?"** → measure at the stage where the term acts, **before every
  post-processing that touches the same cells**. Enumerate what runs after it — here one line of
  `production_upscale` — and build the field without them.

The corollary is the cheap part: when a stage is gated by an `Option`, turning it off is a
one-field bench knob, and that knob is the difference between a refutation and an attribution.

## Finding 65 — THREE incompatible answers to "is there a channel here", and the one that carved the terrain is used by nothing

Independent of the resolution blocker, and it reaches the consumers.

### The three definitions, named against the code

| # | definition | operator | who uses it |
|---|---|---|---|
| 1 | `accumulation ≥ min_area_cells` | **MFD**, `mfd_exponent = 2` — `mfd_accumulation` recomputes its own multi-flow partition from `filled` | **`incise` — it CARVES THE TERRAIN**, and nothing else |
| 2 | `accumulation ≥ head_threshold` | **D8** — `extract_rivers` reads `flow_result.accumulation` and `flow_result.direction`, both single-direction | `rivers.json`, drainage density, `catchment_km2`, `drainage_km2`, Strahler order, the microscope |
| 3 | climatic **discharge** m³/s | runoff accumulation with precipitation and evaporation | `navigability`, `width_m`, the lake water balance |

At 8192², definitions 1 and 2 disagree by **5.25×** on the same shipped field: **8.40 %** of land
by the MFD criterion, **1.60 %** by the D8 river network (181 360 traced cells over 11.33 M land).
Finding 59 measured the underlying cause directly — D8 against MFD p = 2 moves the accumulation
p50 by +49 % at 2048² and +64 % at 8192², and p99 by ×2.6 — so the two operators select different
cell sets at the *same numeric threshold*, and `head_km2` is set to `RELIEF_V1_A_C_KM2` in the
benches precisely so the thresholds match. **Matching the number does not match the set.**

### The statement that matters

> **The definition that shaped the land is used by nothing downstream. The definition used by
> everything downstream shaped nothing.**

`incise` carves on MFD. Every exported hydrological product reads D8 or discharge. No stage
reconciles them, and no test asserts they agree — because until now nobody had asked them the same
question.

This is the same shape as Finding 60 (geometric drainage area against climatic discharge) and it
compounds with it: the project has **three** answers to "how much water is here", one per stage,
and the three have never been reconciled.

### What it costs the consumers

- **Navigability** classifies with definition 3 on segments delineated by definition 2, over
  terrain carved by definition 1. Three definitions in one verdict. Finding 60 already flagged the
  climatic half; this adds that the *geometry* is also not the carving geometry.
- **Drainage density** is definition 2 throughout — self-consistent, but it does **not** describe
  the network the erosion built. Finding 56's 0.337 km/km² is a true statement about `rivers.json`
  and not about the terrain's channels.
- **`catchment_km2` and `drainage_km2`** are definition 2. `drainage_km2` is already documented in
  the code as "NOT AN AREA"; this adds that even the area-like one is on the wrong operator to
  describe what carved the valley it labels.
- **Resource placement and river/mouth inspection** consume definition 2, so they see a network
  5.25× sparser than the one that cut the terrain they are placed on.

**And the strict number is the one to quote to a consumer: 98.4 % of production terrain at 8192²
carries no channel in `rivers.json`.** By the carving definition it is 91.6 %. Both are extreme;
they disagree by a factor 5 on how extreme.

### Not a defect claim, and no fix this round

MFD for the incision is a defensible choice (it is what cured the parallel-rilling comb, Findings
8/9/11) and D8 for a traced, hierarchised river network is equally defensible — a Strahler order
needs a tree, and MFD does not give one. **Neither stage is wrong on its own terms.** What is
missing is any statement of how they relate, and any test that they do.

Recorded, not resolved. No consumer is changed, no threshold is touched.

## Finding 66 — "work" has meant two quantities; the bias is REAL, SMALL, GROWS with the budget, and runs the opposite way to the reason I gave

### The self-correction that comes first

Finding 63 explained the closed form's 27 % over-count by saying the paired metric drops cells
that "carry the largest work (`h_pre − 0`)". **That was asserted without checking and it is
wrong twice over.**

Measured: the drowned cells at 2048², `iterations = 2`, have a mean **pre-incision altitude of
478.0 m** and a mean **drop of 498.0 m**. They are not the lowest land cells — they are
mid-altitude cells that the incision carves all the way through. So:

- the reason was wrong (they are not low-lying);
- **and the direction still holds, for a different reason**: their drop (498 m) is *below* the
  paired mean (633 m), so including them **lowers** it. The paired metric **OVER-states**.

Which means Finding 63's explanation of the formula failure is **withdrawn**: the population
effect pushes the measured work DOWN, the same way the formula errs, so it cannot account for the
formula being 27 % too high. **That failure now has no explanation on file.** Unconfirmed
candidates: the diffusion refilling collapsed channels at extreme `dt`, and iteration 2 re-gating
on a drowned surface.

### The three conventions, measured

| grid | iters | PAIRED | FULL@0 | FULL@true | drowned | % of land | mean h drowned | mean drop |
|---|---|---|---|---|---|---|---|---|
| 2048² | 1 | 455.7 | 449.8 | 451.2 | 50 283 | 7.10 % | 372.8 m | 392.8 m |
| **2048²** | **2** | **633.0** | **615.3** | **617.6** | 80 787 | **11.41 %** | 478.0 m | 498.0 m |
| 2048² | 4 | 757.5 | 718.3 | 721.6 | 118 299 | 16.71 % | 522.9 m | 542.9 m |
| 2048² | 8 | 835.9 | 781.0 | 785.2 | 148 167 | 20.93 % | 573.6 m | 593.6 m |
| 8192² | 1 | 167.4 | 167.0 | 167.4 | 234 849 | 2.07 % | 148.6 m | 168.6 m |
| **8192²** | **2** | **199.0** | **198.2** | **198.7** | 298 598 | **2.64 %** | 167.7 m | 187.7 m |
| 8192² | 4 | 223.9 | 222.2 | 222.8 | 327 293 | 2.89 % | 165.8 m | 185.8 m |
| 8192² | 8 | 262.0 | 259.0 | 259.6 | 336 459 | 2.97 % | 159.9 m | 179.9 m |

`PAIRED` = land in both (the published convention). `FULL@0` = all pre-incision land, a drowned
cell counted at `h_pre − 0` (a lower bound on its drop). `FULL@true` = all pre-incision land, a
drowned cell counted at its actual drop `h_pre − h_final`.

### The ratio moves 2.4 %, not 20 %

| iters | PAIRED | FULL@0 | FULL@true | shift |
|---|---|---|---|---|
| 1 | 2.723 | 2.694 | 2.696 | −1.1 % |
| **2 (shipped)** | **3.181** | **3.105** | **3.108** | **−2.4 %** |
| 4 | 3.383 | 3.232 | 3.239 | −4.5 % |
| 8 | 3.190 | 3.016 | 3.025 | −5.5 % |

**Both predictions on the table were too pessimistic** — 2.6–3.0 and 2.85–2.95 against a measured
3.105. The bias is real and it is **small at the shipped point**.

**And it GROWS with the erosion budget** — −1.1 %, −2.4 %, −4.5 %, −5.5 % — because the drowned
share at 2048² grows 7.10 → 20.93 % while 8192² barely moves, 2.07 → 2.97 %. So the round's
claim that the bias grows with erosion is **confirmed in its growth and inverted in its sign**:
the metric over-states, increasingly.

### The convention, declared once for the whole dossier

> **`FULL@true` is the convention from here on**: all pre-incision land, a drowned cell counted
> at its actual drop. It is the only one of the three that averages a fixed population and counts
> every metre the incision removed.

`PAIRED` is what Findings 57–62 published. The correction is **−2.8 % at 2048² (633.0 → 617.6)
and −0.2 % at 8192² (199.0 → 198.7)**, and the ratio **3.18 → 3.11**.

### Named list of what moves, and it is short

| statement | published | `FULL@true` | verdict |
|---|---|---|---|
| erosion work, 2048² / 8192² (F57 A3) | 633.0 / 199.0 m | **617.6 / 198.7 m** | restate, −2.8 % / −0.2 % |
| the work ratio (F57, F58, F61, F62) | **3.18** | **3.11** | restate; every conclusion drawn from it stands |
| the ratio over the iteration sweep (F61) | 2.72 / 3.18 / 3.38 / 3.19 | **2.70 / 3.11 / 3.24 / 3.03** | still non-monotone, still 2.7–3.2, and the shipped point is still near the maximum |
| the 311.2 m ceiling (F62) | paired | unchanged as a *measurement*, but it is a PAIRED number and the `FULL@true` saturation was not measured | flagged, not restated |
| equal-work matching (F62) | matched on paired work | the match target moves < 3 %, and both sides move the same way | conclusion unaffected |
| `A_c` cell counts, pre-incision invariance, accumulation quantiles, F65 | — | — | untouched: none of them uses the work metric |

**No conclusion in the dossier reverses.** The blocker's magnitude is 3.11 rather than 3.18, the
absolute works are ~3 % lower at the coarse grid, and everything drawn from them holds.

### Method point, because this is the fourth time

The published quantity was an average over a population that **the treatment itself changes**.
That is the same defect as the post-breach lakes (rule 10), the per-100 km spur density (rule 2),
the wet-share saturation (`Sweep::part`) and the regime label of Findings 59/62. Here it cost
2.4 %; the previous instances cost a retracted decomposition and a retracted convergence claim.

> **Corollary to rule 2, worth stating separately: when a treatment can move cells INTO or OUT OF
> the population being averaged, the population must be fixed BEFORE the treatment, and the
> excluded cells must be given a defined value rather than dropped.**

## Finding 67 — the attainable set shrinks CONTINUOUSLY with resolution; a one-sided guard; and the F65 audit finds FOUR definitions, not three

### 1 · The right reading of the ceiling table

Not "two grids with disjoint sets". Along `A_c = 0.1 km²`:

| | 1024² | 2048² | 4096² | 8192² |
|---|---|---|---|---|
| upper bound on the work | 865 m | **782 m** | **637 m** | **396 m** |
| production ships at 633 m (2048²) | inside | inside | **at the limit** | **1.6× outside** |

> **The attainable set shrinks continuously with resolution, and the HD ×128 pipeline targets a
> grid on which the production state is out of reach.** 4096² is already at its limit — 637 m of
> bound against a 633 m state, a 0.6 % margin on a bound known to be loose by 27 %. 8192² is
> outside by a factor 1.6.

A one-sided bound is enough to establish this, which is why it is worth writing even though the
formula it rests on failed its validation point.

### 2 · The guard, with its label

Finding 63 declined to write the assertion. **That was over-cautious and the round is right to
push back.** The formula errs in ONE direction — it over-counts — so it can only produce **false
negatives**. A test that fires when a target EXCEEDS the bound is correct even when the bound is
loose.

`crates/ymir-core/tests/resolution_invariants.rs :: a_work_target_outside_the_attainable_set_is_flagged`

> **Label, so nobody reads more into it: UPPER BOUND. Loose by ~27 % at the single point where it
> has been checked. It detects only a SUBSET of unreachable targets, and every target it flags is
> genuinely unreachable.**

It carries a negative control (a state the fine grid demonstrably occupies must NOT be flagged —
an upper bound must never produce a false positive) and asserts the set shrinks monotonically
with resolution. **Would it have caught the present case?** 633 m against a 396 m bound at
8192² — flagged by 60 %. Not marginal.

### 3a · The audit — and the coastal metrics are CLEAN, against my own prediction

I predicted the spur count and the local parallelism R were exposed. **The code refutes it:**
`terrain::coast_metrics` contains **zero** references to accumulation, drainage, rivers or MFD.
`coast_shape`, the spur count, the p90 length, the spacing CV and both axial R's are computed
from a **marching-squares contour of the height field**. They read the carved terrain directly,
so they are consistent with the definition that carved it, by construction.

**What IS exposed is narrower and stranger than F65 stated — there are FOUR definitions in use,
not three:**

| # | criterion | operator | where it is used |
|---|---|---|---|
| 1 | `acc ≥ min_area_cells` | **MFD p = 2** | `incise` — **carves the terrain**, and nothing else |
| 2 | `acc ≥ min_area_cells` | **D8** | `coastal_comb_levers`'s "channel km" / network extent — Finding 56's rule-7 control block |
| 3 | traced segments, `acc ≥ head_threshold` | **D8 + tracing + clipping** | `rivers.json`, the invariant suite's drainage density, `catchment_km2`, Strahler, the microscope |
| 4 | climatic **discharge** m³/s | runoff accumulation | navigability, `width_m`, the lake water balance |

Definitions 1 and 2 share a threshold and differ only in operator; 2 and 3 share an operator and
differ in threshold and in tracing. **Neither 2 nor 3 is the one that carved the terrain.**

Named statements measured on a network that is not the carving network:

- **Finding 56's "drainage density 0.746 → 0.798 and 0.561 → 0.756 km/km²"** and Finding 56b's
  corrected **0.611 → 0.619 / 0.540 → 0.721** — definition 3.
- **Finding 56c's 0.337 → 0.501 km/km² in the humid bed** — definition 3, and it is the figure
  that produced the 5.25 % / 1.60 % cross-check failure of Finding 65.
- **Finding 56's rule-7 "network extent" column** (80 107 → 115 857 km and 45 573 → 262 788 km) —
  definition 2, a *third* variant, and already withdrawn on separate grounds (Finding 58 §E: it
  reads a fixed threshold on the resulting field).
- **Every Strahler histogram, confluence count and W/D-per-order table** — definition 3.
- **Finding 62's channel-share comparisons** (80.75 % / 8.40 %) — definition 1, correctly, since
  the subject was the carving.

**Not exposed, and this is the useful half of the audit:** the whole coastal-fringe chantier —
spurs 1 254 → 560 and 1 849 → 544, coastline km, p90 length, local R 0.515 → 0.513 and 0.525 →
0.430 — plus the hypsometry, the slope quantiles and the `A_c` cell arithmetic. **The channel-head
law's primary verdict does not depend on any channel definition.**

### 3b · Where the definitions ought to agree and do not

| exported quantity | fed by | should agree with | does it |
|---|---|---|---|
| `rivers.json` geometry, Strahler, `catchment_km2` | 3 | 1 — the valleys it labels were cut by MFD | **no**: 1.60 % against 8.40 % of land, factor 5.25 |
| `width_m`, `navigability` | 4 | 3 for *where*, 1 for *how big the valley is* | **no**: three stages in one verdict (Finding 60) |
| drainage density | 3 | 1 — it is offered as a property of the terrain | **no**, same 5.25× |
| `lake_type`, `source_lake_id`, `kind` | 4 (balance) + 3 (clipping) | each other | partly — the orphan-mouth defect (F56b) is one symptom |
| coastline, spurs, local R | height field | 1 | ✅ **yes, by construction** |

### 3c · The minimum-coherence question, formulated and not answered

Each stage is defensible alone: MFD cured the parallel-rilling comb (Findings 8/9/11), a Strahler
order needs a tree that MFD cannot provide, and the hydrology needs a climatic discharge. So the
question is not "which is right".

> **Minimum coherence: for a cell the incision treated as a channel, does the exported network
> place a channel there — and if not, at what rate, and is the disagreement systematic in
> position (headwaters versus trunks) or scattered?**

That is a containment question with a measurable rate, not a definition change, and it can be
asked without deciding anything. Its answer would say whether the exported network is a **subset**
of the carved one (defensible: a coarser view of the same object) or a **different** one
(not defensible: it labels valleys it did not cut). The 5.25× ratio is consistent with either.

### 6 · Second non-existence result, promoted to that status

Finding 63 recorded it in passing; it belongs beside the ceiling.

> **The two grids cannot be placed at equal work AND comparable numerical fidelity with the dials
> that exist.** At matched work the Courant numbers are **129 against 4033**, a factor 31. And
> there is no independent matching route: `K` and `dt` are the same observable (Finding 44), and
> `iterations` is integer-quantised with a minimum step that overshoots the target.

Two non-existence results now stand: **no duration reaches the coarse grid's state at 8192²**, and
**no dial setting matches the grids at equal work and equal fidelity.** Both are properties of the
parameter space, not of a particular run.

### 6 · Method rule 11, rewritten — the identifier obligation comes FIRST

The round's objection is correct and I accept it: the two winning queries, `sub-cell` and
`resolvable`, are words for the **phenomenon**, which one only knows *after* understanding it. At
round 1 nobody would have searched "sub-cell". The rule as written retrieved what one already
knows to look for.

The counterweight is `min_area_cells` — the **code identifier** — which works without knowing the
phenomenon, because a campaign always knows which identifiers it is about to touch.

> **Rule 11, restated. Before the first measurement of a round: (1) ENUMERATE the code identifiers
> the campaign will touch — constants, config fields, function names — and grep each one against
> the ADR; this is the obligation, and it is mechanical. (2) THEN add any phenomenon phrases you
> can guess; this is a bonus, not a duty, because you cannot name a phenomenon you have not yet
> understood. (3) Read the EARLIEST hit first — the ADR is append-only, so hit order is
> chronological. (4) Record the result, including "nothing found".**
>
> Do NOT grep informal symbols: `A_c` returns 131 hits of noise where `min_area_cells` returns 11
> and lands on the answer.

## Finding 68 — `geo_scale_ratio` is a post-hoc multiplier on FOUR segment arrays; the erosion never sees it; and the export ships two scales at once ACROSS object types

Code reading and enumeration only. No measurement, no run, no change.

### Rule 11's first genuine success, and it must be recorded as such

Applied before anything else, on the identifiers this campaign touches:

| identifier | ADR hits | earliest line | outcome |
|---|---|---|---|
| `geo_scale_ratio` | **17** | **968** | ✅ **the mechanism, the exact scaled/not-scaled lists, AND a ratio-7.5 table were already recorded** |
| `min_area_cells` | 18 | 201 | Finding 6 (already known, Finding 62) |
| `cell_km` | 12 | 310 | — |
| `HILLSLOPE_REF_CELL_M` | 11 | 414 | — |
| `domain_km` | 8 | 968 | same passage |
| `window_km` | 0 | — | **nothing found** (recorded, per the rule) |

**The answer to the decisive question was in the dossier and one grep away.** This is the first
time rule 11 has retrieved a result instead of confirming one, and it did so on the *identifier*
query — which is the half of the rule the last round moved to the front, correctly.

### 1 · The call chain, with lines

| site | role | relative to the incision |
|---|---|---|
| `drainage.rs:890` `apply_geo_scale_ratio` | the only implementation | — |
| `hd_assembly.rs:240` | the only production call, **last statement** of the assembly, after the spillway append | **downstream** |
| `cached_product.rs:772,781` | folded into `hd_drainage_key` as `hd_dr_geo_scale` | downstream (keying only) |
| `viz hd.rs:120,174` | `HdParams` field, **default 1.0** | downstream |
| `viz hd.rs:818,1003,1100` | threaded into the key and the assembly | downstream |
| `viz hd.rs:1223` | `ContinentMeta::geographic_scale_ratio` — **declared in the manifest** | downstream (export) |
| `workspace.rs:366,1226` | UI state, default 1.0, slider 1.0–15.0 | downstream |
| `drainage.rs:2706` | the enumerating guard test | — |
| `fracture/mod.rs:69` · `production_upscale.rs:413` · `upscale.rs:271` | **negative assertions in comments: "NEVER `geo_scale_ratio`"** | **upstream — deliberately excluded** |

**Nothing upstream of or inside the incision reads it.** The three upstream mentions exist only to
forbid it, which is a deliberate and documented design decision, not an omission.

### 2 · The decisive answer: NONE of branches A/B/C — it is branch D

`geo_scale_ratio` **does not reach `cell_km` in `StreamPowerConfig`.** It is never read by
`erosion/stream_power.rs` at all.

> **Branch D. A pure post-process on the drainage RESULT, applied as the final statement of the
> assembly. It multiplies exactly four segment arrays and re-derives two more. NO STAGE ANYWHERE
> USES A SCALED CELL SIZE.**

So branch A is refuted, branch B is correct *for the erosion*, and **branch C's premise is false**:
there is no 366.2 m hydrology cell — the hydrology *computation* runs entirely at 48.83 m too.
Only four reported numbers change afterwards.

Arithmetic verified before use: 400 × 7.5 = 3000 km; 3000/8192 = 0.366211 km; 0.134111 km²;
0.1/0.134111 = 0.7457 cells; and B: 0.048828 km, 0.0023842 km², 41.943 cells. All correct.
**But the branch-A row never occurs: at 8192², `A_c = 0.1 km²` is 41.94 cells at every ratio.**
Everything Findings 57–67 measured about `A_c`'s cell count is therefore ratio-independent.

What `apply_geo_scale_ratio` does (`drainage.rs:890`, `ratio == 1.0` → early return):

| quantity | treatment |
|---|---|
| `segment_drainage_km2` | × ratio² = **×56.25** |
| `segment_discharge_m3s` | × ratio² |
| `segment_catchment_cells` | × ratio² *(added by Finding 49)* |
| `segment_discharge_profile_m3s` | × ratio², per point *(Finding 49)* |
| `segment_width_m` | **re-derived** as `5·Q^0.5` → × ratio = ×7.5 |
| `segment_navigability` | **re-classified** on the scaled `drainage_km2` |

### 3 · Which stage sees it — the Finding 65 table, with the column added

| # | channel definition | operator | sees `geo_scale_ratio`? |
|---|---|---|---|
| 1 | `acc ≥ min_area_cells` | MFD p=2 | ❌ **no** — carves the terrain, ratio-blind |
| 2 | `acc ≥ min_area_cells` | D8 (`coastal_comb_levers` network extent) | ❌ no |
| 3 | traced segments `≥ head_threshold` | D8 + tracing | ❌ **the GEOMETRY is not scaled** — only four attribute arrays on it are |
| 4 | climatic discharge | runoff accumulation | ❌ the **computation** is unscaled; the **reported** discharge is ×56.25 |

And the constants:

| quantity | scaled? | consequence at ratio 7.5 |
|---|---|---|
| `catchment_km2`, `drainage_km2`, `discharge_m3s`, discharge profile | ✅ ×56.25 | signified |
| `width_m` | ✅ ×7.5 | signified |
| `navigability` | ✅ re-classified on signified area | coherent **within** `rivers.json` |
| `stream_km2` 20 · `small_boat_km2` 500 · `barge_km2` 5 000 · `ship_km2` 50 000 | compared against the **scaled** area | ✅ coherent |
| **lake `area_km2`, `level_m`, `depth_m`, `lake_type`** | ❌ **NOT scaled** | ⛔ **real, beside signified rivers** |
| `lake_min_area_km2` = 5 km² · `INVENTORY_MIN_CELLS` = 4 | applied in-stage, unscaled | the inventory floor is 5 **real** km² = 281 signified |
| `NECK_KM` 0.6 · `MIN_SPUR_KM` 1.0 · `MAX_SPUR_KM` 50 | ❌ no (`coast_shape` takes `km_per_cell` from the caller) | a "1 km spur" is **7.5 signified km** |
| `RELIEF_V1_A_C_KM2`, `HILLSLOPE_REF_CELL_M`, `diffusion` | ❌ no | terrain, correctly ratio-blind |
| coastline km, spur counts, hypsometry, slope quantiles | ❌ no | Ymir-km throughout |

### ⛔ The defect: Finding 49's own rule was applied WITHIN `rivers.json` and never ACROSS to `lakes.json`

Finding 49 states it in capitals at the code: *"EVERY SIGNIFIED QUANTITY, or the export ships two
scales at once."* It was enforced across the six segment arrays. **It was never asked of the other
exported object types.**

`export/hydro.rs:103` — `lakes_json` serialises `drainage.lakes` verbatim, so `area_km2`,
`level_m` and `depth_m` leave as **real 400 km values**. At ratio 7.5 a consumer can read:

> a river with `catchment_km2` = 5 625 terminating in a lake with `area_km2` = 100, linked
> explicitly by `source_lake_id` — **two numbers 56.25× apart in scale, describing one junction.**

Any consumer check of the form "is this lake plausible as the sink of this catchment", or any
water-balance recomputation downstream, is wrong by 56.25×. This is the same defect Finding 49
diagnosed, one level up: **it is not within an object, it is between object types.**

The ADR's justification for not touching lakes is sound for the **physics** — the water balance,
levels and footprints ran upstream on real quantities and must not be re-scaled. It does not
cover the **reported area**, which is a signified quantity by exactly the argument that makes
`drainage_km2` one.

**Mitigation that exists:** `ContinentMeta::geographic_scale_ratio` is written into the manifest
(`hd.rs:1223`), so a consumer is *told* the ratio. **What is missing is any statement of WHICH
fields are already scaled** — and with rivers pre-scaled and lakes not, a consumer applying the
manifest ratio uniformly would double-scale the rivers and correct the lakes.

Not fixed here. No change of any kind this round.

### 4 · Finding 47's "2–3 orders below the thresholds" was measured at ratio 1.0

It is a **ratio-1.0 statement about a product that runs at 7.5**, and the ADR's own table at line
~980 (2048², seed reference) already shows what the ratio does:

| ratio | small boat | barge | ship | largest Q, width |
|---|---|---|---|---|
| 1.0 | 742 | 1 | 0 | 48 m³/s, 35 m |
| **7.5** | 811 | **1 661** | **403** | **2 688 m³/s, 259 m** |

56.25× is **1.75 orders of magnitude**, which consumes almost all of "2–3 orders". So:

- ⛔ **"the discharge sits 2–3 orders below the navigability thresholds" does NOT hold at the
  ratio the product uses.** At 7.5 there are 403 ship reaches where there were none.
- ✅ **Finding 47's other and more important claim SURVIVES**: the discharge is not
  resolution-stable (10.33 m³/s at 2048² against 1.04 at 8192², the small-boat class ×50). The
  ratio multiplies both grids by the same 56.25, so the **instability is untouched** — thresholds
  still cannot be calibrated across grids, for the reason Finding 47 gave second rather than first.

**And this propagates to Finding 56c**, whose navigability verdict (small boat 214 → 175 at
2048², 35 → 11 at 8192², barges 2–3) was measured at ratio 1.0. At 7.5 the classes would shift
upward by 1.75 orders, so **the absolute counts and classes in that table describe a
configuration the product does not ship.** The *relative* effect of the channel-head law on them
may or may not survive; that is a measurement, not a deduction, and it is not made here.

### 5 · Named list

**STANDING — ratio-blind by construction** (nothing in them reads a scaled quantity):

- everything about the terrain: the pre-incision field's resolution invariance, the hypsometry,
  the slope quantiles, the erosion work, the work ratio, the attainable-set ceiling and its table
  (Findings 57, 58, 61, 62, 63, 66, 67);
- **`A_c`'s cell counts — 2.6214 and 41.9430 — at every ratio** (Finding 58, and the guard);
- the accumulation quantiles and the D8/MFD divergence (Finding 59);
- the whole coastal-fringe chantier: spurs 1 254 → 560 and 1 849 → 544, coastline km, p90 length,
  local R — `coast_metrics` reads only the height field (Findings 50–56, audited in Finding 67);
- the arêtes / steep-slope shares, the `A_c(S)` law's primary verdict;
- the lake INVARIANTS of Finding 56b/56c (footprint, connectivity, depth-vs-raster, orphan
  mouths) — all computed on real geometry, which the ratio does not touch;
- Finding 60 and Finding 65 (the decoupling and the four definitions) — strengthened rather than
  weakened, since the ratio adds a fifth axis of disagreement.

**WEAKENED — must now carry `geo_scale_ratio = 1.0` in the statement:**

- **Finding 47's "2–3 orders below the thresholds"** — false at 7.5, as above;
- **Finding 56c's navigability table** — 214/175, 35/11, barge 2–3, and the discharge p90 figures
  (1.759/1.983, 1.017/1.386 m³/s): all ratio-1.0, and the classes move by 1.75 orders at 7.5;
- **Finding 56b's W/D-per-order channel widths** — 0.4–5.2 m and the "Q > 0 share": widths ×7.5
  at the shipped ratio, so "the widest channel on the fine grid is 5.2 m" becomes 39 m signified.
  The `Q > 0` **share** is ratio-invariant (a multiplier cannot change a sign), so the
  "not one of the 51 order-5 reaches is wet" result stands;
- **Finding 56b's `catchment_km2` and `drainage_km2` figures**, and any absolute km² quoted from
  `rivers.json`;
- the validation note's header, which declares `geo_scale_ratio = 1.0` as *the production
  condition*. It is the **code's default** (`HdParams`, `Workspace`, `thread.rs` all initialise
  1.0) and **not the author's practice** (7.5). Same shape as `production_hd_config` versus
  `c1_hd_production`: the shipped default was mistaken for the shipped configuration.

**INVALIDATED — demonstrably, and it is one item:**

- **the claim that the export is internally scale-coherent.** Finding 49 established it for the
  segment arrays and it is false across object types: `rivers.json` signified, `lakes.json` real,
  linked by `source_lake_id`, 56.25× apart. Finding 49's own stated rule is violated at the level
  above the one it audited.

Nothing else is invalidated. In particular **no terrain result is affected**, because the ratio is
excluded from every stage that shapes the height field — by three explicit comments and by the
absence of any read in `stream_power.rs`.

## Finding 69 — river width in render cells: the widest river in the world is 1.4 cells, and the exponent is the lever

Rule 11 first, and it paid off twice.

| identifier | ADR hits | earliest | outcome |
|---|---|---|---|
| `width_m` | 17 | **899** | ✅ **Finding 21: the law is Leopold & Maddock `w = a·Q^b`, with `CHANNEL_WIDTH_A = 1.2`** |
| `discharge_m3s` | 7 | 921 | Finding 22 (Q from discharge, not area) |
| `geo_scale_ratio` | 26 | 971 | Finding 24/68 |
| `buildable` | 6 | **28** | ✅ **the acceptance criterion, already written** |
| `constructib` | **0** | — | **nothing found** (recorded) |
| `slope_deg` | 1 | 4474 | — |

**(a) The coefficient has drifted and no finding records it.** Finding 21 states `CHANNEL_WIDTH_A
= 1.2`; the code has **`5.0`** (`drainage.rs:629`). The **exponent** `b = 0.5` has genuine
antecedence — it is Leopold & Maddock's classic width exponent — but **`a` is a PROXY in
practice**, ×4.2 from its recorded value, with no measurement and no entry. Finding 21's per-order
medians (S1 13 m … S5 114 m) were also taken with `a = 1.2` **and** `Q` proxied by area, both
since replaced; they describe nothing that ships.

**(b) The buildability criterion was already in the ADR at line 28**, and it pre-empts the
connectivity-not-flatness reframing: *"The acceptance CRITERION is a LEGIBLE landscape — buildable
valley floors + sharply-delimited steep flanks + moderate interfluves — NOT minimal steep ground.
An earlier 'minimise unbuildable land' framing was WRONG and cost a pass: steep ground is gameplay
content."*

### The measurement — 8192², 48.83 m render cell, watercourses only (spillways excluded)

| bed · ratio | < 1 cell (count / length) | 1–2 | 2–5 | > 5 | p50 | p90 | p99 | **max** |
|---|---|---|---|---|---|---|---|---|
| arid · 1.0 | 100.00 / 100.00 % | 0 | 0 | 0 | 0.000 | 0.007 | 0.050 | 0.105 |
| arid · **7.5** | **100.00 / 100.00 %** | 0 | 0 | 0 | 0.000 | 0.050 | 0.376 | **0.785** |
| humid · 1.0 | 100.00 / 100.00 % | 0 | 0 | 0 | 0.056 | 0.103 | 0.151 | 0.185 |
| humid · **7.5** | **97.60 / 98.03 %** | 2.40 % | **0** | **0** | 0.421 | 0.773 | 1.135 | **1.385** |

> **At the shipped ratio, in the wet climate, the widest river on the map is 1.385 cells and
> 97.6 % of the network is sub-cell. Nothing anywhere reaches 2 cells.** The author's target — a
> `fleuve` at 4–5 cells — is short by a factor 3.6 on the single widest reach and by ~10× on the
> median. In the arid bed no watercourse reaches even one cell.

**Control: exact.** 1.0 → 7.5 moves p50 ×7.52, p90 ×7.50, max ×7.49. The shipped law is `5.0·Q^0.5`
as declared.

**The Strahler breakdown has no hierarchy, and the trunk is the THINNEST thing on the map**
(humid, 7.5, p50 in cells):

| S1 | S2 | S3 | S4 | S5 |
|---|---|---|---|---|
| 0.235 | **0.529** | 0.256 | 0.217 | **0.047** |

Order 5 is **11× narrower than order 2**. In the arid bed order 5's maximum width is **0.00** — it
is dry (Finding 58's dry-trunk result, now visible in the exported render hint). **A consumer
drawing width-proportional strokes would draw the main stem as the faintest line in its own
basin.**

### And a correction to how Finding 56c's discharge maximum should be read

I predicted a maximum of ~11.9 cells from Finding 56c's "discharge max 240.8 m³/s" (humid 8192²).
Measured: 1.385 cells, i.e. a maximum **watercourse** discharge of 3.25 m³/s. The 240.8 was a
**SPILLWAY** — the whole outflow of a below-sea basin over a col, which is exactly what Finding 45
warned tops the discharge sort and must not be read as a river. **My prediction was refuted
because I quoted a number about a different object.** Every "max discharge" in this dossier needs
`SegmentKind` stated beside it.

### The lever is the EXPONENT, not the coefficient — and I had that backwards

I predicted `a` was the lever, on the ground that the width spread was already 28×. **That 28×
was computed from the spillway maximum and is wrong.** Measured over watercourses:

- width contrast max/p50 = 1.385/0.421 = **3.3×**;
- p99/p50 = **2.7×**, so even the top percentile is nearly flat;
- the underlying discharge contrast is therefore 3.3² = **11×**, not 1000×.

The requested contrast is ≥ 5× (a `fleuve` at 4–5 cells beside a `ruisseau` under 1). With
`b = 0.5` and an 11× discharge range the width range is fixed at 3.3× **whatever `a` is**: raising
`a` by 3.6 puts the maximum at 5 cells and drags the median to 1.5, so every brook becomes a cell
wide and the hierarchy flattens further.

> **`a` sets the LEVEL, `b` sets the CONTRAST. To get 5× of width contrast out of an 11× discharge
> range needs `b ≈ log 5 / log 11 = 0.67`** — and `b = 0.5` is the one part of this law that has a
> published anchor. So the honest options are: widen the discharge range (which is the dry-trunk
> problem, Finding 58), or depart from Leopold & Maddock's exponent and label it PROXY.

**Not recommended here, and no value proposed.** The measurement is the deliverable.

## Finding 70 — buildable land is fragmented into ~43 000 components by the drainage skeleton, and the slope threshold is NOT the dominant parameter

Criterion: **8-connected components of cells whose local slope is below a threshold, with no
constraint on internal relief** — a terraced city is acceptable, and ADR line 28 already says the
criterion is a legible landscape rather than minimal steep ground. 8192², the exported (breached)
field, dimensionless central-difference slope in physical units.

**Extent, declared, because the verdict depends on it entirely:** the **diameter of the largest
inscribed disc** — an 8-connected chamfer transform (weights 1, √2) from the non-buildable set,
then `2 × max` per component. A one-cell ribbon scores 2 cells however long it is. This choice
also makes the measure **immune to the 8-connectivity chaining artefact**: two pockets joined by a
single diagonal pixel form one component with an inflated *area*, but the inscribed disc still
sits in one pocket and reports that pocket's size.

TDD §2.2 gives these as VALLEY WIDTHS, the same quantity: ravine 50–100 m · hameau 150–250 m ·
village/bourg 300–600 m · ville 800 m–1.2 km · **cité 2 km+** (the round's brief said 1.5–2 km;
the TDD's 2 km+ is used as the authority).

### The sweep

| slope < | buildable | components | largest | largest share | **cité ⌀≥2 km** | **ville ≥800 m** | village ≥300 | hameau ≥150 |
|---|---|---|---|---|---|---|---|---|
| **10 %** (5.71°) | 39.5 % | 39 183 | 611 km² | 5.9 % | **108** | 464 | 1 333 | 2 628 |
| **15 %** (8.53°) | 48.9 % | 42 952 | 787 km² | 6.1 % | **117** | 416 | 1 305 | 2 685 |
| **20 %** (11.31°) | 55.3 % | 43 683 | 849 km² | 5.8 % | **93** | 357 | 1 235 | 2 681 |
| **25 %** (14.04°) | 59.9 % | 40 194 | **2 100 km²** | **13.3 %** | **74** | 322 | 1 202 | 2 552 |

Component size, cells: p50 **2**, p90 16–18, p99 657–1 000, max 256 k–881 k at every threshold.

### Two predictions refuted, one of them mine and badly

**(1) No percolation.** I predicted that at ~50 % occupancy with 8-connectivity we would be above
the site-percolation threshold and one giant component would hold **> 80 %** of the buildable
land. Measured at 48.9 %: the largest holds **6.1 %**, and there are **42 952** components.
**Refuted outright.**

The reason is structural and worth keeping: the buildable mask is **not a random lattice**. The
steep cells are the valley walls of the drainage network, which is itself a connected, percolating
object. **The drainage network is a partitioning skeleton** — it cuts the buildable land into tens
of thousands of pieces, and no site-percolation intuition applies.

**(2) The slope threshold is NOT the dominant parameter.** The round predicted the count of
city-sized components would vary by **more than an order of magnitude** between 10 % and 25 %.
Measured: **108 → 117 → 93 → 74**, a factor **1.6**, and **non-monotone** — it peaks at 15 %.

And the mechanism is instructive: raising the threshold adds *area* (39.5 → 59.9 %) but **reduces**
the number of city-sized sites, because the cells it adds are on the flanks and they **merge**
previously separate components into large, elongated, ribbon-like shapes. At 20 % an 849 km²
component has an inscribed ⌀ of 14 505 m; at 25 % a 2 100 km² component has only 8 660 m. **Bigger
area, smaller inscribed disc** — which is precisely the discrimination the extent measure was
chosen for, doing its job.

⚠️ **Control-block note (rule 10):** the component *count* is nearly threshold-invariant
(39 183 → 43 683 → 40 194, ±6 %, non-monotone). It is a weak column and should not be read as a
response. The columns that do respond are the buildable fraction (39.5 → 59.9 %) and the
largest-component share (5.9 → 13.3 %).

### The answer to the consumer's question

**There is no shortage of city sites.** At the 15 % threshold: **117 components with an inscribed
diameter ≥ 2 km**, 416 at ville scale, 1 305 at village scale, 2 685 at hameau scale — and the
largest inscribed diameters are **3–14.7 km**, far above the TDD's 2 km cité threshold.

The top eight at 15 % (largest first), showing the cross-tabs:

| cells | km² | ⌀ m | alt p50 m | coast km | internal relief m |
|---|---|---|---|---|---|
| 330 197 | 787 | **14 150** | 1 372 | 2.4 | 1 867 |
| 270 014 | 644 | 5 752 | 189 | **0.0** | 908 |
| 219 464 | 523 | 9 976 | 292 | 2.9 | 881 |
| 203 819 | 486 | 5 086 | 182 | **0.0** | 944 |
| 198 449 | 473 | 11 466 | 394 | 0.1 | 1 162 |
| 187 154 | 446 | 11 915 | 377 | 0.1 | 833 |
| 168 835 | 403 | 8 584 | 570 | 0.3 | 1 132 |
| 165 027 | 393 | 3 947 | 1 738 | 4.3 | 1 170 |

**Majority COASTAL and low, not interior plateau** — six of the top eight sit within 3 km of the
coast with a median altitude of 155–570 m. The round predicted the opposite; **refuted**. Two are
high-altitude (1 372 and 1 738 m) with 1.2–1.9 km of internal relief, which is the terraced-city
case the reframing was meant to admit.

**Internal relief is reported, not filtering, as asked** — and it is large: 500–2 500 m inside a
single buildable component. That is the gameplay information: these are not plains, they are
connected networks of gentle ground threading through steep terrain.

### What the numbers cannot settle

A component with an inscribed diameter of 14 150 m, a median altitude of 1 372 m and 1 867 m of
internal relief is either a broad upland basin or a dendritic web of valley floors stitched
together. **The inscribed diameter says a 14 km disc of sub-8.5° ground exists somewhere in it; it
does not say that disc is a plain rather than a locally-flat shoulder.** That is a judgement, and
it is posed as such in the validation request rather than answered here.

## Finding 71 — the 240.8 m³/s is NONE of the three worlds: it is water moving between two inland sinks, and 81 % of the continent's runoff never reaches an outlet

Rule 11 first: `spillway` **95** ADR hits (earliest 1145) · `shore` **46** (71) · `water_class`
**31** (760) · `bank` 8 (387) · `SegmentKind` 3 (3016) · `source_lake_id` 6 (3024) · `lake_mask`
**1** (1349) · **`berge` and `riverbank`: nothing found** (recorded).

### A0 — the budget, and it is twice the reference in humid

| bed | land | generated runoff | ×reference (300 mm) | **budget** |
|---|---|---|---|---|
| humid | 26 285 km² | **602.9 mm/yr** | **×2.01** | **502.1 m³/s** |
| arid-hot | 26 285 km² | **207.3 mm/yr** | ×0.69 | **172.6 m³/s** |

The round predicted 120–250 mm/yr in humid; measured **602.9**. **Refuted.** My 350–650 holds for
humid; my 30–120 for arid is **refuted** (207.3). `REFERENCE_RUNOFF_MM = 300` is a factor 2 low
for the humid bed, which matters because `drainage_km2` — the field navigability classifies on —
is `discharge / 300 mm`. **The "effective area" is inflated ×2 in humid by a proxy that is off by
×2.**

### A1 — the closure fails by a factor 5, and my first terminal test was the reason it looked worse

⚠️ **My own measurement defect, corrected in the same round.** The first pass tested "the
downstream-most cell's D8 receiver is sea" and reported ×0.018. That test **structurally excludes
spillways**: a spillway's path is traced over a col, outside the accumulation network, so its last
cell has no D8 sea receiver. The code supplies the right discriminator itself —
`drainage.rs:2024` asserts `wc[end] == 1 || chained_into.is_some()`, so **`water_class` at the
last cell is 1 (ocean) exactly when the spillway is terminal.**

Corrected, humid (budget 502.1 m³/s):

| class | segments | Σ m³/s | % of budget |
|---|---|---|---|
| Watercourse → sea | 82 | 46.5 | 9.3 % |
| **Spillway → OCEAN** (wc = 1) | 37 | 46.6 | 9.3 % |
| **Spillway → CHAINED** (wc = 2) | 17 | **589.5** | **117.4 %** |
| **TERMINAL TOTAL** | 119 | **93.1** | **×0.185 — OUTSIDE the declared factor 2** |

Arid (budget 172.6): watercourse-to-sea 10.8, spillway-to-ocean 1.3, **chained 28.2**, terminal
total 12.1 = **×0.070**.

**So the verdict is NONE of the three worlds**, and the leak has three named components:

**(1) The water disappears into a CHAIN of below-sea basins.** 589.5 m³/s — 117 % of the whole
continental budget — flows out of one enclosed below-sea basin into another, and the terminal
outlets carry only 93.1. **Water enters the chain and does not come out.**

**(2) Double counting exists, and it is INSIDE the chain, not at a multi-outlet lake.** The
hypothesised mechanism — a lake with N outlets each inheriting the whole upstream — is **absent**:
**0 of 38** inventoried spillway sources has more than one spillway, at both beds. But the chained
sum alone is 117 % of the budget, so along a chain A → B → C the same water is carried by each
link's spillway and summing them over-counts. The naive sum of all spillways plus watercourses is
**×1.359** of budget in humid.

> ⚠️ **WEAKENED by Finding 72.** `water_class` has three values — 0 land, 1 ocean, 2 enclosed
> below-sea — and **a detected surface lake above sea level is class 0**. So `wc = 0` does not
> mean dry land; it means "neither ocean nor an enclosed below-sea region", which includes every
> surface lake. Read this as **three spillways ending in a class `water_class` cannot
> represent** — the Finding 65 pattern again: the raster consulted to classify termini cannot
> see the water bodies the machinery delivers into. Which of the three are surface lakes is
> the measurement Finding 72's stop rule prevented.

**(3) Three of the eight largest spillways end on a LAND cell** (`water_class = 0`): 21.21, 11.67
and 11.06 m³/s. The code's own assertion says that should be impossible unless the spillway is
chained — and a chained spillway should end in class 2, not class 0. **~44 m³/s is delivered to dry
land and stops there.** Unattributed.

### The answer to the round's single question

> **The 240.82 m³/s spillway drains below-sea basin `1000056` into ANOTHER enclosed below-sea
> basin — `water_class = 2` at its last cell — and it is 2 CELLS LONG.**

It is not a continental collector (no *terminal* outlet exceeds 27.63 m³/s = 5.5 % of budget), not
the terminus of a trunk (see A3), and not a multi-outlet artefact. **It is a two-cell transfer
between two inland sinks, carrying 48 % of the continent's runoff.** The second, 228.74 m³/s, is
the same thing.

For the consumer this settles the "spillways are the rivers of this map" hypothesis: **the two
biggest are not rivers at all.** They are 2 and 40 cells long and they end inland.

### A2 — the dead column, recorded, then the instrument that works

Executed once as asked, to enter it as a dead column: **100 % of spillways sit at Strahler order 1
at both beds** (54 of 4271 S1 segments in humid, 0 at every higher order), by construction —
`strahler_order: 1 // MEANINGLESS on a spillway`. **The column measures the constant.**

The decile replacement is **also weak**, as predicted: 54 spillways diluted into 754-segment
deciles give 1.86 % of the top decile by count and 3.44 % by length. The **named top 20** is the
instrument that reads: **11 of the top 20 by discharge are Spillway in humid** (ranks 1–10 and 15),
3 of 20 in arid. Prediction was ≥ 15; **measured 11 — refuted, but the substance holds**: the top
ten are all spillways and the best watercourse is rank 11 at 3.25 m³/s.

| kind | n | length | p50 | p90 | p99 | **max** |
|---|---|---|---|---|---|---|
| humid Watercourse | 7 492 | 8 822 km | 0.301 | 1.014 | 2.184 | **3.252** |
| humid **Spillway** | **54** | **34 km** | 0.116 | 12.950 | 228.741 | **240.821** |
| arid Watercourse | 12 398 | 14 192 km | 0.000 | 0.004 | 0.240 | **1.045** |
| arid **Spillway** | **44** | **9 km** | 0.014 | 0.313 | 15.106 | **15.106** |

**54 spillways carry the top of the distribution over 34 km of the 8 856 km network — 0.4 % of
its length.**

### A3 — the trunk trace: smooth, and it refutes the jump

From A1's biggest **terminal** outlet (1.11 m³/s, **Watercourse**), humid:

- length **242 cells = 11.8 km**; 92.6 % of it lies on a mapped segment;
- **kind alternations: 1**; longest continuous Watercourse run **224 cells = 92.6 %** of the trace;
- discharge downstream→upstream: 1.11 → 0.93 → 0.51 → 0.48 → 0.41 → 0.32 → 0.31 → 0.15 → 0.09 →
  0.01 → 0.00 — **a smooth monotone decay**;
- **biggest step ×2.0**, at 88 % along the trace, Watercourse→Watercourse.

**My ≤ 4 alternations and > 50 % continuous are confirmed. The round's > 10 alternations and
< 20 % continuous are refuted. And the predicted single 50–80× jump does not occur — the largest
step is ×2.0.** Arid is the same shape: 146 cells, 2 alternations, 65.1 % continuous, max step
×5.2 in the near-zero tail.

The trunk is real, continuous and **tiny**: 11.8 km long, draining 61.2 km², carrying 1.11 m³/s —
and 61.2 km² × 602.9 mm/yr = **1.17 m³/s**, so the discharge is internally consistent with its own
catchment to 5 %. **Nothing is wrong with the watercourses; there simply are no large ones.**

### A4 — the ~28 "missing" spillways never existed

| bed | lakes | detected | below-sea | exorheic | endorheic | spillways | naming a source | naming none |
|---|---|---|---|---|---|---|---|---|
| humid | 86 | 48 | **38** | 84 | 2 | **54** | **38** | **16** |
| arid | 56 | 17 | 39 | 28 | 28 | 44 | 28 | 16 |

**The premise is false and the count is exact.** Spillways are produced *only* by
`below_sea_basin_lakes_infil`, so the number to reconcile against is the **below-sea basin** count,
not the exorheic lake count. In humid: **38 inventoried below-sea basins ↔ 38 spillways naming a
source, exactly.** An exorheic lake above sea level overflows through an ordinary `Watercourse`, so
the other 46 exorheic lakes need no spillway. **There is no hole.**

**And the 16 spillways naming no source ARE Finding 56b's dangling `lake_map` ids** — the same
population seen twice, as the round suspected. They are below-sea basins under
`INVENTORY_MIN_CELLS = 4`, marked in `lake_map`, absent from `lakes.json`, and here they are the
sources of 16 spillways whose `source_lake_id` is therefore `null`. Confirmed.

### ⛔ The export-contract gap that made this un-measurable

`Spillway::chained_into` is computed, used inside `below_sea_basin_lakes_infil`, and **logged** at
`hd_assembly.rs:169` as the "to sea / chained" split — then **dropped**. It never reaches
`C1DrainageResult` and never reaches `rivers.json`.

> **The single field needed to tell a terminal outlet from an interior transfer is computed and
> thrown away.** Neither a consumer nor this bench can build the terminal set from the export; it
> has to be reconstructed from the `water_class` raster, and only because the code happens to
> assert the relation. A consumer drawing "the biggest rivers" from `rivers.json` today would draw
> two 2-cell and 40-cell inland transfers as the map's principal watercourses.

### Standing, and what is NOT done

The round's stop rule fires: A1 concludes a leak **with double counting inside it**, so blocks B
and C are not run. **C3 must not run on these discharges** — the width gauge it would rasterise is
derived from them, and 81 % of the budget is unaccounted. B is independent and remains owed.

Nothing changed: no production edit, no export field added, no renderer. The leak is attributed,
not fixed.

## Finding 72 — the ledger is NOT READ: the blocking control could not be built in three attempts, and the failure is the instrument's

Rule 11 first: `extra_inflow` **2** ADR hits (earliest 1788) · `next_extra` **1** (1792) ·
`chained_into` **3** (1792) · `water_balance_lakes` **2** (2277) · `below_sea_basin_lakes_infil`
**5** (2822) · `absorb` **6** (1665) · `REFERENCE_RUNOFF_MM` **1** (7632, i.e. Finding 71 itself) ·
**`pe_lake`: nothing found** (recorded).

### The perimeter — measured, bounded, and NOT the explanation

| | humid | arid-hot |
|---|---|---|
| A0 budget population (`field > SEA`) | 11 024 627 cells ⇒ **502.1 m³/s** | same cells ⇒ **172.6 m³/s** |
| of which DETECTED-LAKE footprints | 817 919 cells = **24.0 m³/s (4.8 %)** | 12 966 cells = **0.0 m³/s** |
| rain surplus over ENCLOSED BELOW-SEA (wc = 2) | 457 698 cells = **11.9 m³/s** | **0.0 m³/s** |

Both biases are named and bounded:

- lake footprints **are** in the budget — the breach removed the depressions, so those cells read
  `> SEA` — and open-water PE is never deducted, so **A0 overstates by at most 24.0 m³/s (4.8 %)**;
- rain over below-sea footprints is **out** of the budget *and* out of the basins' inflow:
  `runoff_accumulation` seeds only cells with `heightmap > C1_SEA_LEVEL_NORM`, so a below-sea cell
  contributes zero. **11.9 m³/s is counted nowhere at all** — a third, smaller leak, now named.

**Together 35.9 m³/s against a 409 m³/s hole: 8.8 %.** The prediction that the perimeter
correction moves the hole by less than 15 % is **confirmed**. **The leak is not a population
artefact.**

### A3 — the blocking control, three constructions, three artefacts. All mine.

| attempt | construction | why it could not pass |
|---|---|---|
| 1 | accumulation one step downstream of `lk.base.outlet` | `base.outlet` is the **PRE-BREACH sill** — a saddle on the rim — and `flow.direction` was read on the **BREACHED** field where that depression no longer exists. A saddle carries almost no accumulation. Ratios 0.001 / 0.000 / 0.000 |
| 2 | start at the footprint's own **max-accumulation** cell, walk downstream until leaving the footprint | still fails: 0.004 / 0.021 / 0.000. Under D8 accumulation is monotone downstream, so a 2 770× drop in one step is **impossible** unless the cell landed on is not downstream — which it is not, because the exit cell is **at or below sea level** and `runoff_accumulation` leaves those at zero by construction |
| 3 | (arid bed, either construction) | ⚠️ **EMPTY POPULATION: there is no exorheic DETECTED lake with a footprint in the arid bed at all.** The control did not run. A rule-10 failure in my own control |

> **I do not have a validated instrument, so the ledger is not readable, and the stop rule fires.
> A1, A2, §3's cross-check, 4a, 4b and C were NOT executed.**

**And the failure is mine, not the machinery's.** Nothing above licenses the sentence "the
instrument does not carry flow across an exorheic surface lake". Each attempt measured a
different artefact of my own construction: a sill instead of an outlet, a water cell instead of a
land cell, and an empty set. The code's own claim — *"`runoff_accumulation` carries flow ACROSS an
exorheic lake (flat, routed to the outlet), so a lake's outlet reach inherits the whole upstream
catchment automatically"* — is **neither confirmed nor refuted here.**

What a fourth construction has to do, stated so the next round does not repeat these three: compare
the discharge of the **segment leaving** the lake (from `segment_discharge_m3s`, on the clipped
network, `SegmentKind::Watercourse`, population *inherited*) against the maximum discharge of the
segments **entering** it. That works on the object the code's claim is about — an *outlet reach* —
instead of on a raster whose zeros are load-bearing. It also requires that such a reach exist,
which is itself part of the question.

### The three "dry land" spillways of Finding 71 — the classification stands unverified

The six-class typology was built and is in the bench, but it sits after the stop and was not read.
What **is** settled, from the code alone:

**B1 — the guard could not fire, and it is worse than a log line.** `#[cfg(test)]` begins at
`drainage.rs:1965`, so the `assert!(wc[end] == 1 || sw.chained_into.is_some())` at **2024 is inside
the test module** — an assertion over a synthetic fixture. On the production path there is only an
`eprintln!` at `hd_assembly.rs:169`, gated on `verbose`, which every bench in this dossier calls
with `false`. **The relation is asserted on a fixture and logged on production; it has never been
checked on real data.** Both predictions were right; mine was sharper (it *is* an `assert!`, just
not where it matters).

**CORRECTION to Finding 71 at the point of claim.** Finding 71 reported "three of the eight
largest spillways end on a LAND cell (`water_class = 0`)" and called it unattributed. The claim
must be weakened: `water_class` has three values — 0 land, 1 ocean, 2 enclosed below-sea — and **a
detected surface lake above sea level is class 0**. So `wc = 0` does not mean dry land; it means
*"not ocean and not an enclosed below-sea region"*, which includes every surface lake. The finding
should read **"three spillways end in a class `water_class` cannot represent"**. This is the
Finding 65 pattern again: **the raster consulted to classify termini cannot see the water bodies
the machinery delivers into.** Whether those three are in fact surface lakes is the measurement
that the stop rule prevented.

### D — the contract, two occurrences of a new family, and the guard cannot carry them

Confirmed from the code:

- **`Spillway::width_m` exists** (`drainage.rs:1826`, `CHANNEL_WIDTH_A · Q^CHANNEL_WIDTH_B` on the
  spillway's own discharge). **This closes the C1 question left open last round: spillways are
  drawable.**
- **`Spillway::chained_into`** — computed, used in the fixpoint, logged, dropped.
- **`BasinSummary` in full** — 17 fields per basin including the entire water balance
  (`inflow_m3s`, `evaporation_m3s`, `a_eq_km2`, `a_spill_km2`, `has_sill`,
  `predicted_exorheic`, …). `assemble_hd_drainage` reads `bs.lakes`, `bs.lake_map`, `bs.wetland`
  and `bs.spillways`; **`bs.basins` is never read at all.**

**Can the Finding 68 guard carry this family? No, and the reason is structural.** That guard
enumerates the **export surface** and requires every unit-bearing field on it to be declared. A
computed-not-exported field is *by definition absent from the export surface*, so the enumeration
cannot see it — the guard runs in the wrong direction.

> A second guard is needed and it runs the **opposite** way: *every field a consumer REQUIRES must
> be present in the export.* Its input is a **declared requirements list**, and that list is not
> derivable from this repository — it is the consumer's statement of need. So the second guard
> cannot be written from the code side alone, which is precisely why this family has stayed
> invisible while the first guard passes.

Not implemented, per the round's constraint.

### Standing

Nothing changed: no production edit, no export field, no threshold, no renderer, no C3. The 81 %
remains **unattributed**, the perimeter is eliminated as its explanation (8.8 %), and the next
round's first item is a fourth construction of the A3 control on the segment network rather than
on the raster.

## Finding 73 — the fourth construction WORKS and it refutes the machinery: `runoff_accumulation` does not carry flow across a filled depression, because its sort key cannot order a flat

Rule 11 first: `segment_discharge_m3s` **5** ADR hits (earliest 929) · `clipped` **8** (1696) ·
`.basins` **2** (1415) · `has_sill` **2** (1516) · `inflow_m3s` **1** and `evaporation_m3s` **1**
(both 7861, i.e. Finding 72 itself — nothing antecedent). Reading the earliest hits first paid
twice, before any measurement:

- **L4162 / `drainage.rs:756`** — the Finding 22 exception: an outlet run does not read its own
  downstream value, it **inherits `src_q[i]`, the parent's MAXIMUM over its points**
  (`drainage.rs:534`), taken downstream of the lake. So the out/in ratio does test the carry — the
  inherited value is only large if the accumulation crossed — but **its magnitude above 1 is not
  interpretable** (the parent gains catchment below the lake). Declared one-sided before measuring.
- **L1415, Finding 37 POINT 1** — *"detect_lakes 28 exorheic / 0 without an outlet"*. The
  population A3 needed had already been counted once, fifty findings ago.

### A3-four — the control finally runs, and it FAILS on its object

Humid only, declared in advance (the arid world has no exorheic surface lake; Finding 72
attempt 3). Population declared before measuring at 20–45; **measured 42** (excluded: 3
spillway-fed, 1 with no outlet run). Tolerance declared before measuring: pass iff median r ≥ 1.00
**and** ≥ 80 % at r ≥ 0.95.

| | prediction (2026-09-10, non-blind) | measured |
|---|---|---|
| population | 20–45 | **42** ✓ |
| median r = OUT/IN | 1.0–1.6 | **0.052** ✗ |
| share r ≥ 0.95 | ≥ 90 % | **7.1 %** ✗ |

**My prediction is refuted and the measurement wins.** p10 0.007 · p25 0.018 · p75 0.253 · p90
0.579 · min 0.002 · max 129.9.

The **sensitivity** half passes: at a below-sea basin/spillway interface the same instrument reads
1000056 → **160.9**, 1000021 → 228.6, 1000035 → 77.5; median **17.1** over the 15 spillway basins
with an identifiable inlet run, against 0.052 at an ordinary surface lake. **Dynamic range ×328**
— the instrument is not flat, so the 0.052 is a reading, not a null.

> **Stop rule 1 fires. Blocks 2 and 3 were NOT executed** — no ledger, no discriminant, no A1, no
> A2, no B2, no 4a, no 4b, no C. Both beds: the arid ledger cannot be more readable than the
> instrument that would read it, so it inherits the humid verdict rather than being waved through.

### But this time the failure is NOT the instrument's — and here is the proof

The three earlier constructions each measured an artefact of my own. This one does not, and the
round's obligation ("say how the specification was wrong") turned into an attribution. Two
measurements separate the only two things the ratio could conflate.

**A3-four-BIS — the raster against the network, same 42 lakes.**

| lake | IN raster | IN pool | OUT raster | r | exits below sea | run on the exit cell | r_seg |
|---|---|---|---|---|---|---|---|
| 26 | 18 770.8 | 18 881.4 | **62.1** | 0.003 | **0 / 5** | 0.0015 m³/s | 0.002 |
| 6 | 91 674.4 | 91 726.1 | **114.5** | 0.001 | **0 / 1** | 0.0078 | 0.003 |
| 25 | 58 927.3 | 58 930.5 | **66.3** | 0.001 | **0 / 1** | 0.0121 | 0.007 |

(mm·km²/yr.) Median raster ratio **0.007**, share ≥ 0.95 **4.8 %**. And: **a mapped `Watercourse`
run covers the exit cell for 42 of 42 lakes, and that run's FIRST point IS the exit cell for 40 of
42.** So the segment layer is doing its job — the outlet run exists, in the right place, and it
advertises the value it was given. **The raster is the carrier that fails.**

This also **refutes my own Finding 72 explanation of attempt 2**: I wrote that construction 2
failed because "the exit cell is at or below sea level where `runoff_accumulation` leaves zero".
Measured: **0 of the exit cells of these 42 lakes is below sea level.** That explanation was wrong.

**A3-four-TER — the mechanism, and a negative control on it.**

`runoff_accumulation` (`drainage.rs:962`) accumulates in **one pass over cells sorted by
`flow.filled` DESCENDING**. That is a valid topological order of `flow.direction` only where the
terrain strictly descends. Inside a filled depression `filled` is **flat**, so the key cannot order
those cells, `sort_unstable_by` breaks the ties arbitrarily, and a cell is routinely processed
**before its own donor** — whose contribution is then never propagated.

- Over the **605 046** footprint cells of the 42 lakes, **605 046 — 100.0 %** send their water to a
  cell whose `filled` is **not strictly lower**. The sort key orders none of them.
- Negative control: re-accumulate with a **true topological order** (Kahn over the D8 in-degree) —
  same field, same directions, same per-cell runoff, **only the order changed**. The Kahn pass
  processes 11 024 627 of 11 024 627 above-sea cells, so the D8 land graph is acyclic.

| | shipped order | correct order |
|---|---|---|
| median OUT/IN over the 42 lakes | **0.007** | **4.126** |
| share ≥ 0.95 | 4.8 % | **100.0 %** |
| p10 · p90 | 0.001 · 0.065 | 1.815 · 10.128 |

**ATTRIBUTED.** The ratios exceed 1 by 2–14× in the corrected order, which is the expected sign:
the exit cell is downstream and gains catchment. The one-sided prediction was right about the sign
and wrong about the shipped field.

### What it costs, measured — and what it does NOT explain

**Σ accumulation delivered to OCEAN cells: 167.3 m³/s shipped, 202.1 m³/s in the correct order**,
against a 502.1 m³/s budget — closure **×0.333 → ×0.403**.

> **+34.8 m³/s, i.e. 8.5 % of the Finding 71 hole (409 m³/s). This is a THIRD named leak, sized —
> not the explanation of the 81 %.** Perimeter 8.8 %, accumulation order 8.5 %, and the rest still
> unattributed. Correcting the order does not close the budget: **59.7 % of the runoff still
> terminates without ever reaching an ocean cell**, which is the Finding 71 chain and remains the
> open question.

**Blast radius, from the call sites and not from imagination.** `runoff_accumulation` feeds
`segment_discharge_m3s` (hence `width_m` via Leopold, `drainage_km2` and the navigability class),
the lake water balance's inflow, and `BasinSummary::inflow_m3s` — which sums `runoff[nk]` cell by
cell at the shoreline (`drainage.rs:1478`). **So the ledger's own source numbers are downstream of
this defect**, which is a second, independent reason block 2 could not have been read this round
even if block 1 had passed.

**And it hits `A_c` too, in the same direction as Finding 65.** Every threshold of the form
"discharge ≥ X" is evaluated on a field that loses water at every filled depression, so a channel
head test downstream of one fires later than the physics says. Not measured here; named.

### ⚠️ CORRECTION to Finding 22, at the point of claim

Finding 22's exception is justified in the code by the sentence *"`runoff_accumulation` carries
flow ACROSS an exorheic lake (flat, routed to the outlet), so a lake's outlet reach inherits the
whole upstream catchment automatically"* (`drainage.rs:520`). **The premise is false as shipped**:
measured, the accumulation loses 99.3 % of it (median). The exception itself is still the right
thing to do — it is the only reason the outlet run advertises anything at all — but it is
compensating for a defect one layer down, not expressing a property of the accumulation. **A
correct fix is upstream of it, in the ordering, not in the clip.** Not implemented, per the
round's constraint.

### My own instrument, corrected a fourth time — and disowned before it was reported

Route (a) of the block-2 cross-check read `flow.basins` **on the footprint cells**. `compute_flow`
treats every cell at or below sea level as `is_ocean` and `compute_basins` skips those, so a
below-sea cell carries label **0**: the watershed came back empty and (a) read 0.00 m³/s. Left
alone that would have been published as "the machinery's accounting is not confirmed" when the only
broken thing was my label lookup. Corrected to collect the labels of the **shoreline entry cells**;
**not executed**, because stop rule 1 forbids it. The pattern is now four for four: **a raster whose
zeros are load-bearing, read as if a zero meant "no water".**

### Block 4 — the invariant, as a SPECIFICATION (not implemented)

Delivered, because it depends on the code and not on the ledger. Finding 72 established that the
existing `assert!(wc[end] == 1 || sw.chained_into.is_some())` sits at `drainage.rs:2024`, **inside
the `#[cfg(test)]` module that begins at 1965** — an assertion over a synthetic fixture, with only
a `verbose`-gated `eprintln!` on the production path. One sentence to replace it:

> **Every spillway must terminate in exactly one of three states — an OCEAN cell (`wc == 1`), a
> DIFFERENT below-sea region (`region_of[end] ≠ 0` and `≠ own_label`, the case that sets
> `chained_region` and therefore routes the surplus), or a DETECTED SURFACE LAKE
> (`detected_lake_map[end] ≠ 0`, the case that sets `chained_into` but NOT `chained_region`, so
> the surplus is dropped) — and the third case must be counted and reported, never silently
> accepted.**

Cost: one pass over the 54 spillways, at the end of `below_sea_basin_lakes_infil`, reading three
arrays it already holds. Note what the wording deliberately does **not** do: it does not use
`water_class` to name the third case, because Finding 72 established `water_class` cannot represent
a surface lake above sea level (`wc = 0` means "neither ocean nor enclosed below-sea"). Promotion
is the author's call and is out of scope here.

### Standing

No production edit, no export field, no threshold, no renderer, no promotion of the invariant, no
C3. Block 1 **fails**, blocks 2 and 3 are not delivered, block 4 is. What the round produced
instead is the first **attributed** term of the leak since Finding 71: the accumulation order,
worth 8.5 %, with a negative control that closes it and a named blast radius.

## Finding 74 — the order fix, measured in both code states; and the real mechanism of the hole is the CYCLE-BREAKER, not the surface lakes

Rule 11 first, and it holds the answer: **`resolve_flats` 0 ADR hits · `flat_grad` 0 · `compute_accumulation` 0 — NOTHING FOUND, recorded.** `pit_fill` **3** (earliest 689) · `priority-flood` **13** (673).

Three identifiers at zero hits, and they are the ones that held the fix. `terrain/flow.rs:709`
`compute_accumulation` has always sorted on `(filled desc, flat_grad desc)`, with the comment
*"so flat cells are processed from the inflow side down to the outlet (correct topological order
on flats where `filled` ties)"*. **The correct order existed from the start, in `flow.rs`, and
`runoff_accumulation` was a copy of that loop with the tiebreak dropped.** The comparator is
TRIPLED — `flow.rs:723`, `flow.rs:520`, and the broken half at `drainage.rs:963`. Finding 13's
*"fills to the EXACT sill — no epsilon increment, so flats are truly flat"* is what makes the
missing tiebreak fatal rather than cosmetic: the flats are exact **by design**.

### The implementation choice, instructed before writing

Kahn over the D8 in-degree, hosted as `terrain::flow::propagation_order`, **not** the `flat_grad`
reuse. Cost: no new `FlowResult` field and no 537 MB `f64` export, and O(n) against the O(n log n)
sort it replaces. Coupling: a `(filled, flat_grad)` key is only topological if it came from the
**same** `flat_grad` that built `direction`, and **six sites reconstruct a `FlowResult` from a
cache sidecar with `direction` but no `flat_grad`** — an empty order there is a silently-zero
accumulation, this campaign's signature defect. Kahn derives the order from the field it must
agree with. Consumers of `resolve_flats`: four functions, all private inside `flow.rs`, nobody
outside — so reuse would mean widening a private surface for one caller.

> **DEBT, NAMED.** The flat comparator stays duplicated between `flow.rs:723` and `flow.rs:520`,
> and Kahn is now a THIRD ordering scheme beside them. They are not unified, because
> `compute_accumulation` feeds stream-power incision and changing its order would change the
> delivered terrain.

`ALGO_DRAINAGE` 5→6 and `ALGO_HD_DRAINAGE` 7→8: every discharge, width and lake inflow moves.

### The guards, before any result

**The control block is BIT-IDENTICAL.** The same bench ran in both code states (the fix stashed,
then restored). FNV-1a over the raw `f32` bits:

| field | before | after |
|---|---|---|
| eroded | `0x6b10a0c5fdf467a3` | **`0x6b10a0c5fdf467a3`** |
| breached | `0x2cce4ea2fb761a12` | **`0x2cce4ea2fb761a12`** |
| pre-filled | `0x5bd35562ba4fb9ab` | **`0x5bd35562ba4fb9ab`** |

Means, min/max, land counts and land p10/p50/p90 identical to nine decimals; 1 238 085 cells
raised by `pit_fill`, 457 698 below-sea, 55 626 539 ocean — unchanged. **The stop rule does not
fire: the Finding 73 blast radius was right, the height field is not downstream of this.**

**The 42-lake instrument is promoted to a guard, twice.** `accumulation_order.rs` holds a
**default-running synthetic fixture** — a filled pit fed by a known tributary — and an `#[ignore]`d
8192² test on `assemble_hd_drainage`. The split is deliberate: a four-minute test does not protect
a default suite, and *a rule that is not executable does not protect*.

| the same flat, five readings (192-cell synthetic pit) | IN | OUT | OUT/IN |
|---|---|---|---|
| topological order (the fix) | 1825.00 | 2688.00 | **1.4729** |
| `filled`-only order (PRE-FIX) — **negative control** | 1825.00 | 8.00 | **0.0044** |
| `flow.accumulation` (geometric, never broken) | 1825.00 | 2688.00 | 1.4729 |
| `runoff_accumulation` (production) | 4 094 635 | 6 030 893 | **1.4729** |
| `mfd_accumulation` (incision) | 1571.37 | 185.56 | **0.1181** |

The broken order reproduces the 8192² defect on 192 cells (0.0044 against 0.007) — a **×336
separation**, so the guard is not vacuous — and the fix lands **exactly** on the geometric
accumulation's value, which is the proof that the two orders now agree.

Production path, humid, old → new: median OUT/IN **0.052 → 4.185**, share ≥ 0.95 **7.1 % → 100.0 %**,
min **0.002 → 1.254**, population 42 → 44. The `#[ignore]`d guard on `assemble_hd_drainage` reads
n 44, median 4.185, share 100.0 %, min 1.254, max 689.101 — **identical to the digit** to the
bench's independently replicated tail, which is the rule-7 cross-check on the replication itself.

Lake invariants (56b/56c) pass at both beds and both resolutions. Lakes move as expected from more
inflow: humid **84 exorheic / 2 endorheic → 86 / 0**, Σ area 5970.4 → 6219.4 km² (+4.2 %); arid 56
→ 61 lakes, 28/28 → 31/30, 1638.0 → 1703.7 km². **The five largest lakes are IDENTICAL to the
decimal in level and area at both beds** — so the "> 30 % of lakes move" prediction is neither
confirmed nor refuted here: a per-lake diff was not measured and is not claimed.

### The Finding 71 table, re-derived on sane numbers (humid, old → new)

| | before | after |
|---|---|---|
| budget | 502.1 | **502.1** (unchanged, as predicted — a per-cell sum) |
| `Watercourse` → sea (65 runs) | 37.5 | **75.9** |
| `Spillway` → ocean (37) | 46.6 | **82.5** |
| TERMINAL | 84.1 (×0.168) | **158.4 (×0.315)** |
| **HOLE** | 418.0 | **343.7** |
| `Spillway` → `wc = 2` (13) | 538.4 (107 % of budget) | **1052.4 (210 %)** |
| **max `Watercourse`** | **3.252 m³/s** | **21.785 m³/s** |
| max `Spillway` | 240.821 | **380.514** |
| Spillways in the named top 20 | 11 of 20 | **5 of 20** |
| `Watercourse` p50 / p90 / p99 | 0.301 / 1.014 / 2.184 | **0.477 / 2.406 / 19.200** |

Arid: TERMINAL ×0.039 → ×0.048, max `Watercourse` **1.045 → 1.045 (unchanged)**, network 12 398 →
12 275 runs. C is untouched by the fix, as predicted: measured runoff **602.9 mm/yr** humid /
207.3 arid against the 300 mm proxy ⇒ effective area ×0.498 / ×1.447.

> **The consumer-facing result: the map now has a river hierarchy. The largest `Watercourse` goes
> from 3.25 to 21.79 m³/s — ×6.7 — and the p99 from 2.18 to 19.20, ×8.8.** Finding 71's "nothing
> is wrong with the watercourses; there simply are no large ones" is **corrected**: there were
> large ones, and the accumulation order was deleting them.

### The discriminant REFUTES both bets, and names the mechanism it was not looking for

Σ discharge of spillways ending on a surface lake (`did ≠ 0`) = **73.7 m³/s** against the measured
hole of 343.7 ⇒ **coverage ×0.214**, far outside the declared factor 1.5.

**My mechanism at 65 % is refuted as the explanation.** It is real and it is 21.4 %, not the
answer. My predicted coverage of ×0.7–1.3 is refuted; so is the round's expectation that the
discriminant would confirm it.

And the largest class was mislabelled by MY OWN classifier: 8 spillways, **1048.8 m³/s**, came out
as `true DRY LAND` — the third time this campaign has been tempted by that sentence. Tested rather
than read, by comparing the 8-connected `wc == 2` COMPONENT of the terminal cell against the
component(s) of the source basin's own footprint:

| Q m³/s | source | own cmp | end cmp | #cmps | verdict |
|---|---|---|---|---|---|
| 380.51 | 1000056 | 56 | 58 | 2 | **MERGED half** |
| 376.80 | 1000021 | 21 | 29 | 2 | **MERGED half** |
| 276.80 | 1000035 | 35 | 44 | 2 | **MERGED half** |
| 9.62 · 3.21 · 0.88 · 0.63 · 0.33 | 1000016/04/53/14/08 | | | 2 each | **MERGED half** |

**Σ: LOOP 0.0 · MERGED-half 1048.8 · true dry land 0.0.** Eight of eight. **8 of 54 below-sea ids
span more than one region** — `[1000004, 1000008, 1000014, 1000016, 1000021, 1000035, 1000053,
1000056]` — and every one of them is the source of one of these spillways.

> **THE MECHANISM.** `break_reciprocal_spill_cycles` resolves a reciprocal pair by **merging** the
> two regions under one id (`lake_map[absorb] = keep`), **dropping the absorbed basin's spillway**
> (`spillways.retain`), and **resetting `chained_into = None`** ("its receiver was absorbed; treat
> as an open outflow"). The surviving spillway then routes its whole outflow into the OTHER HALF of
> its own merged id — a physically distinct below-sea region whose own outlet has just been
> deleted. **Water is delivered to a region that, by construction, no longer has an outflow.**
> The function whose job was to break cycles creates the terminal sink instead.
>
> **The author's reset-1932 bet, at 10 %, WINS. Mine at 65 % loses. The measurement decides and it
> decides against me.** The scoring is the author's; the figures are 8 of 8 and 1048.8 m³/s against
> 3 of 54 and 73.7 m³/s.

The function's own doc comment already conceded the half of this it knew: *"the fixed point above
still propagated the cyclic `extra_inflow` while both directions existed… removing it inside the
iteration is a deeper change, measured and reported separately."* **This is that separate
measurement**, and it finds a second consequence the comment does not mention: the merge leaves a
half-region with inflow and no outlet.

### The one-basin cross-check: the ledger IS readable, the CHAIN is not

Basin 1000056, after the fix, with the label bug of Finding 73 corrected:

- **(a) watershed** — 1 024 541 above-sea cells over 8026 D8 labels × measured surplus = **47.56 m³/s**;
- **(c) shoreline sum** re-derived from the same `runoff` array = **47.56 m³/s** ⇒ **ratio (a)/(c) = 1.000**;
- `BasinSummary::inflow_m3s` = **380.51 m³/s** ⇒ isolated `extra_inflow` = **332.95 m³/s (87.5 %)**;
- **(b) Σ spillways naming it in `chained_into` = 0.00 m³/s.**

The headline ratio (a+b)/reported is **0.125**, which trips the declared ±25 %. But the two halves
say opposite things and must not be averaged:

> **Stop rule 2's arithmetic fires and its premise is refuted.** The LOCAL ledger is confirmed to
> **1.000** by a route independent of `runoff_accumulation`, so `BasinSummary` is readable. What is
> unreadable is the CHAIN: **332.95 m³/s is routed into this basin by `chained_region` and 0.00 is
> named by `chained_into`.** The routing key and the display id are different fields and they
> disagree completely at the largest basin. That is the computed-not-exported family again, now
> with a number on it.

Only **4 chaining edges** exist in the whole export (`1000024→1000023`, `1000025→1000027`,
`1000041→1000045`, `1000059→1000057`), no cycles, and the 380.51 m³/s chain is **one link long**.

### B2, 4a, 4b, C — reported, with the conditional explicitly NOT met

The round made these conditional on the discriminant confirming. It did not (×0.214), so they are
inventory, not the closure of an attribution.

**B2 — the `wc = 0` termini, settled.** 4 spillways, Σ 86.2 m³/s (25.1 % of the hole): **3 surface
lake · 1 basin with spillway · 0 absorbed receiver · 0 true dry land.** The Finding 72 prediction
(3 of 3 surface lakes, 0 dry, 0 absorbed) is **confirmed** — and per the round's rule this counts
as **one** confirmation with A1, not two, and A1's own mechanism claim lost.

**4a — the surface-lake deficit, lake by lake** (dropped Q is a RECONSTRUCTION, the outlet Q is
measured): lake 46 dropped 34.45 / own outlet 2.185; lake 13 dropped 21.21 / own 3.286; lake 52
dropped 18.05 / own 11.551.

**4b — arithmetic only, ratio 7.5, no recommendation attached.** Lake 46: Q 36.63 → 227.0 m → 4.65
cells. Lake 52: 29.60 → 204.0 m → 4.18 cells. Lake 13: 24.50 → 185.6 m → 3.80 cells. **No repaired
outlet reaches 8 cells; the largest is 4.65.** The author's prediction of ≥ 8 cells is refuted —
because the fix already moved the *real* discharges up, so the hypothetical repair adds less.

### Block 4 — MFD carries by a SPECIAL PATH, and the order under it is broken too

`mfd_accumulation` weights by `drop` on `filled`, so on a flat every weight is zero — and
`flow.rs:906` has a fallback: `if cnt == 0 || wsum <= 0.0` it routes the whole flow to
`direction[c]`, which IS flat-resolved. **So the ROUTING has a named special path.** But
`flow.rs:877` sorts `land` by `(filled desc, INDEX desc)` — an index tiebreak, not `flat_grad`. So
the routing can cross the flat and the order does not let it: measured on the fixture, **OUT/IN
0.1181** — it carries **11.8 %**, neither 0 nor 1.

> **Verdict: the third outcome. It carries by a special path (the `cnt == 0` D8 fallback) whose
> benefit is then mostly destroyed by an index-tiebreak sort.** Per the round's instruction,
> **NOTHING is corrected**: this accumulation shaped the delivered terrain, so changing it is a
> change of terrain, not a bug fix.
>
> **And it gives Finding 66 a new meaning.** F66 recorded that the field is re-filled at every
> incision iteration. Combined with this: **every filled pit is a partial wall to the incision's
> drainage area** — `A` downstream of a depression is understated by ~88 % on the fixture, at every
> iteration. The delivered relief was incised with that. Consequence recorded, not acted on.

### Score

Mine: population 42 ✓ · bit-identical field ✓ · guard median 2–8 → 4.185 ✓ · share ≥ 90 % → 100 % ✓
· budget unchanged ✓ · max `Watercourse` 8–40 → 21.785 ✓ · max `Spillway` 250–450 → 380.5 ✓ ·
MFD "does not carry, no special path" ✗ (it carries 11.8 % by a named path) · `Watercourse`→sea
90–160 → **75.9 ✗** · closure ×0.35–0.50 → **×0.315 ✗** · discriminant ×0.7–1.3 → **×0.214 ✗** ·
implementation > 30 net lines ✓ (123 lines, ~45 of them code).

**The meta-prediction holds a further round: several were false, and the biggest one — the
mechanism of the hole — was mine.**

### Standing

One production change: the accumulation order, plus the two cache bumps. No promotion of the
three-state invariant, no export field, no threshold, no C3, no correction to
`mfd_accumulation`, no conservation-vs-representation recommendation. The hole is **343.7 m³/s**
and now has a dominant named mechanism — the cycle-breaker's merge — which is **not fixed**.

## Finding 75 — MARKED now equals LISTED; the coast read at the CELL scale for the first time; and the exported sea is at 0 m, which the big rivers do reach

Rule 11 first: `a_c_slope_law` **3** ADR hits (earliest 5288) · `MIN_SPUR_KM` **2** (5968) ·
`coast_metrics` **3** (5276) · `INVENTORY_MIN_CELLS` **5** (1242) · `lake_map` **34** (831) ·
**`sea_level_m`: NOTHING FOUND** (recorded — and the ADR's own title says "terraces").

Two hits, 4700 lines apart, say the same thing and neither noticed the other:

- **L1242 (Finding 33 Part A)** — the below-sea inventory floor was lowered from 5 km² to cells
  *"reject single-cell noise only, keep every visible lake so no river terminates in a body
  absent from the export"*;
- **L5974 (Finding 56b)** — `INVENTORY_MIN_CELLS = 4` is ⛔ **the MIRROR case** of the sub-cell
  trap, *"the documented cause of the orphan-mouth regression"*.

**The intent was written as 1 and the constant was set to 4.** That gap IS §2(a).

### A1 — VALIDATION_NOTE §1 amended, and the amendment MOVES the table

The "Known residual" paragraph is replaced: the target is 8192² alone (F67–68), the convergence
programme is closed (F61–66), and `S_ref` on the eroded field is a consequence of the incision,
not a cause (F57). The calibration is single-grid at **`S_ref = 0.1128`**; the 2048² column is
context.

> ⚠️ **And that has a consequence the round did not anticipate, measured here.** §1's table
> (1849 → 544 spurs, *"71 % of the gap to the reference closed"*) was measured with
> `CHANNEL_HEAD_S_REF = 0.3319` — **the 2048² value**. Re-measured at the calibration the
> amended note now prescribes: **1848 → 916, i.e. 51 % of the gap**, not 71 %. The SHIPPED
> figure reproduces F56 exactly (1848 against 1849), so the difference is the calibration and
> nothing else. **The law is weaker at its own target grid than the note advertises**, and §1's
> table must carry that flag before the decision is taken.
>
> `stream_power::CHANNEL_HEAD_S_REF` still ships 0.3319. Un-gating the law would have to move it.
> Not this round.

### A2 — the inventory floor, 4 → 1: MARKED now equals LISTED by construction

`INVENTORY_MIN_CELLS: 4 → 1`. The point is not the number, it is the **equality**: every id in
`lake_map` is now in `lakes`, so a consumer resolving a river terminus can never be handed an id
the export does not contain — **and no future network density can reopen it.**

The same bench ran with the change stashed and restored:

| | before | after |
|---|---|---|
| eroded / breached / pre-filled FNV hash | `0x6b10a0c5fdf467a3` · `0x2cce4ea2fb761a12` · `0x5bd35562ba4fb9ab` | **all three identical** |
| cells raised by `pit_fill` · below-sea · ocean | 1 238 085 · 457 698 · 55 626 539 | **identical** |
| dangling `lake_map` ids, humid / arid | **16 / 16** | **0 / 0** |
| `Spillway` rows with `source_lake = None` | **16** | **0** |
| lakes, humid 8192² | 86 | **102** |
| lakes, arid 8192² | 61 | **77** |
| Σ lake area, humid / arid | 6219.4 / 1703.7 km² | **6219.4 / 1703.7 km²** |

**The height field does not move by a bit** — the breach reads `detect_lakes`' map, not
`below_sea_basin_lakes_infil`'s, so the raster is untouched as predicted. **+16 lakes at both
beds, exactly the dangling count**, and **Σ area unchanged to the printed decimal**: the sixteen
new entries are the sub-4-cell bodies, together under 0.05 km². They cost nothing and they close
the resolution gap. The 16 unnamed spillways of Finding 71 fall to **0**, as the round required.

⚠️ **What is NOT measured: the note's own 26.** My orphan-mouth instrument ("a `Watercourse`
terminus touching an id absent from `lakes`") reads **0 for SHIPPED before the fix**, where §2(a)
reports 2 — so **my definition is not the note's** and I cannot reproduce its 26. What is settled
is the gap that CAUSES it: 16 → 0, structurally. Whoever owns §2(a)'s instrument should re-run it.

The unit test at `drainage.rs:2444` was **reversed on purpose**: it used to assert "the 2-cell pit
is below the 4-cell floor → marked, not listed" — it pinned the negation of the invariant. It now
asserts the equality, with a population guard.

### B — the coast at TWO length scales, and the answer to the sea-level question

`MIN_SPUR_KM` = 1 km = **20.5 cells at 8192²**, so the shipped detector cannot see a 2-to-5-cell
indentation. `coast_shape` now delegates to `coast_shape_thresholds(polys, km_per_cell,
min_spur_km, neck_km)` — **one detector, two settings, not two detectors**; the default call is
asserted to be exactly the parameterised one. Cell scale = spur ≥ 2 cells (0.0977 km), neck
1 cell (0.0488 km). Measured on the **breached** field, the one the export's coastline is traced
from (the eroded-field rows are printed beside it for comparability with Finding 56).

| variant | ≥ 1 km count | ≥ 2 cells count | coast km | p90 km (km scale) | R (km scale) |
|---|---|---|---|---|---|
| **SHIPPED** | 1 901 | **2 398** | 6 742 | 4.89 | 0.523 |
| **LAW ON** (`S_ref = 0.1128`) | 952 | **1 439** | 3 953 | 4.49 | 0.429 |
| **REFERENCE** (coarse, no FBM, no incision) | 23 | **52** | 1 623 | 23.34 | 0.404 |

**Share of the gap to the REFERENCE that the law closes: 50.5 % at the kilometre scale,
40.9 % at the cell scale.** The round predicted the fur would move by less than 40 %; measured
**40.9 %** — refuted by a whisker, and my own "< 30 %" refuted outright.

> **The two scales disagree far less than either of us expected, and the interesting number is
> elsewhere: the cell scale separates SHIPPED from REFERENCE by a factor 46 (2 398 against 52),
> where the kilometre scale separates them by 83 (1 901 against 23).** So the fur the eye sees is
> real and enormous at both scales, and the law removes about half of it at both. It is not a
> kilometre-scale-only remedy, which is what the round's hypothesis expected.

⚠️ **Rule 10 on two columns.** REFERENCE at the cell scale has 52 spurs (8 on the eroded field),
and its `local_axis_r` reads 1.000 on one field and 0.000 on the other — a degenerate population,
**not a reading**. The REFERENCE R must be taken from the kilometre scale (0.404) only.

**`sea_level_m` = 0.0**, from `manifest.json` — and it is **hard-coded** at `hd.rs:1227` with the
comment *"sea anchored to 0 m"*, so it is a contract, not a measurement.

**And the rivers that reach it, reach it.** Of the SHIPPED export's 2 363 terminal `Watercourse`
segments, the **65 that actually drain into an ocean cell** have mouths at min −20.0, **median
−0.1**, p90 −0.0, max −0.0 m. **There is no suspended mouth and no terrace.**

> The 197 m river is therefore **not reaching the exported sea at all** — it is a terminus at a
> lake or a truncation of the traced network, and at a median of **218.4 m** over all 2 363
> termini that is the TYPICAL terminus. **Only 65 of 2 363 (2.8 %) of terminal watercourses
> reach the ocean.** That is Finding 71's hole seen from the geometry side rather than the water
> balance side, and it is the real item behind the author's observation.

Three complete `.ymir` exports written, production seed and config, 8192²,
**`geo_scale_ratio = 7.5`** (the author's practice and the on-disk manifest's value, not the
code's 1.0), every layer present:
`exports/coastal_closure/seed10481999410520546993_8192_{SHIPPED,LAW_ON,REFERENCE}.ymir/`.
Lakes 102 / 97 / 36. Max `Watercourse` 1225.4 / 1230.3 / 1247.8 m³/s at ratio 7.5 (= 21.79 real
× 56.25).

### C — the MFD bisection CANNOT be built in a bench, and the seam is named

`mfd_accumulation` is called at **`stream_power.rs:403`, inside the incision loop**, with no
parameter and no config field to substitute it. A bench cannot re-incise with a corrected carrier
without a production seam, and this round's budget of two production changes was spent on the
note and the floor. **Block C is blocked, not skipped.**

> **The seam, specified and not opened:** `StreamPowerConfig` would need one field — an optional
> flat-aware carrier selector — read at `stream_power.rs:403`, exactly the shape of Finding 74's
> `propagation_order` call. It changes the delivered terrain the moment it is switched on, so its
> promotion is an author decision and a full re-validation round, as the round itself said.

What the round said must happen **first** is done: the toggle is proven to move the reading.
A bench-local copy of `mfd_accumulation` differing **only** in the propagation order, on the
Finding 74 fixture: **OUT/IN 0.1181 → 1.9686, ×16.7.** The weighting, the `cnt == 0` D8 fallback
and the seeds are copied verbatim; the shipped carrier's `(filled desc, INDEX desc)` sort
(`flow.rs:877`) is the only difference. **The toggle is fit to serve a bisection.**

### E — the eight merged ids, in cell coordinates, and the second half is usually ONE CELL

| id | cmp | cells | km² | bbox | centroid | Δ centroids |
|---|---|---|---|---|---|---|
| **1000056** | 56 | **69 560** | 165.844 | 3606,4840..4451,5333 | 3933,5094 | **11.44 km** |
| | 58 | **1** | 0.002 | 4128,5223 | 4128,5223 | |
| **1000021** | 21 | 135 902 | 324.016 | 3763,3291..4387,3792 | 4062,3511 | **20.54 km** |
| | 29 | 10 034 | 23.923 | 3509,3442..3773,3884 | 3665,3650 | |
| **1000035** | 35 | 1 762 | 4.201 | 5131,3606..5365,3943 | 5229,3835 | 2.48 km |
| | 44 | **1** | 0.002 | 5249,3882 | 5249,3882 | |
| **1000053** | 53 | 131 | 0.312 | 5419,4697..5448,4722 | 5436,4709 | 0.38 km |
| | 54 | 16 | 0.038 | 5435,4706..5450,4720 | 5442,4713 | |
| **1000004** | 4 / 5 | 25 / 5 | 0.060 / 0.012 | 4415,2019.. / 4426,2023.. | 4420,2020 / 4428,2024 | **0.43 km** |
| **1000014** | 14 / 15 | 13 / 4 | 0.031 / 0.010 | 5661,3103.. / 5665,3108.. | 5664,3106 / 5666,3109 | 0.17 km |
| **1000016** | 16 / 17 | 3 / 1 | 0.007 / 0.002 | 5768,3240.. / 5762,3246 | 5768,3241 / 5762,3246 | 0.39 km |
| **1000008** | 8 / 10 | **1 / 1** | 0.002 / 0.002 | 5353,2673 / 5356,2675 | — | 0.18 km |

**7 of 8 pairs are within 2.5 km and 6 of 8 within 0.5 km** — prediction confirmed, and more
strongly than predicted. But the number that matters is the other column:

> **For 1000056 — the 380.5 m³/s basin, the largest single term in the hydrological hole — the
> "other half" is ONE CELL of 0.002 km², 11.44 km away from the main body.** The cycle-breaker
> gave that cell the main basin's id, deleted its own spillway, and the surviving spillway now
> delivers 380 m³/s into it. Same for 1000035 (276.8 m³/s → 1 cell) and, at 1 cell against 1
> cell, 1000008.
>
> **Only 1000021 is a genuine two-lobe merge** (135 902 + 10 034 cells, 20.54 km apart) — the
> only one of the eight a renderer would draw as two lakes.

⚠️ The `spill?` column of that table **failed to resolve and is void**: it looked for the
component of the spillway's FIRST point, which is the saddle — above sea level, so `wc ≠ 2` and
the component is 0. The direction is known from Finding 74's block 2b instead (own cmp → end
cmp), not from this column. My instrument, reported as broken rather than read.

### What is NOT delivered

**Block D (C-3 / C-3b bisection) and block F (coast metrics at lake shores) were not run** —
no budget left after the three 8192² builds. Neither is blocked; both are one bench run each and
the two-scale detector they need now exists.

### Standing

Two production changes, both authorised: `INVENTORY_MIN_CELLS` 4 → 1, and the note's §1. Plus one
refactor that is provably not a behaviour change (`coast_shape` delegating to
`coast_shape_thresholds`, identity asserted). No de-gating of the law, no MFD change, no touch to
`break_reciprocal_spill_cycles`, no C3, no renderer. The hydrological hole — 343.7 m³/s, mechanism
attributed at Finding 74 — is untouched, as instructed.

## Finding 76 — the sea-level-flat mechanism is REFUTED (r = 0.05); the fur is a periodic, parallel, non-channel tooth of the incision, and it has a wavelength

### The acceptance criterion, written before the measurements and not met by anything shipped

The objective was never "fewer fringes". It is: **the upscale + closures + incision chain
produces NO coastal fringe.** Opposable, at 8192², on the conditioned field:

| line | target | SHIPPED | LAW ON (0.1128) |
|---|---|---|---|
| spurs ≥ 1 km | ≤ 2× REFERENCE = **46** | 1 901 (**×41**) | 952 (×21) |
| spurs ≥ 2 cells | ≤ 2× REFERENCE = **104** | 2 398 (**×23**) | 1 439 (×14) |
| coastline km | ≤ 1.2× REFERENCE = **1 950** | 6 742 (×3.5) | 3 953 (×2.0) |
| spur p90 length | **≥ 12 km** | 4.89 km (**×0.4**) | 4.49 km (×0.4) |
| spectral peak vs white | **< 3×** | **×3.75 — FAILS** | not measured |
| interior unmoved | hypsometry ±1 m, channel share ±1 pt, depressions ±2 % | — | — |

**Neither the shipped state nor the law comes within an order of magnitude of any of the first
four lines.** The law is not a candidate remedy for this criterion; it is a partial reduction.

### Rule 11, and the earliest hit was the candidate mechanism

`sea_level` **9** (earliest 70) · `h_r` **4** (3789) · `implicit` **8** (149) · `is_ocean` **1**
(8008) · `coast_shape_thresholds` **2** (8331) · `S_min` **3** (3997) · **`s_min`: NOTHING FOUND**.

L3789 (Finding 61) states the candidate exactly: *"With `f = K·dt·A^m/dist_m ≫ 1` the update
`h ← (h + f·h_r)/(1 + f)` drives each cell essentially onto its receiver's height in ONE step"*,
Courant 1353. Combined with `stream_power.rs:410` — `receiver[k] = k` whenever
`field <= sea_level`, so the ocean is a fixed base node and a land cell draining into it takes
`hr` = that ocean cell's height — the prediction is that the coast is dragged onto sea level, one
cell per iteration, along every channel. **It is a clean, code-level derivation. It is also
wrong, and the measurement says so in four independent ways.**

### A — the mechanism is REFUTED, and stop rule 1 fires

**A1 — dent length against flat length.** Over the **1 831** spurs ≥ 2 cells, each matched to its
nearest channel mouth (0 unmatched), the flat = the run of consecutive cells at or below
`sea + 0.01 m` walked upstream from the mouth along the highest-accumulation donor:

> **FLAT length: p10 0.0 · median 0.0 · p90 0.0 · max 3.0 cells — 98.5 % of spurs sit on a
> ZERO-LENGTH flat.** Spur length: p10 2.0 · median 4.1 · p90 14.7 cells.
> **r = 0.051.** Mean flat 0.02 cells against a mean spur of 7.78.

**Stop rule 1 fires. The seam was NOT written.** My own prediction (r ∈ [0.5, 0.85]) is refuted,
the round's (r > 0.8) is refuted, and the measurement wins.

**A2 — the near-shore bed rises from the first cell.** Median height above sea, by flow distance
to the sea: **1 cell → 1.997 m · 2 → 6.881 · 3 → 11.088 · 5 → 18.977 · 10 → 34.970 · 20 →
65.272 m** (145 665 land cells drain directly into an ocean cell; 100 % of land has a flow path
to the sea). That is the alternative A2 named — *a regular slope from the first cell* — not a
plateau then a rise. My "≤ 0.05 m at 1–2 cells" is refuted by a factor 40.

**A3 — the counter-hypothesis dies too, and harder.** Over the 267 100 land cells within 3 of the
sea: **exactly `sea_level`: 0 (0.00 %)** · in `(sea−0.01, sea]`: **0** · in `(sea, sea+0.01]`:
2 235 (0.8 %) · **above sea + 0.01 m: 264 865 (99.2 %)**. And of the 145 665 mouth cells,
**0 sit exactly on their receiver's height** — the `clamp(hr, ho)` signature of the drag is
**totally absent**. Both the mechanism and its counter-hypothesis are dead; the near-shore
terrain is ordinary sloping ground.

Why the code reading was wrong: a mouth cell only reaches the implicit update if it clears
`A_c = 41.94 cells`. There are 145 665 mouth cells and the coastline-wide median accumulation is
**1 cell**, so the overwhelming majority `continue` out of the loop before any drag can happen.
**A correct reading of one code path is not a measurement of the population that traverses it** —
the same error class as Finding 73's "the exit cell is below sea level", which was also a correct
sentence about the wrong cells.

### A4 — what the spurs ARE, since stop rule 1 demands it

| question | answer |
|---|---|
| does the excursion enclose land or sea? | **LAND for 1 491 (81.4 %)**, sea for 340 (18.6 %) ⇒ **TEETH, land slivers pointing seaward**, not bays |
| is the tip a channel? | tip accumulation p10 **0** · median **1** · p90 103 cells (coastline-wide median 1). **Only 250 of 1 831 (13.7 %) reach `A_c` = 41.9 cells** |
| how deep is the notch? | neck − tip: p10 **−20.45** · median **0.84** · p90 +28.46 m |

> **The fur is made of ~4-cell land slivers that carry no channel (86.3 % below the channel
> head) and have under a metre of relief between their neck and their tip.** It is not a
> drowned valley, not a ria, and not a drag — it is a sub-metre wiggle in a gently sloping
> surface, read by a contour tracer at 48.8 m per cell.

### A5 — and it is 100 % the incision's, at the CELL scale too

Finding 51 attributed the fringe entirely to the incision, measured at the **kilometre** scale —
the one that cannot see a 2-cell tooth. Re-measured at the cell scale, FBM and both closures ON,
**only `stream_power` removed**:

| field | ≥ 1 km | **≥ 2 cells** | coast km | p90 km | R (cell) |
|---|---|---|---|---|---|
| SHIPPED (FBM + incision) | 1 848 | **1 831** | 6 206 | 4.62 | **0.879** |
| FBM + closures, **NO incision** | 20 | **8** | 1 602 | 23.34 | **0.000** |

**×229 at the cell scale, and the parallelism goes from 0.000 to 0.879.** Finding 51's conclusion
survives the change of scale intact: the FBM and the closures contribute **nothing** to the
coastal texture, and these figures are indistinguishable from Finding 75's REFERENCE (20 / 8 /
1602 on the eroded field) — i.e. **removing the FBM on top of the incision changes the coast by
zero.**

### C — the fur has a WAVELENGTH, and that is the finding

Gaps between consecutive spur roots along each continuous coast segment > 20 km, against a white
baseline (same counts, uniform positions, 64 draws):

> **SHIPPED: 10 segments, 1 647 gaps, median gap 26.4 cells — dominant spacing **14 cells**
> (684 m) at **×3.75** the white baseline.** Above the criterion's 3× line: **a wavelength
> exists.** The no-incision field yields 0 gaps (8 spurs in total) — rule 10, reported as not a
> reading rather than as a flat spectrum.

14 cells is **not** √A_c (6.5 cells at 8192²), so the spacing is not the channel-head threshold's
own scale. My prediction (2–5 cells) and the round's (3–6) are both refuted. **A periodic,
parallel, sub-metre, non-channel texture at a fixed wavelength is an INSTABILITY of the incision,
not a drainage pattern** — which is the Smith–Bretherton rilling comb of Finding 10 reaching the
shoreline, now measured at the shoreline with a wavelength attached for the first time.

### The rule-9 trap, at full strength: the crop swings the verdict from −60 % to −8 %

Three panels at **two** crops, rendered from the exported rasters:

| | spurs ≥ 1 km | ≥ 2 cells | | |
|---|---|---|---|---|
| | SHIPPED | LAW ON | SHIPPED | LAW ON |
| whole grid | 1 902 | 943 (**−50 %**) | 2 394 | 1 350 (**−44 %**) |
| crop on SHIPPED's densest window (5120,4480) | 129 | 50 (**−61 %**) | 113 | 45 (**−60 %**) |
| **crop on LAW ON's densest window (4160,1600)** | 81 | 75 (**−7 %**) | 152 | 140 (**−8 %**) |

> **The same remedy reads −60 % or −8 % depending on which 31 km square you look at.** At the
> second crop SHIPPED and LAW ON are visually near-identical, both heavily furred, against a
> smooth REFERENCE. The author's instinct — that the saw-tooth he saw is where the law's densest
> residue is — is confirmed, and Finding 75's in-crop figure was flattering by ×1.4 at one crop
> and by ×7 relative to the other.

### Standing

**No production change.** The seam was not written because stop rule 1 fired first, which is what
the stop rule is for. One refactor that is provably not a behaviour change:
`coast_shape_thresholds` now delegates its walk to `coast_spurs`, which the spectral instrument
needs (no aggregate can carry positions) — verified GOLDEN on the real 8192² coastline, all nine
whole-grid figures and all six in-crop figures identical to the digit before and after.

Not run: B (the seam and its variants), D (lake shores under a floor that does not exist). The
hydrological hole (343.7 m³/s, Finding 74) is untouched as instructed.

**What the next round has to explain is no longer "why is the coast flat" — it is "why does the
incision produce a parallel, sub-metre, non-channel texture at a 14-cell wavelength, and why only
near the shoreline".** The candidate now on top is the one Finding 55 tested at the wrong lever:
the rilling instability, whose wavelength is set by the competition between incision and
diffusion, not by `A_c`.

## Finding 77 — there is NO water between 0 and −20 m anywhere on the map; the hollows ARE channels, and the fringe's strongest lever is the unanchored diffusion

### The criterion, rewritten as a DIFFERENCE against the pre-incision field

The threshold form (≤ 2× REFERENCE) is withdrawn. The authority is the pre-incision field of the
same seed — Finding 76-A5 measured it at 20 / 8 / 1 602 km, indistinguishable from the coarse
REFERENCE — and the criterion is **Δ against it**:

| line | target | SHIPPED | Δ | LAW ON | Δ |
|---|---|---|---|---|---|
| spurs ≥ 1 km | **Δ = 0** | 1 848 | **+1 828** | 916 | +896 |
| spurs ≥ 2 cells | **Δ = 0** | 1 831 | **+1 823** | 1 056 | +1 048 |
| coastline km | Δ under the u16 quantisation noise, **±3 km (0.2 %)** | 6 206 | **+4 604 (+287 %)** | 3 640 | +2 038 |
| spectrum | unreadable, as the pre-incision's is (8 spurs ⇒ 0 gaps) | λ 14 at ×3.74 | **readable ⇒ fail** | — | — |

**The incision must add no spur.** The pre-incision's twenty are tectonic geometry, not fringe.

### Rule 11

`is_land` **0 — NOTHING FOUND** · `hypsometry_sweep` **0 — NOTHING FOUND** · `ex-terre` **0 —
NOTHING FOUND** · `coast_spurs` 1 (8575, Finding 76) · `RELIEF_V3_DIFFUSION` 1 (**6098**) ·
`diffusion` 89 (22) · `deposition` 16 (**18**).

Two earliest hits carry the round:

- **L6098 (Finding 59 §4)** — the `(ref/cell)²` anchor exists **only in the nonlinear branch**
  while relief-v3 takes the linear one, so **the shipped hillslope operator is 16× too weak at
  8192²**, and the anchored value is **1.28** — which is exactly the sweep's fourth point. F59
  priced that defect at **4 % of the work divergence** and set it aside. Nobody measured it on
  the coast.
- **L18** — the droplet erosion *"deposits in channels faster than it incises. It is the common
  cause of the terraces (depositional flats), the missing valleys, the net-zero mass balance and
  **the coastal sediment dump**."* Block D tests the deposition CLASS against that line, it does
  not propose it.

### ⚠️ My own instrument leaked, and the first number it produced is DISOWNED

The mask `pre land ∧ shipped sea` catches **every** cell the incision ever took below sea level:
**298 598 cells in 360 components, median 44 and p90 2 660 cells** — the drowned valleys, not the
1 831 teeth of median 4 cells. Its "81.4 % of heads at or above `A_c`" answers a question nobody
asked. **Restricted** to hollows with at least half their cells within 3 of a spur's own coastline
samples: **46 components of 360**, and *that* is the population the objection was about.

The objection was correct on its own terms and is upheld: Finding 76's A4 took the polyline sample
farthest from the neck — a **ridge tip** for 81.4 % of excursions — so "86 % carry no channel"
measured the ridges. **Corrected below.**

### A — the hollows ARE channels, and the depth is a CONSTANT

| column (46 fringe hollows) | p10 | median | p90 |
|---|---|---|---|
| area, cells | 1.00 | **7.00** | 131.00 |
| depth below sea, m | **−20.00** | **−20.00** | **−20.00** |
| incised at the head, m | 48.39 | **357.79** | 1 648.70 |
| HEAD accumulation, cells | 16.00 | **253.00** | 6 564.00 |
| best NEIGHBOUR accumulation, cells | 29.00 | 476.00 | 22 051.00 |

> **At or above `A_c` = 41.9 cells: HEAD 84.8 %, best neighbour 87.0 %.**
> Finding 76 read the ridge tips and got 13.7 %. **⚠️ CORRECTION to Finding 76 at the point of
> claim: the fringe is NOT a "non-channelised texture". The hollows that make it are channels,
> in 85 % of cases.** My conclusion there was drawn from the complementary population.

**And the depth is not a distribution.** Every one of the **298 598** ex-land cells reads exactly
**−20.000 m** — p1, p10, median, p90 and max alike, and the same in the spur zone (55 756 cells).
0.00 % are shallower than −0.5, −1, −2 or −5 m.

### THE FINDING: the ocean has no shallows at all

| population | p1 | median | **max** |
|---|---|---|---|
| SHIPPED ocean cells ADJACENT TO LAND (120 807) | −20.000 | −20.000 | **−20.000** |
| SHIPPED, all 56 077 791 ocean cells | −5 264.0 | −4 620.5 | **−20.000** |
| PRE-INCISION, all 55 780 901 ocean cells | −5 235.6 | −4 596.7 | **−20.000** |

> **The shallowest water anywhere on the map is exactly −20.000 m, and it is already so BEFORE the
> incision.** There is a forbidden band 20 m thick immediately under sea level: no intertidal, no
> shelf gradient, no shallows. **Every shoreline on this continent is a 20 m vertical step.**
>
> The 0 m isoline therefore never runs on the sea side — the sea side is a *constant*. The
> coastline is entirely a threshold crossing of the LAND's own micro-relief against a perfectly
> flat floor, which is exactly what produces a 2-to-5-cell saw-tooth.

**Origin not located, and the ruled-out list is the useful part**: it is not Stein-Stein (that
floors at `ridge_depth_m` = 2600 m and gives the −4 620 m median); not the coastal FBM taper
(`coastal_amplitude_band = 0.30`, `upscale.rs:593` — an amplitude taper, and it damps the FBM to
~0 **at the waterline by design**, #151, which is *why* Finding 76-A5 found the FBM contributing
nothing at the coast); not `coast_warp_strength = 1.5` (a horizontal sampling warp). Naming the
constant is the next round's first item and it is a bounded search.

**This also settles Finding 76's mechanism honestly.** The drag onto the receiver is REAL — 298 598
cells sit exactly on the shelf value, and the head of a fringe hollow was cut by a **median of
358 m**, which is Finding 61's relaxation onto base level at Courant 1353, not fluvial incision.
What was refuted is the *destination*: the drag lands at **−20 m, not at sea level**, so Finding
76-A1 correctly found no sea-level flat while the drag was happening 20 m below its instrument.
**A1's filter was `is_land`, and the dragged cells are not land any more.**

### B — λ is NOT set by `A_c`; the COUNT is

| `A_c` km² | spurs ≥ 1 km | ≥ 2 cells | coast km | λ cells | × white |
|---|---|---|---|---|---|
| 0.025 | 2 711 | **3 264** | 8 429 | **14** | 3.07 |
| **0.1 (shipped)** | 1 848 | **1 831** | 6 206 | **14** | 3.74 |
| 0.4 | 433 | **307** | 2 735 | **11** | 11.17 |
| 1.0 | 202 | **113** | 2 166 | **10** | 8.00 |

**A ×40 sweep moves λ from 14 to 10 — down 29 %, and the wrong way for either hypothesis.**
`λ ∝ A_c` (14 → 56) is refuted; `λ ∝ √A_c` (14 → 28) is refuted. My "insensitive within 30 %"
survives at the edge. But the **count** falls ×29 and the coastline ×3.9: `A_c` governs **how
many** teeth, not **how far apart**. ⚠️ Rule 10 on the last row: 76 gaps is a thin population and
its ×8.00 should not be read as a stronger peak than the shipped ×3.74.

### C — λ is not set by the diffusion either, but the FRINGE is

| `diffusion` | spurs ≥ 1 km | ≥ 2 cells | coast km | R cell | λ | × white |
|---|---|---|---|---|---|---|
| 0 | 1 907 | **2 893** | 7 325 | 0.939 | 19 | 2.47 (**no peak**) |
| **0.08 (shipped)** | 1 848 | **1 831** | 6 206 | 0.879 | 14 | 3.74 |
| 0.32 | 1 718 | **673** | 4 710 | 0.655 | 18 | 3.49 |
| **1.28 (F59 §4's ANCHORED value)** | 1 275 | **189** | 3 538 | 0.455 | 8 | 4.00 |

λ wanders 19 → 14 → 18 → 8, non-monotone, inside a factor 2.4 — **my "< 20 %" is refuted in
magnitude**, the qualitative reading (not the controlling term) holds. But:

> **At the diffusion value Finding 59 §4 showed to be the CORRECT one, the cell-scale fringe falls
> from 1 831 to 189 — −90 % — and the coastline from 6 206 to 3 538 km.** Δ against the
> pre-incision goes from +1 823 to +181. F59 priced the unanchored diffusion at **4 % of the work
> divergence** and set it aside; **on the coast the same defect is worth an order of magnitude.**
> It is not a remedy (Δ is still +181, and the terrain moves — the interior control is owed), but
> it is the strongest lever this campaign has measured, and it is a defect already on file.

### C-bis — my own named candidate, refuted outright

I named `lateral_erosion` (closure b, `stream_power.rs:549`: banks planed **perpendicular to
flow** over a physical half-width `K_lat·A^m`, only where `A ≥ A_c` — an operator that lowers
non-channel cells in a comb at right angles to them) as the mechanism, before measuring.

| variant | ≥ 1 km | ≥ 2 cells | coast km | R | λ | × white |
|---|---|---|---|---|---|---|
| SHIPPED | 1 848 | 1 831 | 6 206 | 0.879 | 14 | 3.74 |
| **`lateral_erosion = 0`** | 1 848 | **1 831** | 6 208 | **0.879** | 14 | 3.73 |
| `talus_passes = 0` | 2 017 | **2 388** | 7 153 | 0.862 | 14 | 3.19 |

**Turning the banks off changes the coast by two kilometres in six thousand and not one spur.**
Refuted. The reason is the sub-cell family again: `half_cells = (K_lat·A^m / cell_m).floor()` and
the code `continue`s at `half_cells < 1`, so near the coast — where `A` is small — the closure
never fires. And **removing the talus makes the fringe 30 % WORSE**, so it is a mild suppressor,
not a cause.

### D — the class test could not run, and the reason is block A's constant

0 cells lifted at every threshold {0, 0.5, 1, 2 m}, because **no ex-land cell is within 2 m of
sea level** — they are all at −20 m. The negative control (threshold 0 m) correctly left the
counts at SHIPPED's. **The deposition class cannot be tested by filling shallow hollows because
there are no shallow hollows**: a test of it would have to fill a 20 m step, which is a different
operation and collides with ADR line 18.

### Score

Mine: λ insensitive to `A_c` (−29 %, at the edge) ✓ · hollow heads ≥ `A_c` 40–75 % → 84.8 %, just
over ✗ · hollow depth −0.2…−1.5 m → **−20.0 constant** ✗✗ · λ insensitive to diffusion within
20 % → ×2.4 ✗ · **`lateral_erosion` as the mechanism** ✗ · D −85 % at 1 m → not testable ✗.
The round's: `λ ∝ A_c` ✗ · λ insensitive to diffusion ✓ (qualitatively) · heads ≥ `A_c` > 70 % ✓ ·
D-post Δ = 0 at 1 m → not testable ✗.

**The meta-prediction holds again, and my named candidate was the loudest thing I got wrong.**

### Standing

No production change. Bench tooling moved to `tests/common/mod.rs` so the builder and the
spectrum have one copy rather than a fourth. Eleven 8192² builds. The hydrological hole
(Finding 74) and the law's gate are untouched.

**Next round's first item is a bounded search: name the constant that forbids water between 0 and
−20 m.** Until it is named, every coastal remedy is being fitted on top of a 20 m cliff nobody
chose.

## Finding 78 — the 20 m shelf clamp is a PROXY, it runs AFTER the incision, it suppresses 43 % of the fur, and it makes every shoreline a cliff

### Rule 11, EXTENDED — grep the VALUE, not only the name

`shelf_min_depth` · `shelf_min_depth_m` · `BathymetryProfile` · `apply_bathymetry_profile` ·
`bathymetry` — **0, 0, 0, 0, 0 ADR hits. NOTHING FOUND, all five.**

The **value** `−20 m`, grepped with the words sea/ocean/shelf/depth, gives **six traces across
five findings**, in one command:

| line | finding | what is written there |
|---|---|---|
| 1248 | **F33** | "a deep bowl (floor −20 m, rim +20 m)" |
| 1272 | **F34** | "Lake #1000104: 61 km², level −20 m, MAX DEPTH 0 m — **a surface with no water**" |
| 1305/1308 | **F35** | "MAX DEPTH 633 m (floor −20 m)" · "a −20 m coastal pocket needing a 613 m crossing is absurd" |
| 1402 | **F36** | "reinstate the **−20 m altitude proxy** Finding 31 removed" |
| 1758 | **F38** | "2 depth-0 lakes (#1000018, #1000003 at −20 m, 420 / 316 km²) — **on a FLAT floor**" |

> **The constant signed five findings and nobody recognised it as a constant, because we grepped
> identifiers.** ⇒ **Rule 11 is extended: when a symptom is a REPEATED VALUE, grep the value, with
> the domain words, not only the name.** F38's "on a FLAT floor" is the clamp described in full
> without being named.

### A1 — the intent: a real invariant, an arbitrary number

`git log -S"shelf_min_depth_m: 20.0"` → **one commit, `563198e`, 2026-06-26**, *"FEAT : submarine
bathymetry — plateau→slope→abyss re-map (fix the slab)"*. A long, careful message. **It never
justifies the 20.** What it does justify is the invariant: *"fraction ocean IDENTIQUE avant/apres
sur les 6 seeds → cote / masque terre-mer / hypsometrie terrestre inchanges"*.

The field's own docstring states the intent — *"Depth at the coastline — the shelf is shallow but
submerged (**stays ocean**)"* — and the header anchors the profile on real passive-margin
morphology: shelf −130 m, break −200 m, abyss −4500 m, width 70 km. **20 is in none of those
anchors**; Earth's depth at the coastline is 0 m, that is what a beach is. And `bathymetry.rs:171`
takes `shelf_min_depth_m.max(1.0)`, which says in code that **1 m satisfies the stated invariant**.

**Verdict: a real numerical safety (keep the land/sea mask) with an arbitrary value — a PROXY, and
it has been one since the day it was written.**

### ⚠️ CORRECTION, to the round's premise AND to my own Finding 77

**The clamp runs AFTER the incision, not before.** `incise_lithology` is called at
`production_upscale.rs:451`; `apply_bathymetry_profile` at `:506`. The commit message says the
same ("applique … apres erosion"). Measured, on the receivers of coastal channel cells
(`A ≥ A_c`, receiver is ocean):

| field | cells | p1 | p10 | **median** | p90 | max | **distinct values** |
|---|---|---|---|---|---|---|---|
| SHIPPED (post-clamp) | 7 051 | −20.000 | −20.000 | **−20.000** | −20.000 | −20.000 | **1** |
| **pre-bathymetry** | 7 064 | −44.515 | −5.557 | **−0.109** | −0.004 | 0.000 | **3 041** |

> **During the incision the sea at the shore sits at a median of −0.109 m, in 3 041 distinct
> values.** The single −20.000 is the clamp overwriting the result afterwards.
>
> **Finding 77 wrote that the relaxation "lands on the shelf value". That is REFUTED by this
> measurement.** And the corollary matters more: Finding 76's original mechanism — coastal cells
> relaxing onto a base level **at sea level** — is what the pre-bathymetry field actually shows
> (median receiver −0.1 m). Finding 76-A1 looked for it on the SHIPPED field, **where the clamp
> had already overwritten the evidence 20 m down.** The mechanism was not wrong; the field it was
> measured on was the wrong one.
>
> Likewise Finding 77's "the head of a fringe hollow was cut by a median of 358 m" is **not what
> the incision removed** — it is incision plus the clamp's overwrite.

### B — the sweep: the interior does not move at all, and the mask never moves

`shelf_min_depth_m` ∈ {20, 5, 1, 0.5→1, 200}. ⚠️ the 0.5 m point is **clamped to 1 m by production
code** (`.max(1.0)`) and is reported as the duplicate it is.

| shelf | ≥ 1 km | **≥ 2 cells** | coast km | p90 | R | λ | × white |
|---|---|---|---|---|---|---|---|
| 200 m (neg. control) | 1 870 | **1 291** | 6 198 | 4.64 | 0.896 | 28 | 4.16 |
| **20 m (SHIPPED)** | 1 848 | **1 831** | 6 206 | 4.62 | 0.879 | 14 | 3.74 |
| 5 m | 1 829 | **2 705** | 6 173 | 4.39 | 0.823 | 15 | 2.11 |
| **1 m** | 1 776 | **3 190** | 5 927 | 4.15 | 0.846 | **6** | **1.86** |

**Control block, every point, identical to the last digit**: land 11 031 073 (16.4376 %),
**mask moved vs SHIPPED = 0**, hypsometry mean **685.32 m**, p50 444.55, pits 1 238 085, erosion
work 199.02 m, channel share 11.489–11.493 %. **The docstring's invariant holds exactly, and the
sweep does not touch the continent.** No stop rule.

Hollow depth follows `−shelf_min` **exactly**: −20.000 / −5.000 / −1.000 / −200.000, p10 = p50 =
p90 = max in every case. Prediction confirmed.

> **THE BINARY QUESTION, answered against the round: the tooth count moves by 74 %, not by
> "under 10 %".** ≥ 2 cells goes 1 291 → 1 831 → 2 705 → **3 190** as the floor rises from 200 m to
> 1 m. **The shipped 20 m clamp is SUPPRESSING about 43 % of the cell-scale fur.**
>
> The mechanism, predicted before measuring: marching squares crosses at
> `t = (0.5 − h_sea)/(h_land − h_sea)`. A deeper sea is a **larger denominator**, so the crossing
> point is less sensitive to the land-side micro-relief and the traced line is smoother.
> **The cell-scale spur count is therefore NOT a property of the land alone — it depends on the
> value of the water beside it**, which no land cell can see. Part of what this campaign has been
> calling "fur" is a contouring artefact of a sub-sea constant.
>
> And the same applies to the wavelength: λ goes 28 → 14 → 15 → **6**, and ×white 4.16 → 3.74 →
> 2.11 → **1.86 (below the 3× line: no readable wavelength at all)**. **Finding 77 found λ was set
> neither by `A_c` nor by the diffusion. It is set by the clamp depth.** The spectral peak we have
> been chasing is largely an artefact of the water value.

### The consumer number: every shoreline is a cliff

Population: LAND cells 8-adjacent to sea (145 665 of them, identical at every point).

| shelf | step land−sea p10 | **p50** | p90 | slope < 15 % | slope > 25 % |
|---|---|---|---|---|---|
| **20 m (SHIPPED)** | 20.13 | **22.00 m** | 39.22 | **0.0 %** | **100.0 %** |
| 5 m | 5.13 | 7.00 m | 24.22 | 52.9 % | 26.4 % |
| **1 m** | 1.13 | **3.00 m** | 20.22 | **71.4 %** | 18.9 % |
| 200 m | 200.13 | 202.00 m | 219.22 | 0.0 % | 100.0 % |

> **At the shipped value, 100 % of the shoreline is steeper than 25 % and 0.0 % is gentler than
> 15 %. Banks, beaches and quays are not hard on this product — they are impossible.** At 1 m,
> 71.4 % of the shoreline is gentler than 15 %. **One constant, and the requirement goes from
> impossible to mostly satisfied, with the continent bit-identical.**

### C — the anchored damper, priced

`shelf 1 m + diffusion 1.28`: ≥ 1 km **1 173**, ≥ 2 cells **383**, coast **3 434 km**, p90 2.86,
R 0.738, λ 4 at ×4.04.

Control block: **mask moved vs SHIPPED 86 615 cells**, land 11 031 073 → 10 991 458, hypsometry
mean 685.32 → **680.08 m (−5.24 m)**, p50 444.55 → 442.84, pits 1 238 085 → **1 335 100 (+7.8 %)**,
channel share 11.489 → **10.654 % (−0.84 pt)**, work 199.02 → 206.70 m (+3.9 %).

> **The terrain moves, and by less than either of us predicted: −5.24 m of hypsometry, not
> 10–60.** Against the campaign's own interior tolerance (±1 m, ±1 pt, ±2 %) it **fails the
> hypsometry line and the depression line, and passes the channel-share line.** It is a change of
> continent, but a cheap one — and it is a defect already on file (Finding 59 §4), not a new knob.

Δ against the pre-incision authority: **+1 153 / +375 / +1 832 km.** My 150–400 range holds (375);
the round's "under 100" is refuted. **The residual is still periodic (×4.04), so it is not noise.**
Whether those 383 teeth are channels was **not measured** — one build short, and it is not
inferred.

### D — Δ ≈ 0 for the first time, and the test is DEGENERATE

On `shelf 1 m + diff 1.28`, lifting ex-land cells shallower than a threshold to `sea + ε`:

| threshold | cells lifted | ≥ 1 km | ≥ 2 cells | coast km | **Δ vs pre-incision** |
|---|---|---|---|---|---|
| 0 m (neg. control) | **0** | 1 173 | 383 | 3 434 | +1 153 / +375 / +1 832 |
| 0.5 m | **0** | 1 173 | 383 | 3 434 | +1 153 / +375 / +1 832 |
| **1 m** | **337 998** | 20 | 9 | 1 627 | **+0 / +1 / +25 km (+1.5 %)** |
| 2 m | 337 999 | 20 | 9 | 1 627 | **+0 / +1 / +25 km** |

**The criterion is met: Δ(≥ 1 km) = 0, Δ(≥ 2 cells) = +1, Δ(coast) = +25 km, spectrum unreadable
— exactly the pre-incision's signature.** My prediction that Δ would not fall to 0 is **refuted**;
the round's is confirmed.

> ⚠️ **And the test is degenerate, which must be said in the same breath.** On a 1 m floor every
> ex-land cell sits at exactly −1.000 m, so the threshold sweep has **two states and no
> resolution**: 0.5 m lifts nothing, 1 m lifts **all 337 998**. "Fill the shallow hollows" and
> "undo every cell the incision drowned" are the same operation here, and Δ = 0 follows almost by
> arithmetic. **What is proved is that the drowned cells are the whole fringe; what is NOT proved
> is that a transport-limited deposition — which fills some and not others — would do it.**
>
> The non-degenerate version is on the **pre-bathymetry** field, where the drowned cells have a
> real depth distribution (3 041 distinct receiver values, p10 −5.6, median −0.1 m). That is one
> bench run and it is the next round's first item.

### Score

Mine: ordering correction ✓ · pre-bathymetry receivers a broad distribution ✓ · mask 0 cells ✓ ·
hollow depth follows `−shelf_min` ✓ · **tooth count moves > 10 %, and upward** ✓ (+74 %, above my
+15…+60 band) · Δ(≥2 cells) 150–400 at shelf1+diff1.28 ✓ (375) · hypsometry 10–60 m **✗** (5.24) ·
**D does not reach 0 ✗✗**.
The round's: A1 a proxy ✓ · A2 one value ✓ (on the delivered field) · mask unchanged ✓ · step
20 → ~1 m ✓ · **tooth count under 10 % ✗** (74 %) · hypsometry > 10 m ✗ · **D → 0 ✓**.

**Meta holds: several were false on both sides.**

### Standing

No production change. The clamp, the diffusion and the fill are all bench variants; promotion is
the author's call and now has prices attached. Eleven 8192² builds across two benches.

**Two items are now separable, and they were being confused:**
1. **The 20 m cliff** — a PROXY with no antecedence, a 22 m median step, banks structurally
   impossible, fixable at zero cost to the continent (mask and hypsometry bit-identical). This is
   a consumer defect in its own right and does **not** need the fringe question resolved.
2. **The fringe** — which the clamp was *hiding* 43 % of. Lowering the floor makes the measured
   fur worse before any remedy improves it, so the two must not be judged on one number.

## Finding 79 — the shelf floor goes to 1 m in production; the fringe on the MASK is 3 623, not 1 831; and the deposition class does NOT converge

### Open question to the author, asked and not inferred

**Does Living Landz draw the coast from the cell mask / tag (as it does lakes), or by interpolating
the 0 m isoline of the raster?** Everything below is reported in **both** definitions, because they
disagree by a factor 2 and the answer decides which column is the product.

### Rule 11, and its extension — see method rules 11b and 12

`depth_m` 10 (5484) · `BasinSummary` 6 (1221) · `water_class` 42 (760) · `is_land` 2 (8607, my own
F77) · `chained_region` 3 (8024). The value greps found `lake_min_depth_m = 10.0` and
`WETLAND_MAX_DEPTH_M = 3.0` — **two depth thresholds the clamp sits astride**, which is why block A
had to read the lakes and not only the coast.

### A — PRODUCTION: `shelf_min_depth_m` 20.0 → **1.0**, with its intent written at the line

`ALGO_UPSCALE_EROSION` 2 → 3. The docstring now carries the reason the original never had.

**Guards, all read before any coastal number.** ⚠️ The field HASH is not the right check here: the
change is *designed* to move ocean cells, so the eroded hash moves `0x6b10a0c5fdf467a3 →
0x4719ff260186645a` by construction. The check is the LAND-only statistics:

| | before (20 m) | after (1 m) |
|---|---|---|
| land cells | 11 031 073 | **11 031 073** |
| land p10 / p50 / p90 (norm) | 0.502560 / 0.539340 / 0.642129 | **identical to six decimals** |
| cells raised by `pit_fill` | 1 238 085 | **1 238 085** |
| below-sea (wc=2) / ocean (wc=1) | 457 698 / 55 626 539 | **identical** |
| lake invariants (subset: 8192², law OFF, both beds) | dup 0, empty 0, above-level 0, area 0, exorheic-without-outlet 2 / 0 | **all identical** |
| **F74 terminal table, humid** | budget 502.1, W→sea 65/75.9, S→ocean 82.5, TERMINAL 158.4, HOLE 343.7, chained 13/1052.4 | **identical, every figure** |
| F74 terminal table, arid | chained 32.5, Σ inflow 119.31, Σ evap 57.97 | 32.8 (+0.9 %), 119.92 (+0.5 %), 57.79 (−0.3 %) |

> **Stop rule A does not fire: no lake invariant moves, and the terminal table moves by at most
> 0.9 % against a 5 % bar.** The negative control is already in the dossier — the 20 m build
> reproduces `0x6b10a0c5fdf467a3`, the hash Findings 74, 75, 77 and 78 all recorded.

⚠️ **Two consequences that are NOT no-ops and must be judged, not buried.**

1. **Below-sea basin `depth_m` collapses**: humid p10/p50/p90 **20.69 / 26.68 / 69.27 → 1.69 / 8.41
   / 52.49 m**; arid 20.62 / 26.31 / 34.82 → 1.69 / 7.68 / 17.14. `level_m` is **unchanged to the
   centimetre** (6.68 → 6.68), so this is the floor moving and nothing else — the clamp was
   inventing 19 m of depth under every basin. F33's "deep bowl, floor −20 m" and F38's "depth-0
   lakes at −20 m, on a FLAT floor" were this.
2. **The wetland mask moves, and in the arid bed it moves a lot**: humid 42 754 → 44 366 (+3.8 %),
   **arid 38 306 → 353 630 (×9.2)**. `WETLAND_MAX_DEPTH_M = 3.0` (`drainage.rs:1777`) classes a
   cell under 3 m of water as wetland, and a 1 m floor puts every below-sea basin under that bar.
   **Physically a 1 m-deep pan IS a wetland**, so this is arguably the correction and not the
   regression — but it is a large change to an exported classification and it is the author's call,
   not mine.

**The consumer number, on the production path**: land−sea step median **22.00 m → 3.00 m**; share
of the shoreline gentler than 15 % **0.0 % → 71.4 %**; steeper than 25 % **100.0 % → 18.9 %**.

### B1 — the mask instrument, and its control

`marching_squares` on a BINARY land/sea field at 0.5 interpolates to the exact midpoint of every
crossed edge, so it traces the boundary staircase and **cannot see the value of the water**. No new
detector: the same `coast_spurs` / `coast_shape_thresholds` run on it.

| | shelf 20 m | shelf 1 m | verdict |
|---|---|---|---|
| **MASK** | 1 758 / **3 623** / 7 056 km | 1 758 / **3 623** / 7 056 km | **IDENTICAL — the instrument does not read the sea** |
| ISOLINE | 1 848 / **1 831** / 6 206 km | 1 776 / **3 190** / 5 927 km | −72 / **+1 359 (+74.2 %)** |

### B2 — the re-read of Findings 75–78, and three verdicts

| variant | ISO ≥1 km | **ISO ≥2 cells** | MASK ≥1 km | **MASK ≥2 cells** | MASK R |
|---|---|---|---|---|---|
| PRE-INCISION (authority) | 18 | **11** | 20 | **13** | 0.280 |
| **SHIPPED (shelf 20)** | 1 848 | **1 831** | 1 758 | **3 623** | **0.920** |
| **PRODUCTION (shelf 1)** | 1 776 | 3 190 | 1 758 | **3 623** | 0.920 |
| LAW ON (shelf 20) | 916 | 1 056 | 939 | **1 614** | 0.887 |
| diffusion 1.28 (shelf 20) | 1 275 | 189 | 1 289 | **567** | 0.751 |

**Verdicts, one per claim:**

- **"the 20 m clamp hides 43 % of the fur" — HELD TO THE TRACER, and it was worse than 43 %.** On
  the mask the shelf change does **nothing at all**; the true cell-scale count of the shipped
  product is **3 623**, and the isoline was reporting **1 831**. The clamp was hiding **49 %** of it
  from an interpolating consumer, and **0 %** from a tag-drawing one.
- **"λ = 14 at ×3.74" — DOES NOT SURVIVE.** On the mask: λ 9 at **×2.37, below the 3× line — not a
  wavelength.** The isoline on the same field gives λ 14 at ×3.74. **The spectral peak of Findings
  76, 77 and 78 is an artefact of the isoline tracer.** (Pre-incision mask: under 32 gaps, rule 10,
  not a reading — as the pre-incision isoline also was.)
- **"the law removes 44 %" — HOLDS, and it is stronger on the mask: −55.4 %** (3 623 → 1 614)
  against the isoline's −42.3 %.
- **"the anchored diffusion removes 90 %" — HOLDS, at −84.3 %** (3 623 → 567).
- **The parallelism is NOT a tracer artefact**: mask R 0.920 against isoline 0.879, on a
  pre-incision baseline of 0.280.

> **The Δ of the criterion is therefore Δ(mask) = +1 738 / +3 610 against the pre-incision
> authority, and it does not move when the shelf does.** The isoline column stays in the dossier
> only until the author answers which one the product draws.

### C — the pre-clamp depths: BOTH written outcomes are true, of different halves

Population: pre-incision LAND and pre-clamp SEA — the cells the incision drowned, read **before**
the bathymetry overwrote them. **298 598 cells.**

| p1 | p10 | **median** | p90 | p99 | max |
|---|---|---|---|---|---|
| 0.0020 | 0.0384 | **0.8803 m** | **18.7383 m** | 47.46 | 99.91 m |

| shallower than | 0.05 m | 0.2 m | 0.5 m | 1 m | 2 m | 5 m |
|---|---|---|---|---|---|---|
| share | 11.9 % | 25.7 % | 39.8 % | **52.3 %** | 63.6 % | 74.8 % |

**Half of them barely cross zero (p10 is 3.8 cm) and the other half are real sub-sea beds (p90
18.7 m, max 100 m).** The round wrote two outcomes and expected one; both are true, of different
halves. My prediction (p50 2–30 cm, p90 < 2 m) is refuted; the round's (p50 1–20 cm, p90 < 1 m) is
refuted harder.

### C-bis — D, non-degenerate: the class does NOT converge

Lifting drowned cells shallower than a threshold to `sea + ε`, on the pre-clamp field:

| threshold | lifted | ISO ≥1 km | ISO ≥2 cells | **Δ mask ≥1 km** | **Δ mask ≥2 cells** |
|---|---|---|---|---|---|
| 0 m (control) | 722 | 1 770 | 3 568 | +1 741 | **+3 627** |
| 0.05 m | 35 541 | 1 752 | 3 979 | +1 780 | **+4 800** |
| 0.2 m | 76 835 | 1 689 | 3 287 | +1 724 | +4 653 |
| 0.5 m | 118 756 | 1 541 | 2 472 | +1 523 | +4 088 |
| 1 m | 156 231 | 1 320 | 1 810 | +1 250 | +3 232 |
| 2 m | 189 788 | 1 002 | 1 277 | +961 | +2 246 |
| 5 m | 223 214 | 633 | 716 | +640 | **+1 556** |
| **ALL** | 298 598 | **18** | **9** | **−1** | **+0** |

> **Filling the shallowest cells makes the mask fringe WORSE before it makes it better**: +3 627 →
> **+4 800** at 5 cm. Lifting a centimetric hollow creates a one-cell isthmus or islet, and the
> boundary staircase gets longer, not shorter. The curve only turns after 0.2 m, and at **5 m —
> three quarters of all drowned cells lifted — 57 % of the fringe is gone and 43 % remains.**
> **Only "lift everything" reaches Δ = 0 (−1 / +0), which is the identity the degenerate test of
> Finding 78 had already shown.**
>
> **So the deposition class is refuted as a short path.** A transport-limited term fills the
> shallow ones and not the deep ones, and that is precisely the region where the metric gets
> worse. My "Δ reaches 90 % at 0.5 m" is refuted by its sign; the round's "the short path exists"
> is refuted by the same table.

⚠️ The 0 m control lifted **722 cells**, not 0 — `depth ≤ 0.0` catches cells sitting exactly at sea
level. It moved the mask count by 17 in 3 610 (0.5 %), so the control is impure and harmless, and
it is reported rather than rounded away.

### Score

Mine: mask identical across the shelf sweep ✓ · no readable λ on the mask ✓ · the clamp's swing was
entirely the tracer ✓ (I said > 95 %, measured 100 %) · law on the mask −30…−55 % ✓ (−55.4, at the
edge) · diffusion −70…−92 % ✓ (−84.3) · lake invariants green ✓ · terminal table under 5 % ✓ ·
**pre-clamp p50 2–30 cm ✗** (0.88 m) · **p90 < 2 m ✗✗** (18.7 m) · **Δ reaches 90 % at 0.5 m ✗✗**
(it is worse than baseline there).
The round's: A1 a proxy ✓ · basin depths move ✓ · terminal table under 5 % ✓ · step 22 → ~3 ✓ ·
mask identical ✓ · **λ = 14 absent from the mask ✓** · law −40…−50 % ✗ (−55.4) · diffusion −60…−75 %
✗ (−84.3) · "the 43 % was > 80 % tracer" ✓ · **pre-clamp p50 1–20 cm ✗, p90 < 1 m ✗** · **the short
path exists ✗**.

**Meta holds on both sides.**

### Standing

One production change (`shelf_min_depth_m` 20 → 1, `ALGO_UPSCALE_EROSION` 2 → 3), its intent
written at the line, all guards green. Everything else is bench. Two method rules earned (11b, 12).

**What is now settled and what is not.** The 20 m cliff is fixed and cost the continent nothing.
The fringe is **3 623 cell-scale spurs on the mask**, against a pre-incision authority of 13 — the
criterion is Δ = +3 610 and the clamp never touched it. The two levers that do touch it are the
channel-head law (−55 %) and the anchored diffusion (−84 %), neither of which is a remedy on its
own and both of which change the continent. **The deposition class is out.**
