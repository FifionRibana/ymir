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
**S1 13 m · S2 21 m · S3 37 m · S4 68 m · S5 114 m** — sanity-matching a Thames-scale trunk
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
