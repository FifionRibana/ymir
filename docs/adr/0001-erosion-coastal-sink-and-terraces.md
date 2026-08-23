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
