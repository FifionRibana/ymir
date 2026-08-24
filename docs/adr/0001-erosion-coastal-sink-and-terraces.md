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
