# ADR 0001 — Erosion: coastal sediment sink, and the terrace / dendritic-valley findings

Status: accepted (2026-08). Scope: C1 HD relief (erosion, FBM upscale, isostasy).
Diagnostic harness: `crates/ymir-core/tests/terrace_diagnosis.rs` (all `#[ignore]`),
commits `019113f`, `76ae09f`, `94fe11b`. Re-run:
`cargo test -p ymir-core --test terrace_diagnosis --release -- --ignored --nocapture`.

This ADR records four findings from a measurement campaign on why the C1 HD relief
shows (a) terraces following isolines and (b) no carved dendritic valleys under a
fully-computed river network. Each is stated with the numbers, because the
conclusions are counter-intuitive and will otherwise be re-litigated.

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

**Measurement (real field, incision = fraction of top-1 % flow cells in a local
altitude minimum; FBM floor 7 %).** No `(sink, density)` pair satisfies "incision
rises AND the Strahler histogram deepens":

| f | droplets/cell | net % | carved (Δ) | maxStrahler | confluences |
|---|---|---|---|---|---|
| 1.0 | 0.95 | +1 | 11 % (+4) | 6 | 2366 |
| 0.0 | 0.95 | +59 | 11 % (+4) | 6 | 2183 |
| 0.25 | 4.0 | +76 | 18 % (+11) | **3** | 649 |
| 0.0 | 4.0 | +87 | 17 % (+11) | **3** | 556 |

The sink leaves incision at 11 %; more droplets raise it to ~18 % only by fragmenting
the network (maxStrahler 6→3, confluences 2366→649 — grain, not drainage). On a
synthetic smooth cone at production parameters, surface roughness barely changed
(0.00005 → 0.00006): **no incision even on an ideal input.**

**Root cause.** Droplets are stochastic and UNCORRELATED — nothing makes neighbouring
rills converge into a shared, deepening channel, so a hierarchical valley network
cannot emerge and heavy droplets just pit the surface.

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

---

## Finding 4 — H3 terraces originate in the coarse isostasy / equilibrium-height, not erosion

**Context.** The relief shows terraces PARALLEL to contours (concentric loops around
hills), with jumps of several hundred metres, at both 8192² and 2048².

**Measurement.**
- u16 quantisation REFUTED: 0.1–0.14 m/unit, vs observed jumps of 100s of m.
- Step height ~**120–176 m**, roughly resolution-INDEPENDENT (512²/2048²) → set by a
  resolution-independent source (coarse tectonic/isostasy field).
- Terrace-source disentangle (flat fraction of a transect): **pure bilinear
  (coarse only) 13 % → FBM 6 % → after erosion 24 %.** So FBM ROUGHENS (not a source);
  erosion DEPOSITION re-flattens (secondary). The flat plateaux originate in the
  coarse post-isostasy field.

The step height matches the C1 equilibrium-height closures producing flat crustal
levels. Note the scale: a 120–176 m step exceeds the relief of a hamlet-scale valley
(~150–250 m wide), so this is not cosmetic.

**Decision.** Since erosion (Finding 3) does not dissect these steps, they must be
addressed at the SOURCE — the isostasy / equilibrium-height closures — a separate
job. Not touched here.

**Consequences.** Terraces remain visible until the closure work is done. See
Finding 3's successor and the recommended order in the campaign summary. When the
equilibrium-height closure is touched, re-check the pre-existing Picard
non-convergence in the same tectonic layer —
[docs/issues/picard-nonconvergence-rectangular-smoke.md](../issues/picard-nonconvergence-rectangular-smoke.md).
