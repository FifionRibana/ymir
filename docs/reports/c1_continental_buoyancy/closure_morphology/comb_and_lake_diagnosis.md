# Diagnostic: the river "comb" artefact, and lake abundance

Read-only, on the SHIPPED production export (`seed10481999410520546993_8192.ymir`, humid,
lithology + C-3b ON). Benches: `tests/comb_diagnosis.rs`, `tests/lake_abundance_diagnosis.rs`
(both `#[ignore]`, parse the export files). Two distinct subjects.

## Subject 1 — the comb (a defect): mono-D8 planform quantization, not flat-ground

**Verdict: the striking parallel rake is in the POLYLINES (D8 direction quantization of
headwater tributaries on planar hillslopes), riding real-but-SHALLOW terrain grooves. It is
NOT a flat-ground tie-break, and none of the three flat-ground suspects is the cause.**

### The first discriminant — terrain vs polyline (cross-section ⊥ the segment)

Mean Δelevation vs the channel cell, across ±12 cells (×49 m), comb (order 1, axiality > 0.9)
vs trunk (order ≥ 3):

| offset (cells) | −12 | −8 | −4 | −2 | −1 | 0 | +1 | +2 | +4 | +8 | +12 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| COMB (o1) | 14 | 13 | 13 | 9 | 6 | 0 | 5 | 8 | 11 | 11 | 5 |
| TRUNK (o≥3) | 112 | 91 | 58 | 29 | 18 | 0 | 17 | 26 | 53 | 83 | 103 |

Groove depth (mean flank − channel): **COMB 12 m, TRUNK 87 m.** So the comb channels are REAL
but shallow (a headwater rill), riding on ~12 m grooves — the channel rises immediately at ±1
(no wide flat floor). Not spurious polylines on smooth ground; not lateral-widened flats.

### Grid vs flow, and the straightness

- Per-segment axiality (concentration R of its own step directions; 1.0 = perfectly straight):
  p10/p50/p90/p99 = **0.55 / 0.79 / 0.90 / 1.00**. Many segments are PERFECTLY straight
  (p99 = 1.00) — unnatural; real channels meander. This is the comb fingerprint.
- Comb candidates (order 1, axiality > 0.9, ≥ 5 pts): **589 = 8.2 %** of segments.
- Slope UNDER comb cells: **median 6.0°, only 14 % < 1°** (trunk 5.4°, 19 % < 1°). The comb is
  NOT on flat ground — it is on planar hillslopes.
- Step-direction histogram is (trivially) only on the 4 grid axes — every D8 path is; that
  alone does not distinguish the comb. The discriminator is the STRAIGHTNESS + the slope.

### Why the three suspects are NOT it (measurement + config, not argument)

The code SELF-DOCUMENTS the mechanism: `terrain/flow.rs:45-51` — the D8 trace "collapses to
CARDINAL combs… grid-rendered parallelism is a floor of drawing rivers as grid cells." Routing
is mono-D8 (`compute_d8`, `flow.rs:296`), tie-break by fixed neighbour order N-first.

1. **Upstream extension (head_km2, full_tree)** — production uses `full_tree = false`
   (`hd.rs:669`, main stem only), so it does NOT export a ramified comb tree. The order-1
   segments are the real leaf tributaries (one, `basin 3042`, carries 31 km² — a substantial
   channel, not a sub-A_c tail). Head extension only lengthens the teeth slightly.
2. **Lateral widening** (`stream_power.rs:448-488`, relief-v3 `K_lat = 4.0`) planes flat floors
   on TRUNKS (A ~ 10⁴), not order-1 headwaters — and the comb cross-section rises immediately
   (no flat floor). Not the order-1 comb cause.
3. **breach_monotone** carves axial grooves across FLATS — but `FlatPerturbation` is ON by
   default (`drainage.rs:193`), mitigating the flat tie-break, and the comb cells are on 6°
   slopes, not flats.

**The real cause is more fundamental than the three suspects: D8 has only 8 directions, so
parallel steepest-descent fall-lines on a uniform (FBM-textured) hillslope all round to the
SAME grid axis → a regular parallel rake (Pass-1 steepest descent, which no flat-ground
mitigation touches).** It is a drainage-EXTRACTION artefact (the polyline planform), not
erosion; the terrain is only mildly grooved. The fix direction (NOT implemented): D∞ routing
(`compute_dinf` exists, off by default) for continuous directions, or planform de-quantization
of the exported polylines. The erosion is not the problem.

## Subject 2 — lake abundance (a design question): the frozen sill

**Verdict: in humid climate 99 % of lakes are exorheic and filled to their sill, and NO sill
incision exists, so every basin is a PERMANENT lake. The France↔Scotland dial is sill
incision, not the noise.**

### Population (lakes.json, 92 lakes)

- By type: **91 Exorheic + 1 CraterAcidic; 0 Endorheic.** 44 below-sea basin ids; 5 shallow.
- Total lake area **6220 km²**; area p50/p90/max = 7.3 / 256 / 868 km²; depth p50/p90/max =
  46 / 357 / 1194 m.
- **The largest 10 lakes hold 74 % of the lake area** (4628 / 6220 km²) — the "Scotland" look is
  a few deep sill-perched basins, not scattered FBM ponds (only 5 shallow).

### Why they fill and stay (recon, Finding 39)

- `detect_lakes` marks a cell as water iff `filled − eroded > 1e-6` — EVERY pit-filled
  depression is a lake. The FBM creates the hollows (16 → 90 682), C-1 conditions them, so the
  lake count tracks the depression count.
- `water_balance_lakes` (`drainage.rs:768`): `a_eq = inflow / pe_lake`; if `a_eq ≥ a_sill` →
  Exorheic, `level = surface_elevation` (the pit-fill sill), UNCHANGED. Humid net_evap ≈ 0 →
  `a_eq` huge → essentially everything exorheic, filled to its col.
- **No sill incision anywhere** — the outlet level is the geometric sill, never lowered. In
  nature the outlet CARVES its sill and drains the lake (post-glacial lakes vanish this way);
  the model freezes it → permanent lake.

### The two dials, measured

- **Climate sets the exorheic share.** Across comparable exports: arid-HOT → 22 % endorheic
  (11/49); humid / arid-cold / tropical → 0 % endorheic. So high PE (hot + arid), not aridity
  alone, flips lakes endorheic (they shrink/dry). Humid → 0 % → every basin exorheic.
- **Sill incision decides permanence.** An exorheic lake HAS an outlet → its sill would incise
  and drain it; the frozen-sill model keeps it full. Estimate: a sill-incision pass would drain
  **~100 % of the lake area** here (91/91 exorheic), leaving only the endorheic + crater residue
  (1 lake). THIS is the France↔Scotland dial.

### Blast radius of sill incision (do NOT implement)

Partially draining exorheic lakes turns each outlet into a river, shrinks/removes lake
footprints, recomputes levels and regimes, moves river mouths and the coastline, and shifts
biomes — it reopens the whole hydro chain stabilised over ~15 passes (drainage + lakes + rivers
+ coast + biomes). The France/Scotland dial belongs here, as an explicit outlet-incision stage,
NOT in the FBM noise.
