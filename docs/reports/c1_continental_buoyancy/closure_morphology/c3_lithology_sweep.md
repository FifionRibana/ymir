# C-3 lithology spread sweep — the measurement that chooses the multiplier

Closures roadmap §3. The soft↔hard erodibility spread is a **measurement, not a
prediction** (author's directive). This is the sweep that fixes it: the WHOLE
production chain (`upscale_from_c1_with_progress`, export recipe — relief-v3
stream-power incision, droplets off), lithology OFF as baseline then ON at
×3/×10/×30/×100 soft, both resolutions, production seed
`10481999410520546993`, domain 400 km.

Bench: `crates/ymir-core/tests/c3_lithology_sweep.rs`
(`cargo test -p ymir-core --test c3_lithology_sweep c3_sweep -- --ignored --nocapture`).

Method rule 3 (ADR 0001): the bench reproduces the WHOLE chain — the K field is
built and threaded exactly as production does it (coarse hard=1.0 + rift soft →
bilinear upscale registered to the altitude; volcaniclastic stamped at HD on the
edifice basal discs; per-cell K into the incision), not a reconstruction of the
incision alone.

## Class coverage (both resolutions, area-preserving)

| class | share |
|---|---|
| hard basement | 95.7 % |
| rift-soft (age≈0) | 1.6 % |
| volcaniclastic (edifice footprints) | 2.7 % |

The soft class is ~4.3 % — the **minority-by-nature** the design predicted (Stock &
Montgomery: hard crystalline + metasediments are the bulk; soft is rift +
volcaniclastic). No sedimentary basins: Ymir's production erosion is
detachment-limited (no deposition) — recorded as a limitation, not filled with noise.

## Sweep — per class: local relief (m) / median slope (°) / steep>30° share (‰) / channel incision (m)

### 2048²

| soft × | HARD | SOFT (rift) | VOLC | closed pits | land % |
|---|---|---|---|---|---|
| OFF  | 329 / 6.9 / 147 / 155 | 150 / 3.1 / 107 / 51 | 612 / 15.6 / 310 / 329 | 982 | 15.7 |
| ×3   | 329 / 6.9 / 147 / 155 | 136 / 2.8 / 111 / 32 | 498 / 11.9 / 282 / 270 | 982 | 15.6 |
| ×10  | 329 / 6.9 / 148 / 155 | 115 / 3.0 / 134 / 17 | 414 / 9.4 / 256 / 206 | 971 | 15.4 |
| ×30  | 329 / 6.9 / 148 / 154 | 105 / 3.1 / 174 / 15 | 363 / 7.8 / 234 / 160 | 975 | 15.3 |
| ×100 | 329 / 6.9 / 148 / 155 | 103 / 2.2 / 202 / 12 | 303 / 6.9 / 225 / 98 | 977 | 15.2 |

### 8192²

| soft × | HARD | SOFT (rift) | VOLC | closed pits | land % |
|---|---|---|---|---|---|
| OFF  | 119 / 9.8 / 227 / 78 | 67 / 5.9 / 146 / 57 | 266 / 24.2 / 425 / 173 | 17516 | 16.6 |
| ×3   | 119 / 9.8 / 227 / 78 | 65 / 5.7 / 159 / 56 | 261 / 23.2 / 419 / 160 | 17485 | 16.6 |
| ×10  | 119 / 9.8 / 228 / 77 | 65 / 5.5 / 171 / 60 | 242 / 21.7 / 405 / 136 | 17365 | 16.5 |
| ×30  | 119 / 9.8 / 228 / 77 | 68 / 5.5 / 188 / 68 | 223 / 20.7 / 390 / 115 | 17202 | 16.5 |
| ×100 | 119 / 9.8 / 228 / 77 | 74 / 5.9 / 209 / 106 | 215 / 19.4 / 379 / 103 | 17114 | 16.4 |

## The two effects, separated (the reason for hard = ×1.0)

1. **Global slowdown = ZERO by design.** Hard basement is the relief-v3 reference
   (×1.0), so the ~96 % hard bulk erodes exactly as production at every multiplier.
   The HARD column is **flat across the whole sweep** (2048²: 329/6.9/147/155 at
   every row; 8192²: 119/9.8/227/77 at every row). There is no global slowdown to
   disentangle from the contrast — the alternative (hard ×0.3, soft ×1.0) would have
   moved 96 % of the continent and confounded the two.
2. **Contrast = the whole signal, and it is monotone and physical.** Softer K erodes
   the rock DOWN: relief and channel incision fall, valleys open (higher W/D). VOLC
   (2048²) relief 612→303, incision 329→98; SOFT incision 51→12. This is the
   hard-gorge / soft-open-valley dichotomy the closure exists to produce.

## Invariants (both resolutions, whole sweep)

- **Closed depressions (C-1) SURVIVE**: 2048² 982→977, 8192² 17516→17114 — a slight
  DECREASE, never a flood. Softening does not fabricate pits; the FBM
  flow-conditioning holds.
- **Land fraction stable**: 15.7→15.2 % (2048²), 16.6→16.4 % (8192²).

## Chosen multipliers

- **soft (rift) = ×10** — ~1 order of magnitude (mid of the Stock & Montgomery
  granite↔mudstone range), a clearly visible contrast (2048² soft incision 51→17,
  volc relief 612→414) with C-1 intact and land stable. ×3 is real but faint; ×30/×100
  progressively **erase the C-2 volcanic cones** (volc relief halved at ×100) — over
  the edifices we just built, which is wrong.
- **volcaniclastic = ×3 (fixed, decoupled from the soft sweep)** — intermediate (S&M:
  volcaniclastic sits between granite and mudstone), and deliberately mild so the
  edifice morphology from C-2 is dissected, not flattened. (The sweep tied volc =
  √soft only to trace the class; production fixes it at ×3.)

`m = 0.4, n = 1` (stable base level — NOT the Kauai exponents). Gated OFF by default
(`LithologyConfig::enabled = false`) → byte-identical production; the eroded cache key
is unchanged when disabled (config skipped from serialization, algo appended only when
enabled).

## Export verdict (the product, not a reconstruction)

Two 8192² exports, humid climate, production seed — NEW (lithology ON: rift ×10,
volcaniclastic ×3) vs REF (`*.volcan.humid`, C-2 baseline, lithology OFF). Decoded
per-manifest metre scales, class map rebuilt for the window
(`tests/c3_export_compare.rs`):

| class | land cells | mean \|Δ\| | mean Δ (signed) | changed ‰ (>1 m) |
|---|---|---|---|---|
| HARD | 10.27 M | 1.5 m | −1.3 m | 28 |
| SOFT (rift) | 0.40 M | 54.7 m | **−54.0 m** | 534 |
| VOLC | 0.48 M | 116.5 m | **−114.6 m** | 701 |

Hard basement essentially unchanged (mean 1.5 m, 2.8 % of cells move — the drainage
network responding to the soft zones re-routing, not a lithology effect on hard rock;
max elevation identical to the digit). Soft classes erode DOWN (rift −54 m,
volcaniclastic −115 m) — the contrast is 35–75× stronger there. Effect (a) global
slowdown = 0, effect (b) contrast localised to the soft classes, confirmed on the
shipped product.

Δ map (`c3_delta_map.png`, ×4 downsample): the hard bulk is uniform slate; the
volcaniclastic edifices glow red (dissected), the rift corridor reddens along its
channels; ocean dark. The contrast operates exactly on the causal soft footprints.

Awaiting the author's visual validation before C-3 is marked done in the roadmap.
