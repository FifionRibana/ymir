# F143 — instruments, declared BEFORE the measurements (2026-10-03)

**The worlds** are the témoin, with `GorgeRetreat::v4` plus the exact scope (`scope_lows`, below), and k = F121_AGE_K ×
age (as F140 / F142). "The construction" is `skeleton` + `carve` + `gorge_plain` on S1. The other stages are built by
`build_world` with the knobs of F142:
- (b) the light pass: `erosion_off`, `bathymetry_off`;
- (d) the full world: the breach's input.

## S — the exact scope (bench `f143_l`)

- ON ×1's 14 D8 lakes (F139's set) are matched to the v4-proxy bodies at ×1, p 0.5, each to the body with its
  largest footprint overlap. Those bodies' lowest cells are the list.
- **The list is a per-world measurement.** The construction knows no final lake: F132-P4's circularity.
- The dropped bodies are printed with their class flags (why the proxy kept them) and their fate in ON: the ring's
  minimum against L_in on ON's construction and on ON's world (d), and their final lake.

## L — localisation (bench `f143_l`)

- **G-ring** (GORGE ×1, p 0.5, the construction, `carve_diag` for the attribution). Each violating pair is a clamped
  ring cell c (S1 > rf + 5 cm, |zb − rf| < 5 cm) with an outside neighbour m more than 28° lower. For m:
  - cut or not (S1 against zb);
  - its laying sample's line: **the outlet line** (the line laying the col) or **another line**;
  - and the count again on the light pass's world (b).
- **G-drained** (×1.2 p 0.5 and ×1.4 p 0, the construction): the cells upstream of a drained body (its skeleton
  D8) under L_floor − 1 m. For each:
  - the laying sample (`carve_diag`'s `who`): the body its cell is in, its zf, and whether its line runs to the same
    body;
  - counted by category, with examples;
  - at `carve` and after the plain.
- **G-tag** (×1 p 0.5, the construction): each 2-cell drop > H_f along the gorge path beyond the lip. For each:
  - its distance from the col;
  - whether its cells are on the carved floor;
  - the lines laying them (the outlet line or another);
  - its body.
- **G-sea** (×1.4 p 0.5): the land cells ≤ sea after the breach that ON ×1 lacks. Through F141's instrumented copy of
  the breach (re-copied, checked bit for bit against production):
  - each cell's pit, floor, spill and body, and the ramp's base;
  - the cells' height on the construction and at (d);
  - their distance from the original coast (the sea cells of (d)), to say whether the sea advances inland.

## θ — the rebased links (bench `f143_l`)

On OFF's trunk links (the construction's mask, without lake / depression links, as F142):
- the rebased links: a cell whose GORGE ×1 p 0.5 base is > 0.5 m under ON's;
- control: the unchanged links.

For each, the local slope S = Δz / length against the law's slope S_law = k·max(A, a_c)^−0.5 (A km², the GORGE
skeleton's area at the donor), as the ratio S / S_law:
- on the construction, and on the light pass's world (b);
- binned by the distance upstream from the nearest kept body: 0–1, 1–2, 2–5, 5–10 and > 10 km;
- medians and the share of links with the ratio outside [0.5, 2].

## H — the head falls (bench `f143_l`)

- At ×1, p 0.5 on the 14: A, D_g, φ, the head fall, S and H_f, from the skeleton.
- **Caps**: 50 m, 100 m, 200 m, none, and **H_cap(A) = 100 m · (A / 418.7 km²)^−0.5**. Provenance: PROXY. The shape
  mirrors S_loi ∝ A^−0.5; the 100 m is the order of the large rivers' falls named by the reviewer, not yet sourced.
- Per cap: the tagged count (min(head, cap) > H_f), the median and the max of the capped tagged heads.
- **The lip drop observed in C0** (the ring's minimum at (b) − L(r), per body) against the head fall.

## C — the controls (bench `f143_c`), declared before it runs

- **The worlds**: ×1 (p 0.5), ×1.4 (p 0.5), ×1.4 (p 0), each under:
  - C0 (v4 + the exact scope);
  - C1 (+ φ = 0);
  - C2 (+ `freeze_design`);
  - C3 (C1 + C2).
  - ON ×1 first (the references); ON and OFF ×1 last (the controls).
- **`freeze_design`** is a bench option. The designed geometry is restored after the light pass: the kept bodies'
  footprints, their rings (8-neighbours), and the gorge's path cells with their Chebyshev radius-2 neighbourhood.
- Each world runs isolated (`catch_unwind`, twice as in F142).
- **On the construction** (C0, C1; C2 = C0 and C3 = C1 by construction, not re-run):
  - G-pits;
  - G-levels and G-area (the construction's own lakes: climate-free drainage of the construction);
  - G-rim;
  - G-ring, G-drained and G-tag (as F142);
  - G-slope100 (the construction's own outlets);
  - G-sea (the breach on the construction against the breach on ON's construction at the same age);
  - θ.
- **On the final world**: every F142 gate as in `f142_g`, with G-ring, G-drained and G-tag also on the eroded world.
  Also:
  - canyons, the coast (the extras), θ;
  - rule 14 (removal and raise; the plain's volume);
  - G-sea against ON ×1's 292;
  - F38.
