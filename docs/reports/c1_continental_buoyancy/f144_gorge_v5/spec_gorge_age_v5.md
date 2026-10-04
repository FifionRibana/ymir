# Specification v5 — the gorge's retreat by age (Finding 144, BUILT GATED, OFF by default)

**v5 = v4** (`f142_gorge_v4/spec_gorge_age_v4.md`) **+ the author's decisions of 2026-10-02 + the parts S, B, H and P
of F144.** Everything in v4 not restated stands. Tags as before: MEASURED, ANCHORED, PROXY, DECISION.

## 0. The decisions (recorded as given, ADR Finding 144)

- **The gorge's retreat applies to every input lake outside the below-sea basins (17)** (DECISION).
- **The outlet falls' height: a cap by the river's size (~70–130 m)** (DECISION).

## 1. The parameters

`GorgeRetreat::v5(r, p)` = v4 (m = 10, A_ref = 418.7 km², `plain`) with three gated fields, each absent from the key
when off:
- `input_scope` (S);
- `bound_b` (B);
- `head_cap` (H).

The bench option `light_mode` (P0–P3, with `light_dt_m` for P2) is absent from the key at 0. `scope_lows` (F143's
list) still takes precedence when given; v5 does not set it.

## 2. What changes against v4

### S — the scope, computed on the input

A body is kept iff:
- it is an **input lake** (a lake of the pre-drainage);
- **none of its cells lies in a `basin_base` closed depression** (`spill > bf + 1 cm` on the breached input, which
  includes the merged below-sea water bodies of F38);
- and **none lies under the sea** (`!land`).

This replaces v4's ring proxy and F143's per-world list. Checked against **the gorge world's own lakes** (G-bodies),
not ON's.

### B — L(r) bounded by the outlet's base

`L = min(L_in, max(L(r), B))`. B is the first base on the outlet path from the col (the skeleton's D8), as the χ walk
reads it:
- the sea's `base_m`;
- an inland basin's spill;
- another kept body's level, or another input lake's surface;
- a closed depression's spill;
- or a land terminal's height.

It is iterated, so that a chain of lakes sees its downstream levels.
- **A lake bounded by the sea (B = base_m)** stays a lake at the sea's level, joined to it by its gorge, if its
  L_floor lies below the sea. It is declared as such, **not made an arm of the sea** without the author's decision.

### H — the capped head fall

- **The head fall is min(φ·D_g, H_cap(A))**, with H_cap(A) = 100 m·(A / A_ref)^−0.5 (PROXY). The shape mirrors S_loi
  ∝ A^−0.5. The 100 m is the order of the large rivers' falls (Victoria, Iguazú, Niagara), **to be sourced before any
  use**. The gorge takes the rest.
- **H_f** = max(1.5·S·100 m, 10 m) with the gorge's real slope S = max(m·S_loi, S_req), capped. **The code already
  used it** (read in F144 Partie 0).
- **G-tag's instrument is corrected** (F143-L): a gorge step is a two-cell drop exceeding the designed slope's drop
  over the same span by more than H_f.

### P — the light pass and the designed geometry (bench candidates)

The design mask is C2's (F143):
- the kept bodies' footprints;
- their rings (8-neighbours);
- the gorge path cells with their Chebyshev radius-2 neighbourhood.

The candidates:
- **P0**: the light pass free (the reference, F142 / F143-C0).
- **P1, the floor**: after the light pass, `z ← max(z, construction)` on the mask. The rim does not go under L(r), the
  corridor not under z_gorge(s), the plain not under its filled level, a lake's bed not under its input level.
  Deposition stays allowed.
- **P2, the transition**: where the light pass lowers a cell, the drop is weighted by w(d) = min(1, d / d_t), with d
  the cell's distance (m) from the mask. It is 0 on the design and 1 from d_t on. **d_t** is declared in
  `f144_declared.md`, measured from C2's walls.
- **P3, the order**: the light pass runs on the construction **without** the gorge's invariant, the rim clamp and the
  plain (`bare`). After it, the design is laid back as `z ← max(z, construction)` on the mask, before the droplets
  and the breach. The droplets and the breach can still undo it (measured).
