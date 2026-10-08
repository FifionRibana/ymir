# The gorge's retreat by age — PAUSED (2026-10-08), the state for whoever resumes it

*Written at F157. One page. The author's decision (2026-10-08): « Pause, avec bilan de clôture ». Everything below is
gated and off; the viz toggle is hidden (`GORGE_TOGGLE_VISIBLE = false`). The production world does not carry it.*

## What the author asked, and what stopped it

- The criterion (2026-09-29 → 10-02):
  - the lakes shrink with age;
  - a present lake is a base level; an emptied one is not;
  - the below-sea basins' lakes are present lakes;
  - « Chute ou gorge, les deux sont valides »;
  - « C'est bien d'avoir l'âge qui fixe jusqu'où la gorge a reculé ».
- **The stop** (F157), after the hillshades of T / T-all / P0 / NP: « Je ne suis pas convaincu. Je ne vois pas où
  s'arrête le lac, mais je vois un rebord qui apparaît en T-all par rapport à P0 (et en T d'ailleurs). »
- **The reviewer's reading** (F157): the visible rim is the ring's clamp (G-ring). P0 erases it because the light pass
  erodes the ring freely; once the designed geometry is protected (T, T-all), the edge stands. No setting of the light
  pass gives a held rim without an edge: the defect is in the construction.

## What is built (gated, toggle hidden)

- **`GorgeRetreat`**, `tectonics_c1/valley_construction.rs`:
  - **v3** (F140): the retreat per lake r_lake = min(2, r_world·(A/A_ref)^p), p ∈ {0, 0.25, 0.5} (default 0.5), map
    (a) of the age selector; the level L(r); the rim clamp; the gorge's slope with m = 10 and the « raidir » rule;
    the head and shortage falls, tagged.
  - **v4** (F142): the D8 scope and the minimal plain (`gorge_plain`).
  - **v5** (F144): the input scope (**18 bodies**); L bounded by the outlet's base; the head fall capped by the
    river's size (H_cap(A) = 100 m·(A/A_ref)^−0.5, PROXY).
  - **v6** (F155): the falls' rule, the current φ draw except a soft lip (no fall). The soft cells come from F154's
    substratum.
- **The light pass's candidates** (bench options of `light_mode`): 1 P1, 2 P2, 3 P3, 4 P4, 5 T, 6 T-all; plus
  `freeze_design` (C2) and `phi_zero`.
- **Instruments**:
  - `GorgeBody::head_for_phi`;
  - `gorge_design_mask`, `gorge_design_distance_m`, `gorge_catchment_mask`;
  - `transition_weight`, `blend_light_pass`.
- **The falls' export**: `UpscaleResult::gorge_falls` → `rivers_ll`'s falls.

## What holds

- **The construction leaves no closed hollow.** The minimal plain equals the bench's (1 cm). The bodies are paired
  1:1. F38 never fires.
- **The falls' rule** is exactly the current draw on the témoin (no soft lip; F155).
- **Protecting the design holds G-pits = 0 and G-rim = 0 in three worlds**: C2 (F143) and T-all (F156).

## What fails, with its cause

1. **The ring's clamp (G-ring), the rim the author sees.**
   - The ring cells are clamped to ≥ L(r) after the cones laid them.
   - Once the light pass no longer erodes them, they stand as an edge above the eroded slope: G-ring 99–105 (P0),
     168–245 (T, T-all), 310–588 (C2).
   - **24–31 % of its violators lie in the head falls' footprint, the majority beyond it** (F156).
2. **The light pass undoes the design.**
   - Its erosion lowers the lips and the gorges: G-pits 43 190 at ×1 under P0.
   - **Its DEPOSITION fills the design**: T, which weights only the erosion, keeps F144-P1's G-pits exactly (1 792 /
     2 330 / 2 487; F156).
   - Any hard protection leaves a wall at the design's edge: C2 4 858–6 277 edge cells, T-all 1 120–1 259, against
     P0's 500–573.
3. **The cones in the drained bowls (G-drained).**
   - `carve_diag` lays a cell from its nearest sample, even across a divide (F144-Z: 94–99 % of G-drained's violators
     laid from outside their catchment).
   - At ×1.4, 2 190 / 2 318 at the construction; 7 176–14 366 after the light pass.

## What was tried and set aside (one line each)

| candidate | round | why it was set aside |
|---|---|---|
| P1 — a floor at the construction on the design | F144 | G-pits 1 792–2 487 (the deposition); walls 4 936–6 367 |
| P2 — the erosion weighted linearly over d_t | F144 | d_t's instrument was blind (49 m); identical to P1 |
| P3 — the light pass first, the design laid after | F144 | walls up to 80°; +54 to +189 s per world (a second construction) |
| P4 — no light pass on the kept bodies' catchments | F154 | G-pits 2–18, not 0; divide walls 1 309–1 423; +8.5 s |
| C2 — the design frozen exactly | F143 / F156 | G-pits 0, but walls 4 858–6 277 at the mask's edge |
| no light pass at all (NP), with the gorge or everywhere | F155 | the author: « Inacceptable » on the hillshades (the construction's facets bare); θ 0.495; bundles × 2.5 |
| T — the true transition (the lowering weighted, d_t = 146 m) | F156 | G-pits = P1's (the deposition); walls × 2.1–2.6 P0 |
| T-all — the whole change weighted | F156 | G-pits 0 and G-rim 0, but walls × 2.0–2.4 P0, and the author sees the rim |

## The track for a resumption (hypothesis to check, not measured)

- **Protect the rim before the cones, not after.**
  - Bound the floor of every line whose cones reach the ring, so that no cone lays a ring cell under L(r).
  - Today the ring is clamped after the cones have cut into it, which makes the edge.
- The light pass would then have nothing to undo on the rim. Whether the deposition still fills the lips (T's
  lesson) is to measure first.
- **Before any build**, measure how many ring cells a cone lays under L(r), per body, and from which line.

## The queue (unchanged by the pause)

- **The breach's ramp anchored on each pit's floor** (F141 γ).
  - The spill-anchored patch is ready, NOT applied: `docs/reports/lakes_gorges/f157_close/breach_spill_anchor.patch`.
    It holds `RampAnchor`, the skeleton's explicit floor anchor, and `ALGO_BREACH` / `ALGO_CLIMATE` /
    `ALGO_HD_DRAINAGE` bumps.
  - Its doc comments still say "confirmed blind", which is false, and it was not compiled.
  - **Its confirmation criterion is to be rewritten.** The lakes also lose the pits the old ramps drained (F156).
- **The breach's second defect**: livré's coastal pits whose spill lies 0.5–3 m above the sea leave 6 294 land cells
  under it, even with the spill anchor (F155).
- **Seed 20261008002's F85 non-convergence**: the production world does not generate. The level / outlet / balance
  fixed point does not converge in 16 passes (`drainage.rs:2831`), in both breach variants (F156). The author: « Le
  laisser dans la file ».
- **Lake (basin) 1000001**: the humid below-sea basin whose outflow is stuck behind its outlet (F86–F87).

## Where things are

- **The specs**:
  - `docs/reports/c1_continental_buoyancy/f140_gorge_build/spec_gorge_age_v3.md`;
  - `…/f142_gorge_v4/spec_gorge_age_v4.md`;
  - `…/f144_gorge_v5/spec_gorge_age_v5.md`.
- **The states**:
  - `…/f145_cones/etat_au_F144.md` (to F144);
  - this page (to F157).
- **The rounds**:
  - `docs/reports/c1_continental_buoyancy/f140_*` … `f145_*`;
  - `docs/reports/lakes_gorges/f154_resume`, `f155_round2`, `f156_round3`, `f157_close`.
- **The benches**: `crates/ymir-core/tests/f126_coast.rs` (`f140_*` … `f156_t`), `f154_resume.rs`, `f155_*.rs`,
  `f156_*.rs`.
