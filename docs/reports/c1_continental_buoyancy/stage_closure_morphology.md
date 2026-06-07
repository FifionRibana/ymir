# Closure morphology validation on the rigid transport (#145 follow-up)

**Question (orthogonal to the fix):** the buoyancy fix proved continents
HOLD (global form). It does NOT answer whether the closures produce
credible *fine* morphology on the corrected transport. "Continents hold"
≠ "closures produce good morphology." This stage measures it.

**Protocol:** rigid production (`rigid_continental_crust=true`), gallery
path (single `run_with_closures`, the authoritative visual reference —
no reclassify), 64²×300 steps. **Leave-one-out (LOO)** ablation isolates
each closure's morphological contribution. Metric = SPATIAL
(`LandMorphology`: land%, perim/area, n_components, largest) on the land
mask (altitude>0, the bipolar render sea level) **+ VISUAL** (8× PNGs,
read by eye). Multi-seed {42 (561 craton cells), 1337 (0), 99 (0)} —
covers the craton-diversity flagged in FOLLOWUPS.md. EH not ablated (it
is the load-bearing regulator; disabling = known +201% runaway).

Harness: `crates/ymir-core/tests/c1_closure_morphology.rs`
(`#[ignore]`, invoke `--ignored --nocapture`). PNGs under
`closure_morphology/` (not committed, gallery convention).

## Scalar table (final state, step 300)

```
config / seed              land%  perim/A  n_comp  largest  cont% | S̃mean  S̃max
--- seed 42 (cratonic cells = 561) ---
full                        27.7   0.254     3      0.937   27.7  | 0.362  2.181
no_davis_suppe              27.7   0.254     3      0.937   27.7  | 0.151  2.121
no_erosion                  27.6   0.254     3      0.936   27.6  | 0.518  2.180
no_subduction_accretion     27.4   0.256     3      0.936   27.4  | 0.410  2.181
no_stein_stein             100.0   0.062     1      1.000   27.7  | 0.363  2.177
no_track_d                  27.4   0.256     3      0.936   27.4  | 0.427  2.181
--- seed 1337 (cratonic = 0) ---
full                        28.3   0.225     3      0.980   28.3  | 0.267  2.181
no_stein_stein             100.0   0.062     1      1.000   28.3  | ...
--- seed 99 (cratonic = 0) ---
full                        12.1   0.379     2      0.792   12.1  | 0.178  2.180
no_subduction_accretion     11.4   0.404     2      0.778   11.4  | 0.284  2.181
no_track_d                  11.4   0.404     2      0.778   11.4  | 0.288  2.181
```

### Extended seed set {2, 2026, 1988, 4138}

```
config / seed              land%  perim/A  n_comp  largest  cont% | S̃mean  S̃max
--- seed 2    (craton 714) full  36.3  0.232   4   0.997  36.3 | 0.409  2.181 ; no_DS S̃mean 0.200 ; no_SS land 100%
--- seed 2026 (craton 472) full  34.2  0.199   4   0.957  34.2 | 0.317  2.129 ; no_DS S̃max 1.791 (DS dominant) ; no_SS 100%
--- seed 1988 (craton 699) full  30.1  0.275   3   0.663  30.1 | 0.295  2.161 ; no_DS S̃max 1.934 ; no_subacc land 28.3%
--- seed 4138 (craton 202) full  21.6  0.298   6   0.763  21.6 | 0.200  2.180 ; no_subacc 19.5% ; 6 compact islands
```

Same pattern on all four: `land% == cont%`; `no_stein_stein` → 100 % land;
DS the dominant thickening source (S̃mean/S̃max drop when off, e.g. seed
2026 S̃max 2.13→1.79); subduction/Track D nibble 0.5–2 % margin land.
Visual: seed 1988 = large central continent with a well-developed brown
orogenic interior (the most pronounced DS massif of the set) + compact
neighbours; seed 4138 (21 % land, ocean-dominated) = compact rounded
islands/continents, **no filaments**. The F3 age-discontinuity ocean line
is slightly more visible on 4138 (a curved plate-boundary arc) but stays
minor. **Verdict unchanged across all 7 seeds (craton 0–714).**

## Findings

**F1 — STRUCTURAL: the land mask ≡ plate_type, not closure relief.**
`land% == cont%` to 0.1 % in *every* row, and DS / erosion / subduction /
Track D leave land% essentially unchanged (Δ ≤ 0.7 %). The land/ocean
boundary (coastline) is set entirely by **plate_type** (Continental =
land) **+ Stein-Stein** (creates ocean depth). The closures sculpt
S̃/altitude *inside* the continents; they do **not** move the coast. This
is Architecture C working as designed, and the Viz-0 "coast = plate_type"
contract confirmed on the rigid transport.

**F2 — the coast is COMPACT MASSES, not filaments.** perim/A 0.25–0.40,
largest 0.79–0.98, n_comp 2–3 across all seeds incl. craton-less. The
#141 filament failure does **not** recur on the rigid transport — the
rigid masses ARE the credible coastline (visual: seed 99 = two compact
continents at only 12 % land, no fragments).

**F3 — Stein-Stein is load-bearing for the ocean.** `no_stein_stein` →
**100 % land** (isostasy alone puts everything ≥ 0). S-S single-handedly
produces the ocean via an age-keyed depth gradient — coherent deep→shallow,
credible. Lone artifact: faint thin dark lines in the ocean = age
discontinuity at plate boundaries (minor; ties to the Track-D age-pile-up
follow-up).

**F4 — Davis-Suppe: credible orogens in S̃, marginal in the render.** In
the S̃ field DS adds brown thickening ridges roughly linear along
convergent boundaries (credible chains; S̃mean 0.36→0.15 when off proves
it is the dominant thickening source). But on the *rendered altitude* its
imprint is marginal — the thick continental interior is already brown via
isostasy. DS rides with the bounded **curtain** oscillation (blue
grid-aligned mesh on the sharp boundary) — already a registered follow-up,
bounded/cosmetic (S̃max stays ≈ 2.18 capped by equilibrium).

**F5 — erosion shapes broad relief, no dendritic valleys at 64².** Acts as
a mass sink (S̃mean 0.36→0.52 when off). At the 64² tectonic grid there
are NO resolvable fluvial valleys — fine dendritic morphology is a later
upscale+detail-phase product (4096²), not a C1 feature. Broad relief
credible; the deposition half-fix produces no visible artifact.

**F6 — subduction/accretion + Track D: small, clean margin contribution.**
Nibble ~0.5–0.7 % margin land (full→no_track_d), **no fingers** (the ≥2
convergent-neighbour fix holds — no 1 px false-land spikes visible).
Topology quasi-static in a gallery run (merge/split rarely fire; matches
the #141 registered Phase-3 prerequisite). Credible margins.

**F7 — Track D trajectory credible.** seed 42 cycles 0→300: overthick
white plateau (init) settles to differentiated green/brown relief via
equilibrium-height. A credible thickness-relaxation trajectory; footprint
near-static (quasi-static topology).

## Verdict — closures SOUND on the rigid transport → LIGHT re-validation confirmed

No closure produces doubtful/artefact morphology that blocks downstream
work. **Phase 3 / piste 4 may proceed.** Two cosmetic artefacts, both
bounded/minor and already registered: the DS curtain (F4) and the S-S
age-discontinuity ocean line (F3).

**Caveat to carry forward (not a defect):** the closures' contribution to
the *rendered* terrain is smaller than one might expect — the C1 64²
gallery render is dominated by plate_type + isostasy + Stein-Stein. Rich
terrain morphology (orogen chains, valleys) lives in S̃ and is only weakly
expressed in altitude at 64². By architecture, fine morphology is the
later upscale+detail+erosion phase's job at 4096². C1's job is the
plate-scale land/ocean geography + broad relief — which is credible.
