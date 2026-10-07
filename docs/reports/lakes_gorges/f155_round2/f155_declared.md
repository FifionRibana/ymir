# F155 — instruments, declared BEFORE the measurements (2026-10-07)

**The worlds**:
- **The six guard states**, with `f134_lake_guard`'s Knobs (the viz's origin; 45° / 40°): livré, A1+B2, C1 nue,
  C2 /10 col (the témoin), C2 /3, C2 /10 niveau mer.
- **ON** = the témoin with `LakeBase::InputLakesAndBasins`.
- **The gorge** = `GorgeRetreat::v5` (or v6 below) on `on_at(age)`.
- Benches:
  - `f155_br` and `f155_f` in `crates/ymir-core/tests/f155_round2.rs`;
  - `f155_p` in `crates/ymir-core/tests/f126_coast.rs` (it needs that file's θ, coast and gate helpers).

## Br — the breach, blind confirmation then production

**The copy is F154's amendment 2, untouched**: the ramp's anchor is `max(height[nb], filled[nb]) − EPS`. Per state,
production then the copy, through the same run_hd tail (`viz_hd_lakes_on` / `viz_hd_lakes_with`):
- **the land cells the breach leaves ≤ the sea** (eroded > sea, conditioned ≤ sea): their count, and those NEW under
  the copy (not under production);
- **the lakes**:
  - the count;
  - `lake_fingerprint`'s two halves: **the mask's hash = the footprint (the first half)**, and the list's hash
    (levels and areas, the second half);
  - production's fingerprint is also checked against the guard's stored `lakes_hash` (the bench is the guard);
- **the conditioned cells** raised / lowered by > 1 cm;
- **the river segments**: the count, and its variation.

**The stop rule (the brief's, verbatim)**: « si un état change un lac (compte ou empreinte), ou crée une cellule de
terre sous la mer qui n'existait pas, ne pas passer en production. Rapporter et localiser. »
- "A lake changes" = the count or the first half of the fingerprint differs.
- "A new below-sea land cell" = a land cell ≤ sea under the copy and not under production.

**Production, if the rule does not fire**:
- **The scope**: ONLY the HD conditioning breach, the one `run_hd` caches as `conditioned` (`hd.rs`), and its bench
  mirror `common::viz_hd_lakes_on`.
  - It goes through `terrain::flow::breach_monotone_anchored(.., RampAnchor::Spill)`.
  - `breach_monotone` / `breach_monotone_protected` keep the floor anchor. **The skeleton's breach
    (`valley_construction.rs:642`, F120's A0 instrument) is NOT changed**: changing it would move the construction
    and the eroded field. It is a separate decision, named for the author.
- **The cache**:
  - `ALGO_BREACH` 1 → 2 (the conditioned field);
  - `ALGO_CLIMATE` 1 → 2 (the climate reads the conditioned field but is keyed on the eroded key);
  - `ALGO_HD_DRAINAGE` 12 → 13 (the bundle reads the conditioned field and the climate, keyed on the eroded key).
  - `ALGO_BREACH` alone would serve a stale climate and stale lakes and rivers.
- **The permanent test** (rule 13): a row of pits beside a lake, behind a long plateau just above the sea. The spill
  anchor sends no ramp below the sea level. Negative control: the floor anchor does.
- **The guard**:
  - **the field guard must stay identical** (the eroded field: no ALGO in its key moves, and the breach is after
    it);
  - **the lake guard's hashes must stay identical** where the confirmation holds. Its digests change (the bumped
    `ALGO_HD_DRAINAGE`), so it is regenerated (`f134_lake_guard`), and the diff must show the six digests only;
  - the conditioned field and the rivers change: they are not guarded.
- The ADR's queue entry for the defect (F142, `etat_au_F144.md` defect 3) is marked resolved.

## F — the soft-lip rule (the gorge, gated)

- **The rule**: `GorgeRetreat::soft_lip`, in `GorgeRetreat::v6` = v5 + `soft_lip`. When on, φ = 0 for a kept body
  whose col lies on a **soft** substratum cell (hardness 0, F154's layer); otherwise φ is the current draw.
- **The soft cells**:
  - `rocks::build_substratum`, read on the construction's input, at hardness 0;
  - its sources are C1's rift mask (`derive_tectonic_labels`) and the edifices;
  - the craton and belt masks are passed empty. Both are hard and are tested AFTER the volcanic and rift tests, so
    the soft cells do not depend on them.
  - They go to the skeleton through `skeleton_with_soft_lips`; `skeleton` passes `None`.
- **Permanent test**: on F154's synthetic world, a soft col gives φ = 0. Negative control: a hard col keeps the draw.
- **The témoin check** (bench `f155_f`), v5 against v6 at ×1 p .5:
  - every φ, the tagged falls and the constructed field (FNV-1a), which must be identical;
  - the soft mask is also compared with `tectonic_sources(..).substratum(S1)`'s hardness-0 cells.

## P — what the light pass brings (measurement only)

**Grep first** (rules 11 / 11c): "light pass", "light incision", "light_k_time_fraction", "k_time", "Finding 121",
"texture", "A1+B2". The original reason is quoted before any number.

**The worlds**: the témoin (C2 /10 col, gorge off) WITH the light pass (`Some(0.1)`) and WITHOUT it (`None`, the
same knobs otherwise). Measured after the Br production change, if it ships.

**The instruments**, on the final world (the eroded field, the run_hd tail):
- **relief p50**: F95's paired relief (`f95_criteria`);
- **σ local p50 / p90**: 3 × 3 in metres, stride 37 (F121);
- **R8**:
  - terrain (land, 16 bins);
  - network (the watercourse segments' 32-cell chords);
  - the comb tile's terrain R8 (F121's `COMB_TILE`);
- **θ**: the trunk links on the construction's carved mask, against OFF ×1's eroded CI widened by 0.005 (F154);
- **canyons**: `over_dug_depression` (F144's extras);
- **the parallel bundles**: F146's instrument; pairs of unjoined `rivers_ll` rivers running within 2 cells for > 2 km;
- **rivers**: the segment count, `rivers_ll`'s count, and the Strahler ≥ 2 count;
- **lakes**: the count and the total area;
- **biomes**: the share of land cells whose biome differs (`c1_biomes_classified_wet`);
- **the cost**: `build_world`'s time, twice per world (single timings carry ±30 s);
- **what appears without it**:
  - **planar walls**: the land share at 28° ± 0.5°, and the wall slope p50 on the construction's wall mask;
  - **sharp crests**: a land cell higher than both neighbours along one of the 4 axes, with both sides steeper than
    28°;
  - **combs**: the comb tile's R8 above.

**The images**:
- the hillshade of `run_hd`'s Relief → Ombrage (azimuth 315°, altitude 45°, flat = 1), north up, 512 × 512 cells,
  with and without the pass;
- **three crops, chosen on the WITH world before the WITHOUT world is built**, by these rules:
  - **montagne**: the 512² tile (stride 256) with the highest p90 altitude;
  - **plateau**: among tiles whose land p50 lies above the land's p75, the one with the lowest σ p50;
  - **vallée construite**: the tile with the most construction wall cells.

**The gorge without the pass**:
- v5 (soft_lip off) at ×1 p .5, ×1.4 p .5 and ×1.4 p 0;
- WITH the pass (P0) and WITHOUT it (NP), with F144's gates (F154's code, the same definitions and tolerances);
- the cost per world.

**The two options, costed, none chosen**:
- **(a)** no light pass when the gorge is on;
- **(b)** no light pass at all.
- For each: the seconds per world against 249.8 s, and the gates and production measures above.

## Amendment after the Br run (NOT blind): the stop rule's localisation

- **Why**: the stop rule fired on « livré » (the lake mask's hash differs). The brief asks: « Rapporter et
  localiser. »
- **The instrument** (bench `f155_br_loc` in `f155_soft.rs`, on the states the rule fired on):
  - **the lakes**: `lake_listing` production against spill, line by line (id, type, level, area, cells, the mask's
    hash), and per lake id the cells lost and gained;
  - **the land cells each anchor leaves below the sea**: the count, the share taken by a ramp, the ramp's anchor
    height above the sea and its length (p0 / p50 / p90), the original heights, and the 1024² tiles holding most.
- **Production stays as it is**: the breach is not changed (the stop rule).
- **P is measured on the current production** (the floor anchor), as is.

