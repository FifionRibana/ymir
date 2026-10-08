# F157 — instruments, declared BEFORE the measurements (2026-10-08)

**The states**:
- **OFF** = the production témoin, C2 /10 col (gorge off).
- **ON** = OFF with `lake_base: Some(LakeBase::InputLakesAndBasins)` (the extended lake base, F133), gorge off.

## C — the gorge's closing summary

- `etat_gorge_pause.md` in this folder.
- An ADR closing entry.
- `accepted_defects.md` only if a defect is accepted. The gorge is paused, not accepted, so none is expected.

## V — the lakes' outline (viz)

- On Relief → Ombrage, the outline of `run_hd`'s final lake mask (`drainage.lake_map`, after the crater pass):
  - a cyan line on the lake cells with a 4-neighbour outside their lake;
  - a dark navy line on the shore cells 4-adjacent to a lake.
- A « Contour des lacs » box in the relief menu (on by default) hides it.
- A view only. Permanent test `the_outline_rings_each_lake_inside_and_out`, with its negative control.
- The guard must stay 6 / 6.

## B — the lake base alone

**B1, F133's table replayed** (bench `f135_w3m`, unchanged, on the current code):
- OFF / F132 / extended: canyons, coast, R8 terrain, relief p50, the Δz between flat resolutions, the lakes, and θ
  with its CI.
- Compared with F133 and F135-W (26 / 34) **to the digit**. Every difference is attributed before going on.
- **The cost**: `build_world` OFF and ON, plus the run_hd tail (`viz_hd_lakes_on`), timed in `f157_b`.

**B2, the crops** (viz bench `f157_viz`, on `run_hd` itself, the guard's production `HdParams`):
- **The new lakes** (F133v's rule): the ON final lakes with < 50 % of their cells under any OFF lake. Each crop is
  615 × 615 cells centred on the lake's cells' centroid, in viz coordinates (data rows, south-first).
- **The control crop**:
  - the 615² window (stride 64) with zero conditioned cells changed OFF → ON and the most land;
  - **if none exists**, the window whose cells all have **|Δz| ≤ 1 m (the declared threshold)**, with the most land;
  - if none either, it is said.
- **Per crop**: the lakes OFF / ON with ≥ 1 cell in the crop.

**B3, the images**:
- `layer_color_image(Relief, Shade { outline: true })`, the viz's own buffer, cropped, north up, OFF and ON;
- one PNG per crop and state.

**B4, the steps below the based lakes** (`f157_b`, ON, F136's amended instrument, completed). For every ON final lake
≥ 1 km² that is not a crater lake and has a D8 receiver at `Lake::outlet`, along the ON D8 path from the col:
- **the drop**: the lake's level minus the first graded reach ≥ 2 km (slope over 2 km < 2 %) within 10 km (K6b);
- **the length**: the distance to that reach;
- **the mean slope**: atan(drop / length);
- **the max slope over 100 m**: the steepest 100 m window between the col and that reach;
- the crop holding the lake's centroid.
- The steps > 200 m are listed.

**B5, the Living Landz exports**:
- `run_hd` with `export_dir`, the production `HdParams` (the guard's C2 /10 col literal), OFF and ON;
- into `exports/f157/off` and `exports/f157/on` (repo root; `exports/` is git-ignored).

**B6, what the promotion would change** (without making it):
- the production field's hash (ON against the guard's OFF entry);
- the guard entries that would move;
- the lake references;
- the exports;
- the cost per world (B1's timings).

**B7, the production defects** (OFF; the 11 crops plus the control):
- **comb teeth** (F124-P4's definition): watercourse segments of the final OFF network with ≥ 80 % of their cells on
  the construction's wall cells (`carve(S1)`'s carved and not floor);
- **planar walls** (F155): the land share at 28° ± 0.5°;
- **axis alignment** (F121 / F126):
  - the terrain R8 (`aniso`, 16) on the crop;
  - the network R8 (c8 of 32-cell chords) of the watercourse segments with an area of 1–3 km².
- **The stage creating each**, measured on the same crops at each stage:
  - **S1** (tectonics + FBM, no incision, no erosion);
  - **the construction** (`carve(S1)`);
  - **C1 bare** (the construction plus everything after, without the light pass);
  - **OFF** (with the light pass).
  - Planar walls and terrain R8 on each. The teeth and the network R8 on C1 bare's and OFF's networks.
- **No remedy.** An inventory for a possible later work.

## The author's questions (copied as given into the report, for the look)

Per crop, in the viz then in Living Landz:
1. Amont : la vallée arrive-t-elle au lac à hauteur de sa surface ?
2. Aval : en sortie de lac, la descente te semble-t-elle acceptable (gorge, rapides), ou vois-tu un mur ou une marche
   artificielle ?
3. Non-régression : vois-tu en ON un défaut absent en OFF ?
4. Les lacs ajoutés te conviennent-ils ?

## Amendment after `f157_b` (NOT blind): B7 split out for memory

- **Why**: `f157_b` measured the cost and B4, then aborted in B7 on a failed 512 MB allocation. It held the OFF world,
  its tail and C1 bare's at once.
- **The change**: B7 runs in its own bench (`tests/f157_b7.rs`), one stage at a time, each freed before the next. The
  crops' origins are those `f157_viz` printed (the same rule, the same 12 crops). The instruments are unchanged.

