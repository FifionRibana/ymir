# Finding 177 — closing the accretion: the guard has no stale cache; per A1 the active share is 27.4 % and every continent has an active margin, but the active margins get no Andean uplift — the B2 stop fires

**Status: measured (no closure modified, no solver); the declared B2 stop fired after tier 1, so tiers 2–4 (B3's attribution,
the formations, the unclipped transect) were not run; nothing adopted, nothing committed without the feu vert.** F176 stays
uncommitted (no feu vert in this brief). Production unchanged; no cache key, no `ALGO_*` touched.

- Declared before any measurement and before re-reading F176's data: `f177_declared.md` (A1, A2, the unit convention, the
  author's criterion verbatim, the rules). Predictions: `f177_predictions.md`.
- The bench: `crates/ymir-core/tests/f177_cascade.rs` (F176's helpers + F177's observer); Part 0's test `f177_fresh_guard` in
  `crates/ymir-viz/src/bridge/c1/hd.rs` (test-only). Each world's lines in `results/`, assembled in `f177_bench_output.txt`.
- **Units**: world km (the 400 km world), the terrestrial equivalent (× 7.5) in parentheses; a C1 cell is 6.25 km (46.9 km).

## The budget (announced: ~20 min wall)

| tier | content | wall |
|---|---|---|
| 0 | « livré » recomputed at 8 192² from an empty cache | 544 s |
| 1 | C1: 15 seeds under F176's selection + the panel at the fixed continent (A1, B1, B2, B4, B5), parallel | ~3 min |
| — | B2's mechanism (C1, 15 seeds, parallel) | ~1 min |
| 2–4 | not run (the stop) | — |

## Part 0 — the guard: there is no stale cache

- **The code**: the guard looks its reference up by digest (`bench_guard::check_against`, `check_lakes_against`) and returns
  `NoReference` when no entry has the run's digest; `run_hd`'s digest chains `tectonic_key`, which folds `ALGO_TECTONICS`. So
  F175's bump gave `NoReference` on every state, and `f123_viz_guard` prints « ≠ / non gardé » for `NoReference` and `Mismatch`
  alike (`workspace.rs:7045`); the panic line with the variant was filtered out by the checks' grep.
- **The test**: « livré » recomputed at unchanged digests from an empty cache (the test binary run from an empty directory):
  **field `Match`, lakes `Match`** (`results/part0_fresh_guard.txt`, 544 s).
- **The cause**: F175 read the guard's `NoReference` as a mismatch. **There is no stale downstream cache, and no commit to
  blame.** F175's « all six fields and lakes different » is withdrawn.
- **The proposal** (nothing changed without the author's feu vert; the references untouched): `f123_viz_guard` prints the
  variant (« non gardé » apart from « ≠ »), and the checks keep the panic line; the five other states fresh (~45 min) only if
  the author wants them.

## A2 — the terrestrial reference (comparison, not a target)

- **Leat & Larter 2003** (p. 1): intra-oceanic arcs « comprise some 40% of the subduction margins of the Earth »; « Some 17 000
  km, or nearly 40%, of the global length of volcanic arc is not situated on the continental margins » out of « about 43 500 km
  (von Huene & Scholl 1991) » (39 %).
- **Bradley 2008** (p. 2): passive margins 105,000 km, « the convergent plate boundaries (53,000 km) » — **no source and no
  definition** given. **Bird 2003** (Table 3, p. 46): subduction zones 51,310 km; continental convergent 23,003; oceanic
  convergent 17,449. Bradley's total is 3 % from Bird's subduction zones and 42 % short of all three → **the reading (an
  inference): collisions are not included.**
- **The recalculation**: 0.6 × 53,000 = 31,800 km of active continental margins → **23.2 %** of the continental margins; with
  Leat & Larter's own totals 20.2 %, with Bird's subduction zones 22.7 % → **~20–23 %**.
- **Variant B closed by the sources**: the Izu-Ogasawara arc crust is « about 20 km » (Nishizawa et al. 2006, p. e33) and
  intra-oceanic arcs are « typically mostly submerged » (Leat & Larter p. 2); under F167's Airy, 20 km is −2.6 km (emergence
  needs 32 km). F176's « no arc island » agrees with them.

## B1 — the active share per A1 (15 seeds, F176's selection)

| seed | 8 | 18 | 19 | 23 | 29 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 9 | 10 | 11 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A1 share, % | 9.0 | 37.4 | 16.4 | 28.3 | 27.4 | 13.5 | 21.3 | 28.1 | 31.5 | 15.0 | 25.2 | 54.4 | 34.8 | 23.2 | 45.6 |

- **Median 27.4 %** (9.0–54.4 %); ≥ 10 % on 14; above 50 % on one (seed 7). Bradley-Leat's ~20–23 % lies inside the spread.
- **Every continent (15 / 15) has an active margin** (≥ 1 group of ≥ 3 active coast cells); 15 / 15 also by F176's « newly
  active » reading.
- **The coast → trench distance** (the active walks, 743 cells): median **3 cells = 18.75 km world (141 km terrestrial)**; 1–8
  cells (6–50 km world, 47–375 km terrestrial); the histogram 28 / 65 / 310 / 236 / 73 / 24 / 4 / 3 for 1…8 cells. The 3 × 3
  cases (the coast on a boundary): 110 of 853 active cells.
- **Mixed ties**: 288 of the 743 walk-active cells reach, at the first depth, both an active and a passive boundary cell (the
  declared « any » rule calls them active).
- Maps: `images/a1_panel_selection.png` (red the active walks, dark red the 3 × 3 actives, green passive, black the trench cells
  reached). The trenches lie inside the seaward bands, ~3 cells off the coasts.

## B2 — the uplift at the coasts, and the stop

The band: the continental cells within Davis-Suppe's reach (30 cells at 64²) of a coast cell, each assigned to its nearest one.

| seed | 8 | 18 | 19 | 23 | 29 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 9 | 10 | 11 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Δh_src active / passive | 1.02 | 1.72 | 1.57 | 1.31 | 0.49 | 2.57 | 0.70 | 1.06 | 2.42 | 1.32 | 2.62 | 0.70 | 1.51 | 0.38 | 1.46 |
| active band: Andean (C-O), m | 817 | 1188 | 145 | **0** | **0** | 2206 | 783 | 157 | 468 | 1018 | **0** | **0** | 130 | 24 | 12 |
| active band: collision (C-C), m | 1330 | 1305 | 1186 | 1238 | 978 | 920 | 491 | 1552 | 2524 | 2412 | 2261 | 930 | 2410 | 423 | 1929 |
| active band: O-O arcs, m | 67 | 4 | 0 | 4 | 3 | 47 | 14 | 11 | 0 | 2 | 40 | 4 | 12 | 3 | 13 |

- **The active band's uplift comes from collision** on every seed. **The continent-ocean wedge gives nothing on 4 seeds (23, 29,
  6, 7)** and ≤ 157 m on 5 more (19, 3, 9, 10, 11); the O-O arcs give ≤ 67 m.
- **The declared rule fires** (« a defect … if the active band has no continent-ocean wedge uplift (mean ≤ 0) on a seed with
  active coasts: the round stops »). **The round stops here.**
- Active ≥ 2 × passive on 3 / 15 seeds (1, 4, 6), on 0 of the 5 panel seeds under the selection (1 / 5 at the fixed continent).

### The mechanism in the code

Of the 853 active coast cells, 110 sit on a boundary (3 × 3), 403 have their trench within the arc's reach (≤ 3 cells), and
**340 beyond it**; only 198 received any continent-ocean wedge or arc.
1. **Davis-Suppe's O-C seed** (`oc_override_seed_mask`, `boundary_classification.rs:383`) needs a convergent **continental** cell
   with an oceanic neighbour of another plate: the continent must touch the subducting plate. With the seaward margin, the
   continent's own oceanic band lies between, so the convergent boundary is ocean-ocean.
2. **Option 5's A7** (`clear_oo_seeds`, `:420`) cuts Davis-Suppe's O-O seeds (it removed F175's dome), so the overriding band
   gets no wedge either.
3. **The O-O arc** (`distribute_arc_mass_plate`, `subduction/source_term.rs:534`) goes to the overriding plate's continental cells
   within `arc_distance` = 3 cells of the consumed cell (`subduction/params.rs:73`), else to its ocean; with the trench at 3–8
   cells most of it lands on the band, and its volume (half the consumed mass) is small.

So an A1-active coast is Andean by kinematics — its plate overrides an active subduction — but no relief closure acts on it: the
relief closures act only at the plate boundary.

### The physical source of a correction (for the author's next finding; nothing changed)

**Hypotheses to verify — the sources are not in `docs/refs`**:
- Volcanic arcs stand where the slab reaches ~100 km depth, typically 150–300 km behind the trench (the arc-trench gap:
  England, Engdahl & Thatcher 2004; Syracuse & Abers 2006; Gill 1981).
- Andean crustal shortening and thickening happen in the overriding plate behind the arc (Isacks 1988; Oncken et al. 2006).
- Our median trench distance is 141 km terrestrial, inside that gap, so the arc front and the wedge would sit at or behind the
  coast.

A possible correction, sourced by these papers once supplied: seed the O-C wedge on the overriding plate's continental cells
within the arc-trench gap of an active trench (not only on a continental cell that touches the subducting plate), and place the
arc at the arc-trench distance rather than within 3 cells of the consumed cell. These are closure changes, for the next finding.

## B4 — the channels (at the end of C1, the cascade's coast)

- **15 seeds (the selection)**:
  - 252 channels ≤ 2 cells wide: rift 157 (62 %), band remnant 84 (33 %), trapped by a merge 6 (2 %), other 5;
  - **one-cell-wide groups 278: rift 202 (73 %)**, band remnant 57 (21 %), trapped 4 (1 %), other 15 (5 %).
- **The priority straits** (the panel at the fixed continent, F176's 256² young montage, `images/channels_panel_fixed.png`):

| seed | channel (x, y) | size | class |
|---|---|---|---|
| 29 | (39, 16) | 62 cells, 194 km world (1 455 km terrestrial) | rift |
| 18 | (49, 37) | 36 cells, 106 km | rift |
| 23 | (10, 28) | 40 cells, 106 km | rift |
| 19 | (53, 20) and (54, 13) | 17 + 8 cells | rift |

  Seed 19's sinuous strait also has an « other » piece at its north end ((47, 25), 8 cells) and two band remnants. The strait
  that separates the white highlands from the continent is a rift on 4 / 4 (its relief's origin was B3's question, not run).
- **Their mechanism**: the split (`apply_rifting_split`, `rifting/split.rs:212`) makes **the rift strip — the continental cells
  of plate a 4-adjacent to plate b (`find_rift_strip`, `:291`), one cell wide — a new plate**. The strip and its neighbours keep
  thinning at a divergent boundary (`apply_rifting_thinning`, `rifting/source_term.rs:57`) until they go below sea level. The
  channel's width is the strip's single cell, by construction.

## B5 — the plates

- 15 worlds: **212 plates** (10–24 per world): < 1 % 82 (39 %), 1–5 % 81 (38 %), > 5 % 49 (23 %).
- **157 merges and 137 splits** (0–23 splits per world; seed 7: 23 splits, 24 plates).
- Bird 2003 (Table 1, p. 6; 52 plates): 39 (75 %) / 6 (12 %) / 7 (13 %). Our worlds have fewer small plates, more mid-size
  ones and more large ones in proportion. Not a target.

## B6 — the rifting (inventory)

- **The trigger**: a negative source on **continental cells classified divergent** (`rifting/source_term.rs:57`, `:82`), Δs =
  `thinning_rate` (1.0) × |v_rel · n̂| × dt.
- **The split**: `rifting/split.rs:212`. It fires when a pair has diverged for ≥ 75 consecutive steps (`DivergenceTracker`,
  `:94`) and the strip's minimum s is < 0.7. The new plate is the one-cell strip, with its age reset and the parent's velocity
  + 0.005 perpendicular.
- **The sources**: McKenzie 1978 for β = 1.4 (s 0.7); the rates are a « Stage E2 W7 analytical first-pass » (`params.rs`
  docs).
- **Was it set assuming the merges compensate? Not stated anywhere.** The Track D report declares the splits « rare by design
  (0–3 / seed at 300 steps) » (README §4) beside merges that collapse 8 plates to 1–2 (« Pangaea narrative accepted », §3). The
  split threshold (75) was set above the merge's (50) « because rifting needs MORE sustained extension than collision needs
  sustained compression » (`params.rs:17`). So the plate count's balance did rest on the merges in practice. Here the splits
  are 9 per world on average, not rare.

## B7 — the transect of seed 18 (from the code; not rerun)

`transect_png` clamps the altitude to [−5 500, 2 500] m (`f176_cascade.rs:1730`, `v.clamp(-5500.0, 2500.0)`): **the flat top
at 2.5 km is the plot's clip.** The unclipped image (tier 4, `transect_png_free`, ~1.5 min) was not run because of the stop.

## The rules

- **The stop fired (B2)**; tiers 2–4 were not run (B3's attribution maps, the formations and condition (v), the unclipped
  transect). They are written and resumable (~10 min) if the author wants them despite the stop.
- **Adoption: not proposed.**
  - (i) the A1 median 27.4 % ✓;
  - (ii) 100 % ✓;
  - (iii) 15 one-cell groups are « other » (no class) ✗;
  - (iv) B3's maps were not produced;
  - (v) not measured, and B2 shows that the continent-ocean subduction builds no relief on the active margins of 4 seeds.

## Predictions

**The reviewer's**:

| | prediction | outcome |
|---|---|---|
| T1 | A1 median in 15–35 % | held (27.4 %) |
| T2 | median coast → trench distance ≥ 2 cells | held (3 cells) |
| T3 | collision gives more than half of the mountain gain on ≥ 3 panel seeds | not measured (stop). B2 shows collision dominates the active bands' uplift on every seed. |
| T4 | active ≥ 2 × passive uplift on ≥ 4 panel seeds | refuted (0 / 5 under the selection, 1 / 5 at the fixed continent) |
| T5 | ≥ half of the one-cell channels trapped by merges rather than rifts | refuted (rift 73 %, trapped 1 %) |
| T6 | the flat top is the plot's clip | held, from the code |
| T7 | Andean cordilleras on ≥ 5 seeds, the separated highland never alone | not measured |

**Mine**:
- held:
  - Part 0 `Match`;
  - T1's range (20–35 %), T2 (2–3 cells), 15 / 15 continents;
  - T4 refuted; T5 refuted; T6;
  - ≥ 30 % of the plates below 1 %; Bird's 7 / 6 / 39;
  - B6 (not stated; the merges carried the balance).
- refuted:
  - « no seed above 50 % » (seed 7: 54.4 %);
  - **« B2's defect not found »** (found on 4 seeds);
  - « band remnants and O-O reassignments ≥ 60 % of the one-cell channels » (rifts 73 %);
  - « 1–3 plates above 5 % » (2–6).

## Open, for the author (nothing chosen)

- (a) **B2's defect**: the relief closures only act at a plate boundary, so an Andean margin with its trench offshore (in the
  continent's own band) gets no wedge. The correction is the author's to write (the arc-trench gap; the sources to supply).
- (b) The rift strip is one cell wide by construction. Splits run at ~9 per world against the « rare by design » note. Do the
  channels it makes count as physical (a rift sea, 47 km terrestrial) or as an artefact of the strip rule?
- (c) The guard: print the variant. Should the five other states be run fresh (~45 min)?
- (d) Tiers 2–4 (B3, formations, transect, ~10 min) on request.
- (e) F176's commit (no feu vert in this brief).

**At commit (2026-10-18), decided**: the author's words: « Feu vert pour commiter le F176 puis le F177 (config.rs et workspace.rs restent dehors ; le test ignoré de hd.rs peut entrer). »
- **F177 is committed as measures only**: the `f177_cascade` bench, and the ignored test `f177_fresh_guard` in
  `crates/ymir-viz/src/bridge/c1/hd.rs`, which is test-only. No closure was modified; production is unchanged and the guard
  stays « banc » 6/6.
- **The reviewer's errors, recorded (the author's)**:
  1. **F174 B′**: « la bande vers le large a été conçue sans vérifier la condition de la fermeture continent-océan (le continent
     doit toucher la plaque qui plonge). Avec la fosse au large, cette fermeture ne s'allume plus : c'est la cause du défaut
     B2 ».
  2. **F175**: « l'affichage « ≠ / non gardé » a été lu comme un cache périmé, sans demander le verdict exact de la garde ».
- **Question (b), the author's answer (verbatim)**: « la largeur d'une cellule (47 km terrestres) est plausible pour un fossé de
  rift à son ouverture. En revanche, la règle « le ruban devient une plaque » est suspecte. Elle est mise en file d'attente pour
  un finding à part, après les cordillères. »
- **Queued, in this order**:
  - the cordilleras: B2's defect, the next finding, written by the author;
  - then the rift rule « the strip becomes a plate » (`find_rift_strip`), in a finding of its own;
  - deposition stays F178's subject as announced in F177's brief.
- **Still open**: the guard's printout (show the exact verdict); the five other states recomputed fresh (~45 min); tiers 2–4
  of F177 (~10 min) on request.
