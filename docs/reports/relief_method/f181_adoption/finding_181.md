# Finding 181 — the adoption of the active margins in the continent profile: the rules hold at the intermediate and old ages, but the cordilleras lower and narrow with age, so the plates go to the author before the adoption is made effective

**Status:**
- Part B and Part C measured. The adoption rule holds at both ages: (v), no rampart, the circumnavigation.
- **The cordilleras' aspect changes with age** (lower, narrower, the mountain share falling), so, as the author's complement
  asks, **the plates are shown to the author before Part A is made effective**.
- Part A is prepared, not applied: no default changed, no code changed outside the bench. Production unchanged.
- Nothing committed without the feu vert.

- Declared before any measurement: `f181_declared.md`. Predictions: `f181_predictions.md`.
- The bench: `crates/ymir-core/tests/f181_cascade.rs` (F180's helpers). Results in `results/`, assembled in
  `f181_bench_output.txt`.
- Budget: announced ~38 min wall; actual: équerres ~1 min, tier 1 1 453 s, the young widths ~8 min, tier 2 78 s.

## The author's words (condition v, verbatim)

> « Non rien ne semble me gêner pour l'instant. Les graines montrent de la variabilité. »

They concern F180's diagnostic maps (the young profile). « Pour l'instant »: the validation holds for the young profile; the
adoption depends on Part B.

## The decision recorded

- F180's coverage rule fired, and F180 adopted nothing.
- The adoption is decided here, for two reasons:
  - that rule measured the presence of the Andean wedge, not the rampart (the reviewer's error, recorded at F180);
  - the direct measures (W₁₀, the active-passive contrast, the formations) pass.
- It becomes effective only if Part B holds, and, the aspect having changed, after the author's look at the plates below.

## Part B — the ages (256², the 15 seeds of F176's selection; F176 the reference)

| | F176 interm. | F176 old | F180 young (F180) | **F180 interm.** | **F180 old** |
|---|---|---|---|---|---|
| seeds with ≥ 1 Andean cordillera | 9 | 9 | 15 | **15** | **15** |
| collision the only formation | 0 | **1** | 0 | **0** | **0** |
| W₁₀ panel, global (max) | 0.31 (seed 18) | 0.27 | 0.28 | **0.21** | **0.19** |
| active ≥ 2 × passive (W₁₀), panel | 2 | 2 | 4 | 4 | 4 |
| mountain %, median (range) | 47.4 (7.5–79.6) | 40.0 (0.3–75.8) | ~37 (25–55) | **23.1 (10.8–47.7)** | **15.6 (5.8–40.5)** |
| median altitude | 303 m | 202 m | — | 115 m | 53 m |
| coast dimension, median | 1.136 | 1.140 | — | 1.180 | 1.181 |
| cordillera width, median of the seeds | 9.0 km world (68 terr.) | 4.8 (36) | **42.4 (318)** | **38.8 (291)** | **25.0 (188)** |
| 256² circumnavigation lost | 0 | 0 | — | 1 | 1 (seed 23) |

- **The adoption rule holds at both ages**:
  - (v): 15 / 15 seeds with a cordillera, collision never alone;
  - no rampart: the panel's global W₁₀ is ≤ 0.21, far under 0.43;
  - the circumnavigation is lost on 1 seed (fewer than 2).
- **The aspect changes**:
  - the mountain share falls (F180: young ~37 % → intermediate 23 % → old 16 %);
  - the cordilleras narrow (width 42 → 39 → 25 km world, 318 → 291 → 188 km terrestrial; per seed at old 0.6–95 km world);
  - their heights fall (2 048², seed 7, below).
  - The active coasts stay ≥ 2 × steeper than the passive ones on 4 panel seeds at every age (active W₁₀ 0.35–0.45
    intermediate, 0.37–0.40 old on seeds 8, 23, 29).
- **Against F176 at the same ages**, F180 keeps cordilleras on 15 / 15 (F176: 9). F176 at old has collision as the only
  formation on 1 seed; F180 on none.

## Part B — 2 048², old

- **Seed 7** (also the longest cordillera at the old age, so its whole map): the transect spans −5 025 … 1 903 m. The active
  coast peaks at **1 903 m** (young: 4 693 m), the passive coast at 83 m; the world's maximum is 4 675 m (young: 5 964 m).
- `images/transect_F180_old_seed7.png`, `images/map2048_F180_old_seed7.png`.

## Part C — the « équerres » (diagnostic, no correction)

- **The majority-active coast cells at sea** (their C1 cell and its 8 neighbours sea): **190 of 644 (29.5 %)**, in 75 groups on
  15 seeds.
- **Two kinds**:
  - **the right-angle lines** (seeds 18 and 2, `images/equerres_15.png`, red around pink): they border **blocks of cells re-typed
    continental during the run** (57 and 64 cells; s at the end ≈ 0.001). The closure is the continent-ocean floor-triggered
    reassignment, **`subduction/source_term.rs:266–268`**: a consumed oceanic cell below s = 0.05 with ≥ 2 convergent
    continental neighbours joins the continental plate and is re-typed continental, with almost no crust. Its edges follow the
    consumed plate's straight boundaries. They are 18 % of the offshore active cells, about 5 % of the active coast cells;
  - **the scattered offshore cells** (82 %): the seaward band's continental cells (F174 B′, continental while the crust exceeds
    10 km), at sea with s ≈ 0.2 at the end.
- Not the 4-connected rift strip. Nothing corrected.

## Part A — prepared, not applied

When the author confirms after the plates:
- the continent profile's defaults take `Margin::seaward_oo_collision_f180()`: the seaward margin, option 5, the merge at
  collision, the active-coast seed at the arc (`andean_active = (4, 10)`) and Davis-Suppe's mass bound (`ds_mass_bound = 1`);
- the viz boxes are checked for the continent profile;
- the old controls keep their flags explicitly off.
- **The caches**: the continent profile runs only in the cascade (`cascade/history.rs`, `cascade/profile.rs`), which no
  `ALGO_*` keys. Production's `run_hd` does not use it (F179 Part 0), so no `ALGO_*` needs a bump. The guard (production) is
  untouched, which the checks will verify.

## Predictions

**The reviewer's**:

| | prediction | outcome |
|---|---|---|
| E1 | ≥ 10 / 15 keep a cordillera at the old age | held (15) |
| E2 | old W₁₀ < young on ≥ 4 panel seeds, under 0.30 | held (5 / 5 lower: 0.13 / 0.12 / 0.03 / 0.19 / 0.17 against 0.19 / 0.18 / 0.06 / 0.28 / 0.27) |
| E3 | old median mountain share ≥ 20 % | refuted (15.6 %) |
| E4 | collision never the only formation, any age | held for F180 (F176 old: 1 seed) |
| E5 | the équerres from the rift strip, < 5 % of the active length | refuted (the continent-ocean re-typing; the right-angle lines ~5 %, all offshore active cells 29.5 %) |

**Mine**:
- held: E2, E3 refuted (5–15 % predicted, 15.6 % measured: just above), E4, E5 refuted on the origin (the re-typing), « the plates
  to the author ».
- refuted: E1 « 6–10 seeds at old » (15), the young width « 15–30 km » (42).

## For the author

- **The plates to look at before the adoption is made effective**:
  - `images/formations_15_F180_age2.png` and `images/formations_15_F180_age3.png`;
  - `images/relief_active_15_F180_age2.png` and `images/relief_active_15_F180_age3.png`;
  - `images/transect_F180_old_seed7.png` and `images/map2048_F180_old_seed7.png`;
  - against F180's young maps.
- The équerres: two blocks re-typed by the continent-ocean reassignment. A queued item, not corrected.

**At commit (2026-10-18), decided**: the author's words: « Feu vert pour commiter le F181 (le rapport et le banc). La Partie A n'est pas encore appliquée : elle attend la décision de l'auteur sur l'âge par défaut. »
- F181 is committed as the report and the bench. **Part A is not applied**: it waits for the author's decision on the default
  age.
- **The reviewer's errors, recorded (the author's)**:
  1. **F179**: « the relief accepted since F166 does not rest on pumping » was concluded from a 4–7 % fall in s, without the
     conversion by Airy. But 0.07 in s is about 400 m. At the old age, the bound takes the median from 202 to 53 m and the
     mountain from 40 to 16 %.
  2. **E3**: the prediction « at least 20 % of mountain at the old age » rested on the same error.
- **Queued**: the continent-ocean reassignment (`subduction/source_term.rs:266–268`) creates a continental type with almost no
  crust.
  - **Production passes through this code**: `C1RunSpec::default()` takes `C1Closures::default()`
    (`crates/ymir-viz/src/bridge/c1/spec.rs:58`), whose subduction is enabled (`subduction/params.rs:70`, `enabled: true`), so
    `apply_subduction_step` runs in `run_hd`'s C1.
  - Nothing corrected.
- **Next, while the author decides**: the attribution of the fall between F176 and F180, the bound against the seed
  (« F176 + the bound alone », 256², intermediate and old ages).
