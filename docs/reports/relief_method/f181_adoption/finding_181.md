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

## F181 — attribution (while the author decides on the default age; nothing else changed)

**« F176 + the bound alone »** (no Andean seed): `Margin { ds_mass_bound: Some(1.0), ..seaward_oo_collision() }`, 256², the 15
seeds, the intermediate and old ages. The budget was announced at ~15 min; it took 1 612 s (27 min), over the estimate and under
45 min.

| | F176 | F176 + bound | F180 (+ bound + seed) |
|---|---|---|---|
| intermediate: cordilleras / collision alone | 9 / 0 | **4 / 5** | 15 / 0 |
| intermediate: mountain % median (range) | 47.4 (7.5–79.6) | **12.6 (1.3–39.7)** | 23.1 (10.8–47.7) |
| intermediate: median altitude | 303 m | **51 m** | 115 m |
| old: cordilleras / collision alone | 9 / 1 | **4 / 4** | 15 / 0 |
| old: mountain % median (range) | 40.0 (0.3–75.8) | **7.0 (0.0–30.7)** | 15.6 (5.8–40.5) |
| old: median altitude | 202 m | **32 m** | 53 m |
| old: land median | 84 541 km² | 76 321 km² | 82 288 km² |

**The split of the change F176 → F180** (per seed, medians; the bound's part = (F176 + bound) − F176, the seed's part =
F180 − (F176 + bound)):

| | F180 − F176 | the bound's part | the seed's part |
|---|---|---|---|
| intermediate, mountain (points) | −15.2 | **−28.2** | **+13.8** |
| intermediate, median altitude | −119 m | −234 m | +64 m |
| intermediate, W₁₀ global | −0.016 | −0.131 | +0.103 |
| old, mountain (points) | −15.5 | **−27.6** | **+12.7** |
| old, median altitude | −76 m | −154 m | +49 m |
| old, W₁₀ global | −0.031 | −0.121 | +0.093 |

- **The whole fall comes from the bound, and more**: alone it cuts ~1.7–1.8 × the total change of the mountain share (the
  summed shares 172–179 %) and ~1.4 × the median altitude. **The seed gives back about half of it**: +13 points of mountain,
  +50–65 m of median altitude, and most of W₁₀ (its part +0.09–0.10 against the bound's −0.12–0.13).
- **The bound alone would fail condition (v)**: cordilleras on 4 / 15, the collisional highland the only formation on 4–5
  seeds, the mountain share 7–13 %, the median altitude 32–51 m. The bound removes the pumped relief everywhere, the
  collisions included (the reviewer's error recorded at commit: 0.07 in s is about 400 m). The Andean seed puts relief back
  where the kinematics supplies it, on the active coasts.
- So F176's relief was largely pumped (the bound takes the old median from 202 to 32 m), and F180's relief is the bounded
  relief plus the Andean cordilleras.

**The author's decision (verbatim)**: « Je penche pour B aussi. On pourra inspecter plus tard pourquoi C. »
- **B: the adoption is effective, and the continent profile's default age becomes intermediate (A = 100 Myr).** The age slider
  stays available.
- **From « on part sur du vieux » (F168) to intermediate.** Under Davis-Suppe's mass bound, the old profile gives 16 % of
  mountain and a median altitude of 53 m, against 24.8 % and 195 m for the Europe reference; the intermediate profile gives
  23 % and 115 m. (« On part sur du vieux » had been a design intent; no code carried it as a default.)
- **Applied**:
  - `Margin::continent_default()` = `seaward_oo_collision_f180()`;
  - `select_continent_default` and `select_continent_plates` take it;
  - `cascade::profile::CONTINENT_DEFAULT_AGE_YEARS = 100 Myr`;
  - in the viz, entering a continent profile (L1, L2 or the selection) checks the adopted boxes and sets the age slider to
    100 Myr, and going back to production clears the boxes;
  - the old benches (F170–F174) call the `_margin` entries with `Margin::off()` explicitly, so their sense is unchanged;
  - the ignored F174 control now checks the default selection against `Margin::continent_default()`.
- **No bench or test assumed a default age** (they pass their ages explicitly); `PhysOpts`'s 5.62 Myr is the P0 level's
  (F159–F162) and stays.
- **No `ALGO_*` bump**: the continent profile runs only in the cascade, which no cache key covers. The checks below verify
  production's hash and the guard.
- **Queued, in this order**:
  1. F182, the deposition (its reference profile becomes the new default: F180 at the intermediate age);
  2. the trenches' coherence (polarity per segment) and the 2-D wedge (Dahlen p. 60–61), together;
  3. the rift strip that becomes a plate;
  4. the continent-ocean reassignment that creates a continental type without crust (`subduction/source_term.rs:266–268`;
     production passes through it, checked at commit);
  5. C: why are the old continents lower than the Europe reference? A sourced reference on the hypsometry of old continents is
     needed (the attribution « F176 + the bound alone » above serves it);
  6. the usual queue.
