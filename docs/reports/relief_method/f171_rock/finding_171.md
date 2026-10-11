# Finding 171 — a rock that depends on depth: a soft cover on a hard basement (K1), then the heterogeneous K (K2)

**Status: built, measured, gated. Production unchanged. Nothing committed without the feu vert.**

**Verdict**:
- **The fallback works**: the témoin qualifies at 20 plates (r* 2); 42, 1 and 9 stay at 16.
- **K1 removes the peneplain on 4 / 4 seeds.**
  - The old profile's median land rises from 3–95 m (K0) to 132–423 m; Europe is 195 m.
  - 88–91 % of the old land is basement.
  - The plains standing on it lie at 123–260 m, near the reviewer's ~180 m estimate.
- The mountain share now decreases with A gradually instead of collapsing at 20 Myr.
- **The macro success is 1 / 4** (seed 42, K1 and K2). It needed 3.
  - The old mountains spread from 6.5 to 37 % across the seeds.
  - Seeds 1 and 9 fail the coast clause: their K0 reference was the drowned peneplain's intricate coast.
- **No stop fires.**
- K2 adds little: it lowers R2_g, does not raise σ_θ by 2°, and moves the classes by up to 10 points.

**Files**:
- Declared before any measurement: `f171_declared.md`. Predictions: `f171_predictions.md` (written first).
- The bench `crates/ymir-core/tests/f171_cascade.rs`:
  - `f171_rock` → `bench_output.txt`;
  - `f171_river` → `bench_river.txt`, the river-profile image redone.
- `images/`.
- Checks: `checks_before.txt`, `checks_after.txt`.

## Part 0

- F170 committed (`f21c5c2`). Its ADR carries the torus default and its cost, the 16-plate selection, the reviewer's error on
  the sutures, the peneplain with the reviewer's reading, and the coast and confluence notes.
- The author's two answers (« oui et oui »: K first and the sea-level rise at F172; 16 plates by default with the 20 / 24
  fallback and a slider) are added under it.

## A — the plates by default

| seed | N 16 | N 20 | retained | the chosen continent | the selection's cost |
|---|---|---|---|---|---|
| témoin | none (best r* 0) | compact k 9, start 0, 51.5 %, **r* 2** | **20** | 51.5 %, r* 2 | 145 s CPU |
| 42 | compact k 7, r* 2 | — | 16 | 51.3 % | 45 s |
| 1 | BFS k 10, r* 1 | — | 16 | 50.4 % | 53 s |
| 9 | compact k 7, r* 3 | — | 16 | 51.2 % | 52 s |

- **The code**: `profile::select_continent_default` (16 → 20 → 24) and `k_range`.
- **The viz**:
  - a slider « Plaques » (8–24) drives « sélection ». At 16 it falls back; at any other value it selects at that count alone
    and shows the result in the progress bar.
  - A « Roche : K0 / K1 / K2 » selector feeds the rock law to the history physics level and K2 to the chain.
- **The negative control**: 8 plates without the selection is F168 bit for bit (the existing test).

## B, C — K1 and K2 (`cascade::Rock`, `physics_history_rock`, `amplify::KHetero`)

- **K1**: K = K_couv while E_cum < d_c, K_couv / ρ_K beyond.
  - d_c = 0 on C1's shields.
  - K is evaluated at the start of each step; the crossing step uses the fraction-of-step rule.
  - ρ_K = 10 and d_c = 1.5 km, declared: Stock & Montgomery 1999 p. 4983 and Table 2 p. 4986; d_c is a design value.
- **K2**: K1 × exp(ln M).
  - ln M is an fBm of 7 octaves, 25 km down to 0.39 km, σ_lnK = 0.5 (Stock & Montgomery's within-lithology spreads).
  - Each octave is born at its level and carried by the upsample and the warp, so it is fixed to the rock.
  - It applies to the physics level and to every amplification level's k.
- **The negative controls** (permanent; `the_rock_law_reduces_to_f170_and_k1_bit_for_bit`):
  - ρ_K = 1 is F170 bit for bit; σ_lnK = 0 is K1 bit for bit, at the physics level and at an amplification level;
  - ρ_K = 10 changes the field and exhausts the cover.

## D — the measures

Europe: plain 51.6 / plateau 4.7 / hill 18.9 / mountain 24.8 %; median 195 m.

### The mountain share against A

| seed | K0 (B0 → 20 → 100 → 300 Myr) | **K1** | K2 |
|---|---|---|---|
| témoin | 80.2 → 34.5 → 21.0 → 18.4 % | 81.3 → 75.6 → 62.6 → **37.0** % | 81.8 → 75.3 → 62.9 → 40.5 % |
| 42 | 72.8 → 11.6 → 7.1 → 6.2 % (F170) | 68.0 → 60.2 → 38.2 → **20.3** % | 67.6 → 59.0 → 44.5 → 21.0 % |
| 1 | 52.2 → 0.6 → 0.3 → 0.3 % (F170) | 59.7 → 58.6 → 31.7 → **6.5** % | 62.4 → 56.8 → 32.4 → 9.0 % |
| 9 | 59.9 → 8.1 → 2.6 → 2.4 % (F170) | 76.2 → 67.9 → 38.7 → **13.9** % | 76.3 → 66.7 → 44.1 → 21.3 % |

### The old profile (300 Myr)

| seed | variant | plain / plateau / hill / mountain | median | peak (2 048²) | cover / basement / shield | plains on the basement |
|---|---|---|---|---|---|---|
| témoin | K0 | 55.0 / 0.2 / 26.4 / 18.4 | 95 m | 2.0 km | 90 / 0 / 10 | — |
| | K1 | 16.4 / 9.0 / 37.6 / 37.0 | 423 m | 5.2 km | 1 / 88 / 10 | 260 m |
| | K2 | 17.2 / 5.4 / 36.8 / 40.5 | 348 m | 5.6 km | 1 / 88 / 10 | 151 m |
| 42 | K0 | 76.0 / 0.0 / 17.8 / 6.2 | 50 m | 2.7 km | 91 / 0 / 9 | — |
| | K1 | 55.9 / 0.5 / 23.3 / 20.3 | 172 m | 3.9 km | 2.5 / 88.5 / 9 | 126 m |
| | K2 | 49.8 / 0.8 / 28.4 / 21.0 | 155 m | 4.2 km | 2.5 / 88.5 / 9 | 97 m |
| 1 | K0 | 96.6 / 0.0 / 3.1 / 0.3 | 3 m | 0.8 km | 92 / 0 / 8 | — |
| | K1 | 49.9 / 0.2 / 43.4 / 6.5 | 132 m | 1.8 km | 2.5 / 90 / 7.5 | 123 m |
| | K2 | 48.5 / 0.1 / 42.4 / 9.0 | 108 m | 2.0 km | 1.3 / 91 / 7.5 | 97 m |
| 9 | K0 | 91.4 / 0.0 / 6.2 / 2.4 | 9 m | 1.5 km | 91 / 0 / 9 | — |
| | K1 | 44.8 / 2.9 / 38.3 / 13.9 | 190 m | 4.3 km | 1.6 / 89 / 9.3 | 126 m |
| | K2 | 34.9 / 2.7 / 41.1 / 21.3 | 179 m | 4.3 km | 1.6 / 89 / 9.3 | 75 m |

### The rock at the surface (K1)

- B0: cover 34–47 %, basement 45–56 %.
- Young: cover 10–18 %, basement 73–81 %.
- Intermediate: cover 2–4 %, basement 87–89 %.
- Old: cover 1–2.5 %, basement 88–91 %. The shields stay 7.5–10.5 %.
- The 1.5 km cover is gone from half the land already at B0 (5.62 Myr), and from three quarters by 20 Myr.
- Map: `images/f171_rock.png`.

### Peaks and the slope guard

- **K1 raises the peaks.**
  - B0: 5.2–6.4 km at 2 048², against F170's 2.0–3.7 km and the témoin's K0 4.1 km.
  - Young: 4.8–5.7 km, against 1.0–2.9 km.
  - Never above 7 km.
- **The slope guard holds** on every K1 / K2 world: the mountain-class slope p50 is 0.25–0.45 against Corsica's 0.318, within
  ×1.5. The one failure of the whole bench is K0 old, seed 1 (0.128).
- **But the facets and the walls fail in the mountain class** at B0 and young, on 4 / 4 seeds, for K1 and K2. With K0 only the
  facets failed there.

### The lithological knickpoints (the declared detector, 256², Strahler ≥ 3)

| | B0 | young | intermediate | old |
|---|---|---|---|---|
| count (4 seeds) | 31–55 | 57–78 | 5–13 | 0–8 |
| median height | 184–315 m | 30–82 m | 20–86 m | 3–19 m |

- They are most numerous while the contact is still crossing the land (B0, young). They fade once the cover is gone.
- **The river profile** (`f171_river_contact.png`; the témoin, K1 intermediate):
  - an upstream basement reach, a flat crossing a closed depression (a lake), a short cover reach;
  - then a 1.1 km drop on the basement to the coastal cover (blue tick: the knickpoint).

### Regularity (391 m), K1 → K2 at the same age (16 pairs)

- **σ_θ** changes by −0.1 to +1.2°, never +2°. It stays 23–26° against Corsica's 29.9°.
- **R2_g falls on 15 / 16 pairs**, by up to 0.155 (seed 1, intermediate: 0.797 → 0.642). Corsica is 0.48.
- **CV_L**: 0.96–2.17. **Hovius R** (A ≥ 10 km²): 2.3–8.7, and 65 on one B0 world with few basins.
- Crops: `f171_crops_25km.png` (Corsica | K1 | K2 at B0 and old, the témoin).

### The coast on the torus (old)

| seed | K0 | K1 | K2 |
|---|---|---|---|
| témoin | 1.00 / 0.97 / 0.31 / 1.20 | 1.01 / 0.95 / 0.30 / 0.77 | 1.01 / 0.95 / 0.30 / 0.82 |
| 42 | 0.99 / 0.95 / 0.31 / 0.96 | 1.00 / 0.96 / 0.30 / 0.86 | 1.00 / 0.96 / 0.29 / 0.91 |
| 1 | 1.06 / 1.10 / 0.34 / 1.18 | 1.01 / 0.97 / 0.31 / 1.12 | 1.00 / 0.96 / 0.31 / 1.18 |
| 9 | 1.07 / 1.18 / 0.34 / 1.04 | 1.01 / 1.00 / 0.30 / 0.83 | 1.01 / 1.01 / 0.30 / 0.82 |

- The columns are dimension / length ratio / aligned / block index. Corsica: 1.12 / 1.37 / 0.27 / 0.78.
- Seeds 1 and 9's K0 coast was the drowned peneplain's, more intricate and closer to Corsica's length ratio. K1 brings back
  ≈ 1.0, which the clause counts as « degraded > 10 % ».

### Sensitivity (the témoin, old profile)

| run | mountain | plain + plateau | median | peak |
|---|---|---|---|---|
| ρ_K 3 | 26.8 % | 52.1 % | 183 m | 3.3 km |
| **ρ_K 10 (central)** | 37.0 % | 25.4 % | 423 m | 5.2 km |
| ρ_K 30 | 58.1 % | 22.5 % | 649 m | 5.9 km |
| d_c 0.75 km | 38.1 % | 28.7 % | 432 m | 5.4 km |
| d_c 3 km | 31.2 % | 37.4 % | 278 m | 4.9 km |
| K2 σ 0.25 | 38.5 % | 23.1 % | 395 m | 5.4 km |
| K2 σ 0.5 (central) | 40.5 % | 22.6 % | 348 m | 5.6 km |
| K2 σ 1.0 | 42.5 % | 20.8 % | 307 m | 5.9 km |

- **ρ_K dominates.** At ρ_K = 3 the témoin's old profile meets every macro clause (26.8 / 52.1 % / 183 m).
- That is a measured sensitivity, **not a choice**: ρ_K stays at its declared 10.

### Cost (idle)

- 434–490 s CPU per world, K0, K1 and K2 alike: the rock law and the carried ln M cost nothing visible.
- The selection: 45–53 s CPU per seed at 16 plates, 145 s for the témoin (16 then 20).

## R

| rule | result |
|---|---|
| the témoin's fallback | passes at 20 plates; the « ≥ 3 » rules stay « ≥ 3 of 4 » |
| macro success, old, K1 or K2 | **1 / 4** (seed 42); seed 9 K2 fails the coast clause only; the témoin fails plain + plateau (25 %) and its K1 median (423 m); seed 1 fails the mountains (6.5 %) and the coast |
| a persistent peneplain (median < 50 m, old, K1) | 0 / 4 |
| peaks off the scale (> 7 km or the slope guard failing, K1 B0 / young) | 0 / 4 (peaks 4.8–6.4 km) |

No stop fires as written; the macro success fails.

## The reading (hypotheses to check)

1. **A depth-dependent rock gives an old continent a relief** (the reviewer's reading holds on the median: 3 / 4 within
   100–300 m). It keeps a basement plain at 100–260 m with residual ranges.
2. **The cover is thin against the model's erosion**: half of it is gone at B0. The hard basement then carries the relief at
   every age, which explains the high B0 / young peaks (5–6 km) and the new walls.
3. **The spread between the seeds is the tectonics'**: the témoin (20 plates) keeps 37 % of mountains, seed 1 only 6.5 %, at
   the same rock.
4. **The coast clause measures against a flooded reference**: K0 old's intricate coast came from its sea-level plains. The
   sea-level rise (F172) is the natural test of the coast.
5. **K2, as declared, does not reach Corsica's confluence spread**: σ_θ is still 23–26°, against 29.9°. It does decorrelate the
   gullies (R2_g down).

## Predictions

### The reviewer's

| held | refuted |
|---|---|
| X1 (the témoin at 20) | X4's second clause (more knickpoints on the young than on the intermediate) |
| X2 (median 100–300 m on 3 / 4) | X6 (σ_θ never +2°; R2_g does fall, 15 / 16) |
| X3 (old mountains in band on 2; peaks > 1.5 km on 4 / 4) | X7 (up to 9.9 points: seed 9 old plain 44.8 → 34.9 %) |
| X4's first clause (≥ 3 seeds) | |
| X5 (+46 to +380 % against F170 at B0 and young, 3 / 3; the témoin +54 % against K0) | |
| the meta | |

### Mine

| held | refuted |
|---|---|
| X1 at 20 | X2 refuted (it held) |
| the peneplain stop not firing | X3's first clause refuted (it held) |
| X3's second clause | X4's « more on the intermediate » |
| X4 on ≥ 3 | X5 refuted (it held) |
| no macro success | X6 (σ_θ +1–4°: no; R2_g unchanged: no) |
| ρ_K 30 ×1.5 the median (×1.53) | X7 (it was refuted) |
| σ ×2 < 5 points | the rock shares (the basement 88–91 %, not 50–85 %; the cover 1–2.5 %, not 5–40 %) |
| the costs | d_c < 30 % (3 km: −34 %) |
| the meta | |

## Open, for the author (nothing chosen)

- **(a)** K1 / K2 as the continent profile's rock (the median and the plains now Europe-like on 3 / 4).
- **(b)** The spread of the old mountains (6.5–37 %) and ρ_K's sensitivity: the témoin passes at ρ_K = 3, not at 10. A choice
  of ρ_K on our worlds would be a calibration; the literature range is 10–10⁵.
- **(c)** The coast clause's reference, before F172 (the sea-level rise) changes the coasts anyway.
- **(d)** The walls and facets that the hard basement brings at B0 and young.

## At commit (2026-10-11)
- **(a)** The continent profile's rock becomes **K2** (K1 + the heterogeneity, σ_lnK = 0.5), with 16 plates and the fallback to
  20 then 24.
- **(b) A change of criterion, recorded as such**: the seeds are no longer judged one by one against Europe.
  - A mechanism is judged on its own effect, which must appear on nearly every seed.
  - The macro layout is judged on a **distribution** of seeds: the share of interesting seeds, the variety.
  - The age A is the author's slider. Its default is **provisional**, to be set when the work items are finished.
  - ρ_K stays at 10.
- **(c)** From F172, the coast is judged against Corsica in absolute value.
- **(d)** The walls, the facets and the coastal rampart are measured at F172 (its parts A and D).

**The author's view of F171, verbatim**: « Assez intéressant. Donc maintenant j'ai un moyen de régler l'age. Le K1/K2 sur le
crop 25km est très intéressant. Les montagnes semble réalistes et irégulières. Sur la 2eme ligne un peu moins. Donc à B0 ça
va, en old un peu moins. »

**The author's position on the seeds, verbatim**: « Age par défaut? Je ne sais pas. Il bougera certainement lorsque tous les
chantiers seront terminés. Ensuite, toutes les seeds ne seront pas intéressantes. C'est ce qu'il faut voir et on le sait.
Car toutes les seeds ne vont pas donner des structures convergentes aux bordures systématiquement. C'est bien évidemment le
"jeu" ici. Il faut trouver des seeds intéressantes. Elles allient à la fois forme des plaques, nombre de plaques de terres,
direction et vitesse de convergence / divergence / rift / etc... C'est multifactoriel et donc c'est pour ça que les seeds 1
à 9 peuvent ne pas forcément donner un résultat intéressant (typiquement continent avec mur, fini avec quasi que des
plaines, etc). »
