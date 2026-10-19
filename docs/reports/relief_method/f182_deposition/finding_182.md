# Finding 182 — the sediment deposition at the physics level: built behind the flag; the declared stop fires before the tiers, because on the closed torus the deposits fill the ocean

**Status: Part A is built and its controls pass. Part B stopped at the cost probe and the panel diagnostic, by the brief's stop
« la circumnavigation perdue sur ≥ 2 graines »: with the deposition the sea vanishes entirely on 9 of the panel's 10 worlds
(young and intermediate), in 14–36 Myr.** The tiers, the reference (Yuan's iterative coupling) and the sensitivity have not run.
Nothing is committed. The reference profile is F181's default (`Margin::continent_default()`, the intermediate age 100 Myr).

Declared before any measurement: `f182_declared.md` (Part 0 and Part A). Predictions: `f182_predictions.md`.

## Part 0 (earlier, unchanged)

- **0.1, today**:
  - the physics level is detachment-limited; the eroded matter disappears (`stream_power.rs:716 ff.`);
  - the lakes exist only in the routing (priority flood, `stream_power.rs:645`);
  - the sea cells are fixed base nodes (`:664`), reset after the hillslope step (`cascade.rs:670`).
- **0.2, the sources**: Yuan et al. 2019 JGR (the fluvial deposition, eq. 8–15; a Gauss-Seidel iteration, eq. 19–21) and EPSL
  (the marine diffusion, eq. 5–6, Table 1).
- **0.3, the scale**: F162-0's convention.
- **0.4, the stability**: the explicit marine diffusion needs sub-steps of ≤ 1.53 kyr at K_M = 400 m²/yr.

## The resumption (the author's complement)

- **The reviewer's error, recorded (the author's)**: the reviewer took « Harris ÷ 7.5 » as the reference for a shelf built by
  the deposition, without rereading F162-0's convention.
- **The principle, confirmed by the author**: the tectonic structures read at the signified scale (÷ 7.5); the surface
  processes (erosion, diffusion, deposition) run at the real scale of the 400 km map, with no conversion. (The author's
  verbatim answer was not included in the complement.)
- **G = 1.55**: Guerit et al. 2019, *Geology* 47(9), p. 855, natural systems « ~0 to 3.1 », the middle as the value, the bounds
  as the sensitivity.
- **The route and the marine scheme** are the author's decisions (the explicit splitting; one silt class, explicit sub-steps).

## Part A — what is built (`crates/ymir-core/src/cascade/deposition.rs`, behind `Drive::dep`)

- **The fluvial pass**: the incision is unchanged. Then one pass runs upstream → downstream on the incised surface's routing,
  depositing D_i = min((G / Ñ_i) Q_in, Q_in) (Yuan JGR eq. 13–15 with the step's quantities; the author's cap).
- **The lakes**: a lake takes no fluvial deposit. At its outlet the sub-basin fills from the bottom up by volume, at most to
  the spill; the remainder continues.
- **The marine step**: one class, the silt, K_M = 400 m²/yr and φ = 0.63 (Yuan EPSL), with the sand omitted (declared).
  - Explicit sub-steps of ≤ Δx² / (4 K_M); no iteration, no linear system.
  - Fluxes run between sea cells only, and the transport moves only the deposited sediment (the one declared departure from
    eq. 5).
- **The hillslope's loss to the sea** is now a marine source, not deleted.
- **The deposits' load is not flexed**; the rebound reads the incision alone, as before.
- **The mass control**: e + H = land + lakes + stranded + sea solid + exported, per step in f64. The sea solid is measured on
  the field.
- **The reference** (Yuan's iterative coupling, eq. 19–21) is in `tests/f182_cascade.rs` only, passed through
  `Drive::dep_hook`. The cascade's code holds no iteration.
- **Entry**: `physics_history_dep` (`physics_history_sea` = it with `None`).
- **Viz**: the box « Dépôt (niveau de physique) » (history mode) and the layer « Épaisseur de sédiments ».

### The controls (A5, permanent: `cascade::history::f169_tests::the_deposition_controls`, and the module's two unit tests)

- `Deposition::nothing()` (G = 0, no lake fill, nothing deposited at sea) is the level without deposition **bit for bit** (the
  field and the rebound).
- G = 0 deposits nothing fluvial.
- The declared deposition changes the field; the mass control holds < 10⁻⁹; nothing is exported on the torus.
- The pass conserves at 10⁻¹² for G = 0, 1.55, 3.1, with and without lakes.
- The marine step conserves the solid volume and never lowers the bathymetry below its start.

## Part B — the cost, then the stop

### The cost per world (seed 8, 256², `f182_cost`, the marine step then parallel)

| age | without | with | marine sub-steps | marine wall | the mass control's worst step |
|---|---|---|---|---|---|
| young 20 Myr | 16.3 s wall / 22.9 s CPU | 42.0 s / 171.9 s | 11 712 | 12.2 s | 4.1 × 10⁻¹³ |
| intermediate 100 Myr | 16.3 s / 23.2 s | 48.7 s / 201.5 s | 22 261 | 16.8 s | 5.8 × 10⁻¹⁴ |
| old 300 Myr | 16.0 s / 22.7 s | 79.2 s / 504.8 s | 57 448 | 46.9 s | 9.3 × 10⁻¹⁴ |

- The marine step has since been made serial (rayon's overhead dominated a ~10⁵-flop sub-step: CPU ≈ 8 × wall). The per-cell
  arithmetic is unchanged; the serial run counts the same 22 261 sub-steps on seed 8.
- **The sub-step counts were the clue.** The intermediate age should need ~65 500 sub-steps, but the probe counted 22 261,
  because the marine step stops running when there is no sea left. On seed 8 at the intermediate age: 40 774 sea cells at
  step 0, 13 618 at step 24, 6 at step 27, 0 from step 28 (33.9 Myr) to the end.

### The stop: the deposits fill the ocean (`f182_sea_fill`, the panel, young and intermediate)

| seed | age | eroded without the deposition (km³) | ocean at t₀ (km³) | eroded / (1 − φ) / ocean₀ | with it: sea under half | sea gone | land at the end, with / without |
|---|---|---|---|---|---|---|---|
| 8 | young | 1.40 × 10⁵ | 4.29 × 10⁵ | 0.88 | 12.1 Myr | 17.7 Myr | 100 % / 51 % |
| 18 | young | 1.06 × 10⁵ | 4.01 × 10⁵ | 0.71 | 16.9 Myr | never | 87 % / 49 % |
| 19 | young | 1.22 × 10⁵ | 3.36 × 10⁵ | 0.98 | 8.7 Myr | 13.8 Myr | 100 % / 53 % |
| 23 | young | 1.48 × 10⁵ | 2.98 × 10⁵ | 1.35 | 8.3 Myr | 13.6 Myr | 100 % / 62 % |
| 29 | young | 1.45 × 10⁵ | 3.39 × 10⁵ | 1.15 | 11.0 Myr | 14.8 Myr | 100 % / 56 % |
| 8 | intermediate | 1.93 × 10⁵ | 4.29 × 10⁵ | 1.21 | 28.9 Myr | 33.9 Myr | 100 % / 51 % |
| 18 | intermediate | 1.31 × 10⁵ | 4.01 × 10⁵ | 0.88 | 30.9 Myr | 35.9 Myr | 100 % / 49 % |
| 19 | intermediate | 1.59 × 10⁵ | 3.36 × 10⁵ | 1.28 | 19.9 Myr | 25.6 Myr | 100 % / 53 % |
| 23 | intermediate | 1.84 × 10⁵ | 2.98 × 10⁵ | 1.67 | 18.7 Myr | 24.5 Myr | 100 % / 62 % |
| 29 | intermediate | 1.87 × 10⁵ | 3.39 × 10⁵ | 1.49 | 22.2 Myr | 26.7 Myr | 100 % / 56 % |

- **The criterion is independent of the scheme.** Over the history, the incision alone erodes 0.71–1.67 times the ocean's volume
  at t₀ once inflated by the surface porosity (eroded / (1 − φ)). On a closed torus that keeps everything (the brief's A4: « le
  tore ne perd rien »), with no load, no compaction and no recycling, the ocean must fill.
- **The filling is faster than the ratio says.** Even at ratio 0.71 (seed 18, young) the sea falls to 13 %. Two things drive
  it:
  - the sediment builds the shallow margins first (the transport moves only the deposits, downslope);
  - every emerged cell erodes in turn, so the eroded volume grows as the sea retreats.
- **The mass control holds throughout** (≤ 4 × 10⁻¹³). The stop is the physics of the closed budget, not a leak.
- The circumnavigation is lost on 9 / 10 worlds, by the brief's rule a stop. (v), W₁₀, the lakes, the shelf and the reference were
  not measured: the tiers did not run.
- **The orders of magnitude**:
  - the incision's 1.3–1.9 × 10⁵ km³ over 100 Myr on ~80 000 km² of land is 1.6–2.4 km of mean denudation, plausible under
    C1's tectonic sources;
  - the ocean's mean depth is ~5 km over half the map.
  - The Earth is not closed in this way: its sediment is recycled by subduction, it loads the lithosphere, and it compacts.

## The predictions

- P4 (the mass control holds): held. M1: held (≤ 4 × 10⁻¹³).
- P5 and M7 (the protection): refuted. The circumnavigation is lost on 9 / 10 panel worlds.
- M9 (the cost: the marine step ~10–40 s at 100 Myr): held at 16.8 s, but only because the sea was gone after step 28.
- P1–P3, P6, P7, M2–M6, M8: not measured.

## For the author

**What a closed budget lacks against the Earth, each sourced or not:**

1. **The sediments' isostatic load.** A column of thickness h in water subsides by h (ρ_s − ρ_w) / (ρ_m − ρ_w). The cascade's
   rebound reads only the unloading (max(0, e)). The densities need a source; the order is ~40–60 % of the rise taken back.
2. **The compaction**: Yuan EPSL eq. 7, Table 1 (the silt's e-folding depth 1 960 m, Sclater & Christie 1980), in the
   brief's source. The surface porosity 0.63 throughout overstates a km-thick column's thickness by up to ~1.7×.
3. **The recycling by subduction.** C1 renews the ocean floor, but the cascade's sea cells only follow C1's U; nothing carries
   the deposits to the trenches. This contradicts « le tore ne perd rien »: a decision, not a parameter.
4. **The erosion itself.** The cascade's incision runs on C1's tectonic sources, which C1's own erosion term already balanced
   once. This is out of F182's scope; it is noted because it sets the eroded volume.

**Measured: none of 1–3 alone obviously suffices.** At the intermediate age the ratio reaches 1.67 (seed 23). Even halving it
(1 + 2) leaves 0.8, and seed 18 young shows that 0.71 already drains the sea to 13 %. Recycling (3) is the only one that removes
mass, and it is the author's decision against the brief's A4.

**Not changed, waiting for the decision**: the code stays behind the flag (off by default, production unchanged); the tiers,
the reference and the sensitivity are ready in `tests/f182_cascade.rs` and resumable.

## The checks (after; the deposition off by default)

- `cargo check --workspace`: clean.
- 665 lib tests pass (F181's 662 plus the three deposition controls); 33 viz tests pass.
- The guard reads Match on all six states (field and lakes).
- Production's C1 hash (`tests/f175_c1_hash.rs`) is unchanged.

**At commit** -- the author's words (verbatim): « Feu vert pour commiter le F182 tel quel (le dépôt derrière son drapeau, coupé
par défaut ; contrôles A5 ; rapport avec l'arrêt). Tu peux supprimer target_f182/ si tu n'en as plus besoin, sinon garde-le
pour le F183. Ajoute à l'ADR une erreur du relecteur : la règle A4 « le tore ne perd rien » ignorait le recyclage des sédiments
par subduction, et le rapport entre volume érodé et volume de l'océan n'a pas été calculé à l'avance, alors qu'il annonçait
l'arrêt (0,71 à 1,67 sur le panel). »
- **The reviewer's error, recorded (the author's)**:
  - the rule A4 « le tore ne perd rien » ignored the recycling of the sediments by subduction;
  - the ratio of the eroded volume to the ocean's volume was not computed beforehand, although it announced the stop
    (0.71–1.67 on the panel).
- **Committed as is**:
  - the deposition behind its flag (`Drive::dep`, `None` by default; the viz box unchecked);
  - the A5 controls;
  - the report with the stop.
- Production is unchanged: the C1 hash, the guard 6/6 Match, 665 lib and 33 viz tests.
