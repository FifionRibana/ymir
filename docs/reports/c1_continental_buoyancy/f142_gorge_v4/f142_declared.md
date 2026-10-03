# F142 — instruments, declared BEFORE the measurements (2026-10-03)

## R — re-tabulated as built, the input alone (bench `f142_r`)

- **The bodies**: `gorge_bodies` of the GORGE skeleton (r 1, p 0.5) on S1.
  - The 14 D8 lakes are F139's set (ON ×1's final lakes ≥ 1 km², not craters, with a D8 receiver at the outlet).
  - Each lake's body is the one with its largest footprint overlap. Both shares are printed: the body's and the
    lake's.
- **Every body's class**, for the scope's criterion: input lake or depression component, whether its ring touches
  another body's depression cell or a cell under the sea, and its final lakes by overlap.
- **T3**: the as-built L_in, L_bed and L_floor (S1). L(r) at r = 0 / 0.5 / 1 / 1.5 / 2. The area under L(r) grows
  8-connected from the body's lowest cell on S1.
  - Compared with F138-T3's totals (1 308.5 / 1 223.6 / 1 145.7 / 367.3 / 0 km²).
  - Per lake: L_bed against today's level (deviation 2).
- **R**: r_lake = min(2, r_world·(A/418.7)^p) with A = the col's area as built. For p ∈ {0, 0.25, 0.5} and r_world ∈
  {0, 0.5, 1, 1.5, 2}: the lakes ≥ 1 km² left and the total area, against F139-R's table. The emptying age by map (a),
  per lake, against F139-R. The as-built median of A is printed against A_ref.
- **T4**: the construction's own head falls, slope and shortage (steepen, m = 10), at r_world 1 with p = 0
  (r_lake = 1, F139-C4's level) and p = 0.5. H_f = max(1.5·S·100 m, 10 m). Compared with F138-T4 / F139-C4's D8
  lines (`f139_t.txt`).

## G — the gates (bench `f142_g`), declared before it runs

- **The worlds**: GORGE v4 (`GorgeRetreat::v4`, k = F121_AGE_K × age as F140) at:
  - ×1, ×0.7, ×0.85, ×1.2, ×1.4 (p = 0.5);
  - ×1 (p = 0, 0.25);
  - ×1.4 (p = 0).
  - ON ×1 first (the references); ON and OFF ×1 last (the controls).
- **Each world runs inside `catch_unwind`.** A panic is reported as a failed gate with its message; the run_hd tail
  is isolated a second time, so a panic there is F38's invariant firing and is named so.
- **G-bodies**: ON ×1's 14 D8 lakes (F139's set). Each lake goes to the kept body with its largest overlap.
  - Gate: 14 matched, 14 distinct bodies, the body's share > 50 %, and no lake overlapping > 1 body.
  - The kept bodies that match no D8 lake are listed with their ON lake.
- **G-pits**: the cells a footprint-restricted fill (`gorge_plain` on a copy, r_lake > 1 bodies) would raise by
  > 0.1 m, at:
  - the bench's construction + plain;
  - the pipeline's (a) construction + plain + rims, (b) + the light pass and (c) + the droplets (built only when a
    body has r_lake > 1);
  - (d) the full world, which is the breach's input.
  - At (a) the pipeline's field is compared bit for bit with the bench's construction + plain on the kept bodies
    (the production wiring of the plain).
  - Gate: 0 at every stage.
- **G-sea**: the land cells (world > sea) ≤ sea after the protected breach on the pre-breach drainage (as run_hd).
  Compared cell for cell with ON ×1's set (added, missing). Gate: identical.
- **F38**: the run_hd tail completes.
- **G-deposit**: `gorge_plain`'s records per body (cells raised > 1 cm, km², deepest, km³).
  - Against the gorge's removal: Σ max(0, ON − GORGE) over land on the constructions at the same age, km³, before
    the plain.
  - Reported, no gate (the reviewer's < 5 % at ×1.4 is a prediction).
- **F140's gates**, as in `f140_g`, on the kept bodies, with the construction + plain as "the construction":
  - G-levels (drained bodies: no final lake over > 50 % of the body);
  - G-rim, G-ring, G-slope100 (negative control: ON's 14 D8 outlets), G-drained, G-tag.
  - G-area: the 14 D8 bodies' final lakes by age at p = 0.5, against F142-R's re-tabulated R (non-increasing is the
    gate).
- **Extras** (canyons, coast, θ with the construction mask, removal and raise in km³ inside / outside the kept
  bodies): at GORGE ×1 p 0.5, GORGE ×1.4 p 0, ON ×1 and OFF ×1 (rule 14, rule 18).
- **Cost**: construction + plain against ON's construction at the same age, per world.

## Amendment to G-pits, after the first run, declared NON-BLIND (2026-10-03)

**The first run** (`f142_g_first_run_defective.txt`, stopped in its first GORGE world) measured G-pits as "the cells
a footprint fill would raise above the ring's minimum".
- At ×1, p 0.5 that read 0 at the construction and at (a), and **159 239 cells, up to 391 m, from (b) on**.
- Body 1 (lake 2) alone had 21 555 of its ~24 300 cells counted, 391 m = its lake's level minus its floor.
- **The measure took the lake to be the cells at the ring's minimum.** Once the light pass lowers the col by a few
  centimetres, the lake is held by a sill INSIDE the footprint, and the whole lake reads as "pits".

**Amended** (the bench's `pits` closure):
- per body with r_lake > 1, the same footprint-restricted fill seeded on the ring;
- **the lake** = the depression cells (fill > z + 0.1 m) at the fill level of the body's lowest cell (none if that
  cell is not in a depression);
- **a pit** = a depression cell at any other level (> 1 cm away);
- per body it prints L(r), the ring's minimum, the lake's level, the pit cells, their distinct levels and the deepest.

**What it no longer catches**: a drained body (r = 2) whose lowest cell holds water reads that water as "the lake".
G-levels catches it ("drained, absent at the end").

The rest of the bench is unchanged. At (a), the bench-against-pipeline difference now prints its bodies and its
largest value.

## Amendment, θ's attribution (bench `f142_theta`), declared after f142_g's x1 world and before it runs

f142_g read theta = 0.264 [0.253, 0.274] at GORGE v4 x1 p 0.5, against OFF's CI 0.493 [0.492, 0.495].
On the construction alone (ON x1 and GORGE v4 x1 + the plain, OFF's trunk links, the construction's mask), theta is
recomputed:
- each world with its own mask;
- GORGE without the links touching a cell whose skeleton base the gorge RAISES (> 0.5 m, the invariant on the
  outlet paths);
- GORGE without the links touching a cell whose base it LOWERS (> 0.5 m, the rebase upstream of the lowered lakes);
- GORGE without both;
- ON on those same links;
- GORGE on each class alone.

## Amendment, G-sea's attribution (bench `f142_sea`), declared after f142_g's x1.4 p 0.5 world (G-sea: 404 against 292)

At GORGE v4 x1.4 with p 0.5 and p 0, the land cells <= sea after the breach that ON x1 does not have are measured:
- whether each lies in a kept body (its r_lake, L_floor and L(r)), or how far from one (cells, by tens);
- inland or ocean-connected (`water_class` on the breached field);
- covered by a final lake, and which;
- its height before the breach.
