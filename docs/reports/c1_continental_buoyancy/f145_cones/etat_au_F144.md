# State at F144 (2026-10-04) — one page, to replace the outdated Part 0 of the v2 note

*Written at F145. The "v2 note" is not in the repository, so this page stands alone; the rest of that note is not
rewritten.*

## What is delivered, what is a candidate

- **Delivered (on by default)**: the témoin is C2 /10 col: the valley construction (Finding 121), `ValleyConstruction::new(k = 0.07183, light pass 0.1 × k_time)`, the absolute slope floor
  S_eq = 0.024·A^−1/2 + A1 (Finding 109), relief-v3 breach, HD drainage, climate, biomes. Guarded: 6 states, field
  and lakes (`bench_field_hashes.json`, `bench_lake_hashes.json`); C2 /10 col = `a8d2d538d692c2f0`.
- **Gated, off by default (candidates):**
  - **the lake base** (F132–F133): `lake_base: InputLakes | InputLakesAndBasins`; the present lakes, and the
    below-sea basins' lakes, are base levels of the χ walk ("ON extended");
  - **the gorge's retreat by age** (F140–F144), `GorgeRetreat` v3 → v5 (v5: the input scope, 18 bodies; the outlet-base
    bound; the capped head falls; the minimal plain), plus the bench options φ = 0, `freeze_design`, `light_mode`
    P1–P3. **ON PAUSE** (the author, 2026-10-04); code committed, viz toggle hidden;
  - the older gated fields (wall profile, confluence band, LTD, off-grid tracé, mur ↔ mer…): see the ADR.

## The gorge on pause (for whoever resumes it)

- **Holds**: the construction leaves no closed hollow; the production plain equals the bench's (1 cm); bodies 1:1;
  F38 never fires.
- **Fails**: G-pits after the light pass; walls at the designed geometry's edge (C2, P1–P3); G-drained; G-sea at ×1.4.
- **Causes, in the base (the « socle »)**: (1) `carve_diag`'s cones laying cells outside their basin (F141-Z, F144-Z);
  (2) the light pass as a second age process (F143-C2); (3) the breach's ramp anchored on each pit's floor (F141 γ).
- P3 is ruled out by its cost (+54 to +189 s per world); P2 was never really tested (d_t degenerated to 1 cell); θ and
  the 82 hanging junctions are to resume.

## The three base defects (the current work, in order)

1. **The cones** (F145): a cell is laid by its nearest sample, even across a divide.
2. **The light pass as a second age**: it makes knickpoints retreat at a grid celerity (F115, Courant 3 699) on top of
   the age the construction already encodes.
3. **The breach's ramp**: from each pit's floor, 0.113 m per cell to the flood's root, without reading the spill
   (`flow.rs:1066–1075`); 292 land cells under the sea in production, under lake 1000016; F38 holds by covering.

## The rules earned since F131

- **Pair a body and a lake by footprint overlap, never by index** (F141).
- **Re-tabulate before building any deviation that touches a tabulated quantity** (F140 / F141).
- **A reviewer's hypothesis enters the ADR labelled "hypothèse à vérifier", never as a fact** (F145; two false lines
  in four rounds, F140 and F143).
- **Every new stage reports its added seconds per world against the reference run_hd 249.8 s**; a second skeleton +
  carve is a red flag before it is proposed (F144).
- Also standing: predictions before measurement, declared non-blind; instruments declared before; rules 11–18; the
  témoin is C2 /10 col; one world one seed is stated as a limitation.

## The author's decisions in force (dated)

- 2026-09-29/30: a present lake is a base level; the below-sea basins' lakes are present lakes; the lakes shrink with
  age; fall or gorge are both valid; falls tagged where relevant.
- 2026-10-02: the gorge's retreat follows the age selector; p ∈ {0, 0.25, 0.5}, default 0.5; steepen when room
  lacks; the minimal plain (« Ok go »); the scope: every input lake outside the below-sea basins; the outlet falls
  capped by the river's size.
- 2026-10-04: **the gorge is paused; the base is treated, starting with the cones.** Not contested: the scope gives
  18 bodies (body 15 meets the criterion on the input); the cap rises to ~240 m on small outlets.
- Standing: **the world's generation time must not rise appreciably** unless a correlation or a feature is added
  (reference 249.8 s).
