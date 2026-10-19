# Finding 183 — declared before any measurement (Part 0, the sources)

**Part 0 is complete.** It first stopped Part A because ρ_g was in none of `docs/refs`. The author deposited Sclater & Christie
1980 (`docs/refs/sclater1980.pdf`), and the resumption below sources ρ_g and declares Part A before any measurement. The reference is F182's deposition (`Deposition::declared()`) on F181's default profile.

**The author's decision (verbatim)**: « Oui ok pour ajouter les 3 points. » It covers the three elements (A1 the load, A2 the
compaction, A3 the recycling by subduction) and the amendment of A4: the mass control now counts the subducted term.

## 0.1 — the densities

| quantity | value | where in the code | source |
|---|---|---|---|
| ρ_m (mantle) | 3 300 kg/m³ | `Rebound::declared()`, `crates/ymir-core/src/cascade.rs:465–467` | the F165 declaration (`f165_continent/f165_declared.md`, the brief's values, no page); the same 3 300 is Whitehead & Clift 2009 p. 7 |
| ρ_c (crust) | 2 700 kg/m³ | the same | the F165 declaration (no page). The cascade's Airy mapping (F167, `history.rs:162–165`) uses 2 800 (Whitehead & Clift 2009 p. 7) |
| ρ_w (sea water) | 1 030 kg/m³ | `Airy::declared()`, `crates/ymir-core/src/cascade/history.rs:165` | **Whitehead & Clift 2009, JGR, p. 7, §[24]**: « the density of ocean seawater equal to ρo = 1030 kg m⁻³ » (already F167's) |
| **ρ_g (the silt's grains)** | — | — | **not found** |

- **The search for ρ_g**: the text of every PDF in `docs/refs` (51 files; `pelletier2010.pdf` does not open) was searched for
  « grain density », « sediment density », « solid density », « bulk density », 2.65 / 2.7 g/cm³, 2 650 / 2 700 kg/m³.
  - **Yuan et al. 2019 EPSL Table 1 (p. 6)** gives φ₁,₀ = 0.63 and z₁,₀ = 1 960 m « (Sclater and Christie, 1980) », but no
    density.
  - **Whitehead & Clift 2009 p. 7, §[24]**: « the densities of both continent material and sediment equal to ρc = 2800
    kg m⁻³ ». This is a model's lumped bulk density for compacted sediment, set equal to the crust's, not a grain density.
  - **Braun & Robert 2005 p. 4**: the flexural load's density depends on the process (sediment, water, crust or zero),
    with no value for the sediment. Its « 2600 and 3200 kg m⁻³ » (p. 3) are the upper crust and the asthenosphere.
  - **Lachenbruch & Morgan 1990**: 2.7 g/cm³ is overthrust crust.
  - Sclater & Christie 1980 (JGR 85(B7), 3711–3739) is not in `docs/refs`.
- **The stop (the brief's Part 0 rule)**: ρ_g is needed for the bulk density ρ_b = (1 − φ̄) ρ_g + φ̄ ρ_w of A1. **The source
  is asked for**: Sclater & Christie 1980, the brief's suggestion, to verify. It is Yuan's own source for φ₀ and z₀.
- **A reviewer's memory, NOT a source (labelled, to verify)**: Sclater & Christie 1980's Table 1 is believed to give for
  « shale » φ₀ = 0.63 and c = 0.51 km⁻¹ (1/c = 1.96 km, Yuan's silt values), with a grain density of about 2.72 g/cm³. If
  so, Yuan's « silt » is Sclater & Christie's shale and ρ_g follows from the same table. Nothing here rests on it.

## 0.2 — the compaction (Yuan et al. 2019 EPSL)

- **Eq. 7, p. 3**: φ(z) = φ₀ e^(−z/z₀), z the depth below the sea floor's surface, « related to the weight of the overburden ».
- **Table 1, p. 6**: φ₁,₀ = 0.63, z₁,₀ = 1 960 m for the silt (« d (Sclater and Christie, 1980) »).
- **Appendix A, p. 9–10**: the solid thickness of a column of height h is ∫₀ʰ (1 − φ₀ e^(−z/z₀)) dz = h − φ₀ z₀ (1 − e^(−h/z₀)).
  Yuan inverts it for the mixture with an approximate average porosity. For one class the integral is exact, so A2 needs only
  the inverse of a monotone function: a table, precomputed once.

## 0.3 — the flexure

- **The same as the rebound**: Ĝ(k) = 1 / (1 + (k α)⁴ / 4), α = 64.4 km (Te 25 km, E 100 GPa, ν 0.25, ρ_m 3 300; F165).
  - `Rebound` and `flexural_alpha_km`: `crates/ymir-core/src/cascade.rs:455–476`.
  - The step's filter is applied in the loop: `cascade.rs`, the rebound block after the fluvial stage (F165-C2; the FFT in the
    terrain's frame on the torus).
- The rebound today: Δz = (ρ_c / ρ_m) · (G ∗ e), e the incision's lowering ≥ 0, so an unloading only.

## 0.4 — the land / sea ratio and the reviewer's order of magnitude

- **Measured at F182** (the panel, without the deposition): land 49–62 % of the torus, a land / sea ratio of 0.96–1.63,
  against 0.41 on the Earth (29 / 71 %).
- **The reviewer's estimate** (« sans recyclage, l'océan terrestre se comblerait en quelques 10⁸ ans »): **a hypothesis to
  verify**. `docs/refs` holds no global sediment flux to check it against.

## What Part A will need once ρ_g is sourced (declared now, not built)

- **A1, the load**: the sediment's weight per cell. Under water it is (ρ_b − ρ_w) g h, on land ρ_b g h, with
  ρ_b = (1 − φ̄) ρ_g + φ̄ ρ_w, where φ̄ is the column's mean porosity from A2.
  - **One load balance**: the step's net load change ΔL = (load of the deposits) − ρ_c (the incision's lowering). It passes
    once through G and is divided by ρ_m. The rebound becomes this single term: no double counting, and with the deposition
    off it is the rebound exactly.
  - The water under a new deposit is displaced. Under water (ρ_b − ρ_w) counts the deposit's excess over the water it
    replaces.
- **A2, the compaction**: a cell keeps its solid volume (norm); h = T⁻¹(solid) by a table of T(h) = h − φ₀ z₀ (1 − e^(−h/z₀)),
  computed once, read by linear interpolation; its accuracy will be declared.
- **A3, the recycling**: the sediment on cells of an oceanic plate other than the continent's is removed at each step, in
  proportion to the area of that plate C1 consumes in the step. It is C1's kinematics, with no parameter; the upper bound
  « all the sediment reaching the trench subducts » is declared.

## The resumption (the author's complement)

### R1 — ρ_g, sourced

- **Sclater & Christie 1980**, JGR 85(B7), 3711–3739, **Table A1a, p. 3732** (the PDF's page 22): « Parameters of the
  Exponential Relation Between Porosity and Depth for Various Normal Pressured North Sea Lithologies and the Density Assumed
  for the Sediment Grains in the Same Lithology ».
- Read on the scanned page (the PDF's text layer garbles the digits). The row **« Shale »: f₀ = 0.63, c = 0.51 × 10⁻⁵ cm⁻¹,
  ρ_sg = 2.72 g/cm³**.
  - c = 0.51 km⁻¹ gives 1 / c = 1 961 m, Yuan's z₁,₀ = 1 960 m; f₀ = 0.63 is Yuan's φ₁,₀. The row matches Yuan's « silt »
    exactly.
  - The next row, « Sand » (0.49, 0.27 → 3 704 m, 2.65), matches Yuan's sand (0.49, 3 700 m).
- **ρ_g = 2 720 kg/m³.** The reviewer's memory (« shale », about 2.72 g/cm³) is confirmed by the PDF, not taken from memory.

### R2 — ρ_c: not harmonised here

- The rebound's ρ_c = 2 700 kg/m³ (F165, no page) and the Airy mapping's 2 800 (Whitehead & Clift 2009 p. 7) differ.
- F183 keeps the current values, so that the deposition with its three new flags off is F182 bit for bit. A1 uses the
  rebound's ρ_c = 2 700 and ρ_m = 3 300 (the load passes through the rebound's flexure).
- **Queued**: « harmoniser ρ_c entre le rebond et Airy, sur une source avec page ».

### R3 — the load balance (A1), as the author specified

- **The incision's lowering of a cell takes its sediment cover first, then the basement.** With the cover's thickness H, the
  lowering e splits into e_c = min(e, H) from the cover and e_b = e − e_c from the basement.
- **ΔL = (the load of the deposits) − ρ_b × e_c − ρ_c × e_b**, in the step's net load change of the cover and the basement:
  - the cover's load is computed from its solid volume S and thickness H, both tracked per cell;
  - under water it is S (ρ_g − ρ_w) (= (ρ_b − ρ_w) H, ρ_b = (S ρ_g + (H − S) ρ_w) / H);
  - on land it is S ρ_g + (H − S) ρ_w (= ρ_b H, the pores water-filled);
  - ΔL_cover = the cover's load at the step's end − at its start (the deposits, the cover's erosion, the marine transport,
    the compaction, the subduction), so the ρ_b × e_c term is inside it;
  - the basement's part is − ρ_c e_b.
- The flexure is the rebound's: Δz = − (G ∗ ΔL) / ρ_m, applied once at the step's end on every cell. It replaces the rebound
  (which is the case H = 0: Δz = (ρ_c / ρ_m) G ∗ e). With A1 off, the rebound is F182's exactly (the incision's e at ρ_c).
- **The cover's removed part leaves its solid volume**:
  - with A2, the solid removed is S − T(H − e_c), the column kept in Athy's equilibrium;
  - without A2 it is S e_c / H, proportional.
  - The amended mass control counts it (R5).
- **The erodibility stays the basement's** (the rock law unchanged). **Queued**: the difference between sediment and rock
  erodibility.

### R4 — Part A, the mechanism (declared before building; each behind its own flag of `Deposition`)

- **A1 the load** (`Deposition::load`, ρ_g = 2 720, ρ_w = 1 030 kg/m³): R3. It needs the rebound (its ρ_c, ρ_m, α).
- **A2 the compaction** (`Deposition::compaction`, φ₀ = 0.63, z₀ = 1 960 m: Yuan EPSL eq. 7, Table 1; Sclater & Christie Table
  A1a):
  - T(H) = H − φ₀ z₀ (1 − e^(−H/z₀)), the column's solid thickness (Yuan EPSL Appendix A, one class);
  - every deposit enters at the surface porosity (thickness = solid / (1 − φ₀)), on land, in the lakes and at sea;
  - at the step's end every cover re-compacts, H = T⁻¹(S), and z moves by the change of H;
  - T⁻¹ is a table computed once (uniform in S, 0.5 m steps to 30 km; beyond, the slope 1 / (1 − φ₀ e^(−H/z₀)) at the end),
    read by linear interpolation. Its accuracy is measured and declared (the round trip |T(T⁻¹(S)) − S|);
  - the lakes fill in thickness at the surface porosity (closed form, as F182). The fill compacts at the step's end, so a lake
    can stand slightly under its spill until the next step: a declared consequence of no iteration.
- **A3 the recycling** (`Deposition::recycle`, C1's `Consumption`):
  - at each physics step, on every cell whose 64² parent belongs to a plate with no continental cell at the interval's start,
    a fraction r = 1 − (1 − F_j)^(share) of the cover's solid and thickness is removed;
  - F_j is C1's consumed fraction of that plate over the interval, and share is the step's duration over the interval's,
    so the interval's removals multiply to F_j;
  - **all the sediment of a consumed area subducts** (the upper bound, declared); no part is accreted, for want of a source.
    The continent's plate (its band and its margin prisms included) keeps everything.
- **A4 the amended mass control**, per step, every quantity in solid volume:
  - in = the incision's eroded solid + the hillslope's loss to the sea;
  - ΔStock = Σ over the cells of the cover's solid, measured, end − start;
  - in = ΔStock + the cover eroded by the incision + the cover lowered by the hillslope + subducted + exported
    (+ without A2, F182's porosity term Σ e_c (1 − S / H): the porous cover's height counted as solid, reported apart);
  - checked at 10⁻⁹ relative; F182's flux control (in = deposits + exported) is kept; the subducted volume is reported.
- **The order in the step** (the three flags on):
  1. uplift;
  2. incision;
  3. the cover's erosion;
  4. the fluvial pass and the lakes;
  5. the hillslope and the cover it lowers;
  6. the marine step;
  7. the recycling;
  8. the compaction;
  9. the load's flexure.
  With A1 off the rebound stays where F182 has it (after the incision).

### R5 — the controls (A5, permanent)

- The three flags off: F182 bit for bit (the field and the rebound).
- ρ_g = ρ_w = 0: the cover's load is zero at every step.
- z₀ = ∞: T is linear, so the compaction changes no thickness.
- No consumption (every F_j = 0): no recycling, bit for bit with A3 off.
- Production's hash and the guard.

### R6 — the queue (added by the author)

1. Harmonise ρ_c between the rebound and Airy, on a source with a page.
2. The erodibility of the sediment against the rock.
