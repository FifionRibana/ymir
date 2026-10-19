# Finding 182 — declared before any measurement (Part 0, the inventory; Part A, the resumption)

**Part 0 stopped Part A because G was not sourced (0.2). The author's complement resumes it: G is sourced (Guerit et al. 2019),
the route and the marine scheme are the author's decisions. Part A below is declared before any measurement.** The reference
profile is F181's default (`Margin::continent_default()`, the intermediate age 100 Myr).

## 0.1 — where the eroded matter goes today (the physics level, 256²)

- **The level's loop** (`crates/ymir-core/src/cascade.rs`, `run_level`):
  - `:545–551`: the sea is held at sea level during the level (`drive.hold_sea`);
  - `:619`: the erosion is `incise_with_floor_torus` (one iteration per step, the talus and the diffusion off inside);
  - `:657–676`: the hillslope sub-step (`talus_sweep_torus` `:657`, then `linear_diffusion_torus` `:669`), after which
    **every sea cell is reset** to sea level (`drive.hold_sea`), or to its pre-hillslope value: « the talus may have shed onto
    it » (`:670`).
    Whatever the hillslope moved onto the sea is deleted.
- **The incision** (`crates/ymir-core/src/erosion/stream_power.rs`, `incise_with_floor_torus`, `:626`):
  - `:645`: the routing on the current surface, `compute_flow` (`terrain/flow.rs`): **priority-flood depression filling**
    (`flow.filled`, the spill surface), D8 receivers on the filled surface, accumulation. **The real field is not filled**:
    the lake exists only in the routing;
  - `:664`: a cell with no receiver (`DIR_NONE`) or at or below the sea level is a **fixed base node** (`receiver[k] =
    k`). The sea cells take no incision and receive nothing;
  - `:682`: the topological stack (receivers before donors) via a donor CSR tree;
  - `:716 ff.`: the **implicit Braun & Willett incision** in stack order, from base level up (E = K A^m S^n). The base level
    bound (F80) and the closed-depression option (F96 A1, `depression_floor`, `false` in the cascade's solver) can stop the
    incision, never deposit;
  - **the eroded matter is not tracked**: detachment-limited. What a step lowers disappears, as the C-3 module's doc says
    (F150). At the mouth nothing happens: the river's last land cell grades to its sea receiver, which is fixed.
- **The lakes**: a closed depression is only the filled routing surface. Its floor keeps its height, so nothing fills it.
  Hence the large closed lakes and the « dust » of small lakes at the physics level.

## 0.2 — the sources (both in `docs/refs`)

- **Yuan, Braun, Guerit, Rouby & Cordonnier 2019** (JGR Earth Surface, doi:10.1029/2018JF004867, `yuan2019_2.pdf`, the
  accepted manuscript):
  - **the fluvial deposition** (p. 5, eq. 8, after Davy & Lague 2009): ∂h/∂t = U − K′ q_w^m′ S^n′ + (d* v_s / q_w) q_s. Eq. 9:
    the transport length L_f = q_w / (d* v_s) = p̃ A G / W, with **G = d* v_s / p₀** a dimensionless deposition coefficient.
    In discrete form (eq. 13–15, p. 7) the deposition term is (G / p̃) / Ã_i × Σ_upstream (h^t − h^{t+Δt} + U Δt);
  - **the numerical method**: implicit, with a **Gauss-Seidel iteration** (eq. 19–21, p. 7–8) until max|Δh| < 10⁻³ U Δt;
    **up to ~12 iterations** (Fig. 2, p. 8), the count depending on G, not on the resolution. It is proven convergent for
    G/p̃ ≤ 1 (p. 8). With G = 0 it is FastScape with no iteration (p. 9). **It does not hold in one pass when G > 0**;
  - **the values of G**: « for typical values of the various parameters (p = 10⁻⁷ m/s, v_s ∈ [10⁻⁶ − 10⁻¹] m/s), G/p̃ is of
    order 1 or greater, in good agreement with estimates from natural sedimentary systems **(Guerit et al. (2018), pers.
    comm.)** » (p. 6). The simulations use G = 0, 0.1, 1, 10 (p. 13–15) as a sensitivity, not a calibrated value.
  - **So G is not sourced as a value here**: only an order of magnitude (≥ 1, from Davy & Lague's ranges), and an unpublished
    personal communication.
- **Yuan, Braun, Guerit, Simon, Bovy, Rouby, Robin & Jiao 2019** (EPSL 524, 115728, `yuan2019_3.pdf`):
  - p. 3: « we recently implement this procedure in a new erosion model that accounts for sediment deposition in the
    continental system (**Guerit et al., 2019**; Yuan et al., 2019). However, in this contribution, we do not include
    continental deposition »;
  - **the marine model** (p. 3, eq. 5): ∂h/∂t = K_M ∇²h + Q_S for h < h_sea. With two grain sizes (eq. 6, after Rivenæs
    1992): silt K₁, sand K₂, an active layer L = 100 m, the silt fraction of the source f = 0.5, the surface porosities
    φ₁,₀ = 0.63 (silt) and φ₂,₀ = 0.49 (sand), and the compaction e-folding depths 1 960 and 3 700 m (Sclater & Christie
    1980), all in Table 1, p. 6;
  - **the coefficients**: Table 1 gives the ranges K₁ = 10–1 000 m²/yr and K₂ = 1–1 000 m²/yr. The Ogooué inversion gives
    « the best-fit values of the transport coefficients for silt in the marine domain are in the range of 300−500 m²/yr, in
    agreement with previous studies on offshore diffusion » (abstract, p. 1);
  - **its method**: an ADI implicit scheme with iterations on the silt fraction (p. 4, « Iterations are performed until
    convergence is reached »).
- **The verdict of 0.2**:
  - **G is not sourced.** The deposition coefficient rests on Guerit et al. (2018, pers. comm.), which the EPSL paper cites
    as Guerit et al. 2019 (the continental deposition model). **By the brief's rule, Part A is stopped until the source is
    deposited** (Guerit et al. 2019, *Geology*, presumably « Fluvial landscape evolution controlled by the sediment deposition
    coefficient »: to verify).
  - The marine silt coefficient is sourced (300–500 m²/yr, the Ogooué inversion); the sand coefficient only as a range.
- **The shelf break, sourced** (for Part B, when it runs): Harris et al. 2014 (Marine Geology 352, 4–24, p. 6), after IHO
  2008: the shelf extends « from the low water line to a depth at which there is usually a marked increase of slope towards
  oceanic depths ». The shelf break is « the line along which there is marked increase of slope at the seaward margin of a
  shelf ».
  - No depth is given. Declared: along each seaward path from a coast cell, the break is the cell of **the largest increase of
    the seaward gradient**; the shelf width is the distance from the coast to it.

## 0.3 — the ×7.5 compression

- **The cascade's convention** is written in F162-0: « the relief is at the real 400 km (F67); the ×7.5 is for the exported
  hydrology only; every cascade area is real ». The permanent test `the_cascade_reads_no_signified_quantity` enforces it.
  The cascade's process coefficients apply **unconverted** on the 400 km world:
  - the stream power K (m^{1−2m}/yr) on real areas;
  - the hillslope diffusivity D (m²/yr, `cascade.rs:146–148`);
  - the flexure α (Te 25 km, `cascade.rs:463`) on the real cell (`cascade.rs:636–640`).
- **Declared**: the marine K_M (m²/yr) and the fluvial deposition follow the same convention: **no conversion**, applied on the
  256² cell of 1.5625 km.
- **A tension to flag** (not resolved here): the tectonic geometry since F174 (the margin's widths, D_gap, the shelf reference
  « Harris ÷ 7.5 ») treats the world as compressed ×7.5. Under F162-0's convention the shelf a real K_M builds would be compared
  with Harris unconverted (the real widths), not ÷ 7.5. This is the author's call.

## 0.4 — the stability

- **Explicit diffusion**: Δt ≤ Δx² / (4D), with Δx = 1 562.5 m at 256²:

  | D (m²/yr) | Δt limit |
  |---|---|
  | 1 (the sand range's bottom) | 610 kyr |
  | 300 | 2.03 kyr |
  | 500 | 1.22 kyr |
  | 1 000 (the ranges' top) | 0.61 kyr |

- **The history's step durations** (`history::step_durations`, T = 5.623 Myr, 300 steps):

  | age | shortest step | longest step |
  |---|---|---|
  | young, 20 Myr | 18.7 kyr | 163 kyr |
  | intermediate, 100 Myr | 18.7 kyr | 1.47 Myr |
  | old, 300 Myr | 18.7 kyr | 5.70 Myr |

- **The verdict**: one explicit step per history step is **unstable at every age** for the sourced silt coefficient. The
  ratio is 9–15 at the shortest step and ~1 200 (intermediate) to ~4 700 (old) at the longest. An explicit scheme stays a
  closure only with sub-steps at the stability limit: Σ steps / Δt_lim ≈ 82 000 sub-steps at the intermediate age, 245 000 at
  the old one, on the marine cells. That is stable and explicit, but it is not « one pass per step ». The implicit route (a
  linear system) would be a solver. **This is the author's call, together with G.**

## The route (to be chosen once G is sourced)

- Yuan's complete method iterates (Gauss-Seidel fluvial, ADI marine), so **route 1 is not one-pass**.
- The preferred route (operator splitting: the implicit incision as now, then one upstream-to-downstream pass depositing
  G·Q_s/(p̃A) with the current step's quantities) treats Yuan's coupling explicitly. Yuan's stated reason for the implicit
  coupling is stability at large Δt (p. 2). Whether the explicit splitting contradicts the source's physics at our steps
  (18.7 kyr – 5.7 Myr) is to be shown once G is known.

## Part A — the resumption (declared before any measurement)

### A0 — the scale (the principle confirmed by the author)

The tectonic structures read at the signified scale (÷ 7.5); the surface processes (erosion, diffusion, deposition) run at the
real scale of the 400 km map, with no conversion. Every coefficient below is applied unconverted on the 256² cell of
1 562.5 m. (The author's verbatim answer was not included in the complement; the principle is recorded as stated there.)

### A1 — G, sourced

- **Guerit, Yuan, Carretier, Rouby, Braun et al. 2019**, *Geology* 47(9), 853–856 (`docs/refs/guerit2019.pdf`, the published
  version), p. 855, on 68 natural systems: « Excluding these systems, G ranges from ~0 to 3.1 ». The excluded systems are three
  with G → ∞. About one third of the natural systems have G < 0.4. G was estimated with n = 2 ± 0.5: « A more linear
  relationship between slope and erosion would lead to higher values of G but in the same order of magnitude » (p. 855).
- **The rule (the author's)**: a range gives the median as the value and the bounds as the sensitivity. **G = 1.55**,
  sensitivity **G = 0** and **G = 3.1**. (The experiments' 0.4–1.9 is not the natural range and is not used.)
- p̃ = 1: a uniform precipitation equal to the reference p₀ (Yuan's p̃ = p / p₀, JGR p. 5).

### A2 — the mechanism (`crates/ymir-core/src/cascade/deposition.rs`, behind `Drive::dep`)

The step at the physics level, with the deposition on (each item is the code's order):

1. **The uplift**, unchanged.
2. **The incision**, unchanged: the implicit Braun-Willett step with the rock law (`incise_with_floor_torus`, then the
   cover/basement correction). Its eroded matter per cell is e_i = before_i − after_i (≥ 0).
3. **The fluvial pass** (A1–A2, one pass, no iteration):
   - the routing is the incised surface's: `compute_flow` (the priority flood, D8 steepest descent on the filled surface,
     periodic), the same routing the incision builds, read on the incised field. The incision's MFD area (p = 6) is the
     incision's only; the pass's Ñ_i is Yuan's eq. 14 count, the cells of i's D8 catchment, i included;
   - upstream → downstream (the reverse of the receivers-first stack): Q_in(i) = Σ over i's donors of their outflow;
     **D_i = min((G / Ñ_i) · Q_in(i), Q_in(i))** (Yuan JGR eq. 13–15 with the step's quantities, and the author's cap: a cell's
     deposit never exceeds the flux entering it); the outflow Q_in + e_i − D_i goes to the receiver;
   - **the lakes (A2)**: a lake cell is a land cell with `filled > z` on the routing surface. It takes no fluvial deposit; its
     inflow and its own e pass along the lake's routing to its outlet (the lake cell whose receiver leaves the lake). At the
     outlet, the sub-basin draining to it fills **from the bottom up by volume**: the level L holding the arrived volume, at
     most the spill (`filled`). The remainder continues downstream;
   - a sea cell receives the flux (the marine source). A land cell with no receiver keeps what reaches it (« stranded »).
4. **The rebound** reads the incision alone (e_i), as before F182. **The deposits' load is not flexed**: a declared
   simplification, consistent with the rebound's existing max(0, e).
5. **The hillslope** (the talus, the linear diffusion) and **the sea reset** are unchanged. **What the hillslope loses to
   the sea is no longer deleted**: H = Σ_land (z before − z after the hillslope) becomes a marine source. It is spread over the
   coastal sea cells in proportion to Σ over their land 4-neighbours of (z_land − z_sea)⁺, the gradient that drives both
   operators across the coast.
6. **The marine step** (A3), on the step's sea cells:
   - one class, the silt, **K_M = 400 m²/yr** (sensitivity 300 and 500), surface porosity **φ = 0.63** (Yuan et al. 2019 EPSL,
     abstract and Table 1). A solid volume V deposits a thickness V / (1 − φ);
   - **the sand is omitted** for want of a sourced value (Table 1 gives only the range 1–1 000 m²/yr): a declared
     simplification, so eq. 6 reduces to one class (F = 1), ∂h/∂t = K_M ∇²h + Q_S;
   - **explicit sub-steps**: n_sub = ⌈4 K_M Δt / Δx²⌉, Δt_sub = Δt / n_sub ≤ Δx² / (4 K_M) = 1.53 kyr at K_M = 400
     (λ = K_M Δt_sub / Δx² ≤ 1/4: the explicit scheme's stability bound). The source enters in equal parts at each sub-step.
     No iteration, no linear system;
   - fluxes only between sea cells (no flux across the coast); periodic on the torus;
   - **the transport moves only the deposited sediment** (each donor's outgoing flux is limited to its sediment thickness).
     This is the one departure from eq. 5, which diffuses the whole surface h: on C1's bathymetry an unlimited K_M would
     erode the tectonic relief (the trench, the slope) over √(K_M A) ≈ 200 km at 100 Myr, half the map. Yuan's Ogooué setup
     does not discuss it. Declared, so it can be reversed;
   - compaction (Table 1's e-folding depths) is ignored; a marine sediment later emerged and eroded is eroded as solid
     (its volume reported).
7. **The sediment thickness** (the viz layer): the deposits of 3 and 6, less what the incision and the hillslope lowered on
   the cell since. Land and lake deposits carry no porosity (Yuan JGR applies none).
8. **The mass control (A4)**, per step in f64: e (the incision) + H = fluvial + lakes + stranded + marine solid + exported.
   Exported is 0 on the torus with the marine part on (nothing leaves the map). The marine solid is measured on the field
   (Σ of the sediment change over the sea cells × (1 − φ)), not booked. The worst per-step relative residual must stay
   < 10⁻⁹. The f32 field's residual (the deposits as written in z against the pass's f64) is reported apart.

### A3 — the reference, in the bench only (Yuan's iterative coupling)

- **Yuan JGR eq. 19–21**, coded in `tests/f182_cascade.rs` and passed through `Drive::dep_hook` (a replacement of the fluvial
  stage; the cascade's code holds no iteration).
- Each iterate is eq. 20's form: the implicit incision (the cascade's own, the rock law included) of the field plus the
  previous iterate's deposit; the deposit from that iterate's erosion by the cascade's pass (the same cap, the same lakes,
  routed on the iterate's surface without its lake fill). The first iterate is the splitting.
- **The stop** is eq. 21: max |h^{k+1} − h^k| < 10⁻³ U Δt, with U Δt the step's mean |U Δt| (norm). At most 50 iterates; the
  count is reported per step. Yuan proves convergence for G/p̃ ≤ 1 only (p. 8); G = 1.55 is above, so a non-convergence is
  reported, not hidden.
- The output is made mass-exact: z = before − e^K + D^K (it differs from Yuan's h^K by D^K − D^{K−1}, within the tolerance).
- **Run** on 2 panel seeds at the intermediate age. **Compared, with no threshold**:
  - the cumulative land deposit maps (fluvial + lakes): Pearson correlation and the mass difference;
  - the masses in the land, the lakes and the sea;
  - the lake area;
  - the median altitude;
  - P7: the gap per step against the step's duration (the longest 1.47 Myr, the shortest 18.7 kyr).

### A4 — the cost

The marine sub-steps scale with the history's duration: ~65 500 sub-steps at 100 Myr, ~13 100 at 20 Myr, ~196 600 at 300 Myr
(K_M = 400). The measured cost per world at the three ages is announced before the tiers run.

### A5 — the measures (Part B)

- **Lakes** (on the final 256² field): a lake is an 8-connected component (on the torus) of land cells with `filled > z`.
  Size classes: **dust 1–2 cells** (≤ 4.9 km²), 3–15, 16–63, **large ≥ 64 cells** (≥ 156 km²). Count, area, distribution.
- **The shelf** (the definition of 0.2): from each coastal sea cell (a sea cell of the main ocean with a land 4-neighbour),
  the seaward path follows the sea floor's D8 steepest descent until no neighbour is lower (at most 128 cells). Along it the
  gradient of each link is (z_j − z_{j+1}) / its length. **The break is the start of the link with the largest increase of
  gradient over the previous link**; the width is the path length from the coast to it. The median per seed, for passive and
  active coasts (A1's classes of the 64² parent cell or its 8 neighbours).
  - Reported against **Harris unconverted** (Table 2, p. 11: passive 88.2 km for all oceans, 34.9–123.0 km per ocean; active
    31 km) **and Harris ÷ 7.5** (4.65–16.4 km passive, 4.1 km active). **No criterion on the width**: the author judges on the
    plates.
- **Kapos**: plains, plateaus, hills, mountains (`macro_stats`), against Europe's plains 51.6 %; the median and the deciles.
- **Coasts**: the dimension (`coast_torus`, box size 4), the length ratio (with / without deposition), the deltas: components
  of land gained at the coast (land with deposition, sea without), their count and area.
- **The sediment**: the thickness map; the shares (land, lakes, sea) from the budget.
- **The protection**: condition (v) (collision is not the only formation), W₁₀ global / active / passive, the circumnavigation.
- **The mass control** per world.

### A6 — the stops (the brief's)

The mass control fails; any iteration or solver without consent (the reference is the bench's, with consent); an unsourced
constant or an undeclared scale conversion; an unstable explicit scheme; the circumnavigation lost on ≥ 2 seeds; (v) fails.
