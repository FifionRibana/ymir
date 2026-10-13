# F174, Part D — a bibliographic note for option 5 (an ocean-ocean subduction closure)

**Nothing here is implemented.** The note gathers what the four sources the author placed in `docs/refs` say, and proposes a
closure rule with calibration targets and validation criteria, for the author's decision. The author: « Ca demanderait de
faire de la biblio pour ça + calage sur base des publications et expériences. Mais pourquoi pas. »

**The sources read** (all four were readable; none is missing):
- **[S04]** Stern, R. J. (2004). Subduction initiation: spontaneous and induced. *Earth and Planetary Science Letters* 226,
  275–292. doi:10.1016/j.epsl.2004.08.007 (`stern2004.pdf`).
- **[N10]** Nikolaeva, K., Gerya, T. V., Marques, F. O. (2010). Subduction initiation at passive margins: Numerical modeling.
  *JGR* 115, B03406. doi:10.1029/2009JB006549 (`nikolaeva2010.pdf`; cited by paragraph [n], figure and table).
- **[C20]** Crameri, F., et al. (2020). A transdisciplinary and community-driven database to unravel subduction zone
  initiation. *Nature Communications* 11, 3750. doi:10.1038/s41467-020-17522-9 (`crameri2020.pdf`, open access).
- **[SG18]** Stern, R. J., Gerya, T. (2018). Subduction initiation in nature and models: A review. *Tectonophysics* (accepted
  manuscript, doi:10.1016/j.tecto.2017.10.014; `stern2018.pdf`; pages are the manuscript's own numbers).
- Works these four cite (Hall et al. 2003, Toth & Gurnis 1998, Cloetingh et al. 1989, Gerya et al. 2008, Faccenda et al. 2008)
  are quoted **through** them and were not read; they are marked « via ».

## 1. At an ocean-ocean convergent boundary, which plate sinks?

- **Age and density**: oceanic lithosphere « becomes denser than the underlying asthenosphere within ~20–50 million years
  after it forms », and density and mass excess keep growing after that [S04, p. 276]; N10 gives 10–50 Ma [N10, [2]].
- **But not indefinitely, and the strength grows too**: the plate « does not become significantly denser after 70–100 Ma »,
  and older plates are harder to bend [N10, [24]]; « oceanic lithosphere gets stronger as well as denser with age » [S04,
  p. 281].
- **The observed polarity at fracture zones**: the old, dense side sinks beneath the young, buoyant one (Izu-Bonin-Mariana:
  old Pacific lithosphere sinking next to the young Philippine-basin lithosphere) [SG18, p. 26; Gerya et al. 2008 via SG18,
  p. 44–45].
- **Buoyant crust does not sink**: around a plume-built plateau, new subduction zones always dip beneath the plateau [SG18,
  p. 46–47]; a block of buoyant crust entering a trench jams it [S04, p. 278].
- **Reading for C1**: the **older** oceanic side sinks; buoyant (continental or plateau) crust overrides. C1 carries an age
  field (`state.age`), so the rule needs no new state.

## 2. Initiating subduction at the edge of a passive margin

- **Rare in nature**:
  - « Passive margin collapse is a corollary of the Wilson cycle but no Cenozoic examples are known; furthermore, the expected
    strength of the lithosphere makes this mode unlikely » [S04, p. 275]; « extreme skepticism is warranted » [S04, p. 282];
  - passive margin collapse is one of the two SI modes « not well supported by natural observations » [SG18, p. 2]; no
    Cenozoic new subduction zone formed at a passive margin, with one possible exception [SG18, p. 36–37];
  - of the 13 events of the SZI database, only two (Anatolia, Oman) lie close to — and none directly along — a passive margin
    [C20, p. 6].
- **The strength argument** (Cloetingh et al. 1989 via SG18, p. 37): if a passive margin does not turn convergent in its first
  **20 Ma**, the growing strength of its lithosphere makes it increasingly unlikely ever to do so.
- **N10's numerical regimes** (2-D, ~100 experiments, Table 1) [N10, [13]–[16], Fig. 2]:
  1. **stable margin**; 2. **overthrusting** (the continent thrusts over the oceanic plate, which stays attached); 3.
     **subduction** (the mantle lithosphere breaks and a retreating slab sinks).
  - Stable → overthrusting is controlled by the ductile strength of the lower continental crust; overthrusting → subduction by
    the ductile strength of the subcontinental mantle and its chemical density contrast with the suboceanic mantle [N10,
    abstract; [18], Table 4, Fig. 7].
  - **The age of the oceanic plate is secondary** (20–90 Ma tested): it matters only when other parameters are critical [N10,
    [19], Table 5, Figs. 6 and 8].
  - **Only thin, hot continental margins subduct**: Moho temperature > 660 °C [N10, abstract; SG18, p. 42, Fig. 15D] — a
    situation « occasionally created by superimposed external processes such as rifting and/or thermal-chemical plume
    activity » [N10, abstract].
  - **Timing**: subduction is always preceded by overthrusting (0.3–3.6 Ma in the fast mode, Fig. 4); in the slow mode the
    sinking starts after ≥ 15 Ma [N10, [16], Figs. 4–5]. Initiation is dated when the oceanic crust reaches 175 km depth
    [N10, [16]].
  - **Weakening**: stronger plastic crust delays fast-mode initiation by 0.3–2.5 Ma and speeds slow-mode initiation by 4–7 Ma
    [N10, [20], Figs. 13–14]; fluids weaken the slab interface; stronger interplate coupling shortens the thrust sheet and
    postpones initiation [N10, [21], Table 8].
- **Induced initiation** (where plates already converge): a fracture zone becomes a self-sustaining subduction zone after
  **~100 km of convergence**, within 6.8 Myr at **2 cm/yr** with a 40 Myr age offset (Hall et al. 2003 via SG18, p. 44);
  « unless one plate is forced down at a rate of 1 cm/year or more, subduction is unlikely to become self-sustaining »
  (Toth & Gurnis via S04, p. 280); transform collapse needs ~2 cm/yr of convergence to overcome the plate strength [S04,
  p. 286]; induced nucleation begins with compression and uplift, spontaneous with extension [S04, p. 275, 281].

## 3. Frequencies and durations (the SZI database, C20)

- **13 well-constrained events since ~105 Ma** (Table 1, p. 11): Anatolia and Oman 104 ± 1, Aleutian 53 ± 2, Izu-Bonin-Mariana
  52 ± 0.5, Sunda-Java 50 ± 10, Lesser Antilles 49 ± 10, Tonga-Kermadec 50 ± 5, Cascadia 48 ± 5, South Sandwich 34 ± 5,
  Halmahera 15 ± 4, New Hebrides-New Britain 10 ± 2, Philippine 9 ± 3, Ryukyu 6 ± 1 Ma — **about one new subduction zone
  per 8 Myr** over the globe (13 / 105 Myr), clustered at ~55–40 and ~16–6 Ma [C20, p. 5, Fig. 1].
- **Where**: ≥ 9 of 13 at or near a pre-existing convergent boundary; **all within 1 500 km of an active subduction zone,
  mostly within 500 km**; mostly in oceanic plate settings near pre-existing arcs; collisions along a trench often precede
  them; « subduction breeds subduction » [C20, p. 5–6]. Horizontally forced initiation dominates the last 100 Ma [C20, p. 1,
  6].
- **The types** (Table 1): episodic subduction (trench jump after accretion) 5, polarity reversal 4, newly destructive 4.
- **Duration**: from lithosphere-scale faulting to a slab at ~100 km depth, « this process may take up to 5–10 Myr »
  [C20, p. 6–7]; SI propagates along a weakness at ~100 mm/yr [SG18, p. 27].

## 4. When the continent reaches the trench: an Andean margin, or a collision

- **Buoyant crust jams the subduction zone** [S04, p. 275, 278]. Then, with continued convergence (induced nucleation)
  [S04, p. 278–280, Fig. 2; SG18, p. 16, 18]:
  - **polarity reversal**: a new subduction zone forms behind the arc, dipping the other way — the ocean now sinks beneath the
    arrived continent, i.e. an **Andean-type margin** on it (Solomon arc after the Ontong-Java plateau → the New Hebrides
    trench [SG18, p. 18; S04, p. 275]); in continental collisions it accompanies narrow, two-sided collision zones with a weak
    lower crust (Faccenda et al. 2008 via SG18, p. 18, Fig. 5);
  - **transference**: the new subduction zone forms outboard of the collision (Mussau trench) — but « there are no Cenozoic
    examples » that are certain; India has kept colliding with Asia for 50 Myr instead [SG18, p. 18; S04, p. 278].
- **Reading**: polarity reversal is the likely outcome; a continent-continent arrival is a collision (C1's existing C-C
  Davis-Suppe path).

## The proposal: an O-O subduction closure for C1 (a closure, no solver)

It extends `apply_subduction_step` behind a flag of the continent profile; production unchanged.

1. **What is consumed**: at an **ocean-ocean** convergent boundary cell, the plate whose boundary cell is **older**
   (`state.age`) sinks [§1]; ties go to the existing velocity rule. Its oceanic cell loses Δs = rate × convergence × dt, as in
   the C-O rule. The overriding side receives the arc mass on its cells within `arc_distance` — **oceanic cells too** when no
   continental cell is reached (an **island arc**, « un peu de terre dans l'océan »), instead of losing it.
2. **The boundary migrates**: a consumed cell whose s falls under the existing floor threshold is **reassigned to the
   overriding plate** (type unchanged, oceanic), with the existing ≥ 2-neighbour rule against 1-pixel fingers. So a seaward
   band of a continental plate (B′) is eaten from its outer edge when it is the older, sinking side, and the trench walks to
   the continent.
3. **What triggers activity**: a plate pair's O-O boundary segment becomes subducting only after **D_init of cumulative
   convergence** (C1 already tracks convergent plate pairs for the accretion merge). Literature: ~100 km at 2 cm/yr (§2), so
   D_init ≈ 100 km ÷ 7.5 ≈ 13 km on the map (~2 cells of 64²) under the declared layout compression, and a minimum rate of
   ~1 cm/yr below which nothing is consumed.
4. **When the continent arrives**: the band consumed, the continental cells meet the other plate's ocean across the boundary
   — the **existing C-O rule** then consumes that ocean beneath the continent (an Andean-type margin, the polarity
   reversal of §4). A continental cell on the sinking plate never subducts (C1's rigidity, #145) — a collision.
5. **Passive margins stay passive by default**: no spontaneous collapse rule (§2: no Cenozoic example, N10's thin-hot
   condition, Cloetingh's 20 Ma). A margin becomes active only by item 2 (its own older ocean consumed by an induced O-O
   subduction).

**Calibration targets** (from the sources; to be scaled by the declared ×7.5 layout and C1's time mapping):
- D_init ≈ 100 km of convergence, the minimum rate ~1 cm/yr, the initiation duration 5–10 Myr [SG18, p. 44; S04, p. 280;
  C20, p. 6–7].
- The consumption rate: C1's existing C-O `consumption_rate`, unchanged (no new constant).
- Note: C1's history lasts T = 5.6 Myr at B0 (F166's calibration) and up to 300 Myr with the age mapping; a 5–10 Myr initiation
  is the whole history at B0 — the delay matters only at intermediate and old ages.

**Validation criteria** (measured, not tuned on our worlds):
- **(V1) Rarity at passive margins**: ≤ 15 % of the initiations at a passive margin (C20: 2 of 13 near, 0 along).
- **(V2) Proximity**: most initiations within 1 500 km (÷ 7.5 → 200 km on the map) of an existing subduction (C20).
- **(V3) Frequency**: of the order of one initiation per ~8 Myr per Earth-size world (C20: 13 in 105 Myr), at the old age.
- **(V4) Polarity**: the older side sinks on every O-O pair by construction; the check is that arcs land on the younger side.
- **(V5) No artefacts**: no 1-pixel fingers (the ≥ 2-neighbour rule), the arc mass conserved (`SubductionStats`).
- **(V6) The outcome**: at old, a share of the seaward margins turned Andean, the rest passive; W₁₀ high only there.

**What the note does not settle** (for the author): whether a C1 closure change is acceptable at all (it is a change of an
existing closure, behind the flag, not an initial condition), and the time scaling of D_init in C1's step units.
