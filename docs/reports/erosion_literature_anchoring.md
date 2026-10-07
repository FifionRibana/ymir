# Empirical anchoring of Ymir's erosion constants against four primary publications

**Date:** 2026-09-10 · **Scope:** extraction and comparison only. No code was modified, no run was
launched, no new production value is recommended.

This report answers one question per Ymir constant: **is there a published measurement behind it,
on what terrain, with what uncertainty — or is there not?** Every value carries a page. Where a
publication does not answer the question asked of it, that is recorded as a result rather than
filled in by inference.

---

## 1. Provenance and method

| Publication | Edition read | Pagination cited |
|---|---|---|
| Harel, Mudd & Attal (2016), `10.1016/j.geomorph.2016.05.035` | **ACCEPTED MANUSCRIPT** (52 pp.), *not* the journal version | `PDF p.X / ms p.X−1` — the manuscript prints its own page numbers |
| Montgomery & Foufoula-Georgiou (1993), `10.1029/93WR02463` | Journal version, WRR 29(12) | Journal pages 3925–3934 |
| Whipple & Tucker (1999), `10.1029/1999JB900120` | Journal version, JGR 104(B8) | Journal pages 17,661–17,674 |
| Stock & Dietrich (2003), `10.1029/2001WR001057` | Journal version, WRR 39(4), art. 1089 | `ESG 1-N` page marks **and** AGU paragraph numbers `[n]` |

**Reading method.** Text layers were extracted with `pypdf`. Montgomery 1993 is an OCR'd scan and
Stock & Dietrich's Table 3 is rotated 90° in the PDF; both were re-rendered as images with
`pymupdf` and read visually to verify every number quoted here. Two OCR corruptions in Montgomery
were caught and corrected this way. One extraction error of my own was caught the same way: Whipple
& Tucker's Table 1 is in units of **10⁵ m²**, not 10⁶.

**Values read off figures** are labelled as such. Values absent from a publication are written
*absent from the publication* and are not interpolated.

---

## 2. Headline findings for design

1. **`A_c = 0.1 km²` is not one boundary — the code uses one number for three distinct physical
   transitions**, and the literature places them two orders of magnitude apart:
   - fluvial/debris-flow break in slope–area data: **0.059–0.14 km²** (Whipple & Tucker Table 1)
     and **≈0.1 km²** (Montgomery, Tennessee Valley) → Ymir is *anchored*;
   - systematic departure from the fluvial power law: **1–21 km²** (Stock & Dietrich Table 3,
     28 valleys) → Ymir is *10–200× too small*;
   - field-mapped channel heads, which is what the constant is *named* after:
     **0.001–0.05 km²** (Montgomery Fig. 6) → Ymir is *2–7× too large*.

2. **`K` now has a dimensionally identical published counterpart, and Ymir is 2.25 × 10⁵ above it.**
   Whipple & Tucker Table 2 gives `K (n = 1) = 2.00 × 10⁻⁵ [m^(1−2m) yr⁻¹]` with `m/n = 0.50`,
   i.e. exactly Ymir's exponents and units. Ymir v3 is `4.5 yr⁻¹` in SI. The authors explicitly
   disown their own value as a measurement ("the actual values used are inconsequential").

3. **Pinning `K` does not pin the duration.** Depending on which published `K` is adopted, the two
   shipped erosion steps represent **4.5 × 10⁵**, **1.2 × 10⁷**, or **10⁷–10⁸ years**. Three orders
   of magnitude. The duration is not recoverable from `K` alone.

4. **The stream-power law has a published *upper* slope bound that Ymir does not implement.**
   Stock & Dietrich: the fluvial power law breaks down at `S = 0.03–0.10`, and headwater reaches
   above `S = 0.03` "should be excluded from stream power law analysis". Ymir's `A_c(S)` law runs
   from `S_min = 0.0087` up to the ridges, i.e. straight through the excluded band. Ymir's only
   guard is a *lower* bound, and it is a proxy.

5. **`(m, n) = (0.5, 1)` is a deliberate physical statement, not a numerical convenience** — it is
   the unit-stream-power case with `a = 1`, i.e. erosion strictly linear in dissipation per unit
   bed area. It is also the case that Harel's 1457 ¹⁰Be samples reject (`⟨n⟩ = 2.6 ± 0.4`).

6. **`stream_km2 = 20 km²` is 20–200× above published blue-line first-order channel areas**
   (0.1–1 km², Merritts & Vincent 1989 via Stock & Dietrich).

---

## 3. Comparison table

Verdicts: **anchored** · **in range but at the edge** · **out of range** · **PROXY, no antecedent
found** · **not comparable**.

### 3.A — Comparisons in km² or dimensionless: **not affected** by `geo_scale_ratio`

| Constant | Ymir value | Published range | Source + page | Measurement terrain | Verdict |
|---|---|---|---|---|---|
| `(m, n)` pair | (0.5, 1.0) | exact derivation: unit stream power with `a = 1`, `b = 0.5`, `c = 1` → **`n = a = 1`**, **`m = ac(1−b) = 0.5`** | W&T eq. (10b)(10c), p. 17,664 | theoretical, no terrain | **anchored** theoretically — this is the primary-source derivation |
| `n` | 1.0 | `⟨n⟩ = 2.6 ± 0.4`; median `2.43 ± 0.15`. Theory: `n = a`, with `a` "argued to range from 1 to as much as 7/2" | Harel PDF p.25/ms p.24, p.36/ms p.35; W&T p. 17,663 | Harel: 33 accepted zones of 59, 1457 ¹⁰Be samples worldwide, SRTM 2000 | **out of range** empirically (4σ below Harel's mean); at the bottom edge of the theoretical range of `a` |
| `m` | 0.5 | `m` alone is **absent from Harel** (he publishes only `θ` and `n`); my derivation `θ·n ≈ 1.24`. In W&T, `m = ac(1−b)` | Harel PDF p.23 + p.36; W&T p. 17,664 | — | **out of range (by derivation)** empirically; **anchored** theoretically |
| `θ = m/n` | 0.5 | `0.51 ± 0.12`; `0.35 ≤ m/n ≤ 0.6` predicted; measured **`0.40 ± 0.10`, `0.49 ± 0.10`, `0.41 ± 0.10`**. **But** 0.08–1.40 across 28 valleys, and "*there is no single m/n value that can be extracted from topographic data for a general stream power law*" | Harel PDF p.23/ms p.22; W&T p. 17,664 and Table 1 p. 17,663; **S&D [56], ESG 1-23** | W&T: King Range CA (×2), Central Range Taiwan. S&D: 28 valleys, 14 countries | **anchored**, but the "narrow range" claim is contested by the most recent source |
| `K` ([`RELIEF_V1_K`](../../crates/ymir-core/src/erosion/stream_power.rs)) | 1500 (v1) / 4500 (v3) in km² units → **1.5 / 4.5 yr⁻¹** in SI | **`2.00 × 10⁻⁵ yr⁻¹`** (n = 1, m/n = 0.50); `9.28 × 10⁻⁵` (n = 2/3); `2.00 × 10⁻⁵` (n = 2) | **W&T Table 2, p. 17,666** | *none* — "SS denotes Ks chosen to yield equivalent steady state profiles for U = 2 mm yr⁻¹", and p. 17,666: "the actual values used are inconsequential" | **out of range — factor 2.25 × 10⁵**; but the only published value in the right units is explicitly disowned as a measurement |
| `RELIEF_V1_A_C_KM2` as **fluvial/debris-flow break** | 0.1 km² | **`0.059 ± 0.020` / `0.072 ± 0.024` / `0.140 ± 0.048` km²** measured; W&T's own model uses **10⁵ m² = 0.1 km²** | **W&T Table 1, p. 17,663** (units 10⁵ m²) and **Table 2, p. 17,666** (`x_c = 320 m`, note "10⁵ m² source area") | King Range CA high/low uplift; Central Range Taiwan | **anchored** — inside the measured range, between King Range and Taiwan, and equal to the authors' own modelling choice |
| `RELIEF_V1_A_C_KM2` as **departure from the fluvial power law** | 0.1 km² | scaling transitions measured at **1 to 21 km²** ("Scaling Transitions / Area, km²" column) | **S&D Table 3, ESG 1-7 and 1-8**; method in **[36]** | 28 valleys: 18 US basins at 1:24,000; 10 basins at 1:50,000 (Ecuador, Peru, N. Korea, Argentina, Greece, Kenya, Guatemala, Philippines, Vietnam, Chile) | **out of range — 10 to 200× too small.** S&D note in **[5]** that the literature places the break at 0.1–1 km² and measure something different |
| `RELIEF_V1_A_C_KM2` as **channel head** (its documented name and its use as `head_km2`) | 0.1 km² | field-mapped channel heads at **10³–5 × 10⁴ m² = 0.001–0.05 km²** (read off Fig. 6); the fitted law gives 0.0136 km² at Ymir's own `S_ref` | Montgomery p. 3930, Fig. 6 and eq. (9) | Tennessee Valley CA, field mapping [Montgomery & Dietrich 1992] | **out of range** — 2 to 7× too large |
| exponent of `A_c(S) ∝ S⁻²` | 2 | 2 (theoretical, eq. 7); **1.84 fitted**, `r² = 0.68` | Montgomery p. 3927 eq. (7); p. 3930 eq. (9) | eq. (9): Tennessee Valley channel heads [Montgomery 1991] | **anchored** |
| `stream_km2` (minimal mapped channel) | 20 km² | "*blue-line first order channels occur at drainage areas of ~0.1 to 1 km²*" [Merritts & Vincent 1989] | **S&D [53], ESG 1-22** | northern California coast | **out of range — 20 to 200× too large** |
| navigability thresholds 500 / 5 000 / 50 000 km² | — | none of the four publications addresses navigability | — | — | **not comparable** — outside the subject of all four |
| `CHANNEL_WIDTH_B` | 0.5 | **`b = 1/2`** [Leopold & Maddock 1953]; similar values apply to partially alluviated bedrock channels [Hack 1973]; typical **`b ≈ 0.4–0.6`** | W&T eq. (6) and p. 17,664 | alluvial channels (hydraulic geometry literature) | **anchored** |
| `CHANNEL_WIDTH_A` | 5.0 m·(m³/s)^−0.5 | `k_w` appears in the notation section p. 17,673 with its dimensions `[L^−3b T^b]` and **no numerical value anywhere**. W&T: "*besides the pioneering work of Suzuki [1982], no comprehensive study of the controls on bedrock channel width has been done*" | W&T p. 17,664 and p. 17,673 | — | **PROXY, no antecedent found** — the coefficient is absent from the publication |
| `REFERENCE_RUNOFF_MM = 300` uniform | ⟺ **`c = 1`** in `Q = k_q A^c` | `0.7 ≤ c ≤ 1` | W&T eq. (7) and p. 17,664 | — | **in range but at the edge** (high) — uniform runoff depth is the `c = 1` limiting case |
| `floor/local-ridge ≈ 0.5` (via `iterations = 2`) | ≈0.5 | no valley-floor-to-local-ridge ratio is published in any of the four | — | — | **PROXY, no antecedent found** (confirmed) |

### 3.B — Comparisons that depend on cell size: **conditional** on `geo_scale_ratio`

If `geo_scale_ratio` reaches `cell_km`, the production cell goes from 48.8 m to 366 m. Every row
below compares a grid gradient (`Δz/Δx`) or a per-cell weight and **must be redone**. `K` in
section 3.A is partly conditional too, since its calibration runs through `S`.

| Constant | Ymir value | Published range | Source + page | Measurement terrain | Verdict |
|---|---|---|---|---|---|
| `C = A_c·S_ref²` | `1.10 × 10⁴ m²` | 500 m²; 1000–2000 m²; **2000 m²** working value; 8000–16000 m²; 16000 m²; Fig. 8 sweep 2000–64000 m² | Montgomery p. 3930 ("`A_cr ≥ 2000 tan θ⁻²`"); p. 3931; Fig. 8 p. 3932 | Tennessee Valley CA (steep); Schoharie Creek NY (low gradient) | **in range but at the edge** (high) — 5.5× the working value, and it lands in the *low-gradient* catchment's band |
| `CHANNEL_HEAD_S_REF` | 0.3319 (2048²); 0.1128 (8192²) | average colluvial slope **`0.54 ± 0.11` / `0.36 ± 0.05` / `0.63 ± 0.26`**; channel-head cloud 0.15–0.9; transition slopes 0.02–0.60 (mostly 0.03–0.15) | **W&T Table 1, p. 17,663**; Montgomery Fig. 6 p. 3930; S&D Table 3 | King Range CA ×2, Taiwan; Tennessee Valley CA; 28 valleys | **in range but at the edge** (low) — just inside 1σ of the lowest of the three (0.36 − 0.05 = 0.31); well above S&D's transition slopes |
| `CHANNEL_HEAD_S_MIN` | 0.0087 (tan 0.5°) | **no published lower bound.** But an **upper** bound is published: "*headwater reaches above 0.03 slope should be excluded from stream power law analysis*", and the power law departs at **`S = 0.03–0.10`** | **S&D [55], ESG 1-23**; **[37]**; conclusions **[57]**; abstract ESG 1-1 | 28 valleys + field observation, western US and Taiwan | **PROXY, no antecedent found** for the lower bound — and the publication supplies an *upper* bound Ymir does not implement |
| `RELIEF_V3_DIFFUSION` | 0.08 | Montgomery eq. (3) p. 3926 names `k₁` "an erosional diffusivity" and gives **no value**. W&T and S&D both explicitly exclude hillslopes | Montgomery p. 3926 eq. (3) | — | **PROXY, no antecedent found**, and **not comparable**: the code itself declares this a bare dimensionless weight whose implied `κΔt = dsub·dx²` varies with cell size (`stream_power.rs`, the `DO NOT "FIX" THIS` note) |

---

## 4. The three hard questions

### 4.1 Does `K` mean anything outside steady state?

The answer splits, and both halves are now documented.

**Measuring `K` without steady state: yes.** Harel does it. PDF p.18 / ms p.17, §2.3: "*the
steady-state hypothesis frequently adopted in this framework … is not necessary and has not been
used in the present work. Royden and Perron (2013) demonstrated that even in cases of transient
channel incision, the slope in χ-space will reflect local erosion rates.*" Equilibrium is required
only **per profile segment**, through the piecewise fit of Mudd et al. (2014). Repeated PDF p.35 /
ms p.34.

What the method does require is an **independently measured `E`** — `K̄ = (1/j)Σ Eₖ /
(A₀^{m/n}·M_χ^n̄)`, eq. (10), PDF p.20 / ms p.19. Ymir has no observed `E` and no uplift term, so
the blocker is the absence of `E`, not the absence of equilibrium.

**Every published value, however, is steady-state.** W&T Table 2 footnote b: "SS denotes Ks chosen
to yield equivalent steady state profiles for U = 2 mm yr⁻¹ for all n (for convenience only)."
W&T Table 1 footnote b: "The condition θ = m/n holds if and only if channels are in equilibrium and
both U and K are constants." S&D [56]: "We have not established that each basin has uniform
lithology or is at steady state."

**The number.** `K_Ymir(v3) = 4.5 yr⁻¹` against `K_W&T(n=1) = 2.00 × 10⁻⁵ yr⁻¹` — same units, same
`(m, n)`. **Factor 2.25 × 10⁵.**

**The duration follows the choice of `K`, not the other way round:**

| `K` pinned at | implied duration for Ymir's two shipped steps |
|---|---|
| W&T `2.00 × 10⁻⁵ yr⁻¹` (Table 2, p. 17,666) | **4.5 × 10⁵ yr** |
| Harel's global mean denudation `⟨E⟩ = 242 ± 59 mm/ky` (at an assumed `A = 1085 km²`, `S = 0.01`) | **1.2 × 10⁷ yr** |
| Stock & Montgomery hard-rock end, as cited in `stream_power.rs` | 10⁷–10⁸ yr |

Three orders of magnitude. `k_time = 9000` cannot be factorised into a physical `(K, dt)` pair by
appeal to the literature — the literature does not narrow it enough. The middle row is my own
arithmetic and is linear in the assumed `S`; the assumption is stated so it can be redone at a
measured slope.

### 4.2 Does the break in Flint's law mark the end of fluvial erosion, or the onset of another process?

**The onset of another process. Confirmed by the primary source, with a quantified disagreement on
the slope at which it happens.**

Stock & Dietrich, abstract ESG 1-1: "*we find that this inverse power law rarely extends to slopes
greater than ~0.03 to 0.10, values below which debris flows rarely travel*". Conclusions [57]:
"*This topographic signature is consistent with a fundamentally different valley incision law by
debris flows.*" Paragraph [5] sets out the debate explicitly — hillslope processes (Ijjasz-Vasquez
& Bras, Moglen & Bras, Lague et al.) against debris flows (Seidl & Dietrich, Montgomery &
Foufoula-Georgiou, Sklar & Dietrich) — and the paper comes down on debris flows. Montgomery p. 3932
had already said the same: "*We suggest that the inflection observed in the drainage area–slope
relation derived from DEMs reflects this transition in valley incision processes.*"

**Hergarten's reading is correct.**

The two sources disagree on where: **0.20–0.30** in Montgomery (p. 3929, p. 3932, corroborated by
Seidl & Dietrich 1992 at ~20% in the Oregon Coast Range) against **0.03–0.10** in Stock & Dietrich
([37] and abstract). A factor of ~5. S&D [4] add that they have seen no field evidence of long-term
bedrock river incision (a strath terrace) above valley slopes of 0.05–0.10.

**Consequence for Ymir.** The `A_c(S)` law applies from `S_min = 0.0087` up to the ridges, which
includes the whole 0.03–0.10 band that S&D ask to be excluded from stream-power analysis. The
published guard is an *upper* bound on `S`; Ymir has only a lower one, and it is a proxy.

Montgomery makes a second point that bears directly on the naming of `RELIEF_V1_A_C_KM2`
(abstract p. 3925, conclusions p. 3932): a constant critical support area "*is more appropriate for
depicting the hillslope/valley transition than for identifying channel heads*", and "*the extent of
the contemporary channel network cannot be directly determined from drainage area–slope relations
extracted from DEMs*".

### 4.3 Is `A_c = 0.1 km²` at the bottom edge of the range, or inside it?

**Inside it for two sources, two orders of magnitude below the third — and the dissenting source is
the most direct measurement of exactly this boundary.**

- **Inside, and strikingly so.** W&T Table 1 measures `0.059 ± 0.020`, `0.072 ± 0.024` and
  `0.140 ± 0.048 km²` across three active orogens; 0.1 km² falls between King Range and Taiwan.
  W&T Table 2 uses **10⁵ m² for their own model** — exactly Ymir's value.
- **Central**, in Montgomery: `≈10⁻¹ km²` at both DEM resolutions in Tennessee Valley (p. 3929,
  p. 3932), `≈1 km²` at South Fork Smith River (p. 3928).
- **Ten to two hundred times too small**, in Stock & Dietrich: Table 3 measures the transition at
  **1–21 km²** over 28 valleys in 14 countries, by the most direct method of the four (bracketing
  by the last drainage area in the debris-flow region and the first in the power-law region, [36]).

The three do not measure the same thing. W&T and Montgomery read a **break in slope on the area–slope
data**; S&D read the **onset of systematic departure from the power law**, which occurs much further
downstream. Ymir implements one of them, and the constant's name designates a third — field-mapped
channel heads, at 0.001–0.05 km².

**Three distinct physical boundaries, one number.** This is the substantive design finding.

**Provenance note on the received "0.1–5 km²" range.** It is not in Montgomery 1993 — "5 km²"
appears nowhere in that paper. Harel PDF p.8 / ms p.7 attributes it to **Wobus et al. (2006)**, GSA
Special Papers 398, 55–74, with Stock & Dietrich (2003) cited only for what lies below it. The
second-hand chain via Hergarten (2020) stops one link short of its actual source.

**Scale caveat on Harel.** Only ~3% of his basins are below 5 km², and "*these 45 samples actually
have an area ranging between 0.7 and 1 km²*" (PDF p.34 / ms p.33). Harel says nothing about the
regime below 0.7 km², which is where `A_c` lives.

---

## 5. Width plausibility check

| | value |
|---|---|
| Ymir exponent `CHANNEL_WIDTH_B` | 0.5 |
| published exponent `b` | **1/2** [Leopold & Maddock 1953], typical 0.4–0.6 — W&T eq. (6), p. 17,664 |
| composite width–**area** exponent `W ∝ A^{bc}` | **0.35–0.50** (`b = 0.4–0.6`, `c = 0.7–1`) — W&T p. 17,664 |
| Ymir coefficient `CHANNEL_WIDTH_A` | 5.0 |
| published coefficient `k_w` | **absent from the publication** — dimensions given p. 17,673, value nowhere |
| Ymir max width, 8192² arid | **5.1 m** → `Q ≈ 1.04 m³/s`, i.e. a runoff depth of ≈30 mm/yr over the largest basin |
| same basin (1085 km², measured at 2048²) at the reference runoff of 300 mm/yr | **16.1 m** |
| width predicted by W&T for that basin | **cannot be computed — the published relation has no coefficient** |

The two numbers side by side cannot be produced, and this is not a limitation of the extraction:
W&T say so themselves on p. 17,664 — no comprehensive study of bedrock channel width existed at the
time beyond Suzuki [1982]. What the comparison does establish: **Ymir's exponent is anchored, its
coefficient is a proxy.** And through eq. (10a), `K ∝ k_w^(−a)` — the width coefficient enters `K`
directly, so the two proxies are not independent.

Caveat on the two Ymir figures: the 5.1 m is at 8192² and the 1085 km² at 2048².

**On "dissipation per unit bed area".** W&T p. 17,663 eq. (2b) defines unit stream power as
"*stream power per unit area of channel bed (the product of shear stress and mean velocity V)*",
with `ε = k_b(τ_b V)^a`. Eq. (10c) gives **`n = a`**. Ymir's `n = 1` therefore states, literally,
that **erosion is strictly linear in dissipation per unit bed area**. That is an explicit physical
commitment, not a numerical convenience — and it is the commitment Harel's 1457 ¹⁰Be samples reject.

---

## 6. Open gaps

1. **Harel's supplementary material** holds the per-zone `n` and `K` values. They are the only way
   to tell whether `n = 1` exists at all in the measured population; the accepted manuscript gives
   only the global mean and median. Not in the PDF.
2. **Wobus et al. (2006)**, GSA Special Papers 398, 55–74 — the actual primary source of the
   "0.1–5 km²" range, not read here.
3. **`docs/refs/stock1999.pdf` is present and readable** (Stock & Montgomery 1999). It is the paper
   Harel cites for the `10⁻²–10⁻⁷ m^0.2/yr` erodibility range (PDF p.27 / ms p.26) and the one
   `stream_power.rs` already names for pinning `K`. Note the units imply `m = 0.4`, not Ymir's 0.5.
   Deliberately not read — outside the declared scope of this exercise.
4. **`docs/refs/pelletier2010.pdf` is still a proxy block page** (2985 bytes of HTML), not a PDF.

---

## 7. Reproduction notes

- Native PDF rendering is unavailable in this environment (poppler absent). `pypdf` and `pymupdf`
  were installed into `.venv` for text extraction and page rasterisation.
- Montgomery 1993 is an OCR'd scan: its text layer corrupts `≈10⁻¹ km²` to `10- km2` in two places
  (p. 3929 and p. 3932). Both were verified by rendering the pages at 210 dpi.
- Stock & Dietrich's Table 3 spans two pages rotated 90°; it was rendered at 300 dpi and read
  visually. The "Scaling Transitions" block is two columns: `Area, km²` then `Slope`.
- Whipple & Tucker's Table 1 units are **10⁵ m²**. Reading them as 10⁶ m² inflates every `A_c` by
  ten and inverts the verdict.
