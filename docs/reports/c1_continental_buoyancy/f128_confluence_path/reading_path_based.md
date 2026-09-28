# Finding 128-B — reading, no code: the path-based nondispersive directions (D8-LTD, D8-LAD)

**Sources.** They are in `docs/refs`, **not committed** (rights, as at Finding 124). The text was
extracted with PyMuPDF; the figures were read on page renders. "p." is the journal page, with the PDF
page in brackets.

- **[O03]** Orlandini, Moretti, Franchini, Aldighieri & Testa 2003, *Path-based methods for the
  determination of nondispersive drainage directions in grid-based digital elevation models*, WRR
  39(6), 1144, doi:10.1029/2002WR001639. Pages TNN 1-1 to 1-8.
- **[O14]** Orlandini, Moretti & Gavioli 2014, *Analytical basis for determining slope lines in grid
  digital elevation models*, WRR 50, 526–539, doi:10.1002/2013WR014606.
- **[P08]** Paik 2008, *Global search algorithm for nondispersive flow path extraction*, JGR 113,
  F04001, doi:10.1029/2007JF000964. Pages 1–9.

## 1. The methods, as defined

- **The theoretical drainage direction (TDD)** is Tarboton's (1997) steepest direction over the 8
  triangular facets around the cell: a continuous angle, clamped to the facet's edges ([O03] §2.1.1,
  p. 1-2; [O14] §2.1, p. 527).
- **The two candidates** are the facet's cardinal pointer p1 and diagonal pointer p2, with angular
  deviations α1 = r and α2 = π/4 − r ([O03] §2.1.2, p. 1-2).
  - **LAD** keeps the least angular deviation. With no memory it **is** D8 ([O03] p. 1-2; [O14] §2.2,
    p. 528: "exactly equivalent to … the classical D8 method").
  - **LTD** keeps the least TRANSVERSAL deviation: δ1 = d·sin α1 and δ2 = √2·d·sin α2, the distance
    from the drained cell's centre to the TDD line through the draining cell. The two tie at r = 26.56°,
    not 22.5° ([O03] §2.1.2, p. 1-3).
- **The path-based part.** The signed deviations are ACCUMULATED down the path, and the pointer that
  minimises |cumulative deviation| is chosen ([O03] §2.1.3, eqs. (1)–(4), p. 1-3; [O14] §2.3,
  eqs. (1)–(6), p. 528).
- **Networks** ([O03] §2.2, p. 1-3). Cells are processed in descending elevation. Where paths
  converge, "the cumulative deviation is calculated by considering the path with the largest upstream
  drainage area". [O14] treats a single slope line and defers networks to Orlandini & Moretti 2009b
  (p. 528), which is **not** in `docs/refs`.

## 2. Their parameters — listed, with pages (the addendum's point 1)

| parameter | where | value | status for us |
|---|---|---|---|
| **λ, "dampening factor"** on the inherited cumulative deviation, 0 ≤ λ ≤ 1 | [O03] §2.1.3, eqs. (1)–(2), p. 1-3 (not in the abstract) | **λ = 1** "advocated": the D8-LTD method *is* λ = 1 ([O03] §4 [17], p. 1-8), chosen numerically ([O03] §3.1 [12], Fig. 3, pp. 1-4/1-5) | **not a PROXY: the method's definition.** [O14] §2.3 [9], p. 528: "no boundary conditions nor parameters need to be specified in the D8-LTD and D8-LAD methods" (its eqs. (3)–(4) have no λ) |
| a flat / depression pre-treatment | [O03] §2.2, p. 1-3 ("a recursive procedure … raise the elevations … a small positive slope"); [O14] §3 [12], p. 529 (∇z ≠ 0 required; "algorithms that identify and treat all flat areas and close depressions [Martz and Garbrecht, 1992; Grimaldi et al., 2007]") | not specified by either paper | **OURS** — a declared choice, and the weak point (§5) |

**The reviewer's "damping parameter" exists; the addendum's caution was right to ask for the page.**
It is λ in [O03] only. It is not in the abstract, and [O14] removes it. There is no other parameter.

## 3. What they prove on planar surfaces (the grid bias)

- **[O03]** §2.1.3 [9], p. 1-3: the benefit "can be demonstrated geometrically only if the D8-LTD method
  with λ = 1 and the simple case of a planar slope are considered". LAD with λ = 1 "produces nonlocally
  biased drainage paths as these paths become sufficiently long".
- **[O03] on the synthetics** (Fig. 3, p. 1-5, read on the page render). For the planar slope, from
  λ = 0 (D8) to λ = 1 (LTD), MAE of the relative area error goes ≈ 0.50 → ≈ 0.22 and RMSE ≈ 1.9 → ≈ 0.25
  (ME ≈ +0.13 → ≈ −0.21). For the parabolic valley, MAE ≈ 0.57 → ≈ 0.29. **The spherical mountain and
  crater show no average improvement** (MAE flat in λ; [12], p. 1-4): the mountain is "intrinsically
  dispersive". Along arcs of the mountain the improvement is significant (Fig. 4, [13], p. 1-5).
- **[O14]**, the analytical basis. **Theorem 3.1** (p. 530): for any C² surface whose gradient never
  vanishes, LTD's polygonal path lies in any ε-neighbourhood of the exact slope line for h small enough.
  **Lemma 3.3** (p. 530): the cumulative transverse deviation is bounded, |δ⁺(k, h)| ≤ h, and it is "the
  key factor" (Remark 3.7). Remark 3.8: the same proof does not carry to LAD. Numerically (Table 2,
  p. 536), mean |deviation| is D8 7.4 m against LTD 0.00–1.58 m for h ≤ 5 m; LTD outperforms D8 up to
  h = 20 m.

## 4. Their declared limits

1. **Not local** ([O03] [3], p. 1-1): "this strategy does not eliminate the bias at the local level, it
   provides nonlocally constrained drainage paths". Every step is still one of 8 directions. **At a
   1-step chord LTD and D8 are identical by definition (R8 = 1.000 for both)**; any effect is at longer
   chords. The addendum's point 2 is right, and a gate read at one chord is worthless.
2. **Divergent terrain** ([O03] [12], [13]): no average gain on the spherical mountain and crater.
3. **DEM errors** ([O03] [14], Fig. 5a, p. 1-6): with σz/Δx < ~5 % on the spherical mountain, D8 gives
   BETTER average areas than LTD; above ~5 % they are comparable.
4. **Resolution** ([O14] §5 [52], p. 536): the grid must describe the topography (h ≤ 20 m on their
   surface). "Further research is needed".

## 5. Paik 2008's critique, and whether 2014 answers it

- **(a) No analytical basis for λ** ([P08] §1 [8], p. 2): D8-LTD's "complicated algorithms require users
  to specify … the 'dampening factor' … recommended … 1 on the basis of several empirical tests. However,
  there is no analytical background". **[O14] answers it**: it removes λ and gives Theorem 3.1. It cites
  Paik for exactly this in [3], p. 527: "as pointed out by Paik [2008], an analytical basis … has not been
  provided so far".
- **(b) Asymmetry and false confluences on divergent terrain** ([P08] §3 [31], p. 6; Figs. 4f and 6e).
  On an outward cone D8-LTD's pattern is "radially asymmetric". Even the paths from the nine highest
  cells, which "should clearly follow eight straight directions", deviate. LTD "suffers the occurrence of
  unrealistic confluences". **[O14] does NOT answer (b)**: it treats one slope line on a C² surface with
  non-vanishing gradient, and a cone's tip is neither. Networks are deferred to 2009b.
- **Paik's numbers** ([P08] Table 1, p. 7): cumulative lateral deviation, LD8 = 100:

  | surface | ED8 2nd | ED8 3rd | GD8 | D8-LTD |
  |---|---|---|---|---|
  | plane (1:4) | 29 | 16 | 4 | **4** |
  | cone | 64 | 55 | 51 | **60** |
  | inward cone | 70 | 61 | 59 | **63** |

- **Paik's competing method, GD8** ([P08] §2, pp. 3–4). Still one receiver per cell. It switches to the
  "secondary" direction when the deviation is one-sided and the global slope to a "start cell" is
  steeper. It has no parameter, but it has a search extent.
- **Paik's own caveat** ([P08] [9], p. 2): "the problems associated with D8 are inherent in the nature of
  digitization; hence, no absolute solution exists".

## 6. `FlowConfig::dinf` — its history, read whole, and whether LTD falls into the same trap

**Commit `4f47495`** (2026-06-29): "D-infini (Tarboton) sur les plats — ESSAYÉ, gated OFF, REPLI D8 (le
tracé rendu re-quantise)". The doc at `flow.rs:45` says the same. What was tried: D∞ **on FLATS only**,
for the ROUTING and the RENDERED river. The continuous angle came from the gradient of `flat_grad`, and
the flow was split between the two bracketing D8 neighbours. Result (`probe_dinf_compare`): WORSE.
Diagonal share on large flats 67 % → 1 %; local R2 0.40 → 0.91. Two causes, both written down:
1. **Re-quantisation.** The rendered river follows the primary receiver, one D8 neighbour per cell:
   "re-quantise quoi qu'il arrive". The continuous benefit was in the ACCUMULATION, not in the trace.
2. **A cardinal flat gradient.** `flat_grad` is the Garbrecht–Martz distance to the outlet (8-connected
   BFS ≈ Chebyshev) plus a bounded noise (`FlatPerturbation`). Its gradient is CARDINAL, so the
   continuous angle points cardinally.

**Would LTD fall into the same trap?**
- **(1) No.** LTD never projects a continuous angle onto a cell: it IS a one-receiver D8 path, and its
  gain is the memory along the path (Lemma 3.3). This is the part of the reviewer's "D∞ re-projected
  without deviation memory — which LTD corrects" that holds.
- **(2) Yes, if flats are resolved as they are now.** On a flat the TDD is undefined, and [O14]'s
  theorem requires ∇z ≠ 0. If flats are resolved by `flat_grad`, LTD's TDD there points cardinally, and
  its memory makes it follow that cardinal direction faithfully: the same comb. LTD corrects cause (1),
  not cause (2).
- **So, for C:** on cells where the resolved surface has no strictly descending facet, the direction is
  the existing D8 of `compute_flow` (Garbrecht–Martz + perturbation), with the cumulative deviation
  RESET. Those cells are counted, declared, and not claimed to be unbiased.

## 7. The choice, motivated

**D8-LTD, in its parameter-free form ([O14] eqs. (1)–(6), i.e. [O03] with λ = 1), for the SKELETON'S
directions only.** The routing of the hydrology stays D8.
- It is the only one of the two with an analytical basis: [O14] Theorem 3.1, with Remark 3.8 excluding
  LAD.
- LAD with memory is biased on long planar paths ([O03] [9]).
- The sources give no support for the reviewer's reason ("LAD too sensitive on flats"). Their reason is
  the transverse deviation's geometry ([O03] [16], p. 1-8; [O14] Remark 3.8).
- **λ = 1 is the method, not a PROXY.** OURS, labelled: the flat threshold (a facet counts as descending
  if its slope is > 0), the flat fallback, and the pit-free surface. The pit-free surface is the breached
  field the D8 skeleton already uses, so the three trees of C share one surface.

**Do these methods do what the round needs? Yes, with two named caveats, so B does not stop the round.**
They remove the D8 path-level grid bias on planar and C² surfaces ([O03] Figs. 2–3; [O14] Theorem 3.1,
Table 2). The caveats:
- **(i) flats**, the §6 trap: resolved as the D8 skeleton resolves them, and counted;
- **(ii) divergent terrain**: [P08] predicts asymmetry and false confluences on a cone. Finding 126's
  synthetic is a DOME, which is a cone, so C's synthetic is exactly Paik's hard case, and his
  cone / plane contrast (60 against 4) is the expectation to read against.

**C0's reproduction target, declared here (the addendum's point 3).** The planar slope of [O03] and
Paik's plane.
- **Geometry.** [O03] does not state its plane's orientation or size in the text. I take a 1:4
  lateral:longitudinal gradient (TDD at arctan ¼ = 14.04° from the cardinal axis, Paik's plane), draining
  to one edge, on 30 × 30 cells (the cell count visible in [O03] Fig. 2a) and on 60 × 60 (Paik: the value
  depends on the domain size).
- **Theoretical areas.** The analytic strip of each cell's drainage lines through its corners, as [O03]
  §3.1 [11] defines them.
- **Pass (the order of magnitude of their improvement):**
  - LTD/D8 MAE ratio in [0.2, 0.9] ([O03]: ≈ 0.44);
  - RMSE ratio in [0.05, 0.4] (≈ 0.13);
  - Paik's cumulative lateral deviation ratio in [0.01, 0.15] (0.04).
- **Otherwise stop: the code, not the method.**
