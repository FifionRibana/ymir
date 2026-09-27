# Finding 120 — reading (Génevaux et al. 2013) and the laws a valley construction would use

Written for the author before any primitive exists. Every number carries its page; where a source
could not supply a number, the gap is stated instead of filled from memory.

> ⛔ **Status (2026-09-25).** Section 2 was read WHILE block A0 was running, before A0 returned —
> a deviation from the round's order, declared. A0 then fired the round's stop rule (the trunks move
> between A1+B2 OFF and ON: p90 37–537 cells at A ≥ 100–1000 km²), so **these laws are read and not
> used**: no primitive has been written. The author's resolution amendment (physical units only) is
> applied throughout: every law below is in km², m, m/m or degrees.

Sources read in full or in the cited sections, all stored in `docs/refs/`:

| file | reference |
| --- | --- |
| `genevaux2013.pdf` | Génevaux, Galin, Guérin, Peytavie, Beneš 2013. *Terrain Generation Using Procedural Models Based on Hydrology.* ACM Trans. Graph. 32(4), Article 143, 10 pp. Open copy from Purdue CGVLab. |
| `clubb2022.pdf` | Clubb, Weir, Mudd 2022. *Continuous measurements of valley floor width in mountainous landscapes.* Earth Surf. Dynam. 10, 437–456. Open access. |
| `harel2016.pdf` | Harel, Mudd, Attal 2016. *Global analysis of the stream power law parameters based on worldwide ¹⁰Be denudation rates.* Geomorphology 268 (accepted manuscript; page numbers below are PDF pages). Already in the dossier. |
| `whipple1999.pdf` | Whipple, Tucker 1999. *Dynamics of the stream-power river incision model.* JGR 104(B8), 17,661–17,674. Already in the dossier. |

## 1. Génevaux et al. 2013

**What it is.** A procedural generator that builds the river network FIRST and the terrain SECOND,
as a continuous analytic function stored in a construction tree (143:1, abstract; 143:2–143:3, §3).

**Inputs** (143:2, §3; 143:3 Fig. 2):

- a domain contour, river mouths, optional sketched river parts;
- two user maps: a **river slope map** (drives the network's dendritic shape and the Rosgen type)
  and a **terrain slope map** (where mountains and flatlands go);
- growth parameters: edge length `e`, the elevation window `ζ` for node selection (143:4), rule
  probabilities `Pc, Ps, Pa` (143:4), distance limits `η = σ = 3/4 e = 1500 m` (143:5), a Lipschitz
  slope bound `κ` (143:5).

**Pipeline.**
1. **Network growth from the mouths upstream** (143:3–143:5): candidate nodes seeded on the contour,
   expanded by a grammar whose priority index IS the Horton–Strahler number (143:4, "Index priority
   and Horton-Strahler's number"; rules Table 1, 143:5). Each new node is higher than its parent, and
   the Lipschitz condition `|p_z − p_zi| < κ(p)·d(p, p_i)` "prevents the creation of huge cliffs"
   (143:5).
2. **Watersheds** from the Voronoi cells of the river nodes (143:3; 143:5 §5.1). Flow from an
   empirical law `φ = 0.42·A^0.69`, `A` in m², `φ` in m³/s, after Dunne & Leopold 1978 (143:6).
3. **Crests**: `q_z = max(a_z, b_z, c_z) + λ(q)·d`, with `λ ∈ [0, 0.25]` a terrain-slope magnitude
   (143:6, Fig. 11). **The crest rises at most 0.25 m/m (14°) above its three river nodes.**
4. **Rosgen classification** of each reach (types A+ … G), which fixes its trajectory and its bed
   profile (143:6, §5.2, Fig. 9).
5. **Primitives and construction tree** (143:6–143:8): river primitives `h(p) = u_z(p) + δ(d(p))`
   with `δ` a 1-D cross-profile per Rosgen type (143:7); terrain primitives on a Poisson-disk
   sampling, 50 per cell, whose elevation interpolates between the river projection and the ridge
   projection (143:7); Perlin noise with amplitude set by the distance to the river and the relief
   (143:7). Blend and replace operators (143:8).

**Outputs.** A continuous function `h(p)` with unlimited level of detail; domains of **970–3 386
km²** with 399–2 106 km of network (Table 3, 143:9); "several hundreds of square kilometers in a few
seconds" (143:8).

**Declared limitations** (143:9, "Limitations"):

- "it can generate only terrains that were once subject to hydraulic erosion";
- the network subdivides only upstream, so **no deltas and no oxbow lakes**; lakes appear nowhere
  else in the paper (the word occurs once, in this sentence);
- the greedy construction "can lead to meandering crests and outcroppings of unnatural-looking
  hills";
- **"the generated river network cannot easily adapt to large mountains with clearly articulated
  valleys. If the slope maps have strong variations flat rivers can cut through high mountains"**;
- many parameters, which allow unnatural terrains.

**Laws borrowed**: Horton–Strahler ordering (143:4), Rosgen river types (143:6), a flow–area power
law from Dunne & Leopold 1978 (143:6). **No valley-width law** (the word "width" does not occur)
and **no hillslope law** beyond the crest bound `λ ≤ 0.25`.

### What transfers to Ymir, and what does not

- **Inverse problem.** Génevaux grows a network and builds the relief to fit it. Ymir HAS a relief
  (tectonics + FBM) and wants valleys that fit it. The paper's fourth limitation is exactly Ymir's
  case — a network that must adapt to "large mountains with clearly articulated valleys" — and the
  paper says it does not do that well.
- **Scale.** Ymir's domain is 160 000 km², **47 to 165 times** the paper's largest and smallest
  examples.
- **What transfers**: the two-function idea (a river profile along a skeleton, a cross-profile
  perpendicular to it), and the crest as the place where two slopes from two rivers meet at a
  bounded gradient.

## 2. The laws, anchored

| law | form | value | source, page | domain | status |
| --- | --- | --- | --- | --- | --- |
| **valley floor profile** | `z(x) = z(x_b) + (E / K·A0^m)^{1/n}·χ`, `χ = ∫ (A0/A)^{m/n} dx` | `m/n = 0.5` (the model's) | Harel 2016 eqs. (4)–(7), PDF p. 15; method of Perron & Royden 2013 (ESPL 38, 570–576), cited PDF p. 14 | channel reaches that obey the stream-power law at uniform `E` | **anchored**. ADR Finding 96: valid as a PROFILE along one channel, invalid as a FIELD |
| **age parameter** `k = (E/K·A0^m)^{1/n}` | one number | Finding 96 measured `k_s = 0.0646` on the delivered field (Flint intercept, A0 = 1 km²) | ADR Finding 96, block B1 | the delivered world | **a calibration, not a literature value**: the round's "a relief target, not a duration" |
| **valley floor width** | `W_v = K_v·A^{c_v}` | `c_v = 0.3 ± 0.06` (Appalachian Plateau) | Clubb 2022 eq. (1) p. 438; abstract p. 437 | bedrock valleys, tectonically quiescent, homogeneous lithology, A 7–24 020 km² (Table 2, p. 452) | **exponent anchored** |
| | | large rivers `c_v = 0.22–0.37` | Clubb 2022 Table 2, p. 452 | Cumberland, Kentucky, Licking, Guyandotte, Little Kanawha, A 3 365–24 020 km² | anchored |
| | | model value `c_v = 0.46`, `K_v = 0.16` (undercutting–slump mechanism, soft bedrock) | Clubb 2022 p. 438–439 | numerical model | anchored as a bracket |
| | coefficient `K_v` | Table 2: 0.01–5.68 (tributaries), 0.08–2.71 (large rivers), "units m km^−2cv, drainage area in km²" | Clubb 2022 Table 2 caption, p. 452 | | ⛔ **NOT anchored — the units are inconsistent with the paper's own widths.** With `A` in km², Crane Creek (A 15.06 km², `c_v` 0.24, `K_v` 2.09) gives **4 m**, and the Cumberland River (24 020 km², 0.37, 0.08) gives **3 m**, while the widths the paper measures on its validation valleys are **344 ± 281 m** and **500–700 m** (pp. 443, 445; other valleys than Table 2's). With `A` in m² Table 2 gives 110 m and 552 m. I read the fits as made with `A` in m², **a reading the text does not state** |
| **hillslope gradient** | planar wall at a bounded slope | "average colluvial slope" **0.54 ± 0.11, 0.36 ± 0.05, 0.63 ± 0.26** m/m (28°, 20°, 32°) | Whipple & Tucker 1999 Table 1, p. 17,663 | King Range CA (high and low uplift), Central Range Taiwan: active orogens | **anchored** for an active orogen's colluvial slopes |
| **angle of repose 30–35°** | | `relief_v3` ships `talus_slope = tan(33°)` (ADR L14344) | — | — | ⛔ **no page in the dossier.** The anchoring report (`erosion_literature_anchoring.md`) does not cover it; the round calls it anchored, the dossier does not |
| **crest bound (Génevaux)** | `q_z = max(river nodes) + λ·d` | `λ ≤ 0.25` (14°) | Génevaux 2013 p. 143:6 | the paper's own generator | a design value, not a measurement |

### Which hillslope, and why

The round offers two: a wall at the angle of repose, or a mature diffusion profile (convex at the
top, straight, concave at the foot). **For the prototype I choose the planar wall at a bounded
slope, and I state the reason as a measurement to make rather than a belief.** A convex–concave
profile is what diffusion PRODUCES; laying it by hand is laying a second model's output. The
prototype's light incision pass and the shipped diffusion then act on the wall, and the σ and R8
columns say whether that is enough. The slope value is the one with a page: Whipple & Tucker's
colluvial range, **28°** (the King Range high-uplift value and the middle of the three).

⚠️ **Finding 96's warning carries over directly.** A field built from a **D8 path** stamps D8's
eight directions into the relief. The valley cross-section is therefore laid by the distance to a
**smoothed** trunk polyline, not to the D8 cells.
