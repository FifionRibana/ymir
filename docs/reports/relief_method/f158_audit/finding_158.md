# Finding 158 — changing the relief method: the current chain audited from the code, the literature gathered, and a multi-scale cascade (64² → 8192²) sketched with three costed paths and a first prototype for F159 (no production code)

**Status: an audit and a bibliography. No production code changes; nothing measured on a world. Nothing committed.**
**ESTIMATE** marks every cost not measured.

**Process:**
- My predictions were written before reading the code for this round (`f158_predictions.md`, the non-blind items
  declared).
- The references are in `docs/refs` (PDFs, never committed). The texts were read from the PDFs themselves;
  paraphrased, with short quotes only.

## The author's criterion, recorded as given

- On the lake base's validation crops (2026-10-08): upstream « oui », downstream « non », non-regression « pas
  spécialement », the added lakes « oui globalement ».
- « En sortie de lac (et à certains endroits le long d'une rivière), on voit un promontoire en demi-cercle. La rivière
  le contourne. C'est très étrange. » « Les crêtes sont horribles ! Il y a aussi des murs systématiquement de
  quasiment 900 m entre le haut et le bas de la chaîne de montagne. » « Les incisions dans le relief sont toujours
  comme un pinceau qu'on a utilisé pour creuser. Ça ne fait en aucun cas naturel. »
- « À mon avis, la méthode utilisée a résolu des problèmes mais donne un relief qui est trop artificiel. Je pense
  qu'il faut qu'on fasse de la biblio et qu'on trouve quelque chose de mieux et plus adapté. Avant, on avait
  l'avantage d'avoir un relief plus naturel, mais on galérait avec les rivières et les lacs, et le vieillissement. De
  plus, on n'a pas la sédimentation, qu'il faudrait qu'on déduise à un moment. Il doit bien exister quelque chose qui
  peut caractériser le relief d'un continent vieux avec ses lacs et rivières. »
- **The author's approach (2026-10-08)**: « Peut-on avoir une approche incrémentale ? Typiquement, les reliefs
  répondent à des structures à différentes échelles. On pourrait les construire petit à petit, raffinant une carte
  d'abord à 64², puis 128², puis 256², jusqu'au 8192². Chaque niveau produirait la forme ou le relief répondant à la
  résolution, plutôt que de ne faire que la tectonique à 64² puis passer directement à 8192². Typiquement le principe
  du LOD. »
- **The objective**: a terrain without artefacts; otherwise the method or the closure is not the right one.
- **The cost principle**: a world in minutes (reference 249.8 s).

## Partie 0

- **F157 committed** (`a8723db`), its ADR entry completed at commit:
  - the author's verdict; the lake base is **not promoted** and stays gated;
  - **the reviewer's correction**: F145's closing of the cones (« the normal geometry of a dissected range ») was
    false. `accepted_defects.md` § 1 is **reopened** (the original entry kept below it, struck), and the line in ADR
    F145's entry marked "[CORRECTED in F158]";
  - the attribution of the author's three defects, **hypotheses to check**: the promontory = the cone of a
    polyline's last sample; the crests = the intersection of the 28° walls; the ~900 m walls = flat-floored valleys
    carved down to the law; the "brush" = polyline × width × cones.
- **No production code changes this round.**
- **Checks before** = F157's after-checks (`a8723db`; no `src` change since).
- **Checks after**: `cargo check` clean; lib 621; viz 32; **the guard 6 / 6, field and lakes = banc** (C2 /10 col
  `a8d2d538d692c2f0`). The geology stage logged 5–6 s, against ~2 s before, with the machine loaded (the code is
  unchanged).

## A — the audit of the current chain (the code, lines cited)

The chain, in `tectonics_c1/production_upscale.rs`, `upscale_from_c1_with_progress` (432):
1. the coarse altitude at 64² (455);
2. `upscale_with_fbm` to 8192² (459);
3. the C-2 edifices (472);
4. [the témoin: the valley construction (525)];
5. the stream-power block (580; the témoin: one light pass, 625);
6. droplets (706, **off** in production);
7. the active rims (736);
8. the bathymetry (759);
9. then `run_hd`'s breach, drainage, climate, biomes, geology and export.

| # | question | what the code does | lines |
|---|---|---|---|
| 1 | **stages and resolutions** | **The tectonics at 64² (C1), then ONE jump to 8192².** `c1_coarse_normalized_altitude` (64²: isostasy + Stein-Stein, normalised) → `upscale_with_fbm` straight to the target. **No intermediate computation resolution.** The other grids are side products: the geology's 1024² smoothing grid (F152), the K field upscaled to HD (`upscale_k_to_hd`), the coarse land report (64²). | `production_upscale.rs:287, 455, 459` |
| 2 | **the jump's interpolation** | **Bilinear** (`sample_bilinear_periodic`) of the 64² field, sampled at a coast-warped position (coast warp 1.5 coarse cells). | `upscale.rs:662–668` |
| 2b | **the detail between 64² and 8192²** | **Noise only**: an anisotropic FBM (7 octaves, oriented along the coarse slope, max anisotropy 3, domain warp 0 in production), its frequency fixed in coarse-cell coordinates (#151). Its amplitude is capped by C-1's relief budget (β = 0.1 × the coarse slope), and tapered at the coast. **The processes act only at 8192²**: the construction and the stream power. | `upscale.rs:29–125, 455–490, 536–700` |
| 3 | **the incision** | **IMPLICIT, Braun & Willett (2013) "FastScape"**: a topological stack of the D8 receivers, updated from base level upward; closed form for n = 1, « UNCONDITIONALLY STABLE (no CFL timestep limit) ». E = K·A_km²^0.5·S. Delivered: relief-v3, K = 1500 × 3 = 4500, dt 1, **2 iterations** (k_time 9 000). Témoin: **one light pass** at 0.1 k_time. At 8192² only. F115's Courant 3 699 is therefore stable but not integrating ("stability is not accuracy", 1013–1016). **Cost**: one 8192² pass ≈ 30 s (F155: removing the light pass saves 31.5 s per world). | `stream_power.rs:1–15, 352, 410–420, 584–720, 1011–1016`; `production_upscale.rs:620–626` |
| 4 | **routing** | **D8** receivers and stack (`compute_flow` on the FILLED surface, Garbrecht–Martz flats); **MFD (p = 2) for the incision's drainage area only**. The skeleton: D8 on the breached field (D8-LTD optional, off). | `stream_power.rs:597–604`; `flow.rs:190–330, 878` |
| 5 | **depressions during erosion** | **Priority-flood FILLING** before each pass: the receivers come from the filled surface. **A1** (`depression_floor`): a cell inside a filled depression is not incised. **No basin-graph routing.** After the chain, `run_hd`'s breach (priority-flood breach plus a fill mop-up). | `flow.rs:201–222, 1018`; `stream_power.rs:703–710` |
| 6 | **deposition** | **None.** A detachment-limited law; the base-level bound « can stop erosion and can never deposit » (49). The only raises are the talus's mass-conserving transfer and the construction's minimal plain (gated). Droplet deposition exists but is **off** (`cfg.erosion = None`). F156's "light-pass deposition" on the design is the talus and the diffusion moving matter, not sediment. | `stream_power.rs:23–60, 837–870`; `upscale.rs:410` |
| 7 | **hillslope diffusion** | **Yes, linear and explicit**: D = 0.08 (relief-v3), 4 sub-steps, applied everywhere (`diffuse_channels`). **Plus a talus sweep** at tan 33° (4 passes, factor 0.5, mass-conserving), and lateral bank planing. The nonlinear Roering-type closure exists (relief-v2: S_c = tan 33°, implicit Picard) but is off in v3. | `stream_power.rs:163–168, 236–290, 410–436, 800–870` |
| 8 | **uplift and time** | **No uplift during the erosion.** The isostatic altitude is the initial condition, and time enters only as k_time = K·dt·iterations, a dimensionless budget. The doc: « the duration is NOT observable ». | `stream_power.rs:131–160, 1051–1053` |
| 9 | **why the construction (F121) was added** | Against the delivered world's defects (ADR F121 table, OFF delivered): coast spurs **+115** (F76–F83); canyon / over-dug class (F102: 29.6 % of the depressions); R8 terrain **0.0924**, stripes (F121); lakes 54 with D_L > 5 at 35 %; transverse thalweg offset **−5.06 m**, entrenchment (F113); relief p50 424 m against the oracle's 488 m. **The oracle that fixed most of it cost 300 passes, 6 210 s** (F95); the χ floor « is a profile, not a field » (F96). The construction lays the χ law along the pre-incision trunks (F120–F121) instead. | ADR F95, F96, F102, F120, F121 |
| 10 | **resolution dependence** | **The law is physical**: A in km², S in m/m, so K does not depend on the cell. The diffusion is scaled by (195.3 m / cell)². **But**: the channels are 1 cell wide; the D8 / MFD areas and the A_c = 0.1 km² threshold sit on the grid; the base-level bound is 0.5 m; k_time is a budget, not a duration. **The ×7.5 geographic scale ratio acts on the hydrology only** (areas × 56.25, discharge, widths, navigability), never on the terrain's K or slopes. | `stream_power.rs:486–493`; `drainage.rs:1281–1310` |

**Read**:
- **Ymir already holds a FastScape-type implicit solver, MFD, diffusion and a talus**, all at 8192² only, run 1–2
  times, from a bilinear + FBM jump of ×128.
- **What it lacks**:
  - any intermediate level;
  - any uplift and time;
  - any sediment transport or deposition;
  - any basin-graph routing of the depressions.

## B — the bibliography (the files in `docs/refs`)

### Multi-scale and amplification (the author's track, first)

**Schott, Galin, Guérin, Peytavie & Paris 2024**, « Terrain Amplification using Multi-scale Erosion », ACM TOG /
SIGGRAPH (`2024-MultiScaleHydro-Author.pdf`).
- **What**: amplifies a low-resolution terrain into a high-resolution, « hydrologically consistent » one. The resolution
  is raised step by step (A_k = D_k ∘ T_k ∘ E_k ∘ U_k):
  - **U**: a **bicubic** upsampling, ×2 or ×4 per step, that adds no detail;
  - **E**: an iterative stream power, **without uplift**, with clamps s_max and a_max on the slope and the area (« more
    uniformly distributed erosion features »), and K(p) = k(1 − ρ(p)) with a hardness map ρ, often a fractal noise
    to break self-similarity and the « axis-aligned artifacts »;
  - **T**: thermal stabilisation to a talus angle, whose critical slope is **perturbed by noise** (a constant one
    « yields uniform slopes »);
  - **D**: a separate sediment creation, transport and deposition step, not mass-conserving by design.
  - At the end, a diffusion-based **re-targeting** restores the coarse elevations of peaks, crests and saddles. A
    **multi-scale breaching** (radii halving) guarantees drainage without one-cell canyons.
- **Routing**: an iterative parallel MFD-like drainage (exponent p = 1.3, « to avoid sharp fluvial incision »), one
  routing iteration per erosion step.
- **Scales**: inputs of 128² or 256².
  - Table 1 lists the amplification steps up to **4096²**, with 300–7 400 erosion iterations per level and 300–800 at
    4096².
  - **8192²** is shown in Fig. 1 (256² → 8192², ×32), and the paper states it as its experiments' maximum, « on a
    single up-to-date GPU ».
  - **[CORRECTED in F159, precision added: the paper gives no total time at 8192², only per-iteration timings.]**
- **Cost**: GPU, per iteration 19.2 ms (erosion), 3.5 (thermal), 41.5 (deposition) at 8192². The re-targeting and the
  breaching take seconds on a CPU at 4096².
- **What it brings to Ymir**: the author's cascade, exactly. The coarse terrain is kept (bounded carving + re-
  targeting). Details come from processes at each scale, not from noise. Drainage is guaranteed without one-cell
  trenches.
- **Its limits**:
  - a design tool, not a physical model: no uplift, no time, no mass conservation;
  - it is controlled by per-level maximum carvings;
  - its realism relies on the input's large forms;
  - and on a GPU (Ymir is CPU-only).

**Schott, Paris, Fournier, Guérin & Galin 2023**, « Large-scale terrain authoring through interactive erosion
simulation », ACM TOG (`schott2023.pdf`; MIT code, `github.com/H-Schott/StreamPowerErosion`).
- **What**: the terrain emerges from an **uplift** field under the stream power, with « a fast yet accurate
  approximation of drainage area » parallel on a GPU. It also does inverse modelling: the uplift reconstructed from a
  DEM.
- **For Ymir**: Ymir's C1 tectonics could supply the uplift. The relief would then be a steady state of uplift and
  erosion, which gives the « vieux » (old) relief its time.
- **Limits**: a single resolution; uniform fractal patterns (« exaggerated self similarities across scales », as the
  2024 paper says of it).

**Cordonnier, Braun, Cani, Benes, Galin, Peytavie & Guérin 2016**, « Large Scale Terrain Generation from Tectonic
Uplift and Fluvial Erosion », CGF (`cordonnier2016.pdf`).
- **What**: a stream graph under the stream power with uplift, on a coarse vector representation. **Then the DEM is
  made « by blending landform feature kernels »** whose parameters come from the graph.
- **For Ymir**: its first half (uplift + stream power on a graph) is relevant. **Its second half is the family Ymir's
  valley construction belongs to**: a profile laid along lines and blended. That is what the author sees as the
  "brush".

### The process pieces

**Braun & Willett 2013**, Geomorphology 180–181, 170–179 (DOI 10.1016/j.geomorph.2012.10.008).
- **NOT in `docs/refs`**: no open copy was found (no HAL or author copy; the publisher's page only). This sheet is
  written from the papers that restate it (Cordonnier 2019, Yuan 2019, Schott 2024) and from Ymir's own
  implementation.
- **What**: the stream power E = K A^m S^n solved implicitly over a stack of the receivers, O(n), unconditionally
  stable for n = 1. Accuracy in transients depends on the step (its Appendix B, quoted by Landlab), but it is exact at
  steady state.
- **For Ymir**: Ymir already has it (`stream_power.rs`). **Its use is what differs**: Ymir runs 1–2 huge steps
  without uplift. FastScape runs many steps with uplift towards a steady state.

**Cordonnier, Bovy & Braun 2019**, « A versatile, linear complexity algorithm for flow routing in topographies with
depressions », ESurf 7, 549–562 (`cordonnier2019.pdf`).
- **What**: instead of a priority-flood fill, it builds **a graph connecting the adjacent drainage basins** and
  computes the flow paths within and across the depressions explicitly, in linear time. The user chooses filling or
  carving inside the depressions. The graph « has the potential to be reused … such as the simulation of erosion ».
- **For Ymir**: the lakes as depressions that the routing crosses without filling the terrain. It fits Ymir's lakes
  (levels, spills) better than fill + breach, and it is the basis for depositing sediment into lakes.
- **Limits**: sequential; its gain over priority flood is mostly speed inside an LEM.

**Yuan, Braun, Guerit, Rouby & Cordonnier 2019**, « A new efficient method to solve the stream power law model taking
into account sediment deposition », JGR Earth Surface 124(6) (`yuan2019.pdf`, HAL author version).
- **What**: an **implicit O(N)** solution of Davy & Lague's (2009) erosion–deposition model. The deposition depends on
  the upstream sediment flux, and the solution is « unconditionally stable even when large time steps are used ». It
  covers the transition from detachment-limited to transport-limited behaviour. It is shown on an uplifting region
  beside a foreland basin: the foreland aggrades, with « autogenic aggradation and incision cycles ».
- **For Ymir**: **the sedimentation the author asks for**, inside the same implicit stack Ymir already has: plains,
  forelands, the filling of basins.
- **Limits**: lakes and depressions need the routing above to receive sediment; a G (deposition) coefficient to
  calibrate.

**Perron, Kirchner & Dietrich 2009**, « Formation of evenly spaced ridges and valleys », Nature 460, 502–505
(`perron2009.pdf`).
- **What**: a Péclet-like number compares advective stream incision with diffusive creep. Pe = 1 gives a length
  L_c, and **the valley spacing is proportional to L_c**. It is checked against five landscapes.
- **For Ymir**: the scale of valleys is set by D / K, not by the grid. **A level whose cell is larger than L_c cannot
  carry its first-order valleys**: they belong to the finer levels. This is the physical basis of "each level carries
  the forms its resolution allows".

**Salles, Mallard & Zahirovic 2020**, « gospl: Global Scalable Paleo Landscape Evolution », JOSS (`salles2020.pdf`).
- **What**: a parallel global LEM (unstructured meshes, variable resolution). Its pieces: « inland river deposition in
  depressions computed using priority-flood techniques »; marine deposition by diffusion at river mouths; hillslope
  processes on land and at sea.
- **For Ymir**: a reference implementation of the full set (fluvial, depressions filled by sediment, marine, hillslope)
  at continental scale.
- **Limits**: Python / PETSc, made for geological times and coarse meshes.

### What these papers cite on the three asked topics

- **Post-orogenic decay (an old continent's relief)**:
  - **Baldwin, Whipple & Tucker 2003**, JGR 108, 2158. Abstract only, not in `docs/refs`: a detachment-limited stream
    power decays the relief in ~1–10 Myr. Isostasy, transport-limited behaviour and an erosion threshold all
    lengthen it; a threshold ×~20.
  - **Braun & Robert 2005**, ESPL 30, 1203 (`braun2005.pdf`): in the Dabie Shan, the mean relief fell by « a factor
    of 2·5 to 4·5 during the last 60–80 Ma », at 0.01–0.04 km/Myr, with flexural isostatic rebound.
  - ⇒ **An old continent is a decaying relief with isostatic rebound and thresholds**, not a fresh steady state. Ymir
    has the isostasy (C1) but no time.
- **Lakes filled by sediment**: goSPL's depression deposition (above), and Yuan 2019's deposition. Cordonnier 2019's
  basin graph is the routing they rely on.
- **Grid artefacts and resolution dependence of the stream power**:
  - **Kwang & Parker 2017** (`kwang2017.pdf`): with m/n = 0.5 and no hillslope diffusion, the steady landscape is
    « invariant to horizontal stretching ». A 10 km² domain can be stretched into the 1 000 km² one. ⇒ **diffusion
    (a length scale) is needed for scales to differ.**
  - **Armitage 2019**, ESurf 7, 67 (`armitage2019.pdf`): routing node to node with multiple directions « significantly
    reduces » the resolution dependence. « LEMs need to capture processes at a sub-grid-scale. »
  - **Perron et al. 2008** (cited by Armitage): a sub-grid flow width removes the resolution dependence of the valley
    spacing, at the price of a response time that depends on that width. Not fetched.
  - **Pelletier 2010** (`pelletier2010.pdf`, already present; its text did not extract): minimising the grid dependence
    of flow routing.
  - Schott 2023 / 2024: axis-aligned artefacts reduced by noisy hardness and by a softer MFD exponent (p = 1.3).
- **Coarse-to-fine and multigrid for LEMs**: **no geoscience paper found** that builds a landscape coarse-to-fine with
  an LEM at each level. Schott 2024 is the only multi-scale erosion method found. The LEM literature handles
  resolution through sub-grid processes (Armitage, Perron) and variable meshes (goSPL, FastScapeLib on irregular
  meshes). **This is a gap**: the cascade would be new as a physical model.

## S — the synthesis for the author

### S1 — a sketch of the cascade (domain 400 km; the cell = 400 km / n)

| level | cell | the forms it can carry | the processes that run there | what it hands on, what is frozen | cost (ESTIMATE, CPU) |
|---|---|---|---|---|---|
| **64²** | 6.25 km | plates, cratons, belts, the continent's outline | C1 tectonics, isostasy (as now) | the altitude, the uplift field (new: U from the tectonics), the masks | as now (≈ 0.3 s) |
| **128²** | 3.1 km | mountain ranges, major basins, the main divides | **uplift + implicit stream power + diffusion to a steady state or a decay age** (FastScape-type, many steps); flexural isostasy | **the major drainage tree, the big lakes (basin graph), the sediment budget** | < 1 s |
| **256²** | 1.56 km | the large rivers' valleys, forelands, big lakes | the same, starting from the upsampled level; **deposition** (Yuan) into forelands and lakes | the trunks (A ≥ ~1 000 km²) frozen as routing | 1–2 s |
| **512²** | 780 m | tributary valleys, plains, lake shorelines | stream power + deposition + diffusion; lake routing (basin graph) | trunks ≥ ~100 km²; lake levels | 2–5 s |
| **1024²** | 390 m | second-order valleys, terraces, alluvial fans | the same; the talus with a noisy critical slope | valleys ≥ ~10 km² | 5–15 s |
| **2048²** | 195 m | first-order valleys (the scale of L_c), crests, scree | stream power (bounded, à la Schott) + diffusion + talus | the river network as exported | 15–40 s |
| **4096²** | 98 m | gullies, crest detail | a few bounded erosion and talus iterations; no new tree | — | 20–60 s |
| **8192²** | 49 m | the finest texture | a few bounded iterations, or noise only; a multi-scale breach; climate, geology, export as now | the export | 30–90 s |

- **Total ESTIMATE: 2–4 min per world** on a CPU, against 249.8 s today.
- **The basis**: Ymir's own light pass costs ≈ 30 s per 8192² pass (F155), so ~8 s at 4096², ~2 s at 2048², < 0.5 s at
  1024² (cost ∝ cells). The step counts per level follow Schott 2024's (thousands at 256², a few hundred above 2048²);
  the CPU costs come from this round's ratios. **Not measured.**

### S2 — the pieces

| piece | Ymir has it | missing | the reference | reusable as is |
|---|---|---|---|---|
| tectonics, isostasy (64²) | yes (C1) | an **uplift rate** U(x, t) rather than an altitude | Schott 2023, Cordonnier 2016 (uplift-driven) | **yes** |
| implicit stream power | yes (`stream_power.rs`) | uplift; time; many steps per level | Braun & Willett 2013 | **yes** (the solver) |
| routing D8 + MFD | yes | MFD in the stack; the node-to-node distributive routing | Armitage 2019; Schott 2024 (p = 1.3) | partly |
| depressions | fill + breach | **the basin graph** (lakes crossed, not filled) | Cordonnier 2019 | no (to build) |
| deposition, sedimentation | none (droplets off) | **erosion–deposition** (plains, forelands, lake infill) | Yuan 2019; goSPL | no (to build) |
| hillslope | linear D + talus 33° | a noisy critical slope; the Péclet / L_c scale per level | Perron 2009; Schott 2024 | yes, with changes |
| intermediate levels | none | **upsampling + the per-level loop** | Schott 2024 | no (to build) |
| re-targeting, multi-scale breach | the breach (single scale) | the multi-scale breach | Schott 2024 | partly |
| climate, biomes | yes | — | — | **yes** (they read the final field) |
| geology (substratum / surface) | yes | — (the rock could feed K's hardness ρ) | Schott 2024 (hardness) | **yes** |
| rivers export (rivers_ll), lakes, guard, `.ymir` | yes | — (they read the final drainage) | — | **yes**; the guard to regenerate |
| valley construction (F121), gorge, cones | yes | — | — | **dropped** by the cascade |

### S3 — the hard points, named

1. **The resolution dependence.** The law is physical (A in km², S in m/m), but:
   - the channels are a cell wide;
   - the D8 / MFD areas depend on the grid;
   - without diffusion the stream power is scale-free (Kwang & Parker).
   - **Each level needs a length scale (D / K, L_c) and a channel-head area that are physical, not "1 cell".**
     Armitage and Perron 2008 show sub-grid routing or width is the known fix.
2. **Keeping the coarse level.**
   - A finer level must not move the coarse forms.
   - Schott 2024 bounds the carving per level and re-targets the peaks and saddles.
   - A physical cascade could instead pass **the coarse level's uplift and base levels** (lakes, trunks) as
     constraints and let the fine level reach its own steady state.
   - **Which constraint keeps the coarse level without freezing it into a template is the central open question.**
3. **The interpolation.**
   - Bilinear ×128 today (a fan of planes); Schott uses bicubic ×2–×4 per level.
   - Each upsample creates new pits and flats that the next level's routing must resolve: the basin graph, not a fill.
4. **The share of noise.**
   - Today all the detail between 6.25 km and 49 m is FBM.
   - In the cascade, noise only perturbs the processes: the hardness / K (Schott), the critical slope, the uplift.
   - **Possibly no additive noise at all below 1 024².**
   - The geology's rock grid is a natural, non-random hardness field.

### S4 — three paths, costed, none chosen

| | (a) the multi-scale cascade (the author's track) | (b) one full LEM at a coarse resolution, then the current upscale | (c) the delivered chain completed, without the construction |
|---|---|---|---|
| **what it is** | 64² → 128² → … → 8192², an LEM (uplift, implicit stream power, deposition, diffusion, basin-graph lakes) at each level, bicubic ×2 between | an LEM with uplift and deposition at 512² or 1024² to a decay age, then the bilinear + FBM jump to 8192² as now | the delivered relief-v3 (implicit stream power at 8192²) plus uplift / time, the basin graph and deposition; the construction removed |
| **the artefacts** (promontories, wall crests, 900 m walls, brush) | **removed at the source**: no construction, no cones; valleys emerge from the processes at their scale | removed down to the LEM's cell (~400–800 m); below it, FBM (the "natural but plain" look) | removed (no construction); the old defects come back (stripes R8 0.092, entrenchment, coast spurs), to fix in the processes |
| **lakes** | the basin graph at every level; lakes arise and fill with sediment | at the coarse level; then the upscale makes pits again (fill / breach) | the basin graph at 8192² |
| **ageing** | uplift and time per level (decay ages, Baldwin / Braun & Robert) | the same, at the coarse level only | needs many steps at 8192² (F95: 300 passes = 6 210 s) unless the uplift / steady-state form converges fast |
| **sedimentation** | yes (Yuan / goSPL), at each level | yes, at the coarse level | yes, at 8192² (costly) |
| **what it breaks** | the construction, the gorge, the cones, the χ laid law, F133's lake base (replaced); the guard (all states); the exports' content (not their format) | the same minus the multi-level loop; the upscale stays | the construction and its dependent work |
| **cost per world (ESTIMATE)** | **2–4 min** (S1) | **1–3 min** (an LEM at 1024² for ~10³ steps ≈ 0.5 s each ≈ 1–2 min, plus the current upscale) | **5–60 min** (many 8192² steps at ≈ 30 s each: 10 steps = 5 min, 100 = 50 min) |
| **its own risk** | new as a physical model (S3-2); the most work | the FBM below 1 km stays "natural but plain"; the jump ×8–16 is still a jump | the cost; F95's oracle already showed it |

### S5 — a first prototype for F159 (proposed, not built)

**"Does the cascade hold?"** — the smallest test: **64² → 128² → 256² → 512² on the témoin's seed**, no export, a
bench only.
- **The pieces**:
  - the current C1 output at 64² as the initial altitude;
  - an uplift field derived from C1 (a declared PROXY: the isostatic altitude's excess as U);
  - Ymir's implicit stream power with uplift and **many steps** per level, to a declared age;
  - linear diffusion with a declared D (L_c ≈ 2–4 cells at 512²);
  - bicubic ×2 between levels;
  - priority-flood filling for now (the basin graph and deposition come later, separately).
- **The measures that would decide**, at 512², against the current témoin downsampled to 512²:
  1. **the artefacts**: the planar-wall share (28° ± 0.5°), sharp crests, comb teeth (F157-B7's instruments; the
     construction gives 17 % planar walls);
  2. **the stability across levels**: the 128² divides and trunks still at the same place at 512² (the share of
     trunks within 2 cells; F120-A0's instrument);
  3. **the forms per level**: the valley spacing against L_c (Perron), the R8 terrain (axis alignment), θ against the
     declared CI;
  4. **the resolution dependence**: the same cascade stopped at 256² and at 512², compared on the shared scales
     (hypsometry, relief p50, drainage density at A ≥ 10 km²);
  5. **the cost per level**.
  6. **Images for the author**: the hillshade (with the lake outline) at 512² of the cascade, the current témoin, and
     the old delivered (relief-v3 without the construction), on 3 declared crops.
- **The stop rule, to write before measuring**: if the 128² trunks move by more than 2 cells at 512², or the
  planar-wall share exceeds the delivered's, the cascade as built does not hold, and F160 says why.

## Predictions

**The reviewer's (hypotheses to check; judged by the reading):**
- "the current chain has no intermediate computation resolution between the tectonic grid and the HD grid (at most
  diagnostic grids)": **held** (one bilinear jump 64² → 8192²; the 1024² geology grid and the K field are side
  products).
- "the delivered incision is explicit, with neither deposition nor hillslope diffusion": **refuted**.
  - It is **implicit** (Braun & Willett's stack, « unconditionally stable »).
  - It has **linear hillslope diffusion** (D = 0.08) and a talus.
  - "No deposition" **held**.
- "its depressions by filling or breaching, not a basin graph": **held**.
- "tectonics, isostasy, climate, geology and the rivers' export reusable as they are": **held**. The climate, the
  geology and rivers_ll read the final field and its drainage.
  - The tectonics would gain an uplift output.
  - The guard must be regenerated.
- Meta: **held**.

**Mine (`f158_predictions.md`):**
- **A**:
  - P-A1 (no intermediate level; bilinear; the 1024² geology grid): **held**;
  - P-A2 (noise + processes at 8192² only): **held**;
  - P-A3 (implicit; 2 passes delivered, one light pass; tens of seconds): **held** (≈ 30 s per pass);
  - P-A4 (D8 + MFD for the area; D8 skeleton): **held**;
  - P-A5 (filling, A1, final breach; no basin graph): **held**;
  - P-A6 (no deposition; the raises are fills): **held in part**. A1 does not fill: it **skips** the incision in the
    depressions. The raises are the talus's transfer and the gated plain;
  - **P-A7 (diffusion off or negligible in relief-v3): refuted.** D = 0.08 explicit, plus the talus, on everywhere;
  - P-A8 (no uplift; k_time a budget): **held**;
  - P-A9 (the construction's reasons): **held** (canyons, stripes, lakes, the oracle's cost, χ), plus the coast spurs
    and the entrenchment;
  - P-A10 (×7.5 on the hydrology only): **held**.
- **B**:
  - P-B1 (Schott 2024 keeps the coarse elevations, GPU seconds at 4K–8K): **held**; the "coarse kept" is a bound plus
    a re-targeting, not a hard constraint;
  - P-B2 (resolution dependence through the area and the slope; fixed by K scaling or χ forms): **refuted in its
    fix**. The literature fixes it by a hillslope length scale (Kwang & Parker; Perron's L_c) and sub-grid routing or
    width (Armitage; Perron 2008), not by scaling K.
- **Meta**: **held**.

## Limitations, stated

1. **No world was measured.** Every cost in S is an ESTIMATE from this round's ratios.
2. **Braun & Willett 2013 and Baldwin et al. 2003 are not in `docs/refs`** (no open copy found). Their sheets come
   from the papers restating them and from abstracts. Perron 2008 was not fetched.
3. **The coarse-to-fine LEM** has no geoscience precedent found here. Absence in this search is not proof of
   absence.
4. **The attribution of the author's three defects** (Partie 0) is not measured. It is labelled as hypotheses.

## State

**Uncommitted**, awaiting the feu vert.
- ADR Finding 158. (The F157 completion, `accepted_defects.md` § 1's reopening and the F145 correction mark were
  committed with F157, `a8723db`.)
- This folder.
- `docs/refs`: new PDFs, not committed: `cordonnier2019`, `yuan2019`, `salles2020`, `schott2023`, `perron2009`,
  `braun2005`, `armitage2019`.
- No code.
