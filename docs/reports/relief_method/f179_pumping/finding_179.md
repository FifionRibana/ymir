# Finding 179 — the Davis-Suppe pumping in the cascade's source (a diagnostic), and the Andean seed placed at the arc

**At the head of the report, as the rules ask:**
- **The sources do not contradict the reviewer's reading.** Dahlen 1990 (p. 59–60) and Willett, Slingerland & Hovius 2001
  (p. 456, 459) both make the orogenic influx kinematic (hV; U = V_c h / L) and let the wedge adjust to the erosion
  (W = hV / ē). Willett adds that the balance is not pointwise.
- **The cascade's relief accepted since F166 does not rest mostly on pumping, by this bound.**
  - On the collision belts the pumping factor is ~3 (Σ Davis-Suppe ≈ 3 × C1's net change of s), but the kinematic bound would
    cut their median s_src by only 4 %. P3 is false.
  - The pumping that matters is at the continent-ocean seeds: the bound cuts the active coasts' s_src by 13 % under F176 and by
    34–36 % under F178–F179.

**Status: Part A measured (nothing bounded, nothing fed back); Part B built behind the continent profile's flag (one closure
changed: where the Andean wedge starts, at the arc, the seed skipped like a boundary seed) and measured on C1; no adoption in
this finding; nothing committed without the feu vert.** Production unchanged (`f175_c1_hash` passes; production does not use the
sources mode, below).

- Declared before any measurement: `f179_declared.md` (the bound's formula from Dahlen and Willett with pages, the mapping to
  C1, d_arc). Predictions: `f179_predictions.md`.
- The bench: `crates/ymir-core/tests/f179_cascade.rs`, with an instrumented observer (`observe179`) that replays the loop's
  seeds and records, per Davis-Suppe cell and step, the family and the kinematic supply of the nearest seed. Results in
  `results/`, assembled in `f179_bench_output.txt`. C1 only, the 15 seeds of F176's selection × F175, F176, F178, F179.

## Part 0

- **`f165`**: nothing was removed or renamed. `tests/f165/main.rs` (added by `ee841a8`, F165) is a submodule of
  `f165_cascade.rs`, but cargo also builds every `tests/<dir>/main.rs` as a test target of its own, where its `super::` paths
  fail. **Repaired, test-only**: renamed `tests/f165/f165_main.rs`, with the `#[path]` updated. Every `ymir-core` test now
  compiles.
- **Production reads C1's own s.**
  - `run_with_closures` excludes no term (`tectonics_c1/cached_product.rs:281`);
  - `upscale_from_c1_with_progress(&state, …)` (`:305`) → `c1_coarse_normalized_altitude(state, …)`
    (`production_upscale.rs:455`) → `c1_coarse_raw_altitude(state, …)` (`:293`).
  - The sources mode (`cascade/history.rs`, `run_c1_sources*`) serves only the cascade. **The pumping does not reach
    production.**

## Part A — the pumping (measured, nothing bounded)

### A1 — the pumping factor (Σ Davis-Suppe / C1's net change of s, cells with a positive net, pooled over the 15 seeds)

| family | F175 | F176 | F178 | F179 |
|---|---|---|---|---|
| collision | 2.92 (p90 20.0) | **3.55** (p90 24.2) | 3.28 (20.1) | 3.06 (22.1) |
| O-C boundary | 5.07 (39.5) | 6.34 (41.7) | 6.44 (43.4) | 5.25 (39.2) |
| Andean | — | — | **9.71** (84.8) | 9.54 (103.8) |

The maxima run to 10³–10⁵ (cells whose net change is near 0). The cells with a net ≤ 0 (unbounded factor) are counted in
the per-seed lines.

### A2 — the kinematic supply (Dahlen p. 59–60, Willett p. 459): Δs_kin = (v_rel · n̂) · dt · 64 · s_in / L, L = 30 cells

| family | share of the steps above the supply (F175 / F176 / F178 / F179) | excess per family cell, F176 | F178 | F179 |
|---|---|---|---|---|
| collision | 38 / 35 / 36 / 38 % | 0.087 s (3.0 km of crust) | 0.077 s (2.7 km) | 0.082 s (2.9 km) |
| O-C boundary | 44 / **54** / 38 / 32 % | 0.685 s (24.0 km) | 0.694 s (24.3 km) | 0.899 s (31.5 km) |
| Andean | — / — / **49** / 46 % | — | **2.79 s (97.7 km)** | 1.89 s (66.0 km) |

### A3 — the bounded s_src (an estimate of what F180 would change; nothing fed back)

| | F175 | F176 | F178 | F179 |
|---|---|---|---|---|
| collision belts, pooled median s_src → bounded | 1.04 → 1.00 (−4 %) | **1.11 → 1.06 (−4 %)** | 1.20 → 1.15 (−4 %) | 1.19 → 1.14 (−4 %) |
| active 0–3 band, median of the seeds' medians | 0.82 → 0.72 | 0.87 → 0.76 (−13 %) | **4.07 → 2.70 (−34 %)** | 3.11 → 1.98 (−36 %) |

- **The reading**: on the collision belts Davis-Suppe gives ~3 × the net change, but its increments mostly stay under the
  kinematic supply. A continental incoming s (~1) gives a generous supply, so the bound removes little.
- At the continent-ocean seeds (an oceanic incoming s of ~0.2, the boosted O-C rate) the increments exceed the supply on half
  the steps, and the excess is tens of km of crust per cell. That is where the cascade's runaway sits (F178), and where F180's
  bound would act.
- The supply is generous by construction: h is the incoming cell's whole s, an upper bound of the accreted thickness. So these
  are lower estimates of the excess.

## Part B — the Andean seed at the arc (d_arc = 4, D_gap = 10; skipped like a boundary seed)

| | F176 | F178 (front) | **F179 (arc)** |
|---|---|---|---|
| land at the end of C1 | 50.0–57.5 % (median 51.1) | 54.8–70.4 % (65.3) | 53.3–65.7 % (61.5) |
| circumnavigation lost | 0 | 4 (23, 3, 4, 10) | **1 (seed 10)** |
| r* = 0 (still circumnavigable) | 0 | 6 | 5 |
| land added over F176, median | — | 522 cells (99 % in the band, 82 % continental at t₀) | 339 cells (98 % in the band, 96 % continental at t₀) |
| A1, any / majority | 27.4 / 17.9 % | 25.0 / 18.3 % | 24.9 / 18.3 % |
| active 0–6 cells: Δh_src median (range) | 2.4 km (0.4–4.7) | 18.3 km (12.3–22.2) | **12.0 km (9.0–14.6)** |
| of which the C-O wedge | 0.4 km | 18.7 km | 11.8 km |
| passive 0–6 cells | 1.6 km | 9.4 km | 4.9 km |
| active ≥ 2 × passive (0–6) | 4 / 15 | 4 / 15 | 10 / 15 |
| Andean pumping factor (median) | — | 9.71 | 9.54 |
| Andean seed instances on cells emerged at t₀ | — | **0 %** | **0 %** |

- **The seeds' position** (`images/seeds_F178_vs_F179_panel.png`):
  - F179's seeds sit at walk depth 4 (222 411 of 249 451 instances; 3 when another trench is nearer). They are a cell or two
    inland of F178's, but still on the band's submerged continental shelf: **0 % on cells emerged at t₀**.
  - The seaward margin puts 32 km of crust (s = 0.914) on Part A's edge, a hair under the emergence (s ≥ 0.915). A depth-4 seed
    reaches the shelf, not the emerged continent.
- **So placing the seed at the arc reduces F178's runaway but does not stop it.**
  - The active coasts still rise by 9–15 km in the cascade's source (the pumping, Part A);
  - the band's shelf still emerges: +339 cells median, 98 % in the band;
  - the circumnavigation is lost on 1 seed, and r* is 0 on 5 more.
  - The bound (F180) is what Part A says acts there.

## Predictions

**The reviewer's**:

| | prediction | outcome |
|---|---|---|
| P1 | collision belts' median factor > 2 | held (3.55 under F176) |
| P2 | bound exceeded on more than half the steps of F176's O-C and F178's Andean seeds | refuted narrowly (O-C 54 %, Andean 49 %) |
| P3 | the bound cuts the collision belts' median s_src by ≥ 30 % | refuted (−4 %) |
| P4 | the circumnavigation lost on ≤ 1 seed with the arc seed | held (1) |
| P5 | ≥ 80 % of the arc seeds on cells emerged at t₀ | refuted (0 %) |
| P6 | the arc seed does not halve the Andean pumping | held (factor 9.54 against 9.71; excess −32 %) |

**Mine**:
- held:
  - Part 0;
  - P1 (2–5);
  - P5 refuted;
  - P6;
  - the A1 shares within ± 3 points.
- refuted:
  - P2 « ≥ 80 % » (32–54 %);
  - P3 « 20–50 % » (4 %);
  - P4 « lost on 1–3 » held at the low end (1);
  - the land « 52–62 % » (53–66 %).

## Open, for the author (nothing chosen; F180 will build the bound)

- (a) The bound's terms for F180:
  - h = the incoming cell's whole s (generous);
  - L = Davis-Suppe's reach (30 cells, wide: the collision family covers most continental interiors,
    `images/pumping_panel_F176_F178_F179.png`). Dahlen's W is the wedge's actual width; a narrower L gives a tighter bound.
- (b) Under the seaward margin even the arc's depth falls on the submerged shelf: 0 % of seeds on emerged cells, which the
  brief chose not to use as a criterion.
- (c) The `f165` repair is test-only, and enters with this finding if the author agrees.

**At commit (2026-10-18), decided**: the author's words: « Feu vert pour commiter le F179 : drapeau coupé, réparation test-only de f165 incluse (renommée en f165_main.rs). »
- **F179 is committed with its flag off**: `Margin::andean_arc = None` by default, so production, F176, F177 and F178 are
  unchanged bit for bit (`the_arc_seed_controls`; production's C1 hash). `Margin::seaward_oo_collision_arc()` is only used on
  request.
- Committed: the bench (`f179_cascade.rs`), the arc seed (`andean_seeds_arc`, `C1Ext::andean_arc`, the Davis-Suppe skip mask
  `apply_davis_suppe_step_routed_skip`), the control, and the test-only `f165` repair (`tests/f165/main.rs` → `f165_main.rs`).
  `config.rs` and `workspace.rs` stay out.
- **The reviewer's errors, recorded (the author's)**:
  1. **F178**: D_gap was applied as the radius of a walk leaving the trench in every direction. The seeds then reach coasts that
     do not face the trench, and the Andean family nearly encircles the continents in C1 (the pumping images, rows F178 and
     F179). The rampart comes back in C1 unseen, since W₁₀ at 256² was not measured.
  2. **F179, P5**: the reviewer predicted that a seed at the arc's depth would fall on emerged cells. He had not seen that the
     B′ ramp puts the continent's edge at s = 0.914, just under the emergence (0.915).
- **Next**: F180 builds the kinematic bound (the author's brief).
