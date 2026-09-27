# Finding 123 (continued) — the identity guard, the teeth under rule 18, and the trunk identity

Draft for the ADR. Seed 1, 8192². Follows `finding_123_A.md`.

## 1 · The identity guard (`tectonics_c1/bench_guard.rs`)

**What it compares.** The viz's eroded field (pre-breach) against the bench's field AT THE SAME
SETTINGS, the settings being the eroded product's cache digest (`eroded_key_full`, the function the
viz already uses for its cache). Reference: `crates/ymir-core/data/bench_field_hashes.json`, written by
`f123_guard`, embedded at compile time. Verdicts: `Match` / `Mismatch` (the viz then shows "≠ banc"
in place of every microscope number) / `NoReference` ("non gardé"). A badge states the verdict at the
top of the microscope.

**Controls, all measured.**
- **Inertness since Finding 109**: the delivered field at the benches' framing hashes to
  `316ca8f727ca6f62`, Finding 109's recorded value (ADR L17300), bit for bit. Every seam added since
  (a1_exempt, valley_construction, basin_base, the benches' config refactor) is inert by default.
- **End to end through `run_hd`** (`f123_viz_guard`, the viz's own pipeline with the workspace's
  HdParams for seed 1 at 8192²): **`Match` for "livré" and for "C2 /10 col"**, digests
  `f61ef23f4f4e9d94` and `ba7a1a85a85c6b92`.

⛔ **A framing is NOT a translation.** The viz auto-frames seed 1 at torus offset x = 0.09375; the
benches of Findings 107–123 use 0.0. The viz-framed delivered field against the benches' field rolled
by 768 cells: **22.3 % of cells bit-identical, max |Δ| 382.9 m.** The world the author has been
looking at is NOT the world those Findings measured. The reference therefore records both framings,
labelled "cadrage viz" and "cadrage des Findings 107-123"; the Findings' numbers apply only to the
second. To read them in the viz, pan the framing to x offset 0.

**The "S2, 580 km², 654 km" entry is still not found**, at either framing (the viz-framed C2/10's
closest S2 entries: 459 km² / 300 km, 798 km² / 299 km). The framing does not explain it; the run's
resolution and the entry's mouth are needed.

**"C2 /10 col" in the panel.** The valley states are a combo box now (a five-button row could clip to
three in a narrow panel). They sit in EXPERT mode, beside the age box.

## 2 · Rule 18 on the teeth (`f123_rule18`)

The 1 787 C2/10 teeth paired with the A1+B2 segment whose own mouth lies within 3 cells and whose
own area is within a factor 1.5:

| population | sinuosity p50 | p90 | R8 at chord 8 cells |
| --- | --- | --- | --- |
| C2/10 teeth, paired (82) | **1.000** | 1.068 | **0.684** |
| their A1+B2 pairs (82) | 1.044 | 1.149 | 0.375 |
| every A1+B2 segment | 1.055 | 1.137 | 0.475 |

**Only 82 of 1 787 teeth (4.6 %) have a counterpart in A1+B2.** The rest do not exist there: the
walls CREATE the teeth, at places where the incised world has no drawn segment. Where a counterpart
exists, it is sinuous (1.044) and the tooth is dead straight (1.000), with an 8-fold anisotropy nearly
double. By the round's reading, **the teeth are the plane, not the tectonic D8** — and B2 (hierarchy,
which dissects the walls) is the remedy it names, B1 (profile + noise) the one to compare it with.

## 3 · The trunk identity (C2), applied

The microscope's trunk is now the path of **maximum drained area** from the mouth to the source, at
every confluence and at every exorheic-lake crossing (`TRUNK_BY_MAX_AREA` in `workspace.rs`), with the
same rule in the benches' port.

**Two corrections the control forced.** Plain max-A lost the lake crossing of 3 delivered trunks:
each was an EXACT tie of inherited areas (5 906 / 5 906, 46 093 / 46 093, 24 604 / 24 604 km² — a
clipped fragment carries its parent's area, Findings 42/93), which `max_by` broke by position. Ties
to the longest path left 2: near-ties of 0.05 % and 0.004 % (5 742 / 5 745, 44 645 / 44 647 km²).
**Two areas within 1 % (PROXY) now tie**, and the tie goes to the longest-path rule.

| world | rule | one-reach | length p50 / p90 / max (km, signified) | lake crossings Σ | trunks that stop crossing a lake |
| --- | --- | --- | --- | --- | --- |
| delivered | longest path (was) | 69.7 % | 10.3 / 103.3 / 272.5 | 7 | — |
| delivered | **max-A, 1 % tie** | 69.7 % | 10.3 / 94.8 / 313.8 | **8** | **0** |
| C2/10 | longest path (was) | 58.7 % | 20.5 / 137.0 / 1 522.7 | 0 | — |
| C2/10 | **max-A, 1 % tie** | 58.7 % | 20.1 / 132.9 / 1 172.6 | 0 | **0** |

Finding 93's control does not regress. The viz builds and its 29 unit tests pass.

## Waiting

- **D** — the author's eye on "C2 /10 col" (now guarded: `= banc`). B and E wait on it.
- **The "654 km" entry** — the run's resolution and the mouth.
