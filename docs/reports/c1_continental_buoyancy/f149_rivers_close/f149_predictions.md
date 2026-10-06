# Finding 149 — predictions, written 2026-10-06 BEFORE any measurement of this round

**Reading declaration, unfavourable.** Non-blind on Findings 73–148 and the reviewer's predictions, in particular:
- **F148's head-score numbers**:
  - AUC 0.784 at σ = 2 and 0.857 at σ = 1, the head being the first 10 trace cells;
  - the head medians (a) −0.12 · (b) −0.07 · (c) +1.26 km⁻¹;
  - (c)'s max-score quantiles (p0.5 = −1.24).
- So **N = 10 is chosen KNOWING it gave 0.784**. σ = 1 is non-blind, as the brief says.
- **F148-K**: lake 1000011's spillway 16217 carries `source_lake = 1000011` and starts 2 cells off the shore. In the
  témoin it is the only Unresolved lake.
- **F148-J**: 10 246 rivers. 873 fusion candidates and 705 micro-rivers into a lake < 1 km outside the control; 726
  micro-lake rivers in all; (a) ∩ (b) = 0.

## Predictions

- **P-H (the head score, σ = 2, N = 10)**:
  - **the stop rule fires again**;
  - at the θ losing ≤ 1 % of (c), (a) ∪ (b) loses **5–20 %** (the reviewer: 15–30 %);
  - the head AUC is the same population statistic as F148's (0.784), but the rule now judges a river once, at its
    head, so the control's low tail sets θ;
  - σ = 1 (non-blind) does better but **also fires** (≤ 35 % at ≤ 1 %).
- **P-H, the control rivers lost at that θ**:
  - mostly SHORT rivers (< 3 km) lying ≥ 50 % on a construction floor, that is the control's floor criterion, not its
    ≥ 10 km one;
  - **fewer than a third are pair members** (`parallel_of`), against the reviewer's ≥ 1/3 parallel traces on the
    built floors;
  - most are tributaries whose head sits on the construction's planar floor or wall.
- **P-K**:
  - on the témoin, **only lake 1000011** changes type (Unresolved → Exorheic), with the reviewer;
  - in the other five guard states, **0–2 lakes each** change (the same rim-cell mechanism), so their lake hashes
    change too, listed;
  - the field guard is untouched (6/6 field = banc).
- **P-T**: the toggle hides about **1 650 rivers**:
  - ~940 fusion candidates + 726 micro-lake rivers, with an overlap of a few tens;
  - about **3 800 km** (the candidates' ~3 450 km plus 726 short micro-rivers of ~0.5 km); with the reviewer's
    ~1 600 in count;
  - a cost of ~0.
- **P-cost**: no world-time change (the flag is a filter; the K fix is a set lookup).

## Meta

My disagreements with the reviewer: H's range and the composition of the lost control rivers. **At least one of my
predictions is wrong.**
