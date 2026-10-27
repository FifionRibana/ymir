# F189 Part 2 -- the predictions (written before the computation)

**Part 2: the brief and the predictions** (written before any computation of Part 2):
- **The author's decision** (relayed by the reviewer): the corrections of `flow.rs` apply everywhere, production included;
  the guard will move.
- **The reviewer's predictions**:
  - **V2 bis redone** (V2's setup, the lakes on, no fluvial deposit; fluvial steps 160 / 40 / 10 kyr):
    - no land cell above 5 × the land's median, at any step;
    - on the quantities above E3's noise floor, r = change(40, 10) / change(160, 40) between 0.1 and 0.5, with the same sign.
  - **If V2 bis passes, E3's triplet 20 / 10 / 5 kyr redone** on the corrected adopted configuration:
    r = change(10, 5) / change(20, 10) between 0.25 and 1, with the same sign, on the quantities above the floor. Production
    is added against the 5 kyr step.
  - **The reading**:
    - if both pass, the step dependence came from the stranding. Production's change against the fine step is documented
      as an accepted limit, and the hillslope's splitting in sub-steps is not needed;
    - if V2 bis fails, stop: the hillslope's splitting is decided next.
- **CC's predictions**:
  - the controls:
    - the stranded volume (fluvial and marine) is 0 by construction;
    - no land cell without a receiver outside the lakes at any sub-step on seed 8;
    - the mass < 10⁻⁹ with the « exported » term;
  - **the « exported » volume** (the ocean filled with no deeper cell): below 1 % of the marine sources on every panel
    seed. A hypothesis: the open water of the panel's worlds stays deep enough at the intermediate age;
  - **the guard**: the six states' field and lake digests all move. `run_hd`'s routing goes through `compute_flow`; the
    per-cell sea makes (iii) impossible there, but the flats' tie (ii) touches every flat;
  - **V2 bis**: a hypothesis. G1 bis traced V2's spikes to the fluvial stranding (the talus shed 93 km in one step), so the
    spikes go and the 5 × criterion passes; r is not predicted beyond the reviewer's range;
  - **the step dependence (E3)**: a hypothesis, no number. Most of the fluvial stranding sat on two cells (71 %), so part of
    E3's dependence may go with it.
