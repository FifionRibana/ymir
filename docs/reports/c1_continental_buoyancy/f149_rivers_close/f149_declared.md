# F149 — instruments and methods, declared BEFORE the measurements (2026-10-06)

**The world**: the témoin (C2 /10 col), `viz_hd_lakes_on`, as F146–F148.

**Before / after**:
- the guard 6/6 (`ymir-viz f123_viz_guard`, field and lakes through `run_hd`);
- C2 /10 col field `a8d2d538d692c2f0`;
- lib and viz tests; `cargo check --workspace`.
- "Before" = F148's after-checks (`d937aec`, the same tree).

**Rule 15, recorded**: F148's stop rule is not relaxed after the fact. Part H uses it unchanged.

## T — the masking toggle (display only)

- **`micro_lake_inflow`** = `end` is a lake or an endorheic lake, AND `length_km < 1.0`. This is F148-M (b)'s
  definition, without the control exclusion: a flag describes the river, not a population.
- Format 0.5.0.
- **The viz**: « Masquer faisceaux et micro-rivières » hides `fusion_candidate || micro_lake_inflow`. A separate box
  hides `retraces_spillway` (F147's). The counts and the km hidden are shown.
- **The gates**:
  - `rivers.json` byte-identical (the existing test);
  - the field and lakes guards unchanged by T (T touches only the export and the viz);
  - the bench reports the hidden count and km on the témoin.

## K — the core fix

- `exorheic_lakes_missing_outlet`: a lake HAS an outlet if a segment carries `segment_source_lake = Some(id)`, in
  addition to the existing 1-cell border test. No wider distance.
- **The type changes are attributed in one run.** The only change is in that function, so the lakes whose type
  changes are exactly the Exorheic lakes whose only outlet is a tagged segment with no border source.
  - `f134_lake_guard` (the writer of `bench_lake_hashes.json`, run as declared) reports them per state:
    `id · area · the tagged segment's distance to the footprint`.
  - **The regeneration is declared.** The diff of `bench_lake_hashes.json` must touch only the states where that list
    is non-empty, and only their `lakes_hash`. Each lake count is unchanged.
- **Gates**:
  - on the témoin, only lake 1000011 changes (any other is listed);
  - a permanent test: a tagged spillway 2 cells off the shore keeps its lake Exorheic. The negative controls are the
    same segment untagged, and a lake with no segment at all, which both go Unresolved;
  - the F86 tests stay green;
  - the viz guard 6/6 after the regeneration.

## H — one attempt: the head score (the stop rule unchanged)

- **The score**: the mean of the convergence C (F148's index, σ = 2 cells) over the **first N = 10 cells** of a
  candidate head's trace.
  - The reason: F148's head score (the same N) is where the separation was measured, AUC 0.784. N = 10 is ~0.5 km,
    half of the micro-rivers' 1 km bound.
  - **NOT blind**: F148 printed that AUC.
- **σ = 1 is reported beside it, labelled NOT blind** (seen in F148). The decision uses σ = 2.
- **The rule simulated** (`HeadRule` gains a mode; bench only), per segment upstream before downstream:
  - a reach leaving a lake's shore, or a spillway: channel from its start (F148);
  - a segment with a channel reach joining it at index j: channel from j. Its own head part (above j, when j > 0) is
    channel from 0 if its head score passes;
  - a segment with no channel input: channel from 0 if its head score (its first min(N, len) cells) is ≥ θ, else not
    a channel. **No sliding**: a head is judged once.
  - Then the downstream closure as in F148.
- **The populations** (a), (b), (c) are F148's definitions on the same build. **θ** = 0, and the quantiles 0.5 / 1 /
  2 / 5 / 10 / 20 % of (c)'s head score. The table is F148's.
- **The stop rule, unchanged**: a θ must remove ≥ 50 % of (a) ∪ (b) while losing ≤ 1 % of (c).
  - **If it fires**: nothing is built. The bundles and micro-rivers become an accepted defect, masked by T.
  - **If it does not**: the rule is built gated in `rivers_ll` only, with F148-P's gates.
- **The lost control rivers**, at the θ nearest the 1 % bound (the largest θ losing ≤ 1 %, else the smallest
  quantile):
  - length;
  - the floor share;
  - pair members (`parallel_of`) or not;
  - their head's position: on a construction floor cell / a carved wall cell (`carved && !floor`) / elsewhere;
  - their end type.

## The cost

T and K's seconds per world (expected ~0).

## Amendment after the first after-check (2026-10-06, NOT blind)

**What happened**: the viz guard read "lakes != banc" for exactly the four regenerated states, while the field read
"= banc" everywhere (`checks_after_run1_STALE_CACHE.txt`).

**The cause, my error**: K changes the HD drainage bundle's code, and `cache.rs` requires `ALGO_HD_DRAINAGE` to be
bumped on any such change. It was not. So `run_hd` served the pre-fix bundle from the drainage cache (lake 1000011
still Unresolved), while the bench recomputed.

**The amendment**:
- `ALGO_HD_DRAINAGE` 11 -> 12;
- the bump changes the HD drainage key, which is the lake guard's `digest`, so the reference is regenerated a second
  time (`f149_lake_guard_raw_regen2.txt`);
- **expected diff against the first regeneration: the six `digest`s only, every `lakes_hash` identical** (kept:
  `bench_lake_hashes_REGEN1.json`);
- then the viz guard is run again.

