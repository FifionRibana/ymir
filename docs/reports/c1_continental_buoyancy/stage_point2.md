# Issue #145 — Point 2 result: rigid continent WITH closures ON (Track D)

First run of the rigid continent with the **full tectonics** (all 7 closures incl Track D). RUN config: `run_with_closures`, `C1Closures::default()`, rigidity ON vs OFF, seeds 42 + 2, 64², 300 steps. Acceptance is **three-part** (scalar + spatial + visual) — never the boundary scalar alone.

## Measurements

| seed | rigid | craton area | land perim/area | n_comp | largest | bound edge S̃ | mass Δ |
|---|---|---|---|---|---|---|---|
| 42 | OFF | 0/561 (0%) | 1.33 | 23 | 0.67 | 0.20 | −20.3% |
| 42 | **ON** | **461/561 (82%)** | **0.58** | 17 | 0.63 | **0.21** | **+25.7%** |
| 2 | OFF | 116/714 (16%) | 0.76 | 27 | 0.69 | 0.73 | −28.9% |
| 2 | **ON** | **625/714 (88%)** | **0.53** | 43 | 0.93 | **0.79** | **+16.8%** |

## (a) Boundary scalar GATE — PASS

The closures-OFF prototype's edge pile-up (0.60→1.99) **does not occur with subduction on**: rigidON edge S̃ ≈ rigidOFF baseline (seed 42: 0.21 vs 0.20; seed 2: 0.79 vs 0.73). **Subduction drains the oceanic inflow at the rigid margin.** Stage S's S3 hypothesis confirmed in production — the caveat self-resolves.

## (b) Spatial bar — PASS (one soft spot)

Craton area preserved 82% / 88% (vs 0% / 16% OFF); land perim/area compact 0.58 / 0.53 (vs 1.33 / 0.76 OFF — OFF "land" is convergence pile-up, not continents); largest-component 0.63 / 0.93 (dominant mass). Soft spot: seed 2 `n_components = 43` — many tiny specks alongside the dominant mass (P95-cap threshold on a noisy field), not a structural failure.

## (c) Visual — PASS

`point2/seed{42,02}_{rigidON,rigidOFF}_land.png` (binary land/sea; also `_altitude` P95-cap render + `_s` thickness). rigidON (both seeds): credible continents — one dominant compact landmass, coast reasonably sharp, **no dentelures / holes / false-coast artefacts** from the rigidity+subduction interaction. rigidOFF: fragmented; its "land" is pile-up artefacts. Clear ON-vs-OFF improvement.

## Decisive signal for Point 6 (closure behaviour)

**Mass: rigidON +25.7% / +16.8% vs rigidOFF −20.3% / −28.9%.** Closures-OFF rigid was mass-exact (−0.00%); with closures ON the rigid continent shows a large **mass GAIN**. The closures **net-produce** on correct transport — they were calibrated to compensate the broken transport's destruction (≈−20%), so on rigid transport that compensation becomes over-production (≈+20%). This is the **Point-6 "calibration is a transport artefact" hypothesis confirmed in direction** → expect the **heavier closure re-foundation** branch of the conditional roadmap, not light re-validation. To quantify per-closure at Point 6.

## Verdict

Point 2 **PASSES all three volets**; S3 (subduction drains the boundary) confirmed in production. The integration is sound. The mass-gain signal pre-stages Point 6.
