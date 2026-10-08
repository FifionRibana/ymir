# F158 — my predictions, written BEFORE reading the code for this round (2026-10-08)

**Not blind, declared.** In earlier rounds I had read:
- `production_upscale.rs`'s stream-power block: `incise_with_floor`, the light pass, A1+B2, the P-modes;
- its receiver array and topological stack (« receivers before donors, via the donor CSR tree »), and the MFD
  accumulation used « for the incision ONLY »;
- `production_hd_config` (`cfg.erosion = None`: relief-v3 has no droplet pass);
- F115's « Courant 3 699 »;
- `breach_monotone_protected` (a priority-flood breach plus a fill mop-up);
- `compute_flow`'s priority-flood fill;
- the valley construction (F120–F121);
- the bench `Knobs::passes(2)`.

## A — the audit

- **P-A1**: no intermediate COMPUTATION resolution between the 64² tectonics and the 8192² HD grid.
  - The jump is one interpolation (bicubic or bilinear) of the coarse normalised altitude, then FBM.
  - The only intermediate grids are a geology smoothing grid (1024²) and diagnostics.
- **P-A2**: the detail between 64² and 8192² is noise (FBM, anisotropic, domain warp), plus processes only at 8192²
  (the construction and the stream power).
- **P-A3**: **the incision is IMPLICIT** (a Braun–Willett receiver stack, O(n)). It is stable at any Courant, hence
  F115's 3 699.
  - It runs at 8192² only.
  - The delivered world takes 2 passes at k_time / 2; the témoin, one light pass at 0.1 k_time.
  - Tens of seconds per pass.
- **P-A4**: D8 for the receivers and the stack; MFD (p = 2) for the incision's drainage area only; D8 for the
  skeleton (D8-LTD optional, off).
- **P-A5**: the depressions are FILLED before each pass (priority flood), with A1's depression floor and the final
  breach. **No basin-graph routing** (Cordonnier et al. 2019).
- **P-A6**: **no sediment transport and no deposition** (detachment-limited). The only "deposition" is A1's filling
  of the depressions, which is not mass-conserving. F156's light-pass deposition is that fill.
- **P-A7**: a linear hillslope diffusion term exists in the code, but it is off or negligible in relief-v3.
- **P-A8**: no uplift during the erosion. The isostatic altitude is the initial condition, and time enters only as
  the K·t budget (k_time).
- **P-A9**: the construction (F121) was added against the delivered world's:
  - canyons (the class at 29.6 %, F102);
  - stripes and anisotropy (R8 0.092);
  - drowned / brush lakes (F122);
  - the oracle's cost (300 passes, 6 210 s, F95);
  - θ / the χ profile not matching the law (F96).
- **P-A10**: the ×7.5 geographic scale ratio acts on the hydrology only (discharge, widths, signified areas), never on
  the terrain's K or slopes. Areas are in km² by the cell size, and K·t is per pass.

## B — the literature

- **P-B1**: Schott et al. 2024 amplifies a coarse terrain with stream power and hillslope processes at several
  scales, keeping the coarse elevations as a constraint. It runs in seconds on a GPU at 4K–8K.
- **P-B2**: the stream power's resolution dependence comes through the drainage area and the slope's grid scale.
  The literature keeps the forms across levels by scaling K with the cell size, or by working in χ / steady-state
  forms.

## The reviewer's predictions (hypotheses to check)

- "no intermediate resolution of computation": I expect **held** (P-A1);
- "the incision is explicit, with neither deposition nor hillslope diffusion": I expect **refuted on "explicit"**
  (P-A3), and held on deposition. Diffusion: probably held in relief-v3 (P-A7);
- "the depressions by filling or breaching, not a basin graph": I expect **held**;
- "tectonics, isostasy, climate, geology and the rivers' export reusable as they are": I expect **held**, with the
  climate, the geology and rivers_ll reading the final field and its drainage (downstream, so unchanged);
- meta: held.

## Meta

- At least one of my predictions is refuted.
