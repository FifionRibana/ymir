# #155 relief — critical-wedge / high-mountain ceiling diagnostic

Target: high mountains (Himalaya/Tibet). Current orogens render ~2265-2596 m (the
metric scale, `probe_sea_unification_acceptance`); the goal is the Tibetan-plateau
scale. This finishes terrestrial relief BEFORE climate (which depends on orography).
Diagnostic only — the mechanism follows the verdict.

## 1. EH exonerated (re-confirmed on current terrain)

`probe_orogen_equilibrium` (current code): margin orogen S̃ p95 = **1.30 / 1.47 / 1.35**
(seeds 42/1988/2026), erosion-OFF **1.59 / 1.73 / 1.59**. `probe_prominence_attribution`:
**0 %** of margin cells above `h_eq` (0/481, 0/909, 0/461). EH (`apply_equilibrium_height_step`)
collapses S̃ quadratically ONLY above `h_eq=2.0` — a no-op below 2.0, so it does
NOTHING to orogens stuck at 1.5. **Raising `h_eq` would have zero effect.** (The old
"EH-bound at 2.0" was the oceanic advective spike, not the continental orogen.)

## 2. The thickening mechanism + what limits it

Davis-Suppe source (`apply_davis_suppe_step`):
`∂S̃/∂t = coupling·|v|·max(0, h_crit(d) − S̃)·exp(−d/L_decay)`, defaults `coupling=2.0`,
`h_max=2.5`, `l_taper=4`, `l_decay=6`; O-C margin-peaked profile (1b-i) `h_crit =
h_max·exp(−d/l_taper)` peaks ~1.95 at the margin. The source deposits TOWARD `h_crit`
and halts as `S̃ → h_crit`. **Achieved S̃ (1.3-1.5) sits well below the target (1.95)**:
the equilibrium is where the source balances (a) the in-loop stream-power EROSION sink
(removes ~0.2-0.25: 1.5 baseline vs 1.7 erosion-OFF) and (b) ADVECTIVE smear (erosion-OFF
still 1.7 < 1.95 target) over the 300-step run. **The limiter is the deposition RATE
(coupling) vs erosion+advection removal — NOT the EH ceiling.**

## 3. Critical-wedge — the anchored target is S̃≈2.0 (= h_eq)

Himalayan crust ~70 km vs normal ~35 km = **2×** (TDD §11: S̃=1.0 = 35 km normal
continental). So **S̃=2.0 = 70 km = the Tibetan-plateau thickness — which is EXACTLY
`h_eq=2.0`.** EH=2.0 is therefore the CORRECT Himalayan-plateau cap; the orogens should
reach ~2.0 but the force only gets them to 1.5. **Critical-wedge fix = stronger
convergent-margin deposition to reach S̃≈2.0** (the EH cap = 70 km Himalayan), NOT
raising EH. The unused headroom 1.5→2.0 IS the missing thickening. Levers (anchored,
for the fix maillon): raise DS `coupling` (deposition rate) at active convergent margins
and/or reduce the margin erosion sink, targeting S̃≈2.0 (= 70 km, the existing h_eq). The
ratio is the craton precedent (craton 1.25×; orogen 2×), physically anchored, not a knob.

## 4. The codomain tension (decide BEFORE coding)

Even at S̃=2.0, the altitude is bounded:
- **Isostasy** `peak_altitude_m ≤ max_elevation_m = 4000 m` — a hard codomain ceiling.
  An S̃=2.0 (70 km) orogen renders ≤4000 m, but Airy of 70 km vs 35 km ≈ **5.4 km** — so
  `max_elevation_m=4000` UNDER-represents the 70 km thickness. Raise to ~5500 m (Airy of
  70 km).
- **Vertical contract** `norm→m = (norm−0.5)·11300`, norm 1.0 = **+5650 m**. Land
  currently caps at norm 0.854 (the 4000 m max_elev). So raising `max_elevation_m` to
  ~5500 uses the EXISTING contract headroom (5500 < 5650) — **no re-pin needed**.
- **8000 m (Everest peaks) > 5650 m contract max** → would require RE-PINNING the
  vertical contract (raise `depth_scale`/half-range), which is COUPLED to the S-S ocean
  scale (half-range = 5651/5000). TENSION.

**VERDICT on the target:** C1 is a 64²→2048² MEAN field. The right C1 target is the
**Tibetan-PLATEAU scale ~5-5.5 km** (Airy of 70 km crust, S̃=2.0) — reachable by
(force → S̃≈2.0) + (raise `max_elevation_m` to ~5500, within the contract, NO re-pin).
**8000 m Everest POINT peaks are sub-grid erosion-sharpening beyond C1's mean-field
resolution AND beyond the vertical contract — out of scope / defer.** Targeting the
plateau (~5.5 km) avoids re-pinning the vertical contract (and its ocean coupling).

## 5. Orographic effect (for the climate ordering)

The critical-wedge AMPLIFIES the existing convergent-margin orogens (raises their S̃ →
higher coastal cordilleras), it does not create new chains. So higher orographic
barriers at the existing margins → confirms doing the critical-wedge BEFORE climate is
the right order (the barriers will shape precipitation).

## Output / scope for the fix maillon
1. **Force:** raise DS convergent-margin deposition (coupling / margin-erosion balance)
   to lift orogen S̃ from ~1.5 toward **S̃≈2.0** (= 70 km Himalayan = the existing h_eq).
   Do NOT raise `h_eq` (2.0 is correct). Anchored on the 2× crustal ratio.
2. **Codomain:** raise `max_elevation_m` 4000 → ~5500 (Airy of 70 km crust), within the
   vertical contract (no re-pin).
3. **Target = Tibetan plateau ~5.5 km** (mean field), NOT 8000 m Everest peaks (sub-grid
   + beyond contract → defer). This resolves the codomain tension without a vertical re-pin.
4. EH untouched; the vertical/horizontal coordinate contract untouched.
