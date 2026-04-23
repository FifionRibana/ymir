# Ymir — Tectonic Solver Scaling & Nondimensionalization

**Status**: analysis note, not yet implemented
**Scope**: Phase 1 (tectonics) — Stokes solve, viscosity closure, advection of $S$, boundary source terms, basal drag, mantle forcing, geological age field
**Out of scope**: the six-phase pipeline downstream of tectonics (isostasy, upscale+FBM, hydraulic erosion, rivers/lakes, climate) — each will get its own note.

---

## 1. Problem statement

The tectonic solver as it stands (MAC-staggered Stokes, matrix-free operator, JFNK or Picard with BiCGSTAB/CG inner, adaptive sub-stepping, plastic yielding, cratonic rigidity, basal drag, boundary sources, mantle convection proxy) is **functionally capable of producing continents that look broadly plausible**, but three pain points keep resurfacing.

**P1 — Time cost at working grid sizes.** Development iterations at $64^2$ and $128^2$ should be near-interactive (sub-second per tectonic step at $64^2$, seconds per step at $128^2$). Production at $512^2$ is accepted to be slow (tens of seconds to minutes per step) but must terminate reliably. Today, Newton often exhausts its iteration budget and sub-stepping shrinks $\Delta t$ aggressively, which multiplies the cost beyond what the grid size alone justifies.

**P2 — Nonlinear solver robustness.** JFNK classification distinguishes `Stagnation`, `Oscillation`, `Divergence`, `ConvergedOnState`, `ConvergedOnResidual`, `MaxIterations`. The taxonomy is itself a sound design decision — early bail-out with a smaller sub-step is cheaper than plowing through Newton iterations that won't converge — but the *frequency* of non-`ConvergedOnResidual` outcomes indicates that the problem we hand to Newton is harder than it needs to be.

**P3 — Physical incoherence between effects.** The solver currently exposes several independent knobs — `subduction_rate = 0.5`, `volcanic_arc_rate = 0.15`, `spreading_rate = 0.3`, `slab_pull_factor = 0.05`, `basal_friction = 0.05`, `cratonic.max_factor = 5.0`, `yielding.yield_stress = 50.0`, `eta_min = 1e-3`, `eta_max = 1e4`, `gravity_factor = 1.0`, `max_plate_velocity = 5.0`, … — each tuned in isolation. None of these numbers is expressed against a common physical reference, so their *ratios* carry no intuitive meaning and there is no way to answer questions like "is slab pull currently stronger or weaker than GPE spreading at a 70 km-thick collision front?" without instrumenting the code. Incoherent parameter cocktails feed back into P1–P2 because they produce viscosity contrasts $\eta_\text{max}/\eta_\text{min}$ much larger than what the target physics actually requires, which wrecks the conditioning of the Stokes operator.

The three pain points are linked: a well-posed, well-scaled problem is cheaper to solve, more robustly, with solver parameters that can be set once from physical reasoning rather than tuned empirically per run.

**Target outcomes of this analysis.**

1. A consistent set of physical scales $(L^\*, S^\*, \tau^\*, v^\*, \eta^\*, \rho^\*)$ against which every equation and every user-facing parameter is expressed.
2. A rewritten set of equations where every surviving parameter is a **dimensionless number** with a known admissible range.
3. An inventory of those dimensionless numbers with their physical meaning and their **qualitative effect on continental character** — not on specific mountain types, see §3.3.
4. A validation strategy based on **stochastic sampling over seeds** rather than hand-picked benchmarks.
5. An implementation roadmap that introduces a `Scales` struct, migrates the solver to internal non-dimensional variables, adds a geological age field, and preserves backward compatibility at the I/O boundary.
6. An explicit list of open questions — choices we are *not* making in this note.

---

## 2. Root cause analysis

### 2.1 Numerical conditioning of the Stokes operator

On the MAC grid with periodic BCs, the Stokes operator $A$ for the momentum balance has, at each velocity DOF, a diagonal contribution scaling as

$$
A_\text{diag}(i,j) \;\sim\; \frac{\eta(i,j)}{(\Delta x)^2} \;+\; C_b \cdot S(i,j) \cdot \frac{1}{(\Delta x)^2}.
$$

The condition number of $A$ for a Stokes problem with variable viscosity on a periodic grid scales, to leading order, as

$$
\kappa(A) \;\sim\; \underbrace{\frac{\eta_\text{max}}{\eta_\text{min}}}_{\text{viscosity contrast}} \;\times\; \underbrace{\left(\frac{L}{\Delta x}\right)^2}_{\text{grid resolution}}.
$$

With the current defaults $\eta_\text{min} = 10^{-3}$, $\eta_\text{max} = 10^4$, the viscosity ratio is $10^7$. At $128^2$, $L/\Delta x = 128$, so $\kappa(A) \sim 10^7 \cdot 10^4 = 10^{11}$. At $512^2$, $\kappa(A) \sim 10^{12}$. BiCGSTAB iteration count grows roughly as $\sqrt{\kappa}$, so we pay a $\sim\!10^{5.5}$ iteration penalty *from conditioning alone*, before any physical nonlinearity enters.

The Jacobi preconditioner (currently with a `1e-20` diagonal floor to preserve near-singular information) only normalizes the diagonal; it does not help with the off-diagonal coupling that dominates at high viscosity contrast. SSOR is marginally better but still a point preconditioner. Multigrid or block preconditioners would be the textbook fix, but they come with significant implementation cost — **nondimensionalization attacks the same problem much more cheaply first, because the physically necessary viscosity contrast for this application is much smaller than $10^7$** (see §4.2 and §5). Multigrid remains the strategic direction for $256^2$–$512^2$ performance (see §6 Phase F).

### 2.2 Parametric incoherence across physical effects

Every physical mechanism currently has its own dimensional scale, implicitly fixed by the choice of numerical value in the config. Let us inventory the state of play.

**GPE driving force.** `compute_rhs` emits $-\nabla(\rho(1 - \rho/\rho_m) S^2)$. The `gravity_factor = 1.0` is a dimensionless scalar prefactor, but the density-corrected GPE prefactor $\rho(1 - \rho/\rho_m)$ is in $\text{kg/m}^3$, so the RHS carries hidden SI units that the viscous LHS $\eta \nabla^2 v$ does not match — $[\eta]$ is unspecified and $[\Delta x] = 1/n$ in grid units. **The RHS and LHS do not live in the same unit system.** They happen to produce plausible numbers because magic numbers were tuned to match, but no scaling argument protects the solver against a change in $n$, $\rho_m$, or $\eta_\text{min}$.

**Viscosity closure.** $\eta = (\dot\varepsilon_{II} + \varepsilon_\text{min})^{1/n - 1}$ with `strain_rate_min = 1e-3`, `eta_max = 1e4`, `eta_min = 1e-3`. The strain-rate floor $\varepsilon_\text{min}$ has hidden units $[1/\tau]$ where $\tau$ is the dimensional time of the simulation, itself not explicitly defined — $\Delta t$ is picked by a CFL condition on $\Delta x / v_\text{max}$ in grid units. So $\varepsilon_\text{min}$ means different things at different grid sizes, and the viscosity field it feeds into is implicitly resolution-dependent.

**Source terms.** `subduction_rate = 0.5`, `volcanic_arc_rate = 0.15`, `spreading_rate = 0.3`, `collision_volcanism_rate = 0.05`, etc., are multiplied by `convergence_rate = (velocity difference) / dx`. The product has units of $[v]/[L] = [1/\tau]$, which is added to $\partial_t S$. Fine so far — but the *ratios* between these rates are arbitrary in the current config. In reality these ratios are physically constrained (arc volcanism is always a fraction of slab mass consumed; spreading at ridges must balance cooling so oceanic crust stabilises near 7 km) but those constraints are not encoded anywhere.

**Slab pull.** `slab_pull_factor = 0.05` multiplied by `plate.subducted_mass` gives an extra plate traction. `subducted_mass` grows monotonically (no units documented). Over long simulations it pushes plate velocities to `max_plate_velocity = 5.0` (hard clamp). The clamp is a safeguard, but it is also a symptom that slab-pull growth is not inherently bounded by a physical time scale.

**Basal drag.** `basal_friction = 0.05`, coupled as $C_b S v / \Delta x^2$. The `1/dx²` normalization (added in commit `3e26890`) is a *resolution-independence fix*, not a physical scaling. Physically, basal drag scales with $\eta_\text{asthenosphere} / h_\text{asthenosphere}$ — a real ratio with real units. Here it is an opaque magic number.

**Cratonic rigidity.** `cratonic.max_factor = 5.0` multiplies $\eta$ near plate centers. Physically cratons are $10^2$–$10^3 \times$ stiffer than mobile belts. A factor of 5 is on the low end, which is fine as a *choice* but the motivation must be explicit. Empirical observation: raising it to 20+ causes Newton to stall or fail to converge at current scaling — a direct symptom of the $\kappa(A)$ budget being already saturated. §4.10 proposes decoupling the "crust resists deformation" behaviour from the viscosity contrast.

**Yielding.** `yield_stress = 50` in solver units. Solver units for stress $\tau = 2\eta\dot\varepsilon$ are $[\eta][1/\tau]$, which given the implicit unit system means the yield stress is neither in Pa nor in anything physically labelled. A Drucker-Prager brittle-ductile transition in the real lithosphere occurs around 100–500 MPa; we have no way to know whether `50.0` is close to that or off by three orders of magnitude.

The issue is not that any of these knobs is *wrong* — a skilled user can tune them empirically to produce decent output. The issue is that **without a common physical scaling, you cannot tell whether the solver is hard because the physics is hard, or because the numbers are mutually incoherent and produce a problem the physics itself would never produce.**

### 2.3 Stiffness of $\eta(\dot\varepsilon)$

Power-law viscosity with $n = 3$ gives $\eta \propto \dot\varepsilon^{-2/3}$. Near the strain-rate floor $\varepsilon_\text{min}$, the derivative $\mathrm{d}\eta/\mathrm{d}\dot\varepsilon$ diverges as $\dot\varepsilon \to 0$. This is a genuine physical property of dislocation creep, not a bug — but it means the Jacobian of the nonlinear Stokes residual is highly sensitive in low-strain regions (plate interiors, cratons). The `smooth_saturate` upper bound and the `soft_min_harmonic` yielding blend (both introduced in #51) smooth out the *upper* kink at $\eta_\text{max}$ and the plastic corner; the *lower* end is still regularized by $\varepsilon_\text{min}$ additively.

Choosing $\varepsilon_\text{min}$ too small makes the Jacobian stiff in quiet regions (bad for Newton convergence); choosing it too large makes the rheology artificially linear everywhere (loss of localisation). The right value is tied to the physical velocity scale — specifically, $\varepsilon_\text{min} \ll v^*/L^* = 1/\tau^*$ but not by more than two orders of magnitude. Today the choice is disconnected from any scale.

### 2.4 Null space handling

Stokes with periodic BCs has a rank-2 null space: adding a constant to $v_x$ and/or $v_y$ does not change the residual. The code projects this out of the RHS and out of the iterate at every Newton step (and at the end of the BiCGSTAB solve). Good so far. **But the preconditioner (Jacobi or SSOR) is not null-space-aware.** It maps $r \mapsto z = M^{-1} r$ including any residual component in the null direction. This residual null component is then amplified by the Krylov basis construction and contaminates the search directions. On grids with very smooth initial conditions, this is a minor noise floor; on grids with sharp tectonic discontinuities (which is most realistic configurations), it can measurably slow BiCGSTAB. A cheap fix is to apply the null-space projector inside the preconditioner wrapper, before and after $M^{-1}$. Cost: two array-wide means per preconditioner application, i.e. $O(N)$ — negligible.

---

## 3. Characteristic scales and modelling intent

### 3.1 Primary scales

The goal is to pick a small set of physical scales that together span the whole tectonic problem, then derive everything else from them. We want scales that are (a) physically meaningful, (b) directly related to Living Landz gameplay requirements, and (c) chosen so that nondimensional fields are naturally $O(1)$ in the target regime.

| Symbol | Name | Value | Source |
|---|---|---|---|
| $L^\*$ | Horizontal length | $350$ km | Target continent diameter (TDD §2.1) |
| $S^\*$ | Crustal thickness | $35$ km | Reference continental crust (Turcotte & Schubert 2002) |
| $\tau^\*$ | Tectonic time scale | $25$–$30$ Myr | Single orogenic episode duration (calibration scale, not a run budget — see §3.3) |
| $\rho^\*$ | Density | $3300$ kg/m³ | Mantle density — buoyancy reference |

These four scales are *chosen*; everything else is *derived*.

### 3.2 Derived scales

$$
\begin{aligned}
v^\* &= L^\* / \tau^\* \approx 1.2\ \text{cm/yr} &&\text{velocity scale} \\
\dot\varepsilon^\* &= v^\*/L^\* = 1/\tau^\* \approx 10^{-15}\ \text{s}^{-1} &&\text{strain-rate scale} \\
\eta^\* &= \rho^\* g \tau^\* S^\* \approx 10^{24}\ \text{Pa·s} &&\text{viscosity scale} \\
\sigma^\* &= \eta^\* \dot\varepsilon^\* \approx 10^9\ \text{Pa} &&\text{stress scale} \\
p^\* &= \rho^\* g S^\* &&\text{lithostatic pressure scale}
\end{aligned}
$$

The definition $\eta^\* = \rho^\* g \tau^\* S^\*$ matches the standard convention of England & McKenzie (1982) and aligns with the accepted order of magnitude for effective lithospheric viscosity at long wavelengths. Combined with $v^\* = L^\*/\tau^\*$, it has a direct consequence for the Argand number derived in §4.1:

$$
\mathrm{Ar} \;=\; \frac{\rho^\* g (S^\*)^2}{\eta^\* v^\*} \;=\; \frac{\rho^\* g (S^\*)^2}{\rho^\* g \tau^\* S^\* \cdot L^\*/\tau^\*} \;=\; \frac{S^\*}{L^\*}.
$$

In thin-sheet geometry ($S^\* \ll L^\*$), $\mathrm{Ar}$ is therefore **necessarily small** — with the default primary scales, $\mathrm{Ar} = 35/350 = 0.1$. The statement that $\mathrm{Ar}$ is "$O(1)$ by construction" is only true if one adopts non-thin-sheet scales ($S^\* \sim L^\*$), which contradicts the model's own geometric assumption. This document's earlier drafts had an inconsistency between §3.2 (scales that force $\mathrm{Ar} = S^\*/L^\*$) and §5.1 (target range $\mathrm{Ar} \in [1, 5]$); §5.1 has been corrected.

$v^\* \approx 1.2$ cm/yr is deliberately at the low end of observed plate velocities (1–10 cm/yr). The continent is quiescent on average; peaks of $\tilde v \sim 5$–$10$ are allowed in active boundary zones. Choosing $v^\*$ at the low end keeps $\tilde v$ near unity most of the time and preserves the physical meaning of "one" as "typical motion".

### 3.3 Modelling intent: character vs geography, and the freedom of history

A point that is central to the whole nondimensionalization effort and easy to miss:

**The dimensionless numbers set the character of the world, not the geography of any particular continent.**

- **Character** — the *quality* of the output: how sharp boundaries are, how tall mountains get, how stable cratons are, whether deformations localize or distribute, whether old scars are reused or forgotten. This is what the dimensionless numbers control.
- **Geography** — the *specific* output: where mountains sit, how many there are, what shape the coastline takes, which plates collide with which. This is determined by the stochastic Voronoï initialisation (number of plates, types, velocity vectors, seed).

A single parametric configuration $(\mathrm{Ar}, \mathrm{Bi}, \mathrm{De}_p, K, \mathrm{Br}, \mathrm{Sp}, \mathrm{Mf}, \ldots)$ produces **many different continents** when sampled over different seeds — but all those continents share the same *character*. One seed may produce a convergent orogen along a northern margin; another seed with the same parameters but different Voronoï produces a rift opening through the interior. Both are valid outcomes of the same physics.

This has two consequences:

- There is no such thing as a parameter setting that "produces type X mountain chain". The mountain typology that emerges depends on the *local boundary configuration* (two continental plates converging → collision-type orogen; oceanic plate subducting under continental → arc-type orogen; divergent continental → rift), which is a geometric accident of the Voronoï seed.
- Validation must be **stochastic** (§3.4): we sample seeds, not single configurations.

#### Simulation duration is a user choice, not a model property

The total simulated geological time $T$ is **not** fixed by the physics or by the scaling. $\tau^\*$ is a *calibration scale* (a reference unit in which dimensionless numbers are expressed) — it is not a budget. The user chooses $T$ in multiples of $\tau^\*$ depending on the narrative intent:

- **Young continent** ($T \approx 1$–$3\,\tau^\*$, ~30–90 Myr). Single active orogenic phase still in progress at output time. One dominant boundary type. Topography is raw, not yet relaxed.
- **Mature continent** ($T \approx 4$–$7\,\tau^\*$, ~100–200 Myr). Multiple boundary types coexist at output time: one margin is active (arc or collision), another has already relaxed, interior may show older sutures. This is the "standard realistic" continent.
- **Fragmenting continent** ($T \approx 8$–$12\,\tau^\*$, ~250–350 Myr). Long enough for a rift to propagate across the continent, opening a new ocean and producing multiple island masses with coherent geological histories. The Voronoï topology evolves substantially (thanks to `dynamic_boundaries`), possibly producing two or more landmasses from an initial single continent.
- **Very old world** ($T \geq 15\,\tau^\*$, ~450+ Myr). Multiple Wilson-cycle-like episodes on the same footprint, deep plastic memory (if $\mathrm{De}_p$ permits healing, else saturated scar patterns). Expected to be numerically demanding — candidate for the hardest validation case.

What changes with $T$ is **how much history accumulates**, not the qualitative character of the physics. Long $T$ with low $\mathrm{De}_p$ (no healing) saturates the continent with scars; long $T$ with moderate $\mathrm{De}_p$ produces belts of ages with mobile corridors; long $T$ with high $\mathrm{De}_p$ keeps the continent perpetually "young-looking" regardless of actual age. The age field (§4.11) records the history in any case so downstream phases can use it for differential erosion, climate, resource placement, or visualisation.

The Phase E stochastic validation protocol (§3.4) runs at multiple $T$ values to ensure numerical stability across the full range of intended narrative durations.

### 3.4 Stochastic validation strategy

Given §3.3, the relevant validation question is not *"does this parameter set produce a particular mountain type?"* but *"does this parameter set produce a coherent diversity of continents across seeds, all numerically stable, at the intended range of simulation durations?"*.

The proposed validation protocol:

1. Fix a "healthy" reference parameter set $\mathcal{P}_\text{ref}$ (see §5 for candidate values).
2. Sample $N = 20$ Voronoï seeds with varying plate counts ($5$–$15$), continental ratios ($20\%$–$50\%$), and random velocity vectors (magnitudes $\tilde v_0 \in [0.5, 2.0]$).
3. Run each seed through the full tectonic phase at multiple durations $T \in \{3, 6, 10\}\,\tau^\*$ on both $128^2$ and $512^2$ grids.
4. Collect per-seed metrics:
   - **Numerical health**: Newton outcomes distribution, mean BiCGSTAB iterations, estimated $\kappa(A)$, fraction of sub-steps that had to shrink.
   - **Physical plausibility**: peak $\tilde S$, fraction of cells in plastic regime, land/ocean ratio, presence of distinct boundary types.
   - **Diversity**: spread across seeds — we want the metrics to *vary* between seeds (different continents) but stay in physical bounds.
   - **Duration scaling**: metrics should evolve smoothly with $T$ (no blow-up, no collapse to triviality at long $T$).
5. The parameter set passes if: all seeds×durations converge; metrics distribution is unimodal without outliers; at least 80% of seeds exhibit two or more distinct boundary typologies at $T \geq 6\,\tau^\*$.

No individual seed is a "correct" benchmark. A **regression** is a shift in the *distribution* of metrics, not a change in any single run.

This protocol replaces single-config benchmarks. Single-config tests remain in the suite as smoke tests (verify the solver runs, no NaN, converges on a trivial case) but they are no longer the validation standard.

---

## 4. Nondimensionalized equations

We introduce tildes for nondimensional variables throughout:

$$
\tilde x = x/L^\*, \quad \tilde t = t/\tau^\*, \quad \tilde v = v/v^\*, \quad \tilde S = S/S^\*, \quad \tilde\eta = \eta/\eta^\*, \quad \tilde\rho = \rho/\rho^\*.
$$

Tildes drop on operators where convenient: $\tilde\nabla = L^\* \nabla$, $\tilde\nabla^2 = (L^\*)^2 \nabla^2$. Tildes also drop on the dimensionless strain rate $\tilde{\dot\varepsilon} = \dot\varepsilon/\dot\varepsilon^\* = \dot\varepsilon\,\tau^\*$.

### 4.1 Stokes momentum balance

The dimensional momentum equation (neglecting inertia — Reynolds $\ll 1$):

$$
\nabla \cdot \boldsymbol\tau - \nabla \Phi = \mathbf{f}_\text{ext},
$$

where $\boldsymbol\tau = 2\eta \boldsymbol{\dot\varepsilon}$ is the deviatoric stress, $\Phi$ is the gravitational potential energy per unit area (depth-integrated), and $\mathbf{f}_\text{ext}$ gathers external tractions (plate driving, slab pull, basal drag). Substituting scales:

$$
\frac{\eta^\* \dot\varepsilon^\*}{L^\*} \tilde\nabla \cdot (2\tilde\eta \,\tilde{\boldsymbol{\dot\varepsilon}}) \;-\; \frac{\rho^\* g (S^\*)^2}{L^\*} \tilde\nabla \tilde\Phi \;=\; \frac{[\mathbf{f}_\text{ext}]}{L^\*}.
$$

Dividing through by $\eta^\* \dot\varepsilon^\* / L^\*$:

$$
\boxed{\;\tilde\nabla \cdot (2\tilde\eta \,\tilde{\boldsymbol{\dot\varepsilon}}) \;-\; \mathrm{Ar}\, \tilde\nabla \tilde\Phi \;=\; \tilde{\mathbf{f}}_\text{ext}.\;}
$$

The **Argand number** is

$$
\mathrm{Ar} \;=\; \frac{\rho^\* g (S^\*)^2 / L^\*}{\eta^\* \dot\varepsilon^\*} \;=\; \frac{\rho^\* g (S^\*)^2}{\eta^\* v^\*}.
$$

**Numerical value under the default scales.** With the derived $\eta^\* = \rho^\* g \tau^\* S^\*$ and $v^\* = L^\*/\tau^\*$ of §3.2, this expression simplifies to $\mathrm{Ar} = S^\*/L^\*$, i.e. the aspect ratio of the thin viscous sheet. For the default primary scales ($S^\* = 35$ km, $L^\* = 350$ km), $\mathrm{Ar} = 0.1$. This is **not $O(1)$ as earlier drafts claimed** — in thin-sheet geometry $\mathrm{Ar}$ is necessarily small. The consequences for the time-scale hierarchy of the full system are discussed in §5.4.

$\mathrm{Ar}$ measures the competition between gravitational spreading of thickened crust and viscous resistance. Crucially, $\mathrm{Ar}$ **sets the equilibrium thickness** at which GPE spreading balances convergent input: a higher $\mathrm{Ar}$ means crust spreads more readily and equilibrates at lower $\tilde S_\text{eq}$; a lower $\mathrm{Ar}$ allows thicker equilibrium orogens. The value of $\tilde S_\text{eq}$ emerging from the solver is therefore a direct diagnostic of whether $\mathrm{Ar}$ is correctly calibrated — see §7 Q1.

The GPE potential nondimensionalizes as $\tilde\Phi = \tilde\rho (1 - \tilde\rho/\tilde\rho_m) \tilde S^2$ (or simply $\tilde S^2$ in the uniform-density limit).

### 4.2 Power-law rheology

Dimensional:

$$
\eta = B \,\dot\varepsilon^{\,1/n - 1},
$$

with $B$ a prefactor carrying units of $\text{Pa} \cdot \text{s}^{1/n}$. Nondimensional:

$$
\tilde\eta = \tilde B \,\tilde{\dot\varepsilon}^{\,1/n - 1}, \qquad \tilde B = B (\dot\varepsilon^\*)^{1/n - 1} / \eta^\*.
$$

Setting $\tilde B = 1$ fixes the prefactor: $B = \eta^\* (\dot\varepsilon^\*)^{1 - 1/n}$. The rheology becomes parameter-free once scales are chosen, modulo the choice of $n$ and the regularization.

**Regularization.** The additive floor $\varepsilon_\text{min}$ becomes $\tilde\varepsilon_\text{min}$ in non-dim. Target: $\tilde\varepsilon_\text{min} \in [10^{-3}, 10^{-2}]$. Below $10^{-4}$ the Jacobian stiffness of §2.3 makes Newton painful; above $10^{-2}$ the rheology is linearised in most of the domain.

**Saturation.** The upper saturation at $\tilde\eta_\text{max}$ is a numerical safeguard reflecting the fact that real rocks become effectively rigid below a temperature-dependent strain rate. Target: $\tilde\eta_\text{max} \in [10^2, 10^3]$, not $10^4$.

**Floor.** $\tilde\eta_\text{min} \sim 0.1$–$1$. Active plastic yield zones can locally soften, but there is no physical regime below $\eta^\*/10$.

**Resulting contrast.** $\tilde\eta_\text{max} / \tilde\eta_\text{min} \in [10^2, 10^4]$, versus $10^7$ today. At $128^2$, this drops $\kappa(A)$ from $\sim 10^{11}$ to $\sim 10^{8}$ — a factor $10^{1.5}$ in BiCGSTAB iterations before any preconditioner improvement.

### 4.3 Plastic yielding and plastic memory

Drucker-Prager-like yielding is currently

$$
\eta_\text{plastic} = \frac{\tau_Y}{2 \dot\varepsilon_{II}}, \qquad \eta = \mathrm{soft\_min}(\eta_\text{visc}, \eta_\text{plastic}).
$$

Nondimensional:

$$
\tilde\eta_\text{plastic} = \frac{\mathrm{Bi}}{2 \tilde{\dot\varepsilon}_{II}}, \qquad \mathrm{Bi} = \frac{\tau_Y}{\eta^\* \dot\varepsilon^\*} = \frac{\tau_Y}{\sigma^\*}.
$$

The **Bingham number** $\mathrm{Bi}$ sets the yield threshold relative to typical viscous stress. Target: $\mathrm{Bi} \in [0.05, 0.5]$. Real lithosphere: $\tau_Y \approx 100$ MPa, $\sigma^\* \approx 10^9$ Pa, so $\mathrm{Bi} \approx 0.1$.

**Plastic memory.** Each cell accumulates $\int \dot\varepsilon_{II}\,\mathrm{d}t$ while in the plastic regime (field `plastic_strain`). This accumulation has two effects: it weakens the local yield stress ($\tau_Y \to \tau_Y(1 - w)$ with $w$ proportional to accumulated strain), so a cell that has plastified once plastifies more easily next time — this is a *scar* in the rock. The healing term $-\mathrm{d}\epsilon_p/\mathrm{d}t = r_\text{healing}$ reduces this memory over time. The nondimensional combination

$$
\mathrm{De}_p = \tau^\* \cdot r_\text{healing}
$$

is a Deborah number for plastic memory: it answers "how many tectonic time scales to erase a scar?".

- $\mathrm{De}_p = 0$: no healing, scars forever (current default `healing_rate = 0`). Successive deformation episodes preferentially reuse the same weak zones → a single dominant suture concentrates everything.
- $\mathrm{De}_p \gg 1$: scars erase immediately, plastic memory plays no role → successive episodes pick new locations each time → diffusely deformed continent.
- $\mathrm{De}_p \sim 0.1$–$0.5$: partial reuse → belts of successive ages, stable cores alternating with mobile zones — the realistic case.

**This is one of the highest-leverage knobs on the qualitative diversity of the output** and deserves dedicated exploration. Plastic memory is also the mechanism used in §4.10 to implement cratonic stability without piling onto the viscosity contrast.

### 4.4 Basal drag (mantle drag)

Dimensional:

$$
\mathbf{f}_\text{drag} = -C_b \,\rho^\* g\, S\, \mathbf{v}.
$$

The current code uses $C_b / \Delta x^2$ scaling for resolution independence; that is a discretisation choice, not a physical one. The physical coefficient $C_b$ has units of $1/\text{length}$ (asthenospheric viscosity over a reference depth).

Nondimensional:

$$
\tilde{\mathbf{f}}_\text{drag} = -\mathrm{Br} \,\tilde S\, \tilde{\mathbf{v}}, \qquad \mathrm{Br} = \frac{C_b \rho^\* g (S^\*)^2}{\eta^\* v^\* / L^\*}.
$$

$\mathrm{Br}$ measures the ratio of basal traction to internal viscous dissipation. Physical estimate: asthenospheric viscosity $\eta_a \sim 10^{20}$ Pa·s over a $h_a \sim 100$ km reference depth gives $\mathrm{Br} \approx 0.01$–$0.1$. Target: $\mathrm{Br} \in [0.01, 0.3]$.

The `S >= 0.3` threshold that skips friction on thin oceanic crust (commit `ab1ac63`) is a fix for a symptom — in the current scheme, thin oceanic plates over-decelerate. A cleaner formulation is drag proportional to $\tilde S^2$ or a smooth transition rather than a hard on/off. This eliminates the kink and is physically closer to the truth.

### 4.5 Gravitational potential energy (GPE) spreading

Already covered in §4.1 via $\tilde\Phi$. For completeness, the current code uses

$$
\Phi = \rho_c \left(1 - \frac{\rho_c}{\rho_m}\right) g S^2 / 2 \;\;\Rightarrow\;\; \tilde\Phi = \tilde\rho (1 - \tilde\rho/\tilde\rho_m) \tilde S^2.
$$

Oceanic cells have $\tilde\Phi/\tilde S^2 \approx 0.08$, continental $\approx 0.14$: a factor $\sim 1.75$ less spreading pressure for oceanic than continental, matching the physical intent. This is fine as-is; the Argand number $\mathrm{Ar}$ absorbs the scalar prefactor.

### 4.6 Crustal thickness advection

Dimensional:

$$
\frac{\partial S}{\partial t} + \nabla \cdot (S \mathbf{v}) = Q(x, t).
$$

Nondimensional (divide by $S^\*/\tau^\*$):

$$
\boxed{\;\frac{\partial \tilde S}{\partial \tilde t} + \tilde\nabla \cdot (\tilde S \tilde{\mathbf{v}}) = \tilde Q.\;}
$$

Parameter-free by construction. $\tilde Q$ is the nondimensional source rate (§4.7).

**CFL.** The advection stability bound becomes

$$
\Delta \tilde t \leq C_\text{CFL} \cdot \frac{\Delta \tilde x}{\max |\tilde{\mathbf{v}}|},
$$

with $C_\text{CFL} \in [0.3, 0.5]$ for first-order upwind (current `cfl_factor = 0.5` is at the edge, $0.3$ is safer with nonlinear source feedback). In physical units, a macro step at $128^2$ with typical $\tilde v \sim 1$ gives $\Delta \tilde t \approx 0.004$, so $\Delta t \approx 0.12$ Myr. Number of steps required for a run of duration $T$ is then $T/\Delta t$, which for $T = 6\,\tau^\*$ gives roughly $1500$ macro steps. The adaptive sub-stepping already targets this regime; what is new is that the budget is now expressed in Myr, not in arbitrary solver ticks.

### 4.7 Source/sink terms at plate boundaries

These are the opaquest parameters. Let us nondimensionalize each.

**Subduction consumption.** $Q_\text{sub} = -k_\text{sub} \cdot \|\Delta\mathbf{v}\| / \Delta x$ on the subducting oceanic cell. The rate is already dimensionless when expressed as *fraction of column consumed per unit of convergent motion $L^\*$*. Target: $\tilde k_\text{sub} \in [0.3, 1.0]$. Current default $0.5$ lands in this range.

**Volcanic arc resurfacing.** Fraction of subducted mass that resurfaces. Physical estimate from arc magma flux: $\sim 10$–$30\%$ of trench consumption. Target: $\tilde k_\text{arc} \in [0.05, 0.3]$. Current $0.15$ is centre-of-range.

**Oceanic spreading.** $\tilde k_\text{spread}$ is chosen such that the steady-state ridge thickness is $\tilde S_\text{oceanic} \approx 0.2 \pm 10\%$. This is a *solver calibration test*, not a user-facing parameter: pick $\tilde k_\text{spread}$ by closure with the thickness field, don't treat it as an independent knob.

**Collision volcanism.** Minor source at continental collision zones. Target: $\tilde k_\text{coll-v} \in [0.02, 0.1]$, clearly smaller than arc volcanism.

**Rift volcanism.** Minor source at rift zones with thick crust. Target: $\tilde k_\text{rift-v} \in [0.01, 0.05]$.

The **relative ordering** is physically constrained: $\tilde k_\text{sub} > \tilde k_\text{arc} > \tilde k_\text{coll-v}, \tilde k_\text{rift-v}$. Any run violating this ordering is physically inconsistent. This gives a cheap *parameter sanity check* at solver initialisation.

### 4.8 Slab pull

$$
\mathbf{f}_\text{slab} = k_\text{slab} \cdot m_\text{subducted} \cdot \hat{\mathbf{n}}_\text{convergence}.
$$

Nondimensional:

$$
\tilde{\mathbf{f}}_\text{slab} = \mathrm{Sp} \cdot \tilde m_\text{subducted} \cdot \hat{\mathbf{n}}_\text{convergence}, \qquad \mathrm{Sp} = \frac{k_\text{slab} M^\*}{\eta^\* \dot\varepsilon^\* / L^\*},
$$

where $M^\*$ is a reference subducted mass scale.

**Bounded-growth issue.** The current formulation has $m_\text{subducted}$ growing monotonically, so $\mathbf{f}_\text{slab}$ grows without bound and is only stopped by the hard `max_plate_velocity = 5.0` clamp. The hard clamp is a red flag — the physics should self-limit.

Proposed fix: exponential decay representing slab detachment and mantle residence,

$$
\frac{\partial \tilde m}{\partial \tilde t} = \tilde Q_\text{sub-conv} - \frac{\tilde m}{\tilde\tau_\text{slab}}, \qquad \tilde\tau_\text{slab} \approx 0.3\text{–}1.0 \;(\sim 10\text{–}30\ \text{Myr}).
$$

This eliminates the velocity clamp entirely. Moved to §6.

**Activation regime requires pre-existing non-quiescent flow.** The nondimensional $\mathrm{Sp}$ defined above assumes $\eta \sim O(1)$. In the floor-dominated regime where $\dot\varepsilon_{II} < \dot\varepsilon_\text{min}$ everywhere, the effective viscosity is $\eta_\text{newton} = \dot\varepsilon_\text{min}^{1/n-1} \approx 100$, and the closed-loop gain $G = \mathrm{Sp}\cdot k_\text{slab-accum}\cdot \tau_\text{slab}\,/\,(\eta\cdot L)$ stays $\ll 1$ over the $\mathrm{Sp} \in [0.5, 3]$ range. The quiescent fixed point is linearly stable — slab-pull alone cannot bootstrap the system out of floor-domination. Activation requires an external initiator (mantle forcing §4.9) that imposes flow independently of the local closed-loop gain. Consequence: **slab-pull is an amplifier, not an initiator**; it transforms pre-existing convergence into traction but cannot create convergence where none exists. This was demonstrated empirically at the Step 7 baseline (see `step7_physics_report.md §Yielding checkpoint`) — `peak|v|` unchanged from Step 6 across the full $\mathrm{Sp}$ band, matching the loop-gain prediction $G \sim 7\times 10^{-3}$.

### 4.9 Mantle convection forcing

The mantle convection proxy applies a continuous velocity bias to plates:

$$
\tilde{\mathbf{v}}_\text{mantle}(\tilde x, \tilde t) = \mathrm{Mf} \cdot \tilde{\mathbf{v}}_\text{pattern}(\tilde x, \tilde t),
$$

where $\tilde{\mathbf{v}}_\text{pattern}$ is an $O(1)$ flow pattern and $\mathrm{Mf}$ is the forcing magnitude. Target: $\mathrm{Mf} \in [0.3, 2.0]$.

### 4.10 Cratonic rigidity — split-mechanism design

Current implementation: multiplicative on $\tilde\eta$,

$$
\tilde\eta_\text{cratonic}(\tilde x) = K \cdot \tilde\eta_\text{base}(\tilde x), \qquad K \in [1, 10^3]\ \text{in principle}.
$$

Physically, real cratonic roots are $10^2$–$10^3$ stiffer than surrounding mobile belts. Current `max_factor = 5.0` is very low by geological standards, but empirical observation confirms that raising $K$ substantially causes Newton to fail — the $\kappa(A)$ budget is already saturated.

**Decision**: decouple "cratons resist deformation" from "cratons have high viscosity". The primary mechanism for cratonic stability becomes **plastic memory immunity**, with viscosity contrast as a secondary contributor.

**Primary mechanism — plastic memory template.** Cratons are regions where the plastic weakening never accumulates: either `plastic_strain` is pinned to zero, or `weakening_strain_ref` is set high enough that realistic stresses cannot weaken the yield threshold. The yield stress stays at its full $\mathrm{Bi}$ value throughout the run, so cratons do not yield under stresses that would soften a mobile belt. They can still flex elastically (through the viscous branch), but they do not localize, do not develop shear zones, and do not accumulate damage.

- Pro: no $\kappa(A)$ penalty — plastic memory is a diagonal-in-memory mechanism, not a multiplicative modification of $\eta$.
- Pro: matches physical intuition better. Real cratons are strong *because* their cold, thick roots prevent ductile localisation, not only because their mean viscosity is higher.
- Pro: interacts naturally with $\mathrm{De}_p$ — a craton is a region where, in addition to initial zero strain, the healing rate is effectively infinite (or equivalently, `weakening_strain_ref` is infinite).

**Secondary mechanism — modest $K$.** A small viscosity contrast ($K \in [3, 8]$) is retained as a bulk stiffening that keeps cratons from participating in wide-wavelength flow. This does cost some $\kappa(A)$ but the cost is manageable because $K$ stays well below $10$.

**Admissible combinations** (to be validated in Phase E):
- $K = 5$, plastic immunity on → strong cratonic character, $\kappa(A)$ only slightly degraded.
- $K = 3$, plastic immunity on → softer cratons, still recognisable.
- $K = 1$, plastic immunity on → cratons differ from mobile belts only by plastic behaviour.

**Thermal coupling (§7 Q3)** is the physically complete answer — cratons are cold, therefore Arrhenius viscosity is naturally high there, and there is no need to encode rigidity explicitly. But thermal coupling is a larger project and is deferred.

### 4.11 Geological age field

The tectonic solver maintains a scalar field $A(\tilde x, \tilde t)$ representing the **age of the crustal material** at each cell. This is a new field, not currently in the code, added because it enables several downstream phases (erosion, climate, resource placement) and costs almost nothing.

**Evolution.**

$$
\frac{\partial A}{\partial \tilde t} + \tilde{\mathbf{v}} \cdot \tilde\nabla A = \tilde\Gamma(A, \text{events}),
$$

where the advective term transports age with the crust (a cell moved by plate motion keeps its age) and $\tilde\Gamma$ is a *resetting* term triggered by boundary events:

- **New oceanic crust at ridge.** Cell age reset to $A = 0$.
- **Volcanic arc resurfacing.** Cell age reset to $A = 0$ (freshly emplaced volcanic material).
- **Collision thickening.** The resulting cell inherits the *maximum* age of the contributing cells (the scar is as old as the older protolith). Alternative: weighted by mass contribution — to discuss, probably not worth the complexity at this stage.
- **Consumed cells (subduction).** Removed entirely; their age disappears with them.
- **Quiescent cells.** $\tilde\Gamma = 0$ — age grows linearly with $\tilde t$ through the advection equation alone (reduced to $\partial A/\partial \tilde t = 1$ in Lagrangian terms).

**Initialisation.** At $\tilde t = 0$, cells are assigned an age representing "basement age": continental cells start at a large value (representing pre-simulation cratonic age, e.g. $A_0 = 5$–$10\,\tau^\*$), oceanic cells at a smaller value reflecting ridge proximity. This is a user-configurable distribution and does not need great precision — the simulation overwrites it through events.

**Transport scheme.** Same upwind scheme as $S$ advection (§4.6), with the same CFL limit. Source terms are applied per-cell after advection, triggered by the boundary classification already computed by `BoundaryField`.

**Output.** $A(\tilde x, T)$ is exported alongside the heightmap. Units: dimensionless multiples of $\tilde t$ at export, converted to Myr using $\tau^\*$ if a physical output is desired.

**Downstream uses.**

- **Differential erosion (Phase 4).** Older chains are more eroded at a given altitude. Concretely, the erosion rate on a cell can be weighted by a factor $f(A)$ that is larger for older material — or the full erosion simulation can be replayed for a duration proportional to age, starting from the age-adjusted initial condition. To be specified in the erosion phase note.
- **Weathering and soil maturity (Phase 5 climate).** Older substrates have more developed soils, affect vegetation biomes.
- **Geological classification (Phase 6 export).** The current `geology.raw` export (TDD §9.1) is based on tectonic *event class* (shield, collision, subduction, rift). Adding age lets the consumer distinguish "young collision" from "ancient collision".

The age field is a universal descriptor that unifies "this zone is a young Alpine analogue" and "this zone is a worn Appalachian analogue" into a single quantitative parameter, with the typology (arc, collision, rift) carried separately by the boundary classification.

---

## 5. Dimensionless numbers and continental character

### 5.1 Reference table

| Symbol | Name | Definition | Target range |
|---|---|---|---|
| $\mathrm{Ar}$ | Argand | $\rho^\* g (S^\*)^2 / (\eta^\* v^\*) = S^\*/L^\*$ | derived, thin-sheet $\ll 1$ |
| $\mathrm{Bi}$ | Bingham | $\tau_Y / \sigma^\*$ | $0.05$–$0.5$ |
| $\mathrm{De}_p$ | Plastic Deborah | $\tau^\* \cdot r_\text{healing}$ | $0.1$–$0.5$ |
| $\mathrm{Br}$ | Basal drag | $C_b \rho^\* g (S^\*)^2/(\eta^\* v^\*/L^\*)$ | $0.01$–$0.3$ |
| $K$ | Cratonic viscous contrast | $\eta_\text{craton}/\eta^\*$ | $3$–$8$ (secondary, see §4.10) |
| $\mathrm{Sp}$ | Slab pull | $k_\text{slab} M^\*/(\eta^\* \dot\varepsilon^\*/L^\*)$ | $0.5$–$3$ |
| $\mathrm{Mf}$ | Mantle forcing | velocity bias amplitude | $0.3$–$2.0$ |
| $\tilde\eta_\text{max}/\tilde\eta_\text{min}$ | Viscosity contrast | ratio | $\leq 10^3$ |
| $\tilde\varepsilon_\text{min}$ | Strain-rate floor | $\varepsilon_\text{min} \tau^\*$ | $10^{-3}$–$10^{-2}$ |
| $\tilde k_\text{sub}$ | Subduction rate | fraction consumed per $L^\*$ of convergence | $0.3$–$1.0$ |
| $\tilde k_\text{arc}$ | Arc resurfacing | fraction of subducted mass returned | $0.05$–$0.3$ |
| $\tilde k_\text{spread}$ | Spreading rate | calibrated to give $\tilde S_\text{ocean} \approx 0.2$ | calibrated |
| $n$ | Power-law exponent | — | $1 \to 3$ via continuation |
| $T/\tau^\*$ | Simulation duration | user choice | $1$–$15+$ |

$\mathrm{Ar}$ is **derived** from the four primary scales, not a free knob: $\mathrm{Ar} = S^\*/L^\*$. With default scales ($S^\* = 35$ km, $L^\* = 350$ km) this gives $\mathrm{Ar} = 0.1$, and any thin-sheet geometry ($S^\* \ll L^\*$) forces $\mathrm{Ar} \ll 1$. It must be reported as a diagnostic at solver startup; the user adjusts it only indirectly, by changing the primary scales. See §5.4 for the consequences of $\mathrm{Ar} \ll 1$ on the relative dynamics of the system.

### 5.2 Qualitative effects on continental character

For each of the highest-leverage numbers, what varying it does to the *character* of continents produced by the solver (averaged over seeds — individual outputs always vary with the Voronoï seed):

**Argand $\mathrm{Ar}$** — "how fast does thickened crust spread?"

With $\mathrm{Ar} = S^\*/L^\*$ fixed by the primary scales (§3.2), Ar is not a direct character knob but a diagnostic of the GPE time scale relative to the tectonic one: $\tau_\text{GPE}/\tau^\* = 1/\mathrm{Ar}$. In thin-sheet geometry $\mathrm{Ar} \ll 1$, meaning GPE spreading is the **slow** mechanism of the system (see §5.4).

To tune "how rapidly crust spreads" within the model, adjust the primary scales rather than $\mathrm{Ar}$ itself:
- **Thinner crust or wider domain** ($S^\*$ down, $L^\*$ up): $\mathrm{Ar}$ smaller → GPE slower, other mechanisms dominate more clearly.
- **Thicker crust or narrower domain**: $\mathrm{Ar}$ larger → GPE faster relative to tectonic motion, spreading more visible on short runs.

Earlier drafts of this note listed $\mathrm{Ar} \in [1, 5]$ as a target range; this was inconsistent with the thin-sheet geometry fixed by §3.2 and has been retired.

**Bingham $\mathrm{Bi}$** — "does deformation localize or distribute?"
- Low ($0.05$–$0.1$): easy yielding, deformation localizes in faults and shear zones → **sharp tectonic boundaries**, narrow deformed belts separating rigid blocks.
- High ($0.3$–$0.5$): yielding rare, deformation stays viscous → **diffuse deformation**, gentle folds, no sharply delimited provinces.

**Plastic Deborah $\mathrm{De}_p$** — "are old scars remembered?"
- $0$ (current): scars forever → **single dominant suture** that absorbs all subsequent deformation → bipolar continent (one big orogen + large undeformed domains).
- $0.1$–$0.3$: partial memory → **belts of successive ages**, mobile corridors alongside stable cores — the canonical geological map.
- $> 1$: scars vanish quickly → **diffusely accidented** continent, no coherent long-term structural pattern.

**Basal drag $\mathrm{Br}$** — "how tightly does the mantle hold the plates?"
- Low ($0.01$): plates move freely → **fast tectonic evolution**, initial conditions quickly forgotten, dynamic geography.
- High ($0.3$): plates strongly dragged → **slow evolution**, initial plate config persists long, more static geography.

**Cratonic system (plastic immunity + modest $K$)** — "are there stable ancient cores?"
- Disabled (no immunity, $K = 1$): no protected zones → **uniform deformation**, no recognisable cratonic provinces.
- Enabled ($K = 5$, plastic immunity on): cratons push deformation outward → **mobile belts at craton margins**, stable old interiors, scars concentrate at craton boundaries.

**Slab pull $\mathrm{Sp}$** — "do plates self-accelerate toward convergence?"
- Low ($0.5$): plates drift by initial velocities only → **static geography**, initial Voronoï setup dominates outcome.
- High ($3$): convergence self-reinforces → **convergence-dominated continent** with sustained subduction zones, few rifts.

**Mantle forcing $\mathrm{Mf}$** — "is the continent continuously driven?"
- Low: plates settle into a steady state after a few episodes → **tectonic quiescence** at late times.
- High: continuous activity → **sustained tectonics** throughout the run, more episodes within a given time budget.

**Duration $T/\tau^\*$** — "how long is the story?"
- Short ($1$–$3\,\tau^\*$): **single-phase continent**, raw topography, one dominant boundary type.
- Medium ($4$–$7$): **mature continent** with coexisting active and relaxed regions.
- Long ($8$–$12$): **fragmenting or multi-cycle continent**, possibly multiple landmasses, complex age stratification.
- Very long ($\geq 15$): **saturated history**, numerical stress test.

### 5.3 Character presets (hypotheses, to be validated)

Three reference configurations worth empirically validating. Each is a hypothesis about a coherent combination — §6 Phase E is the experiment that confirms or refutes them. Duration $T/\tau^\*$ is kept as a separate dimension: each preset should be tested at short, medium, and long durations.

**Preset A — "Dynamic accidented" (default candidate)**
Mobile, varied, visually rich. Mountains of mixed ages, active tectonics at output time, a few cratonic cores.

| $\mathrm{Ar}$ | $\mathrm{Bi}$ | $\mathrm{De}_p$ | $\mathrm{Br}$ | $K$ | plastic immunity | $\mathrm{Sp}$ | $\mathrm{Mf}$ |
|---|---|---|---|---|---|---|---|
| 2 | 0.15 | 0.3 | 0.05 | 5 | on | 1.5 | 1.0 |

**Preset B — "Stable shield with active margins"**
Older feel: a large stable interior, deformation concentrated at edges.

| $\mathrm{Ar}$ | $\mathrm{Bi}$ | $\mathrm{De}_p$ | $\mathrm{Br}$ | $K$ | plastic immunity | $\mathrm{Sp}$ | $\mathrm{Mf}$ |
|---|---|---|---|---|---|---|---|
| 1.5 | 0.1 | 0.15 | 0.1 | 8 | strong | 2.0 | 0.5 |

**Preset C — "Soft planet" (exploratory)**
Low-relief world: broad orogens, smoother topography, gentle diffuse deformation.

| $\mathrm{Ar}$ | $\mathrm{Bi}$ | $\mathrm{De}_p$ | $\mathrm{Br}$ | $K$ | plastic immunity | $\mathrm{Sp}$ | $\mathrm{Mf}$ |
|---|---|---|---|---|---|---|---|
| 4 | 0.3 | 0.4 | 0.03 | 3 | off | 0.8 | 1.2 |

**Important caveat.** These tables are *hypotheses* about which parameter combinations produce coherent mental pictures. They are not validated. §6 Phase E is the experiment that turns hypotheses into known working presets. If a preset fails to produce the claimed character across its seed sample (or produces numerical instability), the table is wrong and gets revised.

### 5.4 Time-scale hierarchy and the "slow GPE" regime

A practical consequence of $\mathrm{Ar} = S^\*/L^\* \approx 0.1$ deserves to be made explicit because it governs the balance of the whole model.

Each physical mechanism in the solver has a characteristic time. Dividing them by $\tau^\*$ gives their dimensionless clock:

| Mechanism | Dimensional time | $/\tau^\*$ | Default value |
|---|---|---|---|
| Tectonic baseline | $\tau^\*$ (definition) | $1$ | 30 Myr |
| GPE spreading | $\tau^\*/\mathrm{Ar}$ | $1/\mathrm{Ar}$ | 300 Myr at $\mathrm{Ar} = 0.1$ |
| Plastic yielding onset | $\tau^\* \cdot \mathrm{Bi}$ | $\mathrm{Bi}$ | 1.5–15 Myr at $\mathrm{Bi} \in [0.05, 0.5]$ |
| Scar healing | $\tau^\*/\mathrm{De}_p$ | $1/\mathrm{De}_p$ | 60–300 Myr at $\mathrm{De}_p \in [0.1, 0.5]$ |
| Basal drag response | $\tau^\*/\mathrm{Br}$ | $1/\mathrm{Br}$ | 100–3000 Myr at $\mathrm{Br} \in [0.01, 0.3]$ |
| Slab-pull acceleration | $\tau^\*/\mathrm{Sp}$ | $1/\mathrm{Sp}$ | 10–60 Myr at $\mathrm{Sp} \in [0.5, 3]$ |
| Mantle forcing cycle | $\tau^\*/\mathrm{Mf}$ | $1/\mathrm{Mf}$ | 15–100 Myr at $\mathrm{Mf} \in [0.3, 2]$ |

**Key observation**: GPE spreading at $\mathrm{Ar} = 0.1$ is one of the slowest mechanisms (alongside basal drag and scar healing). It is **dominated** on short runs ($T < 3\,\tau^\*$) by yielding, slab-pull, and mantle forcing, which all operate on time scales close to or shorter than $\tau^\*$ itself.

**Why this matters:**

- Earlier versions of the solver (before this milestone) produced rapid GPE-driven collapse of thick crust in roughly 20 macro steps. That dynamic required an implicit $\mathrm{Ar} \sim 1$–$10$, which is inconsistent with thin-sheet scaling. The fast spreading that "looked right" was symptomatic of a miscalibrated force balance in which GPE overwhelmed subduction, collision, and slab-pull. Correcting $\mathrm{Ar}$ to its scale-derived value $0.1$ is part of the reconstruction milestone's objective.

- With the corrected $\mathrm{Ar}$, a thick-crust anomaly requires $O(10\,\tau^\*) \approx 300$ Myr to diffuse away by GPE alone. Over a typical mature continent run ($T = 6\,\tau^\*$), GPE does about $60\%$ of one spreading time scale — significant smoothing of long-wavelength features, but not a collapse. This matches geological observations of persistent orogenic roots over hundreds of Myr.

- **Dynamics on short runs must come from the fast mechanisms**: yielding localizes, slab-pull accelerates subducting plates, mantle forcing introduces continental motion, boundary sources create and destroy crust. GPE operates as a long-timescale relaxation over all of this.

**Implications for calibration:**

- If $T$ in a run is short and visible GPE dynamics are desired, the honest move is to extend $T$ rather than inflate $\mathrm{Ar}$.
- If, after all mechanisms are in place (Steps 3–10), the visible dynamics of runs at $T = 6\,\tau^\*$ are judged insufficient for the narrative purpose (Living Landz "playability"), the options are, in order of physical honesty:
  1. Shorten $\tau^\*$ (e.g., from 30 Myr to 10 Myr), which mechanically increases all dimensionless numbers that depend on it ($\mathrm{Ar}$, $\mathrm{De}_p$, etc.). Most honest — it acknowledges the user is running a "fast-forward" world.
  2. Decouple the GPE coefficient from $\mathrm{Ar}$ by introducing an explicit amplification factor. Less honest — it breaks the nondimensional self-consistency — but may be necessary if time compression alone does not suffice.
  3. Accept long $T$ (e.g., $T = 20\,\tau^\*$) as the default and design the pipeline to handle it. Cleanest physically; costs wallclock.

This decision is deferred to post-Step 4, when all fast mechanisms except slab-pull and mantle are in place and the short-run dynamics can be evaluated empirically.

---

## 6. Implementation roadmap

The goal is to migrate the solver to internal nondimensional variables **without breaking the current test suite** and **without changing user-visible output** (initially). User-facing scaling of parameters — i.e. replacing `subduction_rate = 0.5` with a dimensionless parameter expressing the same quantity — comes in a second pass once numerical equivalence is confirmed.

### Phase A — `Scales` infrastructure (non-breaking)

Add a `Scales` struct in `ymir-core/src/tectonics/scales.rs` exposing the four primary scales and derived ones:

```rust
pub struct Scales {
    // Primary (user-configurable)
    pub length: f64,       // L* in km
    pub thickness: f64,    // S* in km
    pub time: f64,         // τ* in Myr
    pub density: f64,      // ρ* in kg/m³

    // Derived (computed)
    pub velocity: f64,          // v* = L*/τ*
    pub strain_rate: f64,       // ε̇* = 1/τ*
    pub viscosity: f64,         // η* = ρ*·g·τ*·S*
    pub stress: f64,            // σ* = η*·ε̇*
    pub argand: f64,            // diagnostic
}

impl Scales {
    pub fn from_primary(...) -> Self { ... }
    pub fn report(&self) -> String { ... }  // logs at solver startup
}
```

`Scales::from_primary(L=350, S=35, τ=30, ρ=3300)` reproduces §3. The struct is threaded into `TectonicsConfig` and made available to all solver modules. This phase adds no behavioural change.

### Phase B — Internal nondimensionalization and new fields

Convert solver modules one at a time to operate in nondimensional units, with conversion to/from physical units only at the I/O boundary (input config, output heightmap).

Order of migration:

1. **Advection of $S$ (§4.6)** — already parameter-free in non-dim form; mostly a no-op renaming.
2. **Stokes + GPE (§4.1, §4.5)** — `apply_stokes`, `compute_rhs`. Absorb `gravity_factor * (rho_c/rho_m) * ...` into a single $\mathrm{Ar}$ coefficient.
3. **Power-law rheology (§4.2)** — `compute_viscosity`. Replace independent `eta_min/max`, `strain_rate_min` by dimensionless `visc_contrast` and `strain_rate_floor_rel`.
4. **Yielding + plastic memory (§4.3)** — Replace `yield_stress` by $\mathrm{Bi}$, add $\mathrm{De}_p$ for healing. Raise healing rate default from $0$ to a non-zero value once $\mathrm{De}_p$ is calibrated.
5. **Basal drag (§4.4)** — Replace hard `S >= 0.3` threshold by a smooth $\tilde S^2$ scaling; replace `basal_friction` by $\mathrm{Br}$.
6. **Source terms (§4.7)** — Enforce the ordering constraint as soft validation at config load.
7. **Slab pull (§4.8)** — Replace monotonic growth by exponential decay; remove `max_plate_velocity` clamp.
8. **Mantle forcing (§4.9)** — Replace raw velocity amplitude by $\mathrm{Mf}$.
9. **Cratonic rigidity (§4.10)** — Split-mechanism design: reduce $K$ to $[3, 8]$, add plastic immunity template for craton cells. Validate that $\kappa(A)$ stays within budget.
10. **Age field (§4.11)** — Add the scalar field $A$, its advection, its event-driven resets. Wire into export.

Each migration is a single PR with a numerical equivalence test: same seed, same config translated to the new variables, produces byte-identical output within floating-point tolerance (items 1–8). Items 9 and 10 are additive and do not have a byte-equivalent reference — they are validated against the Phase E stochastic protocol.

### Phase C — Null-space-aware preconditioner

In `apply_jacobi` and `apply_ssor`, wrap the output $z$ with a mean-projection: $z \leftarrow z - (\bar z_x, \bar z_y)$. Cost: two reductions per application, $O(N)$. Expected gain: $\sim 20$–$30\%$ reduction in BiCGSTAB iterations on realistic tectonic fields.

### Phase D — Diagnostics

At every tectonic step, log:

- Measured $\tilde\eta_\text{max}/\tilde\eta_\text{min}$ (warn if it exceeds the target contrast).
- Measured $\mathrm{Ar}$ (should stay near chosen value; drift indicates a scale error).
- Estimated $\kappa(A)$ via cheap power iteration, every $N$ steps.
- Measured $\tilde S_\text{eq}$ (mean thickness of active orogens) — expected near $1.8$–$2.0$; deviation flags $\mathrm{Ar}$ miscalibration.
- Fraction of cells in plastic regime.
- Distribution of boundary types on active boundaries (how many subduction vs collision vs rift cells).
- Distribution of cell ages (min, mean, max, histogram).
- Newton outcome classification — already logged, keep.

This gives a per-run fingerprint that makes "is this run in a healthy regime?" answerable at a glance.

### Phase E — Stochastic validation (per §3.4)

Implement the 20-seed × 3-duration sampling protocol as a test harness:

```rust
#[test]
#[ignore = "stochastic validation, run on demand"]
fn stochastic_validation_preset_a() { ... }
```

Run on demand (not in CI — too expensive). Produces a report (JSON + plots) with distributions of Newton outcomes, $\kappa(A)$, peak $\tilde S$, boundary type diversity, age field statistics, etc. Regression is a *shift in the distribution*, not a change in any single run.

Individual smoke tests (single seed, single config, short $T$) remain in the CI suite to catch obvious breakage.

This phase also validates (or refutes) the presets of §5.3 and the split-rigidity design of §4.10. A preset that fails its stochastic validation gets removed or revised.

### Phase F — Performance targets and multigrid preconditioner

Immediate targets:

- $64^2$: < 0.1 s per tectonic step, < 30 s for a 300-step run.
- $128^2$: < 1 s per step, < 5 min for a full simulation at $T = 6\,\tau^\*$.
- $256^2$: < 10 s per step.
- $512^2$: < 1 min per step, < 4 h for full simulation at $T = 6\,\tau^\*$.

If Phases C and D bring $128^2$ into the $< 1$ s/step target, the primary $128^2$ workflow is comfortable. **$256^2$ and $512^2$ will require a multigrid preconditioner** to reach acceptable times — point Jacobi and SSOR, even with null-space projection, scale poorly beyond $\sim 10^4$ unknowns.

Multigrid for variable-viscosity Stokes on a staggered grid is a known hard problem — a pointwise smoother does not converge reliably in the presence of strong viscosity contrasts. The accepted solution is a **block-staggered smoother** (Vanka-type or incomplete-LU variants) applied at each grid level, with geometric coarsening and matrix-free restriction/prolongation operators. Good references: Trottenberg et al. (2001) chapter on mixed problems, Kaus et al. (2016) for the tectonic-specific implementation.

**Do not attempt multigrid before Phases A–E are complete.** The current high $\kappa(A)$ is primarily a consequence of the unit mismatch described in §2.2, not a fundamental inadequacy of point preconditioners. Fixing the scaling first reduces the problem to a size where simple preconditioners are adequate for $128^2$, which is enough to validate the scaling approach. Multigrid becomes the *next* project once the scaled solver is working and profiled.

---

## 7. Open questions

**Q1. Is `s_max` still needed, and if so at what value?** The current hard clamp at `s_max = 2.5` is a safety net. Once $\mathrm{Ar}$ is calibrated, the physical equilibrium thickness $\tilde S_\text{eq}$ is an output of the solver — probably near $1.8$–$2.0$, but this is to be measured. The right question post-calibration is not "what value for `s_max`" but "does `s_max` still fire in practice":

- If post-scaling simulations show $\tilde S_\text{eq} \approx 1.8$ and peaks around $2.1$–$2.2$, then `s_max = 2.5` with smooth saturation is a fine safety net that rarely activates. Keep it.
- If peaks routinely reach or exceed $2.5$, then either $\mathrm{Ar}$ is mis-calibrated (the GPE-driven equilibrium is not happening), or the boundary sources are producing unphysical thickening faster than GPE can respond. Either diagnosis points to a scaling problem, not a clamp problem.

**Defer this decision to post-Phase-D**, when we have $\tilde S_\text{eq}$ diagnostics running. If needed, lower `s_max` to $2.2$ and use smooth saturation; if $\tilde S_\text{eq}$ sits at $1.8$–$2.0$ cleanly, leave at $2.5$.

**Q2. Cratonic rigidity via viscosity or plastic memory?** Decided: primary mechanism is plastic memory immunity (§4.10), with modest $K \in [3, 8]$ as secondary stiffening. This preserves the $\kappa(A)$ budget. Empirical observation supports this — the current $K = 5$ default is already near the practical upper limit for Newton convergence, and pushing higher breaks the solver without giving more realistic cratons.

**Q3. Thermal coupling — worth the cost?** Recommendation: implement (a) static temperature template first, as a simple spatial modulator of $\eta_0$ and $\tau_Y$. This is almost free and captures "cratons are cold and rigid, mobile belts are warm and weak" with no additional nonlinearity. If (a) plus the plastic memory mechanism of §4.10 give satisfactory craton behaviour, (b) and (c) are unnecessary — and that is the likely outcome.

Retain (b) and (c) in the backlog for the case where cratonic stability is insufficient even with (a) + plastic immunity:
- (b) = "interaction heats" feedback via a thermal-like field coupled to plastic strain. Adds a scalar field and a new nondimensional number. Moderate cost.
- (c) = full thermal coupling with advection, diffusion, Arrhenius viscosity, shear heating. Significant project, separate note.

**Q4. Do we keep both Picard and JFNK?** Keep both during Phase B migration — Picard serves as a robust fallback if Newton fails on specific seeds. After Phase E stochastic validation, if Newton converges reliably across the seed sample, Picard becomes a debug-only path. If Newton fails on > 5% of seeds even post-scaling, Picard remains the default and JFNK is the fast-path for easy cases.

**Q5. Healing rate coupled to crustal thickness?** Plausible as a future refinement — physically, thicker crust has a hotter base, so healing should be faster there. Implementation cost is low (a spatial modulation of $r_\text{healing}$). The concern is not complexity but *model rigidity*: adding couplings creates hard-to-debug interactions. Recommendation: implement as an opt-in feature (default off), validate that it does not destabilise the solver or narrow the range of producible continents, then enable by default if it improves realism without cost.

**Q6. Multigrid preconditioner for $256^2$–$512^2$.** Deferred to Phase F+ as explained in §6. Confirmed as the strategic direction for production-grade performance at larger grids. Implementation will require a dedicated design note covering smoother choice (Vanka vs ILU vs block Jacobi), coarsening strategy (geometric vs algebraic), and handling of the null space across levels.

**Q7. Erosion-induced feedback on tectonics.** The current pipeline runs tectonics to completion then erosion. Real tectonics has feedback: erosion removes topographic load, which lets GPE push up more material from depth, sustaining orogeny. Ymir does not model this — the resulting underestimation of orogenic lifetime is compensated implicitly by the choice of $\tau^\*$. Flagged for future "coupled pipeline" exploration, out of scope here.

---

## 8. References

- England, P., & McKenzie, D. (1982). A thin viscous sheet model for continental deformation. *Geophys. J. R. Astron. Soc.*, 70(2), 295–321. — Canonical source for Argand number and the sheet model.
- Houseman, G., & England, P. (1986). Finite strain calculations of continental deformation. *J. Geophys. Res.*, 91(B3), 3651–3663. — Nondimensional form and solution method.
- Turcotte, D. L., & Schubert, G. (2002). *Geodynamics* (2nd ed.). Cambridge University Press. — Physical scales, rheology, isostasy.
- Gerya, T. (2010). *Introduction to Numerical Geodynamics Modelling*. Cambridge University Press. — Staggered grid Stokes, nonlinear solvers, pedagogical reference.
- Kaus, B. J. P., Popov, A. A., et al. (2016). Forward and inverse modelling of lithospheric deformation on geological timescales. *NIC Symposium*. — Matrix-free Stokes with variable viscosity, preconditioning strategies, multigrid for staggered-grid Stokes.
- Moresi, L., & Solomatov, V. S. (1995). Numerical investigation of 2D convection with extremely large viscosity variations. *Phys. Fluids*, 7(9), 2154–2162. — Handling of $10^6+$ viscosity contrasts in Stokes solvers.
- Trottenberg, U., Oosterlee, C., & Schüller, A. (2001). *Multigrid*. Academic Press. — Reference for Phase F multigrid preconditioner, chapter on mixed problems.

---

*Working document. §5.1 Ar target range corrected 2026-04-21 after discovery of inconsistency with §3.2 thin-sheet scales (see §5.4 for the derivation and its implications). §5.3 presets remaining values will be refined empirically through Phase E. §6 migration order may be adjusted based on practical coupling between modules. §7 Q1 and Q3 explicitly wait for Phase D diagnostics before being decided.*
