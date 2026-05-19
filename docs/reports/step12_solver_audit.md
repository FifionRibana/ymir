# Step 12 R5b D0 — Solver Audit

**Date** : 2026-05-14
**Branche** : `112-step-12-interleaved-tectonic-erosion-workflow`
**Contexte** : R5.0 + R5.0.1 ont confirmé que CG cap = 2000 est atteint sur 100 % des configs dt testées (de 0.06 à 0.0006, factor 100) au régime mantle nominal mf=1.0. Kappa = 1.3-2.1 × 10⁴. Decision tree node 3 — diagnostic structurel solveur AVANT toute calibration.

Audit en lecture pure, pas de modification de code. Findings agrégés depuis 8 fichiers source + recherche git history.

## A — Architecture AMG actuelle

| Aspect | Valeur | Source |
|---|---|---|
| **Smoother** | Red-Black Gauss-Seidel (RBGS), 1 fwd + 1 bwd sweep / level | `amg/smoother.rs` |
| **Pre/post sweeps** | 1 / 1 (configurable) | `amg/mod.rs` |
| **Strong connections θ** | 0.25 (Ruge-Stüben classique) | `amg/strong_connections.rs:33-80` |
| **Cycle type** | V-cycle uniquement (par défaut) | `amg/vcycle.rs:28-92` |
| **FMG disponible** | Oui (`amg/fmg.rs:28-86`) mais **non utilisé** dans `AmgPreconditioner::apply` | `amg/mod.rs:210-211` |
| **Mode** | Preconditioner-only — PCG avec AMG comme M⁻¹ | `amg/mod.rs:148-216` |
| **Coupling blocks u-v / v-u** | **Ignorés** dans le préconditionneur (Option B'), conservés dans matvec | `amg/mod.rs:148-216` |
| **Coarse min unknowns** | 50 | `amg/mod.rs:106` |
| **Max hierarchy depth** | 7 levels | `amg/mod.rs` |
| **Coarse solve** | Dense LU (Doolittle, partial pivoting) | `amg/coarse_solve.rs` |
| **Benchmarks 32²/64²/128²** | **AUCUN** dédié AMG dans testsuite | `benches/`, `tests/` |

**Default actuel** : `LinearSolverConfig::default() = JacobiCG` (`stokes/solver.rs:226-230`). AMG est **opt-in** via `LinearSolverConfig::AmgCG(AmgConfig::default())`.

## B — Critères convergence Newton/CG actuels (tectonics_v2)

| Critère | Statut | Détail |
|---|---|---|
| **Résidu** | PRÉSENT | `resid ≤ abs_tol_eff` OR `resid ≤ rel_tol · r0` (`nonlinear_solver.rs:556-566`) |
| **abs_tol_eff** | PRÉSENT | `abs_tol.max(10 · linear_tol)` floor pour éviter Newton de chasser inatteignable (`:333`) |
| **Stall (résidu)** | PRÉSENT | `<1 % reduction × 3 iter` quand `far_from_convergence` (`:375-387`) |
| **Critère état (`\|Δu\|/\|u\|`)** | **ABSENT** | Aucun champ `state_tolerance` dans `NewtonConfig` |
| **Détection oscillation (cosinus)** | **ABSENT** | Aucune piste de Newton steps dans la boucle |
| **Outcome `Oscillation`** | **ABSENT** | Variants : `Converged`, `Stalled`, `Diverged`, `CappedIters` (`:85-107`) |
| **Line search** | PRÉSENT | Armijo c₁=1e-4, max 10 backtracks (`:498-543`) |
| **Adaptive dt sur stall** | **ABSENT** | `AdaptiveDtConfig` (Step 4) existe mais **non liée** au statut Newton |

**Defaults** : `rel_tol = 1e-6`, `abs_tol = 1e-10`, `max_outer_iters = 20`, `linear_max_iter = 2000`.

## C — Historique critères pré-milestone

L'utilisateur a mentionné une version pré-milestone avec critère état + oscillation + dt réduction. **Findings** :

`tectonics_v2/` est un **clean rewrite** (Step 0 commit `a8c8f3a`), pas une migration de `tectonics/`. Les critères hybrides #49 existaient dans la **branche `49-hybrid-convergence-criterion-...`** mais **jamais mergés** dans `main` ni dans `112-step-12-...`.

**Commits identifiés** (branche #49 archivée) :
- `654a8e1` — FEAT : Add state_tolerance and oscillation detection config fields
- `e4a7ddb` — FEAT : State-based convergence criterion in Newton solve
- `45bdfc0` — FEAT : Oscillation detection via consecutive step alignment
- `7bf23b0` — TEST : near-tolerance stagnation convergence test

### Snippet — Convergence sur l'état (e4a7ddb, `tectonics/solver/newton.rs:375-395` pré-rewrite)

```rust
let actual_step: Vec<f64> =
    ws.jfnk_delta_v.iter().map(|x| final_alpha * x).collect();

if k >= newton_config.min_iterations_before_classification
    && residual_history.len() > newton_config.trend_window
{
    let v_state_norm: f64 =
        v_old.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-30);
    let step_norm: f64 =
        actual_step.iter().map(|x| x * x).sum::<f64>().sqrt();
    let relative_step = step_norm / v_state_norm;

    if relative_step < newton_config.state_tolerance && trend_descending {
        return NewtonResult {
            outcome: NewtonOutcome::ConvergedOnState,
            ...
        };
    }
}
```

### Snippet — Détection oscillation (45bdfc0, `tectonics/solver/newton.rs:365-401` pré-rewrite)

```rust
let dot: f64 =
    actual_step.iter().zip(prev.iter()).map(|(a, b)| a * b).sum();
let denom = n_curr * n_prev;
if denom > 1e-30 {
    let cos_theta = dot / denom;
    if cos_theta < newton_config.oscillation_cosine_threshold {
        consecutive_anti_aligned += 1;
        if consecutive_anti_aligned >= 2 {
            return NewtonResult {
                outcome: NewtonOutcome::Oscillation,
                ...
            };
        }
    }
}
```

### Paramètres cibles (snapshot pré-rewrite #49)

- `state_tolerance: 1e-4` (relative velocity increment floor)
- `trend_window: 3` (iter window pour résidu descendante)
- `oscillation_cosine_threshold: -0.5` (angle anti-aligné)
- `min_iterations_before_classification: 3` (gate anti-bruit early)

## D — Gaps identifiés (vs ce que user mentionne)

1. **State-based convergence absente** dans `tectonics_v2`. Le solveur chasse un résidu plateau alors que l'état physique est stable → `Stalled` qui n'aurait pas dû se produire si on regardait `|Δu|/|u|`.

2. **Oscillation detection absente**. Le solveur brûle ses 20 outer-iters sur des steps anti-alignés sans sortir explicitement.

3. **AMG opt-in et NON activé par défaut**. 2812 lignes d'implémentation AMG complète mais default = JacobiCG saturé.

4. **Adaptive dt non lié au solver state**. `AdaptiveDtConfig` existe (Step 4) mais ne consomme pas `NonlinearOutcome::Stalled/Oscillating` pour réduire dt.

5. **V-cycle seul** dans `AmgPreconditioner::apply`. FMG implémenté (~2× speedup sur Poisson dans benchmarks unitaires) mais non utilisé.

6. **Strong connections θ = 0.25 fixe**, ne s'adapte pas aux régimes très ill-conditioned (κ ~ 10⁴).

7. **Couplage u-v ignoré dans préconditionneur** (Option B') — peut être insuffisant pour les régimes où couplage Stokes domine.

## E — Recommandations classées par effort

### Effort faible (< 1 h)

- **R5b D1** : test `LinearSolverConfig::AmgCG(AmgConfig::default())` vs `JacobiCG` sur 10 steps × 32² × mf=1.0 workflow OFF, contrôler CG iter. Cf. brief — non-spéculatif, lecture directe des métriques.
- **Logs détaillés Newton/CG** dans le harness (résidu per-iter, alpha, CG iters) pour identifier où le solver stalle exactement.

### Effort moyen (1–4 h)

- **R5b D2 — Réintégrer state-based convergence** : porter `e4a7ddb` dans `tectonics_v2/stokes/nonlinear_solver.rs`. Ajouter `state_tolerance`, `trend_window`, `min_iterations_before_classification` à `NewtonConfig`. Implémenter le critère ConvergedOnState + flag `trend_descending`.
- **R5b D2 — Réintégrer oscillation detection** : porter `45bdfc0`. Tracker `prev_delta_v` + cosine. Ajouter variant `NonlinearOutcome::Oscillating { outer_iters }`.
- **Wire adaptive dt** : faire que `Oscillating` ou `Stalled` chez `nonlinear_solver` déclenche `dt /= 2` chez le harness.

### Effort élevé (> 4 h)

- **R5b D3 — Reformulation Stokes** (conditionnel D1+D2 KO) : AMG point-based 2×2 blocks (couplage explicite), SA-AMG (smoothed aggregation), ou refonte physique velocity-pressure. Pas démarrer sans confirmation user explicite.

## F — Finding rétroactif important

> **Step 8 calibration solveur s'est faite avec régime mantle OFF ou faible.** Le bug "CG cap saturé sur mf=1.0 nominal" a été masqué jusqu'à Step 12 workflow ON, qui multiplie les cycles macro et expose la non-convergence. Step 11 standalone runtime (~50 min/64² rapporté par user) est probablement victime du même bug, moins visible parce que pas multiplié par 5 cycles. **Step 8 + Step 11 validation à reprendre rétroactivement** après fix solveur.

## Conclusion D0

Le rewrite `tectonics_v2` (Step 0) a **délibérément simplifié** le solveur Newton — résidu seulement, pas d'état, pas d'oscillation. Cohérent avec la philosophie de la milestone (rebâtir le cœur Stokes propre avant d'empiler des features). Mais cela a **perdu les correctifs hybrides #49** qui auraient évité la saturation actuelle.

**Le fix n'est pas dans AMG par défaut seul** (proba qu'il désature à 32² faible — note user "AMG inefficace à 32²/64²"). Le fix probable est **D2 (criteres pré-#49)** réintégré dans `tectonics_v2`, possiblement combiné avec D1 (AMG activé sur workflow ON régime).

Plan D1 + D2 reste dans la fourchette d'effort raisonnable (< 5h). D3 (reformulation Stokes) en dernier recours.
