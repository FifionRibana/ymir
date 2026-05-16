# Step 12 R6.1 — Mantle evolution_rate physics decision

**Branche** : `112-step-12-interleaved-tectonic-erosion-workflow`
**Statut** : décision physique validée par l'utilisateur. Implémentation R6.2 enchaînée.
**Référence amont** : finding R5b evolution_rate sweep (3 configs bit-identical) → `MantleConfig::Enabled.evolution_rate` exposé mais non câblé dans le solveur depuis Step 8 (D6, "Out of scope").

## 1. Contexte

Le sweep `mf ∈ {0.5, 0.6, 0.7, 0.8, 0.9}` à 32² × 5 cycles × 20 steps a démontré un **régime binaire** :

| mf | preservation cratons | dynamique soutenue (peak |v| > 0.1) |
|---|---|---|
| 0.5 | PASS (98 %) | FAIL (peak |v| ≈ 1.2e-3) |
| 0.7 | PASS (84 %) | FAIL (1 cycle / 5) |
| 0.8 | FAIL — transition pathologique |
| 0.9, 1.0 | FAIL (dissolution) | PASS partiel |

Aucune valeur de mf (seule) ne produit simultanément **continents préservés + dynamique tectonique soutenue**. Le pattern mantle statique tire le système vers un fixed point : soit le bias est trop faible (système quiescent post-équilibre) soit trop fort (continents dissous).

**Hypothèse** : un pattern mantle **time-varying** maintient le flux tectonique en permanence en redéplaçant les zones de convergence/divergence, débloquant la fenêtre multi-dim (mf, evolution_rate).

## 2. Décision : Phys.A — Phase drift

### Mécanisme

Les modes Fourier du stream function `ψ = Σ_k a_k · sin(kx·2π·x + φx_k) · sin(ky·2π·y + φy_k)` voient leurs phases dériver linéairement avec le temps :

```
φx_k(t) = φx_k(0) + ω · t
φy_k(t) = φy_k(0) + ω · t
```

avec `ω = evolution_rate · TAU` (unité : une rotation complète en phase espace par temps non-dim quand `evolution_rate = 1`).

### Rationale (4 arguments)

**1. Div-freeness préservée exactement.**

La cancellation algébrique des 4 termes nodaux dans `pattern.rs:37-43` (`(ψ[i+1,j+1] − ψ[i+1,j] − ψ[i,j+1] + ψ[i,j])` etc.) ne dépend que de la propriété structurelle "ψ stocké aux nodes du grid". Décaler les phases change les **valeurs** des sines mais préserve `div(curl(ψ)) ≡ 0`. Acceptance Step 8 critique `div_v_mantle_max < 1e-10` survit à tout `(φx, φy)`.

**2. Implémentation minimale.**

Seul `draw_modes(...)` est étendu par un `apply_phase_drift(t, evolution_rate)`. Le reste de la chaîne (`sample_nodal`, normalisation, staggered curl `pattern.rs:107-133`, MAC layout, force assembly `mantle_force.rs`) reste intact. Surface : ~50 LOC production + 4 tests régression.

**3. Hypothèse falsifiable.**

Après `t·ω = π/2` (≈ 5–6 steps à 32² avec dt non-dim ≈ 0.02 et evolution_rate=0.1), les sines sont quasi-orthogonales aux originales → vraie évolution mesurable. Si Phys.A débloque le sweet spot, c'est la preuve formelle que **mantle dynamique** est la pièce manquante du puzzle Living Landz. Si Phys.A échoue R6.3, signal négatif clair pour escalader vers Phys.C.

**4. Pas de scope creep.**

- Phys.B (amplitudes a_k(t)) demande 6 fréquences ω_k indépendantes + phases → 12 paramètres au lieu d'1. Spec physique à dérouler.
- Phys.C (advection `ψ_0(x-vt, y)`) demande direction de drift + magnitude + gestion CFL. Plus physique mais plus lourd.

Aucun des deux ne se justifie tant que Phys.A n'est pas falsifié.

### Caveat plumes lifecycle accepté avec nuance

Phys.A est une **rotation rigide en espace de phase**, pas une naissance/mort de modes. Les zones de convergence/divergence se **déplacent** dans l'espace mais ne se créent ni ne disparaissent. Conséquence physique :

Sur Terre, les chaînes orogéniques se forment en 50+ Ma. Avec `dt = 1.8 Ma/step × 100 steps = 180 Ma simulés` et `evolution_rate=0.10` → rotation phase `0.6 · TAU = 216°` cumulée. Les zones migrent **lentement** : les chaînes ont le temps de se former avant migration substantielle. C'est l'argument empirique d'utilisateur, à valider sur R6.3 visuel.

Si R6.3 montre que les zones migrent trop vite (chaînes pas le temps de se former), `evolution_rate` doit baisser ou Phys.C devient nécessaire.

## 3. Convention de mise à l'échelle

```
ω = evolution_rate · TAU
phase_offset(t_nondim) = ω · t_nondim = evolution_rate · TAU · t_nondim
```

Calibration cible :

| `evolution_rate` | Cumul phase sur 180 Ma (= 100 steps × 1.8 Ma) | Régime visuel attendu |
|---|---|---|
| 0.0  | 0 (statique) | Step 8 baseline — quiescence après équilibre |
| 0.05 | 0.3 · TAU ≈ 108° | Lente migration, chaînes faciles à former |
| 0.10 | 0.6 · TAU ≈ 216° | Migration modérée — sweet spot probable |
| 0.20 | 1.2 · TAU ≈ 432° | Plus d'un tour — chaînes éphémères |
| 1.0  | 6 · TAU | Pattern méconnaissable d'un step à l'autre |

Cette mise à l'échelle est **dimensionnelle**. La vérification empirique en R6.3 trouvera la fenêtre utile pour Living Landz.

## 4. Caveats actés en R6.1

### Normalisation post-drift : **alternative t=0 frozen**

Décision (déléguée par utilisateur, tranchée R6.2) : **figer la normalisation à `init_norm = max|ψ(t=0)|`** plutôt que renormaliser à chaque step.

Raisons :
1. Préférence utilisateur ("dynamique sinusoïdale plus propre")
2. Évite jitter sur la position de l'argmax step-à-step
3. Coût `max|ψ|` économisé à chaque step (négligeable, mais propre)
4. Cohérent avec le découplage modes ↔ rendu : les modes sont figés à init via `seed`, leur normalisation l'est aussi

Conséquence : l'amplitude pic `max|ψ(t)|` peut fluctuer ±20 % autour de 1.0 selon (φx, φy). Le `peak_v_mantle_pattern` reporté dans `newton_agg` sera l'amplitude **à t=0** (cohérent avec le contrat actuel) — mais le pattern effectif à `t>0` aura une amplitude légèrement différente.

**Test R6.2 acté** : `evolution_rate_pattern_amplitude_bounded` — vérifier que `max|ψ(t)| ∈ [0.5, 1.5]` sur 100 steps (encadrement raisonnable).

### Pattern reconstruit en place, pas réalloué

`MantlePattern::rebuild_from_psi(...)` est ajouté pour réécrire les buffers existants sans `Field2D::new(nx, ny)` à chaque step. Coût alloc évité (~32 KB par step à 64²).

## 5. Hypothèse falsifiable R6.3

**Statement** : Si une config `(mf ∈ {0.5, 0.7, 0.8, 1.0}, evolution_rate ∈ {0.05, 0.10, 0.20})` avec Phys.A active produit les **6 critères** R4.1–R4.6 (continents émergés + cratons préservés + bordures déformées + conservation + drainage actif + dynamique soutenue), alors **Phys.A suffit pour Living Landz MVP** et Step 12 est fermable.

**Falsifiable** : si **aucune** des 12 configs (16 − 4 déjà fait) ne passe les 6 critères, preuve formelle que **phase drift est inadéquat** comme mécanisme d'évolution → pivot **Phys.C** (stream function advection) en Step 12.X. Pas de tergiversation, pas de tuning marginal pour faire passer 5/6 critères.

## 6. Alternative Phys.C planifiée

Si R6.3 verdict EVO.D (aucune config ne passe), Step 12.X ouvert avec spec Phys.C :

```
ψ(x, y, t) = ψ_0(x - v_drift_x · t, y - v_drift_y · t)
```

Implémentation : décaler les arguments `x → (x - v_drift_x · t) mod 1`, `y → (y - v_drift_y · t) mod 1` dans `sample_nodal`. Modes inchangés, phases inchangées. Plumes/slabs persistent et défilent.

**Paramètres additionnels** : `(v_drift_x, v_drift_y)` ∈ `[-1.0, 1.0]²` non-dim. Convention sweep à définir au moment Step 12.X.

## 7. Plan implémentation R6.2

**Fichiers modifiés** :

1. `crates/ymir-core/src/tectonics_v2/mantle/stream_function.rs`
   - Ajouter `StreamFunctionBuilder { base_modes, init_norm }` avec `new(nx, ny, config)` + `sample_at_time(nx, ny, t_nondim, evolution_rate)`
   - Garder `generate_stream_function` comme wrapper appellant `builder.sample_at_time(0.0, 0.0)` → bit-identical avec pré-R6
   - Helper privé `apply_phase_drift(modes, t_nondim, evolution_rate) -> Vec<Mode>`

2. `crates/ymir-core/src/tectonics_v2/mantle/pattern.rs`
   - Ajouter `MantlePattern::rebuild_from_psi(&mut self, psi, dx, dy, idx_x, idx_y)` pour reconstruction in-place sans alloc

3. `crates/ymir-core/src/tectonics_v2/diagnostics/harness.rs`
   - Stocker `mantle_psi_builder: Option<StreamFunctionBuilder>` alongside `mantle_pattern`
   - Marquer `mantle_pattern` comme `mut`
   - Insérer rebuild juste avant `MantleForce::new(...).accumulate(...)` ligne 1332 :
     ```rust
     if evolution_rate > 0.0 {
         let t_nondim = step as f64 * dt_target;
         let psi = builder.sample_at_time(nx, ny, t_nondim, evolution_rate);
         pattern.rebuild_from_psi(&psi, dx, dy, &idx_x, &idx_y);
     }
     ```
   - Mettre à jour le commentaire obsolète ligne 977 ("static at Step 8 (evolution_rate = 0, D6)")

**Tests régression ajoutés** (`crates/ymir-core/tests/v2_mantle_evolution_rate.rs`) :

| Test | Vérifie |
|---|---|
| `evolution_rate_zero_pattern_constant_multistep` | `evolution_rate=0` → ψ(t=0) ≡ ψ(t=10·dt) bit-identical |
| `generate_stream_function_unchanged` | wrapper public `generate_stream_function` produit même output pre/post-refactor (seed=42, multi-tailles) |
| `evolution_rate_nonzero_evolves_measurably` | `evolution_rate=0.1` → `||ψ(t=10·dt) - ψ(t=0)||_inf > 0.1` (vraie évolution) |
| `evolution_rate_preserves_div_freeness_multistep` | `evolution_rate=0.1` → `div_v_mantle_max < 1e-10` à chaque step pour 20 steps |

## 8. Anti-patterns écartés

- ❌ Renormaliser `max|ψ|` runtime à chaque step → jitter argmax
- ❌ Modifier `MantleConfig::Enabled` wire-format → break sérialisation V2Spec / preset JSON
- ❌ Implémenter Phys.A en hardcoded values → param `evolution_rate` existe déjà, le câbler suffit
- ❌ Skip tests régression Disabled bit-identical → tous les tests pré-R6 doivent rester verts
- ❌ Tuner `evolution_rate` ad-hoc pour faire passer R6.3 → si 6/6 critères impossibles, EVO.D et remontée
- ❌ Étendre vers Step 13.6 (bordures angulaires) dans R6 → scope séparé

## 9. Critère de succès R6.2

Tous les 4 tests régression verts **et** suite `cargo test -p ymir-core` reste verte (tests pré-R6 pas cassés). Mesure empirique `||ψ(t=10) - ψ(t=0)||_inf > 0.1` confirmée. Div-freeness multi-step `< 1e-10` confirmée. **Pause obligatoire** ensuite : validation utilisateur avant lancer R6.3 sweep.

## 10. Caveat structurel découvert pendant R6.2 — period field = π (pas TAU)

Pendant l'écriture du test `evolution_rate_nonzero_evolves_measurably`, j'ai pris `evolution_rate = 0.5, t = 1.0` espérant `phase_offset = π` produirait `ψ(t=1) = −ψ(t=0)` (chaque sin flippe). Test FAIL : `||·||_inf = 1.9e-15` au lieu de `~2.0`. Cause structurelle de la formulation Phys.A.

### Cause

La forme produit `sin(arg_x + φ_x) · sin(arg_y + φ_y)` est invariante au shift simultané de `φ_x` et `φ_y` par `π` :

```
sin(arg_x + φ_x + π) · sin(arg_y + φ_y + π)
= (−sin(arg_x + φ_x)) · (−sin(arg_y + φ_y))
= sin(arg_x + φ_x) · sin(arg_y + φ_y)
```

Les deux signes se compensent. Conséquence : le **champ ψ a une période effective de π en phase**, pas TAU comme je l'avais supposé.

### Impact sur la convention de mise à l'échelle

Recalcul correct (avec `cfg.total_time_nondim ≈ 2.0` pour les presets active_medley / single_continent, dt_nondim ≈ 0.02 sur 100 steps) :

| `evolution_rate` | Cumul `phase_offset = evo · TAU · t_end` | En units field-period (π) | Régime |
|---|---|---|---|
| 0.0  | 0 | 0 | statique (Step 8 baseline) |
| 0.05 | 0.05·TAU·2 = 0.1·TAU ≈ 36°  | 0.2·π | très lente — moins d'1/5 de période |
| 0.10 | 0.10·TAU·2 = 0.2·TAU ≈ 72°  | 0.4·π | lente — pas de retour identité |
| 0.20 | 0.20·TAU·2 = 0.4·TAU ≈ 144° | 0.8·π | modérée — proche de l'identité-π |
| 0.50 | 0.50·TAU·2 = 1.0·TAU = 360°  | 2·π    | passe par 1 identité-π mid-run |
| 1.0  | 1.0·TAU·2 = 2·TAU = 720°    | 4·π    | 2 identités-π pendant le run |

**Conséquence pour R6.3 sweep** : la fenêtre `evolution_rate ∈ {0.05, 0.10, 0.20}` initialement planifiée stays dans `[0, π]` côté phase cumulée. Aucune `evolution_rate` ne fait revenir le pattern à `ψ(0)` pendant un run de 100 steps. **Bonne nouvelle** — pas d'effet "yo-yo" non-intentionnel sur l'évolution du pattern.

Si R6.3 EVO.D (aucune config ne passe) et qu'on veut tester un régime "yo-yo" volontaire pour cartographier le comportement avant pivot Phys.C, étendre le sweep à `evolution_rate ∈ {0.50, 1.0}` avant d'abandonner Phys.A.

### Mitigation alternative non retenue ici

Désymétriser le drift entre x et y : `φx += ω·t`, `φy += ω·t · sqrt(2)` (ou tout incommensurable). La compensation `(−1)·(−1) = 1` disparaît, la période effective redevient TAU, et la dynamique est plus riche. Pas implémenté en R6.2 — préserver la simplicité MVP. À retenir si Phys.A se révèle trop pauvre.

### Test correspondant

Le test 3 `evolution_rate_nonzero_evolves_measurably` :
- Probe primaire : `phase_offset = π/2` (quadrature) → `||·||_inf > 0.1` ✓
- Probe secondaire : `phase_offset = π` → `||·||_inf < 1e-10` (documente l'identité-π comme propriété structurelle, pas comme bug)

Les deux probes coexistent dans le même test pour ancrer le caveat dans la suite régression.
