# Step 12 R7.A.2.1 — Composite continental profile formula spec

**Branche** : `112-step-12-interleaved-tectonic-erosion-workflow`
**Tête** : `5cf5380` (R7.A.1 + altitude blur fix)
**Statut** : décisions formule documentées avant impl. R7.A.2.2 enchaîne sous validation.

## 1. Contexte amont

R7.A.1 a livré l'infrastructure `InitMode::Orogenic` mais le visuel init révèle que **profil unique par plate est insuffisant** :
- `RadialPeak` seul → continent "blob" homogène (dôme uniforme)
- `Orogenic` seul → "plaine + spikes ponctuels" (crête 3 cellules à 64²)

Pattern méthodologique R7.A.1 documenté en mémoire ([feedback_recursive_tuning_signals_structural.md](C:\Users\obruneau\.claude\projects\d--Personnel-Project-ymir\memory\feedback_recursive_tuning_signals_structural.md)) : 5 tweaks paramétriques consécutifs cette session → cause root structurelle, pas paramétrique.

R7.A.2 valide l'hypothèse **composition additive** : un dôme général + une crête superposée reproduisent la morphologie type Amérique du Sud (continent élevé + Andes côte ouest). C'est la dernière approche "intra-plate" avant pivot R7.B (Voronoï hiérarchique).

## 2. Architecture : variant sibling `InitMode::Composite`

Cohérent avec la stratégie R7.A.1 : un nouveau variant `InitMode::Composite { … }` à côté de `RadialProfile`, `RadialProfileWithFBM`, `Orogenic`. Pas de refactor de l'enum existant.

```rust
pub enum InitMode {
    Checkerboard,
    Uniform { ... },
    Gaussian { ... },
    Convolution { ... },
    RadialProfile { ... },
    RadialProfileWithFBM { ... },
    Orogenic { ... },           // R7.A.1
    Composite { ... },          // R7.A.2 (nouveau)
}
```

`RadialProfile`, `Orogenic`, et tous les autres variants restent **bit-identical** par construction (signatures inchangées, dispatch dans `init_s_field` intact, modules `radial_profile.rs` et `orogenic_profile.rs` inchangés).

### Structure de paramètres : nested structs

Le user avait esquissé `RadialPeakParams` / `OrogenicParams` — ces types n'existent pas encore. Trois options :

| Approche | Avantage | Inconvénient |
|---|---|---|
| Inline flat (~10 champs dans Composite) | Pas de nouveau type | Verbeux, JSON plat encombré |
| Nested structs | Lisible, JSON clean | 2-3 nouveaux types |
| Réutiliser InitMode::RadialProfile fields | Pas de duplication | Couplage variant ↔ struct |

**Décision** : **nested structs**. Définit `CompositeRadialParams` et `CompositeOrogenicRidgeParams` dans `composite_profile.rs`. Surface +2 structs, mais JSON et code restent lisibles.

```rust
pub struct CompositeRadialParams {
    pub continental_value: f64,
    pub profile_shape: ProfileShape,
}

pub struct CompositeOrogenicRidgeParams {
    pub peak_value: f64,
    pub base_continental_value: f64,
    pub half_length_ratio: f64,
    pub width_sigma_ratio: f64,
    pub orientation: OrogenicOrientation,
    pub offset_along_axis_ratio: f64,   // R7.A.2 — Q.f.3
}

pub enum CompositeCap {
    UsePeakOrogenic,         // default — Q.f.2
    Fixed { value: f64 },
}

InitMode::Composite {
    radial: CompositeRadialParams,
    orogenic_ridge: CompositeOrogenicRidgeParams,
    oceanic_value: f64,
    cap: CompositeCap,
}
```

## 3. Formule composite (Q.f.1 + Q.f.2 + Q.f.3 + Q.f.4)

### Q.f.1 — Formule `ridge_only` : **option (a)** ✅

```
ridge_only(cell) = (peak_value - base_continental_value) · long_mod · trans_profile
```

Où :
- `long_mod = smoothstep((1 - |d_along_adjusted| / half_length).clamp(0, 1))`
- `trans_profile = exp(-(d_perp / width_sigma)²)`
- `d_along_adjusted = d_along - offset_along_axis_ratio · half_length` (offset Q.f.3)

**Pourquoi (a) vs (b)** :
- (a) borne par construction `∈ [0, peak - base]`, jamais négative
- (a) découplée de la valeur absolue `orogenic.sample(cell)` → indépendante de la convention de `base_continental_value` dans le contexte de l'orogenic seul vs composite
- (b) aurait besoin d'un clamp à 0 si `orogenic.sample < base`, ajoute une branche inutile

**Garanties** :
- `ridge_only ≥ 0` strictement
- `ridge_only ≤ peak_value - base_continental_value` (typiquement 0.35 avec defaults)
- Continu et C¹ partout (smoothstep + gaussienne, comme R7.A.1)
- À l'axe de crête et au centroïde de la plate : `ridge_only = peak - base = 0.35`
- Loin de l'axe (`|d_perp| > 3σ`) ou hors plateau (`|d_along| > half_length`) : `ridge_only → 0`

### Q.f.2 — Cap : `cap = peak_orogenic` ✅

```
cap_value = match cap {
    CompositeCap::UsePeakOrogenic => orogenic_ridge.peak_value,
    CompositeCap::Fixed { value } => *value,
}
```

**Trade-off** : avec `peak_radial = 0.95` et `ridge_only_max = 0.35`, la somme uncapped = 1.30. Le cap à 1.20 (peak_orogenic) coupe les 0.10 du sommet — petit écrêtage local au centre exact de la crête.

Alternative `cap = 1.30` (sum max) préserverait la totalité du ridge_only. Mais le user a tranché pour `cap = peak_orogenic` (cohérence : composite max = orogenic-only max, comportement prévisible).

Conséquence : visuel composite aura **même peak amplitude** que orogenic-only à 64² (~1.17 dans la pratique car centroïde discret), mais **élévation continentale beaucoup plus élevée autour** (dôme), donc contraste S̃ entre base continentale (0.85) et crête (1.20) reste 0.35 — comme orogenic-only — **MAIS** le contraste vs océan (0.20) sur les piémonts devient 0.95 au lieu de 0.85. Le rendu altitude doit être plus net car la dynamique S̃ visible s'étire.

### Q.f.3 — Offset crête vs centre plate : **(ii) configurable, default 0** ✅

`offset_along_axis_ratio: f64` (default `0.0`) déplace la crête le long de l'axe PCA depuis le centroïde :

```
d_along_adjusted = d_along - offset_along_axis_ratio · half_length
```

- `0.0` : crête centrée sur le centroïde (R7.A.1 behavior)
- `+0.3` : crête déplacée de 30 % de half_length le long de l'axe positif
- `-0.5` : crête près du bord négatif de la plate

Pour MVP R7.A.2 : laisser à `0.0` (cohérent avec R7.A.1, valide formule). Si Andes-côte-ouest visuellement nécessaire : R7.A.2.bis ajustement vers `0.3-0.5`.

**Caveat** : `offset_along_axis_ratio > 0.5` met la crête au-delà de `half_length` (crête sort du plateau longitudinal). Pas blocant numériquement (`long_mod` se clamp à 0 toujours), mais visuellement la crête disparaît hors du plateau. À documenter dans l'API.

### Q.f.4 — Calibration radial vs orogenic ✅

Défauts cohérents avec l'analyse user :

| Param | Valeur | Origine |
|---|---|---|
| `radial.continental_value` | 0.95 | CONTINENTAL_VALUE_DEFAULT (Step 13) |
| `radial.profile_shape` | Smoothstep | défault Step 13 |
| `orogenic_ridge.peak_value` | 1.20 | OROGENIC_PEAK_VALUE_DEFAULT (R7.A.1) |
| `orogenic_ridge.base_continental_value` | 0.85 | OROGENIC_BASE_VALUE_DEFAULT (R7.A.1) |
| `orogenic_ridge.half_length_ratio` | 0.40 | défaut R7.A.1 |
| `orogenic_ridge.width_sigma_ratio` | 0.10 | **R7.A.1 σ=0.10 (Himalaya-like)** — variant validée empiriquement |
| `orogenic_ridge.orientation` | PlateMainAxisPca | défaut R7.A.1 |
| `orogenic_ridge.offset_along_axis_ratio` | 0.0 | Q.f.3 default |
| `oceanic_value` | 0.20 | OROGENIC_OCEANIC_VALUE_DEFAULT |
| `cap` | UsePeakOrogenic | Q.f.2 |

**S̃ aux points cardinaux d'une plate continentale typique (L_plate ≈ 12 cells à 64²)** :

| Position | radial | ridge_only | composite | capped |
|---|---|---|---|---|
| Boundary plate (t=0) | 0.20 | 0 | 0.20 | 0.20 |
| Mid-radius hors crête (t=0.5) | 0.58 | 0 | 0.58 | 0.58 |
| Centroïde hors crête | 0.95 | 0 | 0.95 | 0.95 |
| Centroïde axe crête | 0.95 | 0.35 | 1.30 | **1.20** (cap) |
| Centroïde 1σ off-axis | 0.95 | 0.13 | 1.08 | 1.08 |
| Centroïde 2σ off-axis | 0.95 | 0.006 | 0.96 | 0.96 |
| Pic crête à 0.7 half_length | 0.92 | 0.27 | 1.19 | 1.19 |
| Bout crête à half_length | 0.85 | 0 | 0.85 | 0.85 |

**Lecture visuelle attendue** :
- Continent entier élevé entre 0.20 (bord) et 0.95 (centroïde) — dôme général
- Crête superposée monte de 0.95 (centroïde hors axe) à 1.20 (axe), sur ~2-3 cellules transversales et ~5-6 cellules longitudinales
- Piémont continu entre 0.95 et 1.20 (pas de cliff)
- Bordures plate transition continue 0.20 → 0.95 (smoothstep RadialProfile, comme Step 13)

C'est exactement la signature « dôme + ridge » type Andes/Amérique du Sud.

## 4. Cellules océaniques

Comme R7.A.1, les cellules océaniques **bypassent** la formule composite et retournent uniformément `oceanic_value` (0.20 par défaut). Le composite ne s'applique qu'aux cellules `PlateType::Continental`.

## 5. Caveats acceptés

### 5.1 Le cap à 1.20 écrête le sommet

Au centre exact de la crête, `radial + ridge_only = 1.30` est coupé à 1.20. C'est ~7 % de la dynamique perdue au pic, mais le visuel reste lisible (le pic est local à 1-2 cellules, capper là n'aplatit pas le dôme général).

### 5.2 Pas d'élargissement de la crête vs Orogenic-seul

Le ridge_amount est le **même** que dans `Orogenic` mode. Composite ne rend pas la crête plus large — il ajoute un dôme **autour**. Donc :
- La crête à 64² σ=0.10 reste ~3 cellules FWHM
- MAIS elle est sur un plateau à 0.85-0.95 au lieu de plaine à 0.85 → la chaîne est moins "isolée" visuellement
- Après érosion 5 cycles, hypothèse de travail : le plateau dome préserve la signature visuelle même si le pic se dégrade

Si après simulation Run C la crête s'érode quand même (résolution discrete trop limitée), pivot R7.B.

### 5.3 Hot reload / régression bit-identical strict

Tous les modes pré-R7.A.2 (Checkerboard, Uniform, Gaussian, Convolution, RadialProfile, RadialProfileWithFBM, Orogenic) doivent rester **byte-equal** post-impl. Test à chaque commit.

## 6. Tests régression prévus pour R7.A.2.2

| # | Test | Vérifie |
|---|---|---|
| 1 | `composite_oceanic_uniform` | Oceanic cells retournent exactement `oceanic_value` |
| 2 | `composite_dome_visible_without_ridge` | Avec `peak = base` (ridge_only ≡ 0) → S̃ === RadialProfile output |
| 3 | `composite_ridge_visible_without_dome` | Avec `radial.continental = oceanic` (dôme = 0.20 partout) + ridge → S̃ ≈ 0.20 + ridge (capé) |
| 4 | `composite_cap_respected` | max(S̃) ≤ `cap_value` exact |
| 5 | `composite_determinism` | Même seed → byte-equal output |
| 6 | `composite_peak_at_centroid_on_axis` | Cellule au centroïde sur l'axe PCA → S̃ proche de cap (large grid) |

Pas de test "byte-equal pré-R7.A.2 pour RadialProfile/Orogenic" — c'est garanti par construction (signatures + modules existants intacts) et déjà couvert par les tests d'origine.

## 7. Hors scope explicite

- Pas de modification Voronoï (R7.B reserved)
- Pas de FBM superposé (R7.A.2.1 garde le composite simple ; éventuel R7.A.3 si justifié)
- Pas de multi-crête par plate (1 ridge max)
- Pas de cratonic specialisation (factor field reste calculé par BFS, indépendant du profil S̃)
- Pas de yielding feedback sur la crête (rheology inchangée)

## 8. Hypothèse falsifiable R7.A.2

**Statement** : Si Run C (Composite, 64² × mf=1.0 × evo=0.10 × craton_amp=3) passe **5/6 ou 6/6 critères R4.1-R4.6** (incluant visuel R4.3 chaînes formées + bordures déformées) ET domine Run A (Radial) et Run B (Orogenic-seul) sur les critères visuels, alors **composition additive suffit pour MVP Living Landz** et on continue vers R7.B pour les bordures fractales.

**Falsifiable** : si Run C ≈ Run A (composite ne diffère pas significativement de Radial à cause de la crête écrasée par érosion) OU dissolution catastrophique comme R6.3 mf=1.0, alors **composition ne résout pas la richesse continentale** → pivot R7.B (Voronoï hiérarchique) prioritaire, OU lecture utilisateur 4 (limites 2D thin-sheet confirmées).

## 9. Demande de validation utilisateur

Tu valides :
- (a) Architecture nested structs (`CompositeRadialParams`, `CompositeOrogenicRidgeParams`, `CompositeCap`) ?
- (b) Q.f.1 formule ridge_only option (a) ?
- (c) Q.f.2 cap = peak_orogenic ?
- (d) Q.f.3 offset configurable, default 0 ?
- (e) Q.f.4 defaults (peak=1.20, base=0.85, σ=0.10, etc.) ?
- (f) 6 tests régression listés § 6 ?

Si tous OK → je commit R7.A.2.1 puis enchaîne R7.A.2.2 (impl + 6 tests + V2 spec + preset + UI). Sinon, précise les ajustements et je retravaille.
