# Step 12 R7.A.1.1 — Orogenic continental profile geometric spec

**Branche** : `112-step-12-interleaved-tectonic-erosion-workflow`
**Statut** : décisions géométriques documentées avant impl. R7.A.1.2 enchaîne sous validation.

## 1. Contexte amont

R6.3 a empiriquement établi (Option β + C.4) que **aucune amplification paramétrique simple** (mf, evolution_rate, craton_amp) ne débloque le sweet spot Living Landz tout en partant de l'IC actuelle (active_medley + RadialProfile + Voronoï rectiligne). Le verdict EVO.C identifie l'IC comme cause structurelle.

R7.A adresse une dimension de cette cause : **profil radial pic central** est inadapté pour produire des chaînes orogéniques persistantes. Sur Terre, les chaînes sont **linéaires** (axe long, perpendiculaire à la convergence), pas radiales. R7.A.1 introduit un profil orogenic complémentaire à `RadialProfile`.

R7.A.1 n'adresse PAS la dimension Voronoï rectiligne → R7.B.

## 2. Architecture : variant sibling, pas refactor

Le user a proposé un nouvel enum `ContinentalProfile { RadialPeak, Orogenic }` nested dans une struct `InitConfig`. Le code existant utilise déjà `InitMode` (enum unique) avec une variant `RadialProfile { … }` qui est précisément le "RadialPeak" du user.

**Décision** : ajouter `InitMode::Orogenic { … }` comme variant sibling. Surface minimale, dispatch identique au pattern existant, **aucune signature publique cassée**. `RadialProfile` reste comme variant.

Conséquences :
- `InitMode::default()` reste `Uniform` (pas de régression contre Step 13)
- Tous les presets référant `RadialProfile` ou `RadialProfileWithFBM` sont bit-identical
- Le nouveau preset orogenic instancie `InitMode::Orogenic { … }` explicitement

## 3. Décisions géométriques

### Q.geom.1 — Orientation de la crête : **PCA principal axis (option a)**

Pour chaque plate continentale, l'orientation de la crête est déterminée par le **premier vecteur propre** de la matrice de covariance des coordonnées de ses cellules.

**Pourquoi** :
- Déterministe (même seed → même PCA → même orientation)
- "Naturel" : suit la dimension longue de la plate, ce qui matérialise un alignement émergent de la géométrie Voronoï elle-même
- Ne dépend pas de la connaissance des voisins (qui appartient à R7.B)
- Robuste pour plates rectangulaires/allongées (covariance bien définie)

**Caveats techniques** :

1. **Périodicité torique** : les cellules d'une plate peuvent traverser les bords du domaine (wrap autour `i = 0` ↔ `i = nx-1`). Une PCA naïve sur coordonnées `(i, j)` brutes mettrait ces cellules artificiellement loin → centroïde faux + axe PCA faux.

   **Solution** : calcul du centroïde via **moyenne circulaire** (project sur le cercle unité, moyenne vectorielle, arctan), puis unwrap chaque cellule en choisissant l'image périodique la plus proche du centroïde. PCA s'applique ensuite sur ces coordonnées unwrappées.

2. **Plate dégénérée (peu de cellules ou cellules colinéaires)** : si une plate a < 5 cellules ou si sa covariance est ~rank-1 (toutes colinéaires), la PCA donne un axe peu fiable ou indéterminé.

   **Solution fallback** : si `det(cov) < epsilon` ou `n_cells < 5`, utiliser une orientation fixe `(1, 0)` (axe x). Documenter le fallback ; l'utilisateur peut alors choisir `OrogenicOrientation::Fixed { angle_rad }` pour le preset.

3. **Sign ambiguity** : le PCA donne un axe (direction non orientée). C'est suffisant pour le profil orogène (symétrique en `d_along` autour du centroïde).

### Q.geom.2 — Position de la crête : **centroïde de la plate (option i)**

La crête est centrée sur le centroïde géométrique de la plate (calcul périodique-aware comme en Q.geom.1).

**Pourquoi** :
- Cohérent avec la philosophie "la plate elle-même définit son centre de gravité"
- Pas besoin de connaître les voisins (option iii bord de plate appartient à R7.B après que les bordures aient été redéfinies)
- Géométriquement bien défini pour toutes les plates

**Caveats** :
- Pour des plates très allongées, le "centroïde" n'est pas forcément le point physique le plus naturel pour ancrer une chaîne. Acceptable pour MVP.
- Si la plate est concave (rare avec Voronoï), le centroïde géométrique peut tomber hors plate. Le profil restera valide mathématiquement (distance dans le plan, pas BFS) mais le pic de la crête sera "atténué" car les cellules les plus hautes seront sur la frange.

### Q.geom.3 — Formule mathématique du profil

Pour chaque cellule continentale `(x, y)` (coordonnées `i + 0.5, j + 0.5` en unités de cellules), de la plate `p` :

```
1. Centroïde périodique-aware : (cx_p, cy_p)
2. PCA axis périodique-aware : direction unitaire (ux_p, uy_p)
3. Vecteur de la cellule au centroïde (en coordonnées unwrappées)
   dx = unwrap(x - cx_p,  nx)
   dy = unwrap(y - cy_p,  ny)
4. Projection :
   d_along = dx · ux_p + dy · uy_p     (signed, le long de l'axe)
   d_perp  = dx · (-uy_p) + dy · ux_p  (signed, perpendiculaire)
5. Échelle : L_plate = max BFS distance to inter-plate boundary pour la plate p
   half_length = half_length_ratio · L_plate    (default ratio = 0.40)
   width_sigma = width_sigma_ratio · L_plate    (default ratio = 0.08)
6. Modulation longitudinale (smoothstep clamp décroissant aux extrémités) :
   t_along = (1 - |d_along| / half_length).clamp(0, 1)
   long_mod = smoothstep(t_along) = t_along² · (3 - 2 · t_along)
7. Profil transversal (Gauss) :
   trans_profile = exp(-(d_perp / width_sigma)²)
8. Combinaison :
   ridge_amount = long_mod · trans_profile        ∈ [0, 1]
   S̃ = base_continental_value
       + (peak_value - base_continental_value) · ridge_amount
9. Cellules océaniques : S̃ = oceanic_value (uniforme, comme RadialProfile)
```

**Garanties** :
- `ridge_amount ∈ [0, 1]` strictly → `S̃ ∈ [base, peak]` strictement (ou exactement `base` au centroïde si `width_sigma → 0`, donc clamp en pratique)
- Symétrique en `±d_along` et `±d_perp` (la crête est un segment, pas une demi-droite)
- Continu et C¹ partout (smoothstep + gaussienne)
- À `d_perp = 0` et `d_along = 0` : `S̃ = peak_value` exact
- Loin de l'axe (`|d_perp| > 3 · width_sigma` ou `|d_along| > half_length`) : `S̃ → base_continental_value`

**Pourquoi cette formule** :
- La gaussienne en `d_perp` donne un flanc symétrique et naturellement abrupt en `width_sigma`. Plus simple qu'un smoothstep transversal ; le décroissance Gaussien à `exp(-9)` à `3σ` rend l'effet largement local.
- Le smoothstep longitudinal en `d_along` donne une crête plateau au centre et des extrémités progressivement basses (smoothstep est C¹ continu aux bornes ; à `t_along=0` la dérivée est 0, donc pas de discontinuité au bout de la crête).
- Le L_plate-relative scaling garantit que chaque plate a une crête à l'échelle de sa propre taille, sans paramètre absolu (cohérent avec l'esprit de `RadialProfile`).

## 4. Paramètres par défaut MVP

| Paramètre | Valeur | Sens physique |
|---|---|---|
| `peak_value` | **1.20** | S̃ équivalent altitude ~5-8 km (crête Himalaya / Andes proxy) |
| `base_continental_value` | **0.85** | S̃ équivalent plaine continentale (~28 km crust, descend d'une couche reference) |
| `oceanic_value` | **0.20** | S̃ équivalent océanique (~7 km) — match RadialProfile default |
| `half_length_ratio` | **0.40** | Crête longe 80% de la plate (`= 2 × 0.40 × L_plate` en diamètre projeté) |
| `width_sigma_ratio` | **0.08** | Largeur σ ~ 8% de L_plate, soit ~5-10 cellules à 32-64² |
| `orientation` | `PlateMainAxisPca` | PCA fallback `Fixed { angle_rad: 0.0 }` si dégénéré |

Aucun de ces défauts n'est tuné. Ils sont des **points de départ** dérivés des proxies physiques (Himalaya 250 km / 2500 km ≈ 0.1 = 2·width_sigma, half_length 80% de la plate ~= ratio 0.4). Le verdict R7.A.1.3 dira si ces défauts produisent un visuel acceptable ou s'ils demandent un ajustement (avec documentation explicite — pas de tuning silencieux V2 vigilance).

## 5. Critères d'acceptance R7.A.1.2 (régression bit-identical)

Tests régression à exécuter sans modification existante (Step 13 baselines) :

1. **`radial_profile_default_unchanged`** : `InitMode::RadialProfile { … }` produit même Field2D bit-identical qu'au commit pré-R7.A.1
2. **`orogenic_continental_fraction_correct`** : `InitMode::Orogenic { … }` ne change pas le nombre de cellules continentales vs RadialProfile (même Voronoï, même classif)
3. **`orogenic_oceanic_uniform`** : océaniques retournent exactement `oceanic_value` (pas de bruit)
4. **`orogenic_peak_at_centroid`** : la cellule la plus proche du centroïde d'une plate continentale a `S̃` proche de `peak_value`
5. **`orogenic_decays_perpendicular`** : à distance > `3 · width_sigma` perpendiculaire à l'axe, `S̃ → base_continental_value`
6. **`orogenic_determinism`** : seed identique → byte-equal output (deterministic PCA)

## 6. Hors scope explicite

- **Pas de modification Voronoï** (R7.B)
- **Pas de profil multi-crête** (1 seule crête par plate continentale)
- **Pas de yielding/crust dynamics** sur la crête au-delà de ce que la chaîne tectonique fait déjà
- **Pas de craton spécial sur la crête** (cratonic factor reste calculé comme avant via distance à la bordure inter-plate)
- **Pas de FBM superposé** (analogue à `RadialProfileWithFBM` pourrait suivre dans R7.A.2 si pertinent ; pas dans R7.A.1)
- **Pas de support Step 8/Step 9 yielding feedback** (orogenic c'est juste IC, la rheologie/yielding reste inchangée)

## 7. Anti-patterns à éviter

- ❌ Modifier signature `InitMode::RadialProfile` → break regression Step 13
- ❌ Tuner paramètres defaults pour faire passer le visuel sans documenter
- ❌ Implémenter PCA naïf sans gestion périodicité (résultats faux aux bords)
- ❌ Implémenter avant validation user des décisions ci-dessus

## 8. Demande de validation

Tu valides :
- (a) Architecture : variant sibling `InitMode::Orogenic` au lieu de refactor en `ContinentalProfile` ?
- (b) PCA + centroïde périodique-aware (Q.geom.1 + Q.geom.2) ?
- (c) Formule smoothstep × gaussienne (Q.geom.3) ?
- (d) Defaults peak=1.20 / base=0.85 / half=0.40 / sigma=0.08 ?
- (e) Fallback `Fixed { angle_rad: 0.0 }` pour plates < 5 cellules ou PCA dégénérée ?
- (f) Q.config validation : (ii) single shot mf=1.0, evo=0.10, craton_amp=3 ?

Si tous OK → je commit R7.A.1.1 DOC + enchaîne R7.A.1.2 (impl + 6 tests régression).
