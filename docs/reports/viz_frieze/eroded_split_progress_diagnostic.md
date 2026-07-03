# Diagnostic — découper le nœud « Relief/Eroded » + progression fine (frise, suite e)

**But** : la frise (e) a UN nœud « Relief » = la phase `Eroded` (tectonique + upscale + érosion
+ bathymétrie en bloc). La maquette voulait 3 nœuds (Tectonique → Relief → Érosion) et
l'utilisateur veut une **barre continue** (l'érosion à 20 % → segment rempli à 20 %). Ce
diagnostic vérifie ce qui est **réellement** découpable et mesurable AVANT de coder — ne pas
promettre une barre dont on n'a pas la donnée. Lecture seule.

## Ce que fait réellement la phase Eroded

`cached_c1_eroded` ([cached_product.rs:98](../../../crates/ymir-core/src/tectonics_c1/cached_product.rs#L98)),
dans sa closure de calcul (MISS uniquement) :

```
init_c1_state_phase_2_r7            (init 64², instantané)
run_with_closures(…, |_,_| {})      (TECTONIQUE coarse — callback step JETÉ)
upscale_from_c1(…)                  (UPSCALE HD → ÉROSION → BATHYMÉTRIE)
```

et `upscale_from_c1` ([production_upscale.rs:256](../../../crates/ymir-core/src/tectonics_c1/production_upscale.rs#L256)) :

```
c1_production_altitude_craton       (altitude coarse : isostasie + Stein-Stein — instantané)
upscale_with_fbm(…)                 (RELIEF HD : bicubique + FBM — opaque, pas de callback)
run_erosion(…, |_,_,_| true)        (ÉROSION HD — callback batch JETÉ)   ← le point clé
apply_bathymetry_profile(…)         (BATHYMÉTRIE — 1 passe + BFS, quasi instantané)
```

## Partie A — découpage en nœuds

**Séparables ?** Oui, mais à deux niveaux :
- **Tectonique** (`run_with_closures`) est un appel DISTINCT de `upscale_from_c1` dans la closure
  de `cached_c1_eroded` → séparable directement.
- **Relief (upscale)**, **Érosion**, **Bathymétrie** sont des sous-appels INLINE dans
  `upscale_from_c1` → il faut lui passer un callback (ou le découper en sous-fonctions) pour
  émettre un event entre chacun.

**Poids temporel (estimé, 2048²)** :

| sous-étape | ~durée | progression interne |
|---|---|---|
| tectonique coarse (300 pas @64²) | ~0.25 s | step N/300 (déterminé, mais quasi instantané) |
| upscale FBM (2048²) | ~quelques s | aucune (opaque) |
| **érosion (4 M gouttes)** | **~100 s (dominante)** | **callback par batch (40 batches)** |
| bathymétrie | < 1 s | aucune (opaque) |

**Mapping maquette Tectonique → Relief → Érosion** : réalisable —
Tectonique=`run_with_closures`, Relief=`upscale_with_fbm` (relief HD avant érosion),
Érosion=`run_erosion` ; **bathymétrie fondue** dans Érosion (négligeable). Les 2 « nœuds
différés » de (d2) sont donc récupérables ici.

**⚠️ Caveat CACHE (structurant)** : `cached_c1_eroded` est **UNE** entrée de cache
(« eroded » = tecto+upscale+érosion+bathy en un seul `.raw`). Sur un **HIT**, aucune sous-étape
ne tourne → rechargement instantané. Donc le découpage en 3 nœuds + la barre n'ont de sens
**qu'à froid (MISS)** ; à chaud, les 3 nœuds s'allument ensemble instantanément (cache). Ce
n'est PAS un découpage en unités de cache séparées — juste un affichage de progression pendant
un build à froid.

## Partie B — progression fine (la barre à 20 %)

4. **`run_erosion` a un callback exploitable ?** OUI.
   [hydraulic.rs:113](../../../crates/ymir-core/src/erosion/hydraulic.rs#L113) :
   `progress_callback(batch_end, total, &heightmap) -> bool`, appelé **après chaque batch**
   ([:165](../../../crates/ymir-core/src/erosion/hydraulic.rs#L165)) ; `false` = annuler. Config
   HD : `num_droplets = 4 M·(target/2048)²`, `batch_size = 100 k` → **40 batches @2048** (≈2.5 %
   par tick), 10 @1024, ~2-3 @512. → **vraie barre % continue** sur l'érosion. Actuellement le
   callback est **jeté** (`|_,_,_| true`). Bonus : renvoyer `false` = **annulation de l'érosion
   en plein calcul** (le cancel actuel n'agit qu'entre phases).

5. **`cached_c1_eroded` peut-il PROPAGER ce callback ?** OUI, coût MOYEN. Il faut faire remonter
   le callback de `run_erosion` à travers `upscale_from_c1` puis `cached_c1_eroded` jusqu'au
   worker. Ripple des signatures : `upscale_from_c1` a ~7 sites d'appel, `cached_c1_eroded` ~9
   (surtout des tests) → **contournable par des variantes `_with_progress(...)`** (les fonctions
   actuelles deviennent des wrappers no-op) → **zéro churn de tests**. Calcul byte-identique, **pas
   de changement de clé de cache**. Touche le core mais de façon **additive** (nouvelles fns).

6. **Tectonique — vraie barre N/300 ?** OUI techniquement (`run_with_closures` a le callback
   `|step, state|`, actuellement jeté dans `cached_c1_eroded`). MAIS ~0.25 s à 64² → la barre se
   remplit en un clin d'œil (valeur visuelle faible ; un flash). 

7. **Upscale / bathymétrie** : aucun callback interne (opaques), courtes → **waiter** (ou fondre).

## Verdict

- **Découper « Relief » en Tectonique / Relief / Érosion** : **OUI** (sous-étapes séparables ;
  bathy fondue dans Érosion). Mais **cosmétique à froid uniquement** (cache « eroded » monolithique
  → à chaud, tout instantané). Pas d'unités de cache séparées.
- **Barre de progression FINE** :
  - **Érosion → OUI, vraie barre %** (callback par batch propagé) — **c'est LA phase longue
    (~100 s), la barre y a une vraie valeur**. + annulation mid-érosion en bonus.
  - **Tectonique → barre N/300 possible mais flash** (~0.25 s) — optionnel, faible valeur.
  - **Upscale / bathymétrie → waiter** (opaques, courtes).
  - Donc la « barre continue entre nœuds » se justifie surtout sur le **segment Érosion** ; les
    autres segments *sautent* (opaques) ou *flashent* (tectonique).
- **Ampleur du travail core** : **MOYEN**. Additif via variantes `_with_progress` de
  `upscale_from_c1` + `cached_c1_eroded` (wrappers ⇒ pas de churn), propagation du callback
  d'érosion (+ step tectonique), nouveaux events sous-phase, la frise les consomme. Pas de
  changement d'algorithme ni de clé de cache. Risque faible (wrappers), byte-identique.

**Recommandation** : découper en 3 nœuds (Tecto/Relief/Érosion) + **barre % réelle sur
l'Érosion** (le gros gain) via la propagation du callback ; Tectonique en barre N/300 (ou
waiter) ; Upscale/bathy en waiter. Tout ceci ne s'anime qu'à froid (MISS) — à chaud (HIT) la
frise saute (cache), cohérent avec (e). Diagnostic seulement — le fix suit ce verdict.
