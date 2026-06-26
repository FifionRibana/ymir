# FIX — bimodal bouclier/plateforme cratonique (#165 relief)

**Code :** `Phase2InitParams::craton_shield_fraction: Option<f64>` (gated, None = pré-#165
byte-identique) + `select_shield_mask` dans `init_r7`. `c1_default` (= `Phase2InitParams::
default()`) → `Some(0.15)`. **Validé** par `probe_craton_bimodal` (64²) + `probe_hypsometry`
+ `probe_climate_biomes` (2048², cache).

## Le mécanisme

Le masque cratonique pose l'AIRE (ancienne lithosphère, ~50-70 % des continents — réaliste).
Le modèle rendait TOUTE cette aire en HAUT bouclier. `select_shield_mask` la rend BIMODALE :
- **bouclier** (~15 % de l'aire, ancré sur la part réelle de bouclier exposé ~10-20 %) :
  garde le traitement haut (thick + resist + dense) → boucliers usés hauts (~1100 m).
- **plateforme** (~85 %) : sort du masque → croûte continentale normale BASSE (déjà
  Terre-like). Sélection = mosaïque de value-noise cohérente seedée (patches contigus, pas
  poivre-et-sel), proportion par seuil-quantile → ~15 % quel que soit le seed.

Pourquoi cet axe (et pas l'altitude uniforme, qui a buté 2× — borne densité + couplage
niveau-marin) : on ne baisse pas tout le bloc, on en RECLASSE 85 % en bas. La minorité haute
ne domine plus le percentile-mer → le couplage ne se déclenche pas.

## Validation hypsométrie (2048², frac<500 m toute terre)

| seed | craton aire | avant | **après bimodal** | médiane après |
|---|---|---|---|---|
| 42 | 49 % | 30 % | **59 %** | 387 m |
| 1988 | 57 % | 18 % | **40 %** | 688 m |
| 4138 | 23 % | — | 44 % | 569 m |
| 2026 | 34 % | — | 43 % | 632 m |
| 99 | 0 % | 47 % | 47 % | 545 m |
| 1337 | 0 % | 46 % | 46 % | 557 m |

**Agrégat ~34 % → ~46 %** (cible Terre 52 %). Les seeds catastrophiques (42, 1988) réparés ;
les seeds sans craton (99, 1337) intacts (byte-identique en pratique — pas de craton à
reclasser). Couplage dissous (plancher plateforme/non-craton stable).

## Validation biomes (2048²) — forêt tempérée (forêt + rainforest)

| seed | forêt tempérée tot. | boréal/taïga | steppe (grassland) |
|---|---|---|---|
| 42 | 30 % | 12 % | 56 % |
| 99 | 36 % | 27 % | 36 % |
| 1337 | 20 % | 26 % | 51 % |
| 4138 | 27 % | 18 % | 53 % |
| 1988 | 25 % | 26 % | 38 % |
| 2026 | 23 % | 29 % | 48 % |

- ✅ Le **boréal/taïga** (qui dominait les intérieurs HAUTS) est RÉDUIT — les ex-cratons
  hauts boréaux sont maintenant bas + tempérés.
- ✅ La **forêt tempérée** a grossi : 20-36 % total (vs « mince 1-16 % » avant le fix).
- ⚠️ Le nouveau biome dominant est la **steppe tempérée (grassland) 36-56 %** : les
  plateformes basses sont SÈCHES (intérieur, faible précip) → grassland. **Géologiquement
  correct** (steppe eurasienne / Grandes Plaines / Pampa = intérieurs tempérés secs).

## Verdict

Le **fix relief est un succès** : hypsométrie 34 → 46 % (vers 52 %), intérieurs
boréaux → tempérés, forêt tempérée substantielle (20-36 %). Le couplage niveau-marin est
contourné (axe étendue, pas altitude). Gated, lib verte (468), cache invalidé par le param
(pas de bump ALGO).

## Résidus (sujets distincts, tracké)

1. **Plus de forêt vs steppe = axe CLIMAT** : les plateformes basses sont sèches → steppe
   (correct). Davantage de forêt exigerait plus de précipitation intérieure (le fond frontal
   / pénétration d'humidité), PAS le relief. Distinct.
2. **Structure d'intérieur continental non-craton** : les seeds sans craton (1988 a un
   non-craton plus haut) gardent un léger excès — résidu ~10 pts non-cratonique, à
   diagnostiquer séparément.
3. **Texture intra-cratonique** : le bouclier est un value-noise ; un raffinement
   (distance-au-cœur, sédiments de plateforme) viendrait après si « mosaïque » devient le
   manque dominant.
