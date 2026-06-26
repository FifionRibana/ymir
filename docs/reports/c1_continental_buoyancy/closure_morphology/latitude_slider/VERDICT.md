# Verdict — le slider de latitude est-il expressif ? (#165 jonctions latitudinales)

**Probe :** `c1_closure_morphology::probe_latitude_slider` (`#[ignore]`).
**Run :** seed 42, terrain érodé PARTAGÉ entre latitudes (la latitude ne change que le climat
dérivé, pas le relief — cache HIT). Latitudes 0/15/30/45/60/75°.

## Mesure (seed 42)

| lat | vent | belt | T_sea | T méd terre | précip mm/yr p10/méd/p90 | biomes dominants |
|---|---|---|---|---|---|---|
| 0° | **E→W alizés** | 1.00 | 27.0° | 24.4° | **2054**/2054/4073 | trop. rainforest 85 %, temp. rainforest 15 % |
| 15° | E→W alizés | 0.33 | 24.9° | 22.5° | 599/599/2392 | **savane 60 %**, temp. forest 21 % |
| 30° | **W→E westerlies** | 0.21 | 19.8° | 17.4° | 279/279/2376 | **steppe 67 %**, temp. forest 11 %, temp. rainforest 13 % |
| 45° | W→E westerlies | 0.55 | 12.1° | 9.7° | 448/448/1753 | **steppe 56 %**, temp. forest 19 %, boréal 12 % |
| 60° | **E→W polaires** | 0.53 | 1.9° | −0.4° | 213/213/617 | **desert 62 %**, boréal 24 %, toundra 14 % |
| 75° | E→W polaires | 0.20 | −10.5° | −12.8° | **32**/32/191 | **toundra 100 %** |

## Ce qui MARCHE déjà (largement câblé — ne pas réinventer)

1. **Direction du vent — CÂBLÉE & fonctionne.** `wind_zonal_dir(lat)` bascule par bande et le
   transport scanne dans le bon sens : E→W alizés (<30°), W→E westerlies (30-60°), E→W
   polaires (≥60°). L'orographie est calculée au bon vent par bande.
2. **Température zonale — CÂBLÉE.** `sea_level_temperature(lat)` = 27 °C équateur → −25 °C
   pôle (`(lat/90)^1.8`). T médiane terre 24.4 → −12.8 °C.
3. **Profil de précip de base zonal — CÂBLÉ.** `belt_factor(lat)` × `e_sat(T_sea(lat))` donne
   le vrai profil méridien : ITCZ très humide (2054 mm @0°), déclin vers les subtropiques,
   pôle sec (32 mm @75°).
4. **Biomes DISTINCTS par latitude — le slider EST expressif :** jungle (0°) → savane (15°)
   → steppe (30°) → tempéré (45°) → froid-sec (60°) → toundra (75°). Signatures lisibles.

## Le GAP (2 raffinements de PROFIL, PAS de câblage de direction)

1. **Le désert subtropical ~30° est TROP FAIBLE.** Médiane 279 mm → **steppe**, pas désert
   (ligne 250). Le minimum subtropical de `belt_factor(30°)=0.21` ne descend pas l'intérieur
   sous 250. Réel : la ceinture saharienne (~25-30°, subsidence) est désertique. → creuser le
   minimum subtropical (belt_factor un peu plus bas vers 25-30°, ou un terme de sécheresse de
   subsidence), pour que l'intérieur 30° passe < 250 = désert.
2. **Le front polaire ~60° n'est PAS humide.** Médiane 213 mm → lu comme **desert froid**,
   pas « front polaire humide ». Cause : T_sea(60°)=1.9 °C → `e_sat` faible (l'air froid
   tient peu d'humidité) écrase le `frontal_base` malgré `belt_factor(60°)=0.53` correct. Le
   60° continental EST sec (légitime), mais l'enrichissement frontal des rails de tempête
   (storm-track maritime) n'est pas modélisé → pas de signature « humide ». À trancher selon
   la cible produit (un boost de front polaire si « 60° humide » est voulu).

## Observation transversale

L'intérieur est **PLAT au frontal_base** à chaque latitude (p10 = médiane partout) : le profil
zonal pose un PLANCHER uniforme par bande, la variété orographique n'apparaît qu'au p90.
Donc la variété vient (a) de la latitude (le profil zonal) + (b) de l'orographie (p90) — PAS
de variété intra-bande au niveau médian. C'est le motif de saturation frontale déjà connu
(cf les commentaires `k_frontal`).

## Verdict

Le slider est **LARGEMENT CÂBLÉ et EXPRESSIF** — direction du vent, température, profil de
précip zonal et biomes changent tous avec la latitude, produisant des signatures distinctes.
**Le câblage de direction ne manque RIEN.** Le seul GAP est 2 raffinements de PROFIL de
précip : (1) creuser le désert subtropical ~30° (279 → < 250 mm), (2) optionnellement
enrichir le front polaire ~60° (sinon il reste froid-sec, pas humide). Diagnostic seulement —
le fix (ajuster `belt_factor` / un terme subsidence / un boost front polaire) suit ce verdict.
