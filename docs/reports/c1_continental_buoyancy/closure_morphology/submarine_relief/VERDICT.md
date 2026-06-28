# Verdict — qu'y a-t-il sous le niveau marin aujourd'hui ? (sous-marin)

**Probe :** `c1_closure_morphology::probe_submarine_relief` (`#[ignore]`).
**Run :** 6 seeds, 2048² HD (terrain partagé via cache), cellules sous le niveau marin
(`n ≤ 0.5`), profondeur via `c1_altitude_norm_to_metres` (océan = norm [0,0.5] → [−5650, 0] m).

## Mesure

**Hypsométrie océanique agrégée (6 seeds) :**

| bande | fraction de l'océan |
|---|---|
| 0-200 m (plateau / shelf) | **1.5 %** |
| 200-1000 m (talus) | 3.2 % |
| 1000-3000 m | **95.2 %** |
| 3000-5000 m (abysse) | **0.1 %** |
| 5000+ m (fosses) | 0.0 % |

**Par seed :** océan 73-92 % de la carte ; profondeur moyenne ~**2400-2540 m**, médiane
~2600 m, **la plus profonde seulement ~2850-3180 m**, std ~430-620 m.

**Profil côte→large (profondeur moyenne par distance) :** 0-25 km ~1700-2100 m, 25-100 km
~2500 m, 100-300 km ~2630-2690 m, 300+ km ~2620-2640 m (ou vide — petits bassins).

## Lecture (réponses aux 3 questions)

1. **L'échelle est SOUS-UTILISÉE.** 95 % du fond est dans UNE seule bande (1000-3000 m), la
   plus profonde ~3000 m, alors que le contrat descend à −5650 m. **Aucune plaine abyssale**
   (3000-5000 m = 0.1 %), **aucune fosse**. Le fond ne descend jamais aux profondeurs
   abyssales réelles (−4000/−5000 m).
2. **Structure = DALLE mi-profonde quasi-uniforme (~−2600 m), PAS bimodale.** std ~500 m sur
   une moyenne de 2600 m → faible variation. 95 % à une profondeur ~uniforme. Le réel est
   BIMODAL (plateau peu profond + abysse profond, peu entre) ; le nôtre est une dalle plate
   à mi-profondeur.
3. **Côte→large : gradient présent mais COMPRIMÉ.** La profondeur augmente bien avec la
   distance (1800 → 2700 m) — donc PAS totalement plat (le gradient d'âge Stein-Stein
   existe). MAIS : **quasi aucun plateau continental** (0-200 m = 1.5 % ; à 0-25 km de la côte
   on est déjà à ~1800-2100 m, pas 0-200 m), et le gradient **sature à ~2700 m** (jamais
   d'abysse). La côte plonge presque direct à ~1800 m, sans plateau.

## Cause racine

Le fond océanique est posé par **Stein-Stein** (`depth = 2500 + 350·√âge`, subsidence
thermique). Avec les âges JEUNES de la run C1 courte (pas d'expansion océanique / seafloor
spreading pour vieillir la croûte), `√âge` reste petit → profondeur plafonnée à
~2500-3000 m partout → la dalle uniforme. Deux mécanismes manquent :
- **Plateau continental** : la marge continentale immergée (croûte continentale amincie près
  des côtes) devrait être PEU profonde (0-200 m). Stein-Stein ne s'applique qu'aux cellules
  océaniques → la transition continent→océan saute directement à la profondeur océanique
  (~1800 m), sans plateau.
- **Plaine abyssale** (−4000/−5000 m) : croûte océanique vieille → absente faute d'âge.

## Verdict — structuré-mais-MAUVAISE-FORME (à corriger, pas à générer de zéro)

Une BASE existe (Stein-Stein pose un gradient côte→large par âge — ne pas réinventer), mais
sa FORME est fausse : une dalle uniforme à ~−2600 m au lieu du profil bimodal réel
**plateau → talus → plaine abyssale**. Ce qui manque :
1. **Plateau continental** (0-200 m près des côtes) — quasi absent (1.5 %). **Fonctionnel pour
   Living Landz** (eaux côtières jouables, forme des côtes) → priorité.
2. **Plaine abyssale** (−4000/−5000 m) — absente (l'échelle dispo est inutilisée). Réalisme
   cartographique.
3. **Bimodalité** — le fond est mono-modal (dalle) au lieu de bimodal (plateau + abysse).
4. Fosses / dorsales — absentes (réalisme, basse priorité).

Pistes (cadrage APRÈS ce verdict) : un **profil bathymétrique plateau→talus→abysse** basé sur
la distance à la côte (le `dist_to_coast` existe), OU vieillir la croûte océanique (âge) pour
que Stein-Stein descende à l'abysse + un terme de plateau sur la marge continentale immergée.
Diagnostic seulement — le cadrage de ce qu'on génère (et jusqu'où : plateau fonctionnel vs
abysse cartographique) suit ce verdict.

## FIX appliqué — re-map bathymétrique plateau→talus→abysse (chemin 2)

`FbmUpscaleConfig::bathymetry: Option<BathymetryProfile>` (gated, `None` = byte-identique),
appliqué dans `upscale_from_c1` APRÈS l'érosion sur les cellules sous-marines uniquement
([`terrain::bathymetry`]). Re-map chaque cellule océan vers l'enveloppe `dist_to_coast`
(plateau ~30 km → talus → abysse ~−4500 m), MODULÉE par sa déviation relative existante
(`env·(1+texture·dev)`) → la texture FBM/Stein-Stein survit (pas d'océan en oignon).

**Mesuré (probe, 6 seeds), avant → après :**

| bande | avant (dalle) | après |
|---|---|---|
| 0-200 m plateau | 1.5 % | **11.6 %** |
| 200-1000 m talus | 3.2 % | 5.5 % |
| 1000-3000 m | 95.2 % | 6.8 % |
| 3000-5000 m abysse | 0.1 % | **64.7 %** |
| 5000+ m | 0.0 % | 11.4 % |

Le plus profond ~3000 → **~5000-5600 m** (échelle utilisée), std ~500 → **~2000 m** (varié),
distribution **BIMODALE** (plateau + abysse, talus étroit). **Invariants tenus** : fraction
océan IDENTIQUE avant/après sur les 6 seeds (79/92/77/83/76/73 %) → côte / masque terre-mer /
hypsométrie terrestre inchangés (le re-map ne touche que sous la mer et garde les cellules
immergées). Tests unitaires (`envelope_is_monotone`, `remap_keeps_coastline_and_deepens_offshore`)
+ lib verte (471). Cache invalidé par le champ (FbmUpscaleConfig serde), pas de bump ALGO.

**Calibration** : plateau 11.6 % (un peu au-dessus du ~7-8 % Terre — généreux/fonctionnel,
eaux côtières jouables sur ces cartes régionales à côtes denses ; `shelf_width_km` = le knob
si on veut plus serré). **Différé** : chemin 1 (spreading tectonique → abysse/dorsales/fosses
émergents de la physique), fosses/dorsales réelles (ne suivent pas `dist_to_coast`).
