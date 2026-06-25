# Verdict — diagnostic du plancher continental (amont #165)

**Probe :** `c1_closure_morphology::probe_floor_diagnosis` (`#[ignore]`).
**Run :** 2048² HD, seeds {42, 1337}, levier ×3 sur seed 42. ~800 s.
**Question :** le socle à 1000-2000 m vient-il (1) d'une naissance haute (flottabilité/
isostasie), (2) d'une érosion insuffisante, ou (3) d'un seuil mer mal placé ? Les trois ont
des fix OPPOSÉS — mesurer avant de coder.

## Données

| Mesure | Seed 42 | Seed 1337 | Terre |
|---|---|---|---|
| frac<500m **pré-érosion** (FBM nu) | 27.3 % | 42.4 % | ~52 % |
| frac<500m **post-érosion** | 30.0 % | 46.5 % | ~52 % |
| effet érosion | +2.7 pts | +4.1 pts | — |
| bande 1000-2000m pré→post | 46.1 → 46.4 % | 25.0 → 22.9 % | ~15 % |
| bande 3000+m pré→post | 0.5 → 0.2 % | 2.7 → 1.8 % | ~6 % |
| **levier ×3 gouttes** frac<500m | 30.0 → **30.4 %** | — | — |
| masse déposée sur terre réelle | 10 % | 18 % | — |
| masse déposée sous la mer réelle | 90 % | 82 % | — |
| Δalt depuis 1000-2000m (terre) | −76 m | −111 m | — |
| Δalt depuis 0-500m (terre) | −94 m | −34 m | — |
| bande 0-250m pré→post | 16.2 → 16.0 % | — | — |

CSV : `floor_prepost_bands.csv`, `floor_mass_balance.csv`.

## Verdict : CAUSE 1 — le socle naît HAUT (flottabilité / isostasie)

1. **Né haut — CONFIRMÉ.** Le bourrelet 1000-2000 m est présent dans le terrain
   PRÉ-érosion (FBM, zéro goutte) : 46 % @seed 42. Le déficit de plaine basse existe à la
   sortie isostasie+FBM, avant toute érosion.
2. **Érosion insuffisante — RÉFUTÉ.** L'érosion ne déplace la distribution que de
   +2.7/+4.1 pts et ne réduit PAS la bande 1000-2000 m. Le **test du levier** est décisif :
   ×3 gouttes ne gagne que **+0.4 pt** — l'érosion est saturée, ce n'est pas le levier.
   (L'indice « signature d'érosion insuffisante » du handoff est réfuté par la mesure, cf
   k_oro : l'intuition sur le levier était fausse.)
3. **Seuil mer — ÉCARTÉ.** Bande 0-250 m stable à ~16 % : pas de sliver artificiel à la
   côte, le contrat vertical (`c1_altitude_norm_to_metres`, sea=0.5) coupe au bon endroit.

## Mécanisme secondaire (réel, mais pas le levier)

82-90 % de la masse érodée se dépose **sous la mer réelle** (10-18 % sur terre). Cause :
l'érosion HD utilise `sea_level = 0.1` alors que la mer réelle est `0.5` — une goutte ne
dépose qu'en passant sous 0.1, donc au-delà de la côte réelle ET de tout le plateau
[0.1, 0.5]. L'érosion strippe (Δalt négatif partout) sans construire de plaine côtière.
C'est un défaut à corriger (aligner le `sea_level` de l'érosion sur le contrat vertical),
mais le test du levier prouve qu'il ne suffirait pas : le plancher est dominé par sa
naissance haute, pas par le devenir des sédiments.

## Levier du fix

→ **Isostasie** : abaisser le plancher continental (futur `continental_floor_m` gated dans
`IsostasyConfig`, None = byte-identique), conformément au handoff §4. PAS l'érosion, PAS le
seuil mer. Ne pas raboter les montagnes (3000+ déjà sous le réel : 0.2-1.8 % vs ~6 %).

## Sous-mesure recommandée AVANT de tuner

Décomposer le « né haut » : l'altitude isostatique brute (pré-FBM, à 64²/2048² coarse) vs
l'apport du FBM. Le FBM est censé être ~zéro-moyenne (perturbation de détail) ; si le
plancher haut vient de l'isostasie, le bouton est dans `IsostasyConfig` / l'ancre verticale.
À mesurer avant de choisir le paramètre exact.
