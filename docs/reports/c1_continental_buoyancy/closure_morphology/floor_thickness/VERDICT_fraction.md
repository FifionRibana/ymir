# Verdict — la FRACTION (étendue) high-craton est le vrai levier, pas l'altitude seule

**Probe :** `c1_closure_morphology::probe_craton_fraction` (`#[ignore]`).
**Question :** le socle haut = cratons trop HAUTS (altitude, qui a buté sur 2 murs : borne
densité + couplage niveau-marin) ou trop NOMBREUX (fraction) ?

## Origine de la fraction (mesurée)

`build_phase_1_1_cratonic_mask` ([init.rs](../../../../crates/ymir-core/src/tectonics_c1/init.rs))
flague une plaque cratonique SSI son seed Voronoï est dans la **moitié gauche** du domaine
(`x < nx/2`). **Non ancré** sur aucun processus physique — un accident géométrique.

Conséquence : la fraction cratonique est **follement variable selon le seed** :

| seed | craton/continental | craton/terre-émergée |
|---|---|---|
| 99, 1337 | 0 % | 0 % |
| 4138 | 23 % | 27 % |
| 2026 | 34 % | 40 % |
| 42 | 49 % | 58 % |
| 1988 | 57 % | 65 % |

0 % à 65 % selon où tombent les seeds de plaques. Les seeds à forte fraction (1988, 42)
portent le déficit de plaine basse ; les seeds à 0 (99, 1337) sont Terre-like. C'est ce qui
fait l'agrégat ~34 % <500 m du handoff (moyenne sur 0-65 %).

## Ancrage réel (Terre)

Le craton précambrien sous-tend ~50-70 % de la croûte continentale — donc une grande AIRE
cratonique est réaliste. MAIS le **bouclier exposé haut** ≈ 10-20 % seulement ; le reste de
l'aire cratonique est **plateforme BASSE** (couverte de sédiments phanérozoïques). Notre
modèle rend TOUTE l'aire cratonique en HAUT bouclier → c'est ça l'erreur, pas l'aire en soi.

## Test du levier (seed 42) — décisif

Reconstruction avec le traitement cratonique mis à NORMAL (thickness_ratio 1.0 +
craton_resist 1.0 + craton_rho_crust None = la limite « pas de high-craton ») :

| variant | ALL land<500 | craton méd | non-craton méd | non-craton<500 |
|---|---|---|---|---|
| spécial (c1_default) | 25 % | 1814 m | 429 m | 54 % |
| croûte normale | **42 %** | 833 m | 467 m | 50 % |

→ Retirer le high-craton : **ALL<500 25 → 42 %** (+17 pts) ET **non-craton stable**
(429 → 467 m) → **le couplage niveau-marin se dissout** (les cratons ne dominent plus le
percentile-mer). C'est le levier le plus fort ET le plus propre — il évite le mur du
couplage qui plafonnait l'approche altitude à 38 %.

## Verdict : le levier est l'ÉTENDUE HIGH-CRATON (fraction × rendu tout-haut)

FRACTION et ALTITUDE sont deux faces du MÊME problème : **trop d'aire rendue en HAUT
bouclier.** L'axe FRACTION est le bon levier (l'axe altitude pur a buté 2 fois) car :
1. Le masque est non ancré ET follement seed-variable (0-65 %) — manifestement faux.
2. Réduire l'étendue high-craton → +17 pts (25→42 %) ET dissout le couplage niveau-marin
   (le mur de l'approche altitude).
3. Le non-craton est déjà Terre-like (médiane ~430-470 m, ~50-54 % <500 m) — le modèle SAIT
   faire la plaine basse ; il faut juste moins de high-craton par-dessus.

**Nuance :** « croûte normale » atteint 42 %, pas 52 % — les cellules ex-craton (plaques de
gauche) restent un peu hautes (833 m médiane) même en croûte normale (contexte d'intérieur
continental, pas l'identité cratonique). Donc l'étendue high-craton est le levier DOMINANT
mais pas l'unique ; un résidu ~10 pts est de la structure d'intérieur continental.

## Fix (espace, pas appliqué)

Réduire l'étendue HIGH-craton, deux voies complémentaires :
- **Ancrer/contrôler la fraction du masque** : remplacer la règle « moitié gauche »
  (non ancrée, 0-65 % aléatoire) par une fraction cratonique ANCRÉE et stable.
- **Bimodal bouclier/plateforme** (le « suite » déjà tracké, en fait le cœur du fix) :
  garder une minorité (~10-20 %) en haut bouclier, rendre la majorité de l'aire cratonique
  en BASSE plateforme. Ajoute de la plaine basse + dissout le couplage (la minorité haute ne
  domine plus le percentile).

→ PAS l'altitude cratonique uniforme (mauvais axe : borne densité + couplage). PAS l'érosion,
PAS un plancher global. Le levier est l'ÉTENDUE high-craton (masque + bimodal).
