# SIGNAL 2 — le couplage niveau-marin bloque l'abaissement isolé des cratons

**Probe :** `c1_closure_morphology::probe_craton_density_sweep` (variante keel + densité).
**Run :** seed 42, état tectonique figé, balayage du modèle de compensation cratonique.

## Contexte

Voie A choisie : abaisser les cratons via une isostasie compositionnelle (Jordan). Deux
formes testées : (1) **densité effective de colonne** (proportionnelle à S̃, le knob
`craton_rho_crust` réinterprété) ; (2) **double couche additive** (un sink de freeboard
constant = keel froid sous la croûte, params `craton_root_thickness_km` / `_density_excess`).

## Mesure — la double couche additive (keel constant)

| keel (km × Δρ) | craton méd | craton<500 | non-craton méd | non-craton<500 | ALL land<500 |
|---|---|---|---|---|---|
| no root (2900) | 1814 m | 5 % | **429 m** | 54 % | **25 %** |
| 150 km +35 | 1438 m | 6 % | 698 m | 35 % | 21 % |
| 180 km +45 | 1234 m | 6 % | 946 m | 16 % | 11 % |
| 200 km +45 | 1178 m | 6 % | 1026 m | 7 % | 7 % |
| 220 km +55 | 978 m | 8 % | 1263 m | 1 % | 4 % |
| 250 km +66 | 708 m | 18 % | **1541 m** | 0 % | **9 %** |

→ Le keel constant abaisse bien les cratons MAIS fait exploser le non-craton (429 → 1541 m)
et DÉGRADE le global (25 → 9 %). **Contre-productif.**

## Cause racine (tracée)

Le sink CONSTANT s'applique à TOUS les cratons, y compris les cratons FINS (érodés, S̃≈0.01)
→ il les tire en très négatif → effondre `h_min` → effondre le niveau marin
(`PercentileCapped 0.92` : `h_sea = h_min + 0.4·(h_cap − h_min)`) → **tout le non-craton se
retrouve plus haut au-dessus d'une mer abaissée.** La forme DENSITÉ (proportionnelle à S̃)
ne souffre pas de ça (un craton fin reçoit un petit ajustement, reste près de 0, `h_min`
préservé) — d'où la dérive non-craton modérée du sweep densité (429 → 547 m).

**Mathématiquement** : densité constante ρ_eff ⟺ keel d'épaisseur PROPORTIONNELLE à la
croûte (bien comporté) ; keel d'épaisseur CONSTANTE ⟺ densité effective explosant sur les
cratons fins (casse `h_min`). La forme proportionnelle (densité) est la bonne ; l'additive
constante est à rejeter.

## Le couplage niveau-marin — le signal plus profond (vaut AUSSI pour la densité)

Même avec la forme densité (bien comportée), le niveau marin = percentile 0.92 de
l'altitude. Les cratons occupent le haut de la distribution (58 % de la terre, hauts). Les
abaisser baisse le percentile → baisse la mer → REMONTE le non-craton. Mesure (sweep
densité) : non-craton 429 → 547 m de 2900 à 3100. Conséquence : viser la cible globale
(ALL<500 ≈ 52 %) via la densité exige ρ≈3200 où les cratons tombent à 124 m (SOUS le worn-
shield) ; à ρ≈3100 les cratons sont à 573 m (worn-shield correct) mais ALL<500 plafonne à
38 % (le gain craton est mangé par la remontée non-craton).

| densité | craton méd | non-craton méd | ALL land<500 |
|---|---|---|---|
| 2900 (actuel) | 1814 m | 429 m | 25 % |
| 3100 | 573 m (worn-shield ✓) | 547 m | 38 % |
| 3200 | 124 m (trop bas) | 541 m | 51 % |

→ **On peut amener les cratons à la hauteur worn-shield (~573 m) de façon physique
(densité de colonne effective ~3100, ancrable croûte+keel) — gain réel ALL<500 25 → 38 % —
MAIS le couplage niveau-marin empêche d'atteindre la cible 52 % sans sur-couler les cratons.
La cible pleine exige de traiter le couplage niveau-marin (un 2nd chantier de formulation),
pas seulement la densité cratonique.**

## Décision (en attente user)

- A1. Appliquer le fix craton physique (densité de colonne effective ~3100, cratons →
  worn-shield ~573 m, ALL<500 25 → 38 %, gain réel) et MESURER les biomes (forêts
  partiellement débloquées ?). Tracker le couplage niveau-marin comme 2nd signal.
- A2. Traiter d'abord le couplage niveau-marin (découpler le percentile de l'abaissement
  cratonique : calculer la mer sur une référence stable / hors-craton) PUIS le fix craton —
  pour atteindre la cible pleine.
- Rejeter la double couche additive (keel constant) : casse `h_min`.
