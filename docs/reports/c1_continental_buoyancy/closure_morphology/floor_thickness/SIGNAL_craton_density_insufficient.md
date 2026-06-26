# SIGNAL — le réglage densité cratonique ne suffit pas dans la borne physique

**Probe :** `c1_closure_morphology::probe_craton_density_sweep` (`#[ignore]`).
**Run :** seed 42 (porteur de cratons), état tectonique 64² construit une fois, balayage de
`craton_rho_crust` sur le mapping d'altitude (la densité n'entre pas dans l'évolution de S̃).

## Contexte

Le verdict ([VERDICT.md](VERDICT.md)) a établi que le socle haut = les cratons (médiane
émergée 1814 m). Le fix demandé : abaisser la freeboard cratonique vers la hauteur
worn-shield (~400-600 m) via `craton_rho_crust`, ANCRÉ dans la borne physique (croûte
continentale ≤ ~3000 kg/m³ ; 3300 = manteau). Discipline : si la borne ne suffit pas,
REMONTER un signal de formulation, ne pas forcer un paramètre hors borne.

## Mesure (seed 42)

| rho (kg/m³) | craton médiane | craton <500m | non-craton médiane | non-craton <500m | **ALL land <500m** | borne |
|---|---|---|---|---|---|---|
| 2750 (penalty off) | 2711 m | 7 % | 349 m | 60 % | 24 % | ok |
| **2900 (c1_default actuel)** | 1814 m | 5 % | 429 m | 54 % | 25 % | ok |
| 2950 | 1508 m | 5 % | 450 m | 53 % | 26 % | ok |
| **3000 (borne haute)** | **1197 m** | 5 % | 476 m | 52 % | **26 %** | ok |
| 3050 | 881 m | 8 % | 520 m | 49 % | 26 % | OUT |
| 3100 | 573 m | 30 % | 547 m | 47 % | **38 %** | OUT |
| 3200 | 124 m | 86 % | 541 m | 47 % | **51 %** | OUT |

Cible Terre : ~52 % de terre < 500 m ; worn shield ~300-600 m.

## Verdict du signal

1. **Le knob densité est INSUFFISANT dans la borne.** À `rho ≤ 3000`, la médiane cratonique
   plancher à **1197 m** — loin de la cible 400-600 m. La terre < 500 m ne monte que de
   25 → 26 % : les forêts ne se débloquent pas (il faut ~52 %).
2. **Atteindre la cible exige `rho ≈ 3100-3200`** — c.-à-d. une densité de croûte ≥ manteau
   (3300), **non physique**. À 3200 la terre < 500 m atteint 51 % (≈ Terre) : la mécanique
   FONCTIONNE, mais seulement hors borne. Cela CONFIRME que les cratons sont bien le levier
   (les descendre comble le déficit), mais que **la densité crustale n'est pas le bon
   véhicule** pour le faire.
3. **L'épaisseur n'est pas une alternative** : `craton_thickness_ratio = 1.25` est ancré au
   ratio crustal réel (40/32 km) et explicitement « NOT a knob » dans
   [init_r7 docstring](../../../../crates/ymir-core/src/tectonics_c1/init_r7/mod.rs) ; le
   baisser le rendrait non-cratonique.
4. **Effet de bord** : le balayage densité n'est pas parfaitement isolé au craton — la
   médiane non-craton dérive (429 → 476 m à 3000) via le percentile de sea-level partagé.
   Mineur mais l'invariant « non-craton intact » n'est pas strictement tenu.

## Cause racine (formulation)

L'isostasie d'Airy sur l'épaisseur CRUSTALE donne à une croûte cratonique épaisse une
freeboard haute. Les cratons réels ont une freeboard BASSE malgré une croûte épaisse, grâce
à une **racine lithosphérique mantellique froide et dense** (compensation isostatique
profonde — l'isostasie compositionnelle de Jordan, déjà citée dans le docstring init_r7
comme le raffinement de modèle prévu). `craton_rho_crust` (croûte plus dense) est un proxy
grossier de cet effet, MAIS il est plafonné par la densité crustale physique (≤3000) — donc
incapable de représenter la racine dense.

→ La hauteur worn-shield est un **RAFFINEMENT DE MODÈLE** (isostasie compositionnelle /
racine cratonique / worn-init — « 300 pas ≠ éons »), PAS un réglage de paramètre dans les
bornes. Conformément à la discipline, ce signal est remonté plutôt que forcé.

## Pistes (pour décision)

- **A. Raffinement de modèle (le vrai fix)** : modéliser la freeboard cratonique via une
  densité de COLONNE effective (croûte + racine dense) ou une isostasie compositionnelle —
  ce qui légitime un `rho` effectif > 3000 car il représente la colonne compensée, pas la
  croûte seule. Re-documenter/renommer le paramètre en conséquence (`craton_rho_effective`).
- **B. Worn-init** : initialiser les cratons déjà usés (S̃ abaissé représentant l'usure sur
  des éons que 300 pas ne simulent pas), au lieu d'épais-puis-résistant.
- **C. Accepter le partiel** : pousser `craton_rho_crust` à 3000 (borne) — gain marginal
  (1814 → 1197 m, terre < 500 m 25 → 26 %), ne débloque PAS les forêts. Probablement pas
  suffisant pour justifier le changement seul.
