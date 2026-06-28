# Verdict — routage des rivières : rectiligne sur le plat + rivières fantômes

**Probe :** `c1_closure_morphology::probe_river_routing_audit` (`#[ignore]`).
**Run :** 3 seeds (42, 1988, 2026), terrain + drainage cachés. Déclencheur : le fix bilan
hydrique a asséché les lacs sur-remplis (18 % → 0.7 %) → en retirant les plans d'eau qui
MASQUAIENT les plats, il a re-révélé deux défauts du routage.

## Mesure

| seed | river cells | sur plats pit-remplis | pas même-dir (droit) | cardinaux | navigable en déficit (précip<PE) | endorhéiques ignorés |
|---|---|---|---|---|---|---|
| 42 | 85 224 | 47 % | 54 % | 41 % | **92 %** | 113 |
| 1988 | 111 360 | 46 % | 54 % | 44 % | **81 %** | 186 |
| 2026 | 119 169 | 39 % | 53 % | 45 % | **79 %** | 84 |

## Partie A — motif RECTILIGNE (priorité 1)

**Cause (code-confirmée)** : `compute_d8` route en **D8 (8 directions discrètes)** ; sur les
**plats pit-remplis** (dépressions remplies au seuil par le priority-flood), il suit le
gradient **Garbrecht-Martz** (`resolve_flats`, un transform de distance vers l'exutoire).
D8 sur un champ de distance plat → toutes les cellules pointent dans la même direction
cardinale/diagonale vers l'exutoire → **canaux droits + affluents à ~90°**.

**Ampleur** : **39-47 % des cellules de rivière sont sur des plats pit-remplis** ; **53-54 %
des pas continuent la même direction D8** (longues droites), **~43 % cardinaux**. Le motif
est sur les plats — les plaines basses + les plateformes cratoniques (qu'on a créées) + les
bassins drainés. G-M avait corrigé l'ancien eps-fill (barres-parallèles/éventails 45°) mais
le D8 sur le plat reste cardinal-aligné.

**Approches (surfacées, pas choisies)** :
- **D∞ (Tarboton)** : direction de flux CONTINUE (pas 8 discrètes) → pas de verrouillage
  cardinal, angles naturels. Le plus direct sur le défaut, mais change le cœur du routage.
- **Micro-perturbation des plats** : ajouter un bruit léger au `filled` sur les plats → des
  micro-pentes réelles → D8 non dégénéré. Le moins invasif (local au pit-fill).
- **Méandrage stochastique** : post-traiter les tracés de canaux (perturber les polylignes).
  Cosmétique, n'attaque pas la cause.

## Partie B — rivières FANTÔMES (priorité 2)

**Cause (code-confirmée)** : `extract_rivers` tourne sur l'accumulation **GÉOMÉTRIQUE**
(priority-flood, compte de cellules) calculée AVANT le bilan hydrique. Elle ne connaît NI
les bassins endorhéiques (le bilan hydrique vit dans `c1_drainage`, après) NI le runoff
réel. Donc une rivière sort d'un bassin fermé vers la mer (fantôme), et la navigabilité
(cell-count) sur-estime les rivières partout.

**Ampleur** : **79-92 % des cellules de rivière NAVIGABLE sont en déficit hydrique
(précip < PE)** → ne portent aucune eau réelle. (Le monde à 45° est semi-aride : précip ~450
< PE ~854 presque partout → quasi pas de surplus → les rivières cell-count sont massivement
fantômes.) 84-186 bassins endorhéiques sont ignorés par le routage. `runoff_accum` (local au
bilan hydrique, sans reset endorhéique) n'aide pas encore.

## Verdict : DEUX fix SÉPARÉS (causes différentes)

- **A = tracé sur le plat** (D8 sur gradient plat → rectiligne) — fix de GÉOMÉTRIE du routage
  (D∞ / micro-perturbation / méandrage). Priorité 1 (visuel).
- **B = cohérence discharge + bassins fermés** (l'accumulation géométrique ignore le bilan
  hydrique) — fix de DISCHARGE : router/navigabilité sur le **runoff** (eau réelle) avec
  **reset aux bassins endorhéiques** → déserts sans rivières navigables, rivières mortes aux
  bassins fermés. Priorité 2. (Plus large que « juste l'endorhéique » : la navigabilité
  entière devrait dériver du runoff, pas du compte de cellules — 92 % de fantômes le montre.)

Ils peuvent partager un refactor du routage mais sont conceptuellement distincts. Diagnostic
seulement — le(s) fix suivent ce verdict.
