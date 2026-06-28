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

---

## FIX B appliqué — débit sur le runoff réel + bassins fermés en puits

**Quoi** : la navigabilité (et l'aire de drainage par segment) ne lisent plus le **compte de
cellules** géométrique mais le **DÉBIT runoff réel** (`runoff_accumulation` : `max(0, précip−PE)·
aire`, accumulé en aval le long du D8). Les **bassins endorhéiques** (du bilan hydrique) sont des
**PUITS** : le runoff y meurt (évaporation) → aucun débit fantôme ne continue vers la mer en aval
d'un lac fermé. L'aire effective d'un segment = `débit / runoff de référence (300 mm)`, donc les
seuils de navigabilité en km² s'appliquent inchangés à la profondeur de runoff de référence (une
rivière humide classe exactement comme l'ancien cell-count ; plus sèche → déclassée).

**Ordre** : la navigabilité est calculée APRÈS les lacs (il faut connaître les bassins
endorhéiques). `extract_rivers` garde les **tracés** géométriques (le défaut A, rectiligne,
n'est pas touché ici) ; seul le **débit/navigabilité** porté par ces tracés change.

**Allochtones préservés** (la nuance critique) : le débit est ACCUMULÉ depuis l'amont. Une rivière
née dans un massif humide qui traverse ensuite un terrain localement aride GARDE son débit (le
Nil) → reste navigable malgré le déficit local. On ne coupe jamais sur le déficit LOCAL.

**Mesure (3 seeds, 42/1988/2026), avant → après** :

| seed | navigable en déficit (avant) | navigable en déficit (après) | navigable / cellules rivière (après) |
|---|---|---|---|
| 42   | **92 %** | **48 %** | 4 % (3 729 cell.) |
| 1988 | **81 %** | **38 %** | 3 % (2 934 cell.) |
| 2026 | **79 %** | **43 %** | 2 % (2 917 cell.) |

Les fantômes purs (eau inexistante) sont tués : la navigabilité passe à **2-4 %** des cellules de
rivière (monde semi-aride à 45° → navigable seulement là où un bassin humide produit du surplus).
Le **résidu ~40 %** « navigable en déficit local » n'est PAS un fantôme : ce sont les rivières
**ALLOCHTONES** (débit accumulé amont traversant un terrain aride) — exactement ce qu'on doit
garder. La fraction de cellules sur plats (47/46/39 %) et la rectilinéarité (53-54 %) sont
inchangées : c'est le périmètre du fix A, distinct.

**Gating** : `c1_drainage` sans `climate` (None) → chemin géométrique byte-identique (logique
inchangée). `ALGO_DRAINAGE` 1→2 (le code a changé sans changer d'input → le cache se ré-adresse).
Lib verte (472). `runoff_accumulation` est factorisé et partagé avec le bilan hydrique des lacs.

**Reste** : fix A (tracé rectiligne sur les plats) — géométrie du routage, à traiter séparément.

---

## FIX A appliqué — micro-perturbation du gradient G-M sur les plats

**Cause rappel** : le pit-fill rend les dépressions PARFAITEMENT plates ; Garbrecht-Martz y impose
un gradient distance-à-l'exutoire UNIFORME → D8 route en peignes droits cardinaux/diagonaux.

**Le fix (la version retenue après itération)** : on ajoute un **bruit cohérent (value-noise)** au
**gradient G-M `flat_grad`** sur les plats — PAS à la surface `filled`. Le bruit est **borné à
`amplitude·(fhmax+1)` avec `amplitude < 0.5`** (où `fhmax+1` = le poids d'un pas `tl` du gradient).

**Pourquoi cette borne = pas de pit, garanti mathématiquement** : un pas vers l'exutoire fait baisser
`flat_grad` de `(fhmax+1)` ; l'écart de bruit entre deux cellules est `< 2·0.5·(fhmax+1) = (fhmax+1)`,
donc le voisin-vers-l'exutoire reste **strictement plus bas** après bruit → drainage garanti, aucun
minimum intérieur, réseau connecté. Sous cette borne le bruit domine le choix LATÉRAL (le tie-break
`fh` et les petits écarts `tl`) → la rivière **serpente**. `filled` n'est jamais touché → niveaux de
lac, bilan hydrique, bassins endorhéiques **inchangés** (B intact).

**Première tentative écartée** : perturber la surface `filled` puis re-pit-filler. Le re-fill
**re-aplatissait** les pits creusés par le bruit → re-créait des plats → G-M cardinal (le cardinal
remontait même). Le re-fill contrecarrait la perturbation. La perturbation du gradient (ci-dessus)
n'a pas ce problème : pas de re-fill, pas de pit par construction.

**Itération visuelle** (probe `probe_flat_tracing_compare`, crop des plats, OFF vs ON) :
`flat_tracing/seed00042_tracing_OFF.png` montre les **peignes diagonaux parallèles** (l'artefact) ;
`_ON.png` montre des **réseaux dendritiques organiques** (arborescences naturelles, peignes dissous).
Réglage retenu : `amplitude = 0.45` (juste sous la borne 0.5 → serpentement fort), `frequency = 0.07`
(longueur d'onde ≈ 14 cellules), `octaves = 4`.

**Mesure (seed 42)** :

| | valeur |
|---|---|
| straightness SUR LES PLATS (pas même-dir) | **OFF 60 % → ON 44 %** |
| aire de lacs (invariant) | OFF 0.16 % = ON 0.16 % ✓ |
| lacs endorhéiques (invariant) | OFF 113 ≈ ON 112 ✓ |

L'audit global (3 seeds) : la rectilinéarité totale 53-54 % → 46-47 % (diluée — 60 % des cellules de
rivière sont sur de vraies pentes, déjà bien routées) ; endorhéiques 113/186/84 → 112/185/85
(inchangés) ; navigabilité (fix B) préservée.

**Gating** : `flat_perturbation: None` (sur `FlowConfig` / `C1DrainageConfig`) → routage légataire
byte-identique. Plié dans la clé de cache via la config (pas de re-bump ALGO). Plats SEULEMENT (les
pentes passent par la passe 1 de `compute_d8`, intouchées). Lib verte (472 ; l'échec
`rectangular_simulation_smoke_test` est pré-existant, tectonique, sans rapport).

**Verdict** : A + B faits. Le réseau hydrographique est désormais physique (débit réel, bassins
fermés respectés, allochtones préservés) ET son tracé sur les plats est naturel (dendritique).

---

## DIAGNOSTIC du rectiligne RÉSIDUEL (fix A suite)

La perturbation a amélioré le tracé (straightness-plats 60 % → 44 %) mais **insuffisamment** : il
reste des rivières droites. 44 % est une MOYENNE — ce diagnostic localise et caractérise le résidu
pour choisir le bon fix (D∞ de fond vs réglage de perturbation). Probe `probe_residual_rectilinear`
(2 seeds, cache, perturbation ON) : longueur des runs de pas droits, straightness par TAILLE de plat
et par DÉBIT, + crops colorés par longueur de run (vert = run court, rouge = run long).

### Mesure

| | seed 42 | seed 2026 |
|---|---|---|
| pas droits en runs COURTS (2-8, escalier) | **71 %** | **72 %** |
| pas droits en runs LONGS (9+, peigne) | **29 %** | **28 %** |
| straightness petits plats (<100 c.) | 34 % | 34 % |
| straightness GRANDS plats (>5k c.) | **59 %** | **57 %** |
| straightness par débit (<50 / 50-1k / >1k km²) | 49 / 42 / 47 % | 49 / 41 / 18 % |

### Le résidu a DEUX natures distinctes

1. **DOMINANT (~71 %) — escalier D8 DIFFUS.** La majorité des pas droits sont dans des runs COURTS
   (2-8 pas) : la rivière serpente mais D8 (8 directions discrètes) la rend en marches d'escalier.
   C'est le **plancher de discrétisation**, partout (le vert dendritique des crops). Visuellement
   plutôt naturel. Seul **D∞** (direction de flux continue) l'enlève.

2. **SECONDAIRE (~28 %) — peignes CONCENTRÉS sur les GRANDS plats.** Les runs LONGS (9+) sont des
   **traînées diagonales parallèles** (les lignes rouges des crops, surtout `seed02026_residual_runlen.png`),
   concentrées sur les GRANDS plats : straightness **59 % sur les grands plats vs 34 % sur les
   petits**. Cause : la fréquence du bruit (0.07, longueur d'onde ≈ 14 cellules) est trop FINE pour
   les grands plats — sur des centaines de cellules le bruit **moyenne à ~0** et le gradient G-M
   uniforme (diagonal vers l'exutoire lointain) **redomine** → longs peignes diagonaux. L'amplitude
   (0.45) ne peut pas monter (borne 0.5 de garantie anti-pit).

   Ce n'est **PAS** lié au DÉBIT (la straightness ne croît pas avec le débit ; le bucket >1k km² est
   trop petit — 67-187 pas — pour conclure). Le driver est la **TAILLE du plat**, pas le débit.

### Verdict — les deux fix sont distincts, le réglage d'abord

- **Les peignes (rouge, ~28 %, le plus VISIBLE — longues droites parallèles)** = perturbation MAL
  CALÉE pour les grands plats (fréquence trop fine). Fix = **régler la perturbation** : abaisser la
  fréquence (longueur d'onde plus grande) ou ajouter une octave BASSE fréquence qui varie à travers
  le grand plat. **Cheap, ciblé, attaque l'artefact le plus visible. À faire en premier.**
- **L'escalier (vert, ~71 %, diffus)** = plancher de discrétisation D8. Fix = **D∞** (mono → flux
  fractionnaire, touche `compute_accumulation` / `extract_rivers` — chantier de FOND). À n'engager
  que si l'escalier reste jugé trop rectiligne UNE FOIS les peignes dissous.

Donc : **réglage de la perturbation d'abord** (les peignes des grands plats), **D∞ ensuite si
nécessaire** (l'escalier diffus). Diagnostic seulement — le fix suit.
