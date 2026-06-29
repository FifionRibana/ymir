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

---

## FIX peignes appliqué — la FRÉQUENCE relevée (pas une octave basse)

**Hypothèse de départ (le brief)** : une octave BASSE fréquence (λ ~ taille des grands plats)
donnerait au grand plat une variation à son échelle. **RÉFUTÉE par la mesure.** Une composante
basse fréquence est une **pente lisse à grande échelle** → les rivières la descendent tout droit
sur une longue distance → les peignes EMPIRENT : straightness grands plats 59 % → **72 %**, runs de
17+ pas 6 % → 33 % (`seed02026_residual_runlen.png` de cet essai : traînées diagonales plus longues).

**Ce que la donnée disait** : pour casser un long run droit en escalier court, le bruit doit
**basculer le choix latéral D8 plus SOUVENT** → il faut une fréquence **plus HAUTE** (longueur
d'onde plus courte), l'inverse de l'hypothèse. La basse fréquence flippe rarement (lisse) → runs
longs ; la haute fréquence flippe souvent → runs courts (escalier = acceptable).

**Le fix** : `frequency` 0.07 → **0.17** (λ ≈ 14 → 6 cellules), `amplitude` inchangée (0.45, sous la
borne anti-pit). Pas de composante basse fréquence (la machinerie multi-échelle a été retirée — elle
ne servait pas). Une seule manette : la fréquence, relevée.

**Mesure (seeds 42 / 2026), avant (λ14) → après (λ6)** :

| | λ14 (0.07) | λ6 (0.17) |
|---|---|---|
| runs LONGS (9+, peignes) | 28-29 % | **9-10 %** |
| runs 17+ (longues droites) | 5-6 % | **1 %** |
| straightness GRANDS plats | 57-59 % | **48 %** |
| straightness petits plats | 34 % | 34 % (inchangé) |
| straightness-plats globale OFF→ON | 60 → 44 % | 60 → **38 %** |

**Visuel** (`flat_tracing/seed02026_residual_runlen.png`) : les traînées diagonales rouges (peignes)
sont **dissoutes** — réseau quasi entièrement vert dendritique. Les petits plats sont inchangés
(la haute fréquence marche toujours). **Garde-fous** : aire de lacs 0.16 % = 0.16 % (pas de pit —
la fréquence ne touche pas l'amplitude/la borne), endorhéiques 113 ≈ 111, réseau connecté. Gated
(None byte-identique), plié dans la clé de cache via la config. Lib verte (472).

**Jugement de l'escalier diffus restant** : une fois les peignes partis, le résidu est l'escalier D8
court (90-91 % des pas droits en runs 2-8), qui rend un tracé **dendritique d'aspect naturel** (crops).
Les grands plats à 48 % et les petits à 34 % sont du même ordre — le motif n'est plus pathologique.
→ **On s'arrête là. D∞ reste en RÉSERVE** (chantier de fond : flux fractionnaire, touche
`compute_accumulation`/`extract_rivers`) — à n'engager que si l'escalier est un jour jugé gênant.

---

## DIAGNOSTIC directionnel — l'impression de lignes = quantisation DIAGONALE D8

L'utilisateur voit ENCORE des traînées diagonales (« plus courtes mais on a l'impression de longues
lignes »). La métrique des RUNS (longueur des segments droits) ne capte pas ça. Hypothèse : l'œil
lit le **parallélisme directionnel** (rivières voisines toutes orientées pareil), pas la longueur.
Probe `probe_direction_parallelism` (2 seeds, cache) : histogramme des 8 directions D8 + paramètre
d'ordre d'orientation LOCAL R2 + crop teinté (1 teinte/direction).

### Mesure

| classe | diagonal % (50 % = neutre) | maxbin | entropie /3 | R2 local (1=parallèle) |
|---|---|---|---|---|
| petits plats <100 c. | 38-41 % | 18-20 % | 2.94 | 0.40-0.41 |
| **GRANDS plats >5k c.** | **67-68 %** | 19-20 % | 2.87-2.90 | 0.39-0.40 |
| pente (non plat) | 46-47 % | 15-16 % | 2.97 | 0.45-0.51 |

### Lecture

1. **Biais DIAGONAL fort et SPÉCIFIQUE aux grands plats** : 67-68 % des pas y sont diagonaux, contre
   46-47 % sur les pentes (et ~50 % neutre). C'est le cœur de l'artefact. Sur une pente réelle, D8
   pénalise les diagonales (le `D8_DIST` √2 en pas 1) → diag < 50 %. Sur un plat, la passe 2 choisit
   le `flat_grad` voisin minimal **sans pénalité de distance** → les diagonales (qui « avancent » plus
   vite vers l'exutoire dans le champ de distance) sont sur-représentées. Un gradient continu pointant
   vers l'exutoire lointain est **snappé à 45°** → traînées diagonales parallèles de MÊME direction.

2. **Pas un peigne « serré » uniforme** : le R2 local (0.40 sur les grands plats) n'est PAS supérieur
   aux pentes (0.45-0.51) et l'entropie reste haute (2.87/3) — donc pas une seule direction partout.
   Mais le biais diagonal + la quantisation produisent les **streaks mono-teinte diagonaux** bien
   visibles sur `flat_tracing/seed02026_direction_hue.png` (grands plats), absents de
   `seed00042_direction_hue.png` (petits plats, dendritique multi-teinte). Le parallélisme est
   **concentré sur les grands plats, en diagonale** — exactement là où l'œil voit les lignes.

3. **D8-quantisé par construction** : les directions ne prennent que 8 valeurs. La perturbation ne
   change QUE *quelle* cellule bascule, jamais le fait qu'elle bascule vers une des 8 → elle ne peut
   pas étaler le biais diagonal. C'est pourquoi le tracé reste « en lignes » malgré la perturbation.

### Verdict — D∞ justifié (avec une piste plus cheap à tester d'abord)

L'impression de lignes vient bien de la **quantisation directionnelle D8** (biais diagonal sur les
grands plats), PAS de la longueur des runs (déjà optimisée). **D∞** (directions continues) est le
fix de fond justifié : il étale le biais diagonal sur un continuum d'angles → les streaks diagonaux
se dissolvent. Confirmé par la mesure (pas une supposition).

### Piste cheap ESSAYÉE puis RÉFUTÉE — normaliser la passe plat par la distance

L'hypothèse : une grande part du biais diagonal vient de la passe 2 de `compute_d8` qui compare
`flat_grad` **sans** diviser par `D8_DIST` (contrairement à la passe 1 des pentes). Fix tenté :
prendre le plus grand **drop / distance** (÷√2 pour les diagonales) dans la passe plat, gated.

**RÉFUTÉ par la mesure — surcorrection catastrophique** :

| classe | diag % avant | diag % après ÷dist | R2 avant | R2 après |
|---|---|---|---|---|
| GRANDS plats | 67-68 % | **2 %** | 0.39-0.40 | **0.90-0.91** |

Le biais bascule de diagonal à **cardinal** (diag 2 %) et le parallélisme EXPLOSE (R2 0.40 → **0.90**)
— des **peignes cardinaux pleins** bien pires que les diagonaux (`seed02026_direction_hue_NORMALIZED_WORSE.png` :
blocs mono-teinte, longues droites horizontales/verticales). En prime, ~3× plus de cellules de rivière
sur les grands plats (le routage cardinal collapse en quelques gros canaux).

**Pourquoi** : `flat_grad` est une **distance-graphe (BFS 8-connexe)**, PAS une hauteur euclidienne.
Un voisin cardinal ET un voisin diagonal vers l'exutoire ont TOUS DEUX `tl−1` → le MÊME drop
(`fhmax+1`). Diviser par la distance euclidienne √2 pénalise alors les diagonales à tort → cardinal
pur. L'analogie avec la passe 1 (qui, elle, opère sur de vraies hauteurs euclidiennes où ÷distance
est physiquement correct) NE TIENT PAS. Code reverté.

### Verdict final — D∞ est le seul fix de fond

La piste cheap est morte (la normalisation euclidienne ne s'applique pas à une distance-graphe). Le
biais diagonal n'est donc PAS un simple oubli de pénalité : c'est l'interaction **D8 × transformée de
distance** (quantisation à 8 directions d'un champ de distance octogonal). Seul **D∞** (direction de
flux CONTINUE, Tarboton) l'enlève — il étale le biais sur un continuum d'angles. Chantier de fond
(flux fractionnaire — touche `compute_accumulation`/`extract_rivers`/la cohérence runoff de B), à
n'engager que si l'utilisateur juge les streaks diagonaux résiduels gênants. **Rappel** : le tracé
est DÉJÀ dendritique et d'aspect naturel (perturbation + fréquence) ; D∞ est un polish, pas une
correction de défaut bloquant. Diagnostic seulement.

---

## FIX D∞ — ESSAYÉ puis REPLI D8 (le tracé rendu re-quantise)

Tentative (repli D8 gratuit, gating). D∞ (Tarboton) implémenté sur les plats : flux à direction
CONTINUE depuis le gradient de `flat_grad`, réparti fractionnairement entre les 2 voisins D8
encadrants ; `direction` = la primaire (plus grande fraction) → lacs/bassins/B gardent une mono
direction valide ; `accumulation` devient fractionnaire (eau conservée, `frac1 + (1−frac1) = 1`).
Probe `probe_dinf_compare` (seed 2026, grands plats), D8 vs D∞ (perturbation ON dans les deux).

**RÉSULTAT — D∞ est PIRE, même pathologie que la normalisation** :

| | D8 (actuel) | D∞ |
|---|---|---|
| diagonal grands plats | 67 % | **1 %** (bascule cardinal) |
| R2 local (1=parallèle) | 0.40 | **0.91** |
| aire de lacs | 0.14 % | 0.15 % (préservée) |
| endorhéiques | 86 | 85 (préservés) |

`seed02026_dinf_Dinf.png` : **blocs cardinaux pleins** (cyan/jaune/rouge), bien pires que le
diagonal de `seed02026_dinf_D8.png`. (L'hydrologie de B, elle, EST préservée — lacs, endorhéiques :
le garde-fou tient. Mais le tracé est pire.)

**Pourquoi — la cause est FONDAMENTALE** : le tracé RENDU suit la primaire (UN voisin D8 par
cellule) → il **re-quantise** quoi qu'il arrive. Et le champ `flat_grad` (distance BFS 8-connexe ≈
Chebyshev) a un gradient **CARDINAL** → l'angle continu pointe cardinal → la primaire se verrouille
en cardinal → peignes cardinaux. C'est la MÊME bascule que la normalisation par distance, par la même
raison (la géométrie du champ de distance). Le bénéfice continu de D∞ est dans l'ACCUMULATION
(fractionnaire), pas dans le TRACÉ (qui reste une suite de cellules 8-connexes).

**Conclusion** : le parallélisme rendu est un **plancher du dessin des rivières comme cellules de
grille** sur un grand plat à gradient uniforme. AUCUN réglage du routage (perturbation, normalisation,
D∞) ne l'enlève — le D8 actuel (diagonal, R2 0.40) est le MOINS parallèle des options testées. Le
dissoudre demanderait un rendu VECTORIEL/sous-pixel des rivières (un changement de RENDU viz, pas de
routage), hors périmètre.

**Décision : REPLI D8** (assumé, pas un échec — le tracé D8 + perturbation + fréquence est déjà
dendritique et naturel). D∞ reste une **option gated documentée** (`dinf`, défaut `false`, mesurée
inférieure) pour reproductibilité ; D8 jamais cassé (byte-identique, lib verte 472). Fin du fix A :
les leviers de routage sont épuisés ; tout gain supplémentaire est côté rendu, pas routage.
