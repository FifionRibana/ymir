# Ymir — feuille de route : du bruit aux closures

Document de décision, arrêté à la clôture de la phase hydrologie.

---

## Le constat qui structure la suite

Le diagnostic des dépressions a désigné un coupable unique et inattendu : **le FBM crée
90 682 cuvettes fermées** là où la tectonique en produit 16. Ni MFD, ni talus, ni
diffusion n'en créent — tous les trois en réduisent le compte. Et un curseur de maturité
ne suffit pas : 16 itérations d'incision ne descendent qu'à ~45 000, soit encore 3 000×
le compte tectonique, au prix d'un rabotage des chenaux qu'on avait délibérément borné.

Ce chiffre éclaire une question plus large, posée en fin de phase : **les plaques ne se
sont pas solidifiées avec des irrégularités qui persisteraient.** La croûte océanique a
moins de 200 millions d'années, les cratons ont été érodés et réenfouis plusieurs fois.
Aucune rugosité originelle n'a survécu ; le relief est **entièrement** le produit de
processus actifs. Il n'y a donc rien de « donné » à représenter par du bruit.

Le FBM a été introduit pour une raison technique légitime — combler les 128× d'upscale
entre le champ grossier 64² et le HD — **pas** pour représenter un phénomène. Il tient la
place d'une physique manquante, et c'est pourquoi il fabrique des artefacts.

**La démarche par closures est donc cohérente ; c'est le FBM qui est l'intrus.**

## La trajectoire

Deux lectures du même chantier, l'une immédiate, l'autre de fond :

- **court terme** — empêcher le FBM de créer des dépressions ;
- **fond** — **ajouter de la STRUCTURE CAUSALE aux échelles ≥ cellule grossière**, ce
  qu'aucune closure ne faisait avant. Formulation corrigée en C-3b : « remplacer
  progressivement le FBM » était inatteignable — aucune closure ne détient d'information
  SOUS la cellule grossière (~6 km), donc aucune ne peut le remplacer dans sa propre bande
  (~6 km → 49 m), que seule l'érosion comble à partir d'une graine de brisure de symétrie.
  Ce qui reste vrai et acquis : les closures retirent au FBM la charge de représenter une
  physique aux échelles qu'elles couvrent, et l'empêchent de fabriquer des artefacts (C-1).
  Le FBM ne disparaît pas ; il redevient une **graine de détail sous-grossier**, pas le
  porteur d'une structure qu'il n'a jamais eu les moyens de représenter. (Ce critère vaut
  aussi pour C-4 : l'érosion côtière sera jugée sur l'ajout de structure littorale causale,
  pas sur un rétrécissement du plancher FBM qu'elle ne peut pas produire non plus.)

Règle qui accompagne cette trajectoire : à mesure que les closures arrivent, **le FBM doit
rétrécir, pas cohabiter**. Deux sources de détail — une causale, une aléatoire — qui se
superposent sans se coordonner reproduiraient le problème des stries, qui a coûté
plusieurs passes à diagnostiquer.

**Correction (établie en C-3b) — jusqu'où le FBM peut rétrécir.** Le FBM comble le 128×
d'upscale : des longueurs d'onde de la cellule grossière (~6 km) à la cellule HD (~49 m).
**Aucune closure tectonique ne détient d'information SOUS la cellule grossière.** Les
closures ajoutent donc de la structure aux échelles ≥ cellule grossière ; elles **ne
peuvent pas remplacer le bruit DANS SA PROPRE BANDE** (le détail sous-grossier ne peut venir
que de l'érosion, qui exige une graine de brisure de symétrie — le plancher de dégénérescence
de C-1). Faire *baisser le plancher d'amplitude FBM* est donc **inatteignable par principe**,
pour C-3b comme pour C-4. Ce que les closures font, et doivent faire, c'est empêcher le FBM
de fabriquer des artefacts (C-1) et lui retirer la charge de représenter une physique aux
échelles qu'elles couvrent — pas le supprimer sous la cellule grossière. Consigné à l'ADR
(section « C-3b »).

## Ce qui a été écarté, et pourquoi

**Le remplissage de cuvettes à seuil.** Peu coûteux (un flood), il se compose proprement
avec le breach et donne un cadran direct sur le nombre de lacs. Écarté parce qu'il
**masque** un bruit non physique au lieu de le supprimer, et que la carte actuelle est
acceptable — donc le palliatif n'est pas nécessaire. Décision et raison consignées à
l'ADR, car c'est le raccourci évident que quelqu'un retentera.

**La saisonnalité.** Hors périmètre pour l'instant. Conséquence assumée : un bassin
sahélien ressortira en **lac permanent peu profond** plutôt qu'en playa. Le niveau moyen
est mathématiquement correct ; c'est la dimension temporelle qui manque. Note pour plus
tard : un bassin endoréique dont l'apport est très inférieur à l'évaporation potentielle
**est** une playa — il suffirait de l'étiqueter, sans simuler le cycle.

## Les quatre closures, dans l'ordre retenu

### 1. Conditionnement du FBM — ✅ FAIT (voir ADR 0001, section « C-1 »)

**Livré :** budget de relief `flow_conditioning = β` (plafond d'amplitude ∝ pente locale),
associé à un étirement aval fixe ×8 du bruit (features allongées le long de l'écoulement,
jamais comprimées → pas d'aliasing). β = 0.1 en production. Cuvettes post-FBM : 90 682 → 6 999 à
8192², 6 070 → 220 à 2048² ; les cuvettes profondes (≥ 50 m) tombent à l'ordre tectonique
(23) ; morphologie de montagne préservée (pentes > 30°/> 45° tenues ou accentuées).
Recherche littérature : aucune formulation nommée n'existe (domain warp, bruit anisotrope,
retrait de minima en post-traitement sont voisins mais résolvent un autre problème) — la
dérivation est consignée à l'ADR. Le résidu sub-métrique et les artefacts de l'incision
restent le travail de fond ci-dessous.

**Pourquoi en premier :** chaque closure suivante sera calibrée **sur le champ que le FBM
produit**. Le modifier après obligerait à tout recalibrer — la même leçon que l'échelle
métrique posée avant les seuils physiques.

**Principe :** un bruit respecte l'écoulement s'il **ne crée pas de minimum local**.
Techniquement, moduler la perturbation par la direction du gradient plutôt que de
l'ajouter isotropiquement — un bruit qui déplace latéralement le lit sans créer de
contre-pente. Proche de `amplitude_slope_factor`, mais avec la **monotonie** comme
critère, pas un simple facteur d'amplitude.

**Critère de succès, déjà mesuré :** le compte de cuvettes après FBM doit rester du même
ordre que les **16** tectoniques. Net, falsifiable, référence en main.

### 2. Volcanisme — ✅ FAIT (voir ADR 0001, section « C-2 »)

**Livré :** édifices dérivés de la tectonique C1 (arcs insulaires sur marge de subduction,
chaînes de points chauds d'âge croissant le long du mouvement de plaque, rifts) ; géométrie
ancrée (Wood 1978 strato `H=0.122·Wb+0.45`, Grosse & Kervyn 2018 boucliers H/Wb 0.01–0.1,
Wb ≥ 2 km, cratère D/W ≤ 0.25) ; injection HD après FBM, avant érosion, avec reconstruction
du rebord des cratères **actifs** (construction ≥ érosion) protégée du breach relief-v3 ;
`LakeType::CraterAcidic` (Varekamp 2000, pH < 2 sur cratère actif), portée à l'export,
détectée par une passe dédiée contournant le plancher générique de 5 km². Confirmé **sur
l'export** : minorité climat-dépendante de lacs de cratère acides (0 en désert chaud → 2–3
en humide/tropical, échelle Kawah Ijen/Pavin), pas « partout ». Marqueurs actif/éteint dans
la viz (toggle Symboles). `CraterNeutral` inatteignable par construction (documenté).

**Pourquoi tôt :** bien délimité — il ajoute du relief construit **sans toucher au
drainage**, donc il ne rouvre pas la chaîne hydrologique stabilisée sur une dizaine de
passes. Risque de régression faible.

**Ce qu'il apporte :** des lacs de cratère à eau acide, un type de plan d'eau qu'aucune
autre closure ne produira. Il complète `lake_type` (aujourd'hui exoréique/endoréique)
avec une troisième nature **dont la chimie découle du mécanisme**, pas d'un réglage.

**Emplacement causal plutôt qu'aléatoire :** plaques et zones de subduction existent déjà
dans C1 — arcs insulaires le long des subductions, points chauds intraplaques, rifts.
C'est ce qui distingue un volcanisme structuré d'un semis de cônes.

**Point de vigilance :** un édifice volcanique **crée une dépression fermée** — le
cratère. Après le conditionnement du FBM, une cuvette sera légitime ou artefactuelle ;
celles du volcanisme sont **volontaires** et doivent être reconnues comme telles, non
filtrées avec le bruit. Faire le conditionnement d'abord met ce critère en place avant
d'en avoir besoin.

### 3. Hétérogénéité lithologique — ✅ FAIT (voir ADR 0001, section « C-3 »)

**La seule des quatre qui REMPLACE du bruit par de la structure.** Volcanisme et érosion
côtière *ajoutent* du relief ; la lithologie change la façon dont **tout le reste**
s'érode. C'est donc elle qui fera le plus rétrécir le FBM — et qui validera le recadrage
jusqu'au bout.

**Ce qu'elle apporte :** la variété de vallées demandée — gorges dans la roche dure,
vallées larges dans la tendre —, les crêtes et corniches que le FBM ne sait pas
structurer, et le placement des minerais côté Living Landz.

**État :** livré. `K` par cellule = multiplicateur d'érodabilité CAUSAL (jamais du bruit) :
socle dur ×1 (référence relief-v3), rift tendre ×10, footprints volcaniclastiques ×3 —
tous des signaux déjà présents dans l'état tectonique (rift `age=0`, édifices C-2), donc
aucun nouveau champ advecté. Le spread a été MESURÉ, pas prédit : balayage ×3/×10/×30/×100
sur la chaîne de production complète, deux résolutions
(`tests/c3_lithology_sweep.rs`, rapport `docs/reports/.../closure_morphology/c3_lithology_sweep.md`).
Deux effets séparés par design : socle dur = référence intacte (colonnes HARD plates sur
tout le balayage → pas de ralentissement global), seul le contraste bouge. C-1 survit
(dépressions 982→977 à 2048², 17516→17114 à 8192²). Désactivé par défaut → production
byte-identique. Toggle viz « Lithologie (C-3) » (mode Expert, avec la stream-power).
Reste : la validation visuelle de l'auteur sur l'export, puis marquer ✅ FAIT.

La déposition manque (érosion production détachment-limited) → pas de bassins sédimentaires
causaux ; enregistré comme limitation dans l'ADR, non comblé par du bruit ou de la géométrie.

### 3b. Structure héritée (fracturation) — ✅ FAIT (densité seule ; voir ADR 0001, section « C-3b »)

**Pourquoi.** C-3 a établi que le socle est lithologiquement uniforme et dur — de la physique,
pas un manque. La structure d'un socle mûr n'est donc pas lithologique mais **tectonique** :
la même roche, mais **découpée** par des fractures. La densité de fracturation contrôle
l'érodabilité (Molnar 2007 : la tectonique érode surtout en fracturant → plucking).

**Ce qui est livré — la densité seule, isotrope.** `K = 1 + amplitude·densité`, densité =
proximité aux contacts tectoniques (frontières **convergentes + transformantes** de la
classification dynamique — jamais `cratonic_mask` ni le placeholder géométrique). Le craton
intact **émerge** loin des contacts à `K=1` (référence, ralentissement global nul par
construction) ; les ceintures orogéniques fracturées s'érodent/disséquent plus. Mesuré sur la
chaîne complète, deux résolutions : à 8192² (l'export) le craton tient exactement (105→104 m)
et la ceinture gagne de la dissection (398→538 m à ×16), **C-1 s'améliore** (pits 17382→15639).
Défauts : amplitude ×4, ceinture étroite (decay 25 km, craton = 53 % majoritaire). Gated OFF →
byte-identique.

**Ce qui a été écarté — l'orientation, mesurée puis abandonnée.** Le closure directionnel
(vallées alignées sur la fabrique via incision anisotrope) a été **implémenté et mesuré** : il
ne marche pas. Le taux d'incision ne peut pas réorienter un récepteur figé par la topo
(alignement 0.639→0.536 à relief préservé, pits 1001→2146). Le verrou est le **routage**, mais
son rayon d'impact couvre C-1/rivières/lacs et le champ directionnel de C1 est trop pauvre
(vitesses constantes par plaque, pas de strain, pas d'histoire) → grain uniforme = artefact.
Et la premisse est faible (le treillis appalachien = plis stratifiés, hors-scope). Limitation
caractérisée + spécifiée à l'ADR.

**Sutures fossiles** différées avec spécification (l'accrétion n'enregistre pas *où* elle soude ;
la densité-frontières capte déjà les collisions actives).

Reste : validation visuelle de l'auteur sur l'export, puis ✅ FAIT.

### H-1 / H-2. Chaîne hydro : infiltration puis incision de seuil — 🔬 EN COURS (AVANT C-4)

**Pourquoi elles passent AVANT l'érosion côtière (C-4).** L'incision de seuil (H-2) vidange
les bassins exoréiques : l'exutoire creuse son col, le lac se retire, l'eau ressort en
rivière → **le trait de côte bouge** (embouchures, deltas, niveau de base littoral). Faire
C-4 avant reviendrait à sculpter un littoral que H-2 déplacerait ensuite — on reprendrait le
travail. Diagnostic à l'appui (voir `closure_morphology/comb_and_lake_diagnosis.md`) : 92
lacs, 91 exoréiques, 0 endoréique, 6220 km² d'eau, les 10 plus grands tenant 74 % de l'aire —
des bassins remplis au col, pas des mares de bruit. La chaîne hydro doit d'abord assécher ce
qui doit l'être, PUIS on sculpte la côte qui en résulte.

**H-1 — infiltration (la plus légère, en premier).** Le bilan lit `runoff = max(0, precip −
PE)` : tout le surplus devient ruissellement de surface. En réalité une fraction **s'infiltre**
et n'atteint jamais un lac par écoulement de surface. Ymir n'a **aucune eau souterraine**. H-1
ajoute ce premier terme : un champ de perméabilité CAUSAL (classes lithologiques C-3 + densité
de fracturation C-3b — la roche fracturée/volcaniclastique infiltre plus que le socle
cristallin intact ; jamais du bruit) module la fraction infiltrée. Rayon d'impact **petit** :
change le BILAN, pas la géométrie ni le routage. Des bassins basculent en endoréique, des
niveaux baissent — un pas correct qui ne casse rien. Gated OFF → byte-identique.

**H-1 — ce qu'elle a RÉELLEMENT livré (voir ADR « H-1 »).** La mesure a invalidé la prémisse.
Le vrai défaut était que **les lacs de surface ne voyaient JAMAIS le bilan hydrique** : le
pré-breach tourne avec `climate = None` et sa classification purement géométrique écrase
celle calculée avec le climat. Finding 39 ne valait donc que pour les bassins **sous-marins**.
Correctif livré **par défaut** (ce n'est pas une option : une classification fausse ne se
gate pas) → **53 à 100 % des lacs de surface sont endoréiques** une fois le bilan appliqué,
dont 45 % en humide. ⚠️ change `lake_type`, donc biomes et rendu, sur toute carte existante.
L'**infiltration** elle-même est un terme **secondaire mesuré** (+0 à 3 lacs), conservé car
physiquement juste et **gated OFF**.

**H-2 — incision de seuil : périmètre EFFONDRÉ par H-1.** Un lac exoréique a un exutoire →
son seuil s'inciserait et le vidangerait, mais le modèle le gèle. Or H-1 a montré que le
reliquat exoréique n'est que de **3 à 27 lacs selon le climat, pas 91**, et que les
endoréiques se rétractent d'eux-mêmes (mesuré 1798 → 483 km² en tropical, **rapporté et non
appliqué**). Le rayon d'impact est donc bien plus petit qu'anticipé — **résultat direct de
l'ordre H-1 puis H-2**.

**⚠️ CORRECTION AU CADRAN DE H-2 — le « cadran temporel » France↔Écosse ne peut pas exister
comme durée** (ADR Finding 44, mesuré). `K` et la durée ne sont **pas observables séparément**
dans le modèle : à `k_time = K·dt·iterations` fixé, choisir n'importe quel `T` donne
`K = k_time/T` et `dt = T/iterations`, donc `k·dt = k_time/iterations` — tout ce que lit
l'incision est inchangé et **`T` s'annule**. 10⁴ ans avec un `K` cent fois plus grand produit
**exactement le même terrain**. Épinglé par `duration_cancels_out_of_the_incision` à
T = 10⁴, 10⁶ et 10⁸ ans.

Conséquence, énoncée **visiblement** parce qu'elle invalide la conception initiale du cadran
(portée sur plusieurs tours, y compris le nom de la branche) : **le contrôle exposé au joueur
n'est pas « combien d'années » mais « combien d'érosion cumulée »** — un produit
érodabilité×temps intégré. Traduire `k_time` en années exige de fixer `K` indépendamment, ce
que fournissent les multiplicateurs par lithologie de C-3 plus **un** ancrage absolu Stock &
Montgomery (ordre 10⁷–10⁸ ans à érodabilité de roche dure, réserve m = 0,4 vs m = 0,5 énoncée).

**Ce que ça ne détruit pas :** l'incision de seuil reste le **seul** mécanisme qui vide les
lacs, et le périmètre de H-2 est intact. Seule la **nature du paramètre** change.

**Et une question structurelle ouverte, plus grande qu'un calibrage — à ne PAS ouvrir avec
H-2** (ADR Finding 44) : le nombre de Courant de la config livrée est **1353 à 2048² et 3699 à
8192²**. Le schéma implicite est *stable* à n'importe quel Courant, mais la stabilité n'est pas
l'intégration : avec `f = K·dt·A^m/dist_m ≫ 1`, la mise à jour pose chaque cellule sur la
hauteur de son récepteur en **un** pas. Le terrain livré est donc **le point fixe de deux
relaxations locales, pas le résultat d'un épisode d'érosion** — le modèle n'intègre aucun temps,
à aucune des deux résolutions. C'est une façon légitime de fabriquer du terrain ; ce n'est pas
une façon d'exprimer une durée. Le compte de pas que la borne CFL exige serait de 2706 à 2048²
et 7399 à 8192² (contre 2 livrés), soit ~7400 recalculs de flux sur 67 M cellules : dériver le
compte honnêtement ne corrige pas le modèle, ça révèle qu'il ne fait pas ce qu'on croyait.

**Le critère de jugement de H-2 est prêt** (et il a fallu le réparer d'abord) : la métrique
naïve de plaine contiguë n'était pas une propriété du continent mais du pas de grille
(×4.97 entre 2048² et 8192²). La métrique retenue — **fermeture morphologique à 200 m de
distance de pontage physique, puis composantes connexes** — converge à **×1.03**. La variante
« au pas hexagonal » a été **écartée par la mesure** (×3.25). H-2 se jugera sur la plus grande
plaine et l'aire totale de plaine (3024/3129 km² et 10643/9989 km² aujourd'hui), les deux
seuls agrégats qui convergent.

### 4. Érosion côtière — APRÈS H-1/H-2

**Pourquoi en dernier :** elle travaille sur un littoral dont la forme dépend de tout
l'amont — y compris de H-2, qui déplace le trait de côte (voir ci-dessus). La faire après
lithologie, volcanisme et la chaîne hydro lui donne un trait de côte déjà
structuré à sculpter, plutôt qu'à reprendre ensuite.

**Ce qu'elle apporte :** falaises, plages, plateformes d'abrasion — la morphologie
littorale qui a préoccupé plusieurs fois pendant la phase hydrologie (embouchures
perchées, falaise contre plage au rendu).

**Note :** `erosion/coastal.rs` est un fichier vide marqué M5, comme `thermal.rs`,
`aeolian.rs` et `glacial.rs`.


## ⛔ CRITÈRE D'ACCEPTATION DES FRANGES CÔTIÈRES — opposable, écrit avant mesure (ADR Finding 76)

L'objectif n'a jamais été « moins de franges ». Il est : **la chaîne upscale + closures +
incision ne produit AUCUNE frange côtière.** Mesuré à 8192², sur le champ conditionné, avec
`coast_shape_thresholds` aux deux échelles (≥ 1 km et ≥ 2 cellules) :

⚠️ **Le seuil « ≤ 2× RÉFÉRENCE » est RETIRÉ (Finding 77).** Le critère est une **différence
contre le champ PRÉ-INCISION de la même seed** — l'autorité, mesurée au F76-A5 : 20 / 8 / 1 602 km,
indiscernable de la RÉFÉRENCE. **L'incision ne doit ajouter aucun éperon.**

| ligne | cible | SHIPPED | Δ | LAW ON |
|---|---|---|---|---|
| éperons ≥ 1 km | **Δ = 0** | 1 848 | **+1 828** | +896 |
| éperons ≥ 2 cellules | **Δ = 0** | 1 831 | **+1 823** | +1 048 |
| km de côte | Δ sous le bruit u16, **±3 km (0,2 %)** | 6 206 | **+4 604 (+287 %)** | +2 038 |
| spectre | illisible, comme celui du pré-incision (8 éperons ⇒ 0 écart) | λ 14 à ×3,74 | **lisible ⇒ ÉCHEC** | — |
| et l'intérieur ne bouge pas | hypsométrie ±1 m, part chenal ±1 pt, dépressions ±2 % | — | — | — |

Les vingt éperons du pré-incision sont la géométrie tectonique, pas des franges.

> ⛔ **PRÉALABLE NOMMÉ (Finding 77)** — **il n'y a aucune eau entre 0 et −20 m sur toute la carte**,
> avant comme après incision : le maximum des 56 077 791 cellules océan vaut exactement
> **−20,000 m**, et les 120 807 cellules océan bordant la terre y sont toutes. Chaque trait de côte
> est une **marche verticale de 20 m**, sans estran ni pente de plateau. L'isoligne 0 m ne court
> donc jamais côté mer — le côté mer est une constante — et la côte est entièrement un
> franchissement de seuil du micro-relief terrestre contre un plancher parfaitement plat.
> **LOCALISÉE (Finding 78)** : `BathymetryProfile::shelf_min_depth_m = 20.0`
> (`terrain/bathymetry.rs:49`), clampée ligne 178, posée par le commit `563198e` du 2026-06-26
> qui **ne la justifie nulle part** — et `bathymetry.rs:171` prend `.max(1.0)`, donc **1 m suffit
> à l'invariant énoncé** (« le masque terre/mer ne bouge pas »). C'est un **PROXY**.
>
> ⛔ **DEUX ITEMS, à ne plus confondre :**
>
> **(1) La falaise de 20 m est un défaut consommateur à part entière.** Marche terre−mer p50
> **22,00 m** ; **0,0 %** du trait de côte est sous 15 % de pente et **100,0 %** est au-dessus de
> 25 %. Berges, plages et quais ne sont pas difficiles : ils sont **impossibles**. À
> `shelf_min_depth_m = 1`, la marche p50 tombe à **3,00 m** et **71,4 %** du trait passe sous 15 %
> — **avec le continent bit-identique** (masque 0 cellule déplacée, hypsométrie 685,32 m, puits,
> travail d'érosion : identiques au dernier chiffre). **Ce point ne dépend pas de la question des
> franges.**
>
> **(2) Le clamp CACHE 43 % de la fourrure.** Éperons ≥ 2 cellules : 1 291 (200 m) → **1 831
> (20 m livré)** → 2 705 (5 m) → **3 190 (1 m)**. Marching-squares croise à
> `t = (0,5 − h_mer)/(h_terre − h_mer)` : une mer plus profonde est un dénominateur plus grand,
> donc un tracé moins sensible au micro-relief terrestre. **Le compte d'éperons à l'échelle
> cellulaire n'est pas une propriété de la terre seule.** Idem pour λ : 28 → 14 → 15 → **6**, et le
> pic passe sous la ligne des 3× à 1 m. **Baisser le plancher AGGRAVE la mesure avant que tout
> remède ne l'améliore ; les deux items ne se jugent pas sur un seul chiffre.**
>
> ✅ **CORRIGÉ EN PRODUCTION (Finding 79)** : `shelf_min_depth_m` **20 → 1,0 m**,
> `ALGO_UPSCALE_EROSION` 2 → 3, intention écrite à la ligne. Marche terre−mer p50 **22,00 → 3,00 m**,
> part du trait sous 15 % de pente **0,0 % → 71,4 %**. Gardes verts : masque terre/mer, hypsométrie
> terrestre, puits, invariants de lacs et **table terminale du F74 identiques** (aride ≤ 0,9 %).
> Deux conséquences à juger, pas des no-ops : la profondeur des cuvettes sous-marines tombe
> (p50 26,7 → 8,4 m — le clamp inventait 19 m sous chaque cuvette) et le masque de marécages
> monte ×9,2 en aride (`WETLAND_MAX_DEPTH_M = 3 m`).
>
> ⛔ **ET LE CRITÈRE CHANGE DE COLONNE (Finding 79).** Tracée depuis le **masque** (ce qu'un
> consommateur qui dessine des tags voit), la frange livrée vaut **3 623 éperons ≥ 2 cellules**,
> pas 1 831 : l'isoligne interpolée en cachait **49 %**. Le masque **ne bouge pas d'un chiffre**
> quand le plancher passe de 20 m à 1 m. Δ(masque) contre le pré-incision = **+1 738 / +3 610**.
> **Et la longueur d'onde λ = 14 des F76–F78 n'existe pas sur le masque** (λ 9 à ×2,37, sous la
> ligne des 3×) : c'était le traceur. **Question ouverte à l'auteur : Living Landz trace-t-il la
> côte depuis le masque ou depuis l'isoligne ?**
>
> ⛔ **ET LA BORNE N'EST PAS NEUTRE HYDROLOGIQUEMENT (Finding 81).** Toutes les lignes de la
> table terminale du F74 bougent bien au-delà de 5 % en humide : **budget 502,1 → 474,7 (−5,5 %,
> le climat est calculé SUR le terrain)**, terminal 158,4 → 211,9 (**+33,8 %**, ×0,315 → ×0,446),
> somme chaînée 1 052,4 → 522,1 (**−50,4 %**), trou 343,7 → 262,7 (−23,6 %). Et surtout :
> **cuvettes sous-marines 62 → 20** (aride 62 → 17), aire médiane 0,021 → **31,3 km²** — les
> rebords rendus à la terre font cesser d'exister les petites cuvettes, qui étaient des artefacts
> de rebords noyés. Réconciliation F71-A4 **tenue en forme** (spillway sans source nommée = 0
> partout), sur un compte divisé par trois. Invariants de lacs **verts** aux deux lits et aux deux
> variantes ; exorhéique-sans-exutoire 2 → 1.
>
> **⇒ Le F74 doit être re-dérivé SOUS la borne avant que son remède ne soit conçu.** Et le
> résidu de noyage (7 738 cellules) est **côtier à 87,4 %**, pas intérieur : une seconde couture,
> si elle est un jour écrite, serait côtière — et devra répondre au L502 (`diffuse_channels = true`
> est un choix LEM-correct délibéré).
>
> ⛔ **ET LE TROU DU F74 N'EST PAS UNE ERREUR DE COMPTABILITÉ (Finding 84).** La couture gatée
> `C1DrainageConfig::merged_basin_outlet` re-dérive l'exutoire d'une cuvette FUSIONNÉE (le
> déversoir survivant traversait le col INTERNE de la paire, donc il terminait dans le corps
> qu'il drainait : **3 déversoirs, 501,99 m³/s** en humide). Avec la porte ouverte : déversoirs
> dans leur propre corps **3 → 0**, chaîné **522,1 → 156,4 (−70 %)** — et **terminal ×0,446 et
> trou 262,7 INCHANGÉS**, aux deux lits.
>
> La raison : **Σ évaporation = 0,00 m³/s** sur les vingt cuvettes sous-marines en humide
> (`net_evap = max(0, PE − précip) = 0`), et la plus grande union **n'a aucun rebord par lequel
> déborder**. **365,62 m³/s entrent dans une cuvette fermée, dans un climat qui ne peut pas les
> évaporer, par une géométrie qui ne peut pas les drainer.** Les trois familles de remède
> (a)/(b)/(c) portaient toutes sur *où pointe le déversoir* ; le pointer juste ne déplace rien.
>
> **⇒ Ce qui ferme le trou est H-2 vu du côté des lacs : remonter l'union à son NOUVEAU seuil et
> re-inonder.** Elle atteint alors le rebord supérieur et déborde, ou sa surface croît jusqu'à
> l'équilibre évaporatoire — les deux terminent l'eau. **Stop rule tirée au F84** (exorhéique
> sans exutoire 1 → 2 : la classification du lac et son exutoire sont décidés en deux passes et
> la seconde ne peut plus atteindre la première), donc **rien n'est promu** ; la porte reste
> `None`. Et le retracé referme un 2-cycle explicite (136,38 aller, 19,37 retour) parce qu'il
> tourne APRÈS le détecteur de réciprocité : une couture complète itère les deux.
>
> ⛔ **ET LE POINT FIXE LIVRÉ NE CONVERGE PAS (Finding 85).** `C1DrainageConfig::merged_union_relevel`
> (gatée `None`) fait d'une paire réciproque **une seule région** pour la passe suivante, donc le
> bilan du F39 voit l'union. Résultat mesuré à 8192², deux lits :
>
> * **passes du point fixe 16 → 8 (humide) et 16 → 11 (aride).** Seize est la borne de sécurité
>   de la boucle : **le produit livré en sort par la borne, pas par la convergence.** Tous les
>   chiffres sous-marins des F71 à F84 ont été lus dans un état non convergé ;
> * déversoirs dans leur propre corps **3 → 0** (humide) et **2 → 0** (aride), sans rien retracer
>   et sans supprimer un seul déversoir — **le F85 subsume le F84** et laisse
>   exorhéique-sans-exutoire à 1 au lieu de 2 ;
> * **le trou du F74 ne ferme pas** : humide 262,7 → 182,2 (−31 %), **aride 123,8 → 123,8** ;
> * mais **le BILAN par cuvette ferme** : résidu Σ humide **520,7 → 65,2 m³/s (−87 %)**, aride
>   **28,1 → −1,5** — l'aride se ferme **par évaporation**. **Une grande part du « trou » est un
>   artefact de métrique** : `TERMINAL` ne compte ni l'évaporation ni la terminaison endorhéique ;
> * **le prix est nul** : niveaux +**0,29 m au plus** (5 cuvettes sur 17), terrain noyé **1,7 km²**,
>   `footprint above level` 0 → 0, masque océan inchangé, max `Watercourse` **bit-identique**
>   (19,243 m³/s, 3,37 cases). Le plus gros `Spillway` **baisse** : 360,34 → 92,97 m³/s, **14,58 →
>   7,40 cases** — à porte fermée la plus large ligne bleue du produit est un déversoir à **trois
>   fois** la cible 4–5 cases du F69, et personne ne l'a regardée ;
> * **le F29 (H2) est vivant dans le produit** : à porte fermée le résidu est porté par des
>   cuvettes marquées `exorheic = true` avec `out = 0,00` — dont une de **0,007 km²** avec une aire
>   d'équilibre de **917,5 km²**. Le F39 disait cette classe fermée « par construction ».
>
> **⇒ Rien n'est promu (stop rule : exorhéique-sans-exutoire = 1, ROUGE sur l'état livré, aux deux
> positions de porte).** Reste 65,2 m³/s en humide dans **deux** cuvettes dont `a_eq` égale
> `a_spill` à la troisième décimale, qui ne peuvent ni déborder (aucun exutoire ne se trace) ni
> évaporer (il pleut plus qu'il n'évapore). **Aucune couture comptable ne fermera cela.**
>
> ✅ **PROMU, ET LE TROU EST RETIRÉ COMME MÉTRIQUE DE CONSERVATION (Finding 86).** Deux
> changements de production, aucun gaté :
>
> 1. **La ligne H2** — `LakeType::Unresolved` + `resolve_exorheic_without_outlet`, en fin de
>    `assemble_hd_drainage`. Le F29 avait nommé le défaut, le F30 avait **interdit de le
>    relabelliser endorhéique** (« ce serait affirmer un bassin qui gagne plus qu'il ne perd »),
>    le F39 l'avait fermé **sur la moitié sous-marine seulement** : `apply_lake_water_balance`
>    (`drainage.rs:1078`, H-1c, écrite après) pose encore `Exorheic` depuis `a_eq ≥ a_sill`
>    **sans tracé**. Un seul lac viole : **id 55, un lac de SURFACE** de 19,2 km² à 51,05 m.
>    exorhéique-sans-exutoire **1 → 0** (humide), 0 → 0 (aride). ⚠️ **Item de contrat** : la
>    variante est sérialisée dans `lakes.json` et a cassé **quatre** `match` exhaustifs, dont
>    trois dans le viz — qui affichent désormais « ⚠ sans exutoire (F86) » en rouge.
> 2. **Le relevel du F85, promu** — `ALGO_DRAINAGE` 6 → 7, `ALGO_HD_DRAINAGE` 8 → 9. Le point
>    fixe converge (**8** et **11** passes contre la borne de 16) et épuiser la borne est
>    désormais une assertion. **La non-convergence livrée est attribuée** : une **alternance de
>    période 2** entre deux classes, d'amplitude constante **23,472 m³/s** (humide) et **4,149**
>    (aride) — une paire réciproque qui se nourrit mutuellement à chaque passe, exactement ce que
>    le L2911 avait nommé (« l'hypothèse DAG du F40 était ASSERTÉE, jamais vérifiée ») en ne
>    corrigeant que la SORTIE, jamais la boucle.
>
> **⇒ `TERMINAL` / `TROU` ne sont plus des chiffres de conservation** : ils mesurent « la part du
> budget qui atteint l'océan **par un tronçon terminal** ». L'instrument de conservation est le
> **budget à cinq termes**, et il ferme : **102,5 % en aride**, **92,5 % en humide** contre
> **70,6 %** pour le contrôle A/B. Le terme qui le ferme n'était dans aucune de nos deux listes :
> **le ruissellement qui ARRIVE physiquement à la côte — 180,5 m³/s (38 %) humide et 90,6 (66 %)
> aride — contre 81,0 (17 %) et 13,1 (9,5 %) crédités par la somme des tronçons terminaux.**
> Les F71 et F74 mesuraient la résolution du réseau et l'appelaient un trou.
>
> **Reste UNE cuvette** : 45,04 m³/s (apport LOCAL ; les 65,2 du F85 étaient un chiffre chaîné)
> dans la cuvette humide 1000001, fond à −10,60 m, pleine à ras d'un seuil à 76,25 m dont la
> première cellule aval est **douze pas u16 plus HAUT**. Elle ne peut ni déborder ni évaporer.
> **Nommé, non corrigé** : le côté SORTIE est encore chaîné (le système sous-marin sur-ferme à
> 131,5 % humide), et `a_eq_km2` rapporte `a_spill` quand `a_eq` est infini — l'ambiguïté qui a
> trompé le F85.
>
> ✅ **LE BUDGET FERME À ±3,3 % DANS LES QUATRE CAS (Finding 87).** Un changement de production
> survit, un a été écrit puis **retiré**.
>
> 1. **F, livré** — `BasinSummary::local_outflow_m3s` = `max(0, apport local − évaporation)`, et
>    `a_eq_km2` rapporte enfin **l'infini comme un infini**, avec `a_eq_is_infinite` à côté :
>    l'ambiguïté qui avait trompé le F85 est sortie du type. Le système sous-marin, local **des
>    deux côtés**, ferme à **100,0 %** (humide, les deux positions de porte). Et le **budget à six
>    termes** (arrivée côtière + `Spillway→océan` **local** + évaporation sous la mer + évaporation
>    de surface + retenu + lacs `Unresolved`) ferme à **103,3 / 98,4 / 103,2 / 103,1 %** — deux
>    lits climatiques, deux positions de porte. Les **7,5 %** résiduels du F86 n'étaient pas un
>    terme manquant : c'était `Spillway → océan` qui sur-créditait de **80,5 m³/s** parce qu'il
>    était chaîné (la même table en chaîné lit **120,3 %**).
> 2. **D, retiré** — tracer avant d'étiqueter sur le chemin des lacs *détectés* a échoué à la garde
>    du tour même (« aucun autre lac ne change ») : **19** lacs changent d'état, pas 1, et **les
>    dix-neuf** pour la même raison — *« la descente revient dans sa propre empreinte »*. La
>    surface d'un lac détecté est **plate**, donc une marche D8 partie de l'intérieur n'a aucun
>    gradient pour en sortir. Ce chemin a besoin du **col/échappée du F37c**, qui est un mécanisme.
>    `surface_lake_outlet_trace` reste comme l'instrument qui l'a prouvé.
>
> **Trois nombres qui remplacent trois impressions.**
>
> - **La plus large ligne bleue est cohérente** : les 92,97 m³/s vont **à la MER**, **115,33 m de
>   chute sur 61 cellules (≈ 3,0 km)**, pente **3,36·10⁻²** — au-dessus d'un seuil rocheux (10⁻²),
>   en dessous d'une cascade (10⁻¹). **Aucun déversoir livré n'a une chute ≤ 0** (0 sur 11 humide,
>   0 sur 3 aride) ; le contrôle A/B en donne **4** (et non 3 : le quatrième est un lien entre deux
>   corps *distincts* au même niveau — un détroit est plus large que « un corps qui se déverse en
>   lui-même »). Deux déversoirs livrés **sont** de vraies chutes d'eau : 34,29 m³/s sur 337,8 m
>   (2,10·10⁻¹) et 17,60 m³/s sur 418,4 m (1,44·10⁻¹).
> - **Le lac de 623,6 m est INCISÉ**, aux deux lits : coupe médiane de **165,80 m** (humide) et
>   **161,30 m** (aride) sur l'empreinte, fond livré à **−8,67 m** contre **+12,63 m** avant
>   incision. ⚠️ **Correction à la description du F83 tenue depuis** : le plancher borne la CIBLE de
>   **relaxation**, pas la **diffusion de versant** ni le talus — donc une cellule peut encore être
>   portée sous le niveau de la mer après le F83, et c'est ce que le F80-B1b mesurait à **88,7 %**.
> - **L'arrivée côtière diffuse est de la physique, sa queue ne l'est pas** : **87,7 %** des
>   cellules d'arrivée non jaugées sont sous `head_threshold` — du ruissellement de versant —, et
>   ce chiffre est **identique dans les deux lits** : la question était climatique, la réponse est
>   géométrique. Les **12,3 %** restants, ~**4 850 cellules**, portent des bassins jusqu'à
>   **19,6 km²** sans rivière tracée.
>
> **Nommé, non corrigé** : l'**évaporation est la dernière grandeur chaînée** du bilan (le
> sous-marin aride ferme à **112,1 %** pour cette seule raison) ; la cuvette 1000001 est bloquée
> derrière une marche de **0,014 m** dont l'autre côté n'a **aucune direction D8** (un item de
> routage sur plat) ; et le col que le F86 avait opposé au F36 était **mon propre instrument** — le
> plus bas voisin de l'empreinte est au niveau **par construction** (F36 PART B, L1250).
>
> ⛔ **LE BILAN DU F87 EST UNE PARTITION, PAS UNE CONSERVATION (Finding 88).** Deux changements de
> production, aucun gaté, et **l'un des deux était faux à la première écriture — c'est la mesure qui
> l'a trouvé.**
>
> 1. **A, livré** — `BasinSummary::local_evaporation_m3s = min(évaporation, apport local)`. Le
>    sous-marin aride passe de 112,1 % à **100,0 %** ; l'excès chaîné mesuré vaut **5,70 m³/s**.
>    ⚠️ Et la docstring le dit **en production** : `local_out + local_evap == local_in` est une
>    **identité algébrique**, donc ce 100,0 % *se déduit*. Le « 100,0 % » humide du F87-F était déjà
>    cette identité, énoncée comme une mesure (règle 10).
> 2. **C, livré** — le **col/échappée du F37c porté aux lacs de surface** : `surface_lake_escape`,
>    `surface_lake_escape_trace`, un `walk_to_sink` partagé, et `apply_lake_water_balance` qui
>    **trace avant d'étiqueter**. **57 des 59** lacs se résolvent ; les 2 restants nomment
>    `ReturnsIntoOwnFootprint`. `UnresolvedReason` (enum à sept variantes) voyage dans `C1Lake`,
>    dans `lakes.json` et jusqu'aux trois rendus du viz. `ALGO_DRAINAGE` 7 → 8, `ALGO_HD_DRAINAGE`
>    9 → 10. ⚠️ **Première écriture : 35 lacs sur 59 en échec**, tous sur
>    `SaddleHasNoLowerNeighbour`, parce que je prenais *la plus basse cellule du rebord* au lieu de
>    *la plus basse cellule du rebord depuis laquelle l'eau peut sortir*. Sur un champ breaché la
>    première est souvent un puits dans un chenal de percée. Le F37c le disait et l'instrument du
>    F87-E l'avait déjà implémenté correctement.
>
> **Le résultat qui déclasse le bilan du F87.** `arrivée côtière + Σ apport LOCAL sous-marin =
> **474,7** contre un budget de **474,7**, **Δ −0,00, aux deux lits**. Les deux premiers termes
> épuisent le budget **exactement** — c'est ce qu'une accumulation de ruissellement *est*. Donc ce
> tableau **ne mesure pas la conservation le long du routage** : il redit la partition, et tout
> chiffre au-dessus de 100 % y est un terme compté deux fois. Le sixième terme en était un :
> **tous** les lacs `Unresolved` routent, en traversant leur propre surface plate, **jusqu'à la
> mer** — leur apport était déjà dans les 180,5. **Le point de claim du F87 est réécrit sur place
> en « partition vérifiée ».**
>
> **Et lu à la FIN de chaque chaîne** (ce que le premier saut cachait : **98,1 m³/s**), le routage
> ferme à **96,3 %** (humide) et **98,4 %** (aride), et le reste n'est plus un pourcentage mais un
> objet : **un déversoir par lit qui termine « nulle part d'identifiable »** — 19,5 m³/s humide,
> 4,1 aride, `SpillwayTermination::to_nothing == 1`, la famille du F38.
>
> **Trois corrections au dossier, toutes à mes propres tours.** Le F87-F ci-dessus ; le **F87-E
> retiré** (l'échappée de 1000001 porte une direction D8 sur **les deux** champs et descend à la mer
> en **44 pas** — le défaut est dans le chemin qui **émet** le déversoir, pas dans le terrain) ; et
> le « 165,8 m » du F87-B qui est une **différence de médianes**, la statistique appariée valant
> **158,31 m**.
>
> **Règle 12, écrite en code ce tour** : `flow_field_hash` + `declared_flow` — tout banc qui trace
> déclare son champ et le prouve par empreinte. **Elle a tiré à sa première sortie** : `dr.flow` et
> le `compute_flow` du banc du F87 diffèrent sur **512 130 directions (0,763 %)**. ⚠️ Et le chiffre
> du F87 y survit : l'arrivée côtière vaut **180,5 sur les deux** (Δ 0,0). Le garde a raison sur le
> principe, la conclusion était intacte. **Elle ne couvre ni les unités, ni les populations, ni les
> cartes** — et ce sont eux qui ont coûté trois des quatre re-runs (une pente en unités normalisées,
> 11 300× trop petite ; un p50 d'érodabilité pris sur 55,8 M de cellules de mer ; un classifieur de
> chaîne aveugle aux lacs détectés).
>
> **Nommé, non corrigé** : `FlatPerturbation` **n'a aucun antécédent dans ce dossier** et décide
> 7,1 % des cellules de terre (L733) ; les 45,0 m³/s de 1000001 ont une échappée tracée mais aucun
> déversoir émis ; `stream_km2 = 20` est la **seule** cause (100,0 %) des 4 843 cellules côtières à
> bassin jusqu'à 19,65 km² sans rivière exportée — décision consommateur, non touchée.
>
> ⛔ **ET LE LAC DE 623,6 M N'EST PAS UN OBJET, C'EST UNE CLASSE (Finding 88-D4).** **16 des 53**
> corps livrés ≥ 1 km² ont une coupe médiane > 50 m **et** une paroi p50 > 30° — et le lac promu
> n'est pas le cas extrême : **le corps 1 a perdu 862,2 m avec une paroi à 74,9°** sur 33,4 km². La
> roche y est ordinaire (**0,92×** l'érodabilité médiane des terres) et le volcan le plus proche est
> à **61,6 km** : ni C-2 ni C-3. C'est `E = K·A^m·S^n` à Courant 1353 dans des cuvettes fermées,
> atterrissant en **une** passe (la passe 1 porte **93,7 %** de la coupe du corps, 55,5 % de celle
> du fond). **Si l'œil dit que seize est trop, l'item est l'échelle de temps du F44** — spécifiée,
> jamais implémentée.
>
> ✅ **CRITÈRE ATTEINT (Finding 80) — pour la première fois du chantier.**
> `StreamPowerConfig::base_level_floor`, **gatée `None`, octet-identique**, borne la CIBLE de
> relaxation : `h_r_eff = max(h_r, min(h_o, sea + ε))`. À ε = 0,5 m (3,4 pas u16) :
> **Δ(masque, post-u16) = +0 à l'échelle cellulaire sur les quatre grilles** (f32 8192, u16 8192,
> terrain 2048, océan 1024), Δ(≥ 1 km) = −2 / −1 / −1 / +0, contour 1 634 → 1 639 km (+0,3 %).
> Le F6 avait nommé le défaut (*« No incision bound — floors planed to base level »*) et choisi un
> cadran de DURÉE ; `base_level` n'avait jamais été écrit.
>
> **Le prix, mesuré au bon étage** (pré-clamp) : incision refusée p50 **1,08 m**, p90 5,58 m, sur
> **14 554** cellules, en biefs de **1 à 2 cellules** — pas de banquette plane à +ε. Intérieur :
> hypsométrie **appariée +0,30 m** (la −17,26 m non appariée est l'effet de population des
> 290 903 cellules rendues à la terre, à 2,32 m de moyenne), dépressions +0,41 %, part chenal
> +0,36 pt. Shader : composantes terrestres de 4–8 cellules **504 → 63** (pré-incision 50),
> îlots d'un pixel 5 → 0, bande de sable 16 426 → 9 900 px (pré-incision 9 874).
>
> **Résidu attribué** : 7 738 cellules franchissent encore zéro (−97,4 %), portées à **88,7 % par
> la diffusion de versant** (`diffuse_channels: true`, elle ne saute que la mer) et à 3,2 % par le
> talus ; 11,3 % non attribués. **Non mesuré, non revendiqué** : la table terminale du F74 sous la
> borne.
>
> ⛔ **LA CLASSE « DÉPÔT » EST ÉCARTÉE (Finding 79).** Sur le champ pré-clamp, les cellules noyées
> ont une médiane de **0,88 m** mais un p90 de **18,7 m** : moitié centimétrique, moitié lits
> sous-marins réels. Combler les plus faibles **aggrave** la mesure (Δ masque +3 627 → **+4 800**
> à 5 cm : un creux comblé devient un isthme d'une cellule), et à 5 m — trois quarts des cellules
> relevées — **43 % de la frange subsiste**. Seul « tout relever » atteint Δ = 0.
>
> **Le levier le plus fort mesuré à ce jour n'est pas côtier** : à la valeur ANCRÉE de la
> diffusion (1,28 — le défaut du F59 §4, chiffré là-bas à 4 % de la divergence de travail), la
> fourrure à l'échelle cellulaire passe de **1 831 à 189 (−90 %)** et la côte de 6 206 à 3 538 km.

> ✅ **PROMU EN PRODUCTION (Finding 83), 2026-09-15.** `relief_v3` embarque désormais
> `base_level_floor: Some(BaseLevelFloor { epsilon_m: 0,5 })`, `free_above_km2: None` ;
> `ALGO_UPSCALE_EROSION` 3 → 4. Le mot du critère, fixé par l'auteur après lecture des deux
> mondes dans Living Landz toutes closures actives : **aucun éperon AJOUTÉ**, le recul diffusif
> (6 764 cellules centimétriques, +5 km de contour) est accepté comme une côte qui vit.
>
> **Les seize gardes rejouées sur le CHEMIN DE PRODUCTION reproduisent le banc du Finding 81
> au chiffre**, et le champ au BIT : FNV pré-incision `0xc359c555ec22175c`, borne désactivée
> `0x4719ff260186645a` (le monde d'avant, octet pour octet), production `0x1a022c8785c47eed`.
> Table terminale F74 humide 474,7 / 82:81,0 / 5:130,9 / 8:522,1 / 211,9 (×0,446) / trou 262,7 ;
> aride 137,7 / ×0,101 / 123,8 ; invariants de lacs verts ; `BasinSummary` 20 ; fusions
> `[121234, 6245] [1,1] [2,2]`. **Le toggle viz devient le contrôle A/B** (coché = production ;
> décocher = `base_level_off`, le monde d'avant).


## ⛔ CRITÈRE « CÔTE ORGANIQUE » — écrit avant mesure (ADR Finding 83)

La cible du chantier suivant n'est pas « indentée », elle est **organique** : une côte organique
n'a **ni longueur d'onde ni générateur visible**. Sur le masque CONSOMMATEUR (u16 post-export,
`marching_squares` à 0,5), instrument = **courbe de Richardson au compas** (distance de corde,
8 origines moyennées, 13 règles log-espacées de 100 m à 10 km, polygones ≥ 20 km) :

1. **linéarité** — `log L` contre `log(règle)` droite sur les deux décades, résidu déclaré ;
2. **dimension** — la pente donne **D ∈ [1,10 ; 1,30]** (Grande-Bretagne ≈ 1,25 ; côte lisse 1,00) ;
3. **pas de générateur** — le spectre des positions d'excursions sans pic au-dessus de **3× le blanc** ;
4. **rien de fabriqué** — **Δ(éperons ≥ 2 cellules) contre le pré-incision reste 0** hors estuaire,
   et toute excursion ajoutée porte à sa tête un chenal d'aire ≥ A_est ;
5. **et l'œil de l'auteur**, dans Living Landz, avec les comptes à côté.

**L'instrument est étalonné, pas supposé** : cercle rastérisé D vrai 1,000 → mesuré **1,002** ;
flocon de Koch g7 D vrai 1,2619 → mesuré **1,244**. Le plancher d'instrument est nul (le
`marching_squares` interpole, le contour n'est pas un escalier). ⚠️ Sur une courbe quasi plate le
R² est petit pour des raisons étrangères à la linéarité : **la linéarité se juge au résidu %**.

| étage (u16 8192²) | D | **étalonné** | résidu | D < 1 km | D > 1 km |
|---|---|---|---|---|---|
| pré-incision | 1,051 | **1,053** | 8,16 % | 1,013 | 1,115 |
| **livré (pré-83, la fourrure)** | **1,273** | **1,293** | **4,74 %** | **1,292** | 1,229 |
| **production (bornée)** | 1,051 | **1,053** | 8,07 % | **1,013** | **1,116** |

> ⚠️ **ET LE RÉSULTAT EST INCONFORTABLE, il est écrit tel quel.** La fourrure passe les pattes
> 1, 2 et 3 : loi de puissance **propre** sur deux décades (le plus petit résidu des trois), D
> **1,29** — dans la fenêtre écrite pour une vraie côte — et spectre à 2,43× seulement. Elle
> n'échoue qu'à la patte 4 (Δ ≥ 2 cellules = +3 597). **Tout l'appareil fractal ne distingue pas
> la frange d'une côte plausible ; seul le Δ contre le pré-incision le fait**, et c'est le
> critère de l'auteur, pas celui du géomètre.
> Et la côte que la borne nous donne **échoue à la patte 2** : D = 1,05, et sous 1 km **1,013**.
> Le pré-incision aussi : **l'autorité que le chantier utilise depuis le F76 n'a jamais été une
> côte organique.**

**⇒ B3 : le manque est SUB-KILOMÉTRIQUE, et c'est tout le manque.** Au-dessus d'un kilomètre la
côte bornée a déjà D = 1,116 (1,137 étalonné), dans la fenêtre ; en dessous, 1,013 — plate. Le
premier suspect est nommé et il est délibéré : **`coastal_amplitude_band = 0.30`,
`upscale.rs:593`**, qui amortit le FBM à ~0 **au trait d'eau par conception** (#151). C'est là
que C-4 commence.

**⛔ LE LEVIER ESTUAIRE EST MESURÉ ET IL ÉCHOUE (Finding 83-B2).** `free_above_km2` (PROXY,
hors production) libère la borne au-dessus de A_est. Balayage 50 / 200 / 1 000 km² :
D bouge de **+0,003 au plus** ; à 50 km² **70,5 % des 44 excursions ajoutées n'ont aucun chenal
à leur tête** (c'est la fourrure sous un autre nom, par la diffusion de versant du F80-B1b) ; à
200 km² la tête est propre à **100 %** mais l'effet est **3 entailles de 0,2 km de long, 1
cellule de large, 24 cm de fond** ; à 1 000 km² **rien ne bouge du tout**. Contre la Rance
(≈ 20 km × 0,5–2 km) c'est **10 à 40× trop étroit et 14× trop court**. **Une ria est une vallée
noyée par une mer qui est montée ; la mer d'Ymir n'a jamais bougé, et aucune fonction de l'aire
drainée ne remplace une histoire.**

**Ni l'état livré ni la loi n'approchent d'un ordre de grandeur les quatre premières lignes.**
Un remède qui passe la côte en changeant le continent n'a rien passé.

**Ce qui est ÉLIMINÉ** (Finding 76) : le plat au niveau de la mer (r = 0,051, 98,5 % des éperons
sur un plat de longueur nulle), la classification du zéro (0 cellule exactement à `sea_level`),
et le FBM/les closures (×229 entre « sans incision » et « livré » à l'échelle cellulaire).
**Ce qui reste** : une texture parallèle (R 0,879 contre 0,000), sub-métrique (encoche médiane
0,84 m), non chenalisée (86,3 % des pointes sous `A_c`), **à longueur d'onde fixe de 14 cellules
(684 m) à ×3,75 le blanc**. C'est une instabilité de rilling (Finding 10), pas un réseau de
drainage — et sa longueur d'onde est fixée par la compétition incision/diffusion, pas par `A_c`.

## ⛔ L'HYPSOMÉTRIE EST LE VERROU — priorité 1, devant tout le reste

**Quatre chantiers indépendants sont venus buter dessus**, chacun par un chemin différent, et
c'est la raison de l'inscrire ici plutôt que de la redécouvrir une cinquième fois :

1. **la navigabilité** (ADR Finding 47) — le débit n'est pas stable en résolution (max 10,3 m³/s
   à 2048² contre 1,04 à 8192², la classe « petite embarcation » ×50), donc les seuils en m³/s
   ne peuvent pas être calibrés ;
2. **le mécanisme 1 du Finding 43** — l'incision laisse 90 % des terres à grille fine à une
   branche de versant inerte ;
3. **le cadran temporel de H-2** (Finding 44) — il lui faut une échelle de temps explicite, et
   celle-ci n'a de sens qu'une fois le terme de versant refait ;
4. **le seuil de tête de chenal** (Finding 56) — `S_ref` vaut 0,3319 à 2048² contre 0,1128 à
   8192², donc **aucune calibration de tête de chenal ne peut tenir aux deux grilles** tant que
   l'hypsométrie ne converge pas. La loi comble 56 % / 71 % de l'écart de franges côtières ; le
   reste est hérité de ce défaut, pas de la loi ;
5. **la FORME DES VALLÉES** (Finding 56b) — le W/D de vallée mesuré sur un transect
   perpendiculaire à l'écoulement donne 35–98 à 2048² contre 17–42 à 8192², et **la loi le
   resserre à la grille grossière tout en l'élargissant à la grille fine**. Les deux grilles ne
   sont pas dans le même régime morphologique, et aucune calibration de forme de vallée ne peut
   être ancrée tant que c'est le cas. Cinquième chemin indépendant vers le même verrou.

Ce n'est donc plus un chantier parmi d'autres : **c'est le nœud qui bloque la convergence de
tout ce qui dépend de la pente.**

> ⛔ **BLOQUANT, Finding 61 : la question de convergence en résolution est MAL POSÉE à compte
> d'itérations fixe.** `iterations` est un cadran de DURÉE (`E = K·A^m·S^n` n'a pas de terme de
> soulèvement, donc pas d'état stationnaire : l'attracteur est le niveau de base). Mesuré : à
> 2048² le travail croît 455,7 → 633,0 → 757,5 → 835,9 m avec un Δ décroissant, et l'altitude
> moyenne tombe 447 → 107 m ; à **8192² le Δ décroît puis REMONTE** (+31,6, +24,9, **+38,1**),
> donc aucune limite n'y est extrapolable. Et le rapport inter-grilles vaut 2,72 / 3,18 / 3,38 /
> 3,19 selon le budget — **24 % d'amplitude**.
>
> Conséquence : comparer les deux grilles à compute égal n'est pas les comparer à érosion égale.
> Il faut les apparier sur un **état** (travail égal, ou moyenne hypsométrique égale), ce qui
> exige l'échelle de temps explicite que le Finding 44 a spécifiée sans l'implémenter. Et toute
> reformulation candidate du stream power **qui suppose un état stationnaire est hors domaine
> ici**.
>
> ⛔ **ET L'APPARIEMENT À ÉROSION ÉGALE NE SUFFIT PAS NON PLUS (Finding 62).** Mesuré à travail
> égal (201,5 contre 197,6 m, `iterations` = 2 des deux côtés, seul `dt` bougé) : la **masse**
> s'apparie — terre 16,325 % contre 16,443 %, moyenne 685,7 contre 686,5 m — et la **forme** non :
> part chenal **80,75 % contre 8,40 %**, un facteur 9,6, avec des quantiles d'altitude de formes
> différentes (p10 51 contre 29, p90 1507 contre 1607).
>
> **Et 8192² a un PLAFOND de travail de 311,2 m** : `dt` sur 40 000× d'amplitude ne dépasse pas
> cette valeur, alors que l'état livré à 2048² vaut 633,0 m. **Aucune durée n'apparie les deux
> grilles** — elles n'ont pas le même ensemble atteignable, et le plafond est le compte de
> cellules d'`A_c` vu de l'autre bout du cadran. Donc **une échelle de temps ne peut pas être le
> remède** : elle ne peut que placer une grille dans son propre ensemble atteignable.
>
> Corollaire à porter dans toute comparaison future : à 8192² la part chenal tombe de 53,62 % en
> entrée à 8,40 % au point livré et 5,71 % à huit itérations — **le réseau s'éteint**. Un rapport
> entre un réseau vivant et un réseau qui s'éteint n'est pas un écart de magnitude sur un même
> objet.

**Un corollaire mesuré (Finding 56b), qui bloque aussi la largeur des chenaux :** dans le banc
d'essai aride, **aucun des 51 tronçons d'ordre 5 à 8192² ne porte de débit** et le chenal le plus
large de toute la grille fine fait 5,2 m. Le tronc de la hiérarchie est à sec pendant que l'eau
reste dans les tronçons côtiers d'ordre bas. C'est cohérent pour un continent désertique, mais
cela veut dire que **`w = 5·Q^0,5` n'a rien à calibrer ici** — même verrou, par le débit.

### Ce qui est DÉJÀ SU — à ne pas re-dériver

- **Attribution 63 / 18 / 18** (Finding 43) : 63 % de l'excès d'altitude à 8192² vient de la
  **partition de régime fluvial/versant**, ~18 % de la dispersion MFD, et **~50 % du résidu
  (73 m) reste NON ATTRIBUÉ** — il survit à `A_c = 0` **et** à MFD éteint.
- ⚠️ **CORRIGÉ par le bloc A — deux étiquettes périmées dans les lignes qui suivent.**
  1. **`S_ref` n'est pas une trace de la cause.** 0,3319 / 0,1128 est mesuré sur le champ
     **ÉRODÉ**, aux têtes de chenal. Le champ **avant incision** a la même distribution de pente
     D8 aux deux grilles à 3–4 % près sur les queues. L'écart de `S_ref` est donc une
     **conséquence** de la divergence, et la circularité de la loi `A_c(S)` est **démontrée** :
     elle est calibrée sur une grandeur produite par le défaut qu'elle contourne.
  2. **« Branche inerte » est trop fort : l'opérateur est CONSERVATIF.** 66,5 % des terres
     communes déplacées de plus d'1 m à 8192² pour ~3 m de moyenne. Il fait beaucoup, tout
     à masse nulle.
  3. **La rugosité de l'upscale ×128 n'est PAS la cause racine** : le champ avant incision est
     invariant en résolution en moyenne (865,4 / 865,5 m), en quantiles d'altitude et en
     fraction émergée (16,88 % / 16,88 %). Ajouter une octave ne change rien (< 0,5 %) et
     `amplitude_base ×4` ne change **rien du tout** — le budget de relief C-1 borne la
     rugosité à une valeur invariante en résolution, par construction.

- **Le mécanisme 1 est « un seuil correct alimentant une branche inerte »** : `A_c = 0,1 km²` est
  juste en km², mais la fraction de terres qu'il couvre passe de 55,4 % à 2048² à 10,1 % à 8192²
  (loi aire-fréquence sous-linéaire), et les 90 % restants sont confiés à un terme de versant qui
  ne fait rien (`diffusion = 0` change le résultat de 3 m).
- **Le remède n'est PAS un laplacien plus fort.** Mesuré et réfuté (Finding 45b) : le laplacien
  est **conservatif** — il remplit les vallées autant qu'il abaisse les crêtes — et le renforcer
  a fait passer les bassins sous-marins de 43 à 1994 à 8192² tout en n'élevant l'altitude
  moyenne que de 17 m. Il faut un terme **transport-limited à flux sédimentaire explicite**, le
  seul type d'agent qui **retire** de la masse des versants et la **livre** au réseau.
- **Ce terme exige l'échelle de temps explicite** que le Finding 44 a **spécifiée sans
  l'implémenter** : une diffusivité en m²/an n'a rien à multiplier tant que `dt = 1.0` est un
  simple placeholder d'unité. Et le Finding 44 a montré au passage que `K` et la durée ne sont
  **pas observables séparément** — l'observable est `k_time = K·dt·iterations`.

### Ordre imposé par ces dépendances — ⛔ RÉVISÉ AU FINDING 89

`échelle de temps explicite (F44)` → `terme de versant transport-limited (F43 méc. 1)` →
`hypsométrie convergente` → puis, redevenus calibrables : les seuils de navigabilité en m³/s,
le cadran de H-2, et la calibration de tête de chenal aux deux résolutions.

> ⛔ **LE PREMIER MAILLON EST DÉJÀ LIVRÉ, ET LE VRAI BLOQUEUR EST AILLEURS (Finding 89).** Aucun
> changement de production ce tour ; trois lectures et deux mesures.
>
> 1. **Le F44 n'est pas « spécifié sans être implémenté ».** Ses items 1 (`dt` en années,
>    `T = iterations·dt`) et 4 (la borne CFL, `dt_max = cell_m / c`) **ont été livrés** —
>    `SHIPPED_K_TIME`, `k_time()`, `dt_max_yr()`, `cfl_iterations()`, `courant()`,
>    `timescale_plan()` dans `erosion/stream_power.rs`, avec deux tests d'épinglage. Seuls les
>    items 2 et 3 (dériver le nombre de passes d'une durée et d'une borne) manquent. **Le chantier
>    citait une phrase de la note de reprise, pas le finding.**
> 2. **Le bloqueur est le terme de SOULÈVEMENT absent, et le F61 l'avait prouvé** : sans lui, le
>    seul point fixe de `h ← (h + f·h_r)/(1+f)` est le niveau de base, donc *« "convergence" ne peut
>    signifier que planation »*. **Tout critère qui demande un relief conservé demande à l'équation
>    ce qu'elle ne peut pas produire.** Nommé, non décidé.
> 3. **L'objection de coût du F44 est RÉFUTÉE.** Une passe d'incision coûte **22,1 s** à 8192²
>    (l'incision est **94 %** du pipeline : 44,2 s sur 47,0), le coût marginal est **constant à
>    4,6 %** sur 2→3→4, et les **7 399 passes** que le F44 dérivait font **1,9 jour**. Mille passes :
>    **6,3 h**. Son objection d'identifiabilité tient, elle : la littérature ne borne `T` qu'entre
>    **4,5·10⁵ et 3,6·10⁸ ans**.
> 4. **La boucle temporelle n'a qu'UN étage.** Le F60 (*« l'érosion ne voit jamais le climat »*) rend
>    le climat, le bilan des lacs, les rivières et la bathymétrie tous en queue, une fois :
>    `N × 22,1 s + 12 s`. La question « tous les N pas ? » n'a pas lieu.
>
> **Le temps zéro, chiffré (lit humide, livré).** **56** plans d'eau, **6 487 km² = 24,05 %** des
> 26 969 km² de terre, **586,3 km³** — le 23 % de l'auteur confirmé, soit **24×** le critère
> (< 1 %) et **80×** la France. ⚠️ **Mais le même champ donne 5,28 % en aride** : le terrain est
> climato-indépendant, donc **les quatre cinquièmes du 24 % sont le bilan hydrique humide** (F39 :
> `net_evap = 0` ⇒ `a_eq = ∞` ⇒ remplissage au seuil **par construction**), pas le relief. Et
> **70,4 % de l'aire et 86,8 % du volume** sont dans des corps de **plus de 100 m** de fond ; les
> corps sous 30 m ne portent que **2,9 %** de l'aire. Les **16** corps de la classe canyon (F88-D4,
> reproduit) portent **30,0 % de l'aire et 65,4 % du volume**.
>
> ⛔ **ET LE +115 EST `breach_monotone`, PAS L'INCISION (Finding 90).** Bissection par étage,
> aucun changement de production.
>
> | étage, Δ(≥ 1 km) vs pré-incision | 8192² | 4096² | 2048² | 1024² |
> |---|---|---|---|---|
> | **ÉRODÉ** (sortie d'incision) | **−1** | 0 | −1 | 0 |
> | **BRÈCHÉ** (`hd_assembly`) | **+110** | +176 | 0 | +1 |
> | la ligne du F89-C5 (étages mixtes) | **+115** | +184 | 0 | +1 |
>
> **L'incision n'ajoute aucun éperon** — la jambe 4 est tenue à l'étage sur lequel elle a été
> écrite. Les 115 se décomposent en **+110** (la brèche) et **+5** (l'autorité lue **non brèchée**
> face à un livré brèché : le pré-incision brèché a **25** éperons, pas 20). **117 ajoutés, dont
> 113 — 96,6 % — ont à moins de 4 cellules une cellule TERRE dans l'érodé et MER dans le brèché.**
> Le « −2 » du F80 n'était jamais un Δ contre pré-incision : c'était un Δ de Δ entre deux variantes
> **livrées**.
>
> **Et ce qui a déplacé le chiffre est le plancher du F83, pas le relevel** : plancher coupé
> **1 816** éperons (le +1 828 du F80 reproduit à **+1 791**), plancher mis **135** — soit **×13,5**.
> Relevel coupé : **135, identique à l'unité**. La chaîne relevel → brèche est **réfutée**.
>
> **`breach_monotone` n'a aucun plancher, par construction** : `is_base` arrête la marche **à** une
> cellule déjà sous la mer et ne borne **jamais** la cible, donc tout chemin d'exutoire d'un puits
> sous-marin y descend. Elle **noie 10 375 cellules** (24,74 km², p50 **0,44 m**, max 11,91 m — le
> −8,67 du F89-D est au **99,4ᵉ centile**) et n'en relève aucune. ⚠️ **Mais elle ne crée que 1,6 %
> des cuvettes sous-marines encloses** (6 100 sur 387 710) : **aucune ligne n'est due au point de
> claim des F71–F88**, la machinerie des chaînes travaillait sur de vraies dépressions.
>
> ⛔ **ET LE CHAMP LIVRÉ EST SUR-ÉRODÉ PAR SA PROPRE NUMÉRIQUE (Finding 90-C).** À `k_time` **tenu
> à 9 000** — même budget d'érosion, seule la distribution change — le F44 item 2 mesuré :
>
> | jalon | Courant | côte érodée | côte brèchée | noyées | relief p50 | lacs | classe canyon |
> |---|---|---|---|---|---|---|---|
> | **2** (livré) | 3 699 | 19 (Δ −1) | 135 (**Δ +115**) | 10 375 | 424,2 m | 24,05 % | **16 / 53** |
> | **10** | 740 | 18 (Δ −2) | 106 (**Δ +86**) | 15 549 | **503,4 m (+18,7 %)** | 18,99 % | **5 / 30** |
> | **100** | 74 | 17 (Δ −3) | **47 (Δ +27)** | **2 880** | **505,2 m (+19,1 %)** | **17,87 %** | **3 / 22** |
>
> **Le relief MONTE de 19 %**, et nos deux prédictions le voyaient descendre de 15 à 35 %. À budget
> tenu, intégrer fidèlement **retire moins de masse** : c'est le « stable mais pas intégrant » du
> F44 lu par l'autre bout. **Donc la classe canyon est NUMÉRIQUE** (16 → 3, −81 %, à budget
> constant) et le soulèvement ne sert pas à la corriger — il sert à augmenter `T`. Cela ne réfute
> pas le F61, qui parle de `T → ∞`. Et **deux items ouverts se ferment seuls** à 100 passes :
> `to_nothing` **1 → 0** (la dernière écharde du F38) et `Unresolved` **3 → 0**.
> ⚠️ **La fraction de lacs, elle, ne bouge presque pas** (24,05 → 17,87 %, encore **18×** le
> critère) : c'est un bilan hydrique que l'incision ne lit jamais (F39, F60).
>
> **Le devis D dit que le correctif de la brèche est le mauvais** : un plancher à `sea + 0,5 m`
> empêche **152 des 169** tranchées d'atteindre leur cible, qui retombent alors sur le FILL et
> deviennent des **mares plates** — l'option que le F13 avait proposée et que l'auteur avait
> refusée. Et **tous les gardes de la fonction resteraient verts** (ils vérifient la monotonie, que
> le fill garantit). **Le correctif est le bloc C**, pas D : à 100 passes la brèche noie 2 880
> cellules au lieu de 10 375, sans toucher à `breach_monotone`.
>
> **Nommé, non fait** : les 7 299 passes restantes jusqu'à Courant 1 (~47 h hors session), qui
> diraient si le critère côtier est tenu par la seule intégration.

> ⛔ **ET LE CRITÈRE CÔTIER A SA PREMIÈRE MESURE SUR LE PRODUIT PROMU : Δ(éperons ≥ 1 km) = +115**
> à 8192² (135 livrés contre **20** pré-incision — l'instrument reproduit au chiffre les « vingt »
> du F80). Contre les **+1 828** du F80 sur le même instrument, c'est une réduction **×16** ; ce
> n'est pas zéro, et la jambe 4 (« l'incision ne doit ajouter aucun éperon ») **n'est pas tenue**.
> Invisible à 2048² (Δ = 0), pire à 4096² (+184) qu'à 8192².
>
> **Trois corrections au dossier, toutes à mes tours.** Le **F87-B est retiré** : le fond à
> −8,67 m est l'œuvre de **`breach_monotone`**, pas de la diffusion — le champ **érodé** est à
> **+0,50 m**, soit exactement le plancher du F83, et à **−1,00 m** plancher coupé : **le plancher
> mord et fait ce pour quoi il a été promu**. Le **F88-D1 mélangeait deux populations** dans sa
> ligne du fond (passes 0-1 érodé, passe 2 breaché) : la part de l'incision est 12,13 m dont la
> passe 1 porte **97,4 %**. Et **« Courant 1353 » est la valeur à 2048²** — le livré à 8192² est
> **3699**, et cinq lignes du dossier attachent le chiffre grossier au monde fin.
>
> **Et une passe de plus supprime le trou** : à `iterations = 3` le fond livré de (3487, 5902) est
> **+0,50 m** au lieu de −8,67 m.
>
> **Nommé, non corrigé** : le soulèvement ; la cellule `to_nothing` (3203, 3586) à −1,000 m,
> `water_class 2`, dans aucun plan d'eau, qui reçoit **19,5 m³/s** (famille du F38, un survivant
> des 68) ; la conséquence du F66 (`A` sous-estimé de ~88 % en aval de chaque dépression, **à chaque
> itération** — une boucle temporelle la multiplie par N) ; et les cibles des critères B1/B2 qui
> restent **PROXY** faute de source primaire, là où `K` en a trois avec leurs pages.

## Reste en attente (hors closures)

- **Séparation tronc/affluents** façon Azgaar : un tronc nommé portant un profil ordonné
  source→embouchure, chaque affluent maximal comme cours d'eau propre avec son lien
  {rejoint le tronc S au point P}. La structure actuelle (liste plate de tronçons,
  confluence implicite dans le graphe) est décrite à l'ADR.
- **Érosion glaciaire** (`glacial.rs`, vide, M6) : c'est elle qui donnerait l'Écosse **pour
  de bon** — de vrais surcreusements plutôt qu'un bruit qui les imite.
- **Côté Living Landz**, exporté mais ignoré : `width_m` (rivières rendues à largeur
  constante), `lake_type` (tous les plans d'eau aplatis en une classe), et le biome
  `Wetland` qui a désormais une source de données mais aucun consommateur.

## Les deux règles de méthode acquises

Elles ont chacune coûté plusieurs passes et valent pour tous les chantiers à venir.

**Mesurer en configuration de production, aux deux résolutions.** Un compteur validé à
2048² ne dit rien de 8192² — les micro-cuvettes survivent quand une cellule est 16 fois
plus petite.

**Un terrain reconstruit n'est pas le produit.** Il a trompé le diagnostic **six fois**,
y compris un `merge_verify` masqué par la quantification u16. L'export fait foi.