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

### Ordre imposé par ces dépendances

`échelle de temps explicite (F44)` → `terme de versant transport-limited (F43 méc. 1)` →
`hypsométrie convergente` → puis, redevenus calibrables : les seuils de navigabilité en m³/s,
le cadran de H-2, et la calibration de tête de chenal aux deux résolutions.

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