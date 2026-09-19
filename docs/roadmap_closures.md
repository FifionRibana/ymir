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
> 2. **Le relevel du F85, promu** — ⛔ **ET IL A EMBARQUÉ UNE RÉGRESSION (Finding 91)** :
>    avec le relevel, **2 des 20 composantes sous-marines encloses du livré ne portent AUCUNE
>    étiquette de plan d'eau** (13 corps) ; sans lui, **0** (17 corps). La fusion est par
>    **réciprocité, pas par adjacence** : la classe s'inonde d'un seul fond vers une seule surface
>    et marque `seen` sur toutes les régions qu'elle absorbe — une région au-dessus de cette surface
>    n'est jamais traitée ni revendiquée. **34,8 m³/s de rivière et de déversoir terminent sur ces
>    cellules**, et la plus grande des deux fait **14,89 km²**, plus que huit des treize corps
>    restants. Les gardes du F86 vérifiaient le **champ** (inchangé — le relevel ne déplace aucune
>    hauteur) et `exorheic_lakes_missing_outlet` (0) ; **aucun n'énumère l'invariant du F38 : toute
>    composante `wc == 2` est couverte par un plan d'eau.** Fermé en 2024 pour 68 embouchures
>    orphelines, jamais épinglé. — `ALGO_DRAINAGE` 6 → 7, `ALGO_HD_DRAINAGE` 8 → 9.
>    ✅ **CORRIGÉ AU FINDING 92 (variante β, choix de l'auteur).**
>    `MergedUnionRelevel::separate_unclaimed_regions`, défaut `true` — un correctif, pas une porte.
>    Deux lignes : l'absorption pousse la cellule dans le `comp` de la classe **sans la marquer
>    `seen`**, et le balayage extérieur saute les cellules qu'une inondation a déjà revendiquées.
>    Une région atteignable fusionne encore ; une région inatteignable devient **son propre corps**.
>    La chaîne est intacte (`own_label` pilote toujours `extra_inflow`), donc la convergence du F85
>    et la fermeture du F86 ne dépendaient pas de l'absorption des cellules.
>
>    | | garde F92-B | corps | passes | `to_own_body` | `to_nothing` |
>    |---|---|---|---|---|---|
>    | **ROUGE** (`separate_unclaimed_regions: false`) | **2 non couvertes** : `(3656,3691)` 6 245 c., `(3203,3586)` 1 c. | 13 | 8, converge | 0 | **1** |
>    | **VERT** (livré) | **0** | **15** | **8, converge** | **0** | **0** |
>    | contrôle (relevel OFF) | 0 | 17 | **16, ne converge pas** | **3 / 501,99 m³/s** | 0 |
>
>    **Le garde d'abord, rouge, puis le correctif, vert** — le seul ordre qui prouve quelque chose.
>    Et β ferme aussi le **`to_nothing` 1 → 0** : les 19,5 m³/s du bassin 1000004 atterrissent enfin
>    dans un plan d'eau réel (l'item nommé-non-corrigé du F89-C3).
>
>    ⚠️ **Le prix, dit et non enfoui** : la composante de 14,89 km² devient **son propre lac de
>    173,19 km², niveau 76,25 m, profondeur 83,20 m, `Endorheic`** — parce que seule, elle se remplit
>    jusqu'à son propre seuil. L'eau sous-marine passe de **4 775,9 à 4 949,1 km² (+3,6 %)** et la
>    fraction de lacs du F89-C1 va de 24,05 % à ≈ **24,69 %** : la mauvaise direction pour le critère
>    < 1 %, pour une bonne raison — cette eau était là et n'était pas comptée. **Ma prédiction
>    (+14,89 km²) est fausse d'un ordre de grandeur.**
>
>    ⛔ **Et le garde du F38 est enfin épinglé** : `uncovered_below_sea_components` + une
>    **assertion de production** en fin de `assemble_hd_drainage` (pas `#[cfg(test)]`, F72), avec son
>    **contrôle négatif** en test unitaire permanent. **Règle de méthode 13 : un invariant qu'on
>    ferme reçoit une assertion permanente, ou il n'est pas fermé.**
>
>    Et le viz dit enfin ce que la condition voit : « → cuvette sous-marine SANS exutoire tracé »
>    remplace « puits sous-marin (évaporatif) », qui affirmait un processus que
>    `if !sea && wc == 2` ne peut pas lire (F39 : en humide `net_evap = 0 ⇒ a_eq = ∞`). Le point
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

## ⛔ LE TEMPS INTÉGRÉ : LES CANYONS PARTENT, LA CÔTE PASSE, LE RELIEF DESCEND ENCORE (Finding 95)

Le F44 item 2 mesuré au bout, à `k_time` **tenu à 9 000** — même budget d'érosion, seule la
distribution change. `cfl_iterations` **7 399**, Courant livré **3 699**, `dt_max` **2,703·10⁻⁴ an**
(et non « 60 à 120 ans » : `dt_max = cell_m / c` avec `c = 180 618 m/yr` — au K livré un knickpoint
traverse 400 km en deux ans).

| | **2 (livré)** | **10** | **100** | **300** | **1 000** |
|---|---|---|---|---|---|
| Courant | 3 699 | 740 | 74 | 24,7 | **7,4** |
| fraction de lacs | **24,05 %** | 18,99 | 17,87 | 15,92 | **18,05 %** |
| **classe canyon** | **16 / 53** | 5 / 30 | 3 / 22 | **0 / 24** | **0 / 24** |
| `to_nothing` / `Unresolved` | 1 / 3 | 1 / 1 | 0 / 0 | 0 / 1 | **0 / 0** |
| relief apparié p50 | 424,2 m | +18,7 % | **+19,1 %** | +15,1 % | **+10,2 %** |
| **Δ(≥ 1 km), champ ÉRODÉ** | −1 | −2 | −3 | −2 | **+0** |
| Δ(≥ 1 km), brèché + u16 | +115 | +86 | +27 | +29 | **+26** |

> ⛔ **Trois des quatre critères sont atteints, et le levier était numérique pour les trois.** La
> **classe canyon tombe à 0** et y reste ; **Δ(≥ 1 km) sur le champ érodé atteint exactement 0** —
> le stade sur lequel les F76–83 ont écrit le critère (F90-A) ; `to_nothing` **et** `Unresolved`
> tombent à **0**. Le résidu de brèche du F90, lui, **survit à l'intégration** (+26 éperons sur le
> brèché) : son attribution tient et c'est un défaut séparé.
>
> ⛔ **Mais le relief n'a pas convergé : il culmine vers Courant 74 puis redescend** — +18,7 → +19,1
> → +15,1 → **+10,2 %**, soit **0,2 point** hors du ±10 % déclaré, et **par le HAUT**, pas en
> planation. **Aucune des trois lignes de décision ne s'applique telle qu'écrite** : le quantité
> n'est pas un état, elle est en descente. La coupe médiane par cellule, elle, monte de façon
> monotone (69,5 → 88,3 m) pendant que l'altitude p50 baisse — la signature d'un champ travaillé
> plus uniformément, non plus fort.
>
> **Donc la décision honnête est plus étroite que la table ne le demandait** : les canyons et le
> critère côtier **n'avaient pas besoin du soulèvement**, ils avaient besoin de l'intégration.
> Reste **une** grandeur, avec une trajectoire connue et un coût connu pour la trancher : les
> **6 399 passes restantes jusqu'à Courant 1**, soit **36 h d'un cœur éveillé**.
>
> ⚠️ **Et le ±10 % est suspect pour une raison neuve** : il a été posé symétrique en supposant que
> l'intégration ne pouvait que retirer du relief. La mesure dit qu'elle en **retient 19 % de plus**
> d'abord, puis les rend. **Une bande symétrique autour du champ livré teste la conformité à une
> référence sous-intégrée** — celle que ce tour montre être l'artefact.
>
> ⚠️ **La fraction de lacs ne suit pas l'intégration** (24,05 → 18,99 → 17,87 → 15,92 → **18,05 %**,
> non monotone) et reste **18×** le critère : c'est un bilan hydrique que l'incision ne lit jamais
> (F39, F60). Prédit des deux côtés, vérifié.
>
> ⛔ **Et le devis du F89-D est amendé : la boucle est SÉRIELLE.** Parallélisme mesuré **0,99 cœur
> sur 24**, soutenu. Les 7 399 passes ne font pas « 1,9 jour de machine » mais **≈ 42 h d'un seul
> cœur**, et **ajouter des cœurs n'y change rien** tel que le code est écrit (le priority-flood de
> `compute_flow` est un tas binaire, l'accumulation D8 une descente ordonnée). Coût propre mesuré à
> 1 000 passes : **20 387,7 s = 20,4 s/passe**, ce qui confirme les 22,1 s/passe du F89-D.
> ⚠️ Le « 28 335 s » imprimé au jalon 300 est **contaminé** : `Instant::elapsed()` a compté la veille
> de la machine — ne pas le citer.

## ⛔ LE CLIP N'EST PAS LE DÉFAUT — l'inconsistance rivière/lac l'est (Finding 93)

Le log `rivers clipped to lakes: 14708 -> 7265` n'est **pas** un compte de suppressions : le F20 a
bâti le clip pour **découper** chaque parent en ses suites maximales de points hors lac.

| sort de la suite | compte | part de 7 780 |
|---|---|---|
| **GARDÉE** (≥ 2 points hors lac) | **7 265** | **93,4 %** |
| rejetée : suite de **< 2 points** | **515** | 6,6 % — **417,1 m³/s** de débit parent |
| rejetée : exutoire d'un lac **ENDORHÉIQUE** | **0** | 0,0 % |
| rejetée : exutoire d'une cuvette **SOUS-MARINE** | **0** | 0,0 % |
| parents n'émettant **aucune** suite | **7 448** | leur polyligne entière est dans des lacs |

Le recensement reproduit le clip **exactement** (7 265 = le compte gardé), ce qui est son propre
contrôle.

> ⛔ **ZÉRO faux positif du type supposé** : aucune suite n'est coupée parce que sa destination est
> un plan d'eau couvert — les deux clauses qui le pourraient tirent **0 fois**. **La réduction, ce
> sont 7 448 parents entièrement à l'intérieur des lacs** : c'est le défaut du F20 lui-même
> (*« 44 segments traversant un polygone de lac »* à 2048²) devenu, à 8192², **44,6 % de tous les
> points de rivière dans de l'eau stagnante** — les lacs couvrent 6 660 km² et le réseau est tracé
> sur le champ **breaché**, qui draine toutes les cuvettes.
>
> ⛔ **Le fleuve fantôme n'est pas clippé.** La plus grande portion de ≤ 4 points post-clip est la
> 2612 : **2 022 km², 19,2 m³/s, 4 points = 0,15 km**, `downstream = Some(2613)`, et **son parent
> 5166 faisait déjà 4 points — le clip en a gardé 100 %**, avec **0** point dans un lac. Avec la
> signification géographique (ratio 7,5, aires ×56,25) cela fait **113 762 km², 1 081 m³/s,
> 1,10 km** — à **1,2 %** près des « 1 095 m³/s, 115 181 km², 1 km » de l'écran. **#381 est une
> portion de quatre points d'un tronc correctement chaîné, affichée en unités signifiées.**
>
> **Les deux vrais ingrédients étaient déjà au dossier.** L'aire est **héritée** du parent — défaut
> du **F42**, dont le remède est écrit et chiffré depuis 2024 (*« quelques lignes dans
> `clip_rivers_to_lakes` : lire l'aire à la cellule aval de la suite elle-même au lieu de l'hériter,
> SAUF pour une suite d'exutoire exorhéique où le F22 exige l'héritage, plus un bump
> `ALGO_DRAINAGE`. Coût faible, effet cosmétique »*) et **toujours pas fait**. Et le panneau montre
> une **portion** là où l'œil attend un **système** — la sémantique du **F45**, dont le chaînage
> exclut **délibérément** les cuvettes sous-marines (*« son écoulement est un `Spillway` typé sans
> hiérarchie »*), ce qui est aussi la raison des connexions « coupées » de #64/#212.
>
> **Et la famille tracé-de-lac est TROIS mécanismes, pas un** : les lacs **22** et **41** ont une
> échappée à **0,00 m** et **0,01 m** sous leur col — le port du F37c rencontre un rebord sans
> gradient et la descente y retombe ; le lac **55** voit son tracé **réussir** (échappée 3,58 m plus
> bas, déjà dans un autre plan d'eau) et ne manque que d'un **tronçon exporté** (test réseau du F86,
> l'écart ×5,25 du L7031) ; et **1000001** (660,79 km², 45,04 m³/s) et **1000002** (173,19 km²,
> 21,57 m³/s) n'ont **aucun déversoir émis** — et **aucune suite de carve rejetée** non plus, donc
> il n'y avait rien à couper.
>
> ⚠️ **Deux erreurs d'instrument de ma part, corrigées et consignées** : un rattachement de parent
> par égalité d'aire (ambigu — plusieurs parents partagent 2 022 km²) qui donnait la conclusion
> inverse, et l'**oubli de l'ajout des déversoirs** avant le test du F86, qui faisait lire **9** lacs
> `Unresolved` et 7 relabels au lieu de **3** et **1**. Le F92 n'est pas régressé.
>
> **Nommé, non fait** : les 515 suites de moins de deux points (**417,1 m³/s**), et **le 44,6 %** —
> le chiffre qui compte, et qu'aucun tour n'a encore demandé de changer.

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

## ⛔ LA FRANGE CÔTIÈRE EST DE LA SUR-INCISION LE LONG DES CHENAUX — ET LE PLANCHER QUI LE PROUVE IMPRIME D8 DANS LE TERRAIN (Finding 96)

Recherche d'une closure bon marché pour l'âge du continent, après le F95 (6 210 s d'intégration).
**Deux coutures ajoutées en production, additives et éteintes par défaut** (`depression_floor:
false`, `incision_floor: None`), épinglées par un test d'inertie byte-à-byte (règle 13). Tout le
reste en banc.

**La construction χ est neuve au dossier** (`chi`, `Flint`, `Perron`, `Royden`, `stationnaire` :
**0 occurrence**). Elle est **valide comme profil et invalide comme champ** :

| | |
|---|---|
| l'intégrale `χ = ∫(A₀/A)^{m/n} dx`, une descente ordonnée | **0,37 s** sur 11 327 660 cellules, 0 non atteinte |
| **contrôle de Flint** (666 293 points) | `θ = 0,4998` pour 0,5 · l'ordonnée redonne `(U/K)^{1/n}` ⇒ **passe** |
| cuvettes | **0** contre 1 240 677 au livré |
| **mais** σ local 3×3 | **44,73 m** contre 6,28 au livré ⇒ **7,1× plus RUGUEUX** |
| part de terre > 30° | **79,07 %** contre 28,53 |

> ⛔ **La cause n'est pas un bug, elle est arithmétique** : le gradient du profil est
> `(U/K)^{1/n}(A₀/A)^{m/n}`, soit **52,9° à une cellule de bassin** avec A₀ = 1 km². J'appliquais la
> loi fluviale **trois décades sous son propre seuil de chenal** (`A_c = 0,1 km²`, F43/F65, ancré).
> Brider A à `A_c` ramène ce gradient à 12,8° et la part > 30° à 60,87 % — **un tiers de l'excès,
> pas plus.**
>
> ⛔ **Le reste a une autre cause, que la mesure permet de nommer : χ est DISCONTINU EN TRAVERS DES
> LIGNES DE PARTAGE.** Deux cellules voisines de part et d'autre d'une crête ont des intégrales
> qui diffèrent de centaines de mètres, et `z = z_base + C·χ` en fait une falaise (à C = 0,0719 et
> 48,8 m de cellule, 45° ne demande que Δχ = 679 m). **χ n'a aucun mécanisme pour faire s'accorder
> deux bassins voisins sur leur crête commune** — la propriété même qui fait du χ-mapping un
> diagnostic de capture dans la littérature. Donc χ fournit une **borne**, jamais une surface — ce
> qui est exactement où pointait le F6 : *« la tectonique d'Ymir a déjà fait le soulèvement, le
> correctif est de LIMITER l'incision totale, pas d'ajouter U »*. ⛔ **Et « une borne, jamais une
> surface » ne suffit pas non plus : voir la validation visuelle plus bas — une borne construite
> sur une intégrale de chemin D8 imprime les axes D8 dans le terrain.**

**La table de décision, en coût marginal (la porte « 30 s » recale le champ livré lui-même, 45,5 s) :**

| | **livré** | **A1** exclusion | **B3** plancher χ | **A1+B3** | **ORACLE** (300 p.) |
|---|---|---|---|---|---|
| classe canyon, **en taux** | 29,6 % | **21,9 %** | 29,8 % | 26,8 % | **0 %** |
| **côte, brèché + u16** | **Δ +115** | Δ +104 | **Δ +3** | Δ +4 | Δ +29 |
| relief apparié p50 | 424,2 m | **447,6 (−8,3 %)** | 535,7 (+9,8 %) | 546,6 (+12,0 %) | **488,1 m** |
| fraction de lacs | 24,70 % | 22,67 | 20,19 | **19,72** | **15,92** |
| `to_nothing` / `Unresolved` | 0 / 3 | 0 / **1** | 0 / 4 | 0 / 4 | **0 / 1** |
| **coût marginal** | — | +7,2 s | **+2,2 s** | +3,0 s | **6 210 s** |

> ⛔ **Les deux défauts sont séparables et n'ont pas le même propriétaire.** La **frange côtière est
> de la sur-incision le long des chenaux** : B3 la ferme (+115 → **+3**) et A1 ne l'effleure pas
> (+104). La **classe canyon est de l'incision dans les cuvettes fermées** : A1 la mord (29,6 →
> **21,9 %**) et B3 ne la bouge pas (29,8 %) — χ est construit sur le champ pré-incision **brèché**,
> où les futures empreintes de lac ont été percées, donc à l'intérieur χ est bas et le plancher
> autorise l'entaille qu'il devrait arrêter.
>
> ⛔ **Et elles ne se composent pas** : l'union prend la côte de B3 mais seulement une partie des
> canyons d'A1 (26,8 %), et c'est le seul candidat hors ±10 %. Les deux closures se disputent le
> même budget de relief.
>
> ⛔ **B3 bat l'oracle sur le critère côtier d'un ordre de grandeur, pour 0,04 % de son coût** :
> Δ +3 contre Δ +29 à 300 passes (+26 à 1 000). **Le résidu de brèche du F90, que le F95 avait
> explicitement laissé ouvert comme « survivant à l'intégration », se ferme ici — et par le côté
> incision, pas par le côté brèche.**
>
> ⚠️ **Ce qu'aucun candidat ne fait** : atteindre la classe canyon 0 de l'oracle, ni sa fraction de
> lacs (meilleur 19,72 %, encore **20×** le critère — le bilan hydrique que l'incision ne lit jamais,
> F39/F60/F95). **La réponse est partielle et doit être dite ainsi.**
>
> ⛔ **ET LA VALIDATION VISUELLE RÉFUTE B3 COMME CHAMP LIVRABLE — ce qu'aucun chiffre de la table
> ne voyait.** Sur une tuile de 1024², le champ livré est dendritique ; le champ B3 porte des
> **striations horizontales et verticales et des terrasses en blocs qui suivent les axes D8**, et
> son quadrant inférieur gauche a **perdu son réseau** au profit d'une rampe lisse. La cause est
> structurelle : **χ est une intégrale de chemin D8**, donc ses iso-niveaux héritent de la
> quantification à huit directions, et brider le terrain dessus **imprime ces axes dans le relief**.
> L'union est **visuellement identique** : A1 n'y change rien.
>
> ⚠️ **Et le point de méthode est la moitié la plus tranchante.** La table dit que la texture de B3
> est **meilleure** que celle du livré (p90 38,1° contre 44,6 ; > 30° 23,07 % contre 28,53). Les
> deux chiffres sont justes, et les deux sont des **statistiques de MAGNITUDE isotropes** : elles
> mesurent la raideur, jamais les directions ni la survie du réseau. **Un champ rayé et
> dé-dendritisé marque MIEUX sur toutes**, parce qu'il est vraiment plus lisse — mal lisse. Toutes
> les mesures de texture du dossier (σ local, quantiles de pente, part > 30°, percentiles
> hypsométriques) partagent cet aveuglement. **Un regard sur une tuile a trouvé ce que neuf tours
> de chiffres n'avaient pas soulevé.**
>
> ⇒ **Donc B3 est un INSTRUMENT D'ATTRIBUTION valide et une CLOSURE invalide.** Ce qu'il prouve
> tient (la frange côtière est de la sur-incision le long des chenaux) ; ce qu'il ne peut pas, c'est
> être livré. **Le candidat suivant est moins cher que χ et n'a besoin d'aucune intégrale de
> chemin** : le F6 lu littéralement, un **plafond sur l'ENTAILLE en fonction de l'aire drainée**,
> `h_r_eff = max(h_r, min(h_pre − c·A^q, h_o))` — même couture (`incision_floor` existe déjà et est
> déjà inerte), **isotrope par construction**, un ou deux paramètres libres. **Non mesuré ce
> tour**, posé comme candidat et pas comme résultat.

**Deux voies fermées, chiffrées :**

- ⛔ **`talus_passes` ne supprime pas les murs.** 4 → 16 → 64 passes : p90 44,6 → 37,1 → **33,7°**,
  mais **la part > 30° MONTE** (28,82 → 32,96 %) et **le maximum ne bouge pas** (85,2 → 86,3°), pour
  **268,7 s** — 6,5× l'incision livrée. Le F6 avait raison : *« les faces quasi-verticales sont les
  murs des fentes de 1 px »*, et un balayage local ne relaxe pas une fente d'une cellule sans la
  combler.
- ⛔ **Un champ sans dépression ne rachète pas la brèche** : **−6,8 / −8,7 / −11,3 s** sur trois
  runs, pas les ~155 s espérés. Le coût de la brèche est le **parcours**, qui visite chaque cellule
  qu'il y ait quelque chose à creuser ou non.
- ⛔ **Le dépôt intra-boucle (transport-limited) est hors domaine pour un tour « bon marché »**, sur
  les chiffres du dossier lui-même : *« le plus gros changement (un budget sédimentaire) »* (L487),
  *« consigné comme spécification, pas construit »* (L2514), et surtout *« cette loi doit être posée
  dans un régime où le modèle intègre, **ce que la configuration livrée ne fait pas** »* (L4129) —
  à Courant 3 699. **Sa propre précondition est le régime intégrant que ce tour refuse.**

**Règle de méthode 14, gagnée ce tour** : *une closure qui BORNE l'érosion doit publier l'érosion
qu'elle a permise.* La première forme d'A1 marquait **4 des 6 critères du F95** — canyons 0/29,
lacs 17,51 %, côte +2, `to_nothing` 0 — **en éteignant l'incision partout** (`min(filled, h_o) =
h_o` puisque `filled ≥ h`). Seules deux grandeurs non-critères l'ont dénoncée : `median cut 0,0 m`
et un p50 égal à celui du champ non incisé. **Le dossier avait déjà enregistré ce mode d'échec une
fois, en production** (`a_c_slope_law` : `A_c × 100` *« retire le peigne — et éteint l'incision avec
lui »*). Un mode d'échec qui revient sur deux paramètres différents est une règle, pas une anecdote.

## ⛔ LE PEIGNE D8 EST CELUI DU PRODUIT LIVRÉ — ET TOUTE CLOSURE BON MARCHÉ PAIE SON CRITÈRE AVEC LA MOITIÉ DE L'ÉROSION (Finding 97)

Tour d'instrument. **Une couture ajoutée en production** (`slope_floor_uk`, additive, `None` par
défaut, byte-identique, épinglée par le test d'inertie — règle 13). Le reste en banc.

**L'instrument**, neuf au dossier (`structure tensor` : 0 occurrence) : tenseur de structure sur
fenêtres disjointes, et surtout **R8 = |⟨e^{8iθ}⟩|**, la huitième harmonique circulaire de
l'orientation. ⛔ **La statistique d'isotropie que le dossier possède déjà, `axis R` (F84), est la
DEUXIÈME harmonique — elle est aveugle par symétrie à un motif D8 à quatre axes** : {0°, 45°, 90°,
135°} donne ⟨e^{2iθ}⟩ = 0, zéro par symétrie et non par isotropie. Échelle fixée aux deux bouts
avant tout jugement : **0,0402** sur un synthétique prouvé isotrope, **0,9999** sur un rayé à 45° et
14 cellules — la longueur d'onde même que le F76 a mesurée pour la fourrure.

| R8, tuile (5120, 3072) — celle que l'œil avait vue au F96 | w=4 | w=8 | w=16 |
|---|---|---|---|
| **livré** | **0,0710** | **0,0990** | **0,1382** |
| **B3, le plancher χ** | **0,0079** | 0,0313 | 0,0779 |
| pré-incision | 0,0020 | 0,0071 | 0,0324 |

> ⛔ **L'ATTRIBUTION VISUELLE DU F96 EST RETIRÉE.** J'avais écrit *« χ imprime les axes D8 dans le
> relief »*. Sur la tuile même, **B3 est 9× MOINS aligné que le livré**, et l'ordre est
> pré-incision (0,002) < B3 (0,008) < **livré (0,071)**. Les striations appartiennent au **produit
> livré** : c'est **le peigne du F76** — périodique, parallèle, non-chenal, 14 cellules — et le
> plancher χ en **retire** la majeure partie. Ce que l'œil avait vu, c'était l'artefact préexistant
> rendu lisible parce que le plancher avait enlevé la texture dendritique qui le masquait.
>
> Ce qui reste vrai est plus étroit, et la mesure le donne : **l'IMPRINT** (B3 − livré) a R8 =
> **0,22**, contre 0,092 pour le livré. **L'action du plancher est anisotrope ; le champ qu'elle
> produit l'est moins que celui qu'elle corrige.** Le F96 affirmait la première et en déduisait la
> seconde, ce qui ne suit pas.
>
> ⚠️ **Et l'hypothèse de repli que j'avais posée — des bandes de quantification u8 — est fausse
> aussi** : la tuile couvre 1–3 831 m, les **237 codes sur 237** sont utilisés, une rampe locale
> n'est pas plus fine. C'est un **hillshade** qui a tranché, et il montre le champ livré saturé de
> striations strictement verticales et horizontales. **Première image du peigne du F76 sur le
> terrain** — le F76 l'avait mesuré sur le trait de côte et lui avait donné une longueur d'onde ;
> personne n'était allé le regarder à l'intérieur des terres.

**La table, avec la colonne qui manquait, et la règle 14 qui mord à son premier emploi :**

| | livré | **B1** plafond d'aire | **A1** exclusion | **B2** pente locale | **A1+B2** | ORACLE 300 p. |
|---|---|---|---|---|---|---|
| classe canyon (taux) | 29,6 % | **0,0 %** | 21,9 % | 26,3 % | **3,4 %** | 0,0 % |
| **côte brèchée + u16** | **+115** | +58 | +104 | **+3** | **+3** | +29 |
| relief apparié vs oracle | −13,1 % | +17,1 % | **−8,3 %** | +10,1 % | +13,2 % | — |
| **érosion permise (règle 14)** | 100 % | ⛔ **43,7 %** | **88,4 %** | ⛔ 59,0 % | ⛔ 51,2 % | **> 100 %** |
| **R8 (w16)** | 0,0924 | 0,1048 | 0,0870 | **0,0367** | **0,0304** | n/a |
| fraction de lacs | 24,70 % | 22,09 | 22,67 | **18,90** | 19,21 | 15,92 |

> ⛔ **B2 récupère exactement la côte de χ — Δ +3 — sans aucune grandeur accumulée.** La contrainte
> différentielle s'intègre d'elle-même le long de la descente ordonnée (le solveur visite les
> récepteurs avant les donneurs). **Et B2 est le champ le plus isotrope mesuré, R8 = 0,0367, sous
> le livré ET sous le synthétique prouvé isotrope.** Le clamp gradé n'avait jamais été testé : le
> F76 dit *« Not run: B (the seam and its variants) »*, sa règle d'arrêt 1 avait tiré avant.
>
> ⛔ **B1 atteint la classe canyon de l'oracle — 0 sur 37, exactement 0 % — en supprimant 56 % de
> l'érosion.** C'est le mode d'échec pour lequel la règle 14 a été écrite au F96, revenu un tour
> plus tard sur un autre paramètre. Contre le seuil du tour (> 70 %), **seul A1 passe (88,4 %)**, et
> A1 ne ferme aucun des deux défauts.
>
> ⛔ **DONC, à ce `k_time`, les défauts ET l'érosion sont la même grandeur.** Aucune borne sur
> l'incision n'est assez forte pour fermer la côte ou les canyons sans retirer la moitié de la
> coupe. **L'oracle fait l'inverse** : le F95 mesure sa coupe médiane qui MONTE (69,5 → 88,3 m), donc
> son « érosion permise » est **au-dessus** de 100 % — il retire les défauts en **redistribuant** le
> même budget sur plus de passes, pas en le retenant.
>
> ⇒ **Le mur, nommé.** La seule route mesurée qui retire les défauts *sans* retirer l'érosion est
> l'intégration à `k_time` tenu (F44 item 2, F95), à 6 399 passes. Une closure qui **redistribue**
> au lieu de retenir serait un terme de dépôt transport-limited — que le F96-A3 a écarté sur les
> chiffres du dossier (L4129 : il doit être posé dans un régime où le modèle intègre, et la
> configuration livrée tourne à Courant 3 699). **Ces deux faits ensemble sont le mur contre lequel
> ce chantier se trouve maintenant, et il vaut mieux le dire que le contourner une fois de plus.**

⚠️ **La colonne coût de ce run est contaminée** (le livré y prend 116,5 s contre 45,5 s au F96 :
construit en premier, en contention avec une compilation). Les quatre candidats tiennent en **1,8 s
les uns des autres**, donc ils sont équivalents entre eux ; le marginal contre le livré n'est pas
lisible ici. Les mesures du F96 tiennent : +2,2 s pour le plancher χ, +7,2 s pour A1.

**Règle 14b amendée** : une revendication de texture demande **un rendu CALIBRÉ (hillshade, pas une
rampe d'altitude) ET une statistique directionnelle étalonnée aux deux bouts** — et quand les deux
divergent, aucune ne gagne d'office, on cherche laquelle ment. Le F97 a essayé quatre statistiques ;
**une seule marche** : ni `axis R` (aveugle par symétrie), ni la cohérence C (plus haute pour le
champ lisse non incisé que pour le livré), ni l'entropie d'orientation (0,996–0,999 partout).

## ⛔ LE MUR N'EN EST PAS UN : LE CADRAN TIENT TOUS LES DÉFAUTS ET RAMÈNE LE RELIEF DANS LA BANDE DE L'ORACLE — MAIS LA CLOSURE TAXE L'ÉROSION DE MOITIÉ PARTOUT, ET LE CADRAN SATURE (Finding 98)

Aucun changement de production. Neuf champs, six chaînes de critères, 1 h 50 d'un cœur.

**Le contrôle tire, et fort.** Monter `k_time` sans closure aggrave **tous** les défauts :

| | canyons | côte brèchée | relief vs oracle | R8 |
|---|---|---|---|---|
| livré ×1 | 29,6 % | +115 | −13,1 % | 0,0924 |
| **CONTRÔLE livré ×1,5** | **40,4 %** | **+147** | **−20,3 %** | **0,1017** |

**Et la closure tient à travers le cadran** — c'est le résultat du tour :

| A1+B2 | canyons | côte | relief vs oracle | lacs % | **érosion permise (même `k_time`)** | R8 |
|---|---|---|---|---|---|---|
| ×1 | 3,4 % | **+3** | +13,2 % | 19,21 | 51,1 % | 0,0303 |
| ×1,5 | **3,4 %** | **+3** | **+10,3 %** | 19,44 | 50,0 % | 0,0319 |
| ×2 | 6,7 % | **+3** | **+8,6 % — DANS ±10 %** | 19,60 | 49,2 % | 0,0317 |

> ⛔ **LE MUR NOMMÉ AU F97 EST RETIRÉ.** Le F97 écrivait *« à ce `k_time`, les défauts et l'érosion
> sont la même grandeur »*. **Ils ne le sont pas** : A1+B2 à ×2 tient chacun de ses défauts à sa
> valeur de ×1 **et** entre dans la bande ±10 % de l'oracle, ce qu'aucun candidat n'avait fait. La
> closure **découple le où du combien**, et `k_time` redevient le cadran de l'âge.
>
> ⚠️ **Ce qui survit du mur est plus étroit et stable : la taxe.** Contre le livré **au même
> `k_time`**, l'érosion permise vaut **51,1 → 50,0 → 49,2 %** — plate. Contre le livré ×1 elle
> *monte* (51,1 → 56,3 → 59,4 %), uniquement parce que le budget a grossi : **c'est le piège
> comptable, et les deux dénominateurs sont imprimés pour qu'on ne puisse pas lire le flatteur.**
>
> ⛔ **Et c'est la SATURATION qui rend la taxe chère.** `k_time` ×2 n'achète que **×1,21** de coupe
> réelle (2,241e9 → 2,705e9 m). Rattraper la moitié manquante par le seul cadran demanderait, à ce
> taux de change, de l'ordre de **×6 à ×8**, à un Courant déjà 3 699 fois au-delà de la borne
> d'onde. **Le mur n'est pas « défauts contre érosion » : c'est le TAUX DE CHANGE du cadran.**

**La carte de l'érosion retirée — et la réponse est « partout ».**

| | retiré | (i) chenal sous χ | (ii) cuvette | **(iii) reste** | **porte : AUCUNE** |
|---|---|---|---|---|---|
| A1+B2 | 1,101e9 m | 7,8 % | 3,3 % | **88,9 %** | **91,5 %** |
| B1 plafond d'aire | 1,265e9 m | 7,2 % | 2,0 % | **90,7 %** | **94,0 %** |
| B2 seul | 9,369e8 m | 8,4 % | 1,8 % | **89,8 %** | **93,2 %** |

> ⛔ **Le discriminant du tour ne discrimine pas — et le mien non plus.** Les trois candidats
> tiennent en **1,2 point** sur (iii) et **2,5** sur « aucune porte ». Ma partition par les portes
> ne les sépare pas mieux que celle par χ.
>
> ⛔ **Normalisé par le nombre de cellules, le résultat est plus fort que ce que les deux
> partitions cherchaient.** Les classes de porte couvrent 6,4 % (cuvette), 5,8 % (le plancher
> mordrait) et 87,9 % (aucune) des terres, pour 3,4 / 5,0 / 91,5 % du volume retiré : **intensité
> par cellule 0,53 · 0,87 · 1,04**. **Le retrait est UNIFORME sur le paysage, et plutôt plus FAIBLE
> là où les closures agissent.** En clair : le livré coupe **198 m par cellule** en moyenne, A1+B2
> en retire **97 m, par cellule, partout**.
>
> **Le mécanisme est la transmission.** Les deux closures tiennent le CHENAL en hauteur ; le chenal
> est le niveau de base local de chaque versant au-dessus, donc la diffusion et le talus — qui
> agissent sur toutes les cellules et portent l'essentiel du volume — déplacent beaucoup moins de
> matière. **Une closure branchée sur 12 % des cellules divise par deux l'érosion sur 100 %
> d'entre elles.** Ce n'est ni « empêcher » ni « supprimer » : c'est relever le plancher sous tout
> le paysage.

**Le bloc C ne mesurait pas ce qu'il croyait**, et le code le disait d'avance : `anti-peigne` a
**0 occurrence** au dossier ; le levier du F78 est le **clamp de plateau**, il tourne **après**
l'incision (`:451` puis `:506`), et le **F79 l'a ramené de 20 m à 1 m en production** avec un bump
`ALGO_UPSCALE_EROSION`. Mesuré quand même, comme contrôle négatif : **identique sur chaque
colonne** (canyons 29,6 %, côte +115, relief 424,2 m, érosion 100 %, R8 0,0924 → 0,0926). La côte
ne bouge pas non plus, et ce n'est pas une contradiction avec les 43 % de fourrure du F78 : le F78
compte les **éperons ≥ 2 cellules** sur le champ érodé, cette table compte les **≥ 1 km sur le
brèché u16** que le F90 attribue à la traînée de brèche. **Séparation d'échelle, pas désaccord.**

**Règle 14 amendée, deux fois** : (a) le dénominateur doit être la référence **au même budget** —
sinon une closure « permet plus d'érosion » dès qu'on monte le cadran tout en en permettant la même
fraction ; (b) la carte doit être une **intensité par cellule**, pas une part de volume — une part
de volume ne fait que dire où sont les cellules.

⚠️ **La colonne coût est inutilisable pour le deuxième tour de suite**, et c'est désormais un point
de méthode : les temps de construction tombent de façon monotone au fil d'un run (146,4 → 128,5 →
74,1 → 64,8 s pour un travail de taille identique). **Ils ne sont comparables que mesurés côte à
côte.**
