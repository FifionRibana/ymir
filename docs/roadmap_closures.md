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
- **fond** — **remplacer progressivement le FBM par les closures manquantes**. Chaque
  closure ajoutée réduit la part de bruit nécessaire, jusqu'à ce qu'il redevienne ce
  qu'il aurait dû être : une simple **graine de brisure de symétrie**.

Règle qui accompagne cette trajectoire : à mesure que les closures arrivent, **le FBM doit
rétrécir, pas cohabiter**. Deux sources de détail — une causale, une aléatoire — qui se
superposent sans se coordonner reproduiraient le problème des stries, qui a coûté
plusieurs passes à diagnostiquer.

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

### 3. Hétérogénéité lithologique — 🔬 IMPLÉMENTÉ, en attente de validation visuelle (voir ADR 0001, section « C-3 »)

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

### 4. Érosion côtière

**Pourquoi en dernier :** elle travaille sur un littoral dont la forme dépend de tout
l'amont. La faire après lithologie et volcanisme lui donne un trait de côte déjà
structuré à sculpter, plutôt qu'à reprendre ensuite.

**Ce qu'elle apporte :** falaises, plages, plateformes d'abrasion — la morphologie
littorale qui a préoccupé plusieurs fois pendant la phase hydrologie (embouchures
perchées, falaise contre plage au rendu).

**Note :** `erosion/coastal.rs` est un fichier vide marqué M5, comme `thermal.rs`,
`aeolian.rs` et `glacial.rs`.

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