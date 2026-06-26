# Conception — cache des steps de génération C1 (repartir de l'amont inchangé)

**Statut :** conception (décisions arrêtées avant code). Fil *infrastructure*, distinct de #165.
**Principe :** re-calculabilité appliquée à l'exécution — si une couche amont n'a pas
changé, son résultat est réutilisable ; sinon il (et tout son aval) se recalcule.

---

## 0. La chaîne et son profil de coût (mesuré / inféré)

```
init_c1_state_phase_2_r7 (64²)          ~ s        cheap        [état 64²]
run_with_closures (64², 300 steps)      ~ s        cheap        [état 64²]
upscale_from_c1 (2048²)                 ~ MINUTES  LE COÛT      [eroded HD]  ← isostasy+Stein-Stein+FBM+ÉROSION
c1_drainage (2048²)                     ~ s–10s    modéré       [drainage]
c1_climate (T, P) (2048²)               ~ s        cheap        [climat]
c1_biomes (2048²)                       ~ ms       trivial      [biomes]
```

Ancrage : le run hypsométrie #165 = 663 s pour 6 seeds ≈ **110 s/seed à 2048²**, dominé par
le FBM + l'érosion hydraulique HD *à l'intérieur de `upscale_from_c1`*. Le reste est du
per-cellule (climat, biomes) ou un parcours de grille (drainage).

---

## Décision 1 — quels steps cacher

**Verdict : cache PRINCIPAL = le terrain érodé (sortie de `upscale_from_c1`). Cache
SECONDAIRE = le drainage. Climat + biomes = recalculés à chaque fois (cheap).**

| frontière | produit | coût recalcul | taille (2048² f32) | cache ? |
|---|---|---|---|---|
| état 64² | `C1State` post-`run_with_closures` | cheap | ~0.5 MB | optionnel (petit ; cache-le pour redémarrer l'aval sans la boucle tectonique) |
| **eroded HD** | `upscale_from_c1().heightmap` | **minutes** | 16 MB | **OUI — le step qui justifie tout** |
| drainage | `c1_drainage()` | modéré | ~50 MB (filled+dir+accum+basins) | OUI — modéré, et la calib climat/biome n'y touche pas |
| climat | `c1_climate()` (T, P) | cheap | 32 MB | non — recalcul depuis l'eroded |
| biomes | `c1_biomes()` | trivial | 4 MB | non |

### ⚠️ Caveat honnête sur la boucle « fix relief » (le cas d'usage cité)
Le déficit de plaine basse se corrige **en amont** du FBM+érosion (plancher continental =
`IsostasyConfig` / l'altitude de production). Donc **chaque itération qui change le relief
invalide légitimement le cache eroded et refait l'érosion** — le cache n'accélère PAS le step
qu'on est en train de modifier. C'est l'invalidation qui fait son travail, pas un échec.

Ce que le cache accélère réellement :
- **(a) la calibration EN AVAL d'un relief figé** — faire varier `PrecipParams`, le seuil
  Whittaker, la latitude, ou re-mesurer l'hypsométrie / la température / les biomes **sans
  refaire l'érosion**. C'est le gros gain, et c'est fréquent dans la boucle relief→climat→biomes
  (les sous-itérations climat/biomes à relief constant).
- **(b)** ré-examiner plusieurs diagnostics aval sur le même terrain.
- **(c)** les runs multi-seed où certains seeds n'ont pas changé.

Le cadrage juste : **le cache matérialise la chaîne de dépendances** — seul le step modifié
*et son aval* recalculent. Toucher le climat ⇒ climat+biomes seulement (énorme). Toucher le
relief ⇒ érosion+aval (incompressible). Le **cache par-frontière** (l'option « idéale »)
est ce qui rend ce comportement propre quel que soit le step touché.

---

## Décision 2 — la clé d'invalidation (le cœur)

**La clé encode TOUTES les entrées du step caché ET de tout son amont. On sérialise en JSON
canonique puis on hashe (content-addressing), parce qu'aucun config ne dérive `Hash`/`Eq`
(champs flottants) — c'est déjà l'approche de `metadata.json`.**

### Composants de la clé du cache `eroded`
- `seed: u64`
- `grid` / `target_size` (2048)
- `Phase2InitParams` (init de l'état)
- `C1TimeLoopConfig` : `n_steps`, `rigid_continental_crust`, `dx`, `dy`, `drainage_max_distance`
- identité des `C1Closures` actives (la boucle tectonique)
- `IsostasyConfig`  *(le plancher continental du fix relief vit ici si c'est un paramètre)*
- `SteinSteinParams`  *(⚠ ajouter `#[derive(Serialize, Deserialize)]`)*
- `FbmUpscaleConfig` — **embarque `Option<ErosionConfig>`** ⇒ couvre FBM + érosion d'un coup
- **`algo_version` par step** : un entier explicite, bumpé à CHAQUE changement de CODE du step
  ou d'un step amont *sans* changement de paramètre. C'est le levier anti-péremption.
  `CARGO_PKG_VERSION` est trop grossier (ne bouge pas par édition de code).

### Content-addressing (le mécanisme qui rend la lecture-périmée IMPOSSIBLE)
```
digest = hash(canonical_json(key))           // sha256/blake3 hex tronqué
path   = cache_dir / format!("{step}_{digest}.raw")
```
Une clé différente ⇒ nom de fichier différent ⇒ **miss automatique** (jamais de lecture d'un
fichier dont la clé ne correspond pas). Pas de « comparer puis décider » faillible : le
système de fichiers fait l'index. Un sidecar `{step}_{digest}.json` stocke la clé lisible
(debug) + dims + longueur (backstop de cohérence, déjà vérifié par `load_raw`).

### Chaînage des clés = matérialisation de la chaîne de dépendances
Chaque step compose la clé de son parent :
```
k_state   = root(seed, init, run_cfg, closures, ALGO_TECTONICS)
k_eroded  = k_state.then(iso, ss, upscale_cfg /*incl. erosion*/, ALGO_UPSCALE_EROSION)
k_drain   = k_eroded.then(drainage_cfg, ALGO_DRAINAGE)
```
⇒ **amont changé ⇒ digest amont change ⇒ tous les digests aval changent ⇒ tout l'aval miss**,
automatiquement. C'est exactement le comportement voulu.

### Le cas du fix relief, traité explicitement
- Fix = **paramètre** (p.ex. un `continental_floor_m` dans `IsostasyConfig`) → il est dans la
  clé → digest change → **miss du cache eroded → l'érosion refait sur le NOUVEAU relief**. ✓
- Fix = **changement de code** (sans paramètre, p.ex. réécrire `c1_production_altitude`) →
  **bumper `ALGO_UPSCALE_EROSION`** (et tout step dont le code change) → digest change → miss.
  Discipline contractuelle : *toute* modif de code d'un step caché ou de son amont exige de
  bumper son `algo_version`. (Un hash du source automatiserait ça mais est lourd ; l'entier
  explicite + checklist est le choix pragmatique.)

### Anti-piège (explicite)
**En cas de doute, invalider.** Un recalcul redondant coûte des minutes ; une lecture périmée
coûte la validation de résultats FAUX. Si on hésite à savoir si un changement touche un step
caché → bumper la version. Le content-addressing penche déjà du bon côté (clé trop sensible →
recalcul, jamais périmé).

---

## Décision 3 — format binaire + rapport à l'export §9

**Format : raw f32 little-endian, row-major, sans header — exactement `GridF32::save_raw` /
`load_raw` + `save_raw_f32/u8/u32` qui EXISTENT déjà dans `export/mod.rs`. On réutilise, on
n'écrit pas un second sérialiseur.**

### Partagé OU distinct ? → **codec partagé, politiques distinctes**
Le cache et l'export §9 sont les *mêmes octets sur disque* mais deux *politiques* opposées :

| | cache (interne) | export §9 (livrable) |
|---|---|---|
| but | reprise de calcul | produit Living Landz |
| nommage | content-addressed `{step}_{digest}.raw` | répertoire humain `seed2048/` |
| métadonnées | sidecar mince (clé + dims) | `PipelineMetadata` complet (contrat §9.3) |
| durée de vie | jetable, gitignored (`.ymir_cache/`) | stable, versionné |
| index | le digest | le nom + metadata.json |

**Action concrète :** activer le module déjà anticipé `// pub mod raw;` ([export/mod.rs:19](crates/ymir-core/src/export/mod.rs#L19))
— y extraire le codec bas-niveau (`save_raw_f32`/`load_raw_f32`/u8/u32 + `GridF32::save_raw`/
`load_raw`), consommé par **les deux** couches. Le cache ajoute son sidecar mince ; l'export
garde `PipelineMetadata`. **Codec partagé, pas de sur-couplage** (le cache ne dépend pas de
`PipelineMetadata`, l'export ne dépend pas du digest).

---

## Décision 4 — scope du maillon (insertion + API)

### Où il s'insère
Un wrapper de cache fin autour de chaque frontière cachable, dans le **chemin produit C1**
(le harness / futur entrypoint produit qui enchaîne run→upscale→drainage→climat→biomes).
Les fonctions `upscale_from_c1`, `c1_drainage`, … **ne changent pas** : le cache les *enrobe*.

### API — lire-si-présent-valide / calculer-et-écrire-sinon
```rust
/// Content-addressed : HIT si le .raw du digest existe et que dims/len collent ; sinon
/// compute() puis écrit. Jamais de lecture d'un fichier de clé non concordante.
fn cached<T: RawCodec>(dir: &Path, step: &str, key: &CacheKey, compute: impl FnOnce() -> T) -> T;
```
```rust
let dir = cache_dir(); // .ymir_cache/, gitignored

let k0 = CacheKey::root(seed, &init, &run_cfg, &closures, ALGO_TECTONICS);
let state  = cached(&dir, "state",  &k0, || run_c1(seed, &init, &run_cfg, &closures)); // optionnel

let k1 = k0.then(&iso, &ss, &upscale_cfg, ALGO_UPSCALE_EROSION);   // upscale_cfg embarque l'érosion
let eroded = cached(&dir, "eroded", &k1, || upscale_from_c1(&state, &iso, &ss, &seed, &upscale_cfg).heightmap); // LE gain

let k2 = k1.then(&drainage_cfg, ALGO_DRAINAGE);
let drainage = cached(&dir, "drainage", &k2, || c1_drainage(&eroded, &drainage_cfg, &ss)); // secondaire

// aval cheap — recalcul systématique, pas de cache
let clim   = c1_climate(&eroded, &ss, lat, &precip);
let biomes = c1_biomes(&eroded, &clim);
```

`RawCodec` = trait implémenté via le codec partagé (Décision 3) pour `GridF32` et les sorties
multi-champ (`C1DrainageResult` → plusieurs `.raw` + un sidecar de dims). `CacheKey::then`
compose le digest parent (Décision 2).

---

## Prérequis avant code (checklist)
1. `SteinSteinParams` + `PrecipParams` : ajouter `#[derive(Serialize, Deserialize)]`.
2. Activer `export::raw`, y déplacer le codec bas-niveau (sans changer le comportement de l'export §9).
3. Définir les constantes `ALGO_*` (une par step caché) + documenter la règle « code modifié ⇒ bump ».
4. `.ymir_cache/` dans `.gitignore`.
5. Choisir le hash (blake3 conseillé : rapide ; sinon sha256). Digest hex tronqué (16 hex suffisent).

## Ce qu'on NE fait pas
- Pas de second sérialiseur binaire (réutiliser `export::raw`).
- Pas de cache du climat/biomes (cheap ; et c'est souvent ce qu'on calibre).
- Pas de couplage cache↔`PipelineMetadata` (sidecar mince distinct).
- Pas de clé basée sur `std::hash::Hash` (flottants) — JSON canonique → digest.
