# Handoff — C1 climat (#165) + cache, reprise sur nouveau PC

**Pour :** un agent d'exécution Claude Code SANS contexte préalable (nouvelle session,
nouveau PC). Ce document + le repo suffisent à reprendre. Factuel : état du code,
conventions, prochain travail. Pas d'historique de raisonnement.

**Date de rédaction :** 2026-06-25.

---

## 1. État du projet

- **Repo :** Ymir — générateur de continents physiquement fondé, Rust, workspace 2 crates
  (`ymir-core` = logique ; `ymir-viz` = Bevy). Licence propriétaire.
- **Branche active :** `165-c1-climate-precipitation-field-and-climate-derived-layers`.
  Convention : **une branche par issue**, tous les maillons de l'issue dessus.
  L'intégration se fait sur `YmirDvlp0.X`, **pas** `main`.
- **Issues :**
  - **#155 (relief, MERGÉ)** — contrats de coordonnées, drainage, upscale HD, érosion HD,
    cratons. Le pipeline relief→altitude métrique existe et est stable.
  - **#165 (climat, EN COURS)** — champ de précipitation + couches dérivées du climat
    (température, biomes Whittaker). C'est le travail courant.
- **HEAD :** voir `git log --oneline -10`. Les derniers commits #165 portent sur le climat
  (re-ancrage T_sea, fond frontal abaissé, verdict-grid) PUIS le **maillon cache** (commit
  de sauvegarde de fin de session — s'il n'est pas dans `git log`, c'est qu'il est resté en
  working-tree non commité ; vérifier `git status`).
- ⚠️ **Le compte de commits « vs milestone » vient TOUJOURS du REMOTE** (`git fetch` puis
  comparer), jamais du local périmé. Vérifier à CHAQUE PR.

---

## 2. Maillons livrés — le code qui EXISTE (ne pas réinventer)

### Climat — `crates/ymir-core/src/climate/`
- **`mod.rs`** : `c1_climate(heightmap, ss, latitude_deg, params) -> ClimateResult { temperature, precipitation }` (fonction pure du relief + latitude) ; `c1_biomes(heightmap, &climate) -> Vec<Biome>`.
- **`temperature.rs`** : `compute_temperature`. T_sea ré-ancré sur la vraie courbe
  T-vs-latitude (`(lat/90)^1.8`, pas `sin²`) ; lapse altitudinal **6.5 °C/km**.
  `SEA_LEVEL_NORM = 0.5`.
- **`precipitation.rs`** : `compute_precipitation_with_budget` / `precip_mm_per_year` ;
  transport conservatif de l'humidité (Clausius-Clapeyron / Smith-Barstad) + **fond frontal
  `k_frontal`** (terme de base advecté). Le fond frontal a été **abaissé 1281→450 mm** pour
  démasquer le gradient orographique (il saturait l'intérieur). `PrecipParams` (défaut 45°
  westerlies, vent W→E). `SEA_LEVEL_NORM = 0.5`.
- **Biomes** : classification Whittaker sur seuils thermiques −5 / +5 / +20 °C ; légende
  bandes de précip mm/an.

### Cache des steps de génération (infra, ce session) — voir `docs/design/c1_generation_cache.md`
- **`crates/ymir-core/src/cache.rs`** : `CacheKey` (JSON canonique trié via `BTreeMap` +
  chaînage `derived_from` ; `with` serde / `with_debug` pour les configs sans serde),
  `cached(dir, step, key, compute)` **content-addressed** (`digest = blake3(canonical_json)`,
  `path = {step}_{digest}.raw` + sidecar `.json` ; lecture de périmé STRUCTURELLEMENT
  impossible). Trait `RawCodec` (impl pour `GridF32`). Constantes **`ALGO_TECTONICS` /
  `ALGO_UPSCALE_EROSION` / `ALGO_DRAINAGE`** avec règle de bump commentée.
- **`crates/ymir-core/src/export/raw.rs`** : codec binaire bas niveau f32/u8/u32 LE sans
  header — **PARTAGÉ** par le cache ET l'export §9 (`export/mod.rs`). `GridF32::save_raw/
  load_raw` y délèguent. Un seul sérialiseur.
- **`crates/ymir-core/src/tectonics_c1/cached_product.rs`** :
  `cached_c1_eroded(dir, seed, grid, init, run, closures, ss, upscale_cfg) -> GridF32`
  (enrobe `init → run_with_closures → upscale_from_c1`, le ~108 s d'érosion HD) ;
  `cached_c1_drainage(...)` (chaîné sur la clé érodée ; `impl RawCodec for C1DrainageResult`
  = rasters en `.raw` + réseau/lacs/stats en sidecar JSON) ; `tectonic_key` / `eroded_key` /
  `drainage_key`. Les fonctions de génération `upscale_from_c1`/`c1_drainage` sont
  **INCHANGÉES** (pures) — le cache est un wrapper opt-in.
- **`.ymir_cache/`** gitignored, jetable, content-addressed.
- **Mesuré :** MISS 108 s → HIT ~10 ms (×~10000), round-trip byte-identique, transparence
  CSV bout-en-bout vérifiée (hypsométrie identique cache vs direct).

### Probes #165 — `crates/ymir-core/tests/c1_closure_morphology.rs` (`#[ignore]`)
**Toutes** passent par le helper unique **`c165_eroded(seed, &iso)`** (clé identique → SHARING :
l'érosion est calculée une fois par seed, réutilisée par toutes les analyses aval) :
- `probe_climate_acceptance` (conservation/magnitude/structure), `probe_climate_maps`
  (cartes T/P), `probe_climate_biomes` (histogrammes + cartes), `probe_climate_verdict_grid`
  (relief+biomes+précip+température+transect), `probe_frontal_oro_ratio` (diagnostic),
  `probe_hypsometry` (distribution d'altitude vs Terre).
- `probe_cache_loop_gain` : démonstrateur MISS/HIT/invalidation (`#[ignore]`).
- Invocation : `cargo test --release -p ymir-core --test c1_closure_morphology <nom> -- --ignored --nocapture`.
- Sorties : `docs/reports/c1_continental_buoyancy/closure_morphology/` (+ scripts Python de
  grille `image_grid_*.py`, `hypsometry_plot.py` à la racine).
- Probes 512²/1024² (`probe_precip_structure`, `probe_climate_control_survey`,
  `probe_oro_depletion`) **NON** branchées (résolution ≠ 2048², ne partagent pas le terrain ;
  laissées en appel direct ; 5 warnings cosmétiques `unused grid` dans les probes converties,
  sans conséquence).

### Contrats de coordonnées — `crates/ymir-core/src/tectonics_c1/production_upscale.rs`
- **Vertical :** `c1_altitude_norm_to_metres(norm, ss) = (norm − 0.5) · 11300` (mer 0.5 → 0 m ;
  défaut `ALTITUDE_NORM_HALF_RANGE=1.13`, `depth_scale_m=5000`). Inverse `c1_metres_to_altitude_norm`.
- **Horizontal :** `C1_DOMAIN_KM = 1024` ; `c1_km_per_cell(grid)` = `1024/grid` (0.5 km @2048²) ;
  `c1_cell_area_km2`.
- **Ancre physique unique :** `S̃ = 2.0 = 70 km = h_eq = plafond EH` (équilibre-hauteur).

---

## 3. Conventions du projet (discipline d'exécution)

- **Gated / byte-identical :** toute extension est un opt-in (`Option`, `#[serde(default)]`,
  flag `enabled`) dont l'état par défaut/None reproduit le comportement antérieur **au bit
  près**. Une régression byte-identity se prouve en énumérant TOUS les échecs préexistants
  (lib + intégration), pas en `tail`.
- **Un sujet par commit.** Messages : `TYPE : Issue #N <sujet> [#N]`, **sans** trailer
  Co-Authored-By (voir les commits récents pour le format).
- **Ancré, pas tuné :** chaque paramètre repose sur une grandeur physique réelle (lapse
  6.5 °C/km, seuils Whittaker, navigabilité km², Stein-Stein GDH1…). Ne pas re-tuner pour
  rattraper une cible quand seule une partie de la physique est active — documenter +
  anticiper.
- **Cache sur nouveau PC :** `.ymir_cache/` sera **VIDE** → premier run par seed = MISS
  (~108 s) qui le re-remplit (normal, jetable, reconstructible à l'identique — transparence
  prouvée). ⚠️ **Bump `ALGO_*`** à TOUT changement de code d'un step caché (érosion,
  upscale, drainage, production-altitude) : le content-addressing protège les changements
  de **config** (infaillible) mais le **code** repose sur ce bump manuel — l'oublier = servir
  un terrain périmé (bug silencieux). En cas de doute, invalider (recalcul = minutes ;
  périmé = résultats faux validés, bien pire).
- **Format de sortie :** raw binary §9 (`export/raw.rs`, f32 LE sans header).
- **Build :** toujours `--release` (érosion 10-20× plus lente en debug). `cargo test
  --workspace`, `cargo fmt --all` (max_width 100), `cargo clippy --workspace` (note : une
  erreur clippy PRÉEXISTANTE `Ord`/`PartialOrd` dans `tectonics_c1/distance_field.rs`, sans
  rapport avec le cache).
- Windows / PowerShell : pour l'édition de fichiers non-ASCII en masse, `Set-Content
  -Encoding utf8` corrompt (BOM + mojibake) — utiliser `[System.IO.File]::WriteAllText`
  avec `UTF8Encoding($false)`.

---

## 4. Prochain travail (en attente)

### IMMÉDIAT — le fix relief (le levier courant)
**Constat (mesuré, `probe_hypsometry`) :** les continents manquent de **plaine basse** —
~34 % des terres sous 500 m vs ~52 % réel ; la masse s'accumule en **1000–2000 m** ; ce
n'est **PAS** un excès de montagnes (3000+ sous le réel). Conséquence aval : trop peu de
**forêt tempérée** (les intérieurs hauts sont thermiquement boréaux/taïga).

**Fix = abaisser le PLANCHER continental** (ramener de la plaine basse), **pas raboter les
montagnes**. Le paramètre vivra dans **`IsostasyConfig`** (chemin `crate::tectonics::isostasy`,
réexposé par `IsostasyConfig::c1_default()`). ⚠️ **`continental_floor_m` N'EXISTE PAS encore** —
c'est le nom anticipé du futur paramètre (référencé seulement dans des commentaires de test).
Il faut le créer (gated : `Option`, None = byte-identique).

**Diagnostiquer AVANT de tuner :** établir POURQUOI le plancher est à 1000–2000 m
(flottabilité de base de la croûte continentale dans l'isostasie ? érosion qui n'abaisse pas
assez l'intérieur ? ancre verticale ?), pas raboter à l'aveugle.

**La boucle est maintenant gratuite via le cache :** changer `iso` → la clé invalide
(prouvé : digest différent) → l'érosion se refait **une fois** sur le nouveau relief →
calibration climat/biomes/hypsométrie **gratuite** (HIT) sur ce relief figé. Itérer ainsi.

### Différés #165
- Jonctions latitudinales / slider (terme de base latitudinal — partiellement fait via le
  fond frontal) ; viz log-scale des précip.

### Fils distincts (hors #165)
- **viz-rewire** : brancher la viz interactive sur `c1_drainage` (chemin produit) au lieu du
  cache d'érosion v2.
- **export §9** : le codec existe (`export/raw.rs`) ; reste le `metadata.json` complet + le
  rebouclage vers Living Landz.
- **sous-marin** : océans plats (bathymétrie / détail submarin).

---

## 5. Verdict en cours (état du jugement)

- Le **climat → biomes est VALIDÉ EN TANT QUE CLIMAT** : structure spatiale correcte
  (ombre pluviométrique W→E), précip et thermique cohérents ; la taïga dominante est un
  **froid d'altitude légitime** (vérifié par overlay biome × température × altitude — la
  taïga est bien en altitude, lapse correct).
- **MAIS il manque des forêts tempérées**, et la cause est le **RELIEF** (intérieurs trop
  hauts), **pas le climat**. Le climat répond correctement à un relief qui manque de plaine
  basse.
- → **Le fix relief (§4 IMMÉDIAT) est le prochain levier.** Ne pas re-toucher le climat pour
  « fabriquer » de la forêt : ce serait masquer le vrai problème (le relief).

---

## Démarrage rapide (nouvelle session)
1. `git fetch && git status && git log --oneline -10` — confirmer la branche et que le commit
   cache est présent (sinon il est en working-tree).
2. `cargo build --release && cargo test -p ymir-core --lib` — doit être vert (incl.
   `cache::tests`, `tectonics_c1::cached_product::tests`).
3. Lire `docs/design/c1_generation_cache.md` (décisions cache) puis ce handoff §4.
4. Commencer par DIAGNOSTIQUER le plancher continental (§4 IMMÉDIAT) avant tout tuning.
