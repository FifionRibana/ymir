# Verdict — sous-diagnostic du plancher : réglage ou closure manquante ?

**Probe :** `c1_closure_morphology::probe_floor_thickness` (`#[ignore]`).
**Run :** tectonique 64², seeds {42 (avec cratons), 1337 (sans craton)}, rifting ON/OFF,
décompose coarse-iso vs FBM 2048². ~2 s.
**Question :** le socle haut (né haut, cf [VERDICT floor_diagnosis](../floor_diagnosis/VERDICT.md))
est-il un RÉGLAGE trop haut (baisser un plancher → tout descend) ou une CLOSURE
d'amincissement MANQUANTE (croûte épaissie quelque part, amincie nulle part → tout au
niveau « épais »/haut) ?

## ÉTAPE 1 — le doc C1 décrit-il l'amincissement ?

OUI, abondamment (`docs/c1_lightweight_dynamic_tectonics.md`) :
- **Rifting / thinning (§4.5, §5.2)** : closure dédiée `closures/rifting/`, formulée comme
  un **sink de S̃ sur les frontières continentales divergentes** (« negative s_field source,
  mirror of Davis-Suppe positive on convergent »), ancrée **McKenzie 1978 / Buck 1991**
  (β = 1.4, 30 % d'amincissement). Statut doc : Track D (#132) ⏳ « in progress ».
- **Foreland basin (flexure, Beaumont 1981)** : explicitement « produces fertile plains —
  prime city-placement terrain ». Statut : **§7.3 Phase 3, NON implémenté.**
- Code : la closure rifting EST câblée ([time_loop.rs] `apply_rifting_thinning`) et **activée
  par défaut** (`C1Closures::default()`, donc dans le chemin `c165_eroded`).

→ Le mécanisme d'amincissement n'est pas une découverte : il est prévu ET partiellement
codé. La question devient « tire-t-il ? ».

## ÉTAPE 2 — mesures

| Mesure (coarse 64², continental) | Seed 42 (cratons) | Seed 1337 (0 craton) |
|---|---|---|
| S̃ non-craton : moyenne / std / **CV** | 0.41 / 0.29 / **0.69** | 0.49 / 0.34 / **0.70** |
| S̃ non-craton min..max | 0.01 .. 1.61 | 0.00 .. 1.86 |
| % non-craton à S̃≈1.0 (|S̃-1|<0.1) | 3 % | 6 % |
| craton : fraction de la terre émergée | **58 %** | 0 % |
| **émergé craton** : médiane, frac<500m | **1814 m, 5 %** | — |
| **émergé non-craton** : médiane, frac<500m | **429 m, 54 %** | 773 m, 37 % |
| rifting : mass_removed / splits | 35.1 / **0** | 35.2 / **0** |
| rifting OFF : Δ S̃ moyenne / Δ frac<500m | +0.006 / 0.0 pt | +0.006 / +0.7 pt |
| DECOMPOSE coarse→FBM (frac<500m, terre) | 25.3 % → 27.3 % | — |

CSV : `floor_thickness.csv`.

## Verdict : RÉGLAGE de la buoyancy CRATONIQUE (ni closure manquante, ni plancher global)

1. **Croûte VARIÉE, pas uniforme** (CV 0.70 ; S̃ de 0.01 à 1.86 ; seulement 3-6 % à la
   baseline 1.0) → RÉFUTE « socle uniforme = closure d'amincissement manquante ». Il y a
   déjà énormément de croûte fine.
2. **L'amincissement tire déjà** — mais via l'**érosion macro** (boucle tectonique), pas le
   rifting : le non-craton est aminci à S̃≈0.41-0.49. Le **rifting est inactif** (0 split sur
   300 pas, le désactiver ne change rien : Δ S̃ +0.006) — c'est le risque « Track C : events
   trop rares » du doc §7.2 réalisé. Activer/renforcer le rifting n'est PAS le levier.
3. **Le socle haut = les CRATONS.** Décomposition de la terre émergée (seed 42) : cratons =
   **58 % de la terre à médiane 1814 m (5 % <500 m)** = le bourrelet 1000-2000 m+. Le
   non-craton est **Terre-like (médiane 429 m, 54 % <500 m)**. Sans cratons (seed 1337) le
   plancher est bien plus bas (médiane 773 m). Le déficit de plaine basse vient des cratons
   qui dominent l'hypsométrie émergée et l'écrasent vers le haut.
4. **Fixé par l'isostasie, pas le FBM** : coarse 1247 m → +FBM 1176 m (frac<500m
   25.3 %→27.3 %). Le FBM préserve ; le bouton est dans l'isostasie.

## Levier du fix (à valider, pas encore appliqué)

La buoyancy cratonique : `IsostasyConfig::c1_default` a DÉJÀ `craton_rho_crust: Some(2900)`
(croûte cratonique plus dense → flotte plus bas) — mais c'est **insuffisant** : les cratons
culminent encore à médiane 1814 m. Réel : les boucliers cratoniques usés sont à ~300-600 m,
et la majeure partie de l'aire cratonique est **plateforme** (couverte de sédiments, plaine
basse), pas bouclier exposé haut. Le modèle élève TOUTE la croûte cratonique en haut bouclier.

→ Le fix est un **RÉGLAGE** : abaisser la freeboard cratonique (`craton_rho_crust` plus dense
vers le mantle, et/ou baisser `craton_thickness_ratio`). Cela descend le socle haut SANS
noyer les plaines non-craton (déjà correctes) → crée le contraste (cratons = boucliers
modérés, non-craton = plaines, orogènes = hauts). PAS une closure d'amincissement (déjà
présente/active via érosion), PAS un plancher continental global (noierait les plaines OK).

## Nuance (résiduel secondaire)

Le non-craton émergé est légèrement haut sur seed 1337 (médiane 773 m, 37 % <500 m). Après le
fix craton, un résidu mineur peut subsister (mapping `land_ref_thickness`/`max_elevation_m`,
ou la dynamique « l'érosion noie les bas au lieu de bâtir des plaines », cf le défaut
`sea_level=0.1` du floor_diagnosis). À re-mesurer après le levier craton. Premier ordre
= cratons.
