# Grille d'intégration — la chaîne complète sur un terrain partagé (viz)

**Probe :** `c1_closure_morphology::probe_integration_grid` (`#[ignore]`).
**Assemblage :** `python image_grid_integration.py`.

## But

La **vue d'ensemble** : tout le produit actuel visible en UN endroit, par seed, sur le MÊME
terrain caché — l'instrument de **validation visuelle globale** avant l'export Living Landz
(la cohérence inter-maillons que les probes par maillon ne voient pas). C'est de l'assemblage
de l'existant (rendus climat/biomes/drainage déjà écrits), pas de la réinvention.

## Ce qu'elle montre (6 colonnes, 6 seeds, échelle commune)

1. **relief + bathymétrie** — terre hypsométrique (vert→brun→blanc) + océan par profondeur
   (plateau cyan → abysse navy) : le profil plateau/talus/abysse du sous-marin se lit.
2. **drainage** — `c1_drainage` : rivières navigables (small-boat/barge/ship) + lacs typés
   (exo/endo) en overlay sur le relief assombri.
3. **précip** (mm/an, bandes communes), 4. **température** (°C, seuils Whittaker), 5. **biomes**
   (palette commune), 6. **transect** O→E (altitude + précip + bande biome).

## Le cache paie ici

Toutes les vues lisent le MÊME `c165_eroded` (HIT) et le MÊME `c165_drainage` (caché, chaîné
sur la clé érodée) → **une seule érosion par seed alimente toutes les colonnes**. La grille
6 seeds se génère en ~7 s (érodé en HIT depuis le cache, seul le drainage calculé une fois).
Avant, le rendu drainage reconstruisait l'érosion en direct (lent + incohérent) ;
`probe_c1_drainage_acceptance` est aussi rebranchée sur le cache (1,1 s vs minutes).

## Validation d'intégration (cohérence inter-maillons)

- **Bathymétrie autour des côtes** : halo plateau cyan → abysse navy sur chaque continent. ✓
- **Rivières → mer** : réseaux dendritiques drainant vers les côtes ; lacs **exorhéiques**
  dominants (189/308/142 exo, 0 endo) → l'eau rejoint la mer. ✓
- **Biomes cohérents** relief × climat (forêt en plaine humide, toundra/steppe ailleurs). ✓
- Ensemble cohérent, aucune incohérence flagrante entre maillons. ✓

## Suite tracké

L'**interactif `ymir-viz`** est sur le vieux chemin v2 (pipeline phases Stokes + stubs CLI/BIO
qui prétendent faussement que le climat est un placeholder — faux depuis #165) ou le moteur
c1 tectonique-coarse. À rediriger vers le produit C1 (`cached_c1_eroded` HD → `cached_c1_drainage`
→ `c1_climate` → `c1_biomes`) quand on voudra un viewer live. DIFFÉRÉ (confort, pas bloquant
pour l'export). Les PNG (tuiles + grille) sont régénérables (non commités).
