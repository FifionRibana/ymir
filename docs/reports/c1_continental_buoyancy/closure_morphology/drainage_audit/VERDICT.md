# Verdict — les lacs sont-ils sur-remplis ? (drainage fill-and-spill)

**Probe :** `c1_closure_morphology::probe_drainage_lake_audit` (`#[ignore]`).
**Run :** 6 seeds, terrain + drainage cachés partagés (`c165_eroded` + `c165_drainage`) +
`c1_climate`. Déclencheur : la grille d'intégration montrait beaucoup de gros lacs, tous les
bassins ~remplis, tous avec exutoire (0 endorhéique) — la signature d'un fill-and-spill.

## Mesure

| seed | lacs | aire lacs % terre | exo/endo | gros lacs (1000-5000 / 5000+ km²) | aire en semi-aride <500mm |
|---|---|---|---|---|---|
| 42 | 189 | 19.8 % | 189 / 0 | 10 / 1 | 76 % |
| 99 | 67 | 7.9 % | 67 / 0 | 2 / 0 | 65 % |
| 1337 | 203 | 14.6 % | 203 / 0 | 10 / 1 | 86 % |
| 4138 | 117 | 26.8 % | 117 / 0 | 6 / 1 | 79 % |
| 1988 | 308 | 16.5 % | 308 / 0 | 4 / 3 | 75 % |
| 2026 | 142 | 18.3 % | 142 / 0 | 10 / 4 | 77 % |
| **agg** | — | **18.0 %** | **1026 / 0** | — | **78 %** |

## Verdict : SUR-REMPLI — fill-and-spill géométrique, l'évaporation est ignorée

L'hypothèse (3 observations = 1 mécanisme) est **confirmée sur les 5 angles** :

1. **Aire de lacs = 18 % de la terre** (par seed 8-27 %) vs Terre **~2 %** → **~9× trop**.
2. **Distribution de taille** : beaucoup de petits (<100 km²) MAIS une queue de **gros lacs**
   (1000-5000 km² : ~10/seed ; 5000+ : 1-4/seed) qui DOMINE l'aire — ce sont les bassins
   remplis jusqu'au débordement.
3. **L'algorithme = priority-flood FILL-AND-SPILL, confirmé par le CODE.** `c1_drainage`
   ([drainage.rs]) ne prend QUE le heightmap (pas de champ précip/évaporation). Le docstring
   `LakeType` le dit lui-même : *« Priority-flood pit-filling routes every depression to an
   overflow sill that ultimately spills to the ocean, so PURE GEOMETRY yields Exorheic for
   all lakes. A TRUE endorheic lake — where evaporation balances inflow — is a CLIMATE
   phenomenon, not geometry; this field will carry Endorheic only once a hydroclimate layer
   couples. »* C'est un placeholder géométrique assumé.
4. **Lien climat (le test décisif) : 78 % de l'aire de lacs est en semi-aride (<500 mm)** —
   précisément là où la Terre a ses bassins FERMÉS (Caspienne, Tchad, Aral). Mais TOUS
   débordent (exorhéiques). **Le climat est IGNORÉ** : un bassin aride déborde exactement
   comme un bassin humide. (0 % en <250 mm car le climat à 45° pose un fond frontal ~450 mm ;
   à 30° — le désert subtropical — la part aride serait plus forte. Incohérence inter-maillon
   confirmée : on a créé des intérieurs arides, ils n'ont aucun bassin fermé.)
5. **0 endorhéique = STRUCTUREL**, pas « pas de bassins fermés » : priority-flood remplit
   CHAQUE dépression jusqu'à son seuil de débordement → `outlet_reaches_sea` toujours vrai →
   toujours exorhéique. Les bassins fermés existent (topographiquement) mais sont sur-remplis.

## Cause racine + fix (le fix suit, pas dans cet audit)

Le drainage est **purement géométrique** : `surface_elevation = outlet sill` (remplit à ras
bord), aucun bilan hydrique. La Terre : un bassin déborde SEULEMENT si apport > évaporation ;
sinon il se stabilise SOUS le seuil (lac endorhéique à l'équilibre apport = évaporation).

→ **Fix = coupler un bilan hydrique** (le « hydroclimate layer » que le code annonce) :
par bassin, comparer l'apport (précip × aire drainée − évaporation sur le lac) au volume
jusqu'au seuil. Si apport net > 0 au seuil → déborde (exorhéique, niveau = seuil). Sinon →
endorhéique, niveau d'équilibre SOUS le seuil (apport = évaporation). Effets attendus :
l'aire de lacs chute vers ~2 %, des bassins fermés apparaissent dans les intérieurs arides
(cohérence avec le désert 30° / les steppes), les gros lacs sur-remplis se rétractent.
Tous les champs nécessaires existent déjà (précip + température → évaporation potentielle).

Audit seulement — le fix bilan-hydrique suit ce verdict.
