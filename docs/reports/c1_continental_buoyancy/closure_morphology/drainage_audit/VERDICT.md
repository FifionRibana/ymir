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

## FIX appliqué — bilan hydrique (lacs endorhéiques)

`c1_drainage` prend désormais un `Option<&DrainageClimate>` (précip + température) — le
« hydroclimate layer » que le placeholder attendait ; `None` = géométrie pure (byte-identique).
`cached_c1_drainage` calcule le climat en interne et le fold dans la clé (nouvelle dépendance :
changer le climat réinvalide le drainage). Activé par défaut via `c165_drainage` (45°).

**Modèle** : runoff par cellule = `max(0, précip − PE)·aire` (PE = évaporation potentielle
`61·e_sat(T)`, ancrée : ~850 mm à 12 °C, ~2200 à 27 °C) ; apport = runoff accumulé en aval
jusqu'au lac ; un bassin déborde (exorhéique) si `aire_équilibre = apport/PE_lac ≥ aire_seuil`,
sinon ENDORHÉIQUE au niveau où l'aire couvre `aire_équilibre` (lu sur l'hypsométrie du lac,
direct — pas d'itération). Les cellules au-dessus de l'équilibre sont drainées.

**Mesuré (6 seeds), avant → après :**

| | géométrique | bilan hydrique |
|---|---|---|
| aire lacs / terre (agg) | 18.0 % | **0.7 %** |
| exo / endo (agg) | 1026 / 0 | **258 / 601** |
| par seed | 8-27 % | 0.5-1.0 % |

Le sur-remplissage 9× est éliminé ; **601 endorhéiques** apparaissent en zone aride/semi-aride
(85 % de l'aire de lacs), exorhéiques préservés (258). 0.7 % est **sous** le ~2 % terrestre
mais **cohérent** : le monde à 45° est semi-aride (précip ~450 < PE ~854 → terre en déficit
hydrique → peu de lacs, surtout fermés, comme les steppes) — un résultat ANCRÉ, pas calé sur
2 %. Re-validé sur la grille d'intégration (les gros lacs sur-remplis ont disparu). Lib verte ;
None byte-identique ; clé cache inclut le climat. Gated.

## SUITE tracké — les rivières (discharge depuis le runoff)

Les rivières utilisent encore l'accumulation géométrique (compte de cellules), pas le runoff :
une « rivière » dans un désert (runoff ~0) reste navigable à tort, et les rivières en aval d'un
bassin endorhéique persistent (fantômes — l'eau s'est évaporée dans le lac). FIX suivant
(parallèle, le runoff_accum existe déjà) : navigabilité depuis le runoff (discharge réel) +
reset du runoff aux lacs endorhéiques → rivières sèches dans les déserts, mortes aux bassins
fermés. Différé (le défaut majeur — 18 % de lacs — est corrigé ; les rivières fantômes sont
secondaires).
