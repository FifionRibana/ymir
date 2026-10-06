# Geology v1 — the catalogue of candidate rocks and resources (Finding 150; no decision)

**What this is**: the list the author chooses from in round 2.
- Nothing here is built or decided.
- Each item says what Ymir would need to place it, and how far the existing fields go.
- The field names refer to `finding_150.md`'s inventory.

**The author's criterion** (2026-10-04): "Chutes + ressources … on va faire du zoning : ici, il y a de l'or possible ;
ici, du fer, etc." Hence **possibility zones**, not deposits.

**Sources.** The deposit contexts below are textbook-level geology, cited **from memory and to check** (none of these
is in `docs/refs`):
- **[Robb]** L. Robb 2005, *Introduction to Ore-Forming Processes*, Blackwell;
- **[Evans]** A.M. Evans 1993, *Ore Geology and Industrial Minerals*, 3rd ed., Blackwell;
- **[C&S]** D.P. Cox & D.A. Singer (eds.) 1986, *Mineral Deposit Models*, USGS Bulletin 1693;
- **[Clifford]** T.N. Clifford 1966, "Tectono-metallogenic units and metallogenic provinces of Africa", *EPSL* 1
  (the "Clifford's rule": diamondiferous kimberlites on Archean cratons).

**The existing fields this catalogue leans on** (see the inventory for whether each is filled, exported and read):

| field | grid | what it says |
|---|---|---|
| `plate_type` | 64² (6.25 km cells) | continental / oceanic |
| `cratonic_mask` | 64² | old, stable continental interior (built at init, static) |
| boundary labels | 64² | subduction upper plate (Andean margin), subducting slab, continental collision, divergent; transform via `classify_boundaries` |
| rift-soft (`age < 1` on continental crust) | 64² | young rifted crust (C-3's soft class) |
| C-3b fracture density | 64² → HD | proximity to convergent and transform contacts |
| volcanic edifices | list | centre, size, `kind` (stratocone / shield), `setting` (arc / hotspot / rift), `active`, age along a hotspot chain |
| relief, slope, flow accumulation, rivers, lakes (type: exorheic / endorheic), valley floors | HD (49 m) | the surface |
| climate (temperature, precipitation), biomes, wetland mask | HD | the surface environment |

**Missing everywhere**: a deposit / sediment-thickness field. Production erosion is detachment-limited (C-3's
module doc), so Ymir has **no sedimentary basin**.

## C1 — candidate rock classes

Hardness is relative, for the falls at lake sills: a hard sill keeps a fall, a soft one lets it retreat into a rapid.

| class | geological context | relative hardness | placed by (existing / missing) | artefact risks |
|---|---|---|---|---|
| **Ancient basement / craton** (granite, gneiss, greenstone belts) | stable Precambrian continental interiors | **hard** | existing: `cratonic_mask`, `plate_type` | the 6.25 km cells give blocky contacts. The mask is static from init, and an FBM-refined version was already rejected by C-3 as noise-made. |
| **Fold-and-thrust / orogenic belt** (metamorphic rocks and deformed sediments) | collision and active-margin belts | hard to medium | existing: collision labels, subduction upper plate, C-3b density. Missing: the deformation HISTORY (C1 records only current boundaries, so a past belt reads as craton: the fracture module's stated limit) | belts drawn along today's boundaries only; a uniform width |
| **Sedimentary basin** (sandstone, limestone, shale; also coal, evaporites) | subsiding basins: forelands, intracratonic sags, rift fills, passive margins | soft to medium | **missing**: no deposit field. Proxies: the rift-soft class (rift fill); continental lowland far from contacts and off the craton; the foreland beside a collision belt | proxies from elevation become circular (low ground called "basin" because it is low) |
| **Volcanic** (basalt shields; andesite stratocones; tuffs) | arcs, hotspots, rifts | basalt hard; andesite medium; tuff soft | existing: edifices with setting and kind; C-3's volcaniclastic discs | **perfect circles** (C-3's discs); no lava flows beyond the edifice |
| **Recent unconsolidated deposits** (alluvium, lake beds, slope deposits / colluvium) | valley floors, floodplains, lake margins and drained lake floors, the foot of steep slopes | **soft** (no fall survives) | existing, at HD: flow accumulation, the valley construction's floor width, `lake_map`, slope, the wetland mask | none of scale (49 m cells); only thresholds to declare |
| *(proposed addition)* **Evaporites / salt flats** | arid endorheic basins, playas | soft | existing: endorheic lakes, climate (aridity) | — |
| *(proposed addition)* **Oceanic crust on land / ophiolites** | slivers obducted at sutures | hard | missing: accretion sutures (C-3b's Phase B is `None`) | — |

**The C-3 lesson, to keep**: 6 km grids and perfect circles read as artefacts. Three ways to draw a contact (C3
below):
- the bilinear field's iso-line;
- a contact snapped to the terrain (a valley or a divide);
- a declared low-amplitude warp. C-3 rejected noise-made boundaries, so this one is the last resort.

## C2 — candidate resources

**Grade**:
- **A** = zonable with existing fields;
- **B** = existing fields plus a field derived from them (e.g. a downstream propagation);
- **C** = needs a missing field.

| resource | real deposit context (to check) | Ymir fields for a possibility zone | grade |
|---|---|---|---|
| **gold** | orogenic lode gold in metamorphic belts, including the greenstone belts of old cratons [Robb]; epithermal gold in volcanic arcs [C&S]; **placers** in alluvium downstream of a source [Evans] | craton, collision, subduction upper plate, arc edifices; placers = downstream of a lode zone along rivers, on low-gradient reaches | **A** (lode), **B** (placer) |
| **iron** | banded iron formations of Precambrian cratons [Robb]; **bog iron** in wetlands and peatlands [Evans] | `cratonic_mask` (BIF possibility); the wetland mask, cool humid biomes, flat lowlands (bog iron) | **A** |
| **copper** | porphyry copper above subduction zones (C&S model 17); sediment-hosted stratiform copper in rift basins [Robb]; VMS at spreading ridges and arcs | subduction upper plate, arc edifices (porphyry); rift-soft (stratiform) | **A** |
| **tin** (and tungsten) | granites of collision belts (S-type, greisens) [Robb]; cassiterite placers | collision labels; placers downstream | **A** (belt), **B** (placer) |
| **silver, lead, zinc** | veins in arcs and orogens [C&S]; SEDEX in rift or sedimentary basins; **Mississippi-Valley type in carbonate platforms** [Robb] | veins: subduction upper plate, collision; SEDEX: rift-soft. MVT: carbonate basins, missing | **A** (veins), **C** (MVT) |
| **coal** | peat swamps in subsiding sedimentary basins (forelands, intracratonic) [Evans] | **missing**: no basin. A foreland proxy (beside a collision belt, continental lowland) is possible but not causal | **C** |
| **salt** | evaporites in arid endorheic basins and playas; rift basins with marine incursions; coastal sabkhas [Evans] | endorheic lakes plus aridity (Ymir already has them); arid coastal flats | **A** |
| **building stone** | granite and gneiss (craton), basalt (volcanic), limestone and sandstone (basins) [Evans] | craton, edifices. Limestone and sandstone need the missing basin | **A** (granite, basalt), **C** (limestone, sandstone) |
| **clay** | floodplains, lake beds; kaolin from deep weathering of granite in humid tropics [Evans] | flow accumulation, valley floors, `lake_map`; craton plus a humid tropical climate (kaolin) | **A** |
| **sulphur, obsidian** | sulphur at active volcanoes (fumaroles, solfataras); obsidian from rhyolitic volcanism in arcs and continental rifts [Evans] | edifices: `active` (sulphur); `setting` arc or rift with `kind` stratocone (obsidian) | **A** |
| **gemstones** | diamonds in kimberlites on Archean cratons [Clifford]; corundum and other gems in metamorphic belts (marbles, collision); pegmatites; gem gravels (placers) [Evans] | craton (kimberlite possibility), collision (metamorphic gems), placers downstream. Pegmatites need granite bodies (missing) | **A** / **B** (partial) |
| *(proposed)* **bauxite, nickel laterite** | laterites from deep tropical weathering; nickel on ultramafic rock [Robb] | climate (humid tropical, long stable surfaces: low relief, off the active belts) | **A** (bauxite), **C** (nickel: no ultramafic rock) |
| *(proposed)* **peat** | cool, waterlogged lowlands [Evans] | wetland mask, cool humid biomes, flat ground | **A** |
| *(proposed)* **sand and gravel** | river channels and terraces, beaches [Evans] | rivers (bed width, discharge), coastline | **A** |
| *(proposed)* **hot springs / geothermal** | active volcanism, rifts | active edifices, rift-soft | **A** |

**The count over the reviewer's 11 starting resources** (a split resource counts by its main context):
- **A or B: 9** (gold, iron, copper, tin, salt, building stone, clay, sulphur / obsidian, gemstones partly);
- **C: 2** (coal; silver-lead-zinc's carbonate deposits).
- Building stone and Ag-Pb-Zn are mixed: their basin-hosted kinds are C.

### "Zonable" is not "present": the témoin (F150-I, measured)

The témoin's settled tectonic state has:
- **2 plates left** (8 at init, merged by accretion);
- **0 collision cells**, **7 subduction upper-plate cells**, 142 divergent cells;
- 57 craton cells (6.9 % of the continental cells);
- 15 edifices.

On this world, the zones that follow **collision** (tin, orogenic gold, metamorphic gems) would be **empty**, and those
that follow the **active margin** (porphyry copper, arc veins) would be a few coastal cells.
- The craton, the rifts, the edifices and every surface-driven resource (placers, clay, salt, bog iron, peat, sand
  and gravel) would be populated.
- **Whether a world has belts at all depends on the tectonic run.** On worlds with fewer merges the balance would
  differ (not measured here).

## C3 — the zoning's form (options, no decision)

1. **A possibility per cell, aggregated per hex by Living Landz.**
   - Each resource is a raster (u8 classes or a 0–255 "possibility") on the 49 m grid, exported like the biomes.
   - Living Landz takes the max or the share per hex.
   - It is simple, it follows the terrain, and the contact problem shows only where the source field is coarse.
2. **Zone polygons** (GeoJSON, like the coastline): fewer bytes, explicit names ("ceinture aurifère du nord").
   - But a polygon's edge is a decision; it inherits the 6.25 km blockiness of the coarse fields unless they are
     smoothed first.
3. **A rock-class raster plus rules.** Ymir exports the rock classes (C1), and the resource possibilities are rules
   on them.
   - It is the closest to the author's "la roche en découle", and it keeps the resources editable without
     regenerating.
   - The cost: the rules live in two places.

**Avoiding artefacts at the contacts, any option:**
- never draw a contact on the 64² grid directly (bilinear, then an iso-line);
- prefer a contact the terrain explains (a divide, a valley floor, a lake shore);
- if a warp is used, declare it as non-causal;
- a placer zone follows the rivers downstream of its source, so it is never a disc.
