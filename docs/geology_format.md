# The geology export (format 0.1.0, `.ymir` container 1.1.0)

Two layers of the `.ymir` container, added by ADR Finding 152. The existing layers are unchanged.

- **`geologie_roches`**: the rock class of every cell.
- **`geologie`**: the manifest `geologie.json`, which names the favourability files written beside it.

**What the values are.** The favourabilities are **relative** (0–100 per resource): « Favorabilité dans Ymir,
rareté dans Living Landz ».
- They are not probabilities and not deposits.
- Living Landz applies a global rarity per resource and draws the deposits from the favourabilities at each reset.

## Orientation

Every PNG is **row-major with row 0 = the SOUTHERNMOST row** (the container's invariant, `y = 0` = south). Viewed in
an image viewer, the map appears upside down.

## `geologie_roches.png`

8-bit grey PNG (lossless), the grid's size (e.g. 8192 × 8192, 49 m cells). One class id per cell:

| id | key | class | hardness |
|---|---|---|---|
| 0 | `none` | water (sea, exorheic lake) | — |
| 1 | `craton` | craton (granite, gneiss) | 2 hard |
| 2 | `belt` | belt (metamorphic; fossil or active) | 2 hard |
| 3 | `basaltic_volcanic` | basaltic volcanic (shields: hotspots, rifts) | 2 hard |
| 4 | `arc_volcanic` | arc volcanic (andesite, tuff) | 1 medium |
| 5 | `rift_fill` | rift fill | 0 soft |
| 6 | `loose_deposits` | recent loose deposits (alluvium, lake margins, slope feet) | 0 soft |
| 7 | `evaporites` | evaporites (an endorheic lake's bed and its arid margin) | 0 soft |
| 8 | `basement` | undifferentiated continental basement (the default) | 2 hard |

**The hardness** is for falls at lake sills: hard keeps a fall, soft lets it retreat. The contrast is anchored on
Stock & Montgomery 1999; the values are a PROXY.

**The colours** are in `geologie.json`.

## `favorabilite_<id>.png`

- 8-bit grey PNG, values 0–100, on the **¼ grid** (e.g. 2048 × 2048, 195 m cells, ≈ 9 Living Landz hexes per cell).
- One file per resource of the rules file.
- 0 = not favourable.

## `geologie.json`

```json
{
  "format_version": "0.1.0",
  "value": "relative favourability 0-100 per resource (not a probability, not a deposit)",
  "orientation": "row-major, y = 0 = the SOUTH edge: row 0 of each PNG is the southernmost row",
  "rules": { "format": "ymir-geology-rules", "version": "0.2.0", "sha256": "…" },
  "rocks": { "file": "geologie_roches.png", "width": 8192, "height": 8192,
             "classes": [ { "id": 1, "key": "craton", "name_fr": "Craton (granite, gneiss)", "color": "#C9A55A", "hardness": 2 }, … ] },
  "favourability": { "width": 2048, "height": 2048, "factor": 4,
                     "resources": [ { "id": "gold", "name_fr": "Or", "color": "#E6C229", "file": "favorabilite_gold.png" }, … ] }
}
```

- **`rules.sha256`** is the SHA-256 of the rules file's bytes. The same world with the same rules gives the same
  favourabilities.

## The rules file

- **Location**: `geology_rules.toml` (TOML), read by Ymir. When the file is absent, Ymir uses the shipped
  `crates/ymir-core/data/geology_rules_v0.toml`.
- **The viz**:
  - the « ⛏ Géologie » menu reloads the file and redoes the zoning only;
  - a file that fails to parse is refused with its error, and the shipped rules are used.

**The semantics**:
- a rule's conditions all hold (AND), and the rule gives `chance`;
- a resource's favourability is the **max** of its rules;
- placers then travel **downstream along the rivers** from the resource's own zones (`source_min_chance`):
  - they start at `start_chance` and lose `decay_per_km` per km of flow path;
  - they are set where the slope is < `max_slope` and the drained area ≥ `river_min_area_km2` (≥ 1 km²).
- `modulate = "structural"` multiplies a rule's favourability by
  `g(d) = g_min + (1 − g_min)·min(d / d_ref, 1)` (`[modulation.structural]`):
  - d is the structural density: fracturing near the plate contacts (C-3b), the sutures and the past collisions;
  - vein deposits form districts.

**The conditions**:

| condition | meaning |
|---|---|
| `rock = ["craton", …]` | the cell's rock class key is one of these |
| `belt = "fossil" \| "active"` | a fossil belt (past collision or continental suture) or an active one (a convergent boundary today) |
| `arc_within_km` | within d km of an active-margin (subduction upper plate) cell |
| `volcano_within_km`, `volcano_active`, `volcano_setting` | within d km of an edifice's base (active; arc / hotspot / rift) |
| `wetland`, `lake`, `valley_floor` | the wetland mask; a lake; a floodplain (slope < 0.02, A ≥ 1 km²) |
| `slope_min`, `slope_max` | m/m |
| `endorheic_lake_within_km`, `coast_within_km` | distances |
| `river_min_area_km2` | drained area (real, km²) |
| `precip_min_mm`, `precip_max_mm`, `temp_min_c`, `temp_max_c` | climate |

- Unknown keys, rocks, belts or settings are refused.
- Every rule may carry a free `source` string.
