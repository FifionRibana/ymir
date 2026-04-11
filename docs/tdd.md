# Ymir — World Generator
## Technical Design Document (TDD)

**Version**: 0.1 — Initial Draft
**Date**: April 2026
**Project**: Ymir — Physically-grounded continent generator
**Language**: Rust

---

## 1. Purpose & Scope

### 1.1 What this project is

A standalone Rust tool that generates continents for games and worldbuilding. It produces structured terrain data — heightmaps, drainage networks, biome maps, geological maps — through physically-grounded simulation rather than artistic presets. The tool runs offline; a developer iterates on parameters until the continent is satisfactory, then exports the data for consumption by a game engine or other application.

### 1.2 What this project is NOT

This is not a real-time terrain engine, not a game component, and not a map editor. It does not handle political structures, settlements, roads, or any player-created content. Those are the responsibility of the consuming application. This tool produces the natural world that exists before any human activity shapes it.

### 1.3 Why a custom generator

The current pipeline uses Azgaar's Fantasy Map Generator to produce PNG images (heightmap, binary map, biome map, lake map) that are consumed downstream. This approach has fundamental limitations that have consumed significant development time in workarounds.

The source images are raster data at low resolution (~1920×1006 pixels). Upscaling ×100 creates aliasing artifacts at every boundary: pixelated coastlines requiring blur and noise to mask, biome boundaries with anti-aliased RGB values that don't match any known biome, height quantization in u8 (256 levels for 2500m of range) causing visible terracing, and no smooth coastal gradient requiring a synthetic coastal slope to be injected. Each correction introduces new artifacts requiring further corrections — a spiraling complexity.

Beyond artifacts, the data itself is limited. Azgaar provides no bathymetry (underwater depth is faked from coastal SDF), no drainage network at meso scale (only major rivers as polylines), no geological data for mineral placement, and no wind/precipitation data for climate modeling. The terrain lacks meso-scale coherence: FBM noise adds hills but they are randomly oriented, creating geologically impossible configurations like ridges perpendicular to drainage or wetlands on cliffs.

A physically-grounded generator produces terrain where the meso structures (valleys, ridges, drainage basins) are causally coherent. Mountains form because plates converge. Valleys form because water erodes. Biomes form because rain falls on windward slopes and doesn't fall on leeward slopes. Minerals occur because specific tectonic processes concentrate specific elements. Everything has a cause, and the causes are computable.

### 1.4 Integration with downstream consumers

The generator outputs a set of structured data files (see §9). These can be consumed by any game engine, simulation tool, or rendering pipeline that needs continent-scale terrain data. The output format is designed to be self-describing (JSON metadata + raw binary rasters) and independent of any specific game engine or framework.

---

## 2. Gameplay Requirements

The generated continent must satisfy specific gameplay needs for a medieval sandbox game (the primary use case driving this generator's design).

### 2.1 Scale and traversal

A medieval pedestrian covers approximately 30-40 km per day. Exploration of the continent should represent a significant investment — weeks of in-game time to map the coastline, months to explore the interior. This implies a continent of 280-350 km across. At 35-40m per simulation pixel on the erosion grid, this maps to a grid of approximately 8192² pixels.

### 2.2 Settlement placement and valley sizing

The urban classification from the GDD defines the settlement hierarchy:

| Classification | Population | Approximate footprint | Minimum valley width |
|---------------|------------|----------------------|---------------------|
| Hameau | up to 200 | ~12 hex residential + fields, 150-200m core | 150m (small valleys OK) |
| Village | up to 1,500 | ~20 hex residential, extends along river/valley | 200-300m |
| Bourg | up to 5,000 | ~25-30 hex residential, compact | 400-600m |
| Ville | up to 20,000 | ~50-70 hex residential | 800m-1.2km |
| Cité | 100,000+ | ~190 hex residential, multi-Voronoï-cell | 1.5-2km+ (large valley or coastal plain) |

A key gameplay insight: hameaux can squeeze into almost any valley, including narrow mountain ones. This means hameaux can access the most interesting mountain resources (rare ores in subduction zones, gems in rift zones) but at the cost of limited growth potential — the valley can't support a larger settlement. A player who wants to exploit mountain resources must choose: stay small in a rich narrow valley, or settle in a wide valley with growth potential but fewer exotic resources. This creates a natural strategic tension that the terrain generator must support.

The erosion simulation must produce valleys at all these scales with a natural hierarchy: ravines too narrow for any settlement (50-100m, adds visual texture), hameau valleys (150-250m, common in foothills and mountain flanks), village/bourg valleys (300-600m, along secondary rivers), ville valleys (800m-1.2km, along major rivers or in coastal plains), and cité-scale floodplains (2km+, at river confluences or in large coastal embayments).

### 2.3 Geographic diversity

A single continent must contain multiple distinct geographic regions that feel different to explore and offer different gameplay: coastal lowlands with ports and fishing, agricultural plains, forested hills with timber and game, mountain ranges with mines and strategic passes, river valleys as natural corridors and trade routes.

The tectonics simulation should produce 2-4 distinct mountain systems (collision zones, volcanic arcs, rift shoulders) rather than a single uniform range. The climate model should create rain shadows and continental gradients that differentiate windward forests from leeward grasslands.

### 2.4 Continent specialization

Each continent corresponds to a latitude band and has a dominant climate character: temperate, tropical, cold, arid. This simplifies the climate model (temperature is primarily a function of altitude rather than latitude) while creating distinct gameplay per continent. New continents appear narratively when the Mist concentrates and recedes elsewhere.

### 2.5 River networks and commerce

Rivers are both geographic features and gameplay infrastructure: barriers to cross (requiring bridges or fords), corridors to follow (valleys are the flattest paths through mountains), and trade routes (water transport is cheaper than overland). The generated drainage network must include major navigable rivers (flux above a threshold), secondary rivers suitable for small boats, and tributaries that structure the landscape without being individually navigable.

### 2.6 Mineral resources

The game features mining as a core economic activity. Different ores and stone types must be placed in geologically plausible locations. The tectonic simulation provides the data for this: subduction zones concentrate copper, gold, and silver; old continental shields host iron and nickel; rift zones produce gemstones; sedimentary basins contain coal, limestone, salt, and clay. The generator must output a geological classification per region that downstream systems use for resource spawning.

---

## 3. Pipeline Overview

The generation proceeds in six sequential phases, each building on the output of the previous one. The pipeline is designed so that each phase can be run independently during development, and the first phases (which are fastest) can be iterated rapidly while the later phases (which are most expensive) run less frequently.

```
Phase 1: Tectonics (128²–512²)
    Crustal thickness field from plate simulation
        ↓
Phase 2: Isostasy (same grid)
    Altitude = f(crustal thickness)
        ↓
Phase 3: Upscale + Anisotropic FBM (4096²–8192²)
    High-resolution heightmap with directionally coherent detail
        ↓
Phase 4: Hydraulic Erosion (same grid)
    Carved valleys, ridges, drainage network, sediment deposits
        ↓
Phase 5: Climate & Biomes (same grid)
    Temperature, precipitation, wind effect → Whittaker biome classification
        ↓
Phase 6: Export
    Structured data files for downstream consumption
```

Total estimated generation time on a modern multi-core CPU: 3-8 minutes for the full pipeline at 8192² erosion resolution. Development iterations on 1024² take under 30 seconds.

---

## 4. Phase 1 — Tectonics

### 4.1 Model: Thin Viscous Sheet

The lithosphere is modeled as a 2D thin viscous sheet (England & McKenzie, 1982; Houseman & England, 1986). Instead of simulating 3D mantle convection, this model solves a 2D flow equation for crustal thickness, which captures the essential deformation behavior at a fraction of the computational cost.

The primary field is **crustal thickness S(x,y,t)**, measured in dimensionless units (1.0 = reference continental thickness ≈ 35 km, 0.2 = oceanic thickness ≈ 7 km). The evolution equation is:

```
∂S/∂t = -∇·(S·v) + source_terms
```

where **v(x,y)** is the velocity field, solved at each timestep from force balance.

### 4.2 Plate setup

The domain is a 2D grid with periodic boundary conditions (torus topology), representing a section of planetary surface large enough to contain one continent and surrounding ocean. Grid size: 128² for development, 256²–512² for production.

Plates are initialized as Voronoï regions (5-15 plates). Each plate is assigned a type (continental: S=1.0, or oceanic: S=0.2) and a velocity vector. The continental/oceanic ratio should be approximately 30%/70%. Velocity magnitudes are in the range 1-5 cm/year (dimensionless units scaled to the grid).

The plate configuration is the primary creative input. Different configurations produce different continent shapes: two converging continental plates create a Himalaya-type collision; an oceanic plate subducting under a continent creates an Andes-type arc; a diverging continental plate creates an East Africa-type rift. The developer adjusts plate count, types, and velocities until the resulting geography is interesting.

### 4.3 Velocity field

At each timestep, the velocity field is computed from the balance of three forces:

**Plate traction**: the imposed velocity from plate motion. Each grid cell belongs to a plate and inherits its velocity. At plate boundaries, velocities are interpolated or discontinuous (depending on the boundary type).

**Gravitational pressure gradient**: thick crust "wants" to spread laterally (gravitational potential energy). This force is proportional to the gradient of S²: `F_grav = -∇(S²)`. It prevents unlimited crustal buildup at convergent boundaries and drives collapse of overthickened orogens.

**Viscous resistance**: the lithosphere resists deformation. In the thin sheet approximation, this is modeled as a 2D viscous stress proportional to the strain rate (velocity gradient). The effective viscosity can be power-law (non-Newtonian: n=3 gives more realistic localized deformation than n=1 linear viscosity).

The velocity is solved by discretizing the force balance on the grid (finite differences) and solving the resulting linear system iteratively (Gauss-Seidel or conjugate gradient). This is the most computationally expensive step per timestep, but on a 512² grid it takes milliseconds.

### 4.4 Boundary processes

**Convergent (continental-continental)**: crustal thickness increases (S grows). The collision zone thickens to S=1.5-2.0, creating a mountain belt. The gravitational pressure gradient limits the maximum thickness.

**Convergent (oceanic-continental / subduction)**: the oceanic plate's thickness is removed (subducted into the mantle). A volcanic source term is added ~100-200 km behind the subduction front: `source_volcanic = rate × exp(-distance/decay)`. This creates an arc of thickened crust parallel to the coast.

**Divergent**: the crust thins (S decreases). If S drops below ~0.3, new oceanic crust forms (S resets to 0.2). This creates rifts and new ocean basins.

**Transform**: plates slide past each other. Minimal crustal thickness change, but the velocity discontinuity creates a narrow zone of deformation.

### 4.5 Volcanism

Volcanic activity is modeled as a source term in the thickness equation, not as individual cone placement. At subduction zones, the source creates a broad belt of thickened crust (volcanic arc). At rifts, the source creates localized thickening where the crust is thinnest (shield volcanoes on rift floors). At hotspots (optional), a fixed point source creates a chain of thickened spots as the plate moves over it (Hawaiian-type chains).

### 4.6 Simulation duration and timestep

The simulation runs for N timesteps representing tens of millions of years of geological time. The number of timesteps is a parameter (typically 200-500). More timesteps produce more evolved, mature landscapes; fewer timesteps produce younger, more active geographies. The CFL condition constrains the timestep size relative to grid spacing and maximum velocity.

### 4.7 Output

Phase 1 produces a 2D field of crustal thickness S(x,y) on the tectonic grid (128²–512²). It also produces the final velocity field (useful for determining tectonic stress orientation for mineral placement) and a classification of each grid cell's tectonic history (stable continental, collision zone, subduction zone, rift zone, oceanic) for geological resource mapping.

### 4.8 References

England, P., & McKenzie, D. (1982). "A thin viscous sheet model for continental deformation." Geophysical Journal of the Royal Astronomical Society, 70(2), 295-321.

Houseman, G., & England, P. (1986). "Finite strain calculations of continental deformation: 1. Method and general results for convergent zones." Journal of Geophysical Research, 91(B3), 3651-3663.

Flesch, L. M., et al. (2001). "Dynamics of the India-Eurasia collision zone." Journal of Geophysical Research, 106(B8), 16435-16460.

---

## 5. Phase 2 — Isostasy

### 5.1 Altitude from crustal thickness

The conversion from crustal thickness to surface altitude follows Airy isostasy. The crust "floats" on the denser mantle; thicker crust protrudes higher above the reference level (like a thicker iceberg).

```
altitude(x,y) = (S(x,y) - S_ref) × (ρ_mantle - ρ_crust) / ρ_mantle × scale_factor
```

where S_ref is the reference oceanic thickness (~0.2), ρ_mantle ≈ 3.3, ρ_crust ≈ 2.7. The scale_factor maps dimensionless units to meters, calibrated so that a collision zone (S=1.8) produces mountains of ~3000-4000m and normal continental crust (S=1.0) sits at ~200-400m.

Sea level is defined such that the oceanic crust (S=0.2) is submerged. Altitude values below sea level represent ocean floor (bathymetry).

### 5.2 Smoothing

The raw isostatic altitude has the resolution of the tectonic grid (128²–512²). Before upscaling, a gentle Gaussian blur (sigma=2-3 grid cells) smooths the sharpest transitions while preserving the large-scale structure. This prevents aliasing when the next phase upscales to the erosion grid.

### 5.3 Output

A 2D altitude field at tectonic grid resolution, in meters. Positive = above sea level, negative = below. This is the "initial" landscape before erosion — geologically young, with broad plateaus, sharp collision fronts, and no drainage carving.

---

## 6. Phase 3 — Upscale & Anisotropic FBM

### 6.1 Interpolation

The tectonic altitude grid (128²–512²) is upscaled to the erosion grid (4096²–8192²) via bicubic interpolation. This produces a smooth surface that preserves the large-scale features without introducing artifacts.

### 6.2 Anisotropic noise

Standard FBM noise is isotropic — it creates bumps in all directions equally. Real mountain terrain is strongly anisotropic: ridges run parallel to the collision front, valleys run perpendicular (draining toward the lowlands). The anisotropic FBM uses the gradient of the tectonic thickness field to orient the noise.

At each pixel of the high-resolution grid, the gradient of S (interpolated from the tectonic grid) gives the direction of steepest ascent. Ridged noise is generated along the perpendicular direction (parallel to the mountain chain), creating ridge-and-valley texture aligned with the geological structure.

```
ridge_direction(x,y) = rotate_90(normalize(∇S(x,y)))
noise_value = ridged_fbm(position projected onto ridge_direction, octaves=4-6)
```

The amplitude of the noise scales with the local gradient magnitude: flat plains get very little noise (they're already flat), steep mountain fronts get strong ridging. This naturally creates the transition from smooth plains to corrugated mountains.

### 6.3 Noise parameters

Meso layer (hills, ridges): frequency scaled to produce features of 500m-2km wavelength at the erosion grid resolution. Amplitude 50-200m, modulated by tectonic gradient. 4-5 octaves of ridged noise.

Micro layer (surface texture): frequency scaled for 100-300m features. Amplitude 10-30m. 2-3 octaves of standard FBM. This layer is isotropic — surface roughness doesn't need directional coherence.

### 6.4 Output

A high-resolution heightmap (4096²–8192²) at 35-40m per pixel, combining the tectonic macro structure with directionally coherent meso detail. This is the input to hydraulic erosion.

---

## 7. Phase 4 — Hydraulic Erosion

### 7.1 Algorithm

The erosion simulation follows the particle-based approach of Beyer (2015), adapted for large grids. Individual water droplets are simulated: each drops at a random position, follows the terrain gradient downhill, erodes material from steep slopes, carries sediment, and deposits when velocity decreases.

For each droplet:
1. Place at random position on the terrain
2. Calculate gradient (bilinear interpolation of neighboring heights)
3. Update velocity: accelerate downhill, decelerate by friction
4. Move in the direction of the gradient (with inertia for smooth paths)
5. Erode: remove material proportional to (carrying_capacity - current_sediment), capped by erosion_rate
6. Deposit: if current_sediment > carrying_capacity, deposit the excess
7. Evaporate: reduce water volume slightly each step
8. Repeat until the droplet evaporates or exits the grid

Carrying capacity is proportional to velocity × water_volume × slope. This naturally creates deeper erosion on steep slopes and deposition on flat areas.

### 7.2 Parameters

The erosion parameters control the character of the resulting landscape:

Erosion rate: 0.3-0.5 (higher = more aggressive carving, sharper valleys). Deposition rate: 0.3-0.5. Evaporation rate: 0.01-0.02 per step. Inertia: 0.05-0.1 (higher = smoother, more meandering rivers). Gravity: 4.0-10.0 (higher = faster water, deeper cuts). Min slope: 0.01 (prevents water from pooling on perfectly flat terrain). Max lifetime: 100-200 steps per droplet.

Number of droplets: typically 2-5 million for a 4096² grid, 10-20 million for 8192². This is the main performance knob — more droplets = more refined terrain but longer computation.

### 7.3 Parallelization

Individual droplets are independent except for shared heightmap access. The grid can be divided into tiles, with droplets processed in parallel within tiles. Boundary regions require synchronization. Rayon's parallel iterators work well for this: process droplets in batches of 10,000-50,000, synchronize the heightmap after each batch.

### 7.4 Flow accumulation

As a post-process after erosion, compute the flow accumulation: for each cell, count how many cells drain through it. This is computed by sorting all cells by decreasing altitude, then for each cell adding its accumulated flow to its downhill neighbor.

Flow accumulation directly gives river size: cells with flow > threshold_1 are navigable rivers, cells with flow > threshold_2 are major navigable waterways. The thresholds are calibrated to match the desired river density.

### 7.5 Sediment map

The erosion process produces a secondary output: the sediment thickness at each cell (material deposited by water). Sediment-rich areas correspond to alluvial plains and river deltas — the best agricultural land. This map is exported for downstream use in fertility calculations.

### 7.6 Lake formation

After hydraulic erosion, closed depressions in the terrain naturally fill with water to form lakes. Detection uses a priority flood fill algorithm: starting from the ocean (the global minimum), propagate inward, and any cell lower than all its outflow paths is filled to the level of its lowest outlet. The filled volume is the lake.

Each lake is characterized by its surface altitude, maximum depth, surface area, and outlet river (if any — endorheic lakes have no outlet). Lake type is inferred from context:

- **Glacial lake**: in a U-shaped valley at high altitude (if glacial erosion is implemented) or in an over-deepened erosion basin in mountains.
- **Rift lake**: in a tectonic graben (low crustal thickness relative to surroundings), typically long and narrow. Deep.
- **Volcanic crater lake**: within a volcanic caldera (identified from the tectonic volcanic source term). Can be hostile — acidic if the volcano is active (recent volcanic source term), neutral if dormant.
- **Oxbow lake**: in a floodplain near a meandering river (low altitude, near high-flux river, in sediment-rich zone). Small and shallow.
- **Frozen lake**: any lake where T_mean at lake altitude < 0°C. Permanently frozen or seasonally frozen depending on T_amplitude.

The lake mask and depth map are exported alongside the heightmap.

### 7.7 Coastal erosion (optional post-process)

After hydraulic erosion and lake formation, an optional coastal erosion phase sculpts the shoreline into organic shapes. Wave erosion acts on the coastline from the ocean side, controlled by exposure to prevailing wind and rock resistance.

For each coastal cell, wave energy is computed as a function of fetch (distance of open ocean in the wind direction) and exposure angle. High-energy coasts (facing the wind, long fetch) erode faster. Rock resistance is derived from the geological classification: sedimentary coasts (basins) are soft and erode into smooth bays; crystalline coasts (shields, collision zones) are hard and form resistant headlands; volcanic coasts are intermediate.

The algorithm iterates: erode exposed soft-rock coastline cells, check if newly exposed cells are also vulnerable, repeat. This naturally creates the alternating bay-headland pattern seen on real coastlines.

### 7.8 Output

The eroded heightmap (same resolution as input), the flow accumulation map (for river extraction), the sediment deposit map (for agricultural fertility), the drainage direction map (for each cell, which neighbor it drains toward), the lake mask with depth and type classification, and optionally the coastally-eroded shoreline.

### 7.7 References

Beyer, H.T. (2015). "Implementation of a Simple Erosion Model." Technische Universität München.

Cordonnier, G., et al. (2016). "Large Scale Terrain Generation from Tectonic Uplift and Fluvial Erosion." Computer Graphics Forum, 35(2), 165-175.

---

## 8. Phase 5 — Climate & Biomes

### 8.1 Temperature

Temperature at each grid cell is computed from three components.

Base temperature is a continent-level parameter representing the mean annual temperature at sea level at the continent's central latitude. Values: tropical ~25°C, temperate ~12-15°C, cold ~3-5°C, polar ~-5°C.

Altitude lapse: temperature decreases with altitude at the adiabatic lapse rate of approximately 6.5°C per 1000m. `T_altitude = -6.5 × altitude_km`.

Continentality: coastal areas have buffered temperatures (mild winters, cool summers), inland areas have extreme temperatures (cold winters, hot summers). The amplitude of seasonal variation scales with distance from the ocean: coastal ~5-8°C amplitude, deep inland ~15-20°C. Distance to ocean is computed as a distance transform from the coastline (cells with altitude ≤ 0).

```
T_mean(x,y) = T_base - 6.5 × altitude(x,y)/1000
T_amplitude(x,y) = 5 + 12 × (distance_to_ocean(x,y) / max_distance)
T(x,y,day) = T_mean + T_amplitude × sin(2π × (day - offset) / 365)
```

An optional north-south gradient (0.5-1.0°C per 100km) adds a directional temperature trend across the continent, creating slightly warmer southern and cooler northern regions (or vice versa).

### 8.2 Precipitation

Precipitation is the most important factor for biome diversity. The model has three components.

**Oceanic moisture**: precipitation decreases exponentially with distance from the coast. `P_base = P_coastal × exp(-distance / moisture_decay)` where P_coastal ≈ 800-1200mm/year and moisture_decay ≈ 300-500 km.

**Orographic effect**: this is the dominant source of precipitation variability. Moist air carried by the prevailing wind rises when it encounters terrain, cools, and drops rain. The windward side of mountains receives enhanced rainfall; the leeward side is in rain shadow.

The prevailing wind is defined as a single direction vector per continent (e.g., from the west-southwest for a temperate continent, matching Earth's mid-latitude westerlies). For each grid cell, the "windward exposure" is computed as the dot product of the terrain gradient with the wind direction. Positive exposure (slope facing the wind) enhances precipitation; negative (slope facing away) reduces it.

```
exposure(x,y) = dot(∇altitude(x,y), wind_direction)
P_orographic = P_base × (1 + orographic_factor × exposure)
```

The orographic_factor is typically 2-4, meaning a windward slope can receive 3-5× the base precipitation while the leeward slope receives 0.3-0.5× the base.

**Altitude reduction**: very high altitudes receive less precipitation (the air is already dry). `P_altitude_factor = exp(-altitude / 3000)`.

Final precipitation: `P(x,y) = P_base × P_altitude_factor × (1 + orographic_factor × clamp(exposure, -0.8, 1.0))`

### 8.3 Biome classification (Whittaker diagram)

Each grid cell has a mean annual temperature and annual precipitation. These two values are mapped to a biome via a simplified Whittaker diagram:

| Biome | T_mean (°C) | Precipitation (mm/year) |
|-------|-------------|------------------------|
| Ice/Glacier | < -5 | any |
| Tundra | -5 to 0 | < 500 |
| Taiga | 0 to 5 | > 300 |
| Cold Desert | -5 to 5 | < 300 |
| Temperate Rain Forest | 5 to 15 | > 1500 |
| Temperate Deciduous Forest | 5 to 15 | 600-1500 |
| Grassland | 5 to 15 | 300-600 |
| Desert | > 15 | < 300 |
| Savanna | > 15 | 300-800 |
| Tropical Seasonal Forest | > 20 | 800-2000 |
| Tropical Rain Forest | > 20 | > 2000 |
| Wetland | any (special) | (special: flat + high water table) |

Wetlands are classified separately: cells that are flat (slope < threshold), at low altitude (< 50m above nearest river), and have moderate-to-high precipitation are classified as wetland regardless of the Whittaker lookup. This prevents the issue from the current pipeline where wetlands appear on cliffs.

The exact temperature and precipitation thresholds are parameters to be calibrated during development. The Whittaker diagram is stored as a lookup texture or a decision tree, easily adjustable.

### 8.4 Seasonal data for gameplay

Beyond biome classification, the climate model exports seasonal data for gameplay mechanics in downstream applications.

Growing season length: the number of days per year where T > 5°C. This directly impacts agriculture (the GDD's "warm time" / "cold time" tempo).

Frost risk: the probability of sub-zero temperatures in each month. Determines crop vulnerability and building heating requirements.

Precipitation seasonality: whether rain is evenly distributed (oceanic climate), summer-dry (Mediterranean), or winter-dry (continental). Affects irrigation needs and flood risk.

### 8.5 Output

Per-cell: biome ID, mean temperature, temperature amplitude, annual precipitation, growing season length, precipitation seasonality index. Per-continent: prevailing wind direction, base temperature, climate type.

---

## 9. Phase 6 — Export

### 9.1 Geological map for mineral placement

The tectonic simulation classifies each region of the continent by its geological history. This classification maps directly to mineral resource families.

| Tectonic Context | Geological Character | Associated Resources |
|------------------|---------------------|---------------------|
| Stable continental interior | Old shield / craton | Iron, nickel, gold (alluvial) |
| Collision zone (active) | Folded sedimentary + metamorphic | Marble, slate, precious gems |
| Collision zone (mature) | Granite intrusions | Tin, tungsten, rare earths |
| Subduction zone / volcanic arc | Volcanic + hydrothermal | Copper, gold, silver, sulfur |
| Rift zone | Basalt + alkaline intrusions | Gemstones, soda minerals |
| Sedimentary basin | Thick sediment layers | Coal, limestone, salt, clay, sandstone |
| Alluvial plain (from erosion) | River deposits | Sand, gravel, alluvial gold |
| Coastal | Marine sediments | Salt, chalk, clay |

The consuming application uses this geological map combined with noise-based placement to spawn specific resource nodes. The generator provides the probability distribution; the consumer instantiates individual deposits.

### 9.2 River network

The flow accumulation map is thresholded and vectorized into a river network graph. Each river segment has: start/end coordinates, upstream/downstream connections, average flux (used for width and navigability classification), gradient (used for flow speed), and the basin ID it belongs to.

Navigability classes based on flux thresholds: non-navigable (streams, < flux_1), small boat (flux_1 to flux_2), barge (flux_2 to flux_3), ship (> flux_3). The thresholds are calibrated to produce a plausible network density.

### 9.3 Output file format

The generator produces a directory of files for each continent:

```
continent_name/
├── metadata.json          # Dimensions, scale, climate params, generation params, seed
├── heightmap.raw          # u16 LE, eroded altitude, 8192² (or configured size)
├── bathymetry.raw         # u16 LE, underwater depth (0 = surface, 65535 = max depth)
├── geology.raw            # u8, tectonic classification per cell (tectonic grid resolution)
├── biomes.raw             # u8, biome ID per cell (erosion grid resolution)
├── temperature.raw        # i16 LE, mean annual temperature × 100 (for 0.01°C precision)
├── precipitation.raw      # u16 LE, annual precipitation in mm
├── flow_accumulation.raw  # u32 LE, drainage flux per cell
├── drainage_direction.raw # u8, direction to downhill neighbor (0-7 for 8 directions)
├── sediment.raw           # u16 LE, sediment deposit thickness
├── lakes.json             # Lake catalog: id, surface altitude, depth, area, type, outlet
├── lake_mask.raw          # u8, lake ID per cell (0 = not lake)
├── lake_depth.raw         # u16 LE, depth below lake surface per cell
├── rivers.json            # Vectorized river network graph
└── crustal_thickness.raw  # f32, tectonic thickness field (tectonic grid resolution)
```

All raster files are raw binary at their documented resolution, no headers. The metadata.json contains everything needed to interpret them: grid dimensions, meters-per-pixel, max elevation, max depth, sea level reference, climate parameters, and generation seed for reproducibility.

### 9.4 Integration with existing server pipeline

### 9.4 Integration path

The export format is designed for progressive adoption by a consuming application.

**Phase 1 (minimal)**: convert the generator's output to PNG images matching common fantasy map formats (heightmap, binary map, biome map, lake map). An existing pipeline that consumes Azgaar-style PNGs can switch to Ymir output without code changes. This allows immediate testing.

**Phase 2 (proper)**: the consuming application reads Ymir's native format directly. The richer data (drainage, geology, sediment, climate) becomes available to downstream systems. The PNG conversion step is removed.

**Phase 3 (full integration)**: downstream systems directly consume the generator's data for resource placement, agriculture fertility, climate effects, river navigation. Any enrichment pipeline in the consumer can be simplified or removed, since Ymir already produces a high-resolution eroded heightmap with coherent drainage.

---

## 10. Development Visualization

### 10.1 Requirements

During development, the generator needs real-time visualization to iterate on parameters. This is a minimal Bevy application that displays the terrain being generated and allows parameter adjustment.

### 10.2 Views

**Tectonic view**: crustal thickness as a colormap (thin/blue → thick/red), plate boundaries highlighted, velocity vectors shown as arrows. Updated each tectonic timestep.

**Altitude view**: heightmap as a hillshaded terrain with hypsometric coloring (green lowlands → brown hills → white peaks → blue ocean). This is the primary view for evaluating the overall continent shape.

**Erosion view**: animated during erosion simulation, showing valleys being carved in real-time. Water flow paths visible as blue traces.

**Climate view**: temperature and precipitation as colormaps. Biome classification as a flat-color map matching the game's biome palette.

**Geological view**: tectonic classification colored by type (shield, collision, subduction, rift, sedimentary). This is the map that determines mineral placement.

**River view**: drainage network overlaid on the terrain, colored by flux (thin blue → thick blue for increasing navigability).

### 10.3 Controls

Camera: pan and zoom on the 2D map. Parameter sliders for the current phase (plate velocities, erosion rates, climate thresholds). "Run phase" buttons to execute each phase individually. "Run all" button for the full pipeline. Export button to save the current state.

---

## 11. Grid Sizing Reference

### 11.1 Development grid (fast iteration)

| Phase | Grid | m/pixel | Coverage | Est. time |
|-------|------|---------|----------|-----------|
| Tectonics | 64² | 540 | ~35 km | ~1s |
| Isostasy | 64² | 540 | ~35 km | <1s |
| Upscale + FBM | 512² | 68 | ~35 km | ~2s |
| Erosion | 1024² | 34 | ~35 km | ~15s |
| Climate | 1024² | 34 | ~35 km | ~2s |
| **Total** | | | | **~20s** |

### 11.2 Production grid (final continent)

| Phase | Grid | m/pixel | Coverage | Est. time |
|-------|------|---------|----------|-----------|
| Tectonics | 512² | 680 | ~350 km | ~30s |
| Isostasy | 512² | 680 | ~350 km | <1s |
| Upscale + FBM | 4096² | 85 | ~350 km | ~10s |
| Erosion | 8192² | 43 | ~350 km | ~3-5 min |
| Climate | 4096² | 85 | ~350 km | ~5s |
| **Total** | | | | **~4-6 min** |

Memory peak: 8192² × 4 bytes (f32 working heightmap) + secondary buffers ≈ 400 MB.

### 11.3 Continent sizing (for reference)

At 8192² erosion grid with 40m/pixel = 328 km side. A continental landmass occupying ~40% of the grid (rest is ocean) gives approximately 43 000 km², comparable to Switzerland + Belgium combined. This provides adequate space for the gameplay requirements described in §2.

---

## 12. Design Decisions (resolved)

**Thermal erosion**: not in the initial implementation. Added as a refinement phase once hydraulic erosion is functional and validated. Thermal erosion sharpens ridges and creates talus slopes — it improves realism but is not structurally necessary for the first version.

**Lake formation**: in scope. Lakes are a gameplay resource (fresh water, fishing, potentially hostile environments). The erosion phase must detect closed depressions and fill them to form lakes. Lake classification is derived from geological context: glacial lakes (eroded valleys in mountains), rift lakes (tectonic depressions), volcanic crater lakes (subduction zone calderas — can be hostile, e.g. acidic/sulfuric like Kawah Ijen), and oxbow lakes (river meanders, post-erosion). Frozen lakes are computed from the climate model (T_mean < 0°C at lake altitude). The lake type affects gameplay resources and hazards. Implementation: priority flood fill as a post-process after hydraulic erosion, with lake type inferred from tectonic classification and altitude.

**Coastal erosion**: in scope as an optional post-processing phase after hydraulic erosion. Wave erosion creates the organic coastal features (bays, headlands, sea cliffs, arches) that are currently missing. The model is exposure-based: coastline segments facing the prevailing wind receive more wave energy and erode faster. Rock resistance (from the geological classification) modulates the rate — sedimentary coasts erode into gentle bays, granitic coasts resist and form headlands. This is simpler than the hydraulic erosion model and can be implemented as a separate, independent phase.

**Multi-continent consistency**: yes, shared global parameters. Each world has a global configuration (sea level, prevailing wind direction per latitude band, global temperature gradient) that all continents inherit. Individual continents are generated independently but within this shared framework. This ensures that a tropical continent and a temperate continent in the same world have consistent climate physics — they differ because of latitude, not because of different rules.

**Seed reproducibility**: required. Same seed + same parameters = identical output, byte-for-byte. This is essential for iterative development (tweak one parameter, keep everything else identical) and for reproducible world generation. Implementation constraints: all parallel operations must be order-independent or explicitly ordered. The erosion simulation processes droplets in deterministic batches — within each batch, droplets are sorted by initial position (derived from seed) to ensure consistent processing order regardless of thread scheduling. Floating-point operations use consistent rounding (no fast-math optimizations that reorder operations). The seed cascades through phases: each phase derives its own sub-seed from the master seed, so changing Phase 5 parameters doesn't affect Phase 1-4 output.

**Integration timeline**: Ymir will be developed and validated independently. When it produces satisfactory results, the consuming application will switch from Azgaar PNG input to Ymir's native format. During development, Ymir can export PNG images matching the current Azgaar format (§9.4 Phase 1) for testing without modifications to the consumer. No formal transition period with dual support — the switch is a single commit on the consumer side.

**Aeolian (wind) erosion**: in scope as an additional erosion mode for arid continents. Hydraulic erosion assumes significant rainfall and produces V-shaped valleys; wind erosion produces different landforms (sand dunes, deflation basins, yardangs, desert pavements). For desert-specialized continents, aeolian erosion runs after hydraulic erosion (which still operates, just weakly, since even deserts have occasional flash floods). The model is exposure-based: loose sediment (from the sediment map) is picked up by wind and deposited downwind, creating dune fields in sediment-rich areas and exposing bedrock elsewhere. Not in the first milestone, but designed into the pipeline as a pluggable erosion phase.

**Glacial erosion**: in scope as an additional erosion mode for cold continents. Glaciers carve U-shaped valleys (wide, flat-bottomed) rather than the V-shaped valleys of water erosion, and they create distinctive features: cirques (amphitheater-shaped headwalls), arêtes (sharp ridges between two glacial valleys), fjords (drowned glacial valleys at the coast), moraines (sediment ridges at glacier termini), and hanging valleys. For cold-specialized continents, glacial erosion runs on cells where T_mean < 0°C and precipitation is sufficient for ice accumulation. The glacier flows downhill (like water, but much slower and wider), carving a wider trough. Not in the first milestone, but the temperature model from Phase 5 provides the necessary input.

**River deltas**: in scope, integrated into the hydraulic erosion phase rather than as a separate system. Deltas form naturally when erosion droplets are allowed to continue past the coastline (altitude ≤ 0) for a few steps before dying. Sediment carried by the droplet deposits in the shallow coastal zone, building up a fan of new land. This is implemented as a parameter in the erosion loop (coastal_deposition_range: number of steps a droplet survives below sea level). Included from the first milestone since it affects coastline shape significantly and is trivial to implement (a small change to the droplet termination condition).

## 13. Open Questions (remaining)

**Cave systems**: out of scope. The game currently has no notion of underground exploration or sub-surface levels. The generator focuses exclusively on surface terrain. Caves and underground networks could be a future extension if the game design evolves to include them, but they are not considered in this pipeline.

---

## 13. References

England, P., & McKenzie, D. (1982). "A thin viscous sheet model for continental deformation." Geophysical Journal of the Royal Astronomical Society, 70(2), 295-321.

Houseman, G., & England, P. (1986). "Finite strain calculations of continental deformation." Journal of Geophysical Research, 91(B3), 3651-3663.

Flesch, L. M., Haines, A. J., & Holt, W. E. (2001). "Dynamics of the India-Eurasia collision zone." Journal of Geophysical Research, 106(B8), 16435-16460.

Turcotte, D. L., & Schubert, G. (2002). "Geodynamics." Cambridge University Press.

Beyer, H. T. (2015). "Implementation of a Simple Erosion Model." Technische Universität München.

Cordonnier, G., Braun, J., Cani, M.-P., Benes, B., Galin, E., Peytavie, A., & Guérin, E. (2016). "Large Scale Terrain Generation from Tectonic Uplift and Fluvial Erosion." Computer Graphics Forum, 35(2), 165-175.

Whittaker, R. H. (1975). "Communities and Ecosystems." Macmillan.