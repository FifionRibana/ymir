# M0 Issues

---

## Issue #1 — Project scaffolding and workspace setup

### Description

Set up the Ymir project structure as a Rust workspace with the core library and visualization binary. The workspace must support future growth (additional binaries for batch export, CLI tools) without restructuring.

### Tasks

- Initialize the Cargo workspace with two members: `ymir-core` (library crate for all generation logic) and `ymir-viz` (binary crate for the Bevy visualization app).
- `ymir-core` depends on: `rayon` (parallelism), `image` (PNG I/O), `rand` and `rand_chacha` (seeded deterministic RNG), `serde` + `serde_json` (configuration and metadata serialization).
- `ymir-viz` depends on: `ymir-core` and `bevy` (visualization).
- Create the module structure in `ymir-core`: `erosion/` (hydraulic, later thermal/aeolian/glacial), `tectonics/` (plate simulation, isostasy), `climate/` (temperature, precipitation, biomes), `terrain/` (heightmap grid, noise functions, upscale), `export/` (output file writers), `config.rs` (generation parameters), `seed.rs` (deterministic seed management).
- Add a `GridF32` type (or similar) wrapping a `Vec<f32>` with width/height and indexed access, used throughout the pipeline as the core heightmap representation. Include basic operations: sample with bilinear interpolation, gradient (Sobel), min/max/mean statistics, load from PNG, save to PNG.
- Add a global `WorldSeed` struct that derives sub-seeds per phase (erosion seed, tectonic seed, noise seed) from a master u64 seed, ensuring phase independence (changing erosion parameters doesn't change tectonic output).
- README.md with project description, build instructions, and link to the TDD.
- `.gitignore` for Rust (target/, *.swp, etc.) and generated output files.

### Acceptance criteria

`cargo build` compiles both crates. `cargo run -p ymir-viz` opens an empty Bevy window. `cargo test -p ymir-core` runs and passes (even if there are only placeholder tests). The module structure exists even if modules are mostly empty. `GridF32` can load a PNG heightmap and save it back without data loss.

### Labels

`chore`, `M0`

---

## Issue #2 — Noise functions (FBM, ridged, rotated)

### Description

Implement the procedural noise functions needed across the pipeline: hash, value noise, FBM (standard and rotated), and ridged FBM. These functions must be deterministic (same input = same output, no thread-dependent variation) and performant (they are called millions of times during terrain enrichment).

### Tasks

- Implement `hash21(x: f32, y: f32) -> f32`: 2D to 1D hash function.
- Implement `noise2d(x: f32, y: f32) -> f32`: 2D value noise with smooth Hermite interpolation.
- Implement `fbm(x: f32, y: f32, octaves: u32) -> f32`: standard fractional Brownian motion, summing octaves with halving amplitude and doubling frequency.
- Implement `fbm_rotated(x: f32, y: f32, octaves: u32) -> f32`: FBM with ~37° rotation between octaves (matrix [0.8, 0.6; -0.6, 0.8] × 2.0) to reduce grid-alignment artifacts.
- Implement `ridged_fbm(x: f32, y: f32, octaves: u32) -> f32`: ridged multifractal noise where `abs(noise - 0.5) * 2.0` creates sharp ridge lines. This is the key noise type for mountain chains in the anisotropic FBM (M2).
- All functions take f32 and return f32 in approximately [0, 1] range.
- Unit tests: determinism (same input = same output across calls), range (output stays within expected bounds), visual plausibility (generate a 256×256 noise image and verify it looks correct).

### Acceptance criteria

All noise functions compile, pass unit tests for determinism and range. A test binary can generate a 512×512 PNG of each noise type for visual inspection.

### Labels

`feat`, `M0`

---

## Issue #3 — Hydraulic erosion engine

### Description

Implement the particle-based hydraulic erosion algorithm (Beyer 2015). This is the core M0 deliverable — the algorithm that carves valleys, creates drainage networks, and deposits sediment on the heightmap.

### Context

The algorithm simulates individual water droplets: each drops at a random position, follows the terrain gradient downhill, erodes material from steep slopes, transports sediment, and deposits when velocity decreases. Millions of droplets collectively carve a coherent drainage network.

### Tasks

- Define `ErosionConfig` struct with all tunable parameters: `erosion_rate` (0.3-0.5), `deposition_rate` (0.3-0.5), `evaporation_rate` (0.01-0.02), `inertia` (0.05-0.1), `gravity` (4.0-10.0), `min_slope` (0.01), `max_lifetime` (100-200 steps), `num_droplets` (millions), `coastal_deposition_range` (steps a droplet survives below sea level, for delta formation). Provide `Default` impl with recommended starting values.
- Implement the droplet simulation loop: spawn at random position (seeded RNG), compute bilinear gradient, update velocity (inertia-weighted blend of old direction and gradient direction), move, compute carrying capacity (velocity × water × slope), erode or deposit based on capacity vs current sediment load, evaporate water, terminate on evaporation or grid exit.
- The erosion modifies the heightmap in-place. A separate sediment accumulation map tracks total deposited material per cell.
- Implement coastal deposition: when a droplet reaches altitude ≤ 0 (sea level), it continues for `coastal_deposition_range` steps before dying, depositing its remaining sediment in the coastal zone. This creates river deltas.
- Parallelization strategy: process droplets in deterministic batches (batch size ~10,000-50,000). Within each batch, droplets run sequentially on the shared heightmap (to maintain determinism). Batches are processed with rayon but in a deterministic order (sorted by batch index). The heightmap is synchronized between batches.
- Sediment output: after all erosion, the sediment map records where material accumulated. High values = alluvial plains and deltas.

### Acceptance criteria

Running erosion on a 512×512 test heightmap (a simple cone or Perlin noise mountain) produces visible valleys that converge into a drainage network. The output is deterministic: same seed + same config = identical heightmap byte-for-byte. The sediment map shows deposits at the base of slopes and at the coastline. Performance: 1 million droplets on a 1024×1024 grid completes in under 30 seconds.

### Labels

`feat`, `M0`

---

## Issue #4 — Flow accumulation and river extraction

### Description

After erosion, compute the flow accumulation map (how much water passes through each cell) and extract the river network as a structured graph. This is the foundation for river rendering in the game and for navigability classification.

### Tasks

- Implement drainage direction computation: for each cell, determine which of the 8 neighbors is the lowest (steepest descent). Store as a u8 direction map (0-7 for the 8 cardinal/diagonal directions). Handle flat areas and local minima (cells with no lower neighbor) by routing to the nearest cell that does have a lower neighbor (priority flood fill or similar).
- Implement flow accumulation: sort all cells by decreasing altitude. For each cell (from highest to lowest), add its accumulated flow (starting at 1.0 for rainfall) to its downstream neighbor (from the direction map). After processing all cells, each cell's flow value represents the total upstream area draining through it.
- River thresholding: cells with flow > `river_threshold` are classified as river. Multiple thresholds define navigability classes: stream (small boats), river (barges), major river (ships).
- River vectorization: trace connected river cells into polyline segments. Each segment has: list of coordinates, average flux, average gradient, upstream/downstream connections to other segments, and basin ID (all cells draining to the same ocean outlet share a basin ID).
- Export the river network as a JSON structure (`rivers.json` per the TDD format).

### Acceptance criteria

Flow accumulation on the eroded heightmap produces a dendritic (tree-like) network that visually matches the valleys carved by erosion. Rivers start in mountains and end at the ocean. Tributaries merge into larger rivers (flux increases downstream). No rivers split (except at deltas, which are acceptable). The vectorized river graph is acyclic (no loops). Basin IDs correctly partition the continent into distinct drainage basins.

### Labels

`feat`, `M0`

---

## Issue #5 — Bevy visualization (M0)

### Description

Create the interactive Bevy application for visualizing erosion results and tuning parameters. This is the primary development tool — without it, tuning erosion parameters is guesswork.

### Tasks

- Bevy app setup: 2D camera with pan (drag) and zoom (scroll wheel), sized to display the full heightmap.
- Heightmap rendering: display the GridF32 as a colored quad using hypsometric tinting (blue below sea level, green lowlands 0-200m, yellow-green hills 200-600m, brown mountains 600-1500m, gray-white peaks 1500m+, white snow 2500m+). The color ramp is a configurable gradient.
- Hillshading: compute shading from the heightmap gradient (Sobel) with a configurable light direction. Multiply over the hypsometric colors for a relief effect.
- River overlay: display the flow accumulation as blue lines (alpha proportional to log(flux)) over the terrain. Togglable on/off.
- Erosion progress: run the erosion on a background thread (or in Bevy's async task system). Update the displayed heightmap periodically (every N batches of droplets) so the developer sees the terrain being carved in real-time.
- Parameter panel: a simple UI panel (Bevy UI or egui integration) with sliders for the main ErosionConfig parameters. A "Run erosion" button that starts/restarts the erosion with current parameters. A "Reset" button that reloads the original heightmap.
- PNG export: a button or keyboard shortcut (Ctrl+S) that saves the current heightmap and flow accumulation as PNG files.
- View modes: keyboard shortcuts to cycle between views: hypsometric altitude, hillshade only, slope magnitude, flow accumulation heatmap.

### Acceptance criteria

The developer can load a heightmap PNG, see it rendered with hypsometric colors and hillshading, adjust erosion parameters via sliders, run erosion and watch valleys appear in real-time, inspect the river network overlay, and export the result. The app runs at interactive frame rates (30+ FPS) even while erosion is computing in the background.

### Labels

`feat`, `M0`

---

## Issue #6 — PNG I/O and Azgaard heightmap compatibility

### Description

Implement loading of the existing Azgaard heightmap PNG and export of results in formats compatible with the current Living Landz server pipeline. This enables immediate testing of erosion results in-game without any server modifications.

### Tasks

- Load Azgaard heightmap PNG (grayscale u8, ~1920×1006) into GridF32, normalizing to a configurable elevation range (0-2500m). Handle the land/ocean distinction: pixels at value 0 in the binary map are sea level, land pixels start at a minimum altitude. Apply the same normalization as the current server pipeline (subtract min_land, divide by range).
- Load the binary map and lake map PNGs to mask the heightmap (ocean = 0, lakes = 0 or special value).
- Export the eroded heightmap as a u8 PNG matching the Azgaard format (for drop-in replacement in the current server) and as a u16 raw file (for the enriched heightmap pipeline).
- Export a flow accumulation visualization as PNG (log-scaled, blue channel).
- Export a "effective binary map" PNG derived from the eroded heightmap (altitude > 0 = white, else black) — this replaces the Azgaard binary map with coastlines modified by erosion and delta deposition.

### Acceptance criteria

Loading the Gaulyia heightmap + binary map + lake map produces a GridF32 where ocean is at 0, lakes are at 0, and land ranges from near-0 (coast) to ~1.0 (peaks). Running erosion and exporting as PNG produces files that the Living Landz server can load via the existing `WorldMaps::load` path without code changes. The exported binary map reflects any coastline changes from erosion (deltas extending the coast, etc.).

### Labels

`feat`, `M0`

---

---

# M1 Issues

---

## Issue #7 — Plate initialization and Voronoï partitioning

### Description

Generate the initial plate configuration: a set of tectonic plates as Voronoï regions on a 2D grid, each with a type (continental or oceanic) and a velocity vector. This is the creative input to the tectonic simulation — different plate configurations produce different continent shapes.

### Tasks

- Implement Voronoï partitioning on a 2D grid with periodic boundary conditions (torus topology). Place N seed points (5-15) randomly (seeded RNG), then assign each grid cell to the nearest seed using a distance metric that wraps around the grid edges.
- Each plate gets: an ID, a type (Continental or Oceanic), and a velocity vector (Vec2). The ratio of continental to oceanic plates should be configurable (default ~30% continental).
- Define a `PlateConfig` struct: `num_plates`, `continental_ratio`, `velocity_range` (min/max speed), `seed`. Provide presets for common configurations: "single continent" (one large continental plate surrounded by oceanic), "collision" (two continental plates converging), "archipelago" (several small continental plates), "rift" (one continental plate with divergent velocities pulling it apart).
- Initialize the crustal thickness field: continental plate cells get S=1.0, oceanic cells get S=0.2. Apply a gentle Gaussian blur (sigma=2-3 cells) at plate boundaries to smooth the sharp transitions.
- Visualization: color each cell by plate ID, overlay velocity arrows, highlight continental vs oceanic plates with different color saturation.

### Acceptance criteria

Voronoï partitioning works correctly with periodic boundaries (no edge artifacts, plates can wrap around the grid). The same seed produces identical plate configurations. Preset configurations produce recognizable patterns (the "single continent" preset gives one large landmass, the "collision" preset shows two converging masses). The crustal thickness field is smooth at plate boundaries.

### Labels

`feat`, `M1`

---

## Issue #8 — Thin viscous sheet solver

### Description

Implement the core tectonic simulation: the thin viscous sheet model (England & McKenzie 1982) that solves for the velocity field from force balance and evolves the crustal thickness over time. This is the most mathematically challenging component of the entire project.

### Context

The equation to solve at each timestep is a force balance: plate traction + gravitational pressure gradient (from crustal thickness variations) + viscous resistance = 0. This gives a velocity field. The crustal thickness is then advected by this velocity field: `∂S/∂t = -∇·(S·v) + sources`.

### Tasks

- Implement the gravitational pressure gradient force: `F_grav = -∇(S²)`. This force drives thick crust to spread laterally, preventing unlimited buildup.
- Implement the viscous stress term. Start with linear viscosity (n=1, simpler) for the first iteration: `F_visc = η∇²v` where η is the viscosity parameter. The power-law version (n=3, more realistic) can be added later as a refinement.
- Implement the velocity solver: discretize the force balance on the grid using finite differences. The resulting linear system `A·v = b` is solved iteratively (Gauss-Seidel with SOR, or conjugate gradient). Periodic boundary conditions mean the grid wraps. The solver runs until convergence (residual < threshold) or max iterations.
- Implement the advection step: `S_new = S_old - dt * ∇·(S·v)`. Use an upwind scheme for stability (the flux at each cell face uses the S value from the upwind side).
- Implement the CFL condition for timestep stability: `dt < dx / max(|v|)` where dx is the grid spacing.
- Define `TectonicsConfig`: `viscosity`, `gravity_factor` (strength of the gravitational spreading), `num_timesteps` (200-500), `convergence_threshold` for the velocity solver, `power_law_exponent` (1 for linear, 3 for realistic).
- Logging: at each timestep, log max velocity, max/min crustal thickness, solver iterations to convergence. This is essential for debugging numerical instabilities.

### Acceptance criteria

The velocity solver converges within a reasonable number of iterations (< 100 for a 128² grid). The advection is stable (no NaN, no negative crustal thickness, no unbounded growth). Two converging continental plates produce thickened crust at the collision front. Two diverging plates produce thinned crust at the rift. The simulation runs for 300 timesteps on a 128² grid in under 30 seconds.

### Labels

`feat`, `M1`

### References

England & McKenzie (1982), Houseman & England (1986), Flesch et al. (2001). See TDD §4.

---

## Issue #9 — Boundary processes (subduction, volcanism, rifting)

### Description

Implement the specific physical processes that occur at plate boundaries: subduction with volcanic arc formation, oceanic crust creation at divergent boundaries, and boundary type detection. These processes create the distinctive features (volcanic arcs, island chains, rift valleys) that make the tectonic output geologically interesting.

### Dependencies

Issue #8 (thin sheet solver) must be functional. This issue adds source/sink terms to the thickness evolution equation.

### Tasks

- Implement automatic boundary type detection: for each cell, compare the velocity of its plate with the velocities of neighboring plates. Convergent = velocities point toward each other. Divergent = velocities point away. Transform = velocities are parallel but in different directions. The boundary type + plate types (continental/oceanic) determine the geological process.
- Implement subduction: at convergent oceanic-continental boundaries, the oceanic plate's thickness is consumed (decreases toward the boundary). A volcanic source term is added at a distance of ~100-200 grid cells (configurable) behind the subduction front on the continental side: `source = volcanic_rate × exp(-distance_from_front / volcanic_decay)`. This creates a belt of thickened crust parallel to the coast — the volcanic arc.
- Implement oceanic crust creation: at divergent boundaries where both plates are oceanic (or where continental crust thins below a threshold), new oceanic crust is generated. Cells that thin below `S_ocean_min` (≈0.15) are reset to `S_ocean_new` (≈0.2), representing new ocean floor formation at a mid-ocean ridge.
- Implement rift volcanism: at divergent continental boundaries, where the crust is actively thinning, a small volcanic source term creates localized thickening (representing rift volcanism, like the East African volcanic centers). The source is proportional to the thinning rate.
- All source terms are added to the advection equation: `S_new = S_old - dt * ∇·(S·v) + dt * sources`.

### Acceptance criteria

A configuration with an oceanic plate subducting under a continental plate produces a volcanic arc 100-200 km inland from the coast. A divergent configuration produces a rift valley with thinned crust and occasional volcanic thickening. A mid-ocean divergent boundary creates new oceanic crust. These features are visible in the crustal thickness visualization and translate to recognizable altitude features after isostasy (issue #10).

### Labels

`feat`, `M1`

---

## Issue #10 — Isostasy and altitude conversion

### Description

Convert the crustal thickness field from the tectonic simulation into a surface altitude field using the Airy isostasy model. This transforms the abstract thickness values into a heightmap that can be visualized as terrain and eventually fed into the erosion pipeline.

### Tasks

- Implement the isostatic altitude formula: `altitude = (S - S_ref) × (ρ_mantle - ρ_crust) / ρ_mantle × scale_factor`. Parameters: `S_ref` (reference oceanic thickness, ~0.2), `ρ_mantle` (3.3), `ρ_crust` (2.7), `scale_factor` (calibrated so collision zones at S=1.8 produce ~3000-4000m peaks).
- Sea level determination: compute the altitude value that separates land from ocean. All cells below sea level are ocean (negative altitude = bathymetry). The sea level can be calibrated to achieve a target land/ocean ratio (~30% land).
- Apply a gentle Gaussian blur (sigma=2-3 grid cells) on the resulting altitude field to smooth the sharpest tectonic transitions. This prevents aliasing when the field is later upscaled.
- Export the altitude field as a GridF32, with positive values for land and negative for ocean floor. Also export a tectonic classification map: each cell gets a label (StableContinental, CollisionZone, SubductionArc, RiftZone, OceanicCrust, VolcanicArc) based on its tectonic history (thickness change patterns and boundary proximity).
- Define `IsostasyConfig`: density values, scale_factor, sea_level_target_land_ratio, smoothing_sigma.

### Acceptance criteria

Continental plates stand above sea level (positive altitude). Oceanic plates are below sea level (negative altitude). Collision zones produce the highest peaks (2000-4000m depending on calibration). Rift zones produce below-average altitude (grabens). The tectonic classification map correctly labels boundary zones. The land/ocean ratio is approximately 30% (configurable).

### Labels

`feat`, `M1`

---

## Issue #11 — Bevy visualization (M1 — tectonic views)

### Description

Extend the Bevy visualization app to display the tectonic simulation, including real-time animated plate evolution and altitude output.

### Dependencies

Issue #5 (M0 Bevy viz) provides the base app. This issue adds tectonic-specific views.

### Tasks

- Tectonic thickness view: display the crustal thickness field as a colormap (thin/blue 0.1-0.3 → normal/green 0.8-1.2 → thick/red 1.5-2.0). Update in real-time as the simulation advances.
- Plate overlay: display plate boundaries as colored lines. Show velocity vectors as arrows on each plate. Color plates by ID with continental plates in warm tones and oceanic in cool tones.
- Altitude view: after isostasy, display the altitude field with the same hypsometric coloring and hillshading as the erosion view (issue #5). This allows direct visual comparison between the tectonic altitude (before erosion) and the eroded altitude (after M2).
- Tectonic classification view: color cells by their tectonic label (stable=gray, collision=red, subduction=orange, rift=purple, volcanic=yellow, oceanic=blue).
- Simulation controls: play/pause/step buttons. Speed slider (timesteps per frame). Reset button to restart with new seed or plate configuration.
- Configuration panel: plate count slider, continental ratio slider, velocity range slider. A "randomize" button that generates a new plate configuration. Preset buttons for common configurations (from issue #7).
- Timeline: a progress bar or counter showing the current timestep out of total.

### Acceptance criteria

The developer can watch plates collide and mountains form in real-time. The simulation can be paused, stepped, and resumed. Switching between thickness, altitude, and classification views is instantaneous. The developer can randomize plate configurations and quickly evaluate whether the resulting continent is interesting before running the full simulation.

### Labels

`feat`, `M1`

---

## Issue #12 — Seed determinism infrastructure and testing

### Description

Implement and validate the seed determinism system across both M0 and M1 components. The same master seed must produce identical output across runs, which is essential for iterative parameter tuning ("change one thing, keep everything else identical").

### Tasks

- Implement `WorldSeed` struct: takes a master u64 seed and derives independent sub-seeds for each pipeline phase using a deterministic hash (e.g., `phase_seed = hash(master_seed, phase_name)`). Phases: "plates", "tectonics", "noise_meso", "noise_micro", "erosion", "climate".
- Audit all random operations in M0 (erosion droplet placement) and M1 (plate seed placement, velocity assignment) to ensure they use the seeded RNG from WorldSeed, not thread-local or system RNG.
- Audit the erosion batch processing to ensure deterministic ordering: batches are indexed, processed in index order, and each batch's droplets are generated from the batch-specific seed offset.
- Implement a determinism test: run the full M0 pipeline (load heightmap → erode → flow accumulation) twice with the same seed. Compare output byte-for-byte. They must be identical.
- Implement the same test for M1: run tectonics twice with the same seed. Identical output.
- Add a CI-friendly test that catches regressions: save a reference output hash for a known seed, verify it matches after code changes.

### Acceptance criteria

All determinism tests pass. Running `ymir-viz` twice with `--seed 42` produces visually and numerically identical terrain. Changing the erosion config while keeping the same tectonic seed produces identical tectonic output (phase independence).

### Labels

`test`, `M0`, `M1`