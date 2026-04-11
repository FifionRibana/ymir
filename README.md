# Ymir

Physically-grounded continent generator. Produces terrain through tectonic plate simulation, hydraulic erosion, and climate modeling — from crustal forces to river valleys.

Named after the Norse primordial giant whose body became the world: his flesh the earth, his blood the oceans, his bones the mountains.

## What Ymir does

Ymir generates continents for games and worldbuilding. Rather than painting terrain with noise functions or artistic presets, it simulates the physical processes that create real landscapes. Tectonic plates collide to form mountain ranges. Water erodes the mountains into valleys. Rain falls on windward slopes and creates forests; leeward slopes stay dry and become grasslands. Rivers carve drainage networks that converge into navigable waterways. Minerals concentrate where the geology demands it.

The result is terrain where every feature has a cause. Valleys exist because water carved them. Mountains exist because plates collided. Biomes exist because the climate dictates them. This causal coherence produces terrain that feels natural at every scale, from continental geography down to individual valley floors where a village might settle.

## Pipeline overview

The generation proceeds in six phases, each building on the previous one.

**Phase 1 — Tectonics.** A thin viscous sheet model (England & McKenzie, 1982) simulates crustal deformation from plate motion on a low-resolution grid (128²–512²). Continental plates collide and thicken, oceanic plates subduct and create volcanic arcs, divergent boundaries rift apart. The output is a crustal thickness field.

**Phase 2 — Isostasy.** The crustal thickness is converted to surface altitude via the Airy isostasy formula. Thick crust stands high (mountains), thin crust sits low (basins), and the sea level is calibrated to achieve the desired land/ocean ratio.

**Phase 3 — Upscale and detail.** The tectonic altitude (128²–512²) is upscaled to the erosion grid (4096²–8192² at 35–40 meters per pixel) with bicubic interpolation. Anisotropic FBM noise adds directionally coherent detail: ridges parallel to mountain chains, perpendicular drainage texture on flanks.

**Phase 4 — Erosion.** Millions of simulated water droplets carve the terrain via hydraulic erosion (Beyer, 2015). Valleys form, ridges sharpen, sediment accumulates in floodplains and deltas. The flow accumulation map gives the river network directly. Lakes are detected and classified by geological context (glacial, volcanic, rift).

**Phase 5 — Climate and biomes.** Temperature (from altitude and distance to ocean) and precipitation (from orographic effects with prevailing wind) are computed per cell. A Whittaker diagram maps temperature × precipitation to biome type. Seasonal data (growing season, frost patterns) is exported for downstream use.

**Phase 6 — Export.** The generator outputs structured data files: heightmap, bathymetry, biomes, river network, geological classification, and climate data. The output format is documented in the [Technical Design Document](docs/tdd.md).

For full technical details, see the [Technical Design Document](docs/tdd.md) and the [Milestone Roadmap](docs/milestones.md).

## Project structure

```
ymir/
├── crates/
│   ├── ymir-core/       # Library: all generation logic (tectonics, erosion, climate)
│   └── ymir-viz/        # Binary: Bevy visualization app for interactive development
├── docs/                # TDD, milestones, research references
├── assets/              # Input data and plate configuration presets
└── output/              # Generated terrain data (gitignored)
```

Ymir is a Rust workspace with two crates. `ymir-core` contains the generation pipeline as a library with no UI dependencies — it can be used as a dependency by other tools or integrated into a game server. `ymir-viz` is a Bevy application that provides interactive visualization for development: real-time display of tectonic evolution, erosion progress, and terrain output, with parameter sliders for iterative tuning.

## Building

Ymir requires Rust 1.80+ (2024 edition). No external system dependencies beyond the Rust toolchain.

```bash
git clone https://github.com/FifionRibana/ymir.git
cd ymir
cargo build --release
```

To run the visualization app:

```bash
cargo run -p ymir-viz --release
```

To run the tests:

```bash
cargo test --workspace
```

Release mode (`--release`) is strongly recommended for the visualization app and for any actual terrain generation. The erosion simulation processes millions of particles and is 10–20× slower in debug mode.

## Usage

Ymir is an offline tool. The developer generates a continent, iterates on parameters until the result is satisfactory, then exports the data for consumption by a game engine or other application.

A typical workflow looks like this. Launch the visualization app. Choose or randomize a plate configuration. Run the tectonic simulation and watch the continent form. When the macro shape is satisfying, run the erosion phase and watch valleys carve into the mountains. Adjust parameters (erosion intensity, wind direction for climate, etc.) and re-run as needed. When everything looks good, export the data files.

## References

The generation pipeline is grounded in published geophysics and computer graphics research.

England, P., & McKenzie, D. (1982). "A thin viscous sheet model for continental deformation." *Geophysical Journal of the Royal Astronomical Society*, 70(2), 295–321.

Beyer, H. T. (2015). "Implementation of a Simple Erosion Model." Technische Universität München.

Cordonnier, G., et al. (2016). "Large Scale Terrain Generation from Tectonic Uplift and Fluvial Erosion." *Computer Graphics Forum*, 35(2), 165–175.

Whittaker, R. H. (1975). *Communities and Ecosystems.* Macmillan.

## License

Copyright © 2026 Olivier BRUNEAU. All rights reserved. See [LICENSE](LICENSE) for details.