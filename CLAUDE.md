# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Ymir is a physically-grounded continent generator written in Rust. It simulates geological processes (tectonics, erosion, climate) to produce terrain where every feature has a causal explanation. Proprietary license — no redistribution.

## Commands

```bash
cargo build --release                # Always use --release (erosion is 10-20x slower in debug)
cargo test --workspace               # Run all tests
cargo run -p ymir-viz --release      # Launch Bevy visualization app
cargo test -p ymir-core              # Test core library only
cargo fmt --all                      # Format (edition 2021, max_width 100)
cargo clippy --workspace             # Lint
```

## Architecture

Rust workspace (edition 2024, requires rustc 1.80+) with two crates:

- **ymir-core** — library containing all generation logic, no UI dependencies. Pipeline phases are independent modules that read/write `GridF32` (the core 2D heightmap type in `grid.rs`).
- **ymir-viz** — Bevy 0.18.1 binary for interactive visualization during development.

### Generation pipeline (six phases)

1. **Tectonics** (`tectonics/`) — thin viscous sheet simulation at 128²–512²
2. **Isostasy** (`tectonics/isostasy.rs`) — crustal thickness → altitude
3. **Upscale + detail** (`terrain/upscale.rs`, `terrain/noise.rs`) — bicubic interpolation + anisotropic FBM to 4096²–8192²
4. **Erosion** (`erosion/`) — hydraulic, thermal, coastal, aeolian, glacial
5. **Climate** (`climate/`) — temperature, precipitation, Whittaker biome classification
6. **Export** (`export/`) — PNG heightmaps and raw data files

### Key types

- **`GridF32`** — row-major `Vec<f32>` grid with interpolation, Sobel gradients, statistics, PNG I/O
- **`WorldSeed`** — master seed → phase-specific RNGs via `ChaCha8Rng`. Changing parameters in one phase does not affect other phases' randomness.
- **`GenerationConfig`** — serializable master config grouping all pipeline parameters

### Design invariants

- **Deterministic**: same seed + config = identical output. Rayon batches are sorted before processing.
- **Phase-independent seeding**: `WorldSeed` derives sub-seeds per phase so parameter tweaks don't cascade.
- **Horizontal layers**: each pipeline phase is a standalone module communicating via `GridF32`.

## Key dependencies

rayon (parallelism), rand + rand_chacha (deterministic RNG), serde + serde_json (config serialization), image (PNG I/O), bevy (viz only).

## Docs

- `docs/tdd.md` — Technical Design Document with detailed specifications for each phase
- `docs/milestones.md` — Issue roadmap (M0–M6)
- `docs/thin_viscous_sheet_resolution.md` — Thin viscous sheet resolution plan
- `docs/adr/0001-erosion-coastal-sink-and-terraces.md` — the erosion/hydrology ADR: every
  finding, its measurement, and the method rules earned along the way

## Tooling

**BANNED: `Get-Content -Raw` piped into `Add-Content` (PowerShell) for any file containing
non-ASCII text.** PowerShell 5.1's `Get-Content` reads UTF-8 as cp1252 and `Add-Content
-Encoding utf8` re-encodes it, so every non-ASCII character is DOUBLE-ENCODED. It corrupted
275 lines across the ADR and two source files, and it is silent — the code still compiled and
`grep` matched nothing unusual. Worse, cp1252 leaves `0x81/0x8D/0x8F/0x90/0x9D` undefined, so
those bytes are DESTROYED rather than transformed: `∝`, `↔` and superscripts cannot be
recovered by inverting the encoding and must be retyped. Use the Write/Edit tools, or Python
with an explicit `encoding='utf-8'`, to append to any such file.
