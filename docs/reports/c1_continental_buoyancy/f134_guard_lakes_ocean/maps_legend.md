# Legend — `d_lake{1,2,11}_dz.png` (Finding 134-D)

Written by `f134_d` (`crates/ymir-core/tests/f126_coast.rs`).

- **Crop**: centred on the lake's ON centroid, ± 307 cells (615 × 615 cells = 30 km at 48.8 m per cell),
  viz cell coordinates. The lake's centre is the image centre.
  - Lake 1: (3491, 2500).
  - Lake 2: (4019, 2839).
  - Lake 11: (1429, 4287).
- **Orientation: SOUTH UP.** The rows are the DATA rows, and the field is stored south-first (y = 0 is the south,
  Finding 27). The image is the viz's north-up view flipped vertically: east is right, south is up.
  *Corrected in Finding 135: this file first said "north up".*
- **Colour**: Δz = z_ON − z_OFF on the CONDITIONED field (`HdResult.eroded`). Its pre-breach lakes sit at
  their water surface.
  - **Red**: ON higher than OFF. **Blue**: ON lower. **White**: equal.
  - The scale is linear from 0 to ± 1 000 m. Beyond ± 1 000 m it is **saturated** (full red / full blue):
    a saturated red region means "ON is at least 1 km higher".
  - No colour bar is drawn in the image; this file is the scale.
- **Outlines**, drawn on the cells whose lake id differs from a 4-neighbour's:
  - **green**: the ON lake the map is about (the based lake);
  - **grey**: every other ON lake;
  - **black**: the OFF lakes.

  Where an ON and an OFF outline coincide, ON is drawn over OFF.
- **The worlds**:
  - OFF = `ValleyConstruction::new(F121_AGE_K, Some(0.1))` with slope floor 0.024, the viz's
    "C2 /10 col (défaut)".
  - ON = the same with `lake_base = InputLakesAndBasins`.
