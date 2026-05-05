//! Step 13 Phase 4 — sanity visual checkpoint for the new init modes.
//!
//! Generates four S̃ initialisation fields on the same Voronoï layout
//! (active_medley preset: seed=42, 8 plates, 30 % continental) at 64²
//! so the reviewer can validate visually that:
//!
//! 1. `Uniform` (reference) shows rigid Voronoï borders + uniform
//!    interior thickness — the limitation Step 13 addresses.
//! 2. `RadialProfile { Smoothstep, default }` produces gradual
//!    margins (interior peaks at `continental_value`, transitions
//!    smoothly to `oceanic_value` at inter-plate boundaries).
//! 3. `RadialProfile { Pow { 2.0 } }` keeps more cells closer to
//!    `oceanic_value` (sharper transition near the interior end).
//! 4. `RadialProfileWithFBM { default amplitude 0.10 }` adds intra-
//!    plate speckle on top of the radial profile.
//!
//! Marked `#[ignore]` because it touches the filesystem and is
//! invoked explicitly:
//!
//! ```text
//! cargo test --release -p ymir-core --test v2_step13_visual_checkpoint -- --ignored --nocapture
//! ```
//!
//! Output PNGs land in `docs/reports/step13_visual_checkpoint/`:
//! - `<mode>_64sq.png` — each individual field with its own
//!   `[min, max]` colour range (via `save_heightmap`) for stats.
//! - `patchwork_init_modes_64sq.png` — 1×4 patchwork at fixed
//!   `[0, 1]` range so the four modes are directly visually
//!   comparable.
//! - `plate_type_64sq.png`, `plate_id_64sq.png` — context maps
//!   (where the continental plates are, how the Voronoï cells lay
//!   out).
//!
//! Acceptance is human-eyeball: (2) must show gradient margins vs
//! (1)'s sharp Voronoï outlines; (3) must look "blockier" than (2);
//! (4) must show intra-plate variation on top of (2)/(3).

use std::path::PathBuf;

use ymir_core::grid::GridF32;
use ymir_core::tectonics_v2::boundaries::PlateType;
use ymir_core::tectonics_v2::diagnostics::heightmap::save_heightmap;
use ymir_core::tectonics_v2::field::Field2D;
use ymir_core::tectonics_v2::init::{
    init_s_field, FBM_AMPLITUDE_DEFAULT, FBM_LACUNARITY_DEFAULT, FBM_OCTAVES_DEFAULT,
    FBM_PERSISTENCE_DEFAULT, FBM_SCALE_DEFAULT, FBM_SEED_DEFAULT, InitContext, InitMode,
    PlateInitData, ProfileShape,
};
use ymir_core::tectonics_v2::voronoi::{VoronoiConfig, generate_voronoi};

const NX: usize = 64;
const NY: usize = 64;
const SEED: u64 = 42; // active_medley.json
const NUM_PLATES: usize = 8;
const CONTINENTAL_RATIO: f64 = 0.3;
const CONTINENTAL_VALUE: f64 = 0.95;
const OCEANIC_VALUE: f64 = 0.20;

#[test]
#[ignore]
fn dump_step13_init_modes_64sq() {
    let plates = generate_voronoi(
        NX,
        NY,
        &VoronoiConfig {
            num_plates: NUM_PLATES,
            continental_ratio: CONTINENTAL_RATIO,
        },
        SEED,
    );

    let plate_data = PlateInitData {
        plate_id: &plates.plate_id,
        plate_type: &plates.plate_type,
        seed_coords: Some(&plates.seed_coords),
    };
    let ctx = InitContext {
        nx: NX,
        ny: NY,
        seed: SEED,
        // Legacy Checkerboard amplitude — unused by the modes we
        // exercise here, but the context API requires a value.
        amplitude: 0.2,
        plate_data: Some(plate_data),
    };

    let modes: Vec<(&str, InitMode)> = vec![
        (
            "uniform",
            InitMode::Uniform { boundary_smoothing_width: 1.0 },
        ),
        (
            "radial_smoothstep",
            InitMode::RadialProfile {
                continental_value: CONTINENTAL_VALUE,
                oceanic_value: OCEANIC_VALUE,
                profile_shape: ProfileShape::Smoothstep,
            },
        ),
        (
            "radial_pow_2_0",
            InitMode::RadialProfile {
                continental_value: CONTINENTAL_VALUE,
                oceanic_value: OCEANIC_VALUE,
                profile_shape: ProfileShape::Pow { exponent: 2.0 },
            },
        ),
        (
            "radial_fbm_default",
            InitMode::RadialProfileWithFBM {
                continental_value: CONTINENTAL_VALUE,
                oceanic_value: OCEANIC_VALUE,
                profile_shape: ProfileShape::Smoothstep,
                fbm_amplitude: FBM_AMPLITUDE_DEFAULT,
                fbm_octaves: FBM_OCTAVES_DEFAULT,
                fbm_persistence: FBM_PERSISTENCE_DEFAULT,
                fbm_lacunarity: FBM_LACUNARITY_DEFAULT,
                fbm_scale: FBM_SCALE_DEFAULT,
                fbm_seed: FBM_SEED_DEFAULT,
            },
        ),
    ];

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step13_visual_checkpoint");
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    println!();
    println!(
        "Step 13 Phase 4 visual checkpoint — 64² Voronoï \
         (active_medley layout: seed={}, {} plates, {:.0}% continental)",
        SEED,
        NUM_PLATES,
        CONTINENTAL_RATIO * 100.0
    );
    println!(
        "Output dir: {}",
        out_dir.canonicalize().unwrap_or(out_dir.clone()).display()
    );
    println!();

    let mut fields: Vec<(String, Field2D)> = Vec::new();
    for (label, mode) in &modes {
        let s = init_s_field(*mode, &ctx);
        let png_path = out_dir.join(format!("{}_64sq.png", label));
        let meta = save_heightmap(&s, &png_path).expect("save individual png");

        // Continental-only stats: mean and std-dev. These quantify
        // (a) how much the radial profile pulls the average toward
        // oceanic_value (Pow{2} → lower mean than Smoothstep →
        // lower than Uniform's flat continental_value); and (b)
        // the intra-plate variability that FBM introduces (std-dev
        // should be visibly larger for radial_fbm_default than for
        // the FBM-free modes).
        let mut cont_sum = 0.0_f64;
        let mut cont_count = 0usize;
        for j in 0..NY {
            for i in 0..NX {
                if matches!(plates.plate_type.get(i, j), PlateType::Continental) {
                    cont_sum += s.get(i, j);
                    cont_count += 1;
                }
            }
        }
        let cont_mean = cont_sum / cont_count.max(1) as f64;
        let mut cont_var = 0.0_f64;
        for j in 0..NY {
            for i in 0..NX {
                if matches!(plates.plate_type.get(i, j), PlateType::Continental) {
                    let d = s.get(i, j) - cont_mean;
                    cont_var += d * d;
                }
            }
        }
        let cont_std = (cont_var / cont_count.max(1) as f64).sqrt();
        println!(
            "  {:<22} : range=[{:.4}, {:.4}], mean={:.4}, continental_mean={:.4}, continental_std={:.4}",
            label, meta.min, meta.max, meta.mean, cont_mean, cont_std
        );

        fields.push((label.to_string(), s));
    }

    // Context maps.
    let pt_meta = save_heightmap(
        &plates.plate_type.to_heightmap(),
        &out_dir.join("plate_type_64sq.png"),
    )
    .expect("save plate type");
    let pid_meta = save_heightmap(
        &plates.plate_id.to_heightmap(),
        &out_dir.join("plate_id_64sq.png"),
    )
    .expect("save plate id");
    println!();
    println!("Context maps:");
    println!("  plate_type : {}", pt_meta.png_path.display());
    println!("  plate_id   : {}", pid_meta.png_path.display());

    // 1×4 patchwork at fixed [0, 1] range so the four init modes
    // are directly visually comparable. 1-pixel mid-grey separator
    // between tiles. Tiles are 64×64; final image 4·64 + 3·1 = 259
    // wide × 64 high.
    let tile_w = NX;
    let tile_h = NY;
    let sep = 1usize;
    let pw = tile_w * fields.len() + sep * fields.len().saturating_sub(1);
    let ph = tile_h;
    let mut patch = vec![0.5_f32; pw * ph];
    for (k, (_, s)) in fields.iter().enumerate() {
        let x_off = k * (tile_w + sep);
        for j in 0..tile_h {
            for i in 0..tile_w {
                let v = s.get(i, j) as f32;
                patch[j * pw + (x_off + i)] = v;
            }
        }
    }
    let patch_grid = GridF32::from_vec(pw, ph, patch);
    let patch_path = out_dir.join("patchwork_init_modes_64sq.png");
    patch_grid.save_png_u16(&patch_path).expect("save patchwork");
    println!();
    println!("Patchwork (fixed [0, 1] range, 1×4 layout, mid-grey separators):");
    println!("  {}", patch_path.display());
    println!(
        "  Layout (left → right): Uniform | RadialProfile{{Smoothstep}} | \
         RadialProfile{{Pow 2.0}} | RadialProfileWithFBM{{default}}"
    );
}

/// Phase 7 — multi-preset gallery. Render the four init modes on
/// two contrasting Voronoï layouts and assemble a 2×4 patchwork
/// (rows = presets, cols = modes). Single artefact, directly
/// embeddable in the physics report.
///
/// Presets:
/// - **`single_continent`** (seed=12, 4 plates, 50 % continental):
///   pedagogical — few large continental plates so each gradient
///   has many cells.
/// - **`convergence`** (seed=23, 6 plates, 40 % continental):
///   multi-plate — exercises the same machinery on a layout with
///   plate-plate contacts and smaller continents, closer to a
///   real run.
///
/// `active_medley`, `subduction`, and `divergence` are intentionally
/// out of scope for Step 13 — multi-preset coverage of Step 13
/// init modes is deferred to Step 12 (where the erosion pipeline
/// will exercise the same fields with more visual signal).
#[test]
#[ignore]
fn dump_step13_galerie_multi_preset() {
    struct PresetCfg {
        label: &'static str,
        seed: u64,
        num_plates: usize,
        continental_ratio: f64,
    }

    let presets = [
        PresetCfg {
            label: "single_continent",
            seed: 12,
            num_plates: 4,
            continental_ratio: 0.5,
        },
        PresetCfg {
            label: "convergence",
            seed: 23,
            num_plates: 6,
            continental_ratio: 0.4,
        },
    ];

    let modes: Vec<(&str, InitMode)> = vec![
        (
            "uniform",
            InitMode::Uniform { boundary_smoothing_width: 1.0 },
        ),
        (
            "radial_smoothstep",
            InitMode::RadialProfile {
                continental_value: CONTINENTAL_VALUE,
                oceanic_value: OCEANIC_VALUE,
                profile_shape: ProfileShape::Smoothstep,
            },
        ),
        (
            "radial_pow_2_0",
            InitMode::RadialProfile {
                continental_value: CONTINENTAL_VALUE,
                oceanic_value: OCEANIC_VALUE,
                profile_shape: ProfileShape::Pow { exponent: 2.0 },
            },
        ),
        (
            "radial_fbm_default",
            InitMode::RadialProfileWithFBM {
                continental_value: CONTINENTAL_VALUE,
                oceanic_value: OCEANIC_VALUE,
                profile_shape: ProfileShape::Smoothstep,
                fbm_amplitude: FBM_AMPLITUDE_DEFAULT,
                fbm_octaves: FBM_OCTAVES_DEFAULT,
                fbm_persistence: FBM_PERSISTENCE_DEFAULT,
                fbm_lacunarity: FBM_LACUNARITY_DEFAULT,
                fbm_scale: FBM_SCALE_DEFAULT,
                fbm_seed: FBM_SEED_DEFAULT,
            },
        ),
    ];

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step13_visual_checkpoint");
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    println!();
    println!("Step 13 Phase 7 — multi-preset gallery (64²)");
    println!(
        "Output dir: {}",
        out_dir.canonicalize().unwrap_or(out_dir.clone()).display()
    );

    let tile_w = NX;
    let tile_h = NY;
    let sep = 1usize;
    let cols = modes.len();
    let rows = presets.len();
    let pw = tile_w * cols + sep * cols.saturating_sub(1);
    let ph = tile_h * rows + sep * rows.saturating_sub(1);
    let mut patch = vec![0.5_f32; pw * ph];

    for (row, preset) in presets.iter().enumerate() {
        let plates = generate_voronoi(
            NX,
            NY,
            &VoronoiConfig {
                num_plates: preset.num_plates,
                continental_ratio: preset.continental_ratio,
            },
            preset.seed,
        );
        let plate_data = PlateInitData {
            plate_id: &plates.plate_id,
            plate_type: &plates.plate_type,
            seed_coords: Some(&plates.seed_coords),
        };
        let ctx = InitContext {
            nx: NX,
            ny: NY,
            seed: preset.seed,
            amplitude: 0.2,
            plate_data: Some(plate_data),
        };

        println!();
        println!(
            "  {} (seed={}, {} plates, {:.0}% continental):",
            preset.label,
            preset.seed,
            preset.num_plates,
            preset.continental_ratio * 100.0
        );
        for (col, (mode_label, mode)) in modes.iter().enumerate() {
            let s = init_s_field(*mode, &ctx);

            // Per-(preset, mode) PNG with dynamic range for stats.
            let png_path = out_dir
                .join(format!("galerie_{}_{}_64sq.png", preset.label, mode_label));
            let meta = save_heightmap(&s, &png_path).expect("save individual png");

            // Continental-only mean + std-dev for the report table.
            let mut cont_sum = 0.0_f64;
            let mut cont_count = 0usize;
            for j in 0..NY {
                for i in 0..NX {
                    if matches!(plates.plate_type.get(i, j), PlateType::Continental) {
                        cont_sum += s.get(i, j);
                        cont_count += 1;
                    }
                }
            }
            let cont_mean = cont_sum / cont_count.max(1) as f64;
            let mut cont_var = 0.0_f64;
            for j in 0..NY {
                for i in 0..NX {
                    if matches!(plates.plate_type.get(i, j), PlateType::Continental) {
                        let d = s.get(i, j) - cont_mean;
                        cont_var += d * d;
                    }
                }
            }
            let cont_std = (cont_var / cont_count.max(1) as f64).sqrt();
            println!(
                "    {:<22} : cont_mean={:.4}, cont_std={:.4}, range=[{:.4}, {:.4}]",
                mode_label, cont_mean, cont_std, meta.min, meta.max
            );

            // Place tile in patchwork.
            let x_off = col * (tile_w + sep);
            let y_off = row * (tile_h + sep);
            for j in 0..tile_h {
                for i in 0..tile_w {
                    let v = s.get(i, j) as f32;
                    patch[(y_off + j) * pw + (x_off + i)] = v;
                }
            }
        }
    }

    let patch_grid = GridF32::from_vec(pw, ph, patch);
    let patch_path = out_dir.join("galerie_multi_preset_64sq.png");
    patch_grid.save_png_u16(&patch_path).expect("save patchwork");
    println!();
    println!("Multi-preset patchwork (fixed [0, 1] range, 2×4 layout):");
    println!("  {}", patch_path.display());
    println!("  Rows  (top → bottom): single_continent | convergence");
    println!(
        "  Cols  (left → right): Uniform | RadialProfile{{Smoothstep}} | \
         RadialProfile{{Pow 2.0}} | RadialProfileWithFBM{{default}}"
    );
}
