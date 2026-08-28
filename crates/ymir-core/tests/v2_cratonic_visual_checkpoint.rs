//! Step 9 Phase 3 sanity checkpoint — generate PNGs of the
//! `cratonic_factor` field on a typical 64² Voronoï so a human can
//! eyeball that cratons land at the centres of the large continental
//! plates. Marked `#[ignore]` because it touches the filesystem and
//! is meant to be invoked explicitly:
//!
//! ```text
//! cargo test --release -p ymir-core --test v2_cratonic_visual_checkpoint -- --ignored --nocapture
//! ```
//!
//! Output PNGs land in `docs/reports/step9_visual_checkpoint/`:
//! - `cratonic_factor_64sq.png` — the factor field itself.
//! - `plate_type_64sq.png`     — companion plate-type map (oceanic 0,
//!                               continental 1) for context.
//! - `plate_id_64sq.png`       — companion plate-id map so the
//!                               cratonic centres can be matched to
//!                               the right plates.
//! - one colour-bar PNG per snapshot (saved automatically).
//!
//! The cratons must visibly land at the *interiors* of the large
//! continental plates and stay clear of the oceanic regions.

use std::path::PathBuf;

use ymir_core::tectonics_v2::cratonic::{
    CratonicConfigEnabled, factor::build_cratonic_factor_field,
};
use ymir_core::tectonics_v2::diagnostics::heightmap::save_heightmap;
use ymir_core::tectonics_v2::voronoi::{VoronoiConfig, generate_voronoi};

fn dump_for_grid(nx: usize, ny: usize, label: &str) {
    let vcfg = VoronoiConfig::default(); // num_plates=8, continental_ratio=0.3
    let plates = generate_voronoi(nx, ny, &vcfg, 42);
    let crcfg = CratonicConfigEnabled::default();
    let factor = build_cratonic_factor_field(&plates, &crcfg);

    // `cargo test` runs with CWD = crate root (`crates/ymir-core`).
    // Place outputs at workspace-root `docs/reports/...` so the
    // report markdown can reference them with the same relative
    // path used by the rest of the milestone.
    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/reports/step9_visual_checkpoint");
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    let factor_meta =
        save_heightmap(&factor, &out_dir.join(format!("cratonic_factor_{}.png", label)))
            .expect("save factor png");
    let pt_meta = save_heightmap(
        &plates.plate_type.to_heightmap(),
        &out_dir.join(format!("plate_type_{}.png", label)),
    )
    .expect("save plate-type png");
    let pid_meta = save_heightmap(
        &plates.plate_id.to_heightmap(),
        &out_dir.join(format!("plate_id_{}.png", label)),
    )
    .expect("save plate-id png");

    // Stat summary so the human running this test sees what they
    // are looking at without having to open the PNGs.
    let mut cratonic_count = 0usize;
    let mut nonzero_count = 0usize;
    for v in factor.data() {
        if *v > 0.0 {
            nonzero_count += 1;
        }
        if *v > 0.5 {
            cratonic_count += 1;
        }
    }
    let total = (nx * ny) as f64;
    let cratonic_frac = cratonic_count as f64 / total;
    let nonzero_frac = nonzero_count as f64 / total;

    let continental_count: usize = plates
        .plate_type
        .data()
        .iter()
        .filter(|&&t| matches!(t, ymir_core::tectonics_v2::boundaries::PlateType::Continental))
        .count();
    let continental_frac = continental_count as f64 / total;
    let expected_cratonic_frac = crcfg.cr * continental_frac;

    println!();
    println!("Step 9 cratonic_factor visual checkpoint — {}, seed=42", label);
    println!("  factor png            : {}", factor_meta.png_path.display());
    println!("  plate_type png        : {}", pt_meta.png_path.display());
    println!("  plate_id png          : {}", pid_meta.png_path.display());
    println!(
        "  factor range          : [{:.6}, {:.6}], mean = {:.6}",
        factor_meta.min, factor_meta.max, factor_meta.mean
    );
    println!("  cells with factor > 0   : {} ({:.2} %)", nonzero_count, 100.0 * nonzero_frac);
    println!(
        "  cells with factor > 0.5 : {} ({:.2} %) — cratonic_cell_fraction",
        cratonic_count,
        100.0 * cratonic_frac
    );
    println!("  continental fraction    : {:.2} %", 100.0 * continental_frac);
    println!(
        "  expected craton frac    : Cr * continental = {:.2} %",
        100.0 * expected_cratonic_frac
    );
    let rel_diff = if expected_cratonic_frac > 0.0 {
        (cratonic_frac - expected_cratonic_frac).abs() / expected_cratonic_frac
    } else {
        0.0
    };
    println!(
        "  relative diff vs expected : {:.1} % (acceptance #8 tolerates 20 %)",
        100.0 * rel_diff
    );
}

#[test]
#[ignore]
fn dump_cratonic_factor_64sq() {
    dump_for_grid(64, 64, "64sq");
}

/// Step 9 Phase 8 — companion `cratonic_factor` PNG at 32², the
/// resolution used for the Section 2 immunity demonstration on
/// Step 8 shape. Same Voronoï seed (42) so the plate layout is
/// recognisable across the 32² and 64² visuals.
#[test]
#[ignore]
fn dump_cratonic_factor_32sq() {
    dump_for_grid(32, 32, "32sq");
}
