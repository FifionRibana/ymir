//! Integration test for the dynamic-remap heightmap saver.

use std::path::PathBuf;

use ymir_core::tectonics_v2::diagnostics::heightmap::save_heightmap;
use ymir_core::tectonics_v2::field::Field2D;

#[test]
fn dynamic_range_captures_wide_signal() {
    let mut f = Field2D::new(16, 16);
    for j in 0..16 {
        for i in 0..16 {
            let x = i as f64 / 15.0;
            let y = j as f64 / 15.0;
            f.set(i, j, 0.2 + 3.0 * x * y);
        }
    }
    let tmp: PathBuf = std::env::temp_dir().join("v2_heightmap_dyn.png");
    let md = save_heightmap(&f, &tmp).unwrap();
    assert!((md.min - 0.2).abs() < 1e-12);
    assert!((md.max - 3.2).abs() < 1e-12);
    // Colourbar PNG created next to the main PNG.
    assert!(md.colorbar_path.is_file(), "colourbar at {:?}", md.colorbar_path);
}

#[test]
fn constant_field_does_not_panic() {
    let f = Field2D::filled(8, 8, 1.0);
    let tmp: PathBuf = std::env::temp_dir().join("v2_heightmap_const.png");
    let md = save_heightmap(&f, &tmp).unwrap();
    assert_eq!(md.min, md.max);
}
