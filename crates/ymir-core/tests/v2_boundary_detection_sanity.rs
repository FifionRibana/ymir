//! Step 6 — dynamic boundary detection on a prescribed velocity field.
//!
//! Velocity `v = (sin(2πx), 0)` → `div(v) = 2π cos(2πx)`. On an
//! all-oceanic domain, cells with `cos > threshold/2π` get `Rift`,
//! cells with `cos < -threshold/2π` get `OceanicSubduction`, cells
//! near zero stay `None`. This is the unit-level sanity check for
//! the flag classification rule.

use std::f64::consts::PI;
use ymir_core::tectonics_v2::boundaries::boundary_flag::{BoundaryFlag, BoundaryFlagField};
use ymir_core::tectonics_v2::boundaries::plate_type::{PlateType, PlateTypeField};
use ymir_core::tectonics_v2::boundary_detection::{detect_boundaries, DetectionConfig};
use ymir_core::tectonics_v2::field::PeriodicIndex;
use ymir_core::tectonics_v2::voronoi::PlateIdField;

#[test]
fn sinusoidal_vx_produces_mixed_rift_and_subduction() {
    let nx = 64;
    let ny = 8;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);
    let mut vx = vec![0.0; nx * ny];
    let vy = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let x_face = i as f64 * dx;
            vx[j * nx + i] = (2.0 * PI * x_face).sin();
        }
    }
    let plate_type = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
    let pid = PlateIdField::new(nx, ny);
    let mut out = BoundaryFlagField::filled(nx, ny, BoundaryFlag::None);
    let cfg = DetectionConfig { threshold: 1e-4 };
    detect_boundaries(
        nx, ny, dx, dy, &idx_x, &idx_y, &vx, &vy, &plate_type, &pid, &cfg, &mut out,
    );

    let mut rift = 0usize;
    let mut sub = 0usize;
    let mut none = 0usize;
    for i in 0..nx {
        match out.get(i, 0) {
            BoundaryFlag::Rift => rift += 1,
            BoundaryFlag::OceanicSubduction => sub += 1,
            BoundaryFlag::None => none += 1,
            _ => panic!("unexpected flag on oceanic-only domain"),
        }
    }
    assert!(rift > nx / 8, "expected many rift cells, got {}", rift);
    assert!(sub > nx / 8, "expected many subduction cells, got {}", sub);
    assert_eq!(rift + sub + none, nx);
    let _ = dy;
}

#[test]
fn zero_velocity_means_all_none() {
    let nx = 16;
    let ny = 16;
    let dx = 1.0 / nx as f64;
    let dy = 1.0 / ny as f64;
    let idx_x = PeriodicIndex::new(nx);
    let idx_y = PeriodicIndex::new(ny);
    let vx = vec![0.0; nx * ny];
    let vy = vec![0.0; nx * ny];
    let plate_type = PlateTypeField::filled(nx, ny, PlateType::Oceanic);
    let pid = PlateIdField::new(nx, ny);
    let mut out = BoundaryFlagField::filled(nx, ny, BoundaryFlag::Rift);
    let cfg = DetectionConfig::default();
    detect_boundaries(
        nx, ny, dx, dy, &idx_x, &idx_y, &vx, &vy, &plate_type, &pid, &cfg, &mut out,
    );
    for &f in out.data() {
        assert!(matches!(f, BoundaryFlag::None));
    }
    let _ = dy;
}
