//! ADR Finding 80 — the four gate tests for `StreamPowerConfig::base_level_floor`.
//!
//! A seam that leaks makes every measurement after it meaningless, so these run before any
//! 8192² build: byte-identical at `None`, the parameter reaches the cache key (with the
//! negative control that catches a key folding a timestamp), and a two-cell fixture where the
//! bound either holds or does not.
//!
//! Run: cargo test -p ymir-core --release --test base_level_floor

use ymir_core::cache::CacheKey;
use ymir_core::erosion::stream_power::{BaseLevelFloor, StreamPowerConfig, incise};
use ymir_core::grid::GridF32;
use ymir_core::tectonics_c1::cached_product::{drainage_key_windowed, eroded_key_full};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::closures::volcanism::VolcanismConfig;
use ymir_core::tectonics_c1::drainage::C1DrainageConfig;
use ymir_core::terrain::upscale::FbmUpscaleConfig;

const SEA: f32 = 0.5;
/// 1 norm unit in metres on the production vertical contract (2 · 1.1302 · 5000).
const NORM_TO_M: f32 = 11302.0;

fn upscale_cfg(floor: Option<f32>) -> FbmUpscaleConfig {
    let mut c = FbmUpscaleConfig::default();
    let mut sp = c.stream_power.unwrap_or_default();
    sp.base_level_floor = floor.map(|epsilon_m| BaseLevelFloor { epsilon_m });
    c.stream_power = Some(sp);
    c
}

fn keys(floor: Option<f32>) -> (String, String) {
    let ss = SteinSteinParams::default();
    let volc = VolcanismConfig::default();
    let tectonic = CacheKey::root().with("seed", &10_481_999_410_520_546_993u64);
    let ek = eroded_key_full(&tectonic, &ss, &upscale_cfg(floor), &volc);
    let dk = drainage_key_windowed(&ek, &C1DrainageConfig::default(), &ss, None, 400.0);
    (ek.digest().to_string(), dk.digest().to_string())
}

/// A ramp of land running into one column of sea, so every land cell's chain of receivers ends
/// on an ocean cell well below sea level. This is the coastal situation in miniature.
fn ramp(w: usize, h: usize, sea_depth_m: f32) -> GridF32 {
    let mut g = GridF32::new(w, h, 0.0);
    for y in 0..h {
        for x in 0..w {
            // land descending eastward from +40 m to +0.5 m, then the sea column
            g.data[y * w + x] = if x + 1 >= w {
                SEA - sea_depth_m / NORM_TO_M
            } else {
                // +40 m inland down to +2.0 m at the last land column, so the cell whose
                // receiver IS the sea starts at the round's +2 m.
                SEA + (40.0 - 38.0 * (x as f32 / (w - 2) as f32)) / NORM_TO_M
            };
        }
    }
    g
}

/// The round's fixture, dialled so `f = K*dt*A_km2^m / dist_m` is about 1000 -- the Courant
/// regime Finding 61 measured (1353 on the production seed), where the implicit update drives a
/// cell essentially onto its receiver in ONE step. `cell_km = 0.0488` is the 8192 cell.
fn cfg(floor: Option<f32>) -> StreamPowerConfig {
    let mut c = StreamPowerConfig::default();
    c.sea_level = SEA;
    c.cell_km = 400.0 / 8192.0;
    c.depth_scale_m = 5000.0;
    c.k = 4500.0;
    c.dt = 1.0;
    c.m = 0.5;
    c.n = 1.0;
    c.iterations = 2;
    c.threshold = 0.0;
    c.min_area_cells = 1.0; // every cell is a channel: the bound is what is under test
    c.mfd_exponent = None; // D8 accumulation, so `f` is the plain single-receiver figure
    c.diffusion = 0.0; // isolate the implicit update from the other lowering operators
    c.diffusion_substeps = 1;
    c.talus_passes = 0;
    c.lateral_erosion = 0.0;
    c.base_level_floor = floor.map(|epsilon_m| BaseLevelFloor { epsilon_m });
    c
}

#[test]
fn none_is_byte_identical() {
    // The gate's whole promise: with `None`, not one bit moves. Compared against a config
    // built WITHOUT ever touching the field, so a default that silently changed would show.
    let g = ramp(64, 16, 0.1);
    let legacy = incise(&g, &{
        let mut c = cfg(None);
        c.base_level_floor = None;
        c
    });
    let explicit_none = incise(&g, &cfg(None));
    assert_eq!(legacy.data.len(), explicit_none.data.len());
    let moved = (0..legacy.data.len())
        .filter(|&k| legacy.data[k].to_bits() != explicit_none.data[k].to_bits())
        .count();
    assert_eq!(moved, 0, "`None` is not byte-identical: {moved} cells differ");
}

#[test]
fn the_floor_reaches_the_cache_key() {
    let (e_off, d_off) = keys(None);
    let (e_on, d_on) = keys(Some(0.5));
    eprintln!("   eroded  {e_off} -> {e_on}\n   drainage {d_off} -> {d_on}");
    assert_ne!(e_off, e_on, "the eroded key does not see `base_level_floor`: a stale field ships");
    assert_ne!(d_off, d_on, "the drainage key does not see it");
    // epsilon itself must reach the key, or a recalibration is served stale
    let (e_small, _) = keys(Some(0.2));
    assert_ne!(e_on, e_small, "`epsilon_m` does not reach the key");
    // NEGATIVE CONTROL — without this the test passes on a key that folds a timestamp.
    assert_eq!(keys(None).0, e_off, "the key is not stable across identical calls");
    assert_eq!(keys(Some(0.5)).0, e_on, "the key is not stable across identical calls");
}

#[test]
fn the_bound_holds_and_its_negative_control_drowns() {
    let eps = 0.5f32;
    let g = ramp(64, 16, 0.1);
    let floor_norm = SEA + eps / NORM_TO_M;

    let bounded = incise(&g, &cfg(Some(eps)));
    let free = incise(&g, &cfg(None));

    // population: cells that were LAND before (every column but the sea one)
    let land: Vec<usize> = (0..g.data.len()).filter(|&k| g.data[k] > SEA).collect();
    let drowned_free = land.iter().filter(|&&k| free.data[k] <= SEA).count();
    let drowned_bounded = land.iter().filter(|&&k| bounded.data[k] <= SEA).count();
    let below_floor = land.iter().filter(|&&k| bounded.data[k] < floor_norm - 1e-9).count();
    let lifted = land.iter().filter(|&&k| bounded.data[k] > g.data[k] + 1e-9).count();
    eprintln!(
        "   {} land cells | free: {drowned_free} drowned | bounded: {drowned_bounded} drowned, \
         {below_floor} under the floor, {lifted} LIFTED above their start",
        land.len()
    );

    // NEGATIVE CONTROL FIRST: if the unbounded run does not drown anything, the fixture does
    // not exercise the defect and the assertion below would pass for the wrong reason.
    assert!(
        drowned_free > 0,
        "the fixture does not drown a single cell without the bound — it tests nothing"
    );
    assert_eq!(drowned_bounded, 0, "the bound let {drowned_bounded} land cells cross sea level");
    assert_eq!(below_floor, 0, "{below_floor} cells sit under `sea + epsilon`");
    // the bound may stop erosion; it must never deposit
    assert_eq!(lifted, 0, "the bound LIFTED {lifted} cells — it is depositing, not bounding");
}
