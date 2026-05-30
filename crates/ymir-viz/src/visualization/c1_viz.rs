//! C1 raster visualization — field enum + render-from-snapshot
//! function. Bevy sprite / Image plumbing is wired in Stage E5
//! alongside the UI panels.
//!
//! ## C1Field — 6 modes (Issue #137 Q-S.5)
//!
//! Mirrors `V2Field` shape but for the C1 raw-fields snapshot:
//!
//! - `SThickness` — `S̃` crustal thickness, `[0, 3.0]` palette
//!   (Track A/B/D gallery continuity per Q-V.3).
//! - `Age` — cell age, auto-normalised per-frame.
//! - `PlateId` — categorical hash-mod hue palette (12-color HSV
//!   ramp; deterministic so same id renders the same color
//!   frame-to-frame; adjacent plate ids `0..7` map to distinct
//!   hues for the Phase 2 R7 default 8-plate init).
//! - `PlateType` — 2-color (cyan oceanic / beige continental).
//! - `Altitude` — Architecture C derived `compute_isostasy` +
//!   `apply_stein_stein_bathymetry`. **Stage E3 ships a stub
//!   showing S̃** (gallery-anchored derivation lands at Stage E4).
//! - `Cratonic` — binary grayscale (0/1) MVP standalone view.

use ymir_core::tectonics_c1::stats::C1StepStats;

use crate::bridge::c1::snapshot::C1Snapshot;
use crate::visualization::colormap::{age_colormap, cratonic_grayscale, hypsometric_colormap};

/// Field selector for the C1 raster view. 6 modes per Q-S.5.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum C1Field {
    #[default]
    SThickness,
    Age,
    PlateId,
    PlateType,
    /// Stage E3 stub renders S̃; Stage E4 wires Architecture C
    /// derivation (`compute_isostasy` + Stein-Stein).
    Altitude,
    Cratonic,
}

impl C1Field {
    pub const ALL: &'static [C1Field] = &[
        C1Field::SThickness,
        C1Field::Age,
        C1Field::PlateId,
        C1Field::PlateType,
        C1Field::Altitude,
        C1Field::Cratonic,
    ];

    pub fn label(self) -> &'static str {
        match self {
            C1Field::SThickness => "S̃ (crustal thickness)",
            C1Field::Age => "Age field",
            C1Field::PlateId => "Plate ID",
            C1Field::PlateType => "Plate type",
            C1Field::Altitude => "Altitude (post-isostasy)",
            C1Field::Cratonic => "Cratonic mask",
        }
    }

    pub fn legend_caption(self) -> &'static str {
        match self {
            C1Field::SThickness => "S̃: 0 → 3.0 (Track A/B/D palette)",
            C1Field::Age => "Age: 0 (ridge / rift-spawned) → run max",
            C1Field::PlateId => "Plate ID: 12-color hash-mod hue",
            C1Field::PlateType => "Oceanic (cyan) / Continental (beige)",
            C1Field::Altitude => {
                "Altitude (Stage E3 stub — Stage E4 wires Architecture C)"
            }
            C1Field::Cratonic => "Cratonic: false (black) / true (white)",
        }
    }
}

/// Track A/B/D gallery palette fixed range for `S̃`.
const S_VIZ_MAX: f64 = 3.0;

/// Number of distinct hues in the PlateId hash-mod palette.
/// Picked so adjacent plate ids `0..7` (default 8-plate init) and
/// any spawned by rifting (typically ≤ 12 total over a 300-step
/// run per Track D Stage V) map to visibly distinct hues.
const PLATE_ID_PALETTE_SIZE: usize = 12;

/// Render a snapshot to a row-major RGBA8 buffer (`nx × ny × 4`
/// bytes). Caller is responsible for the Bevy Image update —
/// that wiring lives in `C1VisualizationPlugin` (Stage E5).
///
/// The function is pure / deterministic: given the same snapshot
/// + field, it always produces the same bytes (load-bearing for
/// the view-switch-during-pause behaviour — Stage E3 W4).
pub fn field_to_rgba(snapshot: &C1Snapshot, field: C1Field) -> Vec<u8> {
    let nx = snapshot.nx;
    let ny = snapshot.ny;
    let n_cells = nx * ny;
    let mut rgba = vec![0_u8; n_cells * 4];

    match field {
        C1Field::SThickness => render_s_thickness(snapshot, &mut rgba),
        C1Field::Age => render_age(snapshot, &mut rgba),
        C1Field::PlateId => render_plate_id(snapshot, &mut rgba),
        C1Field::PlateType => render_plate_type(snapshot, &mut rgba),
        // Stage E3 stub: alias to SThickness. Stage E4 replaces
        // with Architecture C derivation.
        C1Field::Altitude => render_s_thickness(snapshot, &mut rgba),
        C1Field::Cratonic => render_cratonic(snapshot, &mut rgba),
    }

    rgba
}

fn render_s_thickness(snapshot: &C1Snapshot, rgba: &mut [u8]) {
    let nx = snapshot.nx;
    let ny = snapshot.ny;
    debug_assert_eq!(rgba.len(), nx * ny * 4);
    for j in 0..ny {
        for i in 0..nx {
            let v = snapshot.s[j * nx + i];
            let t = (v / S_VIZ_MAX).clamp(0.0, 1.0);
            let [r, g, b, a] = hypsometric_colormap(t);
            let k = (j * nx + i) * 4;
            rgba[k] = r;
            rgba[k + 1] = g;
            rgba[k + 2] = b;
            rgba[k + 3] = a;
        }
    }
}

fn render_age(snapshot: &C1Snapshot, rgba: &mut [u8]) {
    let nx = snapshot.nx;
    let ny = snapshot.ny;
    // Per-frame auto-normalise (age range varies wildly with the
    // density-advection pile-up factor; per Track D Stage V seed
    // 42 reaches ≈ 6958 cumulative max).
    let mut vmax = 0.0_f64;
    for &v in &snapshot.age {
        if v > vmax {
            vmax = v;
        }
    }
    let denom = vmax.max(1e-12);
    for j in 0..ny {
        for i in 0..nx {
            let v = snapshot.age[j * nx + i];
            let t = (v / denom).clamp(0.0, 1.0);
            let [r, g, b, a] = age_colormap(t);
            let k = (j * nx + i) * 4;
            rgba[k] = r;
            rgba[k + 1] = g;
            rgba[k + 2] = b;
            rgba[k + 3] = a;
        }
    }
}

fn render_plate_id(snapshot: &C1Snapshot, rgba: &mut [u8]) {
    let nx = snapshot.nx;
    let ny = snapshot.ny;
    for j in 0..ny {
        for i in 0..nx {
            let pid = snapshot.plate_id[j * nx + i];
            let [r, g, b, a] = plate_id_color(pid);
            let k = (j * nx + i) * 4;
            rgba[k] = r;
            rgba[k + 1] = g;
            rgba[k + 2] = b;
            rgba[k + 3] = a;
        }
    }
}

fn render_plate_type(snapshot: &C1Snapshot, rgba: &mut [u8]) {
    let nx = snapshot.nx;
    let ny = snapshot.ny;
    // 2-color palette: oceanic cyan / continental beige.
    const OCEANIC: [u8; 4] = [80, 160, 200, 255];
    const CONTINENTAL: [u8; 4] = [220, 200, 160, 255];
    for j in 0..ny {
        for i in 0..nx {
            let t = snapshot.plate_type[j * nx + i];
            let c = if t == 0 { OCEANIC } else { CONTINENTAL };
            let k = (j * nx + i) * 4;
            rgba[k..k + 4].copy_from_slice(&c);
        }
    }
}

fn render_cratonic(snapshot: &C1Snapshot, rgba: &mut [u8]) {
    let nx = snapshot.nx;
    let ny = snapshot.ny;
    for j in 0..ny {
        for i in 0..nx {
            let t = snapshot.cratonic_mask[j * nx + i] as f64;
            let [r, g, b, a] = cratonic_grayscale(t);
            let k = (j * nx + i) * 4;
            rgba[k] = r;
            rgba[k + 1] = g;
            rgba[k + 2] = b;
            rgba[k + 3] = a;
        }
    }
}

/// Deterministic 12-color HSV hash-mod palette for `PlateId`
/// view. Same `pid` always returns the same color (load-bearing
/// for the view-switch / no-scintillation contract — Stage E3 W1).
///
/// Hue computed via `(pid as f32 * 360.0 / 12.0) % 360.0`. Saturation
/// 0.7, value 0.8 — visually distinct adjacent hues at the
/// 30°-step granularity. Adjacent plate ids `0..7` (default
/// 8-plate init) and the new-rift ids spawned by Track D (up to
/// ~12 total in a 300-step run) all map to distinct hues.
fn plate_id_color(pid: u16) -> [u8; 4] {
    let h = (pid as f32 * (360.0 / PLATE_ID_PALETTE_SIZE as f32))
        % 360.0;
    hsv_to_rgba(h, 0.7, 0.8)
}

/// HSV → RGB conversion in `[0, 1]` saturation/value, hue in
/// degrees `[0, 360)`. Helper for [`plate_id_color`].
fn hsv_to_rgba(h: f32, s: f32, v: f32) -> [u8; 4] {
    let c = v * s;
    let h_prime = h / 60.0;
    let x = c * (1.0 - (h_prime % 2.0 - 1.0).abs());
    let (r1, g1, b1) = match h_prime as u32 {
        0 => (c, x, 0.0),
        1 => (x, c, 0.0),
        2 => (0.0, c, x),
        3 => (0.0, x, c),
        4 => (x, 0.0, c),
        _ => (c, 0.0, x),
    };
    let m = v - c;
    let r = ((r1 + m) * 255.0).round().clamp(0.0, 255.0) as u8;
    let g = ((g1 + m) * 255.0).round().clamp(0.0, 255.0) as u8;
    let b = ((b1 + m) * 255.0).round().clamp(0.0, 255.0) as u8;
    [r, g, b, 255]
}

/// Live diagnostics summary derived from a snapshot — helper for
/// the Stage E5 stats panel. Returns `(live_plates, subduction_cells,
/// accretion_merges, rifting_splits, rifting_cells_thinned)`.
pub fn snapshot_event_summary(
    snapshot: &C1Snapshot,
) -> (usize, usize, usize, usize, usize) {
    let stats: &C1StepStats = &snapshot.stats;
    (
        snapshot.live_plate_count,
        stats.subduction.cells_consumed,
        stats.accretion.merges_count,
        stats.rifting_split.splits_count,
        stats.rifting_thinning.cells_thinned,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dummy_snapshot(nx: usize, ny: usize) -> C1Snapshot {
        let n = nx * ny;
        C1Snapshot {
            step: 0,
            nx,
            ny,
            dx: 1.0 / nx as f64,
            dy: 1.0 / ny as f64,
            s: vec![0.5; n],
            age: vec![1.0; n],
            plate_id: vec![0; n],
            plate_type: vec![1; n],
            cratonic_mask: vec![0; n],
            num_plates: 1,
            live_plate_count: 1,
            plate_velocities: vec![(0.0, 0.0)],
            stats: C1StepStats::default(),
        }
    }

    #[test]
    fn render_produces_rgba8_buffer() {
        let snap = dummy_snapshot(8, 8);
        for &field in C1Field::ALL {
            let rgba = field_to_rgba(&snap, field);
            assert_eq!(
                rgba.len(),
                8 * 8 * 4,
                "field {:?}: rgba len {} != expected {}",
                field,
                rgba.len(),
                8 * 8 * 4
            );
        }
    }

    #[test]
    fn plate_id_palette_deterministic() {
        // Same pid always maps to the same color.
        assert_eq!(plate_id_color(0), plate_id_color(0));
        assert_eq!(plate_id_color(7), plate_id_color(7));
    }

    #[test]
    fn plate_id_palette_adjacent_ids_distinct() {
        // Adjacent ids in the default 8-plate init should be
        // visibly distinct (any RGB channel differs by ≥ 20).
        for pid in 0..PLATE_ID_PALETTE_SIZE as u16 - 1 {
            let a = plate_id_color(pid);
            let b = plate_id_color(pid + 1);
            let d = (a[0] as i32 - b[0] as i32).abs()
                + (a[1] as i32 - b[1] as i32).abs()
                + (a[2] as i32 - b[2] as i32).abs();
            assert!(
                d >= 20,
                "adjacent ids {pid}/{} too close: {a:?} vs {b:?} (sum-diff = {d})",
                pid + 1
            );
        }
    }

    #[test]
    fn plate_type_2_color_distinct() {
        let nx = 4;
        let ny = 4;
        let mut snap = dummy_snapshot(nx, ny);
        // First half oceanic, second half continental.
        for k in 0..nx * ny / 2 {
            snap.plate_type[k] = 0;
        }
        for k in nx * ny / 2..nx * ny {
            snap.plate_type[k] = 1;
        }
        let rgba = field_to_rgba(&snap, C1Field::PlateType);
        // Ocean cell + continental cell differ on R/G/B.
        let oceanic = &rgba[0..4];
        let continental = &rgba[(nx * ny / 2) * 4..(nx * ny / 2) * 4 + 4];
        assert_ne!(oceanic, continental);
    }

    #[test]
    fn altitude_stub_matches_s_thickness() {
        let snap = dummy_snapshot(4, 4);
        let alt = field_to_rgba(&snap, C1Field::Altitude);
        let s = field_to_rgba(&snap, C1Field::SThickness);
        assert_eq!(
            alt, s,
            "Stage E3 Altitude stub must alias SThickness exactly (Stage E4 will diverge)"
        );
    }
}
