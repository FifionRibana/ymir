//! Coastline shape metrics — the instrument ADR Findings 52–54 built, promoted out of a test
//! file so every sweep uses ONE definition.
//!
//! **What it measures and why the obvious metrics do not.** A turn count is a function of vertex
//! spacing (Finding 52: the same curve reads 20 % or 30 % depending only on how finely it is
//! sampled) and it cannot see how FAR a spur runs. A global axial concentration is blind to
//! parallelism arranged in bundles (Finding 54: 0.05 global against 0.52 local, on a coast the
//! renders show as an unmistakable comb).
//!
//! **The primary quantities, and the anchored target** — measured on the coarse field, which the
//! author confirmed reads as a credible coastline:
//!
//! | | production (shipped) | REFERENCE / target |
//! |---|---|---|
//! | indentation DENSITY | 25.8–29.8 per 100 km | **~1 per 100 km** |
//! | length p90 | 3.5–4.6 km | **~24 km** |
//! | LOCAL parallelism | 0.52 | **~0.33** |
//!
//! ⚠️ **Spurs must become LONGER and ~25× FEWER.** A remedy that SHORTENS them moves away from
//! the target — the opposite of the intuition this chantier started with.

use crate::grid::GridF32;
use crate::terrain::contour::Polyline;

/// Search window for a spur, in km. ⚠️ This is an INSTRUMENT CEILING: a spur longer than this
/// cannot be found, so every length statistic saturates at it. At 8 km the `max` column read
/// 7.99–8.00 at every setting INCLUDING the reference, whose true p90 is 24.68 km — the
/// saturation hid the target itself (ADR Finding 54).
pub const MAX_SPUR_KM: f32 = 50.0;
/// A spur must run at least this far out and back to count.
pub const MIN_SPUR_KM: f32 = 1.0;
/// …and close through a neck no wider than this. Wider and it is a bay, not a spur.
pub const NECK_KM: f32 = 0.6;
/// Neighbours used for the LOCAL parallelism measure.
const LOCAL_K: usize = 8;

/// Coastline shape, in physical units.
#[derive(Debug, Clone, Copy)]
pub struct CoastShape {
    pub coast_km: f32,
    pub count: usize,
    /// PRIMARY — indentations per 100 km of coastline. Target ≈ 1.
    pub density_per_100km: f32,
    pub med_len_km: f32,
    /// PRIMARY — the long spurs are what is visible. Target ≈ 24 km.
    pub p90_len_km: f32,
    pub max_len_km: f32,
    /// Share of the coastline lying inside a spur.
    pub share_pct: f32,
    /// Reported only: CV of the spacing between consecutive spurs (Poisson reference 1.0).
    pub spacing_cv: f32,
    /// Reported only: GLOBAL axial concentration. Kept to document that it is blind — bundles
    /// pointing different ways cancel. Use `local_axis_r`.
    pub global_axis_r: f32,
    /// Median over spurs of the axial concentration among each spur and its 8 nearest
    /// neighbours. This is what "parallel fringes" means. Target ≈ 0.33.
    pub local_axis_r: f32,
}

fn pct(sorted: &[f32], f: f32) -> f32 {
    if sorted.is_empty() {
        return 0.0;
    }
    let i = (((sorted.len() - 1) as f64 * f as f64).round() as usize).min(sorted.len() - 1);
    sorted[i]
}

/// Measure a coastline. A SPUR is an excursion that returns on itself: an arc of at least
/// [`MIN_SPUR_KM`] whose endpoints lie within [`NECK_KM`] of each other.
#[must_use]
pub fn coast_shape(polys: &[Polyline], km_per_cell: f32) -> CoastShape {
    coast_shape_thresholds(polys, km_per_cell, MIN_SPUR_KM, NECK_KM)
}

/// [`coast_shape`] with the two length scales given explicitly, so ONE detector can be read at
/// more than one scale instead of a second detector being written beside it.
///
/// ## Why this exists (ADR 0001 Finding 75)
///
/// [`MIN_SPUR_KM`] is 1 km, which is **20.5 cells at 8192²** — so the default instrument cannot
/// see an indentation 2 to 5 cells across, which is the texture the eye reads as "fur" on the
/// rendered coast. Reading the same coastline at the kilometre scale and at the cell scale is
/// two settings of one detector, not two detectors: the spur definition, the arc walk, the
/// axis and the parallelism are shared, so the two columns are commensurable.
///
/// `coast_shape(p, k)` is exactly `coast_shape_thresholds(p, k, MIN_SPUR_KM, NECK_KM)` — asserted
/// in `coast_shape_default_is_the_parameterised_call`, because a refactor that silently moves the
/// shipped numbers would invalidate every coastal figure in this dossier.
/// ONE detected spur, with the geometry the aggregates are built from — exposed because the
/// spectral instrument (ADR Finding 76 block C) needs the spurs' POSITIONS along the curvilinear
/// abscissa, which no aggregate can carry. [`coast_shape_thresholds`] is computed from exactly
/// this list, so a spectrum and a count can never describe different populations.
#[derive(Debug, Clone, Copy)]
pub struct Spur {
    /// Index into the `polys` slice.
    pub poly: usize,
    /// First and last sample of the excursion, into that polyline.
    pub i: usize,
    pub j: usize,
    /// Curvilinear abscissa of the spur's root along its polyline, km.
    pub root_km: f32,
    pub len_km: f32,
    /// Midpoint of the neck, in CELL coordinates.
    pub mid: (f32, f32),
    /// Axial angle of the excursion (neck midpoint -> farthest sample), radians.
    pub axis_rad: f32,
}

/// The spur walk itself. Returns the spurs and the total coastline length in CELLS.
///
/// A SPUR is an excursion that returns on itself: an arc of at least `min_spur_km` whose
/// endpoints lie within `neck_km` of each other, capped at [`MAX_SPUR_KM`].
#[must_use]
pub fn coast_spurs(
    polys: &[Polyline],
    km_per_cell: f32,
    min_spur_km: f32,
    neck_km: f32,
) -> (Vec<Spur>, f32) {
    let neck_cells = neck_km / km_per_cell;
    let min_arc_cells = min_spur_km / km_per_cell;
    let mut out: Vec<Spur> = Vec::new();
    let mut coast_cells = 0.0f32;
    for (pi, pl) in polys.iter().enumerate() {
        let mut arc = vec![0.0f32; pl.len()];
        for i in 1..pl.len() {
            arc[i] = arc[i - 1]
                + ((pl[i].0 - pl[i - 1].0).powi(2) + (pl[i].1 - pl[i - 1].1).powi(2)).sqrt();
        }
        coast_cells += *arc.last().unwrap_or(&0.0);
        if pl.len() < 8 {
            continue;
        }
        let mut i = 0usize;
        while i < pl.len() {
            let mut best: Option<usize> = None;
            let mut j = i + 1;
            while j < pl.len() && (arc[j] - arc[i]) * km_per_cell < MAX_SPUR_KM {
                if arc[j] - arc[i] >= min_arc_cells {
                    let d = ((pl[j].0 - pl[i].0).powi(2) + (pl[j].1 - pl[i].1).powi(2)).sqrt();
                    if d <= neck_cells {
                        best = Some(j);
                    }
                }
                j += 1;
            }
            match best {
                Some(j) => {
                    let mid = ((pl[i].0 + pl[j].0) * 0.5, (pl[i].1 + pl[j].1) * 0.5);
                    let mut far = (0.0f32, (0.0f32, 0.0f32));
                    for p in &pl[i..=j] {
                        let d = ((p.0 - mid.0).powi(2) + (p.1 - mid.1).powi(2)).sqrt();
                        if d > far.0 {
                            far = (d, *p);
                        }
                    }
                    let th = ((far.1.1 - mid.1) as f64).atan2((far.1.0 - mid.0) as f64) as f32;
                    out.push(Spur {
                        poly: pi,
                        i,
                        j,
                        root_km: arc[i] * km_per_cell,
                        len_km: (arc[j] - arc[i]) * km_per_cell,
                        mid,
                        axis_rad: th,
                    });
                    i = j;
                }
                None => i += 1,
            }
        }
    }
    (out, coast_cells)
}

#[must_use]
pub fn coast_shape_thresholds(
    polys: &[Polyline],
    km_per_cell: f32,
    min_spur_km: f32,
    neck_km: f32,
) -> CoastShape {
    // ONE walk, shared with `coast_spurs`, so the spectrum and the count are the same
    // population by construction (ADR Finding 76).
    let (spurs, coast_cells) = coast_spurs(polys, km_per_cell, min_spur_km, neck_km);
    let mut lens: Vec<f32> = spurs.iter().map(|s| s.len_km).collect();
    let roots: Vec<f32> = spurs.iter().map(|s| s.root_km).collect();
    let axes: Vec<(f32, f32, f32)> = spurs.iter().map(|s| (s.mid.0, s.mid.1, s.axis_rad)).collect();
    let in_spur_cells: f32 = spurs.iter().map(|s| s.len_km / km_per_cell).sum();
    let (mut c2, mut s2) = (0.0f64, 0.0f64);
    for sp in &spurs {
        c2 += (2.0 * sp.axis_rad as f64).cos();
        s2 += (2.0 * sp.axis_rad as f64).sin();
    }
    lens.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let mut gaps: Vec<f32> = Vec::new();
    for w in roots.windows(2) {
        if w[1] > w[0] {
            gaps.push(w[1] - w[0]);
        }
    }
    let spacing_cv = if gaps.is_empty() {
        0.0
    } else {
        let m = gaps.iter().sum::<f32>() / gaps.len() as f32;
        let var = gaps.iter().map(|g| (g - m) * (g - m)).sum::<f32>() / gaps.len() as f32;
        if m > 0.0 { var.sqrt() / m } else { 0.0 }
    };

    // LOCAL parallelism: each spur against its nearest neighbours by midpoint.
    let mut local: Vec<f32> = Vec::new();
    for (i, a) in axes.iter().enumerate() {
        if axes.len() <= LOCAL_K {
            break;
        }
        let mut d: Vec<(f32, usize)> = axes
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .map(|(j, b)| (((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt(), j))
            .collect();
        d.select_nth_unstable_by(LOCAL_K - 1, |x, y| {
            x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal)
        });
        let (mut cs, mut sn) = ((2.0 * a.2).cos(), (2.0 * a.2).sin());
        for (_, j) in d.iter().take(LOCAL_K) {
            cs += (2.0 * axes[*j].2).cos();
            sn += (2.0 * axes[*j].2).sin();
        }
        let n = (LOCAL_K + 1) as f32;
        local.push(((cs / n).powi(2) + (sn / n).powi(2)).sqrt());
    }
    local.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let coast_km = coast_cells * km_per_cell;
    let naxis = axes.len();
    CoastShape {
        coast_km,
        count: lens.len(),
        density_per_100km: 100.0 * lens.len() as f32 / coast_km.max(1e-6),
        med_len_km: pct(&lens, 0.5),
        p90_len_km: pct(&lens, 0.9),
        max_len_km: lens.last().copied().unwrap_or(0.0),
        share_pct: 100.0 * in_spur_cells / coast_cells.max(1e-6),
        spacing_cv,
        global_axis_r: if naxis == 0 {
            0.0
        } else {
            (((c2 / naxis as f64).powi(2) + (s2 / naxis as f64).powi(2)).sqrt()) as f32
        },
        local_axis_r: if local.is_empty() { 0.0 } else { local[local.len() / 2] },
    }
}

/// Land/sea render with the contour on top — sea 0, land 0.65, coastline 1.0.
///
/// ADR rule 9: one panel per setting on any visual sweep, FIXED crop, plus a reference panel of
/// what "absent" looks like. Three wrong attributions and two remedies evaluated against blind
/// metrics were each caught by a picture and none by the summary statistics.
pub fn render_coast_crop(
    field: &GridF32,
    polys: &[Polyline],
    sea: f32,
    ox: usize,
    oy: usize,
    crop: usize,
    path: &std::path::Path,
) -> Result<(), String> {
    let mut img = GridF32::new(crop, crop, 0.0);
    for y in 0..crop {
        for x in 0..crop {
            let (sx, sy) = (ox + x, oy + y);
            if sx < field.width && sy < field.height && field.data[sy * field.width + sx] > sea {
                img.set(x, y, 0.65);
            }
        }
    }
    for pl in polys {
        for &(px, py) in pl.iter() {
            let (ix, iy) = (px.round() as isize - ox as isize, py.round() as isize - oy as isize);
            if ix >= 0 && iy >= 0 && (ix as usize) < crop && (iy as usize) < crop {
                img.set(ix as usize, iy as usize, 1.0);
            }
        }
    }
    img.save_png_u8(path)
}

/// The window holding the most coastline vertices — choose ONCE, then hold fixed across panels.
#[must_use]
pub fn densest_window(polys: &[Polyline], w: usize, h: usize, crop: usize) -> (usize, usize) {
    let step = crop / 2;
    let (mut best, mut bxy) = (0usize, (0usize, 0usize));
    let mut oy = 0;
    while oy + crop <= h {
        let mut ox = 0;
        while ox + crop <= w {
            let mut n = 0usize;
            for pl in polys {
                for &(px, py) in pl.iter() {
                    let (x, y) = (px as usize, py as usize);
                    if x >= ox && x < ox + crop && y >= oy && y < oy + crop {
                        n += 1;
                    }
                }
            }
            if n > best {
                best = n;
                bxy = (ox, oy);
            }
            ox += step;
        }
        oy += step;
    }
    bxy
}

#[cfg(test)]
mod tests {
    use super::*;

    /// ADR Finding 75 — the parameterised detector must BE the shipped one at the shipped
    /// thresholds. Every coastal figure in the dossier was measured through `coast_shape`; a
    /// refactor that moved them silently would invalidate the lot, so the identity is pinned
    /// rather than assumed. Rule 10: the fixture must actually produce spurs, or the comparison
    /// is two zeroes agreeing.
    #[test]
    fn coast_shape_default_is_the_parameterised_call() {
        // a ragged closed loop: a circle with a radial spur every few samples
        let mut pl: Polyline = Vec::new();
        for i in 0..720 {
            let t = i as f32 * std::f32::consts::TAU / 720.0;
            let r = 120.0 + if i % 17 < 4 { 28.0 } else { 0.0 };
            pl.push((256.0 + r * t.cos(), 256.0 + r * t.sin()));
        }
        pl.push(pl[0]);
        let polys = vec![pl];
        let km_per_cell = 400.0 / 8192.0;
        let a = coast_shape(&polys, km_per_cell);
        let b = coast_shape_thresholds(&polys, km_per_cell, MIN_SPUR_KM, NECK_KM);
        assert!(a.count > 0, "population guard: the fixture must produce spurs, got {}", a.count);
        assert_eq!(a.count, b.count);
        assert_eq!(a.coast_km.to_bits(), b.coast_km.to_bits(), "coast_km moved");
        assert_eq!(a.p90_len_km.to_bits(), b.p90_len_km.to_bits(), "p90 moved");
        assert_eq!(a.local_axis_r.to_bits(), b.local_axis_r.to_bits(), "local R moved");
        // NEGATIVE CONTROL: the thresholds must matter, or the identity above is vacuous.
        let cell = coast_shape_thresholds(&polys, km_per_cell, 2.0 * km_per_cell, km_per_cell);
        assert_ne!(cell.count, a.count, "the cell-scale setting changed nothing — not a detector");
    }
}
