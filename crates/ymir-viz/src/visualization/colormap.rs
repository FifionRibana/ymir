//! Colormaps for crustal thickness visualization.

/// Hypsometric colormap for crustal thickness.
/// `t` is normalized in [0, 1]. Returns [r, g, b, a] in u8.
pub fn hypsometric_colormap(t: f64) -> [u8; 4] {
    const STOPS: &[(f64, [u8; 3])] = &[
        (0.0, [20, 50, 120]),   // deep blue — thin oceanic
        (0.2, [40, 120, 160]),  // blue-green — thick oceanic
        (0.4, [60, 160, 80]),   // green — continental plains
        (0.6, [180, 140, 60]),  // light brown — hills
        (0.8, [140, 80, 40]),   // dark brown — mountains
        (1.0, [240, 240, 240]), // white — max thickening
    ];

    let t = t.clamp(0.0, 1.0);

    // Find the segment
    let mut i = 0;
    while i < STOPS.len() - 2 && t > STOPS[i + 1].0 {
        i += 1;
    }

    let (t0, c0) = STOPS[i];
    let (t1, c1) = STOPS[i + 1];
    let frac = ((t - t0) / (t1 - t0)).clamp(0.0, 1.0);

    let r = (c0[0] as f64 + frac * (c1[0] as f64 - c0[0] as f64)) as u8;
    let g = (c0[1] as f64 + frac * (c1[1] as f64 - c0[1] as f64)) as u8;
    let b = (c0[2] as f64 + frac * (c1[2] as f64 - c0[2] as f64)) as u8;

    [r, g, b, 255]
}

/// Slope heatmap: red (steep) to green-yellow (flat).
/// `gx`, `gy` are gradient components from `GridF32::gradient_at`.
pub fn slope_color(gx: f32, gy: f32) -> [u8; 4] {
    let mag = (gx * gx + gy * gy).sqrt();
    let t = (mag / 0.5).clamp(0.0, 1.0);
    [(t * 220.0) as u8, ((1.0 - t) * 180.0 + 40.0) as u8, 30, 255]
}

// ── Step 8.6 v2 colormaps ──────────────────────────────────────────────
//
// These colormaps consume a normalized parameter `t ∈ [0, 1]` produced by
// the field-specific normaliser in `visualization::v2_viz`. Linear scales
// (S̃, age, cratonic) take `t = (v − vmin) / (vmax − vmin)`; log scales
// (ε̇_II, |v|) take `t = (log(v) − log(lo)) / (log(hi) − log(lo))` with a
// cosmetic floor so `log(0)` cannot fire.

fn lerp_stops(stops: &[(f64, [u8; 3])], t: f64) -> [u8; 4] {
    let t = t.clamp(0.0, 1.0);
    let mut i = 0;
    while i < stops.len() - 2 && t > stops[i + 1].0 {
        i += 1;
    }
    let (t0, c0) = stops[i];
    let (t1, c1) = stops[i + 1];
    let frac = ((t - t0) / (t1 - t0)).clamp(0.0, 1.0);
    let r = (c0[0] as f64 + frac * (c1[0] as f64 - c0[0] as f64)) as u8;
    let g = (c0[1] as f64 + frac * (c1[1] as f64 - c0[1] as f64)) as u8;
    let b = (c0[2] as f64 + frac * (c1[2] as f64 - c0[2] as f64)) as u8;
    [r, g, b, 255]
}

/// Age colormap. Young (`t = 0`) renders teal, intermediate ages run
/// through olive / brown, old age (`t = 1`) saturates near orange-red.
/// Choice: D5 says "younger to older gradient, e.g. greens to browns".
pub fn age_colormap(t: f64) -> [u8; 4] {
    const STOPS: &[(f64, [u8; 3])] = &[
        (0.00, [40, 100, 90]),  // teal — freshly reset
        (0.25, [120, 150, 60]), // olive
        (0.50, [180, 150, 70]), // tan
        (0.75, [180, 100, 50]), // ochre
        (1.00, [200, 60, 30]),  // burnt-orange — oldest
    ];
    lerp_stops(STOPS, t)
}

/// Cratonic factor. `t = 0` → black (mobile), `t = 1` → white (cratonic
/// core). Linear grayscale; matches the §9 immunity intuition.
pub fn cratonic_grayscale(t: f64) -> [u8; 4] {
    let v = (t.clamp(0.0, 1.0) * 255.0) as u8;
    [v, v, v, 255]
}

/// Bipolar altitude colormap — replicates the local `hypsometric`
/// helper used by the C1 Track A / B / D visual galleries
/// (`c1_phase_2_*_visual_gallery.rs::hypsometric(h, sea_norm)`)
/// verbatim. Single source of truth for the C1 live-altitude
/// render path (Issue #137 Stage E4) — calling this function with
/// the same `(t, sea_norm)` as the gallery produces byte-identical
/// RGB.
///
/// `t ∈ [0, 1]` is the **normalised** altitude (callers apply
/// `t = (alt + half_range) / (2 * half_range)` for the bipolar
/// `[-half_range, +half_range]` altitude range; the gallery uses
/// `half_range = 1.13`). `sea_norm ∈ [0, 1]` is the position of
/// sea level in the normalised axis (gallery uses `0.5`, so
/// `alt = 0` sits on the blue→green transition).
///
/// Four-segment ramp:
///   - `[0, sea_norm/2)`         deep ocean — `[10, 20, 60]` → `[40, 80, 160]`
///   - `[sea_norm/2, sea_norm)`  shallow ocean — `[40, 80, 160]` → `[120, 180, 230]`
///   - `[sea_norm, mid)`         low land — `[60, 130, 60]` → `[140, 100, 50]`
///   - `[mid, 1]`                high land — `[140, 100, 50]` → `[245, 245, 245]`
///
/// where `mid = (sea_norm + 1.0) / 2.0` (= 0.75 at `sea_norm = 0.5`).
pub fn hypsometric_bipolar(t: f32, sea_norm: f32) -> [u8; 4] {
    let mid = (sea_norm + 1.0) * 0.5;
    let lerp = |t: f32, a: [u8; 3], b: [u8; 3]| -> [u8; 3] {
        let t = t.clamp(0.0, 1.0);
        [
            (a[0] as f32 + t * (b[0] as f32 - a[0] as f32)).round() as u8,
            (a[1] as f32 + t * (b[1] as f32 - a[1] as f32)).round() as u8,
            (a[2] as f32 + t * (b[2] as f32 - a[2] as f32)).round() as u8,
        ]
    };
    let rgb = if t <= sea_norm * 0.5 {
        let tt = t / (sea_norm * 0.5).max(1e-6);
        lerp(tt, [10, 20, 60], [40, 80, 160])
    } else if t <= sea_norm {
        let tt = (t - sea_norm * 0.5) / (sea_norm * 0.5).max(1e-6);
        lerp(tt, [40, 80, 160], [120, 180, 230])
    } else if t <= mid {
        let tt = (t - sea_norm) / (mid - sea_norm).max(1e-6);
        lerp(tt, [60, 130, 60], [140, 100, 50])
    } else {
        let tt = (t - mid) / (1.0 - mid).max(1e-6);
        lerp(tt, [140, 100, 50], [245, 245, 245])
    };
    [rgb[0], rgb[1], rgb[2], 255]
}

/// Inferno-like ramp for log-scaled fields (ε̇_II, |v|). Deep purple →
/// red → orange → yellow → near-white. `t ∈ [0, 1]` already sits on
/// the log axis; this function does no additional scaling.
pub fn log_hot(t: f64) -> [u8; 4] {
    const STOPS: &[(f64, [u8; 3])] = &[
        (0.00, [10, 0, 30]),     // dark purple — quiescent
        (0.25, [80, 10, 80]),    // violet
        (0.50, [200, 40, 40]),   // red
        (0.75, [240, 150, 30]),  // orange
        (1.00, [250, 250, 200]), // pale yellow — saturated yielding / fast flow
    ];
    lerp_stops(STOPS, t)
}

/// Map a value to `t ∈ [0, 1]` on a logarithmic axis bounded
/// `[lo, hi]`. `lo` is a cosmetic floor (e.g. `1e-3` for ε̇_II,
/// `1e-5` for |v|) so `log(0)` can never fire.
pub fn log_normalize(value: f64, lo: f64, hi: f64) -> f64 {
    let v = value.max(lo);
    let log_lo = lo.ln();
    let log_hi = hi.ln();
    let span = (log_hi - log_lo).max(1e-12);
    ((v.ln() - log_lo) / span).clamp(0.0, 1.0)
}
