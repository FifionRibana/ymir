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
