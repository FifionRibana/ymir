//! ADR Finding 152-B8 — the two views' colours, shared by the viz and the benches (so a bench's capture is the view).
//! Both write into a south-first RGB(A) buffer (y = 0 = south), like every layer of the viz before its north-up flip.

use super::rocks::ROCK_CLASSES;
use super::zoning::Zoning;

/// The « Roches » view: one colour per class (water keeps the base colour).
pub fn rocks_rgba(rocks: &[u8], base: &mut [u8]) {
    for (k, &r) in rocks.iter().enumerate() {
        if r == 0 {
            continue;
        }
        let c = ROCK_CLASSES[r as usize].color;
        base[k * 4] = c[0];
        base[k * 4 + 1] = c[1];
        base[k * 4 + 2] = c[2];
    }
}

/// The « Chances » view: at each HD cell, the resource `only` (or, with `None`, the resource of highest
/// favourability; ties: the rules file's order) blended over the base with opacity = favourability %.
pub fn chances_rgba(z: &Zoning, only: Option<usize>, w: usize, h: usize, base: &mut [u8]) {
    let f = (w / z.w4.max(1)).max(1);
    for y in 0..h {
        let y4 = (y / f).min(z.h4 - 1);
        for x in 0..w {
            let q = y4 * z.w4 + (x / f).min(z.w4 - 1);
            let pick = match only {
                Some(i) => z.resources.get(i).map(|r| (r, r.fav[q])),
                None => z.resources.iter().map(|r| (r, r.fav[q])).fold(None, |acc: Option<(&_, u8)>, cur| match acc {
                    Some(a) if a.1 >= cur.1 => Some(a),
                    _ => Some(cur),
                }),
            };
            let Some((r, v)) = pick else { continue };
            if v == 0 {
                continue;
            }
            let a = v as f32 / 100.0;
            let k = (y * w + x) * 4;
            for c in 0..3 {
                base[k + c] = (base[k + c] as f32 * (1.0 - a) + r.color[c] as f32 * a) as u8;
            }
        }
    }
}
