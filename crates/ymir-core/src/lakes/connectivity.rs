//! Water-body connectivity: split below-sea cells into edge-connected OCEAN vs
//! enclosed INLAND basins via a border flood-fill.
//!
//! This is pure topology over an altitude field at a given sea level — it does
//! NOT re-run drainage. It answers a question the coastline isoline alone
//! cannot: a below-sea cell surrounded by land (an endorheic depression, a
//! below-sea inland sea) is NOT the ocean. A consumer needs that split to shade
//! ocean vs inland water differently.
//!
//! Encoding (see [`WATER_CLASS_LAND`]/`_OCEAN`/`_INLAND`):
//! - `0` land — altitude `> sea_level`.
//! - `1` ocean — altitude `<= sea_level` AND reachable from a grid edge through
//!   other below-sea cells (8-connected — a DIAGONAL contact is a real hydrological
//!   connection, ADR 0001 Finding 35; a coastal pocket touching the sea only at a
//!   corner IS sea, and 4-connectivity mis-classed it as enclosed inland).
//! - `2` inland — altitude `<= sea_level` but NOT edge-reachable (enclosed).
//!
//! The caller passes the SAME `sea_level` the coastline isoline is traced at, so
//! the two layers agree by construction (no second sea-level constant here).

use std::collections::VecDeque;

use crate::grid::GridF32;

/// Above sea level.
pub const WATER_CLASS_LAND: u8 = 0;
/// Below sea level and connected to a grid edge (the open ocean).
pub const WATER_CLASS_OCEAN: u8 = 1;
/// Below sea level but enclosed by land (an inland below-sea basin).
pub const WATER_CLASS_INLAND: u8 = 2;

/// Classify every cell as land / ocean / inland at `sea_level` (row-major
/// `Vec<u8>`, length `width*height`). Deterministic: BFS from a fixed edge-cell
/// order over 8-neighbours (Finding 35 — diagonal contact is a water connection);
/// the result depends only on the reachable set.
pub fn water_class(height: &GridF32, sea_level: f32) -> Vec<u8> {
    let w = height.width;
    let h = height.height;
    let n = w * h;
    if n == 0 {
        return Vec::new();
    }

    let below = |k: usize| height.data[k] <= sea_level;
    let mut ocean = vec![false; n]; // edge-reachable below-sea cells
    let mut queue: VecDeque<usize> = VecDeque::new();

    // Seed the BFS from every below-sea cell on the four borders (deterministic
    // order: top row, bottom row, then left/right columns).
    let seed = |k: usize, ocean: &mut Vec<bool>, queue: &mut VecDeque<usize>| {
        if below(k) && !ocean[k] {
            ocean[k] = true;
            queue.push_back(k);
        }
    };
    for x in 0..w {
        seed(x, &mut ocean, &mut queue); // top (y=0)
        seed((h - 1) * w + x, &mut ocean, &mut queue); // bottom
    }
    for y in 0..h {
        seed(y * w, &mut ocean, &mut queue); // left (x=0)
        seed(y * w + (w - 1), &mut ocean, &mut queue); // right
    }

    // Flood over 8-connected below-sea neighbours (diagonal contact = water connection).
    while let Some(k) = queue.pop_front() {
        let (x, y) = ((k % w) as i32, (k / w) as i32);
        for (dx, dy) in [
            (-1i32, 0i32), (1, 0), (0, -1), (0, 1),
            (-1, -1), (-1, 1), (1, -1), (1, 1),
        ] {
            let (nx, ny) = (x + dx, y + dy);
            if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                let nk = ny as usize * w + nx as usize;
                if below(nk) && !ocean[nk] {
                    ocean[nk] = true;
                    queue.push_back(nk);
                }
            }
        }
    }

    (0..n)
        .map(|k| {
            if !below(k) {
                WATER_CLASS_LAND
            } else if ocean[k] {
                WATER_CLASS_OCEAN
            } else {
                WATER_CLASS_INLAND
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Edge sea → ocean(1); an enclosed below-sea pit → inland(2); rim → land(0).
    #[test]
    fn edge_sea_is_ocean_enclosed_pit_is_inland() {
        let (w, h) = (7usize, 7usize);
        let sea = 0.5f32;
        // Border ring below sea (0.3), interior land (0.7), one enclosed pit at
        // the centre (0.2) surrounded entirely by land.
        let mut data = vec![0.7f32; w * h];
        for x in 0..w {
            data[x] = 0.3; // top
            data[(h - 1) * w + x] = 0.3; // bottom
        }
        for y in 0..h {
            data[y * w] = 0.3; // left
            data[y * w + (w - 1)] = 0.3; // right
        }
        data[3 * w + 3] = 0.2; // enclosed pit at centre
        let grid = GridF32::from_vec(w, h, data);

        let class = water_class(&grid, sea);
        assert_eq!(class[0], WATER_CLASS_OCEAN, "corner border is edge-connected ocean");
        assert_eq!(class[3], WATER_CLASS_OCEAN, "top-edge cell is ocean");
        assert_eq!(class[3 * w + 3], WATER_CLASS_INLAND, "enclosed pit is inland");
        assert_eq!(class[w + 1], WATER_CLASS_LAND, "interior high ground is land");
        // Determinism.
        assert_eq!(class, water_class(&grid, sea));
    }

    /// Finding 35 — a below-sea pocket touching the edge sea ONLY at a DIAGONAL corner is OCEAN
    /// (8-connectivity), not enclosed inland. Under 4-connectivity it was mis-classed as inland,
    /// which starved the priority-flood of a seed and produced aberrant sills.
    #[test]
    fn diagonal_contact_is_ocean_not_inland() {
        let (w, h) = (5usize, 5usize);
        let sea = 0.5f32;
        let mut data = vec![0.7f32; w * h]; // land everywhere
        data[0] = 0.3; // top-left corner = edge sea (ocean seed)
        data[w + 1] = 0.3; // (1,1) touches the corner ONLY diagonally
        let grid = GridF32::from_vec(w, h, data);
        let class = water_class(&grid, sea);
        assert_eq!(class[0], WATER_CLASS_OCEAN, "corner is edge sea");
        assert_eq!(class[w + 1], WATER_CLASS_OCEAN, "diagonally-connected pocket IS ocean (8-conn)");
    }
}
