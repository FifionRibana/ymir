//! #165 C1 climate — Whittaker biome classification from temperature +
//! precipitation. The READABLE judgement instrument: a biome ("steppe",
//! "temperate forest") is directly judgeable where "600 mm/yr" is not.
//!
//! Each cell is classified by (mean annual temperature °C, annual precipitation
//! mm/yr) per the Whittaker diagram. Thresholds anchored on the real Whittaker
//! biome boundaries, NOT arbitrary. Derived from (T, P) → re-runnable like
//! `c1_climate`: if the precipitation changes, the biomes recompute.
//!
//! # Wire contract — `ymir.WhittakerBiome@v1`
//!
//! [`Biome::to_u8`] / [`Biome::from_u8`] define the STABLE on-disk id used by
//! the `.ymir` `biome.u8` raster (see [`crate::export::container`]). `0` is a
//! reserved sentinel (Undefined, no variant); the ten Whittaker biomes are
//! `1..=10` in the order below. These numbers are FROZEN — a consumer decodes
//! against them, so a variant reorder must NOT change them (a test pins each).

/// Reserved sentinel id in the `ymir.WhittakerBiome@v1` contract (no variant).
pub const BIOME_UNDEFINED_U8: u8 = 0;

use crate::grid::GridF32;

use super::precipitation::{SEA_LEVEL_NORM, precip_mm_per_year};

/// Whittaker biomes (+ Ocean for sub-sea cells).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Biome {
    Ocean,
    Tundra,
    BorealForest,       // taiga
    TemperateGrassland, // steppe / prairie
    TemperateForest,
    TemperateRainforest,
    Desert,
    Savanna,
    TropicalSeasonalForest,
    TropicalRainforest,
    /// Inland water body (lake) — distinct from `Ocean`. Assigned from `water_class`
    /// (enclosed water) / the drainage lake map, NOT from altitude, so a below-sea
    /// enclosed basin's water surface reads as a lake, not the sea (ADR 0001 Finding 18).
    Lake,
}

impl Biome {
    /// Display name.
    pub fn name(self) -> &'static str {
        match self {
            Biome::Ocean => "ocean",
            Biome::Tundra => "tundra",
            Biome::BorealForest => "boreal forest (taiga)",
            Biome::TemperateGrassland => "temperate grassland (steppe)",
            Biome::TemperateForest => "temperate forest",
            Biome::TemperateRainforest => "temperate rainforest",
            Biome::Desert => "desert",
            Biome::Savanna => "savanna",
            Biome::TropicalSeasonalForest => "tropical seasonal forest",
            Biome::TropicalRainforest => "tropical rainforest",
            Biome::Lake => "lake (inland water)",
        }
    }

    /// Stable `ymir.WhittakerBiome@v1` id (`1..=10`; see module header). FROZEN
    /// wire format — never renumber. `0` is reserved ([`BIOME_UNDEFINED_U8`]).
    pub fn to_u8(self) -> u8 {
        match self {
            Biome::Ocean => 1,
            Biome::Tundra => 2,
            Biome::BorealForest => 3,
            Biome::TemperateGrassland => 4,
            Biome::TemperateForest => 5,
            Biome::TemperateRainforest => 6,
            Biome::Desert => 7,
            Biome::Savanna => 8,
            Biome::TropicalSeasonalForest => 9,
            Biome::TropicalRainforest => 10,
            Biome::Lake => 11,
        }
    }

    /// Inverse of [`Biome::to_u8`]. `None` for the reserved `0` sentinel or any
    /// unknown id (a reader must tolerate ids it does not recognise).
    pub fn from_u8(id: u8) -> Option<Biome> {
        Some(match id {
            1 => Biome::Ocean,
            2 => Biome::Tundra,
            3 => Biome::BorealForest,
            4 => Biome::TemperateGrassland,
            5 => Biome::TemperateForest,
            6 => Biome::TemperateRainforest,
            7 => Biome::Desert,
            8 => Biome::Savanna,
            9 => Biome::TropicalSeasonalForest,
            10 => Biome::TropicalRainforest,
            11 => Biome::Lake,
            _ => return None,
        })
    }

    /// Categorical render colour (RGB) — distinct per biome for a readable map.
    pub fn color(self) -> [u8; 3] {
        match self {
            Biome::Ocean => [30, 50, 90],
            Biome::Tundra => [200, 205, 215],
            Biome::BorealForest => [70, 110, 90],
            Biome::TemperateGrassland => [200, 195, 110],
            Biome::TemperateForest => [80, 160, 80],
            Biome::TemperateRainforest => [40, 110, 70],
            Biome::Desert => [225, 200, 140],
            Biome::Savanna => [190, 175, 90],
            Biome::TropicalSeasonalForest => [120, 175, 70],
            Biome::TropicalRainforest => [20, 110, 50],
            Biome::Lake => [60, 110, 175],
        }
    }
}

/// Classify a LAND cell by Whittaker (mean annual T in °C, annual P in mm/yr).
/// Thresholds anchored on the Whittaker diagram: cold → tundra/taiga; the 250 mm
/// arid line → desert; temperate split grassland/forest/rainforest at ~500/1500
/// mm; warm split savanna/seasonal/rainforest at ~1000/2000 mm.
pub fn classify(t_c: f32, p_mm: f32) -> Biome {
    if t_c < -5.0 {
        Biome::Tundra
    } else if p_mm < 250.0 {
        Biome::Desert
    } else if t_c < 5.0 {
        Biome::BorealForest
    } else if t_c < 20.0 {
        if p_mm < 500.0 {
            Biome::TemperateGrassland
        } else if p_mm < 1500.0 {
            Biome::TemperateForest
        } else {
            Biome::TemperateRainforest
        }
    } else if p_mm < 1000.0 {
        Biome::Savanna
    } else if p_mm < 2000.0 {
        Biome::TropicalSeasonalForest
    } else {
        Biome::TropicalRainforest
    }
}

/// Biome map from CONNECTIVITY + drainage water, not altitude (ADR 0001 Finding 18).
/// `water_class` (0 land / 1 ocean / 2 inland below-sea) decides sea-vs-inland membership;
/// `lake_map` (per-cell lake id, 0 = none) marks the drainage water surface (real lakes AND
/// the water-balanced below-sea basins at their computed level). Per cell:
///   - `water_class == 1` (edge-connected below-sea) → `Ocean`;
///   - `lake_map != 0` (under an inland water surface) → `Lake`;
///   - otherwise → Whittaker land — EVEN where altitude < 0 (exposed dry land below the sea
///     level around an endorheic basin: the level came from `min(spill, evaporative)`, so we
///     do NOT flood it to 0 m).
/// Pass empty `water_class`/`lake_map` slices to fall back to the legacy altitude rule.
pub fn compute_biomes(
    heightmap: &GridF32,
    temperature: &GridF32,
    precipitation: &GridF32,
    water_class: &[u8],
    lake_map: &[u32],
) -> Vec<Biome> {
    let n = heightmap.width * heightmap.height;
    let has_wc = water_class.len() == n;
    let has_lakes = lake_map.len() == n;
    (0..n)
        .map(|k| {
            if has_wc {
                if water_class[k] == 1 {
                    return Biome::Ocean; // edge-connected sea
                }
                if has_lakes && lake_map[k] != 0 {
                    return Biome::Lake; // inland water surface (incl. below-sea basins)
                }
                // land (class 0) OR exposed below-sea land (class 2 above the lake level)
                return classify(temperature.data[k], precip_mm_per_year(precipitation.data[k]));
            }
            // Legacy fallback (no connectivity info): altitude membership.
            if heightmap.data[k] <= SEA_LEVEL_NORM {
                Biome::Ocean
            } else {
                classify(temperature.data[k], precip_mm_per_year(precipitation.data[k]))
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `ymir.WhittakerBiome@v1` — the on-disk ids are FROZEN. This pins every
    /// number so a future variant reorder can't silently change the wire format,
    /// and checks the `from_u8` round-trip + the reserved `0` sentinel.
    #[test]
    fn biome_u8_ids_are_frozen() {
        let table = [
            (Biome::Ocean, 1u8),
            (Biome::Tundra, 2),
            (Biome::BorealForest, 3),
            (Biome::TemperateGrassland, 4),
            (Biome::TemperateForest, 5),
            (Biome::TemperateRainforest, 6),
            (Biome::Desert, 7),
            (Biome::Savanna, 8),
            (Biome::TropicalSeasonalForest, 9),
            (Biome::TropicalRainforest, 10),
            (Biome::Lake, 11),
        ];
        for (biome, id) in table {
            assert_eq!(biome.to_u8(), id, "{biome:?} id must stay {id}");
            assert_eq!(Biome::from_u8(id), Some(biome), "id {id} must decode to {biome:?}");
        }
        assert_eq!(Biome::from_u8(BIOME_UNDEFINED_U8), None, "0 is the reserved sentinel");
        assert_eq!(Biome::from_u8(12), None, "unknown ids decode to None");
    }

    /// The full climate contract: a hand-built (T, P) classifies to the biome
    /// whose FROZEN id the `biome.u8` raster would carry.
    #[test]
    fn classify_to_u8_matches_expected_id() {
        assert_eq!(classify(12.0, 900.0).to_u8(), 5, "temperate forest → 5");
        assert_eq!(classify(25.0, 100.0).to_u8(), 7, "hot+dry desert → 7");
        assert_eq!(classify(-12.0, 600.0).to_u8(), 2, "cold tundra → 2");
        assert_eq!(classify(26.0, 2500.0).to_u8(), 10, "tropical rainforest → 10");
    }

    #[test]
    fn whittaker_known_cases() {
        // temperate + moderate rain → temperate forest.
        assert_eq!(classify(12.0, 900.0), Biome::TemperateForest);
        // temperate + dry → steppe.
        assert_eq!(classify(12.0, 350.0), Biome::TemperateGrassland);
        // hot + dry → desert.
        assert_eq!(classify(25.0, 100.0), Biome::Desert);
        // cold → tundra (regardless of P).
        assert_eq!(classify(-12.0, 600.0), Biome::Tundra);
        // cool + wet → boreal forest.
        assert_eq!(classify(2.0, 600.0), Biome::BorealForest);
        // hot + very wet → tropical rainforest.
        assert_eq!(classify(26.0, 2500.0), Biome::TropicalRainforest);
        // temperate + very wet → temperate rainforest.
        assert_eq!(classify(10.0, 1800.0), Biome::TemperateRainforest);
    }

    /// #165 — biome↔climate COHERENCE invariant. Locks the manual cross-check
    /// from the final climate re-validation ("temperate forest must NOT sit on
    /// dry/cold ground") into an automatic guard. EXHAUSTIVE over the (T, P)
    /// plane: every cell `classify` calls a temperate FOREST must lie in the
    /// temperate thermal niche [+5, +20) °C AND be humid enough — anchored on
    /// the SAME Whittaker thresholds `classify` uses (500 mm forest floor,
    /// 1500 mm rainforest floor). By contraposition this forbids temperate
    /// forest on steppe-dry (< 500 mm), on desert (< 250), on taiga-cold
    /// (< +5 °C), or on warm/tropical (≥ +20 °C). A future change that placed a
    /// temperate forest off this niche (the false-alarm the re-validation
    /// caught by eye) breaks this test.
    #[test]
    fn temperate_forest_implies_temperate_and_humid() {
        let mut t = -40.0f32;
        while t <= 50.0 {
            let mut p = 0.0f32;
            while p <= 5000.0 {
                match classify(t, p) {
                    Biome::TemperateForest => {
                        assert!(
                            (5.0..20.0).contains(&t),
                            "temperate forest off the thermal niche: t={t} °C (p={p} mm)"
                        );
                        assert!(
                            (500.0..1500.0).contains(&p),
                            "temperate forest outside its humid band: p={p} mm (t={t} °C)"
                        );
                    }
                    Biome::TemperateRainforest => {
                        assert!(
                            (5.0..20.0).contains(&t),
                            "temperate rainforest off the thermal niche: t={t} °C (p={p} mm)"
                        );
                        assert!(
                            p >= 1500.0,
                            "temperate rainforest below its humid floor: p={p} mm (t={t} °C)"
                        );
                    }
                    _ => {}
                }
                p += 5.0;
            }
            t += 0.25;
        }
    }
}
