//! #165 C1 climate — Whittaker biome classification from temperature +
//! precipitation. The READABLE judgement instrument: a biome ("steppe",
//! "temperate forest") is directly judgeable where "600 mm/yr" is not.
//!
//! Each cell is classified by (mean annual temperature °C, annual precipitation
//! mm/yr) per the Whittaker diagram. Thresholds anchored on the real Whittaker
//! biome boundaries, NOT arbitrary. Derived from (T, P) → re-runnable like
//! `c1_climate`: if the precipitation changes, the biomes recompute.

use crate::grid::GridF32;

use super::precipitation::{precip_mm_per_year, SEA_LEVEL_NORM};

/// Whittaker biomes (+ Ocean for sub-sea cells).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Biome {
    Ocean,
    Tundra,
    BorealForest,        // taiga
    TemperateGrassland,  // steppe / prairie
    TemperateForest,
    TemperateRainforest,
    Desert,
    Savanna,
    TropicalSeasonalForest,
    TropicalRainforest,
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
        }
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

/// Biome map for a relief heightmap given its climate (temperature °C +
/// precipitation in INTERNAL units). Sub-sea cells → `Ocean`. Derived & gated:
/// recompute when T/P change. Row-major `Vec<Biome>` (length nx·ny).
pub fn compute_biomes(heightmap: &GridF32, temperature: &GridF32, precipitation: &GridF32) -> Vec<Biome> {
    let n = heightmap.width * heightmap.height;
    (0..n)
        .map(|k| {
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
