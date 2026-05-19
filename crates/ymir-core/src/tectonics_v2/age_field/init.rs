//! Step 10 — initial `A` field construction from the initial `S̃`
//! field. The construction uses a static threshold (`S̃ > 0.5`) to
//! classify continental vs oceanic per D2 / D7; subsequent
//! transitions between continental and oceanic during the run
//! happen via S̃ advection + recycling and do NOT trigger
//! re-initialisation of A. The `A` field at later times reflects
//! the cell's advected age modulated by boundary-event resets.
//!
//! The construction lives mostly in `AgeFieldState::from_initial_thickness`
//! (parent module). This file collects classification helpers
//! used by the metrics layer.

use crate::tectonics_v2::field::Field2D;

/// `S̃ > 0.5` classification used at init **and** at run-time
/// metrics. Returns true for continental, false for oceanic.
/// Centralised here so the same threshold is enforced wherever
/// the age-field pipeline classifies a cell.
#[inline]
pub fn is_continental_thickness(s: f64) -> bool {
    s > 0.5
}

/// Count of continental cells in `s` per the [`is_continental_thickness`]
/// threshold. Used by the metrics layer for `age_at_continental_cells_mean`.
pub fn count_continental(s: &Field2D) -> usize {
    s.data().iter().filter(|&&v| is_continental_thickness(v)).count()
}

/// Count of oceanic cells in `s`.
pub fn count_oceanic(s: &Field2D) -> usize {
    s.data().iter().filter(|&&v| !is_continental_thickness(v)).count()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classification_threshold_is_half() {
        assert!(!is_continental_thickness(0.5)); // boundary excluded
        assert!(!is_continental_thickness(0.4));
        assert!(is_continental_thickness(0.51));
        assert!(is_continental_thickness(1.0));
    }

    #[test]
    fn counts_partition_the_grid() {
        let mut s = Field2D::filled(4, 4, 0.2);
        s.set(0, 0, 1.0);
        s.set(1, 1, 1.0);
        assert_eq!(count_continental(&s), 2);
        assert_eq!(count_oceanic(&s), 14);
        assert_eq!(count_continental(&s) + count_oceanic(&s), 16);
    }
}
