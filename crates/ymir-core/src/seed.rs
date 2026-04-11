//! Deterministic seed management for reproducible world generation.
//!
//! [`WorldSeed`] takes a master `u64` seed and derives independent sub-seeds
//! for each pipeline phase. This ensures two properties:
//!
//! **Reproducibility**: the same master seed always produces the same world.
//!
//! **Phase independence**: changing parameters in one phase (e.g. erosion rate)
//! does not affect the output of other phases (e.g. plate configuration). Each
//! phase gets its own deterministic RNG stream derived from the master seed and
//! the phase name.

use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

/// Master seed for an entire world generation run.
///
/// Derive phase-specific RNGs with [`WorldSeed::rng_for`]. The phase name is
/// hashed together with the master seed to produce a unique sub-seed per phase,
/// so "erosion" and "tectonics" get completely independent random streams.
///
/// # Example
///
/// ```
/// use ymir_core::seed::WorldSeed;
///
/// let world = WorldSeed::new(42);
///
/// // Each call with the same phase name returns an identically-seeded RNG.
/// let mut rng_a = world.rng_for("erosion");
/// let mut rng_b = world.rng_for("erosion");
/// // rng_a and rng_b will produce the same sequence.
///
/// // Different phase names produce different sequences.
/// let mut rng_tecto = world.rng_for("tectonics");
/// // rng_tecto diverges from rng_a immediately.
/// ```
#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct WorldSeed {
    master: u64,
}

impl WorldSeed {
    /// Create a new world seed from a master value.
    pub fn new(master: u64) -> Self {
        Self { master }
    }

    /// The master seed value. Useful for logging and metadata export.
    pub fn master(&self) -> u64 {
        self.master
    }

    /// Derive a deterministic RNG for a named pipeline phase.
    ///
    /// The phase name is hashed together with the master seed using a
    /// deterministic hasher (SipHash via `DefaultHasher`). The resulting u64
    /// is expanded to a 32-byte seed for ChaCha8, which is a fast
    /// cryptographic-quality PRNG that guarantees identical output across
    /// platforms and Rust versions.
    ///
    /// The hash function is deterministic within a Rust version. Since Ymir
    /// pins its toolchain, this is sufficient for reproducibility. If
    /// cross-version reproducibility becomes a requirement, replace
    /// `DefaultHasher` with a fixed hash algorithm.
    pub fn rng_for(&self, phase: &str) -> ChaCha8Rng {
        let sub_seed = self.derive_seed(phase);
        ChaCha8Rng::seed_from_u64(sub_seed)
    }

    /// Derive a u64 sub-seed for a named phase.
    ///
    /// Useful when you need the raw seed value (e.g. for logging or for
    /// passing to a third-party function that takes a u64 seed).
    pub fn derive_seed(&self, phase: &str) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.master.hash(&mut hasher);
        phase.hash(&mut hasher);
        hasher.finish()
    }

    /// Derive a sub-seed with an additional numeric discriminator.
    ///
    /// Useful for generating multiple independent streams within a single
    /// phase, e.g. one per erosion batch:
    /// ```
    /// # use ymir_core::seed::WorldSeed;
    /// let world = WorldSeed::new(42);
    /// let batch_0_rng = world.rng_for_indexed("erosion_batch", 0);
    /// let batch_1_rng = world.rng_for_indexed("erosion_batch", 1);
    /// ```
    pub fn rng_for_indexed(&self, phase: &str, index: u64) -> ChaCha8Rng {
        let mut hasher = DefaultHasher::new();
        self.master.hash(&mut hasher);
        phase.hash(&mut hasher);
        index.hash(&mut hasher);
        ChaCha8Rng::seed_from_u64(hasher.finish())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    #[test]
    fn test_same_seed_same_output() {
        let seed = WorldSeed::new(12345);
        let mut rng_a = seed.rng_for("erosion");
        let mut rng_b = seed.rng_for("erosion");

        let vals_a: Vec<f64> = (0..100).map(|_| rng_a.random()).collect();
        let vals_b: Vec<f64> = (0..100).map(|_| rng_b.random()).collect();

        assert_eq!(vals_a, vals_b, "Same phase + same seed must give identical output");
    }

    #[test]
    fn test_different_phases_diverge() {
        let seed = WorldSeed::new(12345);
        let mut rng_erosion = seed.rng_for("erosion");
        let mut rng_tecto = seed.rng_for("tectonics");

        let val_e: f64 = rng_erosion.random();
        let val_t: f64 = rng_tecto.random();

        assert_ne!(val_e, val_t, "Different phases should produce different streams");
    }

    #[test]
    fn test_different_seeds_diverge() {
        let seed_a = WorldSeed::new(1);
        let seed_b = WorldSeed::new(2);

        let mut rng_a = seed_a.rng_for("erosion");
        let mut rng_b = seed_b.rng_for("erosion");

        let val_a: f64 = rng_a.random();
        let val_b: f64 = rng_b.random();

        assert_ne!(val_a, val_b, "Different master seeds should produce different streams");
    }

    #[test]
    fn test_indexed_determinism() {
        let seed = WorldSeed::new(42);
        let mut rng_a = seed.rng_for_indexed("batch", 7);
        let mut rng_b = seed.rng_for_indexed("batch", 7);

        let vals_a: Vec<f64> = (0..50).map(|_| rng_a.random()).collect();
        let vals_b: Vec<f64> = (0..50).map(|_| rng_b.random()).collect();

        assert_eq!(vals_a, vals_b);
    }

    #[test]
    fn test_indexed_different_indices_diverge() {
        let seed = WorldSeed::new(42);
        let mut rng_0 = seed.rng_for_indexed("batch", 0);
        let mut rng_1 = seed.rng_for_indexed("batch", 1);

        let val_0: f64 = rng_0.random();
        let val_1: f64 = rng_1.random();

        assert_ne!(val_0, val_1);
    }
}
