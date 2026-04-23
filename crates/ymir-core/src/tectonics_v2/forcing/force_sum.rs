//! `ForceSum` — compose multiple `BodyForce` terms.
//!
//! Each term is a `Box<dyn BodyForce>`. `accumulate` walks the
//! vector in order and calls each term's `accumulate` on the shared
//! output. Because `BodyForce::accumulate` is documented to **add**
//! (never overwrite), this produces the exact sum of the terms.
//!
//! # Canonical order (documentation, not yet enforced)
//!
//! Floating-point addition is not associative, so the order of
//! terms matters for bit-identical regression. At Step 2 only one
//! or two terms are ever used, so no canonical order is enforced.
//! The canonical order is documented below; Step 7 ships slab-pull
//! as a separate assembly step (harness adds `SlabPullForce`
//! contributions after `ForceSum` has laid down GPE), so this
//! ordering is not yet materialised in code.
//!
//! ```text
//! GPE -> basal drag -> slab-pull -> mantle flow
//! ```
//!
//! Step 8 (mantle flow) may introduce the typed-tuple refactor
//! mentioned in the original spec, which would also bring
//! slab-pull into this container.

use super::body_force::{BodyForce, SimulationState, VectorField};

pub struct ForceSum {
    terms: Vec<Box<dyn BodyForce>>,
}

impl ForceSum {
    pub fn new() -> Self {
        Self { terms: Vec::new() }
    }

    pub fn push(&mut self, f: Box<dyn BodyForce>) -> &mut Self {
        self.terms.push(f);
        self
    }

    pub fn with(mut self, f: Box<dyn BodyForce>) -> Self {
        self.push(f);
        self
    }

    pub fn len(&self) -> usize { self.terms.len() }
    pub fn is_empty(&self) -> bool { self.terms.is_empty() }

    /// List of term names, in the order they were pushed. Reported
    /// as `term_names` in the diagnostic markdown.
    pub fn term_names(&self) -> Vec<&'static str> {
        self.terms.iter().map(|t| t.name()).collect()
    }
}

impl Default for ForceSum {
    fn default() -> Self { Self::new() }
}

impl BodyForce for ForceSum {
    fn accumulate(&self, state: &SimulationState, out: &mut VectorField) {
        for term in &self.terms {
            term.accumulate(state, out);
        }
    }
    fn name(&self) -> &'static str { "ForceSum" }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::super::field::{Field2D, PeriodicIndex};
    use super::super::sinusoidal::SinusoidalForce;

    fn env(nx: usize, ny: usize) -> (PeriodicIndex, PeriodicIndex, Field2D) {
        (
            PeriodicIndex::new(nx),
            PeriodicIndex::new(ny),
            Field2D::filled(nx, ny, 1.0),
        )
    }

    #[test]
    fn empty_sum_leaves_output_untouched() {
        let nx = 4;
        let ny = 4;
        let (idx_x, idx_y, s) = env(nx, ny);
        let st = SimulationState { nx, ny, dx: 0.25, dy: 0.25, idx_x: &idx_x, idx_y: &idx_y, s: &s };
        let mut fx = Field2D::filled(nx, ny, 3.0);
        let mut fy = Field2D::filled(nx, ny, -1.0);
        ForceSum::new().accumulate(
            &st,
            &mut VectorField { fx: &mut fx, fy: &mut fy },
        );
        for v in fx.data().iter() {
            assert_eq!(*v, 3.0);
        }
        for v in fy.data().iter() {
            assert_eq!(*v, -1.0);
        }
    }

    #[test]
    fn single_term_equals_term_alone() {
        let nx = 8;
        let ny = 8;
        let (idx_x, idx_y, s) = env(nx, ny);
        let st = SimulationState { nx, ny, dx: 0.125, dy: 0.125, idx_x: &idx_x, idx_y: &idx_y, s: &s };
        let sin = SinusoidalForce::new(10.0, 1.0);
        let mut fx_direct = Field2D::new(nx, ny);
        let mut fy_direct = Field2D::new(nx, ny);
        sin.accumulate(&st, &mut VectorField { fx: &mut fx_direct, fy: &mut fy_direct });

        let mut fx_sum = Field2D::new(nx, ny);
        let mut fy_sum = Field2D::new(nx, ny);
        let mut sum = ForceSum::new();
        sum.push(Box::new(sin));
        sum.accumulate(&st, &mut VectorField { fx: &mut fx_sum, fy: &mut fy_sum });

        for k in 0..(nx * ny) {
            assert_eq!(fx_direct.data()[k], fx_sum.data()[k]);
            assert_eq!(fy_direct.data()[k], fy_sum.data()[k]);
        }
    }

    #[test]
    fn two_terms_equal_sum_of_each() {
        let nx = 8;
        let ny = 8;
        let (idx_x, idx_y, s) = env(nx, ny);
        let st = SimulationState { nx, ny, dx: 0.125, dy: 0.125, idx_x: &idx_x, idx_y: &idx_y, s: &s };
        let a = SinusoidalForce::new(3.0, 1.0);
        let b = SinusoidalForce::new(-2.0, 1.0);

        let mut fx_direct = Field2D::new(nx, ny);
        let mut fy_direct = Field2D::new(nx, ny);
        a.accumulate(&st, &mut VectorField { fx: &mut fx_direct, fy: &mut fy_direct });
        b.accumulate(&st, &mut VectorField { fx: &mut fx_direct, fy: &mut fy_direct });

        let mut fx_sum = Field2D::new(nx, ny);
        let mut fy_sum = Field2D::new(nx, ny);
        let sum = ForceSum::new().with(Box::new(a)).with(Box::new(b));
        sum.accumulate(&st, &mut VectorField { fx: &mut fx_sum, fy: &mut fy_sum });

        for k in 0..(nx * ny) {
            assert!((fx_direct.data()[k] - fx_sum.data()[k]).abs() < 1e-14);
            assert!((fy_direct.data()[k] - fy_sum.data()[k]).abs() < 1e-14);
        }
    }

    #[test]
    fn term_names_preserves_insertion_order() {
        let sum = ForceSum::new()
            .with(Box::new(SinusoidalForce::default()))
            .with(Box::new(crate::tectonics_v2::forcing::ZeroForce));
        assert_eq!(sum.term_names(), vec!["SinusoidalForce", "ZeroForce"]);
    }
}
