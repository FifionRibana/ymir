//! Null-space projection for the periodic thin-sheet momentum operator.
//!
//! On a fully periodic torus the discrete operator `A = -∇·(2η ε̇(·))`
//! has a **2-D null space**: constant velocity per component
//! (`(a, 0)` and `(0, b)`) — rigid-body translation of the entire
//! torus, which the momentum balance does not penalise.
//!
//! There is **no pressure null space**: pressure is not an unknown of
//! the thin-sheet formulation.
//!
//! Every preconditioner application must kill the velocity null-space
//! modes BEFORE and AFTER `M⁻¹`. Two array-wide means per application
//! is O(N) and negligible against the stencil cost.
//!
//! Step 8.5b: the mean is computed with [`par_sum`] (chunk-sequential
//! reduction, bit-identical across thread counts) and the subtract
//! step is `par_iter_mut` (cell-local, order-independent).

use rayon::prelude::*;

use super::parallel_reduce::par_sum;

/// Subtract the arithmetic mean of `data` from every element.
pub fn subtract_mean(data: &mut [f64]) {
    if data.is_empty() {
        return;
    }
    let mean = par_sum(data) / data.len() as f64;
    data.par_iter_mut().for_each(|v| *v -= mean);
}

/// Return the arithmetic mean of `data`.
pub fn mean(data: &[f64]) -> f64 {
    if data.is_empty() {
        0.0
    } else {
        par_sum(data) / data.len() as f64
    }
}

/// Project a velocity pair (vx, vy) onto the zero-mean-per-component
/// subspace. The two components are orthogonal in the null-space
/// decomposition and project independently.
pub fn project_velocity(vx: &mut [f64], vy: &mut [f64]) {
    subtract_mean(vx);
    subtract_mean(vy);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn subtract_mean_zeros_the_mean() {
        let mut v = vec![1.0, 2.0, 3.0, 4.0];
        subtract_mean(&mut v);
        let m: f64 = v.iter().sum::<f64>() / v.len() as f64;
        assert!(m.abs() < 1e-14);
    }

    #[test]
    fn subtract_mean_is_idempotent() {
        let mut v = vec![3.7, -1.2, 0.0, 42.0, 100.0];
        subtract_mean(&mut v);
        let first = v.clone();
        subtract_mean(&mut v);
        for (a, b) in first.iter().zip(v.iter()) {
            assert!((a - b).abs() < 1e-14);
        }
    }

    #[test]
    fn project_velocity_zeros_both_components() {
        let mut vx = vec![1.0, 2.0, 3.0, 4.0];
        let mut vy = vec![10.0, 20.0, 30.0, 40.0];
        project_velocity(&mut vx, &mut vy);
        assert!(mean(&vx).abs() < 1e-14);
        assert!(mean(&vy).abs() < 1e-14);
    }

    #[test]
    fn empty_input_is_safe() {
        let mut v: Vec<f64> = vec![];
        subtract_mean(&mut v);
        assert_eq!(mean(&v), 0.0);
    }
}
