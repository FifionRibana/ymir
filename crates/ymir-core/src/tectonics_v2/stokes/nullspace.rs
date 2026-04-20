//! Null-space projection for the periodic Stokes operator.
//!
//! On a fully periodic torus the continuous Stokes system has:
//! - a 1D pressure null space (the constant mode) — any p + c is a
//!   solution because only ∇p enters the momentum balance;
//! - a 2D velocity null space (constant per component: (a, 0) and
//!   (0, b)) — net rigid-body translation of the entire torus, which
//!   the incompressibility constraint does not forbid since ∇·const = 0.
//!
//! Every preconditioner application must kill these modes BEFORE and
//! AFTER applying M⁻¹. Killing them only at the end of each outer
//! Krylov step is insufficient because intermediate search directions
//! contaminated by the null space degrade orthogonality.
//!
//! These projectors are used by [`super::precond`] to wrap preconditioner
//! applications, and by [`super::solver`] to clean the final iterates.

/// Subtract the arithmetic mean of `data` from every element.
///
/// For periodic pressure/cell-centered fields, this is the orthogonal
/// projector onto the zero-mean subspace.
pub fn subtract_mean(data: &mut [f64]) {
    if data.is_empty() {
        return;
    }
    let mean = data.iter().sum::<f64>() / data.len() as f64;
    for v in data.iter_mut() {
        *v -= mean;
    }
}

/// Return the arithmetic mean of `data`.
pub fn mean(data: &[f64]) -> f64 {
    if data.is_empty() {
        0.0
    } else {
        data.iter().sum::<f64>() / data.len() as f64
    }
}

/// Project a pressure iterate onto the zero-mean subspace.
pub fn project_pressure(p: &mut [f64]) {
    subtract_mean(p);
}

/// Project a velocity pair (vx, vy) onto the zero-mean-per-component subspace.
///
/// The two components are orthogonal in the periodic null space
/// decomposition, so they project independently.
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
