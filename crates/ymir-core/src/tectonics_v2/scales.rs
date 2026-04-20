//! Physical scales and dim ↔ nondim conversions.
//!
//! Primary scales (length, thickness, time, density) are chosen; the
//! remaining scales are derived. The reference choice is documented in
//! `docs/solver-scaling.md` §3 and yields a viscosity scale
//! η* = ρ*·g·τ*·S* ≈ 1.07×10²⁴ Pa·s with the default primary values.

use std::f64::consts::PI;

/// Seconds in one million years (Julian year × 1e6).
pub const SECONDS_PER_MYR: f64 = 365.25 * 24.0 * 3600.0 * 1.0e6;

/// Standard gravitational acceleration (m/s²).
pub const GRAVITY: f64 = 9.81;

/// Physical scales used to non-dimensionalize the tectonic solver.
///
/// All fields are in SI. Primary scales are the four the user sets;
/// derived scales are computed from them at construction time.
#[derive(Clone, Copy, Debug)]
pub struct Scales {
    // ---- Primary (configurable) ----
    /// Horizontal length scale L* (m). Default 350 km.
    pub length: f64,
    /// Crustal thickness scale S* (m). Default 35 km.
    pub thickness: f64,
    /// Tectonic time scale τ* (s). Default 30 Myr.
    pub time: f64,
    /// Reference density ρ* (kg/m³). Default mantle 3300.
    pub density: f64,

    // ---- Derived ----
    /// Velocity scale v* = L*/τ* (m/s).
    pub velocity: f64,
    /// Strain-rate scale ε̇* = 1/τ* (1/s).
    pub strain_rate: f64,
    /// Viscosity scale η* = ρ*·g·τ*·S* (Pa·s).
    pub viscosity: f64,
    /// Stress scale σ* = η*·ε̇* (Pa).
    pub stress: f64,
    /// Lithostatic pressure scale p* = ρ*·g·S* (Pa).
    pub pressure: f64,
    /// Body-force scale (per unit volume) f* = σ*/L* (N/m³).
    pub body_force: f64,
    /// Argand number (diagnostic): Ar = ρ*·g·S*²/(η*·v*).
    pub argand: f64,
}

impl Default for Scales {
    fn default() -> Self {
        Self::from_primary(350.0e3, 35.0e3, 30.0 * SECONDS_PER_MYR, 3300.0)
    }
}

impl Scales {
    /// Build a scale set from the four primary scales.
    pub fn from_primary(length: f64, thickness: f64, time: f64, density: f64) -> Self {
        let velocity = length / time;
        let strain_rate = 1.0 / time;
        let viscosity = density * GRAVITY * time * thickness;
        let stress = viscosity * strain_rate;
        let pressure = density * GRAVITY * thickness;
        let body_force = stress / length;
        let argand = density * GRAVITY * thickness * thickness / (viscosity * velocity);
        Self {
            length,
            thickness,
            time,
            density,
            velocity,
            strain_rate,
            viscosity,
            stress,
            pressure,
            body_force,
            argand,
        }
    }

    /// Nondim length x/L*.
    pub fn to_nondim_length(&self, x: f64) -> f64 { x / self.length }
    pub fn to_dim_length(&self, xt: f64) -> f64 { xt * self.length }

    pub fn to_nondim_thickness(&self, s: f64) -> f64 { s / self.thickness }
    pub fn to_dim_thickness(&self, st: f64) -> f64 { st * self.thickness }

    pub fn to_nondim_time(&self, t: f64) -> f64 { t / self.time }
    pub fn to_dim_time(&self, tt: f64) -> f64 { tt * self.time }

    pub fn to_nondim_velocity(&self, v: f64) -> f64 { v / self.velocity }
    pub fn to_dim_velocity(&self, vt: f64) -> f64 { vt * self.velocity }

    pub fn to_nondim_viscosity(&self, eta: f64) -> f64 { eta / self.viscosity }
    pub fn to_dim_viscosity(&self, etat: f64) -> f64 { etat * self.viscosity }

    pub fn to_nondim_stress(&self, sigma: f64) -> f64 { sigma / self.stress }
    pub fn to_dim_stress(&self, st: f64) -> f64 { st * self.stress }

    pub fn to_nondim_pressure(&self, p: f64) -> f64 { p / self.pressure }
    pub fn to_dim_pressure(&self, pt: f64) -> f64 { pt * self.pressure }

    pub fn to_nondim_body_force(&self, f: f64) -> f64 { f / self.body_force }
    pub fn to_dim_body_force(&self, ft: f64) -> f64 { ft * self.body_force }

    /// Emit a human-readable summary used as a solver-startup fingerprint
    /// and in the Step 0 diagnostics report.
    pub fn report(&self) -> String {
        let myr = self.time / SECONDS_PER_MYR;
        let cm_per_yr = self.velocity * 100.0 * 365.25 * 24.0 * 3600.0;
        format!(
            "Scales:\n\
             - L* = {:.3e} m ({:.0} km)\n\
             - S* = {:.3e} m ({:.0} km)\n\
             - τ* = {:.3e} s ({:.2} Myr)\n\
             - ρ* = {:.1} kg/m³\n\
             - v* = {:.3e} m/s ({:.3} cm/yr)\n\
             - ε̇* = {:.3e} 1/s\n\
             - η* = {:.3e} Pa·s\n\
             - σ* = {:.3e} Pa\n\
             - p* = {:.3e} Pa\n\
             - f* = {:.3e} N/m³\n\
             - Ar = {:.3}\n\
             - 2π check (not used, informational): {:.6}",
            self.length,
            self.length / 1e3,
            self.thickness,
            self.thickness / 1e3,
            self.time,
            myr,
            self.density,
            self.velocity,
            cm_per_yr,
            self.strain_rate,
            self.viscosity,
            self.stress,
            self.pressure,
            self.body_force,
            self.argand,
            2.0 * PI,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, rel: f64) -> bool {
        (a - b).abs() <= rel * a.abs().max(b.abs()).max(1.0)
    }

    #[test]
    fn default_matches_design_note() {
        let s = Scales::default();
        // η* = ρ*·g·τ*·S*
        //    = 3300 × 9.81 × (30·3.156e13) × 3.5e4
        //    ≈ 1.07e24 Pa·s
        assert!(approx_eq(s.viscosity, 1.07e24, 5e-2), "η* = {:.3e}", s.viscosity);
        // v* = L*/τ* ≈ 3.7e-10 m/s ≈ 1.17 cm/yr
        assert!(approx_eq(s.velocity, 3.7e-10, 5e-2), "v* = {:.3e}", s.velocity);
        // Ar should be O(1) by construction
        assert!(s.argand > 0.05 && s.argand < 20.0, "Ar = {}", s.argand);
    }

    #[test]
    fn length_roundtrip() {
        let s = Scales::default();
        for x in [1.0, 1.23e4, 9.87e7] {
            assert!(approx_eq(s.to_dim_length(s.to_nondim_length(x)), x, 1e-12));
        }
    }

    #[test]
    fn all_roundtrips_are_exact() {
        let s = Scales::default();
        let x = 1.234567e8;
        assert!(approx_eq(s.to_dim_length(s.to_nondim_length(x)), x, 1e-12));
        assert!(approx_eq(s.to_dim_thickness(s.to_nondim_thickness(x)), x, 1e-12));
        assert!(approx_eq(s.to_dim_time(s.to_nondim_time(x)), x, 1e-12));
        assert!(approx_eq(s.to_dim_velocity(s.to_nondim_velocity(x)), x, 1e-12));
        assert!(approx_eq(s.to_dim_viscosity(s.to_nondim_viscosity(x)), x, 1e-12));
        assert!(approx_eq(s.to_dim_stress(s.to_nondim_stress(x)), x, 1e-12));
        assert!(approx_eq(s.to_dim_pressure(s.to_nondim_pressure(x)), x, 1e-12));
        assert!(approx_eq(s.to_dim_body_force(s.to_nondim_body_force(x)), x, 1e-12));
    }

    #[test]
    fn derived_scales_consistent() {
        let s = Scales::default();
        assert!(approx_eq(s.velocity, s.length / s.time, 1e-14));
        assert!(approx_eq(s.strain_rate, 1.0 / s.time, 1e-14));
        assert!(approx_eq(s.viscosity, s.density * GRAVITY * s.time * s.thickness, 1e-14));
        assert!(approx_eq(s.stress, s.viscosity * s.strain_rate, 1e-14));
        assert!(approx_eq(s.pressure, s.density * GRAVITY * s.thickness, 1e-14));
        assert!(approx_eq(s.body_force, s.stress / s.length, 1e-14));
    }

    #[test]
    fn report_contains_key_values() {
        let r = Scales::default().report();
        assert!(r.contains("η*"));
        assert!(r.contains("Ar"));
    }
}
