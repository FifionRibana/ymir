//! Boundary-flag enum + static boundary-condition config for Step 5.
//!
//! [`BoundaryFlag`] tags each cell with the boundary mechanism it is
//! subject to (if any). Flags are **prescribed statically** at
//! simulation start — dynamic detection from the velocity field is
//! Step 6.
//!
//! [`BoundaryConfig`] is the enable/disable switch wired through
//! [`crate::tectonics_v2::diagnostics::harness::BaselineConfig`]. When
//! `Disabled`, the source/sink pipeline is **structurally bypassed**
//! (no Q evaluation, no clamp, no tracking): the `S̃` evolution
//! reduces to plain advection, same as Step 0–4. This is what the
//! Step 5 regression's zero-cost invariant requires.

use super::super::field::Field2D;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BoundaryFlag {
    /// Cell is not on an active boundary. `Q_sub, Q_arc, Q_coll-v,
    /// Q_rift-v` all evaluate to zero on it.
    None,
    /// Generic subduction flag. Oceanic subducting under an adjacent
    /// plate; consumes mass via `Q_sub`.
    Subduction,
    /// Mid-ocean ridge / continental rift. Produces mass via
    /// `Q_spread` (oceanic) or `Q_rift-v` (continental).
    Rift,
    /// Continental collision zone (active orogen). Produces mass via
    /// `Q_coll-v`.
    ContinentalCollision,
    /// Explicit oceanic-subduction variant. Functionally equivalent
    /// to `Subduction` for `Q_sub`, kept as a distinct tag so layout
    /// generators can advertise the physical intent (e.g., an
    /// oceanic plate subducting under a continent vs. a generic
    /// trench).
    OceanicSubduction,
}

impl BoundaryFlag {
    /// Returns `true` for flags that mark a subducting cell
    /// (consumes mass via `Q_sub`).
    pub fn is_subduction(&self) -> bool {
        matches!(self, BoundaryFlag::Subduction | BoundaryFlag::OceanicSubduction)
    }
}

/// Cell-centred boundary-flag field, shape `nx × ny`.
#[derive(Clone, Debug)]
pub struct BoundaryFlagField {
    nx: usize,
    ny: usize,
    data: Vec<BoundaryFlag>,
}

impl BoundaryFlagField {
    pub fn filled(nx: usize, ny: usize, f: BoundaryFlag) -> Self {
        Self {
            nx,
            ny,
            data: vec![f; nx * ny],
        }
    }

    pub fn nx(&self) -> usize { self.nx }
    pub fn ny(&self) -> usize { self.ny }

    #[inline]
    pub fn get(&self, i: usize, j: usize) -> BoundaryFlag {
        self.data[j * self.nx + i]
    }

    #[inline]
    pub fn set(&mut self, i: usize, j: usize, f: BoundaryFlag) {
        self.data[j * self.nx + i] = f;
    }

    pub fn data(&self) -> &[BoundaryFlag] { &self.data }

    /// Encoding for the layout PNG: None=0, Rift=1, Subduction=2,
    /// OceanicSubduction=3, ContinentalCollision=4.
    pub fn to_heightmap(&self) -> Field2D {
        let mut out = Field2D::new(self.nx, self.ny);
        for j in 0..self.ny {
            for i in 0..self.nx {
                let v = match self.get(i, j) {
                    BoundaryFlag::None => 0.0,
                    BoundaryFlag::Rift => 1.0,
                    BoundaryFlag::Subduction => 2.0,
                    BoundaryFlag::OceanicSubduction => 3.0,
                    BoundaryFlag::ContinentalCollision => 4.0,
                };
                out.set(i, j, v);
            }
        }
        out
    }
}

/// Rate coefficients for the five source/sink terms. See issue #89
/// D3 and `solver-scaling.md` §5.1 for target ranges.
#[derive(Clone, Copy, Debug)]
pub struct BoundaryRates {
    pub k_sub: f64,
    pub k_arc: f64,
    pub k_spread: f64,
    pub k_coll_v: f64,
    pub k_rift_v: f64,
}

impl BoundaryRates {
    /// Baseline values from issue #89 D3. `k_spread` is a **placeholder**;
    /// the real value must come from
    /// [`super::calibration::calibrate_k_spread`]. Using these defaults
    /// without calibration gives a run that will not hit
    /// `s_oceanic_mean ∈ [0.18, 0.22]`.
    pub fn baseline_uncalibrated() -> Self {
        Self {
            k_sub: 0.5,
            k_arc: 0.15,
            k_spread: 0.5,
            k_coll_v: 0.05,
            k_rift_v: 0.02,
        }
    }

    pub fn with_k_spread(mut self, k: f64) -> Self {
        self.k_spread = k;
        self
    }

    pub fn with_k_sub(mut self, k: f64) -> Self {
        self.k_sub = k;
        self
    }

    /// Physical-ordering sanity check from `solver-scaling.md` §4.7:
    /// `k_sub > k_arc > k_coll_v, k_rift_v`. Violating this ordering
    /// is a sign of a mis-specified preset, not a runtime bug —
    /// callers may treat a `false` return as grounds to refuse to
    /// start a run.
    pub fn ordering_is_physical(&self) -> bool {
        self.k_sub > self.k_arc
            && self.k_arc > self.k_coll_v
            && self.k_arc > self.k_rift_v
    }
}

/// Boundary-mechanism enable/disable flag.
///
/// `Disabled` by-passes the source/sink pipeline structurally — no
/// Q evaluation, no clamp, no tracking. Callers feed
/// [`BoundaryConfig::Disabled`] through [`BaselineConfig`] for the
/// Step 5 regression and for Steps 0–4 compatibility (default). The
/// `Enabled` variant carries the plate-type field, the boundary-flag
/// field, and the rate coefficients.
///
/// Boxed because the two fields together can grow to 128²·(1 + 1) =
/// ~32 KiB of enum data per field; at 512² that's ~500 KiB — still
/// Stack-heavy but tolerable. We keep the enum un-boxed for API
/// simplicity; callers that hold long-lived `BoundaryConfig` values
/// can wrap it themselves.
/// Step 6 — recycling mode for the source/sink pipeline.
///
/// - `Open` (Step 5 behaviour) applies each rate per-cell:
///   `Q_arc = k_arc · Σ_voisins_subducting |Q_sub_voisin|`, etc. No
///   buffer, no cross-step flux accounting beyond the clamp flux.
/// - `Closed` (Step 6) drives all creation from the subduction
///   budget via distributive fractions, routing the spread portion
///   through a delayed ring buffer with rollover. `k_spread` is no
///   longer a rate — it disappears from [`BoundaryRates`] for this
///   mode and is replaced by [`super::super::recycling::RecyclingConfig::spread_fraction`].
#[derive(Clone, Debug)]
pub enum RecyclingModeInit {
    Open,
    Closed(super::super::recycling::RecyclingConfig),
}

impl RecyclingModeInit {
    pub fn label(&self) -> &'static str {
        match self {
            RecyclingModeInit::Open => "open",
            RecyclingModeInit::Closed(_) => "closed",
        }
    }

    pub fn describe(&self) -> String {
        match self {
            RecyclingModeInit::Open => "Open (Step 5 rate-based)".to_string(),
            RecyclingModeInit::Closed(cfg) => format!(
                "Closed (arc={:.3}, coll_v={:.3}, rift_v={:.3}, spread={:.3}, mantle_loss={:.3}, delay={} steps)",
                cfg.arc_fraction, cfg.coll_v_fraction, cfg.rift_v_fraction,
                cfg.spread_fraction, cfg.mantle_loss_fraction, cfg.mantle_delay_steps,
            ),
        }
    }
}

#[derive(Clone, Debug)]
pub enum BoundaryConfig {
    Disabled,
    Enabled {
        /// Initial crust geometry (plate_type, plate_id, initial
        /// boundary_flag, detection_config, geometry_kind). Shared
        /// via `Arc` so sweep runners don't deep-copy the fields.
        geometry: super::crust_geometry::CrustGeometryShared,
        /// Rate coefficients for the Step 5 Open-mode source/sink
        /// pipeline. In Closed mode, only `k_sub` is used; the other
        /// rates are ignored in favour of the recycling fractions.
        rates: BoundaryRates,
        /// Recycling mode — Open (Step 5 per-cell rates) or Closed
        /// (Step 6 budget-distributive fractions).
        recycling_mode: RecyclingModeInit,
    },
}

/// Shared (`Arc`) wrapper around [`super::plate_type::PlateTypeField`].
///
/// The harness clones `BaselineConfig` into per-scenario runs and
/// passes `&BaselineConfig` into inner helpers; `Arc` avoids deep
/// copies of the plate-type field (~nx·ny enum bytes) on each clone.
pub type PlateTypeFieldShared = std::sync::Arc<super::plate_type::PlateTypeField>;

/// Shared (`Arc`) wrapper around [`BoundaryFlagField`].
pub type BoundaryFlagFieldShared = std::sync::Arc<BoundaryFlagField>;

impl BoundaryConfig {
    /// Short stable label for reports and logs.
    pub fn label(&self) -> &'static str {
        match self {
            BoundaryConfig::Disabled => "disabled",
            BoundaryConfig::Enabled { .. } => "enabled",
        }
    }

    /// Step 5 back-compat builder: accept a plate_type field, a
    /// boundary_flag field, rates, and a layout name; wrap into a
    /// static `CrustGeometry` and build an Enabled config with
    /// `RecyclingModeInit::Open`.
    pub fn enabled_static(
        plate_type: super::plate_type::PlateTypeField,
        boundary_flag: BoundaryFlagField,
        layout_name: impl Into<String>,
        rates: BoundaryRates,
    ) -> Self {
        let geometry = super::crust_geometry::CrustGeometry::from_static(
            plate_type,
            boundary_flag,
            layout_name,
        );
        Self::Enabled {
            geometry: std::sync::Arc::new(geometry),
            rates,
            recycling_mode: RecyclingModeInit::Open,
        }
    }

    /// Step 6 builder: Voronoi tessellation with Closed-mode
    /// recycling. Validates the recycling config at construction —
    /// caller receives an `Err` if the fractions don't sum to 1
    /// or any other invariant is violated.
    pub fn enabled_voronoi_closed(
        nx: usize,
        ny: usize,
        voronoi_config: &super::super::voronoi::VoronoiConfig,
        seed: u64,
        rates: BoundaryRates,
        recycling_config: super::super::recycling::RecyclingConfig,
    ) -> Result<Self, super::super::recycling::RecyclingConfigError> {
        recycling_config.validate()?;
        let geometry =
            super::crust_geometry::CrustGeometry::from_voronoi(nx, ny, voronoi_config, seed);
        Ok(Self::Enabled {
            geometry: std::sync::Arc::new(geometry),
            rates,
            recycling_mode: RecyclingModeInit::Closed(recycling_config),
        })
    }

    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "disabled" | "off" => Ok(BoundaryConfig::Disabled),
            "enabled" | "on" => Err(
                "--boundary-config=enabled requires a --layout argument; \
                 the CLI builds the enabled variant from a named layout"
                    .to_string(),
            ),
            other => Err(format!(
                "unknown --boundary-config value '{}'; expected disabled|enabled",
                other,
            )),
        }
    }

    /// Human-readable descriptor threaded into `SolverConfigDump`.
    pub fn describe(&self) -> String {
        match self {
            BoundaryConfig::Disabled => "Disabled".to_string(),
            BoundaryConfig::Enabled {
                rates,
                recycling_mode,
                geometry,
            } => format!(
                "Enabled [{}] (layout='{}', k_sub={:.3}, k_arc={:.3}, k_spread={:.3}, k_coll-v={:.3}, k_rift-v={:.3})",
                recycling_mode.describe(),
                geometry.layout_name,
                rates.k_sub, rates.k_arc, rates.k_spread, rates.k_coll_v, rates.k_rift_v,
            ),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn subduction_flags_detect_is_subduction() {
        assert!(BoundaryFlag::Subduction.is_subduction());
        assert!(BoundaryFlag::OceanicSubduction.is_subduction());
        assert!(!BoundaryFlag::None.is_subduction());
        assert!(!BoundaryFlag::Rift.is_subduction());
        assert!(!BoundaryFlag::ContinentalCollision.is_subduction());
    }

    #[test]
    fn rates_baseline_has_physical_ordering() {
        let r = BoundaryRates::baseline_uncalibrated();
        assert!(r.ordering_is_physical());
    }

    #[test]
    fn describe_shows_all_rates_when_enabled() {
        let plate_type = super::super::plate_type::PlateTypeField::filled(
            2, 2, super::super::plate_type::PlateType::Continental,
        );
        let flags = BoundaryFlagField::filled(2, 2, BoundaryFlag::None);
        let cfg = BoundaryConfig::enabled_static(
            plate_type,
            flags,
            "test",
            BoundaryRates::baseline_uncalibrated(),
        );
        let desc = cfg.describe();
        assert!(desc.contains("k_sub=0.500"));
        assert!(desc.contains("k_arc=0.150"));
        assert!(desc.contains("k_spread=0.500"));
    }

    #[test]
    fn label_matches_variant() {
        assert_eq!(BoundaryConfig::Disabled.label(), "disabled");
    }
}
