//! Step 8.6 — translate a flat `V2RunSpec` into a fully-populated
//! `BaselineConfig` ready for `run_baseline`.

use ymir_core::tectonics_v2::age_field::{AgeFieldConfig, AgeFieldConfigEnabled};
use ymir_core::tectonics_v2::basal_drag::{BasalDragConfig, BasalDragLaw};
use ymir_core::tectonics_v2::boundaries::{BoundaryConfig, BoundaryRates};
use ymir_core::tectonics_v2::cratonic::CratonicConfig;
use ymir_core::tectonics_v2::cratonic::CratonicConfigEnabled;
use ymir_core::tectonics_v2::diagnostics::harness::{
    BaselineConfig, ForceKind, NonlinearChoice, build_force,
};
use ymir_core::tectonics_v2::mantle::MantleConfig;
use ymir_core::tectonics_v2::presets::{Preset, YieldingConfig};
use ymir_core::tectonics_v2::recycling::RecyclingConfig;
use ymir_core::tectonics_v2::rheology::YieldingLaw;
use ymir_core::tectonics_v2::scales::Scales;
use ymir_core::tectonics_v2::slab::SlabPullConfig;
use ymir_core::tectonics_v2::stokes::nonlinear_solver::NewtonConfig;
use ymir_core::tectonics_v2::stokes::picard::PicardConfig;
use ymir_core::tectonics_v2::stokes::solver::LinearSolverConfig;
use ymir_core::tectonics_v2::voronoi::VoronoiConfig;

use super::spec::{
    V2AgeFieldSpec, V2CratonicSpec, V2ForceKind, V2LinearSolverSpec, V2MantleSpec, V2RunSpec,
};

/// Build a `BaselineConfig` from the supplied `V2RunSpec`. Mirrors the
/// `build_step10_config` helper from the v2 integration tests, but
/// driven by a serialisable spec rather than hard-coded constants so
/// presets / UI inputs can flow through unchanged.
///
/// Boundary defaults: `k_sub = 0.5`, all other rates `0.0` (matches
/// the Step 8/9/10 baselines). The Voronoï layout is closed-mode (the
/// only mode validated through Steps 6-10).
pub fn build(spec: &V2RunSpec) -> BaselineConfig {
    let scales = Scales::default();
    let preset = Preset::by_name("dynamic-accidented").expect("dynamic-accidented preset exists");

    let vcfg = VoronoiConfig {
        num_plates: spec.num_plates,
        continental_ratio: spec.continental_ratio,
    };
    let rates = BoundaryRates {
        k_sub: 0.5,
        k_arc: 0.0,
        k_spread: 0.0,
        k_coll_v: 0.0,
        k_rift_v: 0.0,
    };
    let boundary = BoundaryConfig::enabled_voronoi_closed(
        spec.grid_nx,
        spec.grid_ny,
        &vcfg,
        spec.seed,
        rates,
        RecyclingConfig::default(),
    )
    .expect("recycling config valid");

    let (force_kind, sin_amplitude) = match spec.force {
        V2ForceKind::Gpe => (ForceKind::Gpe, 0.0),
        V2ForceKind::Sinusoidal { amplitude } => (ForceKind::Sinusoidal, amplitude),
    };
    let force = build_force(force_kind, &scales, sin_amplitude.max(1e-12), 1.0);

    let mantle = match spec.mantle {
        V2MantleSpec::Off => MantleConfig::Disabled,
        V2MantleSpec::On { mf, coupling, num_modes, seed, evolution_rate } => {
            MantleConfig::Enabled { mf, coupling, num_modes, seed, evolution_rate }
        }
    };

    let cratonic = match spec.cratonic {
        V2CratonicSpec::Off => CratonicConfig::Disabled,
        V2CratonicSpec::On {
            cr,
            k_viscous,
            b_factor,
            smoothing_width,
            plate_area_min,
        } => CratonicConfig::Enabled(CratonicConfigEnabled {
            cr,
            k_viscous,
            b_factor,
            smoothing_width,
            plate_area_min,
            ..CratonicConfigEnabled::default()
        }),
    };

    let age_field = match spec.age_field {
        V2AgeFieldSpec::Off => AgeFieldConfig::Disabled,
        V2AgeFieldSpec::On { continental_age_init, oceanic_age_init } => {
            AgeFieldConfig::Enabled(AgeFieldConfigEnabled {
                continental_age_init,
                oceanic_age_init,
            })
        }
    };

    let linear_solver = match spec.linear_solver {
        V2LinearSolverSpec::Jacobi => LinearSolverConfig::default(),
        V2LinearSolverSpec::Amg => LinearSolverConfig::default(), // Phase 5 — wire AmgCG; default until then
    };

    let heightmap_fractions = if spec.capture_endpoints { vec![0.0, 1.0] } else { vec![] };

    BaselineConfig {
        seed: spec.seed,
        grid_nx: spec.grid_nx,
        grid_ny: spec.grid_ny,
        domain_lx: 1.0,
        domain_ly: 1.0,
        steps: spec.steps,
        cfl_factor: spec.cfl_factor,
        total_time_nondim: spec.total_time_nondim,
        preset,
        nonlinear: NonlinearChoice::Newton,
        newton_cfg: NewtonConfig::default(),
        picard_cfg: PicardConfig::default(),
        heightmap_fractions,
        output_dir: spec.output_dir.clone(),
        force,
        force_kind,
        sinusoidal_amplitude: sin_amplitude,
        s_perturbation_amplitude: spec.s_perturbation_amplitude,
        yielding: YieldingConfig::Enabled(YieldingLaw {
            bi: spec.bi,
            ..YieldingLaw::default()
        }),
        basal_drag: BasalDragConfig::Enabled(BasalDragLaw {
            br: spec.br,
            ..BasalDragLaw::default()
        }),
        boundary,
        boundary_layout_name: format!("v2bridge_{}_seed{}", spec.preset_label, spec.seed),
        slab_pull: if spec.slab_enabled {
            // Step 8 left slab co-calibration as a follow-up; in v2
            // bridge we keep it Disabled by default until the
            // calibration lands. The toggle stays in the spec for
            // forward-compat — when slab support arrives the only
            // line that needs editing is this match arm.
            SlabPullConfig::Disabled
        } else {
            SlabPullConfig::Disabled
        },
        mantle,
        cratonic,
        age_field,
        capture: None,
        linear_solver,
        // Step 8.6 Phase 8d — `init_mode` flows through `V2RunSpec`,
        // exposed in the UI via the "Initialisation" section of the
        // parameter panel. Existing preset JSON files that predate the
        // field deserialise to `V2InitModeSpec::default()` (= Uniform)
        // via `#[serde(default)]`.
        init_mode: spec.init_mode.into_core(),
        // Step 8.6 follow-up — `None` for fresh runs; the bridge's
        // `ContinueRun` command path overrides this after `build()`
        // returns.
        continuation: None,
    }
}
