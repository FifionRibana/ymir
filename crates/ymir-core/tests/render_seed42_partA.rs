//! PART A verification — render seed 42 at the RESTORED defaults (tlf None,
//! 8 plates, window 328 km) and report HD emerged fraction + a coastline speckle
//! metric, to confirm the island (and crisp coasts) are back.
//! Run: cargo test -p ymir-core --test render_seed42_partA --release -- --ignored --nocapture

use std::time::Instant;

use ymir_core::tectonics::isostasy::IsostasyConfig;
use ymir_core::tectonics_c1::cached_product::{c1_land_centroid, cached_c1_eroded};
use ymir_core::tectonics_c1::closures::oceanic_bathymetry::params::SteinSteinParams;
use ymir_core::tectonics_c1::init_r7::Phase2InitParams;
use ymir_core::tectonics_c1::production_upscale::C1_DOMAIN_KM;
use ymir_core::tectonics_c1::time_loop::{C1Closures, C1TimeLoopConfig};
use ymir_core::terrain::upscale::FbmUpscaleConfig;

#[test]
#[ignore]
fn render_seed42_partA() {
    let (seed, grid, target) = (42u64, 64usize, 2048usize);
    let window_km = 328.0f32;
    let init = Phase2InitParams::default(); // 8 plates / cc 1 (restored default)
    let run = C1TimeLoopConfig {
        rigid_continental_crust: true,
        n_steps: 300,
        dx: 1.0 / grid as f64,
        dy: 1.0 / grid as f64,
        iso_config: IsostasyConfig::c1_default(),
        drainage_max_distance: 30,
    };
    let clo = C1Closures::default();
    let ss = SteinSteinParams::default();

    let wf = (window_km / C1_DOMAIN_KM) as f64;
    let cache = std::env::temp_dir().join("ymir_seed42_partA_cache");
    let _ = std::fs::create_dir_all(&cache);
    let sea = 0.5f32;

    // Compare the RESTORED default (tlf None) against the regressed tlf 0.08, same
    // seed/window, so speckle is judged relative — not against an absolute guess.
    for tlf in [None, Some(0.08f32)] {
        let mut upscale = FbmUpscaleConfig::c1_hd_production(target);
        upscale.target_land_fraction = tlf;
        let c = c1_land_centroid(seed, grid, &init, &run, &clo, &ss, tlf);
        let origin =
            [(c[0] - wf * 0.5).clamp(0.0, 1.0 - wf), (c[1] - wf * 0.5).clamp(0.0, 1.0 - wf)];
        upscale.sample_origin = origin;
        upscale.sample_size = wf;

        let t = Instant::now();
        let eroded = cached_c1_eroded(
            &cache,
            seed,
            grid,
            &init,
            &run,
            &clo,
            &ss,
            &upscale,
            &ymir_core::tectonics_c1::closures::volcanism::VolcanismConfig::default(),
        )
        .unwrap()
        .heightmap;
        let (w, h) = (eroded.width, eroded.height);
        let n = w * h;
        let d = &eroded.data;
        let is_land = |k: usize| d[k] > sea;
        let emerged = (0..n).filter(|&k| is_land(k)).count() as f32 / n as f32;

        // Speckle: interior cells whose 4-neighbourhood is ≥3 of the OPPOSITE type
        // (isolated spits / pinholes), per 1000 coastline cells (comparable across
        // different coastline lengths).
        let (mut speckle, mut coast) = (0usize, 0usize);
        for y in 1..h - 1 {
            for x in 1..w - 1 {
                let k = y * w + x;
                let land = is_land(k);
                let opp =
                    [k - 1, k + 1, k - w, k + w].iter().filter(|&&nk| is_land(nk) != land).count();
                if opp > 0 {
                    coast += 1;
                }
                if opp >= 3 {
                    speckle += 1;
                }
            }
        }
        let per_1000 = if coast > 0 { speckle as f32 / coast as f32 * 1000.0 } else { 0.0 };
        let tlf_s = tlf.map(|f| format!("{f:.2}")).unwrap_or_else(|| "None".into());
        eprintln!(
            "seed 42 @ {w}² tlf={tlf_s:>4} in {:.0} s: emerged {:>5.1}%  coast {coast:>6}  \
             speckle {speckle:>5} ({per_1000:>5.1} / 1000 coast)",
            t.elapsed().as_secs_f32(),
            emerged * 100.0,
        );
    }
}
