# References — the relief method (F158 onward)

The PDFs marked **committed** carry an open licence (CC BY 4.0) and are in the repository.
- The others are cited by their reference only.
  - Some are under the publisher's copyright.
  - Others are author versions (HAL, an author's site) with no stated licence.
- A local copy may sit in `docs/refs` (git-ignored by not being added).

| file | reference | in the repo |
|---|---|---|
| `2024-MultiScaleHydro-Author.pdf` | Schott, H., Galin, E., Guérin, E., Peytavie, A., Paris, A. (2024). Terrain Amplification using Multi-scale Erosion. ACM Transactions on Graphics (SIGGRAPH 2024). Author version (the copy the author provided, 2026-10-09; read for F160). | no (author version, no licence stated) |
| `schott2023.pdf` | Schott, H., Paris, A., Fournier, L., Guérin, E., Galin, E. (2023). Large-scale terrain authoring through interactive erosion simulation. ACM Transactions on Graphics. HAL hal-04049125. Code (MIT): github.com/H-Schott/StreamPowerErosion | no (HAL, no licence stated) |
| `cordonnier2016.pdf` | Cordonnier, G., Braun, J., Cani, M.-P., Benes, B., Galin, É., Peytavie, A., Guérin, É. (2016). Large Scale Terrain Generation from Tectonic Uplift and Fluvial Erosion. Computer Graphics Forum 35(2) (Eurographics 2016). | no (publisher's copyright) |
| — | Braun, J., Willett, S. D. (2013). A very efficient O(n), implicit and parallel method to solve the stream power equation governing fluvial incision and landscape evolution. Geomorphology 180–181, 170–179. doi:10.1016/j.geomorph.2012.10.008 | no (no open copy found) |
| `cordonnier2019.pdf` | Cordonnier, G., Bovy, B., Braun, J. (2019). A versatile, linear complexity algorithm for flow routing in topographies with depressions. Earth Surface Dynamics 7, 549–562. doi:10.5194/esurf-7-549-2019 | **committed** (CC BY 4.0) |
| `yuan2019.pdf` | Yuan, X. P., Braun, J., Guerit, L., Rouby, D., Cordonnier, G. (2019). A new efficient method to solve the stream power law model taking into account sediment deposition. JGR Earth Surface 124(6), 1346–1365. doi:10.1029/2018JF004867. HAL hal-02136641 (author version). | no (HAL, no licence stated) |
| `perron2009.pdf` | Perron, J. T., Kirchner, J. W., Dietrich, W. E. (2009). Formation of evenly spaced ridges and valleys. Nature 460, 502–505. doi:10.1038/nature08174 | no (all rights reserved) |
| `salles2020.pdf` | Salles, T., Mallard, C., Zahirovic, S. (2020). gospl: Global Scalable Paleo Landscape Evolution. Journal of Open Source Software 5(56), 2804. doi:10.21105/joss.02804 | **committed** (CC BY 4.0) |
| `armitage2019.pdf` | Armitage, J. J. (2019). Short communication: flow as distributed lines within the landscape. Earth Surface Dynamics 7, 67–75. doi:10.5194/esurf-7-67-2019 | **committed** (CC BY 4.0) |
| `braun2005.pdf` | Braun, J., Robert, X. (2005). Constraints on the rate of post-orogenic erosional decay from low-temperature thermochronological data. Earth Surface Processes and Landforms 30, 1203–1225. doi:10.1002/esp.1271 | no (publisher's copyright) |
| — | Baldwin, J. A., Whipple, K. X., Tucker, G. E. (2003). Implications of the shear stress river incision model for the timescale of postorogenic decay of topography. JGR 108(B3), 2158. doi:10.1029/2001JB000550 | no (no open copy found; abstract only) |
| `kwang2017.pdf` | Kwang, J. S., Parker, G. (2017). Landscape evolution models using the stream power incision model show unrealistic behavior when m/n equals 0.5. Earth Surface Dynamics 5, 807–820. | committed earlier |
| `pelletier2010.pdf` | Pelletier, J. D. (2010). Minimizing the grid-resolution dependence of flow-routing algorithms for geomorphic applications. Geomorphology 122, 91–98. | committed earlier |

## Data (F161 onward)

| data | source | in the repo |
|---|---|---|
| Copernicus DEM GLO-30, five 1-arcsecond tiles over Corsica (N41–N43, E008–E009), the texture reference of F161 | Copernicus Digital Elevation Model (DEM) was accessed on 2026-10-09 from https://registry.opendata.aws/copernicus-dem (`copernicus-dem-30m`, COG tiles). Licence: the Copernicus DEM GLO-30 public licence (free of charge; attribution and a no-liability sentence required) | **no**: the tiles stay outside the repository; the derived grids live in `data/corsica/` (git-ignored); the preparation script is `docs/reports/relief_method/f161_corsica/prep_corse.py` |

**The notice every derived figure carries** (the licence's « modified data » form): « produced using Copernicus
WorldDEM-30 © DLR e.V. 2010-2014 and © Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS by the
European Union and ESA; all rights reserved ».

**The licence's no-liability sentence**: « The organisations in charge of the Copernicus programme by law or by
delegation do not incur any liability for any use of the Copernicus WorldDEM-30. »

**What the data is**: a surface model (DSM) that includes the vegetation and buildings. At ≥ 98 m cells, the canopy
(10–20 m) is small against the relief. Declared, not corrected.
