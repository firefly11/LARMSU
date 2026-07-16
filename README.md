# LARMSU

**LARMSU** (LAndslide Runout Model based on Slope Units) is a semi-empirical model for simulating landslide runout extent and susceptibility, particularly in terrain with pronounced slope transitions.

The model combines a slope-unit coordinate system, an energy-conservation-based distance controller, and a diffusion-angle-based direction controller with Monte Carlo simulation. It distinguishes downhill and uphill movement when a landslide crosses different slope units.

## Software and source code

This repository provides:

- a user-friendly Windows GUI in `exe/dist/`;
- source code for reference, reproducibility, and further development; and
- `Dinf` and `aspect` propagation-direction methods.

Repository: [https://github.com/firefly11/LARMSU](https://github.com/firefly11/LARMSU)

## Requirements

- Microsoft Windows
- ArcGIS Desktop 10.x with ArcPy
- ArcGIS Spatial Analyst extension

The supplied executable was built for **ArcGIS Desktop 10.2 with 32-bit ArcPy**.

## Quick start

1. Download the complete `exe/dist` directory. The executable cannot run without the adjacent runtime files.
2. Open `exe/dist/site-packages/102_32bit.pth` and update the ArcGIS and ArcPy paths if necessary.
3. Run `exe/dist/larmsu.exe`.
4. Select a source-area shapefile (`.shp`) and a DEM (`.tif`).
5. Enter the model parameters, select an output directory and project name, and click **Execute**.

The input shapefile and DEM should use a compatible projected coordinate reference system.

## Main parameters

- **Travel angle:** controls frictional energy loss; the general recommended range is 5-30 degrees.
- **Diffusion angle:** controls lateral spreading; the general recommended range is 15-30 degrees.
- **Maximum velocity:** use 100 m/s when no reliable constraint is available.
- **MCMC times:** controls the number of Monte Carlo simulations.
- **ToCF:** controls the scale of slope-unit extraction.

Parameter values should be adjusted according to landslide material, mobility, scale, and terrain. See the accompanying paper for detailed sensitivity analysis and material-specific recommendations.

## Main outputs

- `resultShp.shp`: simulated landslide impact-area boundary.
- `resultf.tif`: normalized susceptibility raster with values from 0 to 1.

Intermediate terrain and slope-unit products are also written to the selected project directory.

## Citation

Liu, J., Gao, X., & Wu, Y. (2026). LARMSU: A semi-empirical landslide runout model based on slope units. *Engineering Geology, 372*, 108964. [https://doi.org/10.1016/j.enggeo.2026.108964](https://doi.org/10.1016/j.enggeo.2026.108964)

## License

Copyright (C) 2026 Liu Jia. The original LARMSU source code is licensed under the GNU General Public License v3.0 only (`GPL-3.0-only`). Third-party components remain subject to their respective licenses.
