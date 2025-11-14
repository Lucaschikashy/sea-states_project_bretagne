# Sea States Project - Bretagne

Numerical and statistical study of the sea states that govern an offshore wind development zone off Southern Brittany. The repository combines (1) a data-driven analysis of 27 years of met-ocean hindcasts plus buoy observations, and (2) a physics-based wave-ray tracing solver that explores how incoming Atlantic swell refracts over the regional bathymetry.

## Highlights

- **Mean and Extreme Value Analysis (Part 1):** Jupyter notebook that evaluates trends with OLS and Mann-Kendall tests, classifies wave regimes, and delivers both Block Maxima (BM) and Peak-over-Threshold (POT) return levels using `pyextremes` and `mhkit`.
- **Wave Ray Tracer (Part 2):** Stand-alone Python package (`ocean_wave_tracing.py` with helper modules) that advects hundreds of rays with Runge-Kutta 4, supports time-evolving currents, and exports energy density grids plus xarray datasets.
- **Datasets:** ResourceCode hindcast extracts, Candhis buoy data (`Part_1_Mean_EVA_Waves/data/`), and a cropped GMRT bathymetry grid for ray tracing.
- **Figures:** Generated plots are stored in `Part_1_Mean_EVA_Waves/fig`, `graphs`, and `Part_2_Ray_Tracing/figures`.

## Repository layout

```
.
|- Part_1_Mean_EVA_Waves/
|  |- index_project_EVA_1.ipynb      # Notebook with analysis narrative + code
|  |- data/                          # Hindcast, buoy, and EVA CSV exports
|  |- fig/, graphs/                  # Saved figures, QC plots, tables
|  \- LYue_Sea_States_Intermediate_Report.pdf # Report analysing the results of Part 1
|- Part_2_Ray_Tracing/
|  |- ocean_wave_tracing.py          # Wave_tracing class + RK4 solver
|  |- util_methods.py, util_solvers.py
|  |- demo_CA.py                     # Example configuration & plotting script
|  |- GMRTv4_4_0_20251107topo.asc    # Bathymetry grid (ASCII)
|  \- figures/                       # Output PNGs
|- data/, fig/, graphs/              # Empty hooks for project-level outputs
|- requirements.txt
\- .venv/ (optional local virtual environment)
```

## Getting started

1. **Create an environment (recommended)**
   ```powershell
   python -m venv .venv
   .\.venv\Scripts\Activate
   ```
2. **Install dependencies** (Python 3.10+)
   ```powershell
   pip install --upgrade pip
   pip install -r requirements.txt
   ```
3. **Launch Jupyter Lab or Notebook** for Part 1, or run scripts for Part 2 with `python`.

Tip: the repo already contains large data artifacts (CSV, ASC, PDF). Keep them in place so that relative paths inside the notebook and scripts resolve without extra configuration.

## Part 1 - Mean sea states and EVA workflow

The notebook at `Part_1_Mean_EVA_Waves/index_project_EVA_1.ipynb` follows the ResourceCode wind-turbine coursework brief:

- **Site definition:** Bretagne Sud 1 (47.3236 deg N, -3.5522 deg E), depth derived from GMRT and Candhis buoy metadata.
- **Time series analysis:** Uses `pandas`, `numpy`, `matplotlib`, `resourcecode`, and `pymannkendall` to calculate monthly statistics and check long-term trends in Hs, Tm02, and directional stability.
- **Wave regime classification:** Automatically tags deep, transitional, and shallow conditions from Hs, Tm02, and bathymetry.
- **Seasonal roses:** Produces wind and current roses to describe mean vectors and variability.
- **Buoy vs. hindcast validation:** Root-mean-square differences between Candhis station 05602 (2020) and the hindcast.
- **Extreme value assessment:**
  - BM method on yearly or seasonal blocks.
  - POT method with interactively tuned threshold, bootstrap confidence bounds, and direction checks.
  - Summaries stored in `Part_1_Mean_EVA_Waves/data/*.csv`.
- **Reporting:** Intermediate summary PDF plus ready-made plots in `fig/` and `graphs/`.

### Reproducing the notebook

```powershell
cd Part_1_Mean_EVA_Waves
jupyter lab
```

Run the cells sequentially. Any regenerated outputs will be written back into the `fig/`, `graphs/`, and `data/` subfolders. The notebook expects the CSV files already present; to update them from fresh ResourceCode queries, adapt the corresponding download cells and save to the same filenames.

## Part 2 - Wave ray tracing module

`Part_2_Ray_Tracing` implements a modular ray tracer aimed at surf-beat or WW3-style studies:

- **`Wave_tracing` class:**
  - Handles static or time-varying current fields (`U`, `V`) and bathymetry validation.
  - Provides CFL diagnostics, initial-condition helpers (either along a domain edge or pointwise), and dispersion utilities (intrinsic phase and group velocity).
  - Numerically integrates ray position, wave number, and energy using RK4; stores complete trajectories plus collocated bathymetry and current snapshots.
  - Exports trajectories to `xarray.Dataset` (`to_ds`), georeferences rays (`to_latlon`), and collapses track density into gridded energy (`ray_density`).
- **Utility modules:** `util_methods.py` (xarray helpers, field validation) and `util_solvers.py` (advection and wave-number right-hand sides plus RK4 implementation).
- **`demo_CA.py`:** Example experiment that propagates 300 rays with a 15 s swell over the supplied GMRT grid, writing `wave_rays.png` and `relative_wave_energy_flux_density.png` into `figures/`.

### Running the demo

```powershell
cd Part_2_Ray_Tracing
python demo_CA.py
```

By default the demo assumes quiescent currents and launches rays from the western edge. Edit the script (or build your own) to:

- Swap in current fields (for example an `xarray.DataArray`) and enable `temporal_evolution=True`.
- Change `incoming_wave_side`, `wave_period`, number of rays, or simulation time.
- Provide geographic coordinates via `Wave_tracing.to_latlon()` for map-ready outputs.

`ocean_wave_tracing.log` records solver diagnostics; delete it between runs if you want a clean log.

## Data sources

- **ResourceCode hindcast** - Hourly wave, wind, and current data extracted via the `resourcecode` Python client; see notebook cells for download snippets.
- **Candhis buoy 05602** - CSV included under `Part_1_Mean_EVA_Waves/data/Candhis_05602_2020_arch.csv`.
- **GMRT v4** - Cutout bathymetry (`GMRTv4_4_0_20251107topo.asc`) bundled for convenience. Replace with your own ASCII raster to study another site (ensure spacing metadata is present in the header).

## Troubleshooting and tips

- If figures fail to render in Jupyter, ensure you are using a kernel backed by the virtual environment that has `matplotlib`, `cmocean`, and `mhkit` installed.
- NetCDF inputs for currents should be converted to `xarray.DataArray` objects on `x`, `y`, and optional `time` dimensions; see `check_velocity_field` in `util_methods.py` for the expected shape.
- Large ray counts times time steps accumulate quickly. Watch RAM usage or lower `nt` and `nb_wave_rays` when experimenting interactively.
- Regenerate CSV summaries before sharing results to make sure the BM and POT statistics reflect the latest filtering choices.
