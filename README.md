# KNIME Workflow for Surface Urban Heat Island (SUHI) Analysis — Phoenix Case Study

This repository contains the complete and reproducible **KNIME workflow**, associated **Python and R scripts**, and documentation used in the study:

> **"Harmonizing Multi-Resolution Earth Observation Data through Scalable KNIME Workflows:
A Surface Urban Heat Island Case Study in Phoenix, Arizona"**

The workflow demonstrates how **KNIME** can serve as a visual, modular, and reproducible platform for **Earth Observation (EO)** data processing — integrating MODIS NDVI and LST datasets, and applying both R and Python analytics for SUHI computation and zonal statistics.

⸻

## 🌍 Overview

This end-to-end KNIME workflow automates the main steps for Surface Urban Heat Island (SUHI) analysis in Phoenix, Arizona.  
The workflow includes:

- **Importing MODIS datasets:**  
  - LST (MOD11A2, 8-day, 1km resolution)  
  - NDVI (MOD13Q1, 16-day, 250m resolution)

- **Data preprocessing:** QC filtering, masking, reprojection, and resampling  
- **Median compositing:** Generation of seasonal (May–Aug 2024) LST and NDVI composites  
- **SUHI computation:** Urban–rural LST differentials  
- **NDVI–LST regression (global & stratified) analysis:** Assessing vegetation–temperature relationships  
- **Zonal statistics computation:** Aggregating LST, NDVI, and SUHI by Phoenix urban villages (Python-based)  
- **Visualization:** Histograms, scatter plots, and summary tables  

⸻

## 🧩 Repository Structure

```bash
📦 Phoenix_SUHI_KNIME_Workflow
│
├── Data_Samples/                  # Small sample inputs for trying the nodes
│   ├── City_Limit_Light_Outline.geojson
│   ├── Sample_LST.tif
│   ├── Sample_NDVI.tif
│   └── Villages.geojson
│
├── Docs/                          # Example outputs
│   ├── NDVI_LST_scatter.png
│   ├── Phoenix_Workflow.svg
│   ├── SUHI_map.tif
│   └── villages_zonal_stats_py.csv
│
├── Environments/                  # Python and R Conda environments
│   ├── R_environment.yml
│   └── environment.yml
│
├── KNIME_Workflow/
│   └── Phoenix_SUHI_Workflow.knwf # Main KNIME workflow file
│
├── scripts/                       # Scripts embedded in the KNIME nodes
│   ├── NDVI_preprocessing.R / LST_preprocessing.R      # QC filtration and scaling
│   ├── NDVI_compositing.R / LST_compositing.R          # Seasonal median compositing
│   ├── LST_resampling.R                                # TsHARP sharpening (tsharp_disaggregate) + diagnose_sharpening_quality
│   ├── LST_resampling_diagnosis.R                      # Diagnostics for the LST resampling step
│   ├── Clipping&Masking.R / Clipping&Masking_diagnosis.R
│   ├── UrbanVsRural_Mask.R                             # Elevation-matched rural reference zone and baseline
│   ├── rural_reference_zone_redefinition.R             # Standalone Census Urbanized Area + ESA WorldCover mask
│   ├── rural_reference_zone_sensitivity_analysis.R     # Standalone sensitivity analysis of the rural baseline (run in RStudio)
│   ├── SUHI_computation.R
│   ├── Global_NDVI-LST_regression_analysis.R           # City-wide (pooled) regression
│   ├── Stratified_NDVI-LST_regression_analysis.R       # Per-village regression
│   └── Zonal_statistics.py                             # Python zonal statistics for villages
│
├── .gitignore
├── LICENSE
└── README.md
```


⸻


## ⚙️ Dependencies

**KNIME**
- KNIME Analytics Platform ≥ 5.5
- R Integration Extension
- Python Integration Extension

**Environment Setup**

- **Python** (Set under Preferences → KNIME → Python)

Install with Conda:
```bash
conda install -c conda-forge rasterio geopandas rasterstats pandas numpy
```
Or use the provided environment file:
```bash
conda env create -f Environments/environment.yml
```

- **R** (Linked via the Conda Environment Propagation node)
```bash
install.packages(c("terra", "raster", "ggplot2", "dplyr", "tigris", "elevatr"))
```
Or use the provided environment file:
```bash
conda env create -f Environments/R_environment.yml
```

This workflow uses two separate Conda environments — one for Python and one for R — to support different scripting nodes in KNIME.
Each environment can be recreated from the provided .yml files.

⸻


## 📥 Data Acquisition and Local Paths

Full-extent inputs are not stored in this repository. To reproduce the reported values:

1. **MODIS** – submit an AppEEARS area-sample request (https://appeears.earthdatacloud.nasa.gov/task/area) for the Phoenix extent, 1 May – 31 Aug 2024, for MOD11A2.061 (`LST_Day_1km` with `QC_Day`) and MOD13Q1.061 (`250m_16_days_NDVI` with `250m_16_days_VI_Quality`), GeoTIFF output.
2. **Boundaries** – City of Phoenix city limit and urban villages (https://www.phoenixopendata.com); US Census 2020 Urbanized Areas (retrieved with the `tigris` R package).
3. **Land cover** – ESA WorldCover v200 tile N33W114 (https://esa-worldcover.org).
4. **Elevation** – SRTM-derived DEM retrieved with the `elevatr` R package.

Scripts read a data root from the environment variable `PHX_ROOT` (default: a `PhoenixData2` folder in the working directory) and expect sub-folders such as `Boundary/`, `LST/`, `Preprocessed/`, and `WorldCover/`. Set it before running, e.g. `Sys.setenv(PHX_ROOT = "/path/to/PhoenixData2")` in R, or `export PHX_ROOT=/path/to/PhoenixData2` for the Python node. Values inside KNIME nodes that arrive as flow variables (e.g. `lst_wide_path`) are unchanged.

The sensitivity analysis of the rural reference zone (`scripts/rural_reference_zone_sensitivity_analysis.R`) runs outside KNIME in R/RStudio and writes its outputs to `<PHX_ROOT>/Sensitivity/`.


⸻


## 🧮 Key Equations

(1) LST Median Composite:
```bash
LSTcomposite(x,y) = median(LSTt(x,y)),		∀ t ∈ May - Aug 2024
```

(2) NDVI Median Composite:
```bash
NDVIcomposite(x,y) = median(NDVIt(x,y)),		∀ t ∈ May - Aug 2024
```

(3) Surface Urban Heat Island (SUHI):
```bash
SUHI(x,y) = LST(x,y) - mean(LSTrural)     # rural = non-built-up land outside the Census Urbanized Area, elevation-matched to the urban core
```

(4) NDVI–LST Regression:
```bash
LST = 𝜶 + 𝜷 · NDVI
```


⸻


## 📊 Outputs	
- **Composite maps:** LST and NDVI median composites
- **SUHI GeoTIFF and histogram of intensity distribution**
- **NDVI–LST scatterplot and regression statistics**
- **Village-level zonal statistics (CSV + GeoJSON)**
- **Table View summaries of LST, NDVI, and SUHI metrics**


⸻


## 🔍 Notes on Data Quality
- **Missing zonal statistics in some villages (e.g., Encanto) are linked to missing raster cells caused by:**
  - Masking and QC filtering steps
  - Resampling mismatches between 1 km and 250 m resolutions
  - Edge effects during reprojection
- **These spatial gaps appear as empty regions in the composite or SUHI maps.**

**Documenting and managing such artifacts ensures reproducibility and transparency.**


⸻


## 🌍 Reproducibility and Reuse

**This workflow is designed for:**
  - Teaching and capacity building in EO data analytics
  - Rapid prototyping of urban climate monitoring pipelines
  - Extensibility to other cities or sensors (e.g., Landsat, Sentinel-3)


⸻


## 📦 Citation

If you use this repository, please cite:

Rop, Alex; Jain, Devika (2026). Harmonizing Multi-Resolution Earth Observation Data through Scalable KNIME Workflows: A Surface Urban Heat Island Case Study in Phoenix, Arizona.

[GitHub Repository](https://github.com/roplex/Phoenix_SUHI_KNIME_Workflow)


⸻


## ✉️ Contact

Author: Rop K. Alex | Geospatial Engineer

Email: ropalex44@gmail.com

Affiliation: Spatial Data Lab, Center for Geographic Analysis - Harvard University

GitHub: https://github.com/roplex

Location: Nairobi, Kenya


⸻


## 🧠 License

This project is distributed under the MIT License.

You are free to reuse and adapt with proper attribution.
