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
├── Data_Samples/
│   ├── City_Limit_Light_Outline.geojson
│   ├── Sample_LST.tif
│   ├── Sample_NDVI.tif
│   ├── Villages.tif
│
├── Docs/
│   ├── NDVI_LST_scatter.png        
│   ├── Phoenix_Workflow.svg
│   ├── SUHI_Map.png
│   ├── villages_zonal_stats_py.csv
│
├── Environments                   # Python and R environments
│   ├── R_environment.yml        
│   ├── environment.yml
│
├── KNIME_Workflow/
│   └── Phoenix_SUHI_Workflow.knwf     # Main KNIME workflow file
│
├── Scripts/
│   ├── Clipping&Masking.R            		             # Clipping each to the city's bounding box.
│   ├── Clipping&Masking_diagnosis.R                   # Diagnosing the clipping of the rasters.
│   ├── Global_NDVI-LST_regression_analysis.R          # NDVI-LST regression analysis per urban villages.
│   ├── LST_compositing.R             		             # Median compositing.
│   ├── LST_preprocessing.R           		             # QC filtration and scaling.
│   ├── LST_resampling.R             	 	               # Upscaling to the 250m grid.
│   ├── LST_resampling_diagnosis.R                     # Diagnosing the resampling of LST to the NDVI grid.
│   ├── NDVI_compositing.R            		             # Median compositing.
│   ├── NDVI_preprocessing.R          		             # QC filtration and scaling.
│   ├── SUHI_computation.R            		             # Quantifying surface urban heat island (SUHI) effect.
│   ├── Stratified_NDVI-LST_regression_analysis.R      # NDVI-LST regression analysis for the whole of Phoenix.
│   ├── UrbanVsRural_Mask.R           		             # Contrasting urban temperatures to the rural reference zone.
│   ├── Zonal_statistics.py           		             # Python zonal statistics for villages.
│
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
conda env create -f environment.yml
```

- **R** (Linked via the Conda Environment Propagation node)
```bash
install.packages(c("terra", "raster", "ggplot2", "dplyr"))
```
Or use the provided environment file:
```bash
conda env create -f R_environment.yml
```

This workflow uses two separate Conda environments — one for Python and one for R — to support different scripting nodes in KNIME.
Each environment can be recreated from the provided .yml files.

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
SUHI(x,y) = LST(x,y) - mean(LSTrural)
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

Rop, Alex; Jain, Devika (2025). Designing Scalable Workflows for Spatial Analysis of Earth Observation Data using KNIME: A Case Study of Phoenix (May–August 2024).

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
