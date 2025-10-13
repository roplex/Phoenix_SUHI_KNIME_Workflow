# Phoenix_SUHI_KNIME_Workflow
A reproducible Earth Observation (EO) analysis workflow built in KNIME integrating R and Python for LST, NDVI, and SUHI mapping over Phoenix, AZ (May–Aug 2024).

# KNIME Workflow for Surface Urban Heat Island (SUHI) Analysis — Phoenix Case Study

This repository contains the complete and reproducible **KNIME workflow**, associated **Python and R scripts**, and documentation used in the study:

> **"Modular EO Workflow for Surface Urban Heat Island Analysis Using KNIME: A Case Study of Phoenix (May–Aug 2024)"**

The workflow demonstrates how **KNIME** can serve as a visual, modular, and reproducible platform for **Earth Observation (EO)** data processing — integrating MODIS NDVI and LST datasets, and applying both R and Python analytics for SUHI computation and zonal statistics.

---

## 🌍 Overview

This end-to-end KNIME workflow automates the main steps for Surface Urban Heat Island (SUHI) analysis in Phoenix, Arizona.  
The workflow includes:

- **Importing MODIS datasets:**  
  - LST (MOD11A2, 1 km, 8-day composites)  
  - NDVI (MOD13Q1, 250 m, 16-day composites)

- **Data preprocessing:** QC filtering, masking, reprojection, and resampling  
- **Median compositing:** Generation of seasonal (May–Aug 2024) LST and NDVI composites  
- **SUHI computation:** Urban–rural LST differentials  
- **NDVI–LST regression:** Assessing vegetation–temperature relationships  
- **Zonal statistics:** Aggregating LST, NDVI, and SUHI by Phoenix urban villages (Python-based)  
- **Visualization:** Histograms, scatter plots, and summary tables  

---

## 🧩 Repository Structure

```bash
📦 KNIME_Phoenix_UHI_Workflow
│
├── KNIME_Workflow/
│   └── Phoenix_UHI_Workflow.knwf     # Main KNIME workflow file
│
├── Scripts/
│   ├── lst_composite.R               # Median compositing in R
│   ├── ndvi_lst_regression.R         # NDVI-LST regression analysis
│   ├── zonal_stats.py                # Python zonal statistics for villages
│
├── Data_Samples/
│   ├── Phoenix_AOI.geojson
│   ├── Sample_LST.tif
│   ├── Sample_NDVI.tif
│
├── Docs/
│   ├── workflow_diagram.png          # Overview diagram of KNIME workflow
│   ├── suhi_map_example.png
│   ├── ndvi_lst_scatter.png
│   ├── table_summary_example.png
│
├── environment.yml                   # Conda environment for Python dependencies
├── LICENSE
└── README.md
