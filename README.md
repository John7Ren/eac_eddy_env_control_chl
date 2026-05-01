# Environmental controls on eddy-induced chlorophyll-a anomalies in the East Australian Current System

This repository contains code for the manuscript:

**Environmental controls on eddy-induced chlorophyll-a anomalies in the East Australian Current System**

The workflow processes satellite chlorophyll-a fields, eddy trajectory data, and background ocean climatologies to diagnose how eddies influence chlorophyll-a anomalies across the East Australian Current (EAC) system.

## Repository status

This repository is being prepared as the reproducible code archive for Paper 1.

The code currently includes:

- MATLAB scripts for MODIS chlorophyll-a processing, LOESS background calculation, and eddy collocation.
- Python modules for loading data, preparing background fields, building regional climatologies, and plotting manuscript figures.
- Figure scripts for Figures 5–8.

Large data files and generated outputs are not stored in this GitHub repository. They are archived separately on Zenodo.

## Repository structure

```text
eac_eddy_env_control_chl/
├── README.md
├── environment.yml
├── pyproject.toml
├── config/
│   ├── paths_template.yml
│   └── paths_local.yml              # local only; not committed
├── src/
│   ├── matlab/
│   │   └── utils/
│   └── python/
│       └── eac_eddy_chl/
│           ├── io.py
│           ├── time_utils.py
│           ├── grid_utils.py
│           ├── regions.py
│           ├── eddy_utils.py
│           ├── background_fields.py
│           ├── climatology.py
│           └── stats_utils.py
├── scripts/
│   ├── matlab/
│   │   ├── step01_download_modis_chl_and_calculate_loess.m
│   │   └── step02_collocate_eddy_chl_donut.m
│   └── python/
│       ├── 01_load_inputs.py
│       ├── 02_prepare_background_fields.py
│       ├── 04_build_climatologies.py
│       └── figures/
│           ├── make_fig05_07_depth_climatology.py
│           └── make_fig08_lifecycle.py
└── docs/
