# Workflow map

This document maps the manuscript Methods to the code and data files that need to be sorted into the GitHub repository and Zenodo archive.

## 00. Repository rule

GitHub stores:
- reusable Python functions
- workflow scripts
- lightweight notebooks
- documentation
- config templates

Zenodo stores:
- processed data
- figure data
- regional masks
- large NetCDF/Zarr/pickle files
- final archived code snapshot

Raw third-party datasets are cited, not re-uploaded.

# Workflow 01: MODIS CHL field, LOESS background, and eddy collocation

## Step 1: Download/process MODIS CHL

Script:

`scripts/matlab/01_download_modis_chl_and_calculate_loess.m`

Purpose:

- Download/load MODIS-Aqua 8-day chlorophyll-a fields.
- Save the full gridded time series.
- Calculate 6° × 6° LOESS background fields.
- Save `Chl_fields.mat`.

Main outputs:

- `Full_Time_Series.mat`
- `Chl_fields.mat`

These outputs are archived on Zenodo, not stored in GitHub.

## Step 2: Eddy collocation and eddy-induced CHL-a anomaly

Script:

`scripts/matlab/02_collocate_eddy_chl_donut.m`

Purpose:

- Load AE and CE AVISO META eddy trajectories.
- Remove virtual eddies.
- Load raw CHL and LOESS background fields.
- Calculate CHL anomaly as raw CHL minus LOESS background.
- Average raw CHL and CHL anomaly inside normalized eddy-contour rings.

Main output:

- `donut_struct`, saved as a `.mat` file.

This output is archived on Zenodo, not stored in GitHub.

## Required MATLAB helper functions

- `Swath_OSU.m`
- `smooth2d_loess.m`
- `omitvirtualeddy4mMETA.m`
- `nc2mat.m`
- `mat2col.m`
