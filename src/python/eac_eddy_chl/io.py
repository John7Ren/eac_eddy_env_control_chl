"""Input/output utilities for the EAC eddy chlorophyll-a workflow."""

from __future__ import annotations

import os
import pickle
from glob import glob
from pathlib import Path
from typing import Any

import numpy as np
import xarray as xr

try:
    import mat73
except ImportError:  # pragma: no cover
    mat73 = None

try:
    import scipy.io
except ImportError:  # pragma: no cover
    scipy = None

from .background_fields import calc_woa_density_n2, k490_to_zeu
from .eddy_utils import add_day_of_life, donut_velsshr_chl, subset_dict_by_time
from .time_utils import argo_juld_to_month, argo_juld_to_timestamp, argo_juld_to_tstep8

FilePath = str | os.PathLike[str]


def load_pickle(filepath: FilePath) -> Any:
    """Load a Python pickle file."""
    with open(filepath, "rb") as f:
        return pickle.load(f)


def save_pickle(obj: Any, filepath: FilePath) -> None:
    """Save an object to a Python pickle file."""
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "wb") as f:
        pickle.dump(obj, f)


def load_mat_file(filepath: FilePath) -> dict[str, Any]:
    """Load a MATLAB .mat file using mat73 first, then scipy.io fallback."""
    filepath = os.fspath(filepath)
    if mat73 is not None:
        try:
            return mat73.loadmat(filepath)
        except Exception:
            pass
    if scipy is None:
        raise ImportError("Install mat73 or scipy to load .mat files.")
    return scipy.io.loadmat(filepath, squeeze_me=True, struct_as_record=False)


def load_xarray_dataset(filepath: FilePath, **kwargs) -> xr.Dataset:
    """Load a NetCDF-compatible xarray dataset."""
    return xr.open_dataset(filepath, **kwargs)


def load_zarr_dataset(filepath: FilePath, **kwargs) -> xr.Dataset:
    """Load a Zarr dataset."""
    return xr.open_zarr(filepath, **kwargs)


def load_region_mask(filepath: FilePath) -> Any:
    """Load a region-mask dictionary, usually mask15n15_region_dict_*.pkl."""
    return load_pickle(filepath)


def load_eddy_pickles(data_root: FilePath, ce_name="cycvar.pickle", ae_name="acycvar.pickle"):
    """Load CE and AE eddy dictionaries from pickle files."""
    data_root = Path(data_root)
    ce = load_pickle(data_root / ce_name)
    ae = load_pickle(data_root / ae_name)
    return ce, ae


def prepare_eddy_dicts(
    ce,
    ae,
    donut_file: FilePath | None = None,
    start_date="2003-01-01",
    end_date="2022-12-31",
):
    """Add lifecycle, month, tstep8, and optional donut CHL to CE/AE eddy dictionaries."""
    ce = add_day_of_life(ce)
    ae = add_day_of_life(ae)
    for eddy in (ce, ae):
        eddy["month"] = argo_juld_to_month(eddy["time"])
        eddy["dayofyear"] = argo_juld_to_timestamp(eddy["time"]).dayofyear
        eddy["tstep8"] = argo_juld_to_tstep8(eddy["time"])

    if donut_file is not None:
        ce["chl_r_range_VHR"] = donut_velsshr_chl("CE", donut_file)
        ae["chl_r_range_VHR"] = donut_velsshr_chl("AE", donut_file)

    ce_subset = subset_dict_by_time(ce, start_date=start_date, end_date=end_date, time_key="time", time_converter=argo_juld_to_timestamp)
    ae_subset = subset_dict_by_time(ae, start_date=start_date, end_date=end_date, time_key="time", time_converter=argo_juld_to_timestamp)
    return ce_subset, ae_subset


def load_modis_chl(chl_dir: FilePath):
    """Load MODIS CHL Full_Time_Series.mat and Chl_fields.mat from the MATLAB workflow."""
    chl_dir = Path(chl_dir)
    chl_ts_mat = load_mat_file(chl_dir / "Full_Time_Series.mat")
    chl_field_mat = load_mat_file(chl_dir / "Chl_fields.mat")
    return {
        "chl_ts": chl_ts_mat["Chl"],
        "chl_loess": chl_field_mat.get("Chl_loess"),
        "Lat_c": chl_ts_mat["Lat"],
        "Lon_c": chl_ts_mat["Lon"],
        "raw": {"time_series": chl_ts_mat, "fields": chl_field_mat},
    }


def load_modis_k490_zeu(k490_dir: FilePath):
    """Load MODIS K490 MATLAB file and calculate KPAR/Zeu."""
    k490_dir = Path(k490_dir)
    k490_ts_mat = load_mat_file(k490_dir / "Full_Time_Series.mat")
    k490_ts = k490_ts_mat["K490"]
    zeu_ts, kd_par = k490_to_zeu(k490_ts)
    return {
        "k490_ts": k490_ts,
        "kd_par": kd_par,
        "zeu_ts": zeu_ts,
        "Lat_k": k490_ts_mat["Lat"],
        "Lon_k": k490_ts_mat["Lon"],
        "raw": k490_ts_mat,
    }


def load_mld_holte(mld_file: FilePath, lat_min=-50, lat_max=-20, lon_min=145, lon_max=175):
    """Load Holte Argo monthly MLD climatology and subset to the EAC region."""
    ds = xr.open_dataset(mld_file)
    ds_subset = ds.where((ds.lat >= lat_min) & (ds.lat <= lat_max) & (ds.lon >= lon_min) & (ds.lon <= lon_max), drop=True)
    return {"dataset": ds_subset, "mld_holte": ds_subset["mld_da_mean"].data}


def load_woa23_monthly(fn_dir: FilePath, lat_min=-50, lat_max=-20, lon_min=145, lon_max=175):
    """Load monthly WOA23 salinity, temperature, and nitrate fields into NumPy arrays."""
    fn_dir = Path(fn_dir)
    time_woa = np.arange(0, 12, 1) + 0.5
    sp_woa = t_woa = n_woa = None
    for i in range(1, 13):
        sp_fn = glob(str(fn_dir / "WOA2023_salinity_20052014" / f"woa23*_s{i:02d}_01.nc"))
        t_fn = glob(str(fn_dir / "WOA2023_temperature_20052014" / f"woa23*_t{i:02d}_01.nc"))
        n_fn = glob(str(fn_dir / "WOA2023_nitrate_all" / f"woa23*_n{i:02d}_01.nc"))
        if not sp_fn or not t_fn or not n_fn:
            raise FileNotFoundError(f"Missing WOA file for month {i:02d} under {fn_dir}")
        sp_ds = xr.open_dataset(sp_fn[0], decode_times=False)
        t_ds = xr.open_dataset(t_fn[0], decode_times=False)
        n_ds = xr.open_dataset(n_fn[0], decode_times=False)
        ix_int = np.where((sp_ds["lon"] >= lon_min) & (sp_ds["lon"] <= lon_max))[0]
        iy_int = np.where((sp_ds["lat"] >= lat_min) & (sp_ds["lat"] <= lat_max))[0]
        lon_ar = sp_ds.lon[ix_int].data
        lat_ar = sp_ds.lat[iy_int].data
        depth_ts_ar = sp_ds.depth.data
        depth_n_ar = n_ds.depth.data
        sp_ar = sp_ds.s_an[0, :, iy_int, ix_int].data[..., np.newaxis]
        t_ar = t_ds.t_an[0, :, iy_int, ix_int].data[..., np.newaxis]
        n_ar = n_ds.n_an[0, :, iy_int, ix_int].data[..., np.newaxis]
        if i == 1:
            sp_woa = sp_ar.copy()
            t_woa = t_ar.copy()
            n_woa = n_ar.copy()
            lat_woa = lat_ar.copy()
            lon_woa = lon_ar.copy()
            depth_ts_woa = depth_ts_ar.copy()
            depth_n_woa = depth_n_ar.copy()
        else:
            sp_woa = np.concatenate((sp_woa, sp_ar), axis=3)
            t_woa = np.concatenate((t_woa, t_ar), axis=3)
            n_woa = np.concatenate((n_woa, n_ar), axis=3)
    density = calc_woa_density_n2(sp_woa, t_woa, depth_ts_woa, lat_woa, lon_woa)
    return {
        "time_woa": time_woa,
        "SP_woa": sp_woa,
        "T_woa": t_woa,
        "N_woa": n_woa,
        "lat_woa": lat_woa,
        "lon_woa": lon_woa,
        "depth_TS_woa": depth_ts_woa,
        "depth_N_woa": depth_n_woa,
        **density,
    }
