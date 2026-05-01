"""Eddy-track helper functions for the EAC eddy chlorophyll-a workflow."""

from __future__ import annotations

import os
from typing import Any

import numpy as np
import pandas as pd

try:
    import mat73
except ImportError:  # pragma: no cover
    mat73 = None

try:
    from haversine import Unit, haversine
except ImportError:  # pragma: no cover
    Unit = None
    haversine = None


def add_day_of_life(eddy: dict[str, Any]) -> dict[str, Any]:
    """Add normalized day-of-life, stage, and total age to an eddy dictionary."""
    track_unique_ar = np.unique(eddy["track"])
    eddy["dayoflife_norm"] = np.nan * np.ones_like(eddy["track"], dtype=float)
    eddy["stage"] = np.nan * np.ones_like(eddy["track"], dtype=float)
    eddy["age"] = np.nan * np.ones_like(eddy["track"], dtype=float)

    for track in track_unique_ar:
        mask_set = eddy["track"] == track
        idx_mask = np.where(mask_set)[0]
        stage_ed = np.nan * np.ones_like(idx_mask, dtype=float)
        dol = eddy["observation_number"][mask_set]
        age = np.nanmax(dol)
        if not np.isfinite(age) or age == 0:
            continue
        dol_norm = dol / age
        if age <= 60:
            stage_ed[np.where(dol_norm <= 0.4)[0]] = 0
            stage_ed[np.where((dol_norm > 0.4) & (dol_norm <= 0.5))[0]] = 1
            stage_ed[np.where(dol_norm > 0.5)[0]] = 2
        else:
            stage_ed[np.where(dol_norm <= 0.25)[0]] = 0
            stage_ed[np.where((dol_norm > 0.25) & (dol_norm <= 0.75))[0]] = 1
            stage_ed[np.where(dol_norm > 0.75)[0]] = 2
        eddy["dayoflife_norm"][idx_mask] = dol_norm
        eddy["stage"][idx_mask] = stage_ed
        eddy["age"][idx_mask] = age
    return eddy

def calculate_track_top_n_mean(values, tracks, n: int = 5):
    """
    Calculate a track-level top-n mean and return it for each realization.

    For each eddy track, this calculates the mean of the n largest finite
    values, then assigns that track-level value back to all realizations
    belonging to the same track.

    Parameters
    ----------
    values
        Realization-level values, for example eddy amplitude.
    tracks
        Track ID for each realization.
    n
        Number of largest values to average per track.

    Returns
    -------
    numpy.ndarray
        Same length as values, containing the track-level top-n mean.
    """
    import numpy as np
    import pandas as pd

    values = np.asarray(values, dtype=float)
    tracks = np.asarray(tracks)

    df = pd.DataFrame(
        {
            "track": tracks,
            "value": values,
        }
    )

    def _top_n_mean(series):
        arr = series.to_numpy(dtype=float)
        arr = arr[np.isfinite(arr)]

        if arr.size == 0:
            return np.nan

        n_use = min(n, arr.size)
        return np.nanmean(np.sort(arr)[-n_use:])

    return (
        df.groupby("track")["value"]
        .transform(_top_n_mean)
        .to_numpy(dtype=float)
    )
def top_n_mean(arr, n: int = 5):
    """Return mean, count, and values for the top ``n`` finite values."""
    arr = np.asarray(arr)
    valid_vals = arr[np.isfinite(arr)]
    sorted_vals = np.sort(valid_vals)[::-1]
    count = len(sorted_vals)
    if count >= n:
        top_values = sorted_vals[:n]
        mean_value = np.mean(top_values)
    else:
        top_values = sorted_vals
        mean_value = np.nan
    return mean_value, count, top_values


def calculate_bearing(lon0, lat0, lon1, lat1):
    """Calculate bearing from point 0 to point 1 in degrees clockwise from north."""
    lon0_rad, lat0_rad = np.radians(lon0), np.radians(lat0)
    lon1_rad, lat1_rad = np.radians(lon1), np.radians(lat1)
    dlon = lon1_rad - lon0_rad
    x = np.sin(dlon) * np.cos(lat1_rad)
    y = np.cos(lat0_rad) * np.sin(lat1_rad) - np.sin(lat0_rad) * np.cos(lat1_rad) * np.cos(dlon)
    return (np.degrees(np.arctan2(x, y)) + 360) % 360


def calculate_distance_speed_direction(eddy: dict[str, Any]):
    """Calculate track distance, direction, speed, top amplitude, radius, and birthday."""
    if haversine is None or Unit is None:
        raise ImportError("Install haversine to use calculate_distance_speed_direction().")

    track_ar = eddy["track"]
    track_unique = np.unique(track_ar)
    dist_accum_ar = np.full_like(track_ar, np.nan, dtype=np.float32)
    dist_vector_ar = np.full_like(track_ar, np.nan, dtype=np.float32)
    dist_track_ar = np.full_like(track_ar, np.nan, dtype=np.float32)
    dir_track_ar = np.full_like(track_ar, np.nan, dtype=np.float32)
    vector_track_ar = np.full_like(track_ar, np.nan, dtype=np.float32)
    speed_envel_ar = np.full_like(track_ar, np.nan, dtype=np.float32)
    speed_vector_ar = np.full_like(track_ar, np.nan, dtype=np.float32)
    dir_envel_ar = np.full_like(track_ar, np.nan, dtype=np.float32)
    amp_track_ar = np.full_like(track_ar, np.nan, dtype=np.float32)
    rad_track_ar = np.full_like(track_ar, np.nan, dtype=np.float32)
    birthday_ar = np.full_like(track_ar, np.nan, dtype=np.float32)

    for track in track_unique:
        track_index = np.where(track_ar == track)[0]
        if len(track_index) == 0:
            continue
        lon_ar = (eddy["longitude"][track_index] + 180) % 360 - 180
        lat_ar = eddy["latitude"][track_index]
        obs_ar = eddy["observation_number"][track_index]
        amp_realz = eddy["amplitude"][track_index]
        amp_track_ar[track_index] = top_n_mean(amp_realz, n=5)[0]
        rad_realz = eddy["speed_radius"][track_index] / 1e3
        rad_track_ar[track_index] = top_n_mean(rad_realz, n=5)[0]
        dist_accum_ar[track_index[0]] = 0
        speed_envel_ar[track_index[0]] = 0
        dir_envel_ar[track_index[0]] = 0
        birthday_ar[track_index] = eddy["dayofyear"][track_index[np.argmin(obs_ar)]]
        for i in range(1, len(track_index)):
            dist_vector_ar[track_index[i]] = haversine(
                (lat_ar[0], lon_ar[0]), (lat_ar[i], lon_ar[i]), unit=Unit.KILOMETERS
            )
            speed_vector_ar[track_index[i]] = dist_vector_ar[track_index[i]] / (obs_ar[i] - obs_ar[0])
            dist_realz = haversine(
                (lat_ar[i - 1], lon_ar[i - 1]), (lat_ar[i], lon_ar[i]), unit=Unit.KILOMETERS
            )
            speed_envel_ar[track_index[i]] = dist_realz / (obs_ar[i] - obs_ar[i - 1])
            dir_envel_ar[track_index[i]] = calculate_bearing(lon_ar[i - 1], lat_ar[i - 1], lon_ar[i], lat_ar[i])
            dist_accum_ar[track_index[i]] = dist_accum_ar[track_index[i - 1]] + dist_realz
        last_i = len(track_index) - 1
        dist_track_ar[track_index] = dist_accum_ar[track_index[last_i]]
        vector_track_ar[track_index] = haversine(
            (lat_ar[0], lon_ar[0]), (lat_ar[last_i], lon_ar[last_i]), unit=Unit.KILOMETERS
        )
        if np.any(np.isfinite(dist_vector_ar[track_index])):
            ii = np.nanargmax(dist_vector_ar[track_index])
            dir_track_ar[track_index] = calculate_bearing(lon_ar[0], lat_ar[0], lon_ar[ii], lat_ar[ii])

    return (
        dist_track_ar,
        dir_track_ar,
        dist_accum_ar,
        vector_track_ar,
        dist_vector_ar,
        speed_envel_ar,
        dir_envel_ar,
        speed_vector_ar,
        amp_track_ar,
        rad_track_ar,
        birthday_ar,
    )


def chl_from_donut_struct(donut_struct: dict[str, Any], polarity: str):
    """Calculate core/velocity/extended CHL anomaly summaries from donut_struct."""
    if polarity in {"C", "CE", "ce"}:
        key = "CE"
    elif polarity in {"A", "AE", "ae"}:
        key = "AE"
    else:
        key = polarity

    ds = donut_struct[key]
    nobs_ivel = ds["Nobs_r_range"][:, :5]
    nobs_dr = ds["Nobs_dr_velssh"]
    nobs_core = ds["Nobs_r_range"][:, :2]
    chl_ivel = ds["chlanom_r_range"][:, :5].copy()
    chl_dr = ds["chlanom_dr_velssh"].copy()
    chl_core = ds["chlanom_r_range"][:, :2].copy()
    chl_ivel[np.isnan(chl_ivel)] = 0
    chl_dr[np.isnan(chl_dr)] = 0
    chl_core[np.isnan(chl_core)] = 0
    chl_r0 = np.nansum(chl_core * nobs_core, axis=1) / np.nansum(nobs_core, axis=1)
    chl_r1 = np.nansum(chl_ivel * nobs_ivel, axis=1) / np.nansum(nobs_ivel, axis=1)
    chl_r2 = ((chl_dr * nobs_dr) + np.nansum(chl_ivel * nobs_ivel, axis=1)) / (
        np.nansum(nobs_ivel, axis=1) + nobs_dr
    )
    return np.column_stack([chl_r0, chl_r1, chl_r2])


def donut_velsshr_chl(polarity: str, donut_file: str | os.PathLike[str]):
    """Load donut_struct from a .mat file and calculate CHL anomaly summaries."""
    if mat73 is None:
        raise ImportError("Install mat73 to read MATLAB v7.3 .mat files.")
    ds_mat = mat73.loadmat(donut_file)
    return chl_from_donut_struct(ds_mat["donut_struct"], polarity)


def subset_dict_by_time(data_dict, start_date, end_date, time_key="time", time_converter=None):
    """Subset a dictionary of arrays/Series/Index by a time range."""
    raw_time = data_dict[time_key]
    times = time_converter(raw_time) if time_converter is not None else pd.to_datetime(raw_time)
    mask = (times >= pd.to_datetime(start_date)) & (times <= pd.to_datetime(end_date))
    subset = {}
    for key, val in data_dict.items():
        if isinstance(val, np.ndarray) and val.shape and val.shape[0] == mask.size:
            subset[key] = val[mask]
        elif isinstance(val, np.ndarray) and val.ndim >= 2 and val.shape[1] == mask.size:
            subset[key] = val[:, mask]
        elif isinstance(val, pd.Series) and len(val) == mask.size:
            subset[key] = val[mask].reset_index(drop=True)
        elif isinstance(val, pd.Index) and len(val) == mask.size:
            subset[key] = val[mask]
        else:
            subset[key] = val
    return subset


# Backward-compatible names from the notebook.
addDayOfLife = add_day_of_life
Donut_VelSSHR_Chl = donut_velsshr_chl
