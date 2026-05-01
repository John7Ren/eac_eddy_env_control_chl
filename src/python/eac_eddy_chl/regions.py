"""Region and background-climatology helper functions."""

from __future__ import annotations

import numpy as np
from scipy.interpolate import interp1d

from .grid_utils import make_2deg_bins


def mask_from_region_grid(region_mask, lat, lon, lat_bins=None, lon_bins=None):
    """Return a boolean lat/lon mask by expanding a 2°x2° region mask onto a dataset grid."""
    if lat_bins is None or lon_bins is None:
        lat_bins, lon_bins = make_2deg_bins()
    LAT, LON = np.meshgrid(lat, lon, indexing="ij")
    mask_loc = np.zeros_like(LAT, dtype=bool)
    y_ar, x_ar = np.where(region_mask)
    for iy, ix in zip(y_ar, x_ar):
        ilat = lat_bins[iy]
        ilat_1 = lat_bins[iy + 1]
        ilon = lon_bins[ix]
        ilon_1 = lon_bins[ix + 1]
        mask_ilat = (LAT >= ilat_1) & (LAT < ilat)
        mask_ilon = (LON >= ilon) & (LON < ilon_1)
        mask_loc |= mask_ilat & mask_ilon
    return mask_loc


def return_index_dstype_region(
    region_name,
    ds_type,
    region_masks,
    *,
    lat_woa=None,
    lon_woa=None,
    depth_woa=None,
    lat_wd=None,
    lon_wd=None,
    depth_max=800,
):
    """Return vertical, y, x indices for a named region and dataset type."""
    region_mask = region_masks[region_name]
    if ds_type == "woa":
        if lat_woa is None or lon_woa is None or depth_woa is None:
            raise ValueError("lat_woa, lon_woa, and depth_woa are required for ds_type='woa'.")
        lat = lat_woa
        lon = lon_woa
        depth = depth_woa
    elif ds_type == "wind_stress":
        if lat_wd is None or lon_wd is None:
            raise ValueError("lat_wd and lon_wd are required for ds_type='wind_stress'.")
        lat = lat_wd
        lon = lon_wd
        depth = None
    else:
        raise ValueError("ds_type must be 'woa' or 'wind_stress'.")

    mask_loc = mask_from_region_grid(region_mask, lat, lon)
    iy_ar, ix_ar = np.where(mask_loc)
    if ds_type == "wind_stress":
        id_ar = 0
    else:
        id_ar = np.where(depth <= depth_max)[0]
    return id_ar, iy_ar, ix_ar


def max_depth_below_threshold(data, depths, threshold=1, max_depth=200):
    """Return deepest depth where data < threshold and depth <= max_depth."""
    data = np.asarray(data)
    depths = np.asarray(depths)
    depths_expand = np.expand_dims(depths, axis=tuple(range(1, len(data.shape))))
    mask = (data < threshold) & (depths_expand <= max_depth)
    index_matrix = np.expand_dims(np.arange(len(depths)), axis=tuple(range(1, len(data.shape))))
    valid_indices = np.where(mask, index_matrix, -1)
    max_indices = valid_indices.max(axis=0).astype(int)
    return np.where(max_indices >= 0, depths[max_indices], 0.0)


def interpolate_cyclic_time_series(time, data, t_interp, cycle_length=360):
    """Cubic interpolation for a cyclic 1D climatological time series."""
    time = np.asarray(time)
    data = np.asarray(data)
    time_extended = np.append(time, time[0:3] + cycle_length)
    data_extended = np.append(data, data[0:3])
    interp_func = interp1d(time_extended, data_extended, kind="cubic", bounds_error=False, fill_value="extrapolate")
    return interp_func(t_interp)


def interpolate_N2(time, data, t_interp, cycle_length=360):
    """Cyclic cubic interpolation for a 2D array with shape (depth, month)."""
    time = np.asarray(time)
    data = np.asarray(data)
    time_extended = np.append(time, time[0:3] + cycle_length)
    data_extended = np.hstack((data, data[:, [0, 1, 2]]))
    data_interp = np.full((data_extended.shape[0], len(t_interp)), np.nan)
    for i in range(data_extended.shape[0]):
        idata = data_extended[i, :]
        nanmask = np.isnan(idata)
        if (~nanmask).sum() > 2:
            interp_func = interp1d(
                time_extended[~nanmask], idata[~nanmask], kind="cubic", bounds_error=False, fill_value="extrapolate"
            )
            data_interp[i, :] = interp_func(t_interp)
    return data_interp


def find_crossings(t, y1, y2):
    """Find times/values where y1 crosses y2."""
    d = y1 - y2
    crossings = []
    for i in range(len(d) - 1):
        if d[i] == 0:
            crossings.append((t[i], y1[i]))
        elif d[i] * d[i + 1] < 0:
            frac = -d[i] / (d[i + 1] - d[i])
            t_cross = t[i] + frac * (t[i + 1] - t[i])
            y_cross = y1[i] + frac * (y1[i + 1] - y1[i])
            crossings.append((t_cross, y_cross))
    return crossings
