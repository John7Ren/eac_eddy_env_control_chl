"""Climatology builders for Paper 1 diagnostics."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

from .grid_utils import make_2deg_bins, resample_spatial
from .regions import (
    interpolate_N2,
    interpolate_cyclic_time_series,
    mask_from_region_grid,
    max_depth_below_threshold,
    return_index_dstype_region,
)
from .time_utils import get_8day_bins

def _smooth_1d_nan_safe(y, window_length=5, polyorder=2):
    """
    Smooth a 1D seasonal series after linearly filling internal NaNs.

    Returns all-NaN if the input has no finite values.
    Returns the original/fill-interpolated series if there are too few points for Savitzky-Golay.
    """
    y = np.asarray(y, dtype=float)
    out = np.full_like(y, np.nan, dtype=float)

    finite = np.isfinite(y)
    if finite.sum() == 0:
        return out

    x = np.arange(y.size)

    # Fill NaNs by linear interpolation so savgol_filter does not propagate NaNs.
    y_fill = y.copy()
    if finite.sum() == 1:
        y_fill[:] = y[finite][0]
    else:
        y_fill[~finite] = np.interp(x[~finite], x[finite], y[finite])

    # Savitzky-Golay requires an odd window and polyorder < window_length.
    if y.size < window_length or finite.sum() <= polyorder:
        return y_fill

    return savgol_filter(y_fill, window_length=window_length, polyorder=polyorder)

def build_background_map_climatology(
    *,
    chl_ts,
    zeu_ts,
    mld_holte,
    N_woa,
    N2_woa,
    lat_woa,
    lon_woa,
    depth_N_woa,
    Lat_c,
    Lon_c,
    Lat_k,
    Lon_k,
    lat_bins=None,
    lon_bins=None,
):
    """Build 2°x2° background metrics for CHL, Zeu, MLD, N2, and nitracline depth."""
    if lat_bins is None or lon_bins is None:
        lat_bins, lon_bins = make_2deg_bins()
    lx = len(lon_bins) - 1
    ly = len(lat_bins) - 1
    _, tsteps = get_8day_bins()
    m_ar = np.arange(1, 13)

    climt_bg_dict = {
        "Chl": {k: np.full((ly, lx), np.nan) for k in ["timing", "mag", "mean", "std"]},
        "MLD": {k: np.full((ly, lx), np.nan) for k in ["timing", "mag", "mean", "std"]},
        "N2": {k: np.full((ly, lx), np.nan) for k in ["timing", "depth", "value"]},
        "Zeu": {k: np.full((ly, lx), np.nan) for k in ["timing", "mag", "mean", "std"]},
        "Ncline": {k: np.full((ly, lx), np.nan) for k in ["mean", "std"]},
    }
    climt2d_bg_dict = {"N2": {"maxvalue": np.full((12, ly, lx), np.nan), "maxdepth": np.full((12, ly, lx), np.nan)}}

    LAT_WOA, LON_WOA = np.meshgrid(lat_woa, lon_woa, indexing="ij")
    for iy, ilat in enumerate(lat_bins[:-1]):
        for ix, ilon in enumerate(lon_bins[:-1]):
            ilat_1 = lat_bins[iy + 1]
            ilon_1 = lon_bins[ix + 1]
            mask_LOC = (LAT_WOA >= ilat_1) & (LAT_WOA < ilat) & (LON_WOA >= ilon) & (LON_WOA < ilon_1)
            iy_ar_woa, ix_ar_woa = np.where(mask_LOC)
            id_ar_woa = np.where(depth_N_woa <= 800)[0]
            depth_ar_woa = depth_N_woa[id_ar_woa.squeeze()]
            id_broad_N = id_ar_woa[:, None, None, None]
            iy_broad_N = iy_ar_woa[None, :, None, None]
            ix_broad_N = ix_ar_woa[None, None, :, None]
            incline = np.full(12, np.nan)
            in2_depth = np.full(12, np.nan)
            in2_value = np.full(12, np.nan)
            for it in range(12):
                iiN_woa = N_woa[id_broad_N, iy_broad_N, ix_broad_N, it]
                if not np.isnan(np.nanmean(iiN_woa)):
                    incline[it] = np.nanmean(max_depth_below_threshold(iiN_woa, depth_ar_woa))
                iiN2_woa = np.nanmean(N2_woa[id_broad_N, iy_broad_N, ix_broad_N, it], axis=(1, 2, 3))
                if np.any(np.isfinite(iiN2_woa)):
                    iidx_N2 = np.nanargmax(iiN2_woa)
                    in2_depth[it] = depth_N_woa[iidx_N2]
                    in2_value[it] = iiN2_woa[iidx_N2]
            climt_bg_dict["Ncline"]["mean"][iy, ix] = np.nanmean(incline)
            climt_bg_dict["Ncline"]["std"][iy, ix] = np.nanstd(incline)
            climt2d_bg_dict["N2"]["maxvalue"][:, iy, ix] = in2_value.copy()
            climt2d_bg_dict["N2"]["maxdepth"][:, iy, ix] = in2_depth.copy()

            if np.any(np.isfinite(in2_value)):
                idx_N2 = np.nanargmax(in2_value)
                climt_bg_dict["N2"]["timing"][iy, ix] = m_ar[idx_N2]
                climt_bg_dict["N2"]["depth"][iy, ix] = in2_depth[idx_N2]
                climt_bg_dict["N2"]["value"][iy, ix] = in2_value[idx_N2]

            imld_climt = np.nanmean(mld_holte[mask_LOC, :], axis=0)
            if np.any(np.isfinite(imld_climt)):
                ipk = np.nanargmax(imld_climt)
                climt_bg_dict["MLD"]["timing"][iy, ix] = m_ar[ipk]
                climt_bg_dict["MLD"]["mag"][iy, ix] = imld_climt[ipk]
                climt_bg_dict["MLD"]["mean"][iy, ix] = np.nanmean(imld_climt)
                climt_bg_dict["MLD"]["std"][iy, ix] = np.nanstd(imld_climt)

            mask_LOC_k = (Lat_k >= ilat_1) & (Lat_k < ilat) & (Lon_k >= ilon) & (Lon_k < ilon_1)
            izeu_climt = np.nanmean(zeu_ts[:, :, mask_LOC_k], axis=(0, 2))
            if np.any(np.isfinite(izeu_climt)):
                izeu_climt_smooth = _smooth_1d_nan_safe(izeu_climt, window_length=5, polyorder=2)

                if np.any(np.isfinite(izeu_climt_smooth)):
                    ipk = np.nanargmin(izeu_climt_smooth)
                    climt_bg_dict["Zeu"]["mag"][iy, ix] = izeu_climt_smooth[ipk]
                    climt_bg_dict["Zeu"]["mean"][iy, ix] = np.nanmean(izeu_climt)
                    climt_bg_dict["Zeu"]["std"][iy, ix] = np.nanstd(izeu_climt)
                    climt_bg_dict["Zeu"]["timing"][iy, ix] = pd.to_datetime(
                        f"2001-{tsteps[ipk]:03d}", format="%Y-%j"
                    ).month

            mask_LOC_c = (Lat_c >= ilat_1) & (Lat_c < ilat) & (Lon_c >= ilon) & (Lon_c < ilon_1)
            ichl_climt = np.nanmean(chl_ts[:, :, mask_LOC_c], axis=(0, 2))
            if np.any(np.isfinite(ichl_climt)):
                ichl_climt_smooth = _smooth_1d_nan_safe(ichl_climt, window_length=5, polyorder=2)

                if np.any(np.isfinite(ichl_climt_smooth)):
                    ipk = np.nanargmax(ichl_climt_smooth)
                    climt_bg_dict["Chl"]["mag"][iy, ix] = ichl_climt_smooth[ipk]
                    climt_bg_dict["Chl"]["mean"][iy, ix] = np.nanmean(ichl_climt)
                    climt_bg_dict["Chl"]["std"][iy, ix] = np.nanstd(ichl_climt)
                    climt_bg_dict["Chl"]["timing"][iy, ix] = pd.to_datetime(
                        f"2001-{tsteps[ipk]:03d}", format="%Y-%j"
                    ).month
    return climt_bg_dict, climt2d_bg_dict


def build_depth_climatology_for_region(
    *,
    region_name,
    region_masks,
    ce_subset,
    ae_subset,
    chl_ts,
    zeu_ts,
    N_woa,
    N2_woa,
    mld_holte,
    time_woa,
    depth_TS_woa,
    depth_N_woa,
    lat_woa,
    lon_woa,
    yrls=None,
    lat_bins=None,
    lon_bins=None,
):
    """Build regional 8-day depth/environment and eddy CHL anomaly climatologies."""
    if lat_bins is None or lon_bins is None:
        lat_bins, lon_bins = make_2deg_bins()
    lx = len(lon_bins) - 1
    ly = len(lat_bins) - 1
    tbins_8day, tsteps = get_8day_bins()
    lt = len(tsteps)
    if yrls is None:
        yrls = np.arange(2002, 2023)

    id_ar_woa, iy_ar_woa, ix_ar_woa = return_index_dstype_region(
        region_name,
        "woa",
        region_masks,
        lat_woa=lat_woa,
        lon_woa=lon_woa,
        depth_woa=depth_TS_woa,
    )
    id_broad_woa = id_ar_woa[:, None, None, None]
    iy_broad_woa = iy_ar_woa[None, :, None, None]
    ix_broad_woa = ix_ar_woa[None, None, :, None]
    iy_broad_mld_woa = iy_ar_woa[:, None, None]
    ix_broad_mld_woa = ix_ar_woa[None, :, None]

    ncline_ave_woa = np.full(12, np.nan)
    mld_ave_holte = np.full(12, np.nan)
    n_contour_ls = np.arange(0.1, 1.0, 0.1)
    n_matrix = np.full((len(n_contour_ls), 12), np.nan)
    N2_ave_woa = np.full((len(id_ar_woa.squeeze()), 12), np.nan)
    depth_ar_woa = depth_TS_woa[id_ar_woa.squeeze()]

    for it in range(12):
        idata_woa = N_woa[id_broad_woa, iy_broad_woa, ix_broad_woa, it]
        ncline_ave_woa[it] = np.nanmean(max_depth_below_threshold(idata_woa, depth_ar_woa, threshold=1))
        for iiN, iN in enumerate(n_contour_ls):
            n_matrix[iiN, it] = np.nanmean(max_depth_below_threshold(idata_woa, depth_ar_woa, threshold=iN))
        mld_ave_holte[it] = np.nanmean(mld_holte[iy_broad_mld_woa, ix_broad_mld_woa, it], axis=(0, 1))
        N2_ave_woa[:, it] = np.nanmean(N2_woa[id_broad_woa, iy_broad_woa, ix_broad_woa, it], axis=(1, 2)).squeeze()

    time_inday = time_woa * 30
    ncline_ave_woa_8d = interpolate_cyclic_time_series(time_inday, ncline_ave_woa, tsteps)
    mld_ave_holte_8d = interpolate_cyclic_time_series(time_inday, mld_ave_holte, tsteps)
    N2_ave_woa_8d = interpolate_N2(time_inday, N2_ave_woa, tsteps)
    n_matrix_8d = interpolate_N2(time_inday, n_matrix, tsteps)

    chl_ts_resampled = np.full((len(yrls), lt, ly, lx), np.nan)
    zeu_ts_resampled = np.full((len(yrls), lt, ly, lx), np.nan)
    for iyr, _ in enumerate(yrls):
        chl_ts_resampled[iyr, :, :, :] = resample_spatial(np.squeeze(chl_ts[iyr, :, :, :]), ly, lx)
        zeu_ts_resampled[iyr, :, :, :] = resample_spatial(np.squeeze(zeu_ts[iyr, :, :, :]), ly, lx)

    depth_climt_dict = {
        "chl_bg": np.full(lt, np.nan),
        "zeu_bg": np.full(lt, np.nan),
        "ncline_bg": ncline_ave_woa_8d.copy(),
        "mld_holte": mld_ave_holte_8d.copy(),
        "N2_woa_2d": N2_ave_woa_8d.copy(),
        "n_matrix_8d": n_matrix_8d.copy(),
        "chl_ce": np.full(lt, np.nan),
        "chl_ae": np.full(lt, np.nan),
    }

    region_mask = region_masks[region_name]
    region_mask_tile = np.tile(region_mask, (len(yrls), 1, 1))
    for ist in range(lt):
        depth_climt_dict["chl_bg"][ist] = np.nanmean(chl_ts_resampled[:, ist, ::-1, :][region_mask_tile])
        depth_climt_dict["zeu_bg"][ist] = np.nanmean(zeu_ts_resampled[:, ist, ::-1, :][region_mask_tile])

    eddy_climt = build_region_eddy_climatology(region_name, region_masks, ce_subset, ae_subset, lat_bins=lat_bins, lon_bins=lon_bins)
    depth_climt_dict["chl_ce"] = eddy_climt["chl_anom_climt_dict"]["C"].copy()
    depth_climt_dict["chl_ae"] = eddy_climt["chl_anom_climt_dict"]["A"].copy()
    return depth_climt_dict, eddy_climt


def eddy_region_mask(eddy, region_mask, lat_bins=None, lon_bins=None):
    """Return boolean mask selecting eddy realizations inside a 2°x2° region mask."""
    if lat_bins is None or lon_bins is None:
        lat_bins, lon_bins = make_2deg_bins()
    ilat_ar = eddy["latitude"]
    ilon_ar = eddy["longitude"]
    y_ar, x_ar = np.where(region_mask)
    mask_loc = np.zeros_like(ilat_ar, dtype=bool)
    for iy, ix in zip(y_ar, x_ar):
        ilat = lat_bins[iy]
        ilat_1 = lat_bins[iy + 1]
        ilon = lon_bins[ix]
        ilon_1 = lon_bins[ix + 1]
        mask_loc |= (ilat_ar >= ilat_1) & (ilat_ar < ilat) & (ilon_ar >= ilon) & (ilon_ar < ilon_1)
    return mask_loc


def build_region_eddy_climatology(region_name, region_masks, ce_subset, ae_subset, lat_bins=None, lon_bins=None):
    """Build 8-day eddy CHL anomaly climatology and metadata arrays for a region."""
    if lat_bins is None or lon_bins is None:
        lat_bins, lon_bins = make_2deg_bins()
    tbins_8day, tsteps = get_8day_bins()
    lt = len(tsteps)
    out = {
        "tstep_realz_dict": {},
        "chl_anom_dict": {},
        "amp_realz_dict": {},
        "amp_track_dict": {},
        "age_realz_dict": {},
        "birthday_track_dict": {},
        "doy_realz_dict": {},
        "age_norm_dict": {},
        "life_span_dict": {},
        "chl_anom_climt_dict": {"C": np.full(lt, np.nan), "A": np.full(lt, np.nan)},
    }
    region_mask = region_masks[region_name]
    for ie, eddy in {"C": ce_subset, "A": ae_subset}.items():
        mask_loc = eddy_region_mask(eddy, region_mask, lat_bins=lat_bins, lon_bins=lon_bins)
        doy = eddy["dayofyear"]
        out["tstep_realz_dict"][ie] = eddy["tstep8"][mask_loc]
        out["chl_anom_dict"][ie] = eddy["chl_r_range_VHR"][mask_loc, 1]
        out["amp_realz_dict"][ie] = eddy["amplitude"][mask_loc]
        out["amp_track_dict"][ie] = eddy.get("amp_track", np.full_like(eddy["amplitude"], np.nan))[mask_loc]
        out["age_realz_dict"][ie] = eddy["observation_number"][mask_loc]
        out["birthday_track_dict"][ie] = eddy.get("birthday_track", np.full_like(eddy["dayofyear"], np.nan))[mask_loc]
        out["doy_realz_dict"][ie] = eddy["dayofyear"][mask_loc]
        out["age_norm_dict"][ie] = eddy["dayoflife_norm"][mask_loc]
        out["life_span_dict"][ie] = eddy["age"][mask_loc]
        for ist in range(lt):
            st = tbins_8day[ist]
            ed = tbins_8day[ist + 1]
            mask_doy = (doy >= st) & (doy < ed)
            out["chl_anom_climt_dict"][ie][ist] = np.nanmedian(eddy["chl_r_range_VHR"][mask_loc & mask_doy, 1])
    return out
