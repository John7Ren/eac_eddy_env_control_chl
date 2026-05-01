"""Time conversion utilities for the EAC eddy chlorophyll-a workflow."""

from __future__ import annotations

import numpy as np
import pandas as pd


def argo_juld_to_timestamp(argo_juld):
    """Convert Argo-style days since 1950-01-01 to pandas timestamps."""
    day0_ts = pd.Timestamp("1950-01-01T00:00:00")
    day0_juld = pd.Timestamp.to_julian_date(day0_ts)
    return pd.to_datetime(np.asarray(argo_juld) + day0_juld, origin="julian", unit="D")


def argo_juld_to_dayofyear(argo_juld):
    """Convert Argo-style days since 1950-01-01 to day-of-year."""
    return argo_juld_to_timestamp(argo_juld).dayofyear


def argo_juld_to_year(argo_juld):
    """Convert Argo-style days since 1950-01-01 to calendar year."""
    return argo_juld_to_timestamp(argo_juld).year


def argo_juld_to_month(argo_juld):
    """Convert Argo-style days since 1950-01-01 to calendar month."""
    return argo_juld_to_timestamp(argo_juld).month.astype(float)


def argo_juld_to_season(argo_juld):
    """Return Australian seasons: 1 summer, 2 autumn, 3 winter, 4 spring."""
    month = argo_juld_to_timestamp(argo_juld).month
    season = np.zeros_like(month, dtype=float)
    season_groups = np.array([[12, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]])
    for i, months in enumerate(season_groups):
        season[np.isin(month, months)] = i + 1
    return season


def get_8day_bins():
    """Return 8-day day-of-year bin edges and bin labels used in Paper 1."""
    tbins_8day = np.arange(1, 362, 8)
    tbins_8day = np.append(tbins_8day, 367)
    tsteps = tbins_8day[:-1] + 4
    return tbins_8day, tsteps


def dayofyear_to_tstep8(dayofyear):
    """Map day-of-year values to the 8-day climatological time-step labels."""
    dayofyear = np.asarray(dayofyear)
    tbins_8day, tlabels = get_8day_bins()
    tstep8 = np.full_like(dayofyear, np.nan, dtype=float)
    for i, dt in enumerate(tlabels):
        st = tbins_8day[i]
        ed = tbins_8day[i + 1]
        tstep8[np.isin(dayofyear, np.arange(st, ed, 1))] = dt
    return tstep8.astype(float)


def argo_juld_to_tstep8(argo_juld):
    """Convert Argo-style days since 1950-01-01 to 8-day climatological labels."""
    return dayofyear_to_tstep8(argo_juld_to_dayofyear(argo_juld))


def julian_to_tstep8(juld):
    """Convert absolute Julian day values to 8-day climatological labels."""
    dayofyear = pd.to_datetime(juld, origin="julian", unit="D").dayofyear
    return dayofyear_to_tstep8(dayofyear)


# Backward-compatible names from the original notebook.
argoJULD2dayofyear = argo_juld_to_dayofyear
argoJULD2year = argo_juld_to_year
argoJULD2TS = argo_juld_to_timestamp
argoJULD2season = argo_juld_to_season
argoJULD2month = argo_juld_to_month
argoJULD2tstep8 = argo_juld_to_tstep8
JULD2tstep8 = julian_to_tstep8
