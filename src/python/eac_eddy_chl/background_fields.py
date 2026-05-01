"""Background-field calculations for WOA, MLD, K490/Zeu, and stratification."""

from __future__ import annotations

import numpy as np

try:
    import gsw
except ImportError:  # pragma: no cover
    gsw = None


def k490_to_zeu(k490):
    """Convert MODIS K490 to KPAR and euphotic depth using the paper workflow."""
    k490 = np.asarray(k490, dtype=float)
    kd_par = 0.0665 + 0.874 * k490 - 0.00121 / k490
    zeu = 4.6 / kd_par
    return zeu, kd_par


def calc_woa_density_n2(sp_woa, t_woa, depth_ts_woa, lat_woa, lon_woa, g=9.81, rho0=1025):
    """Calculate SA, CT, sigma0, density, and N2 from WOA salinity/temperature fields."""
    if gsw is None:
        raise ImportError("Install gsw to calculate density and N2 from WOA fields.")
    pressure_woa = gsw.p_from_z(-depth_ts_woa[:, None], lat_woa[None, :])
    sa_woa = gsw.SA_from_SP(
        sp_woa, pressure_woa[:, :, None, None], lon_woa[None, None, :, None], lat_woa[None, :, None, None]
    )
    ct_woa = gsw.CT_from_t(sa_woa, t_woa, pressure_woa[:, :, None, None])
    sigma0_woa = gsw.sigma0(sa_woa, ct_woa)
    rho = 1000 + sigma0_woa
    drho_dz = np.gradient(rho, depth_ts_woa, axis=0)
    n2_woa = g / rho0 * drho_dz
    return {
        "pressure_woa": pressure_woa,
        "SA_woa": sa_woa,
        "CT_woa": ct_woa,
        "Sigma0_woa": sigma0_woa,
        "rho_woa": rho,
        "N2_woa": n2_woa,
    }


def calc_energy_diagnostics(rho, n2_woa, depth_ts_woa, g=9.81):
    """Calculate MPE, EPE, and depth-integrated N2 diagnostics from density/N2."""
    idx_500 = depth_ts_woa <= 500
    z = depth_ts_woa.copy()
    dz = np.gradient(depth_ts_woa)
    mpe = np.sum(rho[idx_500, :, :, :] * g * z[idx_500, None, None, None] * dz[idx_500, None, None, None], axis=0)
    rho_mean = np.nanmean(rho, axis=3, keepdims=True)
    rho_anom = rho - rho_mean
    n2_safe = np.where(n2_woa <= 0, np.nan, n2_woa)
    epe = 0.5 * g**2 / n2_safe * rho_anom**2
    epe_integrated = np.nansum(epe * dz[:, None, None, None], axis=0)
    n2_integrated = np.nanmean(n2_woa * dz[:, None, None, None], axis=0)
    return {"MPE": mpe, "EPE": epe, "EPE_integrated": epe_integrated, "N2_integrated": n2_integrated}
