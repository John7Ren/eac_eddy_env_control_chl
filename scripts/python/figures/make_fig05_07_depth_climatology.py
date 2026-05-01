"""
Make paper-style Figures 5–7: regional depth climatology and amplitude-group
eddy-induced chlorophyll-a anomalies.

This script reproduces the structure of the paper plotting code:
(a) stratification and key depths
(b) background CHL plus eddy-induced anomaly
(c) eddy-induced CHL anomaly grouped by amplitude quartiles
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter
from scipy.stats import wilcoxon

from eac_eddy_chl.io import load_pickle
from eac_eddy_chl.regions import find_crossings


# =============================================================================
# Paths
# =============================================================================

DATA_ROOT = Path(
    "/Users/pearl/Library/CloudStorage/OneDrive-UniversityofTasmania/Work/"
    "eac_eddy_env_control_chl_data"
)

PROCESSED_DIR = DATA_ROOT / "processed"
SAVE_DIR = DATA_ROOT / "figures" / "paper_style"
SAVE_DIR.mkdir(parents=True, exist_ok=True)

CLIM_FILE = PROCESSED_DIR / "step04_climatology_outputs.pkl"
WOA_FILE = PROCESSED_DIR / "step02_woa_background_bundle.pkl"


# =============================================================================
# Plot settings
# =============================================================================

lbsz = 28
tksz = 24
lgsz = 18
lw = 4

y_tt = 1.05
x_tt = 0.0

percentile_ls = [0, 25, 50, 75, 100]

c_dict = {
    "C": ["tab:cyan", "tab:blue", "midnightblue", "Purple"],
    "A": ["moccasin", "tab:orange", "tab:red", "maroon"],
}

REGION_LABELS = {
    "ED": "ED",
    "Open": "TS",
    "Low-N": "NL",
}

FIGURE_NUMBERS = {
    "ED": 5,
    "Open": 6,
    "Low-N": 7,
}

# You can manually tune these later to exactly match the submitted paper.
REGION_LIMITS = {
    "ED": {
        "zmin": -0.12,
        "zmax": 0.08,
        "bg_max": 0.60,
        "regime_t0": None,
    },
    "Open": {
        "zmin": -0.04,
        "zmax": 0.03,
        "bg_max": 0.35,
        "regime_t0": None,
    },
    "Low-N": {
        "zmin": -0.015,
        "zmax": 0.015,
        "bg_max": 0.18,
        # Change this if your paper used a different vertical boundary.
        "regime_t0": 181,
    },
}


# 46 8-day bins: centers are 5, 13, ..., 365.
tbins_8day = np.append(np.arange(1, 362, 8), 367)
tsteps = tbins_8day[:-1] + 4
xe = 0.5 * (tbins_8day[1:] + tbins_8day[:-1])

xticks = np.array([15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349])
xticklabels = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]


# =============================================================================
# Helper functions
# =============================================================================

def _safe_wilcoxon_p(x: np.ndarray) -> float:
    """
    Wilcoxon signed-rank test against zero.

    Returns NaN if the test cannot be performed.
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]

    if x.size < 5:
        return np.nan

    if np.allclose(x, 0):
        return np.nan

    try:
        stat, p = wilcoxon(x, alternative="two-sided", zero_method="zsplit")
        return p
    except ValueError:
        return np.nan


def _get_depth_axis(woa_bundle: dict, n_depth: int) -> np.ndarray:
    """
    Get the WOA depth axis matching the saved regional N2 section.
    """
    depth = np.asarray(woa_bundle["depth_TS_woa"], dtype=float)
    return depth[:n_depth]


def _get_regime_lines(region: str, depth_climt_dict: dict) -> tuple[np.ndarray, np.ndarray]:
    """
    Get vertical regime-boundary lines.
    For ED/Open, use MLD-Zeu crossings.
    For Low-N, use the manually specified t0 in REGION_LIMITS.
    """
    imld_holte = np.asarray(depth_climt_dict["mld_holte"], dtype=float)
    izeu_bg = np.asarray(depth_climt_dict["zeu_bg"], dtype=float)

    if region != "Low-N":
        crossings = find_crossings(tsteps, imld_holte, izeu_bg)
        tc = np.full(len(crossings), np.nan)
        yc = np.full(len(crossings), np.nan)

        for i, crossing in enumerate(crossings):
            tc[i], yc[i] = crossing

        return tc, yc

    t0 = REGION_LIMITS[region]["regime_t0"]
    return np.array([t0], dtype=float), np.array([np.nan], dtype=float)


def _choose_amplitude_group_variable(region_eddy: dict, ie: str) -> np.ndarray:
    """
    Use track-level amplitude if available.
    If all NaN, fall back to realization-level amplitude so the script can run.

    For exact paper reproduction, we should later fix Step 4 so amp_track_dict is populated.
    """
    amp_track = np.asarray(region_eddy["amp_track_dict"][ie], dtype=float)

    if np.isfinite(amp_track).sum() > 0:
        return amp_track

    print(
        f"WARNING: amp_track_dict[{ie}] is all NaN. "
        f"Using amp_realz_dict[{ie}] as temporary fallback."
    )
    return np.asarray(region_eddy["amp_realz_dict"][ie], dtype=float)


def _calculate_amplitude_group_chl(region_eddy: dict):
    """
    Calculate median/IQR eddy-induced CHL anomaly by amplitude quartile and 8-day bin.
    """
    chl_anom_sub = {"C": {}, "A": {}}
    plot_records = {"C": [], "A": []}

    for ie in ["C", "A"]:
        x = np.asarray(region_eddy["tstep_realz_dict"][ie], dtype=float)
        z = np.asarray(region_eddy["chl_anom_dict"][ie], dtype=float)
        amp_group_data = _choose_amplitude_group_variable(region_eddy, ie)

        df = pd.DataFrame(
            {
                "tstep_realz": x,
                "chl_anom": z,
                "amp_group": amp_group_data,
            }
        )

        df = df[np.isfinite(df["tstep_realz"]) & np.isfinite(df["chl_anom"]) & np.isfinite(df["amp_group"])]
        df["x_bin"] = pd.cut(df["tstep_realz"], bins=tbins_8day)
        df["x_bin_center"] = df["x_bin"].apply(lambda b: b.mid if pd.notnull(b) else np.nan)

        for isub in range(len(percentile_ls) - 1):
            p0 = percentile_ls[isub]
            p1 = percentile_ls[isub + 1]

            sub0 = np.nanpercentile(df["amp_group"], p0)
            sub1 = np.nanpercentile(df["amp_group"], p1)

            if isub == 0:
                mask_sub = (df["amp_group"] >= sub0) & (df["amp_group"] <= sub1)
            else:
                mask_sub = (df["amp_group"] > sub0) & (df["amp_group"] <= sub1)

            df_sub = df[mask_sub]

            val_ar = np.full(len(xe), np.nan)
            q1_ar = np.full(len(xe), np.nan)
            q3_ar = np.full(len(xe), np.nan)
            insufficient_ar = np.full(len(xe), False, dtype=bool)

            chl_anom_sub[ie][isub] = np.full(len(xe), np.nan)

            grouped = df_sub.groupby("x_bin_center", observed=False)

            for dx, group in grouped:
                if pd.isnull(dx):
                    continue

                ix = int(np.argmin(np.abs(xe - float(dx))))
                valid_vals = group["chl_anom"].dropna()

                if len(valid_vals) < 5:
                    insufficient_ar[ix] = True
                    continue

                p = _safe_wilcoxon_p(valid_vals.to_numpy())

                val_ar[ix] = valid_vals.median()
                q1_ar[ix] = valid_vals.quantile(0.25)
                q3_ar[ix] = valid_vals.quantile(0.75)
                insufficient_ar[ix] = (not np.isfinite(p)) or (p > 0.05)

                chl_anom_sub[ie][isub][ix] = val_ar[ix]

            plot_records[ie].append(
                {
                    "isub": isub,
                    "sub0": sub0,
                    "sub1": sub1,
                    "median": val_ar,
                    "q1": q1_ar,
                    "q3": q3_ar,
                    "insufficient": insufficient_ar,
                }
            )

    return chl_anom_sub, plot_records


def _add_regime_lines(ax, tc: np.ndarray):
    for itc in tc:
        if np.isfinite(itc):
            ax.axvline(x=itc, color="k", linestyle="--", linewidth=lw)


def _add_regime_arrows(ax, region: str):
    """
    Add the arrow annotations used in the submitted-style panels.
    These are approximate and can be manually tuned.
    """
    hdl = 1.8
    hdw = 0.9
    y_arw = 0.92

    if region == "ED":
        arrows = [
            ((0.3, y_arw), (0.02, y_arw), "->"),
            ((0.98, y_arw), (0.78, y_arw), "<-"),
            ((0.75, y_arw), (0.33, y_arw), "<->"),
        ]
    elif region == "Open":
        arrows = [
            ((0.33, y_arw), (0.02, y_arw), "->"),
            ((0.98, y_arw), (0.76, y_arw), "<-"),
            ((0.73, y_arw), (0.36, y_arw), "<->"),
        ]
    else:
        arrows = [
            ((0.45, y_arw), (0.02, y_arw), "<->"),
            ((0.98, y_arw), (0.50, y_arw), "<->"),
        ]

    for xy, xytext, style in arrows:
        ax.annotate(
            "",
            xy=xy,
            xycoords=ax.transAxes,
            xytext=xytext,
            textcoords=ax.transAxes,
            arrowprops=dict(
                arrowstyle=f"{style},head_length={hdl},head_width={hdw}",
                color="black",
                linewidth=4,
            ),
        )


# =============================================================================
# Main plotting function
# =============================================================================

def plot_region(region: str, bundle: dict, woa_bundle: dict) -> None:
    region_label = REGION_LABELS[region]
    fig_num = FIGURE_NUMBERS[region]

    depth_climt_dict = bundle["depth_climatologies"][region]
    region_eddy = bundle["eddy_climatologies"][region]

    iN2_woa = np.asarray(depth_climt_dict["N2_woa_2d"], dtype=float)
    imld_holte = np.asarray(depth_climt_dict["mld_holte"], dtype=float)
    izeu_bg = np.asarray(depth_climt_dict["zeu_bg"], dtype=float)
    incline_bg = np.asarray(depth_climt_dict["ncline_bg"], dtype=float)
    ichl_bg = np.asarray(depth_climt_dict["chl_bg"], dtype=float)
    n_matrix_8d = np.asarray(depth_climt_dict["n_matrix_8d"], dtype=float)

    depth_ar_woa = _get_depth_axis(woa_bundle, iN2_woa.shape[0])
    n_contour_ls = np.arange(0.1, 1.0, 0.1)

    tc, yc = _get_regime_lines(region, depth_climt_dict)
    chl_anom_sub, plot_records = _calculate_amplitude_group_chl(region_eddy)

    zmin = REGION_LIMITS[region]["zmin"]
    zmax = REGION_LIMITS[region]["zmax"]
    bg_max = REGION_LIMITS[region]["bg_max"]

    xlim = [5, 365]

    fig = plt.figure(figsize=(20, 22))
    gs = gridspec.GridSpec(
        7,
        1,
        height_ratios=[0.5, 5, 0.5, 5, 0.5, 5, 0.5],
        hspace=0.3,
    )

    # -------------------------------------------------------------------------
    # (c) Eddy-induced CHL anomaly
    # -------------------------------------------------------------------------
    ax = fig.add_subplot(gs[5])
    lines = []

    for ie in ["C", "A"]:
        for rec in plot_records[ie]:
            isub = rec["isub"]
            ic = c_dict[ie][isub]

            label = f"{ie}E, {rec['sub0']:.2f}-{rec['sub1']:.2f}"

            line, = ax.plot(
                xe,
                rec["median"],
                color=ic,
                linewidth=lw,
                label=label,
            )

            ax.fill_between(
                xe,
                rec["q1"],
                rec["q3"],
                color=ic,
                alpha=0.1,
            )

            insufficient = rec["insufficient"]
            if insufficient.any():
                ax.scatter(
                    xe[insufficient],
                    rec["median"][insufficient],
                    s=100,
                    color=ic,
                    marker="x",
                    linewidths=lw,
                    zorder=5,
                )

            lines.append(line)

    ax.set_xlim(xlim)
    ax.set_ylim([zmin, zmax])
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
    ax.set_ylabel(r"$\mathrm{CHL_e^{\mathbf{\prime}}\ (mg\ m^{-3})}$", fontsize=lbsz)
    ax.tick_params(labelsize=tksz)
    ax.axhline(y=0, color="black", linestyle="--", linewidth=lw * 0.5)
    ax.grid(True)

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)

    # Reorder legend: C quartile 1, A quartile 1, C quartile 2, A quartile 2...
    if len(lines) == 8:
        handles_ordered = [lines[0], lines[4], lines[1], lines[5], lines[2], lines[6], lines[3], lines[7]]
    else:
        handles_ordered = lines

    ax.legend(
        handles=handles_ordered,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.1),
        frameon=False,
        fontsize=lgsz,
        ncol=4,
        columnspacing=5.5,
    )

    ax.text(
        x_tt,
        y_tt,
        "(c) Eddy-induced chlorophyll-a anomaly",
        transform=ax.transAxes,
        fontsize=lbsz * 1.2,
        bbox=dict(facecolor="white", alpha=0.4, edgecolor="none", boxstyle="square,pad=0.1"),
        zorder=6,
    )

    _add_regime_lines(ax, tc)

    # -------------------------------------------------------------------------
    # (a) Stratification and key depths
    # -------------------------------------------------------------------------
    ax = fig.add_subplot(gs[1])

    levels = np.arange(0, 2.1e-4, 0.2e-4)
    cf = ax.contourf(
        tsteps,
        depth_ar_woa,
        iN2_woa,
        cmap="turbo",
        levels=levels,
        extend="both",
    )

    ax.plot(tsteps, imld_holte, color="black", linewidth=lw * 1.5, label="Mixed layer depth")
    ax.plot(tsteps, izeu_bg, color="white", linewidth=lw * 1.5, label="Euphotic depth")
    ax.plot(tsteps, incline_bg, color="magenta", linewidth=lw * 1.5, label="Nitrate isopleth")

    for iiN, _ in enumerate(n_contour_ls):
        ax.plot(tsteps, n_matrix_8d[iiN, :], color="grey", linewidth=lw * 0.5)

    ax.grid()
    ax.set_ylabel(r"Depth ($\mathrm{m}$)", fontsize=lbsz)
    ax.text(
        x_tt,
        y_tt,
        "(a) Stratification and key depths",
        transform=ax.transAxes,
        fontsize=lbsz * 1.2,
        bbox=dict(facecolor="white", alpha=0.4, edgecolor="none", boxstyle="square,pad=0.1"),
        zorder=6,
    )

    ax.set_xlim(xlim)
    ax.set_ylim(0, 200)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    ax.tick_params(labelsize=tksz)

    ax.legend(
        loc="lower right",
        ncol=1,
        bbox_to_anchor=(1.005, 0),
        frameon=True,
        fontsize=lgsz,
        columnspacing=1,
    )

    ax.invert_yaxis()

    cbar_ax = fig.add_axes([0.905, 0.646, 0.01, 0.199])
    cb = plt.colorbar(cf, cax=cbar_ax, orientation="vertical")
    cb.ax.tick_params(labelsize=tksz)

    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-2, 2))
    cb.ax.yaxis.set_major_formatter(formatter)
    cb.ax.yaxis.get_offset_text().set_visible(False)
    cb.ax.set_title(r"  $N^{2}$", pad=12, fontsize=tksz)
    cb.set_ticks([0, 1e-4, 2e-4])
    cb.set_ticklabels([r"0", r"$1 \times 10^{-4}$", r"$2 \times 10^{-4}$"])

    if region != "Low-N":
        ax.scatter(tc, yc, c="k", marker="^", s=400, zorder=6)

    _add_regime_lines(ax, tc)

    # -------------------------------------------------------------------------
    # (b) Background CHL plus anomaly
    # -------------------------------------------------------------------------
    ax = fig.add_subplot(gs[3])

    ax.plot(tsteps, ichl_bg, color="black", linestyle="-", linewidth=lw, label="background")

    for ie in ["C", "A"]:
        for isub in range(4):
            ic = c_dict[ie][isub]
            ichl_anom = chl_anom_sub[ie][isub]
            ichl_bg_anom = ichl_bg + ichl_anom
            ax.plot(tsteps, ichl_bg_anom, color=ic, linestyle="-", linewidth=lw * 0.5)

    ax.grid()
    ax.set_ylabel(r"CHL ($\mathrm{mg\ m^{-3}}$)", fontsize=lbsz)
    ax.text(
        x_tt,
        y_tt,
        "(b) Background chlorophyll-a with anomaly",
        transform=ax.transAxes,
        fontsize=lbsz * 1.2,
        bbox=dict(facecolor="white", alpha=0.4, edgecolor="none", boxstyle="square,pad=0.1"),
        zorder=6,
    )

    ax.set_xlim(xlim)
    ax.set_ylim(0, bg_max)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    ax.tick_params(labelsize=tksz)

    ax.legend(
        loc="lower right",
        ncol=1,
        frameon=True,
        fontsize=lgsz,
        columnspacing=1,
    )

    _add_regime_lines(ax, tc)
    _add_regime_arrows(ax, region)

    # -------------------------------------------------------------------------
    # Save
    # -------------------------------------------------------------------------
    filename_pdf = SAVE_DIR / f"figure_0{fig_num}.pdf"
    filename_png = SAVE_DIR / f"figure_0{fig_num}.png"

    fig.savefig(filename_pdf)
    fig.savefig(filename_png, dpi=300)
    plt.close(fig)

    print(f"Saved: {filename_pdf}")
    print(f"Saved: {filename_png}")


def main() -> None:
    bundle = load_pickle(CLIM_FILE)
    woa_bundle = load_pickle(WOA_FILE)

    for region in ["ED", "Open", "Low-N"]:
        plot_region(region, bundle, woa_bundle)


if __name__ == "__main__":
    main()