"""
Make paper-style Figure 8: lifecycle-climatology diagrams of eddy-induced
chlorophyll-a anomaly by climatological day and eddy age.

Inputs
------
- step04_climatology_outputs.pkl

Outputs
-------
- figure_08.pdf
- figure_08.png
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import gridspec
from scipy.stats import wilcoxon

from eac_eddy_chl.io import load_pickle


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


# =============================================================================
# Settings
# =============================================================================

REGION_ORDER = ["ED", "Open", "Low-N"]
POLARITY_ORDER = ["C", "A"]

REGION_TITLES = {
    "ED": "Eddy District (0.1 $\\mathbf{mg\\ m^{-3}}$)",
    "Open": "Tasman Sea (0.03 $\\mathbf{mg\\ m^{-3}}$)",
    "Low-N": "Northern Low (0.01 $\\mathbf{mg\\ m^{-3}}$)",
}

REGION_LABELS = {
    "ED": "ED",
    "Open": "TS",
    "Low-N": "NL",
}

# Match your paper limits.
LIM_DICT = {
    "Low-N": {"C": {"2d": 0.01}, "A": {"2d": 0.01}},
    "ED": {"C": {"2d": 0.10}, "A": {"2d": 0.10}},
    "Open": {"C": {"2d": 0.03}, "A": {"2d": 0.03}},
}

# 8-day bins: x-axis = climatological day; y-axis = eddy age.
TBINS_8DAY = np.append(np.arange(1, 362, 8), 367)
TSTEPS = TBINS_8DAY[:-1] + 4
XE = 0.5 * (TBINS_8DAY[1:] + TBINS_8DAY[:-1])

YBIN = np.arange(0, 8 * 47, 8)
YE = 0.5 * (YBIN[1:] + YBIN[:-1])

XLIM = [5, 365]
YLIM = [0, 46 * 8]

LETTERS = "abcdefghijklmnopqrstuvwxyz"
LW = 3

# Month-like ticks, using every second month label as in your paper code.
MONTH_TICK_DOY = np.array([15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349])
MONTH_LABELS = np.array(["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"])
XTICKS = MONTH_TICK_DOY[1::2]
XTICKLABELS = MONTH_LABELS[1::2]


# =============================================================================
# Helpers
# =============================================================================

def safe_wilcoxon_p(values: pd.Series | np.ndarray) -> float:
    """Return Wilcoxon signed-rank p-value against zero, or NaN if unavailable."""
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]

    if vals.size < 5:
        return np.nan

    if np.allclose(vals, 0):
        return np.nan

    try:
        _, p = wilcoxon(vals, alternative="two-sided", zero_method="zsplit")
        return float(p)
    except ValueError:
        return np.nan


def find_crossings(x: np.ndarray, y1: np.ndarray, y2: np.ndarray) -> list[tuple[float, float]]:
    """
    Find approximate crossings between two 1D curves by linear interpolation.
    Used here for MLD-Zeu regime boundaries.
    """
    x = np.asarray(x, dtype=float)
    y1 = np.asarray(y1, dtype=float)
    y2 = np.asarray(y2, dtype=float)
    diff = y1 - y2

    crossings: list[tuple[float, float]] = []

    for i in range(len(x) - 1):
        if not np.all(np.isfinite([diff[i], diff[i + 1], x[i], x[i + 1], y1[i], y1[i + 1], y2[i], y2[i + 1]])):
            continue

        if diff[i] == 0:
            crossings.append((x[i], y1[i]))
            continue

        if diff[i] * diff[i + 1] < 0:
            frac = abs(diff[i]) / (abs(diff[i]) + abs(diff[i + 1]))
            xc = x[i] + frac * (x[i + 1] - x[i])
            yc = y1[i] + frac * (y1[i + 1] - y1[i])
            crossings.append((xc, yc))

    return crossings


def get_regime_boundaries(region: str, depth_dict: dict) -> list[float]:
    """Get vertical dashed-line positions for Figure 8."""
    if region == "Low-N":
        # Manual fallback matching the older plotting workflow. Adjust if needed.
        return [181.0]

    crossings = find_crossings(
        TSTEPS,
        np.asarray(depth_dict["mld_holte"], dtype=float),
        np.asarray(depth_dict["zeu_bg"], dtype=float),
    )

    return [float(c[0]) for c in crossings]


def build_lifecycle_matrix(region_eddy: dict, polarity: str) -> dict[str, np.ndarray]:
    """
    Bin eddy-induced CHL-a anomaly by eddy age and climatological day.
    """
    age_region = np.asarray(region_eddy["age_realz_dict"][polarity], dtype=float)
    tstep_region = np.asarray(region_eddy["tstep_realz_dict"][polarity], dtype=float)
    chl_region = np.asarray(region_eddy["chl_anom_dict"][polarity], dtype=float)

    df = pd.DataFrame(
        {
            "age_realz": age_region,
            "tstep_realz": tstep_region,
            "chl_anom": chl_region,
        }
    )
    df = df[np.isfinite(df["age_realz"]) & np.isfinite(df["tstep_realz"]) & np.isfinite(df["chl_anom"])]

    df["y_bin"] = pd.cut(df["age_realz"], bins=YBIN)
    df["x_bin"] = pd.cut(df["tstep_realz"], bins=TBINS_8DAY)
    df["y_bin_center"] = df["y_bin"].apply(lambda b: b.mid if pd.notnull(b) else np.nan)
    df["x_bin_center"] = df["x_bin"].apply(lambda b: b.mid if pd.notnull(b) else np.nan)

    med_matrix = np.full((len(YBIN) - 1, len(TBINS_8DAY) - 1), np.nan)
    count_matrix = np.zeros_like(med_matrix)
    insufficient_matrix = np.full_like(med_matrix, False, dtype=bool)

    grouped = df.groupby(["y_bin_center", "x_bin_center"], observed=False)

    for (dy, dx), group in grouped:
        if pd.isnull(dy) or pd.isnull(dx):
            continue

        iy = int(np.argmin(np.abs(YE - float(dy))))
        ix = int(np.argmin(np.abs(XE - float(dx))))

        valid_vals = group["chl_anom"].dropna()
        count_matrix[iy, ix] = len(valid_vals)

        if len(valid_vals) < 5:
            insufficient_matrix[iy, ix] = True
            continue

        med_matrix[iy, ix] = valid_vals.median()
        p = safe_wilcoxon_p(valid_vals.to_numpy())
        insufficient_matrix[iy, ix] = (not np.isfinite(p)) or (p > 0.05)

    return {
        "med_matrix": med_matrix,
        "count_matrix": count_matrix,
        "insufficient_matrix": insufficient_matrix,
    }


def build_all_lifecycle_matrices(bundle: dict) -> dict:
    """Build lifecycle-climatology matrices for all regions and polarities."""
    out = {}

    for region in REGION_ORDER:
        out[region] = {}
        region_eddy = bundle["eddy_climatologies"][region]

        for polarity in POLARITY_ORDER:
            out[region][polarity] = build_lifecycle_matrix(region_eddy, polarity)

    return out


def add_insufficient_markers(ax, insufficient_matrix: np.ndarray) -> None:
    """Add x markers to bins that are insufficient or not significantly different from zero."""
    for iy in range(insufficient_matrix.shape[0]):
        for ix in range(insufficient_matrix.shape[1]):
            if insufficient_matrix[iy, ix]:
                ax.text(
                    XE[ix],
                    YE[iy],
                    "x",
                    color="gray",
                    ha="center",
                    va="center",
                    fontsize=10,
                    fontweight="bold",
                )


def plot_figure08(bundle: dict, lifecycle_dict: dict) -> None:
    """Plot Figure 8."""
    fig = plt.figure(figsize=(24, 12))
    gs = gridspec.GridSpec(
        2,
        9,
        figure=fig,
        width_ratios=[1, 0.05, 0.2, 1, 0.05, 0.2, 1, 0.05, 0.2],
        wspace=0.1,
        hspace=0.1,
    )

    for i_region, region in enumerate(REGION_ORDER):
        tcs = get_regime_boundaries(region, bundle["depth_climatologies"][region])

        for i_pol, polarity in enumerate(POLARITY_ORDER):
            zlim = LIM_DICT[region][polarity]["2d"]
            zmin2d = -zlim
            zmax2d = zlim

            ax = fig.add_subplot(gs[i_pol, i_region * 3])

            z = lifecycle_dict[region][polarity]["med_matrix"]
            z_insig = lifecycle_dict[region][polarity]["insufficient_matrix"]

            pc = ax.pcolormesh(
                XE,
                YE,
                z,
                cmap="RdBu_r",
                vmin=zmin2d,
                vmax=zmax2d,
                shading="auto",
            )

            ax.set_xlim(XLIM)
            ax.set_ylim(YLIM)
            ax.set_xticks(XTICKS)
            ax.set_xticklabels(XTICKLABELS)
            ax.tick_params(labelsize=18)

            for tc in tcs:
                if np.isfinite(tc):
                    ax.axvline(x=tc, c="k", linestyle="--", linewidth=LW)

            add_insufficient_markers(ax, z_insig)

            if i_pol == 0:
                cax = fig.add_subplot(gs[:, i_region * 3 + 1])
                fig.colorbar(pc, cax=cax, ticks=np.linspace(zmin2d, zmax2d, 9))
                ax.set_title(REGION_TITLES[region], fontsize=22, fontweight="bold")

            if i_region == 0:
                ax.set_ylabel("Eddy age (d)", fontsize=20)
            else:
                ax.set_yticklabels([])

            if i_pol == 1:
                ax.set_xlabel("Month", fontsize=20)
            else:
                ax.set_xticklabels([])

            ax.text(
                0.05,
                0.9,
                f"({LETTERS[i_region * 2 + i_pol]})",
                transform=ax.transAxes,
                fontsize=20,
                fontweight="bold",
                ha="left",
                va="center",
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", boxstyle="square,pad=0.1"),
            )

            if (region == "Low-N") and (polarity == "A"):
                ax.text(
                    0.5,
                    0.9,
                    "eddy aging \nin time",
                    fontsize=20,
                    fontweight="bold",
                    ha="left",
                    va="center",
                    transform=ax.transAxes,
                    bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", boxstyle="square,pad=0.1"),
                )
                x0, y0 = 260, 260
                dx, dy = 80, 80
                ax.annotate(
                    "",
                    xy=(x0 + dx, y0 + dy),
                    xytext=(x0, y0),
                    arrowprops=dict(facecolor="black", edgecolor="black", shrink=0, width=4, headwidth=15),
                    annotation_clip=False,
                )

    fig.text(0.09, 0.7, "CE", ha="center", va="center", fontsize=32, fontweight="bold")
    fig.text(0.09, 0.3, "AE", ha="center", va="center", fontsize=32, fontweight="bold")

    out_pdf = SAVE_DIR / "figure_08.pdf"
    out_png = SAVE_DIR / "figure_08.png"

    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_png}")


def main() -> None:
    bundle = load_pickle(CLIM_FILE)
    lifecycle_dict = build_all_lifecycle_matrices(bundle)
    plot_figure08(bundle, lifecycle_dict)


if __name__ == "__main__":
    main()
