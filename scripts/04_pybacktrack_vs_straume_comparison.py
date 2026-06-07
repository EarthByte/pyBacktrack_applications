#!/usr/bin/env python3
"""
Figure 4 - Frame-invariant statistical comparison of the
pyBacktrack 1.5 paleobathymetry (Zahirovic et al. 2022 GDH1 release;
mantle frame by default) against the independent Straume et al. (2020)
Cenozoic paleobathymetry + paleotopography grids.

IMPORTANT -- the two products are reconstructed in DIFFERENT plate
reference frames, so the same (lon, lat) at a past time does NOT
correspond to the same physical seafloor in both models.  Per-pixel
difference fields would therefore be meaningless.  This script
deliberately reports only **frame-invariant** summary statistics
that depend on each model's *distribution* of bathymetric values at
time t, not on spatial registration.

Top panel: time-series of the **median** ocean depth at every common
time slice plus a shaded **median +/- 1 MAD** band (MAD = median
absolute deviation, the robust analogue of standard deviation).  The
user asked for the time-series panel to focus on these two robust
statistics only -- mean markers and p10/p90 envelopes were removed
(2026-05-31) because the p10/p90 bands of the two products are
visually almost identical and add clutter without information.

Bottom panel: depth-distribution histograms at
SPATIAL_HISTOGRAM_TIME_MA (default 40 Ma; one of the times shown in
Fig 2) so the reader can cross-reference the maps with the actual
shape of the bathymetric distribution at a single time slice.  Per-
product median is overlaid as a dashed vertical line.

Time coverage:
    pyBacktrack : 0..170 Ma (1 Myr)   [Zahirovic 2022 grids]
    Straume     : 1..65  Ma (1 Myr)
    -> common window = 1..65 Ma.

Outputs:
    figures/output/fig04_pybacktrack_vs_straume_stats.csv
    figures/output/fig04_pybacktrack_vs_straume.png
    figures/output/fig04_pybacktrack_vs_straume.pdf
"""
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (OUTPUT_DIR, PALEO_BATHY_FMT, STRAUME_GRID_DIR,
                    STRAUME_GRID_FMT)


# ----------------------------------------------------------------------------
# Tunables
# ----------------------------------------------------------------------------
# Time to overlay depth-distribution histograms at, in the bottom panel.
# 40 Ma is one of the times shown in Fig 2 of this paper, so the reader
# can cross-reference the maps with the distributions.
SPATIAL_HISTOGRAM_TIME_MA = 40


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def discover_common_times():
    """Sorted list of integer ages (Ma) where both grid sets have a
    file.
    """
    if not os.path.isdir(STRAUME_GRID_DIR):
        sys.exit(
            "Straume 2020 grid directory not found:\n"
            f"    {STRAUME_GRID_DIR}\n"
            "Download the Straume et al. (2020) paleobathymetry release\n"
            "from https://zenodo.org/records/4193576 and place the .nc\n"
            "files into one of the locations config.py looks at "
            "(see comment block above STRAUME_GRID_DIR in config.py),\n"
            "or point $PYBT_STRAUME_GRID_DIR at the directory you have.")
    straume_re = re.compile(r"paleobathy-topo_([0-9]+(?:\.[0-9]+)?)Ma_")
    s_times = set()
    for fn in os.listdir(STRAUME_GRID_DIR):
        m = straume_re.search(fn)
        if m:
            s_times.add(int(float(m.group(1))))
    a_times = set()
    for t in s_times:
        if os.path.exists(PALEO_BATHY_FMT.format(time=float(t))):
            a_times.add(t)
    return sorted(s_times & a_times)


def ocean_values_pybacktrack(t_ma):
    """1-D array of finite (= ocean) depth values from the pyBacktrack
    paleobathymetry grid at age t_ma.  Native 0.2 deg, no regridding.
    """
    fa = PALEO_BATHY_FMT.format(time=float(t_ma))
    if not os.path.exists(fa):
        return None
    arr = xr.open_dataset(fa).z.values
    return arr[np.isfinite(arr)]


def ocean_values_straume(t_ma):
    """1-D array of depth values from the Straume grid at age t_ma
    masked to z < 0 (= ocean).  Native 0.1 deg, no regridding.
    """
    fs = STRAUME_GRID_FMT.format(time=float(t_ma))
    if not os.path.exists(fs):
        return None
    arr = xr.open_dataset(fs).z.values
    return arr[arr < 0.0]


def ocean_area_fraction_pybacktrack(t_ma):
    """Fraction of the present-day-grid surface that is finite (ocean)
    in the pyBacktrack paleobathymetry at age t_ma.  Native 0.2 deg.
    """
    fa = PALEO_BATHY_FMT.format(time=float(t_ma))
    arr = xr.open_dataset(fa).z.values
    return float(np.isfinite(arr).mean())


def ocean_area_fraction_straume(t_ma):
    """Fraction of the Straume grid cells with z < 0 at age t_ma.
    Native 0.1 deg.
    """
    fs = STRAUME_GRID_FMT.format(time=float(t_ma))
    arr = xr.open_dataset(fs).z.values
    return float((arr < 0.0).mean())


def mad(values):
    """Median absolute deviation (no normalisation -- the raw MAD, in
    metres, matches what the y-axis is going to plot).
    """
    return float(np.median(np.abs(values - np.median(values))))


def per_time_stats(t_ma):
    A = ocean_values_pybacktrack(t_ma)
    S = ocean_values_straume(t_ma)
    if A is None or S is None or len(A) == 0 or len(S) == 0:
        return None
    return dict(
        time_Ma          = t_ma,
        n_A              = int(len(A)),
        n_S              = int(len(S)),
        median_A         = float(np.median(A)),
        mad_A            = mad(A),
        median_S         = float(np.median(S)),
        mad_S            = mad(S),
        ocean_frac_A     = ocean_area_fraction_pybacktrack(t_ma),
        ocean_frac_S     = ocean_area_fraction_straume(t_ma),
    )


# ----------------------------------------------------------------------------
# CSV writer
# ----------------------------------------------------------------------------
def write_csv(rows, path):
    fields = ("time_Ma", "n_A", "n_S",
              "median_A", "mad_A",
              "median_S", "mad_S",
              "ocean_frac_A", "ocean_frac_S")
    with open(path, "w") as fh:
        fh.write("# Frame-invariant per-time statistics.  PyBacktrack "
                 "(Zahirovic 2022 grids) and Straume 2020 grids are in "
                 "different plate reference frames so per-pixel "
                 "comparison is NOT meaningful; values below are "
                 "independent summaries of each product's depth "
                 "distribution at time t.  MAD = median absolute "
                 "deviation (robust analogue of standard deviation), "
                 "in metres.\n")
        fh.write(",".join(fields) + "\n")
        for r in rows:
            fh.write(",".join(
                str(r[f]) if isinstance(r[f], int) else f"{r[f]:.4f}"
                for f in fields
            ) + "\n")
    print(f"  wrote {path}")


# ----------------------------------------------------------------------------
# Figure: TOP -- time-series of median + median+/-1 MAD envelope per product.
#         BOTTOM -- depth-distribution histograms at hist_time_ma.
# ----------------------------------------------------------------------------
def plot_summary(rows, hist_time_ma, out_base):
    times = np.array([r["time_Ma"] for r in rows])
    med_A = np.array([r["median_A"] for r in rows])
    mad_A = np.array([r["mad_A"]    for r in rows])
    med_S = np.array([r["median_S"] for r in rows])
    mad_S = np.array([r["mad_S"]    for r in rows])

    A_color = "#1f77b4"     # blue  = PyBacktrack 1.5 (Z22 grids)
    S_color = "#d62728"     # red   = Straume 2020

    # Font sizes (user 2026-06-02): bumped up to be readable when the
    # figure is reproduced at single-column width in the docx.
    AXIS_LABEL_SIZE = 16
    TICK_LABEL_SIZE = 14
    LEGEND_SIZE     = 13
    TITLE_SIZE      = 16

    fig, (ax_ts, ax_hist) = plt.subplots(
        nrows=2, ncols=1, figsize=(8.5, 8.0),
        gridspec_kw=dict(height_ratios=[1.1, 1.0]),
    )

    # ---- Top: time-series of median + median +/- 1 MAD --------------------
    # Legend text uses the GENERIC "PyBacktrack 1.5" label (rather than
    # naming the specific plate-model release used as input) so the
    # figure does not need to be edited whenever the user swaps the
    # paleobathymetry source.
    ax_ts.fill_between(times, med_A - mad_A, med_A + mad_A,
                       color=A_color, alpha=0.20, linewidth=0,
                       label="PyBacktrack 1.5 (median ± MAD)")
    ax_ts.plot(times, med_A, color=A_color, lw=2.0,
               label="PyBacktrack 1.5 median")
    # Straume
    ax_ts.fill_between(times, med_S - mad_S, med_S + mad_S,
                       color=S_color, alpha=0.20, linewidth=0,
                       label="Straume 2020 (median ± MAD)")
    ax_ts.plot(times, med_S, color=S_color, lw=2.0,
               label="Straume 2020 median")

    # Tighter axis limits, calibrated to the data
    ymin = min((med_A - mad_A).min(), (med_S - mad_S).min())
    ymax = max((med_A + mad_A).max(), (med_S + mad_S).max())
    pad = 200.0
    ax_ts.set_xlim(times.max(), times.min())       # 65 -> 0 Ma
    ax_ts.set_ylim(np.floor((ymin - pad)/100)*100,
                   np.ceil((ymax + pad)/100)*100)
    ax_ts.set_xlabel("Age (Ma)", fontsize=AXIS_LABEL_SIZE)
    ax_ts.set_ylabel("Ocean depth (m)", fontsize=AXIS_LABEL_SIZE)
    ax_ts.tick_params(axis="both", labelsize=TICK_LABEL_SIZE)
    ax_ts.grid(True, alpha=0.3, linestyle=":")
    ax_ts.legend(loc="lower right", fontsize=LEGEND_SIZE,
                 framealpha=0.9, ncol=2)
    # Panel tag (a) in the upper-left corner.
    ax_ts.text(0.015, 0.97, "(a)", transform=ax_ts.transAxes,
               fontsize=TITLE_SIZE, fontweight="bold",
               va="top", ha="left",
               bbox=dict(facecolor="white", edgecolor="black",
                         linewidth=0.5, pad=3, alpha=0.85))

    # ---- Bottom: depth-distribution histograms at hist_time_ma -----------
    A = ocean_values_pybacktrack(hist_time_ma)
    S = ocean_values_straume(hist_time_ma)
    bins = np.arange(-6500, 1, 100)
    ax_hist.hist(A, bins=bins, density=True, alpha=0.45,
                 color=A_color, edgecolor="none",
                 label=f"PyBacktrack 1.5 ({hist_time_ma} Ma)")
    ax_hist.hist(S, bins=bins, density=True, alpha=0.45,
                 color=S_color, edgecolor="none",
                 label=f"Straume 2020 ({hist_time_ma} Ma)")
    ax_hist.axvline(np.median(A), color=A_color, ls="--", lw=1.5,
                    label=f"PyBacktrack 1.5 median ({np.median(A):+5.0f} m)")
    ax_hist.axvline(np.median(S), color=S_color, ls="--", lw=1.5,
                    label=f"Straume 2020 median ({np.median(S):+5.0f} m)")
    ax_hist.set_xlabel("Depth (m)", fontsize=AXIS_LABEL_SIZE)
    ax_hist.set_ylabel("Probability density (per 100 m bin)",
                       fontsize=AXIS_LABEL_SIZE)
    ax_hist.set_title(
        f"Seafloor depth distributions at {hist_time_ma} Ma ",
        fontsize=TITLE_SIZE,
    )
    ax_hist.tick_params(axis="both", labelsize=TICK_LABEL_SIZE)
    ax_hist.grid(True, alpha=0.3, linestyle=":")
    # Legend moved to upper-RIGHT so the (b) panel tag has room in the
    # upper-left corner.
    ax_hist.legend(loc="upper right", fontsize=LEGEND_SIZE, framealpha=0.9)
    # Tight x-limits so the plot field ends at 0 m on the right
    # (no white space past the y-axis 0 mark) and at -6500 m on the
    # left (crops the empty deep tail).
    ax_hist.set_xlim(-6500, 0)
    # Panel tag (b) in the upper-left corner.
    ax_hist.text(0.015, 0.97, "(b)", transform=ax_hist.transAxes,
                 fontsize=TITLE_SIZE, fontweight="bold",
                 va="top", ha="left",
                 bbox=dict(facecolor="white", edgecolor="black",
                           linewidth=0.5, pad=3, alpha=0.85))

    fig.tight_layout()
    fig.savefig(out_base + ".png", dpi=200, bbox_inches="tight")
    fig.savefig(out_base + ".pdf", bbox_inches="tight")
    print(f"  wrote {out_base}.png")
    print(f"  wrote {out_base}.pdf")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    times = discover_common_times()
    if not times:
        sys.exit("No common-time slices between pyBacktrack and Straume grids.")
    print(f"\nComputing frame-invariant per-time stats over "
          f"{len(times)} common slices ({times[0]}..{times[-1]} Ma) ...")
    rows = []
    for t in times:
        s = per_time_stats(t)
        if s is None:
            continue
        rows.append(s)
        if t % 5 == 0:
            print(f"  t={t:>3} Ma | "
                  f"A: med={s['median_A']:+5.0f} mad={s['mad_A']:4.0f} | "
                  f"S: med={s['median_S']:+5.0f} mad={s['mad_S']:4.0f} | "
                  f"ocean frac A={s['ocean_frac_A']:.3f} "
                  f"S={s['ocean_frac_S']:.3f}")

    out_csv = os.path.join(OUTPUT_DIR, "fig04_pybacktrack_vs_straume_stats.csv")
    write_csv(rows, out_csv)

    out_base = os.path.join(OUTPUT_DIR, "fig04_pybacktrack_vs_straume")
    plot_summary(rows, SPATIAL_HISTOGRAM_TIME_MA, out_base)
    print("Done.")


if __name__ == "__main__":
    main()
