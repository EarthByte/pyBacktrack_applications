#!/usr/bin/env python3
"""
Figure 9 - NW Shelf tectonic-subsidence rate time-series.

Reads the gridded rate fields produced by
`07a_backstrip_all_nwshelf.py` and computes the spatial mean and
standard deviation across the NW Shelf at every 1 Myr time slice.
Three curves are plotted on a single panel:

    A : no sea level, no dynamic topography
    C : Haq2024 long-term hybrid SL + D10_gmcm9 DT (Braz et al. 2021)
    D : C - A   (the correction the SL + DT model imposes on the
                 inferred rate)

Shaded bands show +/- 1 sigma; a horizontal zero line is drawn for
reference (positive = subsiding, negative = uplifting).  Time axis
runs 150 Ma -> 0 Ma left to right.

Output:
    figures/output/fig09_nwshelf_rate_timeseries.png
    figures/output/fig09_nwshelf_rate_timeseries.pdf
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import OUTPUT_DIR, DYNAMIC_TOPOGRAPHY_MODEL, MAX_ANALYSIS_AGE

# ----------------------------------------------------------------------------
# Tunables
# ----------------------------------------------------------------------------
OUT_BASE = os.path.join(OUTPUT_DIR, "nwshelf_subsidence")
TIME_STEP = 1                                # Myr

# Display labels relabel C -> "B" and D -> "C" so the inset legend
# reads as A / B / C instead of the internal A / C / D (which still
# apply to directory names and the rest of the pipeline).
PANELS = [
    ("A", os.path.join(OUT_BASE, "A_no_sl_no_dt_rate"),
     "A: no SL, no DT", "#1f77b4"),           # blue
    ("C", os.path.join(OUT_BASE, "C_sl_and_dt_rate"),
     f"B: Haq2024 SL + {DYNAMIC_TOPOGRAPHY_MODEL} DT", "#2ca02c"),  # green
    ("D", os.path.join(OUT_BASE, "D_rate_difference"),
     "C = B - A", "#d62728"),                 # red
]

FIG_SIZE = (10, 5.0)


# ----------------------------------------------------------------------------
# Aggregation
# ----------------------------------------------------------------------------
def aggregate_from_grids(rate_dir):
    """Walk rate_dir, returning (times, means, stds) at every 1 Myr slice
    for which a rate_<t>.nc grid exists in 0..MAX_ANALYSIS_AGE."""
    times, means, stds = [], [], []
    for t in range(0, int(MAX_ANALYSIS_AGE) + 1, int(TIME_STEP)):
        path = os.path.join(rate_dir, f"rate_{int(t)}.nc")
        if not os.path.exists(path):
            continue
        arr = xr.open_dataset(path).z.values
        v = arr[np.isfinite(arr)]
        if v.size == 0:
            continue
        times.append(t)
        means.append(float(v.mean()))
        stds.append(float(v.std()))
    return (np.array(times), np.array(means), np.array(stds))


# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------
def main():
    if not os.path.isdir(OUT_BASE):
        sys.exit(f"{OUT_BASE} not found -- run 05_backstrip_all_nwshelf.py first.")

    series = []
    for tag, rate_dir, label, colour in PANELS:
        if not os.path.isdir(rate_dir):
            print(f"  ! {rate_dir} missing - run 05_backstrip_all_nwshelf.py")
            continue
        t, m, s = aggregate_from_grids(rate_dir)
        if len(t) == 0:
            print(f"  ! no rate_*.nc grids found in {rate_dir}")
            continue
        print(f"  {label}: {len(t)} time slices "
              f"(rate range [{m.min():.1f}, {m.max():.1f}] m/Myr, "
              f"max sigma {s.max():.1f} m/Myr)")
        series.append((tag, t, m, s, label, colour))

    if not series:
        sys.exit("No data series to plot.")

    fig, ax = plt.subplots(figsize=FIG_SIZE)
    # Time axis: 150 Ma on the left, 0 Ma on the right
    ax.set_xlim(MAX_ANALYSIS_AGE, 0)
    # Y-axis fixed at -50 to +75 m/Myr.  This clips a portion of the
    # +/- 1 sigma envelope at recent times where rates briefly spike
    # outside this band, in exchange for a much more legible plot for
    # the bulk of the record (the 150-10 Ma curves were too compressed
    # against the very-recent extremes with auto-scaling).
    ax.set_ylim(-50, 75)

    # Envelopes first so the mean lines sit on top.
    for tag, t, m, s, label, colour in series:
        ax.fill_between(t, m - s, m + s,
                        color=colour, alpha=0.18, linewidth=0)
    for tag, t, m, s, label, colour in series:
        ax.plot(t, m, color=colour, lw=2.0, label=label)

    # Font sizes (user 2026-06-02): bumped up to match Fig 4.
    AXIS_LABEL_SIZE = 16
    TICK_LABEL_SIZE = 14
    LEGEND_SIZE     = 13

    ax.axhline(0, color="0.4", lw=0.7, ls=":")
    ax.set_xlabel("Age (Ma)", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel("Tectonic-subsidence rate (m/Myr)",
                  fontsize=AXIS_LABEL_SIZE)
    ax.tick_params(axis="both", labelsize=TICK_LABEL_SIZE)
    # Title intentionally removed -- the manuscript figure caption
    # carries the full description.
    ax.grid(True, alpha=0.3, linestyle=":")
    ax.legend(loc="upper left", framealpha=0.9, fontsize=LEGEND_SIZE)

    fig.tight_layout()
    base = os.path.join(OUTPUT_DIR, "fig09_nwshelf_rate_timeseries")
    fig.savefig(base + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(base + ".pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {base}.png")
    print(f"wrote {base}.pdf")


if __name__ == "__main__":
    main()
