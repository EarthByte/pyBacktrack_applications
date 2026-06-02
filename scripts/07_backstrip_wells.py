#!/usr/bin/env python3
"""
Figure 7 - Tectonic-subsidence at the featured NW Shelf well
(asteras), three backstripping configurations overlaid on a single
panel with their min/max paleo-water-depth uncertainty envelopes:

    (a)  no sea level, no dynamic topography           -- baseline
    (b)  Haq2024 long-term hybrid sea level only       -- isolates SL effect
    (c)  Haq2024 long-term sea level + D10_gmcm9       -- full correction
         dynamic topography (Braz et al. 2021)

Mirrors the bundled `pybacktrack/notebooks/backstrip.ipynb` workflow
(cells 12, 14, 16/18, 20) but draws all three runs on one shared
axis instead of three separate panels.

Time axis runs 150 Ma -> 0 Ma left to right.

Output:
    figures/output/fig07_backstripping_3configs.png
    figures/output/fig07_backstripping_3configs.pdf
"""
import os
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pybacktrack

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (WELLS, OUTPUT_DIR, DYNAMIC_TOPOGRAPHY_MODEL,
                    SEA_LEVEL_MODEL, MAX_ANALYSIS_AGE)

warnings.filterwarnings(
    "ignore",
    "Well thickness .* is larger than the total sediment thickness")
warnings.filterwarnings("ignore", "Dynamic topography model .* cannot")

# ----------------------------------------------------------------------------
# Tunables
# ----------------------------------------------------------------------------
FIG_SIZE = (10, 5.2)            # inches
# Y-axis depth maximum is auto-scaled from the deepest envelope across
# the three configurations and snapped up to a multiple of
# DEPTH_AXIS_STEP, with DEPTH_AXIS_PAD m of headroom past the data.
# The top of the (inverted) axis extends at least to -200 m so the
# SL+DT trace is not clipped where it briefly excurses above zero.
DEPTH_AXIS_PAD = 200            # m of headroom past the data envelope
DEPTH_AXIS_STEP = 500           # m, snap-to grid for the rounded max

# Each config: (short_id, label, kwargs to backstrip_well, subtract DT?, colour).
CONFIGS = [
    ("no_sl_no_dt",
     "(a) no SL, no DT",
     dict(),
     False,
     "#7f7f7f"),                 # neutral grey
    ("sl_only",
     "(b) Haq2024 SL only",
     dict(sea_level_model=SEA_LEVEL_MODEL),
     False,
     "#1f77b4"),                 # blue
    ("sl_and_dt",
     f"(c) Haq2024 SL + {DYNAMIC_TOPOGRAPHY_MODEL} DT",
     dict(sea_level_model=SEA_LEVEL_MODEL),
     True,
     "#d62728"),                 # red
]


# ----------------------------------------------------------------------------
# Backstripping
# ----------------------------------------------------------------------------
def backstrip_one(well_file, extra_kwargs, subtract_dyn_topo):
    """Return ages (Ma), mean tectonic subsidence (m), and min/max envelope.

    `subtract_dyn_topo`: if True, sample the configured DT model at the
    well location and subtract the elevation anomaly (relative to
    present) from the tectonic-subsidence trace -- the convention used
    in the bundled backstrip.ipynb.
    """
    kwargs = dict(lithology_filenames=pybacktrack.BUNDLE_LITHOLOGY_FILENAMES)
    kwargs.update(extra_kwargs)
    well, decomp = pybacktrack.backstrip_well(well_file, **kwargs)

    ages = np.array([d.get_age() for d in decomp])
    tect = np.array([d.get_tectonic_subsidence() for d in decomp])
    mn_mx = np.array([d.get_min_max_tectonic_subsidence() for d in decomp])
    tmin, tmax = mn_mx[:, 0], mn_mx[:, 1]

    if subtract_dyn_topo:
        dt = pybacktrack.DynamicTopography.create_from_bundled_model(
            DYNAMIC_TOPOGRAPHY_MODEL, well.longitude, well.latitude
        )
        samples = np.array([dt.sample(float(a)) for a in ages], dtype=float)
        samples = np.where(np.isnan(samples), 0.0, samples)
        dyn_anom = samples - samples[0]
        tect = tect - dyn_anom
        tmin = tmin - dyn_anom
        tmax = tmax - dyn_anom

    return ages, tect, tmin, tmax, well


# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------
def main():
    if not WELLS:
        sys.exit("No wells configured in config.WELLS")
    well = WELLS[0]
    print(f"=== {well['name']}  ({well['file']}) ===")

    runs = []
    for cfg_short, cfg_label, kwargs, sub_dt, colour in CONFIGS:
        print(f"  backstripping {cfg_short} ...")
        ages, tect, tmin, tmax, _ = backstrip_one(
            well["file"], kwargs, sub_dt)
        runs.append((cfg_short, cfg_label, colour, ages, tect, tmin, tmax))

    fig, ax = plt.subplots(figsize=FIG_SIZE)
    # X-axis: 150 Ma on the left, 0 Ma on the right (calendar direction).
    ax.set_xlim(MAX_ANALYSIS_AGE, 0)
    # Y-axis: subsidence is positive *downward*; invert so deeper -> lower
    # on the page (matches the convention of every other subsidence plot).
    # Auto-scale to the deepest envelope across all three configurations
    # so the data are never clipped.
    deepest = max(
        float(tmax.max()) for _, _, _, _, _, _, tmax in runs
    )
    depth_axis_max = (
        float(np.ceil((deepest + DEPTH_AXIS_PAD) / DEPTH_AXIS_STEP)
              * DEPTH_AXIS_STEP)
    )
    # Upper end (visually at the TOP of the inverted axis): show down to
    # at least -200 m (i.e. 200 m above present sea level) so the
    # SL+DT-corrected (config c) curve has headroom when it briefly
    # excurses above zero in early history.  Extend further if any of
    # the three envelopes actually crosses below that floor.
    shallowest = min(
        float(tmin.min()) for _, _, _, _, _, tmin, _ in runs
    )
    depth_axis_min = min(-200.0, shallowest - DEPTH_AXIS_PAD)
    depth_axis_min = (
        float(np.floor(depth_axis_min / DEPTH_AXIS_STEP) * DEPTH_AXIS_STEP)
    )
    print(f"  envelope: {shallowest:.0f} m (top) -> {deepest:.0f} m (bottom)")
    print(f"  y-axis  : {depth_axis_min:.0f} m (top) -> "
          f"{depth_axis_max:.0f} m (bottom)")
    ax.set_ylim(depth_axis_max, depth_axis_min)

    for cfg_short, cfg_label, colour, ages, tect, tmin, tmax in runs:
        # Min/max envelope.
        ax.fill_between(ages, tmin, tmax,
                        color=colour, alpha=0.18, linewidth=0)
        # Mean curve.
        ax.plot(ages, tect, color=colour, lw=2.0, label=cfg_label)
        # Sample points so the user can see where the stratigraphy is
        # actually defined (asteras has only 20 strat horizons).
        ax.scatter(ages, tect, s=18, color=colour, edgecolor="white",
                   linewidth=0.5, zorder=3)

    # Font sizes (user 2026-06-02): bumped up to match Fig 4.
    AXIS_LABEL_SIZE = 16
    TICK_LABEL_SIZE = 14
    LEGEND_SIZE     = 13
    INSET_SIZE      = 14

    ax.set_xlabel("Age (Ma)", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel("Tectonic subsidence (m)", fontsize=AXIS_LABEL_SIZE)
    ax.tick_params(axis="both", labelsize=TICK_LABEL_SIZE)
    ax.grid(True, alpha=0.3, linestyle=":")
    ax.legend(loc="lower left", framealpha=0.9, fontsize=LEGEND_SIZE)
    # Inset annotation in the upper-right corner (replaces the figure
    # title that was previously drawn above the axes).
    title_text = (
        f"Tectonic subsidence of {well['name']}\n"
        f"Three backstripping configurations\n"
        f"Envelopes = Paleo-water depth uncertainty"
    )
    ax.text(
        0.98, 0.98, title_text,
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=INSET_SIZE,
        bbox=dict(boxstyle="round,pad=0.4",
                  facecolor="white", edgecolor="0.4", alpha=0.9),
    )

    fig.tight_layout()
    base = os.path.join(OUTPUT_DIR, "fig07_backstripping_3configs")
    fig.savefig(base + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(base + ".pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {base}.png")
    print(f"wrote {base}.pdf")


if __name__ == "__main__":
    main()
