#!/usr/bin/env python3
"""
Figure 11 - Box plots through time of the NW-Shelf-wide spatial
distribution of tectonic-subsidence rates.

Three stacked panels (one per configuration), sharing the same x-axis
of time (Ma), with one box per 10 Myr time slice in 0..150 Ma:

    A: rate without sea level / dynamic topography corrections
    C: rate with Haq2024 long-term SL + D10_gmcm9 DT (Braz et al. 2021)
    D: C - A (the "correction" the SL + DT model imposes on the
       inferred subsidence rate)

Time runs LEFT to RIGHT, so 150 Ma is on the left and 0 Ma on the right
of each panel.

Inputs (produced by 07a_backstrip_all_nwshelf.py):
    figures/output/nwshelf_subsidence/A_no_sl_no_dt_rate/rate_<t>.nc
    figures/output/nwshelf_subsidence/C_sl_and_dt_rate/rate_<t>.nc
    figures/output/nwshelf_subsidence/D_rate_difference/rate_<t>.nc

Outputs:
    figures/output/fig11_nwshelf_rate_boxplots.png
    figures/output/fig11_nwshelf_rate_boxplots.pdf
"""
import os
import sys

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import OUTPUT_DIR, MAX_ANALYSIS_AGE

# ----------------------------------------------------------------------------
# Tunables
# ----------------------------------------------------------------------------
OUT_BASE = os.path.join(OUTPUT_DIR, "nwshelf_subsidence")
RATE_INTERVAL = 5                        # Myr per bin (width of interval)
# Bin centres: 2.5, 7.5, 12.5, ... up to MAX_ANALYSIS_AGE - RATE_INTERVAL/2.
# Each bin spans [centre - RATE_INTERVAL/2, centre + RATE_INTERVAL/2)
# of age, so the 2.5 Ma box pools rate_0..rate_4, the 7.5 Ma box
# pools rate_5..rate_9, etc.  30 bins total over 0-150 Ma.
BIN_CENTRES = np.arange(RATE_INTERVAL / 2,
                        int(MAX_ANALYSIS_AGE),
                        RATE_INTERVAL)

# Inputs, in panel order.  Fields:
#   tag         -- subdirectory tag ("A", "C", "D")
#   rate_dir    -- where the rate_<t>.nc files live
#   inset       -- per-panel "(x) ..." annotation drawn in the inset
#                  corner of the panel (replaces the per-panel title)
#   colour      -- box facecolour
#   ylim        -- per-panel y-axis range (m/Myr)
#   inset_side  -- "left" or "right": which side of the panel hosts
#                  the inset box.
PANELS = [
    ("A", os.path.join(OUT_BASE, "A_no_sl_no_dt_rate"),
     "(a) Without sea level fluctuations and no dynamic topography "
     "anomalies",
     "#1f77b4",                           # blue
     (-10, 50),
     "left"),
    ("C", os.path.join(OUT_BASE, "C_sl_and_dt_rate"),
     "(b) With sea level fluctuations and no dynamic topography "
     "anomalies",
     "#2ca02c",                           # green
     (-10, 50),
     "left"),
    ("D", os.path.join(OUT_BASE, "D_rate_difference"),
     "(c) Difference (b) - (a)",
     "#d62728",                           # red
     (-20, 30),
     "left"),
]

# Box-plot cosmetics.  Box width capped to ~80% of RATE_INTERVAL so
# adjacent boxes don't visually collide.
BOX_WIDTH = 0.8 * RATE_INTERVAL          # data units (Myr) = 4.0 with 5 Myr bins
BOX_WHIS = (5, 95)                       # whiskers at 5/95 percentile
SHOW_OUTLIERS = False

FIG_SIZE = (12, 8)


# ----------------------------------------------------------------------------
# Data loaders
# ----------------------------------------------------------------------------
def gather_box_data(rate_dir):
    """For each BIN_CENTRE, pool every 1-Myr rate grid whose age falls
    inside the interval [centre - RATE_INTERVAL/2,
                         centre + RATE_INTERVAL/2)
    and concatenate all finite values into one sample.  Each box-plot
    box therefore summarises the spatial AND temporal variability of
    the rate field across the full 10 Myr interval.

    Returns:
        positions : list of bin centres (Myr) we actually have data for
        samples   : list of 1-D arrays (one concatenated sample per bin)
    """
    positions, samples = [], []
    half = RATE_INTERVAL / 2
    for centre in BIN_CENTRES:
        t_lo = int(round(centre - half))            # inclusive
        t_hi = int(round(centre + half))            # exclusive
        bucket = []
        for t in range(t_lo, t_hi):
            path = os.path.join(rate_dir, f"rate_{int(t)}.nc")
            if not os.path.exists(path):
                continue
            arr = xr.open_dataset(path).z.values
            v = arr[np.isfinite(arr)]
            if v.size:
                bucket.append(v)
        if not bucket:
            continue
        positions.append(float(centre))
        samples.append(np.concatenate(bucket))
    return positions, samples


# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------
# Font sizes (user 2026-06-02): bumped up to match Fig 4 / Figs 6-9.
AXIS_LABEL_SIZE = 16
TICK_LABEL_SIZE = 14
INSET_SIZE      = 14

# Anchor age (Ma) for the "Subsidence" / "Uplift" sign-convention
# labels on every panel.  Must be inside the visible x-range.
SIGN_X_MA       = 35

# Bin centres (Ma) to highlight on panel (c).  Used to mark the four
# 5-Myr intervals with the largest mean |D| (the same four intervals
# whose rate / DT maps are drawn as rows of Fig 10 and panels of
# Fig 12).
HIGHLIGHT_C_CENTRES_MA = [12.5, 27.5, 112.5, 132.5]
HIGHLIGHT_PEN = dict(edgecolor="black", linewidth=1.8)


def main():
    if not os.path.isdir(OUT_BASE):
        sys.exit(f"{OUT_BASE} not found -- run 07a_backstrip_all_nwshelf.py "
                 "first.")

    fig, axes = plt.subplots(
        nrows=len(PANELS), ncols=1, figsize=FIG_SIZE, sharex=True
    )
    if len(PANELS) == 1:
        axes = [axes]

    for ax, (tag, rate_dir, inset_text, colour, ylim, inset_side) in zip(
            axes, PANELS):
        if not os.path.isdir(rate_dir):
            ax.text(0.5, 0.5,
                    f"{rate_dir} missing -- run 07a first",
                    transform=ax.transAxes, ha="center",
                    fontsize=INSET_SIZE)
            continue

        positions, samples = gather_box_data(rate_dir)
        # Diagnostic so we can see what the gather step actually returned.
        if positions:
            sample_sizes = [len(s) for s in samples]
            print(f"  [{tag}] {len(positions)} boxes  "
                  f"positions {positions[0]:.1f} .. {positions[-1]:.1f}  "
                  f"sample sizes range {min(sample_sizes)}..{max(sample_sizes)}")
        if not positions:
            ax.text(0.5, 0.5,
                    f"No rate_*.nc grids in\n{rate_dir}",
                    transform=ax.transAxes, ha="center",
                    fontsize=INSET_SIZE)
            continue

        bp = ax.boxplot(
            samples,
            positions=positions,
            widths=BOX_WIDTH,
            whis=BOX_WHIS,
            showfliers=SHOW_OUTLIERS,
            patch_artist=True,
            medianprops=dict(color="black", linewidth=1.2),
        )
        for box in bp["boxes"]:
            box.set(facecolor=colour, alpha=0.55, edgecolor=colour,
                    linewidth=0.8)
        for whisker in bp["whiskers"]:
            whisker.set(color=colour, linewidth=0.8)
        for cap in bp["caps"]:
            cap.set(color=colour, linewidth=0.8)

        # Panel (c) only: thicken the box edges of the four highlighted
        # intervals (the same four 5-Myr windows that Fig 10 maps and
        # Fig 12 plots DT for).
        if tag == "D":
            highlight = set(HIGHLIGHT_C_CENTRES_MA)
            for pos, box in zip(positions, bp["boxes"]):
                if pos in highlight:
                    box.set(**HIGHLIGHT_PEN)

        ax.axhline(0, color="0.4", lw=0.7, ls=":")
        ax.set_ylabel("Rate (m/Myr)", fontsize=AXIS_LABEL_SIZE)
        ax.tick_params(axis="both", labelsize=TICK_LABEL_SIZE)
        ax.grid(True, axis="y", alpha=0.3, linestyle=":")
        # Pin the x-range PER PANEL to 0..MAX_ANALYSIS_AGE and lock
        # x-autoscale OFF so later calls (tight_layout, savefig)
        # don't reset it back to the boxplot's data range.
        ax.set_xlim(MAX_ANALYSIS_AGE, 0)
        ax.set_autoscalex_on(False)
        # Reset xticks per-axis.  ax.boxplot() above called
        # set_xticks(positions) AND set_xticklabels(str(positions)).
        # set_xticks alone doesn't overwrite the labels, so we pass
        # explicit string labels too.  Done on every axis (not just
        # axes[-1]) because sharex=True doesn't propagate the label
        # override reliably through tight_layout.
        major = np.arange(0, int(MAX_ANALYSIS_AGE) + 1, 10)
        minor = np.arange(0, int(MAX_ANALYSIS_AGE) + 1, RATE_INTERVAL)
        ax.set_xticks(major)
        ax.set_xticklabels([str(int(m)) for m in major])
        ax.set_xticks(minor, minor=True)
        # Per-panel y-window: top + middle (0..50 m/Myr) zoom on the
        # positive distribution; bottom (-40..+40 m/Myr) is the
        # symmetric correction field.  Whisker/box tails that extend
        # past the limit are clipped; medians and the bulk of the
        # inter-quartile box stay clearly readable.
        ax.set_ylim(*ylim)
        # Per-panel inset annotation (replaces the per-panel title).
        # Position chosen per panel by inset_side from PANELS so the
        # box overlays the quieter part of each distribution.
        if inset_side == "left":
            inset_x, inset_ha = 0.02, "left"
        else:
            inset_x, inset_ha = 0.98, "right"
        ax.text(
            inset_x, 0.95, inset_text,
            transform=ax.transAxes,
            ha=inset_ha, va="top",
            fontsize=INSET_SIZE,
            bbox=dict(boxstyle="round,pad=0.35",
                      facecolor="white", edgecolor="0.4", alpha=0.9),
        )

        # "Subsidence" / "Uplift" sign-convention labels.  Both labels
        # anchored at x = sign_x_ma (data coords, passed in per
        # figure so the labels stay inside whichever age window the
        # current figure covers); y is in axes coords (0.93 = near
        # top for "Subsidence", 0.04 = near bottom for "Uplift").
        # Uplift only labelled on panels whose y-range extends below
        # zero.
        SIGN_FONT  = "italic"
        SIGN_SIZE  = INSET_SIZE
        SIGN_COLOR = "0.25"
        trans_sign = mtransforms.blended_transform_factory(
            ax.transData, ax.transAxes)
        ax.text(SIGN_X_MA, 0.93, "Subsidence",
                transform=trans_sign, ha="left", va="top",
                fontsize=SIGN_SIZE, fontstyle=SIGN_FONT,
                color=SIGN_COLOR)
        if ylim[0] < 0:
            ax.text(SIGN_X_MA, 0.04, "Uplift",
                    transform=trans_sign, ha="left", va="bottom",
                    fontsize=SIGN_SIZE, fontstyle=SIGN_FONT,
                    color=SIGN_COLOR)

    # Shared x-axis: 150 Ma on the left, 0 Ma on the right.
    # MAJOR ticks (with labels) every 10 Myr -- keeps the axis
    # readable instead of cramming 31 labels on it.
    # MINOR ticks (no labels) every RATE_INTERVAL Myr (5) so the
    # bin boundaries are still marked, just unlabelled.
    axes[-1].set_xlim(MAX_ANALYSIS_AGE, 0)
    axes[-1].set_xticks(np.arange(0, int(MAX_ANALYSIS_AGE) + 1, 10))
    axes[-1].set_xticks(np.arange(0, int(MAX_ANALYSIS_AGE) + 1,
                                  RATE_INTERVAL), minor=True)
    axes[-1].set_xlabel("Age (Ma)", fontsize=AXIS_LABEL_SIZE)
    # No fig.suptitle -- the manuscript caption carries the description.

    fig.tight_layout()
    # Diagnostic: confirm the actual xlim of each axis at save time.
    for i, ax in enumerate(axes):
        print(f"  axis[{i}] final xlim = {ax.get_xlim()}")
    base = os.path.join(OUTPUT_DIR, "fig11_nwshelf_rate_boxplots")
    # NOTE: dropped bbox_inches="tight" -- it was cropping the saved
    # figure to the rendered boxplot's data range (which, due to a
    # matplotlib quirk involving boxes at sparse old-age positions,
    # ended at ~77.5 Ma) instead of the axis xlim (0..150 Ma).
    fig.savefig(base + ".png", dpi=300)
    fig.savefig(base + ".pdf")
    plt.close(fig)
    print(f"\nwrote {base}.png")
    print(f"wrote {base}.pdf")


if __name__ == "__main__":
    main()
