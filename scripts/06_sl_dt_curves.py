#!/usr/bin/env python3
"""
Figure 6 - Sea level and dynamic topography time series at the featured
NW Shelf well (asteras).

Single-panel figure with a dual y-axis:
    left  y: long-term sea level (Haq & Ogg 2024 hybrid long-term curve, m)
    right y: dynamic topography sampled at the well location, expressed
             as elevation anomaly relative to present (m), from the
             D10_gmcm9 model (Braz et al. 2021).

Time axis runs 150 Ma -> 0 Ma left to right (conventional chronological
direction).  150 Ma is the D10_gmcm9 oldest time slice; the asteras
record extends back to 146 Ma.

Output:
    figures/output/fig06_sl_dt_curves.png
    figures/output/fig06_sl_dt_curves.pdf
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

warnings.filterwarnings("ignore", "Dynamic topography model .* cannot")

# ----------------------------------------------------------------------------
# Tunables
# ----------------------------------------------------------------------------
FIG_SIZE = (10, 4.2)            # inches; wide single panel
DT_SAMPLE_STEP = 1.0            # Myr, sampling cadence for the DT curve
SL_COLOR = "#1f77b4"            # matplotlib blue
DT_COLOR = "#d62728"            # matplotlib red
SL_LINEWIDTH = 1.6
DT_LINEWIDTH = 1.6


# ----------------------------------------------------------------------------
# Data loaders
# ----------------------------------------------------------------------------
def load_sea_level_curve(model_name, age_max):
    """Read the bundled sea-level model file (2 columns: age_Ma, sl_m)
    and clip to [0, age_max].

    pybacktrack.bundle_data.BUNDLE_SEA_LEVEL_MODELS maps name -> path,
    so we use that lookup rather than hard-coding a path here.
    """
    sl_path = pybacktrack.BUNDLE_SEA_LEVEL_MODELS[model_name]
    data = np.loadtxt(sl_path)
    ages = data[:, 0]
    sl = data[:, 1]
    keep = ages <= age_max
    return ages[keep], sl[keep]


def load_dynamic_topography_at_well(model_name, lon, lat, age_max, step):
    """Sample bundled DT model at (lon, lat) for ages 0..age_max in Myr
    steps of `step`.  Returns anomalies relative to the present-day
    sample (matches the backstripping convention)."""
    dt = pybacktrack.DynamicTopography.create_from_bundled_model(
        model_name, lon, lat
    )
    ages = np.arange(0.0, age_max + 0.5 * step, step)
    samples = np.array([dt.sample(float(a)) for a in ages], dtype=float)
    # If samples ran off the model's age range we'd get NaN; report it.
    n_nan = int(np.isnan(samples).sum())
    if n_nan:
        print(f"  ! {n_nan} DT samples are NaN (off model coverage); "
              "they will appear as gaps in the curve.")
    anomaly = samples - samples[0]
    return ages, anomaly


def well_lat_lon(well_file):
    """Read the SiteLatitude / SiteLongitude header from a NWSHELF .txt."""
    lat = lon = None
    with open(well_file) as fh:
        for _ in range(20):
            line = fh.readline()
            if not line:
                break
            if "SiteLatitude" in line:
                try: lat = float(line.split("=")[-1])
                except ValueError: pass
            if "SiteLongitude" in line:
                try: lon = float(line.split("=")[-1])
                except ValueError: pass
    if lat is None or lon is None:
        raise ValueError(f"Could not read lat/lon from {well_file}")
    return lon, lat


# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------
def main():
    if not WELLS:
        sys.exit("No wells configured in config.WELLS")
    well = WELLS[0]                          # featured single-well example
    lon, lat = well_lat_lon(well["file"])
    print(f"=== {well['name']}  (lon={lon:.3f}, lat={lat:.3f}) ===")

    sl_ages, sl_m = load_sea_level_curve(SEA_LEVEL_MODEL, MAX_ANALYSIS_AGE)
    print(f"  SL ({SEA_LEVEL_MODEL}): {len(sl_ages)} samples "
          f"over 0..{MAX_ANALYSIS_AGE} Ma")
    dt_ages, dt_anom = load_dynamic_topography_at_well(
        DYNAMIC_TOPOGRAPHY_MODEL, lon, lat, MAX_ANALYSIS_AGE, DT_SAMPLE_STEP)
    print(f"  DT ({DYNAMIC_TOPOGRAPHY_MODEL}): {len(dt_ages)} samples "
          f"at ({lon:.3f}, {lat:.3f})")

    fig, ax_sl = plt.subplots(figsize=FIG_SIZE)
    ax_dt = ax_sl.twinx()

    # Font sizes (user 2026-06-02): bumped up to match Fig 4.
    AXIS_LABEL_SIZE = 16
    TICK_LABEL_SIZE = 14
    LEGEND_SIZE     = 13

    # X-axis: 150 Ma on the left, 0 on the right (calendar direction).
    ax_sl.set_xlim(MAX_ANALYSIS_AGE, 0)

    # Sea level (left axis, blue).
    ln_sl, = ax_sl.plot(sl_ages, sl_m,
                        color=SL_COLOR, lw=SL_LINEWIDTH,
                        label="Long-term eustatic sea level fluctuations")
    ax_sl.set_ylabel("Sea level (m)", color=SL_COLOR,
                     fontsize=AXIS_LABEL_SIZE)
    ax_sl.tick_params(axis="y", colors=SL_COLOR, labelsize=TICK_LABEL_SIZE)
    ax_sl.axhline(0, color=SL_COLOR, lw=0.5, ls=":", alpha=0.5)

    # Dynamic topography (right axis, red).
    ln_dt, = ax_dt.plot(dt_ages, dt_anom,
                        color=DT_COLOR, lw=DT_LINEWIDTH,
                        label="Dynamic topography anomaly")
    ax_dt.set_ylabel("Dynamic topography anomaly (m)",
                     color=DT_COLOR, fontsize=AXIS_LABEL_SIZE)
    ax_dt.tick_params(axis="y", colors=DT_COLOR, labelsize=TICK_LABEL_SIZE)
    ax_dt.axhline(0, color=DT_COLOR, lw=0.5, ls=":", alpha=0.5)

    ax_sl.set_xlabel("Age (Ma)", fontsize=AXIS_LABEL_SIZE)
    ax_sl.tick_params(axis="x", labelsize=TICK_LABEL_SIZE)
    # No title -- the manuscript figure caption carries the attribution.
    ax_sl.grid(True, alpha=0.3, linestyle=":")
    ax_sl.legend(handles=[ln_sl, ln_dt], loc="lower left",
                 framealpha=0.9, fontsize=LEGEND_SIZE)

    fig.tight_layout()
    base = os.path.join(OUTPUT_DIR, "fig06_sl_dt_curves")
    fig.savefig(base + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(base + ".pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {base}.png")
    print(f"wrote {base}.pdf")


if __name__ == "__main__":
    main()
