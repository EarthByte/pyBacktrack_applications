#!/usr/bin/env python3
"""
Figure 8 - Geohistory analysis (Wheeler diagram) for the featured
NW Shelf well (asteras).

Reproduces the bundled `geohistory_analysis.ipynb` Wheeler diagram
in the standard geological positive-up sign convention:

    - stacked decompacted stratigraphic units, stepping downward
      from the time-dependent sea floor
    - sea-floor curve (sea level minus water depth)
    - sea-level curve (Haq2024 long-term hybrid; positive-up,
      highstand above zero, lowstand below)
    - basement / tectonic-subsidence curve (large negative)
    - dynamic-topography elevation anomaly (D10_gmcm9, Braz et al.
      2021) relative to present (positive = uplift)

Time axis runs 150 Ma -> 0 Ma left to right.  Y-axis is depth below
present-day sea level in metres: the zero datum marks present-day
sea level, positive values plot above it, negative below.

Output:
    figures/output/fig08_geohistory_<wellname>.png
    figures/output/fig08_geohistory_<wellname>.pdf
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

# Stratigraphic-unit fill colours, cycled if there are more units than colours.
UNIT_COLOURS = [
    "mistyrose", "bisque", "lemonchiffon", "lightgreen", "aquamarine",
    "powderblue", "lavender", "thistle", "khaki", "palegreen",
    "lightsteelblue", "wheat", "peachpuff", "mintcream", "lavenderblush",
    "honeydew", "linen", "antiquewhite",
]


# ----------------------------------------------------------------------------
# Decompaction history (port of bundled notebook's helper)
# ----------------------------------------------------------------------------
def get_decompacted_history(well, decomp, dt_model):
    """Return ages, sea levels, surface (sea-floor) elevations,
    tectonic subsidence, dynamic-topography anomalies, and per-unit
    cumulative elevations.  All in m relative to present-day sea level.

    Sign convention: **positive UP** (the standard geological geohistory
    convention, also what every sea-level curve in the literature uses).

      * Sea level: positive when above present sea level (highstand),
        negative when below (lowstand).
      * Sea floor / surface elevation: typically negative (below present
        sea level); briefly positive if the shelf was exposed.
      * Tectonic subsidence (basement): large negative (deep).
      * Dynamic topography anomaly: positive = uplift relative to today.
      * Stratigraphic-unit boundaries: stack downward from the sea floor
        (each younger-to-older transition subtracts the unit's decompacted
        thickness), so all unit elevations are negative.
    """
    n_units = len(well.stratigraphic_units)
    ages, sea_levels, surface_elevations = [], [], []
    tect_elevations, dyn_topo = [], []
    all_unit_elevations = [[] for _ in range(n_units)]

    for ds in decomp:
        age = ds.get_age()
        ages.append(age)

        # Sea level: native positive-up sign from pybacktrack.
        sea_level = ds.get_sea_level() if hasattr(ds, "get_sea_level") else 0.0
        sea_levels.append(sea_level)

        # Sea-floor elevation = sea level (positive up) - water depth
        # (positive down).  Result is negative for any submerged site.
        water_depth = ds.get_water_depth()
        surface_elev = sea_level - water_depth
        surface_elevations.append(surface_elev)

        # Tectonic subsidence: pybacktrack returns it positive-down, so
        # we negate to get the basement *elevation* (a large negative).
        tect_elevations.append(-ds.get_tectonic_subsidence())

        # Dynamic topography: pybacktrack samples as elevation
        # (positive = uplift).  Keep that sign.
        if dt_model is not None:
            dt_val = dt_model.sample(float(age))
            dt = (0.0 if (dt_val is None or np.isnan(dt_val))
                  else float(dt_val))
        else:
            dt = 0.0
        dyn_topo.append(dt)

        # Stratigraphic-unit cumulative elevations: stack downward from
        # the sea floor by subtracting each decompacted unit thickness.
        # Already-stripped units are pinned to the current sea floor.
        num_still_present = len(ds.decompacted_stratigraphic_units)
        ui = 0
        while ui < n_units - num_still_present:
            all_unit_elevations[ui].append(surface_elev)
            ui += 1
        cum_elev = surface_elev
        for ds_unit in ds.decompacted_stratigraphic_units:
            cum_elev -= ds_unit.decompacted_thickness
            all_unit_elevations[ui].append(cum_elev)
            ui += 1

    return (
        np.array(ages),
        np.array(sea_levels),
        np.array(surface_elevations),
        np.array(tect_elevations),
        np.array(dyn_topo),
        [np.array(u) for u in all_unit_elevations],
    )


# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------
def plot_geohistory(name, ages, sea, surface, tect, dynt, unit_elevations,
                    out_path):
    """Single-panel matplotlib geohistory plot.

    Time axis: 150 Ma on the left, 0 Ma on the right.
    Y axis: **positive UP** (standard geological convention).  Sea-level
    highstands above zero; sea floor and basement below.  Default
    window is +300 m (top) to -4500 m (bottom), auto-extended on either
    end to fit the actual data.
    """
    # Auto-extend the default window if the data exceeds it.
    default_top = 300.0
    default_bot = -4500.0
    # Highest visible value: include sea-level highstand + a little
    # headroom; include surface-elevation excursions too (rare but
    # possible during regression).
    candidates_top = [default_top, float(sea.max()) + 50.0,
                      float(surface.max()) + 50.0]
    if np.any(dynt != 0):
        candidates_top.append(float((dynt - dynt[0]).max()) + 50.0)
    y_top = max(candidates_top)
    # Lowest visible value: include basement + a little headroom.
    deepest = float(tect.min())
    if unit_elevations and len(unit_elevations[-1]) > 0:
        deepest = min(deepest, float(unit_elevations[-1].min()))
    y_bot = min(default_bot, deepest - 200.0)
    # Snap to round multiples for clean tick spacing.
    y_top = float(np.ceil(y_top / 100.0) * 100.0)
    y_bot = float(np.floor(y_bot / 500.0) * 500.0)

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.set_xlim(MAX_ANALYSIS_AGE, 0)
    ax.set_ylim(y_bot, y_top)                       # NOT inverted

    # Fill stratigraphic units oldest-first so the youngest end up on top.
    # In the new positive-up convention each unit's "upper" boundary is at
    # a more positive elevation than its "lower" boundary.
    n_units = len(unit_elevations)
    for i in reversed(range(n_units)):
        upper = unit_elevations[i - 1] if i > 0 else surface
        lower = unit_elevations[i]
        col = UNIT_COLOURS[i % len(UNIT_COLOURS)]
        ax.fill_between(ages, lower, upper,
                        color=col, edgecolor="gray", linewidth=0.3)

    # Reference curves.
    ax.plot(ages, surface, color="royalblue", lw=1.5, label="Sea floor")
    ax.plot(ages, sea, color="navy", lw=1.0, ls="--",
            label=f"Sea level ({SEA_LEVEL_MODEL})")
    ax.plot(ages, tect, color="red", lw=2.0,
            label="Tectonic subsidence (basement)")
    if np.any(dynt != 0):
        ax.plot(ages, dynt - dynt[0],
                color="orange", lw=1.25, ls=":",
                label=f"{DYNAMIC_TOPOGRAPHY_MODEL} DT anomaly")

    # Horizontal zero line: marks present-day sea level.  Useful visual
    # anchor now that the y-axis crosses zero.
    ax.axhline(0, color="black", lw=0.7, alpha=0.6)

    # Font sizes (user 2026-06-02): bumped up to match Fig 4.
    AXIS_LABEL_SIZE = 16
    TICK_LABEL_SIZE = 14
    LEGEND_SIZE     = 13

    ax.set_xlabel("Age (Ma)", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel("Depth below present-day sea level (m)",
                  fontsize=AXIS_LABEL_SIZE)
    ax.tick_params(axis="both", labelsize=TICK_LABEL_SIZE)
    # Title intentionally removed -- the manuscript figure caption
    # carries the well name and full attribution.
    ax.grid(True, alpha=0.3, linestyle=":")
    # Legend on the LOWER LEFT.  In the new positive-up convention, the
    # ages 150-146 Ma (left edge of the plot, before asteras' first
    # stratigraphic horizon at 146 Ma) are still empty AND the lower
    # half of the panel is well below any active curve / stratigraphic
    # polygon, so this corner doesn't overlap any data.
    ax.legend(loc="lower left", framealpha=0.9, fontsize=LEGEND_SIZE)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    fig.savefig(os.path.splitext(out_path)[0] + ".pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    if not WELLS:
        sys.exit("No wells configured in config.WELLS")
    for w in WELLS:
        print(f"\n=== Geohistory analysis: {w['name']} ({w['file']}) ===")
        well, decomp = pybacktrack.backstrip_well(
            w["file"],
            lithology_filenames=pybacktrack.BUNDLE_LITHOLOGY_FILENAMES,
            sea_level_model=SEA_LEVEL_MODEL,
        )
        dt = pybacktrack.DynamicTopography.create_from_bundled_model(
            DYNAMIC_TOPOGRAPHY_MODEL, well.longitude, well.latitude
        )
        ages, sea, surface, tect, dynt, units = get_decompacted_history(
            well, decomp, dt)

        slug = w["name"].lower().replace(" ", "_")
        out = os.path.join(OUTPUT_DIR, f"fig08_geohistory_{slug}.png")
        plot_geohistory(w["name"], ages, sea, surface, tect, dynt, units, out)


if __name__ == "__main__":
    main()
