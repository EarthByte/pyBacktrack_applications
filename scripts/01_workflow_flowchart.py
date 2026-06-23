#!/usr/bin/env python3
"""
Figure 1 - pyBacktrack 1.5 paleobathymetry-gridding workflow flowchart.

Renders the data flow of the `paleo_bathymetry` module as a labelled
node-and-arrow diagram:

    Input data (peach)     -> present-day age, sediment thickness,
                              bathymetry; rift start/end age grids;
                              plate model; dynamic-topography model;
                              sea-level curve.
    Processing steps (blue) -> decompaction; oceanic/continental
                              subsidence; plate reconstruction; grid
                              interpolation.
    Decision (yellow)      -> per-grid-point oceanic-vs-continental
                              branch based on the age grid.
    External input (lavender) -> optional traditional-paleobathymetry
                              input for subducted crust (merging).
    Output (green)         -> per-Myr NetCDF paleobathymetry grids.

Canvas is 12 x N units of axes coordinates; box sizes are tuned so
that text fits cleanly inside every rounded rectangle.

Source material:
    pyBacktrack-master/docs/pybacktrack_paleo_bathymetry.rst
    pyBacktrack-master/pybacktrack/paleo_bathymetry.py
    pyBacktrack-master/pybacktrack/notebooks/paleobathymetry.ipynb

Output:
    figures/output/fig01_paleobathy_workflow.png
    figures/output/fig01_paleobathy_workflow.pdf   (vector, for editing)
    figures/output/fig01_paleobathy_workflow.svg   (vector, for Illustrator)
"""
import os
import sys

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import OUTPUT_DIR


# Colours by node category.
COLOURS = {
    "input":    "#FFE5B4",   # peach - input data
    "process":  "#C8E1FF",   # light blue - core algorithm step
    "decision": "#FFEFA1",   # yellow - branching point
    "output":   "#C4EFC4",   # green - final product
    "ext":      "#E1D5FF",   # lavender - optional external input
}


def box(ax, x, y, w, h, label, kind="process", fontsize=9):
    """Draw a rounded rectangle with a centred label."""
    p = FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.06",
        linewidth=1.2, edgecolor="black",
        facecolor=COLOURS.get(kind, "white"),
        zorder=2,
    )
    ax.add_patch(p)
    ax.text(
        x, y, label, ha="center", va="center",
        fontsize=fontsize, zorder=3,
    )


def arrow(ax, x1, y1, x2, y2, label=None, label_offset=(0.0, 0.0)):
    a = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle="-|>", mutation_scale=12,
        linewidth=1.0, color="black", zorder=1,
    )
    ax.add_patch(a)
    if label:
        ax.text(
            0.5 * (x1 + x2) + label_offset[0],
            0.5 * (y1 + y2) + label_offset[1],
            label, fontsize=8, ha="center", va="center", style="italic",
            zorder=4,
        )


def main():
    fig, ax = plt.subplots(figsize=(12, 13), dpi=200)
    # ylim extends below 0 so the bottom-row "output" box (centred at
    # y=0.3 with h=0.8 -> bottom edge at y=-0.1) is fully inside the
    # axes rather than clipped at the bottom, AND so there is room
    # underneath the diagram for the category legend.
    ax.set_xlim(0, 12); ax.set_ylim(-1.6, 14)
    ax.set_aspect("equal"); ax.axis("off")

    # Central column x-coordinate (was 5.0 on the old narrow canvas).
    CX = 6.0

    # -------------------- Input layer (top row) --------------------
    box(ax, 2.0, 13, 3.0, 0.9,
        "Present-day age grid\n(Seton et al., 2020)", kind="input")
    box(ax, 6.0, 13, 3.0, 0.9,
        "Total sediment thickness\n(GlobSed; Straume et al., 2019)",
        kind="input", fontsize=8)
    box(ax, 10.0, 13, 3.0, 0.9,
        "Present-day bathymetry\n(ETOPO1)", kind="input")

    box(ax, 2.0, 11.7, 3.0, 0.9,
        "Built-in rift start/end\nage grids", kind="input")
    box(ax, 6.0, 11.7, 3.0, 0.9,
        "Plate model\n(rotations, static polygons)",
        kind="input", fontsize=9)
    box(ax, 10.0, 11.7, 3.0, 0.9,
        "Dynamic topography model",
        kind="input", fontsize=10)

    # Row 3 -- widened (was w=4.0) so the label no longer touches the
    # rounded-rectangle edges.
    box(ax, CX, 10.5, 6.0, 0.8,
        "Sea-level curve",
        kind="input", fontsize=10)

    # -------------------- Per-grid-point loop (middle) --------------------
    # Row 4 -- widened (was w=5.0) and h bumped from 0.7 to 0.8.
    box(ax, CX, 9.3, 7.0, 0.8,
        "For every present-day grid point (default 6 arc-min spacing) ...",
        kind="process", fontsize=10)

    # Row 5 -- decision + sediment thickness sampler.
    box(ax, 2.6, 8.0, 3.6, 0.9,
        "Inside age grid?\n(oceanic vs continental crust)",
        kind="decision", fontsize=9)

    box(ax, 8.8, 8.0, 3.4, 0.9,
        "Sample sediment thickness\n+ single base lithology\n(Average ocean-floor sed.)",
        kind="process", fontsize=8)

    # Row 6 -- THREE process boxes.
    # Oceanic box was at x=0.7 (left edge -0.3, OFF the canvas).  Now at
    # x=1.3 (left edge 0.3, well inside).
    box(ax, 1.3, 6.6, 2.0, 1.0,
        "Oceanic\nsubsidence:\nage-depth model\n(RHCW18 / GDH1)",
        kind="process", fontsize=8)
    box(ax, 4.1, 6.6, 2.8, 1.0,
        "Continental subsidence:\nMcKenzie stretching\nfrom rift start/end ages",
        kind="process", fontsize=8)
    box(ax, 8.8, 6.6, 3.4, 1.0,
        "Decompact sediment column\nbackward through time\n(porosity-depth)",
        kind="process", fontsize=8)

    # Row 7 -- combine.  Widened (was w=5.0, h=0.9) so the longer of the
    # two text lines no longer overruns the box edge.
    box(ax, CX, 5.1, 7.0, 1.0,
        "Combine tectonic subsidence + decompacted thickness\n"
        "+ dynamic topography + sea-level correction",
        kind="process", fontsize=10)

    # Row 8 -- reconstruct.  Widened to match row 7.
    box(ax, CX, 3.8, 7.0, 1.0,
        "Reconstruct each grid point to its paleo-position\n"
        "using rotation file + static polygons",
        kind="process", fontsize=10)

    # Row 9 -- interpolate (centre) + traditional bathymetry (left).
    # Interpolate widened to w=6.5 (was 5.0).  Traditional box moved a
    # touch further left and slightly narrower so the two no longer
    # visually overlap.
    box(ax, 1.5, 2.6, 2.6, 0.9,
        "Forward-modelled\nbathymetry on\nsubducted crust",
        kind="ext", fontsize=8)
    box(ax, 6.5, 2.6, 6.5, 0.9,
        "Interpolate (paleo-lon, paleo-lat, paleo-depth) tuples\n"
        "onto a regular output grid at every output time",
        kind="process", fontsize=10)

    # Row 10 -- merge.
    box(ax, CX, 1.4, 7.0, 1.0,
        "Optionally merge: pyBacktrack on extant crust\n"
        "+ forward-modelled bathymetry on subducted crust",
        kind="process", fontsize=10)

    # -------------------- Output --------------------
    # Row 11 -- widened (was w=5.0) and h bumped 0.6 -> 0.8.
    box(ax, CX, 0.3, 7.0, 0.8,
        "Time-dependent paleobathymetry NetCDF grids (per Myr)",
        kind="output", fontsize=11)

    # -------------------- Arrows --------------------
    # Every input box feeds the per-grid-point loop directly.  Earlier
    # revisions had the row-2 inputs converge on the sea-level-curve
    # box, which incorrectly implied they were inputs to the sea-level
    # model (John Cannon review, 2026-06).  Each of the seven input
    # boxes now arrows straight into a distinct x along the top edge
    # of the loop-entry blue box ("For every present-day grid point
    # ..."), spanning x in [2.5, 9.5] at y = 9.70.  Arrow lines pass
    # behind intermediate boxes (peach facecolors are opaque, box
    # zorder=2 > arrow zorder=1), so the diagram reads as a clean
    # funnel of inputs converging on the loop entry.
    LOOP_TOP_Y = 9.70
    # Row 1 (y=13, h=0.9 -> bottom 12.55)
    arrow(ax, 2.0, 12.55, 3.0, LOOP_TOP_Y)   # Age grid
    arrow(ax, 6.0, 12.55, 5.5, LOOP_TOP_Y)   # Sediment thickness
    arrow(ax, 10.0, 12.55, 8.5, LOOP_TOP_Y)  # Bathymetry
    # Row 2 (y=11.7, h=0.9 -> bottom 11.25)
    arrow(ax, 2.0, 11.25, 3.5, LOOP_TOP_Y)   # Rift age grids
    arrow(ax, 6.0, 11.25, 6.0, LOOP_TOP_Y)   # Plate model
    arrow(ax, 10.0, 11.25, 9.0, LOOP_TOP_Y)  # Dynamic topography
    # Row 3 (y=10.5, h=0.8 -> bottom 10.10)
    arrow(ax, CX, 10.10, 6.5, LOOP_TOP_Y)    # Sea-level curve

    # Row 4 loop -> branching (decision + sediment sample)
    arrow(ax, CX, 8.90, 2.6, 8.50)
    arrow(ax, CX, 8.90, 8.8, 8.50)

    # Decision -> branches
    arrow(ax, 2.3, 7.55, 1.5, 7.10, label="yes",
          label_offset=(-0.55, 0.10))
    arrow(ax, 2.9, 7.55, 3.9, 7.10, label="no",
          label_offset=(0.55, 0.10))

    # Sample sediment thickness (row 5 R) -> Decompact sediment column
    # (row 6 R).  Was missing in earlier revisions (John Cannon review,
    # 2026-06); without this, Decompact appeared to start from nowhere
    # and the sediment-sampling step had no visible outflow.
    arrow(ax, 8.8, 7.55, 8.8, 7.10)

    # Row 6 -> combine (row 7)
    arrow(ax, 1.3, 6.10, 4.0, 5.65)
    arrow(ax, 4.1, 6.10, 5.5, 5.65)
    arrow(ax, 8.8, 6.10, 7.5, 5.65)

    # Combine -> reconstruct -> merge -> output, centre column
    arrow(ax, CX, 4.60, CX, 4.30)        # row 7 -> row 8
    arrow(ax, CX, 3.30, CX, 3.05)        # row 8 -> row 9 (interpolate)
    arrow(ax, CX, 2.15, CX, 1.90)        # row 9 -> row 10 (merge)
    arrow(ax, 2.8, 2.6, 3.5, 1.90)       # traditional bathy -> merge
    arrow(ax, CX, 0.90, CX, 0.65)        # merge -> output

    # -------------------- Legend (UNDER the diagram) --------------------
    # Placed beneath the bottom output box.  Anchored at axes-x=0.5
    # (horizontal centre) and axes-y just inside the extended lower
    # ylim, with loc="upper center" so the legend's TOP edge sits at
    # the anchor and the legend grows downward into the spare ylim
    # space.
    legend_items = [
        ("Input data", "input"),
        ("Process step", "process"),
        ("Decision", "decision"),
        ("Optional external", "ext"),
        ("Output", "output"),
    ]
    handles = [
        mpatches.Patch(facecolor=COLOURS[kind], edgecolor="black",
                       label=label)
        for label, kind in legend_items
    ]
    ax.legend(
        handles=handles, loc="upper center",
        bbox_to_anchor=(0.5, 0.06), ncol=5,
        frameon=False, fontsize=11,
    )

    # No suptitle: the manuscript's figure caption supplies the title;
    # an in-figure title was duplicating that and pushing the diagram
    # down the page (Dietmar 2026-05-31).

    # ---- Output: PNG (raster) + PDF + SVG (both vector) ----
    # Tight crop -- pad_inches=0.02 keeps just a hairline of whitespace
    # around the diagram so the box outlines aren't flush with the
    # canvas edge, but no more.  Default pad_inches=0.1 was leaving
    # noticeable margins (Dietmar 2026-05-31).
    base = os.path.join(OUTPUT_DIR, "fig01_paleobathy_workflow")
    fig.savefig(base + ".png", dpi=300, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(base + ".pdf",          bbox_inches="tight", pad_inches=0.02)
    fig.savefig(base + ".svg",          bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print(f"wrote {base}.png")
    print(f"wrote {base}.pdf")
    print(f"wrote {base}.svg")


if __name__ == "__main__":
    main()
