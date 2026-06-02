#!/usr/bin/env python3
"""
Figure 5 - NW Shelf well location map.

Vertical-gravity-gradient basemap (Sandwell et al. v31.1, +/-60
Eotvos colour range, Crameri `vik` diverging colourmap) overlaid
with the 109 industry wells used in the backstripping and gridded
rate-map analyses.  The featured comprehensive-example well, asteras
(used in Figs 5-7), is emphasised with a red square and label; all
other wells are small black dots.

The VGG basemap highlights short-wavelength changes in subsurface
density and exposes the NE-SW Mesozoic rift trends, depocentre
boundaries and transform-margin lineaments that motivate the
distribution of NW Shelf industry drilling.

Output:
    figures/output/fig05_nwshelf_map.png
    figures/output/fig05_nwshelf_map.pdf
"""
import glob
import os
import sys

import numpy as np
import pygmt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import WELLS, NWSHELF_WELL_DIR, OUTPUT_DIR, VGG_GRID

# ----------------------------------------------------------------------------
# Tunables
# ----------------------------------------------------------------------------
# Map region (lon_min, lon_max, lat_min, lat_max).  Framed around the
# active NW Shelf depocentre cluster and the featured well (asteras).
# Keep this in sync with VGG_REGION in tools/populate_data.py -- that's
# the window the regional VGG grid was cropped to.
REGION = [114, 130, -20, -9]

# Symmetric colour range, capped to keep the dynamic range readable
# (the global VGG distribution has long tails).  60 E (was 80) on
# 2026-05-28: lets the inner-shelf structural fabric stand out more
# clearly at the tighter [114, 130, -20, -9] window.
VGG_CLIP = 60.0                          # Eotvos
VGG_CMAP = "vik"                         # Crameri diverging blue-white-red
VGG_SHADING = "+a45+nt0.5"               # gentle hillshade for relief feel

# Highlighted-well markers and labels.
# Plain Helvetica, non-bold, black -- no halo, no red fill -- so the
# label sits cleanly beside the smaller red square without competing
# with the VGG colour-bar reading.
HIGHLIGHT_FILL = "red"
HIGHLIGHT_PEN = "0.75p,black"
HIGHLIGHT_MARKER_SIZE = "0.28c"          # square side; was 0.45c
HIGHLIGHT_LABEL_FONT = "14p,Helvetica,black"

# Geographic labels: physiographic features (basins, abyssal plain, sea)
# in italic, land features (country, island) in bold roman.  Each label
# gets a semi-transparent white box backdrop for legibility over the
# colourful VGG basemap.  Coordinates are tuned to (a) sit on the
# correct feature in the basemap, (b) avoid the asteras well marker at
# (124.12, -13.15), and (c) stay clear of the densely-drilled Browse
# Basin core cluster around (123-125 E, 13-14 S).
GEO_LABELS = [
    # (text, lon, lat, font, justify)
    ("Northern Carnarvon Basin", 115.3, -18.0, "11p,Helvetica-Oblique,black", "ML"),
    ("Browse Basin",             120.8, -14.2, "12p,Helvetica-Oblique,black", "MC"),
    ("Bonaparte Basin",          128.0, -12.5, "11p,Helvetica-Oblique,black", "MC"),
    ("Argo Abyssal Plain",       117.0, -13.0, "11p,Helvetica-Oblique,black", "MC"),
    ("Arafura Sea",              129.4, -10.0, "11p,Helvetica-Oblique,black", "MR"),
    ("AUSTRALIA",                126.7, -18.5, "14p,Helvetica-Bold,black",    "MC"),
    ("Timor",                    124.5, -9.35, "12p,Helvetica-Bold,black",    "MC"),
]


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def read_all_wells():
    """Walk NWSHELF_WELL_DIR and return a list of (key, lon, lat) tuples.

    `key` is the well filename without its extension -- matches the
    `file=` paths in config.WELLS, which is how the highlighted-set
    lookup is done.
    """
    out = []
    for f in (sorted(glob.glob(os.path.join(NWSHELF_WELL_DIR, "*.txt"))) +
              sorted(glob.glob(os.path.join(NWSHELF_WELL_DIR, "*.dat")))):
        lat = lon = None
        with open(f) as fh:
            for _ in range(10):
                line = fh.readline()
                if not line:
                    break
                if "SiteLatitude" in line:
                    try: lat = float(line.split("=")[-1])
                    except ValueError: pass
                if "SiteLongitude" in line:
                    try: lon = float(line.split("=")[-1])
                    except ValueError: pass
        if lat is not None and lon is not None:
            key = os.path.splitext(os.path.basename(f))[0]
            out.append((key, lon, lat))
    return out


def find_vgg_grid():
    """Resolve the regional VGG basemap path (see config.VGG_GRID).

    Shipped in-repo as ``data/grids/vgg_nwshelf.nc`` after running
    ``tools/populate_data.py``; raise informatively if missing.
    """
    if os.path.exists(VGG_GRID):
        return VGG_GRID
    raise FileNotFoundError(
        f"Could not locate the regional VGG grid: {VGG_GRID}\n"
        "Run `python tools/populate_data.py` from the repo root to "
        "generate it from the global Sandwell V31.1 source."
    )


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    wells = read_all_wells()
    print(f"Read {len(wells)} wells from {NWSHELF_WELL_DIR}")

    highlighted = {}
    for cfg in WELLS:
        key = os.path.splitext(os.path.basename(cfg["file"]))[0]
        highlighted[key] = cfg["name"]
    print(f"Highlighting: {highlighted}")

    vgg_path = find_vgg_grid()
    print(f"VGG basemap:  {vgg_path}")

    fig = pygmt.Figure()
    pygmt.config(FONT_TITLE="14p,Helvetica-Bold",
                 FONT_LABEL="11p", FONT_ANNOT="10p",
                 COLOR_NAN="220/220/220")

    # Symmetric diverging VGG palette, capped at +/- VGG_CLIP Eotvos.
    pygmt.makecpt(cmap=VGG_CMAP,
                  series=[-VGG_CLIP, VGG_CLIP, 5.0],
                  continuous=True)

    fig.basemap(
        region=REGION, projection="M16c",
        frame=['WSne',
               "xa5f1+lLongitude",
               "ya5f1+lLatitude"],
    )
    # VGG basemap.  Cropped to REGION by pyGMT automatically.
    # NB: GMT 6.5 rejects `nan_transparent=True` combined with `shading=`
    # ("Cannot specify a transparent color for grids when intensities are
    # also used") -- and crucially the rejection is SILENT for grdimage:
    # the basemap simply doesn't draw, but the colourbar below still does.
    # Drop nan_transparent; NaN cells (if any) take COLOR_NAN (light grey).
    fig.grdimage(
        grid=vgg_path,
        region=REGION,
        cmap=True,
        shading=VGG_SHADING,
    )
    # Coastlines and country borders, drawn on top of the VGG so the
    # land-water split is still legible.
    fig.coast(shorelines="0.4p,black", resolution="i",
              borders="2/0.25p,gray40")

    # All wells: small black dots.
    lons = np.array([w[1] for w in wells])
    lats = np.array([w[2] for w in wells])
    fig.plot(x=lons, y=lats, style="c0.15c",
             fill="black", pen="0.4p,white")

    # Highlighted well(s): small red square + plain Helvetica label.
    for key, lon, lat in wells:
        if key in highlighted:
            fig.plot(x=[lon], y=[lat], style=f"s{HIGHLIGHT_MARKER_SIZE}",
                     fill=HIGHLIGHT_FILL, pen=HIGHLIGHT_PEN)
            fig.text(
                x=lon, y=lat,
                text=highlighted[key],
                font=HIGHLIGHT_LABEL_FONT,
                justify="ML",
                offset="0.25c/0.0c",
            )

    # Geographic labels (basins, abyssal plain, sea, land/island).
    # Semi-transparent white box backdrop keeps them legible over the
    # VGG colour swing.
    for text, lon, lat, font, justify in GEO_LABELS:
        fig.text(
            x=lon, y=lat,
            text=text,
            font=font,
            justify=justify,
            fill="white@30",
            clearance="0.08c/0.05c",
        )

    # Colourbar for VGG.  Annotation every 20 E, finer ticks every 5 E.
    fig.colorbar(
        frame=["a20f5", "x+lVertical gravity gradient (Eotvos)"],
        position="JBC+w12c/0.35c+h+o0/1.3c+ma",
    )

    base = os.path.join(OUTPUT_DIR, "fig05_nwshelf_map")
    fig.savefig(base + ".png", dpi=300)
    fig.savefig(base + ".pdf")
    print(f"\nwrote {base}.png")
    print(f"wrote {base}.pdf")


if __name__ == "__main__":
    main()
