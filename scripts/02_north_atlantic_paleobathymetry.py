#!/usr/bin/env python3
"""
Figure 2 - Lambert-azimuthal-equal-area paleobathymetry centred in
the mid-North Atlantic, at 120, 80, 60 and 40 Ma.

Default centre (CENTRAL_LON=0, CENTRAL_LAT=45) with a HORIZON_DEG=50
cap covers the western Tethys (Mediterranean and Caribbean realm),
the full North Atlantic, the Arctic Ocean (well past the pole) and
extends south to roughly the equator.  Equal-area, so ocean basins
can be compared visually without size distortion.

Rendering helpers (gplately driver, paleobathymetry rendering with
hillshade, plate-model overlay) come from `paleobathy_render.py`.

Output:
    figures/output/fig02_paleobathymetry_<time>Ma.png
    figures/output/fig02_paleobathymetry_4panel.png
    (and .pdf vector siblings)
"""
import os
import sys

import pygmt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import OUTPUT_DIR  # noqa

import paleobathy_render as pbr

# ============================================================================
# USER-TUNABLE HILLSHADE CONTROLS
# ----------------------------------------------------------------------------
# Edit these to change the bathymetry relief shading.  Applied to the
# shared rendering module by _apply_hillshade_settings() in main().
#
# Plumbing: pbr.gradient_for() calls
#     pygmt.grdgradient(azimuth=HILLSHADE_AZIMUTH,
#                       normalize=HILLSHADE_NORMALIZE,
#                       outgrid=<cached grad_<N>Ma.nc>)
# and pbr.draw_paleobathymetry() then passes the cached .nc as
# `shading=...` to fig.grdimage.  Changing HILLSHADE_NORMALIZE without
# wiping the cache is a silent no-op for ages whose gradient is already
# on disk -- _apply_hillshade_settings() takes care of the wipe.
# ============================================================================
HILLSHADE_ON        = True       # set False to disable shading entirely
HILLSHADE_AZIMUTH   = "315"      # deg clockwise from north (315 = NW)
HILLSHADE_NORMALIZE = "e0.3+a0.3"  # soft preset: gentle exponential + 30 %
                                  # ambient illumination.  Same as the
                                  # video so the static panels and the
                                  # animation look consistent.


def _apply_hillshade_settings():
    """Push the user-tunable HILLSHADE_* constants above into the
    paleobathy_render module and wipe its cached gradients for this
    figure's TIMES_MA so the new settings actually take effect on
    the very next pygmt.grdgradient call.

    Set PYBT_KEEP_HILLSHADE_CACHE=1 in the environment to SKIP the
    wipe step -- useful when you've only changed the CPT, fonts, or
    overlay settings and want to keep the (slow) hillshade gradients
    from the last run.
    """
    pbr.HILLSHADE_ON        = HILLSHADE_ON
    pbr.HILLSHADE_AZIMUTH   = HILLSHADE_AZIMUTH
    pbr.HILLSHADE_NORMALIZE = HILLSHADE_NORMALIZE
    if os.environ.get("PYBT_KEEP_HILLSHADE_CACHE"):
        print("  PYBT_KEEP_HILLSHADE_CACHE set -- keeping existing "
              f"gradient cache for ages {TIMES_MA}")
        return
    pbr.wipe_gradient_cache(TIMES_MA)


# Convenience re-exports so the rest of the script reads cleanly.
make_gplately_plot   = pbr.make_gplately_plot
_draw_paleobathymetry = pbr.draw_paleobathymetry
BATHY_CMAP           = pbr.BATHY_CMAP
BATHY_CPT_SERIES     = pbr.BATHY_CPT_SERIES
NAN_COLOR            = pbr.NAN_COLOR
FONT_TITLE           = pbr.FONT_TITLE
FONT_LABEL           = pbr.FONT_LABEL
FONT_ANNOT           = pbr.FONT_ANNOT

# ----------------------------------------------------------------------------
# Colour palette: GMT's built-in `bathy` (sequence of blues, classical
# oceanographic bathymetry look).  Same palette is now used by Fig 3
# (south-pole) and the deprecated global Fig 2.  Rebinding here pins
# the cpt locally even if `02_global_paleobathymetry.py:BATHY_CMAP`
# changes later.
# ----------------------------------------------------------------------------
BATHY_CMAP = "bathy"


# Times for this figure -- the Cretaceous-Cenozoic interval chosen to
# show the wide-open Cretaceous Tethys (120 Ma), the maximum extent of
# the Late Cretaceous oceans (80 Ma), the early Cenozoic Atlantic
# opening (60 Ma) and the mid-Eocene establishment of the modern
# Southern Ocean (40 Ma).  Edit if a different time list is more
# illuminating for the Atlantic / Arctic / Tethys story.
TIMES_MA = [80, 60, 40, 20]

# Lambert Azimuthal Equal Area centred at (CENTRAL_LON, CENTRAL_LAT)
# with a circular cap of HORIZON_DEG angular radius.  The defaults
# below put the centre in the mid-North Atlantic at ~45 deg N so the
# visible cap covers:
#   - the entire western Tethys (Mediterranean + Caribbean realm),
#   - the full North Atlantic ocean,
#   - the Arctic Ocean (well past the pole),
#   - extends south to roughly the equator,
#   - and reaches west into the eastern Pacific / east into Asia.
# Equal-area, so ocean basins can be compared visually without size
# distortion.
#
# Tune these to taste:
#   - move CENTRAL_LON east (e.g. 30) to shift toward Eurasia / Indian
#     Ocean; west (e.g. -90) to shift toward the Pacific
#   - lower CENTRAL_LAT (e.g. 30) to weight more toward the tropics,
#     raise (e.g. 60) to weight more toward the Arctic
#   - lower HORIZON_DEG (e.g. 40) for a tighter cap; raise (e.g. 60)
#     for closer to the full hemisphere
#
# Angular alternative: if a non-round fan / rectangular layout reads
# better, swap the projection in plot_one_time / plot_four_panel to
# either:
#   - Albers equal-area conic:
#       projection = f"B{CENTRAL_LON}/{CENTRAL_LAT}/20/70/{w}c"
#       region     = [-120, 120, 0, 85]
#   - Lambert conformal conic:
#       projection = f"L{CENTRAL_LON}/{CENTRAL_LAT}/20/70/{w}c"
#       region     = [-120, 120, 0, 85]
# Both give a fan-shaped map with straight meridians and curved
# parallels (rectangular at standard parallels, narrowing toward the
# pole), good for showing N Pacific + N Atlantic + Arctic + equator
# side-by-side with minimal high-latitude distortion.
CENTRAL_LON = -20.0      # mid-Atlantic centre (geographic identity:
                          #   340 also valid; lat/lon labels are drawn
                          #   manually below so the sign no longer
                          #   affects rendering)
CENTRAL_LAT = 45.0
HORIZON_DEG = 40.0
PROJ_WIDTH_CM_SINGLE = 14
PROJ_WIDTH_CM_PANEL = 14


# ============================================================================
# Manual lon/lat labels.
# ----------------------------------------------------------------------------
# GMT 6.5's `-B WSne` axis annotations DO NOT render on the circular
# boundary of an oblique LAEA projection (verified by side-by-side
# comparison with Fig 3, where the same code path with a polar centre
# DOES produce labels).  The frame box, fancy ring, age tag and
# graticule all draw fine -- only the WSne axis annotations are
# silently dropped.
#
# Workaround: compute label positions GEOMETRICALLY from the cap
# boundary.  The cap boundary is the set of points at exactly
# HORIZON_DEG great-circle distance from (CENTRAL_LON, CENTRAL_LAT),
# so for any chosen meridian/parallel we can solve for the (lon, lat)
# at which that graticule line exits the cap.  Labels are placed at
# those exit points, which means they automatically follow the cap
# whenever CENTRAL_LON / CENTRAL_LAT / HORIZON_DEG are adjusted.
#
# Tunables:
#   LABEL_LON_MERIDIANS  -- list of (meridian_deg, lon_nudge, lat_nudge).
#                           Each lon label is independently tunable.
#   LABEL_LAT_PARALLELS  -- list of (parallel_deg, lon_nudge, lat_nudge).
#                           Each lat label is independently tunable.
# ============================================================================
import math as _math

_LABEL_FONT      = "16p,Helvetica,black"
_LABEL_FILL      = "white@30"
_LABEL_CLEARANCE = "0.06c/0.04c"

# Per-label tuning.  Each entry is a tuple:
#     (graticule_value_deg, lon_nudge_deg, lat_nudge_deg)
# where the nudges are added to the computed cap-edge intersection.
# Each label is INDEPENDENT -- editing one label's nudges has no
# effect on the others.  Drop a label by removing/commenting its
# tuple.  Add more labels by appending a tuple with the desired
# graticule value (the cap-edge computation handles any meridian /
# parallel that actually intersects the cap).
#
# Sign convention:
#     lon_nudge < 0  ->  shift WEST  (left)
#     lon_nudge > 0  ->  shift EAST  (right)
#     lat_nudge < 0  ->  shift SOUTH (down)
#     lat_nudge > 0  ->  shift NORTH (up)
#
# Meridian labels are anchored at the SOUTH arc; parallel labels at
# the WEST arc.

#  (meridian_deg,  lon_nudge,  lat_nudge)
LABEL_LON_MERIDIANS = [
    (-60,          -3.5,       -4.0),     # 60 deg W label
    (  0,           0.5,       -3.25),     # 0 deg    label
]

#  (parallel_deg,  lon_nudge,  lat_nudge)
LABEL_LAT_PARALLELS = [
    (40,           -2.25,        0.0),     # 40 deg N label
    (60,           -3.75,        0.0),     # 60 deg N label
]


def _cap_intersect_meridian(meridian_lon, southern=True):
    """Return the (lon, lat) where the meridian at `meridian_lon`
    crosses the LAEA cap boundary.  `southern=True` picks the southern
    of the two intersection points (suitable for "south edge" labels).
    Returns None if the meridian misses the cap entirely.
    """
    h  = _math.radians(HORIZON_DEG)
    cl = _math.radians(CENTRAL_LAT)
    dl = _math.radians(meridian_lon - CENTRAL_LON)
    # Solve cos(h) = sin(cl)*sin(lat) + cos(cl)*cos(lat)*cos(dl)
    a = _math.sin(cl)
    b = _math.cos(cl) * _math.cos(dl)
    R = _math.hypot(a, b)
    if R == 0:
        return None
    target = _math.cos(h) / R
    if abs(target) > 1:
        return None
    phi    = _math.atan2(b, a)
    angle  = _math.asin(target)
    lat_a  = _math.degrees(angle      - phi)
    lat_b  = _math.degrees(_math.pi - angle - phi)
    lat    = min(lat_a, lat_b) if southern else max(lat_a, lat_b)
    return (meridian_lon, lat)


def _cap_intersect_parallel(parallel_lat, western=True):
    """Return the (lon, lat) where the parallel at `parallel_lat`
    crosses the LAEA cap boundary.  `western=True` picks the western
    intersection.  Returns None if the parallel misses the cap.
    """
    h  = _math.radians(HORIZON_DEG)
    cl = _math.radians(CENTRAL_LAT)
    pl = _math.radians(parallel_lat)
    num = _math.cos(h) - _math.sin(cl) * _math.sin(pl)
    den = _math.cos(cl) * _math.cos(pl)
    if den == 0 or abs(num / den) > 1:
        return None
    dlon = _math.degrees(_math.acos(num / den))
    lon  = CENTRAL_LON - dlon if western else CENTRAL_LON + dlon
    return (lon, parallel_lat)


def _format_lon(deg):
    deg = ((deg + 180) % 360) - 180          # wrap to -180..180
    if deg < 0:
        return f"{int(round(-deg))}\\260W"
    if deg > 0:
        return f"{int(round(deg))}\\260E"
    return "0\\260"


def _format_lat(deg):
    if deg > 0:
        return f"{int(round(deg))}\\260N"
    if deg < 0:
        return f"{int(round(-deg))}\\260S"
    return "0\\260"


def _add_manual_lonlat_labels(fig):
    """Compute and draw geometrically-correct lon/lat labels at the
    cap boundary, with PER-LABEL nudges applied from the tuples
    LABEL_LON_MERIDIANS / LABEL_LAT_PARALLELS above.  Adapts
    automatically to CENTRAL_LON / CENTRAL_LAT / HORIZON_DEG.
    """
    # Lon labels at south arc of cap.
    for meridian, lon_nudge, lat_nudge in LABEL_LON_MERIDIANS:
        hit = _cap_intersect_meridian(meridian, southern=True)
        if hit is None:
            continue
        lon, lat = hit
        fig.text(x=lon + lon_nudge, y=lat + lat_nudge,
                 text=_format_lon(meridian), font=_LABEL_FONT,
                 justify="MC", fill=_LABEL_FILL,
                 clearance=_LABEL_CLEARANCE, no_clip=True)
    # Lat labels at west arc of cap.
    for parallel, lon_nudge, lat_nudge in LABEL_LAT_PARALLELS:
        hit = _cap_intersect_parallel(parallel, western=True)
        if hit is None:
            continue
        lon, lat = hit
        fig.text(x=lon + lon_nudge, y=lat + lat_nudge,
                 text=_format_lat(parallel), font=_LABEL_FONT,
                 justify="MR", fill=_LABEL_FILL,
                 clearance=_LABEL_CLEARANCE, no_clip=True)


def plot_one_time(time_ma, gplot):
    print(f"  rendering {time_ma} Ma")
    fig = pygmt.Figure()
    pygmt.config(
        FONT_TITLE=FONT_TITLE,
        FONT_LABEL="24p,Helvetica",          # bigger colorbar LABEL (user 2026-05-31)
        FONT_ANNOT="16p,Helvetica",          # lat/lon labels at panel rim
                                              #   (colorbar ticks re-bumped to
                                              #   18p just before fig.colorbar
                                              #   below)
        COLOR_NAN=NAN_COLOR,
        MAP_FRAME_TYPE="fancy",               # match Fig 3.  If labels still
                                              #   fail to render with the dash
                                              #   removed from CENTRAL_LON
                                              #   above, switch back to "plain"
                                              #   and we'll have isolated the
                                              #   cause to fancy + oblique LAEA.
        MAP_GRID_PEN_PRIMARY="1p,gray30",     # visible 60/30 deg graticule
    )
    # background=True clamps out-of-range cells (depths > 5000 m) to the
    # endpoint colour instead of rendering them white.  No output= here:
    # the downstream paleobathy_render.draw_paleobathymetry passes
    # cmap=True to grdimage, which means "use the session-current CPT",
    # and adding output= breaks that hand-off (grdimage falls back to
    # GMT's default rainbow palette).
    pygmt.makecpt(cmap=BATHY_CMAP, series=BATHY_CPT_SERIES,
                  continuous=True, background=True)

    # Thin lon/lat graticule: lines + labels every 60 deg lon, 30 deg
    # lat.  Boundary annotations via "WSne".  Lambert Azimuthal
    # Equal-Area projection (A) with a horizon angle that crops the
    # visible cap to a circle of HORIZON_DEG angular radius around
    # the centre.  `region="g"` is global; the horizon in the
    # projection string does the cropping.  The centred +t title was
    # dropped (user 2026-05-31): the age is now shown only as a
    # single big boxed tag in the upper-left of the panel.
    fig.basemap(
        region="g",
        projection=(f"A{CENTRAL_LON}/{CENTRAL_LAT}/{HORIZON_DEG}/"
                    f"{PROJ_WIDTH_CM_SINGLE}c"),
        # WSne brings the lon/lat annotations back to the outside of the
        # circular frame.  Fancy frame's 60 deg divisions follow the
        # xa60g60 annotation interval, not the WSne corners.
        frame=["xa60g60", "ya20g20", "WSne"],
    )
    _draw_paleobathymetry(fig, gplot, time_ma)
    # Re-draw the basemap ON TOP so the 60/20 deg graticule isn't
    # painted over by grdimage and the coastline polygons.
    fig.basemap(frame=["xa60g60", "ya20g20"])

    # Manual lat/lon labels (GMT WSne doesn't render on oblique LAEA --
    # see comment block at the top of this file).
    _add_manual_lonlat_labels(fig)

    # Big age tag in the top-left corner of the panel.
    fig.text(
        text=f"{int(time_ma)} Ma",
        position="TL", offset="0.3c/-0.05c", no_clip=True,
        font="28p,Helvetica-Bold,black",
        fill="white", pen="0.5p,black",
        clearance="0.18c/0.10c",
    )

    # Bump FONT_ANNOT back up just for the colorbar tick numbers (the
    # 14p above kept the lat/lon labels readable but small).
    pygmt.config(FONT_ANNOT="18p,Helvetica")
    fig.colorbar(
        frame=['a1000f200+lPaleobathymetry (m)'],
        position="JBC+w12c/0.4c+h+o0/1.6c+e+ma",
    )

    base = os.path.join(
        OUTPUT_DIR, f"fig02_paleobathymetry_{time_ma}Ma")
    fig.savefig(base + ".png", dpi=300)
    fig.savefig(base + ".pdf")          # vector
    print(f"  wrote {base}.png")


def plot_four_panel(gplot):
    print("\n  rendering 4-panel figure")
    fig = pygmt.Figure()
    pygmt.config(
        FONT_TITLE=FONT_TITLE,
        FONT_LABEL="24p,Helvetica",          # bigger colorbar LABEL
        FONT_ANNOT="16p,Helvetica",          # lat/lon labels at panel rim
                                              #   (colorbar ticks re-bumped to
                                              #   18p just before fig.colorbar
                                              #   below)
        COLOR_NAN=NAN_COLOR,
        MAP_FRAME_TYPE="fancy",               # match Fig 3.  If labels still
                                              #   fail to render with the dash
                                              #   removed from CENTRAL_LON
                                              #   above, switch back to "plain"
                                              #   and we'll have isolated the
                                              #   cause to fancy + oblique LAEA.
        MAP_GRID_PEN_PRIMARY="1p,gray30",     # visible 60/30 deg graticule
    )
    # background=True clamps out-of-range cells (depths > 5000 m) to the
    # endpoint colour instead of rendering them white.  No output= here:
    # the downstream paleobathy_render.draw_paleobathymetry passes
    # cmap=True to grdimage, which means "use the session-current CPT",
    # and adding output= breaks that hand-off (grdimage falls back to
    # GMT's default rainbow palette).
    pygmt.makecpt(cmap=BATHY_CMAP, series=BATHY_CPT_SERIES,
                  continuous=True, background=True)

    # Manual 2x2 layout via shift_origin instead of fig.subplot().
    # Reason: pyGMT subplot consistently DROPS the per-panel lat/lon
    # annotations on LAEA / circular projections, regardless of
    # WSne, frame, sharex/sharey, or margins (verified through
    # several iterations on 2026-05-31).  Treating each panel as an
    # independent basemap call keeps the WSne behaviour identical
    # to plot_one_time(), so the labels render correctly.
    PANEL_W = 14.0    # cm per panel
    GAP_X   = 2.5     # cm gap between columns
    GAP_Y   = 2.5     # cm gap between rows

    def _render(t):
        """Render one panel at the current origin."""
        fig.basemap(
            region="g",
            projection=(f"A{CENTRAL_LON}/{CENTRAL_LAT}/"
                        f"{HORIZON_DEG}/{PANEL_W}c"),
            frame=["xa60g60", "ya20g20", "WSne"],
        )
        _draw_paleobathymetry(fig, gplot, t)
        # Redraw 60/20 deg graticule on top of bathymetry.
        fig.basemap(frame=["xa60g60", "ya20g20"])
        # 4 manual lat/lon labels (GMT WSne doesn't render on oblique
        # LAEA -- see comment block at top of file).
        _add_manual_lonlat_labels(fig)
        fig.text(
            text=f"{int(t)} Ma",
            position="TL", offset="0.3c/-0.05c", no_clip=True,
            font="28p,Helvetica-Bold,black",
            fill="white", pen="0.5p,black",
            clearance="0.18c/0.10c",
        )

    # Layout: TIMES_MA = [oldest, ..., youngest] -> reading order
    # (a) top-left  (b) top-right
    # (c) bottom-left  (d) bottom-right
    # Render bottom row first (starting at default origin), then
    # shift up to render the top row, so shift_origin offsets are
    # always positive in y.
    _render(TIMES_MA[2])                                      # bottom-left
    fig.shift_origin(xshift=f"{PANEL_W + GAP_X}c")
    _render(TIMES_MA[3])                                      # bottom-right
    fig.shift_origin(xshift=f"-{PANEL_W + GAP_X}c",
                     yshift=f"{PANEL_W + GAP_Y}c")
    _render(TIMES_MA[0])                                      # top-left
    fig.shift_origin(xshift=f"{PANEL_W + GAP_X}c")
    _render(TIMES_MA[1])                                      # top-right

    # Now origin = top-right panel.  Shift back to bottom-left of the
    # layout so the colorbar position math is anchored cleanly.
    fig.shift_origin(xshift=f"-{PANEL_W + GAP_X}c",
                     yshift=f"-{PANEL_W + GAP_Y}c")

    # Colorbar centred under the 2-column layout.  Layout total width
    # is 2*PANEL_W + GAP_X; centre x = PANEL_W + GAP_X/2.  Colorbar
    # placed via x/y coords (relative to current origin = bottom-left
    # panel), 16 cm wide, anchored at its top-centre 1.6 cm below the
    # bottom-left panel.
    pygmt.config(FONT_ANNOT="18p,Helvetica")
    cb_centre_x = PANEL_W + GAP_X / 2
    fig.colorbar(
        frame=['a1000f200+lPaleobathymetry (m)'],
        position=(f"x{cb_centre_x}c/-1.6c+w16c/0.4c+h+e+ma+jTC"),
    )

    base = os.path.join(
        OUTPUT_DIR, "fig02_paleobathymetry_4panel")
    fig.savefig(base + ".png", dpi=300)
    fig.savefig(base + ".pdf")          # vector
    print(f"  wrote {base}.png")


def main():
    # Apply the user-tunable HILLSHADE_* constants from the top of this
    # file (and clear any cached gradients so the new settings take
    # effect on the very next pygmt.grdgradient call).
    _apply_hillshade_settings()

    # Lambert Azimuthal Equal-Area centred at lon=CENTRAL_LON -- pass
    # that to gplately so coastline / continent polygons get split at
    # the antimeridian and don't wrap as fans across the visible cap.
    gplot = make_gplately_plot(central_meridian=CENTRAL_LON)
    for t in TIMES_MA:
        plot_one_time(t, gplot)
    plot_four_panel(gplot)


if __name__ == "__main__":
    main()
