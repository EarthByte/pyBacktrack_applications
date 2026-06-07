#!/usr/bin/env python3
"""
Figure 3 - Lambert-azimuthal-equal-area paleobathymetry centred on the
South Pole, at 35, 25, 15 and 0 Ma.

Same projection family as Fig 2 (LAEA `A<lon>/<lat>/<horizon>/<width>`)
but centred at the South Pole.  HORIZON_DEG=50 puts the edge of the
visible cap at latitude -40 deg, i.e. the Southern Ocean and Antarctic
margin out to ~Cape Town / southern tip of South America / Tasmania.
Equal-area, so ocean basins (Weddell, Ross, Australian-Antarctic) can
be compared visually without size distortion.

Rendering helpers (gplately driver, paleobathymetry rendering with
hillshade, plate-model overlay) come from `paleobathy_render.py`.

Output:
    figures/output/fig03_southern_ocean_paleobathymetry_<time>Ma.png
    figures/output/fig03_southern_ocean_paleobathymetry_4panel.png
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


# ============================================================================
# Projection (LAEA, south-pole centred, 50 deg horizon -> cap edge at -40 S)
# ============================================================================
TIMES_MA              = [35, 25, 15, 0]
CENTRAL_LON           = 0.0
CENTRAL_LAT           = -90.0    # south pole
HORIZON_DEG           = 50.0     # cap edge at lat -40 S (90 - 50)
PROJ_WIDTH_CM_SINGLE  = 14
PROJ_WIDTH_CM_PANEL   = 14


# ============================================================================
# Present-day acronym labels (rendered only on the 0 Ma map).
# ----------------------------------------------------------------------------
# White bold text drawn at explicit (lon, lat) geographic coordinates.
# Each entry is (lon_deg, lat_deg, acronym).  Edit any line to nudge
# an individual label; comment out a line to drop a label; append a
# new tuple to add one.
# ============================================================================
PRESENT_DAY_ACRONYMS = [
    # (lon,  lat,  acronym -- full name)
    ( 146, -58,  "TG"),    # Tasman Gateway
    ( -59, -58,  "DP"),    # Drake Passage
    ( 146, -50,  "STR"),   # South Tasman Rise
    (  78, -50,  "KP"),    # Kerguelen Plateau
    ( 180, -52,  "CP"),    # Campbell Plateau
]
_ACRONYM_FONT = "16p,Helvetica-Bold,white"


def _add_present_day_acronyms(fig):
    """Place the white-text acronym labels (TG, DP, STR, KP, CP) on
    the present-day map only.  Called only when time_ma == 0."""
    for lon, lat, text in PRESENT_DAY_ACRONYMS:
        fig.text(x=lon, y=lat, text=text,
                 font=_ACRONYM_FONT, justify="MC", no_clip=True)


def plot_one_time(time_ma, gplot):
    print(f"  rendering {time_ma} Ma")
    fig = pygmt.Figure()
    pygmt.config(
        FONT_TITLE=FONT_TITLE,
        FONT_LABEL="24p,Helvetica",          # bigger colorbar LABEL
        FONT_ANNOT="16p,Helvetica",          # lat/lon labels at panel rim
                                              #   (colorbar ticks re-bumped to
                                              #   18p just before fig.colorbar
                                              #   below)
        COLOR_NAN=NAN_COLOR,
        MAP_FRAME_TYPE="fancy",               # fancy ring; the 90 deg divisions
                                              #   came from the "WSne" axis-letters
                                              #   in the frame= list -- WSne is now
                                              #   dropped (irrelevant on circular
                                              #   projections), and fancy follows
                                              #   the xa60g60 annotation interval.
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

    # LAEA projection (A) with south-pole centre and 50 deg horizon cap.
    # 60 deg lon / 30 deg lat graticule; no centred title (the age is
    # drawn as a boxed tag in the top-left after the basemap).
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

    # Present-day acronym labels (TG, DP, STR, KP, CP) only on 0 Ma.
    if int(time_ma) == 0:
        _add_present_day_acronyms(fig)

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
        OUTPUT_DIR, f"fig03_southern_ocean_paleobathymetry_{time_ma}Ma")
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
        MAP_FRAME_TYPE="fancy",               # fancy ring; the 90 deg divisions
                                              #   came from the "WSne" axis-letters
                                              #   in the frame= list -- WSne is now
                                              #   dropped (irrelevant on circular
                                              #   projections), and fancy follows
                                              #   the xa60g60 annotation interval.
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
        # Redraw graticule on top of the bathymetry.
        fig.basemap(frame=["xa60g60", "ya20g20"])
        # Present-day acronym labels only on the 0 Ma panel.
        if int(t) == 0:
            _add_present_day_acronyms(fig)
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
        OUTPUT_DIR, "fig03_southern_ocean_paleobathymetry_4panel")
    fig.savefig(base + ".png", dpi=300)
    fig.savefig(base + ".pdf")          # vector
    print(f"  wrote {base}.png")


def main():
    # Apply the user-tunable HILLSHADE_* constants from the top of this
    # file (and clear any cached gradients so the new settings take
    # effect on the very next pygmt.grdgradient call).
    _apply_hillshade_settings()

    # LAEA centred at (CENTRAL_LON, CENTRAL_LAT) -- pass to gplately so
    # coastline / continent polygons get split at the antimeridian and
    # don't wrap as fans across the visible cap.
    gplot = make_gplately_plot(central_meridian=CENTRAL_LON)
    for t in TIMES_MA:
        plot_one_time(t, gplot)
    plot_four_panel(gplot)


if __name__ == "__main__":
    main()
