#!/usr/bin/env python3
"""
paleobathy_render
=================

Shared rendering helpers used by both static paleobathymetry figures
(Fig 2 -- North-Atlantic-centred LAEA, Fig 3 -- Southern-Ocean
orthographic) and by the per-Myr video companion scripts.

This module is the SOLE source of truth for:

  * the plate-model + gplately driver setup        (make_gplately_plot)
  * gridded paleobathymetry rendering + hillshade  (draw_paleobathymetry)
  * the cached grdgradient pipeline                (gradient_for)
  * the colour palette, fonts and overlay style    (BATHY_*, COASTLINE_*, ...)

Consumer scripts import it as a normal Python module:

    import paleobathy_render as pbr

    gplot = pbr.make_gplately_plot(central_meridian=0)
    pbr.draw_paleobathymetry(fig, gplot, time_ma)

User-tunable knobs:
  * HILLSHADE_ON / HILLSHADE_AZIMUTH / HILLSHADE_NORMALIZE
    -- relief shading; edit at the top of each consumer script by
    mutating the module attribute (e.g. ``pbr.HILLSHADE_NORMALIZE = "e0.6"``)
    BEFORE the first call into draw_paleobathymetry / gradient_for.
  * DRAW_CONTINENTS / DRAW_COASTLINES / DRAW_TRENCHES_AND_TEETH
    -- overlay toggles.

NOTE: draw_paleobathymetry does NOT add the per-panel age label any
more.  The consumer script is expected to add its own fig.text() after
calling draw_paleobathymetry, so each figure can choose its own font,
size, position and box style without fighting a default.
"""
import logging
import os
import sys

# Silence fiona's flood of "Skipping field because of invalid value:
# key='feature_type', value=<pygplates.pygplates.FeatureType ...>"
# warnings (fired by gplately's get_* accessors).  fiona drops only
# that one unserialisable column; the geometries + plate IDs survive.
logging.getLogger("fiona").setLevel(logging.ERROR)

import pygmt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (PALEO_BATHY_FMT, OUTPUT_DIR,
                    PLATE_MODEL_NAME, PLATE_MODEL_ANCHOR_PLATE,
                    PALEO_BATHY_SOURCE)


# ----------------------------------------------------------------------------
# Cache directories.  Source-aware so the Z22 mantle-frame and Z22
# paleomagnetic-frame variants don't clobber each other's downloaded
# plate model or gradient files.
# ----------------------------------------------------------------------------
PLATE_MODEL_CACHE = os.path.join(
    OUTPUT_DIR, f"plate-model-repo_{PALEO_BATHY_SOURCE}")
os.makedirs(PLATE_MODEL_CACHE, exist_ok=True)

GRAD_DIR = os.path.join(
    OUTPUT_DIR, f"paleobathy_gradients_{PALEO_BATHY_SOURCE}")
os.makedirs(GRAD_DIR, exist_ok=True)


# ----------------------------------------------------------------------------
# Colour palette + cpt range.  GMT's "bathy" is the classic blue ramp;
# VMIN=-5000 / VMAX=0 covers the bathymetric half of typical seafloor.
# ----------------------------------------------------------------------------
VMIN = -5000.0
VMAX = 0.0
BATHY_CMAP = "bathy"
BATHY_CPT_SERIES = [VMIN, VMAX, 100]


# ----------------------------------------------------------------------------
# Hillshade.  DEFAULTS only -- consumers are encouraged to override via
# `pbr.HILLSHADE_NORMALIZE = "..."` at the top of their script.
#
# Mechanics: gradient_for(time_ma) calls
#     pygmt.grdgradient(grid=..., azimuth=HILLSHADE_AZIMUTH,
#                       normalize=HILLSHADE_NORMALIZE, outgrid=<cached.nc>)
# and caches the resulting per-age intensity grid under GRAD_DIR.
# draw_paleobathymetry then passes that .nc as `shading=...` to
# fig.grdimage.  Changing HILLSHADE_NORMALIZE without wiping the cache
# is a silent no-op for ages whose gradient is already on disk.
# ----------------------------------------------------------------------------
HILLSHADE_ON         = True
HILLSHADE_AZIMUTH    = "315"     # deg clockwise from north (315 = NW)
HILLSHADE_NORMALIZE  = "e0.6"    # "t<sigma>" tangent or "e<sigma>" exp.
                                 # Smaller sigma -> sharper relief.
                                 # "t1"   = very gentle tangent
                                 # "t0.5" = moderate tangent
                                 # "e0.6" = sharper exponential
HILLSHADE_INTENSITY  = None      # None = use cached gradient unchanged;
                                 # float in (0, 1] = dim further at
                                 # grdimage time via "+a<val>" modifier.


# ----------------------------------------------------------------------------
# Overlay toggles + styles.
# ----------------------------------------------------------------------------
DRAW_CONTINENTS         = False  # continent polygons look like COB lines;
                                 #   off by default
DRAW_COASTLINES         = True
DRAW_ALL_BOUNDARIES     = False  # OFF by default (user 2026-05-31): the
                                 #   thin grey underlay of the full
                                 #   topological-section network read as
                                 #   thin black lines on the published
                                 #   plots.  Subduction zones still show
                                 #   via DRAW_TRENCHES_AND_TEETH below.
DRAW_TRENCHES_AND_TEETH = True

NAN_COLOR        = "200/200/200"
CONTINENT_FILL   = "#c0c0c0"
# Suppress the coastline-polygon outline.  GMT 6.5 on this install
# silently ignores both "0p,white@100" and "0p" (the parsed pen
# falls back to the 0.25p black default), so we make the pen match
# the fill colour -- any line GMT chooses to draw is then visually
# indistinguishable from the silver coastline interior.
COASTLINE_PEN    = f"default,{CONTINENT_FILL}"
ALL_BOUNDARY_PEN = "0.35p,#666666"
TRENCH_PEN       = "0.8p,red"
TEETH_GAP        = "0.5c"
TEETH_HEIGHT     = "0.18c"
TEETH_FILL       = "red"
TEETH_PEN        = "0.4p,red"


# ----------------------------------------------------------------------------
# Default fonts used by the colourbar in single-panel layouts.  Consumer
# scripts typically override these via pygmt.config(FONT_LABEL=...,
# FONT_ANNOT=...) inside their plot_one_time() / plot_four_panel().
# ----------------------------------------------------------------------------
FONT_TITLE = "22p,Helvetica-Bold"
FONT_LABEL = "16p,Helvetica"
FONT_ANNOT = "14p,Helvetica"


# ----------------------------------------------------------------------------
# Plate-model + gplately driver setup.
# ----------------------------------------------------------------------------
def make_gplately_plot(central_meridian=0.0):
    """Build the gplately.PlotTopologies driver for the configured plate
    model.

    `central_meridian` -- the central meridian of the map projection the
    rendered polygons will be plotted on.  gplately uses this to split
    coastline / continent polygons at the antipodal meridian
    (`central_meridian +/- 180`) so they render correctly in projections
    that cannot natively handle a polygon crossing the antimeridian.
    Without this, Antarctic polygons in particular wrap as long fans
    across the visible hemisphere at certain reconstruction times
    (e.g. 38 Ma in the Southern Ocean orthographic view).

    Version-tolerant: the constructor signature for `central_meridian`
    was added to PlotTopologies in newer gplately releases.  On older
    versions the kwarg raises TypeError; we fall back to setting the
    attribute and forwarding `central_meridian` per-call through
    `get_coastlines() / get_continents()` (see `_get_polys`).  The
    value is stashed on the returned gplot as
    `_pyback_central_meridian`.
    """
    import gplately
    from plate_model_manager import PlateModelManager

    pm = PlateModelManager()
    plate_model = pm.get_model(PLATE_MODEL_NAME, data_dir=PLATE_MODEL_CACHE)
    rotation_model = plate_model.get_rotation_model()
    topology_features = plate_model.get_topologies()
    static_polygons = plate_model.get_static_polygons()
    coastlines = plate_model.get_coastlines()
    continents = plate_model.get_continental_polygons()

    plate_reconstruction = gplately.PlateReconstruction(
        rotation_model, topology_features, static_polygons,
        anchor_plate_id=PLATE_MODEL_ANCHOR_PLATE,
    )
    try:
        gplot = gplately.PlotTopologies(
            plate_reconstruction,
            coastlines=coastlines,
            continents=continents,
            central_meridian=central_meridian,
        )
    except TypeError:
        gplot = gplately.PlotTopologies(
            plate_reconstruction,
            coastlines=coastlines,
            continents=continents,
        )
        try:
            gplot.central_meridian = central_meridian
        except Exception:
            pass
    gplot._pyback_central_meridian = central_meridian
    return gplot


def _get_polys(gplot_method, time_ma, central_meridian):
    """Call a gplately `get_*` accessor and (when supported) forward
    `central_meridian` so the returned GeoDataFrame has polygons
    pre-split at the antimeridian.  Falls back silently to the
    no-kwarg call on older gplately versions.
    """
    if central_meridian is not None:
        try:
            return gplot_method(central_meridian=central_meridian)
        except TypeError:
            pass
    return gplot_method()


# ----------------------------------------------------------------------------
# Trenches + subduction teeth.
# ----------------------------------------------------------------------------
def _split_by_polarity(gdf):
    """Return (left_gdf, right_gdf) by sniffing whichever polarity column
    gplately's current release decided to use."""
    if gdf is None:
        return None, None
    for col in ("polarity", "subduction_polarity", "POLARITY",
                "SUBDUCTION_POLARITY", "tooth_side", "side",
                "subduction_side"):
        if col in gdf.columns:
            vals = gdf[col].astype(str).str.lower().str.strip()
            return gdf[vals.isin(("left", "l"))], gdf[vals.isin(("right", "r"))]
    return gdf, gdf.iloc[0:0]


def _coerce_left_right(result):
    """Normalise gplately's get_subduction_direction() return value to
    (left_gdf, right_gdf) across the (tuple | single-GeoDataFrame)
    variation seen across gplately releases."""
    if isinstance(result, tuple) and len(result) == 2:
        return result[0], result[1]
    return _split_by_polarity(result)


_ALL_BOUNDARY_GETTERS = (
    "get_all_topological_sections",
    "get_all_topologies",
    "get_topological_plate_boundaries",
)


def _plot_all_plate_boundaries(fig, gplot, time_ma):
    """Plot every topological boundary segment as one thin grey pen so the
    plate-boundary network is visibly closed regardless of how each
    individual segment is typed in the plate model."""
    for attr in _ALL_BOUNDARY_GETTERS:
        getter = getattr(gplot, attr, None)
        if getter is None:
            continue
        try:
            gdf = getter()
        except Exception as exc:
            print(f"  ({time_ma} Ma) {attr}() raised: {exc}")
            continue
        if gdf is None or len(gdf) == 0:
            continue
        try:
            fig.plot(data=gdf, pen=ALL_BOUNDARY_PEN)
            return
        except Exception as exc:
            print(f"  ({time_ma} Ma) {attr} plot failed: {exc}")

    # Fallback: union whatever individual boundary-type getters exist.
    for attr in ("get_ridges_and_transforms", "get_ridges",
                 "get_transforms", "get_trenches", "get_misc_boundaries"):
        getter = getattr(gplot, attr, None)
        if getter is None:
            continue
        try:
            gdf = getter()
        except Exception:
            continue
        if gdf is None or len(gdf) == 0:
            continue
        try:
            fig.plot(data=gdf, pen=ALL_BOUNDARY_PEN)
        except Exception as exc:
            print(f"  ({time_ma} Ma) {attr} (fallback) skipped: {exc}")


def _plot_trenches_with_teeth(fig, gplot, time_ma):
    """Plot trenches as red polylines plus subduction teeth on the
    overriding-plate side."""
    trench_line_gdf = None
    try:
        trench_line_gdf = gplot.get_trenches()
    except Exception as exc:
        print(f"  ({time_ma} Ma) get_trenches() failed: {exc}")
    if trench_line_gdf is not None and len(trench_line_gdf) > 0:
        try:
            fig.plot(data=trench_line_gdf, pen=TRENCH_PEN)
        except Exception as exc:
            print(f"  ({time_ma} Ma) trench line skipped: {exc}")

    try:
        result = gplot.get_subduction_direction()
    except Exception as exc:
        print(f"  ({time_ma} Ma) get_subduction_direction() failed: {exc}")
        return
    if result is None:
        print(f"  ({time_ma} Ma) no subduction-direction data")
        return

    left, right = _coerce_left_right(result)

    if left is not None and len(left) > 0:
        try:
            fig.plot(
                data=left,
                style=f"f{TEETH_GAP}/{TEETH_HEIGHT}+l+t",
                fill=TEETH_FILL, pen=TEETH_PEN,
            )
        except Exception as exc:
            print(f"  ({time_ma} Ma) left teeth skipped: {exc}")

    if right is not None and len(right) > 0:
        try:
            fig.plot(
                data=right,
                style=f"f{TEETH_GAP}/{TEETH_HEIGHT}+r+t",
                fill=TEETH_FILL, pen=TEETH_PEN,
            )
        except Exception as exc:
            print(f"  ({time_ma} Ma) right teeth skipped: {exc}")


# ----------------------------------------------------------------------------
# Hillshade gradient (cached per age).
# ----------------------------------------------------------------------------
def gradient_for(time_ma):
    """Return (grid_path, gradient_path) for `time_ma`, computing the
    gradient via pygmt.grdgradient if it isn't already cached."""
    grid_path = PALEO_BATHY_FMT.format(time=float(time_ma))
    if not os.path.exists(grid_path):
        raise FileNotFoundError(
            f"paleobathy grid for {time_ma} Ma not found at {grid_path}")
    grad_path = os.path.join(GRAD_DIR, f"grad_{int(time_ma)}Ma.nc")
    if not os.path.exists(grad_path):
        pygmt.grdgradient(
            grid=grid_path,
            azimuth=HILLSHADE_AZIMUTH,
            normalize=HILLSHADE_NORMALIZE,
            outgrid=grad_path,
        )
    return grid_path, grad_path


def wipe_gradient_cache(times_ma):
    """Delete cached gradient .nc files for the listed ages.  Call this
    from a consumer's main() whenever HILLSHADE_* settings have just been
    changed -- otherwise gradient_for() reuses any existing cached file
    and the new settings appear to have no effect."""
    if not os.path.isdir(GRAD_DIR):
        return 0
    n = 0
    for t in times_ma:
        p = os.path.join(GRAD_DIR, f"grad_{int(t)}Ma.nc")
        if os.path.exists(p):
            os.remove(p)
            n += 1
    if n:
        print(f"  cleared {n} cached gradient(s) from {GRAD_DIR} "
              f"(HILLSHADE_NORMALIZE={HILLSHADE_NORMALIZE!r})")
    return n


# ----------------------------------------------------------------------------
# Bathymetry + overlay rendering for a single time slice.
#
# No basemap call here -- the caller has already set region + projection
# (either via fig.basemap or via the subplot context with projection=).
# No age-label fig.text either -- the caller adds its own so each figure
# can pick its own font / size / position.
# ----------------------------------------------------------------------------
def draw_paleobathymetry(fig, gplot, time_ma):
    """Draw the paleobathymetry grid + (optional) hillshade + (optional)
    plate-model overlay (coastlines, full boundary network, trenches +
    teeth) for `time_ma` into the current pyGMT panel."""
    grid_path = PALEO_BATHY_FMT.format(time=float(time_ma))
    if not os.path.exists(grid_path):
        raise FileNotFoundError(
            f"paleobathy grid for {time_ma} Ma not found at {grid_path}")

    grdimage_kwargs = dict(grid=grid_path, cmap=True)
    if HILLSHADE_ON:
        _, grad_path = gradient_for(time_ma)
        if HILLSHADE_INTENSITY is None:
            grdimage_kwargs["shading"] = grad_path
        else:
            grdimage_kwargs["shading"] = f"{grad_path}+a{HILLSHADE_INTENSITY}"
    fig.grdimage(**grdimage_kwargs)

    gplot.time = float(time_ma)
    cm = getattr(gplot, "_pyback_central_meridian", None)

    if DRAW_CONTINENTS:
        try:
            gdf = _get_polys(gplot.get_continents, time_ma, cm)
            if gdf is not None and len(gdf) > 0:
                fig.plot(data=gdf, fill=CONTINENT_FILL, pen="0p")
        except Exception as exc:
            print(f"  ({time_ma} Ma) continents skipped: {exc}")

    if DRAW_COASTLINES:
        try:
            gdf = _get_polys(gplot.get_coastlines, time_ma, cm)
            if gdf is not None and len(gdf) > 0:
                if COASTLINE_PEN is None:
                    fig.plot(data=gdf, fill=CONTINENT_FILL)
                else:
                    fig.plot(data=gdf, fill=CONTINENT_FILL,
                             pen=COASTLINE_PEN)
        except Exception as exc:
            print(f"  ({time_ma} Ma) coastlines skipped: {exc}")

    if DRAW_ALL_BOUNDARIES:
        _plot_all_plate_boundaries(fig, gplot, time_ma)
    if DRAW_TRENCHES_AND_TEETH:
        _plot_trenches_with_teeth(fig, gplot, time_ma)
