#!/usr/bin/env python3
"""
Figure 10 - Combined 3 x 4 matrix of NW Shelf tectonic-subsidence-rate
maps at four well-separated times where the SL+DT correction has the
largest spatial impact (max |C - A|).

Columns (configurations):
    A: no sea level, no dynamic topography                    cpt: dem1
    C: Haq2024 long-term SL + D10_gmcm9 DT (Braz et al. 2021) cpt: dem1
    D = C - A                                                 cpt: vik

Rows: four times in 0-150 Ma, picked greedily from the per-time
spatial-mean of |D| with `DT_TIME_SEP` Myr minimum separation so the
four times spread across the record rather than clustering around
one peak.  Selected times + the full per-time stats are written to
a sidecar text file for the paper's methods section.

Rate-map pipeline per panel:
    blockmedian (0.05 deg) -> surface (tension 0.5) ->
    mask by the concave hull (alpha shape) of the 109 well points,
    dilated by HULL_BUFFER_KM, with KEEP_LARGEST_POLYGON_ONLY
    enforcing one connected polygon.  Well points are over-plotted
    as small white-filled circles on every panel.

Inputs (produced by 07a_backstrip_all_nwshelf.py):
    figures/output/nwshelf_subsidence/A_no_sl_no_dt_rate/rate_<t>.nc
    figures/output/nwshelf_subsidence/C_sl_and_dt_rate/rate_<t>.nc
    figures/output/nwshelf_subsidence/D_rate_difference/rate_<t>.nc
    figures/output/nwshelf_subsidence/all_wells_locations.txt
    figures/output/nwshelf_subsidence/cutoff_time_Ma.txt

Outputs:
    figures/output/fig10_nwshelf_rate_maps_combined.png
    figures/output/fig10_nwshelf_rate_maps_combined.pdf
    figures/output/nwshelf_subsidence/rate_delta_stats.csv
    figures/output/nwshelf_subsidence/rate_delta_picked_times.txt
"""
import os
import sys
import time

import numpy as np
import pygmt
import xarray as xr
from scipy.spatial import Delaunay
import shapely
from shapely.geometry import Polygon
from shapely.ops import unary_union

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import OUTPUT_DIR, DYNAMIC_TOPOGRAPHY_MODEL, MAX_ANALYSIS_AGE

# ----------------------------------------------------------------------------
# Tunables
# ----------------------------------------------------------------------------
OUT_BASE = os.path.join(OUTPUT_DIR, "nwshelf_subsidence")
REGION = [113, 132, -22, -8]

A_DIR = os.path.join(OUT_BASE, "A_no_sl_no_dt_rate")
C_DIR = os.path.join(OUT_BASE, "C_sl_and_dt_rate")
D_DIR = os.path.join(OUT_BASE, "D_rate_difference")

# CPT settings.  Rate columns A and C use GMT's `dem1` topo palette
# (warm-cool ramp); the difference column D uses Crameri's `vik`
# diverging palette centred on zero so positive (under-estimated
# subsidence without SL/DT) and negative (over-estimated) corrections
# are immediately distinguishable.  The two ranges are independent so
# the small-amplitude D field is not washed out by the larger A/C
# values.  Out-of-range A/C values are clamped to the end colours and
# indicated by the +e side-arrow triangles on the colourbar.
RATE_CMAP = "dem4"                   # for A and C
RATE_SERIES = (0, 50, 2.5)           # m/Myr; clamped + arrowed out of range
DELTA_CMAP = "vik"                   # for delta (diverging blue-white-red)
DELTA_SERIES = (-50, 50, 5)          # m/Myr, symmetric around zero

# How to pick the 4 representative times.
N_TIMES = 4
DT_TIME_SEP = 25.0                   # Myr, minimum spacing between picks

# Grid layout.  CRITICAL: with projection "M?" each Mercator panel's
# height is determined by the lat range and the panel WIDTH, not by
# the cell height -- so any cell taller than the projection becomes
# whitespace INSIDE the cell, which reads as a gap between rows.
# For REGION=[113,132,-22,-8] and a 6.77c column width, each Mercator
# panel is ~5.15c tall; four rows + bottom-axis labels need ~22c
# total figure height.
FIG_FIGSIZE = ("21c", "22c")         # 3 columns x 4 rows
SUBPLOT_MARGINS = ["0.35c", "0c"]    # [horizontal, vertical] -- 0c
                                     # vertical so the panel cells butt
                                     # right up against each other; the
                                     # bottom-row sharex hides any
                                     # tick-label overlap.
PANEL_TITLE_FONT = "16p,Helvetica-Bold"

# Continuous-gridding parameters.  Each input rate field is re-gridded
# via blockmedian -> surface to fill 07a's per-well gaps into a
# continuous field, then masked against the concave hull of the well
# locations (see the alpha-shape parameters below).
REGRID_SPACING = 0.05                # degrees (~5.5 km in lon)
REGRID_TENSION = 0.50                # GMT surface tension.  T=0 is the
                                     # pure minimum-curvature (natural)
                                     # spline; T=1 is the pure harmonic
                                     # (Laplacian) solution.  T=0.5 is
                                     # a balanced compromise -- enough
                                     # tension to suppress overshoot
                                     # wiggles near sparse data
                                     # without flattening the field
                                     # into the harmonic limit.
# Concave-hull (alpha-shape) parameters.  We mask via the concave
# envelope of the 109 well locations, NOT a union of circles, so the
# resulting field does not have the visible blob/scallop boundary that
# any radius-based mask produces.  The hull is built by:
#   1. Delaunay-triangulating the well points in local equirectangular
#      km coordinates.
#   2. Keeping only triangles where every edge is shorter than
#      `ALPHA_MAX_EDGE_KM`.  Smaller threshold = tighter concave fit;
#      larger threshold = closer to the convex hull.
#   3. `shapely.unary_union` of the kept triangles -> a single polygon
#      that follows the actual outline of the well cluster.
ALPHA_MAX_EDGE_KM = 400              # km, edges longer than this drop
                                     # out of the alpha shape.  Tuned so
                                     # the 109 NW Shelf wells produce a
                                     # single connected polygon (smaller
                                     # values leave isolated outlier
                                     # wells as their own micro-hull).
HULL_BUFFER_KM = 80                  # km of dilation added around the
                                     # concave hull so the smooth
                                     # surface reaches well beyond the
                                     # outermost wells.  Also helps fuse
                                     # any residual disconnected pieces
                                     # into a single connected polygon.
KEEP_LARGEST_POLYGON_ONLY = True     # if the alpha shape + buffer still
                                     # produces a MultiPolygon, keep
                                     # only the largest piece.  Set to
                                     # False to render every piece.

# Per-time metric for "how big is the SL+DT correction here?".
# We use spatial-mean(|D|) which gives equal weight to all pixels in
# the masked NW Shelf area.  Alternatives: spatial std(D), or quantile.
METRIC_NAME = "mean_abs_delta"


# ----------------------------------------------------------------------------
# Stats & time selection
# ----------------------------------------------------------------------------
def collect_per_time_stats():
    """Walk D_DIR and build a per-time stats table.

    Returns
    -------
    times : np.ndarray of int, sorted ascending
    stats : dict mapping name -> np.ndarray aligned with `times`.
            Keys: mean_A, std_A, mean_C, std_C, mean_D, std_D,
            mean_abs_delta, n_pixels
    """
    times = []
    for fn in sorted(os.listdir(D_DIR)):
        if not fn.startswith("rate_") or not fn.endswith(".nc"):
            continue
        try:
            t = int(fn[len("rate_"):-len(".nc")])
        except ValueError:
            continue
        if 0 <= t <= int(MAX_ANALYSIS_AGE):
            times.append(t)
    times = np.array(sorted(times), dtype=int)
    print(f"Found {len(times)} delta-rate grids in {D_DIR}")

    out = {k: np.full(len(times), np.nan, dtype=float)
           for k in ("mean_A", "std_A", "mean_C", "std_C",
                     "mean_D", "std_D", "mean_abs_delta", "n_pixels")}
    for i, t in enumerate(times):
        for name, sub in (("A", A_DIR), ("C", C_DIR), ("D", D_DIR)):
            path = os.path.join(sub, f"rate_{int(t)}.nc")
            if not os.path.exists(path):
                continue
            arr = xr.open_dataset(path).z.values
            finite = np.isfinite(arr)
            if not finite.any():
                continue
            arr_v = arr[finite]
            out[f"mean_{name}"][i] = float(np.mean(arr_v))
            out[f"std_{name}"][i] = float(np.std(arr_v))
            if name == "D":
                out["mean_abs_delta"][i] = float(np.mean(np.abs(arr_v)))
                out["n_pixels"][i] = int(finite.sum())
    return times, out


def write_stats_csv(times, stats):
    path = os.path.join(OUT_BASE, "rate_delta_stats.csv")
    with open(path, "w") as f:
        f.write(
            "# Per-time NW-Shelf-wide tectonic-subsidence-rate statistics.\n"
            "# A = no SL, no DT.  C = Haq2024 SL + "
            f"{DYNAMIC_TOPOGRAPHY_MODEL} DT.  D = C - A.  All in m/Myr.\n"
            "time_Ma,mean_A,std_A,mean_C,std_C,mean_D,std_D,"
            "mean_abs_delta,n_pixels\n"
        )
        for i, t in enumerate(times):
            row = [str(int(t))]
            for k in ("mean_A", "std_A", "mean_C", "std_C",
                      "mean_D", "std_D", "mean_abs_delta", "n_pixels"):
                v = stats[k][i]
                row.append("" if np.isnan(v) else f"{v:.3f}")
            f.write(",".join(row) + "\n")
    print(f"  wrote {path}")
    return path


def pick_max_delta_times(times, metric, n, min_sep):
    """Greedy non-maximum suppression to pick `n` times where `metric`
    is large but separated by at least `min_sep` Myr.

    Returns the picks sorted oldest-to-youngest (since rows of the
    figure go top-down = oldest at top).
    """
    available = np.ones(len(times), dtype=bool)
    available &= np.isfinite(metric)
    picks = []
    for _ in range(n):
        if not available.any():
            break
        cand = metric.copy()
        cand[~available] = -np.inf
        idx = int(np.argmax(cand))
        picks.append(idx)
        t0 = times[idx]
        # Mask out everything within min_sep
        for j, t in enumerate(times):
            if abs(t - t0) < min_sep:
                available[j] = False
    picks_t = sorted([int(times[i]) for i in picks], reverse=True)
    return picks_t


# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------
def _grid_is_drawable(path):
    """Quick check: file exists AND at least one finite cell."""
    if not os.path.exists(path):
        return False
    try:
        arr = xr.open_dataset(path).z.values
        return bool(np.isfinite(arr).any())
    except Exception as exc:
        print(f"    ! could not read {path}: {exc}")
        return False


# RGB endpoints of the diverging master CPTs we support for hinged
# (asymmetric) ranges.  Mapping: cmap_name -> (low_end, midpoint,
# high_end).  All three are RGB triples interpreted as the colours of
# the master CPT at z=0, z=0.5 and z=1 respectively.  When `reverse=True`
# is requested we swap low_end and high_end.
#
# We hand-code these (instead of round-tripping through `makecpt -G`)
# because the truncate-then-concatenate path silently produced a
# white->red gradient over the full range in practice -- the half-CPT
# files appear to be written empty or have a format that doesn't
# concatenate cleanly into GMT 6.x.  Writing the CPT bytes directly
# from Python is robust and gives us bit-exact control.
_DIVERGING_CPT_ENDPOINTS = {
    "polar": ((0, 0, 255), (255, 255, 255), (255, 0, 0)),
}

# Cmaps with an INTRINSIC HARD HINGE at z=0 baked into the master CPT
# (e.g. topo's sea-level break, earth, etoposl).  For these, GMT 6.x
# auto-anchors the hinge to z=0 when -T spans across 0 -- the two
# halves of the master CPT are stretched independently to [lo, 0] and
# [0, hi].  No custom Python construction is needed; plain makecpt
# with the requested series is correct.
_INTRINSIC_HARD_HINGE_CMAPS = {"topo", "earth", "etoposl", "geo"}


def _save_hinged_diverging_cpt(cmap, lo, hi, hinge, step, out_path,
                               reverse=False):
    """Write a custom hinged-diverging CPT directly to ``out_path``.

    Interpolates linearly in RGB between the master CPT's documented
    endpoints (see ``_DIVERGING_CPT_ENDPOINTS``) so that:

    * z in [lo, hinge] runs low_end -> midpoint
    * z in [hinge, hi] runs midpoint -> high_end

    with the colour transition pinned at ``z = hinge`` regardless of
    how asymmetric ``[lo, hi]`` is around ``hinge``.

    Step size is shared across both halves.  Both halves may end with
    a short final cell if (hinge - lo) or (hi - hinge) is not an
    integer multiple of ``step``.
    """
    try:
        c_low_end, c_mid, c_high_end = _DIVERGING_CPT_ENDPOINTS[cmap]
    except KeyError as exc:
        raise NotImplementedError(
            f"_save_hinged_diverging_cpt: no hard-coded endpoints for "
            f"cmap '{cmap}'.  Add an entry to "
            f"_DIVERGING_CPT_ENDPOINTS or call _save_cpt without "
            f"a hinge."
        ) from exc
    if reverse:
        c_low_end, c_high_end = c_high_end, c_low_end

    def _interp(c0, c1, t):
        return tuple(int(round(c0[i] + t * (c1[i] - c0[i])))
                     for i in range(3))

    def _fmt(c):
        return f"{c[0]}/{c[1]}/{c[2]}"

    def _emit_half(fh, z_start, z_end, c_start, c_end):
        """Emit step-spaced data rows from z_start to z_end (z_start
        is inclusive, z_end is the last z that should be present)."""
        span = z_end - z_start
        z = z_start
        eps = step * 1e-6
        while z < z_end - eps:
            z_next = min(z + step, z_end)
            t_lo = (z - z_start) / span
            t_hi = (z_next - z_start) / span
            fh.write(
                f"{z:g} {_fmt(_interp(c_start, c_end, t_lo))} "
                f"{z_next:g} {_fmt(_interp(c_start, c_end, t_hi))}\n"
            )
            z = z_next

    with open(out_path, "w") as fh:
        fh.write(f"# Hinged diverging CPT, cmap={cmap}, "
                 f"range=[{lo:g}, {hi:g}], hinge={hinge:g}, "
                 f"step={step:g}, reverse={reverse}\n")
        fh.write("# COLOR_MODEL = RGB\n")
        _emit_half(fh, lo, hinge, c_low_end, c_mid)
        _emit_half(fh, hinge, hi, c_mid, c_high_end)
        fh.write(f"B {_fmt(c_low_end)}\n")
        fh.write(f"F {_fmt(c_high_end)}\n")
        fh.write("N 200/200/200\n")
    return out_path


def _save_cpt(cmap, series, out_path, reverse=False, background=False,
              hinge=None):
    """Build a CPT once and persist it to disk so the rendering loop can
    pass a *filename* to grdimage and avoid issuing makecpt inside the
    subplot context (which has caused hangs on macOS Apple-Silicon).

    Pass ``reverse=True`` to flip the colour progression of the master
    CPT (pygmt's documented `+i`-equivalent flag).

    Pass ``background=True`` (GMT `-D` flag) to clamp the back- and
    fore-ground colours of the resulting CPT to the bottom / top
    colour of the master CPT.  Use this for "clamp" semantics on
    one-sided ranges (e.g. rate 0..50 m/Myr): values above 50 then
    render as the 50-m/Myr colour instead of GMT's default white.

    Pass ``hinge=<z>`` to anchor the midpoint of a diverging master CPT
    at the data value ``z`` (typically 0).  Routes to
    ``_save_hinged_diverging_cpt`` which writes the CPT directly in
    Python, bypassing ``makecpt`` entirely.  Only the cmaps listed in
    ``_DIVERGING_CPT_ENDPOINTS`` are supported in this mode.
    """
    lo, hi, step = float(series[0]), float(series[1]), float(series[2])
    if hinge is not None and lo < float(hinge) < hi:
        if cmap in _INTRINSIC_HARD_HINGE_CMAPS:
            # Master CPT carries its own hinge tag; trust GMT to
            # anchor it.  Fall through to plain makecpt below.
            pass
        else:
            return _save_hinged_diverging_cpt(
                cmap, lo, hi, float(hinge), step, out_path,
                reverse=reverse,
            )
    pygmt.makecpt(
        cmap=cmap, series=list(series), continuous=True,
        reverse=reverse, background=background, output=out_path,
    )
    return out_path


def _lon_lat_arrays(ds_or_da):
    """Return ``(lon_values, lat_values)`` for a NetCDF that uses either
    ``lon``/``lat`` or ``x``/``y`` as coord names.

    pyGMT's ``Figure.surface`` and other tabular-to-grid tools write
    the output NetCDF with ``x``/``y`` coordinate variables (the GMT
    convention).  The bundled paleobathymetry grids, by contrast, use
    ``lon``/``lat``.  We probe both, then fall back to the last two
    dims of the data variable.

    Accepts either an xarray.Dataset (with a 'z' variable) or an
    xarray.DataArray directly.
    """
    if isinstance(ds_or_da, xr.DataArray):
        da = ds_or_da
        coords = da.coords
        dims = da.dims
    else:
        da = ds_or_da["z"]
        coords = ds_or_da.coords
        dims = da.dims
    if "lon" in coords:
        lon = coords["lon"].values
    elif "x" in coords:
        lon = coords["x"].values
    else:
        lon = coords[dims[-1]].values
    if "lat" in coords:
        lat = coords["lat"].values
    elif "y" in coords:
        lat = coords["y"].values
    else:
        lat = coords[dims[-2]].values
    return lon, lat


def _local_km_projector(region):
    """Return a (lon, lat) -> (x_km, y_km) converter for the region.

    Uses an equirectangular projection anchored on the region's
    mid-latitude.  Adequate for NW-Shelf-scale geometry: distortion at
    the corners is well under 1 percent compared with the haversine
    distance, far below our 30-60 km tolerance for mask boundaries.
    """
    lat0 = 0.5 * (region[2] + region[3])
    coslat = float(np.cos(np.radians(lat0)))

    def project(lon, lat):
        lon = np.asarray(lon, dtype=float)
        lat = np.asarray(lat, dtype=float)
        return lon * 111.0 * coslat, lat * 111.0

    return project, coslat


def _alpha_shape_polygon(well_lon, well_lat, region,
                         max_edge_km=ALPHA_MAX_EDGE_KM,
                         buffer_km=HULL_BUFFER_KM):
    """Concave hull of the well cluster -- single polygon, no circles.

    Algorithm
    ---------
    1. Project the wells into local km coordinates.
    2. Delaunay-triangulate.
    3. Keep only triangles whose longest edge is shorter than
       `max_edge_km` (the alpha-shape "longest-edge" filter).  Bigger
       triangles bridge gaps in the well distribution and would
       otherwise add the convex-hull tail back in.
    4. `shapely.unary_union` of the kept triangles -> the concave hull
       as a (possibly multi-) polygon.
    5. Optional dilation by `buffer_km` so the smooth surface can run
       a bit beyond the outermost wells.
    6. Project back to (lon, lat) for the masking step.

    Returns
    -------
    shapely.geometry.(Multi)Polygon in (lon, lat) coordinates.
    """
    project, coslat = _local_km_projector(region)
    wx, wy = project(well_lon, well_lat)
    pts_km = np.column_stack([wx, wy])
    tri = Delaunay(pts_km)
    triangles = pts_km[tri.simplices]            # (n_tri, 3, 2)
    # Longest edge of each triangle.
    a = np.linalg.norm(triangles[:, 0] - triangles[:, 1], axis=1)
    b = np.linalg.norm(triangles[:, 1] - triangles[:, 2], axis=1)
    c = np.linalg.norm(triangles[:, 2] - triangles[:, 0], axis=1)
    longest = np.maximum.reduce([a, b, c])
    keep = longest < max_edge_km
    kept_polys = [Polygon(triangles[i]) for i in np.flatnonzero(keep)]
    if not kept_polys:
        # Fall back to convex hull rather than crash; user can lift
        # max_edge_km to recover the alpha-shape behaviour.
        print(f"    ! alpha shape produced 0 triangles at "
              f"max_edge={max_edge_km} km -- falling back to convex hull")
        from scipy.spatial import ConvexHull
        hull = ConvexHull(pts_km)
        hull_km = Polygon(pts_km[hull.vertices])
    else:
        hull_km = unary_union(kept_polys)
    # Dilate by buffer_km so the surface reaches well beyond the
    # outermost wells.  Buffer is in the same km units as the hull.
    # A large buffer also fuses any residual disconnected sub-hulls
    # (small clusters of outlier wells) into one polygon.
    if buffer_km > 0:
        hull_km = hull_km.buffer(buffer_km)
    # If the result is still a MultiPolygon (rare with a big buffer,
    # but possible when an outlier well sits more than 2 * buffer_km
    # from the main cluster), keep only the largest piece -- the user
    # explicitly wants a single hull.
    if (KEEP_LARGEST_POLYGON_ONLY and
            getattr(hull_km, "geom_type", "") == "MultiPolygon"):
        biggest = max(hull_km.geoms, key=lambda p: p.area)
        n_dropped = len(hull_km.geoms) - 1
        print(f"    alpha-shape after buffer is a MultiPolygon "
              f"({len(hull_km.geoms)} pieces); keeping the largest "
              f"({biggest.area:.0f} km^2) and dropping {n_dropped} smaller "
              f"piece(s)")
        hull_km = biggest
    # Project the hull boundary back to (lon, lat).
    def unproject(x_km, y_km):
        return x_km / (111.0 * coslat), y_km / 111.0

    def transform_polygon(poly):
        if poly.is_empty:
            return poly
        # poly may be Polygon or MultiPolygon
        if poly.geom_type == "MultiPolygon":
            return shapely.geometry.MultiPolygon(
                [transform_polygon(p) for p in poly.geoms]
            )
        ext_x, ext_y = poly.exterior.coords.xy
        ext = np.column_stack(unproject(np.asarray(ext_x), np.asarray(ext_y)))
        interiors = []
        for ring in poly.interiors:
            ix, iy = ring.coords.xy
            interiors.append(
                np.column_stack(unproject(np.asarray(ix), np.asarray(iy)))
            )
        return Polygon(ext, holes=interiors)

    return transform_polygon(hull_km)


def _concave_hull_mask(grid_da, well_lon, well_lat, region):
    """Build a 1.0/NaN mask aligned with ``grid_da`` from the
    concave hull (alpha shape) of the well points -- a single
    connected polygon, not a union of circles.
    """
    poly = _alpha_shape_polygon(
        well_lon, well_lat, region,
        max_edge_km=ALPHA_MAX_EDGE_KM,
        buffer_km=HULL_BUFFER_KM,
    )
    # Vectorised inside-polygon test via shapely.contains_xy (shapely
    # >=2.0) -- avoids building one Point per cell.
    lon, lat = _lon_lat_arrays(grid_da)
    cell_lon, cell_lat = np.meshgrid(lon, lat)
    flat_lon = cell_lon.ravel()
    flat_lat = cell_lat.ravel()
    inside = shapely.contains_xy(poly, flat_lon, flat_lat)
    mask_flat = np.where(inside, 1.0, np.nan)
    return mask_flat.reshape(cell_lon.shape)


# Cached well-location array.  Loaded once on first call, reused for
# every panel.  Reading from disk per panel is cheap but the cache makes
# the intent explicit and avoids 12 redundant disk hits.
_WELL_LOCS_CACHE = None


def _well_locations():
    """Return the (109, 2) array of (lon, lat) well locations written by
    `07a_backstrip_all_nwshelf.py`.  Cached on first call.
    """
    global _WELL_LOCS_CACHE
    if _WELL_LOCS_CACHE is None:
        loc_path = os.path.join(OUT_BASE, "all_wells_locations.txt")
        if not os.path.exists(loc_path):
            raise FileNotFoundError(
                f"{loc_path} not found -- run 07a_backstrip_all_nwshelf.py "
                "first."
            )
        _WELL_LOCS_CACHE = np.loadtxt(loc_path, usecols=(1, 2), skiprows=1)
        print(f"  loaded {len(_WELL_LOCS_CACHE)} well locations from "
              f"{os.path.basename(loc_path)} (cached for the rest of the run)")
    return _WELL_LOCS_CACHE


def make_continuous_masked_grid(rate_path, region, spacing=REGRID_SPACING,
                                tension=REGRID_TENSION):
    """Return a continuous, well-coverage-masked rate grid for one panel.

    Steps
    -----
    1. Read the (gappy) ``rate_<t>.nc`` and extract every finite
       ``(lon, lat, rate)`` triplet.  Script 07a produced these by
       gridding per-well rates with `surface` + a per-well buffer mask,
       so cells outside that buffer are NaN here.  Coord names may be
       ``lon``/``lat`` or ``x``/``y`` depending on which pyGMT tool
       wrote them.
    2. ``pygmt.blockmedian`` consolidates the points at the regrid cell
       spacing (drops the per-input-cell noise that would otherwise pin
       the surface solver into wiggles).
    3. ``pygmt.surface`` interpolates a smooth, continuous grid over the
       entire `region` at `spacing` with `tension`.  This grid is
       defined everywhere -- but it has extrapolated values far from
       the data.
    4. Mask = concave hull (alpha shape) of the 109 well locations.  A
       single connected polygon following the outline of the well
       cluster -- explicitly NOT a union of circles around individual
       wells, which would reproduce the gappy "blob" appearance of the
       earlier nearest-neighbour mask.
    5. Return surface * mask -- continuous over the shelf, NaN outside.

    Returns
    -------
    xarray.DataArray | None
        The masked continuous grid, or None if there were no finite
        input cells to interpolate from.
    """
    ds = xr.open_dataset(rate_path)
    z = ds["z"].values
    lon, lat = _lon_lat_arrays(ds)
    lon2, lat2 = np.meshgrid(lon, lat)
    finite = np.isfinite(z)
    if not finite.any():
        return None
    pts = np.column_stack([lon2[finite].ravel(),
                           lat2[finite].ravel(),
                           z[finite].ravel()])
    # blockmedian consolidates duplicates that fall in the same regrid cell.
    bm = pygmt.blockmedian(data=pts, region=region, spacing=spacing)
    # surface fills the whole region (smooth minimum-curvature spline).
    cont = pygmt.surface(data=bm, region=region, spacing=spacing,
                         tension=tension)
    # Build the coverage mask from the concave hull of the well
    # locations -- a single connected polygon, not a union of circles.
    wells = _well_locations()
    mask = _concave_hull_mask(cont, wells[:, 0], wells[:, 1], region)
    masked = cont.copy()
    masked.values = cont.values * mask
    return masked


def plot_matrix(picked_times):
    """Render the 4-row x 3-col rate-map matrix.

    Refactored for reliability on macOS Apple-Silicon:
      * The two CPTs (imola for A/C, vik for delta) are pre-built and
        saved to disk *before* the subplot context; grdimage receives a
        filename rather than relying on the session-default CPT inside
        the subplot.
      * Each panel emits a short progress print so a hang can be located
        precisely without waiting for the whole 12-panel matrix.
      * `fig.coast(resolution="c")` (crude) is used because the bundled
        GMT "low" / "intermediate" GSHHG dataset can trigger a one-off
        download/cache build that visibly hangs on a fresh conda env;
        "c" is in-binary and instant.  This is plenty of detail given
        the map scale.
      * Per-panel `frame=` is dropped because the subplot already sets
        the shared axes; redundant frame= calls inside the subplot
        have been observed to confuse GMT's frame state on M3.
    """
    overall_t0 = time.time()
    print(f"\nRendering rate-map matrix for times {picked_times} Ma")
    print(f"  output base: {os.path.join(OUTPUT_DIR, 'fig10_nwshelf_rate_maps_combined')}")

    # --- 1. CPTs ------------------------------------------------------------
    cpt_dir = os.path.join(OUT_BASE, "cpts")
    os.makedirs(cpt_dir, exist_ok=True)
    rate_cpt = os.path.join(cpt_dir, "fig10_rate.cpt")
    delta_cpt = os.path.join(cpt_dir, "fig10_delta.cpt")
    print(f"  building {RATE_CMAP} CPT  {RATE_SERIES} -> {rate_cpt} "
          "(background clamped so rates > 50 m/Myr render as the 50-m/Myr "
          "colour rather than white)")
    _save_cpt(RATE_CMAP, RATE_SERIES, rate_cpt, background=True)
    print(f"  building {DELTA_CMAP} CPT {DELTA_SERIES} -> {delta_cpt}"
          " (reversed)")
    _save_cpt(DELTA_CMAP, DELTA_SERIES, delta_cpt, reverse=True)

    # --- 2. Figure + global config ------------------------------------------
    fig = pygmt.Figure()
    pygmt.config(
        FONT_TITLE="18p,Helvetica-Bold",
        FONT_LABEL="16p",
        FONT_ANNOT="13p",
        MAP_FRAME_TYPE="plain",
        COLOR_NAN="240/240/240",
    )

    # Pre-load the 109 well locations once so they can be over-plotted
    # on every panel as small white-filled circles (so the reader can
    # see where the data actually constrain the gridded rate field).
    wells = _well_locations()
    well_lons = wells[:, 0]
    well_lats = wells[:, 1]

    # --- 3. 4x3 subplot loop ------------------------------------------------
    with fig.subplot(
        nrows=4, ncols=3,
        figsize=FIG_FIGSIZE,
        margins=SUBPLOT_MARGINS,
        sharex="b", sharey="l",
        frame=["WSne", "xa5f1+lLongitude", "ya5f1+lLatitude"],
    ):
        for row, t_ma in enumerate(picked_times):
            # Display tags relabel C -> "B" and D -> "C" so the figure
            # reads as A / B / C instead of the internal A / C / D
            # (which still apply to directory names and the rest of
            # the pipeline).
            for col, (sub_dir, cpt_path, tag) in enumerate((
                (A_DIR, rate_cpt, "A"),
                (C_DIR, rate_cpt, "B"),
                (D_DIR, delta_cpt, "C = B - A"),
            )):
                idx = row * 3 + col
                path = os.path.join(sub_dir, f"rate_{t_ma}.nc")
                panel_t0 = time.time()
                print(f"  panel [{row + 1}/4 {tag:<9}] {t_ma:>4} Ma -> ",
                      end="", flush=True)

                with fig.set_panel(idx):
                    # Always lay down the basemap so the panel is sized
                    # consistently even when the grid is missing/empty.
                    fig.basemap(region=REGION, projection="M?")

                    if _grid_is_drawable(path):
                        # Continuous gridding step: blockmedian -> surface
                        # -> concave-hull mask of the 109 well locations
                        # so the smooth surface is shown only inside the
                        # actual outline of the well cluster.
                        grid = make_continuous_masked_grid(path, REGION)
                        if grid is not None:
                            fig.grdimage(
                                grid=grid, region=REGION,
                                cmap=cpt_path, nan_transparent=True,
                            )
                        else:
                            fig.text(
                                x=(REGION[0] + REGION[1]) / 2,
                                y=(REGION[2] + REGION[3]) / 2,
                                text="no finite cells",
                                font="12p,Helvetica", justify="MC",
                            )
                    else:
                        fig.text(
                            x=(REGION[0] + REGION[1]) / 2,
                            y=(REGION[2] + REGION[3]) / 2,
                            text=f"{os.path.basename(path)} missing/empty",
                            font="12p,Helvetica", justify="MC",
                        )

                    # Crude shorelines: in-binary, no network/cache work.
                    fig.coast(shorelines="0.3p,black", resolution="c")
                    # Well markers: tiny white-filled circles with a thin
                    # black outline so they read against any background
                    # colour in the rate field.  Drawn before the panel
                    # title so the title box sits on top of any wells
                    # underneath it.
                    fig.plot(
                        x=well_lons, y=well_lats,
                        style="c0.10c", fill="white", pen="0.25p,black",
                    )
                    fig.text(
                        x=REGION[0] + 0.4, y=REGION[3] - 0.4,
                        text=f"{tag}  -  {t_ma} Ma",
                        font=PANEL_TITLE_FONT,
                        justify="TL", fill="white@30", pen="0.25p,black",
                    )

                print(f"{time.time() - panel_t0:5.1f} s")

    print(f"  all 12 panels rendered in {time.time() - overall_t0:.1f} s")

    # --- 4. Two independent colourbars, stacked, full figure width ---------
    # Rates (A, C) use batlow; delta (D) uses vik.  The two ranges are
    # NOT the same (delta is typically much smaller than the absolute
    # rate) so a single shared bar would compress one of them.  Both
    # colourbars are now drawn at the *full figure width* (20 cm, leaving
    # ~0.5 cm gutter on each side of the 21 cm panel matrix) and stacked
    # vertically below the figure -- this makes the colour ramps and
    # tick marks readable at print size instead of being squeezed into
    # two 10 cm half-width bars.
    #
    # The rate bar (top of the stack) puts its annotations ABOVE the bar
    # (+ma) so they sit in the gap between the panel matrix and the
    # rate bar; the delta bar (bottom of the stack) keeps its default
    # below-the-bar annotation position so its labels never overlap with
    # the rate bar above.
    #
    # +e on each position string adds side-arrow triangles indicating
    # values beyond the CPT range are clamped to the end colours.
    CBAR_WIDTH = "20c"                  # cm; figure is 21 cm wide
    CBAR_HEIGHT = "0.5c"
    RATE_CBAR_Y = "1.6c"                # offset below the figure
    DELTA_CBAR_Y = "3.6c"               # ~2 cm gap below the rate bar
    print("  drawing rate (cubhelix) colourbar -- top of stack, full width")
    fig.colorbar(
        cmap=rate_cpt,
        position=(f"JBC+w{CBAR_WIDTH}/{CBAR_HEIGHT}+h"
                  f"+o0c/{RATE_CBAR_Y}+ma+e"),
        # No literal double quotes around the +l label -- GMT renders
        # them as part of the text on this build.  The pyGMT string is
        # already passed atomically so the quotes are not needed to
        # escape spaces.
        frame=["a10f5", "x+lSubsidence rate (m/Myr)"],
    )
    print("  drawing delta (vik) colourbar -- bottom of stack, full width")
    fig.colorbar(
        cmap=delta_cpt,
        position=(f"JBC+w{CBAR_WIDTH}/{CBAR_HEIGHT}+h"
                  f"+o0c/{DELTA_CBAR_Y}+e"),
        frame=[
            f"a{int((DELTA_SERIES[1] - DELTA_SERIES[0]) / 4)}f"
            f"{int((DELTA_SERIES[1] - DELTA_SERIES[0]) / 20)}",
            "x+lRate difference B - A (m/Myr)",
        ],
    )

    # --- 5. Save ------------------------------------------------------------
    base = os.path.join(OUTPUT_DIR, "fig10_nwshelf_rate_maps_combined")
    save_t0 = time.time()
    print("  saving PNG ...")
    fig.savefig(base + ".png", dpi=300)
    print(f"    wrote {base}.png ({time.time() - save_t0:.1f} s)")
    save_t0 = time.time()
    print("  saving PDF ...")
    fig.savefig(base + ".pdf")
    print(f"    wrote {base}.pdf ({time.time() - save_t0:.1f} s)")
    print(f"  TOTAL plot_matrix time: {time.time() - overall_t0:.1f} s")


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Hand-picked four times that correspond to the four 5-Myr intervals
# identified in Fig 11 panel (c) as having the largest mean |D| across
# the NW Shelf.  Each pick is the middle integer age of its 5-Myr
# interval (e.g. 12 Ma = centre of the 10..14 Ma window).
#
# Set to None to fall back to the original greedy non-maximum
# suppression picker (pick_max_delta_times) acting on the per-time
# stats.
# ----------------------------------------------------------------------------
HARDCODED_PICKED_TIMES = [12, 27, 112, 132]   # Ma; youngest -> oldest


def main():
    if not os.path.isdir(D_DIR):
        sys.exit(
            f"{D_DIR} not found -- run 07a_backstrip_all_nwshelf.py first."
        )
    times, stats = collect_per_time_stats()
    if len(times) == 0:
        sys.exit(f"No rate-difference grids found under {D_DIR}.")

    write_stats_csv(times, stats)

    if HARDCODED_PICKED_TIMES is not None:
        picked = sorted(HARDCODED_PICKED_TIMES, reverse=True)  # old -> young
        # Plot rows go oldest at the top; the matrix code reverses
        # this list when iterating, so honour the same convention.
        print(f"\nUsing HARDCODED_PICKED_TIMES: {picked} Ma "
              "(four 5-Myr intervals from Fig 11 panel (c) "
              "with largest mean |D|)")
    else:
        picked = pick_max_delta_times(
            times, stats[METRIC_NAME], n=N_TIMES, min_sep=DT_TIME_SEP)
        print(f"\nPicked {N_TIMES} max-{METRIC_NAME} times: {picked} Ma  "
              f"(min separation = {DT_TIME_SEP} Myr)")

    # Persist the picks so build_paper.js, Fig 12 (12_dt_maps.py) and
    # reviewers can read them.
    picks_path = os.path.join(OUT_BASE, "rate_delta_picked_times.txt")
    with open(picks_path, "w") as f:
        if HARDCODED_PICKED_TIMES is not None:
            f.write(
                "# Hand-picked from the four 5-Myr intervals with the\n"
                "# largest mean |D| in Fig 11 panel (c) (see\n"
                "# 11_rate_boxplots.py).\n"
                "# Used as the rows of fig10_nwshelf_rate_maps_combined.*\n"
                "# and the four time slices of fig12_dt_maps_2x2.*.\n"
            )
        else:
            f.write(
                f"# Times (Ma) at which spatial-mean(|D|) is maximal,\n"
                f"# with a minimum separation of {DT_TIME_SEP} Myr.\n"
                "# Used as the rows of fig10_nwshelf_rate_maps_combined.*\n"
                "# and the four time slices of fig12_dt_maps_2x2.*.\n"
            )
        for t in picked:
            try:
                i = int(np.where(times == t)[0][0])
                metric_val = stats[METRIC_NAME][i]
                f.write(f"{t}\t{metric_val:.3f}\n")
            except (IndexError, ValueError):
                # Time isn't in `times` (hand-picked time outside the
                # per-time stats sample); just write the time.
                f.write(f"{t}\n")
    print(f"  wrote {picks_path}")

    plot_matrix(picked)


if __name__ == "__main__":
    main()
