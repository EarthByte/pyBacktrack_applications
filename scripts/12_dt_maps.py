#!/usr/bin/env python3
"""
Figure 12 - Raw dynamic-topography maps over the NW Shelf at the same
four times picked for Fig. 10.

For each picked time t the field plotted is the bundled D10_gmcm9
dynamic-topography elevation anomaly (Braz et al., 2021), positive-
up relative to present-day sea level:

    positive (warm) = location was higher than today (uplifted)
    negative (cool) = location was lower than today (subsided down)
    zero            = same as present

The bundled D10_gmcm9 grids are sampled at 1 deg resolution and at
roughly 5 Myr time slices.  For each picked time t we:

  1. Linearly interpolate the 5 Myr grids in time to obtain the DT
     field at exactly that age.
  2. Extract the (lon, lat, DT) point cloud over the NW Shelf bounding
     box and feed it to GMT `blockmedian` (0.05 deg cells) -> `surface`
     (tension 0.5) -- the same upsampling pipeline that produces the
     well-derived rate maps of Fig. 10, so the visual resolution and
     smoothness match column-for-column.
  3. Mask with the same concave hull of the 109 well locations used
     in Fig. 10 so the spatial footprint of the panels matches exactly.

No well markers are over-plotted: the field is a pure model output and
does not depend on the well distribution.

Inputs:
    pybacktrack bundle:
        bundle_data/dynamic_topography/models/Braz2021/D10_gmcm9/<t>.00.nc
    figures/output/nwshelf_subsidence/all_wells_locations.txt
    figures/output/nwshelf_subsidence/rate_delta_picked_times.txt
        (the four times used as panels of Fig. 10 -- if missing, falls
         back to the hard-coded default [130, 100, 65, 10] Ma)

Outputs:
    figures/output/fig12_dt_maps_2x2.png
    figures/output/fig12_dt_maps_2x2.pdf
"""
import importlib.util
import os
import re
import sys
import time

import numpy as np
import pygmt
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (OUTPUT_DIR, DYNAMIC_TOPOGRAPHY_MODEL, MAX_ANALYSIS_AGE,
                    DT_GRIDS_DIR)

# Re-use Fig 10's helpers so the spatial-resolution pipeline and the
# concave-hull mask match exactly.
_spec = importlib.util.spec_from_file_location(
    "_fig10",
    os.path.join(os.path.dirname(__file__), "10_rate_maps.py"),
)
_fig10 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_fig10)
_concave_hull_mask = _fig10._concave_hull_mask
_well_locations = _fig10._well_locations
_save_cpt = _fig10._save_cpt
_lon_lat_arrays = _fig10._lon_lat_arrays
REGRID_SPACING = _fig10.REGRID_SPACING       # 0.05 deg (~5.5 km in lon)
REGRID_TENSION = _fig10.REGRID_TENSION       # GMT surface tension 0.5

# ----------------------------------------------------------------------------
# Tunables
# ----------------------------------------------------------------------------
OUT_BASE = os.path.join(OUTPUT_DIR, "nwshelf_subsidence")
REGION = [113, 132, -22, -8]                  # match Fig 9

# Buffer added when subsetting the global DT grids so surface has
# data points outside the masked region (avoids edge ringing).
GRID_BUFFER_DEG = 4.0

# GMT's `polar` diverging palette: blue = negative DT (subsided down
# relative to today), white = 0 (same as today), red = positive DT
# (uplifted relative to today).  polar has a hinge baked into its
# master CPT at z=0, so the white midpoint anchors to z=0 regardless
# of how asymmetric the series range is -- exactly what we want for
# the masked NW Shelf DT field, which at the picked times is
# dominantly negative (down to ~-350 m) but has a moderate positive
# band at ~112 Ma (up to ~+50 m).  Cells outside [-350, +50] clamp to
# the endpoint colours, with triangular arrows on BOTH ends of the
# colourbar (`+ebf`) flagging the clamp.
#
# Tried polar first (custom Python-built CPT to pin white at z=0); the
# user asked to switch to GMT's hypsometric `topo` palette, which has
# a HARD HINGE built into its master CPT at sea level (z=0).  GMT
# 6.5/6.6 honours intrinsic hard hinges automatically when `-T` spans
# across z=0 -- the two halves of the master CPT are stretched
# independently to [lo, 0] and [0, hi], so the sea-level colour break
# lands exactly at z=0 regardless of how asymmetric the range is.
# Semantically nice for DT: negative DT renders as ocean blues,
# positive DT as land greens/browns.  No custom Python path needed.
DT_CMAP = "topo"
DT_CMAP_REVERSE = False
DT_SERIES = (-350, 50, 10)                    # m
DT_CMAP_HINGE = 0                             # auto-honoured by topo's intrinsic hard hinge
DT_CBAR_TICKS = "a50f25"                      # major every 50 m
DT_CBAR_END_ARROWS = "+ebf"                   # both-end arrows

PANEL_TITLE_FONT = "16p,Helvetica-Bold"

# Layout (2 cols x 2 rows).  Per the same Mercator-aspect-driven sizing
# argument as Fig 9, the figure height must equal n_rows * panel_h
# (no dead vertical space inside cells).  For REGION=[113,132,-22,-8]
# with 2 cols and ~8 c column width, each panel is ~6.1c tall;
# 2 rows + bottom-axis labels -> ~13c.
FIG_FIGSIZE = ("17c", "13c")
SUBPLOT_MARGINS = ["0.35c", "0c"]
CBAR_WIDTH = "16c"
CBAR_HEIGHT = "0.5c"
CBAR_Y_OFFSET = "1.6c"

# Fallback in case the picked-times sidecar is missing.
DEFAULT_PICKED_TIMES = [130, 100, 65, 10]


# ----------------------------------------------------------------------------
# Locate the regional D10_gmcm9 grids in the repo
# ----------------------------------------------------------------------------
def _bundled_dt_dir(model_name=DYNAMIC_TOPOGRAPHY_MODEL):
    """Return the directory holding the per-time DT NetCDFs.

    The repo ships a region-cropped copy of the Braz et al. (2021)
    ``D10_gmcm9`` model at ``data/grids/dynamic_topography/D10_gmcm9/``
    (written by ``tools/populate_data.py`` from the pybacktrack bundle).
    All slice filenames match the upstream convention -- ``<time>.00.nc``
    -- so the existing per-time interpolation logic is unchanged.

    If the directory is missing we fall back to the bundled global copy
    inside the installed pybacktrack so a stock install still renders.
    """
    if os.path.isdir(DT_GRIDS_DIR) and any(
            fn.endswith(".nc") for fn in os.listdir(DT_GRIDS_DIR)):
        return DT_GRIDS_DIR

    # Fallback: search the bundled pybacktrack copy.
    try:
        import pybacktrack
    except ImportError:
        raise FileNotFoundError(
            f"No regional DT grids at {DT_GRIDS_DIR} and pybacktrack is "
            "not installed.  Run `python tools/populate_data.py` from the "
            "repo root, or install pybacktrack and try again.")
    base = os.path.join(
        os.path.dirname(pybacktrack.__file__),
        "bundle_data", "dynamic_topography", "models",
    )
    for sub in (("Braz2021", model_name), (f"{model_name}.grids",),
                (model_name,)):
        cand = os.path.join(base, *sub)
        if os.path.isdir(cand) and any(
                fn.endswith(".nc") for fn in os.listdir(cand)):
            return cand
    raise FileNotFoundError(
        f"could not find DT grids for model {model_name} at {DT_GRIDS_DIR} "
        f"or under the pybacktrack bundle ({base})")


# ----------------------------------------------------------------------------
# Build a (time, lat, lon) DT stack covering the NW Shelf region
# ----------------------------------------------------------------------------
def _load_dt_stack():
    """Read every <age>.00.nc grid in [0, MAX_ANALYSIS_AGE], subset to the
    NW Shelf bbox + buffer, and return:

        grid_times : sorted ndarray of grid ages (Ma)
        dt_stack   : (n_times, n_lat, n_lon) float32 array, m
        lon_vec    : (n_lon,) ndarray of lon centres
        lat_vec    : (n_lat,) ndarray of lat centres
    """
    dt_dir = _bundled_dt_dir()
    print(f"  bundled D10_gmcm9 grids: {dt_dir}")

    times = []
    for fn in sorted(os.listdir(dt_dir)):
        m = re.match(r"^([0-9]+\.[0-9]+)\.nc$", fn)
        if m is None:
            continue
        t = float(m.group(1))
        if 0.0 <= t <= MAX_ANALYSIS_AGE:
            times.append(t)
    times = sorted(times)
    if not times:
        raise RuntimeError(f"no DT grids found in [0, {MAX_ANALYSIS_AGE}] Ma "
                           f"under {dt_dir}")
    print(f"  loaded {len(times)} DT time slices: "
          f"{times[0]:.0f} .. {times[-1]:.0f} Ma "
          f"(median spacing "
          f"{float(np.median(np.diff(times))):.1f} Myr)")

    lon_min, lon_max = REGION[0] - GRID_BUFFER_DEG, REGION[1] + GRID_BUFFER_DEG
    lat_min, lat_max = REGION[2] - GRID_BUFFER_DEG, REGION[3] + GRID_BUFFER_DEG

    stack = []
    lon_vec = lat_vec = None
    for t in times:
        path = os.path.join(dt_dir, f"{t:.2f}.nc")
        ds = xr.open_dataset(path)
        sub = ds.z.sel(lon=slice(lon_min, lon_max),
                       lat=slice(lat_min, lat_max))
        if lon_vec is None:
            lon_vec = sub.lon.values
            lat_vec = sub.lat.values
        stack.append(sub.values.astype(np.float32))
    return (np.array(times, dtype=float),
            np.stack(stack, axis=0),
            lon_vec, lat_vec)


def _dt_at_time(t, grid_times, dt_stack):
    """Linearly interpolate the DT field at age t (Ma) from the
    bundled grids.
    """
    if t <= grid_times[0]:
        return dt_stack[0]
    if t >= grid_times[-1]:
        return dt_stack[-1]
    i_hi = int(np.searchsorted(grid_times, t, side="left"))
    if grid_times[i_hi] == t:
        return dt_stack[i_hi]
    i_lo = i_hi - 1
    f = (t - grid_times[i_lo]) / (grid_times[i_hi] - grid_times[i_lo])
    return dt_stack[i_lo] * (1.0 - f) + dt_stack[i_hi] * f


# ----------------------------------------------------------------------------
# Picked-times sidecar
# ----------------------------------------------------------------------------
def _picked_times():
    """Read the same four times Fig 9 picked.  Falls back to a
    hard-coded default if the sidecar is missing.
    """
    path = os.path.join(OUT_BASE, "rate_delta_picked_times.txt")
    if not os.path.exists(path):
        print(f"  ! {path} missing -- using fallback "
              f"picked times {DEFAULT_PICKED_TIMES}")
        return list(DEFAULT_PICKED_TIMES)
    picks = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                picks.append(int(float(line.split()[0])))
            except (ValueError, IndexError):
                continue
    picks = sorted(set(picks), reverse=True)
    if len(picks) != 4:
        print(f"  ! found {len(picks)} picked times in {path}; "
              f"falling back to {DEFAULT_PICKED_TIMES}")
        return list(DEFAULT_PICKED_TIMES)
    return picks


# ----------------------------------------------------------------------------
# Upsample one time slice through Fig 9's blockmedian -> surface pipeline
# ----------------------------------------------------------------------------
def upsample_dt_field(dt_field, lon_vec, lat_vec):
    """Re-grid the 1 deg DT field at one time slice through GMT
    `blockmedian` (REGRID_SPACING) -> `surface` (REGRID_TENSION) so the
    rendered map has the same visual resolution as the well-based rate
    maps in Fig 9.

    Returns an xarray.DataArray over REGION with the same target
    spacing as Fig 9's continuous rate field.
    """
    lon2, lat2 = np.meshgrid(lon_vec, lat_vec)
    finite = np.isfinite(dt_field)
    if not finite.any():
        return None
    pts = np.column_stack([
        lon2[finite].ravel(),
        lat2[finite].ravel(),
        dt_field[finite].ravel().astype(float),
    ])
    bm = pygmt.blockmedian(data=pts, region=REGION, spacing=REGRID_SPACING)
    cont = pygmt.surface(
        data=bm, region=REGION,
        spacing=REGRID_SPACING, tension=REGRID_TENSION,
    )
    return cont


def make_masked_dt_grid(t, grid_times, dt_stack, lon_vec, lat_vec, mask):
    """Return a continuous, well-coverage-masked DT field at age t.

    Steps:
        1. Linearly interpolate the 1-deg bundled DT grids at age t.
        2. Upsample with blockmedian -> surface to Fig 9's 0.05-deg grid.
        3. Multiply by the concave-hull mask (already 0.05-deg aligned).
    """
    raw = _dt_at_time(t, grid_times, dt_stack)            # (n_lat, n_lon)
    cont = upsample_dt_field(raw, lon_vec, lat_vec)
    if cont is None:
        return None
    masked = cont.copy()
    masked.values = cont.values * mask
    return masked


# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------
def plot_matrix(picked_times, grid_times, dt_stack, lon_vec, lat_vec, mask):
    """Render the 2 x 2 DT matrix at picked_times."""
    overall_t0 = time.time()
    print(f"\nRendering raw-DT 2x2 matrix at times {picked_times} Ma")

    cpt_dir = os.path.join(OUT_BASE, "cpts")
    os.makedirs(cpt_dir, exist_ok=True)
    dt_cpt = os.path.join(cpt_dir, "fig12_dt.cpt")
    print(f"  building {DT_CMAP} CPT {DT_SERIES} -> {dt_cpt}")
    _save_cpt(DT_CMAP, DT_SERIES, dt_cpt,
              reverse=DT_CMAP_REVERSE, hinge=DT_CMAP_HINGE)

    fig = pygmt.Figure()
    pygmt.config(
        FONT_TITLE="18p,Helvetica-Bold",
        FONT_LABEL="16p",
        FONT_ANNOT="13p",
        MAP_FRAME_TYPE="plain",
        COLOR_NAN="240/240/240",
    )

    with fig.subplot(
        nrows=2, ncols=2,
        figsize=FIG_FIGSIZE,
        margins=SUBPLOT_MARGINS,
        sharex="b", sharey="l",
        frame=["WSne", "xa5f1+lLongitude", "ya5f1+lLatitude"],
    ):
        for idx, t_ma in enumerate(picked_times):
            panel_t0 = time.time()
            print(f"  panel [{idx + 1}/4] {t_ma:>4} Ma -> ", end="",
                  flush=True)
            with fig.set_panel(idx):
                fig.basemap(region=REGION, projection="M?")
                grid = make_masked_dt_grid(
                    float(t_ma), grid_times, dt_stack,
                    lon_vec, lat_vec, mask,
                )
                if grid is not None:
                    fig.grdimage(
                        grid=grid, region=REGION,
                        cmap=dt_cpt, nan_transparent=True,
                    )
                else:
                    fig.text(
                        x=(REGION[0] + REGION[1]) / 2,
                        y=(REGION[2] + REGION[3]) / 2,
                        text="no finite DT cells",
                        font="12p,Helvetica", justify="MC",
                    )
                # Intermediate-resolution coastline -- the NW Shelf
                # outline reads with too few vertices at "c" (crude).
                fig.coast(shorelines="0.3p,black", resolution="i")
                fig.text(
                    x=REGION[0] + 0.4, y=REGION[3] - 0.4,
                    text=f"DT  -  {t_ma} Ma",
                    font=PANEL_TITLE_FONT,
                    justify="TL", fill="white@30", pen="0.25p,black",
                )
            print(f"{time.time() - panel_t0:5.1f} s")

    print("  drawing colourbar")
    fig.colorbar(
        cmap=dt_cpt,
        position=(f"JBC+w{CBAR_WIDTH}/{CBAR_HEIGHT}+h"
                  f"+o0c/{CBAR_Y_OFFSET}+ma{DT_CBAR_END_ARROWS}"),
        frame=[DT_CBAR_TICKS, "x+lDynamic topography (m)"],
    )

    base = os.path.join(OUTPUT_DIR, "fig12_dt_maps_2x2")
    print("  saving PNG ...")
    fig.savefig(base + ".png", dpi=300)
    print(f"    wrote {base}.png")
    print("  saving PDF ...")
    fig.savefig(base + ".pdf")
    print(f"    wrote {base}.pdf")
    print(f"  TOTAL: {time.time() - overall_t0:.1f} s")


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    picked = _picked_times()
    print(f"  picked times = {picked} Ma (oldest first)")

    grid_times, dt_stack, lon_vec, lat_vec = _load_dt_stack()

    # Mask aligned with Fig 9's blockmedian/surface output spacing.
    n_lat = int(round((REGION[3] - REGION[2]) / REGRID_SPACING)) + 1
    n_lon = int(round((REGION[1] - REGION[0]) / REGRID_SPACING)) + 1
    template = xr.DataArray(
        np.zeros((n_lat, n_lon), dtype=np.float32),
        coords={
            "lat": np.linspace(REGION[2], REGION[3], n_lat),
            "lon": np.linspace(REGION[0], REGION[1], n_lon),
        },
        dims=("lat", "lon"),
    )
    wells = _well_locations()
    mask = _concave_hull_mask(template, wells[:, 0], wells[:, 1], REGION)
    print(f"  mask: {int(np.isfinite(mask).sum())} / {mask.size} cells active")

    plot_matrix(picked, grid_times, dt_stack, lon_vec, lat_vec, mask)


if __name__ == "__main__":
    main()
