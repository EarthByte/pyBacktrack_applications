#!/usr/bin/env python3
"""
Figure 12 video - 1 Myr animation of raw dynamic topography over the
NW Shelf.

Per-Myr companion to `12_dt_maps.py`: renders one single-panel frame
per Myr in [0, MAX_ANALYSIS_AGE] showing the bundled D10_gmcm9
dynamic-topography elevation anomaly (Braz et al., 2021), linearly
interpolated in time between the 5 Myr bundled grids and upsampled
spatially through the same `blockmedian -> surface` pipeline (0.05
deg, tension 0.5) and concave-hull mask used in Fig 10 and Fig 12.

Frame styling -- topo hypsometric colourmap, -350 to +50 m range
with sea-level hinge at z=0 and `+ebf` clamp arrows on the
colourbar, shoreline, panel label -- matches `12_dt_maps.py` so the
animation is visually continuous with the static figure.

Outputs:
    figures/output/dt_field/frames/frame_NNNN.png
    figures/output/fig12_dt_video.mp4
        ^ geological-time backward (MAX_ANALYSIS_AGE -> 0)
    figures/output/fig12_dt_video_forward.mp4
        ^ chronological forward (0 -> MAX_ANALYSIS_AGE)

Usage:
    python3 12a_dt_video.py                  # full 0..MAX_ANALYSIS_AGE
    python3 12a_dt_video.py --time-step 5    # quick coarse check
    python3 12a_dt_video.py --max-age 100    # truncate
"""
import argparse
import importlib.util
import os
import shutil
import subprocess
import sys
import time

import numpy as np
import pygmt
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import OUTPUT_DIR, MAX_ANALYSIS_AGE


def _import_sibling(name, filename):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(os.path.dirname(__file__), filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_fig12 = _import_sibling("_fig12", "12_dt_maps.py")
_fig10 = _import_sibling("_fig10", "10_rate_maps.py")

REGION = _fig12.REGION
DT_CMAP = _fig12.DT_CMAP
DT_CMAP_REVERSE = _fig12.DT_CMAP_REVERSE
DT_CMAP_HINGE = _fig12.DT_CMAP_HINGE
DT_SERIES = _fig12.DT_SERIES
DT_CBAR_TICKS = _fig12.DT_CBAR_TICKS
DT_CBAR_END_ARROWS = _fig12.DT_CBAR_END_ARROWS
PANEL_TITLE_FONT = _fig12.PANEL_TITLE_FONT
REGRID_SPACING = _fig12.REGRID_SPACING
_save_cpt = _fig12._save_cpt
_load_dt_stack = _fig12._load_dt_stack
make_masked_dt_grid = _fig12.make_masked_dt_grid
_concave_hull_mask = _fig10._concave_hull_mask
_well_locations = _fig10._well_locations


# ----------------------------------------------------------------------------
# Tunables
# ----------------------------------------------------------------------------
OUT_BASE = os.path.join(OUTPUT_DIR, "dt_field")
FRAME_DIR = os.path.join(OUT_BASE, "frames")
FRAMERATE = 6
FRAME_WIDTH = "12c"                   # Mercator NW Shelf panel
CBAR_WIDTH = "11c"
CBAR_HEIGHT = "0.4c"
CBAR_Y_OFFSET = "1.2c"


# ----------------------------------------------------------------------------
# Args
# ----------------------------------------------------------------------------
def _parse_args():
    p = argparse.ArgumentParser(
        description=("1 Myr raw-DT animation over the NW Shelf, "
                     "stitched into MP4 via ffmpeg.  Frames already "
                     "on disk are reused unless --force is passed."),
    )
    p.add_argument("--time-step", type=int, default=1,
                   help="Render every N-th Myr (default 1).")
    p.add_argument("--max-age", type=float, default=MAX_ANALYSIS_AGE,
                   help=f"Stop at this age (default {MAX_ANALYSIS_AGE}).")
    p.add_argument("--force", action="store_true",
                   help="Wipe the frame cache and re-render every frame.")
    p.add_argument("--skip-video", action="store_true",
                   help="Render frames but don't run ffmpeg.")
    p.add_argument("--keep-frames", action="store_true",
                   help="Do not delete the per-frame PNGs after MP4 stitch.")
    return p.parse_args()


# ----------------------------------------------------------------------------
# Frame renderer
# ----------------------------------------------------------------------------
def _render_frame(t_ma, frame_path, grid_times, dt_stack,
                  lon_vec, lat_vec, mask, dt_cpt):
    """Render a single-panel PNG for one time slice."""
    fig = pygmt.Figure()
    pygmt.config(
        FONT_TITLE="14p,Helvetica-Bold",
        FONT_LABEL="13p",
        FONT_ANNOT="10p",
        MAP_FRAME_TYPE="plain",
        COLOR_NAN="240/240/240",
    )
    fig.basemap(region=REGION, projection=f"M{FRAME_WIDTH}",
                frame=["WSne", "xa5f1+lLongitude", "ya5f1+lLatitude"])
    grid = make_masked_dt_grid(
        float(t_ma), grid_times, dt_stack, lon_vec, lat_vec, mask)
    if grid is not None:
        fig.grdimage(grid=grid, region=REGION,
                     cmap=dt_cpt, nan_transparent=True)
    fig.coast(shorelines="0.3p,black", resolution="i")
    fig.text(
        x=REGION[0] + 0.4, y=REGION[3] - 0.4,
        text=f"DT  -  {int(t_ma)} Ma",
        font=PANEL_TITLE_FONT, justify="TL",
        fill="white@30", pen="0.25p,black",
    )
    fig.colorbar(
        cmap=dt_cpt,
        position=(f"JBC+w{CBAR_WIDTH}/{CBAR_HEIGHT}+h"
                  f"+o0c/{CBAR_Y_OFFSET}+ma{DT_CBAR_END_ARROWS}"),
        frame=[DT_CBAR_TICKS, "x+lDynamic topography (m)"],
    )
    fig.savefig(frame_path, dpi=150)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    args = _parse_args()

    if shutil.which("ffmpeg") is None:
        sys.exit("ffmpeg not found in PATH; install ffmpeg first.")

    grid_times, dt_stack, lon_vec, lat_vec = _load_dt_stack()

    # Mask aligned with Fig 9 / Fig 11 spacing.
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

    cpt_dir = os.path.join(OUTPUT_DIR, "nwshelf_subsidence", "cpts")
    os.makedirs(cpt_dir, exist_ok=True)
    dt_cpt = os.path.join(cpt_dir, "fig12_dt.cpt")
    _save_cpt(DT_CMAP, DT_SERIES, dt_cpt,
              reverse=DT_CMAP_REVERSE, hinge=DT_CMAP_HINGE)

    os.makedirs(FRAME_DIR, exist_ok=True)

    times = list(range(0, int(args.max_age) + 1, max(1, args.time_step)))
    # Pre-compute the output paths so the cache audit and the render
    # loop both refer to the same filenames.  Convention here matches
    # the bathymetry videos (02a/02b/03a/03b/10a): frame_0000 = the
    # present (t = 0 Ma) and frame_<MAX_AGE> = the oldest slice.  So
    # the direct-frame-order ffmpeg stitch plays present -> past (the
    # geological "rewind" direction; written to `<base>.mp4`) and the
    # `-vf reverse` stitch plays past -> present (chronological
    # forward; written to `<base>_forward.mp4`).
    frame_paths = [
        os.path.join(FRAME_DIR, f"frame_{i:04d}.png")
        for i in range(len(times))
    ]

    if args.force:
        wiped = 0
        for fn in os.listdir(FRAME_DIR):
            if fn.endswith(".png"):
                os.remove(os.path.join(FRAME_DIR, fn))
                wiped += 1
        if wiped:
            print(f"  --force: wiped {wiped} cached frame(s) from "
                  f"{FRAME_DIR}")

    existing = sum(1 for p in frame_paths if os.path.exists(p))
    if args.force:
        print(f"  --force: rendering {len(times)} frames fresh at "
              f"{args.time_step} Myr step, 0..{int(args.max_age)} Ma")
    elif existing == 0:
        print(f"  cache empty: rendering {len(times)} frames at "
              f"{args.time_step} Myr step, 0..{int(args.max_age)} Ma")
    elif existing == len(times):
        print(f"  cache hit: all {existing} frames already on disk "
              f"(re-stitch only).  Pass --force to regenerate.")
    else:
        print(f"  partial cache: {existing} frame(s) exist, "
              f"{len(times) - existing} will be rendered.  "
              f"Pass --force to regenerate everything.")

    overall_t0 = time.time()
    regenerated = 0
    for i, t in enumerate(times):
        fpath = frame_paths[i]
        if os.path.exists(fpath) and not args.force:
            if i % 25 == 0:
                print(f"  [{i+1:3d}/{len(times)}] t={t:>4} Ma "
                      f"-> {os.path.basename(fpath)} (cached)")
            continue
        frame_t0 = time.time()
        _render_frame(t, fpath, grid_times, dt_stack,
                      lon_vec, lat_vec, mask, dt_cpt)
        regenerated += 1
        if i % 25 == 0:
            print(f"  [{i+1:3d}/{len(times)}] t={t:>4} Ma "
                  f"-> {os.path.basename(fpath)} "
                  f"({time.time() - frame_t0:.1f} s)")
    print(f"  frames ready in {time.time() - overall_t0:.1f} s "
          f"({regenerated} freshly rendered, "
          f"{len(times) - regenerated} reused)")

    if args.skip_video:
        print("  --skip-video: stopping after frame render")
        return

    backward_mp4 = os.path.join(OUTPUT_DIR, "fig12_dt_video.mp4")
    forward_mp4 = os.path.join(OUTPUT_DIR, "fig12_dt_video_forward.mp4")
    glob = os.path.join(FRAME_DIR, "frame_%04d.png")
    common = ["-y", "-framerate", str(FRAMERATE), "-i", glob,
              "-c:v", "libx264", "-pix_fmt", "yuv420p",
              "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2"]
    # frame_0000 = present day, so direct frame-order stitch plays
    # present -> past (geological rewind = `<base>.mp4`); the
    # `<base>_forward.mp4` is the same content reversed so it plays
    # past -> present (chronological forward).
    print("  stitching geological backward / rewind MP4 ...")
    subprocess.run(["ffmpeg", *common, backward_mp4], check=True)
    print("  stitching chronological forward MP4 ...")
    subprocess.run(
        ["ffmpeg", *common[:5],
         "-vf", "reverse,scale=trunc(iw/2)*2:trunc(ih/2)*2",
         "-c:v", "libx264", "-pix_fmt", "yuv420p",
         forward_mp4],
        check=True,
    )
    print(f"  wrote {backward_mp4}")
    print(f"  wrote {forward_mp4}")

    if not args.keep_frames:
        for fn in os.listdir(FRAME_DIR):
            os.remove(os.path.join(FRAME_DIR, fn))


if __name__ == "__main__":
    main()
