#!/usr/bin/env python3
"""
1 Myr animations of NW Shelf subsidence-rate fields.

Per-Myr companion to Fig. 10 (the 3 x 4 rate-map matrix): renders one
single-panel frame per Myr in [0, cutoff] for each of three
configurations, and stitches the frames into MP4s via ffmpeg.

    rate_A : no SL, no DT                              (dem1 cmap,
                                                        0-50 m/Myr)
    rate_C : Haq2024 SL + D10_gmcm9 DT                 (dem1 cmap,
             (Braz et al., 2021)                        0-50 m/Myr)
    rate_D : difference  rate_A - rate_C               (vik diverging,
                                                       -50..+50 m/Myr)

Frame styling matches the per-panel look of Fig. 10 exactly so the
animations read as a continuous extension of the static figure:

    - blockmedian -> surface upsampling (0.05 deg, tension 0.5)
      via Fig. 10's `make_continuous_masked_grid`,
    - concave-hull (alpha-shape) mask of the 109 well locations,
    - white-filled circles over-plotted at every well,
    - dem1 + e clamp arrows for the absolute-rate panels,
      vik symmetric for the difference panel.

The time range is auto-clipped to the defensible-time cutoff written
by `07a_backstrip_all_nwshelf.py` (`cutoff_time_Ma.txt`).  A forward
(0 Ma -> cutoff) and a geological-time backward (cutoff -> 0 Ma) MP4
are produced for each configuration.

Outputs (in figures/output/):
    fig10_nwshelf_rate_video_A_no_sl_no_dt.mp4
    fig10_nwshelf_rate_video_C_sl_and_dt.mp4
    fig10_nwshelf_rate_video_D_difference.mp4
        (plus the matching *_forward.mp4 trio)
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import OUTPUT_DIR


def _parse_args():
    p = argparse.ArgumentParser(
        description=("Render 1 Myr NW Shelf subsidence-rate frames for "
                     "configs A, C and A - C, and stitch MP4s.  By "
                     "default any frame PNG already on disk is reused; "
                     "pass --force to regenerate everything."),
    )
    p.add_argument(
        "--force", action="store_true",
        help="Wipe the per-config frame caches and re-render every frame.",
    )
    p.add_argument(
        "--skip-video", action="store_true",
        help="Render frames but don't run ffmpeg.",
    )
    return p.parse_args()


_args = _parse_args()


# Reuse Fig 10 helpers via importlib so the per-frame look matches
# exactly (concave-hull mask, blockmedian->surface upsampling, CPT
# builder, cached well loader).
_spec = importlib.util.spec_from_file_location(
    "_fig10",
    os.path.join(os.path.dirname(__file__), "10_rate_maps.py"),
)
_fig10 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_fig10)

REGION = _fig10.REGION
RATE_CMAP = _fig10.RATE_CMAP
RATE_SERIES = _fig10.RATE_SERIES
DELTA_CMAP = _fig10.DELTA_CMAP
DELTA_SERIES = _fig10.DELTA_SERIES
PANEL_TITLE_FONT = _fig10.PANEL_TITLE_FONT
_save_cpt = _fig10._save_cpt
make_continuous_masked_grid = _fig10.make_continuous_masked_grid
_well_locations = _fig10._well_locations
_grid_is_drawable = _fig10._grid_is_drawable


# ----------------------------------------------------------------------------
# Tunables
# ----------------------------------------------------------------------------
OUT_BASE = os.path.join(OUTPUT_DIR, "nwshelf_subsidence")
TIME_STEP = 1
FRAMERATE = 6
FRAME_WIDTH = "12c"             # single Mercator NW Shelf panel
CBAR_WIDTH = "11c"
CBAR_HEIGHT = "0.4c"
CBAR_Y_OFFSET = "1.2c"

# (subdir, cfg_label, is_delta, out_basename) -- the cmap/range comes
# from RATE_* / DELTA_* (imported above) so the videos stay in lockstep
# with Fig 9's static panels.
CONFIGS = [
    ("A_no_sl_no_dt_rate", "A: no SL, no DT",                 False,
     "fig10_nwshelf_rate_video_A_no_sl_no_dt"),
    ("C_sl_and_dt_rate",   "C: Haq2024 SL + D10_gmcm9 DT",    False,
     "fig10_nwshelf_rate_video_C_sl_and_dt"),
    ("D_rate_difference",  "D = A - C",                        True,
     "fig10_nwshelf_rate_video_D_difference"),
]


def read_cutoff():
    """Return the defensible-time cutoff written by 07a_backstrip_all_nwshelf.py."""
    p = os.path.join(OUT_BASE, "cutoff_time_Ma.txt")
    if os.path.exists(p):
        try:
            return int(open(p).read().strip())
        except Exception:
            pass
    return 150


# ----------------------------------------------------------------------------
# Per-frame renderer
# ----------------------------------------------------------------------------
def render_frame(subdir, cfg_label, t, cpt_path, cbar_ticks, cbar_label,
                 well_lons, well_lats, frame_dir, force=False):
    """Render a single-panel rate-map frame matching Fig 9's per-panel look.

    If the PNG already exists on disk and `force` is False the frame is
    reused (the path is returned with a flag).  With `force=True` the
    frame is always rendered fresh.

    Returns (path, was_freshly_rendered) or (None, False) if the
    underlying grid is missing.
    """
    grid_path = os.path.join(OUT_BASE, subdir, f"rate_{int(t)}.nc")
    if not _grid_is_drawable(grid_path):
        return None, False
    fname = os.path.join(frame_dir, f"frame_{int(t):04d}.png")
    if os.path.exists(fname) and not force:
        return fname, False
    fig = pygmt.Figure()
    pygmt.config(
        FONT_TITLE="14p,Helvetica-Bold",
        FONT_LABEL="13p",
        FONT_ANNOT="10p",
        MAP_FRAME_TYPE="plain",
        COLOR_NAN="240/240/240",
    )
    fig.basemap(
        region=REGION, projection=f"M{FRAME_WIDTH}",
        frame=["WSne", "xa5f1+lLongitude", "ya5f1+lLatitude"],
    )
    grid = make_continuous_masked_grid(grid_path, REGION)
    if grid is not None:
        fig.grdimage(
            grid=grid, region=REGION,
            cmap=cpt_path, nan_transparent=True,
        )
    fig.coast(shorelines="0.3p,black", resolution="c")
    fig.plot(
        x=well_lons, y=well_lats,
        style="c0.10c", fill="white", pen="0.25p,black",
    )
    fig.text(
        x=REGION[0] + 0.4, y=REGION[3] - 0.4,
        text=f"{cfg_label}  -  {int(t)} Ma",
        font=PANEL_TITLE_FONT, justify="TL",
        fill="white@30", pen="0.25p,black",
    )
    fig.colorbar(
        cmap=cpt_path,
        position=(f"JBC+w{CBAR_WIDTH}/{CBAR_HEIGHT}+h"
                  f"+o0c/{CBAR_Y_OFFSET}+ma+e"),
        frame=[cbar_ticks, f"x+l{cbar_label}"],
    )
    fig.savefig(fname, dpi=150)
    return fname, True


# ----------------------------------------------------------------------------
# Per-configuration driver
# ----------------------------------------------------------------------------
def make_video_for_config(subdir, cfg_label, is_delta, out_basename,
                          cutoff, well_lons, well_lats, cpt_path,
                          force=False):
    """Render all frames for one configuration, then stitch MP4s.

    Frames already present on disk are reused unless `force=True`,
    in which case the per-config frame cache is wiped first and every
    frame is rendered from scratch.
    """
    if is_delta:
        cbar_ticks = "a25f5"
        cbar_label = "Rate difference A - C (m/Myr)"
    else:
        cbar_ticks = "a10f5"
        cbar_label = "Subsidence rate (m/Myr)"

    frame_dir = os.path.join(OUT_BASE, f"frames_{subdir}")
    os.makedirs(frame_dir, exist_ok=True)
    if force:
        wiped = 0
        for fn in os.listdir(frame_dir):
            if fn.endswith(".png"):
                os.remove(os.path.join(frame_dir, fn))
                wiped += 1
        if wiped:
            print(f"  --force: wiped {wiped} cached frame(s) from "
                  f"{frame_dir}")

    times_to_render = list(np.arange(0, cutoff + 1, TIME_STEP))
    existing = sum(
        1 for t in times_to_render
        if os.path.exists(os.path.join(frame_dir, f"frame_{int(t):04d}.png"))
    )
    if force:
        print(f"  --force: rendering {len(times_to_render)} frames fresh")
    elif existing == 0:
        print(f"  cache empty: rendering {len(times_to_render)} frames fresh")
    elif existing == len(times_to_render):
        print(f"  cache hit: all {existing} frames already on disk "
              f"(re-stitch only).  Pass --force to regenerate.")
    else:
        print(f"  partial cache: {existing} frame(s) exist, "
              f"{len(times_to_render) - existing} will be rendered.  "
              f"Pass --force to regenerate everything.")

    print(f"\n=== Rendering frames for {cfg_label} (0-{cutoff} Ma) ===")
    overall_t0 = time.time()
    rendered, regenerated = [], 0
    for t in times_to_render:
        result, was_fresh = render_frame(
            subdir, cfg_label, t, cpt_path, cbar_ticks, cbar_label,
            well_lons, well_lats, frame_dir, force=force,
        )
        if result is not None:
            rendered.append(result)
            if was_fresh:
                regenerated += 1
            if len(rendered) % 25 == 0:
                tag = "(cached)" if not was_fresh else ""
                print(f"  [{len(rendered):3d}] t={int(t):>3} Ma "
                      f"-> {os.path.basename(result)} {tag}")
    print(f"  {len(rendered)} frames available "
          f"({regenerated} freshly rendered, "
          f"{len(rendered) - regenerated} reused) in "
          f"{time.time() - overall_t0:.1f} s")
    if not rendered:
        print("  ! no frames rendered; skipping MP4 stitch")
        return
    if _args.skip_video:
        print("  --skip-video: stopping after frame render")
        return

    forward = os.path.join(OUTPUT_DIR, f"{out_basename}_forward.mp4")
    backward = os.path.join(OUTPUT_DIR, f"{out_basename}.mp4")

    # frame_0000 = present day so the direct frame-order stitch is
    # the geological rewind (present -> past) and `-vf reverse` is
    # the chronological forward (past -> present).
    cmd_backward = [
        "ffmpeg", "-y", "-framerate", str(FRAMERATE),
        "-i", os.path.join(frame_dir, "frame_%04d.png"),
        "-pix_fmt", "yuv420p", "-vcodec", "libx264", "-crf", "23",
        "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
        backward,
    ]
    cmd_forward = [
        "ffmpeg", "-y", "-i", backward,
        "-vf", "reverse,scale=trunc(iw/2)*2:trunc(ih/2)*2",
        "-vcodec", "libx264", "-crf", "23", "-pix_fmt", "yuv420p",
        forward,
    ]
    print(f"  stitching geological backward (rewind) -> {backward}")
    subprocess.run(cmd_backward, check=True)
    print(f"  stitching chronological forward        -> {forward}")
    subprocess.run(cmd_forward, check=True)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    if shutil.which("ffmpeg") is None:
        sys.exit("ffmpeg not found on PATH -- install it "
                 "(brew install ffmpeg or apt install ffmpeg)")

    cutoff = read_cutoff()
    wells = _well_locations()
    well_lons = wells[:, 0]
    well_lats = wells[:, 1]

    # Pre-build the two CPTs once (the same hang-fix as Fig 9: avoid
    # makecpt-inside-loop on macOS Apple-Silicon).
    cpt_dir = os.path.join(OUT_BASE, "cpts")
    os.makedirs(cpt_dir, exist_ok=True)
    rate_cpt = os.path.join(cpt_dir, "fig10_rate.cpt")
    delta_cpt = os.path.join(cpt_dir, "fig10_delta.cpt")
    print(f"  building {RATE_CMAP} CPT {RATE_SERIES} -> {rate_cpt}")
    _save_cpt(RATE_CMAP, RATE_SERIES, rate_cpt)
    print(f"  building {DELTA_CMAP} CPT {DELTA_SERIES} -> {delta_cpt}")
    _save_cpt(DELTA_CMAP, DELTA_SERIES, delta_cpt)

    for subdir, cfg_label, is_delta, out_basename in CONFIGS:
        cpt_path = delta_cpt if is_delta else rate_cpt
        make_video_for_config(
            subdir, cfg_label, is_delta, out_basename, cutoff,
            well_lons, well_lats, cpt_path, force=_args.force,
        )


if __name__ == "__main__":
    main()
