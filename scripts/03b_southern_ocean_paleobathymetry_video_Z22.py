#!/usr/bin/env python3
"""
Figure 3 video (Z22) - 1 Myr south-pole orthographic animation of
global paleobathymetry, hardwired to the Zahirovic et al. (2022) GDH1
grids in the paleomagnetic reference frame.

Hardwired:
    PALEO_BATHY_SOURCE       = "Z22_PMag"
    PLATE_MODEL_NAME         = "Zahirovic2022"
    PLATE_MODEL_ANCHOR_PLATE = 701701       (paleomagnetic frame)
    PALEO_BATHY_DIR          = /Users/dietmar/Documents/GPlates/
                               pyBacktrack1.5/Paleobathymetry_Z22_PMag
                               (fallback: Zahirovic2022_paleobathymetry_grids/
                                under the paper root)

The script sets PYBACKTRACK_PALEO_BATHY_SOURCE before importing
config so the right paths, plate model and anchor are resolved.
No CLI override -- run `03a_southern_ocean_paleobathymetry_video.py` for
the Z22 mantle-frame variant.  A startup banner echoes the resolved
source / dir / plate model / anchor.

Outputs:
    figures/output/southpole_paleobathy/frames_Z22_PMag/frame_NNNN.png
    figures/output/fig03_southern_ocean_paleobathymetry_video_Z22_PMag.mp4
        ^ geological-time backward (MAX_AGE_MA -> 0)
    figures/output/fig03_southern_ocean_paleobathymetry_video_Z22_PMag_forward.mp4
        ^ chronological forward (0 -> MAX_AGE_MA)

Usage:
    python3 03b_southern_ocean_paleobathymetry_video_Z22.py              # reuse cached frames if present
    python3 03b_southern_ocean_paleobathymetry_video_Z22.py --force      # wipe cache and re-render every frame
    python3 03b_southern_ocean_paleobathymetry_video_Z22.py --skip-video # render frames but don't run ffmpeg
"""
import argparse
import os
import shutil
import subprocess
import sys


def _parse_args():
    p = argparse.ArgumentParser(
        description=("Z22 paleomagnetic-frame south-pole paleobathymetry "
                     "animation."),
    )
    p.add_argument(
        "--force", action="store_true",
        help="Re-render frames even if they already exist on disk.",
    )
    p.add_argument(
        "--skip-video", action="store_true",
        help="Render frames but don't run ffmpeg.",
    )
    return p.parse_args()


_args = _parse_args()

# ----------------------------------------------------------------------------
# HARDWIRED: set the env var BEFORE importing config so config.py's
# PALEO_BATHY_SOURCE branch picks Z22 (and therefore Zahirovic2022 plate
# model + anchor 701701).
# ----------------------------------------------------------------------------
os.environ["PYBACKTRACK_PALEO_BATHY_SOURCE"] = "Z22_PMag"

import pygmt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as _cfg
from config import OUTPUT_DIR, PALEO_BATHY_FMT          # noqa: E402

# Shared rendering helpers (gplately driver, hillshade pipeline,
# bathymetry draw + plate-model overlay).
import paleobathy_render as pbr

make_gplately_plot   = pbr.make_gplately_plot
_draw_paleobathymetry = pbr.draw_paleobathymetry
BATHY_CMAP           = pbr.BATHY_CMAP
BATHY_CPT_SERIES     = pbr.BATHY_CPT_SERIES
NAN_COLOR            = pbr.NAN_COLOR
FONT_TITLE           = pbr.FONT_TITLE
FONT_LABEL           = pbr.FONT_LABEL
FONT_ANNOT           = pbr.FONT_ANNOT

# ----------------------------------------------------------------------------
# Hardwired params (no CLI).
# ----------------------------------------------------------------------------
SOURCE = "Z22_PMag"                          # for output filenames + banner
EXPECTED_PLATE_MODEL = "Zahirovic2022"       # belt-and-braces sanity check
EXPECTED_ANCHOR = 701701                     # belt-and-braces sanity check
MAX_AGE_MA = 170                             # match the Z22 grid range
TIME_STEP = 1                                # Myr between frames
FRAMERATE = 8                                # MP4 fps

# LAEA south-pole, HORIZON_DEG=50 -> cap edge at -40 S.  Matches the
# static Fig 3 projection family.
SOUTHPOLE_LON = 0.0
SOUTHPOLE_LAT = -90.0
HORIZON_DEG   = 50.0
PROJ_WIDTH_CM = 14

# ----------------------------------------------------------------------------
# Resolved from config (post env-var override) and used for output paths.
# ----------------------------------------------------------------------------
PLATE_MODEL_NAME = _cfg.PLATE_MODEL_NAME
PLATE_MODEL_ANCHOR = _cfg.PLATE_MODEL_ANCHOR_PLATE
FRAME_LABEL = {0: "mantle", 701701: "paleomagnetic"}.get(
    PLATE_MODEL_ANCHOR, f"anchor={PLATE_MODEL_ANCHOR}")

OUT_BASE = os.path.join(OUTPUT_DIR, "southpole_paleobathy")
FRAME_DIR = os.path.join(OUT_BASE, f"frames_{SOURCE}")
os.makedirs(FRAME_DIR, exist_ok=True)


def _banner():
    bar = "=" * 78
    print(bar)
    print(" SOUTH-POLE PALEOBATHYMETRY VIDEO -- Zahirovic 2022 (hardwired)")
    print(bar)
    print(f"   PALEO_BATHY_SOURCE       = {_cfg.PALEO_BATHY_SOURCE}")
    print(f"   PALEO_BATHY_DIR          = {_cfg.PALEO_BATHY_DIR}")
    print(f"   PLATE_MODEL_NAME         = {PLATE_MODEL_NAME}")
    print(f"   PLATE_MODEL_ANCHOR_PLATE = {PLATE_MODEL_ANCHOR}"
          f"   ({FRAME_LABEL} reference frame)")
    print(f"   projection               = A{SOUTHPOLE_LON}/{SOUTHPOLE_LAT}/"
          f"{HORIZON_DEG}/{PROJ_WIDTH_CM}c (Lambert azimuthal "
          f"equal-area, {HORIZON_DEG:.0f} deg south-pole cap)")
    print(f"   time range               = 0..{MAX_AGE_MA} Ma "
          f"in {TIME_STEP} Myr steps  ->  "
          f"{(MAX_AGE_MA // TIME_STEP) + 1} frames")
    print(f"   framerate                = {FRAMERATE} fps")
    print(f"   frame cache              = {FRAME_DIR}")
    print(f"   output MP4 (geological)  = "
          f"fig03_southern_ocean_paleobathymetry_video_{SOURCE}.mp4")
    print(f"   output MP4 (forward)     = "
          f"fig03_southern_ocean_paleobathymetry_video_{SOURCE}_forward.mp4")
    print(bar)
    # Hard-stop if config didn't actually land on Z22 + Zahirovic2022 +
    # anchor 701701 -- only fires if someone hand-edited the env-var
    # handling or the Z22_PMag branch in config.py.
    if _cfg.PALEO_BATHY_SOURCE != SOURCE:
        sys.exit(
            f"ERROR: env-var override didn't stick -- expected "
            f"PALEO_BATHY_SOURCE={SOURCE!r}, got "
            f"{_cfg.PALEO_BATHY_SOURCE!r}"
        )
    if PLATE_MODEL_NAME != EXPECTED_PLATE_MODEL:
        sys.exit(
            f"ERROR: plate model {PLATE_MODEL_NAME!r} does not match "
            f"the expected {EXPECTED_PLATE_MODEL!r} for source {SOURCE!r}."
        )
    if PLATE_MODEL_ANCHOR != EXPECTED_ANCHOR:
        sys.exit(
            f"ERROR: anchor plate {PLATE_MODEL_ANCHOR} does not match "
            f"the expected anchor {EXPECTED_ANCHOR} for source {SOURCE!r}."
        )


def _frame_path(time_ma):
    return os.path.join(FRAME_DIR, f"frame_{int(time_ma):04d}.png")


def render_frame(time_ma, gplot, force=False):
    out = _frame_path(time_ma)
    if os.path.exists(out) and not force:
        return out

    grid_path = PALEO_BATHY_FMT.format(time=float(time_ma))
    if not os.path.exists(grid_path):
        print(f"  ! {time_ma} Ma: grid missing ({grid_path}), skipping")
        return None

    fig = pygmt.Figure()
    pygmt.config(
        FONT_TITLE=FONT_TITLE,
        FONT_LABEL="18p,Helvetica", FONT_ANNOT="14p,Helvetica",
        COLOR_NAN=NAN_COLOR, MAP_FRAME_TYPE="fancy",
        MAP_GRID_PEN_PRIMARY="0.6p,gray40",   # visible 60/30 deg graticule
    )
    pygmt.makecpt(cmap=BATHY_CMAP, series=BATHY_CPT_SERIES, continuous=True)

    # LAEA south-pole, HORIZON_DEG cap.  60/30 deg graticule, no centred
    # title -- the age is the boxed tag in the top-left.
    fig.basemap(
        region="g",
        projection=(f"A{SOUTHPOLE_LON}/{SOUTHPOLE_LAT}/"
                    f"{HORIZON_DEG}/{PROJ_WIDTH_CM}c"),
        frame=["xa60g60", "ya30g30", "WSne"],
    )
    _draw_paleobathymetry(fig, gplot, time_ma)

    fig.text(
        text=f"{int(time_ma)} Ma",
        position="TL", offset="0.3c/-0.3c", no_clip=True,
        font="22p,Helvetica-Bold,black",
        fill="white", pen="0.5p,black",
        clearance="0.18c/0.10c",
    )

    fig.colorbar(
        frame=['a1000f200+lPaleobathymetry (m)'],
        position="JBC+w12c/0.4c+h+o0/1.2c+e+ma",
    )

    fig.savefig(out, dpi=200)
    return out


def _wipe_gradient_cache_if_forced():
    """If --force is set, also wipe the hillshade gradient cache
    (`pbr.GRAD_DIR`).  gradient_for() keys its cache only on the
    integer age, so stale grad files produced under different
    HILLSHADE_* params silently survive a plain frame-only --force.
    """
    if not _args.force:
        return
    grad_dir = pbr.GRAD_DIR
    if not os.path.isdir(grad_dir):
        return
    wiped = 0
    for fn in os.listdir(grad_dir):
        if fn.endswith(".nc"):
            os.remove(os.path.join(grad_dir, fn))
            wiped += 1
    if wiped:
        print(f"  --force: wiped {wiped} cached gradient(s) from "
              f"{grad_dir} (will recompute with current "
              f"HILLSHADE_* parameters)")


def render_all_frames(gplot):
    print(f"\n=== Rendering {SOURCE} south-pole frames "
          f"({FRAME_LABEL} frame) ===")
    _wipe_gradient_cache_if_forced()
    all_times = list(range(0, MAX_AGE_MA + 1, TIME_STEP))
    existing = sum(1 for t in all_times if os.path.exists(_frame_path(t)))
    if _args.force:
        print(f"  --force: {existing} cached frame(s) will be overwritten; "
              f"rendering {len(all_times)} frames fresh")
    elif existing == 0:
        print(f"  cache empty: rendering {len(all_times)} frames fresh")
    elif existing == len(all_times):
        print(f"  cache hit: all {existing} frames already on disk "
              f"(re-stitch only).  Pass --force to regenerate.")
    else:
        print(f"  partial cache: {existing} frame(s) exist, "
              f"{len(all_times) - existing} will be rendered.  "
              f"Pass --force to regenerate everything.")

    rendered, regenerated = [], 0
    for t in all_times:
        already = os.path.exists(_frame_path(t))
        out = render_frame(t, gplot, force=_args.force)
        if out is not None:
            rendered.append(out)
            if not already or _args.force:
                regenerated += 1
            if t % 10 == 0:
                tag = "(cached)" if (already and not _args.force) else ""
                print(f"  {t:3d} Ma -> {os.path.basename(out)} {tag}")
    print(f"  done: {len(rendered)} frames available "
          f"({regenerated} freshly rendered, "
          f"{len(rendered) - regenerated} reused)")
    return rendered


def stitch_videos(rendered):
    if shutil.which("ffmpeg") is None:
        print("\nffmpeg not found on PATH -- skipping video step")
        return
    if not rendered:
        print("\nno frames rendered, skipping video step")
        return

    forward = os.path.join(
        OUTPUT_DIR,
        f"fig03_southern_ocean_paleobathymetry_video_{SOURCE}_forward.mp4",
    )
    backward = os.path.join(
        OUTPUT_DIR,
        f"fig03_southern_ocean_paleobathymetry_video_{SOURCE}.mp4",
    )

    # frame_0000 = present day -> direct frame-order stitch plays
    # present -> past (geological rewind = `<base>.mp4`); the
    # `<base>_forward.mp4` is obtained via `-vf reverse` so it plays
    # past -> present (chronological forward).
    cmd_backward = [
        "ffmpeg", "-y", "-framerate", str(FRAMERATE),
        "-i", os.path.join(FRAME_DIR, "frame_%04d.png"),
        "-pix_fmt", "yuv420p", "-vcodec", "libx264", "-crf", "23",
        "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
        backward,
    ]
    cmd_forward = [
        "ffmpeg", "-y", "-i", backward, "-vf", "reverse",
        "-vcodec", "libx264", "-crf", "23", "-pix_fmt", "yuv420p",
        forward,
    ]
    print(f"\n=== ffmpeg: geological backward / rewind ({backward}) ===")
    subprocess.run(cmd_backward, check=True)
    print(f"\n=== ffmpeg: chronological forward ({forward}) ===")
    subprocess.run(cmd_forward, check=True)


def main():
    _banner()
    # South-pole orthographic centred at lon=SOUTHPOLE_LON -- pass it
    # to gplately so coastlines get split at the antimeridian.
    gplot = make_gplately_plot(central_meridian=SOUTHPOLE_LON)
    rendered = render_all_frames(gplot)
    if not _args.skip_video:
        stitch_videos(rendered)


if __name__ == "__main__":
    main()
