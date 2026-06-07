#!/usr/bin/env python3
"""
Figure 2 video - 1 Myr Lambert-azimuthal-equal-area animation of
paleobathymetry, centred in the mid-North Atlantic and cropped to a
HORIZON_DEG circular cap that covers western Tethys + North Atlantic
+ Arctic Ocean down to roughly the equator.  Per-Myr companion to
the static `02_north_atlantic_paleobathymetry.py` figure (the script
filename is historic; the projection is no longer strictly polar).
Same rendering pipeline (`grdgradient` hillshading + plate-model
overlay) as the deprecated global Winkel-Tripel video at
`02a_global_paleobathymetry_video.py`.

Plate model + anchor plate are selected by config.py via
PALEO_BATHY_SOURCE:

    PALEO_BATHY_SOURCE = "Z22_mantle" -> Zahirovic 2022, anchor 0
                                          (mantle frame)
    PALEO_BATHY_SOURCE = "Z22_PMag"   -> Zahirovic 2022, anchor 701701
                                          (paleomagnetic frame)

A startup banner echoes the resolved source / dir / plate model / anchor.

Outputs (source-aware so Z22_mantle and Z22_PMag videos coexist;
distinct filename prefix from the canonical Fig 2 global video so this
script never clobbers it):
    figures/output/north_atlantic_paleobathy/frames_<SOURCE>/frame_NNNN.png
    figures/output/fig02_paleobathymetry_video_<SOURCE>.mp4
    figures/output/fig02_paleobathymetry_video_<SOURCE>_forward.mp4

Usage:
    python3 02_north_atlantic_paleobathymetry_video.py --source Z22_mantle
    python3 02_north_atlantic_paleobathymetry_video.py --source Z22_PMag

    # Use the PALEO_BATHY_SOURCE set in config.py
    python3 02_north_atlantic_paleobathymetry_video.py

    # Quick check: render every 5 Myr up to 100 Ma
    python3 02_north_atlantic_paleobathymetry_video.py --time-step 5 --max-age 100
"""
import argparse
import os
import shutil
import subprocess
import sys


# --------------------------------------------------------------------------
# CLI parse FIRST -- before importing config -- so we can set the env var
# that config.py reads.
# --------------------------------------------------------------------------
def _parse_args():
    p = argparse.ArgumentParser(
        description=("Render a 1 Myr north-pole-centred polar-"
                     "stereographic paleobathymetry animation, then "
                     "stitch frames into MP4."),
    )
    p.add_argument(
        "--source",
        choices=("Z22_mantle", "Z22_PMag"),
        default=None,
        help=("Paleobathymetry grid source.  Sets the env var "
              "PYBACKTRACK_PALEO_BATHY_SOURCE before config.py is "
              "imported, which selects both the grid directory AND the "
              "plate model + anchor (Zahirovic2022 model with anchor 0 "
              "for the Z22_mantle source; same model with anchor 701701 "
              "for Z22_PMag).  Leave unset to use the default in "
              "config.py."),
    )
    p.add_argument(
        "--max-age", type=int, default=150,
        help="Oldest frame in Ma (default: 170, matches Z22 grid range).",
    )
    p.add_argument(
        "--time-step", type=int, default=1,
        help="Frame spacing in Myr (default: 1).",
    )
    p.add_argument(
        "--framerate", type=int, default=8,
        help="Output MP4 framerate (default: 8).",
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
if _args.source is not None:
    os.environ["PYBACKTRACK_PALEO_BATHY_SOURCE"] = _args.source

# --------------------------------------------------------------------------
# NOW import config -- whatever source we resolved above will be honoured.
# --------------------------------------------------------------------------
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

# Lambert Azimuthal Equal-Area projection (A<lon>/<lat>/<horizon>/<width>)
# centred in the mid-North Atlantic.  Equal-area + bounded by the
# HORIZON_DEG cap, so visible cap covers: western Tethys
# (Mediterranean + Caribbean realm), North Atlantic, Arctic Ocean,
# down to roughly the equator, plus a slice of N Pacific / N Asia at
# the western and eastern limbs.  Tunables match the static figure
# `02_north_atlantic_paleobathymetry.py`.
CENTRE_LON = 0.0
CENTRE_LAT = 45.0
HORIZON_DEG = 50.0
PROJ_WIDTH_CM = 14

# --------------------------------------------------------------------------
# Run-time constants derived from CLI args + config.
# --------------------------------------------------------------------------
SOURCE = _cfg.PALEO_BATHY_SOURCE
PLATE_MODEL_NAME = _cfg.PLATE_MODEL_NAME
PLATE_MODEL_ANCHOR = _cfg.PLATE_MODEL_ANCHOR_PLATE
FRAME_LABEL = {0: "mantle", 701701: "paleomagnetic"}.get(
    PLATE_MODEL_ANCHOR, f"anchor={PLATE_MODEL_ANCHOR}")

TIME_STEP = _args.time_step
MAX_AGE_MA = _args.max_age
FRAMERATE = _args.framerate

# Source-aware paths so the Z22_mantle and Z22_PMag runs never clobber
# each other's frame caches.  Distinct top-level dir from the south-pole
# video (`southern_ocean_paleobathy/`) so the two views coexist.
OUT_BASE = os.path.join(OUTPUT_DIR, "north_atlantic_paleobathy")
FRAME_DIR = os.path.join(OUT_BASE, f"frames_{SOURCE}")
os.makedirs(FRAME_DIR, exist_ok=True)


def _banner():
    bar = "=" * 78
    print(bar)
    print(" NORTH-POLE PALEOBATHYMETRY VIDEO -- run configuration")
    print(bar)
    print(f"   PALEO_BATHY_SOURCE       = {SOURCE}")
    print(f"   PALEO_BATHY_DIR          = {_cfg.PALEO_BATHY_DIR}")
    print(f"   PLATE_MODEL_NAME         = {PLATE_MODEL_NAME}")
    print(f"   PLATE_MODEL_ANCHOR_PLATE = {PLATE_MODEL_ANCHOR}"
          f"   ({FRAME_LABEL} reference frame)")
    print(f"   projection               = A{CENTRE_LON}/{CENTRE_LAT}/"
          f"{HORIZON_DEG}/{PROJ_WIDTH_CM}c (Lambert azimuthal "
          f"equal-area, {HORIZON_DEG:.0f} deg horizon cap)")
    print(f"   time range               = 0..{MAX_AGE_MA} Ma "
          f"in {TIME_STEP} Myr steps  ->  "
          f"{(MAX_AGE_MA // TIME_STEP) + 1} frames")
    print(f"   framerate                = {FRAMERATE} fps")
    print(f"   frame cache              = {FRAME_DIR}")
    print(f"   output MP4 (geological)  = "
          f"fig02_paleobathymetry_video_{SOURCE}.mp4")
    print(f"   output MP4 (forward)     = "
          f"fig02_paleobathymetry_video_{SOURCE}_forward.mp4")
    print(bar)
    # Bail out loudly if the (source, model, anchor) triple is inconsistent.
    expected = {
        "Z22_mantle": ("Zahirovic2022", 0),
        "Z22_PMag":   ("Zahirovic2022", 701701),
    }.get(SOURCE)
    if expected is not None:
        exp_model, exp_anchor = expected
        if (exp_model, exp_anchor) != (PLATE_MODEL_NAME, PLATE_MODEL_ANCHOR):
            sys.exit(
                f"ERROR: source {SOURCE!r} should pair with plate model "
                f"{exp_model!r} + anchor {exp_anchor}, but got "
                f"{PLATE_MODEL_NAME!r} + anchor {PLATE_MODEL_ANCHOR}.  "
                "Fix the PALEO_BATHY_SOURCE branch in config.py."
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
    # background=True clamps out-of-range cells (depths > 5000 m) to the
    # endpoint colour instead of rendering them white.  No output= here:
    # downstream grdimage uses the session-current CPT (cmap=True), and
    # adding output= breaks that hand-off (GMT default rainbow palette).
    pygmt.makecpt(cmap=BATHY_CMAP, series=BATHY_CPT_SERIES,
                  continuous=True, background=True)

    # LAEA with HORIZON_DEG cap; `region="g"` global, the horizon in
    # the projection string does the cropping.  60/30 deg graticule;
    # no +t centred title -- the age is the boxed tag in the top-left.
    fig.basemap(
        region="g",
        projection=(f"A{CENTRE_LON}/{CENTRE_LAT}/"
                    f"{HORIZON_DEG}/{PROJ_WIDTH_CM}c"),
        frame=["xa60g60", "ya30g30", "WSne"],
    )
    _draw_paleobathymetry(fig, gplot, time_ma)

    # Big age tag in the top-left corner of the panel (matches statics).
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
    (`pbr.GRAD_DIR`, shared with the south-pole video).  gradient_for()
    keys its cache only on the integer age, so stale grad files
    produced under different HILLSHADE_* params silently survive a
    plain frame-only --force.

    Set PYBT_KEEP_HILLSHADE_CACHE=1 in the environment to override
    --force on the gradient side only (frames still get re-rendered,
    but the existing per-age gradient .nc files are kept).  Useful
    when you want to redo frames against the same hillshade.
    """
    if not _args.force:
        return
    if os.environ.get("PYBT_KEEP_HILLSHADE_CACHE"):
        print("  PYBT_KEEP_HILLSHADE_CACHE set -- keeping gradient "
              "cache even with --force (frame PNGs still wiped)")
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
    print(f"\n=== Rendering {SOURCE} north-pole frames "
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
        f"fig02_paleobathymetry_video_{SOURCE}_forward.mp4",
    )
    backward = os.path.join(
        OUTPUT_DIR,
        f"fig02_paleobathymetry_video_{SOURCE}.mp4",
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
    # Lambert Azimuthal Equal-Area centred at lon=CENTRE_LON -- pass
    # it to gplately so coastlines get split at the antimeridian
    # and don't wrap as fans across the visible cap.
    gplot = make_gplately_plot(central_meridian=CENTRE_LON)
    rendered = render_all_frames(gplot)
    if not _args.skip_video:
        stitch_videos(rendered)


if __name__ == "__main__":
    main()
