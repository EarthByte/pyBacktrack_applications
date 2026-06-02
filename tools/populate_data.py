#!/usr/bin/env python3
"""
One-shot data-population helper for pyBacktrack_applications/.

Reads from the *originating* on-disk locations (Sandwell VGG on disk,
D10_gmcm9 grids bundled with pyBacktrack, NW Shelf wells in the
pyBacktrack1.5 paper repo) and writes region-cropped, repo-local copies
under ``../data/``.  Run this once after cloning the repo on a machine
that already has the originating files; from then on every figure
script reads exclusively from ``data/`` and the repo is self-contained.

What it does
------------
1. ``data/wells/`` -- copy 109 NW Shelf wells (asteras.txt, plover.txt, ...)
   from ``pyBacktrack1.5/NWSHELF/backstripping_example/wells/``.
2. ``data/grids/vgg_nwshelf.nc`` -- grdcut the global Sandwell V31.1 VGG
   grid to the Fig 5 region [114, 130, -20, -9].
3. ``data/grids/dynamic_topography/D10_gmcm9/<time>.nc`` -- grdcut every
   per-time slice of the bundled D10_gmcm9 dynamic-topography model
   (Braz et al. 2021) to the NW Shelf region [113, 132, -22, -8] +
   buffer.

Usage
-----
    cd pyBacktrack_applications/tools
    python populate_data.py [--vgg PATH] [--dt-bundle PATH] [--wells PATH]

By default it expects the source files at:
    --vgg          /Users/dietmar/grids/vgg_31.1.nc
    --dt-bundle    <pybacktrack>/bundle_data/dynamic_topography/models/Braz2021/D10_gmcm9/
    --wells        ../../NWSHELF/backstripping_example/wells/

Each flag can also be set via the environment variables PYBT_VGG_SOURCE,
PYBT_DT_SOURCE and PYBT_WELLS_SOURCE.
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Region definitions — keep in sync with the consumers!
#   - VGG     : Fig 5's REGION (05_nwshelf_map.py)
#   - DT      : Fig 12's REGION + 4 deg buffer (12_dt_maps.py)
# ---------------------------------------------------------------------------
VGG_REGION = (114, 130, -20, -9)
DT_REGION = (109, 136, -26, -4)   # Fig 12's [113,132,-22,-8] + 4 deg buffer


# ---------------------------------------------------------------------------
def _here() -> Path:
    return Path(__file__).resolve().parent


def _repo_root() -> Path:
    return _here().parent


def _default_vgg() -> Path:
    return Path(os.environ.get("PYBT_VGG_SOURCE",
                               "/Users/dietmar/grids/vgg_31.1.nc"))


def _default_dt_bundle() -> Path | None:
    """Resolve the bundled D10_gmcm9 dir from the installed pybacktrack."""
    if "PYBT_DT_SOURCE" in os.environ:
        return Path(os.environ["PYBT_DT_SOURCE"])
    try:
        import pybacktrack
    except ImportError:
        return None
    base = (Path(pybacktrack.__file__).resolve().parent
            / "bundle_data" / "dynamic_topography" / "models" / "Braz2021"
            / "D10_gmcm9")
    return base if base.is_dir() else None


def _default_wells() -> Path:
    return Path(os.environ.get(
        "PYBT_WELLS_SOURCE",
        str(_repo_root().parent / "NWSHELF" / "backstripping_example" / "wells")))


# ---------------------------------------------------------------------------
def _grdcut(src: Path, dst: Path, region: tuple[float, float, float, float]):
    """Wrap pygmt.grdcut with a clear error if pygmt isn't installed."""
    try:
        import pygmt
    except ImportError as exc:
        raise SystemExit(
            "pygmt is required for grdcut.  Install via "
            "`conda install -c conda-forge pygmt` and re-run."
        ) from exc
    print(f"  grdcut {src.name}  R={region} -> {dst.relative_to(_repo_root())}")
    pygmt.grdcut(grid=str(src), region=list(region), outgrid=str(dst))


def _do_wells(src_dir: Path, dst_dir: Path) -> int:
    if not src_dir.is_dir():
        print(f"[skip] wells source missing: {src_dir}")
        return 0
    dst_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    for fn in sorted(os.listdir(src_dir)):
        if not (fn.endswith(".txt") or fn.endswith(".dat")):
            continue
        shutil.copy2(src_dir / fn, dst_dir / fn)
        n += 1
    print(f"[wells] copied {n} files into {dst_dir.relative_to(_repo_root())}/")
    return n


def _do_vgg(src: Path, dst: Path) -> bool:
    if not src.exists():
        print(f"[skip] VGG source missing: {src}\n"
              "       Supply a Sandwell V31.1 .nc via --vgg or "
              "$PYBT_VGG_SOURCE.")
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    _grdcut(src, dst, VGG_REGION)
    return True


def _do_dt(src_dir: Path | None, dst_dir: Path) -> int:
    if src_dir is None or not src_dir.is_dir():
        print(f"[skip] DT bundle source missing: {src_dir}\n"
              "       Install pybacktrack (so the bundled D10_gmcm9 grids are "
              "available) or supply --dt-bundle / $PYBT_DT_SOURCE.")
        return 0
    dst_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    for fn in sorted(os.listdir(src_dir)):
        if not fn.endswith(".nc"):
            continue
        m = re.match(r"^([0-9]+\.[0-9]+)\.nc$", fn)
        if m is None:
            continue
        src = src_dir / fn
        dst = dst_dir / fn
        _grdcut(src, dst, DT_REGION)
        n += 1
    print(f"[dt]    cropped {n} D10_gmcm9 slices into "
          f"{dst_dir.relative_to(_repo_root())}/")
    return n


# ---------------------------------------------------------------------------
def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vgg",        type=Path, default=_default_vgg())
    p.add_argument("--dt-bundle",  type=Path, default=_default_dt_bundle())
    p.add_argument("--wells",      type=Path, default=_default_wells())
    p.add_argument("--skip-wells", action="store_true")
    p.add_argument("--skip-vgg",   action="store_true")
    p.add_argument("--skip-dt",    action="store_true")
    args = p.parse_args()

    repo = _repo_root()
    print(f"== pyBacktrack_applications/tools/populate_data.py")
    print(f"   repo root  : {repo}")
    print(f"   VGG source : {args.vgg}")
    print(f"   DT bundle  : {args.dt_bundle}")
    print(f"   wells source: {args.wells}")

    if not args.skip_wells:
        _do_wells(args.wells, repo / "data" / "wells")
    if not args.skip_vgg:
        _do_vgg(args.vgg, repo / "data" / "grids" / "vgg_nwshelf.nc")
    if not args.skip_dt:
        _do_dt(args.dt_bundle,
               repo / "data" / "grids" / "dynamic_topography" / "D10_gmcm9")

    print("done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
