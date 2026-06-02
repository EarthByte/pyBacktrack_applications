#!/usr/bin/env bash
# ============================================================================
# run_all_figures.sh
#
# One-shot batch runner that builds every figure for the pyBacktrack 1.5
# Geoscientific Model Development paper.  It will:
#
#   1) ensure pyBacktrack, pyGMT, xarray, gplately and plate-model-manager
#      are installed (no-op if they already are);
#   2) install the pyBacktrack example data (used by the well files);
#   3) run each figure script in turn, writing PNGs to ../output/.
#
# Usage:
#   ./run_all_figures.sh             # run everything
#   ./run_all_figures.sh --skip-install   # skip the pip install steps
#   ./run_all_figures.sh fig01 fig03       # run only specific figures
#   ./run_all_figures.sh --clean fig02 fig02_video fig03 fig03_video
#                                     # wipe bathymetry outputs FIRST, then
#                                     # rebuild the listed targets from
#                                     # scratch.  --clean delegates to
#                                     # ./clean_bathy.sh and only fires when
#                                     # at least one bathymetry target
#                                     # (fig02, fig03, fig02_video*, fig03_video*)
#                                     # is requested; for non-bathy targets it
#                                     # is silently a no-op.
#
# All paths inside the scripts are relative to this directory; you can move
# the whole repo as long as ../../ still points to the paper root containing
# Zahirovic2022_pybacktrack_merged_paleobathymetry/ and NWSHELF/.
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PYTHON=${PYTHON:-python3}
SKIP_INSTALL=0
CLEAN_BATHY=0
TARGETS=()

for arg in "$@"; do
  case "$arg" in
    --skip-install) SKIP_INSTALL=1 ;;
    --clean)        CLEAN_BATHY=1 ;;
    fig01|fig02|fig02_video|fig03|fig03_video|fig03_video_z22|fig04|fig05|fig06|fig07|fig08|fig09|fig10|fig11|fig12|fig10_video|fig12_video) TARGETS+=("$arg") ;;
    *) echo "Unknown argument: $arg" >&2; exit 2 ;;
  esac
done
if [ ${#TARGETS[@]} -eq 0 ]; then
  # Default: build figures 1-12.  Fig 02 is now the LAEA-on-mid-N-Atlantic
  # 4-panel; the old Winkel-Tripel global script is reachable via
  # `fig02_global_deprecated` for comparison only.  Animations
  # (fig*_video) remain opt-in.
  TARGETS=(fig01 fig02 fig03 fig04 fig05 fig06 fig07 fig08 fig09 fig10 fig11 fig12)
fi

# ---------------------------------------------------------------------------
# 0.  Optional: wipe stale bathymetry outputs before re-rendering.
# ---------------------------------------------------------------------------
# Only triggers if --clean was passed AND at least one of the requested
# TARGETS is a bathymetry target.  Avoids accidentally nuking NW-Shelf
# outputs when the user runs --clean fig06 (which is unrelated to
# clean_bathy.sh's scope).
if [ "$CLEAN_BATHY" -eq 1 ]; then
  CLEAN_TRIGGERED=0
  for tgt in "${TARGETS[@]}"; do
    case "$tgt" in
      fig02|fig02_video|fig02_video_z22|fig03|fig03_video|fig03_video_z22)
        CLEAN_TRIGGERED=1
        break
        ;;
    esac
  done
  if [ "$CLEAN_TRIGGERED" -eq 1 ]; then
    echo "==> --clean: wiping stale bathymetry outputs first"
    "$SCRIPT_DIR/clean_bathy.sh"
    echo
  else
    echo "==> --clean: requested but no bathymetry targets in list; skipping"
  fi
fi

# ---------------------------------------------------------------------------
# 1.  Dependencies
# ---------------------------------------------------------------------------
if [ "$SKIP_INSTALL" -eq 0 ]; then
  echo "==> Installing / refreshing Python dependencies (pyBacktrack 1.5, pyGMT, ...)"
  # --break-system-packages is needed on some macOS / Linux distros that
  # restrict pip writes outside venvs.  Remove the flag if you are using a
  # virtualenv or conda environment.
  $PYTHON -m pip install --upgrade pip --break-system-packages 2>/dev/null || \
      $PYTHON -m pip install --upgrade pip
  $PYTHON -m pip install --upgrade --break-system-packages \
      "pybacktrack>=1.5" pygmt xarray netCDF4 \
      gplately plate-model-manager \
      matplotlib cartopy \
      2>/dev/null || \
      $PYTHON -m pip install --upgrade \
          "pybacktrack>=1.5" pygmt xarray netCDF4 \
          gplately plate-model-manager \
          matplotlib cartopy

  echo "==> Installing pyBacktrack example data (idempotent)"
  $PYTHON -m pybacktrack.install_examples -- examples_data || true
fi

# ---------------------------------------------------------------------------
# 2.  Run each figure script
# ---------------------------------------------------------------------------
# Filter pattern for harmless macOS stderr noise from pyGMT's GMT
# subprocesses (and Qt5/Qt6 mixed loads when running in conda base).
# Real errors are still passed through; only matching lines are silenced.
NOISE_REGEX='MallocStackLogging|MallocStackLoggingNoCompact'

run_step() {
  local name="$1"
  local script="$2"
  echo
  echo "========================================================================"
  echo " Running ${name}: ${script}"
  echo "========================================================================"
  # 2> >(grep ...) pipes stderr through a filter; everything that doesn't
  # match NOISE_REGEX is forwarded back to stderr (>&2).
  time $PYTHON "$script" 2> >(grep -E -v "$NOISE_REGEX" >&2)
}

for tgt in "${TARGETS[@]}"; do
  case "$tgt" in
    # Script numbers now align with figure numbers (2026-05-28 rename).
    fig01) run_step "Fig 01 - Workflow flowchart"                              "01_workflow_flowchart.py" ;;
    fig02)
      run_step "Fig 02 - LAEA paleobathymetry, mid-N-Atlantic centred (4-panel)" "02_north_atlantic_paleobathymetry.py"
      ;;
    fig02_video)
      run_step "Fig 02 (video) - LAEA paleobathymetry 1 Myr animation"           "02_north_atlantic_paleobathymetry_video.py"
      ;;
    fig03) run_step "Fig 03 - LAEA south-pole paleobathymetry (4-panel)"       "03_southern_ocean_paleobathymetry.py" ;;
    fig04)
      # Frame-invariant comparison of pyBacktrack 1.5 vs
      # Straume 2020 forward-model paleobathymetry: median + MAD.
      run_step "Fig 04 - PyBacktrack vs Straume paleobathymetry comparison"          "04_pybacktrack_vs_straume_comparison.py"
      ;;
    fig05) run_step "Fig 05 - NW Shelf map (VGG basemap, asteras highlighted)" "05_nwshelf_map.py" ;;
    fig06) run_step "Fig 06 - SL + DT time series (dual y-axis) at asteras"    "06_sl_dt_curves.py" ;;
    fig07) run_step "Fig 07 - Asteras backstripping (3 configs in one panel)"  "07_backstrip_wells.py" ;;
    fig08) run_step "Fig 08 - Asteras geohistory (Wheeler diagram)"            "08_geohistory_wells.py" ;;
    fig09)
      # Margin-scale: data step (slow!) then rate-time-series figure.
      run_step "Fig 09 (data step)    - Backstrip all NW Shelf wells + grid rates" "07a_backstrip_all_nwshelf.py"
      run_step "Fig 09 (time-series)  - Mean +/- 1 sigma rate, A / C / A-C"        "09_rate_timeseries.py"
      ;;
    fig10)
      # Margin-scale combined 3x4 rate-map matrix (4 max-delta times).
      # Requires 07a_backstrip_all_nwshelf.py to have been run first
      # (fig09 takes care of that).
      run_step "Fig 10 - NW Shelf rate maps, 3 cols x 4 times (combined)"        "10_rate_maps.py"
      ;;
    fig11)
      # Margin-scale box plots of the spatial rate distribution.
      run_step "Fig 11 - NW Shelf rate box plots (10 Myr intervals)"             "11_rate_boxplots.py"
      ;;
    fig12)
      # Raw dynamic-topography 2x2 panel at the same 4 times as Fig 10.
      # Depends on the picked-times sidecar written by 10_rate_maps.py
      # (run fig10 first) and on the well-locations sidecar written by
      # 07a_backstrip_all_nwshelf.py (run fig09 first).
      run_step "Fig 12 - Raw dynamic-topography 2x2 panel"                       "12_dt_maps.py"
      ;;
    fig03_video)
      run_step "Fig 03 (video, default source) - South-pole paleobathymetry 1 Myr animation (0-170 Ma + ffmpeg)" "03a_southern_ocean_paleobathymetry_video.py"
      ;;
    fig03_video_z22)
      run_step "Fig 03 (video, Z22 hardwired) - South-pole paleobathymetry 1 Myr animation in Z22 paleomag frame" "03b_southern_ocean_paleobathymetry_video_Z22.py"
      ;;
    fig10_video)
      # 1-Myr subsidence-rate animations (companion to Fig 10 maps).
      run_step "Fig 10 (video) - Subsidence-rate animations (1 Myr frames + ffmpeg)" "10a_make_subsidence_videos.py"
      ;;
    fig12_video)
      run_step "Fig 12 (video) - Raw dynamic-topography 1 Myr animation"             "12a_dt_video.py"
      ;;
  esac
done

echo
echo "All figures written to: $SCRIPT_DIR/../output/"
ls -1 "$SCRIPT_DIR/../output/" || true
