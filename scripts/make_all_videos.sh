#!/usr/bin/env bash
# ============================================================================
# make_all_videos.sh
#
# Build every figure-script video for the pyBacktrack 1.5 paper FROM
# SCRATCH: wipe any cached frame PNGs, re-render every frame, then
# stitch the MP4s.  Companion to run_all_figures.sh.
#
# Targets:
#   02     - Fig 2 (LAEA on mid-N-Atlantic) per-Myr video
#   03a    - Fig 3 LAEA south-pole paleobathymetry video (default source)
#   03b    - Fig 3 LAEA south-pole paleobathymetry video, Z22 hardwired
#   10a    - NW Shelf well-derived subsidence-rate videos (A, C, D)
#   12a    - NW Shelf raw dynamic-topography video
#
# Usage:
#   ./make_all_videos.sh                  # rebuild ALL five from scratch
#   ./make_all_videos.sh 02 03a           # rebuild a subset from scratch
#   ./make_all_videos.sh --no-force 10a   # honour cached frames for 10a
#   ./make_all_videos.sh --skip-stitch    # render frames but skip ffmpeg
#
# By default every script is invoked with --force, so existing frame
# caches under figures/output/{southpole_paleobathy,northpole_paleobathy,
# nwshelf_subsidence,dt_field}/ are wiped first.
# Pass --no-force to fall back to the default "skip frames that already
# exist" policy baked into each script.
#
# Prerequisites (only flagged if the relevant target is requested):
#   - ffmpeg on PATH
#   - For 10a: figures/output/nwshelf_subsidence/{rate_*.nc,
#     all_wells_locations.txt} produced by 07a_backstrip_all_nwshelf.py
#   - For 11a: same well-locations file, plus
#     figures/output/nwshelf_subsidence/rate_delta_picked_times.txt
#     written by 09_rate_maps.py (the four times used as static-figure
#     panels; the video reuses them only as a cross-check anchor).
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PYTHON=${PYTHON:-python3}
FORCE_FLAG="--force"
SKIP_STITCH=0
TARGETS=()

for arg in "$@"; do
  case "$arg" in
    --no-force)    FORCE_FLAG="" ;;
    --skip-stitch) SKIP_STITCH=1 ;;
    02|03a|03b|10a|12a) TARGETS+=("$arg") ;;
    -h|--help)
      sed -n '1,/^set -euo pipefail/p' "$0" | sed 's/^# \{0,1\}//'
      exit 0
      ;;
    *)
      echo "Unknown argument: $arg (expected: 02, 03a, 03b, 10a, 12a, --no-force, --skip-stitch, --help)" >&2
      exit 2
      ;;
  esac
done
if [ ${#TARGETS[@]} -eq 0 ]; then
  # Default rebuild list: Fig 2 video, both Fig 3 video variants, and the
  # NW Shelf rate + DT videos.
  TARGETS=(02 03a 03b 10a 12a)
fi

# ---------------------------------------------------------------------------
# Preflight
# ---------------------------------------------------------------------------
if [ "$SKIP_STITCH" -eq 0 ] && ! command -v ffmpeg >/dev/null 2>&1; then
  echo "ERROR: ffmpeg not found on PATH.  Install it (brew install ffmpeg" >&2
  echo "       or apt install ffmpeg) or rerun with --skip-stitch." >&2
  exit 3
fi

OUT_DIR="$SCRIPT_DIR/../output"
NWSHELF_DIR="$OUT_DIR/nwshelf_subsidence"

# Per-target prerequisite checks (only fire if the target is requested).
for tgt in "${TARGETS[@]}"; do
  case "$tgt" in
    10a)
      if ! ls "$NWSHELF_DIR/A_no_sl_no_dt_rate"/rate_*.nc >/dev/null 2>&1 \
         || ! ls "$NWSHELF_DIR/C_sl_and_dt_rate"/rate_*.nc >/dev/null 2>&1 \
         || ! ls "$NWSHELF_DIR/D_rate_difference"/rate_*.nc >/dev/null 2>&1 \
         || [ ! -f "$NWSHELF_DIR/all_wells_locations.txt" ]; then
        echo "ERROR: 10a requires NW Shelf rate grids + well-locations file." >&2
        echo "       Run ./run_all_figures.sh fig08 first." >&2
        exit 4
      fi
      ;;
    12a)
      if [ ! -f "$NWSHELF_DIR/all_wells_locations.txt" ]; then
        echo "ERROR: 12a requires $NWSHELF_DIR/all_wells_locations.txt" >&2
        echo "       Run ./run_all_figures.sh fig08 first." >&2
        exit 4
      fi
      # picked-times sidecar is optional (11a falls back to a hardcoded
      # default if missing) but warn so the run-log is self-explanatory.
      if [ ! -f "$NWSHELF_DIR/rate_delta_picked_times.txt" ]; then
        echo "WARNING: $NWSHELF_DIR/rate_delta_picked_times.txt missing;" >&2
        echo "         12a will use its built-in fallback picked times." >&2
      fi
      ;;
  esac
done

# ---------------------------------------------------------------------------
# Per-target invocation
# ---------------------------------------------------------------------------
STITCH_FLAG=""
if [ "$SKIP_STITCH" -eq 1 ]; then
  STITCH_FLAG="--skip-video"
fi

NOISE_REGEX='MallocStackLogging|MallocStackLoggingNoCompact'

run_step() {
  local name="$1"; shift
  local script="$1"; shift
  echo
  echo "========================================================================"
  echo " Video build: ${name}"
  echo "   script : ${script}"
  echo "   flags  : ${FORCE_FLAG} ${STITCH_FLAG} $*"
  echo "========================================================================"
  # Filter macOS-only Malloc* stderr noise but pass real errors through.
  time $PYTHON "$script" $FORCE_FLAG $STITCH_FLAG "$@" \
      2> >(grep -E -v "$NOISE_REGEX" >&2)
}

for tgt in "${TARGETS[@]}"; do
  case "$tgt" in
    02)  run_step "02  Fig 2: LAEA on mid-N-Atlantic" \
                  "02_north_atlantic_paleobathymetry_video.py" ;;
    03a) run_step "03a Fig 3: LAEA south-pole (default source)" \
                  "03a_southern_ocean_paleobathymetry_video.py" ;;
    03b) run_step "03b Fig 3: LAEA south-pole (Z22 hardwired)" \
                  "03b_southern_ocean_paleobathymetry_video_Z22.py" ;;
    10a) run_step "10a NW Shelf subsidence-rate videos (A, C, D)" \
                  "10a_make_subsidence_videos.py" ;;
    12a) run_step "12a NW Shelf raw dynamic-topography video" \
                  "12a_dt_video.py" ;;
  esac
done

echo
echo "All requested video targets finished."
echo "  Outputs under: $OUT_DIR/"
ls -1 "$OUT_DIR"/*.mp4 2>/dev/null || echo "  (no MP4s found yet -- check the per-step logs above)"
