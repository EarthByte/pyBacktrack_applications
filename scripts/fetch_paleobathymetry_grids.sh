#!/usr/bin/env bash
# ============================================================================
# fetch_paleobathymetry_grids.sh
#
# Downloads the pyBacktrack 1.5 released paleobathymetry grid set from the
# GPlates WebDAV repository into ../../Zahirovic2022_paleobathymetry_grids/
# (so that figures/scripts/config.py can find them).
#
# These are the paleomagnetic-frame Zahirovic et al. (2022) +M7+GDH1+merged
# paleobathymetry grids linked from the pyBacktrack 1.5 documentation:
#   https://pybacktrack.readthedocs.io/en/stable/pybacktrack_paleo_bathymetry.html
#
# Usage:
#     ./fetch_paleobathymetry_grids.sh           # download paleomagnetic frame
#     ./fetch_paleobathymetry_grids.sh mantle    # download mantle frame instead
# ============================================================================
set -euo pipefail

FRAME="${1:-paleomagnetic}"
case "$FRAME" in
  paleomagnetic|pmag)
    REMOTE="https://repo.gplates.org/webdav/PlateModel_Age_SR_Grids/Zahirovic_etal_2022_GDJ/02_AgegridsUsingTopologies/PaleomagneticFrame/PaleobathymetryGrids/GDH1/"
    ;;
  mantle|mantleframe)
    REMOTE="https://repo.gplates.org/webdav/PlateModel_Age_SR_Grids/Zahirovic_etal_2022_GDJ/02_AgegridsUsingTopologies/OptimisedMantleFrame/PaleobathymetryGrids/GDH1/"
    ;;
  *)
    echo "Unknown frame: $FRAME (expected 'paleomagnetic' or 'mantle')" >&2
    exit 2 ;;
esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# Match the fallback path tried by scripts/config.py so the grids land
# in the right place for the figure scripts to find them automatically.
case "$FRAME" in
  paleomagnetic|pmag) SUBDIR="Zahirovic2022_PMag"   ;;
  mantle|mantleframe) SUBDIR="Zahirovic2022_mantle" ;;
esac
DEST="$REPO_ROOT/data/grids/paleobathymetry/$SUBDIR"
mkdir -p "$DEST"

echo "==> Downloading $FRAME-frame paleobathymetry grids"
echo "    from: $REMOTE"
echo "    into: $DEST"

# Recursive wget — keeps only the NetCDF files, strips path components, and
# does not re-download files we already have.
wget --recursive --no-parent --no-host-directories --cut-dirs=10 \
     --accept "*.nc" --reject "index.html*" \
     --no-clobber --quiet --show-progress \
     --directory-prefix="$DEST" \
     "$REMOTE"

echo
echo "Done. Grids in $DEST"
ls -1 "$DEST" | head -10
echo "..."
echo "$(ls -1 "$DEST" | wc -l) files total"
