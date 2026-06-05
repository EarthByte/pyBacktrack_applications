"""
Shared configuration for the pyBacktrack_applications figure scripts.

Layout assumed by this file:

    pyBacktrack_applications/
        data/
            wells/                              <- 109 NW Shelf wells
            grids/
                vgg_nwshelf.nc                  <- regional Sandwell V31.1 VGG
                dynamic_topography/D10_gmcm9/   <- regional D10_gmcm9 grids
        scripts/
            config.py                           <- this file
            01_workflow_flowchart.py ... 12_dt_maps.py
        figures/
            output/                             <- written by the scripts

All paths below are resolved relative to ``REPO_ROOT`` (the directory
containing both ``scripts/`` and ``data/``), so the repo can be cloned
anywhere on disk and the scripts will still find their inputs.

The two paleobathymetry grid sets (Zahirovic 2022 mantle frame and
Zahirovic 2022 paleomagnetic frame) are too large to ship in-tree
(~hundreds of MB).  They live outside the repo and are populated by
``scripts/fetch_paleobathymetry_grids.sh``; see ``PALEO_BATHY_DIR``
below.
"""
import os

# ---------------------------------------------------------------------------
# Silence the "MallocStackLogging: can't turn off malloc stack logging
# because it was not enabled" noise that pyGMT's subprocess GMT calls
# emit on macOS hundreds of times during the data step + rate-map +
# video figures.  Two layers of protection:
#
#   1. Pop any inherited Malloc*/MSL* env vars so child processes
#      start with a clean environment.
#   2. Install a file-descriptor-level stderr filter on macOS so that
#      even if a child process emits the warning anyway (which it
#      does, because something in conda's dyld load path -- probably
#      the Qt5/Qt6 mixed install -- triggers it independent of env
#      vars), the line is dropped before reaching the terminal.
#
# Opt out by setting PYBACKTRACK_RAW_STDERR=1 in your environment.
# ---------------------------------------------------------------------------
for _msl_var in ("MallocStackLogging", "MallocStackLoggingNoCompact",
                 "MallocScribble", "MallocGuardEdges",
                 "MallocCheckHeapAbort", "MallocErrorAbort"):
    os.environ.pop(_msl_var, None)


def _install_stderr_filter():
    """Replace fd 2 with a pipe whose reader strips lines that match
    NOISE_REGEX and forwards the rest to the original stderr.

    Works for ALL writes to fd 2, including writes from C subprocesses
    (GMT, ffmpeg, ...) and shared libraries -- this is what pure-Python
    sys.stderr redirection cannot do.
    """
    import re
    import sys
    import threading

    # Patterns that are pure noise.  Add more if other categories of
    # macOS stderr garbage start appearing in your runs.
    NOISE_REGEX = re.compile(
        r"MallocStackLogging|"
        r"MallocStackLoggingNoCompact"
    )

    # Duplicate fd 2 so we can keep writing the real (kept) lines to
    # the original stderr, even after we redirect fd 2 itself.
    saved_fd = os.dup(2)
    r_fd, w_fd = os.pipe()
    os.dup2(w_fd, 2)
    os.close(w_fd)
    # Force Python's sys.stderr to write into the same pipe so its
    # output is filtered too (it would write into fd 2 directly anyway
    # because of the dup2, but updating sys.stderr keeps Python's
    # internal buffering aligned).
    sys.stderr = os.fdopen(2, "w", buffering=1)

    def _pump():
        sink = os.fdopen(saved_fd, "w", buffering=1)
        src = os.fdopen(r_fd, "r")
        for line in src:
            if not NOISE_REGEX.search(line):
                try:
                    sink.write(line)
                    sink.flush()
                except Exception:
                    # Original stderr was closed (process tearing
                    # down).  Just exit the pump thread silently.
                    return

    t = threading.Thread(target=_pump, daemon=True,
                         name="pybacktrack-stderr-filter")
    t.start()


# Install on macOS unless the user explicitly opts out.  No-op on Linux
# (where this noise doesn't exist) and Windows (where pyGMT scripts are
# rarely run for this project).
import sys as _sys
if (_sys.platform == "darwin"
        and not os.environ.get("PYBACKTRACK_RAW_STDERR")):
    try:
        _install_stderr_filter()
    except Exception as _exc:
        # Don't take the whole figure run down if the filter fails to
        # install -- just warn and keep going.
        print(f"[pybacktrack-figs] WARNING: could not install stderr "
              f"filter ({_exc}); MallocStackLogging noise will be "
              "visible.", flush=True)

# Root of the repo (the directory containing scripts/ and data/).
# Kept as PAPER_ROOT for backwards-compat with downstream scripts that
# still import that name.
REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir)
)
PAPER_ROOT = REPO_ROOT     # alias

# Repo-local data tree -- populated by tools/populate_data.py on first run.
DATA_DIR = os.path.join(REPO_ROOT, "data")
WELLS_DIR_INREPO = os.path.join(DATA_DIR, "wells")
GRIDS_DIR_INREPO = os.path.join(DATA_DIR, "grids")

# ---------------------------------------------------------------------------
# Paleobathymetry grid source.
#
# Toggle between the two Zahirovic et al. (2022) GDH1 grid sets we
# support:
#
#   "Z22_PMag"    - the Z22 grids in the *paleomagnetic* reference frame.
#                   Default for this repo (uses John Cannon's 2026
#                   pybacktrack-nearneighbor re-grid when available on disk).
#   "Z22_mantle"  - the Z22 grids in the optimised *mantle* reference frame.
#
# Both share the same topology model (Zahirovic 2022 in the GPlates
# registry); only the anchor plate differs (0 for the mantle frame,
# 701701 for the paleomagnetic frame).
#
# Each branch sets:
#   PALEO_BATHY_DIR             - directory containing the .nc grids
#   PLATE_MODEL_ANCHOR_PLATE    - anchor plate matching the grid frame
# ---------------------------------------------------------------------------
# The env var below lets a script override the source for a single run
# without editing this file (e.g. the global-paleobathy video script
# accepts --source on the command line and sets this env var before
# importing config).  Leave PALEO_BATHY_SOURCE_DEFAULT as the day-to-day
# default; the env var takes precedence.
PALEO_BATHY_SOURCE_DEFAULT = "Z22_PMag"           # "Z22_PMag" | "Z22_mantle"
PALEO_BATHY_SOURCE = os.environ.get(
    "PYBACKTRACK_PALEO_BATHY_SOURCE", PALEO_BATHY_SOURCE_DEFAULT)

if PALEO_BATHY_SOURCE == "Z22_mantle":
    # Zahirovic et al. (2022) GDH1 grids, optimised mantle reference frame.
    # Topology model: `Zahirovic2022` in the GPlates registry; anchor 0
    # because this is a mantle-frame product.
    #
    # Resolution order:
    #   1. $PYBT_PALEO_BATHY_DIR override (if set)
    #   2. data/grids/paleobathymetry/Zahirovic2022_mantle/ in the repo
    #      (populated by `scripts/fetch_paleobathymetry_grids.sh mantle`)
    #   3. legacy absolute path on Dietmar's laptop (back-compat)
    _candidates = [
        os.environ.get("PYBT_PALEO_BATHY_DIR"),
        os.path.join(GRIDS_DIR_INREPO, "paleobathymetry",
                     "Zahirovic2022_mantle"),
        ("/Users/dietmar/Documents/GPlates/pyBacktrack1.5/"
         "Paleobathymetry_Z22_mantle"),
        os.path.join(PAPER_ROOT, "Zahirovic2022_paleobathymetry_grids_mantle"),
    ]
    PALEO_BATHY_DIR = next(
        (p for p in _candidates if p and os.path.isdir(p)),
        _candidates[1])   # fall through to repo path even if missing
    _PALEO_BATHY_PLATE_MODEL = "Zahirovic2022"
    _PALEO_BATHY_ANCHOR = 0          # mantle reference frame

elif PALEO_BATHY_SOURCE == "Z22_PMag":
    # Zahirovic et al. (2022) GDH1 grids, paleomagnetic reference frame.
    # Topology model: `Zahirovic2022` in the GPlates registry; anchor
    # 701701 pins the paleomagnetic frame.
    #
    # Resolution order (first that exists on disk wins):
    #   1. $PYBT_PALEO_BATHY_DIR override (if set)
    #   2. data/grids/paleobathymetry/Zahirovic2022_PMag/ (in-repo;
    #      populated by `scripts/fetch_paleobathymetry_grids.sh paleomagnetic`)
    #   3. John Cannon's 2026 pybacktrack-nearneighbor re-grid on Dietmar's
    #      laptop -- this is what the paper figures use by default
    #   4. legacy mirror of the public Z22 release
    #
    # Public-release source for option (2):
    #   https://repo.gplates.org/webdav/PlateModel_Age_SR_Grids/
    #   Zahirovic_etal_2022_GDJ/02_AgegridsUsingTopologies/
    #   PaleomagneticFrame/PaleobathymetryGrids/GDH1/
    # PAPER_ROOT = pyBacktrack_applications/.  John's nearneighbour
    # re-grid lives in the parent pyBacktrack1.5/ tree, one level up.
    _PYB15_ROOT = os.path.dirname(PAPER_ROOT)
    _candidates = [
        os.environ.get("PYBT_PALEO_BATHY_DIR"),
        os.path.join(GRIDS_DIR_INREPO, "paleobathymetry",
                     "Zahirovic2022_PMag"),
        # John's nearneighbour re-grid in the parent paper tree:
        os.path.join(_PYB15_ROOT,
                     "Zahirovic2022_PMag_paleobathymetry_pybacktrack-nearneighbor"),
        # Legacy laptop mirror of the public Z22 release:
        ("/Users/dietmar/Documents/GPlates/pyBacktrack1.5/"
         "Paleobathymetry_Z22_PMag"),
        os.path.join(_PYB15_ROOT, "Zahirovic2022_paleobathymetry_grids"),
    ]
    PALEO_BATHY_DIR = next(
        (p for p in _candidates if p and os.path.isdir(p)),
        _candidates[1])
    _PALEO_BATHY_PLATE_MODEL = "Zahirovic2022"
    _PALEO_BATHY_ANCHOR = 701701     # paleomagnetic frame

else:
    raise ValueError(
        f"Unknown PALEO_BATHY_SOURCE: {PALEO_BATHY_SOURCE!r}. "
        "Set to 'Z22_mantle' or 'Z22_PMag'."
    )

# Filename template for the released grids.  The default matches the
# convention used by `pybacktrack.reconstruct_paleo_bathymetry_grids` (one
# decimal place for the time).  If the grids in PALEO_BATHY_DIR use a
# different naming convention (e.g. integer Ma, or `paleobathymetry_<time>Ma.nc`),
# we attempt to auto-detect the pattern from any .nc file present.  Edit the
# fallback string by hand if your filenames don't fit one of the candidate
# patterns tried below.
def _detect_paleo_bathy_fmt(directory):
    import glob, re
    candidates = [
        "paleo_bathymetry_{time:.1f}.nc",   # default pyBacktrack output
        "paleo_bathymetry_{time:.0f}.nc",   # integer-Ma variant
        "paleobathymetry_{time:.0f}Ma.nc",  # GPlates repo style A
        "paleobathymetry_{time:.1f}Ma.nc",  # GPlates repo style B
        "paleobathymetry_{time:.0f}.nc",    # bare integer Ma
    ]
    # Probe each candidate with a few common times.
    if os.path.isdir(directory):
        for fmt in candidates:
            for probe in (0, 60, 100):
                if os.path.exists(os.path.join(
                        directory, fmt.format(time=float(probe)))):
                    return os.path.join(directory, fmt)
        # Last resort: inspect any .nc file and derive the format from it.
        ncs = sorted(glob.glob(os.path.join(directory, "*.nc")))
        if ncs:
            base = os.path.basename(ncs[0])
            m = re.search(r"(\d+)(\.\d+)?", base)
            if m:
                if m.group(2):  # has decimal
                    fmt = re.sub(r"\d+\.\d+", "{time:.1f}", base, count=1)
                else:
                    fmt = re.sub(r"\d+", "{time:.0f}", base, count=1)
                return os.path.join(directory, fmt)
    # Fall through to the pyBacktrack default; downstream scripts will
    # error informatively if no matching file is found.
    return os.path.join(directory, candidates[0])


PALEO_BATHY_FMT = _detect_paleo_bathy_fmt(PALEO_BATHY_DIR)

# Plate model + anchor.  Both come from the PALEO_BATHY_SOURCE branch
# above so the overlay model + frame always match the grid model + frame
# (otherwise coastlines drift away from bathymetric features).
#   - "Z22_mantle" source -> Zahirovic2022 model + anchor 0       (mantle)
#   - "Z22_PMag"   source -> Zahirovic2022 model + anchor 701701 (pmag)
PLATE_MODEL_NAME = _PALEO_BATHY_PLATE_MODEL
PLATE_MODEL_ANCHOR_PLATE = _PALEO_BATHY_ANCHOR

# Sanity print on import so every figure run logs which frame is active
# (silent mismatches are the worst kind -- coastlines look "roughly
# right" but are tens to hundreds of km off the bathymetric features).
_ANCHOR_FRAME = {0: "mantle", 701701: "paleomagnetic"}.get(
    PLATE_MODEL_ANCHOR_PLATE, f"anchor={PLATE_MODEL_ANCHOR_PLATE}")
print(
    f"[pybacktrack-figs] PALEO_BATHY_SOURCE = {PALEO_BATHY_SOURCE} "
    f"-> plate model {PLATE_MODEL_NAME}, anchor {PLATE_MODEL_ANCHOR_PLATE} "
    f"({_ANCHOR_FRAME} frame)",
)

# NW Shelf well files (lithology format consumed by pybacktrack.backstrip_well).
# Lives in-repo under data/wells/ after running tools/populate_data.py.
NWSHELF_WELL_DIR = WELLS_DIR_INREPO

# Featured well for the comprehensive single-well example.
#
# We chose "asteras" because all 20 of its stratigraphic entries fall
# inside 0-150 Ma (oldest = 146 Ma, youngest = 24 Ma), which matches the
# coverage of the D10_gmcm9 dynamic topography model exactly.  Dense
# Cretaceous + Cenozoic record, mid-NW-Shelf location.
#
# Previous featured wells (Dillon Shoals, Plover) were dropped on
# 2026-05-27 when the DT model was switched from gld428 to D10_gmcm9 --
# both had Permian-Triassic syn-rift records that extended well beyond
# the D10_gmcm9 0-150 Ma range, making meaningful comparison impossible.
WELLS = [
    {
        "name": "Asteras",
        "file": os.path.join(NWSHELF_WELL_DIR, "asteras.txt"),
        "max_age": 146.0,
        "rift_start_age": 155.0,       # purely Cretaceous-onwards record;
        "rift_end_age": 130.0,         # post-syn-rift in this well
    },
]

# Output directory for all figures (PNG, PDF).
OUTPUT_DIR = os.path.join(REPO_ROOT, "figures", "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Regional vertical-gravity-gradient grid (Sandwell & Smith V31.1, cropped
# to the Fig 5 NW Shelf region by tools/populate_data.py).
VGG_GRID = os.path.join(GRIDS_DIR_INREPO, "vgg_nwshelf.nc")

# Regional dynamic-topography grids (Braz et al. 2021, D10_gmcm9, cropped
# to the Fig 12 NW Shelf region + buffer by tools/populate_data.py).
DT_GRIDS_DIR = os.path.join(
    GRIDS_DIR_INREPO, "dynamic_topography", "D10_gmcm9")

# Straume et al. (2020) paleobathy+topo grids, used by
# `04_pybacktrack_vs_straume_comparison.py` as an independent reference
# to benchmark the released pyBacktrack 1.5 paleobathymetry against.
# 0.1 deg, 1 Myr, 1-65 Ma, full-earth (no NaN, includes topography).
STRAUME_GRID_DIR = os.path.join(PAPER_ROOT, "Straume_paleogeography_grids")
STRAUME_GRID_FMT = os.path.join(
    STRAUME_GRID_DIR, "paleobathy-topo_{time:.2f}Ma_Straume_et_al.nc")

# Dynamic topography model used for the NW Shelf well figures.
#
# We use 'D10_gmcm9' (Braz et al. 2021), bundled with pyBacktrack at
# pybacktrack/bundle_data/dynamic_topography/models/Braz2021/D10_gmcm9/.
# It covers 0-150 Ma in 5 Myr time slices.  Switched from 'gld428' on
# 2026-05-27 because gld428 had too many artefacts at the NW Shelf scale.
# D10_gmcm9's 150 Ma cutoff drove the choice of asteras as the featured
# well (record only goes back to 146 Ma -- fully covered).
DYNAMIC_TOPOGRAPHY_MODEL = "D10_gmcm9"

# Sea level model used for the NW Shelf well figures.
#
# We use the bundled long-term hybrid Phanerozoic curve of Haq and Ogg (2024),
# which covers the full record of the NW Shelf wells with a smoothed,
# highstand-following envelope appropriate for long-period subsidence analysis.
SEA_LEVEL_MODEL = "Haq2024_Hybrid_SealevelCurve_Longterm"

# Hard upper bound for every x-axis on the NW Shelf well + margin-scale
# figures.  Matches the D10_gmcm9 dynamic-topography coverage and is
# also the asteras well's effective record length.  Time runs LEFT to
# RIGHT in every figure, so the x-axis is (MAX_ANALYSIS_AGE, 0) -- 150
# Ma on the left, 0 Ma on the right.
MAX_ANALYSIS_AGE = 150.0
