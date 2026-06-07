#!/usr/bin/env python3
"""
Data step for the margin-scale subsidence-rate figures (Figs 8, 9,
10).  Backstrips every well in ``data/wells/`` (see
``config.NWSHELF_WELL_DIR``) in two configurations and grids both the
tectonic subsidence and the forward 1-Myr subsidence rate at every time
slice, then writes their pixel-wise difference grid.

Configurations:

    A.  no sea level, no dynamic topography
    C.  bundled Haq2024 long-term hybrid sea-level curve + D10_gmcm9
        dynamic topography (Braz et al. 2021; 0-150 Ma in 5 Myr
        steps).  Wells reaching beyond 150 Ma have DT samples set
        to 0 outside the model coverage, so their syn-rift records
        remain in the population but the analysis is meaningfully
        constrained to 0-150 Ma.

For each 1 Myr time slice the per-well subsidence point cloud is
gridded with GMT `surface` and clipped with a buffered hull mask.
The difference grid

    D = C - A                      (full corrections minus no SL/DT)

quantifies how much of the inferred tectonic subsidence (or rate)
at each pixel is absorbed by the SL + DT corrections.

Defensible time cutoff
----------------------
The well population thins back through time (~109 wells at 50 Ma,
~107 at 100 Ma, ~68 at 200 Ma, only a handful at 250 Ma).  The
script automatically picks the oldest defensible time slice -- the
oldest age at which the number of contributing wells is >=
`MIN_WELLS_FOR_GRID` AND sites span at least `MIN_LON_RANGE` deg
of longitude and `MIN_LAT_RANGE` deg of latitude -- and only
writes grids out to that age.  The cutoff is printed at the end
of the run and recorded in `cutoff_time_Ma.txt`.

Outputs (under figures/output/nwshelf_subsidence/)
--------------------------------------------------
A_no_sl_no_dt/   data_<wellname>.txt
A_no_sl_no_dt/   subsidence_<time>.nc
C_sl_and_dt/     data_<wellname>.txt
C_sl_and_dt/     subsidence_<time>.nc
D_difference/    subsidence_<time>.nc          (= C - A)
all_wells_locations.txt
all_wells_subsidence_<config>.csv
well_counts_per_time.csv
cutoff_time_Ma.txt
"""
import glob
import os
import sys
import time
import warnings

import numpy as np
import pygmt
import xarray as xr
import pybacktrack

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    NWSHELF_WELL_DIR, OUTPUT_DIR,
    DYNAMIC_TOPOGRAPHY_MODEL, SEA_LEVEL_MODEL,
)

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
MAX_TIME = 150                # Myr - hard upper bound matching D10_gmcm9
                              # coverage.  Actual cutoff is data-driven
                              # (see MIN_WELLS_FOR_GRID below) so this is
                              # just a ceiling -- it stops the script from
                              # backstripping the Permo-Triassic portion
                              # of long-record wells that the DT model
                              # no longer covers.
TIME_STEP = 1.0
TIMES_FULL = np.arange(0, MAX_TIME + 1, TIME_STEP, dtype=float)

# Map region (degrees).  Includes a buffer beyond the data hull for
# surface gridding.
REGION = [113, 132, -22, -8]
GRID_SPACING = 0.1
SURFACE_TENSION = 0.25
MASK_BUFFER_KM = 80

# Defensible-time-cutoff thresholds.
MIN_WELLS_FOR_GRID = 10       # at least this many wells contribute
                              # (relaxed from 15 to 10 on 2026-06-02
                              # to extend the temporal coverage of
                              # the rate grids)
MIN_LON_RANGE_DEG = 3.0       # spread over at least 3 deg of longitude
MIN_LAT_RANGE_DEG = 3.0       # spread over at least 3 deg of latitude

# Sub-output directory tree.
OUT_BASE = os.path.join(OUTPUT_DIR, "nwshelf_subsidence")
CONFIGS = [
    ("A_no_sl_no_dt", "A. no sea level, no dynamic topography",
     dict(sea_level_model=None, dynamic_topography_model=None)),
    ("C_sl_and_dt",
     f"C. Haq2024 long-term SL + {DYNAMIC_TOPOGRAPHY_MODEL} dynamic topography",
     dict(sea_level_model=SEA_LEVEL_MODEL,
          dynamic_topography_model=DYNAMIC_TOPOGRAPHY_MODEL)),
]
DIFF_DIR = "D_difference"
RATE_SUFFIX = "_rate"          # appended to config dirs for rate grids
RATE_DIFF_DIR = "D_rate_difference"
for cfg_name, _, _ in CONFIGS:
    os.makedirs(os.path.join(OUT_BASE, cfg_name), exist_ok=True)
    os.makedirs(os.path.join(OUT_BASE, cfg_name + RATE_SUFFIX), exist_ok=True)
os.makedirs(os.path.join(OUT_BASE, DIFF_DIR), exist_ok=True)
os.makedirs(os.path.join(OUT_BASE, RATE_DIFF_DIR), exist_ok=True)

warnings.filterwarnings(
    "ignore", "Well thickness .* is larger than the total sediment thickness")
warnings.filterwarnings("ignore", "Dynamic topography model .* cannot")


# --------------------------------------------------------------------------
# Buffer mask helper.
#
# Implementation note: pyGMT 0.18.0 does NOT expose `pygmt.grdmask` as a
# top-level function (there is no pygmt/src/grdmask.py).  Rather than
# wrap a low-level Session.call_module() invocation with on-disk temp
# files, we build the mask in pure numpy.  At our scale -- ~190 x 140
# grid cells and 109 well points -- the broadcast distance test is
# essentially instant (~3 M float operations) and has no Python-package
# dependencies beyond numpy.
# --------------------------------------------------------------------------
def build_buffer_mask(locs, region, spacing, buffer_km):
    """Return a 2-D float numpy array over (region, spacing) with 1.0
    inside `buffer_km` of any (lon, lat) point in `locs` and NaN
    elsewhere.

    Uses an equirectangular-km approximation centred on the region's
    mid-latitude.  Errors at the NW Shelf corners are <1 percent of
    the buffer radius -- negligible for our 80 km buffer.
    """
    lonmin, lonmax, latmin, latmax = region
    lons = np.arange(lonmin, lonmax + 1e-6, spacing)
    lats = np.arange(latmin, latmax + 1e-6, spacing)
    LON, LAT = np.meshgrid(lons, lats)

    lat0_rad = np.deg2rad(0.5 * (latmin + latmax))
    km_per_deg_lat = 111.32
    km_per_deg_lon = 111.32 * np.cos(lat0_rad)

    well_lon = locs[:, 0]
    well_lat = locs[:, 1]
    # Broadcast: (n_lat, n_lon, n_wells)
    dlon_km = (LON[..., None] - well_lon[None, None, :]) * km_per_deg_lon
    dlat_km = (LAT[..., None] - well_lat[None, None, :]) * km_per_deg_lat
    min_dist_sq = (dlon_km ** 2 + dlat_km ** 2).min(axis=-1)
    return np.where(min_dist_sq <= buffer_km ** 2, 1.0, np.nan)


def _apply_buffer_mask(grid, mask_arr):
    """Multiply a pygmt.surface output DataArray by a 2-D numpy mask,
    handling either lon/lat or lat/lon axis order between pyGMT releases.

    Returns a new DataArray with the same coords/dims as `grid`.
    """
    arr = grid.values.astype(float)
    if arr.shape == mask_arr.shape:
        m = mask_arr
    elif arr.shape == mask_arr.T.shape:
        m = mask_arr.T
    else:
        raise ValueError(
            f"mask shape {mask_arr.shape} does not match grid shape "
            f"{arr.shape} (or its transpose)"
        )
    masked_vals = arr * m
    return xr.DataArray(
        masked_vals, coords=grid.coords, dims=grid.dims, name=grid.name,
    )


# --------------------------------------------------------------------------
# 1. Backstrip every well in every configuration
# --------------------------------------------------------------------------
def read_well_metadata(path):
    lat = lon = None
    with open(path) as fh:
        for _ in range(10):
            line = fh.readline()
            if not line:
                break
            if 'SiteLatitude' in line:
                try: lat = float(line.split('=')[-1])
                except ValueError: pass
            if 'SiteLongitude' in line:
                try: lon = float(line.split('=')[-1])
                except ValueError: pass
    return lat, lon


def backstrip_single(well_file, cfg_kwargs):
    """Run backstripping at the 1 Myr cadence requested by TIMES_FULL.

    Mirrors the API used in the bundled ``geohistory_analysis.ipynb``
    notebook (cell 13): pass `times=` to backstrip_well to get decompacted
    sections at every requested time.  Dynamic topography is sampled
    separately and subtracted from the tectonic-subsidence trace, exactly
    as the bundled notebook does.
    """
    kwargs = dict(lithology_filenames=pybacktrack.BUNDLE_LITHOLOGY_FILENAMES)
    if cfg_kwargs.get("sea_level_model"):
        kwargs["sea_level_model"] = cfg_kwargs["sea_level_model"]

    well, decomp = pybacktrack.backstrip_well(
        well_file, times=TIMES_FULL.tolist(), **kwargs)

    dt_sampler = None
    if cfg_kwargs.get("dynamic_topography_model"):
        dt_sampler = pybacktrack.DynamicTopography.create_from_bundled_model(
            cfg_kwargs["dynamic_topography_model"],
            well.longitude, well.latitude,
        )

    rows = []
    for ds in decomp:
        age = ds.get_age()
        ts = ds.get_tectonic_subsidence()
        if dt_sampler is not None:
            dt_val = dt_sampler.sample(float(age))
            if dt_val is not None and not np.isnan(dt_val):
                ts = ts - float(dt_val)
        rows.append((age, well.longitude, well.latitude, ts))
    return rows


def run_all_backstripping():
    well_files = (sorted(glob.glob(os.path.join(NWSHELF_WELL_DIR, "*.txt"))) +
                  sorted(glob.glob(os.path.join(NWSHELF_WELL_DIR, "*.dat"))))
    print(f"Found {len(well_files)} well files in {NWSHELF_WELL_DIR}")

    # Locations file
    loc_path = os.path.join(OUT_BASE, "all_wells_locations.txt")
    with open(loc_path, "w") as out:
        out.write("# name lon lat\n")
        for wf in well_files:
            lat, lon = read_well_metadata(wf)
            if lat is None or lon is None:
                continue
            name = os.path.splitext(os.path.basename(wf))[0]
            out.write(f"{name} {lon:.5f} {lat:.5f}\n")
    print(f"  wrote {loc_path}")

    for cfg_name, cfg_title, cfg_kwargs in CONFIGS:
        print(f"\n=== Configuration {cfg_title} ===")
        csv_path = os.path.join(OUT_BASE, f"all_wells_subsidence_{cfg_name}.csv")
        # Resume-from-CSV: backstripping all 109 wells in two configs takes
        # ~50 minutes per config; if the CSV is already on disk from a
        # previous (interrupted) run we skip the per-well loop entirely.
        # Force a fresh backstripping run by deleting the CSV or by setting
        # PYBACKTRACK_FORCE_BACKSTRIP=1 in your shell.
        if (os.path.exists(csv_path)
                and os.path.getsize(csv_path) > 0
                and not os.environ.get("PYBACKTRACK_FORCE_BACKSTRIP")):
            print(f"  {csv_path} already exists "
                  f"({os.path.getsize(csv_path) / 1024:.0f} KB) -- "
                  "skipping backstripping for this config.  "
                  "Set PYBACKTRACK_FORCE_BACKSTRIP=1 to rerun.")
            continue
        # Per-well heartbeat so a stalled run is visible; flushing the
        # CSV every well so the on-disk size also reflects progress.
        t_start = time.time()
        with open(csv_path, "w") as out_csv:
            out_csv.write("# time_Ma lon lat tectonic_subsidence_m well\n")
            for i, wf in enumerate(well_files, start=1):
                name = os.path.splitext(os.path.basename(wf))[0]
                t_well = time.time()
                try:
                    rows = backstrip_single(wf, cfg_kwargs)
                except Exception as exc:
                    print(f"  [{i:>3}/{len(well_files)}] ! skipping "
                          f"{name}: {exc}", flush=True)
                    continue
                for age, lon, lat, ts in rows:
                    out_csv.write(
                        f"{age:.1f} {lon:.5f} {lat:.5f} {ts:.2f} {name}\n")
                out_csv.flush()
                elapsed = time.time() - t_well
                # Print every well -- it's 1 line/30 s on average so
                # not noisy and gives an immediate "still alive" signal.
                eta_min = (time.time() - t_start) / i * (
                    len(well_files) - i) / 60.0
                print(f"  [{i:>3}/{len(well_files)}] {name} "
                      f"({elapsed:5.1f} s/well, "
                      f"ETA {eta_min:5.1f} min for this config)",
                      flush=True)
        print(f"  wrote {csv_path}")


# --------------------------------------------------------------------------
# 2. Determine the defensible time cutoff
# --------------------------------------------------------------------------
def determine_time_cutoff():
    """Find the oldest time with a defensible well population.

    Returns the integer cutoff time in Ma.  Also writes a per-time
    well-count CSV for the methods section of the paper.
    """
    # Use config A (which always succeeds because no SL/DT is required).
    csv_a = os.path.join(OUT_BASE, "all_wells_subsidence_A_no_sl_no_dt.csv")
    data = np.loadtxt(csv_a, usecols=(0, 1, 2))     # time, lon, lat
    print("\nDetermining defensible time cutoff...")
    counts_path = os.path.join(OUT_BASE, "well_counts_per_time.csv")
    cutoff = 0
    with open(counts_path, "w") as out:
        out.write("# time_Ma n_wells lon_range_deg lat_range_deg defensible\n")
        for t in TIMES_FULL:
            sel = np.abs(data[:, 0] - t) < 0.5
            lons = data[sel, 1]; lats = data[sel, 2]
            n = len(lons)
            lon_rng = (lons.max() - lons.min()) if n else 0.0
            lat_rng = (lats.max() - lats.min()) if n else 0.0
            ok = (n >= MIN_WELLS_FOR_GRID and
                  lon_rng >= MIN_LON_RANGE_DEG and
                  lat_rng >= MIN_LAT_RANGE_DEG)
            out.write(f"{int(t)} {n} {lon_rng:.2f} {lat_rng:.2f} "
                      f"{'Y' if ok else 'N'}\n")
            if ok and t > cutoff:
                cutoff = int(t)
    print(f"  wrote {counts_path}")
    print(f"  defensible cutoff = {cutoff} Ma "
          f"(thresholds: n_wells >= {MIN_WELLS_FOR_GRID}, "
          f"lon range >= {MIN_LON_RANGE_DEG}, "
          f"lat range >= {MIN_LAT_RANGE_DEG})")
    with open(os.path.join(OUT_BASE, "cutoff_time_Ma.txt"), "w") as out:
        out.write(f"{cutoff}\n")
    return cutoff


# --------------------------------------------------------------------------
# 3. Grid the subsidence at every Myr for each configuration
# --------------------------------------------------------------------------
def grid_subsidence(cutoff_time):
    """Grid every time slice from 0 to cutoff_time (inclusive)."""
    loc_path = os.path.join(OUT_BASE, "all_wells_locations.txt")
    locs = np.loadtxt(loc_path, usecols=(1, 2), skiprows=1)

    print(f"\nBuilding hull mask buffer = {MASK_BUFFER_KM} km around "
          f"{len(locs)} wells (numpy broadcast)")
    mask = build_buffer_mask(locs, REGION, GRID_SPACING, MASK_BUFFER_KM)

    times = np.arange(0, cutoff_time + 1, int(TIME_STEP))
    for cfg_name, cfg_title, _ in CONFIGS:
        cfg_dir = os.path.join(OUT_BASE, cfg_name)
        csv_path = os.path.join(OUT_BASE,
                                f"all_wells_subsidence_{cfg_name}.csv")
        print(f"\n=== Gridding {cfg_title} (0-{cutoff_time} Ma) ===")
        data = np.loadtxt(csv_path, usecols=(0, 1, 2, 3))
        for t in times:
            sel = np.abs(data[:, 0] - t) < 0.5
            sub = data[sel]
            if len(sub) < MIN_WELLS_FOR_GRID:
                continue
            tbl = np.column_stack([sub[:, 1], sub[:, 2], sub[:, 3]])
            try:
                # Pre-process with blockmean (GMT best practice -- collapses
                # any same-cell points to their cell mean and silences
                # "You should have pre-processed ..." warnings from surface).
                tbl_bm = pygmt.blockmean(
                    data=tbl, region=REGION, spacing=GRID_SPACING)
                grid = pygmt.surface(
                    data=tbl_bm, region=REGION, spacing=GRID_SPACING,
                    tension=SURFACE_TENSION,
                )
                grid = _apply_buffer_mask(grid, mask)
                out_nc = os.path.join(cfg_dir, f"subsidence_{int(t)}.nc")
                grid.to_netcdf(out_nc)
                if int(t) % 25 == 0:
                    print(f"  {t:>5.0f} Ma -{len(sub):3d} wells -> {out_nc}")
            except Exception as exc:
                print(f"  ! gridding failed at {t} Ma: {exc}")


# --------------------------------------------------------------------------
# 4. Difference grid: D = C - A
# --------------------------------------------------------------------------
def make_difference_grids(cutoff_time):
    import xarray as xr

    print(f"\n=== Computing difference grids D = C - A (0-{cutoff_time} Ma) ===")
    a_dir = os.path.join(OUT_BASE, "A_no_sl_no_dt")
    c_dir = os.path.join(OUT_BASE, "C_sl_and_dt")
    d_dir = os.path.join(OUT_BASE, DIFF_DIR)

    times = np.arange(0, cutoff_time + 1, int(TIME_STEP))
    written = 0
    for t in times:
        a_path = os.path.join(a_dir, f"subsidence_{int(t)}.nc")
        c_path = os.path.join(c_dir, f"subsidence_{int(t)}.nc")
        if not (os.path.exists(a_path) and os.path.exists(c_path)):
            continue
        try:
            a = xr.open_dataset(a_path).z
            c = xr.open_dataset(c_path).z
            d = c - a              # sign convention C - A (user 2026-06-02)
            d_path = os.path.join(d_dir, f"subsidence_{int(t)}.nc")
            d.to_netcdf(d_path)
            written += 1
            if int(t) % 25 == 0:
                print(f"  {t:>5.0f} Ma -> {d_path}")
        except Exception as exc:
            print(f"  ! difference failed at {t} Ma: {exc}")
    print(f"  wrote {written} difference grids in {d_dir}")


# --------------------------------------------------------------------------
# 5. Per-well subsidence RATES (forward 1-Myr difference) and rate grids
# --------------------------------------------------------------------------
def write_rate_csv(cfg_name):
    """Derive a per-well subsidence-rate CSV from the per-well subsidence CSV.

    Rate convention:
        rate(t)  =  subsidence(t-1)  -  subsidence(t)    [m / Myr]
    where ``t`` is the age in Ma and t-1 is 1 Myr younger.  A positive rate
    means the site was subsiding (going deeper) during that 1 Myr interval;
    negative = uplifting.  Rate is reported at the midpoint of the 1 Myr
    interval (i.e. age = t - 0.5 Ma, rounded to the nearest integer for the
    output grid times).
    """
    in_csv = os.path.join(OUT_BASE, f"all_wells_subsidence_{cfg_name}.csv")
    out_csv = os.path.join(OUT_BASE, f"all_wells_rate_{cfg_name}.csv")
    data = np.loadtxt(in_csv, usecols=(0, 1, 2, 3))
    # Re-load names from column 4 (strings).
    names = np.loadtxt(in_csv, usecols=(4,), dtype=str)

    with open(out_csv, "w") as out:
        out.write("# time_Ma lon lat rate_m_per_Myr well\n")
        for w in np.unique(names):
            sel = (names == w)
            t = data[sel, 0]
            lon = data[sel, 1]
            lat = data[sel, 2]
            s = data[sel, 3]
            # Sort by time ascending.
            order = np.argsort(t)
            t = t[order]; s = s[order]
            lon = lon[order]; lat = lat[order]
            if len(t) < 2:
                continue
            # Forward difference: rate at midpoint between t[i] and t[i+1]
            for i in range(len(t) - 1):
                if abs(t[i + 1] - t[i] - 1.0) > 0.1:
                    continue        # only use successive 1 Myr pairs
                mid_t = 0.5 * (t[i] + t[i + 1])
                # rate convention above:  sub(younger) - sub(older)
                # younger time = smaller Ma value = t[i]
                rate = s[i] - s[i + 1]
                out.write(f"{mid_t:.1f} {lon[i]:.5f} {lat[i]:.5f} "
                          f"{rate:.3f} {w}\n")
    print(f"  wrote {out_csv}")


def grid_rates(cutoff_time):
    """Grid the per-well subsidence rates at every integer Myr."""
    loc_path = os.path.join(OUT_BASE, "all_wells_locations.txt")
    locs = np.loadtxt(loc_path, usecols=(1, 2), skiprows=1)
    mask = build_buffer_mask(locs, REGION, GRID_SPACING, MASK_BUFFER_KM)

    times = np.arange(0, cutoff_time + 1, int(TIME_STEP))
    for cfg_name, cfg_title, _ in CONFIGS:
        rate_csv = os.path.join(OUT_BASE, f"all_wells_rate_{cfg_name}.csv")
        rate_dir = os.path.join(OUT_BASE, cfg_name + RATE_SUFFIX)
        print(f"\n=== Gridding RATE for {cfg_title} (0-{cutoff_time} Ma) ===")
        data = np.loadtxt(rate_csv, usecols=(0, 1, 2, 3))
        for t in times:
            # Match midpoints within ±0.6 Myr of integer t.
            sel = np.abs(data[:, 0] - t) < 0.6
            sub = data[sel]
            if len(sub) < MIN_WELLS_FOR_GRID:
                continue
            tbl = np.column_stack([sub[:, 1], sub[:, 2], sub[:, 3]])
            try:
                # Same blockmean pre-processing as in grid_subsidence above.
                tbl_bm = pygmt.blockmean(
                    data=tbl, region=REGION, spacing=GRID_SPACING)
                grid = pygmt.surface(
                    data=tbl_bm, region=REGION, spacing=GRID_SPACING,
                    tension=SURFACE_TENSION,
                )
                grid = _apply_buffer_mask(grid, mask)
                out_nc = os.path.join(rate_dir, f"rate_{int(t)}.nc")
                grid.to_netcdf(out_nc)
                if int(t) % 25 == 0:
                    print(f"  {t:>5.0f} Ma -{len(sub):3d} wells -> {out_nc}")
            except Exception as exc:
                print(f"  ! gridding failed at {t} Ma: {exc}")


def make_rate_difference_grids(cutoff_time):
    """Difference of rate grids: rate_C - rate_A."""
    import xarray as xr

    print(f"\n=== Computing rate-difference grids (rate_C - rate_A) "
          f"(0-{cutoff_time} Ma) ===")
    a_dir = os.path.join(OUT_BASE, "A_no_sl_no_dt" + RATE_SUFFIX)
    c_dir = os.path.join(OUT_BASE, "C_sl_and_dt" + RATE_SUFFIX)
    d_dir = os.path.join(OUT_BASE, RATE_DIFF_DIR)

    times = np.arange(0, cutoff_time + 1, int(TIME_STEP))
    written = 0
    for t in times:
        a_path = os.path.join(a_dir, f"rate_{int(t)}.nc")
        c_path = os.path.join(c_dir, f"rate_{int(t)}.nc")
        if not (os.path.exists(a_path) and os.path.exists(c_path)):
            continue
        try:
            a = xr.open_dataset(a_path).z
            c = xr.open_dataset(c_path).z
            d = c - a              # sign convention C - A (user 2026-06-02)
            d_path = os.path.join(d_dir, f"rate_{int(t)}.nc")
            d.to_netcdf(d_path)
            written += 1
            if int(t) % 25 == 0:
                print(f"  {t:>5.0f} Ma -> {d_path}")
        except Exception as exc:
            print(f"  ! rate-difference failed at {t} Ma: {exc}")
    print(f"  wrote {written} rate-difference grids in {d_dir}")


# --------------------------------------------------------------------------
def main():
    run_all_backstripping()
    cutoff = determine_time_cutoff()
    grid_subsidence(cutoff)
    make_difference_grids(cutoff)

    # -- new: subsidence rates --
    print("\nDeriving per-well subsidence rates (1 Myr forward differences)")
    for cfg_name, _, _ in CONFIGS:
        write_rate_csv(cfg_name)
    grid_rates(cutoff)
    make_rate_difference_grids(cutoff)

    print(f"\nDone. Outputs under {OUT_BASE}")
    print(f"Defensible time cutoff: {cutoff} Ma "
          f"(see cutoff_time_Ma.txt and well_counts_per_time.csv)")


if __name__ == "__main__":
    main()
