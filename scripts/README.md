# scripts/

One Python module per figure, two shared helpers, and three shell
runners.

## Layout

```
scripts/
├── README.md                                   <- this file
├── config.py                                   <- shared paths + constants
├── paleobathy_render.py                        <- shared rendering pipeline (Figs 2, 3)
│
├── 01_workflow_flowchart.py                    <- Fig 1
├── 02_north_atlantic_paleobathymetry.py        <- Fig 2 (statics)
├── 02_north_atlantic_paleobathymetry_video.py  <- Fig 2 (MP4)
├── 03_southern_ocean_paleobathymetry.py        <- Fig 3 (statics)
├── 03a_southern_ocean_paleobathymetry_video.py <- Fig 3 (MP4)
├── 03b_southern_ocean_paleobathymetry_video_Z22.py
├── 04_pybacktrack_vs_straume_comparison.py     <- Fig 4
├── 05_nwshelf_map.py                           <- Fig 5
├── 06_sl_dt_curves.py                          <- Fig 6
├── 07_backstrip_wells.py                       <- Fig 7
├── 07a_backstrip_all_nwshelf.py                <- DATA STEP (writes rate grids)
├── 08_geohistory_wells.py                      <- Fig 8
├── 09_rate_timeseries.py                       <- Fig 9
├── 10_rate_maps.py                             <- Fig 10
├── 10a_make_subsidence_videos.py               <- Fig 10 MP4 companion
├── 11_rate_boxplots.py                         <- Fig 11
├── 12_dt_maps.py                               <- Fig 12
├── 12a_dt_video.py                             <- Fig 12 MP4 companion
│
├── run_all_figures.sh                          <- batch runner (every Fig)
├── make_all_videos.sh                          <- batch runner (every MP4)
└── fetch_paleobathymetry_grids.sh              <- one-time external grid fetch
```

## Figure catalogue

| Fig. | Script | Inputs | What it shows |
|---|---|---|---|
| 1 | `01_workflow_flowchart.py` | none | Paleobathymetry-gridding workflow (matplotlib) |
| 2 | `02_north_atlantic_paleobathymetry.py` | paleobathy grids, plate model | LAEA centred on mid-N-Atlantic, 4 panels |
| 3 | `03_southern_ocean_paleobathymetry.py` | paleobathy grids, plate model | LAEA south-polar, 4 panels, present-day acronyms |
| 4 | `04_pybacktrack_vs_straume_comparison.py` | paleobathy + Straume 2020 grids | Time-series median ± MAD comparison |
| 5 | `05_nwshelf_map.py` | `data/grids/vgg_nwshelf.nc`, wells | Well-location map on regional VGG basemap |
| 6 | `06_sl_dt_curves.py` | bundled SL + DT models | Asteras-well sea-level & DT time series |
| 7 | `07_backstrip_wells.py` | asteras well | Backstripped tectonic subsidence (3 configs) |
| 8 | `08_geohistory_wells.py` | asteras well | Geohistory / Wheeler diagram |
| 9 | `09_rate_timeseries.py` | rate grids (data step) | Margin-scale time series A, B, C = B − A |
| 10 | `10_rate_maps.py` | rate grids (data step) | 3 × 4 rate-map matrix |
| 11 | `11_rate_boxplots.py` | rate grids (data step) | 5-Myr-binned boxplots, 3 panels |
| 12 | `12_dt_maps.py` | `data/grids/dynamic_topography/D10_gmcm9/`, picked times | Raw DT at same four times as Fig 10 |

## The data step (`07a_backstrip_all_nwshelf.py`)

Runs `pybacktrack.backstrip_well` over every well in `data/wells/`,
computes per-Myr subsidence rates for configurations A (no SL/no DT)
and B (SL + DT), gridds them with `blockmedian` → `surface`, and writes:

* `figures/output/nwshelf_subsidence/well_summary.csv`
* `figures/output/nwshelf_subsidence/all_wells_locations.txt`
* `figures/output/nwshelf_subsidence/sub_<well>.csv` (per-well subsidence)
* `figures/output/nwshelf_subsidence/rate_<well>.csv` (per-well rate)
* `figures/output/nwshelf_subsidence/rate_<time>.nc` (per-time grids, 0–150 Ma)
* `figures/output/nwshelf_subsidence/rate_delta_stats.csv`
* `figures/output/nwshelf_subsidence/rate_delta_picked_times.txt`
* `figures/output/nwshelf_subsidence/cutoff_time_Ma.txt`

The data step is invoked automatically by the `fig09` target. Re-run
manually to refresh the per-Myr grids:

```bash
python 07a_backstrip_all_nwshelf.py
```

## Dependency graph

```
[wells]──┬──────────────────────► fig07, fig08
         └─► 07a_data_step ──┬─► fig09 ─► fig10 ──► fig12
                             ├─► fig11
                             └─► 10a_videos
```

Fig 10 writes `rate_delta_picked_times.txt`, which Fig 12 (and the
Fig 11 highlighting logic) consume. Run Fig 10 before Fig 12 the first
time through.

## Runner targets

```bash
./run_all_figures.sh                 # every static figure (01..12)
./run_all_figures.sh fig09           # just Fig 9 (with data step if needed)
./run_all_figures.sh fig10 fig11     # multiple targets

./make_all_videos.sh                 # every MP4 animation
./make_all_videos.sh 10a 12a         # just NW Shelf rate + DT MP4s
./make_all_videos.sh --no-force      # honour cached frames
```

Add `--skip-stitch` to `make_all_videos.sh` to halt before the `ffmpeg`
mp4 assembly step (useful when you only want the per-frame PNGs).

## Key constants (in `config.py`)

| Constant | Default | Notes |
|---|---|---|
| `PALEO_BATHY_SOURCE` | `"Z22_mantle"` | Toggle to `"Z22_PMag"` to switch grid + frame |
| `DYNAMIC_TOPOGRAPHY_MODEL` | `"D10_gmcm9"` | Braz et al. (2021), 0–150 Ma |
| `SEA_LEVEL_MODEL` | `"Haq2024_Hybrid_SealevelCurve_Longterm"` | Bundled long-term hybrid |
| `MAX_ANALYSIS_AGE` | `150.0` | Upper bound on every NW Shelf x-axis |
| `WELLS` | `[{name: "Asteras", ...}]` | Single featured well for Figs 6–8 |
| `VGG_GRID` | `data/grids/vgg_nwshelf.nc` | Regional Sandwell cut |
| `DT_GRIDS_DIR` | `data/grids/dynamic_topography/D10_gmcm9` | Regional Braz 2021 cut |

## Conventions in force

* Subsidence rate sign: **positive = subsiding (deeper)**, negative =
  uplifting.
* Time-axis direction: every x-y plot runs LEFT to RIGHT from 150 Ma to
  0 Ma (use `ax.set_xlim(MAX_ANALYSIS_AGE, 0)`).
* Plotting libraries: pyGMT for everything 1D plus the NW Shelf maps;
  matplotlib + cartopy + gplately for the global / south-pole paleobathy
  views (Figs 2, 3).
* CPT defaults: `dem1` (Crameri) for one-sided rate maps with the
  above-range clamped to the 50 m/Myr colour; `vik` (Crameri) for
  differences and the VGG basemap; `nuuk` for paleobathymetry; `cool`
  reversed for dynamic topography.
