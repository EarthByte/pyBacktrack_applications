# Input data

Everything in this folder is committed to the repo so the figures are
reproducible from a fresh clone (with the two exceptions noted at the
bottom — the global paleobathymetry grids, which are too large for
GitHub). The contents of `data/` are written by
`../tools/populate_data.py`; see that script's docstring for the
originating-on-disk locations.

```
data/
├── README.md      <- this file
├── wells/         <- 109 NW Shelf wells (.txt / .dat)
└── grids/
    ├── vgg_nwshelf.nc                          <- regional Sandwell VGG
    └── dynamic_topography/D10_gmcm9/<t>.00.nc  <- regional D10_gmcm9 slices
```

---

## `wells/` — 109 NW Shelf wells

Lithology-format well files consumed by `pybacktrack.backstrip_well` and
`pybacktrack.paleo_bathymetry`. Each file is a plain-text table with a
header block of metadata (`SiteLongitude`, `SiteLatitude`,
`RiftStartAge`, `RiftEndAge`, ...) followed by one row per stratigraphic
interval.

Origin
:   Australian Geoscience Information Network (AGIN) public-domain
    petroleum-exploration releases, redistributed with the pyBacktrack 1.0
    release as the bundled NW Shelf backstripping example.

Featured well
:   `asteras.txt` (124.12 °E, 13.15 °S; 20 stratigraphic intervals; oldest
    horizon 146 Ma; youngest 24 Ma). Used as the single-well comprehensive
    example in Figs 6–8.

Number of files
:   109 (`.txt` + `.dat`).

Licence
:   Same terms as the pyBacktrack distribution under which these files
    are bundled (see https://pybacktrack.readthedocs.io).

Citation
:   Müller, R. D., Cannon, J., Williams, S., and Dutkiewicz, A. (2018).
    *PyBacktrack 1.0: A Tool for Reconstructing Paleobathymetry on
    Oceanic and Continental Crust.* Geochemistry, Geophysics,
    Geosystems 19, 1898–1909. https://doi.org/10.1029/2017GC007313

---

## `grids/vgg_nwshelf.nc` — regional vertical gravity gradient

A regional crop of the Sandwell & Smith vertical-gravity-gradient (VGG)
release V31.1, cut to the Fig 5 window `[114°E, 130°E, 20°S, 9°S]`.
Units: Eötvös.

Origin
:   `/Users/dietmar/grids/vgg_31.1.nc` (the global Sandwell V31.1
    release). Cut by `tools/populate_data.py` via `pygmt.grdcut` to the
    region above.

Resolution
:   ~1 arc-minute (matches the source release).

Use
:   Basemap for Fig 5 (NW Shelf well-location map). The diverging
    Crameri `vik` colourmap, capped at ±60 E, exposes the NE–SW
    Mesozoic rift trends, depocentre boundaries and transform-margin
    lineaments that motivate the distribution of NW Shelf wells.

Licence
:   Public domain (data product of Scripps Institution of Oceanography).
    Please credit Sandwell & Smith.

Citation
:   Sandwell, D. T., Müller, R. D., Smith, W. H. F., Garcia, E., and
    Francis, R. (2014). *New global marine gravity model from CryoSat-2
    and Jason-1 reveals buried tectonic structure.* Science 346, 65–67.
    Version 31.1 release: https://topex.ucsd.edu/marine_grav

---

## `grids/dynamic_topography/D10_gmcm9/` — regional dynamic topography

Per-time NetCDF slices of the Braz et al. (2021) `D10_gmcm9`
dynamic-topography model, region-cropped to `[109°E, 136°E, 26°S, 4°S]`
(the Fig 12 plotting window `[113, 132, -22, -8]` plus a 4° edge buffer
to keep the GMT `surface` reconstruction free of edge ringing).

Coverage
:   0–150 Ma in 5 Myr time slices (31 files in total). Slice filenames
    match the upstream convention (`<time>.00.nc`) so the time-interp
    logic in `12_dt_maps.py` is unchanged.

Origin
:   pybacktrack bundled data,
    `pybacktrack/bundle_data/dynamic_topography/models/Braz2021/D10_gmcm9/`.
    Cut by `tools/populate_data.py` via `pygmt.grdcut`.

Sign convention
:   Positive = location was higher than today (uplifted); negative =
    location was lower than today (subsided down); zero = same as
    present. Matches the global D10_gmcm9 release.

Use
:   Raw dynamic-topography input for Fig 12 and (interpolated to 1 Myr
    spacing) the optional `12a_dt_video.py` animation. Also the dynamic-
    topography correction term in the backstripping configurations C and
    SL+DT used by Figs 7, 9, 10, 11.

Citation
:   Braz, C., Flament, N., and Müller, R. D. (2021). *Long-wavelength
    dynamic topography from upper-mantle convection models constrained
    by 130 Ma of palaeogeography.* Earth and Planetary Science Letters
    574, 117174. https://doi.org/10.1016/j.epsl.2021.117174

---

## External: global paleobathymetry grids

The Zahirovic et al. (2022) GDH1 paleobathymetry release is needed by
Figs 2, 3 and 4 but is too large (~200–300 MB per reference frame) to
commit to GitHub. Both reference frames (optimised mantle and
paleomagnetic) are fetched on demand by:

```bash
cd ../scripts
./fetch_paleobathymetry_grids.sh mantle           # Z22 mantle frame
./fetch_paleobathymetry_grids.sh paleomagnetic    # Z22 paleomagnetic frame
```

The downloaded grids land under
`data/grids/paleobathymetry/Zahirovic2022_mantle/` or
`data/grids/paleobathymetry/Zahirovic2022_PMag/` (the fallback
resolution paths used by `config.py`). Toggle which frame the figures
consume via `PALEO_BATHY_SOURCE` in `scripts/config.py`
(`"Z22_mantle"` or `"Z22_PMag"`).

Citation
:   Zahirovic, S., Eleish, A., Doss, S., Pall, J., Cannon, J., Müller,
    R. D., Williams, S., Tetley, M., Hassan, R., Bahlburg, H., Boggiani,
    P., and Fox, P. A. (2022). *Subduction and carbonate platform
    interactions.* Geoscience Data Journal 9, 109–128.

---

## External: Straume 2020 paleogeography grids

Required only by Fig 4 (the pyBacktrack 1.5 vs Straume 2020 comparison).
Public Zenodo / Pangaea release; the figure script will error
informatively if the grids are not present under
`Straume_paleogeography_grids/` next to the repo. See
`04_pybacktrack_vs_straume_comparison.py` for the expected filename
pattern.

Citation
:   Straume, E. O., Gaina, C., Medvedev, S., Hochmuth, K., Gohl, K.,
    Whittaker, J. M., Abdul Fattah, R., Doornenbal, J. C., and Hopper,
    J. R. (2020). *GlobSed: Updated total sediment thickness in the
    world's oceans.* Geochemistry, Geophysics, Geosystems 21,
    e2020GC009107. (Paleogeography grids released in the same project
    series.)
