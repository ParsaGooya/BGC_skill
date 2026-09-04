# Observational Data Preparation

> **AI Assistance Notice:** This documentation was produced with the assistance of artificial intelligence (AI) and reviewed by the project author(s).

The `obs_data_prep` directory contains dataset-specific workflows for downloading, extracting, transforming, and standardizing observational and reanalysis products used in the BGC skill analysis.

Unlike `model_data_prep`, which provides a common processing pipeline for multiple CanESM experiments, the observational workflows are organized separately by source dataset because each product has its own native format, variable names, spatial grid, temporal structure, and preprocessing requirements.

The currently supported observational products are:

* **ESA CCI / related chlorophyll products** for surface chlorophyll concentration;
* **GLODAP** bottle observations and mapped climatologies for ocean biogeochemical variables;
* **SODA** ocean reanalysis fields used for physical-ocean validation.

## Directory structure

```text
obs_data_prep/
├── ESACCI_chlos/
│   ├── data_links.txt
│   ├── download_data.py
│   └── regrid_chlos_ESACCI_or_Shujie23.py
│
├── GLODAP/
│   ├── extract_GLODAP.ipynb
│   ├── fill_data_w_climatology.ipynb
│   └── salinity_normalized.ipynb
│
└── SODA/
    ├── download_soda.py
    └── soda_data_regrid.py
```

Each subdirectory represents an independent preprocessing workflow.

## General workflow

Although implementation details differ between datasets, observational preparation generally follows:

```text
raw observational source
        │
        ▼
download or load source data
        │
        ▼
extract required variables
        │
        ▼
rename variables and coordinates
        │
        ▼
quality-control / missing-value handling
        │
        ▼
optional derived-variable calculation
        │
        ▼
spatial aggregation or regridding
        │
        ▼
save in project-standard format
```

The resulting files are organized using variable names consistent with the model-analysis workflows wherever possible.

---

# ESA CCI chlorophyll

The `ESACCI_chlos` workflow prepares remotely sensed surface chlorophyll concentration.

It contains two main stages:

```text
download_data.py
        │
        ▼
raw monthly satellite files
        │
        ▼
regrid_chlos_ESACCI_or_Shujie23.py
        │
        ▼
monthly 1° × 1° chlorophyll dataset
```

## Downloading chlorophyll data

`download_data.py` reads a collection of source URLs and downloads monthly NetCDF files.

The script includes logic for correcting outdated Ocean Colour CCI download links before requesting the corresponding files.

The URL conversion currently updates source identifiers associated with older product versions to the product version expected by the workflow.

Already-downloaded files are detected and skipped.

### Output

Raw monthly chlorophyll files are written to the configured observational-data directory.

The download location is currently specified directly in the script and must therefore be changed when running the workflow in another environment.

## Chlorophyll regridding

`regrid_chlos_ESACCI_or_Shujie23.py` prepares chlorophyll data from either

```text
ESACCI
```

or

```text
Shujie23
```

depending on the configured `source`.

The source variable is expected to be

```python
chlor_a
```

and is converted to the project variable name

```text
chlos
```

during output.

### Spatial aggregation

The native high-resolution chlorophyll grid is reduced to approximately \(1^\circ \times 1^\circ\) resolution by averaging groups of grid cells in the latitude and longitude directions.

For the high-resolution input used by this workflow, blocks of 24 native grid cells are averaged along each horizontal dimension.

The resulting grid uses nominal cell centers at

```text
latitude:  -89.5, -88.5, ..., 89.5
longitude: -179.5, -178.5, ..., 179.5
```

### Time handling

For products whose timestamp is encoded in the filename, the year and month are extracted and assigned as an xarray time coordinate.

Processed monthly files are first written individually and are subsequently concatenated along `time`.

The final product is sorted chronologically and saved as a single NetCDF file containing the full available time range.

---

# GLODAP

The `GLODAP` directory prepares hydrographic observations from the GLODAP merged bottle-data product and related mapped climatologies.

The workflow has three distinct purposes:

```text
extract_GLODAP.ipynb
        │
        ├── standard GLODAP variables
        │
        ├── model-compatible names
        │
        ▼
CSV observational records

salinity_normalized.ipynb
        │
        ▼
salinity-normalized carbonate variables

fill_data_w_climatology.ipynb
        │
        ▼
selected missing measurements filled from
GLODAP mapped climatology
```

## Extracting GLODAP observations

`extract_GLODAP.ipynb` reads the GLODAP merged master CSV and creates smaller, variable-specific observational files.

The notebook currently targets variables including:

```text
silicate
phosphate
nitrite
oxygen
talk
theta
salinity
tco2
```

Additional variables can be included through `var_list`.

For each requested variable, the workflow retains the observation metadata required by the downstream analysis:

```text
year
month
day
hour
minute
latitude
longitude
depth
pressure
variable value
```

A full datetime coordinate is also constructed from the observation date and time.

### Missing values

The GLODAP missing-value indicator

```text
-9999
```

is converted to `NaN`.

### Coordinate standardization

Several columns are renamed to match conventions used elsewhere in the project:

```text
month     → time
latitude  → lat
longitude → lon
depth     → lev
```

### Variable-name standardization

Where possible, GLODAP variable names are converted to their model/CMIP-style equivalents.

Examples include:

```text
salinity  → so
theta     → thetao
oxygen    → o2
tco2      → dissic
nitrate   → no3
phosphate → po4
```

The resulting variable-specific tables are saved as CSV files beneath the corresponding observational-data directory.

The current GLODAP extraction notebook uses the GLODAP v2.2023 merged master file and produces records covering the observational years available in that source dataset.

## Filling missing values with mapped climatology

`fill_data_w_climatology.ipynb` provides an optional preprocessing step for variables whose observational records contain missing measurements.

The workflow uses the GLODAP v2.2016b mapped climatology as a spatial and vertical reference.

For each missing observational value:

1. the nearest valid climatological horizontal grid point is identified;
2. the nearest climatological depth is selected;
3. the corresponding climatological value is retrieved;
4. the missing observational value is replaced with that value.

Longitude coordinates are first normalized so that the observational and climatological grids use compatible longitude conventions.

The notebook currently applies this procedure to selected nutrient variables such as

```text
silicate
NO3
PO4
```

The operation should be interpreted as **gap filling using a climatological estimate**, not as an observed measurement. Downstream analyses should retain this distinction when interpreting filled records.

## Salinity-normalized carbonate variables

`salinity_normalized.ipynb` computes salinity-normalized GLODAP carbonate-system variables.

The current workflow applies this normalization to

```text
talk
tco2
```

using simultaneously extracted GLODAP salinity.

For a variable \(X\), the normalized value is calculated as

$$
X_{\mathrm{norm}}
=
X
\frac{35}{S},
$$

where \(S\) is the corresponding observed salinity.

The reference salinity is therefore

$$
S_{\mathrm{ref}} = 35.
$$

The resulting variables are stored with an `n` prefix, for example:

```text
ntalk
ntco2
```

These normalized quantities are saved independently from the original GLODAP measurements.

---

# SODA

The `SODA` directory prepares ocean physical variables from the SODA reanalysis.

The workflow consists of:

```text
download_soda.py
        │
        ▼
raw SODA NetCDF files
        │
        ▼
soda_data_regrid.py
        │
        ▼
standardized 1° × 1° fields
```

## Downloading SODA

`download_soda.py` identifies available NetCDF files from the SODA data server.

The helper

```python
find_netcdf_links(url)
```

retrieves the dataset webpage, extracts links ending in `.nc`, and resolves relative paths to full download URLs.

The script then selects monthly ocean files identified by

```text
_mn_ocean
```

and downloads files that are not already present locally.

The current workflow targets SODA version 3.15.2.

## Variable preparation

`soda_data_regrid.py` processes several physical-ocean fields and renames them to the project conventions:

```text
temp → thetao
salt → so
wt   → wo
u    → uo
v    → vo
```

This allows the resulting SODA data to be used alongside CanESM fields without additional variable-name translation.

## Dimension standardization

Native SODA dimension names depend on the variable grid.

Vertical coordinates are converted from either

```text
st_ocean
```

or

```text
sw_ocean
```

to

```text
lev
```

and horizontal dimensions are converted from either

```text
yt_ocean, xt_ocean
```

or

```text
yu_ocean, xu_ocean
```

to

```text
lat, lon
```

as appropriate.

## Spatial aggregation

The workflow creates a regular \(1^\circ\) grid by averaging source values falling within each one-degree longitude and latitude interval.

Longitude bins span

```text
0°–360°
```

during aggregation.

The resulting longitudes are subsequently converted to the

```text
-180°–180°
```

convention using `lonfixer`.

The processed fields are saved as NetCDF files under variable-specific observational directories.

---

# Output conventions

The observational workflows aim to produce data that can be consumed by the analysis modules with minimal dataset-specific handling.

Processed products generally follow the directory structure

```text
<data_root>/
└── <variable>/
    └── observation/
        └── <dataset>/
            └── processed files
```

Examples include

```text
dissic/observation/GLODAP/
thetao/observation/SODA/
chlos/observations/ESACCI/
```

Variable names are standardized to match the model data where practical.

Common standardized names include:

| Quantity                   | Project name |
| -------------------------- | ------------ |
| Potential temperature      | `thetao`     |
| Salinity                   | `so`         |
| Dissolved inorganic carbon | `dissic`     |
| Alkalinity                 | `talk`       |
| Oxygen                     | `o2`         |
| Nitrate                    | `no3`        |
| Phosphate                  | `po4`        |
| Zonal velocity             | `uo`         |
| Meridional velocity        | `vo`         |
| Vertical velocity          | `wo`         |
| Surface chlorophyll        | `chlos`      |

---

# Dependencies

The observational preparation workflows use combinations of:

```text
numpy
pandas
xarray
xESMF
requests
BeautifulSoup
wget
tqdm
```

Not all dependencies are required for every dataset.

---

# Environment-specific paths

Many input and output directories are currently specified directly inside the scripts and notebooks and correspond to the ECCC/CCCma computing environment.

For example, several workflows refer directly to paths under

```text
/space/hall5/...
```

or

```text
/space/hall7/...
```

These paths must be changed when reproducing the workflows in another computing environment.

Unlike `model_data_prep`, `obs_data_prep` does not currently provide a shared configuration object or command-line pipeline.

---

# Reproducibility notes

The observational workflows are tied to particular source-product versions.

At present these include products such as:

```text
GLODAP v2.2023 merged observations
GLODAP v2.2016b mapped climatology
SODA 3.15.2
ESA CCI / associated chlorophyll products
```

When updating an observational product, users should verify:

* source variable names;
* missing-value conventions;
* coordinate names;
* spatial resolution;
* longitude convention;
* temporal coverage;
* units; and
* file naming conventions.

Changes in any of these may require updates to the corresponding preprocessing workflow.

