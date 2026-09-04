# Model Data Preparation

> **AI Assistance Notice:** This documentation was produced with the assistance of artificial intelligence (AI) and reviewed by the project author(s).

The `model_data_prep` package prepares CanESM model output for use by the BGC skill-analysis workflows.

It provides a common pipeline for extracting model variables from CMIP6-style archives and related CanESM experiments, selecting vertical levels, combining ensemble members, optionally regridding fields to a regular \(1^\circ \times 1^\circ\) grid, standardizing dimensions, and saving analysis-ready NetCDF files.

The package currently supports:

* DCPP decadal hindcasts and forecasts,
* DCPP assimilation experiments,
* assimilation extensions extracted separately from disk,
* historical and SSP2-4.5 simulations,
* CanOE BGC assimilation experiments, and
* assimilation-sensitivity experiments.

The main user interface is the `DataPrepConfig` configuration object together with the `run_data_prep` pipeline.

## Package structure

| Module            | Purpose                                                                                                                      |
| ----------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| `config.py`       | Defines the `DataPrepConfig` dataclass containing user-configurable options.                                                 |
| `config.yaml`     | Example runtime configuration for the data-preparation pipeline.                                                             |
| `pipeline.py`     | Coordinates the selected preparation workflows for each requested variable.                                                  |
| `hindcast.py`     | Processes DCPP initialized hindcasts and forecasts.                                                                          |
| `assimilation.py` | Processes DCPP assimilation data, assimilation extensions, and sensitivity experiments.                                      |
| `historical.py`   | Processes historical simulations and, when required, their SSP2-4.5 continuation.                                            |
| `canoe.py`        | Processes CanOE BGC assimilation experiments.                                                                                |
| `paths.py`        | Constructs input and output paths used by the package.                                                                       |
| `utils.py`        | Shared utilities for loading, ensemble concatenation, vertical selection, regridding, dimension standardization, and saving. |
| `__main__.py`     | Command-line entry point for running the configured pipeline.                                                                |

## Processing workflow

The data-preparation pipeline follows the general sequence

```text
configuration
    │
    ▼
locate experiment files
    │
    ▼
load requested ensemble members
    │
    ▼
limit vertical extent
    │
    ▼
combine ensemble members
    │
    ▼
optional 1° × 1° regridding
    │
    ▼
standardize dimensions and coordinates
    │
    ▼
save analysis-ready NetCDF files
```

The exact loading procedure depends on the experiment type, but the preprocessing operations are shared wherever possible.

## Configuration

Data preparation is controlled through `DataPrepConfig`, defined in `config.py`.

A typical configuration can be loaded from `config.yaml` and passed to the main pipeline.

```python
from model_data_prep.config import DataPrepConfig
from model_data_prep.pipeline import run_data_prep

cfg = DataPrepConfig(
    user_data_root="/path/to/output",
    var_list=("dissic", "talk"),
    realm_list=("Omon",),
    historical=True,
    assimilation=False,
    hindcast=False,
)

run_data_prep(cfg)
```

The package also provides a YAML-based command-line workflow through `__main__.py`.

### Core configuration

`user_data_root`
: Root directory under which processed data are saved.

`model`
: CMIP model name. The default is `CanESM5`.

`var_list`
: Sequence of variables to process.

`realm_list`
: CMIP realm or table associated with each variable, such as `Omon`.

If only one realm is provided, it is used for every variable. Otherwise, entries in `realm_list` are paired with entries in `var_list`.

`realizations`
: Ensemble realizations to process, for example `r1i1p2f1`. When set to `None`, available realizations are inferred from the archive directory.

`lev_range`
: Maximum vertical coordinate value retained from datasets containing a `lev` dimension. Set to `None` to retain all vertical levels.

`regrid`
: If `True`, fields are bilinearly regridded to a regular global \(1^\circ \times 1^\circ\) grid.

## Selecting experiments

Different workflows are enabled using Boolean configuration options.

### Historical simulations

```yaml
historical: true
historical_initial_year: 1950
historical_final_year: 2025
```

Historical model output is read from the CMIP6 `historical` experiment.

When the requested final year extends beyond the historical experiment, the workflow also reads the corresponding SSP2-4.5 simulation.

When realizations are inferred automatically, the historical workflow retains realizations available in both the historical and SSP2-4.5 experiments when both are required.

### Hindcasts and forecasts

```yaml
hindcast: true
hindcast_initial_year: 1975
hindcast_final_year: 2025
```

Initialized predictions are divided between two DCPP decadal experiments:

* `dcppA-hindcast` for initialization years before 2020;
* `dcppB-forecast` for initialization years from 2020 onward.

Each initialization year is processed independently.

The native forecast time dimension is renamed to `lead_time`, and lead times are represented using one-based integer indices:

```text
1, 2, 3, ..., N
```

The initialization information is stored separately using the `year` coordinate.

Processed initialized predictions use dimensions of the form

```text
ensembles × lead_time × ... × lat × lon
```

where `...` represents optional dimensions such as depth.

### Assimilation

```yaml
assimilation: true
```

Assimilation data are read from the `dcppA-assim` experiment.

Selected realizations are concatenated along the ensemble dimension, optionally regridded, and standardized to use

```text
ensembles × time × ... × lat × lon
```

before being written to disk.

### Assimilation extensions

```yaml
assimilation_extracted_from_disc: true
```

This workflow handles assimilation years that are not contained in the standard CMIP-style archive used by the main assimilation workflow.

The required NetCDF files must first be extracted from the original model output and placed under the expected `extentions` directory.

Files are grouped by year and ensemble member, concatenated through time, vertically limited when requested, and regridded using the same utilities as the primary assimilation workflow.

### CanOE BGC assimilation

```yaml
canoe_assimilation_bgc: true
```

The CanOE workflow processes BGC assimilation experiments stored outside the standard CMIP6 archive.

Individual CanOE assimilation runs are discovered from the configured CanOE root directory, concatenated as ensemble members, vertically subset when needed, optionally regridded, and saved using the standard output conventions.

The output model identifier for these runs is currently

```text
CanESM5-CanOE_1
```

### Assimilation sensitivity experiments

Assimilation-sensitivity runs are handled separately by `process_assimilation_sensitivity`.

The workflow loops over the available BGC sensitivity experiments, reads their variables from the configured sensitivity archive, applies the common vertical-selection and regridding operations, adds an ensemble dimension, and writes each processed experiment to its own output directory.

## Output organization

Processed data are stored beneath `user_data_root`.

The general directory structure is

```text
user_data_root/
└── <variable>/
    └── <collection>/
        └── <model>/
            └── *.nc
```

Examples of collections include

```text
historical
assimilation
forecast
```

Assimilation-sensitivity experiments use experiment-specific directory names.

This organization allows downstream analysis utilities to identify datasets by variable, experiment type, and model.

## Regridding

Regridding is implemented using `xESMF`.

When

```python
cfg.regrid is True
```

the package constructs a regular global \(1^\circ \times 1^\circ\) target grid and applies bilinear interpolation:

```python
xe.Regridder(
    ds,
    ds_out,
    "bilinear",
    ignore_degenerate=True,
    periodic=True,
)
```

Longitude is treated as periodic.

The xESMF output coordinates are subsequently converted from `x` and `y` dimensions into conventional `lon` and `lat` dimensions by `coords_edit`.

If regridding is disabled, the requested variable is returned on its native model grid.

## Ensemble handling

Model realizations can either be specified explicitly or inferred from the archive.

When

```python
cfg.realizations is not None
```

the requested realization IDs are used directly.

When

```python
cfg.realizations is None
```

the package searches the corresponding experiment directory and infers the available realization IDs.

Datasets from individual realizations are concatenated along an ensemble dimension and ultimately standardized to use the dimension name

```text
ensembles
```

for saved products.

## Vertical selection

For three-dimensional ocean variables, `lev_range` can be used to restrict the maximum model depth included during preprocessing.

Conceptually,

```python
ds = ds.where(ds.lev <= lev_range, drop=True)
```

is applied whenever the dataset contains a `lev` dimension.

For example,

```yaml
lev_range: 600
```

retains model levels whose `lev` coordinate is less than or equal to 600.

Setting

```yaml
lev_range: null
```

disables vertical subsetting.

## Shared utilities

The `utils.py` module contains the low-level operations shared among experiment-specific workflows.

### `coords_edit`

Converts the two-dimensional latitude and longitude coordinates generated by xESMF into one-dimensional `lat` and `lon` coordinates and renames the associated `y` and `x` dimensions.

### `available_realizations`

Searches an experiment directory and infers realization identifiers from its subdirectory names.

Some workflows can optionally exclude realizations containing `p1`.

### `resolve_realizations`

Returns explicitly configured realizations when they are provided; otherwise discovers available realizations from disk.

### `limit_level`

Restricts a dataset to model levels at or above the configured maximum `lev` value.

### `concat_members`

Concatenates individually loaded model realizations along an ensemble dimension.

### `open_member_dataset`

Opens all NetCDF files associated with a realization as a single time-concatenated xarray dataset.

### `open_initialized_dataset`

Opens the NetCDF dataset associated with one initialized prediction and applies the configured vertical selection.

### `regrid`

Extracts the requested variable and, when enabled, regrids it to the common global \(1^\circ \times 1^\circ\) grid.

### `standardize_member_time`

Renames the internal member dimension to `ensembles` and places ensemble and time dimensions before the spatial dimensions.

### `save_dataarray`

Creates the output directory when necessary and saves a DataArray as a named NetCDF variable.

### `year_from_time_value`

Extracts the calendar year from an xarray time-coordinate value for use in output filenames.

## Path utilities

`paths.py` centralizes the directory conventions used to navigate the CMIP6 archive.

The main helpers construct locations for

```text
DCPP experiments
historical experiments
SSP2-4.5 experiments
initialized predictions
individual variables
processed output
```

Centralizing these paths prevents the experiment-specific processing modules from duplicating the CMIP6 directory hierarchy.

## Running the pipeline

The intended command-line workflow is:

```bash
python -m model_data_prep
```

The entry point reads the YAML configuration, constructs a `DataPrepConfig`, and passes it to `run_data_prep`.

The pipeline loops over each configured variable, determines its corresponding realm, and executes the enabled experiment workflows.

For example, a configuration requesting historical and assimilation data will conceptually execute

```text
variable 1
    ├── historical
    └── assimilation

variable 2
    ├── historical
    └── assimilation
```

before continuing to the next variable.

## Dependencies

The data-preparation package relies primarily on

```text
numpy
xarray
xESMF
PyYAML
tqdm
```

and the NetCDF/Dask dependencies required by `xarray.open_mfdataset`.

## Notes

The default archive locations in `DataPrepConfig` are specific to the ECCC/CCCma computing environment.

Users running the package elsewhere must provide appropriate values for the CMIP6 archive root, CanOE experiment root, assimilation-sensitivity root, and output root.

The preparation routines also assume the CMIP6 directory conventions used by the CanESM archives. Alternative archives may require changes to the path-construction utilities.
