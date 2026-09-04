# Configuration

> **AI Assistance Notice:** This documentation was produced with the assistance of artificial intelligence (AI) and reviewed by the project author(s).

The `configs` package contains shared configuration data used throughout the BGC skill analysis and plotting workflows. It centralizes variable metadata, unit conversions, geographical regions, plotting conventions, climatological data locations, ENSO event definitions, and CanESM ocean vertical-coordinate information.

Keeping these settings outside the analysis code avoids duplicating project-specific constants across notebooks and modules and provides a common set of conventions for loading, processing, and visualizing model and observational data.

## Usage

Several commonly used configuration files are loaded automatically when the package is imported:

```python
import configs

configs.units
configs.varx
configs.unit_change_dics
configs.boundaries_dict
configs.ENSO
configs.contourf
configs.pcolor_ranges
configs.model_levels
configs.model_lev_bounds
```

Alternatively, individual configuration objects can be imported directly:

```python
from configs import units, boundaries_dict, model_levels
```

Configuration files that are not loaded in `configs.__init__`, such as climatology paths and plotting styles, can be read explicitly when required.

## Configuration files

| File                                 | Purpose                                                                                                |
| ------------------------------------ | ------------------------------------------------------------------------------------------------------ |
| `units.yaml`                         | Display units associated with physical and biogeochemical variables.                                   |
| `unit_changes.yaml`                  | Multiplicative conversion factors used to convert source data to the units adopted by the analysis.    |
| `varx.yaml`                          | Maps alternative or derived variable names to their corresponding canonical variable names.            |
| `regions.yaml`                       | Latitude and longitude boundaries for commonly analysed geographical regions.                          |
| `ENSO.yaml`                          | Timing information used to identify observed ENSO events.                                              |
| `contourf_levels.yaml`               | Explicit contour levels used for selected variables in filled-contour plots.                           |
| `pcolor_ranges.yaml`                 | Recommended plotting ranges for variables, including absolute fields and anomaly/difference fields.    |
| `styles.yaml`                        | Matplotlib marker, line-style, and transparency conventions for different experiment or dataset types. |
| `mplRC.py`                           | Utilities and predefined Matplotlib `rcParams` settings used to maintain consistent figure formatting. |
| `GLODAP_climatology.yaml`            | File location and variable-name mapping for the GLODAP mapped climatology.                             |
| `model_climatology.yaml`             | Locations of model-derived climatological datasets used by selected diagnostics.                       |
| `CanESM_ocean_levels.json`           | Representative depths of the CanESM ocean model levels.                                                |
| `CanESM_ocean_level_boundaries.json` | Vertical boundaries associated with the CanESM ocean model levels.                                     |

## Variable names and units

### `units.yaml`

`units.yaml` maps variable names to the units used when displaying or labelling results.

For example:

```yaml
dissic: '$\mu$mol kg$^{-1}$'
thetao: '$^{o}$ C'
mlotst: 'm'
```

The values are primarily intended for plotting labels and therefore may contain Matplotlib-compatible LaTeX notation.

The keys generally correspond to the variable names used internally by the analysis library rather than necessarily the native names contained in every source dataset.

### `varx.yaml`

Some related variables use different names depending on their source or preprocessing state. `varx.yaml` provides mappings from these alternative names to canonical variable names.

For example:

```yaml
ntalk: talk
ndissic: dissic
no3os: no3
```

This mapping allows downstream code to recover common metadata, such as units or plotting conventions, for related variables.

### `unit_changes.yaml`

`unit_changes.yaml` contains multiplicative conversion factors used to convert variables from their source units into the units adopted by the analysis.

A conversion is generally applied as

```python
converted = original * unit_change_dics[var]
```

For variables that require no conversion, the conversion factor is `1`.

For example, carbonate-system variables such as dissolved inorganic carbon and alkalinity are converted using factors appropriate for their source units, while temperature and salinity currently use unity factors.

These values describe project-specific conversions and should be updated if the units of an input dataset change.

## Geographical regions

### `regions.yaml`

`regions.yaml` defines frequently used geographical domains using latitude and longitude limits:

```yaml
NEP:
    lat_min: 40
    lat_max: 65
    lon_min: -160
    lon_max: -110
```

Each region contains

* `lat_min`: southern latitude boundary,
* `lat_max`: northern latitude boundary,
* `lon_min`: western longitude boundary,
* `lon_max`: eastern longitude boundary.

The configuration includes regional domains such as the Northeast Pacific (`NEP`) and Northwest Atlantic (`NWA`), hemispheric domains, and equatorial Pacific and Atlantic regions.

Some regions cross the longitude discontinuity. For example, the equatorial Pacific configuration may have `lon_min > lon_max`. Code selecting such regions must therefore support longitude wrapping rather than assuming that `lon_min < lon_max`.

## ENSO configuration

### `ENSO.yaml`

`ENSO.yaml` stores observed temporal information used by ENSO-related analyses.

Event boundaries are represented as fractional years, where the fractional component represents the month within a year. These values are used by the ENSO-analysis utilities to identify periods associated with specified events.

The configuration is intended to keep ENSO event definitions independent from the analysis implementation so that event selections can be modified without changing the analysis code.

## Plot configuration

Three configuration files control different aspects of plotting.

### `contourf_levels.yaml`

Contains explicit contour values for variables for which fixed contour intervals are useful.

For example:

```yaml
uo:
    [-1.0, -0.8, -0.6, ..., 0.8, 1.0]
```

Using common contour levels makes figures from different experiments or datasets visually comparable.

### `pcolor_ranges.yaml`

Contains preferred minimum and maximum plotting limits.

Ranges are grouped according to the type of quantity being visualized. For example, separate ranges are provided for absolute fields and difference/anomaly fields and, for some variables, for surface and subsurface values.

A typical entry has the form

```yaml
thetao:
    surface:
        vmin: -2
        vmax: 2
    depth:
        vmin: -1
        vmax: 1
```

These ranges provide project defaults and can be overridden by individual plotting functions when required.

### `styles.yaml`

Defines common plotting styles for dataset or experiment categories such as observations, assimilation runs, historical simulations, hindcasts, and control simulations.

Style entries may specify properties such as

```yaml
marker
linestyle
alpha
```

Centralizing these choices gives the different notebooks and plotting utilities a consistent visual language.

### `mplRC.py`

`mplRC.py` contains Matplotlib configuration utilities used to establish common figure formatting.

The `setRC` function applies attributes from an RC configuration object to Matplotlib settings such as

* figure size and resolution,
* font sizes,
* axis labels and titles,
* legend formatting,
* line widths, and
* marker sizes.

This provides a common figure style without repeating `matplotlib.rc` calls throughout analysis notebooks.

## Climatology configuration

### `GLODAP_climatology.yaml`

Defines the location of the GLODAP mapped climatological dataset and maps the internal variable names used by this repository to their corresponding GLODAP names.

For example:

```yaml
rename_dict:
    dissic: TCO2
    talk: TAlk
    no3: NO3
```

The mapping allows observational climatologies to be loaded using the same internal variable naming conventions as model data.

### `model_climatology.yaml`

Contains locations of precomputed model climatology products required by specific diagnostics.

Unlike general variable metadata, these entries are environment-dependent filesystem paths and may need to be changed when the analysis is run on another system.

## CanESM ocean vertical coordinates

The two JSON files

```text
CanESM_ocean_levels.json
CanESM_ocean_level_boundaries.json
```

describe the vertical ocean grid used by the CanESM datasets analysed in this repository.

`CanESM_ocean_levels.json` contains representative model-level depths, while `CanESM_ocean_level_boundaries.json` contains the corresponding vertical interfaces.

When `configs` is imported, these values are loaded as NumPy arrays:

```python
from configs import model_levels, model_lev_bounds
```

These arrays can be used when selecting model depths, constructing vertical masks, or calculating quantities that depend on layer thickness.

## Adding new configuration

When introducing a new variable, check whether entries should also be added to:

```text
units.yaml
unit_changes.yaml
varx.yaml
contourf_levels.yaml
pcolor_ranges.yaml
```

Not every variable requires entries in every file. For example, an alias is only needed in `varx.yaml` when the variable name differs from its canonical representation, and explicit plotting limits are only required when fixed defaults are useful.

When introducing a new geographical domain, add the latitude and longitude boundaries to `regions.yaml`.

Configuration values should remain data-oriented. Analysis logic should generally remain in the corresponding analysis or plotting modules rather than being embedded in YAML files.

