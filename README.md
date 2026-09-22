# BGC_skill

**BGC_skill** is a Python-based analysis framework for evaluating and comparing **ocean biogeochemical (BGC) and physical properties across climate-model simulations and observational products**.

The repository provides workflows for preparing model and observational datasets, calculating derived quantities and diagnostics, comparing model simulations against observations, and visualizing spatial and temporal variability. Particular emphasis is placed on evaluating **CanESM5 and CanESM5-CanOE**, including historical simulations, assimilation experiments, and seasonal-to-decadal prediction systems. However, the design goal was to keep hardcoding minimal to allow adabtability to future needs and datasets with minimal effort. The hope is for this package to be easily adaptable for when new model versions and runs, obervational datasets or analysis is available.

The analysis tools cover global and regional climatologies, time series, vertical structure, model–observation comparisons, ENSO-related variability, and tropical ocean dynamics and biogeochemistry.

> **AI Assistance Notice:** Documentation in this repository was produced with the assistance of artificial intelligence (AI) and reviewed by the project author(s).

---

## Repository Overview

The repository is organized around three main components:

1. **Data preparation** – Prepare model and observational datasets in consistent formats for subsequent analysis.
2. **Analysis utilities** – Shared routines for loading, processing, masking, aggregating, and comparing datasets.
3. **Analysis notebooks** – Configurable workflows for specific scientific analyses and visualizations.

A typical analysis follows:

**Configuration → Data loading → Data preparation → Analysis → Visualization**

Most notebooks expose their main scientific choices near the beginning of the notebook so that variables, experiments, observational products, time periods, regions, depths, and optional calculations can be changed without modifying the underlying analysis code.

---

## Repository Structure

```text
BGC_skill/
│
├── Configs/
│   ├── units.yaml
│   ├── unit_changes.yaml
│   ├── varx.yaml
│   ├── regions.yaml
│   ├── ENSO.yaml
│   ├── contourf_levels.yaml
│   ├── pcolor_ranges.yaml
│   ├── styles.yaml
│   ├── mplRC.py
│   ├── GLODAP_climatology.yaml
│   ├── model_climatology.yaml
│   └── ...
│
├── Data_prep/ 
│    ├── model_data_prep/
│    |
│    ├── obs_data_prep/
│          ├── ESACCI_chlos/
│          ├── GLODAP/
│          └── SODA/
│    
│
├── Modules/
│
├── Notebooks/
│        ├── global_map_notebooks/
│        │
│        ├── timeseries_notebooks/
│        │
│        ├── GLODAP_notebooks/
│        │
│        ├── tropics_notebooks/
│
└── README.md
```

Each major package or workflow contains additional documentation describing its configuration and intended use.

---

# Configuration

The `Configs/` package provides shared configuration information used throughout the repository.

It centralizes definitions such as:

* variable names and aliases;
* plotting units;
* unit conversions;
* geographic regions;
* ENSO observational temporal boundaries;
* plotting ranges and contour levels;
* plotting styles;
* CanESM5 model vertical levels and level boundaries;
* observational and model climatology locations.

Centralizing these definitions ensures that notebooks and analysis routines use consistent variable naming, units, regions, and visualization settings.

See the `configs/` README for detailed descriptions of individual configuration files.

---

# Model Data Preparation

The `model_data_prep/` package prepares climate-model output for use by the analysis workflows.

The package supports several types of CanESM experiments, including:

* **DCPP hindcasts and forecasts**
* **DCPP assimilation**
* **assimilation extensions**
* **historical and SSP2-4.5 simulations**
* **CanESM5-CanOE biogeochemical assimilation**
* **assimilation sensitivity experiments**

The data-preparation pipeline handles tasks such as:

* discovering available ensemble realizations;
* selecting requested variables and vertical levels;
* concatenating ensemble members;
* standardizing dimensions and coordinates;
* converting prediction time axes to initialization and lead-time representations;
* regridding fields to a common spatial grid;
* writing processed datasets for subsequent analysis.

Experiment-specific processing is separated from shared utilities so that new model experiments can be incorporated without duplicating the common data-processing logic.

See `model_data_prep/README.md` for configuration and workflow details.

---

# Observational Data Preparation

The `obs_data_prep/` directory contains dataset-specific workflows for preparing observational and reanalysis products used by the analysis notebooks. The expectation is that all data prep for future observational sources
be included here.

Currently supported preparation workflows include:

### ESA CCI Chlorophyll

Tools are provided for downloading and spatially aggregating ESA CCI ocean-colour chlorophyll observations onto the analysis grid.

### GLODAP

The GLODAP workflows prepare discrete ocean biogeochemical observations for comparison with model simulations.

Processing includes:

* extracting selected variables from the GLODAP merged dataset;
* standardizing variable names;
* handling missing values;
* salinity normalization of selected carbonate-system variables;
* filling selected observational gaps using mapped climatological fields where required.

Common GLODAP variables include temperature, salinity, oxygen, dissolved inorganic carbon, alkalinity, nitrate, phosphate, and related carbonate-system properties.

### SODA

The SODA workflow downloads and prepares ocean physical fields, including temperature, salinity, vertical velocity, and horizontal currents. Variables and coordinates are standardized and spatially aggregated to the analysis grid.

Observation-preparation workflows are intentionally dataset-specific because the source products have substantially different structures, resolutions, and sampling characteristics.

See `obs_data_prep/README.md` for details.

---

# Analysis Modules

The `Modules/` package contains the shared analysis and plotting functionality used throughout the notebooks.

These utilities support operations such as:

* model and observational data loading;
* temporal and spatial subsetting;
* regional and biome averaging;
* climatology and anomaly calculation;
* detrending;
* ENSO event selection;
* model–observation matching;
* vertical-level processing;
* derived-variable calculations;
* statistical diagnostics;
* map, profile, scatter, and time-series visualization.

Keeping these operations in reusable modules allows the notebooks to remain primarily focused on **scientific configuration and analysis**, rather than duplicating data-processing code.

---

# Analysis Notebooks

The notebooks provide configurable entry points for the main scientific analyses in the repository.

## Global Maps

`global_map_notebooks/`

These workflows examine the climatological spatial structure of ocean biogeochemical and physical properties.

### Climatological Maps

The climatological-map workflow compares spatial climatologies across model simulations and observational products.

It supports annual, seasonal, and individual-month climatologies. Monthly-resolution data can also be retained for event-based analyses, including **El Niño and La Niña composites and their differences**.

GLODAP climatologies are included in a separate notebook since the GLODAP published a full coverage spatial-vertical climatology in netcdf which are gap filled as well as point measurements.

### ONI Correlation Maps

This workflow calculates grid-point relationships between ocean variables anomalies and the **Oceanic Niño Index (ONI)**, allowing the spatial response of physical and biogeochemical variables to ENSO variability to be investigated.


---

## Regional Time Series

`timeseries_notebooks/`

These notebooks examine temporal variability within selected ocean regions and biomes.

### Regional Seasonal Cycles

Calculates and compares climatological seasonal cycles across models and observational products, allowing differences in the timing, magnitude, and amplitude of seasonal variability to be assessed.

### Regional Time Series

Examines the temporal evolution of regionally averaged anomalies, including long-term variability, trends, interannual variability, and differences between experiments or datasets.

Depending on the analysis, optional preprocessing can include climatology calculation, detrending, surface \(pCO_2\) decomposition, and ENSO event selection.

---

## GLODAP Model–Observation Comparisons

`GLODAP_notebooks/`

These workflows provide detailed comparisons between model simulations and discrete **GLODAP observations**.

Rather than simply selecting the nearest individual observation to a model point, GLODAP measurements are **grouped around model-grid locations and depth using a configurable maximum spatial distance**. Model-grid locations can also be required to contain a minimum number of observations before being retained.

This provides a consistent framework for comparing irregularly sampled GLODAP measurements with gridded model output.

### Scatter Maps

Displays the spatial distribution of model–observation comparisons at the retained model-grid locations.

### Scatter Plots

Evaluates pointwise model–observation relationships using scatter diagrams. The plots can include density information, statistical diagnostics such as \(R^2\) and RMSE, geographic highlighting, and coloring by an additional variable.

### Vertical Profiles

Compares the vertical structure of GLODAP observations and model simulations as a function of depth within selected regions or grid cell with configureable inclusion distance.

### Time Series

Examines the temporal evolution of GLODAP observations and corresponding model simulations where the observational sampling permits temporal comparisons.

Optional GLODAP analyses can include:

* calculation of model climatologies;
* derivation of carbonate-system variables;
* calculation of alkalinity minus dissolved inorganic carbon;
* calculation of differences between assimilation and historical experiments.

---

## Tropical Ocean Analysis

`tropics_notebooks/`

The tropical-ocean notebooks focus on the spatial, vertical, and temporal structure of physical and biogeochemical variability in the equatorial oceans.

These analyses are particularly useful for examining **ENSO-related ocean variability and the coupling between physical circulation and biogeochemical responses**.

### Depth–Latitude Cross-Sections

Examines the **vertical and meridional structure** of anomalies using latitude–depth cross-sections.

The workflow can be used to compare the location, magnitude, and vertical extent of large-scale tropical ocean features across models, observations, and experiments.

Biomes used by this analysis include:

* Pacific Equatorial (PEQ)
* Pacific Equatorial East (PEQ-E)
* Pacific Equatorial West (PEQ-W)
* Atlantic Equatorial (AEQ)

### Depth–Longitude Cross-Sections

Examines the **vertical and zonal structure** of tropical ocean anomalies using longitude–depth sections.

This is particularly useful for investigating zonal gradients and the vertical structure of the equatorial Pacific and other tropical basins.

### Depth–Time Series

Uses depth–time diagrams to investigate the **temporal evolution of subsurface physical and biogeochemical anomalies**.

These diagrams provide a direct view of how vertical ocean structure changes through time and can be related to ENSO variability.

### Longitude–Latitude Aerial Views

Examines the **horizontal structure of tropical ocean anomalies** using longitude–latitude maps.

These maps complement the vertical cross-sections by showing the horizontal extent and spatial organization of physical and biogeochemical responses.

---

# Biomes and Regional Analysis

Many workflows support analysis within predefined ocean **regions and biomes**.

Biome masks provide a consistent method for spatially aggregating model and observational fields and allow comparisons between dynamically and biogeochemically distinct ocean environments.

The exact biome definitions used by an analysis are documented within the corresponding notebook.

---

# ENSO Analysis

ENSO variability is an important component of several workflows in this repository.

Where supported, `mask_ENSO_events` can be used to restrict an analysis to selected **El Niño or La Niña events**.

Depending on the workflow, selected events can:

* retain their original time dimensions for subsequent seasonal analysis;
* be averaged by year and month;
* be combined into climatological ENSO composites;
* be used to calculate **La Niña minus El Niño composite differences**.

The latter represents climatological differences calculated across the selected months and ENSO events.

ENSO timing and related configuration information are maintained centrally in the project configuration files.

---

# Derived Physical and Biogeochemical Quantities

Several workflows can optionally calculate additional quantities when the required source variables are available.

Examples include:

### Ekman Transport

`calculate_ekman` derives wind-driven Ekman transport from the required wind-stress fields. This can be used to investigate relationships between atmospheric forcing, upper-ocean circulation, and biogeochemical variability.

### Seawater Density

`calculate_density` derives seawater density from the required temperature and salinity fields, allowing changes in ocean stratification and water-mass structure to be examined alongside biogeochemical variables.

### Surface \(pCO_2\) Decomposition

Where surface \(pCO_2\) and sea-surface temperature are available, `decompose_spco2` can separate \(pCO_2\) variability into a **temperature-driven component** and a **residual component**.

### Detrending

Selected variables can be detrended prior to analysis. Trends can be estimated over the complete time series, independently for individual months, or using **LOWESS (locally weighted scatterplot smoothing)** when nonlinear trend estimation is desired.

---

# Typical Workflow

Most analyses follow the same general pattern:

```text
1. Select variables
        ↓
2. Select observations / models / experiments
        ↓
3. Define temporal and spatial domain
        ↓
4. Configure optional preprocessing
        ↓
5. Load standardized data
        ↓
6. Calculate derived variables / anomalies / composites
        ↓
7. Apply regional or biome averaging
        ↓
8. Calculate diagnostics
        ↓
9. Visualize and compare results
```

This separation between **configuration, reusable processing routines, and visualization** makes it possible to apply the same scientific analysis consistently across multiple variables, datasets, and model experiments.

---

# Documentation

More detailed documentation is provided within the individual components of the repository:

* **`configs/README.md`** – shared project configurations, units, aliases, regions, and plotting settings.
* **`model_data_prep/README.md`** – model-data preparation and experiment-specific processing.
* **`obs_data_prep/README.md`** – observational-data preparation workflows.
* **Notebook introductions** – scientific purpose, configuration options, biomes, and optional calculations specific to each analysis.

Users should consult the notebook-level documentation for the exact configuration options supported by a particular analysis.

---

# Python Dependencies

BGC_skill is built around the scientific Python ecosystem and requires several packages for multidimensional data processing, statistical analysis, oceanographic calculations, and visualization.

The main dependencies include:

### Core scientific computing

* **NumPy** – numerical array operations
* **pandas** – tabular data manipulation, particularly for discrete observational datasets such as GLODAP
* **xarray** – labelled multidimensional arrays and the primary interface for model and gridded observational data
* **SciPy** – scientific and statistical calculations
* **Dask** – parallel and out-of-core processing of large xarray datasets

### Data I/O and processing

* **netCDF4** – reading and writing NetCDF datasets
* **cftime** – handling climate-model calendars and non-standard datetime coordinates
* **PyYAML** – reading the YAML configuration files used throughout the repository

### Oceanographic and geospatial analysis

* **xESMF** – spatial regridding of model and observational datasets
* **GSW (TEOS-10)** – seawater thermodynamic calculations, including quantities used in density calculations
* **statsmodels** – statistical analysis, including LOWESS-based nonlinear trend estimation

### Visualization

* **Matplotlib** – general plotting and figure generation
* **Cartopy** – geographic map projections and coastlines
* **cmocean** – oceanographically appropriate scientific colour maps

### Notebook environment

* **JupyterLab** or **Jupyter Notebook** – required to run the analysis notebooks interactively
* **IPython** – interactive Python utilities used within the notebook workflows

Some workflows may require additional packages depending on the selected dataset or analysis.

## Installation

A typical environment can be created using Conda or Mamba. For example:

```bash
conda create -n bgc_skill \
    --override-channels \
    -c conda-forge \
    python \
    numpy \
    pandas \
    xarray \
    scipy \
    dask \
    netcdf4 \
    cftime \
    pyyaml \
    xesmf \
    gsw \
    statsmodels \
    matplotlib \
    cartopy \
    cmocean \
    jupyterlab
```

Then activate the environment with:

```bash
conda activate bgc_skill
```

Using **conda** is **not** recommended as it required license. However, packages such as Cartopy and xESMF rely on compiled geospatial and regridding libraries that are generally easier to install through Conda/Mamba than individually through `pip`. If using conda, make sure **conda-forge** channel is being used and **completely avoid** using **default** channels. **Do not use** conda unless you are absolutely certain the **default** channel is **not** being accessed. Otherwise, use at your own responsibility.

> **Note:** Exact package versions are not currently specified here. For reproducible analyses, it is recommended to maintain an `environment.yml` or equivalent dependency file containing the versions used for the project.


---

## Citation and Use

If this repository contributes to published research, please cite the associated publication(s) and datasets as appropriate. Dataset-specific citation requirements should also be followed for observational and model products used in an analysis.

## Disclaimer

This repository contains research software developed for scientific analysis. Individual workflows may rely on project-specific directory structures, datasets, and computing environments and may therefore require configuration before being applied in other environments.

**AI Assistance Notice:** Documentation in this repository was produced with the assistance of artificial intelligence (AI) and reviewed by the project author(s).


## Contributors
This work was developed by **Parsa Gooya** in collaboration with the **Ocean Predictions Group** at the **Canadian Centre for Climate Modeling and Analysis**.

## Copyright
© Environment and Climate Change Canada and the contributors, 2025. All rights reserved.  
For inquiries, contact **parsa.gooya@ec.gc.ca**.  
Do not copy or reproduce without proper citation.
