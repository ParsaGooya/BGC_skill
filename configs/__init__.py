
"""
Shared configuration for the BGC skill analysis package.

This package loads commonly used project configuration files and exposes
them as Python objects. Configuration data include variable aliases,
display units, unit-conversion factors, geographical regions, ENSO event
definitions, plotting ranges, and CanESM ocean vertical coordinates.

Attributes
----------
varx : dict
    Mapping from alternative variable names to canonical variable names.
units : dict
    Display units associated with model and observational variables.
unit_change_dics : dict
    Multiplicative factors used to convert variables to the units adopted
    by the analysis.
boundaries_dict : dict
    Latitude and longitude boundaries for predefined geographical regions.
ENSO : dict
    ENSO event definitions used by ENSO-related analyses.
contourf : dict
    Predefined contour levels for selected variables.
pcolor_ranges : dict
    Default plotting ranges for selected variables and diagnostics.
model_lev_bounds : numpy.ndarray
    Vertical boundaries of the CanESM ocean model levels.
model_levels : numpy.ndarray
    Representative depths of the CanESM ocean model levels.

Notes
-----
The configuration files are resolved relative to this package rather than
the current working directory, allowing ``configs`` to be imported from
scripts and notebooks located elsewhere in the repository.
"""



from pathlib import Path
import yaml
import json
import numpy as np
from . import mplRC

CONFIG_DIR = Path(__file__).parent

with open(CONFIG_DIR / "varx.yaml", "r") as f:
    varx = yaml.safe_load(f)

with open(CONFIG_DIR / "units.yaml", "r") as f:
    units = yaml.safe_load(f)

with open(CONFIG_DIR / "unit_changes.yaml", "r") as f:
    unit_change_dics = yaml.safe_load(f)

with open(CONFIG_DIR / "regions.yaml", "r") as f:
    boundaries_dict = yaml.safe_load(f)

with open(CONFIG_DIR / "ENSO.yaml", "r") as f:
    ENSO = yaml.safe_load(f)

with open(CONFIG_DIR / "contourf_levels.yaml", "r") as f:
    contourf = yaml.safe_load(f)


with open(CONFIG_DIR / "pcolor_ranges.yaml", "r") as f:
    pcolor_ranges = yaml.safe_load(f)

with open(CONFIG_DIR / "CanESM_ocean_level_boundaries.json", "r") as f:
    model_lev_bounds = np.array(json.load(f))


with open(CONFIG_DIR / "CanESM_ocean_levels.json", "r") as f:
    model_levels = np.array(json.load(f))