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

with open(CONFIG_DIR / "CanESM_ocean_level_boundaries.json", "r") as f:
    model_lev_bounds = np.array(json.load(f))


with open(CONFIG_DIR / "CanESM_ocean_levels.json", "r") as f:
    model_levels = np.array(json.load(f))