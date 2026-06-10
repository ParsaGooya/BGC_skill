from pathlib import Path
import yaml

from .config import DataPrepConfig
from .pipeline import run_data_prep

'''
To prepare model data for analysis, edit config.yaml and run: 

python -m data_prep.

Note that you need to have yaml package installed:

pip install pyyaml


author: Parsa Gooya (parsa.gooya@ec.gc.ca)

'''

def load_config(path: str | Path) -> DataPrepConfig:
    with open(path, "r") as f:
        data = yaml.safe_load(f)

    return DataPrepConfig(**data)


if __name__ == "__main__":
    cfg = load_config("config.yaml")
    run_data_prep(cfg)