from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import numpy as np



@dataclass
class DataPrepConfig:
    """
    
    Configuration for CanESM data preparation.
    All user-editable choices are kept here.

    """

    regrid: bool = True
    model: str = "CanESM5"
    var_list: Sequence[str] = ("mlotst",)
    realm_list: Sequence[str] = ("Omon",)
    lev_range: float | None = 600

    realizations: Sequence[str] | None = ("r1i1p2f1",)

    hindcast: bool = False
    hindcast_initial_year: int = 1975
    hindcast_final_year: int = 2025

    assimilation: bool = False
    assimilation_sensitivity: bool = False
    assimilation_extracted_from_disc: bool = False
    extension_years: Sequence[int] = field(default_factory=lambda: tuple(np.arange(2021, 2024)))

    historical: bool = True
    historical_initial_year: int = 1950
    historical_final_year: int = 2025

    canoe_assimilation_bgc: bool = False

    # Base directories. Change these once here if the filesystem changes.
    user_data_root: Path = Path("/space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data")
    cmip6_root: Path = Path("/space/hall5/sitestore/eccc/crd/ccrn/model_output/CMIP6/final/CMIP6")
    canoe_bgc_root: Path = Path("/space/hall6/sitestore/eccc/crd/ccrn/users/scrd107/canesm_runs")
    assim_sensitivity_root = Path("/home/scrd107/site6")
    
    version: str = "v20190429"
    grid: str = "gn"

    def realm_for(self, index: int) -> str:
        """Return matching realm for var_list[index]."""
        if len(self.realm_list) == 1:
            return self.realm_list[0]
        return self.realm_list[index]
