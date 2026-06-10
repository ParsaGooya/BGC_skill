from __future__ import annotations

import glob
import numpy as np
import xarray as xr

from .config import DataPrepConfig
from .paths import output_dir
from .utils import limit_level, regrid, save_dataarray, standardize_member_time, year_from_time_value


def process_canoe_bgc_assimilation(cfg: DataPrepConfig, var: str, realm: str) -> None:
    """Process the newly created CanOE BGC assimilation runs."""
    model = "CanESM5-CanOE_1"
    run_dirs = sorted(glob.glob(str(cfg.canoe_bgc_root / "d2k-asm-*a")))

    datasets = []
    for run_dir in run_dirs:
        member = int(run_dir.split("-")[-1][1:-1])
        pattern = f"{run_dir}/data/nc_output/CMIP6/CCCma/CCCma/*/dcppA-assim/*/{realm}/{var}/{cfg.grid}/{cfg.version}/*.nc"
        ds = xr.open_mfdataset(pattern, combine="nested", concat_dim="time")
        datasets.append(ds.assign_coords(member=np.array(member)))

    if not datasets:
        raise FileNotFoundError("No CanOE BGC assimilation runs were found.")

    ds = limit_level(xr.concat(datasets, dim="member"), cfg.lev_range)
    da = standardize_member_time(regrid(ds, var, cfg))

    out_dir = output_dir(cfg, var, "assimilation", model=model)
    start = year_from_time_value(da.time[0].values)
    end = year_from_time_value(da.time[-1].values)
    save_dataarray(da, var, out_dir / f"{var}_{realm}_ensmebles_{start}01_{end}12_1x1_LE.nc")
