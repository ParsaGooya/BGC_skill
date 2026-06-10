from __future__ import annotations

from tqdm import tqdm
import numpy as np
import xarray as xr

from .config import DataPrepConfig
from .paths import dcpp_dir, initialized_variable_dir, output_dir
from .utils import concat_members, regrid, open_initialized_dataset, resolve_realizations, save_dataarray


def _process_initialized_years(
    cfg: DataPrepConfig,
    *,
    var: str,
    realm: str,
    experiment: str,
    first_year: int,
    last_year_exclusive: int,
) -> None:
    base_dir = dcpp_dir(cfg, experiment)
    realizations = resolve_realizations(cfg, base_dir)
    print(f"Available realizations - {experiment}: {realizations}")

    out_dir = output_dir(cfg, var, "forecast")

    for year in tqdm(range(first_year, last_year_exclusive)):
        dsets = {}
        for r in realizations:
            try:
                path = initialized_variable_dir(base_dir, year, r, realm, var, cfg)
                dsets[r] = open_initialized_dataset(path, cfg)
            except Exception as exc:
                print(f"Skipping {experiment} {year} {r}: {exc}")

        ds = concat_members(dsets)
        da = regrid(ds, var, cfg)
        da = da.rename({"time": "lead_time", "member": "ensembles"})
        da = da.assign_coords(lead_time=np.arange(1, da.sizes["lead_time"] + 1), year=year + 1)
        da = da.transpose("ensembles", "lead_time", ..., "lat", "lon")

        out_path = out_dir / f"{var}_{realm}_{year}_ensmebles_{year + 1}01_{year + 10}12_1x1_LE.nc"
        save_dataarray(da, var, out_path)


def process_hindcast_and_forecast(cfg: DataPrepConfig, var: str, realm: str) -> None:
    """Process dcppA-hindcast and, when requested, dcppB-forecast."""
    print(f"Loading Hindcast data: {cfg.hindcast_initial_year} - 2019")
    _process_initialized_years(
        cfg,
        var=var,
        realm=realm,
        experiment="dcppA-hindcast",
        first_year=cfg.hindcast_initial_year,
        last_year_exclusive=2020,
    )

    if cfg.hindcast_final_year > 2019:
        print(f"Loading Forecast data: 2020 - {cfg.hindcast_final_year}")
        _process_initialized_years(
            cfg,
            var=var,
            realm=realm,
            experiment="dcppB-forecast",
            first_year=2020,
            last_year_exclusive=cfg.hindcast_final_year,
        )
