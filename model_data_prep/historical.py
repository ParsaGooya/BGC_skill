from __future__ import annotations

from .config import DataPrepConfig
from .paths import historical_dir, output_dir, ssp245_dir, variable_dir
from .utils import concat_members, regrid, open_member_dataset, resolve_realizations, save_dataarray, standardize_member_time, year_from_time_value


def _common_historical_realizations(cfg: DataPrepConfig) -> list[str]:
    hist_base = historical_dir(cfg)
    realizations = resolve_realizations(cfg, hist_base, exclude_p1=True)

    if cfg.realizations is None and cfg.historical_final_year > 2015:
        ssp_base = ssp245_dir(cfg)
        ssp_realizations = set(resolve_realizations(cfg, ssp_base, exclude_p1=True))
        realizations = [r for r in realizations if r in ssp_realizations]

    return realizations


def _open_members(cfg: DataPrepConfig, base_dir, var: str, realm: str, realizations: list[str]):
    dsets = {}
    for r in realizations:
        try:
            dsets[r] = open_member_dataset(variable_dir(base_dir, r, realm, var, cfg), cfg)
        except Exception as exc:
            print(f"Skipping {base_dir.name} {r}: {exc}")
    return concat_members(dsets)


def _save_time_series(cfg: DataPrepConfig, da, var: str, realm: str) -> None:
    out_dir = output_dir(cfg, var, "historical")
    start = year_from_time_value(da.time[0].values)
    end = year_from_time_value(da.time[-1].values)
    save_dataarray(da, var, out_dir / f"{var}_{realm}_ensmebles_{start}01_{end}12_1x1_LE.nc")


def process_historical(cfg: DataPrepConfig, var: str, realm: str) -> None:
    realizations = _common_historical_realizations(cfg)
    print(f"Available realizations - Historical: {realizations}")

    if cfg.historical_initial_year < 2015:
        end = None if cfg.historical_final_year > 2015 else str(cfg.historical_final_year)
        ds = _open_members(cfg, historical_dir(cfg), var, realm, realizations)
        ds = ds.sel(time=slice(str(cfg.historical_initial_year), end))
        da = standardize_member_time(regrid(ds, var, cfg))
        _save_time_series(cfg, da, var, realm)

    if cfg.historical_final_year > 2015:
        ds = _open_members(cfg, ssp245_dir(cfg), var, realm, realizations)
        ds = ds.sel(time=slice(None, str(cfg.historical_final_year)))
        da = standardize_member_time(regrid(ds, var, cfg))
        _save_time_series(cfg, da, var, realm)
