from __future__ import annotations

import glob
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import xarray as xr
import xesmf as xe

from .config import DataPrepConfig


def coords_edit(da: xr.DataArray) -> xr.DataArray:
    """Convert x/y grid coordinates from xESMF output to lat/lon dimensions."""
    lat = da.lat.values[:, 0]
    lon = da.lon.values[0, :]

    out = (
        xr.DataArray(da.values, dims=da.dims)
        .rename({"x": "lon", "y": "lat"})
        .assign_coords({"lat": lat, "lon": lon})
    )
    for dim in da.dims:
        if dim not in {"x", "y"}:
            out = out.assign_coords({dim: da[dim]})
    return out


def available_realizations(base_dir: Path, *, exclude_p1: bool = False) -> list[str]:
    """Infer available realization IDs from subdirectory names."""
    realizations: list[str] = []
    for path in glob.glob(str(base_dir / "*")):
        name = Path(path).name
        realization = name.split("-")[-1]
        if exclude_p1 and "p1" in realization:
            continue
        realizations.append(realization)
    return sorted(np.unique(realizations).tolist())


def resolve_realizations(cfg: DataPrepConfig, base_dir: Path, *, exclude_p1: bool = False) -> list[str]:
    if cfg.realizations is not None:
        return list(cfg.realizations)
    return available_realizations(base_dir, exclude_p1=exclude_p1)


def limit_level(ds: xr.Dataset, lev_range: float | None) -> xr.Dataset:
    if lev_range is not None and "lev" in ds.dims:
        return ds.where(ds.lev <= lev_range, drop=True)
    return ds


def concat_members(dsets: dict[str, xr.Dataset], member_dim: str = "member") -> xr.Dataset:
    if not dsets:
        raise FileNotFoundError("No datasets were found/opened for the requested realizations.")
    return xr.concat(list(dsets.values()), dim=member_dim).assign_coords({member_dim: list(dsets.keys())})


def open_member_dataset(path: Path, cfg: DataPrepConfig) -> xr.Dataset:
    files = sorted(path.glob("*.nc"))
    if not files:
        raise FileNotFoundError(f"No NetCDF files found under {path}")
    ds = xr.open_mfdataset(str(path / "*.nc"), combine="nested", concat_dim="time")
    return limit_level(ds, cfg.lev_range)


def open_initialized_dataset(path: Path, cfg: DataPrepConfig) -> xr.Dataset:
    files = sorted(path.glob("*.nc"))
    if not files:
        raise FileNotFoundError(f"No NetCDF files found under {path}")
    ds = xr.open_dataset(files[0])
    return limit_level(ds, cfg.lev_range)


def regrid(ds: xr.Dataset, var: str, cfg: DataPrepConfig) -> xr.DataArray:
    da = ds[var]
    if not cfg.regrid:
        return da

    ds_out = xe.util.grid_global(1, 1)
    regridder = xe.Regridder(ds, ds_out, "bilinear", ignore_degenerate=True, periodic=True)
    return coords_edit(regridder(da))


def standardize_member_time(da: xr.DataArray, *, time_dim: str = "time") -> xr.DataArray:
    return da.rename({"member": "ensembles"}).transpose("ensembles", time_dim, ..., "lat", "lon")


def save_dataarray(da: xr.DataArray, var: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    da.to_dataset(name=var).to_netcdf(out_path)


def year_from_time_value(value) -> int:
    return value.item().year
