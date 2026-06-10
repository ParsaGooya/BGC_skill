from __future__ import annotations

import glob
from pathlib import Path

import xarray as xr

from .config import DataPrepConfig
from .paths import dcpp_dir, output_dir, variable_dir
from .utils import concat_members, limit_level, regrid, open_member_dataset, resolve_realizations, save_dataarray, standardize_member_time, year_from_time_value





def sensitivity_out_dir(cfg: DataPrepConfig, var: str, sensitivity_index: int) -> Path:
    return cfg.user_data_root / var / f"assimilation_bgc{sensitivity_index}" / cfg.model


def sensitivity_input_dir(cfg: DataPrepConfig, realm: str, var: str, sensitivity_index: int) -> Path:
    run_name = f"d2a-asm-bgc{sensitivity_index}"
    model_name = f"CanESM5-d2a-asm-bgc{sensitivity_index}"
    return (
        cfg.assim_sensitivity_root
        / run_name
        / "data/nc_output/CMIP6/CCCma/CCCma"
        / model_name
        / "dcppA-assim/r1i1p2f1"
        / realm
        / var
        / cfg.grid
        / cfg.version
    )




def process_assimilation(cfg: DataPrepConfig, var: str, realm: str) -> None:
    '''prepare from dcpp-assim runs from saved CMIP style experiments on the disk'''
    base_dir = dcpp_dir(cfg, "dcppA-assim")
    realizations = resolve_realizations(cfg, base_dir)
    print(f"Available realizations - Assimilation: {realizations}")

    dsets = {}
    for r in realizations:
        try:
            dsets[r] = open_member_dataset(variable_dir(base_dir, r, realm, var, cfg), cfg)
        except Exception as exc:
            print(f"Skipping assimilation {r}: {exc}")

    ds = concat_members(dsets)
    da = standardize_member_time(regrid(ds, var, cfg))

    out_dir = output_dir(cfg, var, "assimilation")
    start = year_from_time_value(da.time[0].values)
    end = year_from_time_value(da.time[-1].values)
    save_dataarray(da, var, out_dir / f"{var}_{realm}_ensmebles_{start}01_{end}12_1x1_LE.nc")


def process_assimilation_extensions(cfg: DataPrepConfig, var: str, realm: str) -> None:
    '''
    
    prepare from dcpp-assim runs not included in the saved CMIP style experiments on the disk.
    Need ro run extract_from_disk.sh first to first create the netcdf files in the same output 
    directory as where you stored CIP6 data for the variable. See  process_assimilation function.
    
    ''' 
    
    out_dir = output_dir(cfg, var, "assimilation")

    yearly = []
    for year in cfg.extension_years:
        files = sorted(glob.glob(str(out_dir / "extentions" / f"*{year}*.nc")))
        if not files:
            print(f"No extension files found for {year}")
            continue
        ensembles = [Path(link).name.split("_")[4] for link in files]
        yearly.append(
            xr.concat([xr.open_dataset(link) for link in files], dim="ensembles")
            .assign_coords(ensembles=ensembles))
        

    if not yearly:
        raise FileNotFoundError("No extension files were found.")

    ds = limit_level(xr.concat(yearly, dim="time"), cfg.lev_range)
    da = regrid(ds, var, cfg)
    if "ensembles" in ds.coords:
        da = da.assign_coords(ensembles=ds.ensembles)

    start = year_from_time_value(da.time[0].values)
    end = year_from_time_value(da.time[-1].values)
    save_dataarray(da, var, out_dir / f"{var}_{realm}_ensmebles_{start}01_{end}12_1x1_LE.nc")





def process_assimilation_sensitivity(cfg: DataPrepConfig) -> None:
    """Process assimilation sensitivity BGC runs.

    This reproduces the original loop over bgc0, bgc1, bgc2, bgc3, but uses
    the shared config, level limiting, regridding, and saving helpers.
    """
    for sensitivity_index in range(4):
        print("=========================================================")
        print(f"Loading assimilation sensitivity run: bgc{sensitivity_index}")

        for index, var in enumerate(cfg.assim_sensitivity_var_list):
            realm = cfg.assim_sensitivity_realm_for(index)
            in_dir = sensitivity_input_dir(cfg, realm, var, sensitivity_index)
            out_dir = sensitivity_out_dir(cfg, var, sensitivity_index)

            print(f"Processing {var=} {realm=} from {in_dir}")

            ds = xr.open_mfdataset(
                str(in_dir / "*.nc"),
                combine="nested",
                concat_dim="time",
            )
            ds = limit_level(ds, cfg.lev_range)

            da = regrid(ds, var, cfg)
            da = da.expand_dims("ensembles", axis=0).transpose(
                "ensembles", "time", ..., "lat", "lon"
            )

            start_year = year_from_time_value(da.time[0].values)
            end_year = year_from_time_value(da.time[-1].values)
            out_path = out_dir / f"{var}_{realm}_ensmebles_{start_year}01_{end_year}12_1x1_LE.nc"

            save_dataarray(da, var, out_path)
            print(f"Assimilation sensitivity {var} data saved to {out_path}\n")
