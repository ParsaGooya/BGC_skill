from __future__ import annotations
import numpy as np
import pandas as pd
import xarray as xr
import seawater as sw

from typing import Sequence

from modules.analysis.module_data_postprocessing import haversine

COMMON_KEYS = ["year", "time", "lat", "lon", "lev"]

def is_depth_range(depth) -> bool:
    return isinstance(depth, (tuple, list))


def rmse(ds1, ds2):
    return np.sqrt(((ds1 - ds2)**2).mean())

def corr(ds1, ds2):
    mask = ~np.isnan(ds1) & ~np.isnan(ds2)
    ds1_valid = ds1[mask]
    ds2_valid = ds2[mask]
    return np.corrcoef(ds1_valid, ds2_valid)[0, 1]


def _filter_depth(df: pd.DataFrame, depth):
    if depth is None:
        return df

    if isinstance(depth, (tuple, list)):
        return df[(df["lev"] >= depth[0]) & (df["lev"] <= depth[1])]

    return df[df["lev"] == depth]


def _reduce_depth_xr(da: xr.DataArray, depth=None) -> xr.DataArray:
    da = _filter_depth_xr(da, depth)

    if "lev" in da.dims:
        if is_depth_range(depth) or depth is None:
            return da.mean("lev")
        return da

    return da



def _filter_depth_xr(da: xr.DataArray, depth=None) -> xr.DataArray:
    if depth is None:
        return da

    if is_depth_range(depth):
        return da.where((da.lev >= depth[0]) & (da.lev <= depth[1]), drop=True)

    return da.sel(lev=depth)

def _filter_bounds_xr(
    da: xr.DataArray,
    boundaries: dict,
) -> xr.DataArray:
    da = da.where(
        (da["lat"] >= boundaries["lat_min"])
        & (da["lat"] <= boundaries["lat_max"])
    )
    da = da.where(
        (da["lon"] >= boundaries["lon_min"])
        & (da["lon"] <= boundaries["lon_max"])
    )
    return da



def _filter_month(df: pd.DataFrame, month: int | None):
    if month is None:
        return df

    return df[df["time"] == month]

def _filter_location(
    df: pd.DataFrame,
    lat: float | None = None,
    lon: float | None = None,
    max_dist_km: float | None = None,
):
    if lat is None or lon is None:
        return df

    df = df.copy()
    df["dist"] = [
        haversine(row_lat, row_lon, lat, lon)
        for row_lat, row_lon in zip(df["lat"], df["lon"])
    ]

    if max_dist_km is None:
        return df[df["dist"] == df["dist"].min()]

    return df[df["dist"] <= max_dist_km * 1000]



def _apply_filters(
    df: pd.DataFrame,
    *,
    depth=None,
    month=None,
    lat=None,
    lon=None,
    max_dist_km=None,
):
    df = _filter_depth(df, depth)
    df = _filter_month(df, month)
    df = _filter_location(df, lat, lon, max_dist_km)
    return df




def _align_dataframes(
    dfs: dict[str, pd.DataFrame],
    keys: Sequence[str] = COMMON_KEYS,
) -> dict[str, pd.DataFrame]:

    keys = list(keys)

    common = dfs[list(dfs.keys())[0]][keys].drop_duplicates()

    for df in list(dfs.values())[1:]:
        common = common.merge(
            df[keys].drop_duplicates(),
            on=keys,
            how="inner",
        )

    return {
        name: df.merge(common, on=keys, how="inner")
        for name, df in dfs.items()
    }

def select_matching(ds_list: list[str], index: int) -> str:
    if len(ds_list) == 1:
        return ds_list[0]
    return ds_list[index]




def _add_density_contours(ax):
    SS, TT = np.meshgrid(
        np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100),
        np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], 100),
    )

    rho = sw.dens(SS, TT, 0)
    sigma0 = rho - 1000

    contour = ax.contour(
        SS,
        TT,
        sigma0,
        levels=15,
        colors="black",
        alpha=0.5,
    )
    ax.clabel(contour, inline=True, fontsize=8)



def get_color(colors_dict, var: str, ds: str, default=None):
    if colors_dict is None:
        return default

    if var in colors_dict and ds in colors_dict[var]:
        return colors_dict[var][ds].get("color", default)

    if ds in colors_dict:
        value = colors_dict[ds]
        if isinstance(value, dict):
            return value.get("color", default)
        return value

    return default



def score_pattern(
    data_frame: pd.DataFrame,
    model_run: str,
    ref: str,
    score_function: str = "MSE",
) -> pd.DataFrame:
    ds = data_frame.copy()

    if score_function.lower() == "mse":
        col = f"{model_run} - {ref}"
        ds[col] = (ds[model_run] - ds[ref]) ** 2

        return (
            ds.groupby(["lat", "lon"])
            .agg(
                mse=(col, "mean"),
                count=(col, "count"),
            )
            .reset_index()
        )

    if score_function.lower() == "r2":
        return (
            ds.groupby(["lat", "lon"])
            .apply(lambda g: g[model_run].corr(g[ref]))
            .reset_index(name="corr")
        )

    raise ValueError(f"Unknown score_function: {score_function}")


def average_on_ref_times(ds: xr.DataArray, ref_dataframe: pd.DataFrame) -> xr.DataArray:
    ref_times = (
        ref_dataframe[["year", "time"]]
        .drop_duplicates()
        .sort_values(["year", "time"])
        .reset_index(drop=True)
    )

    return xr.concat(
        [
            ds.sel(year=row.year, time=row.time)
            for row in ref_times.itertuples(index=False)
        ],
        dim="t",
    )


markers = [
    ".",   # point
    ",",   # pixel
    "o",   # circle
    "v",   # triangle down
    "^",   # triangle up
    "<",   # triangle left
    ">",   # triangle right
    "1",   # tri-down
    "2",   # tri-up
    "3",   # tri-left
    "4",   # tri-right
    "8",   # octagon
    "s",   # square
    "p",   # pentagon
    "P",   # plus (filled)
    "*",   # star
    "h",   # hexagon 1
    "H",   # hexagon 2
    "+",   # plus
    "x",   # x
    "D",   # diamond
    "d",   # thin diamond
    "|",   # vertical line
    "_",   # horizontal line
]
