from __future__ import annotations

from typing import Sequence
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import numpy as np


from modules.plotting.utils import (
    _align_dataframes,
    COMMON_KEYS,
    get_color, 
    _filter_depth, 
    rmse, corr)



def _prepare_aligned_variable_dataframe(
    dataframe: dict,
    *,
    var: str,
    location_ref_var: str | None = None,
    outlier_var=None,
    keys: Sequence[str] = COMMON_KEYS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if location_ref_var is None:
        location_ref_var = var

    location_ref = dataframe[location_ref_var]
    data = dataframe[var]

    if outlier_var is not None:
        location_ref = location_ref[location_ref["var"] < outlier_var]

    aligned = _align_dataframes(
        {
            "data": data,
            "location_ref": location_ref,
        },
        keys=keys,
    )

    return aligned["data"], aligned["location_ref"]




def timeseries_comaprisons(
    dataframe,
    var: str,
    ds_list: list,
    units: dict,
    colors_dict: dict,
    location_ref_var=None,
    ref_ds="obs",
    depth_range: list | tuple = None,
    outlier_var=None,
    figsize=(50, 6),
    fontsize=20,
    calculate_rmse=False,
    calculate_r2=False,
    add_grid=False,
):
    ds_toplot, location_ref = _prepare_aligned_variable_dataframe(
        dataframe,
        var=var,
        location_ref_var=location_ref_var,
        outlier_var=outlier_var,
    )

    if location_ref_var is None:
        location_ref_var = var

    ds_toplot["yearmonth"] = (
        ds_toplot["year"] + (ds_toplot["time"] + 0.5) / 12
    )

    title = (
        f"{var} ({units[var]}) - "
        f"{location_ref['year'].min()} - {location_ref['year'].max()}"
    )

    if location_ref_var != var:
        title += f"\n location reference = {location_ref_var}"

    if depth_range is not None:
        ds_toplot = _filter_depth(ds_toplot, depth_range)
        title += f" {depth_range[0]} - {depth_range[1]} m"

    ds_toplot_spatial_mean = (
        ds_toplot.groupby(["yearmonth"])
        .mean(numeric_only=True)
        .reset_index()
    )

    offsets = np.arange(
        -((len(ds_list) + 1) // 2),
        len(ds_list) + 1 - ((len(ds_list) + 1) // 2),
    )

    datetimes = [
        datetime.datetime(int(row.year), int(row.time) + 1, 15)
        for row in ds_toplot.itertuples(index=False)
    ]

    xtick_dates = [dt.strftime("%Y-%m") for dt in datetimes]

    plt.figure(figsize=figsize)

    for ind, ds in enumerate(ds_list):
        color = get_color(colors_dict, var, ds)

        plt.scatter(
            ds_toplot["yearmonth"] + offsets[ind] * 0.05,
            ds_toplot[ds],
            color=color,
            label=ds,
            alpha=0.5,
        )

        plt.scatter(
            ds_toplot_spatial_mean["yearmonth"] + offsets[ind] * 0.05,
            ds_toplot_spatial_mean[ds],
            marker="_",
            s=100,
            color=color,
        )

    if calculate_rmse:
        rmse_scores = {
            ds: np.round(rmse(ds_toplot[ref_ds], ds_toplot[ds]), 2)
            for ds in ds_list
            if ds != ref_ds
        }
        title += "\nRMSE: " + ", ".join(
            f"{ds}={score}" for ds, score in rmse_scores.items()
        )

    if calculate_r2:
        r2_scores = {
            ds: np.round(corr(ds_toplot[ref_ds], ds_toplot[ds]) ** 2, 2)
            for ds in ds_list
            if ds != ref_ds
        }
        title += "\n$R^2$: " + ", ".join(
            f"{ds}={score}" for ds, score in r2_scores.items()
        )

    plt.scatter(
        ds_toplot["yearmonth"] + offsets[-1] * 0.05,
        ds_toplot[ref_ds],
        marker="x",
        color="k",
        label=ref_ds,
        alpha=0.5,
    )

    plt.scatter(
        ds_toplot_spatial_mean["yearmonth"] + offsets[-1] * 0.05,
        ds_toplot_spatial_mean[ref_ds],
        marker="_",
        color="k",
    )

    xticks = (
        ds_toplot[["yearmonth", "year", "time"]]
        .drop_duplicates()
        .sort_values("yearmonth")
    )

    xtick_labels = [
        f"{int(row.year)}-{int(row.time) + 1:02d}"
        for row in xticks.itertuples(index=False)
    ]

    plt.xticks(
        xticks["yearmonth"],
        xtick_labels,
        rotation=90,
        fontsize=fontsize,
    )

    plt.yticks(fontsize=fontsize)
    plt.ylabel(f"{units[var]}", fontsize=fontsize)
    plt.xlabel("time", fontsize=fontsize)
    plt.title(title, fontsize=fontsize)
    plt.legend(fontsize=fontsize)

    if add_grid:
        plt.grid()