
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Sequence

from modules.plotting.utils import (
    COMMON_KEYS,
    _align_dataframes, 
    _apply_filters, 
    select_matching, 
    _add_density_contours, 
    markers, 
    rmse, 
    corr
)


def _add_reference_lines(ax, add_1_1_line=True, add_zeros_grid=False):
    if add_1_1_line:
        lo = min(ax.get_xlim()[0], ax.get_ylim()[0])
        hi = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([lo, hi], [lo, hi], color="k", alpha=1)

    if add_zeros_grid:
        ax.axhline(0, color="k", alpha=0.5)
        ax.axvline(0, color="k", alpha=0.5)




def prepare_scatter_data(
    dataframe: dict[str, pd.DataFrame],
    *,
    var_y: str,
    var_x: str | None = None,
    location_ref_var: str | None = None,
    color_bar_var: str | None = None,
    outlier_variance: float | None = None,
    depth=None,
    month: int | None = None,
    lat: float | None = None,
    lon: float | None = None,
    max_dist_km: float | None = None,
    title_prefix: str | None = None,
    keys: Sequence[str] = COMMON_KEYS,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None, str]:
    
    var_x = var_y if var_x is None else var_x
    location_ref_var = var_y if location_ref_var is None else location_ref_var

    dfs = {
        "y": dataframe[var_y],
        "x": dataframe[var_x],
        "location_ref": dataframe[location_ref_var],
    }

    if color_bar_var is not None and color_bar_var.lower() != "depth":
        dfs["color"] = dataframe[color_bar_var]

    if outlier_variance is not None:
        dfs["location_ref"] = dfs["location_ref"][
            dfs["location_ref"]["var"] < outlier_variance
        ]

    dfs = _align_dataframes(dfs, keys=keys)

    for name in dfs:
        dfs[name] = _apply_filters(
            dfs[name],
            depth=depth,
            month=month,
            lat=lat,
            lon=lon,
            max_dist_km=max_dist_km,
        )

    y_df = dfs["y"]
    x_df = dfs["x"]
    color_df = dfs.get("color")

    title = f"{y_df['year'].min()} - {y_df['year'].max()}"

    if title_prefix is not None:
        title = f"{title_prefix} - {title}"

    if var_x != var_y or location_ref_var != var_y:
        title += f"\nlocation reference = {location_ref_var}"

    if depth is not None:
        title += f" | lev = {depth}"

    if month is not None:
        title += f" | month = {month}"

    if lat is not None and lon is not None:
        if max_dist_km is None:
            title += f" | nearest: lat={x_df['lat'].iloc[0]}, lon={x_df['lon'].iloc[0]}"
        else:
            title += f"\nwithin {max_dist_km} km from lat={lat}, lon={lon}"

    return x_df, y_df, color_df, title



def scatter_comparison(
    dataframe: dict[str, pd.DataFrame],
    var_y: str,
    ds_list_var_y: list[str],
    units: dict,
    *,
    var_x: str | None = None,
    ds_list_var_x: list[str] = ["obs"],
    location_ref_var: str | None = None,
    outlier_variance: float | None = None,
    color_bar_var: str | None = None,
    ds_list_color_bar: list[str] | None = None,
    title_prefix: str | None = None,
    cmap: str = "coolwarm",
    density_contour: bool = False,
    depth=None,
    month: int | None = None,
    lat: float | None = None,
    lon: float | None = None,
    max_dist_km: float | None = None,
    figsize=(16, 6),
    plot_lim_x=None,
    plot_lim_y=None,
    vmin=None,
    vmax=None,
    month_based_colorbar: bool = False,
    calculate_rmse: bool = True,
    calculate_r2: bool = True,
    add_1_1_line: bool = True,
    add_zeros_grid: bool = False,
):
    """
    Scatter comparison for already-selected dataframes.

    Expected input
    --------------
    dataframe is a dictionary of variable -> dataframe:

        dataframe = {
            "talk": talk_df_for_selected_biome,
            "dissic": dissic_df_for_selected_biome,
            "thetao": thetao_df_for_selected_biome,
        }

    Each dataframe should have common location/time columns:

        year, time, lat, lon, lev

    and data columns such as:

        obs, CanESM5_assimilation, CanESM5_historical, ...

    Notes
    -----
    - Biome selection is assumed to happen before calling this function.
    - If var_x is None, var_x = var_y.
    - If location_ref_var is None, location_ref_var = var_y.
    """
    var_x = var_y if var_x is None else var_x

    if density_contour:
        assert (var_y == "thetao" and
            var_x == "so"), 'density contours are only available for TS diagrams.'

    if len(ds_list_var_x) > 1:
        assert len(ds_list_var_x) == len(ds_list_var_y), 'Either compare all  experiments for y variable against a single experiment'\
                                                        'for var x or provide a corresponding experiment to compare against for each.'


    if color_bar_var is not None and color_bar_var.lower() != "depth":
        assert ds_list_color_bar is not None, 'specify at least one experiment for color bar variable.'
        assert not month_based_colorbar, 'you cannot ask for month based colorbar if a color_bar_var is requested.'

        if len(ds_list_color_bar) > 1:
            assert len(ds_list_color_bar) == len(ds_list_var_y), 'Either compare all experiments for x-y variables against a single experiment'\
                                                        'for colorbar var or provide a corresponding experiment to compare against for each.'
    

    x_df, y_df, color_df, title_base = prepare_scatter_data(
        dataframe,
        var_y=var_y,
        var_x=var_x,
        location_ref_var=location_ref_var,
        color_bar_var=color_bar_var,
        outlier_variance=outlier_variance,
        depth=depth,
        month=month,
        lat=lat,
        lon=lon,
        max_dist_km=max_dist_km,
        title_prefix=title_prefix,
    )

    for i, y_ds in enumerate(ds_list_var_y):
        x_ds = select_matching(ds_list_var_x, i)

        if month_based_colorbar:
            fig, axes = plt.subplots(1, 2, figsize=figsize)
            ax0, ax1 = axes
        else:
            fig, ax0 = plt.subplots(1, 1, figsize=figsize)

        if color_bar_var is None or color_bar_var.lower() == "depth":
            c = x_df["lev"]
            color_label = "depth"
        else:
            color_ds = select_matching(ds_list_color_bar, i)
            c = color_df[color_ds]
            color_label = f"{color_bar_var} ({units[color_bar_var]})\n{color_ds}"

        scatter = ax0.scatter(
            x_df[x_ds],
            y_df[y_ds],
            c=c,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )

        fig.colorbar(scatter, ax=ax0, label=color_label)

        title = f"{y_ds} vs {x_ds} - {title_base}"

        if calculate_r2:
            title += f" | $R^2$ = {np.round(corr(x_df[x_ds], y_df[y_ds]) ** 2, 2)}"

        if calculate_rmse:
            title += f" | RMSE = {np.round(rmse(x_df[x_ds], y_df[y_ds]), 2)}"

        ax0.set_title(title)
        ax0.set_xlabel(f"{var_x} ({units[var_x]})")
        ax0.set_ylabel(f"{var_y} ({units[var_y]})")

        if plot_lim_x is not None:
            ax0.set_xlim(*plot_lim_x)
        if plot_lim_y is not None:
            ax0.set_ylim(*plot_lim_y)

        _add_reference_lines(
            ax0,
            add_1_1_line=add_1_1_line,
            add_zeros_grid=add_zeros_grid,
        )

        if density_contour:
            _add_density_contours(ax0)

        if month_based_colorbar:
            scatter_month = ax1.scatter(
                x_df[x_ds],
                y_df[y_ds],
                c=x_df["time"] + 1,
                cmap="tab10",
            )

            fig.colorbar(scatter_month, ax=ax1, label="month")

            ax1.set_xlabel(f"{var_x} ({units[var_x]})")
            ax1.set_ylabel(f"{var_y} ({units[var_y]})")

            if plot_lim_x is not None:
                ax1.set_xlim(*plot_lim_x)
            if plot_lim_y is not None:
                ax1.set_ylim(*plot_lim_y)

            _add_reference_lines(
                ax1,
                add_1_1_line=add_1_1_line,
                add_zeros_grid=add_zeros_grid,
            )

            if density_contour:
                _add_density_contours(ax1)

        plt.tight_layout()






def scatter_comparison_singlepanel(
    dataframe: dict[str, pd.DataFrame],
    var_y: str,
    ds_list: list[str],
    units: dict,
    *,
    var_x: str | None = None,
    location_ref_var: str | None = None,
    outlier_variance: float | None = None,
    shapes_dict: dict | None = None,
    colors_dict: dict | None = None,
    alphas_dict: dict | None = None,
    color_bar_var: str | None = None,
    title_prefix: str | None = None,
    cmap: str = "coolwarm",
    density_contour: bool = False,
    depth=None,
    month: int | None = None,
    lat: float | None = None,
    lon: float | None = None,
    max_dist_km: float | None = None,
    figsize=(16, 6),
    plot_lim_x=None,
    plot_lim_y=None,
    vmin=None,
    vmax=None,
    month_based_colorbar: bool = False,
    calculate_rmse: bool = True,
    calculate_r2: bool = True,
    add_1_1_line: bool = True,
    add_zeros_grid: bool = False,
):
    
    if shapes_dict is None:
        shapes_dict = {}
        for ind, ds in enumerate(ds_list):
            shapes_dict[ds] = markers[ind]

    var_x = var_y if var_x is None else var_x

    if density_contour:
        assert (var_y == "thetao" and
            var_x == "so"), 'density contours are only available for TS diagrams.'


    if color_bar_var is None:
        assert colors_dict is not None
        assert not month_based_colorbar

    x_df, y_df, color_df, title_base = prepare_scatter_data(
        dataframe,
        var_y=var_y,
        var_x=var_x,
        location_ref_var=location_ref_var,
        color_bar_var=color_bar_var,
        outlier_variance=outlier_variance,
        depth=depth,
        month=month,
        lat=lat,
        lon=lon,
        max_dist_km=max_dist_km,
        title_prefix=title_prefix,
    )

    if month_based_colorbar:
        fig, axes = plt.subplots(1, 2, figsize=figsize)
        ax0, ax1 = axes
    else:
        fig, ax0 = plt.subplots(1, 1, figsize=figsize)

    last_scatter = None

    for i, y_ds in enumerate(ds_list):


        marker = shapes_dict.get(y_ds, "o") if shapes_dict else "o"
        alpha = alphas_dict.get(y_ds, 1) if alphas_dict else 1

        label = y_ds

        if calculate_r2:
            label += f" | $R^2$ = {np.round(corr(x_df[y_ds], y_df[y_ds]) ** 2, 2)}"

        if calculate_rmse:
            label += f" | RMSE = {np.round(rmse(x_df[y_ds], y_df[y_ds]), 2)}"

        if color_bar_var is None:
            last_scatter = ax0.scatter(
                x_df[y_ds],
                y_df[y_ds],
                color=colors_dict[y_ds],
                marker=marker,
                alpha=alpha,
                label=label,
            )

        else:
            if color_bar_var.lower() == "depth":
                c = x_df["lev"]
                color_label = "depth"
            else:
                c = color_df[y_ds]
                color_label = f"{color_bar_var} ({units[color_bar_var]})"

            edgecolors = colors_dict.get(y_ds, None) if colors_dict else None

            last_scatter = ax0.scatter(
                x_df[y_ds],
                y_df[y_ds],
                c=c,
                cmap=cmap,
                edgecolors=edgecolors,
                marker=marker,
                alpha=alpha,
                label=label,
                vmin=vmin,
                vmax=vmax,
            )

        if month_based_colorbar:
            scatter_month = ax1.scatter(
                x_df[y_ds],
                y_df[y_ds],
                c=x_df["time"] + 1,
                cmap="tab10",
                marker=marker,
                alpha=alpha,
            )

    ax0.legend()
    ax0.set_title(title_base)
    ax0.set_xlabel(f"{var_x} ({units[var_x]})")
    ax0.set_ylabel(f"{var_y} ({units[var_y]})")

    if color_bar_var is not None:
        fig.colorbar(last_scatter, ax=ax0, label=color_label)

    if plot_lim_x is not None:
        ax0.set_xlim(*plot_lim_x)
    if plot_lim_y is not None:
        ax0.set_ylim(*plot_lim_y)

    _add_reference_lines(
        ax0,
        add_1_1_line=add_1_1_line,
        add_zeros_grid=add_zeros_grid,
    )

    if density_contour:
        _add_density_contours(ax0)

    if month_based_colorbar:
        fig.colorbar(scatter_month, ax=ax1, label="month")

        ax1.set_xlabel(f"{var_x} ({units[var_x]})")
        ax1.set_ylabel(f"{var_y} ({units[var_y]})")

        if plot_lim_x is not None:
            ax1.set_xlim(*plot_lim_x)
        if plot_lim_y is not None:
            ax1.set_ylim(*plot_lim_y)

        _add_reference_lines(
            ax1,
            add_1_1_line=add_1_1_line,
            add_zeros_grid=add_zeros_grid,
        )

        if density_contour:
            _add_density_contours(ax1)

    plt.tight_layout()




