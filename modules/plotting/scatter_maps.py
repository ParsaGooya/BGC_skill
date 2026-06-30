from __future__ import annotations

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


from modules.plotting.utils import (
    _align_dataframes,
    COMMON_KEYS,
    is_depth_range, 
    _filter_depth, 
    _filter_bounds_xr, 
    _reduce_depth_xr, 
    score_pattern,
    rmse, corr)



def _plot_biome_background(
    ax,
    mask_ocean,
    mask_biome=None,
    *,
    lat_min_to_show=None,
    lat_max_to_show=None,
    lon_min_to_show=None,
    lon_max_to_show=None,
    lat=None,
    lon=None
):
    lat_slice = slice(lat_min_to_show, lat_max_to_show)
    lon_slice = slice(lon_min_to_show, lon_max_to_show)

    ax.pcolormesh(
        mask_ocean.lon[lon_slice],
        mask_ocean.lat[lat_slice],
        mask_ocean.where(mask_ocean == 0).values[lat_slice, lon_slice],
        vmin=0,
        vmax=1,
    )

    if mask_biome is not None:
        ax.pcolormesh(
            mask_biome.lon[lon_slice],
            mask_biome.lat[lat_slice],
            mask_biome.where(mask_biome == 1).values[lat_slice, lon_slice],
            vmin=0,
            vmax=1,
            alpha=0.25,
        )

    if all([lat is not None, lon is not None]):
        plt.scatter(lon, lat, color='blue', s=10, label='My Point', zorder=10)

def _plot_xr_field(
    ax,
    da: xr.DataArray,
    *,
    lat_min_to_show=None,
    lat_max_to_show=None,
    lon_min_to_show=None,
    lon_max_to_show=None,
    cmap="Reds",
    vmin=0,
    vmax=100,
):
    lat_slice = slice(lat_min_to_show, lat_max_to_show)
    lon_slice = slice(lon_min_to_show, lon_max_to_show)

    return ax.pcolormesh(
        da.lon[lon_slice],
        da.lat[lat_slice],
        da.values[lat_slice, lon_slice],
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )


def _prepare_model_field(
    ds_dict_models: dict,
    *,
    var: str,
    ds: str,
    boundaries_dict: dict,
    depth=None,
):
    da = ds_dict_models[var][ds]
    da = _filter_bounds_xr(da, boundaries_dict)
    da = _reduce_depth_xr(da, depth)
    return da


def _spatial_title(
    *,
    var: str,
    units: dict,
    ds: str,
    ref_ds=None,
    years=None,
    depth=None,
    location_ref_var=None,
    is_climatology=False,
    compare=False,
):
    clim = " climatology" if is_climatology else ""

    if compare:
        title = f"{var}{clim} ({units[var]}) {ds} vs {ref_ds}"
    else:
        title = f"{var}{clim} ({units[var]}) {ds}"

    if years is not None:
        title += f" - {years[0]} - {years[1]}"

    if location_ref_var is not None and location_ref_var != var:
        title += f"\n location reference = {location_ref_var}"

    if depth is not None:
        if is_depth_range(depth):
            title += f" - {depth[0]} - {depth[1]} m"
        else:
            title += f" - depth = {np.round(depth, 2)}"

    return title



def spatial_maps(
    mask_ocean,
    mask_biome,
    dataframe,
    var: str,
    ds_list: list,
    units: dict,
    location_ref_var=None,
    ref_ds="obs",
    depth: list | tuple | int = None,
    outlier_var=None,
    figsize=(6, 6),
    fontsize=20,
    lat_min_to_show=None,
    lat_max_to_show=None,
    lon_min_to_show=None,
    lon_max_to_show=None,
    vmin=0,
    vmax=100,
    s=5,
    cmap="Reds",
    plot_rmse=False,
    plot_diff=False,
    plot_r2=False,
    plot_count=False,
    calculate_rmse=False,
    calculate_r2=False,
):

    """Plot spatial maps of observation/model values or pointwise scores.

    This function works with already-selected dataframe inputs for one analysis
    region or biome. It aligns the requested variable dataframe with an optional
    location-reference variable using the common spatiotemporal keys, optionally
    filters by depth and outlier variance, and then creates one map per dataset
    in ``ds_list``.

    Parameters
    ----------
    mask_ocean : xarray-like object
        Ocean mask used as the background layer. Expected to expose ``lat`` and
        ``lon`` coordinates and support ``where``.
    mask_biome : xarray-like object
        Biome or regional mask used as the highlighted overlay. Also expected
        to contain ``lat_min``, ``lat_max``, ``lon_min``, and ``lon_max`` values
        used as default plotting bounds.
    dataframe : dict[str, pandas.DataFrame]
        Dictionary mapping variable names to dataframes. Each dataframe should
        contain the common keys ``year``, ``time``, ``lat``, ``lon``, and ``lev``,
        plus columns for the datasets requested in ``ds_list`` and ``ref_ds``.
    var : str
        Variable to plot.
    ds_list : list
        Dataset or experiment names to plot. One figure is created per entry.
    units : dict
        Mapping from variable name to display units.
    location_ref_var : str, optional
        Variable used to define the valid sampling locations. Defaults to
        ``var``.
    ref_ds : str, default "obs"
        Reference dataset used for difference, RMSE, R2, and optional summary
        score calculations.
    depth : list, tuple, int, optional
        Depth selection. A tuple/list is treated as a depth range; an integer or
        scalar is treated as an exact depth level. Score maps require a range,
        not a single integer depth.
    outlier_var : float, optional
        If given, rows in the location-reference dataframe with ``var`` greater
        than or equal to this value are removed before alignment.
    figsize, fontsize, lat_min_to_show, lat_max_to_show, lon_min_to_show,
    lon_max_to_show, vmin, vmax, s, cmap
        Plot formatting controls. If plotting bounds are omitted, they are read
        from ``mask_biome``.
    plot_rmse, plot_diff, plot_r2, plot_count : bool
        Mutually constrained map modes for plotting RMSE, difference, R2, or
        observation count instead of the mean spatial field.
    calculate_rmse, calculate_r2 : bool
        If true, append whole-sample RMSE and/or R2 values to the plot title.

    Notes
    -----
    Biome selection is assumed to happen before calling this function. This
    function only plots the data it receives.
    """

    lat_min_to_show = mask_biome["lat_min"].values if lat_min_to_show is None else lat_min_to_show
    lat_max_to_show = mask_biome["lat_max"].values if lat_max_to_show is None else lat_max_to_show
    lon_min_to_show = mask_biome["lon_min"].values if lon_min_to_show is None else lon_min_to_show
    lon_max_to_show = mask_biome["lon_max"].values if lon_max_to_show is None else lon_max_to_show

    assert sum([plot_rmse, plot_r2, plot_count]) <= 1, (
        "specify only one of plot_rmse, plot_r2 or plot_count"
    )

    comparing = any([plot_rmse, plot_r2, calculate_rmse, calculate_r2, plot_diff])
    if comparing:
        assert ref_ds is not None

    if any([plot_rmse, plot_r2, plot_count]):
        assert not isinstance(depth, int), "specify a depth range not a depth point for calculation of scores!"

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
        keys=COMMON_KEYS,
    )

    ds_dicts = _filter_depth(aligned["data"], depth)
    location_ref = _filter_depth(aligned["location_ref"], depth)

    years = (
        location_ref["year"].min(),
        location_ref["year"].max(),
    )

    for ds in ds_list:
        plt.figure(figsize=figsize)
        ax = plt.subplot(1, 1, 1)

        _plot_biome_background(
            ax,
            mask_ocean,
            mask_biome,
            lat_min_to_show=lat_min_to_show,
            lat_max_to_show=lat_max_to_show,
            lon_min_to_show=lon_min_to_show,
            lon_max_to_show=lon_max_to_show,
        )

        title = _spatial_title(
            var=var,
            units=units,
            ds=ds,
            ref_ds=ref_ds,
            years=years,
            depth=depth,
            location_ref_var=location_ref_var,
            compare=comparing,
        )

        if not any([plot_rmse, plot_r2, plot_count, plot_diff]):
            if "obs" in ds:
                to_plot = (
                    ds_dicts.groupby(["lat", "lon"])
                    .agg(mean=("obs", "mean"), count=("obs", "count"))
                    .reset_index()
                )

                scatter = ax.scatter(
                    to_plot["lon"],
                    to_plot["lat"],
                    c=to_plot["mean"],
                    s=s,
                    cmap=cmap,
                    vmin=vmin,
                    vmax=vmax,
                )

            else:
                to_plot = (
                    ds_dicts.groupby(["lat", "lon"])
                    .agg(mean=(ds, "mean"), count=(ds, "count"))
                    .reset_index()
                )

                scatter = ax.scatter(
                    to_plot["lon"],
                    to_plot["lat"],
                    c=to_plot["mean"],
                    s=s,
                    cmap=cmap,
                    vmin=vmin,
                    vmax=vmax,
                )

            plt.colorbar(scatter, label=f"{var} ({units[var]})")

        elif plot_rmse:
            scores = score_pattern(ds_dicts, ds, ref_ds, score_function="MSE")
            scatter = ax.scatter(
                scores["lon"],
                scores["lat"],
                c=np.sqrt(scores["mse"]),
                s=s,
                cmap="Reds",
                vmin=vmin,
                vmax=vmax,
            )
            plt.colorbar(scatter, label=f"RMSE ({units[var]})")

        elif plot_diff:
            diff = ds_dicts.copy()
            diff["diff"] = diff[ds] - diff[ref_ds]

            scatter = ax.scatter(
                diff["lon"],
                diff["lat"],
                c=diff["diff"],
                s=s,
                cmap="Reds",
                vmin=vmin,
                vmax=vmax,
            )
            plt.colorbar(scatter, label=f"Diff ({units[var]})")

        elif plot_r2:
            scores = score_pattern(ds_dicts, ds, ref_ds, score_function="R2")
            scatter = ax.scatter(
                scores["lon"],
                scores["lat"],
                c=scores["corr"] ** 2,
                s=s,
                cmap="Reds",
                vmin=0,
                vmax=1,
            )
            plt.colorbar(scatter, label="$R^2$")

        elif plot_count:
            counts = (
                ds_dicts.groupby(["lat", "lon"])
                .agg(count=(ds, "count"))
                .reset_index()
            )

            scatter = ax.scatter(
                counts["lon"],
                counts["lat"],
                c=counts["count"],
                s=s,
                vmin=0,
                vmax=vmax,
                cmap="Blues_r",
            )
            plt.colorbar(scatter, label="count")

        if calculate_r2 or calculate_rmse:
            r2 = np.round(corr(ds_dicts[ref_ds], ds_dicts[ds]) ** 2, 2)
            rmse_score = np.round(rmse(ds_dicts[ref_ds], ds_dicts[ds]), 2)

            if calculate_r2:
                title += f"\n $R^2$ score {r2}"

            if calculate_rmse:
                sep = " - " if calculate_r2 else "\n "
                title += f"{sep}RMSE {rmse_score}"

        ax.set_ylabel("lat")
        ax.set_xlabel("lon")
        ax.set_title(title)


def spatial_maps_climatology(
    mask_ocean,
    mask_biome,
    ds_dict,
    var: str,
    ds_list: list,
    units: dict,
    ref_ds="obs",
    depth: list | tuple | int = None,
    figsize=(6, 6),
    fontsize=20,
    lat_min_to_show=None,
    lat_max_to_show=None,
    lon_min_to_show=None,
    lon_max_to_show=None,
    vmin=0,
    vmax=100,
    cmap="Reds",
    plot_rmse=False,
    plot_diff=False,
    calculate_rmse=False,
):
    
    """Plot spatial climatology maps from gridded xarray fields.

    This function creates one spatial map per dataset in ``ds_list`` using
    gridded fields from ``ds_dict``. It applies the plotting bounds stored on
    ``mask_biome``, optionally filters or averages over depth, and can instead
    plot RMSE or difference maps against ``ref_ds``.

    Parameters
    ----------
    mask_ocean : xarray-like object
        Ocean mask used as the background layer. Expected to expose ``lat`` and
        ``lon`` coordinates and support ``where``.
    mask_biome : xarray-like object
        Biome or regional mask used as the highlighted overlay. Also expected
        to contain ``lat_min``, ``lat_max``, ``lon_min``, and ``lon_max`` values
        used as default plotting bounds.
    ds_dict : dict
        Nested dictionary of gridded fields, accessed as
        ``ds_dict[var][dataset_name]``. Each field is expected to be an
        ``xarray.DataArray`` with latitude and longitude coordinates and, when
        applicable, a ``lev`` dimension.
    var : str
        Variable to plot.
    ds_list : list
        Dataset or experiment names to plot. One figure is created per entry.
    units : dict
        Mapping from variable name to display units.
    ref_ds : str, default "obs"
        Reference dataset used when plotting RMSE, plotting differences, or
        calculating RMSE for the title.
    depth : list, tuple, int, optional
        Depth selection. A tuple/list is treated as a depth range and reduced
        over ``lev`` by the existing depth helper; a scalar is treated as an
        exact depth level. RMSE and difference maps require a range rather than
        a single integer depth.
    figsize, fontsize, lat_min_to_show, lat_max_to_show, lon_min_to_show,
    lon_max_to_show, vmin, vmax, cmap
        Plot formatting controls. If plotting bounds are omitted, they are read
        from ``mask_biome``.
    plot_rmse, plot_diff : bool
        Plot RMSE or reference-minus-model differences instead of the model
        climatology field.
    calculate_rmse : bool
        If true, append the whole-field RMSE against ``ref_ds`` to the title.

    Notes
    -----
    This function assumes climatological gridded fields are already prepared in
    ``ds_dict``. It does not perform biome selection beyond the mask bounds.
    """
    comparing = any([plot_rmse, calculate_rmse, plot_diff])

    boundaries_dict = dict(
        lat_min = mask_biome["lat_min"].values,
        lat_max = mask_biome["lat_max"].values,
        lon_min = mask_biome["lon_min"].values,
        lon_max = mask_biome["lon_max"].values
    )

    lat_min_to_show = boundaries_dict["lat_min"] if lat_min_to_show is None else lat_min_to_show
    lat_max_to_show = boundaries_dict["lat_max"] if lat_max_to_show is None else lat_max_to_show
    lon_min_to_show = boundaries_dict["lon_min"] if lon_min_to_show is None else lon_min_to_show
    lon_max_to_show = boundaries_dict["lon_max"] if lon_max_to_show is None else lon_max_to_show

    if comparing:
        assert ref_ds is not None

    if any([plot_rmse, plot_diff]) and depth is not None:
        assert not isinstance(depth, int), "specify a depth range not a depth point for score calculations!"

    ds_ref = None
    if ref_ds is not None:
        ds_ref = _prepare_model_field(
            ds_dict,
            var=var,
            ds=ref_ds,
            boundaries_dict=boundaries_dict,
            depth=depth,
        )

    for ds in ds_list:
        plt.figure(figsize=figsize)
        ax = plt.subplot(1, 1, 1)

        _plot_biome_background(
            ax,
            mask_ocean,
            mask_biome,
            lat_min_to_show=lat_min_to_show,
            lat_max_to_show=lat_max_to_show,
            lon_min_to_show=lon_min_to_show,
            lon_max_to_show=lon_max_to_show,
        )

        title = _spatial_title(
            var=var,
            units=units,
            ds=ds,
            ref_ds=ref_ds,
            depth=depth,
            is_climatology=True,
            compare=comparing,
        )

        ds_model = _prepare_model_field(
            ds_dict,
            var=var,
            ds=ds,
            boundaries_dict=boundaries_dict,
            depth=depth,
        )

        if calculate_rmse:
            rmse_score = np.round(rmse(ds_ref, ds_model).values, 2)

        if plot_rmse:
            ds_model_to_plot = np.sqrt(((ds_ref - ds_model) ** 2).mean("lev"))

        elif plot_diff:
            ds_model_to_plot = ds_ref.mean("lev") - ds_model.mean("lev")

        else:
            ds_model_to_plot = ds_model

        scatter = _plot_xr_field(
            ax,
            ds_model_to_plot,
            lat_min_to_show=lat_min_to_show,
            lat_max_to_show=lat_max_to_show,
            lon_min_to_show=lon_min_to_show,
            lon_max_to_show=lon_max_to_show,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )

        plt.colorbar(scatter, label=f"{var} ({units[var]})")

        if calculate_rmse:
            title += f"\n RMSE {rmse_score}"

        ax.set_ylabel("lat")
        ax.set_xlabel("lon")
        ax.set_title(title)