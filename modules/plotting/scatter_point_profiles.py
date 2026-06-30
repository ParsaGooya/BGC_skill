
import random
import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt

from modules.plotting.utils import (
    COMMON_KEYS,
    _align_dataframes, 
    _filter_depth,
    _filter_location,
    _filter_depth_xr,
    markers, 
    rmse, 
    corr
)

def _make_shapes_dict(ds_list, shapes_dict=None):
    """Return marker styles for each dataset, using a provided mapping if available."""
    if shapes_dict is not None:
        return shapes_dict

    assert len(ds_list) <= len(markers), (
        "Not enough marker shapes for this many datasets"
    )

    shapes = random.sample(markers, len(ds_list))
    return {ds: shapes[i] for i, ds in enumerate(ds_list)}


def _get_plot_color(colors_dict, ds):
    """Return the plotting color for a dataset, or choose a random CSS color."""
    if colors_dict is None:
        return random.choice(list(colors.CSS4_COLORS.values()))
    return colors_dict[ds]




def _plot_profile_points(
    ax,
    x,
    lev,
    *,
    ds,
    color,
    marker,
    label,
    cmap="tab10",
):
    """Plot one profile as variable values against depth on the given axes."""
    return ax.scatter(
        x,
        lev,
        c=color,
        marker=marker,
        cmap=cmap,
        label=label,
    )


def _profile_title(
    *,
    var,
    time_title,
    lat,
    lon,
    max_dist=None,
    location_ref_var=None,
    mean_lev=False,
):
    """Build the title used for point-profile plots."""
    if max_dist is not None:
        title = (
            f"{var} - {time_title}\n"
            f"within {max_dist} km from lat: {lat} lon: {lon}"
        )
    else:
        title = f"{var} - {time_title} - lat: {lat} lon: {lon}"

    if location_ref_var is not None:
        title += f"\n location reference = {location_ref_var}"

    if mean_lev:
        title += " (mean per lev)"

    return title



def point_profile(
    dataframe,
    var: str,
    ds_list: list,
    units: dict,
    location_ref_var=None,
    ref_ds="obs",
    depth: list | tuple = None,
    outlier_var=None,
    figsize=(6, 6),
    fontsize=20,
    lat=None,
    lon=None,
    xlims=None,
    ylims=None,
    max_dist=None,
    isolate_time=False,
    mean_lev=False,
    legend=True,
    bbox_to_anchor=(1, 0.5),
    shapes_dict=None,
    colors_dict=None,
    calculate_rmse=False,
    calculate_r2=False,
):
    """Plot vertical profiles from matched point-sample dataframes.

    The input dataframe is expected to map variable names to pandas DataFrames.
    Each variable dataframe should contain the shared matching columns in
    COMMON_KEYS, typically year, time, lat, lon, and lev, plus one column per
    dataset or experiment listed in ds_list.

    The function aligns the requested variable dataframe with a location
    reference dataframe before plotting. If location_ref_var is not provided,
    the plotted variable itself is used as the location reference. This keeps
    only samples that are available in both the plotted variable and the
    location reference. The aligned data can be filtered by depth range,
    outlier variance, and geographic location.

    Location selection
    ------------------
    If lat and lon are both None, the function chooses the lat/lon pair with
    the largest number of obs samples. If lat and lon are provided and max_dist
    is None, the nearest available sample location is selected. If max_dist is
    provided, all samples within max_dist kilometers of the requested location
    are retained.

    Time and depth selection
    ------------------------
    If isolate_time is True, the year/time pair with the largest number of
    samples at the selected location is retained. If mean_lev is True, the
    filtered dataframe is averaged by lev before plotting. If depth is provided,
    the dataframe is filtered with _filter_depth before location and plotting
    operations.

    Plot behavior
    -------------
    Each dataset in ds_list is plotted as variable value on the x-axis against
    lev on the y-axis. Observations use an X marker, while other datasets use
    markers from shapes_dict or automatically assigned markers. Colors come
    from colors_dict when provided; otherwise random CSS colors are used.

    Optional scores
    ---------------
    If calculate_rmse and/or calculate_r2 are True, scores are computed against
    ref_ds using the filtered data and appended to legend labels for each
    non-observation dataset.

    Parameters
    ----------
    dataframe : dict
        Mapping from variable name to point-sample dataframe.
    var : str
        Variable to plot.
    ds_list : list
        Dataset or experiment columns to plot.
    units : dict
        Mapping from variable name to display units.
    location_ref_var : str, optional
        Variable used to define common sample locations. Defaults to var.
    ref_ds : str, default "obs"
        Reference dataset used for optional RMSE and R² calculations.
    depth : list or tuple, optional
        Depth range passed to _filter_depth.
    outlier_var : float, optional
        Maximum allowed variance in the location-reference dataframe.
    figsize : tuple, default (6, 6)
        Matplotlib figure size.
    fontsize : int, default 20
        Kept for API compatibility; not directly used in the current body.
    lat, lon : float, optional
        Requested location. If omitted, the most sampled obs location is used.
    xlims, ylims : tuple, optional
        Axis limits passed to ax.set_xlim and ax.set_ylim.
    max_dist : float, optional
        Search radius in kilometers around lat/lon.
    isolate_time : bool, default False
        Whether to keep only the most sampled year/time pair.
    mean_lev : bool, default False
        Whether to average values by lev before plotting.
    legend : bool, default True
        Whether to draw the legend.
    bbox_to_anchor : tuple, default (1, 0.5)
        Legend anchor passed to ax.legend.
    shapes_dict : dict, optional
        Mapping from dataset name to marker style.
    colors_dict : dict, optional
        Mapping from dataset name to color.
    calculate_rmse, calculate_r2 : bool, default False
        Whether to append RMSE and/or R² scores to labels.
    """
    if calculate_rmse or calculate_r2:
        assert ref_ds is not None

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

    if lat is None and lon is None:
        max_dist = None

        counts = (
            ds_dicts.groupby(["lat", "lon"])
            .agg(count=("obs", "count"))
            .reset_index()
        )

        row = counts.loc[counts["count"].idxmax()]
        lat = row["lat"]
        lon = row["lon"]

        ds_dicts_toplot = ds_dicts[
            (ds_dicts["lat"] == lat) & (ds_dicts["lon"] == lon)
        ]

    else:
        ds_dicts_toplot = _filter_location(
            ds_dicts,
            lat=lat,
            lon=lon,
            max_dist_km=max_dist,
        )

        if max_dist is None:
            lat = ds_dicts_toplot["lat"].iloc[0]
            lon = ds_dicts_toplot["lon"].iloc[0]

    if isolate_time:
        counts = (
            ds_dicts_toplot.groupby(["year", "time"])
            .size()
            .reset_index(name="count")
        )

        row = counts.loc[counts["count"].idxmax()]
        year_sel = row["year"]
        time_sel = row["time"]

        ds_dicts_toplot = ds_dicts_toplot[
            (ds_dicts_toplot["year"] == year_sel)
            & (ds_dicts_toplot["time"] == time_sel)
        ]

        time_title = f"year {year_sel} month {time_sel}"

    else:
        time_title = (
            f"{location_ref['year'].min()} - "
            f"{location_ref['year'].max()}"
        )

    if mean_lev:
        ds_dicts_toplot = (
            ds_dicts_toplot.groupby(["lev"], as_index=False)
            .mean(numeric_only=True)
        )

    shapes_dict = _make_shapes_dict(ds_list, shapes_dict)

    fig = plt.figure(figsize=figsize)
    ax = plt.subplot(1, 1, 1)

    title = _profile_title(
        var=var,
        time_title=time_title,
        lat=lat,
        lon=lon,
        max_dist=max_dist,
        location_ref_var=(
            location_ref_var if location_ref_var != var else None
        ),
        mean_lev=mean_lev,
    )

    for ds in ds_list:
        label = ds

        if calculate_rmse or calculate_r2:
            r2 = np.round(corr(ds_dicts_toplot[ref_ds], ds_dicts_toplot[ds]) ** 2, 2)
            rmse_score = np.round(rmse(ds_dicts_toplot[ref_ds], ds_dicts_toplot[ds]), 2)

        if calculate_r2 and ds != "obs":
            label += f" - $R^2$ score (vs {ref_ds}) {r2}"

        if calculate_rmse and ds != "obs":
            if calculate_r2:
                label += f" - RMSE {rmse_score}"
            else:
                label += f" - RMSE (vs {ref_ds}) {rmse_score}"


        color = _get_plot_color(colors_dict, ds)

        marker = "X" if ds == "obs" else shapes_dict[ds]

        scatter = _plot_profile_points(
            ax,
            ds_dicts_toplot[ds],
            ds_dicts_toplot["lev"],
            ds=ds,
            color=color,
            marker=marker,
            label=label,
        )

    if legend:
        ax.legend(loc="center left", bbox_to_anchor=bbox_to_anchor)

    ax.set_ylabel("depth")
    ax.set_xlabel(f"{var} ({units[var]})")
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_title(title)
    ax.invert_yaxis()




def point_profile_climatology(
    ds_dict,
    var: str,
    ds_list: list,
    units: dict,
    lat,
    lon,
    ref_ds="obs",
    depth: list | tuple = None,
    figsize=(6, 6),
    fontsize=20,
    xlims=None,
    ylims=None,
    legend=False,
    bbox_to_anchor=(1, 0.5),
    shapes_dict=None,
    colors_dict=None,
    calculate_rmse=False,
):
    """Plot vertical climatological profiles from gridded xarray fields.

    The input ds_dict is expected to be a nested dictionary of the form
    ds_dict[var][dataset], where each entry is an xarray DataArray with lat,
    lon, and lev coordinates. For each dataset in ds_list, the function
    optionally filters the vertical coordinate, selects the nearest grid point
    to the requested lat/lon, and plots the resulting vertical profile.

    Reference handling
    ------------------
    If calculate_rmse is True, ref_ds must be provided. The reference profile is
    selected from ds_dict[var][ref_ds] using the same depth filtering and
    nearest-location selection as the plotted profiles. RMSE values are appended
    to legend labels for non-observation datasets.

    Plot behavior
    -------------
    Each dataset is plotted as climatological variable value on the x-axis
    against lev on the y-axis. Observations use an X marker, while other
    datasets use markers from shapes_dict or automatically assigned markers.
    Colors come from colors_dict when provided; otherwise random CSS colors are
    used.

    Parameters
    ----------
    ds_dict : dict
        Nested mapping from variable name to dataset name to xarray DataArray.
    var : str
        Variable to plot.
    ds_list : list
        Dataset or experiment names to plot.
    units : dict
        Mapping from variable name to display units.
    lat, lon : float
        Requested location; the nearest grid point is selected.
    ref_ds : str, default "obs"
        Reference dataset used for optional RMSE calculations.
    depth : list or tuple, optional
        Depth range passed to _filter_depth_xr.
    figsize : tuple, default (6, 6)
        Matplotlib figure size.
    fontsize : int, default 20
        Kept for API compatibility; not directly used in the current body.
    xlims, ylims : tuple, optional
        Axis limits passed to ax.set_xlim and ax.set_ylim.
    legend : bool, default False
        Whether to draw the legend.
    bbox_to_anchor : tuple, default (1, 0.5)
        Legend anchor passed to ax.legend.
    shapes_dict : dict, optional
        Mapping from dataset name to marker style.
    colors_dict : dict, optional
        Mapping from dataset name to color.
    calculate_rmse : bool, default False
        Whether to append RMSE scores against ref_ds to labels.
    """
    if calculate_rmse:
        assert ref_ds is not None

    shapes_dict = _make_shapes_dict(ds_list, shapes_dict)

    ds_ref = None
    if ref_ds is not None:
        ds_ref = ds_dict[var][ref_ds]
        ds_ref = _filter_depth_xr(ds_ref, depth)
        ds_ref = ds_ref.sel(lat=lat, lon=lon, method="nearest")

    fig = plt.figure(figsize=figsize)
    ax = plt.subplot(1, 1, 1)

    title = f"{var} - climatology - lat: {lat} lon: {lon}"

    for ds in ds_list:
        ds_model = ds_dict[var][ds]
        ds_model = _filter_depth_xr(ds_model, depth)
        ds_model = ds_model.sel(lat=lat, lon=lon, method="nearest")

        label = ds

        if calculate_rmse and ds != "obs":
            rmse_score = np.round(rmse(ds_model, ds_ref).values, 2)
            label += f" - RMSE (vs {ref_ds}) {rmse_score}"

        color = _get_plot_color(colors_dict, ds)
        marker = "X" if ds == "obs" else shapes_dict[ds]

        _plot_profile_points(
            ax,
            ds_model,
            ds_model["lev"],
            ds=ds,
            color=color,
            marker=marker,
            label=label,
        )

    if legend:
        ax.legend(loc="center left", bbox_to_anchor=bbox_to_anchor)

    ax.set_ylabel("depth")
    ax.set_xlabel(f"{var} ({units[var]})")
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_title(title)
    ax.invert_yaxis()