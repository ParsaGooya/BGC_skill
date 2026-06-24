
import random
import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt

from modules.plotting.utils import (
    COMMON_KEYS,
    _align_dataframes, 
    _filter_depth,
    _filter_location,
    markers, 
    rmse, 
    corr
)

def _make_shapes_dict(ds_list, shapes_dict=None):
    if shapes_dict is not None:
        return shapes_dict

    assert len(ds_list) <= len(markers), (
        "Not enough marker shapes for this many datasets"
    )

    shapes = random.sample(markers, len(ds_list))
    return {ds: shapes[i] for i, ds in enumerate(ds_list)}


def _get_plot_color(colors_dict, ds):
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

        time_title = f"year {year_sel} month {time_sel + 1}"

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




