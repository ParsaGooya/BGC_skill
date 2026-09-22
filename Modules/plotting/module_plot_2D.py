import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import ticker
import cartopy.crs as ccrs
from pathlib import Path
import cartopy
import cartopy.crs as ccrs
from numpy.typing import NDArray
from mpl_toolkits.axes_grid1 import make_axes_locatable
import xarray as xr
from typing import Sequence, Literal

from Modules.analysis.module_global_averages import area_weighted_avg
from Modules.analysis.module_data_postprocessing import (spatial_mask, 
                                                         Metrics, 
                                                         calculate_measure)
from Modules.data_info.module_state_dict import state_dict
from Modules.plotting.utils import *

    

def seasonal_pattern_correlation(
    reference: xr.DataArray,
    target: xr.DataArray,
    season: str,
    ldyr_ini: int = 0,
    ldyr_end: int = 1,
    spatial_dims: tuple[str, ...] = ("lat", "lon"),
) -> xr.DataArray:
    """
    Calculate spatial pattern correlation for a requested season.

    If temporal dimensions are available:
        1. form the seasonal mean independently for each year,
        2. calculate spatial pattern correlation for each year,
        3. average the correlations over years.

    If no month dimension exists, pattern correlation is calculated
    directly from the provided fields.
    """
    reference, target = xr.align(
        reference,
        target,
        join="inner",
    )

    if "month" in reference.dims:
        reference = seasonal_mean(
            reference,
            season=season,
            ldyr_ini=ldyr_ini,
            ldyr_end=ldyr_end,
        )

        target = seasonal_mean(
            target,
            season=season,
            ldyr_ini=ldyr_ini,
            ldyr_end=ldyr_end,
        )

    corr = xr.corr(
        reference,
        target,
        dim=list(spatial_dims),
    )

    if "year" in corr.dims:
        corr = corr.mean("year")

    return corr


def prepare_snapshot(
    data: xr.DataArray,
    season: str = "ANN",
    ldyr: int = 0,
    years_to_plot: int | list[int] | np.ndarray | tuple[int, int] | None = None,
    return_year_string: bool = False,
) -> xr.DataArray:
    """
    Prepare a climatological snapshot.

    Parameters
    ----------
    data
        Input data with dimensions including some combination of
        year, month, lev, and lon.
    season
        Season or month to average. This code be individual months,
        normal seasons (e.g. JFM, AMJ, ...) or shifted seasons (DJF, MAM, ...).
    ldyr
        Lead year.
    years_to_plot
        Year selection.

        - None:
            Average over all available years.
        - int:
            Select one year.
        - list/array:
            Composite mean over selected years.
        - tuple[int, int]:
            Difference between two years, interpreted as
            (year1, year0) -> year1 - year0.
    return_year_string
        whether or not to select the year string for plot title.

    Returns
    -------
    xr.DataArray
        Data reduced to depth x longitude.
    year_string
        string reprenseting the chosen years for plot title.
    """
    year_string = None
    # --------------------------------------------------------------
    # Seasonal mean
    # --------------------------------------------------------------
    if "month" in data.dims:
        data = seasonal_mean(
            data,
            season=season,
            ldyr_ini=ldyr,
            ldyr_end=ldyr + 1,
        )

    # --------------------------------------------------------------
    # Year selection
    # --------------------------------------------------------------
    if "year" not in data.dims:
        year_string = None
        out = data.squeeze()

        if return_year_string:
            return out, year_string
        
        return out

    if years_to_plot is None:
        y0 = data.year.min().item()
        y1 = data.year.max().item()
        year_string = f'{y0} - {y1} mean'
        out = data.mean("year").squeeze()

    # Explicit year difference: (year1, year0)
    elif isinstance(years_to_plot, tuple):
        if len(years_to_plot) != 2:
            raise ValueError(
                "A year-difference tuple must contain exactly two years."
            )

        year1, year0 = years_to_plot
        year_string = f'{year1} minus {year0}'
        out = (
            data.sel(year=year1)
            - data.sel(year=year0)
        ).squeeze()

    # One year
    elif np.isscalar(years_to_plot):
        year_string = f'{years_to_plot}'
    
        out =  data.sel(
            year=years_to_plot,
        ).squeeze()

    else:

        # Composite over several years
        out = (
            data.sel(year=years_to_plot)
            .mean("year")
            .squeeze()
        )

        years = data.sel(year=years_to_plot).year.values

        year_string = (
            f"{years.min():g}-{years.max():g} composite "
            f"(n={years.size})"
        )


    if return_year_string:
        return out, year_string

    return out


def select_depth_range(
    data: xr.DataArray,
    lev_interp: np.ndarray | None = None,
    lev_range: float | tuple[float, float] | None = None,
    return_range_string: bool = False,
) -> xr.DataArray:

    if "lev" not in data.dims:
        raise ValueError(
            "The data must contain a 'lev' dimension."
        )

    if lev_interp is not None:
        data = data.interp(
            lev=lev_interp,
        )

    if lev_range is not None:
        if np.isscalar(lev_range):
            lev_min = 0
            lev_max = lev_range
        else:
            lev_min, lev_max = lev_range

        data = data.where(
            (data.lev >= lev_min)
            & (data.lev <= lev_max),
            drop=True,
        )

    if return_range_string:
        lev_min = data.lev.min().item()
        lev_max = data.lev.max().item()

        range_string = f"{lev_min:.2f}-{lev_max:.2f} m"

        return data, range_string

    return data


def plot_composites(
    ds_list: list[Exp],
    data_dict: dict[Exp, state_dict],
    mask=None,
    specific_years: list | None = None,
    figsize=(12, 12),
    depth_range: float | Sequence[float] = None,
    central_longitude=260,
    ldyr_ini=0,
    ldyr_end=1,
    vmax=2,
    vmin=-2,
    cmap="RdBu_r",
    cbar_label=r"mol m$^{-2}$ yr$^{-1}$",
    std=False,
    seasons_to_plot=("ANN",),
    fontsize=20,
    var_name="",
    show_equator=False,
    dir_name=None,
    file_name=None,
    save=False,
):

    """
    Plot seasonal spatial composites for multiple datasets.

    The function creates a grid of global maps with one row per requested
    season and one column per dataset. For each dataset, data can optionally
    be restricted to specific years and a specified depth range. Data with a
    depth dimension are averaged over depth before plotting.

    When an observational dataset is available under the ``"obs"`` key in
    ``data_dict``, model and observational years are aligned before computing
    seasonal pattern correlations. Each panel reports the area-weighted
    spatial average and, when observations are available, the corresponding
    pattern correlation.

    Parameters
    ----------
    ds_list : list[Exp]
        Dataset or experiment names to plot. Each entry must correspond to a
        key in ``data_dict``.
    data_dict : dict[Exp, state_dict]
        Mapping from dataset or experiment names to state objects containing
        the data to plot. An optional ``"obs"`` entry is used as the
        observational reference for pattern-correlation calculations.
    mask : optional
        Spatial mask used when calculating the area-weighted global average.
        If ``None``, a mask is generated from each plotted field using
        ``spatial_mask``.
    specific_years : list or None, optional
        Years to include in the analysis. If ``None``, all available years
        are used.
    figsize : tuple, optional
        Figure size passed to ``matplotlib.pyplot.figure``.
    depth_range : float or Sequence[float] or None, optional
        Depth or depth range to select before averaging over the ``lev``
        dimension. If ``None``, no explicit depth selection is applied.
    central_longitude : float, optional
        Central longitude of the Robinson map projection.
    ldyr_ini : int, optional
        Initial lead year passed to ``seasonal_mean`` and
        ``seasonal_pattern_correlation``.
    ldyr_end : int, optional
        Final lead year passed to ``seasonal_mean`` and
        ``seasonal_pattern_correlation``.
    vmax : float, optional
        Maximum value of the plotting color scale.
    vmin : float, optional
        Minimum value of the plotting color scale.
    cmap : str, optional
        Matplotlib colormap used for the spatial fields.
    cbar_label : str, optional
        Label for the shared colorbar.
    std : bool, optional
        If ``True``, plot the interannual standard deviation instead of the
        temporal mean.
    seasons_to_plot : Sequence[str], optional
        Seasons to plot. A separate row is created for each season.
    fontsize : float, optional
        Font size used for panel and row titles.
    var_name : str, optional
        Variable name included in row labels.
    show_equator : bool, optional
        If ``True``, draw a dotted line along the equator.
    dir_name : str or Path or None, optional
        Directory in which to save the figure when ``save`` is ``True``.
    file_name : str or None, optional
        Output filename without the ``.png`` extension.
    save : bool, optional
        If ``True``, save the figure as a PNG file.

    Raises
    ------
    ValueError
        If ``std`` is ``True`` and fewer than two ``specific_years`` are
        provided.
    RuntimeError
        If the observational data contain a depth dimension that is averaged
        internally but a corresponding model dataset does not contain a
        ``lev`` dimension.

    Notes
    -----
    If a ``month`` dimension is present, seasonal means are calculated using
    ``seasonal_mean``. If a ``year`` dimension is subsequently present, the
    function averages over years unless ``std=True``.

    Pattern correlations are calculated only when observational data are
    provided. The global average shown in each panel is calculated with
    ``area_weighted_avg``.

    Documentation produced with the assistance of AI.
    """

    if std and specific_years is not None and len(specific_years) < 2:
        raise ValueError(
            "At least two years are required to calculate interannual std."
        )

    has_obs = "obs" in data_dict
    _obs_depth_averaged = False

    fig = plt.figure(figsize=figsize)

    for season_idx, season in enumerate(seasons_to_plot):

        # --------------------------------------------------------------
        # Prepare observational reference
        # --------------------------------------------------------------
        obs_ref = None

        if has_obs:
            obs_ref = data_dict["obs"].data

            if specific_years is not None and "year" in obs_ref.dims:
                obs_ref = obs_ref.sel(year=specific_years)

            if depth_range is not None and "lev" in obs_ref.dims:
                obs_ref, _ = resolve_depth_range(obs_ref, depth_range)
                

            if "lev" in obs_ref.dims:
                obs_ref = obs_ref.mean("lev")
                _obs_depth_averaged = True

        # --------------------------------------------------------------
        # Plot each dataset
        # --------------------------------------------------------------
        for ds_idx, name in enumerate(ds_list):

            ax = plt.subplot(
                len(seasons_to_plot),
                len(ds_list),
                season_idx * len(ds_list) + ds_idx + 1,
                projection=ccrs.Robinson(
                    central_longitude=central_longitude
                ),
            )

            ds = data_dict[name].data

            if _obs_depth_averaged and "lev" not in ds.dims:
                raise RuntimeError(
                    "Observation data has depth dimension but the model data does not."
                )

            # ----------------------------------------------------------
            # Align years with observations
            # ----------------------------------------------------------
            if (
                obs_ref is not None
                and "year" in ds.dims
                and "year" in obs_ref.dims
            ):
                ds, obs_aligned = xr.align(
                    ds,
                    obs_ref,
                    join="inner",
                )
            else:
                obs_aligned = obs_ref

            if specific_years is not None and "year" in ds.dims:
                ds = ds.sel(year=specific_years)

            # ----------------------------------------------------------
            # Depth selection
            # ----------------------------------------------------------
            selected_depth = None

            if depth_range is not None and "lev" in ds.dims:
                ds, selected_depth = resolve_depth_range(ds, depth_range)

            if "lev" in ds.dims:
                ds = ds.mean("lev")
            # ----------------------------------------------------------
            # Build field to plot
            # ----------------------------------------------------------
            if "month" in ds.dims:
                ds_seasonal = seasonal_mean(
                    ds,
                    season=season,
                    ldyr_ini=ldyr_ini,
                    ldyr_end=ldyr_end,
                )
            else:
                ds_seasonal = ds

            if "year" in ds_seasonal.dims:
                if std:
                    ds_toplot = ds_seasonal.std("year")
                else:
                    ds_toplot = ds_seasonal.mean("year")
            else:
                ds_toplot = ds_seasonal

            # ----------------------------------------------------------
            # Plot
            # ----------------------------------------------------------
            if central_longitude != 0:
                plot_data = add_cyclic_point(ds_toplot)
            else:
                plot_data = ds_toplot

            cb = ax.pcolormesh(
                plot_data.lon,
                plot_data.lat,
                plot_data,
                cmap=plt.get_cmap(cmap),
                vmax=vmax,
                vmin=vmin,
                rasterized=True,
                transform=ccrs.PlateCarree(),
            )

            if show_equator:
                ax.plot(
                    ds_toplot.lon,
                    np.zeros(len(ds_toplot.lon)),
                    color="black",
                    linewidth=1,
                    linestyle="dotted",
                    transform=ccrs.PlateCarree(),
                )

            ax.coastlines()
            ax.set_ylabel("")
            ax.set_xlabel("")

            # ----------------------------------------------------------
            # Global average
            # ----------------------------------------------------------
            plot_mask = mask
            if plot_mask is None:
                plot_mask = spatial_mask(ds_toplot)

            glbavg = np.round(
                area_weighted_avg(
                    ds_toplot,
                    mask=plot_mask,
                ).values,
                4,
            )

            # ----------------------------------------------------------
            # Pattern correlation
            # ----------------------------------------------------------
            if obs_aligned is not None:

                corr_pat = seasonal_pattern_correlation(
                    obs_aligned,
                    ds,
                    season=season,
                    ldyr_ini=ldyr_ini,
                    ldyr_end=ldyr_end,
                )

                corr_pat = np.round(corr_pat.values, 2)

                panel_title = f"{glbavg}, {corr_pat}"

            else:
                panel_title = f"{glbavg}"

            if season_idx == 0:
                panel_title = f"{name}\n{panel_title}"

            ax.set_title(
                panel_title,
                fontsize=fontsize,
            )

            # ----------------------------------------------------------
            # Row label
            # ----------------------------------------------------------
            if ds_idx == 0:


                row_title = f"Composite {var_name} {season}"

                if selected_depth is not None:
                    row_title = (
                        f"{var_name} {season} "
                        f"lev: {selected_depth}"
                    )
                
                if _obs_depth_averaged:
                    row_title += " depth average"

                if "year" in ds_toplot.dims:
                    y0 = ds_toplot.year.values[0]
                    y1 = ds_toplot.year.values[-1]

                    row_title += f" ({y0}–{y1})"


                ax.text(
                    -0.25,
                    1.5,
                    row_title,
                    fontsize=fontsize,
                    transform=ax.transAxes,
                )

    # ------------------------------------------------------------------
    # Colorbar
    # ------------------------------------------------------------------
    divider = make_axes_locatable(ax)

    ax_cb = divider.append_axes(
        "bottom",
        size="10%",
        pad=0.1,
        axes_class=plt.Axes,
    )

    cbar = plt.colorbar(
        cb,
        cax=ax_cb,
        orientation="horizontal",
    )

    if std:
        cbar_label = f"std ({cbar_label})"

    cbar.set_label(
        label=cbar_label,
        size=20,
    )
    cbar.ax.tick_params(labelsize=20)

    plt.tight_layout()
    plt.subplots_adjust(
        wspace=0.1,
        hspace=0.3,
    )

    if save:
        Path(dir_name).mkdir(
            parents=True,
            exist_ok=True,
        )

        plt.savefig(
            f"{dir_name}/{file_name}.png"
        )
        
        
def plot_measures(
    ds_list: list[Exp],
    data_dict: dict[Exp, state_dict],
    measure: Metrics = "rmse",
    figsize=(12, 12),
    central_longitude=260,
    ldyr_ini=0,
    ldyr_end=1,
    vmax=2,
    vmin=-2,
    label="",
    cmap="RdBu_r",
    dir_name=None,
    file_name=None,
    var_name="",
    fontsize=20,
    depth_range: float | Sequence[float] = None,
    individual_months=False,
    shifted_seasons=False,
    mask=None,
    save=False,
):

    _has_depth = False

    if individual_months:
        seasons = MONTH_NAMES

    elif shifted_seasons:
        seasons = ("DJF", "MAM", "JJA", "SON", "ANN")

    else:
        seasons = ("JFM", "AMJ", "JAS", "OND", "ANN")

    fig = plt.figure(figsize=figsize)

    obs = data_dict["obs"].data

    for season_idx, season in enumerate(seasons):

        for ds_idx, name in enumerate(ds_list):

            ax = plt.subplot(
                len(seasons),
                len(ds_list),
                season_idx * len(ds_list) + ds_idx + 1,
                projection=ccrs.Robinson(
                    central_longitude=central_longitude
                ),
            )

            target = data_dict[name].data

            obs_aligned, target_aligned = xr.align(
                obs,
                target,
                join="inner",
            )

            selected_depth = None
            if (depth_range is not None and 
             "lev" in obs_aligned.dims and 
             "lev" in target_aligned.dims):
                    obs_aligned, _ = resolve_depth_range(obs_aligned, depth_range)
                    target_aligned, selected_depth = resolve_depth_range(target_aligned, depth_range)
                    
            if "lev" in obs_aligned.dims:
                if "lev" not in target_aligned.dims:
                    raise RuntimeError(
                        "Observation has depth dimension but the model data does not."
                    )
                _has_depth = True

            obs_seasonal = seasonal_mean(
                obs_aligned,
                season=season,
                ldyr_ini=ldyr_ini,
                ldyr_end=ldyr_end,
            )

            target_seasonal = seasonal_mean(
                target_aligned,
                season=season,
                ldyr_ini=ldyr_ini,
                ldyr_end=ldyr_end,
            )

            ds_toplot = calculate_measure(
                obs_seasonal,
                target_seasonal,
                measure=measure,
                dim="year",
            )

            plot_mask = mask

            if plot_mask is None:
                plot_mask = spatial_mask(ds_toplot)

            glbavg = np.round(
                area_weighted_avg(
                    ds_toplot,
                    mask=plot_mask,
                    is_ds=False,
                ).values,
                2,
            )

            if central_longitude != 0:
                plot_data = add_cyclic_point(ds_toplot)
            else:
                plot_data = ds_toplot

            cb = ax.pcolormesh(
                plot_data.lon,
                plot_data.lat,
                plot_data,
                cmap=plt.get_cmap(cmap),
                vmax=vmax,
                vmin=vmin,
                rasterized=True,
                transform=ccrs.PlateCarree(),
            )

            ax.coastlines()
            ax.set_ylabel("")
            ax.set_xlabel("")

            panel_title = f"{glbavg}"

            if season_idx == 0:
                panel_title = f"{name}\n{panel_title}"

            ax.set_title(
                panel_title,
                fontsize=fontsize,
            )

            if ds_idx == 0:

                row_title = f"{var_name} {season}"

                if selected_depth is not None:
                    row_title = (
                        f"{var_name} {season} "
                        f"lev: {selected_depth} "
                    )

                if _has_depth:
                    row_title += " depth average"

                if "year" in obs_aligned.dims:
                    y0 = obs_aligned.year.values[0]
                    y1 = obs_aligned.year.values[-1]

                    row_title += f" ({y0}–{y1})"                    

                ax.text(
                    -0.25,
                    1.1,
                    row_title,
                    fontsize=fontsize,
                    transform=ax.transAxes,
                )

    divider = make_axes_locatable(ax)

    ax_cb = divider.append_axes(
        "bottom",
        size="10%",
        pad=0.1,
        axes_class=plt.Axes,
    )

    cbar = plt.colorbar(
        cb,
        cax=ax_cb,
        orientation="horizontal",
    )

    cbar.set_label(
        label=label,
        size=20,
    )

    cbar.ax.tick_params(labelsize=20)

    plt.tight_layout()
    plt.subplots_adjust(
        wspace=0.1,
        hspace=0.3,
    )

    if save:
        Path(dir_name).mkdir(
            parents=True,
            exist_ok=True,
        )

        plt.savefig(
            f"{dir_name}/{file_name}.png"
        )




def plot_global_map(ds_list: list[Exp],
                    data_dict: dict[Exp, state_dict],
                    data_sig_dict: dict[Exp, state_dict] | None = None,
                    central_longitude=180,
                    gridlines=False,
                    cmap=mpl.cm.RdYlBu,
                    vmin=-1,
                    vmax=1,
                    vals=None,
                    nvals=10,
                    cbar=False,
                    cbar_label='',
                    ticks_rotation=0,
                    title=None,
                    show_mean = True,
                    show_equator = False,
                    fnt_size=12,
                    figsize=None,
                    fig_dir=None,
                    fig_name=None,
                    save=False,
                    **kwargs): 
    """
    Plot global spatial maps for one or more datasets.

    A separate map is created for each dataset listed in ``ds_list``. The
    plotted field is displayed on a Plate Carree projection with configurable
    central longitude, color levels, map annotations, and optional statistical
    significance hatching.

    Parameters
    ----------
    ds_list : list[Exp]
        Dataset or experiment names to plot. Each entry must correspond to a
        key in ``data_dict``.
    data_dict : dict[Exp, state_dict]
        Mapping from dataset or experiment names to state objects containing
        the spatial data and associated ``y0`` and ``y1`` year metadata.
    data_sig_dict : dict[Exp, state_dict] or None, optional
        Mapping containing statistical significance fields corresponding to
        the datasets in ``ds_list``. If provided, significant regions are
        overlaid using hatching.
    central_longitude : float, optional
        Central longitude of the Plate Carree projection.
    gridlines : bool, optional
        If ``True``, draw map gridlines.
    cmap : str or matplotlib colormap, optional
        Colormap used to display the spatial field.
    vmin : float, optional
        Minimum value of the color scale.
    vmax : float, optional
        Maximum value of the color scale.
    vals : array-like or None, optional
        Explicit boundaries for the discrete color scale. If ``None``,
        equally spaced boundaries between ``vmin`` and ``vmax`` are
        generated.
    nvals : int, optional
        Number of color intervals used when ``vals`` is not provided.
    cbar : bool, optional
        If ``True``, add a horizontal colorbar to each figure.
    cbar_label : str, optional
        Label applied to the colorbar.
    ticks_rotation : float, optional
        Rotation angle, in degrees, for colorbar tick labels.
    title : str, optional
        Base figure title. The dataset name and corresponding year range are
        appended automatically.
    show_mean : bool, optional
        If ``True``, append the area-weighted spatial mean of the plotted
        field to the title.
    show_equator : bool, optional
        If ``True``, draw a dotted line along the equator.
    fnt_size : int, optional
        Base font size used for figure text.
    figsize : tuple or None, optional
        Figure size passed to ``matplotlib.pyplot.subplots``.
    fig_dir : str or Path or None, optional
        Directory in which to save figures when ``save`` is ``True``.
    fig_name : str or None, optional
        Filename used when saving the figure.
    save : bool, optional
        If ``True``, save the generated figure to ``fig_dir``.
    **kwargs
        Additional keyword arguments accepted by the function. These are
        currently not used internally.

    Notes
    -----
    A cyclic longitude point is added before plotting to avoid discontinuities
    at the map boundary. When ``central_longitude == 0``, longitude coordinates
    are shifted to the ``[-180, 180)`` range and sorted to remove the visual
    seam at the central longitude.

    If ``data_sig_dict`` is provided, its fields are overlaid with dotted
    hatching to indicate statistical significance.

    Documentation produced with the assistance of AI.
    """
        
    
    mpl.rcParams.update({'font.size': fnt_size})

    for name in ds_list:
        
        ds = add_cyclic_point(data_dict[name].data)
        ds_sig = add_cyclic_point(data_sig_dict[name].data) if data_sig_dict is not None else None
        y0 = data_dict[name].y0
        y1 = data_dict[name].y1

        if central_longitude == 0: # remove white line at central longitude
            ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
            ds = ds.sortby(ds.lon)
            if ds_sig is not None:
                ds_sig.coords['lon'] = (ds_sig.coords['lon'] + 180) % 360 - 180
                ds_sig = ds_sig.sortby(ds_sig.lon)

                
        lat = ds.lat
        lon = ds.lon
        img_extent = [lon[0], lon[-1], lat[0], lat[-1]]
            
        crs = ccrs.PlateCarree(central_longitude=central_longitude)    

        fig, ax = plt.subplots(nrows=1,
                            ncols=1, 
                            figsize=figsize, 
                            subplot_kw={'projection' : crs})                           
        
        if vals is not None:
            nvals = len(vals) - 1
            
        if vals is None:
            scale = (vmax-vmin)/float(nvals)
            vals = vmin + (vmax-vmin)*np.arange(nvals+1)/float(nvals)
        
        
        norm = mpl.colors.BoundaryNorm(vals, plt.cm.get_cmap(cmap).N)
        
        axis = ax
        if gridlines:
            axis.gridlines(draw_labels=False)

        im = axis.imshow(ds, 
                        origin='lower',
                        extent=img_extent,
                        # cmap=cmap,
                        cmap=plt.cm.get_cmap(cmap),
                        norm=norm,
                        interpolation='none',                     
                        transform=ccrs.PlateCarree())

        title_toplot = title + f" {name} {y0}-{y1}"
        if show_mean:
            title_toplot = title_toplot + f' ({np.round(area_weighted_avg(ds).values,2)})'
        im.set_clim(vmin,
                    vmax)
        
                
        if show_equator:
                    axis.plot(ds.lon,  # Longitude range
                            [0] * len(ds.lon),  # Latitude at the equator
                            color='black',  # Choose any color
                            linewidth=1, 
                            linestyle='dotted',  # Dashed line
                            transform=ccrs.PlateCarree())
        axis.coastlines()
        axis.set_title(title_toplot,
                    fontsize=fnt_size)

        if ds_sig is not None:  # statistical significance
            cs = axis.contourf(ds_sig,
                            1,
                            # hatches=['','....'],
                            hatches=['....'],
                            alpha=0,
                            # origin='lower',
                            extent=img_extent,
                            transform=ccrs.PlateCarree())
            
        if cbar:
            clb_x = 0.055 #0.095 
            clb_y = 0.05
            clb_w = 0.9 #0.8
            clb_h = 0.04

            cax = plt.axes([clb_x, # left
                            clb_y, # bottom
                            clb_w, # width
                            clb_h])# height
            cb = mpl.colorbar.ColorbarBase(ax=cax,
                                        cmap=plt.cm.get_cmap(cmap),
                                        # cmap=cmap,
                                        norm=norm,
                                        spacing='uniform',
                                        orientation='horizontal',
                                        extend='both',
                                        ticks=vals)

            cax.tick_params(labelsize=fnt_size-2)
            cb.set_ticks(ticks=vals, 
                        rotation=ticks_rotation,
                        labels=np.round(vals,3))
            
            cb.set_label(label=cbar_label,
                        size=fnt_size-2) 
            
        
        fig.tight_layout()
        if save:
            Path(fig_dir).mkdir(parents=True, exist_ok=True)
            plt.savefig(f'{fig_dir}/{fig_name}',
                        bbox_inches='tight',
                        dpi=300)
            



def plot_depth_vs_time_biomeavg(
    ds_list: list[Exp],
    ds_dicts: dict[Biome, dict[Exp, state_dict]],
    biome: Biome,
    ldyr: int = 0,
    title: str = "",
    figsize: tuple[float, float] = (45, 10),
    contourf_levels=None,
    cmap: str = "viridis",
    dir_name=None,
    file_name=None,
    colorbar_label: str | None = None,
    season: str = "ANN",
    lev_interp: np.ndarray | None = None,
    lev_range: float | tuple[float, float] | None = None,
    monthly_res: bool = False,
    ELNINO_years: NDArray[np.floating] | None = None,
    ELNINO_color: str = "r", 
    LANINA_years: NDArray[np.floating] | None = None,
    LANINA_color: str = "b", 
    return_fig_handles: bool = False,
    font_size: int = 14,
    save: bool = False,
):

    """
    Plot biome-averaged depth-versus-time sections for multiple datasets.

    The function creates one depth-time contour plot for each dataset in
    ``ds_list`` using data from the selected ``biome``. Depending on
    ``monthly_res``, the horizontal axis represents either seasonal/annual
    values indexed by year or monthly values represented on a continuous
    fractional-year axis.

    For each dataset, the function optionally interpolates the vertical
    coordinate to a common set of depth levels and restricts the plotted
    domain to a specified depth range. Filled contours are accompanied by
    labeled white contour lines. El Niño and La Niña years may optionally be
    indicated by vertical dotted lines.

    Parameters
    ----------
    ds_list : list[Exp]
        Dataset or experiment names to plot. Each entry must correspond to a
        key in the dictionary associated with ``biome`` in ``ds_dicts``.
        A separate subplot is created for each entry, in the order provided.
    ds_dicts : dict[Biome, dict[Exp, state_dict]]
        Nested mapping from biome names to dataset or experiment names and
        their associated state objects. Each selected state object must
        contain the data to plot through its ``data`` attribute.
    biome : Biome
        Biome whose horizontally averaged depth-time fields are plotted.
        The value must correspond to a key in ``ds_dicts``.
    ldyr : int, optional
        Lead year used to select the temporal window. In seasonal mode,
        ``seasonal_mean`` is called with ``ldyr_ini=ldyr`` and
        ``ldyr_end=ldyr + 1``. In monthly-resolution mode, the corresponding
        monthly indices are obtained with ``get_season_indices``.
    title : str, optional
        Base title prepended to each subplot title. The dataset name, biome,
        and selected season are appended automatically.
    figsize : tuple[float, float], optional
        Figure size passed to ``matplotlib.pyplot.subplots``.
    contourf_levels : array-like or int or None, optional
        Levels supplied to ``matplotlib.axes.Axes.contourf``. If ``None``,
        Matplotlib determines the contour levels automatically. The same
        resulting levels are also used for the overlaid contour lines.
    cmap : str, optional
        Matplotlib colormap used for the filled contours.
    dir_name : str or Path or None, optional
        Directory in which to save the figure when ``save`` is ``True``.
        The directory is created if it does not already exist.
    file_name : str or None, optional
        Output filename without the ``.png`` extension. Required when
        ``save=True``.
    colorbar_label : str or None, optional
        Label applied to the shared vertical colorbar. If ``None``, no
        colorbar label is added.
    season : str, optional
        Season passed to ``seasonal_mean`` when ``monthly_res=False``.
        When ``monthly_res=True``, this parameter must be ``"ANN"`` because
        the function retains the individual monthly values rather than
        calculating a seasonal average.
    lev_interp : numpy.ndarray or None, optional
        Target depth levels used to interpolate each dataset along the
        ``lev`` dimension. If ``None``, the native depth coordinate of each
        dataset is retained.
    lev_range : float or tuple[float, float] or None, optional
        Depth interval retained for plotting. If a scalar is provided, it is
        interpreted as the maximum depth and the retained interval is
        ``0 <= lev <= lev_range``. If a two-element tuple is provided, it is
        interpreted as ``(minimum_depth, maximum_depth)``. If ``None``, the
        full available depth range is used.
    monthly_res : bool, optional
        If ``False``, calculate the requested seasonal mean and plot the
        result against year. If ``True``, retain monthly values for the
        selected lead year and construct a continuous fractional-year
        coordinate from the ``year`` and ``month`` dimensions.
    ELNINO_years : NDArray[np.floating] or None, optional
        El Niño event years or fractional-year locations to mark with
        vertical dotted lines. Markers outside the plotted temporal range
        are ignored.
    ELNINO_color : str, optional
        Matplotlib color specification used for El Niño markers.
    LANINA_years : NDArray[np.floating] or None, optional
        La Niña event years or fractional-year locations to mark with
        vertical dotted lines. Markers outside the plotted temporal range
        are ignored.
    LANINA_color : str, optional
        Matplotlib color specification used for La Niña markers.
    return_fig_handles : bool, optional
        If ``True``, return the Matplotlib figure and axes objects. If
        ``False``, the function returns nothing.
    font_size : int, optional
        Font size used for subplot titles, axis labels, tick labels, contour
        labels, and the colorbar.
    save : bool, optional
        If ``True``, save the figure as a PNG file in ``dir_name`` using
        ``file_name``.

    Returns
    -------
    tuple[matplotlib.figure.Figure, numpy.ndarray] or None
        If ``return_fig_handles=True``, returns ``(fig, axes)``, where ``fig``
        is the Matplotlib figure and ``axes`` is the two-dimensional array of
        subplot axes created by ``plt.subplots``. Otherwise, returns ``None``.

    Raises
    ------
    ValueError
        If ``monthly_res=True`` and ``season`` is not ``"ANN"``.
    ValueError
        If ``ds_list`` is empty.
    ValueError
        If a selected dataset does not contain a ``lev`` dimension.
    ValueError
        If ``save=True`` and either ``dir_name`` or ``file_name`` is not
        provided.

    Notes
    -----
    In seasonal mode, each dataset is first reduced with ``seasonal_mean``
    over the requested season and lead year. The resulting horizontal
    coordinate is expected to be ``year``.

    In monthly-resolution mode, the function selects all monthly indices
    belonging to the requested lead year, stacks the ``year`` and ``month``
    dimensions into a ``yearmonth`` dimension, and constructs a continuous
    fractional-year coordinate as

    ``year + (mod(month, 12) + 0.5) / 12``.

    This places each monthly value approximately at the midpoint of its
    corresponding month on the time axis.

    Depth interpolation, when requested, is performed before the depth-range
    restriction. The depth axis is inverted so that shallower values appear
    near the top of each subplot and greater depths appear lower on the
    figure.

    All datasets share a single colorbar based on the final ``contourf``
    object created in the plotting loop. For most meaningful visual
    comparisons, a common ``contourf_levels`` specification should therefore
    be used across datasets.

    El Niño and La Niña markers are drawn only when their supplied values lie
    within the minimum and maximum values of the plotted time coordinate.

    Documentation produced with the assistance of AI.
    """


    if monthly_res and season != "ANN":
        raise ValueError(
            "monthly_res=True currently requires season='ANN'."
        )

    if not ds_list:
        raise ValueError("'ds_list' cannot be empty.")

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(
        len(ds_list),
        1,
        figsize=figsize,
        squeeze=False,
    )

    ds_dict = ds_dicts[biome]

    contour_f = None

    # ------------------------------------------------------------------
    # Datasets
    # ------------------------------------------------------------------
    for ds_idx, name in enumerate(ds_list):

        ax = axes[ds_idx, 0]

        data = ds_dict[name].data

        # --------------------------------------------------------------
        # Seasonal / monthly selection
        # --------------------------------------------------------------
        if monthly_res:
            month_indices = get_season_indices(
                season="ANN",
                ldyr_ini=ldyr,
                ldyr_end=ldyr + 1,
            )
            
            ts = data.isel(
                month=month_indices,
            )

        else:
            ts = seasonal_mean(
                data,
                season=season,
                ldyr_ini=ldyr,
                ldyr_end=ldyr + 1,
            )

        # --------------------------------------------------------------
        # Depth handling
        # --------------------------------------------------------------
        if "lev" not in ts.dims:
            raise ValueError(
                f"Dataset {name!r} does not contain a 'lev' dimension."
            )

        if lev_interp is not None:
            ts = ts.interp(
                lev=lev_interp,
            )

        if lev_range is not None:

            if np.isscalar(lev_range):
                lev_min = 0
                lev_max = lev_range

            else:
                lev_min, lev_max = lev_range

            ts = ts.where(
                (ts.lev >= lev_min)
                & (ts.lev <= lev_max),
                drop=True,
            )

        # --------------------------------------------------------------
        # Horizontal coordinate
        # --------------------------------------------------------------
        if monthly_res:
            ts = ts.stack(
                yearmonth=("year", "month"),
            )

            yearmonth = (
                ts.year.values
                + (np.mod(ts.month.values, 12) + 0.5) / 12
            )

            ts = ts.assign_coords(
                yearmonth=yearmonth,
            )

            dim = "yearmonth"

        else:
            dim = "year"

        xx = ts[dim].values

        # --------------------------------------------------------------
        # Plot
        # --------------------------------------------------------------
        contour_f = ax.contourf(
            xx,
            ts.lev.values,
            ts,
            levels=contourf_levels,
            cmap=cmap,
        )

        contours = ax.contour(
            xx,
            ts.lev.values,
            ts,
            levels=contour_f.levels,
            colors="white",
        )

        ax.clabel(
            contours,
            inline=True,
            fontsize=font_size,
            colors="white",
        )

        # --------------------------------------------------------------
        # Axis formatting
        # --------------------------------------------------------------
        years = np.unique(
            np.floor(xx)
        )

        ax.set_xticks(
            years,
            labels=years.astype(int),
            rotation=45,
        )

        ax.set_title(
            f"{title} - {name} - {biome} - {season}",
            fontsize=font_size,
        )

        ax.set_ylabel(
            "depth (m)",
            fontsize=font_size,
        )

        ax.tick_params(
            axis="both",
            labelsize=font_size,
        )

        ax.invert_yaxis()

        # --------------------------------------------------------------
        # ENSO markers
        # --------------------------------------------------------------
        if ELNINO_years is not None:
            for year in ELNINO_years:
                if xx.min() <= year <= xx.max():
                    ax.axvline(
                        x=year,
                        linestyle="dotted",
                        color=ELNINO_color,
                        alpha=0.5,
                    )

        if LANINA_years is not None:
            for year in LANINA_years:
                if xx.min() <= year <= xx.max():
                    ax.axvline(
                        x=year,
                        linestyle="dotted",
                        color=LANINA_color,
                        alpha=0.5,
                    )

    # ------------------------------------------------------------------
    # Shared colorbar
    # ------------------------------------------------------------------
    cbar = fig.colorbar(
        contour_f,
        ax=axes[:, 0],
        orientation="vertical",
    )

    if colorbar_label is not None:
        cbar.set_label(
            colorbar_label,
            fontsize=font_size,
        )

    cbar.ax.tick_params(
        labelsize=font_size,
    )

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    if save:
        if dir_name is None or file_name is None:
            raise ValueError(
                "'dir_name' and 'file_name' are required when save=True."
            )

        Path(dir_name).mkdir(
            parents=True,
            exist_ok=True,
        )

        plt.savefig(
            f"{dir_name}/{file_name}.png",
            bbox_inches="tight",
            dpi=300,
        )

    if return_fig_handles:
        return fig, axes


def snapshot_depth_vs_lon(
    ds_list: list[Exp],
    ds_dicts: dict[Biome, dict[Exp, state_dict]],
    biome: Biome,
    years_to_plot: int | list[int] | np.ndarray | tuple[int, int] | None = None,
    contour_dict: dict[Biome, dict[Exp, state_dict]] | None = None,
    quiver_dicts: dict[Biome, dict[Exp, state_dict]] | None = None,
    quiver_axis: Literal["U", "V"] | None = None,
    ldyr: int = 0,
    title: str = "",
    figsize: tuple[float, float] = (10, 45),
    contourf_levels=None,
    contour_var: str | None = None,
    contour_var_levels=None,
    cmap: str = "viridis",
    dir_name=None,
    file_name=None,
    colorbar_label: str | None = None,
    season: str = "ANN",
    lev_interp: np.ndarray | None = None,
    lev_range: float | tuple[float, float] | None = None,
    return_fig_handles: bool = False,
    headwidth: float = 5,
    headlength: float = 1,
    headaxislength: float = 2,
    width: float = 0.005,
    font_size: int = 14,
    save: bool = False,
):
    

    """
    Plot longitude-depth snapshots or composites for multiple datasets.

    The function creates one longitude-versus-depth cross-section for each
    dataset in ``ds_list`` using data associated with the selected ``biome``.
    The primary field is displayed using filled contours and may optionally be
    supplemented with contours from a second variable and directional quiver
    vectors from an additional field.

    Temporal selection and compositing are handled by ``prepare_snapshot``,
    allowing the plot to represent a specific year, multiple selected years,
    a range of years, or the available climatological period, depending on
    ``years_to_plot``. The resulting data may also be interpolated to common
    depth levels and restricted to a requested depth interval before plotting.

    Parameters
    ----------
    ds_list : list[Exp]
        Dataset or experiment names to plot. Each entry must correspond to a
        key in ``ds_dicts[biome]``. A separate subplot is created for each
        dataset, preserving the order of ``ds_list``.
    ds_dicts : dict[Biome, dict[Exp, state_dict]]
        Nested mapping from biome names to dataset or experiment names and
        their associated state objects. The ``data`` attribute of each state
        object supplies the primary longitude-depth field.
    biome : Biome
        Biome whose longitude-depth cross-sections are plotted. The value must
        correspond to a key in ``ds_dicts`` and, when provided, in
        ``contour_dict`` and ``quiver_dicts``.
    years_to_plot : int, list[int], numpy.ndarray, tuple[int, int], or None, optional
        Year selection passed to ``prepare_snapshot``. The interpretation of
        this argument is handled by that helper function and may represent a
        single year, multiple selected years, or a year range. If ``None``,
        the helper uses its default temporal averaging behavior.
    contour_dict : dict[Biome, dict[Exp, state_dict]] or None, optional
        Optional nested mapping containing a secondary variable to overlay as
        contour lines. The biome and experiment structure is expected to
        correspond to that of ``ds_dicts``. If provided,
        ``contour_var_levels`` must also be specified.
    quiver_dicts : dict[Biome, dict[Exp, state_dict]] or None, optional
        Optional nested mapping containing a field to represent using quiver
        arrows. The quiver field is processed with the same temporal and depth
        selections as the primary field.
    quiver_axis : {"U", "V"} or None, optional
        Direction assigned to values from ``quiver_dicts``. If ``"U"``, the
        supplied field is treated as the horizontal quiver component and the
        vertical component is set to zero. If ``"V"``, the supplied field is
        treated as the vertical component and the horizontal component is set
        to zero. Must be specified as ``"U"`` or ``"V"`` when
        ``quiver_dicts`` is provided.
    ldyr : int, optional
        Lead year passed to ``prepare_snapshot`` for temporal selection and
        seasonal compositing.
    title : str, optional
        Base text included at the beginning of each subplot title. Dataset
        name, biome, temporal description, season, and optional contour
        variable information are appended automatically.
    figsize : tuple[float, float], optional
        Figure size passed to ``matplotlib.pyplot.subplots``.
    contourf_levels : array-like, int, or None, optional
        Levels supplied to ``matplotlib.axes.Axes.contourf`` for the primary
        field. If ``None``, Matplotlib determines the filled-contour levels
        automatically.
    contour_var : str or None, optional
        Descriptive name of the secondary contour variable. When both
        ``contour_dict`` and ``contour_var`` are provided, this name is added
        to the subplot title.
    contour_var_levels : array-like or int or None, optional
        Contour levels used for the secondary field supplied through
        ``contour_dict``. This argument is required when ``contour_dict`` is
        provided.
    cmap : str, optional
        Matplotlib colormap used for the filled contours of the primary field.
    dir_name : str or Path or None, optional
        Directory in which to save the figure when ``save=True``. The
        directory is created if it does not already exist.
    file_name : str or None, optional
        Output filename without the ``.png`` extension. Required when
        ``save=True``.
    colorbar_label : str or None, optional
        Label applied to the shared vertical colorbar. If ``None``, no label
        is added.
    season : str, optional
        Season passed to ``prepare_snapshot`` for construction of the plotted
        snapshot or composite.
    lev_interp : numpy.ndarray or None, optional
        Target depth coordinates used by ``select_depth_range`` to interpolate
        the primary, contour, and quiver fields before plotting. If ``None``,
        each field retains its native vertical coordinates.
    lev_range : float or tuple[float, float] or None, optional
        Depth interval passed to ``select_depth_range``. A scalar typically
        represents a maximum depth, while a two-element tuple specifies an
        explicit minimum and maximum depth. If ``None``, the available depth
        range is retained.
    return_fig_handles : bool, optional
        If ``True``, return the Matplotlib figure and axes objects. If
        ``False``, the function returns nothing.
    headwidth : float, optional
        Width of quiver arrow heads passed to
        ``matplotlib.axes.Axes.quiver``.
    headlength : float, optional
        Length of quiver arrow heads passed to
        ``matplotlib.axes.Axes.quiver``.
    headaxislength : float, optional
        Length of the quiver arrow head along its axis, passed to
        ``matplotlib.axes.Axes.quiver``.
    width : float, optional
        Shaft width of the quiver arrows.
    font_size : int, optional
        Font size used for subplot titles, axis labels, tick labels, contour
        labels, and colorbar text.
    save : bool, optional
        If ``True``, save the generated figure as a PNG file.

    Returns
    -------
    tuple[matplotlib.figure.Figure, numpy.ndarray] or None
        If ``return_fig_handles=True``, returns ``(fig, axes)``, where ``fig``
        is the Matplotlib figure and ``axes`` is the two-dimensional array of
        subplot axes returned by ``plt.subplots``. Otherwise, returns
        ``None``.

    Raises
    ------
    ValueError
        If ``ds_list`` is empty.
    ValueError
        If ``quiver_dicts`` is provided and ``quiver_axis`` is neither
        ``"U"`` nor ``"V"``.
    ValueError
        If ``contour_dict`` is provided without ``contour_var_levels``.
    ValueError
        If ``save=True`` and either ``dir_name`` or ``file_name`` is not
        provided.

    Notes
    -----
    The primary field, optional contour field, and optional quiver field are
    independently passed through ``prepare_snapshot`` using the same
    ``season``, ``ldyr``, and ``years_to_plot`` arguments. They are then
    passed through ``select_depth_range`` using the same interpolation and
    depth-selection settings.

    When a secondary contour field is supplied, it is aligned with the
    primary field using ``xr.align(..., join="inner")`` before plotting.
    Quiver data are aligned with the primary field in the same way. This
    restricts the corresponding fields to common coordinates before the
    overlays are drawn.

    If ``contour_dict`` is provided, the secondary field is shown as
    semi-transparent black contour lines at ``contour_var_levels``. Otherwise,
    contour lines are generated directly from the primary field using the
    same levels produced by ``contourf`` and are displayed in white.

    Quiver values represent only one vector component. For
    ``quiver_axis="U"``, the supplied values define the longitude-direction
    component and the depth-direction component is zero. For
    ``quiver_axis="V"``, the supplied values define the depth-direction
    component and the longitude-direction component is zero. Quiver arrows
    are subsampled every second longitude and depth grid point, and NaN values
    are replaced by zero before plotting.

    The plotted depth axis is inverted so that shallower depths appear near
    the top and greater depths appear toward the bottom of each panel.

    The temporal label included in each subplot title is obtained from
    ``prepare_snapshot`` when available. If no year string is returned, the
    function falls back to the ``y0`` and ``y1`` metadata stored in the
    corresponding primary ``state_dict`` and labels the period as a mean.

    A single vertical colorbar is shared across all subplots and is based on
    the final filled-contour object created in the dataset loop. For direct
    comparison among experiments, supplying common ``contourf_levels`` across
    all datasets is therefore recommended.

    Documentation produced with the assistance of AI.
    """

    if not ds_list:
        raise ValueError(
            "'ds_list' cannot be empty."
        )

    if quiver_dicts is not None and quiver_axis not in ("U", "V"):
        raise ValueError(
            "'quiver_axis' must be either 'U' or 'V' "
            "when quiver data are provided."
        )

    if contour_dict is not None and contour_var_levels is None:
        raise ValueError(
            "'contour_var_levels' must be provided when "
            "'contour_dict' is specified."
        )

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(
        len(ds_list),
        1,
        figsize=figsize,
        squeeze=False,
    )

    ds_dict = ds_dicts[biome]

    contour_f = None

    # ------------------------------------------------------------------
    # Datasets
    # ------------------------------------------------------------------
    for ds_idx, name in enumerate(ds_list):

        ax = axes[ds_idx, 0]

        # --------------------------------------------------------------
        # Primary field
        # --------------------------------------------------------------
        ts, year_string = prepare_snapshot(
            ds_dict[name].data,
            season=season,
            ldyr=ldyr,
            years_to_plot=years_to_plot,
            return_year_string=True,
        )

        ts = select_depth_range(
            ts,
            lev_interp=lev_interp,
            lev_range=lev_range,
        )

        # --------------------------------------------------------------
        # Contour field
        # --------------------------------------------------------------
        contour_ts = None

        if contour_dict is not None:
            contour_ts = prepare_snapshot(
                contour_dict[biome][name].data,
                season=season,
                ldyr=ldyr,
                years_to_plot=years_to_plot,
            )

            contour_ts = select_depth_range(
                contour_ts,
                lev_interp=lev_interp,
                lev_range=lev_range,
            )

        # --------------------------------------------------------------
        # Quiver field
        # --------------------------------------------------------------
        qv_ts = None

        if quiver_dicts is not None:
            qv_ts = prepare_snapshot(
                quiver_dicts[biome][name].data,
                season=season,
                ldyr=ldyr,
                years_to_plot=years_to_plot,
            )

            qv_ts = select_depth_range(
                qv_ts,
                lev_interp=lev_interp,
                lev_range=lev_range,
            )

        # --------------------------------------------------------------
        # Align longitude / depth coordinates
        # --------------------------------------------------------------
        if contour_ts is not None:
            ts, contour_ts = xr.align(
                ts,
                contour_ts,
                join="inner",
            )

        if qv_ts is not None:
            ts, qv_ts = xr.align(
                ts,
                qv_ts,
                join="inner",
            )

        xx = ts.lon.values

        # --------------------------------------------------------------
        # Filled contours
        # --------------------------------------------------------------
        contour_f = ax.contourf(
            xx,
            ts.lev.values,
            ts,
            levels=contourf_levels,
            cmap=cmap,
        )

        # --------------------------------------------------------------
        # Contours
        # --------------------------------------------------------------
        if contour_ts is not None:
            contours = ax.contour(
                contour_ts.lon.values,
                contour_ts.lev.values,
                contour_ts,
                colors="black",
                levels=contour_var_levels,
                alpha=0.5,
            )

            ax.clabel(
                contours,
                inline=True,
                fontsize=font_size,
                colors="black",
            )

        else:
            contours = ax.contour(
                xx,
                ts.lev.values,
                ts,
                colors="white",
                levels=contour_f.levels,
            )

            ax.clabel(
                contours,
                inline=True,
                fontsize=font_size,
                colors="white",
            )

        # --------------------------------------------------------------
        # Quiver
        # --------------------------------------------------------------
        if qv_ts is not None:

            if quiver_axis == "U":
                U = qv_ts.values
                V = np.zeros_like(U)

            else:
                V = qv_ts.values
                U = np.zeros_like(V)

            ax.quiver(
                qv_ts.lon.values[::2],
                qv_ts.lev.values[::2],
                np.nan_to_num(
                    U[::2, ::2],
                    nan=0.0,
                ),
                np.nan_to_num(
                    V[::2, ::2],
                    nan=0.0,
                ),
                alpha=0.5,
                width=width,
                headwidth=headwidth,
                headlength=headlength,
                headaxislength=headaxislength,
            )

        # --------------------------------------------------------------
        # Year fallback
        # --------------------------------------------------------------
        if year_string is None:
            y0 = ds_dict[name].y0
            y1 = ds_dict[name].y1

            if y0 is not None and y1 is not None:
                year_string = f"{y0} - {y1} mean"
            else:
                year_string = ""

        # --------------------------------------------------------------
        # Title
        # --------------------------------------------------------------
        title_parts = [
            title,
            name,
            biome,
            year_string,
            f"{season} composite",
        ]

        if contour_dict is not None and contour_var is not None:
            title_parts.append(
                f"{contour_var} contours"
            )

        title_ = " - ".join(
            part for part in title_parts if part
        )

        ax.set_title(
            title_,
            fontsize=font_size,
        )

        # --------------------------------------------------------------
        # Axis formatting
        # --------------------------------------------------------------
        ax.set_ylabel(
            "depth (m)",
            fontsize=font_size,
        )

        if ds_idx == len(ds_list) - 1:
            ax.set_xlabel(
                "Lon ($^o$ East)",
                fontsize=font_size,
            )
        else:
            ax.set_xlabel("")

        ax.invert_yaxis()

        ax.tick_params(
            axis="both",
            which="major",
            labelsize=font_size,
        )

    # ------------------------------------------------------------------
    # Shared colorbar
    # ------------------------------------------------------------------
    cbar = fig.colorbar(
        contour_f,
        ax=axes[:, 0],
        orientation="vertical",
    )

    if colorbar_label is not None:
        cbar.set_label(
            colorbar_label,
            fontsize=font_size,
        )

    cbar.ax.tick_params(
        labelsize=font_size,
    )

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    if save:
        if dir_name is None or file_name is None:
            raise ValueError(
                "'dir_name' and 'file_name' are required "
                "when save=True."
            )

        Path(dir_name).mkdir(
            parents=True,
            exist_ok=True,
        )

        plt.savefig(
            f"{dir_name}/{file_name}.png",
            bbox_inches="tight",
            dpi=300,
        )

    if return_fig_handles:
        return fig, axes




def snapshot_aerial(
    ds_list: list[Exp],
    ds_dicts: dict[Biome, dict[Exp, state_dict]],
    biome: Biome,
    years_to_plot: int | list[int] | np.ndarray | tuple[int, int] | None = None,
    quiver_dicts: dict[Biome, dict[Exp, state_dict]] | None = None,
    quiver_axis: Literal["U", "V"] | None = None,
    ldyr: int = 0,
    title: str = "",
    figsize: tuple[float, float] = (10, 45),
    vmax: float | None = None,
    vmin: float | None = None,
    cmap: str = "viridis",
    dir_name=None,
    file_name=None,
    colorbar_label: str | None = None,
    season: str = "ANN",
    lev_interp: np.ndarray | None = None,
    lev_range: float | tuple[float, float] | None = None,
    integrate: bool = False,
    return_fig_handles: bool = False,
    width: float = 0.005,
    headwidth: float = 5,
    headlength: float = 1,
    headaxislength: float = 2,
    font_size: int = 14,
    save: bool = False,
):

    """
    Plot aeraial spatial snapshots or composites for multiple datasets.

    The function creates one latitude-longitude panel for each dataset in
    ``ds_list`` using data associated with the selected ``biome``. Temporal
    selection and compositing are handled by ``prepare_snapshot``. If the
    resulting field contains a depth dimension, the selected vertical range can
    optionally be interpolated and then either averaged or integrated before
    plotting.

    An optional secondary field can be overlaid as quiver arrows. The supplied
    quiver variable is interpreted as either a zonal-like or meridional-like
    component according to ``quiver_axis``, while the orthogonal component is set
    to zero.

    Parameters
    ----------
    ds_list : list[Exp]
        Dataset or experiment names to plot. Each entry must correspond to a key
        in ``ds_dicts[biome]``. A separate subplot is created for each dataset in
        the order provided.
    ds_dicts : dict[Biome, dict[Exp, state_dict]]
        Nested mapping from biome names to dataset or experiment names and their
        associated state objects. The ``data`` attribute of each selected state
        object provides the primary field.
    biome : Biome
        Biome whose horizontal fields are plotted. The value must correspond to a
        key in ``ds_dicts`` and, when provided, in ``quiver_dicts``.
    years_to_plot : int, list[int], numpy.ndarray, tuple[int, int], or None, optional
        Year selection passed to ``prepare_snapshot``. Depending on the helper
        implementation, this may represent a single year, a collection of years,
        or a year range. If ``None``, the default temporal averaging behavior of
        ``prepare_snapshot`` is used.
    quiver_dicts : dict[Biome, dict[Exp, state_dict]] or None, optional
        Optional nested mapping containing a secondary field to display using
        quiver arrows. The field is processed using the same temporal and depth
        selections as the primary field.
    quiver_axis : {"U", "V"} or None, optional
        Direction assigned to values from ``quiver_dicts``. If ``"U"``, the
        supplied field is treated as the horizontal U component and the V
        component is set to zero. If ``"V"``, the supplied field is treated as
        the V component and the U component is set to zero. Must be specified as
        ``"U"`` or ``"V"`` when ``quiver_dicts`` is provided.
    ldyr : int, optional
        Lead year passed to ``prepare_snapshot`` for temporal selection and
        seasonal compositing.
    title : str, optional
        Base text prepended to each subplot title. Dataset name, biome, temporal
        description, season, selected depth range, and depth-reduction method are
        appended automatically where applicable.
    figsize : tuple[float, float], optional
        Figure size passed to ``matplotlib.pyplot.subplots``.
    vmax : float or None, optional
        Upper limit of the color scale used by ``pcolormesh``. If ``None``,
        Matplotlib determines the upper bound automatically.
    vmin : float or None, optional
        Lower limit of the color scale used by ``pcolormesh``. If ``None``,
        Matplotlib determines the lower bound automatically.
    cmap : str, optional
        Matplotlib colormap used for the primary spatial field.
    dir_name : str or Path or None, optional
        Directory in which to save the figure when ``save=True``. The directory
        is created if it does not already exist.
    file_name : str or None, optional
        Output filename without the ``.png`` extension. Required when
        ``save=True``.
    colorbar_label : str or None, optional
        Label applied to the shared vertical colorbar. If ``None``, no colorbar
        label is added.
    season : str, optional
        Season passed to ``prepare_snapshot`` when constructing the plotted
        snapshot or composite.
    lev_interp : numpy.ndarray or None, optional
        Target depth coordinates passed to ``select_depth_range`` for vertical
        interpolation. If ``None``, native depth coordinates are retained.
    lev_range : float or tuple[float, float] or None, optional
        Depth interval passed to ``select_depth_range``. A scalar typically
        specifies a maximum depth, while a two-element tuple specifies an
        explicit minimum and maximum depth. If ``None``, the available depth
        range is retained.
    integrate : bool, optional
        Method used to collapse the ``lev`` dimension when one is present. If
        ``True``, integrate over depth using ``xarray.DataArray.integrate``. If
        ``False``, calculate the arithmetic mean over depth.
    return_fig_handles : bool, optional
        If ``True``, return the Matplotlib figure and axes objects. If ``False``,
        the function returns nothing.
    width : float, optional
        Shaft width of quiver arrows.
    headwidth : float, optional
        Width of quiver arrow heads passed to
        ``matplotlib.axes.Axes.quiver``.
    headlength : float, optional
        Length of quiver arrow heads passed to
        ``matplotlib.axes.Axes.quiver``.
    headaxislength : float, optional
        Length of the quiver arrow head along its axis, passed to
        ``matplotlib.axes.Axes.quiver``.
    font_size : int, optional
        Font size used for subplot titles, axis labels, tick labels, and colorbar
        text.
    save : bool, optional
        If ``True``, save the generated figure as a PNG file.

    Returns
    -------
    tuple[matplotlib.figure.Figure, numpy.ndarray] or None
        If ``return_fig_handles=True``, returns ``(fig, axes)``, where ``fig`` is
        the Matplotlib figure and ``axes`` is the two-dimensional array of subplot
        axes returned by ``plt.subplots``. Otherwise, returns ``None``.

    Raises
    ------
    ValueError
        If ``ds_list`` is empty.
    ValueError
        If ``quiver_dicts`` is provided and ``quiver_axis`` is neither ``"U"``
        nor ``"V"``.
    ValueError
        If ``save=True`` and either ``dir_name`` or ``file_name`` is not
        provided.

    Notes
    -----
    The primary field is first processed with ``prepare_snapshot`` using the
    requested ``season``, ``ldyr``, and ``years_to_plot`` settings.

    If the resulting field contains a ``lev`` dimension, ``select_depth_range``
    is applied before the vertical dimension is reduced. When ``integrate=True``,
    the field is vertically integrated using ``integrate("lev")``; otherwise,
    the arithmetic mean over ``lev`` is calculated. If no depth dimension is
    present, no vertical processing is applied.

    The optional quiver field undergoes the same temporal selection and, when
    applicable, the same depth interpolation, depth restriction, and depth
    reduction as the primary field. The two fields are then aligned using
    ``xr.align(..., join="inner")`` so that only common coordinates are retained
    before plotting.

    Quiver data represent a single vector component. For ``quiver_axis="U"``,
    the supplied values define the U component while the V component is zero. For
    ``quiver_axis="V"``, the supplied values define the V component while the U
    component is zero. Quiver arrows are subsampled every second latitude and
    longitude grid point, and NaN values are replaced by zero before plotting.

    The temporal label included in each subplot title is obtained from
    ``prepare_snapshot`` when available. If no year string is returned, the
    function falls back to the ``y0`` and ``y1`` metadata stored in the
    corresponding primary ``state_dict`` and labels that period as a mean.

    When depth selection is applied, the range description returned by
    ``select_depth_range`` is included in the subplot title. The title also
    indicates whether the plotted field is a depth mean or a depth integral.

    A single vertical colorbar is shared across all subplots and is based on the
    final ``pcolormesh`` object created in the dataset loop. Supplying common
    ``vmin`` and ``vmax`` values is therefore recommended when direct comparison
    among experiments is desired.

    Documentation produced with the assistance of AI.
    """


    if not ds_list:
        raise ValueError(
            "'ds_list' cannot be empty."
        )

    if quiver_dicts is not None and quiver_axis not in ("U", "V"):
        raise ValueError(
            "'quiver_axis' must be either 'U' or 'V' "
            "when quiver data are provided."
        )

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(
        len(ds_list),
        1,
        figsize=figsize,
        squeeze=False,
    )

    ds_dict = ds_dicts[biome]

    im = None

    # ------------------------------------------------------------------
    # Datasets
    # ------------------------------------------------------------------
    for ds_idx, name in enumerate(ds_list):

        ax = axes[ds_idx, 0]

        # --------------------------------------------------------------
        # Primary field
        # --------------------------------------------------------------
        ts, year_string = prepare_snapshot(
            ds_dict[name].data,
            season=season,
            ldyr=ldyr,
            years_to_plot=years_to_plot,
            return_year_string=True,
        )

        has_lev = "lev" in ts.dims
        range_string = ""

        if has_lev:
            ts, range_string = select_depth_range(
                ts,
                lev_interp=lev_interp,
                lev_range=lev_range,
                return_range_string=True,
            )

            if integrate:
                ts = ts.integrate("lev")
            else:
                ts = ts.mean("lev")

        # --------------------------------------------------------------
        # Quiver field
        # --------------------------------------------------------------
        qv_ts = None

        if quiver_dicts is not None:
            qv_ts = prepare_snapshot(
                quiver_dicts[biome][name].data,
                season=season,
                ldyr=ldyr,
                years_to_plot=years_to_plot,
            )

            if "lev" in qv_ts.dims:
                qv_ts = select_depth_range(
                    qv_ts,
                    lev_interp=lev_interp,
                    lev_range=lev_range,
                )

                if integrate:
                    qv_ts = qv_ts.integrate("lev")
                else:
                    qv_ts = qv_ts.mean("lev")

        # --------------------------------------------------------------
        # Align coordinates
        # --------------------------------------------------------------
        if qv_ts is not None:
            ts, qv_ts = xr.align(
                ts,
                qv_ts,
                join="inner",
            )

        # --------------------------------------------------------------
        # Plot
        # --------------------------------------------------------------
        im = ax.pcolormesh(
            ts.lon,
            ts.lat,
            ts,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )

        # --------------------------------------------------------------
        # Quiver
        # --------------------------------------------------------------
        if qv_ts is not None:

            if quiver_axis == "U":
                U = qv_ts.values
                V = np.zeros_like(U)

            else:
                V = qv_ts.values
                U = np.zeros_like(V)

            ax.quiver(
                qv_ts.lon.values[::2],
                qv_ts.lat.values[::2],
                np.nan_to_num(
                    U[::2, ::2],
                    nan=0.0,
                ),
                np.nan_to_num(
                    V[::2, ::2],
                    nan=0.0,
                ),
                alpha=0.5,
                width=width,
                headwidth=headwidth,
                headlength=headlength,
                headaxislength=headaxislength,
            )

        # --------------------------------------------------------------
        # Year fallback
        # --------------------------------------------------------------
        if year_string is None:
            y0 = ds_dict[name].y0
            y1 = ds_dict[name].y1

            if y0 is not None and y1 is not None:
                year_string = f"{y0} - {y1} mean"
            else:
                year_string = ""

        # --------------------------------------------------------------
        # Title
        # --------------------------------------------------------------
        title_parts = [
            title,
            name,
            biome,
            year_string,
            f"{season} composite",
            range_string,
        ]

        if has_lev:
            if integrate:
                title_parts.append(
                    "depth integrated"
                )
            else:
                title_parts.append(
                    "depth mean"
                )

        title_ = " - ".join(
            part for part in title_parts if part
        )

        ax.set_title(
            title_,
            fontsize=font_size,
        )

        # --------------------------------------------------------------
        # Axis formatting
        # --------------------------------------------------------------
        ax.set_ylabel(
            "Lat ($^o$ North)",
            fontsize=font_size,
        )

        if ds_idx == len(ds_list) - 1:
            ax.set_xlabel(
                "Lon ($^o$ East)",
                fontsize=font_size,
            )
        else:
            ax.set_xlabel("")

        ax.tick_params(
            axis="both",
            which="major",
            labelsize=font_size,
        )

    # ------------------------------------------------------------------
    # Shared colorbar
    # ------------------------------------------------------------------
    cbar = fig.colorbar(
        im,
        ax=axes[:, 0],
        orientation="vertical",
    )

    if colorbar_label is not None:
        cbar.set_label(
            colorbar_label,
            fontsize=font_size,
        )

    cbar.ax.tick_params(
        labelsize=font_size,
    )

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    if save:
        if dir_name is None or file_name is None:
            raise ValueError(
                "'dir_name' and 'file_name' are required "
                "when save=True."
            )

        Path(dir_name).mkdir(
            parents=True,
            exist_ok=True,
        )

        plt.savefig(
            f"{dir_name}/{file_name}.png",
            bbox_inches="tight",
            dpi=300,
        )

    if return_fig_handles:
        return fig, axes
            

def snapshot_crossection(
    ds_list: list[Exp],
    ds_dicts: dict[Biome, dict[Exp, state_dict]],
    biome: Biome,
    longitude: float,
    years_to_plot: int | list[int] | np.ndarray | tuple[int, int] | None = None,
    contour_dict: dict[Biome, dict[Exp, state_dict]] | None = None,
    quiver_dicts: dict[Biome, dict[Exp, state_dict]] | None = None,
    quiver_axis: Literal["U", "V"] | None = None,
    ldyr: int = 0,
    title: str = "",
    figsize: tuple[float, float] = (10, 45),
    contourf_levels=None,
    contour_var: str | None = None,
    contour_var_levels=None,
    cmap: str = "viridis",
    dir_name=None,
    file_name=None,
    colorbar_label: str | None = None,
    season: str = "ANN",
    lev_interp: np.ndarray | None = None,
    lev_range: float | tuple[float, float] | None = None,
    return_fig_handles: bool = False,
    width: float = 0.005,
    headwidth: float = 5,
    headlength: float = 1,
    headaxislength: float = 2,
    font_size: int = 14,
    save: bool = False,
):

    """
    Plot latitude-depth cross-sections at a selected longitude for multiple datasets.

    The function creates one latitude-versus-depth cross-section for each dataset
    in ``ds_list`` using data associated with the selected ``biome``. The
    requested longitude is mapped to the nearest available longitude coordinate
    in each primary dataset. Temporal selection and compositing are handled by
    ``prepare_snapshot``, while optional depth interpolation and depth-range
    selection are handled by ``select_depth_range``.

    The primary field is displayed using filled contours. An optional secondary
    field may be overlaid as contour lines, and an optional vector component may
    be displayed using quiver arrows. :contentReference[oaicite:0]{index=0}

    Parameters
    ----------
    ds_list : list[Exp]
        Dataset or experiment names to plot. Each entry must correspond to a key
        in ``ds_dicts[biome]``. A separate subplot is created for each dataset in
        the order provided.
    ds_dicts : dict[Biome, dict[Exp, state_dict]]
        Nested mapping from biome names to dataset or experiment names and their
        associated state objects. The ``data`` attribute of each selected state
        object provides the primary field used to construct the cross-section.
    biome : Biome
        Biome whose latitude-depth cross-sections are plotted. The value must
        correspond to a key in ``ds_dicts`` and, when provided, in
        ``contour_dict`` and ``quiver_dicts``.
    longitude : float
        Longitude at which to extract the latitude-depth cross-section. For each
        primary dataset, the nearest available longitude coordinate is selected
        using xarray's nearest-neighbor selection.
    years_to_plot : int, list[int], numpy.ndarray, tuple[int, int], or None, optional
        Year selection passed to ``prepare_snapshot``. Depending on the helper
        implementation, this may specify a single year, multiple selected years,
        or a year range. If ``None``, the default temporal averaging behavior of
        ``prepare_snapshot`` is used.
    contour_dict : dict[Biome, dict[Exp, state_dict]] or None, optional
        Optional nested mapping containing a secondary variable to overlay as
        contour lines. The biome and experiment structure is expected to
        correspond to that of ``ds_dicts``. If provided,
        ``contour_var_levels`` must also be specified.
    quiver_dicts : dict[Biome, dict[Exp, state_dict]] or None, optional
        Optional nested mapping containing a field to display using quiver
        arrows. The field is processed using the same temporal, depth, and
        longitude selections as the primary field.
    quiver_axis : {"U", "V"} or None, optional
        Direction assigned to values from ``quiver_dicts``. If ``"U"``, the
        supplied values form the horizontal quiver component and the vertical
        component is set to zero. If ``"V"``, the supplied values form the
        vertical component and the horizontal component is set to zero. Must be
        ``"U"`` or ``"V"`` when ``quiver_dicts`` is provided.
    ldyr : int, optional
        Lead year passed to ``prepare_snapshot`` for temporal selection and
        seasonal compositing.
    title : str, optional
        Base text included in each subplot title. The selected longitude,
        dataset name, biome, temporal description, season, and optional contour
        variable information are appended automatically.
    figsize : tuple[float, float], optional
        Figure size passed to ``matplotlib.pyplot.subplots``.
    contourf_levels : array-like, int, or None, optional
        Levels supplied to ``matplotlib.axes.Axes.contourf`` for the primary
        field. If ``None``, Matplotlib determines the filled-contour levels
        automatically.
    contour_var : str or None, optional
        Descriptive name of the secondary variable supplied through
        ``contour_dict``. When provided together with ``contour_dict``, the name
        is included in the subplot title.
    contour_var_levels : array-like or int or None, optional
        Contour levels used for the secondary field supplied through
        ``contour_dict``. This argument is required whenever ``contour_dict`` is
        provided.
    cmap : str, optional
        Matplotlib colormap used for the filled contours of the primary field.
    dir_name : str or Path or None, optional
        Directory in which to save the figure when ``save=True``. The directory
        is created if it does not already exist.
    file_name : str or None, optional
        Output filename without the ``.png`` extension. Required when
        ``save=True``.
    colorbar_label : str or None, optional
        Label applied to the shared vertical colorbar. If ``None``, no colorbar
        label is added.
    season : str, optional
        Season passed to ``prepare_snapshot`` when constructing the temporal
        snapshot or composite.
    lev_interp : numpy.ndarray or None, optional
        Target depth coordinates passed to ``select_depth_range`` for vertical
        interpolation. If ``None``, native depth coordinates are retained.
    lev_range : float or tuple[float, float] or None, optional
        Depth interval passed to ``select_depth_range``. A scalar typically
        specifies a maximum depth, while a two-element tuple specifies an
        explicit minimum and maximum depth. If ``None``, the available depth
        range is retained.
    return_fig_handles : bool, optional
        If ``True``, return the Matplotlib figure and axes objects. If ``False``,
        the function returns nothing.
    width : float, optional
        Shaft width of quiver arrows.
    headwidth : float, optional
        Width of quiver arrow heads passed to
        ``matplotlib.axes.Axes.quiver``.
    headlength : float, optional
        Length of quiver arrow heads passed to
        ``matplotlib.axes.Axes.quiver``.
    headaxislength : float, optional
        Length of the quiver arrow head along its axis, passed to
        ``matplotlib.axes.Axes.quiver``.
    font_size : int, optional
        Font size used for subplot titles, axis labels, tick labels, contour
        labels, and colorbar text.
    save : bool, optional
        If ``True``, save the generated figure as a PNG file.

    Returns
    -------
    tuple[matplotlib.figure.Figure, numpy.ndarray] or None
        If ``return_fig_handles=True``, returns ``(fig, axes)``, where ``fig`` is
        the Matplotlib figure and ``axes`` is the two-dimensional array of subplot
        axes returned by ``plt.subplots``. Otherwise, returns ``None``.

    Raises
    ------
    ValueError
        If ``ds_list`` is empty.
    ValueError
        If ``quiver_dicts`` is provided and ``quiver_axis`` is neither ``"U"``
        nor ``"V"``.
    ValueError
        If ``contour_dict`` is provided without ``contour_var_levels``.
    ValueError
        If ``save=True`` and either ``dir_name`` or ``file_name`` is not
        provided.

    Notes
    -----
    For each dataset, the primary field is first processed with
    ``prepare_snapshot`` using the requested ``season``, ``ldyr``, and
    ``years_to_plot``. Depth interpolation and restriction are then applied
    through ``select_depth_range``.

    The requested ``longitude`` does not need to exactly match a grid coordinate.
    The primary field selects the nearest available longitude and stores that
    coordinate as ``selected_lon``. Optional contour and quiver fields are then
    selected at the grid longitude nearest to this selected primary-field
    longitude. The actual selected longitude is reported in the subplot title.

    When a secondary contour field is provided, it is aligned with the primary
    field using ``xr.align(..., join="inner")`` after temporal, depth, and
    longitude selection. Quiver data are aligned in the same manner. This
    restricts the respective fields to common latitude and depth coordinates
    before plotting.

    If ``contour_dict`` is supplied, the secondary variable is shown using
    semi-transparent black contour lines at ``contour_var_levels``. Otherwise,
    contour lines are generated from the primary field using the levels produced
    by the filled-contour plot and are displayed in white.

    Quiver data represent a single vector component. For ``quiver_axis="U"``,
    the supplied values define the horizontal component and the vertical
    component is zero. For ``quiver_axis="V"``, the supplied values define the
    vertical component and the horizontal component is zero. Quiver arrows are
    subsampled every second latitude and depth grid point, and NaN values are
    replaced by zero before plotting.

    The temporal description included in each subplot title is obtained from
    ``prepare_snapshot`` when available. If no year string is returned, the
    function falls back to the ``y0`` and ``y1`` metadata stored in the
    corresponding primary ``state_dict`` and identifies the period as a mean.

    The depth axis is inverted so that shallower depths appear near the top of
    each panel and greater depths appear toward the bottom.

    A single vertical colorbar is shared across all subplots and is based on the
    final filled-contour object created in the dataset loop. For direct
    comparison among experiments, supplying common ``contourf_levels`` across
    all datasets is recommended.

    Documentation produced with the assistance of AI.
    """
    if not ds_list:
        raise ValueError(
            "'ds_list' cannot be empty."
        )

    if quiver_dicts is not None and quiver_axis not in ("U", "V"):
        raise ValueError(
            "'quiver_axis' must be either 'U' or 'V' "
            "when quiver data are provided."
        )

    if contour_dict is not None and contour_var_levels is None:
        raise ValueError(
            "'contour_var_levels' must be provided when "
            "'contour_dict' is specified."
        )

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(
        len(ds_list),
        1,
        figsize=figsize,
        squeeze=False,
    )

    ds_dict = ds_dicts[biome]

    contour_f = None

    # ------------------------------------------------------------------
    # Datasets
    # ------------------------------------------------------------------
    for ds_idx, name in enumerate(ds_list):

        ax = axes[ds_idx, 0]

        # --------------------------------------------------------------
        # Primary field
        # --------------------------------------------------------------
        ts, year_string = prepare_snapshot(
            ds_dict[name].data,
            season=season,
            ldyr=ldyr,
            years_to_plot=years_to_plot,
            return_year_string=True,
        )

        ts = select_depth_range(
            ts,
            lev_interp=lev_interp,
            lev_range=lev_range,
        )

        ts = ts.sel(
            lon=longitude,
            method="nearest",
        )

        selected_lon = ts.lon.item()

        # --------------------------------------------------------------
        # Contour field
        # --------------------------------------------------------------
        contour_ts = None

        if contour_dict is not None:
            contour_ts = prepare_snapshot(
                contour_dict[biome][name].data,
                season=season,
                ldyr=ldyr,
                years_to_plot=years_to_plot,
            )

            contour_ts = select_depth_range(
                contour_ts,
                lev_interp=lev_interp,
                lev_range=lev_range,
            )

            contour_ts = contour_ts.sel(
                lon=selected_lon,
                method="nearest",
            )

        # --------------------------------------------------------------
        # Quiver field
        # --------------------------------------------------------------
        qv_ts = None

        if quiver_dicts is not None:
            qv_ts = prepare_snapshot(
                quiver_dicts[biome][name].data,
                season=season,
                ldyr=ldyr,
                years_to_plot=years_to_plot,
            )

            qv_ts = select_depth_range(
                qv_ts,
                lev_interp=lev_interp,
                lev_range=lev_range,
            )

            qv_ts = qv_ts.sel(
                lon=selected_lon,
                method="nearest",
            )

        # --------------------------------------------------------------
        # Align latitude / depth coordinates
        # --------------------------------------------------------------
        if contour_ts is not None:
            ts, contour_ts = xr.align(
                ts,
                contour_ts,
                join="inner",
            )

        if qv_ts is not None:
            ts, qv_ts = xr.align(
                ts,
                qv_ts,
                join="inner",
            )

        xx = ts.lat.values

        # --------------------------------------------------------------
        # Filled contours
        # --------------------------------------------------------------
        contour_f = ax.contourf(
            xx,
            ts.lev.values,
            ts,
            levels=contourf_levels,
            cmap=cmap,
        )

        # --------------------------------------------------------------
        # Contours
        # --------------------------------------------------------------
        if contour_ts is not None:

            contours = ax.contour(
                contour_ts.lat.values,
                contour_ts.lev.values,
                contour_ts,
                colors="black",
                levels=contour_var_levels,
                alpha=0.5,
            )

            ax.clabel(
                contours,
                inline=True,
                fontsize=font_size,
                colors="black",
            )

        else:
            contours = ax.contour(
                xx,
                ts.lev.values,
                ts,
                colors="white",
                levels=contour_f.levels,
            )

            ax.clabel(
                contours,
                inline=True,
                fontsize=font_size,
                colors="white",
            )

        # --------------------------------------------------------------
        # Quiver
        # --------------------------------------------------------------
        if qv_ts is not None:

            if quiver_axis == "U":
                U = qv_ts.values
                V = np.zeros_like(U)

            else:
                V = qv_ts.values
                U = np.zeros_like(V)

            ax.quiver(
                qv_ts.lat.values[::2],
                qv_ts.lev.values[::2],
                np.nan_to_num(
                    U[::2, ::2],
                    nan=0.0,
                ),
                np.nan_to_num(
                    V[::2, ::2],
                    nan=0.0,
                ),
                alpha=0.5,
                width=width,
                headwidth=headwidth,
                headlength=headlength,
                headaxislength=headaxislength,
            )

        # --------------------------------------------------------------
        # Year fallback
        # --------------------------------------------------------------
        if year_string is None:
            y0 = ds_dict[name].y0
            y1 = ds_dict[name].y1

            if y0 is not None and y1 is not None:
                year_string = f"{y0} - {y1} mean"
            else:
                year_string = ""

        # --------------------------------------------------------------
        # Title
        # --------------------------------------------------------------
        title_parts = [
            title,
            f"lon: {selected_lon:g}°E",
            name,
            '\n',
            biome,
            year_string,
            '\n',
            f"{season} composite",
        ]

        if contour_dict is not None and contour_var is not None:
            title_parts.append(
                f"{contour_var} contours"
            )

        title_ = " - ".join(
            part for part in title_parts if part
        )

        ax.set_title(
            title_,
            fontsize=font_size,
        )

        # --------------------------------------------------------------
        # Axis formatting
        # --------------------------------------------------------------
        ax.set_ylabel(
            "depth (m)",
            fontsize=font_size,
        )

        if ds_idx == len(ds_list) - 1:
            ax.set_xlabel(
                "Latitude (°)",
                fontsize=font_size,
            )
        else:
            ax.set_xlabel("")

        ax.invert_yaxis()

        ax.tick_params(
            axis="both",
            which="major",
            labelsize=font_size,
        )

    # ------------------------------------------------------------------
    # Shared colorbar
    # ------------------------------------------------------------------
    cbar = fig.colorbar(
        contour_f,
        ax=axes[:, 0],
        orientation="vertical",
    )

    if colorbar_label is not None:
        cbar.set_label(
            colorbar_label,
            fontsize=font_size,
        )

    cbar.ax.tick_params(
        labelsize=font_size,
    )

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    if save:
        if dir_name is None or file_name is None:
            raise ValueError(
                "'dir_name' and 'file_name' are required "
                "when save=True."
            )

        Path(dir_name).mkdir(
            parents=True,
            exist_ok=True,
        )

        plt.savefig(
            f"{dir_name}/{file_name}.png",
            bbox_inches="tight",
            dpi=300,
        )

    if return_fig_handles:
        return fig, axes
            

def snapshot_depth_vs_lon_quiver(
    ds_list: list[Exp],
    quiver_dicts_X: dict[Biome, dict[Exp, state_dict]] | None,
    quiver_dicts_Y: dict[Biome, dict[Exp, state_dict]] | None,
    biome: Biome,
    years_to_plot: int | list[int] | np.ndarray | tuple[int, int] | None = None,
    contour_dict: dict[Biome, dict[Exp, state_dict]] | None = None,
    ldyr: int = 0,
    title: str = "",
    figsize: tuple[float, float] = (10, 45),
    contourf_levels=None,
    cmap: str = "viridis",
    dir_name=None,
    file_name=None,
    colorbar_label: str | None = None,
    season: str = "ANN",
    lev_interp: np.ndarray | None = None,
    lev_range: float | tuple[float, float] | None = None,
    return_fig_handles: bool = False,
    save: bool = False,
    arrow_scale: float | None = None,
    quiver_step: int = 1,
    headwidth: float = 5,
    headlength: float = 1,
    headaxislength: float = 2,
    width: float = 0.005,
    font_size: int = 14,
    scale_factor_X: float = 1,
    scale_factor_Y: float = 1,
):

    """
    Plot longitude-depth vector fields for multiple datasets, optionally overlaid
    on a scalar contour field.

    The function creates one longitude-versus-depth panel for each dataset in
    ``ds_list``. Two independently supplied fields may be used as the horizontal
    and vertical components of the quiver vectors. Either component may be
    omitted; when only one is provided, the missing component is replaced with
    zeros so that a valid vector field can still be plotted.

    An optional scalar field from ``contour_dict`` can be displayed beneath the
    vectors using filled contours with overlaid contour lines. Temporal
    selection and compositing are handled by ``prepare_snapshot``, while
    vertical interpolation and depth-range selection are handled by
    ``select_depth_range``.

    Parameters
    ----------
    ds_list : list[Exp]
        Dataset or experiment names to plot. Each entry must correspond to a key
        under the selected ``biome`` in the supplied quiver dictionaries and, when
        provided, ``contour_dict``. A separate subplot is created for each dataset
        in the order given.
    quiver_dicts_X : dict[Biome, dict[Exp, state_dict]] or None
        Nested mapping containing the X-component of the quiver field. The
        ``data`` attribute of each selected state object is processed with
        ``prepare_snapshot`` and ``select_depth_range`` before plotting. At least
        one of ``quiver_dicts_X`` or ``quiver_dicts_Y`` must be provided.
    quiver_dicts_Y : dict[Biome, dict[Exp, state_dict]] or None
        Nested mapping containing the Y-component of the quiver field. The
        ``data`` attribute of each selected state object is processed in the same
        way as the X-component. At least one of ``quiver_dicts_X`` or
        ``quiver_dicts_Y`` must be provided.
    biome : Biome
        Biome whose longitude-depth vector fields are plotted. The value must
        correspond to a key in all supplied nested dictionaries.
    years_to_plot : int, list[int], numpy.ndarray, tuple[int, int], or None, optional
        Year selection passed to ``prepare_snapshot``. Depending on the helper
        implementation, this may specify a single year, multiple selected years,
        or a year range. If ``None``, the default temporal averaging behavior of
        ``prepare_snapshot`` is used.
    contour_dict : dict[Biome, dict[Exp, state_dict]] or None, optional
        Optional nested mapping containing a scalar field to display beneath the
        quiver vectors as filled contours. If ``None``, only the vector field is
        plotted and no colorbar is created.
    ldyr : int, optional
        Lead year passed to ``prepare_snapshot`` for temporal selection and
        seasonal compositing.
    title : str, optional
        Base text included at the beginning of each subplot title. Dataset name,
        biome, temporal description, and season are appended automatically.
    figsize : tuple[float, float], optional
        Figure size passed to ``matplotlib.pyplot.subplots``.
    contourf_levels : array-like, int, or None, optional
        Levels supplied to ``matplotlib.axes.Axes.contourf`` when
        ``contour_dict`` is provided. If ``None``, Matplotlib determines the
        contour levels automatically.
    cmap : str, optional
        Matplotlib colormap used for the optional scalar contour field.
    dir_name : str or Path or None, optional
        Directory in which to save the figure when ``save=True``. The directory
        is created if it does not already exist.
    file_name : str or None, optional
        Output filename without the ``.png`` extension. Required when
        ``save=True``.
    colorbar_label : str or None, optional
        Label applied to the shared vertical colorbar when ``contour_dict`` is
        provided. If ``None``, no label is added.
    season : str, optional
        Season passed to ``prepare_snapshot`` when constructing the temporal
        snapshot or composite.
    lev_interp : numpy.ndarray or None, optional
        Target depth coordinates passed to ``select_depth_range`` for vertical
        interpolation. If ``None``, native depth coordinates are retained.
    lev_range : float or tuple[float, float] or None, optional
        Depth interval passed to ``select_depth_range``. A scalar typically
        specifies a maximum depth, while a two-element tuple specifies an
        explicit minimum and maximum depth. If ``None``, the available depth
        range is retained.
    return_fig_handles : bool, optional
        If ``True``, return the Matplotlib figure and axes objects. If ``False``,
        the function returns nothing.
    save : bool, optional
        If ``True``, save the generated figure as a PNG file.
    arrow_scale : float or None, optional
        Scaling factor passed to ``matplotlib.axes.Axes.quiver`` through its
        ``scale`` argument. Smaller values generally produce longer arrows, while
        larger values produce shorter arrows. If ``None``, Matplotlib determines
        the scale automatically.
    quiver_step : int, optional
        Spatial subsampling interval applied independently to longitude and depth
        when plotting quiver arrows. A value of ``1`` plots every vector, ``2``
        plots every second vector along each dimension, and so on. Must be at
        least 1.
    headwidth : float, optional
        Width of quiver arrow heads passed to
        ``matplotlib.axes.Axes.quiver``.
    headlength : float, optional
        Length of quiver arrow heads passed to
        ``matplotlib.axes.Axes.quiver``.
    headaxislength : float, optional
        Length of the quiver arrow head along its axis, passed to
        ``matplotlib.axes.Axes.quiver``.
    width : float, optional
        Shaft width of quiver arrows.
    font_size : int, optional
        Font size used for subplot titles, axis labels, tick labels, contour
        labels, and colorbar text.
    scale_factor_X : float, optional
        Multiplicative scaling factor applied to the X-component after temporal
        and depth processing and before alignment and plotting.
    scale_factor_Y : float, optional
        Multiplicative scaling factor applied to the Y-component after temporal
        and depth processing and before alignment and plotting.

    Returns
    -------
    tuple[matplotlib.figure.Figure, numpy.ndarray] or None
        If ``return_fig_handles=True``, returns ``(fig, axes)``, where ``fig`` is
        the Matplotlib figure and ``axes`` is the two-dimensional array of subplot
        axes returned by ``plt.subplots``. Otherwise, returns ``None``.

    Raises
    ------
    ValueError
        If ``ds_list`` is empty.
    ValueError
        If both ``quiver_dicts_X`` and ``quiver_dicts_Y`` are ``None``.
    ValueError
        If ``quiver_step`` is less than 1.
    ValueError
        If ``save=True`` and either ``dir_name`` or ``file_name`` is not
        provided.

    Notes
    -----
    Each supplied quiver component is processed independently with
    ``prepare_snapshot`` using the requested ``season``, ``ldyr``, and
    ``years_to_plot`` arguments. The resulting field is then passed through
    ``select_depth_range`` using the same ``lev_interp`` and ``lev_range``
    settings.

    The optional ``scale_factor_X`` and ``scale_factor_Y`` values are applied
    after temporal and vertical preprocessing. These factors can be used to
    rescale components with different units or magnitudes before vector plotting.

    If one quiver component is omitted, it is replaced using
    ``xr.zeros_like`` based on the provided component. For example, if
    ``quiver_dicts_X`` is ``None``, the X-component is set to zero everywhere
    using the processed Y-component as the template. This allows purely
    horizontal or purely vertical vector fields to be represented.

    Before plotting, the X- and Y-components are aligned using
    ``xr.align(..., join="inner")``. If a contour field is also provided, all
    three fields are aligned together. This restricts them to common longitude
    and depth coordinates.

    When ``contour_dict`` is supplied, the scalar field is shown using filled
    contours and white contour lines at the same levels. The contour lines are
    labeled directly on the plot. A shared vertical colorbar is created only in
    this case.

    Quiver vectors are subsampled using ``quiver_step`` along both longitude and
    depth. NaN values in either component are replaced with zero immediately
    before plotting.

    The temporal description used in the subplot title is obtained first from
    the X-component when available, otherwise from the Y-component. If no year
    string is returned by ``prepare_snapshot``, the function falls back to the
    ``y0`` and ``y1`` metadata of whichever quiver source dictionary is
    available and labels that interval as a mean.

    The depth axis is inverted so that shallow depths appear near the top of
    each panel and greater depths appear toward the bottom.

    Documentation produced with the assistance of AI.
    """
    if not ds_list:
        raise ValueError(
            "'ds_list' cannot be empty."
        )

    if quiver_dicts_X is None and quiver_dicts_Y is None:
        raise ValueError(
            "At least one of 'quiver_dicts_X' or "
            "'quiver_dicts_Y' must be provided."
        )

    if quiver_step < 1:
        raise ValueError(
            "'quiver_step' must be at least 1."
        )

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(
        len(ds_list),
        1,
        figsize=figsize,
        squeeze=False,
    )

    contour_f = None

    # ------------------------------------------------------------------
    # Datasets
    # ------------------------------------------------------------------
    for ds_idx, name in enumerate(ds_list):

        ax = axes[ds_idx, 0]

        # --------------------------------------------------------------
        # Quiver X component
        # --------------------------------------------------------------
        tsX = None
        year_string = None

        if quiver_dicts_X is not None:
            tsX, year_string = prepare_snapshot(
                quiver_dicts_X[biome][name].data,
                season=season,
                ldyr=ldyr,
                years_to_plot=years_to_plot,
                return_year_string=True,
            )

            tsX = select_depth_range(
                tsX,
                lev_interp=lev_interp,
                lev_range=lev_range,
            )

            tsX = tsX * scale_factor_X

        # --------------------------------------------------------------
        # Quiver Y component
        # --------------------------------------------------------------
        tsY = None

        if quiver_dicts_Y is not None:
            tsY, year_string_Y = prepare_snapshot(
                quiver_dicts_Y[biome][name].data,
                season=season,
                ldyr=ldyr,
                years_to_plot=years_to_plot,
                return_year_string=True,
            )

            tsY = select_depth_range(
                tsY,
                lev_interp=lev_interp,
                lev_range=lev_range,
            )

            tsY = tsY * scale_factor_Y

            if year_string is None:
                year_string = year_string_Y


        # --------------------------------------------------------------
        # Missing quiver component
        # --------------------------------------------------------------
        if tsX is None:
            tsX = xr.zeros_like(tsY)

        if tsY is None:
            tsY = xr.zeros_like(tsX)

        # --------------------------------------------------------------
        # Contour field
        # --------------------------------------------------------------
        contour_ts = None

        if contour_dict is not None:
            contour_ts = prepare_snapshot(
                contour_dict[biome][name].data,
                season=season,
                ldyr=ldyr,
                years_to_plot=years_to_plot,
            )

            contour_ts = select_depth_range(
                contour_ts,
                lev_interp=lev_interp,
                lev_range=lev_range,
            )

        # --------------------------------------------------------------
        # Align longitude / depth coordinates
        # --------------------------------------------------------------
        if contour_ts is not None:
            tsX, tsY, contour_ts = xr.align(
                tsX,
                tsY,
                contour_ts,
                join="inner",
            )
        else:
            tsX, tsY = xr.align(
                tsX,
                tsY,
                join="inner",
            )

        # --------------------------------------------------------------
        # Filled contours
        # --------------------------------------------------------------
        if contour_ts is not None:
            contour_f = ax.contourf(
                contour_ts.lon.values,
                contour_ts.lev.values,
                contour_ts,
                levels=contourf_levels,
                cmap=cmap,
            )

            contours = ax.contour(
                contour_ts.lon.values,
                contour_ts.lev.values,
                contour_ts,
                colors="white",
                levels=contour_f.levels,
            )

            ax.clabel(
                contours,
                inline=True,
                fontsize=font_size,
                colors="white",
            )

        # --------------------------------------------------------------
        # Quiver
        # --------------------------------------------------------------
        ax.quiver(
            tsX.lon.values[::quiver_step],
            tsX.lev.values[::quiver_step],
            np.nan_to_num(
                tsX.values[::quiver_step, ::quiver_step],
                nan=0.0,
            ),
            np.nan_to_num(
                tsY.values[::quiver_step, ::quiver_step],
                nan=0.0,
            ),
            alpha=0.5,
            scale=arrow_scale,
            width=width,
            headwidth=headwidth,
            headlength=headlength,
            headaxislength=headaxislength,
        )

        # --------------------------------------------------------------
        # Year fallback
        # --------------------------------------------------------------
        if year_string is None:

            source_dict = (
                quiver_dicts_X
                if quiver_dicts_X is not None
                else quiver_dicts_Y
            )

            y0 = source_dict[biome][name].y0
            y1 = source_dict[biome][name].y1

            if y0 is not None and y1 is not None:
                year_string = f"{y0} - {y1} mean"
            else:
                year_string = ""

        # --------------------------------------------------------------
        # Title
        # --------------------------------------------------------------
        title_parts = [
            title,
            name,
            biome,
            year_string,
            f"{season} composite",
        ]

        title_ = " - ".join(
            part for part in title_parts if part
        )

        ax.set_title(
            title_,
            fontsize=font_size,
        )

        # --------------------------------------------------------------
        # Axis formatting
        # --------------------------------------------------------------
        ax.set_ylabel(
            "depth (m)",
            fontsize=font_size,
        )

        if ds_idx == len(ds_list) - 1:
            ax.set_xlabel(
                "Lon ($^o$ East)",
                fontsize=font_size,
            )
        else:
            ax.set_xlabel("")

        ax.invert_yaxis()

        ax.tick_params(
            axis="both",
            which="major",
            labelsize=font_size,
        )

    # ------------------------------------------------------------------
    # Shared colorbar
    # ------------------------------------------------------------------
    if contour_f is not None:
        cbar = fig.colorbar(
            contour_f,
            ax=axes[:, 0],
            orientation="vertical",
        )

        if colorbar_label is not None:
            cbar.set_label(
                colorbar_label,
                fontsize=font_size,
            )

        cbar.ax.tick_params(
            labelsize=font_size,
        )

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    if save:
        if dir_name is None or file_name is None:
            raise ValueError(
                "'dir_name' and 'file_name' are required "
                "when save=True."
            )

        Path(dir_name).mkdir(
            parents=True,
            exist_ok=True,
        )

        plt.savefig(
            f"{dir_name}/{file_name}.png",
            bbox_inches="tight",
            dpi=300,
        )

    if return_fig_handles:
        return fig, axes



def snapshot_crossection_quiver(
    ds_list: list[Exp],
    quiver_dicts_X: dict[Biome, dict[Exp, state_dict]] | None,
    quiver_dicts_Y: dict[Biome, dict[Exp, state_dict]] | None,
    biome: Biome,
    longitude: float,
    years_to_plot: int | list[int] | np.ndarray | tuple[int, int] | None = None,
    contour_dict: dict[Biome, dict[Exp, state_dict]] | None = None,
    ldyr: int = 0,
    title: str = "",
    figsize: tuple[float, float] = (10, 45),
    contourf_levels=None,
    cmap: str = "viridis",
    dir_name=None,
    file_name=None,
    colorbar_label: str | None = None,
    season: str = "ANN",
    lev_interp: np.ndarray | None = None,
    lev_range: float | tuple[float, float] | None = None,
    return_fig_handles: bool = False,
    save: bool = False,
    arrow_scale: float | None = None,
    quiver_step: int = 1,
    headwidth: float = 5,
    headlength: float = 1,
    headaxislength: float = 2,
    width: float = 0.005,
    font_size: int = 15,
    scale_factor_X: float = 1,
    scale_factor_Y: float = 1,
):

    """
    Plot latitude-depth vector cross-sections at a selected longitude for multiple
    datasets, optionally overlaid on a scalar contour field.

    The function creates one latitude-versus-depth panel for each dataset in
    ``ds_list``. Two independently supplied fields may be used as the horizontal
    and vertical components of the quiver vectors. Either component may be
    omitted; when only one is provided, the missing component is replaced with
    zeros so that a valid vector field can still be plotted.

    For each dataset, the requested longitude is mapped to the nearest available
    longitude coordinate using whichever quiver component is available as the
    reference. The same selected longitude is then used for the remaining quiver
    component and optional contour field. Temporal selection and compositing are
    handled by ``prepare_snapshot``, while vertical interpolation and depth-range
    selection are handled by ``select_depth_range``. :contentReference[oaicite:0]{index=0}

    Parameters
    ----------
    ds_list : list[Exp]
        Dataset or experiment names to plot. Each entry must correspond to a key
        under the selected ``biome`` in the supplied quiver dictionaries and, when
        provided, ``contour_dict``. A separate subplot is created for each dataset
        in the order given.
    quiver_dicts_X : dict[Biome, dict[Exp, state_dict]] or None
        Nested mapping containing the X-component of the quiver field. The
        ``data`` attribute of each selected state object is processed with
        ``prepare_snapshot`` and ``select_depth_range`` before longitude
        selection and plotting. At least one of ``quiver_dicts_X`` or
        ``quiver_dicts_Y`` must be provided.
    quiver_dicts_Y : dict[Biome, dict[Exp, state_dict]] or None
        Nested mapping containing the Y-component of the quiver field. The
        ``data`` attribute of each selected state object is processed in the same
        way as the X-component. At least one of ``quiver_dicts_X`` or
        ``quiver_dicts_Y`` must be provided.
    biome : Biome
        Biome whose latitude-depth vector cross-sections are plotted. The value
        must correspond to a key in all supplied nested dictionaries.
    longitude : float
        Requested longitude at which to extract the latitude-depth cross-section.
        The nearest available longitude is selected from whichever processed
        quiver component is available first and is then used as the reference
        longitude for all other plotted fields.
    years_to_plot : int, list[int], numpy.ndarray, tuple[int, int], or None, optional
        Year selection passed to ``prepare_snapshot``. Depending on the helper
        implementation, this may specify a single year, multiple selected years,
        or a year range. If ``None``, the default temporal averaging behavior of
        ``prepare_snapshot`` is used.
    contour_dict : dict[Biome, dict[Exp, state_dict]] or None, optional
        Optional nested mapping containing a scalar field to display beneath the
        quiver vectors as filled contours. The field is processed using the same
        temporal, depth, and longitude selections as the vector components.
    ldyr : int, optional
        Lead year passed to ``prepare_snapshot`` for temporal selection and
        seasonal compositing.
    title : str, optional
        Base text included at the beginning of each subplot title. The selected
        longitude, dataset name, biome, temporal description, and season are
        appended automatically.
    figsize : tuple[float, float], optional
        Figure size passed to ``matplotlib.pyplot.subplots``.
    contourf_levels : array-like, int, or None, optional
        Levels supplied to ``matplotlib.axes.Axes.contourf`` when
        ``contour_dict`` is provided. If ``None``, Matplotlib determines the
        contour levels automatically.
    cmap : str, optional
        Matplotlib colormap used for the optional scalar contour field.
    dir_name : str or Path or None, optional
        Directory in which to save the figure when ``save=True``. The directory
        is created if it does not already exist.
    file_name : str or None, optional
        Output filename without the ``.png`` extension. Required when
        ``save=True``.
    colorbar_label : str or None, optional
        Label applied to the shared vertical colorbar when ``contour_dict`` is
        provided. If ``None``, no colorbar label is added.
    season : str, optional
        Season passed to ``prepare_snapshot`` when constructing the temporal
        snapshot or composite.
    lev_interp : numpy.ndarray or None, optional
        Target depth coordinates passed to ``select_depth_range`` for vertical
        interpolation. If ``None``, native depth coordinates are retained.
    lev_range : float or tuple[float, float] or None, optional
        Depth interval passed to ``select_depth_range``. A scalar typically
        specifies a maximum depth, while a two-element tuple specifies an
        explicit minimum and maximum depth. If ``None``, the available depth
        range is retained.
    return_fig_handles : bool, optional
        If ``True``, return the Matplotlib figure and axes objects. If ``False``,
        the function returns nothing.
    save : bool, optional
        If ``True``, save the generated figure as a PNG file.
    arrow_scale : float or None, optional
        Scaling factor passed to ``matplotlib.axes.Axes.quiver`` through its
        ``scale`` argument. Smaller values generally produce longer arrows, while
        larger values produce shorter arrows. If ``None``, Matplotlib determines
        the scale automatically.
    quiver_step : int, optional
        Spatial subsampling interval applied independently to latitude and depth
        when plotting quiver arrows. A value of ``1`` plots every vector, ``2``
        plots every second vector along each dimension, and so on. Must be at
        least 1.
    headwidth : float, optional
        Width of quiver arrow heads passed to
        ``matplotlib.axes.Axes.quiver``.
    headlength : float, optional
        Length of quiver arrow heads passed to
        ``matplotlib.axes.Axes.quiver``.
    headaxislength : float, optional
        Length of the quiver arrow head along its axis, passed to
        ``matplotlib.axes.Axes.quiver``.
    width : float, optional
        Shaft width of quiver arrows.
    font_size : int, optional
        Font size used for subplot titles, axis labels, tick labels, contour
        labels, and colorbar text.
    scale_factor_X : float, optional
        Multiplicative scaling factor applied to the X-component after temporal
        and depth processing and before longitude selection, alignment, and
        plotting.
    scale_factor_Y : float, optional
        Multiplicative scaling factor applied to the Y-component after temporal
        and depth processing and before longitude selection, alignment, and
        plotting.

    Returns
    -------
    tuple[matplotlib.figure.Figure, numpy.ndarray] or None
        If ``return_fig_handles=True``, returns ``(fig, axes)``, where ``fig`` is
        the Matplotlib figure and ``axes`` is the two-dimensional array of subplot
        axes returned by ``plt.subplots``. Otherwise, returns ``None``.

    Raises
    ------
    ValueError
        If ``ds_list`` is empty.
    ValueError
        If both ``quiver_dicts_X`` and ``quiver_dicts_Y`` are ``None``.
    ValueError
        If ``quiver_step`` is less than 1.
    ValueError
        If ``save=True`` and either ``dir_name`` or ``file_name`` is not
        provided.

    Notes
    -----
    Each supplied quiver component is processed independently with
    ``prepare_snapshot`` using the requested ``season``, ``ldyr``, and
    ``years_to_plot`` arguments. The resulting field is then passed through
    ``select_depth_range`` using the same ``lev_interp`` and ``lev_range``
    settings.

    The optional ``scale_factor_X`` and ``scale_factor_Y`` values are applied
    after temporal and vertical preprocessing. These factors can be used to
    rescale vector components before plotting.

    Longitude selection is based on whichever quiver component is available as
    the reference field. The nearest longitude to ``longitude`` is selected from
    that field and stored as ``selected_lon``. Both quiver components, when
    present, and the optional contour field are then independently selected at
    the grid longitude nearest to ``selected_lon``. The actual selected longitude
    is included in the subplot title.

    If one quiver component is omitted, it is replaced using
    ``xr.zeros_like`` based on the provided component. For example, if
    ``quiver_dicts_X`` is ``None``, the X-component is set to zero everywhere
    using the processed and longitude-selected Y-component as the template.

    Before plotting, the X- and Y-components are aligned using
    ``xr.align(..., join="inner")``. If a contour field is provided, all three
    fields are aligned together. This restricts them to common latitude and depth
    coordinates.

    When ``contour_dict`` is supplied, the scalar field is displayed with filled
    contours and white contour lines using the levels generated by ``contourf``.
    The contour lines are labeled directly on the plot. A shared vertical
    colorbar is created only when this scalar contour field is present.

    Quiver vectors are subsampled using ``quiver_step`` along both latitude and
    depth. NaN values in either component are replaced with zero immediately
    before plotting.

    The temporal description used in the subplot title is obtained first from
    the X-component when available, otherwise from the Y-component. If no year
    string is returned by ``prepare_snapshot``, the function falls back to the
    ``y0`` and ``y1`` metadata of whichever quiver source dictionary is
    available and labels that interval as a mean.

    The depth axis is inverted so that shallow depths appear near the top of
    each panel and greater depths appear toward the bottom.

    Documentation produced with the assistance of AI.
    """

    if not ds_list:
        raise ValueError(
            "'ds_list' cannot be empty."
        )

    if quiver_dicts_X is None and quiver_dicts_Y is None:
        raise ValueError(
            "At least one of 'quiver_dicts_X' or "
            "'quiver_dicts_Y' must be provided."
        )

    if quiver_step < 1:
        raise ValueError(
            "'quiver_step' must be at least 1."
        )

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(
        len(ds_list),
        1,
        figsize=figsize,
        squeeze=False,
    )

    contour_f = None

    # ------------------------------------------------------------------
    # Datasets
    # ------------------------------------------------------------------
    for ds_idx, name in enumerate(ds_list):

        ax = axes[ds_idx, 0]

        # --------------------------------------------------------------
        # Quiver X component
        # --------------------------------------------------------------
        tsX = None
        year_string = None

        if quiver_dicts_X is not None:
            tsX, year_string = prepare_snapshot(
                quiver_dicts_X[biome][name].data,
                season=season,
                ldyr=ldyr,
                years_to_plot=years_to_plot,
                return_year_string=True,
            )

            tsX = select_depth_range(
                tsX,
                lev_interp=lev_interp,
                lev_range=lev_range,
            )

            tsX = tsX * scale_factor_X

        # --------------------------------------------------------------
        # Quiver Y component
        # --------------------------------------------------------------
        tsY = None

        if quiver_dicts_Y is not None:
            tsY, year_string_Y = prepare_snapshot(
                quiver_dicts_Y[biome][name].data,
                season=season,
                ldyr=ldyr,
                years_to_plot=years_to_plot,
                return_year_string=True,
            )

            tsY = select_depth_range(
                tsY,
                lev_interp=lev_interp,
                lev_range=lev_range,
            )

            tsY = tsY * scale_factor_Y

            if year_string is None:
                year_string = year_string_Y


        # --------------------------------------------------------------
        # Select longitude
        #
        # Use whichever quiver component exists as the reference for the
        # actual nearest longitude.
        # --------------------------------------------------------------
        reference_ts = (
            tsX
            if tsX is not None
            else tsY
        )

        reference_ts = reference_ts.sel(
            lon=longitude,
            method="nearest",
        )

        selected_lon = reference_ts.lon.item()

        if tsX is not None:
            tsX = tsX.sel(
                lon=selected_lon,
                method="nearest",
            )

        if tsY is not None:
            tsY = tsY.sel(
                lon=selected_lon,
                method="nearest",
            )

        # --------------------------------------------------------------
        # Missing quiver component
        # --------------------------------------------------------------
        if tsX is None:
            tsX = xr.zeros_like(tsY)

        if tsY is None:
            tsY = xr.zeros_like(tsX)

        # --------------------------------------------------------------
        # Contour field
        # --------------------------------------------------------------
        contour_ts = None

        if contour_dict is not None:
            contour_ts = prepare_snapshot(
                contour_dict[biome][name].data,
                season=season,
                ldyr=ldyr,
                years_to_plot=years_to_plot,
            )

            contour_ts = select_depth_range(
                contour_ts,
                lev_interp=lev_interp,
                lev_range=lev_range,
            )

            contour_ts = contour_ts.sel(
                lon=selected_lon,
                method="nearest",
            )

        # --------------------------------------------------------------
        # Align latitude / depth coordinates
        # --------------------------------------------------------------
        if contour_ts is not None:
            tsX, tsY, contour_ts = xr.align(
                tsX,
                tsY,
                contour_ts,
                join="inner",
            )

        else:
            tsX, tsY = xr.align(
                tsX,
                tsY,
                join="inner",
            )

        # --------------------------------------------------------------
        # Filled contours
        # --------------------------------------------------------------
        if contour_ts is not None:
            contour_f = ax.contourf(
                contour_ts.lat.values,
                contour_ts.lev.values,
                contour_ts,
                levels=contourf_levels,
                cmap=cmap,
            )

            contours = ax.contour(
                contour_ts.lat.values,
                contour_ts.lev.values,
                contour_ts,
                colors="white",
                levels=contour_f.levels,
            )

            ax.clabel(
                contours,
                inline=True,
                fontsize=font_size,
                colors="white",
            )

        # --------------------------------------------------------------
        # Quiver
        # --------------------------------------------------------------
        ax.quiver(
            tsX.lat.values[::quiver_step],
            tsX.lev.values[::quiver_step],
            np.nan_to_num(
                tsX.values[::quiver_step, ::quiver_step],
                nan=0.0,
            ),
            np.nan_to_num(
                tsY.values[::quiver_step, ::quiver_step],
                nan=0.0,
            ),
            scale=arrow_scale,
            width=width,
            headwidth=headwidth,
            headlength=headlength,
            headaxislength=headaxislength,
        )

        # --------------------------------------------------------------
        # Year fallback
        # --------------------------------------------------------------
        if year_string is None:

            source_dict = (
                quiver_dicts_X
                if quiver_dicts_X is not None
                else quiver_dicts_Y
            )

            y0 = source_dict[biome][name].y0
            y1 = source_dict[biome][name].y1

            if y0 is not None and y1 is not None:
                year_string = f"{y0} - {y1} mean"
            else:
                year_string = ""

        # --------------------------------------------------------------
        # Title
        # --------------------------------------------------------------
        title_parts = [
            title,
            f"lon: {selected_lon:g}°E",
            name,
            '\n',
            biome,
            year_string,
            '\n',
            f"{season} composite",
        ]

        title_ = " - ".join(
            part for part in title_parts if part
        )

        ax.set_title(
            title_,
            fontsize=font_size,
        )

        # --------------------------------------------------------------
        # Axis formatting
        # --------------------------------------------------------------
        ax.set_ylabel(
            "depth (m)",
            fontsize=font_size,
        )

        if ds_idx == len(ds_list) - 1:
            ax.set_xlabel(
                "Latitude ($^o$)",
                fontsize=font_size,
            )
        else:
            ax.set_xlabel("")

        ax.invert_yaxis()

        ax.tick_params(
            axis="both",
            which="major",
            labelsize=font_size,
        )

    # ------------------------------------------------------------------
    # Shared colorbar
    # ------------------------------------------------------------------
    if contour_f is not None:
        cbar = fig.colorbar(
            contour_f,
            ax=axes[:, 0],
            orientation="vertical",
        )

        if colorbar_label is not None:
            cbar.set_label(
                colorbar_label,
                fontsize=font_size,
            )

        cbar.ax.tick_params(
            labelsize=font_size,
        )

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    if save:
        if dir_name is None or file_name is None:
            raise ValueError(
                "'dir_name' and 'file_name' are required "
                "when save=True."
            )

        Path(dir_name).mkdir(
            parents=True,
            exist_ok=True,
        )

        plt.savefig(
            f"{dir_name}/{file_name}.png",
            bbox_inches="tight",
            dpi=300,
        )

    if return_fig_handles:
        return fig, axes

            

def snapshot_aerial_quiver(
    ds_list: list[Exp],
    quiver_dicts_X: dict[Biome, dict[Exp, state_dict]] | None,
    quiver_dicts_Y: dict[Biome, dict[Exp, state_dict]] | None,
    biome: Biome,
    years_to_plot: int | list[int] | np.ndarray | tuple[int, int] | None = None,
    background_dict: dict[Biome, dict[Exp, state_dict]] | None = None,
    ldyr: int = 0,
    title: str = "",
    figsize: tuple[float, float] = (10, 45),
    vmax: float | None = None,
    vmin: float | None = None,
    cmap: str = "viridis",
    dir_name=None,
    file_name=None,
    colorbar_label: str | None = None,
    season: str = "ANN",
    lev_interp: np.ndarray | None = None,
    lev_range: float | tuple[float, float] | None = None,
    return_fig_handles: bool = False,
    save: bool = False,
    arrow_scale: float | None = None,
    quiver_step: int = 1,
    headwidth: float = 5,
    headlength: float = 1,
    headaxislength: float = 2,
    width: float = 0.005,
    font_size: int = 15,
    scale_factor_X: float = 1,
    scale_factor_Y: float = 1,
):

    """
    Plot horizontal vector fields for multiple datasets, optionally overlaid on a
    scalar background field.

    The function creates one latitude-longitude panel for each dataset in
    ``ds_list``. Two independently supplied fields may be used as the X and Y
    components of the quiver vectors. Either component may be omitted; when only
    one is provided, the missing component is replaced with zeros so that a valid
    vector field can still be plotted.

    If the quiver or background fields contain a depth dimension, the requested
    depth range can optionally be interpolated and selected before the fields are
    averaged over depth. Temporal selection and compositing are handled by
    ``prepare_snapshot``, while vertical processing is handled by
    ``select_depth_range``.

    Parameters
    ----------
    ds_list : list[Exp]
        Dataset or experiment names to plot. Each entry must correspond to a key
        under the selected ``biome`` in the supplied quiver dictionaries and, when
        provided, ``background_dict``. A separate subplot is created for each
        dataset in the order given.
    quiver_dicts_X : dict[Biome, dict[Exp, state_dict]] or None
        Nested mapping containing the X-component of the quiver field. The
        ``data`` attribute of each selected state object is processed with
        ``prepare_snapshot`` and, when a depth dimension is present,
        ``select_depth_range`` before depth averaging and plotting. At least one
        of ``quiver_dicts_X`` or ``quiver_dicts_Y`` must be provided.
    quiver_dicts_Y : dict[Biome, dict[Exp, state_dict]] or None
        Nested mapping containing the Y-component of the quiver field. The
        ``data`` attribute of each selected state object is processed in the same
        way as the X-component. At least one of ``quiver_dicts_X`` or
        ``quiver_dicts_Y`` must be provided.
    biome : Biome
        Biome whose horizontal vector fields are plotted. The value must
        correspond to a key in all supplied nested dictionaries.
    years_to_plot : int, list[int], numpy.ndarray, tuple[int, int], or None, optional
        Year selection passed to ``prepare_snapshot``. Depending on the helper
        implementation, this may specify a single year, multiple selected years,
        or a year range. If ``None``, the default temporal averaging behavior of
        ``prepare_snapshot`` is used.
    background_dict : dict[Biome, dict[Exp, state_dict]] or None, optional
        Optional nested mapping containing a scalar field to display beneath the
        quiver vectors using ``pcolormesh``. If the field contains a depth
        dimension, the same depth interpolation and range selection are applied
        before averaging over depth.
    ldyr : int, optional
        Lead year passed to ``prepare_snapshot`` for temporal selection and
        seasonal compositing.
    title : str, optional
        Base text included at the beginning of each subplot title. Dataset name,
        biome, temporal description, season, and selected depth range are appended
        automatically where applicable.
    figsize : tuple[float, float], optional
        Figure size passed to ``matplotlib.pyplot.subplots``.
    vmax : float or None, optional
        Upper bound of the color scale used for the optional background field. If
        ``None``, Matplotlib determines the upper limit automatically.
    vmin : float or None, optional
        Lower bound of the color scale used for the optional background field. If
        ``None``, Matplotlib determines the lower limit automatically.
    cmap : str, optional
        Matplotlib colormap used for the optional scalar background field.
    dir_name : str or Path or None, optional
        Directory in which to save the figure when ``save=True``. The directory
        is created if it does not already exist.
    file_name : str or None, optional
        Output filename without the ``.png`` extension. Required when
        ``save=True``.
    colorbar_label : str or None, optional
        Label applied to the shared vertical colorbar when ``background_dict`` is
        provided. If ``None``, no colorbar label is added.
    season : str, optional
        Season passed to ``prepare_snapshot`` when constructing the temporal
        snapshot or composite.
    lev_interp : numpy.ndarray or None, optional
        Target depth coordinates passed to ``select_depth_range`` for vertical
        interpolation when a ``lev`` dimension is present. If ``None``, native
        depth coordinates are retained.
    lev_range : float or tuple[float, float] or None, optional
        Depth interval passed to ``select_depth_range`` when a ``lev`` dimension
        is present. A scalar typically specifies a maximum depth, while a
        two-element tuple specifies an explicit minimum and maximum depth. If
        ``None``, the available depth range is retained.
    return_fig_handles : bool, optional
        If ``True``, return the Matplotlib figure and axes objects. If ``False``,
        the function returns nothing.
    save : bool, optional
        If ``True``, save the generated figure as a PNG file.
    arrow_scale : float or None, optional
        Scaling factor passed to ``matplotlib.axes.Axes.quiver`` through its
        ``scale`` argument. Smaller values generally produce longer arrows, while
        larger values produce shorter arrows. If ``None``, Matplotlib determines
        the scale automatically.
    quiver_step : int, optional
        Spatial subsampling interval applied independently to longitude and
        latitude when plotting quiver arrows. A value of ``1`` plots every
        vector, ``2`` plots every second vector along each dimension, and so on.
        Must be at least 1.
    headwidth : float, optional
        Width of quiver arrow heads passed to
        ``matplotlib.axes.Axes.quiver``.
    headlength : float, optional
        Length of quiver arrow heads passed to
        ``matplotlib.axes.Axes.quiver``.
    headaxislength : float, optional
        Length of the quiver arrow head along its axis, passed to
        ``matplotlib.axes.Axes.quiver``.
    width : float, optional
        Shaft width of quiver arrows.
    font_size : int, optional
        Font size used for subplot titles, axis labels, tick labels, and colorbar
        text.
    scale_factor_X : float, optional
        Multiplicative scaling factor applied to the X-component after temporal
        and optional depth processing and before alignment and plotting.
    scale_factor_Y : float, optional
        Multiplicative scaling factor applied to the Y-component after temporal
        and optional depth processing and before alignment and plotting.

    Returns
    -------
    tuple[matplotlib.figure.Figure, numpy.ndarray] or None
        If ``return_fig_handles=True``, returns ``(fig, axes)``, where ``fig`` is
        the Matplotlib figure and ``axes`` is the two-dimensional array of subplot
        axes returned by ``plt.subplots``. Otherwise, returns ``None``.

    Raises
    ------
    ValueError
        If ``ds_list`` is empty.
    ValueError
        If both ``quiver_dicts_X`` and ``quiver_dicts_Y`` are ``None``.
    ValueError
        If ``quiver_step`` is less than 1.
    ValueError
        If ``save=True`` and either ``dir_name`` or ``file_name`` is not
        provided.

    Notes
    -----
    Each supplied quiver component is processed independently with
    ``prepare_snapshot`` using the requested ``season``, ``ldyr``, and
    ``years_to_plot`` arguments.

    If a processed quiver component contains a ``lev`` dimension,
    ``select_depth_range`` is applied using ``lev_interp`` and ``lev_range``.
    The resulting field is then averaged over ``lev`` using ``mean("lev")``.
    The same vertical processing is applied to the optional background field
    when it contains a depth dimension.

    The depth-range description included in the subplot title is obtained from
    the first available quiver component for which
    ``select_depth_range(..., return_range_string=True)`` is called. If the
    X-component does not provide one, the corresponding range description from
    the Y-component is used.

    The optional ``scale_factor_X`` and ``scale_factor_Y`` values are applied
    after temporal and vertical preprocessing. These factors can be used to
    rescale the two vector components before plotting.

    If one quiver component is omitted, it is replaced using
    ``xr.zeros_like`` based on the provided component. For example, if
    ``quiver_dicts_X`` is ``None``, the X-component is set to zero everywhere
    using the processed Y-component as the template.

    Before plotting, the X- and Y-components are aligned using
    ``xr.align(..., join="inner")``. If a background field is provided, all three
    fields are aligned together. This restricts them to common latitude and
    longitude coordinates.

    When ``background_dict`` is supplied, the scalar field is shown with
    ``pcolormesh`` and a shared vertical colorbar is created. If no background
    field is supplied, only the quiver vectors are plotted and no colorbar is
    added.

    Quiver vectors are subsampled using ``quiver_step`` along both longitude and
    latitude. NaN values in either vector component are replaced with zero
    immediately before plotting.

    The temporal description used in each subplot title is obtained first from
    the X-component when available, otherwise from the Y-component. If no year
    string is returned by ``prepare_snapshot``, the function falls back to the
    ``y0`` and ``y1`` metadata of whichever quiver source dictionary is
    available and labels that interval as a mean.

    Documentation produced with the assistance of AI.
    """
    if not ds_list:
        raise ValueError(
            "'ds_list' cannot be empty."
        )

    if quiver_dicts_X is None and quiver_dicts_Y is None:
        raise ValueError(
            "At least one of 'quiver_dicts_X' or "
            "'quiver_dicts_Y' must be provided."
        )

    if quiver_step < 1:
        raise ValueError(
            "'quiver_step' must be at least 1."
        )

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(
        len(ds_list),
        1,
        figsize=figsize,
        squeeze=False,
    )

    im = None

    # ------------------------------------------------------------------
    # Datasets
    # ------------------------------------------------------------------
    for ds_idx, name in enumerate(ds_list):

        ax = axes[ds_idx, 0]

        # --------------------------------------------------------------
        # Quiver X component
        # --------------------------------------------------------------
        tsX = None
        year_string = None
        range_string = ""

        if quiver_dicts_X is not None:
            tsX, year_string = prepare_snapshot(
                quiver_dicts_X[biome][name].data,
                season=season,
                ldyr=ldyr,
                years_to_plot=years_to_plot,
                return_year_string=True,
            )

            if "lev" in tsX.dims:
                tsX, range_string = select_depth_range(
                    tsX,
                    lev_interp=lev_interp,
                    lev_range=lev_range,
                    return_range_string=True,
                )

                tsX = tsX.mean("lev")

            tsX = tsX * scale_factor_X

        # --------------------------------------------------------------
        # Quiver Y component
        # --------------------------------------------------------------
        tsY = None

        if quiver_dicts_Y is not None:
            tsY, year_string_Y = prepare_snapshot(
                quiver_dicts_Y[biome][name].data,
                season=season,
                ldyr=ldyr,
                years_to_plot=years_to_plot,
                return_year_string=True,
            )

            range_string_Y = ""

            if "lev" in tsY.dims:
                tsY, range_string_Y = select_depth_range(
                    tsY,
                    lev_interp=lev_interp,
                    lev_range=lev_range,
                    return_range_string=True,
                )

                tsY = tsY.mean("lev")

            tsY = tsY * scale_factor_Y

            if year_string is None:
                year_string = year_string_Y

            if not range_string:
                range_string = range_string_Y

        # --------------------------------------------------------------
        # Missing quiver component
        # --------------------------------------------------------------
        if tsX is None:
            tsX = xr.zeros_like(tsY)

        if tsY is None:
            tsY = xr.zeros_like(tsX)

        # --------------------------------------------------------------
        # Background field
        # --------------------------------------------------------------
        background_ts = None

        if background_dict is not None:
            background_ts = prepare_snapshot(
                background_dict[biome][name].data,
                season=season,
                ldyr=ldyr,
                years_to_plot=years_to_plot,
            )

            if "lev" in background_ts.dims:
                background_ts = select_depth_range(
                    background_ts,
                    lev_interp=lev_interp,
                    lev_range=lev_range,
                )

                background_ts = background_ts.mean("lev")

        # --------------------------------------------------------------
        # Align latitude / longitude coordinates
        # --------------------------------------------------------------
        if background_ts is not None:
            tsX, tsY, background_ts = xr.align(
                tsX,
                tsY,
                background_ts,
                join="inner",
            )

        else:
            tsX, tsY = xr.align(
                tsX,
                tsY,
                join="inner",
            )

        # --------------------------------------------------------------
        # Background field
        # --------------------------------------------------------------
        if background_ts is not None:
            im = ax.pcolormesh(
                background_ts.lon,
                background_ts.lat,
                background_ts,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
            )

        # --------------------------------------------------------------
        # Quiver
        # --------------------------------------------------------------
        ax.quiver(
            tsX.lon.values[::quiver_step],
            tsX.lat.values[::quiver_step],
            np.nan_to_num(
                tsX.values[::quiver_step, ::quiver_step],
                nan=0.0,
            ),
            np.nan_to_num(
                tsY.values[::quiver_step, ::quiver_step],
                nan=0.0,
            ),
            scale=arrow_scale,
            width=width,
            headwidth=headwidth,
            headlength=headlength,
            headaxislength=headaxislength,
        )

        # --------------------------------------------------------------
        # Year fallback
        # --------------------------------------------------------------
        if year_string is None:

            source_dict = (
                quiver_dicts_X
                if quiver_dicts_X is not None
                else quiver_dicts_Y
            )

            y0 = source_dict[biome][name].y0
            y1 = source_dict[biome][name].y1

            if y0 is not None and y1 is not None:
                year_string = f"{y0} - {y1} mean"
            else:
                year_string = ""

        # --------------------------------------------------------------
        # Title
        # --------------------------------------------------------------
        title_parts = [
            title,
            name,
            biome,
            year_string,
            f"{season} composite",
            range_string,
        ]

        title_ = " - ".join(
            part for part in title_parts if part
        )

        ax.set_title(
            title_,
            fontsize=font_size,
        )

        # --------------------------------------------------------------
        # Axis formatting
        # --------------------------------------------------------------
        ax.set_ylabel(
            "Lat ($^o$ North)",
            fontsize=font_size,
        )

        if ds_idx == len(ds_list) - 1:
            ax.set_xlabel(
                "Lon ($^o$ East)",
                fontsize=font_size,
            )
        else:
            ax.set_xlabel("")

        ax.tick_params(
            axis="both",
            which="major",
            labelsize=font_size,
        )

    # ------------------------------------------------------------------
    # Shared colorbar
    # ------------------------------------------------------------------
    if im is not None:
        cbar = fig.colorbar(
            im,
            ax=axes[:, 0],
            orientation="vertical",
        )

        if colorbar_label is not None:
            cbar.set_label(
                colorbar_label,
                fontsize=font_size,
            )

        cbar.ax.tick_params(
            labelsize=font_size,
        )

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    if save:
        if dir_name is None or file_name is None:
            raise ValueError(
                "'dir_name' and 'file_name' are required "
                "when save=True."
            )

        Path(dir_name).mkdir(
            parents=True,
            exist_ok=True,
        )

        plt.savefig(
            f"{dir_name}/{file_name}.png",
            bbox_inches="tight",
            dpi=300,
        )

    if return_fig_handles:
        return fig, axes





