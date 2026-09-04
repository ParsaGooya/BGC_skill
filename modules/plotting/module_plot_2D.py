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

from modules.analysis.module_global_averages import area_weighted_avg
from modules.analysis.module_data_postprocessing import (spatial_mask, 
                                                         Metrics, 
                                                         calculate_measure)
from modules.data_info.module_state_dict import state_dict
from modules.plotting.utils import *

    

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
    Prepare a seasonal snapshot.

    Parameters
    ----------
    data
        Input data with dimensions including some combination of
        year, month, lev, and lon.
    season
        Season or month to average.
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

        range_string = f"{lev_min:2}-{lev_max:2} m"

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





