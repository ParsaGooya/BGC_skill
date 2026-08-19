import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import ticker
import cartopy.crs as ccrs
from pathlib import Path
import cartopy
import cartopy.crs as ccrs
from mpl_toolkits.axes_grid1 import make_axes_locatable
import xarray as xr
from typing import Sequence

from modules.analysis.module_global_averages import area_weighted_avg
from modules.analysis.module_data_postprocessing import (spatial_mask, 
                                                         Metrics, 
                                                         calculate_measure)
from modules.data_info.module_state_dict import state_dict
from modules.plotting.utils import add_cyclic_point

MONTH_NAMES = (
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
)

SHIFTED_SEASONS = {"DJF", "MAM", "JJA", "SON"}


def get_season_indices(
    season: str,
    ldyr_ini: int = 0,
    ldyr_end: int = 1,
) -> np.ndarray:
    """
    Return monthly time indices for a requested season.

    JFM/AMJ/JAS/OND use the original January-centered temporal structure.

    DJF/MAM/JJA/SON use the same indices after the data have been shifted
    to a December-centered temporal structure.
    """
    if ldyr_end <= ldyr_ini:
        raise ValueError("'ldyr_end' must be greater than 'ldyr_ini'.")

    if season in MONTH_NAMES:
        month_idx = MONTH_NAMES.index(season)

        return np.array([
            year * 12 + month_idx
            for year in range(ldyr_ini, ldyr_end)
        ])

    seasons = {
        "JFM": (0, 1, 2),
        "AMJ": (3, 4, 5),
        "JAS": (6, 7, 8),
        "OND": (9, 10, 11),
        "DJF": (0, 1, 2),
        "MAM": (3, 4, 5),
        "JJA": (6, 7, 8),
        "SON": (9, 10, 11),
    }

    if season in seasons:
        return np.array([
            year * 12 + month_idx
            for year in range(ldyr_ini, ldyr_end)
            for month_idx in seasons[season]
        ])

    if season == "ANN":
        return np.arange(ldyr_ini * 12, ldyr_end * 12)

    raise ValueError(f"Unsupported season: {season!r}")


def shift_to_december_centered(
    ds: xr.DataArray,
) -> xr.DataArray:
    """
    Shift monthly data by one month so that time=0 corresponds to December
    of the previous year.

    For example, before shifting:

        year=2000, time=0 -> Jan 2000
        year=2000, time=1 -> Feb 2000

    After shifting:

        year=2000, time=0 -> Dec 1999
        year=2000, time=1 -> Jan 2000
        year=2000, time=2 -> Feb 2000

    The original dimension names are preserved.
    """
    if "year" not in ds.dims or "month" not in ds.dims:
        raise ValueError(
            "December-centered shifting requires both 'year' and 'time' dimensions."
        )

    stacked = ds.stack(_year_time=("year", "month"))

    shifted = stacked.shift(_year_time=1)

    return (
        shifted
        .unstack("_year_time")
        .transpose(*ds.dims)
    )


def seasonal_mean(
    ds: xr.DataArray,
    season: str,
    ldyr_ini: int = 0,
    ldyr_end: int = 1,
) -> xr.DataArray:
    """
    Select a season and average over its monthly time indices.

    DJF/MAM/JJA/SON are evaluated after shifting the data to the
    December-centered temporal convention.
    """
    if "month" not in ds.dims:
        return ds

    if season in SHIFTED_SEASONS:
        ds = shift_to_december_centered(ds)

    month_indices = get_season_indices(
        season=season,
        ldyr_ini=ldyr_ini,
        ldyr_end=ldyr_end,
    )

    return ds.isel(month=month_indices).mean("month")


def resolve_depth_range(ds: xr.DataArray | xr.Dataset, depth_range: float | Sequence):

    if isinstance(depth_range, float):
        ds = ds.sel(lev = depth_range, method = 'nearest') 
        selected_depth = np.round(ds.lev.values,2)
        
        return ds, selected_depth

    if len(depth_range) > 2:
        depth_range = slice(np.floor(depth_range.min()), np.ceil(depth_range.max()))

    else:
        depth_range = slice(depth_range[0], depth_range[1])
    
    ds = ds.sel(lev = depth_range)
    selected_depth = (np.round(ds.lev.values.min(),2), np.round(ds.lev.values.max(),2))

    return ds, selected_depth
    

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


def plot_composites(
    ds_list: list[str],
    data_dict: dict[str, state_dict],
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

                if "year" in obs_aligned.dims:
                    y0 = obs_aligned.year.values[0]
                    y1 = obs_aligned.year.values[-1]

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
    ds_list: list[str],
    data_dict: dict[str, state_dict],
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



def plot_maps(ds,
              mask = None,
              lat=None,
              lon=None,
              central_longitude=180,
              gridlines=False,
              cmap='RdBu_r',
              vmin=0,
              vmax=1,
              nvals=12,
              cbar=False,
              cbar_label='',
              ncols=3,
              titles=None,
              figsize=(12,3),
              fig_dir=None,
              fig_name=None,
              show=False,
              save=False,
              **kwargs):
    fnt_size = 14
    mpl.rcParams.update({'font.size': fnt_size})
#   hfont = {'fontname':'Calibri'}

    if lat is None:
        lat = ds.lat
    if lon is None:
        lon = ds.lon
        
    img_extent = [lon[0], lon[-1], lat[0], lat[-1]]
    fig, ax = plt.subplots(nrows=1,
                           ncols=len(ds.time),
                           figsize=figsize, 
                           subplot_kw={'projection': ccrs.PlateCarree(central_longitude=central_longitude)})

    scale = (vmax-vmin)/float(nvals)
    vals = vmin + (vmax-vmin)*np.arange(nvals+1)/float(nvals)
    # norm = mpl.colors.BoundaryNorm(vals, cmap.N)
    norm = mpl.colors.BoundaryNorm(vals, plt.cm.get_cmap(cmap).N)
    
    for i, (axis, ds_lead) in enumerate(zip(ax, ds)):
        if gridlines:
            axis.gridlines(draw_labels=False)
        im = axis.imshow(ds_lead, 
                         origin='lower',
                         extent=img_extent,
                         # cmap=cmap,
                         cmap=plt.cm.get_cmap(cmap),
                         norm=norm,
                         transform=ccrs.PlateCarree())
        im.set_clim(vmin, vmax)
        axis.coastlines()
        if mask is None:
            mask = spatial_mask(ds_lead.to_dataset())
        if titles:
            # area_txt = str(np.round_(area_weighted_avg(ds_lead,
            #                                           is_ds=False),2).values)
            area_txt = str(np.round(area_weighted_avg(ds_lead,
                                                      mask=mask,
                                                      is_ds=False),2).values)
                                                      
                                                      
            axis.set_title(f'{titles}, Yr {ds_lead.time.values + 1}, avg={area_txt}', fontsize=fnt_size)
    
    if cbar:
        if ncols == 2:
            clb_x = 0.2 
            clb_y = -0.1
            clb_w = 0.6
            clb_h = 0.07
        if ncols == 3:
            clb_x = 0.2 
            clb_y = 0.05 
            clb_w = 0.6
            clb_h = 0.05 
        cax = plt.axes([clb_x, # left
                        clb_y, # bottom
                        clb_w, # width
                        clb_h])# height
        cb = mpl.colorbar.ColorbarBase(ax=cax,
                                       cmap=cmap,
                                       norm=norm,
                                       spacing='uniform',
                                       orientation='horizontal',
                                       extend='both',
                                       ticks=vals)
        tick_locator = ticker.MaxNLocator(nbins=nvals/2)
        cb.locator = tick_locator
        cb.update_ticks()

        cb.set_label(label=cbar_label) #, **hfont) 
        
    fig.tight_layout()
    
    if save:
        Path(fig_dir).mkdir(parents=True,
                            exist_ok=True)
        plt.savefig(f'{fig_dir}/{fig_name}',
                    bbox_inches='tight',
                    dpi=300)
    
    if show is False:
        plt.close()


def plot_depth_vs_time_biomeavg(ds_list,
                            ds_dicts,
                            bms_label,
            xx=None,
            ldyr=1-1,
            title='',
            figsize=(10,45),
            dict_label='ts',
            contourf_levels = [0,0.5, 1,2,4,6, 8,10,15,20,25,30,35,40,50],
            cmap = 'viridis',
            dir_name=None,
            file_name=None,
            colorbar_label = None, 
            season = 'ANN',
            lev_interp = None,
            lev_range = None,
            monthly_res = False,
            ENSO_years = True,
            show=False,
            return_fig_handles = False,
            font_size = 14,
            save=False):
    '''
     to plot data written in terms of target years
     on a common period for each lead year
     
    '''


    for jj in range(4):
        sea = [ii*12+3*jj+np.arange(3) for ii in range(ldyr,
                                                        ldyr + 1 )]
        seas = list(np.stack(sea,
                            axis=0).flatten())
        if jj == 0:
            JFM = seas
            print(f"JFM:{seas}")
        if jj == 1:
            AMJ = seas
            print(f"AMJ:{seas}")
        if jj == 2:
            JAS = seas
            print(f"JAS:{seas}")
        if jj == 3:
            OND = seas
            print(f"OND:{seas}")
            print('======')
    print(f"ANN:{np.arange(12) + 12* ldyr}")
    print('======')
    seasons = {'JFM': JFM,
            'AMJ': AMJ,
            'JAS': JAS,
            'OND': OND,
            'ANN': np.arange(12) + 12* ldyr}

    for ind, month in enumerate(['Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']):
            seasons[month] = [ii*12+ ind for ii in range(ldyr, ldyr + 1 )]
    
    seasons['DJF'] = [0,1,2]
    seasons['MAM'] = [2,3,4]
    seasons['JJA'] = [5,6,7]
    seasons['SON'] = [8,9,10]


    if season in ['DJF','MAM','JJA','SON']:
        assert ldyr == 0
        for ds in ds_list:
            assert 'hindcast' not in ds

    if monthly_res:
        assert season == 'ANN', ldyr == 0

    if show:
        ldyr1 = ldyr + 1
    

        fig, ax = plt.subplots(len(ds_list), 1 ,figsize=figsize)
        ds_dict = ds_dicts[bms_label]

        for ind, ds in enumerate(ds_list):
            # if len(ds.split('-')) == 1:
            dstp = ds_dict[ds][dict_label]
            # else:
            #     ds1 = (ds.split('-')[0]).split(' ')[0]
            #     ds0 = (ds.split('-')[1]) .split(' ')[1]
            #     dstp = ds_dict[ds1][dict_label] - ds_dict[ds0][dict_label]

            if season == 'DJF':
                try:
                    ts = DJFy(dstp).sel(time=slice(ldyr*12,
                                                        (ldyr1)*12-1))
                except:
                    ts = dstp.sel(time=slice(ldyr*12,
                                                        (ldyr1)*12-1))                    
                    print('DJFy unseccessfull, make sure the conversion is already done!')
            else:
                ts = dstp.sel(time=slice(ldyr*12,
                                                        (ldyr1)*12-1))
            ts = ts.sel(time = seasons[season] )
            if lev_interp is not None:
                ts = ts.interp(lev = lev_interp )
            
            if lev_range is not None:
                if type(lev_range) is not list:
                    ls = [0, lev_range]
                    lev_range = ls
                    del ls
            
                ts = ts.where((ts.lev <= lev_range[1]) & (ts.lev >= lev_range[0]), drop = True)

            if monthly_res:
                ts = ts.stack(yearmonth = ('year','time'))
                xx = ts.year.values + (ts.time.values + 0.5)/12
                ts['yearmonth'] = ts.year.values + (ts.time.values + 0.5)/12
                dim = 'yearmonth'
            else:
                ts = ts.mean('time')
                xx = ts.year.values     
                dim = 'year' 

            label = f'{ds}'

            if len(ds_list) == 1:
                ax_ = ax
            else:
                ax_ = ax[ind]
            
            contour_f = ax_.contourf(xx, ts.lev.values,
                        ts,
                        levels = contourf_levels,
                        cmap = cmap
                        # vmax = 30,
                        # vmin = 0,
                        # ds_dict[ds]['linestyle'],
                            )
            contours = ax_.contour(xx, ts.lev.values, ts, colors='white', levels=contour_f.levels)
            ax_.clabel(contours, inline=True, fontsize=font_size, colors='white')
            ax_.set_xticks(np.unique(np.floor(xx)), np.unique(np.floor(xx)).astype(int), rotation  = 45)
                        # color=ds_dict[ds]['color'])
            # ax[ind].colorbar()
                
            ax_.set_title(title + ' - ' +  ds +' - ' + bms_label, fontsize = font_size)
            ax_.set_ylabel('depth (m)')
            ax_.invert_yaxis()

            if ENSO_years is not None:
                if ENSO_years == 'obs':
                    if monthly_res:
                        ENSO_years = [1982 + 3/12, 1983 + 6/12, 1986 + 8/12, 1988 + 2/12, 1991 + 4/12, 1992 + 6/12, 1994 + 8/12, 1995 + 4/ 12, 1997 + 4/12, 1998 + 5/12 , 
                                      2002 + 5 /12, 2003 + 3/12,  2004 + 6/12, 2005 + 2/12, 2006 + 8/12, 2006 + 1/12, 2009 + 6/12, 2010 + 3/12, 2014 + 9/12, 2016 + 4/12,  
                                      2018 + 8/12, 2019 + 6/12]#, 2023 + 5/12, 2024  +5/12]
                    else:
                        ENSO_years = [1982.5, 1983.5, 1986.5,1987.5, 1996.5, 1999.5, 2008.5, 2010.5,2014.5, 2016.5, 2022.5, 2024.5]    


                for x in (ENSO_years):
                    if cmap == 'viridis':
                        ax_.axvline(x=x,linestyle = 'dashed', color = 'r', alpha = 0.5)
                    else:
                        ax_.axvline(x=x,linestyle = 'dashed', color = 'g', alpha = 0.5)

        cbar = fig.colorbar(contour_f , ax=ax, orientation='vertical', label = colorbar_label)
        cbar.ax.set_ylabel(colorbar_label ,fontsize=font_size)
        cbar.ax.tick_params(labelsize=font_size)

    if save:
        if show is False:
            print(f"If show = {show} then save = False --change to show = True")
        else:
            Path(dir_name).mkdir(parents=True,
                                     exist_ok=True)
            plt.savefig(f'{dir_name}/{file_name}.png',
                            bbox_inches='tight',
                            dpi=300)   
            

def snapshot_depth_vs_lon(ds_list,
                            ds_dicts,
                            bms_label,
                            yeartoplot = None,
                            contour_dict = None,
                            quiver_dicts = None, 
                            quiver_axis = None,
            ldyr=1-1,
            title='',
            figsize=(10,45),
            contourf_levels = [0,0.5, 1,2,4,6, 8,10,15,20,25,30,35,40,50],
            contour_var = None,
            contour_var_levels = None,
            cmap = 'viridis',
            dir_name=None,
            file_name=None,
            colorbar_label = None, 
            season = 'ANN',
            lev_interp = None,
            lev_range = None,
            show=False,
            return_fig_handles = False,
            headwidth = 5,
            headlength = 1, 
            headaxislength = 2,
            width =  0.005,
            font_size = 14,
            save=False):
    '''
     to plot data written in terms of target years
     on a common period for each lead year
     
    '''
        
    for jj in range(4):
        sea = [ii*12+3*jj+np.arange(3) for ii in range(ldyr,
                                                        ldyr + 1 )]
        seas = list(np.stack(sea,
                            axis=0).flatten())
        if jj == 0:
            JFM = seas
            print(f"JFM:{seas}")
        if jj == 1:
            AMJ = seas
            print(f"AMJ:{seas}")
        if jj == 2:
            JAS = seas
            print(f"JAS:{seas}")
        if jj == 3:
            OND = seas
            print(f"OND:{seas}")
            print('======')
    print(f"ANN:{np.arange(12) + 12* ldyr}")
    print('======')
    seasons = {'JFM': JFM,
            'AMJ': AMJ,
            'JAS': JAS,
            'OND': OND,
            'ANN': np.arange(12) + 12* ldyr}

    for ind, month in enumerate(['Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']):
            seasons[month] = [ii*12+ ind for ii in range(ldyr, ldyr + 1 )]
    
    seasons['DJF'] = [0,1,2]
    seasons['MAM'] = [2,3,4]
    seasons['JJA'] = [5,6,7]
    seasons['SON'] = [8,9,10]
    seasons['ENSO_diff'] = [-1]

    if season in ['DJF','MAM','JJA','SON']:
        assert ldyr == 0
        for ds in ds_list:
            assert 'hindcast' not in ds
    if season not in seasons.keys():
        season = 'ENSO_diff'
        raise Warning('Season was set to ENSO_diff !')
    if season == 'ENSO_diff':
        yeartoplot = ['diff']

    if show:
        ldyr1 = ldyr + 1
    

        fig, ax = plt.subplots(len(ds_list), 1 ,figsize=figsize)
        ds_dict = ds_dicts[bms_label]

        if contour_dict is not None:
            ds_dict_contour = contour_dict[bms_label]

        if quiver_dicts is not None:
            qv_dict = quiver_dicts[bms_label]

        for ind, ds in enumerate(ds_list):
            # if len(ds.split('-')) == 1:
            dstp = ds_dict[ds]

            if quiver_dicts is not None:
                qvtp = qv_dict[ds]
            if contour_dict is not None:
                dstp_contour = ds_dict_contour[ds]
            # else:
            #     ds1 = (ds.split('-')[0]).split(' ')[0]
            #     ds0 = (ds.split('-')[1]) .split(' ')[1]
            #     dstp = ds_dict[ds1] - ds_dict[ds0]

            #     if quiver_dicts is not None:
            #         qvtp = qv_dict[ds1] - qv_dict[ds0]

            if season == 'DJF':
                try:
                    ts = DJFy(dstp).sel(time=slice(ldyr*12,(ldyr1)*12-1))
                                              
                    if quiver_dicts is not None:
                        qv_ts = DJFy(qvtp).sel(time=slice(ldyr*12,(ldyr1)*12-1))
                    if contour_dict is not None:
                        contour_ts = DJFy(dstp_contour).sel(time=slice(ldyr*12,(ldyr1)*12-1))
                            
                except:
                    ts = dstp.sel(time=slice(ldyr*12,(ldyr1)*12-1))

                    if quiver_dicts is not None:
                        qv_ts = qvtp.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                    if contour_dict is not None:
                        contour_ts = dstp_contour.sel(time=slice(ldyr*12,(ldyr1)*12-1))  
                                                        
                    print('DJFy unseccessfull, make sure the conversion is already done!')
            else:
                
                    ts = dstp.sel(time=slice(ldyr*12, (ldyr1)*12-1)) if season not in ['ENSO_diff'] else dstp
                    if quiver_dicts is not None:
                        qv_ts = qvtp.sel(time=slice(ldyr*12, (ldyr1)*12-1))  if season not in ['ENSO_diff'] else qvtp
                    if contour_dict is not None:
                        contour_ts = dstp_contour.sel(time=slice(ldyr*12,(ldyr1)*12-1))  if season not in ['ENSO_diff'] else dstp_contour                               
                                                        
            if yeartoplot is not None:
                if all([type(yeartoplot) == str, 'diff' not in yeartoplot]): 
                
                    y1 = eval(yeartoplot.split('-')[0])
                    y0 = eval(yeartoplot.split('-')[1])
                    ts = (ts.sel(time = seasons[season] ).sel(year = y1 ).mean('time') - ts.sel(time = seasons[season] ).sel(year = y0).mean('time')).squeeze()

                    if quiver_dicts is not None:
                        qv_ts = (qv_ts.sel(time = seasons[season] ).sel(year = y1 ).mean('time') - qv_ts.sel(time = seasons[season] ).sel(year = y0).mean('time')).squeeze()
                    if contour_dict is not None:
                        contour_ts = (contour_ts.sel(time = seasons[season] ).sel(year = y1 ).mean('time') - contour_ts.sel(time = seasons[season] ).sel(year = y0).mean('time')).squeeze()

                else:
                    ts = ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean(['year','time']).squeeze()

                    if quiver_dicts is not None:
                        qv_ts = qv_ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean(['year','time']).squeeze()
                    if contour_dict is not None:
                        contour_ts = contour_ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean(['year','time']).squeeze()
            else:
                ts = ts.sel(time = seasons[season] ).mean(['year','time']).squeeze()

                if quiver_dicts is not None:
                    qv_ts = qv_ts.sel(time = seasons[season] ).mean(['year','time']).squeeze()
                if contour_dict is not None:
                    contour_ts = contour_ts.sel(time = seasons[season] ).mean(['year','time']).squeeze()

            xx = ts.lon.values
            
            if lev_interp is not None:
                ts = ts.interp(lev = lev_interp )
                if quiver_dicts is not None:
                    qv_ts = qv_ts.interp(lev = lev_interp )
                if contour_dict is not None:
                    contour_ts = contour_ts.interp(lev = lev_interp )
                        
            if lev_range is not None:
                if type(lev_range) is not list:
                    ls = [0, lev_range]
                    lev_range = ls
                    del ls
            
                ts = ts.where((ts.lev <= lev_range[1]) & (ts.lev >= lev_range[0]), drop = True)
                if quiver_dicts is not None:
                    qv_ts = qv_ts.where((qv_ts.lev <= lev_range[1]) & (qv_ts.lev >= lev_range[0]), drop = True)
                if contour_dict is not None:
                    contour_ts = contour_ts.where((contour_ts.lev <= lev_range[1]) &  (contour_ts.lev >= lev_range[0]), drop = True)
                

            label = f'{ds}'
            if len(ds_list) == 1:
                ax_ = ax
            else:
                ax_ = ax[ind]
            
            contour_f = ax_.contourf(xx, ts.lev.values,
                        ts,
                        levels = contourf_levels,
                        cmap = cmap
                        # vmax = 30,
                        # vmin = 0,
                        # ds_dict[ds]['linestyle'],
                            )
            if contour_dict is not None:
                contours = ax_.contour(xx, contour_ts.lev.values, contour_ts, colors='black', levels = contour_var_levels, alpha = 0.5) 
                ax_.clabel(contours, inline=True, fontsize=font_size, colors='black')
        
            else:
                contours = ax_.contour(xx, ts.lev.values, ts, colors='white', levels=contour_f.levels)
                ax_.clabel(contours, inline=True, fontsize=font_size, colors='white')

            if quiver_dicts is not None:
                if quiver_axis == 'U':
                    U = qv_ts.values
                    V = np.zeros_like(qv_ts.values)
                else:
                    V = qv_ts.values
                    U = np.zeros_like(qv_ts.values)
                ax_.quiver(xx[::2], qv_ts.lev.values[::2], np.nan_to_num(U[::2,::2], nan=0.0), np.nan_to_num(V[::2,::2], nan=0.0),  alpha = 0.5,width = width, headwidth = headwidth,headlength = headlength, headaxislength = headaxislength)

            title_ = title + ' - ' +  ds +' - ' + bms_label
            # if  type(yeartoplot) == str:
            #     title_ = title_ + f' {yeartoplot}' + f' {season}'         
            # elif all([len(yeartoplot) >1 , type(yeartoplot) == list]):
            title_ = title_ + f' composite' + f' {season}'               
            # else:
            #     if yeartoplot[0] not in states_dict.keys():
            #         states_dict[yeartoplot[0]] = 'NA'
            #     title_ = title_ + f' {yeartoplot[0]}' + f' {season}' + f': ENSO {states_dict[yeartoplot[0]]}'
            
            if contour_dict is not None:
                title_ = title_ + f' - {contour_var} contours'            
            
            
            ax_.set_title( title_, fontsize = font_size)
            ax_.set_ylabel('depth (m)',fontsize = font_size)
            ax_.set_xlabel('Lon ($^o$ East)', fontsize = font_size)
            if ind < len(ds_list) - 1:
                ax_.set_xlabel('', fontsize = font_size)
            ax_.invert_yaxis()
            ax_.tick_params(axis='x', which='major', labelsize=font_size)
            ax_.tick_params(axis='y', which='major', labelsize=font_size)


        cbar = fig.colorbar(contour_f , ax=ax, orientation='vertical', label = colorbar_label)
        cbar.ax.set_ylabel(colorbar_label ,fontsize=font_size)
        cbar.ax.tick_params(labelsize=font_size)
    if save:
        if show is False:
            print(f"If show = {show} then save = False --change to show = True")
        else:
            Path(dir_name).mkdir(parents=True,
                                     exist_ok=True)
            plt.savefig(f'{dir_name}/{file_name}.png',
                            bbox_inches='tight',
                            dpi=300)   

def snapshot_arial(ds_list,
                            ds_dicts,
                            bms_label,
                            yeartoplot = None, 
                            quiver_dicts = None, 
            ldyr=1-1,
            title='',
            figsize=(10,45),
            vmax=None,
            vmin=None,
            cmap = 'viridis',
            dir_name=None,
            file_name=None,
            colorbar_label = None, 
            season = 'ANN',
            lev_interp = None,
            lev_range:list = None,
            integrate = False,
            show=False,
            return_fig_handles = False,
            font_size = 14,
            save=False):
    '''
     to plot data written in terms of target years
     on a common period for each lead year
     
    '''
        
    for jj in range(4):
        sea = [ii*12+3*jj+np.arange(3) for ii in range(ldyr,
                                                        ldyr + 1 )]
        seas = list(np.stack(sea,
                            axis=0).flatten())
        if jj == 0:
            JFM = seas
            print(f"JFM:{seas}")
        if jj == 1:
            AMJ = seas
            print(f"AMJ:{seas}")
        if jj == 2:
            JAS = seas
            print(f"JAS:{seas}")
        if jj == 3:
            OND = seas
            print(f"OND:{seas}")
            print('======')
    print(f"ANN:{np.arange(12) + 12* ldyr}")
    print('======')
    seasons = {'JFM': JFM,
            'AMJ': AMJ,
            'JAS': JAS,
            'OND': OND,
            'ANN': np.arange(12) + 12* ldyr}

    for ind, month in enumerate(['Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']):
            seasons[month] = [ii*12+ ind for ii in range(ldyr, ldyr + 1 )]
    
    seasons['DJF'] = [0,1,2]
    seasons['MAM'] = [2,3,4]
    seasons['JJA'] = [5,6,7]
    seasons['SON'] = [8,9,10]
    seasons['ENSO_diff'] = [-1]
    
    if season in ['DJF','MAM','JJA','SON']:
        assert ldyr == 0
        for ds in ds_list:
            assert 'hindcast' not in ds
    if season not in seasons.keys():
        season = 'ENSO_diff'
        raise Warning('Season was set to ENSO_diff !')
    if season == 'ENSO_diff':
        yeartoplot = ['diff']

    if show:
        ldyr1 = ldyr + 1
    

        fig, ax = plt.subplots(len(ds_list), 1 ,figsize=figsize)
        ds_dict = ds_dicts[bms_label]

        if quiver_dicts is not None:
            qv_dict = quiver_dicts[bms_label]

        for ind, ds in enumerate(ds_list):

            dstp = ds_dict[ds]
            dstp_lev = True if 'lev' in dstp.dims else False 
            if quiver_dicts is not None:
                qvtp = qv_dict[ds]
                qvtp_lev = True if 'lev' in qvtp.dims else False 

            if season == 'DJF':
                try:
                    ts = DJFy(dstp).sel(time=slice(ldyr*12,(ldyr1)*12-1))
                                                        
                    if quiver_dicts is not None:
                        qv_ts = DJFy(qvtp).sel(time=slice(ldyr*12,(ldyr1)*12-1))

                except:
                    ts = dstp.sel(time=slice(ldyr*12, (ldyr1)*12-1))
                                                       
                    if quiver_dicts is not None:
                        qv_ts = qvtp.sel(time=slice(ldyr*12, (ldyr1)*12-1))

                    print('DJFy unseccessfull, make sure the conversion is already done!')
            else:
                
                    ts = dstp.sel(time=slice(ldyr*12, (ldyr1)*12-1))  if season not in ['ENSO_diff'] else dstp

                    if quiver_dicts is not None:
                        qv_ts = qvtp.sel(time=slice(ldyr*12, (ldyr1)*12-1))  if season not in ['ENSO_diff'] else qvtp
                                                       
            if yeartoplot is not None:
                if all([type(yeartoplot) == str, 'diff' not in yeartoplot]): 
                    y1 = eval(yeartoplot.split('-')[0])
                    y0 = eval(yeartoplot.split('-')[1])
                    ts = (ts.sel(time = seasons[season] ).sel(year = y1 ).mean('time') - ts.sel(time = seasons[season] ).sel(year = y0).mean('time')).squeeze()

                    if quiver_dicts is not None:
                        qv_ts = (qv_ts.sel(time = seasons[season] ).sel(year = y1 ).mean('time') - qv_ts.sel(time = seasons[season] ).sel(year = y0).mean('time')).squeeze()

                else:
                    ts = ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean(['year','time']).squeeze()

                    if quiver_dicts is not None:
                        qv_ts = qv_ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean(['year','time']).squeeze()
            else:
                    ts = ts.sel(time = seasons[season] ).mean(['year','time']).squeeze()

                    if quiver_dicts is not None:
                        qv_ts = qv_ts.sel(time = seasons[season] ).mean(['year','time']).squeeze()   

            xx = ts.lon.values
            
            if lev_interp is not None:
                ts = ts.interp(lev = lev_interp ) if dstp_lev else ts
                if quiver_dicts is not None:
                    qv_ts = qv_ts.interp(lev = lev_interp ) if qvtp_lev else qv_ts
            
            if lev_range is not None:
                ts = ts.where(ts.lev >= lev_range[0], drop = True)  if dstp_lev else ts
                ts = ts.where(ts.lev <= lev_range[1], drop = True)  if dstp_lev else ts
                if quiver_dicts is not None:
                    qv_ts = qv_ts.where(qv_ts.lev >= lev_range[0], drop = True)  if qvtp_lev else qv_ts
                    qv_ts = qv_ts.where(qv_ts.lev <= lev_range[1], drop = True)  if qvtp_lev else qv_ts

            if integrate:
                ts = ts.integrate('lev')  if dstp_lev else ts
            else: 
                ts = ts.mean('lev')  if dstp_lev else ts
            if quiver_dicts is not None:
                qv_ts = ts.mean('lev')  if qvtp_lev else qv_ts


            label = f'{ds}'
            if len(ds_list) == 1:
                ax_ = ax
            else:
                ax_ = ax[ind]
            
            # crs = ccrs.PlateCarree(central_longitude=central_longitude)    
            im = ax_.pcolormesh(ts.lon, ts.lat, ts, 
                            cmap=cmap, vmin = vmin, vmax = vmax)
            # ax_.invert_yaxis()

            # if quiver_dicts is not None:
                ##### add #####
            
            # if  type(yeartoplot) == str:
            #     ax_.set_title(title + ' - ' + f'{lev_range[0]}-{lev_range[1]} (m)' + ' - ' +  ds +' - ' + bms_label + f' {yeartoplot}' + f' {season}', fontsize = font_size)
            
            # elif all([len(yeartoplot) >1 , type(yeartoplot) == list]):
            ax_.set_title(title + ' - ' + f'{lev_range[0]}-{lev_range[1]} (m)' + ' - ' +  ds +' - ' + bms_label + f' composite' + f' {season}', fontsize = font_size)
            # else:
            #     if yeartoplot[0] not in states_dict.keys():
            #         states_dict[yeartoplot[0]] = 'NA'
            #     ax_.set_title(title + ' - ' + f'{lev_range[0]}-{lev_range[1]} (m)' + ' - ' +  ds +' - ' + bms_label + f' {yeartoplot[0]}' + f' {season}' + f': ENSO {states_dict[yeartoplot[0]]}', fontsize = font_size)
            ax_.set_ylabel('lat ($^o$ North)', fontsize = font_size)
            ax_.set_xlabel('Lon ($^o$ East)', fontsize = font_size)
            if ind < len(ds_list) - 1:
                ax_.set_xlabel('', fontsize = font_size)
            ax_.tick_params(axis='x', which='major', labelsize=font_size)
            ax_.tick_params(axis='y', which='major', labelsize=font_size)
            # ax_.invert_yaxis()


        cbar = fig.colorbar(im , ax=ax, orientation='vertical', label = colorbar_label)
        cbar.ax.set_ylabel(colorbar_label ,fontsize=font_size)
        cbar.ax.tick_params(labelsize=font_size)
    if save:
        if show is False:
            print(f"If show = {show} then save = False --change to show = True")
        else:
            Path(dir_name).mkdir(parents=True,
                                     exist_ok=True)
            plt.savefig(f'{dir_name}/{file_name}.png',
                            bbox_inches='tight',
                            dpi=300)   
            

def snapshot_crossection(ds_list,
                            ds_dicts,
                            bms_label,
                            longitude,
                            yeartoplot = None,
                            contour_dict = None,
                            quiver_dicts = None, 
                            quiver_axis = 'U',
            ldyr=1-1,
            title='',
            figsize=(10,45),
            contourf_levels = [0,0.5, 1,2,4,6, 8,10,15,20,25,30,35,40,50],
            contour_var = None,
            contour_var_levels = None,
            cmap = 'viridis',
            dir_name=None,
            file_name=None,
            colorbar_label = None, 
            season = 'ANN',
            lev_interp = None,
            lev_range:list = None,
            show=False,
            return_fig_handles = False,
            width =  0.005,
            headwidth = 5,
            headlength = 1, 
            headaxislength = 2,
            font_size = 14,
            save=False):
    '''
     to plot data written in terms of target years
     on a common period for each lead year
     
    '''
        
    for jj in range(4):
        sea = [ii*12+3*jj+np.arange(3) for ii in range(ldyr,
                                                        ldyr + 1 )]
        seas = list(np.stack(sea,
                            axis=0).flatten())
        if jj == 0:
            JFM = seas
            print(f"JFM:{seas}")
        if jj == 1:
            AMJ = seas
            print(f"AMJ:{seas}")
        if jj == 2:
            JAS = seas
            print(f"JAS:{seas}")
        if jj == 3:
            OND = seas
            print(f"OND:{seas}")
            print('======')
    print(f"ANN:{np.arange(12) + 12* ldyr}")
    print('======')
    seasons = {'JFM': JFM,
            'AMJ': AMJ,
            'JAS': JAS,
            'OND': OND,
            'ANN': np.arange(12) + 12* ldyr}

    for ind, month in enumerate(['Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']):
            seasons[month] = [ii*12+ ind for ii in range(ldyr, ldyr + 1 )]
    
    seasons['DJF'] = [0,1,2]
    seasons['MAM'] = [2,3,4]
    seasons['JJA'] = [5,6,7]
    seasons['SON'] = [8,9,10]
    seasons['ENSO_diff'] = [-1]

    if season in ['DJF','MAM','JJA','SON']:
        assert ldyr == 0
        for ds in ds_list:
            assert 'hindcast' not in ds
    if season not in seasons.keys():
        season = 'ENSO_diff'
        raise Warning('Season was set to ENSO_diff !')
    if season == 'ENSO_diff':
        yeartoplot = ['diff']
    if show:
        ldyr1 = ldyr + 1
    

        fig, ax = plt.subplots(len(ds_list), 1 ,figsize=figsize)
        ds_dict = ds_dicts[bms_label]

        if contour_dict is not None:
            ds_dict_contour = contour_dict[bms_label]

        if quiver_dicts is not None:
            qv_dict = quiver_dicts[bms_label]

        for ind, ds in enumerate(ds_list):
            # if len(ds.split('-')) == 1:
            dstp = ds_dict[ds]
            if quiver_dicts is not None:
                qvtp = qv_dict[ds]
            if contour_dict is not None:
                dstp_contour = ds_dict_contour[ds]
            # else:
            #     ds1 = (ds.split('-')[0]).split(' ')[0]
            #     ds0 = (ds.split('-')[1]) .split(' ')[1]
            #     dstp = ds_dict[ds1] - ds_dict[ds0]
            #     if quiver_dicts is not None:
            #         qvtp = qv_dict[ds1] - qv_dict[ds0]

            if season == 'DJF':
                try:
                    ts = DJFy(dstp).sel(time=slice(ldyr*12,(ldyr1)*12-1))
                    if quiver_dicts is not None:
                        qv_ts = DJFy(qvtp).sel(time=slice(ldyr*12,(ldyr1)*12-1))
                    if contour_dict is not None:
                        contour_ts = DJFy(dstp_contour).sel(time=slice(ldyr*12,(ldyr1)*12-1))
                                                        
                except:
                    ts = dstp.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                    if quiver_dicts is not None:
                        qv_ts = qvtp.sel(time=slice(ldyr*12,(ldyr1)*12-1))  
                    if contour_dict is not None:
                        contour_ts = dstp_contour.sel(time=slice(ldyr*12,(ldyr1)*12-1))  

                    print('DJFy unseccessfull, make sure the conversion is already done!')
            else:

                    ts = dstp.sel(time=slice(ldyr*12,(ldyr1)*12-1))  if season not in ['ENSO_diff'] else dstp
                    if quiver_dicts is not None:
                        qv_ts = qvtp.sel(time=slice(ldyr*12,(ldyr1)*12-1))  if season not in ['ENSO_diff'] else qvtp
                    if contour_dict is not None:
                        contour_ts = dstp_contour.sel(time=slice(ldyr*12,(ldyr1)*12-1))  if season not in ['ENSO_diff'] else dstp_contour                               
            if yeartoplot is not None:
                if all([type(yeartoplot) == str, 'diff' not in yeartoplot]): 
                    y1 = eval(yeartoplot.split('-')[0])
                    y0 = eval(yeartoplot.split('-')[1])
                    ts = (ts.sel(time = seasons[season] ).sel(year = y1 ).mean('time') - ts.sel(time = seasons[season] ).sel(year = y0).mean('time')).squeeze()
                    if quiver_dicts is not None:
                        qv_ts = (qv_ts.sel(time = seasons[season] ).sel(year = y1 ).mean('time') - qv_ts.sel(time = seasons[season] ).sel(year = y0).mean('time')).squeeze()
                    if contour_dict is not None:
                        contour_ts = (contour_ts.sel(time = seasons[season] ).sel(year = y1 ).mean('time') - contour_ts.sel(time = seasons[season] ).sel(year = y0).mean('time')).squeeze()

                else:
                    ts = ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean(['year','time']).squeeze()
                    if quiver_dicts is not None:
                        qv_ts = qv_ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean(['year','time']).squeeze()
                    if contour_dict is not None:
                        contour_ts = contour_ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean(['year','time']).squeeze()
            else:
                ts = ts.sel(time = seasons[season] ).mean(['year','time']).squeeze()
                if quiver_dicts is not None:
                    qv_ts = qv_ts.sel(time = seasons[season] ).mean(['year','time']).squeeze()
                if contour_dict is not None:
                    contour_ts = contour_ts.sel(time = seasons[season] ).mean(['year','time']).squeeze()             
            xx = ts.lat.values
            
            if lev_interp is not None:
                ts = ts.interp(lev = lev_interp )
                if quiver_dicts is not None:
                    qv_ts = qv_ts.interp(lev = lev_interp )
                if contour_dict is not None:
                    contour_ts = contour_ts.interp(lev = lev_interp )
            
            if lev_range is not None:
                if type(lev_range) is not list:
                    ls = [0, lev_range]
                    lev_range = ls
                    del ls                
                ts = ts.where((ts.lev <= lev_range[1]) & (ts.lev >= lev_range[0]), drop = True)
                if quiver_dicts is not None:
                    qv_ts = qv_ts.where((qv_ts.lev <= lev_range[1]) &  (qv_ts.lev >= lev_range[0]), drop = True)
                if contour_dict is not None:
                    contour_ts = contour_ts.where((contour_ts.lev <= lev_range[1]) &  (contour_ts.lev >= lev_range[0]), drop = True)
                
            ts = ts.sel(lon = longitude, method = 'nearest')
            if quiver_dicts is not None:
                qv_ts = qv_ts.sel(lon = longitude, method = 'nearest')
            if contour_dict is not None:
                contour_ts = contour_ts.sel(lon = longitude, method = 'nearest')

            label = f'{ds}'
            if len(ds_list) == 1:
                ax_ = ax
            else:
                ax_ = ax[ind]
            
            contour_f = ax_.contourf(xx, ts.lev.values,
                        ts,
                        levels = contourf_levels,
                        cmap = cmap)

            

            if contour_dict is not None:
                contours = ax_.contour(xx, contour_ts.lev.values, contour_ts, colors='black', levels = contour_var_levels, alpha = 0.5) 
                ax_.clabel(contours, inline=True, fontsize=font_size, colors='black')
            else:
                contours = ax_.contour(xx, ts.lev.values, ts, colors='white', levels=contour_f.levels)
                ax_.clabel(contours, inline=True, fontsize=font_size, colors='white')
            

            if quiver_dicts is not None:
                if quiver_axis == 'U':
                    U = qv_ts.values
                    V = np.zeros_like(qv_ts)
                else:
                    V = qv_ts.values
                    U = np.zeros_like(qv_ts)
                ax_.quiver(xx[::2], qv_ts.lev.values[::2], np.nan_to_num(U[::2,::2], nan=0.0), np.nan_to_num(V[::2,::2], nan=0.0),  alpha = 0.5, width = width, headwidth = headwidth,headlength = headlength, headaxislength = headaxislength)

            title = title + ' - ' +  f'lon: {ts.lon.values}$^o$ east' + ' - ' + ds +' - ' + bms_label 
            # if  type(yeartoplot) == str:
            #     title = title + f'\n {yeartoplot}'   + f' {season}'
            
            # elif all([len(yeartoplot) >1 , type(yeartoplot) == list]):
            title = title + f'\n composite'  + f' {season}'
            # else:
            #     if yeartoplot[0] not in states_dict.keys():
            #         states_dict[yeartoplot[0]] = 'NA'
            #     title = title +  f'\n {yeartoplot[0]}' + f' {season}' + f': ENSO {states_dict[yeartoplot[0]]}'
            
            if contour_dict is not None:
                title = title + f' - {contour_var} contours'
            
            ax_.set_title(title , fontsize = font_size)           
            ax_.set_ylabel('depth (m)', fontsize = font_size)
            ax_.set_xlabel('lat (deg)', fontsize = font_size)
            ax_.invert_yaxis()
            ax_.tick_params(axis='x', which='major', labelsize=font_size)
            ax_.tick_params(axis='y', which='major', labelsize=font_size)

        cbar = fig.colorbar(contour_f , ax=ax, orientation='vertical', label = colorbar_label)
        cbar.ax.set_ylabel(colorbar_label ,fontsize=font_size) 
        cbar.ax.tick_params(labelsize=font_size)
    if save:
        if show is False:
            print(f"If show = {show} then save = False --change to show = True")
        else:
            Path(dir_name).mkdir(parents=True,
                                     exist_ok=True)
            plt.savefig(f'{dir_name}/{file_name}.png',
                            bbox_inches='tight',
                            dpi=300)   
            

def snapshot_depth_vs_lon_quiver(ds_list,
                            quiver_dicts_X, 
                            quiver_dicts_Y,
                            bms_label,
                            yeartoplot = None, 
            ldyr=1-1,
            title='',
            figsize=(10,45),
            contour_dict = None,
            contourf_levels = [0,0.5, 1,2,4,6, 8,10,15,20,25,30,35,40,50],
            cmap = 'viridis',
            dir_name=None,
            file_name=None,
            colorbar_label = None, 
            season = 'ANN',
            lev_interp = None,
            lev_range = None,
            show=False,
            return_fig_handles = False,
            save=False,
            arrow_scale = None,
            gapped_quiver = False,
            headwidth = 5,
            headlength = 1, 
            headaxislength = 2,
            width =  0.005,
            font_size = 14,
            scale_factor_X = 1,
            scale_factor_Y = 1):
    '''
     to plot data written in terms of target years
     on a common period for each lead year
     
    '''
        
    for jj in range(4):
        sea = [ii*12+3*jj+np.arange(3) for ii in range(ldyr,
                                                        ldyr + 1 )]
        seas = list(np.stack(sea,
                            axis=0).flatten())
        if jj == 0:
            JFM = seas
            # print(f"JFM:{seas}")
        if jj == 1:
            AMJ = seas
            # print(f"AMJ:{seas}")
        if jj == 2:
            JAS = seas
            # print(f"JAS:{seas}")
        if jj == 3:
            OND = seas
            # print(f"OND:{seas}")
            # print('======')
    # print(f"ANN:{np.arange(12) + 12* ldyr}")
    # print('======')
    seasons = {'JFM': JFM,
            'AMJ': AMJ,
            'JAS': JAS,
            'OND': OND,
            'ANN': np.arange(12) + 12* ldyr}

    for ind, month in enumerate(['Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']):
            seasons[month] = [ii*12+ ind for ii in range(ldyr, ldyr + 1 )]
    
    seasons['DJF'] = [0,1,2]
    seasons['MAM'] = [2,3,4]
    seasons['JJA'] = [5,6,7]
    seasons['SON'] = [8,9,10]
    seasons['ENSO_diff'] = [-1]

    print( f'Scale Factor Y : {scale_factor_Y} \n Scale Factor X : {scale_factor_X} \n ====================')
    assert all([quiver_dicts_Y is None, quiver_dicts_X is  None]) is False

    if season in ['DJF','MAM','JJA','SON']:
        assert ldyr == 0
        for ds in ds_list:
            assert 'hindcast' not in ds
    if season not in seasons.keys():
        season = 'ENSO_diff'
        raise Warning('Season was set to ENSO_diff !')
    if season == 'ENSO_diff':
        yeartoplot = 'diff'

    if show:
        ldyr1 = ldyr + 1
    

        fig, ax = plt.subplots(len(ds_list), 1 ,figsize=figsize)
        if quiver_dicts_X is not None:
            ds_dict_X = quiver_dicts_X[bms_label]
        else:
            ds_dict_X = None
        if quiver_dicts_Y is not None:
            ds_dict_Y = quiver_dicts_Y[bms_label]
        else:
            ds_dict_Y = None
        if contour_dict is not None:
            ds_dict = contour_dict[bms_label]

        for ind, ds in enumerate(ds_list):
            # if len(ds.split('-')) == 1:
            if ds_dict_X is None:
                    dstpY = ds_dict_Y[ds] * scale_factor_Y
                    dstpX = xr.zeros_like(dstpY)
            elif ds_dict_Y is None:
                dstpX = ds_dict_X[ds] * scale_factor_X
                dstpY = xr.zeros_like(dstpX)                
            else:
                dstpY = ds_dict_Y[ds] * scale_factor_Y
                dstpX = ds_dict_X[ds] * scale_factor_X
            if contour_dict is not None:
                dstp = ds_dict[ds]


            if season == 'DJF':
                try:
                    tsX = DJFy(dstpX).sel(time=slice(ldyr*12,(ldyr1)*12-1))
                    tsY = DJFy(dstpY).sel(time=slice(ldyr*12,(ldyr1)*12-1))
                    if contour_dict is not None:
                        ts = DJFy(dstp).sel(time=slice(ldyr*12,(ldyr1)*12-1))
                except:
                    tsX = dstpX.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                    tsY = dstpY.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                    if contour_dict is not None:
                        ts = dstp.sel(time=slice(ldyr*12,(ldyr1)*12-1))                                                        
                    print('DJFy unseccessfull, make sure the conversion is already done!')
            else:
                
                    tsX = dstpX.sel(time=slice(ldyr*12,(ldyr1)*12-1)) if season not in ['ENSO_diff'] else dstpX
                    tsY = dstpY.sel(time=slice(ldyr*12,(ldyr1)*12-1)) if season not in ['ENSO_diff'] else dstpY
                    if contour_dict is not None:
                        ts = dstp.sel(time=slice(ldyr*12,(ldyr1)*12-1))   if season not in ['ENSO_diff'] else dstp                                                   
            if yeartoplot is not None:
                if all([type(yeartoplot) == str, 'diff' not in yeartoplot]): 
                
                    y1 = eval(yeartoplot.split('-')[0])
                    y0 = eval(yeartoplot.split('-')[1])
                    tsX = (tsX.sel(time = seasons[season] ).sel(year = y1 ).mean('time') - tsX.sel(time = seasons[season] ).sel(year = y0).mean('time')).squeeze()
                    tsY = (tsY.sel(time = seasons[season] ).sel(year = y1 ).mean('time') - tsY.sel(time = seasons[season] ).sel(year = y0).mean('time')).squeeze()
                    if contour_dict is not None:
                        ts = (ts.sel(time = seasons[season] ).sel(year = y1 ).mean('time') - ts.sel(time = seasons[season] ).sel(year = y0).mean('time')).squeeze()
                else:
                    tsX = tsX.sel(time = seasons[season] ).sel(year = yeartoplot ).mean(['year','time']).squeeze()
                    tsY = tsY.sel(time = seasons[season] ).sel(year = yeartoplot ).mean(['year','time']).squeeze()
                    if contour_dict is not None:
                        ts = ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean(['year','time']).squeeze()
            else:
                tsX = tsX.sel(time = seasons[season] ).mean(['year','time']).squeeze()
                tsY = tsY.sel(time = seasons[season] ).mean(['year','time']).squeeze()
                if contour_dict is not None:
                    ts = ts.sel(time = seasons[season] ).mean(['year','time']).squeeze()
            xx = tsX.lon.values
            
            if lev_interp is not None:
                tsX = tsX.interp(lev = lev_interp )
                tsY = tsY.interp(lev = lev_interp )
                if contour_dict is not None:
                    ts = ts.interp(lev = lev_interp )

            if lev_range is not None:
                if type(lev_range) is not list:
                    ls = [0, lev_range]
                    lev_range = ls
                    del ls
                tsX = tsX.where((tsX.lev <= lev_range[1]) & (tsX.lev >= lev_range[0]) , drop = True)
                tsY = tsY.where((tsY.lev <= lev_range[1]) & (tsY.lev >= lev_range[0]), drop = True)
                if contour_dict is not None:
                    ts = ts.where((ts.lev <= lev_range[1]) & (ts.lev >= lev_range[0]), drop = True)

            label = f'{ds}'
            if len(ds_list) == 1:
                ax_ = ax
            else:
                ax_ = ax[ind]

            if contour_dict is not None:
                contour_f = ax_.contourf(xx, ts.lev.values,
                            ts,
                            levels = contourf_levels,
                            cmap = cmap
                            # vmax = 30,
                            # vmin = 0,
                            # ds_dict[ds]['linestyle'],
                                )
                contours = ax_.contour(xx, ts.lev.values, ts, colors='white', levels=contour_f.levels)
                ax_.clabel(contours, inline=True, fontsize=font_size, colors='white')

            if gapped_quiver is None:
                step = 1
            else:
                step = gapped_quiver   
            ax_.quiver(xx[::step], tsX.lev.values[::step], np.nan_to_num(tsX[::step,::step], nan=0.0), np.nan_to_num(tsY[::step,::step], nan=0.0),  alpha = 0.5, scale = arrow_scale, headwidth = headwidth,headlength = headlength, headaxislength = headaxislength, width = width)
            
            # if  type(yeartoplot) == str:
            #     ax_.set_title(title + ' - ' +  ds +' - ' + bms_label + f' {yeartoplot}' + f' {season}', fontsize = font_size)
            
            # elif all([len(yeartoplot) >1 , type(yeartoplot) == list]):
            ax_.set_title(title + ' - ' +  ds +' - ' + bms_label + f' composite' + f' {season}', fontsize = font_size)

            # else:
            #     if yeartoplot[0] not in states_dict.keys():
            #         states_dict[yeartoplot[0]] = 'NA'
            #     ax_.set_title(title + ' - ' +  ds +' - ' + bms_label + f' {yeartoplot[0]}' + f' {season}' + f': ENSO {states_dict[yeartoplot[0]]}', fontsize = font_size)
            ax_.set_ylabel('depth (m)', fontsize = font_size)
            ax_.set_xlabel('Lon ($^o$ East)', fontsize = font_size)
            if ind < len(ds_list) - 1:
                ax_.set_xlabel('', fontsize = font_size)
            ax_.set_ylim(lev_range)
            ax_.invert_yaxis()
            ax_.tick_params(axis='x', which='major', labelsize=font_size)
            ax_.tick_params(axis='y', which='major', labelsize=font_size)

        if contour_dict is not None:
            cbar = fig.colorbar(contour_f , ax=ax, orientation='vertical', label = colorbar_label)
            cbar.ax.set_ylabel(colorbar_label ,fontsize=font_size)
            cbar.ax.tick_params(labelsize=font_size)
    if save:
        if show is False:
            print(f"If show = {show} then save = False --change to show = True")
        else:
            Path(dir_name).mkdir(parents=True,
                                     exist_ok=True)
            plt.savefig(f'{dir_name}/{file_name}.png',
                            bbox_inches='tight',
                            dpi=300)   

def snapshot_crossection_quiver(ds_list,
                            quiver_dicts_X, 
                            quiver_dicts_Y,
                            bms_label,
                            longitude,
                            yeartoplot = None,
            ldyr=1-1,
            title='',
            figsize=(10,45),
            contour_dict = None,
            contourf_levels = [0,0.5, 1,2,4,6, 8,10,15,20,25,30,35,40,50],
            cmap = 'viridis',
            dir_name=None,
            file_name=None,
            colorbar_label = None, 
            season = 'ANN',
            lev_interp = None,
            lev_range:list = None,
            show=False,
            return_fig_handles = False,
            save=False            ,
            arrow_scale = None, 
            gapped_quiver = False,
            headwidth = 5,
            headlength = 1, 
            headaxislength = 2,
            width = 0.005,
            font_size = 15,
            scale_factor_X = 1,
            scale_factor_Y = 1):
    '''
     to plot data written in terms of target years
     on a common period for each lead year
     
    '''
        
    for jj in range(4):
        sea = [ii*12+3*jj+np.arange(3) for ii in range(ldyr,
                                                        ldyr + 1 )]
        seas = list(np.stack(sea,
                            axis=0).flatten())
        if jj == 0:
            JFM = seas
            # print(f"JFM:{seas}")
        if jj == 1:
            AMJ = seas
            # print(f"AMJ:{seas}")
        if jj == 2:
            JAS = seas
            # print(f"JAS:{seas}")
        if jj == 3:
            OND = seas
            # print(f"OND:{seas}")
            # print('======')
    # print(f"ANN:{np.arange(12) + 12* ldyr}")
    # print('======')
    seasons = {'JFM': JFM,
            'AMJ': AMJ,
            'JAS': JAS,
            'OND': OND,
            'ANN': np.arange(12) + 12* ldyr}

    for ind, month in enumerate(['Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']):
            seasons[month] = [ii*12+ ind for ii in range(ldyr, ldyr + 1 )]
    
    seasons['DJF'] = [0,1,2]
    seasons['MAM'] = [2,3,4]
    seasons['JJA'] = [5,6,7]
    seasons['SON'] = [8,9,10]
    # seasons['ENSO_diff'] = [-1]

    print( f'Scale Factor Y : {scale_factor_Y} \n Scale Factor X : {scale_factor_X} \n ====================')

    assert all([quiver_dicts_Y is None, quiver_dicts_X is  None]) is False

    if season in ['DJF','MAM','JJA','SON']:
        assert ldyr == 0
        for ds in ds_list:
            assert 'hindcast' not in ds
    if season not in seasons.keys():
        season = 'ENSO_diff'
        raise Warning('Season was set to ENSO_diff !')
    if season == 'ENSO_diff':
        yeartoplot = ['diff']

    if show:
        ldyr1 = ldyr + 1
    

        fig, ax = plt.subplots(len(ds_list), 1 ,figsize=figsize)
        if quiver_dicts_X is not None:
            ds_dict_X = quiver_dicts_X[bms_label]
        else:
            ds_dict_X = None
        if quiver_dicts_Y is not None:
            ds_dict_Y = quiver_dicts_Y[bms_label]
        else:
            ds_dict_Y = None
        if contour_dict is not None:
            ds_dict = contour_dict[bms_label]
        
        for ind, ds in enumerate(ds_list):
            # if len(ds.split('-')) == 1:
            if ds_dict_X is None:
                dstpY = ds_dict_Y[ds] * scale_factor_Y
                dstpX = xr.zeros_like(dstpY)
            elif ds_dict_Y is None:
                dstpX = ds_dict_X[ds] * scale_factor_X
                dstpY = xr.zeros_like(dstpX)                
            else:
                dstpY = ds_dict_Y[ds] * scale_factor_Y
                dstpX = ds_dict_X[ds] * scale_factor_X
            if contour_dict is not None:
                dstp = ds_dict[ds]


            if season == 'DJF':
                try:
                    tsX = DJFy(dstpX).sel(time=slice(ldyr*12,(ldyr1)*12-1))
                    tsY = DJFy(dstpY).sel(time=slice(ldyr*12,(ldyr1)*12-1))
                    if contour_dict is not None:
                        ts = DJFy(dstp).sel(time=slice(ldyr*12,(ldyr1)*12-1))
                                                        
                except:
                    tsX = dstpX.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                    tsY = dstpY.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                    if contour_dict is not None:
                        ts = dstp.sel(time=slice(ldyr*12,(ldyr1)*12-1))

                    print('DJFy unseccessfull, make sure the conversion is already done!')
            else:
                
                    tsX = dstpX.sel(time=slice(ldyr*12,(ldyr1)*12-1))  if season not in ['ENSO_diff'] else dstpX
                    tsY = dstpY.sel(time=slice(ldyr*12,(ldyr1)*12-1))  if season not in ['ENSO_diff'] else dstpY
                    if contour_dict is not None:
                        ts = dstp.sel(time=slice(ldyr*12,(ldyr1)*12-1))  if season not in ['ENSO_diff'] else dstp
                                   
            if yeartoplot is not None:
                if all([type(yeartoplot) == str, 'diff' not in yeartoplot]): 
                    y1 = eval(yeartoplot.split('-')[0])
                    y0 = eval(yeartoplot.split('-')[1])
                    tsX = (tsX.sel(time = seasons[season] ).sel(year = y1 ).mean('time') - tsX.sel(time = seasons[season] ).sel(year = y0).mean('time')).squeeze()
                    tsY = (tsY.sel(time = seasons[season] ).sel(year = y1 ).mean('time') - tsY.sel(time = seasons[season] ).sel(year = y0).mean('time')).squeeze()
                    if contour_dict is not None:
                        ts = (ts.sel(time = seasons[season] ).sel(year = y1 ).mean('time') - ts.sel(time = seasons[season] ).sel(year = y0).mean('time')).squeeze()

                else:
                    tsX = tsX.sel(time = seasons[season] ).sel(year = yeartoplot ).mean(['year','time']).squeeze()
                    tsY = tsY.sel(time = seasons[season] ).sel(year = yeartoplot ).mean(['year','time']).squeeze()
                    if contour_dict is not None:
                        ts = ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean(['year','time']).squeeze()
            else:
                tsX = tsX.sel(time = seasons[season] ).mean(['year','time']).squeeze()
                tsY = tsY.sel(time = seasons[season] ).mean(['year','time']).squeeze()
                if contour_dict is not None:
                    ts = ts.sel(time = seasons[season] ).mean(['year','time']).squeeze()
            xx = tsX.lat.values
            
            if lev_interp is not None:
                tsX = tsX.interp(lev = lev_interp )
                tsY = tsY.interp(lev = lev_interp )
                if contour_dict is not None:
                    ts = ts.interp(lev = lev_interp )
            
            if lev_range is not None:
                if type(lev_range) is not list:
                    ls = [0, lev_range]
                    lev_range = ls
                    del ls
                tsX = tsX.where((tsX.lev <= lev_range[1]) & (tsX.lev >= lev_range[0]) , drop = True)
                tsY = tsY.where((tsY.lev <= lev_range[1]) & (tsY.lev >= lev_range[0]), drop = True)
                if contour_dict is not None:
                    ts = ts.where((ts.lev <= lev_range[1]) & (ts.lev >= lev_range[0]), drop = True)
                
            tsX = tsX.sel(lon = longitude, method = 'nearest')
            tsY = tsY.sel(lon = longitude, method = 'nearest')
            if contour_dict is not None:
                ts = ts.sel(lon = longitude, method = 'nearest')

            label = f'{ds}'
            if len(ds_list) == 1:
                ax_ = ax
            else:
                ax_ = ax[ind]
        
            
            if contour_dict is not None:
                contour_f = ax_.contourf(xx, ts.lev.values,
                            ts,
                            levels = contourf_levels,
                            cmap = cmap)
                contours = ax_.contour(xx, ts.lev.values, ts, colors='white', levels=contour_f.levels)
                ax_.clabel(contours, inline=True, fontsize=font_size, colors='white')

            if gapped_quiver is None:
                step = 1
            else:
                step = gapped_quiver   
            ax_.quiver(xx[::step], tsX.lev.values[::step], np.nan_to_num(tsX.values[::step,::step], nan=0.0), np.nan_to_num(tsY.values[::step,::step], nan=0.0), scale = arrow_scale, headwidth = headwidth,headlength = headlength, headaxislength = headaxislength, width = width)


            # if  type(yeartoplot) == str:
            #     ax_.set_title(title + ' - ' +  f'lon: {tsX.lon.values} deg east' + ' - ' + ds +' - ' + bms_label + f'\n {yeartoplot}' + f' {season}', fontsize = font_size)
            
            # elif all([len(yeartoplot) >1 , type(yeartoplot) == list]):
            ax_.set_title(title + ' - ' +  f'lon: {tsX.lon.values}$^o$ east' + ' - ' + ds +' - ' + bms_label + f'\n composite' + f' {season}', fontsize = font_size)

            # else:
            #     if yeartoplot[0] not in states_dict.keys():
            #         states_dict[yeartoplot[0]] = 'NA'
            #     ax_.set_title(title + ' - ' + f'lon: {tsX.lon.values} deg east' + ' - ' +  ds +' - ' + bms_label + f'\n {yeartoplot[0]}' + f' {season}' + f': ENSO {states_dict[yeartoplot[0]]}', fontsize = font_size)
            ax_.set_ylabel('depth (m)', fontsize = font_size)
            ax_.set_xlabel('lat (deg)', fontsize = font_size)
            ax_.tick_params(axis='x', which='major', labelsize=font_size)
            ax_.tick_params(axis='y', which='major', labelsize=font_size)
            ax_.set_ylim(lev_range)
            ax_.invert_yaxis()
        if contour_dict is not None:
            cbar = fig.colorbar(contour_f , ax=ax, orientation='vertical', label = colorbar_label)
            cbar.ax.set_ylabel(colorbar_label ,fontsize=font_size)
            cbar.ax.tick_params(labelsize=font_size)
    
    if save:
        if show is False:
            print(f"If show = {show} then save = False --change to show = True")
        else:
            Path(dir_name).mkdir(parents=True,
                                     exist_ok=True)
            plt.savefig(f'{dir_name}/{file_name}.png',
                            bbox_inches='tight',
                            dpi=300)   
            

def snapshot_arial_quiver(ds_list,
                            quiver_dicts_X, 
                            quiver_dicts_Y,
                            bms_label,
                            yeartoplot = None,
                            quiver_dicts = None, 
            ldyr=1-1,
            title='',
            figsize=(10,45),
            vmax=None,
            vmin=None,
            cmap = 'viridis',
            dir_name=None,
            file_name=None,
            colorbar_dict = None,
            colorbar_label = None, 
            season = 'ANN',
            lev_interp = None,
            lev_range:list = None,
            show=False,
            return_fig_handles = False,
            save=False,
            arrow_scale = None, 
            gapped_quiver = None,
            headwidth = 5,
            headlength = 1, 
            headaxislength = 2,
            width =  0.005,
            font_size = 15,
            scale_factor_X = 1,
            scale_factor_Y = 1):
    '''
     to plot data written in terms of target years
     on a common period for each lead year
     
    '''
        
    for jj in range(4):
        sea = [ii*12+3*jj+np.arange(3) for ii in range(ldyr,
                                                        ldyr + 1 )]
        seas = list(np.stack(sea,
                            axis=0).flatten())
        if jj == 0:
            JFM = seas
            # print(f"JFM:{seas}")
        if jj == 1:
            AMJ = seas
            # print(f"AMJ:{seas}")
        if jj == 2:
            JAS = seas
            # print(f"JAS:{seas}")
        if jj == 3:
            OND = seas
            # print(f"OND:{seas}")
            # print('======')
    # print(f"ANN:{np.arange(12) + 12* ldyr}")
    # print('======')
    seasons = {'JFM': JFM,
            'AMJ': AMJ,
            'JAS': JAS,
            'OND': OND,
            'ANN': np.arange(12) + 12* ldyr}

    for ind, month in enumerate(['Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']):
            seasons[month] = [ii*12+ ind for ii in range(ldyr, ldyr + 1 )]
    
    seasons['DJF'] = [0,1,2]
    seasons['MAM'] = [2,3,4]
    seasons['JJA'] = [5,6,7]
    seasons['SON'] = [8,9,10]
    seasons['ENSO_diff'] = [-1]

    print( f'Scale Factor Y : {scale_factor_Y} \n Scale Factor X : {scale_factor_X} \n ====================')
    assert all([quiver_dicts_Y is None, quiver_dicts_X is  None]) is False

    if season in ['DJF','MAM','JJA','SON']:
        assert ldyr == 0
        for ds in ds_list:
            assert 'hindcast' not in ds
    if season not in seasons.keys():
        season = 'ENSO_diff'
        raise Warning('Season was set to ENSO_diff !')
    if season == 'ENSO_diff':
        yeartoplot = ['diff']

    if show:
        ldyr1 = ldyr + 1
    

        fig, ax = plt.subplots(len(ds_list), 1 ,figsize=figsize)
        if quiver_dicts_X is not None:
            ds_dict_X = quiver_dicts_X[bms_label]
        else:
            ds_dict_X = None
        if quiver_dicts_Y is not None:
            ds_dict_Y = quiver_dicts_Y[bms_label]
        else:
            ds_dict_Y = None
        if colorbar_dict is not None:
            ds_dict = colorbar_dict[bms_label]

        for ind, ds in enumerate(ds_list):
            # if len(ds.split('-')) == 1:
            if ds_dict_X is None:
                dstpY = ds_dict_Y[ds] * scale_factor_Y
                dstpX = xr.zeros_like(dstpY)
            elif ds_dict_Y is None:
                dstpX = ds_dict_X[ds] * scale_factor_X
                dstpY = xr.zeros_like(dstpX)                
            else:
                dstpY = ds_dict_Y[ds] * scale_factor_Y
                dstpX = ds_dict_X[ds] * scale_factor_X
            if colorbar_dict is not None:
                dstp = ds_dict[ds]
            
            dstpY_lev = True if 'lev' in dstpY.dims else False
            dstpX_lev = True if 'lev' in dstpX.dims else False  
            if colorbar_dict is not None:
                dstp_lev = True if 'lev' in dstp.dims else False 

            if season == 'DJF':
                try:
                    tsX = DJFy(dstpX).sel(time=slice(ldyr*12,(ldyr1)*12-1))
                    tsY = DJFy(dstpY).sel(time=slice(ldyr*12,(ldyr1)*12-1))
                    if colorbar_dict is not None:
                        ts = DJFy(dstp).sel(time=slice(ldyr*12,(ldyr1)*12-1))
                                                        
                except:
                    tsX = dstpX.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                    tsY = dstpY.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                    if colorbar_dict is not None:
                        ts = dstp.sel(time=slice(ldyr*12,(ldyr1)*12-1))

                    print('DJFy unseccessfull, make sure the conversion is already done!')
            else:
                
                    tsX = dstpX.sel(time=slice(ldyr*12,(ldyr1)*12-1))  if season not in ['ENSO_diff'] else dstpX
                    tsY = dstpY.sel(time=slice(ldyr*12,(ldyr1)*12-1))  if season not in ['ENSO_diff'] else dstpY
                    if colorbar_dict is not None:
                        ts = dstp.sel(time=slice(ldyr*12,(ldyr1)*12-1))  if season not in ['ENSO_diff'] else dstp
                                                       
            if yeartoplot is not None:
                if all([type(yeartoplot) == str, 'diff' not in yeartoplot]): 
                    y1 = eval(yeartoplot.split('-')[0])
                    y0 = eval(yeartoplot.split('-')[1])
                    tsX = (tsX.sel(time = seasons[season] ).sel(year = y1 ).mean('time') - tsX.sel(time = seasons[season] ).sel(year = y0).mean('time')).squeeze()
                    tsY = (tsY.sel(time = seasons[season] ).sel(year = y1 ).mean('time') - tsY.sel(time = seasons[season] ).sel(year = y0).mean('time')).squeeze()
                    if colorbar_dict is not None:
                        ts = (ts.sel(time = seasons[season] ).sel(year = y1 ).mean('time') - ts.sel(time = seasons[season] ).sel(year = y0).mean('time')).squeeze()

                else:
                    tsX = tsX.sel(time = seasons[season] ).sel(year = yeartoplot ).mean(['year','time']).squeeze()
                    tsY = tsY.sel(time = seasons[season] ).sel(year = yeartoplot ).mean(['year','time']).squeeze()
                    if colorbar_dict is not None:
                        ts = ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean(['year','time']).squeeze()
            else:
                tsX = tsX.sel(time = seasons[season] ).mean(['year','time']).squeeze()
                tsY = tsY.sel(time = seasons[season] ).mean(['year','time']).squeeze()
                if colorbar_dict is not None:
                    ts = ts.sel(time = seasons[season] ).mean(['year','time']).squeeze()
            xx = tsX.lon.values
            
            if lev_interp is not None:
                tsX = tsX.interp(lev = lev_interp )
                tsY = tsY.interp(lev = lev_interp )
                if colorbar_dict is not None:
                    ts = ts.interp(lev = lev_interp )
            
            if lev_range is not None:
                tsX = tsX.where(tsX.lev >= lev_range[0], drop = True) if dstpX_lev else tsX
                tsX = tsX.where(tsX.lev <= lev_range[1], drop = True) if dstpX_lev else tsX
                tsY = tsY.where(tsY.lev >= lev_range[0], drop = True) if dstpY_lev else tsY
                tsY = tsY.where(tsY.lev <= lev_range[1], drop = True) if dstpY_lev else tsY
                if colorbar_dict is not None:
                    ts = ts.where(ts.lev >= lev_range[0], drop = True) if dstp_lev else ts
                    ts = ts.where(ts.lev <= lev_range[1], drop = True) if dstp_lev else ts
                
            tsX = tsX.mean('lev') if dstpX_lev else tsX
            tsY = tsY.mean('lev') if dstpY_lev else tsY
            if colorbar_dict is not None:
                ts = ts.mean('lev') if dstp_lev else ts


            label = f'{ds}'
            if len(ds_list) == 1:
                ax_ = ax
            else:
                ax_ = ax[ind]
            
            if colorbar_dict is not None:
                im = ax_.pcolormesh(ts.lon, ts.lat, ts, 
                                cmap=cmap, vmin = vmin, vmax = vmax)
            if gapped_quiver is None:
                step = 1
            else:
                step = gapped_quiver                
            ax_.quiver(xx[::step], tsX.lat.values[::step], np.nan_to_num(tsX.values[::step,::step], nan=0.0), np.nan_to_num(tsY.values[::step,::step], nan=0.0), scale = arrow_scale, headwidth = headwidth,headlength = headlength, headaxislength = headaxislength, width = width)

            # if  type(yeartoplot) == str:
            #     ax_.set_title(title + ' - ' + f'{lev_range[0]}-{lev_range[1]} (m)' + ' - ' +  ds +' - ' + bms_label + f' {yeartoplot}' + f' {season}', fontsize = font_size)
            
            # elif all([len(yeartoplot) >1 , type(yeartoplot) == list]):
            ax_.set_title(title + ' - ' + f'{lev_range[0]}-{lev_range[1]} (m)' + ' - ' +  ds +' - ' + bms_label + f' composite' + f' {season}', fontsize = font_size)
            # else:
            #     if yeartoplot[0] not in states_dict.keys():
            #         states_dict[yeartoplot[0]] = 'NA'
            #     ax_.set_title(title + ' - ' + f'{lev_range[0]}-{lev_range[1]} (m)' + ' - ' +  ds +' - ' + bms_label + f' {yeartoplot[0]}' + f' {season}' + f': ENSO {states_dict[yeartoplot[0]]}', fontsize = font_size)
            # ax_.set_ylabel('lat (deg North)', fontsize = font_size)
            if ind < len(ds_list) -1 :
                ax_.set_xlabel('')
            else:
                ax_.set_xlabel('Lon ($^o$ East)', fontsize = font_size)
            ax_.tick_params(axis='x', which='major', labelsize=font_size)
            ax_.tick_params(axis='y', which='major', labelsize=font_size)
            # ax_.invert_yaxis()

        if colorbar_dict is not None:
            cbar = fig.colorbar(im , ax=ax, orientation='vertical', label = colorbar_label)
            cbar.ax.set_ylabel(colorbar_label ,fontsize=font_size)  
            cbar.ax.tick_params(labelsize=font_size)

    if save:
        if show is False:
            print(f"If show = {show} then save = False --change to show = True")
        else:
            Path(dir_name).mkdir(parents=True,
                                     exist_ok=True)
            plt.savefig(f'{dir_name}/{file_name}.png',
                            bbox_inches='tight',
                            dpi=300)   
            



def trend(ds, dim = 'year', return_detrended = False ):
    m = ds.polyfit( dim  = dim, deg = 1).polyfit_coefficients
    out = m[0]*ds[dim] +  m[1]
    if return_detrended:
        return  out, ds - out
    else:
        return out












