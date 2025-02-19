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
from module_global_averages import area_weighted_avg
from module_metrics import rmse, corr, corr_map
from module_data_postprocessing import spatial_mask





def add_cyclic_point(ds):
    add = ds.isel(lon = -1).assign_coords(lon = ds.isel(lon = -1).lon.values + 1)
    return xr.concat([ds,add], dim = 'lon')


def plot_composites(ds_list,
                    data_dict,
                    mask = None,
                    specific_years:list = None,
                    figsize=(12,12),                    
                    central_longitude=260,
                    ldyr_ini=0,            
                    ldyr_end=1,
                    vmax=2,
                    vmin=-2,
                    cmap = 'RdBu_r',
                    cbar_label=r'mol m$^{-2}$ yr$^{-1}$',  
                    std  = False,        
                    individual_months = False,  
                    shifted_seasons = False,   
                    pattern_corr = 1,      
                    dir_name=None,
                    file_name=None,
                    save=False):
    
    print('======')
    print(f"init year starts: {data_dict['obs'].year.values[0] + ldyr_ini}")
    print(f"init year ends  : {data_dict['obs'].year.values[-1]}")
    print(f"ave lead years  : {ldyr_ini+1} to {ldyr_end-1+1}")
    print('======')

    plt.figure(figsize=figsize)    
    
    nseas = 4
    nmnth = 12
    nldyr = ldyr_end - ldyr_ini 
    if individual_months:
        seasons = {}
        for ind, month in enumerate(['Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']):
            seasons[month] = [ii*12+ ind for ii in range(ldyr_ini, ldyr_end)]
            
    elif shifted_seasons:
        assert ldyr_ini ==0, ldyr_end == 1
        for ds in ds_list:
            assert 'hindcast' not in ds 
        seasons = {}
        seasons['DJF'] = [0,1,2] 
        seasons['MAM'] = [2,3,4]
        seasons['JJA'] = [5,6,7]
        seasons['SON'] = [8,9,10]
        seasons['ANN'] = np.arange(12)

    else:
        for jj in range(nseas):
            sea = [ii*nmnth+3*jj+np.arange(3) for ii in range(ldyr_ini,
                                                            ldyr_end)]
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

        print(f"ANN:{np.arange(nldyr*nmnth) + ldyr_ini * 12}")
        print('======')

        seasons = {'JFM': JFM,
                'AMJ': AMJ,
                'JAS': JAS,
                'OND': OND,
                'ANN': np.arange(nldyr*nmnth) + 12 * ldyr_ini}
    
    i = 0

    if specific_years is not None:
        if len(specific_years) >= 1:
                assert std is False, 'Need at least two years to caluculate interannual std '

    for season, inds in seasons.items(): 

        if season == 'DJF':
            obs_ref = DJFy(data_dict['obs'])
        else:
            obs_ref = data_dict['obs']

        if specific_years is not None:
            obs_ref = obs_ref.sel(year = specific_years)

        if std:
            obs_ref = obs_ref.sel(time=inds)
            obs_ref = obs_ref.mean(['time']).std('year')                                                   
        else:
            obs_ref = obs_ref.sel(time=inds).mean(['year',
                                                            'time'])   
        
        for ind, ds in enumerate(ds_list):
            
            ax = plt.subplot(len(seasons),
                             len(ds_list),
                             i*len(ds_list)+ind+1,
                             projection=ccrs.Robinson(central_longitude=central_longitude))
            if season == 'DJF':
                ds_toplot = DJFy(data_dict[ds])
            else:
                ds_toplot = data_dict[ds] 

            if specific_years is not None:
                ds_toplot = ds_toplot.sel(year = specific_years)
            
            if std:
                ds_toplot = ds_toplot.sel(time=inds)
                ds_toplot = ds_toplot.mean(['time']).std('year')  
            else:
                ds_toplot = ds_toplot.sel(time=inds).mean(['year',
                                                                'time'])
            if central_longitude != 0:
                cb = ax.pcolormesh(add_cyclic_point(ds_toplot).lon,
                                add_cyclic_point(ds_toplot).lat,
                                add_cyclic_point(ds_toplot),
                                cmap=plt.cm.get_cmap(cmap),
                                vmax=vmax,
                                vmin=vmin,
                                rasterized=True,
                                transform=ccrs.PlateCarree()) 
            else:
                cb = ax.pcolormesh(ds_toplot.lon,
                                ds_toplot.lat,
                                ds_toplot,
                                cmap=plt.cm.get_cmap(cmap),
                                vmax=vmax,
                                vmin=vmin,
                                rasterized=True,
                                transform=ccrs.PlateCarree())
            
            _ = ax.coastlines()

            ax.set_ylabel('')
            ax.set_xlabel('')
            if mask is None:
                mask = spatial_mask(ds_toplot)
            glbavg = np.round(area_weighted_avg(ds_toplot,
                                                mask=mask.values).values,4)
            
            if pattern_corr == 1:
                if season == 'DJF':
                    corr_pat  = corr_map(DJFy(data_dict['obs']).sel(time=inds).mean(['time']),
                                        DJFy(data_dict[ds]).sel(time=inds).mean(['time']))
                else:
                    corr_pat  = corr_map(data_dict['obs'].sel(time=inds).mean(['time']),
                                        data_dict[ds].sel(time=inds).mean(['time']))
                corr_pat_avg  = np.round(corr_pat.values.mean(),2)
                ax.set_title(f'{str(glbavg)}, {str(corr_pat_avg)}',
                         fontsize=20)
            elif pattern_corr == 2:
                if season == 'DJF':
                    corr_pat  = corr_map(DJFy(data_dict['obs']).sel(time=inds).mean(['time','year']),
                                        DJFy(data_dict[ds]).sel(time=inds).mean(['time','year']))
                else:
                    corr_pat  = corr_map(data_dict['obs'].sel(time=inds).mean(['time','year']),
                                        data_dict[ds].sel(time=inds).mean(['time','year']))
                corr_pat_avg  = np.round(corr_pat.values,2)
                ax.set_title(f'{str(glbavg)}, {str(corr_pat_avg)}$*$',
                         fontsize=20)

            if i == 0:
                ax.set_title(ds + f'\n {str(glbavg)}, {str(corr_pat_avg)}',
                             fontsize=20)
                
            if ind ==0:
                ax.text(-0.25,
                        1.1,
                        season,
                        fontsize=20,
                        transform=ax.transAxes)
        i = i + 1
    
    divider = make_axes_locatable(ax)
    
    ax_cb   = divider.append_axes('bottom',
                                  size="10%",
                                  pad=0.1,
                                  axes_class=plt.Axes)
    
    cbar    = plt.colorbar(cb,
                           cax=ax_cb,
                           orientation="horizontal")
    if std:
        cbar_label = f'std ({cbar_label})'
    cbar.set_label(label=cbar_label, #r'mol m$^{-2}$ yr$^{-1}$',
                   size=20)
    cbar.ax.tick_params(labelsize=20)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.1,
                        hspace=0.3)

    if save:
        Path(dir_name).mkdir(parents=True,
                             exist_ok=True)
        plt.savefig(f'{dir_name}/{file_name}.png')
        
        
def plot_measures(ds_list,
                  data_dict,
                  measure='rmse',
                  figsize=(12,12),                    
                  central_longitude=260,
                  ldyr_ini=0,            
                  ldyr_end=1,
                  vmax=2,
                  vmin=-2,
                  label='',
                  cmap='RdBu_r',
                  dir_name=None,
                  file_name=None,
                  individual_months = False,
                  shifted_seasons = False,
                  mask = None,
                  save=False):
    
    print('======')
    print(f"init year starts: {data_dict['obs'].year.values[0]  + ldyr_ini}")
    print(f"init year ends  : {data_dict['obs'].year.values[-1]}")
    print(f"ave lead years  : {ldyr_ini+1} to {ldyr_end-1+1}")
    print('======')
    
    plt.figure(figsize=figsize)    
    
    nseas = 4
    nmnth = 12
    nldyr = ldyr_end - ldyr_ini 
    if individual_months:
        seasons = {}
        for ind, month in enumerate(['Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']):
            seasons[month] = [ii*12+ ind for ii in range(ldyr_ini, ldyr_end)]

    elif shifted_seasons:
        assert ldyr_ini == 0, ldyr_end == 1
        for ds in ds_list:
            assert 'hindcast' not in ds 
        seasons = {}
        seasons['DJF'] = [0,1,2] 
        seasons['MAM'] = [2,3,4]
        seasons['JJA'] = [5,6,7]
        seasons['SON'] = [8,9,10]
        seasons['ANN'] = np.arange(12)

    else:
        for jj in range(nseas):
            sea = [ii*nmnth+3*jj+np.arange(3) for ii in range(ldyr_ini,
                                                            ldyr_end)]
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
            
        print(f"ANN:{np.arange(nldyr*nmnth) + ldyr_ini * 12}")
        print('======')
        seasons = {'JFM': JFM,
                'AMJ': AMJ,
                'JAS': JAS,
                'OND': OND,
                'ANN': np.arange(nldyr*nmnth) + ldyr_ini * 12}
    
    i = 0
    
    for season, inds in seasons.items():     
        
        for ind, ds in enumerate(ds_list):
            
            ax = plt.subplot(len(seasons),
                             len(ds_list),
                             i*len(ds_list)+ind+1,
                             projection=ccrs.Robinson(central_longitude=central_longitude))

            if season == 'DJF':
                ds_obs = DJFy(data_dict['obs']).sel(time=inds)
            else:
                ds_obs = data_dict['obs'].sel(time=inds) 

            if season == 'DJF':
                ds_targ = DJFy(data_dict[ds]).sel(time=inds)
            else:
                ds_targ = data_dict[ds].sel(time=inds) 

            ds_obs_clim = ds_obs.mean(['time', 'year'])
            ds_obs = ds_obs.mean('time')
            ds_targ = ds_targ.mean('time')
            metric_dim = 'year'


            if measure == 'rmse':
                ds_toplot = rmse(ds_obs,
                                 ds_targ,metric_dim)
            if measure == 'corr':
                ds_toplot = xr.corr(ds_obs,
                                 ds_targ, dim = metric_dim)
            if measure == 'corr_det':
                ds_toplot = xr.corr(trend(ds_obs, dim = metric_dim, return_detrended = True)[1],
                                 trend(ds_targ, dim = metric_dim, return_detrended= True)[1], dim = metric_dim)
            if measure == 'bias':
                ds_toplot =  ds_targ.mean(metric_dim) - ds_obs.mean(metric_dim)
            if measure == 'MSSS':

                nom = ((ds_obs - ds_targ) **2).mean(metric_dim)   
                denom =   ((ds_obs - ds_obs_clim) **2).mean(metric_dim)          
                ds_toplot =  1 - nom/denom
                                 
                
            if mask is None:
                mask = spatial_mask(ds_toplot)
            glbavg = np.round(area_weighted_avg(ds_toplot,
                                               mask=mask,
                                               is_ds=False).values,2)
            
            if central_longitude != 0:
                cb = ax.pcolormesh(add_cyclic_point(ds_toplot).lon,
                                add_cyclic_point(ds_toplot).lat,
                                add_cyclic_point(ds_toplot),
                                cmap=plt.cm.get_cmap(cmap),
                                vmax=vmax,
                                vmin=vmin,
                                rasterized=True,
                                transform=ccrs.PlateCarree())
            else:
                cb = ax.pcolormesh(ds_toplot.lon,
                                ds_toplot.lat,
                                ds_toplot,
                                cmap=plt.cm.get_cmap(cmap),
                                vmax=vmax,
                                vmin=vmin,
                                rasterized=True,
                                transform=ccrs.PlateCarree())
            
            _ = ax.coastlines()

            ax.set_ylabel('')
            ax.set_xlabel('')

            ax.set_title(f'{str(glbavg)}',
                         fontsize=20)
            if i == 0:
                ax.set_title(ds + f'\n {str(glbavg)}',
                             fontsize=20)
                
            if ind ==0:
                ax.text(-0.25,
                        1.1,
                        season,
                        fontsize=20,
                        transform=ax.transAxes)
        i = i + 1
    
    divider = make_axes_locatable(ax)
    
    ax_cb   = divider.append_axes('bottom',
                                  size="10%",
                                  pad=0.1,
                                  axes_class=plt.Axes)
    
    cbar    = plt.colorbar(cb,
                           cax=ax_cb,
                           orientation="horizontal")
    cbar.set_label(label=label,
                   size=20)
    cbar.ax.tick_params(labelsize=20)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.1,
                        hspace=0.3)

    if save:
        Path(dir_name).mkdir(parents=True,
                             exist_ok=True)
        plt.savefig(f'{dir_name}/{file_name}.png')
        

def plot_single_map_wmo(ds,
                        ds_sig=None,
                        lat=None,
                        lon=None,
                        polar_stereo=False, 
                        central_longitude=180,
                        lat_lims=[50,90],
                        gridlines=False,
                        cmap=mpl.cm.RdYlBu,
                        vmin=-1.5,
                        vmax=1.5,
                        vals=None,
                        nvals=10,
                        cbar=False,
                        cbar_label='',
                        ticks_rotation=0,
                        title=None,
                        show_mean = True,
                        fnt_size=12,
                        figsize=None,
                        fig_dir=None,
                        fig_name=None,
                        show=False,
                        save=False,
                        **kwargs): 
    
    import cartopy.feature as cfeature    
    
    mpl.rcParams.update({'font.size': fnt_size})

    
    if central_longitude == 0: # remove white line at central longitude
        ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
        ds = ds.sortby(ds.lon)
        if ds_sig is not None:
            ds_sig.coords['lon'] = (ds_sig.coords['lon'] + 180) % 360 - 180
            ds_sig = ds_sig.sortby(ds_sig.lon)

            
    if lat is None:
        lat = ds.lat
    if lon is None:
        lon = ds.lon

    img_extent = [lon[0], lon[-1], lat[0], lat[-1]]
        
    crs = ccrs.PlateCarree(central_longitude=central_longitude)    
    if polar_stereo:
        crs = ccrs.NorthPolarStereo()    
        if max(lat_lims) < 0:
            crs = ccrs.SouthPolarStereo()    
        
    
    fig, ax = plt.subplots(nrows=1,
                           ncols=1, 
                           figsize=figsize, 
                           subplot_kw={'projection' : crs})                           
    
    if vals is not None:
        nvals = len(vals) - 1
        
    if vals is None:
        scale = (vmax-vmin)/float(nvals)
        vals = vmin + (vmax-vmin)*np.arange(nvals+1)/float(nvals)
    
    # norm = mpl.colors.BoundaryNorm(vals, cmap.N)
    norm = mpl.colors.BoundaryNorm(vals, plt.cm.get_cmap(cmap).N)
    
    axis = ax
    if gridlines:
        axis.gridlines(draw_labels=False)
    if polar_stereo:
        polarCentral_set_latlim(lat_lims, axis)
        axis.add_feature(cfeature.NaturalEarthFeature('physical',
                                                      'land', 
                                                      '50m',
                                                      # edgecolor='face',
                                                      facecolor='grey'))
    im = axis.imshow(ds, 
                     origin='lower',
                     extent=img_extent,
                     # cmap=cmap,
                     cmap=plt.cm.get_cmap(cmap),
                     norm=norm,
                     interpolation='none',                     
                     transform=ccrs.PlateCarree())
    if show_mean:
        title = title + f' ({np.round(area_weighted_avg(ds).values,2)})'
    im.set_clim(vmin,
                vmax)
    axis.coastlines()
    axis.set_title(title,
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
        if polar_stereo:
            clb_y = 0.0
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
        # tick_locator = ticker.MaxNLocator(nbins=nvals/2)
        # cb.locator = tick_locator
        # cb.update_ticks()

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
        
    if show is False:
        plt.close()


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
            contour_levels = [0,0.5, 1,2,4,6, 8,10,15,20,25,30,35,40,50],
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
                        levels = contour_levels,
                        cmap = cmap
                        # vmax = 30,
                        # vmin = 0,
                        # ds_dict[ds]['linestyle'],
                            )
            contours = ax_.contour(xx, ts.lev.values, ts, colors='white', levels=contour_f.levels)
            ax_.clabel(contours, inline=True, fontsize=10, colors='white')
            ax_.set_xticks(np.unique(np.floor(xx)), np.unique(np.floor(xx)).astype(int), rotation  = 45)
                        # color=ds_dict[ds]['color'])
            # ax[ind].colorbar()
                
            ax_.set_title(title + ' - ' +  ds +' - ' + bms_label)
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

        fig.colorbar(contour_f , ax=ax, orientation='vertical', label = colorbar_label)

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
                            yeartoplot,
                            quiver_dicts = None, 
                            quiver_axis = None,
            states_dict = {1997 : 'model negative', 2001: 'model positive', 1987 : 'real positive', 2000: 'real negeative'},
            ldyr=1-1,
            title='',
            figsize=(10,45),
            contour_levels = [0,0.5, 1,2,4,6, 8,10,15,20,25,30,35,40,50],
            cmap = 'viridis',
            dir_name=None,
            file_name=None,
            colorbar_label = None, 
            season = 'ANN',
            lev_interp = None,
            lev_range = None,
            show=False,
            return_fig_handles = False,
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

    if show:
        ldyr1 = ldyr + 1
    

        fig, ax = plt.subplots(len(ds_list), 1 ,figsize=figsize)
        ds_dict = ds_dicts[bms_label]
        if quiver_dicts is not None:
            qv_dict = quiver_dicts[bms_label]

        for ind, ds in enumerate(ds_list):
            # if len(ds.split('-')) == 1:
            dstp = ds_dict[ds]

            if quiver_dicts is not None:
                qvtp = qv_dict[ds]
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

                except:
                    ts = dstp.sel(time=slice(ldyr*12,(ldyr1)*12-1))

                    if quiver_dicts is not None:
                        qv_ts = qvtp.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                                                        
                    print('DJFy unseccessfull, make sure the conversion is already done!')
            else:
                ts = dstp.sel(time=slice(ldyr*12, (ldyr1)*12-1))
                if quiver_dicts is not None:
                    qv_ts = qvtp.sel(time=slice(ldyr*12, (ldyr1)*12-1))
                                                       
            
            if all([type(yeartoplot) == str, 'diff' not in yeartoplot]): 
            
                y1 = eval(yeartoplot.split('-')[0])
                y0 = eval(yeartoplot.split('-')[1])
                ts = (ts.sel(time = seasons[season] ).sel(year = y1 ) - ts.sel(time = seasons[season] ).sel(year = y0)).mean('time').squeeze()

                if quiver_dicts is not None:
                    qv_ts = (qv_ts.sel(time = seasons[season] ).sel(year = y1 ) - qv_ts.sel(time = seasons[season] ).sel(year = y0)).mean('time').squeeze()

            else:
                ts = ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean('year').mean('time').squeeze()

                if quiver_dicts is not None:
                    qv_ts = qv_ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean('year').mean('time').squeeze()

            xx = ts.lon.values
            
            if lev_interp is not None:
                ts = ts.interp(lev = lev_interp )
                if quiver_dicts is not None:
                    qv_ts = qv_ts.interp(lev = lev_interp )
            
            if lev_range is not None:
                if type(lev_range) is not list:
                    ls = [0, lev_range]
                    lev_range = ls
                    del ls
            
                ts = ts.where((ts.lev <= lev_range[1]) & (ts.lev >= lev_range[0]), drop = True)
                if quiver_dicts is not None:
                    qv_ts = qv_ts.where((qv_ts.lev <= lev_range[1]) & (qv_ts.lev >= lev_range[0]), drop = True)


            label = f'{ds}'
            if len(ds_list) == 1:
                ax_ = ax
            else:
                ax_ = ax[ind]
            
            contour_f = ax_.contourf(xx, ts.lev.values,
                        ts,
                        levels = contour_levels,
                        cmap = cmap
                        # vmax = 30,
                        # vmin = 0,
                        # ds_dict[ds]['linestyle'],
                            )
            contours = ax_.contour(xx, ts.lev.values, ts, colors='white', levels=contour_f.levels)
            ax_.clabel(contours, inline=True, fontsize=10, colors='white')

            if quiver_dicts is not None:
                if quiver_axis == 'U':
                    U = qv_ts.values
                    V = np.zeros_like(qv_ts.values)
                else:
                    V = qv_ts.values
                    U = np.zeros_like(qv_ts.values)
                ax_.quiver(xx, qv_ts.lev.values, np.nan_to_num(U, nan=0.0), np.nan_to_num(V, nan=0.0),  alpha = 0.5)

            
            if  type(yeartoplot) == str:
                ax_.set_title(title + ' - ' +  ds +' - ' + bms_label + f' {yeartoplot}' + f' {season}')
            
            elif all([len(yeartoplot) >1 , type(yeartoplot) == list]):
                ax_.set_title(title + ' - ' +  ds +' - ' + bms_label + f' composite' + f' {season}')

            else:
                if yeartoplot[0] not in states_dict.keys():
                    states_dict[yeartoplot[0]] = 'NA'
                ax_.set_title(title + ' - ' +  ds +' - ' + bms_label + f' {yeartoplot[0]}' + f' {season}' + f': ENSO {states_dict[yeartoplot[0]]}')
            ax_.set_ylabel('depth (m)')
            ax_.set_xlabel('Lon (deg East)')
            ax_.invert_yaxis()


        fig.colorbar(contour_f , ax=ax, orientation='vertical', label = colorbar_label)

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
                            yeartoplot,
                            quiver_dicts = None, 
            states_dict = {1997 : 'model negative', 2001: 'model positive', 1987 : 'real positive', 2000: 'real negeative'},
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
            show=False,
            return_fig_handles = False,
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

    if show:
        ldyr1 = ldyr + 1
    

        fig, ax = plt.subplots(len(ds_list), 1 ,figsize=figsize)
        ds_dict = ds_dicts[bms_label]

        if quiver_dicts is not None:
            qv_dict = quiver_dicts[bms_label]

        for ind, ds in enumerate(ds_list):
            # if len(ds.split('-')) == 1:
            dstp = ds_dict[ds]
            dstp_lev = True if 'lev' in dstp.dims else False 
            if quiver_dicts is not None:
                qvtp = qv_dict[ds]
                qvtp_lev = True if 'lev' in qvtp.dims else False 
            # else:
            #     ds1 = (ds.split('-')[0]).split(' ')[0]
            #     ds0 = (ds.split('-')[1]) .split(' ')[1]
            #     dstp = ds_dict[ds1] - ds_dict[ds0]

                # if quiver_dicts is not None:
                #     qvtp = qv_dict[ds1] - qv_dict[ds0]

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
                ts = dstp.sel(time=slice(ldyr*12, (ldyr1)*12-1))

                if quiver_dicts is not None:
                    qv_ts = qvtp.sel(time=slice(ldyr*12, (ldyr1)*12-1))
                                                       

            if all([type(yeartoplot) == str, 'diff' not in yeartoplot]): 
                y1 = eval(yeartoplot.split('-')[0])
                y0 = eval(yeartoplot.split('-')[1])
                ts = (ts.sel(time = seasons[season] ).sel(year = y1 ) - ts.sel(time = seasons[season] ).sel(year = y0)).mean('time').squeeze()

                if quiver_dicts is not None:
                    qv_ts = (qv_ts.sel(time = seasons[season] ).sel(year = y1 ) - qv_ts.sel(time = seasons[season] ).sel(year = y0)).mean('time').squeeze()

            else:
                ts = ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean('year').mean('time').squeeze()

                if quiver_dicts is not None:
                    qv_ts = qv_ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean('year').mean('time').squeeze()

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
            
            if  type(yeartoplot) == str:
                ax_.set_title(title + ' - ' + f'{lev_range[0]}-{lev_range[1]} (m)' + ' - ' +  ds +' - ' + bms_label + f' {yeartoplot}' + f' {season}')
            
            elif all([len(yeartoplot) >1 , type(yeartoplot) == list]):
                ax_.set_title(title + ' - ' + f'{lev_range[0]}-{lev_range[1]} (m)' + ' - ' +  ds +' - ' + bms_label + f' composite' + f' {season}')

            else:
                if yeartoplot[0] not in states_dict.keys():
                    states_dict[yeartoplot[0]] = 'NA'
                ax_.set_title(title + ' - ' + f'{lev_range[0]}-{lev_range[1]} (m)' + ' - ' +  ds +' - ' + bms_label + f' {yeartoplot[0]}' + f' {season}' + f': ENSO {states_dict[yeartoplot[0]]}')
            ax_.set_ylabel('lat (deg North)')
            ax_.set_xlabel('Lon (deg East)')
            # ax_.invert_yaxis()


        fig.colorbar(im , ax=ax, orientation='vertical', label = colorbar_label)

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
                            yeartoplot,
                            longitude,
                            quiver_dicts = None, 
                            quiver_axis = 'U',
            states_dict = {1997 : 'model negative', 2001: 'model positive', 1987 : 'real positive', 2000: 'real negeative'},
            ldyr=1-1,
            title='',
            figsize=(10,45),
            contour_levels = [0,0.5, 1,2,4,6, 8,10,15,20,25,30,35,40,50],
            cmap = 'viridis',
            dir_name=None,
            file_name=None,
            colorbar_label = None, 
            season = 'ANN',
            lev_interp = None,
            lev_range:list = None,
            show=False,
            return_fig_handles = False,
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

    if show:
        ldyr1 = ldyr + 1
    

        fig, ax = plt.subplots(len(ds_list), 1 ,figsize=figsize)
        ds_dict = ds_dicts[bms_label]

        if quiver_dicts is not None:
            qv_dict = quiver_dicts[bms_label]

        for ind, ds in enumerate(ds_list):
            # if len(ds.split('-')) == 1:
            dstp = ds_dict[ds]
            if quiver_dicts is not None:
                qvtp = qv_dict[ds]
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
                                                        
                except:
                    ts = dstp.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                    if quiver_dicts is not None:
                          qv_ts = qvtp.sel(time=slice(ldyr*12,(ldyr1)*12-1))   

                    print('DJFy unseccessfull, make sure the conversion is already done!')
            else:
                ts = dstp.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                if quiver_dicts is not None:
                    qv_ts = qvtp.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                                                        

            if all([type(yeartoplot) == str, 'diff' not in yeartoplot]): 
                y1 = eval(yeartoplot.split('-')[0])
                y0 = eval(yeartoplot.split('-')[1])
                ts = (ts.sel(time = seasons[season] ).sel(year = y1 ) - ts.sel(time = seasons[season] ).sel(year = y0)).mean('time').squeeze()
                if quiver_dicts is not None:
                    qv_ts = (qv_ts.sel(time = seasons[season] ).sel(year = y1 ) - qv_ts.sel(time = seasons[season] ).sel(year = y0)).mean('time').squeeze()

            else:
                ts = ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean('year').mean('time').squeeze()
                if quiver_dicts is not None:
                    qv_ts = qv_ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean('year').mean('time').squeeze()

            xx = ts.lat.values
            
            if lev_interp is not None:
                ts = ts.interp(lev = lev_interp )
                if quiver_dicts is not None:
                    qv_ts = qv_ts.interp(lev = lev_interp )
            
            if lev_range is not None:
                if type(lev_range) is not list:
                    ls = [0, lev_range]
                    lev_range = ls
                    del ls                
                ts = ts.where((ts.lev <= lev_range[1]) & (ts.lev >= lev_range[0]), drop = True)
                if quiver_dicts is not None:
                    qv_ts = qv_ts.where((qv_ts.lev <= lev_range[1]) &  (qv_ts.lev >= lev_range[0]), drop = True)
                
            ts = ts.sel(lon = longitude, method = 'nearest')
            if quiver_dicts is not None:
                qv_ts = qv_ts.sel(lon = longitude, method = 'nearest')


            label = f'{ds}'
            if len(ds_list) == 1:
                ax_ = ax
            else:
                ax_ = ax[ind]
            
            contour_f = ax_.contourf(xx, ts.lev.values,
                        ts,
                        levels = contour_levels,
                        cmap = cmap)
            contours = ax_.contour(xx, ts.lev.values, ts, colors='white', levels=contour_f.levels)
            ax_.clabel(contours, inline=True, fontsize=10, colors='white')

            if quiver_dicts is not None:
                if quiver_axis == 'U':
                    U = qv_ts.values
                    V = np.zeros_like(qv_ts)
                else:
                    V = qv_ts.values
                    U = np.zeros_like(qv_ts)
                ax_.quiver(xx, qv_ts.lev.values, np.nan_to_num(U, nan=0.0), np.nan_to_num(V, nan=0.0),  alpha = 0.5)

            
            if  type(yeartoplot) == str:
                ax_.set_title(title + ' - ' +  f'lon: {ts.lon.values} deg east' + ' - ' + ds +' - ' + bms_label + f'\n {yeartoplot}' + f' {season}')
            
            elif all([len(yeartoplot) >1 , type(yeartoplot) == list]):
                ax_.set_title(title + ' - ' +  f'lon: {ts.lon.values} deg east' + ' - ' + ds +' - ' + bms_label + f'\n composite' + f' {season}')

            else:
                if yeartoplot[0] not in states_dict.keys():
                    states_dict[yeartoplot[0]] = 'NA'
                ax_.set_title(title + ' - ' + f'lon: {ts.lon.values} deg east' + ' - ' +  ds +' - ' + bms_label + f'\n {yeartoplot[0]}' + f' {season}' + f': ENSO {states_dict[yeartoplot[0]]}')
            ax_.set_ylabel('depth (m)')
            ax_.set_xlabel('lat (deg)')
            ax_.invert_yaxis()


        fig.colorbar(contour_f , ax=ax, orientation='vertical', label = colorbar_label)

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
                            yeartoplot,
            states_dict = {1997 : 'model negative', 2001: 'model positive', 1987 : 'real positive', 2000: 'real negeative'},
            ldyr=1-1,
            title='',
            figsize=(10,45),
            contour_dict = None,
            contour_levels = [0,0.5, 1,2,4,6, 8,10,15,20,25,30,35,40,50],
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
            headlengthfloat = 2,
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

    print( f'Scale Factor Y : {scale_factor_Y} \n Scale Factor X : {scale_factor_X} \n ====================')
    assert all([quiver_dicts_Y is None, quiver_dicts_X is  None]) is False

    if season in ['DJF','MAM','JJA','SON']:
        assert ldyr == 0
        for ds in ds_list:
            assert 'hindcast' not in ds

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
                tsX = dstpX.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                tsY = dstpY.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                if contour_dict is not None:
                    ts = dstp.sel(time=slice(ldyr*12,(ldyr1)*12-1))                                                       
            
            if all([type(yeartoplot) == str, 'diff' not in yeartoplot]): 
            
                y1 = eval(yeartoplot.split('-')[0])
                y0 = eval(yeartoplot.split('-')[1])
                tsX = (tsX.sel(time = seasons[season] ).sel(year = y1 ) - tsX.sel(time = seasons[season] ).sel(year = y0)).mean('time').squeeze()
                tsY = (tsY.sel(time = seasons[season] ).sel(year = y1 ) - tsY.sel(time = seasons[season] ).sel(year = y0)).mean('time').squeeze()
                if contour_dict is not None:
                    ts = (ts.sel(time = seasons[season] ).sel(year = y1 ) - ts.sel(time = seasons[season] ).sel(year = y0)).mean('time').squeeze()


            else:
                tsX = tsX.sel(time = seasons[season] ).sel(year = yeartoplot ).mean('year').mean('time').squeeze()
                tsY = tsY.sel(time = seasons[season] ).sel(year = yeartoplot ).mean('year').mean('time').squeeze()
                if contour_dict is not None:
                    ts = ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean('year').mean('time').squeeze()

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
                            levels = contour_levels,
                            cmap = cmap
                            # vmax = 30,
                            # vmin = 0,
                            # ds_dict[ds]['linestyle'],
                                )
                contours = ax_.contour(xx, ts.lev.values, ts, colors='white', levels=contour_f.levels)
                ax_.clabel(contours, inline=True, fontsize=10, colors='white')

            ax_.quiver(xx, tsX.lev.values, np.nan_to_num(tsX, nan=0.0), np.nan_to_num(tsY, nan=0.0),  alpha = 0.5, scale = arrow_scale)
            
            if  type(yeartoplot) == str:
                ax_.set_title(title + ' - ' +  ds +' - ' + bms_label + f' {yeartoplot}' + f' {season}')
            
            elif all([len(yeartoplot) >1 , type(yeartoplot) == list]):
                ax_.set_title(title + ' - ' +  ds +' - ' + bms_label + f' composite' + f' {season}')

            else:
                if yeartoplot[0] not in states_dict.keys():
                    states_dict[yeartoplot[0]] = 'NA'
                ax_.set_title(title + ' - ' +  ds +' - ' + bms_label + f' {yeartoplot[0]}' + f' {season}' + f': ENSO {states_dict[yeartoplot[0]]}')
            ax_.set_ylabel('depth (m)')
            ax_.set_xlabel('Lon (deg East)')
            ax_.set_ylim(lev_range)
            ax_.invert_yaxis()

        if contour_dict is not None:
            fig.colorbar(contour_f , ax=ax, orientation='vertical', label = colorbar_label)
    
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
                            yeartoplot,
                            longitude, 
            states_dict = {1997 : 'model negative', 2001: 'model positive', 1987 : 'real positive', 2000: 'real negeative'},
            ldyr=1-1,
            title='',
            figsize=(10,45),
            contour_dict = None,
            contour_levels = [0,0.5, 1,2,4,6, 8,10,15,20,25,30,35,40,50],
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

    print( f'Scale Factor Y : {scale_factor_Y} \n Scale Factor X : {scale_factor_X} \n ====================')

    assert all([quiver_dicts_Y is None, quiver_dicts_X is  None]) is False

    if season in ['DJF','MAM','JJA','SON']:
        assert ldyr == 0
        for ds in ds_list:
            assert 'hindcast' not in ds

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
                tsX = dstpX.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                tsY = dstpY.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                if contour_dict is not None:
                    ts = dstp.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                                   

            if all([type(yeartoplot) == str, 'diff' not in yeartoplot]): 
                y1 = eval(yeartoplot.split('-')[0])
                y0 = eval(yeartoplot.split('-')[1])
                tsX = (tsX.sel(time = seasons[season] ).sel(year = y1 ) - tsX.sel(time = seasons[season] ).sel(year = y0)).mean('time').squeeze()
                tsY = (tsY.sel(time = seasons[season] ).sel(year = y1 ) - tsY.sel(time = seasons[season] ).sel(year = y0)).mean('time').squeeze()
                if contour_dict is not None:
                    ts = (ts.sel(time = seasons[season] ).sel(year = y1 ) - ts.sel(time = seasons[season] ).sel(year = y0)).mean('time').squeeze()

            else:
                tsX = tsX.sel(time = seasons[season] ).sel(year = yeartoplot ).mean('year').mean('time').squeeze()
                tsY = tsY.sel(time = seasons[season] ).sel(year = yeartoplot ).mean('year').mean('time').squeeze()
                if contour_dict is not None:
                    ts = ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean('year').mean('time').squeeze()

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
                            levels = contour_levels,
                            cmap = cmap)
                contours = ax_.contour(xx, ts.lev.values, ts, colors='white', levels=contour_f.levels)
                ax_.clabel(contours, inline=True, fontsize=10, colors='white')
            ax_.quiver(xx, tsX.lev.values, np.nan_to_num(tsX.values, nan=0.0), np.nan_to_num(tsY.values, nan=0.0), scale = arrow_scale)


            if  type(yeartoplot) == str:
                ax_.set_title(title + ' - ' +  f'lon: {tsX.lon.values} deg east' + ' - ' + ds +' - ' + bms_label + f'\n {yeartoplot}' + f' {season}')
            
            elif all([len(yeartoplot) >1 , type(yeartoplot) == list]):
                ax_.set_title(title + ' - ' +  f'lon: {tsX.lon.values} deg east' + ' - ' + ds +' - ' + bms_label + f'\n composite' + f' {season}')

            else:
                if yeartoplot[0] not in states_dict.keys():
                    states_dict[yeartoplot[0]] = 'NA'
                ax_.set_title(title + ' - ' + f'lon: {tsX.lon.values} deg east' + ' - ' +  ds +' - ' + bms_label + f'\n {yeartoplot[0]}' + f' {season}' + f': ENSO {states_dict[yeartoplot[0]]}')
            ax_.set_ylabel('depth (m)')
            ax_.set_xlabel('lat (deg)')
            ax_.set_ylim(lev_range)
            ax_.invert_yaxis()
        if contour_dict is not None:
            fig.colorbar(contour_f , ax=ax, orientation='vertical', label = colorbar_label)

    
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
                            yeartoplot,
                            quiver_dicts = None, 
            states_dict = {1997 : 'model negative', 2001: 'model positive', 1987 : 'real positive', 2000: 'real negeative'},
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
            headwidth = 1,
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

    print( f'Scale Factor Y : {scale_factor_Y} \n Scale Factor X : {scale_factor_X} \n ====================')
    assert all([quiver_dicts_Y is None, quiver_dicts_X is  None]) is False

    if season in ['DJF','MAM','JJA','SON']:
        assert ldyr == 0
        for ds in ds_list:
            assert 'hindcast' not in ds

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
                tsX = dstpX.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                tsY = dstpY.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                if colorbar_dict is not None:
                    ts = dstp.sel(time=slice(ldyr*12,(ldyr1)*12-1))
                                                       

            if all([type(yeartoplot) == str, 'diff' not in yeartoplot]): 
                y1 = eval(yeartoplot.split('-')[0])
                y0 = eval(yeartoplot.split('-')[1])
                tsX = (tsX.sel(time = seasons[season] ).sel(year = y1 ) - tsX.sel(time = seasons[season] ).sel(year = y0)).mean('time').squeeze()
                tsY = (tsY.sel(time = seasons[season] ).sel(year = y1 ) - tsY.sel(time = seasons[season] ).sel(year = y0)).mean('time').squeeze()
                if colorbar_dict is not None:
                    ts = (ts.sel(time = seasons[season] ).sel(year = y1 ) - ts.sel(time = seasons[season] ).sel(year = y0)).mean('time').squeeze()

            else:
                tsX = tsX.sel(time = seasons[season] ).sel(year = yeartoplot ).mean('year').mean('time').squeeze()
                tsY = tsY.sel(time = seasons[season] ).sel(year = yeartoplot ).mean('year').mean('time').squeeze()
                if colorbar_dict is not None:
                    ts = ts.sel(time = seasons[season] ).sel(year = yeartoplot ).mean('year').mean('time').squeeze()

            xx = ts.lon.values
            
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
                
            ax_.quiver(xx, tsX.lat.values, np.nan_to_num(tsX.values, nan=0.0), np.nan_to_num(tsY.values, nan=0.0), scale = arrow_scale, headwidth = headwidth)

            
            if  type(yeartoplot) == str:
                ax_.set_title(title + ' - ' + f'{lev_range[0]}-{lev_range[1]} (m)' + ' - ' +  ds +' - ' + bms_label + f' {yeartoplot}' + f' {season}')
            
            elif all([len(yeartoplot) >1 , type(yeartoplot) == list]):
                ax_.set_title(title + ' - ' + f'{lev_range[0]}-{lev_range[1]} (m)' + ' - ' +  ds +' - ' + bms_label + f' composite' + f' {season}')

            else:
                if yeartoplot[0] not in states_dict.keys():
                    states_dict[yeartoplot[0]] = 'NA'
                ax_.set_title(title + ' - ' + f'{lev_range[0]}-{lev_range[1]} (m)' + ' - ' +  ds +' - ' + bms_label + f' {yeartoplot[0]}' + f' {season}' + f': ENSO {states_dict[yeartoplot[0]]}')
            ax_.set_ylabel('lat (deg North)')
            ax_.set_xlabel('Lon (deg East)')
            # ax_.invert_yaxis()

        if colorbar_dict is not None:
            fig.colorbar(im , ax=ax, orientation='vertical', label = colorbar_label)

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

def DJFy(ds):
    ds_long = ds.sel(time = np.arange(12)).stack(month = ('year','time')).transpose('month',...)
    ds_shifted = xr.full_like(ds_long, np.nan).transpose('month',...)
    ds_shifted[1:,] = ds_long[:-1,].values
    return ds_shifted.unstack().transpose(*ds.dims)



