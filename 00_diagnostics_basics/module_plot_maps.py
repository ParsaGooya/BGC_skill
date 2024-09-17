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
                    cbar_label=r'mol m$^{-2}$ yr$^{-1}$',  
                    std  = False,        
                    individual_months = False,           
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
        
        if specific_years is not None:
            obs_ref = data_dict['obs'].sel(year = specific_years)
        else:
            obs_ref = data_dict['obs']
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
            if specific_years is not None:
                ds_toplot = data_dict[ds].sel(year = specific_years)
            else:
                ds_toplot = data_dict[ds]
            
            if std:
                ds_toplot = ds_toplot.sel(time=inds)
                ds_toplot = ds_toplot.mean(['time']).std('year')  


            else:
                ds_toplot = ds_toplot.sel(time=inds).mean(['year',
                                                                'time'])
            
            cb = ax.pcolormesh(ds_toplot.lon,
                               ds_toplot.lat,
                               ds_toplot,
                               cmap=plt.cm.get_cmap('RdBu_r'),
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
            
            corr_pat      = corr_map(data_dict['obs'].sel(time=inds).mean(['time']),
                                     data_dict[ds].sel(time=inds).mean(['time']))
            corr_pat_avg  = np.round(corr_pat.values.mean(),2)

            ax.set_title(f'{str(glbavg)}, {str(corr_pat_avg)}',
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
                  individual_months = True,
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
            ds_obs = data_dict['obs'].sel(time=inds)
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


def trend(ds, dim = 'year', return_detrended = False ):
    m = ds.polyfit( dim  = dim, deg = 1).polyfit_coefficients
    out = m[0]*ds[dim] +  m[1]
    if return_detrended:
        return  out, ds - out
    else:
        return out