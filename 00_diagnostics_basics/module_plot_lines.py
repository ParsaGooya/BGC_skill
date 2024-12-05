import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import xarray as xr

def plot_ts_glbavg(ds_list,
            ds_dict,
            xx=None,
            nx=1,
            dict_label='ts',
            ref='obs',
            err=False,
            title='',
            bbox=(.68,.5,.5,.5),
            figsize=(6,3),
            dir_name=None,
            file_name=None,
            ylabel = None,
            show=False,
            save=False):
    '''
     to plot data written in terms of initial year
    '''
    if show:
        xx1 = ds_dict[ds_list[0]][dict_label].year.values
        xx2 = xx1[-1] + 1
        for ix in np.arange(nx):
            plt.figure(figsize=figsize)
            for ind, ds in enumerate(ds_list):
                
                ts = ds_dict[ds][dict_label].sel(time=slice(ix*12,
                                                       (ix+1)*12-1)).mean(dim='time')
                if err:
                    ts = abs(ts - ds_dict[ref][dict_label].sel(time=slice(ix*12,
                                                                          (ix+1)*12-1)).mean(dim='time'))
                xx = xx1
                if ts.year.size == 1:
                    xx = xx2
                    
                plt.plot(ix + xx,
                         ts,
                         ds_dict[ds]['linestyle'],
                         label=f'{ds}',
                         color=ds_dict[ds]['color'])
                
                plt.title(title)

                
            plt.legend(loc='best',
                       bbox_to_anchor=bbox,
                       fontsize='small',
                       handlelength=2,
                       frameon=False)   
                
            plt.ylabel(ylabel)   
        if save:
            Path(dir_name).mkdir(parents=True,
                                exist_ok=True)
            plt.savefig(f'{dir_name}/{file_name}.png',
                        bbox_inches='tight',
                        dpi=300)   

            
def plot_ts_glbavg_on_target(ds_list,
                            ds_dict,
            xx=None,
            ldyr=1-1,
            dict_label='ts',
            ref='obs',
            err=False,
            title='',
            bbox=(.68,.5,.5,.5),
            figsize=(6,3),
            dir_name=None,
            file_name=None,
            ylabel = None, 
            season = 'ANN',
            monthly_res = False,
            Correlations = False,
            Trend = False,
            ENSO = True,
            show=False,
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
            'ANN': np.arange(12) + 12*ldyr}
    
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
        
        plt.figure(figsize=figsize)
        if season == 'DJF':
            ref = DJFy(ds_dict['obs'][dict_label]).sel(time=slice(ldyr*12,
                                                        (ldyr1)*12-1))
        else:
            ref = ds_dict['obs'][dict_label].sel(time=slice(ldyr*12,
                                                        (ldyr1)*12-1))

        ref = ref.sel(time = seasons[season] )

        if monthly_res:
            ref =  ref.stack(yearmonth = ('year','time'))
            ref['yearmonth'] = ref.year.values + (ref.time.values + 0.5)/12
        else:
            ref =  ref.mean('time')

        for ind, ds in enumerate(ds_list):
            if season == 'DJF':
                ts = ds_dict[ds][dict_label].sel(time=slice(ldyr*12,
                                                        (ldyr1)*12-1))
            else:
                ts = ds_dict[ds][dict_label].sel(time=slice(ldyr*12,
                                                        (ldyr1)*12-1))
            ts = ts.sel(time = seasons[season] )

            if monthly_res:
                ts = ts.stack(yearmonth = ('year','time'))
                xx = ts.year.values + (ts.time.values + 0.5)/12
                ts['yearmonth'] = ts.year.values + (ts.time.values + 0.5)/12
                dim = 'yearmonth'
            else:
                ts = ts.mean(dim='time')
                xx = ts.year.values 
                dim = 'year'
            # if err:
            #     y0 = np.max([ds_dict[ref][dict_label].year.values[0],
            #                  ts.year.values[0]])
            #     y1 = np.min([ds_dict[ref][dict_label].year.values[-1],
            #                  ts.year.values[-1]])
            #     xx = np.arange(y0,y1+1)
            #     ts = (ts.sel(year=slice(y0,y1)) - ds_dict[ref][dict_label].sel(time=slice(ldyr*12,
            #                                                                                  ldyr1*12-1)).mean(dim='time').sel(year=slice(y0,y1)))
                
            if Correlations:
                corr = xr.corr(ts, ref, dim = dim).values
                corr_detrend = xr.corr(trend(ts, dim = dim, return_detrended = True)[1], trend(ref, dim = dim, return_detrended = True)[1], dim =dim).values
                label = f'{ds} {np.round(corr,2)} ({np.round(corr_detrend,2)})'
            else:
                label = f'{ds}'
            
            plt.plot(xx,
                         ts,
                         ds_dict[ds]['linestyle'],
                         label=label,
                         color=ds_dict[ds]['color'])

            if Trend:
                ts_trend = trend(ts, dim = dim)
                plt.plot(xx,
                        ts_trend,
                        linestyle = 'dashed',
                        color=ds_dict[ds]['color'])

                
            plt.title(title)
            plt.ylabel(ylabel)
        
        if ENSO:
                ylim = plt.get_ylim()
                
                if monthly_res:
                    ENSO_years = [1982 + 3/12, 1983 + 6/12, 1986 + 8/12, 1988 + 2/12, 1991 + 4/12, 1992 + 6/12, 1994 + 8/12, 1995 + 4/ 12, 1997 + 4/12, 1998 + 5/12 , 
                                      2002 + 5 /12, 2003 + 3/12,  2004 + 6/12, 2005 + 2/12, 2006 + 8/12, 2006 + 1/12, 2009 + 6/12, 2010 + 3/12, 2014 + 9/12, 2016 + 4/12,  
                                      2018 + 8/12, 2019 + 6/12, 2023 + 5/12, 2024  +5/12]

                else:
                    ENSO_years = [1996.5, 1999.5, 2008.5, 2010.5,2014.5, 2016.5, 2022.5, 2024.5]
                for x in ENSO_years:
                    plt.axvline(x=x,linestyle = 'dashed', color = 'r', alpha = 0.5)

                
        plt.legend(loc='best',
                       bbox_to_anchor=bbox,
                       fontsize='small',
                       handlelength=2,
                       frameon=False)   
                
                
    if save:
        if show is False:
            print(f"If show = {show} then save = False --change to show = True")
        else:
            Path(dir_name).mkdir(parents=True,
                                     exist_ok=True)
            plt.savefig(f'{dir_name}/{file_name}.png',
                            bbox_inches='tight',
                            dpi=300)   


def plot_ts_vs_lead(xa_list,
            xa_dict,
            ylim_min=0.,
            ylim_max=0.1,
            xlim_min=0,
            xlim_max=59,
            xticks_step=1,       
            xticks_labels=None,
            color_dict=None,
            title='',
            xlabel='',
            ylabel='',
            labels=None,
            ncol_labels=1,
            bbox=(.68,.5,.5,.5),
            figsize=(6,3),
            fontsize=14,
            show_leg=True,
            dict_label=None,
            dir_name=None,
            file_name=None,
            show=False,
            save=False):
    if show:
        
        if labels is None:
            labels = xa_list
        
        fig, ax = plt.subplots(figsize=figsize)
        
        for ind, xa in enumerate(xa_list):
            if dict_label is None:
                ts = xa_dict[xa]
                ts = ts[np.where(ts != 0.)]
            if dict_label is not None:
                ts = xa_dict[xa][dict_label]
            xx = np.arange(ts.size)
            
            color = color_dict
            if color_dict is not None:
                color = color_dict[xa]#['color']
                
            plt.plot(xx,
                     ts,
                     'o-',
                     markersize=5,
                     color=color,
                     label=labels[ind])
            
        plt.xlim(xlim_min,
                 xlim_max)
        plt.ylim(ylim_min,
                 ylim_max)
        plt.title(title,
                  fontsize=fontsize+2)
        
        plt.xlabel(xlabel,
                   fontsize=fontsize)
        plt.ylabel(ylabel,
                   fontsize=fontsize)
        
            
        # xticks_arr = np.arange(.2+xlim_min,
        #                        xlim_max+.2,
        xticks_arr = np.arange(xlim_min,
                               xlim_max+1,
                               xticks_step)              
        ax.set_xticks(xticks_arr) # <--- set the ticks first
        # ax.set_xticklabels(np.arange(xlim_min_label,
        #                              xlim_max_label,
        #                              xticks_step,
        #                              dtype=int))  
        xticks_labels = ['']+list(np.arange(xlim_min+2,
                                            xlim_max+1,
                                            xticks_step))+['']
        ax.set_xticklabels(xticks_labels)  
        plt.xticks(fontsize=fontsize)        
        plt.yticks(fontsize=fontsize)        
        
        
        if show_leg:
            plt.legend(loc='best',
                       bbox_to_anchor=bbox,
                       # fontsize='small',
                       fontsize=fontsize,
                       handlelength=1,
                       ncol=ncol_labels,
                       frameon=False)   
                
        if save:
            Path(dir_name).mkdir(parents=True,
                                exist_ok=True)
            plt.savefig(f'{dir_name}/{file_name}.png',
                        bbox_inches='tight',
                        dpi=300)                           


def plot_lines(ds,
               dim_x='time',
               title='',
               xlabel='',
               ylabel='',
               ylim_min=0,
               ylim_max=.8,
               xticks_step=1,
               bbox=(.75,.53,.5,.5),
               show_leg=True,
               show_grid=False,
               fig_dir=None,
               fig_name=None,
               save=False):

    fig, axes = plt.subplots(nrows=1,
                             ncols=1,
                             figsize=(6,3))

    for var in ds.data_vars:
        dsx = ds[var].to_dataset()    
        ax = sns.lineplot(dsx,
                          x=dim_x, 
                          y=var,
                          label=var,
                          linewidth=1,
                          marker='o')
        ax.set(xlabel=xlabel)
        ax.set(ylabel=ylabel)
        
        xticks_arr = np.arange(0,
                               ds[dim_x].values.size,
                               xticks_step)              
        ax.set_xticks(xticks_arr) # <--- set the ticks first
        ax.set_xticklabels(np.arange(1,
                                     ds[dim_x].values.size+1,
                                     xticks_step))     
        plt.legend(loc='best',
                       bbox_to_anchor=bbox,
                       fontsize='small',
                       handlelength=1,
                       frameon=False)            
        
        if show_leg is False:
            plt.legend('',
                       frameon=False)
        
        xlim_min = -0.3
        xlim_max = ds[dim_x].values.size-1 + 0.3
        
        a = ax.set(title=title,
                   xlim=(xlim_min,
                         xlim_max),
                   ylim=(ylim_min,
                         ylim_max))
        
        if show_grid:
            plt.grid()
        
    fig.tight_layout()
    if save:
        Path(fig_dir).mkdir(parents=True,
                            exist_ok=True)
        plt.savefig(f'{fig_dir}/{fig_name}.png',
                        bbox_inches='tight',
                        dpi=300)                                                           
            
def plot_ts_vs_lead_biomes(xa_list,
            xa_dict,
            bms_labels,
            mask_biomes = None,
            ylim_min=None,
            ylim_max=None,
            xlim_min=0,
            xlim_max=59,
            xticks_step=1,       
            xticks_labels=None,
            color_dict=None,
            title='',
            xlabel='',
            ylabel='',
            labels=None,
            ncol_labels=1,
            bbox=(.68,.5,.5,.5),
            figsize=(10,30),
            fontsize=10,
            show_leg=True,
            dict_label=None,
            dir_name=None,
            file_name=None,
            return_fig_handles = False,
            show=False,
            save=False):
    if show:
        
        if labels is None:
            labels = xa_list
        if mask_biomes is not None:
            fig, ax = plt.subplots(len(bms_labels), 2 ,figsize=figsize, gridspec_kw={'width_ratios': [2,0.75]})
        else:
            fig, ax = plt.subplots(len(bms_labels), 1 ,figsize=figsize)
        for ii,bm_label in enumerate(bms_labels):
            if mask_biomes is not None:
                ax_ = ax[(ii, 0)]
            else:
                try:
                    ax_ = ax[ii]
                except:
                    ax_ = ax
 
            if ii<17:
                for ind, xa in enumerate(xa_list):
                    if dict_label is None:
                        ts = xa_dict[bm_label][xa]
                        ts = ts[np.where(ts != 0.)]
                    if dict_label is not None:
                        ts = xa_dict[xa][dict_label]
                    xx = np.arange(ts.size)
                    
                    color = color_dict
                    if color_dict is not None:
                        color = color_dict[xa]#['color']
                        
                    ax_.plot(xx,
                            ts,
                            'o-',
                            markersize=5,
                            color=color,
                            label=labels[ind])
                    
                ax_.set_xlim(xlim_min,
                        xlim_max)
                ax_.set_ylim(ylim_min,
                        ylim_max)
                ax_.set_title(title + '-' + bm_label,
                        fontsize=fontsize+2)
                
                ax_.set_xlabel('')
                ax_.set_ylabel(ylabel,
                        fontsize=fontsize)
                
                ax_.set_xticks([])
                ax_.set_xticklabels([])

                if ii == len(bms_labels) - 1:
                    ax_.set_xlabel(xlabel,
                            fontsize=fontsize)
                    # xticks_arr = np.arange(.2+xlim_min,
                    #                        xlim_max+.2,
                    xticks_arr = np.arange(xlim_min,
                                        xlim_max+1,
                                        xticks_step)              
                    ax_.set_xticks(xticks_arr, fontsize=fontsize) # <--- set the ticks first
                    # ax.set_xticklabels(np.arange(xlim_min_label,
                    #                              xlim_max_label,
                    #                              xticks_step,
                    #                              dtype=int))  
                    xticks_labels = ['']+list(np.arange(xlim_min+2,
                                                        xlim_max+1,
                                                        xticks_step))+['']
                    ax_.set_xticklabels(xticks_labels)  
                    # ax[ii].set_xticks(fontsize=fontsize)        
                    # ax[ii].set_yticks(fontsize=fontsize)        
                if mask_biomes is not None:
                    
                    ax[ii,1].pcolormesh(mask_biomes[bm_label].lon, mask_biomes[bm_label].lat, mask_biomes[bm_label].values)
                    
                if show_leg:
                    ax_.legend(loc='best',
                            bbox_to_anchor=bbox,
                            # fontsize='small',
                            fontsize=fontsize,
                            handlelength=1,
                            ncol=ncol_labels,
                            frameon=False)   
        plt.subplots_adjust(wspace=0.55,
                        hspace=0.15)        
        if save:
            Path(dir_name).mkdir(parents=True,
                                exist_ok=True)
            plt.savefig(f'{dir_name}/{file_name}.png',
                        bbox_inches='tight',
                        dpi=300)                           
        if return_fig_handles:
            return fig, ax


            
            
def plot_ts_biomeavg_on_target(ds_list,
                            ds_dicts,
                            bms_labels,
            mask_biomes = None,
            xx=None,
            ldyr=1-1,
            dict_label='ts',
            ref_ds='obs',
            err=False,
            title='',
            bbox=(.68,.5,.5,.5),
            figsize=(10,45),
            dir_name=None,
            file_name=None,
            ylabel = None, 
            season = 'ANN',
            lev_range = None,
            monthly_res = False,
            Correlations = False,
            Trend = False,
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
    
        if mask_biomes is not None:
            fig, ax = plt.subplots(len(bms_labels), 2 ,figsize=figsize, gridspec_kw={'width_ratios': [2,0.75]})
        else:
            fig, ax = plt.subplots(len(bms_labels), 1 ,figsize=figsize)
        for ii,bm_label in enumerate(bms_labels):
            if mask_biomes is not None:
                ax_ = ax[(ii, 0)]
            else:
                try:
                    ax_ = ax[ii]
                except:
                    ax_ = ax
            if ii<16:
                ds_dict = ds_dicts[bm_label]
                if ref_ds is not None:
                    if season == 'DJF':
                        ref = DJFy(ds_dicts[bm_label][ref_ds][dict_label]).sel(time=slice(ldyr*12,
                                                                    (ldyr1)*12-1))
                    else:
                        ref = ds_dicts[bm_label][ref_ds][dict_label].sel(time=slice(ldyr*12,
                                                                    (ldyr1)*12-1))

                    ref = ref.sel(time = seasons[season] )

                    if 'lev' in ref.dims:
                        if lev_range is not None:
                            ref = ref.where(ref.lev <= lev_range, drop = True)
                        ref = ref.mean('lev')
                                                                    
                    if monthly_res:
                        ref =  ref.stack(yearmonth = ('year','time'))
                        ref['yearmonth'] = ref.year.values + (ref.time.values + 0.5)/12
                    else:
                        ref =  ref.mean('time')
                else:
                    assert Correlations == False, 'provide reference dataset for correlations.'

                for ind, ds in enumerate(ds_list):
                    if season == 'DJF':
                        ts = DJFy(ds_dict[ds][dict_label]).sel(time=slice(ldyr*12,
                                                                (ldyr1)*12-1))
                    else:
                        ts = ds_dict[ds][dict_label].sel(time=slice(ldyr*12,
                                                                (ldyr1)*12-1))
                    ts = ts.sel(time = seasons[season] )
                    if 'lev' in ts.dims:
                        if lev_range is not None:
                            ts = ts.where(ts.lev <= lev_range, drop = True)
                        ts = ts.mean('lev')

                    if monthly_res:
                        ts = ts.stack(yearmonth = ('year','time'))
                        xx = ts.year.values + (ts.time.values + 0.5)/12
                        ts['yearmonth'] = ts.year.values + (ts.time.values + 0.5)/12
                        dim = 'yearmonth'
                    else:
                        ts = ts.mean('time')
                        xx = ds_dict[ds_list[ind]][dict_label].year.values     
                        dim = 'year' 
                    # if err:
                    #     y0 = np.max([ds_dict[ref][dict_label].year.values[0],
                    #                 ts.year.values[0]])
                    #     y1 = np.min([ds_dict[ref][dict_label].year.values[-1],
                    #                 ts.year.values[-1]])
                    #     xx = np.arange(y0,y1+1)
                    #     ts = (ts.sel(year=slice(y0,y1)) - ds_dict[ref][dict_label].sel(time=slice(ldyr*12,
                    #                                                                                 ldyr1*12-1)).mean(dim='time').sel(year=slice(y0,y1)))

                    if Correlations:
                        
                        corr = xr.corr(ts, ref, dim = dim).values
                        corr_detrend = xr.corr(trend(ts, dim = dim, return_detrended = True)[1], trend(ref, dim = dim, return_detrended = True)[1], dim =dim).values
                        label = f'{ds} {np.round(corr,2)} ({np.round(corr_detrend,2)})'
                    else:
                        label = f'{ds}'

                    ax_.plot(xx,
                                ts,
                                ds_dict[ds]['linestyle'],
                                label=label,
                                color=ds_dict[ds]['color'])
                    if Trend:
                        
                        ts_trend = trend(ts, dim = dim)
                        ax_.plot(xx,
                                ts_trend,
                                linestyle = 'dashed',
                                color=ds_dict[ds]['color'])

                        
                    ax_.set_title(title + '-' + bm_label)
                    ax_.set_ylabel(ylabel)
                if ii <len(bms_labels) - 1:
                    ax_.set_xlabel('')
                    ax_.set_xticks([])
                    ax_.set_xticklabels([])

                if ENSO_years is not None:
                        ylim = ax_.get_ylim()
                        if ENSO_years == 'obs':
                            if monthly_res:
                                ENSO_years =  [1982 + 3/12, 1983 + 6/12, 1986 + 8/12, 1988 + 2/12, 1991 + 4/12, 1992 + 6/12, 1994 + 8/12, 1995 + 4/ 12, 1997 + 4/12, 1998 + 5/12 , 
                                      2002 + 5 /12, 2003 + 3/12,  2004 + 6/12, 2005 + 2/12, 2006 + 8/12, 2006 + 1/12, 2009 + 6/12, 2010 + 3/12, 2014 + 9/12, 2016 + 4/12,  
                                      2018 + 8/12, 2019 + 6/12, 2023 + 5/12, 2024  +5/12]
                            else:
                                ENSO_years = [1996.5, 1999.5, 2008.5, 2010.5,2014.5, 2016.5, 2022.5, 2024.5]
                        
                        for x in ENSO_years:
                            ax_.axvline(x=x,linestyle = 'dashed', color = 'g', alpha = 0.25)

                if mask_biomes is not None:

                    ax[ii,1].pcolormesh(mask_biomes[bm_label].lon, mask_biomes[bm_label].lat, mask_biomes[bm_label].values)
                        
                ax_.legend(loc='best',
                            bbox_to_anchor=bbox,
                            fontsize='small',
                            handlelength=2,
                            frameon=False)   
                if return_fig_handles:
                    return fig, ax

    
    plt.subplots_adjust(wspace=0.55,
                        hspace=0.25)             
                
    if save:
        if show is False:
            print(f"If show = {show} then save = False --change to show = True")
        else:
            Path(dir_name).mkdir(parents=True,
                                     exist_ok=True)
            plt.savefig(f'{dir_name}/{file_name}.png',
                            bbox_inches='tight',
                            dpi=300)   




def trend(ds, dim = 'year', return_detrended = False , remove_intercept = False):
    if dim == 'ref':
        ref_ds = xr.full_like(ds, np.nan)
        ds_copy = ds.copy()
        ds_copy['ref'] = np.arange(len(ds.ref))
    else:
        ds_copy = ds.copy()

    m = ds_copy.polyfit( dim  = dim, deg = 1).polyfit_coefficients
    if remove_intercept:
        out = m[0]*ds_copy[dim] 
    else:
        out = m[0]*ds_copy[dim] +  m[1]
    if dim == 'ref':
        ref_ds[:] = out.values
        out = ref_ds
    if return_detrended:
        return  out, ds - out
    else:
        return out

def DJFy(ds):
    ds_long = ds.sel(time = np.arange(12)).stack(month = ('year','time')).transpose('month',...)
    ds_shifted = xr.full_like(ds_long, np.nan).transpose('month',...)
    ds_shifted[1:,] = ds_long[:-1,].values
    return ds_shifted.unstack().transpose(*ds.dims)