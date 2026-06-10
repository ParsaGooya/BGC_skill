import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import random
from matplotlib import colors
import seawater as sw
from modules.analysis.module_data_postprocessing import haversine
def rmse(ds1, ds2):
    return np.sqrt(((ds1 - ds2)**2).mean())

def corr(ds1, ds2):
    mask = ~np.isnan(ds1) & ~np.isnan(ds2)
    ds1_valid = ds1[mask]
    ds2_valid = ds2[mask]
    return np.corrcoef(ds1_valid, ds2_valid)[0, 1]


def scatter_comparison(biome_tp, dataframe, var_y:str, ds_list_var_y: list, units:dict, var_x:str = None, ds_list_var_x: list = ['obs'], 
                       location_ref_var = None, 
                       outlier_var = None, 
                       color_bar_var = None,
                       ds_list_color_bar: list = None,
                       cmap = 'coolwarm',
                       density_contour = False, 
                       depth = None, 
                       month = None,
                        lat = None,
                        lon = None,
                        max_dist = None, 
                       figsize = (16,6), 
                       plot_lim_x = None, 
                       plot_lim_y = None, 
                       vmin = None, 
                       vmax = None, 
                       month_based_scatter = False, 
                       calculate_rmse = True,
                       calculate_r2 = True,
                       add_1_1_line = True,
                       add_zeros_grid = False):
    


    if var_x is None:
        var_x = var_y
    
    if density_contour:
        assert var_y == 'thetao'
        assert var_x == 'so'
  
    if len(ds_list_var_x) > 1:
        assert len(ds_list_var_x) == len(ds_list_var_y)
    
    if color_bar_var is not None:
        assert ds_list_color_bar is not None
        assert month_based_scatter is False
        if len(ds_list_color_bar) > 1:
            assert len(ds_list_color_bar) == len(ds_list_var_y)


    if var_x != var_y:
        assert location_ref_var is not None, 'Specify which variable to be used as the location reference'
    else:
        if plot_lim_y is not None:
            plot_lim_x = plot_lim_y
        elif plot_lim_x is not None:
            plot_lim_y =  plot_lim_x 

    if location_ref_var is None:
        location_ref_var = var_y


    location_ref = dataframe[location_ref_var][biome_tp]
    if outlier_var is not None: 
        location_ref = location_ref[location_ref['var']<outlier_var]

    keys = ['year', 'time', 'lat', 'lon' ,'lev']

    if color_bar_var is not None:
        common = (
            location_ref[keys]
            .merge(dataframe[var_y][biome_tp][keys], on=keys, how="inner")
            .merge(dataframe[var_x][biome_tp][keys], on=keys, how="inner")
            .merge(dataframe[color_bar_var][biome_tp][keys], on=keys, how="inner")
            .drop_duplicates()
        )
        ds_dicts_color = dataframe[color_bar_var][biome_tp].merge(common, on=keys, how="inner") 
    else:
        common = (
            location_ref[keys]
            .merge(dataframe[var_y][biome_tp][keys], on=keys, how="inner")
            .merge(dataframe[var_x][biome_tp][keys], on=keys, how="inner")
            .drop_duplicates()
        )
    
    ds_dicts_y = dataframe[var_y][biome_tp].merge(common, on=keys, how="inner")  # talk measured in the lab, Total scale
    ds_dicts_x = dataframe[var_x][biome_tp].merge(common, on=keys, how="inner")   # DIC measured in the lab in μmol/kg-sw

    location_ref = location_ref.merge(common, on=keys, how="inner")
    # ds_dicts_y = pd.merge(dataframe[var_y][biome_tp], location_ref[['year', 'time', 'lat', 'lon' ,'lev']], on=['year', 'time', 'lat', 'lon' ,'lev'], how='inner')
    # ds_dicts_x = pd.merge(dataframe[var_x][biome_tp], location_ref[['year', 'time', 'lat', 'lon' ,'lev']], on=['year', 'time', 'lat', 'lon' ,'lev'], how='inner')
    # if color_bar_var is not None:
    #     ds_dicts_color = pd.merge(dataframe[color_bar_var][biome_tp], location_ref[['year', 'time', 'lat', 'lon' ,'lev']], on=['year', 'time', 'lat', 'lon' ,'lev'], how='inner')
    

    if var_x == var_y:
        title = f' {biome_tp} - {ds_dicts_y["year"].min()} - {ds_dicts_y["year"].max()} - ' ###
    else:
        title = f' {biome_tp} - {location_ref["year"].min()} - {location_ref["year"].max()} - ' 

    if any([var_x != var_y, location_ref_var != var_y]):
        title = title + f'\n location reference = {location_ref_var} '
    
    if depth is not None:
        if any([type(depth) == tuple , type(depth) == list]):
            ds_dicts_y = ds_dicts_y[(ds_dicts_y['lev'] >= depth[0]) &  (ds_dicts_y['lev'] <= depth[1])]
            ds_dicts_x = ds_dicts_x[(ds_dicts_x['lev'] >= depth[0]) &  (ds_dicts_x['lev'] <= depth[1])]
            if  color_bar_var is not None:
                ds_dicts_color = ds_dicts_color[(ds_dicts_color['lev'] >= depth[0]) &  (ds_dicts_color['lev'] <= depth[1])]
            title = title + f'lev = {depth[0]} - {depth[1]}'
        else:
            ds_dicts_y = ds_dicts_y[ds_dicts_y['lev'] == depth]
            ds_dicts_x = ds_dicts_x[ds_dicts_x['lev'] == depth]
            if color_bar_var is not None:
                ds_dicts_color = ds_dicts_color[ds_dicts_color['lev'] == depth]

            title = title + f' lev = {np.round(depth,2)} '

    if month is not None:
            month_based_scatter = False
            ds_dicts_y = ds_dicts_y[ds_dicts_y['time'] +1  == month ]
            ds_dicts_x = ds_dicts_x[ds_dicts_x['time']  + 1 == month  ]
            if color_bar_var is not None:
                ds_dicts_color = ds_dicts_color[ds_dicts_color['time'] + 1 == month]

            title = title + f' month = {month}'

    if all([lat is not None, lon is not None]):

        ds_dicts_x['dist'] = np.array([haversine(ila,ilo,lat, lon) for ila, ilo in zip(ds_dicts_x['lat'].values,ds_dicts_x['lon'].values)])
        ds_dicts_y['dist'] = np.array([haversine(ila,ilo,lat, lon) for ila, ilo in zip(ds_dicts_y['lat'].values,ds_dicts_y['lon'].values)])
        if color_bar_var is not None:
            ds_dicts_color['dist'] = np.array([haversine(ila,ilo,lat, lon) for ila, ilo in zip(ds_dicts_color['lat'].values,ds_dicts_color['lon'].values)])

        if max_dist is None:
            ds_dicts_x = ds_dicts_x[ds_dicts_x['dist'] == ds_dicts_x['dist'].min()]
            ds_dicts_y = ds_dicts_y[ds_dicts_y['dist'] == ds_dicts_y['dist'].min()]
            if color_bar_var is not None:
                ds_dicts_color = ds_dicts_color[ds_dicts_color['dist'] == ds_dicts_color['dist'].min()]
            lat = np.unique(ds_dicts_x['lat'].values)[0]
            lon = np.unique(ds_dicts_x['lon'].values)[0]

            title = title + f' lat : {lat} lon {lon}' 

        else:
            ds_dicts_x = ds_dicts_x[ds_dicts_x['dist'] <= max_dist * 1000] 
            ds_dicts_y = ds_dicts_y[ds_dicts_y['dist'] <= max_dist * 1000]    
            if color_bar_var is not None:
                ds_dicts_color = ds_dicts_color[ds_dicts_color['dist'] <= max_dist * 1000] 
            title = title + f'\n within {max_dist}km from lat : {lat} lon {lon}' 


    for ind, ds in enumerate(ds_list_var_y):

        if month_based_scatter:
            fig, ax = plt.subplots(1, 2 ,figsize=figsize)
            ax0 = ax[0]
            ax1 = ax[1]
        else:
            fig, ax0 = plt.subplots(1, 1 ,figsize=figsize)

        if len(ds_list_var_x) > 1:
            ref_ds = ds_list_var_x[ind]
        else:
            ref_ds = ds_list_var_x[0]

        if color_bar_var is not None:
            if len(ds_list_color_bar) > 1:
                colorbar_ds = ds_list_color_bar[ind]
            else:
                colorbar_ds = ds_list_color_bar[0]

        if ds == ref_ds:
            title0 = f'{ds} - {title}' ###
        else:
            title0 = f' {ds} vs {ref_ds} - {title}' 

        r2 = np.round(corr(ds_dicts_x[ref_ds], ds_dicts_y[ds]) ** 2,2)
        rmse_score = np.round(rmse(ds_dicts_x[ref_ds], ds_dicts_y[ds]),2)
        
        if color_bar_var is not None:
            c = ds_dicts_color[colorbar_ds]
            label = f'{color_bar_var} ({units[color_bar_var]}) \n {colorbar_ds}'
        else:
            c = ds_dicts_x['lev']
            label = 'depth'


        scatter = ax0.scatter(ds_dicts_x[ref_ds], ds_dicts_y[ds], c = c, cmap = cmap, vmin = vmin, vmax = vmax)

        # ax0.set_ylabel(f'{ds} ({units[var_y]})')
        # ax0.set_xlabel(f'{ref_ds} ({units[var_x]})')
        ax0.set_ylabel(f'{var_y} ({units[var_y]})')
        ax0.set_xlabel(f'{var_x} ({units[var_x]})')
        fig.colorbar(scatter, ax=ax0, label = label)
        # ax[ind,0].colorbar(scatter, label = label)

        if calculate_r2:
            title0 = title0 + f'\n $R^{2}$ score {r2}'
        if calculate_rmse:
            if calculate_r2:
                title0 = title0 + f'- RMSE {rmse_score}'
            else:
                title0 = title0 + f'\n RMSE {rmse_score}'



        ax0.set_title(title0)
        
        if plot_lim_y is not None:
            ax0.set_ylim(plot_lim_y[0], plot_lim_y[1])
        if plot_lim_x is not None:
            ax0.set_xlim(plot_lim_x[0], plot_lim_x[1])
            
        
        if add_1_1_line:
            x = np.linspace(ax0.get_xlim()[0], ax0.get_xlim()[1]  , 100)
            y = np.linspace(ax0.get_ylim()[0], ax0.get_ylim()[1]  , 100)
            ax0.plot(x,
                y,
                color='k',
                alpha = 1)
        
        if density_contour:
            SS, TT = np.meshgrid(np.linspace(ax0.get_xlim()[0], ax0.get_xlim()[1]  , 100), np.linspace(ax0.get_ylim()[0], ax0.get_ylim()[1]  , 100))
            rho = sw.dens(SS, TT, 0)
            sigma0  = rho - 1000
            contour = ax0.contour(SS, TT, sigma0, levels=15, colors='black', alpha = 0.5)
            ax0.clabel(contour, inline=True, fontsize=8)


        if add_zeros_grid:
            ax0.hlines(0, xmin = ax0.get_xlim()[0], xmax = ax0.get_xlim()[1],color = 'k', alpha = 0.5)
            ax0.vlines(0, ymin = ax0.get_ylim()[0], ymax = ax0.get_ylim()[1],color = 'k', alpha = 0.5)
        
        if month_based_scatter:
            
            scatter = ax1.scatter(ds_dicts_x[ref_ds], ds_dicts_y[ds], c = ds_dicts_x['time'] + 1, cmap = 'tab10')
            ax1.set_ylabel(f'{ds} ({units[var_y]})')
            ax1.set_xlabel(f'{ref_ds} ({units[var_x]})')
            fig.colorbar(scatter, ax=ax1, label = 'month')
            # ax[ind,1].colorbar(scatter, label = 'month')


            if plot_lim_y is not None:
                ax1.set_ylim(plot_lim_y[0], plot_lim_y[1])
            if plot_lim_x is not None:
                ax1.set_xlim(plot_lim_x[0], plot_lim_x[1])

            if density_contour:
                SS, TT = np.meshgrid(np.linspace(ax1.get_xlim()[0], ax1.get_xlim()[1]  , 100), np.linspace(ax1.get_ylim()[0], ax1.get_ylim()[1]  , 100))
                rho = sw.dens(SS, TT, 0)
                sigma0  = rho - 1000
                contour = ax1.contour(SS, TT, sigma0, levels=15, colors='black', alpha = 0.5)
                ax0.clabel(contour, inline=True, fontsize=8)

            if add_1_1_line:
                x = np.linspace(ax1.get_xlim()[0], ax1.get_xlim()[1] + 1 , 100)
                y = np.linspace(ax1.get_ylim()[0], ax1.get_ylim()[1] + 1 , 100)
                ax1.plot(x,
                    y,
                    color='k',
                    alpha = 1)
            if add_zeros_grid:
                ax1.hlines(0, xmin = ax1.get_xlim()[0], xmax = ax1.get_xlim()[1],color = 'k', alpha = 0.5)
                ax1.vlines(0, ymin = ax1.get_ylim()[0], ymax = ax1.get_ylim()[1],color = 'k', alpha = 0.5)




def scatter_comparison_singlepanel(biome_tp, dataframe, var_y:str, ds_list_var_y: list, units:dict, var_x:str = None, ds_list_var_x: list = ['obs'], 
                       location_ref_var = None, 
                       outlier_var = None, 
                       shapes_dict:dict = None,
                       color_bar_var = 'depth',
                    #    ds_list_color_bar: list = None,
                       cmap = 'coolwarm',
                       colors_dict :dict = None,
                       alphas_dict :dict = None, 
                       density_contour = False, 
                       depth = None, 
                       month = None,
                        lat = None,
                        lon = None,
                        max_dist = None, 
                       figsize = (16,6), 
                       plot_lim_x = None, 
                       plot_lim_y = None, 
                       vmin = None, 
                       vmax = None, 
                       month_based_scatter = False, 
                       calculate_rmse = True,
                       calculate_r2 = True,
                       add_1_1_line = True,
                       add_zeros_grid = False):
    
    if shapes_dict is None:
        shapes_dict = {}
        for ind, ds in enumerate(ds_list):
            shapes_dict[ds] = shapes[ind]

    if color_bar_var is None:
        assert colors_dict is not None
        color_bar_var = 'None'
        assert month_based_scatter is False

    if var_x is None:
        var_x = var_y
    
    if density_contour:
        assert var_y == 'thetao'
        assert var_x == 'so'
  
    if len(ds_list_var_x) > 1:
        assert ds_list_var_x == ds_list_var_y
    
    if all([color_bar_var != 'None',color_bar_var.lower() != 'depth']):
            ds_list_color_bar = ds_list_var_y


    if var_x != var_y:
        assert location_ref_var is not None, 'Specify which variable to be used as the location reference'
    else:
        if plot_lim_y is not None:
            plot_lim_x = plot_lim_y
        elif plot_lim_x is not None:
            plot_lim_y =  plot_lim_x 

    if location_ref_var is None:
        location_ref_var = var_y


    location_ref = dataframe[location_ref_var][biome_tp]
    if outlier_var is not None: 
        location_ref = location_ref[location_ref['var']<outlier_var]

    keys = ['year', 'time', 'lat', 'lon' ,'lev']

    if all([color_bar_var != 'None',color_bar_var.lower() != 'depth']):
        common = (
            location_ref[keys]
            .merge(dataframe[var_y][biome_tp][keys], on=keys, how="inner")
            .merge(dataframe[var_x][biome_tp][keys], on=keys, how="inner")
            .merge(dataframe[color_bar_var][biome_tp][keys], on=keys, how="inner")
            .drop_duplicates()
        )
        ds_dicts_color = dataframe[color_bar_var][biome_tp].merge(common, on=keys, how="inner") 
        clabel = f'{color_bar_var} ({units[color_bar_var]})'
    else:
        common = (
            location_ref[keys]
            .merge(dataframe[var_y][biome_tp][keys], on=keys, how="inner")
            .merge(dataframe[var_x][biome_tp][keys], on=keys, how="inner")
            .drop_duplicates()
        )
        clabel = f'depth'
    
    ds_dicts_y = dataframe[var_y][biome_tp].merge(common, on=keys, how="inner")  # talk measured in the lab, Total scale
    ds_dicts_x = dataframe[var_x][biome_tp].merge(common, on=keys, how="inner")   # DIC measured in the lab in μmol/kg-sw

    location_ref = location_ref.merge(common, on=keys, how="inner")


    if var_x == var_y:
        title = f' {biome_tp} - {ds_dicts_y["year"].min()} - {ds_dicts_y["year"].max()} - ' ###
    else:
        title = f' {biome_tp} - {location_ref["year"].min()} - {location_ref["year"].max()} - ' 

    if any([var_x != var_y, location_ref_var != var_y]):
        title = title + f'\n location reference = {location_ref_var} '
    
    if depth is not None:
        if any([type(depth) == tuple , type(depth) == list]):
            ds_dicts_y = ds_dicts_y[(ds_dicts_y['lev'] >= depth[0]) &  (ds_dicts_y['lev'] <= depth[1])]
            ds_dicts_x = ds_dicts_x[(ds_dicts_x['lev'] >= depth[0]) &  (ds_dicts_x['lev'] <= depth[1])]
            if all([color_bar_var != 'None', color_bar_var.lower() != 'depth']):
                ds_dicts_color = ds_dicts_color[(ds_dicts_color['lev'] >= depth[0]) &  (ds_dicts_color['lev'] <= depth[1])]
            title = title + f'lev = {depth[0]} - {depth[1]}'
        else:
            ds_dicts_y = ds_dicts_y[ds_dicts_y['lev'] == depth]
            ds_dicts_x = ds_dicts_x[ds_dicts_x['lev'] == depth]
            if all([color_bar_var != 'None', color_bar_var.lower() != 'depth']):
                ds_dicts_color = ds_dicts_color[ds_dicts_color['lev'] == depth]

            title = title + f' lev = {np.round(depth,2)} '

    if month is not None:
            month_based_scatter = False
            ds_dicts_y = ds_dicts_y[ds_dicts_y['time'] +1  == month ]
            ds_dicts_x = ds_dicts_x[ds_dicts_x['time']  + 1 == month  ]
            if all([color_bar_var != 'None', color_bar_var.lower() != 'depth']):
                ds_dicts_color = ds_dicts_color[ds_dicts_color['time'] + 1 == month]

            title = title + f' month = {month}'

    if all([lat is not None, lon is not None]):

        ds_dicts_x['dist'] = np.array([haversine(ila,ilo,lat, lon) for ila, ilo in zip(ds_dicts_x['lat'].values,ds_dicts_x['lon'].values)])
        ds_dicts_y['dist'] = np.array([haversine(ila,ilo,lat, lon) for ila, ilo in zip(ds_dicts_y['lat'].values,ds_dicts_y['lon'].values)])
        if all([color_bar_var != 'None', color_bar_var.lower() != 'depth']):
            ds_dicts_color['dist'] = np.array([haversine(ila,ilo,lat, lon) for ila, ilo in zip(ds_dicts_color['lat'].values,ds_dicts_color['lon'].values)])

        if max_dist is None:
            ds_dicts_x = ds_dicts_x[ds_dicts_x['dist'] == ds_dicts_x['dist'].min()]
            ds_dicts_y = ds_dicts_y[ds_dicts_y['dist'] == ds_dicts_y['dist'].min()]
            if all([color_bar_var != 'None', color_bar_var.lower() != 'depth']):
                ds_dicts_color = ds_dicts_color[ds_dicts_color['dist'] == ds_dicts_color['dist'].min()]
            lat = np.unique(ds_dicts_x['lat'].values)[0]
            lon = np.unique(ds_dicts_x['lon'].values)[0]

            title = title + f' lat : {lat} lon {lon}' 

        else:
            ds_dicts_x = ds_dicts_x[ds_dicts_x['dist'] <= max_dist * 1000] 
            ds_dicts_y = ds_dicts_y[ds_dicts_y['dist'] <= max_dist * 1000]    
            if all([color_bar_var != 'None', color_bar_var.lower() != 'depth']):
                ds_dicts_color = ds_dicts_color[ds_dicts_color['dist'] <= max_dist * 1000] 
            title = title + f'\n within {max_dist}km from lat : {lat} lon {lon}' 


    if month_based_scatter:
        fig, ax = plt.subplots(1, 2 ,figsize=figsize)
        ax0 = ax[0]
        ax1 = ax[1]
    else:
        fig, ax0 = plt.subplots(1, 1 ,figsize=figsize)

    for ind, ds in enumerate(ds_list_var_y):

        if len(ds_list_var_x) > 1:
            ref_ds = ds_list_var_x[ind]
        else:
            ref_ds = ds_list_var_x[0]

        if all([color_bar_var != 'None', color_bar_var.lower() != 'depth']):
            if len(ds_list_color_bar) > 1:
                colorbar_ds = ds_list_color_bar[ind]
            else:
                colorbar_ds = ds_list_color_bar[0]

        if ds == ref_ds:
            title0 = title ###
        else:
            title0 = f'vs {ref_ds} - {title}' 

        r2 = np.round(corr(ds_dicts_x[ref_ds], ds_dicts_y[ds]) ** 2,2)
        rmse_score = np.round(rmse(ds_dicts_x[ref_ds], ds_dicts_y[ds]),2)
        label = ds
        if calculate_r2:
            label = label + f'\n $R^{2}$ score {r2}'
        if calculate_rmse:
            if calculate_r2:
                label = label + f'- RMSE {rmse_score}'
            else:
                label = label + f'\n RMSE {rmse_score}'
        if alphas_dict is None:
            alpha = 1
        else:
            alpha = alphas_dict[ds]

        if color_bar_var != 'None':
            if color_bar_var.lower() != 'depth':
                c = ds_dicts_color[colorbar_ds]
            else:
                c = ds_dicts_x['lev']
            if colors_dict is None:
                edgecolor = None
            else:
                edgecolors= colors_dict[ds]
            scatter = ax0.scatter(ds_dicts_x[ref_ds], ds_dicts_y[ds], c = c, cmap = cmap, edgecolors= edgecolors , marker = shapes_dict[ds],alpha = alpha, label = label, vmin = vmin, vmax = vmax)

        else:
            scatter = ax0.scatter(ds_dicts_x[ref_ds], ds_dicts_y[ds], color = colors_dict[ds], marker = shapes_dict[ds], alpha = alpha, label = label)

        ax0.legend()


        
        if month_based_scatter:
            
            scatter = ax1.scatter(ds_dicts_x[ref_ds], ds_dicts_y[ds], c = ds_dicts_x['time'] + 1, cmap = 'tab10')
            ax1.set_ylabel(f'{ds} ({units[var_y]})')
            ax1.set_xlabel(f'{ref_ds} ({units[var_x]})')
            fig.colorbar(scatter, ax=ax1, label = 'month')


    ax0.set_ylabel(f'{var_y} ({units[var_y]})')
    ax0.set_xlabel(f'{var_x} ({units[var_x]})')
    if color_bar_var != 'None':
        fig.colorbar(scatter, ax=ax0, label = clabel)
    ax0.set_title(title0)
    
    if plot_lim_y is not None:
        ax0.set_ylim(plot_lim_y[0], plot_lim_y[1])
    if plot_lim_x is not None:
        ax0.set_xlim(plot_lim_x[0], plot_lim_x[1])
    
    if add_1_1_line:
        x = np.linspace(ax0.get_xlim()[0], ax0.get_xlim()[1]  , 100)
        y = np.linspace(ax0.get_ylim()[0], ax0.get_ylim()[1]  , 100)
        ax0.plot(x,
            y,
            color='k',
            alpha = 1)
    
    if all([density_contour]):
        SS, TT = np.meshgrid(np.linspace(ax0.get_xlim()[0], ax0.get_xlim()[1]  , 100), np.linspace(ax0.get_ylim()[0], ax0.get_ylim()[1]  , 100))
        rho = sw.dens(SS, TT, 0)
        sigma0  = rho - 1000
        contour = ax0.contour(SS, TT, sigma0, levels=15, colors='black', alpha = 0.5)
        ax0.clabel(contour, inline=True, fontsize=8)

    if add_zeros_grid:
        ax0.hlines(0, xmin = ax0.get_xlim()[0], xmax = ax0.get_xlim()[1],color = 'k', alpha = 0.5)
        ax0.vlines(0, ymin = ax0.get_ylim()[0], ymax = ax0.get_ylim()[1],color = 'k', alpha = 0.5)

    if month_based_scatter:
        if plot_lim_y is not None:
            ax1.set_ylim(plot_lim_y[0], plot_lim_y[1])
        if plot_lim_x is not None:
            ax1.set_xlim(plot_lim_x[0], plot_lim_x[1])

        if density_contour:
            SS, TT = np.meshgrid(np.linspace(ax1.get_xlim()[0], ax1.get_xlim()[1]  , 100), np.linspace(ax1.get_ylim()[0], ax1.get_ylim()[1]  , 100))
            rho = sw.dens(SS, TT, 0)
            sigma0  = rho - 1000
            contour = ax1.contour(SS, TT, sigma0, levels=15, colors='black', alpha = 0.5)
            ax0.clabel(contour, inline=True, fontsize=8)

        if add_1_1_line:
            x = np.linspace(ax1.get_xlim()[0], ax1.get_xlim()[1] + 1 , 100)
            y = np.linspace(ax1.get_ylim()[0], ax1.get_ylim()[1] + 1 , 100)
            ax1.plot(x,
                y,
                color='k',
                alpha = 1)
        if add_zeros_grid:
            ax1.hlines(0, xmin = ax1.get_xlim()[0], xmax = ax1.get_xlim()[1],color = 'k', alpha = 0.5)
            ax1.vlines(0, ymin = ax1.get_ylim()[0], ymax = ax1.get_ylim()[1],color = 'k', alpha = 0.5)

import datetime
def timeseries_comaprisons(biome_tp, dataframe, var:str, ds_list: list,  units:dict, colors_dict:dict,
                       location_ref_var = None, 
                       ref_ds = 'obs',
                    #    color_bar_var = None,
                    #    ds_list_color_bar: list = None,
                       depth_range:list | tuple  = None, 
                       outlier_var = None, 
                       figsize = (50,6), 
                        fontsize = 20,
                       calculate_rmse = False,
                       calculate_r2 = False,
                       add_grid = False):
    


    if location_ref_var is None:
        location_ref_var = var


    location_ref = dataframe[location_ref_var][biome_tp]
    if outlier_var is not None: 
        location_ref = location_ref[location_ref['var']<outlier_var]

    keys = ['year', 'time', 'lat', 'lon' ,'lev']
    common = (
            location_ref[keys]
            .merge(dataframe[var][biome_tp][keys], on=keys, how="inner")
            .drop_duplicates()
        )
    
    ds_dicts = dataframe[var][biome_tp].merge(common, on=keys, how="inner")  # talk measured in the lab, Total scale
    location_ref = location_ref.merge(common, on=keys, how="inner")
    ds_dicts['yearmonth'] = ds_dicts['year'] + (ds_dicts['time'] + 0.5)/12

    title = f'{var} ({units[var]}) - {biome_tp} - {location_ref["year"].min()} - {location_ref["year"].max()}'
    offsets = np.arange(- ((len(ds_list) + 1)//2) , len(ds_list) + 1 - (len(ds_list) + 1)//2)

    if location_ref_var != var:
        title = title + f'\n location reference = {location_ref_var} '

    plt.figure(figsize = figsize)
    
    if depth_range is not None:
        ds_toplot = ds_dicts[(ds_dicts['lev'] >= depth_range[0]) &  (ds_dicts['lev'] <= depth_range[1])]
        title = title + f' {depth_range[0]} - {depth_range[1]} m'

    ds_toplot_spatial_mean = ds_toplot.groupby(['yearmonth']).mean().reset_index()
    datetimes = [datetime.datetime(int(ds_toplot['year'].values[ind]), int(ds_toplot['time'].values[ind]) + 1, 15) for ind  in range(len(ds_toplot))]
    xtick_dates = [dt.strftime('%Y-%m') for dt in datetimes]

    for ind, ds in enumerate(ds_list):
 
        plt.scatter(ds_toplot['yearmonth'] + offsets[ind] * 0.05, ds_toplot[ds], color = colors_dict[var][ds]['color'] ,label = ds, alpha = 0.5)
        plt.scatter(ds_toplot_spatial_mean['yearmonth'] + offsets[ind] * 0.05, ds_toplot_spatial_mean[ds], marker = '_', s = 100, color = colors_dict[var][ds]['color'] )

    if calculate_rmse: 
        raise NotImplementedError('RMSE calculation is not implemented yet!')
    if calculate_r2: 
        raise NotImplementedError('R2 calculation is not implemented yet!')
        

    plt.scatter(ds_toplot['yearmonth'] + offsets[-1] * 0.05, ds_toplot[ref_ds], marker = 'x' , color = 'k', label = ref_ds, alpha = 0.5)
    plt.scatter(ds_toplot_spatial_mean['yearmonth'] + offsets[-1] * 0.05, ds_toplot_spatial_mean[ref_ds], marker = '_' , color = 'k')

    plt.xticks(np.unique(ds_toplot['yearmonth'].values), np.unique(xtick_dates), rotation=90, fontsize = fontsize)
    plt.yticks(fontsize = fontsize)
    plt.ylabel(f'{units[var]}', fontsize = fontsize)
    plt.xlabel('time', fontsize = fontsize)
    plt.title(title, fontsize = fontsize)
    plt.legend( fontsize = fontsize)
    if add_grid:
        plt.grid()






def spatial_maps(biome_tp, mask_ocean, mask_biomes, dataframe, var:str, ds_list: list,  units:dict, boundaries_dict:dict, 
                       location_ref_var = None, 
                       ref_ds = 'obs',
                       depth:list | tuple | int = None, 
                       outlier_var = None, 
                       figsize = (6,6), 
                        fontsize = 20,
                        lat_min_to_show = None,
                        lat_max_to_show = None,
                        lon_min_to_show = None,
                        lon_max_to_show = None,
                        vmin = 0,
                        vmax = 100,
                        s = 5,
                        cmap = 'Reds',
                       plot_rmse = False,
                       plot_diff = False, 
                       plot_r2 = False,
                       plot_count = False,
                       calculate_rmse = False,
                       calculate_r2 = False,
                       plot_full_model_grid = False,
                       ds_dict_models :dict = None,
                       biome_boundaries_dict = None):
    
    assert sum([plot_rmse, plot_r2,  plot_count])<=1, 'specify one of plot_rmse, plot_r2 or plot_count'

    if any([plot_rmse,plot_r2, calculate_rmse, calculate_r2, plot_diff ]):
        assert ref_ds is not None

    if any([plot_rmse,plot_r2,  plot_count]):
        assert type(depth) != int, 'specify a depth range not a depth point!'

    if location_ref_var is None:
        location_ref_var = var

    location_ref = dataframe[location_ref_var][biome_tp]
    if outlier_var is not None: 
        location_ref = location_ref[location_ref['var']<outlier_var]

    keys = ['year', 'time', 'lat', 'lon' ,'lev']
    common = (
            location_ref[keys]
            .merge(dataframe[var][biome_tp][keys], on=keys, how="inner")
            .drop_duplicates()
        )
    
    ds_dicts = dataframe[var][biome_tp].merge(common, on=keys, how="inner")  # talk measured in the lab, Total scale
    location_ref = location_ref.merge(common, on=keys, how="inner")

    if depth is not None:
        if any([type(depth) == tuple, type(depth) == list]):
            ds_dicts = ds_dicts[(ds_dicts['lev'] >= depth[0]) &  (ds_dicts['lev'] <= depth[1])]
        else:
            ds_dicts = ds_dicts[ds_dicts['lev'] == depth]

    if plot_full_model_grid:
        assert ds_dict_models is not None, 'specify full model data dictionary'
        assert biome_boundaries_dict is not None, 'specify biome boundaries for model data'

    for i, ds in enumerate(ds_list):

        plt.figure(figsize = figsize)
        ax = plt.subplot(1,1,1)
        ax.pcolormesh(mask_ocean.lon[lon_min_to_show:lon_max_to_show], mask_ocean.lat[lat_min_to_show:lat_max_to_show], mask_ocean.where(mask_ocean==0).values[lat_min_to_show:lat_max_to_show,lon_min_to_show:lon_max_to_show], vmin = 0, vmax =1)
        ax.pcolormesh(mask_biomes[biome_tp].lon[lon_min_to_show:lon_max_to_show], mask_biomes[biome_tp].lat[lat_min_to_show:lat_max_to_show], mask_biomes[biome_tp].where(mask_biomes[biome_tp] == 1).values[lat_min_to_show:lat_max_to_show,lon_min_to_show:lon_max_to_show], vmin = 0, vmax =1, alpha = 0.25)


        if any([plot_rmse,plot_r2, calculate_rmse, calculate_r2, plot_diff ]):
            title = f'{var} ({units[var]}) {ds} vs {ref_ds} - {location_ref["year"].min()} - {location_ref["year"].max()}' 
        else:
            title = f'{var} ({units[var]}) {ds} - {location_ref["year"].min()} - {location_ref["year"].max()}' 

        if location_ref_var != var:
            title = title + f'\n location reference = {location_ref_var} '

        if depth is not None:
            if any([type(depth) == tuple, type(depth) == list]):
                title = title + f'- {depth[0]} - {depth[1]} m'
            else:
                title = title + f'- depth = {np.round(depth,2)}'

        if all([not plot_rmse,not plot_r2,  not plot_count, not plot_diff]):

            ds_topolot_obs = ds_dicts.groupby([ 'lat', 'lon']).agg(
                        mean=('obs' , 'mean') ,
                        count = ('obs' , 'count')).reset_index()

            if all([plot_full_model_grid, 'obs' not in ds]):
                ds_model = ds_dict_models[var][ds]
                ds_model = ds_model.where((ds_model['lat'] >=  boundaries_dict[biome_tp]['lat_min'] ) & (ds_model['lat'] <=  boundaries_dict[biome_tp]['lat_max'] ) )
                ds_model = ds_model.where((ds_model['lon'] >=  boundaries_dict[biome_tp]['lon_min'] ) & (ds_model['lon'] <=  boundaries_dict[biome_tp]['lon_max'] )  )
                if any([type(depth) == tuple, type(depth) == list]):
                    ds_model = ds_model.where((ds_model.lev >=  depth[0]) & (ds_model.lev <=  depth[1]), drop = True).mean('lev')
                else:
                    ds_model = ds_model.sel(lev = depth)

                ds_topolot_model = average_on_ref_times(ds_model, ds_dicts).mean('t')
            else:
                ds_topolot_model =  ds_dicts.groupby([ 'lat', 'lon']).agg(
                        mean=(ds , 'mean') ,
                        count = (ds , 'count')).reset_index()
            
            if 'obs' in ds:
                scatter = ax.scatter(ds_topolot_obs['lon'], ds_topolot_obs['lat'], c =  ds_topolot_obs['mean'], s = s, cmap = cmap, vmin = vmin, vmax = vmax)
            else:
                if plot_full_model_grid:
                    scatter = ax.pcolormesh(ds_topolot_model.lon[lon_min_to_show:lon_max_to_show], ds_topolot_model.lat[lat_min_to_show:lat_max_to_show], ds_topolot_model.values[lat_min_to_show:lat_max_to_show,lon_min_to_show:lon_max_to_show], cmap = cmap, vmin = vmin, vmax = vmax)
                else:
                    scatter = ax.scatter(ds_topolot_model['lon'], ds_topolot_model['lat'], c =  ds_topolot_model['mean'], s = s, cmap = cmap, vmin = vmin, vmax = vmax)


            plt.colorbar(scatter, label = f'{var} ({units[var]})')

        elif plot_rmse:
            scores_dicts =  score_pattern(ds_dicts, ds, ref_ds)
            scatter = ax.scatter(scores_dicts['lon'], scores_dicts['lat'], c =  np.sqrt(scores_dicts['mse']), s = s, cmap = 'Reds', vmin = vmin, vmax = vmax)
            plt.colorbar(scatter, label = f'RMSE ({units[var]})')
        elif plot_diff:
            scores_dicts =  ds_dicts[ds]  - ds_dicts[ref_ds] 
            scatter = ax.scatter(scores_dicts['lon'], scores_dicts['lat'], c =  scores_dicts, s = s, cmap = 'Reds', vmin = vmin, vmax = vmax)
            plt.colorbar(scatter, label = f'Diff ({units[var]})')
        elif plot_r2:
            scores_dicts =  score_pattern(ds_dicts, ds, ref_ds, score_function='R2')
            scatter = ax.scatter(scores_dicts['lon'], scores_dicts['lat'], c =  (scores_dicts['corr'])**2, s = s, cmap = 'Reds', vmin = 0, vmax = 1)
            plt.colorbar(scatter, label = f'R2')
        elif plot_count:
            ds_topolot_count =  ds_dicts.groupby([ 'lat', 'lon']).agg(
                        count = (ds , 'count')).reset_index()
            scatter = ax.scatter(ds_topolot_count['lon'], ds_topolot_count['lat'], c =  (ds_topolot_count['count']), s = s, vmin = 0, vmax = vmax, cmap = 'Blues_r')
            plt.colorbar(scatter, label = 'count')

        if any([calculate_rmse, calculate_r2]):
                r2 = np.round(corr(ds_dicts[ref_ds], ds_dicts[ds]) ** 2,2)
                rmse_score = np.round(rmse(ds_dicts[ref_ds], ds_dicts[ds]),2)
        if calculate_r2:
            title = title + f'\n $R^{2}$ score {r2}'
        if calculate_rmse:
            if calculate_r2:
                title = title + f'- RMSE {rmse_score}'
            else:
                title = title + f'\n RMSE {rmse_score}'

        plt.ylabel('lat')
        plt.xlabel('lon')
        plt.title(title)




def spatial_maps_climatology(biome_tp, mask_ocean, mask_biomes, ds_dict, var:str, ds_list: list,  units:dict, boundaries_dict:dict, 
                       ref_ds = 'obs',
                       depth:list | tuple | int = None, 
                       figsize = (6,6), 
                        fontsize = 20,
                        lat_min_to_show = None,
                        lat_max_to_show = None,
                        lon_min_to_show = None,
                        lon_max_to_show = None,
                        vmin = 0,
                        vmax = 100,
                        cmap = 'Reds',
                       plot_rmse = False,
                       plot_diff = False, 
                       calculate_rmse = False,
                       biome_boundaries_dict = None):
    

    if any([plot_rmse, calculate_rmse, plot_diff ]):
        assert ref_ds is not None

    if any([plot_rmse, plot_diff]):
        if depth is not None:   
            assert type(depth) != int, 'specify a depth range not a depth point!'

    if ref_ds is not None:

            ds_ref = ds_dict[var][ref_ds]
            ds_ref = ds_ref.where((ds_ref['lat'] >=  boundaries_dict[biome_tp]['lat_min'] ) & (ds_ref['lat'] <=  boundaries_dict[biome_tp]['lat_max'] ) )
            ds_ref = ds_ref.where((ds_ref['lon'] >=  boundaries_dict[biome_tp]['lon_min'] ) & (ds_ref['lon'] <=  boundaries_dict[biome_tp]['lon_max'] )  )
            if depth is not None:   
                if any([type(depth) == tuple, type(depth) == list]):
                    ds_ref = ds_ref.where((ds_ref.lev >=  depth[0]) & (ds_ref.lev <=  depth[1]), drop = True)


    for i, ds in enumerate(ds_list):

        plt.figure(figsize = figsize)
        ax = plt.subplot(1,1,1)
        ax.pcolormesh(mask_ocean.lon[lon_min_to_show:lon_max_to_show], mask_ocean.lat[lat_min_to_show:lat_max_to_show], mask_ocean.where(mask_ocean==0).values[lat_min_to_show:lat_max_to_show,lon_min_to_show:lon_max_to_show], vmin = 0, vmax =1)
        ax.pcolormesh(mask_biomes[biome_tp].lon[lon_min_to_show:lon_max_to_show], mask_biomes[biome_tp].lat[lat_min_to_show:lat_max_to_show], mask_biomes[biome_tp].where(mask_biomes[biome_tp] == 1).values[lat_min_to_show:lat_max_to_show,lon_min_to_show:lon_max_to_show], vmin = 0, vmax =1, alpha = 0.25)


        if any([plot_rmse,calculate_rmse, plot_diff ]):
            title = f'{var} climatology ({units[var]}) {ds} vs {ref_ds}' 
        else:
            title = f'{var} climatology ({units[var]}) {ds}' 


        if depth is not None:
            if any([type(depth) == tuple, type(depth) == list]):
                title = title + f'- {depth[0]} - {depth[1]} m'
            else:
                title = title + f'- depth = {np.round(depth,2)}'


        ds_model = ds_dict[var][ds]
        ds_model = ds_model.where((ds_model['lat'] >=  boundaries_dict[biome_tp]['lat_min'] ) & (ds_model['lat'] <=  boundaries_dict[biome_tp]['lat_max'] ) )
        ds_model = ds_model.where((ds_model['lon'] >=  boundaries_dict[biome_tp]['lon_min'] ) & (ds_model['lon'] <=  boundaries_dict[biome_tp]['lon_max'] )  )
        if depth is not None:   
            if any([type(depth) == tuple, type(depth) == list]):
                    ds_model = ds_model.where((ds_model.lev >=  depth[0]) & (ds_model.lev <=  depth[1]), drop = True)

        if any([calculate_rmse]):
                rmse_score = np.round(np.sqrt((( ds_ref - ds_model)**2).mean()).values,2)

        if plot_rmse:
            ds_model_to_plot = np.sqrt((( ds_ref - ds_model)**2).mean('lev'))
            # scatter = ax.pcolormesh(ds_model.lon[lon_min_to_show:lon_max_to_show], ds_model.lat[lat_min_to_show:lat_max_to_show], ds_model.values[lat_min_to_show:lat_max_to_show,lon_min_to_show:lon_max_to_show], cmap = cmap, vmin = vmin, vmax = vmax)
            # plt.colorbar(scatter, label = f'RMSE ({units[var]})')
        elif plot_diff:
            ds_model_to_plot =  ds_ref.mean('lev') - ds_model.mean('lev')
            # scatter = ax.scatter(scores_dicts['lon'], scores_dicts['lat'], c =  scores_dicts, s = s, cmap = 'Reds', vmin = vmin, vmax = vmax)
            # plt.colorbar(scatter, label = f'Diff ({units[var]})')
        else:         
            if depth is not None:   
                if any([type(depth) == tuple, type(depth) == list]):
                    ds_model_to_plot = ds_model.mean('lev')
                else:
                    ds_model_to_plot = ds_model.sel(lev = depth)
            else:
                ds_model_to_plot = ds_model.mean('lev')
            

        scatter = ax.pcolormesh(ds_model_to_plot.lon[lon_min_to_show:lon_max_to_show], ds_model_to_plot.lat[lat_min_to_show:lat_max_to_show], ds_model_to_plot.values[lat_min_to_show:lat_max_to_show,lon_min_to_show:lon_max_to_show], cmap = cmap, vmin = vmin, vmax = vmax)
        plt.colorbar(scatter, label = f'{var} ({units[var]})')


        if calculate_rmse:

            title = title + f'\n RMSE {rmse_score}'

        plt.ylabel('lat')
        plt.xlabel('lon')
        plt.title(title)






def score_pattern(data_frame: pd.DataFrame, model_run, ref, score_function = 'MSE'):
    ds = data_frame.copy()
    if score_function.lower() == 'mse':
        ds[f'{model_run} - {ref}'] = (data_frame[model_run] - data_frame[ref])**2
        return ds.groupby([ 'lat', 'lon']).agg(
        mse=(f'{model_run} - {ref}', 'mean') ,
        count = (f'{model_run} - {ref}', 'count')).reset_index()
    elif score_function.lower() == 'r2':
        return ds.groupby([ 'lat', 'lon']).apply(lambda g: g[model_run].corr(g[ref])).reset_index(name='corr') 
    

def average_on_ref_times(ds: xr.DataArray, ref_dataframe : pd.DataFrame):
    ref = ref_dataframe.groupby(['year','time']).agg(
                    mean=('obs' , 'mean') ,
                    count = ('obs' , 'count')).reset_index()
    return xr.concat([ds.sel(year = ref['year'][i], time = ref['time'][i]) for i in range(len(ref))], dim = 't')




def point_profile(biome_tp,  dataframe, var:str, ds_list: list,  units:dict,
                       location_ref_var = None, 
                       ref_ds = 'obs',
                       depth:list | tuple = None, 
                       outlier_var = None, 
                       figsize = (6,6), 
                        fontsize = 20,
                        lat = None,
                        lon = None,
                        xlims = None,
                        ylims = None,
                        max_dist = None,
                        isolate_time = False,
                        mean_lev = False,
                        legend = True,
                        bbox_to_anchor=(1, 0.5),
                        shapes_dict = None,
                        colors_dict = None ,
                        monthly_colorbar = False,
                       calculate_rmse = False,
                       calculate_r2 = False):


    if any([calculate_rmse, calculate_r2 ]):
        assert ref_ds is not None


    if location_ref_var is None:
        location_ref_var = var

    location_ref = dataframe[location_ref_var][biome_tp]
    if outlier_var is not None: 
        location_ref = location_ref[location_ref['var']<outlier_var]

    keys = ['year', 'time', 'lat', 'lon' ,'lev']
    common = (
            location_ref[keys]
            .merge(dataframe[var][biome_tp][keys], on=keys, how="inner")
            .drop_duplicates()
        )
    
    ds_dicts = dataframe[var][biome_tp].merge(common, on=keys, how="inner")  # talk measured in the lab, Total scale
    location_ref = location_ref.merge(common, on=keys, how="inner")

    if depth is not None:
            ds_dicts = ds_dicts[(ds_dicts['lev'] >= depth[0]) &  (ds_dicts['lev'] <= depth[1])]


    if all([lat is None, lon is None]):
        max_dist = None
        ds_topolot_obs = ds_dicts.groupby([ 'lat', 'lon']).agg(
                    mean=('obs' , 'mean') ,
                    count = ('obs' , 'count')).reset_index()
        lat = ds_topolot_obs.loc[ds_topolot_obs['count'].idxmax()]['lat']
        lon = ds_topolot_obs.loc[ds_topolot_obs['count'].idxmax()]['lon']
        ds_dicts_toplot = ds_dicts[ds_dicts['lat'] == lat]
        ds_dicts_toplot = ds_dicts_toplot[ds_dicts_toplot['lon'] == lon]
  
    else:
        ds_dicts_toplot = ds_dicts.copy()
        ds_dicts_toplot['dist'] = np.array([haversine(ila,ilo,lat, lon) for ila, ilo in zip(ds_dicts_toplot['lat'].values,ds_dicts_toplot['lon'].values)])
        if max_dist is None:
            ds_dicts_toplot = ds_dicts_toplot[ds_dicts_toplot['dist'] == ds_dicts_toplot['dist'].min()]
            lat = np.unique(ds_dicts_toplot['lat'].values)[0]
            lon = np.unique(ds_dicts_toplot['lon'].values)[0]
        else:
            ds_dicts_toplot = ds_dicts_toplot[ds_dicts_toplot['dist'] <= max_dist * 1000]       

    if isolate_time:
  
            counts = (
                ds_dicts_toplot
                .groupby(["year", "time"])
                .size()
                .reset_index(name="count")
            )
            row = counts.loc[counts["count"].idxmax()]
            year_sel = row["year"]
            time_sel = row["time"]
            ds_dicts_toplot = ds_dicts_toplot[(ds_dicts_toplot['year'] == year_sel )& (ds_dicts_toplot['time'] == time_sel )]

    if mean_lev:
            ds_dicts_toplot = (ds_dicts_toplot.groupby(['lev'], as_index=False).mean(numeric_only=True))

                           
    if ylims is None:
        ylim = (0, ds_dicts_toplot['lev'].max())

    if shapes_dict is None:
        assert len(ds_list) <= len(markers), 'Not enough marker shapes for this many number of data to plot'
        shapes_dict = {}
        shapes = random.sample(markers, len(ds_list))
        for ind, ds in enumerate(ds_list):
            shapes_dict[ds] = shapes[ind]

    fig = plt.figure(figsize = figsize)
    ax = plt.subplot(1,1,1)
    if isolate_time:
        time_title = f'year {year_sel} month {time_sel + 1} '
    else:
        time_title = f'{location_ref["year"].min()} - {location_ref["year"].max()}' 


    if all([max_dist is not None, lat is not None, lon is not None]):
        title = f'{var} - {time_title} \n within {max_dist}km from lat : {lat} lon {lon}' 
    else:
        title = f'{var} - {time_title} - lat : {lat} lon {lon}' 

    if location_ref_var != var:
            title = title + f'\n location reference = {location_ref_var} '
    if mean_lev:
            title = title + '(mean per lev)'
    for i, ds in enumerate(ds_list):


        label = ds
        if any([calculate_rmse, calculate_r2]):
                r2 = np.round(corr(ds_dicts_toplot[ref_ds], ds_dicts_toplot[ds]) ** 2,2)
                rmse_score = np.round(rmse(ds_dicts_toplot[ref_ds], ds_dicts_toplot[ds]),2)
        if all([calculate_r2, ds != 'obs']):
            label = label + f'- $R^{2}$ score (vs {ref_ds}) {r2}'
        if all([calculate_rmse, ds != 'obs']):
            if calculate_r2:
                label = label + f'- RMSE {rmse_score}'
            else:
                label = label + f'- RMSE  (vs {ref_ds}){rmse_score}'
        if monthly_colorbar:
            colorbar = ds_dicts_toplot['time'] + 1
        else:
            if colors_dict is None:
                colorbar = random.choice(list(colors.CSS4_COLORS.values()))
            else:
                colorbar = colors_dict[ds]
        if ds == 'obs':
             scatter = plt.scatter( ds_dicts_toplot[ds], ds_dicts_toplot['lev'], c = colorbar, marker = 'X' , cmap = 'tab10', label = label )
        else:
            scatter = plt.scatter( ds_dicts_toplot[ds], ds_dicts_toplot['lev'],c = colorbar, marker = shapes_dict[ds] , cmap = 'tab10', label = label )  
    if monthly_colorbar:
        fig.colorbar(scatter, ax=ax, label = 'month')
    if legend:
        plt.legend(loc='center left', bbox_to_anchor=bbox_to_anchor)
    plt.ylabel('depth')
    plt.xlabel(var + f' ({units[var]})')
    plt.xlim(xlims)
    plt.ylim(ylims)
    plt.title(title)
    plt.gca().invert_yaxis()




def point_profile_climatology( ds_dict, var:str, ds_list: list,  units:dict,lat, lon,
                       ref_ds = 'obs',
                       depth:list | tuple = None, 
                       figsize = (6,6), 
                        fontsize = 20,
                        xlims = None,
                        ylims = None,
                        # max_dist = None,
                        legend = False,
                        bbox_to_anchor=(1, 0.5),
                        shapes_dict = None,
                        colors_dict = None ,
                       calculate_rmse = False):


    if any([calculate_rmse ]):
        assert ref_ds is not None



    if ref_ds is not None:

        ds_ref = ds_dict[var][ref_ds]
        if depth is not None:   
            if any([type(depth) == tuple, type(depth) == list]):
                ds_ref = ds_ref.where((ds_ref.lev >=  depth[0]) & (ds_ref.lev <=  depth[1]), drop = True)

        ds_ref = ds_ref.sel(lat = lat, lon = lon, method = 'nearest')
                

    if shapes_dict is None:
        assert len(ds_list) <= len(markers), 'Not enough marker shapes for this many number of data to plot'
        shapes_dict = {}
        shapes = random.sample(markers, len(ds_list))
        for ind, ds in enumerate(ds_list):
            shapes_dict[ds] = shapes[ind]

    fig = plt.figure(figsize = figsize)
    ax = plt.subplot(1,1,1)



    title = f'{var} - climatology - lat : {lat} lon {lon}' 

    for i, ds in enumerate(ds_list):

        ds_model = ds_dict[var][ds]
        if depth is not None:   
            if any([type(depth) == tuple, type(depth) == list]):
                    ds_model = ds_model.where((ds_model.lev >=  depth[0]) & (ds_model.lev <=  depth[1]), drop = True)
        
        ds_model = ds_model.sel(lat = lat, lon = lon, method = 'nearest')

        if ylims is None:
            ylim = (0, ds_model['lev'].max())

        label = ds
        if any([calculate_rmse]):
                rmse_score = np.round(np.sqrt(((ds_model - ds_ref)**2).mean()).values,2)

        if all([calculate_rmse, ds != 'obs']):

                label = label + f'- RMSE  (vs {ref_ds}){rmse_score}'

        if colors_dict is None:
            colorbar = random.choice(list(colors.CSS4_COLORS.values()))
        else:
            colorbar = colors_dict[ds]
        if ds == 'obs':
             scatter = plt.scatter( ds_model, ds_model['lev'], c = colorbar, marker = 'X' , cmap = 'tab10', label = label )
        else:
            scatter = plt.scatter( ds_model, ds_model['lev'],c = colorbar, marker = shapes_dict[ds] , cmap = 'tab10', label = label )  
    if legend:
        plt.legend(loc='center left', bbox_to_anchor=bbox_to_anchor)
    plt.ylabel('depth')
    plt.xlabel(var + f' ({units[var]})')
    plt.xlim(xlims)
    plt.ylim(ylims)
    plt.title(title)
    plt.gca().invert_yaxis()





markers = [
    ".",   # point
    ",",   # pixel
    "o",   # circle
    "v",   # triangle down
    "^",   # triangle up
    "<",   # triangle left
    ">",   # triangle right
    "1",   # tri-down
    "2",   # tri-up
    "3",   # tri-left
    "4",   # tri-right
    "8",   # octagon
    "s",   # square
    "p",   # pentagon
    "P",   # plus (filled)
    "*",   # star
    "h",   # hexagon 1
    "H",   # hexagon 2
    "+",   # plus
    "x",   # x
    "D",   # diamond
    "d",   # thin diamond
    "|",   # vertical line
    "_",   # horizontal line
]