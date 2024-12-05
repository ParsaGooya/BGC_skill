import numpy as np
import xarray as xr
from pathlib import Path
import re


def coords_edit(ds):
    
    lat = ds.lat.values[:,0]
    lon = ds.lon.values[0,:]
    ds_like = xr. DataArray(ds.values, dims = ds.dims).rename({'x' :'lon', 'y': 'lat'}).assign_coords({'lat': lat, 'lon': lon })

    if 'member' in ds.dims:
        ds_like = ds_like.assign_coords({'member' : ds.member})
    if 'time' in ds.dims:
        ds_like = ds_like.assign_coords({'time' : ds.time})
    return ds_like


def load_biomes(dir_in,
                key='*',
                verbose=True):
    if verbose:
        print("loading global ocean biomes..")
    for f in sorted(Path(dir_in).glob(f"{key}.nc")):
        ds = xr.open_dataset(f)
        ds = ds.transpose('lat',
                          'lon',
                          'year')
    if verbose:
        print("done")
    return ds


# def load_data(dir_in,
#               verbose=True):
#     if verbose:
#         print("loading data..")
#     datasets = []
#     for f in sorted(Path(dir_in).glob(f"*.nc")):
#         ds = xr.open_dataset(f)
#         fname_split = f.stem.split("_")   # ripf info
#         iyr = int(fname_split[-2])
#         ds = ds.assign_coords(year=iyr)
#         # ds['time'] = np.arange(len(ds.time.dt.month))
#         ds['time'] = np.arange(ds.time.size)
#         ds = ds.expand_dims('year',axis=0)
#         datasets.append(ds)
#     ds_combined = xr.concat(datasets, dim='year').sortby('year')
#     if verbose:
#         print("done")
#     return ds_combined

def load_data(dir_in,var,
              verbose=True,
              ensemble_mean = False,
              ensemble_id = None,
              load = True):
    if verbose:
        print("loading data..")

    ds = (xr.open_mfdataset(dir_in, combine = 'nested', concat_dim = 'time')[var]).transpose('time', ...) 
    if ensemble_id is not None:
        ds = ds.sel(ensembles = ensemble_id)
    if ensemble_mean:
        ds = ds.mean('ensembles')

    try:
        ds = xr.concat([ds[i:i+12].expand_dims('year',axis = 0).assign_coords(year =  [int(np.datetime_as_string(ds[i].time.values)[:4])]).assign_coords(time = np.arange(len(ds[i:i+12].time))) for i in range(0,len(ds.time),12)], dim = 'year')
    except:
        ds = xr.concat([ds[i:i+12].expand_dims('year',axis = 0).assign_coords(year =  [ds[i].time.item().year]).assign_coords(time = np.arange(len(ds[i:i+12].time))) for i in range(0,len(ds.time),12)], dim = 'year')
    if verbose:
        print("done")
    if load:
        return ds.load()
    else:
        return ds

def load_ensemble(dir_in,
# def load_forecasts(dir_in,
                   var,
                   y0=1981,
                   y1=2022,
                   verbose=True,
                   ensemble_mean = False):
    """
     read in hindcast/forecast data
    """
    try:
        data = []
        if verbose:
            print('loading forecasts..')
        datasets = []
        for f in sorted(Path(dir_in).glob(f"*.nc")):
            ds = xr.open_dataset(f)[var].mean('ensembles').load() if ensemble_mean else xr.open_dataset(f)[var].load()
            fname_split = f.stem.split("_")   # ripf info
            iyr = int(fname_split[-6])
            ds = ds.assign_coords(year=iyr+1)
            ds = ds.rename({'lead_time' : 'time'})
            ds['time'] = np.arange(ds.time.size)
            ds = ds.expand_dims('year',axis=0)
            datasets.append(ds)
        ds_combined = xr.concat(datasets, dim='year').sortby('year')
        if verbose:
            print("done")
    except IOError:
        print('The input forecast file is missing! --check proper location')
        quit()
    return ds_combined

def load_forecasts(dir_in,
# def load_nadj_forecasts(dir_in,
                        key='*',
                        var_in='nn_adjusted',
                        var_out='fgco2',
                        verbose=True):
    if verbose:
        print("loading data..")
    ds = xr.open_mfdataset(f'{dir_in}/{key}*.nc',
                           combine='nested',
                           concat_dim='year').squeeze().load()[var_in].to_dataset(name=var_out)   
    ds = ds.rename({'lead_time':'time'})
    ds['time'] = np.arange(ds.time.size)
    if verbose:
        print("done")
    return ds


def load_adj_forecasts(dir_in,
                       key='*',
                       var='fgco2',
                       verbose=True):
    if verbose:
        print("loading data..")
    for f in sorted(Path(dir_in).glob(f"{key}.nc")):
        ds = xr.open_dataset(f)   
        ds = ds.rename({'lead_time':'time'})
        ds['time'] = np.arange(ds.time.size)
    if verbose:
        print("done")
    return ds

def load_nadj_ensemble(dir_in,
                       var_in='nn_adjusted',
                       var_out='fgco2',
                       etype='LE_members',
                       EE=10,
                       verbose=True):
    
    if verbose:
        print("loading ensemble..")

    if etype == 'LE_members':
        ds_combined = xr.open_mfdataset(f'{dir_in}/*_LE.nc',
                                        combine='nested',
                                        concat_dim='year').squeeze().load()[var_in].to_dataset(name=var_out)   
        ds_combined = ds_combined.rename({'lead_time':'time'})
        ds_combined['time'] = np.arange(ds_combined.time.size)
        
    if etype == 'LE_seeds':
        datasets = []
        for iE in np.arange(EE):
            iE1 = iE + 1
            ds = load_forecasts(f'{dir_in}/E{iE1}',
            # ds = load_nadj_forecasts(f'{dir_in}/E{iE1}',
                                    verbose=False)
            ds = ds.assign_coords(ensembles=iE1)
            ds = ds.expand_dims('ensembles',axis=1)
            datasets.append(ds)
            ds_combined = xr.concat(datasets,
                                    dim='ensembles').sortby('ensembles')

    if verbose:
        print("done")
        
    return ds_combined


