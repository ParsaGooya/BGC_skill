import warnings
warnings.filterwarnings('ignore') # don't output warnings

import os
import xarray as xr
import cftime
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import numpy as np
import xesmf as xe
from pathlib import Path
import glob
import wget

########################################################## Regrid observations ##########################################################

source = 'Shujie23' # or 'ESACCI'  ## Choose obsearvational source dataset
var_name = 'chlor_a'
out_dir = f'/space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/chlos/observations/{source}/raw/'  #### path to local directory where raw obseravations are downloaded

def regridder(ds, var):  #### averagin into 1x1 resolution
    ds = ds[var].transpose(...,'lat','lon') if var is not None else ds.transpose(...,'lat','lon')
    ds = xr.concat([ds[...,i:i+24,:].mean('lat') for i in range(0, len(ds.lat), 24)], dim = 'lat').assign_coords(lat = np.arange(-89.5,90,1)*-1)
    ds = xr.concat([ds[...,:,i:i+24].mean('lon') for i in range(0, len(ds.lon), 24)], dim = 'lon').assign_coords(lon = np.arange(-179.5,180,1))
    return ds

# ds_out = xe.util.grid_global(1, 1)
for link in glob.glob(out_dir + 'raw/*.nc'):
    print(f'start regridding {link.split("/")[-1]} ...\n')
    ds = xr.open_dataset(link)
    if source == 'Shujie23':
        # ds = ds[var_name].expand_dims('time', axis = 0)
        year, month = np.divmod(int(link.split("/")[-1].split('.')[0][1:]),100)
        ds = ds.assign_coords(time = np.datetime64(f'{year:04d}-{month:02d}-01T00:00:00.000000000', 'ns'))
    # regridder = xe.Regridder(ds, ds_out, 'bilinear', ignore_degenerate=True, periodic=True)
    ds= regridder(ds, var_name)
    ds.to_dataset(name = 'chlos').to_netcdf(out_dir + 'averaged/' +  link.split("/")[-1].split('4km')[0] + '1x1_averaged'+ link.split("/")[-1].split('4km')[1])
    print(f'{link.split("/")[-1]} saved!\n') ### save regridded obsearvations

############### concat individual years and save as one .nc file ######################
ds = xr.concat([xr.open_dataset(link) for link in glob.glob(out_dir + 'averaged/*.nc')], dim = 'time').sortby('time').transpose(...,'time','lat','lon')
min_time = np.datetime_as_string(ds.time.min().values)[:4] + np.datetime_as_string(ds.time.min().values)[5:7]
max_time = np.datetime_as_string(ds.time.max().values)[:4] + np.datetime_as_string(ds.time.max().values)[5:7]
ds.to_netcdf(out_dir.split('/raw/')[0] + f'/{source}-1M_MONTHLY_1x1_{min_time}-{max_time}.nc')