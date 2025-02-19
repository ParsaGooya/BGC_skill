
from pathlib import Path
import glob
import numpy as np
import xarray as xr
import xesmf as xe

def coords_edit(ds):
    lat = ds.lat.values[:,0]
    lon = ds.lon.values[0,:]
    ds_like = xr. DataArray(ds.values, dims = ds.dims).rename({'x' :'lon', 'y': 'lat'}).assign_coords({'lat': lat, 'lon': lon })
    for d in ds.dims:
        if d not in ['x','y']:
            ds_like = ds_like.assign_coords(d = ds[d])
    return ds_like

def lonfixer(ds):
    len_lon = len(ds.lon)
    ds0 = ds.isel(lon = np.arange(0, (len_lon//2)))
    ds1 = ds.isel(lon = np.arange((len_lon//2), len_lon)).assign_coords(lon = ds.lon.values[len_lon//2:] - 360)
    return xr.concat([ds1,ds0], dim = 'lon')

data_dir = '/space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/soda/raw'
rename_dict = {'temp' : 'thetao', 'salt' : 'so', 'wt' : 'wo', 'u':'uo','v':'vo'}
# rename_dict = { 'wt' : 'wo', 'u':'uo','v':'vo'}

ds_out = xe.util.grid_global(1, 1) 


for var in rename_dict.keys():

    ds = xr.open_mfdataset(str(Path(data_dir, "*.nc")), combine='nested', concat_dim='time')[var]
    out_dir = f'/space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/{rename_dict[var]}/observation/SODA/'

    Path(out_dir).mkdir(parents=True, exist_ok=True)
    print(f'start processing {var} ...')
    try:
        ds = ds.rename({'st_ocean' : 'lev'})
    except:
        ds = ds.rename({'sw_ocean' : 'lev'})
    try:
        ds = ds.rename({'yt_ocean' : 'lat','xt_ocean' : 'lon'})
    except:
        ds = ds.rename({'yu_ocean' : 'lat','xu_ocean' : 'lon'})

    ds_merid = []
    for lon in range(0,360):
        ds_merid.append(ds.where((ds.lon>lon) & (ds.lon<lon + 1), drop = True).mean('lon').expand_dims( 'lon', axis = -1).assign_coords(lon = [lon + 0.5]))
    ds_merid = xr.concat(ds_merid, dim = 'lon')
    ds_zonal = []
    for lat in range(-74,90):
        ds_zonal.append(ds_merid.where((ds_merid.lat>lat) & (ds_merid.lat<lat + 1), drop = True).mean('lat').expand_dims( 'lat', axis = -2).assign_coords(lat = [lat + 0.5]))   
    ds = xr.concat(ds_zonal, dim = 'lat')
       
    # regridder = xe.Regridder(ds, ds_out, 'bilinear', ignore_degenerate=True, periodic=True)
    # ds = coords_edit(regridder(ds.load()))

    lonfixer(ds).to_dataset(name = rename_dict[var]).to_netcdf(out_dir + f'{rename_dict[var]}_Omon_obs_{ds.time[0].values.astype(str)[:4]}01_{ds.time[-1].values.astype(str)[:4]}12_1x1.nc')
    print(f'Saved to {out_dir + f"{rename_dict[var]}_Omon_obs_{ds.time[0].values.astype(str)[:4]}01_{ds.time[-1].values.astype(str)[:4]}12_1x1.nc"} \n')
