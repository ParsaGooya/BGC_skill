# supress warnings
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
from tqdm import tqdm


def coords_edit(ds):
    lat = ds.lat.values[:,0]
    lon = ds.lon.values[0,:]
    ds_like = xr. DataArray(ds.values, dims = ds.dims).rename({'x' :'lon', 'y': 'lat'}).assign_coords({'lat': lat, 'lon': lon })
    if 'member' in ds.dims:
        ds_like = ds_like.assign_coords({'member' : ds.member})
    if 'time' in ds.dims:
        ds_like = ds_like.assign_coords({'time' : ds.time})
    return ds_like

var = 'intpp' #chlos
##############################################  Hindcast #####################################

out_dir = f'/space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/{var}/forecast/' ## specify local directory
dir_hindcast = '/space/hall5/sitestore/eccc/crd/ccrn/model_output/CMIP6/final/CMIP6/DCPP/CCCma/CanESM5/dcppA-hindcast' ### dir to dcppA_hindcast

realization  = []    ### extract available realizations
for dir in glob.glob(dir_hindcast + '/*'):  
    name = dir.split('/')[-1]
    realization.append(name.split('-')[-1])
realization = np.unique(realization)
print(f'Available realizations : \n {realization}')

############# clear output directoy ###############
Path(out_dir).mkdir(parents=True, exist_ok=True)
####### load data based on intialization year #####
hindcast_dict = {}
for year in tqdm(range(1998,2020)):
    dsets = {}
    for r in realization:
        try:
            path = glob.glob(dir_hindcast + '/' + f's{year}-' + r + f'/Omon/{var}/gn/v20190429/*nc')
            dsets[r] = xr.open_dataset(path[0])
        except:
            pass
    hindcast_dict[year] = xr.concat([ds for key, ds in dsets.items()], dim = 'member').assign_coords(member = list(dsets.keys()) )
###### regrid data based on initialization year #####
ds_out = xe.util.grid_global(1, 1)  ### definde regrider
for key, ds in tqdm(hindcast_dict.items()):
    regridder = xe.Regridder(ds, ds_out, 'bilinear', ignore_degenerate=True, periodic=True)
    ds = coords_edit(regridder(ds[var]))
    ds = ds.rename({'time': 'lead_time', 'member' : 'ensembles'})
    ds['lead_time'] = np.arange(1,121)
    ds = ds.assign_coords(year = key+1).transpose('ensembles','lead_time','lat','lon')
    ds.to_dataset(name = var).to_netcdf(out_dir + f'{var}_Omon_{key}_ensmebles_{key+1}01_{key+10}12_1x1_LE.nc')


##############################################  Forecast #####################################
dir_hindcast = '/space/hall5/sitestore/eccc/crd/ccrn/model_output/CMIP6/final/CMIP6/DCPP/CCCma/CanESM5/dcppB-forecast'

realization  = []
for dir in glob.glob(dir_hindcast + '/*'):
    name = dir.split('/')[-1]
    realization.append(name.split('-')[-1])
realization = np.unique(realization)
print(f'Available realizations : \n {realization}')

hindcast_dict = {}
for year in tqdm(range(2020,2024)):
    dsets = {}
    for r in realization:
        try:
            path = glob.glob(dir_hindcast + '/' + f's{year}-' + r + f'/Omon/{var}/gn/v20190429/*nc')
            dsets[r] = xr.open_dataset(path[0])
        except:
            pass
    hindcast_dict[year] = xr.concat([ds for key, ds in dsets.items()], dim = 'member').assign_coords(member = list(dsets.keys()) )

ds_out = xe.util.grid_global(1, 1)
for key, ds in tqdm(hindcast_dict.items()):
    
    regridder = xe.Regridder(ds, ds_out, 'bilinear', ignore_degenerate=True, periodic=True)
    ds = coords_edit(regridder(ds[var]))
    ds = ds.rename({'time': 'lead_time', 'member' : 'ensembles'})
    ds['lead_time'] = np.arange(1,121)
    ds = ds.assign_coords(year = key+1).transpose('ensembles','lead_time','lat','lon')
    ds.to_dataset(name = var).to_netcdf(out_dir + f'{var}_Omon_{key}_ensmebles_{key+1}01_{key+10}12_1x1_LE.nc')
##################################################################################################################
########################################### assimilation data ####################################################
out_dir = f'/space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/{var}/assimilation/' ### local directory
dir_assimilation = '/space/hall5/sitestore/eccc/crd/ccrn/model_output/CMIP6/final/CMIP6/DCPP/CCCma/CanESM5/dcppA-assim'

############# clear output directoy ###############
Path(out_dir).mkdir(parents=True, exist_ok=True)
##################################################
realization  = []
for dir in glob.glob(dir_assimilation + '/*'):
    name = dir.split('/')[-1]
    realization.append(name.split('-')[-1])
realization = np.unique(realization)
print(f'Available realizations : \n {realization}')
###### extract data based on realization ######

dsets = {}
for r in realization:
    try:
        dsets[r] = xr.open_mfdataset(str(Path(dir_assimilation + '/'  + r + f'/Omon/{var}/gn/v20190429/', "*.nc")), combine='nested', concat_dim='time')
    except:
        pass
###### concat and regrid ##### 
assim = xr.concat([ds for key, ds in dsets.items()], dim = 'member').assign_coords(member = list(dsets.keys()) )

ds_out = xe.util.grid_global(1, 1) 
regridder = xe.Regridder(assim, ds_out, 'bilinear', ignore_degenerate=True, periodic=True)
assim = coords_edit(regridder(assim[var]))
assim = assim.rename({'member' : 'ensembles'}).transpose('ensembles','time','lat','lon')
assim.to_dataset(name = var).to_netcdf(out_dir + f'{var}_Omon_ensmebles_{assim.time[0].values.item().year}01_{assim.time[-1].values.item().year}12_1x1_LE.nc')


#### assimilation data for 2021-2023 should be extracted from disc first which are saved based on year #####
ls_year = []
for year in [2021,2022,2023]:
    ls = glob.glob(f'/space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/{var}/assimilation/extentions/*{year}*.nc')
    ensembles = [link.split('_')[4] for link in ls] ### extract ensemble id from the name of the file
    ls_year.append(xr.concat([xr.open_dataset(link) for link in ls], dim = 'ensembles').assign_coords(ensembles = ensembles))   
assim_ext = xr.concat(ls_year, dim = 'time')
ds_out = xe.util.grid_global(1, 1) 
regridder = xe.Regridder(assim_ext, ds_out, 'bilinear', ignore_degenerate=True, periodic=True)
assim_ext = coords_edit(regridder(assim_ext[var])).assign_coords(ensembles = ensembles)
assim_ext.to_dataset(name = var).to_netcdf(out_dir + f'{var}_Omon_ensmebles_{assim_ext.time[0].values.item().year}01_{assim_ext.time[-1].values.item().year}12_1x1_LE.nc')

#### assimilation data for 2021-2023 should be extracted from disc first which are saved based on year #####
ls_year = []
for year in [2021,2022,2023]:
    ls = glob.glob(f'/space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/{var}/assimilation/extentions/*{year}*.nc')
    ensembles = [link.split('_')[4] for link in ls] ### extract ensemble id from the name of the file
    ls_year.append(xr.concat([xr.open_dataset(link) for link in ls], dim = 'ensembles').assign_coords(ensembles = ensembles))   
assim_ext = xr.concat(ls_year, dim = 'time')
ds_out = xe.util.grid_global(1, 1) 
regridder = xe.Regridder(assim_ext, ds_out, 'bilinear', ignore_degenerate=True, periodic=True)
assim_ext = coords_edit(regridder(assim_ext[var])).assign_coords(ensembles = ensembles)
assim_ext.to_dataset(name = var).to_netcdf(out_dir + f'{var}_Omon_ensmebles_{assim_ext.time[0].values.item().year}01_{assim_ext.time[-1].values.item().year}12_1x1_LE.nc')

##################################################################################################################
########################################### simulation data ####################################################
out_dir = f'/space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/{var}/simulation/' ### local directory
dir_simmulation = '/space/hall5/sitestore/eccc/crd/ccrn/model_output/CMIP6/final/CMIP6/CMIP/CCCma/CanESM5/historical'

############# clear output directoy ###############
Path(out_dir).mkdir(parents=True, exist_ok=True)
##################################################
realization  = []
for dir in glob.glob(dir_simmulation + '/*'):
    name = dir.split('/')[-1]
    if 'p1' not in name.split('-')[-1]:
        realization.append(name.split('-')[-1])
realization = np.unique(realization)
print(f'Available realizations : \n {realization}')

###### extract data based on realization ######

dsets = {}
for r in realization:
    try:
        dsets[r] = xr.open_mfdataset(str(Path(dir_simmulation + '/'  + r + f'/Omon/{var}/gn/v20190429/', "*.nc")), combine='nested', concat_dim='time')
    except:
        pass
###### concat and regrid ##### 
sim = xr.concat([ds for key, ds in dsets.items()], dim = 'member').assign_coords(member = list(dsets.keys()) )

ds_out = xe.util.grid_global(1, 1) 
regridder = xe.Regridder(sim, ds_out, 'bilinear', ignore_degenerate=True, periodic=True)
sim = coords_edit(regridder(sim[var]))
sim = sim.rename({'member' : 'ensembles'}).transpose('ensembles','time','lat','lon')
sim.to_dataset(name = var).to_netcdf(out_dir + f'{var}_Omon_ensmebles_{sim.time[0].values.item().year}01_{sim.time[-1].values.item().year}12_1x1_LE.nc')
