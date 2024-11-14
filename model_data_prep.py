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


####################################################### specify : ###############################################################################################

model = 'CanESM5' # CanESM5 or CanESM5-CanOE
var = 'thetao' 
lev_range = 600 ## if 3D

assimilation = True
extracted_from_disc = True
extention_years = [2021,2022,2023]

hindcast = False
hindcat_initial_year = 1975
hindcat_final_year = 2025

simulation = False
simulation_initial_year = 1950
simulation_final_year = 2025

print(f'Model : {model}')
print(f'Variable : {var}')
print(f'vertical range : {lev_range} \n\n')


###################################################################################################################################################################
def coords_edit(ds):
    lat = ds.lat.values[:,0]
    lon = ds.lon.values[0,:]
    ds_like = xr. DataArray(ds.values, dims = ds.dims).rename({'x' :'lon', 'y': 'lat'}).assign_coords({'lat': lat, 'lon': lon })
    for d in ds.dims:
        if d not in ['x','y']:
            ds_like = ds_like.assign_coords(d = ds[d])
    return ds_like
#######################################################################  Hindcast ##################################################################################
if hindcast:

    print(f' Loading Hindcast data : {hindcat_initial_year} - 2019 \n')

    out_dir = f'/space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/{var}/forecast/{model}/' ## specify local directory
    dir_hindcast = f'/space/hall5/sitestore/eccc/crd/ccrn/model_output/CMIP6/final/CMIP6/DCPP/CCCma/{model}/dcppA-hindcast' ### dir to dcppA_hindcast

    realization  = []    ### extract available realizations
    for dir in glob.glob(dir_hindcast + '/*'):  
        name = dir.split('/')[-1]
        realization.append(name.split('-')[-1])
    realization = np.unique(realization)
    print(f'Available realizations - Hindcast : \n {realization} \n')

    ############# clear output directoy ###############
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    ####### load data based on intialization year #####
    hindcast_dict = {}
    for year in tqdm(range(hindcat_initial_year,2020)):
        dsets = {}
        for r in realization:
            try:
                path = glob.glob(dir_hindcast + '/' + f's{year}-' + r + f'/Omon/{var}/gn/v20190429/*nc')
                ds = xr.open_dataset(path[0])
                if 'lev' in ds.dims:
                    ds = ds.where(ds.lev <= lev_range, drop =True )
                dsets[r] = ds
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
        ds = ds.assign_coords(year = key+1).transpose('ensembles','lead_time',...,'lat','lon')
        ds.to_dataset(name = var).to_netcdf(out_dir + f'{var}_Omon_{key}_ensmebles_{key+1}01_{key+10}12_1x1_LE.nc')


    ##############################################  Forecast #####################################
    if hindcat_final_year > 2019:

        print(f' Loading Forecast data : 2020 - {hindcat_final_year}  \n')

        dir_hindcast = f'/space/hall5/sitestore/eccc/crd/ccrn/model_output/CMIP6/final/CMIP6/DCPP/CCCma/{model}/dcppB-forecast'

        realization  = []
        for dir in glob.glob(dir_hindcast + '/*'):
            name = dir.split('/')[-1]
            realization.append(name.split('-')[-1])
        realization = np.unique(realization)
        print(f'Available realizations - Forecast : \n {realization}')

        hindcast_dict = {}
        for year in tqdm(range(2020,hindcat_final_year)):
            dsets = {}
            for r in realization:
                try:
                    path = glob.glob(dir_hindcast + '/' + f's{year}-' + r + f'/Omon/{var}/gn/v20190429/*nc')
                    ds = xr.open_dataset(path[0])
                    if 'lev' in ds.dims:
                        ds = ds.where(ds.lev <= lev_range, drop =True )
                    dsets[r] = ds
                except:
                    pass
            hindcast_dict[year] = xr.concat([ds for key, ds in dsets.items()], dim = 'member').assign_coords(member = list(dsets.keys()) )

        ds_out = xe.util.grid_global(1, 1)
        for key, ds in tqdm(hindcast_dict.items()):
            
            regridder = xe.Regridder(ds, ds_out, 'bilinear', ignore_degenerate=True, periodic=True)
            ds = coords_edit(regridder(ds[var]))
            ds = ds.rename({'time': 'lead_time', 'member' : 'ensembles'})
            ds['lead_time'] = np.arange(1,121)
            ds = ds.assign_coords(year = key+1).transpose('ensembles','lead_time',...,'lat','lon')
            ds.to_dataset(name = var).to_netcdf(out_dir + f'{var}_Omon_{key}_ensmebles_{key+1}01_{key+10}12_1x1_LE.nc')


    print(' Hindcast data saved! \n\n')
##################################################################################################################
########################################### assimilation data ####################################################
if assimilation: 
    print(f' Loading assimilation data : 1958 - ? (see what is available on disc) \n runs beyond that should be extrtacted from disc first. see extract_from_disc.sh \n\n ')

    out_dir = f'/space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/{var}/assimilation/{model}/' ### local directory
    dir_assimilation = f'/space/hall5/sitestore/eccc/crd/ccrn/model_output/CMIP6/final/CMIP6/DCPP/CCCma/{model}/dcppA-assim'

    ############# clear output directoy ###############
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    ##################################################
    realization  = []
    for dir in glob.glob(dir_assimilation + '/*'):
        name = dir.split('/')[-1]
        realization.append(name.split('-')[-1])
    realization = np.unique(realization)
    print(f'Available realizations - Assimilation : \n {realization} \n')
    ###### extract data based on realization ######

    dsets = {}
    for r in realization:
        try:
            ds = xr.open_mfdataset(str(Path(dir_assimilation + '/'  + r + f'/Omon/{var}/gn/v20190429/', "*.nc")), combine='nested', concat_dim='time')
            if 'lev' in ds.dims:
                    ds = ds.where(ds.lev <= lev_range, drop =True )
            dsets[r] = ds 
        except:
            pass
    ###### concat and regrid ##### 

    assim = xr.concat([ds for key, ds in dsets.items()], dim = 'member').assign_coords(member = list(dsets.keys()) )

    ds_out = xe.util.grid_global(1, 1) 
    regridder = xe.Regridder(assim, ds_out, 'bilinear', ignore_degenerate=True, periodic=True)
    assim = coords_edit(regridder(assim[var]))
    assim = assim.rename({'member' : 'ensembles'}).transpose('ensembles','time',...,'lat','lon')
    assim.to_dataset(name = var).to_netcdf(out_dir + f'{var}_Omon_ensmebles_{assim.time[0].values.item().year}01_{assim.time[-1].values.item().year}12_1x1_LE.nc')

    print(' Assimilation data saved! \n\n')
######################## assimilation data for 2021-2023 should be extracted from disc first which are saved based on year. see: extract_from_disc.sh ################################
    if extracted_from_disc:
        print(f' Loading assimilation data : {extention_years}\n ')
        ls_year = []
        for year in extention_years:
            ls = glob.glob(f'/space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/{var}/assimilation/{model}/extentions/*{year}*.nc')
            ensembles = [link.split('_')[4] for link in ls] ### extract ensemble id from the name of the file
            ls_year.append(xr.concat([xr.open_dataset(link) for link in ls], dim = 'ensembles').assign_coords(ensembles = ensembles))  
        
        assim_ext = xr.concat(ls_year, dim = 'time')
        if 'lev' in assim_ext.dims:
              assim_ext = assim_ext.where(assim_ext.lev <= lev_range, drop =True )
        ds_out = xe.util.grid_global(1, 1) 
        regridder = xe.Regridder(assim_ext, ds_out, 'bilinear', ignore_degenerate=True, periodic=True)
        assim_ext = coords_edit(regridder(assim_ext[var])).assign_coords(ensembles = ensembles)
        assim_ext.to_dataset(name = var).to_netcdf(out_dir + f'{var}_Omon_ensmebles_{assim_ext.time[0].values.item().year}01_{assim_ext.time[-1].values.item().year}12_1x1_LE.nc')
########################################################################################################################
################################################### simulation data ####################################################
if simulation:
    print(f' Loading simulation data : {simulation_initial_year} - {simulation_final_year} \n\n')

    out_dir = f'/space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/{var}/simulation/{model}/' ### local directory
    dir_simmulation = f'/space/hall5/sitestore/eccc/crd/ccrn/model_output/CMIP6/final/CMIP6/CMIP/CCCma/{model}/historical'
    dir_simmulation_ssp245 = f'/space/hall5/sitestore/eccc/crd/ccrn/model_output/CMIP6/final/CMIP6/ScenarioMIP/CCCma/{model}/ssp245'

    ############# clear output directoy ###############
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    ##################################################

    realization  = []
    for dir in glob.glob(dir_simmulation + '/*'):
        name = dir.split('/')[-1]
        if 'p1' not in name.split('-')[-1]:
            realization.append(name.split('-')[-1])
    realization = np.unique(realization)
    # realization = ['r10i1p2f1' ,'r1i1p2f1', 'r2i1p2f1' ,'r3i1p2f1' ,'r4i1p2f1' ,'r5i1p2f1', 'r6i1p2f1', 'r7i1p2f1' ,'r8i1p2f1' ,'r9i1p2f1']

    if simulation_final_year > 2015:
        realization_ssp245  = []
        for dir in glob.glob(dir_simmulation_ssp245 + '/*'):
            name = dir.split('/')[-1]
            if 'p1' not in name.split('-')[-1]:
                if name.split('-')[-1] in realization:
                    realization_ssp245.append(name.split('-')[-1])
        realization_ssp245 = np.unique(realization_ssp245)

        for element in realization:
            if element not in realization_ssp245:
                realization.remove(element)

    print(f'Available realizations - Historical : \n {realization} \n')

    ###### extract data based on realization ######

    dsets = {}
    for r in realization:
        try:
            ds = xr.open_mfdataset(str(Path(dir_simmulation + '/'  + r + f'/Omon/{var}/gn/v20190429/', "*.nc")), combine='nested', concat_dim='time').sel(time = slice(str(simulation_initial_year), None))
            if 'lev' in ds.dims:
                ds = ds.where(ds.lev <= lev_range, drop = True)
            dsets[r] = ds
        except:
            pass
    ###### concat and regrid ##### 

    sim = xr.concat([ds for key, ds in dsets.items()], dim = 'member').assign_coords(member = list(dsets.keys()) )

    ds_out = xe.util.grid_global(1, 1) 
    regridder = xe.Regridder(sim, ds_out, 'bilinear', ignore_degenerate=True, periodic=True)
    sim = coords_edit(regridder(sim[var]))
    sim = sim.rename({'member' : 'ensembles'}).transpose('ensembles','time',...,'lat','lon')

    ###### condisering saving historical and ssp245 runs separately if the dataset is superlarge (e.g. 3D):#######
    sim.to_dataset(name = var).to_netcdf(out_dir + f'{var}_Omon_ensmebles_{sim.time[0].values.item().year}01_{sim.time[-1].values.item().year}12_1x1_LE.nc')
    ###### extract data based on realization ######

    if simulation_final_year > 2015:
        dsets = {}
        for r in realization:
            try:
                ds = xr.open_mfdataset(str(Path(dir_simmulation_ssp245 + '/'  + r + f'/Omon/{var}/gn/v20190429/', "*.nc")), combine='nested', concat_dim='time')
                if 'lev' in ds.dims:
                    ds = ds.where(ds.lev <= lev_range, drop = True)
                dsets[r] = ds
            except:
                pass
        ###### concat and regrid ##### 
        sim_ssp = xr.concat([ds.sel(time = slice(None, str(simulation_final_year))) for key, ds in dsets.items()], dim = 'member').assign_coords(member = list(dsets.keys()) )

        ds_out = xe.util.grid_global(1, 1) 
        regridder = xe.Regridder(sim_ssp, ds_out, 'bilinear', ignore_degenerate=True, periodic=True)
        sim_ssp = coords_edit(regridder(sim_ssp[var]))
        sim_ssp = sim_ssp.rename({'member' : 'ensembles'}).transpose('ensembles','time',...,'lat','lon')
        ###### condisering saving historical and ssp245 runs separately if the dataset is superlarge (e.g. 3D):#######
        sim_ssp.to_dataset(name = var).to_netcdf(out_dir + f'{var}_Omon_ensmebles_{sim_ssp.time[0].values.item().year}01_{sim_ssp.time[-1].values.item().year}12_1x1_LE.nc')
        ##############################################################################################################
        # sim = xr.concat([sim.sel(ensembles = sim_ssp.ensembles), sim_ssp], dim = 'time')
        # sim.to_dataset(name = var).to_netcdf(out_dir + f'{var}_Omon_ensmebles_{sim.time[0].values.item().year}01_{sim.time[-1].values.item().year}12_1x1_LE.nc')
    else:
        sim.to_dataset(name = var).to_netcdf(out_dir + f'{var}_Omon_ensmebles_{sim.time[0].values.item().year}01_{sim.time[-1].values.item().year}12_1x1_LE.nc')

    print(' Simulation data saved! ')

#################################################### CanOE BGC new assim runs ##################################################################################
# list = glob.glob('/space/hall6/sitestore/eccc/crd/ccrn/users/scrd107/canesm_runs/d2k-asm-*')
# for l in list:
    # member = int(l.split('-')[-1][1:])
    # ds = xr.open_mfdataset(l + '/data/nc_output/CMIP6/CCCma/CCCma/*/dcppA-assim/*' + f'/Omon/{var}/gn/v20190429/' + "*.nc", combine='nested', concat_dim='time')