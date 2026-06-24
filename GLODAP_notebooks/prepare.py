

import sys
sys.path.insert(1, '/home/rpg002/BGC_skill')
from pathlib import Path
import dataclasses
from typing import Sequence

import numpy as np
import sys
import xarray as xr
import pandas as pd
from fnmatch import fnmatch


from modules.data_info.module_state_dict import (get_data_dict, 
                                                state_dict)

from modules.analysis.module_data_postprocessing import (extract_model_grid_within_distance, 
                                                        interp_xarray_to_dataframe, 
                                                        nanmasker, 
                                                        carbonate)





@dataclasses.dataclass
class data_dicts:
    var_list : list[str]
    experiment_list : list[str]  #['observation','assimilation', 'historical']
    model_list : list[str] #['CanESM5', 'CanESM5-CanOE']
    obs_source : dict[str] | str #'GLODAP'  ## this could be a dictionary for each variable as they might have different sources.
    data_directory : str | Path #'/space/hall7/sitestore/eccc/crd/cccma/users/rpg002/data'
    assimilation_BGC_run_id : int = None
    CanOE_assimilation_BGC_run_id : int = 1

    def __post_init__(self):

        self.info_dicts = {}

        for var in self.var_list:
            self.info_dicts[var] = {}
            for exp in self.experiment_list:
                state_dict = get_data_dict(
                                                        self.data_directory,
                                                        var,
                                                        exp,
                                                        self.assimilation_BGC_run_id,
                                                        self.CanOE_assimilation_BGC_run_id)
                if exp == 'assimilation' and self.assimilation_BGC_run_id is not None:
                    exp = f"assimilation_bgc{self.assimilation_BGC_run_id}"
                if state_dict is not None:
                   self.info_dicts[var][exp] = state_dict

    def get_obs_dicts(self) -> dict[str, state_dict]:
        obs_dicts = {}
        if 'observation' in self.experiment_list:
            
            for var in self.var_list:
                    if isinstance(self.obs_source, dict):
                        source = self.obs_source.get('var')
                    else:
                        source = self.obs_source
                        
                    if self.info_dicts[var].get('observation') is not None:
                        if self.info_dicts[var]['observation'].get(source) is not None:
                           obs_dicts[var] = self.info_dicts[var]['observation'].get(source)
                        else:
                            raise ValueError(f'{source} observations for {var} not found')
                    else:
                        print(f'{var} observations do not exist')
        
        return obs_dicts

    def get_model_dicts(self) -> dict[str, dict[str, dict[str, state_dict]]]:       

        model_dicts = {}
        for var in self.var_list:

            model_experiments = [exp for exp in self.info_dicts[var] if exp != 'observation']
            if len(model_experiments) > 0:
                model_dicts[var] = {}

                for exp in model_experiments:
                
                    model_dicts[var][exp] = {}

                    for model in self.info_dicts[var][exp]:
                        if  self.info_dicts[var][exp].get(model) is not None:
                            model_dicts[var][exp][model] = self.info_dicts[var][exp].get(model) 
                        else:
                            print(f'{model} {var} {exp} does not exist')

        return model_dicts


    def get_var_time_ranges(self, model_dicts : dict[str, dict[str, dict[str, state_dict]]],  obs_dicts : dict[str, state_dict]):

        var_ranges = {}
        for var in self.var_list:
            years_min = []
            years_max = []
            if var in model_dicts:
                for exp in model_dicts[var]:
                    for model in model_dicts[var][exp]:
                        years_min.append(model_dicts[var][exp][model].y0)
                        years_max.append(model_dicts[var][exp][model].y1)

            if  var in self.var_list:
                years_min.append(obs_dicts[var].y0)
                years_max.append(obs_dicts[var].y1)

            if len(years_min) > 0:
                var_ranges[var] = max(years_min), min(years_max)

        return var_ranges
    


    
def _load_model_data(model_dicts : dict[str, dict[str, dict[str, state_dict]]], unit_change_dics : dict, varx_dicts : dict = {}, verbose = True):
    return_mask = True
    model_mask = None
    for var in model_dicts:
        for exp in model_dicts[var]:

                varx = varx_dicts.get(var, var)

                for  model in model_dicts[var][exp]:
                    if verbose:
                        print(f'loading {model} {var} {exp}...')

                    ensemble_id = None
                    if any(['piControl' in exp , 'historical' in exp]):
                            ensemble_id = ['r1i1p2f1'] if  'CanESM5' in model else ['r1i1p1f1']
                    
                    model_mask_ = model_dicts[var][exp][model].load_data(varx, ensemble_id = ensemble_id, return_mask = return_mask, unit_change = unit_change_dics[varx])
    
                    return_mask = False 
                    if model_mask_ is not None:
                        model_mask = model_mask_     
                    
                    if verbose:
                        print('done.')


    return model_dicts, model_mask

def _load_obs_data( obs_dicts : dict[str, state_dict], varx_dicts : dict = {}, verbose = True):  
    obs_mask = {}
    for var in obs_dicts:
        varx = varx_dicts.get(var, var)
        if verbose:
            print(f'loading {var}')
        mask = obs_dicts[var].load_data(varx, return_mask = True)
        if mask is not None:
            obs_mask[var] = mask
        if verbose:
            print('done.')


    return obs_dicts, obs_mask


def _combine_model_exp(model_dicts: dict[str, dict[str, dict[str, state_dict]]]):

    model_em_dicts = {}

    for var in model_dicts:
        model_em_dicts[var] = {}
        for exp in model_dicts[var]:
            for model in model_dicts[var][exp]:
                    if 'CanOE' not in model:
                        model_exp = model + '-CMOC' + '_' + exp 
                    else:
                        model_exp = model +  '_' + exp    
                    model_em_dicts[var][model_exp]  =  model_dicts[var][exp][model]   
    
    return model_em_dicts
                                                            # dict(data = model_data[var][exp][model].squeeze(),
                                                            # color = model_dicts[var][exp][model].color,
                                                            # linestyle = model_dicts[var][exp][model].linestyle,
                                                            # marker = model_dicts[var][exp][model].marker,
                                                            # alpha = model_dicts[var][exp][model].alpha)


           

def extract_regional_mask(mask : xr.DataArray, 
                        lat_min : int,
                        lat_max : int,
                        lon_min : int,
                        lon_max : int):
    
    regions = mask.where((mask.lat>= lat_min) & (mask.lat <= lat_max))
    regions = regions.where( (mask.lon <= lon_max) | (mask.lon >= lon_min))
    regions = regions.where(regions == 1, 0)

    regions = regions.assign_coords(lat_min = lat_min)
    regions = regions.assign_coords(lat_max = lat_max)
    regions = regions.assign_coords(lon_min = lon_min)
    regions = regions.assign_coords(lon_max = lon_max)

    return regions



def prepare_data_for_analysis(var_list : list[str],
            experiment_list :list[str],
            model_list:list[str],
            obs_source : dict | str,
            data_directory : str,
            unit_change_dics : dict,
            # biomes_directory :str = None,
            assimilation_BGC_run_id: int = None,
            CanOE_assimilation_BGC_run_id : int = 1,
            nldyr = 1, 
            y0_show_cntrl = 2022,
            y1_show_cntrl = 2062,
            varx_dicts : dict = {'ntalk' : 'talk' ,  'ndissic' :'dissic', 'no3os' : 'no3', 'intno3' : 'no3', 'intchl' : 'chl'},
            verbose = True):


    ## prepare data configs

    data_info = data_dicts(var_list = var_list,
                       experiment_list = experiment_list,
                       model_list = model_list,
                       obs_source = obs_source,
                       data_directory = data_directory,
                       assimilation_BGC_run_id = assimilation_BGC_run_id,
                       CanOE_assimilation_BGC_run_id = CanOE_assimilation_BGC_run_id)



    obs_dicts = data_info.get_obs_dicts()
    model_dicts = data_info.get_model_dicts()
    var_ranges = data_info.get_var_time_ranges(model_dicts, obs_dicts)




    if verbose:
        print('======================================================= \n')
        print('observation directories: \n')
        for var in obs_dicts:
            obs_dicts[var].PrintLoc()
        
        print('\nmodel directories: \n')
        for var  in model_dicts:
            for exp in model_dicts[var]:
                print(f'{var} {exp}: ')
                for model in model_dicts[var][exp]:
                    model_dicts[var][exp][model].PrintLoc()
                    


    ## load data from memory

    model_dicts, model_mask = _load_model_data(model_dicts, unit_change_dics = unit_change_dics, varx_dicts = varx_dicts, verbose = verbose)
    obs_dicts, obs_mask = _load_obs_data(obs_dicts, varx_dicts = varx_dicts, verbose = verbose)

 
    if 'lev' in model_mask.dims:
        mask_ocean_surface  = model_mask.isel(lev = 0).drop('lev').squeeze().load()


    data_em_dicts = _combine_model_exp(model_dicts)

   

    for var in data_em_dicts:        
        for exp in data_em_dicts[var]:
            if obs_mask.get(var, None) is not None:
                data_em_dicts[var][exp].apply_nc_mask(obs_mask[var])

            if 'piControl' in exp:
                time_selection_dict = dict(time = np.arange(1, nldyr * 12 +1 ), year=slice(y0_show_cntrl ,y1_show_cntrl))
            else:
                time_selection_dict = dict(time = np.arange(1, nldyr * 12 +1 ), year=slice(*var_ranges[var]))

            data_em_dicts[var][exp].sel(time_selection_dict)
    
    for var in obs_dicts:
        if var not in data_em_dicts:
            data_em_dicts[var] = {}
        obs_dicts[var].apply_nc_mask(model_mask)
        obs_dicts[var].sel(dict( year=slice(*var_ranges[var])))

        data_em_dicts[var]['obs'] = obs_dicts[var]

    mask_ocean_surface = mask_ocean_surface.fillna(0)
    model_mask = model_mask.squeeze().fillna(0)

    for var in obs_mask:
        obs_mask[var] =  obs_mask[var].fillna(0).squeeze() 


    return data_em_dicts, obs_mask, model_mask, mask_ocean_surface





def write_model_obs_data_to_dataframe(dict_data: dict[str, dict[str, state_dict]],
                                  biomes_dict: dict,  
                                  min_count = 1,
                                  model_lev_bounds: np.ndarray | list = None):
    
    dataframe_dict = {}



    for var in dict_data:
        experiments = [exp for exp in dict_data[var] if exp != 'obs']
        if len(experiments)>0:
            first_var_w_model_data = var
            first_model_exp_for_that =  experiments[0]
            break
    
    for bms_label, mask in biomes_dict.items():
        delete_dummy = False
        print(f'{bms_label}:')
        dataframe_dict[bms_label] = {}

        for var in dict_data:
        
            if var not in dataframe_dict:
                print(f'        {var} ' )
                
                if 'obs' not in dict_data[var]:
                    print(f'no observation exists for {var}')
                    return
                
                ref = dict_data[var]['obs'].data
                ref = ref[(ref['lat'] >=  mask['lat_min'].values ) & (ref['lat'] <=  mask['lat_max'].values ) ]
                ref = ref[(ref['lon'] >=  mask['lon_min'].values ) & (ref['lon'] <=  mask['lon_max'].values ) ]
                ref = ref[~np.isnan(ref['obs'])]


                data = []
                for ds in [exp for exp in dict_data[var] if exp != 'obs']: #['historical_CMOC', 'historical_CanOE', f'{assimilation_CanOE}',f'{assimilation_CMOC}', 'piControl' , 'piControl_CanOE']:
                    data.append(dict_data[var][ds].data.assign_coords(run = ds).load())

                if len(data) > 0:
                    data = xr.concat(data, dim = 'run')
                else:
                    delete_dummy = True
                    data = xr.full_like(dict_data[first_var_w_model_data][first_model_exp_for_that].data, np.nan).expand_dims(run = 1).assign_coords(run = ['dummy']).load()
                    
                gridded = extract_model_grid_within_distance(ref, data, min_count = min_count,mask = mask, lev_bins = model_lev_bounds, tol=2,thresh=200,badval=-999999 )
            
                if model_lev_bounds is None:
                    data = data.where((data['lat'] >=  mask['lat_min'].values  - 1) & (data['lat'] <=  mask['lat_max'].values + 1), drop = True)
                    data = data.where((data['lon'] >=  mask['lon_min'].values - 1) & (data['lon'] <=  mask['lon_max'].values + 1), drop = True)
                    gridded = gridded[gridded['lev'] <= data['lev'].max().values]
                    gridded[list(data.run.values)] = interp_xarray_to_dataframe(data, gridded)
                else:
                    gridded[list(data.run.values)] = [data.sel(year = gridded['year'].values[i], time = gridded['time'].values[i], lat = gridded['lat'].values[i], lon = gridded['lon'].values[i], lev = gridded['lev'].values[i], method = 'nearest').values 
                                                    for i in range(len(gridded))]
                    
                gridded['lev'] = data['lev'].sel(lev = gridded['lev'].values, method = 'nearest').values
                if delete_dummy:
                    delete_dummy = False
                    gridded = gridded.drop(columns=["dummy"])

                dataframe_dict[bms_label][var] = gridded



    return dataframe_dict
            



def infer_carbonate_chemistry(dataframe_dict : dict, carbonate_var_list : list, silicate_climatologes_dirs : dict = None ):

        for var in carbonate_var_list:
                dataframe_dict[var] = {}


        if not all(['talk' in dataframe_dict, 
                'dissic' in dataframe_dict, 
                'po4' in dataframe_dict, 
                'no3' in dataframe_dict, 
                'silicate' in dataframe_dict, 
                'so' in dataframe_dict, 
                'thetao' in dataframe_dict]):

                raise RuntimeError('all of talk, dissic, po4, no3, silicate, so, and thetao should be available for carbonate chemistry calculation.')
        

  
        model_runs  = [i for i in list(dataframe_dict['talk'].columns) if any(['CanOE' in i, 'CMOC' in i])]

        talk = dataframe_dict['talk']
        dissic = dataframe_dict['dissic']
        thetao= dataframe_dict['thetao']
        so = dataframe_dict['so']
        pressure = None #dataframe_dict['pressure'][bms_label] 
        silicate =  infer_model_silicate(dataframe_dict['silicate'], silicate_climatologes_dirs, model_runs) 
        po4 =  infer_model_phosphate(dataframe_dict['po4'] , dataframe_dict['no3'])
        sulfide = 0 
        ammonia = 0 
        output = carbonate(carbonate_var_list, talk, dissic, thetao, so, pressure, silicate , po4 , sulfide , ammonia , temperature_out = None, pressure_out = None )
        for ind, var in enumerate(carbonate_var_list):
            dataframe_dict[var] = output[ind]

        
        return dataframe_dict


def infer_model_silicate(silicate_obs_dataframe : pd.DataFrame, silicate_climatologes_dirs : dict | str, model_runs : list):

    df = silicate_obs_dataframe.copy()

    if isinstance(silicate_climatologes_dirs, dict):  
        for model_run in model_runs:
           
            silicate = xr.open_dataset(silicate_climatologes_dirs[model_run])['silicate'] 

            df[model_run] = np.array([silicate.sel( 
                                        time = silicate_obs_dataframe['time'].values[i], 
                                        lat = silicate_obs_dataframe['lat'].values[i], 
                                        lon = silicate_obs_dataframe['lon'].values[i], 
                                        deptht = silicate_obs_dataframe['lev'].values[i], method = 'nearest').values 
                                                            for i in range(len(silicate_obs_dataframe))] )[:,None]


    else:  
            silicate = xr.open_dataset(silicate_climatologes_dirs)['silicate'] 

            df[model_runs] = np.repeat(np.array([silicate.sel( 
                                        time = silicate_obs_dataframe['time'].values[i], 
                                        lat = silicate_obs_dataframe['lat'].values[i], 
                                        lon = silicate_obs_dataframe['lon'].values[i], 
                                        deptht = silicate_obs_dataframe['lev'].values[i], method = 'nearest').values 
                                                            for i in range(len(silicate_obs_dataframe))] )[:,None], len(model_runs), axis = 1)
        
    return df




def infer_model_phosphate(po4_obs_dataframe : pd.DataFrame, no3_dataframe : pd.DataFrame):
    df = po4_obs_dataframe.copy()
    model_runs  = [i for i in list(no3_dataframe.columns) if any(['CanOE' in i, 'CMOC' in i])]
    df[model_runs] = no3_dataframe[model_runs]/16

    return df


COMMON_KEYS = ["year", "time", "lat", "lon", "lev"]




def get_common_locations(
    dataframe: pd.DataFrame,
    location_ref_dataframe: pd.DataFrame,
    keys: Sequence[str] = COMMON_KEYS,
):
    keys = list(keys)

    common = (
        location_ref_dataframe[keys]
        .drop_duplicates()
        .merge(
            dataframe[keys].drop_duplicates(),
            on=keys,
            how="inner",
        )
    )

    return (
        dataframe.merge(common, on=keys, how="inner"),
        location_ref_dataframe.merge(common, on=keys, how="inner"),
    )



def experiment_finder(dataframe: pd.DataFrame, model_experiment : list ):
    model_runs  = [i for i in list(dataframe.columns) if any(['CanOE' in i, 'CMOC' in i])]
    modelexp_list = []

    for item in model_experiment:
        if 'obs' in item:
            modelexp_list.append('obs')
        else:
            model, exp = tuple(item.split(' '))
            for model_exp in model_runs:
                if fnmatch(model_exp, f"{model}*{exp}*"):
                    modelexp_list.append(model_exp)
                    break

    return modelexp_list