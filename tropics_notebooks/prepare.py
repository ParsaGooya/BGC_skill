

import sys
sys.path.insert(1, '/home/rpg002/BGC_skill')
from pathlib import Path
import dataclasses
import glob
import numpy as np
import sys
import xarray as xr
import pandas as pd
from fnmatch import fnmatch
import yaml
import copy

from modules.data_info.module_state_dict import (get_data_dict, 
                                                state_dict)

from modules.analysis.module_data_postprocessing import (get_climatology_on_base, 
                                                        get_detrended, 
                                                        apply_lowess as apply_lowess_,
                                                        spco2_temp,
                                                        dens0,
                                                        ekman_pumping,
                                                        wind_stress)

from modules.analysis.module_global_averages import area_weighted_avg

from modules.plotting.utils import Var, Exp, Model, Obs, Biome

with open(Path("/home/rpg002/BGC_skill/configs") / "GLODAP_climatology.yaml", "r") as f:
    GLODAP_clim_config = yaml.safe_load(f)

with open(Path("/home/rpg002/BGC_skill/configs") / "model_climatology.yaml", "r") as f:
    model_climatology = yaml.safe_load(f)



@dataclasses.dataclass
class data_dicts:
    var_list : list[Var]
    experiment_list : list[Exp]  #['observation','assimilation', 'historical']
    model_list : list[Model] #['CanESM5', 'CanESM5-CanOE']
    obs_source : dict[Obs] | Obs #'GLODAP'  ## this could be a dictionary for each variable as they might have different sources.
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

    def get_obs_dicts(self) -> dict[Var, state_dict]:
        obs_dicts = {}
        if 'observation' in self.experiment_list:
            
            for var in self.var_list:
 
                if self.info_dicts[var].get('observation') is not None:
                    if isinstance(self.obs_source, dict):        
                        source = self.obs_source[var]
                    else:
                        source = self.obs_source
                
                    if self.info_dicts[var]['observation'].get(source) is not None:
                        obs_dicts[var] = self.info_dicts[var]['observation'].get(source)
                    else:
                        raise ValueError(f'{source} observations for {var} not found')
        
        return obs_dicts

    def get_model_dicts(self) -> dict[Var, dict[Exp, dict[Model, state_dict]]]:       

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


        return model_dicts


    def get_var_time_ranges(self, model_dicts : dict[Var, dict[Exp, dict[Model, state_dict]]],  obs_dicts : dict[Var, state_dict]):

        var_ranges = {}
        for var in self.var_list:
            years_min = []
            years_max = []
            if var in model_dicts:
                for exp in model_dicts[var]:
                    for model in model_dicts[var][exp]:
                        years_min.append(model_dicts[var][exp][model].y0)
                        years_max.append(model_dicts[var][exp][model].y1)

            if  var in obs_dicts:
                years_min.append(obs_dicts[var].y0)
                years_max.append(obs_dicts[var].y1)

            if len(years_min) > 0:
                var_ranges[var] = max(years_min), min(years_max)

        return var_ranges
    


    
def _load_model_data(model_dicts : dict[Var, dict[Exp, dict[Model, state_dict]]], 
                     unit_change_dics : dict[Var, str], 
                     varx_dicts : dict[Var, str] = {}, 
                     target_levels: np.ndarray | xr.DataArray | None = None,
                     verbose = True):
    
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
                    
                    model_mask_ = model_dicts[var][exp][model].load_data(varx, 
                                                                         ensemble_id = ensemble_id, 
                                                                         target_levels = target_levels,
                                                                         return_mask = return_mask, 
                                                                         unit_change = unit_change_dics.get(varx))
    
                    if model_mask is None:
                        model_mask = model_mask_  
                    else:
                        if model_mask_ is not None:
                            if ("lev" in model_mask_.dims
                            and "lev" not in model_mask.dims):
                                model_mask = model_mask_     
                    
                    if verbose:
                        print('done.')

    try:
        model_mask = model_mask.drop_vars("d")
        return model_dicts, model_mask
    except :
        return model_dicts, model_mask
    

    

def _load_obs_data( obs_dicts : dict[Var, state_dict], 
                   varx_dicts : dict[Var, str] = {}, 
                   target_levels: np.ndarray | xr.DataArray | None = None,
                   verbose = True):  
    
    obs_mask = {}
    for var in obs_dicts:
        varx = varx_dicts.get(var, var)
        if verbose:
            print(f'loading {var}')
        mask = obs_dicts[var].load_data(varx, target_levels = target_levels, return_mask = True)
        if mask is not None:
            obs_mask[var] = mask
        if verbose:
            print('done.')


    return obs_dicts, obs_mask


def _combine_model_exp(model_dicts: dict[Var, dict[Exp, dict[Model, state_dict]]]):

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
    regions = regions.where( (mask.lon <= lon_max) & (mask.lon >= lon_min))
    regions = regions.where(regions == 1, 0)

    regions = regions.assign_coords(lat_min = lat_min)
    regions = regions.assign_coords(lat_max = lat_max)
    regions = regions.assign_coords(lon_min = lon_min)
    regions = regions.assign_coords(lon_max = lon_max)

    return regions



def prepare_data_for_analysis(var_list : list[Var],
            experiment_list :list[Exp],
            model_list:list[Model],
            obs_source : dict[Var, Obs] | Obs,
            data_directory : str,
            unit_change_dics : dict,
            assimilation_BGC_run_id: int = None,
            CanOE_assimilation_BGC_run_id : int = 1,
            target_levels: np.ndarray | xr.DataArray | None = None,
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

    model_dicts, model_mask = _load_model_data(model_dicts, unit_change_dics = unit_change_dics, varx_dicts = varx_dicts, verbose = verbose, target_levels = target_levels)
    obs_dicts, obs_mask = _load_obs_data(obs_dicts, varx_dicts = varx_dicts, verbose = verbose, target_levels = target_levels)

 
    if 'lev' in model_mask.dims:
        mask_ocean_surface  = model_mask.isel(lev = 0).drop('lev').squeeze().load()
    else:
        mask_ocean_surface = model_mask


    data_em_dicts = _combine_model_exp(model_dicts)

    tolerance = {
        "lat": 0.5,
        "lon": 0.5,
        "lev": 5.0,
    }

    for var in data_em_dicts:        
        for exp in data_em_dicts[var]:
            if obs_mask.get(var, None) is not None:
                data_em_dicts[var][exp].apply_nc_mask(obs_mask[var], tolerance)

            if 'piControl' in exp:
                time_selection_dict = dict(month = np.arange(1, nldyr * 12 +1 ), year=slice(y0_show_cntrl ,y1_show_cntrl))
            else:
                time_selection_dict = dict(month = np.arange(1, nldyr * 12 +1 ), year=slice(*var_ranges[var]))

            data_em_dicts[var][exp].sel(time_selection_dict)
    
    for var in obs_dicts:
        if var not in data_em_dicts:
            data_em_dicts[var] = {}
        obs_dicts[var].apply_nc_mask(model_mask, tolerance)
        obs_dicts[var].sel(dict( year=slice(*var_ranges[var])))

        data_em_dicts[var]['obs'] = obs_dicts[var]

    mask_ocean_surface = mask_ocean_surface.fillna(0)
    model_mask = model_mask.squeeze().fillna(0)

    for var in obs_mask:
        obs_mask[var] =  obs_mask[var].fillna(0).squeeze() 


    return data_em_dicts, obs_mask, model_mask, mask_ocean_surface

                 


def load_ONI(experiment_list: list[Exp],
            model_list: list[Model],
            assimilation_BGC_run_id: int | None = None,
            CanOE_assimilation_BGC_run_id: int | None = 1,
            obs_source= 'SODA'):
    _path = '/home/rpg002/BGC_skill/ONI'
    ONI_dict = {}
    
    for exp in experiment_list:
        if 'obs' in exp:
            model_exp = 'obs'
            data_path = glob.glob(f'{_path}/{obs_source}*_ONI_*.nc')[0]
            ONI_dict[model_exp] = xr.open_dataarray(data_path)

        else:
            for model in model_list:
                if 'CanOE' not in model:
                    model += '-CMOC' 
                    if 'assim' in exp and assimilation_BGC_run_id is not None:
                        model += f'_{assimilation_BGC_run_id}'
                elif 'assim' in exp and CanOE_assimilation_BGC_run_id is not None:
                    model += f'_{CanOE_assimilation_BGC_run_id}'

                model_exp = model + '_' + exp    
                _phys_model = model.split('-')[0]
                data_path = glob.glob(f'{_path}/{_phys_model}_{exp}*_ONI_*.nc')[0]
                ONI_dict[model_exp] = xr.open_dataarray(data_path)


    return ONI_dict

def get_enso_events(ONI_dict: dict[Exp, xr.DataArray],
                   lower_percentile: float = 25,
                   upper_percentike: float = 75):

    ElNino = {}
    LaNina = {}

    for exp, ONI in ONI_dict.items():
        ONI_stacked = ONI.stack(ref = ('year','month'))
        time_coords = ONI_stacked.year + (ONI_stacked.month - 0.5) /12
        ONI_stacked = ONI_stacked.assign_coords(ref = time_coords.values)

        ElNino[exp] = ONI_stacked.where((ONI_stacked >= np.percentile(ONI.dropna("month"),upper_percentike)), drop = True).ref
        LaNina[exp] = ONI_stacked.where((ONI_stacked <= np.percentile(ONI.dropna("month"),lower_percentile)), drop = True).ref

    return ElNino, LaNina


def spco2_decomposition(dict_em_data: dict[Var, dict[Exp, state_dict]]):

    if any(['spco2' not in dict_em_data and 'tos' not in dict_em_data]):
        print('spco2 decomposition not successful. Both spco2 and tos variables must exist.')
        return
    

    dict_em_data['spco2_temp'] = copy.copy(dict_em_data['spco2'])
    dict_em_data['spco2_res'] = copy.copy(dict_em_data['spco2'])

    dict_em_data['spco2_temp'].data = spco2_temp(dict_em_data['spco2'].data, dict_em_data['tos'].data)
    dict_em_data['spco2_res'].data = dict_em_data['spco2'].data - dict_em_data['spco2_temp'].data
    return dict_em_data


def get_ekman(dict_em_data: dict[Var, dict[Exp, state_dict]]):
    if not all(['uas' in dict_em_data, 
                'vas' in dict_em_data]):
        print("uas or vas does not exist. wE cannot be calculated.")
        return dict_em_data

    dict_em_data['wind_speed'] = {}
    dict_em_data['wE'] = {}

    for exp in dict_em_data['uas']:
        uas = dict_em_data['uas'][exp]
        if exp not in dict_em_data['vas']:
            pass
        vas = dict_em_data['vas'][exp]

        uas_data, vas_data = xr.align(
                uas.data,
                vas.data,
                join="inner",
            )

        wind_speed = np.sqrt((uas_data**2)+ (vas_data ** 2))
        wE = ekman_pumping(wind_stress(uas_data), wind_stress(vas_data))

        dict_em_data["wind_speed"][exp] = state_dict(
            var = "wind_speed",
            experiment = uas.experiment,
            model_key= uas.model_key,
            data = wind_speed,
            assimilation_BGC_run_id= uas.assimilation_BGC_run_id,
            CanOE_assimilation_BGC_run_id= uas.CanOE_assimilation_BGC_run_id
        )

        dict_em_data["wE"][exp] = state_dict(
            var = "wE",
            experiment = uas.experiment,
            model_key= uas.model_key,
            data = wE,
            assimilation_BGC_run_id= uas.assimilation_BGC_run_id,
            CanOE_assimilation_BGC_run_id= uas.CanOE_assimilation_BGC_run_id
        )

    return dict_em_data

def get_density(dict_em_data: dict[Var, dict[Exp, state_dict]]):
    if not all(['so' in dict_em_data, 
                'thetao' in dict_em_data]):
        print("so or thetao does not exist. Density cannot be calculated.")
        return dict_em_data

    dict_em_data['density'] = {}

    for exp in dict_em_data['thetao']:
        thetao = dict_em_data['thetao'][exp]
        if exp not in dict_em_data['so']:
            pass
        so = dict_em_data['so'][exp]

        thetao_data, so_data = xr.align(
                thetao.data,
                so.data,
                join="inner",
            )

        density = dens0(so_data, thetao_data)

        dict_em_data["density"][exp] = state_dict(
            var = "density",
            experiment = thetao.experiment,
            model_key= thetao.model_key,
            data = density,
            assimilation_BGC_run_id= thetao.assimilation_BGC_run_id,
            CanOE_assimilation_BGC_run_id= thetao.CanOE_assimilation_BGC_run_id
        )

    return dict_em_data


def calculate_climatology(
        dict_em_data: dict[Var, dict[Exp, state_dict]],
        y0_base: int | None = None,
        y1_base: int | None = None,
        center_on_zero: bool = False,
        return_anomaly: bool = False,
        verbose: bool = True
):
    dict_clim = {}
    dict_anom = {}

    for var in dict_em_data:
        dict_clim[var] = {}
        dict_anom[var] = {}

        if y0_base is None:
            y0_base_list = [dict_em_data[var][model_exp].y0 for model_exp in dict_em_data[var]]
        if y1_base is None:
            y1_base_list = [dict_em_data[var][model_exp].y1 for model_exp in dict_em_data[var]]      

        y0_base_ = y0_base or max(y0_base_list)
        y1_base_ = y1_base or min(y1_base_list)

        for model_exp in dict_em_data[var]:

            if verbose:
                print(f'calculating {var} {model_exp} climatology over {y0_base_} - {y1_base_} ...')

            dict_clim[var][model_exp] = copy.copy(dict_em_data[var][model_exp])
            dict_clim[var][model_exp].data = get_climatology_on_base( dict_em_data[var][model_exp].data,
                                                        y0_base_,
                                                        y1_base_,
                                                        center_on_zero = center_on_zero).load()
            

            if return_anomaly:
                dict_anom[var][model_exp] = copy.copy(dict_em_data[var][model_exp]) 
                dict_anom[var][model_exp].data = dict_em_data[var][model_exp].data - dict_clim[var][model_exp].data

            if verbose:
                print('done.')
    if return_anomaly:
        return dict_clim, dict_anom

    return dict_clim


def calculate_detrended(
        dict_em_data: dict[Var, dict[Exp, state_dict]],
        variable_list: list[Var] = None,
        y0_base: int | None = None,
        y1_base: int | None = None,
        month_specific_det: bool = True,
        apply_lowess: bool = False,
        lo_pts=120 , 
        lo_delta=0.01, 
        it=3,
        verbose: bool = True
):
    dict_det = {}

    if variable_list is None:
        return

    for var in variable_list:
        if var not in dict_em_data:
            raise ValueError(
                f"{var} does not exist!"
            )
        dict_det[var] = {}

        if y0_base is None:
            y0_base_list = [dict_em_data[var][model_exp].y0 for model_exp in dict_em_data[var]]
        if y1_base is None:
            y1_base_list = [dict_em_data[var][model_exp].y1 for model_exp in dict_em_data[var]]        

        y0_base_ = y0_base or max(y0_base_list)
        y1_base_ = y1_base or min(y1_base_list)

        for model_exp in dict_em_data[var]:

            if verbose:
                print(f'calculating {var} {model_exp} detrended over {y0_base_} - {y1_base_} ...')

            dict_det[var][model_exp] = copy.copy(dict_em_data[var][model_exp])
            data = copy.deepcopy(dict_em_data[var][model_exp].data.sel(year = slice(y0_base_, y1_base_)))
            detrended = get_detrended( data, month_specific_det = month_specific_det).load()
                                                        
            if apply_lowess:
                detrended = apply_lowess_(detrended, lo_pts=lo_pts , lo_delta=lo_delta, it=it)

            dict_det[var][model_exp].data = detrended

            if verbose:
                print('done.')

    return dict_det

def mask_NESO_events(
        dict_em_data: dict[Var, dict[Exp, state_dict]],
        ONI_dict: dict[Exp, xr.DataArray],
        y0_base: int = None,
        y1_base: int = None,
        calculate_mean: bool = True,
        upper_percentile: float = 75,
        lower_percentile: float = 25,
        return_diff: bool = True
):
    dict_LaNina = copy.deepcopy(dict_em_data)
    dict_ElNino = copy.deepcopy(dict_em_data)
    if return_diff:
        dict_diff = copy.deepcopy(dict_em_data)

    for var in dict_em_data:

        if y0_base is None:
            y0_base_list = [dict_em_data[var][model_exp].y0 for model_exp in dict_em_data[var]]
        if y1_base is None:
            y1_base_list = [dict_em_data[var][model_exp].y1 for model_exp in dict_em_data[var]]   


        y0_base_ = y0_base or max(y0_base_list)
        y1_base_ = y1_base or min(y1_base_list)

        for exp in dict_em_data[var]:
            if exp not in ONI_dict:
                raise ValueError(
                    f"No ONI dataset for {exp}."
                )

            data = copy.deepcopy(dict_em_data[var][exp].data.sel(year = slice(y0_base_, y1_base_)))

            ONI, data = xr.align(
                ONI_dict[exp],
                data,
                join="inner",
            )

            lanina = data.where(ONI<=np.percentile(ONI_dict[exp].dropna('month'),lower_percentile))
            elnino = data.where(ONI>=np.percentile(ONI_dict[exp].dropna('month'),upper_percentile))   

            dict_LaNina[var][exp].y0 = lanina.year.min().values
            dict_LaNina[var][exp].y1 = lanina.year.max().values

            dict_ElNino[var][exp].y0 = elnino.year.min().values
            dict_ElNino[var][exp].y1 = elnino.year.max().values

            if calculate_mean:
                lanina = lanina.mean(['year', 'month'])     
                elnino = elnino.mean(['year', 'month'])   

            if return_diff:
                diff = lanina.mean(['year', 'month']) - elnino.mean(['year', 'month'])
                dict_diff[var][exp].data = diff
                dict_diff[var][exp].y0 = lanina.year.min().values
                dict_diff[var][exp].y1 = lanina.year.max().values              

            dict_LaNina[var][exp].data = lanina
            dict_ElNino[var][exp].data = elnino

    if return_diff:
        return dict_LaNina, dict_ElNino, dict_diff
    
    return dict_LaNina, dict_ElNino


def take_lat_average(dict_em_data: dict[Var, dict[Exp, state_dict]],
                      EQ_masks_dict: dict[Biome, xr.DataArray | xr.Dataset]):
    
    dict_em_data_regional = {}
    for var in dict_em_data:
        dict_em_data_regional[var] = {}

        for region, mask in EQ_masks_dict.items():
            dict_em_data_regional[var][region] = {}

            for model_exp in dict_em_data[var]:
                state_dict = copy.deepcopy( dict_em_data[var][model_exp])
                
                state_dict.data = state_dict.data.where(mask != 0, drop = True).mean('lat')  
                if region == 'PEQ':
                    state_dict.data['lon'] = np.mod(state_dict.data['lon'] + 360,360)
                    state_dict.data = state_dict.data.sortby("lon")
                dict_em_data_regional[var][region][model_exp] = state_dict


    return dict_em_data_regional


def mask_EQ(dict_em_data: dict[Var, dict[Exp, state_dict]],
            EQ_masks_dict: dict[Biome, xr.DataArray | xr.Dataset]):
    
    dict_em_data_regional = {}
    for var in dict_em_data:
        dict_em_data_regional[var] = {}

        for region, mask in EQ_masks_dict.items():
            dict_em_data_regional[var][region] = {}

            for model_exp in dict_em_data[var]:
                state_dict = copy.deepcopy( dict_em_data[var][model_exp])
                
                state_dict.data = state_dict.data.where(mask != 0, drop = True)     
                if region == 'PEQ':
                    state_dict.data['lon'] = np.mod(state_dict.data['lon'] + 360,360)
                    state_dict.data = state_dict.data.sortby("lon")         
                dict_em_data_regional[var][region][model_exp] = state_dict

    return dict_em_data_regional




def experiment_finder(data_dict: dict[Exp, state_dict], model_experiment : list[Exp] ):
    model_runs  = [i for i in data_dict if any(['CanOE' in i, 'CMOC' in i])]
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


def create_eq_biomes(biomes: xr.Dataset, lat_max = 2, lat_min = -2):
    

    PEQ = biomes.MeanBiomes.where((biomes.MeanBiomes.lat>= lat_min) & (biomes.MeanBiomes.lat <= lat_max))
    PEQ = PEQ.where( (PEQ.lon <= -59) | (PEQ.lon >= 121))
    PEQ = PEQ.where(np.isnan(PEQ), 1)

    PEQ_E = PEQ.where(PEQ.lon<0)
    PEQ_W = PEQ.where(PEQ.lon>0)

    AEQ = biomes.MeanBiomes.where((biomes.MeanBiomes.lat>= lat_min) & (biomes.MeanBiomes.lat <= lat_max))
    AEQ = AEQ.where( (AEQ.lon >= -55) & (AEQ.lon <= 20))
    AEQ = AEQ.where(np.isnan(AEQ), 1)


    return { 'PEQ' : PEQ.where(~np.isnan(PEQ),0), 
            'AEQ' : AEQ.where(~np.isnan(AEQ),0), 
            'PEQ_E' : PEQ_E.where(~np.isnan(PEQ_E),0), 
            'PEQ_W' : PEQ_W.where(~np.isnan(PEQ_W),0)}