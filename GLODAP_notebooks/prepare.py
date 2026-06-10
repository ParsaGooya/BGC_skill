
import warnings
warnings.filterwarnings('ignore')
import sys
sys.path.insert(1, '/home/rpg002/BGC_skill')
from pathlib import Path
import dataclasses

import numpy as np
import sys
import xarray as xr
from modules.data_info.module_state_dict import get_data_dict, state_dict
from modules.analysis.module_data_load import load_nc_data, load_csv_data, load_biomes


from modules.data_info.module_data_info import biomes_module
from modules.analysis.module_data_load import load_biomes

from modules.analysis.module_data_postprocessing import nanmasker 






@dataclasses.dataclass
class data_dicts:
    var_list : list[str]
    experiment_list : [str]  #['observation','assimilation', 'historical']
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
    return regions



def prepare_data_for_analysis(var_list : list[str],
            experiment_list :list[str],
            model_list:list[str],
            obs_source : dict | str,
            data_directory : str,
            biomes_directory :str,
            unit_change_dics : dict,
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


    bms_info = biomes_module(biomes_directory)
    bms_dict = bms_info.data_dict['biomes']
    bms_dir  = bms_dict['dir']
    bms_file = bms_dict['file']
    bms_dict_plot = bms_info.dict_biomes_plot


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
                    
        print('\nbiomes direcoty: \n')
        print('======================================================= \n')

    ## load data from memory

    model_dicts, model_mask = _load_model_data(model_dicts, unit_change_dics = unit_change_dics, varx_dicts = varx_dicts, verbose = verbose)
    obs_dicts, obs_mask = _load_obs_data(obs_dicts, varx_dicts = varx_dicts, verbose = verbose)

 
    if 'lev' in model_mask.dims:
        mask_ocean_surface  = model_mask.isel(lev = 0).drop('lev').squeeze().load()


    biomes = load_biomes(f'{bms_dir}',
                        f'{bms_file}') * mask_ocean_surface
    nbms = 17
    bms_labels = [bms_info.biomes_dict[ii]['label'] for ii in np.arange(nbms)+1]

    mask_biomes = {}
    for bms_label in bms_labels:
        bm = (biomes.MeanBiomes).where(biomes.MeanBiomes  == bms_labels.index(bms_label)+1, 0)
        mask_biomes[bms_label] = bm.where(bm == 0,1) +  (biomes.MeanBiomes ).where(np.isnan((biomes.MeanBiomes)),0)


    ## prepare datasets


    data_em_dicts = _combine_model_exp(model_dicts)

   

    for var in data_em_dicts:        
        for exp in data_em_dicts[var]:
            if obs_mask.get(var, None) is not None:
                data_em_dicts[var][exp].apply_mask(obs_mask[var])

            if 'piControl' in exp:
                time_selection_dict = dict(time = np.arange(1, nldyr * 12 +1 ), year=slice(y0_show_cntrl ,y1_show_cntrl))
            else:
                time_selection_dict = dict(time = np.arange(1, nldyr * 12 +1 ), year=slice(*var_ranges[var]))

            data_em_dicts[var][exp].sel(time_selection_dict)
    
    for var in obs_dicts:
        if var not in data_em_dicts:
            data_em_dicts[var] = {}
        obs_dicts[var].apply_mask(model_mask)
        obs_dicts[var].sel(dict( year=slice(*var_ranges[var])))

        data_em_dicts[var]['obs'] = obs_dicts[var]

    mask_ocean_surface = mask_ocean_surface.fillna(0)
    model_mask = model_mask.squeeze().fillna(0)

    for var in obs_mask:
        obs_mask[var] =  obs_mask[var].fillna(0).squeeze() 



    return data_em_dicts, obs_mask, model_mask, mask_ocean_surface, mask_biomes