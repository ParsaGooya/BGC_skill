
import glob

parent_data_dir = '/space/hall7/sitestore/eccc/crd/cccma/users/rpg002/data'
biomes_dir = "/space/hall7/sitestore/eccc/crd/cccma/users/rsa001/big_store/AI/FGCO2/data/biomes"
output_dir = "/space/hall7/sitestore/eccc/crd/cccma/users/rpg002/output"


def dir_generator(parent_data_dir , varlist, experiment_list, model_list, verbose = True ):
    dirs = {}
    for model in model_list:
        dirs[model] = {}
        for var in varlist:
            dirs[model][var] = {}
            if var == 'intchl':
                raise RuntimeError('For intchl load CHL at this step and later integrate it when loading the data in the memory! ')
            for exp in experiment_list:
                if exp =='observation':
                    path_ = f'{parent_data_dir}/{var}/{exp}'
                elif 'CanOE' in exp:
                    path_ = f'{parent_data_dir}/{var}/{exp.split("_")[0]}/{model}-CanOE*'
                else:
                    path_ = f'{parent_data_dir}/{var}/{exp.split("_")[0]}/{model}'

                if len(glob.glob(path_ + '/')) > 0:
                    dirs[model][var][exp] = path_
                else:
                    if verbose:
                        print(f'{model} {var} {exp} does not exist')

    dirs['biomes'] = biomes_dir
    dirs['output'] = output_dir
    return dirs


