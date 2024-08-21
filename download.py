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
############################## read addresses to OCCCI data, edit the links and download to local directory #########################

out_dir = '/space/hall5/sitestore/eccc/crd/ccrn/users/rpg002/data/chlos/observations/raw/' ## path to local directory where you want your observational data saved

lines = []
with open('links.txt', 'r') as file:  ## read addresses
    for line in file:
        lines.append(line.strip())

def address_update(line): ### update addresses
    line = line.split('occci-v1.0')[0] + 'occci-v6.0' + line.split('occci-v1.0')[1]
    line = line.split('OC4v6')[0] + 'OCx' + line.split('OC4v6')[1]
    line = line.split('fv1')[0] + 'fv6' + line.split('fv1')[1]
    return line

lines = (list(set(lines)))

for line in lines: ### download and save the data with concise naming
    for m in range(12):
        downloaded = glob.glob(out_dir + '*.nc')
        dl_link = address_update(line)
        dl_link = dl_link.split('-fv6.0.nc')[0].split('OCx-')
        dl_link = dl_link[0] + 'OCx-' + f'{int(dl_link[1]) + m}' + '-fv6.0.nc' 
        print(f'\n downloading {dl_link.split("/")[-1]} .... \n')
        try:
            if (out_dir + dl_link.split('/')[-1]) not in downloaded:
                wget.download(dl_link, out=out_dir + dl_link.split('/')[-1])
            else:
                print('Already downloaded.')
        except:
            print('file not existant')


