import numpy as np
import xarray as xr
from xarray import DataArray



def rmse(ds1,ds2, dim = 'year'):
    return ((ds1 - ds2)**2).mean(dim)**(1/2)

def corr(ds1,ds2, dim = 'year'):
        mean1 = ds1.mean(dim)
        mean2 = ds2.mean(dim)
        std1  = ds1.std(dim)
        std2  = ds2.std(dim)
        covariance = ((ds1 - mean1) * (ds2 - mean2)).mean(dim)
        return (covariance/(std1 * std2))
    
def corr_map(ds1, ds2):
        mean1 = ds1.mean(['lat',
                          'lon'])
        mean2 = ds2.mean(['lat',
                          'lon'])
        std1  = ds1.std(['lat',
                         'lon'])
        std2  = ds2.std(['lat',
                         'lon'])
        covariance = ((ds1 - mean1) * (ds2 - mean2)).mean(['lat',
                                                           'lon'])
        return (covariance/(std1 * std2)) #.values

    
def get_variance(ds,
                 var='tas',
                 dim_ens='ensemble',
                 dim_year='year',
                 check=False):

    nens  = len(ds[dim_ens].values)
    nens1 = nens - 1

    var_emean = ds.mean(dim=dim_ens).var(dim=dim_year,ddof=1)
    residual  = ds - ds.mean(dim=dim_ens)

    var_total = ds.var(dim=dim_year,ddof=1).mean(dim=dim_ens) 
    var_noise = residual.var(dim=dim_ens,ddof=1).mean(dim=dim_year)
    var_predc = var_emean - var_noise/float(nens)

    # var_noise_alternative = nens*(var_total - var_emean)/float(nens1)
    var_total_alternative = var_predc + var_noise

    if check:
        print(f"Min Tvar Method 1: {float(var_total.to_array().min())}")
        print(f"Min TVar Method 2: {float(var_total_alternative.to_array().min())}")
        print(f"Max TVar Method 1: {float(var_total.to_array().max())}")
        print(f"Max TVar Method 2: {float(var_total_alternative.to_array().max())}")
        
        
    var_total = var_total.rename({var:'total'})
    var_noise = var_noise.rename({var:'noise'})
    var_predc = var_predc.rename({var:'predictable'})

    ds_variance = var_total.merge(var_predc.merge(var_noise,
                                                  join='exact'),
                                  join='exact')

    
    return ds_variance


