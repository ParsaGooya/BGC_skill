import numpy as np
import xarray as xr


def align_data_and_targets(data,
                           targets):

    lead_years = data.time.size // 12
    last_target_year = data.year[-1].item() + lead_years - 1
    last_years_diff = data.year[-1].item() - targets.year[-1].item()
    year_diff = targets.year[-1].item() - last_target_year
    if year_diff >= 0:
        if 'ensembles' in data.dims:
            ds = data[:,:,:(lead_years*12),...]
        else:
            ds = data[:,:(lead_years*12),...]
        data_out = targets[:len(targets.year) - year_diff,...]
    elif (year_diff< 0) & (last_years_diff <= 0):
        if 'ensembles' in data.dims:
            ds = data[:,:,:(lead_years*12),...]
        else:
            ds = data[:,:(lead_years*12),...]
        data_out = targets
    else:
        if 'ensembles' in data.dims:
            ds = data[:-last_years_diff,:,:(lead_years*12),...]
        else:
            ds = data[:-last_years_diff,:(lead_years*12),...]
        data_out = targets

    return ds, data_out  



def reshape_data_to_reference(data,
                              data_ref): 
    """
     input : ds     : time      x lat x lon
     ouput : ds_out : year x 12 x lat x lon
    """
    lead_years = int(data_ref.shape[1] / 12)
    ls = [data[y :y + lead_years,...].stack(flat=['year',
                                                 'time']).reset_index('flat',
                                                                      drop = True) for y in range(len(data.year))]
    data_reshaped = xr.concat([ds.assign_coords(flat=np.arange(len(ds.flat))) for ds in ls],
                             dim='year').assign_coords(year=data.year)
    data_reshaped = data_reshaped.rename({'flat':'time'}).transpose('year',
                                                                    'time',
                                                                  ...).where((data_reshaped.year >= data_ref.year.min()) & (data_reshaped.year <= data_ref.year.max()),
                                                                             drop = True)    
    return data_reshaped



def rewrite_data_like_hindcasts(dat_in,
                                hnd_in,
                                check=False):
    '''
       input :  dat_in : year*12 x lat x lon 
                hnd_in : year x time x lat x lon
       output: targets : year x time x lat x lon
    '''
    hnd_var, dat_var = align_data_and_targets(hnd_in,
                                              dat_in)
    dat_out =  reshape_data_to_reference(dat_var,
                                         hnd_var)
    if check: # check it works as intended
        for yr in np.arange(dat_out.year.values[-dat_out.time.values.shape[0]//12],
                            dat_out.year.values[-1]):
            print(yr)
            print(dat_out.sel(lat=0,
                                   method='nearest').sel(lon=0,
                                                         method='nearest').sel(year=yr).values)
            print("====")
    return dat_out

def rewrite_data_like_obs_on_target(data, data_ref):
    time_min = max(data.year[0].values , data_ref.year[0].values)
    time_max = min(data.year[-1].values + np.divmod(data.time[-1].values,12)[0] , data_ref.year[-1].values)
    data_target =  xr.concat([data[:,i].assign_coords(year = data.year.values + np.divmod(i,12)[0] ) for i in range(len(data.time))], dim = 'time').transpose('year','time',...)
    data_target = data_target.sel(year = slice(time_min, time_max))
    data_ref_target = data_ref.sel(year = slice(time_min, time_max))
    return data_target, data_ref_target
##

def xa_empty(lon,
             lat,
             time,
             year):
    data = np.empty((len(year),
                     len(time),
                     len(lat),
                     len(lon)))
    data[:] = np.nan
    xa = xr.Dataset(
                    coords=dict(
                                time=time,
                                year=year,
                                lon=lon,
                                lat=lat,
                                ),
                    data_vars=dict(
                                   data=(["year", "time", "lat", "lon"], data),
                                  ),
                    ).to_array().squeeze()
    return xa


def rewrite_adj_data_like_hindcasts(dat,
                                    hnd,
                                    var = 'chlos',
                                    check=False):
    '''
       input :  dat_in : year_a x time_a x lat x lon 
                hnd_in : year_b x time_b x lat x lon
       output:  dat_out: year_b x time_b x lat x lon (fill in nan)
    '''
    xa = xa_empty(hnd.lon,
                  hnd.lat,
                  hnd.time,
                  hnd.year)
    dt = np.max([0,dat.year[0].values - hnd.year.values[0]]) 
    for iyr in np.arange(dat.year.size):
        if hnd.year.values[0]+iyr+dt <= hnd.year.values[-1]:
            xa[iyr+dt,0:dat.time.size,:,:] = dat.to_array().squeeze()[iyr,:,:,:]
    dat_out = xa.to_dataset(name=var)
    if check: # check it works as intended
        for yr in np.arange(dat_out.year.values[-dat_out.time.values.shape[0]//12],
                            dat_out.year.values[-1]):
            print(yr)
            print(dat_out[var].sel(lat=0,
                                   method='nearest').sel(lon=0,
                                                         method='nearest').sel(year=yr).values)
            print("====")
    return dat_out


def align_data_to_common_base(data,
                              nldyr):
    '''
    rewrites data from initial years to target years
    on the same base period for each lead year
    '''

    datasets = []

    for ildyr in np.arange(nldyr):
        
        ildyr1 = ildyr + 1
            
        inds = np.arange(ildyr*12,(ildyr1)*12)
            
        years_ild = data.year.values - ildyr # common period per lead year

        years_ild = years_ild[years_ild>=data.year.values.min()]    

        ds = data.sel(time=inds).sel(year=years_ild)
            
        ds['year'] = years_ild +  ildyr       # target years

        datasets.append(ds)

    data_aligned = xr.concat(datasets,
                             dim='time').sortby('time')
    
    return data_aligned



def write_monthly_to_annual(ds,
                            time='time',
                            dataset=True):
    
    if dataset:
        
        dim_year = False
        if 'year' in ds.dims:
            dim_year = True
            ds = ds.rename({'year' : 'year_tmp'})
    
        ds.coords['month'] = np.ceil(ds[time] % 12).astype('int') 
        ds.coords['year'] = (ds[time] // 12).astype('int')
        ds_am = ds.groupby('year').mean(dim='time')
    
        if dim_year:
            ds_am = ds_am.rename({'year' : 'time'})
            ds_am = ds_am.rename({'year_tmp':'year'})
            
    if not dataset:
        
        ds_am = []
        for iyr in np.arange(ds.shape[0] // 12):
            ds_time = ds[iyr*12:(iyr+1)*12-1].mean()
            ds_am.append(ds_time)
            

    return ds_am
