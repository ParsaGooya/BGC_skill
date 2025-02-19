import numpy as np
import xarray as xr
from xarray import DataArray


def dens0(S,T):
    """
    Calculates density of seawater at atmospheric pressure and a 
    hydrostatic pressure of 0 dbar (surface).

    Parameters
    ----------
    S : xarray
           practical salinity [(PSS-78 scale)]
    T : xarray
           temperature [degree C (ITS-90)]

    Returns
    -------
    density : array_like
           density of seawater [kg m^-3]

    Usage
    --------
    >>> import eh22tools as eh
    >>> density = eh.dens0(S,T)
    if S and T are not singular they must have same dimensions

    References
    ----------
    .. [1] Unesco 1983. Algorithms for computation of fundamental properties of 
       seawater, 1983. Unesco Tech. Pap. in Mar. Sci., No. 44, 53 pp.

    .. [2] Millero, F.J. and  Poisson, A.
       International one-atmosphere equation of state of seawater.
       Deep-Sea Res. 1981. Vol28A(6) pp625-629.

    """

    # check S,T dimensions and verify they have the same shape or are singular
    if (S.shape!=T.shape):
        raise TypeError('dens0: S & T must have same dimensions or be singular')

    # Convert temperature on ITS-90 scale to IPTS-68 scale for use with density equations
    T = 1.00024 * T
    
    # Calculate density of pure water
    a0 = 999.842594
    a1 =   6.793952e-2
    a2 =  -9.095290e-3
    a3 =   1.001685e-4
    a4 =  -1.120083e-6
    a5 =   6.536332e-9
    dens_0sal = a0 + a1*T + a2*T**2 + a3*T**3 + a4*T**4 + a5*T**5;

    # Correct density for salinity
    b0 =  8.24493e-1
    b1 = -4.0899e-3
    b2 =  7.6438e-5
    b3 = -8.2467e-7
    b4 =  5.3875e-9
    c0 = -5.72466e-3
    c1 =  1.0227e-4
    c2 = -1.6546e-6
    d0 = 4.8314e-4
    return(dens_0sal + (b0 + b1*T + b2*T**2 + b3*T**3 + b4*T**4)*S + (c0 + c1*T + c2*T**2)*S**1.5 + d0*S**2)


def spatial_mask(ds):
    ds_mask = ds.where(ds.to_array().isnull,1)
    ds_mask = ds_mask.where(ds_mask == 1,0)
    return ds_mask.to_array().squeeze()


def create_train_mask(ds,
                      exclude_idx=0):
    
    xds = ds[list(ds.data_vars)[0]]
    
    mask = np.full((xds.shape[0],
                    xds.shape[1]),
                   False,
                   dtype=bool)
    x = np.arange(0,
                  12*xds.shape[0],
                  12)   
    y = np.arange(1,
                  xds.shape[1] + 1)
    idx_array = x[..., None] + y
    mask[idx_array >= idx_array[-1,
                                exclude_idx + 12]] = True
    return mask

def get_climatology(ds,
                    mask=None,
                    axis=0):
    
    xds = ds[list(ds.data_vars)[0]]
    
    if mask is None:
        mask = np.full((xds.shape[0],
                        xds.shape[1]),
                        False,
                        dtype=bool)
        
    preprocessing_mask = np.broadcast_to(mask[...,
                                              None,
                                              None],
                                         xds.shape)
    preprocessing_mask = xr.DataArray(preprocessing_mask,
                                      dims=xds.dims,
                                      coords=xds.coords)        
    ds_clim = (xds.where(~preprocessing_mask)).mean('year').to_dataset()
    
    return ds_clim



def write_monthly_to_annual(ds,
                            time='time',
                            season = 'ANN',
                            dataset=True):
    
    for jj in range(4):
        sea = [ii*12+3*jj+np.arange(3) for ii in range(0,
                                                          0 + 1 )]
        seas = list(np.stack(sea,
                             axis=0).flatten())
        if jj == 0:
            JFM = seas
        if jj == 1:
            AMJ = seas
        if jj == 2:
            JAS = seas
        if jj == 3:
            OND = seas

    seasons = {'JFM': JFM,
               'AMJ': AMJ,
               'JAS': JAS,
               'OND': OND,
               'ANN': np.arange(12)}

    for ind, month in enumerate(['Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']):
            seasons[month] = [ii*12+ ind for ii in range(0, 0 + 1 )]
    
    if dataset:
        
        dim_year = False
        if 'year' in ds.dims:
            dim_year = True
            ds = ds.rename({'year' : 'year_tmp'})
            
        ds.coords['year'] = (ds[time] // 12).astype('int')
        ds.coords['time'] = np.ceil(ds[time] % 12).astype('int') 
        

        ds_am = ds.where((ds.time >= min(seasons[season])) & (ds.time <= max(seasons[season])) , drop = True).groupby('year').mean(dim='time')
        # else:
        #     ls = []
        #     for ds_year in ds.groupby('year'):
        #         temp = ds_year[1].where((ds.time >= min(seasons[season])) & (ds.time <= max(seasons[season])) , drop = True).stack(ref = ['year_tmp', 'time'])
        #         target_time = temp.year_tmp + temp.time/len(temp.time)
        #         ls.append(temp.reset_index('year_tmp','time').drop('year').assign_coords(ref = target_time.values))
        #     ds_am = xr.concat(ls, dim = 'year').rename({'ref' : 'year_tmp'}).assign_coords(year = np.arange(len(ls)))
        if dim_year:
            ds_am = ds_am.rename({'year' : 'time'})
            ds_am = ds_am.rename({'year_tmp':'year'})
            
    if not dataset:
            ds_am = []
            for iyr in np.arange(ds.shape[0] // 12):
                ds_time = ds[iyr*12:(iyr+1)*12-1].mean()
                ds_am.append(ds_time)

          

    return ds_am.transpose('year','time',...)


def broadcast(arr,
              xds,
              arr_dim='2D'):
    
    if arr_dim == '2D':
        arr_broadcast = np.broadcast_to(arr[...,
                                            None,
                                            None],
                                        xds.shape)
    if arr_dim == '3D':
        arr_broadcast = np.broadcast_to(arr[...,
                                            None,
                                            None,
                                            None],
                                        xds.shape)
    if arr_dim == '4D':
        arr_broadcast = np.broadcast_to(arr[...,
                                            None,
                                            None,
                                            None],
                                        xds.shape)
    arr_broadcast = xr.DataArray(arr_broadcast,
                                 dims=xds.dims,
                                 coords=xds.coords).squeeze()  
    return arr_broadcast


def detrend(data,
            trend_dim,
            deg=1,
            mask=None,
            with_intercept=True):
    
    var = list(data.data_vars)[0]    
    
    data[trend_dim] = np.arange(len(data[trend_dim]))

    if mask is not None:
        trend_coefs = data.where(~mask).polyfit(dim=trend_dim,
                                                deg=deg,
                                                skipna=True)
    else:
        trend_coefs = data.polyfit(dim=trend_dim,
                                   deg=deg,
                                   skipna=True)

    slope = trend_coefs[f'{var}_polyfit_coefficients'][0].to_numpy()
    intercept = trend_coefs[f'{var}_polyfit_coefficients'][1].to_numpy()
    
    trend_axis = int(np.where(np.array(data.dims) == trend_dim)[0])
    
    timesteps = np.expand_dims(np.arange(data[var].shape[trend_axis]),
                               axis=[i for i in range(0,data[var].ndim) if i != trend_axis])
    
    slope = np.expand_dims(slope,
                           axis=trend_axis)
    
    intercept = np.expand_dims(intercept, axis=trend_axis)

    if with_intercept:
        trend = timesteps * slope + intercept
    else:
        trend = timesteps * slope
        
    data_detrend = data - trend
    
    
    trend = DataArray(trend,
                      coords=data.coords,
                      dims=data.dims).to_dataset(name=var)
    
    return data_detrend, trend, slope, intercept



def get_climatology_on_base(ds,
                            iyr_base,
                            eyr_base,
                            idate_ds='YYYY0101', # Jan 1
                            freq='mon',
                            ltime_dep=True,
                            full_set=False,
                            remove_mean = False):
    
    if not ltime_dep:  # relative to climatology of year 1
        ds_base = ds.sel(year=slice(iyr_base,
                                      eyr_base)).isel(time=np.arange(12))
        nly = len(ds.time) // 12
        ds_base = xr.concat([ds_base for _ in range(nly)],
                             dim='time').assign_coords(time=ds.time)
        ds_clim = ds_base.mean(dim='year')
        

    # could use get_climatology with proper mask instead..
        
    if  ltime_dep:     # relative to ltime dependent climatology
        
        if not full_set: # if full hindcast set not available
            print('Full data set not available for aligned climatology..')
            print('Period varies for each lead time..')
            ds_clim = ds.sel(year=slice(iyr_base,
                                        eyr_base)).mean(dim='year')
    
        if full_set: # if full hindcast set is available
        
            if freq == 'mon':  # for monthly climatologies
            
                nmth = 12
                nld_mth = ds['time'].size
                nld_yr = nld_mth // nmth

                imth_hind = int(idate_ds[4:6])
    
                ds_clim_xr = []

                for ild_mth in np.arange(nld_mth): # align and take climatology on base

                    ild_yr = ild_mth // nmth 
        
                    iyr_base_ild = iyr_base - ild_yr
                    eyr_base_ild = eyr_base - ild_yr

                    if ild_mth <= nmth - imth_hind + ild_yr*nmth:
                        iyr_base_ild_x = iyr_base_ild
                        eyr_base_ild_x = eyr_base_ild
                
                    if ild_mth > nmth - imth_hind + ild_yr*nmth:
                        iyr_base_ild_x = iyr_base_ild - 1
                        eyr_base_ild_x = eyr_base_ild - 1
                
                    ds_sliced_ild = ds.sel(year=slice(iyr_base_ild_x,
                                                      eyr_base_ild_x)).sel(time=slice(ild_mth,
                                                                                  ild_mth))

                    ds_clim_ild = ds_sliced_ild.mean(dim='year')   

                    ds_clim_xr.append(ds_clim_ild)

                ds_clim = xr.concat(ds_clim_xr,
                                    dim='time').sortby('time')
            
            
            if freq == "ann":  # for annual climatologies
        
                nld_yr = ds['time'].size

                ds_clim_xr = []
        
                for ild_yr in np.arange(nld_yr):

                    iyr_base_ild = iyr_base - ild_yr
                    eyr_base_ild = eyr_base - ild_yr

                    ds_sliced_ild = ds.sel(year=slice(iyr_base_ild,
                                                      eyr_base_ild)).sel(time=slice(ild_yr,
                                                                                ild_yr))
                
                    ds_clim_ild = ds_sliced_ild.mean(dim='year')   
            
                    ds_clim_xr.append(ds_clim_ild)
            
                ds_clim = xr.concat(ds_clim_xr,
                                    dim='time').sortby('time')
    if remove_mean:
             ds_clim = ds_clim - ds_clim.isel(time = np.arange(12)).mean('time')  
           
    return ds_clim