import numpy as np
import xarray as xr
from xarray import DataArray
import pandas as pd
import PyCO2SYS as pyco2


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

def ekman_pumping(tauu, tauv, density = None, epsilon = 1 / (60*60*24)):
    a_E =  6371e3 #m
    f = coriolis_freq(tauu.lat)
    B = Beta(tauu.lat)
    
    dtauu_dy = tauu.differentiate('lat') * 360 / (a_E * 2 * np.pi)
    dtauu_dx = tauu.differentiate('lon') * 360 / (a_E * 2 * np.pi * np.cos(np.deg2rad(tauv.lat)))
    dtauv_dy = tauv.differentiate('lat') * 360 /( a_E * 2 * np.pi )
    dtauv_dx = tauv.differentiate('lon') * 360 /( a_E * 2 * np.pi  * np.cos(np.deg2rad(tauv.lat)))


    if density == None:
        density = 1025 # kg m-3
    else:
        if 'lev' in density.dims:
            densit = density.isel(lev = 0)

    e2f2 = epsilon ** 2 + f **2
    nominator = (e2f2 * (epsilon * dtauu_dx - f * dtauu_dy + f * dtauv_dx + epsilon * dtauv_dy - B * tauu) - 2 * B * f * epsilon * tauv + 2 * (f ** 2) * B * tauu  )
    
    return (nominator / (density * (e2f2 ** 2))).transpose(..., 'lat','lon')
    # return (-dtauu_dy / (density * f)) + ((tauu * B) / (density * (f**2))) + (dtauv_dx / (density * f))

def wind_stress(wind_10m, density_air = None, Cd = None):
    if density_air == None:
        density_air = 1.225
    if Cd == None:
        Cd = 1.2e-3
    return density_air * Cd * (wind_10m * np.abs(wind_10m)) 
    
def Beta(lat):
    omega = 7.292116e-5 #rad s-1
    a_E =  6371e3 #m
    return 2 * omega * np.cos(np.deg2rad(lat))  /  a_E 

def coriolis_freq(lat):
    omega = 7.292116e-5 #rad s-1
    return 2*omega * np.sin(np.deg2rad(lat))


def spatial_mask(ds):
    ds_mask = ds.where(np.isnan(ds),1)
    ds_mask = ds_mask.where(ds_mask == 1,0)
    return ds_mask.squeeze()


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


def spco2_temp(spco2, tos):
    tos_anom = tos - tos.mean('time')
    return np.exp(0.0423* tos_anom) * spco2.mean('time')


def accumulation(ds, uo,  vo, wo, wo_levels = None, conversion = True, return_grads = False):
    if wo_levels is None:
        wo_levels = wo.lev

    a_E =  6371e3 #m

    wo = wo.interp(lev = wo_levels)
    expand_dims  = [0,1,-2,-1] if 'ensembles' not in ds.dims else [0,1,2,-2,-1]
    dndz=(ds[...,1:,:,:].data- ds[...,:-1,:,:].data)/np.expand_dims(ds.lev[1:].data-ds.lev[:-1].data, expand_dims)
    wdndz = xr.DataArray(dndz, dims = ds.dims, coords = wo[...,1:-1,:,:].coords ) * wo[...,1:-1,:,:]

    uo = uo.interp(lon = np.arange(-179,180)) 
    expand_dims  = [0,1,-3,-2] if 'ensembles' not in ds.dims else [0,1,2,-3,-2]
    expand_dims_lat = [0,1,-3,-1] if 'ensembles' not in ds.dims else [0,1,2,-3,-1]
    dndx = (ds[...,1:].data- ds[...,:-1].data) * 360 /np.expand_dims(ds.lon[1:].data-ds.lon[:-1].data, expand_dims)  / np.expand_dims(a_E * 2 * np.pi * np.cos(np.deg2rad(ds.lat)), expand_dims_lat) 
    udndx = xr.DataArray(dndx, dims = ds.dims, coords = uo.coords ) * uo

    vo = vo.interp(lat = np.arange(-89,90)) 
    expand_dims  = [0,1,-3,-1] if 'ensembles' not in ds.dims else [0,1,2,-3,-1]
    dndy = (ds[...,:,1:,:].data- ds[...,:,:-1,:].data) * 360 /np.expand_dims(ds.lat[1:].data-ds.lat[:-1].data, expand_dims) / (a_E * 2 * np.pi) 
    vdndy = xr.DataArray(dndy, dims = ds.dims, coords = vo.coords )  * vo

    if conversion:
        udndx = udndx * 1e6 
        vdndy = vdndy * 1e6
        wdndz = wdndz * 1e6
        print(r'unit coverted to umol m$^{-3}$ s$^{-1}$')
    if return_grads:
        return udndx , vdndy , wdndz , xr.DataArray(dndx, dims = ds.dims, coords = uo.coords ), xr.DataArray(dndy, dims = ds.dims, coords = vo.coords ),  xr.DataArray(dndz, dims = ds.dims, coords = wo[...,1:-1,:,:].coords )
    else:
        return udndx , vdndy , wdndz 



def annual_mean_stacked(data, dim = 'ref'):
    y0 = int(min(data[dim].values))
    y1 = int(max(data[dim].values))
    years = np.arange(y0, y1+1)
    out =  xr.concat( [data[i:i+12].mean(dim).assign_coords(year = years[ind]) for ind, i in enumerate(np.arange(0,len(data[dim]),12))], dim = 'year')
    out = out.rename({'year' : dim})
    return out





import statsmodels.nonparametric.smoothers_lowess

def _lowess_ufunc(data, lo_pts=None, lo_delta=None, it=None):

    lowess = statsmodels.nonparametric.smoothers_lowess.lowess
    
    ### If importing an xr.DataArray make numpy array
    if (type(data)==type(xr.DataArray([]))):
        data = data.values
        
    ### This adds an extra dimension if 1D
    ### Turns DataArray into numpy array
    if (len(data.shape)==1):
        data = np.expand_dims(data, axis=1)
        
    ### Get dimensions
    ndim0 = np.shape(data)[0]
    ndim1 = np.shape(data)[1]

    ### Allocate space to store data
    data_dt = np.ones((ndim0, ndim1))*np.NaN
    #slope = np.ones((ndim1))*np.NaN
    #intercept = np.ones((ndim1))*np.NaN

    ### Loop over the stacked dimension
    #for dim1 in tqdm(range(ndim1)):
    for dim1 in range(ndim1):  
        ### Mask is true if not a NaN
        mask = ~np.isnan(data[:, dim1])
        
        ### If the mask is all false
        ### We will skip that point
        if np.sum(mask)!=0:
            ### apply smoother with parameters 
            trend = lowess(data[:,dim1], 
                    np.array([x for x in range(len(data[:,dim1]))], dtype=np.float64),
                    frac=lo_pts / len(data[:,dim1]),
                    delta=lo_delta * len(data[:,dim1]),
                    it=it)
            
            ### subtract linear trend
            data_dt[mask, dim1] = trend[:,1]
            del trend

    return data_dt

#==============================================
# LOWESS FUNCTION
#=============================================
def lowess(data, dim=None, lo_pts=None, lo_delta=None, it=3):
    ### Get coordinate names
    coords = list(dict(data.coords).keys())

    ### Pop out the coordinate you want to detrend over
    coords.pop(coords.index(dim))
    try:
        coords.pop(coords.index('d'))
    except:
        pass
    ### stack the other dimensions
    if len(coords) > 0:

        data = data.stack({'z':coords})
        ## Apply detrend
        out = xr.apply_ufunc(_lowess_ufunc, data.load(), lo_pts, lo_delta, it)
        ## Unstack it
        return out.unstack('z').transpose(dim, *coords)
    
    else:
        out = xr.full_like(data, np.nan)
        out[:] = _lowess_ufunc(data.load().values, lo_pts=10, lo_delta=0.01 ,it=3).squeeze()
        return out





def fmtVarName(strx):
    """ transform string into one that meets python naming conventions
    arg: str
    returns: str suitable for file name
    """
    vName=re.sub('[^a-zA-Z0-9_\-\s/]','',strx.strip())
    vName=re.sub('[\s/]','_',vName)
    vName=re.sub('-','_',vName)
    if re.match('[0-9]',vName):
        vName='_'+vName
    return vName

def haversine(la0,lo0,la1,lo1):
    """ haversine formula with numpy array handling
    Calculates spherical distance between points on Earth in meters
    Compares elements of (la0,lo0) with (la1,lo1)
    Shapes must be compatible with numpy array broadcasting
    args: lats and lons in decimal degrees
    returns: distance on sphere with volumetric mean Earth radius in meters
    """
    rEarth=6371*1e3 # 
    # convert to radians
    la0=np.radians(la0)
    la1=np.radians(la1)
    lo0=np.radians(lo0)
    lo1=np.radians(lo1)
    theta=2*np.arcsin(np.sqrt(np.sin((la0-la1)/2)**2+np.cos(la0)*np.cos(la1)*np.sin((lo0-lo1)/2)**2))
    d=rEarth*theta
    return d

def nearest_point(la,lo,gridla,gridlo,mask=None,tol=2,thresh=40,badval=-999999):
    """ find nearest point to la,lo on arbitrary horizontal grid with points at gridla, gridlo
    with optional masking (eg for ocean values) 
    args: 
        la,lo: coordinates of point(s) to find
        gridla,gridlo: 2-d array of grid latitude, longitude
        mask: boolean mask of same shape ad gridla and gridlo with ones at valid grid points
        tol: number of degrees lat/lon within which to search; could be problematic near poles
        thresh: only accept grid points within this many km of target 
                default 40 km based on .5 deg grid and half diagonal grid distance:
                        .5*111km/deg*sqrt(2)/2 ~ 39.2 km
    returns:
        j,i model indices along y and x dimensions, respectively
        if mask is provided, j and i are set to badval; note that if negative integer, might index array w/o err
    """
    if not isinstance(la,(int,float,np.int64,np.float64,np.int32,np.float32)):
        try:
            len(la)
        except TypeError as err:
            raise Exception(f'TypeError:{err}\nAdd type to list? '+repr(type(la))) from None
        ji=np.array([nearest_point(ila,ilo,gridla,gridlo,mask,tol,thresh,badval) for ila, ilo in zip(la,lo)])
        return ji[:,0],ji[:,1]
    jjj,iii=np.where(((gridlo>lo-tol)&(gridlo<lo+tol)&(gridla>la-tol)&(gridla<la+tol))|\
                     ((gridlo+360>lo-tol)&(gridlo+360<lo+tol)&(gridla>la-tol)&(gridla<la+tol))|\
                     ((gridlo>lo+360-tol)&(gridlo<lo+360+tol)&(gridla>la-tol)&(gridla<la+tol)))
    #jjj,iii=np.where((gridlo>lo-tol)&(gridlo<lo+tol)&(gridla>la-tol)&(gridla<la+tol))
    ddd=np.minimum(np.minimum(haversine(la,lo+360,gridla[jjj,iii],gridlo[jjj,iii]),
                            haversine(la,lo,gridla[jjj,iii],gridlo[jjj,iii])),
                            haversine(la,lo,gridla[jjj,iii],gridlo[jjj,iii]+360))

    ind=np.argmin(ddd)
    if ddd[ind]<=thresh*1e3:
        j=jjj[ind]
        i=iii[ind]
        if mask is not None:
            if mask[j,i]==0:
                j=badval
                i=badval
    else:
        j=badval
        i=badval
    return   [j,i]

def nearest_valid_point(la,lo,gridla,gridlo,mask,tol=2,thresh=160,badval=-999999):
    """ find nearest point to la,lo on arbitrary horizontal grid with points at gridla, gridlo
    with optional masking (eg for ocean values) 
    args: 
        la,lo: coordinates of point(s) to find
        gridla,gridlo: 2-d array of grid latitude, longitude
        mask: boolean mask of same shape ad gridla and gridlo with ones at valid grid points
        tol: number of degrees lat/lon within which to search; could be problematic near poles
        thresh: only accept grid points within this many km of target 
                for .5 deg grid and half diagonal grid distance:
                        .5*111km/deg/(sqrt(2)/2) ~ 80 km
    returns:
        j,i model indices along y and x dimensions, respectively
        nearest valid (mask=1) point is returned
        if none is found within thresh, badval is returned; note that if negative integer, might index array w/o err

    """
    if not isinstance(la,(int,float,np.int64,np.float64,np.int32,np.float32)):
        try:
            len(la)
        except TypeError as err:
            raise Exception(f'TypeError:{err}\nAdd type to list? '+repr(type(la))) from None
        ji=np.array([nearest_valid_point(ila,ilo,gridla,gridlo,mask,tol,thresh,badval) for ila, ilo in zip(la,lo)])
        return ji[:,0],ji[:,1]
    jjj,iii=np.where((mask==1) & (((gridlo>lo-tol)&(gridlo<lo+tol)&(gridla>la-tol)&(gridla<la+tol))|\
                     ((gridlo+360>lo-tol)&(gridlo+360<lo+tol)&(gridla>la-tol)&(gridla<la+tol))|\
                     ((gridlo>lo+360-tol)&(gridlo<lo+360+tol)&(gridla>la-tol)&(gridla<la+tol))))
    #jjj,iii=np.where((gridlo>lo-tol)&(gridlo<lo+tol)&(gridla>la-tol)&(gridla<la+tol))
    if len(jjj)==0: # none found
        return [badval,badval]
    ddd=np.minimum(np.minimum(haversine(la,lo+360,gridla[jjj,iii],gridlo[jjj,iii]),
                            haversine(la,lo,gridla[jjj,iii],gridlo[jjj,iii])),
                            haversine(la,lo,gridla[jjj,iii],gridlo[jjj,iii]+360))

    ind=np.argmin(ddd)
    if ddd[ind]<=thresh*1e3:
        j=jjj[ind]
        i=iii[ind]
    else:
        j=badval
        i=badval
    return  [j,i]

def carbonate(output_list, talk, dissic, thetao, so, pressure, silicate , po4 , sulfide = 0, ammonia = 0, temperature_out = None, pressure_out = None ):
    pyco2_kws = {}
    if type(talk) == pd.core.frame.DataFrame:
        runs =  [i for i in list(talk.columns) if any(['CanOE' in i, 'CMOC' in i, 'obs' in  i])]

        keys = ['year', 'time', 'lat', 'lon' ,'lev']
        common = (
            talk[keys]
            .merge(dissic[keys], on=keys, how="inner")
            .merge(thetao[keys], on=keys, how="inner")
            .merge(so[keys], on=keys, how="inner")
            .merge(silicate[keys], on=keys, how="inner")
            .merge(po4[keys], on=keys, how="inner")
            .drop_duplicates()
        )

        # Define the known marine carbonate system parameters
        pyco2_kws["par1"] = talk.merge(common, on=keys, how="inner")[runs].values  # talk measured in the lab, Total scale
        pyco2_kws["par2"] = dissic.merge(common, on=keys, how="inner")[runs].values   # DIC measured in the lab in μmol/kg-sw
        pyco2_kws["par1_type"] = 1  # tell PyCO2SYS: "par1 is a talk value"
        pyco2_kws["par2_type"] = 2  # tell PyCO2SYS: "par2 is a DIC value"

        # Define the seawater conditions and add them to the dict
        pyco2_kws["salinity"] = so.merge(common, on=keys, how="inner")[runs].values   # practical salinity
        pyco2_kws["temperature"] = thetao.merge(common, on=keys, how="inner")[runs].values   # lab temperature (input conditions) in °C
        pyco2_kws["total_silicate"]  = silicate.merge(common, on=keys, how="inner")[runs].values  # total silicate in μmol/kg-sw
        pyco2_kws["total_phosphate"] = po4.merge(common, on=keys, how="inner")[runs].values  # total phosphate in μmol/kg-sw     

        if type(pressure) == type(talk):
            pressure = pressure.merge(common, on=keys, how="inner")[runs].values 
        else:
            pressure =  np.repeat(talk.merge(common, on=keys, how="inner")['lev'].values[None,].T, len(runs), axis=1)

        pyco2_kws["pressure"] = pressure # lab pressure (input conditions) in dbar, ignoring the atmosphere

        if type(temperature_out) == type(talk):
            temperature_out = temperature_out.merge(common, on=keys, how="inner")[runs].values
        if type(pressure_out) == type(talk):
            pressure_out = pressure_out.merge(common, on=keys, how="inner")[runs].values
        if type(sulfide) == type(talk):
            sulfide = sulfide.merge(common, on=keys, how="inner")[runs].values
        if type(ammonia) == type(talk):
            ammonia = ammonia.merge(common, on=keys, how="inner")[runs].values

        pyco2_kws["temperature_out"] = 10 if temperature_out is None else temperature_out # in-situ temperature (output conditions) in °C
        pyco2_kws["pressure_out"] = 5 if pressure_out is None else pressure_out # in-situ pressure (output conditions) in dbar, ignoring the atmosphere
        pyco2_kws["total_ammonia"] = ammonia  # total ammonia in μmol/kg-sw
        pyco2_kws["total_sulfide"] = sulfide   # total sulfide in μmol/kg-sw
        results = pyco2.sys(**pyco2_kws)

    else:

        # Define the known marine carbonate system parameters
        pyco2_kws["par1"] = talk  # talk measured in the lab, Total scale
        pyco2_kws["par2"] = dissic  # DIC measured in the lab in μmol/kg-sw
        pyco2_kws["par1_type"] = 1  # tell PyCO2SYS: "par1 is a talk value"
        pyco2_kws["par2_type"] = 2  # tell PyCO2SYS: "par2 is a DIC value"

        # Define the seawater conditions and add them to the dict
        pyco2_kws["salinity"] = so   # practical salinity
        pyco2_kws["temperature"] = thetao   # lab temperature (input conditions) in °C
        pyco2_kws["total_silicate"]  = silicate  # total silicate in μmol/kg-sw
        pyco2_kws["total_phosphate"] = po4  # total phosphate in μmol/kg-sw     
        pyco2_kws["pressure"] = pressure # lab pressure (input conditions) in dbar, ignoring the atmosphere

        pyco2_kws["temperature_out"] = 10 if temperature_out is None else temperature_out # in-situ temperature (output conditions) in °C
        pyco2_kws["pressure_out"] = 5 if pressure_out is None else pressure_out # in-situ pressure (output conditions) in dbar, ignoring the atmosphere
        pyco2_kws["total_ammonia"] = ammonia  # total ammonia in μmol/kg-sw
        pyco2_kws["total_sulfide"] = sulfide   # total sulfide in μmol/kg-sw        
    # Now calculate everything with PyCO2SYS!
        results = xr.apply_ufunc(co2sys, **pyco2_kws, input_core_dims=[[], [], [], [], []],
                                                    output_core_dims=[[]],
                                                    vectorize=True,
                                                    dask="parallelized",
                                                    output_dtypes=[float],
                                                )

    # Now calculate everything with PyCO2SYS!
    
    ls = []
    for output in output_list:
        if type(talk) == pd.core.frame.DataFrame:
            df = talk.merge(common, on=keys, how="inner").copy()
            df[runs] = results[output]
            ls.append(df)
        else:
            ls.append(results[output])
    return tuple(ls)



def co2sys(pyco2_kws):
    out =pyco2.sys(**pyco2_kws)
    return out





from tqdm import tqdm
from scipy.interpolate import griddata
# from module_data_postprocessing import nearest_valid_point


def date_coord(ds):
    data_stacked  = ds.stack(date = ('year', 'time'))
    data_stacked['date'] =  pd.to_datetime({'year': data_stacked.year.values, 'month':data_stacked.time.values , 'day': 15}).values
    return data_stacked.transpose('date', ...)

def interp_xarray_to_dataframe(data: xr.DataArray  , ref:  pd.core.frame.DataFrame):
    interpolated_list = []
    for _, row in tqdm(ref.iterrows()):
        if 'date' in data.dims:
            interpolated_list.append(data.interp(date = row['date'], lat = row['lat'], lon = row['lon'], lev = row['lev'], kwargs={"fill_value": "extrapolate"}).values) 
        else:
            interpolated_list.append(data.interp(year = row['year'], time = row['time'],  lat = row['lat'], lon = row['lon'], lev = row['lev'], kwargs={"fill_value": "extrapolate"}).values) 

    return  interpolated_list


def interp_xarray_to_dataframe_vertorized(data:xr.DataArray, ref:pd.core.frame.DataFrame):
    if 'date' in data.dims:
        epoch = np.datetime64('1980-01-01T00:00:00.000000000')
        T, Z, Y, X = np.meshgrid(
            (data['date'].values - epoch) / np.timedelta64(1, 'D') , 
            data['lev'].values,
            data['lat'].values,
            data['lon'].values,
            indexing='ij')
        data_points = np.column_stack([T.ravel(), Z.ravel(), Y.ravel(), X.ravel()])
        ref_date_pd = pd.to_datetime(ref.date.values)
        ref_date_np = ref_date_pd.values.astype('datetime64[ns]')
        ref_dates = (ref_date_np - epoch) / np.timedelta64(1, 'D')
        interp_points = np.column_stack([ref_dates, ref['lev'].values, ref['lat'].values, ref['lon'].values])

    else:
        
        Yr, T, Z, Y, X = np.meshgrid(
            data['year'].values,
            data['time'].values,
            data['lev'].values,
            data['lat'].values,
            data['lon'].values,
            indexing='ij')

        data_points = np.column_stack([Yr.ravel(), T.ravel(), Z.ravel(), Y.ravel(), X.ravel()])
        interp_points = np.column_stack([ref['year'].values, ref['time'].values, ref['lev'].values, ref['lat'].values, ref['lon'].values])

    data_values =  np.stack([data[i].values.ravel() for i in range(data.shape[0])], axis=1)  
    data_points = data_points[~np.isnan(data_values[:,0])] 
    data_values = data_values[~np.isnan(data_values[:,0])]

    vals = griddata(data_points, data_values, interp_points, method='linear')

    return vals


def grid_obs_dataframe(ref, min_count = 3):

    lat_min = np.floor(ref['lat'].min())
    lat_max = np.ceil(ref['lat'].max())
    bins = np.arange(lat_min  , lat_max + 1, 1)
    ref['lat_bin'] = pd.cut(ref['lat'], bins=bins, right=False) 

    lon_min = np.floor(ref['lon'].min())
    lon_max = np.ceil(ref['lon'].max())
    bins = np.arange(lon_min , lon_max + 1, 1)
    ref['lon_bin'] = pd.cut(ref['lon'], bins=bins, right=False) 
    gridded = ref.groupby(['year','time', 'lat_bin', 'lon_bin', 'lev']).agg(
    obs=('obs', 'mean'),
    var=('obs', 'var'),
    count=('obs', 'count') ).reset_index()

           
    gridded = gridded[gridded['count'] >= min_count]
    lat = [(gridded['lat_bin'].values[i].left +  gridded['lat_bin'].values[i].right) * 0.5 for i in range(len(gridded['lat_bin']) )]
    lon = [(gridded['lon_bin'].values[i].left +  gridded['lon_bin'].values[i].right) * 0.5 for i in range(len(gridded['lon_bin']) )]

    gridded = gridded.rename(columns = {'lat_bin' : 'lat', 'lon_bin' : 'lon'})
    gridded['lat'] = lat
    gridded['lon'] =lon

    return gridded

def extract_model_grid_within_distance(ref, data, min_count = 3,mask = None, lev_bins = None, tol=2,thresh=160,badval=-999999 ):
    
    model_lat = data['lat']
    model_lon = data['lon']

    if len(model_lat.shape) == 1:
        Y, X = np.meshgrid(
                model_lat['lat'].values,
                model_lon['lon'].values,
                indexing='ij')
    if mask is not None:
        assert Y.shape == mask.shape
    else:
        mask =1 
    
    model_lat_ind, model_lon_ind = nearest_valid_point(ref['lat'], ref['lon'] ,Y, X, mask = mask,tol=tol,thresh=thresh,badval=badval )
    ref = ref[model_lat_ind > badval]

    model_lat_ind = model_lat_ind[model_lat_ind>badval]
    model_lon_ind = model_lon_ind[model_lon_ind>badval]
    ref['model_lat'] = Y[model_lat_ind, model_lon_ind]
    ref['model_lon'] = X[model_lat_ind, model_lon_ind]


    if lev_bins is not None:
        ref['lev_bins'] = pd.cut(ref['lev'], bins=lev_bins, right=False) 
    else:
        ref['lev_bins'] = ref['lev']

    gridded = ref.groupby(['year','time', 'model_lat', 'model_lon', 'lev_bins']).agg(
    obs=('obs', 'mean'),
    var=('obs', 'var'),
    count=('obs', 'count') ).reset_index()

    gridded = gridded[gridded['count'] >= min_count]
    if lev_bins is not None:
        lev = [(gridded['lev_bins'].values[i].left +  gridded['lev_bins'].values[i].right) * 0.5 for i in range(len(gridded['lev_bins']) )]
        gridded['lev_bins'] = lev

    gridded = gridded.rename(columns = {'lev_bins' : 'lev', 'model_lat' : 'lat', 'model_lon' : 'lon'})

    return gridded




def nanmasker(ds: xr.DataArray, dim = 'time', return_mask = False, min_valid_Fraction = 0.8):
    fraction = 1 - min_valid_Fraction
    mask = ds.isnull().sum(dim=dim)/len(ds[dim])
    if return_mask:
        return ds * xr.ones_like(mask).where(mask<fraction, np.nan), xr.ones_like(mask).where(mask<fraction, np.nan)
    else:
        return ds * xr.ones_like(mask).where(mask<fraction, np.nan)


def add_cyclic_point(ds):
    add = ds.isel(lon = -1).assign_coords(lon = ds.isel(lon = -1).lon.values + 1)
    return xr.concat([ds,add], dim = 'lon')