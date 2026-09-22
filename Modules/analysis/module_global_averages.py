import numpy as np
import xarray as xr
from numpy import meshgrid, deg2rad, gradient, sin, cos
from xarray import DataArray
from typing import Literal


def earth_radius(lat):
    '''
    calculate radius of Earth assuming oblate spheroid
    defined by WGS84
    
    Input
    ---------
    lat: vector or latitudes in degrees  
    
    Output
    ----------
    r: vector of radius in meters
    
    Notes
    -----------
    WGS84: https://earth-info.nga.mil/GandG/publications/tr8350.2/tr8350.2-a/Chapter%203.pdf
    '''
    # import numpy as np
    # from numpy import deg2rad, sin, cos

    # define oblate spheroid from WGS84
    a = 6378137
    b = 6356752.3142
    e2 = 1 - (b**2/a**2)
    
    # convert from geodecic to geocentric
    # see equation 3-110 in WGS84
    lat = deg2rad(lat)
    lat_gc = np.arctan( (1-e2)*np.tan(lat) )

    # radius equation
    # see equation 3-107 in WGS84
    r = (
        (a * (1 - e2)**0.5) 
         / (1 - (e2 * np.cos(lat_gc)**2))**0.5 
        )

    return r

def area_grid(lat, lon, mask=None):
    """
    Calculate the area of each grid cell
    Area is in square meters
    
    Input
    -----------
    lat: vector of latitude in degrees
    lon: vector of longitude in degrees
    
    Output
    -----------
    area: grid-cell area in square-meters with dimensions, [lat,lon]
    
    Notes
    -----------
    Based on the function in
    https://github.com/chadagreene/CDT/blob/master/cdt/cdtarea.m
    """

    if mask is None:
        mask = np.ones((len(lat),
                        len(lon)))  # new
    
    xlon, ylat = meshgrid(lon, lat)
    R = earth_radius(ylat)

    dlat = deg2rad(gradient(ylat, axis=0))
    dlon = deg2rad(gradient(xlon, axis=1))

    dy = dlat * R
    dx = dlon * R * cos(deg2rad(ylat))

    area = dy * mask * dx    # new
    
    xda = DataArray(
        area,
        dims=["lat", "lon"],
        coords={"lat": lat, "lon": lon},
        attrs={
            "long_name": "area_per_pixel",
            "description": "area per pixel",
            "units": "m^2",
        },
    )
    return xda


def area_weighted_avg(ds: xr.DataArray | xr.Dataset,
                      mask: xr.DataArray | xr.Dataset =None,
                      lat_name='lat',
                      lon_name='lon',
                      Hemisphere: Literal["south", "north"] | None = None,                     
                      integral=False):
    
    if mask is not None:
        mask = mask.sel(lat = ds[lat_name], lon = ds[lon_name])
        
    da_area = area_grid(ds[lat_name],
                        ds[lon_name],
                        mask)
  
    if Hemisphere is not None:
        if Hemisphere.lower() == 'south':
            da_area = da_area.where(da_area.lat<0, drop = True)
            ds = ds.where(ds.lat<0, drop = True)
        elif Hemisphere.lower() == 'north':
            da_area = da_area.where(da_area.lat>0, drop = True)
            ds = ds.where(ds.lat>0, drop = True)
    
    da_area = da_area.fillna(0)
    total_area = da_area.sum([lat_name,
                              lon_name])
    
    ds_weighted = ds*da_area #*12/1000/1e12

    if not integral:
        ds_weighted = (ds_weighted) / total_area

    ds_avg = ds_weighted.sum([lat_name,
                              lon_name])
    

    return ds_avg


