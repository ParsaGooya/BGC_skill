import numpy as np
import xarray as xr
from pathlib import Path
import pandas as pd


def coords_edit(ds):
    
    lat = ds.lat.values[:,0]
    lon = ds.lon.values[0,:]
    ds_like = xr. DataArray(ds.values, dims = ds.dims).rename({'x' :'lon', 'y': 'lat'}).assign_coords({'lat': lat, 'lon': lon })

    if 'member' in ds.dims:
        ds_like = ds_like.assign_coords({'member' : ds.member})
    if 'time' in ds.dims:
        ds_like = ds_like.assign_coords({'time' : ds.time})
    return ds_like


def load_biomes(dir_in,
                key='*',
                verbose=True):
    if verbose:
        print("loading global ocean biomes..")
    for f in sorted(Path(dir_in).glob(f"{key}.nc")):
        ds = xr.open_dataset(f)
        ds = ds.transpose('lat',
                          'lon',
                          'year')
    if verbose:
        print("done")
    return ds



def load_nc_data(files : list[Path], 
              ensemble_mean = True,
              ensemble_id = None,
              rename_dict : dict = None):


    times = (xr.open_mfdataset(files, combine = 'nested', concat_dim = 'time')).transpose('time', ...).time 
    # try:
    #     years = [int(np.datetime_as_string(times[i].values)[:4])  for i in range(0,len(times),12)]
    # except:
    years = [ time.item().year for time in times]

    ds = xr.open_mfdataset(files, combine = 'nested', concat_dim = 'time', decode_times = False).transpose('time', ...) 

    if 'height' in ds.dims:
        ds = ds.drop('height')
    if 'depth'  in ds.dims:
        ds = ds.drop('depth')

    if ensemble_id is not None:
        ds = ds.sel(ensembles = ensemble_id)

    if ensemble_mean:
        ds = ds.mean('ensembles')

    if rename_dict is not None:
        ds = ds.rename(**rename_dict)

    ls = [month_extractor(ds, i).expand_dims('year',axis = 0) for i in range(0,len(ds.time),12)]
    ds = xr.concat([ item.assign_coords(time = np.arange(1, len(item.time) + 1)) for item in ls], dim = 'year').assign_coords(year = np.unique(years))

    return ds

def month_extractor(ds, ind):
    return ds.isel( time = np.arange(ind, ind+12))



def load_csv_data(
    files: list[str | Path],
    *,
    axis: int = 0,
    sort_by: str = "year",
    rename_dict : dict = None,
    ignore_index: bool = True,
    **read_csv_kwargs,
) -> pd.DataFrame:
    """
    Load multiple CSV files, combine them, and sort by year.

    Parameters
    ----------
    files:
        List of CSV file paths.

    axis:
        Axis to concatenate along.
        Use axis=0 to stack rows, which is the usual case.

    sort_by:
        Column name to sort by after concatenation.

    ignore_index:
        Whether to reset the row index after concatenation.

    **read_csv_kwargs:
        Extra keyword arguments passed to pd.read_csv.
    """

    files = [Path(file) for file in files]

    if not files:
        raise ValueError("No CSV files were provided.")

    dfs = [pd.read_csv(file, **read_csv_kwargs) for file in files]

    df = pd.concat(
        dfs,
        axis=axis,
        ignore_index=ignore_index,
    )

    
    if rename_dict is not None:
        df.rename(columns=rename_dict, inplace=True)
    
    if sort_by not in df.columns:
        raise ValueError(f"Column {sort_by!r} not found in combined CSV data.")


    df.sort_values(sort_by).reset_index(drop=True)

    df.loc[df['lon'] > 180, 'lon'] = df.loc[df['lon'] > 180, 'lon'] - 360

    return df

# def load_ensemble(dir_in,
# # def load_forecasts(dir_in,
#                    var,
#                    y0=1981,
#                    y1=2022,
#                    verbose=True,
#                    ensemble_mean = False):
#     """
#      read in hindcast/forecast data
#     """
#     try:
#         data = []
#         if verbose:
#             print('loading forecasts..')
#         datasets = []
#         for f in sorted(Path(dir_in).glob(f"*.nc")):
#             ds = xr.open_dataset(f)[var].mean('ensembles').load() if ensemble_mean else xr.open_dataset(f)[var].load()
#             fname_split = f.stem.split("_")   # ripf info
#             iyr = int(fname_split[-6])
#             ds = ds.assign_coords(year=iyr+1)
#             ds = ds.rename({'lead_time' : 'time'})
#             ds['time'] = np.arange(ds.time.size)
#             ds = ds.expand_dims('year',axis=0)
#             datasets.append(ds)
#         ds_combined = xr.concat(datasets, dim='year').sortby('year')
#         if verbose:
#             print("done")
#     except IOError:
#         print('The input forecast file is missing! --check proper location')
#         quit()
#     return ds_combined






from pathlib import Path
import pandas as pd


def open_mfcsv(
    files: list[str | Path],
    *,
    concat_axis: int = 0,
    ignore_index: bool = True,
    add_source_column: bool = False,
    **read_csv_kwargs,
) -> pd.DataFrame:
    """
    Read multiple CSV files and concatenate them into one DataFrame.

    Similar idea to xr.open_mfdataset, but for CSV files.

    Parameters
    ----------
    files:
        List of CSV file paths.

    concat_axis:
        Axis to concatenate along. Usually 0 for stacking rows.

    ignore_index:
        Whether to reset the row index after concatenation.

    add_source_column:
        If True, adds a column showing which file each row came from.

    **read_csv_kwargs:
        Extra keyword arguments passed to pd.read_csv.
    """

    files = [Path(f) for f in files]

    if not files:
        raise ValueError("No CSV files were provided.")

    dataframes = []

    for file in files:
        df = pd.read_csv(file, **read_csv_kwargs)

        if add_source_column:
            df["source_file"] = str(file)

        dataframes.append(df)

    return pd.concat(
        dataframes,
        axis=concat_axis,
        ignore_index=ignore_index,
    )



