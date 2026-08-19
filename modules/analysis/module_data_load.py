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

    times = ds = xr.open_mfdataset(files, combine = 'nested', concat_dim = 'time').time
    ds = xr.open_mfdataset(files, combine = 'nested', concat_dim = 'time', decode_times = False).rename({'time' : 'year'}).transpose('year', ...)   

    years = times.dt.year 
    months = times.dt.month
    ds = ds.assign_coords({
        "year" : years.values,
    })

    ls = []
    for year in np.unique(years):
        
        ds_year = ds.where(ds.year == year, drop = True)
        month = months[years == year]
        ls.append(ds_year.rename({'year' : 'month'}).assign_coords({'year' : year, 'month' : month.values}))
    
    ds = xr.concat(ls, dim = 'year').assign_coords(year = np.unique(years))

    if 'height' in ds.dims:
        ds = ds.drop('height')
    if 'depth'  in ds.dims:
        ds = ds.drop('depth')

    has_ensemble = 'ensembles' in ds.dims

    if has_ensemble and ensemble_id is not None:
        ds = ds.sel(ensembles = ensemble_id)

    if has_ensemble and ensemble_mean:
        ds = ds.mean('ensembles')

    if rename_dict is not None:
        ds = ds.rename(**rename_dict)


    return ds


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



