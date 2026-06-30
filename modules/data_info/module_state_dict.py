'''

'''
from pathlib import Path
import pandas as pd
import xarray as xr
from modules.analysis.module_data_load import load_nc_data, load_csv_data
from modules.analysis.module_data_postprocessing import nanmasker
import yaml

with open(Path("/home/rpg002/BGC_skill/configs") / "styles.yaml", "r") as f:
    EXPERIMENT_styles = yaml.safe_load(f)


exlude_dirs = [Path('/space/hall7/sitestore/eccc/crd/cccma/users/rpg002/data/no3/observation/monthly_data'),
               Path('/space/hall7/sitestore/eccc/crd/cccma/users/rpg002/data/tas/observations/raw_era5'),
              ]


def model_color(experiment: str, model_key: str) -> str:
    experiment = experiment.lower()

    if experiment == "observation":
        return "black"

    if experiment == "assimilation":
        if model_key == "CanESM5-CanOE":
            return "tab:blue"
        return "tab:green"

    if experiment == "historical":
        if model_key == "CanESM5-CanOE":
            return "tab:purple"
        return "tab:red"

    if experiment == "hindcast":
        return "tab:blue"

    if experiment == "control":
        return "grey"

    return "black"


def infer_realm(var: str) -> str:
    if var in {"uas", "vas", "tauu", "tauv"}:
        return "Amon"
    return "Omon"


def find_data_files(
    path: Path,
    var: str,
    realm: str | None = None,
) -> list[Path]:
    realm = realm or infer_realm(var)

    nc_files = sorted(path.glob(f"{var}_{realm}_*.nc"))
    if not nc_files:
        nc_files = sorted(path.glob(f"{var}_*.nc"))
    if not nc_files:
        nc_files = sorted(path.glob("*.nc"))

    csv_files = sorted(path.glob(f"*{var}*.csv"))
    if not csv_files:
        csv_files = sorted(path.glob("*.csv"))

    if nc_files and csv_files:
        raise ValueError(
            f"Found both NetCDF and CSV files in {path}. "
            "Expected only one file type per model/source directory."
        )

    return nc_files or csv_files



def infer_years_from_nc_files(files: list[Path]) -> tuple[int, int]:
    if not files:
        raise ValueError("No NetCDF files were provided.")

    years = []

    for file in files:
        with xr.open_dataset(file, decode_times=True) as ds:
            if "time" not in ds.coords:
                raise ValueError(f"No time coordinate found in {file}")

            time = ds["time"]
            try:
            
                years.append(int(time.dt.year.isel(time=0).values))
                years.append(int(time.dt.year.isel(time=-1).values))
            except:
                years.append(None)

    return min(years), max(years)


def infer_years_from_csv_files(files: list[Path]) -> tuple[int, int]:
    if not files:
        raise ValueError("No CSV files were provided.")

    years = []

    for file in files:
        df = pd.read_csv(file, usecols=["year"])
        years.append(int(df["year"].min()))
        years.append(int(df["year"].max()))

    return min(years), max(years)


def infer_years(files: list[Path]) -> tuple[int, int]:
    if not files:
        raise ValueError("No files were provided.")

    suffixes = {file.suffix for file in files}

    if suffixes == {".nc"}:
        return (*infer_years_from_nc_files(files), suffixes)

    if suffixes == {".csv"}:
        return (*infer_years_from_csv_files(files), suffixes)

    raise ValueError(
        f"Unsupported or mixed file types: {suffixes}. "
        "Expected only .nc files or only .csv files."
    )





def resolve_experiment_dir(
    data_directory: str | Path,
    var: str,
    experiment: str,
    assimilation_BGC_run_id: int | None = None,
) -> Path:
    data_directory = Path(data_directory)

    if assimilation_BGC_run_id is not None and experiment == "assimilation":
        return data_directory / var / f"assimilation_bgc{assimilation_BGC_run_id}"

    return data_directory / var / experiment



def normalize_model_key(model_name: str) -> str:
    if model_name.startswith("CanESM5-CanOE"):
        return "CanESM5-CanOE"

    return model_name

def include_model_dir(
    source_dir: Path,
    CanOE_assimilation_BGC_run_id: int | None = None,
    exlude_dirs : list[Path] = exlude_dirs) -> bool:

    if source_dir in exlude_dirs:
        return False

    name = source_dir.name

    if name.startswith("CanESM5-CanOE") and 'assimilation' in str(source_dir):
        if CanOE_assimilation_BGC_run_id is None:
            return True
  
        return name == f"CanESM5-CanOE_{CanOE_assimilation_BGC_run_id}"

    return True

import dataclasses
@dataclasses.dataclass
class state_dict:
    var: str
    files: list[Path]
    experiment: str
    model_key: str
    assimilation_BGC_run_id: int | None = None
    CanOE_assimilation_BGC_run_id: int | None = None
    

    def __post_init__(self):
        self.y0, self.y1, self.type = infer_years(self.files)

        self.dir= [str(file) for file in self.files]
        self.color =  model_color(self.experiment, normalize_model_key(self.model_key))  
        self.linestyle =  EXPERIMENT_styles.get(self.experiment).get("linestyle", None)
        self.alpha =  EXPERIMENT_styles.get(self.experiment).get("alpha", 1)
        self.marker =  EXPERIMENT_styles.get(self.experiment).get("marker", None)
        

        if self.assimilation_BGC_run_id is not None:
            self.assimilation_BGC_run_id = self.assimilation_BGC_run_id

        if self.model_key == "CanESM5-CanOE" and self.CanOE_assimilation_BGC_run_id is not None:
            self.CanOE_assimilation_BGC_run_id = self.CanOE_assimilation_BGC_run_id

    def PrintLoc(self):
        print('/'.join(self.dir[0].split('/')[:-1]))

    def load_data(self, varx = None, rename_dict : dict = None,  return_mask = False,unit_change : float = None, **kwargs):
        if unit_change is None:
            unit_change = 1
        mask = None
        if varx is None:
            varx = self.var

        if '.nc' in self.type:
                
            ds = load_nc_data(self.files, rename_dict =rename_dict, **kwargs)[varx].transpose(...,'lat','lon')
            _, mask = nanmasker(ds.stack(ref = ('year','time')), dim = 'ref', return_mask= True)
            self.data = ds.squeeze() * mask * unit_change


        elif '.csv' in self.type:

            ds = load_csv_data(self.files, rename_dict =rename_dict, **kwargs)
            ds.rename(columns={varx : 'obs'}, inplace=True)
            self.data = ds.squeeze() * unit_change

        if return_mask and mask is not None:
            return mask.rename('mask')
    
    def apply_nc_mask(self, mask : xr.DataArray | float):
        if '.nc' in self.type:
            self.data = self.data * mask

    def sel(self, time_selection_dict : dict):
        if '.nc' in self.type:
            self.data = self.data.sel(**time_selection_dict)
        
        elif '.csv' in self.type:
            for key, values in time_selection_dict.items():
                if isinstance(values, slice):
                    min = values.start
                    max = values.stop
                else:
                    min = values.min()
                    max = values.max()

                self.data.loc[(self.data[key] >= min) & (self.data[key] <= max)]


def get_data_dict(
    data_directory: str | Path,
    var: str,
    experiment: str,
    assimilation_BGC_run_id: int | None = None,
    CanOE_assimilation_BGC_run_id: int | None = None,
    realm: str | None = None,
) -> dict:
    """
    Discover available data for a variable and experiment.

    Directory structure assumed
    ---------------------------
    data_directory/
      var/
        experiment/
          model_or_source_1/
            file_1.nc
            file_2.nc

          model_or_source_2/
            file_1.csv
            file_2.csv

    Notes
    -----
    - Files must be inside model/source directories.
    - Each model/source directory may contain multiple files.
    - Each model/source directory must contain either NetCDF files or CSV files, not both.
    - y0/y1 are inferred across all files in that model/source directory.
    - state["dir"] is a list of file paths, suitable for xr.open_mfdataset.
    """

    experiment = experiment.lower()
    realm = realm or infer_realm(var)

    experiment_dir = resolve_experiment_dir(
        data_directory=data_directory,
        var=var,
        experiment=experiment,
        assimilation_BGC_run_id=assimilation_BGC_run_id,
    )

    
    if not experiment_dir.exists():
        print(f"{var} {experiment} does not exist")
        return None


    data_dict = {}

    source_dirs = sorted(p for p in experiment_dir.iterdir() if p.is_dir())

    for source_dir in source_dirs:
        if not include_model_dir(
            source_dir,
            CanOE_assimilation_BGC_run_id=CanOE_assimilation_BGC_run_id,
        ):
            continue

        files = find_data_files(source_dir, var=var, realm=realm)

        if not files:
            continue

        model_key = source_dir.name

        data_dict[model_key] = state_dict(
            var=var,
            files=files,
            experiment=experiment,
            model_key=model_key,
            assimilation_BGC_run_id=assimilation_BGC_run_id,
            CanOE_assimilation_BGC_run_id=CanOE_assimilation_BGC_run_id,
        )

    return data_dict
