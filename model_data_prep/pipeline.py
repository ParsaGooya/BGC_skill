from __future__ import annotations

from .assimilation import process_assimilation, process_assimilation_extensions
from .canoe import process_canoe_bgc_assimilation
from .config import DataPrepConfig
from .hindcast import process_hindcast_and_forecast
from .historical import process_historical


def run_data_prep(cfg: DataPrepConfig) -> None:
    """Run selected data-preparation workflows."""
    for index, var in enumerate(cfg.var_list):
        realm = cfg.realm_for(index)

        print("=========================================================")
        print(f"Extracting data for model={cfg.model}, variable={var}, realm={realm}, lev_range={cfg.lev_range}")

        if cfg.canoe_assimilation_bgc:
            process_canoe_bgc_assimilation(cfg, var, realm)
            continue

        if cfg.hindcast:
            process_hindcast_and_forecast(cfg, var, realm)

        if cfg.assimilation:
            process_assimilation(cfg, var, realm)

        if cfg.assimilation_extracted_from_disc:
            process_assimilation_extensions(cfg, var, realm)

        if cfg.historical:
            process_historical(cfg, var, realm)
