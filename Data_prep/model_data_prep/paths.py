from __future__ import annotations

from pathlib import Path

from .config import DataPrepConfig


def output_dir(cfg: DataPrepConfig, var: str, collection: str, model: str | None = None) -> Path:
    model = model or cfg.model
    return cfg.user_data_root / var / collection / model


def dcpp_dir(cfg: DataPrepConfig, experiment: str, model: str | None = None) -> Path:
    model = model or cfg.model
    return cfg.cmip6_root / "DCPP" / "CCCma" / model / experiment


def historical_dir(cfg: DataPrepConfig, model: str | None = None) -> Path:
    model = model or cfg.model
    return cfg.cmip6_root / "CMIP" / "CCCma" / model / "historical"


def ssp245_dir(cfg: DataPrepConfig, model: str | None = None) -> Path:
    model = model or cfg.model
    return cfg.cmip6_root / "ScenarioMIP" / "CCCma" / model / "ssp245"


def variable_dir(base: Path, realization: str, realm: str, var: str, cfg: DataPrepConfig) -> Path:
    return base / realization / realm / var / cfg.grid / cfg.version


def initialized_variable_dir(base: Path, year: int, realization: str, realm: str, var: str, cfg: DataPrepConfig) -> Path:
    return base / f"s{year}-{realization}" / realm / var / cfg.grid / cfg.version
