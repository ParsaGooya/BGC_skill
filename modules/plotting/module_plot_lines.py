import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import xarray as xr


from modules.analysis.module_global_averages import area_weighted_avg
from modules.analysis.module_data_postprocessing import (trismooth,
                                                         trend)
from modules.data_info.module_state_dict import state_dict
from modules.plotting.utils import *



def prepare_biome_timeseries(
    ds: xr.DataArray,
    season: str,
    ldyr: int = 0,
    lev_range: float | tuple[float, float] | None = None,
    monthly_res: bool = False,
) -> tuple[xr.DataArray, str, float | tuple]:
    """
    Prepare a biome-averaged time series for plotting.

    Returns
    -------
    ts : xr.DataArray
        Prepared time series.
    dim : str
        Temporal dimension used for plotting/correlation.
    """
    if "month" in ds.dims:
        if monthly_res:
            if season != "ANN":
                raise ValueError(
                    "monthly_res=True requires season='ANN'."
                )

            month_indices = get_season_indices(
                season="ANN",
                ldyr_ini=ldyr,
                ldyr_end=ldyr + 1,
            )

            ts = ds.isel(month=month_indices)

        else:
            ts = seasonal_mean(
                ds,
                season=season,
                ldyr_ini=ldyr,
                ldyr_end=ldyr + 1,
            )
    else:
        ts = ds

    selected_range = None
    if "lev" in ts.dims:
        if lev_range is not None:
            ts, selected_range = resolve_depth_range(ts, lev_range) 

        ts = ts.mean("lev")

    if monthly_res:
        ts = ts.stack(
            yearmonth=("year", "month"),
        )

        yearmonth = (
            ts.year.values
            + (np.mod(ts.month.values - 0.5, 12) ) / 12
        )

        ts = ts.assign_coords(
            yearmonth=yearmonth,
        )

        return ts, "yearmonth", selected_range

    if "year" in ts.dims:
        ts = ts.assign_coords(
            year=ts.year.values + 0.5,
        )

    return ts, "year", selected_range


                                                                             
def plot_ts_vs_lead_biomes(
    ds_list: list[Exp],
    ds_dict: dict[Biome, dict[Exp, state_dict]],
    biomes_to_plot: list[Biome],
    mask_biomes: dict[Biome, xr.DataArray]=None,
    depth_range: float | Sequence[float] = None,
    ylim_min=None,
    ylim_max=None,
    xlim_min=-1,
    xlim_max=None,
    xticks_step=1,
    xticks_labels=None,
    var_name="",
    xlabel="",
    ylabel="",
    ncol_labels=1,
    bbox=(0.68, 0.5, 0.5, 0.5),
    figsize=(10, 30),
    fontsize=10,
    show_leg=True,
    dir_name=None,
    file_name=None,
    return_fig_handles=False,
    save=False,
):


    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    if mask_biomes is not None:
        fig, axes = plt.subplots(
            len(biomes_to_plot),
            2,
            figsize=figsize,
            squeeze=False,
            gridspec_kw={"width_ratios": [2, 0.75]},
        )
    else:
        fig, axes = plt.subplots(
            len(biomes_to_plot),
            1,
            figsize=figsize,
            squeeze=False,
        )

    # ------------------------------------------------------------------
    # Biomes
    # ------------------------------------------------------------------
    for biome_idx, biome in enumerate(biomes_to_plot):

        ax = axes[biome_idx, 0]

        # --------------------------------------------------------------
        # Time series
        # --------------------------------------------------------------
        for series_idx, name in enumerate(ds_list):

            ts = ds_dict[biome][name].data
            if "year" in ts.dims:
                ts = ts.mean("year")

            # Preserve the behavior of the old function:
            # remove zero-valued entries.
            ts = ts[np.asarray(ts) != 0.0]
            ts, selected_depth = resolve_depth_range(ts, depth_range)
            if "lev" in ts.dims:
                ts = ts.mean("lev")
            xx = np.arange(ts.size)

            color = ds_dict[biome][name].color
            linestyle = ds_dict[biome][name].linestyle

            ax.plot(
                xx,
                ts,
                linestyle,
                markersize=5,
                color=color,
                label=name,
            )

        # --------------------------------------------------------------
        # Axis formatting
        # --------------------------------------------------------------
        ax.set_xlim(
            xlim_min,
            xlim_max,
        )

        ax.set_ylim(
            ylim_min,
            ylim_max,
        )

        title = f"{var_name} climatology - {biome}"
        if selected_depth is not None:
            title += f" {depth_range}"

        ax.set_title(
            title,
            fontsize=fontsize + 2,
        )

        ax.set_ylabel(
            ylabel,
            fontsize=fontsize,
        )

        ax.tick_params(
            axis="both",
            labelsize=fontsize,
        )

        # Only show x-axis labels on the final row.
        if biome_idx < len(biomes_to_plot) - 1:
            ax.set_xlabel("")
            ax.set_xticks([])

        else:
            ax.set_xlabel(
                xlabel,
                fontsize=fontsize,
            )

            xticks = np.arange(
                xlim_min,
                xlim_max + 1 if xlim_max is not None else xlim_max ,
                xticks_step,
            )

            ax.set_xticks(xticks)

            if xticks_labels is not None:
                if len(xticks_labels) != len(xticks):
                    raise ValueError(
                        "'xticks_labels' must have the same length "
                        "as the generated x ticks."
                    )

                ax.set_xticklabels(
                    xticks_labels,
                    fontsize=fontsize,
                )

        # --------------------------------------------------------------
        # Biome mask
        # --------------------------------------------------------------
        if mask_biomes is not None:
            biome_mask = mask_biomes[biome]
            mask_ax = axes[biome_idx, 1]

            mask_ax.pcolormesh(
                biome_mask.lon,
                biome_mask.lat,
                biome_mask,
            )

        # --------------------------------------------------------------
        # Legend
        # --------------------------------------------------------------
        if show_leg:
            ax.legend(
                loc="best",
                bbox_to_anchor=bbox,
                fontsize=fontsize,
                handlelength=1,
                ncol=ncol_labels,
                frameon=False,
            )

    plt.subplots_adjust(
        wspace=0.55,
        hspace=0.15,
    )

    if save:
        Path(dir_name).mkdir(
            parents=True,
            exist_ok=True,
        )

        plt.savefig(
            f"{dir_name}/{file_name}.png",
            bbox_inches="tight",
            dpi=300,
        )

    if return_fig_handles:
        return fig, axes

            
def plot_ts_biomeavg_on_target(
    ds_list: list[Exp],
    ds_dicts: dict[Biome, dict[Exp, state_dict]],
    biomes_to_plot: list[Biome],
    mask_biomes: dict[Biome, xr.DataArray] | None = None,
    ldyr=0,
    ref_ds: Exp | xr.DataArray = "obs",
    title="",
    bbox=(0.68, 0.5, 0.5, 0.5),
    figsize=(10, 45),
    wspace=0.35,
    hspace=0.35,
    dir_name=None,
    file_name=None,
    ylabel=None,
    season="ANN",
    lev_range=None,
    monthly_res=False,
    correlations=False,
    rolling=None,
    triangular_smoothing=None,
    show_trend=False,
    ELNINO_years: np.ndarray[float] | None = None,
    LANINA_years: np.ndarray[float] | None = None,
    return_fig_handles=False,
    save=False,
):
    if rolling is not None and triangular_smoothing is not None:
        raise ValueError(
            "Specify either rolling or triangular_smoothing, not both."
        )

    if monthly_res and season != "ANN":
        raise ValueError(
            "monthly_res=True currently requires season='ANN'."
        )

    if correlations and ref_ds is None:
        raise ValueError(
            "A reference dataset is required when correlations=True."
        )

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    if mask_biomes is not None:
        fig, axes = plt.subplots(
            len(biomes_to_plot),
            2,
            figsize=figsize,
            gridspec_kw={"width_ratios": [2, 0.75]},
        )
    else:
        fig, axes = plt.subplots(
            len(biomes_to_plot),
            1,
            figsize=figsize,
            squeeze=False,
        )

    # ------------------------------------------------------------------
    # Biomes
    # ------------------------------------------------------------------
    for biome_idx, biome in enumerate(biomes_to_plot):

        if mask_biomes is not None:
            ax = axes[biome_idx, 0]
            mask_ax = axes[biome_idx, 1]
        else:
            ax = axes[biome_idx, 0]

        ds_dict = ds_dicts[biome]

        # --------------------------------------------------------------
        # Reference
        # --------------------------------------------------------------
        ref = None

        if ref_ds is not None:
            if isinstance(ref_ds, str):
                ref_data = ds_dict[ref_ds].data
            else:
                ref_data = ref_ds
                if ("year" not in ref_data.dims or
                    "year" not in ref_data.coords or
                    "month" not in ref_data.dims or
                    "month" not in ref_data.coords):
                    raise ValueError(
                        "The provided ref dataset must have year and month dimensions and coords."
                    )


            ref, ref_dim, _ = prepare_biome_timeseries(
                ref_data,
                season=season,
                ldyr=ldyr,
                lev_range=lev_range,
                monthly_res=monthly_res,
            )

            if rolling is not None:
                ref = ref.rolling(
                    {ref_dim: rolling},
                    center=True,
                ).mean()

            elif triangular_smoothing is not None:
                ref = ref.copy()
                ref[:] = trismooth(
                    ref.values,
                    triangular_smoothing,
                )

        # --------------------------------------------------------------
        # Datasets
        # --------------------------------------------------------------
        for name in ds_list:

            ts, dim, selected_depth = prepare_biome_timeseries(
                ds_dict[name].data,
                season=season,
                ldyr=ldyr,
                lev_range=lev_range,
                monthly_res=monthly_res,
            )

            if rolling is not None:
                ts = ts.rolling(
                    {dim: rolling},
                    center=True,
                ).mean()

            elif triangular_smoothing is not None:
                ts = ts.copy()
                ts[:] = trismooth(
                    ts.values,
                    triangular_smoothing,
                )

            # ----------------------------------------------------------
            # Align reference / target
            # ----------------------------------------------------------
            if ref is not None:
                ts_aligned, ref_aligned = xr.align(
                    ts,
                    ref,
                    join="inner",
                )
            else:
                ts_aligned = ts
                ref_aligned = None

            xx = ts_aligned[dim].values

            # ----------------------------------------------------------
            # Legend label / correlations
            # ----------------------------------------------------------
            if correlations:
                corr = xr.corr(
                    ts_aligned,
                    ref_aligned,
                    dim=dim,
                ).values

                corr_detrended = xr.corr(
                    trend(
                        ts_aligned,
                        dim=dim,
                        return_detrended=True,
                    )[1],
                    trend(
                        ref_aligned,
                        dim=dim,
                        return_detrended=True,
                    )[1],
                    dim=dim,
                ).values

                plot_label = (
                    f"{name} "
                    f"{np.round(corr, 2)} "
                    f"({np.round(corr_detrended, 2)})"
                )

            else:
                plot_label = name

            # ----------------------------------------------------------
            # Plot time series
            # ----------------------------------------------------------
            ax.plot(
                xx,
                ts_aligned,
                ds_dict[name].linestyle,
                label=plot_label,
                color=ds_dict[name].color,
            )

            if show_trend:
                ts_trend = trend(
                    ts_aligned,
                    dim=dim,
                )

                ax.plot(
                    xx,
                    ts_trend,
                    linestyle="dashed",
                    color=ds_dict[name].color,
                )

        # --------------------------------------------------------------
        # Axis formatting
        # --------------------------------------------------------------
        if biome_idx < len(biomes_to_plot) - 1:
            ax.set_xlabel("")
            ax.set_xticks([])

        biome_title = f"{title} - {biome}"

        if selected_depth is not None:
                biome_title += f' {selected_depth}'

        if rolling is not None:
            biome_title += f" rolling mean {rolling}"

        elif triangular_smoothing is not None:
            biome_title += (
                f" triangular smoothing {triangular_smoothing}"
            ) 
                    

        ax.set_title(biome_title)
        ax.set_ylabel(ylabel)

        # --------------------------------------------------------------
        # ENSO markers
        # --------------------------------------------------------------
        if ELNINO_years is not None:

            for year in [year for year in ELNINO_years
                 if (xx.min() <= year and year  <= xx.max())]:
            
                ax.axvline(
                    x=year,
                    linestyle="dotted",
                    color="r",
                    alpha=0.25,
                )

        if LANINA_years is not None:
            for year in [year for year in LANINA_years
                 if (xx.min() <= year and year  <= xx.max())]:
                
                ax.axvline(
                    x=year,
                    linestyle="dotted",
                    color="b",
                    alpha=0.25,
                )

        # --------------------------------------------------------------
        # Biome mask
        # --------------------------------------------------------------
        if mask_biomes is not None:
            biome_mask = mask_biomes[biome]

            mask_ax.pcolormesh(
                biome_mask.lon,
                biome_mask.lat,
                biome_mask,
                vmin=0,
                vmax=1,
            )

        ax.legend(
            loc="best",
            bbox_to_anchor=bbox,
            handlelength=2,
            frameon=False,
        )

    plt.subplots_adjust(
        wspace=wspace,
        hspace=hspace,
    )

    if save:
        Path(dir_name).mkdir(
            parents=True,
            exist_ok=True,
        )

        plt.savefig(
            f"{dir_name}/{file_name}.png",
            bbox_inches="tight",
            dpi=300,
        )

    if return_fig_handles:
        return fig, axes
        # return xx, ts_aligned