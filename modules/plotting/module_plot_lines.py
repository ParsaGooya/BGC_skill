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
    Prepare a biome-averaged time series for temporal plotting and comparison.

    The function converts an input field into a one-dimensional temporal series
    suitable for downstream plotting or correlation analysis. For data containing
    a ``month`` dimension, the series can either be reduced to seasonal means or
    retained at monthly resolution. If a depth dimension is present, an optional
    depth range is selected before averaging over depth.

    For monthly-resolution output, the ``year`` and ``month`` dimensions are
    stacked into a continuous ``yearmonth`` coordinate expressed as fractional
    years. Otherwise, the function returns a yearly time coordinate shifted by
    0.5 years.

    Parameters
    ----------
    ds : xr.DataArray
        Input DataArray containing the biome-averaged data to prepare. The array
        may contain ``year``, ``month``, and ``lev`` dimensions depending on the
        temporal and vertical resolution of the source data.
    season : str
        Season used when reducing monthly data to seasonal means. The value is
        passed to ``seasonal_mean`` when ``monthly_res=False``. When
        ``monthly_res=True``, this must be ``"ANN"`` because the function retains
        the individual monthly values rather than calculating a seasonal mean.
    ldyr : int, optional
        Lead year used to determine the months associated with the requested
        season. For seasonal means, it is passed to ``seasonal_mean`` as
        ``ldyr_ini=ldyr`` and ``ldyr_end=ldyr + 1``. For monthly-resolution
        output, the corresponding annual sequence of month indices is obtained
        from ``get_season_indices``.
    lev_range : float or tuple[float, float] or None, optional
        Depth range to retain before averaging over the ``lev`` dimension. If
        provided, the range is resolved using ``resolve_depth_range``. If
        ``None``, the full available depth dimension is averaged. Ignored when
        the input does not contain a ``lev`` dimension.
    monthly_res : bool, optional
        If ``False``, monthly data are reduced to seasonal means using
        ``seasonal_mean`` and the returned temporal dimension is ``"year"``.
        If ``True``, monthly values are retained, ``season`` must be ``"ANN"``,
        and the ``year`` and ``month`` dimensions are stacked into a
        ``"yearmonth"`` dimension with fractional-year coordinates.

    Returns
    -------
    ts : xr.DataArray
        Prepared time series. Any ``lev`` dimension is removed by averaging over
        depth after optional depth-range selection. For monthly-resolution output,
        the temporal dimension is ``yearmonth``. Otherwise, the temporal
        dimension is expected to be ``year``.
    dim : str
        Name of the temporal dimension to use for downstream plotting or
        correlation calculations. Returns ``"yearmonth"`` when
        ``monthly_res=True`` and ``"year"`` otherwise.
    selected_range : float, tuple, or None
        Depth range selected by ``resolve_depth_range`` when ``lev_range`` is
        provided and the data contain a ``lev`` dimension. Returns ``None`` when
        no explicit depth range is selected.

    Raises
    ------
    ValueError
        If ``monthly_res=True`` while ``season`` is not ``"ANN"``.

    Notes
    -----
    If the input contains a ``month`` dimension and ``monthly_res=False``,
    ``seasonal_mean`` is used to calculate the requested seasonal composite for
    the specified lead year.

    When ``monthly_res=True``, the function retains individual monthly values.
    The relevant month indices are obtained using ``get_season_indices`` with
    ``season="ANN"``, ``ldyr_ini=ldyr``, and ``ldyr_end=ldyr + 1``.

    If a ``lev`` dimension is present and ``lev_range`` is specified,
    ``resolve_depth_range`` is used to restrict the field and return the actual
    selected range. Regardless of whether a range is explicitly specified, the
    remaining depth dimension is subsequently reduced using an arithmetic mean.

    For monthly-resolution output, ``year`` and ``month`` are stacked into a
    single ``yearmonth`` dimension. The coordinate is represented in fractional
    years using the month midpoint, so successive monthly values are separated
    by approximately ``1 / 12`` year.

    For non-monthly output, when a ``year`` dimension is present, 0.5 is added to
    the year coordinate so that each value is positioned at the midpoint of the
    corresponding year.

    Documentation produced with the assistance of AI.
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

    """
    Plot biome-averaged climatologies as a function of lead position.

    The function creates one row per requested biome and plots the selected
    datasets as one-dimensional climatological series along the x-axis. If the
    input data contain a ``year`` dimension, values are first averaged over year.
    Zero-valued entries are then removed to preserve the behavior of the original
    implementation. An optional depth range can be selected before averaging over
    the remaining ``lev`` dimension.

    When ``mask_biomes`` is provided, each biome row also includes a second panel
    showing the corresponding spatial biome mask.

    Parameters
    ----------
    ds_list : list[Exp]
        Dataset or experiment names to plot. Each entry must correspond to a key
        in ``ds_dict[biome]`` for every biome listed in ``biomes_to_plot``.
    ds_dict : dict[Biome, dict[Exp, state_dict]]
        Nested mapping from biome names to dataset or experiment names and their
        associated state objects. The ``data`` attribute provides the series to
        plot, while ``color`` and ``linestyle`` define the visual style of each
        dataset.
    biomes_to_plot : list[Biome]
        Biomes to include in the figure. A separate row of subplots is created for
        each biome in the order provided.
    mask_biomes : dict[Biome, xr.DataArray] or None, optional
        Optional mapping from biome names to spatial masks. When provided, a
        second subplot is created in each biome row and the corresponding mask is
        displayed with ``pcolormesh``.
    depth_range : float or Sequence[float] or None, optional
        Depth range passed to ``resolve_depth_range`` before vertical averaging.
        If the resulting series contains a ``lev`` dimension, that dimension is
        averaged after the requested range has been selected.
    ylim_min : float or None, optional
        Lower y-axis limit for each time-series panel.
    ylim_max : float or None, optional
        Upper y-axis limit for each time-series panel.
    xlim_min : float, optional
        Lower x-axis limit. The default is ``-1``.
    xlim_max : float or None, optional
        Upper x-axis limit. This value is also used when constructing x-axis tick
        locations on the final biome row.
    xticks_step : float, optional
        Spacing between generated x-axis tick locations.
    xticks_labels : Sequence[str] or None, optional
        Optional labels to assign to the generated x-axis tick locations on the
        final biome row. If provided, its length must match the number of
        generated ticks.
    var_name : str, optional
        Variable name included in each subplot title.
    xlabel : str, optional
        Label applied to the x-axis of the final biome row.
    ylabel : str, optional
        Label applied to the y-axis of each time-series panel.
    ncol_labels : int, optional
        Number of columns used for the legend.
    bbox : tuple, optional
        Bounding-box anchor passed to ``Axes.legend``.
    figsize : tuple[float, float], optional
        Figure size passed to ``matplotlib.pyplot.subplots``.
    fontsize : int, optional
        Base font size used for axis labels, tick labels, and legends. Subplot
        titles use ``fontsize + 2``.
    show_leg : bool, optional
        If ``True``, display a legend on each biome time-series panel.
    dir_name : str or Path or None, optional
        Directory in which to save the figure when ``save=True``. The directory
        is created if it does not already exist.
    file_name : str or None, optional
        Output filename without the ``.png`` extension when ``save=True``.
    return_fig_handles : bool, optional
        If ``True``, return the Matplotlib figure and axes objects. If ``False``,
        the function returns nothing.
    save : bool, optional
        If ``True``, save the generated figure as a PNG file.

    Returns
    -------
    tuple[matplotlib.figure.Figure, numpy.ndarray] or None
        If ``return_fig_handles=True``, returns ``(fig, axes)``, where ``fig`` is
        the Matplotlib figure and ``axes`` is the two-dimensional array of subplot
        axes returned by ``plt.subplots``. Otherwise, returns ``None``.

    Raises
    ------
    ValueError
        If ``xticks_labels`` is provided and its length does not match the number
        of generated x-axis tick locations.

    Notes
    -----
    For each dataset, the underlying ``state_dict.data`` is used as the source
    series. If a ``year`` dimension is present, the function first computes the
    mean over year so that the plotted curve represents a climatological lead
    profile.

    Zero-valued entries are removed before depth processing using
    ``ts[np.asarray(ts) != 0.0]``. This intentionally preserves the behavior of
    the earlier implementation and means that exact zeros are treated as values
    to exclude rather than valid climatological values.

    Depth selection is handled through ``resolve_depth_range``. If a ``lev``
    dimension remains after selection, the function calculates the arithmetic
    mean over depth before plotting.

    The x-coordinate is not taken directly from a DataArray coordinate. Instead,
    it is generated as ``np.arange(ts.size)``, so the horizontal axis represents
    the sequential position of the remaining values in the processed series.

    Each dataset is plotted using the ``color`` and ``linestyle`` metadata stored
    in its corresponding state object.

    Only the final biome row displays x-axis ticks and the x-axis label. Earlier
    rows have their x ticks removed to reduce visual clutter.

    When ``mask_biomes`` is supplied, the subplot layout contains two columns.
    The first column contains the climatological series and the second displays
    the corresponding biome mask using its ``lon`` and ``lat`` coordinates. The
    mask panel does not otherwise alter the time-series calculation.

    The subplot title is constructed as ``"<var_name> climatology - <biome>"``.
    If a depth range was selected, the supplied ``depth_range`` value is appended
    to the title.

    Documentation produced with the assistance of AI.
    """

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
    ELNINO_years: np.ndarray | None = None,
    LANINA_years: np.ndarray | None = None,
    return_fig_handles=False,
    save=False,
):

    """
    Plot biome-averaged time series for multiple datasets against a common target
    or reference series.

    The function creates one time-series panel for each requested biome and
    optionally adds a second panel showing the corresponding biome mask. Each
    dataset is prepared with ``prepare_biome_timeseries`` so that seasonal or
    monthly-resolution series, optional depth averaging, and lead-year selection
    are handled consistently.

    An optional reference series can be used to align the plotted datasets and,
    when requested, calculate both raw and detrended correlations. The plotted
    time series can additionally be smoothed with either a centered rolling mean
    or triangular smoothing, and linear trends can be overlaid. El Niño and
    La Niña years may also be marked with vertical reference lines.

    Parameters
    ----------
    ds_list : list[Exp]
        Dataset or experiment names to plot. Each entry must correspond to a key
        in ``ds_dicts[biome]`` for every biome listed in ``biomes_to_plot``.
    ds_dicts : dict[Biome, dict[Exp, state_dict]]
        Nested mapping from biome names to dataset or experiment names and their
        associated state objects. The ``data`` attribute provides the time series,
        while ``linestyle`` and ``color`` define the plotting style for each
        dataset.
    biomes_to_plot : list[Biome]
        Biomes to include in the figure. A separate row is created for each biome
        in the order provided.
    mask_biomes : dict[Biome, xr.DataArray] or None, optional
        Optional mapping from biome names to spatial masks. When provided, a
        second subplot is created in each biome row and the corresponding mask is
        displayed with ``pcolormesh``.
    ldyr : int, optional
        Lead year passed to ``prepare_biome_timeseries`` for seasonal or monthly
        time-series preparation.
    ref_ds : Exp or xr.DataArray, optional
        Reference dataset used for alignment and optional correlation
        calculations. If provided as a string, it is interpreted as a dataset key
        in the current biome's dictionary. If provided as an ``xr.DataArray``,
        the array must contain both ``year`` and ``month`` as dimensions and
        coordinates.
    title : str, optional
        Base text prepended to each biome subplot title.
    bbox : tuple, optional
        Bounding-box anchor passed to ``Axes.legend``.
    figsize : tuple[float, float], optional
        Figure size passed to ``matplotlib.pyplot.subplots``.
    wspace : float, optional
        Horizontal spacing between subplot columns passed to
        ``plt.subplots_adjust``.
    hspace : float, optional
        Vertical spacing between biome rows passed to
        ``plt.subplots_adjust``.
    dir_name : str or Path or None, optional
        Directory in which to save the figure when ``save=True``. The directory
        is created if it does not already exist.
    file_name : str or None, optional
        Output filename without the ``.png`` extension when ``save=True``.
    ylabel : str or None, optional
        Label applied to the y-axis of each time-series panel.
    season : str, optional
        Season passed to ``prepare_biome_timeseries``. When
        ``monthly_res=True``, this must be ``"ANN"``.
    lev_range : float, tuple[float, float], or None, optional
        Depth range passed to ``prepare_biome_timeseries`` before averaging over
        depth. If no explicit range is selected, the helper determines the
        corresponding output behavior.
    monthly_res : bool, optional
        If ``False``, plot seasonal or annual biome-averaged series along the
        ``year`` dimension. If ``True``, retain monthly resolution and plot along
        the fractional-year ``yearmonth`` coordinate produced by
        ``prepare_biome_timeseries``.
    correlations : bool, optional
        If ``True``, calculate the correlation of each dataset with the aligned
        reference series and include both the raw and detrended correlations in
        the legend label. Requires ``ref_ds`` to be provided.
    rolling : int or None, optional
        Window length for centered rolling-mean smoothing. Applied independently
        to the reference series and each dataset along their prepared temporal
        dimension. Cannot be used together with ``triangular_smoothing``.
    triangular_smoothing : int or None, optional
        Smoothing parameter passed to ``trismooth``. Applied independently to the
        reference and dataset values. Cannot be used together with ``rolling``.
    show_trend : bool, optional
        If ``True``, calculate a trend for each aligned dataset using ``trend``
        and overlay it as a dashed line in the same dataset color.
    ELNINO_years : np.ndarray or None, optional
        El Niño years to mark with vertical dotted red lines. Only years falling
        within the plotted temporal range are shown.
    LANINA_years : np.ndarray or None, optional
        La Niña years to mark with vertical dotted blue lines. Only years falling
        within the plotted temporal range are shown.
    return_fig_handles : bool, optional
        If ``True``, return the Matplotlib figure and axes objects. If ``False``,
        the function returns nothing.
    save : bool, optional
        If ``True``, save the generated figure as a PNG file.

    Returns
    -------
    tuple[matplotlib.figure.Figure, numpy.ndarray] or None
        If ``return_fig_handles=True``, returns ``(fig, axes)``, where ``fig`` is
        the Matplotlib figure and ``axes`` contains the subplot axes. Otherwise,
        returns ``None``.

    Raises
    ------
    ValueError
        If both ``rolling`` and ``triangular_smoothing`` are specified.
    ValueError
        If ``monthly_res=True`` while ``season`` is not ``"ANN"``.
    ValueError
        If ``correlations=True`` and ``ref_ds`` is ``None``.
    ValueError
        If ``ref_ds`` is supplied as an ``xr.DataArray`` without both ``year``
        and ``month`` dimensions and coordinates.

    Notes
    -----
    For each biome, the optional reference field is prepared first using
    ``prepare_biome_timeseries`` with the same ``season``, ``ldyr``,
    ``lev_range``, and ``monthly_res`` settings as the target datasets.

    If ``ref_ds`` is a string, the reference field is retrieved from the current
    biome's dataset dictionary. If it is supplied directly as an
    ``xr.DataArray``, the same reference array is used for each biome after
    validation of its ``year`` and ``month`` dimensions and coordinates.

    Rolling and triangular smoothing are mutually exclusive. A rolling mean is
    calculated with ``center=True`` along the prepared temporal dimension.
    Triangular smoothing is applied by copying the DataArray and replacing its
    values with the output from ``trismooth``.

    Each target dataset is aligned with the prepared reference using
    ``xr.align(..., join="inner")``. Consequently, plotting and correlation
    calculations use only temporal coordinates common to both series. If no
    reference is supplied, the prepared target series is plotted without this
    alignment step.

    When ``correlations=True``, the legend label contains two correlation
    coefficients. The first is the direct xarray correlation between the aligned
    dataset and reference. The value in parentheses is calculated after both
    series have been detrended using ``trend(..., return_detrended=True)``.

    When ``show_trend=True``, the trend returned by ``trend`` is plotted as a
    dashed line in the same color as the corresponding dataset.

    Only the final biome row retains x-axis ticks and labels. Earlier rows have
    their x ticks removed to reduce visual clutter.

    If an explicit depth range is selected by ``prepare_biome_timeseries``, the
    resolved range is appended to the biome subplot title. Information about
    rolling or triangular smoothing is also appended when applicable.

    El Niño and La Niña markers are filtered against the minimum and maximum
    values of the plotted temporal coordinate before vertical lines are added.
    El Niño events are shown in red and La Niña events in blue, both using dotted
    lines with partial transparency.

    When ``mask_biomes`` is provided, each biome row contains a second panel
    showing the corresponding spatial mask using its ``lon`` and ``lat``
    coordinates with display limits of 0 and 1.

    A legend is always added to each biome time-series panel using the supplied
    ``bbox`` anchor and without a frame.

    Documentation produced with the assistance of AI.
    """
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