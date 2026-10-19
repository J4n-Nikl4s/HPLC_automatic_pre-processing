# enzymeml_plot.py
from __future__ import annotations

import math
from typing import Any, Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from pyenzyme import EnzymeMLDocument


def _nice_grid(n: int) -> Tuple[int, int]:
    if n <= 0:
        return 1, 1
    c = math.ceil(math.sqrt(n))
    r = math.ceil(n / c)
    return r, c


def _group_measurements(
    measurements: List[Dict[str, Any]],
) -> Dict[str, List[Dict[str, Any]]]:
    groups: Dict[str, List[Dict[str, Any]]] = {}
    for m in measurements:
        gid = m.group_id or "ungrouped"
        groups.setdefault(str(gid), []).append(m)
    return groups


def _collect_species_series(
    group: List[Dict[str, Any]],
) -> Dict[str, Dict[float, List[float]]]:
    # species_id -> { time_value: [replicate_values...] }
    out: Dict[str, Dict[float, List[float]]] = {}
    for meas in group:
        for sd in meas.species_data:
            sid = sd.species_id
            t = sd.time or []
            y = sd.data or []
            if not sid or not t or not y or len(t) != len(y):
                continue
            tmap = out.setdefault(str(sid), {})
            for ti, yi in zip(t, y):
                tmap.setdefault(float(ti), []).append(float(yi))
    return out


def _aggregate_group(group: List[Dict[str, Any]]):
    # Returns: species_id -> (times, means, stds, meta)
    series = _collect_species_series(group)
    meta: Dict[str, Dict[str, str]] = {}
    for meas in group:
        for sd in meas.species_data:
            sid = sd.species_id
            if not sid or sid in meta:
                continue
            tu = sd.time_unit.id or sd.time_unit.name or ""
            du = sd.data_unit.id or sd.data_unit.name or ""
            dt = sd.data_type or ""
            meta[str(sid)] = {"time_unit": tu, "data_unit": du, "data_type": dt}

    agg: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray, Dict[str, str]]] = {}
    for sid, tmap in series.items():
        if not tmap:
            continue
        times = np.array(sorted(tmap.keys()), dtype=float)
        means = np.array([np.mean(tmap[t]) for t in times], dtype=float)
        stds = np.array(
            [np.std(tmap[t], ddof=1) if len(tmap[t]) > 1 else 0.0 for t in times],
            dtype=float,
        )
        agg[sid] = (
            times,
            means,
            stds,
            meta.get(sid, {"time_unit": "", "data_unit": "", "data_type": ""}),
        )
    return agg


def plot_replicates(
    doc: EnzymeMLDocument,
    figsize_per_panel: Tuple[float, float] = (5.0, 4.0),
    sharex: bool = False,
    sharey: bool = False,
    tight: bool = False,
):
    """Takes an EnzymeML-like document (dict). Builds subplots by measurement.group_id and plots mean±std per species."""
    measurements = doc.measurements
    groups = _group_measurements(measurements)
    gids = list(groups.keys()) or ["ungrouped"]

    r, c = _nice_grid(len(gids))
    fig = plt.figure(figsize=(figsize_per_panel[0] * c, figsize_per_panel[1] * r))
    axes: List[plt.Axes] = []

    for i, gid in enumerate(gids, 1):
        ax = fig.add_subplot(
            r,
            c,
            i,
            sharex=axes[0] if (sharex and axes) else None,
            sharey=axes[0] if (sharey and axes) else None,
        )
        axes.append(ax)
        agg = _aggregate_group(groups[gid])

        # Labels from first species present
        xlab = ylab = ""
        for _, (_t, _m, _s, meta) in agg.items():
            tu, du, dt = (
                meta["time_unit"],
                meta["data_unit"],
                meta["data_type"],
            )
            xlab = xlab or (f"time [{tu}]" if tu else "time")
            ylab = ylab or (
                f"{dt.name.lower()} [{du}]" if dt and du else (dt or du or "value")
            )
            break

        for sid, (t, m, s, _meta) in sorted(agg.items()):
            if t.size == 0:
                continue
            ax.errorbar(t, m, yerr=s, fmt="o", capsize=3, label=str(sid))

        ax.set_title(str(gid))
        if xlab:
            ax.set_xlabel(xlab)
        if ylab:
            ax.set_ylabel(ylab)
        ax.legend(title="species_id", loc="best")
        ax.grid(True, linestyle=":", linewidth=0.5)

    if tight:
        fig.tight_layout()
    return fig, axes


def _species_names(doc: EnzymeMLDocument) -> Dict[str, str]:
    # species_id -> name, falls back to the id if no name is set
    names: Dict[str, str] = {}
    for species in [*doc.small_molecules, *doc.proteins, *doc.complexes]:
        names[str(species.id)] = str(getattr(species, "name", None) or species.id)
    return names


def plot_replicates_broken_axis(
    doc: EnzymeMLDocument,
    x_ranges: List[Tuple[float, float]] = [(0, 130), (1430, 1460)],
    color_dict: Dict[str, str] | None = None,
    exclude: Tuple[str, ...] = ("Neg. Kontrolle", "Summe"),
    figsize_per_panel: Tuple[float, float] = (10.0, 5.0),
    breakmark_size: float = 0.015,
    wspace: float = 0.05,
    xlabel: str | None = None,
    ylabel: str | None = None,
    xlabel_pad: float = 25,
):
    """Plots mean±std of the replicates per measurement.group_id with a broken x-axis.

    Uses the `brokenaxes` package: each entry of `x_ranges` becomes one x-axis segment,
    the segments share the y-axis. The last point of a segment is connected to the
    first point of the next one by a dashed line.

    Args:
        doc: EnzymeML document, replicates are grouped by measurement.group_id.
        x_ranges: x-limits of the axis segments, e.g. [(0, 130), (1430, 1460)].
        color_dict: colors per species name or species id, also sets the legend order.
            Species without a color get the default cycle.
        exclude: species whose name or id contains one of these strings are not plotted.
        figsize_per_panel: size of the figure per group.
        breakmark_size: size of the diagonal break marks.
        wspace: gap between the axis segments.
        xlabel: x-axis label, derived from the time unit if None.
        ylabel: y-axis label, derived from the data type and unit if None.
        xlabel_pad: distance of the x-axis label from the axis, larger values move it down.

    Returns:
        fig, baxes: the figure and one BrokenAxes object per group.
    """
    from brokenaxes import brokenaxes
    from matplotlib.gridspec import GridSpec
    from matplotlib.patches import ConnectionPatch
    from matplotlib.ticker import MaxNLocator

    color_dict = color_dict or {}
    names = _species_names(doc)
    groups = _group_measurements(doc.measurements)
    gids = list(groups.keys()) or ["ungrouped"]

    r, c = _nice_grid(len(gids))
    fig = plt.figure(figsize=(figsize_per_panel[0] * c, figsize_per_panel[1] * r))
    grid = GridSpec(r, c, figure=fig, hspace=0.4, wspace=0.25)
    baxes = []

    for k, gid in enumerate(gids):
        bax = brokenaxes(
            xlims=tuple(x_ranges), subplot_spec=grid[k], fig=fig, wspace=wspace, d=breakmark_size
        )
        baxes.append(bax)
        segments = bax.axs
        agg = _aggregate_group(groups[gid])

        # keep the order of color_dict, species without a color afterwards
        def order(sid: str) -> int:
            keys = list(color_dict)
            for key in (names.get(sid, sid), sid):
                if key in keys:
                    return keys.index(key)
            return len(keys)

        species = [
            sid for sid in sorted(agg, key=order)
            if not any(e in names.get(sid, sid) or e in sid for e in exclude)
        ]

        xlab, ylab = xlabel, ylabel
        legend_entries = []
        for i, sid in enumerate(species):
            t, m, s, meta = agg[sid]
            name = names.get(sid, sid)
            color = color_dict.get(name, color_dict.get(sid, f"C{i % 10}"))
            if xlab is None:
                tu = meta["time_unit"]
                xlab = f"Time [{tu}]" if tu else "Time"
            if ylab is None:
                du, dt = meta["data_unit"], meta["data_type"]
                dt_name = getattr(dt, "name", dt)
                ylab = f"{str(dt_name).capitalize()} [{du}]" if dt_name and du else (du or "Value")

            # last point of the previous segment, connected to the first point of the next one
            prev = None
            for j, (ax, (lo, hi)) in enumerate(zip(segments, x_ranges)):
                # skip missing values (NaN), e.g. time points where no peak was found
                mask = (t >= lo) & (t <= hi) & np.isfinite(m)
                if not mask.any():
                    continue
                errorbar = ax.errorbar(t[mask], m[mask], yerr=np.nan_to_num(s[mask]), alpha=0.5,
                                       capsize=3, color=color)
                points = ax.scatter(t[mask], m[mask], zorder=10, clip_on=False, color=color)
                if prev is None:
                    # legend entry shows the line with error bar and the marker on top
                    legend_entries.append(((errorbar, points), name))
                if prev is not None:
                    con = ConnectionPatch(
                        xyA=prev[1], coordsA="data", axesA=prev[0],
                        xyB=(t[mask][0], m[mask][0]), coordsB="data", axesB=ax,
                        color=color, alpha=0.5, linestyle=(0, (4, 6)),
                    )
                    fig.add_artist(con)
                prev = (ax, (t[mask][-1], m[mask][-1]))

        bax.standardize_ticks()
        bax.set_spines()
        bax.set_ylim(bottom=0)
        for ax in segments:
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        bax.set_ylabel(ylab or "Value")
        bax.set_xlabel(xlab or "Time", labelpad=xlabel_pad)
        bax.set_title(str(doc.name) if len(gids) == 1 and doc.name else str(gid))
        if legend_entries:
            handles, labels = zip(*legend_entries)
            bax.legend(list(handles), list(labels), loc="upper right")

    return fig, baxes
