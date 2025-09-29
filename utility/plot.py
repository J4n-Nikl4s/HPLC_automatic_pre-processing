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
