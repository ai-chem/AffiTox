"""Per-target scatter plots: predicted score (x) vs experimental pKi (y)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.stats import linregress, pearsonr, spearmanr

from ..config import AnalysisConfig
from ..data import load_merged, prepare_xy
from .style import FS_ANNOT, FS_AXIS, FS_PANEL, FS_TICK, apply_style, method_color, save_figure


def plot_target_scatter(target: str, cfg: AnalysisConfig) -> None:
    df = load_merged(cfg.merged_dir, target)
    n_methods = len(cfg.methods)
    apply_style(cfg)

    series: list[tuple[str, str, np.ndarray, np.ndarray]] = []
    y_values: list[np.ndarray] = []
    for method_id in cfg.methods:
        label = cfg.method_label(method_id)
        metric = cfg.primary_metric(method_id)
        if metric not in df.columns:
            series.append((method_id, label, np.array([]), np.array([])))
            continue
        pki, score = prepare_xy(df, metric, cfg.exp_col)
        series.append((method_id, label, pki, score))
        if len(pki) >= 3:
            y_values.append(pki)

    if y_values:
        y_all = np.concatenate(y_values)
        y_span = float(y_all.max() - y_all.min()) or 1.0
        y_pad = 0.05 * y_span
        y_lim = (float(y_all.min()) - y_pad, float(y_all.max()) + y_pad)
    else:
        y_lim = None

    fig, axes = plt.subplots(
        1,
        n_methods,
        figsize=(3.7 * n_methods + 0.7, 5.55),
        sharey=True,
        squeeze=False,
        gridspec_kw={"wspace": 0.07},
    )
    axes = axes[0]

    for i, (ax, (method_id, label, pki, score)) in enumerate(zip(axes, series)):
        color = method_color(cfg, method_id)
        ax.set_box_aspect(1)
        ax.set_title(label, fontsize=FS_PANEL, pad=6)

        if len(pki) < 3:
            ax.text(
                0.5, 0.5,
                "insufficient data" if len(pki) else "no data",
                ha="center", va="center", transform=ax.transAxes,
            )
        else:
            ax.scatter(
                score, pki, s=26, alpha=0.55, color=color,
                edgecolors="white", linewidths=0.35,
            )
            slope, intercept, _, _, _ = linregress(score, pki)
            xx = np.linspace(score.min(), score.max(), 100)
            ax.plot(xx, slope * xx + intercept, color="black", lw=1.7, alpha=0.85)
            pr, _ = pearsonr(score, pki)
            sr, _ = spearmanr(score, pki)
            ax.text(
                0.04, 0.97,
                f"r={pr:.3f}\n$\\rho$={sr:.3f}\nn={len(pki)}",
                transform=ax.transAxes, va="top", ha="left", fontsize=FS_ANNOT,
                bbox=dict(
                    boxstyle="round,pad=0.25", facecolor="white",
                    alpha=0.92, edgecolor="#cccccc",
                ),
            )

        ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=4))
        ax.set_xlabel("Predicted score", fontsize=FS_AXIS, labelpad=6)
        ax.xaxis.set_label_coords(0.5, -0.14)
        ax.tick_params(axis="x", labelsize=FS_TICK, pad=3, rotation=0)
        ax.grid(True, alpha=0.25)
        if y_lim is not None:
            ax.set_ylim(y_lim)
        if i == 0:
            ax.set_ylabel("Experimental pKi", fontsize=FS_AXIS, labelpad=6)
            ax.tick_params(axis="y", labelsize=FS_TICK, pad=2, labelleft=True)
        else:
            ax.set_ylabel("")
            ax.tick_params(axis="y", labelleft=False, length=0)

    fig.subplots_adjust(top=0.90, bottom=0.16, left=0.055, right=0.995)
    out = cfg.figures_dir / "correlations" / "scatter" / f"scatter_{target.lower()}"
    save_figure(fig, out)


def plot_all_scatters(cfg: AnalysisConfig) -> None:
    for target in cfg.targets:
        plot_target_scatter(target, cfg)
