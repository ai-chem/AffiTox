#!/usr/bin/env python3
"""Pearson heatmap with square cells and BMC-minimum PNG size (>=1500x1200)."""

from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.collections import QuadMesh
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec


ROOT = Path("/mnt/tank/scratch/okonovalova/AffiTox")
ANALYSIS_ROOT = ROOT / "analysis" / "excluding_2z5x_3mjg"
TARGETS = [
    "1g5m",
    "2v5z",
    "3eyg",
    "3jy9",
    "3lxk",
    "11ue",
    "4ase",
    "4f65",
    "4tz4",
    "4zau",
    "5jkv",
    "5mo4",
    "6gqj",
    "6jok",
    "5lf3",
    "7kk3",
]
METHODS = ["boltz2", "dynamicbind", "gnina", "plapt", "qvina"]
METHOD_LABELS = {
    "boltz2": "Boltz-2",
    "dynamicbind": "DynamicBind",
    "gnina": "GNINA 1.3",
    "plapt": "PLAPT",
    "qvina": "QVina2",
}


def _annotations(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.map(lambda value: "" if pd.isna(value) else f"{value:.2f}")


def _fit_horizontal_ticks(fig: plt.Figure, ax, start_fs: float, min_fs: float = 7.5) -> float:
    """Keep x tick labels horizontal and shrink until neighbours do not overlap."""
    labels = [t for t in ax.get_xticklabels() if t.get_text().strip()]
    fs = float(start_fs)
    while fs >= min_fs:
        for text in labels:
            text.set_fontsize(fs)
            text.set_rotation(0)
            text.set_ha("center")
            text.set_va("top")
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        boxes = [text.get_window_extent(renderer) for text in labels]
        overlap = any(boxes[i].x1 + 4.0 > boxes[i + 1].x0 for i in range(len(boxes) - 1))
        if not overlap:
            return fs
        fs -= 0.4
    return min_fs


def main() -> None:
    table_path = ANALYSIS_ROOT / "tables" / "correlations" / "summary_all_proteins.csv"
    correlations = pd.read_csv(table_path)
    pearson = correlations[correlations["method_id"].isin(METHODS)].pivot(
        index="method_id", columns="target", values="pearson_r"
    )
    pearson = pearson.reindex(index=METHODS, columns=TARGETS).astype(float)
    pearson.index = [METHOD_LABELS[method] for method in pearson.index]
    pearson.columns = [target.upper() for target in pearson.columns]
    median = pearson.median(axis=1, skipna=True).to_frame("All targets")

    cmap = LinearSegmentedColormap.from_list(
        "ultramarine_ivory_vermilion",
        ["#3558C5", "#FFF9EE", "#EF4938"],
        N=256,
    )
    # Same 4-size ladder as article figures, scaled for the wide square-cell canvas.
    # Base: FS_PANEL=12, FS_AXIS=10, FS_TICK=9, FS_ANNOT=8  →  ×1.9
    axis_fs = 19.0   # FS_AXIS * 1.9
    method_fs = 17.1  # FS_TICK * 1.9
    tick_fs = 17.1    # FS_TICK * 1.9
    annot_fs = 15.2   # FS_ANNOT * 1.9
    cbar_tick_fs = 17.1
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": tick_fs,
            "axes.labelsize": axis_fs,
            "axes.titlesize": 22.8,  # FS_PANEL * 1.9
            "xtick.labelsize": tick_fs,
            "ytick.labelsize": method_fs,
            "font.weight": "normal",
            "axes.labelweight": "normal",
            "axes.titleweight": "normal",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )

    # Wide canvas so 5x17 square cells stay >=1200 px tall at 300 DPI.
    fig = plt.figure(figsize=(24.0, 8.8))
    grid = GridSpec(
        1,
        4,
        figure=fig,
        width_ratios=[1.0, 0.253, 16.0, 0.55],
        wspace=0.053,
        left=0.12,
        right=0.95,
        bottom=0.26,
        top=0.93,
    )
    median_ax = fig.add_subplot(grid[0, 0])
    body_ax = fig.add_subplot(grid[0, 2])
    colorbar_ax = fig.add_subplot(grid[0, 3])

    common = {
        "cmap": cmap,
        "vmin": -1.0,
        "vmax": 1.0,
        "center": 0.0,
        "linewidths": 0.4,
        "linecolor": "white",
        "square": True,
        "annot_kws": {"size": annot_fs, "color": "black"},
    }
    sns.heatmap(
        median,
        ax=median_ax,
        annot=_annotations(median),
        fmt="",
        cbar=False,
        **common,
    )
    sns.heatmap(
        pearson,
        ax=body_ax,
        annot=_annotations(pearson),
        fmt="",
        yticklabels=False,
        cbar=True,
        cbar_ax=colorbar_ax,
        cbar_kws={"label": "Correlation coefficient"},
        **common,
    )

    median_ax.set_xlabel("")
    median_ax.set_ylabel("Docking method", fontsize=axis_fs)
    median_ax.set_xticklabels(["All\ntargets"], rotation=0, ha="center", va="top")
    median_ax.tick_params(axis="x", rotation=0, labelsize=tick_fs, pad=6)
    median_ax.tick_params(axis="y", rotation=0, labelsize=method_fs)
    body_ax.set_xlabel("Target", fontsize=axis_fs, labelpad=8)
    body_ax.set_ylabel("")
    body_ax.tick_params(axis="x", rotation=0, labelsize=tick_fs, pad=4)
    body_ax.tick_params(axis="y", left=False, labelleft=False)
    plt.setp(body_ax.get_xticklabels(), rotation=0, ha="center", va="top")
    plt.setp(median_ax.get_xticklabels(), linespacing=0.95)

    cbar = body_ax.collections[0].colorbar
    bottom, top = cbar.ax.get_ylim()
    if bottom < top:
        cbar.ax.invert_yaxis()
    cbar.set_ticks([-1.0, 0.0, 1.0])
    cbar.ax.tick_params(labelsize=cbar_tick_fs)
    cbar.set_label("Correlation coefficient", fontsize=axis_fs)

    _fit_horizontal_ticks(fig, body_ax, start_fs=tick_fs, min_fs=tick_fs)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    mesh = next(c for c in body_ax.collections if isinstance(c, QuadMesh))
    table_bbox = mesh.get_window_extent(renderer)
    inv = fig.transFigure.inverted()
    _, y0 = inv.transform((table_bbox.x0, table_bbox.y0))
    _, y1 = inv.transform((table_bbox.x1, table_bbox.y1))
    cpos = colorbar_ax.get_position()
    colorbar_ax.set_position([cpos.x0, y0, cpos.width, y1 - y0])

    output_base = (
        ANALYSIS_ROOT
        / "figures"
        / "correlations"
        / "heatmaps"
        / "pearson_heatmap_ultramarine_ivory_vermilion_square_cells"
    )
    output_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(
        output_base.with_suffix(".png"),
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)

    for variant in ("by_protein", "by_pdb"):
        dest = (
            ANALYSIS_ROOT
            / "figures"
            / variant
            / "correlations"
            / "heatmaps"
            / output_base.name
        )
        dest.parent.mkdir(parents=True, exist_ok=True)
        for suffix in (".png", ".svg"):
            src = output_base.with_suffix(suffix)
            if src.exists():
                shutil.copy2(src, dest.with_suffix(suffix))
    print(output_base.with_suffix(".png"))


if __name__ == "__main__":
    main()
