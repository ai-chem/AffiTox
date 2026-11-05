#!/usr/bin/env python3
"""Pass-all PoseBusters heatmap: methods × targets, Pearson square-cell style."""

from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.collections import QuadMesh
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec


ROOT = Path("/nfs/home/okonovalova/AffiTox")
ANALYSIS_ROOT = ROOT / "analysis" / "excluding_2z5x_3mjg"
TB = ROOT / "Toxicity_Benchmark"
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
METHODS = ["boltz2", "dynamicbind", "gnina", "qvina"]  # PLAPT has no PoseBusters
METHOD_LABELS = {
    "boltz2": "Boltz-2",
    "dynamicbind": "DynamicBind",
    "gnina": "GNINA 1.3",
    "qvina": "QVina2",
}


def _annotations(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.map(lambda value: "" if pd.isna(value) else f"{value:.1f}")


def _fit_horizontal_ticks(fig: plt.Figure, ax, start_fs: float, min_fs: float = 7.5) -> float:
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
    table_path = ANALYSIS_ROOT / "tables" / "posebusters" / "pass_rates_by_target_method.csv"
    summary = pd.read_csv(table_path)
    summary = summary[summary["method_id"].isin(METHODS)].copy()

    rates = summary.pivot(index="method_id", columns="target", values="pass_rate_all")
    rates = rates.reindex(index=METHODS, columns=TARGETS).astype(float) * 100.0
    rates.index = [METHOD_LABELS[m] for m in rates.index]
    rates.columns = [t.upper() for t in rates.columns]

    # Weighted "All targets" column (same as combined panel B).
    all_targets = []
    for method_id in METHODS:
        part = summary[summary["method_id"] == method_id]
        n = part["n_poses"].sum()
        if n == 0:
            all_targets.append(float("nan"))
        else:
            all_targets.append(100.0 * (part["pass_rate_all"] * part["n_poses"]).sum() / n)
    all_col = pd.DataFrame({"All targets": all_targets}, index=rates.index)

    # Full article palette, then remap so 0↔old(15%) and 100↔old(90%).
    base_cmap = LinearSegmentedColormap.from_list(
        "ultramarine_ivory_vermilion",
        ["#3558C5", "#FFF9EE", "#EF4938"],
        N=256,
    )
    lo, hi = 0.15, 0.90
    cmap = LinearSegmentedColormap.from_list(
        "pass_rate_compressed",
        [base_cmap(lo + (hi - lo) * i / 255.0) for i in range(256)],
        N=256,
    )
    # Same 4-size ladder as Pearson square-cell figure (×1.9).
    axis_fs = 19.0
    method_fs = 17.1
    tick_fs = 17.1
    annot_fs = 15.2
    cbar_tick_fs = 17.1
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": tick_fs,
            "axes.labelsize": axis_fs,
            "axes.titlesize": 22.8,
            "xtick.labelsize": tick_fs,
            "ytick.labelsize": method_fs,
            "font.weight": "normal",
            "axes.labelweight": "normal",
            "axes.titleweight": "normal",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )

    # 4 methods vs Pearson's 5 → scale height ~4/5.
    fig = plt.figure(figsize=(24.0, 8.8 * 4 / 5))
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
    summary_ax = fig.add_subplot(grid[0, 0])
    body_ax = fig.add_subplot(grid[0, 2])
    colorbar_ax = fig.add_subplot(grid[0, 3])

    common = {
        "cmap": cmap,
        "vmin": 0.0,
        "vmax": 100.0,
        "linewidths": 0.4,
        "linecolor": "white",
        "square": True,
        "annot_kws": {"size": annot_fs, "color": "black"},
    }
    sns.heatmap(
        all_col,
        ax=summary_ax,
        annot=_annotations(all_col),
        fmt="",
        cbar=False,
        **common,
    )
    sns.heatmap(
        rates,
        ax=body_ax,
        annot=_annotations(rates),
        fmt="",
        yticklabels=False,
        cbar=True,
        cbar_ax=colorbar_ax,
        cbar_kws={"label": "Pass rate (%)"},
        **common,
    )

    summary_ax.set_xlabel("")
    summary_ax.set_ylabel("Docking method", fontsize=axis_fs)
    summary_ax.set_xticklabels(["All\ntargets"], rotation=0, ha="center", va="top")
    summary_ax.tick_params(axis="x", rotation=0, labelsize=tick_fs, pad=6)
    summary_ax.tick_params(axis="y", rotation=0, labelsize=method_fs)
    body_ax.set_xlabel("Target", fontsize=axis_fs, labelpad=8)
    body_ax.set_ylabel("")
    body_ax.tick_params(axis="x", rotation=0, labelsize=tick_fs, pad=4)
    body_ax.tick_params(axis="y", left=False, labelleft=False)
    plt.setp(body_ax.get_xticklabels(), rotation=0, ha="center", va="top")
    plt.setp(summary_ax.get_xticklabels(), linespacing=0.95)

    cbar = body_ax.collections[0].colorbar
    bottom, top = cbar.ax.get_ylim()
    if bottom < top:
        cbar.ax.invert_yaxis()
    cbar.set_ticks([0, 20, 40, 60, 80, 100])
    cbar.ax.tick_params(labelsize=cbar_tick_fs)
    cbar.set_label("Pass rate (%)", fontsize=axis_fs)

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

    name = "pass_all_checks_by_method_target_square_cells"
    output_base = ANALYSIS_ROOT / "figures" / "posebusters" / name
    output_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(
        output_base.with_suffix(".png"),
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    fig.savefig(
        output_base.with_suffix(".pdf"),
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)

    for dest_dir in (
        TB,
        TB / "by_pdb" / "posebusters",
        ANALYSIS_ROOT / "figures" / "by_pdb" / "posebusters",
    ):
        dest_dir.mkdir(parents=True, exist_ok=True)
        for suffix in (".png", ".pdf", ".svg"):
            src = output_base.with_suffix(suffix)
            if src.exists():
                shutil.copy2(src, dest_dir / src.name)

    print(output_base.with_suffix(".png"))


if __name__ == "__main__":
    main()
