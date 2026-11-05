"""Correlation heatmaps (article palette; median summary, not pooled)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kendalltau, pearsonr, spearmanr

from ..config import AnalysisConfig
from ..data import load_merged, prepare_xy
from ..pki_range_audit import run_pki_range_audit
from ..property_baselines import (
    baseline_pearson_matrix,
    baseline_summary_medians,
    run_property_baselines,
)
from .style import (
    FS_ANNOT,
    FS_AXIS,
    FS_PANEL,
    FS_TICK,
    apply_style,
    diverging_cmap,
    method_color,
    plot_combined_heatmaps_with_summary_row,
    plot_heatmap_with_summary_col,
    plot_heatmap_with_summary_row,
    save_figure,
)

# Explicit labels so pooled vs median cannot be confused (W5).
MEDIAN_LABEL = "Median"
POOLED_LABEL = "Pooled"


def _pivot(corr: pd.DataFrame, value_col: str, cfg: AnalysisConfig) -> pd.DataFrame:
    labels = [cfg.method_label(m) for m in cfg.methods]
    idx = [cfg.target_label(t) for t in cfg.targets]
    mat = pd.DataFrame(index=idx, columns=labels, dtype=float)
    for _, row in corr.iterrows():
        if row["method_id"] not in cfg.methods:
            continue
        mat.loc[
            cfg.target_label(str(row["target"])),
            cfg.method_label(row["method_id"]),
        ] = row[value_col]
    return mat


def _pooled_correlations(cfg: AnalysisConfig) -> dict[str, dict[str, float]]:
    """Pooled pKi vs score across all ligands/targets (between-target variance included)."""
    pooled: dict[str, dict[str, float]] = {}
    frames = []
    for target in cfg.targets:
        try:
            frames.append(load_merged(cfg.merged_dir, target))
        except FileNotFoundError:
            continue
    if not frames:
        return pooled

    merged = pd.concat(frames, ignore_index=True)
    for method_id in cfg.methods:
        metric = cfg.primary_metric(method_id)
        if metric not in merged.columns:
            continue
        x, y = prepare_xy(merged, metric, cfg.exp_col)
        label = cfg.method_label(method_id)
        if len(x) < 3:
            pooled[label] = {"pearson_r": np.nan, "spearman_rho": np.nan, "kendall_tau": np.nan}
            continue
        pr, _ = pearsonr(x, y)
        sr, _ = spearmanr(x, y)
        kt, _ = kendalltau(x, y)
        pooled[label] = {
            "pearson_r": float(pr),
            "spearman_rho": float(sr),
            "kendall_tau": float(kt),
        }
    return pooled


def _median_summary(body: pd.DataFrame) -> pd.Series:
    """Cross-target median of within-target coefficients (matches main-text 0.552 claim)."""
    return body.median(axis=0, skipna=True)


def plot_correlation_heatmaps(corr: pd.DataFrame, cfg: AnalysisConfig) -> None:
    apply_style(cfg)
    cmap = diverging_cmap(cfg)
    out_dir = cfg.figures_dir / "correlations" / "heatmaps"
    pooled = _pooled_correlations(cfg)

    specs = [
        ("pearson_r", "Pearson r"),
        ("spearman_rho", "Spearman rho"),
        ("kendall_tau", "Kendall tau"),
    ]
    combined_panels: list[tuple[str, pd.Series, pd.DataFrame]] = []
    median_lookup: dict[str, pd.Series] = {}
    for col, title in specs:
        body = _pivot(corr, col, cfg)
        if body.isna().all().all():
            continue
        median = _median_summary(body.astype(float))
        median_lookup[col] = median
        pooled_series = pd.Series(
            {m: pooled.get(m, {}).get(col, np.nan) for m in body.columns},
            dtype=float,
        )

        # Main heatmaps use MEDIAN (not pooled) — W5
        plot_heatmap_with_summary_row(
            median,
            body.astype(float),
            title=f"{title} by target and method",
            cbar_label="Correlation coefficient",
            cmap=cmap,
            vmin=-1.0,
            vmax=1.0,
            center=0.0,
            fmt=".2f",
            path_base=out_dir / f"heatmap_{col}",
            summary_label=MEDIAN_LABEL,
        )
        if col == "pearson_r":
            plot_heatmap_with_summary_col(
                median,
                body.astype(float).T,
                title=f"{title} by method and target",
                cbar_label="Correlation coefficient",
                cmap=cmap,
                vmin=-1.0,
                vmax=1.0,
                center=0.0,
                fmt=".2f",
                body_annot=True,
                y_label="Docking method",
                figsize=(11.8, 4.6),
                equal_summary_width=True,
                label_panel_ratio=1.35,
                label_summary_wspace=0.005,
                cbar_width_ratio=0.35,
                summary_body_gap_scale=0.5,
                summary_body_wspace=0.14,
                y_tick_rotation=0,
                path_base=out_dir / "heatmap_pearson_r_by_method",
                summary_label=MEDIAN_LABEL,
            )
            # Optional SI companion: pooled column, explicitly labeled
            plot_heatmap_with_summary_col(
                pooled_series,
                body.astype(float).T,
                title=f"{title} by method and target (pooled summary)",
                cbar_label="Correlation coefficient",
                cmap=cmap,
                vmin=-1.0,
                vmax=1.0,
                center=0.0,
                fmt=".2f",
                body_annot=True,
                y_label="Docking method",
                figsize=(11.8, 4.6),
                equal_summary_width=True,
                label_panel_ratio=1.35,
                label_summary_wspace=0.005,
                cbar_width_ratio=0.35,
                summary_body_gap_scale=0.5,
                summary_body_wspace=0.14,
                y_tick_rotation=0,
                path_base=out_dir / "heatmap_pearson_r_by_method_pooled",
                summary_label=POOLED_LABEL,
            )
            _plot_pearson_with_property_baselines(
                cfg, body.astype(float), median, cmap, out_dir
            )
        combined_panels.append((title, median, body.astype(float)))

    # Fill median rows in summary table
    summary_records = []
    for method_id in cfg.methods:
        label = cfg.method_label(method_id)
        med = {
            "method_id": method_id,
            "method_label": label,
            "summary": "median_across_targets",
            "pearson_r": float(median_lookup.get("pearson_r", pd.Series()).get(label, np.nan)),
            "spearman_rho": float(median_lookup.get("spearman_rho", pd.Series()).get(label, np.nan)),
            "kendall_tau": float(median_lookup.get("kendall_tau", pd.Series()).get(label, np.nan)),
        }
        pool = {
            "method_id": method_id,
            "method_label": label,
            "summary": "pooled_all_ligands",
            "pearson_r": pooled.get(label, {}).get("pearson_r", np.nan),
            "spearman_rho": pooled.get(label, {}).get("spearman_rho", np.nan),
            "kendall_tau": pooled.get(label, {}).get("kendall_tau", np.nan),
        }
        summary_records.extend([med, pool])
    pd.DataFrame(summary_records).to_csv(
        cfg.tables_dir / "correlations" / "median_vs_pooled_summary.csv",
        index=False,
    )

    if len(combined_panels) == len(specs):
        plot_combined_heatmaps_with_summary_row(
            combined_panels,
            cbar_label="Correlation coefficient",
            cmap=cmap,
            vmin=-1.0,
            vmax=1.0,
            center=0.0,
            fmt=".2f",
            path_base=out_dir / "heatmap_correlations_combined",
            col_width=1.28,
            row_height=0.38,
            summary_label=MEDIAN_LABEL,
        )

    # W5: r vs pKi spread
    panel = run_pki_range_audit(cfg, corr)
    plot_r_vs_pki_spread(panel, cfg)


def _plot_pearson_with_property_baselines(
    cfg: AnalysisConfig,
    body_targets_by_method: pd.DataFrame,
    method_summary: pd.Series,
    cmap,
    out_dir: Path,
) -> None:
    """Figure-2 style heatmap: method rows + HAC/MW/cLogP/permutation-floor rows."""
    baseline_path = cfg.tables_dir / "correlations" / "property_null_baselines.csv"
    if baseline_path.exists():
        baselines = pd.read_csv(baseline_path)
    else:
        baselines = run_property_baselines(cfg)

    base_mat = baseline_pearson_matrix(baselines, cfg)
    base_summary = baseline_summary_medians(baselines)

    method_body = body_targets_by_method.T.astype(float)
    base_mat = base_mat.reindex(columns=list(method_body.columns))
    combined_body = pd.concat([method_body, base_mat], axis=0)
    combined_summary = pd.concat([method_summary.astype(float), base_summary])
    combined_summary = combined_summary.reindex(combined_body.index)

    n_rows = len(combined_body.index)
    plot_heatmap_with_summary_col(
        combined_summary,
        combined_body,
        title="Pearson r by method/target with property null baselines",
        cbar_label="Correlation coefficient",
        cmap=cmap,
        vmin=-1.0,
        vmax=1.0,
        center=0.0,
        fmt=".2f",
        body_annot=True,
        y_label="Method / null baseline",
        figsize=(12.4, max(5.2, 0.42 * n_rows + 1.2)),
        equal_summary_width=True,
        label_panel_ratio=1.55,
        label_summary_wspace=0.005,
        cbar_width_ratio=0.35,
        summary_body_gap_scale=0.5,
        summary_body_wspace=0.14,
        y_tick_rotation=0,
        path_base=out_dir / "heatmap_pearson_r_with_baselines",
        summary_label=MEDIAN_LABEL,
    )


def plot_r_vs_pki_spread(panel: pd.DataFrame, cfg: AnalysisConfig) -> None:
    """Scatter: per-target Pearson r vs pKi SD (and IQR) for primary methods."""
    apply_style(cfg)
    out_dir = cfg.figures_dir / "correlations" / "pki_range"
    out_dir.mkdir(parents=True, exist_ok=True)

    methods = [(mid, cfg.method_label(mid)) for mid in cfg.methods]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), sharey=True)
    for ax, spread_col, xlab in zip(
        axes,
        ("pki_sd", "pki_iqr"),
        (r"Within-target p$K_i$ SD", r"Within-target p$K_i$ IQR"),
    ):
        for method_id, label in methods:
            col = f"pearson_r_{method_id}"
            if col not in panel.columns:
                continue
            x = panel[spread_col].to_numpy(dtype=float)
            y = panel[col].to_numpy(dtype=float)
            mask = np.isfinite(x) & np.isfinite(y)
            color = method_color(cfg, method_id)
            ax.scatter(
                x[mask],
                y[mask],
                s=36,
                alpha=0.85,
                color=color,
                label=label,
                edgecolors="white",
                linewidths=0.4,
                zorder=3,
            )
            # Annotate Boltz-2 points with PDB ids
            if method_id == "boltz2":
                for _, row in panel.loc[mask].iterrows():
                    ax.annotate(
                        cfg.target_label(str(row["target"])),
                        (row[spread_col], row[col]),
                        textcoords="offset points",
                        xytext=(3, 3),
                        fontsize=FS_ANNOT,
                        color=color,
                        alpha=0.9,
                    )
        ax.axhline(0.0, color="#888888", lw=0.8, zorder=1)
        ax.set_xlabel(xlab)
        ax.set_xlim(left=0)
        ax.grid(True, alpha=0.25)
    axes[0].set_ylabel(r"Per-target Pearson $r$")
    axes[1].legend(loc="lower right", fontsize=FS_ANNOT, frameon=False)
    fig.suptitle("Scoring power vs within-target affinity spread", y=1.02, fontsize=FS_PANEL)
    fig.tight_layout()
    save_figure(fig, out_dir / "pearson_r_vs_pki_spread")

    # Boltz-only focused panel (main-text / SI companion)
    if "pearson_r_boltz2" not in panel.columns:
        return
    fig2, ax = plt.subplots(figsize=(5.6, 4.4))
    x = panel["pki_sd"].to_numpy(dtype=float)
    y = panel["pearson_r_boltz2"].to_numpy(dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    color = method_color(cfg, "boltz2")
    ax.scatter(x[mask], y[mask], s=55, color=color, edgecolors="white", linewidths=0.5, zorder=3)
    for _, row in panel.loc[mask].iterrows():
        ax.annotate(
            cfg.target_label(str(row["target"])),
            (row["pki_sd"], row["pearson_r_boltz2"]),
            textcoords="offset points",
            xytext=(4, 4),
            fontsize=FS_ANNOT,
        )
    if mask.sum() >= 4:
        sr, sp = spearmanr(x[mask], y[mask])
        ax.text(
            0.04,
            0.96,
            rf"Spearman $\rho$($r$, SD)$=${sr:.2f}" + "\n" + rf"$p=${sp:.3f}",
            transform=ax.transAxes,
            va="top",
            fontsize=FS_ANNOT,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#cccccc", alpha=0.9),
        )
    ax.axhline(0.0, color="#888888", lw=0.8)
    ax.set_xlabel(r"Within-target p$K_i$ SD")
    ax.set_ylabel(r"Boltz-2 Pearson $r$")
    ax.set_title("Is target-dependence just affinity-range attenuation?")
    ax.grid(True, alpha=0.25)
    fig2.tight_layout()
    save_figure(fig2, out_dir / "boltz2_pearson_r_vs_pki_sd")
