"""PoseBusters heatmaps with All targets summary row/column."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec

from ..config import AnalysisConfig
from ..constants import POSEBUSTERS_EXCLUDE
from ..posebusters import collect_pooled_pass_counts, pass_count_distribution
from .style import (
    ALL_TARGETS_LABEL,
    FS_ANNOT,
    FS_AXIS,
    FS_PANEL,
    FS_TICK,
    SUMMARY_BODY_HSPACE,
    apply_style,
    method_color,
    plot_combined_heatmaps_with_summary_row,
    plot_heatmap_with_summary_row,
    plot_heatmap_with_summary_row_and_col,
    save_figure,
    sequential_cmap,
)
from .style import _format_annot  # shared heatmap cell labels

ALL_TESTS_HEATMAP_VMIN = 76.0  # 40–100 base; color at 75 → now at 90 (vmax=100)
ALL_TESTS_HEATMAP_VMAX = 100.0

ARTICLE_CHECK_LABELS = {
    "mol_pred_loaded": "Molecule loaded",
    "mol_cond_loaded": "Condition loaded",
    "sanitization": "Sanitization",
    "inchi_convertible": "InChI conversion",
    "all_atoms_connected": "All atoms connected",
    "bond_lengths": "Bond lengths",
    "bond_angles": "Bond angles",
    "internal_steric_clash": "Internal steric clash",
    "aromatic_ring_flatness": "Aromatic ring flatness",
    "non-aromatic_ring_non-flatness": "Non-aromatic ring non-flatness",
    "double_bond_flatness": "Double bond flatness",
    "internal_energy": "Internal energy",
    "protein-ligand_maximum_distance": "Protein–ligand max. distance",
    "minimum_distance_to_protein": "Min. distance to protein",
    "minimum_distance_to_organic_cofactors": "Min. distance to organic cofactors",
    "minimum_distance_to_inorganic_cofactors": "Min. distance to inorganic cofactors",
    "minimum_distance_to_waters": "Min. distance to waters",
    "volume_overlap_with_protein": "Volume overlap with protein",
    "volume_overlap_with_organic_cofactors": "Volume overlap with organic cofactors",
    "volume_overlap_with_inorganic_cofactors": "Volume overlap with inorganic cofactors",
    "volume_overlap_with_waters": "Volume overlap with waters",
}


def _style_heatmap_text_contrast(
    ax: plt.Axes, vmin: float, vmax: float, *, fontsize: float = FS_ANNOT
) -> None:
    """White/dark annotation text tuned for the 76–100% pass-rate palette."""
    for text in ax.texts:
        raw = text.get_text().strip()
        if not raw:
            continue
        try:
            val = float(raw)
        except ValueError:
            continue
        t = (val - vmin) / (vmax - vmin) if vmax > vmin else 0.5
        text.set_color("#ffffff" if t >= 0.50 else "#222222")
        text.set_fontsize(fontsize)


def _add_panel_letter(ax: plt.Axes, letter: str, fontsize: float) -> None:
    """Place A/B the same number of points above the axes top (not axes-fraction)."""
    ax.annotate(
        letter,
        xy=(0.0, 1.0),
        xycoords=ax.transAxes,
        xytext=(-2, 6),
        textcoords="offset points",
        fontsize=fontsize,
        fontweight="bold",
        ha="right",
        va="bottom",
        clip_on=False,
        annotation_clip=False,
    )


def _draw_all_tests_heatmap_panel(
    fig: plt.Figure,
    gs_slot,
    *,
    summary_df: pd.DataFrame,
    body: pd.DataFrame,
    n_body: int,
    cmap,
    title: str,
    adaptive_text: bool,
    vmin: float = ALL_TESTS_HEATMAP_VMIN,
    vmax: float = ALL_TESTS_HEATMAP_VMAX,
    cbar_ticks: list[float] | None = None,
    panel_letter: str | None = None,
    title_fontsize: float | None = None,
    annot_size: float = FS_ANNOT,
    show_cbar: bool = True,
) -> plt.Axes:
    if show_cbar:
        inner = gs_slot.subgridspec(
            2,
            2,
            height_ratios=[1.0, max(n_body, 1)],
            width_ratios=[1.0, 0.052],
            hspace=SUMMARY_BODY_HSPACE / 4,
            wspace=0.05,
        )
        ax_top = fig.add_subplot(inner[0, 0])
        ax_bot = fig.add_subplot(inner[1, 0], sharex=ax_top)
        cbar_ax = fig.add_subplot(inner[:, 1])
    else:
        inner = gs_slot.subgridspec(
            2,
            1,
            height_ratios=[1.0, max(n_body, 1)],
            hspace=SUMMARY_BODY_HSPACE / 4,
        )
        ax_top = fig.add_subplot(inner[0, 0])
        ax_bot = fig.add_subplot(inner[1, 0], sharex=ax_top)
        cbar_ax = None

    if cbar_ticks is None:
        cbar_ticks = [76, 80, 85, 90, 95, 100]
    heatmap_kw = dict(
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        linewidths=0.5,
        linecolor="white",
        annot_kws={"size": annot_size, "color": "black"},
    )
    sns.heatmap(
        summary_df,
        ax=ax_top,
        annot=_format_annot(summary_df, ".1f"),
        fmt="",
        cbar=False,
        **heatmap_kw,
    )
    sns.heatmap(
        body,
        ax=ax_bot,
        annot=_format_annot(body, ".1f"),
        fmt="",
        cbar=show_cbar,
        cbar_ax=cbar_ax,
        cbar_kws={"label": "Pass rate (%)", "ticks": cbar_ticks} if show_cbar else None,
        **heatmap_kw,
    )
    if adaptive_text:
        _style_heatmap_text_contrast(ax_top, vmin, vmax, fontsize=annot_size)
        _style_heatmap_text_contrast(ax_bot, vmin, vmax, fontsize=annot_size)

    title_fs = FS_PANEL if title_fontsize is None else title_fontsize
    if title:
        ax_top.set_title(title, fontsize=title_fs, pad=6)
    ax_top.set_xlabel("")
    ax_top.set_ylabel("")
    ax_top.tick_params(axis="x", bottom=False, labelbottom=False, labelsize=FS_TICK)
    ax_top.tick_params(axis="y", rotation=0, labelsize=FS_TICK)

    ax_bot.set_xlabel("")
    ax_bot.set_ylabel("Target", fontsize=FS_AXIS)
    ax_bot.tick_params(axis="x", rotation=0, labelsize=FS_TICK)
    ax_bot.tick_params(axis="y", rotation=0, pad=4, labelsize=FS_TICK)
    ax_top.tick_params(axis="y", pad=4, labelsize=FS_TICK)
    if cbar_ax is not None:
        cbar_ax.tick_params(labelsize=FS_TICK)
    if panel_letter:
        _add_panel_letter(ax_top, panel_letter, title_fs)
    return ax_top


def _prepare_at_least_all_tests_panels(
    summary: pd.DataFrame,
    checks: pd.DataFrame,
    cfg: AnalysisConfig,
) -> tuple[int, dict[str, pd.DataFrame], pd.DataFrame, pd.DataFrame, int] | None:
    methods = [m for m in cfg.methods if m not in POSEBUSTERS_EXCLUDE]
    n_tests, pooled = collect_pooled_pass_counts(cfg)
    if n_tests == 0 or not any(pooled.values()):
        return None

    at_least: dict[str, pd.DataFrame] = {}
    for method_id, counts in pooled.items():
        if counts:
            at_least[method_id] = pass_count_distribution(counts, n_tests, mode="at_least")

    summary_row = _all_tests_summary_row(checks, summary, cfg, methods)
    body = _all_tests_body(checks, cfg, methods)
    if body.isna().all().all():
        return None

    summary_df = summary_row.to_frame().T
    summary_df.index = [ALL_TARGETS_LABEL]
    return n_tests, at_least, summary_df, body.astype(float), len(body)

def _pretty_check(name: str) -> str:
    if name in ARTICLE_CHECK_LABELS:
        return ARTICLE_CHECK_LABELS[name]
    return (
        name.replace("minimum_distance_to_", "min_dist_")
        .replace("volume_overlap_with_", "overlap_")
        .replace("protein-ligand_", "prot_lig_")
    )


def _pooled_check_pass_rates(
    checks: pd.DataFrame,
    summary: pd.DataFrame,
    cfg: AnalysisConfig,
    methods: list[str],
) -> pd.DataFrame:
    """Pose-weighted pooled pass rate (%) per check × method (article panel A)."""
    labels = [cfg.method_label(m) for m in methods]
    merged = checks.merge(
        summary[["target", "method_id", "n_poses"]],
        on=["target", "method_id"],
        how="left",
    )
    check_order: list[str] = []
    for check in merged["check"]:
        if check not in check_order:
            check_order.append(check)
    records: dict[str, dict[str, float]] = {}
    for check in check_order:
        pretty = _pretty_check(check)
        records[pretty] = {}
        for method_id in methods:
            part = merged[(merged["check"] == check) & (merged["method_id"] == method_id)]
            label = cfg.method_label(method_id)
            n = float(part["n_poses"].fillna(0).sum())
            if n <= 0:
                records[pretty][label] = np.nan
            else:
                records[pretty][label] = 100.0 * float(
                    (part["pass_rate"] * part["n_poses"]).sum() / n
                )
    return pd.DataFrame.from_dict(records, orient="index")[labels]


def _pass_rate_body(summary: pd.DataFrame, rate_col: str, cfg: AnalysisConfig, methods: list[str]) -> pd.DataFrame:
    labels = [cfg.method_label(m) for m in methods]
    idx = [cfg.target_label(t) for t in cfg.targets]
    body = pd.DataFrame(index=idx, columns=labels, dtype=float)
    for _, row in summary.iterrows():
        if row["method_id"] not in methods:
            continue
        body.loc[
            cfg.target_label(str(row["target"])),
            cfg.method_label(row["method_id"]),
        ] = 100.0 * row[rate_col]
    return body


def _weighted_summary_row(summary: pd.DataFrame, rate_col: str, cfg: AnalysisConfig, methods: list[str]) -> pd.Series:
    sub = summary[summary["method_id"].isin(methods)].copy()
    sub["weighted"] = sub[rate_col] * sub["n_poses"]
    summary_row = pd.Series(dtype=float)
    for method_id in methods:
        label = cfg.method_label(method_id)
        part = sub[sub["method_id"] == method_id]
        if part.empty or part["n_poses"].sum() == 0:
            summary_row[label] = np.nan
        else:
            summary_row[label] = 100.0 * part["weighted"].sum() / part["n_poses"].sum()
    return summary_row


def _all_tests_body(checks: pd.DataFrame, cfg: AnalysisConfig, methods: list[str]) -> pd.DataFrame:
    """Mean per-check pass rate by target and method (All tests row from per-check heatmaps)."""
    labels = [cfg.method_label(m) for m in methods]
    idx = [cfg.target_label(t) for t in cfg.targets]
    body = pd.DataFrame(index=idx, columns=labels, dtype=float)
    for method_id in methods:
        sub = checks[checks.method_id == method_id]
        if sub.empty:
            continue
        label = cfg.method_label(method_id)
        for target in cfg.targets:
            part = sub[sub.target == target.lower()]
            if part.empty:
                continue
            body.loc[cfg.target_label(target), label] = 100.0 * part["pass_rate"].mean()
    return body


def _all_tests_summary_row(checks: pd.DataFrame, summary: pd.DataFrame, cfg: AnalysisConfig, methods: list[str]) -> pd.Series:
    if "mean_frac_pass" in summary.columns:
        sub = summary[summary["method_id"].isin(methods)].copy()
        summary_row = pd.Series(dtype=float)
        for method_id in methods:
            label = cfg.method_label(method_id)
            part = sub[sub["method_id"] == method_id]
            if part.empty or part["n_poses"].sum() == 0:
                summary_row[label] = np.nan
            else:
                summary_row[label] = 100.0 * (
                    part["mean_frac_pass"] * part["n_poses"]
                ).sum() / part["n_poses"].sum()
        return summary_row

    body = _all_tests_body(checks, cfg, methods)
    return body.mean(axis=0, skipna=True)


def plot_posebusters_summary(summary: pd.DataFrame, cfg: AnalysisConfig) -> None:
    apply_style(cfg)
    cmap = sequential_cmap(cfg)
    out_dir = cfg.figures_dir / "posebusters" / "heatmaps"
    methods = [m for m in cfg.methods if m not in POSEBUSTERS_EXCLUDE]

    for rate_col, title in (
        ("pass_rate_all", "Pass all checks (100%)"),
        ("pass_rate_90pct", "Pass >=90% checks"),
        ("pass_rate_50pct", "Pass >=50% checks"),
    ):
        body = _pass_rate_body(summary, rate_col, cfg, methods)
        if body.isna().all().all():
            continue

        summary_row = _weighted_summary_row(summary, rate_col, cfg, methods)
        plot_heatmap_with_summary_row(
            summary_row,
            body.astype(float),
            title=title,
            cbar_label="Pass rate (%)",
            cmap=cmap,
            vmin=0.0,
            vmax=100.0,
            center=None,
            fmt=".1f",
            path_base=out_dir / f"heatmap_{rate_col}",
        )


def plot_posebusters_combined(summary: pd.DataFrame, checks: pd.DataFrame, cfg: AnalysisConfig) -> None:
    """Three-panel figure: pass all, pass >=90%, and mean All-tests pass rate."""
    apply_style(cfg)
    cmap = sequential_cmap(cfg)
    out_dir = cfg.figures_dir / "posebusters" / "heatmaps"
    methods = [m for m in cfg.methods if m not in POSEBUSTERS_EXCLUDE]

    panels: list[tuple[str, pd.Series, pd.DataFrame]] = []
    for rate_col, title in (
        ("pass_rate_all", "Pass all checks (100%)"),
        ("pass_rate_90pct", "Pass >=90% checks"),
    ):
        body = _pass_rate_body(summary, rate_col, cfg, methods)
        if body.isna().all().all():
            continue
        panels.append((title, _weighted_summary_row(summary, rate_col, cfg, methods), body.astype(float)))

    all_tests_body = _all_tests_body(checks, cfg, methods)
    if not all_tests_body.isna().all().all():
        panels.append(
            (
                "All tests (mean check pass rate)",
                _all_tests_summary_row(checks, summary, cfg, methods),
                all_tests_body.astype(float),
            )
        )

    if len(panels) == 3:
        plot_combined_heatmaps_with_summary_row(
            panels,
            cbar_label="Pass rate (%)",
            cmap=cmap,
            vmin=0.0,
            vmax=100.0,
            center=None,
            fmt=".1f",
            path_base=out_dir / "heatmap_posebusters_combined",
            col_width=0.95,
            row_height=0.36,
        )


def plot_posebusters_per_check(checks: pd.DataFrame, cfg: AnalysisConfig) -> None:
    apply_style(cfg)
    cmap = sequential_cmap(cfg)
    out_dir = cfg.figures_dir / "posebusters" / "heatmaps"
    methods = [m for m in cfg.methods if m not in POSEBUSTERS_EXCLUDE]

    for method_id in methods:
        sub = checks[checks.method_id == method_id].copy()
        if sub.empty:
            continue
        label = cfg.method_label(method_id)
        body = sub.pivot(index="check", columns="target", values="pass_rate")
        body.index = [_pretty_check(c) for c in body.index]
        body.columns = [cfg.target_label(str(c)) for c in body.columns]
        body = 100.0 * body

        summary_col = body.mean(axis=1, skipna=True)
        summary_row = body.mean(axis=0, skipna=True)

        plot_heatmap_with_summary_row_and_col(
            summary_row,
            summary_col,
            body.astype(float),
            title=f"PoseBusters per-check pass rate — {label}",
            cbar_label="Pass rate (%)",
            cmap=cmap,
            vmin=0.0,
            vmax=100.0,
            fmt=".1f",
            path_base=out_dir / f"heatmap_checks_{method_id}",
            figsize_scale=(0.90, 0.44),
            y_label="Tests",
        )


def _draw_pass_count_curve_on_ax(
    ax: plt.Axes,
    curves: dict[str, pd.DataFrame],
    cfg: AnalysisConfig,
    n_tests: int,
    *,
    x_label: str,
    title: str,
    panel_letter: str | None = None,
    title_fontsize: float | None = None,
) -> None:
    methods = [m for m in cfg.methods if m not in POSEBUSTERS_EXCLUDE]
    for method_id in methods:
        df = curves.get(method_id)
        if df is None or df.empty:
            continue
        ax.plot(
            df["n_tests_passed"],
            df["pct_ligands"],
            label=cfg.method_label(method_id),
            color=method_color(cfg, method_id),
            linestyle="-",
            linewidth=2.0,
            marker="o",
            markersize=3.5,
        )

    ax.set_xlim(0.5, n_tests + 0.5)
    ax.set_ylim(0, 105)
    ax.set_xticks(range(1, n_tests + 1))
    ax.set_yticks(range(0, 101, 10))
    title_fs = FS_PANEL if title_fontsize is None else title_fontsize
    ax.set_xlabel(x_label, fontsize=FS_AXIS)
    ax.set_ylabel("% ligands", fontsize=FS_AXIS)
    ax.tick_params(axis="both", labelsize=FS_TICK)
    if title:
        ax.set_title(title, fontsize=title_fs, pad=6)
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.legend(loc="upper left", frameon=True, framealpha=0.9, fontsize=FS_TICK)
    if panel_letter:
        _add_panel_letter(ax, panel_letter, title_fs)


def _plot_pass_count_curve(
    curves: dict[str, pd.DataFrame],
    *,
    cfg: AnalysisConfig,
    n_tests: int,
    title: str,
    x_label: str,
    path_base,
) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    _draw_pass_count_curve_on_ax(ax, curves, cfg, n_tests, x_label=x_label, title=title)
    save_figure(fig, path_base)


def plot_posebusters_at_least_all_tests_combined(
    summary: pd.DataFrame,
    checks: pd.DataFrame,
    cfg: AnalysisConfig,
) -> None:
    """Article figure: A per-check pooled rates; B pass-all by target (+ All targets).

    Panel B keeps its heatmap width; height is stretched to match A. Horizontal
    gaps A—(B block)—colorbar are equal; cell annotations use ``FS_ANNOT * 1.25``.
    Entire figure (canvas + fonts + strokes) is scaled by ``SCALE``.
    """
    apply_style(cfg)
    cmap = sequential_cmap(cfg)
    out_dir = cfg.figures_dir / "posebusters"
    methods = [m for m in cfg.methods if m not in POSEBUSTERS_EXCLUDE]

    panel_a = _pooled_check_pass_rates(checks, summary, cfg, methods)
    if panel_a.isna().all().all():
        return

    body = _pass_rate_body(summary, "pass_rate_all", cfg, methods)
    if body.isna().all().all():
        return
    summary_row = _weighted_summary_row(summary, "pass_rate_all", cfg, methods)
    summary_df = summary_row.to_frame().T
    summary_df.index = [ALL_TARGETS_LABEL]
    n_body = len(body)

    scale = 1.2
    font_scale = scale * 1.3
    fs_panel = FS_PANEL * font_scale
    fs_axis = FS_AXIS * font_scale
    fs_tick = FS_TICK * font_scale
    annot_size = FS_ANNOT * 1.25 * font_scale
    lw = 0.5 * scale

    fig = plt.figure(figsize=(17.2 * scale, 7.4 * scale))
    outer = GridSpec(
        1,
        3,
        figure=fig,
        width_ratios=[1.35, 1.0, 0.045],
        wspace=0.28,
        left=0.16,
        right=0.95,
        top=0.88,
        bottom=0.10,
    )

    ax_a = fig.add_subplot(outer[0, 0])
    sns.heatmap(
        panel_a,
        ax=ax_a,
        cmap=cmap,
        vmin=0.0,
        vmax=100.0,
        linewidths=lw,
        linecolor="white",
        annot=_format_annot(panel_a, ".1f"),
        fmt="",
        cbar=False,
        annot_kws={"size": annot_size, "color": "black"},
    )
    ax_a.set_xlabel("")
    ax_a.set_ylabel("Test", fontsize=fs_axis)
    ax_a.tick_params(axis="x", rotation=0, labelsize=fs_tick)
    ax_a.tick_params(axis="y", rotation=0, labelsize=fs_tick, pad=2 * scale)
    plt.setp(ax_a.get_xticklabels(), rotation=0, ha="center")

    ax_b_top = _draw_all_tests_heatmap_panel(
        fig,
        outer[0, 1],
        summary_df=summary_df,
        body=body.astype(float),
        n_body=n_body,
        cmap=cmap,
        title="",
        adaptive_text=False,
        vmin=0.0,
        vmax=100.0,
        cbar_ticks=[0, 20, 40, 60, 80, 100],
        title_fontsize=fs_panel,
        annot_size=annot_size,
        show_cbar=False,
    )

    cbar_ax = fig.add_subplot(outer[0, 2])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0.0, vmax=100.0))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax, ticks=[0, 20, 40, 60, 80, 100])
    cbar.set_label("Pass rate (%)", fontsize=fs_axis)
    cbar_ax.tick_params(labelsize=fs_tick)

    fig.canvas.draw()
    ax_b_bot = max(
        (
            ax
            for ax in fig.axes
            if ax is not ax_b_top
            and ax is not ax_a
            and ax is not cbar_ax
            and abs(ax.get_position().x0 - ax_b_top.get_position().x0) < 1e-4
        ),
        key=lambda a: a.get_position().height,
    )

    fig_w, fig_h = fig.get_size_inches()
    pa = ax_a.get_position()
    bt = ax_b_top.get_position()
    bb = ax_b_bot.get_position()
    pc = cbar_ax.get_position()

    # Both tables +1/3 taller; cells scale proportionally; widths unchanged.
    height_scale = 4.0 / 3.0
    a_w_in = pa.width * fig_w
    a_h_nat = pa.height * fig_h
    a_h_in = a_h_nat * height_scale
    b_w_in = bb.width * fig_w
    cbar_w_in = pc.width * fig_w
    gap_sum_in = max(0.0, (bt.y0 - (bb.y0 + bb.height)) * fig_h)
    new_sum_in = bt.height * fig_h * height_scale
    new_gap_in = gap_sum_in * height_scale
    new_body_in = a_h_in - new_sum_in - new_gap_in

    top_margin_in = (1.0 - pa.y0 - pa.height) * fig_h
    bottom_margin_in = pa.y0 * fig_h
    fig_h = top_margin_in + a_h_in + bottom_margin_in
    fig.set_size_inches(fig_w, fig_h, forward=True)
    y1 = fig_h - top_margin_in
    a_y0 = y1 - a_h_in
    b_top_y0 = y1 - new_sum_in
    b_bot_y0 = a_y0

    # Match panel-B fonts/strokes to the scaled ladder (helper still uses FS_*).
    for ax in (ax_b_top, ax_b_bot):
        ax.tick_params(axis="both", labelsize=fs_tick, pad=4 * scale)
        for coll in ax.collections:
            coll.set_linewidths([lw])
        for t in list(ax.texts):
            t.set_fontsize(annot_size)
    ax_b_bot.set_ylabel("Target", fontsize=fs_axis, labelpad=6 * scale)

    # Pass-rate colorbar digits: same size + family as docking method labels.
    method_fs = fs_tick
    method_family = "DejaVu Sans"
    for t in ax_a.get_xticklabels():
        if t.get_text().strip():
            method_fs = float(t.get_fontsize())
            method_family = t.get_fontfamily()
            if isinstance(method_family, (list, tuple)):
                method_family = method_family[0] if method_family else "DejaVu Sans"
            break

    def _style_cbar_ticks() -> None:
        cbar.set_label("Pass rate (%)", fontsize=method_fs)
        if cbar.ax.yaxis.label is not None:
            cbar.ax.yaxis.label.set_fontsize(method_fs)
            cbar.ax.yaxis.label.set_fontname(method_family)
        cbar_ax.tick_params(labelsize=method_fs)
        for t in cbar_ax.get_yticklabels():
            t.set_fontsize(method_fs)
            t.set_fontname(method_family)

    _style_cbar_ticks()
    fig._affitox_cbar_tick_fs = method_fs
    fig._affitox_cbar_tick_family = method_family
    fig._affitox_style_cbar_ticks = _style_cbar_ticks

    # Place A, then set equal heatmap↔heatmap↔colorbar gaps large enough for Target.
    label_left_in = 2.8 * scale
    a_x = label_left_in
    ax_a.set_position([a_x / fig_w, a_y0 / fig_h, a_w_in / fig_w, a_h_in / fig_h])
    # Temporary B so we can measure how much left margin Target+ticks need.
    b_x = a_x + a_w_in + 1.2 * scale
    ax_b_top.set_position(
        [b_x / fig_w, b_top_y0 / fig_h, b_w_in / fig_w, new_sum_in / fig_h]
    )
    ax_b_bot.set_position(
        [b_x / fig_w, b_bot_y0 / fig_h, b_w_in / fig_w, new_body_in / fig_h]
    )
    cbar_ax.set_position(
        [(b_x + b_w_in + 0.6 * scale) / fig_w, b_bot_y0 / fig_h, cbar_w_in / fig_w, a_h_in / fig_h]
    )

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()

    def _axes_box_in(ax):
        p = ax.get_position()
        return p.x0 * fig_w, (p.x0 + p.width) * fig_w

    def _label_left_extent_in(ax) -> float:
        """Inches from axes left edge to leftmost y-label/tick ink."""
        lefts = []
        for artist in list(ax.get_yticklabels()) + [ax.yaxis.label]:
            if artist is None or not artist.get_visible():
                continue
            try:
                ab = artist.get_window_extent(renderer)
            except Exception:
                continue
            lefts.append(inv.transform((ab.x0, ab.y0))[0] * fig_w)
        if not lefts:
            return 0.55 * scale
        ax_left, _ = _axes_box_in(ax)
        return max(0.45 * scale, ax_left - min(lefts) + 0.08 * scale)

    def _cbar_right_extent_in(ax) -> float:
        rights = []
        for artist in list(ax.get_yticklabels()) + [ax.yaxis.label]:
            if artist is None or not artist.get_visible():
                continue
            try:
                ab = artist.get_window_extent(renderer)
            except Exception:
                continue
            rights.append(inv.transform((ab.x1, ab.y1))[0] * fig_w)
        _, ax_right = _axes_box_in(ax)
        if not rights:
            return 0.7 * scale
        return max(0.55 * scale, max(rights) - ax_right + 0.06 * scale)

    b_label_in = max(
        _label_left_extent_in(ax_b_bot),
        _label_left_extent_in(ax_b_top),
    )
    cbar_right_in = _cbar_right_extent_in(cbar_ax)
    right_margin = 0.08 * scale

    def _apply_x(b_xin: float, c_xin: float) -> None:
        ax_b_top.set_position(
            [b_xin / fig_w, b_top_y0 / fig_h, b_w_in / fig_w, new_sum_in / fig_h]
        )
        ax_b_bot.set_position(
            [b_xin / fig_w, b_bot_y0 / fig_h, b_w_in / fig_w, new_body_in / fig_h]
        )
        cbar_ax.set_position(
            [c_xin / fig_w, b_bot_y0 / fig_h, cbar_w_in / fig_w, a_h_in / fig_h]
        )

    # A–B gap fits Target+ticks; B–colorbar gap is tighter (colorbar closer to B).
    target_clearance_in = 0.45 * scale
    cbar_gap_in = 0.22 * scale
    for _ in range(4):
        gap_ab = b_label_in + target_clearance_in
        fixed = (
            label_left_in
            + a_w_in
            + gap_ab
            + b_w_in
            + cbar_gap_in
            + cbar_w_in
            + cbar_right_in
            + right_margin
        )
        if fixed > fig_w + 1e-6:
            fig_w = fixed + 0.02 * scale
            fig.set_size_inches(fig_w, fig_h, forward=True)
            ax_a.set_position(
                [label_left_in / fig_w, a_y0 / fig_h, a_w_in / fig_w, a_h_in / fig_h]
            )
        b_x = label_left_in + a_w_in + gap_ab
        c_x = b_x + b_w_in + cbar_gap_in
        _apply_x(b_x, c_x)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        b_label_in = max(
            _label_left_extent_in(ax_b_bot),
            _label_left_extent_in(ax_b_top),
            b_label_in,
        )
        cbar_right_in = max(cbar_right_in, _cbar_right_extent_in(cbar_ax))

    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    y_shared = max(
        inv.transform(ax.transAxes.transform((0.0, 1.0)))[1]
        for ax in (ax_a, ax_b_top)
    )
    dy = 8.0 * scale / 72.0 / fig.get_figheight()
    y_label = y_shared + dy
    for ax, letter, heading in (
        (ax_a, "A", "Pass rate of each test"),
        (ax_b_top, "B", "Pass all checks"),
    ):
        x0, _ = inv.transform(ax.transAxes.transform((0.0, 1.0)))
        x1, _ = inv.transform(ax.transAxes.transform((1.0, 1.0)))
        fig.text(
            x0,
            y_label,
            letter,
            fontsize=fs_panel,
            fontweight="bold",
            ha="right",
            va="bottom",
        )
        fig.text(
            0.5 * (x0 + x1),
            y_label,
            heading,
            fontsize=fs_panel,
            ha="center",
            va="bottom",
        )

    from . import style as _plot_style

    _fit_cell = _plot_style._fit_heatmap_cell_text
    _prevent = _plot_style._prevent_overlapping_text
    _old_web_w = _plot_style.WEB_HIGH_RES_WIDTH_PX
    _plot_style._fit_heatmap_cell_text = lambda *args, **kwargs: None
    _plot_style.WEB_HIGH_RES_WIDTH_PX = int(round(_old_web_w * scale))

    def _prevent_keep_cbar_fs(fig_arg):
        _prevent(fig_arg)
        style = getattr(fig_arg, "_affitox_style_cbar_ticks", None)
        if style is not None:
            style()
            return
        fs = getattr(fig_arg, "_affitox_cbar_tick_fs", None)
        family = getattr(fig_arg, "_affitox_cbar_tick_family", "DejaVu Sans")
        if fs is None:
            return
        for ax in fig_arg.axes:
            texts = [t.get_text().strip() for t in ax.get_yticklabels()]
            if texts and all(t in {"0", "20", "40", "60", "80", "100"} for t in texts if t):
                ax.tick_params(labelsize=fs)
                for t in ax.get_yticklabels():
                    t.set_fontsize(fs)
                    t.set_fontname(family)

    _plot_style._prevent_overlapping_text = _prevent_keep_cbar_fs
    try:
        save_figure(fig, out_dir / "pass_count_at_least_with_all_tests")
    finally:
        _plot_style._fit_heatmap_cell_text = _fit_cell
        _plot_style._prevent_overlapping_text = _prevent
        _plot_style.WEB_HIGH_RES_WIDTH_PX = _old_web_w


def plot_posebusters_at_least_all_tests_combined_var2(
    summary: pd.DataFrame,
    checks: pd.DataFrame,
    cfg: AnalysisConfig,
) -> None:
    """Refined layout: same structure as v1, adaptive heatmap labels, explicit colorbar ticks."""
    apply_style(cfg)
    cmap = sequential_cmap(cfg)
    out_dir = cfg.figures_dir / "posebusters"

    prepared = _prepare_at_least_all_tests_panels(summary, checks, cfg)
    if prepared is None:
        return
    n_tests, at_least, summary_df, body, n_body = prepared

    fig = plt.figure(figsize=(15.5, 6.8))
    outer = GridSpec(
        1,
        2,
        figure=fig,
        width_ratios=[1.08, 1.0],
        wspace=0.30,
        left=0.07,
        right=0.97,
        top=0.92,
        bottom=0.12,
    )

    ax_line = fig.add_subplot(outer[0, 0])
    _draw_pass_count_curve_on_ax(
        ax_line,
        at_least,
        cfg,
        n_tests,
        x_label=f"Minimum N tests passed (of {n_tests})",
        title="All targets: ligands passing at least N tests",
    )

    _draw_all_tests_heatmap_panel(
        fig,
        outer[0, 1],
        summary_df=summary_df,
        body=body,
        n_body=n_body,
        cmap=cmap,
        title="All tests (mean check pass rate)",
        adaptive_text=True,
    )

    save_figure(fig, out_dir / "pass_count_at_least_with_all_tests_var2")


def plot_posebusters_pass_count_curves(cfg: AnalysisConfig) -> None:
    """Line plots: exactly / at least N PoseBusters tests passed (All targets pooled)."""
    apply_style(cfg)
    out_dir = cfg.figures_dir / "posebusters"
    out_dir.mkdir(parents=True, exist_ok=True)

    n_tests, pooled = collect_pooled_pass_counts(cfg)
    if n_tests == 0 or not any(pooled.values()):
        return

    exactly: dict[str, pd.DataFrame] = {}
    at_least: dict[str, pd.DataFrame] = {}
    for method_id, counts in pooled.items():
        if not counts:
            continue
        exactly[method_id] = pass_count_distribution(counts, n_tests, mode="exactly")
        at_least[method_id] = pass_count_distribution(counts, n_tests, mode="at_least")

    tables_dir = cfg.tables_dir / "posebusters"
    tables_dir.mkdir(parents=True, exist_ok=True)
    for method_id in exactly:
        label = cfg.method_label(method_id)
        exactly[method_id].assign(method_id=method_id, method_label=label).to_csv(
            tables_dir / f"pass_count_exactly_{method_id}.csv", index=False
        )
        at_least[method_id].assign(method_id=method_id, method_label=label).to_csv(
            tables_dir / f"pass_count_at_least_{method_id}.csv", index=False
        )

    x_label = f"Number of PoseBusters tests passed (of {n_tests})"
    _plot_pass_count_curve(
        exactly,
        cfg=cfg,
        n_tests=n_tests,
        title="All targets: ligands passing exactly N tests",
        x_label=x_label,
        path_base=out_dir / "pass_count_exactly_n_tests",
    )
    _plot_pass_count_curve(
        at_least,
        cfg=cfg,
        n_tests=n_tests,
        title="All targets: ligands passing at least N tests",
        x_label=f"Minimum N tests passed (of {n_tests})",
        path_base=out_dir / "pass_count_at_least_n_tests",
    )
