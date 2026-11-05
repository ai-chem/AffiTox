"""Inferential statistics plots for pairwise Wilcoxon comparisons."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ..config import AnalysisConfig
from .style import FS_ANNOT, FS_AXIS, FS_PANEL, FS_TICK,  apply_style, save_figure, diverging_cmap, sequential_cmap


FAMILY_ORDER: tuple[str, ...] = (
    "Pearson r",
    "nEF active (top 10%)",
    "Spearman rho",
    "Kendall tau",
    "nEF inactive (top 10%)",
)
CORR_FAMILIES: tuple[str, ...] = (
    "Pearson r",
    "Spearman rho",
    "Kendall tau",
)
NEF_FAMILIES: tuple[str, ...] = (
    "nEF active (top 10%)",
    "nEF inactive (top 10%)",
)
PER_TARGET_FAMILY_ORDER: tuple[str, ...] = (
    "Pearson r (per-target)",
    "Spearman rho (per-target)",
    "nEF10 active (per-target)",
    "nEF10 inactive (per-target)",
)


def _method_order(tests: pd.DataFrame, cfg: AnalysisConfig) -> list[str]:
    ordered = [m for m in cfg.methods if m in set(tests["method_a"]) | set(tests["method_b"])]
    extras = sorted((set(tests["method_a"]) | set(tests["method_b"])) - set(ordered))
    return ordered + extras


def _family_matrix(
    tests: pd.DataFrame,
    family: str,
    methods: list[str],
    value_col: str,
) -> pd.DataFrame:
    mat = pd.DataFrame(np.nan, index=methods, columns=methods, dtype=float)
    np.fill_diagonal(mat.values, 0.0)

    sub = tests[tests["family"] == family]
    for _, row in sub.iterrows():
        a = row["method_a"]
        b = row["method_b"]
        v = float(row[value_col])
        mat.loc[a, b] = v
        mat.loc[b, a] = -v if value_col == "median_diff_a_minus_b" else v
    return mat


def _logp_matrix(pm: pd.DataFrame) -> pd.DataFrame:
    return -np.log10(pm.clip(lower=1e-12))


def _logp_vmax(logp_mats: dict[str, pd.DataFrame]) -> float:
    vals = [
        mat.values[np.isfinite(mat.values)]
        for mat in logp_mats.values()
        if mat is not None and mat.size
    ]
    if not vals or not any(v.size for v in vals):
        return 3.0
    return max(3.0, float(np.nanmax(np.concatenate(vals))))


def _annot_with_sig(effect: pd.DataFrame, pvals: pd.DataFrame, alpha: float = 0.05) -> np.ndarray:
    annot = np.empty(effect.shape, dtype=object)
    for i in range(effect.shape[0]):
        for j in range(effect.shape[1]):
            val = effect.iat[i, j]
            if pd.isna(val):
                annot[i, j] = ""
                continue
            if i == j:
                annot[i, j] = "0.00"
                continue
            star = "*" if pd.notna(pvals.iat[i, j]) and pvals.iat[i, j] < alpha else ""
            annot[i, j] = f"{val:.2f}{star}"
    return annot


def plot_inferential_heatmaps(tests: pd.DataFrame, cfg: AnalysisConfig) -> None:
    """Create SI-ready heatmaps of pairwise effect sizes and Holm-adjusted p-values."""
    if tests.empty:
        return

    apply_style(cfg)
    out_dir = cfg.figures_dir / "inferential"
    out_dir.mkdir(parents=True, exist_ok=True)

    methods = _method_order(tests, cfg)
    labels = {m: cfg.method_label(m) for m in methods}
    families = [f for f in FAMILY_ORDER if f in set(tests["family"])]
    if not families:
        families = sorted(set(tests["family"]))

    vmax = 0.0
    effect_mats: dict[str, pd.DataFrame] = {}
    pval_mats: dict[str, pd.DataFrame] = {}
    for family in families:
        em = _family_matrix(tests, family, methods, "median_diff_a_minus_b")
        pm = _family_matrix(tests, family, methods, "p_holm")
        effect_mats[family] = em
        pval_mats[family] = pm
        finite_vals = np.abs(em.values[np.isfinite(em.values)])
        if finite_vals.size:
            vmax = max(vmax, float(np.nanmax(finite_vals)))
    if vmax == 0.0:
        vmax = 0.5

    corr_families = [f for f in CORR_FAMILIES if f in set(tests["family"])]
    nef_families = [f for f in NEF_FAMILIES if f in set(tests["family"])]
    if len(corr_families) < 3:
        corr_families = [f for f in families if "nEF" not in f][:3]
    if len(nef_families) < 2:
        nef_families = [f for f in families if "nEF" in f][:2]

    # (1) Correlation families in one horizontal row, single scale after third panel.
    fig1 = plt.figure(figsize=(16.0, 5.4))
    gs1 = fig1.add_gridspec(1, 4, width_ratios=[1.0, 1.0, 1.0, 0.06], wspace=0.18)
    cmap_eff = diverging_cmap(cfg, name="inferential_effect")
    for idx, family in enumerate(corr_families[:3]):
        ax = fig1.add_subplot(gs1[0, idx])
        em = effect_mats[family].rename(index=labels, columns=labels)
        pm = pval_mats[family].rename(index=labels, columns=labels)
        cbar_ax = fig1.add_subplot(gs1[0, 3]) if idx == 2 else None
        sns.heatmap(
            em,
            ax=ax,
            cmap=cmap_eff,
            center=0.0,
            vmin=-vmax,
            vmax=vmax,
            annot=_annot_with_sig(em, pm),
            fmt="",
            linewidths=0.5,
            linecolor="white",
            cbar=(idx == 2),
            cbar_ax=cbar_ax,
            cbar_kws={"label": "Median paired difference (A - B)"},
            yticklabels=(idx == 0),
            xticklabels=True,
        )
        ax.set_title(family, fontsize=FS_PANEL)
        ax.tick_params(axis="x", rotation=35)
        plt.setp(ax.get_xticklabels(), ha="right")
        if idx == 0:
            ax.tick_params(axis="y", rotation=0)
            plt.setp(ax.get_yticklabels(), ha="right")
        else:
            ax.set_yticklabels([])
            ax.tick_params(axis="y", left=False, labelleft=False)
    fig1.suptitle("Pairwise Wilcoxon effects: Pearson / Spearman / Kendall", y=1.02)
    save_figure(fig1, out_dir / "heatmap_wilcoxon_effects")

    # (2) nEF families side-by-side; shared colorbar on the right of both panels.
    n_methods = max(len(methods), 1)
    cell_in = 1.15
    panel_in = n_methods * cell_in
    label_fs = 8
    panel_title_fs = FS_PANEL
    suptitle_fs = FS_PANEL
    annot_size = 7.5
    fig_nef = plt.figure(figsize=(2 * panel_in + 1.6, panel_in + 3.0))
    gs_nef = fig_nef.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 0.07], wspace=0.22)
    nef_axes: list = []
    cbar_ax = None
    for idx, family in enumerate(nef_families[:2]):
        ax = fig_nef.add_subplot(gs_nef[0, idx])
        nef_axes.append(ax)
        em = effect_mats[family].rename(index=labels, columns=labels)
        pm = pval_mats[family].rename(index=labels, columns=labels)
        if idx == 1:
            cbar_ax = fig_nef.add_subplot(gs_nef[0, 2])
        sns.heatmap(
            em,
            ax=ax,
            cmap=cmap_eff,
            center=0.0,
            vmin=-vmax,
            vmax=vmax,
            annot=_annot_with_sig(em, pm),
            fmt="",
            linewidths=0.5,
            linecolor="white",
            square=True,
            cbar=(idx == 1),
            cbar_ax=cbar_ax,
            cbar_kws={},
            yticklabels=(idx == 0),
            xticklabels=True,
            annot_kws={"size": annot_size, "ha": "center", "va": "center"},
        )
        ax.set_title(family, fontsize=panel_title_fs)
        ax.tick_params(axis="x", rotation=35, labelsize=label_fs)
        plt.setp(ax.get_xticklabels(), ha="right")
        if idx == 0:
            ax.tick_params(axis="y", rotation=0, labelsize=label_fs)
            plt.setp(ax.get_yticklabels(), ha="right")
        else:
            ax.set_yticklabels([])
            ax.tick_params(axis="y", left=False, labelleft=False)
        for text in ax.texts:
            val_txt = text.get_text().rstrip("*")
            try:
                val = float(val_txt)
            except ValueError:
                continue
            text.set_color("#ffffff" if abs(val) >= 0.35 else "#262626")

    if nef_axes and cbar_ax is not None:
        cbar_ax.tick_params(labelsize=label_fs)
        fig_nef.canvas.draw()
        renderer = fig_nef.canvas.get_renderer()
        hm_pos = nef_axes[0].get_position()
        cb_pos = cbar_ax.get_position()
        cbar_ax.set_position([cb_pos.x0, hm_pos.y0, cb_pos.width, hm_pos.height])
        title_bbox = nef_axes[0].title.get_window_extent(renderer).transformed(
            fig_nef.transFigure.inverted()
        )
        font_h_fig = suptitle_fs / 72.0 / fig_nef.get_figheight()
        fig_nef.suptitle(
            r"Pairwise Wilcoxon effects: nEF$_{10,\mathrm{active}}$ / nEF$_{10,\mathrm{inactive}}$",
            y=title_bbox.y1 + font_h_fig,
            va="bottom",
            fontsize=suptitle_fs,
        )
    save_figure(fig_nef, out_dir / "heatmap_wilcoxon_effects_nef10")

    cmap_p = sequential_cmap(cfg, name="inferential_logp")
    logp_label = r"$-\log_{10}(p_{\mathrm{Holm}})$"
    corr_logp_mats = {f: _logp_matrix(pval_mats[f]) for f in corr_families[:3]}
    nef_logp_mats = {f: _logp_matrix(pval_mats[f]) for f in nef_families[:2]}
    vmax_logp_corr = _logp_vmax(corr_logp_mats)
    vmax_logp_nef = _logp_vmax(nef_logp_mats)

    # (3) Correlation families: Holm-adjusted significance, one horizontal row.
    fig_logp_corr = plt.figure(figsize=(16.0, 5.4))
    gs_logp_corr = fig_logp_corr.add_gridspec(1, 4, width_ratios=[1.0, 1.0, 1.0, 0.06], wspace=0.18)
    for idx, family in enumerate(corr_families[:3]):
        ax = fig_logp_corr.add_subplot(gs_logp_corr[0, idx])
        logp = corr_logp_mats[family].rename(index=labels, columns=labels)
        cbar_ax = fig_logp_corr.add_subplot(gs_logp_corr[0, 3]) if idx == 2 else None
        sns.heatmap(
            logp,
            ax=ax,
            cmap=cmap_p,
            vmin=0.0,
            vmax=vmax_logp_corr,
            annot=True,
            fmt=".2f",
            linewidths=0.5,
            linecolor="white",
            cbar=(idx == 2),
            cbar_ax=cbar_ax,
            cbar_kws={"label": logp_label},
            yticklabels=(idx == 0),
            xticklabels=True,
        )
        ax.set_title(family, fontsize=FS_PANEL)
        ax.tick_params(axis="x", rotation=35)
        plt.setp(ax.get_xticklabels(), ha="right")
        if idx == 0:
            ax.tick_params(axis="y", rotation=0)
            plt.setp(ax.get_yticklabels(), ha="right")
        else:
            ax.set_yticklabels([])
            ax.tick_params(axis="y", left=False, labelleft=False)
    fig_logp_corr.suptitle(
        rf"Pairwise Wilcoxon significance ({logp_label}): Pearson / Spearman / Kendall",
        y=1.02,
    )
    save_figure(fig_logp_corr, out_dir / "heatmap_wilcoxon_logp_holm")

    # (4) nEF families: Holm-adjusted significance, side-by-side with shared colorbar.
    fig_logp_nef = plt.figure(figsize=(2 * panel_in + 1.6, panel_in + 3.0))
    gs_logp_nef = fig_logp_nef.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 0.07], wspace=0.22)
    nef_logp_axes: list = []
    logp_cbar_ax = None
    for idx, family in enumerate(nef_families[:2]):
        ax = fig_logp_nef.add_subplot(gs_logp_nef[0, idx])
        nef_logp_axes.append(ax)
        logp = nef_logp_mats[family].rename(index=labels, columns=labels)
        if idx == 1:
            logp_cbar_ax = fig_logp_nef.add_subplot(gs_logp_nef[0, 2])
        sns.heatmap(
            logp,
            ax=ax,
            cmap=cmap_p,
            vmin=0.0,
            vmax=vmax_logp_nef,
            annot=True,
            fmt=".2f",
            linewidths=0.5,
            linecolor="white",
            square=True,
            cbar=(idx == 1),
            cbar_ax=logp_cbar_ax,
            cbar_kws={"label": logp_label},
            yticklabels=(idx == 0),
            xticklabels=True,
            annot_kws={"size": annot_size, "ha": "center", "va": "center"},
        )
        ax.set_title(family, fontsize=panel_title_fs)
        ax.tick_params(axis="x", rotation=35, labelsize=label_fs)
        plt.setp(ax.get_xticklabels(), ha="right")
        if idx == 0:
            ax.tick_params(axis="y", rotation=0, labelsize=label_fs)
            plt.setp(ax.get_yticklabels(), ha="right")
        else:
            ax.set_yticklabels([])
            ax.tick_params(axis="y", left=False, labelleft=False)
        for text in ax.texts:
            val_txt = text.get_text()
            try:
                val = float(val_txt)
            except ValueError:
                continue
            text.set_color("#ffffff" if val >= 0.6 * vmax_logp_nef else "#262626")

    if nef_logp_axes and logp_cbar_ax is not None:
        logp_cbar_ax.tick_params(labelsize=label_fs)
        fig_logp_nef.canvas.draw()
        renderer = fig_logp_nef.canvas.get_renderer()
        hm_pos = nef_logp_axes[0].get_position()
        cb_pos = logp_cbar_ax.get_position()
        logp_cbar_ax.set_position([cb_pos.x0, hm_pos.y0, cb_pos.width, hm_pos.height])
        title_bbox = nef_logp_axes[0].title.get_window_extent(renderer).transformed(
            fig_logp_nef.transFigure.inverted()
        )
        font_h_fig = suptitle_fs / 72.0 / fig_logp_nef.get_figheight()
        fig_logp_nef.suptitle(
            rf"Pairwise Wilcoxon significance ({logp_label}): "
            r"nEF$_{10,\mathrm{active}}$ / nEF$_{10,\mathrm{inactive}}$",
            y=title_bbox.y1 + font_h_fig,
            va="bottom",
            fontsize=suptitle_fs,
        )
    save_figure(fig_logp_nef, out_dir / "heatmap_wilcoxon_logp_holm_nef10")


def plot_per_target_inferential_heatmaps(per_target_tests: pd.DataFrame, cfg: AnalysisConfig) -> None:
    """Create per-target pairwise effect and significance heatmaps."""
    if per_target_tests.empty:
        return

    apply_style(cfg)
    out_dir = cfg.figures_dir / "inferential"
    out_dir.mkdir(parents=True, exist_ok=True)

    targets = [t.lower() for t in cfg.targets]
    fams = [f for f in PER_TARGET_FAMILY_ORDER if f in set(per_target_tests["family"])]
    if not fams:
        fams = sorted(per_target_tests["family"].unique())
    all_comparisons = sorted(per_target_tests["comparison"].unique())

    vmax = float(np.nanmax(np.abs(per_target_tests["effect_obs"].to_numpy(dtype=float))))
    if not np.isfinite(vmax) or vmax == 0.0:
        vmax = 0.5

    cmap_eff = diverging_cmap(cfg, name="per_target_effect")
    cmap_logp = sequential_cmap(cfg, name="per_target_logp")
    vmax_logp = float(
        np.nanmax(
            -np.log10(per_target_tests["perm_p_bh_target"].clip(lower=1e-12).to_numpy(dtype=float))
        )
    )
    if not np.isfinite(vmax_logp):
        vmax_logp = 3.0
    vmax_logp = max(3.0, vmax_logp)

    for target in targets:
        sub_t = per_target_tests[per_target_tests["target"] == target].copy()
        if sub_t.empty:
            continue

        eff = sub_t.pivot(index="family", columns="comparison", values="effect_obs")
        pbh = sub_t.pivot(index="family", columns="comparison", values="perm_p_bh_target")
        eff = eff.reindex(index=fams, columns=all_comparisons)
        pbh = pbh.reindex(index=fams, columns=all_comparisons)
        logp = -np.log10(pbh.clip(lower=1e-12))

        fig_e, ax_e = plt.subplots(1, 1, figsize=(max(9.0, 0.40 * len(all_comparisons)), 3.8))
        sns.heatmap(
            eff,
            ax=ax_e,
            cmap=cmap_eff,
            center=0.0,
            vmin=-vmax,
            vmax=vmax,
            annot=True,
            fmt=".2f",
            linewidths=0.4,
            linecolor="white",
            cbar=True,
            cbar_kws={"label": "Effect size (A - B)"},
        )
        ax_e.set_title(
            f"{cfg.target_label(target)}: per-target pairwise effects",
            fontsize=FS_AXIS,
        )
        ax_e.set_xlabel("Method pair")
        ax_e.set_ylabel("Family")
        ax_e.tick_params(axis="x", rotation=45)
        ax_e.tick_params(axis="y", rotation=0)
        save_figure(fig_e, out_dir / f"per_target_{target}_wilcoxon_effects")

        fig_p, ax_p = plt.subplots(1, 1, figsize=(max(9.0, 0.40 * len(all_comparisons)), 3.8))
        sns.heatmap(
            logp,
            ax=ax_p,
            cmap=cmap_logp,
            vmin=0.0,
            vmax=vmax_logp,
            annot=True,
            fmt=".2f",
            linewidths=0.4,
            linecolor="white",
            cbar=True,
            cbar_kws={"label": r"$-\log_{10}(p_{\mathrm{BH,target}})$"},
        )
        ax_p.set_title(
            f"{cfg.target_label(target)}: per-target pairwise significance",
            fontsize=FS_AXIS,
        )
        ax_p.set_xlabel("Method pair")
        ax_p.set_ylabel("Family")
        ax_p.tick_params(axis="x", rotation=45)
        ax_p.tick_params(axis="y", rotation=0)
        save_figure(fig_p, out_dir / f"per_target_{target}_wilcoxon_logp")
