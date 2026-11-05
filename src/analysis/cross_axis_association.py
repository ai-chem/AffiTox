#!/usr/bin/env python3
"""Part A: target×method cross-axis association of pose validity vs affinity/screening."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from analysis.constants import PRIMARY_METRICS
from analysis.data import load_merged, resolve_posebusters_csv
from analysis.inferential import benjamini_hochberg
from analysis.robustness_common import (
    ALL_METHODS,
    DEFAULT_N_BOOT,
    DEFAULT_N_PERM,
    DEFAULT_SEED,
    FIGURES,
    MERGED_DIR,
    POSE_DIR,
    POSE_METHODS,
    ROBUST_FIGURES,
    ROBUST_TABLES,
    TABLES,
    TARGETS,
    RobustnessConfig,
    cluster_bootstrap_spearman,
    collapse_posebusters,
    inventory_row,
    kendall_safe,
    load_nodubl_index,
    method_label,
    residualize_two_way,
    restricted_permutation_p,
    spearman_safe,
    target_label,
)
from analysis.plots.style import DEFAULT_METHOD_COLORS, save_figure

POSE_ENDPOINTS = (
    "strict_pass_end_to_end",
    "strict_pass_conditional",
    "mean_fraction_checks_passed",
)
AFFINITY_ENDPOINTS = (
    "pearson_r",
    "spearman_rho",
    "kendall_tau",
    "nEF10_active",
    "nEF10_low",
)
PRIMARY_PAIR = ("strict_pass_end_to_end", "pearson_r")
RANK_AXES = ("pearson_r", "strict_pass_end_to_end", "nEF10_active", "nEF10_low")


def _load_published_affinity() -> pd.DataFrame:
    corr = pd.read_csv(TABLES / "correlations" / "correlations_with_ci_perm.csv")
    ef = pd.read_csv(TABLES / "enrichment" / "ef_summary_all_proteins.csv")
    corr = corr.loc[corr["metric"].isin(PRIMARY_METRICS.values())].copy()
    corr = corr.loc[corr["method_id"].isin(POSE_METHODS)].copy()
    ef = ef.loc[ef["method_id"].isin(POSE_METHODS)].copy()
    if corr.duplicated(["target", "method_id"]).any():
        raise AssertionError("Duplicate target×method rows in correlations_with_ci_perm.csv")
    if ef.duplicated(["target", "method_id"]).any():
        raise AssertionError("Duplicate target×method rows in ef_summary_all_proteins.csv")
    keep_c = [
        "target",
        "method_id",
        "n_points",
        "pearson_r",
        "spearman_rho",
        "kendall_tau",
    ]
    keep_e = [
        "target",
        "method_id",
        "nEF10_active",
        "nEF10_low",
        "N_10_active",
        "N_10_low",
        "A_10_active",
        "A_10_low",
    ]
    out = corr[keep_c].merge(ef[keep_e], on=["target", "method_id"], how="inner", validate="1:1")
    out["target"] = out["target"].astype(str).str.lower()
    return out


def _pose_panel_for_target(target: str) -> dict[str, dict]:
    merged = load_merged(MERGED_DIR, target)
    # Repeated Ki rows share a canonical SMILES and are not extra ligands.
    # The end-to-end pose denominator is one compound, matching the correlation unit.
    n_input = int(merged["canonical_smiles"].nunique())
    nodubl = load_nodubl_index(target)
    per_method: dict[str, dict] = {}
    for method_id in POSE_METHODS:
        csv_path = resolve_posebusters_csv(POSE_DIR, POSE_DIR, target, method_id)
        rec = {
            "n_input_ligands": n_input,
            "n_generated_pose_rows": 0,
            "n_unique_pb_tokens": 0,
            "n_unmatched_pb_ligands": 0,
            "n_generated_poses": 0,
            "n_evaluable_poses": 0,
            "n_pass_all_posebusters": 0,
            "mean_fraction_checks_passed": np.nan,
            "pose_generation_coverage": np.nan,
            "strict_pass_conditional": np.nan,
            "strict_pass_end_to_end": np.nan,
            "n_pose_end_to_end": n_input,
            "n_pose_conditional": 0,
            "pb_csv": str(csv_path) if csv_path else "",
        }
        if csv_path is None or not csv_path.exists():
            rec["pose_generation_coverage"] = 0.0
            rec["strict_pass_end_to_end"] = 0.0
            rec["strict_pass_conditional"] = np.nan
            per_method[method_id] = rec
            continue
        raw = pd.read_csv(csv_path)
        rec["n_generated_pose_rows"] = int(len(raw))
        collapsed = collapse_posebusters(raw, nodubl, merged)
        if collapsed.empty:
            per_method[method_id] = rec
            continue
        rec["n_unique_pb_tokens"] = int(collapsed["token"].nunique())
        rec["n_unmatched_pb_ligands"] = int((~collapsed["matched_to_merged"]).sum())
        matched = collapsed.loc[collapsed["matched_to_merged"]].copy()
        if matched.duplicated("canonical_smiles").any():
            matched = matched.drop_duplicates("canonical_smiles", keep="first")
        rec["n_generated_poses"] = int(len(matched))
        evaluable = matched.loc[matched["evaluable"]]
        rec["n_evaluable_poses"] = int(len(evaluable))
        rec["n_pass_all_posebusters"] = int(evaluable["pass_all"].sum()) if len(evaluable) else 0
        rec["mean_fraction_checks_passed"] = (
            float(evaluable["frac_checks_passed"].mean()) if len(evaluable) else np.nan
        )
        rec["n_pose_conditional"] = rec["n_evaluable_poses"]
        rec["pose_generation_coverage"] = rec["n_evaluable_poses"] / n_input if n_input else np.nan
        rec["strict_pass_conditional"] = (
            rec["n_pass_all_posebusters"] / rec["n_evaluable_poses"]
            if rec["n_evaluable_poses"]
            else np.nan
        )
        rec["strict_pass_end_to_end"] = rec["n_pass_all_posebusters"] / n_input if n_input else np.nan
        per_method[method_id] = rec
    return per_method


def build_long_table() -> pd.DataFrame:
    published = _load_published_affinity()
    rows = []
    for target in TARGETS:
        pose = _pose_panel_for_target(target)
        merged = load_merged(MERGED_DIR, target)
        n_input = int(merged["canonical_smiles"].nunique())
        for method_id in POSE_METHODS:
            aff = published.loc[
                (published["target"] == target) & (published["method_id"] == method_id)
            ]
            if aff.empty:
                raise AssertionError(f"Missing published affinity row for {target} {method_id}")
            a = aff.iloc[0]
            p = pose[method_id]
            if p["n_input_ligands"] != n_input:
                raise AssertionError("n_input mismatch")
            rows.append(
                {
                    "target": target,
                    "target_label": target_label(target),
                    "method_id": method_id,
                    "method": method_label(method_id),
                    "n_input_ligands": p["n_input_ligands"],
                    "n_generated_pose_rows": p["n_generated_pose_rows"],
                    "n_unique_pb_tokens": p["n_unique_pb_tokens"],
                    "n_unmatched_pb_ligands": p["n_unmatched_pb_ligands"],
                    "n_generated_poses": p["n_generated_poses"],
                    "n_evaluable_poses": p["n_evaluable_poses"],
                    "n_pass_all_posebusters": p["n_pass_all_posebusters"],
                    "pose_generation_coverage": p["pose_generation_coverage"],
                    "strict_pass_conditional": p["strict_pass_conditional"],
                    "strict_pass_end_to_end": p["strict_pass_end_to_end"],
                    "mean_fraction_checks_passed": p["mean_fraction_checks_passed"],
                    "pearson_r": float(a["pearson_r"]),
                    "spearman_rho": float(a["spearman_rho"]),
                    "kendall_tau": float(a["kendall_tau"]),
                    "nEF10_active": float(a["nEF10_active"]),
                    "nEF10_low": float(a["nEF10_low"]),
                    "n_pearson": int(a["n_points"]),
                    "n_spearman": int(a["n_points"]),
                    "n_kendall": int(a["n_points"]),
                    "n_nef10_active": int(a["N_10_active"]),
                    "n_nef10_low": int(a["N_10_low"]),
                    "n_pose_end_to_end": p["n_pose_end_to_end"],
                    "n_pose_conditional": p["n_pose_conditional"],
                    "n_mean_frac_checks": p["n_evaluable_poses"],
                    "A_10_active": int(a["A_10_active"]),
                    "A_10_low": int(a["A_10_low"]),
                }
            )
    long = pd.DataFrame(rows)
    if long.duplicated(["target", "method_id"]).any():
        raise AssertionError("Duplicate target×method rows in cross-axis long table")
    expected = 16 * 4
    if len(long) != expected:
        raise AssertionError(f"Expected {expected} target×method rows, got {len(long)}")
    return long


def _association_row(
    long: pd.DataFrame,
    pose_col: str,
    aff_col: str,
    cfg: RobustnessConfig,
    *,
    family: str,
) -> dict:
    x = long[pose_col].to_numpy(dtype=float)
    y = long[aff_col].to_numpy(dtype=float)
    targets = long["target"].to_numpy()
    methods = long["method_id"].to_numpy()
    mask = np.isfinite(x) & np.isfinite(y)
    rho_raw = spearman_safe(x, y)
    lo_raw, hi_raw, mean_raw = cluster_bootstrap_spearman(
        x, y, targets, cfg.n_boot, cfg.seed
    )
    x_adj = residualize_two_way(x, targets, methods)
    y_adj = residualize_two_way(y, targets, methods)
    rho_adj, p_adj = restricted_permutation_p(
        x, y, targets, methods, cfg.n_perm, cfg.seed + 1, adjusted=True
    )
    lo_adj, hi_adj, mean_adj = cluster_bootstrap_spearman(
        x_adj, y_adj, targets, cfg.n_boot, cfg.seed + 2
    )
    _, p_raw = restricted_permutation_p(
        x, y, targets, methods, cfg.n_perm, cfg.seed + 3, adjusted=False
    )
    return {
        "family": family,
        "pose_endpoint": pose_col,
        "affinity_endpoint": aff_col,
        "pose_axis": "docking_power",
        "affinity_axis": {
            "pearson_r": "scoring_power",
            "spearman_rho": "ranking_power",
            "kendall_tau": "ranking_power",
            "nEF10_active": "screening_power",
            "nEF10_low": "screening_power",
        }[aff_col],
        "n_observations": int(mask.sum()),
        "n_targets": int(pd.unique(targets[mask]).size),
        "n_methods": int(pd.unique(methods[mask]).size),
        "spearman_raw": rho_raw,
        "spearman_raw_ci_low": lo_raw,
        "spearman_raw_ci_high": hi_raw,
        "spearman_raw_boot_mean": mean_raw,
        "spearman_two_way_adjusted": rho_adj,
        "spearman_adj_ci_low": lo_adj,
        "spearman_adj_ci_high": hi_adj,
        "spearman_adj_boot_mean": mean_adj,
        "p_perm_raw_restricted": p_raw,
        "p_perm_adjusted_restricted": p_adj,
        "n_boot": cfg.n_boot,
        "n_perm": cfg.n_perm,
        "seed": cfg.seed,
        "is_primary": int((pose_col, aff_col) == PRIMARY_PAIR),
    }


def association_table(long: pd.DataFrame, cfg: RobustnessConfig) -> pd.DataFrame:
    rows = []
    rows.append(_association_row(long, PRIMARY_PAIR[0], PRIMARY_PAIR[1], cfg, family="primary"))
    for pose_col in POSE_ENDPOINTS:
        for aff_col in AFFINITY_ENDPOINTS:
            if (pose_col, aff_col) == PRIMARY_PAIR:
                continue
            rows.append(
                _association_row(long, pose_col, aff_col, cfg, family="secondary")
            )
    out = pd.DataFrame(rows)
    sec = out["family"] == "secondary"
    out["p_perm_adjusted_bh"] = np.nan
    if sec.any():
        out.loc[sec, "p_perm_adjusted_bh"] = benjamini_hochberg(
            out.loc[sec, "p_perm_adjusted_restricted"].to_numpy(dtype=float)
        )
    return out


def rank_agreement(long: pd.DataFrame) -> pd.DataFrame:
    rows = []
    disagreements = []
    for target, sub in long.groupby("target"):
        ranks = {}
        tops = {}
        for axis in RANK_AXES:
            vals = sub[axis].to_numpy(dtype=float)
            # Higher is better on every listed axis. Average ranks for ties.
            rnk = rankdata(-vals, method="average")
            ranks[axis] = dict(zip(sub["method_id"], rnk))
            best = float(np.nanmin(rnk))
            method_ids = sub["method_id"].to_numpy()
            tops[axis] = set(method_ids[np.isclose(rnk, best)])
        for i, a in enumerate(RANK_AXES):
            for b in RANK_AXES[i + 1 :]:
                ra = np.array([ranks[a][m] for m in sub["method_id"]])
                rb = np.array([ranks[b][m] for m in sub["method_id"]])
                tau = kendall_safe(ra, rb)
                mad = float(np.mean(np.abs(ra - rb)))
                same_top = int(len(tops[a] & tops[b]) > 0 and tops[a] == tops[b])
                rows.append(
                    {
                        "target": target,
                        "target_label": target_label(target),
                        "axis_a": a,
                        "axis_b": b,
                        "same_top_method": same_top,
                        "top_a": "|".join(sorted(tops[a])),
                        "top_b": "|".join(sorted(tops[b])),
                        "mean_abs_rank_diff": mad,
                        "kendall_tau_ranks": tau,
                        "n_methods": int(len(sub)),
                    }
                )
                if not same_top:
                    disagreements.append(
                        {
                            "target": target,
                            "axis_a": a,
                            "axis_b": b,
                            "top_a": "|".join(sorted(tops[a])),
                            "top_b": "|".join(sorted(tops[b])),
                        }
                    )
    detail = pd.DataFrame(rows)
    summary_rows = []
    for (a, b), g in detail.groupby(["axis_a", "axis_b"]):
        summary_rows.append(
            {
                "axis_a": a,
                "axis_b": b,
                "n_targets": int(g["target"].nunique()),
                "proportion_same_top_method": float(g["same_top_method"].mean()),
                "mean_of_mean_abs_rank_diff": float(g["mean_abs_rank_diff"].mean()),
                "median_kendall_tau_ranks": float(g["kendall_tau_ranks"].median()),
                "n_targets_discordant_top": int((g["same_top_method"] == 0).sum()),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary["block"] = "pair_summary"
    detail["block"] = "per_target"
    return pd.concat([summary, detail], ignore_index=True, sort=False)


def leave_one_target_out(long: pd.DataFrame, cfg: RobustnessConfig) -> pd.DataFrame:
    rows = []
    pose_col, aff_col = PRIMARY_PAIR
    for held in TARGETS:
        sub = long.loc[long["target"] != held]
        x = sub[pose_col].to_numpy(dtype=float)
        y = sub[aff_col].to_numpy(dtype=float)
        t = sub["target"].to_numpy()
        m = sub["method_id"].to_numpy()
        rho_raw = spearman_safe(x, y)
        x_adj = residualize_two_way(x, t, m)
        y_adj = residualize_two_way(y, t, m)
        rho_adj = spearman_safe(x_adj, y_adj)
        rows.append(
            {
                "held_out_target": held,
                "held_out_label": target_label(held),
                "n_observations": int(len(sub)),
                "n_targets": int(sub["target"].nunique()),
                "spearman_raw": rho_raw,
                "spearman_two_way_adjusted": rho_adj,
            }
        )
    return pd.DataFrame(rows)


def plot_association_matrix(assoc: pd.DataFrame, path_base: Path) -> None:
    mat = assoc.pivot(index="pose_endpoint", columns="affinity_endpoint", values="spearman_two_way_adjusted")
    mat = mat.reindex(index=list(POSE_ENDPOINTS), columns=list(AFFINITY_ENDPOINTS))
    fig, ax = plt.subplots(figsize=(8.2, 3.6))
    im = ax.imshow(mat.to_numpy(dtype=float), cmap="coolwarm", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=30, ha="right")
    ax.set_yticks(range(len(mat.index)))
    ax.set_yticklabels(mat.index)
    ax.set_title("Two-way-adjusted Spearman (pose vs affinity/screening endpoints)")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iat[i, j]
            if pd.notna(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Spearman ρ (adjusted)")
    fig.tight_layout()
    save_figure(fig, path_base)


def plot_pose_vs_pearson(long: pd.DataFrame, path_base: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    for method_id, sub in long.groupby("method_id"):
        ax.scatter(
            sub["strict_pass_end_to_end"],
            sub["pearson_r"],
            label=method_label(method_id),
            color=DEFAULT_METHOD_COLORS.get(method_id, "#555555"),
            s=42,
            edgecolors="k",
            linewidths=0.3,
        )
    ax.set_xlabel("Strict PoseBusters pass (end-to-end)")
    ax.set_ylabel("Within-target Pearson r (scoring power)")
    ax.set_title("Pose validity vs scoring power (16 targets × 4 pose methods)")
    ax.legend(frameon=False, fontsize=8)
    ax.axhline(0.0, color="#999999", lw=0.6)
    fig.tight_layout()
    save_figure(fig, path_base)


def write_cross_axis_report(long: pd.DataFrame, assoc: pd.DataFrame, ranks: pd.DataFrame, loto: pd.DataFrame) -> None:
    prim = assoc.loc[assoc["is_primary"] == 1].iloc[0]
    same_top = ranks.loc[
        (ranks["block"] == "pair_summary")
        & (ranks["axis_a"] == "pearson_r")
        & (ranks["axis_b"] == "strict_pass_end_to_end")
    ]
    prop = float(same_top["proportion_same_top_method"].iloc[0]) if len(same_top) else np.nan
    text = f"""# Cross-axis association analysis

Unit of analysis: one target × method observation.
Primary panel: 16 targets × 4 pose-generating methods (Boltz-2, DynamicBind, GNINA 1.3, QVina2).
PLAPT is excluded because it does not generate poses.

Pose denominators are ligand-level after collapsing duplicate PoseBusters rows and joining
to the merged AffiTox panel. Missing pose generation is **not** treated as a PoseBusters
failure in `strict_pass_conditional`. It is counted as unsuccessful only in
`strict_pass_end_to_end`.

## Terminology

These results describe **distinct and potentially nonredundant evaluation axes**
and **partially associated but sometimes discordant endpoints**.
They do **not** establish statistical independence of docking power and scoring power.

## Primary prespecified association

`strict_pass_end_to_end` vs Pearson r (scoring power):

- raw Spearman ρ = {prim['spearman_raw']:.3f} (95% target-cluster bootstrap CI {prim['spearman_raw_ci_low']:.3f} to {prim['spearman_raw_ci_high']:.3f})
- two-way-adjusted Spearman ρ = {prim['spearman_two_way_adjusted']:.3f} (CI {prim['spearman_adj_ci_low']:.3f} to {prim['spearman_adj_ci_high']:.3f})
- restricted permutation p (adjusted) = {prim['p_perm_adjusted_restricted']:.4g}
- n = {int(prim['n_observations'])} target–method observations ({int(prim['n_targets'])} targets, {int(prim['n_methods'])} methods)
- B_boot = {int(prim['n_boot'])}, B_perm = {int(prim['n_perm'])}, seed = {int(prim['seed'])}

The secondary family (remaining pose × affinity/screening pairs) uses Benjamini–Hochberg
adjustment on restricted-permutation p-values and is not used as a confirmatory claim.

## Rank nonredundancy

Proportion of targets with the same unique top method on Pearson r and strict end-to-end
validity: {prop:.3f}. Target-specific disagreements are in `cross_axis_rank_agreement.csv`.

## Leave-one-target-out

Primary-association sensitivity after dropping each target is in
`cross_axis_leave_one_target_out.csv` (n={len(loto)}).

## What this does not support

- Independence of pose validity and affinity estimation.
- That a high PoseBusters pass rate implies informative pKi ranking, or the converse.
"""
    (ROBUST_TABLES / "cross_axis_analysis.md").write_text(text)


def write_n_audit(long: pd.DataFrame) -> Path:
    corr = pd.read_csv(TABLES / "correlations" / "correlations_with_ci_perm.csv")
    rows = []
    for _, r in long.iterrows():
        published = corr.loc[
            (corr["target"] == r["target"]) & (corr["method_id"] == r["method_id"])
        ]
        merged = load_merged(MERGED_DIR, r["target"])
        rows.append(
            {
                "target": r["target"],
                "method_id": r["method_id"],
                "n_input_ligands": r["n_input_ligands"],
                "n_merged_rows": int(len(merged)),
                "n_unique_smiles": int(merged["canonical_smiles"].nunique()),
                "n_pearson": r["n_pearson"],
                "n_points_published": int(published["n_points"].iloc[0]) if len(published) else np.nan,
                "n_input_matches_unique_smiles": int(
                    r["n_input_ligands"] == merged["canonical_smiles"].nunique()
                ),
                "n_pearson_matches_published": int(
                    len(published) == 1 and int(published["n_points"].iloc[0]) == int(r["n_pearson"])
                ),
                "n_evaluable_le_input": int(r["n_evaluable_poses"] <= r["n_input_ligands"]),
                "end_to_end_uses_input_denom": int(
                    np.isclose(
                        r["strict_pass_end_to_end"] * r["n_input_ligands"],
                        r["n_pass_all_posebusters"],
                    )
                ),
            }
        )
    audit = pd.DataFrame(rows)
    path = ROBUST_TABLES / "cross_axis_n_audit.csv"
    audit.to_csv(path, index=False)
    if not bool(audit["n_input_matches_unique_smiles"].all()):
        raise AssertionError("n_input_ligands does not match unique canonical SMILES")
    if not bool(audit["n_pearson_matches_published"].all()):
        raise AssertionError("n_pearson does not match published n_points")
    return path


def run_cross_axis(cfg: RobustnessConfig) -> dict[str, Path]:
    cfg.ensure_dirs()
    long = build_long_table()
    assoc = association_table(long, cfg)
    ranks = rank_agreement(long)
    loto = leave_one_target_out(long, cfg)
    long_path = ROBUST_TABLES / "cross_axis_long.csv"
    assoc_path = ROBUST_TABLES / "cross_axis_associations.csv"
    rank_path = ROBUST_TABLES / "cross_axis_rank_agreement.csv"
    loto_path = ROBUST_TABLES / "cross_axis_leave_one_target_out.csv"
    long.to_csv(long_path, index=False)
    assoc.to_csv(assoc_path, index=False)
    ranks.to_csv(rank_path, index=False)
    loto.to_csv(loto_path, index=False)
    write_n_audit(long)
    plot_association_matrix(assoc, ROBUST_FIGURES / "cross_axis_association_matrix")
    plot_pose_vs_pearson(long, ROBUST_FIGURES / "pose_vs_pearson_scatter")
    write_cross_axis_report(long, assoc, ranks, loto)
    return {
        "long": long_path,
        "associations": assoc_path,
        "ranks": rank_path,
        "loto": loto_path,
    }


if __name__ == "__main__":
    run_cross_axis(RobustnessConfig())
