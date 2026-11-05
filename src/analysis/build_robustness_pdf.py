#!/usr/bin/env python3
"""Build a standalone PDF report of the AffiTox robustness analysis (code, methods, results)."""

from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "analysis" / "tables" / "robustness"
FIGDIR = ROOT / "analysis" / "figures" / "robustness"
SRC = ROOT / "src" / "analysis"
OUT = TABLES / "AffiTox_robustness_report.pdf"

A4 = (8.27, 11.69)
MONO = "DejaVu Sans Mono"
SANS = "DejaVu Sans"

CODE_FILES = [
    SRC / "robustness_common.py",
    SRC / "cross_axis_association.py",
    SRC / "training_exposure_extension.py",
    SRC / "temporal_validation.py",
    SRC / "screening_balanced_sensitivity.py",
    SRC / "run_robustness.py",
    ROOT / "tests" / "test_robustness.py",
]


def _page(pdf: PdfPages, title: str, lines: list[str], *, fontsize: float = 9.0) -> None:
    fig = plt.figure(figsize=A4)
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0.07, 0.05, 0.86, 0.90])
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.0, 0.98, title, fontsize=13, fontweight="bold", fontfamily=SANS, va="top")
    y = 0.93
    for line in lines:
        wrapped = textwrap.wrap(line, width=98) if line else [""]
        for w in wrapped:
            if y < 0.04:
                pdf.savefig(fig, dpi=150)
                plt.close(fig)
                fig = plt.figure(figsize=A4)
                fig.patch.set_facecolor("white")
                ax = fig.add_axes([0.07, 0.05, 0.86, 0.90])
                ax.axis("off")
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                ax.text(0.0, 0.98, title + " (continued)", fontsize=12, fontweight="bold", fontfamily=SANS, va="top")
                y = 0.93
            ax.text(0.0, y, w if w else " ", fontsize=fontsize, fontfamily=SANS, va="top")
            y -= 0.018 + (0.004 if fontsize >= 10 else 0.0)
    pdf.savefig(fig, dpi=150)
    plt.close(fig)


def _table_page(pdf: PdfPages, title: str, df: pd.DataFrame, col_width: float | None = None) -> None:
    fig = plt.figure(figsize=A4)
    fig.patch.set_facecolor("white")
    fig.suptitle(title, fontsize=11, fontfamily=SANS, fontweight="bold", y=0.97)
    ax = fig.add_axes([0.04, 0.04, 0.92, 0.90])
    ax.axis("off")
    show = df.copy()
    for c in show.columns:
        if pd.api.types.is_float_dtype(show[c]):
            show[c] = show[c].map(lambda v: "" if pd.isna(v) else f"{v:.3f}")
        else:
            show[c] = show[c].astype(str).str.slice(0, 42)
    tbl = ax.table(
        cellText=show.values.tolist(),
        colLabels=list(show.columns),
        loc="upper center",
        cellLoc="left",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(6.2)
    tbl.scale(1.0, 1.25)
    pdf.savefig(fig, dpi=150)
    plt.close(fig)


def _code_pages(pdf: PdfPages, path: Path) -> None:
    text = path.read_text(errors="replace").splitlines()
    rel = str(path.relative_to(ROOT))
    chunk = 58
    for i in range(0, len(text), chunk):
        block = text[i : i + chunk]
        fig = plt.figure(figsize=A4)
        fig.patch.set_facecolor("white")
        ax = fig.add_axes([0.04, 0.03, 0.94, 0.93])
        ax.axis("off")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        start = i + 1
        end = i + len(block)
        ax.text(
            0.0,
            0.99,
            f"Code: {rel}  lines {start}–{end} / {len(text)}",
            fontsize=8,
            fontweight="bold",
            fontfamily=SANS,
            va="top",
        )
        body = "\n".join(f"{start + j:4d}  {line[:108]}" for j, line in enumerate(block))
        ax.text(0.0, 0.95, body, fontsize=5.4, fontfamily=MONO, va="top")
        pdf.savefig(fig, dpi=140)
        plt.close(fig)


def _figures(pdf: PdfPages) -> None:
    long = pd.read_csv(TABLES / "cross_axis_long.csv")
    assoc = pd.read_csv(TABLES / "cross_axis_associations.csv")

    fig, ax = plt.subplots(figsize=A4)
    fig.subplots_adjust(top=0.72, bottom=0.38)
    ax.set_title("Pose validity vs scoring power (16 targets × 4 pose methods)")
    colors = {
        "boltz2": "#E5885F",
        "dynamicbind": "#FDBBAA",
        "gnina": "#AFC6F2",
        "qvina": "#6584E1",
    }
    for mid, sub in long.groupby("method_id"):
        ax.scatter(
            sub["strict_pass_end_to_end"],
            sub["pearson_r"],
            label=sub["method"].iloc[0],
            color=colors.get(mid, "#555"),
            s=36,
            edgecolors="k",
            linewidths=0.3,
        )
    ax.set_xlabel("Strict PoseBusters pass (end-to-end)")
    ax.set_ylabel("Within-target Pearson r")
    ax.axhline(0, color="#999", lw=0.6)
    ax.legend(frameon=False, fontsize=8)
    fig.text(
        0.08,
        0.12,
        "Each point is one target×method. Missing poses are unsuccessful only in the end-to-end rate.",
        fontsize=8,
        fontfamily=SANS,
    )
    pdf.savefig(fig, dpi=150)
    plt.close(fig)

    mat = assoc.pivot(index="pose_endpoint", columns="affinity_endpoint", values="spearman_two_way_adjusted")
    order_p = ["strict_pass_end_to_end", "strict_pass_conditional", "mean_fraction_checks_passed"]
    order_a = ["pearson_r", "spearman_rho", "kendall_tau", "nEF10_active", "nEF10_low"]
    mat = mat.reindex(index=order_p, columns=order_a)
    fig, ax = plt.subplots(figsize=A4)
    fig.subplots_adjust(top=0.72, bottom=0.32)
    im = ax.imshow(mat.to_numpy(dtype=float), cmap="coolwarm", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=25, ha="right", fontsize=8)
    ax.set_yticks(range(len(mat.index)))
    ax.set_yticklabels(mat.index, fontsize=8)
    ax.set_title("Two-way-adjusted Spearman (pose vs affinity/screening)")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iat[i, j]
            if pd.notna(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03, label="adjusted Spearman")
    pdf.savefig(fig, dpi=150)
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    assoc = pd.read_csv(TABLES / "cross_axis_associations.csv")
    long = pd.read_csv(TABLES / "cross_axis_long.csv")
    ranks = pd.read_csv(TABLES / "cross_axis_rank_agreement.csv")
    loto = pd.read_csv(TABLES / "cross_axis_leave_one_target_out.csv")
    ev = pd.read_csv(TABLES / "training_exposure_evidence.csv")
    bym = pd.read_csv(TABLES / "training_exposure_by_method.csv")
    exp = pd.read_csv(TABLES / "exposure_performance_sensitivity.csv")
    delta = pd.read_csv(TABLES / "temporal_within_target_delta_r.csv")
    feas = pd.read_csv(TABLES / "external_validation_feasibility.csv")
    bal = pd.read_csv(TABLES / "balanced_class_nef_sensitivity.csv")
    prim = assoc.loc[assoc["is_primary"] == 1].iloc[0]
    pair = ranks.loc[
        (ranks["block"] == "pair_summary")
        & (ranks["axis_a"] == "pearson_r")
        & (ranks["axis_b"] == "strict_pass_end_to_end")
    ].iloc[0]

    with PdfPages(OUT) as pdf:
        _page(
            pdf,
            "AffiTox robustness analysis — code, results, conclusions",
            [
                "Repository: https://github.com/ai-chem/AffiTox",
                "Local clone: /mnt/tank/scratch/okonovalova/AffiTox",
                "Seed 42. Bootstrap/permutation B = 10,000. Balanced nEF B = 2,000.",
                "Two identical analysis runs (CSV MD5 00e68d0c74b689ae4bf2845961257903).",
                "The main manuscript TeX was not edited. Proposed wording is in manuscript_review_text.md.",
                "",
                "This PDF documents the new robustness modules (Parts A–D), how denominators and",
                "score directions were defined, what was computable, what was not, and the numerical",
                "results. Terminology: distinct and potentially nonredundant evaluation axes;",
                "partially associated but sometimes discordant endpoints. Statistical independence",
                "was not tested and is not claimed.",
            ],
            fontsize=10,
        )

        _page(
            pdf,
            "1. Scientific objective",
            [
                "Three limitations of the published AffiTox evaluation were addressed:",
                "1. No formal cross-axis association between pose validity and affinity/screening.",
                "2. Incomplete characterization of potential training-data exposure.",
                "3. Insufficient evidence for prospective or temporally novel transfer.",
                "",
                "Unit of analysis for Part A: one target × method observation, 16 targets × 4",
                "pose-generating methods (Boltz-2, DynamicBind, GNINA 1.3, QVina2). PLAPT is",
                "excluded because it does not generate poses.",
                "",
                "Scoring, ranking, screening and docking power are kept as separate estimands.",
                "Missing poses are not imputed as PoseBusters failures in the conditional rate.",
                "They are unsuccessful only in the end-to-end rate.",
                "",
                "pose_generation_coverage = n_evaluable / n_input",
                "strict_pass_conditional = n_pass_all / n_evaluable   (missing poses excluded)",
                "strict_pass_end_to_end = n_pass_all / n_input         (missing = unsuccessful)",
            ],
        )

        _page(
            pdf,
            "2. Code map and what each module does",
            [
                "src/analysis/run_robustness.py — orchestrator: inventory, Parts A–D, second run,",
                "MD5 identity check, manuscript_review_text.md, robustness_validation_report.md.",
                "",
                "src/analysis/robustness_common.py — paths, input-inventory helper, PoseBusters",
                "ligand collapse (token CHEMBL/ligand_N → merged SMILES), two-way residualization",
                "y_ij − mean_i − mean_j + grand mean, cluster bootstrap Spearman, restricted",
                "within-target permutation of the pose endpoint with residualization inside each",
                "permutation, evidence-level vocabulary, Morgan ECFP4 helper.",
                "",
                "src/analysis/cross_axis_association.py — long table 64 rows; primary pair",
                "end-to-end pass vs Pearson r; secondary family with Benjamini–Hochberg;",
                "rank agreement; leave-one-target-out; figures.",
                "",
                "src/analysis/training_exposure_extension.py — method-specific evidence catalog",
                "from overlap_fractions_by_method.csv; DynamicBind ECFP4 similarity proxy vs",
                "processed-table crystal ligands; descriptive PDB-flag performance; NOT_COMPUTABLE.md.",
                "Proxies are never relabelled as exact membership. QVina2 overlap is N/A, not zero.",
                "",
                "src/analysis/temporal_validation.py — method-specific cutoffs (Boltz-2 year split",
                "is not copied onto PLAPT/DynamicBind/GNINA); eligibility n≥10 and pKi range ≥1.0;",
                "within-target Δr = r_post − r_pre with independent ligand bootstrap; candidate",
                "in-sample post-2022 manifest; feasibility simulation.",
                "",
                "src/analysis/screening_balanced_sensitivity.py — within-target equal-n resampling",
                "of actives and inactives (B=2000). Changes the estimand; labelled sensitivity.",
                "Active enrichment: stronger predicted binding first. Inactive: reversed direction.",
                "",
                "tests/test_robustness.py — formulas, score direction, residualization, BH,",
                "reproducibility of the 16×4 long table.",
            ],
        )

        _page(
            pdf,
            "3. Statistical methods (implemented as specified)",
            [
                "Primary association: Spearman of strict_pass_end_to_end vs Pearson r on 64 rows.",
                "Uncertainty: 10,000 target-cluster bootstraps (resample 16 targets with replacement).",
                "Two-way adjustment: residualize pose and affinity on C(target)+C(method), then Spearman.",
                "Permutation: shuffle the pose endpoint among the four methods within each target;",
                "repeat residualization inside every permutation; p = (#|ρ_perm|≥|ρ_obs|+1)/(B+1).",
                "Secondary pairs: same machinery; BH on restricted-permutation p-values only.",
                "",
                "Score conventions (unchanged from AffiTox): Pearson uses SIGN_FLIP for QVina2 and",
                "Boltz-2 affinity (lower-better). Enrichment does not use that flip; detect_score_direction",
                "ranks stronger binding first for actives and the opposite end for nEF10_low.",
                "GNINA primary score is gnina_cnn_affinity_bestpose (higher-better).",
                "",
                "Duplicate ligand–target keys are rejected in the merged panel and in the long table.",
                "Target-level rows are never treated as ligand-level observations.",
                "n_input is checked against merged row counts; n_pearson against published n_points",
                "(cross_axis_n_audit.csv: all 64 rows matched).",
            ],
        )

        _page(
            pdf,
            "4. Results — Part A cross-axis association",
            [
                f"Panel: {len(long)} observations, {long['target'].nunique()} targets, "
                f"{long['method'].nunique()} methods.",
                f"Primary raw Spearman ρ = {prim['spearman_raw']:.3f} "
                f"(95% cluster CI {prim['spearman_raw_ci_low']:.3f} to {prim['spearman_raw_ci_high']:.3f}).",
                f"Two-way-adjusted ρ = {prim['spearman_two_way_adjusted']:.3f} "
                f"(CI {prim['spearman_adj_ci_low']:.3f} to {prim['spearman_adj_ci_high']:.3f}).",
                f"Restricted permutation p (adjusted) = {prim['p_perm_adjusted_restricted']:.4g}.",
                "",
                "Interpretation: unadjusted, pose end-to-end validity and scoring power are positively",
                "associated. After removing additive target and method main effects, the residual",
                "association is near zero and the permutation test does not reject no residual association.",
                "This is evidence of partially associated but sometimes discordant endpoints, not of",
                "statistical independence.",
                "",
                f"Same unique top method on Pearson r vs end-to-end validity: "
                f"{100 * float(pair['proportion_same_top_method']):.0f}% of targets "
                f"(mean |rank difference| {float(pair['mean_of_mean_abs_rank_diff']):.2f}; "
                f"{int(pair['n_targets_discordant_top'])} discordant tops).",
                "Active vs inactive screening tops agree on only 31% of targets (rank table).",
                "",
                "Leave-one-target-out: raw Spearman stays between 0.40 and 0.47; adjusted Spearman",
                "stays near zero (range about −0.14 to +0.05). No single target drives the primary result.",
                "",
                "Secondary BH-adjusted permutation p-values are all large (smallest BH ≈ 0.72).",
                "Pose vs nEF10_active is essentially unassociated even before adjustment (raw ρ ≈ 0.02).",
            ],
        )

        show_assoc = assoc[
            [
                "family",
                "pose_endpoint",
                "affinity_endpoint",
                "spearman_raw",
                "spearman_raw_ci_low",
                "spearman_raw_ci_high",
                "spearman_two_way_adjusted",
                "p_perm_adjusted_restricted",
                "p_perm_adjusted_bh",
            ]
        ].copy()
        _table_page(pdf, "Table. Cross-axis Spearman associations (B=10,000)", show_assoc)

        _figures(pdf)

        _table_page(
            pdf,
            "Table. Leave-one-target-out (primary pair)",
            loto.rename(columns={"held_out_label": "held_out", "spearman_two_way_adjusted": "rho_adj"}),
        )

        _page(
            pdf,
            "5. Results — Part B training exposure",
            [
                "Evidence levels used: exact manifest intersection; processed-source-table;",
                "temporal eligibility proxy; source-domain exposure; similarity proxy;",
                "indeterminate; not applicable.",
                "",
                "Boltz-2: source-domain 100% and year≤2022 ≈ 88.6% are proxies. Exact pair",
                "membership is indeterminate (no official manifest).",
                "PLAPT: exact first-100k membership indeterminate (no frozen snapshot/row order).",
                "DynamicBind: processed-table intersections are exact for that released CSV",
                "(PDB 3/16; pairs 44/10,703). ECFP4 vs crystal ligands is a similarity proxy.",
                "GNINA 1.3: exact PDB in CrossDocked v1.3 training types 8/16. Ligand/pair",
                "membership indeterminate (no molcache).",
                "QVina2: supervised bioactivity overlap is not applicable, not zero.",
                "",
                "Descriptive target-level median Pearson r (PDB flags only, not leakage tests):",
                "GNINA types-exposed 8 targets median r ≈ 0.13 vs not-exposed ≈ 0.30;",
                "DynamicBind processed-PDB exposed 3 targets ≈ 0.27 vs 13 not-exposed ≈ 0.29.",
                "Sample sizes are too small for a causal leakage claim. These strata must not be",
                "compared with the Boltz-2 88.6% year proxy as if they were the same estimand.",
            ],
        )
        _table_page(pdf, "Table. Training-exposure evidence by method", bym)
        _table_page(
            pdf,
            "Table. Evidence rows (truncated fields)",
            ev[["method", "quantity", "n_items", "n_overlap", "frac", "evidence_level", "exact_or_proxy"]],
        )
        _table_page(pdf, "Table. Descriptive performance vs PDB exposure flags", exp)

        both = delta.loc[delta["both_eligible"] == True]  # noqa: E712
        _page(
            pdf,
            "6. Results — Part C temporal sensitivity",
            [
                "Label: temporal sensitivity / external temporal validation — not prospective validation.",
                "A retrospective post-cutoff split cannot become prospective validation by reanalysis.",
                "True prospective validation needs frozen models and rules before new labels exist.",
                "",
                "Eligibility: n_scored ≥ 10 and pKi range ≥ 1.0. Four AffiTox targets have both",
                "publication-year strata eligible: JAK1 (3EYG), JAK2 (3JY9), JAK3 (3LXK), CRBN (4TZ4).",
                "A 16-target pre-cutoff median is not compared with this four-target post-cutoff set.",
                "",
                "Boltz-2 Δr = r(≥2023) − r(≤2022), independent ligand bootstrap B=10,000:",
                "JAK1: −0.42 (CI −1.04 to 0.07); JAK2: −0.20 (−0.87 to 0.12);",
                "JAK3: +0.45 (0.22 to 0.67); CRBN: −0.29 (−0.53 to −0.06).",
                "Signs differ across targets; several intervals include zero. This is not a powered",
                "prospective test.",
                "",
                "Feasibility: observed Boltz-2 dual-target mean Δr ≈ −0.11 (SD ≈ 0.39). Even 32",
                "targets × 200 ligands/stratum still leaves a wide CI on the cross-target median",
                "(width ≈ 0.32 in the simulation). The current four dual-eligible targets are not",
                "enough for a precise multi-target median difference.",
                "",
                "The candidate manifest lists in-sample AffiTox records with document year ≥2023",
                "(already inside the panel). New Ki after model freeze are not in this repository.",
            ],
        )
        _table_page(
            pdf,
            "Table. Within-target Δr (both strata eligible)",
            both[
                [
                    "target_label",
                    "method",
                    "n_pre",
                    "n_post",
                    "r_pre",
                    "r_post",
                    "delta_r",
                    "delta_r_ci_low",
                    "delta_r_ci_high",
                ]
            ],
        )
        _table_page(
            pdf,
            "Table. Feasibility of a cross-target median Δr (simulation, not a sample-size rule)",
            feas,
        )

        _page(
            pdf,
            "7. Results — Part D balanced-class nEF10 sensitivity",
            [
                "Equal-n resampling of actives and inactives without replacement, B=2,000 per",
                "target×method, changes the population and the estimand. It is not a replacement",
                "for confirmatory unbalanced AffiTox nEF10.",
                "",
                "Example BCL-2 / Boltz-2: unbalanced active nEF10 = 0.96; balanced mean ≈ 0.95.",
                "Inactive nEF10_low moves more (0.52 unbalanced vs ≈ 0.61 balanced), as expected",
                "when prevalence is altered.",
                "Directions in the table match detect_score_direction (Boltz-2/QVina2 active=asc;",
                "GNINA CNN-affinity and DynamicBind/PLAPT active=desc).",
            ],
        )
        _table_page(
            pdf,
            "Table. Balanced vs unbalanced nEF10 (first 16 rows)",
            bal[
                [
                    "target_label",
                    "method",
                    "n_active_unbalanced",
                    "n_inactive_unbalanced",
                    "unbalanced_nEF10_active",
                    "balanced_nEF10_active_mean",
                    "unbalanced_nEF10_low",
                    "balanced_nEF10_low_mean",
                ]
            ].head(16),
        )

        _page(
            pdf,
            "8. What was computable / not computable",
            [
                "Computable: ligand-level PoseBusters panel joined to merged tables; 16×4 cross-axis",
                "Spearman with cluster bootstrap and restricted permutation; rank nonredundancy;",
                "LOTO; method-specific exposure catalog plus DynamicBind ECFP4 proxy; temporal Δr",
                "on four dual-eligible targets; in-sample post-2022 manifest; feasibility table;",
                "balanced nEF sensitivity.",
                "",
                "Not computable (no guessed numbers; see NOT_COMPUTABLE.md):",
                "— Boltz-2 exact affinity-training pair membership (no official manifest).",
                "— PLAPT exact first-100k membership (no frozen parquet + row-order hash).",
                "— GNINA ligand/pair types membership (no molcache/SDF for types poses).",
                "— ECFP4 vs Boltz-2 / PLAPT / GNINA training SMILES.",
                "— Maximum protein sequence identity vs each engine's training receptors.",
                "— New post-cutoff Ki records outside AffiTox.",
            ],
        )

        _page(
            pdf,
            "9. Supported vs unsupported conclusions",
            [
                "Supported:",
                "— Pose validity and scoring/screening can be jointly tabulated on AffiTox.",
                "— They are partially associated in the raw 16×4 panel (primary raw ρ ≈ 0.44)",
                "  and sometimes discordant (adjusted ρ ≈ 0; 12% of targets disagree on the top",
                "  method for Pearson vs end-to-end validity; screening tops disagree more).",
                "— Exposure evidence is method-specific and must not be pooled as one leakage rate.",
                "— Post-cutoff AffiTox slices support temporal sensitivity, not prospective validation.",
                "",
                "Unsupported:",
                "— Statistical independence of docking power and scoring power.",
                "— Exact training-set membership for Boltz-2 or PLAPT.",
                "— That a high PoseBusters pass rate implies informative pKi ranking, or the converse.",
                "— That a retrospective year split is prospective validation.",
                "— That balanced nEF10 is the confirmatory screening estimand.",
            ],
        )

        _page(
            pdf,
            "10. Manuscript-ready wording (author review; TeX not edited)",
            [
                "Pose validity and affinity estimation remain distinct and potentially nonredundant",
                "evaluation axes on AffiTox. Across 16 targets and the four pose-generating methods,",
                "the prespecified association between ligand-level strict PoseBusters end-to-end pass",
                "rate and within-target Pearson r is raw Spearman ρ = 0.438 (95% target-cluster CI",
                "0.272 to 0.629); two-way-adjusted ρ = −0.059 (CI −0.387 to 0.264); restricted-",
                "permutation p = 0.6756; n = 64. The endpoints are partially associated but sometimes",
                "discordant: the same unique top method on Pearson r and strict end-to-end validity",
                "in 88% of targets (mean absolute rank difference 0.66). We do not interpret these",
                "results as statistical independence.",
                "",
                "A retrospective post-cutoff analysis cannot become a true prospective validation",
                "merely through statistical reanalysis. True prospective validation requires frozen",
                "models and analysis rules before new experimental labels are obtained.",
            ],
        )

        _page(
            pdf,
            "11. Reproduce",
            [
                "cd /mnt/tank/scratch/okonovalova/AffiTox",
                "export PYTHONPATH=src PYTHONUNBUFFERED=1 PYTHONHASHSEED=0",
                "/mnt/tank/scratch/okonovalova/miniconda3/envs/docking/bin/python -u \\",
                "  src/analysis/run_robustness.py --seed 42 --n-boot 10000 --n-perm 10000 --n-balanced 2000",
                "",
                "This PDF:",
                "python -u src/analysis/build_robustness_pdf.py",
                "",
                "On this cluster submit via sbatch -p aichem (see build_robustness_pdf.sh).",
                "Outputs live only under analysis/tables/robustness/ and analysis/figures/robustness/.",
            ],
            fontsize=9,
        )

        _page(
            pdf,
            "Appendix. Full source of the new modules",
            [
                "The following pages list every new analysis file used for Parts A–D, with line",
                "numbers. Existing AffiTox modules (enrichment.py, correlations.py, leakage_stratification.py,",
                "etc.) were reused, not rewritten. Existing published CSVs and manuscript TeX were",
                "not modified.",
            ],
        )
        for path in CODE_FILES:
            if path.exists():
                _code_pages(pdf, path)

    print(f"wrote {OUT}  size={OUT.stat().st_size} bytes")


if __name__ == "__main__":
    main()
