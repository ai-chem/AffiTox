#!/usr/bin/env python3
"""Sensitivity statistics for the canonical 16-target AffiTox panel.

Recomputes, from the merged score tables already on this panel:

- Pearson r and active nEF10 on method-available, pairwise-common, and
  five-method-common ligand masks, then Wilcoxon/Holm across targets;
- scaffold-cluster and assay-cluster bootstrap intervals for Pearson r;
- the same cluster bootstrap for the cross-target median Pearson difference
  on the five-method mask;
- publication-year common-mask Pearson r and within-target Δr.

Training-corpus intersections are not recomputed here: the DynamicBind
processed table and the GNINA CrossDocked types archive are not in this tree.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from numpy import log10
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import pearsonr

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from analysis.constants import PRIMARY_METRICS, SIGN_FLIP_METRICS
from analysis.data import ensure_pki, load_merged
from analysis.enrichment import compute_ef_for_metric
from analysis.inferential import holm, wilcoxon_paired
from analysis.panel import CANONICAL_TARGETS
from analysis.training_overlap_audit import LIGANDS, PDB_TO_CSV

RDLogger.DisableLog("rdApp.*")

ROOT = PROJECT_ROOT / "analysis" / "excluding_2z5x_3mjg" / "tables"
OUT = ROOT / "sensitivity"
PUBLISHED = ROOT / "correlations" / "correlations_with_ci_perm.csv"
N_BOOT = 2000
SEED = 42
METHODS = list(PRIMARY_METRICS.keys())


def _murcko(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return ""
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or ""
    except Exception:
        return ""


def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.std(x, ddof=1) == 0 or np.std(y, ddof=1) == 0:
        return np.nan
    r, _ = pearsonr(x, y)
    return float(r) if np.isfinite(r) else np.nan


def _nef10_active(df: pd.DataFrame, metric: str) -> float:
    res = compute_ef_for_metric(df, metric)
    if not res:
        return np.nan
    return float(res["nEF10_active"])


def _cluster_ids(labels: np.ndarray) -> np.ndarray:
    """Integer cluster ids. Empty labels stay as singletons."""
    ids = np.empty(len(labels), dtype=int)
    next_id = 0
    known: dict[str, int] = {}
    for i, lab in enumerate(labels):
        if not lab:
            ids[i] = next_id
            next_id += 1
            continue
        if lab not in known:
            known[lab] = next_id
            next_id += 1
        ids[i] = known[lab]
    return ids


def _resample_indices(cluster_ids: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    uniq = np.unique(cluster_ids)
    drawn = rng.choice(uniq, size=len(uniq), replace=True)
    parts = [np.flatnonzero(cluster_ids == c) for c in drawn]
    return np.concatenate(parts) if parts else np.array([], dtype=int)


def _one_row_per_smiles(df: pd.DataFrame) -> pd.DataFrame:
    """Same pre-dock SMILES collapse as analysis.ligand_collapse.group_by_smiles."""
    from analysis.ligand_collapse import group_by_smiles

    return group_by_smiles(df)


def load_panel() -> dict[str, dict]:
    panel = {}
    for target in CANONICAL_TARGETS:
        df = _one_row_per_smiles(load_merged(ROOT, target))
        pki = ensure_pki(df).to_numpy(dtype=float)
        smiles = df["canonical_smiles"].astype(str).to_numpy()
        assay = df["assay_chembl_id"].astype(str).fillna("").to_numpy()
        scaffolds = np.array([_murcko(s) for s in smiles], dtype=object)
        scores = {}
        for method_id, col in PRIMARY_METRICS.items():
            y = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
            if col in SIGN_FLIP_METRICS:
                y = -y
            scores[method_id] = y
        years = _publication_years(target, smiles)
        panel[target] = {
            "df": df,
            "pki": pki,
            "smiles": smiles,
            "assay": assay,
            "scaffold": scaffolds,
            "scaffold_id": _cluster_ids(scaffolds),
            "assay_id": _cluster_ids(assay),
            "scores": scores,
            "year": years,
        }
    return panel


def _publication_years(target: str, smiles: np.ndarray) -> np.ndarray:
    stem = PDB_TO_CSV[target]
    nodubl = pd.read_csv(LIGANDS / f"{stem}_nodubl.csv", sep=";")
    nodubl["document_year"] = pd.to_numeric(nodubl["document_year"], errors="coerce")
    year = (
        nodubl.groupby(nodubl["canonical_smiles"].astype(str))["document_year"].min()
    )
    return pd.Series(smiles).map(year).to_numpy(dtype=float)


def _finite_mask(pki: np.ndarray, cols: list[np.ndarray]) -> np.ndarray:
    mask = np.isfinite(pki)
    for col in cols:
        mask &= np.isfinite(col)
    return mask


def point_estimates(panel: dict[str, dict]) -> tuple[pd.DataFrame, pd.DataFrame]:
    coef_rows = []
    for target, rec in panel.items():
        pki = rec["pki"]
        df = rec["df"]
        all_scores = [rec["scores"][m] for m in METHODS]
        common = _finite_mask(pki, all_scores)
        for method_id, col in PRIMARY_METRICS.items():
            own = _finite_mask(pki, [rec["scores"][method_id]])
            for mask_name, mask in (("method_available", own), ("common_all_five", common)):
                coef_rows.append(
                    {
                        "target": target,
                        "method_id": method_id,
                        "mask": mask_name,
                        "N": int(mask.sum()),
                        "pearson_r": _pearson(pki[mask], rec["scores"][method_id][mask]),
                        "nEF10_active": _nef10_active(df.loc[mask], col),
                    }
                )
        for i, m1 in enumerate(METHODS):
            for m2 in METHODS[i + 1 :]:
                both = _finite_mask(pki, [rec["scores"][m1], rec["scores"][m2]])
                for method_id in (m1, m2):
                    coef_rows.append(
                        {
                            "target": target,
                            "method_id": method_id,
                            "mask": f"pairwise_{m1}_{m2}",
                            "N": int(both.sum()),
                            "pearson_r": _pearson(pki[both], rec["scores"][method_id][both]),
                            "nEF10_active": _nef10_active(
                                df.loc[both], PRIMARY_METRICS[method_id]
                            ),
                        }
                    )
    coef = pd.DataFrame(coef_rows)
    tests = _wilcoxon_table(coef)
    return coef, tests


def _wilcoxon_table(coef: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for mask, sub in coef.groupby("mask"):
        present = [m for m in METHODS if m in set(sub["method_id"])]
        for endpoint in ("pearson_r", "nEF10_active"):
            raw = []
            for i, m1 in enumerate(present):
                for m2 in present[i + 1 :]:
                    a = sub.loc[sub.method_id == m1].set_index("target")[endpoint]
                    b = sub.loc[sub.method_id == m2].set_index("target")[endpoint]
                    targets = sorted(set(a.index) & set(b.index))
                    va = [float(a.loc[t]) for t in targets]
                    vb = [float(b.loc[t]) for t in targets]
                    p, n_nz, n_tgt, med, mean = wilcoxon_paired(va, vb)
                    raw.append(
                        {
                            "mask": mask,
                            "endpoint": endpoint,
                            "method_a": m1,
                            "method_b": m2,
                            "n_targets": n_tgt,
                            "n_nonzero_diff": n_nz,
                            "median_diff_a_minus_b": med,
                            "mean_diff_a_minus_b": mean,
                            "p_raw": p,
                        }
                    )
            pvals = np.array([r["p_raw"] if np.isfinite(r["p_raw"]) else 1.0 for r in raw])
            adj = holm(pvals)
            for rec, p_holm in zip(raw, adj):
                rec["p_holm"] = float(p_holm)
                rows.append(rec)
    return pd.DataFrame(rows)


def _cluster_cis(panel: dict[str, dict], id_key: str, cluster_name: str) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    rows = []
    for target, rec in panel.items():
        pki = rec["pki"]
        cids = rec[id_key]
        common = _finite_mask(pki, [rec["scores"][m] for m in METHODS])
        for method_id in METHODS:
            own = _finite_mask(pki, [rec["scores"][method_id]])
            for mask_name, mask in (("method_available", own), ("common_all_five", common)):
                idx0 = np.flatnonzero(mask)
                point = _pearson(pki[idx0], rec["scores"][method_id][idx0])
                boots = np.empty(N_BOOT, dtype=float)
                local_ids = cids[idx0]
                # Remap so resampling uses only clusters present in the mask.
                local_ids = _cluster_ids(local_ids.astype(str))
                for b in range(N_BOOT):
                    take = _resample_indices(local_ids, rng)
                    sel = idx0[take]
                    boots[b] = _pearson(pki[sel], rec["scores"][method_id][sel])
                finite = boots[np.isfinite(boots)]
                rows.append(
                    {
                        "target": target,
                        "method_id": method_id,
                        "mask": mask_name,
                        "cluster": cluster_name,
                        "N": int(len(idx0)),
                        "n_clusters": int(len(np.unique(local_ids))),
                        "pearson_r": point,
                        "ci_low": float(np.percentile(finite, 2.5)) if len(finite) else np.nan,
                        "ci_high": float(np.percentile(finite, 97.5)) if len(finite) else np.nan,
                        "n_boot": N_BOOT,
                        "n_boot_finite": int(len(finite)),
                    }
                )
    return pd.DataFrame(rows)


def _ordering_bootstrap(panel: dict[str, dict], id_key: str, cluster_name: str) -> pd.DataFrame:
    """Bootstrap the cross-target median Pearson difference on the five-method mask."""
    rng = np.random.default_rng(SEED + 7)
    targets = list(panel)
    pair_deltas = {
        (m1, m2): np.empty(N_BOOT, dtype=float)
        for i, m1 in enumerate(METHODS)
        for m2 in METHODS[i + 1 :]
    }
    for b in range(N_BOOT):
        rs = {m: [] for m in METHODS}
        for target in targets:
            rec = panel[target]
            pki = rec["pki"]
            common = _finite_mask(pki, [rec["scores"][m] for m in METHODS])
            idx0 = np.flatnonzero(common)
            local_ids = _cluster_ids(rec[id_key][idx0].astype(str))
            sel = idx0[_resample_indices(local_ids, rng)]
            for method_id in METHODS:
                rs[method_id].append(_pearson(pki[sel], rec["scores"][method_id][sel]))
        for i, m1 in enumerate(METHODS):
            for m2 in METHODS[i + 1 :]:
                delta = np.asarray(rs[m1], dtype=float) - np.asarray(rs[m2], dtype=float)
                pair_deltas[(m1, m2)][b] = float(np.nanmedian(delta))
    rows = []
    for (m1, m2), draws in pair_deltas.items():
        finite = draws[np.isfinite(draws)]
        rows.append(
            {
                "cluster": cluster_name,
                "mask": "common_all_five",
                "endpoint": "pearson_r",
                "method_a": m1,
                "method_b": m2,
                "median_diff_point_boot_median": float(np.median(finite)) if len(finite) else np.nan,
                "ci_low": float(np.percentile(finite, 2.5)) if len(finite) else np.nan,
                "ci_high": float(np.percentile(finite, 97.5)) if len(finite) else np.nan,
                "fraction_median_diff_positive": float(np.mean(finite > 0)) if len(finite) else np.nan,
                "n_boot": N_BOOT,
            }
        )
    return pd.DataFrame(rows)


def temporal(panel: dict[str, dict]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(SEED + 11)
    coef_rows = []
    delta_rows = []
    for target, rec in panel.items():
        pki = rec["pki"]
        year = rec["year"]
        common = _finite_mask(pki, [rec["scores"][m] for m in METHODS])
        strata = {
            "pubyear_le_2022": common & np.isfinite(year) & (year <= 2022),
            "pubyear_ge_2023": common & np.isfinite(year) & (year >= 2023),
        }
        stored = {}
        for stratum, mask in strata.items():
            idx = np.flatnonzero(mask)
            for method_id in METHODS:
                r = _pearson(pki[idx], rec["scores"][method_id][idx]) if len(idx) >= 10 else np.nan
                coef_rows.append(
                    {
                        "target": target,
                        "method_id": method_id,
                        "stratum": stratum,
                        "N": int(len(idx)),
                        "pearson_r": r,
                        "eligible": int(len(idx) >= 10),
                    }
                )
            stored[stratum] = idx
        pre = stored["pubyear_le_2022"]
        post = stored["pubyear_ge_2023"]
        both = len(pre) >= 10 and len(post) >= 10
        for method_id in METHODS:
            r_pre = _pearson(pki[pre], rec["scores"][method_id][pre]) if len(pre) >= 10 else np.nan
            r_post = _pearson(pki[post], rec["scores"][method_id][post]) if len(post) >= 10 else np.nan
            boots = np.empty(N_BOOT, dtype=float)
            if both:
                for b in range(N_BOOT):
                    ip = pre[rng.integers(0, len(pre), len(pre))]
                    io = post[rng.integers(0, len(post), len(post))]
                    boots[b] = _pearson(pki[io], rec["scores"][method_id][io]) - _pearson(
                        pki[ip], rec["scores"][method_id][ip]
                    )
                finite = boots[np.isfinite(boots)]
                lo = float(np.percentile(finite, 2.5)) if len(finite) else np.nan
                hi = float(np.percentile(finite, 97.5)) if len(finite) else np.nan
            else:
                lo = hi = np.nan
            delta = r_post - r_pre if np.isfinite(r_pre) and np.isfinite(r_post) else np.nan
            delta_rows.append(
                {
                    "target": target,
                    "method_id": method_id,
                    "n_pre": int(len(pre)),
                    "n_post": int(len(post)),
                    "r_pre": r_pre,
                    "r_post": r_post,
                    "delta_r_ge2023_minus_le2022": delta,
                    "delta_r_ci_low": lo,
                    "delta_r_ci_high": hi,
                    "both_eligible": int(both),
                    "n_boot": N_BOOT if both else 0,
                }
            )
    return pd.DataFrame(coef_rows), pd.DataFrame(delta_rows)


def check_against_published(coef: pd.DataFrame) -> None:
    published = pd.read_csv(PUBLISHED)
    published = published.loc[published.metric.isin(PRIMARY_METRICS.values())].copy()
    col_to_method = {v: k for k, v in PRIMARY_METRICS.items()}
    published["method_id"] = published["metric"].map(col_to_method)
    own = coef.loc[coef["mask"] == "method_available", ["target", "method_id", "pearson_r"]]
    merged = own.merge(
        published[["target", "method_id", "pearson_r"]],
        on=["target", "method_id"],
        suffixes=("_recomputed", "_published"),
    )
    if len(merged) != 16 * 5:
        raise AssertionError(f"Expected 80 published comparisons, got {len(merged)}")
    delta = (merged["pearson_r_recomputed"] - merged["pearson_r_published"]).abs()
    worst = float(delta.max())
    print(
        "method-available Pearson vs the earlier correlation table "
        f"(computed before the current merged files): max |Δ| = {worst:.3e}"
    )
    merged["abs_delta"] = delta
    merged.to_csv(OUT / "method_available_vs_published_pearson.csv", index=False)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"Loading {len(CANONICAL_TARGETS)} targets from {ROOT}")
    panel = load_panel()
    print("Point estimates and Wilcoxon...")
    coef, tests = point_estimates(panel)
    check_against_published(coef)
    coef.to_csv(OUT / "mask_coefficients.csv", index=False)
    tests.to_csv(OUT / "wilcoxon_by_mask.csv", index=False)
    print("Scaffold-cluster Pearson intervals...")
    scaff = _cluster_cis(panel, "scaffold_id", "murcko_scaffold")
    scaff.to_csv(OUT / "scaffold_cluster_pearson_ci.csv", index=False)
    print("Assay-cluster Pearson intervals...")
    assay = _cluster_cis(panel, "assay_id", "assay_chembl_id")
    assay.to_csv(OUT / "assay_cluster_pearson_ci.csv", index=False)
    print("Ordering bootstrap on the five-method mask...")
    order = pd.concat(
        [
            _ordering_bootstrap(panel, "scaffold_id", "murcko_scaffold"),
            _ordering_bootstrap(panel, "assay_id", "assay_chembl_id"),
        ],
        ignore_index=True,
    )
    order.to_csv(OUT / "ordering_bootstrap_common_mask.csv", index=False)
    print("Publication-year common-mask Pearson and Δr...")
    stratum, delta = temporal(panel)
    stratum.to_csv(OUT / "temporal_common_mask_pearson.csv", index=False)
    delta.to_csv(OUT / "temporal_delta_r.csv", index=False)
    note = OUT / "exposure_corpus_not_available.txt"
    note.write_text(
        "DynamicBind processed table (d3_with_clash_info.csv), the timesplit "
        "lists, LP-PDBBind, and the GNINA CrossDocked2020 v1.3 types archive "
        "are not in this tree. Training-set intersection counts were not "
        "recomputed. Publication-year rows above use only the local ChEMBL "
        "document_year joined to the merged scores.\n"
    )
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
