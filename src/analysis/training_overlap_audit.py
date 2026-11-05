#!/usr/bin/env python3
"""Training-corpus overlap audit tables for SI Section 20."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "analysis" / "tables"
LIGANDS = ROOT / "input" / "ligands_nodubl"

PDB_TO_CSV = {
    "1g5m": "BCL2_Ki_WT_ChEMBL_252",
    "2v5z": "MAO-B_Ki_WT_ChEMBL_246",
    "3eyg": "JAK1_Ki_WT_ChEMBL_2255",
    "3jy9": "JAK2_Ki_WT_ChEMBL_2027",
    "3lxk": "JAK3_Ki_WT_ChEMBL_786",
    "11ue": "PDGFRB_Ki_WT_ChEMBL_275",
    "4ase": "VEGFR2_Ki_WT_ChEMBL_875",
    "4f65": "FGFR1_Ki_WT_ChEMBL_134",
    "4tz4": "CRBN_Ki_WT_ChEMBL_127",
    "4zau": "EGFR_Ki_WT_curated_251",
    "5jkv": "CYP19A1_Aromatase_Ki_WT_ChEMBL_548",
    "5mo4": "ABL1_BCR-ABL_Ki_WT_ChEMBL_693",
    "6gqj": "KIT_Ki_WT_curated_1298",
    "6jok": "PDGFRA_Ki_WT_curated_250",
    "5lf3": "PSMB5_Ki_WT_ChEMBL_88",
    "7kk3": "PARP1_Ki_WT_ChEMBL_1075",
}

PDB_DEPOSITION_YEAR = {
    "1G5M": 2001,
    "2V5Z": 2007,
    "3EYG": 2009,
    "3JY9": 2009,
    "3LXK": 2010,
    "11UE": 2026,
    "4ASE": 2012,
    "4F65": 2012,
    "4TZ4": 2014,
    "4ZAU": 2015,
    "5JKV": 2017,
    "5MO4": 2017,
    "6GQJ": 2018,
    "6JOK": 2020,
    "5LF3": 2016,
    "7KK3": 2021,
}

TARGET_LABEL = {
    "1G5M": "BCL-2",
    "2V5Z": "MAO-B",
    "3EYG": "JAK1",
    "3JY9": "JAK2",
    "3LXK": "JAK3",
    "11UE": "PDGFRB",
    "4ASE": "VEGFR2",
    "4F65": "FGFR1",
    "4TZ4": "CRBN",
    "4ZAU": "EGFR",
    "5JKV": "CYP19A1",
    "5MO4": "ABL1",
    "6GQJ": "c-KIT",
    "6JOK": "PDGFRA",
    "5LF3": "PSMB5",
    "7KK3": "PARP1",
}

# Canonical human UniProt accessions for the AffiTox assay targets
# (affinity-pair identity). 5LF3 is the PSMB5 surrogate; the ChEMBL assay
# target is PSMB5 (P28074).
PDB_TO_UNIPROT = {
    "1G5M": "P10415",
    "2V5Z": "P27338",
    "3EYG": "P23458",
    "3JY9": "O60674",
    "3LXK": "P52333",
    "11UE": "P09619",
    "4ASE": "P35968",
    "4F65": "P11362",
    "4TZ4": "Q96SW2",
    "4ZAU": "P00533",
    "5JKV": "P11511",
    "5MO4": "P00519",
    "6GQJ": "P10721",
    "6JOK": "P16234",
    "5LF3": "P28074",
    "7KK3": "P09874",
}

METHODS = [
    "boltz2_affinity_pred_value",
    "dynamicbind_affinity_bestpose",
    "gnina_cnn_affinity_bestpose",
    "plapt_affinity",
    "qvina_affinity_bestpose",
]


def dynamicbind_structure_era(year: int) -> str:
    return "Training-era (<2019 deposition)" if year < 2019 else "Holdout-like (2019+ deposition)"


def pdbbind_v2020_likelihood(pdb: str, year: int) -> str:
    if year >= 2021:
        return "Unlikely (deposited after PDBbind v2020 compilation)"
    if year == 2020:
        return "Borderline (2020 deposition; outside DynamicBind train split)"
    if pdb in {"3EYG", "3JY9", "3LXK", "11UE", "4ASE", "4F65", "4TZ4", "4ZAU", "5MO4", "6GQJ", "1G5M"}:
        return "High (frequent PDBbind/CrossDocked kinase/oncoprotein entry)"
    return "Moderate (target class represented in PDBbind v2020)"


def load_benchmark_ligands() -> pd.DataFrame:
    frames = []
    for pdb, stem in PDB_TO_CSV.items():
        path = LIGANDS / f"{stem}_nodubl.csv"
        df = pd.read_csv(path, sep=";")
        df["pdb"] = pdb.upper()
        df["document_year"] = pd.to_numeric(df["document_year"], errors="coerce")
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def target_structure_audit() -> pd.DataFrame:
    rows = []
    for pdb, year in sorted(PDB_DEPOSITION_YEAR.items(), key=lambda item: item[1]):
        rows.append(
            {
                "pdb": pdb,
                "target": TARGET_LABEL[pdb],
                "deposition_year": year,
                "dynamicbind_structure_era": dynamicbind_structure_era(year),
                "pdbbind_v2020_likelihood": pdbbind_v2020_likelihood(pdb, year),
                "boltz2_structure_pool": "Included (deposited before 2023-06-01 cutoff)",
            }
        )
    return pd.DataFrame(rows)


def ligand_publication_audit(ligands: pd.DataFrame) -> pd.DataFrame:
    uniq = ligands.drop_duplicates(["pdb", "canonical_smiles"])
    summary = (
        uniq.groupby("pdb")
        .agg(
            n_ligands=("canonical_smiles", "count"),
            median_publication_year=("document_year", "median"),
            pct_publication_le_2018=("document_year", lambda s: 100 * (s <= 2018).mean()),
            pct_publication_ge_2019=("document_year", lambda s: 100 * (s >= 2019).mean()),
        )
        .reset_index()
    )
    return summary


def structure_holdout_correlations() -> pd.DataFrame:
    corr = pd.read_csv(TABLES / "correlations/correlations_with_ci_perm.csv")
    corr = corr[corr.metric.isin(METHODS)].copy()
    corr["structure_era"] = corr["target"].str.upper().map(
        lambda t: "pre-2019 deposition" if PDB_DEPOSITION_YEAR[t] < 2019 else "2019+ deposition"
    )
    rows = []
    for (era, method), sub in corr.groupby(["structure_era", "method_label"]):
        rows.append(
            {
                "structure_era": era,
                "method": method,
                "n_targets": len(sub),
                "median_pearson_r": float(sub["pearson_r"].median()),
                "mean_pearson_r": float(sub["pearson_r"].mean()),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    out = TABLES / "overlap"
    out.mkdir(parents=True, exist_ok=True)

    ligands = load_benchmark_ligands()
    uniq_smiles = ligands["canonical_smiles"].nunique()
    uniq_pairs = ligands.drop_duplicates(["pdb", "canonical_smiles"])

    target_structure_audit().to_csv(out / "target_structure_era_audit.csv", index=False)
    ligand_publication_audit(ligands).to_csv(out / "ligand_publication_year_by_target.csv", index=False)
    structure_holdout_correlations().to_csv(out / "structure_era_correlation_summary.csv", index=False)

    summary = {
        "unique_smiles": int(uniq_smiles),
        "unique_target_ligand_pairs": int(len(uniq_pairs)),
        "pct_pairs_publication_le_2018": float((uniq_pairs["document_year"] <= 2018).mean()),
        "pct_pairs_publication_ge_2019": float((uniq_pairs["document_year"] >= 2019).mean()),
        "n_targets_pre2019_deposition": int(sum(y < 2019 for y in PDB_DEPOSITION_YEAR.values())),
        "n_targets_post2019_deposition": int(sum(y >= 2019 for y in PDB_DEPOSITION_YEAR.values())),
    }
    pd.DataFrame([summary]).to_csv(out / "benchmark_overlap_summary.csv", index=False)
    print(f"Wrote training overlap audit to {out}")
    print(summary)


if __name__ == "__main__":
    main()
