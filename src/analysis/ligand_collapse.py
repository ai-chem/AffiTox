"""Pre-dock collapse of duplicate SMILES to one ligand–target row.

Experimental Ki values are averaged; the *first-occurring* existing ligand_id
in table order is kept (never renumbered) so docking artifacts for that id
remain the canonical files for the molecule.
"""

from __future__ import annotations

import re
from typing import Hashable

import numpy as np
import pandas as pd

_EXPERIMENTAL_MEAN_COLS = (
    "standard_value",
    "value",
    "pchembl_value",
    "pKi",
    "pValue",
    "exp_value",
    "exp_standard_value",
    "exp_pchembl_value",
)
_DOCK_PREFIXES = ("qvina_", "gnina_", "plapt_", "dynamicbind_", "boltz2_")
_LIGAND_NUM_RE = re.compile(r"(\d+)")


def parse_ligand_number(token: object) -> int | None:
    if token is None or (isinstance(token, float) and np.isnan(token)):
        return None
    m = re.fullmatch(r"ligand_(\d+)", str(token).strip())
    if not m:
        m = _LIGAND_NUM_RE.search(str(token))
        if not m:
            return None
    n = int(m.group(1))
    return n if n >= 1 else None


def _recompute_pki_activity(df: pd.DataFrame) -> pd.DataFrame:
    if "standard_value" not in df.columns:
        return df
    ki = pd.to_numeric(df["standard_value"], errors="coerce")
    pki = 9.0 - np.log10(ki.where(ki > 0))
    if "pKi" in df.columns:
        df["pKi"] = pki
    if "pValue" in df.columns:
        df["pValue"] = pki
    if "is_active" in df.columns:
        df["is_active"] = ki < 1000.0
    if "activity_class" in df.columns:
        df["activity_class"] = np.where(
            ki <= 100.0, "high", np.where(ki < 1000.0, "medium", "low")
        )
    return df


def duplicate_collapse_plan(df: pd.DataFrame, smiles_col: str = "canonical_smiles") -> pd.DataFrame:
    """Return one row per dropped ligand: keep_ligand_id, drop_ligand_id, smiles, …"""
    if smiles_col not in df.columns or "ligand_id" not in df.columns:
        return pd.DataFrame()

    work = df.copy().reset_index(drop=True)
    work["_ord"] = np.arange(len(work))
    group_cols: list[Hashable] = [smiles_col]
    if "molecule_chembl_id" in work.columns:
        group_cols.append("molecule_chembl_id")

    rows = []
    for key, g in work.groupby(group_cols, sort=False):
        if len(g) < 2:
            continue
        g = g.sort_values("_ord", kind="mergesort")
        keep = str(g.iloc[0]["ligand_id"])
        keep_chembl = g.iloc[0].get("molecule_chembl_id", "")
        smi = g.iloc[0][smiles_col]
        for _, r in g.iloc[1:].iterrows():
            rows.append(
                {
                    "canonical_smiles": smi,
                    "molecule_chembl_id_keep": keep_chembl,
                    "molecule_chembl_id_drop": r.get("molecule_chembl_id", ""),
                    "keep_ligand_id": keep,
                    "drop_ligand_id": str(r["ligand_id"]),
                    "keep_ligand_number": parse_ligand_number(keep),
                    "drop_ligand_number": parse_ligand_number(r["ligand_id"]),
                }
            )
    return pd.DataFrame(rows)


def group_by_smiles(df: pd.DataFrame, smiles_col: str = "canonical_smiles") -> pd.DataFrame:
    """Collapse duplicate SMILES *before* docking (mean Ki; first-occurring ligand_id).

    - Group by canonical_smiles and molecule_chembl_id (when present).
    - Arithmetic mean of experimental Ki-related columns only.
    - Keep the first-occurring ligand_id in table order — never invent ids.
    - If docking columns exist, copy them from that representative row only.
    - Recompute pKi / is_active / activity_class from the mean Ki.
    """
    if smiles_col not in df.columns or len(df) == 0:
        return df

    mask = df[smiles_col].notna() & (
        df[smiles_col].astype(str).str.strip() != ""
    ) & (df[smiles_col].astype(str).str.lower() != "nan")
    work = df.loc[mask].copy().reset_index(drop=True)
    if len(work) == 0:
        return df

    group_cols: list[Hashable] = [smiles_col]
    if "molecule_chembl_id" in work.columns:
        group_cols.append("molecule_chembl_id")

    if work.groupby(group_cols).ngroups == len(work):
        return _recompute_pki_activity(work.reset_index(drop=True))

    dock_cols = [c for c in work.columns if c.startswith(_DOCK_PREFIXES)]
    exp_mean_cols = [
        c
        for c in _EXPERIMENTAL_MEAN_COLS
        if c in work.columns and pd.api.types.is_numeric_dtype(work[c])
    ]

    if "ligand_id" not in work.columns:
        agg: dict = {}
        for col in work.columns:
            if col in group_cols:
                continue
            if pd.api.types.is_numeric_dtype(work[col]):
                agg[col] = "mean"
            else:
                agg[col] = "first"
        out = work.groupby(group_cols, as_index=False).agg(agg)
        return _recompute_pki_activity(out)

    work["_ord"] = np.arange(len(work))
    rows: list[dict] = []
    for _, g in work.groupby(group_cols, sort=False):
        g = g.sort_values("_ord", kind="mergesort")
        rep = g.iloc[0]  # first occurring in table order
        row = rep.drop(labels=["_ord"]).to_dict()
        for col in exp_mean_cols:
            row[col] = pd.to_numeric(g[col], errors="coerce").mean()
        row["ligand_id"] = rep["ligand_id"]
        for col in dock_cols:
            row[col] = rep[col]
        rows.append(row)

    out = pd.DataFrame(rows)
    # Restore approximate original order by first occurrence
    out["_ord"] = out["ligand_id"].map(
        {str(r["ligand_id"]): i for i, r in work.iterrows()}
    )
    out = out.sort_values("_ord", kind="mergesort").drop(columns=["_ord"])
    ordered = [c for c in df.columns if c in out.columns]
    ordered += [c for c in out.columns if c not in ordered]
    out = out[ordered]
    return _recompute_pki_activity(out.reset_index(drop=True))
