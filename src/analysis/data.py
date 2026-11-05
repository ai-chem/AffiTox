"""Load merged docking tables and PoseBusters CSVs."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .constants import POSEBUSTERS_METHOD_KEYS, SIGN_FLIP_METRICS


def merged_csv_path(merged_dir: Path, target: str) -> Path:
    return merged_dir / f"merged_ligands_docking_{target.lower()}.csv"


def load_merged(merged_dir: Path, target: str) -> pd.DataFrame:
    path = merged_csv_path(merged_dir, target)
    if not path.exists():
        raise FileNotFoundError(f"Missing merged table: {path}")
    return pd.read_csv(path)


def ensure_pki(df: pd.DataFrame, exp_col: str = "pValue") -> pd.Series:
    """Return pKi series; compute from Ki if pValue column absent."""
    if exp_col in df.columns:
        return pd.to_numeric(df[exp_col], errors="coerce")
    if "standard_value" in df.columns:
        ki = pd.to_numeric(df["standard_value"], errors="coerce")
        return 9.0 - np.log10(ki.clip(lower=1e-12))
    raise ValueError("Need pValue or standard_value column for experimental axis")


def prepare_xy(
    df: pd.DataFrame,
    metric_col: str,
    exp_col: str = "pValue",
) -> tuple[np.ndarray, np.ndarray]:
    """Paired finite (pKi, score) arrays with sign flip applied."""
    pki = ensure_pki(df, exp_col)
    score = pd.to_numeric(df[metric_col], errors="coerce")
    mask = pki.notna() & score.notna() & np.isfinite(pki) & np.isfinite(score)
    x = pki[mask].values.astype(float)
    y = score[mask].values.astype(float)
    if metric_col in SIGN_FLIP_METRICS:
        y = -y
    return x, y


def classify_activity(ki_nm: float) -> str | None:
    if ki_nm is None or (isinstance(ki_nm, float) and np.isnan(ki_nm)):
        return None
    v = float(ki_nm)
    if v <= 100.0:
        return "high"
    if v < 1000.0:
        return "medium"
    return "low"


def activity_mask(df: pd.DataFrame, set_name: str) -> np.ndarray:
    ki = pd.to_numeric(df.get("standard_value"), errors="coerce")
    classes = ki.map(classify_activity)
    if set_name == "high_active":
        return (classes == "high").fillna(False).values
    if set_name == "active":
        return classes.isin(["high", "medium"]).fillna(False).values
    if set_name == "low":
        return (classes == "low").fillna(False).values
    raise ValueError(set_name)


def resolve_posebusters_csv(
    cfg_pose_dir: Path,
    cfg_boltz_dir: Path,
    target: str,
    method_id: str,
    dynamicbind_label: str = "dynamicbind",
) -> Path | None:
    key = POSEBUSTERS_METHOD_KEYS.get(method_id, method_id)
    if method_id == "boltz2":
        path = cfg_boltz_dir / f"posebusters_results_{target.lower()}_{key}.csv"
        return path if path.exists() else None
    if method_id == "dynamicbind":
        for label in (dynamicbind_label, "dynamicbind", "dynamicbind_new", "dynamicbind_old"):
            path = cfg_pose_dir / f"posebusters_results_{target.lower()}_{label}.csv"
            if path.exists():
                return path
        return None
    path = cfg_pose_dir / f"posebusters_results_{target.lower()}_{key}.csv"
    return path if path.exists() else None
