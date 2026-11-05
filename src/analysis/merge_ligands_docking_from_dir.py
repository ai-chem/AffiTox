#!/usr/bin/env python3
"""
Merge nodubl/curated ligand tables with docking metrics.

Pre-dock curation (enforced here before the score join):
  within each target, rows that share canonical_smiles (+ molecule_chembl_id
  when present) are collapsed by the arithmetic mean of experimental Ki fields.
  The surviving ligand_id is an *existing* id from that group (never renumbered),
  chosen deterministically as the minimum ligand_id. Docking scores are then
  joined only on that ligand_id so the merge cannot attach another molecule's
  poses/scores.

Assay-level nodubl rows may still exist on disk for provenance; the merged
analysis tables are unique ligand–target pairs after this collapse.
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path
from typing import Dict, Optional, List, Tuple

import pandas as pd

_SRC = Path(__file__).resolve().parent.parent
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from analysis.ligand_collapse import group_by_smiles
from analysis.ligand_identity import chembl_to_ligand_ids, native_to_ligand_id
from analysis.panel import (
    CANONICAL_TARGETS,
    LEGACY_PDB_TO_NODUBL_STEM,
    nodubl_csv_name,
)
from docking_benchmark2.ligand_ids import (
    format_ligand_id,
    normalize_ligand_id,
    parse_dynamicbind_idx,
    parse_ligand_number,
)

# Маппинг PDB кодов на имена CSV файлов в ligands_nodubl (без _nodubl.csv)
PDB_TO_CSV_MAPPING = {
    "1ere": "Estrogen_receptor_alpha_Ki_WT_ChEMBL_494",
    "1g5m": "BCL2_Ki_WT_ChEMBL_252",
    "1g5m_wrong": "BCL2L1_BCL-XL_Ki_WT_ChEMBL_287",
    "1tqn": "CYP3A4_Ki_WT_ChEMBL_1",
    "2v5z": "MAO-B_Ki_WT_ChEMBL_246",
    "2z5x": "MAO-B_Ki_WT_ChEMBL_246",
    "3eyg": "JAK1_Ki_WT_ChEMBL_2255",
    "3jy9": "JAK2_Ki_WT_ChEMBL_2027",
    "3lxk": "JAK3_Ki_WT_ChEMBL_786",
    "3mjg": "PDGFRB_Ki_WT_ChEMBL_275",
    "4ase": "VEGFR2_Ki_WT_ChEMBL_875",
    "4f65": "FGFR1_Ki_WT_ChEMBL_134",
    "4k7a": "Androgen_receptor_Ki_WT_ChEMBL_516",
    "4tz4": "CRBN_Ki_WT_ChEMBL_127",
    "4wnv": "HTR2B_5HT2B_IC50_WT_ChEMBL_61",
    "4zau": "EGFR_Ki_WT_curated_251",
    "5jkv": "CYP19A1_Aromatase_Ki_WT_ChEMBL_548",
    "5mo4": "ABL1_BCR-ABL_Ki_WT_ChEMBL_693",
    "6gqj": "KIT_Ki_WT_curated_1298",
    "6jok": "PDGFRA_Ki_WT_curated_250",
    "5lf3": "PSMB5_Ki_WT_ChEMBL_88",
    "7awe": "PSMB5_Ki_WT_ChEMBL_88",
    "7kk3": "PARP1_Ki_WT_ChEMBL_1075",
    "8zyq": "hERG_Ki_WT_curated_417",
    "11ue": "PDGFRB_Ki_WT_ChEMBL_275",
}

# Столбцы для извлечения из ligands_nodubl / ligands_curated CSV
LIGANDS_COLUMNS = [
    "ligand_id",
    "molecule_chembl_id",
    "canonical_smiles",
    "assay_chembl_id",
    "standard_value",
    "value",
    "pchembl_value",
    "pKi",
    "pValue",
    "is_active",
    "activity_class",
    "assay_type",
    "document_chembl_id",
    "type",
    "units",
    "uo_units",
]

# Столбцы boltz для добавления
BOLTZ_COLUMNS = [
    "boltz2_affinity_pred_value",
    "boltz2_affinity_probability_binary",
    "boltz2_affinity_pred_value1",
    "boltz2_affinity_probability_binary1",
    "boltz2_affinity_pred_value2",
    "boltz2_affinity_probability_binary2",
]

# Столбцы экспериментальных данных для добавления
EXP_COLUMNS = [
    "exp_value",
    "exp_standard_value",
    "exp_pchembl_value",
]

# Path to Boltz affinity JSON outputs (override via BOLTZ_RESULTS_DIR env or --boltz-results-dir)
BOLTZ_RESULTS_DIR = os.environ.get("BOLTZ_RESULTS_DIR", "")


def extract_chembl_id_from_path(json_path: str) -> Optional[str]:
    """Извлекает CHEMBL ID из пути к JSON файлу."""
    match = re.search(r'/(CHEMBL\d+)/', json_path)
    if match:
        return match.group(1)
    match = re.search(r'affinity_(CHEMBL\d+)\.json', json_path)
    if match:
        return match.group(1)
    return None


def parse_boltz_affinity_json(json_path: str) -> Optional[Dict[str, float]]:
    """Разбор одного Boltz affinity JSON файла."""
    try:
        with open(json_path, "r") as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError):
        return None

    result: Dict[str, float] = {}
    fields = [
        "affinity_pred_value",
        "affinity_probability_binary",
        "affinity_pred_value1",
        "affinity_probability_binary1",
        "affinity_pred_value2",
        "affinity_probability_binary2",
    ]
    
    for field in fields:
        if field in data:
            try:
                result[f"boltz2_{field}"] = float(data[field])
            except (ValueError, TypeError):
                pass
    
    return result if result else None


def collect_boltz_metrics(
    protein: str, boltz_results_dir: str = BOLTZ_RESULTS_DIR
) -> Dict[str, Dict[str, float]]:
    """
    Собирает метрики Boltz для данного белка.
    Возвращает словарь: {molecule_chembl_id: {boltz2_...: value, ...}, ...}
    
    Поддерживает две структуры директорий:
    1. Старая: {boltz_results_dir}/{protein}/docking/boltz2_*/predictions/{ligand_dataset}/boltz_results_{ligand_dataset}/predictions/{CHEMBL_ID}/affinity_{CHEMBL_ID}.json
    2. Новая: {boltz_results_dir}/{protein}/docking/boltz2_*/boltz_results_{CHEMBL_ID}/predictions/{CHEMBL_ID}/affinity_{CHEMBL_ID}.json
    """
    result: Dict[str, Dict[str, float]] = {}
    
    protein_dir = os.path.join(boltz_results_dir, protein, "docking")
    if not os.path.isdir(protein_dir):
        return result
    
    boltz_docking_dir = None
    for item in os.listdir(protein_dir):
        if item.startswith("boltz2_") and "aff" in item:
            boltz_docking_dir = os.path.join(protein_dir, item)
            break
    
    if not boltz_docking_dir or not os.path.isdir(boltz_docking_dir):
        return result
    
    # НОВАЯ СТРУКТУРА: boltz_results_{CHEMBL_ID}/predictions/{CHEMBL_ID}/affinity_{CHEMBL_ID}.json
    for item in os.listdir(boltz_docking_dir):
        if item.startswith("boltz_results_") and item.startswith("boltz_results_CHEMBL"):
            boltz_results_chembl_dir = os.path.join(boltz_docking_dir, item)
            if not os.path.isdir(boltz_results_chembl_dir):
                continue
            
            # Извлекаем CHEMBL ID из имени директории
            match = re.search(r'boltz_results_(CHEMBL\d+)', item)
            if not match:
                continue
            chembl_id = match.group(1)
            
            # Ищем файл affinity_{CHEMBL_ID}.json
            predictions_dir = os.path.join(boltz_results_chembl_dir, "predictions", chembl_id)
            if os.path.isdir(predictions_dir):
                json_file = os.path.join(predictions_dir, f"affinity_{chembl_id}.json")
                if os.path.isfile(json_file):
                    metrics = parse_boltz_affinity_json(json_file)
                    if metrics:
                        result[chembl_id] = metrics
    
    # СТАРАЯ СТРУКТУРА: predictions/{ligand_dataset}/boltz_results_{ligand_dataset}/predictions/{CHEMBL_ID}/affinity_{CHEMBL_ID}.json
    predictions_dir = os.path.join(boltz_docking_dir, "predictions")
    if os.path.isdir(predictions_dir):
        for ligand_dataset in os.listdir(predictions_dir):
            ligand_dataset_dir = os.path.join(predictions_dir, ligand_dataset)
            if not os.path.isdir(ligand_dataset_dir):
                continue
            
            boltz_results_subdir = os.path.join(
                ligand_dataset_dir, f"boltz_results_{ligand_dataset}"
            )
            if not os.path.isdir(boltz_results_subdir):
                continue
            
            # Ищем JSON файлы с метриками в поддиректории predictions
            predictions_subdir = os.path.join(boltz_results_subdir, "predictions")
            if os.path.isdir(predictions_subdir):
                for chembl_dir in os.listdir(predictions_subdir):
                    chembl_path = os.path.join(predictions_subdir, chembl_dir)
                    if not os.path.isdir(chembl_path):
                        continue
                    
                    # Ищем файл affinity_{CHEMBL_ID}.json
                    json_file = os.path.join(chembl_path, f"affinity_{chembl_dir}.json")
                    if os.path.isfile(json_file):
                        chembl_id = extract_chembl_id_from_path(json_file) or chembl_dir
                        # Пропускаем, если уже добавили из новой структуры
                        if chembl_id in result:
                            continue
                        metrics = parse_boltz_affinity_json(json_file)
                        if metrics:
                            result[chembl_id] = metrics
    
    return result


# group_by_smiles: see analysis.ligand_collapse (pre-dock mean-Ki collapse).


# Импортируем функции парсинга из extract_raw_docking_metrics.py
# (копируем их сюда, чтобы не зависеть от другого скрипта)

def parse_qvina_log(log_path: str) -> Optional[Dict[str, float]]:
    """Разбор одного qvina .log файла."""
    try:
        with open(log_path, "r") as f:
            lines = f.readlines()
    except OSError:
        return None

    data_rows: List[Tuple[int, float, float, float]] = []
    in_table = False

    for line in lines:
        line = line.strip()
        if "mode |" in line and "affinity" in line:
            in_table = True
            continue
        if not in_table:
            continue
        if line.startswith("-----"):
            continue
        if re.match(r"^\s*\d+\s+[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+$", line):
            parts = re.split(r"\s+", line)
            if len(parts) >= 4:
                try:
                    mode = int(parts[0])
                    affinity = float(parts[1])
                    rmsd_lb = float(parts[2])
                    rmsd_ub = float(parts[3])
                    data_rows.append((mode, affinity, rmsd_lb, rmsd_ub))
                except (ValueError, IndexError):
                    continue

    if not data_rows:
        return None

    min_aff_row = min(data_rows, key=lambda x: x[1])
    max_aff_row = max(data_rows, key=lambda x: x[1])
    bestpose_row = data_rows[0] if data_rows else None

    if bestpose_row is None:
        return None

    _, aff_min, rmsd_lb_min, rmsd_ub_min = min_aff_row
    _, aff_max, _, _ = max_aff_row
    _, aff_bp, rmsd_lb_bp, rmsd_ub_bp = bestpose_row

    lb_candidates = [row for row in data_rows if row[2] > 0]
    lb_aff = None
    if lb_candidates:
        lb_best = min(lb_candidates, key=lambda x: x[2])
        lb_aff = lb_best[1]

    ub_candidates = [row for row in data_rows if row[3] > 0]
    ub_aff = None
    if ub_candidates:
        ub_best = min(ub_candidates, key=lambda x: x[3])
        ub_aff = ub_best[1]

    return {
        "affinity_min": aff_min,
        "rmsd_lb_min": rmsd_lb_min,
        "rmsd_ub_min": rmsd_ub_min,
        "affinity_bestpose": aff_bp,
        "rmsd_lb_bestpose": rmsd_lb_bp,
        "rmsd_ub_bestpose": rmsd_ub_bp,
        "max_affinity": aff_max,
        "lb_affinity": lb_aff,
        "ub_affinity": ub_aff,
    }


def parse_gnina_log(log_path: str) -> Optional[Dict[str, float]]:
    """Разбор одного gnina .log файла."""
    try:
        with open(log_path, "r") as f:
            lines = f.readlines()
    except OSError:
        return None

    data_rows: List[Tuple[int, float, float, float, float]] = []

    for line in lines:
        line = line.strip()
        if re.match(r"^\d+\s+[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+$", line):
            parts = re.split(r"\s+", line)
            if len(parts) != 5:
                continue
            try:
                mode = int(parts[0])
                affinity = float(parts[1])
                intramol = float(parts[2])
                cnn_pose = float(parts[3])
                cnn_affinity = float(parts[4])
            except ValueError:
                continue
            data_rows.append((mode, affinity, intramol, cnn_pose, cnn_affinity))

    if not data_rows:
        return None

    min_aff_row = min(data_rows, key=lambda x: x[1])
    bestpose_row = max(data_rows, key=lambda x: x[3])

    _, aff_min, intr_min, cnn_pose_min, cnn_aff_min = min_aff_row
    _, aff_bp, intr_bp, cnn_pose_bp, cnn_aff_bp = bestpose_row

    return {
        "affinity_min": aff_min,
        "intramol_min": intr_min,
        "cnn_pose_score_min": cnn_pose_min,
        "cnn_affinity_min": cnn_aff_min,
        "affinity_bestpose": aff_bp,
        "intramol_bestpose": intr_bp,
        "cnn_pose_score_bestpose": cnn_pose_bp,
        "cnn_affinity_bestpose": cnn_aff_bp,
    }


def parse_plapt_json(json_path: str) -> Optional[Dict[str, float]]:
    """Разбор одного plapt .json файла."""
    try:
        with open(json_path, "r") as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError):
        return None

    result: Dict[str, float] = {}
    if "affinity" in data:
        result["affinity"] = float(data["affinity"])
    if "affinity_uM" in data:
        result["affinity_uM"] = float(data["affinity_uM"])
    return result if result else None


def parse_dynamicbind_filename(fname: str) -> Optional[Tuple[float, float]]:
    """Разбор имени sdf-файла DynamicBind."""
    m = re.search(r"_lddt([0-9.]+)_affinity([0-9.\-]+)\.sdf$", fname)
    if not m:
        return None
    try:
        lddt = float(m.group(1))
        aff = float(m.group(2))
    except ValueError:
        return None
    return lddt, aff


def dynamicbind_ligand_from_folder(folder_name: str) -> str:
    """По имени папки index470_idx_470 делаем идентификатор лиганда idx_470."""
    m = re.search(r"index(\d+)_idx_(\d+)", folder_name)
    if m:
        return f"idx_{m.group(2)}"
    return folder_name


def dynamicbind_ligand_from_csv_row(row: pd.Series) -> Optional[str]:
    """Return the explicit source ligand ID from a DynamicBind summary row.

    DynamicBind can omit failed ligands from affinity_prediction.csv while
    preserving their original ``idx_N`` names.  Therefore the DataFrame row
    number is not a ligand identifier and must never be used as one.
    """
    for column in ("name", "ligand_id", "ligand", "id"):
        if column not in row.index or pd.isna(row[column]):
            continue
        value = str(row[column]).strip()
        idx_match = re.fullmatch(r"idx_(\d+)", value)
        if idx_match:
            return f"idx_{int(idx_match.group(1))}"
        ligand_match = re.fullmatch(r"ligand_(\d+)", value)
        if ligand_match:
            return format_ligand_id(int(ligand_match.group(1)))
        if re.fullmatch(r"\d+", value):
            return f"idx_{int(value)}"
    return None


def collect_metrics_from_docking_dir(
    protein: str, docking_dir: Path, dynamicbind_subdir: str = "dynamicbind"
) -> Dict[str, Dict[str, float]]:
    """
    Собирает метрики всех методов докинга из директории docking.
    
    Returns:
        Словарь: {ligand_id: {method_metric: value, ...}, ...}
        где ligand_id = "ligand_1", "ligand_2", ... или "idx_0", "idx_1", ...
    """
    from typing import Tuple
    
    all_metrics: Dict[str, Dict[str, float]] = {}
    
    # QVina
    qvina_dir = docking_dir / "qvina"
    if qvina_dir.exists():
        print(f"   Собираю метрики qvina...")
        count = 0
        for fname in os.listdir(qvina_dir):
            if not fname.endswith(".log") or not fname.startswith("ligand_"):
                continue
            ligand = os.path.splitext(fname)[0]
            log_path = qvina_dir / fname
            parsed = parse_qvina_log(str(log_path))
            if parsed:
                if ligand not in all_metrics:
                    all_metrics[ligand] = {}
                for k, v in parsed.items():
                    all_metrics[ligand][f"qvina_{k}"] = v
                count += 1
        print(f"     Найдено {count} метрик qvina")
    
    # GNina
    gnina_dir = docking_dir / "gnina"
    if gnina_dir.exists():
        print(f"   Собираю метрики gnina...")
        count = 0
        for fname in os.listdir(gnina_dir):
            if not fname.endswith(".log") or not fname.startswith("ligand_"):
                continue
            ligand = os.path.splitext(fname)[0]
            log_path = gnina_dir / fname
            parsed = parse_gnina_log(str(log_path))
            if parsed:
                if ligand not in all_metrics:
                    all_metrics[ligand] = {}
                for k, v in parsed.items():
                    all_metrics[ligand][f"gnina_{k}"] = v
                count += 1
        print(f"     Найдено {count} метрик gnina")
    
    # PLAPT
    plapt_dir = docking_dir / "plapt"
    if plapt_dir.exists():
        print(f"   Собираю метрики plapt...")
        count = 0
        for sub in os.listdir(plapt_dir):
            sub_dir = plapt_dir / sub
            if not sub_dir.is_dir():
                continue
            for fname in os.listdir(sub_dir):
                if not fname.endswith(".json") or not fname.startswith("ligand_"):
                    continue
                ligand = os.path.splitext(fname)[0]
                jpath = sub_dir / fname
                parsed = parse_plapt_json(str(jpath))
                if parsed:
                    if ligand not in all_metrics:
                        all_metrics[ligand] = {}
                    for k, v in parsed.items():
                        all_metrics[ligand][f"plapt_{k}"] = v
                    count += 1
        print(f"     Найдено {count} метрик plapt")
    
    # DynamicBind
    dynamicbind_dir = docking_dir / dynamicbind_subdir
    if dynamicbind_dir.exists():
        print(f"   Собираю метрики dynamicbind...")
        count = 0
        # Some result trees contain a stale nested copy plus a newer direct
        # rerun (notably 4tz4).  Prefer the shallowest summary deterministically;
        # never combine multiple panels with the same idx_N identifiers.
        summary_candidates = sorted(dynamicbind_dir.glob("**/affinity_prediction.csv"))
        if summary_candidates:
            depths = {
                path: len(path.relative_to(dynamicbind_dir).parts)
                for path in summary_candidates
            }
            min_depth = min(depths.values())
            preferred = [path for path in summary_candidates if depths[path] == min_depth]
            if len(preferred) > 1:
                raise RuntimeError(
                    "Ambiguous DynamicBind summaries at equal depth: "
                    + ", ".join(str(path) for path in preferred)
                )
            csv_path = preferred[0]
            ignored = [path for path in summary_candidates if path != csv_path]
            print(f"     DynamicBind summary: {csv_path}")
            if ignored:
                print(
                    "     WARNING: ignored deeper stale DynamicBind summaries: "
                    + ", ".join(str(path) for path in ignored)
                )
            df_csv = pd.read_csv(csv_path)
            aff_col = next(
                (
                    column
                    for column in ("affinity", "predicted_affinity", "affinity_pred")
                    if column in df_csv.columns
                ),
                None,
            )
            if aff_col is None:
                raise ValueError(f"No affinity column in DynamicBind summary {csv_path}")
            missing_explicit_id = 0
            duplicate_ids: set[str] = set()
            seen_ids: set[str] = set()
            for _, row in df_csv.iterrows():
                try:
                    aff = float(row[aff_col])
                except (TypeError, ValueError):
                    continue
                ligand = dynamicbind_ligand_from_csv_row(row)
                if ligand is None:
                    missing_explicit_id += 1
                    continue
                if ligand in seen_ids:
                    duplicate_ids.add(ligand)
                    continue
                seen_ids.add(ligand)
                if ligand not in all_metrics:
                    all_metrics[ligand] = {}
                all_metrics[ligand]["dynamicbind_affinity_maxaff"] = aff
                all_metrics[ligand]["dynamicbind_affinity_bestpose"] = aff
                all_metrics[ligand]["dynamicbind_lddt_maxaff"] = float("nan")
                all_metrics[ligand]["dynamicbind_lddt_bestpose"] = float("nan")
                count += 1
            if missing_explicit_id:
                raise ValueError(
                    f"{missing_explicit_id} DynamicBind rows lack explicit ligand IDs "
                    f"in {csv_path}"
                )
            if duplicate_ids:
                raise ValueError(
                    "Duplicate DynamicBind ligand IDs in summary: "
                    + ", ".join(sorted(duplicate_ids))
                )
        print(f"     Найдено {count} метрик dynamicbind")
    
    return all_metrics


def _resolve_canonical_ligand_id(
    native: str, lookup: dict[str, str]
) -> Optional[str]:
    """Map an engine filename/folder token onto ligand_0001. Never use glob order."""
    if native in lookup:
        return lookup[native]
    upper = native.upper()
    if upper in lookup:
        return lookup[upper]
    padded = normalize_ligand_id(native)
    if padded and padded in lookup:
        return padded
    number = parse_ligand_number(native)
    if number is not None:
        candidate = format_ligand_id(number)
        if candidate in lookup:
            return candidate
        return candidate
    db_number = parse_dynamicbind_idx(native)
    if db_number is not None:
        candidate = format_ligand_id(db_number)
        if candidate in lookup:
            return candidate
        return candidate
    return None


def _remap_metrics_to_ligand_id(
    metrics_by_native: Dict[str, Dict[str, float]], lookup: dict[str, str]
) -> Dict[str, Dict[str, float]]:
    remapped: Dict[str, Dict[str, float]] = {}
    unmatched: list[str] = []
    for native, metrics in metrics_by_native.items():
        ligand_id = _resolve_canonical_ligand_id(native, lookup)
        if ligand_id is None:
            unmatched.append(native)
            continue
        remapped.setdefault(ligand_id, {}).update(metrics)
    if unmatched:
        print(
            f"   WARNING: {len(unmatched)} engine tokens not in ligand_id map "
            f"(first 8: {unmatched[:8]})"
        )
    return remapped


def _load_ligand_map(base_dir: Path, protein: str) -> pd.DataFrame:
    path = base_dir / "processed" / "id_maps" / f"{protein}_ligand_id_map.csv"
    if path.is_file():
        return pd.read_csv(path)
    nodubl = base_dir / "input" / "ligands_nodubl" / nodubl_csv_name(protein)
    from analysis.ligand_identity import map_rows_from_nodubl

    return map_rows_from_nodubl(protein, nodubl)


def merge_ligands_with_docking(
    ligands_csv: Path,
    protein: str,
    base_dir: Path,
    boltz_results_dir: str = BOLTZ_RESULTS_DIR,
    use_metrics_csv: bool = True,
    dynamicbind_subdir: str = "dynamicbind",
    results_dir: Path | None = None,
    example_n: int | None = None,
) -> pd.DataFrame:
    """
    Объединяет данные лигандов с метриками докинга (boltz2 и exp_*).
    
    Args:
        ligands_csv: Путь к CSV файлу из ligands_nodubl
        protein: PDB ID белка
        base_dir: Базовая директория проекта
        boltz_results_dir: Директория с результатами Boltz
    
    Returns:
        Объединенный DataFrame
    """
    print(f"📖 Читаю лиганды из: {ligands_csv.name}")

    try:
        sample = ligands_csv.read_text(encoding="utf-8", errors="replace")[:2048]
        sep = ";" if sample.count(";") >= sample.count(",") else ","
        ligands_df = pd.read_csv(ligands_csv, sep=sep, low_memory=False)
    except Exception as e:
        print(f"❌ Ошибка при чтении {ligands_csv}: {e}")
        sys.exit(1)

    missing_cols = [col for col in LIGANDS_COLUMNS if col not in ligands_df.columns]
    if missing_cols:
        print(f"⚠️  Отсутствуют столбцы в ligands CSV: {missing_cols}")
        available_cols = [col for col in LIGANDS_COLUMNS if col in ligands_df.columns]
        ligands_selected = ligands_df[available_cols].copy()
    else:
        ligands_selected = ligands_df[LIGANDS_COLUMNS].copy()

    if "molecule_chembl_id" not in ligands_selected.columns and "molecule_chembl_id" in ligands_df.columns:
        ligands_selected["molecule_chembl_id"] = ligands_df["molecule_chembl_id"]

    id_map = _load_ligand_map(base_dir, protein)
    lookup = native_to_ligand_id(id_map)
    if "ligand_id" not in ligands_selected.columns:
        ligands_selected.insert(
            0,
            "ligand_id",
            [format_ligand_id(i) for i in range(1, len(ligands_selected) + 1)],
        )
    else:
        ligands_selected["ligand_id"] = ligands_selected["ligand_id"].map(
            lambda x: normalize_ligand_id(str(x)) or x
        )

    if example_n:
        keep = {format_ligand_id(i) for i in range(1, int(example_n) + 1)}
        ligands_selected = ligands_selected[ligands_selected["ligand_id"].isin(keep)].copy()
        print(f"   EXAMPLE_N={example_n}: keeping {len(ligands_selected)} ligand_id rows")

    print(f"   Найдено {len(ligands_selected)} записей лигандов (assay-level / nodubl)")

    # Pre-dock collapse: mean experimental Ki per SMILES; keep existing ligand_id.
    # Docking artifacts remain keyed by that representative ligand_id (no renumbering).
    n_before = len(ligands_selected)
    ligands_selected = group_by_smiles(ligands_selected)
    print(
        f"   Pre-dock mean-Ki collapse by SMILES: {n_before} → {len(ligands_selected)} "
        f"(ligand_id preserved from each group)"
    )

    if results_dir is None:
        results_dir = base_dir / "results"
    docking_dir = results_dir / protein / "docking"
    print("🔍 Собираю метрики докинга из исходных файлов...")
    all_metrics_by_native = collect_metrics_from_docking_dir(
        protein, docking_dir, dynamicbind_subdir=dynamicbind_subdir
    )
    all_metrics_by_ligand = _remap_metrics_to_ligand_id(all_metrics_by_native, lookup)
    # Restrict to surviving ligand_ids so left-merge cannot pick up orphan scores
    keep_ids = set(ligands_selected["ligand_id"].astype(str))
    all_metrics_by_ligand = {
        lid: m for lid, m in all_metrics_by_ligand.items() if str(lid) in keep_ids
    }
    print(f"   Найдено метрик для {len(all_metrics_by_ligand)} ligand_id (после collapse)")

    print(f"🔍 Собираю метрики Boltz для {protein}...")
    boltz_metrics = collect_boltz_metrics(protein, boltz_results_dir)
    boltz_by_ligand: Dict[str, Dict[str, float]] = {}
    chembl_lookup = chembl_to_ligand_ids(id_map)
    unmatched_boltz = []
    for chembl_id, metrics in boltz_metrics.items():
        lids = chembl_lookup.get(chembl_id) or chembl_lookup.get(str(chembl_id).upper(), [])
        if not lids:
            unmatched_boltz.append(chembl_id)
            continue
        for lid in lids:
            if str(lid) not in keep_ids:
                continue
            boltz_by_ligand.setdefault(lid, {}).update(metrics)
    if unmatched_boltz:
        print(
            f"   WARNING: {len(unmatched_boltz)} Boltz CHEMBL ids not in map "
            f"(first 8: {unmatched_boltz[:8]})"
        )
    print(f"   Найдено {len(boltz_metrics)} Boltz native / {len(boltz_by_ligand)} ligand_id")

    ligands_selected["exp_value"] = ligands_selected.get("value", None)
    ligands_selected["exp_standard_value"] = ligands_selected.get("standard_value", None)
    ligands_selected["exp_pchembl_value"] = ligands_selected.get("pchembl_value", None)

    print("🔗 Объединяю данные по (pdb_id, ligand_id) — только представительские id...")
    merged_df = ligands_selected.copy()
    if all_metrics_by_ligand:
        metrics_df = pd.DataFrame(
            [{"ligand_id": lid, **metrics} for lid, metrics in all_metrics_by_ligand.items()]
        )
        merged_df = merged_df.merge(metrics_df, on="ligand_id", how="left", validate="one_to_one")
    if boltz_by_ligand:
        boltz_df = pd.DataFrame(
            [{"ligand_id": lid, **metrics} for lid, metrics in boltz_by_ligand.items()]
        )
        merged_df = merged_df.merge(boltz_df, on="ligand_id", how="left", validate="one_to_one")

    merged_df["protein_pdb_id"] = protein.upper()
    leading = ["protein_pdb_id", "ligand_id"]
    rest = [c for c in merged_df.columns if c not in leading]
    merged_df = merged_df[leading + rest]
    
    # Статистика объединения
    total_ligands = len(merged_df)
    boltz_cols_present = [col for col in BOLTZ_COLUMNS if col in merged_df.columns]
    if boltz_cols_present:
        has_boltz_data = merged_df[boltz_cols_present[0]].notna().sum()
    else:
        has_boltz_data = 0
    
    print(f"   Всего лигандов: {total_ligands}")
    print(f"   Лигандов с данными Boltz: {has_boltz_data}")
    print(f"   Лигандов без данных Boltz: {total_ligands - has_boltz_data}")
    for col, label in (
        ("qvina_affinity_bestpose", "qvina"),
        ("gnina_cnn_affinity_bestpose", "gnina"),
        ("plapt_affinity", "plapt"),
        ("dynamicbind_affinity_bestpose", "dynamicbind"),
    ):
        if col in merged_df.columns:
            print(f"   {label} scored: {int(merged_df[col].notna().sum())}")

    return merged_df


def process_single_protein(
    pdb_id: str,
    base_dir: Path,
    output_dir: Optional[Path] = None,
    boltz_results_dir: str = BOLTZ_RESULTS_DIR,
    dynamicbind_subdir: str = "dynamicbind_new",
    results_dir: Path | None = None,
    example_n: int | None = None,
) -> bool:
    """
    Обрабатывает один белок: объединяет лиганды с метриками докинга.
    
    Args:
        pdb_id: PDB ID белка (например, "3eyg")
        base_dir: Базовая директория проекта
        output_dir: Директория для сохранения результатов
        boltz_results_dir: Директория с результатами Boltz
    
    Returns:
        True если успешно, False иначе
    """
    pdb_id_lower = pdb_id.lower()
    if pdb_id_lower not in LEGACY_PDB_TO_NODUBL_STEM:
        print(f"❌ PDB ID {pdb_id} не найден в маппинге")
        return False

    curated = base_dir / "input" / "ligands_curated" / f"{pdb_id_lower}_ligands.csv"
    ligands_csv = (
        curated
        if curated.is_file()
        else base_dir / "input" / "ligands_nodubl" / nodubl_csv_name(pdb_id_lower)
    )
    if not ligands_csv.exists():
        print(f"❌ Файл лигандов не найден: {ligands_csv}")
        return False

    if output_dir is None:
        output_dir = base_dir / "analysis" / "excluding_2z5x_3mjg" / "tables"
    
    # Создаем выходной файл
    output_csv = output_dir / f"merged_ligands_docking_{pdb_id_lower}.csv"
    
    # Объединяем данные
    try:
        merged_df = merge_ligands_with_docking(
            ligands_csv,
            pdb_id_lower,
            base_dir,
            boltz_results_dir,
            dynamicbind_subdir=dynamicbind_subdir,
            results_dir=results_dir,
            example_n=example_n,
        )
        
        # Сохраняем результат
        print(f"💾 Сохраняю результат в: {output_csv}")
        merged_df.to_csv(output_csv, index=False, sep=",")
        print(f"✅ Готово! Сохранено {len(merged_df)} записей")
        return True
    except Exception as e:
        print(f"❌ Ошибка при обработке {pdb_id}: {e}")
        import traceback
        traceback.print_exc()
        return False


def process_all_proteins(
    base_dir: Path,
    output_dir: Optional[Path] = None,
    boltz_results_dir: str = BOLTZ_RESULTS_DIR,
    dynamicbind_subdir: str = "dynamicbind_new",
    results_dir: Path | None = None,
    example_n: int | None = None,
) -> None:
    """
    Обрабатывает все белки из маппинга.
    
    Args:
        base_dir: Базовая директория проекта
        output_dir: Директория для сохранения результатов
        boltz_results_dir: Директория с результатами Boltz
    """
    print("=" * 80)
    print("ОБЪЕДИНЕНИЕ ДАННЫХ ЛИГАНДОВ С МЕТРИКАМИ ДОКИНГА (ИЗ ДИРЕКТОРИИ DOCKING)")
    print("=" * 80)
    print()
    
    if output_dir is None:
        output_dir = base_dir / "analysis" / "excluding_2z5x_3mjg" / "tables"

    success_count = 0
    fail_count = 0

    for pdb_id in CANONICAL_TARGETS:
        print(f"\n{'=' * 80}")
        print(f"Обработка белка: {pdb_id.upper()}")
        print(f"{'=' * 80}")
        
        if process_single_protein(
            pdb_id,
            base_dir,
            output_dir,
            boltz_results_dir,
            dynamicbind_subdir=dynamicbind_subdir,
            results_dir=results_dir,
            example_n=example_n,
        ):
            success_count += 1
        else:
            fail_count += 1
    
    print(f"\n{'=' * 80}")
    print("ИТОГИ")
    print(f"{'=' * 80}")
    print(f"✅ Успешно обработано: {success_count}")
    print(f"❌ Ошибок: {fail_count}")
    print(f"📁 Результаты сохранены в: {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Объединяет данные лигандов из ligands_nodubl с метриками докинга из директории docking"
    )
    parser.add_argument(
        "pdb_id",
        nargs="?",
        type=str,
        help="PDB ID белка для обработки (например, 3eyg). Если не указан, обрабатываются все белки.",
    )
    parser.add_argument(
        "--base-dir",
        type=str,
        default=None,
        help="Базовая директория проекта (по умолчанию: директория скрипта)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Директория для сохранения результатов (по умолчанию: analysis/tables)",
    )
    parser.add_argument(
        "--results-dir",
        type=str,
        default=None,
        help="Корневая директория сырых результатов докинга (по умолчанию: <base-dir>/results)",
    )
    parser.add_argument(
        "--boltz-results-dir",
        type=str,
        default=BOLTZ_RESULTS_DIR or None,
        help="Директория с результатами Boltz (env: BOLTZ_RESULTS_DIR)",
    )
    parser.add_argument(
        "--dynamicbind-subdir",
        type=str,
        default="dynamicbind_new",
        help="Подпапка DynamicBind внутри docking (по умолчанию: dynamicbind_new)",
    )
    parser.add_argument(
        "--example-n",
        type=int,
        default=None,
        help="Keep only ligand_0001..ligand_000N (GitHub example). Default: all rows.",
    )

    args = parser.parse_args()
    
    if args.base_dir:
        base_dir = Path(args.base_dir)
    else:
        base_dir = Path(__file__).resolve().parent.parent.parent
    
    output_dir = Path(args.output_dir) if args.output_dir else None
    results_dir = Path(args.results_dir) if args.results_dir else None
    boltz_dir = args.boltz_results_dir or os.environ.get("BOLTZ_RESULTS_DIR", "")
    example_n = args.example_n
    if example_n is None and os.environ.get("EXAMPLE_N"):
        example_n = int(os.environ["EXAMPLE_N"])

    if args.pdb_id:
        success = process_single_protein(
            args.pdb_id,
            base_dir,
            output_dir,
            boltz_dir,
            dynamicbind_subdir=args.dynamicbind_subdir,
            results_dir=results_dir,
            example_n=example_n,
        )
        sys.exit(0 if success else 1)
    else:
        process_all_proteins(
            base_dir,
            output_dir,
            boltz_dir,
            dynamicbind_subdir=args.dynamicbind_subdir,
            results_dir=results_dir,
            example_n=example_n,
        )


if __name__ == "__main__":
    main()

