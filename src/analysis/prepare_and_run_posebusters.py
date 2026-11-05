#!/usr/bin/env python3
"""
Prepare PoseBusters inputs from docking results and optionally run checks.

Pipeline:
1) Collect poses from docking outputs:
   - gnina/qvina: convert ligand_*_out.pdbqt -> top-1 SDF via Open Babel
   - dynamicbind: rank1 ligand SDF gets canonical SMILES bond orders, paired
     with the rank1 receptor PDB from the same index*_idx_* folder
2) Build a PoseBusters table CSV with columns:
   mol_pred,mol_cond,protein,method,ligand
3) Optionally run:
   bust -t <table.csv> --outfmt csv --output <results.csv>

Notes:
- PoseBusters supports ligand inputs in .sdf/.mol/.mol2. .pdbqt is converted here.
- This script runs "dock-like" checks (mol_pred + mol_cond). It does not require mol_true.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    from .posebusters_runner import run_posebusters_chunked
    from .pose_chemistry import (
        TEMPLATE_MARKER,
        pdbqt_to_template_mol,
        sdf_to_template_mol,
        write_pose_molecule,
    )
except ImportError:
    from posebusters_runner import run_posebusters_chunked
    from pose_chemistry import (
        TEMPLATE_MARKER,
        pdbqt_to_template_mol,
        sdf_to_template_mol,
        write_pose_molecule,
    )


PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_POSES_DIR = PROJECT_ROOT / "results"
DEFAULT_PROTEINS_DIR = PROJECT_ROOT / "processed" / "proteins"
DEFAULT_OUTPUT_DIR = (
    PROJECT_ROOT / "analysis" / "excluding_2z5x_3mjg" / "tables" / "posebuster"
)
PDB_TO_LIGAND_CSV = {
    "1g5m": "BCL2_Ki_WT_ChEMBL_252_nodubl.csv",
    "2v5z": "MAO-B_Ki_WT_ChEMBL_246_nodubl.csv",
    "2z5x": "MAO-B_Ki_WT_ChEMBL_246_nodubl.csv",
    "3eyg": "JAK1_Ki_WT_ChEMBL_2255_nodubl.csv",
    "3jy9": "JAK2_Ki_WT_ChEMBL_2027_nodubl.csv",
    "3lxk": "JAK3_Ki_WT_ChEMBL_786_nodubl.csv",
    "3mjg": "PDGFRB_Ki_WT_ChEMBL_275_nodubl.csv",
    "4ase": "VEGFR2_Ki_WT_ChEMBL_875_nodubl.csv",
    "4f65": "FGFR1_Ki_WT_ChEMBL_134_nodubl.csv",
    "4tz4": "CRBN_Ki_WT_ChEMBL_127_nodubl.csv",
    "4zau": "EGFR_Ki_WT_curated_251_nodubl.csv",
    "5jkv": "CYP19A1_Aromatase_Ki_WT_ChEMBL_548_nodubl.csv",
    "5mo4": "ABL1_BCR-ABL_Ki_WT_ChEMBL_693_nodubl.csv",
    "6gqj": "KIT_Ki_WT_curated_1298_nodubl.csv",
    "6jok": "PDGFRA_Ki_WT_curated_250_nodubl.csv",
    "5lf3": "PSMB5_Ki_WT_ChEMBL_88_nodubl.csv",
    "7awe": "PSMB5_Ki_WT_ChEMBL_88_nodubl.csv",
    "7kk3": "PARP1_Ki_WT_ChEMBL_1075_nodubl.csv",
    "11ue": "PDGFRB_Ki_WT_ChEMBL_275_nodubl.csv",
}


def _is_valid_sdf_file(sdf_path: Path) -> bool:
    """Quick integrity check for cached SDF files."""
    if not sdf_path.exists() or sdf_path.stat().st_size == 0:
        return False
    try:
        # SDF record separator; enough for a lightweight sanity check.
        with sdf_path.open("r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
        return "$$$$" in content and "M  END" in content
    except OSError:
        return False


def _find_protein_pdb(proteins_dir: Path, protein: str) -> Optional[Path]:
    """Find protein pdb by cleaned-chainA and fallback names."""
    candidates = [
        proteins_dir / f"{protein.lower()}_chainA.pdb",
        proteins_dir / f"{protein.upper()}_chainA.pdb",
        proteins_dir / f"{protein.lower()}_chainK.pdb",
        proteins_dir / f"{protein.upper()}_chainK.pdb",
        proteins_dir / f"{protein.lower()}.pdb",
        proteins_dir / f"{protein.upper()}.pdb",
    ]
    for path in candidates:
        if path.exists():
            return path
    return None


def _discover_proteins(poses_dir: Path) -> List[str]:
    proteins: List[str] = []
    for sub in poses_dir.iterdir():
        if not sub.is_dir():
            continue
        if (sub / "docking").is_dir():
            proteins.append(sub.name.lower())
    return sorted(proteins)


def _load_smiles_by_ligand(protein: str) -> Dict[str, str]:
    from docking_benchmark2.ligand_ids import (
        format_ligand_id,
        parse_ligand_number,
        vina_legacy_id,
    )

    curated = PROJECT_ROOT / "input" / "ligands_curated" / f"{protein.lower()}_ligands.csv"
    csv_name = PDB_TO_LIGAND_CSV.get(protein.lower())
    csv_path = curated if curated.is_file() else (
        PROJECT_ROOT / "input" / "ligands_nodubl" / csv_name if csv_name else None
    )
    if not csv_path or not csv_path.is_file():
        return {}
    sample = csv_path.read_text(encoding="utf-8", errors="replace")[:2048]
    delimiter = ";" if sample.count(";") >= sample.count(",") else ","
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter=delimiter))
    out: Dict[str, str] = {}
    for row_index, row in enumerate(rows, start=1):
        smiles = row.get("canonical_smiles") or row.get("smiles")
        if not smiles:
            continue
        ligand_id = row.get("ligand_id") or format_ligand_id(row_index)
        number = parse_ligand_number(ligand_id)
        if number is None:
            continue
        # Key only by the ligand's own number. Row order is not the ligand number
        # after duplicate collapse, and must not overwrite another ligand's SMILES.
        out[format_ligand_id(number)] = smiles
        out[vina_legacy_id(number)] = smiles
    return out


def _is_template_sdf_file(sdf_path: Path) -> bool:
    if not _is_valid_sdf_file(sdf_path):
        return False
    try:
        from rdkit import Chem

        supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        mol = supplier[0] if supplier and len(supplier) else None
        return bool(
            mol is not None
            and mol.HasProp(TEMPLATE_MARKER)
            and mol.GetProp(TEMPLATE_MARKER) == "1"
        )
    except Exception:
        return False


def _convert_pdbqt_to_sdf_top1(
    pdbqt_path: Path,
    sdf_path: Path,
    canonical_smiles: str,
    force: bool,
) -> bool:
    """Transfer canonical chemistry onto top-pose PDBQT coordinates."""
    if sdf_path.exists() and not force and _is_template_sdf_file(sdf_path):
        return True
    sdf_path.unlink(missing_ok=True)
    sdf_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        fixed = pdbqt_to_template_mol(pdbqt_path, canonical_smiles)
        write_pose_molecule(fixed, sdf_path)
    except Exception as exc:
        print(f"Template chemistry conversion failed [{pdbqt_path}]: {exc}")
        sdf_path.unlink(missing_ok=True)
        return False
    return _is_template_sdf_file(sdf_path)


def _convert_dynamicbind_sdf(
    src_sdf: Path,
    sdf_path: Path,
    canonical_smiles: str,
    force: bool,
) -> bool:
    """Transfer canonical chemistry onto a DynamicBind rank-1 SDF."""
    if sdf_path.exists() and not force and _is_template_sdf_file(sdf_path):
        return True
    sdf_path.unlink(missing_ok=True)
    sdf_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        fixed = sdf_to_template_mol(src_sdf, canonical_smiles)
        write_pose_molecule(fixed, sdf_path)
    except Exception as exc:
        print(f"Template chemistry conversion failed [{src_sdf}]: {exc}")
        sdf_path.unlink(missing_ok=True)
        return False
    return _is_template_sdf_file(sdf_path)


def _parse_dynamicbind_ligand_id(index_dir_name: str) -> str:
    """index470_idx_470 is 0-based → ligand_0471. Fallback to folder name."""
    from docking_benchmark2.ligand_ids import format_ligand_id

    m = re.search(r"index(\d+)_idx_(\d+)", index_dir_name)
    if m:
        return format_ligand_id(int(m.group(2)) + 1)
    return index_dir_name


def _collect_dynamicbind_rank1_pair(index_dir: Path) -> Optional[Tuple[Path, Path]]:
    """Pick rank1 ligand+receptor pair; fallback to first valid rank* pair."""
    rank1 = sorted(index_dir.glob("rank1_ligand_*.sdf"))
    if rank1:
        ligand = rank1[0]
        receptor = index_dir / ligand.name.replace("_ligand_", "_receptor_").replace(".sdf", ".pdb")
        if receptor.exists():
            return ligand, receptor
    rank_any = sorted(index_dir.glob("rank*_ligand_*.sdf"))
    for ligand in rank_any:
        receptor = index_dir / ligand.name.replace("_ligand_", "_receptor_").replace(".sdf", ".pdb")
        if receptor.exists():
            return ligand, receptor
    return None


def _prepare_rows_for_protein(
    protein: str,
    poses_dir: Path,
    proteins_dir: Path,
    output_sdf_dir: Path,
    methods: List[str],
    obabel_bin: str,
    force: bool,
    dynamicbind_subdir: str,
    dynamicbind_method_label: str,
) -> Tuple[Dict[str, List[Dict[str, str]]], Dict[str, int]]:
    rows_by_method: Dict[str, List[Dict[str, str]]] = {"gnina": [], "qvina": [], "dynamicbind": []}
    stats = {"gnina": 0, "qvina": 0, "dynamicbind": 0, "failed_conversions": 0}
    smiles_by_ligand = _load_smiles_by_ligand(protein)

    docking_dir = poses_dir / protein / "docking"
    if not docking_dir.is_dir():
        return rows_by_method, stats
    protein_pdb = _find_protein_pdb(proteins_dir, protein)

    if "gnina" in methods:
        if protein_pdb is None:
            print(f"Skip GNINA [{protein.upper()}]: protein PDB not found in {proteins_dir}")
        else:
            gnina_dir = docking_dir / "gnina"
            for pdbqt_path in sorted(gnina_dir.glob("ligand_*_out.pdbqt")):
                ligand = pdbqt_path.stem.replace("_out", "")
                sdf_path = output_sdf_dir / protein / "gnina" / f"{ligand}.sdf"
                smiles = smiles_by_ligand.get(ligand)
                ok = bool(
                    smiles
                    and _convert_pdbqt_to_sdf_top1(
                        pdbqt_path, sdf_path, smiles, force
                    )
                )
                if not ok:
                    stats["failed_conversions"] += 1
                    continue
                rows_by_method["gnina"].append(
                    {
                        "mol_pred": str(sdf_path),
                        "mol_cond": str(protein_pdb),
                        "protein": protein.upper(),
                        "method": "gnina",
                        "ligand": ligand,
                    }
                )
                stats["gnina"] += 1

    if "qvina" in methods:
        if protein_pdb is None:
            print(f"Skip QVINA [{protein.upper()}]: protein PDB not found in {proteins_dir}")
        else:
            qvina_dir = docking_dir / "qvina"
            for pdbqt_path in sorted(qvina_dir.glob("ligand_*_out.pdbqt")):
                ligand = pdbqt_path.stem.replace("_out", "")
                sdf_path = output_sdf_dir / protein / "qvina" / f"{ligand}.sdf"
                smiles = smiles_by_ligand.get(ligand)
                ok = bool(
                    smiles
                    and _convert_pdbqt_to_sdf_top1(
                        pdbqt_path, sdf_path, smiles, force
                    )
                )
                if not ok:
                    stats["failed_conversions"] += 1
                    continue
                rows_by_method["qvina"].append(
                    {
                        "mol_pred": str(sdf_path),
                        "mol_cond": str(protein_pdb),
                        "protein": protein.upper(),
                        "method": "qvina",
                        "ligand": ligand,
                    }
                )
                stats["qvina"] += 1

    if "dynamicbind" in methods:
        dynamicbind_dir = docking_dir / dynamicbind_subdir
        # DynamicBind outputs occur in both flat and nested dataset layouts.
        candidates: Dict[str, Path] = {}
        for index_dir in sorted(dynamicbind_dir.glob("**/index*_idx_*")):
            if not index_dir.is_dir():
                continue
            chosen_pair = _collect_dynamicbind_rank1_pair(index_dir)
            if chosen_pair is None:
                continue
            ligand = _parse_dynamicbind_ligand_id(index_dir.name)
            previous = candidates.get(ligand)
            if previous is None or len(index_dir.relative_to(dynamicbind_dir).parts) < len(
                previous.relative_to(dynamicbind_dir).parts
            ):
                candidates[ligand] = index_dir

        for ligand, index_dir in sorted(
            candidates.items(),
            key=lambda item: int(item[0].split("_", 1)[1]),
        ):
            chosen_pair = _collect_dynamicbind_rank1_pair(index_dir)
            if chosen_pair is None:
                continue
            chosen_ligand, chosen_receptor = chosen_pair
            sdf_path = output_sdf_dir / protein / dynamicbind_method_label / f"{ligand}.sdf"
            smiles = smiles_by_ligand.get(ligand)
            ok = bool(
                smiles
                and _convert_dynamicbind_sdf(chosen_ligand, sdf_path, smiles, force)
            )
            if not ok:
                stats["failed_conversions"] += 1
                continue
            rows_by_method["dynamicbind"].append(
                {
                    "mol_pred": str(sdf_path),
                    # For dynamicbind use the receptor from the same rank folder.
                    "mol_cond": str(chosen_receptor),
                    "protein": protein.upper(),
                    "method": dynamicbind_method_label,
                    "ligand": ligand,
                }
            )
            stats["dynamicbind"] += 1

    return rows_by_method, stats


def _flatten_rows(rows_by_method: Dict[str, List[Dict[str, str]]]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for method_rows in rows_by_method.values():
        rows.extend(method_rows)
    return rows


def _prepare_rows_for_protein_legacy(
    protein: str,
    poses_dir: Path,
    proteins_dir: Path,
    output_sdf_dir: Path,
    methods: List[str],
    obabel_bin: str,
    force: bool,
    dynamicbind_subdir: str,
    dynamicbind_method_label: str,
) -> Tuple[List[Dict[str, str]], Dict[str, int]]:
    """Compatibility wrapper that returns a single row list."""
    rows_by_method, stats = _prepare_rows_for_protein(
        protein=protein,
        poses_dir=poses_dir,
        proteins_dir=proteins_dir,
        output_sdf_dir=output_sdf_dir,
        methods=methods,
        obabel_bin=obabel_bin,
        force=force,
        dynamicbind_subdir=dynamicbind_subdir,
        dynamicbind_method_label=dynamicbind_method_label,
    )
    return _flatten_rows(rows_by_method), stats


def _write_table(rows: List[Dict[str, str]], table_csv: Path) -> None:
    table_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["mol_pred", "mol_cond", "protein", "method", "ligand"]
    with table_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _run_posebusters(
    table_csv: Path,
    output_csv: Path,
    bust_bin: str,
    max_workers: int,
    posebusters_config: Optional[str] = None,
    chunk_size: int = 100,
) -> None:
    run_posebusters_chunked(
        table_csv=table_csv,
        output_csv=output_csv,
        bust_bin=bust_bin,
        max_workers=max_workers,
        posebusters_config=posebusters_config,
        chunk_size=chunk_size,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare SDF poses and run PoseBusters for docking-benchmark-2 results."
    )
    parser.add_argument(
        "--poses-dir",
        type=str,
        default=str(DEFAULT_POSES_DIR),
        help="Directory with per-protein docking outputs (default: data/results_bindingdb).",
    )
    parser.add_argument(
        "--proteins-dir",
        type=str,
        default=str(DEFAULT_PROTEINS_DIR),
        help="Directory with protein pdb files (default: data/input/proteins).",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory for prepared SDFs, table and PoseBusters outputs.",
    )
    parser.add_argument(
        "--methods",
        type=str,
        default="gnina,qvina,dynamicbind",
        help="Comma-separated methods to include (gnina,qvina,dynamicbind).",
    )
    parser.add_argument(
        "--proteins",
        type=str,
        default="",
        help="Optional comma-separated protein ids. If empty, auto-discover from poses-dir.",
    )
    parser.add_argument(
        "--obabel-bin",
        type=str,
        default="obabel",
        help="Path/name of Open Babel executable for pdbqt->sdf conversion.",
    )
    parser.add_argument(
        "--bust-bin",
        type=str,
        default="bust",
        help="Path/name of PoseBusters CLI executable.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=8,
        help="PoseBusters max workers when running checks.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=100,
        help="Rows per resumable PoseBusters chunk; use 0 to disable chunking.",
    )
    parser.add_argument(
        "--posebusters-config",
        type=str,
        default="",
        help=(
            "Optional PoseBusters config file passed to '--config' "
            "(e.g. .../posebusters/config/dock_fast.yml)."
        ),
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Only prepare SDF + table CSV(s); do not run PoseBusters.",
    )
    parser.add_argument(
        "--single-csv",
        action="store_true",
        help=(
            "Use legacy single combined CSV mode "
            "(posebusters_input_table.csv / posebusters_results.csv). "
            "By default, one protein+method = one CSV."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force reconvert/copy SDF files even if they already exist.",
    )
    parser.add_argument(
        "--dynamicbind-subdir",
        type=str,
        default="dynamicbind",
        help="DynamicBind subdirectory inside docking/ (default: dynamicbind).",
    )
    parser.add_argument(
        "--dynamicbind-method-label",
        type=str,
        default="dynamicbind",
        help="Method label for dynamicbind rows/output files (default: dynamicbind).",
    )

    args = parser.parse_args()
    poses_dir = Path(args.poses_dir)
    proteins_dir = Path(args.proteins_dir)
    output_dir = Path(args.output_dir)
    output_sdf_dir = output_dir / "sdf"
    table_csv = output_dir / "posebusters_input_table.csv"
    output_csv = output_dir / "posebusters_results.csv"
    per_method_tables_dir = output_dir / "tables_by_method"
    per_method_results_dir = output_dir / "results_by_method"

    methods = [m.strip().lower() for m in args.methods.split(",") if m.strip()]
    allowed = {"gnina", "qvina", "dynamicbind"}
    unknown = [m for m in methods if m not in allowed]
    if unknown:
        raise ValueError(f"Unknown methods: {unknown}. Allowed: {sorted(allowed)}")

    if args.proteins:
        proteins = [p.strip().lower() for p in args.proteins.split(",") if p.strip()]
    else:
        proteins = _discover_proteins(poses_dir)

    all_rows: List[Dict[str, str]] = []
    rows_by_protein_method: Dict[str, Dict[str, List[Dict[str, str]]]] = {}
    summary: Dict[str, int] = {"gnina": 0, "qvina": 0, "dynamicbind": 0, "failed_conversions": 0}

    print(f"Preparing PoseBusters inputs from: {poses_dir}")
    print(f"Protein structures from: {proteins_dir}")
    print(f"Methods: {methods}")
    print(f"Proteins: {len(proteins)}")

    for protein in proteins:
        rows_by_method, stats = _prepare_rows_for_protein(
            protein=protein,
            poses_dir=poses_dir,
            proteins_dir=proteins_dir,
            output_sdf_dir=output_sdf_dir,
            methods=methods,
            obabel_bin=args.obabel_bin,
            force=args.force,
            dynamicbind_subdir=args.dynamicbind_subdir,
            dynamicbind_method_label=args.dynamicbind_method_label,
        )
        rows_by_protein_method[protein] = rows_by_method
        all_rows.extend(_flatten_rows(rows_by_method))
        for k in summary:
            summary[k] += stats[k]

    print(f"Prepared rows: {len(all_rows)}")
    print(
        "By method: "
        f"gnina={summary['gnina']}, "
        f"qvina={summary['qvina']}, "
        f"dynamicbind={summary['dynamicbind']}"
    )
    print(f"Failed conversions: {summary['failed_conversions']}")

    # Default mode: one protein+method = one table/results CSV
    if not args.single_csv:
        per_method_tables_dir.mkdir(parents=True, exist_ok=True)
        if not args.prepare_only:
            per_method_results_dir.mkdir(parents=True, exist_ok=True)

        created = 0
        for protein in proteins:
            for method in methods:
                rows = rows_by_protein_method.get(protein, {}).get(method, [])
                if not rows:
                    continue
                created += 1
                method_tag = (
                    args.dynamicbind_method_label if method == "dynamicbind" else method
                )
                method_table = per_method_tables_dir / f"posebusters_input_{protein.lower()}_{method_tag}.csv"
                _write_table(rows, method_table)
                print(
                    f"Table [{protein.upper()} {method.upper()}]: "
                    f"{method_table} (rows={len(rows)})"
                )

                if args.prepare_only:
                    continue
                method_output = per_method_results_dir / f"posebusters_results_{protein.lower()}_{method_tag}.csv"
                _run_posebusters(
                    table_csv=method_table,
                    output_csv=method_output,
                    bust_bin=args.bust_bin,
                    max_workers=args.max_workers,
                    posebusters_config=(args.posebusters_config.strip() or None),
                    chunk_size=args.chunk_size,
                )
                print(f"PoseBusters [{protein.upper()} {method.upper()}]: {method_output}")

        if args.prepare_only:
            print(f"Prepare-only mode complete. Protein+method tables created: {created}")
        else:
            print(f"PoseBusters complete. Protein+method result CSVs created: {created}")
        return

    # Legacy mode: one combined table/results CSV
    _write_table(all_rows, table_csv)
    print(f"Table: {table_csv}")
    if args.prepare_only:
        print("Prepare-only mode: skipping PoseBusters run.")
        return
    _run_posebusters(
        table_csv=table_csv,
        output_csv=output_csv,
        bust_bin=args.bust_bin,
        max_workers=args.max_workers,
        posebusters_config=(args.posebusters_config.strip() or None),
        chunk_size=args.chunk_size,
    )
    print(f"PoseBusters results: {output_csv}")


if __name__ == "__main__":
    main()

