#!/usr/bin/env python3
"""Apply pre-dock SMILES collapse to curated ligand tables and prune docking files.

1. Backup assay-level curated CSVs.
2. Collapse each curated table: mean Ki, keep first-occurring ligand_id.
3. Move docking artifacts for non-representative ligand_ids into an archive
   (qvina/gnina/plapt files, DynamicBind idx dirs; Boltz CHEMBL dirs only if
   that CHEMBL is not retained by any kept row).

Does not delete — only shutil.move into archive/.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

_SRC = Path(__file__).resolve().parent.parent
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from analysis.ligand_collapse import duplicate_collapse_plan, group_by_smiles, parse_ligand_number
from analysis.panel import CANONICAL_TARGETS, project_root

LIGAND_STEM_RE = re.compile(r"^ligand_(\d+)(?:$|[._\-])")
IDX_DIR_RE = re.compile(r"^index(\d+)_idx_\d+$")
CHEMBL_RE = re.compile(r"CHEMBL\d+", re.IGNORECASE)


def _read_table(path: Path) -> pd.DataFrame:
    sample = path.read_text(encoding="utf-8", errors="replace")[:2048]
    sep = ";" if sample.count(";") >= sample.count(",") else ","
    return pd.read_csv(path, sep=sep, low_memory=False)


def _move(src: Path, dst: Path, moves: list[dict]) -> None:
    if not src.exists():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        return
    shutil.move(str(src), str(dst))
    moves.append({"src": str(src), "dst": str(dst)})


def _ligand_number_from_name(name: str) -> int | None:
    m = LIGAND_STEM_RE.match(name)
    return int(m.group(1)) if m else None


def prune_docking_tree(
    docking_root: Path,
    archive_root: Path,
    drop_numbers: set[int],
    drop_chembl: set[str],
    keep_chembl: set[str],
    moves: list[dict],
) -> None:
    if not docking_root.is_dir() or not drop_numbers:
        return

    to_move: list[Path] = []

    # qvina / gnina / plapt: ligand_N*
    for method in ("qvina", "gnina", "plapt"):
        method_dir = docking_root / method
        if not method_dir.is_dir():
            continue
        for path in method_dir.rglob("*"):
            n = _ligand_number_from_name(path.name)
            if n is not None and n in drop_numbers:
                to_move.append(path)

    # DynamicBind: **/index{k}_idx_{k}  (0-based → ligand k+1)
    for method in ("dynamicbind", "dynamicbind_old", "dynamicbind_new"):
        method_dir = docking_root / method
        if not method_dir.is_dir():
            continue
        for path in method_dir.rglob("index*_idx_*"):
            if not path.is_dir():
                continue
            m = IDX_DIR_RE.match(path.name)
            if not m:
                continue
            n = int(m.group(1)) + 1
            if n in drop_numbers:
                to_move.append(path)

    # Boltz: only under boltz* dirs; archive CHEMBL exclusive to dropped rows
    for boltz_root in docking_root.glob("boltz*"):
        if not boltz_root.is_dir():
            continue
        for path in boltz_root.rglob("*"):
            if not path.is_dir():
                continue
            name = path.name
            if not (
                name.startswith("boltz_results_CHEMBL")
                or (CHEMBL_RE.fullmatch(name) and path.parent.name == "predictions")
            ):
                continue
            m = CHEMBL_RE.search(name)
            if not m:
                continue
            chembl = m.group(0).upper()
            if chembl in keep_chembl or chembl not in drop_chembl:
                continue
            to_move.append(path)

    # Deepest paths first so parents move after children if both listed
    for path in sorted(to_move, key=lambda p: len(p.parts), reverse=True):
        if not path.exists():
            continue
        try:
            rel = path.relative_to(docking_root)
        except ValueError:
            continue
        _move(path, archive_root / "docking" / rel, moves)


def prune_processed_ligands(
    processed_ligands: Path,
    archive_root: Path,
    drop_numbers: set[int],
    moves: list[dict],
) -> None:
    if not processed_ligands.is_dir() or not drop_numbers:
        return
    to_move: list[Path] = []
    for path in processed_ligands.rglob("*"):
        n = _ligand_number_from_name(path.name)
        if n is not None and n in drop_numbers:
            to_move.append(path)
    for path in sorted(to_move, key=lambda p: len(p.parts), reverse=True):
        if not path.exists():
            continue
        try:
            rel = path.relative_to(processed_ligands)
        except ValueError:
            continue
        _move(path, archive_root / "processed_ligands" / rel, moves)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base-dir", type=Path, default=None)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--skip-docking-prune", action="store_true")
    ap.add_argument("--skip-processed-prune", action="store_true")
    args = ap.parse_args()

    base = args.base_dir or project_root()
    curated_dir = base / "input" / "ligands_curated"
    results_dir = base / "results"
    processed_dir = base / "processed" / "ligands"
    processed_sdf = base / "processed" / "ligands_sdf"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_curated = base / "input" / f"ligands_curated_assay_level_backup_{stamp}"
    archive_root = (
        base
        / "analysis"
        / "excluding_2z5x_3mjg"
        / f"docking_archive_nonrepresentative_{stamp}"
    )
    plan_dir = (
        base
        / "analysis"
        / "excluding_2z5x_3mjg"
        / "tables"
        / "pre_dock_collapse_plan"
    )
    plan_dir.mkdir(parents=True, exist_ok=True)

    if not args.dry_run:
        backup_curated.mkdir(parents=True, exist_ok=True)
        archive_root.mkdir(parents=True, exist_ok=True)

    summary = []
    all_moves: list[dict] = []

    for pdb in CANONICAL_TARGETS:
        pdb = pdb.lower()
        src = curated_dir / f"{pdb}_ligands.csv"
        if not src.exists():
            print(f"SKIP missing {src}", flush=True)
            continue

        raw = _read_table(src)
        plan = duplicate_collapse_plan(raw)
        if len(plan):
            plan = plan.copy()
            plan["pdb"] = pdb.upper()
            plan.to_csv(plan_dir / f"drop_plan_{pdb}.csv", index=False)
        else:
            pd.DataFrame(
                columns=[
                    "canonical_smiles",
                    "molecule_chembl_id_keep",
                    "molecule_chembl_id_drop",
                    "keep_ligand_id",
                    "drop_ligand_id",
                    "keep_ligand_number",
                    "drop_ligand_number",
                    "pdb",
                ]
            ).to_csv(plan_dir / f"drop_plan_{pdb}.csv", index=False)

        n_before = len(raw)
        collapsed = group_by_smiles(raw)
        n_after = len(collapsed)

        drop_numbers: set[int] = set()
        if len(plan) and "drop_ligand_number" in plan.columns:
            drop_numbers = {
                int(x)
                for x in plan["drop_ligand_number"].dropna().astype(int).tolist()
            }

        keep_chembl = {
            str(x).upper()
            for x in collapsed.get("molecule_chembl_id", pd.Series(dtype=str)).dropna()
            if str(x).strip()
        }
        drop_chembl: set[str] = set()
        if len(plan) and "molecule_chembl_id_drop" in plan.columns:
            drop_chembl = {
                str(x).upper()
                for x in plan["molecule_chembl_id_drop"].dropna()
                if str(x).strip() and str(x).upper() not in keep_chembl
            }

        print(
            f"{pdb}: curated {n_before} → {n_after}; "
            f"drop_ligand_files={len(drop_numbers)}; "
            f"exclusive_drop_chembl={len(drop_chembl)}",
            flush=True,
        )

        if not args.dry_run:
            shutil.copy2(src, backup_curated / src.name)
            collapsed.to_csv(src, index=False)

            if not args.skip_docking_prune and drop_numbers:
                prune_docking_tree(
                    results_dir / pdb / "docking",
                    archive_root / pdb,
                    drop_numbers,
                    drop_chembl,
                    keep_chembl,
                    all_moves,
                )
            if not args.skip_processed_prune and drop_numbers:
                prune_processed_ligands(
                    processed_dir / pdb, archive_root / pdb, drop_numbers, all_moves
                )
                prune_processed_ligands(
                    processed_sdf / pdb, archive_root / pdb, drop_numbers, all_moves
                )

        summary.append(
            {
                "pdb": pdb.upper(),
                "n_before": n_before,
                "n_after": n_after,
                "n_drop_ids": len(drop_numbers),
                "n_exclusive_drop_chembl": len(drop_chembl),
            }
        )

    pd.DataFrame(summary).to_csv(plan_dir / "collapse_summary.csv", index=False)
    if not args.dry_run:
        (archive_root / "moves.jsonl").write_text(
            "\n".join(json.dumps(m) for m in all_moves) + ("\n" if all_moves else ""),
            encoding="utf-8",
        )
        (archive_root / "README.md").write_text(
            f"""# Non-representative docking archive ({stamp})

Moved docking/processed files for ligand_ids that are **not** the first-occurring
row in a duplicate-SMILES group after pre-dock mean-Ki collapse.

Curated assay-level backup: `{backup_curated}`

Kept ligand tables: `input/ligands_curated/<pdb>_ligands.csv` (collapsed).
Drop plans: `{plan_dir}`
""",
            encoding="utf-8",
        )
    print("SUMMARY", flush=True)
    print(pd.DataFrame(summary).to_string(index=False), flush=True)
    print(f"moves recorded: {len(all_moves)}", flush=True)
    print(f"plan_dir={plan_dir}", flush=True)
    if not args.dry_run:
        print(f"curated_backup={backup_curated}", flush=True)
        print(f"archive={archive_root}", flush=True)


if __name__ == "__main__":
    main()
