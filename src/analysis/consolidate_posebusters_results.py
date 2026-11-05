#!/usr/bin/env python3
"""Consolidate PoseBusters CSVs and audit them against available pose files."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
from pathlib import Path


TARGETS = (
    "1g5m",
    "2v5z",
    "3eyg",
    "3jy9",
    "3lxk",
    "11ue",
    "4ase",
    "4f65",
    "4tz4",
    "4zau",
    "5jkv",
    "5mo4",
    "6gqj",
    "6jok",
    "5lf3",
    "7kk3",
)
METHOD_TAGS = ("gnina", "qvina", "dynamicbind", "boltz2")
REFRESH_TARGETS = {"2v5z", "4tz4", "11ue", "7kk3"}


def _csv_rows(path: Path) -> int:
    if not path.is_file() or path.stat().st_size == 0:
        return 0
    with path.open(newline="", encoding="utf-8") as handle:
        return sum(1 for _ in csv.DictReader(handle))


def _available_pose_count(results_dir: Path, target: str, method: str) -> int:
    docking = results_dir / target / "docking"
    if method in {"gnina", "qvina"}:
        return sum(
            1
            for path in (docking / method).glob("ligand_*_out.pdbqt")
            if path.stat().st_size > 0
        )
    if method in ("dynamicbind", "dynamicbind_new"):
        ligand_ids: set[str] = set()
        sub = "dynamicbind" if (docking / "dynamicbind").exists() else "dynamicbind_new"
        for index_dir in (docking / sub).glob("**/index*_idx_*"):
            match = re.fullmatch(r"index(\d+)_idx_\d+", index_dir.name)
            if not match:
                continue
            if list(index_dir.glob("rank1_ligand_*.sdf")) and list(
                index_dir.glob("rank1_receptor_*.pdb")
            ):
                ligand_ids.add(match.group(1))
        return len(ligand_ids)
    if method == "boltz2":
        ligand_ids = {
            path.name.removesuffix("_model_0.cif")
            for path in docking.glob("boltz2*/**/*_model_0.cif")
            if path.stat().st_size > 0
        }
        return len(ligand_ids)
    raise ValueError(method)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--copy-existing", action="store_true")
    args = parser.parse_args()

    source_dir = args.root / "analysis" / "tables" / "posebuster"
    output_dir = args.analysis_root / "tables" / "posebuster"
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.copy_existing:
        for target in TARGETS:
            if target in REFRESH_TARGETS:
                continue
            for method in METHOD_TAGS:
                source = source_dir / f"posebusters_results_{target}_{method}.csv"
                destination = output_dir / source.name
                if source.is_file() and not destination.is_file():
                    shutil.copy2(source, destination)

    records: list[dict[str, str | int]] = []
    for target in TARGETS:
        for method in METHOD_TAGS:
            source_poses = _available_pose_count(args.root / "results", target, method)
            input_path = (
                output_dir
                / "tables_by_method"
                / f"posebusters_input_{target}_{method}.csv"
            )
            evaluable_poses = _csv_rows(input_path)
            result_path = output_dir / f"posebusters_results_{target}_{method}.csv"
            rows = _csv_rows(result_path)
            records.append(
                {
                    "target": target,
                    "method": method,
                    "source_pose_rows": source_poses,
                    "evaluable_pose_rows": evaluable_poses,
                    "excluded_unconvertible_rows": max(
                        source_poses - evaluable_poses, 0
                    ),
                    "posebusters_rows": rows,
                    "missing_posebusters_rows": max(evaluable_poses - rows, 0),
                    "extra_posebusters_rows": max(rows - evaluable_poses, 0),
                    "status": (
                        "complete"
                        if rows == evaluable_poses and evaluable_poses > 0
                        else "incomplete"
                    ),
                    "result_csv": str(result_path),
                }
            )

    report = args.analysis_root / "tables" / "posebusters_input_output_audit.csv"
    report.parent.mkdir(parents=True, exist_ok=True)
    with report.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    incomplete = [row for row in records if row["status"] != "complete"]
    print(f"PoseBusters inventory: {len(records) - len(incomplete)}/{len(records)} complete")
    for row in incomplete:
        print(
            f"INCOMPLETE {row['target']}/{row['method']}: "
            f"evaluable={row['evaluable_pose_rows']} "
            f"results={row['posebusters_rows']}"
        )
    print(report)


if __name__ == "__main__":
    main()
