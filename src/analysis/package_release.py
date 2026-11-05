#!/usr/bin/env python3
"""Build GitHub example (16×5×5) and six Zenodo packaging trees.

Does not re-download ChEMBL/PDB. Copies local snapshots. Does not copy the
full docking dump into git; record 03 gets a packing list plus a script.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
from pathlib import Path

LEGACY_EXAMPLE_PDBS = ("2z5x", "3mjg", "7awe")
_EXAMPLE_LIGAND_FILE_RE = re.compile(
    r"(?:^|/)(?:ligand_000[1-5]|ligand_[1-5]|idx_[0-4])(?!\d)"
)


def is_legacy_dynamicbind_posebusters(path: Path) -> bool:
    """Superseded 20-check DynamicBind CSVs (renamed *_dynamicbind_old.csv).

    After filter_posebusters_to_10711.py, the active 21-check series is
    ``*_dynamicbind.csv`` (renamed from ``*_dynamicbind_new.csv``).
    """
    name = path.name.lower()
    return name.endswith("_dynamicbind_old.csv")

_SRC = Path(__file__).resolve().parent.parent
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from analysis.ligand_identity import write_maps  # noqa: E402
from analysis.panel import (  # noqa: E402
    CANONICAL_METHODS,
    CANONICAL_TARGETS,
    EXAMPLE_LIGAND_COUNT,
    PDB_TO_RAW_CSV,
    project_root,
)
from docking_benchmark2.ligand_ids import format_ligand_id, vina_legacy_id  # noqa: E402
from docking_benchmark2.preprocessing.curate_chembl import (  # noqa: E402
    curate_from_frozen_nodubl,
    write_csv,
)
from analysis.panel import nodubl_csv_name  # noqa: E402

RECORD_NAMES = (
    "01_raw",
    "02_prepared",
    "03_docking",
    "04_merged",
    "05_posebusters",
    "06_metrics_figures",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_sha256_manifest(root: Path, manifest: Path) -> None:
    rows = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        if path.name in ("SHA256SUMS.txt", "NOTICE"):
            continue
        rel = path.relative_to(root).as_posix()
        rows.append(f"{sha256_file(path)}  {rel}")
    manifest.write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")


def copy_file(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)


def copy_notice(dest_dir: Path, base: Path) -> None:
    src = base / "NOTICE"
    if src.is_file():
        copy_file(src, dest_dir / "NOTICE")


def copy_raw_ligands(base: Path) -> None:
    src_dir = Path(
        "/mnt/tank/scratch/okonovalova/AffiTox/docking-benchmark-2/data/input/ligands"
    )
    dest = base / "input" / "ligands"
    dest.mkdir(parents=True, exist_ok=True)
    for pdb, name in PDB_TO_RAW_CSV.items():
        src = src_dir / name
        if not src.is_file():
            print(f"WARNING: missing raw ligand CSV for {pdb}: {src}")
            continue
        copy_file(src, dest / name)
        print(f"raw ligands {pdb} -> {dest / name}")


def build_curated_and_maps(base: Path) -> int:
    maps = write_maps(base, list(CANONICAL_TARGETS), EXAMPLE_LIGAND_COUNT)
    curated_dir = base / "input" / "ligands_curated"
    curated_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    for pdb in CANONICAL_TARGETS:
        src = base / "input" / "ligands_nodubl" / nodubl_csv_name(pdb)
        curated = curate_from_frozen_nodubl(src)
        write_csv(curated, curated_dir / f"{pdb}_ligands.csv")
        n += len(curated)
        print(f"curated {pdb} {len(curated)}")
    (curated_dir / "N_panel.txt").write_text(f"{n}\n", encoding="utf-8")
    print(f"N_panel={n} (curated / merged rows)")
    return n


def stage_zenodo_01(base: Path, out: Path) -> None:
    rec = out / "01_raw"
    rec.mkdir(parents=True, exist_ok=True)
    proteins = rec / "proteins"
    ligands = rec / "ligands"
    proteins.mkdir(exist_ok=True)
    ligands.mkdir(exist_ok=True)
    for pdb in CANONICAL_TARGETS:
        src = base / "input" / "proteins" / f"{pdb}.pdb"
        if src.is_file():
            copy_file(src, proteins / f"{pdb}.pdb")
        raw_name = PDB_TO_RAW_CSV[pdb]
        src_l = base / "input" / "ligands" / raw_name
        if src_l.is_file():
            copy_file(src_l, ligands / raw_name)
    (rec / "download_manifest.tsv").write_text(
        "pdb_id\tprotein_url\tligand_file\tsource\n"
        + "\n".join(
            f"{pdb}\thttps://files.rcsb.org/download/{pdb.upper()}.pdb\t"
            f"{PDB_TO_RAW_CSV[pdb]}\tChEMBL API snapshot (local copy; not a live query)"
            for pdb in CANONICAL_TARGETS
        )
        + "\n",
        encoding="utf-8",
    )
    (rec / "README.md").write_text(
        _readme_01(),
        encoding="utf-8",
    )
    copy_notice(rec, base)


def stage_zenodo_02(base: Path, out: Path) -> None:
    rec = out / "02_prepared"
    rec.mkdir(parents=True, exist_ok=True)
    for pdb in CANONICAL_TARGETS:
        src = base / "input" / "ligands_curated" / f"{pdb}_ligands.csv"
        if src.is_file():
            copy_file(src, rec / "ligands" / f"{pdb}_ligands.csv")
        src_map = base / "processed" / "id_maps" / f"{pdb}_ligand_id_map.csv"
        if src_map.is_file():
            copy_file(src_map, rec / "id_maps" / f"{pdb}_ligand_id_map.csv")
        box = base / "processed" / "boxes" / f"{pdb}.json"
        if box.is_file():
            copy_file(box, rec / "boxes" / f"{pdb}.json")
        seq = base / "processed" / "plapt_sequences" / f"{pdb}.txt"
        if seq.is_file():
            copy_file(seq, rec / "plapt_sequences" / f"{pdb}.txt")
        for name in (f"{pdb}.pdb", f"{pdb}.pdbqt"):
            prot = base / "processed" / "proteins" / name
            if prot.is_file():
                copy_file(prot, rec / "proteins" / name)
    panel_map = base / "processed" / "id_maps" / "panel_ligand_id_map.csv"
    if panel_map.is_file():
        copy_file(panel_map, rec / "id_maps" / "panel_ligand_id_map.csv")
    n_file = base / "input" / "ligands_curated" / "N_panel.txt"
    if n_file.is_file():
        copy_file(n_file, rec / "N_panel.txt")
    (rec / "weights").mkdir(parents=True, exist_ok=True)
    (rec / "weights" / "README.md").write_text(_readme_weights(), encoding="utf-8")
    (rec / "README.md").write_text(_readme_02(), encoding="utf-8")
    copy_notice(rec, base)


def _docking_inventory(base: Path) -> list[dict]:
    rows = []
    results = base / "results"
    for pdb in CANONICAL_TARGETS:
        docking = results / pdb / "docking"
        for method in ("qvina", "gnina", "plapt", "dynamicbind_new", "dynamicbind"):
            path = docking / method
            exists = path.is_dir()
            rows.append(
                {
                    "pdb_id": pdb,
                    "method": method,
                    "path": str(path),
                    "exists": exists,
                }
            )
        boltz = list(docking.glob("boltz2_*")) if docking.is_dir() else []
        rows.append(
            {
                "pdb_id": pdb,
                "method": "boltz2",
                "path": str(boltz[0]) if boltz else str(docking / "boltz2_*"),
                "exists": bool(boltz),
            }
        )
    return rows


def stage_zenodo_03(base: Path, out: Path) -> None:
    rec = out / "03_docking"
    rec.mkdir(parents=True, exist_ok=True)
    inventory = _docking_inventory(base)
    (rec / "inventory.json").write_text(
        json.dumps(inventory, indent=2), encoding="utf-8"
    )
    pack = rec / "pack_from_workspace.sh"
    pack.write_text(
        "#!/usr/bin/env bash\n"
        "# Zip full docking dumps from the workspace results/ tree.\n"
        "# GPU re-runs of Boltz-2 / DynamicBind are not bit-identical; this zip is the record.\n"
        "set -euo pipefail\n"
        "ROOT=\"${TOXAFFINITY_ROOT:-$(cd \"$(dirname \"$0\")/../../..\" && pwd)}\"\n"
        "DEST=\"${1:-$ROOT/release/zenodo/03_docking/payload}\"\n"
        "mkdir -p \"$DEST\"\n"
        "for pdb in 1g5m 2v5z 3eyg 3jy9 3lxk 11ue 4ase 4f65 4tz4 4zau 5jkv 5mo4 6gqj 6jok 5lf3 7kk3; do\n"
        "  zip -r -0 \"$DEST/results_${pdb}.zip\" \"$ROOT/results/$pdb\" || true\n"
        "done\n",
        encoding="utf-8",
    )
    pack.chmod(0o755)
    (rec / "README.md").write_text(_readme_03(), encoding="utf-8")
    copy_notice(rec, base)


def stage_zenodo_04(base: Path, out: Path) -> None:
    rec = out / "04_merged"
    rec.mkdir(parents=True, exist_ok=True)
    src_dir = base / "analysis" / "excluding_2z5x_3mjg" / "tables"
    n = 0
    coverage = []
    for pdb in CANONICAL_TARGETS:
        src = src_dir / f"merged_ligands_docking_{pdb}.csv"
        if not src.is_file():
            coverage.append({"pdb_id": pdb, "n": 0, "missing": True})
            continue
        copy_file(src, rec / src.name)
        # Count without pandas dependency here if needed — pandas is available.
        import pandas as pd

        df = pd.read_csv(src)
        n += len(df)
        coverage.append(
            {
                "pdb_id": pdb,
                "n": int(len(df)),
                "has_ligand_id": "ligand_id" in df.columns,
                "qvina": int(df["qvina_affinity_bestpose"].notna().sum())
                if "qvina_affinity_bestpose" in df.columns
                else None,
                "gnina": int(df["gnina_cnn_affinity_bestpose"].notna().sum())
                if "gnina_cnn_affinity_bestpose" in df.columns
                else None,
                "plapt": int(df["plapt_affinity"].notna().sum())
                if "plapt_affinity" in df.columns
                else None,
                "dynamicbind": int(df["dynamicbind_affinity_bestpose"].notna().sum())
                if "dynamicbind_affinity_bestpose" in df.columns
                else None,
                "boltz2": int(df["boltz2_affinity_pred_value"].notna().sum())
                if "boltz2_affinity_pred_value" in df.columns
                else None,
            }
        )
    (rec / "merge_audit.json").write_text(
        json.dumps({"N_panel": n, "by_target": coverage}, indent=2), encoding="utf-8"
    )
    (rec / "N_panel.txt").write_text(f"{n}\n", encoding="utf-8")
    (rec / "README.md").write_text(_readme_04(n), encoding="utf-8")
    copy_notice(rec, base)


def stage_zenodo_05(base: Path, out: Path) -> None:
    rec = out / "05_posebusters"
    rec.mkdir(parents=True, exist_ok=True)
    src_dir = base / "analysis" / "excluding_2z5x_3mjg" / "tables" / "posebuster"
    dest = rec / "tables"
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True, exist_ok=True)
    if src_dir.is_dir():
        for csv_path in sorted(src_dir.glob("*.csv")):
            name = csv_path.name.lower()
            if is_legacy_dynamicbind_posebusters(csv_path):
                continue
            if any(old in name for old in LEGACY_EXAMPLE_PDBS):
                continue
            copy_file(csv_path, dest / csv_path.name)
    (rec / "README.md").write_text(_readme_05(), encoding="utf-8")
    copy_notice(rec, base)


def _copy_png_tree(src: Path, dest: Path) -> int:
    n = 0
    if not src.is_dir():
        return 0
    for png in src.rglob("*.png"):
        rel = png.relative_to(src).as_posix().lower()
        if any(old in rel for old in LEGACY_EXAMPLE_PDBS):
            continue
        copy_file(png, dest / png.relative_to(src))
        n += 1
    return n


def stage_zenodo_06(base: Path, out: Path, n_panel: int) -> None:
    rec = out / "06_metrics_figures"
    rec.mkdir(parents=True, exist_ok=True)
    for sub in ("tables", "figures"):
        old = rec / sub
        if old.exists():
            shutil.rmtree(old)
    tables = base / "analysis" / "excluding_2z5x_3mjg" / "tables"
    figures = base / "analysis" / "excluding_2z5x_3mjg" / "figures"
    if tables.is_dir():
        for csv_path in tables.rglob("*.csv"):
            name = csv_path.name.lower()
            parts_l = {p.lower() for p in csv_path.parts}
            # Raw PoseBusters CSVs belong in record 05; keep summary tables under
            # tables/posebusters/. Skip local pre-filter backups.
            if any(
                p == "posebuster" or p.startswith("posebuster_backup") for p in parts_l
            ):
                continue
            if any(old in name for old in LEGACY_EXAMPLE_PDBS):
                continue
            copy_file(csv_path, rec / "tables" / csv_path.relative_to(tables))
    n_png = 0
    for variant in ("by_protein", "by_pdb", ""):
        src = figures / variant if variant else figures
        dest = rec / "figures" / (variant or "root")
        n_png += _copy_png_tree(src, dest)
    (rec / "N_panel.txt").write_text(f"{n_panel}\n", encoding="utf-8")
    (rec / "README.md").write_text(_readme_06(n_panel, n_png), encoding="utf-8")
    copy_notice(rec, base)


def _copy_example_docking(base: Path, dest_root: Path, pdb: str, n: int) -> None:
    docking = base / "results" / pdb / "docking"
    if not docking.is_dir():
        return
    for i in range(1, n + 1):
        legacy = vina_legacy_id(i)
        padded = format_ligand_id(i)
        for method in ("qvina", "gnina"):
            src_dir = docking / method
            if not src_dir.is_dir():
                continue
            for suffix in (".log", ".pdbqt", "_out.pdbqt"):
                src = src_dir / f"{legacy}{suffix}"
                if src.is_file():
                    copy_file(
                        src,
                        dest_root / method / pdb / f"{padded}{suffix}",
                    )
        plapt = docking / "plapt"
        if plapt.is_dir():
            for src in plapt.rglob(f"{legacy}.json"):
                copy_file(
                    src,
                    dest_root / "plapt" / pdb / f"{padded}.json",
                )
        db = docking / "dynamicbind_new"
        if db.is_dir():
            idx = i - 1
            for folder in db.glob(f"**/index{idx}_idx_{idx}"):
                if folder.is_dir():
                    shutil.copytree(
                        folder,
                        dest_root / "dynamicbind" / pdb / padded,
                        dirs_exist_ok=True,
                    )
    # Boltz: copy CHEMBL folders for the first n map rows
    import pandas as pd

    map_path = base / "processed" / "id_maps" / f"{pdb}_ligand_id_map.csv"
    if not map_path.is_file():
        return
    mp = pd.read_csv(map_path).head(n)
    boltz_dirs = list(docking.glob("boltz2_*"))
    if not boltz_dirs:
        return
    for _, row in mp.iterrows():
        chembl = str(row.get("boltz_native_id") or "").strip()
        ligand_id = str(row["ligand_id"])
        if not chembl:
            continue
        hits = list(boltz_dirs[0].rglob(chembl))
        dirs = [p for p in hits if p.is_dir()]
        if not dirs:
            continue
        shutil.copytree(
            dirs[0],
            dest_root / "boltz2" / pdb / ligand_id,
            dirs_exist_ok=True,
        )


def stage_example(
    base: Path, n: int = EXAMPLE_LIGAND_COUNT, copy_docking: bool = True
) -> None:
    import pandas as pd

    ex = base / "example"
    for pdb in CANONICAL_TARGETS:
        raw_name = PDB_TO_RAW_CSV[pdb]
        src = base / "input" / "ligands" / raw_name
        if src.is_file():
            copy_file(src, ex / "01_raw" / "ligands" / raw_name)
        prot = base / "input" / "proteins" / f"{pdb}.pdb"
        if prot.is_file():
            copy_file(prot, ex / "01_raw" / "proteins" / f"{pdb}.pdb")
        curated = base / "input" / "ligands_curated" / f"{pdb}_ligands.csv"
        if curated.is_file():
            df = pd.read_csv(curated).head(n)
            dest = ex / "02_prepared" / "ligands" / f"{pdb}_ligands.csv"
            dest.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(dest, index=False)
        src_map = base / "processed" / "id_maps" / f"{pdb}_example_ligand_id_map.csv"
        if src_map.is_file():
            copy_file(src_map, ex / "02_prepared" / "id_maps" / src_map.name)
        box = base / "processed" / "boxes" / f"{pdb}.json"
        if box.is_file():
            copy_file(box, ex / "02_prepared" / "boxes" / f"{pdb}.json")
        for name in (f"{pdb}.pdb", f"{pdb}.pdbqt"):
            p = base / "processed" / "proteins" / name
            if p.is_file():
                copy_file(p, ex / "02_prepared" / "proteins" / name)
        seq = base / "processed" / "plapt_sequences" / f"{pdb}.txt"
        if seq.is_file():
            copy_file(seq, ex / "02_prepared" / "plapt_sequences" / f"{pdb}.txt")
        for i in range(1, n + 1):
            legacy = vina_legacy_id(i)
            padded = format_ligand_id(i)
            lig_dir = base / "processed" / "ligands" / pdb
            if lig_dir.is_dir():
                matches = list(lig_dir.rglob(f"{legacy}.pdbqt"))
                if matches:
                    copy_file(
                        matches[0],
                        ex / "02_prepared" / "ligand_pdbqt" / pdb / f"{padded}.pdbqt",
                    )
        if copy_docking:
            _copy_example_docking(base, ex / "03_docking", pdb, n)
        merged = (
            base
            / "analysis"
            / "excluding_2z5x_3mjg"
            / "tables"
            / f"merged_ligands_docking_{pdb}.csv"
        )
        if merged.is_file() and "ligand_id" in pd.read_csv(merged, nrows=0).columns:
            dfm = pd.read_csv(merged)
            keep = {format_ligand_id(i) for i in range(1, n + 1)}
            dfm = dfm[dfm["ligand_id"].isin(keep)]
            dest = ex / "04_merged" / merged.name
            dest.parent.mkdir(parents=True, exist_ok=True)
            dfm.to_csv(dest, index=False)
        pb = base / "analysis" / "excluding_2z5x_3mjg" / "tables" / "posebuster"
        if pb.is_dir():
            for csv_path in pb.glob(f"*{pdb}*"):
                if is_legacy_dynamicbind_posebusters(csv_path):
                    continue
                _slice_posebusters_csv(
                    csv_path, ex / "05_posebusters" / "tables" / csv_path.name, n=n
                )
    prune_legacy_example(ex)
    (ex / "README.md").write_text(_readme_example(n), encoding="utf-8")
    copy_notice(ex, base)


def _slice_posebusters_csv(src: Path, dest: Path, n: int) -> None:
    import pandas as pd

    df = pd.read_csv(src)
    if "file" in df.columns:
        sliced = df[df["file"].astype(str).map(lambda s: bool(_EXAMPLE_LIGAND_FILE_RE.search(s)))]
        if sliced.empty:
            sliced = df.head(n)
        df = sliced
    dest.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(dest, index=False)


def prune_legacy_example(ex: Path) -> None:
    """Drop 2z5x/3mjg/7awe leftovers and old DynamicBind PoseBusters CSVs."""
    if not ex.is_dir():
        return
    for path in sorted(ex.rglob("*"), reverse=True):
        name = path.name.lower()
        drop = any(old in name for old in LEGACY_EXAMPLE_PDBS) or (
            path.is_file() and is_legacy_dynamicbind_posebusters(path)
        )
        if not drop:
            continue
        if path.is_dir():
            shutil.rmtree(path, ignore_errors=True)
        elif path.is_file():
            path.unlink(missing_ok=True)


def _readme_01() -> str:
    return """# AffiTox Zenodo record 01 — raw ChEMBL + PDB

Snapshot of the ChEMBL API extracts and RCSB PDB files used for the 16-target panel.
Do not query the live ChEMBL API to reproduce the paper; use these files.

Targets (this order): 1g5m, 2v5z, 3eyg, 3jy9, 3lxk, 11ue, 4ase, 4f65, 4tz4, 4zau, 5jkv, 5mo4, 6gqj, 6jok, 5lf3, 7kk3.

ChEMBL data remain under EMBL-EBI / ChEMBL terms. PDB structures remain under wwPDB terms. See NOTICE.
"""


def _readme_02() -> str:
    return """# AffiTox Zenodo record 02 — prepared inputs

Curated ligand tables with stable `ligand_id` (ligand_0001, …) in docked-table row order,
Meeko receptors, search boxes, PLAPT sequences, and engine-native ID maps.

Join key for every later record is `(pdb_id, ligand_id)`. Historical files named `ligand_1`
or `idx_0` are aliases listed in `id_maps/`.

Place Boltz-2, DynamicBind, and PLAPT weight files under `weights/` (see weights/README.md).
"""


def _readme_weights() -> str:
    return """# Model weights (deposit in this record)

Copy the pinned checkpoints here before upload. Analysis from merged tables does not need them.

| Method | Typical local path | Files |
|--------|--------------------|-------|
| Boltz-2 | `~/.boltz/` | `boltz2_conf.ckpt`, `boltz2_aff.ckpt` |
| DynamicBind | DynamicBind `workdir/big_score_model_sanyueqi_with_time/` | score model checkpoints |
| PLAPT | WELP-PLAPT `models/affinity_predictor.onnx` | ONNX (MD5 `4a340ab3a3179417ce41f3f44328ee38` in SOFTWARE_REGISTRY.md) |

Record SHA256 in this folder after copy. Third-party weights keep their original licenses (NOTICE).
"""


def _readme_03() -> str:
    return """# AffiTox Zenodo record 03 — docking outputs

Frozen poses, logs, and native engine JSON/SDF for the 16-target panel.

Re-running Boltz-2 or DynamicBind on a GPU is **not bit-identical**. This record is the
scientific snapshot. Layers 4–6 must be rebuilt from these files, not from a new GPU run,
if you need the same numbers.

Run `pack_from_workspace.sh` on the cluster to zip `results/<pdb>/`.
"""


def _readme_04(n: int) -> str:
    return f"""# AffiTox Zenodo record 04 — merged tables

One CSV per target: `merged_ligands_docking_<pdb>.csv`.

Join key: `(protein_pdb_id, ligand_id)` only. Missing method scores stay empty (never imputed).
Panel N (curated / merged rows): **{n}**.

Seed: 42 (documented for docking; QVina2 logs may print a different engine seed — see SOFTWARE_REGISTRY.md).
"""


def _readme_05() -> str:
    return """# AffiTox Zenodo record 05 — PoseBusters and pose conversion

Full-panel PoseBusters CSVs belong here, not on GitHub.

GitHub `example/05_posebusters/` holds only the 16 × 5 ligand slice.
PLAPT has no poses and is excluded.

DynamicBind tables are `posebusters_results_<pdb>_dynamicbind.csv` (21 boolean
checks, including `internal_energy`; renamed from ``*_dynamicbind_new.csv``).
Superseded 20-check ``*_dynamicbind_old.csv`` files are not included. 2z5x, 3mjg,
and 7awe are not in this panel.
"""


def _readme_06(n: int, n_png: int) -> str:
    return f"""# AffiTox Zenodo record 06 — metrics and figures

Correlation, nEF, inferential tables, and PNG figures for the canonical panel
(`analysis/excluding_2z5x_3mjg`, `figures/by_protein` as the publication axis).

Heatmap style: ultramarine / ivory / vermilion, square cells
(`pearson_heatmap_ultramarine_ivory_vermilion_square_cells.png` and the same style
for Spearman-only, Kendall-only, combined, and intermediate heatmaps). PNG only.

PoseBusters pass-count curves use 21 boolean checks (including
``internal_energy``); the at-least-21 endpoint matches ``pass_rate_all``.

Panel N: **{n}**. PNG files packaged in this pass: {n_png}.
"""


def _readme_example(n: int) -> str:
    return f"""# AffiTox GitHub example (16 targets × {n} ligands × 5 methods)

Same six information blocks as the Zenodo records, truncated to `ligand_0001` … `ligand_000{n:04d}`
for every target and every method. Seed 42. If a method failed one of these five ligands, the
row stays with a missing score; do not substitute `ligand_0006`.

```bash
source config/project.env.sh
# EXAMPLE_N={n}  # first five ligand_id rows; omit the variable for the full panel
bash pipeline/prepare/run_curate.sh
PYTHONPATH=src python src/analysis/ligand_identity.py
bash pipeline/postprocess/run_merge.sh   # full panel → analysis/excluding_2z5x_3mjg/tables
EXAMPLE_N={n} bash pipeline/postprocess/run_merge.sh --output-dir example/04_merged
```

Full-panel PoseBusters conversion and docking dumps are on Zenodo records 03 and 05.
"""


def reslice_example_posebusters(base: Path, n: int = EXAMPLE_LIGAND_COUNT) -> None:
    src_dir = base / "analysis" / "excluding_2z5x_3mjg" / "tables" / "posebuster"
    dest_dir = base / "example" / "05_posebusters" / "tables"
    dest_dir.mkdir(parents=True, exist_ok=True)
    if not src_dir.is_dir():
        return
    for pdb in CANONICAL_TARGETS:
        for csv_path in src_dir.glob(f"*{pdb}*"):
            if is_legacy_dynamicbind_posebusters(csv_path):
                continue
            _slice_posebusters_csv(csv_path, dest_dir / csv_path.name, n=n)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", type=Path, default=None)
    parser.add_argument("--skip-example-docking", action="store_true")
    parser.add_argument("--prune-example-only", action="store_true")
    args = parser.parse_args()
    base = args.base_dir or project_root()
    if args.prune_example_only:
        reslice_example_posebusters(base)
        prune_legacy_example(base / "example")
        print(f"Pruned example at {base / 'example'}")
        return 0
    zenodo = base / "release" / "zenodo"
    zenodo.mkdir(parents=True, exist_ok=True)

    copy_raw_ligands(base)
    n_panel = build_curated_and_maps(base)
    stage_zenodo_01(base, zenodo)
    stage_zenodo_02(base, zenodo)
    stage_zenodo_03(base, zenodo)
    stage_zenodo_04(base, zenodo)
    stage_zenodo_05(base, zenodo)
    stage_zenodo_06(base, zenodo, n_panel)
    stage_example(base, copy_docking=not args.skip_example_docking)
    for rec in RECORD_NAMES:
        write_sha256_manifest(zenodo / rec, zenodo / rec / "SHA256SUMS.txt")
    (zenodo / "README.md").write_text(
        "# Six related Zenodo records (information blocks, not per-protein DOIs)\n\n"
        "Upload each of `01_raw` … `06_metrics_figures` as its own record, then link them "
        "from one Zenodo community/project. Panel: 16 PDB codes listed in 01_raw/README.md.\n"
        f"Frozen ligand_id N (curated rows): {n_panel}.\n",
        encoding="utf-8",
    )
    print(f"Done. N_panel={n_panel}. Trees: {zenodo} and {base / 'example'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
