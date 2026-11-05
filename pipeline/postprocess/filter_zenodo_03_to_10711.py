#!/usr/bin/env python3
"""Filter AffiTox Zenodo record-03 docking zips to the published N=10711 panel.

Reads existing ``release/zenodo/payloads/03/results_<pdb>.zip`` (frozen dump) and
writes panel-filtered copies under ``payloads/03_panel10711/``. Keep sets come from
``processed/id_maps/<pdb>_ligand_id_map.csv`` (already filtered to curated 10711).

Per-ligand files are dropped when their vina / DynamicBind / CHEMBL token is not
in the keep set. Shared non-ligand paths are kept. Aggregate affinity CSVs that
list ``molecule_chembl_id`` / ``CHEMBL`` columns are row-filtered when possible.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import time
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
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
BUCKETS = [
    ("03a", ["3eyg", "4f65", "5lf3", "2v5z"]),
    ("03b", ["3jy9", "3lxk"]),
    ("03c", ["5mo4", "6gqj", "4zau"]),
    ("03d", ["7kk3", "5jkv", "4tz4"]),
    ("03e", ["4ase", "6jok", "11ue", "1g5m"]),
]

# vina / plapt style: ligand_12.log, ligand_12_out.pdbqt (not ligand_0012 alone ambiguity —
# zero-padded ligand_00012 still yields int 12).
LIGAND_NUM_RE = re.compile(r"(?:^|/)ligand_(\d+)(?:_|\.|/|$)")
DB_INDEX_RE = re.compile(r"(?:^|/)index(\d+)_idx_(\d+)(?:/|$)")
DB_IDX_RE = re.compile(r"(?:^|/)idx_(\d+)(?:_|\.|/|$)")
CHEMBL_RE = re.compile(r"(CHEMBL\d+)", re.IGNORECASE)


@dataclass
class KeepSet:
    pdb: str
    vina_nums: set[int] = field(default_factory=set)
    db_idxs: set[int] = field(default_factory=set)
    chembls: set[str] = field(default_factory=set)
    n_map: int = 0


def load_keep(maps_dir: Path, pdb: str) -> KeepSet:
    path = maps_dir / f"{pdb}_ligand_id_map.csv"
    if not path.is_file():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    keep = KeepSet(pdb=pdb, n_map=len(df))
    for _, row in df.iterrows():
        vina = str(row.get("vina_legacy_id") or "")
        m = re.search(r"(\d+)", vina)
        if m:
            keep.vina_nums.add(int(m.group(1)))
        db = str(row.get("dynamicbind_legacy_id") or "")
        m = re.search(r"(\d+)", db)
        if m:
            keep.db_idxs.add(int(m.group(1)))
        for col in ("boltz_native_id", "molecule_chembl_id"):
            val = str(row.get(col) or "").strip().upper()
            if val.startswith("CHEMBL"):
                keep.chembls.add(val)
        lid = str(row.get("ligand_id") or "")
        m = re.search(r"(\d+)", lid)
        if m:
            keep.vina_nums.add(int(m.group(1)))
    return keep


def member_decision(name: str, keep: KeepSet) -> str:
    """Return 'keep' or 'drop' for a zip member path."""
    # DynamicBind per-pose folders first (most specific).
    db = DB_INDEX_RE.search(name)
    if db:
        idx = int(db.group(2))
        return "keep" if idx in keep.db_idxs else "drop"
    db2 = DB_IDX_RE.search(name)
    if db2 and "/index" not in name:
        # bare idx_N tokens outside index*_idx_* already handled
        idx = int(db2.group(1))
        # only treat as DynamicBind id when under a dynamicbind path
        if "dynamicbind" in name.lower():
            return "keep" if idx in keep.db_idxs else "drop"

    chembls = {c.upper() for c in CHEMBL_RE.findall(name)}
    if chembls:
        if chembls & keep.chembls:
            return "keep"
        return "drop"

    ligs = [int(x) for x in LIGAND_NUM_RE.findall(name)]
    if ligs:
        # keep only if every captured ligand number is in the panel
        if all(n in keep.vina_nums for n in ligs):
            return "keep"
        return "drop"

    return "keep"


def _filter_csv_bytes(raw: bytes, keep: KeepSet) -> bytes | None:
    """Row-filter affinity-like CSVs; return None to keep original bytes."""
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        return None
    try:
        reader = csv.DictReader(io.StringIO(text))
        if not reader.fieldnames:
            return None
        fields = [f.lower() for f in reader.fieldnames]
        chembl_col = None
        for cand in (
            "molecule_chembl_id",
            "chembl_id",
            "boltz_native_id",
            "ligand_id",
            "name",
            "id",
        ):
            if cand in fields:
                chembl_col = reader.fieldnames[fields.index(cand)]
                break
        if chembl_col is None:
            return None
        out = io.StringIO()
        writer = csv.DictWriter(out, fieldnames=reader.fieldnames, lineterminator="\n")
        writer.writeheader()
        kept = 0
        for row in reader:
            val = str(row.get(chembl_col) or "").strip()
            up = val.upper()
            ok = False
            if up.startswith("CHEMBL") and up in keep.chembls:
                ok = True
            else:
                m = re.search(r"(\d+)", val)
                if m and int(m.group(1)) in keep.vina_nums:
                    ok = True
                if m and int(m.group(1)) in keep.db_idxs and "idx" in val.lower():
                    ok = True
            if ok:
                writer.writerow(row)
                kept += 1
        if kept == 0:
            return None  # safer: leave original if we matched nothing
        return out.getvalue().encode("utf-8")
    except Exception:
        return None


def filter_zip(
    src: Path,
    dest: Path,
    keep: KeepSet,
    *,
    dry_run: bool = False,
) -> dict:
    t0 = time.time()
    with zipfile.ZipFile(src, "r") as zin:
        infos = zin.infolist()
        n_in = len(infos)
        bytes_in = sum(i.file_size for i in infos)
        keep_infos: list[zipfile.ZipInfo] = []
        drop_infos: list[zipfile.ZipInfo] = []
        for info in infos:
            if info.is_dir():
                keep_infos.append(info)
                continue
            if member_decision(info.filename, keep) == "keep":
                keep_infos.append(info)
            else:
                drop_infos.append(info)
        bytes_keep = sum(i.file_size for i in keep_infos if not i.is_dir())
        bytes_drop = sum(i.file_size for i in drop_infos)
        report = {
            "pdb": keep.pdb,
            "n_map": keep.n_map,
            "members_in": n_in,
            "members_keep": len(keep_infos),
            "members_drop": len(drop_infos),
            "bytes_in": bytes_in,
            "bytes_keep": bytes_keep,
            "bytes_drop": bytes_drop,
            "src": str(src),
            "dest": str(dest),
            "dry_run": dry_run,
        }
        if dry_run:
            report["seconds"] = round(time.time() - t0, 1)
            return report

        dest.parent.mkdir(parents=True, exist_ok=True)
        tmp = dest.with_suffix(dest.suffix + ".partial")
        if tmp.exists():
            tmp.unlink()
        with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_STORED, allowZip64=True) as zout:
            for info in keep_infos:
                if info.is_dir():
                    continue
                data = zin.read(info)
                name_l = info.filename.lower()
                if name_l.endswith(".csv") and (
                    "affinity" in name_l or "prediction" in name_l or "score" in name_l
                ):
                    filtered = _filter_csv_bytes(data, keep)
                    if filtered is not None:
                        data = filtered
                # rebuild ZipInfo to avoid copy of external attrs quirks
                out_info = zipfile.ZipInfo(filename=info.filename, date_time=info.date_time)
                out_info.compress_type = zipfile.ZIP_STORED
                out_info.external_attr = info.external_attr
                zout.writestr(out_info, data)
        tmp.replace(dest)
        report["bytes_out"] = dest.stat().st_size
        report["seconds"] = round(time.time() - t0, 1)
        return report


def write_buckets(zip_dir: Path, *, rel_prefix: str = "03") -> Path:
    manifest = {"records": [], "panel_n": 10711}
    for name, pdbs in BUCKETS:
        files = []
        total = 0
        for pdb in pdbs:
            z = zip_dir / f"results_{pdb}.zip"
            size = z.stat().st_size if z.is_file() else 0
            files.append(
                {
                    "pdb": pdb,
                    "zip": f"{rel_prefix}/results_{pdb}.zip",
                    "bytes": size,
                }
            )
            total += size
        manifest["records"].append({"id": name, "files": files, "bytes": total})
    path = zip_dir.parent / "03_buckets.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    parser.add_argument(
        "--src-dir",
        type=Path,
        default=None,
        help="Directory with results_<pdb>.zip (default: release/zenodo/payloads/03)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory (default: release/zenodo/payloads/03_panel10711)",
    )
    parser.add_argument("--pdb", action="append", default=None, help="Limit to PDB(s)")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--write-buckets", action="store_true")
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="Write beside source then replace it (saves a full second copy of 03).",
    )
    parser.add_argument(
        "--skip-existing-marker",
        action="store_true",
        help="Skip zips that already have a sibling .panel10711_done marker.",
    )
    args = parser.parse_args()
    root = args.root
    src_dir = args.src_dir or (root / "release" / "zenodo" / "payloads" / "03")
    out_dir = args.out_dir or (root / "release" / "zenodo" / "payloads" / "03_panel10711")
    maps_dir = root / "processed" / "id_maps"
    pdbs = args.pdb or list(TARGETS)
    reports = []
    out_dir.mkdir(parents=True, exist_ok=True)
    for pdb in pdbs:
        src = src_dir / f"results_{pdb}.zip"
        if not src.is_file():
            print(f"WARNING: missing {src}", file=sys.stderr, flush=True)
            continue
        marker = src.with_suffix(src.suffix + ".panel10711_done")
        if args.skip_existing_marker and marker.is_file():
            print(f"[{pdb}] skip (marker {marker.name})", flush=True)
            continue
        keep = load_keep(maps_dir, pdb)
        if args.in_place:
            dest = src.with_suffix(src.suffix + ".filtered")
        else:
            dest = out_dir / f"results_{pdb}.zip"
        print(
            f"[{pdb}] map={keep.n_map} filter {src.name} "
            f"({src.stat().st_size / 1024**3:.2f} GiB) -> {dest}",
            flush=True,
        )
        rep = filter_zip(src, dest, keep, dry_run=args.dry_run)
        if args.in_place and not args.dry_run:
            src.unlink()
            dest.replace(src)
            marker.write_text(
                json.dumps(
                    {
                        "pdb": pdb,
                        "n_map": keep.n_map,
                        "members_keep": rep["members_keep"],
                        "members_drop": rep["members_drop"],
                        "bytes_drop": rep["bytes_drop"],
                    },
                    indent=2,
                )
                + "\n"
            )
            rep["dest"] = str(src)
            rep["in_place"] = True
        reports.append(rep)
        print(
            f"  keep {rep['members_keep']}/{rep['members_in']} members; "
            f"drop {rep['bytes_drop']/1024**3:.2f} GiB of {rep['bytes_in']/1024**3:.2f} GiB "
            f"in {rep['seconds']}s",
            flush=True,
        )
    report_path = out_dir / "filter_report.json"
    # merge with prior reports if present
    prior = []
    if report_path.is_file():
        try:
            prior = json.loads(report_path.read_text())
        except json.JSONDecodeError:
            prior = []
    by_pdb = {r.get("pdb"): r for r in prior if isinstance(r, dict)}
    for r in reports:
        by_pdb[r["pdb"]] = r
    report_path.write_text(json.dumps(list(by_pdb.values()), indent=2) + "\n")
    print(f"wrote {report_path}", flush=True)
    if args.write_buckets and not args.dry_run:
        buckets = write_buckets(src_dir, rel_prefix="03")
        print(f"wrote {buckets}", flush=True)
        for rec in json.loads(buckets.read_text())["records"]:
            print(f"  {rec['id']}: {rec['bytes']/1024**3:.2f} GiB", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
