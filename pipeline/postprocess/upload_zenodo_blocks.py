#!/usr/bin/env python3
"""Upload AffiTox information-block records (01–06). Record 03 is sharded.

Does not upload the historical per-protein zips 10.5281/zenodo.20825057–067.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# Reuse HTTP helper from the legacy uploader.
_POST = Path(__file__).resolve().parent
if str(_POST) not in sys.path:
    sys.path.insert(0, str(_POST))
from upload_zenodo import Client, load_token, save_manifest, load_manifest  # noqa: E402

API = os.environ.get("ZENODO_API", "https://zenodo.org/api")
COMMUNITY = os.environ.get("ZENODO_COMMUNITY", "affitox")
GITHUB = "https://github.com/ai-chem/AffiTox"
MAX_RECORD_GB = 50.0

BLOCK_META = [
    {
        "id": "01_raw",
        "title": "AffiTox: raw ChEMBL extracts and PDB structures (record 01)",
        "zip": "01_raw.zip",
        "description": (
            "Snapshot of the ChEMBL API extracts and RCSB PDB files for the canonical "
            "16-target AffiTox panel (1g5m, 2v5z, 3eyg, 3jy9, 3lxk, 11ue, 4ase, 4f65, "
            "4tz4, 4zau, 5jkv, 5mo4, 6gqj, 6jok, 5lf3, 7kk3). Do not query live ChEMBL "
            f"to reproduce the paper. Code: {GITHUB}."
        ),
    },
    {
        "id": "02_prepared",
        "title": "AffiTox: prepared inputs, ligand_id maps, and method weights (record 02)",
        "zip": "02_prepared.zip",
        "description": (
            "Curated ligand tables with stable ligand_id, Meeko receptors, boxes, "
            "PLAPT sequences, ID maps, and pinned Boltz-2 / DynamicBind / PLAPT weights. "
            f"Join key is (pdb_id, ligand_id). Code: {GITHUB}."
        ),
    },
    {
        "id": "04_merged",
        "title": "AffiTox: merged ligand_id tables (record 04)",
        "zip": "04_merged.zip",
        "description": (
            "One merged_ligands_docking_<pdb>.csv per target after joining on ligand_id. "
            "Panel N = 10711 curated / merged rows. Rebuild from record 03 "
            f"to match these numbers. Code: {GITHUB}."
        ),
    },
    {
        "id": "05_posebusters",
        "title": "AffiTox: PoseBusters tables (record 05)",
        "zip": "05_posebusters.zip",
        "description": (
            "Full-panel PoseBusters CSVs for Boltz-2, DynamicBind, GNINA, and QVina2. "
            f"PLAPT has no poses and is excluded. Code: {GITHUB}."
        ),
    },
    {
        "id": "06_metrics_figures",
        "title": "AffiTox: metrics tables and PNG figures (record 06)",
        "zip": "06_metrics_figures.zip",
        "description": (
            "Correlation, nEF, inferential tables and publication PNGs for "
            "analysis/excluding_2z5x_3mjg (by_protein and by_pdb). "
            f"Code: {GITHUB}."
        ),
    },
]


def metadata_for(title: str, description: str) -> dict:
    return {
        "title": title,
        "upload_type": "dataset",
        "description": description,
        "creators": [
            {"name": "Konovalova, Olga A.", "orcid": "0000-0003-0391-3211"}
        ],
        "license": "cc-by-4.0",
        "keywords": ["docking", "toxicity", "benchmark", "ChEMBL", "AffiTox"],
        "communities": [{"identifier": COMMUNITY}],
        "related_identifiers": [
            {
                "identifier": GITHUB,
                "relation": "isDocumentedBy",
                "resource_type": "software",
            }
        ],
        "notes": (
            "Information-block data record for AffiTox v1.0.0. "
            "Historical per-protein zips 10.5281/zenodo.20825057–067 are not this release."
        ),
    }


def create_deposition(client: Client, metadata: dict) -> dict:
    body = json.dumps({"metadata": metadata}).encode()
    _, raw = client.request(
        "POST",
        f"{API}/deposit/depositions",
        data=body,
        headers={"Content-Type": "application/json"},
    )
    return json.loads(raw)


def _curl_upload(url: str, path: Path, header_file: str) -> int:
    """PUT a file with retries. Token lives in header_file, not argv."""
    cmd = [
        "curl",
        "-fS",
        "--retry",
        "40",
        "--retry-delay",
        "30",
        "--retry-all-errors",
        "--retry-max-time",
        "86400",
        "--connect-timeout",
        "60",
        "--max-time",
        "0",
        "--upload-file",
        str(path),
        "-H",
        f"@{header_file}",
        url,
    ]
    proc = subprocess.run(cmd, check=False)
    return proc.returncode


def upload_file(client: Client, bucket_url: str, path: Path) -> None:
    url = f"{bucket_url}/{path.name}"
    size = path.stat().st_size
    print(f"    uploading {path.name} ({size / 1024**3:.2f} GiB) ...", flush=True)
    fd, header_file = tempfile.mkstemp(prefix="zenodo_hdr_", text=True)
    try:
        os.chmod(header_file, 0o600)
        with os.fdopen(fd, "w") as fh:
            fh.write(f"Authorization: Bearer {client.token}\n")
        last_rc = 1
        t0 = time.time()
        for attempt in range(1, 9):
            last_rc = _curl_upload(url, path, header_file)
            if last_rc == 0:
                elapsed = time.time() - t0 or 1.0
                print(
                    f"    done {path.name} in {elapsed / 60:.1f} min "
                    f"({size / elapsed / 1024**2:.1f} MiB/s)",
                    flush=True,
                )
                return
            wait = min(60 * attempt, 300)
            print(
                f"    curl rc={last_rc} for {path.name} "
                f"(attempt {attempt}/8); sleep {wait}s",
                flush=True,
            )
            time.sleep(wait)
        raise RuntimeError(
            f"upload failed for {path.name} after 8 attempts (last curl rc={last_rc})"
        )
    finally:
        try:
            os.unlink(header_file)
        except OSError:
            pass


def upload_record(
    client: Client,
    manifest: dict,
    record_id: str,
    title: str,
    description: str,
    files: list[Path],
    *,
    dry_run: bool,
    create_only: bool,
    manifest_path: Path,
) -> None:
    missing = [p for p in files if not p.is_file()]
    if missing:
        raise SystemExit(f"{record_id}: missing {missing}")
    size_gb = sum(p.stat().st_size for p in files) / 1024**3
    if size_gb > MAX_RECORD_GB:
        raise SystemExit(f"{record_id}: {size_gb:.1f} GiB exceeds {MAX_RECORD_GB} GiB")
    existing = next((r for r in manifest["records"] if r.get("id") == record_id), None)
    if existing and existing.get("upload_complete"):
        print(f"{record_id}: already complete ({existing.get('prereserved_doi')})")
        return
    print(f"\n=== {record_id} ({size_gb:.2f} GiB, {len(files)} files) ===")
    if dry_run:
        for p in files:
            print(f"  would upload {p}")
        return
    if existing and existing.get("deposition_id"):
        dep_id = existing["deposition_id"]
        bucket = existing["bucket_url"]
        doi = existing.get("prereserved_doi")
        print(f"  resume {dep_id} DOI {doi}")
    else:
        dep = create_deposition(client, metadata_for(title, description))
        dep_id = dep["id"]
        bucket = dep["links"]["bucket"]
        doi = dep["metadata"].get("prereserve_doi", {}).get("doi")
        entry = {
            "id": record_id,
            "title": title,
            "deposition_id": dep_id,
            "draft_html": dep["links"]["latest_draft_html"],
            "prereserved_doi": doi,
            "bucket_url": bucket,
            "files": [p.name for p in files],
            "uploaded_files": [],
            "upload_complete": False,
        }
        manifest["records"] = [r for r in manifest["records"] if r.get("id") != record_id]
        manifest["records"].append(entry)
        save_manifest(manifest_path, manifest)
        print(f"  created {dep_id} DOI {doi}")
    if create_only:
        return
    uploaded = set(existing.get("uploaded_files", []) if existing else [])
    for path in files:
        if path.name in uploaded:
            print(f"    skip {path.name}")
            continue
        upload_file(client, bucket, path)
        uploaded.add(path.name)
        for rec in manifest["records"]:
            if rec.get("id") == record_id:
                rec["uploaded_files"] = sorted(uploaded)
                rec["upload_complete"] = len(uploaded) == len(files)
        save_manifest(manifest_path, manifest)
    for rec in manifest["records"]:
        if rec.get("id") == record_id:
            rec["upload_complete"] = len(uploaded) == len(files)
    save_manifest(manifest_path, manifest)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--packages-dir",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "release" / "zenodo" / "payloads",
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--create-only", action="store_true")
    parser.add_argument("--skip-03", action="store_true")
    args = parser.parse_args()
    packages = args.packages_dir
    manifest_path = packages / "zenodo_blocks_manifest.json"
    manifest = load_manifest(manifest_path)
    manifest.setdefault("community", COMMUNITY)
    manifest.setdefault("community_url", f"https://zenodo.org/communities/{COMMUNITY}/")
    if not args.dry_run:
        client = Client(load_token())
    else:
        client = Client(token="dry-run")

    for block in BLOCK_META:
        zip_path = packages / block["zip"]
        upload_record(
            client,
            manifest,
            block["id"],
            block["title"],
            block["description"],
            [zip_path],
            dry_run=args.dry_run,
            create_only=args.create_only,
            manifest_path=manifest_path,
        )
        save_manifest(manifest_path, manifest)

    if not args.skip_03:
        buckets_path = packages / "03_buckets.json"
        if not buckets_path.is_file():
            raise SystemExit(f"missing {buckets_path}; run pack_zenodo_blocks.sh first")
        buckets = json.loads(buckets_path.read_text())
        for bucket in buckets["records"]:
            files = [packages / item["zip"] for item in bucket["files"]]
            upload_record(
                client,
                manifest,
                f"03_docking_{bucket['id']}",
                f"AffiTox: docking outputs ({bucket['id']} of record 03)",
                (
                    "Frozen docking dumps for a subset of the 16-target panel. "
                    "Record 03 is split across sibling depositions because Zenodo "
                    "limits a dataset to 50 GiB. GPU re-runs of Boltz-2 / DynamicBind "
                    f"are not bit-identical. Code: {GITHUB}."
                ),
                files,
                dry_run=args.dry_run,
                create_only=args.create_only,
                manifest_path=manifest_path,
            )
            save_manifest(manifest_path, manifest)

    manifest["updated"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    save_manifest(manifest_path, manifest)
    print(f"\nManifest: {manifest_path}")
    for rec in manifest["records"]:
        status = "COMPLETE" if rec.get("upload_complete") else "IN PROGRESS"
        print(f"  {rec.get('id')}: {rec.get('prereserved_doi')} [{status}]")


if __name__ == "__main__":
    main()
