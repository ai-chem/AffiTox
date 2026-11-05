#!/usr/bin/env python3
"""Upload AffiTox Zenodo packages as draft depositions (no analysis zip)."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

API = os.environ.get("ZENODO_API", "https://zenodo.org/api")
COMMUNITY = os.environ.get("ZENODO_COMMUNITY", "affitox")
GITHUB = "https://github.com/ai-chem/AffiTox"
MAX_RECORD_GB = 50.0

# Bundles stay under Zenodo 50 GiB/record limit (Jun 2025 packages).
RECORD_BUNDLES: list[dict] = [
    {
        "part": 1,
        "title": "AffiTox: raw input and docking results (part 1/6)",
        "description": (
            "Raw BindingDB input snapshot and per-target raw docking outputs for "
            "7awe, 3mjg, 1g5m, 4f65, 4zau, 4tz4, 2z5x. Part 1/6 of the AffiTox "
            f"16-target toxicity docking benchmark. Analysis tables/figures: {GITHUB}."
        ),
        "files": [
            "toxdock-input-raw.zip",
            "results_by_target/results_7awe.zip",
            "results_by_target/results_3mjg.zip",
            "results_by_target/results_1g5m.zip",
            "results_by_target/results_4f65.zip",
            "results_by_target/results_4zau.zip",
            "results_by_target/results_4tz4.zip",
            "results_by_target/results_2z5x.zip",
        ],
    },
    {
        "part": 2,
        "title": "AffiTox: raw docking results (part 2/6)",
        "description": (
            "Per-target raw docking outputs for 6jok, 3lxk, 5jkv. "
            f"Part 2/6 of the AffiTox benchmark; analysis: {GITHUB}."
        ),
        "files": [
            "results_by_target/results_6jok.zip",
            "results_by_target/results_3lxk.zip",
            "results_by_target/results_5jkv.zip",
        ],
    },
    {
        "part": 3,
        "title": "AffiTox: raw docking results (part 3/6)",
        "description": (
            "Per-target raw docking outputs for 4ase, 6gqj. "
            f"Part 3/6 of the AffiTox benchmark; analysis: {GITHUB}."
        ),
        "files": [
            "results_by_target/results_4ase.zip",
            "results_by_target/results_6gqj.zip",
        ],
    },
    {
        "part": 4,
        "title": "AffiTox: raw docking results (part 4/6)",
        "description": (
            "Per-target raw docking outputs for 7kk3, 5mo4. "
            f"Part 4/6 of the AffiTox benchmark; analysis: {GITHUB}."
        ),
        "files": [
            "results_by_target/results_7kk3.zip",
            "results_by_target/results_5mo4.zip",
        ],
    },
    {
        "part": 5,
        "title": "AffiTox: raw docking results (part 5/6, 3jy9)",
        "description": (
            "Per-target raw docking output for PDB target 3jy9. "
            f"Part 5/6 of the AffiTox benchmark; analysis: {GITHUB}."
        ),
        "files": ["results_by_target/results_3jy9.zip"],
    },
    {
        "part": 6,
        "title": "AffiTox: raw docking results (part 6/6, 3eyg)",
        "description": (
            "Per-target raw docking output for PDB target 3eyg. "
            f"Part 6/6 of the AffiTox benchmark; analysis: {GITHUB}."
        ),
        "files": ["results_by_target/results_3eyg.zip"],
    },
]


@dataclass
class Client:
    token: str

    def request(
        self,
        method: str,
        url: str,
        data: bytes | None = None,
        headers: dict | None = None,
    ) -> tuple[int, bytes]:
        hdrs = {"Authorization": f"Bearer {self.token}"}
        if headers:
            hdrs.update(headers)
        req = urllib.request.Request(url, data=data, headers=hdrs, method=method)
        try:
            with urllib.request.urlopen(req, timeout=3600) as resp:
                return resp.status, resp.read()
        except urllib.error.HTTPError as exc:
            body = exc.read()
            raise RuntimeError(f"{method} {url} -> {exc.code}: {body[:500]!r}") from exc


def load_token() -> str:
    if os.environ.get("ZENODO_ACCESS_TOKEN"):
        return os.environ["ZENODO_ACCESS_TOKEN"]
    path = Path.home() / ".zenodo" / "token"
    if path.is_file():
        return path.read_text().strip()
    raise SystemExit("Zenodo token not found (~/.zenodo/token or ZENODO_ACCESS_TOKEN)")


def bundle_size_gb(packages: Path, files: list[str]) -> float:
    total = sum((packages / f).stat().st_size for f in files)
    return total / (1024**3)


def base_metadata(bundle: dict) -> dict:
    return {
        "title": bundle["title"],
        "upload_type": "dataset",
        "description": bundle["description"],
        "creators": [{"name": "Koont, Olga"}],
        "license": "cc-by-4.0",
        "keywords": [
            "docking",
            "toxicity",
            "benchmark",
            "BindingDB",
            "AffiTox",
        ],
        "communities": [{"identifier": COMMUNITY}],
        "related_identifiers": [
            {
                "identifier": GITHUB,
                "relation": "isDocumentedBy",
                "resource_type": "software",
            }
        ],
        "notes": (
            f"Draft upload for the AffiTox Zenodo community ({COMMUNITY}). "
            "Analysis artifacts (tables/figures) are distributed via GitHub, not this record."
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


def upload_file(client: Client, bucket_url: str, path: Path) -> None:
    name = path.name
    url = f"{bucket_url}/{name}"
    size = path.stat().st_size
    print(f"    uploading {name} ({size / 1024**3:.2f} GiB) ...", flush=True)
    t0 = time.time()
    subprocess.run(
        [
            "curl",
            "-fS",
            "--upload-file",
            str(path),
            "-H",
            f"Authorization: Bearer {client.token}",
            url,
        ],
        check=True,
    )
    elapsed = time.time() - t0
    rate = size / elapsed / (1024**2) if elapsed else 0
    print(f"    done {name} in {elapsed / 60:.1f} min ({rate:.1f} MiB/s)", flush=True)


def load_manifest(path: Path) -> dict:
    if path.is_file():
        return json.loads(path.read_text())
    return {"community": COMMUNITY, "community_url": f"https://zenodo.org/communities/{COMMUNITY}/", "records": []}


def save_manifest(path: Path, manifest: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--packages-dir",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "zenodo_packages",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="JSON manifest path (default: <packages-dir>/zenodo_upload_manifest.json)",
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--create-only",
        action="store_true",
        help="Create draft depositions without uploading files",
    )
    parser.add_argument(
        "--part",
        type=int,
        action="append",
        help="Upload only selected part number(s); repeatable",
    )
    args = parser.parse_args()

    packages = args.packages_dir
    manifest_path = args.manifest or (packages / "zenodo_upload_manifest.json")
    manifest = load_manifest(manifest_path)
    client = Client(load_token())

    selected = set(args.part) if args.part else None

    for bundle in RECORD_BUNDLES:
        part = bundle["part"]
        if selected and part not in selected:
            continue

        files = bundle["files"]
        missing = [f for f in files if not (packages / f).is_file()]
        if missing:
            raise SystemExit(f"Part {part}: missing files: {missing}")

        size_gb = bundle_size_gb(packages, files)
        if size_gb > MAX_RECORD_GB:
            raise SystemExit(f"Part {part}: {size_gb:.1f} GiB exceeds {MAX_RECORD_GB} GiB Zenodo limit")

        existing = next((r for r in manifest["records"] if r.get("part") == part), None)
        if existing and existing.get("upload_complete"):
            print(f"Part {part}: already complete (deposition {existing['deposition_id']}), skip")
            continue

        print(f"\n=== Part {part}/6 ({size_gb:.2f} GiB, {len(files)} files) ===")

        if args.dry_run:
            for f in files:
                print(f"  would upload: {f}")
            continue

        if existing and existing.get("deposition_id"):
            dep_id = existing["deposition_id"]
            bucket = existing["bucket_url"]
            draft_html = existing["draft_html"]
            doi = existing.get("prereserved_doi")
            print(f"  resume deposition {dep_id} ({draft_html})")
        else:
            dep = create_deposition(client, base_metadata(bundle))
            dep_id = dep["id"]
            bucket = dep["links"]["bucket"]
            draft_html = dep["links"]["latest_draft_html"]
            doi = dep["metadata"].get("prereserve_doi", {}).get("doi")
            entry = {
                "part": part,
                "title": bundle["title"],
                "deposition_id": dep_id,
                "conceptrecid": dep.get("conceptrecid"),
                "draft_html": draft_html,
                "prereserved_doi": doi,
                "bucket_url": bucket,
                "files": files,
                "upload_complete": False,
            }
            if existing:
                manifest["records"] = [r for r in manifest["records"] if r.get("part") != part]
            manifest["records"].append(entry)
            manifest["records"].sort(key=lambda r: r["part"])
            save_manifest(manifest_path, manifest)
            print(f"  created draft {dep_id} -> {draft_html}")
            if doi:
                print(f"  prereserved DOI: {doi}")

        if args.create_only:
            continue

        uploaded = set(existing.get("uploaded_files", []) if existing else [])
        for rel in files:
            if rel in uploaded:
                print(f"    skip (already uploaded): {Path(rel).name}")
                continue
            upload_file(client, bucket, packages / rel)
            uploaded.add(rel)
            for rec in manifest["records"]:
                if rec.get("part") == part:
                    rec["uploaded_files"] = sorted(uploaded)
            save_manifest(manifest_path, manifest)

        for rec in manifest["records"]:
            if rec.get("part") == part:
                rec["upload_complete"] = len(uploaded) == len(files)
        save_manifest(manifest_path, manifest)

    manifest["updated"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    save_manifest(manifest_path, manifest)
    print(f"\nManifest: {manifest_path}")
    print(f"Community: {manifest['community_url']}")
    for rec in sorted(manifest["records"], key=lambda r: r["part"]):
        status = "COMPLETE" if rec.get("upload_complete") else "IN PROGRESS"
        print(f"  part {rec['part']}: {rec['draft_html']} [{status}]")


if __name__ == "__main__":
    main()
