#!/usr/bin/env python3
"""Publish sibling Zenodo deposits for AffiTox record 03 shards filtered to N=10711.

Creates new depositions for 03a–03e, marks each current parent isObsoletedBy the
new DOI, and retargets umbrella 22007993 hasPart.

Token from ZENODO_ACCESS_TOKEN or ~/.zenodo/token. Never prints the token.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_POST = Path(__file__).resolve().parent
ROOT = _POST.resolve().parents[1]
if str(_POST) not in sys.path:
    sys.path.insert(0, str(_POST))
from upload_zenodo import GITHUB, Client, load_token  # noqa: E402
from upload_zenodo_blocks import (  # noqa: E402
    COMMUNITY,
    create_deposition,
    metadata_for,
    upload_file,
)
from publish_zenodo_v1 import (  # noqa: E402
    UMBRELLA_DOI,
    UMBRELLA_ID,
    doi_rel,
    editable_metadata,
    enter_edit,
    get_dep,
    merge_related,
    put_and_publish,
)

SHARDS = {
    "03a": {
        "parent_id": 23224051,
        "pdbs": ["3eyg", "4f65", "5lf3", "2v5z"],
        "title": "AffiTox: docking outputs (03a of record 03: 3eyg, 4f65, 5lf3, 2v5z)",
    },
    "03b": {
        "parent_id": 23225395,
        "pdbs": ["3jy9", "3lxk"],
        "title": "AffiTox: docking outputs (03b of record 03: 3jy9, 3lxk)",
    },
    "03c": {
        "parent_id": 23225788,
        "pdbs": ["5mo4", "6gqj", "4zau"],
        "title": "AffiTox: docking outputs (03c of record 03: 5mo4, 6gqj, 4zau)",
    },
    "03d": {
        "parent_id": 23226603,
        "pdbs": ["7kk3", "5jkv", "4tz4"],
        "title": "AffiTox: docking outputs (03d of record 03: 7kk3, 5jkv, 4tz4)",
    },
    "03e": {
        "parent_id": 23227176,
        "pdbs": ["4ase", "6jok", "11ue", "1g5m"],
        "title": "AffiTox: docking outputs (03e of record 03: 4ase, 6jok, 11ue, 1g5m)",
    },
}

ZIP_DIR = ROOT / "release" / "zenodo" / "payloads" / "03"
REPORT_DIR = ROOT / "release" / "zenodo" / "payloads" / "03_panel10711"


def description_html(shard: str, parent_doi: str, pdbs: list[str]) -> str:
    pdb_s = ", ".join(pdbs)
    return (
        f"<p>AffiTox Zenodo record {shard} — frozen docking dumps for "
        f"<code>{pdb_s}</code>, filtered to the published panel "
        f"<strong>N = 10711</strong> (unique target–SMILES / curated ligand_id).</p>"
        "<p>Per-ligand poses/logs/JSON outside the curated id_maps are removed. "
        "GPU re-runs of Boltz-2 / DynamicBind are not bit-identical.</p>"
        f"<p>This deposit replaces "
        f"<a href='https://doi.org/{parent_doi}'>{parent_doi}</a>.</p>"
        f"<p>Index: <a href='https://doi.org/{UMBRELLA_DOI}'>{UMBRELLA_DOI}</a>. "
        f"Code: <a href='{GITHUB}'>{GITHUB}</a>.</p>"
    )


def restore_parent_if_editing(client: Client, parent_id: int) -> None:
    dep = get_dep(client, parent_id)
    if dep.get("state") != "inprogress":
        return
    print(f"republishing locked parent {parent_id} ...", flush=True)
    meta = editable_metadata(dict(dep.get("metadata") or {}))
    meta.pop("communities", None)
    put_and_publish(client, parent_id, meta)


def patch_parent_obsoleted(client: Client, parent_id: int, new_doi: str, shard: str) -> None:
    enter_edit(client, parent_id)
    dep = get_dep(client, parent_id)
    meta = editable_metadata(dict(dep.get("metadata") or {}))
    desc = str(meta.get("description") or "")
    banner = (
        "<p><strong>Superseded fileset.</strong> Use "
        f"<a href='https://doi.org/{new_doi}'>{new_doi}</a> for record {shard} "
        "(panel N=10711 docking dumps).</p>"
    )
    if new_doi not in desc:
        meta["description"] = banner + desc
    meta["notes"] = f"Superseded AffiTox {shard}. Current files: {new_doi}."
    meta["related_identifiers"] = merge_related(
        meta.get("related_identifiers") or [], doi_rel(new_doi, "isObsoletedBy")
    )
    meta.pop("communities", None)
    print(f"marking {parent_id} isObsoletedBy {new_doi} ...", flush=True)
    put_and_publish(client, parent_id, meta)


def patch_umbrella(client: Client, parent_doi: str, new_doi: str, shard: str) -> None:
    enter_edit(client, UMBRELLA_ID)
    dep = get_dep(client, UMBRELLA_ID)
    meta = editable_metadata(dict(dep.get("metadata") or {}))
    desc = str(meta.get("description") or "").replace(parent_doi, new_doi)
    old_id = parent_doi.rsplit(".", 1)[-1]
    new_id = new_doi.rsplit(".", 1)[-1]
    desc = desc.replace(old_id, new_id)
    if "10711" not in desc and "11180" in desc:
        desc = desc.replace("11180", "10711")
    meta["description"] = desc
    related = []
    replaced = False
    for item in meta.get("related_identifiers") or []:
        ident = str(item.get("identifier") or "")
        if item.get("relation") == "hasPart" and ident in {
            parent_doi,
            f"https://doi.org/{parent_doi}",
        }:
            related.append(doi_rel(new_doi, "hasPart"))
            replaced = True
            continue
        related.append(item)
    if not replaced:
        related = merge_related(related, doi_rel(new_doi, "hasPart"))
    meta["related_identifiers"] = related
    meta.pop("communities", None)
    print(f"updating umbrella hasPart {shard} -> {new_doi} ...", flush=True)
    put_and_publish(client, UMBRELLA_ID, meta)


def publish_shard(client: Client, shard: str) -> dict:
    cfg = SHARDS[shard]
    parent_id = int(cfg["parent_id"])
    parent_doi = f"10.5281/zenodo.{parent_id}"
    zips = [ZIP_DIR / f"results_{pdb}.zip" for pdb in cfg["pdbs"]]
    missing = [str(z) for z in zips if not z.is_file()]
    if missing:
        raise SystemExit(f"{shard}: missing zips: {missing}")
    restore_parent_if_editing(client, parent_id)

    meta = metadata_for(cfg["title"], description_html(shard, parent_doi, cfg["pdbs"]))
    meta["related_identifiers"] = [
        {
            "identifier": GITHUB,
            "relation": "isDocumentedBy",
            "resource_type": "software",
        },
        doi_rel(UMBRELLA_DOI, "isPartOf"),
        doi_rel(parent_doi, "isNewVersionOf"),
    ]
    meta["notes"] = (
        f"Current AffiTox {shard} docking dumps (panel N=10711). "
        f"Replaces {parent_doi}. Community {COMMUNITY}."
    )
    print(f"[{shard}] creating deposition ...", flush=True)
    dep = create_deposition(client, meta)
    dep_id = int(dep["id"])
    bucket = (dep.get("links") or {}).get("bucket")
    print(f"  draft {dep_id}", flush=True)
    if not bucket:
        raise SystemExit(f"{dep_id} has no bucket")
    for zpath in zips:
        print(
            f"  uploading {zpath.name} ({zpath.stat().st_size / 1024**3:.2f} GiB) ...",
            flush=True,
        )
        upload_file(client, bucket, zpath)
    dep = get_dep(client, dep_id)
    if not (dep.get("files") or []):
        raise SystemExit(f"{dep_id} has no files after upload")
    print(f"  publishing {dep_id} ...", flush=True)
    published = put_and_publish(
        client, dep_id, editable_metadata(dict(dep.get("metadata") or {}))
    )
    pub_meta = published.get("metadata") or {}
    doi = pub_meta.get("doi") or f"10.5281/zenodo.{dep_id}"
    conceptdoi = published.get("conceptdoi") or pub_meta.get("conceptdoi")
    print(f"  published {doi}", flush=True)
    patch_parent_obsoleted(client, parent_id, doi, shard)
    patch_umbrella(client, parent_doi, doi, shard)
    return {
        "shard": shard,
        "parent_id": parent_id,
        "parent_doi": parent_doi,
        "new_id": dep_id,
        "doi": doi,
        "conceptdoi": conceptdoi,
        "files": [z.name for z in zips],
        "bytes": sum(z.stat().st_size for z in zips),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "shards",
        nargs="*",
        default=None,
        help="Shards to publish (default: all 03a–03e)",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    shards = args.shards or list(SHARDS)
    bad = [s for s in shards if s not in SHARDS]
    if bad:
        raise SystemExit(f"unknown shards: {bad}; choose from {sorted(SHARDS)}")
    if args.dry_run:
        for shard in shards:
            cfg = SHARDS[shard]
            zips = [ZIP_DIR / f"results_{pdb}.zip" for pdb in cfg["pdbs"]]
            total = sum(z.stat().st_size for z in zips if z.is_file())
            print(
                f"{shard}: parent={cfg['parent_id']} "
                f"zips={sum(z.is_file() for z in zips)}/{len(zips)} "
                f"bytes={total/1024**3:.2f} GiB"
            )
        return 0
    client = Client(load_token())
    results = [publish_shard(client, shard) for shard in shards]
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    report = REPORT_DIR / "zenodo_03_panel_10711.json"
    report.write_text(json.dumps(results, indent=2) + "\n")
    print(f"wrote {report}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
