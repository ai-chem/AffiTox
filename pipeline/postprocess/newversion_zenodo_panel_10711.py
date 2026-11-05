#!/usr/bin/env python3
"""Publish sibling Zenodo deposits for AffiTox panel N=10711 (records 02/04/05/06).

Creates new depositions (bucket-locked parents cannot use newversion), marks each
current parent isObsoletedBy the new DOI, and retargets umbrella 22007993 hasPart.

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

# Current published parents (as of docs/zenodo/README.md before this refresh).
RECORDS = {
    "02": {
        "parent_id": 23213816,
        "zip": "02_prepared.zip",
        "title": "AffiTox: prepared inputs, ligand_id maps, and method weights (record 02)",
        "why": "Panel curated to N=10711 (unique target–SMILES); id_maps filtered to match.",
    },
    "04": {
        "parent_id": 23214418,
        "zip": "04_merged.zip",
        "title": "AffiTox: merged ligand_id tables (record 04)",
        "why": "Merged docking tables for the published N=10711 panel.",
    },
    "05": {
        "parent_id": 23979584,
        "zip": "05_posebusters.zip",
        "title": "AffiTox: PoseBusters tables (record 05)",
        "why": "PoseBusters CSVs filtered/reaggregated to the N=10711 panel.",
    },
    "06": {
        "parent_id": 23979728,
        "zip": "06_metrics_figures.zip",
        "title": "AffiTox: metrics tables and PNG figures (record 06)",
        "why": "Metrics/figures rebuilt for the N=10711 panel (21-check PoseBusters curves).",
    },
}


def description_html(rec_id: str, parent_doi: str, why: str) -> str:
    return (
        f"<p>AffiTox Zenodo record {rec_id} for the canonical 16-target panel "
        f"(<code>analysis/excluding_2z5x_3mjg</code>), panel N = <strong>10711</strong>.</p>"
        f"<p>{why}</p>"
        f"<p>This deposit replaces "
        f"<a href='https://doi.org/{parent_doi}'>{parent_doi}</a>.</p>"
        f"<p>Index: <a href='https://doi.org/{UMBRELLA_DOI}'>{UMBRELLA_DOI}</a>. "
        f"Code: <a href='{GITHUB}'>{GITHUB}</a>.</p>"
    )


def restore_parent_if_editing(client: Client, parent_id: int) -> None:
    dep = get_dep(client, parent_id)
    if dep.get("state") != "inprogress":
        return
    print(f"republishing locked parent {parent_id} (old zip kept) ...", flush=True)
    meta = editable_metadata(dict(dep.get("metadata") or {}))
    meta.pop("communities", None)
    put_and_publish(client, parent_id, meta)
    print("  parent restored", flush=True)


def patch_parent_obsoleted(client: Client, parent_id: int, new_doi: str, rec_id: str) -> None:
    enter_edit(client, parent_id)
    dep = get_dep(client, parent_id)
    meta = editable_metadata(dict(dep.get("metadata") or {}))
    desc = str(meta.get("description") or "")
    banner = (
        "<p><strong>Superseded fileset.</strong> Use "
        f"<a href='https://doi.org/{new_doi}'>{new_doi}</a> for record {rec_id} "
        "(panel N=10711).</p>"
    )
    if new_doi not in desc:
        meta["description"] = banner + desc
    meta["notes"] = f"Superseded AffiTox record {rec_id} zip. Current files: {new_doi}."
    meta["related_identifiers"] = merge_related(
        meta.get("related_identifiers") or [], doi_rel(new_doi, "isObsoletedBy")
    )
    meta.pop("communities", None)
    print(f"marking {parent_id} isObsoletedBy {new_doi} ...", flush=True)
    put_and_publish(client, parent_id, meta)
    print("  parent metadata published", flush=True)


def patch_umbrella(client: Client, parent_doi: str, new_doi: str, rec_id: str) -> None:
    enter_edit(client, UMBRELLA_ID)
    dep = get_dep(client, UMBRELLA_ID)
    meta = editable_metadata(dict(dep.get("metadata") or {}))
    desc = str(meta.get("description") or "").replace(parent_doi, new_doi)
    old_id = parent_doi.rsplit(".", 1)[-1]
    new_id = new_doi.rsplit(".", 1)[-1]
    desc = desc.replace(old_id, new_id)
    # Explicit N mention
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
    print(f"updating umbrella hasPart {rec_id} -> {new_doi} ...", flush=True)
    put_and_publish(client, UMBRELLA_ID, meta)
    print(f"  umbrella published {UMBRELLA_DOI}", flush=True)


def publish_one(client: Client, rec_id: str) -> dict:
    cfg = RECORDS[rec_id]
    parent_id = int(cfg["parent_id"])
    parent_doi = f"10.5281/zenodo.{parent_id}"
    zip_path = ROOT / "release" / "zenodo" / "payloads" / cfg["zip"]
    if not zip_path.is_file():
        raise SystemExit(f"missing {zip_path}; pack first")
    restore_parent_if_editing(client, parent_id)

    meta = metadata_for(cfg["title"], description_html(rec_id, parent_doi, cfg["why"]))
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
        f"Current AffiTox record {rec_id} (panel N=10711). "
        f"Replaces {parent_doi}. Community {COMMUNITY}."
    )
    print(f"[{rec_id}] creating new deposition ...", flush=True)
    dep = create_deposition(client, meta)
    dep_id = int(dep["id"])
    bucket = (dep.get("links") or {}).get("bucket")
    doi_pre = (dep.get("metadata") or {}).get("prereserve_doi", {}).get("doi")
    print(f"  draft {dep_id} prereserved {doi_pre}", flush=True)
    if not bucket:
        raise SystemExit(f"{dep_id} has no bucket")
    print(
        f"  uploading {zip_path.name} ({zip_path.stat().st_size / 1024**2:.1f} MiB) ...",
        flush=True,
    )
    upload_file(client, bucket, zip_path)
    dep = get_dep(client, dep_id)
    if not (dep.get("files") or []):
        raise SystemExit(f"{dep_id} has no files after upload; refuse to publish")
    print(f"  publishing {dep_id} ...", flush=True)
    published = put_and_publish(
        client, dep_id, editable_metadata(dict(dep.get("metadata") or {}))
    )
    pub_meta = published.get("metadata") or {}
    doi = pub_meta.get("doi") or f"10.5281/zenodo.{dep_id}"
    conceptdoi = published.get("conceptdoi") or pub_meta.get("conceptdoi")
    print(f"  published {doi} (concept {conceptdoi})", flush=True)

    patch_parent_obsoleted(client, parent_id, doi, rec_id)
    patch_umbrella(client, parent_doi, doi, rec_id)

    return {
        "record": rec_id,
        "parent_id": parent_id,
        "parent_doi": parent_doi,
        "new_id": dep_id,
        "doi": doi,
        "conceptdoi": conceptdoi,
        "zip": zip_path.name,
        "zip_bytes": zip_path.stat().st_size,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "records",
        nargs="*",
        default=["02", "04", "05", "06"],
        choices=sorted(RECORDS),
        help="Which blocks to publish (default: 02 04 05 06)",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.dry_run:
        for rid in args.records:
            z = ROOT / "release" / "zenodo" / "payloads" / RECORDS[rid]["zip"]
            print(f"{rid}: parent={RECORDS[rid]['parent_id']} zip={z} exists={z.is_file()}")
        return 0
    client = Client(load_token())
    report_path = ROOT / "release" / "zenodo" / "payloads" / "zenodo_panel_10711.json"
    results = []
    for rid in args.records:
        results.append(publish_one(client, rid))
    report_path.write_text(json.dumps(results, indent=2) + "\n")
    print(f"wrote {report_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
