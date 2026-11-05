#!/usr/bin/env python3
"""Publish AffiTox v1.0.0 drafts, restrict historical records, fix umbrella 22007993.

Token from ZENODO_ACCESS_TOKEN or ~/.zenodo/token. Never prints the token.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

_POST = Path(__file__).resolve().parent
ROOT = _POST.resolve().parents[1]
if str(_POST) not in sys.path:
    sys.path.insert(0, str(_POST))
from upload_zenodo import API, GITHUB, Client, load_token  # noqa: E402
from upload_zenodo_blocks import upload_file  # noqa: E402

UMBRELLA_ID = 22007993
UMBRELLA_DOI = "10.5281/zenodo.22007993"
CONCEPT_DOI = "10.5281/zenodo.22007992"
HISTORICAL_IDS = (20825057, 20825059, 20825061, 20825063, 20825065, 20825067)
PANEL = (
    "1g5m, 2v5z, 3eyg, 3jy9, 3lxk, 11ue, 4ase, 4f65, "
    "4tz4, 4zau, 5jkv, 5mo4, 6gqj, 6jok, 5lf3, 7kk3"
)
NEW_RECORDS = [
    {
        "id": "01_raw",
        "deposition_id": 23213812,
        "title": "AffiTox: raw ChEMBL extracts and PDB structures (record 01)",
        "pdbs": PANEL,
    },
    {
        "id": "02_prepared",
        "deposition_id": 23979597,
        "title": "AffiTox: prepared inputs, ligand_id maps, and method weights (record 02)",
        "pdbs": PANEL,
    },
    {
        "id": "03_docking_03a",
        "deposition_id": 23981992,
        "title": "AffiTox: docking outputs (03a of record 03: 3eyg, 4f65, 5lf3, 2v5z; panel N=10711)",
        "pdbs": "3eyg, 4f65, 5lf3, 2v5z",
    },
    {
        "id": "03_docking_03b",
        "deposition_id": 23982808,
        "title": "AffiTox: docking outputs (03b of record 03: 3jy9, 3lxk; panel N=10711)",
        "pdbs": "3jy9, 3lxk",
    },
    {
        "id": "03_docking_03c",
        "deposition_id": 23983375,
        "title": "AffiTox: docking outputs (03c of record 03: 5mo4, 6gqj, 4zau; panel N=10711)",
        "pdbs": "5mo4, 6gqj, 4zau",
    },
    {
        "id": "03_docking_03d",
        "deposition_id": 23984015,
        "title": "AffiTox: docking outputs (03d of record 03: 7kk3, 5jkv, 4tz4; panel N=10711)",
        "pdbs": "7kk3, 5jkv, 4tz4",
    },
    {
        "id": "03_docking_03e",
        "deposition_id": 23984506,
        "title": "AffiTox: docking outputs (03e of record 03: 4ase, 6jok, 11ue, 1g5m; panel N=10711)",
        "pdbs": "4ase, 6jok, 11ue, 1g5m",
    },
    {
        "id": "04_merged",
        "deposition_id": 23979579,
        "title": "AffiTox: merged ligand_id tables (record 04)",
        "pdbs": PANEL,
    },
    {
        "id": "05_posebusters",
        "deposition_id": 23979721,
        "title": "AffiTox: PoseBusters tables (record 05)",
        "pdbs": PANEL,
    },
    {
        "id": "06_metrics_figures",
        "deposition_id": 23979767,
        "title": "AffiTox: metrics tables and PNG figures (record 06)",
        "pdbs": PANEL,
    },
]
ACCESS_CONDITIONS = (
    "<p>These files are a <strong>historical, superseded</strong> AffiTox dump "
    "(old per-protein zips, panel included 7awe/2z5x/3mjg). They are not the "
    f"v1.0.0 dataset. Cite <a href='https://doi.org/{UMBRELLA_DOI}'>{UMBRELLA_DOI}</a> "
    f"and <a href='{GITHUB}'>{GITHUB}</a>. Access to these files is granted on request "
    "to the depositor for provenance only.</p>"
)


def _json(client: Client, method: str, url: str, payload: dict | None = None) -> dict:
    body = None
    headers = {}
    if payload is not None:
        body = json.dumps(payload).encode()
        headers["Content-Type"] = "application/json"
    _, raw = client.request(method, url, data=body, headers=headers)
    return json.loads(raw) if raw else {}


def get_dep(client: Client, dep_id: int) -> dict:
    return _json(client, "GET", f"{API}/deposit/depositions/{dep_id}")


def editable_metadata(meta: dict) -> dict:
    out = dict(meta)
    for key in ("doi", "prereserve_doi"):
        out.pop(key, None)
    return out


def merge_related(existing: list, extra: dict) -> list:
    items = list(existing or [])
    key = (extra.get("identifier"), extra.get("relation"))
    if key not in {(i.get("identifier"), i.get("relation")) for i in items}:
        items.append(extra)
    return items


def drop_related(existing: list, predicate) -> list:
    return [i for i in (existing or []) if not predicate(i)]


def doi_rel(identifier: str, relation: str, resource_type: str = "dataset") -> dict:
    return {
        "identifier": identifier,
        "relation": relation,
        "resource_type": resource_type,
        "scheme": "doi",
    }


def enter_edit(client: Client, dep_id: int) -> None:
    try:
        _json(client, "POST", f"{API}/deposit/depositions/{dep_id}/actions/edit")
    except Exception as exc:
        if "400" not in str(exc) and "409" not in str(exc):
            raise


def put_and_publish(client: Client, dep_id: int, meta: dict) -> dict:
    _json(client, "PUT", f"{API}/deposit/depositions/{dep_id}", {"metadata": meta})
    return _json(client, "POST", f"{API}/deposit/depositions/{dep_id}/actions/publish")


def refresh_06(client: Client) -> None:
    zip_path = ROOT / "release" / "zenodo" / "payloads" / "06_metrics_figures.zip"
    if not zip_path.is_file():
        raise SystemExit(f"missing {zip_path}")
    dep = get_dep(client, 23979767)
    if dep.get("submitted"):
        print("06 already published; skip zip replace")
        return
    bucket = (dep.get("links") or {}).get("bucket")
    if not bucket:
        raise SystemExit("06 draft has no bucket URL")
    print(f"uploading refreshed 06 ({zip_path.stat().st_size / 1024**2:.1f} MiB)", flush=True)
    upload_file(client, bucket, zip_path)


def publish_new(client: Client) -> None:
    part_of = doi_rel(UMBRELLA_DOI, "isPartOf")
    github = {
        "identifier": GITHUB,
        "relation": "isDocumentedBy",
        "resource_type": "software",
    }
    for rec in NEW_RECORDS:
        dep_id = rec["deposition_id"]
        dep = get_dep(client, dep_id)
        files = dep.get("files") or []
        if not files:
            raise SystemExit(f"{dep_id} has no files; refuse to publish empty")
        if dep.get("submitted") and dep.get("state") == "done":
            print(f"{rec['id']}: already published")
            continue
        meta = editable_metadata(dict(dep.get("metadata") or {}))
        meta["title"] = rec["title"]
        meta["related_identifiers"] = merge_related(
            merge_related(meta.get("related_identifiers") or [], github), part_of
        )
        if rec["id"].startswith("03_"):
            meta["description"] = (
                "<p>Frozen docking dumps for "
                f"{rec['pdbs']}. Record 03 of AffiTox v1.0.0 is split across "
                "sibling depositions 03a–03e because Zenodo limits a dataset to "
                "50 GiB. GPU re-runs of Boltz-2 / DynamicBind are not bit-identical. "
                f"Index: <a href='https://doi.org/{UMBRELLA_DOI}'>{UMBRELLA_DOI}</a>. "
                f"Code: <a href='{GITHUB}'>{GITHUB}</a>.</p>"
            )
        print(f"publishing {rec['id']} {dep_id} ({len(files)} files) ...", flush=True)
        put_and_publish(client, dep_id, meta)
        print(f"  published 10.5281/zenodo.{dep_id}", flush=True)


def patch_umbrella(client: Client) -> None:
    enter_edit(client, UMBRELLA_ID)
    dep = get_dep(client, UMBRELLA_ID)
    meta = editable_metadata(dict(dep.get("metadata") or {}))
    rows = "".join(
        f"<tr><td>{rec['id']}</td><td><a href='https://doi.org/10.5281/zenodo.{rec['deposition_id']}'>"
        f"10.5281/zenodo.{rec['deposition_id']}</a></td><td>{rec['title']}</td></tr>"
        for rec in NEW_RECORDS
    )
    meta["title"] = "AffiTox v1.0.0: data release manifest"
    meta["upload_type"] = "dataset"
    meta["access_right"] = "open"
    meta["license"] = meta.get("license") or "cc-by-4.0"
    meta["creators"] = [
        {"name": "Konovalova, Olga A.", "orcid": "0000-0003-0391-3211"}
    ]
    meta["keywords"] = ["docking", "toxicity", "benchmark", "ChEMBL", "AffiTox"]
    meta["notes"] = (
        "Citable index for AffiTox v1.0.0 data. Concept DOI "
        f"{CONCEPT_DOI}. Historical per-protein zips 10.5281/zenodo.20825057–067 "
        "are superseded and restricted; do not cite them as the paper dataset."
    )
    meta["description"] = (
        "<p>Index record for the AffiTox v1.0.0 dataset (canonical 16-target panel: "
        f"{PANEL}). Cite this DOI in the paper. Files of the six information blocks "
        "live in the related records below; record 03 is sharded 03a–03e (Zenodo "
        "50 GiB/record limit). GPU re-runs of Boltz-2 / DynamicBind are not "
        "bit-identical; record 03 is the frozen snapshot.</p>"
        f"<p>Code: <a href='{GITHUB}'>{GITHUB}</a>.</p>"
        "<table><thead><tr><th>Block</th><th>DOI</th><th>Title</th></tr></thead>"
        f"<tbody>{rows}</tbody></table>"
        "<p><strong>Historical.</strong> Do not use "
        "10.5281/zenodo.20825057–067 (old per-protein zips including 7awe/2z5x/3mjg). "
        "Those records are restricted and marked isObsoletedBy this DOI.</p>"
    )
    related = drop_related(
        meta.get("related_identifiers") or [],
        lambda i: i.get("relation") in {"hasPart", "obsoletes", "isDocumentedBy"}
        or str(i.get("identifier") or "").startswith("10.5281/zenodo.208250"),
    )
    related = merge_related(
        related,
        {
            "identifier": GITHUB,
            "relation": "isDocumentedBy",
            "resource_type": "software",
        },
    )
    for rec in NEW_RECORDS:
        related = merge_related(
            related, doi_rel(f"10.5281/zenodo.{rec['deposition_id']}", "hasPart")
        )
    for hid in HISTORICAL_IDS:
        related = merge_related(related, doi_rel(f"10.5281/zenodo.{hid}", "obsoletes"))
    meta["related_identifiers"] = related
    # published umbrella: communities often 400 on update
    meta.pop("communities", None)
    print(f"updating umbrella {UMBRELLA_ID} ...", flush=True)
    put_and_publish(client, UMBRELLA_ID, meta)
    print(f"  umbrella published {UMBRELLA_DOI}", flush=True)


def restrict_historical(client: Client) -> None:
    github = {
        "identifier": GITHUB,
        "relation": "isDocumentedBy",
        "resource_type": "software",
    }
    obsolete = doi_rel(UMBRELLA_DOI, "isObsoletedBy")
    for dep_id in HISTORICAL_IDS:
        enter_edit(client, dep_id)
        dep = get_dep(client, dep_id)
        meta = editable_metadata(dict(dep.get("metadata") or {}))
        desc = str(meta.get("description") or "")
        if "HISTORICAL" not in desc:
            desc = (
                "<p><strong>HISTORICAL — not AffiTox v1.0.0.</strong> "
                "Old per-protein docking zips (panel included 7awe/2z5x/3mjg). "
                f"Cite <a href='https://doi.org/{UMBRELLA_DOI}'>{UMBRELLA_DOI}</a>.</p>"
                + desc
            )
        meta["description"] = desc
        meta["access_right"] = "restricted"
        meta["access_conditions"] = ACCESS_CONDITIONS
        meta["notes"] = (
            "HISTORICAL. Files restricted. Not the AffiTox v1.0.0 data release. "
            f"Superseded by {UMBRELLA_DOI}. Current code: {GITHUB}."
        )
        related = merge_related(meta.get("related_identifiers") or [], github)
        related = merge_related(related, obsolete)
        meta["related_identifiers"] = related
        meta.pop("communities", None)
        print(f"restricting {dep_id} ...", flush=True)
        put_and_publish(client, dep_id, meta)
        print(f"  restricted 10.5281/zenodo.{dep_id}", flush=True)


def write_outputs(client: Client) -> None:
    out_dir = ROOT / "release" / "zenodo" / "payloads"
    rows = []
    published = []
    for rec in NEW_RECORDS:
        dep = get_dep(client, rec["deposition_id"])
        meta = dep.get("metadata") or {}
        doi = meta.get("doi") or f"10.5281/zenodo.{rec['deposition_id']}"
        published.append(
            {
                "id": rec["id"],
                "deposition_id": rec["deposition_id"],
                "doi": doi,
                "conceptdoi": meta.get("conceptdoi"),
                "state": dep.get("state"),
                "access_right": meta.get("access_right"),
                "title": meta.get("title"),
                "files": [
                    {
                        "filename": f.get("filename") or f.get("key"),
                        "size": f.get("filesize") or f.get("size"),
                    }
                    for f in dep.get("files") or []
                ],
            }
        )
        for f in dep.get("files") or []:
            rows.append(
                {
                    "record": rec["id"],
                    "zenodo_record": rec["deposition_id"],
                    "doi": doi,
                    "record_title": meta.get("title"),
                    "filename": f.get("filename") or f.get("key"),
                    "size_bytes": f.get("filesize") or f.get("size"),
                }
            )
    umb = get_dep(client, UMBRELLA_ID)
    umeta = umb.get("metadata") or {}
    report = {
        "umbrella_doi": umeta.get("doi") or UMBRELLA_DOI,
        "umbrella_conceptdoi": umeta.get("conceptdoi") or CONCEPT_DOI,
        "umbrella_access_right": umeta.get("access_right"),
        "records": published,
        "historical": [],
    }
    for hid in HISTORICAL_IDS:
        dep = get_dep(client, hid)
        meta = dep.get("metadata") or {}
        report["historical"].append(
            {
                "id": hid,
                "doi": meta.get("doi") or f"10.5281/zenodo.{hid}",
                "state": dep.get("state"),
                "access_right": meta.get("access_right"),
                "title": meta.get("title"),
            }
        )
    (out_dir / "zenodo_published_v1.json").write_text(json.dumps(report, indent=2) + "\n")
    tsv_path = ROOT / "docs" / "zenodo" / "affitox_data_manifest.tsv"
    lines = ["record\tzenodo_record\tdoi\trecord_title\tfilename\tsize_bytes"]
    for row in rows:
        lines.append(
            f"{row['record']}\t{row['zenodo_record']}\t{row['doi']}\t"
            f"{row['record_title']}\t{row['filename']}\t{row['size_bytes']}"
        )
    tsv_path.write_text("\n".join(lines) + "\n")
    print(f"wrote {out_dir / 'zenodo_published_v1.json'}")
    print(f"wrote {tsv_path}")
    for rec in published:
        print(f"  {rec['id']}: {rec['doi']} [{rec['state']}] {rec['access_right']}")
    for rec in report["historical"]:
        print(f"  historical {rec['doi']}: {rec['access_right']} [{rec['state']}]")


def main() -> int:
    client = Client(load_token())
    refresh_06(client)
    publish_new(client)
    patch_umbrella(client)
    restrict_historical(client)
    write_outputs(client)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
