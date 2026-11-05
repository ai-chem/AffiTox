#!/usr/bin/env python3
"""Rename Zenodo community/deposits from ToxDock-Bench to AffiTox."""

from __future__ import annotations

import json
import os
import urllib.error
import urllib.request
from pathlib import Path

GITHUB = "https://github.com/ai-chem/AffiTox"
COMMUNITY = "affitox"

TITLE_MAP = {
    20825057: "AffiTox: raw input and docking results (part 1/6)",
    20825059: "AffiTox: raw docking results (part 2/6)",
    20825061: "AffiTox: raw docking results (part 3/6)",
    20825063: "AffiTox: raw docking results (part 4/6)",
    20825065: "AffiTox: raw docking results (part 5/6, 3jy9)",
    20825067: "AffiTox: raw docking results (part 6/6, 3eyg)",
}

DESC_MAP = {
    20825057: (
        "Raw BindingDB input snapshot and per-target raw docking outputs for "
        "7awe, 3mjg, 1g5m, 4f65, 4zau, 4tz4, 2z5x. Part 1/6 of the AffiTox "
        f"16-target toxicity docking benchmark. Analysis tables/figures: {GITHUB}."
    ),
    20825059: (
        "Per-target raw docking outputs for 6jok, 3lxk, 5jkv. "
        f"Part 2/6 of the AffiTox benchmark; analysis: {GITHUB}."
    ),
    20825061: (
        "Per-target raw docking outputs for 4ase, 6gqj. "
        f"Part 3/6 of the AffiTox benchmark; analysis: {GITHUB}."
    ),
    20825063: (
        "Per-target raw docking outputs for 7kk3, 5mo4. "
        f"Part 4/6 of the AffiTox benchmark; analysis: {GITHUB}."
    ),
    20825065: (
        "Per-target raw docking output for PDB target 3jy9. "
        f"Part 5/6 of the AffiTox benchmark; analysis: {GITHUB}."
    ),
    20825067: (
        "Per-target raw docking output for PDB target 3eyg. "
        f"Part 6/6 of the AffiTox benchmark; analysis: {GITHUB}."
    ),
}


def load_token() -> str:
    if os.environ.get("ZENODO_ACCESS_TOKEN"):
        return os.environ["ZENODO_ACCESS_TOKEN"]
    path = Path.home() / ".zenodo" / "token"
    return path.read_text().strip()


def api_request(token: str, method: str, url: str, data: bytes | None = None) -> dict:
    headers = {"Authorization": f"Bearer {token}"}
    if data is not None:
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, method=method, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            return json.loads(resp.read())
    except urllib.error.HTTPError as exc:
        body = exc.read()[:800]
        raise RuntimeError(f"{method} {url} -> {exc.code}: {body!r}") from exc


def update_depositions(token: str) -> None:
    for dep_id, new_title in TITLE_MAP.items():
        meta_patch = {
            "title": new_title,
            "description": DESC_MAP[dep_id],
            "keywords": ["docking", "toxicity", "benchmark", "BindingDB", "AffiTox"],
            "communities": [{"identifier": COMMUNITY}],
            "notes": (
                "Draft upload for the AffiTox Zenodo community (affitox). "
                "Analysis artifacts (tables/figures) are distributed via GitHub, not this record."
            ),
        }
        updated = api_request(
            token,
            "PUT",
            f"https://zenodo.org/api/deposit/depositions/{dep_id}",
            json.dumps({"metadata": meta_patch}).encode(),
        )
        print(f"OK {dep_id}: {updated['metadata']['title']} ({len(updated.get('files', []))} files)")


def main() -> None:
    token = load_token()
    update_depositions(token)
    recs = api_request(
        token,
        "GET",
        f"https://zenodo.org/api/communities/{COMMUNITY}/records?size=20",
    )
    print(f"\nCommunity: https://zenodo.org/communities/{COMMUNITY}/")
    print(f"Records in community: {recs['hits']['total']}")
    for hit in recs["hits"]["hits"]:
        print(f"  - {hit['metadata']['title']}")


if __name__ == "__main__":
    main()
