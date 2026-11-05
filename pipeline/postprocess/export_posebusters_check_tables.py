#!/usr/bin/env python3
"""Export PoseBusters check-level pass rates as CSV and SI landscape tables."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path("/mnt/tank/scratch/okonovalova/AffiTox")
SRC = (
    ROOT
    / "analysis"
    / "excluding_2z5x_3mjg"
    / "tables"
    / "posebusters"
    / "pass_rates_by_check_target_method.csv"
)
OUT = (
    ROOT
    / "analysis"
    / "excluding_2z5x_3mjg"
    / "tables"
    / "posebusters"
    / "check_level"
)

TARGETS = [
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
]
METHODS = [
    ("boltz2", "Boltz-2"),
    ("dynamicbind", "DynamicBind"),
    ("gnina", "GNINA 1.3"),
    ("qvina", "QVina2"),
]


def pretty_check(name: str) -> str:
    return (
        name.replace("minimum_distance_to_", "min_dist_")
        .replace("volume_overlap_with_", "overlap_")
        .replace("protein-ligand_", "prot_lig_")
    )


def latex_escape(text: str) -> str:
    return text.replace("_", r"\_")


def wide_frame(sub: pd.DataFrame, checks: list[str]) -> pd.DataFrame:
    body = sub.pivot(index="check", columns="target", values="pass_rate_pct")
    body = body.reindex(index=checks, columns=TARGETS)
    body["pooled"] = body.mean(axis=1, skipna=True)
    all_tests = body.mean(axis=0, skipna=True).to_frame().T
    all_tests.index = ["All tests"]
    wide = pd.concat([body, all_tests])
    wide.index.name = "check"
    return wide.reset_index()


def write_tex(path: Path, wide: pd.DataFrame, method_label: str) -> None:
    cols = TARGETS + ["pooled"]
    headers = [t.upper() for t in TARGETS] + ["Pooled"]
    lines = [
        r"\setlength{\tabcolsep}{2.6pt}",
        r"\renewcommand{\arraystretch}{1.05}",
        r"\resizebox{\linewidth}{!}{%",
        r"\begin{tabular}{l" + "r" * len(cols) + "}",
        r"\toprule",
        "Check & " + " & ".join(headers) + r" \\",
        r"\midrule",
    ]
    for _, row in wide.iterrows():
        check = str(row["check"])
        values = " & ".join(f"{float(row[c]):.1f}" for c in cols)
        shown = "All tests" if check == "All tests" else pretty_check(check)
        prefix = r"\midrule" + "\n" if check == "All tests" else ""
        lines.append(f"{prefix}{latex_escape(shown)} & {values} \\\\")
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}%",
            r"}",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(SRC)
    raw["pass_rate_pct"] = 100.0 * raw["pass_rate"]
    first = raw[raw["method_id"] == "boltz2"]
    checks = first.drop_duplicates("check")["check"].tolist()

    long_rows = []
    for method_id, method_label in METHODS:
        sub = raw[raw["method_id"] == method_id].copy()
        wide = wide_frame(sub, checks)
        wide.insert(1, "method_id", method_id)
        wide.insert(2, "method_label", method_label)
        csv_path = OUT / f"check_level_{method_id}.csv"
        wide.to_csv(csv_path, index=False, float_format="%.4f")
        write_tex(OUT / f"check_level_{method_id}.tex", wide, method_label)
        for check in checks + ["All tests"]:
            row = wide.loc[wide["check"] == check].iloc[0]
            for target in TARGETS + ["pooled"]:
                long_rows.append(
                    {
                        "method_id": method_id,
                        "method_label": method_label,
                        "check": check,
                        "check_label": pretty_check(check) if check != "All tests" else "All tests",
                        "target": target,
                        "pass_rate_pct": float(row[target]),
                    }
                )
        print(csv_path)

    long_df = pd.DataFrame(long_rows)
    long_path = OUT / "check_level_long.csv"
    long_df.to_csv(long_path, index=False, float_format="%.4f")
    print(long_path)


if __name__ == "__main__":
    main()
