#!/usr/bin/env python3
"""
Standalone script to reproduce AffiTox nEF10.png as SVG.

Does not import or modify existing analysis modules. Logic mirrors
src/analysis/plots/enrichment.py::plot_nef_violins for frac=10.

Default input:
  analysis/tables/enrichment/ef_summary_all_proteins.csv
Default output:
  nEF10.svg (next to this repo's nEF10.png)

Usage:
  python scripts/plot_nef10_svg.py
  python scripts/plot_nef10_svg.py --input path/to/ef_summary_all_proteins.csv --output path/to/nEF10.svg
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

# Keep SVG text as editable text (not paths).
plt.rcParams["svg.fonttype"] = "none"

# ---------------------------------------------------------------------------
# Article style constants (mirrored from src/analysis/constants.py / style.py)
# ---------------------------------------------------------------------------

METHOD_ORDER: list[tuple[str, str]] = [
    ("boltz2", "Boltz-2"),
    ("dynamicbind", "DynamicBind"),
    ("gnina", "GNINA 1.3"),
    ("plapt", "PLAPT"),
    ("qvina", "QVina2"),
]

POSITIVE_SEQUENTIAL_STOPS: list[tuple[float, str]] = [
    (0.0, "#6584E1"),
    (0.35, "#B7CAEA"),
    (0.55, "#DEDAD7"),
    (0.75, "#FDBBAB"),
    (1.0, "#E5885F"),
]

RANDOM_SEED = 42
FRAC = 10


def positive_sequential_cmap(name: str = "nef_point_gradient") -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list(name, POSITIVE_SEQUENTIAL_STOPS, N=256)


def panel_plot(
    ax: plt.Axes,
    df: pd.DataFrame,
    col: str,
    title: str,
    panel_letter: str,
    order: list[str],
    cmap: LinearSegmentedColormap,
    seed: int,
    frac: int,
) -> None:
    """One panel: violins + jittered points + median diamonds + low-nEF labels."""
    rng = np.random.default_rng(seed)
    x_positions = np.arange(len(order))
    values_by_method = [
        df.loc[df["method_label"] == method, col].astype(float).dropna().to_numpy()
        for method in order
    ]

    violin = ax.violinplot(
        values_by_method,
        positions=x_positions,
        widths=0.82,
        showmeans=False,
        showmedians=False,
        showextrema=False,
    )
    for body in violin["bodies"]:
        body.set_facecolor("#DEDAD7")
        body.set_edgecolor("#2A2A2A")
        body.set_linewidth(0.6)
        body.set_alpha(0.22)

    label_thr = 0.6
    for i, vals in enumerate(values_by_method):
        if len(vals) == 0:
            continue
        method = order[i]
        df_method = (
            df.loc[df["method_label"] == method, ["target", col]]
            .dropna()
            .copy()
            .reset_index(drop=True)
        )
        proteins = df_method["target"].astype(str).str.upper().to_numpy()
        jitter_x = x_positions[i] + rng.uniform(-0.18, 0.18, size=len(vals))
        jitter_y = vals + rng.uniform(-0.015, 0.015, size=len(vals))
        jitter_y = np.clip(jitter_y, 0.0, 1.0)
        low_mask = vals < label_thr
        jitter_x[low_mask] = x_positions[i]
        jitter_y[low_mask] = vals[low_mask]
        ax.scatter(
            jitter_x,
            jitter_y,
            s=34,
            c=[cmap(float(np.clip(v, 0.0, 1.0))) for v in vals],
            alpha=0.82,
            edgecolors="white",
            linewidths=0.55,
            zorder=3,
        )
        label_idx = sorted(
            [j for j, yi in enumerate(vals) if yi < label_thr],
            key=lambda j: vals[j],
        )
        last_y_by_side = {"right": -10.0, "left": -10.0}
        min_gap = 0.03
        for rank, j in enumerate(label_idx):
            side = "right" if rank % 2 == 0 else "left"
            ha = "left" if side == "right" else "right"
            x_text = (jitter_x[j] + 0.08) if side == "right" else (jitter_x[j] - 0.08)
            y_text = float(jitter_y[j])
            if (y_text - last_y_by_side[side]) < min_gap:
                y_text = last_y_by_side[side] + min_gap
            y_text = float(np.clip(y_text, 0.0, 1.0))
            last_y_by_side[side] = y_text
            ax.annotate(
                proteins[j],
                (jitter_x[j], jitter_y[j]),
                textcoords="data",
                xytext=(x_text, y_text),
                ha=ha,
                va="center",
                fontsize=7,
                alpha=0.9,
                zorder=5,
            )
        ax.scatter(
            [x_positions[i]],
            [float(np.median(vals))],
            s=50,
            marker="D",
            c=["#222222"],
            edgecolors="white",
            linewidths=0.7,
            zorder=4,
        )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(order, fontsize=9)
    ax.set_ylim(-0.03, 1.03)
    ax.set_ylabel(f"nEF{frac}", fontsize=11)
    ax.grid(axis="y", linestyle="--", alpha=0.35, linewidth=0.7)
    ax.axhline(0.5, color="#666666", linestyle=":", linewidth=1.0, alpha=0.7)
    ax.set_title(title, fontsize=12, pad=10)
    ax.text(
        0.03,
        0.98,
        panel_letter,
        transform=ax.transAxes,
        fontsize=18,
        fontweight="bold",
        va="top",
        ha="left",
    )


def load_ef_summary(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    required = {"target", "nEF10_active", "nEF10_low"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in {csv_path}: {sorted(missing)}")

    # Ensure method_label exists (article pipeline CSV has it; older CSVs may not).
    if "method_label" not in df.columns:
        if "method_id" in df.columns:
            id_to_label = dict(METHOD_ORDER)
            df = df.copy()
            df["method_label"] = df["method_id"].map(id_to_label)
        elif "metric" in df.columns:
            metric_to_label = {
                "boltz2_affinity_pred_value": "Boltz-2",
                "dynamicbind_affinity_bestpose": "DynamicBind",
                "gnina_cnn_affinity_bestpose": "GNINA 1.3",
                "gnina_cnn_affinity_min": "GNINA 1.3",
                "plapt_affinity": "PLAPT",
                "qvina_affinity_bestpose": "QVina2",
            }
            df = df.copy()
            df["method_label"] = df["metric"].map(metric_to_label)
        else:
            raise ValueError(
                f"{csv_path} needs method_label, method_id, or metric to map methods."
            )

    df = df.dropna(subset=["method_label"]).copy()
    return df


def plot_nef10_svg(ef: pd.DataFrame, output_path: Path, seed: int = RANDOM_SEED) -> Path:
    present_labels = set(ef["method_label"].astype(str))
    order = [label for _, label in METHOD_ORDER if label in present_labels]
    if not order:
        raise ValueError("None of the expected method labels found in the EF summary.")

    cmap = positive_sequential_cmap()
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.6), sharey=True)
    panel_plot(
        axes[0],
        ef,
        "nEF10_active",
        f"nEF{FRAC} (actives: Ki < 1000 nM)",
        "A",
        order,
        cmap,
        seed=seed,
        frac=FRAC,
    )
    panel_plot(
        axes[1],
        ef,
        "nEF10_low",
        f"nEF{FRAC},low (low-affinity: Ki ≥ 1000 nM)",
        "B",
        order,
        cmap,
        seed=seed + 1,
        frac=FRAC,
    )
    plt.tight_layout()

    output_path = output_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    return output_path


def default_paths() -> tuple[Path, Path]:
    repo_root = Path(__file__).resolve().parent.parent
    input_csv = repo_root / "analysis" / "tables" / "enrichment" / "ef_summary_all_proteins.csv"
    output_svg = repo_root / "nEF10.svg"
    return input_csv, output_svg


def main() -> None:
    default_in, default_out = default_paths()
    parser = argparse.ArgumentParser(
        description="Reproduce nEF10.png as an SVG (active + inactive violin panels)."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=default_in,
        help=f"Path to ef_summary_all_proteins.csv (default: {default_in})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=default_out,
        help=f"Output SVG path (default: {default_out})",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=RANDOM_SEED,
        help=f"Jitter random seed (default: {RANDOM_SEED})",
    )
    args = parser.parse_args()

    if not args.input.exists():
        raise FileNotFoundError(f"Input CSV not found: {args.input}")

    ef = load_ef_summary(args.input)
    out = plot_nef10_svg(ef, args.output, seed=args.seed)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
