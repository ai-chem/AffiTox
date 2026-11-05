"""Runtime configuration parsed from CLI / environment."""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass, field
from pathlib import Path

from .constants import METHOD_LABELS, PRIMARY_METRICS, TARGET_LABELS


def _split_list(value: str) -> list[str]:
    return [x.strip() for x in value.replace(",", " ").split() if x.strip()]


def _parse_method_colors(raw: str) -> dict[str, str]:
    out: dict[str, str] = {}
    if not raw:
        return out
    for part in raw.split(","):
        part = part.strip()
        if not part or ":" not in part:
            continue
        mid, color = part.split(":", 1)
        out[mid.strip()] = color.strip()
    return out


@dataclass
class AnalysisConfig:
    analysis_root: Path
    merged_dir: Path
    posebusters_dir: Path
    posebusters_boltz2_dir: Path
    dynamicbind_posebusters_label: str
    targets: list[str]
    methods: list[str]
    exp_col: str
    color_low: str
    color_mid: str
    color_high: str
    method_colors: dict[str, str] = field(default_factory=dict)
    figure_dpi: int = 700
    font_family: str = "DejaVu Sans"
    target_label_mode: str = "pdb"
    figure_variant: str = ""
    n_bootstrap: int = 10_000
    n_permutation: int = 10_000
    random_seed: int = 42
    run_correlations: bool = True
    run_enrichment: bool = True
    run_posebusters: bool = True
    run_inferential: bool = True
    run_figures: bool = True

    @property
    def tables_dir(self) -> Path:
        return self.analysis_root / "tables"

    @property
    def figures_dir(self) -> Path:
        base = self.analysis_root / "figures"
        return base / self.figure_variant if self.figure_variant else base

    def method_label(self, method_id: str) -> str:
        return METHOD_LABELS.get(method_id, method_id)

    def target_label(self, target: str) -> str:
        target_id = target.lower()
        if self.target_label_mode == "protein":
            return TARGET_LABELS.get(target_id, target.upper())
        return target.upper()

    def primary_metric(self, method_id: str) -> str:
        if method_id not in PRIMARY_METRICS:
            raise KeyError(f"Unknown method: {method_id}")
        return PRIMARY_METRICS[method_id]

    def ensure_dirs(self) -> None:
        for sub in (
            "correlations/per_target",
            "correlations",
            "enrichment/per_target",
            "enrichment",
            "posebusters",
            "inferential",
            "correlations/heatmaps",
            "correlations/scatter",
            "enrichment/figures",
            "posebusters/heatmaps",
            "posebusters/summary",
        ):
            (self.tables_dir / sub).mkdir(parents=True, exist_ok=True)
            (self.figures_dir / sub.replace("per_target", "")).mkdir(parents=True, exist_ok=True)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="AffiTox article analysis pipeline")
    p.add_argument("--analysis-root", type=Path, required=True)
    p.add_argument("--merged-dir", type=Path, required=True)
    p.add_argument("--posebusters-dir", type=Path, required=True)
    p.add_argument("--posebusters-boltz2-dir", type=Path, required=True)
    p.add_argument("--dynamicbind-posebusters-label", default="dynamicbind")
    p.add_argument("--targets", default=os.environ.get("TARGETS", ""))
    p.add_argument("--methods", default=os.environ.get("METHODS", ""))
    p.add_argument("--exp-col", default=os.environ.get("EXP_COL", "pValue"))
    p.add_argument("--color-low", default=os.environ.get("COLOR_LOW", "#3558C5"))
    p.add_argument("--color-mid", default=os.environ.get("COLOR_MID", "#FFF9EE"))
    p.add_argument("--color-high", default=os.environ.get("COLOR_HIGH", "#EF4938"))
    p.add_argument("--method-colors", default=os.environ.get("METHOD_COLORS", ""))
    p.add_argument("--figure-dpi", type=int, default=int(os.environ.get("FIGURE_DPI", "300")))
    p.add_argument("--font-family", default=os.environ.get("FONT_FAMILY", "DejaVu Sans"))
    p.add_argument(
        "--target-label-mode",
        choices=("pdb", "protein"),
        default=os.environ.get("TARGET_LABEL_MODE", "pdb"),
    )
    p.add_argument("--figure-variant", default=os.environ.get("FIGURE_VARIANT", ""))
    p.add_argument("--n-bootstrap", type=int, default=int(os.environ.get("N_BOOTSTRAP", "10000")))
    p.add_argument("--n-permutation", type=int, default=int(os.environ.get("N_PERMUTATION", "10000")))
    p.add_argument("--random-seed", type=int, default=int(os.environ.get("RANDOM_SEED", "42")))
    p.add_argument("--skip-correlations", action="store_true")
    p.add_argument("--skip-enrichment", action="store_true")
    p.add_argument("--skip-posebusters", action="store_true")
    p.add_argument("--skip-inferential", action="store_true")
    p.add_argument("--skip-figures", action="store_true")
    return p


def config_from_args(args: argparse.Namespace) -> AnalysisConfig:
    targets = _split_list(args.targets)
    methods = _split_list(args.methods)
    if not targets:
        raise ValueError("No targets specified (--targets or TARGETS env)")
    if not methods:
        raise ValueError("No methods specified (--methods or METHODS env)")

    return AnalysisConfig(
        analysis_root=args.analysis_root,
        merged_dir=args.merged_dir,
        posebusters_dir=args.posebusters_dir,
        posebusters_boltz2_dir=args.posebusters_boltz2_dir,
        dynamicbind_posebusters_label=args.dynamicbind_posebusters_label,
        targets=[t.lower() for t in targets],
        methods=[m.lower() for m in methods],
        exp_col=args.exp_col,
        color_low=args.color_low,
        color_mid=args.color_mid,
        color_high=args.color_high,
        method_colors=_parse_method_colors(args.method_colors),
        figure_dpi=args.figure_dpi,
        font_family=args.font_family,
        target_label_mode=args.target_label_mode,
        figure_variant=args.figure_variant.strip(),
        n_bootstrap=args.n_bootstrap,
        n_permutation=args.n_permutation,
        random_seed=args.random_seed,
        run_correlations=not args.skip_correlations,
        run_enrichment=not args.skip_enrichment,
        run_posebusters=not args.skip_posebusters,
        run_inferential=not args.skip_inferential,
        run_figures=not args.skip_figures,
    )
