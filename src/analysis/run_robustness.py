#!/usr/bin/env python3
"""Orchestrate AffiTox robustness analyses (Parts A–D). Does not edit the manuscript."""

from __future__ import annotations

import argparse
import hashlib
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from analysis.cross_axis_association import run_cross_axis
from analysis.robustness_common import (
    DEFAULT_N_BALANCED,
    DEFAULT_N_BOOT,
    DEFAULT_N_PERM,
    DEFAULT_SEED,
    MERGED_DIR,
    POSE_DIR,
    ROBUST_TABLES,
    TABLES,
    TARGETS,
    RobustnessConfig,
    file_md5,
    inventory_row,
)
from analysis.screening_balanced_sensitivity import run_balanced_screening
from analysis.temporal_validation import run_temporal
from analysis.training_exposure_extension import run_training_exposure


def build_inventory() -> Path:
    files = [
        (TABLES / "correlations" / "correlations_with_ci_perm.csv", "target×method correlation", ["target", "method_id"]),
        (TABLES / "enrichment" / "ef_summary_all_proteins.csv", "target×method screening", ["target", "method_id"]),
        (TABLES / "posebusters" / "pass_rates_by_target_method.csv", "target×method pose-row summary", ["target", "method_id"]),
        (TABLES / "overlap" / "overlap_fractions_by_method.csv", "method×quantity exposure", ["method_id", "overlap_dimension"]),
        (TABLES / "overlap" / "pair_membership_flags.csv", "ligand–target pair", ["pdb", "canonical_smiles"]),
        (TABLES / "overlap" / "affinity_pairs_uniprot_inchikey.csv", "UniProt×InChIKey pair", ["uniprot", "parent_inchikey"]),
        (TABLES / "overlap" / "pearson_paired_delta_dual_stratum.csv", "target×method Δr (legacy B=2000)", ["pdb", "method_id"]),
        (TABLES / "overlap" / "external" / "dynamicbind" / "d3_with_clash_info.csv", "DynamicBind processed complex", ["pdb"]),
        (TABLES / "overlap" / "external" / "gnina" / "crossdocked_v1.3_train_pdbs.txt", "GNINA training PDB line", None),
        (TABLES / "overlap" / "external" / "LP_PDBBind_ligand_features.csv", "LP-PDBBind crystal ligand", ["pdb"]),
    ]
    for t in TARGETS:
        files.append(
            (MERGED_DIR / f"merged_ligands_docking_{t}.csv", "ligand (merged panel)", ["canonical_smiles"])
        )
        for method in ("boltz2", "gnina", "qvina", "dynamicbind"):
            files.append(
                (
                    POSE_DIR / f"posebusters_results_{t}_{method}.csv",
                    "pose row (PoseBusters)",
                    None,
                )
            )
    rows = [inventory_row(path, unit=unit, key_cols=keys) for path, unit, keys in files]
    inv = pd.DataFrame(rows)
    path = ROBUST_TABLES / "input_inventory.csv"
    inv.to_csv(path, index=False)
    return path


def _hash_outputs() -> str:
    h = hashlib.md5()
    for path in sorted(ROBUST_TABLES.glob("*.csv")):
        h.update(path.name.encode())
        h.update(file_md5(path).encode())
    return h.hexdigest()


def write_validation_report(cfg: RobustnessConfig, elapsed: float, hash1: str, hash2: str) -> None:
    text = f"""# Robustness analysis validation report

## 1. What was computable

- Ligand-level PoseBusters panel joined to the merged AffiTox tables (16×4).
- Cross-axis Spearman associations with target-cluster bootstrap (B={cfg.n_boot})
  and restricted within-target permutation (B={cfg.n_perm}), including two-way
  residualization against C(target)+C(method).
- Rank nonredundancy and leave-one-target-out for the primary pair.
- Method-specific training-exposure evidence catalog (existing intersections +
  DynamicBind ECFP4 similarity proxy vs processed-table crystal ligands).
- Temporal sensitivity Δr with independent ligand bootstrap B={cfg.n_boot},
  explicit eligibility (n≥{cfg.min_n} and pKi range ≥{cfg.min_pki_range}),
  candidate post-2022 manifest, and a feasibility table.
- Balanced-class nEF10 sensitivity (B={cfg.n_balanced}); this changes the estimand.

## 2. What was not computable

See `NOT_COMPUTABLE.md`. Exact Boltz-2/PLAPT membership, GNINA ligand/pair types
membership, and training-set ECFP4/sequence identity for engines without SMILES
manifests were **not** replaced by guessed numbers.

## 3. Conclusions that are supported

- Pose validity and scoring/screening endpoints can be jointly tabulated and
  their **association** (or discordance) estimated on the 16×4 panel.
- Training-exposure evidence is **method-specific** and must not be compared as
  interchangeable membership rates.
- Post-cutoff AffiTox slices support **temporal sensitivity**, not prospective
  validation.

## 4. Conclusions that remain unsupported

- Statistical independence of docking power and scoring power.
- Exact training-set membership for Boltz-2 or PLAPT.
- That a retrospective year split is prospective validation.
- That balanced nEF10 is the confirmatory screening estimand.

## 5–6. Manuscript-ready wording

See `manuscript_review_text.md`. The main manuscript was **not** edited.

## 7. Reproduce

```bash
export PYTHON=/mnt/tank/scratch/okonovalova/miniconda3/envs/docking/bin/python
cd {PROJECT_ROOT}
$PYTHON src/analysis/run_robustness.py --seed {cfg.seed} --n-boot {cfg.n_boot} --n-perm {cfg.n_perm} --n-balanced {cfg.n_balanced}
```

Tests:

```bash
$PYTHON -m pytest tests/test_robustness.py tests/test_scientific_invariants.py -q
```

## Reproducibility

- First-run output MD5 over robustness CSVs: `{hash1}`
- Second-run output MD5: `{hash2}`
- Identical: `{hash1 == hash2}`
- Elapsed (both runs): {elapsed:.1f} s
- Seed: {cfg.seed}
"""
    (ROBUST_TABLES / "robustness_validation_report.md").write_text(text)


def write_manuscript_text() -> None:
    assoc_path = ROBUST_TABLES / "cross_axis_associations.csv"
    rank_path = ROBUST_TABLES / "cross_axis_rank_agreement.csv"
    delta_path = ROBUST_TABLES / "temporal_within_target_delta_r.csv"
    prim_txt = "see cross_axis_associations.csv (not yet written)"
    rank_txt = "see cross_axis_rank_agreement.csv"
    delta_txt = "see temporal_within_target_delta_r.csv"
    if assoc_path.exists():
        assoc = pd.read_csv(assoc_path)
        prim = assoc.loc[assoc["is_primary"] == 1].iloc[0]
        prim_txt = (
            f"raw Spearman ρ = {prim['spearman_raw']:.3f} "
            f"(95% target-cluster CI {prim['spearman_raw_ci_low']:.3f} to {prim['spearman_raw_ci_high']:.3f}); "
            f"two-way-adjusted ρ = {prim['spearman_two_way_adjusted']:.3f} "
            f"(CI {prim['spearman_adj_ci_low']:.3f} to {prim['spearman_adj_ci_high']:.3f}); "
            f"restricted-permutation p = {prim['p_perm_adjusted_restricted']:.4g}; "
            f"n = {int(prim['n_observations'])} target–method observations "
            f"({int(prim['n_targets'])} targets × {int(prim['n_methods'])} methods)."
        )
    if rank_path.exists():
        ranks = pd.read_csv(rank_path)
        same = ranks.loc[
            (ranks["block"] == "pair_summary")
            & (ranks["axis_a"] == "pearson_r")
            & (ranks["axis_b"] == "strict_pass_end_to_end")
        ]
        if len(same):
            rank_txt = (
                f"the same unique top method on Pearson r and strict end-to-end validity "
                f"in {100 * float(same['proportion_same_top_method'].iloc[0]):.0f}% of targets "
                f"(mean absolute rank difference {float(same['mean_of_mean_abs_rank_diff'].iloc[0]):.2f})."
            )
    if delta_path.exists():
        delta = pd.read_csv(delta_path)
        both = delta.loc[delta["both_eligible"] == True]  # noqa: E712
        n_dual = int(both["target"].nunique()) if len(both) else 0
        delta_txt = (
            f"{n_dual} targets were inferentially eligible in both publication-year strata "
            f"(n≥10 and pKi range ≥1.0). Within-target Δr = r_post − r_pre uses independent "
            f"ligand bootstrap B=10,000. These are temporal-sensitivity estimates, not "
            f"prospective validation. Do not compare a 16-target pre-cutoff median with a "
            f"smaller post-cutoff median."
        )
    text = f"""# Proposed manuscript / SI wording (author review only)

The main article and SI TeX files were **not** modified. If the numbers below
survive author checks, they can be pasted into SI §13/§20 and a short main-text
pointer.

## Main text (Discussion, pose vs affinity)

Pose validity and affinity estimation remain **distinct and potentially
nonredundant evaluation axes** on AffiTox. Across 16 targets and the four
pose-generating methods, the prespecified association between ligand-level
strict PoseBusters end-to-end pass rate and within-target Pearson *r* is
{prim_txt}
The endpoints are **partially associated but sometimes discordant**:
{rank_txt}
Methods can rank compounds by *pK~i~* with invalid poses, and geometrically
plausible poses can carry little affinity information. We do not interpret
these results as statistical independence.

Training-exposure evidence and temporal sensitivity remain in SI §20. A
retrospective post-cutoff split is **temporal sensitivity**, not prospective
validation.

## SI Methods (new robustness subsection)

Cross-axis analyses use one target×method observation on the 16×4 pose panel
(PLAPT excluded). Pose generation coverage is
*n*_evaluable / *n*_input; conditional strict pass is
*n*_pass-all / *n*_evaluable (missing poses excluded); end-to-end strict pass
is *n*_pass-all / *n*_input (missing poses unsuccessful). Associations are
Spearman correlations with 10,000 target-cluster bootstrap CIs. Two-way
adjustment residualizes each endpoint on C(target)+C(method); uncertainty
uses the same cluster bootstrap; permutation shuffles the pose endpoint among
methods within each target and repeats residualization (B=10,000). The
primary pair is end-to-end strict pass vs Pearson *r*; remaining pairs are a
Benjamini–Hochberg secondary family. Balanced-class nEF~10~ resampling
(B=2,000) is a sensitivity analysis that changes the estimand.

Temporal Δ*r* uses independent ligand bootstrap (B=10,000) on strata defined
by ChEMBL document year. A target–method–stratum cell is eligible only if
*n*≥10 and *pK~i~* range ≥1.0. Medians are not compared across unequal
eligible target sets. The Boltz-2 cutoff is not assigned to PLAPT, DynamicBind,
or GNINA.

## SI Results

Primary association: {prim_txt}

Rank nonredundancy: {rank_txt}

Temporal sensitivity: {delta_txt}

Exact numeric tables: `analysis/tables/robustness/cross_axis_associations.csv`,
`cross_axis_rank_agreement.csv`, and `temporal_within_target_delta_r.csv`.
Published SI Table S24 remains unchanged until the authors elect to update it.

## Explicit statement

A retrospective post-cutoff analysis cannot become a true prospective
validation merely through statistical reanalysis. True prospective validation
requires frozen models and analysis rules before new experimental labels are
obtained.
"""
    (ROBUST_TABLES / "manuscript_review_text.md").write_text(text)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p.add_argument("--n-boot", type=int, default=DEFAULT_N_BOOT)
    p.add_argument("--n-perm", type=int, default=DEFAULT_N_PERM)
    p.add_argument("--n-balanced", type=int, default=DEFAULT_N_BALANCED)
    p.add_argument("--skip-second-run", action="store_true")
    args = p.parse_args()
    cfg = RobustnessConfig(
        seed=args.seed,
        n_boot=args.n_boot,
        n_perm=args.n_perm,
        n_balanced=args.n_balanced,
    )
    cfg.ensure_dirs()
    t0 = time.time()
    print("[0] input inventory")
    inv = build_inventory()
    print(f"    -> {inv}")
    print("[1] Part A cross-axis")
    run_cross_axis(cfg)
    print("[2] Part B training exposure")
    run_training_exposure(cfg)
    print("[3] Part C temporal / external readiness")
    run_temporal(cfg)
    print("[4] Part D balanced nEF sensitivity")
    run_balanced_screening(cfg)
    write_manuscript_text()
    h1 = _hash_outputs()
    if args.skip_second_run:
        h2 = h1
    else:
        print("[5] second run for reproducibility")
        run_cross_axis(cfg)
        run_training_exposure(cfg)
        run_temporal(cfg)
        run_balanced_screening(cfg)
        h2 = _hash_outputs()
    write_manuscript_text()
    write_validation_report(cfg, time.time() - t0, h1, h2)
    print(f"done. identical_outputs={h1 == h2}  hash={h2}")
    print(f"tables: {ROBUST_TABLES}")


if __name__ == "__main__":
    main()
