# GitHub to Zenodo software release (do this after CI is green)

Do **not** publish GitHub `v1.0.0` until Actions is green, including the analysis
smoke test. `analysis/config/defaults.sh` honors `PYTHON` (no cluster interpreter
hard-coded). `DEFAULT_METHOD_COLORS` lives in `src/analysis/constants.py`.

## 1. Release commit (on a clean, reviewed diff)

```bash
git status
# stage only reproducibility-release files, not unrelated WIP
git commit -m "Prepare AffiTox v1.0.0 reproducibility release"
git push origin main
```

Wait for GitHub Actions to pass, including:

```text
python -c "from analysis.plots.style import apply_style"
python src/analysis/run_pipeline.py --help
N_BOOTSTRAP=10 N_PERMUTATION=10 TARGETS=1g5m  analysis smoke
```

## 2. Enable Zenodo GitHub integration **before** publishing the release

https://help.zenodo.org/docs/github/enable-repository/

Profile → GitHub → Sync now → enable `ai-chem/AffiTox`.

Keep `CITATION.cff`; do not add `.zenodo.json`.

## 3. Draft then publish GitHub release

Tag: `v1.0.0` on `main`. Title: `AffiTox v1.0.0 J Cheminform reproducibility release`.

Suggested notes:

```markdown
## AffiTox v1.0.0

Reproducibility release accompanying the AffiTox *Journal of Cheminformatics* manuscript.

### Included

- curated 16-target affinity panel;
- manuscript-level merged prediction tables;
- scoring, ranking, screening, and PoseBusters analysis;
- publication figures in SVG and PNG;
- scientific invariant tests;
- installation and reproduction instructions.

### Raw docking outputs

Do **not** cite `10.5281/zenodo.20825057`–`067` as the current data. Those are
historical per-protein zips (including 7awe/2z5x/3mjg). The v1.0.0 data are the
information-block records in `release/zenodo/` (01 raw, 02 prepared+weights,
03 docking shards, 04 merged, 05 PoseBusters, 06 metrics/PNG). Record 03 is
split across sibling Zenodo depositions because of the 50 GiB/record limit.

### Scope

The repository supports analysis-level reproduction from deposited tables.
Full docking additionally requires separately installed third-party methods,
model checkpoints, and compatible computational hardware.
```

## 4. Check the auto-created Zenodo Software record

Title AffiTox, version 1.0.0, type Software, MIT, ORCID, GitHub URL, community affitox.

Use the **version-specific software DOI** in the paper, not a raw-data part DOI.

## 5. Umbrella data record

Published: [10.5281/zenodo.22007993](https://doi.org/10.5281/zenodo.22007993)
(`AffiTox v1.0.0: data release manifest`), Has-part 01 / 02 / 03a–e / 04 / 05 / 06.
See `docs/zenodo/README.md`. The two small files on that record are the August 2026
stub (old TSV); the landing-page description and related identifiers are the
authoritative index.
