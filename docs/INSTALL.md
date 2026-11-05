# AffiTox installation

Linux x86_64. Windows is not supported.

There is no Docker image, lock file, or `requirements.txt`. Versions for the shared chemistry environment are [`environment.yml`](../environment.yml) (env name `docking`). The Python package is [`setup.py`](../setup.py) (import name `docking_benchmark2`).

`environment.yml` does **not** install QVina2 (`qvina02`), GNINA, DynamicBind, Boltz-2, PLAPT, or PoseBusters. [`environment.cpu.yml`](../environment.cpu.yml) is the same chemistry stack without `cuda-toolkit` and `mgltools` (prepare uses Meeko).

The first-page walkthrough is in the repository [README](../README.md). This page is the longer install note.

## What the `docking` env is for

Prepare, ligand PDBQT, boxes, merge, and analysis all use this stack: Python 3.10, RDKit, OpenBabel, Meeko, pandas, SciPy, matplotlib. AutoDock Vina 1.2.6 is in the YAML because it ships with that conda recipe; AffiTox docking still calls the **QVina2** binary `qvina02`, which you must provide separately.

Checked against a live env named `docking` (not a fresh `conda env create` in this pass):

| Package | `environment.yml` | Seen in live `docking` env |
|---------|-------------------|----------------------------|
| python | 3.10.18 | 3.10.18 |
| numpy | 1.26.4 | 1.26.4 |
| pandas | 2.3.3 | 2.3.3 |
| scipy (pip) | 1.15.3 | 1.15.3 |
| matplotlib-base | 3.10.8 | 3.10.8 |
| seaborn (pip) | 0.13.2 | 0.13.2 |
| meeko (pip) | 0.7.1 | 0.7.1 |
| biopython (pip) | 1.86 | 1.86 |
| rdkit | 2025.09.3 | conda-meta `rdkit-2025.09.3` |
| openbabel | 3.1.1 | 3.1.1 |
| vina | 1.2.6 | binary present |
| mgltools | 1.5.7 | listed in YAML |
| cuda-toolkit | 12.6.2 | listed in YAML |

## Create the env

```bash
git clone https://github.com/ai-chem/AffiTox.git
cd AffiTox
conda env create -f environment.yml   # not re-run in this documentation pass; CUDA + mgltools pins
# fallback without NVIDIA / mgltools:
# conda env create -f environment.cpu.yml
conda activate docking
export TOXAFFINITY_ROOT="$(pwd)"
export PYTHON="$(command -v python)"
bash run_pipeline.sh install          # pip install -e ".[analysis]" (fallback: .[rdkit] or -e .)
```

`TOXAFFINITY_ROOT` is required by the shell wrappers. Do not rename the variable.

## Method binaries

Copy [`config/methods_config.hpc.example.yaml`](../config/methods_config.hpc.example.yaml), replace the placeholders, and `export METHODS_CONFIG` to that file. Or set `PLAPT_PATH`, `DYNAMICBIND_PATH`, `GNINA_BIN`, `BOLTZ_ROOT`.

CLI recovered from deposited logs: [`SOFTWARE_REGISTRY.md`](SOFTWARE_REGISTRY.md).

PoseBusters: put `bust` on `PATH` before `bash run_pipeline.sh posebusters`.

Boltz-2: not launched by `run_pipeline.sh dock`. After `boltz predict`, export `BOLTZ_RESULTS_DIR`.

## Tests

```bash
PYTHONPATH=src pytest tests -q
```

## Data

Dataset index: [10.5281/zenodo.22007993](https://doi.org/10.5281/zenodo.22007993).
What is inside each zip and where to copy it in this clone:
[`docs/zenodo/README.md`](zenodo/README.md).
Do not use the restricted historical DOIs `10.5281/zenodo.20825057`–`067`.

Checksums: [`docs/zenodo/affitox_data_manifest.tsv`](zenodo/affitox_data_manifest.tsv).
