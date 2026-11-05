# Target-Resolved Multi-Axis Benchmarking of Docking and Affinity Methods on Safety-Relevant Proteins

*Journal of Cheminformatics* (in review)

**Abstract.** Evaluating computational binding methods for
safety-relevant proteins requires understanding how
well they rank affinities, identify strong and weak
binders, and predict plausible binding poses.
We introduce AffiTox, a retrospective benchmark
of five docking and affinity-prediction methods
comprising 10711 curated ligand-target records
with experimental $K_i$ measurements for 16 proteins.

Boltz-2 had the highest within-target Pearson
correlation with experimental $\mathrm{p}K_i$
($r=0.55$), although its advantage varied across
targets and evaluation tasks.
High active enrichment occurred even when affinity
correlations were weak, particularly in sets dominated
by active compounds. Low-affinity enrichment revealed
ranking limitations that active enrichment alone
did not capture. Pose validity provided complementary
information: methods with high pose-validity rates
could still produce scores that correlated poorly
with experimental affinity. These findings support
target-specific evaluation of the prediction task
for which a method will be used.

Scientific Contribution:
AffiTox applies screening power to assess how well
methods recover experimentally defined low-affinity
compounds at the low-affinity end of their predicted
rankings. This assessment complements within-target
scoring power, ranking power, and pose validity.
The resource provides data and protocols for
evaluating and developing binding-prediction methods
for safety-relevant proteins.


<p align="center">
  <img src="docs/assets/affitox_panel.png"
       alt="Overview of the AffiTox benchmark dataset, evaluated methods, and analysis axes"
       width="100%">
</p>

## Dataset

Full-panel data: **[10.5281/zenodo.22007993](https://doi.org/10.5281/zenodo.22007993)**
(`AffiTox v1.0.0: data release manifest`). That is the DOI to put in the paper.
Community page (no DOI): [zenodo.org/communities/affitox](https://zenodo.org/communities/affitox/).

What to download, what is inside each zip, and which folder in this clone to
copy into: **[`docs/zenodo/README.md`](docs/zenodo/README.md)**.

This repository is code, 16-target inputs, `processed/` receptors, and a
16×5 `example/` slice. It does **not** contain `results/` (~187 GiB docking
dumps) or full-panel `analysis/tables` and `analysis/figures` (gitignored).
Those are Zenodo records 03 and 04–06.

## Repository layout

```text
AffiTox/
├── config/                 paths, target list, method settings
├── input/                  proteins and curated Ki tables
├── processed/              prepared receptors, ligands, docking boxes
├── results/                raw docking and affinity outputs (not in git)
├── pipeline/               one wrapper per stage
├── src/docking_benchmark2/ prepare and dock adapters
├── src/analysis/           merge, PoseBusters, statistics, figures
├── analysis/tables/        merged scores and manuscript statistics
├── analysis/figures/       manuscript figures
├── docs/                   install, pipeline, methods, data
├── tests/                  metric invariant tests
└── run_pipeline.sh         stage driver
```

Read the folders in the order you actually run the work.

1. **[`config/`](config/project.env.sh)**  
   Sets the repository root (`TOXAFFINITY_ROOT`), the 16 PDB codes, and default YAML files. Source this before any stage. Method flags live in [`config/methods_config.yaml`](config/methods_config.yaml); protein preparation in [`config/protein_settings_keep_cofactors_v2.yaml`](config/protein_settings_keep_cofactors_v2.yaml). Details: [`docs/PIPELINE.md`](docs/PIPELINE.md).

2. **[`input/`](input/README.md)**  
   Starting data. Structures in `input/proteins/`. Curated $K_i$ tables in `input/ligands_nodubl/`. The pipeline reads these; it does not recurate ChEMBL. Column meanings: [`docs/DATA_DICTIONARY.md`](docs/DATA_DICTIONARY.md).

3. **`processed/`**  
   Output of **prepare**. Shared PDBQT receptors, ligand files, search boxes, and protein sequences for PLAPT. Every docking method is supposed to see this same prepared input. How it is built: [`pipeline/prepare/README.md`](pipeline/prepare/README.md).

4. **`results/`**  
   Output of **dock**. One tree per target and method (`results/<pdb>/docking/<method>/`). Git ignores this directory because it is large. Download the Zenodo zips listed under **What to download** if you want to re-run merge or PoseBusters without docking. What each engine writes: [`pipeline/dock/README.md`](pipeline/dock/README.md).

5. **`analysis/excluding_2z5x_3mjg/tables/`**  
   Output of **merge** and **posebusters**, then of **analysis**. Full-panel CSVs are gitignored; unpack Zenodo **04** / **05** / **06** ([`docs/zenodo/README.md`](docs/zenodo/README.md)). The 16×5 slice is in [`example/`](example/README.md).

6. **`analysis/excluding_2z5x_3mjg/figures/`**  
   PNG from record **06**. Which script makes which figure: [`docs/FIGURES.md`](docs/FIGURES.md).

7. **Code that moves data along this path**  
   [`run_pipeline.sh`](run_pipeline.sh) is the public switch. Wrappers: [`pipeline/README.md`](pipeline/README.md). Prepare/dock Python package: `src/docking_benchmark2/` (CLI names `affitox-pipeline` and `toxdock-pipeline`). Manuscript statistics: `src/analysis/`. Tests: [`tests/`](tests/test_scientific_invariants.py).

## Installation

Linux x86_64. Windows is not supported.

`environment.yml` is the version list for the shared chemistry stack (env name `docking`). It does **not** install QVina2, GNINA, DynamicBind, Boltz-2, PLAPT, or PoseBusters. Those are separate programs. After conda, `run_pipeline.sh install` only does `pip install -e .` so Python can import `docking_benchmark2`.

Pins below were checked against a live conda env named `docking` (Python 3.10.18, NumPy 1.26.4, pandas 2.3.3, SciPy 1.15.3, matplotlib 3.10.8, seaborn 0.13.2, meeko 0.7.1, Biopython 1.86, OpenBabel 3.1.1, AutoDock Vina 1.2.6, RDKit conda package 2025.09.3). `conda env create` itself was not re-run in this documentation pass.

### 1. Clone and conda stack (declared versions)

```bash
git clone https://github.com/ai-chem/AffiTox.git   # public code, inputs, and analysis tables
cd AffiTox

conda env create -f environment.yml   # env name: docking; includes CUDA toolkit + mgltools pins
# If that solver fails (no NVIDIA, bioconda mgltools, compilers), use the CPU file instead:
#   conda env create -f environment.cpu.yml
conda activate docking                 # Python 3.10.18 plus RDKit, OpenBabel, Vina, …

export TOXAFFINITY_ROOT="$(pwd)"      # path contract used by every wrapper
export PYTHON="$(command -v python)"  # use this env’s interpreter, not a random python3

bash run_pipeline.sh install          # pip install -e ".[analysis]" into the active env
```

Exact conda/pip versions are in [`environment.yml`](environment.yml): `python=3.10.18`, `numpy=1.26.4`, `pandas=2.3.3`, `rdkit=2025.09.3`, `openbabel=3.1.1`, `vina=1.2.6`, `mgltools=1.5.7`, `cuda-toolkit=12.6.2`, and pip pins `meeko==0.7.1`, `scipy==1.15.3`, `seaborn==0.13.2`, `biopython==1.86`, `pyyaml==6.0.3`, `gemmi==0.7.3`. [`environment.cpu.yml`](environment.cpu.yml) is the same Python/RDKit/OpenBabel/Vina/pip pins without NVIDIA toolkit, mgltools, or compilers. Prepare uses Meeko, not mgltools. GPU methods (Boltz-2, DynamicBind, PLAPT) have their own environments. Package metadata: [`setup.py`](setup.py). More context: [`docs/INSTALL.md`](docs/INSTALL.md).

### 2. If `conda env create -f environment.yml` fails or CUDA is unwanted

Use [`environment.cpu.yml`](environment.cpu.yml) (same env name `docking`). That file was **not** executed in this documentation pass; it is the portable subset of the live `docking` pins.

```bash
conda env create -f environment.cpu.yml
conda activate docking
export PYTHON="$(command -v python)"
bash run_pipeline.sh install
```

A further reduced line (also **not** re-run here) if you already have RDKit and OpenBabel:

```bash
conda create -n docking python=3.10.18 rdkit=2025.09.3 openbabel=3.1.1 vina=1.2.6 -c conda-forge -c bioconda
conda activate docking
python -m pip install biopython==1.86 meeko==0.7.1 scipy==1.15.3 seaborn==0.13.2 pyyaml==6.0.3 gemmi==0.7.3
export PYTHON="$(command -v python)"
bash run_pipeline.sh install
```

Prefer `environment.yml` only when you want the file-for-file match with the HPC env, including CUDA toolkit pins.

### 3. Docking engines (required only to rerun prepare/dock)

Install each tool in its own prefix, then point [`config/methods_config.yaml`](config/methods_config.yaml) or the environment variables `PLAPT_PATH`, `DYNAMICBIND_PATH`, `GNINA_BIN`, `BOLTZ_ROOT` at **your** copies. A path template: [`config/methods_config.hpc.example.yaml`](config/methods_config.hpc.example.yaml). Recovered CLI pins: [`docs/SOFTWARE_REGISTRY.md`](docs/SOFTWARE_REGISTRY.md).

| Program | What AffiTox calls | Notes |
|---------|--------------------|-------|
| QVina2 | binary `qvina02` | Not provided by `vina=1.2.6` in `environment.yml`. That pin is AutoDock Vina. |
| GNINA 1.3 | binary `gnina` | Published runs used `gnina v1.3 master:97fa6bc+` with `--no_gpu`. |
| DynamicBind | checkout with `run_single_protein_inference.py` | GPU. Separate conda env in typical setups. |
| Boltz-2 | `boltz predict` outside `run_pipeline.sh dock` | GPU. Set `BOLTZ_RESULTS_DIR` before merge. Package seen in a Boltz env: `boltz==2.2.0`. |
| PLAPT | WELP-PLAPT (`plapt.py` + ONNX weights) | GPU or CPU. No 3D poses. |
| PoseBusters | CLI `bust` on `PATH` | Needed for the posebusters stage. A current env reports `bust 0.4.6`; deposited CSVs do not store `bust --version`. |

### 4. Check that the Python package imported

```bash
python -c "import docking_benchmark2; print('ok')"   # prepare/dock import path
PYTHONPATH=src pytest tests -q                       # Ki to pKi and nEF invariants
```

## What to download

Index DOI (cite this): [10.5281/zenodo.22007993](https://doi.org/10.5281/zenodo.22007993).
Unpack paths, zip contents, and `curl` examples: [`docs/zenodo/README.md`](docs/zenodo/README.md).

| Stage | Already in the clone | Download from Zenodo |
|-------|----------------------|----------------------|
| tests, 16×5 example | `example/` | nothing |
| `run_pipeline.sh analysis` on the **full panel** | not in git | **04** + **06** → `analysis/excluding_2z5x_3mjg/{tables,figures}/` |
| `prepare` | `input/proteins/`, `input/ligands_nodubl/` | **01** only if you want the deposited ChEMBL/PDB snapshot |
| `dock` | `processed/` | **02** for pinned weights; no docking zips. Install engines (table above) |
| `merge` / re-PoseBusters | ligand CSVs, `processed/` | **03a–e** `results_<pdb>.zip` unzipped at the repo root → `results/<pdb>/` |
| PoseBusters CSVs without re-running `bust` | `example/05_posebusters/` (16×5) | **05** → `analysis/excluding_2z5x_3mjg/tables/posebuster/` |

```bash
# full-panel merged tables (record 04, ~5 MiB)
curl -L -o 04_merged.zip 'https://zenodo.org/records/23979579/files/04_merged.zip?download=1'
unzip 04_merged.zip
mkdir -p analysis/excluding_2z5x_3mjg/tables
cp 04_merged/merged_ligands_docking_*.csv analysis/excluding_2z5x_3mjg/tables/
```

Do not cite `10.5281/zenodo.20825057`–`067`.

## Pipeline

AffiTox takes protein structures and experimental $K_i$ tables and produces, for five methods, per-ligand scores (and poses where the method has them), then the manuscript tables and figures. Scoring asks how well a score tracks $pK_i$. Ranking asks whether the order of compounds is preserved. Screening asks whether strong binders (and, separately, weak binders) rise to the top of a sorted list. Pose validity asks whether a predicted pose passes PoseBusters checks. Those four questions are different; the pipeline keeps them in separate tables.

Driver:

```bash
source config/project.env.sh     # TOXAFFINITY_ROOT, TARGETS, YAML paths
bash run_pipeline.sh --help      # stages listed below; default with no args is analysis
```

Settings for each stage are **not** copied into this page. Follow the link under the stage.

### install

Puts this repository on `PYTHONPATH` via an editable pip install. It does not download GNINA, Boltz, or Zenodo archives. Needs: the clone only.

```bash
bash run_pipeline.sh install     # same as pipeline/install/setup_environment.sh
```

Details: [`pipeline/install/README.md`](pipeline/install/README.md), [`docs/INSTALL.md`](docs/INSTALL.md).

### prepare

Builds one prepared receptor, ligand set, and docking box per target so QVina2, GNINA, DynamicBind, and PLAPT do not each clean the protein differently. Needs: `input/` and [`config/interaction_protein_ligand_16target.json`](config/interaction_protein_ligand_16target.json) from the clone. Prepared outputs are already in `processed/` if you skip this stage.

```bash
bash run_pipeline.sh prepare     # writes processed/
```

The wrapper calls `python -m docking_benchmark2.cli.run_benchmark --stage preparation`. Pairing of PDB codes to ligand CSV names is [`config/interaction_protein_ligand_16target.json`](config/interaction_protein_ligand_16target.json), set as `interaction_config_file` in [`config/toxdock_config.yaml`](config/toxdock_config.yaml). Pass `--interaction-config` to use another JSON with the same schema (`src/docking_benchmark2/utils/settings.py`).

Details: [`pipeline/prepare/README.md`](pipeline/prepare/README.md), [`config/protein_settings_keep_cofactors_v2.yaml`](config/protein_settings_keep_cofactors_v2.yaml).

### dock

Runs the methods listed in [`config/toxdock_config.yaml`](config/toxdock_config.yaml): QVina2, GNINA, PLAPT, DynamicBind. Writes logs, poses, and extracted metric CSVs under `results/`. Needs: `processed/` from the clone, plus the engine binaries. No Zenodo download. If a listed engine is missing, this stage stops and prints what to install. Boltz-2 is not started here.

```bash
bash run_pipeline.sh dock        # does not run Boltz-2
```

**Boltz-2 is not started by this stage.** Run `boltz predict` in a Boltz workspace (sampling used for the paper tag: 80 structure steps, 10 diffusion samples, affinity 80 / 1, seed 42). Then:

```bash
export BOLTZ_RESULTS_DIR=/path/to/boltz/results   # merge and Boltz PoseBusters read this
```

To dock only some engines, pass the CLI flag (the `METHODS` env var is for analysis, not this wrapper):

```bash
bash pipeline/dock/run_docking.sh --methods gnina qvina
```

Details: [`pipeline/dock/README.md`](pipeline/dock/README.md), [`config/methods_config.yaml`](config/methods_config.yaml), [`docs/SOFTWARE_REGISTRY.md`](docs/SOFTWARE_REGISTRY.md).

### merge

Joins each target’s ligand table with per-method scores and adds `pValue` ($pK_i = 9 - \log_{10}(K_i)$ with $K_i$ in nM). Needs: `results/<pdb>/` from your own dock run, or unpacked `results_<pdb>.zip` from Zenodo. The published merged CSVs are already in `analysis/tables/`; skip this stage unless you are rebuilding them.

```bash
bash run_pipeline.sh merge       # writes analysis/tables/merged_ligands_docking_<pdb>.csv
```

Details: [`pipeline/postprocess/README.md`](pipeline/postprocess/README.md). Scripts: `src/analysis/merge_ligands_docking_from_dir.py`, `src/analysis/add_pvalue_column.py`.

### posebusters

Runs PoseBusters (`bust`) on predicted poses for QVina2, GNINA, and DynamicBind. PLAPT has no poses and is skipped. Boltz-2 is included only when `BOLTZ_RESULTS_DIR` is set. Needs: poses under `results/` and `bust` on `PATH` to recompute. If `bust` or `results/` is missing, the wrapper keeps the deposited CSVs in `analysis/tables/posebuster/` and exits 0.

```bash
bash run_pipeline.sh posebusters
```

Details: [`pipeline/postprocess/README.md`](pipeline/postprocess/README.md), [`docs/METHODS.md`](docs/METHODS.md).

### analysis

Computes scoring, ranking, screening, pose-validity summaries, and figures from the merged tables. For the **full panel**, unpack Zenodo records **04** and **06** into `analysis/excluding_2z5x_3mjg/` first ([`docs/zenodo/README.md`](docs/zenodo/README.md)). The 16×5 slice is already in `example/`. Docking is not required.

```bash
bash run_pipeline.sh analysis    # default if you call run_pipeline.sh with no arguments
# same entry point:
bash pipeline/postprocess/run_article_analysis.sh
```

Defaults (`analysis/config/defaults.sh`): all 16 targets, five methods, `N_BOOTSTRAP=10000`, `N_PERMUTATION=10000`, `RANDOM_SEED=42`. That resampling is slow. A short check:

```bash
N_BOOTSTRAP=10 N_PERMUTATION=10 TARGETS=1g5m \
  bash pipeline/postprocess/run_article_analysis.sh
```

Full SI profile (forces 10k / 10k and writes a log):

```bash
bash pipeline/postprocess/run_si_full.sh
```

Details: [`analysis/README.md`](analysis/README.md), [`docs/METHODS.md`](docs/METHODS.md), [`docs/FIGURES.md`](docs/FIGURES.md).

### all

```bash
bash run_pipeline.sh all         # prepare, dock, merge, posebusters, analysis
```

This skips `install` and still skips Boltz-2 inference. Needs: engine binaries for dock. Prepare uses the tracked 16-target JSON. PoseBusters without `bust` keeps deposited CSVs. `all` will still stop at dock if QVina2, GNINA, PLAPT, or DynamicBind are not installed; use `bash run_pipeline.sh analysis` for manuscript tables.

Licence: [MIT](LICENSE) for this repository’s source. BindingDB, ChEMBL, PDB records and third-party weights follow their own terms ([NOTICE](NOTICE)).

## How to cite

*Journal of Cheminformatics* (BMC, Vancouver). Cite the **umbrella data DOI**
once; do not list all ten part records in the article reference list.

**Availability of data and materials.** See the paste-ready paragraph in
[`docs/DATA_AVAILABILITY.md`](docs/DATA_AVAILABILITY.md).

Dataset:

Konovalova OA, Orlova A, Telepov A, Khrabrov K, Karpushkina I, Shestun P, Kadurin A, Vinogradov V, Tsypin A, Dmitrenko A. AffiTox v1.0.0: data release manifest [dataset]. Zenodo; 2026. https://doi.org/10.5281/zenodo.22007993.

Code:

Konovalova OA, Orlova A, Telepov A, Khrabrov K, Karpushkina I, Shestun P, Kadurin A, Vinogradov V, Tsypin A, Dmitrenko A. AffiTox. GitHub; 2026. https://github.com/ai-chem/AffiTox.

Machine-readable: [`CITATION.cff`](CITATION.cff).
