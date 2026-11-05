# AffiTox v1.0.0 data — what to download and where it goes

**Cite this index in the paper:** [10.5281/zenodo.22007993](https://doi.org/10.5281/zenodo.22007993)
(`AffiTox v1.0.0: data release manifest`). Concept DOI:
[10.5281/zenodo.22007992](https://doi.org/10.5281/zenodo.22007992).
Community: https://zenodo.org/communities/affitox/

Code: https://github.com/OlgaKoont/AffiTox

Do **not** cite the restricted historical per-protein zips
`10.5281/zenodo.20825057`–`067` (old panel including 7awe/2z5x/3mjg).

All commands below are from the clone root (`export TOXAFFINITY_ROOT="$(pwd)"`).
Each zip already contains a folder prefix; unzip, then copy into the repo paths
in the right-hand column. File lists: [`affitox_data_manifest.tsv`](affitox_data_manifest.tsv).

Canonical analysis root after unpack: `analysis/excluding_2z5x_3mjg/`.
Panel PDB order: 1g5m, 2v5z, 3eyg, 3jy9, 3lxk, 11ue, 4ase, 4f65, 4tz4, 4zau,
5jkv, 5mo4, 6gqj, 6jok, 5lf3, 7kk3. Panel N (curated / merged rows): **10711**.

GPU re-runs of Boltz-2 / DynamicBind are not bit-identical. Record 03 is the
frozen snapshot; rebuild 04–06 from those files if you need the same numbers.

## Which record for which job

| You want | Download | Unzip then copy into |
|----------|----------|----------------------|
| Paper citation only | nothing (this DOI) | — |
| Tests / 16×5 example | nothing | already in `example/` |
| ChEMBL + PDB snapshot used in the paper | **01** | `input/ligands/`, `input/proteins/` |
| Prepared receptors, boxes, `ligand_id` maps, method weights | **02** | `processed/…`, `input/ligands_curated/`, weights as noted |
| Raw docking poses/logs (re-merge / re-PoseBusters) | **03a–e** | `results/<pdb>/` (prefix is already in the zip) |
| Merged score tables (`merged_ligands_docking_<pdb>.csv`) | **04** | `analysis/excluding_2z5x_3mjg/tables/` |
| Full-panel PoseBusters CSVs | **05** | `analysis/excluding_2z5x_3mjg/tables/posebuster/` |
| Correlation / nEF / inferential tables + PNG | **06** | `analysis/excluding_2z5x_3mjg/tables/` and `…/figures/` |

GitHub does **not** carry full-panel `analysis/tables/` or `analysis/figures/`
(gitignored). Use record 04 + 06 for manuscript tables and PNGs. The clone
already has `input/proteins/`, `input/ligands_nodubl/`, and `processed/` for
the 16 targets, so 01–02 are optional unless you want the exact deposited
snapshot or the pinned weights.

## Records, DOI, zip, what is inside

### Index (cite this)

- DOI: [10.5281/zenodo.22007993](https://doi.org/10.5281/zenodo.22007993)
- Inside: landing page + Has-part links to 01–06. Two small stub files from
  an earlier upload; the HTML description and related identifiers are the
  authoritative index.

### 01 — raw ChEMBL extracts and PDB structures (~36 MiB)

- DOI: [10.5281/zenodo.23213812](https://doi.org/10.5281/zenodo.23213812)
- File: `01_raw.zip` → folder `01_raw/`
- Inside: `proteins/<pdb>.pdb` (16 structures); `ligands/*.csv` (ChEMBL API
  snapshots, original names); `download_manifest.tsv`; NOTICE / README /
  SHA256SUMS.
- Copy:

```bash
curl -L -o 01_raw.zip 'https://zenodo.org/records/23213812/files/01_raw.zip?download=1'
unzip 01_raw.zip
mkdir -p input/proteins input/ligands
cp 01_raw/proteins/*.pdb input/proteins/
cp 01_raw/ligands/*.csv input/ligands/
```

Do not query live ChEMBL to reproduce the paper; use these files.

### 02 — prepared inputs, ligand_id maps, weights (~4.6 GiB)

- DOI: [10.5281/zenodo.23979597](https://doi.org/10.5281/zenodo.23979597)
  (replaces [10.5281/zenodo.23213816](https://doi.org/10.5281/zenodo.23213816))
- File: `02_prepared.zip` → folder `02_prepared/`
- Inside: `ligands/<pdb>_ligands.csv` (stable `ligand_id`); `id_maps/`;
  `proteins/<pdb>.pdb` and `.pdbqt` (Meeko); `boxes/<pdb>.json`;
  `plapt_sequences/<pdb>.txt`; `weights/` (Boltz-2 / DynamicBind / PLAPT
  checkpoints); `N_panel.txt`.
- Copy:

```bash
curl -L -o 02_prepared.zip 'https://zenodo.org/records/23979597/files/02_prepared.zip?download=1'
unzip 02_prepared.zip
mkdir -p input/ligands_curated processed/proteins processed/boxes \
         processed/id_maps processed/plapt_sequences
cp 02_prepared/ligands/*.csv input/ligands_curated/
cp 02_prepared/N_panel.txt input/ligands_curated/
cp 02_prepared/proteins/* processed/proteins/
cp 02_prepared/boxes/*.json processed/boxes/
cp 02_prepared/id_maps/*.csv processed/id_maps/
cp 02_prepared/plapt_sequences/*.txt processed/plapt_sequences/
# weights: see 02_prepared/weights/README.md (Boltz-2 → ~/.boltz/, etc.)
```

Join key for every later record is `(pdb_id, ligand_id)`.

### 03 — frozen docking dumps (sharded; Zenodo 50 GiB/record)

Unzip **from the clone root**. Each archive already starts with `results/<pdb>/`.

| Shard | DOI | Zip files | PDB |
|-------|-----|-----------|-----|
| 03a | [10.5281/zenodo.23981992](https://doi.org/10.5281/zenodo.23981992) (replaces [10.5281/zenodo.23224051](https://doi.org/10.5281/zenodo.23224051)) | `results_{3eyg,4f65,5lf3,2v5z}.zip` | 3eyg, 4f65, 5lf3, 2v5z |
| 03b | [10.5281/zenodo.23982808](https://doi.org/10.5281/zenodo.23982808) (replaces [10.5281/zenodo.23225395](https://doi.org/10.5281/zenodo.23225395)) | `results_{3jy9,3lxk}.zip` | 3jy9, 3lxk |
| 03c | [10.5281/zenodo.23983375](https://doi.org/10.5281/zenodo.23983375) (replaces [10.5281/zenodo.23225788](https://doi.org/10.5281/zenodo.23225788)) | `results_{5mo4,6gqj,4zau}.zip` | 5mo4, 6gqj, 4zau |
| 03d | [10.5281/zenodo.23984015](https://doi.org/10.5281/zenodo.23984015) (replaces [10.5281/zenodo.23226603](https://doi.org/10.5281/zenodo.23226603)) | `results_{7kk3,5jkv,4tz4}.zip` | 7kk3, 5jkv, 4tz4 |
| 03e | [10.5281/zenodo.23984506](https://doi.org/10.5281/zenodo.23984506) (replaces [10.5281/zenodo.23227176](https://doi.org/10.5281/zenodo.23227176)) | `results_{4ase,6jok,11ue,1g5m}.zip` | 4ase, 6jok, 11ue, 1g5m |

Inside each `results/<pdb>/docking/`: method trees (`qvina`, `gnina`, `plapt`,
`dynamicbind` / `dynamicbind_new`, `boltz2_*`) with poses, logs, and native
JSON/SDF.

```bash
# example: 1g5m only (~2.8 GiB), record 03e
curl -L -o results_1g5m.zip \
  'https://zenodo.org/records/23984506/files/results_1g5m.zip?download=1'
unzip results_1g5m.zip
ls results/1g5m/docking
```

You do not need all of record 03 unless you re-merge or re-run PoseBusters on
every target.

### 04 — merged ligand_id tables (~5.3 MiB)

- DOI: [10.5281/zenodo.23979579](https://doi.org/10.5281/zenodo.23979579)
  (replaces [10.5281/zenodo.23214418](https://doi.org/10.5281/zenodo.23214418))
- File: `04_merged.zip` → folder `04_merged/`
- Inside: one `merged_ligands_docking_<pdb>.csv` per target; `merge_audit.json`;
  `N_panel.txt` (10711).
- Copy:

```bash
curl -L -o 04_merged.zip 'https://zenodo.org/records/23979579/files/04_merged.zip?download=1'
unzip 04_merged.zip
mkdir -p analysis/excluding_2z5x_3mjg/tables
cp 04_merged/merged_ligands_docking_*.csv analysis/excluding_2z5x_3mjg/tables/
cp 04_merged/N_panel.txt analysis/excluding_2z5x_3mjg/tables/
```

### 05 — PoseBusters tables (~13 MiB)

- DOI: [10.5281/zenodo.23979721](https://doi.org/10.5281/zenodo.23979721)
  (replaces [10.5281/zenodo.23979584](https://doi.org/10.5281/zenodo.23979584))
- File: `05_posebusters.zip` → folder `05_posebusters/`
- Inside: `tables/posebusters_results_<pdb>_<method>.csv` for Boltz-2,
  DynamicBind (`*_dynamicbind.csv`, 21 checks including `internal_energy`),
  GNINA, QVina2 (no PLAPT). The superseded 20-check `*_dynamicbind.csv` files
  are not included.
- Copy:

```bash
curl -L -o 05_posebusters.zip 'https://zenodo.org/records/23979721/files/05_posebusters.zip?download=1'
unzip 05_posebusters.zip
mkdir -p analysis/excluding_2z5x_3mjg/tables/posebuster
cp 05_posebusters/tables/*.csv analysis/excluding_2z5x_3mjg/tables/posebuster/
```

### 06 — metrics tables and PNG figures (~138 MiB)

- DOI: [10.5281/zenodo.23979767](https://doi.org/10.5281/zenodo.23979767)
  (replaces [10.5281/zenodo.23979728](https://doi.org/10.5281/zenodo.23979728))
- File: `06_metrics_figures.zip` → folder `06_metrics_figures/`
- Inside: `tables/` (merged CSVs plus `correlations/`, `enrichment/`,
  `inferential/`, PoseBusters pass-count curves on **21** checks);
  `figures/by_protein/` (publication axis) and `figures/by_pdb/`.
- Copy:

```bash
curl -L -o 06_metrics_figures.zip \
  'https://zenodo.org/records/23979767/files/06_metrics_figures.zip?download=1'
unzip 06_metrics_figures.zip
mkdir -p analysis/excluding_2z5x_3mjg/tables analysis/excluding_2z5x_3mjg/figures
cp -a 06_metrics_figures/tables/. analysis/excluding_2z5x_3mjg/tables/
cp -a 06_metrics_figures/figures/. analysis/excluding_2z5x_3mjg/figures/
```

Publication heatmap example:
`analysis/excluding_2z5x_3mjg/figures/by_protein/correlations/heatmaps/pearson_heatmap_ultramarine_ivory_vermilion_square_cells.png`

## Journal of Cheminformatics (BMC) — how to cite

Cite **one** data DOI (the umbrella), not the ten part DOIs and not the
historical 20825057–067 records. Cite the GitHub repository for code until a
Zenodo *Software* DOI exists for tag `v1.0.0`.

**Availability of data and materials** (Declarations):

> The AffiTox v1.0.0 dataset is available from Zenodo
> (https://doi.org/10.5281/zenodo.22007993). It indexes six related records:
> raw ChEMBL and PDB snapshots, prepared inputs and method weights, frozen
> docking outputs (split across five deposits because of the 50 GiB/record
> limit), merged ligand tables, PoseBusters tables, and metrics/PNG figures.
> Source code is available at https://github.com/OlgaKoont/AffiTox (MIT).
> ChEMBL and PDB files remain under their original terms. GPU re-runs of
> Boltz-2 and DynamicBind are not bit-identical; the docking record is the
> frozen snapshot used in this study.

In-text: cite the dataset as a numbered reference, e.g. “Data are archived on
Zenodo [n].”

**Vancouver reference** (BMC / *J Cheminform*):

Konovalova OA, Orlova A, Telepov A, Khrabrov K, Karpushkina I, Shestun P, Kadurin A, Vinogradov V, Tsypin A, Dmitrenko A. AffiTox v1.0.0: data release manifest [dataset]. Zenodo; 2026. https://doi.org/10.5281/zenodo.22007993.

Code (until the software Zenodo archive exists):

Konovalova OA, Orlova A, Telepov A, Khrabrov K, Karpushkina I, Shestun P, Kadurin A, Vinogradov V, Tsypin A, Dmitrenko A. AffiTox. GitHub; 2026. https://github.com/OlgaKoont/AffiTox.

BibTeX is in [`CITATION.cff`](../../CITATION.cff) and
[`DATA_AVAILABILITY.md`](../DATA_AVAILABILITY.md).
