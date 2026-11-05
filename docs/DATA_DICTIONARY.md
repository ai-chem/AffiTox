# AffiTox data dictionary

Layer 1 is the ChEMBL API snapshot in `input/ligands/`. Frozen docked tables are
`input/ligands_nodubl/`. Public curated tables with `ligand_id` are
`input/ligands_curated/<pdb>_ligands.csv`. BindingDB dumps in `input/bindingdb/`
are provenance only (gitignored). Curation CLI:
`bash pipeline/prepare/run_curate.sh` (`--mode frozen-nodubl` keeps docked row
order; `--mode from-raw` applies the filters below).

## Curation order (ChEMBL snapshot → panel)

1. Keep experimentally measured \(K_i\) only.
2. Keep nM, or convert M / mM / µM / pM → nM; drop other units.
3. Deduplicate on `molecule_chembl_id` (else SMILES); median \(K_i\) if replicates disagree.
4. Drop RDKit MW \(> 500\) Da.
5. Require parseable SMILES.
6. Assign `ligand_id` (`ligand_0001`, …) in remaining row order **before docking**.
7. Convert to \(\mathrm{p}K_i = 9 - \log_{10}(K_i\,[\mathrm{nM}])\).
8. Label actives \(K_i < 1000\) nM and inactives \(K_i \ge 1000\) nM (screening only; inactives stay in the panel).

## `ligands_nodubl` columns

| Column | Meaning |
|--------|---------|
| `canonical_smiles` | RDKit canonical SMILES |
| `standard_type` | Endpoint label (`Ki`) |
| `standard_value` | \(K_i\) in nM |
| `pchembl_value` | p-scale potency when present |

`duplicates_report.json` records removed replicate IDs.

## Panel target map (from `src/analysis/merge_ligands_docking_from_dir.py`)

| PDB | Ligand CSV stem (`input/ligands_nodubl/`) |
|-----|-------------------------------------------|
| 1g5m | BCL2_Ki_WT_ChEMBL_252 |
| 2v5z | MAO-B_Ki_WT_ChEMBL_246 |
| 3eyg | JAK1_Ki_WT_ChEMBL_2255 |
| 3jy9 | JAK2_Ki_WT_ChEMBL_2027 |
| 3lxk | JAK3_Ki_WT_ChEMBL_786 |
| 11ue | PDGFRB_Ki_WT_ChEMBL_275 |
| 4ase | VEGFR2_Ki_WT_ChEMBL_875 |
| 4f65 | FGFR1_Ki_WT_ChEMBL_134 |
| 4tz4 | CRBN_Ki_WT_ChEMBL_127 |
| 4zau | EGFR_Ki_WT_curated_251 |
| 5jkv | CYP19A1_Aromatase_Ki_WT_ChEMBL_548 |
| 5mo4 | ABL1_BCR-ABL_Ki_WT_ChEMBL_693 |
| 6gqj | KIT_Ki_WT_curated_1298 |
| 6jok | PDGFRA_Ki_WT_curated_250 |
| 5lf3 | PSMB5_Ki_WT_ChEMBL_88 |
| 7kk3 | PARP1_Ki_WT_ChEMBL_1075 |

The merge script also maps extra PDB IDs (for example 8zyq/hERG, 1ere/ERalpha) that are **not** in `TARGETS`.

## Merged analysis tables

`analysis/excluding_2z5x_3mjg/tables/merged_ligands_docking_<pdb>.csv` (join key: `ligand_id`)

| Column | Meaning |
|--------|---------|
| `ligand_id` | Stable id (`ligand_0001`, …); merge key with `protein_pdb_id` |
| `molecule_chembl_id` | Attribute only; not the file join key |
| `standard_value` | \(K_i\) nM |
| `pValue` / `pKi` | \(\mathrm{p}K_i = 9-\log_{10}(K_i\,\mathrm{nM})\) |
| `is_active` | \(K_i < 1000\) nM |
| `boltz2_affinity_pred_value` | Boltz-2 primary |
| `dynamicbind_affinity_bestpose` | DynamicBind primary |
| `gnina_cnn_affinity_bestpose` | GNINA 1.3 primary |
| `plapt_affinity` | PLAPT primary |
| `qvina_affinity_bestpose` | QVina2 primary |

Effective \(N\) is complete-case per method (missing docking rows are not imputed).

Frozen panel \(N\) (curated / merged rows) is **10 711**.
Do not freeze 10 497. The ligand_id merge job writes the same count into
`input/ligands_curated/N_panel.txt` and Zenodo record 04.

## Activity classes (code)

`src/analysis/data.py` / `add_pvalue_column.py`:

- `high`: \(K_i \le 100\) nM
- `medium`: \(100 < K_i < 1000\) nM
- `low` (inactive): \(K_i \ge 1000\) nM
- screening **active** = high ∪ medium
