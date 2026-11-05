# AffiTox input data

One parent directory holds proteins and ligands:

- `input/proteins/<pdb>.pdb` — RCSB PDB snapshots for the public 16-target panel
- `input/ligands/` — raw ChEMBL API extracts (layer 1)
- `input/ligands_nodubl/` — frozen curated tables (row order used for docking)
- `input/ligands_curated/<pdb>_ligands.csv` — same rows plus `ligand_id`, `pKi`, `is_active`

Do not ship 7awe, 2z5x, or 3mjg as panel members. Extra PDBs may remain on disk for provenance.

```bash
source config/project.env.sh
bash run_pipeline.sh curate    # default: freeze nodubl order and assign ligand_0001…
# bash pipeline/prepare/run_curate.sh --mode from-raw   # apply Ki/MW filters to layer 1
```

Join key for docking files and merge: `ligand_id` (`ligand_0001`, …). `molecule_chembl_id` and SMILES are attributes. Sidecar maps: `processed/id_maps/`.

Panel PDB IDs: `TARGETS` in `config/project.env.sh`.

| Path / field | Meaning |
|--------------|---------|
| `proteins/<pdb>.pdb` | target structure |
| `ligands/*.csv` | ChEMBL API snapshot (semicolon) |
| `ligands_curated/<pdb>_ligands.csv` `ligand_id` | stable id assigned before docking |
| `standard_value` | \(K_i\) nM |
| `pKi` / `pValue` | \(9 - \log_{10}(K_i[\mathrm{nM}])\) |
| `is_active` | \(K_i < 1000\) nM |

N is the curated / merged panel size. See `input/ligands_curated/N_panel.txt` (10711).

License: MIT for code; ChEMBL/PDB terms for source records (`NOTICE`).
