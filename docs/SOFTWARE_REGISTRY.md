# AffiTox software pins for full docking reproduction

Analysis-level reproduction from deposited tables does **not** require these binaries.
Full docking does. The table below mixes three evidence sources:

1. **Adapter / YAML**: what AffiTox *intends* to call (`src/docking_benchmark2/docking/`, `config/methods_config.yaml`).
2. **Published logs**: what a deposited run actually printed (local `results/` mirrors Zenodo raw parts).
3. **Current cluster checkout**: git HEAD / `gnina --version` / conda *today*. This can drift after the paper run.

Do not treat (3) as the paper pin unless it matches (2).

Set install locations with `PLAPT_PATH`, `DYNAMICBIND_PATH`, `BOLTZ_ROOT`, `GNINA_BIN`, `BUST_BIN`, `POSEBUSTERS_CONFIG`, or copy `config/methods_config.hpc.example.yaml` and export `METHODS_CONFIG`.

## Summary

| Method | Recovered from published logs / deposited CSVs | Current cluster snapshot (may have drifted) | Still open |
|--------|-------------------------------------------------|---------------------------------------------|------------|
| QVina2 | Log banner = Quick Vina 2 (Alhossary 2015). Example `4f65/ligand_1.log`: 9 modes, exhaustiveness described as "low", **printed seed is not 42**. Logs do **not** store the argv. | `qvina02` is **not** in conda env `docking` (that env has AutoDock Vina **1.2.6**). | Confirm which `qvina02` binary wrote the BindingDB logs; QVina2 does not print a version string. |
| GNINA 1.3 | Log header **`gnina v1.3 master:97fa6bc+` Built Oct 3 2024**. Full argv with `--no_gpu`, no `--cnn`. | Authors' cluster binary named `gnina` (needs `libcudart` from a GNINA env). Default `--cnn_scoring=rescore`. | Named CNN weights are the GNINA 1.3 built-in default (not passed on the CLI). |
| DynamicBind | Log argv: `run_single_protein_inference.py`, `--inference_steps 20 --seed 42 --no_relax`. **`--samples_per_complex 5 --savings_per_complex 1`** (not 10). No `--rigid_protein` in the logged command. CWD = DynamicBind checkout. | git `abdcd83` (`v1.0-6-gabdcd83`), origin `luwei0917/DynamicBind`. Checkpoints in `workdir/big_score_model_sanyueqi_with_time/`. | YAML in the repo currently lists `samples_per_complex: 10`; **logs are the paper run**. |
| Boltz-2 | Run tag `boltz2_s80_d10_seed42` / `boltz2_s80_d10_seed42_aff_s80_d1`. SLURM script: `--model boltz2 --seed 42 --sampling_steps 80 --diffusion_samples 10 --sampling_steps_affinity 80 --diffusion_samples_affinity 1 --affinity_mw_correction`. YAML `version: 1` with protein MSA + ligand SMILES + CIF template + affinity binder. | Package **boltz 2.2.0**, git `ac33fe0` (`v2.2.0-18-gac33fe0`). Weights `~/.boltz/boltz2_conf.ckpt` and `boltz2_aff.ckpt`. | Confirm those `.ckpt` files were not overwritten after the run. |
| PLAPT | Per-ligand JSON with `affinity` / `affinity_uM` (no version field). Adapter loads `WELP-PLAPT/plapt.py` + `models/affinity_predictor.onnx`. | git `ebc6391` (`trrt-good/WELP-PLAPT`). ONNX MD5 `4a340ab3a3179417ce41f3f44328ee38`. | JSON does not record the git commit. |
| PoseBusters | Current record 05 ([10.5281/zenodo.23979721](https://doi.org/10.5281/zenodo.23979721)) CSVs have **21** dock-style boolean columns including `internal_energy` and `non-aromatic_ring_non-flatness`. No version column. CLI used: `bust -t <table.csv> --outfmt csv --output <results.csv>` (config optional). | Current env `posebuster`: **bust 0.4.6**. Configs `dock.yml` / `dock_fast.yml` present. | CSV has no `bust --version`. Older 05 deposits mixed 20- and 21-column files; use the current DOI above. |

---

## QVina2: example log (no argv)

Source: `results/4f65/docking/qvina/ligand_1.log` (same files ship in Zenodo raw `results_4f65.zip`).

The engine identifies itself only by the Quick Vina 2 citation banner. There is no `Commandline:` line (unlike GNINA).

```text
# Fast, Accurate, and Reliable Molecular Docking with QuickVina 2,
# Bioinformatics (2015), doi: 10.1093/bioinformatics/btv082
WARNING: The search space volume > 27000 Angstrom^3 (See FAQ)
WARNING: at low exhaustiveness, it may be impossible to utilize all CPUs
Using random seed: 1102758008
mode |   affinity | dist from best mode
   1         -7.5      0.000      0.000
   ...
   9         -6.1     21.482     23.697
```

Adapter reconstruction (`src/docking_benchmark2/docking/qvina.py` + YAML):

```text
qvina02 --receptor <protein.pdbqt> --ligand <ligand.pdbqt>
  --center_x/y/z <box> --size_x/y/z <box>
  --out <ligand>_out.pdbqt --log <ligand>.log
  --exhaustiveness 8 --seed 42 --num_modes 10 --energy_range 10
```

**Caveat:** deposited logs for 4f65/1g5m print seeds other than 42 and list 9 modes. Either the original run did not pass `--seed 42`, or this `qvina02` build ignores/rewrites the seed. Reconstruct from logs, not from the current YAML alone.

---

## GNINA 1.3: example command line (authoritative)

Source: `results/4f65/docking/gnina/ligand_1.log`.

```text
gnina v1.3 master:97fa6bc+   Built Oct  3 2024.
gnina is based on smina and AutoDock Vina.
Commandline: gnina
  --receptor .../proteins/4f65.pdbqt
  --ligand   .../ligands/4f65/FGFR1_Ki_WT_ChEMBL_134_nodubl/ligand_1.pdbqt
  --out      .../gnina/ligand_1_out.pdbqt
  --log      .../gnina/ligand_1.log
  --exhaustiveness 8 --num_modes 9 --seed 42
  --center_x 17.688 --center_y -3.447 --center_z -5.069
  --size_x 77.834 --size_y 60.176 --size_z 63.062
  --no_gpu
Using random seed: 42
mode | affinity | intramol | CNN pose score | CNN affinity
   1    -9.85     -0.58        0.9756           9.037
```

No `--cnn` flag: GNINA 1.3 still writes CNN pose/affinity columns (default `cnn_scoring=rescore`, built-in default models). AffiTox primary score is `gnina_cnn_affinity_bestpose`.

---

## DynamicBind: example command line (authoritative)

Source: `results/4f65/docking/dynamicbind_new/FGFR1_Ki_WT_ChEMBL_134_nodubl/4f65_FGFR1_Ki_WT_ChEMBL_134_nodubl_dynamicbind.log`  
(same pattern for `1g5m` / BCL2).

```text
Command: .../envs/dynamicbind/bin/python
  <DynamicBind_checkout>/run_single_protein_inference.py
  <protein.pdb> <ligands.csv>
  --header 4f65_FGFR1_Ki_WT_ChEMBL_134_nodubl
  --results <dataset_dir>
  --samples_per_complex 5
  --savings_per_complex 1
  --inference_steps 20
  --seed 42
  --device 0
  --python .../envs/dynamicbind/bin/python
  --relax_python .../envs/relax/bin/python
  --num_workers 20
  --no_relax
Working directory: <DynamicBind_checkout>
```

Checkpoints expected under that checkout:

```text
workdir/big_score_model_sanyueqi_with_time/ema_inference_epoch314_model.pt
workdir/big_score_model_sanyueqi_with_time/pro_ema_inference_epoch138_model.pt
workdir/big_score_model_sanyueqi_with_time/model_parameters.yml
```

MD5 on the current cluster copies: `24eea910a7cf815b3d93c69010fc2782` and `f304ebf5e394eae9170151b53ea5ee30`. Git commit of that checkout today: `abdcd83f313cd20d50c3917e04615e989a8f63e5`.

---

## Boltz-2: example YAML + CLI

Production tree: `boltz/data/results/4f65/docking/boltz2_s80_d10_seed42/`  
(also mirrored under AffiTox `results/4f65/docking/boltz2_s80_d10_seed42_aff_s80_d1/`).

Example YAML (`yaml_inputs/FGFR1_Ki_WT_ChEMBL_134_nodubl/CHEMBL3893914.yaml`):

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: ELPEDPRWELPRDRLVLGKPLGEGAFGQVVLAEAIGL...
      msa: <msa_cache>/precomputed/a3m/4f65_chainA.a3m
  - ligand:
      id: "134"
      smiles: 'Cc1ccc(OCCNc2nc(Nc3cc(C)[nH]n3)cc(C)c2C#N)cn1'
templates:
  - cif: <boltz_input>/proteins/4f65.cif
    chain_id: A
properties:
  - affinity:
      binder: "134"
```

CLI from `boltz/run_boltz2_docking_bindingdb_slurm.sh` (this is the batch that names the `s80_d10_seed42` tag):

```text
boltz predict <yaml_dir> \
  --model boltz2 \
  --out_dir <pred_dir> \
  --accelerator gpu --devices 1 \
  --seed 42 \
  --sampling_steps 80 \
  --diffusion_samples 10 \
  --sampling_steps_affinity 80 \
  --diffusion_samples_affinity 1 \
  --affinity_mw_correction \
  --max_parallel_samples 1 \
  --num_workers 2
```

Package on this cluster today: `boltz==2.2.0`. Weights: `~/.boltz/boltz2_conf.ckpt`, `~/.boltz/boltz2_aff.ckpt`.

---

## PLAPT: example output record

Source: `results/4f65/docking/plapt/FGFR1_Ki_WT_ChEMBL_134_nodubl/ligand_100.json`.

```json
{
  "protein": "4f65",
  "ligand_id": "ligand_100",
  "smiles": "COC(=O)c1ccc2c(c1)NC(=O)/C2=C(\\Nc1ccc(...",
  "affinity": 6.64884988472049,
  "affinity_uM": 0.22446576628216017,
  "prediction_time": 0.0171463187061139
}
```

No engine version in the JSON. Adapter identity: WELP-PLAPT `plapt.py` + `models/affinity_predictor.onnx` (current checkout `ebc6391`, ONNX MD5 `4a340ab3a3179417ce41f3f44328ee38`).

---

## PoseBusters: deposited CSV schema (no version field)

Source: `analysis/tables/posebuster/posebusters_results_4f65_qvina.csv` (and gnina / dynamicbind_new / boltz2). These tables are in git and in the analysis bundle; they are **not** in the 202 GB raw Zenodo zips alone.

Header (20 boolean checks after `file,molecule,position`):

```text
mol_pred_loaded, mol_cond_loaded, sanitization, inchi_convertible,
all_atoms_connected, bond_lengths, bond_angles, internal_steric_clash,
aromatic_ring_flatness, non-aromatic_ring_non-flatness, double_bond_flatness,
protein-ligand_maximum_distance, minimum_distance_to_protein,
minimum_distance_to_organic_cofactors, minimum_distance_to_inorganic_cofactors,
minimum_distance_to_waters, volume_overlap_with_protein,
volume_overlap_with_organic_cofactors, volume_overlap_with_inorganic_cofactors,
volume_overlap_with_waters
```

That column set matches PoseBusters **dock** configs (`dock.yml` / `dock_fast.yml`) from the 0.3+ series (`non-aromatic_ring_non-flatness` is not in 0.2.x). The CSV does **not** record `bust --version`.

Invocation in `src/analysis/prepare_and_run_posebusters.py`:

```text
bust -t <posebusters_input_*.csv> --outfmt csv --output <posebusters_results_*.csv> --max-workers N
# optional: --config .../posebusters/config/dock_fast.yml
```

Current deposited record 05 tables have 21 boolean checks (`internal_energy` from bust 0.4.x). The CSV still does not store `bust --version`. Do not cite 0.4.6 as a recovered paper pin from logs; it is the schema of the current deposit.

Raw docking archives (Zenodo parts 1 to 6) are listed in `docs/zenodo/affitox_data_manifest.tsv`.
