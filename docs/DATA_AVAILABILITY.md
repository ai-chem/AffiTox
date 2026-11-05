# Data availability (*Journal of Cheminformatics* / BMC)

Paste the Declarations paragraph into **Availability of data and materials**.
Cite the umbrella dataset as one numbered Vancouver reference. Do not cite
`10.5281/zenodo.20825057`–`067`. The software Zenodo DOI is filled after GitHub
tag `v1.0.0` (CI green). Unpack map: [`zenodo/README.md`](zenodo/README.md).

## Declarations paragraph

The AffiTox v1.0.0 dataset is available from Zenodo
(https://doi.org/10.5281/zenodo.22007993). That index record links six related
deposits: raw ChEMBL and PDB snapshots; prepared inputs, ligand identifiers and
method weights; frozen docking outputs (split across five records because Zenodo
limits a dataset to 50 GiB); merged ligand tables; PoseBusters tables; and
metrics with PNG figures. Source code is available at
https://github.com/ai-chem/AffiTox under the MIT license. ChEMBL and PDB
redistributions remain under their original terms. GPU re-runs of Boltz-2 and
DynamicBind are not bit-identical; the docking deposit is the frozen snapshot
used in this study.

## Vancouver references

Dataset:

Konovalova OA, Orlova A, Telepov A, Khrabrov K, Karpushkina I, Shestun P, Kadurin A, Vinogradov V, Tsypin A, Dmitrenko A. AffiTox v1.0.0: data release manifest [dataset]. Zenodo; 2026. https://doi.org/10.5281/zenodo.22007993.

Code (replace with the Zenodo Software DOI when `v1.0.0` exists):

Konovalova OA, Orlova A, Telepov A, Khrabrov K, Karpushkina I, Shestun P, Kadurin A, Vinogradov V, Tsypin A, Dmitrenko A. AffiTox. GitHub; 2026. https://github.com/ai-chem/AffiTox.

In-text example: “Data are archived on Zenodo [n].”

## BibTeX

```bibtex
@misc{affitox_data_v1,
  author       = {Konovalova, Olga A. and Orlova, Anastasia and Telepov, Alexander
                  and Khrabrov, Kuzma and Karpushkina, Irina and Shestun, Pavel
                  and Kadurin, Artur and Vinogradov, Vladimir and Tsypin, Artem
                  and Dmitrenko, Andrei},
  title        = {{AffiTox} v1.0.0: data release manifest},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.22007993},
  url          = {https://doi.org/10.5281/zenodo.22007993},
  note         = {Dataset. Concept DOI: 10.5281/zenodo.22007992}
}

@misc{affitox_code,
  author       = {Konovalova, Olga A. and Orlova, Anastasia and Telepov, Alexander
                  and Khrabrov, Kuzma and Karpushkina, Irina and Shestun, Pavel
                  and Kadurin, Artur and Vinogradov, Vladimir and Tsypin, Artem
                  and Dmitrenko, Andrei},
  title        = {{AffiTox}},
  year         = {2026},
  publisher    = {GitHub},
  url          = {https://github.com/ai-chem/AffiTox},
  note         = {Version 1.0.0 source. Software Zenodo DOI to be added after the GitHub release.}
}
```

## Part DOIs (supplement only, not the paper citation)

| Record | DOI |
|--------|-----|
| 01 raw ChEMBL+PDB | https://doi.org/10.5281/zenodo.23213812 |
| 02 prepared + weights | https://doi.org/10.5281/zenodo.23979597 |
| 03a docking | https://doi.org/10.5281/zenodo.23981992 |
| 03b docking | https://doi.org/10.5281/zenodo.23982808 |
| 03c docking | https://doi.org/10.5281/zenodo.23983375 |
| 03d docking | https://doi.org/10.5281/zenodo.23984015 |
| 03e docking | https://doi.org/10.5281/zenodo.23984506 |
| 04 merged tables | https://doi.org/10.5281/zenodo.23979579 |
| 05 PoseBusters | https://doi.org/10.5281/zenodo.23979721 |
| 06 metrics + PNG | https://doi.org/10.5281/zenodo.23979767 |
