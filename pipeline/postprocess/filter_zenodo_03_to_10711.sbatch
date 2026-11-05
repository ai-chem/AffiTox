#!/bin/bash
#SBATCH --job-name=zenodo03_10711
#SBATCH --output=/nfs/home/okonovalova/AffiTox/release/zenodo/payloads/03_panel10711/filter_%j.out
#SBATCH --error=/nfs/home/okonovalova/AffiTox/release/zenodo/payloads/03_panel10711/filter_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH -p aichem

set -euo pipefail
ROOT=/nfs/home/okonovalova/AffiTox
PY=/mnt/tank/scratch/okonovalova/miniconda3/envs/docking/bin/python
mkdir -p "${ROOT}/release/zenodo/payloads/03_panel10711"
cd "${ROOT}"

# Mark already-filtered zips so --skip-existing-marker works if we requeue.
for pdb in 1g5m 2v5z 3eyg 3jy9 5lf3; do
  marker="${ROOT}/release/zenodo/payloads/03/results_${pdb}.zip.panel10711_done"
  if [[ -f "${ROOT}/release/zenodo/payloads/03/results_${pdb}.zip" && ! -f "${marker}" ]]; then
    echo "{\"pdb\": \"${pdb}\", \"prefiltered\": true}" > "${marker}"
  fi
done

# Remaining targets only (in-place to stay under quota).
"${PY}" pipeline/postprocess/filter_zenodo_03_to_10711.py \
  --in-place \
  --skip-existing-marker \
  --write-buckets \
  --pdb 3lxk --pdb 11ue --pdb 4ase --pdb 4f65 \
  --pdb 4tz4 --pdb 4zau --pdb 5jkv --pdb 5mo4 \
  --pdb 6gqj --pdb 6jok --pdb 7kk3

echo "DONE $(date -Is)"
