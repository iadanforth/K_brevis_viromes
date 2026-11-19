#!/bin/bash
#SBATCH --job-name=download_refs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/download_refs_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/download_refs_%j.err

module load miniforge3 
echo "$(conda info --base)"
eval "$(conda shell.bash hook)" 
conda activate datasets-env
# Record the solved environment next to the other env files. This previously wrote
# `datasets-env.yml` into the CURRENT directory - the repository root - leaving a stray
# duplicate of envs/datasets-env.yml in the repo after every run, and the `prefix:` line
# it emits is machine-specific, so it is stripped.
conda env export --no-builds | grep -v '^prefix:' \
    > /sciclone/home/iadanforth/K_brevis_viromes/envs/datasets-env.lock.yml

# NOTE: this wrapper calls the worker script by a RELATIVE path, so the job must be
# submitted from the repository root (sbatch starts the job in the submitting directory).
bash ./scripts/03_download_refs.sh 
