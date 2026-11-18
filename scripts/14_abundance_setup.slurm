#!/bin/bash
#SBATCH --job-name=14_abundance_setup
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=4:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/14_abundance_setup_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/14_abundance_setup_%j.err

# One-time setup for step 14 (vOTU abundance through time):
#   1. coverm-env  - CoverM 0.7.0 + samtools, for read mapping and coverage
#   2. r-plot-env  - R 4.4 + ggplot2 and friends, for the figure
# Runs on a compute node because conda environment solves are too heavy for dev
# (AGENTS.md section 1). Safe to re-run: existing environments are skipped.
#
# Submit:
#   ( module purge; module load slurm/kuro; sbatch -p serial scripts/14_abundance_setup.slurm )

set -euo pipefail
REPO="/sciclone/home/iadanforth/K_brevis_viromes"

module load miniforge3
eval "$(conda shell.bash hook)"

for ENV in coverm-env r-plot-env; do
    if conda env list | awk '{print $1}' | grep -qx "$ENV"; then
        echo "[setup] $ENV already exists"
    else
        echo "[setup] creating $ENV from envs/$ENV.yml"
        conda env create -y -f "$REPO/envs/$ENV.yml"
    fi
    conda env export -n "$ENV" --no-builds > "$REPO/envs/${ENV}.lock.yml"
done

# Record the exact CoverM interface, because the spelling of its filter flags has changed
# between releases (percent vs fraction) and the mapping script must match the installed one.
conda activate coverm-env
echo "[setup] coverm version: $(coverm --version 2>&1)"
echo "[setup] ---- coverm genome --help (filter flags) ----"
coverm genome --help 2>&1 | grep -iE 'min-read-percent-identity|min-read-aligned-percent|min-covered-fraction|--methods|-m,|genome-fasta' | head -20
conda deactivate

conda activate r-plot-env
Rscript -e 'cat("[setup] R", R.version.string, "\n"); for (p in c("ggplot2","readr","dplyr","tidyr","scales","viridisLite","patchwork")) cat(sprintf("  %-12s %s\n", p, as.character(packageVersion(p))))'
conda deactivate

echo "[setup] done"
