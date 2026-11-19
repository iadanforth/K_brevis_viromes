#!/bin/bash
#SBATCH --job-name=00_create_envs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=8:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/00_create_envs_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/00_create_envs_%j.err

# Create the conda environments that steps 01-06 need, from the pinned envs/*.yml files.
#
# Why this exists: steps 07 and 14 build their own environments, but steps 01-06 assumed
# the environments already existed, because on this account they were created by hand
# before the pipeline was scripted. A new user cloning the repository had no single command
# to build them, which made the first six steps un-runnable without reading section B of
# the README and reconstructing the commands. This closes that gap.
#
# Runs on a compute node: solving environments uses several cores and is explicitly not
# allowed on dev (AGENTS.md section 1). Safe to re-run - existing environments are skipped.
#
# Submit (from the repository root):
#   ( module purge; module load slurm/astral; sbatch -p batch --exclude=as01 scripts/00_create_envs.slurm )

set -euo pipefail
REPO="/sciclone/home/iadanforth/K_brevis_viromes"

module load miniforge3
eval "$(conda shell.bash hook)"

# Only the environments steps 01-06 use. gvmag-bin-env, checkv-env, genomad-env,
# aniclustermap-env, phylo-env and tree-env are built by 07_gvmag_setup.slurm; coverm-env
# and r-plot-env by 14_abundance_setup.slurm. Rebuilding them here would duplicate that
# work and risk the two scripts disagreeing about versions.
#
# NOT included, deliberately: binning-env, bbmpa-env and readcorr-env, which section B of
# the README records as belonging to earlier attempts that are no longer part of the
# pipeline.
ENVS=(sra-toolkit-env datasets-env seqtk-env host-screen-env metaspades-env seqkit-env)

for ENV in "${ENVS[@]}"; do
    YML="$REPO/envs/${ENV}.yml"
    [[ -s "$YML" ]] || { echo "[envs] MISSING $YML" >&2; exit 1; }
    if conda env list | awk '{print $1}' | grep -qx "$ENV"; then
        echo "[envs] $ENV already exists"
    else
        echo "[envs] creating $ENV from envs/${ENV}.yml"
        conda env create -y -f "$YML"
    fi
done

echo
echo "[envs] verifying the tools each step actually calls"
check () {   # check <env> <tool> [tool...]
    local env="$1"; shift
    conda activate "$env"
    for t in "$@"; do
        if command -v "$t" >/dev/null 2>&1; then
            printf "  %-18s %-22s OK\n" "$env" "$t"
        else
            printf "  %-18s %-22s ** MISSING **\n" "$env" "$t"; MISSING=1
        fi
    done
    conda deactivate
}
MISSING=0
check sra-toolkit-env fasterq-dump
check datasets-env    datasets
check seqtk-env       seqtk
check host-screen-env bowtie2 bowtie2-build samtools
check metaspades-env  metaspades.py metaviralspades.py
check seqkit-env      seqkit

# fastp is NOT a conda environment on this account: it is a standalone binary put on PATH
# by ~/.bashrc. A new user must install it themselves - see README section B - because
# sbatch does not source ~/.bashrc and step 02 will fail with "fastp is not available".
if command -v fastp >/dev/null 2>&1; then
    echo "  fastp (standalone)                      OK: $(command -v fastp)"
else
    echo "  fastp (standalone)                      ** NOT ON PATH - step 02 will fail **"
    echo "    install it and make it available to jobs; it is not provided by any envs/*.yml"
    MISSING=1
fi

echo
(( MISSING == 0 )) && echo "[envs] all tools for steps 01-06 are present" \
                   || { echo "[envs] some tools are missing (see above)" >&2; exit 1; }
