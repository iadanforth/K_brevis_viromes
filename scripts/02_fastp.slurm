#!/bin/bash
#SBATCH --job-name=kbrevis_fastp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/fastp_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/fastp_%j.err

# Use the pinned conda environment rather than the standalone ~/programs/fastp binary.
# The binary was on PATH only via ~/.bashrc, and `sbatch` does not source ~/.bashrc - this
# step worked solely because a job inherits the submitting shell's environment. That made
# step 02 unreproducible on any other account, and silently dependent on the submitting
# shell here. fastp-env pins the same version (1.3.7), so the trimmed reads are unchanged.
module load miniforge3
eval "$(conda shell.bash hook)"
conda activate fastp-env
fastp --version

# NOTE: this wrapper calls the worker script by a RELATIVE path, so the job must be
# submitted from the repository root (sbatch starts the job in the submitting directory).
bash ./scripts/02_fastp.sh

