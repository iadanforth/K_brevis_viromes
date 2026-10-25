#!/bin/bash
#SBATCH --job-name=clean_files
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/clean_files_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/clean_files_%j.err

set -ueo pipefail

module load miniforge3/
conda env list
conda activate bbmpa-env

bash ./scripts/clean_files.sh
