#!/bin/bash
#SBATCH --job-name=get_ref_genomes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/get_ref_genomes_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/get_ref_genomes_%j.err

set -euo pipefail

module load miniforge3 
conda activate datasets-env
conda env export --no-builds > datasets-env.yml

bash ./scripts/get_ref_genomes.sh
