#!/bin/bash
#SBATCH --job-name=01_assemble_kbrevis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/01_assemble_kbrevis_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/01_assemble_kbrevis_%j.err

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate fastp-env
conda env export --no-builds > fastp-env.yml

bash ./scripts/01_assemble_kbrevis.sh

conda deactivate
