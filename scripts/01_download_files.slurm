#!/bin/bash
#SBATCH --job-name=kbrevis_download
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/kbrevis_download_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/kbrevis_download_%j.err

set -euo pipefail

# make sure that you're working in the SRA toolkit environment
# Activate sra-toolkit environment
module load miniforge3 
eval "$(conda shell.bash hook)"
conda activate sra-toolkit-env
conda env export --no-builds > sra-toolkit-env.yml

bash ./scripts/01_download_files.sh ~/scr10/data/kbrevis/SraRunTable.csv

conda deactivate 
conda deactivate

rm -r ~/scr10/data/kbrevis/temp
