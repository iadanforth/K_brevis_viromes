#!/bin/bash
#SBATCH --job-name=download_phix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/download_phix_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/download_phix_%j.err

module load miniforge3 
echo "$(conda info --base)"
eval "$(conda shell.bash hook)" 
conda activate datasets-env
conda env export --no-builds > datasets-env.yml

bash ./scripts/03b_download_phix.sh 
