#!/bin/bash
#SBATCH --job-name=download_files
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0-02:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/download_files_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/download_files_%j.err

bash ./scripts/download_files.sh ~/scr10/data/kbrevis/SraRunTable.csv
