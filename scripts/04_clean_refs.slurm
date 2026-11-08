#!/bin/bash
#SBATCH --job-name=clean_metagenomes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/clean_metagenomes_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/clean_metagenomes_%j.err

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate host-screen-env
conda env export --no-builds > host-screen-env.yml

bash ./scripts/04_clean_refs.sh
