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

bash ./scripts/02_fastp.sh 

