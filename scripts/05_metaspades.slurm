#!/bin/bash
#SBATCH --job-name=06_metaSPAdes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/06_metaspades_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/06_metaspades_%j.err

set -eo pipefail

set -u pipefail 

bash ./scripts/06_metaspades.sh
