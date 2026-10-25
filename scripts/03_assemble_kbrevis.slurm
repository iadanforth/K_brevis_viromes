#!/bin/bash
#SBATCH --job-name=kb_assembly
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/03_kbrevis_assemble_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/03_kbrevis_assemble_%j.err

set -euo pipefail

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate metaspades-env

READ_DIR="$HOME/scr10/data/kbrevis/host_filtered"
OUT_DIR="$HOME/scr10/data/kbrevis/reference/kbrevis_assembly"

R1="$READ_DIR/kbrevis_reads_1_hostfree.fastq"
R2="$READ_DIR/kbrevis_reads_2_hostfree.fastq"

[[ -s "$R1" ]] || { echo "Missing or empty input: $R1" >&2; exit 1; }
[[ -s "$R2" ]] || { echo "Missing or empty input: $R2" >&2; exit 1; }

metaspades.py \
    -1 "$R1" \
    -2 "$R2" \
    -o "$OUT_DIR" \
    -t "$SLURM_CPUS_PER_TASK" \
    -m 110
