#!/bin/bash
set -euo pipefail
# APPLY METSPADES TO ASSEMBLE CONTIGS DE NOVO 
module load miniforge3
eval "$(conda shell.bash hook)"
conda activate metaspades-env

command -v metaspades.py

READ_DIR="$HOME/scr10/data/kbrevis/reference/kbrevis_clean"
OUT_DIR="$HOME/scr10/data/kbrevis/reference/kbrevis_assembly/SRR14099984_metaspades"

[[ -s "$READ_DIR/kbrevis_1_CLEAN.fastq" ]]
[[ -s "$READ_DIR/kbrevis_2_CLEAN.fastq" ]]

metaspades.py \
    -1 "$READ_DIR/kbrevis_reads_1.fastq" \
    -2 "$READ_DIR/kbrevis_reads_2.fastq" \
    -o "$OUT_DIR" \
    -t "$SLURM_CPUS_PER_TASK" \
    -m 110
    
conda deactivate
