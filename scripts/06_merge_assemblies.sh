#!/bin/bash
# Submit one 06_merge_assemblies.slurm job per sample that has BOTH a
# metaviralSPAdes assembly (contigs/) and a metaSPAdes assembly (contigs_no_v/).
# Run this after the 05_metaspades_nov.sh jobs have finished.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
MV_DIR="$HOME/scr10/data/kbrevis/contigs"
MS_DIR="$HOME/scr10/data/kbrevis/contigs_no_v"
OUT_DIR="$HOME/scr10/data/kbrevis/contigs_merged"
SLURM_SCRIPT="$SCRIPT_DIR/06_merge_assemblies.slurm"

mkdir -p "$OUT_DIR" "$REPO_DIR/logs"

shopt -s nullglob
for MV in "$MV_DIR"/*_contigs.fasta; do
    FILE_NAME="${MV##*/}"
    NAME="${FILE_NAME%_contigs.fasta}"
    MS="$MS_DIR/${NAME}_contigs.fasta"

    if [[ ! -s "$MS" ]]; then
        echo "[$NAME] no metaSPAdes assembly yet ($MS), skipping"
        continue
    fi
    if [[ -s "$OUT_DIR/${NAME}_merged_contigs.fasta" ]]; then
        echo "[$NAME] already merged, skipping"
        continue
    fi

    echo "[$NAME] submitting merge job"
    sbatch --job-name="merge_${NAME}" "$SLURM_SCRIPT" "$NAME"
done
