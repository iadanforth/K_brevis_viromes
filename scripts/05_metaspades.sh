#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
READ_DIR="$HOME/scr10/data/kbrevis/corrected"
OUT_DIR="$HOME/K_brevis_viromes/output/contigs"
SLURM_SCRIPT="$SCRIPT_DIR/05_metaspades.slurm"

[[ -d "$READ_DIR" ]] || {
    echo "Corrected read directory not found: $READ_DIR" >&2
    exit 1
}
[[ -f "$SLURM_SCRIPT" ]] || {
    echo "Slurm script not found: $SLURM_SCRIPT" >&2
    exit 1
}

mkdir -p "$OUT_DIR" "$REPO_DIR/logs"

shopt -s nullglob
R1_FILES=(
    "$READ_DIR"/*_1_*.fastq
    "$READ_DIR"/*_1.fastq
)

((${#R1_FILES[@]} > 0)) || {
    echo "No paired R1 FASTQ files found in $READ_DIR" >&2
    exit 1
}

for R1 in "${R1_FILES[@]}"; do
    case "$R1" in
        *_1_*.fastq)
            R2="${R1/_1_/_2_}"
            FILE_NAME="${R1##*/}"
            NAME="${FILE_NAME%_1_*}"
            ;;
        *_1.fastq)
            R2="${R1%_1.fastq}_2.fastq"
            FILE_NAME="${R1##*/}"
            NAME="${FILE_NAME%_1.fastq}"
            ;;
    esac

    [[ -s "$R1" ]] || {
        echo "Missing or empty R1 file: $R1" >&2
        exit 1
    }
    [[ -s "$R2" ]] || {
        echo "Missing or empty R2 mate for $R1: $R2" >&2
        exit 1
    }

    if [[ -s "$OUT_DIR/${NAME}_contigs.fasta" ]]; then
        echo "[$NAME] already assembled, skipping"
        continue
    fi

    echo "[$NAME] submitting MetaSPAdes job"
    sbatch \
        --job-name="metaspades_${NAME}" \
        "$SLURM_SCRIPT" "$NAME" "$R1" "$R2"
done
