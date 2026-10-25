#!/bin/bash
set -euo pipefail

REF_DIR="$HOME/scr10/data/kbrevis/reference"
CLEAN_DIR="$HOME/scr10/data/kbrevis/reference/kbrevis_clean"

command -v fastp >/dev/null 2>&1 || {
    echo "fastp is not available; activate the environment that contains it." >&2
    exit 1
}

mkdir -p "$CLEAN_DIR"

shopt -s nullglob
R1_FILES=("$REF_DIR"/*_1.fastq)

if ((${#R1_FILES[@]} == 0)); then
    echo "No *_1.fastq files found in $REF_DIR" >&2
    exit 1
fi

for R1 in "${R1_FILES[@]}"; do
    R2="${R1%_1.fastq}_2.fastq"

    if [[ ! -s "$R2" ]]; then
        echo "Missing or empty mate file: $R2" >&2
        exit 1
    fi

    NAME="${R1##*/}"
    NAME="${NAME%_1.fastq}"
    OUT1="$CLEAN_DIR/${NAME}_1_TRIMMED.fastq"
    OUT2="$CLEAN_DIR/${NAME}_2_TRIMMED.fastq"

    echo "Trimming paired reads for $NAME"
    fastp \
        --in1 "$R1" --in2 "$R2" \
        --out1 "$OUT1" --out2 "$OUT2" \
        --trim_front1 8 --trim_front2 8 \
        --trim_tail1 20 --trim_tail2 20 \
        --n_base_limit 0 \
        --length_required 100 \
        --average_qual 20 \
        --json /dev/null \
        --html /dev/null
done

