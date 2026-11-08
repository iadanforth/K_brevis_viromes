#!/bin/bash
set -euo pipefail

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate metaspades-env
conda env export --no-builds > metaspades-env.yml

command -v metaspades.py >/dev/null 2>&1 || { echo "metaspades.py not found" >&2; exit 1; }

READ_DIR="$HOME/scr10/data/kbrevis/corrected"
OUT_DIR="$HOME/K_brevis_viromes/output/contigs"
mkdir -p "$OUT_DIR"

READ_FILES=("$READ_DIR"/*_CORRECTED_interleaved.fastq)
((${#READ_FILES[@]} > 0)) || { echo "No *_CORRECTED_interleaved.fastq in $READ_DIR" >&2; exit 1; }

for READS in "${READ_FILES[@]}"; do
    NAME="${READS##*/}"
    NAME="${NAME%_CORRECTED_interleaved.fastq}"
    SAMPLE_OUT="$OUT_DIR/$NAME"

    # 2. Warn if file is empty
    [[ -s "$READS" ]] || { echo "Missing or empty input: $READS" >&2; exit 1; }

    # Skip samples that already finished
    if [[ -s "$SAMPLE_OUT/contigs.fasta" ]]; then
        echo "[$NAME] already assembled, skipping"
        continue
    fi

    # 3. Run metaSPAdes
    echo "[$NAME] assembling"
    metaspades.py \
        --12 "$READS" \
        -o "$SAMPLE_OUT" \
        -t "$SLURM_CPUS_PER_TASK" \
        -m 110 \
        # you will probably have to change the 122 based on how long your 
        # files end up being when you run this on real data
        -k 21,33,55,77,99,122

    # Copy contigs out with the sample name, for convenience
    cp "$SAMPLE_OUT/contigs.fasta" "$OUT_DIR/${NAME}_contigs.fasta"
done
