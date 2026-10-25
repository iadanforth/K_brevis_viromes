#!/bin/bash
set -euo pipefail

RAW_DIR="$HOME/scr10/data/kbrevis/raw"
CLEAN_DIR="$HOME/scr10/data/kbrevis/clean"
REF="phix"  # BBDuk's built-in PhiX reference

command -v bbduk.sh >/dev/null 2>&1 || {
    echo "bbduk.sh is not available; activate the BBMap environment first." >&2
    exit 1
}

mkdir -p "$CLEAN_DIR"

shopt -s nullglob
FASTQ_FILES=("$RAW_DIR"/*.fastq)

# first, check that there are actually .fastq files found in the directory
if ((${#FASTQ_FILES[@]} == 0)); then
    echo "No .fastq files found in $RAW_DIR" >&2
    # the >&2 redirects this to standard error instead of standard output 
    exit 1
fi

# then,run a for loop on all SEQs in the RAW directory 
# make sure that all files have a "mate" (i.e. that there is a 1 and 2 file for each of them)
# BASE="${SEQ##*/}" — removes the directory part of the path. For example, /data/raw/SRR123_1.fastq becomes SRR123_1.fastq.
# BASE="${BASE%_1.fastq}" — removes the _1.fastq suffix, leaving SRR123.
# MATE2="$RAW_DIR/${BASE}_2.fastq" — builds the expected path to its second mate, such as /data/raw/SRR123_2.fastq.

for SEQ in "${FASTQ_FILES[@]}"; do
    if [[ "$SEQ" == *_1.fastq ]]; then
        BASE="${SEQ##*/}"
        BASE="${BASE%_1.fastq}"
        MATE2="$RAW_DIR/${BASE}_2.fastq"

        if [[ ! -f "$MATE2" ]]; then
            echo "Missing mate file for $SEQ: $MATE2" >&2
            exit 1
        fi
# if you have confirmed that the mate is there, proceed to running bbduk.sh
        echo "Processing paired run $BASE..."
        bbduk.sh -Xmx30g \
            in1="$SEQ" in2="$MATE2" \
            out1="$CLEAN_DIR/${BASE}_1_clean.fastq" \
            out2="$CLEAN_DIR/${BASE}_2_clean.fastq" \
            ref="$REF" k=27 hdist=1 qtrim=rl trimq=17 \
            cardinality=t mingc=0.05 maxgc=0.95

    elif [[ "$SEQ" == *_2.fastq ]]; then
        MATE1="${SEQ%_2.fastq}_1.fastq"
        [[ -f "$MATE1" ]] || {
            echo "Missing mate file for $SEQ: $MATE1" >&2
            exit 1
        }

    else
        BASE="${SEQ##*/}"
        BASE="${BASE%.fastq}"

        echo "Processing single-end file $BASE..."
        bbduk.sh -Xmx30g \
            in="$SEQ" out="$CLEAN_DIR/${BASE}_clean.fastq" \
            ref="$REF" k=27 hdist=1 qtrim=rl trimq=17 \
            cardinality=t mingc=0.05 maxgc=0.95
    fi
done
