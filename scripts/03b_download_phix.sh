#!/bin/bash
set -euo pipefail

module load miniforge3

REF_DIR="$HOME/scr10/data/kbrevis/reference"
PHIX_FASTA="$REF_DIR/phix174.fna"
mkdir -p "$REF_DIR"

if [[ ! -s "$PHIX_FASTA" ]]; then
    echo "Downloading PhiX174 (NC_001422.1)..."
    TMP=$(mktemp -d)
    datasets download virus genome accession NC_001422.1 \
        --filename "$TMP/phix.zip"
    unzip -q "$TMP/phix.zip" -d "$TMP"
    cp "$TMP/ncbi_dataset/data/genomic.fna" "$PHIX_FASTA"
    rm -rf "$TMP"
    echo "Saved $PHIX_FASTA"
else
    echo "$PHIX_FASTA already exists"
fi
