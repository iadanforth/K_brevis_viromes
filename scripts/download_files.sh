#!/bin/bash
set -euo pipefail

RUN_TABLE="${1:-$HOME/scr10/data/kbrevis/SraRunTable.csv}"
RAW_DIR="$HOME/scr10/data/kbrevis/raw"
TEMP_DIR="$HOME/scr10/data/kbrevis/temp"

mkdir -p "$RAW_DIR" "$TEMP_DIR"

# Skip the CSV header, then read the first column (Run).
{
    IFS=, read -r _ || true
    while IFS=, read -r accession _ || [[ -n "${accession:-}" ]]; do
        accession="${accession//\"/}"
        accession="${accession%$'\r'}"

        [[ -z "$accession" ]] && continue

        echo "Downloading ${accession}..."
        fasterq-dump "$accession" \
            --split-files \
            --outdir "$RAW_DIR" \
            --temp "$TEMP_DIR" \
            --progress
    done
} < "$RUN_TABLE"
