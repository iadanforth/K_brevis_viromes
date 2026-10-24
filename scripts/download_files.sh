#!/bin/bash
set -euo pipefail

RUN_TABLE="${1:-$HOME/scr10/data/kbrevis/SraRunTable.csv}"
RAW_DIR="$HOME/scr10/data/kbrevis/raw"

if [[ ! -f "$RUN_TABLE" ]]; then
    echo "Run table not found: $RUN_TABLE" >&2
    exit 1
fi

if ! command -v fasterq-dump >/dev/null 2>&1; then
    echo "fasterq-dump is not available; load or activate the SRA Toolkit environment first." >&2
    exit 1
fi

mkdir -p "$RAW_DIR"

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
            --progress
    done
} < "$RUN_TABLE"
