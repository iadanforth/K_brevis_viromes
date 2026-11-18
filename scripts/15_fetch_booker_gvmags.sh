#!/bin/bash
# Download the GVMAGs published by the original study (Booker et al. 2026) from GenBank
# and reassemble them into per-genome FASTA files.
#
# THE KEY STRUCTURAL POINT: the accessions the paper lists are CONTIGS, not genomes.
#   PV483752-PV483756   5 contigs
#   PV637767-PV637835  69 contigs
# All 74 carry an `isolate` field in their GenBank title naming the GVMAG they belong to,
# and they resolve to just 10 GVMAGs (their vOTU1-vOTU6, some with two member genomes).
# Comparing our whole GVMAGs against individual contigs would give meaningless ANI and
# alignment fractions, so contigs are grouped by isolate into one multi-FASTA per genome -
# the same shape as our own GVMAGs.
#
# Small download (~6 MB over 74 records), so it runs on dev. No email or API key is sent
# to NCBI: esummary/efetch do not require one, and the user's contact details are not ours
# to hand to an external service. Downloaded data lands in its own directory and is only
# ever read as data.
set -euo pipefail

OUT="$HOME/scr10/data/kbrevis/booker_gvmags"
RAW="$OUT/raw"
GENOMES="$OUT/genomes"
mkdir -p "$RAW" "$GENOMES"

ACC=$(python3 -I -c "
a=[f'PV{n}' for n in range(483752,483757)]
b=[f'PV{n}' for n in range(637767,637836)]
print(','.join(a+b))")

# 1. metadata, to learn which contig belongs to which GVMAG
if [[ ! -s "$RAW/esummary.json" ]]; then
    echo "[booker] fetching metadata for 74 accessions"
    curl -sS -A 'kbrevis-pipeline/1.0' \
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&retmode=json&id=${ACC}" \
      -o "$RAW/esummary.json.part"
    mv "$RAW/esummary.json.part" "$RAW/esummary.json"
fi

# 2. the sequences themselves, in two batches to keep each request modest
if [[ ! -s "$RAW/booker_contigs.fasta" ]]; then
    echo "[booker] fetching sequences"
    : > "$RAW/booker_contigs.fasta.part"
    for RANGE in "483752 483756" "637767 637835"; do
        set -- $RANGE
        IDS=$(python3 -I -c "print(','.join(f'PV{n}' for n in range($1,$2+1)))")
        curl -sS -A 'kbrevis-pipeline/1.0' \
          "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id=${IDS}" \
          >> "$RAW/booker_contigs.fasta.part"
        sleep 1   # be polite to NCBI between requests
    done
    mv "$RAW/booker_contigs.fasta.part" "$RAW/booker_contigs.fasta"
fi
echo "[booker] downloaded $(grep -c '^>' "$RAW/booker_contigs.fasta") contigs"

# 3. group contigs into one FASTA per GVMAG, named by the isolate
python3 -I "$HOME/K_brevis_viromes/scripts/group_booker_contigs.py" \
    --esummary "$RAW/esummary.json" \
    --fasta "$RAW/booker_contigs.fasta" \
    --out-dir "$GENOMES" \
    --manifest "$OUT/booker_gvmag_manifest.tsv"
