#!/bin/bash
set -euo pipefail

# we're gonna want the NCBI datasets 

REF_DIR="$HOME/scr10/data/kbrevis/reference"
# set a directory for where the reference genomes will go
 
mkdir -p "$REF_DIR"
 
# --reference restricts the download to the RefSeq reference genome
# --include genome returns just the genomic FASTA
### DOG ### 
datasets download genome taxon "Canis familiaris" \
    --reference \
    --include genome \
    --filename "${REF_DIR}/dog_ref_genome.zip"
unzip -o "${REF_DIR}/dog_ref_genome.zip" -d "$REF_DIR" 

### CAT #### 
datasets download genome accession GCF_018350175.1 \
    --reference \
    --include genome \
    --filename "${REF_DIR}/cat_ref_genome.zip"
unzip -o "${REF_DIR}/cat_ref_genome.zip" -d "$REF_DIR" 

### MOUSE ### 
datasets download genome accession GCF_000001635.27 \
    --reference \
    --include genome \
	--filename "${REF_DIR}/mouse_ref_genome.zip"
unzip -o "${REF_DIR}/mouse_ref_genome.zip" -d "$REF_DIR" 


### HUMAN ### 
datasets download genome accession GCF_000001405.40 \
    --reference \
    --include genome \
	--filename "${REF_DIR}/human_ref_genome.zip"
unzip -o "${REF_DIR}/human_ref_genome.zip" -d "$REF_DIR" \

### K. BREVIS ### 
fasterq-dump SRR14099984 --split-files --outdir "$REF_DIR" --progress
