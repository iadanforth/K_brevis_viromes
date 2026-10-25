#!/bin/bash
set -euo pipefail

# we're gonna want the NCBI datasets 

#!/bin/bash
set -euo pipefail

REF_DIR="$HOME/scr10/data/kbrevis/reference"
mkdir -p "$REF_DIR"

download_and_extract() {
    local name="$1"
    local accession_type="$2"
    local accession="$3"
    local zip="$REF_DIR/${name}_ref_genome.zip"
    local stage

    stage=$(mktemp -d "$REF_DIR/.${name}.XXXXXX")

    if [[ "$accession_type" == "taxon" ]]; then
        datasets download genome taxon "$accession" \
            --reference \
            --include genome \
            --filename "$zip"
    else
        datasets download genome accession "$accession" \
            --include genome \
            --filename "$zip"
    fi

    unzip -q -o "$zip" -d "$stage"

    local -a fasta_files=()
    while IFS= read -r -d '' file; do
        fasta_files+=("$file")
    done < <(find "$stage" -type f -name '*.fna' -print0)

    if ((${#fasta_files[@]} != 1)); then
        echo "Expected one .fna for $name, found ${#fasta_files[@]}; check the package or use a specific assembly accession." >&2
        return 1
    fi

    mv "${fasta_files[0]}" "$REF_DIR/${name}_ref_genome.fna"
    rm -rf -- "$stage"
}

download_and_extract dog taxon "Canis familiaris"
download_and_extract cat accession "GCF_018350175.1"
download_and_extract mouse accession "GCF_000001635.27"
download_and_extract human accession "GCF_000001405.40"

# lastly, delete the zip file 
rm -rf dog_ref_genome.zip
rm -rf cat_ref_genome.zip
rm -rf mouse_ref_genome.zip
rm -rf human_ref_genome.zip


fasterq-dump SRR14099984 \
    --split-files \
    --outdir "$REF_DIR" \
    --progress

mv "$REF_DIR/SRR14099984_1.fastq" "$REF_DIR/kbrevis_reads_1.fastq"
mv "$REF_DIR/SRR14099984_2.fastq" "$REF_DIR/kbrevis_reads_2.fastq"
