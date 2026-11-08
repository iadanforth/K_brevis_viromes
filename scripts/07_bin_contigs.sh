#/!/bash
set -euo pipefail

# 1. create environment with bbmap, samtools and metabat 
module load miniforge3
eval "$(conda shell.bash hook)"
conda activate binning-env

set -u

# 2. Check that the tools are available
for tool in bbmap.sh samtools jgi_summarize_bam_contig_depths metabat2; do
    command -v "$tool" >/dev/null 2>&1 || { echo "$tool not found" >&2; exit 1; }
done

# 3. Directories
READ_DIR="$HOME/scr10/data/kbrevis/corrected"
CONTIG_DIR="$HOME/K_brevis_viromes/output/contigs"
OUT_DIR="$HOME/K_brevis_viromes/output/bins"
THREADS="$SLURM_CPUS_PER_TASK"
mkdir -p "$OUT_DIR"

# 4. Find the per-sample contig files made by the metaSPAdes script
shopt -s nullglob
CONTIG_FILES=("$CONTIG_DIR"/*_contigs.fasta)
((${#CONTIG_FILES[@]} > 0)) || { echo "No *_contigs.fasta in $CONTIG_DIR" >&2; exit 1; }
 
for CONTIGS in "${CONTIG_FILES[@]}"; do
    NAME="${CONTIGS##*/}"
    NAME="${NAME%_contigs.fasta}"
    READS="$READ_DIR/${NAME}_CORRECTED_interleaved.fastq"
    SAMPLE_OUT="$OUT_DIR/$NAME"
    BAM="$SAMPLE_OUT/${NAME}_sorted.bam"
    DEPTH="$SAMPLE_OUT/${NAME}_depth.txt"
 
    # Skip samples with no contigs (common with tiny test files)
    if [[ ! -s "$CONTIGS" ]]; then
        echo "[$NAME] contig file is empty, skipping" >&2
        continue
    fi
 
    [[ -s "$READS" ]] || { echo "Missing or empty reads: $READS" >&2; exit 1; }
 
    # Skip samples that already finished
    if [[ -f "$SAMPLE_OUT/.done" ]]; then
        echo "[$NAME] already binned, skipping"
        continue
    fi
 
    mkdir -p "$SAMPLE_OUT"
 
    # 5. Map reads back to this sample's contigs and sort the alignments
    #    nodisk=t builds the contig index in memory instead of writing a ref/ folder
    echo "[$NAME] mapping reads to contigs"
    bbmap.sh -Xmx50g \
        ref="$CONTIGS" nodisk=t \
        in="$READS" interleaved=t \
        out=stdout.sam \
        threads="$THREADS" \
        2> "$SAMPLE_OUT/${NAME}_bbmap.log" |
        samtools sort -@ "$THREADS" -o "$BAM" -
 
    samtools index "$BAM"
 
    # 6. Calculate per-contig read depth
    echo "[$NAME] calculating contig depths"
    jgi_summarize_bam_contig_depths --outputDepth "$DEPTH" "$BAM"
 
    # 7. Bin contigs (default minimum contig length is 2500 bp)
    echo "[$NAME] binning"
    metabat2 \
        -i "$CONTIGS" \
        -a "$DEPTH" \
        -o "$SAMPLE_OUT/${NAME}_bin" \
        -t "$THREADS" \
        --seed 1
 
    N_BINS=$(find "$SAMPLE_OUT" -maxdepth 1 -name "${NAME}_bin.*.fa" | wc -l)
    echo "[$NAME] produced $N_BINS bins"
 
    touch "$SAMPLE_OUT/.done"
done
