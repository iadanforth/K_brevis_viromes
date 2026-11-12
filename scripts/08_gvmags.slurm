#!/bin/bash
#SBATCH --job-name=08_gvmags
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/08_gvmags_%A_%a.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/08_gvmags_%A_%a.err

# Build giant virus metagenome-assembled genomes (GVMAGs) for ONE sample, following
# Booker et al. 2026 (Supplemental Methods):
#   1. map the sample's reads to its merged contigs        (bbmap 35.85, samtools 1.15.1)
#   2. per-contig read depth + bin contigs                 (MetaBAT 2.15)
#   3. find the giant virus bins                           (GVClass)
#   4. quality/completeness of giant virus bins            (GVClass + CheckV 1.1.1, db v1.5)
#   5. gene prediction + annotation, bitscore >= 50        (geNomad 1.7.4, which predicts
#                                                           genes with pyrodigal-gv itself)
# Run as an array by 08_gvmags.sh (one task per line of the sample list),
# or for one sample:  sbatch 08_gvmags.slurm SAMPLE_NAME

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:?Submit with sbatch}"
REPO="/sciclone/home/iadanforth/K_brevis_viromes"
DATA="$HOME/scr10/data/kbrevis"
READ_DIR="$DATA/clean"
MERGED_DIR="$DATA/contigs_merged"
OUT_ROOT="$DATA/gvmags"
GVCLASS_DIR="$HOME/data10/programs/gvclass-v1.1.1"
GENOMAD_DB="$HOME/data10/databases/genomad_db"
CHECKV_DB="$HOME/data10/databases/checkv_db/checkv-db-v1.5"
MIN_BITSCORE=50

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${1:?sample list file}")
else
    NAME="${1:?Usage: sbatch $0 SAMPLE_NAME}"
fi

CONTIGS_IN="$MERGED_DIR/${NAME}_merged_contigs.fasta"
R1="$READ_DIR/${NAME}_1_CLEAN.fastq"
R2="$READ_DIR/${NAME}_2_CLEAN.fastq"
OUT="$OUT_ROOT/$NAME"
for f in "$CONTIGS_IN" "$R1" "$R2"; do
    [[ -s "$f" ]] || { echo "[$NAME] missing or empty: $f" >&2; exit 1; }
done
[[ -d "$GVCLASS_DIR" && -d "$GENOMAD_DB" && -d "$CHECKV_DB" ]] || {
    echo "Databases/GVClass not found; run scripts/07_gvmag_setup.slurm first" >&2; exit 1; }
if [[ -f "$OUT/.done" ]]; then
    echo "[$NAME] already finished, skipping"; exit 0
fi
mkdir -p "$OUT"/{binning,bins,gvclass_in,gvmags,checkv,genomad}

module load miniforge3
eval "$(conda shell.bash hook)"

# ---- 1-2. Map reads, compute depth, bin (gvmag-bin-env) -----------------------
conda activate gvmag-bin-env
# Contig names only (drop the "found_by=..." descriptions) so names match the BAM.
CONTIGS="$OUT/binning/${NAME}_contigs.fna"
seqkit seq -i "$CONTIGS_IN" > "$CONTIGS"

BAM="$OUT/binning/${NAME}.sorted.bam"
if [[ ! -s "$BAM.bai" ]]; then
    echo "[$NAME] mapping reads to merged contigs (bbmap)"
    bbmap.sh -Xmx$(( SLURM_MEM_PER_NODE * 8 / 10 / 1024 ))g \
        ref="$CONTIGS" nodisk=t \
        in1="$R1" in2="$R2" \
        out=stdout.sam threads="$THREADS" \
        2> "$OUT/binning/${NAME}_bbmap.log" |
        samtools sort -@ "$THREADS" -m 2G -o "$BAM" -
    samtools index "$BAM"
fi

echo "[$NAME] contig depths + MetaBAT2 binning"
jgi_summarize_bam_contig_depths --outputDepth "$OUT/binning/${NAME}_depth.txt" "$BAM"
# MetaBAT2 defaults as in the paper (contigs >= 2,500 bp, bins >= 200 kb); fixed seed
# so re-runs give identical bins.
metabat2 -i "$CONTIGS" -a "$OUT/binning/${NAME}_depth.txt" \
    -o "$OUT/bins/${NAME}_bin" -t "$THREADS" --seed 1
conda deactivate

# GVClass wants .fna files with no "." in the base name: SAMPLE_bin.3.fa -> SAMPLE_bin_3.fna
shopt -s nullglob
BINS=("$OUT"/bins/"${NAME}"_bin.*.fa)
echo "[$NAME] MetaBAT2 produced ${#BINS[@]} bins"
if ((${#BINS[@]} == 0)); then
    echo "[$NAME] no bins, nothing else to do"; touch "$OUT/.done"; exit 0
fi
for b in "${BINS[@]}"; do
    base=$(basename "$b" .fa); cp "$b" "$OUT/gvclass_in/${base//./_}.fna"
done

# ---- 3-4a. GVClass: taxonomy, completeness, contamination of every bin --------
echo "[$NAME] GVClass"
export PIXI_HOME="$HOME/data10/programs/pixi" PIXI_CACHE_DIR="$HOME/data10/programs/pixi-cache"
export PATH="$PIXI_HOME/bin:$PATH"
rm -rf "$OUT/gvclass"
( cd "$GVCLASS_DIR" && pixi run gvclass "$OUT/gvclass_in" -o "$OUT/gvclass" -t "$THREADS" )

# Giant virus bins = GVClass domain-level call is NCLDV (phylum Nucleocytoviricota),
# the classification shared by all GVMAGs in the paper (Supplemental Data 4).
conda activate checkv-env   # any env with python3 works for the helper
python -I "$REPO/scripts/select_gvmags.py" \
    --gvclass "$OUT/gvclass/gvclass_summary.tsv" \
    --bin-dir "$OUT/gvclass_in" \
    --out-dir "$OUT/gvmags" \
    --sample "$NAME"
conda deactivate

GVMAGS=("$OUT"/gvmags/*.fna)
echo "[$NAME] ${#GVMAGS[@]} giant virus bins (GVMAGs)"
if ((${#GVMAGS[@]} == 0)); then
    touch "$OUT/.done"; exit 0
fi

# ---- 4b. CheckV 1.1.1 with database v1.5 -------------------------------------
# CheckV scores single sequences, so each GVMAG's contigs are joined into one record
# (separated by 100 Ns) to get one completeness/quality value per GVMAG.
conda activate checkv-env
JOINED="$OUT/checkv/gvmags_joined.fna"
: > "$JOINED"
for g in "${GVMAGS[@]}"; do
    awk -v name="$(basename "$g" .fna)" '
        BEGIN { sp = sprintf("%100s", ""); gsub(/ /, "N", sp); printf ">%s\n", name }
        /^>/  { if (n++) printf "%s", sp; next }
              { printf "%s", $0 }
        END   { print "" }' "$g" >> "$JOINED"
done
checkv end_to_end "$JOINED" "$OUT/checkv/out" -t "$THREADS" -d "$CHECKV_DB" --remove_tmp
conda deactivate

# ---- 5. Gene prediction + annotation (geNomad) --------------------------------
# `genomad annotate` predicts the genes with pyrodigal-gv itself (the paper's gene caller)
# and writes them to <gvmag>_annotate/<gvmag>_proteins.faa, so pyrodigal-gv is NOT called
# separately here: that would predict the same genes twice with the same tool. --cleanup
# only deletes geNomad's mmseqs2 scratch directory, so the proteins are kept.
conda activate genomad-env
for g in "${GVMAGS[@]}"; do
    gname=$(basename "$g" .fna)
    genomad annotate --cleanup --threads "$THREADS" "$g" "$OUT/genomad/$gname" "$GENOMAD_DB"
    # Keep annotations with bitscore >= 50, as in the paper.
    GENES_TSV="$OUT/genomad/$gname/${gname}_annotate/${gname}_genes.tsv"
    awk -F'\t' -v min="$MIN_BITSCORE" \
        'NR==1{for(i=1;i<=NF;i++) if($i=="bitscore") c=i; print; next} c && $c!="NA" && $c>=min' \
        "$GENES_TSV" > "$OUT/genomad/${gname}_genes_bitscore${MIN_BITSCORE}.tsv"
done
conda deactivate

# ---- Summary table (like the paper's Supplemental Data 4) ---------------------
conda activate checkv-env
python -I "$REPO/scripts/summarize_gvmags.py" \
    --sample "$NAME" \
    --gvclass "$OUT/gvclass/gvclass_summary.tsv" \
    --gvmag-dir "$OUT/gvmags" \
    --checkv "$OUT/checkv/out/quality_summary.tsv" \
    --genomad-dir "$OUT/genomad" \
    --min-bitscore "$MIN_BITSCORE" \
    --out "$OUT/${NAME}_gvmag_summary.tsv"
conda deactivate

touch "$OUT/.done"
echo "[$NAME] done: $OUT/${NAME}_gvmag_summary.tsv"
