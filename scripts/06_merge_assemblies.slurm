#!/bin/bash
#SBATCH --partition=serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=12G
# Sized from a measured pilot run (job 110082, SRR28004145 - the largest metaSPAdes
# assembly, 2.18M contigs of which 100,406 pass the 1000 bp filter): 18 seconds and
# 1.68 GB MaxRSS on 8 cores. The original 16c/64G/12h was ~40x the memory and ~700x the
# walltime actually needed; blastn only has a few hundred metaviralSPAdes contigs to
# query, so it finishes in seconds. 4 cores / 12 GB / 1 h still leaves ~7x the memory and
# a large time margin for the samples with more metaviralSPAdes contigs.
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/06_merge_assemblies_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/06_merge_assemblies_%j.err

# Compare metaviralSPAdes and metaSPAdes contigs for one sample and merge them.
# Same contig = >=95% ANI over >=85% of the shorter contig (blastn + CheckV anicalc).
#
# MV_ONLY=1 handles a sample with NO metaSPAdes assembly (SRR28004142, SRR28004143 - their
# metaSPAdes jobs were cancelled for time; see README "05b-fix-c"). Nothing can be compared,
# so the metaviralSPAdes contigs become that sample's final contig set unchanged. It has to
# be asked for explicitly, not inferred from a missing file, so that running this step
# before the metaSPAdes jobs finish can never silently throw an assembly away.

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:-4}"
MIN_LEN="${MIN_LEN:-1000}"   # minimum metaSPAdes contig length kept
MIN_ANI="${MIN_ANI:-95}" # Average nucleotide identity threshold for merging contigs
MIN_COV="${MIN_COV:-85}"
MV_ONLY="${MV_ONLY:-0}"

if (($# != 1)); then
    echo "Usage: sbatch $0 SAMPLE_NAME" >&2
    exit 2
fi

NAME="$1"
SCRIPT_DIR="/sciclone/home/iadanforth/K_brevis_viromes/scripts"
MV_FASTA="$HOME/scr10/data/kbrevis/contigs/${NAME}_contigs.fasta"
MS_FASTA="$HOME/scr10/data/kbrevis/contigs_no_v/${NAME}_contigs.fasta"
OUT_DIR="$HOME/scr10/data/kbrevis/contigs_merged"
WORK="$OUT_DIR/$NAME"

[[ -s "$MV_FASTA" ]] || { echo "Missing or empty: $MV_FASTA" >&2; exit 1; }

if ((MV_ONLY)); then
    # Refuse if a metaSPAdes assembly does exist: then a real merge is what is wanted, and
    # passing MV_ONLY=1 by mistake would quietly discard half the data.
    if [[ -s "$MS_FASTA" ]]; then
        echo "MV_ONLY=1 was given but $MS_FASTA exists." >&2
        echo "Drop MV_ONLY to merge both assemblies for $NAME." >&2
        exit 1
    fi
    echo "[$NAME] ==== metaviralSPAdes-only mode ===="
    echo "[$NAME] no metaSPAdes assembly; its contigs become the final set unchanged."
else
    [[ -s "$MS_FASTA" ]] || { echo "Missing or empty: $MS_FASTA" >&2; exit 1; }
fi

mkdir -p "$WORK"

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate seqkit-env

# Tag contig names with their assembler so the two sets can be told apart. The MV_ prefix
# is applied in both modes so contig IDs mean the same thing across all samples.
seqkit replace -p '^' -r 'MV_' "$MV_FASTA" > "$WORK/mv.fasta"
if ((MV_ONLY)); then
    conda deactivate
    conda activate checkv-env
    # No --ms/--ani: merge_assemblies.py passes the contigs through as one cluster each and
    # labels every record assembly_source=metaviralspades_only. No length filter is applied,
    # which matches how metaviralSPAdes contigs are treated in the merged samples too
    # (MIN_LEN only ever applied to the metaSPAdes side).
    python -I "$SCRIPT_DIR/merge_assemblies.py" \
        --sample "$NAME" \
        --mv "$WORK/mv.fasta" \
        --out-dir "$OUT_DIR" | tee "$OUT_DIR/${NAME}_merge_summary.txt"
    exit 0
fi

seqkit seq -m "$MIN_LEN" "$MS_FASTA" | seqkit replace -p '^' -r 'MS_' > "$WORK/ms.fasta"

conda deactivate
conda activate checkv-env

echo "[$NAME] blastn: metaviralSPAdes contigs vs metaSPAdes contigs"
makeblastdb -in "$WORK/ms.fasta" -dbtype nucl -out "$WORK/ms_db" > /dev/null
blastn \
    -query "$WORK/mv.fasta" \
    -db "$WORK/ms_db" \
    -outfmt '6 std qlen slen' \
    -max_target_seqs 10000 \
    -num_threads "$THREADS" \
    -out "$WORK/mv_vs_ms.blast.tsv"

anicalc -i "$WORK/mv_vs_ms.blast.tsv" -o "$WORK/mv_vs_ms.ani.tsv"

python -I "$SCRIPT_DIR/merge_assemblies.py" \
    --sample "$NAME" \
    --mv "$WORK/mv.fasta" \
    --ms "$WORK/ms.fasta" \
    --ani "$WORK/mv_vs_ms.ani.tsv" \
    --out-dir "$OUT_DIR" \
    --min-ani "$MIN_ANI" \
    --min-cov "$MIN_COV" | tee "$OUT_DIR/${NAME}_merge_summary.txt"

rm -f "$WORK"/ms_db.*
