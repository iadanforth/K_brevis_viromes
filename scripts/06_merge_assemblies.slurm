#!/bin/bash
#SBATCH --partition=serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/06_merge_assemblies_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/06_merge_assemblies_%j.err

# Compare metaviralSPAdes and metaSPAdes contigs for one sample and merge them.
# Same contig = >=95% ANI over >=85% of the shorter contig (blastn + CheckV anicalc).

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:-4}"
MIN_LEN="${MIN_LEN:-1000}"   # minimum metaSPAdes contig length kept
MIN_ANI="${MIN_ANI:-95}"
MIN_COV="${MIN_COV:-85}"

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

for f in "$MV_FASTA" "$MS_FASTA"; do
    [[ -s "$f" ]] || { echo "Missing or empty: $f" >&2; exit 1; }
done

mkdir -p "$WORK"

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate seqkit-env

# Tag contig names with their assembler so the two sets can be told apart.
seqkit replace -p '^' -r 'MV_' "$MV_FASTA" > "$WORK/mv.fasta"
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
