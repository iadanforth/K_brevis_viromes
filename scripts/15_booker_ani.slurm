#!/bin/bash
#SBATCH --job-name=15_booker_ani
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=30:00
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/15_booker_ani_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/15_booker_ani_%j.err

# Do any of OUR GVMAGs match the GVMAGs the original study published?
#
# Compares our 26 GVMAGs against the 10 GVMAGs of Booker et al. 2026 (GenBank
# PV483752-PV483756 and PV637767-PV637835, reassembled from 74 contigs by
# 15_fetch_booker_gvmags.sh) with fastANI, and summarises the result as CSV.
#
# fastANI is run in BOTH directions. It reports the aligned fraction relative to the
# QUERY, so to express alignment as a fraction of the SHORTER genome - which is what the
# MIUViG species definition uses, and what avoids a small genome looking like a full match
# to a large one - the run where the shorter genome is the query is the one to read.
#
# Two things to expect in the output:
#   * fastANI reports nothing below ~80% ANI, so an absent pair means "no relative within
#     ~80% ANI", not "0% identity". Our GVMAGs with no match at all still appear in the
#     summary, with NA.
#   * "same species" uses the MIUViG rule (>=95% ANI over >=85% of the shorter genome),
#     the same rule used to build our own vOTUs in step 09, so the answer is consistent
#     with the rest of this project.
#
# Submit:
#   ( module purge; module load slurm/astral; sbatch -p batch --exclude=as01 scripts/15_booker_ani.slurm )

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:?Submit with sbatch}"
REPO="/sciclone/home/iadanforth/K_brevis_viromes"
DATA="$HOME/scr10/data/kbrevis"
MINE="$DATA/votus/gvmags"
THEIRS="$DATA/booker_gvmags/genomes"
OUT="$DATA/booker_gvmags/ani"

[[ -d "$THEIRS" ]] || { echo "No published GVMAGs ($THEIRS); run scripts/15_fetch_booker_gvmags.sh first." >&2; exit 1; }
mkdir -p "$OUT"

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate aniclustermap-env

ls "$MINE"/*.fna   > "$OUT/mine.txt"
ls "$THEIRS"/*.fna > "$OUT/theirs.txt"
echo "[ani] ours: $(wc -l < "$OUT/mine.txt")   published: $(wc -l < "$OUT/theirs.txt")"

echo "[ani] fastANI, ours as query"
fastANI --ql "$OUT/mine.txt" --rl "$OUT/theirs.txt" -t "$THREADS" \
        -o "$OUT/mine_vs_theirs.tsv" 2> "$OUT/fastani_fwd.log"
echo "[ani] fastANI, published as query (for alignment fraction of the shorter genome)"
fastANI --ql "$OUT/theirs.txt" --rl "$OUT/mine.txt" -t "$THREADS" \
        -o "$OUT/theirs_vs_mine.tsv" 2> "$OUT/fastani_rev.log"
echo "[ani] pairs found: $(wc -l < "$OUT/mine_vs_theirs.tsv") forward, $(wc -l < "$OUT/theirs_vs_mine.tsv") reverse"
conda deactivate

conda activate checkv-env   # for its Python 3.12
python -I "$REPO/scripts/summarize_booker_ani.py" \
    --forward "$OUT/mine_vs_theirs.tsv" \
    --reverse "$OUT/theirs_vs_mine.tsv" \
    --mine-dir "$MINE" \
    --theirs-manifest "$DATA/booker_gvmags/booker_gvmag_manifest.tsv" \
    --votus "$REPO/output/votus/gvmags_votus.tsv" \
    --summaries "$DATA"/gvmags/*/*_gvmag_summary.tsv \
    --out-summary "$REPO/output/booker_gvmag_ani_summary.csv" \
    --out-pairs   "$REPO/output/booker_gvmag_ani_all_pairs.csv"
conda deactivate
echo "[ani] done"
