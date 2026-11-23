#!/bin/bash
#SBATCH --job-name=14b_abundance_relaxed
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=3:00:00
#SBATCH --array=1-9%3
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/14b_abundance_relaxed_%A_%a.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/14b_abundance_relaxed_%A_%a.err

# SUPERSEDED 2026-11-23 (prompt 60). Kept for provenance only - do not run it.
#
# This script existed to answer one question: how much of the gap between our result and
# Booker et al. 2026 was caused by our 95% read-identity filter, which their Methods do
# not specify? It re-ran step 14 with CoverM's defaults on the nine Gulf samples against
# the 19 vOTU representatives, writing to abundance/per_sample_relaxed/.
#
# Step 14 now does that comparison itself. scripts/14_gvmag_abundance.slurm runs BOTH
# filter settings on ALL 11 environmental samples against ALL 26 individual GVMAGs, and
# merge_gvmag_abundance.py writes them side by side. That is a strict superset of what
# this script produced, so there is nothing left for it to do.
#
# Results from the original run are in README section 14h.


# Sensitivity test: does our 95% read-identity filter explain why we do not see the
# incubation vOTUs in the environment, when Booker et al. 2026 report four of them there?
#
# Their Methods say only "The relative abundance and coverage of environmental metagenomic
# reads recruited to each GVMAG was calculated using coverm" - no identity cutoff, no
# breadth cutoff, no version. CoverM's own default applies NO identity filter. Ours uses
# >=95% identity and >=75% aligned, the MIUViG species boundary, chosen so that a
# related-but-different virus cannot inflate a vOTU's abundance.
#
# That stringency is the right default for "is THIS virus present", but it is exactly wrong
# for "is something LIKE this virus present". If the Gulf carries a relative at 85-94%
# identity rather than the same strain, our filter discards it and theirs would not.
#
# This run repeats step 14 on the nine EH25 monthly samples with CoverM's defaults and
# writes to a SEPARATE directory, so the main result is untouched and the two can be
# compared directly.
set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:?Submit with sbatch}"
DATA="$HOME/scr10/data/kbrevis"
REF="$DATA/abundance/votu_ref"
OUT="$DATA/abundance/per_sample_relaxed"
CLEAN="$DATA/clean"

# The nine Gulf of Mexico monthly samples only - the set the paper analysed (EH25,
# Dec 2020 to Aug 2021). The two coastal Florida samples are a different study.
SAMPLES=(SRR28004153 SRR28004152 SRR28004151 SRR28004150 SRR28004149
         SRR28004148 SRR28004147 SRR28004146 SRR28004145)
NAME="${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}"
TSV="$OUT/${NAME}_coverm.tsv"
mkdir -p "$OUT"
[[ -s "$TSV" ]] && { echo "[relaxed] $NAME already done"; exit 0; }

export TMPDIR="/local/scr/$USER/coverm_relaxed_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"; trap 'rm -rf "$TMPDIR"' EXIT

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate coverm-env

echo "[relaxed] $NAME - CoverM defaults: no identity filter, no aligned-percent filter"
coverm genome \
    --coupled "$CLEAN/${NAME}_1_CLEAN.fastq" "$CLEAN/${NAME}_2_CLEAN.fastq" \
    --genome-fasta-directory "$REF" --genome-fasta-extension fna \
    --methods mean covered_bases length count rpkm tpm relative_abundance \
    --min-covered-fraction 0 \
    --threads "$THREADS" \
    --output-file "${TSV}.part"
mv "${TSV}.part" "$TSV"
echo "[relaxed] wrote $TSV"
