#!/bin/bash
#SBATCH --job-name=14_votu_abundance
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=2:00:00
#SBATCH --array=1-11%4
# Resources measured, not guessed. Smoke test 73931 mapped 200,000 read pairs against the
# 11 MB reference in 3 s on 4 cores. The environmental samples hold roughly 21-56 M pairs,
# i.e. 100-280x that, so 8 cores should finish each sample in well under 15 min; 2 h is
# generous headroom for the I/O of reading a 14 GB FASTQ off scratch. 24 GB is far more
# than CoverM needs for an 11 MB reference with streamed reads, and %4 keeps four tasks
# running without monopolising a shared node.
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/14_votu_abundance_%A_%a.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/14_votu_abundance_%A_%a.err

# How abundant is each vOTU in the environment, month by month?
#
# Maps each ENVIRONMENTAL sample's cleaned reads against the 19 vOTU representative
# genomes with CoverM, giving per-vOTU coverage, breadth and normalised abundance. One
# array task per environmental sample.
#
# WHAT COVERM DOES, and why it is the right tool: it aligns the reads (minimap2 short-read
# preset, bundled), then summarises coverage per GENOME rather than per contig. That
# matters because a GVMAG is many contigs, and a per-contig table would give 19 vOTUs
# hundreds of uninterpretable rows. It also reports `covered_fraction`, the proportion of
# the genome with any coverage, which is what the MIUViG standard uses to decide whether a
# vOTU is actually PRESENT rather than merely attracting a few stray reads.
#
# Read filters, and why these values:
#   --min-read-percent-identity 95  a read must match at >=95% identity. 95% ANI is the
#       MIUViG species boundary and the same cut-off used to build the vOTUs in step 09, so
#       this keeps reads from a related-but-different virus from inflating a vOTU's
#       abundance. Without it, closely related vOTUs steal reads from each other.
#       NOTE THE SCALE: CoverM 0.7.0 wants a PERCENTAGE here ("e.g. 95 for 95%"), not a
#       fraction. Passing 0.95 is accepted silently and means 0.95% identity - i.e. no
#       filtering whatsoever - which would inflate every abundance through cross-mapping
#       with no error to warn you. Verified against `coverm genome --full-help` for 0.7.0.
#   --min-read-aligned-percent 75   at least 75% of the read's bases must align, which
#       rejects reads matching only through a short conserved patch. Percentage, as above.
#   --min-covered-fraction 0          deliberately DISABLED here. CoverM's default (10)
#       silently rewrites a genome's abundance to 0 when breadth is low, which would bake a
#       presence/absence threshold into the raw numbers. Reporting raw values and applying
#       the >=75% breadth rule later (merge_votu_abundance.py) keeps the threshold explicit,
#       documented and changeable without re-mapping ~200 GB of reads.
#
# Breadth is requested as `covered_bases` plus `length`, NOT as a "covered_fraction"
# method: CoverM 0.7.0 has no such method (its methods are relative_abundance, mean,
# trimmed_mean, coverage_histogram, covered_bases, variance, length, count,
# reads_per_base, rpkm, tpm). merge_votu_abundance.py divides the two, which is the same
# number and makes the arithmetic visible.
#
# Each vOTU is represented by ONE genome (its step 09 representative), not all its members,
# so reads are not split between near-identical genomes of the same vOTU.
#
# Submit (after 14_abundance_setup.slurm, and after running 14_votu_reference.sh):
#   ( module purge; module load slurm/kuro; sbatch -p serial scripts/14_votu_abundance.slurm )

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:?Submit with sbatch}"
DATA="$HOME/scr10/data/kbrevis"
REF="$DATA/abundance/votu_ref"
OUT="$DATA/abundance/per_sample"
CLEAN="$DATA/clean"

# The 11 environmental samples, in collection-date order (oldest first) so the array index
# reads naturally in the logs. Dates come from metadata/SraRunTable.csv (Collection_Date).
#   SRR28004153  2020-12   Gulf of Mexico      SRR28004148  2021-05  Gulf of Mexico
#   SRR28004152  2021-01   Gulf of Mexico      SRR28004147  2021-06  Gulf of Mexico
#   SRR28004151  2021-02   Gulf of Mexico      SRR28004146  2021-07  Gulf of Mexico
#   SRR28004150  2021-03   Gulf of Mexico      SRR28004145  2021-08  Gulf of Mexico
#   SRR28004149  2021-04   Gulf of Mexico      SRR34748664  2021-09-29  coastal Florida
#                                              SRR34748663  2022-03-08  coastal Florida
SAMPLES=(SRR28004153 SRR28004152 SRR28004151 SRR28004150 SRR28004149 SRR28004148
         SRR28004147 SRR28004146 SRR28004145 SRR34748664 SRR34748663)

NAME="${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}"
R1="$CLEAN/${NAME}_1_CLEAN.fastq"
R2="$CLEAN/${NAME}_2_CLEAN.fastq"
TSV="$OUT/${NAME}_coverm.tsv"

mkdir -p "$OUT"
[[ -d "$REF" ]] || { echo "No vOTU reference ($REF); run scripts/14_votu_reference.sh first." >&2; exit 1; }
[[ -s "$R1" && -s "$R2" ]] || { echo "Missing cleaned reads for $NAME" >&2; exit 1; }

# Restartable: a finished sample is skipped, so the array can be resubmitted safely.
if [[ -s "$TSV" ]]; then
    echo "[abund] $NAME already done ($TSV); skipping"
    exit 0
fi

# Node-local scratch for CoverM's temporary BAMs - never /tmp (AGENTS.md section 5).
export TMPDIR="/local/scr/$USER/coverm_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate coverm-env

echo "[abund] sample $NAME (task $SLURM_ARRAY_TASK_ID), $(ls "$REF"/*.fna | wc -l) vOTU genomes, $THREADS threads"
echo "[abund] reads: $(du -h "$R1" | cut -f1) + $(du -h "$R2" | cut -f1)"
coverm --version

# Write to a temporary name and mv into place, so an interrupted run never leaves a
# half-written table that the restart check would mistake for a finished one.
coverm genome \
    --coupled "$R1" "$R2" \
    --genome-fasta-directory "$REF" \
    --genome-fasta-extension fna \
    --methods mean covered_bases length count rpkm tpm relative_abundance variance \
    --min-read-percent-identity 95 \
    --min-read-aligned-percent 75 \
    --min-covered-fraction 0 \
    --threads "$THREADS" \
    --output-file "${TSV}.part"
mv "${TSV}.part" "$TSV"

conda deactivate
echo "[abund] wrote $TSV ($(wc -l < "$TSV") lines)"
