#!/bin/bash
#SBATCH --job-name=14_gvmag_abundance
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=3:00:00
#SBATCH --array=1-11%4
# Resources measured, not guessed. The previous per-vOTU version of this step (job 74276)
# mapped each environmental sample against an 11 MB reference in well under 15 min on
# 8 cores. This version maps against a 14 MB reference (26 genomes, not 19) and runs CoverM
# TWICE per sample, so ~3x the work: 3 h is generous headroom, and `seff` after the run
# should be used to shrink it. 24 GB is far more than CoverM needs for a 14 MB reference
# with streamed reads. %4 keeps four tasks running without monopolising a shared node.
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/14_gvmag_abundance_%A_%a.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/14_gvmag_abundance_%A_%a.err

# Step 14 - how much of each GVMAG is recruited from each environmental sample?
#
# Maps every ENVIRONMENTAL sample's cleaned reads against ALL 26 GVMAGs with CoverM. One
# array task per environmental sample (11 tasks).
#
# WHY THIS REPLACED THE PER-vOTU VERSION (2026-11-23, prompt 60)
# -------------------------------------------------------------
# The first version mapped against 19 vOTU REPRESENTATIVES (one genome per vOTU cluster).
# Booker et al. 2026's Methods instead say: "The relative abundance and coverage of
# environmental metagenomic reads recruited to each GVMAG was calculated using coverm."
# Each GVMAG - so the reference is now all 26 individual GVMAGs, built by
# scripts/14_gvmag_reference.sh. 19 of those came from the incubation experiments (the
# viruses whose presence in the environment is the question) and 7 from environmental
# samples (internal positive controls - see that script's header).
#
# WHAT COVERM DOES, and why it is the right tool: it aligns the reads (minimap2 short-read
# preset, bundled), then summarises coverage per GENOME rather than per contig. That
# matters because a GVMAG is many contigs, and a per-contig table would give 26 genomes
# hundreds of uninterpretable rows. From `covered_bases` and `length` we get BREADTH, the
# proportion of the genome with any coverage, which is what decides whether a genome is
# actually PRESENT rather than merely attracting reads at a few conserved genes.
#
# TWO FILTER SETTINGS, RUN SIDE BY SIDE
# -------------------------------------
# The paper specifies no identity cutoff, no breadth cutoff, no version and no settings, so
# "what they did" can only mean CoverM's defaults. Our previous run used a 95% identity
# filter, and the two give answers that differ by 100-500x in recruited reads (README
# section 14h). Rather than pick one and hide the other, each sample is mapped twice:
#
#   paper   CoverM defaults - NO read-identity filter, NO aligned-length filter.
#           Faithful to Booker et al. This is the primary result.
#   strict  --min-read-percent-identity 95 --min-read-aligned-percent 75
#           95% ANI is the MIUViG species boundary and the same cutoff that defined the
#           vOTUs in step 09, so a read must plausibly come from THIS virus rather than a
#           relative. 75% aligned rejects reads matching only through a short conserved
#           patch.
#
#           NOTE THE SCALE: CoverM 0.7.0 wants a PERCENTAGE here ("e.g. 95 for 95%"), not
#           a fraction. Passing 0.95 is accepted silently and means 0.95% identity - i.e.
#           no filtering at all - with no error to warn you. Verified against
#           `coverm genome --full-help` for 0.7.0.
#
# Comparing the two columns is the whole point: a genome whose recruitment collapses when
# the identity filter is applied was recruiting reads from relatives, not from itself.
#
# --min-covered-fraction 0 in BOTH runs: CoverM's default (10) silently rewrites a genome's
# abundance to 0 when breadth is below 10%, which would bake a presence/absence threshold
# into the raw numbers. Reporting raw values and applying the breadth rule later
# (merge_gvmag_abundance.py) keeps the threshold explicit, documented and changeable
# without re-mapping ~100 GB of reads.
#
# Breadth is requested as `covered_bases` plus `length`, NOT as a "covered_fraction"
# method: CoverM 0.7.0 has no such method (its methods are relative_abundance, mean,
# trimmed_mean, coverage_histogram, covered_bases, variance, length, count,
# reads_per_base, rpkm, tpm). merge_gvmag_abundance.py divides the two, which is the same
# number and makes the arithmetic visible.
#
# CAVEAT, stated plainly: mapping to all 26 GVMAGs individually means the 6 vOTUs that
# contain 2-3 near-identical (>=95% ANI) GVMAGs will SPLIT their reads between those
# members, so each member looks less abundant than the vOTU as a whole. That is a direct
# consequence of following the paper's per-GVMAG design. merge_gvmag_abundance.py therefore
# also writes a per-vOTU rollup that sums read counts across members, so both views exist.
#
# Submit (after 14_abundance_setup.slurm, and after running 14_gvmag_reference.sh):
#   ( module purge; module load slurm/astral; \
#     sbatch -p batch --exclude=as01 scripts/14_gvmag_abundance.slurm )

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:?Submit with sbatch}"
DATA="$HOME/scr10/data/kbrevis"
REF="$DATA/abundance/gvmag_ref"
OUT="$DATA/abundance/per_gvmag"
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

mkdir -p "$OUT/paper" "$OUT/strict"
[[ -d "$REF" ]] || { echo "No GVMAG reference ($REF); run scripts/14_gvmag_reference.sh first." >&2; exit 1; }
NREF=$(ls "$REF"/*.fna 2>/dev/null | wc -l)
[[ "$NREF" -gt 0 ]] || { echo "Reference dir $REF holds no .fna files" >&2; exit 1; }
[[ -s "$R1" && -s "$R2" ]] || { echo "Missing cleaned reads for $NAME" >&2; exit 1; }

# Node-local scratch for CoverM's temporary BAMs - never /tmp (AGENTS.md section 5).
export TMPDIR="/local/scr/$USER/coverm_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate coverm-env

echo "[abund] sample $NAME (task $SLURM_ARRAY_TASK_ID), $NREF GVMAG genomes, $THREADS threads"
echo "[abund] reads: $(du -h "$R1" | cut -f1) + $(du -h "$R2" | cut -f1)"
coverm --version

METHODS=(mean covered_bases length count rpkm tpm relative_abundance variance)

run_coverm () {
    local tag="$1"; shift
    local tsv="$OUT/$tag/${NAME}_coverm.tsv"
    # Restartable: a finished sample/filter pair is skipped, so the array can be
    # resubmitted safely after a failure or a walltime overrun.
    if [[ -s "$tsv" ]]; then
        echo "[abund] $NAME/$tag already done ($tsv); skipping"
        return 0
    fi
    echo "[abund] === $NAME / $tag filter: coverm $* ==="
    local t0=$SECONDS
    # Write to a temporary name and mv into place, so an interrupted run never leaves a
    # half-written table that the restart check would mistake for a finished one.
    coverm genome \
        --coupled "$R1" "$R2" \
        --genome-fasta-directory "$REF" \
        --genome-fasta-extension fna \
        --methods "${METHODS[@]}" \
        --min-covered-fraction 0 \
        --threads "$THREADS" \
        "$@" \
        --output-file "${tsv}.part"
    mv "${tsv}.part" "$tsv"
    echo "[abund] wrote $tsv ($(wc -l < "$tsv") lines) in $((SECONDS - t0))s"
}

# Primary result: CoverM defaults, matching the paper's unspecified settings.
run_coverm paper

# Comparison: the MIUViG-consistent identity and alignment-length filters.
run_coverm strict --min-read-percent-identity 95 --min-read-aligned-percent 75

conda deactivate
echo "[abund] done $NAME"
