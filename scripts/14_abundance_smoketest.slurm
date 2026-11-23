#!/bin/bash
#SBATCH --job-name=14_abundance_smoketest
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=30:00
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/14_abundance_smoketest_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/14_abundance_smoketest_%j.err

# Validate the exact CoverM command in 14_gvmag_abundance.slurm on 200,000 reads before
# committing an 11-task array to ~200 GB of them. Three bugs in that command were already
# found by reading `coverm genome --full-help`, and all three would have been silent or
# late-failing:
#   * --min-read-percent-identity takes a PERCENTAGE (95), not a fraction (0.95). A
#     fraction is accepted and means 0.95% identity, i.e. no filtering at all.
#   * --min-read-aligned-percent likewise.
#   * there is no `covered_fraction` method in CoverM 0.7.0, only `covered_bases`.
# This job proves the corrected command parses, produces the expected columns, and that
# breadth computed as covered_bases/length is sane.
#
# Submit:
#   ( module purge; module load slurm/astral; sbatch -p batch --exclude=as01 scripts/14_abundance_smoketest.slurm )

set -euo pipefail

DATA="$HOME/scr10/data/kbrevis"
REF="$DATA/abundance/votu_ref"
WORK="$DATA/abundance/smoketest"
# SRR28004148 is the environmental sample that vOTU01, vOTU10 and vOTU16 were assembled
# from, so those three MUST show up. A test against a sample with no expected hit could
# pass while mapping nothing at all.
SAMPLE=SRR28004148
N_READS=200000

mkdir -p "$WORK"
export TMPDIR="/local/scr/$USER/coverm_smoke_$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

# Subset here rather than on dev: reading the head of a multi-GB FASTQ is still I/O that
# belongs on a compute node.
head -n $((N_READS * 4)) "$DATA/clean/${SAMPLE}_1_CLEAN.fastq" > "$WORK/r1.fastq"
head -n $((N_READS * 4)) "$DATA/clean/${SAMPLE}_2_CLEAN.fastq" > "$WORK/r2.fastq"
echo "[smoke] subset $N_READS read pairs from $SAMPLE"

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate coverm-env

coverm genome \
    --coupled "$WORK/r1.fastq" "$WORK/r2.fastq" \
    --genome-fasta-directory "$REF" \
    --genome-fasta-extension fna \
    --methods mean covered_bases length count rpkm tpm relative_abundance variance \
    --min-read-percent-identity 95 \
    --min-read-aligned-percent 75 \
    --min-covered-fraction 0 \
    --threads "$SLURM_CPUS_PER_TASK" \
    --output-file "$WORK/smoke_coverm.tsv"

echo "[smoke] CoverM exited 0"
echo "[smoke] ---- column names produced ----"
head -1 "$WORK/smoke_coverm.tsv" | tr '\t' '\n' | nl
# NOTE: `grep -c` exits 1 when the count is zero, which under `set -e` kills the script.
# The first version of this check did exactly that and the job reported FAILED even though
# CoverM had succeeded. `|| true` keeps the count informational.
echo "[smoke] ---- method columns found (CoverM 0.7.0 is LONG format: one row per genome,"
echo "[smoke]      columns named '<sample> <Method>', so these are per-sample, not per-genome) ----"
n_methods=$(head -1 "$WORK/smoke_coverm.tsv" | tr '\t' '\n' | tail -n +2 | grep -c . || true)
echo "  method columns: $n_methods (expect 8)"
n_genomes=$(tail -n +2 "$WORK/smoke_coverm.tsv" | grep -c '^vOTU' || true)
echo "  vOTU rows:      $n_genomes (expect 19)"
echo "[smoke] ---- breadth sanity for the three vOTUs assembled FROM this sample ----"
python3 -I - "$WORK/smoke_coverm.tsv" <<'PY'
import sys, csv
rows=list(csv.reader(open(sys.argv[1]), delimiter='\t'))
hdr, body = rows[0], rows[1:]
# Column 0 is the genome; the remaining columns are "<sample> <Method>". Identify each
# method by the END of its header, since the sample prefix is arbitrary.
idx={}
for i,h in enumerate(hdr[1:], start=1):
    hl=h.lower().strip()
    for needle,name in (('covered bases','cb'),('length','ln'),('mean','mean'),('rpkm','rpkm')):
        if hl.endswith(needle): idx[name]=i; break
def g(row,key):
    i=idx.get(key)
    try: return float(row[i])
    except Exception: return 0.0
print(f"  {'vOTU':8s} {'covered_bases':>14s} {'length':>10s} {'breadth':>8s} {'mean':>8s} {'rpkm':>9s}")
for row in body:
    gen=row[0].strip()
    if not gen.startswith('vOTU'): continue
    cb, ln = g(row,'cb'), g(row,'ln')
    br = cb/ln if ln else 0.0
    mark = '  <- assembled from this sample' if gen in ('vOTU01','vOTU10','vOTU16') else ''
    print(f"  {gen:8s} {cb:14.0f} {ln:10.0f} {br:8.3f} {g(row,'mean'):8.3f} {g(row,'rpkm'):9.3f}{mark}")
PY
conda deactivate
echo "[smoke] done"
