#!/bin/bash
#SBATCH --partition=serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=2-00:00:00
#SBATCH --mem=360G
# NOTE ON MEMORY, learned the hard way on 2026-11-21 (jobs 122119/122120, both OOM-killed
# after 7-8 h). SRR28004142 and SRR28004143 are the two largest read sets (252 and 236 GB
# of cleaned reads) and need far more than this default: 143 exceeded a 324 GB SPAdes
# limit outright. Give them a 1 TB-class node (bora hima hi03, or gulf gz) with
# --mem=900G.
#
# THE TRAP: do not size these from the "peak so far" in a cancelled run's spades.log.
# Those runs were killed part-way through K55 with memory still climbing, so their 171 GB
# and 310 GB high-water marks were floors, not requirements. Trimming the request to just
# above them is what caused the OOM.
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/05_metaspades_no_v_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/05_metaspades_no_v_%j.err

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:?Submit with sbatch}"
export OMP_NUM_THREADS="$THREADS"

# SPAdes' -m is a cap, in GB, that SPAdes applies to its OWN allocator. It used to be
# hardcoded at 340 while #SBATCH --mem said 360G, and the two have to be kept consistent:
# Slurm's cgroup refuses any allocation past --mem, and SPAdes' internal accounting runs
# roughly 10% below its true RSS (job 92248 reported "325G / 340G" at the moment the
# kernel refused it at 360 GiB RSS). So derive -m from the actual allocation and keep 10%
# headroom. That way raising --mem on the sbatch command line is enough, and SPAdes aborts
# with its own "memory limit exceeded" message instead of dying inside mimalloc.
if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
    SPADES_MEM=$(( SLURM_MEM_PER_NODE / 1024 * 9 / 10 ))
elif [[ -n "${SLURM_MEM_PER_CPU:-}" ]]; then
    SPADES_MEM=$(( SLURM_MEM_PER_CPU * THREADS / 1024 * 9 / 10 ))
else
    echo "Cannot read the Slurm memory allocation; submit with --mem." >&2
    exit 1
fi
(( SPADES_MEM >= 16 )) || { echo "Allocation too small: -m would be ${SPADES_MEM}G" >&2; exit 1; }

# Resume a run that died part-way instead of starting over. Set to a SPAdes check-point
# ('k21', 'k33', 'k55', 'last'); 'k55' redoes only the final k-mer iteration and keeps the
# finished K21/K33 output in $SAMPLE_OUT. -o, -t and -m may be changed on a restart, the
# input reads may not, so the same SAMPLE/R1/R2 arguments must be passed again. Example:
#   RESTART_FROM=k55 sbatch --export=ALL,RESTART_FROM --mem=950G ... 05_metaspades_nov.slurm ...
RESTART_FROM="${RESTART_FROM:-}"

if (($# != 3)); then
    echo "Usage: sbatch $0 SAMPLE_NAME R1.fastq R2.fastq" >&2
    exit 2
fi

NAME="$1"
R1="$2"
R2="$3"
OUT_DIR="$HOME/scr10/data/kbrevis/contigs_no_v"
SAMPLE_OUT="$OUT_DIR/$NAME"

mkdir -p "$OUT_DIR"

[[ -s "$R1" ]] || {
    echo "Missing or empty R1 file: $R1" >&2
    exit 1
}
[[ -s "$R2" ]] || {
    echo "Missing or empty R2 file: $R2" >&2
    exit 1
}

if [[ -s "$OUT_DIR/${NAME}_contigs.fasta" ]]; then
    echo "[$NAME] already assembled, skipping"
    exit 0
fi
if [[ -n "$RESTART_FROM" ]]; then
    [[ -d "$SAMPLE_OUT" ]] || {
        echo "RESTART_FROM=$RESTART_FROM needs the previous run's output directory," >&2
        echo "but $SAMPLE_OUT does not exist. Drop RESTART_FROM to start from scratch." >&2
        exit 1
    }
elif [[ -e "$SAMPLE_OUT" ]]; then
    echo "Output directory already exists for $NAME: $SAMPLE_OUT" >&2
    echo "Move or remove it before resubmitting this sample," >&2
    echo "or set RESTART_FROM=<check-point> to resume it." >&2
    exit 1
fi

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate metaspades-env

command -v metaspades.py >/dev/null 2>&1 || {
    echo "metaspades.py not found in metaspades-env" >&2
    exit 1
}

# Assemble paired metagenomic DNA reads with metaSPAdes.
#
# A restart takes a DIFFERENT argument list, not an extra flag. SPAdes reads the original
# command line back out of $SAMPLE_OUT/params.txt and refuses to be handed it a second
# time: passing -1/-2 aborts with "you cannot specify input data (-1, -2, -12, --pe-1,
# --pe-2 ...) with --restart-from option!" and --only-assembler aborts the same way (see
# check_options_for_restart_from() in spades_pipeline/options_parser.py). Only -o, -t and
# -m may be given; SPAdes logs them as "with updated parameters", which is what makes the
# new -m take effect. The reads still come from params.txt, so the R1/R2 arguments to this
# script are only used for the sanity checks above on a restart.
if [[ -n "$RESTART_FROM" ]]; then
    [[ -s "$SAMPLE_OUT/params.txt" ]] || {
        echo "No $SAMPLE_OUT/params.txt, so there is no run to restart from." >&2
        echo "Remove $SAMPLE_OUT and resubmit without RESTART_FROM." >&2
        exit 1
    }
    echo "[$NAME] resuming from check-point $RESTART_FROM with ${THREADS} threads, -m ${SPADES_MEM}G"
    metaspades.py \
        --restart-from "$RESTART_FROM" \
        -o "$SAMPLE_OUT" \
        -t "$THREADS" \
        -m "$SPADES_MEM"
else
    echo "[$NAME] assembling paired reads with ${THREADS} threads, -m ${SPADES_MEM}G"
    metaspades.py \
        --only-assembler \
        -1 "$R1" \
        -2 "$R2" \
        -o "$SAMPLE_OUT" \
        -t "$THREADS" \
        -m "$SPADES_MEM"
fi

cp "$SAMPLE_OUT/contigs.fasta" "$OUT_DIR/${NAME}_contigs.fasta"
