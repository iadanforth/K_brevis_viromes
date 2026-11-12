#!/bin/bash
#SBATCH --partition=serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=2-00:00:00
#SBATCH --mem=360G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/05_metaspades_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/05_metaspades_%j.err
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/05_metaspades_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/05_metaspades_%j.err

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:?Submit with sbatch}"
export OMP_NUM_THREADS="$THREADS"

if (($# != 3)); then
    echo "Usage: sbatch $0 SAMPLE_NAME R1.fastq R2.fastq" >&2
    exit 2
fi

NAME="$1"
R1="$2"
R2="$3"
OUT_DIR="$HOME/scr10/data/kbrevis/contigs"
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
if [[ -e "$SAMPLE_OUT" ]]; then
    echo "Output directory already exists for $NAME: $SAMPLE_OUT" >&2
    echo "Move or remove it before resubmitting this sample." >&2
    exit 1
fi

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate metaspades-env

command -v metaviralspades.py >/dev/null 2>&1 || {
    echo "metaviralspades.py not found in metaspades-env" >&2
    exit 1
}

# Assemble paired metagenomic DNA reads with metaviralSPAdes.
echo "[$NAME] assembling paired reads"
metaviralspades.py \
    --only-assembler \
    -1 "$R1" \
    -2 "$R2" \
    -o "$SAMPLE_OUT" \
    -t "$THREADS" \
    -m 340

cp "$SAMPLE_OUT/contigs.fasta" "$OUT_DIR/${NAME}_contigs.fasta"
