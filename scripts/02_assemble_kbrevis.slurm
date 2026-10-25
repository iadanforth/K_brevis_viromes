#!/bin/bash
#SBATCH --job-name=kbrevis_host_screen
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/02_kbrevis_assemble_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/02_kbrevis_assemble_%j.err

set -euo pipefail

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate host-screen-env

for tool in bowtie2 bowtie2-build bowtie2-inspect samtools; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "$tool is not available in the active environment." >&2
        exit 1
    }
done

REF_DIR="$HOME/scr10/data/kbrevis/reference"
TRIM_DIR="$REF_DIR/kbrevis_clean"
OUT_DIR="$HOME/scr10/data/kbrevis/host_filtered"
INDEX_PREFIX="$REF_DIR/host_index"

REF_FASTAS=(
    "$REF_DIR/dog_ref_genome.fna"
    "$REF_DIR/cat_ref_genome.fna"
    "$REF_DIR/mouse_ref_genome.fna"
    "$REF_DIR/human_ref_genome.fna"
)

for ref in "${REF_FASTAS[@]}"; do
    [[ -s "$ref" ]] || {
        echo "Missing or empty reference FASTA: $ref" >&2
        exit 1
    }
done

mkdir -p "$OUT_DIR"

# Build the combined index once; reuse it on later runs.
if ! bowtie2-inspect "$INDEX_PREFIX" >/dev/null 2>&1; then
    HOST_FASTAS=$(IFS=,; printf '%s' "${REF_FASTAS[*]}")
    bowtie2-build --threads "$SLURM_CPUS_PER_TASK" \
        "$HOST_FASTAS" "$INDEX_PREFIX"
fi

shopt -s nullglob
R1_FILES=("$TRIM_DIR"/*_1_TRIMMED.fastq)

if ((${#R1_FILES[@]} == 0)); then
    echo "No trimmed *_1_TRIMMED.fastq files found in $TRIM_DIR" >&2
    exit 1
fi

for R1 in "${R1_FILES[@]}"; do
    R2="${R1%_1_TRIMMED.fastq}_2_TRIMMED.fastq"

    [[ -s "$R2" ]] || {
        echo "Missing or empty mate file: $R2" >&2
        exit 1
    }

    NAME="${R1##*/}"
    NAME="${NAME%_1_TRIMMED.fastq}"

    OUT1="$OUT_DIR/${NAME}_1_hostfree.fastq"
    OUT2="$OUT_DIR/${NAME}_2_hostfree.fastq"

    if [[ -e "$OUT1" || -e "$OUT2" ]]; then
        echo "Output already exists for $NAME; choose another output directory or move the existing files." >&2
        exit 1
    fi

    echo "Screening $NAME against host genomes..."
    bowtie2 -x "$INDEX_PREFIX" \
        -1 "$R1" -2 "$R2" \
        -p "$SLURM_CPUS_PER_TASK" \
        2> "$OUT_DIR/${NAME}_host_screen.log" |
        samtools view -u - |
        samtools fastq -f 12 -F 2304 \
            -1 "$OUT1" \
            -2 "$OUT2" \
            -0 /dev/null -s /dev/null -n -
done
