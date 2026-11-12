#!/bin/bash
# Submit 08_gvmags.slurm as ONE Slurm array job: one task per sample that has merged
# contigs (step 06) and hasn't finished GVMAG building yet.
#
# Usage (from anywhere; load the Slurm module of the cluster you want first):
#   ( module purge; module load slurm/kuro; bash scripts/08_gvmags.sh -p serial )
#   bash scripts/08_gvmags.sh --dry-run          # show what would be submitted
# Anything after the options is passed straight to sbatch (partition, -t, ...).
# Pick the cluster with:  ~/hpc/best_cluster.sh -c 16 -m 96 -H 12
#
# Env vars: MAX_PARALLEL (default 10) = array tasks allowed to run at once.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
DATA="$HOME/scr10/data/kbrevis"
MERGED_DIR="$DATA/contigs_merged"
READ_DIR="$DATA/clean"
OUT_ROOT="$DATA/gvmags"
SLURM_SCRIPT="$SCRIPT_DIR/08_gvmags.slurm"
MAX_PARALLEL="${MAX_PARALLEL:-10}"

DRY_RUN=0
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=1
    shift
fi

# ---- Check the one-time setup (07_gvmag_setup.slurm) was done -------------------
missing=()
[[ -d "$HOME/data10/programs/gvclass-v1.1.1" ]]  || missing+=("GVClass ($HOME/data10/programs/gvclass-v1.1.1)")
[[ -d "$HOME/data10/databases/genomad_db" ]]     || missing+=("geNomad database ($HOME/data10/databases/genomad_db)")
[[ -s "$HOME/data10/databases/checkv_db/checkv-db-v1.5/genome_db/checkv_reps.dmnd" ]] \
                                                 || missing+=("CheckV database v1.5 ($HOME/data10/databases/checkv_db/checkv-db-v1.5)")
if ((${#missing[@]})); then
    printf 'Setup not finished; missing:\n' >&2
    printf '  - %s\n' "${missing[@]}" >&2
    echo "Run scripts/07_gvmag_setup.slurm first." >&2
    ((DRY_RUN)) || exit 1
fi

# ---- Build the sample list -----------------------------------------------------
[[ -d "$MERGED_DIR" ]] || { echo "No merged contigs yet ($MERGED_DIR); run step 06 first." >&2; exit 1; }
mkdir -p "$OUT_ROOT" "$REPO_DIR/logs"

samples=()
shopt -s nullglob
for f in "$MERGED_DIR"/*_merged_contigs.fasta; do
    name="${f##*/}"
    name="${name%_merged_contigs.fasta}"
    if [[ -f "$OUT_ROOT/$name/.done" ]]; then
        echo "[$name] already finished, skipping"
        continue
    fi
    if [[ ! -s "$f" || ! -s "$READ_DIR/${name}_1_CLEAN.fastq" || ! -s "$READ_DIR/${name}_2_CLEAN.fastq" ]]; then
        echo "[$name] merged contigs or clean reads missing/empty, skipping" >&2
        continue
    fi
    samples+=("$name")
done

if ((${#samples[@]} == 0)); then
    echo "Nothing to submit."
    exit 0
fi

# A new, dated list for every submission: the running array reads it line by line,
# so it must not change while the array runs.
LIST="$OUT_ROOT/sample_list_$(date +%Y%m%d_%H%M%S).txt"
if ((DRY_RUN)); then
    LIST="(dry run: list not written)"
else
    printf '%s\n' "${samples[@]}" > "$LIST"
fi

N=${#samples[@]}
CMD=(sbatch --array="1-${N}%${MAX_PARALLEL}" "$@" "$SLURM_SCRIPT" "$LIST")
echo "${N} samples: ${samples[*]}"
echo "Sample list: $LIST"
if ((DRY_RUN)); then
    echo "Would run: ${CMD[*]}"
else
    "${CMD[@]}"
fi
