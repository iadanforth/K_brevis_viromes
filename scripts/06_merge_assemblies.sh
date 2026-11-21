#!/bin/bash
# Submit one 06_merge_assemblies.slurm job per sample that has BOTH a
# metaviralSPAdes assembly (contigs/) and a metaSPAdes assembly (contigs_no_v/).
# Run this after the 05_metaspades_nov.sh jobs have finished.
#
# Exception: the samples listed in MV_ONLY_SAMPLES have no metaSPAdes assembly and never
# will (their metaSPAdes jobs were cancelled for time - README "05b-fix-c"), so they are
# submitted with MV_ONLY=1 and their metaviralSPAdes contigs are used as-is. Any OTHER
# sample missing a metaSPAdes assembly is still skipped, because that means step 05b has
# not finished for it yet. If one of these samples is ever reassembled with metaSPAdes, the
# slurm script refuses MV_ONLY=1 and the sample should be taken off this list for a normal
# merge. Override for a one-off run with, e.g., MV_ONLY_SAMPLES="" bash 06_merge_assemblies.sh
set -euo pipefail

# Samples with no metaSPAdes assembly, merged from metaviralSPAdes contigs alone.
#
# EMPTY BY DEFAULT since 2026-11-21. SRR28004142 and SRR28004143 were on this list because
# their metaSPAdes jobs were cancelled for time on 2026-11-13, which left them with 364 and
# 828 contigs and only 12 and 11 bins, against 41-95 bins for every sample that got both
# assemblers. Those two assemblies are being completed, after which all 20 samples get a
# real merge and this list should stay empty.
#
# If a sample genuinely has no metaSPAdes assembly, put it back:
#   MV_ONLY_SAMPLES="SRRxxxxxxx" bash scripts/06_merge_assemblies.sh
MV_ONLY_SAMPLES="${MV_ONLY_SAMPLES-}"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
MV_DIR="$HOME/scr10/data/kbrevis/contigs"
MS_DIR="$HOME/scr10/data/kbrevis/contigs_no_v"
OUT_DIR="$HOME/scr10/data/kbrevis/contigs_merged"
SLURM_SCRIPT="$SCRIPT_DIR/06_merge_assemblies.slurm"

mkdir -p "$OUT_DIR" "$REPO_DIR/logs"

# Collected so run_pipeline.sh can chain the next stage after all of these.
SUBMITTED=()

shopt -s nullglob
for MV in "$MV_DIR"/*_contigs.fasta; do
    FILE_NAME="${MV##*/}"
    NAME="${FILE_NAME%_contigs.fasta}"
    MS="$MS_DIR/${NAME}_contigs.fasta"

    MV_ONLY=0
    if [[ " $MV_ONLY_SAMPLES " == *" $NAME "* ]]; then
        MV_ONLY=1
    fi

    if [[ ! -s "$MS" ]] && ((!MV_ONLY)); then
        echo "[$NAME] no metaSPAdes assembly yet ($MS), skipping"
        continue
    fi
    if [[ -s "$OUT_DIR/${NAME}_merged_contigs.fasta" ]]; then
        echo "[$NAME] already merged, skipping"
        continue
    fi

    if ((MV_ONLY)); then
        echo "[$NAME] submitting metaviralSPAdes-only job (no metaSPAdes assembly)"
        JID=$(sbatch --parsable --job-name="merge_${NAME}" --export=ALL,MV_ONLY=1 "$@" "$SLURM_SCRIPT" "$NAME")
        echo "[$NAME] job $JID"; SUBMITTED+=("$JID")
    else
        echo "[$NAME] submitting merge job"
        JID=$(sbatch --parsable --job-name="merge_${NAME}" "$@" "$SLURM_SCRIPT" "$NAME")
        echo "[$NAME] job $JID"; SUBMITTED+=("$JID")
    fi
done

if ((${#SUBMITTED[@]})); then
    echo "SUBMITTED_JOB_IDS=$(IFS=:; echo "${SUBMITTED[*]}")"
fi
