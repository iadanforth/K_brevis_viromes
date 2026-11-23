#!/bin/bash
# Build the mapping reference for step 14: one FASTA per GVMAG.
#
# WHAT CHANGED, AND WHY (2026-11-23, prompt 60 in PROMPT_LOG.md)
# -------------------------------------------------------------
# The first version of this script built a reference of 19 vOTU REPRESENTATIVE genomes -
# one genome per vOTU, the other 7 GVMAGs discarded. That is a defensible design, but it is
# not what Booker et al. 2026 did. Their Methods say:
#
#     "The relative abundance and coverage of environmental metagenomic reads recruited to
#      each GVMAG was calculated using coverm."
#
# "each GVMAG" - not each vOTU. So this script now writes ALL 26 GVMAGs, each as its own
# mapping target, and step 14 maps every environmental sample against all of them.
#
# A GVMAG ("giant-virus metagenome-assembled genome") is one binned, quality-checked viral
# genome recovered from ONE sample. A vOTU is a CLUSTER of GVMAGs that are >=95% identical
# (the MIUViG species boundary, applied in step 09) - so a vOTU is roughly "a viral
# species" and a GVMAG is "one assembled copy of it from one sample". 26 GVMAGs collapsed
# into 19 vOTUs in step 09: 13 vOTUs contain a single GVMAG, 6 contain 2-3.
#
# WHICH GVMAGs GO IN THE REFERENCE
# --------------------------------
# The scientific question is whether the giant viruses that bloomed in the INCUBATION
# experiments are present in the environment, so the 19 incubation-derived GVMAGs are the
# targets of interest. The 7 environmental GVMAGs are included as well, by default, for one
# reason: they are INTERNAL POSITIVE CONTROLS. Each was assembled from an environmental
# sample, so each MUST be recoverable by mapping that sample's reads back to it. If the
# incubation GVMAGs come out at near-zero coverage, the only way to know that is a real
# biological result rather than a broken mapping command is to see the environmental
# GVMAGs come out high in the same run. Including them costs nothing - the whole reference
# is 14 MB and CoverM maps the reads once.
#
# Set GVMAG_SET=incubation to build the 19 incubation GVMAGs only (the literal reading of
# the request); GVMAG_SET=all (the default) adds the 7 controls. Either way the manifest
# labels every genome with its source set, so downstream analysis can split them.
#
# OUTPUT
# ------
#   $HOME/scr10/data/kbrevis/abundance/gvmag_ref/<gvmag>.fna      mapping reference
#   $HOME/scr10/data/kbrevis/abundance/gvmag_reference_manifest.tsv
#
# CoverM's `genome` mode takes each genome's NAME from its FILENAME, so the files are named
# exactly as the GVMAG is named in output/votus/gvmags_votus.tsv
# (e.g. SRR28004154__SRR28004154_bin_22.fna). The manifest is what joins that name back to
# its vOTU, its source sample and its source set.
#
# Tiny (26 file copies, 14 MB), so it is allowed to run on dev (AGENTS.md section 1).
#
# Run:
#   bash scripts/14_gvmag_reference.sh                 # all 26 (default)
#   GVMAG_SET=incubation bash scripts/14_gvmag_reference.sh   # 19 incubation only
set -euo pipefail

REPO="/sciclone/home/iadanforth/K_brevis_viromes"
VOTUS="$REPO/output/votus/gvmags_votus.tsv"
GVMAG_DIR="$HOME/scr10/data/kbrevis/votus/gvmags"
OUT="$HOME/scr10/data/kbrevis/abundance/gvmag_ref"
MANIFEST="$HOME/scr10/data/kbrevis/abundance/gvmag_reference_manifest.tsv"
GVMAG_SET="${GVMAG_SET:-all}"

case "$GVMAG_SET" in
    all|incubation) ;;
    *) echo "GVMAG_SET must be 'all' or 'incubation', not '$GVMAG_SET'" >&2; exit 1 ;;
esac

# The paper's two sample sets. Incubation = the CH1 (May 2021) and CH2 (Aug 2021) bloom
# incubation experiments; environmental = the EH25 Gulf of Mexico monthly series plus two
# coastal Florida samples. Taken from metadata/SraRunTable.csv.
INCUBATION="SRR28004137 SRR28004138 SRR28004139 SRR28004140 SRR28004141 SRR28004142 SRR28004143 SRR28004144 SRR28004154"

[[ -s "$VOTUS" ]] || { echo "Missing $VOTUS (run step 09 first)" >&2; exit 1; }
[[ -d "$GVMAG_DIR" ]] || { echo "Missing $GVMAG_DIR (run step 08 first)" >&2; exit 1; }

mkdir -p "$OUT"
rm -f "$OUT"/*.fna
printf 'genome\tgvmag\tvotu\tvotu_rank\tsource_sample\tsource_set\tlength_bp\tcheckv_completeness\tfasta\n' > "$MANIFEST"

n=0; skipped=0
while IFS=$'\t' read -r votu gvmag sample isrep rep length comp ani af nmem; do
    [[ -n "${votu:-}" ]] || continue
    set_name=environmental
    for i in $INCUBATION; do [[ "$sample" == "$i" ]] && set_name=incubation; done
    if [[ "$GVMAG_SET" == "incubation" && "$set_name" != "incubation" ]]; then
        skipped=$((skipped + 1)); continue
    fi
    src="$GVMAG_DIR/${gvmag}.fna"
    [[ -s "$src" ]] || { echo "missing GVMAG FASTA: $src" >&2; exit 1; }
    cp "$src" "$OUT/${gvmag}.fna"
    rank=member; [[ "$isrep" == "yes" ]] && rank=representative
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$gvmag" "$gvmag" "$votu" "$rank" "$sample" "$set_name" "$length" "$comp" "$OUT/${gvmag}.fna" >> "$MANIFEST"
    n=$((n + 1))
done < <(tail -n +2 "$VOTUS")

echo "[ref] GVMAG_SET=$GVMAG_SET: wrote $n GVMAG reference genomes to $OUT (skipped $skipped)"
echo "[ref] manifest: $MANIFEST"
awk -F'\t' 'NR>1{c[$6]++; v[$3]=1} END{for(k in c) printf "  %-14s %d GVMAGs\n", k, c[k]; n=0; for(k in v) n++; printf "  %-14s %d\n", "distinct vOTUs", n}' "$MANIFEST"
echo "[ref] total reference size: $(du -shc "$OUT"/*.fna | tail -1 | cut -f1)"
