#!/bin/bash
# Build the mapping reference for step 14: one FASTA per vOTU, named by vOTU ID.
#
# CoverM's `genome` mode takes the genome name from the FILENAME, so the step 09
# representatives (named <sample>__<bin>.fna) are copied to vOTU01.fna ... vOTU19.fna.
# Without this the abundance table would be keyed by bin name, which hides which vOTU is
# which and makes the figure unreadable.
#
# A vOTU is represented by ONE genome - its representative - not by all its members. If
# every member were included, reads from a shared virus would be split across near-identical
# genomes and each would look less abundant than it is.
#
# Tiny (19 files, 11 MB), so it runs on dev. Writes a manifest recording which bin each
# vOTU ID maps to and whether that vOTU was discovered in an incubation or an environmental
# sample.
set -euo pipefail

REPO="/sciclone/home/iadanforth/K_brevis_viromes"
VOTUS="$REPO/output/votus/gvmags_votus.tsv"
OUT="$HOME/scr10/data/kbrevis/abundance/votu_ref"

# The paper's two sample sets. The incubation experiment (CH1 May 2021, CH2 Aug 2021) and
# the Gulf of Mexico / coastal Florida environmental series.
INCUBATION="SRR28004137 SRR28004138 SRR28004139 SRR28004140 SRR28004141 SRR28004142 SRR28004143 SRR28004144 SRR28004154"

mkdir -p "$OUT"
rm -f "$OUT"/*.fna
MANIFEST="$HOME/scr10/data/kbrevis/abundance/votu_reference_manifest.tsv"
printf 'votu\tsource_sample\tsource_set\trepresentative_bin\tfasta\n' > "$MANIFEST"

n=0
while IFS=$'\t' read -r votu gvmag sample isrep rep length comp ani af nmem; do
    [[ "$isrep" == "yes" ]] || continue
    src="$HOME/scr10/data/kbrevis/votus/gvmags/${gvmag}.fna"
    [[ -s "$src" ]] || { echo "missing representative FASTA: $src" >&2; exit 1; }
    cp "$src" "$OUT/${votu}.fna"
    set_name=environmental
    for i in $INCUBATION; do [[ "$sample" == "$i" ]] && set_name=incubation; done
    printf '%s\t%s\t%s\t%s\t%s\n' "$votu" "$sample" "$set_name" "$gvmag" "$OUT/${votu}.fna" >> "$MANIFEST"
    n=$((n + 1))
done < <(tail -n +2 "$VOTUS")

echo "[ref] wrote $n vOTU reference genomes to $OUT"
echo "[ref] manifest: $MANIFEST"
awk -F'\t' 'NR>1{c[$3]++} END{for(k in c) printf "  %-14s %d vOTUs\n", k, c[k]}' "$MANIFEST"
echo "[ref] total reference size: $(du -shc "$OUT"/*.fna | tail -1 | cut -f1)"
