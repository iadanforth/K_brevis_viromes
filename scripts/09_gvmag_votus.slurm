#!/bin/bash
#SBATCH --job-name=09_gvmag_votus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/09_gvmag_votus_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/09_gvmag_votus_%j.err

# Group the GVMAGs from ALL samples into viral OTUs ("viral species"), following
# Booker et al. 2026 (Supplemental Methods, "GVMAG comparisons"):
#   1. pool every sample's GVMAGs from step 8
#   2. all-vs-all ANI                               (fastANI, via ANIclustermap's engine)
#   3. cluster into vOTUs at >=95% ANI / >=85% AF    (MIUViG standard; votu_cluster.py)
#   4. clustered ANI heatmap                         (ANIclustermap 1.1.0, the paper's tool)
#   5. optional: nearest relative in GVMAGs V2       (Vasquez et al. 2025, 8,508 genomes)
#
# This runs ONCE for the whole project, after 08_gvmags.sh has finished for every sample.
# Submit (pick the cluster with ~/hpc/best_cluster.sh -c 16 -m 64 -H 2):
#   ( module purge; module load slurm/astral; sbatch -p batch --exclude=as01 scripts/09_gvmag_votus.slurm )
#   ...add --with-refs to also compare against the GVMAGs V2 reference genomes:
#   ( module purge; module load slurm/astral; sbatch -p batch --exclude=as01 scripts/09_gvmag_votus.slurm --with-refs )

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:?Submit with sbatch}"
REPO="/sciclone/home/iadanforth/K_brevis_viromes"
DATA="$HOME/scr10/data/kbrevis"
GVMAG_ROOT="$DATA/gvmags"
OUT="$DATA/votus"
GVMAGS_V2_DIR="$HOME/data10/databases/gvmagsV2"
MIN_ANI="${MIN_ANI:-95}"
MIN_AF="${MIN_AF:-85}"

WITH_REFS=0
[[ "${1:-}" == "--with-refs" ]] && WITH_REFS=1

[[ -d "$GVMAG_ROOT" ]] || { echo "No GVMAGs yet ($GVMAG_ROOT); run step 08 first." >&2; exit 1; }
mkdir -p "$OUT"/{gvmags,ani,heatmap,refs}

module load miniforge3
eval "$(conda shell.bash hook)"

# ---- 1. Pool every sample's GVMAGs, renamed <sample>__<bin> --------------------
# Bin names repeat between samples (every sample has a bin_1), so the sample name is
# prefixed to keep each GVMAG traceable to the sample it came from.
rm -f "$OUT"/gvmags/*.fna
shopt -s nullglob
n_samples=0
for sample_dir in "$GVMAG_ROOT"/*/; do
    sample=$(basename "$sample_dir")
    # Only pool directories that step 08 actually finished. Without this, ANY directory
    # under gvmags/ holding a gvmags/*.fna is pooled as an extra genome - which is how the
    # 08 smoke-test directory (_tailtest, a bacterial bin used to exercise CheckV/geNomad)
    # would have entered the vOTU analysis as a 27th "GVMAG" with no summary row.
    if [[ ! -f "${sample_dir}.done" ]]; then
        echo "[votus] skipping $sample: no .done marker from step 08"
        continue
    fi
    found=0
    for g in "$sample_dir"gvmags/*.fna; do
        cp "$g" "$OUT/gvmags/${sample}__$(basename "$g")"
        found=1
    done
    ((found)) && n_samples=$((n_samples + 1))
done
GVMAGS=("$OUT"/gvmags/*.fna)
echo "[votus] pooled ${#GVMAGS[@]} GVMAGs from $n_samples samples"
if ((${#GVMAGS[@]} == 0)); then
    echo "[votus] no GVMAGs were found in any sample; nothing to cluster."
    exit 0
fi

# One combined table of every sample's per-GVMAG summary, with the same <sample>__<bin>
# names, so completeness can be used to pick each vOTU's representative.
#
# A `sample_set` column is inserted after `sample` (added 2026-11-24, prompt 68) saying
# whether that GVMAG came from an incubation experiment or an environmental sample. The
# whole analysis turns on that split - the question is whether the viruses that bloomed in
# the incubations are present in the environment - and without it every consumer of this
# table has to re-derive it from its own copy of the accession list. Both existing readers
# (votu_cluster.py and plot_incubation_votu_timeseries.R) use csv.DictReader / read_tsv and
# address columns by name, so inserting a column does not disturb them.
#
# These 9 accessions are the CH1 (May 2021) and CH2 (Aug 2021) incubation runs; everything
# else in the study is environmental. This list must stay identical to INCUBATION_SAMPLES
# in pipeline.conf and to the one in scripts/14_gvmag_reference.sh; the count is checked
# below so a typo is caught rather than silently reclassifying a sample.
INCUBATION="SRR28004137 SRR28004138 SRR28004139 SRR28004140 SRR28004141 SRR28004142 SRR28004143 SRR28004144 SRR28004154"
N_INC=$(wc -w <<< "$INCUBATION")
if (( N_INC != 9 )); then
    echo "[votus] INCUBATION holds $N_INC accessions, expected 9 - refusing to label" >&2
    exit 1
fi

COMBINED="$OUT/all_gvmag_summaries.tsv"
: > "$COMBINED"
for s in "$GVMAG_ROOT"/*/*_gvmag_summary.tsv; do
    sample=$(basename "$(dirname "$s")")
    set_name=environmental
    for i in $INCUBATION; do [[ "$sample" == "$i" ]] && set_name=incubation; done
    awk -F'\t' -v s="$sample" -v set="$set_name" \
        -v first="$([[ -s $COMBINED ]] && echo 0 || echo 1)" '
        function emit(col2, col1,    i) {
            printf "%s\t%s", col1, col2
            for (i = 2; i <= NF; i++) printf "\t%s", $i
            printf "\n"
        }
        NR==1 { if (first) emit("sample_set", $1); next }
        { $2 = s "__" $2; emit(set, $1) }' "$s" >> "$COMBINED"
done
# A sample with no GVMAG writes no *_gvmag_summary.tsv at all (step 08 exits early), so
# the glob above simply matches fewer files - that is expected and harmless. The counts
# should still agree, because such a sample contributes no genomes either; if they do not,
# say so, but do not fail: votu_cluster.py ranks a GVMAG with no summary row last
# (completeness defaults to -1) and reports its completeness as NA.
N_ROWS=0
[[ -s "$COMBINED" ]] && N_ROWS=$(( $(wc -l < "$COMBINED") - 1 ))
echo "[votus] combined summary rows: $N_ROWS (pooled GVMAGs: ${#GVMAGS[@]})"
if (( N_ROWS != ${#GVMAGS[@]} )); then
    echo "[votus] WARNING: ${#GVMAGS[@]} GVMAGs pooled but $N_ROWS summary rows found." >&2
    echo "[votus] Clustering still runs; the unmatched genomes get completeness=NA." >&2
fi

# ---- 2. All-vs-all ANI (fastANI) ----------------------------------------------
# fastANI compares whole genomes without aligning them base by base, and handles
# multi-contig genomes, which matters because a GVMAG is usually several contigs.
conda activate aniclustermap-env
printf '%s\n' "${GVMAGS[@]}" > "$OUT/ani/gvmag_list.txt"
fastANI --ql "$OUT/ani/gvmag_list.txt" --rl "$OUT/ani/gvmag_list.txt" \
    -t "$THREADS" -o "$OUT/ani/gvmags_all_vs_all.tsv"
echo "[votus] fastANI pairs: $(wc -l < "$OUT/ani/gvmags_all_vs_all.tsv")"

# ---- 4. Clustered ANI heatmap (the paper's Supplemental Figure 3) -------------
# Run before deactivating the env; failure here must not lose the vOTU results.
if ((${#GVMAGS[@]} >= 2)); then
    ANIclustermap -i "$OUT/gvmags" -o "$OUT/heatmap" \
        --fig_width 15 --fig_height 12 \
        || echo "[votus] WARNING: ANIclustermap failed; vOTU results below are unaffected" >&2
else
    echo "[votus] only one GVMAG: skipping the heatmap"
fi
conda deactivate

# ---- 3. Cluster into vOTUs ----------------------------------------------------
conda activate checkv-env   # used only for its Python 3.12
python -I "$REPO/scripts/votu_cluster.py" \
    --fastani "$OUT/ani/gvmags_all_vs_all.tsv" \
    --gvmag-dir "$OUT/gvmags" \
    --summary "$COMBINED" \
    --out-prefix "$OUT/gvmags" \
    --min-ani "$MIN_ANI" --min-af "$MIN_AF" | tee "$OUT/votu_summary.txt"
conda deactivate

# ---- 5. Optional: nearest relative among the GVMAGs V2 references -------------
# Each vOTU representative is compared to the 8,508 ANI-clustered reference genomes
# of Vasquez et al. 2025. This says how novel each vOTU is: a high ANI to a reference
# means it is a known giant virus, while no hit at all means fastANI found nothing
# above ~80% ANI, i.e. no close relative has been sequenced.
if ((WITH_REFS)); then
    REF_FNA_DIR=$(find "$GVMAGS_V2_DIR" -maxdepth 3 -type d -name 'fna*' | head -1)
    if [[ -z "$REF_FNA_DIR" ]]; then
        echo "[votus] --with-refs needs the full GVMAGs V2 bundle, which is not downloaded." >&2
        echo "[votus] Get it with: sbatch --export=ALL,GVMAGS_V2_FULL=1 ... scripts/07_gvmag_setup.slurm" >&2
    else
        conda activate aniclustermap-env
        find "$REF_FNA_DIR" -name '*.fna' > "$OUT/refs/ref_list.txt"
        echo "[votus] comparing vOTU representatives to $(wc -l < "$OUT/refs/ref_list.txt") references"
        fastANI --ql "$OUT/gvmags_representatives.txt" --rl "$OUT/refs/ref_list.txt" \
            -t "$THREADS" -o "$OUT/refs/votus_vs_gvmagsV2.tsv"
        conda deactivate
        # Best hit per representative, annotated from the GVMAGs V2 metadata table.
        conda activate checkv-env
        python -I "$REPO/scripts/votu_nearest_ref.py" \
            --fastani "$OUT/refs/votus_vs_gvmagsV2.tsv" \
            --votus "$OUT/gvmags_votus.tsv" \
            --metadata "$GVMAGS_V2_DIR/GVDB_stats.tsv" \
            --out "$OUT/votu_nearest_reference.tsv"
        conda deactivate
    fi
fi

echo "[votus] done. Results in $OUT:"
echo "  gvmags_votus.tsv            vOTU of every GVMAG"
echo "  gvmags_representatives.txt  one representative genome per vOTU"
echo "  votu_summary.txt            how many vOTUs were formed"
echo "  heatmap/                    clustered ANI heatmap"
((WITH_REFS)) && echo "  votu_nearest_reference.tsv  closest GVMAGs V2 relative per vOTU"
exit 0
