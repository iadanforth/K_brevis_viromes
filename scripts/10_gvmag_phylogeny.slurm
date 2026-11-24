#!/bin/bash
#SBATCH --job-name=10_gvmag_phylogeny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/10_gvmag_phylogeny_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/10_gvmag_phylogeny_%j.err

# Build a trimmed, concatenated giant virus marker-gene alignment from the step 09 vOTU
# representatives plus GVMAGs V2 reference genomes, following Booker et al. 2026
# (Supplemental Methods, "GVMAG comparisons"):
#   1. predict proteins for each vOTU representative       (pyrodigal-gv, Prodigal headers)
#   2. choose reference genomes from GVMAGs V2             (select_refs.py)
#   3. find, extract, concatenate and align the 7 default
#      giant virus marker genes                            (ncldv_markersearch.py -c)
#   4. trim the concatenated alignment                      (trimAl 1.5.1 -automated1)
# The output *.trimmed.aln is the input for tree inference (IQ-TREE; see
# README step 10).
#
# Needs: step 09 finished, and 07_gvmag_setup.slurm run with GVMAGS_V2_FULL=1 so the
# 1.18 GB reference bundle is present.
#
# Submit (16+ cores REQUIRED: ncldv_markersearch hardcodes `hmmsearch --cpu 16`):
#   ( module purge; module load slurm/kuro;   sbatch -p serial scripts/10_gvmag_phylogeny.slurm )
#   ( module purge; module load slurm/astral; sbatch -p batch --exclude=as01 scripts/10_gvmag_phylogeny.slurm )
#
# Options (all optional, after the script name):
#   --tier mcl|skani     reference tier: 2,606 (default) or 8,507 representatives
#   --orders LIST        'auto' (default: the orders of our own GVMAGs), 'all', or e.g.
#                        Imitervirales,Algavirales
#   --max-refs N         cap in-order references (useful to keep the alignment small;
#                        0 = no cap, the default)
#   --aligner clustalo|muscle   per-marker aligner (default clustalo, as in the paper)
#   --no-refs            our vOTUs only, no reference genomes

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:?Submit with sbatch}"
REPO="/sciclone/home/iadanforth/K_brevis_viromes"
DATA="$HOME/scr10/data/kbrevis"
VOTU_DIR="$DATA/votus"
GVMAG_ROOT="$DATA/gvmags"
OUT="$DATA/phylogeny"
NMS_DIR="$HOME/data10/programs/ncldv_markersearch"
GVMAGS_V2_DIR="$HOME/data10/databases/gvmagsV2"
PROJECT="kbrevis_gvmags"

TIER=mcl; ORDERS=auto; MAX_REFS=0; ALIGNER=clustalo; USE_REFS=1
while (($#)); do
    case "$1" in
        --tier) TIER="$2"; shift 2 ;;
        --orders) ORDERS="$2"; shift 2 ;;
        --max-refs) MAX_REFS="$2"; shift 2 ;;
        --aligner) ALIGNER="$2"; shift 2 ;;
        --no-refs) USE_REFS=0; shift ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

if ((THREADS < 16)); then
    echo "WARNING: ncldv_markersearch hardcodes 'hmmsearch --cpu 16'. With only $THREADS" >&2
    echo "cores allocated, hmmsearch will oversubscribe the node. Submit with -c 16 or more." >&2
fi

REPS="$VOTU_DIR/gvmags_representatives.txt"
[[ -s "$REPS" ]] || { echo "No vOTU representatives ($REPS); run step 09 first." >&2; exit 1; }
[[ -d "$NMS_DIR/hmm" ]] || {
    echo "ncldv_markersearch not installed ($NMS_DIR); run 07_gvmag_setup.slurm." >&2; exit 1; }

mkdir -p "$OUT"/{proteins,work}

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate phylo-env

# ---- 1. Proteins for each vOTU representative ---------------------------------
# ncldv_markersearch requires Prodigal-style protein IDs (contigname_1, contigname_2, ...)
# because it joins split marker genes that sit near each other on the same contig; giant
# virus genes are often fragmented across a MAG, so this is essential, not cosmetic.
# It takes the taxon name from the .faa file name, so names must contain no dots.
rm -f "$OUT"/proteins/*.faa
n_q=0
while read -r fna; do
    [[ -s "$fna" ]] || continue
    base=$(basename "$fna" .fna)
    # name the taxon by its vOTU, so tree tips are readable: vOTU003__SRR28004141_bin_2
    # $5 is is_representative: a sample_set column was inserted at $4 on 2026-11-24.
    votu=$(awk -F'\t' -v g="$base" '$2==g && $5=="yes" {print $1}' \
             "$VOTU_DIR/gvmags_votus.tsv" | head -1)
    name=$(echo "${votu:-noVOTU}__${base}" | tr -c 'A-Za-z0-9_-' '_')
    pyrodigal-gv -i "$fna" -a "$OUT/proteins/${name}.faa" -p meta -j "$THREADS" -o /dev/null
    n_q=$((n_q + 1))
done < "$REPS"
echo "[phylo] predicted proteins for $n_q vOTU representatives"

# ---- 2. Reference genomes from GVMAGs V2 --------------------------------------
if ((USE_REFS)); then
    REF_FAA_DIR=$(find "$GVMAGS_V2_DIR" -maxdepth 3 -type d -name 'faa*' | head -1)
    if [[ -z "$REF_FAA_DIR" ]]; then
        echo "The GVMAGs V2 protein bundle is not downloaded." >&2
        echo "Get it with: ( module purge; module load slurm/astral; sbatch -p batch \\" >&2
        echo "  --exclude=as01 --export=ALL,GVMAGS_V2_FULL=1 scripts/07_gvmag_setup.slurm )" >&2
        exit 1
    fi
    python -I "$REPO/scripts/select_refs.py" \
        --metadata "$GVMAGS_V2_DIR/GVDB_stats.tsv" \
        --faa-dir "$REF_FAA_DIR" \
        --out-dir "$OUT/proteins" \
        --tier "$TIER" --orders "$ORDERS" --max-refs "$MAX_REFS" \
        --gvclass-summaries "$GVMAG_ROOT"/*/*_gvmag_summary.tsv
else
    echo "[phylo] --no-refs: building the alignment from our vOTUs alone"
fi

N_TAXA=$(ls "$OUT"/proteins/*.faa | wc -l)
echo "[phylo] $N_TAXA taxa in the alignment ($n_q ours + $((N_TAXA - n_q)) references)"

# Reference .faa files come from the GVMAGs V2 bundle and were not necessarily produced by
# Prodigal. If their headers are not contigname_N, markersearch cannot join split genes for
# those taxa (it still finds single-copy hits, so the run succeeds but fragmented reference
# genes may be missed). Report rather than fail.
BAD=$(for f in "$OUT"/proteins/*.faa; do
          head -1 "$f" | grep -qE '^>[^ ]+_[0-9]+' || basename "$f"; done | wc -l)
((BAD == 0)) && echo "[phylo] all taxa have Prodigal-style protein IDs" \
             || echo "[phylo] NOTE: $BAD of $N_TAXA .faa files lack Prodigal-style IDs; split marker genes in those genomes may not be joined" >&2

# ---- 3. ncldv_markersearch: find, extract, concatenate, align -----------------
# Two quirks drive the directory juggling below:
#  * the script hardcodes a RELATIVE HMM path ("hmm/gvogs.common.hmm"), so it must run from
#    a directory containing `hmm/` - we symlink the installed copy;
#  * it writes all outputs into the current directory, so that same directory is the output.
# Default marker set is exactly the paper's seven: VLTF3 (GVOGm0890), A32 packaging ATPase
# (0760), DNA topoisomerase II (0461), TFIIB (0172), DNA polymerase B (0054), RNA polymerase
# alpha (0023) and SNF2-like helicase (0013). The major capsid protein (0003) and RNAP beta
# (0022) are available but NOT in the default set, as in the paper.
cd "$OUT/work"
ln -sfn "$NMS_DIR/hmm" hmm
echo "[phylo] running ncldv_markersearch (aligner: $ALIGNER); this is the slow step"
python -I "$NMS_DIR/ncldv_markersearch.py" \
    -i "$OUT/proteins" \
    -n "$PROJECT" \
    -t "$THREADS" \
    -g "$ALIGNER" \
    -c

CONCAT="$OUT/work/${PROJECT}.concat.aln"
[[ -s "$CONCAT" ]] || { echo "No concatenated alignment was produced ($CONCAT)." >&2; exit 1; }
echo "[phylo] concatenated alignment: $(grep -c '^>' "$CONCAT") taxa"

# ---- 4. trimAl -automated1 ----------------------------------------------------
# -automated1 chooses between gap-based and similarity-based trimming using a heuristic
# tuned for maximum-likelihood phylogenetics; it is what the paper used. The HTML report
# shows exactly which columns were dropped, which is worth checking before tree building.
TRIMMED="$OUT/${PROJECT}.trimmed.aln"
trimal -in "$CONCAT" -out "$TRIMMED" -automated1 -fasta \
       -htmlout "$OUT/${PROJECT}.trimal_report.html"
trimal -in "$TRIMMED" -sident > "$OUT/${PROJECT}.trimmed.identity.txt" 2>&1 || true

# ---- summary ------------------------------------------------------------------
cp "$OUT/work/${PROJECT}.table.tsv" "$OUT/" 2>/dev/null || true
cp "$OUT/work/${PROJECT}.full_output.txt" "$OUT/" 2>/dev/null || true
aln_len () {   # length of the first sequence in a FASTA alignment
    python -I -c "
import sys
seq = []
for line in open(sys.argv[1]):
    if line.startswith('>'):
        if seq: break
    else:
        seq.append(line.strip())
print(len(''.join(seq)))" "$1"
}
{
  echo "project:        $PROJECT"
  echo "date:           $(date)"
  echo "our vOTUs:      $n_q"
  echo "total taxa:     $(grep -c '^>' "$TRIMMED")"
  echo "reference tier: $([[ $USE_REFS == 1 ]] && echo "$TIER, orders=$ORDERS, max=$MAX_REFS" || echo none)"
  echo "aligner:        $ALIGNER"
  echo "columns before trimming: $(aln_len "$CONCAT")"
  echo "columns after trimming:  $(aln_len "$TRIMMED")"
} | tee "$OUT/${PROJECT}.phylogeny_summary.txt"

echo
echo "[phylo] done. Trimmed alignment for tree building:"
echo "  $TRIMMED"
echo "Per-marker occupancy per taxon: $OUT/${PROJECT}.table.tsv"
echo "trimAl report (which columns were dropped): $OUT/${PROJECT}.trimal_report.html"
echo "Next: infer a tree with IQ-TREE (see README step 10)."
