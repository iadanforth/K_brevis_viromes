#!/bin/bash
#SBATCH --job-name=11_gvmag_tree
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/11_gvmag_tree_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/11_gvmag_tree_%j.err

# Infer a maximum-likelihood phylogenomic tree from the step 10 trimmed, concatenated
# marker-gene alignment, following Booker et al. 2026 (Supplemental Methods):
#   1. IQ-TREE: model selection + ML tree + ultrafast bootstrap   (paper: -m MFP -bb 1000 -nt AUTO -wbt)
#   2. RAxML -f i: internode certainty (IC) and tree certainty (TC) scores
#
# Flag names below are IQ-TREE 2.x, verified against the installed 2.3.6 help text. The
# paper's `-bb`/`-nt`/`-wbt` are IQ-TREE 1.x spellings; the 2.x equivalents are
# `-B`/`-T`/`--boot-trees` (we use `--wbtl`, which also writes branch lengths because
# RAxML's IC/TC step wants them).
#
# Submit (needs a 72 h partition; kuro caps at 48 h, so add -t 2-00:00:00 there):
#   ( module purge; module load slurm/femto;  sbatch -p batch scripts/11_gvmag_tree.slurm )
#   ( module purge; module load slurm/astral; sbatch -p batch --exclude=as01 scripts/11_gvmag_tree.slurm )
#   ( module purge; module load slurm/kuro;   sbatch -p serial -t 2-00:00:00 scripts/11_gvmag_tree.slurm )
#
# Options:
#   --alignment PATH   alignment to use (default: step 10's kbrevis_gvmags.trimmed.aln)
#   --model STRING     IQ-TREE model (default MFP = ModelFinder Plus, as in the paper).
#                      Use e.g. LG+F+G4 to skip model selection on a large alignment.
#   --bootstraps N     ultrafast bootstrap replicates (default 1000, as in the paper)
#   --threads N|AUTO   default: the Slurm allocation (the paper used AUTO)
#   --alrt N           also run SH-aLRT with N replicates (default 0 = off)
#   --skip-ic          skip the RAxML IC/TC step
#   --mset LIST        restrict ModelFinder to these matrices (e.g. LG,WAG,JTT). Only
#                      meaningful with --model MFP/MF. Default: empty = all matrices.
#   --cmax N           max rate categories ModelFinder tries (IQ-TREE default 10)
#   --safe             IQ-TREE's -safe: slower but numerically stable likelihood scaling

set -euo pipefail

THREADS_ALLOC="${SLURM_CPUS_PER_TASK:?Submit with sbatch}"
REPO="/sciclone/home/iadanforth/K_brevis_viromes"
DATA="$HOME/scr10/data/kbrevis"
PHYLO_DIR="$DATA/phylogeny"
OUT="$DATA/tree"
PROJECT="kbrevis_gvmags"

ALN="$PHYLO_DIR/${PROJECT}.trimmed.aln"
MODEL=MFP; BOOTSTRAPS=1000; THREADS="$THREADS_ALLOC"; ALRT=0; DO_IC=1
MSET=""; CMAX=""; SAFE=0
while (($#)); do
    case "$1" in
        --alignment) ALN="$2"; shift 2 ;;
        --model) MODEL="$2"; shift 2 ;;
        --bootstraps) BOOTSTRAPS="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --alrt) ALRT="$2"; shift 2 ;;
        --skip-ic) DO_IC=0; shift ;;
        --mset) MSET="$2"; shift 2 ;;
        --cmax) CMAX="$2"; shift 2 ;;
        --safe) SAFE=1; shift ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

[[ -s "$ALN" ]] || { echo "No trimmed alignment ($ALN); run step 10 first." >&2; exit 1; }
mkdir -p "$OUT"

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate tree-env

# ---- Pre-flight: look at the alignment before spending days on it -------------
# IQ-TREE refuses to bootstrap fewer than 4 sequences ("It makes no sense to perform
# bootstrap with less than 4 sequences"), and it only says so AFTER model selection has
# run. Catch it up front. Found while smoke-testing on a 3-taxon alignment.
N_TAXA=$(grep -c '^>' "$ALN")
if ((N_TAXA < 4)) && ((BOOTSTRAPS > 0)); then
    echo "Only $N_TAXA taxa in $ALN. IQ-TREE cannot bootstrap fewer than 4 sequences." >&2
    echo "Either add reference genomes (step 10 without --no-refs) or pass --bootstraps 0." >&2
    exit 1
fi

# Taxa made up mostly of X padding (markers that were not found) are placed unreliably and
# slow everything down. This reports them so they can be removed and the step re-run; it
# does not drop them automatically, since that is a scientific decision.
python - "$ALN" <<'PY'
import sys
aln, name, seqs = sys.argv[1], None, {}
for line in open(aln):
    line = line.rstrip("\n")
    if line.startswith(">"):
        name = line[1:].split()[0]; seqs[name] = []
    elif name:
        seqs[name].append(line.strip())
seqs = {k: "".join(v) for k, v in seqs.items()}
if not seqs:
    sys.exit("alignment contains no sequences")
L = len(next(iter(seqs.values())))
bad = {k: v for k, v in seqs.items() if len(v) != L}
print(f"[tree] {len(seqs)} taxa, {L} alignment columns")
if bad:
    sys.exit(f"alignment is not rectangular: {len(bad)} sequences differ in length "
             f"(first: {next(iter(bad))})")
def pad(s):
    return 100.0 * (s.count("X") + s.count("x") + s.count("-")) / len(s)
poor = sorted(((pad(v), k) for k, v in seqs.items()), reverse=True)
over = [(p, k) for p, k in poor if p >= 70]
print(f"[tree] worst padding/gap fraction: {poor[0][0]:.1f}% ({poor[0][1]})")
if over:
    print(f"[tree] WARNING: {len(over)} taxa are >=70% X/gap and are placed unreliably:")
    for p, k in over[:10]:
        print(f"[tree]   {p:5.1f}%  {k}")
    print("[tree] Consider removing them (see phylogeny/*.table.tsv) and re-running step 10.")
PY

# ---- 1. IQ-TREE: ModelFinder + ML tree + ultrafast bootstrap ------------------
# -m MFP  : ModelFinder Plus - tests substitution models (for protein data: LG, WAG, JTT,
#           ... with +F/+G/+R variants), picks the best by BIC, then infers the tree.
# -B 1000 : ultrafast bootstrap (UFBoot), 1000 replicates. UFBoot support values are NOT
#           comparable to classical bootstrap: >=95 is roughly as strong as >=70 classical.
# --bnni  : optimise each UFBoot tree by NNI. Recommended by the IQ-TREE authors because it
#           reduces UFBoot's tendency to overestimate support when the model is misspecified
#           - a real risk with concatenated markers. (Not in the paper; a deliberate addition.)
# --wbtl  : write the bootstrap trees, with branch lengths, to <prefix>.ufboot. The paper's
#           -wbt did the same without branch lengths; RAxML's IC/TC step reads this file.
# -T      : explicit thread count rather than the paper's AUTO, because AUTO probes the
#           machine and can exceed the Slurm allocation.
PREFIX="$OUT/$PROJECT"
# --mem caps IQ-TREE's own memory use at 85% of the Slurm allocation, so it reports a
# shortfall rather than being OOM-killed. (On small alignments IQ-TREE prints harmless
# "Too low -mem, automatically increased to N MB" notes per model; ignore them.)
IQ_ARGS=(-s "$ALN" -m "$MODEL" --bnni
         -T "$THREADS" --prefix "$PREFIX" --seed 12345
         --mem "$(( ${SLURM_MEM_PER_NODE:-131072} * 85 / 100 / 1024 ))G")
((BOOTSTRAPS > 0)) && IQ_ARGS+=(-B "$BOOTSTRAPS" --wbtl)
((ALRT > 0)) && IQ_ARGS+=(--alrt "$ALRT")
# -mset / -cmax bound how many models ModelFinder tries. Measured on this alignment
# (job 111549): optimising the parameters of ONE model took 1,492 s, so an unrestricted
# ModelFinder over every protein matrix x rate-heterogeneity variant would spend days
# before the tree search even starts. Restricting to the standard matrices keeps the
# paper's "pick the best-fit model by BIC" approach affordable.
[[ -n "$MSET" ]] && IQ_ARGS+=(-mset "$MSET")
[[ -n "$CMAX" ]] && IQ_ARGS+=(-cmax "$CMAX")
# -safe switches IQ-TREE to scaled likelihood computation. IQ-TREE itself asked for it on
# this alignment ("WARNING: Numerical underflow for lh-branch, use `-safe` option"), which
# is expected at 14.5% average pairwise identity across 1,564 taxa. It is slower, but an
# underflowing likelihood can mislead both the tree search and the support values.
((SAFE)) && IQ_ARGS+=(-safe)

# No --redo: if a previous run of this job hit the walltime, IQ-TREE finds
# <prefix>.ckp.gz and resumes from it rather than starting over. That is what makes this
# job safe to run on a partition whose time limit is shorter than the whole analysis -
# just resubmit the identical command and it continues. Pass --redo by hand only to
# deliberately discard a previous run.
if [[ -s "$PREFIX.ckp.gz" ]]; then
    echo "[tree] found $PREFIX.ckp.gz - IQ-TREE will resume from the checkpoint"
fi
echo "[tree] iqtree2 ${IQ_ARGS[*]}"
iqtree2 "${IQ_ARGS[@]}"

TREE="$PREFIX.treefile"
[[ -s "$TREE" ]] || { echo "IQ-TREE produced no tree ($TREE)" >&2; exit 1; }

# ---- 2. RAxML -f i: internode certainty / tree certainty ---------------------
# -f i computes IC and TC scores (Salichos & Rokas 2013): how consistently each internal
# branch is supported across the bootstrap trees, which exposes conflicting signal that a
# single bootstrap number can hide.
#
# DELIBERATE DEVIATION FROM THE PAPER: it reports `-f i -m GTRCAT`, but GTRCAT is RAxML's
# GTR *nucleotide* model and this alignment is amino acids. RAxML's protein models are
# PROTCAT<matrix>/PROTGAMMA<matrix>, so PROTCATLG is used here (LG is the modern general
# protein matrix). The IC/TC calculation compares tree topologies rather than optimising a
# likelihood, so the matrix choice has little effect - but the model must at least match the
# data type or RAxML will reject the alignment.
if ((DO_IC)); then
    UFBOOT="$PREFIX.ufboot"
    if [[ -s "$UFBOOT" ]]; then
        IC_DIR="$OUT/ic"; rm -rf "$IC_DIR"; mkdir -p "$IC_DIR"
        # RAxML writes RAxML_* files into -w, which must be an absolute path.
        raxmlHPC-PTHREADS-SSE3 -T "$THREADS" -f i \
            -t "$TREE" -z "$UFBOOT" -s "$ALN" \
            -m PROTCATLG -n ic -w "$IC_DIR" -C \
            || echo "[tree] WARNING: RAxML IC/TC step failed; the ML tree above is unaffected" >&2
        cp "$IC_DIR"/RAxML_info.ic "$OUT/${PROJECT}.IC_TC_info.txt" 2>/dev/null || true
    else
        echo "[tree] no .ufboot file, so the IC/TC step was skipped" >&2
    fi
fi

# ---- Summary -----------------------------------------------------------------
{
  echo "project:            $PROJECT"
  echo "date:               $(date)"
  echo "alignment:          $ALN"
  echo "model requested:    $MODEL"
  echo "model chosen:       $(grep -m1 'Best-fit model' "$PREFIX.log" 2>/dev/null | sed 's/.*: //' || echo 'see log')"
  echo "log-likelihood:     $(grep -m1 'BEST SCORE FOUND' "$PREFIX.log" 2>/dev/null | sed 's/.*: //' || echo 'see log')"
  echo "bootstrap reps:     $( ((BOOTSTRAPS>0)) && echo "$BOOTSTRAPS (ultrafast UFBoot, --bnni)" || echo "none (--bootstraps 0)" )"
  echo "threads:            $THREADS"
  echo "taxa in tree:       $(tr ',' '\n' < "$TREE" | grep -c ':' || true)"
  echo
  if ((BOOTSTRAPS > 0)); then
  echo "UFBoot support distribution (branch support values):"
  python - "$TREE" <<'PY'
import re, sys
vals = [float(x) for x in re.findall(r"\)(\d+(?:\.\d+)?):", open(sys.argv[1]).read())]
if not vals:
    print("  (no support values found in the tree file)")
else:
    vals.sort()
    n = len(vals)
    strong = sum(1 for v in vals if v >= 95)
    print(f"  {n} internal branches; median {vals[n//2]:.0f}")
    print(f"  >=95 (well supported): {strong} ({100*strong/n:.0f}%)")
    print(f"  <70  (poorly supported): {sum(1 for v in vals if v < 70)}")
PY
  fi
} | tee "$OUT/${PROJECT}.tree_summary.txt"

echo
echo "[tree] done. Key outputs in $OUT:"
echo "  ${PROJECT}.treefile          the ML tree with UFBoot support (Newick)"
echo "  ${PROJECT}.iqtree            human-readable report: model, likelihood, parameters"
echo "  ${PROJECT}.log               full IQ-TREE log"
[[ -s "$PREFIX.ufboot" ]] && echo "  ${PROJECT}.ufboot            the bootstrap trees"
[[ -s "$OUT/${PROJECT}.IC_TC_info.txt" ]] && echo "  ${PROJECT}.IC_TC_info.txt    RAxML internode/tree certainty scores"
echo "  ${PROJECT}.tree_summary.txt  this summary"
echo
echo "The tree is UNROOTED. Booker et al. midpoint-rooted theirs; do that in iTOL,"
echo "FigTree or PearTree, or with ape::midpoint in R, before interpreting it."
