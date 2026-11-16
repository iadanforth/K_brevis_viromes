#!/bin/bash
#SBATCH --job-name=11_pre_tree_check
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=30:00
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/11_pre_tree_check_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/11_pre_tree_check_%j.err

# Cheap diagnostic run BEFORE committing a multi-day IQ-TREE job (step 11).
#
# Step 10 produced a 1,564-taxon x 37,739-column amino-acid alignment. IQ-TREE's memory
# use scales with taxa x site-patterns x states x rate-categories, which at this scale can
# land anywhere between tens and hundreds of GB - and if the request is too small the job
# dies, while if it is too large it waits in the queue for nothing and blocks other users.
# IQ-TREE prints the RAM it needs within seconds of reading an alignment, so this job asks
# it, rather than guessing. It also reports how much of the alignment is X/gap padding,
# which is what decides whether any taxa should be dropped before building the tree.
#
# Submit:
#   ( module purge; module load slurm/kuro; sbatch -p serial scripts/11_pre_tree_check.slurm )

set -euo pipefail

ALN="$HOME/scr10/data/kbrevis/phylogeny/kbrevis_gvmags.trimmed.aln"
WORK="$HOME/scr10/data/kbrevis/tree_probe"
[[ -s "$ALN" ]] || { echo "No trimmed alignment ($ALN); run step 10 first." >&2; exit 1; }
mkdir -p "$WORK"

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate tree-env

echo "===== alignment composition ====="
python -I - "$ALN" <<'PY'
import sys, collections
aln = sys.argv[1]; name = None; seqs = {}
for line in open(aln):
    line = line.rstrip("\n")
    if line.startswith(">"):
        name = line[1:].split()[0]; seqs[name] = []
    elif name:
        seqs[name].append(line.strip())
seqs = {k: "".join(v) for k, v in seqs.items()}
L = len(next(iter(seqs.values())))
print(f"taxa: {len(seqs)}  columns: {L}")
print("non-rectangular sequences:", sum(1 for v in seqs.values() if len(v) != L))

def pad(s):
    return 100.0 * (s.count("X") + s.count("x") + s.count("-")) / len(s)

p = {k: pad(v) for k, v in seqs.items()}
buckets = collections.Counter()
for v in p.values():
    buckets[min(int(v // 10) * 10, 90)] += 1
print("padding/gap % histogram (lower bound of 10%-bucket: taxa):",
      dict(sorted(buckets.items())))
print(f">=70% padded: {sum(1 for v in p.values() if v >= 70)} taxa")
print(f">=90% padded: {sum(1 for v in p.values() if v >= 90)} taxa")

ours = {k: v for k, v in p.items() if k.startswith("vOTU")}
print(f"\nour vOTUs present in the alignment: {len(ours)}")
for k, v in sorted(ours.items(), key=lambda x: -x[1]):
    print(f"  {v:5.1f}% padded  {k}")
missing = [n for n in ("vOTU05", "vOTU19") if not any(k.startswith(n) for k in seqs)]
print("\nzero-marker vOTUs absent from the alignment:", missing or "none (both present!)")
PY

echo
echo "===== IQ-TREE memory requirement ====="
# -n 0 stops after initialisation: enough for IQ-TREE to report the RAM it needs and the
# number of distinct site patterns (the real driver of memory), without searching for a
# tree. A fixed model is used because ModelFinder would do far more work than needed here.
cd "$WORK"
iqtree2 -s "$ALN" -m LG+F+G4 -n 0 -T 4 --prefix probe --seed 1 --redo 2>&1 \
  | grep -iE 'RAM|memory|sequences with|columns|patterns|constant|ERROR|WARNING' | head -30 || true

echo
echo "--- RAM / pattern lines from probe.log:"
grep -iE 'RAM|pattern' probe.log 2>/dev/null | head -15 || true
