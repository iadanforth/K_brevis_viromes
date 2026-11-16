#!/bin/bash
#SBATCH --job-name=11_pre_tree_check2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=30:00
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/11_pre_tree_check2_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/11_pre_tree_check2_%j.err

# Second diagnostic on the step 10 alignment, separating the TWO different reasons a
# column can be empty for a taxon. The first check (11_pre_tree_check.slurm) counted them
# together and reported that all 1,564 taxa are >=70% "padded", which is true but useless:
#
#   'X'  = a marker gene that was NOT FOUND in that genome. This is missing DATA, and a
#          taxon made mostly of X cannot be placed in the tree.
#   '-'  = an alignment gap. This says nothing about the genome; it means Clustal Omega
#          inserted a gap because some OTHER sequence has residues there. At 14.5% average
#          pairwise identity across 1,564 divergent genomes, the alignment is mostly gaps
#          for everyone, including complete genomes.
#
# Only the X fraction should influence whether a taxon is dropped. This also reports how
# many columns would survive stricter gap trimming, to judge whether trimAl -automated1
# (which removed only 19% of columns) was lenient enough to be worth revisiting.
#
# Submit:
#   ( module purge; module load slurm/kuro; sbatch -p serial scripts/11_pre_tree_check2.slurm )

set -euo pipefail

ALN="$HOME/scr10/data/kbrevis/phylogeny/kbrevis_gvmags.trimmed.aln"
[[ -s "$ALN" ]] || { echo "No trimmed alignment ($ALN); run step 10 first." >&2; exit 1; }

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate tree-env

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
print(f"taxa: {len(seqs)}  columns: {L}\n")

def frac(s, chars):
    return 100.0 * sum(s.count(c) for c in chars) / len(s)

xs   = {k: frac(v, "Xx") for k, v in seqs.items()}
gaps = {k: frac(v, "-")  for k, v in seqs.items()}
res  = {k: 100.0 - xs[k] - gaps[k] for k in seqs}

print("=== MISSING DATA (X only) - this is what matters for dropping taxa ===")
b = collections.Counter()
for v in xs.values():
    b[min(int(v // 10) * 10, 90)] += 1
print("X% histogram (bucket: taxa):", dict(sorted(b.items())))
for t in (50, 70, 90):
    print(f"  >={t}% X: {sum(1 for v in xs.values() if v >= t)} taxa")

print("\n=== ALIGNMENT GAPS (- only) - affects everyone, not a quality signal ===")
g = sorted(gaps.values())
print(f"  median {g[len(g)//2]:.1f}%   min {g[0]:.1f}%   max {g[-1]:.1f}%")

print("\n=== REAL RESIDUES per taxon ===")
r = sorted(res.values())
print(f"  median {r[len(r)//2]:.1f}%   min {r[0]:.1f}%   max {r[-1]:.1f}%")
print(f"  that is a median of ~{int(r[len(r)//2]/100*L)} amino acids of actual data")

print("\n=== our vOTUs: X% (missing markers) vs residues ===")
print(f"  {'vOTU':10s} {'X%':>7s} {'gap%':>7s} {'residue%':>9s} {'aa':>6s}")
for k in sorted(seqs):
    if k.startswith("vOTU"):
        print(f"  {k.split('.')[0]:10s} {xs[k]:7.1f} {gaps[k]:7.1f} {res[k]:9.1f} "
              f"{int(res[k]/100*L):6d}")

print("\n=== column occupancy: how many columns would survive stricter gap trimming ===")
names = list(seqs); n = len(names)
cols = [0] * L
for s in seqs.values():
    for i, c in enumerate(s):
        if c not in "-Xx":
            cols[i] += 1
for keep in (0.05, 0.10, 0.20, 0.50):
    surv = sum(1 for c in cols if c >= keep * n)
    print(f"  columns with data in >={int(keep*100):3d}% of taxa: {surv:6d} "
          f"({100.0*surv/L:5.1f}% of {L})")
PY
