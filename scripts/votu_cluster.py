#!/usr/bin/env python3
"""Group GVMAGs from all samples into viral operational taxonomic units (vOTUs).

A vOTU is the metagenomics stand-in for a viral "species": a cluster of genomes
similar enough to be treated as the same virus. The community standard (MIUViG,
Roux et al. 2019 Nat Biotechnol) is >=95% average nucleotide identity (ANI) over
>=85% of the shorter genome, which is also the default in CheckV's `aniclust`.
Booker et al. 2026 say only that they grouped GVMAGs into vOTUs "based on their
average nucleotide identity" without giving a cut-off, so the MIUViG values are
used here and can be changed with --min-ani / --min-af.

Clustering is greedy and centroid-based, like `aniclust`: GVMAGs are ranked
best-first, and each unassigned GVMAG in turn becomes a vOTU representative and
absorbs every unassigned GVMAG that matches it. Ranking is by CheckV completeness
(if a summary table is given), then genome length, then name, so the most complete
and longest genome represents each vOTU. That mirrors the paper, where one vOTU
(e.g. vOTU5) was recovered from two samples and reported once.

Inputs:
  --fastani   all-vs-all fastANI output (query, ref, ANI, mapped_frags, total_frags)
  --gvmag-dir directory of <sample>__<bin>.fna, one file per GVMAG
  --summary   optional combined per-GVMAG summary table (needs columns `gvmag`
              and `checkv_completeness`) used only for ranking

Outputs (written next to --out-prefix):
  <prefix>_votus.tsv          one row per GVMAG: its vOTU, whether it represents it,
                              and its ANI / alignment fraction to the representative
  <prefix>_representatives.txt  the representative GVMAGs' .fna paths, one per line
"""
import argparse
import csv
import sys
from pathlib import Path


def genome_length(path):
    return sum(len(line.strip()) for line in open(path) if not line.startswith(">"))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fastani", required=True)
    p.add_argument("--gvmag-dir", required=True)
    p.add_argument("--summary")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--min-ani", type=float, default=95.0)
    p.add_argument("--min-af", type=float, default=85.0,
                   help="minimum aligned fraction of the SHORTER genome, in percent")
    a = p.parse_args()

    lengths = {f.stem: genome_length(f) for f in sorted(Path(a.gvmag_dir).glob("*.fna"))}
    if not lengths:
        sys.exit(f"No .fna files in {a.gvmag_dir}")

    completeness = {}
    if a.summary and Path(a.summary).exists():
        with open(a.summary) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                try:
                    completeness[row["gvmag"]] = float(row.get("checkv_completeness", ""))
                except (ValueError, TypeError):
                    pass

    # fastANI rows: query_path  ref_path  ANI  mapped_fragments  total_query_fragments
    # The aligned fraction is mapped/total, i.e. the fraction of the QUERY that aligned.
    # Running all-vs-all gives both directions, so for each pair take the fraction
    # belonging to the shorter genome, as aniclust's --min_tcov does.
    edges = {n: {} for n in lengths}
    with open(a.fastani) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 5:
                continue
            q, r = Path(f[0]).stem, Path(f[1]).stem
            if q == r or q not in lengths or r not in lengths:
                continue
            # Judge each pair on the row whose QUERY is the shorter genome, so `af`
            # is the aligned fraction of the shorter one. The reciprocal row is
            # skipped; all-vs-all guarantees the row we want exists.
            if lengths[q] > lengths[r]:
                continue
            ani, af = float(f[2]), 100.0 * float(f[3]) / float(f[4])
            if ani < a.min_ani or af < a.min_af:
                continue
            prev = edges[q].get(r)
            if prev is None or ani > prev[0]:
                edges[q][r] = edges[r][q] = (ani, af)

    def rank(n):
        return (-completeness.get(n, -1.0), -lengths[n], n)

    assigned, clusters = set(), []
    for rep in sorted(lengths, key=rank):
        if rep in assigned:
            continue
        assigned.add(rep)
        members = []
        for other in sorted(edges[rep], key=rank):
            if other not in assigned:
                assigned.add(other)
                members.append(other)
        clusters.append((rep, members))

    width = len(str(len(clusters)))
    with open(f"{a.out_prefix}_votus.tsv", "w") as tsv, \
         open(f"{a.out_prefix}_representatives.txt", "w") as reps:
        tsv.write("votu\tgvmag\tsample\tis_representative\trepresentative\t"
                  "length_bp\tcheckv_completeness\tani_to_rep\taf_to_rep\tn_members\n")
        for i, (rep, members) in enumerate(clusters, start=1):
            votu = f"vOTU{str(i).zfill(width)}"
            reps.write(f"{Path(a.gvmag_dir) / (rep + '.fna')}\n")
            for name in [rep] + members:
                ani, af = ("100.00", "100.00") if name == rep else \
                    tuple(f"{v:.2f}" for v in edges[rep][name])
                comp = completeness.get(name)
                tsv.write("\t".join([
                    votu, name, name.split("__")[0],
                    "yes" if name == rep else "no", rep,
                    str(lengths[name]), f"{comp:.2f}" if comp is not None else "NA",
                    ani, af, str(len(members) + 1)]) + "\n")

    multi = sum(1 for _, m in clusters if m)
    print(f"{len(lengths)} GVMAGs -> {len(clusters)} vOTUs "
          f"({multi} found in more than one sample/bin, "
          f"{len(clusters) - multi} singletons) "
          f"at >={a.min_ani}% ANI over >={a.min_af}% of the shorter genome")


if __name__ == "__main__":
    main()
