#!/usr/bin/env python3
"""For each vOTU representative, report its closest relative in GVMAGs V2.

GVMAGs V2 (Vasquez et al. 2025, doi:10.1101/2025.09.26.678796) is a collection of
18,725 giant virus genomes reduced to 8,508 clustered representatives. Comparing our
vOTU representatives to it answers "is this a virus someone has already sequenced,
and where was it found?".

Reading the output:
  * ANI >= 95% over most of the genome  -> the same viral species as a known genome
  * ANI 80-95%                          -> a relative of a known genome, new species
  * no row at all                       -> fastANI found nothing above ~80% ANI, so no
                                           close relative has been sequenced yet. This
                                           is NOT evidence the GVMAG is wrong; most
                                           environmental giant viruses have no close
                                           relative in any database.

Inputs:
  --fastani   fastANI output, query = vOTU representatives, ref = GVMAGs V2 genomes
  --votus     votu_cluster.py's <prefix>_votus.tsv
  --metadata  GVDB_stats.tsv from the GVMAGs V2 release
"""
import argparse
import csv
from pathlib import Path

# Columns copied from GVDB_stats.tsv, as (output name, column in the file).
META_COLS = [("ref_phylum", "Phylum"), ("ref_class", "Class"), ("ref_order", "Order"),
             ("ref_family", "gvclass_family"), ("ref_genus", "new_genus"),
             ("ref_length_bp", "LENbp"), ("ref_gc_percent", "GCperc"),
             ("ref_order_completeness", "order_completeness"),
             ("ref_ecosystem", "Ecosystem")]
NA = "NA"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fastani", required=True)
    p.add_argument("--votus", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    # Best hit per query: highest ANI, then largest aligned fraction.
    best = {}
    with open(a.fastani) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 5:
                continue
            q, r = Path(f[0]).stem, Path(f[1]).stem
            ani, af = float(f[2]), 100.0 * float(f[3]) / float(f[4])
            if q not in best or (ani, af) > (best[q][1], best[q][2]):
                best[q] = (r, ani, af)

    meta = {}
    with open(a.metadata) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            meta[row["genome"]] = row

    with open(a.votus) as fh:
        reps = [r for r in csv.DictReader(fh, delimiter="\t")
                if r["is_representative"] == "yes"]

    hits = 0
    with open(a.out, "w") as out:
        out.write("\t".join(["votu", "representative", "n_members", "length_bp",
                             "checkv_completeness", "nearest_reference",
                             "ani_percent", "aligned_fraction_percent"]
                            + [c for c, _ in META_COLS]) + "\n")
        for r in reps:
            name = r["gvmag"]
            row = [r["votu"], name, r["n_members"], r["length_bp"],
                   r["checkv_completeness"]]
            if name in best:
                ref, ani, af = best[name]
                hits += 1
                row += [ref, f"{ani:.2f}", f"{af:.2f}"]
                m = meta.get(ref, {})
                row += [(m.get(col) or NA) for _, col in META_COLS]
            else:
                row += ["none_above_fastani_threshold", NA, NA] + [NA] * len(META_COLS)
            out.write("\t".join(row) + "\n")

    print(f"{len(reps)} vOTU representatives: {hits} with a GVMAGs V2 relative, "
          f"{len(reps) - hits} with none (no close relative sequenced)")


if __name__ == "__main__":
    main()
