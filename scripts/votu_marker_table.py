#!/usr/bin/env python3
"""Which of the seven NCLDV marker genes does each vOTU carry?

Step 10 already computes this - `ncldv_markersearch` writes an occupancy table as a
by-product - but it lands in scratch, which is purged after 90 days without access, and it
is keyed by the long alignment taxon name. This turns it into a committed, readable table.

WHAT A "MARKER" IS HERE, AND WHY IT HAS NO NAME
-----------------------------------------------
The phylogeny in step 11 is built from seven conserved genes that nearly all giant viruses
carry, so they can be aligned across the whole group. They are identified ONLY by their
GVOG accession (GVOGm0013, GVOGm0023, ...), because that is all the HMM models carry: a
check of `work/hmm/gvogs.common.hmm` found NAME lines holding the accession and no DESC,
ACC or any other free-text field, and the GVMAGs V2 proteomes have bare headers with no
functional annotation. Translating an accession to a familiar gene name (PolB, VLTF3, ...)
needs the Giant Virus Database's own GVOG annotation table, which this project does not
hold. **Do not assume the mapping** - README 10j once named the seven from general
knowledge and that list has never been verified against this run.

(The same directory holds a second, properly named HMM set, `NCLDV.hmm`: A32, D5, SFII,
mcp, mRNAc, PolB, RNAPL, RNAPS, RNR, VLTF3. The phylogeny did NOT use it - it used the
GVOG set - so those names do not map onto these columns.)

ONE GENOME PER vOTU
-------------------
Step 10 put one genome per vOTU into the alignment: the step 09 REPRESENTATIVE. So this
table describes the representative, not the union over the cluster's members. A
non-representative member could in principle carry a marker its representative lacks; for
the 6 multi-member vOTUs that is untested here.

Inputs
  output/tree/marker_occupancy_raw.tsv   ncldv_markersearch occupancy, copied out of
                                         scratch so the repo does not depend on it
  output/votus/gvmags_votus.tsv          vOTU membership, for the representative + sample
Output
  output/votus/votu_marker_genes.csv     one row per vOTU, one column per marker

Run:
    python3 -I scripts/votu_marker_table.py
"""
import argparse
import csv
import re

INCUBATION = {"SRR28004137", "SRR28004138", "SRR28004139", "SRR28004140", "SRR28004141",
              "SRR28004142", "SRR28004143", "SRR28004144", "SRR28004154"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--markers", default="output/tree/marker_occupancy_raw.tsv")
    ap.add_argument("--votus", default="output/votus/gvmags_votus.tsv")
    ap.add_argument("--out", default="output/votus/votu_marker_genes.csv")
    a = ap.parse_args()

    rep, origin = {}, {}
    for r in csv.DictReader(open(a.votus), delimiter="\t"):
        if r["is_representative"] == "yes":
            rep[r["votu"]] = r["gvmag"]
            origin[r["votu"]] = ("incubation" if r["sample"] in INCUBATION
                                 else "environmental")

    rows = list(csv.DictReader(open(a.markers), delimiter="\t"))
    cols = [c for c in rows[0] if c != "genome"]

    # The occupancy table's row name is the alignment taxon, "vOTU01__<sample>__<bin>_".
    # Pull the vOTU off the front rather than trying to rebuild the whole mangled name.
    found = {}
    for r in rows:
        m = re.match(r"(vOTU\d+)", r["genome"])
        if not m:
            continue                       # reference genomes, not ours
        found[m.group(1)] = {c: (float(r[c] or 0) > 0) for c in cols}

    missing = sorted(set(rep) - set(found))
    if missing:
        raise SystemExit(f"vOTUs absent from {a.markers}: {', '.join(missing)}. Every vOTU "
                         f"should have a row even when it carries zero markers; a missing "
                         f"row means the occupancy table does not match step 09's vOTUs.")

    out = []
    for votu in sorted(rep):
        hits = found[votu]
        n = sum(hits.values())
        row = {"votu": votu, "representative_gvmag": rep[votu],
               "sample_set": origin[votu]}
        row.update({c: ("yes" if hits[c] else "no") for c in cols})
        row["n_markers"] = n
        # ncldv_markersearch drops a taxon with no markers, so it never reaches the tree.
        row["in_tree"] = "yes" if n > 0 else "no"
        out.append(row)

    fields = ["votu", "representative_gvmag", "sample_set"] + cols + ["n_markers", "in_tree"]
    with open(a.out, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=fields)
        w.writeheader(); w.writerows(out)

    w1 = max(len(c) for c in cols)
    print(f"| {'vOTU':6s} | {'set':13s} | " + " | ".join(f"{c:{w1}s}" for c in cols)
          + " | n | tree |")
    print(f"|{'-'*8}|{'-'*15}|" + "|".join("-" * (w1 + 2) for c in cols) + "|---|------|")
    for r in out:
        print(f"| {r['votu']:6s} | {r['sample_set']:13s} | "
              + " | ".join(f"{('Y' if r[c]=='yes' else '.'): ^{w1}s}" for c in cols)
              + f" | {r['n_markers']} | {r['in_tree']:4s} |")

    full = sum(1 for r in out if r["n_markers"] == len(cols))
    none = [r["votu"] for r in out if r["n_markers"] == 0]
    print(f"\n{len(out)} vOTUs x {len(cols)} markers")
    print(f"  all {len(cols)} markers : {full}")
    print(f"  zero markers  : {len(none)}" + (f"  ({', '.join(none)}) - not in the tree"
                                              if none else ""))
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
