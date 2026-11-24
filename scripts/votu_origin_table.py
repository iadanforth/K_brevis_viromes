#!/usr/bin/env python3
"""One row per vOTU: which sample it was recovered from, and what kind of sample that was.

"Found in" here means ASSEMBLED FROM - the sample whose assembly produced the GVMAG that
was binned, quality-checked and clustered into this vOTU. That is how a vOTU enters the
study at all. (For the 7 environmental vOTUs it is also where read mapping detects them in
step 14; no vOTU is detected in a sample it was not assembled from.)

Six vOTUs were recovered from more than one sample, because the same virus assembled
independently in several incubations; those samples are listed together in one row, so the
table stays at one row per vOTU. All six are incubation-only, so the sample type is never
mixed within a row - the script checks this and fails loudly if that ever stops being true.

Input : output/votus/gvmags_votus.tsv (step 09)
Output: output/votu_origin.csv, and the same table printed as Markdown

Run:
    python3 -I scripts/votu_origin_table.py
"""
import argparse
import csv
from collections import OrderedDict

# The 9 runs from the bloom incubation experiments (CH1 May 2021, CH2 Aug 2021).
# Everything else in the study is an environmental sample.
INCUBATION = {"SRR28004137", "SRR28004138", "SRR28004139", "SRR28004140", "SRR28004141",
              "SRR28004142", "SRR28004143", "SRR28004144", "SRR28004154"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--votus", default="output/votus/gvmags_votus.tsv")
    ap.add_argument("--out", default="output/votu_origin.csv")
    a = ap.parse_args()

    samples = OrderedDict()
    with open(a.votus) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            samples.setdefault(row["votu"], []).append(row["sample"])

    rows = []
    for votu in sorted(samples):
        srr = sorted(set(samples[votu]))
        kinds = {"incubation" if s in INCUBATION else "environmental" for s in srr}
        if len(kinds) > 1:
            raise SystemExit(
                f"{votu} was assembled from both incubation and environmental samples "
                f"({', '.join(srr)}). This table assumes one kind per vOTU and would "
                f"misreport it; give it a 'mixed' value or split the row.")
        rows.append({"votu": votu, "sample": ", ".join(srr),
                     "sample_type": kinds.pop()})

    with open(a.out, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=["votu", "sample", "sample_type"])
        w.writeheader()
        w.writerows(rows)

    w1 = max(len(r["votu"]) for r in rows)
    w2 = max(len(r["sample"]) for r in rows)
    w3 = max(len(r["sample_type"]) for r in rows)
    print(f"| {'vOTU':<{w1}} | {'Sample':<{w2}} | {'Sample type':<{w3}} |")
    print(f"|{'-' * (w1 + 2)}|{'-' * (w2 + 2)}|{'-' * (w3 + 2)}|")
    for r in rows:
        print(f"| {r['votu']:<{w1}} | {r['sample']:<{w2}} | {r['sample_type']:<{w3}} |")

    n_inc = sum(1 for r in rows if r["sample_type"] == "incubation")
    print(f"\n{len(rows)} vOTUs: {n_inc} incubation, {len(rows) - n_inc} environmental")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
