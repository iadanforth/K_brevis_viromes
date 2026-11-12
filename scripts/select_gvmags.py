#!/usr/bin/env python3
"""Pick the giant virus bins (GVMAGs) out of one sample's MetaBAT2 bins.

A bin is called a GVMAG when GVClass's majority-rule taxonomy places it in the
domain NCLDV (phylum Nucleocytoviricota), i.e. `taxonomy_majority` starts with
"d_NCLDV". Every GVMAG reported by Booker et al. 2026 (Supplemental Data 4) has
this classification. No completeness or contamination cut-off is applied here:
the paper kept GVMAGs from 66 % completeness up; quality columns are reported in
the summary table instead so they can be filtered later.

Selected bins are copied to OUT_DIR as <query>.fna and listed in
OUT_DIR/gvmag_list.tsv.
"""
import argparse
import csv
import shutil
import sys
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gvclass", required=True, help="GVClass gvclass_summary.tsv")
    p.add_argument("--bin-dir", required=True, help="directory of bin .fna files given to GVClass")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--sample", required=True)
    a = p.parse_args()

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    n_bins = n_gv = 0
    with open(a.gvclass) as fh, open(out / "gvmag_list.tsv", "w") as lst:
        lst.write("sample\tgvmag\ttaxonomy_majority\torder_completeness\tcellular_unique\n")
        for row in csv.DictReader(fh, delimiter="\t"):
            n_bins += 1
            tax = row.get("taxonomy_majority", "")
            if not tax.startswith("d_NCLDV"):
                continue
            query = row["query"]
            src = Path(a.bin_dir) / f"{query}.fna"
            if not src.exists():
                sys.exit(f"GVClass query {query} has no matching bin file {src}")
            shutil.copy(src, out / f"{query}.fna")
            lst.write(f"{a.sample}\t{query}\t{tax}\t{row.get('order_completeness', 'NA')}\t"
                      f"{row.get('cellular_unique', 'NA')}\n")
            n_gv += 1
    print(f"[{a.sample}] GVClass classified {n_bins} bins; {n_gv} are NCLDV (GVMAGs)")


if __name__ == "__main__":
    main()
