#!/usr/bin/env python3
"""Summarise one sample's giant virus MAGs (GVMAGs) in one table.

One row per GVMAG, similar to Supplemental Data 4 of Booker et al. 2026:

  assembly statistics  (from the GVMAG .fna)    n_contigs, length_bp, gc_percent,
                                                n50_bp, longest_contig_bp
  CheckV 1.1.1         (quality_summary.tsv)    checkv_quality, checkv_completeness,
                                                checkv_completeness_method,
                                                checkv_contamination, checkv_viral_genes,
                                                checkv_host_genes, checkv_warnings
  geNomad annotate     (<gvmag>_genes.tsv)      genomad_genes, genomad_annotated
                                                (bitscore >= MIN_BITSCORE),
                                                genomad_virus_hallmarks
  GVClass              (gvclass_summary.tsv)    every GVClass column, prefixed
                                                "gvclass_" (copied as-is, so the table
                                                works whatever columns this GVClass
                                                version writes)

No GVMAG is filtered out here; filter on the quality columns afterwards.
Missing values are written as NA. Only the Python standard library is used.

Inputs are the files written by 08_gvmags.slurm:
  --gvmag-dir    OUT/gvmags/           <gvmag>.fna (one per GVMAG)
  --checkv       OUT/checkv/out/quality_summary.tsv  (contig_id = GVMAG name)
  --genomad-dir  OUT/genomad/          <gvmag>/<gvmag>_annotate/<gvmag>_genes.tsv
  --gvclass      OUT/gvclass/gvclass_summary.tsv     (query = GVMAG name)
"""
import argparse
import csv
import sys
from pathlib import Path

NA = "NA"


def read_fasta_lengths(path):
    """Return (list of contig lengths, total G+C count) for a FASTA file."""
    lengths, gc, cur = [], 0, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    lengths.append(cur)
                cur = 0
            else:
                s = line.strip().upper()
                cur = (cur or 0) + len(s)
                gc += s.count("G") + s.count("C")
    if cur is not None:
        lengths.append(cur)
    return lengths, gc


def n50(lengths):
    total, run = sum(lengths), 0
    for length in sorted(lengths, reverse=True):
        run += length
        if run * 2 >= total:
            return length
    return 0


def read_tsv_by_key(path, key):
    """Read a TSV into {row[key]: row}. Missing file -> ({}, [])."""
    if not path or not Path(path).is_file():
        return {}, []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = {row[key]: row for row in reader if row.get(key)}
        return rows, list(reader.fieldnames or [])


def genomad_counts(genes_tsv, min_bitscore):
    """Count genes, annotated genes (bitscore >= min) and virus hallmark genes."""
    if not genes_tsv.is_file():
        return NA, NA, NA
    n = n_ann = n_hall = 0
    with open(genes_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n += 1
            try:
                if float(row.get("bitscore", "NA")) >= min_bitscore:
                    n_ann += 1
            except ValueError:
                pass  # "NA" = no marker hit
            if row.get("virus_hallmark", "0") not in ("0", "", "NA"):
                n_hall += 1
    return n, n_ann, n_hall


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True)
    p.add_argument("--gvclass", required=True, help="GVClass gvclass_summary.tsv")
    p.add_argument("--gvmag-dir", required=True, help="folder of <gvmag>.fna files")
    p.add_argument("--checkv", required=True, help="CheckV quality_summary.tsv")
    p.add_argument("--genomad-dir", required=True, help="folder of geNomad annotate outputs")
    p.add_argument("--min-bitscore", type=float, default=50.0)
    p.add_argument("--out", required=True, help="output TSV")
    a = p.parse_args()

    fnas = sorted(Path(a.gvmag_dir).glob("*.fna"))
    if not fnas:
        sys.exit(f"[{a.sample}] no .fna files in {a.gvmag_dir}")

    checkv, _ = read_tsv_by_key(a.checkv, "contig_id")
    gvclass, gvclass_cols = read_tsv_by_key(a.gvclass, "query")
    gvclass_cols = [c for c in gvclass_cols if c != "query"]
    if not checkv:
        print(f"[{a.sample}] warning: no CheckV results read from {a.checkv}", file=sys.stderr)
    if not gvclass:
        print(f"[{a.sample}] warning: no GVClass results read from {a.gvclass}", file=sys.stderr)

    checkv_cols = [("checkv_quality", "checkv_quality"),
                   ("checkv_completeness", "completeness"),
                   ("checkv_completeness_method", "completeness_method"),
                   ("checkv_contamination", "contamination"),
                   ("checkv_viral_genes", "viral_genes"),
                   ("checkv_host_genes", "host_genes"),
                   ("checkv_warnings", "warnings")]
    bs = int(a.min_bitscore) if a.min_bitscore.is_integer() else a.min_bitscore
    header = (["sample", "gvmag", "n_contigs", "length_bp", "gc_percent", "n50_bp",
               "longest_contig_bp"]
              + [c for c, _ in checkv_cols]
              + ["genomad_genes", f"genomad_annotated_bitscore{bs}", "genomad_virus_hallmarks"]
              + [f"gvclass_{c}" for c in gvclass_cols])

    out = Path(a.out)
    tmp = out.with_name(out.name + ".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for fna in fnas:
            name = fna.stem
            lengths, gc = read_fasta_lengths(fna)
            total = sum(lengths)
            cv = checkv.get(name, {})
            gv = gvclass.get(name, {})
            genes_tsv = Path(a.genomad_dir) / name / f"{name}_annotate" / f"{name}_genes.tsv"
            row = [a.sample, name, len(lengths), total,
                   f"{100 * gc / total:.2f}" if total else NA,
                   n50(lengths), max(lengths) if lengths else 0]
            row += [cv.get(col) or NA for _, col in checkv_cols]
            row += list(genomad_counts(genes_tsv, a.min_bitscore))
            row += [gv.get(c) or NA for c in gvclass_cols]
            w.writerow(row)
    tmp.replace(out)
    print(f"[{a.sample}] wrote {len(fnas)} GVMAG rows to {out}")


if __name__ == "__main__":
    main()
