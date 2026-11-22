#!/usr/bin/env python3
"""Merge CoverM's per-sample tables into one tidy abundance-through-time table.

CoverM 0.7.0 writes ONE ROW PER GENOME, with the columns named after the SAMPLE, not the
genome: `Genome`, `<sample> Mean`, `<sample> Covered Bases`, `<sample> Length`,
`<sample> Read Count`, `<sample> RPKM`, `<sample> TPM`,
`<sample> Relative Abundance (%)`, `<sample> Variance`. (Verified on smoke-test job 73929 -
an earlier version of this script assumed the opposite layout, one row with
`<genome> <method>` columns, and would have silently matched nothing.) The sample prefix is
taken from the R1 FASTQ filename, so it is not relied on here; the sample name comes from
the CoverM output filename instead.

CoverM also emits an `unmapped` row, which is skipped because it is not a vOTU.

This produces one long/tidy table - one row per vOTU per sample - with the collection date
joined on, so R can plot abundance against time directly.

The presence rule is applied HERE rather than during mapping, so it is explicit and can be
changed without re-mapping ~200 GB of reads:

    a vOTU counts as PRESENT in a sample when >=75% of its genome is covered

That is the MIUViG convention (Roux et al. 2019): a genuine detection covers most of the
genome, whereas a few conserved genes attracting stray reads covers only a little. Depth
alone cannot distinguish the two, which is why breadth is the gate and coverage/RPKM are
only interpreted once it passes. `rpkm_detected` is RPKM with below-threshold values set to
0, which is the column to plot; `rpkm` keeps the raw value so a different threshold can be
explored.

Outputs:
  votu_abundance_long.tsv   tidy: one row per vOTU per sample, with dates and flags
  votu_abundance_matrix.tsv vOTU x sample matrix of rpkm_detected, for a quick look
"""
import argparse
import csv
import re
from datetime import date
from pathlib import Path

# CoverM method name -> our column name. CoverM's exact header spelling varies by version,
# so matching is done case-insensitively on a substring.
# CoverM's exact column suffixes (0.7.0) -> our column names. Matched case-insensitively
# on the END of the header, because the header is "<sample> <Method>" and the sample name
# can contain anything. "covered bases" is listed before "bases"-like suffixes and
# "read count" before "count" so the more specific name wins.
METHODS = [
    ("covered bases", "covered_bases"),
    ("relative abundance (%)", "relative_abundance"),
    ("relative abundance", "relative_abundance"),
    ("read count", "read_count"),
    ("trimmed mean", "trimmed_mean"),
    ("mean", "mean_coverage"),
    ("rpkm", "rpkm"),
    ("tpm", "tpm"),
    ("length", "genome_length"),
    ("variance", "variance"),
]


def parse_date(raw):
    """Collection_Date is 'YYYY-MM' for the Gulf series and 'YYYY-MM-DD' for the two
    coastal Florida samples. A month is represented by its first day so both sit on one
    axis; `date_precision` records which it was, so the figure can say so."""
    raw = (raw or "").strip()
    m = re.match(r"^(\d{4})-(\d{2})-(\d{2})$", raw)
    if m:
        y, mo, d = (int(x) for x in m.groups())
        return date(y, mo, d).isoformat(), "day"
    m = re.match(r"^(\d{4})-(\d{2})$", raw)
    if m:
        y, mo = (int(x) for x in m.groups())
        return date(y, mo, 1).isoformat(), "month"
    return "", "unknown"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--coverm-dir", required=True)
    p.add_argument("--manifest", required=True, help="votu_reference_manifest.tsv")
    p.add_argument("--sra-table", required=True, help="metadata/SraRunTable.csv")
    p.add_argument("--out-long", required=True)
    p.add_argument("--out-matrix", required=True)
    p.add_argument("--out-csv", help="tidy CSV holding exactly what the figure needs")
    p.add_argument("--only-source-set", choices=["environmental", "incubation"],
                   help="restrict --out-csv to vOTUs discovered in this sample set. "
                        "'environmental' gives the 7 vOTUs that were assembled from an "
                        "environmental sample, which are also exactly the 7 ever detected "
                        "in one.")
    p.add_argument("--min-covered-fraction", type=float, default=0.75,
                   help="breadth required to call a vOTU present (default 0.75 = MIUViG)")
    a = p.parse_args()

    origin = {}
    with open(a.manifest) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            origin[row["votu"]] = row

    meta = {}
    with open(a.sra_table) as fh:
        for row in csv.DictReader(fh):
            iso, prec = parse_date(row.get("Collection_Date"))
            meta[row["Run"]] = {
                "collection_date": iso,
                "date_precision": prec,
                "isolation_source": (row.get("isolation_source") or "").strip() or "NA",
                "depth": (row.get("Depth") or "").strip() or "NA",
            }

    rows = []
    files = sorted(Path(a.coverm_dir).glob("*_coverm.tsv"))
    if not files:
        raise SystemExit(f"No *_coverm.tsv in {a.coverm_dir}")

    for f in files:
        sample = f.name.replace("_coverm.tsv", "")
        with open(f) as fh:
            rdr = csv.reader(fh, delimiter="\t")
            header = next(rdr)
            body = [r for r in rdr if r]
        # Column 0 is "Genome"; every other column is "<sample> <Method>", so the method
        # is identified by the header's suffix and applies to whichever genome's row we
        # are on.
        idx = {}
        for i, h in enumerate(header[1:], start=1):
            h_l = h.lower().strip()
            for needle, name in METHODS:
                if h_l.endswith(needle):
                    idx[name] = i
                    break
        missing = [n for n in ("covered_bases", "genome_length", "rpkm") if n not in idx]
        if missing:
            raise SystemExit(
                f"{f.name}: could not find column(s) {missing} in CoverM output.\n"
                f"  header was: {header}\n"
                f"  This usually means the CoverM version changed its column names; "
                f"update METHODS at the top of this script.")

        for line in body:
            genome = line[0].strip()
            if genome not in origin:
                continue          # skips CoverM's "unmapped" row and anything unexpected

            def val(key):
                i = idx.get(key)
                if i is None or i >= len(line):
                    return ""
                return line[i]

            def f_(key, default=0.0):
                try:
                    return float(val(key))
                except (TypeError, ValueError):
                    return default

            # CoverM 0.7.0 has no covered_fraction method, so breadth is computed
            # here from the two values it does report. Guard against a zero length,
            # which would otherwise divide by zero for a genome CoverM did not see.
            glen = f_("genome_length")
            cbases = f_("covered_bases")
            cf = (cbases / glen) if glen > 0 else 0.0
            rpkm = f_("rpkm")
            present = cf >= a.min_covered_fraction
            m = meta.get(sample, {})
            o = origin[genome]
            rows.append({
                "votu": genome,
                "sample": sample,
                "collection_date": m.get("collection_date", ""),
                "date_precision": m.get("date_precision", "unknown"),
                "isolation_source": m.get("isolation_source", "NA"),
                "depth_m": m.get("depth", "NA"),
                "votu_source_set": o["source_set"],
                "votu_source_sample": o["source_sample"],
                "genome_length": val("genome_length"),
                "covered_bases": val("covered_bases"),
                "mean_coverage": f"{f_('mean_coverage'):.6g}",
                "covered_fraction": f"{cf:.6g}",
                "rpkm": f"{rpkm:.6g}",
                "rpkm_detected": f"{(rpkm if present else 0.0):.6g}",
                "tpm": f"{f_('tpm'):.6g}",
                "relative_abundance_pct": f"{f_('relative_abundance'):.6g}",
                "read_count": val("read_count"),
                "present": "TRUE" if present else "FALSE",
            })

    rows.sort(key=lambda r: (r["votu"], r["collection_date"]))
    cols_out = ["votu", "sample", "collection_date", "date_precision", "isolation_source",
                "depth_m", "votu_source_set", "votu_source_sample", "genome_length",
                "covered_bases", "mean_coverage", "covered_fraction", "rpkm",
                "rpkm_detected", "tpm", "relative_abundance_pct", "read_count", "present"]
    with open(a.out_long, "w") as out:
        w = csv.DictWriter(out, fieldnames=cols_out, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    votus = sorted({r["votu"] for r in rows})
    samples = sorted({(r["collection_date"], r["sample"]) for r in rows})
    with open(a.out_matrix, "w") as out:
        out.write("votu\tsource_set\t" + "\t".join(s for _, s in samples) + "\n")
        out.write("date\t\t" + "\t".join(d for d, _ in samples) + "\n")
        lookup = {(r["votu"], r["sample"]): r["rpkm_detected"] for r in rows}
        for v in votus:
            src = origin[v]["source_set"]
            out.write(f"{v}\t{src}\t" + "\t".join(lookup.get((v, s), "0") for _, s in samples) + "\n")

    # A self-contained CSV for the figure: one row per vOTU per environmental sample,
    # carrying only the columns plot_votu_abundance.R reads. Kept separate from the long
    # TSV so the figure can be regenerated from one small file with no other inputs.
    if a.out_csv:
        csv_cols = ["votu", "sample", "collection_date", "survey", "votu_source_set",
                    "relative_abundance_pct", "rpkm", "mean_coverage",
                    "covered_fraction", "detected"]
        keep = [r for r in rows
                if not a.only_source_set or r["votu_source_set"] == a.only_source_set]
        if a.only_source_set:
            n_v = len({r["votu"] for r in keep})
            print(f"--only-source-set {a.only_source_set}: {n_v} vOTUs, {len(keep)} rows")
        with open(a.out_csv, "w", newline="") as fo:
            w = csv.DictWriter(fo, fieldnames=csv_cols)
            w.writeheader()
            for r in keep:
                w.writerow({
                    "votu": r["votu"],
                    "sample": r["sample"],
                    "collection_date": r["collection_date"],
                    # The nine Gulf samples are a monthly series; the two coastal Florida
                    # samples come from a separate study and must not be joined to it.
                    "survey": ("Gulf of Mexico (monthly)" if "Gulf" in r["isolation_source"]
                               else "coastal Florida (separate study)"),
                    "votu_source_set": r["votu_source_set"],
                    "relative_abundance_pct": r["relative_abundance_pct"],
                    "rpkm": r["rpkm"],
                    "mean_coverage": r["mean_coverage"],
                    "covered_fraction": r["covered_fraction"],
                    "detected": "TRUE" if r["present"] == "TRUE" else "FALSE",
                })
        print(f"wrote {a.out_csv}")

    n_present = sum(1 for r in rows if r["present"] == "TRUE")
    print(f"{len(rows)} rows ({len(votus)} vOTUs x {len(samples)} samples)")
    print(f"detections at >={a.min_covered_fraction:.0%} breadth: {n_present}")
    print("\nper-vOTU detections:")
    for v in votus:
        rs = [r for r in rows if r["votu"] == v]
        det = [r for r in rs if r["present"] == "TRUE"]
        src = origin[v]["source_set"][:4]
        best = max((float(r["rpkm"]) for r in rs), default=0.0)
        print(f"  {v:8s} {src:5s} detected in {len(det):2d}/{len(rs):2d} samples, max RPKM {best:9.3f}")
    print(f"\nwrote {a.out_long}\nwrote {a.out_matrix}")


if __name__ == "__main__":
    main()
