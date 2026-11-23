#!/usr/bin/env python3
"""Merge CoverM's per-sample tables into tidy abundance tables, for both filter settings.

Replaces merge_votu_abundance.py (2026-11-23, prompt 60 in PROMPT_LOG.md). The change:
step 14 now maps against each individual GVMAG rather than against one representative per
vOTU, following Booker et al. 2026's Methods ("...reads recruited to each GVMAG..."), and
it runs CoverM twice per sample - once at CoverM's defaults ("paper") and once with a 95%
read-identity filter ("strict"). This script reads both and writes them side by side.

WHAT COVERM'S OUTPUT LOOKS LIKE
-------------------------------
CoverM 0.7.0 writes ONE ROW PER GENOME, with the columns named after the SAMPLE, not the
genome: `Genome`, `<sample> Mean`, `<sample> Covered Bases`, `<sample> Length`,
`<sample> Read Count`, `<sample> RPKM`, `<sample> TPM`,
`<sample> Relative Abundance (%)`, `<sample> Variance`. (Verified on smoke-test job 73929 -
an earlier version of this script assumed the opposite layout and would have silently
matched nothing.) The `<sample>` prefix is taken from the R1 FASTQ filename, so it is not
relied on here; the sample name comes from the CoverM output FILENAME instead.

CoverM also emits an `unmapped` row, which is skipped because it is not a genome.

THE PRESENCE RULE
-----------------
Applied HERE rather than during mapping, so it is explicit and can be changed without
re-mapping ~100 GB of reads:

    a genome counts as PRESENT in a sample when >=75% of its length is covered

That is the MIUViG convention (Roux et al. 2019). BREADTH, not depth, is the gate: a
genuine detection covers most of the genome, whereas a few conserved genes attracting
stray reads from a relative covers only a little but can still show high average depth on
the part it does cover. Breadth saturates as 1 - e^-coverage, so >=75% breadth needs only
about 1.4x coverage - it is a weak demand on depth and a strong demand on evenness.

THE PER-vOTU ROLLUP, AND WHY IT IS NEEDED
-----------------------------------------
Mapping to 26 individual GVMAGs means the 6 vOTUs holding 2-3 near-identical (>=95% ANI)
members SPLIT their recruited reads between those members, so each member alone understates
the vOTU. The rollup therefore reports, per vOTU per sample:
  read_count_sum        reads recruited by the whole cluster (the split undone)
  best_covered_fraction breadth of the single best-covered member (breadths cannot be
                        summed - overlapping members would double-count)
  rep_*                 the representative's own values, for a like-for-like comparison
                        with the previous per-representative analysis

Outputs (all written to --out-dir):
  gvmag_abundance_long.tsv     one row per genome per sample per filter setting
  gvmag_abundance.csv          the same, trimmed to what the figures read
  votu_rollup.tsv              one row per vOTU per sample per filter
  gvmag_abundance_matrix.tsv   genome x sample matrix of relative abundance, per filter
"""
import argparse
import csv
import re
from collections import defaultdict
from datetime import date
from pathlib import Path

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

# The two filter settings run by step 14. Written to a column called `filter_set` rather
# than `filter`: in R, a column named `filter` shadows dplyr's verb inside tidy-eval and
# produces confusing failures.
FILTERS = ["paper", "strict"]


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


def read_coverm(path):
    """Return {genome: {metric: raw string}} for one CoverM table."""
    with open(path) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr)
        body = [r for r in rdr if r]
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
            f"{path.name}: could not find column(s) {missing} in CoverM output.\n"
            f"  header was: {header}\n"
            f"  This usually means the CoverM version changed its column names; "
            f"update METHODS at the top of this script.")
    out = {}
    for line in body:
        out[line[0].strip()] = {
            k: (line[i] if i < len(line) else "") for k, i in idx.items()}
    return out


def fnum(d, key, default=0.0):
    try:
        return float(d.get(key, ""))
    except (TypeError, ValueError):
        return default


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--coverm-dir", required=True,
                   help="abundance/per_gvmag, holding paper/ and strict/ subdirectories")
    p.add_argument("--manifest", required=True, help="gvmag_reference_manifest.tsv")
    p.add_argument("--sra-table", required=True, help="metadata/SraRunTable.csv")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--min-covered-fraction", type=float, default=0.75,
                   help="breadth required to call a genome present (default 0.75 = MIUViG)")
    a = p.parse_args()

    out_dir = Path(a.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    origin = {}
    with open(a.manifest) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            origin[row["genome"]] = row
    if not origin:
        raise SystemExit(f"{a.manifest} is empty")

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

    # Refuse to merge a half-finished run. The two filter arms are written by the same
    # array task but one after the other, so mid-run the "paper" directory is ahead of
    # "strict". Merging anyway produces a table that looks complete and is quietly wrong:
    # a statistic like "max breadth across samples" silently drops the samples that have
    # not landed, and a genome's own source sample going missing can turn a detection
    # into an absence. Observed for real on 2026-11-23 during job 74013.
    present = {}
    for filt in FILTERS:
        d = Path(a.coverm_dir) / filt
        present[filt] = {f.name.replace("_coverm.tsv", "") for f in d.glob("*_coverm.tsv")}
        if not present[filt]:
            raise SystemExit(f"No *_coverm.tsv in {d} - has step 14 finished?")
    if present[FILTERS[0]] != present[FILTERS[1]]:
        only = {f: sorted(present[f] - present[g])
                for f, g in (FILTERS, FILTERS[::-1])}
        raise SystemExit(
            "The filter arms cover different samples, so step 14 is still running or a "
            "task failed:\n"
            + "".join(f"  only in {f}: {', '.join(v) or '(none)'}\n"
                      for f, v in only.items())
            + "Wait for the array to finish, or check logs/14_gvmag_abundance_*.err.")
    n_samples_expected = len(present[FILTERS[0]])

    rows = []
    for filt in FILTERS:
        d = Path(a.coverm_dir) / filt
        files = sorted(d.glob("*_coverm.tsv"))
        for f in files:
            sample = f.name.replace("_coverm.tsv", "")
            table = read_coverm(f)
            m = meta.get(sample, {})
            for genome, vals in table.items():
                if genome not in origin:
                    continue      # skips CoverM's "unmapped" row and anything unexpected
                # CoverM 0.7.0 has no covered_fraction method, so breadth is computed here
                # from the two values it does report. Guard against a zero length, which
                # would otherwise divide by zero for a genome CoverM did not see.
                glen = fnum(vals, "genome_length")
                cbases = fnum(vals, "covered_bases")
                cf = (cbases / glen) if glen > 0 else 0.0
                o = origin[genome]
                rows.append({
                    "filter_set": filt,
                    "genome": genome,
                    "votu": o["votu"],
                    "votu_rank": o["votu_rank"],
                    "sample": sample,
                    "collection_date": m.get("collection_date", ""),
                    "date_precision": m.get("date_precision", "unknown"),
                    "isolation_source": m.get("isolation_source", "NA"),
                    "depth_m": m.get("depth", "NA"),
                    "genome_source_set": o["source_set"],
                    "genome_source_sample": o["source_sample"],
                    "genome_length": vals.get("genome_length", ""),
                    "covered_bases": vals.get("covered_bases", ""),
                    "mean_coverage": f"{fnum(vals, 'mean_coverage'):.6g}",
                    "covered_fraction": f"{cf:.6g}",
                    "rpkm": f"{fnum(vals, 'rpkm'):.6g}",
                    "tpm": f"{fnum(vals, 'tpm'):.6g}",
                    "relative_abundance_pct": f"{fnum(vals, 'relative_abundance'):.6g}",
                    "read_count": vals.get("read_count", ""),
                    "detected": "TRUE" if cf >= a.min_covered_fraction else "FALSE",
                })

    rows.sort(key=lambda r: (r["filter_set"], r["votu"], r["genome"], r["collection_date"]))

    def survey(r):
        # The nine Gulf samples are a monthly series; the two coastal Florida samples come
        # from a separate study and must not be joined to it.
        return ("Gulf of Mexico (monthly)" if "Gulf" in r["isolation_source"]
                else "coastal Florida (separate study)")

    long_cols = ["filter_set", "genome", "votu", "votu_rank", "sample", "collection_date",
                 "date_precision", "isolation_source", "depth_m", "genome_source_set",
                 "genome_source_sample", "genome_length", "covered_bases",
                 "mean_coverage", "covered_fraction", "rpkm", "tpm",
                 "relative_abundance_pct", "read_count", "detected"]
    with open(out_dir / "gvmag_abundance_long.tsv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=long_cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    # Self-contained CSV for the figures: everything they read, nothing else, so a figure
    # can be rebuilt from one small file with no other inputs.
    csv_cols = ["filter_set", "genome", "votu", "votu_rank", "sample", "collection_date",
                "survey", "genome_source_set", "genome_source_sample",
                "relative_abundance_pct", "rpkm", "mean_coverage", "covered_fraction",
                "read_count", "detected"]
    with open(out_dir / "gvmag_abundance.csv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=csv_cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: (survey(r) if k == "survey" else r[k]) for k in csv_cols})

    # Per-vOTU rollup - see the module docstring for why each statistic is what it is.
    grouped = defaultdict(list)
    for r in rows:
        grouped[(r["filter_set"], r["votu"], r["sample"])].append(r)
    roll = []
    for (filt, votu, sample), rs in grouped.items():
        rep = next((x for x in rs if x["votu_rank"] == "representative"), rs[0])
        best = max(rs, key=lambda x: float(x["covered_fraction"]))
        roll.append({
            "filter_set": filt,
            "votu": votu,
            "n_members": len(rs),
            "sample": sample,
            "collection_date": rep["collection_date"],
            "survey": survey(rep),
            "votu_source_set": rep["genome_source_set"],
            "read_count_sum": sum(int(float(x["read_count"] or 0)) for x in rs),
            "rpkm_sum": f"{sum(float(x['rpkm']) for x in rs):.6g}",
            "relative_abundance_pct_sum":
                f"{sum(float(x['relative_abundance_pct']) for x in rs):.6g}",
            "best_covered_fraction": best["covered_fraction"],
            "best_member": best["genome"],
            "rep_genome": rep["genome"],
            "rep_covered_fraction": rep["covered_fraction"],
            "rep_relative_abundance_pct": rep["relative_abundance_pct"],
            "rep_rpkm": rep["rpkm"],
            "detected": "TRUE" if float(best["covered_fraction"]) >= a.min_covered_fraction
                        else "FALSE",
        })
    roll.sort(key=lambda r: (r["filter_set"], r["votu"], r["collection_date"]))
    with open(out_dir / "votu_rollup.tsv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(roll[0].keys()),
                           delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(roll)

    # Quick-look matrix: genome x sample relative abundance, one block per filter.
    samples = sorted({(r["collection_date"], r["sample"]) for r in rows})
    genomes = sorted(origin)
    with open(out_dir / "gvmag_abundance_matrix.tsv", "w") as fo:
        for filt in FILTERS:
            lookup = {(r["genome"], r["sample"]): r["relative_abundance_pct"]
                      for r in rows if r["filter_set"] == filt}
            fo.write(f"# filter = {filt}; values are relative abundance (%)\n")
            fo.write("genome\tvotu\tsource_set\t" + "\t".join(s for _, s in samples) + "\n")
            fo.write("date\t\t\t" + "\t".join(d for d, _ in samples) + "\n")
            for g in genomes:
                o = origin[g]
                fo.write(f"{g}\t{o['votu']}\t{o['source_set']}\t"
                         + "\t".join(lookup.get((g, s), "0") for _, s in samples) + "\n")
            fo.write("\n")

    # ---- report ----
    n_s = len(samples)
    print(f"{len(rows)} rows = {len(genomes)} GVMAGs x {n_s} samples x {len(FILTERS)} filters")
    print(f"presence rule: covered_fraction >= {a.min_covered_fraction:.0%}\n")
    for filt in FILTERS:
        rs = [r for r in rows if r["filter_set"] == filt]
        det = [r for r in rs if r["detected"] == "TRUE"]
        inc = [r for r in det if r["genome_source_set"] == "incubation"]
        print(f"[{filt}] {len(det)} detections, {len(inc)} of them incubation GVMAGs")
    print(f"\nper-GVMAG summary (max over the {n_samples_expected} "
          f"environmental samples):")
    print(f"  {'genome':34s} {'vOTU':7s} {'set':6s} "
          f"{'paper reads':>12s} {'breadth':>8s} | {'strict reads':>12s} {'breadth':>8s}  det")
    for g in genomes:
        o = origin[g]
        cells = []
        for filt in FILTERS:
            rs = [r for r in rows if r["genome"] == g and r["filter_set"] == filt]
            cells.append((max((int(float(r["read_count"] or 0)) for r in rs), default=0),
                          max((float(r["covered_fraction"]) for r in rs), default=0.0),
                          sum(1 for r in rs if r["detected"] == "TRUE")))
        (pr, pb, _), (sr, sb, sd) = cells
        print(f"  {g:34s} {o['votu']:7s} {o['source_set'][:4]:6s} "
              f"{pr:12d} {pb:8.3f} | {sr:12d} {sb:8.3f}  {sd}")
    print(f"\nwrote {out_dir}/gvmag_abundance_long.tsv")
    print(f"wrote {out_dir}/gvmag_abundance.csv")
    print(f"wrote {out_dir}/votu_rollup.tsv")
    print(f"wrote {out_dir}/gvmag_abundance_matrix.tsv")


if __name__ == "__main__":
    main()
