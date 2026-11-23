#!/usr/bin/env python3
"""Which vOTUs appeared in which samples?

"Appeared" means two different things here, and conflating them would misrepresent the
data, so both are reported and kept distinct:

  ASSEMBLED  the vOTU was recovered from that sample - a GVMAG was binned out of that
             sample's assembly and clustered into this vOTU. This is available for all
             20 samples.
  DETECTED   that sample's reads were mapped back to the vOTU representative and covered
             >=75% of the genome (MIUViG). This is available for the 11 ENVIRONMENTAL
             samples only, because they are the only ones step 14 mapped.
  mapped     reads were mapped and the vOTU was not detected: a true absence.
  -          not assessed - the 9 incubation samples were never read-mapped, so a blank
             here means "unknown", NOT "absent". That asymmetry is the main thing this
             table has to keep honest.

Outputs a long table (one row per vOTU per sample) and a wide matrix.
"""
import argparse
import csv
import glob
import os
from pathlib import Path

STATUS = {"assembled+detected": "ASSEMBLED+DETECTED", "assembled": "ASSEMBLED",
          "detected": "DETECTED", "mapped": "mapped, not detected", "none": "not assessed"}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--votus", required=True, help="gvmags_votus.tsv from step 09")
    p.add_argument("--coverm-dir", required=True, help="abundance/per_sample")
    p.add_argument("--sra-table", required=True)
    p.add_argument("--incubation", required=True,
                   help="comma-separated incubation sample accessions")
    p.add_argument("--min-breadth", type=float, default=0.75)
    p.add_argument("--out-long", required=True)
    p.add_argument("--out-matrix", required=True)
    a = p.parse_args()

    inc = {s.strip() for s in a.incubation.split(",") if s.strip()}

    # every sample in the study, with its collection date
    samples, date_of = [], {}
    with open(a.sra_table) as fh:
        for row in csv.DictReader(fh):
            samples.append(row["Run"])
            date_of[row["Run"]] = (row.get("Collection_Date") or "").strip()
    samples = sorted(set(samples), key=lambda s: (s not in inc, date_of.get(s, ""), s))

    # which sample each vOTU was assembled from (a vOTU can come from several)
    assembled = {}
    votus = set()
    with open(a.votus) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            votus.add(row["votu"])
            assembled.setdefault(row["votu"], set()).add(row["sample"])

    # read-mapping detections, from the CoverM tables
    detected, mapped_samples = {}, set()
    for f in sorted(glob.glob(os.path.join(a.coverm_dir, "*_coverm.tsv"))):
        samp = os.path.basename(f).replace("_coverm.tsv", "")
        mapped_samples.add(samp)
        rows = list(csv.reader(open(f), delimiter="\t"))
        hdr, body = rows[0], rows[1:]
        idx = {}
        for i, h in enumerate(hdr[1:], 1):
            hl = h.lower().strip()
            for needle, key in (("covered bases", "cb"), ("length", "ln")):
                if hl.endswith(needle):
                    idx[key] = i
                    break
        for r in body:
            g = r[0].strip()
            if not g.startswith("vOTU"):
                continue
            try:
                cb, ln = float(r[idx["cb"]]), float(r[idx["ln"]])
            except Exception:
                continue
            if ln and cb / ln >= a.min_breadth:
                detected.setdefault(g, set()).add(samp)

    votus = sorted(votus)
    long_rows = []
    for v in votus:
        for s in samples:
            asm = s in assembled.get(v, set())
            det = s in detected.get(v, set())
            if asm and det:
                st = "assembled+detected"
            elif asm:
                st = "assembled"
            elif det:
                st = "detected"
            elif s in mapped_samples:
                st = "mapped"
            else:
                st = "none"
            long_rows.append({
                "votu": v, "sample": s, "collection_date": date_of.get(s, ""),
                "sample_set": "incubation" if s in inc else "environmental",
                "read_mapped": "yes" if s in mapped_samples else "no",
                "status": st, "status_label": STATUS[st],
            })

    with open(a.out_long, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(long_rows[0].keys()))
        w.writeheader(); w.writerows(long_rows)

    SHORT = {"assembled+detected": "A+D", "assembled": "A", "detected": "D",
             "mapped": ".", "none": ""}
    look = {(r["votu"], r["sample"]): SHORT[r["status"]] for r in long_rows}
    with open(a.out_matrix, "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["votu"] + samples)
        w.writerow(["sample_set"] + ["incubation" if s in inc else "environmental" for s in samples])
        w.writerow(["collection_date"] + [date_of.get(s, "") for s in samples])
        for v in votus:
            w.writerow([v] + [look[(v, s)] for s in samples])

    n_asm = sum(1 for r in long_rows if r["status"].startswith("assembled"))
    n_det = sum(1 for r in long_rows if "detected" in r["status"])
    print(f"{len(votus)} vOTUs x {len(samples)} samples")
    print(f"  assembled from a sample : {n_asm} cells")
    print(f"  detected by mapping     : {n_det} cells")
    print(f"  read-mapped samples     : {len(mapped_samples)} of {len(samples)}"
          f"  (the other {len(samples)-len(mapped_samples)} are 'not assessed', not 'absent')")
    print(f"\nwrote {a.out_long}\nwrote {a.out_matrix}")


if __name__ == "__main__":
    main()
