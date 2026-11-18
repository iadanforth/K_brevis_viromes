#!/usr/bin/env python3
"""Summarise ANI between our GVMAGs and the GVMAGs published by Booker et al. 2026.

Produces two CSVs:

  <out-summary>  ONE ROW PER OUR GVMAG - the answer to "does this match anything they
                 published?". Every one of our GVMAGs appears, including those with no
                 match, because "no match" is a result.
  <out-pairs>    every fastANI pair, for checking the summary or applying a different
                 threshold.

Reading ANI output correctly needs two things that are easy to get wrong:

1. fastANI reports the aligned fraction relative to the QUERY. A 200 kb genome fully
   contained in a 2 Mb genome aligns over 100% of the small one but 10% of the large one.
   Both runs are therefore read, and the fraction reported is the one belonging to the
   SHORTER genome, which is what the MIUViG definition uses.
2. fastANI reports nothing below roughly 80% ANI. A missing pair means "no relative within
   ~80% ANI", not "0% identity", and it is NOT evidence that the genomes are unrelated at
   a deeper taxonomic level - that is what the step 11 tree is for.

"same species" applies the MIUViG rule (>=95% ANI over >=85% of the shorter genome), the
same rule used to build our vOTUs in step 09, so the answer is consistent with the rest of
the project.
"""
import argparse
import csv
import glob
from pathlib import Path


def genome_length(path):
    n = 0
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                n += len(line.strip())
    return n


def read_fastani(path):
    """-> {(query_stem, ref_stem): (ani, aligned_fraction_of_query_percent)}"""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 5:
                continue
            q, r = Path(f[0]).stem, Path(f[1]).stem
            ani = float(f[2])
            af = 100.0 * float(f[3]) / float(f[4]) if float(f[4]) else 0.0
            prev = out.get((q, r))
            if prev is None or ani > prev[0]:
                out[(q, r)] = (ani, af)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--forward", required=True, help="ours as query")
    p.add_argument("--reverse", required=True, help="published as query")
    p.add_argument("--mine-dir", required=True)
    p.add_argument("--theirs-manifest", required=True)
    p.add_argument("--votus")
    p.add_argument("--summaries", nargs="*", default=[])
    p.add_argument("--out-summary", required=True)
    p.add_argument("--out-pairs", required=True)
    p.add_argument("--out-votu", help="optional vOTU-level correspondence CSV")
    p.add_argument("--min-ani", type=float, default=95.0)
    p.add_argument("--min-af", type=float, default=85.0)
    a = p.parse_args()

    mine = {f.stem: f for f in sorted(Path(a.mine_dir).glob("*.fna"))}
    my_len = {k: genome_length(v) for k, v in mine.items()}

    theirs, their_votu, their_len = {}, {}, {}
    with open(a.theirs_manifest) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            stem = Path(row["fasta"]).stem
            theirs[stem] = row["booker_gvmag"]
            their_votu[stem] = row["booker_votu"]
            their_len[stem] = int(row["total_bp"])

    # our own vOTU assignment and CheckV completeness, for context in the summary
    my_votu, my_sample = {}, {}
    if a.votus and Path(a.votus).exists():
        with open(a.votus) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                my_votu[row["gvmag"]] = row["votu"]
                my_sample[row["gvmag"]] = row["sample"]
    my_comp, my_tax = {}, {}
    for pat in a.summaries:
        for f in glob.glob(pat):
            sample = Path(f).parent.name
            with open(f) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    key = f"{sample}__{row.get('gvmag', '')}"
                    my_comp[key] = row.get("checkv_completeness", "")
                    my_tax[key] = row.get("gvclass_taxonomy_majority", "")

    fwd = read_fastani(a.forward)     # (ours, theirs) -> af of OURS
    rev = read_fastani(a.reverse)     # (theirs, ours) -> af of THEIRS

    pairs = []
    for (q, r), (ani, af_q) in fwd.items():
        af_r = rev.get((r, q), (None, None))[1]
        ani_r = rev.get((r, q), (None, None))[0]
        # the aligned fraction belonging to the shorter genome
        if my_len.get(q, 0) <= their_len.get(r, 0):
            af_short, which = af_q, "ours"
        else:
            af_short, which = (af_r if af_r is not None else af_q), ("published" if af_r is not None else "ours")
        pairs.append({
            "our_gvmag": q, "our_votu": my_votu.get(q, "NA"),
            "our_length_bp": my_len.get(q, 0),
            "booker_gvmag": theirs.get(r, r), "booker_votu": their_votu.get(r, "NA"),
            "booker_length_bp": their_len.get(r, 0),
            "ani_percent": f"{ani:.3f}",
            "aligned_fraction_of_ours_percent": f"{af_q:.3f}",
            "aligned_fraction_of_published_percent": "NA" if af_r is None else f"{af_r:.3f}",
            "aligned_fraction_of_shorter_percent": f"{af_short:.3f}",
            "shorter_genome_is": which,
            "reciprocal_ani_percent": "NA" if ani_r is None else f"{ani_r:.3f}",
            "same_species_miuvig": "yes" if (ani >= a.min_ani and af_short >= a.min_af) else "no",
        })
    # also keep reverse-only pairs (fastANI is not perfectly symmetric)
    for (r, q), (ani, af_r) in rev.items():
        if (q, r) in fwd:
            continue
        af_short = af_r if their_len.get(r, 0) <= my_len.get(q, 0) else af_r
        pairs.append({
            "our_gvmag": q, "our_votu": my_votu.get(q, "NA"),
            "our_length_bp": my_len.get(q, 0),
            "booker_gvmag": theirs.get(r, r), "booker_votu": their_votu.get(r, "NA"),
            "booker_length_bp": their_len.get(r, 0),
            "ani_percent": f"{ani:.3f}",
            "aligned_fraction_of_ours_percent": "NA",
            "aligned_fraction_of_published_percent": f"{af_r:.3f}",
            "aligned_fraction_of_shorter_percent": f"{af_short:.3f}",
            "shorter_genome_is": "published",
            "reciprocal_ani_percent": "NA",
            "same_species_miuvig": "yes" if (ani >= a.min_ani and af_short >= a.min_af) else "no",
            "note": "detected only with the published genome as query",
        })

    pair_cols = ["our_gvmag", "our_votu", "our_length_bp", "booker_gvmag", "booker_votu",
                 "booker_length_bp", "ani_percent", "aligned_fraction_of_ours_percent",
                 "aligned_fraction_of_published_percent",
                 "aligned_fraction_of_shorter_percent", "shorter_genome_is",
                 "reciprocal_ani_percent", "same_species_miuvig", "note"]
    pairs.sort(key=lambda d: (d["our_gvmag"], -float(d["ani_percent"])))
    with open(a.out_pairs, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=pair_cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(pairs)

    # one row per our GVMAG
    rows = []
    for g in sorted(mine):
        hits = [d for d in pairs if d["our_gvmag"] == g]
        hits.sort(key=lambda d: (-float(d["ani_percent"]),
                                 -float(d["aligned_fraction_of_shorter_percent"])))
        best = hits[0] if hits else None
        rows.append({
            "our_gvmag": g,
            "our_votu": my_votu.get(g, "NA"),
            "our_sample": my_sample.get(g, g.split("__")[0]),
            "our_length_bp": my_len.get(g, 0),
            "our_checkv_completeness": my_comp.get(g, "NA") or "NA",
            "our_gvclass_taxonomy": my_tax.get(g, "NA") or "NA",
            "n_published_matches": len(hits),
            "best_match_booker_gvmag": best["booker_gvmag"] if best else "none_above_fastani_threshold",
            "best_match_booker_votu": best["booker_votu"] if best else "NA",
            "best_match_length_bp": best["booker_length_bp"] if best else "NA",
            "ani_percent": best["ani_percent"] if best else "NA",
            "aligned_fraction_of_shorter_percent": best["aligned_fraction_of_shorter_percent"] if best else "NA",
            "same_species_miuvig": best["same_species_miuvig"] if best else "no_match",
        })
    sum_cols = list(rows[0].keys())
    with open(a.out_summary, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=sum_cols)
        w.writeheader()
        w.writerows(rows)

    # vOTU-level rollup: the most interpretable form of the answer, because a vOTU is
    # the species-level unit. Several of our GVMAGs can belong to one vOTU and match the
    # same published genome, which at GVMAG level looks like several results and at vOTU
    # level is one.
    if a.out_votu:
        import collections
        by_votu = collections.defaultdict(lambda: {"gvmags": [], "pub": set(),
                                                   "pubvotu": set(), "ani": []})
        for r in rows:
            v = by_votu[r["our_votu"]]
            v["gvmags"].append(r["our_gvmag"])
            if r["same_species_miuvig"] == "yes":
                v["pub"].add(r["best_match_booker_gvmag"])
                v["pubvotu"].add(r["best_match_booker_votu"])
                v["ani"].append(float(r["ani_percent"]))
        vcols = ["our_votu", "n_our_gvmags", "matches_published", "published_gvmags",
                 "published_votu", "min_ani_percent", "max_ani_percent"]
        with open(a.out_votu, "w", newline="") as fo:
            w = csv.DictWriter(fo, fieldnames=vcols)
            w.writeheader()
            for v in sorted(by_votu):
                d = by_votu[v]
                w.writerow({
                    "our_votu": v,
                    "n_our_gvmags": len(d["gvmags"]),
                    "matches_published": "yes" if d["pub"] else "no",
                    "published_gvmags": "; ".join(sorted(d["pub"])) or "none",
                    "published_votu": "; ".join(sorted(x for x in d["pubvotu"] if x != "NA")) or "NA",
                    "min_ani_percent": f"{min(d['ani']):.3f}" if d["ani"] else "NA",
                    "max_ani_percent": f"{max(d['ani']):.3f}" if d["ani"] else "NA",
                })
        matched = sum(1 for v in by_votu.values() if v["pub"])
        print(f"vOTU-level: {matched} of {len(by_votu)} of our vOTUs match a published GVMAG")
        print(f"wrote {a.out_votu}")

    same = [r for r in rows if r["same_species_miuvig"] == "yes"]
    any_hit = [r for r in rows if r["n_published_matches"] > 0]
    print(f"{len(mine)} of our GVMAGs vs {len(theirs)} published GVMAGs")
    print(f"  with any fastANI hit (>=~80% ANI): {len(any_hit)}")
    print(f"  same species by MIUViG (>={a.min_ani}% ANI over >={a.min_af}% of the shorter): {len(same)}")
    if same:
        print("\n  matches:")
        for r in same:
            print(f"    {r['our_gvmag']:42s} ({r['our_votu']})  ==  {r['best_match_booker_gvmag']:26s} "
                  f"({r['best_match_booker_votu']})  ANI {r['ani_percent']}%  AF {r['aligned_fraction_of_shorter_percent']}%")
    sub = [r for r in any_hit if r["same_species_miuvig"] != "yes"]
    if sub:
        print("\n  related but below the species threshold:")
        for r in sub:
            print(f"    {r['our_gvmag']:42s} -> {r['best_match_booker_gvmag']:26s} "
                  f"ANI {r['ani_percent']}%  AF {r['aligned_fraction_of_shorter_percent']}%")
    print(f"\nwrote {a.out_summary}\nwrote {a.out_pairs}")


if __name__ == "__main__":
    main()
