#!/usr/bin/env python3
"""Choose giant virus reference genomes from GVMAGs V2 for the phylogenomic tree.

GVMAGs V2 (Vasquez et al. 2025, doi:10.1101/2025.09.26.678796) is the largest curated
giant virus genome collection: 18,725 genomes reduced to 8,508 ANI-cluster
representatives, whose protein FASTAs ship in GVMAGS_V2_database.gz. Putting all 8,508
into one concatenated alignment is possible but slow, and far beyond what a Bayesian
phylogenetic program can comfortably handle, so this script picks a defensible subset.

Two nested tiers are available, both read from the release's GVDB_stats.tsv:

  --tier mcl    the coarser MCL clustering (column `mcl_rep`): 2,606 representatives.
                Closest in scale to the published tree, which had 2,037 genomes after
                filtering. This is the default.
  --tier skani  every ANI-cluster representative (column `skani_rep`): 8,507. Finest
                resolution, slowest.

On top of the tier you can restrict by viral order, which is usually what you want: a
tree of Imitervirales genomes is far more informative about Imitervirales vOTUs than one
diluted with distant orders. `--orders auto` reads the orders of your own GVMAGs from the
step 08 GVClass summaries and keeps references from those orders, plus (unless
--no-outgroup) a small, even sample of other orders to root the tree.

--max-refs caps the number of in-order references, sampled evenly across families so the
subset stays diverse rather than dominated by one clade.

Output: the chosen references' .faa files copied into --out-dir with sanitised names
(no dots, because ncldv_markersearch derives the taxon name from the file name), and a
manifest TSV recording what was chosen and why.
"""
import argparse
import csv
import random
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path


def sanitise(name):
    """ncldv_markersearch uses the .faa basename as the taxon name and does regex
    substitutions on '.faa', so dots and whitespace in names must go."""
    return re.sub(r"[^A-Za-z0-9_-]", "_", name)


def order_of(row):
    """The viral order of one step 08 GVMAG summary row, as a bare name ("Imitervirales").

    Read from `gvclass_taxonomy_majority`, which is the only clean source. GVClass's
    `gvclass_order` column is a VOTE STRING, not a label - a real row reads
    "NCLDV__Imitervirales:10(100.00%)", and one with a contested order reads
    "NCLDV__Imitervirales:20(55.56%),BAC__UM-FILTER-46-8:8(22.22%),...". Using that column
    directly (as this script did until 2026-11-16) put those whole strings into the wanted
    set, where they matched none of GVDB_stats.tsv's `Order` values, so NO in-order
    reference was ever selected and the tree was silently built from outgroups alone.

    `gvclass_taxonomy_majority` is the resolved lineage,
    "d_NCLDV;p_Nucleocytoviricota;c_Megaviricetes;o_Imitervirales;f_IM_01;g_;s_singleton",
    so the order is the o_ field. An unresolved order is an empty "o_", which returns None
    and correctly contributes nothing.
    """
    tax = (row.get("gvclass_taxonomy_majority")
           or row.get("taxonomy_majority") or "").strip()
    for field in tax.split(";"):
        field = field.strip()
        if field.startswith("o_"):
            name = field[2:].strip()
            return name or None
    # Fallback: pull the order out of the vote string if the lineage column is absent.
    vote = (row.get("gvclass_order") or "").strip()
    m = re.search(r"NCLDV__([A-Za-z_]+?):", vote)
    return m.group(1) if m else None


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--metadata", required=True, help="GVDB_stats.tsv")
    p.add_argument("--faa-dir", required=True, help="directory of reference .faa files")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--tier", choices=["mcl", "skani"], default="mcl")
    p.add_argument("--orders", default="auto",
                   help="'auto' (from --gvclass-summaries), 'all', or a comma-separated list")
    p.add_argument("--gvclass-summaries", nargs="*", default=[],
                   help="step 08 *_gvmag_summary.tsv files, for --orders auto")
    p.add_argument("--max-refs", type=int, default=0, help="0 = no cap")
    p.add_argument("--outgroup-per-order", type=int, default=3,
                   help="references kept from each OTHER order, to root the tree (0 = none)")
    p.add_argument("--seed", type=int, default=1)
    a = p.parse_args()

    rows = list(csv.DictReader(open(a.metadata), delimiter="\t"))
    by_genome = {r["genome"]: r for r in rows}
    col = "mcl_rep" if a.tier == "mcl" else "skani_rep"
    reps = sorted({r[col] for r in rows if r[col] not in ("NA", "")})
    print(f"[refs] tier '{a.tier}' ({col}): {len(reps)} representatives in the metadata")

    # ---- which orders do we care about? ----
    if a.orders == "all":
        wanted = None
    elif a.orders == "auto":
        wanted = set()
        for f in a.gvclass_summaries:
            if not Path(f).exists():
                continue
            with open(f) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    o = order_of(row)
                    if o:
                        wanted.add(o)
        if not wanted:
            sys.exit("--orders auto found no orders in the GVClass summaries; pass --orders "
                     "explicitly (e.g. --orders Imitervirales) or 'all'")
        print(f"[refs] orders taken from our own GVMAGs: {', '.join(sorted(wanted))}")
    else:
        wanted = {o.strip() for o in a.orders.split(",") if o.strip()}
        print(f"[refs] orders requested: {', '.join(sorted(wanted))}")

    # ---- split representatives into in-order and other-order ----
    in_order, others_by_order = [], defaultdict(list)
    for g in reps:
        r = by_genome.get(g)
        order = (r or {}).get("Order", "NA")
        if wanted is None or order in wanted:
            in_order.append(g)
        else:
            others_by_order[order].append(g)
    print(f"[refs] {len(in_order)} representatives in those orders, "
          f"{sum(len(v) for v in others_by_order.values())} in others")
    if wanted is not None and not in_order:
        sys.exit("No reference genome matched the requested orders "
                 f"({', '.join(sorted(wanted))}). Those names must match GVDB_stats.tsv's "
                 "`Order` column exactly (e.g. 'Imitervirales', no 'o_' prefix). Refusing "
                 "to continue: the alignment would contain only outgroups, which looks "
                 "like a valid tree but is not the one you asked for.")

    # ---- cap, sampling evenly across families to keep the subset diverse ----
    rng = random.Random(a.seed)
    if a.max_refs and len(in_order) > a.max_refs:
        by_family = defaultdict(list)
        for g in in_order:
            by_family[(by_genome.get(g) or {}).get("gvclass_family", "NA")].append(g)
        for v in by_family.values():
            rng.shuffle(v)
        picked, families = [], sorted(by_family)
        while len(picked) < a.max_refs and any(by_family[f] for f in families):
            for f in families:
                if by_family[f] and len(picked) < a.max_refs:
                    picked.append(by_family[f].pop())
        print(f"[refs] capped to {len(picked)} references across {len(families)} families "
              f"(--max-refs {a.max_refs})")
        in_order = sorted(picked)

    # ---- outgroups: a few genomes from each other order, to root the tree ----
    outgroups = []
    if wanted is not None and a.outgroup_per_order > 0:
        for order, gs in sorted(others_by_order.items()):
            if order in ("NA", ""):
                continue
            rng.shuffle(gs)
            outgroups.extend(gs[:a.outgroup_per_order])
        print(f"[refs] plus {len(outgroups)} outgroup references from "
              f"{len(others_by_order)} other orders")

    # ---- copy the .faa files ----
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    faa_dir = Path(a.faa_dir)
    index = {f.stem: f for f in faa_dir.rglob("*.faa")}
    print(f"[refs] {len(index)} .faa files available in {faa_dir}")

    copied, missing = 0, []
    with open(out.parent / "reference_manifest.tsv", "w") as man:
        man.write("taxon_file\tgenome\trole\torder\tfamily\tgenus\tlength_bp\tecosystem\n")
        for g, role in [(x, "in_order") for x in in_order] + \
                       [(x, "outgroup") for x in outgroups]:
            src = index.get(g)
            if src is None:
                missing.append(g); continue
            name = sanitise(g)
            shutil.copy(src, out / f"{name}.faa")
            r = by_genome.get(g, {})
            man.write("\t".join([name, g, role, r.get("Order", "NA"),
                                 r.get("gvclass_family", "NA"), r.get("new_genus", "NA"),
                                 r.get("LENbp", "NA"), r.get("Ecosystem", "NA")]) + "\n")
            copied += 1

    print(f"[refs] copied {copied} reference .faa files into {out}")
    if missing:
        print(f"[refs] WARNING: {len(missing)} selected genomes had no .faa in the bundle "
              f"(first few: {', '.join(missing[:5])})", file=sys.stderr)
    if copied == 0:
        sys.exit("No reference proteins were copied; check --faa-dir points at the "
                 "unpacked GVMAGs V2 bundle.")


if __name__ == "__main__":
    main()
