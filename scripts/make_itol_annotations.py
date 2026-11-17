#!/usr/bin/env python3
"""Build iTOL annotation files for the step 11 giant-virus phylogeny.

The tree has 1,564 tips, of which only 17 are ours. Uploaded bare to iTOL it is an
unreadable hairball: every tip is an anonymous accession like `IMGM3300034190.BIN56`,
and there is no way to see where the Karenia brevis vOTUs sit or which clade is which
viral order. These annotation files fix that. Upload the .treefile to iTOL first, then
drag each annotation file onto the rendered tree.

A naming detail that matters: ncldv_markersearch/Clustal Omega rewrote `_` as `.` in the
FASTA headers, so a genome called `IMGM3300034190_BIN56` in GVDB_stats.tsv and in
reference_manifest.tsv appears as `IMGM3300034190.BIN56` in the tree. Every lookup here
normalises `_` to `.` before matching, or nothing would join.

Outputs (all iTOL "dataset" text files):
  itol_labels.txt            readable tip labels for our vOTUs (references left alone)
  itol_tree_colors.txt       our vOTUs' branches and labels in red and bold
  itol_order_colorstrip.txt  a colour strip giving every tip's viral order
  itol_family_colorstrip.txt the same for family, restricted to the common families
  itol_popup.txt             hover text: completeness, nearest reference, ecosystem

Usage:
  python3 -I scripts/make_itol_annotations.py \
      --tree   ~/scr10/data/kbrevis/tree/kbrevis_gvmags.treefile \
      --manifest ~/scr10/data/kbrevis/phylogeny/reference_manifest.tsv \
      --votus  ~/K_brevis_viromes/output/votus/gvmags_votus.tsv \
      --nearest ~/K_brevis_viromes/output/votus/votu_nearest_reference.tsv \
      --out-dir ~/K_brevis_viromes/output/tree
"""
import argparse
import csv
import re
from collections import Counter
from pathlib import Path

OURS_COLOR = "#d62728"

# Colour-blind-safe qualitative palette (Okabe-Ito, plus greys for the long tail).
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9",
           "#D55E00", "#F0E442", "#999999", "#7F7F7F", "#BCBD22"]


def norm(name):
    """Tree tips carry '.' where the source tables carry '_'."""
    return name.replace("_", ".")


def read_tips(tree_path):
    t = Path(tree_path).read_text()
    return re.findall(r"[(,]([A-Za-z0-9._\-]+):", t)


def header(kind, label, color, extra=()):
    out = [kind, "SEPARATOR TAB", f"DATASET_LABEL\t{label}", f"COLOR\t{color}"]
    out.extend(extra)
    out.append("DATA")
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tree", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--votus")
    p.add_argument("--nearest")
    p.add_argument("--out-dir", required=True)
    a = p.parse_args()

    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    tips = read_tips(a.tree)
    ours = [t for t in tips if t.startswith("vOTU")]
    print(f"{len(tips)} tips, {len(ours)} of them ours")

    # ---- reference metadata, keyed by the tree's spelling ----
    meta = {}
    with open(a.manifest) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            meta[norm(row["taxon_file"])] = row
    matched = sum(1 for t in tips if t in meta)
    print(f"manifest rows matched to tips: {matched} / {len(tips) - len(ours)} references")

    # ---- our vOTUs' own details ----
    nearest = {}
    if a.nearest and Path(a.nearest).exists():
        with open(a.nearest) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                nearest[row["votu"]] = row
    members = Counter()
    if a.votus and Path(a.votus).exists():
        with open(a.votus) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                members[row["votu"]] += 1

    def votu_of(tip):
        return tip.split(".")[0]

    def sample_of(tip):
        m = re.match(r"vOTU\d+\.\.([A-Za-z0-9]+)\.\.(.+?)\.?$", tip)
        return (m.group(1), m.group(2).replace(".", "_")) if m else ("", "")

    # ---- 1. readable labels for our tips ----
    lines = ["LABELS", "SEPARATOR TAB", "DATA"]
    for t in ours:
        v = votu_of(t)
        srr, binname = sample_of(t)
        n = members.get(v, 1)
        tag = f" [{n} genomes]" if n > 1 else ""
        lines.append(f"{t}\t{v} | {srr}{tag}")
    (out / "itol_labels.txt").write_text("\n".join(lines) + "\n")

    # ---- 2. highlight our tips ----
    lines = ["TREE_COLORS", "SEPARATOR TAB", "DATA"]
    for t in ours:
        lines.append(f"{t}\tlabel\t{OURS_COLOR}\tbold\t2")
        lines.append(f"{t}\tbranch\t{OURS_COLOR}\tnormal\t4")
    (out / "itol_tree_colors.txt").write_text("\n".join(lines) + "\n")

    # ---- 3 & 4. colour strips by order and by family ----
    def colorstrip(field, filename, title, our_value):
        values = Counter()
        for t in tips:
            if t in meta:
                values[meta[t].get(field) or "NA"] += 1
        common = [v for v, _ in values.most_common(len(PALETTE)) if v != "NA"]
        cmap = {v: PALETTE[i] for i, v in enumerate(common)}
        legend_shapes = "\t".join(["1"] * (len(common) + 1))
        legend_colors = "\t".join([cmap[v] for v in common] + [OURS_COLOR])
        legend_labels = "\t".join(common + ["K. brevis vOTU"])
        head = header("DATASET_COLORSTRIP", title, PALETTE[0], extra=[
            "COLOR_BRANCHES\t0",
            "STRIP_WIDTH\t30",
            "MARGIN\t5",
            f"LEGEND_TITLE\t{title}",
            f"LEGEND_SHAPES\t{legend_shapes}",
            f"LEGEND_COLORS\t{legend_colors}",
            f"LEGEND_LABELS\t{legend_labels}",
        ])
        body = []
        for t in tips:
            if t.startswith("vOTU"):
                body.append(f"{t}\t{OURS_COLOR}\t{our_value}")
            elif t in meta:
                v = meta[t].get(field) or "NA"
                if v in cmap:
                    body.append(f"{t}\t{cmap[v]}\t{v}")
        (out / filename).write_text("\n".join(head + body) + "\n")
        print(f"{filename}: {len(body)} tips coloured over {len(cmap)} {field} values")

    colorstrip("order", "itol_order_colorstrip.txt", "Viral order", "K. brevis vOTU")
    colorstrip("family", "itol_family_colorstrip.txt", "Family", "K. brevis vOTU")

    # ---- 5. hover pop-ups ----
    lines = ["POPUP_INFO", "SEPARATOR TAB", "DATA"]
    for t in ours:
        v = votu_of(t)
        srr, binname = sample_of(t)
        n = nearest.get(v, {})
        html = (f"<h3>{v}</h3>"
                f"<p><b>Karenia brevis virome GVMAG</b><br>"
                f"sample: {srr}<br>bin: {binname}<br>"
                f"genomes in vOTU: {members.get(v, 1)}<br>"
                f"CheckV completeness: {n.get('checkv_completeness', 'NA')}%<br>"
                f"length: {n.get('length_bp', 'NA')} bp<br>"
                f"nearest GVMAGs V2 relative: {n.get('nearest_reference', 'NA')} "
                f"({n.get('ani_percent', 'NA')}% ANI over {n.get('aligned_fraction_percent', 'NA')}%)<br>"
                f"reference family: {n.get('ref_family', 'NA')}<br>"
                f"reference ecosystem: {n.get('ref_ecosystem', 'NA')}</p>")
        lines.append(f"{t}\t{v}\t{html}")
    for t in tips:
        if t in meta:
            m = meta[t]
            html = (f"<h3>{m['genome']}</h3><p>GVMAGs V2 reference ({m['role']})<br>"
                    f"order: {m['order']}<br>family: {m['family']}<br>"
                    f"genus: {m['genus']}<br>length: {m['length_bp']} bp<br>"
                    f"ecosystem: {m['ecosystem']}</p>")
            lines.append(f"{t}\t{m['genome']}\t{html}")
    (out / "itol_popup.txt").write_text("\n".join(lines) + "\n")
    print(f"itol_popup.txt: {len(lines) - 3} entries")


if __name__ == "__main__":
    main()
