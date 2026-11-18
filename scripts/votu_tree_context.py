#!/usr/bin/env python3
"""Where does each K. brevis vOTU sit in the step 11 tree, and how well supported is it?

A 1,564-tip Newick file does not answer "what is my virus related to?" by inspection. For
each of our vOTU tips this walks up to its parent node and reports:

  * the UFBoot support on the branch that places it (>=95 is well supported; UFBoot values
    are NOT comparable to classical bootstrap, where >=70 is the usual threshold)
  * its sister group - the other tips descending from the same parent - and, where those
    are GVMAGs V2 references, their viral order and family from reference_manifest.tsv

Reading the output: a vOTU whose sister group is a tight set of one family, on a branch
with high support, is confidently placed in that family. A vOTU on a low-support branch
with a taxonomically mixed sister group is not placed by this tree, however long the
sister list is.

Usage:
  python3 -I scripts/votu_tree_context.py \
      --tree output/tree/kbrevis_gvmags.treefile \
      --manifest ~/scr10/data/kbrevis/phylogeny/reference_manifest.tsv \
      --out output/tree/votu_tree_context.tsv
"""
import argparse
import csv
import re
from collections import Counter
from pathlib import Path


class Node:
    __slots__ = ("name", "support", "parent", "children")

    def __init__(self):
        self.name = None
        self.support = None
        self.parent = None
        self.children = []

    def tips(self):
        if not self.children:
            return [self.name]
        out = []
        for c in self.children:
            out.extend(c.tips())
        return out


def parse_newick(text):
    """Minimal Newick reader. Internal node labels are treated as support values."""
    text = text.strip().rstrip(";")
    root = Node()
    cur = root
    i = 0
    n = len(text)
    token = ""

    def finish_label(node, tok):
        tok = tok.strip()
        if not tok:
            return
        # label[:branchlength]
        if ":" in tok:
            lab, _, _bl = tok.partition(":")
        else:
            lab = tok
        lab = lab.strip()
        if not lab:
            return
        if node.children:            # internal -> support
            try:
                node.support = float(lab)
            except ValueError:
                node.support = None
        else:                        # tip -> name
            node.name = lab

    while i < n:
        ch = text[i]
        if ch == "(":
            child = Node()
            child.parent = cur
            cur.children.append(child)
            cur = child
            token = ""
        elif ch == ",":
            finish_label(cur, token)
            token = ""
            sib = Node()
            sib.parent = cur.parent
            cur.parent.children.append(sib)
            cur = sib
        elif ch == ")":
            finish_label(cur, token)
            token = ""
            cur = cur.parent
        else:
            token += ch
        i += 1
    finish_label(cur, token)
    return root


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tree", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--max-sisters", type=int, default=8)
    a = p.parse_args()

    root = parse_newick(Path(a.tree).read_text())

    meta = {}
    with open(a.manifest) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            meta[row["taxon_file"].replace("_", ".")] = row

    # index tips
    tips = {}

    def walk(nd):
        if not nd.children:
            if nd.name:
                tips[nd.name] = nd
        for c in nd.children:
            walk(c)

    walk(root)
    ours = sorted(t for t in tips if t.startswith("vOTU"))
    print(f"{len(tips)} tips, {len(ours)} ours\n")

    rows = []
    for name in ours:
        nd = tips[name]
        parent = nd.parent
        # climb until the parent has a sister group that is not just this tip
        sisters = []
        anc = parent
        while anc is not None:
            sisters = [t for t in anc.tips() if t != name]
            if sisters:
                break
            anc = anc.parent
        support = anc.support if anc is not None else None

        orders = Counter()
        families = Counter()
        for s in sisters:
            m = meta.get(s)
            if m:
                orders[m.get("order") or "NA"] += 1
                families[m.get("family") or "NA"] += 1
            elif s.startswith("vOTU"):
                orders["(our vOTU)"] += 1
                families["(our vOTU)"] += 1

        def top(c):
            return ", ".join(f"{k} x{v}" if v > 1 else k for k, v in c.most_common(3))

        rows.append({
            "votu": name.split(".")[0],
            "tip": name,
            "ufboot_support_of_placement": "NA" if support is None else f"{support:.0f}",
            "n_sister_tips": len(sisters),
            "sister_orders": top(orders) or "NA",
            "sister_families": top(families) or "NA",
            "sister_examples": ", ".join(sisters[:a.max_sisters]),
        })

    cols = ["votu", "tip", "ufboot_support_of_placement", "n_sister_tips",
            "sister_orders", "sister_families", "sister_examples"]
    with open(a.out, "w") as out:
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    print(f"{'vOTU':8s} {'supp':>5s} {'sis':>5s}  {'sister orders':38s} sister families")
    for r in sorted(rows, key=lambda x: x["votu"]):
        print(f"{r['votu']:8s} {r['ufboot_support_of_placement']:>5s} "
              f"{r['n_sister_tips']:>5d}  {r['sister_orders'][:38]:38s} "
              f"{r['sister_families'][:40]}")
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
