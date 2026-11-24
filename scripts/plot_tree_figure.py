#!/usr/bin/env python3
"""Draw the giant-virus phylogeny as a publication figure, from the IQ-TREE output.

WHY THIS SCRIPT EXISTS
----------------------
`output/Kbrevis_virome_tree.png` was originally exported by hand from iTOL, which meant it
could not be regenerated, and its legend text and vOTU branch colouring were fixed at
whatever iTOL's defaults produced - about 4 pt type and hairline red branches, both too
weak for publication. This script redraws the same tree from
`output/tree/kbrevis_gvmags.treefile` with full control over type size and line weight, so
the figure is reproducible and can be tuned without a round trip through a web tool.

WHAT IS DRAWN
-------------
A circular phylogram of all 1,564 taxa in the alignment: 1,547 reference genomes from
GVMAGs V2 plus the 17 K. brevis vOTUs that carried at least one marker gene. Reference
branches are thin and grey, because they are context rather than result. The 17 vOTU
branches are drawn thick, in saturated red, over the top of everything else, with a filled
wedge at the tip and a label outside the ring - they are the point of the figure and should
be findable at a glance. An annulus outside the tips carries each taxon's viral ORDER, the
taxonomic rank at which our genomes can actually be placed.

WHY 17 AND NOT 19
-----------------
vOTU05 and vOTU19 carry none of the seven core Nucleocytoviricota marker genes, so
ncldv_markersearch dropped them before alignment and they are not in the tree at all. See
README section 10j; this script re-derives the count from the tree rather than hard-coding
it, and prints which vOTUs are missing.

INPUTS
  output/tree/kbrevis_gvmags.treefile      IQ-TREE maximum-likelihood tree, Newick
  output/tree/itol_order_colorstrip.txt    taxon -> colour, order  (written by
                                           make_itol_annotations.py)
OUTPUT
  output/Kbrevis_virome_tree.{png,pdf}

Run (matplotlib lives in aniclustermap-env):
    module load miniforge3; eval "$(conda shell.bash hook)"; conda activate aniclustermap-env
    python3 -I scripts/plot_tree_figure.py
"""
import argparse
import math
import re
from collections import namedtuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Wedge, Patch

VOTU_RED = "#D7191C"      # saturated, prints well, distinct from every order colour
REF_GREY = "#9A9A9A"
INK, MUTED, PAPER = "#1A1A1A", "#4A4A4A", "#FFFFFF"


class Node:
    __slots__ = ("name", "length", "children", "parent", "angle", "radius", "is_votu")

    def __init__(self, name="", length=0.0):
        self.name, self.length = name, length
        self.children, self.parent = [], None
        self.angle = self.radius = 0.0
        self.is_votu = False


def parse_newick(text):
    """Minimal Newick reader. Internal node labels are IQ-TREE's UFBoot support values and
    are kept but not used here. Written out rather than taking an ete3/Bio.Phylo
    dependency, because neither is installed in this project's environments and the format
    needed here is small and well defined."""
    text = text.strip().rstrip(";")
    root = Node()
    cur = root
    i, n = 0, len(text)
    token = re.compile(r"[^,():]*")
    while i < n:
        c = text[i]
        if c == "(":
            child = Node(); child.parent = cur; cur.children.append(child)
            cur = child; i += 1
        elif c == ",":
            sib = Node(); sib.parent = cur.parent; cur.parent.children.append(sib)
            cur = sib; i += 1
        elif c == ")":
            cur = cur.parent; i += 1
        elif c == ":":
            i += 1
            m = token.match(text, i)
            cur.length = float(m.group() or 0.0)
            i = m.end()
        else:
            m = token.match(text, i)
            cur.name = m.group().strip()
            i = m.end() if m.end() > i else i + 1
    return root


def layout(root, gap_frac=0.015):
    """Radial coordinates. Leaves are spread evenly over the circle, leaving a small gap so
    the first and last tip do not touch; an internal node sits at the mean angle of its
    children. Radius is the cumulative branch length from the root, i.e. a true phylogram -
    distances are meaningful, not just topology."""
    leaves, stack = [], [root]
    order = []
    while stack:                      # iterative, because 1,564 taxa can nest deeply
        nd = stack.pop()
        order.append(nd)
        stack.extend(nd.children)
    for nd in order:                  # root -> tips: radius accumulates
        nd.radius = (nd.parent.radius + nd.length) if nd.parent else 0.0
    for nd in order:
        if not nd.children:
            leaves.append(nd)
    span = 2 * math.pi * (1 - gap_frac)
    start = math.pi / 2              # 12 o'clock, so the gap sits at the top
    for k, lf in enumerate(leaves):
        lf.angle = start - span * k / max(len(leaves) - 1, 1)
    for nd in reversed(order):        # tips -> root: parents average their children
        if nd.children:
            nd.angle = sum(c.angle for c in nd.children) / len(nd.children)
    return leaves


def arc_points(r, a0, a1, step=math.radians(0.6)):
    k = max(int(abs(a1 - a0) / step) + 1, 2)
    return [(r * math.cos(a0 + (a1 - a0) * j / (k - 1)),
             r * math.sin(a0 + (a1 - a0) * j / (k - 1))) for j in range(k)]


def norm(name):
    """iTOL annotation files and the treefile spell the same taxon differently: IQ-TREE
    rewrites '_' and other punctuation as '.'. Collapsing every run of non-alphanumerics to
    a single '.' and dropping leading/trailing ones makes the two agree."""
    return re.sub(r"[^A-Za-z0-9]+", ".", name).strip(".").lower()


def read_orders(path):
    colour, label = {}, {}
    legend = []
    with open(path) as fh:
        in_data = False
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("LEGEND_COLORS"):
                legend_colours = line.split("\t")[1:]
            elif line.startswith("LEGEND_LABELS"):
                legend_labels = line.split("\t")[1:]
                legend = list(zip(legend_labels, legend_colours))
            elif line.strip() == "DATA":
                in_data = True
            elif in_data and line.strip():
                parts = line.split("\t")
                if len(parts) >= 2:
                    colour[norm(parts[0])] = parts[1]
                    label[norm(parts[0])] = parts[2] if len(parts) > 2 else ""
    return colour, label, legend


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", default="output/tree/kbrevis_gvmags.treefile")
    ap.add_argument("--orders", default="output/tree/itol_order_colorstrip.txt")
    ap.add_argument("--out", default="output/Kbrevis_virome_tree")
    ap.add_argument("--dpi", type=int, default=400)
    a = ap.parse_args()

    root = parse_newick(open(a.tree).read())
    leaves = layout(root)
    order_colour, order_label, legend = read_orders(a.orders)

    votu_re = re.compile(r"(vOTU\d+)", re.I)
    for lf in leaves:
        lf.is_votu = bool(votu_re.search(lf.name))
    votu_leaves = [lf for lf in leaves if lf.is_votu]
    found = sorted({votu_re.search(lf.name).group(1) for lf in votu_leaves})
    print(f"{len(leaves)} taxa in the tree; {len(votu_leaves)} are ours, "
          f"covering {len(found)} vOTUs:")
    print("  " + " ".join(found))

    # Red marks each vOTU's OWN terminal branch, not its lineage back to the root. Tracing
    # ancestors was tried first and is wrong for this figure: the 17 vOTUs sit in many
    # different clades, so their paths merge low down and turn most of the backbone red,
    # which hides the thing the colour is supposed to pick out. Each tip instead gets a
    # radial leader line out to the ring, which is what actually makes it findable.
    on_path = {id(lf) for lf in votu_leaves}

    rmax = max(lf.radius for lf in leaves)

    # Build the branch segments. Each branch is an L in polar space: an arc at the parent's
    # radius spanning its children, then a radial line out to the child.
    ref_seg, votu_seg = [], []
    stack = [root]
    while stack:
        nd = stack.pop()
        if nd.children:
            angs = [c.angle for c in nd.children]
            # the connecting arc, split so the red part covers only vOTU-bearing children
            for c in nd.children:
                ref_seg.append(arc_points(nd.radius, nd.angle, c.angle))
            stack.extend(nd.children)
        if nd.parent is not None:
            x0, y0 = nd.parent.radius * math.cos(nd.angle), nd.parent.radius * math.sin(nd.angle)
            x1, y1 = nd.radius * math.cos(nd.angle), nd.radius * math.sin(nd.angle)
            sel = votu_seg if id(nd) in on_path else ref_seg
            sel.append([(x0, y0), (x1, y1)])

    fig = plt.figure(figsize=(15, 15), dpi=a.dpi, facecolor=PAPER)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_aspect("equal"); ax.set_axis_off()
    ax.set_facecolor(PAPER)

    ax.add_collection(LineCollection(ref_seg, colors=REF_GREY, linewidths=0.35,
                                     capstyle="round", zorder=2))
    # The vOTU path is drawn twice: a white casing underneath so it stays legible where it
    # crosses dense reference branches, then the red line on top.
    ax.add_collection(LineCollection(votu_seg, colors="white", linewidths=3.4,
                                     capstyle="round", zorder=3))
    ax.add_collection(LineCollection(votu_seg, colors=VOTU_RED, linewidths=1.9,
                                     capstyle="round", zorder=4))

    # A thin red leader from each vOTU tip out to the ring. The terminal branches alone are
    # short and buried among 1,547 grey neighbours; the leader is what carries the eye from
    # the label in to the branch.
    r0_lead = rmax * 1.015
    ax.add_collection(LineCollection(
        [[(lf.radius * math.cos(lf.angle), lf.radius * math.sin(lf.angle)),
          (r0_lead * math.cos(lf.angle), r0_lead * math.sin(lf.angle))]
         for lf in votu_leaves],
        colors=VOTU_RED, linewidths=1.1, alpha=0.85, zorder=4))

    # Order annulus just outside the tips.
    r0 = rmax * 1.015
    dr = rmax * 0.045
    half = math.radians(360 / len(leaves) / 2) * 1.08
    wedges = [Wedge((0, 0), r0 + dr, math.degrees(lf.angle - half),
                    math.degrees(lf.angle + half), width=dr) for lf in leaves]
    cols = [order_colour.get(norm(lf.name), "#DDDDDD") for lf in leaves]
    ax.add_collection(PatchCollection(wedges, facecolors=cols, edgecolors="none", zorder=5))

    # A second, taller wedge marks each vOTU so it is visible in the ring as well as in the
    # branches.
    rv = r0 + dr
    dv = rmax * 0.055
    vw = [Wedge((0, 0), rv + dv, math.degrees(lf.angle - half * 2.6),
                math.degrees(lf.angle + half * 2.6), width=dv) for lf in votu_leaves]
    ax.add_collection(PatchCollection(vw, facecolors=VOTU_RED, edgecolors="none", zorder=6))

    # vOTU labels outside the ring, rotated to read outward. Several vOTUs are near
    # neighbours in the tree, so their labels would print on top of each other at the tip
    # angle; nudge them apart to a minimum angular separation, keeping tree order. The
    # leader line still points at the true position, so nothing is misplaced.
    MIN_SEP = math.radians(3.6)
    ordered = sorted(votu_leaves, key=lambda l: l.angle)
    label_angle = {id(l): l.angle for l in ordered}
    for _ in range(200):                 # a few relaxation sweeps settle it
        moved = False
        for i in range(len(ordered) - 1):
            lo, hi = ordered[i], ordered[i + 1]
            gap = label_angle[id(hi)] - label_angle[id(lo)]
            if gap < MIN_SEP:
                push = (MIN_SEP - gap) / 2
                label_angle[id(lo)] -= push
                label_angle[id(hi)] += push
                moved = True
        if not moved:
            break

    rl = rv + dv + rmax * 0.025
    for lf in votu_leaves:
        name = votu_re.search(lf.name).group(1)
        ang = label_angle[id(lf)]
        deg = math.degrees(ang) % 360
        flip = 90 < deg < 270
        ax.text(rl * math.cos(ang), rl * math.sin(ang), name,
                rotation=deg + 180 if flip else deg,
                rotation_mode="anchor", ha="right" if flip else "left", va="center",
                fontsize=13, fontweight="bold", color=VOTU_RED, zorder=7)

    lim = rl + rmax * 0.20
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)

    # Legend: viral order, plus the vOTU key. Set at a size that stays readable when the
    # figure is reduced to a single journal column.
    used = {order_label.get(norm(lf.name), "") for lf in leaves}
    handles = [Patch(facecolor=c, edgecolor="none", label=l)
               for l, c in legend if l in used and "vOTU" not in l]
    handles.append(Patch(facecolor=VOTU_RED, edgecolor="none",
                         label=f"K. brevis vOTU ({len(found)} of 19)"))
    # framealpha=1: matplotlib's default 0.8 let the colour ring show through the legend
    # box where the two overlap, which reads as a printing fault.
    leg = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.005, 0.995),
                    frameon=True, framealpha=1.0, facecolor=PAPER, fontsize=15,
                    title="Viral order", title_fontsize=17, labelspacing=0.55,
                    handlelength=1.6, borderpad=0.9, handleheight=1.1)
    leg.get_frame().set_edgecolor("#CCCCCC"); leg.get_frame().set_linewidth(1.0)
    leg.get_title().set_fontweight("bold"); leg.get_title().set_color(INK)
    for t in leg.get_texts():
        t.set_color(INK)

    # Scale bar: one unit is one substitution per site.
    bar = 0.5
    bx, by = -lim * 0.86, -lim * 0.90
    ax.plot([bx, bx + bar], [by, by], color=INK, linewidth=2.6, zorder=8,
            solid_capstyle="butt")
    for dx in (0, bar):                  # end ticks, so the bar's extent is unambiguous
        ax.plot([bx + dx, bx + dx], [by - lim * 0.012, by + lim * 0.012],
                color=INK, linewidth=2.0, zorder=8)
    ax.text(bx + bar / 2, by + lim * 0.022, f"{bar} substitutions per site",
            ha="center", va="bottom", fontsize=14, color=INK)

    for ext in ("png", "pdf"):
        out = f"{a.out}.{ext}"
        fig.savefig(out, dpi=a.dpi, facecolor=PAPER)
        print(f"wrote {out}")
    plt.close(fig)

    missing = sorted(set(f"vOTU{i:02d}" for i in range(1, 20)) - set(found))
    if missing:
        print(f"\nNOT in the tree: {', '.join(missing)} - they carried none of the seven "
              f"marker genes, so ncldv_markersearch dropped them (README 10j).")


if __name__ == "__main__":
    main()
