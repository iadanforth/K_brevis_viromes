#!/usr/bin/env python3
"""Draw the pipeline diagram: a plain vertical flowchart of what was run.

One box per step, top to bottom, connected by arrows. Each box holds WHAT WAS DONE and
THE PROGRAM THAT DID IT; the text to its right says briefly HOW. No title, no results -
this figure documents the method, not the findings (those live in README.md and in the
abundance / ANI / tree outputs).

Redrawn 2026-11-23 (prompt 61 in PROMPT_LOG.md) in this simpler style at the user's
request, replacing a wider three-column version that also summarised results.

Colours are Okabe-Ito, the same colour-blind-safe palette used by every other figure in
the project, so they read as one set. Colour here encodes the PHASE of the pipeline
(read preparation, assembly, genome recovery, phylogeny, quantification), which is the
only grouping the flowchart itself cannot show.

Run (matplotlib lives in aniclustermap-env):
    module load miniforge3; eval "$(conda shell.bash hook)"; conda activate aniclustermap-env
    python3 -I scripts/plot_pipeline_figure.py --out-dir output
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# Okabe-Ito colour-blind-safe palette, shared with the abundance figures.
BLUE, ORANGE, GREEN, PINK, SKY, VERM, GREY = (
    "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#5A5A5A")
INK, MUTED, PAPER = "#1A1A1A", "#4A4A4A", "#FFFFFF"


def tint(hex_colour, amount):
    """Mix a colour with white. Gives each box a pale fill of its phase colour so the
    outline stays legible and the dark text on top keeps full contrast."""
    h = hex_colour.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    mix = lambda c: int(round(c + (255 - c) * amount))
    return "#%02X%02X%02X" % (mix(r), mix(g), mix(b))


# (step number, what was done, program(s) used, how it was done, phase colour)
STEPS = [
    ("01", "Download sequencing reads", "sra-tools 3.4.1",
     ["20 paired-end metagenomes: 18 from BioProject PRJNA1079807",
      "plus 2 coastal Florida runs (9 incubation, 11 environmental)"], GREY),
    ("02", "Quality-trim and filter reads", "fastp 1.3.7",
     ["trim 8 bp from the front and 20 bp from the tail, drop bases",
      "below Q20, discard reads under 100 bp or containing Ns"], GREY),
    ("03", "Fetch contaminant genomes", "NCBI datasets 16.2.0",
     ["human, mouse, dog, cat, PhiX174 spike-in and the",
      "Karenia brevis host genome, built into one bowtie2 index"], GREY),
    ("04", "Remove host and contaminant reads", "bowtie2 2.5.5",
     ["map every sample against the combined index and keep only",
      "the read pairs that fail to align, so what remains is non-host"], GREY),
    ("05", "Assemble reads into contigs", "SPAdes 4.3.0",
     ["each sample assembled twice: metaviralSPAdes, tuned for viral",
      "genomes, and metaSPAdes, which recovers longer contigs"], BLUE),
    ("06", "Merge the two assemblies", "seqkit 2.8.2",
     ["combine both contig sets per sample and remove duplicates,",
      "keeping contigs ≥5 kb as candidate viral genome fragments"], BLUE),
    ("08", "Bin contigs into genomes\nand keep the giant viruses",
     "bbmap 35.85 · MetaBAT2 2.17\nCheckV 1.1.1 · geNomad 1.8.0 · GVClass 1.1.0",
     ["map reads back to get per-contig coverage, group contigs into",
      "bins by coverage and tetranucleotide signal, then keep only the",
      "bins GVClass assigns to the NCLDV giant viruses: 26 GVMAGs"], GREEN),
    ("09", "Cluster genomes into vOTUs", "fastANI 1.33",
     ["all-vs-all average nucleotide identity; genomes sharing ≥95% ANI",
      "across ≥85% of the shorter one are one vOTU (MIUViG standard)"], GREEN),
    ("10", "Build a marker-gene alignment",
     "ncldv_markersearch · Clustal Ω 1.2.4\ntrimAl 1.5.0",
     ["find 7 universal giant-virus marker genes in each genome and in",
      "the GVMAGs V2 reference set, align each gene, concatenate the",
      "alignments and trim poorly aligned columns"], PINK),
    ("11", "Infer the phylogeny", "IQ-TREE 2.3.6",
     ["ModelFinder picks the substitution model, then maximum-likelihood",
      "tree search with 1000 ultrafast bootstrap replicates for support"], PINK),
    ("14", "Quantify each genome in\nthe environmental samples", "CoverM 0.7.0",
     ["map each environmental sample's reads onto every GVMAG and report",
      "coverage, breadth and relative abundance, run both at CoverM's",
      "defaults and with a 95% read-identity filter"], ORANGE),
    ("15", "Compare with the published genomes", "fastANI 1.33",
     ["ANI between our GVMAGs and those of Booker et al. 2026 from",
      "GenBank, to see which of their viruses we independently recovered"], VERM),
]

# The step number sits in a gutter to the LEFT of the box, not inside it: the longest
# titles ("Remove host and contaminant reads") run the full width of the box, and a number
# in the top-left corner collided with them.
NUM_X = 0.30
BOX_X, BOX_W = 0.62, 3.58          # inches
TEXT_X = BOX_X + BOX_W + 0.45
GAP = 0.40                          # vertical space for the arrow between boxes
PAD_Y = 0.22
FIG_W = 9.45


def box_height(step):
    """Boxes grow with their contents so nothing is cramped: one line for the step's
    description (sometimes two) plus one line per line of the program name."""
    _, what, tool, _, _ = step
    return 0.30 + 0.215 * what.count("\n") + 0.215 * (tool.count("\n") + 1) + 0.30


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="output")
    ap.add_argument("--stem", default="pipeline_overview")
    a = ap.parse_args()

    heights = [box_height(s) for s in STEPS]
    fig_h = sum(heights) + GAP * (len(STEPS) - 1) + 2 * PAD_Y

    fig = plt.figure(figsize=(FIG_W, fig_h), dpi=200, facecolor=PAPER)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, FIG_W)
    ax.set_ylim(0, fig_h)
    ax.set_axis_off()
    ax.set_facecolor(PAPER)

    y = fig_h - PAD_Y
    for i, (step, (num, what, tool, how, colour)) in enumerate(zip(heights, STEPS)):
        top, bot = y, y - step
        ax.add_patch(FancyBboxPatch(
            (BOX_X, bot), BOX_W, step,
            boxstyle="round,pad=0,rounding_size=0.07",
            linewidth=1.6, edgecolor=colour, facecolor=tint(colour, 0.88), zorder=2))

        # Step number, set in the phase colour in the left gutter.
        ax.text(NUM_X, (top + bot) / 2, num, ha="center", va="center",
                fontsize=10, fontweight="bold", color=colour, zorder=3)

        # What was done (bold) and the program that did it (smaller, in colour),
        # stacked and centred in the box.
        n_what, n_tool = what.count("\n") + 1, tool.count("\n") + 1
        block = 0.215 * n_what + 0.06 + 0.195 * n_tool
        cy = (top + bot) / 2
        ax.text(BOX_X + BOX_W / 2, cy + block / 2 - 0.215 * n_what / 2, what,
                ha="center", va="center", fontsize=10.5, fontweight="bold",
                color=INK, linespacing=1.35, zorder=3)
        ax.text(BOX_X + BOX_W / 2, cy - block / 2 + 0.195 * n_tool / 2, tool,
                ha="center", va="center", fontsize=8.6, color=colour,
                linespacing=1.35, zorder=3)

        # How it was done, to the right of the box.
        ax.text(TEXT_X, cy, "\n".join(how), ha="left", va="center",
                fontsize=9.2, color=MUTED, linespacing=1.55, zorder=3)

        if i < len(STEPS) - 1:
            ax.add_patch(FancyArrowPatch(
                (BOX_X + BOX_W / 2, bot), (BOX_X + BOX_W / 2, bot - GAP + 0.015),
                arrowstyle="-|>", mutation_scale=15, linewidth=1.5,
                color=INK, shrinkA=0, shrinkB=0, zorder=1))
        y = bot - GAP

    for ext in ("png", "pdf"):
        out = f"{a.out_dir}/{a.stem}.{ext}"
        fig.savefig(out, facecolor=PAPER, bbox_inches=None)
        print(f"wrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
