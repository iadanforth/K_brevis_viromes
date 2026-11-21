#!/usr/bin/env python3
"""Draw the pipeline diagram: what was run, with what, and what came out.

A pipeline figure is only worth having if it carries the real numbers, so every count
here is the measured result of the run documented in README.md, not a placeholder:
20 samples -> 848 bins -> 26 GVMAGs -> 19 vOTUs -> a 1,564-taxon tree.

Run (matplotlib lives in aniclustermap-env):
    module load miniforge3; eval "$(conda shell.bash hook)"; conda activate aniclustermap-env
    python3 -I scripts/plot_pipeline_figure.py --out-dir output
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# Okabe-Ito, the same colour-blind-safe palette used for the abundance figures, so the
# project's figures read as one set.
BLUE, ORANGE, GREEN, PINK, SKY, VERM, GREY = (
    "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#5A5A5A")
INK, MUTED, PAPER = "#1A1A1A", "#5A5A5A", "#FFFFFF"

# step, title, tool, detail lines, accent colour
STEPS = [
    ("01", "Download reads", "sra-tools 3.4.1",
     ["20 metagenomes (PRJNA1079807 + 2)", "9 incubation  ·  11 environmental", "2.0 TB"], GREY),
    ("02", "Quality trim", "fastp 1.3.7",
     ["trim 8 bp front / 20 bp tail", "Q20, ≥100 bp, no Ns", "1.6 TB"], GREY),
    ("03", "Reference genomes", "NCBI datasets",
     ["human, mouse, dog, cat,", "PhiX174, K. brevis host", "79 GB"], GREY),
    ("04", "Remove host reads", "bowtie2 2.5.5",
     ["map to the combined index,", "keep what does not match", "1.1 TB"], GREY),
    ("05", "Assemble, twice", "SPAdes 4.3.0",
     ["metaviralSPAdes  →  90 GB", "metaSPAdes (18/20) →  322 GB", "two assemblers per sample"], BLUE),
    ("06", "Merge assemblies", "seqkit + CheckV",
     ["dereplicate the two sets", "per sample", "6.6 GB merged contigs"], BLUE),
    ("08", "Bin and classify", "bbmap · MetaBAT2 · GVClass\nCheckV 1.1.1 · geNomad",
     ["848 bins from 20 samples", "→ 26 GVMAGs in 9 samples", "5 high-quality (CheckV)"], GREEN),
    ("09", "Cluster into vOTUs", "fastANI 1.33",
     ["≥95% ANI / ≥85% AF (MIUViG)", "26 GVMAGs → 19 vOTUs", "6 found in >1 sample"], GREEN),
    ("10", "Marker alignment", "ncldv_markersearch\nClustal Ω · trimAl",
     ["7 giant-virus marker genes", "+1,523 GVMAGs V2 references", "1,564 taxa × 37,739 columns"], PINK),
    ("11", "Phylogeny", "IQ-TREE 2.3.6 · RAxML",
     ["WAG+F+I+R4, 1000 UFBoot", "64% of branches ≥95", "tree certainty 0.73"], PINK),
]

# The three questions the tree and tables were built to answer.
RESULTS = [
    ("Abundance through time", ORANGE,
     ["CoverM: 11 environmental samples", "mapped onto the 19 vOTUs",
      "", "7 detections, all of them vOTUs",
      "assembled from that same sample.", "None of the 12 incubation vOTUs",
      "is detectable in the environment."]),
    ("Match to the published study", VERM,
     ["fastANI vs Booker et al. 2026", "(GenBank PV483752-6, PV637767-835)",
      "", "12 of 26 GVMAGs are the same",
      "species at 99.4-100% ANI.", "All 6 published vOTUs were",
      "independently recovered."]),
    ("Novelty", SKY,
     ["vs GVMAGs V2 (8,508 genomes)", "",
      "12 of 19 vOTUs have a relative", "above ~80% ANI, 7 have none.",
      "", "Not one reaches the species", "threshold: all 19 are new."]),
]


def box(ax, x, y, w, h, fc, ec, lw=1.2, r=0.012, z=2, alpha=1.0):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                       linewidth=lw, edgecolor=ec, facecolor=fc, zorder=z, alpha=alpha)
    ax.add_patch(p)
    return p


def arrow(ax, x1, y1, x2, y2, color=MUTED, lw=1.4, style="-|>", rad=0.0):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=11, linewidth=lw, color=color,
                                 connectionstyle=f"arc3,rad={rad}", zorder=1))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out-dir", default="output")
    p.add_argument("--basename", default="pipeline_overview")
    a = p.parse_args()

    fig, ax = plt.subplots(figsize=(15.5, 11.4))
    fig.patch.set_facecolor(PAPER)
    ax.set_facecolor(PAPER)
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

    ax.text(1.5, 97.6, "Karenia brevis bloom virome: giant viruses from metagenome to phylogeny",
            fontsize=16.5, fontweight="bold", color=INK, va="top")
    ax.text(1.5, 94.2,
            "Re-analysis and extension of Booker et al. 2026 (ISME Communications).  "
            "Every count is the measured result of this run; see README.md for the step that produced it.",
            fontsize=9.4, color=MUTED, va="top")

    ax.text(2.0, 90.0, "PROCESSING  ·  steps 01-11", fontsize=9.2, fontweight="bold",
            color="#8A8A8A", va="bottom")
    ax.text(48.5, 90.0, "FINDINGS  ·  steps 12-15", fontsize=9.2, fontweight="bold",
            color="#8A8A8A", va="bottom")

    # ---- left column: the sequential pipeline -----------------------------------------
    BX, BW, BH, GAP = 2.0, 40.0, 7.6, 1.15
    top = 88.2
    centres = []
    for i, (num, title, tool, lines, accent) in enumerate(STEPS):
        y = top - i * (BH + GAP) - BH
        box(ax, BX, y, BW, BH, "#FFFFFF", "#D8D8D8")
        box(ax, BX, y, 0.85, BH, accent, accent)           # accent spine
        ax.text(BX + 2.1, y + BH - 1.5, f"{num}", fontsize=11.5, fontweight="bold",
                color=accent, va="top")
        ax.text(BX + 6.0, y + BH - 1.4, title, fontsize=11.2, fontweight="bold",
                color=INK, va="top")
        ax.text(BX + 6.0, y + BH - 3.5, tool, fontsize=8.2, color=accent,
                va="top", style="italic")
        ax.text(BX + 21.0, y + BH - 1.4, "\n".join(lines), fontsize=8.3, color=MUTED,
                va="top", linespacing=1.5)
        centres.append((y, y + BH))
        if i:
            arrow(ax, BX + BW / 2, centres[i - 1][0], BX + BW / 2, y + BH)

    bot = centres[-1][0]

    # ---- right column: the three questions ---------------------------------------------
    RX, RW, RH = 48.5, 49.5, 20.5
    rtop = 88.2
    for i, (title, accent, lines) in enumerate(RESULTS):
        y = rtop - i * (RH + 2.6) - RH
        box(ax, RX, y, RW, RH, "#FCFCFC", "#D8D8D8")
        box(ax, RX, y + RH - 0.85, RW, 0.85, accent, accent)   # accent header bar
        ax.text(RX + 2.0, y + RH - 2.6, title, fontsize=11.6, fontweight="bold",
                color=INK, va="top")
        ax.text(RX + 2.0, y + RH - 5.6, "\n".join(lines), fontsize=9.0, color=MUTED,
                va="top", linespacing=1.62)


    # ---- caveat strip -------------------------------------------------------------------
    cy = rtop - 3 * (RH + 2.6) - 2.0
    box(ax, RX, cy - 16.6, RW, 16.4, "#FFF8EC", "#E8D5AE")
    ax.text(RX + 2.0, cy - 1.5, "Read alongside the result", fontsize=10.2,
            fontweight="bold", color="#8A6412", va="top")
    ax.text(RX + 2.0, cy - 4.0,
            "Alignment identity averages 14.5%, so deep branches are\n"
            "weakly supported — read support, not topology alone.\n"
            "Non-detection means <1.4× coverage in that sample,\n"
            "not absence from the Gulf of Mexico.\n"
            "vOTU05 is called 100% complete yet has no marker genes,\n"
            "no GVMAGs V2 relative and no published match.",
            fontsize=8.4, color="#6B4E0E", va="top", linespacing=1.5)

    foot = min(bot, cy - 16.6) - 3.4
    ax.axhline(y=foot + 1.6, xmin=0.012, xmax=0.988, color="#E4E4E4", linewidth=0.9)
    ax.text(BX, foot - 0.4,
            "Steps 00 and 07 build the conda environments and reference databases; "
            "12-15 collect results and run the three analyses on the right.",
            fontsize=8.1, color=MUTED, va="top")
    ax.text(BX, foot - 3.2,
            "Scripts, parameters and every command: github.com/iadanforth/K_brevis_viromes   ·   "
            "one-command rerun:  bash run_pipeline.sh all",
            fontsize=8.0, color="#9A9A9A", va="top")
    ax.set_ylim(foot - 6.0, 100)

    import os
    os.makedirs(a.out_dir, exist_ok=True)
    for ext in ("pdf", "png"):
        f = os.path.join(a.out_dir, f"{a.basename}.{ext}")
        fig.savefig(f, dpi=240 if ext == "png" else None,
                    bbox_inches="tight", facecolor=PAPER)
        print(f"wrote {f}")


if __name__ == "__main__":
    main()
