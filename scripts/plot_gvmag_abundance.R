#!/usr/bin/env Rscript
# NOT RUN BY DEFAULT (2026-11-24, prompt 64). Its outputs,
# output/abundance/gvmag_abundance_timeseries.{pdf,png}, were deleted at the user's
# request and the call was removed from run_pipeline.sh and from the README's Commands
# section. The script is kept because it is the only thing that draws the lenient-vs-strict
# and abundance-vs-breadth comparison, and because the data it reads
# (output/abundance/gvmag_abundance.csv) is still produced by step 14. Run it by hand if
# that comparison is wanted again:
#
#   Rscript scripts/plot_gvmag_abundance.R output/abundance/gvmag_abundance.csv output/abundance
#
# The headline figure is scripts/plot_incubation_votu_timeseries.R (step 18).
# GVMAG abundance and breadth through time across the environmental samples.
#
# Replaces plot_votu_abundance.R (2026-11-23, prompt 60 in PROMPT_LOG.md), which plotted
# 19 vOTU representatives under a single filter setting. Step 14 now maps every
# environmental sample against all 26 individual GVMAGs, following Booker et al. 2026
# ("...reads recruited to each GVMAG..."), and does it twice - at CoverM's defaults and
# with a 95% read-identity filter - so this figure shows both.
#
# Input : output/abundance/gvmag_abundance.csv  (written by merge_gvmag_abundance.py).
#         Self-contained: this script needs nothing else, so the figure can be rebuilt
#         from one small CSV.
# Output: gvmag_abundance_timeseries.{pdf,png}
#
# WHY TWO ROWS OF PANELS, NOT ONE
# -------------------------------
# The top row is RELATIVE ABUNDANCE - the fraction of the sample's reads that landed on
# each genome. That is the quantity Booker et al. report. The bottom row is BREADTH, the
# fraction of the genome covered by at least one read. They have to be read together:
# abundance alone cannot tell a real infection from reads borrowed off a relative's
# conserved genes, because a genome can soak up reads at high apparent depth while only a
# few percent of it is covered at all. Breadth is what separates the two, so it is plotted
# directly underneath rather than buried in a table. The dashed line marks the 75% breadth
# rule (MIUViG, Roux et al. 2019) used to call a genome present.
#
# WHY TWO COLUMNS
# ---------------
# Left: CoverM's defaults, which apply NO read-identity filter - the paper specifies no
# settings, so this is the closest reading of what they did. Right: reads must match at
# >=95% identity over >=75% of their length, the MIUViG species boundary and the same
# cutoff that defined the vOTUs in step 09. The difference between the columns IS the
# result: recruitment that vanishes under the identity filter was never this genome's.
#
# A METADATA CORRECTION, AND ONE HONESTY CONSTRAINT
# -------------------------------------------------
# 1. SRR34748664 and SRR34748663 carry isolation_source "coastal Florida" in
#    metadata/SraRunTable.csv, and this figure used to draw them as a separate study -
#    open triangles, not joined to the Gulf series. That label is WRONG: per the user on
#    2026-11-24 they come from the same survey as the other nine and were mislabelled on
#    submission to the SRA. All 11 environmental samples are therefore one series, joined
#    by one line in date order. The correction is applied here rather than upstream, so
#    SraRunTable.csv and the merged tables keep the SRA's own record.
# 2. Genomes that were measured and found absent are drawn, faintly, rather than dropped.
#    A reader must be able to see that they were tested.
#
# Usage:
#   Rscript scripts/plot_gvmag_abundance.R <csv> <out-dir>

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(scales); library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: plot_gvmag_abundance.R <csv> <out-dir>")
in_csv <- args[1]; out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito, colour-blind safe, shared with every other figure in the project.
INC <- "#D55E00"; ENV <- "#0072B2"; INK <- "#1A1A1A"; MUTED <- "#5A5A5A"

d <- read_csv(in_csv, show_col_types = FALSE) %>%
  mutate(collection_date = as.Date(collection_date),
         detected = as.logical(detected),
         filter_set = factor(filter_set, levels = c("paper", "strict"),
                         labels = c("CoverM defaults (no identity filter)",
                                    "≥95% read identity, ≥75% aligned")),
         set = factor(genome_source_set, levels = c("incubation", "environmental"),
                      labels = c("assembled from an incubation sample",
                                 "assembled from an environmental sample")))

# The two metrics need DIFFERENT y scales - relative abundance spans orders of magnitude
# and wants a log axis, breadth is a fraction in [0, 1] and wants a linear one - and
# facet_grid cannot mix transforms across rows. So they are built as two separate plots
# and stacked with patchwork, which also lets them share one legend and one x axis.
#
# Zeros cannot be drawn on a log axis, so they are floored at half the smallest non-zero
# value and the caption says so.
floor_val <- min(d$relative_abundance_pct[d$relative_abundance_pct > 0], na.rm = TRUE) / 2

# The two runs the SRA mislabels as "coastal Florida" - see the header.
MISLABELLED <- c("SRR34748664", "SRR34748663")
SETCOL <- c("assembled from an incubation sample" = INC,
            "assembled from an environmental sample" = ENV)

base_theme <- theme_bw(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(colour = INK, size = 9.5, face = "bold"),
    legend.position = "top",
    legend.text = element_text(colour = INK, size = 9.5),
    plot.title = element_text(colour = INK, size = 10.5, face = "bold"),
    plot.caption = element_text(colour = MUTED, size = 8, hjust = 0),
    axis.title.y = element_text(colour = INK, size = 9.5),
    axis.text = element_text(colour = INK, size = 8.5),
    axis.text.x = element_text(lineheight = 0.9)
  )

# One panel row: `yvar` against time, faceted by filter setting. All 11 environmental
# samples are one series (see the metadata correction in the header), so one line runs
# through all of them. Genomes measured and found absent are drawn faintly rather than
# dropped, so a reader can see they were tested.
panel <- function(df, yvar, title, ylab) {
  ggplot(df, aes(collection_date, .data[[yvar]], group = genome, colour = set)) +
    geom_line(alpha = 0.45, linewidth = 0.4) +
    # Undetected: small and faint. Detected: solid and large. Size alone carries
    # detection, so nothing is encoded twice.
    geom_point(data = df %>% filter(!detected), alpha = 0.35, size = 0.9, shape = 16) +
    geom_point(data = df %>% filter(detected), size = 2.6, shape = 16) +
    facet_wrap(~ filter_set, nrow = 1) +
    scale_colour_manual(values = SETCOL, name = NULL, drop = FALSE) +
    # Without this the legend inherits the last layer's tiny translucent glyph, which is
    # unreadable; force one clear solid key per colour.
    guides(colour = guide_legend(
      override.aes = list(shape = 16, size = 3, alpha = 1, linewidth = 1))) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b\n%Y") +
    labs(title = title, x = NULL, y = ylab) +
    base_theme
}

d_ab <- d %>% mutate(value = pmax(relative_abundance_pct, floor_val))

p_ab <- panel(d_ab, "value",
              "Relative abundance: what share of each sample's reads landed on the genome",
              "% of reads (log scale)") +
  scale_y_log10(labels = label_number(drop0trailing = TRUE))

p_br <- panel(d, "covered_fraction",
              "Breadth: how much of the genome was covered at all - the test of real presence",
              "fraction of genome covered") +
  geom_hline(yintercept = 0.75, linetype = "22", colour = INK, linewidth = 0.45) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     labels = label_number(accuracy = 0.01)) +
  labs(caption = paste0(
    "Each line is one of the 26 GVMAGs across the 11 environmental samples, which are ",
    "one survey (station EH25).\nLarge solid markers pass the 75% breadth rule (dashed ",
    "line, MIUViG); small faint ones do not. ",
    "Abundance zeros floored at ", signif(floor_val, 2), "%."))

p <- (p_ab / p_br) + plot_layout(guides = "collect") &
  theme(legend.position = "top")

ggsave(file.path(out_dir, "gvmag_abundance_timeseries.pdf"), p,
       width = 11, height = 8.2, device = cairo_pdf, bg = "white")
ggsave(file.path(out_dir, "gvmag_abundance_timeseries.png"), p,
       width = 11, height = 8.2, dpi = 200, bg = "white")

cat(sprintf("wrote %s/gvmag_abundance_timeseries.{pdf,png}\n", out_dir))
cat(sprintf("%d GVMAGs x %d samples x %d filters\n",
            n_distinct(d$genome), n_distinct(d$sample), n_distinct(d$filter_set)))
det <- d %>% filter(detected) %>% count(filter_set, genome_source_set)
if (nrow(det)) print(as.data.frame(det)) else cat("no detections at the 75% breadth rule\n")
