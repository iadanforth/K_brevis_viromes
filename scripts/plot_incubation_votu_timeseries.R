#!/usr/bin/env Rscript
# Abundance of the INCUBATION-derived vOTUs in the ENVIRONMENTAL metagenomes, through time.
#
# This is the figure that corresponds to Booker et al. 2026 Fig. 4: environmental
# metagenomes mapped onto the giant-virus genomes recovered from the bloom incubation
# experiments, plotted against collection date.
#
# WHAT IS PLOTTED, AND WHY EACH CHOICE
# ------------------------------------
#   series   The 12 vOTUs assembled from an incubation sample. The 7 vOTUs assembled from
#            environmental samples are deliberately NOT here: the question is whether the
#            viruses that bloomed in the incubations are present in the environment, and an
#            environmental vOTU answering "yes" about its own source sample answers nothing.
#   x        Collection date of the environmental sample.
#   y        Relative abundance - the percentage of that sample's reads that CoverM
#            assigned to the vOTU. This is the quantity Booker et al. report. Log scale,
#            because the values span five orders of magnitude (5e-06% to 0.36%).
#   mapping  CoverM's DEFAULT settings: no read-identity filter and no aligned-length
#            filter. The paper specifies no cutoffs, and a lenient threshold is what lets
#            an environmental read recruit to a genome assembled from a different
#            experiment. This is the `paper` arm of step 14.
#   per vOTU A vOTU can hold 2-3 near-identical GVMAGs, which split the reads between them.
#            The rollup sums each cluster's members back together, so a vOTU's abundance is
#            the abundance of the whole cluster, not of one arbitrary member.
#
# THE TWO SURVEYS
# ---------------
# Nine of the 11 environmental samples are the monthly Gulf of Mexico series
# (Dec 2020 - Aug 2021, station EH25 in the paper) - that is the bloom time course, and it
# is what the lines join. The other two are coastal Florida samples from a separate study
# (2021-09-29, 2022-03-08); they are drawn as open triangles at their real dates and are
# NOT joined to the series, because a line through them would imply a continuous 15-month
# survey that does not exist.
#
# COLOUR
# ------
# Twelve series need twelve distinguishable keys. The Okabe-Ito colour-blind-safe palette
# has eight, so colour is paired with line type: the first eight vOTUs are solid, the last
# four dashed. Both are mapped to `votu`, which merges them into one legend.
#
# Usage:
#   Rscript scripts/plot_incubation_votu_timeseries.R <votu_rollup.tsv> <out-dir> [filter]
#   filter: paper (default, lenient - CoverM defaults) | strict (>=95% read identity)
#
# Outputs (in <out-dir>):
#   incubation_votu_abundance.csv              self-contained data behind the figure
#   incubation_votu_timeseries.{pdf,png}       the figure

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("usage: plot_incubation_votu_timeseries.R <votu_rollup.tsv> <out-dir> [filter]")
in_tsv <- args[1]; out_dir <- args[2]
filt <- if (length(args) >= 3) args[3] else "paper"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

GULF <- "Gulf of Mexico (monthly)"
INK <- "#1A1A1A"; MUTED <- "#4A4A4A"

d <- read_tsv(in_tsv, show_col_types = FALSE) %>%
  filter(filter_set == filt, votu_source_set == "incubation") %>%
  mutate(collection_date = as.Date(collection_date),
         detected = as.logical(detected),
         abundance_pct = relative_abundance_pct_sum)

if (nrow(d) == 0)
  stop("no incubation rows with filter_set == '", filt, "' in ", in_tsv)

# Self-contained CSV: everything the figure draws, nothing else, so it can be rebuilt
# from this one small file with no other input.
out_csv <- file.path(out_dir, "incubation_votu_abundance.csv")
d %>%
  select(votu, sample, collection_date, survey, abundance_pct,
         reads = read_count_sum, breadth = best_covered_fraction,
         n_gvmags_in_votu = n_members, detected) %>%
  arrange(votu, collection_date) %>%
  write_csv(out_csv)

votus <- sort(unique(d$votu))
# Okabe-Ito, recycled across 12 series with line type carrying the second cycle.
okabe <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
           "#E69F00", "#56B4E9", "#8C6D00", "#000000")
pal <- setNames(okabe[(seq_along(votus) - 1) %% 8 + 1], votus)
lty <- setNames(ifelse(seq_along(votus) <= 8, "solid", "22"), votus)

gulf  <- d %>% filter(survey == GULF)
coast <- d %>% filter(survey != GULF)

p <- ggplot(d, aes(collection_date, abundance_pct,
                   colour = votu, linetype = votu, group = votu)) +
  geom_line(data = gulf, linewidth = 0.7, alpha = 0.9) +
  geom_point(data = gulf, size = 1.9) +
  geom_point(data = coast, shape = 2, size = 2.4, stroke = 0.9,
             show.legend = FALSE) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_linetype_manual(values = lty, name = NULL) +
  scale_y_log10(labels = label_number(drop0trailing = TRUE),
                breaks = 10^(-6:0)) +
  annotation_logticks(sides = "l", colour = "grey55", linewidth = 0.3,
                      short = unit(0.05, "cm"), mid = unit(0.1, "cm"),
                      long = unit(0.15, "cm")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
  labs(
    title = paste("Abundance of the incubation-derived vOTUs in the environmental",
                  "metagenomes"),
    subtitle = paste0(
      "The 12 giant-virus vOTUs recovered from the bloom incubation experiments, ",
      "quantified by mapping each environmental\nsample's reads onto them with CoverM at ",
      if (filt == "paper") "its default settings (no read-identity filter)."
      else "≥95% read identity.",
      " Lines join the 9-month Gulf of Mexico\nseries (station EH25); open triangles are ",
      "the two coastal Florida samples, which are a separate study."),
    x = NULL, y = "Relative abundance (% of sample reads)",
    caption = paste0(
      "Data: ", basename(out_csv), "  ·  step 14, `", filt, "` arm.  ",
      "Genome breadth stays below 20% for every point shown, so these are read-recruitment\n",
      "values, not confirmed genome detections - see README section 17d.")) +
  theme_bw(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    legend.position  = "right",
    legend.key.width = unit(1.3, "cm"),
    legend.text      = element_text(colour = INK, size = 9),
    plot.title       = element_text(colour = INK, size = 12, face = "bold"),
    plot.subtitle    = element_text(colour = MUTED, size = 8.8, lineheight = 1.3),
    plot.caption     = element_text(colour = MUTED, size = 7.6, hjust = 0,
                                    lineheight = 1.3),
    axis.title.y     = element_text(colour = INK, size = 9.5),
    axis.text        = element_text(colour = INK, size = 8.5),
    axis.text.x      = element_text(lineheight = 0.9))

ggsave(file.path(out_dir, "incubation_votu_timeseries.pdf"), p,
       width = 11, height = 6.2, device = cairo_pdf, bg = "white")
ggsave(file.path(out_dir, "incubation_votu_timeseries.png"), p,
       width = 11, height = 6.2, dpi = 200, bg = "white")

cat(sprintf("wrote %s\n", out_csv))
cat(sprintf("wrote %s/incubation_votu_timeseries.{pdf,png}\n", out_dir))
cat(sprintf("%d incubation vOTUs x %d environmental samples, filter = %s\n",
            n_distinct(d$votu), n_distinct(d$sample), filt))
cat(sprintf("abundance range: %.3g%% to %.3g%%\n",
            min(d$abundance_pct), max(d$abundance_pct)))
