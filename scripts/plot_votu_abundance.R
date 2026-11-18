#!/usr/bin/env Rscript
# Figure: K. brevis vOTU abundance in the environment through time.
#
# Input : votu_abundance_long.tsv from merge_votu_abundance.py
# Output: votu_abundance_timeseries.{pdf,png}  - abundance vs time, one panel per vOTU
#         votu_abundance_heatmap.{pdf,png}     - all vOTUs x time as a heatmap
#
# Two things about the data drive the design, and both are honesty issues rather than
# aesthetic ones:
#
# 1. The 11 "environmental" samples are NOT one survey. Nine are a monthly Gulf of Mexico
#    series (2020-12 to 2021-08, month precision). Two are from a different study in
#    coastal Florida waters (2021-09-29 and 2022-03-08, day precision, with depths). Drawing
#    a single connecting line through all eleven would imply a continuous time series across
#    two locations and two studies. So the line is drawn through the Gulf series only, and
#    the Florida samples are plotted as distinct points beyond a divider.
#
# 2. A zero is not a small number, it is a non-detection. Points failing the >=75% genome
#    breadth rule are drawn hollow at the axis floor rather than joined into the line, so a
#    vOTU that is genuinely absent for months does not look like one at low abundance.
#
# Usage:
#   Rscript scripts/plot_votu_abundance.R <long.tsv> <out-dir> [metric]
#   metric: rpkm_detected (default) | mean_coverage | relative_abundance_pct | tpm

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(scales); library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: plot_votu_abundance.R <long.tsv> <out-dir> [metric]")
in_tsv  <- args[1]
out_dir <- args[2]
metric  <- if (length(args) >= 3) args[3] else "rpkm_detected"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

d <- read_tsv(in_tsv, show_col_types = FALSE) %>%
  mutate(
    collection_date = as.Date(collection_date),
    survey = ifelse(grepl("Gulf", isolation_source), "Gulf of Mexico (monthly)",
                    "coastal Florida (separate study)"),
    value  = .data[[metric]],
    detected = present
  )

# Gulf series only for the connecting line; the two Florida samples stay as points.
gulf <- d %>% filter(survey == "Gulf of Mexico (monthly)")
fl   <- d %>% filter(survey != "Gulf of Mexico (monthly)")
divider <- mean(c(max(gulf$collection_date), min(fl$collection_date)))

# Colour-blind-safe (Okabe-Ito): blue = discovered in the incubation experiment,
# orange = discovered in an environmental sample (a positive control for the mapping).
pal <- c(incubation = "#0072B2", environmental = "#E69F00")

metric_label <- c(
  rpkm_detected          = "Abundance (RPKM, non-detections set to 0)",
  rpkm                   = "Abundance (RPKM, unfiltered)",
  mean_coverage          = "Mean coverage depth (x)",
  relative_abundance_pct = "Relative abundance (% of reads)",
  tpm                    = "Abundance (TPM)"
)[metric]
if (is.na(metric_label)) metric_label <- metric

# ---------------- Figure 1: one small panel per vOTU ----------------
# EVERY vOTU gets a panel, including those never detected. An earlier version dropped
# never-detected vOTUs to avoid empty panels, which was exactly wrong for the question
# being asked: "how abundant were the incubation vOTUs in the environment?" is answered by
# the twelve flat panels, and hiding them would hide the result.
#
# Raw `rpkm` is plotted rather than `rpkm_detected`, so read recruitment that FAILS the
# breadth test is still visible, with detection shown by whether the point is filled. That
# distinction matters here: vOTU13 reaches 32x mean coverage in SRR34748664 over 0.9% of its
# genome - reads piling onto one small conserved region. Plotting rpkm_detected would show
# a flat zero (true, but mute); plotting raw rpkm with a hollow point shows the recruitment
# AND that it is not a real detection.
keep   <- sort(unique(d$votu))
d_k    <- d    %>% mutate(value = rpkm)
gulf_k <- gulf %>% mutate(value = rpkm)

# Detection is encoded by FILL: a point filled with the series colour is a real
# detection, a white-filled point is not. scale_fill_identity is used because the fill
# has to carry the series colour, which a two-value manual scale cannot do. (An earlier
# version mapped TRUE->NA / FALSE->"white", which rendered detections as transparent
# rather than solid - the opposite of the intent.)
d_k <- d_k %>% mutate(fill_col = ifelse(detected, pal[votu_source_set], "white"))

# Panels for vOTUs that are never detected are labelled as such. Without this the free
# y-axis is actively misleading: vOTU13's panel climbs to 275 RPKM and reads as a rising
# abundance, when every one of those points fails the breadth test and the vOTU is absent.
never_lab <- d_k %>% group_by(votu, votu_source_set) %>%
  summarise(any_det = any(detected), .groups = "drop") %>%
  filter(!any_det) %>%
  mutate(panel = paste0(votu, "  (", substr(votu_source_set, 1, 5), ")"),
         label = "never detected")

d_k      <- d_k      %>% mutate(panel = paste0(votu, "  (", substr(votu_source_set, 1, 5), ")"))
gulf_k   <- gulf_k   %>% mutate(panel = paste0(votu, "  (", substr(votu_source_set, 1, 5), ")"),
                                fill_col = ifelse(detected, pal[votu_source_set], "white"))

p1 <- ggplot(d_k, aes(collection_date, value)) +
  geom_vline(xintercept = as.numeric(divider), linetype = "dotted",
             colour = "grey55", linewidth = 0.4) +
  geom_line(data = gulf_k, aes(colour = votu_source_set, group = votu), linewidth = 0.6) +
  geom_point(aes(colour = votu_source_set, shape = survey, fill = fill_col),
             size = 2.3, stroke = 0.8) +
  geom_text(data = never_lab, aes(label = label), x = Inf, y = Inf,
            hjust = 1.05, vjust = 1.4, size = 2.7, colour = "grey45",
            inherit.aes = FALSE) +
  facet_wrap(~ panel, ncol = 4, scales = "free_y") +
  scale_colour_manual(values = pal, name = "vOTU discovered in") +
  scale_fill_identity() +
  scale_shape_manual(values = c(`Gulf of Mexico (monthly)` = 21,
                                `coastal Florida (separate study)` = 24),
                     name = "Sample origin") +
  scale_x_date(date_breaks = "3 months", date_labels = "%b\n%Y") +
  labs(
    title = "Karenia brevis giant-virus vOTU abundance in the environment through time",
    subtitle = paste0(
      "Environmental reads mapped to vOTU representative genomes (CoverM, >=95% read identity). ",
      "FILLED point = detected\n(>=75% of the genome covered). HOLLOW point = reads recruited but ",
      "the virus is NOT present. Note the free y-axes:\na tall hollow curve is read recruitment ",
      "onto a small conserved region, not abundance. Dotted line separates the monthly\n",
      "Gulf of Mexico series from two samples of a separate coastal Florida study."),
    x = NULL, y = "Read recruitment (RPKM)",
    caption = paste0(
      "Panel label shows where the vOTU was discovered (incub = incubation experiment, ",
      "envir = environmental sample).\nAll 12 incubation-derived vOTUs are undetected in ",
      "every environmental sample.")
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92"),
    strip.text = element_text(face = "bold", colour = "grey20"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", colour = "grey10"),
    plot.subtitle = element_text(colour = "grey30", size = 8, lineheight = 1.2),
    plot.caption = element_text(colour = "grey45", size = 7, hjust = 0)
  )

ggsave(file.path(out_dir, "votu_abundance_timeseries.pdf"), p1, width = 13, height = 9, bg = "white")
ggsave(file.path(out_dir, "votu_abundance_timeseries.png"), p1, width = 13, height = 9, dpi = 200, bg = "white")

# ---------------- Figure 2: heatmap of every vOTU ----------------
# log1p because abundance spans orders of magnitude; grey = not detected, so a true
# absence is visually distinct from a low value rather than merging into the dark end.
# Breadth, not abundance, is the right fill for the heatmap: it is what decides presence,
# and in these data it is sharply bimodal (everything is either <16% or >99%), which a
# depth-based fill would blur.
hm <- d %>%
  mutate(
    plot_val = covered_fraction * 100,
    votu_lab = paste0(votu, ifelse(votu_source_set == "incubation", "", " *"))
  )

p2 <- ggplot(hm, aes(factor(collection_date), reorder(votu_lab, desc(votu_lab)))) +
  geom_tile(aes(fill = plot_val), colour = "white", linewidth = 0.4) +
  geom_text(data = subset(hm, detected), label = "*", colour = "white",
            size = 4, vjust = 0.72) +
  scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, 100),
                       na.value = "grey92", name = "Genome breadth (%)") +
  labs(
    title = "Genome breadth: which vOTUs are actually present, and when",
    subtitle = paste0(
      "Breadth is bimodal - every vOTU-sample pair is either <16% or >99% covered, so the ",
      ">=75% presence rule is not a\njudgement call. White * marks a detection. ",
      "* after a vOTU name = it was assembled from an environmental sample,\nso its own ",
      "detection is a positive control for the mapping."),
    x = "Collection date", y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey30", size = 8, lineheight = 1.2)
  )

ggsave(file.path(out_dir, "votu_abundance_heatmap.pdf"), p2, width = 11, height = 6.5, bg = "white")
ggsave(file.path(out_dir, "votu_abundance_heatmap.png"), p2, width = 11, height = 6.5, dpi = 200, bg = "white")

# ---------------- console summary ----------------
# Figure 1 plots RAW rpkm, so the summary reports raw rpkm too - reporting the
# breadth-gated value here would print 0 for every undetected vOTU and hide the fact that
# reads recruited at all, which is the interesting part.
cat("\nFigure 1 y-axis: raw RPKM (all", n_distinct(d$votu), "vOTUs shown, detected or not)\n")
cat("Figure 2 fill:   genome breadth (%)\n")
det <- d %>% group_by(votu, votu_source_set) %>%
  summarise(detected_in = sum(detected), n = n(),
            max_rpkm = max(rpkm), max_breadth = max(covered_fraction),
            max_mean_cov = max(mean_coverage), .groups = "drop") %>%
  arrange(desc(detected_in), desc(max_rpkm))
cat("\nper vOTU (detected_in = samples passing the >=75% breadth rule):\n")
as.data.frame(det) %>% print(row.names = FALSE)
never <- det %>% filter(detected_in == 0) %>% pull(votu)
cat("\nnever detected in any environmental sample:", length(never), "vOTUs\n  ",
    paste(never, collapse = ", "), "\n")
cat("\nNOTE: a high max_rpkm with a low max_breadth is read recruitment onto a small part\n",
    "of the genome, not presence - vOTU13 reaches 32x mean coverage over 0.9% of its\n",
    "genome in SRR34748664. That is why presence is gated on breadth.\n")
cat("\nwrote figures to", out_dir, "\n")
