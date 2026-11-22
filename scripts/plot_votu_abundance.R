#!/usr/bin/env Rscript
# vOTU abundance through time across the environmental samples.
#
# Input : output/abundance/votu_abundance_environmental.csv  (written by
#         merge_votu_abundance.py --out-csv). That file is self-contained: this script
#         needs nothing else, so the figure can be rebuilt from one small CSV.
# Output: votu_abundance_timeseries.{pdf,png}
#
# Two things about the data drive the design, and both are honesty issues:
#
# 1. The 11 environmental samples are NOT one survey. Nine are a monthly Gulf of Mexico
#    series (Dec 2020 - Aug 2021, station EH25 in Booker et al.); two are from a separate
#    coastal Florida study (2021-09-29, 2022-03-08). A single line through all eleven
#    would imply a continuous 15-month survey that does not exist, so the line is drawn
#    through the Gulf series only.
# 2. Each vOTU is detected in exactly ONE sample. That is the result, not a plotting
#    artifact: the giant virus recovered from a given month is not detectable in any
#    other month. Undetected vOTUs are drawn as faint grey lines so the reader can see
#    that they were measured and were absent, rather than simply omitted.
#
# Usage:
#   Rscript scripts/plot_votu_abundance.R <csv> <out-dir> [metric]
#   metric: relative_abundance_pct (default) | rpkm | mean_coverage

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: plot_votu_abundance.R <csv> <out-dir> [metric]")
in_csv  <- args[1]; out_dir <- args[2]
metric  <- if (length(args) >= 3) args[3] else "relative_abundance_pct"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

d <- read_csv(in_csv, show_col_types = FALSE) %>%
  mutate(collection_date = as.Date(collection_date),
         detected = as.logical(detected),
         value = .data[[metric]])

gulf <- d %>% filter(survey == "Gulf of Mexico (monthly)")
fl   <- d %>% filter(survey != "Gulf of Mexico (monthly)")
divider <- mean(c(max(gulf$collection_date), min(fl$collection_date)))

# vOTUs detected at least once get a colour and a label; the rest stay grey.
seen <- d %>% group_by(votu) %>% summarise(any = any(detected), .groups = "drop") %>%
  filter(any) %>% pull(votu) %>% sort()
# Okabe-Ito, colour-blind safe, matching the project's other figures.
pal <- setNames(c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
                  "#56B4E9", "#D55E00", "#000000")[seq_along(seen)], seen)

d    <- d    %>% mutate(grp = ifelse(votu %in% seen, votu, "other"))
gulf <- gulf %>% mutate(grp = ifelse(votu %in% seen, votu, "other"))
lab  <- d %>% filter(detected)

ylab <- c(relative_abundance_pct = "Relative abundance (% of reads)",
          rpkm = "Abundance (RPKM)",
          mean_coverage = "Mean coverage depth (x)")[metric]
if (is.na(ylab)) ylab <- metric

p <- ggplot() +
  # everything never detected, behind
  geom_line(data = filter(gulf, grp == "other"),
            aes(collection_date, value, group = votu),
            colour = "grey85", linewidth = 0.45) +
  geom_vline(xintercept = as.numeric(divider), linetype = "dotted",
             colour = "grey55", linewidth = 0.45) +
  # the seven that are detected somewhere
  geom_line(data = filter(gulf, grp != "other"),
            aes(collection_date, value, colour = grp, group = votu), linewidth = 0.85) +
  # Two explicit point layers rather than one with a conditional fill: hollow for
  # "reads recruited but not present", solid for a real detection.
  geom_point(data = filter(d, grp != "other", !detected),
             aes(collection_date, value, colour = grp, shape = survey),
             fill = "white", size = 2.5, stroke = 0.9) +
  geom_point(data = filter(d, grp != "other", detected),
             aes(collection_date, value, colour = grp, fill = grp, shape = survey),
             size = 3.4, stroke = 0.9) +
  # Plain labels, nudged right. Every detection sits at a distinct height, so they do
  # not collide and no repel package is needed.
  geom_text(data = lab, aes(collection_date, value, label = votu, colour = grp),
            size = 3.2, fontface = "bold", vjust = -1.1, show.legend = FALSE) +
  scale_colour_manual(values = pal, name = "vOTU detected somewhere", na.value = "grey85") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_shape_manual(values = c(`Gulf of Mexico (monthly)` = 21,
                                `coastal Florida (separate study)` = 24),
                     name = "Sample origin") +
  scale_x_date(date_breaks = "2 months", date_labels = "%b\n%Y",
               expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(
    title = "Giant-virus vOTU abundance in the environmental samples",
    subtitle = paste0(
      "Environmental reads mapped to the 19 vOTU representative genomes (CoverM, >=95% read identity).\n",
      "A FILLED point is a detection: >=75% of the genome covered. A HOLLOW point means reads recruited\n",
      "but the genome is not present. Grey lines are the 12 vOTUs never detected in any sample.\n",
      "The dotted line separates the monthly Gulf of Mexico series (Dec 2020 - Aug 2021) from two samples\n",
      "of a separate coastal Florida study (Sep 2021, Mar 2022)."),
    x = NULL, y = ylab,
    caption = paste0(
      "Each detected vOTU appears in exactly one sample: the giant virus recovered from a given month is not\n",
      "detectable in any other month. Data: votu_abundance_environmental.csv")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92"),
    legend.position  = "bottom",
    legend.box       = "vertical",
    plot.title       = element_text(face = "bold", colour = "grey10"),
    plot.subtitle    = element_text(colour = "grey30", size = 8.6, lineheight = 1.35),
    plot.caption     = element_text(colour = "grey45", size = 7.6, hjust = 0, lineheight = 1.3)
  )

ggsave(file.path(out_dir, "votu_abundance_timeseries.pdf"), p,
       width = 10.5, height = 7.2, bg = "white")
ggsave(file.path(out_dir, "votu_abundance_timeseries.png"), p,
       width = 10.5, height = 7.2, dpi = 220, bg = "white")

cat("metric:", metric, "\n")
cat("vOTUs detected at least once:", length(seen), "of", n_distinct(d$votu), "\n  ",
    paste(seen, collapse = ", "), "\n")
d %>% filter(detected) %>%
  select(votu, sample, collection_date, all_of(metric), covered_fraction) %>%
  arrange(collection_date) %>% as.data.frame() %>% print(row.names = FALSE)
cat("\nwrote figures to", out_dir, "\n")
