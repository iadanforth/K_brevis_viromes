#!/usr/bin/env Rscript
# vOTU abundance against Karenia brevis cell concentration, with a linear model.
#
# Inputs
#   1. votu_rollup.tsv          per-vOTU abundance per environmental sample, written by
#                               merge_gvmag_abundance.py. Updated 2026-11-23 (prompt 60):
#                               step 14 now maps against each individual GVMAG, so the
#                               per-vOTU numbers come from the rollup, which sums each
#                               cluster's members back together. All 19 vOTUs are fitted,
#                               not just the 7 assembled from environmental samples - the
#                               question the paper asks is whether the INCUBATION viruses
#                               track cell counts, so excluding them would beg it.
#   2. karenia_cell_counts.csv  sample, collection_date, karenia_cells_per_L
#
# The cell counts are NOT produced by this pipeline. Booker et al. 2026 plot them in
# Fig. 4A but give no numbers in the text, no units in the caption and no source; the
# only numeric counts in the paper are from lab cultures and one New Pass sample, neither
# of which is station EH25. Fill karenia_cell_counts.csv from a real source - their
# Supplemental Data 14, or the FWC/FWRI HAB Monitoring Database, which publishes counts
# for Florida coastal stations - and record which in the `source` column. Values read by
# eye off a published figure are not data and should not be fitted.
#
# Outputs
#   votu_vs_karenia_cells.{pdf,png}   abundance vs cell count, one panel per vOTU, with fit
#   votu_vs_karenia_lm.csv            slope, intercept, R2, p, n for every vOTU
#
# Usage:
#   Rscript scripts/plot_votu_vs_cells.R <votu_rollup.tsv> <cells.csv> <out-dir> \
#       [filter] [metric]
#   filter: strict (default, >=95% read identity) | paper (CoverM defaults)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
  stop("usage: plot_votu_vs_cells.R <votu_rollup.tsv> <cells.csv> <out-dir> [filter] [metric]")
ab_csv <- args[1]; cell_csv <- args[2]; out_dir <- args[3]
filt   <- if (length(args) >= 4) args[4] else "strict"
metric <- if (length(args) >= 5) args[5] else "relative_abundance_pct_sum"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# votu_rollup.tsv is tab-separated and holds both filter settings stacked; keep one.
reader <- if (grepl("\\.tsv$", ab_csv)) read_tsv else read_csv
ab <- reader(ab_csv, show_col_types = FALSE) %>%
  filter(filter_set == filt) %>%
  mutate(collection_date = as.Date(collection_date), detected = as.logical(detected))
if (nrow(ab) == 0) stop("no rows with filter_set == '", filt, "' in ", ab_csv)
cells <- read_csv(cell_csv, show_col_types = FALSE) %>%
  mutate(collection_date = as.Date(collection_date),
         karenia_cells_per_L = suppressWarnings(as.numeric(karenia_cells_per_L)))

n_have <- sum(!is.na(cells$karenia_cells_per_L))
if (n_have < 3) {
  cat("\n", strrep("-", 72), "\n", sep = "")
  cat("STOPPING: karenia_cell_counts.csv has ", n_have,
      " cell counts filled in; a linear model needs at least 3.\n", sep = "")
  cat("This is a data gap, not a code problem. The counts are not in the paper's text.\n")
  cat("Fill the karenia_cells_per_L column from a real source (Supplemental Data 14, or\n")
  cat("the FWC/FWRI HAB Monitoring Database), note it in `source`, and re-run.\n")
  cat(strrep("-", 72), "\n\n", sep = "")
  quit(status = 0)
}

# Keep only the columns the join needs. Both files carry a `survey` column, and joining
# on them both silently produces survey.x / survey.y, after which `aes(shape = survey)`
# fails with "object 'survey' not found".
d <- ab %>%
  inner_join(cells %>% select(sample, collection_date, karenia_cells_per_L),
             by = c("sample", "collection_date")) %>%
  filter(!is.na(karenia_cells_per_L)) %>%
  mutate(value = .data[[metric]],
         # The strip tells the reader where each vOTU came from, because that is the whole
         # comparison: an incubation vOTU correlating with cell counts in the environment
         # would be the paper's result; an environmental one doing so is near-tautological,
         # since it was assembled from one of these very samples.
         votu_label = paste0(votu, "  \u00b7  ", substr(votu_source_set, 1, 3)))

# One fit per vOTU. lm() needs at least 3 points and some spread in x; report rather than
# silently emit a fit from two points.
fits <- d %>% group_by(votu) %>% group_modify(~{
  if (nrow(.x) < 3 || length(unique(.x$karenia_cells_per_L)) < 3) {
    return(tibble(slope = NA_real_, intercept = NA_real_, r_squared = NA_real_,
                  p_value = NA_real_, n = nrow(.x), note = "too few points to fit"))
  }
  m <- lm(value ~ karenia_cells_per_L, data = .x); sm <- summary(m)
  tibble(slope = coef(m)[2], intercept = coef(m)[1],
         r_squared = sm$r.squared,
         p_value = if (nrow(sm$coefficients) > 1) sm$coefficients[2, 4] else NA_real_,
         n = nrow(.x), note = "")
}) %>% ungroup()

write_csv(fits, file.path(out_dir, "votu_vs_karenia_lm.csv"))

lab <- fits %>% mutate(
  txt = ifelse(is.na(r_squared), note,
               sprintf("R² = %.2f, p = %.3g, n = %d", r_squared, p_value, n)))

p <- ggplot(d, aes(karenia_cells_per_L, value)) +
  geom_smooth(method = "lm", se = TRUE, colour = "#0072B2",
              fill = "#0072B2", alpha = 0.13, linewidth = 0.8, formula = y ~ x) +
  # All 11 environmental samples are one survey: the SRA's "coastal Florida" label on
  # SRR34748664 and SRR34748663 is a mislabelling (user, 2026-11-24), so there is no
  # sample-origin distinction left for shape to carry.
  geom_point(aes(fill = detected), shape = 21, colour = "#1A1A1A",
             size = 2.8, stroke = 0.8) +
  geom_text(data = lab, aes(label = txt), x = Inf, y = Inf, hjust = 1.04, vjust = 1.5,
            size = 3.0, colour = "grey30", inherit.aes = FALSE) +
  facet_wrap(~ votu_label, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "white"),
                    name = "Detected (>=75% genome breadth)") +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  labs(
    title = "vOTU abundance against Karenia brevis cell concentration",
    subtitle = paste0(
      "All 19 vOTUs (inc = assembled from an incubation sample, env = from an environmental one).\n",
      "Linear fit per vOTU; shaded band is the 95% CI.\n",
      "Hollow points are samples where reads recruited but the genome was not present (<75% breadth) -\n",
      "they are genuine zeros for abundance, not missing data."),
    x = expression(italic("Karenia brevis")~"(cells L"^-1*")"),
    y = c(relative_abundance_pct = "Relative abundance (% of reads)",
          rpkm = "Abundance (RPKM)", mean_coverage = "Mean coverage (x)")[metric],
    caption = "Cell counts are external to this pipeline; their provenance is in the `source` column of karenia_cell_counts.csv."
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey30", size = 8.4, lineheight = 1.35),
        plot.caption = element_text(colour = "grey45", size = 7.4, hjust = 0))

ggsave(file.path(out_dir, "votu_vs_karenia_cells.pdf"), p, width = 13, height = 8.5, bg = "white")
ggsave(file.path(out_dir, "votu_vs_karenia_cells.png"), p, width = 13, height = 8.5, dpi = 200, bg = "white")

cat("cell counts available for", n_have, "of", nrow(cells), "samples\n\n")
as.data.frame(fits) %>% print(row.names = FALSE)
cat("\nwrote", file.path(out_dir, "votu_vs_karenia_lm.csv"), "and the figures\n")
