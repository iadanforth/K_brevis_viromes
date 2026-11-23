#!/usr/bin/env Rscript
# Which vOTUs appeared in which samples, as a figure.
#
# Input : output/votu_sample_table.csv   (from scripts/votu_sample_table.py)
# Output: votu_sample_matrix.{pdf,png}
#
# The figure has to carry a distinction the data forces on it. "Appeared" means two
# different things:
#   ASSEMBLED - a GVMAG was binned out of that sample and clustered into this vOTU.
#               Known for all 20 samples.
#   DETECTED  - that sample's reads covered >=75% of the vOTU genome. Known for the 11
#               ENVIRONMENTAL samples only, because they are the only ones mapped.
# So a blank cell in an incubation column means NOT ASSESSED, not absent. Showing those
# in the same colour as a true absence would claim 108 absences that were never tested,
# which is why "not assessed" gets its own grey and sits apart in the legend.

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: plot_votu_sample_matrix.R <table.csv> <out-dir>")
in_csv <- args[1]; out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

d <- read_csv(in_csv, show_col_types = FALSE)

lev <- c("not assessed", "mapped, not detected", "ASSEMBLED", "ASSEMBLED+DETECTED")
d <- d %>% mutate(status_label = factor(status_label, levels = lev))

# samples left to right: the 9 incubation samples, then the 11 environmental by date
ord <- d %>% distinct(sample, sample_set, collection_date) %>%
  arrange(sample_set != "incubation", collection_date, sample)
d <- d %>% mutate(sample = factor(sample, levels = ord$sample),
                  votu = factor(votu, levels = rev(sort(unique(votu)))))
split_at <- sum(ord$sample_set == "incubation") + 0.5

pal <- c("not assessed"         = "#DCDCDC",
         "mapped, not detected" = "#F7F7F7",
         "ASSEMBLED"            = "#56B4E9",
         "ASSEMBLED+DETECTED"   = "#0072B2")

p <- ggplot(d, aes(sample, votu, fill = status_label)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_vline(xintercept = split_at, colour = "grey25", linewidth = 0.7) +
  annotate("text", x = split_at / 2 + 0.5, y = nlevels(d$votu) + 0.9,
           label = "INCUBATION  (not read-mapped)", size = 3.1,
           fontface = "bold", colour = "grey35") +
  annotate("text", x = (split_at + nlevels(d$sample)) / 2 + 0.5,
           y = nlevels(d$votu) + 0.9,
           label = "ENVIRONMENTAL  (read-mapped)", size = 3.1,
           fontface = "bold", colour = "grey35") +
  scale_fill_manual(values = pal, name = NULL, drop = FALSE) +
  coord_cartesian(ylim = c(0.5, nlevels(d$votu) + 1.4), clip = "off") +
  labs(
    title = "Which vOTUs appeared in which samples",
    subtitle = paste0(
      "ASSEMBLED = a GVMAG was binned from that sample and clustered into this vOTU (known for all 20 samples).\n",
      "DETECTED  = that sample's reads covered >=75% of the vOTU genome (known for the 11 environmental samples only).\n",
      "The 9 incubation samples were never read-mapped, so grey means NOT ASSESSED - it does not mean absent."),
    x = NULL, y = NULL,
    caption = "Data: votu_sample_table.csv  ·  wide form: votu_sample_matrix.csv"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid       = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7.6),
    axis.text.y      = element_text(size = 8.4),
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(colour = "grey30", size = 8.2, lineheight = 1.4),
    plot.caption     = element_text(colour = "grey45", size = 7.2, hjust = 0),
    plot.margin      = margin(10, 14, 8, 10)
  )

ggsave(file.path(out_dir, "votu_sample_matrix.pdf"), p, width = 10.5, height = 7, bg = "white")
ggsave(file.path(out_dir, "votu_sample_matrix.png"), p, width = 10.5, height = 7, dpi = 220, bg = "white")

cat("cells by status:\n")
d %>% count(status_label, .drop = FALSE) %>% as.data.frame() %>% print(row.names = FALSE)
cat("\nwrote figures to", out_dir, "\n")
