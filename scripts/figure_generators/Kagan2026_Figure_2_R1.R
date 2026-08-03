# ============================================================================
# Comparative CNE figure: Actinopterygii-specific vs Gnathostome-conserved
#
# Gene assignment uses the `geneId` column directly. In the new table format
# `geneId` is already the correctly-assigned gene symbol for each CNE, so the
# previous flank-based reassignment (column Q / flank_geneIds + primary_flank())
# is no longer needed and has been removed.
#
# `geneId` is taken verbatim (only a defensive 'gene-' prefix strip + trim);
# LOC* gene-model IDs are excluded from the hotspot tables as before.
#
# Input files (expected in working directory):
#   actinopterygii_cne_final_table.tsv
#   gnathostomata_cne_final_table.tsv
#
# Output: cne_comparison_geneId_assignment.pdf
#
# Notes:
#   - Palette is Wong's colorblind-safe 8-colour set (Nature Methods 2011);
#     acti = #CC79A7 (reddish-purple), gnath = #009E73 (bluish-green).
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
  library(ggpubr)
  library(scales)
})

# ---- 0. Input file names ---------------------------------------------------
ACTI_FILE  <- "actinopterygii_cne_final_table_R1_20260731.tsv"
GNATH_FILE <- "gnathostomata_cne_final_table_R1_20260731.tsv"

# ---- 1. Colour palette (Wong 2011, colour-blind safe) ----------------------
ACTI    <- "#CC79A7"   # reddish-purple
GNATH   <- "#009E73"   # bluish-green
NEUTRAL <- "#999999"   # neutral grey for backgrounds

palette_pair <- c(`Actinopterygii-specific` = ACTI,
                  `Gnathostome-conserved`   = GNATH)

# ---- 2. Load and assign genes using the geneId column ----------------------
load_table <- function(path) {
  d <- read_tsv(path, show_col_types = FALSE,
                guess_max = 50000, na = c("", "NA"))
  d |>
    mutate(
      # geneId is already the correct assignment; just normalise lightly.
      gene      = trimws(sub("^gene-", "", as.character(geneId))),
      gene      = ifelse(gene == "" , NA_character_, gene),
      width     = as.integer(width),
      phastcons = suppressWarnings(as.numeric(phastcons)),
      active    = is_active == TRUE     | is_active == "TRUE",
      atac      = in_atac_peak == 1     | in_atac_peak == "1"
    )
}

acti  <- load_table(ACTI_FILE)  |> mutate(set = "Actinopterygii-specific")
gnath <- load_table(GNATH_FILE) |>
           # filter(width >= 25) |>
           mutate(set = "Gnathostome-conserved")

both <- bind_rows(acti, gnath) |>
          mutate(set = factor(set,
                              levels = c("Actinopterygii-specific",
                                         "Gnathostome-conserved")))

cat(sprintf("Acti CNEs:  %d  |  unique genes: %d  |  missing geneId: %d\n",
            nrow(acti),
            dplyr::n_distinct(acti$gene[!is.na(acti$gene)]),
            sum(is.na(acti$gene))))
cat(sprintf("Gnath CNEs: %d  |  unique genes: %d  |  missing geneId: %d\n",
            nrow(gnath),
            dplyr::n_distinct(gnath$gene[!is.na(gnath$gene)]),
            sum(is.na(gnath$gene))))

# ---- 3. Per-set summaries used inside the panels ---------------------------
n_acti  <- nrow(acti)
n_gnath <- nrow(gnath)

# Annotation categories (collapsed)
annot_cat <- function(a) {
  case_when(
    grepl("Intron",     a) ~ "Intron",
    grepl("Promoter",   a) ~ "Promoter",
    grepl("UTR",        a) ~ "UTR",
    grepl("Intergenic|Distal", a, ignore.case = TRUE) ~ "Intergenic",
    TRUE                  ~ "Other"
  )
}
both <- both |> mutate(annot_cat = annot_cat(annotation))

# CNEs per gene, using the geneId-based assignment
cnes_per_gene <- both |>
  filter(!is.na(gene)) |>
  count(set, gene, name = "n_cnes")

# Top 15 hotspots per set
top_hot <- cnes_per_gene |>
  filter(!grepl("^loc", gene, ignore.case = TRUE)) |>
  group_by(set) |>
  slice_max(n_cnes, n = 15, with_ties = FALSE) |>
  arrange(set, desc(n_cnes)) |>
  ungroup()

# ---- 3b. Statistical comparisons for panels A, B, C, E ---------------------
# All tests below treat each CNE as an independent observation, which is the
# standard approach for this kind of figure but is worth flagging: CNEs
# cluster heavily within genes (see panel D/F), so a locus with many CNEs
# contributes many correlated data points. Gene-level sensitivity checks
# (aggregating first, one value per gene) are printed to the console
# alongside the per-CNE tests actually used to annotate the plots.

fmt_p <- function(p) {
  if (is.na(p))   return("p = NA")
  if (p < 0.0001) return("p < 0.0001")
  sprintf("p = %.4f", p)
}
sig_stars <- function(p) {
  if (is.na(p))  return("n.s.")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  "n.s."
}

## --- Panel A: CNE width -----------------------------------------------------
# Continuous, heavily right-skewed (log-x axis) -> Wilcoxon rank-sum, not a t-test.
test_width <- wilcox.test(width ~ set, data = both)

gene_level_width <- both |>
  filter(!is.na(gene)) |>
  group_by(set, gene) |>
  summarise(med_width = median(width), .groups = "drop")
test_width_genelevel <- wilcox.test(med_width ~ set, data = gene_level_width)

cat(sprintf("\n[Panel A: CNE width]\n  Per-CNE Wilcoxon rank-sum: W = %.0f, %s\n",
            test_width$statistic, fmt_p(test_width$p.value)))
cat(sprintf("  Gene-level sensitivity check (one median width per gene): %s\n",
            fmt_p(test_width_genelevel$p.value)))

## --- Panel B: PhastCons score -----------------------------------------------
test_phast <- wilcox.test(phastcons ~ set, data = both)

gene_level_phast <- both |>
  filter(!is.na(gene)) |>
  group_by(set, gene) |>
  summarise(med_phast = median(phastcons, na.rm = TRUE), .groups = "drop")
test_phast_genelevel <- wilcox.test(med_phast ~ set, data = gene_level_phast)

cat(sprintf("\n[Panel B: PhastCons score]\n  Per-CNE Wilcoxon rank-sum: W = %.0f, %s\n",
            test_phast$statistic, fmt_p(test_phast$p.value)))
cat(sprintf("  Gene-level sensitivity check (one median score per gene): %s\n",
            fmt_p(test_phast_genelevel$p.value)))

## --- Panel C: annotation category vs set -----------------------------------
# Categorical x categorical -> chi-square test of independence (omnibus),
# plus per-category post-hoc two-proportion tests, FDR-corrected across the
# 5 categories since we're now doing multiple comparisons.
tab_annot  <- table(both$set, both$annot_cat)
test_annot <- chisq.test(tab_annot)

annot_levels <- levels(factor(both$annot_cat,
                              levels = c("Intron","Promoter","UTR",
                                         "Intergenic","Other")))
annot_post <- lapply(annot_levels, function(cat_lvl) {
  a_yes <- sum(both$set == "Actinopterygii-specific" & both$annot_cat == cat_lvl)
  a_tot <- sum(both$set == "Actinopterygii-specific")
  g_yes <- sum(both$set == "Gnathostome-conserved"   & both$annot_cat == cat_lvl)
  g_tot <- sum(both$set == "Gnathostome-conserved")
  pt <- prop.test(c(a_yes, g_yes), c(a_tot, g_tot))
  data.frame(annot_cat = cat_lvl, p_raw = pt$p.value)
}) |> bind_rows() |>
  mutate(p_adj = p.adjust(p_raw, method = "fdr"),
         stars = vapply(p_adj, sig_stars, character(1)))

cat(sprintf("\n[Panel C: annotation category]\n  Chi-square test of independence: X2 = %.1f, df = %d, %s\n",
            test_annot$statistic, as.integer(test_annot$parameter), fmt_p(test_annot$p.value)))
cat("  Per-category post-hoc two-proportion tests (FDR-corrected):\n")
for (i in seq_len(nrow(annot_post))) {
  cat(sprintf("    %-12s raw %s, FDR-adjusted %s (%s)\n",
              annot_post$annot_cat[i], fmt_p(annot_post$p_raw[i]),
              fmt_p(annot_post$p_adj[i]), annot_post$stars[i]))
}

## --- Panel E: activity / ATAC accessibility --------------------------------
# Binary outcome x set -> chi-square (2x2) for each metric separately.
test_active <- chisq.test(table(both$set, both$active))
test_atac   <- chisq.test(table(both$set, both$atac))

cat(sprintf("\n[Panel E: activity & accessibility]\n  gene active: X2 = %.1f, %s\n",
            test_active$statistic, fmt_p(test_active$p.value)))
cat(sprintf("  ATAC peak:   X2 = %.1f, %s\n",
            test_atac$statistic, fmt_p(test_atac$p.value)))

# ---- 4. Panel A — CNE width (log-x histogram) ------------------------------


medians_width <- both |> group_by(set) |>
  summarise(med = median(width, na.rm = TRUE), .groups = "drop")
ann_width <- sprintf("Median CNE width: Acti-spec = %d  |  Gnath-cons = %d  \u00b7  Wilcoxon %s",
                     medians_width$med[1], medians_width$med[2], fmt_p(test_width$p.value))

pA <- ggplot(both, aes(x = width, fill = set)) +
  geom_histogram(bins = 30, alpha = 0.6, color = "white",
                 linewidth = 0.2, position = "identity",
                 aes(y = after_stat(density))) +
  geom_vline(data = medians_width, aes(xintercept = med, color = set),
             linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  scale_x_log10(labels = label_comma()) +
  scale_fill_manual(values = palette_pair,
                    labels = c(sprintf("Acti-spec (n=%s)",  comma(n_acti)),
                               sprintf("Gnath-cons (n=%s)", comma(n_gnath)))) +
  scale_color_manual(values = palette_pair, guide = "none") +
  labs(title = "CNE length statistics",
       subtitle = ann_width,
       x = "CNE width (bp, log scale)", y = "density", fill = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey30"))

# ---- 5. Panel B — PhastCons score ------------------------------------------
ann_phast <- sprintf("Wilcoxon rank-sum test: %s", fmt_p(test_phast$p.value))

pB <- ggplot(both, aes(x = phastcons, fill = set)) +
  geom_histogram(bins = 40, alpha = 0.6, color = "white",
                 linewidth = 0.2, position = "identity",
                 aes(y = after_stat(density))) +
  scale_x_continuous(limits = c(0, 120)) +
  scale_fill_manual(values = palette_pair,
                    labels = c(sprintf("Acti-spec (n=%s)",  comma(n_acti)),
                               sprintf("Gnath-cons (n=%s)", comma(n_gnath)))) +
  labs(title = "Sequence conservation",
       subtitle = ann_phast,
       x = "PhastCons score", y = "density") +
  theme_minimal(base_size = 10) +
  theme(legend.position = c(0.75, 0.85),
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_blank(),
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey30"))

# ---- 6. Panel C — annotation category (grouped bar) ------------------------
annot_df <- both |>
  count(set, annot_cat) |>
  group_by(set) |>
  mutate(pct = n / sum(n) * 100) |>
  ungroup() |>
  mutate(annot_cat = factor(annot_cat,
                            levels = c("Intron","Promoter","UTR",
                                       "Intergenic","Other")))

star_df_C <- annot_df |>
  group_by(annot_cat) |>
  summarise(y = max(pct) + 4, .groups = "drop") |>
  left_join(annot_post, by = "annot_cat")

pC <- ggplot(annot_df, aes(x = annot_cat, y = pct, fill = set)) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.65, color = "white", linewidth = 0.2) +
  geom_text(data = star_df_C, aes(x = annot_cat, y = y, label = stars),
            inherit.aes = FALSE, size = 4, fontface = "bold", color = "grey20") +
  scale_fill_manual(values = palette_pair) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Genomic position",
       subtitle = sprintf("Chi-square test of independence: \u03c7\u00b2 = %.1f (df=%d), %s  \u00b7  stars = FDR-adjusted per-category",
                          test_annot$statistic, as.integer(test_annot$parameter), fmt_p(test_annot$p.value)),
       x = NULL, y = "% of CNEs", fill = NULL) +
  theme_minimal(base_size = 10) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        legend.position  = "none")

# ---- 7. Panel D — CNEs per gene (log-y) ------------------------------------
medians_cnes <- cnes_per_gene |>
  group_by(set) |>
  summarise(med = median(n_cnes), .groups = "drop")

# Pre-compute density per set per bin
binned <- cnes_per_gene |>
  group_by(set, n_cnes) |>
  summarise(n = n(), .groups = "drop_last") |>
  mutate(density = n / sum(n)) |>
  ungroup() |>
  filter(n_cnes <= 30, density > 0)

# Pick a sensible visual baseline well below the smallest density
y_floor <- min(binned$density) / 5

# Use geom_segment to draw bars from y_floor to density, dodged left/right
binned <- binned |>
  mutate(x_offset = ifelse(set == "Actinopterygii-specific", -0.2, 0.2))

pD <- ggplot(binned) +
  geom_segment(aes(x = n_cnes + x_offset, xend = n_cnes + x_offset,
                   y = y_floor, yend = density, color = set),
               linewidth = 2.5, alpha = 0.8, lineend = "butt") +
  geom_vline(data = medians_cnes, aes(xintercept = med, color = set),
             linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  #scale_y_log10(labels = label_percent(accuracy = 0.01)) +
  scale_color_manual(values = palette_pair,
                     name = NULL,
                     guide = guide_legend(override.aes = list(linewidth = 4))) +
  coord_cartesian(xlim = c(0.5, 30), ylim = c(y_floor, max(binned$density) * 1.3)) +
  labs(title = "CNEs per gene",
       subtitle = sprintf("Median CNEs per gene: Acti-spec = %d | Gnath-cons = %d",
                   medians_cnes$med[medians_cnes$set == "Actinopterygii-specific"],
                   medians_cnes$med[medians_cnes$set == "Gnathostome-conserved"]),
       x = "CNEs per gene",
       y = "fraction of genes (log scale)") +
  theme_minimal(base_size = 10) +
  theme(plot.title       = element_text(face = "bold"),
        plot.subtitle    = element_text(size = 9, color = "grey30"),
        legend.position  = "none")

# ---- 8. Panel E — activity & accessibility --------------------------------
act_df <- both |>
  group_by(set) |>
  summarise(
    `gene\nactive`   = mean(active, na.rm = TRUE) * 100,
    `ATAC\npeak`     = mean(atac,   na.rm = TRUE) * 100,
    .groups = "drop"
  ) |>
  pivot_longer(-set, names_to = "metric", values_to = "pct")

star_df_E <- act_df |>
  group_by(metric) |>
  summarise(y = max(pct) + 8, .groups = "drop") |>
  mutate(stars = c(sig_stars(test_active$p.value), sig_stars(test_atac$p.value))[
    match(metric, c("gene\nactive", "ATAC\npeak"))])

pE <- ggplot(act_df, aes(x = metric, y = pct, fill = set)) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.65, color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_dodge(width = 0.75),
            vjust = -0.4, size = 3, color = "grey20") +
  geom_text(data = star_df_E, aes(x = metric, y = y, label = stars),
            inherit.aes = FALSE, size = 5, fontface = "bold", color = "grey20") +
  scale_fill_manual(values = palette_pair,
                    labels = c(sprintf("Acti-spec (n=%s)",  comma(n_acti)),
                               sprintf("Gnath-cons (n=%s)", comma(n_gnath)))) +
  scale_y_continuous(limits = c(0, 112)) +
  labs(title = "Activity & accessibility",
       subtitle = sprintf("\u03c7\u00b2 test: active %s  \u00b7  ATAC peak %s",
                          fmt_p(test_active$p.value), fmt_p(test_atac$p.value)),
       x = NULL, y = "% of CNEs") +
  theme_minimal(base_size = 10) +
  theme(legend.position  = "none",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey30"))

# ---- 9. Panel F — top-15 hotspot tables ------------------------------------
hot_a <- top_hot |> filter(set == "Actinopterygii-specific") |>
  mutate(
    rank      = row_number(),
    n_text    = sprintf("%2d", n_cnes),
    gene_text = gene,
  )
hot_g <- top_hot |> filter(set == "Gnathostome-conserved") |>
   mutate(
    rank      = row_number(),
    n_text    = sprintf("%2d", n_cnes),
    gene_text = gene,
  )
pF <- ggplot() +
  annotate("text", x = 0.05,  y = 16.5, label = "Acti-spec",  fontface = "bold",
           color = ACTI,  hjust = 0, size = 4) +
  annotate("text", x = 0.55,  y = 16.5, label = "Gnath-cons", fontface = "bold",
           color = GNATH, hjust = 0, size = 4) +
    # Column headers under each set
  annotate("text", x = 0,    y = 15, label = "CNE no.", fontface = "bold",
           color = "grey30", hjust = 0, size = 3) +
  annotate("text", x = 0.15, y = 15, label = "Gene name", fontface = "bold",
           color = "grey30", hjust = 0, size = 3) +
  annotate("text", x = 0.5,  y = 15, label = "CNE no.", fontface = "bold",
           color = "grey30", hjust = 0, size = 3) +
  annotate("text", x = 0.65, y = 15, label = "Gene name", fontface = "bold",
           color = "grey30", hjust = 0, size = 3) +
  # Acti panel: count in mono, gene name in italic
  geom_text(data = hot_a, aes(x = 0.05, y = 15 - rank, label = n_text),
            hjust = 0, size = 3.2, family = "mono") +
  geom_text(data = hot_a, aes(x = 0.15, y = 15 - rank, label = gene_text),
            hjust = 0, size = 3.2, fontface = "italic") +
  # Gnath panel: same pattern
  geom_text(data = hot_g, aes(x = 0.55, y = 15 - rank, label = n_text),
            hjust = 0, size = 3.2, family = "mono") +
  geom_text(data = hot_g, aes(x = 0.65, y = 15 - rank, label = gene_text),
            hjust = 0, size = 3.2, fontface = "italic") +
  scale_x_continuous(limits = c(-0.05, 1.0)) +
  scale_y_continuous(limits = c(-1, 17)) +
  labs(title = "Top 15 CNE hotspot genes ") +
  theme_void(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0))

# ---- 10. Compose ----------------------------------------------------------
fig <- ggarrange(
  pA, pB, pC, pD, pE, pF,
  ncol = 2, nrow = 3,
  labels = c("A", "B", "C", "D", "E", "F"),
  font.label = list(size = 18, face = "bold", color = "black"),
  label.x = 0.01,    # offset from left edge of each panel
  label.y = 1.00,    # offset from top
  align = "hv"       # align horizontal+vertical axes across panels
)

# Add the overall title separately
#fig <- annotate_figure(
#  fig,
#  top = text_grob(
#    sprintf(
#      "Comparison of actinopterygii-specific (n=%s) and gnathostome-conserved (n=%s) CNEs",
#      comma(n_acti), comma(n_gnath)),
#    face = "bold", size = 13
#  )
#)

ggsave("cne_comparison_geneId_assignment_R1_20260731.pdf", fig,
       width = 8, height = 12, dpi = 160, bg = "white")
