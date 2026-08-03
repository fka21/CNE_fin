# ============================================================================
# Differential TF usage in lineage-specific CNEs — combined three-panel figure
#
#   Panel A : Volcano of TF motif enrichment in Actinopterygii-specific CNEs
#   Panel B : Volcano of TF motif enrichment in Gnathostomata-specific CNEs
#   Panel C : MA-style differential plot (Actinopterygii vs Gnathostomata)
#
# Geometric relationship between the panels:
#   A and B each plot log2FE vs -log10 FDR for one clade's CNE set, colored by
#   TF mechanism category. C collapses the two log2FE axes onto a single
#   (mean, Δ) coordinate system — a 45° rotation that puts the clade-bias
#   axis on a dedicated vertical dimension. Panels are colored on different
#   schemes by design: A/B by TF mechanism (biological prior), C by
#   differential class (data-driven, post-hoc).
#
# Inputs:
#   actinopterygii_vs_random_genomic.tf_enrichment_real_vs_random_genomic.tsv
#   gnathostomata_vs_random_genomic.tf_enrichment_real_vs_random_genomic.tsv
#
# Outputs:
#   TF_enrichment_combined_figure.pdf
#   TF_enrichment_combined_figure.png
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(scales)
})

# ----------------------------------------------------------------------------
# 0. CONFIG
# ----------------------------------------------------------------------------

FILE_ACTIN <- "actinopteriigy_vs_random_genomic.tf_enrichment_real_vs_random_genomic_R1_20260731.tsv"
FILE_GNATH <- "gnathostomata_vs_random_genomic.tf_enrichment_real_vs_random_genomic_R1_20260731.tsv"

FDR_CUTOFF      <- 0.05
FOLD_CUTOFF     <- 1.5
LOG2FC_CUTOFF   <- log2(FOLD_CUTOFF)

## Volcano panel A/B auto-label selection (both panels use the same
## top-by-FDR / top-by-effect selection; see build_volcano()).
N_GNATH_LABEL_BY_FDR    <- 12
N_GNATH_LABEL_BY_LOG2FE <- 8

OUTPUT_STEM <- "SupplFig_TF_motif_enrichment"

# ----------------------------------------------------------------------------
# 1. HELPERS, PALETTES
# ----------------------------------------------------------------------------

## Clade anchor colors (consistent with other Acti-spec/Gnath-cons figures).
ACTI  <- "#CC79A7"   # reddish-purple
GNATH <- "#009E73"   # bluish-green

## Panel A/B: TF mechanism categories.
palette_mech <- c(
  "Chromatin/Epigenetic"     = "#7F77DD",
  "Cell cycle/Proliferation" = "#1D9E75",
  "Stress/Homeostasis"       = "#D85A30",
  "Nuclear receptor"         = "#BA7517",
  "Hematopoietic"            = "#D4537E",
  "Myogenic"                 = "#378ADD",
  "Neural crest/EMT"         = "#A05A9C",
  "Signaling effector"       = "#6B8E23",
  "ETS family"               = "#5DA9C4",
  "Tissue patterning"        = "#C8A030",
  "General/KLF/SP"           = "#888780",
  "Other"                    = "#BBBBBB"
)

## Shared theme for both panels — tweak once, apply everywhere.
panel_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.background  = element_rect(fill = "#FAFAF8", color = NA),
      plot.background   = element_rect(fill = "white",   color = NA),
      panel.grid.major  = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor  = element_blank(),
      axis.line         = element_line(color = "#666666", linewidth = 0.5),
      axis.ticks        = element_line(color = "#666666", linewidth = 0.4),
      axis.ticks.length = unit(3, "pt"),
      plot.title        = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle     = element_text(size = 9.5, color = "#555555", hjust = 0),
      legend.title      = element_text(size = 9),
      legend.text       = element_text(size = 8.5),
      legend.key        = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = "#CCCCCC",
                                       linewidth = 0.3),
      axis.title        = element_text(size = 10)
    )
}

# ----------------------------------------------------------------------------
# 2. TF MECHANISM CLASSIFIER (used for panels A & B coloring)
# ----------------------------------------------------------------------------
# Order matters: more specific patterns are checked first.

classify_tf <- function(tf_name) {
  t <- tolower(tf_name)
  if (any(str_detect(t, c("mta3","mecp2","mbd2","mbd1","cxxc1",
                          "rcor1","prdm2","prdm4","dnmt"))))
    return("Chromatin/Epigenetic")
  if (str_detect(t, "zfx") || str_detect(t, "znf711"))
    return("Chromatin/Epigenetic")
  if (any(str_detect(t, c("e2f","tfdp","mycn","myca","mych",
                          "mnt","maz","znf143","hinfp"))))
    return("Cell cycle/Proliferation")
  if (any(str_detect(t, c("hsf","mtf1","nfkb","rest"))))
    return("Stress/Homeostasis")
  if (any(str_detect(t, c("ar","esr","pgr","rara","rarg","rarb",
                          "thra","thrb","vdr","pparab","ppar",
                          "nr5a","nr3c","nr2c","nr2f","nr4a",
                          "esrra","nr1i"))))
    return("Nuclear receptor")
  if (any(str_detect(t, c("gata1","gata2","tal1","runx1","runx3",
                          "spi1","ikzf","cebp","mafb","maff","mafg",
                          "klf1","fli1","erg","etv2","mecom","gfi1",
                          "lyl","lmo2","hhex","meis1","myb","rbpj"))))
    return("Hematopoietic")
  if (any(str_detect(t, c("myod1","myog","myf5","myf6","mef2","tcf12"))))
    return("Myogenic")
  if (any(str_detect(t, c("sox10","mitf","foxd3","tfec","snai","tfap2",
                          "pax3","pax7","sox5","sox9","her6","her8",
                          "her9","foxg","pax6","zeb"))))
    return("Neural crest/EMT")
  if (any(str_detect(t, c("smad","tcf3","tcf4","tcf7","lef1","gli",
                          "stat","foxo"))))
    return("Signaling effector")
  if (any(str_detect(t, c("ets1","ets2","etv4","etv5","etv6","etv7",
                          "elf","erf","gabpa","spdef"))))
    return("ETS family")
  if (any(str_detect(t, c("hox","pax2","pax5","pax8","pax1","wt1",
                          "lhx","hnf","ptf1","eomes","tbx","gata3",
                          "gata4","gata5","gata6","isl1","hand","nkx",
                          "six","eya","pitx","irx","dlx","msx","atoh",
                          "pou4f","prox","sp7","runx2","sox17","sox2"))))
    return("Tissue patterning")
  if (any(str_detect(t, c("yy1","klf","sp1","sp2","sp3"))))
    return("General/KLF/SP")
  "Other"
}

# ----------------------------------------------------------------------------
# 3. LOAD + ANNOTATE BOTH TSVs
# ----------------------------------------------------------------------------

read_enrichment <- function(path) {
  read_tsv(path, show_col_types = FALSE) |>
    mutate(
      log2fold      = log2_enrichment_pc,
      fdr           = pmax(FDR_BH, 1e-300),
      real_n        = real_cnes_with_tf,
      neg_log10_fdr = -log10(fdr),
      significant   = FDR_BH < FDR_CUTOFF & fold_enrichment_pc >= FOLD_CUTOFF
    ) |>
    rowwise() |>
    mutate(category = classify_tf(tf)) |>
    ungroup()
}

tf_actin <- read_enrichment(FILE_ACTIN)
tf_gnath <- read_enrichment(FILE_GNATH)

n_sig_actin <- sum(tf_actin$significant)
n_sig_gnath <- sum(tf_gnath$significant)
cat(sprintf("Loaded %d TFs (Actinopterygii: %d significant), %d TFs (Gnathostomata: %d significant)\n",
            nrow(tf_actin), n_sig_actin,
            nrow(tf_gnath), n_sig_gnath))
if (n_sig_actin <= 2 || n_sig_gnath <= 2 ||
    max(n_sig_actin, n_sig_gnath) >= 5 * max(1, min(n_sig_actin, n_sig_gnath))) {
  cat(sprintf(paste0(
    "  NOTE: strong asymmetry / very few significant TFs in one set ",
    "(Acti = %d, Gnath = %d).\n",
    "        Panel C (differential MA) may be near-empty on one side; ",
    "consider a 2-panel figure or a one-sided reframing.\n"),
    n_sig_actin, n_sig_gnath))
}

# ----------------------------------------------------------------------------
# 4. VOLCANO BUILDER (Panels A & B)
# ----------------------------------------------------------------------------
# `curated_labels` may be NULL (auto-label via ggrepel) or a tibble with
# columns `tf`, `lx`, `ly`, `hj` for manual placement with segment connectors.

build_volcano <- function(tf_data, panel_title, curated_labels = NULL,
                          x_limits = NULL) {
  sig_pts   <- tf_data |> filter( significant)
  nosig_pts <- tf_data |> filter(!significant)

  ## Data-driven x-limits unless explicitly supplied. The previous default
  ## (c(-1, 5.2)) was tuned to an earlier dataset; with the revised CNE sets
  ## the point cloud is much smaller, so let the axis frame the actual data.
  if (is.null(x_limits)) {
    xr <- range(tf_data$log2fold[is.finite(tf_data$log2fold)], na.rm = TRUE)
    x_limits <- c(floor(xr[1]), ceiling(xr[2]))
  }

  p <- ggplot() +
    geom_point(data = nosig_pts,
               aes(x = log2fold, y = neg_log10_fdr),
               color = "#DDDDDD", alpha = 0.5, size = 1.0) +
    geom_vline(xintercept = LOG2FC_CUTOFF, linetype = "dashed",
               color = "#999999", linewidth = 0.3, alpha = 0.6) +
    geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed",
               color = "#999999", linewidth = 0.3, alpha = 0.6) +
    geom_point(data = sig_pts,
               aes(x = log2fold, y = neg_log10_fdr,
                   size = real_n, fill = category),
               shape = 21, color = "black", stroke = 0.3, alpha = 0.78)

  ## ---- labels ----
  if (!is.null(curated_labels)) {
    lab_df <- curated_labels |>
      inner_join(tf_data |> select(tf, log2fold, neg_log10_fdr), by = "tf") |>
      rename(x = log2fold, y = neg_log10_fdr) |>
      mutate(display = sub("\\+.*$", "", tf))

    p <- p +
      geom_segment(data = lab_df,
                   aes(x = x, y = y, xend = lx, yend = ly),
                   color = "#666666", linewidth = 0.25) +
      geom_label(data = lab_df,
                 aes(x = lx, y = ly, label = display, hjust = hj),
                 size = 2.9, color = "#222222",
                 fill = "white", label.size = 0.2,
                 label.padding = unit(0.12, "lines"))
  } else {
    ## Auto-select: union of top-by-FDR and top-by-effect among significant.
    ## Guarded so an empty significant set (possible under the revised data,
    ## e.g. the Actinopterygii panel) does not error.
    if (nrow(sig_pts) == 0) {
      auto_lab <- sig_pts |> mutate(display = character(0))
      p <- p  # nothing to label
    } else {
      by_fdr <- sig_pts |> arrange(desc(neg_log10_fdr)) |>
                slice_head(n = N_GNATH_LABEL_BY_FDR)
      by_eff <- sig_pts |> arrange(desc(log2fold)) |>
                slice_head(n = N_GNATH_LABEL_BY_LOG2FE)
      auto_lab <- bind_rows(by_fdr, by_eff) |>
                  distinct(tf, .keep_all = TRUE) |>
                  mutate(display = sub("\\+.*$", "", tf))

      p <- p +
        geom_label_repel(
          data            = auto_lab,
          aes(x = log2fold, y = neg_log10_fdr, label = display),
          size            = 2.9,
          color           = "#222222",
          fill            = "white",
          label.size      = 0.2,
          label.padding   = unit(0.12, "lines"),
          box.padding     = 0.35,
          point.padding   = 0.25,
          segment.color   = "#666666",
          segment.size    = 0.25,
          min.segment.length = 0,
          max.overlaps    = Inf,
          seed            = 1
        )
    }
  }

  p +
    scale_fill_manual(values = palette_mech, name = "TF mechanism category",
                      drop = FALSE) +
    scale_size_continuous(range = c(1.2, 9), guide = "none") +
    scale_x_continuous(limits = x_limits) +
    labs(
      title = panel_title,
      x = expression(log[2] ~ "fold enrichment vs random genomic background"),
      y = expression(-log[10] ~ "FDR")
    ) +
    guides(fill = guide_legend(override.aes = list(size = 4), ncol = 1)) +
    panel_theme()
}

# ----------------------------------------------------------------------------
# 5. PANEL A LABELLING
# ----------------------------------------------------------------------------
# The previous version hand-placed ~21 curated labels tuned to an earlier
# dataset in which many TFs were significant in the Actinopterygii-specific
# set. Under the revised CNE sets that panel has very few significant points,
# so those curated positions no longer correspond to significant hits (and
# would draw callouts to non-significant grey points). Panel A now uses the
# same automatic top-by-FDR / top-by-effect labelling as Panel B.
curated_actin <- NULL

# ----------------------------------------------------------------------------
# 6. BUILD PANELS A & B
# ----------------------------------------------------------------------------

p_volcano_actin <- build_volcano(
  tf_data        = tf_actin,
  panel_title    = "Actinopterygii-specific CNEs",
  curated_labels = curated_actin   # now NULL -> automatic labelling
)

p_volcano_gnath <- build_volcano(
  tf_data        = tf_gnath,
  panel_title    = "Gnathostomata-specific CNEs",
  curated_labels = NULL
)

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# 7. COMPOSE: A | B
# ----------------------------------------------------------------------------
# The differential MA-style panel (previously Panel C) has been removed.
# Under the revised, exon-filtered CNE sets, only one TF motif reaches
# significance in the Actinopterygii-specific set versus sixteen in the
# Gnathostome-conserved set (with a single motif, Six6a/Six6b, significant
# in both) - too lopsided for a meaningful differential contrast between
# the two clades. This figure is now a two-panel supplementary figure
# showing the two per-clade volcano plots side by side; the asymmetry
# itself, and the shared/distinct hits, are described in the main text.

combined <- (p_volcano_actin | p_volcano_gnath) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "TF motif enrichment in Acti-spec and Gnath-cons CNEs",
    subtitle = "Per-clade volcano plots (color = TF mechanism category).",
    caption  = paste0(
      "Significance thresholds: FDR < ", FDR_CUTOFF,
      " AND fold enrichment >= ", FOLD_CUTOFF,
      " (log2FE >= log2(", FOLD_CUTOFF, ")).   ",
      "Point size proportional to # CNEs carrying the motif."
    ),
    tag_levels = "A",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 15, hjust = 0),
      plot.subtitle = element_text(size = 10.5, color = "#555555", hjust = 0),
      plot.caption  = element_text(size = 9, color = "#555555", hjust = 0),
      plot.background = element_rect(fill = "white", color = NA)
    )
  ) &
  theme(
    legend.position   = "right",
    plot.tag          = element_text(face = "bold", size = 14)
  )

# ----------------------------------------------------------------------------
# 8. EXPORT
# ----------------------------------------------------------------------------

ggsave(paste0(OUTPUT_STEM, ".pdf"), combined,
       width = 13, height = 6.5, device = cairo_pdf)

ggsave(paste0(OUTPUT_STEM, ".png"), combined,
       width = 13, height = 6.5, dpi = 160, bg = "white")

cat("\nWrote:\n",
    "  ", OUTPUT_STEM, ".pdf\n",
    "  ", OUTPUT_STEM, ".png\n", sep = "")
