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
  library(colorspace)
})

# ----------------------------------------------------------------------------
# 0. CONFIG
# ----------------------------------------------------------------------------

FILE_ACTIN <- "actinopterygii_vs_random_genomic.tf_enrichment_real_vs_random_genomic.tsv"
FILE_GNATH <- "gnathostomata_vs_random_genomic.tf_enrichment_real_vs_random_genomic.tsv"

FDR_CUTOFF      <- 0.05
FOLD_CUTOFF     <- 1.5
LOG2FC_CUTOFF   <- log2(FOLD_CUTOFF)
SHARED_ABS_DIFF <- 0.5            # |Δlog2FE| below this and sig-in-both ⇒ shared

## MA panel label selection.
N_TOP_DIFF_BY_ABS    <- 14
N_TOP_DIFF_BY_VOLUME <- 12

## Volcano panel B auto-label selection (panel A uses curated positions below).
N_GNATH_LABEL_BY_FDR    <- 12
N_GNATH_LABEL_BY_LOG2FE <- 8

OUTPUT_STEM <- "TF_enrichment_combined_figure"

# ----------------------------------------------------------------------------
# 1. HELPERS, PALETTES
# ----------------------------------------------------------------------------

cap_first   <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
clean_label <- function(x) cap_first(sub("\\+.*$", "", x))

## Clade anchor colors — used in panel C's differential classes.
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

## Panel C: differential classes between the two clades.
palette_diff <- c(
  "Shared"                       = "#9A9A9A",
  "Stronger in Actinopterygii"   = colorspace::lighten(ACTI,  0.35),
  "Actinopterygii-specific"      = ACTI,
  "Stronger in Gnathostomata"    = colorspace::lighten(GNATH, 0.35),
  "Gnathostomata-specific"       = GNATH
)
diff_class_levels <- names(palette_diff)

## Shared theme for all three panels — tweak once, apply everywhere.
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

cat(sprintf("Loaded %d TFs (Actinopterygii: %d significant), %d TFs (Gnathostomata: %d significant)\n",
            nrow(tf_actin), sum(tf_actin$significant),
            nrow(tf_gnath), sum(tf_gnath$significant)))

# ----------------------------------------------------------------------------
# 4. VOLCANO BUILDER (Panels A & B)
# ----------------------------------------------------------------------------
# `curated_labels` may be NULL (auto-label via ggrepel) or a tibble with
# columns `tf`, `lx`, `ly`, `hj` for manual placement with segment connectors.

build_volcano <- function(tf_data, panel_title, curated_labels = NULL,
                          x_limits = c(-1, 5.2)) {
  sig_pts   <- tf_data |> filter( significant)
  nosig_pts <- tf_data |> filter(!significant)

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
# 5. CURATED LABELS FOR PANEL A (Actinopterygii)
# ----------------------------------------------------------------------------
# Hand-tuned positions from the original standalone volcano script.

curated_actin <- tibble::tribble(
  ~tf,             ~lx,   ~ly,    ~hj,
  "rest",          4.50,  2.7,    0,
  "snai1a",        3.10,  18.0,   0,
  "prdm2a",        2.35,  60.5,   0,
  "mta3",          1.95,  36.5,   0,
  "myf6+ttf1.4",   2.45,  45.5,   0,
  "myog",          2.55,  18.5,   0,
  "etv2",          2.55,  15.5,   0,
  "tal1",          0.55,  39.4,   1,
  "pax2a+pax5",    0.50,  35.0,   1,
  "yy1a+yy1b",     0.50,  23.2,   1,
  "smad1",         0.50,  21.0,   1,
  "smad2",         1.55,  22.0,   0,
  "mafb+nrl",      1.50,  26.2,   0,
  "nkx2.1",        0.50,  17.5,   1,
  "klf12b",        0.50,  16.0,   1,
  "ar",            2.30,  11.5,   0,
  "VDR+nr1i2",     2.55,  4.0,    0,
  "tfap2b",        2.40,  5.6,    0,
  "e2f4",          2.65,  7.5,    0,
  "hsf1+hsf4",     2.45,  8.5,    0,
  "mecp2",         1.45,  9.0,    0
)

# ----------------------------------------------------------------------------
# 6. BUILD PANELS A & B
# ----------------------------------------------------------------------------

p_volcano_actin <- build_volcano(
  tf_data        = tf_actin,
  panel_title    = "Actinopterygii-specific CNEs",
  curated_labels = curated_actin
)

p_volcano_gnath <- build_volcano(
  tf_data        = tf_gnath,
  panel_title    = "Gnathostomata-specific CNEs",
  curated_labels = NULL
)

# ----------------------------------------------------------------------------
# 7. PANEL C — MA-STYLE DIFFERENTIAL PLOT
# ----------------------------------------------------------------------------

diff_table <- inner_join(
  tf_actin |> select(tf,
                     log2fold_actin = log2fold, fdr_actin = fdr,
                     sig_actin      = significant, real_n_actin = real_n),
  tf_gnath |> select(tf,
                     log2fold_gnath = log2fold, fdr_gnath = fdr,
                     sig_gnath      = significant, real_n_gnath = real_n),
  by = "tf"
) |>
  mutate(
    ## MA-layout coordinates: mean and difference.
    mean_log2fold  = (log2fold_actin + log2fold_gnath) / 2,
    delta_log2fold = log2fold_actin - log2fold_gnath,
    abs_diff       = abs(delta_log2fold),
    max_real_n     = pmax(real_n_actin, real_n_gnath, na.rm = TRUE),
    diff_class     = case_when(
      sig_actin &  sig_gnath & abs_diff < SHARED_ABS_DIFF ~ "Shared",
      sig_actin &  sig_gnath & delta_log2fold > 0         ~ "Stronger in Actinopterygii",
      sig_actin &  sig_gnath & delta_log2fold < 0         ~ "Stronger in Gnathostomata",
      sig_actin & !sig_gnath                              ~ "Actinopterygii-specific",
      sig_gnath & !sig_actin                              ~ "Gnathostomata-specific",
      TRUE                                                ~ "n.s. in both"
    )
  )

cat(sprintf("Joined %d TFs for differential panel.\n", nrow(diff_table)))
print(diff_table |> count(diff_class, sort = TRUE), n = Inf)

sig_only <- diff_table |>
  filter(diff_class != "n.s. in both") |>
  mutate(diff_class = factor(diff_class, levels = diff_class_levels))
noise <- diff_table |> filter(diff_class == "n.s. in both")

## --- label selection: high-skew ∪ high-volume -------------------------------
by_diff <- sig_only |> arrange(desc(abs_diff))   |> slice_head(n = N_TOP_DIFF_BY_ABS)
by_vol  <- sig_only |> arrange(desc(max_real_n)) |> slice_head(n = N_TOP_DIFF_BY_VOLUME)
labels_df <- bind_rows(by_diff, by_vol) |>
  distinct(tf, .keep_all = TRUE) |>
  mutate(display = clean_label(tf))

## --- axes -------------------------------------------------------------------
x_rng <- range(sig_only$mean_log2fold, na.rm = TRUE)
x_lo  <- floor(x_rng[1] * 2) / 2 - 0.3
x_hi  <- ceiling(x_rng[2] * 2) / 2 + 0.3
y_abs <- max(abs(sig_only$delta_log2fold), na.rm = TRUE)
y_hi  <- ceiling(y_abs * 2) / 2 + 0.3
y_lo  <- -y_hi

## --- per-class counts for corner annotations --------------------------------
cl_counts <- sig_only |> count(diff_class, .drop = FALSE)
get_n <- function(lbl) {
  v <- cl_counts$n[cl_counts$diff_class == lbl]
  if (length(v) == 0) 0L else v
}
n_actin_specific <- get_n("Actinopterygii-specific")
n_actin_stronger <- get_n("Stronger in Actinopterygii")
n_gnath_specific <- get_n("Gnathostomata-specific")
n_gnath_stronger <- get_n("Stronger in Gnathostomata")
n_shared         <- get_n("Shared")

corner_anno <- tibble::tribble(
  ~x,            ~y,           ~hj, ~vj, ~label,
  x_lo + 0.15,   y_hi - 0.15,  0,   1,
    sprintf("Actinopterygii-leaning\nspecific: %d   stronger: %d",
            n_actin_specific, n_actin_stronger),
  x_lo + 0.15,   y_lo + 0.15,  0,   0,
    sprintf("Gnathostomata-leaning\nspecific: %d   stronger: %d",
            n_gnath_specific, n_gnath_stronger)
)
shared_anno <- tibble::tibble(
  x = x_hi - 0.15, y = 0,
  label = sprintf("Shared (n = %d)", n_shared)
)

p_diff_ma <- ggplot() +
  geom_point(data = noise,
             aes(x = mean_log2fold, y = delta_log2fold),
             color = "#E2E2E2", size = 0.6, alpha = 0.45, stroke = 0) +
  geom_hline(yintercept = 0,
             linetype = "dashed", color = "#777777", linewidth = 0.45) +
  geom_hline(yintercept = c(-SHARED_ABS_DIFF, SHARED_ABS_DIFF),
             linetype = "dotted", color = "#BBBBBB", linewidth = 0.3) +
  geom_vline(xintercept = LOG2FC_CUTOFF,
             linetype = "dotted", color = "#BBBBBB", linewidth = 0.3) +
  geom_point(
    data = sig_only |> arrange(diff_class),
    aes(x = mean_log2fold, y = delta_log2fold,
        fill = diff_class, size = max_real_n),
    shape = 21, color = "white", stroke = 0.35, alpha = 0.92
  ) +
  geom_text_repel(
    data            = labels_df,
    aes(x = mean_log2fold, y = delta_log2fold, label = display),
    size            = 3.0,
    color           = "#222222",
    box.padding     = 0.45,
    point.padding   = 0.30,
    segment.color   = "#888888",
    segment.size    = 0.25,
    min.segment.length = 0,
    max.overlaps    = Inf,
    seed            = 1
  ) +
  geom_text(data = corner_anno,
            aes(x = x, y = y, label = label, hjust = hj, vjust = vj),
            size = 3.2, color = "#555555", fontface = "italic",
            lineheight = 0.95) +
  geom_text(data = shared_anno,
            aes(x = x, y = y, label = label),
            hjust = 1, vjust = -0.6,
            size = 3.2, color = "#555555", fontface = "italic") +
  scale_fill_manual(values = palette_diff,
                    name   = "Differential class",
                    drop   = FALSE) +
  scale_size_continuous(range  = c(2, 11),
                        name   = "Max # CNEs with motif",
                        breaks = c(50, 200, 500, 1000, 1500)) +
  coord_cartesian(xlim = c(x_lo, x_hi), ylim = c(y_lo, y_hi)) +
  labs(
    title    = "Differential view (MA layout)",
    subtitle = "x: mean enrichment across CNE sets   |   y: clade-bias axis",
    x = expression("Mean " * log[2] * " fold enrichment, (Actinopterygii + Gnathostomata) / 2"),
    y = expression(Delta * log[2] * "FE  (Actinopterygii − Gnathostomata)")
  ) +
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(size = 5, alpha = 1,
                          shape = 21, color = "white", stroke = 0.35)
    ),
    size = guide_legend(
      order = 2,
      override.aes = list(shape = 21,
                          fill  = "grey55",
                          color = "grey25",
                          stroke = 0.4,
                          alpha = 1)
    )
  ) +
  panel_theme()

# ----------------------------------------------------------------------------
# 8. COMPOSE: (A | B) / C
# ----------------------------------------------------------------------------

combined <- (p_volcano_actin | p_volcano_gnath) / p_diff_ma +
  plot_layout(heights = c(1, 1.05), guides = "collect") +
  plot_annotation(
    title    = "Differential TF motif usage in lineage-specific CNEs",
    subtitle = paste0(
      "Panels A & B: per-clade volcanos (color = TF mechanism category).   ",
      "Panel C: differential view in MA coordinates (color = differential class)."
    ),
    caption  = paste0(
      "Significance thresholds: FDR < ", FDR_CUTOFF,
      " AND fold enrichment ≥ ", FOLD_CUTOFF,
      " (log2FE ≥ log2(", FOLD_CUTOFF, ")).   ",
      "Point size ∝ # CNEs carrying the motif.   ",
      "Panel C |Δlog2FE| threshold for ‘Shared’ vs ‘Stronger’: ", SHARED_ABS_DIFF, "."
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
# 9. EXPORT
# ----------------------------------------------------------------------------

ggsave(paste0(OUTPUT_STEM, ".pdf"), combined,
       width = 15, height = 13, device = cairo_pdf)

ggsave(paste0(OUTPUT_STEM, ".png"), combined,
       width = 15, height = 13, dpi = 160, bg = "white")

cat("\nWrote:\n",
    "  ", OUTPUT_STEM, ".pdf\n",
    "  ", OUTPUT_STEM, ".png\n", sep = "")
