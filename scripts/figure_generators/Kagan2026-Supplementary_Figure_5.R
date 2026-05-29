# ============================================================================
# CNE abundance vs gene structure — ALL genes in the two raw CNE tables.
#
# Unlike the top-N hotspot script, this analyses every gene that appears in
# each CNE table (i.e. every gene with >= 1 CNE), removing the range
# restriction that flattened the earlier top-50 correlations.
#
# For each gene:
#   n_cnes        = number of CNEs assigned to it (rows in the raw table)
#   gene_length   = genomic span of the locus, bp        } from the
#   exon_number   = exons in the longest transcript      } GRCz12tu GTF
#
# It then (1) tests whether CNE count correlates with gene length / exon
# number per set (Spearman, rank-based, robust to the heavy skew), and
# (2) visualises the relationship with binned barplots: genes are grouped
# into gene-length / exon-number bins and the bar height is the mean CNE
# count per gene in that bin. With thousands of genes a per-gene scatter is
# an unreadable hairball; binned means show the trend cleanly.
#
# Inputs (set paths in CONFIG):
#   actinopterygii_cne_final_table.tsv
#   gnathostomata_cne_final_table.tsv
#   GRCz12tu GTF  (NCBI RefSeq, GCF_049306965.1)
#
# Outputs:
#   cne_vs_gene_structure_allgenes.tsv   per-gene table
#   cne_vs_gene_structure_allgenes.pdf   two-panel binned barplot
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ============================================================================
# 0. CONFIG
# ============================================================================
ACTI_FILE  <- "actinopterygii_cne_final_table.tsv"
GNATH_FILE <- "gnathostomata_cne_final_table.tsv"
GTF_FILE   <- "Danio_rerio.GRCz12tu.gtf"

EXCLUDE_LOC            <- FALSE  # FALSE = keep LOC* gene models (analyse all data)
INCLUDE_ZERO_CNE_GENES <- FALSE  # TRUE = add every GTF gene with 0 CNEs to each
                                 # set (a genome-wide background; the raw files
                                 # themselves contain only genes with >= 1 CNE)


# bin edges (edit freely). cut(): right-closed, lowest included.
LENGTH_BREAKS <- c(0, 5e3, 1e4, 2.5e4, 5e4, 1e5, 2.5e5, 5e5, Inf)
LENGTH_LABELS <- c("<5kb","5-10kb","10-25kb","25-50kb",
                   "50-100kb","100-250kb","250-500kb",">500kb")
EXON_BREAKS   <- c(0, 1, 2, 5, 10, 20, 50, 100, Inf)
EXON_LABELS   <- c("1","2","3-5","6-10","11-20","21-50","51-100",">100")

OUT_TSV <- "cne_vs_gene_structure_allgenes.tsv"
OUT_PDF <- "cne_vs_gene_structure_allgenes.pdf"

ACTI_COL  <- "#CC79A7"
GNATH_COL <- "#009E73"
palette_pair <- c(`Actinopterygii-specific` = ACTI_COL,
                  `Gnathostome-conserved`   = GNATH_COL)

# ============================================================================
# 1. CNE count per gene, per set (all genes in the raw tables)
# ============================================================================
count_cnes <- function(path, set_label) {
  d <- readr::read_tsv(path, show_col_types = FALSE,
                       guess_max = 100000, na = c("", "NA"))
  d |>
    dplyr::mutate(gene = tolower(sub("^gene-", "", as.character(geneId)))) |>
    dplyr::filter(!is.na(gene), gene != "") |>
    { \(x) if (EXCLUDE_LOC)
             dplyr::filter(x, !grepl("^loc", gene, ignore.case = TRUE))
           else x }() |>
    dplyr::count(gene, name = "n_cnes") |>
    dplyr::mutate(set = set_label)
}
counts <- dplyr::bind_rows(
  count_cnes(ACTI_FILE,  "Actinopterygii-specific"),
  count_cnes(GNATH_FILE, "Gnathostome-conserved")
)
cat(sprintf("Genes with >=1 CNE: acti %d, gnath %d.\n",
            sum(counts$set == "Actinopterygii-specific"),
            sum(counts$set == "Gnathostome-conserved")))

# ============================================================================
# 2. Gene length + exon number from the GRCz12tu GTF
# ============================================================================
if (!file.exists(GTF_FILE))
  stop("GTF file not found: ", GTF_FILE, " — set GTF_FILE in the CONFIG block.")

make_txdb <- function(gtf) {
  if (requireNamespace("txdbmaker", quietly = TRUE))
    txdbmaker::makeTxDbFromGFF(gtf, format = "auto")
  else
    GenomicFeatures::makeTxDbFromGFF(gtf, format = "auto")
}
cat("Parsing GTF (this can take a minute) ...\n")
txdb <- make_txdb(GTF_FILE)

g_gr     <- GenomicFeatures::genes(txdb)
gene_len <- data.frame(gene_id        = names(g_gr),
                       gene_length_bp = BiocGenerics::width(g_gr),
                       stringsAsFactors = FALSE)

txl <- GenomicFeatures::transcriptLengths(txdb)   # cols incl. gene_id, tx_len, nexon
gene_exon <- txl |>
  dplyr::filter(!is.na(gene_id)) |>
  dplyr::group_by(gene_id) |>
  dplyr::slice_max(tx_len, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(gene_id, exon_number = nexon)

gene_struct <- gene_len |>
  dplyr::left_join(gene_exon, by = "gene_id") |>
  dplyr::mutate(key = tolower(sub("^gene-", "", gene_id)))

# ============================================================================
# 3. Optionally add zero-CNE genes, then attach gene structure
# ============================================================================
if (INCLUDE_ZERO_CNE_GENES) {
  zero <- expand.grid(
    gene = unique(gene_struct$key),
    set  = c("Actinopterygii-specific", "Gnathostome-conserved"),
    stringsAsFactors = FALSE
  ) |>
    dplyr::anti_join(counts, by = c("gene", "set")) |>
    dplyr::mutate(n_cnes = 0L)
  counts <- dplyr::bind_rows(counts, zero)
  cat("Added zero-CNE genes: background = all GTF genes per set.\n")
}

per_gene <- counts |>
  dplyr::left_join(gene_struct, by = c("gene" = "key")) |>
  dplyr::mutate(
    found_in_gtf = !is.na(gene_length_bp),
    length_bin   = cut(gene_length_bp, breaks = LENGTH_BREAKS,
                       labels = LENGTH_LABELS, include.lowest = TRUE),
    exon_bin     = cut(exon_number,    breaks = EXON_BREAKS,
                       labels = EXON_LABELS, include.lowest = TRUE)
  ) |>
  dplyr::select(set, gene, n_cnes, gene_id, gene_length_bp, exon_number,
                length_bin, exon_bin, found_in_gtf)

n_missing <- sum(!per_gene$found_in_gtf)
if (n_missing > 0)
  cat(sprintf("WARNING: %d gene(s) not matched in the GTF (dropped from stats).\n",
              n_missing))

readr::write_tsv(per_gene, OUT_TSV)
cat(sprintf("Wrote %s\n", OUT_TSV))

matched <- per_gene |> dplyr::filter(found_in_gtf)

# ============================================================================
# 4. Correlation: CNE count vs gene length / exon number (Spearman, per set)
# ============================================================================
cat("\n=== Spearman correlations (CNE count vs gene structure) ===\n")
corr_tbl <- matched |>
  dplyr::group_by(set) |>
  dplyr::group_modify(function(d, key) {
    ct_len  <- suppressWarnings(cor.test(d$n_cnes, d$gene_length_bp,
                                         method = "spearman"))
    ct_exon <- suppressWarnings(cor.test(d$n_cnes, d$exon_number,
                                         method = "spearman"))
    data.frame(n_genes      = nrow(d),
               rho_length   = unname(ct_len$estimate),
               p_length     = ct_len$p.value,
               rho_exon     = unname(ct_exon$estimate),
               p_exon       = ct_exon$p.value)
  }) |>
  dplyr::ungroup()
print(as.data.frame(corr_tbl), row.names = FALSE, digits = 3)

rho_txt <- function(set_name, col)
  sprintf("%s rho=%.2f", sub("-.*", "", set_name),
          corr_tbl[[col]][corr_tbl$set == set_name])


# ============================================================================
# 5. Two-panel binned barplot
# ============================================================================
box_panel <- function(df, bin_col, bin_levels, x_title, subtitle) {
  d <- df |>
    dplyr::filter(!is.na(.data[[bin_col]])) |>
    dplyr::mutate(bin = factor(.data[[bin_col]], levels = bin_levels))
  ggplot(d, aes(x = bin, y = n_cnes, fill = set)) +
    geom_boxplot(position = position_dodge2(preserve = "single"),
                 width = 0.7, linewidth = 0.3, colour = "grey30",
                 outlier.size = 0.5, outlier.alpha = 0.25) +
    scale_fill_manual(values = palette_pair, name = NULL) +
    scale_y_log10() +
    labs(title = x_title, subtitle = subtitle,
         x = NULL, y = "CNEs per gene (log scale)") +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, colour = "grey30"),
          axis.text.x   = element_text(angle = 35, hjust = 1),
          panel.grid.major.x = element_blank())
}

pA <- box_panel(matched, "length_bin", LENGTH_LABELS, "By gene length",
                sprintf("Spearman:  %s,  %s",
                        rho_txt("Actinopterygii-specific", "rho_length"),
                        rho_txt("Gnathostome-conserved",   "rho_length")))
pB <- box_panel(matched, "exon_bin", EXON_LABELS, "By exon number",
                sprintf("Spearman:  %s,  %s",
                        rho_txt("Actinopterygii-specific", "rho_exon"),
                        rho_txt("Gnathostome-conserved",   "rho_exon")))

fig <- (pA / pB) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "CNE abundance vs gene structure — all CNE-associated genes",
    theme = theme(plot.title = element_text(face = "bold", size = 13))) &
  theme(legend.position = "bottom")

ggsave(OUT_PDF, fig, width = 9, height = 8.5, bg = "white")
cat(sprintf("Wrote %s\n", OUT_PDF))

# ============================================================================
# 6. Density scatter + trend line + hotspot overlay
#    geom_rect with manual log-space binning handles overplotting; the lm line
#    shows the weak positive trend per set; the labelled points are the top-10
#    CNE hotspots, which sit well ABOVE the trend — i.e. the hotspots are not
#    explained by gene structure.
#
#    Two palette modes (USE_SET_PALETTE):
#      FALSE: shared viridis fill (one legend)
#      TRUE : grey -> ACTI_COL / grey -> GNATH_COL per set; two gradient
#             legends are attached via dummy plots and patchwork.
# ============================================================================
suppressPackageStartupMessages(library(cowplot))

OUT_PDF2 <- "cne_vs_gene_structure_density.pdf"

USE_SET_PALETTE <- TRUE     # TRUE = per-set grey->colour ramps
GREY_LOW        <- "#888888"

top10 <- matched |>
  dplyr::group_by(set) |>
  dplyr::slice_max(n_cnes, n = 10, with_ties = FALSE) |>
  dplyr::ungroup()

# ---- hotspot labels: a function so x/y can vary by panel -------------------
lab_layer <- function(xv) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(
      data = top10,
      aes(x = .data[[xv]], y = n_cnes, label = gene),
      size = 2.4, max.overlaps = 20, segment.size = 0.2,
      colour = "black", fontface = "italic")
  } else {
    geom_text(
      data = top10,
      aes(x = .data[[xv]], y = n_cnes, label = gene),
      size = 2.4, vjust = -0.7, colour = "black", fontface = "italic")
  }
}

# ---- main panel builder ----------------------------------------------------
density_panel <- function(xvar, x_lab, rho_col) {
  if (USE_SET_PALETTE) {
    N_BINS   <- 40
    x_breaks <- seq(min(log10(matched[[xvar]])),
                    max(log10(matched[[xvar]])), length.out = N_BINS + 1)
    y_breaks <- seq(min(log10(matched$n_cnes)),
                    max(log10(matched$n_cnes)), length.out = N_BINS + 1)

    binned <- matched |>
      dplyr::filter(.data[[xvar]] > 0, n_cnes > 0) |>
      dplyr::mutate(
        xi = findInterval(log10(.data[[xvar]]), x_breaks, all.inside = TRUE),
        yi = findInterval(log10(n_cnes),        y_breaks, all.inside = TRUE)) |>
      dplyr::group_by(set, xi, yi) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
      dplyr::mutate(
        xmin = 10 ^ x_breaks[xi],     xmax = 10 ^ x_breaks[xi + 1],
        ymin = 10 ^ y_breaks[yi],     ymax = 10 ^ y_breaks[yi + 1])

    ramp <- function(hi) colorRampPalette(c(GREY_LOW, hi))(100)
    pal  <- list(`Actinopterygii-specific` = ramp(ACTI_COL),
                 `Gnathostome-conserved`   = ramp(GNATH_COL))

    binned <- binned |>
      dplyr::group_by(set) |>
      dplyr::mutate(
        idx  = as.integer(ceiling(100 * rank(n, ties.method = "average") /
                                  dplyr::n())),
        idx  = pmin(pmax(idx, 1L), 100L),
        fill = pal[[as.character(set[1])]][idx]) |>
      dplyr::ungroup()

    base <- ggplot(binned) +
      geom_rect(aes(xmin = xmin, xmax = xmax,
                    ymin = ymin, ymax = ymax, fill = fill)) +
      scale_fill_identity()
  } else {
    base <- ggplot(matched, aes(x = .data[[xvar]], y = n_cnes)) +
      geom_bin2d(bins = 40) +
      scale_fill_viridis_c(name = "genes", trans = "log10")
  }

  base +
    geom_smooth(data = matched, aes(x = .data[[xvar]], y = n_cnes),
                method = "lm", se = FALSE, linewidth = 0.7,
                colour = "#D55E00") +
    geom_point(data = top10, aes(x = .data[[xvar]], y = n_cnes),
               colour = "black", size = 1.1) +
    lab_layer(xvar) +
    scale_x_log10(labels = scales::label_comma()) +
    scale_y_log10() +
    facet_wrap(~ set) +
    labs(x = x_lab, y = "CNEs per gene (log scale)",
         subtitle = sprintf("Spearman:  %s,  %s",
                            rho_txt("Actinopterygii-specific", rho_col),
                            rho_txt("Gnathostome-conserved",   rho_col))) +
    theme_minimal(base_size = 11) +
    theme(plot.subtitle = element_text(size = 9, colour = "grey30"),
          strip.text    = element_text(face = "bold"))
}

# ---- two dummy-plot legends for the per-set gradients ----------------------
legend_for <- function(label, hi_col, n_vec) {
  p <- ggplot(data.frame(x = 1, y = 1, n = n_vec),
              aes(x, y, fill = n)) +
    geom_tile() +
    scale_fill_gradient(name = label, low = GREY_LOW, high = hi_col,
                        guide = guide_colourbar(direction = "horizontal",
                                                title.position = "top",
                                                barwidth  = grid::unit(4, "cm"),
                                                barheight = grid::unit(0.4, "cm"))) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  if (utils::packageVersion("ggplot2") >= "4.0.0") {
    g  <- ggplotGrob(p)
    gb <- g$grobs[[which(g$layout$name == "guide-box-bottom")]]
    if (is.null(gb) || inherits(gb, "zeroGrob"))
      gb <- g$grobs[[which(g$layout$name == "guide-box")]]
    gb
  } else if (utils::packageVersion("cowplot") >= "1.1.3") {
    cowplot::get_plot_component(p, "guide-box-bottom", return_all = TRUE)[[1]]
  } else {
    cowplot::get_legend(p)
  }
}

# ---- build and combine the two panels --------------------------------------
pL <- density_panel("gene_length_bp", "gene length (bp, log scale)", "rho_length")
pE <- density_panel("exon_number",    "exon number (log scale)",     "rho_exon")

fig2 <- (pL / pE) +
  plot_annotation(
    title = "CNE abundance vs gene structure - gene density, with hotspots labelled",
    theme = theme(plot.title = element_text(face = "bold", size = 13)))

if (USE_SET_PALETTE) {
  n_acti  <- (matched |> dplyr::filter(set == "Actinopterygii-specific"))$n_cnes
  n_gnath <- (matched |> dplyr::filter(set == "Gnathostome-conserved"))$n_cnes
  leg_acti  <- legend_for("Actinopterygii CNEs per gene", ACTI_COL,  n_acti)
  leg_gnath <- legend_for("Gnathostome CNEs per gene",    GNATH_COL, n_gnath)

  legend_row <- cowplot::plot_grid(leg_acti, leg_gnath, nrow = 1)
  fig2 <- cowplot::plot_grid(fig2, legend_row,
                             ncol = 1, rel_heights = c(10, 1))
}

ggsave(OUT_PDF2, fig2,
       width  = 10,
       height = if (USE_SET_PALETTE) 10.5 else 9, bg = "white")
cat(sprintf("Wrote %s\n", OUT_PDF2))

n_gnath <- (matched |> dplyr::filter(set == "Gnathostome-conserved"))$n_cnes
leg_gnath <- legend_for("Gnathostome\nCNEs per gene", GNATH_COL, n_gnath)
legends   <- patchwork::wrap_elements(leg_acti) /
             patchwork::wrap_elements(leg_gnath)
print(legends)
combined  <- (pL / pE | legends) + patchwork::plot_layout(widths = c(8, 1))
print(combined)
ggsave("test_combined.pdf", combined, width = 12, height = 9)