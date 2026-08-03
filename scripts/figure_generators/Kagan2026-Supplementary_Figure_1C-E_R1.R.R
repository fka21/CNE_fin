#!/usr/bin/env Rscript
# =============================================================================
# Visualization of Actinopterygii-specific AND Gnathostome-conserved CNEs
# around selected hotspot genes (default: sall1a, zfhx4, ebf3a)
# -----------------------------------------------------------------------------
# Produces one panel per gene, each showing three tracks sharing an x-axis:
#   * "genes"      : gene models in the window, one gene per row (stacked
#                    "under each other", not packed). The target gene is
#                    highlighted and labelled; other genes are shown in grey
#                    for context.
#   * "Gnath-cons" : Gnathostome-conserved CNEs, drawn as boxes (green).
#   * "Acti-spec"  : Actinopterygii-specific CNEs, drawn as boxes (pink).
#
# UNLIKE plot_pcdh_cne_panels.R, this script does NOT require you to look up
# and hardcode chromosome accessions/windows by hand. You give it gene
# names; it finds each gene's coordinates directly in the GTF and builds a
# window of [gene_start - flank_bp, gene_end + flank_bp] automatically. This
# avoids the accession/window typos that came up with the manually-specified
# pcdh regions table.
#
# INPUT FILES (edit the paths below, or run the script from the folder that
# contains them):
#   Danio_rerio.GRCz12tu.gtf                          -> full gene annotation
#   actinopterygii_cne_final_table_R1_20260731.bed    -> Acti-spec CNEs (BED6)
#   gnathostomata_cne_final_table_R1_20260731.bed      -> Gnath-cons CNEs (BED6)
# (Rename these variables below if your local BED files are named
# differently - e.g. if you exported them yourselves under other names.)
#
# OUTPUT: hotspot_gene_panels.pdf  and  .png
#
# *** CHROMOSOME NAMING ***
# This GTF's column 1 uses NCBI RefSeq accessions WITH the version suffix
# (e.g. "NC_133182.1"), while the BED files use the accession WITHOUT the
# suffix (e.g. "NC_133182") - same convention already established for
# plot_pcdh_cne_panels.R. The script strips the suffix internally so the
# two sources match; no manual chromosome mapping is needed for this
# gene-name-based lookup approach.
# =============================================================================

## ---- 0. Packages -----------------------------------------------------------
required <- c("ggplot2", "dplyr", "patchwork", "scales", "ggrepel", "data.table")
missing  <- required[!vapply(required, requireNamespace,
                             FUN.VALUE = logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
  library(scales);  library(ggrepel); library(data.table)
})

## ---- 1. Configuration ------------------------------------------------------
gtf_file    <- "Danio_rerio.GRCz12tu.gtf"
acti_bed    <- "actinopterygii_specific_drer_R1_20260731.bed"   # Acti-spec CNEs
gnath_bed   <- "gnathostomata_conserved_drer_R1_20260731.bed"    # Gnath-cons CNEs
out_prefix  <- "hotspot_gene_panels"

target_genes <- c("sall1a", "zfhx4", "ebf3a")   # <- edit this list as needed
flank_bp     <- 60000                            # bp of context on each side
                                                  #    of the gene body

# Colours (kept consistent with plot_pcdh_cne_panels.R)
acti_fill    <- "#CC79A7"; acti_outline  <- "#9E5A82"  # Actinopterygii-specific
gnath_fill   <- "#009E73"; gnath_outline <- "#00755A"  # Gnathostome-conserved
target_fill  <- "#E69F00"                              # highlighted target gene
other_fill   <- "#BFBFBF"                              # non-target genes (context)

## ---- 2. Read CNEs (BED files) ----------------------------------------------
read_cne_bed <- function(path, label) {
  if (!file.exists(path)) stop("BED file not found: ", path,
                               "\n  Set the corresponding path variable or ",
                               "setwd() to its folder.")
  d <- read.table(path, sep = "\t", header = FALSE, quote = "",
                   fill = TRUE, stringsAsFactors = FALSE,
                   col.names = c("chrom", "start", "end",
                                 "name", "score", "strand"))
  d$chrom <- sub("\\.[0-9]+$", "", d$chrom)   # drop RefSeq version suffix
  d$set   <- label
  d
}
cne_acti  <- read_cne_bed(acti_bed,  "Actinopterygii-specific")
cne_gnath <- read_cne_bed(gnath_bed, "Gnathostome-conserved")

## ---- 3. Read gene models (GTF file) ----------------------------------------
if (!file.exists(gtf_file)) stop("GTF file not found: ", gtf_file,
                                 "\n  Set 'gtf_file' or setwd() to its folder.")

message("Reading GTF (this may take a little while for a >1 GB file)...")

gtf_raw <- fread(gtf_file, header = FALSE, sep = "\t", quote = "",
                 select = c(1, 3, 4, 5, 7, 9),
                 col.names = c("chrom", "feature", "start", "end",
                               "strand", "attribute"),
                 showProgress = TRUE)

gtf_genes_raw <- gtf_raw[feature == "gene"]
rm(gtf_raw); invisible(gc())
message(sprintf("Retained %d 'gene' feature rows genome-wide.", nrow(gtf_genes_raw)))

gtf_attr <- function(attr_str, key) {
  pat <- paste0(key, ' "([^"]*)"')
  m   <- regexpr(pat, attr_str)
  out <- rep(NA_character_, length(attr_str))
  hit <- m > -1
  out[hit] <- sub(paste0('.*', key, ' "([^"]*)".*'), "\\1", attr_str[hit])
  out
}
gtf_genes_raw[, gene_id   := gtf_attr(attribute, "gene_id")]
gtf_genes_raw[, gene_name := gtf_attr(attribute, "gene_name")]
gtf_genes_raw[, attribute := NULL]
gtf_genes_raw[, chrom_noversion := sub("\\.[0-9]+$", "", chrom)]

## ---- 4. Locate each target gene and build its window -----------------------
locate_gene <- function(gene) {
  hit <- gtf_genes_raw[tolower(gene_name) == tolower(gene) |
                        tolower(gene_id)   == tolower(gene)]
  if (nrow(hit) == 0) {
    warning("Gene '", gene, "' was not found (by gene_name or gene_id) in ",
            gtf_file, " - it will be skipped. Check spelling/case, or ",
            "whether this GTF uses a different symbol for it.")
    return(NULL)
  }
  if (nrow(hit) > 1) {
    message("NOTE: '", gene, "' matched ", nrow(hit),
            " GTF rows (e.g. a gene plus a duplicate/paralog record); ",
            "using the widest span across all matches.")
  }
  data.frame(
    cluster   = gene,
    chrom     = hit$chrom_noversion[1],
    win_start = max(0, min(hit$start) - flank_bp),
    win_end   = max(hit$end) + flank_bp,
    stringsAsFactors = FALSE
  )
}

regions <- do.call(rbind, lapply(target_genes, locate_gene))
if (is.null(regions) || nrow(regions) == 0) {
  stop("None of the requested target_genes were found in the GTF - nothing to plot.")
}
message("Resolved gene windows:")
for (i in seq_len(nrow(regions))) {
  message(sprintf("  %-10s %s: %.0f-%.0f bp (%.0f kb window)",
                  regions$cluster[i], regions$chrom[i],
                  regions$win_start[i], regions$win_end[i],
                  (regions$win_end[i] - regions$win_start[i]) / 1e3))
}

## ---- 5. Helper functions ---------------------------------------------------

gene_arrow_polys <- function(g, head_bp, body_hh, head_hh) {
  out <- vector("list", nrow(g))
  for (i in seq_len(nrow(g))) {
    x1 <- g$gstart_d[i]; x2 <- g$gend_d[i]; y <- g$row_y[i]
    hw <- min(head_bp, x2 - x1)
    if (g$strand[i] == "+") {
      bx <- x2 - hw
      xs <- c(x1, bx, bx, x2, bx, bx, x1)
    } else {
      bx <- x1 + hw
      xs <- c(x2, bx, bx, x1, bx, bx, x2)
    }
    ys <- c(y - body_hh, y - body_hh, y - head_hh, y,
            y + head_hh, y + body_hh, y + body_hh)
    out[[i]] <- data.frame(grp = i, x = xs, y = ys,
                           category = g$category[i], stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

## ---- 6. Build one panel -----------------------------------------------------
build_panel <- function(r) {
  target_chrom <- r$chrom; ws <- r$win_start; we <- r$win_end
  span <- we - ws
  row_step <- 1.15
  body_hh  <- 0.22; head_hh <- 0.40

  ## --- CNEs in the window, both sets ---
  get_cne <- function(d) {
    d <- d[d$chrom == target_chrom & d$end >= ws & d$start <= we, , drop = FALSE]
    minw <- span * 0.0022
    mid  <- (d$start + d$end) / 2
    d$xmin <- pmax(ws, pmin(d$start, mid - minw / 2))
    d$xmax <- pmin(we, pmax(d$end,   mid + minw / 2))
    d
  }
  acti  <- get_cne(cne_acti)
  gnath <- get_cne(cne_gnath)

  # Two stacked CNE rows below the gene track: Gnath-cons above Acti-spec
  gnath$ymin <- -0.75; gnath$ymax <- -0.15; gnath$category <- "Gnathostome-conserved"
  acti$ymin  <- -1.55; acti$ymax  <- -0.95; acti$category  <- "Actinopterygii-specific"
  cne <- rbind(
    data.frame(xmin = gnath$xmin, xmax = gnath$xmax,
               ymin = gnath$ymin, ymax = gnath$ymax,
               category = gnath$category),
    data.frame(xmin = acti$xmin,  xmax = acti$xmax,
               ymin = acti$ymin,  ymax = acti$ymax,
               category = acti$category)
  )

  ## --- Gene models in the window (from the GTF) ---
  genes <- gtf_genes_raw[chrom_noversion == target_chrom & end >= ws & start <= we]
  genes <- as.data.frame(genes)
  genes$gstart  <- genes$start
  genes$gend    <- genes$end
  genes$strand  <- ifelse(genes$strand %in% c("+", "-"), genes$strand, "+")
  genes$label   <- ifelse(!is.na(genes$gene_name) & genes$gene_name != "",
                          genes$gene_name, genes$gene_id)
  genes$is_target <- tolower(genes$label) == tolower(r$cluster) |
                       tolower(genes$gene_id) == tolower(r$cluster)
  genes$category  <- ifelse(genes$is_target, "target gene", "other gene")
  genes$gstart_d <- pmax(genes$gstart, ws)
  genes$gend_d   <- pmin(genes$gend,   we)

  genes       <- genes[order(genes$gstart, genes$gend), , drop = FALSE]
  genes$row   <- seq_len(nrow(genes))
  genes$row_y <- genes$row * row_step
  n_rows <- if (nrow(genes) > 0) max(genes$row) else 0L

  gene_poly <- if (nrow(genes) > 0) {
    gene_arrow_polys(genes, head_bp = span * 0.016,
                     body_hh = body_hh, head_hh = head_hh)
  } else {
    data.frame(grp = integer(0), x = numeric(0), y = numeric(0),
               category = character(0))
  }

  lab <- genes
  lab$lab_x <- (lab$gstart_d + lab$gend_d) / 2
  lab$lab_y <- lab$row_y + 0.58

  gene_center  <- if (n_rows > 0) (row_step + n_rows * row_step) / 2 else row_step / 2
  gnath_center <- (gnath$ymin[1] + gnath$ymax[1]) / 2
  acti_center  <- (acti$ymin[1]  + acti$ymax[1])  / 2
  if (length(gnath_center) == 0 || is.na(gnath_center)) gnath_center <- -0.45
  if (length(acti_center)  == 0 || is.na(acti_center))  acti_center  <- -1.25

  p <- ggplot() +
    geom_hline(yintercept = 0.30, colour = "grey85", linewidth = 0.3) +
    geom_rect(data = cne,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                  fill = category),
              colour = NA, linewidth = 0.18) +
    geom_polygon(data = gene_poly,
                 aes(x = x, y = y, group = grp, fill = category),
                 colour = "grey25", linewidth = 0.28) +
    geom_text_repel(data = lab,
                    aes(x = lab_x, y = lab_y, label = label),
                    size = 2.6, fontface = "italic", colour = "grey15",
                    segment.size = 0.2, segment.colour = "grey70",
                    min.segment.length = 0.3, box.padding = 0.18,
                    max.overlaps = Inf, seed = 42) +
    scale_fill_manual(
      name   = NULL,
      values = c("target gene" = target_fill,
                 "other gene"  = other_fill,
                 "Gnathostome-conserved"    = gnath_fill,
                 "Actinopterygii-specific"  = acti_fill),
      limits = c("target gene", "other gene",
                 "Gnathostome-conserved", "Actinopterygii-specific"),
      drop   = FALSE) +
    scale_x_continuous(
      name   = paste0("Position on ", target_chrom, " (Mb)"),
      breaks = scales::pretty_breaks(n = 6),
      labels = function(b) sprintf("%.2f", b / 1e6),
      expand = expansion(mult = c(0.012, 0.012))) +
    scale_y_continuous(
      breaks = c(acti_center, gnath_center, gene_center),
      labels = c("Acti-spec", "Gnath-cons", "genes"),
      expand = expansion(mult = c(0.04, 0.04))) +
    coord_cartesian(xlim = c(ws, we),
                    ylim = c(-1.95, max(n_rows, 1) * row_step + 1.05),
                    clip = "off") +
    labs(title = paste0(r$cluster, "  \u00b7  ", target_chrom, ": ",
                        sprintf("%.2f", ws / 1e6), "\u2013",
                        sprintf("%.2f", we / 1e6), " Mb"),
         subtitle = paste0(nrow(gnath), " Gnath-cons CNEs  \u00b7  ",
                           nrow(acti), " Acti-spec CNEs  \u00b7  ",
                           nrow(genes), " gene models")) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
          axis.text.y   = element_text(face = "bold", colour = "grey30"),
          axis.ticks.y  = element_blank(),
          axis.title.y  = element_blank(),
          axis.line.x   = element_line(colour = "grey60", linewidth = 0.3),
          plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, colour = "grey40"),
          plot.title.position = "plot",
          legend.position = "bottom")

  list(plot = p, n_rows = n_rows,
       n_gnath = nrow(gnath), n_acti = nrow(acti), n_gene = nrow(genes))
}

## ---- 7. Build all panels and assemble the figure ---------------------------
panels <- lapply(seq_len(nrow(regions)), function(i) build_panel(regions[i, ]))

for (i in seq_along(panels)) {
  pi <- panels[[i]]
  message(sprintf("%s: %d Gnath-cons CNEs, %d Acti-spec CNEs, %d gene models, %d gene rows",
                  regions$cluster[i], pi$n_gnath, pi$n_acti, pi$n_gene, pi$n_rows))
}

h <- vapply(panels, function(p) max(p$n_rows, 1) + 5, numeric(1))

plot_list <- lapply(panels, function(p) p$plot)
figure <- wrap_plots(plot_list, ncol = 1, heights = h) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Actinopterygii-specific and Gnathostome-conserved CNEs around selected hotspot genes",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 14))) &
  theme(legend.position = "bottom")

## ---- 8. Save ----------------------------------------------------------------
total_rows <- sum(vapply(panels, function(p) max(p$n_rows, 1), numeric(1)))
fig_w <- 10
fig_h <- 3.0 * length(panels) + total_rows * 0.24

ggsave(paste0(out_prefix, ".pdf"), figure,
       width = fig_w, height = fig_h, limitsize = FALSE)
ggsave(paste0(out_prefix, ".png"), figure,
       width = fig_w, height = fig_h, dpi = 300, bg = "white",
       limitsize = FALSE)

message("Saved: ", out_prefix, ".pdf and ", out_prefix, ".png")
