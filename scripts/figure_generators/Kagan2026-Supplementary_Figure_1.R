#!/usr/bin/env Rscript
# =============================================================================
# Visualization of Actinopterygii-specific CNEs and pcdh gene clusters
# -----------------------------------------------------------------------------
# Produces a two-panel figure:
#   Panel A : pcdh1 cluster  - chr10 (NC_133185), 26.43-26.68 Mb
#   Panel B : pcdh2 cluster  - chr14 (NC_133189),  2.35-2.85 Mb
#
# Each panel shows two tracks sharing an x-axis:
#   * "pcdh genes" : gene models drawn as direction-aware arrows, packed into
#                    non-overlapping rows. pcdh genes are highlighted and
#                    labelled; other genes are shown in grey for context.
#   * "CNEs"       : Actinopterygii-specific CNEs, drawn as boxes (#CC79A7).
#
# INPUT FILES (edit the paths below, or run the script from the folder that
# contains them):
#   actinopterygii_specific_drer.bed   -> CNE coordinates (BED6)
#   actinopteriigy_cne_final_table.tsv -> annotation table, used for gene models
#
# OUTPUT: actinopterygii_cne_pcdh_panels.pdf  and  .png
# =============================================================================

## ---- 0. Packages -----------------------------------------------------------
required <- c("ggplot2", "dplyr", "patchwork", "scales", "ggrepel")
missing  <- required[!vapply(required, requireNamespace,
                             FUN.VALUE = logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
  library(scales);  library(ggrepel)
})

## ---- 1. Configuration ------------------------------------------------------
bed_file   <- "gnathostomata_conserved_drer.bed"      # CNE coordinates
tsv_file   <- "gnathostomata_cne_final_table.tsv"    # annotation / gene models
out_prefix <- "gnathostomata_cne_panels_zfhx4"        # output file prefix

# Colours
# cne_fill    <- "#CC79A7"   # Actinopterygii-specific CNEs (as requested)
# cne_outline <- "#9E5A82"   # darker shade of the CNE colour, for box borders

cne_fill    <- "#009E73"   # Gnathostome-conserved CNEs (as requested)
cne_outline <- "#00755A"   # darker shade of the CNE colour, for box borders
pcdh_fill   <- "#2C6E9B"   # pcdh gene models
other_fill  <- "#BFBFBF"   # non-pcdh gene models (context)

# Regions to plot. The chromosome accessions follow the zebrafish RefSeq
# assembly in which chr1 = NC_133176 ... chr25 = NC_133200, hence
# chr10 = NC_133185 and chr14 = NC_133189. Adjust here if your assembly differs.
#regions <- data.frame(
#  cluster   = c("pcdh1 cluster",  "pcdh2 cluster"),
#  chr_label = c("chr10",          "chr14"),
#  accession = c("NC_133185",      "NC_133199"),
#  win_start = c(26.43e6,          2.35e6),
#  win_end   = c(26.68e6,          2.85e6),
#  stringsAsFactors = FALSE
#)

regions <- data.frame(
  cluster   = c("ebf3a",  "zfhx4"),
  chr_label = c("chr12",          "chr24"),
  accession = c("NC_133187",      "NC_133199"),
  win_start = c(44.44e6,          27.1e6),
  win_end   = c(44.7e6,          27.3e6),
  stringsAsFactors = FALSE
)

## ---- 2. Read data ----------------------------------------------------------
if (!file.exists(bed_file)) stop("BED file not found: ", bed_file,
                                 "\n  Set 'bed_file' or setwd() to its folder.")
if (!file.exists(tsv_file)) stop("TSV file not found: ", tsv_file,
                                 "\n  Set 'tsv_file' or setwd() to its folder.")

# CNEs from the BED file (BED6: chrom, start, end, name, score, strand)
cne_all <- read.table(bed_file, sep = "\t", header = FALSE, quote = "",
                       fill = TRUE, stringsAsFactors = FALSE,
                       col.names = c("chrom", "start", "end",
                                     "name", "score", "strand"))
# Drop the RefSeq version suffix (".1") so BED and TSV chrom names match
cne_all$chrom <- sub("\\.[0-9]+$", "", cne_all$chrom)

# Annotation table -> gene models
tab <- read.delim(tsv_file, header = TRUE, sep = "\t", quote = "",
                   comment.char = "", check.names = FALSE,
                   stringsAsFactors = FALSE)
tab$chrom <- sub("\\.[0-9]+$", "", tab$seqnames)

## ---- 3. Helper functions ---------------------------------------------------

# Numeric strand code (1/2) -> "+"/"-"
strand_sym <- function(s) ifelse(s == 1, "+", ifelse(s == 2, "-", NA_character_))

# Greedy interval packing: assign each [start,end] interval to the lowest row
# on which it does not overlap (within 'pad' bp) anything already placed.
pack_intervals <- function(starts, ends, pad = 0) {
  ord      <- order(starts, ends)
  row_last <- numeric(0)            # last occupied end-coordinate per row
  rows     <- integer(length(starts))
  for (i in ord) {
    placed <- FALSE
    for (r in seq_along(row_last)) {
      if (starts[i] > row_last[r] + pad) {
        rows[i] <- r; row_last[r] <- ends[i]; placed <- TRUE; break
      }
    }
    if (!placed) { row_last <- c(row_last, ends[i]); rows[i] <- length(row_last) }
  }
  rows
}

# Build geom_polygon vertices for direction-aware gene arrows.
# Returns a long data frame (one block of 7 vertices per gene).
gene_arrow_polys <- function(g, head_bp, body_hh, head_hh) {
  out <- vector("list", nrow(g))
  for (i in seq_len(nrow(g))) {
    x1 <- g$gstart_d[i]; x2 <- g$gend_d[i]; y <- g$row_y[i]
    hw <- min(head_bp, x2 - x1)             # arrowhead width (capped for short genes)
    if (g$strand[i] == "+") {               # arrow points right (5'->3')
      bx <- x2 - hw
      xs <- c(x1, bx, bx, x2, bx, bx, x1)
    } else {                                 # arrow points left
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

## ---- 4. Build one panel ----------------------------------------------------
build_panel <- function(r) {
  acc <- r$accession; ws <- r$win_start; we <- r$win_end
  span <- we - ws
  row_step <- 1.15
  body_hh  <- 0.22; head_hh <- 0.40           # arrow body / head half-heights

  ## --- CNEs in the window (from the BED file) ---
  cne <- cne_all[cne_all$chrom == acc &
                 cne_all$end >= ws & cne_all$start <= we, , drop = FALSE]
  # Enforce a minimum drawn width so very short CNEs stay visible
  minw <- span * 0.0022
  mid  <- (cne$start + cne$end) / 2
  cne$xmin <- pmax(ws, pmin(cne$start, mid - minw / 2))
  cne$xmax <- pmin(we, pmax(cne$end,   mid + minw / 2))
  cne_ymin <- -1.05; cne_ymax <- -0.40
  cne$ymin <- cne_ymin; cne$ymax <- cne_ymax
  cne$category <- "Gnathostome-conserved CNE"

  ## --- Gene models (from the TSV) ---
  # Each protocadherin isoform is annotated as a separate gene that extends to
  # the shared constant exons, so collapse to one model per geneId.
  genes <- tab %>%
    filter(chrom == acc, !is.na(geneId), geneId != "", geneId != "NA") %>%
    group_by(geneId) %>%
    summarise(gstart = min(geneStart, na.rm = TRUE),
              gend   = max(geneEnd,   na.rm = TRUE),
              strand = strand_sym(geneStrand[1]),
              .groups = "drop") %>%
    filter(gend >= ws, gstart <= we) %>%          # keep genes overlapping window
    as.data.frame()

  genes$label    <- sub("^gene-", "", genes$geneId)
  genes$is_pcdh  <- grepl("pcdh", genes$geneId, ignore.case = TRUE)
  genes$category <- ifelse(genes$is_pcdh, "pcdh gene", "other gene")
  # Clip to the window so arrowheads never spill into the margin
  genes$gstart_d <- pmax(genes$gstart, ws)
  genes$gend_d   <- pmin(genes$gend,   we)

  # Pack genes into non-overlapping rows (the nested-cluster "staircase")
  genes$row   <- pack_intervals(genes$gstart, genes$gend, pad = span * 0.006)
  genes$row_y <- genes$row * row_step
  n_rows <- max(genes$row)

  gene_poly <- gene_arrow_polys(genes, head_bp = span * 0.016,
                                body_hh = body_hh, head_hh = head_hh)

  # Labels: pcdh genes only (other genes shown for context but left unlabelled)
  # lab <- genes[genes$is_pcdh, , drop = FALSE]
  lab <- genes
  lab$lab_x <- (lab$gstart_d + lab$gend_d) / 2
  lab$lab_y <- lab$row_y + 0.58

  ## --- Track centres for the y-axis labels ---
  gene_center <- (row_step + n_rows * row_step) / 2
  cne_center  <- (cne_ymin + cne_ymax) / 2

  ## --- Assemble the panel ---
  p <- ggplot() +
    # faint separator between the two tracks
    geom_hline(yintercept = -0.05, colour = "grey85", linewidth = 0.3) +
    # CNE feature track
    geom_rect(data = cne,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                  fill = category),
              colour = cne_outline, linewidth = 0.18) +
    # gene-model arrows
    geom_polygon(data = gene_poly,
                 aes(x = x, y = y, group = grp, fill = category),
                 colour = "grey25", linewidth = 0.28) +
    # gene labels (pcdh only)
    geom_text_repel(data = lab,
                    aes(x = lab_x, y = lab_y, label = label),
                    size = 2.6, fontface = "italic", colour = "grey15",
                    segment.size = 0.2, segment.colour = "grey70",
                    min.segment.length = 0.3, box.padding = 0.18,
                    max.overlaps = Inf, seed = 42) +
    scale_fill_manual(
      name   = NULL,
      values = c("pcdh gene" = pcdh_fill,
                 "other gene" = other_fill,
                 "Gnathostome-conserved CNE" = cne_fill),
      limits = c("pcdh gene", "other gene", "Gnathostome-conserved CNE"),
      drop   = FALSE) +
    scale_x_continuous(
      name   = paste0("Position on ", r$chr_label, " (Mb)"),
      breaks = scales::pretty_breaks(n = 6),
      labels = function(b) sprintf("%.2f", b / 1e6),
      expand = expansion(mult = c(0.012, 0.012))) +
    scale_y_continuous(
      breaks = c(cne_center, gene_center),
      labels = c("CNEs", "pcdh genes"),
      expand = expansion(mult = c(0.04, 0.04))) +
    coord_cartesian(xlim = c(ws, we),
                    ylim = c(-1.45, n_rows * row_step + 1.05),
                    clip = "off") +
    labs(title = paste0(r$cluster, "  \u00b7  ", r$chr_label, ": ",
                        sprintf("%.2f", ws / 1e6), "\u2013",
                        sprintf("%.2f", we / 1e6), " Mb"),
         subtitle = paste0(nrow(cne), " Actinopterygii-specific CNEs  \u00b7  ",
                           nrow(genes), " gene models (",
                           sum(genes$is_pcdh), " pcdh)")) +
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

  list(plot = p, n_rows = n_rows, n_cne = nrow(cne),
       n_gene = nrow(genes), n_pcdh = sum(genes$is_pcdh))
}

## ---- 5. Build both panels and assemble the figure --------------------------
panels <- lapply(seq_len(nrow(regions)), function(i) build_panel(regions[i, ]))

for (i in seq_along(panels)) {
  pi <- panels[[i]]
  message(sprintf("%s (%s): %d CNEs, %d genes (%d pcdh), %d gene rows",
                  regions$cluster[i], regions$chr_label[i],
                  pi$n_cne, pi$n_gene, pi$n_pcdh, pi$n_rows))
}

# Panel heights scale with the number of gene rows so neither panel is squashed
h <- vapply(panels, function(p) p$n_rows + 5, numeric(1))

figure <- (panels[[1]]$plot / panels[[2]]$plot) +
  plot_layout(heights = h, guides = "collect") +
  plot_annotation(
    title    = "Actinopterygii-specific CNEs across the zebrafish pcdh clusters",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 14))) &
  theme(legend.position = "bottom")

## ---- 6. Save ---------------------------------------------------------------
total_rows <- sum(vapply(panels, function(p) p$n_rows, numeric(1)))
fig_w <- 10
fig_h <- 4.5 + total_rows * 0.24

ggsave(paste0(out_prefix, ".pdf"), figure,
       width = fig_w, height = fig_h, limitsize = FALSE)
ggsave(paste0(out_prefix, ".png"), figure,
       width = fig_w, height = fig_h, dpi = 300, bg = "white",
       limitsize = FALSE)

message("Saved: ", out_prefix, ".pdf and ", out_prefix, ".png")
