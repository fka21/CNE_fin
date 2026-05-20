### LIBRARIES ###
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(
  GenomicRanges,
  tidyverse,
  rtracklayer,
  ChIPseeker,
  GenomicFeatures,
  patchwork,
  ComplexHeatmap,
  circlize,
  GenomeInfoDb,
  CNEr,
  Gviz
)

setwd(this.path::here())
source("custom_functions.R")

### --- Load preprocessed objects --------------------------------------------
gene_activity <- readRDS("../output/preprocessed/gene_activity.rds")
peak_anno_list_pre <- readRDS("../output/preprocessed/peak_anno_list.rds")

# Assemble the GRanges set the rest of the script expects, from explicit paths.
# Globbing was fragile — if any other script writes into granges_rds/ the
# naming clashed.
gr <- list(
  actinopteriigy_cne_gr = readRDS(
    "../output/preprocessed/actinopteriigy_cne_gr.rds"
  ),
  gnathostomata_cne_gr = readRDS(
    "../output/preprocessed/gnathostomata_cne_gr.rds"
  ),
  atac_peaks_gr = readRDS("../output/granges_rds/atac_peaks_gr.rds")
)

### --- TxDb (rebuilt; serialisation of TxDb objects is brittle) -------------
drer_anno <- txdbmaker::makeTxDbFromGFF(
  "../ancilliary_files/drer.gff",
  organism = "Danio rerio"
)
gr$genome_genes_gr <- genes(drer_anno)
gr$atac_peaks_overlapping_actinopteriigy_cne_gr <- subsetByOverlaps(
  gr$actinopteriigy_cne_gr,
  gr$atac_peaks_gr
)

### --- Build overlap list for all subsets
overlap_list <- list(
  `Actinopteriigy specific CNE` = gr$actinopteriigy_cne_gr,
  `Gnathostomata specific CNE` = gr$gnathostomata_cne_gr,
  `Actinopteriigy CNE with\nATACseq peak overlap` = subsetByOverlaps(
    gr$actinopteriigy_cne_gr,
    gr$atac_peaks_gr
  ),
  `Gnathostomata CNE with\nATACseq peak overlap` = subsetByOverlaps(
    gr$gnathostomata_cne_gr,
    gr$atac_peaks_gr
  ),
  `Actinopteriigy CNE with\nactive gene nearby` = gr$actinopteriigy_cne_gr[
    gr$actinopteriigy_cne_gr$is_active == T
  ],
  `Gnathostomata CNE with\nactive gene nearby` = gr$gnathostomata_cne_gr[
    gr$gnathostomata_cne_gr$is_active == T
  ],
  `Actinopteriigy CNE with\nATACseq peak overlap and\nactive gene nearby` = subsetByOverlaps(
    gr$actinopteriigy_cne_gr[gr$actinopteriigy_cne_gr$is_active == T],
    gr$atac_peaks_gr
  ),
  `Gnathostomata CNE with\nATACseq peak overlap and\nactive gene nearby` = subsetByOverlaps(
    gr$gnathostomata_cne_gr[gr$gnathostomata_cne_gr$is_active == T],
    gr$atac_peaks_gr
  )
)

###############################################################################
### Circos plot
###############################################################################
genome <- read_tsv("../ancilliary_files/sequence_report.tsv")

overlap_list <- lapply(overlap_list, relabel_seqlevels, genome)
gr <- lapply(gr, relabel_seqlevels, genome)

bin_size <- 2e6
autosomes <- paste0("chr", 1:25)

genes <- gr$genome_genes_gr
genes <- genes[seqnames(genes) %in% autosomes]
genes <- keepSeqlevels(genes, autosomes, pruning.mode = "coarse")

bins <- tileGenome(
  seqlengths(genes)[autosomes],
  tilewidth = bin_size,
  cut.last.tile.in.chrom = TRUE
)

gene_counts <- countOverlaps(bins, genes)
gene_density <- cbind(
  as.data.frame(bins)[, c("seqnames", "start", "end")],
  density = as.numeric(gene_counts)
) |>
  dplyr::rename(chromosome_name = seqnames) |>
  mutate(
    chromosome_name = factor(
      chromosome_name,
      levels = autosomes,
      ordered = TRUE
    )
  )

atac_bins <- bin_granges(
  overlap_list$`Actinopteriigy CNE with\nATACseq peak overlap`,
  bins,
  autosomes
)
teleost_bins <- bin_granges(
  overlap_list$`Actinopteriigy specific CNE`,
  bins,
  autosomes
)
vertebrate_bins <- bin_granges(
  overlap_list$`Gnathostomata specific CNE`,
  bins,
  autosomes
)

chrom_lengths <- as.numeric(seqlengths(genes)[autosomes])
xlim <- cbind(start = rep(1, length(autosomes)), end = chrom_lengths)

teleost_line <- "#CC79A7"
verte_line <- "#009E73"
atac_line <- "gray50"
teleost_fill <- teleost_line
verte_fill <- verte_line
atac_fill <- atac_line

pdf("../output/circos_element_density.pdf", width = 9, height = 9)
circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = c(rep(1, length(autosomes) - 1), 4),
  points.overflow.warning = FALSE
)
circos.initialize(factors = autosomes, xlim = xlim)

# track 1 — actinopteriigy
circos.trackPlotRegion(
  factors = teleost_bins$`Chromosome name`,
  x = teleost_bins$start,
  y = teleost_bins$hit_count,
  track.height = 0.12,
  panel.fun = function(x, y) {
    circos.lines(x, y, col = teleost_line, area = TRUE, lwd = 1, type = "s")
    circos.segments(
      x0 = x,
      y0 = 0,
      x1 = x,
      y1 = y,
      col = adjustcolor("grey70", alpha.f = 0.3)
    )
    circos.xaxis(
      labels.facing = "clockwise",
      labels.niceFacing = TRUE,
      major.at = c(0, 20e6, 40e6, 60e6),
      labels = c("0 Mb", "20 Mb", "40 Mb", "60 Mb"),
      labels.cex = 0.5
    )
  }
)

# track 2 — gnathostomata
circos.trackPlotRegion(
  factors = vertebrate_bins$`Chromosome name`,
  x = vertebrate_bins$start,
  y = vertebrate_bins$hit_count,
  track.height = 0.12,
  panel.fun = function(x, y) {
    circos.lines(x, y, col = verte_line, area = TRUE, lwd = 1, type = "s")
    circos.segments(
      x0 = x,
      y0 = 0,
      x1 = x,
      y1 = y,
      col = adjustcolor("grey70", alpha.f = 0.3)
    )
  }
)

# track 3 — ATAC
circos.trackPlotRegion(
  factors = atac_bins$`Chromosome name`,
  x = atac_bins$start,
  y = atac_bins$hit_count,
  track.height = 0.12,
  panel.fun = function(x, y) {
    circos.lines(x, y, col = atac_line, area = TRUE, lwd = 1, type = "s")
    circos.segments(
      x0 = x,
      y0 = 0,
      x1 = x,
      y1 = y,
      col = adjustcolor("grey70", alpha.f = 0.3)
    )
  }
)

for (chr in autosomes) {
  mid <- as.numeric(seqlengths(gr$genome_genes_gr)[chr]) / 2
  dens_chr <- gene_density$density[gene_density$chromosome_name == chr]
  if (!length(dens_chr)) {
    next
  }
  circos.text(
    sector.index = chr,
    track.index = 1,
    x = mid,
    y = -310,
    labels = chr,
    facing = "bending.inside",
    cex = 0.6
  )
}

circos.yaxis(sector.index = "chr1", track.index = 1, labels.cex = 0.50)
circos.yaxis(sector.index = "chr1", track.index = 2, labels.cex = 0.40)
circos.yaxis(sector.index = "chr1", track.index = 3, labels.cex = 0.35)

par(fig = c(0, 1, 0, 1), new = TRUE)
plot.new()
legend(
  "center",
  legend = c("Actinopteriigy", "Gnathostomata", "ATAC-seq peaks"),
  fill = c(teleost_fill, verte_fill, atac_fill),
  border = c(teleost_line, verte_line, atac_line),
  bty = "n",
  cex = 0.9
)
circos.clear()
dev.off()

###############################################################################
### CNEr-based visualisations
###############################################################################
gr$actinopteriigy_cne_gr <- gr$actinopteriigy_cne_gr[
  seqnames(gr$actinopteriigy_cne_gr) %in% autosomes
]
gr$gnathostomata_cne_gr <- gr$gnathostomata_cne_gr[
  seqnames(gr$gnathostomata_cne_gr) %in% autosomes
]

pdf("../output/power_law_like_fits.pdf", width = 11, height = 7)
par(mfrow = c(1, 2))
plot_cne_width_powerlaw(
  gr$actinopteriigy_cne_gr,
  main = "Power-law like distribution of\nActinopteriigy specific CNE widths"
)
plot_cne_width_powerlaw(
  gr$gnathostomata_cne_gr,
  main = "Power-law like distribution of\nGnathostomata specific CNE widths"
)
par(mfrow = c(1, 1))
dev.off()

p1 <- plotCNEDistribution(gr$actinopteriigy_cne_gr, chrs = c("chr5", "chr14")) +
  theme_bw(base_size = 15) +
  labs(title = "Actinopteriigy specific CNEs")
p2 <- plotCNEDistribution(gr$gnathostomata_cne_gr, chrs = c("chr5", "chr14")) +
  theme_bw(base_size = 15) +
  labs(title = "Gnathostomata specific CNEs")
p1 / p2
ggsave("../output/cne_clustering.pdf", width = 9, height = 8, units = "in")

###############################################################################
### Gviz — fixed zoom regions
###############################################################################
zooms <- tribble(
  ~chr    , ~start   , ~end     , ~file                           ,
  "chr5"  ,   32.5e6 ,   33.5e6 , "../output/chr5_gviz_zoom.pdf"  ,
  "chr14" ,  2e6     ,  3e6     , "../output/chr14_gviz_zoom.pdf" ,
  "chr18" ,   39.5e6 ,   40.5e6 , "../output/chr18_gviz_zoom.pdf" ,
  "chr9"  ,   39.5e6 ,   40.5e6 , "../output/chr9_gviz_zoom.pdf"  ,
  "chr19" , 39e6     , 41e6     , "../output/chr19_gviz_zoom.pdf"
)

for (i in seq_len(nrow(zooms))) {
  plot_gviz_zoom(
    gr = gr,
    chr = zooms$chr[i],
    start = zooms$start[i],
    end = zooms$end[i],
    filepath = zooms$file[i],
    height = 10,
    gene_activity = gene_activity,
    active_col = "royalblue",
    inactive_col = "skyblue",
    teleost_fill = teleost_line,
    verte_fill = verte_line,
    atac_fill = atac_line
  )
}

###############################################################################
### Gviz — GO:0007224 (smoothened signaling pathway) gene panels
###############################################################################
go_genes_gr <- readRDS("../output/great/gviz/GO0007224_genes.rds")
assoc_smo_actino <- readRDS("../output/great/gviz/GO0007224_cnes_actino.rds")
assoc_smo_gnatho <- readRDS("../output/great/gviz/GO0007224_cnes_gnatho.rds")

# Relabel the GO gene/CNE objects to the same chr* scheme as the rest of `gr`.
go_genes_gr <- relabel_seqlevels(go_genes_gr, genome)
assoc_smo_actino <- relabel_seqlevels(assoc_smo_actino, genome)
assoc_smo_gnatho <- relabel_seqlevels(assoc_smo_gnatho, genome)

# Filter to autosomes so they line up with the precomputed `gr` set.
go_genes_gr <- go_genes_gr[as.character(seqnames(go_genes_gr)) %in% autosomes]
assoc_smo_actino <- assoc_smo_actino[
  as.character(seqnames(assoc_smo_actino)) %in% autosomes
]
assoc_smo_gnatho <- assoc_smo_gnatho[
  as.character(seqnames(assoc_smo_gnatho)) %in% autosomes
]

dir.create("../output/gviz_GO0007224", showWarnings = FALSE, recursive = TRUE)

# Each plot centres on a gene body + 200 kb flank so the basal regulatory
# domain is visible. plot_gviz_zoom already knows how to render `gr`'s tracks.
flank_bp <- 200e3
for (i in seq_along(go_genes_gr)) {
  g <- go_genes_gr[i]
  chr <- as.character(seqnames(g))
  start <- max(1L, start(g) - flank_bp)
  end <- end(g) + flank_bp
  name <- mcols(g)$gene_id
  if (is.null(name) || is.na(name) || !nzchar(name)) {
    name <- paste0(chr, "_", start(g))
  }
  out <- file.path("../output/gviz_GO0007224", paste0(name, ".pdf"))

  tryCatch(
    plot_gviz_zoom(
      gr = gr,
      chr = chr,
      start = start,
      end = end,
      filepath = out,
      height = 10,
      gene_activity = gene_activity,
      active_col = "royalblue",
      inactive_col = "skyblue",
      teleost_fill = teleost_line,
      verte_fill = verte_line,
      atac_fill = atac_line
    ),
    error = function(e) {
      message("Gviz failed for ", name, ": ", conditionMessage(e))
    }
  )
}
