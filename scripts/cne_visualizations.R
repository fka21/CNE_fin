### LIBRARIES ###
# Install and load necessary libraries
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

# Load all required libraries using pacman
pacman::p_load(
  GenomicRanges,
  tidyverse,
  rtracklayer,
  gprofiler2,
  rGREAT,
  regioneR,
  patchwork,
  RColorBrewer,
  Biostrings,
  ChIPseeker,
  GenomicFeatures,
  UpSetR,
  ComplexHeatmap,
  GenometriCorr,
  circlize,
  GenomeInfoDb
)
### SET WORKING DIRECTORY ###
setwd(this.path::here())

source('custom_functions.R')

### READ IN DATA ###

gene_activity <- readRDS('../output/gene_activity.RDS')

files <- list.files("../output/granges_rds/", full.names = TRUE)
files <- files[str_detect(files, 'actino|gnath')]

gr <- setNames(
  lapply(files, readRDS),
  stringr::str_remove(basename(files), "\\.rds$")
)

######################
### CNE annotation ###
######################

drer_anno <- txdbmaker::makeTxDbFromGFF(
  "../ancilliary_files/drer.gff",
  organism = "Danio rerio"
)

# Annotate peaks and create visualizations
peak_anno_list <- lapply(
  list(
    gr$actinopteriigy_cne_gr,
    gr$gnathostomata_cne_gr,
    gr$actinopteriigy_cne_overlapping_atac_gr,
    gr$gnathostomata_cne_overlapping_atac_gr,
    gr$actinopteriigy_cne_active_with_atac_gr
  ),
  annotatePeak,
  TxDb = drer_anno,
  tssRegion = c(-3000, 3000),
  genomicAnnotationPriority = c(
    "Intergenic",
    "Downstream",
    "Promoter",
    "5UTR",
    "3UTR",
    "Intron",
    "Exon"
  )
)

names(peak_anno_list) <- c(
  "Actinopteriigy specific CNE",
  "Gnathostomata specific CNE",
  "Actinopteriigy CNE with\nATACseq peak overlap",
  "Gnathostomata CNE with\nATACseq peak overlap",
  "Actinopteriigy CNE with\nATACseq peak overlap and\nactive gene nearby"
)

p1 <- plotAnnoBar(peak_anno_list) + ggtitle(NULL)
p2 <- plotDistToTSS(peak_anno_list) + ggtitle(NULL) + ylab("CNE (%) (5' -> 3')")

p1 / p2
ggsave(
  '../output/subset_cne_chipseekr_annotated.pdf',
  height = 9,
  width = 9,
  units = 'in'
)

#############################
### Preparing circos plot ###
#############################

genome <- read_tsv('../ancilliary_files/sequence_report.tsv')
gr$genome_genes_gr <- genes(drer_anno)

for (nm in names(gr)) {
  print(nm)
  x <- gr[[nm]]
  # Map to new seqlevels first (drops NAs)
  new_levels <- genome$`Sequence name`[match(
    seqlevels(x),
    genome$`RefSeq seq accession`
  )]
  seqlevels(x) <- new_levels
  # Now set lengths (match to new levels)
  seqlengths(x) <- genome$`Seq length`[match(
    seqlevels(x),
    genome$`Sequence name`
  )]
  gr[[nm]] <- x
}

bin_size <- 2e6


# Keep only autosomes 1–25 for circos (optionally drop MT)
genes <- gr$genome_genes_gr
autosomes <- paste0("chr", 1:25)
genes <- genes[seqnames(genes) %in% autosomes]
genes <- keepSeqlevels(genes, autosomes, pruning.mode = "coarse")

# Make 2 Mb tiles across the genome
bins <- tileGenome(
  seqlengths(genes)[autosomes],
  tilewidth = bin_size,
  cut.last.tile.in.chrom = TRUE
)

# Count how many genes fall into each bin
gene_counts <- countOverlaps(bins, genes)

gene_density <- cbind(
  as.data.frame(bins)[, c("seqnames", "start", "end")],
  density = as.numeric(gene_counts)
) |>
  dplyr::rename(chromosome_name = seqnames)

# Ensure chromosome factor order is chr1..chr25
gene_density <- gene_density |>
  mutate(
    chromosome_name = factor(
      chromosome_name,
      levels = autosomes,
      ordered = TRUE
    )
  )


# example: bin ATAC peaks and one CNE set
atac_bins <- bin_granges(
  gr$atac_peaks_overlapping_actinopteriigy_cne_gr,
  bins,
  autosomes
)
teleost_bins <- bin_granges(gr$actinopteriigy_cne_gr, bins, autosomes)
vertebrate_bins <- bin_granges(gr$gnathostomata_cne_gr, bins, autosomes)

# xlim: one chrom per row, from 1 to seqlength
chromosome_levels <- autosomes
chrom_lengths <- as.numeric(seqlengths(genes)[chromosome_levels])

xlim <- cbind(
  start = rep(1, length(chromosome_levels)),
  end = chrom_lengths
)

# ---- colors ----
teleost_line <- "#1F78B4"
verte_line <- "skyblue"
atac_line <- "#E31A1C"

teleost_fill <- "#1F78B4"
verte_fill <- "skyblue"
atac_fill <- "#E31A1C"

pdf("../output/circos_element_density.pdf", width = 9, height = 9)

circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = c(rep(1, length(chromosome_levels) - 1), 4),
  points.overflow.warning = FALSE
)

circos.initialize(
  factors = chromosome_levels,
  xlim = xlim
)

## Gene density track (analogous to the TE script)
circos.trackPlotRegion(
  factors = teleost_bins$`Chromosome name`,
  x = teleost_bins$start,
  y = teleost_bins$hit_count,
  track.height = 0.12,
  panel.fun = function(x, y) {
    # line + area for density
    circos.lines(x, y, col = teleost_line, area = T, lwd = 1, type = "s")
    circos.segments(
      x0 = x,
      y0 = 0,
      x1 = x,
      y1 = y,
      col = adjustcolor("grey70", alpha.f = 0.3)
    )

    # x‑axis only once per sector, with Mb labels

    circos.xaxis(
      labels.facing = "clockwise",
      labels.niceFacing = TRUE,
      major.at = c(0, 20e6, 40e6, 60e6),
      labels = c("0 Mb", "20 Mb", "40 Mb", "60 Mb"),
      labels.cex = 0.5
    )
  }
)

circos.trackPlotRegion(
  factors = vertebrate_bins$`Chromosome name`,
  x = vertebrate_bins$start,
  y = vertebrate_bins$hit_count,
  track.height = 0.12,
  panel.fun = function(x, y) {
    # line + area for density
    circos.lines(x, y, col = verte_line, area = T, lwd = 1, type = "s")
    circos.segments(
      x0 = x,
      y0 = 0,
      x1 = x,
      y1 = y,
      col = adjustcolor("grey70", alpha.f = 0.3)
    )
  }
)

circos.trackPlotRegion(
  factors = atac_bins$`Chromosome name`,
  x = atac_bins$start,
  y = atac_bins$hit_count,
  track.height = 0.12,
  panel.fun = function(x, y) {
    # line + area for density
    circos.lines(x, y, col = atac_line, area = T, lwd = 1, type = "s")
    circos.segments(
      x0 = x,
      y0 = 0,
      x1 = x,
      y1 = y,
      col = adjustcolor("grey70", alpha.f = 0.3)
    )
  }
)


for (chr in chromosome_levels) {
  # explicitly coerce the scalar to numeric, but preserve name lookup
  mid <- as.numeric(seqlengths(gr$genome_genes_gr)[chr]) / 2

  dens_chr <- gene_density$density[gene_density$chromosome_name == chr]
  if (!length(dens_chr)) {
    next
  } # skip chromosomes without density bins

  circos.text(
    sector.index = chr,
    track.index = 1, # or the track where your gene density is
    x = mid,
    y = -310,
    labels = chr,
    facing = "bending.inside",
    cex = 0.6
  )
}


circos.yaxis(sector.index = "chr1", track.index = 1, labels.cex = 0.5)
circos.yaxis(sector.index = "chr1", track.index = 2, labels.cex = 0.4)
circos.yaxis(sector.index = "chr1", track.index = 3, labels.cex = 0.35)

# ---- legend overlay in the center ----
par(fig = c(0, 1, 0, 1), new = TRUE)
plot.new()

legend(
  "center",
  legend = c(
    "Actinopteriigy",
    "Gnathostomata",
    "ATAC-seq peaks"
  ),
  fill = c(teleost_fill, verte_fill, atac_fill),
  border = c(teleost_line, verte_line, atac_line),
  bty = "n",
  cex = 0.9
)

circos.clear()
dev.off()

#################################
### CNEr based visualizations ###
#################################

library(CNEr)

gr_plot <- GRangePairs(
  first = gr$actinopteriigy_cne_gr,
  second = gr$actinopteriigy_cne_gr
)


gr$actinopteriigy_cne_gr <- gr$actinopteriigy_cne_gr[
  seqnames(gr$actinopteriigy_cne_gr) %in% autosomes
]
genes <- keepSeqlevels(genes, autosomes, pruning.mode = "coarse")


par(mfrow = c(1, 2))
pdf('../output/power_law_like_fits.pdf', width = 9, height = 7)
plot_cne_width_powerlaw(
  gr$actinopteriigy_cne_gr,
  main = 'Power-law like distribution of\nTeleost specific CNE widths'
)
plot_cne_width_powerlaw(
  gr$gnathostomata_cne_gr,
  main = 'Power-law like distribution of\nCNEs shared with Vertebrata widths'
)
dev.off()
par(mfrow = c(1, 1))

p1 <- plotCNEDistribution(gr$actinopteriigy_cne_gr, chrs = c('chr5', 'chr14')) +
  theme_bw(base_size = 15) +
  labs(title = 'Teleost specific CNEs')

p2 <- plotCNEDistribution(gr$gnathostomata_cne_gr, chrs = c('chr5', 'chr14')) +
  theme_bw(base_size = 15) +
  labs(title = 'CNEs shared with Vertebrata')

p1 / p2
ggsave('../output/cne_clustering.pdf', width = 9, height = 8, units = 'in')

##################
### Gviz plots ###
##################

library(Gviz)

plot_gviz_zoom(
  gr = gr,
  chr = "chr5",
  start = 32.5e6,
  end = 33.5e6,
  gene_activity = gene_activity,
  active_col = "royalblue",
  inactive_col = "skyblue",
  filepath = "../output/chr5_gviz_zoom.pdf",
  height = 10,
  teleost_fill = teleost_line,
  verte_fill = verte_line,
  atac_fill = atac_line
)

plot_gviz_zoom(
  gr = gr,
  chr = "chr14",
  start = 2e6,
  end = 3e6,
  filepath = '../output/chr14_gviz_zoom.pdf',
  height = 10,
  gene_activity = gene_activity,
  active_col = "royalblue",
  inactive_col = "skyblue",
  teleost_fill = teleost_line,
  verte_fill = verte_line,
  atac_fill = atac_line
)
