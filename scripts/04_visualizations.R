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

preproc_dir <- "../output/preprocessed"

### --- Shared palette -------------------------------------------------------
# Okabe-Ito, matching Fig. 1A: reddish purple for the actinopterygian-specific
# set, bluish green for the gnathostome-conserved set, yellow for the
# sarcopterygian-specific set, grey for ATAC-seq peaks.
teleost_line <- "#CC79A7"
verte_line <- "#009E73"
sarco_line <- "#F0E442"
atac_line <- "gray50"
teleost_fill <- teleost_line
verte_fill <- verte_line
atac_fill <- atac_line

### --- Load preprocessed objects --------------------------------------------
gene_activity <- readRDS(file.path(preproc_dir, "gene_activity.rds"))
peak_anno_list_pre <- readRDS(file.path(preproc_dir, "peak_anno_list.rds"))

# Assemble the GRanges set the rest of the script expects, from explicit paths.
# Globbing was fragile — if any other script writes into granges_rds/ the
# naming clashed.
gr <- list(
  actinopteriigy_cne_gr = readRDS(
    file.path(preproc_dir, "actinopteriigy_cne_gr.rds")
  ),
  gnathostomata_cne_gr = readRDS(
    file.path(preproc_dir, "gnathostomata_cne_gr.rds")
  ),
  atac_peaks_gr = readRDS(file.path(preproc_dir, "atac_peaks_gr.rds"))
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
    which(gr$actinopteriigy_cne_gr$is_active)
  ],
  `Gnathostomata CNE with\nactive gene nearby` = gr$gnathostomata_cne_gr[
    which(gr$gnathostomata_cne_gr$is_active)
  ],
  `Actinopteriigy CNE with\nATACseq peak overlap and\nactive gene nearby` = subsetByOverlaps(
    gr$actinopteriigy_cne_gr[which(gr$actinopteriigy_cne_gr$is_active)],
    gr$atac_peaks_gr
  ),
  `Gnathostomata CNE with\nATACseq peak overlap and\nactive gene nearby` = subsetByOverlaps(
    gr$gnathostomata_cne_gr[which(gr$gnathostomata_cne_gr$is_active)],
    gr$atac_peaks_gr
  )
)

###############################################################################
### Annotation plot
###############################################################################
# ChIPseeker default genomicAnnotationPriority. The distance-to-TSS panel is
# now produced by rGREAT in script 2 from the regulatory-domain associations, so
# only the annotation composition is drawn here.

temp_list <- lapply(
  overlap_list,
  annotatePeak,
  overlap = "all",
  TxDb = drer_anno,
  tssRegion = c(-3000, 3000)
)

p1 <- plotAnnoBar(temp_list, title = NULL)
ggsave('../output/cne_annotations.pdf', p1, width = 8, height = 6)

###############################################################################
### Circos plot — zebrafish coordinates
###############################################################################
genome <- read_tsv("../ancilliary_files/sequence_report.tsv")

overlap_list <- lapply(overlap_list, relabel_seqlevels, genome)
gr <- lapply(gr, relabel_seqlevels, genome)

bin_size <- 2e6
autosomes <- paste0("chr", 1:25)

genome_genes <- gr$genome_genes_gr
genome_genes <- genome_genes[seqnames(genome_genes) %in% autosomes]
genome_genes <- keepSeqlevels(
  genome_genes,
  autosomes,
  pruning.mode = "coarse"
)

drer_chrom_sizes <- setNames(
  genome$`Seq length`[match(autosomes, genome$`Sequence name`)],
  autosomes
)

plot_cne_circos(
  tracks = list(
    overlap_list$`Actinopteriigy specific CNE`,
    overlap_list$`Gnathostomata specific CNE`,
    overlap_list$`Actinopteriigy CNE with\nATACseq peak overlap`
  ),
  chrom_sizes = drer_chrom_sizes,
  filepath = "../output/circos_element_density.pdf",
  track_colours = list(teleost_line, verte_line, atac_line),
  legend_labels = c(
    "Actinopterygii-specific",
    "Gnathostomata-conserved",
    "ATAC-Seq peaks (DANIO-CODE)"
  ),
  bin_size = bin_size,
  label_gap = 1,
  legend_title = "GRCz11"
)

###############################################################################
### Circos plot — human coordinates (sarcopterygii-specific set)
###############################################################################
# Symmetric counterpart of the zebrafish panel. The sarcopterygian-specific set
# and the gnathostome-conserved set called on the human reference share one
# coordinate space, so they are plotted together; there is no ATAC track because
# the accessibility data are zebrafish.

hsap_cne_path <- file.path(preproc_dir, "hsap_cne_gr.rds")
hsap_report <- "../ancilliary_files/sequence_report_hsap.tsv"


hsap_cne_gr <- readRDS(hsap_cne_path)
genome_hsap <- read_tsv(hsap_report, show_col_types = FALSE)

hsap_cne_gr <- lapply(hsap_cne_gr, function(gr) {
  gr <- keepSeqlevels(
    gr,
    genome_hsap$`RefSeq seq accession`,
    pruning.mode = 'coarse'
  )
  seqlevels(gr) <- genome_hsap$`UCSC style name`[
    match(seqlevels(gr), genome_hsap$`RefSeq seq accession`)
  ]
  return(gr)
})

hsap_autosomes <- paste0("chr", 1:22)
hsap_chrom_sizes <- setNames(
  genome_hsap$`Seq length`[
    match(hsap_autosomes, genome_hsap$`UCSC style name`)
  ],
  hsap_autosomes
)

# Same bin width as the zebrafish panel, so counts on the two circos plots are
# in identical units (CNEs per 2 Mb) and can be compared directly.
plot_cne_circos(
  tracks = list(
    hsap_cne_gr$sarcopterygii,
    hsap_cne_gr$gnathostomata
  ),
  chrom_sizes = hsap_chrom_sizes,
  filepath = "../output/circos_element_density_hsap.pdf",
  track_colours = list(sarco_line, verte_line),
  legend_labels = c("Sarcopterygii-specific", "Gnathostomata-conserved"),
  bin_size = bin_size,
  track_height = 0.14,
  label_gap = 0.8,
  legend_title = "GRCh38"
)

###########################################################################
### Width distributions — human sets
###########################################################################
hsap_cne_auto <- lapply(hsap_cne_gr, function(x) {
  x[as.character(seqnames(x)) %in% hsap_autosomes]
})

pdf("../output/power_law_like_fits_hsap.pdf", width = 11, height = 7)
par(mfrow = c(1, 2))
plot_cne_width_powerlaw(
  hsap_cne_auto$sarcopterygii,
  main = "Power-law like distribution of\nSarcopterygii specific CNE widths"
)
plot_cne_width_powerlaw(
  hsap_cne_auto$gnathostomata,
  main = "Power-law like distribution of\nGnathostomata conserved CNE widths\n(human reference)"
)
par(mfrow = c(1, 1))
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
  "chr12" , 44e6     , 45e6     , "../output/chr12_gviz_zoom.pdf" ,
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
