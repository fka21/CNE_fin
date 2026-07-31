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
    gr$atac_peaks_gr
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
  legend_title = "GRCz12"
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

dir.create("../output/reports", showWarnings = FALSE, recursive = TRUE)

# Gene models for the human locus panels. Both packages are optional: without
# them the CNE density and element tracks are still drawn, only the gene track
# is dropped. UCSC seqlevels match the relabelled CNE objects.
hsap_txdb <- if (
  requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)
) {
  TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
} else {
  message(
    "TxDb.Hsapiens.UCSC.hg38.knownGene not installed; ",
    "human panels will omit the gene track."
  )
  NULL
}

hsap_orgdb <- if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  org.Hs.eg.db::org.Hs.eg.db
} else {
  NULL
}

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

###########################################################################
### Gviz — Sarcopterygii-specific density peaks (human coordinates)      ###
###########################################################################
# The circos panel shows sarcopterygii-specific elements concentrating on
# chr1 and chr2 and towards the ends of chr10. Rather than hard-coding those
# coordinates, the densest 1 Mb window is located in each target interval and
# the plotting window is centred on it, so the panels track the element set if
# it is regenerated. Chosen windows are written out so they can be fixed for
# the final figure once the set is frozen.

hsap_gviz_dir <- "../output/gviz_hsap"
dir.create(hsap_gviz_dir, showWarnings = FALSE, recursive = TRUE)

chr10_len <- hsap_chrom_sizes[["chr10"]]
telomere_span <- 10e6

spike_targets <- list(
  list(label = "chr1", chrom = "chr1", region = NULL),
  list(label = "chr2", chrom = "chr2", region = NULL),
  list(
    label = "chr10_p_telomere",
    chrom = "chr10",
    region = c(1, telomere_span)
  ),
  list(
    label = "chr10_q_telomere",
    chrom = "chr10",
    region = c(chr10_len - telomere_span, chr10_len)
  )
)

spike_windows <- purrr::map_dfr(spike_targets, function(tg) {
  top <- find_cne_spikes(
    hsap_cne_gr$sarcopterygii,
    chrom = tg$chrom,
    bin_size = 1e6,
    region = tg$region,
    n = 3,
    chrom_length = hsap_chrom_sizes[[tg$chrom]]
  )
  if (!nrow(top)) {
    return(tibble())
  }
  top %>%
    mutate(target = tg$label, rank = row_number()) %>%
    relocate(target, rank)
})

write_tsv(
  spike_windows,
  "../output/reports/sarcopterygii_density_peaks_hsap.tsv"
)
print(spike_windows)

# Plot the top window per target, padded so the surrounding gene context is
# visible rather than only the peak itself.
plot_span <- 4e6

purrr::walk(unique(spike_windows$target), function(tg) {
  win <- spike_windows %>% filter(target == tg, rank == 1)
  if (!nrow(win)) {
    return(invisible(NULL))
  }

  centre <- (win$start + win$end) / 2
  chrom_len <- hsap_chrom_sizes[[win$chrom]]
  from <- max(1, round(centre - plot_span / 2))
  to <- min(chrom_len, round(centre + plot_span / 2))

  plot_hsap_cne_locus(
    cne_list = list(
      `Sarcopterygii-specific` = hsap_cne_gr$sarcopterygii,
      `Gnathostomata-conserved` = hsap_cne_gr$gnathostomata
    ),
    chrom = win$chrom,
    start = from,
    end = to,
    filepath = file.path(hsap_gviz_dir, paste0(tg, "_cne_density.pdf")),
    track_colours = list(
      `Sarcopterygii-specific` = sarco_line,
      `Gnathostomata-conserved` = verte_line
    ),
    txdb = hsap_txdb,
    orgdb = hsap_orgdb,
    density_bin = 50e3
  )
})


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
