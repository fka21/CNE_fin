### LIBRARIES ###
pacman::p_load(
  GenomicRanges,
  tidyverse,
  rtracklayer,
  readxl,
  ComplexHeatmap,
  circlize,
  GenomicFeatures
)

setwd(this.path::here())
source("custom_functions.R")

### --- Load preprocessed objects -------------------------------------------
actinopteriigy_cne_gr <- readRDS(
  "../output/preprocessed/actinopteriigy_cne_gr.rds"
)
gnathostomata_cne_gr <- readRDS(
  "../output/preprocessed/gnathostomata_cne_gr.rds"
)
peak_anno_list <- readRDS("../output/preprocessed/peak_anno_list.rds")
gene_activity <- readRDS("../output/preprocessed/gene_activity.rds")

### --- Load rGREAT-derived gene lists --------------------------------------
fg <- readRDS("../output/great/filtered_genes_for_downstream.rds")
filtered_genes_fin <- fg$fin
filtered_genes_skeletal_fin <- fg$fin_skeletal
filtered_genes <- fg$all_signif_actino # for the "any signif GO" filters

great_tbl_actino <- read_tsv("../output/great/great_actinopterygii_GObp.tsv")

drer_anno <- txdbmaker::makeTxDbFromGFF(
  "../ancilliary_files/drer.gff",
  organism = "Danio rerio"
)

### --- YueSong overlap (liftover-based) ------------------------------------
sheet7_data <- read_excel(
  "../input/YueSong-et-al_2025.xlsx",
  sheet = 7,
  col_names = TRUE,
  skip = 1
)
alias <- read_tsv(
  "https://hgdownload.soe.ucsc.edu/hubs/GCF/000/002/035/GCF_000002035.6/GCF_000002035.6.chromAlias.txt",
  comment = "#",
  show_col_types = FALSE,
  col_names = FALSE
)
sheet7_data$Chromosome <- alias$X5[match(sheet7_data$Chromosome, alias$X1)]

sheet7_gr <- GRanges(
  seqnames = sheet7_data[[1]],
  ranges = IRanges(start = sheet7_data[[2]], end = sheet7_data[[3]]),
  mcols = sheet7_data[, 4:ncol(sheet7_data)]
)

liftover <- read_tsv(
  "../input/ucsc_GRCz11-GRCz12_liftover_actinopteriigy_cne.bed",
  col_names = FALSE
)
liftover_gr <- GRanges(
  seqnames = liftover$X1,
  ranges = IRanges(start = liftover$X2 + 1, end = liftover$X3),
  strand = "*",
  cne_name = liftover$X4
)

hits <- findOverlaps(liftover_gr, sheet7_gr, ignore.strand = TRUE)
overlapping_cne_names <- unique(mcols(liftover_gr)$cne_name[queryHits(hits)])

# Need the raw CNE table for chrom/start/end strings used by the UpSet plot
actinopteriigy_cne <- as.data.frame(actinopteriigy_cne_gr) |>
  transmute(
    chromosome = as.character(seqnames),
    start,
    end,
    cne_name = paste0(chromosome, ".1")
  ) # legacy concat for sheet7 key
actinopteriigy_cne_ov <- actinopteriigy_cne |>
  filter(cne_name %in% overlapping_cne_names)

### --- GO-derived gene filtering on peak annotations -----------------------
anno_actino <- non_exon(peak_anno_list$actinopteriigy_CNE@anno)
anno_gnatho <- non_exon(peak_anno_list$gnathostomata_CNE@anno)

actinopteriigy_filtered_skeletal_fin <- as.data.frame(
  anno_actino[
    (anno_actino$geneId %in% filtered_genes_skeletal_fin)
  ]
)
actinopteriigy_filtered_fin <- as.data.frame(
  anno_actino[
    (anno_actino$geneId %in% filtered_genes_fin)
  ]
)

write.table(
  actinopteriigy_filtered_fin,
  "../output/actinopteriigy_GO-fin_specific_cne.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  actinopteriigy_filtered_skeletal_fin,
  "../output/actinopteriigy_GO-fin-skeletal_specific_cne.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

### --- Active-gene filtered exports ----------------------------------------
export_active <- as.data.frame(
  anno_actino[
    anno_actino$is_active == TRUE
  ]
)
write.table(
  export_active,
  "../output/actinopteriigy_cne_near_active_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  filter(actinopteriigy_filtered_fin, is_active == TRUE),
  "../output/actinopteriigy_GO-fin_specific_cne_active_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  filter(actinopteriigy_filtered_skeletal_fin, is_active == TRUE),
  "../output/actinopteriigy_GO-fin-skeletal_specific_cne_active_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  as.data.frame(anno_gnatho[!str_detect(anno_gnatho$annotation, "Exon")]),
  "../output/gnathostomata_specific_cne.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

### --- ATAC-seq overlap ----------------------------------------------------
atac_peaks <- read_tsv(
  "../input/consensus_peaks.mLb.clN.bed",
  col_names = FALSE
)
colnames(atac_peaks) <- c(
  "chromosome",
  "start",
  "end",
  "peak_name",
  "score",
  "strand"
)
atac_peaks_gr <- GRanges(
  seqnames = atac_peaks$chromosome,
  ranges = IRanges(atac_peaks$start, atac_peaks$end),
  peak_name = atac_peaks$peak_name,
  score = atac_peaks$score
)

anno_atac_overlaps <- findOverlaps(anno_actino, atac_peaks_gr)
cne_atac_annotated <- anno_actino[queryHits(anno_atac_overlaps)]
cne_atac_annotated <- cne_atac_annotated[
  !(str_detect(cne_atac_annotated$annotation, "Exon"))
]
cne_atac_annotated_df <- as.data.frame(cne_atac_annotated)

write.table(
  cne_atac_annotated_df,
  "../output/actinopteriigy_cne_atac_overlaps_annotated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  filter(cne_atac_annotated_df, is_active == TRUE),
  "../output/actinopteriigy_cne_atac_overlaps_active_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ATAC × annotation bar chart
overlap_by_annotation <- cne_atac_annotated_df |>
  as_tibble() |>
  mutate(annotation = str_remove(annotation, " \\(.*$")) |>
  dplyr::count(annotation) |>
  mutate(perc = n / sum(n) * 100) |>
  arrange(desc(n))

p <- ggplot(
  overlap_by_annotation,
  aes(reorder(annotation, -n), n, fill = annotation)
) +
  geom_col(color = "black", alpha = 0.8) +
  geom_text(
    aes(label = paste0(n, "\n(", round(perc, 1), "%)")),
    vjust = -0.5,
    size = 3
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x = "Annotation Type",
    y = "Count",
    title = "Distribution of ATAC-overlapping CNEs by Annotation"
  )
ggsave(
  "../output/cne_atac_overlaps_by_annotation.pdf",
  p,
  width = 8,
  height = 6,
  dpi = 300
)

# fin-specific ATAC subsets
write.table(
  filter(cne_atac_annotated_df, geneId %in% filtered_genes_fin),
  "../output/actinopteriigy_cne_atac_overlaps_fin_specific.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  filter(
    cne_atac_annotated_df,
    geneId %in% filtered_genes_fin,
    is_active == TRUE
  ),
  "../output/actinopteriigy_cne_atac_overlaps_fin_specific_active_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

### --- actinopteriigy-unique CNEs ∩ ATAC -----------------------------------
ov_tel_verte <- findOverlaps(actinopteriigy_cne_gr, gnathostomata_cne_gr)
actinopteriigy_only_gr <- actinopteriigy_cne_gr[setdiff(
  seq_along(actinopteriigy_cne_gr),
  unique(queryHits(ov_tel_verte))
)]
ov_tel_atac <- findOverlaps(actinopteriigy_only_gr, atac_peaks_gr)
tel_only_with_atac <- actinopteriigy_only_gr[unique(queryHits(ov_tel_atac))]
ov_anno_telonly <- findOverlaps(anno_actino, tel_only_with_atac)
final_df_tel <- as.data.frame(anno_actino[unique(queryHits(ov_anno_telonly))])

write.table(
  final_df_tel,
  "../output/actinopteriigy_unique_cne_with_atac_annotated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  filter(final_df_tel, is_active == TRUE),
  "../output/actinopteriigy_unique_cne_with_atac_annotated_active_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

### --- GRanges RDS exports for GenometriCorr etc. --------------------------
export_rds_dir <- "../output/granges_rds"
dir.create(export_rds_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(
  actinopteriigy_cne_gr,
  file.path(export_rds_dir, "actinopteriigy_cne_gr.rds")
)
saveRDS(
  gnathostomata_cne_gr,
  file.path(export_rds_dir, "gnathostomata_cne_gr.rds")
)
saveRDS(atac_peaks_gr, file.path(export_rds_dir, "atac_peaks_gr.rds"))
saveRDS(peak_anno_list, file.path(export_rds_dir, "peak_annotations.rds"))
saveRDS(genes(drer_anno), file.path(export_rds_dir, "genome_genes_gr.rds"))
saveRDS(
  rtracklayer::import("../ancilliary_files/drer.gff"),
  file.path(export_rds_dir, "genome_annotation_gr.rds")
)

### --- UpSet sets ----------------------------------------------------------
en <- read_tsv("../ancilliary_files/enhancer.grcz12.bed", col_names = FALSE)
enh_gr <- GRanges(
  seqnames = en$X1,
  ranges = IRanges(en$X2 + 1, en$X3),
  ep_id = en$X4
)

yuesong_gr <- GRanges(
  seqnames = str_remove_all(actinopteriigy_cne_ov$chromosome, "\\.1$"),
  ranges = IRanges(
    start = actinopteriigy_cne_ov$start,
    end = actinopteriigy_cne_ov$end
  )
)

base_gr <- unique(anno_actino)

SLACK <- 0L # adjust as needed — 0 for exact, 50-200 for lifted coords

in_atac <- overlapping_idx(base_gr, atac_peaks_gr, slack = SLACK)
in_active <- overlapping_idx(base_gr, base_gr[base_gr$is_active], 0)
in_fin <- which(base_gr$geneId %in% filtered_genes_fin)
in_sheet7 <- overlapping_idx(
  base_gr,
  yuesong_gr,
  slack = SLACK
)
in_chan_enh <- overlapping_idx(base_gr, enh_gr, slack = SLACK)

res_actino <- analyse_cne_universe(
  anno_gr = unique(anno_actino),
  label = "actinopteriigy",
  atac_peaks_gr = atac_peaks_gr,
  fin_geneIds = filtered_genes_fin,
  slack = 0L
)

res_gnatho <- analyse_cne_universe(
  anno_gr = unique(anno_gnatho),
  label = "gnathostomata",
  atac_peaks_gr = atac_peaks_gr,
  fin_geneIds = filtered_genes_fin,
  slack = 0L
)


### --- Final tables for the Shiny app --------------------------------------

write.table(
  res_actino$final,
  "../output/actinopteriigy_cne_final_table.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  res_gnatho$final,
  "../output/gnathostomata_cne_final_table.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
